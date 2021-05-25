### Nicole Bedford, Hoekstra Lab
### September 2019
### Social Network Graphs:
### SNA of Overlap Index and Genetic Relatedness
#### http://kateto.net/network-visualization
#### https://briatte.github.io/ggnet/
library(plyr)
library(dplyr)
library(tidyr)
library(data.table) # setorder
library(igraph)
library(ggraph)
library(ggplot2)
library(ggforce)
library(Rmisc)

## Set working directory:
setwd("~/Dropbox (Personal)/Fieldwork_Paper_NLB/Good/Code/ForGitHub")

#### 1. Read in Overlap Index dataframe: #### --------------------------------------------------------------------------------------------
OI <- read.csv("OI_1.5h_real.csv", header=T, sep=",")
OI[1] <- NULL

## Make nodes dataframe:
nodes <- OI %>% select(mouse1, sex1) %>% distinct()

## Make links dataframe:
links <- OI %>% select(mouse1, mouse2, OI)

## Save copy of bidirectional links dataframe (each mouse appears in both "from" and "to" columns)
links_bi <- links %>% filter(OI > 0)
colnames(links_bi) <- c("from", "to", "weight") # 148 links

#### 2. Make an undirected links dataframe (i.e., take mean OI for each pair) : #### --------------------------------------------------------------------------------------------
## First make mouse1 and mouse2 into integers:
links$mouse1 <- gsub("uk2", 2, links$mouse1)
links$mouse2 <- gsub("uk2", 2, links$mouse2)
links$mouse1 <- as.integer(links$mouse1)
links$mouse2 <- as.integer(links$mouse2)
links$pair <- ifelse(links$mouse1 < links$mouse2,
                     paste(links$mouse1, links$mouse2, sep="."),
                     paste(links$mouse2, links$mouse1, sep="."))
links_un <- links %>% group_by(pair) %>% summarise(meanOI = mean(OI)) # calculate mean OI for each pair
links_un <- separate(links_un, pair, sep="\\.", into = c("mouse1", "mouse2"))
links_un$mouse1 <- gsub("\\<2\\>", "uk2", links_un$mouse1)
links_un$mouse2 <- gsub("\\<2\\>", "uk2", links_un$mouse2)
links_un$mouse1 <- as.factor(links_un$mouse1) 
links_un$mouse2 <- as.factor(links_un$mouse2)
colnames(links_un) <- c("from", "to", "weight")
links_un <- links_un %>% filter(weight > 0) # 74 links

#### 3. Make SNA graph with genetic relatedness data: #### -------------------------------------------------------------------------------------------------
## Add genetic relatedness to links dataframe: from Brock's vcfTools results
kin <- read.csv("kinship.ACfilter.csv", header=T, sep=",")
decode <- read.csv("SampleDecode.csv", header=T, sep=",", colClasses = "character")
decode$ID <- paste(decode$Sex, decode$MouseID, sep="_") # make new ID variable that includes mouse sex and ID number
info <- decode[ ,c(1,5,6,11,14)] # keep only essential columns
colnames(info) <- c("INDV1", "MouseID_1", "Sex_1", "Age_1", "ID_1")
kin1 <- full_join(info, kin[ ,c(1,2,7)]) # Add Individual 1 information
colnames(info) <- c("INDV2", "MouseID_2", "Sex_2", "Age_2", "ID_2")
kin2 <- full_join(info, kin1)
kin2 <- kin2[ ,c(7,2,6,1,9,4,8,3,11)] # re-order the kin2 dataframe
kin2[,9] <- kin2[,9]*2 # convert kinship coefficient to relatedness (2X)
colnames(kin2) <- c("mouse1", "mouse2", "sample1", "sample2", "age1", "age2", "sex1", "sex2", "relatedness")
kin2 <- kin2[complete.cases(kin2), ] # remove rows containing NA

## Make links dataframe (complete with edges for Overlap Index and Genetic Relatedness):
## First make mouse1 and mouse2 into integers:
kin3 <- kin2 # make copy of kinship dataframe
kin3$mouse1 <- as.integer(kin3$mouse1)
kin3$mouse2 <- as.integer(kin3$mouse2)
kin3$pair <- ifelse(kin3$mouse1 < kin3$mouse2,
                     paste(kin3$mouse1, kin3$mouse2, sep="."),
                     paste(kin3$mouse2, kin3$mouse1, sep="."))
links <- left_join(links, kin3[ ,c(9:10)]) %>% distinct
links$mouse1 <- as.factor(links$mouse1) 
links$mouse2 <- as.factor(links$mouse2)
links_gr <- links[ ,c(4,5)] %>% distinct()
links_gr <- separate(links_gr, pair, sep="\\.", into = c("mouse1", "mouse2"))
links_gr$mouse1 <- gsub("\\<2\\>", "uk2", links_gr$mouse1)
links_gr$mouse2 <- gsub("\\<2\\>", "uk2", links_gr$mouse2)
colnames(links_gr) <- c("from", "to", "weight")
links_gr$weight <- ifelse(links_gr$weight < 0, 0, links_gr$weight)
links_gr[is.na(links_gr)] <- 0
links_gr <- links_gr %>% filter(!from %in% c("126", "uk2"), !to %in% c("126", "uk2"), !weight == 1)

## What is the average genetic relatedness in the population?
range(links_gr$weight)
summarySE(links_gr, measurevar = "weight") # 435 links

## Remove zero links (i.e., unrelated pairs):
links_gr <- links_gr %>% filter(!weight == 0) # 272 links

## Make network graph object using genetic relatedness links:
nodes_gr <- nodes %>% filter(!mouse1 %in% c("126", "uk2")) # subset of nodes dataframe including all mice with genotype information
gr_graph <- graph_from_data_frame(links_gr %>% filter(!from %in% c("126", "uk2"), !to %in% c("126", "uk2")), vertices = nodes_gr)
layout_gr <- create_layout(gr_graph, layout="fr")
coords <- as.data.frame(layout_gr[ ,c(1:2)])

## Make the genetic relatedness graph:
ggraph(gr_graph, layout = "manual", node.positions = coords, circular = FALSE) + 
  geom_edge_link(edge_colour="black", edge_alpha=0.2) +
  geom_node_label(aes(label=nodes_gr$mouse1)) +
  geom_node_point(aes(color=nodes_gr$sex1)) +
  theme_void() +
  theme(aspect.ratio=1)
ggsave("Plots/Fig4a.pdf", width=20, height=16, units="cm", dpi=1500, useDingbats=FALSE)

## Make network graph object using 1.5h OI links:
un_graph <- graph_from_data_frame(links_un %>% filter(!from %in% c("126", "uk2"), !to %in% c("126", "uk2")), vertices = nodes_gr) 

## Make the Overlap Index graph:
ggraph(un_graph, layout = "manual", node.positions = coords, circular = FALSE) + 
  geom_edge_link(edge_colour="black", edge_alpha=0.2, aes(width = weight)) +
  geom_node_label(aes(label=nodes_gr$mouse1)) +
  geom_node_point(aes(color=nodes_gr$sex1)) +
  theme_void() +
  theme(aspect.ratio=1)
ggsave("Plots/Fig4b.pdf", width=20, height=16, units="cm", dpi=1500, useDingbats=FALSE)

## Calculate degree centrality for OI network: How many direct, ‘one hop’ connections each node has to other nodes within the network?
nodes2 <- nodes %>% filter(!mouse1 %in% c("uk2"))
un_graph2 <- graph_from_data_frame(links_un %>% filter(!from %in% c("uk2"), !to %in% c("uk2")), vertices = nodes2) 
nodes2$cd <- centr_degree(un_graph2)$res
colnames(nodes2)[1:2] <- c("MouseID", "Sex")

## Calculate betweenness centrality for OI network: which nodes act as ‘bridges’ between nodes in a network?
nodes2$bc <- betweenness(un_graph2)

## Add in age info:
nodes3 <- decode %>% select(MouseID, Sex, Age) %>% distinct()
nodes3 <- left_join(nodes2, nodes3)
nodes3$Age2 <- ifelse(nodes3$Age == "Juvenile" | nodes3$Age == "Subadult", "Juvenile", "Adult")
summary(lm(cd ~ Sex, data = nodes3)) # 0.309
summary(lm(cd ~ Age2, data = nodes3)) # 0.472
summary(lm(bc ~ Sex, data = nodes3)) # 0.445
summary(lm(bc ~ Age2, data = nodes3)) # 0.416

mean(nodes3$cd) # 4.2

## Calculate degree centrality for GR network: How many direct, ‘one hop’ connections each node has to other nodes within the network? 
nodes2 <- nodes %>% filter(!mouse1 %in% c("uk2"))
gr_graph2 <- graph_from_data_frame(links_gr %>% filter(!from %in% c("uk2"), !to %in% c("uk2")), vertices = nodes2) 
nodes2$cd <- centr_degree(gr_graph2)$res
colnames(nodes2)[1:2] <- c("MouseID", "Sex")

## Calculate betweenness centrality for GR network: which nodes act as ‘bridges’ between nodes in a network?
nodes2$bc <- betweenness(gr_graph2)

## Add in age info:
nodes3 <- decode %>% select(MouseID, Sex, Age) %>% distinct()
nodes3 <- left_join(nodes2, nodes3)
nodes3$Age2 <- ifelse(nodes3$Age == "Juvenile" | nodes3$Age == "Subadult", "Juvenile", "Adult")
summary(lm(cd ~ Sex, data = nodes3)) # 0.325
summary(lm(cd ~ Age2, data = nodes3)) # 0.142
summary(lm(bc ~ Sex, data = nodes3)) # 0.183
summary(lm(bc ~ Age2, data = nodes3)) # 0.953

mean(nodes3$cd) # 17.5

#### 4. Run linear regression with genetic relatedness versus Overlap Index: #### -------------------------------------------------------------------------------------------------
## Plot non-zero overlap indices against relatedness:
reg <- links %>% group_by(pair, relatedness) %>% summarise(meanOI = mean(OI)) # calculate mean OI for each pair
reg <- reg %>% filter(meanOI > 0, relatedness > 0)
r <- lm(meanOI ~ relatedness, data=reg)
summary(r) # p = 0.003322 **, Adjusted R-squared: 0.2284

ggplot(reg, aes(x=relatedness, y=meanOI)) +
  geom_point(size=0.5) +
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE) +
  theme_classic()
ggsave("Plots/Fig4c.pdf", width=10, height=8, units="cm", dpi=1500, useDingbats=FALSE)
