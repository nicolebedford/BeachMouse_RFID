#### Nicole Bedford, Hoekstra Lab 2019
#### Code for Florida Fieldwork Paper Analyses
#### Make Figures 1, 2, & 5 using RFID data

## Read-in the data:
## /Users/nlbedford/Dropbox/Fieldwork_Paper_NLB/RFID/MasterDataframe/RFID_Master_dataframe.R
## Run the whole script from RFID_Master_dataframe.R to get the "nov" and "may" data frames:
## Clear the rest of the workspace by switching List to Grid, select all, remove nov and may, then Clear Workspace

## Load libraries:
library(plyr)
library(dplyr)
library(Rmisc)
library(ggplot2)
library(lme4)
library(lmerTest)
library(grid)
library(tibble)
library(factoextra)
library(gdata)
library(cluster)

## Set working directory:
setwd("~/Dropbox (Personal)/Fieldwork_Paper_NLB/Good/Code/ForGitHub")

#### Start with November 2017 Dataset #### --------------------------------------------------------------------------------------------------------------------------------
## Read-in the RFID data: (run RFID_Master_dataframe.R)
copy <- nov # make a copy of the November dataframe

#### 1. Remove "ping trains" #### ------------------------------------------------------------------------------------------------------------------------------------------
## If pings occur at 1 Hz, keep first and last in the "train"
copy <- copy %>% filter(!is.na(Night)) # keep only Nights 1 to 11 (no grid traps in place)
copy <- copy %>% filter(!Mouse_ID %in% c("uk1", "uk3", "uk4")) # these are "starter" RFID tags --> not real mice!
pingSum <- as.data.frame(table(copy$Mouse_ID)) # number of pings per animal before train filtering
preFilt_mouseSummary <- copy %>% group_by(Mouse_ID) %>% 
  summarise(numNight=n_distinct(Night), numBur=n_distinct(Burrow), pingSum=n())
low <- pingSum %>% filter(Freq < 4) # find mice with fewer than 4 pings (120, 122, 110, 113, 119, 123)
copy_high <- copy %>% filter(!Mouse_ID %in% low$Var1) # remove "low" ping mice (for-loop doesn't work for mice with < 4 pings)
highMice <- as.data.frame(factor(unique(copy_high$Mouse_ID))) # "high" ping mice (4 or more pings)
colnames(highMice) <- "Mouse_ID"
ord <- order(highMice)
m <- as.data.frame(cbind(ord, highMice))
colnames(m) <- c("ord", "Mouse_ID")
copy_high <- full_join(copy_high, m)

## Per night summary of RFID data for Supplementary Table 2:
resultlist <- list()
for(i in 1:11) {
  p <- copy %>% filter(Night == i) %>% summarise(burs = n_distinct(Burrow), mice = n_distinct(Mouse_ID), pings = n()) 
  p$Night <- i  
  resultlist[[i]] = p 
}
result = do.call(rbind, resultlist)
copy %>% summarise(burs = n_distinct(Burrow), mice = n_distinct(Mouse_ID), pings = n()) # Summary for all 11 Nights

## Nested for-loop for removing ping trains for each individual:
resultlist <- list()
for (i in 1:nrow(highMice)) {
  p <- copy_high %>% filter(ord %in% i)
  p$diff <- NA # 22nd column
  p$remove <- NA # 23rd column
  p <- p[order(p$Timestamp),]
  for (j in 2:(nrow(p)-1)) {
    p[j,22] <- ifelse(p[j+1,4] - p[j,4] == 1, 1, 0)
    p[j,23] <- ifelse(p[j-1,22] == 1 & p[j,22] == 1, "remove", "keep")
    resultlist[[i]] = p 
  }}
result = do.call(rbind, resultlist) # make sure result has same number of rows as copy!
result <- full_join(copy, result)
table(result$remove) # 7452 keep, 794 remove
trim <- result %>% filter(!remove %in% "remove") %>% droplevels() # good dataframe for all downstream analyses
length(unique(trim$Mouse_ID)) # 32
length(unique(trim$Burrow)) # 40

#### 2. Read in the mark-recapture data: #### ----------------------------------------------------------------------------------------------------
capt <- read.table("capt17.txt")
GridTrap <- read.table("GridTrap17.txt")
colnames(capt)[4] <- "Detector"
colnames(GridTrap)[1] <- "Detector"
traps <- join(capt, GridTrap, by="Detector", type="left")
traps <- traps[,c(2,4:7)] # keep only necessary columns
colnames(traps) <- c("ID", "Detector", "Sex", "X", "Y")
traps <- unique(traps[c("ID", "Detector", "Sex", "X", "Y")]) # for each individual, throw out duplicate trap locations
n.per.detector <- traps %>% group_by(Detector, X, Y) %>% summarise(n=n())
n.per.detector$n <- as.factor(n.per.detector$n)

#### 3a. Summarise: How many pings per burrow over 11 nights? #### ------------------------------------------------------------------------------------------------------
burTot <- trim %>% group_by(Burrow, X, Y) %>% summarise(pingSum = n())
burTot$scale <- log2(as.numeric(burTot$pingSum) + 1) # if pingSum = 1, scale = 0

### Figure 1:
ggplot() +
  geom_point(data=burTot, aes(x=X, y=Y, size=scale)) + 
  geom_point(data=n.per.detector, aes(x=X, y=Y, color=n)) + # successful trap locations (n=69)
  xlab(NULL) + ylab(NULL) + 
  scale_x_continuous(breaks=seq(0,1300,100), limits=c(-5,1300)) + 
  scale_y_continuous(breaks=seq(0,100,20), limits=c(-20,125)) +
  coord_fixed() +
  theme_classic() + guides(fill = FALSE, color = FALSE)
ggsave("Plots/Fig1d.pdf", height = 12, units="cm", dpi=1500, useDingbats=FALSE)

#### 3b. Summarise: How many burrows does a mouse visit each night? #### ----------------------------------------------------------------------------------------------------
## Clean up dataframe:
df <- trim %>% 
  select(Night, Burrow, Mouse_ID, Sex, Age)
df$Mouse_ID <- as.factor(df$Mouse_ID)
df$Burrow <- as.factor(df$Burrow)

# Nice summary of specs per burrow
bp <- df %>% group_by(Burrow) %>% summarise(totPings=n())
bn <- df %>% group_by(Burrow) %>% summarise(totNights=n_distinct(Night))
tm <- df %>% group_by(Burrow) %>% summarise(totMice=n_distinct(Mouse_ID))
Burrow <- full_join(bp, bn, by="Burrow") %>% full_join(., tm, by="Burrow")

# Nice summary of specs per mouse
mp <- df %>% group_by(Mouse_ID, Sex, Age) %>% summarise(totPings=n())
mn <- df %>% group_by(Mouse_ID, Sex, Age) %>% summarise(totNights=n_distinct(Night))
tb <- df %>% group_by(Mouse_ID, Sex, Age) %>% summarise(totBurs=n_distinct(Burrow))
Mice <- full_join(mp, mn, by=c("Mouse_ID", "Sex", "Age")) %>% full_join(., tb, by=c("Mouse_ID", "Sex", "Age")) 

df_noDup <- df[!duplicated(df), ]
length(unique(df_noDup$Burrow)) # 40 burrows
length(unique(df_noDup$Mouse_ID)) # 32 mice 

#### 4. Consistency: How repeatable across nights are visition patterns for mice and burrows? ------------------------------------------------------------------------------
## a. How consistent are mice in their burrow usage across nights?
df_noDup <- with(df_noDup, df_noDup[order(Mouse_ID, Night, Burrow),])
burTally <- df_noDup %>% group_by(Mouse_ID, Burrow) %>% tally()
burTally$n <- burTally$n -1 # subtract first obs. (i.e., only count repeated burrows)
burTally <- burTally %>% group_by(Mouse_ID) %>% summarise(burTally = sum(n))
Mice <- full_join(Mice, burTally)
Mice$consistency <- Mice$burTally / ((Mice$totNights -1) * Mice$totBurs)
Mice_comp <- Mice[complete.cases(Mice), ]
mean(Mice_comp$consistency) # 47%
fo <- Mice_comp %>% filter(totNights == 11)
mean(fo$consistency) # 55%
median(fo$consistency) # 56%

ggplot(Mice %>% filter(totNights == 11), aes(x=consistency)) +
  geom_histogram(binwidth = 0.1, colour="black", fill="grey") +
  ylab("no. mice") + xlab("Burrow use consistency")

## b. How consistent are burrows in their mouse visitation patterns across nights?
mouseTally <- df_noDup %>% group_by(Burrow, Mouse_ID) %>% tally()
mouseTally$n <- mouseTally$n -1 # subtract first obs. (i.e., only count repeated mice)
mouseTally <- mouseTally %>% group_by(Burrow) %>% summarise(mouseTally = sum(n))
Burrow <- full_join(Burrow, mouseTally)
Burrow$consistency <- Burrow$mouseTally / ((Burrow$totNights -1) * Burrow$totMice)
Burrow[is.na(Burrow)] <- 0
mean(Burrow$consistency) 
fo <- Burrow %>% filter(totNights == 11)
mean(fo$consistency) 

ggplot(Burrow %>% filter(totNights == 11), aes(x=consistency)) +
  geom_histogram(binwidth = 0.1, colour="black", fill="grey") +
  ylab("no. burrows") + xlab("Mouse visitation consistency")

#### 5a. Burrow summaries: How many mice per burrow per night? How many mice per burrow total? -----------------------------------------------------------------------------------
## Mice per night
mpn <- df_noDup %>% group_by(Burrow, Night) %>% summarise(numMice=n_distinct(Mouse_ID), 
                                                          numFem=sum(Sex=="F"), numMal=sum(Sex=="M"),
                                                          numJuv=sum(Age=="Juvenile"), numAd=sum(Age=="Adult")) 

mpn2 <- mpn %>% group_by(Burrow) %>% summarise(micePerNight=mean(numMice))
range(mpn2$micePerNight) # 1 - 5.14
mean(mpn2$micePerNight) # 2.21
hist(mpn2$micePerNight)

ggplot(Burrow, aes(x=totMice)) +
  geom_histogram(binwidth = 1, colour="black", fill="grey") +
  scale_x_continuous(breaks=seq(0,10,2)) + 
  scale_y_continuous(breaks=seq(0,10,2), limits=c(0,10)) + 
  ylab("no. burrows") + xlab("total mice per burrow") +
  theme_classic()
ggsave("Plots/Fig2e.pdf", width=8, height=6, units="cm", dpi=1500, useDingbats=FALSE)

range(Burrow$totMice) # 1-10
median(Burrow$totMice) # 3

#### 5b. Mouse summaries: how many burrows per mouse per night? how many burrows per mouse total? ----------------------------------------------------------------------------
## Burrows per night
bpn <- df_noDup %>% group_by(Mouse_ID, Night) %>% summarise(numBurs=n_distinct(Burrow))

bpn2 <- bpn %>% group_by(Mouse_ID) %>% summarise(burrowsPerNight=mean(numBurs))
range(bpn2$burrowsPerNight) # 1 - 5.6
mean(bpn2$burrowsPerNight) # 2.45
hist(bpn2$burrowsPerNight)
                                                         
ggplot(Mice, aes(x=totBurs)) +
  geom_histogram(binwidth = 1, colour="black", fill="grey") +
  scale_x_continuous(breaks=seq(0,10,2)) + 
  scale_y_continuous(breaks=seq(0,8,2), limits=c(0,8)) + 
  ylab("no. mice") + xlab("total burrows per mouse") +
  theme_classic()
ggsave("Plots/Fig2f.pdf", width=8, height=6, units="cm", dpi=1500, useDingbats=FALSE)

range(Mice$totBurs) # 1-9
median(Mice$totBurs) # 5

#### 6a. Burrow Membership by ping portion #### -----------------------------------------------------------------------------------------------------------------------------
perNight_perMouse_pingSummary <- trim %>% group_by(Burrow, X, Y, Night, Mouse_ID, Sex) %>% dplyr::summarise(pingCount = n())
perNight_pingSummary <- trim %>% group_by(Burrow, X, Y, Night) %>% dplyr::summarise(pingCount = n())

## Averge number of pings per night (per mouse)?
mouseNightly <- perNight_perMouse_pingSummary %>% group_by(Mouse_ID, Sex, Night) %>% dplyr::summarise(pingSum = sum(pingCount))
mouseNightly <- mouseNightly %>% group_by(Mouse_ID, Sex) %>% summarise(meanPings_perNight = mean(pingSum))
hist(mouseNightly$meanPings_perNight)
mean(mouseNightly$meanPings_perNight) # 24

# What portion of a burrow's total pings are made up by indiv. mice?
pingSummary <- perNight_perMouse_pingSummary %>% group_by(Burrow, X, Y, Mouse_ID, Sex) %>% dplyr::summarise(pingSum = sum(pingCount))

#### 6b. Equal Usage? Make plot of activity level by burrow rank #### -------------------------------------------------------------------------------------------------------
# pingSum = no. pings for that mouse @ that burrow
# pingTot = no. pings for that mouse total (all burrows, all nights)
pingSummary <- pingSummary[order(pingSummary$Mouse_ID, -pingSummary$pingSum), ] # order by mouse, then pingSum descending
pingSummary <- pingSummary %>% group_by(Mouse_ID) %>% dplyr::mutate(id = row_number()) # add burrow rank for each Mouse_ID
colnames(pingSummary)[7] <- "burRank"
pingTot <- pingSummary %>% group_by(Mouse_ID) %>% summarise(pingTot=sum(pingSum))
pingSummary <- full_join(pingSummary, pingTot)
pingSummary$activityPortion <- pingSummary$pingSum / pingSummary$pingTot 

## Calculate cumulative activity:
resultlist <- list()
for (i in 1:9) {
  c <- pingSummary %>% filter(burRank %in% c(1:i)) %>%
    group_by(Mouse_ID) %>% summarise(cum_activity = sum(activityPortion))
  c$rank <- i
  resultlist[[i]] = c }
cum_activity = do.call(rbind, resultlist)

## Summarise cumulative activity:
cum_activity_sum <- cum_activity %>% group_by(rank) %>% summarise(mean=mean(cum_activity), n=n(), sd=sd(cum_activity), se=sd/sqrt(n))

scree <- cum_activity_sum %>% mutate(diff = mean - lag(mean))
scree[1,6] <- scree[1,2]

ggplot(scree, aes(x=rank, y=diff)) +
  geom_col(colour="black", fill="grey") +
  geom_errorbar(aes(ymin = diff - se, ymax = diff + se)) +
  scale_x_continuous(breaks=seq(0,10,2)) + 
  scale_y_continuous(breaks=seq(0, 0.6, 0.2), limits=c(0, 0.7005)) +
  ylab("Activity portion") + xlab("Burrow rank") +
  theme_classic() +
ggsave("Plots/Fig2g.pdf", width=8, height=6, units="cm", dpi=1500, useDingbats=FALSE)

#### 6c. Hierarchical clustering: cluster burrows based on similarity in membership #### -----------------------------------------------------------------------------------
# totPings = no. of pings @ that burrow (all mice, all nights)
pie <- full_join(pingSummary, Burrow[,c(1:4)])
pie$portion <- pie$pingSum / pie$totPings
pie$scale <- log2(pie$totPings + 1)
pie$Mouse_ID <- as.factor(pie$Mouse_ID)
pie$Burrow <- as.factor(pie$Burrow)
pie <- pie %>%
  dplyr::arrange(Sex, Mouse_ID) %>% ungroup() %>% # sort dataframe by Sex, then Mouse_ID
  dplyr::mutate(Mouse_ID = factor(Mouse_ID, unique(Mouse_ID))) # apply order to mice
levels(pie$Mouse_ID) # 32 mice
levels(pie$Burrow) # 40 burrows
pie_trim <- pie[ ,c(1:5,13)]
#write.csv(pie, "pie.csv")

## Make visitation matrix (expressed as portion of total burrow activity)
pie_trim_wide <- pie_trim[ ,c(1:4,6)] %>% spread(Mouse_ID, portion)
pie_trim_wide[is.na(pie_trim_wide)] <- 0
pie_trim_wide$sum <- rowSums(pie_trim_wide[ ,4:ncol(pie_trim_wide)]) # make sure this adds up to 1!
pie_trim_wide <- pie_trim_wide[ ,-c(1:3,ncol(pie_trim_wide))]
pie_trim_wide <- as.matrix(pie_trim_wide)
row.names(pie_trim_wide) <- Burrow$Burrow
pie_trim_wideScaled <- scale(pie_trim_wide) # scale/standardize the data (mean=0, sd=1)

## Which Agglomerative clustering method is best?
m <- c("average", "single", "complete", "ward")
names(m) <- c("average", "single", "complete", "ward")

## Function to compute coefficients (higher = more cluster structure)
ac <- function(x) {
  agnes(pie_trim_wideScaled, method = x)$ac
  }
map_dbl(m, ac) # ward is best

d <- dist(pie_trim_wideScaled, method = "euclidean") # compute the dissimilarity matrix
clusters <- hclust(d, method="ward.D2") # ward.D minimizes the total within-cluster variance
heatmap(pie_trim_wideScaled, distfun=dist, hclustfun=function(d) hclust(d, method="ward.D2"))

## Define Clusters : What is the optimal number of clusters?
#https://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
# fviz_nbclust(pie_trim_wideScaled, FUN = hcut, method = "silhouette", k.max = 39) # 19
# fviz_nbclust(pie_trim_wideScaled, FUN = hcut, method = "wss", k.max = 39) # unclear
# gap_stat <- clusGap(pie_trim_wideScaled, FUN = hcut, K.max = 39, B = 50)
# fviz_gap_stat(gap_stat) # 39

k=7 # define number of clusters
clusterCut <- cutree(clusters, k=k)
cl <- as.data.frame(table(clusterCut, Burrow$Burrow))
cl <- cl %>% filter(Freq > 0) # which burrows belong to which clusters?
colnames(cl)[2] <- "Burrow"
pie2 <- join(pie, cl[,c(1,2)])
clustPerMouse <- pie2 %>% group_by(Mouse_ID) %>% summarise(numClust = n_distinct(clusterCut))
table(clustPerMouse$numClust)
bursPerClust <- cl %>% group_by(clusterCut) %>% summarise(n=n())
mean(bursPerClust$n)

## Visualize the dendrogram:
pdf("Plots/Fig2a_top.pdf", width=18, height=5)
plot(clusters, hang = -1)
rect.hclust(clusters, k=k, border=1:length(unique(cl$clusterCut))) # draw boxes around the clusters on the dendrogram
dev.off()

## What is the compostition of each cluster? How many burrows? How many mice?
clustSum <- pie2 %>% group_by(clusterCut) %>% 
  summarise(numBurs=n_distinct(Burrow), numMice=n_distinct(Mouse_ID))

## Who is the top-ranked mouse per cluster (i.e. mouse with the largest share of pings?)
topMouse <- pie2 %>% group_by(clusterCut, Mouse_ID) %>% summarise(clusterPings = sum(pingSum)) # top mouse by portion
topMouse <- topMouse %>% group_by(clusterCut) %>% filter(clusterPings == max(clusterPings))
topMouse <- full_join(topMouse, clustSum)
topMouse2 <- pie2 %>% group_by(clusterCut, Mouse_ID) %>% summarise(numBur = n())

## Plot burrow clusters on trapping grid: 
clustBur <- pie2 %>% select(Burrow, X, Y, clusterCut, scale)
clustBur <- clustBur[!duplicated(clustBur), ]
clustBur$clusterCut <- as.factor(clustBur$clusterCut) 

ggplot() +
  geom_point(data=clustBur, aes(x=X, y=Y, size=scale, color=clusterCut)) +
  theme_classic() + guides(fill = FALSE, color = FALSE)
ggsave("Plots/Fig2b.pdf", height = 6, width = 24, units="cm", dpi=1500, useDingbats=FALSE)

#### 6d. Make Pie Charts #### -----------------------------------------------------------------------------------------------------------------------------------------------
c <- as.dendrogram(clusters)
order <- as.data.frame(labels(c))
order$x <- rownames(order)
colnames(order)[1] <- "Burrow"
pie <- join(pie, order)
pie$x <- as.integer(pie$x)
pie$y <- 1

ID <- unique(pie[ ,4:5])
table(ID$Sex) # 13 Female, 18 Male, 1 Unknown
fem <- c("#BC49BC", "#CD0074", "#9E0059", "#E40045", "#EB5B87", "#FF0000", "#C50000",
         "#FF4900", "#FF9063", "#FF7400", "#C55900", "#FF9200", "#FFBC63")
mal <- c("#3914AF", "#2B0E88", "#1B1BB3", "#5F5FC6", "#1240AB", "#0D4184", "#0B61A4", "#508CBB", "#009999", "#007676",
         "#00AF64", "#4CC390", "#00CC00", "#009E00", "#67E300", "#9CEA5B", "#9FEE00", "#7BB800")
pal <- c(fem, mal, "grey") # 32 mice, 32 colours

ggplot(pie, aes(x = scale/2, y = portion, fill = Mouse_ID, width = scale)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = pal) +
  facet_grid(y ~ x) + # manually assigned coordinates
  coord_polar("y") +
  xlab(NULL) + ylab(NULL) +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.background = element_blank())
ggsave("Plots/Fig2a_bottom.pdf", width=40, height=2, units="cm", dpi=1500, useDingbats=FALSE)

# Extract the pie chart legend:
p <- ggplot(pie, aes(x = scale/2, y = portion, fill = Mouse_ID, width = scale)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = pal) +
  facet_grid(y ~ x) +
  coord_polar("y")
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
mylegend <- g_legend(p)
grid.draw(mylegend)

#### 7. Which burrows are mice sleeping in during the day? ### -----------------------------------------------------------------------------------------------------------------------------------
df <- trim %>% dplyr::select(Burrow, Mouse_ID, ord, Sex, Timestamp, Night) %>% distinct()

firstLast <- function(i) {
  fo <- df %>% filter(ord %in% i)
  fo <- fo[order(fo$Timestamp),]
  result <- fo %>% group_by(Night) %>% slice(c(1, n())) %>% ungroup()
  result$time <- rep(c("first", "last"), len = nrow(result))
  return(data.frame(result))
}
result <- lapply(ord, firstLast) %>% bind_rows()
result <- result[order(result$Mouse_ID, result$Night), ] # order result dataframe first by Mouse_ID, then by Night

## Do animals emerge the next day (at dusk) from the same burrow they were last seen at the night before (at dawn)?
result$Day <- as.factor(ifelse(result$time == "last", result$Night + 0.5, result$Night - 0.5)) # make new column "Day" (e.g., Day 1.5 is the daytime period between Nights 1 and 2)
lastFirst <- result[-6] %>% group_by(Mouse_ID, Day, Burrow) %>% spread(time, Timestamp) # remove "Night", make into wide format
lastFirst_diff <- lastFirst[!complete.cases(lastFirst), ] # when do mice exit from a different burrow than the one they entered? (218 times)
lastFirst_same <- lastFirst[complete.cases(lastFirst), ] # when do mice exit from the same burrow as the one they entered? (202 times)
colnames(lastFirst_same)[6:7] <- c("exit", "enter")

## When are mice exiting and entering their home burrows?
lastFirst_same$inBurrow <- as.numeric((lastFirst_same$exit - lastFirst_same$enter)/60) # how many hours are mice spending in the burrow?
hist(lastFirst_same$inBurrow, breaks=100)

lastFirst_same$Time_exit <- format(as.POSIXct(strptime(lastFirst_same$exit, "%Y-%m-%d %H:%M:%S", tz="")), format = "%H:%M:%S")
lastFirst_same$Time_exit <- as.POSIXct(lastFirst_same$Time_exit, format="%H:%M:%S")
hist(lastFirst_same$Time_exit, breaks=100)
lastFirst_same$Time_enter <- format(as.POSIXct(strptime(lastFirst_same$enter, "%Y-%m-%d %H:%M:%S", tz="")), format = "%H:%M:%S")
lastFirst_same$Time_enter <- as.POSIXct(lastFirst_same$Time_enter, format="%H:%M:%S")
hist(lastFirst_same$Time_enter, breaks=100)

HB <- as.data.frame(lastFirst_same)
HB <- join(HB, Mice[ ,c(1:6)])
HB <- HB %>% filter(Time_exit > paste(Sys.Date(), "12:00:00", sep=" ") & Time_enter < paste(Sys.Date(), "12:00:00", sep=" ")) # Sunrise/Sunset on Nov. 25th 2017 = 6:20AM - 4:46PM: 10:26 daylight (Update to Today's Date!!)
HomeBurrows <- HB %>% group_by(Burrow) %>% summarise(days = n_distinct(Day)) # List of "home burrows" and number of days for which it served as "home" for 1 or more mice
write.csv(HomeBurrows, "HomeBurrows.csv")

HB_good <- HB %>% filter(totNights == 11) %>% drop.levels() # keep only mice that were recorded on all 11 nights
clustOrd <- c("108", "97", "104", "105", "125", "124", "116", "92", "106", "129", "118", "95") # put top 12 mice in good order
HB_good$Mouse_ID <- reorder.factor(HB_good$Mouse_ID, new.order=clustOrd) # apply cluster order to Mouse_ID
mean(HB_good$Time_exit) # "2019-10-22 17:48:11 EDT"
mean(HB_good$Time_enter) # "2019-10-22 05:13:46 EDT"

times <- gather(HB_good[ ,c(9:10)], key = "Entry", value = "Time")
times$Time <- format(as.POSIXct(strptime(times$Time, "%Y-%m-%d %H:%M:%S", tz="")), format = "%H:%M:%S")
times$Time <- as.POSIXct(times$Time, format="%H:%M:%S")
times$Time_num <- as.numeric(times$Time)
times$Time_num <- times$Time_num - min(times$Time_num)
times$Time_num <- ifelse(times$Time_num > 43200, times$Time_num - 86400, times$Time_num) # Noon 12:00:00 = 43200 ; Midnight 00:00:00 = 86400
times$Time_num <- times$Time_num + (60*60 + 49*60 + 35) # 0 is 01:49:35 --> fix
24*(max(times$Time_num) - min(times$Time_num)) / 86400 # 12.96 hours

## Sunrise/Sunset on Nov. 25th 2017 = 6:20AM - 4:46PM
sunrise17 = 0 + 6*60*60 + 20*60 # midnight plus 6h20m
sunset17 = -43200 + 4*60*60 + 46*60 # noon plus 4h46m

ggplot(data = times, aes(x = Time_num, group = Entry, fill = Entry)) +
  geom_histogram(colour="black", bins=24) +
  geom_vline(xintercept = sunrise17) +
  geom_vline(xintercept = sunset17) +
  scale_x_continuous(breaks=seq(-43200, 43200, 7200), limits=c(-43200, 43200)) + # 2h increments
  theme_classic()
ggsave("Plots/FigS4b.pdf", width=18, height=12, units="cm", dpi=1500, useDingbats=FALSE)

# All 5 burrows in this plot were recorded on all 11 nights (check Burrow dataframe)
ggplot(HB_good, aes(x=Day, y=Mouse_ID, fill=Burrow)) +
  geom_tile(color="white") +
  coord_fixed(ratio=0.5) +
  theme_classic()
ggsave("Plots/Fig4d.pdf", width=10, height=6, units="cm", dpi=1500, useDingbats=FALSE)

## When do mice use a home burrow vs. an alternate burrow?
HB_good_long <- gather(HB_good[ ,c(1:2, 4, 11, 6:7)], key = "entry", value = "Timestamp", exit:enter)
foo <- full_join(HB_good_long, trim[ ,c(4:6, 17:18)])
foo <- foo %>% filter(Mouse_ID %in% clustOrd) %>% drop.levels() # keep only "home burrow" mice
foo$home <- ifelse(foo$Mouse_ID %in% c("95", "118", "129", "106") & foo$Burrow == "H15c", "home",
                           ifelse(foo$Mouse_ID %in% c("92", "116") & foo$Burrow == "G6a", "home",
                                  ifelse(foo$Mouse_ID %in% c("124", "125") & foo$Burrow == "E70b", "home",
                                         ifelse(foo$Mouse_ID %in% c("104", "105") & foo$Burrow == "F62", "home",
                                                ifelse(foo$Mouse_ID %in% c("97", "108") & foo$Burrow == "B22", "home", "alternate")))))

foo$Time <- format(as.POSIXct(strptime(foo$Timestamp, "%Y-%m-%d %H:%M:%S", tz="")), format = "%H:%M:%S")
foo$Time <- as.POSIXct(foo$Time, format="%H:%M:%S")
foo$Time_num <- as.numeric(foo$Time)
foo$Time_num <- foo$Time_num - min(foo$Time_num)
foo$Time_num <- ifelse(foo$Time_num > 43200, foo$Time_num - 86400, foo$Time_num) # Noon 12:00:00 = 43200 ; Midnight 00:00:00 = 86400
foo$Time_num <- foo$Time_num + 40 # 0 is 00:00:40 --> fix
24*(max(foo$Time_num) - min(foo$Time_num)) / 86400 # 12.99 hours

ggplot(data = foo, aes(x = Time_num, group = home, fill = home)) +
  geom_histogram(colour="black", bins=24) +
  geom_vline(xintercept = sunrise17) +
  geom_vline(xintercept = sunset17) +
  scale_x_continuous(breaks=seq(-43200, 43200, 7200), limits=c(-43200, 43200)) + # 2h increments
  theme_classic()
ggsave("Plots/FigS4c.pdf", width=18, height=12, units="cm", dpi=1500, useDingbats=FALSE)

foob <- foo %>% group_by(Mouse_ID, Sex, Age, home) %>% summarise(n=n())
foob <- foob %>% spread(key = "home", value = "n")
foob$ratio <- 1 - (foob$home / (foob$home + foob$alternate))
summarySE(foob, measurevar = "ratio") # what % of activity is seen at alternate burrow? 52 +/- 7

## Activity rasters for Top 12 mice:
for (i in clustOrd) {
  print(ggplot(data = foo %>% filter(Mouse_ID == i), aes(x=Timestamp, y=0, color=home)) +
    geom_tile() + ggtitle(i) + theme_classic())
}

#### 8. Genetic relatedness #### ---------------------------------------------------------------------------------------------------------------------------
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

## Look only at home burrow mice:
kin3 <- kin2 %>% filter(mouse1 %in% unique(HB_good$Mouse_ID), mouse2 %in% unique(HB_good$Mouse_ID))
str(kin3)
kin3$mouse1 <- as.integer(kin3$mouse1)
kin3$mouse2 <- as.integer(kin3$mouse2)
kin3$pair <- ifelse(kin3$mouse1 < kin3$mouse2,
                    paste(kin3$mouse1, kin3$mouse2, sep="."),
                    paste(kin3$mouse2, kin3$mouse1, sep="."))
kin3 <- kin3 %>% filter(relatedness < 1) # remove relatedness with self (r = 1.0)
kin3$relatedness <- ifelse(kin3$relatedness < 0, 0, kin3$relatedness) # convert negative relatedness values to zero

burs <- HB_good[ ,1:2]
burs <- burs[!duplicated(burs), ]
burs$Mouse_ID <- as.character(burs$Mouse_ID)
colnames(burs)[2] <- "mouse1"
kin3$mouse1 <- as.character(kin3$mouse1)
kin3 <- full_join(kin3, burs)
colnames(kin3)[11] <- "bur1"
colnames(burs)[2] <- "mouse2"
kin3$mouse2 <- as.character(kin3$mouse2)
kin3 <- full_join(kin3, burs)
colnames(kin3)[12] <- "bur2"

## Are mice in the same or different burrows?
kin3$shared <- as.factor(ifelse(kin3$bur1 == kin3$bur2, "same", "diff"))
kin4 <- kin3[ ,c(9,10,13)]
kin4 <- kin4[!duplicated(kin4), ]

ggplot(kin4, aes(x=relatedness, fill=shared)) + geom_density(alpha=0.3) +
  theme_classic()
ggsave("Plots/Fig4e.pdf", width=8, height=6, units="cm", dpi=1500, useDingbats=FALSE)

kin4_wide <- kin4 %>% spread(key="shared", value="relatedness")
ks.test(kin4_wide$diff, kin4_wide$same) # D = 0.62857, p-value = 0.00245 
  
kin4_same <- kin4 %>% filter(shared == "same")
kin4_same$kinship_coeff <- kin4_same$relatedness / 2 # kinship coefficient is relatedness divided by two
kin4_same$burrow <- c("G6a", "H15c", "H15c", "H15c", "B22", "F62", "H15c", "H15c", "H15c", "E70b")
kin4_kc <- kin4_same %>% group_by(burrow) %>% summarise(mean_kc = mean(kinship_coeff))

#### 9. Home ranges: Minimum convex hull for each mouse #### ---------------------------------------------------------------------------------------------------------------
library(pracma)

par(mfrow=c(2,2))
resultlist <- list()
for (i in unique(nov$Mouse_ID)) {
  s <- pie %>% filter(Mouse_ID == i) # subset the dataframe by Mouse_ID (look only at HB mice)
  X <- data.matrix(s[ ,c(2:3)]) # turn burrow coordinates into matrix
  plot(pie$X, pie$Y, pch = 20, cex = 0.5) # plot all burrow locations
  points(X, col="red", pch = 20, cex = 0.5) # plot burrow locations for Mouse i
  hpts <- chull(X) # calculate the convex hull
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ]) # plot the outline of the convex hull
  a <- as.data.frame(X[hpts, ]) # vertices of the convex hull
  area <- polyarea(a$X, a$Y) # calculate the area of the MCP
  title(paste(i, round(area))) # add to plot title
  a$Mouse_ID <- i # save Mouse_ID
  a$MCP_area <- abs(area) # absolute value of MCP area
  resultlist[[i]] = a
}
mcp = do.call(rbind, resultlist) 
mcp <- join(mcp, Mice[ ,c(1:3,6)])
mcp <- mcp[!duplicated(mcp), ]
mcp <- mcp %>% filter(MCP_area > 0)
mcp <- mcp %>% group_by(Mouse_ID, Sex, Age, MCP_area, totBurs) %>% summarise(MCP_vert = n())

## Manually count the number of burrows located inside the MCP: how many used (red) and unused (black)?
useCount <- read.csv("burrowUse_inside_MCP.csv", sep=",", header = TRUE)
mcp2 <- full_join(mcp, useCount)
mcp2$Burs_available <- mcp2$MCP_vert + mcp2$inMCP
mcp2$frac_used <- mcp2$Burs_used / mcp2$Burs_available 
hist(mcp2$frac_used)
range(mcp2$frac_used) # 0.5 - 1.0 (50% - 100%)
mean(mcp2$frac_used, na.rm = T) # On average, mice use 89% +/- 5% of the burrows in their home range
sd(mcp2$frac_used, na.rm = T)/sqrt(length(mcp2))
plot(mcp$MCP_area, mcp$frac_used) # The larger the home range, the more unused burrows!

#### 10. Distance between burrows belonging to same cluster #### ---------------------------------------------------------------------------------------------------------------

# Set up a fake key to join on (just a constant):
clustBur2 <- clustBur[ ,c(1:3)] %>% mutate(k = 1) 

# Perform the join, remove the key, then calculate the distance:
clustBur2 <- clustBur2 %>% 
  full_join(clustBur2, by = "k") %>% 
  filter(Burrow.x != Burrow.y) %>%
  mutate(dist = sqrt((X.x - X.y)^2 + (Y.x - Y.y)^2)) %>%
  select(-k)

# Split cluster 3 (far left in Fig. 2a dendrogram)
info <- clustBur
info$clusterCut <- ifelse(info$Burrow == "D9a" | info$Burrow == "A9", 
                   gsub("3", "8", info$clusterCut),
                   info$clusterCut)
info$clusterCut <- ifelse(info$Burrow == "C49", 
                          gsub("3", "9", info$clusterCut),
                          info$clusterCut)
info$clusterCut <- ifelse(info$Burrow == "I31" | info$Burrow == "H31a", 
                          gsub("3", "10", info$clusterCut),
                          info$clusterCut)

# Add back in clusterCut information:
colnames(info)[1] <- "Burrow.x"
clustBur2 <- full_join(clustBur2, info[ ,c(1,4)], by = "Burrow.x")
colnames(info)[1] <- "Burrow.y"
clustBur2 <- full_join(clustBur2, info[ ,c(1,4)], by = "Burrow.y")
clustBur2$group <- ifelse(clustBur2$clusterCut.x == clustBur2$clusterCut.y, "same", "diff")

ggplot(clustBur2, aes(x=dist, fill=group)) + geom_density(alpha=0.3) +
  theme_classic()
ggsave("Plots/Fig2c.pdf", width=8, height=6, units="cm", dpi=1500, useDingbats=FALSE)

clustBur3 <- clustBur2
clustBur3$pair <- paste(clustBur3$Burrow.x, clustBur3$Burrow.y, sep = ".")
clustBur3 <- clustBur3[ ,c(11,10,7)]
clustBur3_wide <- clustBur3 %>% spread(key="group", value="dist")
ks.test(clustBur3_wide$diff, clustBur3_wide$same) # D = 0.96087, p-value < 2.2e-16