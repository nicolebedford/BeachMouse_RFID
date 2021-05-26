# Nicole Bedford, Hoekstra Lab 2021

## Another approach for visualizing mouse/burrow visitation: 
## https://cran.r-project.org/web/packages/bipartiteD3/vignettes/bipartiteD3_Intro.html

## Interactive plot doesn't work without old version of r2d3:
#packageurl <- "https://cran.r-project.org/src/contrib/Archive/r2d3/r2d3_0.2.3.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")

library(r2d3)
library(bipartite )
library(purrr) 
library(dplyr) 
library(tidyr) 
library(stringr)
library(tibble)
library(RColorBrewer)
library(bipartiteD3)

setwd("~/Dropbox (Personal)/Fieldwork_Paper_NLB/Good/Code/ForGitHub")

pie <- read.csv("pie.csv")
pie_trim <- pie[ ,c(2,5,7)]
mice <- pie[ ,c(5,9)] %>% distinct()
burrows <- pie[ ,c(2,11)] %>% distinct()

bur <- data.frame(higher = pie_trim$Burrow,
                  lower = pie_trim$Mouse_ID,
                  webID = "Visitation",
                  freq = pie_trim$pingSum)

bipartite::frame2webs(bur) -> BurrOcc
bipartite::plotweb(BurrOcc$Visitation)

fo <- bur[,c(1,2,4)] %>% spread(key = higher, value = freq) %>% column_to_rownames(var = "lower")
fo <- fo %>% replace(is.na(.), 0)
bur_order <- OrderByCrossover(fo)

ManualColours <- c('114'="#508CBB", '110'="#1240AB", '123'="#C50000", '112'="#EB5B87", '113'="#0B61A4", '126'="#67E300",
                   '125'="#009E00", '124'="#FF4900", '122'="#FF0000", '111'="#0D4184", '120'="#4CC390", '127'="#FF9063",
                   '104'="#3914AF", '105'="#CD0074", '119'="#00AF64", '109'="#5F5FC6", '93'="#7BB800", '94'="#BC49BC",
                   '121'="#00CC00", 'uk2'="lightgrey", '103'="#FF7400", '115'="#E40045", '107'="#9E0059", '106'="#2B0E88",
                   '129'="#9CEA5B", '118'="#007676", '95'="#C55900", '108'="#1B1BB3", '96'="#FF9200", '97'="#FFBC63",
                   '92'="#9FEE00", '116'="#009999")

all <- bipartite_D3(BurrOcc, SortPrimary = bur_order[[1]], SortSecondary = bur_order[[2]],
             PrimaryLab = "Mouse", SecondaryLab = "Burrow",
             colouroption = "manual",
             NamedColourVector = ManualColours, 
             ColourBy = 1)

all
save_d3_html(all, "Visitations_allMice_allBurs.html")

## Try z-scoring pingSum (number of pings per mouse per burrow)
pie_trim <- pie_trim %>% mutate(zscore = (pingSum - mean(pingSum))/sd(pingSum) + 1) # equivalent to using the function scale() and adding 1

bur <- data.frame(higher = pie_trim$Burrow,
                  lower = pie_trim$Mouse_ID,
                  webID = "Visitation",
                  freq = pie_trim$zscore)

bipartite::frame2webs(bur) -> BurrOcc
bipartite::plotweb(BurrOcc$Visitation)
fo <- bur[,c(1,2,4)] %>% spread(key = higher, value = freq) %>% column_to_rownames(var = "lower")
fo <- fo %>% replace(is.na(.), 0)
bur_order <- OrderByCrossover(fo)

all_zscore <- bipartite_D3(BurrOcc, SortPrimary = bur_order[[1]], SortSecondary = bur_order[[2]],
             PrimaryLab = "Mouse", SecondaryLab = "Burrow",
             colouroption = "manual",
             NamedColourVector = ManualColours, 
             ColourBy = 1) ## Looks much better / less sparse!

all_zscore
save_d3_html(all_zscore, "Visitations_allMice_allBurs_zscore.html")

bipartite_D3(BurrOcc, SortPrimary = bur_order[[1]], SortSecondary = bur_order[[2]],
             PrimaryLab = "Mouse", SecondaryLab = "Burrow",
             colouroption = "manual",
             NamedColourVector = ManualColours, 
             ColourBy = 1,
             Orientation = "horizontal", IncludePerc = FALSE, 
             IndivFigSize = c(175, 600))

## Try a Sankey or Chord Diagram: https://www.data-to-viz.com/graph/sankey.html
library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(chorddiag)
library(networkD3)

nodes <- data.frame(name=c(as.character(pie_trim$Mouse_ID), as.character(pie_trim$Burrow)) %>% unique())

pie_trim$Mouse_ID_num = match(pie_trim$Mouse_ID, nodes$name) -1
pie_trim$Burrow_num = match(pie_trim$Burrow, nodes$name) -1

## Using raw pingSum
sankeyNetwork(Links = pie_trim, Nodes = nodes,
              Source = "Mouse_ID_num", Target = "Burrow_num",
              Value = "pingSum", NodeID = "name")

## Using z-scored pingSum
sankeyNetwork(Links = pie_trim, Nodes = nodes,
              Source = "Mouse_ID_num", Target = "Burrow_num",
              Value = "zscore", NodeID = "name")

node_colors <- 'd3.scaleOrdinal() .domain(["114", "110", "123", "112", "113", "126", "125", "124", 
                                  "122", "111", "120", "127", "104", "105", "119", "109",
                                  "93", "94", "121", "uk2", "103", "115", "107", "106",
                                  "129", "118", "95", "108", "96","97", "92", "116",
                                  "D9a", "A9",
                                  "C49", 
                                  "I31", "H31a",
                                  "G60", "F62", "B59a", "H58b", "D57", "B59b", "H58a",
                                  "A67", "A66a", "F70a", "E70b", "F70b", "E69", "C68", "E70a",
                                  "B22", "B21", "F23a",
                                  "H41b",
                                  "G6a", "G4b", "G4a", "G2", "G6b", "H2",
                                  "D9b", "A2a", "A2b", "A2c",
                                  "I17", "F16", "H15c", "H15a", "C19", "E16"])
.range(["#508CBB","#1240AB","#C50000","#EB5B87","#0B61A4","#67E300","#009E00","#FF4900",
"#FF0000","#0D4184","#4CC390","#FF9063","#3914AF","#CD0074","#00AF64","#5F5FC6", 
"#7BB800","#BC49BC","#00CC00","lightgrey","#FF7400","#E40045","#9E0059","#2B0E88",
"#9CEA5B","#007676","#C55900","#1B1BB3","#FF9200","#FFBC63","#9FEE00","#009999",
"lightgrey", "lightgrey",
"#508CBB",
"#EB5B87", "#EB5B87",
"#3914AF", "#3914AF", "#3914AF", "#3914AF", "#3914AF", "#3914AF", "#3914AF",
"#FF4900", "#FF4900", "#FF4900", "#FF4900", "#FF4900", "#FF4900", "#FF4900", "#FF4900",
"#FFBC63", "#FFBC63", "#FFBC63", 
"#C50000",
"#009999", "#009999", "#009999", "#009999", "#009999", "#009999",
"#BC49BC", "#BC49BC", "#BC49BC", "#BC49BC",
"#C55900", "#C55900", "#C55900", "#C55900", "#C55900", "#C55900"])'

sn <- sankeyNetwork(Links = pie_trim, Nodes = nodes,
              Source = "Mouse_ID_num", Target = "Burrow_num",
              Value = "zscore", NodeID = "name",
              colourScale = node_colors,
              width = 300, height = 1000)

save_d3_html(sn, file = "sn.html", background = "transparent")

#############################################
#### Bipartite plot for home burrow mice ####
#############################################

pie <- pie %>% filter(Mouse_ID %in% c('95', '118', '129', '106', '92', '116',
                                       '124', '125', '105', '105', '104', '97', '108'))

ManualColours_mice <- c('104'="#3914AF", '105'="#CD0074", '124'="#FF4900", '125'="#009E00", '116'="#009999", '92'="#9FEE00",
                       '118'="#007676", '129'="#9CEA5B", '106'="#2B0E88", '95'="#C55900", '108'="#1B1BB3", '97'="#FFBC63")
