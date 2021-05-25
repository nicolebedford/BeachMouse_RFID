## Code for generating Overlap Index Data (Real and Shuffled)
## Read-in the data:
## Run the whole script from RFID_Master_dataframe.R to get the "nov" and "may" data frames:
## Clear the rest of the workspace by switching List to Grid, select all, remove nov and may, then Clear Workspace

library(plyr)
library(dplyr)
library(ggplot2)
library(data.table) # setorder
library(gtools) # permutations
library(evobiR) # SlidingWindow
library(gdata) # reorder.factor
library(forcats) # fct_rev
library(viridis)
library(colormap)
library(grid) # rasterGrob
library(multcomp)
library(tibble)

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

# Nested for-loop for removing ping trains for each individual:
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

#### 2a. Real data distributions of mouse and burrow pings #### ---------------------------------------------------------------------------------------------------------------------------------------------
df <- trim %>% 
  dplyr::select(Burrow, Mouse_ID, Sex, Timestamp)
df$Mouse_ID <- as.factor(df$Mouse_ID)
df$Timestamp <- as.POSIXct(df$Timestamp)
df <- df %>%
  dplyr::select(Burrow, Mouse_ID, Timestamp)
Mice <- df %>% group_by(Mouse_ID) %>% summarise(n=n()) # how many pings per mouse?
hist(Mice$n) # non-normal
Burrow <- df %>% group_by(Burrow) %>% summarise(n=n()) # how many pings per burrow?
hist(Burrow$n) # non-normal

# #### 2b. Optional Shuffle #### ---------------------------------------------------------------------------------------------------------------------------------------------
# df$Mouse_ID <- base::sample(Mice$Mouse_ID, size=nrow(df), replace=TRUE) # randomly draw mouse labels from Mouse_ID list
# df$Burrow <- base::sample(Burrow$Burrow, size=nrow(df), replace=TRUE) # randomly draw burrow labels from Burrow list
# Mice <- df %>% group_by(Mouse_ID) %>% summarise(n=n()) # how many pings per mouse?
# hist(Mice$n) # normal!
# Burrow <- df %>% group_by(Burrow) %>% summarise(n=n()) # how many pings per burrow?
# hist(Burrow$n) # normal!
# #### End of shuffle section #### ---------------------------------------------------------------------------------------------------------------------------------------------

## Make new column "MouseBur"
df$MouseBur <- factor(paste(df$Mouse_ID, df$Burrow, sep="."))
df$ts_num <- as.numeric(df$Timestamp) # convert to seconds 

#### 3. Sliding Window Analysis: window = w ; step = w/2 #### ----------------------------------------------------------------------------------------------------------------
n <- max(df$ts_num) - min(df$ts_num) + 1 # what is the range of time stamps?
all_ts <- data.frame(matrix(0, nrow=n))
all_ts$ts_num <- seq(from = min(df$ts_num), to = max(df$ts_num), by = 1)
colnames(all_ts)[1] <- "count"
mb <- as.character(unique(df$MouseBur))

## Define the window size: try the following time bins: 15m, 30m, 45m, 1h, 1.5h, 2h, 3h, 4h, 6h, 8h, 10h, 12h, 14h
w <- 60*60*1.5

countPings <- function(i) {
  s <- df %>% filter(MouseBur %in% i)
  fo <- full_join(s, all_ts, by="ts_num")
  fo$count <- ifelse(is.na(fo$Mouse_ID), 0, 1)
  fo <- setorder(fo, ts_num)
  sw <- SlidingWindow(FUN=sum, data=fo$count, window = w, step = w/2)
  result <- data.frame(sw)
  result$bin <- seq(1, nrow(result), 1)
  result$window <- w
  result$g <- i
  return(data.frame(result))
}
good <- lapply(mb, countPings) %>% bind_rows()

## Make good dataframe of no. of pings per placeTime bin per animal
x <- data.frame(do.call('rbind', strsplit(as.character(good$g), '.', fixed=TRUE)))
new <- cbind(x, good[1:3])
colnames(new) <- c("Mouse_ID", "Burrow", "count", "bin", "window")
new$placeTime <- as.factor(paste(new$Burrow, new$bin, sep="."))
plot(new$bin, new$count) # you should see one clump per night

## Compute function for all pairwise combinations of mouse IDs
choose(length(unique(df$Mouse_ID)), 2) # 496 possible dyads 
mice <- as.character(unique(new$Mouse_ID))
perms <- as.data.frame(permutations(n=length(unique(df$Mouse_ID)), r=2, v=mice)) # Get all permutations of mice (order does matter)
colnames(perms) <- c("mouse1", "mouse2")
perms$pair <- as.factor(paste(perms$mouse1, perms$mouse2, sep="."))
all_pt_pair <- expand.grid(placeTime = unique(new$placeTime), pair = unique(perms$pair)) # all possible place-times and mouse pairs
pt <- data.frame(do.call('rbind', strsplit(as.character(all_pt_pair$placeTime), '.', fixed=TRUE)))
colnames(pt) <- c("Burrow", "bin")
pair <- data.frame(do.call('rbind', strsplit(as.character(all_pt_pair$pair), '.', fixed=TRUE)))
colnames(pair) <- c("mouse1", "mouse2")
all_pt_pair <- cbind(all_pt_pair, pt, pair) # make dataframe of all possible place-times and mouse pairs, make same columns as new 
all_pt_pair$Mouse_ID <- all_pt_pair$mouse1
all_pt_pair$window <- w

## Add ping count data for mouse1
new$mouse1 <- new$Mouse_ID
new$bin <- as.factor(new$bin) # convert bin to factor so we can join with all_pt_pair
M1 <- full_join(new, all_pt_pair)
M1 <- M1 %>% 
  dplyr::select(mouse1, mouse2, pair, Burrow, window, bin, count, placeTime)
colnames(M1) <- c("mouse1", "mouse2", "pair", "Burrow", "window", "bin", "count1", "placeTime")
M1[is.na(M1)] <- 0 # Replace NAs with zeros

## Add ping count data for mouse2
new$mouse2 <- new$Mouse_ID
M2 <- full_join(M1, new[ ,c(3,6,8)]) # merge with count, placeTime, and mouse2
M2[is.na(M2)] <- 0 # Replace NAs with zeros
colnames(M2)[9] <- "count2"

## Get total number pings per mouse (over entire study, using sliding window method --> will be inflated over actual ping count)
pingTot <- M1 %>% group_by(mouse1) %>% summarise(pingTot = sum(count1))
M2 <- full_join(M2, pingTot, by="mouse1") # add ping total for mouse1
colnames(pingTot)[1] <- "mouse2"
M2 <- full_join(M2, pingTot, by="mouse2") # add ping total for mouse2
colnames(M2)[10:11] <- c("pingTot1", "pingTot2")

#### 4a. Calculate "Overlap Index" #### ------------------------------------------------------------------------------------------------------------------------------------
## OI = for each mouse i, what portion of its total pings occur in same placeTime bin as mouse j?
M2$normCount <- M2$count1 / M2$pingTot1 # First normalize counts per placeTime
total <- M2 %>% group_by(mouse1) %>% summarise(total=sum(normCount)) # each animal's normCount should add up to one!
nz <- subset(M2, M2$count1 > 0 & M2$count2 > 0) # look only at rows with overlap between mice
pair_count <- nz %>% group_by(pair) %>% summarise(count = length(pair)) # how many instances of overlap per pair?

## Calculate the Overlap Index:
OI <- nz %>% group_by(mouse1, mouse2, pair) %>% summarise(OI = sum(normCount)) 
OI <- full_join(perms, OI) # this dataframe should be the same length as perms
OI[is.na(OI)] <- 0 # Replace NAs with zeros

## Add sex information
info <- copy[ ,c(5,17)] # pull out sex
info <- info[!duplicated(info),] # remove duplicates
colnames(info)[1] <- "mouse1"
OI <- full_join(OI, info, by="mouse1") 
colnames(info)[1] <- "mouse2"
OI <- full_join(OI, info, by="mouse2") 
colnames(OI)[5:6] <- c("sex1", "sex2")
OI$type <- paste(OI$sex1, OI$sex2, sep="")
OI$type <- sub("FM", "MF", OI$type)
OI$type <- factor(OI$type, levels=c("FF", "MM", "MF"))
OI$window <- w
#write.csv(OI[,c(1,2,4,8)], "Nov17_OI_5400.csv") # Save a different file for each time bin!

#### 4b. Which time bin is most appropriate? #### --------------------------------------------------------------------------------------------------------------------------
## Calculate Network Load at different OI thresholds for different timebins:
load_data <- function(path) {
  files <- dir(path, pattern = '\\.csv', full.names = TRUE)
  tables <- lapply(files, read.csv)
  do.call(rbind, tables)
}
all <- load_data("windowSize_OI") # load all OI data for different windows
all[1] <- NULL
all[all == 0] <- NA # replace zeroes with NA

## Network Load (i.e. fraction of observed links (M) in the network over all possible pair combinations (from Psorakis et al. 2012)
numPairs <- choose(length(unique(all$mouse1)), 2)
all_trim <- all[complete.cases(all), ] # remove rows containing NAs
linkCount <- all_trim %>% dplyr::count(window) # how many non-zero OI values are observed?
plot(linkCount$window, linkCount$n)
linkCount$timeBin <- linkCount$window / (60*60) # convert to hours
linkCount$NetworkLoad <- linkCount$n / (numPairs*2)
all_trim$timeBin <- all_trim$window / (60*60)
hist(all_trim$OI)

## Find the elbow of the zero threshold curve: https://stackoverflow.com/questions/41518870/finding-the-elbow-knee-in-a-curve
nl <- all_trim %>% dplyr::count(timeBin)
nl$NetworkLoad <- nl$n / (numPairs)

get.elbow.points.indices <- function(x, y, threshold) {
  d1 <- diff(y) / diff(x) # first derivative
  d2 <- diff(d1) / diff(x[-1]) # second derivative
  indices <- which(abs(d2) > threshold)
  return(indices)
}

## First approximate the function (i.e. the curve), since we have only a few points:
x = nl$timeBin
y = nl$NetworkLoad
ap <- approx(x, y, n=1000, yleft=min(y), yright=max(y))
x <- ap$x
y <- ap$y
indices <- get.elbow.points.indices(x, y, 2) # threshold for jump
x[indices] # 0.483984 1.488739
plot(x, y, pch=19)
points(x[indices], y[indices], pch=19, col='red') # which points have the maximum second derivative? 0.5 and 1.5 --> use 1.5 (90 min)

ggplot(data=nl, aes(x=timeBin, y=NetworkLoad/2)) +
  geom_line(size=0.2) +
  geom_point(size=0.2) +
  xlab("Time bin") + ylab("Network Load") +
  scale_x_continuous(breaks=seq(0,14,2)) +
  scale_y_continuous(limits = c(0.09, 0.25)) +
  theme_classic()

## Make a histogram of all timestamps:
df$Time <- format(as.POSIXct(strptime(df$Timestamp, "%Y-%m-%d %H:%M:%S", tz="")), format = "%H:%M:%S")
df$Time <- as.POSIXct(df$Time, format="%H:%M:%S")
hist(df$Time, breaks=24)
df$Time_num <- as.numeric(df$Time)
df$Time_num <- df$Time_num - min(df$Time_num)
hist(df$Time_num, breaks=24)
range(df$Time_num) # 0 to 86349 (Nov 2017)
df$Time_num <- ifelse(df$Time_num > 43200, df$Time_num - 86400, df$Time_num) # Noon 12:00:00 = 43200 ; Midnight 00:00:00 = 86400
hist(df$Time_num, breaks=24)

df$Time_num <- df$Time_num + 40 # Nov17: midnight is 40 seconds ahead --> fix

24*(max(df$Time_num) - min(df$Time_num)) / 86400 # 13.45 hours (Nov 2017)
colnames(info)[1] <- "Mouse_ID"
df <- inner_join(df, info)
df$Sex <- factor(df$Sex, levels = c("M", "F", "U"))

## Sunrise/Sunset on Nov. 25th 2017 = 6:20AM - 4:46PM
sunrise17 = 0 + 6*60*60 + 20*60 # midnight plus 6h20m
sunset17 = -43200 + 4*60*60 + 46*60 # noon plus 4h46m

ggplot(data = df, aes(x = Time_num)) +
  geom_histogram(colour="black", fill="grey", bins=24) +
  geom_vline(xintercept = sunrise17) +
  geom_vline(xintercept = sunset17) +
  scale_x_continuous(breaks=seq(-43200, 43200, 7200), limits=c(-43200, 43200)) + # 2h increments
  theme_classic()
ggsave("Plots/S4a.pdf", width=18, height=12, units="cm", dpi=1500, useDingbats=FALSE)

#### 5. Make a heatmap figure using Overlap Index #### ---------------------------------------------------------------------------------------------------------------------
## Make rows corresponding to diagonal in the heatmap
colnames(info)[1] <- "mouse2"
info$mouse2 <- factor(info$mouse2)
dia <- as.data.frame(matrix(NA, ncol=ncol(OI), nrow=length(unique(OI$mouse1))))
dia[1] <- factor(unique(OI$mouse1))
dia[2:3] <- info[order(info$mouse2), ]
dia$V4 <- dia$V3
dia$V5 <- factor(paste(dia$V3, dia$V4, sep=""))
dia$V6 <- factor(paste(dia$V1, dia$V2, sep="."))
dia$V7 <- NA
colnames(dia) <- c("mouse1", "mouse2", "sex1", "sex2", "type", "pair", "OI", "window")
dia <- dia %>% 
  dplyr::select(mouse1, mouse2, pair, OI, sex1, sex2, type, window) # put diagonal in same order as OI dataframe
OI <- rbind(OI, dia)
OI$mouse1 <- factor(OI$mouse1)
OI$mouse2 <- factor(OI$mouse2)
str(OI)
#write.csv(OI, "OI_1.5h_real.csv")
#write.csv(OI, "OI_1.5h_shuff.csv")

## Turn the data into a matrix: Repeat with both real and shuffled data
#OI <- read.csv("OI_1.5h_real.csv", header=T, sep=",")
#OI <- read.csv("OI_1.5h_shuff.csv", header=T, sep=",")
OI$mouse1 <- as.factor(OI$mouse1)
OI$mouse2 <- as.factor(OI$mouse2)

dataOrdered <- setorder(OI, mouse1, mouse2)
dataMatrix <- as.data.frame(matrix(dataOrdered$OI, nrow=length(unique(dataOrdered$mouse1)), ncol=length(unique(dataOrdered$mouse2))))
colnames(dataMatrix) <- unique(dataOrdered$mouse2)            
rownames(dataMatrix) <- unique(dataOrdered$mouse1)   

## Order the data by cluster: Apply realOrd to shuffle data as well!
realOrd <- hclust(dist(dataMatrix, method = "euclidean"), method = "ward.D")$order # find euclidean distance between all points, then cluster hierarchically using ward.D
OI$mouse1 <- reorder.factor(OI$mouse1, new.order=realOrd) # apply cluster order to mouse1
OI$mouse2 <- reorder.factor(OI$mouse2, new.order=realOrd) # apply cluster order to mouse2
OI$mouse2 <- fct_rev(OI$mouse2) # reverse factor order for mouse2 --> diagonal now runs top to bottom

#### 6. Make Figures (Overlap Index clustered heatmap) #### ----------------------------------------------------------------------------------------------------------------
OI_nz <- OI %>% filter(OI > 0) # only non-zero values of OI
max(OI_nz$OI) # 0.02419355
hist(OI_nz$OI)

ggplot(OI, aes(x=mouse2, y=mouse1, fill=OI)) +
  geom_tile() +
  scale_fill_viridis(limits=c(0, 0.02419355)) +
  ggtitle("Nov17_real_1.5h") +
  theme_classic() + theme(aspect.ratio=1)
ggsave("Plots/Fig3a.pdf", width=10, height=10, units="cm", dpi=1500, useDingbats=FALSE)

## Make a colour gradient for the histogram background palette
grad <- colormap(colormaps$viridis, nshades = 100) # Generate color codes
gradHorz <- t(grad) # Transpose to make a horizontal gradient
g <- rasterGrob(gradHorz, height=unit(1,"npc"),
                width=unit(1,"npc"), interpolate=TRUE)

ggplot(OI_nz, aes(x=OI)) +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + # add gradient background
  geom_line(stat="density", colour="white", size=2) +
  expand_limits(y=0) +
  xlab("Overlap Index") +
  xlim(0, 0.02419355) +
  geom_vline(xintercept=mean(OI_nz$OI, na.rm=T), colour="white", linetype=2, size=1) +
  theme(aspect.ratio=0.5)
ggsave("Plots/Fig3b.pdf", width=10, height=10, units = "cm")

## Network Load, Max, Mean, Median (non-zero) for each year: real, shuff
nrow(OI_nz) # 148/992 = 0.1491935 (real), 992/992 = 1.00 (shuff)
max(OI_nz$OI) # 0.02419355 (real), 0.004198447/0.004411977 (shuff) 
mean(OI_nz$OI, na.rm=T) # 0.005361158 (real), 0.002449568/0.002469737 (shuff) 
