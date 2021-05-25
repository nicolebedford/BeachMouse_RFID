#### Nicole Bedford, Hoekstra Lab 2021
#### Code for Florida Fieldwork Paper Analyses
#### Explore Burrow details

## Read-in the data:
## Run the whole script from RFID_Master_dataframe.R to get the "nov" and "may" data frames:
## Clear the rest of the workspace by switching List to Grid, select all, remove nov and may, then Clear Workspace

## Load libraries:
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(Rmisc)
library(agricolae)
library(lubridate)

## Set working directory:
setwd("~/Dropbox (Personal)/Fieldwork_Paper_NLB/Good/Code/ForGitHub")

#### Start with November 2017 Dataset #### --------------------------------------------------------------------------------------------------------------------------------
## Read-in the RFID data: (run RFID_Master_dataframe.R)
copy <- nov # make a copy of the November dataframe

#### 1a. Remove "ping trains" #### ------------------------------------------------------------------------------------------------------------------------------------------
## If pings occur at 1 Hz, keep first and last in the "train"
copy <- copy %>% filter(!is.na(Night)) # keep only Nights 1 to 11 (no grid traps in place)
copy <- copy %>% filter(!Mouse_ID %in% c("uk1", "uk3", "uk4")) # these are "starter" RFID tags --> not real mice!
pingSum <- as.data.frame(table(copy$Mouse_ID)) # number of pings per animal before train filtering
preFilt_mouseSummary <- copy %>% group_by(Mouse_ID) %>% 
  dplyr::summarise(numNight=n_distinct(Night), numBur=n_distinct(Burrow), pingSum=n())
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
  p <- copy %>% filter(Night == i) %>% dplyr::summarise(burs = n_distinct(Burrow), mice = n_distinct(Mouse_ID), pings = n()) 
  p$Night <- i  
  resultlist[[i]] = p 
}
result = do.call(rbind, resultlist)
copy %>% dplyr::summarise(burs = n_distinct(Burrow), mice = n_distinct(Mouse_ID), pings = n()) # Summary for all 11 Nights

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

#### 1b. Do long burrows have more juvenile pings per night? #### ------------------------------------------------------------------------------------------------------------------------------------------

j <- trim %>% filter(Age == "Juvenile") # 3044
j <- j %>% group_by(Burrow, Length, Cover, Night) %>% dplyr::summarise(juv_n=n())
a <- trim %>% filter(Age == "Adult") # 4372
a <- a %>% group_by(Burrow, Length, Cover, Night) %>% dplyr::summarise(ad_n=n())
all <- full_join(j,a)
all[ ,c(5:6)][is.na(all[ ,c(5:6)])] <- 0
all$juv_portion <- all$juv_n / (all$juv_n + all$ad_n)
all <- all %>% group_by(Burrow, Length, Cover) %>% dplyr::summarise(mean_juvPor = mean(juv_portion))

ggplot(all, aes(x=Length, y=mean_juvPor)) +
  geom_point(size=0.5) +
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE) +
  theme_classic() +
  xlim(c(0,100)) +
  xlab("burrow length") + ylab("mean portion juvenile pings")
ggsave("Plots/Fig1f.pdf", width=8, height=6, units="cm", dpi=1500, useDingbats=FALSE)

hist((all$mean_juvPor)) # not normally distributed :(
summary(lm(mean_juvPor ~ Length, data = all)) # Adjusted R-squared:  0.1971, p-value: 0.006352
cor.test(all$Length, all$mean_juvPor, use="complete.obs", method="spearman") # p-value = 0.005539

#### 1c. Make dataframe of Burrow Specs #### ------------------------------------------------------------------------------------------------------------------------------------------
spec <- copy %>% group_by(Burrow, X, Y, Slope, Heading, Elev, Lat, Lon, Length, Cover) %>% 
  dplyr::summarise(pingCount = n(), mouseCount = n_distinct(Mouse_ID), nightCount = n_distinct(Night))
spec <- as.data.frame(spec)

## Read-in notes data: is burrow collapsed or crabby-looking?
notes <- read.csv("BurrowLength_notes.csv", header=T, sep=",")
spec <- join(spec, notes)

## Make dataframe of per night Burrow summaries:
perNight <- copy %>% group_by(Burrow, Night) %>% dplyr::summarise(pings = n(), mice = n_distinct(Mouse_ID))
perNight <- perNight %>% group_by(Burrow) %>% dplyr::summarise(nightCount = n(), pings_per_night = mean(pings), mice_per_night = mean(mice))

spec <- join(spec, perNight)

#### 2. Explore Correlations #### ------------------------------------------------------------------------------------------------------------------------------------------
## Remove collapsed burrows:
# spec <- spec %>% filter(!Notes %in% "collapsed") # doesn't impact conclusions below!

## Test for associations between burrow length and aspects of mouse behavior: none!
cor.test(spec$Length, spec$pingCount, use="complete.obs", method="spearman") # p-value = 0.2996
cor.test(spec$Length, spec$mouseCount, use="complete.obs", method="spearman") # p-value = 0.7104
cor.test(spec$Length, spec$pings_per_night, use="complete.obs", method="spearman") # p-value = 0.3311
cor.test(spec$Length, spec$mice_per_night, use="complete.obs", method="spearman") # p-value = 0.5645
cor.test(spec$Length, spec$nightCount, use="complete.obs", method="spearman") # p-value = 0.5621

## Test for associations between burrow length and other specs: none!
cor.test(spec$Length, spec$Slope, use="complete.obs", method="spearman") # p-value = 0.4002
cor.test(spec$Length, spec$Heading, use="complete.obs", method="spearman") # p-value = 0.5179
cor.test(spec$Length, spec$Elev, use="complete.obs", method="spearman") # p-value = 0.07851
cor.test(spec$Length, spec$Lat, use="complete.obs", method="spearman") # p-value = 0.9143
cor.test(spec$Length, spec$Lon, use="complete.obs", method="spearman") # p-value = 0.6892
cor.test(spec$Length, spec$Cover, use="complete.obs", method="spearman") # p-value = 0.2675

## Read-in home burrow information:
## days = no. of days on which 1 or more mice slept in the burrow, NA = no daytime occupancy inferred
HB <- read.table("HomeBurrows.csv", header = T, sep=",")
HB[1] <- NULL
spec <- join(spec, HB)
spec$homeBurrow <- ifelse(is.na(spec$days) == TRUE, "other", "home")

## Do home burrows have different attributes than other burrows? 
wilcox.test(Length ~ homeBurrow, data=spec) # p-value = 0.6598
wilcox.test(Slope ~ homeBurrow, data=spec) # p-value = 0.2276
wilcox.test(Heading ~ homeBurrow, data=spec) # p-value = 0.07421
wilcox.test(Elev ~ homeBurrow, data=spec) # p-value = 0.1153
wilcox.test(Lat ~ homeBurrow, data=spec) # p-value = 1
wilcox.test(Lon ~ homeBurrow, data=spec) # p-value = 0.9006
wilcox.test(Cover ~ homeBurrow, data=spec) # p-value = 0.01157 *

## Percent Cover relationships:
cor.test(spec$Cover, spec$pingCount, use="complete.obs", method="spearman") # p-value = 0.01982 *
cor.test(spec$Cover, spec$mouseCount, use="complete.obs", method="spearman") # p-value = 0.4301
cor.test(spec$Cover, spec$pings_per_night, use="complete.obs", method="spearman") # p-value = 0.0123 *
cor.test(spec$Cover, spec$mice_per_night, use="complete.obs", method="spearman") # p-value = 0.5801
cor.test(spec$Cover, spec$nightCount, use="complete.obs", method="spearman") # p-value = 0.1947

#### 3. Plot burrow spec histograms #### ------------------------------------------------------------------------------------------------------------------------------------------
par(mfrow=c(2,2))
hist(spec$Length, breaks=15, col="blue", xlab="entrance tunnel length (cm)", main=NULL)
hist(spec$Slope, breaks=15, col="orange", xlab="burrow site slope (degrees)", main=NULL)
hist(spec$Cover, breaks=15, col="forestgreen", xlab="percent cover (%)", main=NULL)
hist(spec$Elev, breaks=10, col="turquoise", xlab="burrow site elevation (feet)", main=NULL)

ggplot(spec, aes(x=Length)) +
  geom_histogram(binwidth = 5, colour="black", fill="grey") +
  scale_x_continuous(breaks=seq(0,100,25), limits = c(0,100)) + 
  theme_classic()
ggsave("Plots/FigS2a.pdf", width=8, height=6, units="cm", dpi=1500, useDingbats=FALSE)

ggplot(spec, aes(x=Cover)) +
  geom_histogram(binwidth = 5, colour="black", fill="grey") +
  scale_x_continuous(breaks=seq(0,100,20)) + 
  theme_classic()
ggsave("Plots/FigS2b.pdf", width=8, height=6, units="cm", dpi=1500, useDingbats=FALSE)

ggplot(spec, aes(x=Slope)) +
  geom_histogram(binwidth = 2.5, colour="black", fill="grey") +
  scale_x_continuous(breaks=seq(0,50,10)) + 
  theme_classic()
ggsave("Plots/FigS2c.pdf", width=8, height=6, units="cm", dpi=1500, useDingbats=FALSE)

ggplot(spec, aes(x=Elev)) +
  geom_histogram(binwidth = 2.5, colour="black", fill="grey") +
  scale_x_continuous(breaks=seq(0,50,10)) + 
  theme_classic()
ggsave("Plots/FigS2d.pdf", width=8, height=6, units="cm", dpi=1500, useDingbats=FALSE)

ggplot(spec)  +
  geom_histogram(aes(Heading), binwidth = 20, colour="black", fill="grey") +
  scale_x_continuous(breaks=seq(0,360,45), limits = c(0,360), oob = scales::squish) +
  coord_polar() +
  theme(axis.text.x = element_text(size = 16)) +
  theme_classic()
ggsave("Plots/FigS2e.pdf", width=8, height=6, units="cm", dpi=1500, useDingbats=FALSE)

ggplot(spec, aes(x=Cover, y=log2(pings_per_night))) +
  geom_point(size=0.5) +
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE) +
  theme_classic() +
  xlab("Percent Cover") + ylab("mean pings per night")
ggsave("Plots/Fig1e.pdf", width=8, height=6, units="cm", dpi=1500, useDingbats=FALSE)

hist(log2(spec$pings_per_night)) # tranformed data looks better!
summary(lm(log2(pings_per_night) ~ Cover, data = spec)) # Adjusted R-squared:  0.1439, p-value: 0.00993
summary(lm(log2(pings_per_night) ~ Length, data = spec)) # p-value: 0.2962

mean(spec$Length, na.rm = T) # 36.6
mean(spec$Slope, na.rm = T) # 21.6
mean(spec$Elev, na.rm = T) # 19.7
mean(spec$Heading, na.rm = T) # 178
mean(spec$Cover, na.rm = T) # 57.7

#### 4a. HOBO Temperature Data #### ------------------------------------------------------------------------------------------------------------------------------------------
## Artificial Burrow Data from 2015:
nov15 <- read.csv("Nov15_artificial_burrows.csv", header=T, sep=",")
nov15$Date_time <- as.POSIXct(nov15$Date_time) 
decode <- read.csv("Nov15_decode.csv", header = T, sep=",")

ggplot(nov15, aes(x=Date_time, y=Temp, colour = Burrow)) +
  geom_line()

## Remove timepoints when HOBOs were being moved between sites:
## There are 7 sessions (7 unique burrow sites)
s1 <- nov15 %>% filter(Date_time > "0015-11-02 21:00:00", Date_time < "0015-11-04 18:45:00") 
s2 <- nov15 %>% filter(Date_time > "0015-11-04 20:15:00", Date_time < "0015-11-06 10:20:00") # start 1 hour later
s3 <- nov15 %>% filter(Date_time > "0015-11-06 20:25:00", Date_time < "0015-11-08 01:10:00") # start 1.5 hours later
s4 <- nov15 %>% filter(Date_time > "0015-11-08 14:20:00", Date_time < "0015-11-10 18:30:00") # start 1 hour later
s5 <- nov15 %>% filter(Date_time > "0015-11-10 20:00:00", Date_time < "0015-11-12 19:50:00") # start 1 hour later
s6 <- nov15 %>% filter(Date_time > "0015-11-12 21:05:00", Date_time < "0015-11-15 18:05:00") # start 1 hour later
s7 <- nov15 %>% filter(Date_time > "0015-11-15 19:20:00", Date_time < "0015-11-17 08:00:00") # start 1 hour later

## Add session info:
s1$session <- "s1"
s2$session <- "s2"
s3$session <- "s3"
s4$session <- "s4"
s5$session <- "s5"
s6$session <- "s6"
s7$session <- "s7"

## Combine all sessions:
nov15_filt <- rbind(s1, s2, s3, s4, s5, s6, s7)
nov15_filt <- plyr::join(nov15, nov15_filt)
nov15_filt$Temp <- ifelse(is.na(nov15_filt$session) == T, NA, nov15_filt$Temp)

max(nov15_filt$Date_time) - min(nov15_filt$Date_time) # Time difference of 14.40278 days
range(nov15_filt$Temp, na.rm = T) # 7.870 36.824

## Plot the filtered data:
ggplot(nov15_filt %>% filter(Date_time < "0015-11-17 08:00:00"), aes(x=Date_time, y=Temp, colour = Burrow)) +
  geom_line() +
  scale_x_datetime(date_breaks = "1 day") +
  scale_y_continuous(breaks = seq(0, 40, 5)) +  
  theme_classic()
ggsave("Plots/FigS3a.pdf", width=14.4, height=6, units="in", dpi=1500, useDingbats=FALSE)

## Align Ambient with other timestamps (every 15 min, vs. every 10 min):
nov15_filt <- nov15_filt[ ,-c(3,5)]
nov15_filt <- nov15_filt %>% distinct()
nov15_wide <- spread(nov15_filt, key = "Burrow", value = "Temp") # convert to wide format
a <- nov15_wide[ ,c(1,6)]
a <- a[complete.cases(a), ]
a.15 <- a %>% filter(minute(Date_time) == 15) # change all 15:00 timepoints to 10:00
a.15$Date_time <- a.15$Date_time - minutes(5) # should now read 10:00
a.45 <- a %>% filter(minute(Date_time) == 45) # change all 45:00 timepoints to 50:00
a.45$Date_time <- a.45$Date_time + minutes(5) # should now read 50:00
amb <- a %>% filter(minute(Date_time) == 00 | minute(Date_time) == 30)
amb <- rbind(amb, a.15, a.45)
nov15_wide <- full_join(nov15_wide[ ,c(1:5)], amb)  
nov15_wide <- nov15_wide[complete.cases(nov15_wide), ]

## Calculate difference from Ambient:
nov15_wide$diff10 <- abs(nov15_wide$`10` - nov15_wide$A)
nov15_wide$diff40 <- abs(nov15_wide$`40` - nov15_wide$A)
nov15_wide$diff60 <- abs(nov15_wide$`60` - nov15_wide$A)
nov15_wide$diff80 <- abs(nov15_wide$`80` - nov15_wide$A)
nov15_long <- gather(nov15_wide[ ,c(7:10)], "Burrow", "diff_from_Ambient")

## Run the ANOVA:
m <- aov(diff_from_Ambient ~ Burrow, data = nov15_long)
summary(m)
hsd <- HSD.test(m, trt = "Burrow")
hsd$groups

## Make the plot:
summary <- summarySE(data = nov15_long, measurevar = "diff_from_Ambient", groupvars = "Burrow")
summary$upper <- summary$diff_from_Ambient + summary$ci
summary$lower <- summary$diff_from_Ambient - summary$ci

ggplot() +
  geom_pointrange(data=summary, mapping=aes(x=Burrow, y=diff_from_Ambient, ymin=lower, ymax=upper)) +
  theme_classic() +
  scale_y_continuous(breaks=seq(0,5,1), limits=c(0,5)) 
ggsave("Plots/FigS3b.pdf", width=6, height=6, units="in", dpi=1500, useDingbats=FALSE)

#### 4b. HOBO Temperature Data #### ------------------------------------------------------------------------------------------------------------------------------------------
## Natural Burrow Data from 2017:
nov17 <- read.csv("Nov17_natural_burrows.csv", header=T, sep=",")
nov17$Date_time <- as.POSIXct(paste(nov17$Date, nov17$Time, sep = " "))
nov17$HOBO <- factor(nov17$HOBO)
decode <- read.csv("Nov17_decode.csv", header = T, sep=",")

ggplot(nov17, aes(x=Date_time, y=Temp, colour = HOBO)) +
  geom_line()

## There are 3 sessions:
s1 <- nov17 %>% filter(Date_time > "0017-11-04 18:00:00", Date_time < "0017-11-10 08:00:00") # good
s2 <- nov17 %>% filter(Date_time > "0017-11-10 17:00:00", Date_time < "0017-11-21 15:00:00") # start 1 hour later
s3 <- nov17 %>% filter(Date_time > "0017-11-24 11:00:00", Date_time < "0017-11-29 12:00:00") # start 1 hour later, end much earlier

## Add session info:
s1$session <- "s1"
s2$session <- "s2"
s3$session <- "s3"

## Combine all sessions:
nov17_filt <- rbind(s1, s2, s3)
nov17_filt <- plyr::join(nov17, nov17_filt)
nov17_filt$Temp <- ifelse(is.na(nov17_filt$session) == T, NA, nov17_filt$Temp)
nov17_filt <- plyr::join(nov17_filt, decode[ ,c(1,8,9)]) # HOBO, length, session
nov17_filt.noNA <- nov17_filt[complete.cases(nov17_filt), ] 
max(nov17_filt.noNA$Date_time) - min(nov17_filt.noNA$Date_time) # Time difference of 24.72917 days
range(nov17_filt.noNA$Temp) # 2.717 37.096

## Plot the filtered data:
ggplot(nov17_filt.noNA %>% filter(Date_time < "0017-11-29 12:00:00"), aes(x=Date_time, y=Temp, colour = HOBO)) +
  geom_line() +
  scale_x_datetime(date_breaks = "2 days") +
  scale_y_continuous(breaks = seq(0, 40, 5)) +  
  theme_classic()
ggsave("Plots/FigS3d.pdf", width=24.7, height=6, units="in", dpi=1500, useDingbats=FALSE)

#### 4c. HOBO Temperature Data #### ------------------------------------------------------------------------------------------------------------------------------------------
## Natural Burrow Data from May 2016:
library(plotly)
may16 <- read.csv("May16_natural_burrows.csv", header=T, sep=",")
may16$Date_time <- as.POSIXct(paste(may16$Date, may16$Time, sep = " "))
may16$HOBO <- factor(may16$HOBO)
decode <- read.csv("May16_decode.csv", header = T, sep=",")

plot <- ggplot(may16, aes(x=Date_time, y=Temp, colour = HOBO)) +
  geom_line()
ggplotly(plot) # make interactive plot

## There are 2 sessions:
s1 <- may16 %>% filter(Date_time > "0016-05-04 20:00:00", Date_time < "0016-05-11 19:00:00") 
s2 <- may16 %>% filter(Date_time > "0016-05-12 19:00:00", Date_time < "0016-05-17 08:30:00") 

## Add session info:
s1$session <- "s1"
s2$session <- "s2"

## Combine all sessions:
may16_filt <- rbind(s1, s2)
may16_filt <- plyr::join(may16, may16_filt)
may16_filt$Temp <- ifelse(is.na(may16_filt$session) == T, NA, may16_filt$Temp)
may16_filt <- plyr::join(may16_filt, decode[ ,c(1,2,5)]) # HOBO, length, session
may16_filt.noNA <- may16_filt[complete.cases(may16_filt), ] 
max(may16_filt.noNA$Date_time) - min(may16_filt.noNA$Date_time) # Time difference of 12.52083 days
range(may16_filt.noNA$Temp) # 17.272 28.692

## Plot the filtered data:
ggplot(may16_filt.noNA, aes(x=Date_time, y=Temp, colour = HOBO)) +
  geom_line() +
  scale_x_datetime(date_breaks = "2 days") +
  theme_classic() + theme(legend.position = "none") +
  scale_y_continuous(breaks = seq(10, 30, 2))
ggsave("Plots/FigS3c.pdf", width=12.5, height=6, units="in", dpi=1500, useDingbats=FALSE)

#### 4d. Burrow temperature summaries #### ------------------------------------------------------------------------------------------------------------------------------------------
nov17_filt.noNA_noAmb <- nov17_filt.noNA %>% filter(Length != "ambient")
nov17_filt.noNA_noAmb$Length <- as.numeric(as.character(nov17_filt.noNA_noAmb$Length))
nov17_filt.noNA_noAmb$month <- "nov17"
may16_filt.noNA$month <- "may16"
colnames(may16_filt.noNA)[8] <- "Length"
str(nov17_filt.noNA_noAmb)
str(may16_filt.noNA)
allTemps <- rbind(nov17_filt.noNA_noAmb[ ,c(1,2,8,4,9,7)], may16_filt.noNA[ ,c(6,3,8,5,9,7)])
allTemps$Burrow_ID <- paste(allTemps$HOBO, allTemps$Length, sep = "_")
length(unique(allTemps$Burrow_ID)) # 32 burrow recordings (14 from may 2016, 18 from nov 2017)

summary <- allTemps %>% group_by(Burrow_ID, Length, month, session) %>% 
  dplyr::summarise(mean = mean(Temp), median = median(Temp), min = min(Temp), max = max(Temp))
summary$range <- summary$max - summary$min

ggplot(summary, aes(x=log(Length), y=range)) +
  geom_point(size=0.5) +
  geom_smooth(method=lm, se=TRUE) +
  theme_classic() +
  xlab("log (burrow length)") + ylab("temperature range")

summary(lmer(log(Length) ~ range + (1|month), data = summary)) # Adjusted R-squared:  0.6166
summary(lm(log(Length) ~ range, data = summary))

#### 4e. Try Fourier Transform to analyze data #### -------------------------------------------------------------------------------------------------------------------------
allTemps$time <- as.numeric(allTemps$Date_time)
w <- allTemps %>% filter(Burrow_ID == "7_30") %>% dplyr::select(Burrow_ID, time, Temp) # practice with burrow 7_30
w <- w[order(w$time), ] # make sure Temp is in chronological order
traj <- w$Temp
time <- seq(1, length(traj)) # in seconds
plot(traj, type="l")

# First de-trend the data:
trend <- lm(traj ~ time)
detrend.traj <- trend$residuals
plot(detrend.traj, type="l") # plot should now look "flatter"

## Fourier transform ## From: http://www.di.fc.ul.pt/~jpn/r/fourier/fourier.html
plot.frequency.spectrum <- function(X.k, xlimits=c(0,length(X.k))) {
  plot.data  <- cbind(0:(length(X.k)-1), Mod(X.k))
  plot.data[2:length(X.k),2] <- 2*plot.data[2:length(X.k),2] 
  plot(plot.data, t="h", lwd=2, main="", 
       xlab="Frequency (Hz)", ylab="Strength", 
       xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))))
}

X.k <- fft(detrend.traj) # find all harmonics
plot.frequency.spectrum(X.k, xlimits = c(0,50))
max(Mod(X.k)) 

# Hilbert transform
library(seewave)
h <- hilbert(wave = detrend.traj, f = 1, fftw = FALSE)
env(wave = detrend.traj, f = 1, envt = "hil", plot = TRUE) # show plot
ampEnv <- env(wave = detrend.traj, f = 1, envt = "hil", plot = FALSE) # return a matrix
mean(ampEnv) # 1.609319

# Run Hilbert transform for all Burrows
allTemps$burrowNum <- as.numeric(as.factor(allTemps$Burrow_ID))

par(mfrow=c(4,4))
resultlist <- list()
for (i in 1:length(unique(allTemps$burrowNum))) {
  w <- allTemps %>% filter(burrowNum == i) %>% dplyr::select(burrowNum, time, Temp)
  w <- w[order(w$time), ] # make sure Temp is in chronological order
  traj <- w$Temp # save ordered temperature readings as a trajectory
  time <- seq(1, length(traj)) # save vector of time stamps
  trend <- lm(traj ~ time) # regress out time (i.e. correct for shifting baseline)
  detrend.traj <- trend$residuals # save residuals
  print(max(detrend.traj) - min(detrend.traj))
  plot(detrend.traj)
  ampEnv <- env(wave = detrend.traj, f = 1, envt = "hil", plot = FALSE) # get the amplitude envelope from the Hilbert transform
  #env(wave = detrend.traj, f = 1, envt = "hil", plot = TRUE)
  resultlist[[i]] = mean(ampEnv)
}
result = as.data.frame(do.call(rbind, resultlist))
result$burrowNum <- rownames(result) 
result <- join(result, allTemps[ ,c(3:7,9)], type = "left") %>% distinct()
colnames(result)[1] <- "amplitudeEnvelope"

ggplot(result, aes(x=Length, y=amplitudeEnvelope)) +
  geom_point(size=0.5) +
  geom_smooth(method=lm, se=TRUE) +
  theme_classic() +
  xlab("burrow length") + ylab("amplitude envelope")
ggsave("Plots/Fig1g.pdf", width=8, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summary(lm(amplitudeEnvelope ~ Length, data = result)) # Adjusted R-squared:  0.6381

#### 5. Draw a map of Florida #### ------------------------------------------------------------------------------------------------------------------------------------------
library(maps)
library(mapdata)
library(ggspatial)

FL <- map_data("state", region="florida")
site <- data.frame(longitude = -86.723619, latitude = 30.394862) # This is the location of the F23 RFID reader (roughly centre of grid)

ggplot(FL, aes(x=long, y=lat)) +
  geom_polygon(color="grey22", fill="antiquewhite") +
  geom_point(data=site, aes(x=longitude, y=latitude), size=8, shape=23, fill="cyan") +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.25, "in"), pad_y = unit(0.25, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_fixed(1) +
  theme(panel.grid.major=element_line(color=gray(0.5), 
                                      linetype="dashed", size=0.15), panel.background=element_rect(fill="aliceblue"))
ggsave("Plots/FigS1a.pdf", height=6, units="in", dpi=1500, useDingbats=FALSE)

## Zoom in on Santa Rosa Island:
hr <- map_data("worldHires")
usa <- hr %>% filter(region %in% "USA") %>% droplevels()
ggplot() +
  geom_polygon(data=usa, aes(x=long, y=lat, group=group), color="grey22", fill="antiquewhite") +
  geom_point(data=site, aes(x=longitude, y=latitude), size=8, shape=23, fill="cyan") +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.25, "in"), pad_y = unit(0.25, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_fixed(ylim = c(30.2,30.7), xlim = c(-87.5,-86.3), expand = FALSE) +
  theme(panel.grid.major=element_line(color=gray(0.5), 
                                      linetype="dashed", size=0.15), panel.background=element_rect(fill="aliceblue"))
ggsave("Plots/FigS1b.pdf", height=6, units="in", dpi=1500, useDingbats=FALSE)