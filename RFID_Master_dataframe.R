## Nicole Bedford
## Hoekstra Lab
## February - November 2018
## Make a useable Master dataframe of RFID pings per burrow per night ##
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(purrr)
library(sp)
library(rgdal)

# Read in the data:
setwd("~/Dropbox (Personal)/Fieldwork_Paper_NLB/Good/Code/ForGitHub")

pings <- read.table("RFID_rawPings_2017.csv", header=T, sep=",")  # 2014x32=64448
decode <- read.table("RFID_Decode2017.txt", header=T, colClasses = "factor")
burrow <- read.table("uniqueBur_2017.csv", header=T, sep=",")
specs <- read.table("BurrowSpecs_Nov17.csv", header=T, sep=",")

# Convert Burrow information to long format:
burrow[burrow==""] <- NA                                          # convert blank cells to NA
b <- gather(burrow, "Date", "Reader", 2:23, na.rm=TRUE)           # convert to long format
b$Date <- gsub("[.]", ":", b$Date)                                # replace . with : 
b$Date <- gsub("X", "", b$Date)                                   # remove X
b$Date <- mdy(b$Date)                                             # Date format now good!
b <- left_join(b, specs, by="Burrow")                             # add extra burrow features

# Add Night Info:
b$Night <- NA
b$Night <- ifelse(b$Date == "2017-11-20", 1,
                  ifelse(b$Date == "2017-11-21", 2,
                         ifelse(b$Date == "2017-11-22", 3,
                                ifelse(b$Date == "2017-11-23", 4,
                                       ifelse(b$Date == "2017-11-24", 5,
                                              ifelse(b$Date == "2017-11-25", 6,
                                                     ifelse(b$Date == "2017-11-26", 7,
                                                            ifelse(b$Date == "2017-11-27", 8,
                                                                   ifelse(b$Date == "2017-11-28", 9,
                                                                          ifelse(b$Date == "2017-11-29", 10,
                                                                                 ifelse(b$Date == "2017-11-30", 11, NA)))))))))))

# Convert RFID Ping information to long format:
pings[pings==""] <- NA                                            # convert blank cells to NA (43599)
sum(is.na(pings))                                                 # 43599 NAs
p <- gather(pings, "Reader", "info", 1:32, na.rm=TRUE)            # 20849 non-blank cells (64448-43599)
length(unique(p$Reader))                                          # there should be info for 32 Readers
p <- p[!grepl("New Recording Session", p$info), ]                 # 19662 cells contain real data
p <- p %>% separate(info, c("RFID","Time","Date"), " ")           # split info into separate columns
p$Timestamp <- paste(p$Date, p$Time)
p$Timestamp <- mdy_hms(p$Timestamp)
p$num <- as.numeric(p$Timestamp)                                  # converts timestamp into seconds since "origin"
hist(p$Timestamp, breaks=50)                                      # Timestamps for R19 and R30
p$Date <- mdy(p$Date)
p <- p[ ,c(1,2,4,5)]

## Fix timestamps for R19 & R30
# Reader 19
R19 <- p %>% filter(Reader %in% "R19")
start19 <- as.data.frame(as.POSIXct(c("2017-11-14 16:30:00", "2017-11-15 16:10:00", "2017-11-16 16:25:00", "2017-11-19 17:00:00", "2017-11-20 17:00:00", "2017-11-21 16:22:00", "2017-11-22 16:44:00", 
                                      "2017-11-23 16:40:00", "2017-11-24 16:49:00", "2017-11-25 16:41:00", "2017-11-26 16:36:00", "2017-11-27 16:29:00", "2017-11-28 16:21:00", "2017-11-29 17:04:00")))
colnames(start19) <- "Timestamp"
start19$num <- as.numeric(start19$Timestamp)
R19a <- read.table("R19_adjust.csv", header=T, sep=",")
R19a$goodDate <- as.POSIXct(R19a$Good_num, origin = "1940-01-02 00:00:00") # Convert numeric back to Date/Time format
R19a <- cbind(R19, R19a)
R19a <- R19a[ ,c(1,2,16)]
colnames(R19a)[3] <- "Timestamp"
R19a <- R19a[complete.cases(R19a), ]                              # Put this dataframe back into full set, 336 pings
R19a$Date <- as.character(R19a$Timestamp)
R19a$Date <- map(strsplit(R19a$Date, " "), 1)
R19a$Date <- as.Date(unlist(R19a$Date))
R19a <- R19a[ ,c(1,2,4,3)] # 336

# Reader 30
R30 <- p %>% filter(Reader %in% "R30")
start30 <- as.data.frame(as.POSIXct(c("2017-11-17 17:00:00", "2017-11-18 17:04:00", "2017-11-19 16:39:00", "2017-11-20 16:39:00", "2017-11-21 16:20:00", "2017-11-22 17:04:00", "2017-11-23 16:57:00", 
                                      "2017-11-24 16:51:00", "2017-11-25 16:56:00", "2017-11-26 16:57:00", "2017-11-27 16:51:00", "2017-11-28 16:45:00", "2017-11-29 16:45:00", "2017-11-30 16:07:00")))
colnames(start30) <- "Timestamp"
start30$num <- as.numeric(start30$Timestamp)
R30a <- read.table("R30_adjust.csv", header=T, sep=",")
R30a$goodDate <- as.POSIXct(R30a$Good_num, origin = "1940-01-02 00:00:00") # Convert numeric back to Date/Time format
R30a <- cbind(R30, R30a)
R30a <- R30a[ ,c(1,2,16)]
colnames(R30a)[3] <- "Timestamp"
R30a <- R30a[complete.cases(R30a), ]                              # Put this dataframe back into full set, 606 pings
R30a$Date <- as.character(R30a$Timestamp)
R30a$Date <- map(strsplit(R30a$Date, " "), 1)
R30a$Date <- as.Date(unlist(R30a$Date))
R30a <- R30a[ ,c(1,2,4,3)] # 606

## Add new timestamp adjusted R19 & R30 into ping dataframe:
cut <- c("R19", "R30")
p <- p[!p$Reader %in% cut, ]                                      # 18649
p <- rbind(p, R19a, R30a)                                         # 19591
p$Mouse_ID <- NA
p$Mouse_ID <- decode$Mouse_ID[match(unlist(p$RFID), decode$RFID)]
p <- p[!duplicated(p), ]                                          # 15402
p <- p[complete.cases(p), ]                                       # 15360

# Add Night Info:
p$Night <- NA
p$Night <- ifelse(p$Timestamp >= as.POSIXct('2017-11-20 12:00') & p$Timestamp < as.POSIXct('2017-11-21 12:00'), 1,
                   ifelse(p$Timestamp >= as.POSIXct('2017-11-21 12:00') & p$Timestamp < as.POSIXct('2017-11-22 12:00'), 2,
                          ifelse(p$Timestamp >= as.POSIXct('2017-11-22 12:00') & p$Timestamp < as.POSIXct('2017-11-23 12:00'), 3,
                                 ifelse(p$Timestamp >= as.POSIXct('2017-11-23 12:00') & p$Timestamp < as.POSIXct('2017-11-24 12:00'), 4,
                                        ifelse(p$Timestamp >= as.POSIXct('2017-11-24 12:00') & p$Timestamp < as.POSIXct('2017-11-25 12:00'), 5,
                                               ifelse(p$Timestamp >= as.POSIXct('2017-11-25 12:00') & p$Timestamp < as.POSIXct('2017-11-26 12:00'), 6,
                                                      ifelse(p$Timestamp >= as.POSIXct('2017-11-26 12:00') & p$Timestamp < as.POSIXct('2017-11-27 12:00'), 7,
                                                             ifelse(p$Timestamp >= as.POSIXct('2017-11-27 12:00') & p$Timestamp < as.POSIXct('2017-11-28 12:00'), 8,
                                                                    ifelse(p$Timestamp >= as.POSIXct('2017-11-28 12:00') & p$Timestamp < as.POSIXct('2017-11-29 12:00'), 9,
                                                                           ifelse(p$Timestamp >= as.POSIXct('2017-11-29 12:00') & p$Timestamp < as.POSIXct('2017-11-30 12:00'), 10,
                                                                                  ifelse(p$Timestamp >= as.POSIXct('2017-11-30 12:00') & p$Timestamp < as.POSIXct('2017-12-01 12:00'), 11, NA)))))))))))

# Make good dataframe
df <- merge(p, b, by=c("Reader", "Night"), all=TRUE)              # 49570
df <- df[,-8]                                                     # Remove extra Date column
colnames(df)[4] <- "Date"
length(unique(df$Mouse_ID))                                       # 44 mice
table(df$Mouse_ID)
starters <- c("st1", "st2", "st3", "st4")
df <- df[!df$Mouse_ID %in% starters, ]                            # 40486
df <- df[complete.cases(df[,1:6]),]                               # 8314
data <- inner_join(df, decode, by=c("Mouse_ID", "RFID"))
table(data$Burrow) # some burrows have no pings, check
hist(data$Timestamp, breaks=200)

## Convert Lat/Long to UTM
# Transform the Coordinate Reference System (CRS) to UTM by setting the EPSG to 32616 for WGS 84, UTM Zone 16N
# http://spatialreference.org/ref/epsg/wgs-84-utm-zone-16n/
cord.dec = SpatialPoints(cbind(df$Lon, -df$Lat), proj4string=CRS("+proj=longlat"))  
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:32616"))
par(mfrow = c(1, 2))
plot(cord.dec, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(cord.UTM, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
utm <- as.data.frame(cord.UTM)
data <- cbind(data, utm)
nov <- data # 8314 obs. of 19 var.
nov <- nov[ ,c(1,3:16,2,17:20)] # reorder columns

##### Put together the May 2016 dataframe #### -------------------------------------
# Read in the data
pings <- read.table("RFID_rawPings_2016.csv", header=T, sep=",")  # 1156x25=28900
decode <- read.table("RFID_Decode2016.txt", header=T, colClasses = "factor")
burrow <- read.table("uniqueBur_2016.csv", header=T, sep=",")
specs <- read.table("BurrowSpecs_May16.csv", header=T, sep=",")

# Convert Burrow information to long format:
burrow[burrow==""] <- NA                                          # convert blank cells to NA
b <- gather(burrow, "Date", "Reader", 2:13, na.rm=TRUE)           # convert to long format
b$Date <- gsub("[.]", ":", b$Date)                                # replace . with : 
b$Date <- gsub("X", "", b$Date)                                   # remove X
b$Date <- mdy(b$Date)                                             # Date format now good!
b <- left_join(b, specs, by="Burrow")                             # add extra burrow features

# Add Night Info:
b$Night <- NA
b$Night <- ifelse(b$Date == "2016-05-11", 1,
                  ifelse(b$Date == "2016-05-12", 2,
                         ifelse(b$Date == "2016-05-13", 3,
                                ifelse(b$Date == "2016-05-14", 4,
                                       ifelse(b$Date == "2016-05-15", 5,
                                              ifelse(b$Date == "2016-05-16", 6,
                                                     ifelse(b$Date == "2016-05-17", 7,
                                                            ifelse(b$Date == "2016-05-18", 8, NA))))))))

# Convert RFID Ping information to long format:
pings[pings==""] <- NA                                            # convert blank cells to NA
sum(is.na(pings))                                                 # 23755 NAs
p <- gather(pings, "Reader", "info", 1:25, na.rm=TRUE)            # 5145 non-blank cells (28900-23755)
length(unique(p$Reader))                                          # there should be info for 25 Readers
p <- p[!grepl("New Recording Session", p$info), ]                 # 4862 cells contain real data
p <- p %>% separate(info, c("RFID","Time","Date"), " ")           # split info into separate columns
p$Timestamp <- paste(p$Date, p$Time)
p$Timestamp <- mdy_hms(p$Timestamp)
p$num <- as.numeric(p$Timestamp)                                  # converts timestamp into seconds since "origin"
hist(p$Timestamp, breaks=50)
p$Date <- mdy(p$Date)
p <- p[ ,c(1,2,4,5)]
p$Mouse_ID <- NA
p$Mouse_ID <- decode$Mouse_ID[match(unlist(p$RFID), decode$RFID)] # 4862
p <- p[!duplicated(p), ]                                          # 3741
p <- p[complete.cases(p), ]                                       # 3741

# Add Night Info:
p$Night <- NA
p$Night <- ifelse(p$Timestamp >= as.POSIXct('2016-05-11 12:00') & p$Timestamp < as.POSIXct('2016-05-12 12:00'), 1,
                   ifelse(p$Timestamp >= as.POSIXct('2016-05-12 12:00') & p$Timestamp < as.POSIXct('2016-05-13 12:00'), 2,
                          ifelse(p$Timestamp >= as.POSIXct('2016-05-13 12:00') & p$Timestamp < as.POSIXct('2016-05-14 12:00'), 3,
                                 ifelse(p$Timestamp >= as.POSIXct('2016-05-14 12:00') & p$Timestamp < as.POSIXct('2016-05-15 12:00'), 4,
                                        ifelse(p$Timestamp >= as.POSIXct('2016-05-15 12:00') & p$Timestamp < as.POSIXct('2016-05-16 12:00'), 5,
                                               ifelse(p$Timestamp >= as.POSIXct('2016-05-16 12:00') & p$Timestamp < as.POSIXct('2016-05-17 12:00'), 6,
                                                      ifelse(p$Timestamp >= as.POSIXct('2016-05-17 12:00') & p$Timestamp < as.POSIXct('2016-05-18 12:00'), 7, NA)))))))

# Make good dataframe
df <- merge(p, b, by=c("Reader", "Night"), all=TRUE)              # 3138
df <- df[,-8]                                                     # Remove extra Date column
colnames(df)[4] <- "Date"
length(unique(df$Mouse_ID))                                       # 32 mice
table(df$Mouse_ID)
df <- df[complete.cases(df[,1:6]),]                               # 3030
data <- inner_join(df, decode, by=c("Mouse_ID", "RFID"))
table(data$Burrow) # some burrows have no pings, check
hist(data$Timestamp, breaks=200)

## Convert Lat/Long to UTM
# Transform the Coordinate Reference System (CRS) to UTM by setting the EPSG to 32616 for WGS 84, UTM Zone 16N
# http://spatialreference.org/ref/epsg/wgs-84-utm-zone-16n/
cord.dec = SpatialPoints(cbind(df$Lon, -df$Lat), proj4string=CRS("+proj=longlat"))  
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:32616"))
par(mfrow = c(1, 2))
plot(cord.dec, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(cord.UTM, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)
utm <- as.data.frame(cord.UTM)
data <- cbind(data, utm)
may <- data # 3030 obs. of 20 var.
may <- may[ ,c(1,3:16,2,17:20)] # reorder columns