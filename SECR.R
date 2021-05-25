## Nicole Bedford, Hoekstra Lab
## Mark-Recapture analyis for Florida beach mouse fieldwork (Nov 2017)
## Analyse mark-recapture data in R using the secr package
## secr is "spatially explicit capture-recapture in R" based on models and papers by Murray Efford (2016)
## One system of units is used throughout secr, distances are in metres and areas are in hectares (ha)
## The unit of density for 2-dimensional habitat is animals per hectare
## 1 ha = 10000 m2 = 0.01 km2. To convert density to animals per km2, multiply by 100
## g0 is the probability of capture at the homerange centre
## sigma is the scaling function (in metres)

## Note: 7 mice caught ONLY in Accessory traps: 116, 121, 128, 129, 130, 131, 132
## 116, 121, 129 appear in RFID dataset; 128, 130, 131, 132 do not

library(dplyr)
library(plyr)
library(reshape)
library(ggplot2)
library(secr)

setwd("~/Dropbox (Personal)/Fieldwork_Paper_NLB/Good/Code/ForGitHub")

## Make "good" capt files. Save these and read-in below ## -----------------------------------------------------------------------------------------------------------------------
## November 2017 capt file (raw data):
capt <- read.table("Capt_Nov17_raw.csv", sep=",", header=T) # Read-in 2017 mark-capture information
table(capt$Type, capt$Recap)
colnames(capt)[2] <- "Detector"
capt$Session <- "BeachMouse"
capt$Occasion <- as.integer(as.factor(capt$Date))
length(unique(capt$Detector)) # 75 unique traps (Grid + Accessory)
length(unique(capt$ID)) # 43 unique mice (Grid + Accessory)
length(unique(capt$Date)) # 18 trap nights (Grid + Accessory)

## Make plot of mark-recapture information by night (Grid Only):
Gcapt <- capt %>% filter(Type == "G")
length(unique(Gcapt$Detector)) # 69 unique traps (Grid only)
length(unique(Gcapt$ID)) # 36 unique mice (Grid only)
length(unique(Gcapt$Date)) # 16 trap nights (Grid only)

markRecap <- Gcapt %>% group_by(Occasion, Recap) %>% summarise(n=n())

ggplot(markRecap, aes(x = Occasion, y = n, fill = Recap)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("grey", "orange")) +
  xlab("Occasion") +
  ylab("Number of Mice") +
  theme_classic()
ggsave("Plots/FigS1d.pdf", width=12, height=6, units="cm", dpi=1500, useDingbats=FALSE)

## Clean up capt dataframe for SECR analysis:
capt[61,16] <- -1 # Animals sometimes die on capture: mark these detections with a minus sign before the occasion number
capt[90,16] <- -4 # Animals sometimes die on capture: mark these detections with a minus sign before the occasion number
capt <- capt %>% filter(Type=="G") %>% select(Session, ID, Occasion, Detector, Sex) # remove accessory traps
length(unique(capt$Detector)) # 69 unique traps (Grid only)
length(unique(capt$ID)) # 36 unique mice (Grid only)

## Read-in grid trap location information: (this is necessary to provide rownames for usage matrices)
gridTrap17 <- read.table("GridTrap17.txt", row.names = 1)

#### Make a "capthist" object for November 2017 data:
mouse17 <- read.capthist(captfile = "capt17.txt", trapfile = "GridTrap17.txt", detector = "multi") # full dataset (all nights)

pdf("Plots/FigS1c.pdf", width=18, height=5)
plot(mouse17, type = "n.per.detector", ncap = TRUE) # does not include accessory traps
dev.off()

## 1. Try secr.fit with full dataset: (all nights, uniform usage matrix) ## ----------------------------------------------------------------------------------------------------
suggest.buffer(mouse17) # 80, Warning message: using automatic 'detectpar' g0 = 0.01187, sigma = 23.93 
fit17_full <- secr.fit(mouse17, model=list(D~1, g0~1, sigma~1), buffer=80, biasLimit=NA, detectfn='HN', trace=TRUE)
predict(fit17_full)
# link    estimate SE.estimate          lcl         ucl
# D       log  2.77356015 0.491052480  1.965602278  3.91362790
# g0    logit  0.01065237 0.001513292  0.008060414  0.01406598
# sigma   log 28.98020987 2.077410312 25.186183888 33.34576480

## Check buffer width:
esa.plot(fit17_full)
suggest.buffer(fit17_full) # 99; if I re-run the model with buffer=99, I get the same answers!

## 2a. Include variable usage matrix for reporting in paper: ## -------------------------------------------------------------------------------------------------------------------
## See "Grid Notes" in Google Docs
mat17 <- matrix(1, nrow=630, ncol=16) # numTraps x numNights
rownames(mat17) <- rownames(gridTrap17)
mat17[c(1:2,13:72,83:142,153:212,223:282,293:352,363:422,433:492,503:630), 1] <- 0
mat17[c(1:2,44:72,114:142,184:212,254:282,324:352,394:422,464:492,534:630), 2] <- 0
mat17[c(47:70,117:140,187:210,257:280,327:350,397:420,467:490,537:560,607:630), c(3:5)] <- 0
mat17[c(56:70,126:140,196:210,266:280,336:350,406:420,476:490,546:560,616:630), c(6:7)] <- 0
mat17[c(1:15,71:85,141:155,211:225,281:295,351:365,421:435,491:505,561:575), c(8:14)] <- 0
mat17[c(1:43,71:113,141:183,211:253,281:323,351:393,421:463,491:533,561:603), c(15)] <- 0
mat17[c(1:55,71:125,141:195,211:265,281:335,351:405,421:475,491:545,561:615), c(16)] <- 0

## Apply usage matrix to full dataset (all 16 nights)
usage(traps(mouse17)) <- mat17 # add usage information from binary matrix above
summary(mouse17) # Supplementary Table 1 for Paper

## Secr model doesn't fit well for nights with < 400 traps --> remove nights 1,2,15,16
mouse17_mat <- read.capthist(captfile = "capt17_sub.txt", trapfile = "GridTrap17.txt", detector = "multi") # make new capthist object with only nights 3-14
usage(traps(mouse17_mat)) <- mat17[ ,3:14] # add usage information from binary matrix above (only nights 3-14)

## 2b. Try secr.fit with partial dataset: (12 nights, correct usage matrix) 
suggest.buffer(mouse17_mat) # 80, Warning message: using automatic 'detectpar' g0 = 0.01187, sigma = 23.93 
fit17_mat <- secr.fit(mouse17_mat, model=list(D~1, g0~1, sigma~1), buffer=80, biasLimit=NA, detectfn='HN', trace=TRUE)
predict(fit17_mat)
# link    estimate SE.estimate         lcl        ucl
# D       log  2.53994987 0.479621617  1.75993890  3.6656644
# g0    logit  0.01888389 0.003036172  0.01376763  0.0258516
# sigma   log 29.96809816 2.463566207 25.51542748 35.1977998

## Make SECR half-normal function plot:
pdf("Plots/FigS1e.pdf", width = 6, height = 6)
plot(fit17_mat, sigmatick=TRUE, limits=TRUE, alpha=0.5, xval=0:120)
dev.off()

## Calculate home-range size from sigma:
r <- circular.r(p=0.95, detectfn="HN", sigma=29.96809816) # 67
pi*(r^2) # 14212.7 m2 home-range size estimate

## Calculate population density from D (mice/ha):
2.53994987 * 100 # 253.995 mice/km2

## How many mice do we expect in the trap grid? (include 80m buffer)
h = 80 + (2*80) # grid height
w = 690 + (2*80) # grid width
(h*w)/1e+6 # 0.204 km2 trap grid area
253.995 * 0.204 # 51.81498 mice (we caught 36 mice in grid traps, 43 in grid + accessory)
36/52 # 69% (G)
43/52 # 83% (G+A)
32/52 # 62% (RFID)