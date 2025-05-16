########################
# Import libraries
########################
library(dplyr)
library(ggplot2)
library(ggpubr)
#library(PairedData)
library(tidyr)
#library(ggstatsplot)
library(janitor)
library(cowplot)
library(ggplot2)
library(ggforce)
library(here)
library(tidyverse)
#library(reshape2)
library(viridis)
#library(extrafont)
library(lmerTest)
library(EnvStats)


### NEED TO CORRECT ISOLATE NUMBERS TO BE DIFFERNET BETWEEN RECTAL AND VAGINAL!!!
########################
# Import and clean data
########################
#1. Core analyses from multiple DDAs
#pH4 <- read_csv(here("phenotypes", "DDA", "DDA_data_in", "231204_pH4_DDA_df_all.csv"))
# That file includes 24 isolates for YST6 and YST7. To be consistent with the rest of the dataset, removed the isolates that weren't sequenced.

pH4 <- read_csv(here("phenotypes", "DDA", "DDA_data_in", "240527_pH4_DDA_df_all.csv"))
pH4$RAD20[pH4$RAD80==0] <- 1
pH4$FoG20[pH4$RAD80==0] <- NA
#pH4 <- separate(pH4, col = "name", c(NA, "line", "PlaceIsol", NA, "drug"), sep="_", remove=FALSE)
pH4$type <- substr(pH4$PlaceIsol, 1, 1)
pH4$replicate <- as.numeric(substr(pH4$PlaceIsol, 2, 3))
write.csv(pH4, here("phenotypes", "DDA", "DDA_data_in", "240527pH4_DDA_df_all.csv"))

pH4$RAD20[pH4$name == "211216_YST6_V17_B_FLC_48hr"] <- 1
pH4$FoG20[pH4$name == "211216_YST6_V17_B_FLC_48hr"] <- NA
DDA_avg<- pH4 %>%
  group_by(line, type, PlaceIsol, drug) %>%
  summarise(avgRAD20 = mean(RAD20, na.rm=TRUE), avgRAD50 = mean(RAD50, na.rm=TRUE), avgFoG20 = mean(FoG20, na.rm=TRUE), avgFoG50 = mean(FoG50, na.rm=TRUE))

##par(mfrow=c(2,1), mar=c(1, 1, 1, 1))
#plot(DDA_avg$avgRAD20, DDA_avg$avgRAD50)
#plot(DDA_avg$avgFoG20, DDA_avg$avgFoG50)

DDA_avg$cols <-c()
DDA_avg$cols[DDA_avg$line=="YST6" & DDA_avg$type=="V"] <- "steelblue"
DDA_avg$cols[DDA_avg$line=="YST6" & DDA_avg$type=="R"] <- "skyblue"
DDA_avg$cols[DDA_avg$line=="YST7" & DDA_avg$type=="V"] <- "#b8860b"
DDA_avg$cols[DDA_avg$line=="YST7" & DDA_avg$type=="R"] <- "#EEDD82"
DDA_avg$cols[DDA_avg$line=="TVY4" & DDA_avg$type=="V"] <- "#5E478B"
DDA_avg$cols[DDA_avg$line=="TVY4" & DDA_avg$type=="R"] <- "#9270DC"
DDA_avg$cols[DDA_avg$line=="TVY10" & DDA_avg$type=="V"] <- "#EC3194"
DDA_avg$cols[DDA_avg$line=="TVY10" & DDA_avg$type=="R"] <- "#E37AB3"

DDA_avg$pos[DDA_avg$line=="YST6" & DDA_avg$type=="V"] <- 2
DDA_avg$pos[DDA_avg$line=="YST6" & DDA_avg$type=="R"] <- 2.2
DDA_avg$pos[DDA_avg$line=="TVY4" & DDA_avg$type=="V"] <- 3
DDA_avg$pos[DDA_avg$line=="TVY4" & DDA_avg$type=="R"] <- 3.2
DDA_avg$pos[DDA_avg$line=="TVY10" & DDA_avg$type=="V"] <- 4
DDA_avg$pos[DDA_avg$line=="TVY10" & DDA_avg$type=="R"] <- 4.2
DDA_avg$pos[DDA_avg$line=="YST7" & DDA_avg$type=="V"] <- 5
DDA_avg$pos[DDA_avg$line=="YST7" & DDA_avg$type=="R"] <- 5.2

DDA_avg$pch[DDA_avg$type=="V"] <- 2
DDA_avg$pch[DDA_avg$type=="R"] <- 0

DDA_avg_FLC <- subset(DDA_avg, drug == "FLC")
DDA_avg_BA <- subset(DDA_avg, drug == "BA")

pdf("figures_tables/Figure7-DDAall12.pdf", width=7, height=5.5, family = "ArialMT")
par(mfrow=c(2, 2), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
plot(DDA_avg_FLC$pos, DDA_avg_FLC$avgRAD20, col=DDA_avg_FLC$cols, pch = DDA_avg_FLC$pch, ylim=c(30, 0), xlab="", xaxt="n", yaxt="n")
axis(2, las=2)
axis(1, at=c(2.1, 3.1, 4.1, 5.1), labels=FALSE)
mtext("Fluconazole", side=3, adj=0.1)
mtext("Resistance (RAD20)", side=2, line=3)

plot(DDA_avg_BA$pos, DDA_avg_BA$avgRAD20, col=DDA_avg_BA$cols, pch = DDA_avg_BA$pch, ylim=c(30, 0), xlab="", xaxt="n", yaxt="n")
axis(2, las=2, labels=FALSE)
axis(1, at=c(2.1, 3.1, 4.1, 5.1), labels=FALSE)
mtext("Boric acid", side=3, adj=0.1)

plot(DDA_avg_FLC$pos, DDA_avg_FLC$avgFoG20, col=DDA_avg_FLC$cols, pch=DDA_avg_FLC$pch, ylim=c(0, 1), xlab="", xaxt="n", yaxt="n")
axis(2, las=2)
axis(1, at=c(2.1, 3.1, 4.1, 5.1), labels=c("YST6", "TVY4V", "TVY10","YST7"), cex.axis=0.8)
mtext("Tolerance (FoG20)", side=2, line=3)

plot(DDA_avg_BA$pos, DDA_avg_BA$avgFoG20, col=DDA_avg_BA$cols, pch = DDA_avg_BA$pch, ylim=c(0, 1), xlab="", xaxt="n", yaxt="n")
axis(1, at=c(2.1, 3.1, 4.1, 5.1), labels=c("YST6", "TVY4V", "TVY10","YST7"), cex.axis=0.8)
legend("topright", legend=c("vaginal", "rectal"), pch=c(2, 0), cex=0.9, bty="n")
dev.off()

#Fortalk - just resistance
par(mfrow=c(1, 2), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
plot(DDA_avg_FLC$pos, DDA_avg_FLC$avgRAD20, col=DDA_avg_FLC$cols, pch = DDA_avg_FLC$pch, ylim=c(30, 0), xlab="", xaxt="n", yaxt="n")
axis(2, las=2)
axis(1, at=c(2.1, 3.1, 4.1, 5.1), labels=c("YST6", "TVY4V", "TVY10","YST7"), cex.axis=0.8)
mtext("Fluconazole", side=3, adj=0.1)
mtext("Resistance (ZOI)", side=2, line=3)

plot(DDA_avg_BA$pos, DDA_avg_BA$avgRAD20, col=DDA_avg_BA$cols, pch = DDA_avg_BA$pch, ylim=c(30, 0), xlab="", xaxt="n", yaxt="n")
axis(2, las=2, labels=FALSE)
axis(1, at=c(2.1, 3.1, 4.1, 5.1), labels=c("YST6", "TVY4V", "TVY10","YST7"), cex.axis=0.8)
mtext("Boric acid", side=3, adj=0.1)
legend("bottomright", legend=c("vaginal", "rectal"), pch=c(2, 0), cex=0.9, bty="n")


#Separate Datasets by line
YST6_avg <- subset(DDA_avg, line== "YST6")
YST7_avg <- subset(DDA_avg, line== "YST7")
TVY4_avg <- subset(DDA_avg, line== "TVY4")
TVY10_avg <- subset(DDA_avg, line== "TVY10")

#Separate Datasets by line and drug
YST6_F <- subset(YST6_avg, drug== "FLC") #no variation
YST6_B <- subset(YST6_avg, drug== "BA")
YST7_F <- subset(YST7_avg, drug== "FLC")
YST7_B <- subset(YST7_avg, drug== "BA")
TVY4_F <- subset(TVY4_avg, drug== "FLC")
TVY4_B <- subset(TVY4_avg, drug== "BA")
TVY10_F <- subset(TVY10_avg, drug== "FLC")
TVY10_B <- subset(TVY10_avg, drug== "BA")


##########################################
# T-test for vaginal to rectal
##########################################
t.test(subset(YST6_B, type == "R")$avgRAD20, subset(YST6_B, type == "V")$avgRAD20)
#t = 0.80403, df = 21.823, p-value = 0.43

t.test(subset(TVY4_B, type == "R")$avgRAD20, subset(TVY4_B, type == "V")$avgRAD20)
#t = 0.83737, df = 20.992, p-value = 0.4118

t.test(subset(TVY10_B, type == "R")$avgRAD20, subset(TVY10_B, type == "V")$avgRAD20)
#t = -0.12487, df = 17.535, p-value = 0.902

t.test(subset(YST7_B, type == "R")$avgRAD20, subset(YST7_B, type == "V")$avgRAD20)
#t = -1.1541, df = 21.201, p-value = 0.2613

t.test(subset(YST6_B, type == "R")$avgFoG20, subset(YST6_B, type == "V")$avgFoG20)
#t = 2.7665, df = 18.235, p-value = 0.01261 (R = 0.0278, V = 0.281)

t.test(subset(TVY4_B, type == "R")$avgFoG20, subset(TVY4_B, type == "V")$avgFoG20)
#t = -0.33911, df = 17.979, p-value = 0.7385

t.test(subset(TVY10_B, type == "R")$avgFoG20, subset(TVY10_B, type == "V")$avgFoG20)
#t = 1.12, df = 19.799, p-value = 0.2761

t.test(subset(YST7_B, type == "R")$avgFoG20, subset(YST7_B, type == "V")$avgFoG20)
#t = 1.3814, df = 15.992, p-value = 0.1861

#t.test(subset(YST6_F, type == "R")$avgRAD20, subset(YST6_F, type == "V")$avgRAD20)

t.test(subset(TVY4_F, type == "R")$avgRAD20, subset(TVY4_F, type == "V")$avgRAD20)
#t = 1.6696, df = 18.828, p-value = 0.1115

t.test(subset(TVY10_F, type == "R")$avgRAD20, subset(TVY10_F, type == "V")$avgRAD20)
#t = -1.4883, df = 21.034, p-value = 0.1515

t.test(subset(YST7_F, type == "R")$avgRAD20, subset(YST7_F, type == "V")$avgRAD20)
#t = 0.5148, df = 21.568, p-value = 0.6119

#t.test(subset(YST6_F, type == "R")$avgFoG20, subset(YST6_F, type == "V")$avgFoG20)

t.test(subset(TVY4_F, type == "R")$avgFoG20, subset(TVY4_F, type == "V")$avgFoG20)
#t = -0.24957, df = 21.663, p-value = 0.8053

t.test(subset(TVY10_F, type == "R")$avgFoG20, subset(TVY10_F, type == "V")$avgFoG20)
#t = -1.3112, df = 21.425, p-value = 0.2037

t.test(subset(YST7_F, type == "R")$avgFoG20, subset(YST7_F, type == "V")$avgFoG20)
#t = -1.3808, df = 16.344, p-value = 0.1859

##########################################
# outlier test
##########################################
#start with k = 1, if one identified, increase k
Y7VF_out_RAD <- rosnerTest(YST7_F$avgRAD20, k = 1) #0
T4VF_out_RAD <- rosnerTest(TVY4_F$avgRAD20, k = 1) #0
T10VF_out_RAD <- rosnerTest(TVY10_F$avgRAD20, k = 1) #0

Y6VB_out_RAD <- rosnerTest(YST6_B$avgRAD20, k = 1) #0
Y7VB_out_RAD <- rosnerTest(YST7_B$avgRAD20, k = 1) #0
T4VB_out_RAD <- rosnerTest(TVY4_B$avgRAD20, k = 1) #0
T10VB_out_RAD <- rosnerTest(TVY10_B$avgRAD20, k = 1) #0

Y7VF_out_FoG <- rosnerTest(YST7_F$avgFoG20, k = 1) #0
T4VF_out_FoG <- rosnerTest(TVY4_F$avgFoG20, k = 1) #0
T10VF_out_FoG <- rosnerTest(TVY10_F$avgFoG20, k = 1) #0

Y6VB_out_FoG <- rosnerTest(YST6_B$avgFoG20, k = 2) #1 -> 0.5 (obs 23)
Y7VB_out_FoG <- rosnerTest(YST7_B$avgFoG20, k = 2) #1 -> 0.45 (obs 6)
T4VB_out_FoG <- rosnerTest(TVY4_B$avgFoG20, k = 1) #0
T10VB_out_FoG <- rosnerTest(TVY10_B$avgFoG20, k = 1) #0
