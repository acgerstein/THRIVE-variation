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
YST6_YST7_pH4 <- read_csv(here("DDA", "DDA_data_in", "pH4_DDA_df_edited.csv"))
YST6_YST7_pH4$RAD20[YST6_YST7_pH4$RAD80==0] <- 1
YST6_YST7_pH4$FoG20[YST6_YST7_pH4$RAD80==0] <- NA
YST6_YST7_pH4 <- separate(YST6_YST7_pH4, col = "name", c(NA, "line", "PlaceIsol", NA, "drug"), sep="_", remove=FALSE)
YST6_YST7_pH4$type <- substr(YST6_YST7_pH4$PlaceIsol, 1, 1)
YST6_YST7_pH4$replicate <- as.numeric(substr(YST6_YST7_pH4$PlaceIsol, 2, 3))
write.csv(YST6_YST7_pH4, here("DDA", "DDA_data_in", "pH4_DDA_df_allCol.csv"))

YST6_YST7_pH4$RAD20[YST6_YST7_pH4$name == "211216_YST6_V17_B_FLC_48hr"] <- 1
YST6_YST7_pH4$FoG20[YST6_YST7_pH4$name == "211216_YST6_V17_B_FLC_48hr"] <- NA
DDA_YST6_YST7_pH4_avg<- YST6_YST7_pH4 %>%
  group_by(line, type, replicate, drug) %>%
  summarise(avgRAD20 = mean(RAD20, na.rm=TRUE), avgRAD50 = mean(RAD50, na.rm=TRUE), avgFoG20 = mean(FoG20, na.rm=TRUE), avgFoG50 = mean(FoG50, na.rm=TRUE))
DDA_YST6_YST7_pH4_avg$enviro <- "pH4"

YST6_YST7_pH7 <- read_csv(here("DDA", "DDA_data_in", "pH7_DDA_df_edited.csv"))
YST6_YST7_pH7$RAD20[YST6_YST7_pH7$RAD80==0] <- 1
YST6_YST7_pH7$FoG20[YST6_YST7_pH7$RAD80==0] <- NA
YST6_YST7_pH7 <- separate(YST6_YST7_pH7, col = "name", c( "line", "type", "replicate", NA, "drug"), sep="_", remove=FALSE)
YST6_YST7_pH7$replicate <- as.numeric(YST6_YST7_pH7$replicate)
write.csv(YST6_YST7_pH7, here("DDA", "DDA_data_in", "pH7_DDA_df_allCol.csv"))

DDA_YST6_YST7_pH7_avg<- YST6_YST7_pH7 %>%
  group_by(line, type, replicate, drug) %>%
  summarise(avgRAD20 = mean(RAD20, na.rm=TRUE), avgRAD50 = mean(RAD50, na.rm=TRUE), avgFoG20 = mean(FoG20, na.rm=TRUE), avgFoG50 = mean(FoG50, na.rm=TRUE))
DDA_YST6_YST7_pH7_avg$enviro <- "pH7"
DDA_YST6_YST7_pH7_avg$replicate <- as.numeric(DDA_YST6_YST7_pH7_avg$replicate)
DDA_avg <- rbind(DDA_YST6_YST7_pH4_avg, DDA_YST6_YST7_pH7_avg)

plot(DDA_avg$avgRAD20, DDA_avg$avgRAD50)
plot(DDA_avg$avgFoG20, DDA_avg$avgFoG50)

#New column - Line, Type, Drug
# DDA_avg$STED <- paste0(DDA_avg$line, DDA_avg$type, DDA_avg$enviro, DDA_avg$drug)
#
# #New column - Line and Type
# DDA_avg$ST <- paste0(DDA_avg$line, DDA_avg$type)
#
# #New column - Drug, Species, Type
# DDA_avg$DST <- paste0(DDA_avg$drug, DDA_avg$line, DDA_avg$type)

DDA_avg$cols[DDA_avg$line=="YST6" & DDA_avg$type=="V"] <- "steelblue"
DDA_avg$cols[DDA_avg$line=="YST6" & DDA_avg$type=="R"] <- "skyblue"
DDA_avg$cols[DDA_avg$line=="YST7" & DDA_avg$type=="V"] <- "darkred"
DDA_avg$cols[DDA_avg$line=="YST7" & DDA_avg$type=="R"] <- "salmon"

#Separate Datasets by line
YST6_avg <- subset(DDA_avg, line!= "YST7")
YST7_avg <- subset(DDA_avg, line!= "YST6")

# Subset by type, environment and drug
Y6V4F <- subset(YST6_avg, type == "V" & enviro =="pH4" & drug == "FLC")
Y6R4F <- subset(YST6_avg, type == "R" & enviro =="pH4" & drug == "FLC")
Y6V7F <- subset(YST6_avg, type == "V" & enviro =="pH7" & drug == "FLC")
Y6R7F <- subset(YST6_avg, type == "R" & enviro =="pH7" & drug == "FLC")
Y6V4B <- subset(YST6_avg, type == "V" & enviro =="pH4" & drug == "BA")
Y6R4B <- subset(YST6_avg, type == "R" & enviro =="pH4" & drug == "BA")
Y6V7B <- subset(YST6_avg, type == "V" & enviro =="pH7" & drug == "BA")
Y6R7B <- subset(YST6_avg, type == "R" & enviro =="pH7" & drug == "BA")

Y7V4F <- subset(YST7_avg, type == "V" & enviro =="pH4" & drug == "FLC")
Y7R4F <- subset(YST7_avg, type == "R" & enviro =="pH4" & drug == "FLC")
Y7V7F <- subset(YST7_avg, type == "V" & enviro =="pH7" & drug == "FLC")
Y7R7F <- subset(YST7_avg, type == "R" & enviro =="pH7" & drug == "FLC")
Y7V4B <- subset(YST7_avg, type == "V" & enviro =="pH4" & drug == "BA")
Y7R4B <- subset(YST7_avg, type == "R" & enviro =="pH4" & drug == "BA")
Y7V7B <- subset(YST7_avg, type == "V" & enviro =="pH7" & drug == "BA")
Y7R7B <- subset(YST7_avg, type == "R" & enviro =="pH7" & drug == "BA")

# Subset by environment and drug
Y6F4 <- subset(YST6_avg, enviro =="pH4" & drug == "FLC")
Y6F7 <- subset(YST6_avg, enviro =="pH7" & drug == "FLC")
Y6B4 <- subset(YST6_avg, enviro =="pH4" & drug == "BA")
Y6B7 <- subset(YST6_avg, enviro =="pH7" & drug == "BA")

Y7F4 <- subset(YST7_avg, enviro =="pH4" & drug == "FLC")
Y7F7 <- subset(YST7_avg, enviro =="pH7" & drug == "FLC")
Y7B4 <- subset(YST7_avg, enviro =="pH4" & drug == "BA")
Y7B7 <- subset(YST7_avg, enviro =="pH7" & drug == "BA")

Y6F7_out_RAD <- rosnerTest(Y6F7$avgRAD20, k = 3) # 1 (8)
Y6B4_out_RAD <- rosnerTest(Y6B4$avgRAD20, k = 3) # 2 (15.5, 7)
Y6B7_out_RAD <- rosnerTest(Y6B7$avgRAD20, k = 3)
Y7F4_out_RAD <- rosnerTest(Y7F4$avgRAD20, k = 3)
Y7F7_out_RAD <- rosnerTest(Y7F7$avgRAD20, k = 3) # 2 (7, 9.3)
Y7B4_out_RAD <- rosnerTest(Y7B4$avgRAD20, k = 3)
Y7B7_out_RAD <- rosnerTest(Y7B7$avgRAD20, k = 3)

Y6F7_out_FoG <- rosnerTest(Y6F7$avgFoG20, k = 3) # 1 (0.51)
Y6B4_out_FoG <- rosnerTest(Y6B4$avgFoG20, k = 3) # 1 (0.5)
Y6B7_out_FoG <- rosnerTest(Y6B7$avgFoG20, k = 3)
Y7F4_out_FoG <- rosnerTest(Y7F4$avgFoG20, k = 3)
Y7F7_out_FoG <- rosnerTest(Y7F7$avgFoG20, k = 3)
Y7B4_out_FoG <- rosnerTest(Y7B4$avgFoG20, k = 3) # 1 (0.45)
Y7B7_out_FoG <- rosnerTest(Y7B7$avgFoG20, k = 3)

Y6F7_out_RAD <- rosnerTest(Y6F7$avgRAD50, k = 3)
Y6B4_out_RAD <- rosnerTest(Y6B4$avgRAD50, k = 3) # 1 (8.15)
Y6B7_out_RAD <- rosnerTest(Y6B7$avgRAD50, k = 3)
Y7F4_out_RAD <- rosnerTest(Y7F4$avgRAD50, k = 3)
Y7F7_out_RAD <- rosnerTest(Y7F7$avgRAD50, k = 3)
Y7B4_out_RAD <- rosnerTest(Y7B4$avgRAD50, k = 3)
Y7B7_out_RAD <- rosnerTest(Y7B7$avgRAD50, k = 3)

Y6F7_out_FoG <- rosnerTest(Y6F7$avgFoG50, k = 3)
Y6B4_out_FoG <- rosnerTest(Y6B4$avgFoG50, k = 3) # 2 (0.57, 0.52)
Y6B7_out_FoG <- rosnerTest(Y6B7$avgFoG50, k = 3)
Y7F4_out_FoG <- rosnerTest(Y7F4$avgFoG50, k = 3)
Y7F7_out_FoG <- rosnerTest(Y7F7$avgFoG50, k = 3)
Y7B4_out_FoG <- rosnerTest(Y7B4$avgFoG50, k = 3) # 1 (0.22)
Y7B7_out_FoG <- rosnerTest(Y7B7$avgFoG50, k = 3)

#####################
## Statistics
#####################

YST6FLC <- subset(YST6_avg, drug == "FLC")
YST6BA <- subset(YST6_avg, drug == "BA")
YST7FLC <- subset(YST7_avg, drug == "FLC")
YST7BA <- subset(YST7_avg, drug == "BA")

YST6FLC$replicate[YST6FLC$type == "R"] <- YST6FLC$replicate+24
YST6BA$replicate[YST6BA$type == "R"] <- YST6BA$replicate+24
YST7FLC$replicate[YST7FLC$type == "R"] <- YST7FLC$replicate+24
YST7BA$replicate[YST7BA$type == "R"] <- YST7BA$replicate+24

#Type III Analysis of Variance Table with Satterthwaite's method
YST6FLC.lmer <- lmer(avgRAD20 ~ type * enviro + (1|replicate), data = YST6FLC)
anova(YST6FLC.lmer, type = 3)
# Sum Sq Mean Sq NumDF DenDF   F value Pr(>F)
# type           0.0     0.0     1 55.195    0.0074 0.9319
# enviro      8008.7  8008.7     1 55.199 1207.3457 <2e-16 ***
#   type:enviro    0.0     0.0     1 55.199    0.0074 0.9319

YST6BA.lmer <- lmer(avgRAD20 ~ enviro*type + (1|replicate), data = YST6BA)
anova(YST6BA.lmer, type = 3)
# Sum Sq Mean Sq NumDF DenDF  F value Pr(>F)
# enviro      634.65  634.65     1 45.999 305.7975 <2e-16 ***
# type          1.68    1.68     1 45.999   0.8098 0.3729
# enviro:type  11.34   11.34     1 45.999   5.4658 0.0238 *

YST7FLC.lmer <- lmer(avgRAD20 ~ enviro*type + (1|replicate), data = YST7FLC)
anova(YST7FLC.lmer, type = 3)
Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)
# enviro      98.213  98.213     1    46 70.0261 8.526e-11 ***
#   type         4.651   4.651     1    46  3.3159   0.07512 .
# enviro:type  0.113   0.113     1    46  0.0809   0.77738

YST7BA.lmer <- lmer(avgRAD20 ~ enviro * type + (1|replicate), data = YST7BA)
anova(YST7BA.lmer, type = 3)
#              Sum Sq Mean Sq NumDF  DenDF  F value Pr(>F)
# enviro      299.862 299.862     1    46 295.9676 <2e-16 ***
#type          0.439   0.439     1    46   0.4336 0.5135
#enviro:type   2.894   2.894     1    46   2.8559 0.0978 .

#YST6FLC.lmer_FoG <- lmer(avgFoG20 ~ type * enviro + (1|replicate), data = YST6FLC)
#anova(YST6FLC.lmer_FoG, type = 3)
YST6FLC.T_FoG <- t.test(subset(Y6F7, type == "R")$avgFoG20, subset(Y6F7, type == "V")$avgFoG20, paired=TRUE)
YST6FLC.T_FoG
#t = 0.33913, df = 23, p-value = 0.7376


YST6BA.lmer_FoG <- lmer(avgFoG20 ~ enviro * type + (1|replicate), data = YST6BA)
anova(YST6BA.lmer_FoG, type = 3)
# Sum Sq   Mean Sq NumDF DenDF F value   Pr(>F)
# enviro      0.0139362 0.0139362     1    92  7.9242 0.005966 **
#   type        0.0002779 0.0002779     1    92  0.1580 0.691914
# enviro:type 0.0013625 0.0013625     1    92  0.7747 0.381047

YST7FLC.lmer_FoG <- lmer(avgFoG20 ~ enviro * type + (1|replicate), data = YST7FLC)
anova(YST7FLC.lmer_FoG, type = 3)
# Sum Sq   Mean Sq NumDF DenDF F value   Pr(>F)
# enviro      0.0247416 0.0247416     1    92  8.5207 0.004414 **
#   type        0.0063891 0.0063891     1    92  2.2003 0.141401
# enviro:type 0.0009281 0.0009281     1    92  0.3196 0.573199

YST7BA.lmer_FoG <- lmer(avgFoG20 ~ enviro * type + (1|replicate), data = YST7BA)
anova(YST7BA.lmer_FoG, type = 3)
# Sum Sq Mean Sq NumDF DenDF  F value Pr(>F)
# enviro      0.35183 0.35183     1    46 323.1165 <2e-16 ***
#   type        0.00177 0.00177     1    46   1.6259 0.2087
# enviro:type 0.00232 0.00232     1    46   2.1283 0.1514


###################
#Figures
###################
# Main figure

pdf("figures/FigureX-DDA_20.pdf", width=8, height=6, font = "Times")
par(mfrow=c(2, 2), oma=c(4, 4, 1, 1), mar=c(1, 1, 1, 1))
plot(rep(1, nrow(Y6V4F)), Y6V4F$avgRAD20, xlim=c(0.7, 5.1), ylim=c(25, 0), yaxt="n", xaxt="n", col=Y6V4F$cols)
arrows(0.8, median(Y6V4F$avgRAD20), 1.2, median(Y6V4F$avgRAD20), lwd=2, col=Y6V4F$cols[1], length=0)
points(rep(1.4, nrow(Y6R4F)), Y6R4F$avgRAD20, col=Y6R4F$cols)
arrows(1.2, median(Y6R4F$avgRAD20), 1.6, median(Y6R4F$avgRAD20), lwd=2, col=Y6R4F$cols[1], length=0)
points(rep(2, nrow(Y6V7F)), Y6V7F$avgRAD20, col=Y6V7F$cols)
arrows(1.8, median(Y6V7F$avgRAD20), 2.2, median(Y6V7F$avgRAD20), lwd=2, col=Y6V7F$cols[1], length=0)
points(rep(2.4, nrow(Y6R7F)), Y6R7F$avgRAD20, col=Y6R7F$cols)
arrows(2.2, median(Y6R7F$avgRAD20), 2.6, median(Y6R7F$avgRAD20), lwd=2, col=Y6R7F$cols[1], length=0)
points(rep(3.4, nrow(Y6V4B)), Y6V4B$avgRAD20, col=Y6V4B$cols)
arrows(3.2, median(Y6V4B$avgRAD20), 3.6, median(Y6V4B$avgRAD20), lwd=2, col=Y6V4B$cols[1], length=0)
points(rep(3.8, nrow(Y6R4B)), Y6R4B$avgRAD20, col=Y6R4B$cols)
arrows(3.6, median(Y6R4B$avgRAD20), 4, median(Y6R4B$avgRAD20), lwd=2, col=Y6R4B$cols[1], length=0)
points(rep(4.4, nrow(Y6V7B)), Y6V7B$avgRAD20, col=Y6V7B$cols)
arrows(4.2, median(Y6V7B$avgRAD20), 4.6, median(Y6V7B$avgRAD20), lwd=2, col=Y6V7B$cols[1], length=0)
points(rep(4.8, nrow(Y6R7B)), Y6R7B$avgRAD20, col=Y6R7B$cols)
arrows(4.6, median(Y6R7B$avgRAD20), 5, median(Y6R7B$avgRAD20), lwd=2, col=Y6R7B$cols[1], length=0)
abline(v=2.9, lty=2)
axis(2, las=2)
axis(1, at=c(1, 1.4, 2, 2.4, 3.4, 3.8, 4.4, 4.8), labels=FALSE)
mtext("FLC", side=3, adj=0.22)
mtext("BA", side=3, adj=0.75)
mtext("Resistance (RAD20)", side=2, line=3)

plot(rep(1, nrow(Y7V4F)), Y7V4F$avgRAD20, xlim=c(0.7, 5.1), ylim=c(25, 0), yaxt="n", xaxt="n", col=Y7V4F$cols)
arrows(0.8, median(Y7V4F$avgRAD20), 1.2, median(Y7V4F$avgRAD20), lwd=2, col=Y7V4F$cols[1], length=0)
points(rep(1.4, nrow(Y7R4F)), Y7R4F$avgRAD20, col=Y7R4F$cols)
arrows(1.2, median(Y7R4F$avgRAD20), 1.6, median(Y7R4F$avgRAD20), lwd=2, col=Y7R4F$cols[1], length=0)
points(rep(2, nrow(Y7V7F)), Y7V7F$avgRAD20, col=Y7V7F$cols)
arrows(1.8, median(Y7V7F$avgRAD20), 2.2, median(Y7V7F$avgRAD20), lwd=2, col=Y7V7F$cols[1], length=0)
points(rep(2.4, nrow(Y7R7F)), Y7R7F$avgRAD20, col=Y7R7F$cols)
arrows(2.2, median(Y7R7F$avgRAD20), 2.6, median(Y7R7F$avgRAD20), lwd=2, col=Y7R7F$cols[1], length=0)
points(rep(3.4, nrow(Y7V4B)), Y7V4B$avgRAD20, col=Y7V4B$cols)
arrows(3.2, median(Y7V4B$avgRAD20), 3.6, median(Y7V4B$avgRAD20), lwd=2, col=Y7V4B$cols[1], length=0)
points(rep(3.8, nrow(Y7R4B)), Y7R4B$avgRAD20, col=Y7R4B$cols)
arrows(3.6, median(Y7R4B$avgRAD20), 4, median(Y7R4B$avgRAD20), lwd=2, col=Y7R4B$cols[1], length=0)
points(rep(4.4, nrow(Y7V7B)), Y7V7B$avgRAD20, col=Y7V7B$cols)
arrows(4.2, median(Y7V7B$avgRAD20), 4.6, median(Y7V7B$avgRAD20), lwd=2, col=Y7V7B$cols[1], length=0)
points(rep(4.8, nrow(Y7R7B)), Y7R7B$avgRAD20, col=Y7R7B$cols)
arrows(4.6, median(Y7R7B$avgRAD20), 5, median(Y7R7B$avgRAD20), lwd=2, col=Y7R7B$cols[1], length=0)
abline(v=2.9, lty=2)
axis(2, las=2, labels=FALSE)
axis(1, at=c(1, 1.4, 2, 2.4, 3.4, 3.8, 4.4, 4.8), labels=FALSE)
mtext("FLC", side=3, adj=0.22)
mtext("BA", side=3, adj=0.75)

plot(rep(1, nrow(Y6V4F)), Y6V4F$avgFoG20, xlim=c(0.7, 5.1), ylim=c(0, 1), yaxt="n", xaxt="n", col=Y6V4F$cols)
arrows(0.8, median(Y6V4F$avgFoG20), 1.2, median(Y6V4F$avgFoG20), lwd=2, col=Y6V4F$cols[1], length=0)
points(rep(1.4, nrow(Y6R4F)), Y6R4F$avgFoG20, col=Y6R4F$cols)
arrows(1.2, median(Y6R4F$avgFoG20), 1.6, median(Y6R4F$avgFoG20), lwd=2, col=Y6R4F$cols[1], length=0)
points(rep(2, nrow(Y6V7F)), Y6V7F$avgFoG20, col=Y6V7F$cols)
arrows(1.8, median(Y6V7F$avgFoG20), 2.2, median(Y6V7F$avgFoG20), lwd=2, col=Y6V7F$cols[1], length=0)
points(rep(2.4, nrow(Y6R7F)), Y6R7F$avgFoG20, col=Y6R7F$cols)
arrows(2.2, median(Y6R7F$avgFoG20), 2.6, median(Y6R7F$avgFoG20), lwd=2, col=Y6R7F$cols[1], length=0)
points(rep(3.4, nrow(Y6V4B)), Y6V4B$avgFoG20, col=Y6V4B$cols)
arrows(3.2, median(Y6V4B$avgFoG20), 3.6, median(Y6V4B$avgFoG20), lwd=2, col=Y6V4B$cols[1], length=0)
points(rep(3.8, nrow(Y6R4B)), Y6R4B$avgFoG20, col=Y6R4B$cols)
arrows(3.6, median(Y6R4B$avgFoG20), 4, median(Y6R4B$avgFoG20), lwd=2, col=Y6R4B$cols[1], length=0)
points(rep(4.4, nrow(Y6V7B)), Y6V7B$avgFoG20, col=Y6V7B$cols)
arrows(4.2, median(Y6V7B$avgFoG20), 4.6, median(Y6V7B$avgFoG20), lwd=2, col=Y6V7B$cols[1], length=0)
points(rep(4.8, nrow(Y6R7B)), Y6R7B$avgFoG20, col=Y6R7B$cols)
arrows(4.6, median(Y6R7B$avgFoG20), 5, median(Y6R7B$avgFoG20), lwd=2, col=Y6R7B$cols[1], length=0)
abline(v=2.9, lty=2)
axis(2, las=2)
Xlabels <- c("V\npH4", "R\npH4", "V\npH7", "R\npH7", "V\npH4", "R\npH4", "V\npH7", "R\npH7")
axis(1, at=c(1, 1.4, 2, 2.4, 3.4, 3.8, 4.4, 4.8), labels=Xlabels, gap.axis = 0.5, tick=FALSE)
mtext("Tolerance (FoG20)", side=2, line=3)
mtext("YST6", side=1, line=2.2, cex=1.2)

plot(rep(1, nrow(Y7V4F)), Y7V4F$avgFoG20, xlim=c(0.7, 5.1), ylim=c(0, 1), yaxt="n", xaxt="n", col=Y7V4F$cols)
arrows(0.8, median(Y7V4F$avgFoG20), 1.2, median(Y7V4F$avgFoG20), lwd=2, col=Y7V4F$cols[1], length=0)
points(rep(1.4, nrow(Y7R4F)), Y7R4F$avgFoG20, col=Y7R4F$cols)
arrows(1.2, median(Y7R4F$avgFoG20), 1.6, median(Y7R4F$avgFoG20), lwd=2, col=Y7R4F$cols[1], length=0)
points(rep(2, nrow(Y7V7F)), Y7V7F$avgFoG20, col=Y7V7F$cols)
arrows(1.8, median(Y7V7F$avgFoG20), 2.2, median(Y7V7F$avgFoG20), lwd=2, col=Y7V7F$cols[1], length=0)
points(rep(2.4, nrow(Y7R7F)), Y7R7F$avgFoG20, col=Y7R7F$cols)
arrows(2.2, median(Y7R7F$avgFoG20), 2.6, median(Y7R7F$avgFoG20), lwd=2, col=Y7R7F$cols[1], length=0)
points(rep(3.4, nrow(Y7V4B)), Y7V4B$avgFoG20, col=Y7V4B$cols)
arrows(3.2, median(Y7V4B$avgFoG20), 3.6, median(Y7V4B$avgFoG20), lwd=2, col=Y7V4B$cols[1], length=0)
points(rep(3.8, nrow(Y7R4B)), Y7R4B$avgFoG20, col=Y7R4B$cols)
arrows(3.6, median(Y7R4B$avgFoG20), 4, median(Y7R4B$avgFoG20), lwd=2, col=Y7R4B$cols[1], length=0)
points(rep(4.4, nrow(Y7V7B)), Y7V7B$avgFoG20, col=Y7V7B$cols)
arrows(4.2, median(Y7V7B$avgFoG20), 4.6, median(Y7V7B$avgFoG20), lwd=2, col=Y7V7B$cols[1], length=0)
points(rep(4.8, nrow(Y7R7B)), Y7R7B$avgFoG20, col=Y7R7B$cols)
arrows(4.6, median(Y7R7B$avgFoG20), 5, median(Y7R7B$avgFoG20), lwd=2, col=Y7R7B$cols[1], length=0)
abline(v=2.9, lty=2)
axis(2, las=2, labels=FALSE)
axis(1, at=c(1, 1.4, 2, 2.4, 3.4, 3.8, 4.4, 4.8), labels=Xlabels, gap.axis = 0.5, tick=FALSE)
legend("topright", legend=c("YST6V", "YST6R", "YST7V", "YST7R"), col=c("steelblue", "skyblue", "darkred", "salmon"), pch=21)
mtext("YST7", side=1, line=2.2, cex=1.2)
dev.off()

##############
# RAD50/FoG50
###############
#Type III Analysis of Variance Table with Satterthwaite's method
YST6FLC.lmer <- lmer(avgRAD50 ~ type * enviro + (1|replicate), data = YST6FLC)
anova(YST6FLC.lmer, type = 3)
# Sum Sq Mean Sq NumDF  DenDF  F value Pr(>F)
# type           0.42    0.42     1 55.170   0.2844  0.596
# enviro      1214.82 1214.82     1 55.174 826.8487 <2e-16 ***
#   type:enviro    0.42    0.42     1 55.174   0.2844  0.596

YST6BA.lmer <- lmer(avgRAD50 ~ enviro*type + (1|replicate), data = YST6BA)
anova(YST6BA.lmer, type = 3)
#Type III Analysis of Variance Table with Satterthwaite's method
#              Sum Sq Mean Sq NumDF DenDF  F value Pr(>F)
# enviro      279.598 279.598     1    92 320.6904 <2e-16 ***
# type          0.907   0.907     1    92   1.0408 0.3103
# enviro:type   1.260   1.260     1    92   1.4457 0.2323

YST7FLC.lmer <- lmer(avgRAD50 ~ enviro*type + (1|replicate), data = YST7FLC)
anova(YST7FLC.lmer, type = 3)
# Sum Sq Mean Sq NumDF DenDF  F value Pr(>F)
# enviro      357.60  357.60     1    46 832.0827 <2e-16 ***
#   type          0.09    0.09     1    46   0.2175 0.6432
# enviro:type   0.19    0.19     1    46   0.4447 0.5082

YST7BA.lmer <- lmer(avgRAD50 ~ enviro * type + (1|replicate), data = YST7BA)
anova(YST7BA.lmer, type = 3)
#             Sum Sq Mean Sq NumDF DenDF   F value  Pr(>F)
# enviro      523.06  523.06     1    46 1108.9249 < 2e-16 ***
# type          0.06    0.06     1    46    0.1167 0.73422
# enviro:type   2.75    2.75     1    46    5.8316 0.01977 *

#YST6FLC.lmer_FoG <- lmer(avgFoG20 ~ type * enviro + (1|replicate), data = YST6FLC)
#anova(YST6FLC.lmer_FoG, type = 3)

YST6FLC.T_FoG <- t.test(subset(Y6F7, type == "R")$avgFoG50, subset(Y6F7, type == "V")$avgFoG50, paired=TRUE)
YST6FLC.T_FoG
#t = 0.28762, df = 23, p-value = 0.7762

YST6BA.lmer_FoG <- lmer(avgFoG50 ~ enviro * type + (1|replicate), data = YST6BA)
anova(YST6BA.lmer_FoG, type = 3)
# Sum Sq  Mean Sq NumDF DenDF  F value Pr(>F)
# enviro      0.312817 0.312817     1    92 122.8127 <2e-16 ***
#   type        0.001709 0.001709     1    92   0.6708 0.4149
# enviro:type 0.002570 0.002570     1    92   1.0088 0.3178

YST7FLC.lmer_FoG <- lmer(avgFoG50 ~ enviro * type + (1|replicate), data = YST7FLC)
anova(YST7FLC.lmer_FoG, type = 3)
#              Sum Sq Mean Sq NumDF DenDF  F value Pr(>F)
# enviro      1.07696 1.07696     1    92 138.4393 <2e-16 ***
# type        0.00442 0.00442     1    92   0.5681 0.4530
# enviro:type 0.00130 0.00130     1    92   0.1675 0.6833

YST7BA.lmer_FoG <- lmer(avgFoG50 ~ enviro * type + (1|replicate), data = YST7BA)
anova(YST7BA.lmer_FoG, type = 3)
#               Sum Sq  Mean Sq NumDF DenDF F value    Pr(>F)
# enviro      0.061890 0.061890     1    92 245.400 < 2.2e-16 ***
# type        0.001564 0.001564     1    92   6.202 0.0145553 *
# enviro:type 0.003374 0.003374     1    92  13.380 0.0004235 ***

pdf("figures/FigureX-DDA_20.pdf", width=8, height=6, font = "Times")
par(mfrow=c(2, 2), oma=c(4, 4, 1, 1), mar=c(1, 1, 1, 1))
plot(rep(1, nrow(Y6V4F)), Y6V4F$avgRAD50, xlim=c(0.7, 5.1), ylim=c(25, 0), yaxt="n", xaxt="n", col=Y6V4F$cols)
arrows(0.8, median(Y6V4F$avgRAD50), 1.2, median(Y6V4F$avgRAD50), lwd=2, col=Y6V4F$cols[1], length=0)
points(rep(1.4, nrow(Y6R4F)), Y6R4F$avgRAD50, col=Y6R4F$cols)
arrows(1.2, median(Y6R4F$avgRAD50), 1.6, median(Y6R4F$avgRAD50), lwd=2, col=Y6R4F$cols[1], length=0)
points(rep(2, nrow(Y6V7F)), Y6V7F$avgRAD50, col=Y6V7F$cols)
arrows(1.8, median(Y6V7F$avgRAD50), 2.2, median(Y6V7F$avgRAD50), lwd=2, col=Y6V7F$cols[1], length=0)
points(rep(2.4, nrow(Y6R7F)), Y6R7F$avgRAD50, col=Y6R7F$cols)
arrows(2.2, median(Y6R7F$avgRAD50), 2.6, median(Y6R7F$avgRAD50), lwd=2, col=Y6R7F$cols[1], length=0)
points(rep(3.4, nrow(Y6V4B)), Y6V4B$avgRAD50, col=Y6V4B$cols)
arrows(3.2, median(Y6V4B$avgRAD50), 3.6, median(Y6V4B$avgRAD50), lwd=2, col=Y6V4B$cols[1], length=0)
points(rep(3.8, nrow(Y6R4B)), Y6R4B$avgRAD50, col=Y6R4B$cols)
arrows(3.6, median(Y6R4B$avgRAD50), 4, median(Y6R4B$avgRAD50), lwd=2, col=Y6R4B$cols[1], length=0)
points(rep(4.4, nrow(Y6V7B)), Y6V7B$avgRAD50, col=Y6V7B$cols)
arrows(4.2, median(Y6V7B$avgRAD50), 4.6, median(Y6V7B$avgRAD50), lwd=2, col=Y6V7B$cols[1], length=0)
points(rep(4.8, nrow(Y6R7B)), Y6R7B$avgRAD50, col=Y6R7B$cols)
arrows(4.6, median(Y6R7B$avgRAD50), 5, median(Y6R7B$avgRAD50), lwd=2, col=Y6R7B$cols[1], length=0)
abline(v=2.9, lty=2)
axis(2, las=2)
axis(1, at=c(1, 1.4, 2, 2.4, 3.4, 3.8, 4.4, 4.8), labels=FALSE)
mtext("FLC", side=3, adj=0.22)
mtext("BA", side=3, adj=0.75)
mtext("Resistance (RAD20)", side=2, line=3)

plot(rep(1, nrow(Y7V4F)), Y7V4F$avgRAD50, xlim=c(0.7, 5.1), ylim=c(25, 0), yaxt="n", xaxt="n", col=Y7V4F$cols)
arrows(0.8, median(Y7V4F$avgRAD50), 1.2, median(Y7V4F$avgRAD50), lwd=2, col=Y7V4F$cols[1], length=0)
points(rep(1.4, nrow(Y7R4F)), Y7R4F$avgRAD50, col=Y7R4F$cols)
arrows(1.2, median(Y7R4F$avgRAD50), 1.6, median(Y7R4F$avgRAD50), lwd=2, col=Y7R4F$cols[1], length=0)
points(rep(2, nrow(Y7V7F)), Y7V7F$avgRAD50, col=Y7V7F$cols)
arrows(1.8, median(Y7V7F$avgRAD50), 2.2, median(Y7V7F$avgRAD50), lwd=2, col=Y7V7F$cols[1], length=0)
points(rep(2.4, nrow(Y7R7F)), Y7R7F$avgRAD50, col=Y7R7F$cols)
arrows(2.2, median(Y7R7F$avgRAD50), 2.6, median(Y7R7F$avgRAD50), lwd=2, col=Y7R7F$cols[1], length=0)
points(rep(3.4, nrow(Y7V4B)), Y7V4B$avgRAD50, col=Y7V4B$cols)
arrows(3.2, median(Y7V4B$avgRAD50), 3.6, median(Y7V4B$avgRAD50), lwd=2, col=Y7V4B$cols[1], length=0)
points(rep(3.8, nrow(Y7R4B)), Y7R4B$avgRAD50, col=Y7R4B$cols)
arrows(3.6, median(Y7R4B$avgRAD50), 4, median(Y7R4B$avgRAD50), lwd=2, col=Y7R4B$cols[1], length=0)
points(rep(4.4, nrow(Y7V7B)), Y7V7B$avgRAD50, col=Y7V7B$cols)
arrows(4.2, median(Y7V7B$avgRAD50), 4.6, median(Y7V7B$avgRAD50), lwd=2, col=Y7V7B$cols[1], length=0)
points(rep(4.8, nrow(Y7R7B)), Y7R7B$avgRAD50, col=Y7R7B$cols)
arrows(4.6, median(Y7R7B$avgRAD50), 5, median(Y7R7B$avgRAD50), lwd=2, col=Y7R7B$cols[1], length=0)
abline(v=2.9, lty=2)
axis(2, las=2, labels=FALSE)
axis(1, at=c(1, 1.4, 2, 2.4, 3.4, 3.8, 4.4, 4.8), labels=FALSE)
mtext("FLC", side=3, adj=0.22)
mtext("BA", side=3, adj=0.75)

plot(rep(1, nrow(Y6V4F)), Y6V4F$avgFoG50, xlim=c(0.7, 5.1), ylim=c(0, 1), yaxt="n", xaxt="n", col=Y6V4F$cols)
arrows(0.8, median(Y6V4F$avgFoG50), 1.2, median(Y6V4F$avgFoG50), lwd=2, col=Y6V4F$cols[1], length=0)
points(rep(1.4, nrow(Y6R4F)), Y6R4F$avgFoG50, col=Y6R4F$cols)
arrows(1.2, median(Y6R4F$avgFoG50), 1.6, median(Y6R4F$avgFoG50), lwd=2, col=Y6R4F$cols[1], length=0)
points(rep(2, nrow(Y6V7F)), Y6V7F$avgFoG50, col=Y6V7F$cols)
arrows(1.8, median(Y6V7F$avgFoG50), 2.2, median(Y6V7F$avgFoG50), lwd=2, col=Y6V7F$cols[1], length=0)
points(rep(2.4, nrow(Y6R7F)), Y6R7F$avgFoG50, col=Y6R7F$cols)
arrows(2.2, median(Y6R7F$avgFoG50), 2.6, median(Y6R7F$avgFoG50), lwd=2, col=Y6R7F$cols[1], length=0)
points(rep(3.4, nrow(Y6V4B)), Y6V4B$avgFoG50, col=Y6V4B$cols)
arrows(3.2, median(Y6V4B$avgFoG50), 3.6, median(Y6V4B$avgFoG50), lwd=2, col=Y6V4B$cols[1], length=0)
points(rep(3.8, nrow(Y6R4B)), Y6R4B$avgFoG50, col=Y6R4B$cols)
arrows(3.6, median(Y6R4B$avgFoG50), 4, median(Y6R4B$avgFoG50), lwd=2, col=Y6R4B$cols[1], length=0)
points(rep(4.4, nrow(Y6V7B)), Y6V7B$avgFoG50, col=Y6V7B$cols)
arrows(4.2, median(Y6V7B$avgFoG50), 4.6, median(Y6V7B$avgFoG50), lwd=2, col=Y6V7B$cols[1], length=0)
points(rep(4.8, nrow(Y6R7B)), Y6R7B$avgFoG50, col=Y6R7B$cols)
arrows(4.6, median(Y6R7B$avgFoG50), 5, median(Y6R7B$avgFoG50), lwd=2, col=Y6R7B$cols[1], length=0)
abline(v=2.9, lty=2)
axis(2, las=2)
Xlabels <- c("V\npH4", "R\npH4", "V\npH7", "R\npH7", "V\npH4", "R\npH4", "V\npH7", "R\npH7")
axis(1, at=c(1, 1.4, 2, 2.4, 3.4, 3.8, 4.4, 4.8), labels=Xlabels, gap.axis = 0.5, tick=FALSE)
mtext("Tolerance (FoG20)", side=2, line=3)
mtext("YST6", side=1, line=2.2, cex=1.2)

plot(rep(1, nrow(Y7V4F)), Y7V4F$avgFoG50, xlim=c(0.7, 5.1), ylim=c(0, 1), yaxt="n", xaxt="n", col=Y7V4F$cols)
arrows(0.8, median(Y7V4F$avgFoG50), 1.2, median(Y7V4F$avgFoG50), lwd=2, col=Y7V4F$cols[1], length=0)
points(rep(1.4, nrow(Y7R4F)), Y7R4F$avgFoG50, col=Y7R4F$cols)
arrows(1.2, median(Y7R4F$avgFoG50), 1.6, median(Y7R4F$avgFoG50), lwd=2, col=Y7R4F$cols[1], length=0)
points(rep(2, nrow(Y7V7F)), Y7V7F$avgFoG50, col=Y7V7F$cols)
arrows(1.8, median(Y7V7F$avgFoG50), 2.2, median(Y7V7F$avgFoG50), lwd=2, col=Y7V7F$cols[1], length=0)
points(rep(2.4, nrow(Y7R7F)), Y7R7F$avgFoG50, col=Y7R7F$cols)
arrows(2.2, median(Y7R7F$avgFoG50), 2.6, median(Y7R7F$avgFoG50), lwd=2, col=Y7R7F$cols[1], length=0)
points(rep(3.4, nrow(Y7V4B)), Y7V4B$avgFoG50, col=Y7V4B$cols)
arrows(3.2, median(Y7V4B$avgFoG50), 3.6, median(Y7V4B$avgFoG50), lwd=2, col=Y7V4B$cols[1], length=0)
points(rep(3.8, nrow(Y7R4B)), Y7R4B$avgFoG50, col=Y7R4B$cols)
arrows(3.6, median(Y7R4B$avgFoG50), 4, median(Y7R4B$avgFoG50), lwd=2, col=Y7R4B$cols[1], length=0)
points(rep(4.4, nrow(Y7V7B)), Y7V7B$avgFoG50, col=Y7V7B$cols)
arrows(4.2, median(Y7V7B$avgFoG50), 4.6, median(Y7V7B$avgFoG50), lwd=2, col=Y7V7B$cols[1], length=0)
points(rep(4.8, nrow(Y7R7B)), Y7R7B$avgFoG50, col=Y7R7B$cols)
arrows(4.6, median(Y7R7B$avgFoG50), 5, median(Y7R7B$avgFoG50), lwd=2, col=Y7R7B$cols[1], length=0)
abline(v=2.9, lty=2)
axis(2, las=2, labels=FALSE)
axis(1, at=c(1, 1.4, 2, 2.4, 3.4, 3.8, 4.4, 4.8), labels=Xlabels, gap.axis = 0.5, tick=FALSE)
legend("topright", legend=c("YST6V", "YST6R", "YST7V", "YST7R"), col=c("steelblue", "skyblue", "darkred", "salmon"), pch=21)
mtext("YST7", side=1, line=2.2, cex=1.2)
dev.off()

# Combined figure

# Control figure
DDA_controls$de <- paste(DDA_controls$drug, DDA_controls$enviro, sep="_")
yAG003 <- subset(DDA_controls, line == "yAG003")
yAG054 <- subset(DDA_controls, line == "yAG054")
yAG116 <- subset(DDA_controls, line == "yAG116")
yAG126 <- subset(DDA_controls, line == "yAG126")

pdf("figures/manuscript/FigureSX-DDA_controls.pdf", width=7.5, height=4, font = "Times")
par(mfrow=c(2, 2), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
plot(as.numeric(as.factor(yAG003$de)), yAG003$RAD20, xaxt="n", ylim=c(30, 0), xlim=c(0.5, 4.5))
points(c(1, 2), c(DDA_controls_avg$avgRAD20[1], DDA_controls_avg$avgRAD20[3]), type="l")
points(c(3, 4), c(DDA_controls_avg$avgRAD20[2], DDA_controls_avg$avgRAD20[4]), type="l")
axis(1, at = 1:4, labels=FALSE)
mtext("SC5314", side=3, adj=0.01)
mtext("Resistance (RAD20)", side=2, line=2)
plot(as.numeric(as.factor(yAG054$de)), yAG054$RAD20, xaxt="n", ylim=c(30, 0), xlim=c(0.5, 4.5), yaxt="n")
points(c(1, 2), c(DDA_controls_avg$avgRAD20[5], DDA_controls_avg$avgRAD20[7]), type="l")
points(c(3, 4), c(DDA_controls_avg$avgRAD20[6], DDA_controls_avg$avgRAD20[8]), type="l")
axis(2, labels=FALSE)
axis(1, at = 1:4, labels=FALSE)
mtext("T101", side=3, adj=0.01)

plot(as.numeric(as.factor(yAG003$de)), yAG003$FoG20, xaxt="n", ylim=c(0,1), xlim=c(0.5, 4.5))
points(c(1, 2), c(DDA_controls_avg$avgFoG20[1], DDA_controls_avg$avgFoG20[3]), type="l")
points(c(3, 4), c(DDA_controls_avg$avgFoG20[2], DDA_controls_avg$avgFoG20[4]), type="l")
axis(1, at = 1:4, labels=c("BA pH4", "BA pH7", "FLC pH4", "FLC pH7"))
mtext("SC5314", side=3, adj=0.01)
mtext("Tolerance (FoG20)", side=2, line=2)
plot(as.numeric(as.factor(yAG054$de)), yAG054$FoG20, xaxt="n", ylim=c(0,1), xlim=c(0.5, 4.5), yaxt="n")
points(c(1, 2), c(DDA_controls_avg$avgFoG20[5], DDA_controls_avg$avgFoG20[7]), type="l")
points(c(3, 4), c(DDA_controls_avg$avgFoG20[6], DDA_controls_avg$avgFoG20[8]), type="l")
axis(2, labels=FALSE)
axis(1, at = 1:4, labels=c("BA pH4", "BA pH7", "FLC pH4", "FLC pH7"))
dev.off()

#unused
par(mfrow=c(2, 1))
plot(c(1, 3, 2, 4), DDA_controls_avg$avgRAD20[1:4], ylim=c(25, 0), col="red", pch=c(19, 19, 21, 21), xlim=c(0.5, 4.5), xaxt="n")
points(c(1, 3, 2, 4), DDA_controls_avg$avgRAD20[5:8], col="blue", pch=c(19, 19, 21, 21))
points(c(1, 3, 2, 4), DDA_controls_avg$avgRAD20[9:12], col="orange", pch=c(19, 19, 21, 21))
points(c(1, 3, 2, 4), DDA_controls_avg$avgRAD20[13:16], col="purple", pch=c(19, 19, 21, 21))
abline(v=2.5)
axis(1, at=c(1:4), labels=FALSE)

plot(c(1, 3, 2, 4), DDA_controls_avg$avgFoG20[1:4], ylim=c(0, 1), col="red", pch=c(19, 19, 21, 21), xlim=c(0.5, 4.5), xaxt="n")
points(c(1, 3, 2, 4), DDA_controls_avg$avgFoG20[5:8], col="blue", pch=c(19, 19, 21, 21))
points(c(1, 3, 2, 4), DDA_controls_avg$avgFoG20[9:12], col="orange", pch=c(19, 19, 21, 21))
points(c(1, 3, 2, 4), DDA_controls_avg$avgFoG20[13:16], col="purple", pch=c(19, 19, 21, 21))
abline(v=2.5)
axis(1, at=c(1:4), labels=c("pH4", "pH7", "pH4", "pH7"))

### HERE
#Plot YST6 RAD80
YST6_RAD80 <- ggplot(data=YST6_avg,mapping=aes(x = STED,y=avgRAD80, col=type))+
  geom_sina(alpha=0.8,size=3)+
  xlab("Drug")+
  ylab("Resistance (mm)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=16), axis.title = element_text(size=80))+
 ylim(23,0)+
 scale_color_manual(values = c("Salmon","Red"),name="",labels=c("C. krusei Rectal","C krusei Vaginal"))+
  stat_summary(fun = median,
               fun.min = median,
               fun.max = median,
               geom = "crossbar",
              col= "black",
               width = 0.5)
#  annotate("text",x=c(1,4,5,6),y=c(9.5,8.5,11,9.5),label = "*")
#annotate("segment",x=c(0.75,3.75,4.75,5.75),xend = c(1.25,4.25,5.25,6.25),y = c(10,9,11.5,10), yend =c(10,9,11.5,10))

YST6_FoG80 <- ggplot(data=DDA_YST6_avg,mapping=aes(x= Drug_Type,y=avgFoG80,col=type))+
  geom_sina(alpha=0.7,size=3)+
  xlab("Drug")+
  ylab("Tolerance")+
 theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=16), axis.title = element_text(size=80))+
  scale_color_manual(values = c("Salmon","Red"),name="",labels=c("C. krusei Rectal","C. krusei Vaginal"))+
  ylim(0,1)+
  stat_summary(fun = median,
               fun.min = median,
              fun.max = median,
               geom = "crossbar",
              col= "black",
              width = 0.5)
# annotate("text",x=c(1,3,4,5,6,7),y=c(0.67,0.62,0.7,0.64,0.71,0.51),label = "*")
#annotate("segment",x=c(0.75,2.75,3.75,4.75,5.75,6.75),xend = c(1.25,3.25,4.25,5.25,6.25,7.25),y = c(0.66,0.61,0.69,0.63,0.7,0.5), yend =c(0.66,0.61,0.69,0.63,0.7,0.5))

#Saving YST6 resistance and tolerance
DDA_YST6_80 <- plot_grid(YST6_RAD80,YST6_FoG80,nrow = 2,align = "hv")
DDA_YST6_80
ggsave("280524DDA_YST6.pdf",DDA_YST6_80, width = 7, height=10)

######################################################################################################
#Pivot data - YST7 - Drug and type
DDA_YST7_avg$Drug_Type <- paste(DDA_YST7_avg$drug, DDA_YST7_avg$type)

#Plot YST7 RAD80
YST7_RAD80 <- ggplot(data=DDA_YST7_avg,mapping=aes(x = Drug_Type,y=avgRAD80,col=type))+
  geom_sina(alpha=0.8,size=3)+
  xlab("Drug")+
  ylab("Resistance (mm)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=16), axis.title = element_text(size=80))+
  ylim(23,0)+
  scale_color_manual(values = c("steelblue3", "darkred"),name="",labels="")+
  stat_summary(fun = median,
              fun.min = median,
               fun.max = median,
               geom = "crossbar",
               col= "black",
               width = 0.5)
#  annotate("text",x=c(1,4,5,6),y=c(9.5,8.5,11,9.5),label = "*")
#annotate("segment",x=c(0.75,3.75,4.75,5.75),xend = c(1.25,4.25,5.25,6.25),y = c(10,9,11.5,10), yend =c(10,9,11.5,10))

YST7_FoG80 <- ggplot(data=DDA_YST7_avg,mapping=aes(x= Drug_Type,y=avgFoG80,col=type))+
  geom_sina(alpha=0.7,size=3)+
  xlab("Drug")+
  ylab("Tolerance")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=16), axis.title = element_text(size=80))+
  scale_color_manual(values = c("steelblue3", "darkred"),name="",labels="")+
  ylim(0,1)+
  stat_summary(fun = median,
               fun.min = median,
               fun.max = median,
               geom = "crossbar",
               col= "black",
             width = 0.5)
# annotate("text",x=c(1,3,4,5,6,7),y=c(0.67,0.62,0.7,0.64,0.71,0.51),label = "*")
#annotate("segment",x=c(0.75,2.75,3.75,4.75,5.75,6.75),xend = c(1.25,3.25,4.25,5.25,6.25,7.25),y = c(0.66,0.61,0.69,0.63,0.7,0.5), yend =c(0.66,0.61,0.69,0.63,0.7,0.5))

#Saving YST7 resistance and tolerance
DDA_YST7_80 <- plot_grid(YST7_RAD80,YST7_FoG80,nrow = 2,align = "hv")
DDA_YST7_80
ggsave("DDA_YST7.pdf",DDA_YST7_80, width = 7, height=10)

##########################################################################3
#Graphing YST6 and YST7 resistance on one graph
#RAD80
YST6_YST7_RAD80 <- ggplot(data=DDA_YST6_YST7_avg, mapping=aes(x = Drug_Species_Type, y=avgRAD80,col= Species_Type))+
  geom_sina(alpha=0.8,size=3)+
  xlab("Drug")+
  ylab("Resistance (mm)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=16), axis.title = element_text(size=80))+
  ylim(18,0)+
  scale_color_manual(values = c("Salmon","Red", "skyblue2", "mediumblue"),name="",labels=c("C. krusei Rectal","C. krusei Vaginal", "C. albicans Rectal", "C. albicans Vaginal"))+
  stat_summary(fun = median,
               fun.min = median,
               fun.max = median,
               geom = "crossbar",
               col= "black",
               width = 0.5)
YST6_YST7_RAD80

#Fog80
YST6_YST7_FoG80 <- ggplot(data=DDA_YST6_YST7_avg, mapping=aes(x = Species_Type_Drug, y=avgFoG80,col= Species_Type))+
  geom_sina(alpha=0.8,size=3)+
  xlab("Drug")+
  ylab("Tolerance")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=16), axis.title = element_text(size=80))+
  ylim(0,1)+
  scale_color_manual(values = c("Salmon","Red", "skyblue2", "mediumblue"),name="",labels=c("C. krusei Rectal","C. krusei Vaginal", "C. albicans Rectal", "C. albicans Vaginal"))+
  stat_summary(fun = median,
               fun.min = median,
               fun.max = median,
               geom = "crossbar",
               col= "black",
               width = 0.5)
YST6_YST7_FoG80



