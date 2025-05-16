library(here)
library(beeswarm)
library(tidyverse)

# read data
# TVY10 = TVY1 here so can keep the code at 4 characters for the strain name

#d <- read.csv(here("phenotypes", "invasive_growth", "data_in", "230412_IG_YST6_YST7_TVY4_TVY10.csv"))
# to be consistent with other phenotyping and sequencing, just use the 12 isolates that were sequenced
d <- read.csv(here("phenotypes", "invasive_growth", "data_in", "240527_IG_YST6_YST7_TVY4_TVY10.csv"))

d$Strain <- substr(d$name, 1, 4)
d$Place <- substr(d$name, 5, 5)
d$Replicate <- substr(d$name, 6, 7)
d$SP <- substr(d$name, 1, 5)
d$SP <- factor(d$SP, levels = c("YAG00", "YST6V", "YST6R", "TVY4V", "TVY4R", "TVY1V", "TVY1R", "YST7V", "YST7R"))

d$pos[d$Strain=="YAG0" & d$Place=="0"] <- 1
d$pos[d$Strain=="YST6" & d$Place=="V"] <- 2
d$pos[d$Strain=="YST6" & d$Place=="R"] <- 2.2
d$pos[d$Strain=="TVY4" & d$Place=="V"] <- 3
d$pos[d$Strain=="TVY4" & d$Place=="R"] <- 3.2
d$pos[d$Strain=="TVY1" & d$Place=="V"] <- 4
d$pos[d$Strain=="TVY1" & d$Place=="R"] <- 4.2
d$pos[d$Strain=="YST7" & d$Place=="V"] <- 5
d$pos[d$Strain=="YST7" & d$Place=="R"] <- 5.2

d$cols[d$Strain=="YST6" & d$Place=="V"] <- "steelblue"
d$cols[d$Strain=="YST6" & d$Place=="R"] <- "skyblue"
d$cols[d$Strain=="TVY4" & d$Place=="V"] <- "#5E478B"
d$cols[d$Strain=="TVY4" & d$Place=="R"] <- "#9270DC"
d$cols[d$Strain=="TVY1" & d$Place=="V"] <- "#EC3194"
d$cols[d$Strain=="TVY1" & d$Place=="R"] <- "#E37AB3"
d$cols[d$Strain=="YST7" & d$Place=="V"] <- "#b8860b"
d$cols[d$Strain=="YST7" & d$Place=="R"] <- "#EEDD82"
d$cols[d$Strain=="YAG0" & d$Place=="0"] <- "black"

t <- jitter(as.numeric(as.factor(d$SP)))

pdf("figures_tables/FigureX-IV_12.pdf", width=7, height=4.5)
plot(t, d$max, ylim =c(0, 1), yaxt="n", xaxt="n", ylab = "Invasive growth", xlab="", col=d$cols)
axis(2, las=2, at=c(0, 0.25, 0.5, 0.75, 1))
axis(1, at=1:9, labels=c("SC5314", "YST6V", "YST6R", "TVY4V", "TVY4R", "TVY10V", "TVY10R", "YST7V", "YST7R"), cex.axis=0.8)
points(t[97], 1, pch=19)
arrows(4, 0.85, 5, 0.85, length=0, lwd=2)
arrows(8, 0.85, 9, 0.85, length=0, lwd=2)
text(4.75, 0.88, "*", cex=1.5)
text(8.25, 0.88, "*", cex=1.5)
dev.off()

wilcox.test(subset(d, SP == "YST6R")$max, subset(d, SP == "YST6V")$max ) #W = 60, p-value = 0.3041
wilcox.test(subset(d, SP == "TVY4R")$max, subset(d, SP == "TVY4V")$max ) #W = 102, p-value = 0.04242 # R > V
wilcox.test(subset(d, SP == "TVY1R")$max, subset(d, SP == "TVY1V")$max ) #W = 76.5, p-value = 0.9414
wilcox.test(subset(d, SP == "YST7R")$max, subset(d, SP == "YST7V")$max ) #W = 28, p-value = 0.003583 # V > R

wilcox.test(subset(d, SP == "YST6R")$mean, subset(d, SP == "YST6V")$mean ) #W = 60, p-value = 0.3041
wilcox.test(subset(d, SP == "TVY4R")$mean, subset(d, SP == "TVY4V")$mean ) #W = 100, p-value = 0.08513
wilcox.test(subset(d, SP == "TVY1R")$mean, subset(d, SP == "TVY1V")$mean ) #W = 71.5, p-value = 0.677
wilcox.test(subset(d, SP == "YST7R")$mean, subset(d, SP == "YST7V")$mean ) #W = 52.5, p-value = 0.2292

###########
d$overlap <- ifelse(d$scaleB1==d$scaleB2, "purple", NA)

YST6V <- subset(d, SP == "YST6V")
YST6R <- subset(d, SP == "YST6R")
YST7V <- subset(d, SP == "YST7V")
YST7R <- subset(d, SP == "YST7R")
TVY4V <- subset(d, SP == "TVY4V")
TVY4R <- subset(d, SP == "TVY4R")

#YST6R$overlap <- ifelse(YST6R$scaleB1==YST6R$scaleB2, "purple", NA)

pdf("Figures_Tables/FigureSX-IV_12.pdf", width=6, height=7, font = "Times")
par(mfrow=c(3, 2), mar=c(1, 1, 1, 1), oma=c(3, 4.5, 1, 1))
plot(YST6V$Replicate, YST6V$scaleB1, col="blue", pch=19, ylim = c(0,1), yaxt="n", xaxt="n", ylab = "Invasive growth", xlab="", cex=0.6)
axis(2, las=2, at=c(0, 0.25, 0.5, 0.75, 1))
points(YST6V$Replicate, YST6V$scaleB2, col = "red", pch=19, cex=0.6)
points(YST6V$Replicate, YST6V$scaleB2, col = YST6V$overlap, pch=19, cex=0.6)
mtext("A) YST6V", side = 3, outer=FALSE, adj=0.01)

plot(YST6R$Replicate, YST6R$scaleB1, col="blue", pch=19, ylim = c(0,1), yaxt="n", xaxt="n", ylab = "Invasive growth", xlab="", cex=0.6)
axis(2, las=2, at=c(0, 0.25, 0.5, 0.75, 1), labels=FALSE)
points(YST6R$Replicate, YST6R$scaleB2, col = "red", pch=19, cex=0.6)
points(YST6R$Replicate, YST6R$scaleB2, col = YST6R$overlap, pch=19, cex=0.6)
mtext("B) YST6R", side = 3, outer=FALSE, adj=0.01)
legend("topright", c("Bioreplicate 1", "Bioreplicate 2", "Both bioreplicates"), col=c("blue", "red", "purple"), pch=19)

plot(TVY4V$Replicate, TVY4V$scaleB1, col="blue", pch=19, ylim = c(0,1), yaxt="n", xaxt="n", ylab = "Invasive growth", xlab="", cex=0.6)
axis(2, las=2, at=c(0, 0.25, 0.5, 0.75, 1))
points(TVY4V$Replicate, TVY4V$scaleB2, col = "red", pch=19, cex=0.6)
points(TVY4V$Replicate, TVY4V$scaleB2, col = TVY4V$overlap, pch=19, cex=0.6)
mtext("C) TVY4V", side = 3, outer=FALSE, adj=0.01)

plot(TVY4R$Replicate, TVY4R$scaleB1, col="blue", pch=19, ylim = c(0,1), yaxt="n", xaxt="n", ylab = "Invasive growth", xlab="", cex=0.6)
axis(2, las=2, at=c(0, 0.25, 0.5, 0.75, 1), labels=FALSE)
points(TVY4R$Replicate, TVY4R$scaleB2, col = "red", pch=19, cex=0.6)
points(TVY4R$Replicate, TVY4R$scaleB2, col = TVY4R$overlap, pch=19, cex=0.6)
mtext("D) TVY4R", side = 3, outer=FALSE, adj=0.01)

plot(YST7V$Replicate, YST7V$scaleB1, col="blue", pch=19, ylim = c(0,1), yaxt="n", xaxt="n", ylab = "Invasive growth", xlab="", cex=0.6)
axis(2, las=2, at=c(0, 0.25, 0.5, 0.75, 1))
points(YST7V$Replicate, YST7V$scaleB2, col = "red", pch=19, cex=0.6)
points(YST7V$Replicate, YST7V$scaleB2, col = YST7V$overlap, pch=19, cex=0.6)
mtext("E) YST7V", side = 3, outer=FALSE, adj=0.01)
axis(1, 1:24, labels=FALSE)

plot(YST7R$Replicate, YST7R$scaleB1, col="blue", pch=19, ylim = c(0,1), yaxt="n", xaxt="n", ylab = "Invasive growth", xlab="", cex=0.6)
axis(2, las=2, at=c(0, 0.25, 0.5, 0.75, 1), labels=FALSE)
points(YST7R$Replicate, YST7R$scaleB2, col = "red", pch=19, cex=0.6)
points(YST7R$Replicate, YST7R$scaleB2, col = YST7R$overlap, pch=19, cex=0.6)
mtext("F) YST7R", side = 3, outer=FALSE, adj=0.01)
axis(1, 1:24, labels=FALSE)
mtext("Invasive growth", side=2, outer=TRUE, line = 3, cex=1.25)
mtext("Replicate 1-24", side=1, outer=TRUE, line = 1, cex=1.25)
dev.off()
