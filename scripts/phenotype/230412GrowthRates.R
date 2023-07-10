#Growth rate stuffs

#functions
nderiv <- function(fit, x, eps=1e-5)
  (predict(fit, x + eps) - predict(fit, x - eps))/(2 * eps)

spline.slope <- function(x, y, n=101, eps=1e-5)
  max(nderiv(loess(log(y) ~ x, degree=1, span=0.1),
             seq(min(x), max(x), length=n)), na.rm=TRUE)

spline.slope.datO <- function(d, ...)
  sapply(d[-1], spline.slope, x=d$t, ...)

spline.slope.dat <- function(d, t, ...)
    sapply(d, spline.slope, x=t, ...)

library(here)
library(beeswarm)
library(EnvStats)



#########################
#RPMI
#########################
#Load data
dRPMI_1_12 <- read.csv(here("growthRate", "data_in", "220622GC-YST6_YST701_12-RPMI.csv"))
wellsdRPMI_1 <- read.csv(here("growthRate", "data_in", "220622GC-YST6_YST701-12-RPMI_Platemap.csv"))
dRPMI_13_24 <- read.csv(here("growthRate", "data_in", "220601GC-YST6_YST713_24-RPMI.csv"))
wellsdRPMI_13 <- read.csv(here("growthRate", "data_in", "220601GC-YST6_YST713_24-RPMI_Platemap.csv"))
dRPMI_25_36 <- read.csv(here("growthRate", "data_in", "220701GC-YST6_YST725_36-RPMI.csv"))
wellsdRPMI_25 <- read.csv(here("growthRate", "data_in", "220701GC-YST6_YST725-36-RPMI_Platemap.csv"))
dRPMI_37_48 <- read.csv(here("growthRate", "data_in", "220706GC-YST6_YST737-48-RPMI.csv"))
wellsdRPMI_37 <- read.csv(here("growthRate", "data_in", "220706GC-YST6_YST737-48-RPMI_Platemap.csv"))

names(dRPMI_1_12)[1] <- "t"
names(dRPMI_13_24)[1] <- "t"
names(dRPMI_25_36)[1] <- "t"
names(dRPMI_37_48)[1] <- "t"

# YST7V - darkred
# YST7R - salmon
# YST6V - steelblue3
# YST6R - skyblue

#calculate growth rates
lgrRPMI1 <- spline.slope.dat(dRPMI_1_12[,2:97], dRPMI_1_12$t)
lgrRPMI13 <- spline.slope.dat(dRPMI_13_24[,2:97], dRPMI_13_24$t)
lgrRPMI25 <- spline.slope.dat(dRPMI_25_36[,2:97], dRPMI_25_36$t)
lgrRPMI37 <- spline.slope.dat(dRPMI_37_48[,2:97], dRPMI_37_48$t)

grRPMI1 <- data.frame(wellsdRPMI_1, lgrRPMI1)
grRPMI1$Pos[grRPMI1$Strain=="YST6" & grRPMI1$Place=="V"] <- "1"
grRPMI1$Pos[grRPMI1$Strain=="YST6" & grRPMI1$Place=="R"] <- "1.2"
grRPMI1$Pos[grRPMI1$Strain=="YST7" & grRPMI1$Place=="V"] <- "2"
grRPMI1$Pos[grRPMI1$Strain=="YST7" & grRPMI1$Place=="R"] <- "2.2"
grRPMI13 <- data.frame(wellsdRPMI_13, lgrRPMI13)
grRPMI13$Pos[grRPMI13$Strain=="YST6" & grRPMI13$Place=="V"] <- "1"
grRPMI13$Pos[grRPMI13$Strain=="YST6" & grRPMI13$Place=="R"] <- "1.2"
grRPMI13$Pos[grRPMI13$Strain=="YST7" & grRPMI13$Place=="V"] <- "2"
grRPMI13$Pos[grRPMI13$Strain=="YST7" & grRPMI13$Place=="R"] <- "2.2"
grRPMI25 <- data.frame(wellsdRPMI_25, lgrRPMI25)
grRPMI25$Pos[grRPMI25$Strain=="YST6" & grRPMI25$Place=="V"] <- "1"
grRPMI25$Pos[grRPMI25$Strain=="YST6" & grRPMI25$Place=="R"] <- "1.2"
grRPMI25$Pos[grRPMI25$Strain=="YST7" & grRPMI25$Place=="V"] <- "2"
grRPMI25$Pos[grRPMI25$Strain=="YST7" & grRPMI25$Place=="R"] <- "2.2"
grRPMI37 <- data.frame(wellsdRPMI_37, lgrRPMI37)
grRPMI37$Pos[grRPMI37$Strain=="YST6" & grRPMI37$Place=="V"] <- "1"
grRPMI37$Pos[grRPMI37$Strain=="YST6" & grRPMI37$Place=="R"] <- "1.2"
grRPMI37$Pos[grRPMI37$Strain=="YST7" & grRPMI37$Place=="V"] <- "2"
grRPMI37$Pos[grRPMI37$Strain=="YST7" & grRPMI37$Place=="R"] <- "2.2"

grRPMI1_m <- aggregate(grRPMI1[c("lgrRPMI1")], grRPMI1[c("Strain", "Place", "Replicate", "Pos")], mean)
grRPMI1_m <- grRPMI1_m[-12,]
names(grRPMI1_m)[5] <- "lgr"
grRPMI13_m <- aggregate(grRPMI13[c("lgrRPMI13")], grRPMI13[c("Strain", "Place", "Replicate", "Pos")], mean)
names(grRPMI13_m)[5] <- "lgr"
grRPMI25_m <- aggregate(grRPMI25[c("lgrRPMI25")], grRPMI25[c("Strain", "Place", "Replicate", "Pos")], mean)
names(grRPMI25_m)[5] <- "lgr"
grRPMI37_m <- aggregate(grRPMI37[c("lgrRPMI37")], grRPMI37[c("Strain", "Place", "Replicate", "Pos")], mean)
names(grRPMI37_m)[5] <- "lgr"

grRPMI_m <- rbind(grRPMI1_m, grRPMI13_m, grRPMI25_m, grRPMI37_m)

grRPMI_m$cols[grRPMI_m$Strain=="YST6" & grRPMI_m$Place=="V"] <- "steelblue"
grRPMI_m$cols[grRPMI_m$Strain=="YST6" & grRPMI_m$Place=="R"] <- "skyblue"
grRPMI_m$cols[grRPMI_m$Strain=="YST7" & grRPMI_m$Place=="V"] <- "darkred"
grRPMI_m$cols[grRPMI_m$Strain=="YST7" & grRPMI_m$Place=="R"] <- "salmon"

subset(grRPMI_m, Strain == "YST6" & Place == "V")
subset(grRPMI_m, Strain == "YST6" & Place == "R")
subset(grRPMI_m, Strain == "YST7" & Place == "V")
subset(grRPMI_m, Strain == "YST7" & Place == "R")

grRPMI_m_sub <- subset(grRPMI_m, Replicate < 25) # replicates above 25 were done with a different plate reader program so we are not using them

plot(jitter(as.numeric(grRPMI_m_sub$Pos)), grRPMI_m_sub$lgr, col=grRPMI_m_sub$col, xlim=c(0.75, 2.35), ylim=c(0, 0.5), cex=1.2, yaxt="n", xaxt="n", ylab="Growth rate (/h)", xlab="")
axis(2, las=2)
axis(1, at=1:2, labels=FALSE)

#stats
grRPMI_YST6V <- subset(grRPMI_m_sub, Strain == "YST6" & Place == "V")
grRPMI_YST6R <- subset(grRPMI_m_sub, Strain == "YST6" & Place == "R")
grRPMI_YST7V <- subset(grRPMI_m_sub, Strain == "YST7" & Place == "V")
grRPMI_YST7R <- subset(grRPMI_m_sub, Strain == "YST7" & Place == "R")

t.test(grRPMI_YST6V$lgr, grRPMI_YST6R$lgr) #t = -0.60445, df = 44.985, p-value = 0.5486
t.test(grRPMI_YST7V$lgr, grRPMI_YST7R$lgr) #t = -1.7271, df = 34.399, p-value = 0.09313

rosnerTest(grRPMI_YST6V$lgr, k = 1) # 1 (0.1692283)
rosnerTest(grRPMI_YST6R$lgr, k = 1) # 0
rosnerTest(grRPMI_YST7V$lgr, k = 1) # 0
rosnerTest(grRPMI_YST7R$lgr, k = 1) # 0

#########################
#VSM
#########################
#Load data
dVSM <- read.csv(here("growthRate", "data_in", "211008GC_YST6YST7_V1_12-VSM.csv"))
wellsVSM <- read.csv(here("growthRate", "data_in", "211008_PlateMap-VSM.csv"))

names(dVSM)[1] <- "t"

# YST7V - darkred
# YST7R - salmon
# YST6V - steelblue3
# YST6R - skyblue

#calculate growth rates
lgrVSM <- spline.slope.dat(dVSM[,2:97], dVSM$t)

grVSM <- data.frame(wellsVSM, lgrVSM)
grVSM_m <- aggregate(grVSM[c("lgrVSM")], grVSM[c("Strain", "Place", "Replicate", "Pos")], mean)

grVSM_m$col[grVSM_m$Strain=="YST6" & grVSM_m$Place=="V"] <- "steelblue"
grVSM_m$col[grVSM_m$Strain=="YST6" & grVSM_m$Place=="R"] <- "skyblue"
grVSM_m$col[grVSM_m$Strain=="YST7" & grVSM_m$Place=="V"] <- "darkred"
grVSM_m$col[grVSM_m$Strain=="YST7" & grVSM_m$Place=="R"] <- "salmon"

# stats
grVSM_YST6V <- subset(grVSM_m, Strain == "YST6" & Place == "V")
grVSM_YST6R <- subset(grVSM_m, Strain == "YST6" & Place == "R")
grVSM_YST7V <- subset(grVSM_m, Strain == "YST7" & Place == "V")
grVSM_YST7R <- subset(grVSM_m, Strain == "YST7" & Place == "R")

t.test(grVSM_YST6V$lgr, grVSM_YST6R$lgr) #t = -0.50844, df = 21.966, p-value = 0.6162
t.test(grVSM_YST7V$lgr, grVSM_YST7R$lgr) #t = 0.23908, df = 21.989, p-value = 0.8133

rosnerTest(grVSM_YST6V$lgr, k = 1) # 0
rosnerTest(grVSM_YST6R$lgr, k = 1) # 0
rosnerTest(grVSM_YST7V$lgr, k = 1) # 0
rosnerTest(grVSM_YST7R$lgr, k = 1) # 0

#########################
#One plot
#########################

pdf("figures/FigureX-GC.pdf", width=3.5, height=5.5, font = "Times")
par(mfrow=c(2, 1), mar=c(1, 1, 1, 1), oma=c(4, 4, 1, 1))
plot(jitter(as.numeric(grRPMI_m_sub$Pos)), grRPMI_m_sub$lgr, col=grRPMI_m_sub$col, xlim=c(0.75, 2.35), ylim=c(0, 0.5), cex=1.2, yaxt="n", xaxt="n", ylab="Growth rate (/h)", xlab="")
axis(2, las=2)
axis(1, at=c(1.1, 2.1), labels=FALSE)
mtext("A) RPMI", side = 3, outer=FALSE, adj=0.01)

plot(jitter(grVSM_m$Pos), grVSM_m$lgr, col=grVSM_m$col, xlim=c(0.75, 2.35), ylim=c(0, 0.5), cex=1.2, yaxt="n", xaxt="n", ylab="Growth rate (/h)", xlab="")
axis(2, las=2)
axis(1, at=c(1.1, 2.1), labels=c("YST6", "YST7"))
legend("bottomright", legend=c("YST6V", "YST6R", "YST7V", "YST7R"), col=c("steelblue", "skyblue", "darkred", "salmon"), pch=21, cex=0.9, bty="n")
mtext("B) VSM", side = 3, outer=FALSE, adj=0.01)
mtext("Growth rate (/h)", side = 2, outer=TRUE, line = 2.5, cex=1.5)
dev.off()
