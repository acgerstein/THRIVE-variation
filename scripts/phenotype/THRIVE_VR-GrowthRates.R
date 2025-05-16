library(here)
library(beeswarm)
library(EnvStats)
library(lmerTest)

#Growth rate stuffs
#functions to calculate growth rate
nderiv <- function(fit, x, eps=1e-5)
  (predict(fit, x + eps) - predict(fit, x - eps))/(2 * eps)

spline.slope <- function(x, y, n=101, eps=1e-5)
  max(nderiv(loess(log(y) ~ x, degree=1, span=0.1),
             seq(min(x), max(x), length=n)), na.rm=TRUE)

spline.slope.datO <- function(d, ...)
  sapply(d[-1], spline.slope, x=d$t, ...)

spline.slope.dat <- function(d, t, ...)
    sapply(d, spline.slope, x=t, ...)



#########################
#Load OD kinetic data
#########################

RPMI_1 <- read.csv(here("phenotypes","growthRate", "data_in", "231103_THRIVE_YST6_YST7_TVY4_TVY10_48h_RPMI.csv"), header=FALSE)
RPMI_1 <- RPMI_1[1:193,]
RPMI_2 <- read.csv(here("phenotypes","growthRate", "data_in", "231117_THRIVE_YST6_YST7_TVY4_TVY10_48h_RPMI.csv"), header=FALSE)
RPMI_2 <- RPMI_2[1:193,]
VSM_1 <- read.csv(here("phenotypes","growthRate", "data_in", "231117_THRIVE_YST6_YST7_TVY4_TVY10_48h_VSM.csv"), header=FALSE)
VSM_1 <- VSM_1[1:193,]
VSM_2 <- read.csv(here("phenotypes","growthRate", "data_in", "231122_THRIVE_YST6_YST7_TVY4_TVY10_48h_VSM.csv"), header=FALSE)
VSM_2 <- VSM_2[1:193,]

# change the name of the first column to t
names(RPMI_1)[1] <- "t"
names(RPMI_2)[1] <- "t"
names(VSM_1)[1] <- "t"
names(VSM_2)[1] <- "t"

# read in wells information
wells <- read.csv(here("phenotypes","growthRate", "data_in", "231201_PlateMap-RPMIVSM.csv"), header=TRUE)

# desired colour info
# YST7V - darkred
# YST7R - salmon
# YST6V - steelblue3
# YST6R - skyblue

#calculate growth rates (round calculation to 3 decimal places)
lgrRPMI1 <- round(spline.slope.dat(RPMI_1[,2:97], RPMI_1$t), 2)
lgrRPMI2 <- round(spline.slope.dat(RPMI_2[,2:97], RPMI_2$t), 2)
lgrVSM1 <- round(spline.slope.dat(VSM_1[,2:97], VSM_1$t), 2)
lgrVSM2 <- round(spline.slope.dat(VSM_2[,2:97], VSM_2$t), 2)

# attach growth rate calculation to the wells
grRPMI1 <- data.frame(wells, lgrRPMI1)
grRPMI1$BioRep <- "A"
names(grRPMI1)[6] <- "lgr"
grRPMI2 <- data.frame(wells, lgrRPMI2)
grRPMI2$BioRep <- "B"
names(grRPMI2)[6] <- "lgr"
grVSM1 <- data.frame(wells, lgrVSM1)
grVSM1$BioRep <- "A"
names(grVSM1)[6] <- "lgr"
grVSM2 <- data.frame(wells, lgrVSM2)
grVSM2$BioRep <- "B"
names(grVSM2)[6] <- "lgr"

# join two RPMI dataset together
grRPMI12 <- rbind(grRPMI1, grRPMI2)
# join two RPMI dataset together
grVSM12 <- rbind(grVSM1, grVSM2)

# calculate average amount technical and biologial replicates
grRPMI_m <- aggregate(grRPMI12[c("lgr")], grRPMI12[c("Strain", "Place", "Replicate", "Pos")], mean)
grVSM_m <- aggregate(grVSM12[c("lgr")], grVSM12[c("Strain", "Place", "Replicate", "Pos")], mean)

# offset V and R replicates
grRPMI_m$pch[grRPMI_m$Place=="V"] <- 2
grRPMI_m$pch[grRPMI_m$Place=="R"] <- 0

grVSM_m$pch[grVSM_m$Place=="V"] <- 2
grVSM_m$pch[grVSM_m$Place=="R"] <- 0

# manually add colour info by strain and place
grRPMI_m$cols[grRPMI_m$Strain=="YST6" & grRPMI_m$Place=="V"] <- "steelblue"
grRPMI_m$cols[grRPMI_m$Strain=="YST6" & grRPMI_m$Place=="R"] <- "skyblue"
grRPMI_m$cols[grRPMI_m$Strain=="TVY4" & grRPMI_m$Place=="V"] <- "#5E478B"
grRPMI_m$cols[grRPMI_m$Strain=="TVY4" & grRPMI_m$Place=="R"] <- "#9270DC"
grRPMI_m$cols[grRPMI_m$Strain=="TVY10" & grRPMI_m$Place=="V"] <- "#EC3194"
grRPMI_m$cols[grRPMI_m$Strain=="TVY10" & grRPMI_m$Place=="R"] <- "#E37AB3"
grRPMI_m$cols[grRPMI_m$Strain=="YST7" & grRPMI_m$Place=="V"] <- "#b8860b"
grRPMI_m$cols[grRPMI_m$Strain=="YST7" & grRPMI_m$Place=="R"] <- "#EEDD82"

grVSM_m$cols[grVSM_m$Strain=="YST6" & grVSM_m$Place=="V"] <- "steelblue"
grVSM_m$cols[grVSM_m$Strain=="YST6" & grVSM_m$Place=="R"] <- "skyblue"
grVSM_m$cols[grVSM_m$Strain=="TVY4" & grVSM_m$Place=="V"] <- "#5E478B"
grVSM_m$cols[grVSM_m$Strain=="TVY4" & grVSM_m$Place=="R"] <- "#9270DC"
grVSM_m$cols[grVSM_m$Strain=="TVY10" & grVSM_m$Place=="V"] <- "#EC3194"
grVSM_m$cols[grVSM_m$Strain=="TVY10" & grVSM_m$Place=="R"] <- "#E37AB3"
grVSM_m$cols[grVSM_m$Strain=="YST7" & grVSM_m$Place=="V"] <- "#b8860b"
grVSM_m$cols[grVSM_m$Strain=="YST7" & grVSM_m$Place=="R"] <- "#EEDD82"

######################
#stats
######################

# subset datasets by strain and site
grRPMI_YST6V <- subset(grRPMI_m, Strain == "YST6" & Place == "V")
grRPMI_YST6R <- subset(grRPMI_m, Strain == "YST6" & Place == "R")
grRPMI_YST7V <- subset(grRPMI_m, Strain == "YST7" & Place == "V")
grRPMI_YST7R <- subset(grRPMI_m, Strain == "YST7" & Place == "R")
grRPMI_TVY4V <- subset(grRPMI_m, Strain == "TVY4" & Place == "V")
grRPMI_TVY4R <- subset(grRPMI_m, Strain == "TVY4" & Place == "R")
grRPMI_TVY10V <- subset(grRPMI_m, Strain == "TVY10" & Place == "V")
grRPMI_TVY10R <- subset(grRPMI_m, Strain == "TVY10" & Place == "R")

grVSM_YST6V <- subset(grVSM_m, Strain == "YST6" & Place == "V")
grVSM_YST6R <- subset(grVSM_m, Strain == "YST6" & Place == "R")
grVSM_YST7V <- subset(grVSM_m, Strain == "YST7" & Place == "V")
grVSM_YST7R <- subset(grVSM_m, Strain == "YST7" & Place == "R")
grVSM_TVY4V <- subset(grVSM_m, Strain == "TVY4" & Place == "V")
grVSM_TVY4R <- subset(grVSM_m, Strain == "TVY4" & Place == "R")
grVSM_TVY10V <- subset(grVSM_m, Strain == "TVY10" & Place == "V")
grVSM_TVY10R <- subset(grVSM_m, Strain == "TVY10" & Place == "R")


# t-tests comparing each V x R
t.test(grRPMI_YST6V$lgr, grRPMI_YST6R$lgr) #t = -0.27735, df = 21.82, p-value = 0.7841
t.test(grRPMI_TVY4V$lgr, grRPMI_TVY4R$lgr) #t = -2.5017, df = 16.142, p-value = 0.02348, R higher by 0.03
t.test(grRPMI_TVY10V$lgr, grRPMI_TVY10R$lgr) #t = 0.51907, df = 17.78, p-value = 0.6101
t.test(grRPMI_YST7V$lgr, grRPMI_YST7R$lgr) #t = -1.1687, df = 19.718, p-value = 0.2565


t.test(grVSM_YST6V$lgr, grVSM_YST6R$lgr) #t = 0.60441, df = 21.068, p-value = 0.552
t.test(grVSM_TVY4V$lgr, grVSM_TVY4R$lgr) #t = -2.3624, df = 18.125, p-value = 0.02953, R higher by 0.03
t.test(grVSM_TVY10V$lgr, grVSM_TVY10R$lgr) # t = -1.6739, df = 19.476, p-value = 0.1101
t.test(grVSM_YST7V$lgr, grVSM_YST7R$lgr) #t = -3.6258, df = 18.773, p-value = 0.001828, R higher by 0.01

rosnerTest(c(grRPMI_YST6V$lgr, grRPMI_YST6R$lgr), k=1) #0
rosnerTest(c(grRPMI_YST7V$lgr, grRPMI_YST7R$lgr), k=1) #1
rosnerTest(c(grRPMI_YST7V$lgr, grRPMI_YST7R$lgr), k=2) #1 = 0.43, obs 22 (R rep high)
rosnerTest(grRPMI_TVY4V$lgr, k=1) #0
rosnerTest(grRPMI_TVY4R$lgr, k=1) #1
rosnerTest(grRPMI_TVY4R$lgr, k=2) #1 = 0.395, obs 10 (R)
rosnerTest(c(grRPMI_TVY10V$lgr, grRPMI_TVY10R$lgr), k=1) #0

rosnerTest(c(grVSM_YST6V$lgr, grVSM_YST6R$lgr), k=1)
rosnerTest(c(grVSM_YST6V$lgr, grVSM_YST6R$lgr), k=2) #1 =0.43, obs 17 (R rep high)
rosnerTest(grVSM_YST7V$lgr, k=1)
rosnerTest(grVSM_YST7R$lgr, k=1)
rosnerTest(grVSM_TVY4V$lgr, k=1)
rosnerTest(grVSM_TVY4R$lgr, k=1)
rosnerTest(c(grVSM_TVY10V$lgr, grVSM_TVY10R$lgr), k=1)
rosnerTest(c(grVSM_TVY10V$lgr, grVSM_TVY10R$lgr), k=2) #1 = 0.285, obs 22 (R rep high)

#########################
#One plot
#########################

pdf("figures_tables/FigureX-GC_all.pdf", width=4, height=6, font = "Times")
par(mfrow=c(2, 1), mar=c(1, 1, 1, 1), oma=c(4, 4, 1, 1))
plot(jitter(as.numeric(grRPMI_m$Pos)), grRPMI_m$lgr, col=grRPMI_m$col, pch=grRPMI_m$pch, xlim=c(0.75, 4.35), ylim=c(0, 0.5), cex=1.2, yaxt="n", xaxt="n", ylab="Growth rate (/h)", xlab="")
axis(2, las=2)
axis(1, at=c(1.1, 2.1, 3.1, 4.1), labels=FALSE)
mtext("A) RPMI", side = 3, outer=FALSE, adj=0.01)

plot(jitter(grVSM_m$Pos), grVSM_m$lgr, col=grVSM_m$col, pch = grVSM_m$pch, xlim=c(0.75, 4.35), ylim=c(0, 0.5), cex=1.2, yaxt="n", xaxt="n", ylab="Growth rate (/h)", xlab="")
axis(2, las=2)
axis(1, at=c(1.1, 2.1, 3.1, 4.1), labels=c("YST6", "YST7", "TVY4", "TVY10"))
legend("bottomright", legend=c("vaginal", "rectal"), pch=c(2, 0), cex=0.9, bty="n")
mtext("B) VSM", side = 3, outer=FALSE, adj=0.01)
mtext("Growth rate (/h)", side = 2, outer=TRUE, line = 2.5, cex=1.5)
dev.off()
