library(here)
library(data.table)

#Read data
Coverage <- fread(here("genomics", "sequencing_coverage", "C.alb_Ro_coverage.txt"), sep="\t", header=T)

#Look at YST7 chromosome 3R
YST7_chr3 <- subset(Coverage, Chr == "Ca22chr3A_C_albicans_SC5314")[,c(1:2, 19:26)]
YST7_chr3_norm <- subset(Coverage, Chr == "Ca22chr3A_C_albicans_SC5314")[1:1400000,c(1:2, 19:26)]
YST7_chr3_sub <- subset(Coverage, Chr == "Ca22chr3A_C_albicans_SC5314")[1400000:1460000,c(1:2, 19:26)]
YST7_chr3_sub_scaled <- data.frame(YST7_chr3_sub[,1:2], YST7_chr3_sub[,3:10]/median(YST7_chr3_norm$YST7R20))

plot(YST7_chr3_sub_scaled$locus, YST7_chr3_sub_scaled$YST7R2, cex=0.1, xlab="position", ylab="Relative read depth", yaxt="n")
axis(2, las=2)
abline(h=1, col="red")

subCrazy <- subset(YST7_chr3_sub_scaled, YST7R2 > 2.5) # where does it jump? Consistent in all
#1410255-1410472 This is ALS7, consistent in all (>3)
#1447887-1448079 This is in ALS8 (>3)
#1449558-1456181 (>2.5)

YST7_chr3_subMinor <- subset(Coverage, Chr == "Ca22chr3A_C_albicans_SC5314")[1451600:1460000,c(1:2, 19:26)]
plot(YST7_chr3_subMinor$locus, YST7_chr3_subMinor$YST7V12, cex=0.1)


# turn into data frame
cov.df <- as.data.frame(Coverage)
colnames(cov.df)
# split into a list by chromosome
cov.ddn <- split(cov.df, cov.df$Chr)
tail(cov.ddn)
head(cov.ddn[[2]])

#Exclude MitoDNA
cov.ddn2 <- cov.ddn[names(cov.ddn) !="Ca22chrM_C_albicans_SC5314"]



# how long is each chromosome?
lapply(cov.ddn2, function(x) length(t(x[,1])))
bin <- function(ddn, line, len=1000){
  #divide the large list into a list that contains read #s for only one line (based on column number), chr, locus
  line.ddn <-  lapply(ddn, "[",  c(1, 2, line))
  #set up an empty list
  bin_cov <- list()
  # for each of the 8 chromosomes (note there are 8 elements/dataframes stored in the list)
  # loop over all chromosomes
  for (i in 1:8){
    # set up another empty list
    blanksliding <- list()
    bin1000 <- rep(seq(1,length(line.ddn[[i]]$locus),len), each=len)
    #how many extra 'positions' did we create from the above sequence
    extra <- length(bin1000)-length(line.ddn[[i]]$locus)
    #correct for the extra bins by subtracting off the extra, so we have the same number as actual positions
    bins <- bin1000[1:(length(bin1000)-extra)]
    # add as a column to the dataframe that contains data for one line and one chromosome
    line.ddn[[i]]$bin <- bins
    # now split into a new list for the focal line and focal chromosome by bin
    blanksliding <- split(line.ddn[[i]][,3], line.ddn[[i]]$bin)
    # now calculate the coverage of each bin/sliding window as the sum divided by the length and store it in a list; this contains the average # of reads for each of our bins, whicgh contains 5000 positions for the focal line and focal chromsoome
    bin_cov[[i]] <- sapply(blanksliding, function(x) sum(x)/length(x))
  }
  names(bin_cov) <- names(ddn)
  return(bin_cov)
}

head(cov.ddn)
#colnames(cov.ddn)
# then run the function for each line
TVY10R11	<-	bin(cov.ddn2,3)
TVY10R12	<-	bin(cov.ddn2,4)
TVY10R13	<-	bin(cov.ddn2,5)
TVY10R14	<-	bin(cov.ddn2,6)
TVY10V10	<-	bin(cov.ddn2,7)
TVY10V11	<-	bin(cov.ddn2,8)
TVY10V12	<-	bin(cov.ddn2,9)
TVY10V13	<-	bin(cov.ddn2,10)
TVY04R01	<-	bin(cov.ddn2,11)
TVY04R02	<-	bin(cov.ddn2,12)
TVY04R04	<-	bin(cov.ddn2,13)
TVY04R05	<-	bin(cov.ddn2,14)
TVY04V02	<-	bin(cov.ddn2,15)
TVY04V03	<-	bin(cov.ddn2,16)
TVY04V04	<-	bin(cov.ddn2,17)
TVY04V05	<-	bin(cov.ddn2,18)
YST07R20	<-	bin(cov.ddn2,19)
YST07R02	<-	bin(cov.ddn2,20)
YST07R03	<-	bin(cov.ddn2,21)
YST07R06	<-	bin(cov.ddn2,22)
YST07V12	<-	bin(cov.ddn2,23)
YST07V02	<-	bin(cov.ddn2,24)
YST07V05	<-	bin(cov.ddn2,25)
YST07V08	<-	bin(cov.ddn2,26)
#CEC3627
SRR6669882<-	bin(cov.ddn2,27)
#CEC4883
SRR6670022<-	bin(cov.ddn2,28)
#######
#The plot function
#The plot function. The variable that affect the ploidy plot is the stand. I also tweek that for each samples till I get it right.
SWline <- function(bin){
  nm <-deparse(substitute(bin))
  # calculate the number of positions in the first chromosome
  x0 <- length(bin[[1]])
  #stand <- mean(unlist(bin[[1]]))
  stand <- mean(unlist(bin[[1]]))
  #stand <- mean(subset(Flu, coverage<500)$coverage)
  plot(1:    x0,  2*unlist(bin[[1]])/stand,  type="l", col="red", xlim=c(0, length(unlist(bin))), xaxt="n", yaxt="n", xlab="", ylab="", ylim=c(0, 4))
  for (i in 2:8){
    x1 <- length(bin[[i]])+x0-1
    x <- x0:x1
    stand <- median(unlist(bin[[i]]))
    #stand <- median(median(unlist(bin[[1]])),median(unlist(bin[[2]])),median(unlist(bin[[3]])),median(unlist(bin[[4]])),median(unlist(bin[[5]])),median(unlist(bin[[6]])),median(unlist(bin[[7]])),median(unlist(bin[[8]])),median(unlist(bin[[9]])),median(unlist(bin[[10]])),median(unlist(bin[[11]])),median(unlist(bin[[12]])),median(unlist(bin[[13]])))
    #stand <- mean(mean(unlist(bin[[1]])),mean(unlist(bin[[2]])),mean(unlist(bin[[3]])),mean(unlist(bin[[4]])),mean(unlist(bin[[5]])),mean(unlist(bin[[6]])),mean(unlist(bin[[7]])),mean(unlist(bin[[8]])),mean(unlist(bin[[9]])),mean(unlist(bin[[10]])),mean(unlist(bin[[11]])),mean(unlist(bin[[12]])),mean(unlist(bin[[13]])))
    points(c(x0:x1), 2*unlist(bin[[i]]/stand), type="l", col=col[i])
    x0 <- x1
  }
  mtext(nm, side=3, adj=0.01, cex=1)
  abline(h=1, lty=2)
  abline(h=2, lty=2)
  abline(h=3, lty=2)
  abline(h=4, lty=2)
  #axis(1,at=c(49,148,254,375,509,671,862,1067,1282,1511,1761,2037,2323),cex.axis=0.5,labels=FALSE)
}
col <- c(rep(c("red", "blue"), 7))

#####Plotting based on author

#pdf("plots/2302010_THRIVE_coverage_analysis.pdf", width=8, height=10, font = "Times")
par(mfrow=c(14, 2), mar=c(0.5, 1, 1,1), oma=c(1, 1, 1, 1))
SWline(TVY04R01)
axis(2,at=c(1,2,3,4),cex.axis=0.8,las = 2)
SWline(TVY04V02)
SWline(TVY04R02)
axis(2,at=c(1,2,3,4),cex.axis=0.8,las = 2)
SWline(TVY04V03)
SWline(TVY04R04)
axis(2,at=c(1,2,3,4),cex.axis=0.8,las = 2)
SWline(TVY04V04)
SWline(TVY04R05)
axis(2,at=c(1,2,3,4),cex.axis=0.8,las = 2)
SWline(TVY04V05)
SWline(YST07R02)
axis(2,at=c(1,2,3,4),cex.axis=0.8,las = 2)
SWline(YST07V02)
SWline(YST07R03)
axis(2,at=c(1,2,3,4),cex.axis=0.8,las = 2)
SWline(YST07V05)
SWline(YST07R06)
axis(2,at=c(1,2,3,4),cex.axis=0.8,las = 2)
SWline(YST07V12)
SWline(YST07R20)
axis(2,at=c(1,2,3,4),cex.axis=0.8,las = 2)
SWline(YST07V08)
SWline(SRR6669882)
axis(2,at=c(1,2,3,4),cex.axis=0.8,las = 2)
SWline(SRR6670022)
SWline(TVY10R11)
axis(2,at=c(1,2,3,4),cex.axis=0.8,las = 2)
SWline(TVY10V10)
SWline(TVY10R12)
axis(2,at=c(1,2,3,4),cex.axis=0.8,las = 2)
SWline(TVY10V11)
SWline(TVY10R13)
axis(2,at=c(1,2,3,4),cex.axis=0.8,las = 2)
SWline(TVY10V12)
SWline(TVY10R14)
axis(2,at=c(1,2,3,4),cex.axis=0.8,las = 2)
axis(1,at=c(1750,4250,6250,8000,9350,10650,11600,13000),labels=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","ChrR"),cex.axis=1,las = 1)
SWline(TVY10V13)
#axis(1,at=c(350,850,1250,1600,1870,2130,2320,2600),labels=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","ChrR"),cex.axis=1,las = 1)
axis(1,at=c(1750,4250,6250,8000,9350,10650,11600,13000),labels=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","ChrR"),cex.axis=1,las = 1)
dev.off()


pdf("genomics/figures/230220THRIVE_albicans_coverage_sample_analysis.pdf", width=8, height=10, font = "Times")
par(mfrow=c(10, 1), mar=c(0.5, 1, 1,1), oma=c(1, 1, 1, 1))
TVY04 <- TVY04R01
YST07 <- YST07V08
TVY10 <- TVY10R14
SWline(TVY04)
axis(2,at=c(1,2,3,4),cex.axis=1.5,las = 2)
SWline(YST07)
axis(2,at=c(1,2,3,4),cex.axis=1.5,las = 2)
SWline(TVY10)
axis(2,at=c(1,2,3,4),cex.axis=1.5,las = 2)
SWline(SRR6670022)
axis(2,at=c(1,2,3,4),cex.axis=1.5,las = 2)
axis(1,at=c(1750,4250,6250,8000,9350,10650,11600,13000),labels=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","ChrR"),cex.axis=1.5,las = 1)
dev.off()

#Want to look at the CNV region on chr3 in YST7 -
#somewhere between 1400001-1460001 based on previous analysis with large bin sizes
TVY10V11_chr3 <- subset(cov.df, Chr == "Ca22chr3A_C_albicans_SC5314")[,c(1, 2, 8 )]
TVY10V11_chr3_sub <- subset(TVY10V11_chr3, locus > 1440001 & locus < 1470001)

TVY4R1_chr3 <- subset(cov.df, Chr == "Ca22chr3A_C_albicans_SC5314")[,c(1, 2, 11 )]
TVY4R1_chr3_sub <- subset(TVY4R1_chr3, locus > 1440001 & locus < 1470001)

YST7V5_chr3 <- subset(cov.df, Chr == "Ca22chr3A_C_albicans_SC5314")[,c(1, 2, 25 )]
YST7V5_chr3_sub <- subset(YST7V5_chr3, locus > 1440001 & locus < 1470001)

YST07R02_chr3 <- subset(cov.df, Chr == "Ca22chr3A_C_albicans_SC5314")[,c(1, 2, 20 )]
YST07R02_chr3_sub <- subset(YST07R02_chr3, locus > 1440001 & locus < 1470001)

SRR6670022_chr3 <- subset(cov.df, Chr == "Ca22chr3A_C_albicans_SC5314")[,c(1, 2, 28)]
SRR6670022_chr3_sub <- subset(SRR6670022_chr3, locus > 1440001 & locus < 1470001)

SRR6669882_chr3 <- subset(cov.df, Chr == "Ca22chr3A_C_albicans_SC5314")[,c(1, 2, 27)]
SRR6669882_chr3_sub <- subset(SRR6669882_chr3, locus > 1440001 & locus < 1470001)


#stand for YST07R02 to get average coverage to 100:
#2/mean(YST07R02_chr3_sub$YST7R2) = 0.01395793
plot(YST07R02_chr3_sub$locus, YST07R02_chr3_sub$YST7R2*2/mean(YST07R02_chr3_sub$YST7R2), cex=0.2)
points(TVY4R1_chr3_sub$locus, TVY4R1_chr3_sub$TVY4R1*2/mean(TVY4R1_chr3_sub$TVY4R1), cex=0.2, col="red")
points(YST7V5_chr3_sub$locus, YST7V5_chr3_sub$YST7V5*2/mean(YST7V5_chr3_sub$YST7V5), cex=0.2, col="grey")
points(TVY10V11_chr3_sub$locus, TVY10V11_chr3_sub$TVY10V11*2/mean(TVY10V11_chr3_sub$TVY10V11), cex=0.2, col="orange")
points(SRR6670022_chr3_sub$locus, SRR6670022_chr3_sub$SRR6670022*2/mean(SRR6670022_chr3_sub$SRR6670022), cex=0.2, col="blue")
points(SRR6669882_chr3_sub$locus, SRR6669882_chr3_sub$SRR6669882*2/mean(SRR6669882_chr3_sub$SRR6669882), cex=0.2, col="purple")
legend("topright", legend = c("YST07R02","TVY4R1","YST7V5","TVY10V11","SRR6670022","SRR6669882"), col = c("black","red","grey","orange","blue","purple"), pch=1)


#stand for YST07R02 to get average coverage to 100:
#2/mean(YST07R02_chr3_sub$YST7R2) = 0.01395793
par(mfrow=c(3, 2), mar=c(1, 1, 2, 2), oma=c(1, 1, 1, 1))
plot(YST07R02_chr3_sub$locus, YST07R02_chr3_sub$YST7R2*2/mean(YST07R02_chr3_sub$YST7R2), cex=0.2, ylab="Copy number", xlab="", main="YST07R2", adj=0, ylim=c(0,6), xaxt = "n", yaxt = "n")
axis(2,las=2)
abline(h=2, col="blue")
plot(YST7V5_chr3_sub$locus, YST7V5_chr3_sub$YST7V5*2/mean(YST7V5_chr3_sub$YST7V5), cex=0.2, ylab="", xlab="", main="YST7V5", adj=0, ylim=c(0,6), xaxt = "n", yaxt = "n")
abline(h=2, col="blue")
axis(1,labels = F)
axis(2,labels = F)
plot(TVY4R1_chr3_sub$locus, TVY4R1_chr3_sub$TVY4R1*2/mean(TVY4R1_chr3_sub$TVY4R1), cex=0.2, ylab="Copy number", xlab="", main="TVY4R1", adj=0, ylim=c(0,6), xaxt = "n")
abline(h=2, col="blue")
plot(TVY10V11_chr3_sub$locus, TVY10V11_chr3_sub$TVY10V11*2/mean(TVY10V11_chr3_sub$TVY10V1), cex=0.2, ylab="", xlab="", main="TVY10V11", adj=0, ylim=c(0,6), xaxt = "n", yaxt = "n")
abline(h=2, col="blue")
plot(SRR6670022_chr3_sub$locus, SRR6670022_chr3_sub$SRR6670022*2/mean(SRR6670022_chr3_sub$SRR6670022), cex=0.2, ylab="Copy number", xlab="Position", main="SRR6670022", adj=0, ylim=c(0,6))
abline(h=2, col="blue")
plot(SRR6669882_chr3_sub$locus, SRR6669882_chr3_sub$SRR6669882*2/mean(SRR6669882_chr3_sub$SRR6669882), cex=0.2, ylab="", xlab="Position", main="SRR6669882", adj=0, ylim=c(0,6), yaxt = "n")
abline(h=2, col="blue")
