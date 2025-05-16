#libraries
library(tidyverse)
library(ggforce)

# load in snp_diff files
TVY4_snp_diff_df <-read_csv(here("genomics", "pairwise_comparisons","pairwise_df", "TVY4_snp_diff_df.csv"))
TVY10_snp_diff_df <-read_csv(here("genomics","pairwise_comparisons", "pairwise_df", "TVY10_snp_diff_df.csv"))
YST7_snp_diff_df <-read_csv(here("genomics", "pairwise_comparisons", "pairwise_df", "YST7_snp_diff_df.csv"))
YST6_snp_diff_df <-read_csv(here("genomics", "pairwise_comparisons", "pairwise_df", "YST6_snp_diff_df.csv"))

# Visualize the data frame


####Plot
#colours
#YST7
rgb(254,126, 121, maxColorValue = 255) #"#FE7E79"
rgb(181,23, 0, maxColorValue = 255) #"#B51700"
rgb(217.5, 74.5, 60.5, maxColorValue = 255) # "#D94A3C"

#png("~/Nextcloud/THRIVE/Research/THRIVE_yeast-VR-cloned/genomics/figures/240606_YST6_SNP_diff.png")
pdf(here("genomics", "figures", "240607_YST7_SNP_diff.pdf"), width = 3.1, height = 3.9)
#ggplot(YST6_snp_diff_df, aes(x=Compare, y=Pairwise_SNP_diff, fill=Compare, color = Compare))  +
ggplot(YST7_snp_diff_df, aes(x=Compare, y=Pairwise_SNP_diff, color = Compare))  +
  geom_violin(trim=F,alpha = 0.2) +
  scale_color_manual(values=c("#EEDD82","goldenrod","#b8860b")) +
  geom_sina(size = 0.75) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2) +
  stat_summary(fun = "mean",geom = "crossbar", width = 0.5, colour = "Black") +
  ylab("Pairwise SNP Difference")  +
  xlab("Paired groups") +
  theme(legend.title = element_blank()) +
  ylim(0, 2500) +
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())  +
  theme(legend.position = "none")
dev.off()

#TVY4
rgb(146,112, 220, maxColorValue = 255) #"#9270DC"
rgb(94,71, 139, maxColorValue = 255) #"#5E478B"
rgb(120, 91.5, 179.5, maxColorValue = 255) # "#785BB3"

pdf(here("genomics", "figures", "240607_TVY4_SNP_diff.pdf"), width = 3.1, height = 3.9)
#ggplot(YST6_snp_diff_df, aes(x=Compare, y=Pairwise_SNP_diff, fill=Compare, color = Compare))  +
ggplot(TVY4_snp_diff_df, aes(x=Compare, y=Pairwise_SNP_diff, color = Compare))  +
  geom_violin(trim=F,alpha = 0.2) +
  scale_color_manual(values=c("#9270DC", "#785BB3", "#5E478B")) +
  geom_sina(size = 0.75) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2) +
  stat_summary(fun = "mean",geom = "crossbar", width = 0.5, colour = "Black") +
  ylab("Pairwise SNP Difference")  +
  xlab("Paired groups") +
  theme(legend.title = element_blank()) +
  ylim(0, 2500) +
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())  +
  theme(legend.position = "none")
dev.off()

#TVY10
rgb(227,122, 179, maxColorValue = 255) #"#E37AB3"
rgb(236,49, 148, maxColorValue = 255) #"#EC3194"
rgb(231.5, 74.5, 163.5, maxColorValue = 255) # "#E74AA3"

pdf(here("genomics", "figures", "240607_TVY10_SNP_diff.pdf"), width = 3.1, height = 3.9)
ggplot(TVY10_snp_diff_df, aes(x=Compare, y=Pairwise_SNP_diff, color = Compare))  +
  geom_violin(trim=F,alpha = 0.2) +
  scale_color_manual(values=c("#E37AB3", "#E74AA3", "#EC3194")) +
  geom_sina(size = 0.75) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2) +
  stat_summary(fun = "mean",geom = "crossbar", width = 0.5, colour = "Black") +
  ylab("Pairwise SNP Difference")  +
  xlab("Paired groups") +
  theme(legend.title = element_blank()) +
  ylim(0, 2500) +
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())  +
  theme(legend.position = "none")
dev.off()

#YST6
rgb(129,197,225, maxColorValue = 255) #"#81C5E1"
rgb(75,133, 182, maxColorValue = 255) #"#4B85B6"
rgb(102, 165, 203.5, maxColorValue = 255) # "#66A5CB"

pdf(here("genomics", "figures", "240607_YST6_SNP_diff.pdf"), width = 3.1, height = 3.9)
ggplot(YST6_snp_diff_df, aes(x=Compare, y=Pairwise_SNP_diff, color = Compare))  +
  geom_violin(trim=F,alpha = 0.2) +
  scale_color_manual(values=c("#81C5E1", "#66A5CB", "#4B85B6")) +
  geom_sina(size = 0.75) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2) +
  stat_summary(fun = "mean",geom = "crossbar", width = 0.5, colour = "Black") +
  ylab("Pairwise SNP Difference")  +
  xlab("Paired groups") +
  theme(legend.title = element_blank()) +
  ylim(0, 2500) +
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())  +
  theme(legend.position = "none")
dev.off()



######Combined########



### STATISTICS
TVY4aov <- aov(Pairwise_SNP_diff~Compare, data = TVY4_snp_diff_df)
summary(TVY4aov)
#Df Sum Sq Mean Sq F value   Pr(>F)
#Compare       2  57672   28836   8.599 0.000244 ***
#Residuals   250 838327    3353

TukeyHSD(TVY4aov)
aggregate(Pairwise_SNP_diff~Compare, data=TVY4_snp_diff_df, mean)
# Compare Complement_SNP_Sum
# 1 TVY4R-TVY4R           576.1091
# 2 TVY4R-TVY4V           601.9167
# 3 TVY4V-TVY4V           619.8636


#TVY4R-TVY4V-TVY4R-TVY4R 25.80758  3.895241 47.71991 0.0162175
#TVY4V-TVY4V-TVY4R-TVY4R 43.75455 18.827228 68.68186 0.0001411
#TVY4V-TVY4V-TVY4R-TVY4V 17.94697 -2.636085 38.53002 0.1013171


summary(aov(Pairwise_SNP_diff~Compare, data = TVY10_snp_diff_df))
#             Df   Sum Sq Mean Sq F value Pr(>F)
#Compare       2   395702  197851    1.66  0.192
#Residuals   250 29799243  119197

summary(aov(Pairwise_SNP_diff~Compare, data = YST7_snp_diff_df))
# Df  Sum Sq Mean Sq F value Pr(>F)
# Compare       2   54476   27238   2.164  0.117
# Residuals   250 3146688   12587

summary(aov(Pairwise_SNP_diff~Compare, data = YST6_snp_diff_df))
# Df  Sum Sq Mean Sq F value Pr(>F)
# Compare       2   11340    5670   0.785  0.457
# Residuals   273 1972985    7227

