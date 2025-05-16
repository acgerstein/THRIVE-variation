
library(readr)
library(dplyr)
library(ggplot2)
library(ggplot2)
library(ggbreak)
library(here)
thrive <- read.table(here("genomics/pixy-pi","240708_num-sample_Pi_summary.txt"), sep = "\t",header = TRUE)

thrive$pop <- factor(thrive$pop, levels=c("TVY10", "TVY4", "YST7","YST6"))

pdf(here("figures_tables", "240707_pi_vrs_num-sam.pdf"), width = 3.1, height = 3.9)

thrive %>% ggplot(aes(x=num_samples, y=MeanDiversity,color=pop)) +
  geom_point(aes(color=pop), size =4) + geom_line(size =2) +
  geom_errorbar(aes(ymin=MeanDiversity-StdDev, ymax=MeanDiversity+StdDev), width=.1) +
  scale_colour_manual(values=c("#E37AB3","#785BB3","goldenrod","#4B85B6"),name = "Participant") +# scale_y_break(c(0.00062, 0.002), scales = 0.2) +
  theme_classic(base_size = 20) +  xlab("Number of samples") +
  ylab(expression("Average nucleotide diversity " (~pi ))) + #geom_vline(xintercept = 6, linetype="dotted", color = "grey", size=1) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5), labels = function(n){format(n, scientific = FALSE)},
    sec.axis = sec_axis(
      ~ . / 1000
    )
  ) +
  theme_classic(base_size = 20) +
  ggbreak::scale_y_break(c(0.000062, 0.002), scales = c(1, 1000))

