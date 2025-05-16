# Load required libraries
library(here)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(data.table)

# Functions

# Function to extract SNP positions from a BED file
# Function to extract concatenated chromosome, start, and end positions from a BED file
extract_snps <- function(bed_file) {
  # Read BED file
  snps <- fread(bed_file)
  # Extract chromosome, start, and end positions and concatenate them into a single column
  snp_positions <- paste0(snps$V1,sep = "-",snps$V2,sep = "-",snps$V3)
  return(snp_positions)
}

# List of file paths to your 24 BED files
TVY4_bed_files <- list.files(path = here("genomics", "pairwise_comparisons", "TVY4"), pattern = "*.bed", full.names = TRUE)
TVY10_bed_files <- list.files(path =here("genomics", "pairwise_comparisons", "TVY10"), pattern = "*.bed", full.names = TRUE)
YST7_bed_files <- list.files(path = here("genomics", "pairwise_comparisons", "YST7"), pattern = "*.bed", full.names = TRUE)
YST6_bed_files <- list.files(path = here("genomics", "pairwise_comparisons", "YST6"), pattern = "*.bed", full.names = TRUE)

# Function to extract SNP positions from a BED file
# Extract sample names from file names (assuming file names are in the format "sample_name.bed")
TVY4_sample_names <- gsub(".bed$", "", basename(TVY4_bed_files))
TVY10_sample_names <- gsub(".bed$", "", basename(TVY10_bed_files))
YST7_sample_names <- gsub(".bed$", "", basename(YST7_bed_files))
YST6_sample_names <- gsub(".bed$", "", basename(YST6_bed_files))

# Matrix to store pairwise complement SNP positions (i.e., putt out what is in sample A but not B)
TVY4_pairwise_complement_sum <- matrix(NA, nrow = length(TVY4_sample_names), ncol = length(TVY4_sample_names),
                              dimnames = list(TVY4_sample_names, TVY4_sample_names))

TVY10_pairwise_complement_sum <- matrix(NA, nrow = length(TVY10_sample_names), ncol = length(TVY10_sample_names),
                                       dimnames = list(TVY10_sample_names, TVY10_sample_names))


YST7_pairwise_complement_sum <- matrix(NA, nrow = length(YST7_sample_names), ncol = length(YST7_sample_names),
                                        dimnames = list(YST7_sample_names, YST7_sample_names))

YST6_pairwise_complement_sum <- matrix(NA, nrow = length(YST6_sample_names), ncol = length(YST6_sample_names),
                                       dimnames = list(YST6_sample_names, YST6_sample_names))

# Loop through each pair of BED files
for (i in 1:length(TVY4_bed_files)) {
  for (j in i:length(TVY4_bed_files)) { # To make it symmetric and avoid duplicate calculations
    # Extract SNP positions from each BED file
    snp_positions_i <- extract_snps(TVY4_bed_files[i])
    snp_positions_j <- extract_snps(TVY4_bed_files[j])
    # Determine complement SNP positions
    complement_positions <- length(setdiff(snp_positions_i, snp_positions_j)) + length(setdiff(snp_positions_j, snp_positions_i))
    # Store the sum of complement SNP positions
    TVY4_pairwise_complement_sum[i, j] <- complement_positions
    TVY4_pairwise_complement_sum[j, i] <- NA # Set symmetric cell to NA to keep only one of the symmetric cells
  }
}
#View(TVY4_pairwise_complement_sum)

for (i in 1:length(TVY10_bed_files)) {
  for (j in i:length(TVY10_bed_files)) { # To make it symmetric and avoid duplicate calculations
    # Extract SNP positions from each BED file
    snp_positions_i <- extract_snps(TVY10_bed_files[i])
    snp_positions_j <- extract_snps(TVY10_bed_files[j])
    # Determine complement SNP positions
    complement_positions <- length(setdiff(snp_positions_i, snp_positions_j)) + length(setdiff(snp_positions_j, snp_positions_i))
    # Store the sum of complement SNP positions
    TVY10_pairwise_complement_sum[i, j] <- complement_positions
    TVY10_pairwise_complement_sum[j, i] <- NA # Set symmetric cell to NA to keep only one of the symmetric cells
  }
}
#View(TVY10_pairwise_complement_sum)

for (i in 1:length(YST7_bed_files)) {
  for (j in i:length(YST7_bed_files)) { # To make it symmetric and avoid duplicate calculations
    # Extract SNP positions from each BED file
    snp_positions_i <- extract_snps(YST7_bed_files[i])
    snp_positions_j <- extract_snps(YST7_bed_files[j])
    # Determine complement SNP positions
    complement_positions <- length(setdiff(snp_positions_i, snp_positions_j)) + length(setdiff(snp_positions_j, snp_positions_i))
    # Store the sum of complement SNP positions
    YST7_pairwise_complement_sum[i, j] <- complement_positions
    YST7_pairwise_complement_sum[j, i] <- NA # Set symmetric cell to NA to keep only one of the symmetric cells
  }
}
#View(YST7_pairwise_complement_sum)

for (i in 1:length(YST6_bed_files)) {
  for (j in i:length(YST6_bed_files)) { # To make it symmetric and avoid duplicate calculations
    # Extract SNP positions from each BED file
    snp_positions_i <- extract_snps(YST6_bed_files[i])
    snp_positions_j <- extract_snps(YST6_bed_files[j])
    # Determine complement SNP positions
    complement_positions <- length(setdiff(snp_positions_i, snp_positions_j)) + length(setdiff(snp_positions_j, snp_positions_i))
    # Store the sum of complement SNP positions
    YST6_pairwise_complement_sum[i, j] <- complement_positions
    YST6_pairwise_complement_sum[j, i] <- NA # Set symmetric cell to NA to keep only one of the symmetric cells
  }
}
#View(YST6_pairwise_complement_sum)

# Convert the matrix to a data frame
TVY4_snp_diff_df <- reshape2::melt(TVY4_pairwise_complement_sum, na.rm = TRUE)
TVY10_snp_diff_df <- reshape2::melt(TVY10_pairwise_complement_sum, na.rm = TRUE)
YST7_snp_diff_df <- reshape2::melt(YST7_pairwise_complement_sum, na.rm = TRUE)
YST6_snp_diff_df <- reshape2::melt(YST6_pairwise_complement_sum, na.rm = TRUE)

# Rename columns
colnames(TVY4_snp_diff_df) <- c("Sample1", "Sample2","Pairwise_SNP_diff")
colnames(TVY10_snp_diff_df) <- c("Sample1", "Sample2","Pairwise_SNP_diff")
colnames(YST7_snp_diff_df) <- c("Sample1", "Sample2","Pairwise_SNP_diff")
colnames(YST6_snp_diff_df) <- c("Sample1", "Sample2","Pairwise_SNP_diff")

# Add column names indicating which samples were combined
#pairwise_complement_df$Samples_Combined <- paste(pairwise_complement_df$Sample1,pairwise_complement_df$Sample2, sep = "vs")
row.names(TVY4_snp_diff_df) <- c(1:length(TVY4_snp_diff_df$Sample1))
row.names(TVY10_snp_diff_df) <- c(1:length(TVY10_snp_diff_df$Sample1))
row.names(YST7_snp_diff_df) <- c(1:length(YST7_snp_diff_df$Sample1))
row.names(YST6_snp_diff_df) <- c(1:length(YST6_snp_diff_df$Sample1))

#sample_groups <- c(rep("Rectal", 66), rep("Rectal_Vaginal", 83), rep("Vaginal", 127))  # Example: 3 groups of 24 samples each
#snp_diff_df <- cbind(TVY4_snp_diff_df, Group = sample_groups)
#TVY4_snp_diff_df$log_SNP_Difference <- log(TVY4_snp_diff_df$Complement_SNP_Sum)
TVY4_snp_diff_df$S1 <- str_sub(TVY4_snp_diff_df$Sample1, 5, 5)
TVY4_snp_diff_df$S2 <- str_sub(TVY4_snp_diff_df$Sample2, 5, 5)
TVY4_snp_diff_df$Compare <- paste0(TVY4_snp_diff_df$S1,sep = "-",TVY4_snp_diff_df$S2)

TVY10_snp_diff_df$S1 <- str_sub(TVY10_snp_diff_df$Sample1, 6, 6)
TVY10_snp_diff_df$S2 <- str_sub(TVY10_snp_diff_df$Sample2, 6, 6)
TVY10_snp_diff_df$Compare <- paste0(TVY10_snp_diff_df$S1,sep = "-",TVY10_snp_diff_df$S2)

YST7_snp_diff_df$S1 <- str_sub(YST7_snp_diff_df$Sample1, 5, 5)
YST7_snp_diff_df$S2 <- str_sub(YST7_snp_diff_df$Sample2, 5, 5)
YST7_snp_diff_df$Compare <- paste0(YST7_snp_diff_df$S1,sep = "-",YST7_snp_diff_df$S2)

YST6_snp_diff_df$S1 <- str_sub(YST6_snp_diff_df$Sample1, 5, 5)
YST6_snp_diff_df$S2 <- str_sub(YST6_snp_diff_df$Sample2, 5, 5)
YST6_snp_diff_df$Compare <- paste0(YST6_snp_diff_df$S1,sep = "-",YST6_snp_diff_df$S2)

# save the files
write.csv(TVY4_snp_diff_df, here("genomics", "pairwise_comparisons", "pairwise_df", "TVY4_snp_diff_df.csv"), row.names=FALSE)
write.csv(TVY10_snp_diff_df, here("genomics", "pairwise_comparisons", "pairwise_df", "TVY10_snp_diff_df.csv"), row.names=FALSE)
write.csv(YST7_snp_diff_df, here("genomics", "pairwise_comparisons", "pairwise_df", "YST7_snp_diff_df.csv"), row.names=FALSE)
write.csv(YST6_snp_diff_df, here("genomics", "pairwise_comparisons", "pairwise_df", "YST6_snp_diff_df.csv"), row.names=FALSE)
