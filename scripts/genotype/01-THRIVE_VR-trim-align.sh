  #!/bin/bash
#SBATCH --account=def-acgerste
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=0-24:02:00
#SBATCH --mem=0
#SBATCH --mail-user=adamubua@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=Alignement_TVY10_4
#SBATCH --output=%x-%j.out


#Load modules


module load nixpkgs/16.09
module load StdEnv/2020 trimmomatic/0.39  bwa picard


#Make output directories for each sample in the alignment and trimming 
directory

cat filenames.txt | parallel 'mkdir 
/home/abdul/scratch/C.albicans/data_out/trimmed_reads/{}'


cat filenames.txt | parallel 'mkdir 
/home/abdul/scratch/C.albicans/data_out/aligned_reads/{}'


#Run Trimmomatic

cat filenames.txt | parallel 'java -jar 
$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads $SLURM_CPUS_PER_TASK 
-trimlog /home/abdul/scratch/C.albicans/data_out/trimmed_reads/{}/{}.log 
/home/abdul/scratch/seqcoast/230117_Order_2442/renamed_raw_reads/{}_R1.fastq.gz 
/home/abdul/scratch/seqcoast/230117_Order_2442/renamed_raw_reads/{}_R2.fastq.gz 
/home/abdul/scratch/C.albicans/data_out/trimmed_reads/{}/{}_1.trimmed_PE.fastq.gz 
/home/abdul/scratch/C.albicans/data_out/trimmed_reads/{}/{}_1.trimmed_SE.fastq.gz 
/home/abdul/scratch/C.albicans/data_out/trimmed_reads/{}/{}_2.trimmed_PE.fastq.gz 
/home/abdul/scratch/C.albicans/data_out/trimmed_reads/{}/{}_2.trimmed_SE.fastq.gz 
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33'

#Run Alignment

cat filenames.txt | parallel 'time bwa mem -t $SLURM_CPUS_PER_TASK 
/home/abdul/scratch/C.albicans/data_in/reference/Candida_albicans_index 
/home/abdul/scratch/C.albicans/data_out/trimmed_reads/{}/{}_1.trimmed_PE.fastq.gz 
/home/abdul/scratch/C.albicans/data_out/trimmed_reads/{}/{}_2.trimmed_PE.fastq.gz 
-o /home/abdul/scratch/C.albicans/data_out/aligned_reads/{}/{}.sam' 


#Sort the Alignment 

cat filenames.txt | parallel 'time java -jar $EBROOTPICARD/picard.jar 
SortSam I=/home/abdul/scratch/C.albicans/data_out/aligned_reads/{}/{}.sam 
O=/home/abdul/scratch/C.albicans/data_out/aligned_reads/{}/{}.sorted.sam 
SORT_ORDER=coordinate'



## Collect alignment statistics
cat filenames.txt | parallel 'time java -jar $EBROOTPICARD/picard.jar 
CollectAlignmentSummaryMetrics 
R=/home/abdul/scratch/C.albicans/data_in/reference/Candida_albicans.fa.gz 
I=/home/abdul/scratch/C.albicans/data_out/aligned_reads/{}/{}.sorted.sam 
O=/home/abdul/scratch/C.albicans/data_out/aligned_reads/{}/{}.alignment_summary.txt'


##Convert SAM to BAM files

cat filenames.txt | parallel 'time java -jar $EBROOTPICARD/picard.jar 
SamFormatConverter 
I=/home/abdul/scratch/C.albicans/data_out/aligned_reads/{}/{}.sorted.sam 
O=/home/abdul/scratch/C.albicans/data_out/aligned_reads/{}/{}.sorted.bam 
R=/home/abdul/scratch/C.albicans/data_in/reference/Candida_albicans.fa.gz'


##Validate BAM files

cat filenames.txt | parallel 'time java -jar $EBROOTPICARD/picard.jar 
ValidateSamFile 
I=/home/abdul/scratch/C.albicans/data_out/aligned_reads/{}/{}.sorted.bam 
R=/home/abdul/scratch/C.albicans/data_in/reference/Candida_albicans.fa.gz'
