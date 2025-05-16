#!/bin/bash
#SBATCH --account=def-acgerste
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=0-42:02:00
#SBATCH --mem=375G
#SBATCH --mail-user=adamubua@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=Clean_gVCF
#SBATCH --output=%x-%j.out


#Load modules


module load nixpkgs/16.09
module load StdEnv/2020 bwa picard gatk




#Define input and output folders
output="/home/abdul/scratch/C.albicans/data_out/trimmed_reads"
input="/home/abdul/scratch/C.albicans/data_in/raw_files"


#Create directories for each sample
cat filenames.txt | parallel 'mkdir 
/home/abdul/scratch/C.albicans/data_out/clean_alignment/{}'



#Add read groups
cat filenames.txt | parallel 'time java -jar $EBROOTPICARD/picard.jar 
AddOrReplaceReadGroups 
I=/home/abdul/scratch/C.albicans/data_out/aligned_reads/{}/{}.sorted.bam 
R=/home/abdul/scratch/C.albicans/data_in/reference/Candida_albicans.fa.gz 
O=/home/abdul/scratch/C.albicans/data_out/clean_alignment/{}/{}.sorted_RG.bam 
RGID=SeqCoast RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM={}'


#Mark potential reads duplicates


cat filenames.txt | parallel 'time java -jar $EBROOTPICARD/picard.jar 
MarkDuplicates 
I=/home/abdul/scratch/C.albicans/data_out/clean_alignment/{}/{}.sorted_RG.bam 
O=/home/abdul/scratch/C.albicans/data_out/clean_alignment/{}/{}.sorted_RG_dedup.bam 
R=/home/abdul/scratch/C.albicans/data_in/reference/Candida_albicans.fa.gz 
M=/home/abdul/scratch/C.albicans/data_out/clean_alignment/{}/{}.sorted_RG_dedup.txt 
CREATE_INDEX=true READ_NAME_REGEX=null'


cat filenames.txt | parallel 'time java -jar $EBROOTPICARD/picard.jar 
ValidateSamFile 
I=/home/abdul/scratch/C.albicans/data_out/clean_alignment/{}/{}.sorted_RG_dedup.bam 
R=/home/abdul/scratch/C.albicans/data_in/reference/Candida_albicans.fa.gz 
IGNORE=INVALID_TAG_NM'


#echo " Mark duplicates done!"


#Correct possible info differences in the alignmed paired end reads. 
Readmore on this step




#cat filenames.txt | parallel 'time java -jar $EBROOTPICARD/picard.jar 
FixMateInformation 
I=/home/abdul/scratch/C.albicans/data_out/clean_alignment/{}/{}.sorted_RG_dedup.bam 
R=/home/abdul/scratch/C.albicans/data_in/reference/Candida_albicans.fa.gz 
O=/home/abdul/scratch/C.albicans/data_out/clean_alignment/{}/{}.sorted_RG_dedup_Fixmate.bam 
ADD_MATE_CIGAR=true CREATE_INDEX=true'




cat filenames.txt | parallel 'time java -jar $EBROOTPICARD/picard.jar 
ValidateSamFile 
I=/home/abdul/scratch/C.albicans/data_out/clean_alignment/{}/{}.sorted_RG_dedup_Fixmate.bam 
R=/home/abdul/scratch/C.albicans/data_in/reference/Candida_albicans.fa.gz 
IGNORE=INVALID_TAG_NM'



#echo " Fixmate information done!"


cat filenames.txt | parallel 'gatk BaseRecalibrator -I 
/home/abdul/scratch/C.albicans/data_out/clean_alignment/{}/{}.sorted_RG_dedup_Fixmate.bam 
-R /home/abdul/scratch/C.albicans/data_in/reference/Candida_albicans.fa.gz 
--known-sites 
/home/abdul/scratch/C.albicans/data_in/reference/A22_Jones_PMID_15123810_Polymorphisms_hapA.vcf.gz 
-O 
/home/abdul/scratch/C.albicans/data_out/clean_alignment/{}/{}.recal_data.table'


cat filenames.txt | parallel 'gatk ApplyBQSR -R 
/home/abdul/scratch/C.albicans/data_in/reference/Candida_albicans.fa.gz -I 
/home/abdul/scratch/C.albicans/data_out/clean_alignment/{}/{}.sorted_RG_dedup_Fixmate.bam 
--bqsr-recal-file 
/home/abdul/scratch/C.albicans/data_out/clean_alignment/{}/{}.recal_data.table 
-O 
/home/abdul/scratch/C.albicans/data_out/clean_alignment/{}/{}.sorted_RG_dedup_Fixmate_bqsr.bam'


#Run Haplotypecaller in g.vcf mode


cat filenames.txt | parallel 'gatk --java-options "-Xmx4g" HaplotypeCaller 
--native-pair-hmm-threads $SLURM_CPUS_PER_TASK  -R 
/home/abdul/scratch/C.albicans/data_in/reference/Candida_albicans.fa.gz 
-ploidy 2 -I 
/home/abdul/scratch/C.albicans/data_out/clean_alignment/bam_files/{}.sorted_RG_dedup_Fixmate_bqsr.bam 
-O 
/home/abdul/scratch/C.albicans/data_out/varaint_calling/{}/{}.raw.snps.indels.g.vcf.gz 
-ERC GVCF'


#echo " HaplotypeCaller done!"
