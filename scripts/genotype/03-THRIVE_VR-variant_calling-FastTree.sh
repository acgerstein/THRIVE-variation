#!/bin/bash
#SBATCH --account=def-acgerste
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=45:00:00
#SBATCH --mem=0
#SBATCH --mail-user=adamubua@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=Variant_call_THRIVE
#SBATCH --output=%x-%j.out


module load StdEnv/2020
module load picard gatk

output="/home/abdul/scratch/C.albicans/data_out/varaint_calling/combined_variant_calls"
ref="/home/abdul/scratch/C.albicans/data_in/reference/Candida_albicans.fa.gz"

#Combine all gvcf files

gatk CombineGVCFs -R $ref --variant THRIVE_global_gvcfs.list -O $output/230131_THRIVE_calb_combined.g.vcf.gz

#Run genotype vcf to generate final vcf

gatk --java-options "-Xmx4g" GenotypeG VCFs --dbsnp /home/abdul/scratch/C.albicans/data_in/reference/A22_Jones_PMID_15123810_Polymorphisms_hapA.vcf.gz -R $ref -V $output/230131_THRIVE_calb_combined.g.vcf.gz -O $output/230131_THRIVE_calb_joint_called_raw.vcf.gz

#Select only SNPs from the file for filtration

gatk SelectVariants \
    -V $output/230131_THRIVE_calb_joint_called_raw.vcf.gz \
    -select-type SNP \
    -R $ref \
    -O $output/230131_THRIVE_calb_joint_called_raw.snps.vcf.gz

#Filter SNPs

gatk VariantFiltration \
    -R $ref \
    -V $output/230131_THRIVE_calb_joint_called_raw.snps.vcf.gz \
    -O $output/230131_THRIVE_calb_joint_called_filtered.snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" 



###exclude mtDNA
gatk SelectVariants -V $output/230131_THRIVE_calb_joint_called_filtered.snps.vcf.gz  -XL Ca22chrM_C_albicans_SC5314 -R $ref -O $output/230131_THRIVE_calb_joint_called_filtered.nomtd.snps.vcf.gz


#Select InDel only from the file for filtration
gatk SelectVariants -V $output/230131_THRIVE_calb_joint_called_raw.vcf.gz  -R $ref -select-type INDELs -O $output/230131_THRIVE_calb_joint_called_raw.indel.vcf.gz


###Filter InDels

gatk VariantFiltration -R $ref -V $output/230131_THRIVE_calb_joint_called_raw.indel.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O $output/220805_calb_joint_called_filtered.indel.vcf.gz


###exclude mtDNA
gatk SelectVariants -V $output/230131_THRIVE_calb_joint_called_filtered.indel.vcf.gz  -XL Ca22chrM_C_albicans_SC5314 -R $ref -O $output/230131_THRIVE_calb_joint_called_filtered.nomtd.indel.vcf.gz


#Combine/Merge the SNPs and InDels filtered files

java -jar $EBROOTPICARD/picard.jar MergeVcfs \
          I=$output/230131_THRIVE_calb_joint_called_filtered.nomtd.snps.vcf.gz \
          I=output/230131_THRIVE_calb_joint_called_filtered.nomtd.indel.vcf.gz \
          O=$output/230131_THRIVE_calb_merged_filtered_variants.nomtd.vcf.gz


# Select only varaints that passed. NB: the combined variant file has varants marked to be filtered. These are normally not used in subsequent analyses but they contribute to file sizes
gatk SelectVariants -R $ref --exclude-filtered -V $output/230131_THRIVE_calb_merged_filtered_variants.nomtd.vcf.gz -O $output/230131_THRIVE_calb_merged_passed_variants.nomtd.vcf.gz

gatk CollectVariantCallingMetrics --DBSNP /home/abdul/scratch/C.albicans/data_in/reference/A22_Jones_PMID_15123810_Polymorphisms_hapA.vcf.gz -I $output/230131_THRIVE_calb_merged_passed_variants.nomtd.vcf.gz -O $output/230131_THRIVE_variant_stats

echo " Varaint Calling  done!"


#Phylogeny

python vcf2phylip/vcf2phylip.py --input $output/230131_THRIVE_calb_joint_called_filtered.nomtd.snps.vcf.gz --output-folder phylogeny --nexus --nexus-binary --fasta


./FastTreeDbl -nt -gtr phylogeny/*230131_THRIVE*.fasta > phylogeny/230131_THRIVE_calb_tree.nwk

