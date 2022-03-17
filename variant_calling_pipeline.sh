#!/bin/bash
#SBATCH --job-name=bam_pipeline
#SBATCH --output=%parallel_bam_pipeline_%j.out
#SBATCH --error=bam_pipeline.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tpaisie@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --mem=80gb
#SBATCH --time=96:00:00
#SBATCH --account=salemi
#SBATCH --qos=salemi

pwd; hostname; date

# module load intel 
# module load gcc
# module load parallel
# module load fastqc
# module load trimmomatic
# module load bowtie2
module load samtools
module load picard
module load gatk/3.8

export _JAVA_OPTIONS="-Xmx50g"

# REF=/ufrc/data/reference/bowtie2/v_cholerae_o1_2010el_1786
# REF=/ufrc/salemi/tpaisie/vvulnificus/ref_seqs/M06-24.fasta
#REF=/ufrc/data/reference/bowtie2/vibrChol1

# ls -1 | parallel 'fastqc -f fastq {}_1.fastq.gz {}_2.fastq.gz' ::: $(ls *_*.fastq.gz | rev | cut -c 12- | rev | uniq)
# parallel 'trimmomatic PE {}_1.fastq.gz {}_2.fastq.gz {}_pair_1.fastq.gz U_{}_1.fastq.gz {}_pair_2.fastq.gz U_{}_2.fastq.gz ILLUMINACLIP:/apps/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36' ::: $(ls *_*.fastq.gz | rev | cut -c 12- | rev | uniq)
# ls -1 | parallel 'fastqc -f fastq {}_pair_1.fastq.gz {}_pair_2.fastq.gz' ::: $(ls *_pair_*.fastq.gz | rev | cut -c 17- | rev | uniq)

#bowtie2
parallel 'bowtie2 -p 8 -x /blue/salemi/tpaisie/goss_project/ref_seq/NZ_CP018728 -1 {}_1.fastq.gz -2 {}_2.fastq.gz -I 0 -X 1200 --un-conc-gz {}_unmapped | samtools view -bS - | samtools sort -O bam -T tmp_{}.tmp > {}_sorted.bam' ::: $(ls *_*.fastq.gz | rev | cut -c 12- | rev | uniq)

# # Add or replace read groups picard command
# # Replace read groups in a BAM file. 
# # This tool enables the user to replace all read groups in the INPUT file with a single new read group and assign all reads to this read group in the OUTPUT BAM file.
# parallel 'picard AddOrReplaceReadGroups I={}_sorted.bam O={}_re.bam RGID=ID_{} RGLB=LB_{} RGPL=ILLUMINA RGPU=ILLUMINA RGSM=SM_{}' ::: $(ls *_sorted.bam | rev | cut -c 12- | rev | uniq)

# parallel 'picard AddOrReplaceReadGroups -I {}_sorted.bam -O {}_re.bam -RGID ID_{} -RGLB LB_{} -RGPL ILLUMINA -RGPU ILLUMINA -RGSM SM_{}' ::: $(ls *_sorted.bam | rev | cut -c 12- | rev | uniq)


# # # mark duplicates - Identifies duplicate reads
# # This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. 
# parallel 'picard MarkDuplicates I={}_re.bam O={}_nodups.bam M=marked_dup_metrics_{}.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true READ_NAME_REGEX="[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*" OPTICAL_DUPLICATE_PIXEL_DISTANCE=100' ::: $(ls *_re.bam | rev | cut -c 8- | rev | uniq)

# for i in $(ls *_re.bam | rev | cut -c 8- | rev | uniq)
# 	do
# 		picard MarkDuplicates I=${i}_re.bam O=${i}_nodups.bam M=marked_dup_metrics_${i}.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true READ_NAME_REGEX="[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*" OPTICAL_DUPLICATE_PIXEL_DISTANCE=100
# 	done

# index the output bam file from MarkDuplicates
# Must be indexed for the following GATK command

#################################################################################################
# merge together bam tools for 1st round of variant calling with freebayes

FILE=$1
OUTPUT=$2

samtools merge -b $FILE $OUTPUT


FILE1=$1
FILE2=$2

freebayes -L ${FILE1} -f $REF -T 0.001 -p 1 -i -X -n 0 -E 3 --min-repeat-size 5 -m 1 -q 20 -R 0 -Y 0 -e 1000 -F 0.5 -C 2 -3 0 -G 1 -! 0 > ${FILE2}


export _JAVA_OPTIONS="-Xmx10g"

# VCF=$1

# # # # step 2 - extract SNPs and Indels from each vcf file
# # # # extracts SNPS
# gatk SelectVariants -R $ref -V ${VCF} --select-type-to-include SNP -O snps_${VCF}

# gatk SelectVariants -R $ref -V ${VCF} --select-type-to-include INDEL -O indels_${VCF}

# # # # # # # step 3 - filter the SNP and Indel files
# # # # # # # filters the SNP vcf
# gatk VariantFiltration -R $ref -V snps_${VCF} --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filter-name "basic_snp_filter" -O filtered_snps_${VCF}

# gatk VariantFiltration -R $ref -V indels_${VCF} --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filter-name "basic_indel_filter" -O filtered_indels_${VCF}

###########################################################################################################################

# BQSR #1
# parallel 'gatk BaseRecalibrator -R /blue/salemi/tpaisie/cholera/ref_seq/vibrChol1.fa -I {}.bam --known-sites filtered_snps_new_drc_uncal_snps_051121.vcf --known-sites filtered_indels_new_drc_uncal_snps_051121.vcf -O {}.table' ::: $(ls *.bam | rev | cut -c 5- | rev | uniq)

# for i in $(ls *.bam | rev | cut -c 5- | rev | uniq)
# 	do
# 		gatk BaseRecalibrator -R $ref -I ${i}.bam --known-sites filtered_snps_${VCF} --known-sites filtered_indels_${VCF} -O ${i}.table
# 	done

# Analyze the BQSR reports from the base recalibration steps (steps 4 & 5)
# parallel 'gatk ApplyBQSR -R /ufrc/salemi/tpaisie/cholera/ref_seq/v_cholerae_o1_2010el_1786.fa -I {}.bam --bqsr-recal-file {}.table -O recal_{}.bam' ::: $(ls *.bam | rev | cut -c 5- | rev | uniq)

# gatk ApplyBQSR -R $REF -I ${FILE}.bam --bqsr-recal-file ${FILE}.table -O recal_${FILE}.bam

# for i in $(ls *.bam | rev | cut -c 5- | rev | uniq)
# 	do
# 		gatk ApplyBQSR -R $ref -I ${i}.bam --bqsr-recal-file ${i}.table -O recal_${i}.bam
# 	done


# LIST=$1
 # #calling variants on recalibrated bam files (will have only one vcf file as the output)
#freebayes -L ${LIST}.txt -v ${LIST}.vcf -f $REF 0.001 -p 1 -i -X -n 0 -E 3 --min-repeat-size 5 -m 1 -q 20 -R 0 -Y 0 -e 1000 -F 0.5 -C 2 -3 0 -G 1 -! 0


# # # step 9 - filter vcf file from freebayes for SNPs only
# vcffilter -f "TYPE = snp" ${LIST}.vcf > snps_${LIST}.vcf


# # # # step 10 - compressing and indexing the snp only vcf file & variant normalization of the snp only vcf file & decompressing the vcf file

# # # #compress vcf
# bgzip snps_${LIST}.vcf

# # #index vcf
# tabix -p vcf snps_${LIST}.vcf.gz

# # # normalizing the variant vcf file
# bcftools norm -f $ref -o norm_${LIST}.vcf.gz snps_${LIST}.vcf.gz

# # #decompresses the output variant normalized file
# bgzip -d norm_${LIST}.vcf.gz


# # step 11 - filter normalized vcf file and create SNP alignment fasta file
# bash vcflib_pipeline_HPC2.sh norm_${LIST}.vcf



