# Steps for calling variants
1. Trim fastq files - Trimmomatic
2. Reference based mapping - Bowtie2
3. Add or replace read groups in bam file - Picard
	- This tool enables the user to replace all read groups in the INPUT file with a single new read group and assign all reads to this read group in the OUTPUT BAM file
5. Mark duplicates - Picard
	- Identifies duplicate reads
	- This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA
	- Bam file must be indexed for the following GATK command
7. Merge bam files - samtools merge
	- this command is necessary to run Freebayes
	- Freebayes requires one bam file, but it can be merged
9. Run Freebayes to call variants on bam files - Freebayes
10. Filter SNPs and Indels from VCF file output from Freebayes
	- gatk SelectVariants
	- gatk Variant Filtration
12. Recalibrate base quality score - gatk BaseRecalibrator
13. Apply BQSR reports into the bam files - gatk ApplyBQSR
14. Merge recalibrated bam files - samtools merge
15. Run Freebayes to call variants on recalibrated bam files
16. Filter vcf file from freebayes for SNPs only - vcffilter
17. Compressing and indexing the snp only vcf & variant normalization of the snp only vcf & decompressing the vcf
	- compress vcf - bgzip
	- index vcf - tabix
	- normalize snp vcf - bcftools norm
	- decompress the normalized vcf - bgzip
18. Filter normalized vcf file and create SNP alignment fasta file - used custom made scripts
	- vcflib_pipeline_HPC2.sh script
	- uses vcf_fa_extractor.py to filter and make fasta file from vcf