#!/bin/bash

#Script to call germline variants in a tumor WGS paired end reads 2 x 100bp
#Following GATK4 best practices workflow
#I did not independently write this script, I followed this tutorial: https://www.youtube.com/watch?v=iHkiQvxyr5c

if false
then
#download data

wget -P /Users/lina_01/Desktop/variant_calling/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz

wget -P /Users/lina_01/Desktop/variant_calling/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz

#download reference files
wget -P /Users/lina_01/Desktop/variant_calling/supporting_files/hg38/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
 

gunzip /Users/lina_01/Desktop/variant_calling/supporting_files/hg38/hg38.fa



#index ref - .fai file before running haplotype caller
samtools faidx /Users/lina_01/Desktop/variant_calling/supporting_files/hg38/hg38.fa


#ref dict - .dict file before running the haplotype caller
gatk CreateSequenceDictionary R=/Users/lina_01/Desktop/variant_calling/supporting_files/hg38/hg38.fa O=~/Desktop/variant_calling/supporting_files/hg38/hg38.dict 



#download known sites files for BQSR from GATK resource bundle
wget -P /Users/lina_01/Desktop/variant_calling/supporting_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P /Users/lina_01/Desktop/variant_calling/supporting_files/hg39/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

fi

#directories
ref="/Users/lina_01/Desktop/variant_calling/supporting_files/hg38/hg38.fa"
known_sites="/Users/lina_01/Desktop/variant_calling/supporting_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="/Users/lina_01/Desktop/variant_calling/aligned_reads"
reads="/Users/lina_01/Desktop/variant_calling/reads"
results="/Users/lina_01/Desktop/variant_calling/results"
data="/Users/lina_01/Desktop/variant_calling/data"

#----------
#STEP 1: QC - Run FastQC
#----------

echo "STEP 1:QC - Run fastqc"

fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/

#(check for flagged/failed stats)
#(don't worry about GC content if not failed)
#(don't trim poor bases/reduce sequence length)
#(check adapter sequences)

#NO TRIMMING REQUIRED, QUALITY LOOKS OK

#-----------
#STEP 2: Map to reference using BWA-MEM
#----------
echo "STEP 2: Map to reference using BWA-MEM"

#Use to BWA to generate index for reference
bwa index ${ref}

#BWA alignment- 4 threads - we know that the aligned reads will be missing a read group
bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam

#(took over 60 mins for BWA index, about 50 mins for BWA MEM, 32 GB RAM machine)

#-------------
#STEP 3: Mark duplicates and sort - GATK4
#-------------

echo "STEP 3: Mark duplicates and sort - GATK4"

gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam -O ${aligned_reads}/SRR062634_sorted_dedup_reads.bam
#Took 5 mins on 32 GB RAM machine

#--------------
#STEP 4: Base quality recalibration
#--------------

echo "STEP 4: Base quality recalibration"

#1. Build the model
gatk BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table
#took 13 mins

#2. Apply the model to adjust base quality scores
gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${data}/recal_data.table -O ${aligned_reads}/SRR062634_sorted_dedup_reads_BQSR.bam
#took 6 mins

#-------------
#STEP 5: Collect alignment and insert size metrics
#-------------

echo "STEP 5: Collect alignment and insert size metrics"

gatk CollectAlignmentSummaryMetrics -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_reads_BQSR.bam -O ${aligned_reads}/alignment_metrics.txt
#took 3.15 mins

gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR062634_sorted_dedup_reads_BQSR.bam OUTPUT= ${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf
#took 1.6 mins

#----------------------
#STEP 6: Call variants - gatk haplotype caller
#----------------------

echo "STEP 6: Call variants - gatk haplotype caller"

gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_reads_BQSR.bam -O ${results}/raw_variants.vcf
#took 7-8 hours

#EXTRACT SNPS AND INDELS

gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf
#very quick, less than 10 seconds










