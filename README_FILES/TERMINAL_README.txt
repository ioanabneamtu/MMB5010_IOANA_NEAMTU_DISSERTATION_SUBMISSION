TERMINAL SCRIPT IN BASH AND DOCUMENTATION - README

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

This README file contains scripts and descriptions of the pipeline used for RNA-Seq data analysis. The pipeline processes raw FASTQ files, performs quality control, aligns reads to a reference genome, and generates gene expression count matrices.

The following steps are covered in this README file:

	1. FASTQ Merging and Preprocessing
    	2. Quality Control
    	3. Read Alignment
    	4. SAM/BAM File Handling
    	5. Gene Expression Quantification

The pipeline was run on the TITAN server terminal via a VPN connection to the University of Malta, with access provided by Prof. Joseph Borg. All operations were performed in an Ubuntu 22.04.3 LTS environment using Conda 24.5.0.

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

1. FASTQ Merging and Preprocessing

The initial set of raw sequencing data consisted of 92 FASTQ files (46 media-treated and 46 E. coli LPS-treated samples). Each file was split and compressed, then recombined using the following script:

#!/bin/bash
# Getting list of all unique prefixes
prefixes=($(ls -1 | grep -o '^[0-9]*[A-Z]*_S[0-9]*_R1_001.fastq.gz' | sed 's/_S[0-9]*_R1_001.fastq.gz//' | sort -u))

# Looping through all prefixes and zcat matching files
for prefix in "${prefixes[@]}"; do
	zcat "${prefix}"_S*_R1_001.fastq.gz > "${prefix}.fastq.gz"
done

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

2. Quality Control

Trimming of the raw reads was performed using Trimmomatic (downloadable at https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip; v0.39) with Illumina adapters. The following script was used for trimming:

#!/bin/bash

# Looping through all .fastq files and running Trimmomatic to trim reads
for file in *.fastq; do
	java -jar /DISK1/research/bscic04/software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 12 \
	"$file" "${file%.fastq}_trim.fastq" \
	ILLUMINACLIP:/DISK1/research/bscic04/software/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 \
	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

Quality control for both raw and trimmed FASTQ files was performed using FastQC (downloadable at https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip; v0.12.1), and summary reports were generated with MultiQC (v1.15):

bash

# FastQC quality control
fastqc *

# MultiQC installation and report generation
conda install multiqc==1.15
multiqc .

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

3. Read Alignment

Reads were aligned to the human reference genome GRCh38.14 using HISAT2 (realease 7/24/2020 downloadable at 	https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download; v2.2.1). The genome and annotation files were downloaded from Ensembl (Release 110):

bash

# Downloading genome and annotation patch 14
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz--2024-04-11 11:41:40-- https:/ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz


The reference genome index was built, and trimmed reads were aligned using HISAT2 with strand specificity and Phred33 quality scores:

bash

# Building HISAT2 index
/DISK1/research/bscic04/software/hisat2/hisat2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa hisat2index

# Aligning reads to the reference genome
for file in *_trim.fastq; do
	output_prefix=$(basename "$file" _trim.fastq)
	/DISK1/research/bscic04/software/hisat2/hisat2 --threads 16  -q -x /home/inea0001/inea0001/HISAT2/hisat2index \
	-U "$file" -S "${output_prefix}.sam" --phred33 --rna-strandness RF -t
done

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

4. SAM/BAM File Handling

The SAM output files were processed using SAMTools (v1.9; downloadable at https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2) for sorting and compression:

bash

#!/bin/bash

# Sorting SAM files and converting to BAM
for file in *.sam.gz; do
	SAMTools sort -o "${file%.sam.gz}.sorted.bam" "$file"
done

# Compressing BAM files
gzip *bam

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

5. Gene Expression Quantification

Gene expression was quantified using featureCounts (v2.0.6; downloadable at https://sourceforge.net/projects/subread/files/subread-2.0.6/subread-2.0.6-Linux-x86_64.tar.gz/download), which generated count matrices directly from the sorted BAM files:

bash

# Generating gene expression counts
featureCounts -a /home/ioana/…/Homo_sapiens.GRCh38.110.gtf.gz -o subread_counts -T 12 -g gene_id *.sorted.bam


