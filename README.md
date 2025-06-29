# MMB5010 Dissertation Project: Investigating the Role of Lipopolysaccharide-Stimulated Monocytes in the Development of Immune Tolerance in Community-Acquired Pneumonia Patients

## Project Overview

This repository contains scripts, documentation, and analysis notebooks for the processing and differential expression analysis of RNA-Seq data derived from monocyte samples treated with **E. coli LPS** and **control media**. The pipeline consists of two core components:

* **Terminal-based preprocessing and alignment pipeline** (Bash + Conda environment)
* **R-based differential expression and downstream functional analysis** (DESeq2 workflow)

---

## Experimental Context

* **Samples:** 92 monocyte RNA samples (46 media, 46 LPS-treated)
* **Platform:** RStudio (v4.0.3) and TITAN server (Ubuntu 22.04.3 LTS via VPN)
* **Reference Genome:** GRCh38 (Ensembl Release 110)

---

## Terminal Pipeline (Bash Scripts)

### 1. FASTQ Merging & Preprocessing

Recombination of split, compressed FASTQ files:

```bash
for prefix in $(ls *_R1_001.fastq.gz | sed 's/_S[0-9]*_R1_001.fastq.gz//' | sort -u); do
  zcat ${prefix}_S*_R1_001.fastq.gz > ${prefix}.fastq.gz
done
```

### 2. Read Trimming & Quality Control

Trimmomatic v0.39 with Illumina adapters:

```bash
for file in *.fastq; do
  java -jar trimmomatic.jar SE -threads 12 "$file" "${file%.fastq}_trim.fastq" \
  ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
```

Quality control tools: FastQC v0.12.1, MultiQC v1.15

```bash
fastqc *
multiqc .
```

### 3. Alignment to Reference Genome

HISAT2 v2.2.1:

```bash
hisat2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa hisat2index
for file in *_trim.fastq; do
  output_prefix=$(basename "$file" _trim.fastq)
  hisat2 -x hisat2index -U "$file" -S "$output_prefix.sam" \
  --phred33 --rna-strandness RF -t -p 16
done
```

### 4. SAM/BAM Processing

SAMTools v1.9:

```bash
for file in *.sam; do
  samtools sort -o "${file%.sam}.sorted.bam" "$file"
  gzip "${file%.sam}.sorted.bam"
done
```

### 5. Gene Expression Quantification

featureCounts v2.0.6:

```bash
featureCounts -a Homo_sapiens.GRCh38.110.gtf.gz -o subread_counts.txt -T 12 -g gene_id *.sorted.bam
```

---

## Differential Expression Analysis (R + DESeq2)

### Setup

Required libraries:

```r
BiocManager::install(c("DESeq2", "AnnotationDbi", "org.Hs.eg.db", "clusterProfiler", "EnhancedVolcano", "pheatmap"))
install.packages(c("tidyverse", "ggrepel", "stringr"))
```

### Files & Input

* `LPSvsMEDIA_counts.csv`, `TNFalpha_counts.csv`
* Corresponding sample metadata CSVs

### Steps

1. Environment & data load
2. Count filtering and normalisation (size factors, VST)
3. Exploratory quality control: PCA, scree plots, heatmaps
4. Differential expression analysis using `DESeq()` and shrinkage via `ashr`
5. Visualization: MA plots, volcano plots, p-value histograms
6. Gene Ontology enrichment using `clusterProfiler`
7. Pathway and STRING analysis with ENSEMBL mapping
8. Output generation and export of tables and figures

### Output Files

Includes:

* Differential expression tables (up- and downregulated genes)
* Normalised and unnormalised count matrices
* PCA loadings
* GO enrichment results
* Figures in `.png` format
  
---

## GSEA Analysis

Manual generation of input files:

* `.gct` files: `expression_GSEA_Media_vs_LPS.txt.gct`, etc.
* `.cls` phenotype files

To run:

```bash
./gsea.sh
```

Using GSEA v4.3.3 (Linux)

---

## Notes

* All R code is included in `DGE_Analysis.Rmd`
* The notebook output is `DGE_Analysis.nb.html` and `Differential Gene Expression_Notebook.pdf`
* Manual adjustment of the working directory path is required
* Reproducibility is ensured through pinned versions and session info

---

## Author

Ioana Beatrice Neamțu — University of Malta — MMB5010 Dissertation Project

---
