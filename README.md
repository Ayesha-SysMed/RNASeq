This repository provides a comprehensive pipeline for processing RNA-Seq data, from raw FASTQ files to differential expression analysis and pathway enrichment. The pipeline uses various bioinformatics tools like FASTQC, Trimmomatic, STAR, DESeq2, and clusterProfiler.

Requirements

Operating System: Linux (recommended)

Tools:

FASTQC - Quality control for raw reads
Trimmomatic - Adapter trimming and quality filtering
STAR - RNA-Seq read alignment
featureCounts - Read counting
DESeq2 - Differential expression analysis (R package)
clusterProfiler - Pathway enrichment analysis (R package)
biomaRt - Gene annotation (R package)

Libraries (R)

DESeq2
pheatmap
ggplot2
dplyr
tidyverse
AnnotationDbi
clusterProfiler
org.Ce.eg.db (for C. elegans data)

Pipeline Overview

This pipeline performs the following steps:

Quality Control with FASTQC: Check the quality of raw RNA-Seq reads.
Trimming with Trimmomatic: Trim adapter sequences and low-quality bases from the reads.
Read Alignment with STAR: Align the reads to a reference genome.
Counting with featureCounts: Count the number of reads aligned to each gene.
Differential Expression Analysis with DESeq2: Identify differentially expressed genes.
Gene Annotation with biomaRt: Annotate differentially expressed genes with gene symbols.
Visualization: Generate a volcano plot to visualize differential expression results.
Pathway Analysis with clusterProfiler: Perform enrichment analysis to identify enriched biological pathways.

Directory Structure
Ensure that your working directory has the following structure:


Copy
/home/usr/Desktop/Ayesha
├── SampleA
├── SampleB
├── SampleC
├── SampleD
├── SampleE
├── SampleF
├── wbcel235.fa (Reference genome)
├── wbcel235.gtf (GTF file for annotation)
├── N2 VS CX11262.csv (Count matrix)
└── result.csv (Final differential expression results)
