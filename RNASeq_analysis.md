# RNA-Seq Analysis Pipeline

This guide outlines the steps to process RNA-Seq data from raw FASTQ files to differential expression analysis and pathway enrichment using various tools like FASTQC, Trimmomatic, STAR, and DESeq2.

## Initial Setup

### Sample Information
- **Operating System:** Linux
- **Directory:** `/home/usr/Desktop/Ayesha`
- **Samples:** Paired-end reads, 12 reads total, with one folder for each sample (SampleA to SampleF)

### Tools Used
- **FASTQC:** For quality control of raw reads
- **Trimmomatic:** For trimming adapters and low-quality bases
- **STAR:** For aligning reads to the reference genome
- **featureCounts:** For counting reads mapped to genes
- **DESeq2:** For differential expression analysis
- **Pathway Analysis:** Using clusterProfiler

## Step-by-Step Instructions

### 1. Quality Control with FASTQC

#### Download and Install FASTQC
Download FASTQC from the browser.

#### Move to Directory and Unzip FASTQ Files
```bash
cd /home/usr/Desktop/Ayesha
gunzip *.fastq.gz
```

#### Run FASTQC
```bash
fastqc *.fastq --threads 6 --outdir .
```
- `*.fastq`: Specifies all FASTQ files in the directory.
- `--outdir .`: Output directory (current directory).

### 2. Trimming with Trimmomatic

#### Install Trimmomatic
```bash
sudo apt-get install trimmomatic
```

#### Find Installed Packages
```bash
dpkg -L trimmomatic
```

#### Run Trimmomatic for Each Sample
```bash
java -jar /usr/share/java/trimmomatic-0.39.jar PE A_R1.fq.gz A_R2.fq.gz Aoutput_forward_paired.fq.gz Aoutput_forward_unpaired.fq.gz Aoutput_reverse_paired.fq.gz Aoutput_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
```
- `PE`: Paired-End mode.
- `ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True`: Adapter trimming parameters.
- `LEADING:3`: Remove leading low-quality bases with Phred score < 3.
- `TRAILING:3`: Remove trailing low-quality bases with Phred score < 3.
- `MINLEN:36`: Discard reads shorter than 36 bases.

### 3. Alignment with STAR

#### Download Reference Genome and GTF Files
Download from [WormBase Parasite FTP](https://parasite.wormbase.org/ftp.html) and unzip.
```bash
gunzip *.gz
```

#### Build Index
```bash
STAR --runMode genomeGenerate --genomeDir Star_index/ --genomeFastaFiles wbcel235.fa --sjdbGTFfile wbcel235.gtf --runThreadN 6
```

#### Align Reads
```bash
STAR --runMode alignReads --genomeDir Star_index/ --outSAMtype BAM SortedByCoordinate --readFilesIn Aoutput_forward_paired.fastq Aoutput_reverse_paired.fastq --runThreadN 6
```

### 4. Counting Reads with featureCounts

#### Install featureCounts
```bash
sudo apt-get install subread
```

#### Count Reads
```bash
featureCounts -a wbcel235.gtf -o count.out -T 8 *.bam
```

### 5. Differential Expression Analysis with DESeq2

#### Load Required Libraries
```r
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyverse)
```

#### Load and Preprocess Counts Data
```r
counts <- read.delim("N2 VS CX11262.csv", header = TRUE, row.names = 1, sep = ",")
counts <- counts[which(rowSums(counts) > 50),]

conditions <- factor(c("N2", "N2", "CX11262", "CX11262"))
coldata <- data.frame(row.names = colnames(counts), condition = conditions)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
```

#### Normalization and Transformation
```r
dds <- estimateSizeFactors(dds)
exp_normalised_counts <- counts(dds, normalized = TRUE)
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)
```

#### Clustering and PCA
```r
vsd_cor <- cor(vsd_mat)
pheatmap(vsd_cor, annotation = select(coldata, condition))
plotPCA(vsd, intgroup="condition")
```

#### Run DESeq
```r
dds <- DESeq(dds)
```

#### Mean vs Variance Plot
```r
mean_counts_control <- apply(counts[,1:2], 1, mean)
variance_counts_control <- apply(counts[,1:2], 1, var)
df_control <- data.frame(mean_counts_control, variance_counts_control)
ggplot(df_control) + geom_point(aes(x=mean_counts_control, y=variance_counts_control)) + scale_y_log10() + scale_x_log10() + xlab("Mean Count per gene in control") + ylab("Variance per gene in control") + ggtitle("FOR N2")

mean_counts_treated <- apply(counts[,3:4], 1, mean)
variance_counts_treated <- apply(counts[,3:4], 1, var)
df_treated <- data.frame(mean_counts_treated, variance_counts_treated)
ggplot(df_treated) + geom_point(aes(x=mean_counts_treated, y=variance_counts_treated)) + scale_y_log10() + scale_x_log10() + xlab("Mean Count per gene in treated") + ylab("Variance per gene in treated") + ggtitle("FOR CX11262")
```

#### Dispersion Plot
```r
plotDispEsts(dds)
```

#### Wald Test Results
```r
exp_res <- results(dds, alpha = 0.05)
plotMA(exp_res, ylim=c(-8, 8))
```

#### Log Fold Change Shrinkage
```r
BiocManager::install("apeglm")
exp_res <- lfcShrink(dds, coef="condition_N2_vs_CX11262", type="apeglm")
plotMA(exp_res, ylim=c(-8, 8))
```

#### Save Significant Results
```r
exp_res <- results(dds, contrast = c("condition", "N2", "CX11262"), alpha = 0.05, lfcThreshold = 0.00)
exp_res_sig <- subset(exp_res, padj < 0.05)
write.csv(exp_res_sig, "result.csv")
```

### 6. Annotation with biomaRt

#### Install Required Packages
```r
BiocManager::install("org.Ce.eg.db", force = TRUE)
BiocManager::install("AnnotationDbi", force = TRUE)
```

#### Load Libraries and Annotate Results
```r
library(org.Ce.eg.db)
library(AnnotationDbi)

deg_results <- read.csv(file.choose())
gene_names <- AnnotationDbi::mapIds(org.Ce.eg.db, keys = deg_results$wormbase_id, column = "SYMBOL", keytype = "WORMBASE", multiVals = "list", use.names = TRUE)
deg_results$gene_name <- gene_names[as.character(deg_results$wormbase_id)]
write.csv(deg_results, file = "updated_deg_results.csv", row.names = FALSE)
```

### 7. Visualization with Volcano Plot

#### Load Final Results and Visualize
```r
deg_results <- read.csv("updated_deg_results.csv")
deg_results$regulation <- ifelse(deg_results$padj < 0.05 & deg_results$log2FoldChange > 0, "Upregulated", ifelse(deg_results$padj < 0.05 & deg_results$log2FoldChange < 0, "Downregulated", "Not Significant"))

highlighted_genes <- c("srap-1", "C05D12.3", "Y9D1A.1", "fbxa-61", "fbxa-59")
highlight_genes <- deg_results %>% filter(gene_name %in% highlighted_genes)

volcano_plot <- ggplot(deg_results, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = regulation), alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "blue", "Downregulated" = "red", "Not Significant" = "gray")) +
  geom_text_repel(data = highlight_genes, aes(label = gene_name), size = 3, color = "black", fontface = "bold", box.padding = 0.3, point.padding = 0.2, bg.color = "white", bg.r = 0.15) +
  labs(title = "N2 vs CX11262", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  theme_minimal() +
  annotate("text", x = max(deg_results$log2FoldChange), y = max(-log10(deg_results$padj)) * 0.9, label = paste("Upregulated:", sum(deg_results$regulation == "Upregulated")), hjust = 1, color = "blue", size = 5) +
  annotate("text", x = max(deg_results$log2FoldChange), y = max(-log10(deg_results$padj)) * 0.85, label = paste("Downregulated:", sum(deg_results$regulation == "Downregulated")), hjust = 1, color = "red", size = 5)

print(volcano_plot)
```

### 8. Pathway Analysis

#### Load Required Libraries
```r
library(DESeq2)
library(clusterProfiler)
library(org.Ce.eg.db)
library(AnnotationDbi)
```

#### Perform Pathway Analysis
```r
deseq2_results <- read.csv(file.choose())
sig_genes <- deseq2_results[deseq2_results$padj < 0.05, ]
gene_ids <- sig_genes$wormbase_id
entrez_ids <- select(org.Ce.eg.db, keys = gene_ids, keytype = "WORMBASE", columns = "ENTREZID")
reactome_enrich <- enrichPathway(entrez_ids$ENTREZID, organism = "celegans")
dotplot(reactome_enrich)
```

This completes the RNA-Seq analysis pipeline from raw FASTQ files to pathway enrichment analysis.