# RNA-seq Analysis with STAR, Salmon, and DESeq2

This repository contains a full RNA-seq analysis pipeline using both alignment-based (STAR) and alignment-free (Salmon) methods, followed by downstream differential expression analysis using DESeq2.

---

## 📁 Repository Structure

```
RNAseq_analysis/
├── rawdata/                   # SRA or FASTQ input files
├── fastq/                    # FASTQ files extracted by fasterq-dump
├── genome/                   # Genome FASTA, GTF, STAR/salmon indices
├── star_index/              # STAR genome index
├── salmon_index/            # Salmon transcriptome index
├── salmon_quant/            # Quantification outputs from Salmon
├── tx2gene.tsv              # Transcript-to-gene mapping file
├── deseq2_analysis.Rmd      # RMarkdown pipeline using DESeq2
├── DESeq2_results_with_symbols.csv
├── GO_all_DEGs.csv
├── KEGG_all_DEGs.csv
└── figures/                 # Volcano plot, enrichment plots, heatmaps
```

---

## 🔧 Preprocessing: SRA Download & FASTQ Generation

```bash
# Connect to HPC cluster
ssh xinwang@xanadu-submit-ext.cam.uchc.edu

# Load SRA Toolkit
module load sratoolkit

# Download datasets
prefetch SRR17518171~SRR17518178
vdb-validate SRR17518171~SRR17518178

# Convert to FASTQ
fasterq-dump --split-files --threads 4 -O fastq/ rawdata/SRRxxxxxxx
```

---

## 🧬 Genome Preparation

```bash
# Download genome and annotation
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz

gunzip *.gz
```

### ➤ STAR Index
```bash
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir star_index \
     --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbGTFfile Homo_sapiens.GRCh38.113.gtf \
     --sjdbOverhang 100
```

### ➤ Salmon Index
```bash
salmon index \
  -t Homo_sapiens.GRCh38.cdna.all.fa \
  -i salmon_index \
  -k 31
```

---

## 🧪 Alignment & Quantification

### ➤ STAR Alignment
```bash
STAR --genomeDir star_index \
     --readFilesIn sample_1.fastq.gz sample_2.fastq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix star_output/sample_ \
     --outSAMtype BAM SortedByCoordinate
```

### ➤ Salmon Quantification
```bash
salmon quant \
  -i salmon_index \
  -l A \
  -1 sample_1.fastq.gz \
  -2 sample_2.fastq.gz \
  -o salmon_quant/sample_quant \
  -p 4 --validateMappings
```

---

## 📊 Differential Expression with DESeq2 (RMarkdown)

See `deseq2_analysis.Rmd` for a fully annotated pipeline, including:

- Sample import with `tximport`
- Differential gene expression via `DESeq2`
- Volcano plot and PCA
- GO & KEGG enrichment using `clusterProfiler`
- Heatmap of top DEGs

---

## 🧾 Transcript-to-Gene Mapping
Generated from GTF with:
```bash
awk '$3 == "transcript"' Homo_sapiens.GRCh38.113.gtf | \
awk '{ match($0, /transcript_id "([^"]+)"/, a); match($0, /gene_id "([^"]+)"/, b); if (a[1] && b[1]) print a[1] "\t" b[1]; }' > tx2gene.tsv
```

---

## 📬 Contact
For questions, contact: **xinwang@uchc.edu**
