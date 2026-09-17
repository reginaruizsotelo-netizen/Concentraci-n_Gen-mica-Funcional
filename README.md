<h1 align="center">Functional_Genomic_Concentration</h1>

<p align="center">
  Comparative transcriptomic analysis of <i>Coffea arabica</i> cultivars Catuaí and CR95 in response to <i>Xylella fastidiosa</i> infection.
</p>

<p align="center">
  <a href="#about">About</a> ·
  <a href="#tech-stack">Tech Stack</a> ·
  <a href="#methodology">Methodology</a> ·
  <a href="#files-in-this-repository">Files</a> ·
  <a href="#key-results">Results</a>
</p>

---

## About

This repository contains the analysis pipeline used to study the differential transcriptional response of two coffee cultivars — **Catuaí** (susceptible) and **CR95** (resistant) — to *Xylella fastidiosa* infection. The workflow covers RNA-seq preprocessing, differential expression analysis, functional enrichment, and regulatory motif discovery in promoter regions of candidate defense-related genes.

**Dataset:** 15 RNA-seq samples
- Catuaí: n = 8 (infected = 5, control = 3)
- CR95: n = 7 (infected = 3, control = 4)

## Tech Stack

The analysis was run using:

- R 4.6.1
- DESeq2 v1.52.0
- edgeR v4.10.5
- limma v3.68.5
- clusterProfiler v4.20.0
- rtracklayer v1.72.0
- Rsamtools v2.28.0
- Biostrings v2.80.2
- memes v1.20.0 (MEME v5.5.9)
- HISAT2 v2.2.3
- featureCounts (Subread) v2.1.1
- Trimmomatic v0.39
- FastQC v0.12.1

## Methodology

### 1. Quality control
- **FastQC** v0.12.1 — pre- and post-trimming QC (Figures S1, S2)

### 2. Read trimming
- **Trimmomatic** v0.39, paired-end mode (CaC10 processed single-end, missing forward read)
  - `ILLUMINACLIP`: TruSeq3 adapter reference
  - `HEADCROP`: 15 bp
  - `TRAILING`: Q30
  - `SLIDINGWINDOW`: 4:30

### 3. Alignment
- **HISAT2** v2.2.3 (splice-aware)
  - Reference: *Coffea arabica* genome `GCF_036785885.1`

### 4. Read quantification
- **featureCounts** v2.1.1
  - Input: HISAT2 BAM files + GTF/GFF annotation

### 5. Exploratory PCA
📄 `Coffee_PCA.R`
- Low-expression filter: counts-per-million (CPM)
- Variance-stabilizing transformation (VST) — **DESeq2** v1.52.0
- Guided design: model treatment separately per cultivar

### 6. Differential expression analysis (DEGs)
📄 `DEGs_Heatmaps_Networks_Enrichment.R`
📄 `Evidence_1_Report_Supplementary_material.pdf` — contains gene co-expression networks per cultivar (Pearson correlation, threshold > 0.8), showing distinct topological patterns between Catuaí (downregulation-dominated) and CR95 (upregulation-dominated)
- **DESeq2** v1.52.0, run independently per cultivar
- Low-expression filter: `filterByExpr` — **edgeR** v4.10.5
- Design: `~ treatment`
- Shrunken log2FC: **ashr** method
- DEG criteria (both required): `padj` (FDR) < 0.05 and `|log2FC|` > 1

### 7. GO enrichment — competitive test
📄 `DEGs_Heatmaps_Networks_Enrichment.R`
- **CAMERA** function — **limma** v3.68.5
- FDR < 0.05
- Minimum 5 genes per term

### 8. GO enrichment — hypergeometric test
📄 `DEGs_Heatmaps_Networks_Enrichment.R`
📄 `CR95_padj0.05_log2FC1.xlsx`, `Catuai_padj0.05_log2FC1.xlsx` — filtered DEG lists used as input
- **clusterProfiler** v4.20.0, `enricher` function
- Run separately for upregulated and downregulated genes per cultivar
- Thresholds: `|log2FC|` > 1, `padj` < 0.05

### 9. Regulatory motif analysis
📄 `MEME.R`
- 8 candidate genes (ethylene-pathway related, CR95 upregulated)
- Coordinates extracted with **rtracklayer** v1.72.0 (`GCF_036785885.1` annotation)
- Promoter region: **-1000 bp / +200 bp** relative to TSS (strand-aware)
- Sequence extraction: **Rsamtools** v2.28.0 + **Biostrings** v2.80.2
- Motif discovery: **MEME** v5.5.9 (via **memes** v1.20.0), width 6–15 bp, max 3 motifs
- Motif comparison: **TomTom** vs. JASPAR Plants database

## Files in this repository

| File | Used in |
|---|---|
| `Coffee_PCA.R` | Step 5 — Exploratory PCA |
| `DEGs_Heatmaps_Networks_Enrichment.R` | Steps 6, 7, 8 — DEG analysis, heatmaps, GO/hypergeometric enrichment |
| `CR95_padj0.05_log2FC1.xlsx` | Step 8 — CR95 DEG input |
| `Catuai_padj0.05_log2FC1.xlsx` | Step 8 — Catuaí DEG input |
| `MEME.R` | Step 9 — Promoter motif discovery |
| `Evidence_1_Report_Supplementary_material.pdf` | Step 6 — co-expression networks per cultivar |

## Key results

- Catuaí: 741 DEGs (padj<0.05) → 38 after |log2FC|>1 filter (7 up / 31 down)
- CR95: 363 DEGs (padj<0.05) → 225 after |log2FC|>1 filter (216 up / 9 down)
- CR95 mounts a defense response centered on ethylene and fatty acid signaling; Catuaí shows a broader but low-magnitude, less specific transcriptional response
- Motif analysis: top hits matched MYB-family TFs (MYB23, TCX3, ATMYB31) — no direct match to classical ethylene-pathway elements
