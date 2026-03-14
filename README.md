# Single-Cell RNA-seq Analysis: Manual vs Azimuth Cell Annotation

## Introduction

Single-cell RNA sequencing (scRNA-seq) enables the characterization of cellular heterogeneity by profiling gene expression at the individual cell level. Accurate identification of cell types is a critical step in scRNA-seq analysis.

In this project, I performed a complete scRNA-seq analysis workflow using the **Seurat framework** and compared **manual cell type annotation** with **reference-based automated annotation using Azimuth**. In addition, I explored biological pathways enriched in specific cell clusters using **Gene Set Enrichment Analysis (GSEA)**.

This project demonstrates a **fully reproducible single-cell RNA-seq analysis pipeline**, including preprocessing, clustering, annotation, and functional interpretation.

---

# Methods

## Data Processing

The raw dataset was processed using the Seurat pipeline. The following steps were performed:

* Data loading from **10x Genomics format**
* Quality control filtering
* Normalization
* Identification of highly variable genes
* Principal Component Analysis (PCA)

---

## Clustering

Cells were clustered using a **graph-based clustering method**, followed by **UMAP dimensionality reduction** for visualization of cellular structure in low-dimensional space.

---

## Marker Gene Identification

Cluster-specific marker genes were identified using **differential expression analysis**, enabling characterization of gene signatures that distinguish each cluster.

---

## Cell Type Annotation

Two approaches were used to identify cell types:

### Manual Annotation

Cluster marker genes were examined and compared with known **PBMC marker genes** to infer cell identities.

### Automated Annotation

Cell identities were predicted using **Azimuth reference mapping**, which maps query cells to a curated reference atlas to assign cell type labels.

---

## Pathway Analysis

To investigate biological functions associated with specific clusters, **Gene Set Enrichment Analysis (GSEA)** was performed on cluster-specific differential expression results.

---

# Outputs

Figures generated during analysis are saved in the following directory:

`figures/`

These include:

* Quality control violin plots
* Variable feature plots
* PCA visualization
* UMAP clustering plots
* Marker gene dot plots
* Azimuth annotation comparison
* GSEA enrichment visualization

Tables and exported results are saved in:

`Results/`

---

# Results

## Comparison Between Manual and Azimuth Annotations

![alt text](<figures/Comparison between Manual and Azimuth Annotations.png>)

---

## GSEA Analysis

![alt text](<figures/GSEA analysis plot cluster_0 vs others.png>)



---

# Acknowledgment

This analysis and visualization were performed using **R**. Knowledge and skills applied in this project were developed through **DataCamp courses**, hands-on computational biology practice, and publicly available **10x Genomics datasets**.

