# Regenerate figures from research paper
- **Title: Single-cell multi-omic analysis of the vestibular schwannoma ecosystem uncovers a nerve injury-like state** 
- Journal: Nature Communications, 2023
- https://doi.org/10.1038/s41467-023-42762-w
- https://www.nature.com/articles/s41467-023-42762-w

# Data availability
All scRNA-seq, scATAC-seq, and new bulk RNA-seq data is available through the Gene Expression Omnibus with GEO accession “GSE216784”. All WES data is available through the database of Genotypes and Phenotypes (dbGaP) with accession “phs003318.v1.p1”.

**Source data to generate figures: https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-42762-w/MediaObjects/41467_2023_42762_MOESM6_ESM.xlsx**

# About this paper
## Objective
The paper aims to deeply characterize the cellular and molecular ecosystem of vestibular schwannoma (VS) — a benign tumor of the Schwann cells that form the myelin sheath around nerves — by using single-cell RNA-seq (scRNA-seq) and single-nucleus ATAC-seq (snATAC-seq).

The main goal was to understand the heterogeneity within tumor and surrounding microenvironment cells, and to reveal how VS resembles a nerve injury-like state, which might explain how it interacts with neurons and immune cells.

## Key Methods
- Collected VS tumor tissues from patients and performed:
- scRNA-seq to profile gene expression at single-cell resolution.
- snATAC-seq to profile chromatin accessibility and infer regulatory states.
- Integrated these data to identify cell types, subtypes, and regulatory programs.
- Compared VS cells to normal peripheral nerve and injury conditions.

## Main Findings
1. Schwannoma Cells Mimic Nerve Injury:
- Tumor Schwann cells show gene expression and chromatin states similar to those seen in injured peripheral nerves, including upregulation of repair-associated genes.
- They adopt a “repair Schwann cell-like” phenotype.
2. Immune Microenvironment:
- The tumor microenvironment has an enriched population of immune cells, including macrophages with an activated phenotype, which likely interact with tumor Schwann cells.
- Signals such as cytokines and chemokines suggest cross-talk between tumor and immune cells.
3. Regulatory Programs:
- They identified key transcription factors (e.g., c-Jun) and regulatory elements driving this injury-like state.
- The chromatin accessibility profiles confirm reactivation of injury repair programs.
4. Therapeutic Implication: The discovery of an injury-like program suggests possible targets to modulate the repair-like state or immune interactions in VS, potentially opening new avenues for non-surgical therapies.

## Takeaway:
Vestibular schwannomas reactivate an injury-response program, reprogramming Schwann cells and shaping a supportive immune microenvironment. The study highlights the power of single-cell multi-omics to uncover these complex cellular states and interactions.

# 📊 Figure 1 Overview (excluding 1a and 1b)
1. Figure 1c: scRNA-seq UMAP
2. Figure 1d: scATAC-seq UMAP
3. Figure 1e: Dot plot of marker genes
4. Figure 1f: Heatmap

# Figure 2 Overview


## Fig 1c and d use 
- UMAP_1 and 2 data
- Label of different cell types/clusters
## Figure 1e shows:
- Y-axis: major cell types/clusters (grouped together by Group)
- X-axis: marker genes
- Dot size: % of cells expressing the gene
- Dot color: scaled average expression
## Figure 1f shows (note: ComplexHeatmap wants a matrix):
- Rows = Mouse Peripheral Nerve Cell Type
- Columns = Clusters
- Cell values = Mean Module Score
