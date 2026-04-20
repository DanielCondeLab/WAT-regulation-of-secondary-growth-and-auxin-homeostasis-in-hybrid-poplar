# WAT-regulation-of-secondary-growth-and-auxin-homeostasis-in-hybrid-poplar Bioinformatic Pipeline 

---

## Overview

This repository contains the complete computational pipeline to analyze gene expression and functional enrichment in **Populus tremula × alba HAP2** of (poner paper y doi)

The workflow integrates:

- Bulk RNA-seq and deconvolved differential expression (**edgeR**)
- Orthology mapping (**P. tremula × alba HAP2 v5.1 → P. trichocarpa v4.1**)
- Cell-type deconvolution (**BayesPrism**)
- Functional enrichment (**ORA GO (topGO); GSEA KEGG (clusterProfiler)**)

---
## Requirements

- R (v 4.4.2)
- Conda enviroment files are provided for each script of the pipeline (.yaml)

**Note**: BayesPrism R package must be installed via R

---

## Pipeline Structure

### 1. Bulk Differential Expression

**Scripts (DEGs_Bulk):**
- `C1_ZT_Aditive_Bulk_edgeR.R`
- `C9_ZT_Aditive_Bulk_edgeR.R`

**Description:**
This scripts perform differential expression analysis for:
- C1 vs WT
- C9 vs WT

**Key steps:**
- Filtering low expressed genes
- TMM normalization
- Additive model: ~ ZT + Genotype
- Fitted to Genewise Negative Binomial Generalized Linear Models

**Inputs:**
- Bulk raw count matrix for 30 samples:
  - `Counts_30samples.csv`

**Outputs:**
- Full results for each genotype:
  - `Separated_*_Aditive_Model.csv`
- Filtered DEGs for each genotype (FDR < 0.05):
  - `Separated_*_FDR_005_Aditive_Model.csv`

---

### 2. ORA Enrichment of DEGs from Bulk 

**Scripts (ORA_Bulk):**
- `ORA_topGO_Bulk.R`

**Description:**
ORA GO enrichment of shared DEGs from C1 vs WT and C9 vs WT 

for up and downregulated respectively

**Features:**
- Biological Process ontology
- Algorithm: `weight01 + Fisher`
- BH p-value correction

**Inputs:**
- Filtered DEGs from Bulk (FDR < 0.05):
  - `Separated_C1_FDR_005_Aditive_Model.csv`
  - `Separated_C9_FDR_005_Aditive_Model.csv`
- Poplar–Arabidopsis GO annotation database:
  - `2026_Poplar_Arab_GOs_TAIR.csv`
- TAIR functional descriptions:
  - `TAIR10_functional_descriptions.csv`

**Outputs:**
-  Enriched GO terms with contributing genes and Arabidopsis functional information
    - `*_All_Significant_GO_Enrichment_With_Arab_Functional_Info.csv`
- Dotplots (top 20 enriched terms) for each direction
  - `*_Top_20_ORA_topGO_Bulk.svg`

---

### 3. GSEA of Bulk (KEGG Pathways)

**Scripts (GSEA_Bulk):**
- `GSEA_KEGG_Bulk.R`

**Description:**
GSEA of KEGGs using custom ranking: -log10(p_fisher) * sign(log2FC_avg)

**Features:**
- Fisher method to combine p-values and average of logFC (C1 + C9)
- One-tailed GSEA:
  - `pos` → upregulated
  - `neg` → downregulated

**Outputs:**
- Enriched pathways:
  - `*_All_Enriched_KEGG_Pathways.csv`
- Dotplots
- Network plots (cnetplot)

### 4. Orthology Mapping + Deconvolution Input File

**Scripts (Orthologues_Poplars):**
- `1_Orthologue_Dictionary_Generation.R`

**Description:**
- Maps genes: P. tremula × alba HAP2 → P. trichocarpa
- - Generates mixture file for deconvolution

**Strategy:**
- 1:1 → keep all
- 1:N → keep all
- N:1 → prioritize DEGs + highest WT expression
- N:N → select first alphabetically

**Inputs:**
- Bulk raw count matrix:
  - `Counts_30samples.csv`
- Differential expression results for C1 and C9:
  - `Separated_*_FDR_005_Aditive_Model.csv`
- Phytozome Biomart poplar ortholog file:
  - `PtaHAP2vsTrico.txt`

**Outputs:**
- Poplar orthologue dictionary:
  - `PtremxalbaHAP2_to_Ptricho.csv`
- RNA-seq input for Deconvolution:
  - `Tricho_WT_Mutant_BULK_Counts.tsv`

---

### 5. Deconvolution of Bulk Celltypes 

**Scripts (Deconvolution_Analysis):**
- `BayesnPrism.R`

**Description:**
Deconvolves bulk RNA-seq into **cell-type-specific expression**

**Inputs:**
- P.trichocarpa scRNA-seq Seurat object (DOI: 10.1186/s13059-025-03728-x):
  - `Integrated_dataSnRNAseqSTEMFinalclustering.rds`
- Seurat marker file (optional, if already computed for faster):
  - `Filtered_Ptricho_FindAllMarkers_Seurat.csv`
- RNA-Seq input for Deconvolution:
  - `Tricho_WT_Mutant_BULK_Counts.tsv`

**Outputs:**
- Mean cell fractions plots for each genotype (`.svg`)
- Tissue-specific expression matrices (`.csv`)

---

### 6. Deconvolved Celltype Differential Expression 

**Scripts (Deconvolution_Analysis):**
- `BayesPrism_DEGs.R`

**Description:**
This scripts perform differential expression analysis for the deconvolved celltypes:
- C1 vs WT
- C9 vs WT

**Key steps:**
- Filtering low expressed genes
- TMM normalization
- Additive model: ~ ZT + Genotype
- Fitted to Genewise Negative Binomial Generalized Linear Models

**Inputs:**
- Tissue-specific expression matrices (`.csv`):
- Poplar orthologue dictionary:
  - `PtremxalbaHAP2_to_Ptricho.csv`
    
**Outputs:**
- Celltype Full results for each genotype:
  - `*_Separated_*_Aditive_Model.csv`
- Celltype filtered DEGs for each genotype (FDR < 0.05) :
  - `*_Separated_*_FDR_005_Aditive_Model.csv`

---

### 7. Deconvolved Celltype Marker Identification 

**Scripts (Deconvolution_Analysis):**
- `BayesPrism_One_vs_Rest.R`

**Description:**
This scripts identifies highly celltype-specific genes for the deconvolved celltypes

using an One vs Rest aproach:
- Gene expression in X vs mean gene expression in the rest of celltypes

**Key steps:**
- Filtering low expressed genes
- TMM normalization
- Moodel: ~ 0 + celltype
- Fitted to Genewise Negative Binomial Generalized Linear Models

**Inputs:**
- Tissue-specific expression matrices (`.csv`)
  
**Outputs:**
- Marker genes for all the celltypes:
  - `Resumed_General_Markers_FDR_005.csv`

---

### 8. ORA Enrichment of Cell-Type Genes

**Script:**
- `ORA_topGO_Deconvolution.R`

**Description:**
GO enrichment per tissue using **topGO**

**Features:**
- Algorithm: `weight01 + Fisher`
- BH correction

**Inputs:**
- Celltype filtered DEGs for each genotype (FDR < 0.05) :
  - `*_Separated_*_FDR_005_Aditive_Model.csv`
- Marker genes fo all the celltypes:
  - `Resumed_General_Markers_FDR_005.csv`

**Outputs:**
- Enriched GO terms (`.csv` and `.xlsx`)
- Enriched GO terms + contributing genes (`.csv`)
- Enriched GO terms + contributing genes + functional annotations (`.csv`)
- Dotplots (top 50) of most enriched processes (`.svg`)

---


## Workflow

RNA-seq<br>
↓<br>
DEGs analysis edgeR (C1 / C9)<br>
↓<br>
Functional analysis:<br>
├── ORA (topGO; Shared Up or Dwn DEGs)<br>
└── GSEA (KEGG)<br>
<br>

Poplars orthology mapping<br>
↓<br>
Deconvolution of RNA-seq with BayesPrism<br>
↓<br>
DEGs analysis and cell type specific identification (edgeR; C1 / C9 + One vs Rest)<br>
↓<br>
Shared Up or Dwn DEGs C1 + C9 ∩ tissue DEGs<br>
↓<br>
ORA (topGO)

