# WAT-regulation-of-secondary-growth-and-auxin-homeostasis-in-hybrid-poplar Pipeline (Bulk RNA-seq + Deconvolution + Functional Analysis)

---

## Overview

This repository contains a complete computational pipeline to analyze gene expression and functional enrichment in **Populus tremula × alba**.

The workflow integrates:

- Bulk RNA-seq differential expression (**edgeR**)
- Orthology mapping (**P. tremula × alba → P. trichocarpa**)
- Cell-type deconvolution (**BayesPrism**)
- Functional enrichment (**topGO, GSEA, KEGG**)

---

## Pipeline Structure

### 1. Bulk Differential Expression (edgeR)

**Scripts:**
- `Bulk_edgeR_C1_additive_EN.R`
- `Bulk_edgeR_C9_additive_EN.R`

**Description:**
Differential expression analysis for:
- WT vs C1
- WT vs C9

**Key steps:**
- Filtering (`filterByExpr`)
- TMM normalization
- Additive model: ~ ZT + treatment

- - Robust dispersion (`glmQLFit`)

**Outputs:**
- Full results:
- `Separado_C1_Aditivo_Model.csv`
- `Separado_C9_Aditivo_Model.csv`
- Filtered (FDR < 0.05):
- `Separado_C1_FDR_005_Aditivo_Model.csv`
- `Separado_C9_FDR_005_Aditivo_Model.csv`

---

### 2. Orthology Mapping + Mixture File

**Script:**
- `1_Mixture_File_Generation_EN.R`

**Description:**
- Maps genes: P. tremula × alba HAP2 → P. trichocarpa
- - Generates mixture file for deconvolution

**Strategy:**
- 1:1 → keep
- 1:N → keep all
- N:1 → prioritize DEGs + highest WT expression
- N:N → select first alphabetically

**Outputs:**
- `PtremxalbaHAP2_to_Ptricho.csv`
- `Tricho_WT_Mutant_BULK_Counts.tsv`

---

### 3. Deconvolution (BayesPrism)

**Script:**
- `2_BayesnPrism_EN.R`

**Description:**
Deconvolves bulk RNA-seq into **cell-type-specific expression**

**Inputs:**
- Bulk mixture file
- scRNA-seq reference (Seurat)
- Marker genes

**Outputs:**
- Cell fractions (SVG)
- Tissue-specific matrices: Z_Matrix_Pseudocounts_Deconvoluted_Tissues/*.csv
- 
---

### 4. ORA Enrichment (Cell-Type Specific)

**Script:**
- `ORA_topGO_EN.R`

**Description:**
GO enrichment per tissue using **topGO**

**Features:**
- Algorithm: `weight01 + Fisher`
- BH correction

**Outputs:**
- `.csv` and `.xlsx` enriched GO terms
- GO + genes
- GO + genes + annotations
- Dotplots (top 50)

---

### 5. ORA Enrichment (Bulk Shared DEGs)

**Script:**
- `ORA_topGO_bulk_shared_EN.R`

**Description:**
Enrichment of DEGs shared between C1 and C9

**Outputs:**
- Dotplot (top 20 GO terms)
- CSV with:
- GO terms
- genes
- Arabidopsis annotation

---

### 6. GSEA (KEGG Pathways)

**Script:**
- `GSEA_KEGG_one_tailed_EN.R`

**Description:**
GSEA using custom ranking: ranking = -log10(p_fisher) * sign(log2FC_avg)


**Features:**
- Fisher p-value (C1 + C9)
- One-tailed:
  - `pos` → upregulated
  - `neg` → downregulated

**Outputs:**
- Enriched pathways:
  - `*_All_Enriched_KEGG_Pathways.csv`
- Dotplots
- Network plots (cnetplot)

---

## 🔁 Workflow
Bulk RNA-seq
↓
edgeR (C1 / C9)
↓
Functional analysis:
├── ORA (topGO; Shared DEGs)
└── GSEA (KEGG)

↓
Orthology mapping
↓
BayesPrism
↓
edgeR (C1 / C9 + One vs Rest)
↓
Shared DEGs C1 + C9 ∩ tissue DEGs
↓
ORA (topGO)



