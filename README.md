# WAT-regulation-of-secondary-growth-and-auxin-homeostasis-in-hybrid-poplar Bioinformatic Pipeline 

---

## Overview

This repository contains the complete computational pipeline to analyze gene expression and functional enrichment in **Populus tremula × alba HAP2** 

The workflow integrates:

- Bulk RNA-seq and deconvolved celltype differential expression (**edgeR**)
- Orthology mapping (**P. tremula × alba HAP2 v5.1 → P. trichocarpa v4.1** | **P. tremula × alba HAP2 v5.1 → P. Alba (KEGG)**)
- Cell-type deconvolution (**BayesPrism**)
- Functional enrichment (**ORA GO (topGO); GSEA KEGG (clusterProfiler)**)

poner paper y doi
---

## Pipeline Scheme (Incluir esquema)

---
## Software and Databases Requirements
- Databases
  - Populus tremula x alba HAP2 from Phytozome (https://phytozome-next.jgi.doe.gov/info/  PtremulaxPopulusalbaHAP2_v5_1)
  - Populus alba from KEGG (https://www.kegg.jp/kegg-bin/show_organism?org=palz)
  - Biomart Phytozome Tool (https://phytozome-next.jgi.doe.gov/biomart/martview/a5f9a612e8d5ed5ca96db2f9713cb466)
  - Arabidopsis thaliana TAIR10 from Phytozome (https://phytozome-next.jgi.doe.gov/info/Athaliana_TAIR10)
  - Arabidopsis thaliana TAIR10 functional info (https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_functional_descriptions)
  - Populus trichocarpa scRNA-seq reference (https://doi.org/10.1186/s13059-025-03728-x)
- Software
  - Conda enviroment files are provided for each script of the pipeline (.yaml)
  - OrthoFinder configuration file and 1:1 parser script
  - R (v 4.4.2) **Note**: BayesPrism R package must be installed via R
---

## Pipeline Worflow

### 1. Bulk Differential Expression (DE)

**Scripts:**
- `C1_ZT_Aditive_Bulk_edgeR.R`
- `C9_ZT_Aditive_Bulk_edgeR.R`

**Description:**

This scripts perform differential expression analysis for:
- C1 vs WT
- C9 vs WT

**Features:**
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

**Scripts:**
- `ORA_topGO_Bulk.R`

**Description:**

ORA GO enrichment of shared significant DEGs from C1 vs WT and C9 vs WT 
for up and downregulated respectively

**Features:**
- Biological Process ontology
- Algorithm: `weight01 + Fisher`
- BH p-value correction

**Inputs:**
- Filtered DEGs from Bulk (FDR < 0.05):
  - `Separated_*_FDR_005_Aditive_Model.csv`
- Poplar–Arabidopsis GO annotation database:

  Created by intersecting Phytozome P.tremula x alba HAP2 and A.thaliana 
- TAIR functional descriptions:
  - `TAIR10_functional_descriptions.csv`

**Outputs:**
-  Enriched GO terms with contributing genes and Arabidopsis functional information
    - `*_All_Significant_GO_Enrichment_With_Arab_Functional_Info.csv`
- Dotplots (top 20 enriched terms) for each direction
  - `*_Top_20_ORA_topGO_Bulk.svg`

---
### 3. Orthologue Mapping Between P.tremula x alba HAP2 - P.alba (KEGG)

**Description:**

This scripts identifies the orthologues between P.tremula x alba HAP2 from Phytozome and
P.alba from KEGG

**Scripts:**
- Activate orthofinder.yaml in the Terminal and copy config.json parameters
- Create and move to a folder (e.g. Poplar_Proteomes) and add P.tremula x alba and P.alba proteomes
- Execute: `orthofinder -M msa -f Poplar_Proteomes/`
- Enter results folder. Extract 1:1 pairs using `Poplars_parser.py`

**Features:**
- DIAMOND for protein homolgy searches
- FAMSA for multiple sequence alignment
- IQ-TREE 3 for phylogenetic inference 
  - ModelFinder automatic model selection
  - Ultrafast bootstrap with 1000 replicates  

**Inputs:**
- Poplars proteomes (principal transcripts)

**Outputs:**
- 1:1 Poplar Orthologues:
  - `one2one_pairs.tsv`

---

### 4. GSEA of Bulk (KEGG Pathways)

**Scripts:**
- `GSEA_KEGG_Bulk.R`

**Description:**

GSEA of DE results from bulk using KEGG database 

**Features:**
- Fisher method to combine p-values and average of logFC (C1 + C9)
- Custom ranking: -log10(p_fisher) * sign(log2FC_avg)
- One-tailed GSEA:
  - `pos` → upregulated
  - `neg` → downregulated

**Inputs:**
- Annotation info of P. tremula x alba HAP2 from Phytozome
- DE analysis results for each genotype:
  - `Separated_*_Aditive_Model.csv`

**Outputs:**
- Enriched pathways:
  - `*_All_Enriched_KEGG_Pathways.csv`
- Dotplots
- Network plots (cnetplot)

---
### 5. Visualization of KEGG routes with Pathview

**Scripts:**
- `GSEA_Pathview.R`

**Description:**

Visualization of enriched KEGG routes from GSEA with Pathview software

**Features:**
- Only significant genes are selected (FDR < 0.05)
- logFC direction is represented (magnitude is not considered)

**Inputs:**
- Annotation info of P. tremula x alba HAP2 from Phytozome
- DE analysis results for each genotype:
  - `Separated_*_Aditive_Model.csv`
- Significant genes of DE analysis for each genotype:
  - `Separated_*_FDR_005_Aditive_Model.csv`
- 1:1 Poplar Orthologues:
  - `one2one_pairs.tsv`

**Outputs:**
- Enriched map of interest (e.g. palz00400):
  - `Map_00400.csv`
  - `palz00400.png`

---
### 6. Orthology Mapping + Deconvolution Input File

**Scripts:**
- `Orthologue_Dictionary_Generation.R`

**Description:**

Maps genes between P. tremula × alba HAP2 → P. trichocarpa and Generates input file for deconvolution

**Features:**
- 1:1 → keep all
- 1:N → keep all
- N:1 → prioritize DEGs + highest WT expression
- N:N → select first alphabetically

**Inputs:**
- Bulk raw count matrix:
  - `Counts_30samples.csv`
- Significant genes of DE analysis for each genotype:
  - `Separated_*_FDR_005_Aditive_Model.csv`
- Phytozome Biomart P.tremula x alba - P.trichocarpa ortholog file:

**Outputs:**
- Poplar orthologue dictionary:
  - `PtremxalbaHAP2_to_Ptricho.csv`
- Input matrix for deconvolution:
  - `Tricho_WT_Mutant_BULK_Counts.tsv`

---

### 7. Deconvolution of Bulk Celltypes with BayesPrism 

**Scripts:**
- `BayesnPrism.R`

**Description:**

Deconvolves bulk RNA-seq data into **cell-type-specific expression**

**Features:**
- Although BayesPrism’s marker calculation is faster, it identifies far fewer genes than Seurat. For this reason, we recommend using Seurat’s FindAllMarkers.
- Sieve Elements (SE) celltype from P.trichocarpa scRNA-seq reference generates aberrant data.
  It is integrated into the Unknown (???) celltype. 


**Inputs:**
- P.trichocarpa scRNA-seq Seurat object (DOI: 10.1186/s13059-025-03728-x):
  - `Integrated_dataSnRNAseqSTEMFinalclustering.rds`
- Seurat marker file (optional, if already computed for faster performance):
  - `Filtered_Ptricho_FindAllMarkers_Seurat.csv`
- RNA-Seq input for Deconvolution:
  - `Tricho_WT_Mutant_BULK_Counts.tsv`

**Outputs:**
- Mean cell fractions plots for each genotype (`.svg`)
- Tissue-specific expression matrices (`.csv`)

---

### 8. Deconvolved Celltype Differential Expression Analysis 

**Scripts (Deconvolution_Analysis):**
- `BayesPrism_DEGs.R`

**Description:**

This scripts perform differential expression analysis for the deconvolved celltypes:
- C1 vs WT
- C9 vs WT

Parameters are the same as in the bulk DE analysis.

**Features:**
- Filtering low expressed genes
- TMM normalization
- Additive model: ~ ZT + Genotype
- Fitted to Genewise Negative Binomial Generalized Linear Models

**Inputs:**
- Celltype-specific expression matrices from BayesPrism(`.csv`):
- Poplar orthologue dictionary:
  - `PtremxalbaHAP2_to_Ptricho.csv`
    
**Outputs:**
- Celltype Full results for each genotype:
  - `*_Separated_*_Aditive_Model.csv`
- Celltype filtered DEGs for each genotype (FDR < 0.05) :
  - `*_Separated_*_FDR_005_Aditive_Model.csv`

---

### 9. Deconvolved Celltype Marker Identification 

**Scripts (Deconvolution_Analysis):**
- `BayesPrism_One_vs_Rest.R`

**Description:**

This scripts identifies highly celltype-specific genes for the deconvolved celltypes
using an One vs Rest aproach:
- Gene expression vs mean gene expression in the rest of celltypes for each gene

**Features:**
- Filtering low expressed genes
- TMM normalization
- Model: ~ 0 + celltype
- Fitted to Genewise Negative Binomial Generalized Linear Models

**Inputs:**
- Tissue-specific expression matrices (`.csv`)
  
**Outputs:**
- Marker genes for all the celltypes:
  - `Resumed_General_Markers_FDR_005.csv`

---

### 10. ORA Enrichment of Cell-Type Genes

**Script:**
- `ORA_topGO_Deconvolution.R`

**Description:**

ORA GO enrichment of shared significant DEGs from C1 vs WT and C9 vs WT per celltype 
for up and downregulated respectively

**Features:**
- Biological Process ontology
- Algorithm: `weight01 + Fisher`
- BH p-value correction

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
