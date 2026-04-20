# WAT-regulation-of-secondary-growth-and-auxin-homeostasis-in-hybrid-poplar Pipeline (Bulk RNA-seq + Deconvolution + Functional Analysis)
Overview
This repository contains the full computational pipeline used to analyze gene expression and functional enrichment in Populus tremula × alba under drought-related conditions.
The workflow integrates:
Bulk RNA-seq differential expression analysis (edgeR)
Orthology mapping between P. tremula × alba and P. trichocarpa
Cell-type deconvolution using BayesPrism
Functional enrichment analysis (topGO, GSEA, KEGG)
The pipeline is modular and organized into independent scripts that can be executed sequentially.
📂 Pipeline Structure
1. Bulk Differential Expression (edgeR)
Scripts:
Bulk_edgeR_C1_additive_EN.R
Bulk_edgeR_C9_additive_EN.R
Purpose:
Perform differential expression analysis separately for:
WT vs C1
WT vs C9
Key features:
TMM normalization
Filtering using filterByExpr
Additive model controlling for ZT (time effect):
~ ZT + treatment
Robust dispersion estimation (estimateDisp, glmQLFit)
Outputs:
Full DE results:
Separado_C1_Aditivo_Model.csv
Separado_C9_Aditivo_Model.csv
Filtered DEGs (FDR < 0.05):
Separado_C1_FDR_005_Aditivo_Model.csv
Separado_C9_FDR_005_Aditivo_Model.csv
2. Orthology Mapping + Mixture File Generation
Script:
1_Mixture_File_Generation_EN.R
Purpose:
Convert gene IDs:
P. tremula × alba → P. trichocarpa
Generate the mixture file required for deconvolution
Orthology strategy:
1:1 → direct mapping
1:N → keep all orthologs
N:1 → prioritize DEGs + highest WT expression
N:N → select first alphabetically
Outputs:
Orthology dictionary:
PtremxalbaHAP2_to_Ptricho.csv
Mixture file:
Tricho_WT_Mutant_BULK_Counts.tsv
3. Deconvolution (BayesPrism)
Script:
2_BayesnPrism_EN.R
Purpose:
Decompose bulk RNA-seq into cell-type-specific expression profiles
Inputs:
Bulk mixture file
scRNA-seq reference (Seurat object)
Marker genes (Seurat or BayesPrism)
Key steps:
Marker selection (Seurat FindAllMarkers)
Bayesian deconvolution
Extraction of:
Cell fractions
Tissue-specific pseudocount matrices
Outputs:
Cell-type fractions (SVG plot)
Expression matrices per tissue:
Z_Matrix_Pseudocounts_Deconvoluted_Tissues/*.csv
4. ORA Enrichment (Cell-Type Specific)
Script:
ORA_topGO_EN.R
Purpose:
Perform GO enrichment (topGO) on:
Upregulated genes
Downregulated genes
per cell type (deconvolved)
Key features:
Uses Arabidopsis GO annotations
Algorithm: weight01 + Fisher
BH correction
Outputs:
Enriched GO terms:
.csv + .xlsx
GO + genes:
_With_Genes_Enriched_topGo.csv
GO + genes + functional annotation:
_With_Genes_Enriched_And_Arab_Functional_Info_topGo.csv
Dotplots (top 50 terms)
5. ORA Enrichment (Bulk Shared DEGs)
Script:
ORA_topGO_bulk_shared_EN.R
Purpose:
Enrichment of shared DEGs between C1 and C9:
UP (both mutants)
DOWN (both mutants)
Differences vs previous ORA:
No tissue resolution
Single gene list (intersection)
Outputs:
Dotplot (top 20 GO terms)
CSV:
GO + genes + Arabidopsis annotation
6. GSEA (KEGG Pathways)
Script:
GSEA_KEGG_one_tailed_EN.R
Purpose:
Perform GSEA (clusterProfiler) using a custom ranking:
ranking=−log 
10
​	
 (p 
Fisher
​	
 )⋅sign( 
log2FC
​	
 )
Key features:
Fisher method combines C1 + C9 p-values
One-tailed enrichment:
"pos" → upregulated pathways
"neg" → downregulated pathways
KEGG mapping via:
KEGGREST
bitr_kegg
Outputs:
Enriched KEGG pathways:
*_All_Enriched_KEGG_Pathways.csv
Visualizations:
Dotplots
Network plots (cnetplot)
🔁 Workflow Summary
Bulk RNA-seq
     ↓
edgeR (C1 & C9 separately)
     ↓
Shared / tissue-specific DEGs
     ↓
Orthology mapping (Populus → Ptrichocarpa)
     ↓
BayesPrism deconvolution
     ↓
Cell-type expression matrices
     ↓
Functional analysis:
   ├── ORA (topGO)
   └── GSEA (KEGG)
⚙️ Requirements
R (≥ 4.3 recommended)
Packages:
edgeR
limma
tidyverse
Seurat
BayesPrism
topGO
clusterProfiler
KEGGREST
svglite
openxlsx
🧠 Key Design Decisions
ZT included in model → avoids confounding circadian effects
Separate C1 and C9 models → accounts for different variance structures
Orthology mapping to P. trichocarpa → required for compatibility with scRNA reference
Use of pseudocounts → enables tissue-specific DE analysis
Combined Fisher p-value (GSEA) → integrates both mutants
📌 Notes
Some scripts are modular but not fully automated (manual parameter control required)
KEGG queries depend on external API (may require retry if rate-limited)
Deconvolution results depend strongly on marker selection quality
🚀 Suggested Improvements
Wrap pipeline into a Snakemake / Nextflow workflow
Add:
logging
config files
Standardize:
paths
naming conventions
Add reproducibility layer (e.g. renv)
