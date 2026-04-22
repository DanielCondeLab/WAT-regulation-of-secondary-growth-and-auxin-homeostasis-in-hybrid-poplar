#########################################################
# CONDA ENVIRONMENT: Bayes.yaml
# DEGs from BayesnPrism Z-Matrices. zApplied to ALL celltypes (.csv)
# Same additive model (~zt + genotype) as bulk RNA-seq DEGS is used.
# Since C9 replicates have much more variability than C1, two parallel designs and fits are generated.
# Additionally, the mean pseudocounts of WT and Mutants (C1 and C9) are calculated for the Heatmaps plots.
# File routes MUST be changed in order to use this code.
# Input files:
# - Expression Matrices (pseudocounts) per celltype
# - P.tremula x alba - P.trichocarpa Dictionary (1:1). Obtained from Orthologue_Dictionary_Generation.R
# Output files:
# - DEGs (FDR < 0.05) of C1 for each celltype
# - DEGs (FDR < 0.05) of C9 for each celltype
# edgeR Manual: https://bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
#########################################################

# =============================================================================
# 0. IMPORT LIBRARIES
# =============================================================================

library(BayesPrism)
library(Seurat)
library(tidyverse)
library(edgeR)

# =======================================================================================================================
# 1. IMPORT EXPRESSION MATRIX DATA PER CELLTYPE (PSEUDOCOUNTS) AND THE P.TREMULA X ALBA - P.TRICHOCARPA DICTIONARY
# =======================================================================================================================

# Path with pseudocount matrices from each BayesPrism celltype
ruta_datos <- "/Users/danielconde/Desktop/WATs/BayesnPrism/Results/Z_Matrix_Pseudocounts_Deconvoluted_Celltypes/"

# P.Trem x Alba HAP2 - P.tricho Dictionary
PAlba_to_Potri <- read_csv(
  "/Users/danielconde/Desktop/WATs/BayesnPrism/Data/Mapeo_Phytozome_P.tremHAP2_P.tricho/PtremxalbaHAP2_to_Ptricho.csv"
)

dict_poplars <- setNames(
  PAlba_to_Potri$Alba,
  PAlba_to_Potri$Trichocarpa
)

archivos <- list.files(
  ruta_datos,
  pattern = "\\.csv$",
  full.names = TRUE
)

# =============================================================================
# 2. MAIN LOOP → APPLIED TO ALL PSEUDOCOUNTS
# =============================================================================

for (archivo in archivos) {
  
  mensaje <- paste0("\n========== Processing file: ", basename(archivo), " ==========")
  message(mensaje)
  nombre_tejido <- tools::file_path_sans_ext(basename(archivo))
  
  # Special condition for XMC/JX due to its name having an "/" in the name
  if (nombre_tejido == "XMC_JX") {
    nombre_tejido <- "XMC/JX" }
  
  # ---------------------------------------------------------------------------
  # 2.1 Load pseudocounts and clean
  # ---------------------------------------------------------------------------
  df <- read.csv(archivo, row.names = 1)
  
  df_numeric <- as.data.frame(lapply(df, function(x) as.numeric(trimws(x)))) # Transform to numeric and remove whitespace
  rownames(df_numeric) <- rownames(df)
  #df_numeric <- round(df_numeric) # With edgeR, rounding is NOT necessary, it can work with decimals
  
  df_numeric <- df_numeric %>%
    select(where(~ !all(.x == 0)))
  
  transposed_df_numeric <- t(df_numeric)
  
  # Remove Replicates (columns) with sum 0 of expression for all their genes (only happens in one celltype)
  transposed_df_numeric <- transposed_df_numeric[, colSums(transposed_df_numeric) > 0]
  
  # ---------------------------------------------------------------------------
  # 2.2 Define main condition treatment and time factor
  # ---------------------------------------------------------------------------
  
  # Extract Treatment (WT and Mutants)
  treatment <- factor(substring(colnames(transposed_df_numeric),1,2))
  treatment <- relevel(treatment, ref = "WT") # Set WT as reference for dropping later
  
  # Extract ZT + numbers (e.g. ZT8, ZT23)
  zt <- factor(sub(".*(ZT[0-9]+).*", "\\1", colnames(transposed_df_numeric)))
  
  genes <- rownames(transposed_df_numeric)
  
  # ---------------------------------------------------------------------------
  # 2.3 Create edgeR object
  # ---------------------------------------------------------------------------
  y <- DGEList(
    counts = transposed_df_numeric,
    genes  = genes,
    group  = treatment
  )
  
  # ---------------------------------------------------------------------------
  # 2.4 Filtering by expression and Normalization
  # ---------------------------------------------------------------------------
  
  keep <- filterByExpr(y) # Uses group by default
  table(keep) # Number of genes retained
  
  y <- y[keep, , keep.lib.sizes = F]
  
  ### TMM Normalization ###  
  y <- normLibSizes(y)
  #y$samples Normalization factor applied to each sample
  
  plotMDS(y, col= rep(1:6, each = 5)) # Observe that C9 is not completely grouped, variability will increase
  
  # Calculate mean cpms (wt and mutants) for heatmap plot
  cpms <- cpm(y,normalized.lib.sizes = T)
  cpms <- as.data.frame(cpms)
  is_wt <- grepl("WT", colnames(cpms))
  
  cpms$CPM_TMM_WT <- rowSums(cpms[, is_wt]) / sum(is_wt)
  cpms$CPM_TMM_MUTANT <- rowSums(cpms[, !is_wt]) / sum(!is_wt)
  cpms <- cpms %>% 
    mutate(Gene_Potri = rownames(cpms))
  
  # ---------------------------------------------------------------------------
  # 2.5 Define Additive Model Design for C1 and C9 respectively
  # ---------------------------------------------------------------------------
  
  # Design for C1
  
  y_C1 <- y[, treatment %in% c("WT","C1")]
  
  treatment_C1 <- droplevels(treatment[treatment %in% c("WT","C1")])
  zt_C1 <- droplevels(zt[treatment %in% c("WT","C1")])
  
  treatment_C1 <- relevel(treatment_C1, ref="WT")
  
  design_C1 <- model.matrix(~zt_C1 + treatment_C1)
  
  # Design for C9
  
  y_C9 <- y[, treatment %in% c("WT","C9")]
  
  treatment_C9 <- droplevels(treatment[treatment %in% c("WT","C9")])
  zt_C9 <- droplevels(zt[treatment %in% c("WT","C9")])
  
  treatment_C9 <- relevel(treatment_C9, ref="WT")
  
  design_C9 <- model.matrix(~zt_C9 + treatment_C9)
  
  
  # ---------------------------------------------------------------------------
  # 2.6 Dispersion Estimation for C1 and C9
  # ---------------------------------------------------------------------------
  
  y_C1 <- estimateDisp(y_C1, design_C1, robust=TRUE)
  fit_C1 <- glmQLFit(y_C1, design_C1, robust=TRUE)
  
  y_C9 <- estimateDisp(y_C9, design_C9, robust=TRUE)
  fit_C9 <- glmQLFit(y_C9, design_C9, robust=TRUE)
  
  # ---------------------------------------------------------------------------
  # 2.7 Differential Expression
  # ---------------------------------------------------------------------------
  if (nombre_tejido == "XMC/JX") { # To save correctly
    nombre_tejido <- "XMC_JX" }
  
  ### ZT and C1 ###
  zt_23_c1_qlf <- glmQLFTest(fit_C1, coef="treatment_C1C1")
  topTags(zt_23_c1_qlf)
  
  zt_genes_names <- zt_23_c1_qlf$genes
  zt_genes_data  <- zt_23_c1_qlf$table
  
  # Filter by FDR < 0.05
  zt_final_data <- merge(zt_genes_names, 
                         zt_genes_data, 
                         by ="row.names")
  
  zt_final_data <- zt_final_data %>% 
    dplyr::select(-Row.names) %>% 
    mutate(FDR = p.adjust(zt_final_data$PValue, method="BH")) %>% 
    filter(FDR < 0.05)
  
  # Apply poplars dictionary
  zt_final_data$Gene_Potri <- dict_poplars[zt_final_data$genes]
  
  # Add mean CPM calculations for heatmap (defined only once as it's the same for C1 and C9)
  cpms <- cpms %>% 
    dplyr::select(Gene_Potri,CPM_TMM_MUTANT,CPM_TMM_WT) %>% 
    dplyr::rename(genes = Gene_Potri) 
  
  zt_final_data <- inner_join(cpms,zt_final_data)
  
  write_csv(zt_final_data, 
            file = paste0("/Users/danielconde/Desktop/WATs/BayesnPrism/Results/Expression_Data_Bayes_Seurat_C1yC9/",
                          nombre_tejido,
                          "_Separated_C1_FDR_005_Aditive_Model.csv"))
  
  ### ZT and C9 ###
  zt_23_c9_qlf <- glmQLFTest(fit_C9, coef="treatment_C9C9")
  topTags(zt_23_c9_qlf)
  
  # Save genes affected by ZTs
  zt_genes_names <- zt_23_c9_qlf$genes
  zt_genes_data  <- zt_23_c9_qlf$table
  
  # Filter by FDR < 0.05
  zt_final_data <- merge(zt_genes_names, 
                         zt_genes_data, 
                         by ="row.names")
  
  zt_final_data <- zt_final_data %>% 
    dplyr::select(-Row.names) %>% 
    mutate(FDR = p.adjust(zt_final_data$PValue, method="BH")) %>% 
    filter(FDR < 0.05)
  
  # Apply poplars dictionary
  zt_final_data$Gene_Potri <- dict_poplars[zt_final_data$genes]
  
  # Add mean CPM calculations for heatmap
  zt_final_data <- inner_join(cpms,zt_final_data)
  
  write_csv(zt_final_data, 
            file = paste0("/Users/danielconde/Desktop/WATs/BayesnPrism/Results/Expression_Data_Bayes_Seurat_C1yC9/",
                          nombre_tejido,
                          "_Separated_C9_FDR_005_Aditive_Model.csv"))
  
  
  message("✓ Saved: ", nombre_tejido)
}

message("\n=== PROCESSING COMPLETED FOR ALL FILES ===\n")
