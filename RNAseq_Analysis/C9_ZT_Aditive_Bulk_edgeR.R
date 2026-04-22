###########################################################
# CONDA ENVIRONMENT: Bayes.yaml
# This script performs bulk RNA-seq differential expression analysis with edgeR
# separately for C1 and C9 (common filtering and normalization for both mutants).
# We have WAT mutants (C1 and C9) vs WT and two time points: ZT8 and ZT23.
# We do NOT want to study temporal effects directly, but because many genes vary significantly in expression
# across time, we must account for the effect of ZT when performing differential expression analysis.
# This specific script section contains the analysis and output generation for WT vs C9.
# File routes MUST be changed in order to use this code.
# Input files:
# - Bulk raw count matrix for 30 samples
# Output files:
# - One .csv file with all WT vs C9 differential expression results from the additive model
# - One .csv file with WT vs C9 differential expression results filtered at FDR < 0.05
# edgeR manual: https://bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
############################################################

# =============================================================================
# IMPORT LIBRARIES
# =============================================================================

library(edgeR)
library(limma)
library(tidyverse)

# =============================================================================
# LOAD DATA
# =============================================================================

rawdata <- read_csv("/Users/danielconde/Desktop/WATs/Bulk_edgeR/Data/Counts_30samples.csv") %>%
  column_to_rownames(var = colnames(.)[1])

genes = rownames(rawdata)

# =============================================================================
# CREATE GROUPINGS BY TREATMENT (WT AND C9) AND ZTs (ZT8 AND ZT23)
# =============================================================================

# Extract treatment (WT and mutants)
treatment <- factor(substring(colnames(rawdata),1,2))
treatment <- relevel(treatment, ref = "WT") # Set WT as reference to drop it later

# Extract ZT + numbers (e.g. ZT8, ZT23)
zt <- factor(sub(".*(ZT[0-9]+).*", "\\1", colnames(rawdata)))

# =============================================================================
# MAKE DGEList OBJECT
# =============================================================================

y <- DGEList(counts = rawdata, 
             genes = genes, 
             group = treatment)

# =============================================================================
# FILTERING AND NORMALIZATION
# =============================================================================

keep <- filterByExpr(y) # We include design; include the matrix in design
table(keep) # Number of retained genes

y <- y[keep, , keep.lib.sizes = FALSE]

### TMM NORMALIZATION ###

y <- normLibSizes(y)
# y$samples contains the normalization factor applied to each sample

# =============================================================================
# DATA EXPLORATION
# =============================================================================

plotMDS(y, col = rep(1:6, each = 5)) # We see that C9 is not completely grouped, which will increase variability

### model.matrix design
interaction_design <- model.matrix(~zt + zt:treatment) # Rows are samples
logFC <- predFC(y,
                interaction_design,
                prior.count = 1,
                dispersion = 0.05)

cor(logFC[,3:6])

# In this case, the observed correlations indicate that C9 shows relatively consistent behavior across time points,
# whereas C1 shows greater variability, which reinforces the convenience of a design model, especially for C1.
# However, drastic changes in the biological direction of the results are not expected compared with the basic exactTest analysis.

# =============================================================================
### DESIGN MATRIX ###
# =============================================================================

y_C9 <- y[, treatment %in% c("WT","C9")]

treatment_C9 <- droplevels(treatment[treatment %in% c("WT","C9")])
zt_C9 <- droplevels(zt[treatment %in% c("WT","C9")])

treatment_C9 <- relevel(treatment_C9, ref = "WT")

design_C9 <- model.matrix(~zt_C9 + treatment_C9)

# =============================================================================
# DISPERSION ESTIMATION
# =============================================================================
y_C9 <- estimateDisp(y_C9, design_C9, robust = TRUE)
fit_C9 <- glmQLFit(y_C9, design_C9, robust = TRUE)

plotBCV(y_C9)
plotQLDisp(fit_C9)

# The common dispersion is the square root of the BCV. Very, very low values indicate that the experiment is not reliable because there is
# almost no variability. Very high values indicate high variance between replicates, which reduces the power to detect DE genes.
# If “Squeezed” (red) does not change much relative to “Raw” (black), it usually indicates that there were no genes with extreme QL
# dispersions requiring substantial shrinkage.
# The goal of the plot is not for them to be identical, but for the squeeze to reduce noise without distorting the global structure.

# =============================================================================
# DIFFERENTIAL EXPRESSION
# =============================================================================

#####################  ZT AND C9 #####################
zt_23_c9_qlf <- glmQLFTest(fit_C9, coef = "treatment_C9C9")
topTags(zt_23_c9_qlf)

# Store genes affected by ZT (2 files: unfiltered and FDR-filtered)
zt_genes_names <- zt_23_c9_qlf$genes
zt_genes_data  <- zt_23_c9_qlf$table

# Unfiltered
zt_final_data <- merge(zt_genes_names, zt_genes_data, by = "row.names")
zt_final_data <- zt_final_data %>% 
  dplyr::select(-Row.names) %>% 
  mutate(FDR = p.adjust(zt_final_data$PValue, method = "BH")) 

write_csv(zt_final_data, 
          file = '/Users/danielconde/Desktop/WATs/Bulk_edgeR/Results/C9/Aditive/Separated_C9_Aditive_Model.csv')

# Filtered at FDR < 0.05
zt_final_data <- merge(zt_genes_names, 
                       zt_genes_data, 
                       by = "row.names")

zt_final_data <- zt_final_data %>% 
  dplyr::select(-Row.names) %>% 
  mutate(FDR = p.adjust(zt_final_data$PValue, method = "BH")) %>% 
  filter(FDR < 0.05)

write_csv(zt_final_data, 
          file = '/Users/danielconde/Desktop/WATs/Bulk_edgeR/Results/C9/Aditive/Separated_C9_FDR_005_Aditive_Model.csv')
