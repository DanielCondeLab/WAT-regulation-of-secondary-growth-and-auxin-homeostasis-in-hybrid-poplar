################################################
# CONDA ENVIRONMENT: Pathways_analysis.yaml
# This script performs Gene Set Enrichment Analysis (GSEA) using differential expression data and 
# applying one-tailed tests (positive and negative part).
# Step 1: Generate a new p-value using Fisher's method and the mean log2FC of of C1 vs WT and C9 vs WT data.
# Step 2: Prepare the ranking statistic: ranking = -log10(p_value) * sign(log2FC_avg).
# Step 3: Sort from highest to lowest based on the ranking statistic and apply GSEA.
# Step 4: Extract the genes associated with the pathways of interest.
# File routes MUST be changed in order to use this code.
# Input files:
# - Differential expression data for ALL genes from C1 vs WT and C9 vs WT
# - Best Reciprocal Hit (BRH) of Phytozome P.tremula x alba HAP2 v5.1 - Populus Alba from KEGG
# Output files:
# - Populus Dicctionary Mapping (P.tremula x alba HAP2 - P.Alba from KEGG - KEGG ID)
#   - KEGG_Mapping_Info.csv
# - .csv files with all enriched KEGG pathways and their core enrichment genes, separately for positive and negative one-tailed tests
# - DotPlots of the top 20 enriched KEGG pathways filtered by adjusted p-value , separately for positive and negative one-tailed tests
# - Network plots (cnetplot) of the top 20 enriched KEGG pathways filtered by adjusted p-value, separately for positive and negative one-tailed tests
# - Simplified DotPlots with cleaned pathway descriptions, separately for positive and negative one-tailed tests
################################################

# =============================================================================
# 1. IMPORT DATA AND LIBRARIES
# =============================================================================

# Import libraries
library(pathview)
library(tidyverse)
library(clusterProfiler)
library(KEGGREST)
library(svglite)

# Parameters for GSEA
set.seed(15)
seed = 15
minGSSize = 10  # Minimun number of genes required for the category
maxGSSize = 500 # Maximum number of genes required for the category
pvalueCutoff = 0.05 

### Paths for file input ###
ruta_directorio <- "/Users/danielconde/Library/CloudStorage/GoogleDrive-lab171@intranet.cbgp.upm.es/Mi unidad/Drive_Juan/WATs/Pathways_Analysis_Mutantes_ZT8_ZT23/"
ruta_destino <- paste0(ruta_directorio,"Results/KEGGs/Nuevo_GSEA_KEGG/")

### Mutant data ZT8_C1 vs WT & ZT8_C9 vs WT 
c1 <- read.csv(paste0(ruta_directorio,"Data/Datos edgeR/Datos_Bulk_edgeR_Zt_como_Covariable/Separated_C1_Aditive_Model.csv"))
c9 <- read.csv(paste0(ruta_directorio,"Data/Datos edgeR/Datos_Bulk_edgeR_Zt_como_Covariable/Separated_C9_Aditive_Model.csv"))

setwd('/Users/danielconde/Library/CloudStorage/GoogleDrive-lab171@intranet.cbgp.upm.es/Mi unidad/Drive_Juan/WATs/Pathways_Analysis_Mutantes_ZT8_ZT23/Results/KEGGs/Nuevo_GSEA_KEGG')

# Rename column "genes" to "Gene"
# Get object names from the global environment
objs <- ls(envir = .GlobalEnv)

for (obj in objs) {
  # Get the object
  x <- get(obj, envir = .GlobalEnv)
  # Only act if it is a data.frame and contains the column
  if ("genes" %in% colnames(x)) {
    x <- x %>% rename(Gene = genes)
    # Overwrite the original object
    assign(obj, x, envir = .GlobalEnv)
  }
}

### Load Datasets With BRH of P.tremula x alba HAP2 Phytozome - P.alba KEGG
poplar_anotation_route <- '/Users/danielconde/Library/CloudStorage/GoogleDrive-lab171@intranet.cbgp.upm.es/Mi unidad/Drive_Juan/WATs/Pathways_Analysis_Mutantes_ZT8_ZT23/Results/Alineamiento_P.Alba_x_P.trem_x_alba/PtXaAlbH_Palba_RBH.csv'
poplar <- read_csv(poplar_anotation_route)
poplar <- poplar %>% 
  dplyr::select(1,2) %>% 
  rename(Gene = PtXaAlbH)

# =============================================================================
# 2. PREPARATION OF edgeR EXPRESSION DATA AND RANKING STATISTIC
# =============================================================================

# Generation of the combined dataset (Gene Name, log2FC, p_value)
c1_c9 <- inner_join(
  c1 %>% dplyr::select(Gene, logFC, PValue),
  c9 %>% dplyr::select(Gene, logFC, PValue),
  by = "Gene",
  suffix = c("_C1", "_C9")
)

# Generation of p-value using Fisher's method (independence of p-values)
c1_c9  <- c1_c9  %>%
  mutate(
    # Fisher's chi-square statistic: -2 sum log(p)
    chisq = -2 * (log(PValue_C1) + log(PValue_C9)),
    # combined p-value with df = 2*k = 4
    p_fisher = pchisq(chisq, df = 4, lower.tail = FALSE),
    # FDR adjustment
    adj_p_fisher = p.adjust(p_fisher, method = "BH")
  )

# Data for GSEA
c1_c9 <- c1_c9 %>%
  mutate(
    log2FC_avg = (logFC_C1 + logFC_C9) / 2,  # Mean log2FC
    ranking    = -log10(p_fisher) * sign(log2FC_avg) # Ranking statistic
  ) %>%
  distinct(Gene, .keep_all = TRUE)     # no duplicates per gene

# =============================================================================
# 3. MAPPING Palba KEGG Gene - KEGG Gene Number ID  - PATHWAYS MAPPs FOR THE SELECTED SPECIES
# =============================================================================

# Merge annotations and genes (keep only genes with annotations)
gsea_data  <- inner_join(poplar,
                         c1_c9,
                         by = "Gene")

# Keep P.alba KEGG Gene Name and Ranking Statistic
gsea_data <- gsea_data %>% 
  dplyr::select(Palba, 
                ranking)

# Download KEGG P.alba Pathways Info
kegg_palz_info <- keggLink("palz","pathway")
kegg_palz_info_df <- data.frame(
  Path = names(kegg_palz_info),
  Gene = unname(kegg_palz_info),
  stringsAsFactors = FALSE
)

# Eliminate palz prefix
palz_genes <- unique(sub("^palz:", "", kegg_palz_info_df$Gene))

# Map Genes -> KEGG ID Proteina ID Gen 
palz_gene_paths_and_kegg_ids <- bitr_kegg(palz_genes, fromType="kegg", toType="ncbi-proteinid", organism="palz")
palz_gene_paths_and_kegg_ids <- palz_gene_paths_and_kegg_ids %>% 
  rename(Palba = `ncbi-proteinid`)

# Join P.alba from KEGG - KEGG ID - Ranking Statistic
final_gsea_dataset <- inner_join(gsea_data, palz_gene_paths_and_kegg_ids)

# Saving Gene Mapping Info
mapping_all_info <- inner_join(final_gsea_dataset, poplar) %>% 
  relocate(Gene, Palba, kegg)

write_csv(mapping_all_info, paste0(ruta_destino, "KEGG_Mapping_Info.csv"))


# =============================================================================
# 4. PREPARATION OF GSEA INPUT: GENELIST
# =============================================================================
# Keep P.alba KEGG Id 
final_gsea_dataset <- final_gsea_dataset %>% 
  select(kegg, ranking)

# Generate geneList, term2name and term2gene (as indicated in the clusterProfiler documentation)
# GeneList
geneList <- final_gsea_dataset$ranking
names(geneList) <- final_gsea_dataset$kegg

geneList <- geneList[ names(geneList) != "" ]               # remove empty names
geneList <- geneList[ !duplicated(names(geneList)) ]        # remove duplicates (by name)
geneList <- sort(geneList, decreasing = TRUE)               # sort (mandatory)

# =============================================================================
# 5. APPLY GSEA TO THE UPPER (POS) AND LOWER (NEG) TAIL SEPARATELY
# =============================================================================

for (scoreType in c("neg", "pos")) {
  
  ego <-  gseKEGG(
    geneList,
    organism = "palz",
    keyType = "kegg",
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    pvalueCutoff = 1,
    pAdjustMethod = "fdr",
    verbose = T,
    use_internal_data = F,
    seed = seed,
    scoreType = scoreType
  )
  
  
  # Save results
  ego_csv <- as.data.frame(ego)
  ego_csv <- ego_csv %>% 
    separate_rows(core_enrichment, sep = "/") %>% 
    rename(kegg = core_enrichment) 
  
  ego_csv <- inner_join(ego_csv, mapping_all_info %>% select(-ranking))
  
  write_csv(ego_csv,
            file = paste0(ruta_destino,"/", toupper(scoreType), "_All_Enriched_KEGG_Pathways.csv" ))
  
  ### Plots ###
  
  # DotPlot P.adjust
  df <- as.data.frame(ego)
  df_filtered <- df[df$p.adjust <= 0.05, ]   # or qvalues
  ego_q <- ego
  ego_q@result <- df_filtered
  
  q_dot_plot <- dotplot(ego_q,
                        x = "Count",
                        color = "p.adjust",
                        showCategory = 20)
  
  ggsave(filename = paste0(ruta_destino,"/", scoreType, "_P.adjust_Top_20_Enriched_KEGG_Pathways.svg"),
         plot = q_dot_plot,
         width = 14, 
         height = 14)
  
  
  # Network Plot P.adjust
  ego_network <- ego
  ego_network@result <- df_filtered
  ego_network@result$Description <- ego_network@result$Description %>% 
    str_remove("\\s*- Populus alba \\(white poplar\\)") %>%
    str_trim()
  
  ego_network@result$Description <- paste0(
    ego_network@result$ID, "\n", ego_network@result$Description
  )
  
  network_plot <- cnetplot(ego_network,
                           showCategory = 20,
                           circular = FALSE,
                           colorEdge = TRUE,
                           size_item = 0.4,
                           color_category = "red",
                           node_label = "category")
  
  ggsave(filename = paste0(ruta_destino,"/", scoreType, "_NetworkPlot_P.adjust_Top_20_Enriched_KEGG_Pathways.svg"),
         plot = network_plot,
         width = 14, 
         height = 14,
         dpi = "retina")
  
  # DotPlot P.adjust with simplified Populus alba text
  adjusted_dotplot  <- dotplot(ego_network,
                               x = "Count",
                               color = "p.adjust",
                               showCategory = 20)
  
  ggsave(filename = paste0(ruta_destino,"/", scoreType, "_Simplified_P.adjust_Top_20_Enriched_KEGG_Pathways.svg"),
         plot = adjusted_dotplot,
         width = 14, 
         height = 14)
}