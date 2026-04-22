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
# - Phytozome P.tremula x alba HAP2 v5.1 annotations (to obtain KOs)
# Output files:
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

### Path for file input ###
ruta_directorio <- "/Users/juanmurillomurillo/Library/CloudStorage/GoogleDrive-lab171@intranet.cbgp.upm.es/Mi unidad/Drive_Juan/WATs/Pathways_Analysis_Mutantes_ZT8_ZT23/"

### Mutant data ZT8_C1 vs WT & ZT8_C9 vs WT 
ruta_destino <- paste0(ruta_directorio,"Results/KEGGs/Nuevo_GSEA_KEGG/")
c1 <- read.csv(paste0(ruta_directorio,"Data/Datos edgeR/Datos_Bulk_edgeR_Zt_como_Covariable/Separated_C1_Aditive_Model.csv"))
c9 <- read.csv(paste0(ruta_directorio,"Data/Datos edgeR/Datos_Bulk_edgeR_Zt_como_Covariable/Separated_C9_Aditive_Model.csv"))

setwd('/Users/juanmurillomurillo/Library/CloudStorage/GoogleDrive-lab171@intranet.cbgp.upm.es/Mi unidad/Drive_Juan/WATs/Pathways_Analysis_Mutantes_ZT8_ZT23/Results/KEGGs/Nuevo_GSEA_KEGG')

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

# SELECTION OF KO TERMS TO USE
# Phytozome file with Populus tremula x alba HAP2 annotations
ruta <- paste0(ruta_directorio,"Data/Populus Tremula x Alba HAP2/Phytozome/PhytozomeV14/PtremulaxPopulusalbaHAP2/v5.1/annotation/PtremulaxPopulusalbaHAP2_716_v5.1.P14.annotation_info.txt/PtremulaxPopulusalbaHAP2_716_v5.1.P14.annotation_info.txt")
poplar <- read_tsv(ruta)
poplar <- poplar %>% 
  dplyr::rename(Gene = locusName)

# =============================================================================
# 2. PREPARATION OF edgeR EXPRESSION DATA
# =============================================================================

# Generation of the combined dataset (Gene Name, log2FC, p_value)
# Perhaps they should have kept the same sign for logFC
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
# 3. MAPPING Gene - KOs - PATHWAY FOR THE SELECTED SPECIES
# =============================================================================

# Merge annotations and genes (keep only genes with annotations)
poplar_c1_c9 <- inner_join(poplar,
                           c1_c9,
                           by = "Gene")

gsea_data <- poplar_c1_c9 %>% 
  dplyr::select(Gene, 
                ranking)

# Obtain K number mappings (KXXXXX -> KO:XXXXX; KEGG ORTHOLOGY)
Gene_KO <- poplar_c1_c9 %>%
  select(Gene, KO) %>%
  mutate(KO = str_split(KO, "[,;\\\\s]+")) %>%
  unnest(KO) %>%
  filter(!is.na(KO), KO != "") %>%
  distinct(Gene, KO)

KOs <- Gene_KO$KO %>% unique()          

# 1. Map KOs to species genes (Populus alba "palz" in this case)
# Helper function to split into batches of 10
chunk_vec <- function(v, n = 10) {
  split(v, ceiling(seq_along(v) / n))
}

# 1. Map KOs to species genes (Populus alba = "palz")
ko_chunks <- chunk_vec(KOs, 10)

ko_to_gene <- lapply(ko_chunks, function(chunk) {
  Sys.sleep(0.1)  # small pause to avoid 403 error
  keggLink("palz", paste0("ko:", chunk))
})

palz_genes <- unique(sub("^palz:", "", unname(unlist(ko_to_gene))))

# 2. Map Genes -> pathways (remember that a generic mapping exists; organism = "ko")
palz_gene_paths <- bitr_kegg(palz_genes, fromType="kegg", toType="Path", organism="palz")

# 3. Obtain unique pathway names
palz_unique_paths <- unique(palz_gene_paths$Path)

# 4. Obtain KO:xxxxx definitions
# Split into batches because this is the maximum accepted by keggGet
chunks <- split(palz_unique_paths, ceiling(seq_along(palz_unique_paths)/1))

# Apply the function
path_annots <- lapply(chunks, function(p) {
  Sys.sleep(0.1) 
  keggGet(paste0("path:", p))
})

path_annots_2 <- unlist(path_annots, recursive = FALSE)

path_df <- do.call(bind_rows, lapply(path_annots_2, function(xx) {
  data.frame(
    Path             = if (!is.null(xx$ENTRY)) xx$ENTRY else NA_character_,
    Path_name        = if (length(xx$NAME)) xx$NAME[1] else NA_character_,
    Path_description = if (!is.null(xx$DESCRIPTION)) paste(xx$DESCRIPTION, collapse = " ") else NA_character_,
    stringsAsFactors = FALSE
  )
}))

# path_df comes from keggList("pathway"): Path = mapXXXXX, Path_name
kegg_ko_desc <- merge(palz_gene_paths, path_df,
                      by = "Path", 
                      all.x = F)

ko_to_gene_df <- do.call(rbind, lapply(ko_to_gene, function(x) {
  if (length(x) == 0) return(NULL)
  data.frame(
    KO = sub("^ko:", "", names(x)),          # remove 'ko:' prefix
    palz_gene_full = unname(x),              # returned value: 'palz:XXXXX'
    kegg = sub("^palz:", "", unname(x)) # without prefix
  )
}))

gene_kegg_ko_desc <- merge(ko_to_gene_df, 
                           kegg_ko_desc,
                           by = "kegg", 
                           all.x = F) %>% dplyr::select(-NAME)

# =============================================================================
# 4. PREPARATION OF GSEA INPUT: GENELIST | TERM2GENE | TERM2NAME
# =============================================================================

# Generate geneList, term2name and term2gene (as indicated in the clusterProfiler documentation)
# GeneList
geneList <- gsea_data$ranking
names(geneList) <- gsea_data$Gene

geneList <- geneList[ names(geneList) != "" ]               # remove empty names
geneList <- geneList[ !duplicated(names(geneList)) ]        # remove duplicates (by name)
geneList <- sort(geneList, decreasing = TRUE)               # sort (mandatory)

# 1) TERM2GENE: pathway (koXXXXX) ↔ gene (your IDs)
term2gene <- merge(Gene_KO,gene_kegg_ko_desc, by = "KO") %>% 
  dplyr::select(Path, Gene)

# 2) TERM2NAME: pathway (koXXXXX) ↔ name
term2name <- merge(Gene_KO,gene_kegg_ko_desc, by = "KO") %>% 
  dplyr::select(Path,Path_name)

# =============================================================================
# 5. APPLY GSEA TO THE UPPER (POS) AND LOWER (NEG) TAIL SEPARATELY
# =============================================================================

for (scoreType in c("pos", "neg")) {
  
  ego <- GSEA(
    geneList   = geneList,      
    TERM2GENE  = term2gene,
    TERM2NAME  = term2name,
    minGSSize  = minGSSize, 
    maxGSSize = maxGSSize,
    pvalueCutoff = 1, 
    verbose = FALSE, 
    eps = 0,
    scoreType = scoreType,
    seed = seed
  )
  
  # Save results
  ego_csv <- as.data.frame(ego)
  ego_csv <- ego_csv %>% 
    separate_rows(core_enrichment, sep = "/") %>% 
    rename(Alba = core_enrichment)
  
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