############################################################
# CONDA ENVIRONMENT: Pathways_analysis.yaml
# Pathview visualization of KEGG pathways of interest for mutants C1 and C9
# Visualize gene expression changes (logFC direction) for significant genes (FDR < 0.05) on selected KEGG pathway maps.
# NOTE: Pathways of interest should be selected based on prior GSEA results.
# File routes MUST be changed in order to use this code.
# Input files:
# - Differential expression results for ALL genes:
#   - Separated_C1_Aditive_Model.csv
#   - Separated_C9_Aditive_Model.csv
# - Differential expression results filtered (FDR < 0.05):
#   - Separated_C1_FDR_005_Aditive_Model.csv
#   - Separated_C9_FDR_005_Aditive_Model.csv
# - Phytozome annotation file (Populus tremula x alba HAP2)
# - Ortholog mapping file (OrthoFinder one-to-one pairs using P.alba from KEGG)
# Output files:
# - Table with P.tremula x alba HAP2 gene → P.alba KEGG ID mapping:
#   - P.trem_alba_hap2_P.alba_kegg_zt8.csv
# - Pathview KEGG pathway image (PNG):
#   - Generated automatically by pathview (e.g., pathway 00400)
# - Table linking genes, KEGG IDs, and pathway visualization results (e.g. for Phenylalanine, tyrosine and tryptophan biosynthesis):
#   - Map_00400.csv
############################################################

# =============================================================================
# 1. LOAD DATA AND LIBRARIES
# =============================================================================

library(pathview)
library(tidyverse)
library(clusterProfiler)
library(KEGGREST)
library(svglite)

# Select KO species
Species_KOs <- "Poplar" # Options: Poplar / Arab

### Base directory ###
ruta_directorio <- "/Users/danielconde/Desktop/WATs/Pathways_Analysis_Mutantes_ZT8_ZT23/"

### Input differential expression data ###
ruta_destino <- paste0(ruta_directorio,"Results/KEGGs/Nuevo_GSEA_KEGG/")

c1 <- read.csv(paste0(ruta_directorio,"Data/Datos edgeR/Datos_Bulk_edgeR_Zt_como_Covariable/Separated_C1_Aditive_Model.csv"))
c9 <- read.csv(paste0(ruta_directorio,"Data/Datos edgeR/Datos_Bulk_edgeR_Zt_como_Covariable/Separated_C9_Aditive_Model.csv"))

# Filtered (FDR < 0.05)
c1_fdr <- read.csv(paste0(ruta_directorio,"Data/Datos edgeR/Datos_Bulk_edgeR_Zt_como_Covariable/Separated_C1_FDR_005_Aditive_Model.csv")) %>% 
  rename(logFC_C1 = logFC, FDR_C1 = FDR)

c9_fdr <- read.csv(paste0(ruta_directorio,"Data/Datos edgeR/Datos_Bulk_edgeR_Zt_como_Covariable/Separated_C9_FDR_005_Aditive_Model.csv")) %>% 
  rename(logFC_C9 = logFC, FDR_C9 = FDR)

setwd(ruta_destino)

# Rename "genes" column to "Gene" if present
objs <- ls(envir = .GlobalEnv)

for (obj in objs) {
  x <- get(obj, envir = .GlobalEnv)
  if ("genes" %in% colnames(x)) {
    x <- x %>% rename(Gene = genes)
    assign(obj, x, envir = .GlobalEnv)
  }
}

# =============================================================================
# 2. KO ANNOTATION SELECTION
# =============================================================================

if (Species_KOs == "Poplar") {
  
  ruta <- paste0(ruta_directorio,"Data/Populus Tremula x Alba HAP2/Phytozome/PhytozomeV14/PtremulaxPopulusalbaHAP2/v5.1/annotation/PtremulaxPopulusalbaHAP2_716_v5.1.P14.annotation_info.txt/PtremulaxPopulusalbaHAP2_716_v5.1.P14.annotation_info.txt")
  
  poplar <- read_tsv(ruta) %>% 
    rename(Gene = locusName)
  
} else if (Species_KOs == "Arab") {
  
  poplar <- read_csv("/Users/danielconde/Desktop/WATs/Pathways_Analysis_Mutantes_ZT8_ZT23/Data/Populus-Arabidposis-KO/Populus-Arabidposis-KO.csv") %>% 
    rename(Gene = pop_gene)
  
} else {
  stop("Please provide a valid species: 'Poplar' or 'Arab'.")
}

# =============================================================================
# 3. PREPARE EXPRESSION DATA
# =============================================================================

c1_c9 <- inner_join(
  c1 %>% select(Gene, logFC, PValue),
  c9 %>% select(Gene, logFC, PValue),
  by = "Gene",
  suffix = c("_C1", "_C9")
)

# Fisher combined p-value
c1_c9 <- c1_c9 %>%
  mutate(
    chisq = -2 * (log(PValue_C1) + log(PValue_C9)),
    p_fisher = pchisq(chisq, df = 4, lower.tail = FALSE),
    adj_p_fisher = p.adjust(p_fisher, method = "BH")
  )

# Ranking statistic
c1_c9 <- c1_c9 %>%
  mutate(
    log2FC_avg = (logFC_C1 + logFC_C9) / 2,
    ranking    = -log10(p_fisher) * sign(log2FC_avg)
  ) %>%
  distinct(Gene, .keep_all = TRUE)

# =============================================================================
# 4. PATHVIEW VISUALIZATION (SIGNIFICANT GENES ONLY)
# =============================================================================
# Import 1:1 pairs from OrthoFinder
orthologues_poplar <- read_tsv('/Users/danielconde/Desktop/WATs/Pathways_Analysis_Mutantes_ZT8_ZT23/Results/Alineamiento_P.Alba_x_P.trem_x_alba/OrthoFinder/Results_Dec29/one2one_pairs.tsv')

# Remove isoform suffixes
orthologues_poplar <- orthologues_poplar %>%
  mutate(
    HAP2_id = gsub("\\.1(\\.p)?$", "", HAP2_id),
    ALBA_id = gsub("\\.1(\\.p)?$", "", ALBA_id)
  ) %>% 
  dplyr::select(HAP2_id, ALBA_id) %>% 
  rename(Gene = HAP2_id)

# Download KEGG -> NCBI protein mapping for all 'palz'
map_all <- keggConv("ncbi-proteinid", "palz") # map_all: named vector; names = "palz:NNNNNNN", values = "ncbi-proteinid:XP_..."

# Reverse mapping: RefSeq protein -> KEGG
inv_map <- setNames(names(map_all), sub("^ncbi-proteinid:", "", unname(map_all))) # inv_map["XP_073263267"] -> "palz:105153908" (example)

# Apply to ZT8 and ZT23 data
c1_c9_fdr <- inner_join(
  c1_fdr  %>% dplyr::select(Gene, logFC_C1, FDR_C1),
  c9_fdr %>% dplyr::select(Gene, logFC_C9, FDR_C9),
  by = "Gene"
)

# Calculate average logFC for every significant gene
c1_c9_fdr <- c1_c9_fdr %>% 
  mutate(logFC_Average = (logFC_C1 + logFC_C9) /2)

filtered_c1_c9 <- c1_c9_fdr %>% 
  dplyr::filter(FDR_C1 <= 0.05 & FDR_C9 <= 0.05) 

pathview_data <- inner_join(orthologues_poplar, filtered_c1_c9, by = "Gene") %>% 
  dplyr::select(Gene, ALBA_id, logFC_Average)

palba_ids <- unique(na.omit(pathview_data$ALBA_id))

map_df <- tibble(
  ALBA_id = palba_ids,
  KEGG_id = sub("^palz:", "", unname(inv_map[palba_ids]))
)

pathview_ready <- pathview_data %>% 
  left_join(map_df, by = "ALBA_id")

write_csv(pathview_ready, 
          file = "/Users/danielconde/Desktop/WATs/Pathways_Analysis_Mutantes_ZT8_ZT23/Results/KEGGs/Nuevo_GSEA_KEGG/P.trem_alba_hap2_P.alba_kegg_zt8.csv")

gene_vec_pref <- setNames(pathview_ready$logFC_Average, pathview_ready$KEGG_id)

### Generate KEGG native map (PNG)
pv <- pathview (
  gene.data   = gene_vec_pref,
  pathway.id  =  "00400",   # Code of the KEGG pathway         
  species     = "palz",     # Code for Populus Alba
  gene.idtype = "kegg",
  out.dir     = ruta_destino
)

pv_alba <-as.data.frame(pv$plot.data.gene)

#### FINAL MERGE: Poplar gene - KEGG gene name - Pathview results

join_results_pathview <- pv_alba

join_results_pathview <- join_results_pathview %>% 
  separate_rows(all.mapped,sep = ",") %>% 
  dplyr::rename(KEGG_id = all.mapped)

final_join <- inner_join(join_results_pathview, pathview_ready, by ="KEGG_id")

final_join <- final_join %>% 
  relocate(Gene, ALBA_id, KEGG_id)

write_csv(final_join, file = '/Users/danielconde/Desktop/WATs/Pathways_Analysis_Mutantes_ZT8_ZT23/Results/KEGGs/Nuevo_GSEA_KEGG/Map_00400.csv')