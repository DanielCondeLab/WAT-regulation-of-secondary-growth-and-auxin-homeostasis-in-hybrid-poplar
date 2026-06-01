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
# - Populus Dicctionary Mapping (P.tremula x alba HAP2 - P.Alba from KEGG - KEGG ID)
#   - KEGG_Mapping_Info.csv
# Output files:
# - Pathview KEGG pathway image (PNG):
#   - Generated automatically by pathview (e.g., pathway 00400)
############################################################

# =============================================================================
# 1. LOAD DATA AND LIBRARIES
# =============================================================================

library(pathview)
library(tidyverse)
library(clusterProfiler)
library(KEGGREST)
library(svglite)

# Seleccion de mapa a visualizar
map_to_plot <- "00400" # Code of the KEGG pathway -> e.g. (00562 Inositol phosphate metabolism)

### Base directory ###
ruta_directorio <- "/Users/danielconde/Library/CloudStorage/GoogleDrive-lab171@intranet.cbgp.upm.es/Mi unidad/Drive_Juan/WATs/Pathways_Analysis_Mutantes_ZT8_ZT23/"
ruta_destino <- paste0(ruta_directorio,"Results/KEGGs/Pathview/")
setwd(ruta_destino)

### Input differential expression data ###
c1 <- read.csv(paste0(ruta_directorio,"Data/Datos edgeR/Datos_Bulk_edgeR_Zt_como_Covariable/Separated_C1_Aditive_Model.csv"))
c9 <- read.csv(paste0(ruta_directorio,"Data/Datos edgeR/Datos_Bulk_edgeR_Zt_como_Covariable/Separated_C9_Aditive_Model.csv"))

# Filtered (FDR < 0.05)
c1_c9_FDR_005 <- inner_join(c1,c9, by = "genes", suffix = c("_1", "_9")) %>% 
  filter(FDR_1 < 0.05 & FDR_9 < 0.05) %>% 
  mutate(logFC_Mean = (logFC_1 + logFC_9) /2 )

### Archivo Diccionario KEGG
poplar_dictionary <- read_csv(paste0(ruta_directorio,'Results/KEGGs/Nuevo_GSEA_KEGG/KEGG_Mapping_Info.csv'))

### Genes to Visualize
joined_data <- inner_join(c1_c9_FDR_005 %>% rename(Gene = genes),poplar_dictionary)
gene_vec_pref <- setNames(joined_data$logFC_Mean, joined_data$kegg)

# =============================================================================
# 4. PATHVIEW VISUALIZATION (SIGNIFICANT GENES ONLY)
# =============================================================================

### Generate KEGG native map (PNG)
pv <- pathview (
  gene.data   = gene_vec_pref,
  pathway.id  =  map_to_plot,        
  species     = "palz",     # Code for Populus Alba
  gene.idtype = "kegg",
  kegg.native = T # T KEGG format | F Graphviz Format
)
pv_alba <-as.data.frame(pv$plot.data.gene)
