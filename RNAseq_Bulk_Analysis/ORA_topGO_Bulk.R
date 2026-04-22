############################################
# CONDA ENVIRONMENT: Pathways_analysis.yaml
# Over Representation Analysis (ORA) enrichment of UP- and DOWN-regulated DEGs (logFC > 0; logFC < 0; analyzed separately) 
# from significant bulk genes (FDR < 0.05) shared between C1 vs WT and C9 vs WT
# using Biological Process Gene Ontology (BP; GO)
# File routes MUST be changed in order to use this code.
# Input files:
# - Differential expression results for C1 vs WT from bulk RNA-seq (used to identify shared UP and DOWN DEGs)
# - Differential expression results for C9 vs WT from bulk RNA-seq (used to identify shared UP and DOWN DEGs)
# - Phytozome Populus tremula x alba HAP2 Database: P. tremula x alba HAP2 gene | Arabidopsis GO terms
# - Arabidopsis gene functional descriptions from TAIR database
# Output files:
# - One .svg DotPlot showing the top 20 enriched GO Biological Process terms for the selected direction (up or down)
# - One .csv file containing enriched GO terms with contributing genes and Arabidopsis functional information for the selected direction (up or down)
# topGO manual: https://bioconductor.org/packages/release/bioc/manuals/topGO/man/topGO.pdf
############################################

# =============================================================================
# 0. IMPORT LIBRARIES
# =============================================================================

library(topGO)
library(tidyverse)

# =============================================================================
# 1. IMPORT DATA AND SELECT UP OR DOWN
# =============================================================================

### Select data direction: UP or DOWN
up_or_down <- "down" # up or down
###

### Path for file input ###
ruta_directorio <- "/Users/danielconde/Desktop/WATs/Pathways_Analysis_Mutantes_ZT8_ZT23/"

df_fdr_c1 <- read.csv(paste0(ruta_directorio,"Data/Datos edgeR/Datos_Bulk_edgeR_Zt_como_Covariable/Separated_C1_FDR_005_Aditive_Model.csv")) %>% 
  dplyr::rename(Gene = genes,
                logFC_C1 = logFC,
                FDR_C1 = FDR)

df_fdr_c9 <- read.csv(paste0(ruta_directorio,"Data/Datos edgeR/Datos_Bulk_edgeR_Zt_como_Covariable/Separated_C9_FDR_005_Aditive_Model.csv")) %>% 
  dplyr::rename(Gene = genes,
                logFC_C9 = logFC,
                FDR_C9 = FDR)

# Common genes between C1 and C9
df_fdr_c1_c9 <- inner_join(df_fdr_c1, df_fdr_c9, by = "Gene")

# UP in both mutants (filtered by logFC)
up_df_fdr_c1_c9  <- df_fdr_c1_c9 %>% 
  dplyr::filter(logFC_C1 > 0 & logFC_C9 > 0) %>% 
  dplyr::select(Gene)

# DOWN in both mutants (filtered by logFC)
down_df_fdr_c1_c9 <- df_fdr_c1_c9 %>% 
  dplyr::filter(logFC_C1 < 0 & logFC_C9 < 0) %>% 
  dplyr::select(Gene)

# Populus tremula x alba HAP2 - Arabidopsis - Arabidopsis GO file 
GOs_2026 <- read.csv('/Users/danielconde/Desktop/WATs/Pathways_Analysis_Mutantes_ZT8_ZT23/Data/Datos Pop x Arab/2026_Poplar_Arab_GOs_TAIR.csv')
GOs_2026  <- GOs_2026  %>% 
  dplyr::select(-Arab) %>% 
  dplyr::rename(Gene = locusName,
                GO = Go_Arab_TAIR)
  
all_genes_with_go <- GOs_2026  

# Arabidopsis gene functional descriptions (TAIR)
arab_gene_description_tair <- read_tsv("/Users/danielconde/Desktop/WATs/Pathways_Analysis_Mutantes_ZT8_ZT23/Data/Mapeo_Phytozome_P.tremHAP2_P.tricho/TAIR10_functional_descriptions.csv")
colnames(arab_gene_description_tair) <- c(
  "genes",
  "Type",
  "Short_description",
  "Curator_summary",
  "Computational_description"
)

arab_gene_description_tair <- arab_gene_description_tair %>% 
  mutate(genes = str_remove(genes, "\\\\..*$"))

# =============================================================================
# 2. PROCESS UP AND DOWN SEPARATELY (SPECIFY AT THE BEGINNING)
# =============================================================================

if (up_or_down  == "up"){
  genes_of_interest <- up_df_fdr_c1_c9$Gene
} else if (up_or_down  == "down"){
  genes_of_interest <- down_df_fdr_c1_c9$Gene  
} else {
  print("Choose up or down")
}

# =============================================================================
# 3. topGO SCRIPT
# =============================================================================

# Define the gene universe as all the genes with GO annotations
geneUniverse <- all_genes_with_go$Gene

# Identify the genes in the list that have assigned GOs with a binary classification
geneList <- factor(as.integer(geneUniverse %in% genes_of_interest))
names(geneList) <- geneUniverse

# topGO object creation 
GOdata <- new("topGOdata",
              ontology = "BP",  # 'BP' for Biological Process, 'MF' for Molecular Function, 'CC' for Cellular Component
              allGenes = geneList,
              annot = annFUN.gene2GO, # We use this because we use our OWN go anotations since poplar does not have an R database like Arab
              gene2GO = split(all_genes_with_go$GO, all_genes_with_go$Gene)) # Mandatory formating when using own annotations

# Definition of test and statistic used
result <- runTest(GOdata, 
                  algorithm = "weight01", 
                  statistic = "fisher")

total_GO_terms <- length(usedGO(GOdata))

# Result table 
GOresult <- GenTable(GOdata, 
                     weight01Fisher = result,
                     orderBy = "weight01Fisher", 
                     ranksOf = "weight01Fisher", 
                     topNodes = total_GO_terms)

GOresult$weight01Fisher <- as.numeric(GOresult$weight01Fisher)

# Sometimes the p-value is NA when it is lower than 1e-30. We fix it with this by fixing the value in NA cases
GOresult$weight01Fisher[is.na(GOresult$weight01Fisher)] <- 1e-30

# FDR adjustment and filtering
GOresult$adj_p_value <- p.adjust(GOresult$weight01Fisher, method = "BH") 

GOresult_filtered <- GOresult[GOresult$adj_p_value < 0.05, ]

# GeneRatio calculation (Significant / Total Annotated in the Gene Universe)
GOresult_filtered$GeneRatio <- as.numeric(GOresult_filtered$Significant) / as.numeric(GOresult_filtered$Annotated)

### DotPlot ###

# Plot the top 20 terms
GOresults_plot <- GOresult_filtered[1:20,]
plot <- ggplot(GOresults_plot, aes(x = GeneRatio, y = reorder(Term, GeneRatio), 
                                   size = GeneRatio, 
                                   color = adj_p_value)) +
  geom_point() +
  scale_color_gradient(limits = c(0, 0.05), low = "purple", high = "orange", 
                       name = "p.adjusted") + 
  scale_size_continuous(name = "GeneRatio", range = c(2, 6)) +
  labs(x = "GeneRatio", y = NULL,
       title = paste0(toupper(up_or_down), " GO Enrichment TOP 50")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),  # Adjust font size for x-axis
    axis.text.y = element_text(size = 7, hjust = 0, vjust = 0.5),  # Adjust font size and alignment for y-axis
    axis.title.x = element_text(size = 14),  # Adjust font size for x-axis label
    axis.title.y = element_text(size = 14),  # Adjust font size for y-axis label
    plot.title = element_text(size = 16, face = "bold"),  # Adjust title size and style
    legend.title = element_text(size = 12),  # Adjust font size for legend title
    legend.text = element_text(size = 10),  # Adjust font size for legend text
    plot.margin = margin(1, 1, 1.5, 1.5, "cm")  # Increase plot margins to prevent clipping
  )

plot

ggsave(plot, 
       filename = paste0("/Users/danielconde/Desktop/WATs/Pathways_Analysis_Mutantes_ZT8_ZT23/Results/Nuevo_Resultados_topGo_Bulk/",
                         up_or_down,
                         "_Top_20_ORA_topGO_Bulk.svg"))

# =============================================================================
# 4. OBTAIN GO - GENE ASSOCIATION LIST
# =============================================================================

significantGO <- function(result_test, alpha = 0.05) {
  if (!inherits(result_test, "topGOresult")) {
    stop("result_test must be a topGOresult object.")
  }
  
  significant_terms <- names(score(result_test)[score(result_test) < alpha])
  
  if (length(significant_terms) == 0) {
    message("No significant GO terms found at alpha = ", alpha)
  }
  
  return(significant_terms)
}

significant_GO_terms <- significantGO(result, alpha = 1) # Keep all, we are only interested in GO - Gene associations
significant_genes_list <- genesInTerm(GOdata, whichGO = significant_GO_terms)

# Format nicely to Gene - GO Term pairs
significant_genes_df <- stack(significant_genes_list)
colnames(significant_genes_df) <- c("Gene", "GO.ID")
significant_genes_df <- significant_genes_df %>% filter(Gene %in% names(geneList[geneList == 1]))

go_collapsed <- significant_genes_df %>%
  group_by(GO.ID) %>%
  summarise(
    genes = list(sort(unique(Gene))),
    n_genes = n_distinct(Gene),
    .groups = "drop"
  )

# Keep only GOs that are significantly enriched
significant_df <- inner_join(go_collapsed, GOresult_filtered, by = "GO.ID")

# Convert the list into a single cell with all genes separated by commas
significant_df <- significant_df %>%
  mutate(genes = sapply(genes, function(x) paste(x, collapse = ", ")))

# Separate into individual rows and remove blank spaces in cells
splited_significant_df <- significant_df %>% 
  separate_rows(genes, sep = ",") 

splited_significant_df$genes <- trimws(splited_significant_df$genes)  

# Join with Arabidopsis functional annotations
splited_significant_df <- inner_join(splited_significant_df, arab_gene_description_tair, by = "genes")

# Save
write_csv(splited_significant_df,
          file = paste0("/Users/danielconde/Desktop/WATs/Pathways_Analysis_Mutantes_ZT8_ZT23/Results/Nuevo_Resultados_topGo_Bulk/",
                        up_or_down,
                        "_All_Significant_GO_Enrichment_With_Arab_Functional_Info.csv"))
