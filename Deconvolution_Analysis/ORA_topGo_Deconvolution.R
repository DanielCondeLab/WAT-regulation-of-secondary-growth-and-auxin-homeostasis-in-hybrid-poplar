########################################################################
# CONDA ENVIRONMENT: Pathways_analysis.yaml
# Over Representation Analysis (ORA) enrichment of UP- and DOWN-regulated DEGs (logFC > 0; logFC < 0; analyzed separately) 
# for every celltype from BayesPrism DEGs filtered by being celltype markers (BayesPrism_One_vs_Rest.R)
# File routes MUST be changed in order to use this code.
# Input files:
# - UPregulated DEGs (FDR < 0.05) per celltype intersected with their celltype markers
# - DOWNregulated DEGs (FDR < 0.05) per celltype intersected with their celltype markers
# - Phytozome P.trichocarpa v14 annotations (to obtain Arabidopsis orthologs)
# - Phytozome Populus tremula x alba HAP2 Database: P. tremula x alba HAP2 gene | Arabidopsis GO terms
# - Arabidopsis gene functional descriptions from TAIR database
# Output files:
# - .csv and .xlsx files with enriched processes (GOs). Excel files are used for combined celltype DotPlot
# - DotPlots of the top 50 most enriched processes (GOs) per celltype
# - .csv files with enriched processes (GOs) and contributing genes
# - .csv files with P. tricho genes, enriched process, Arabidopsis ortholog, and Arabidopsis gene name and TAIR functional info
########################################################################

# =============================================================================
# 0. IMPORT LIBRARIES AND PATHS
# =============================================================================

library(tidyverse)
library(topGO)
library(GO.db)
library(svglite)
library(openxlsx)

# =============================================================================
# 1. IMPORT INPUT FILES
# =============================================================================

analyze_up_or_down <- "DWN"

if (analyze_up_or_down == "UP") {
  data_files <- list.files(
    "/Users/danielconde/Desktop/WATs/BayesnPrism/Results/Expression_Data_Bayes_Seurat_C1yC9/Cruces_C1_C9/UP/Cruzado_con_edgeR/",
    pattern = "\\\\.csv$",
    full.names = TRUE
  )
  plots_output_path <- "/Users/danielconde/Desktop/WATs/BayesnPrism/Results/ORA_topGo_Celltypes/UP/"
} else if (analyze_up_or_down == "DWN"){
  data_files <- list.files(
    "/Users/danielconde/Desktop/WATs/BayesnPrism/Results/Expression_Data_Bayes_Seurat_C1yC9/Cruces_C1_C9/DWN/Cruzado_con_edgeR/",
    pattern = "\\\\.csv$",
    full.names = TRUE
  )
  plots_output_path <- "/Users/danielconde/Desktop/WATs/BayesnPrism/Results/ORA_topGo_Celltypes/DWN/"
}
rm(df)

# Read each CSV and assign it to a variable named after the celltype
for (file in data_files) {
  celltype_name <- tools::file_path_sans_ext(basename(file))
  df <- read.csv(file)
  assign(celltype_name, df)
}

# Phytozome P.trichocarpa v14 annotations (keep P.tricho gene and Arabidopsis ortholog)
phytozome_ptricho_arab <- read_tsv('/Users/danielconde/Desktop/WATs/BayesnPrism/Results/ORA_topGo_Celltypes/Data_Phytozome_P.tricho_&_Arab/P.tricho_4.1_Phytozome/Phytozome/PhytozomeV14/Ptrichocarpa/v4.1/annotation/Ptrichocarpa_533_v4.1.P14.annotation_info.txt')
ptricho_arab <- phytozome_ptricho_arab %>% 
  dplyr::select(locusName,"Best-hit-arabi-name") 

dicc_ptricho_arab <- setNames(ptricho_arab$`Best-hit-arabi-name`, ptricho_arab$locusName)

# Database Arabidopsis GO slim
database_arab_static_go  <- read_tsv('/Users/danielconde/Desktop/WATs/BayesnPrism/Results/ORA_topGo_Celltypes/Data_Phytozome_P.tricho_&_Arab/GOs_Arab/ATH_GO_GOSLIM 2.txt',
                                     skip = 4, col_names = F) 

database_arab_static_go  <- database_arab_static_go %>% 
  dplyr::select(X1, X6) %>% 
  unique()

database_arab_static_go <- database_arab_static_go %>% 
  dplyr::rename(Arabidopsis = 1,
                goid = 2)

# Arabidopsis gene functional descriptions (TAIR)
arab_gene_description_tair <- read_tsv("/Users/danielconde/Desktop/WATs/BayesnPrism/Data/Mapeo_Phytozome_P.tremHAP2_P.tricho/TAIR10_functional_descriptions.csv")
colnames(arab_gene_description_tair) <- c(
  "Arabidopsis",
  "Type",
  "Short_description",
  "Curator_summary",
  "Computational_description"
)

arab_gene_description_tair <- arab_gene_description_tair %>% 
  mutate(Arabidopsis = str_remove(Arabidopsis, "\\\\..*$"))

# Merge dictionary (P.trichocarpa - Arabidopsis) with TAIR annotations
join_ptricho_arab <- ptricho_arab %>% 
  dplyr::rename(Arabidopsis = "Best-hit-arabi-name") 

final_dicc_poplar <- inner_join(arab_gene_description_tair, join_ptricho_arab, by ="Arabidopsis", relationship = "many-to-many") %>% 
  unique()

final_dicc_poplar <- final_dicc_poplar %>% 
  relocate(locusName)

write_csv(final_dicc_poplar,
          file = "/Users/danielconde/Desktop/WATs/BayesnPrism/Data/Mapeo_Phytozome_P.tremHAP2_P.tricho/Ptricho_With_Arab_Functional_Info.csv")

# Automatically select celltypes
celltypes <- tools::file_path_sans_ext(basename(data_files))

# =============================================================================
# 3. PREPARE OBJECTS FOR TOPGO ENRICHMENT 
# =============================================================================

for (ct in celltypes) {
  
  df <- get(paste0(ct), envir = .GlobalEnv)
  genes_of_interest <- df$Gene_Potri
  
  ptricho_arab <- phytozome_ptricho_arab %>% 
    dplyr::select(locusName,"Best-hit-arabi-name") %>% 
    dplyr::rename(Arabidopsis = "Best-hit-arabi-name") 
  
  # Merge to obtain dataset with all P.tricho genes and Arabidopsis GO terms
  all_genes_with_go <- inner_join(ptricho_arab, 
                                  database_arab_static_go, relationship = "many-to-many")
  
  all_genes_with_go <- all_genes_with_go %>% 
    dplyr::select(-Arabidopsis) %>% 
    unique()
  
  # Define the gene universe as all the genes with GO annotations
  geneUniverse <- all_genes_with_go$locusName
 
  # Identify the genes in the list that have assigned GOs with a binary classification
  geneList <- factor(as.integer(geneUniverse %in% genes_of_interest))
  names(geneList) <- geneUniverse
  
  # topGO object creation 
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                annot = annFUN.gene2GO, 
                gene2GO = split(all_genes_with_go$goid, all_genes_with_go$locusName))
  
  # Definition of test and statistic used
  result <- runTest(GOdata, 
                    algorithm = "weight01", 
                    statistic = "fisher")
  
  total_GO_terms <- length(usedGO(GOdata))
  
  # Result table  
  GOresult <- GenTable(GOdata, weight01Fisher = result,
                       orderBy = "weight01Fisher", 
                       ranksOf = "weight01Fisher", 
                       topNodes = total_GO_terms)
  
  GOresult$Term_full <- Term(GOTERM[GOresult$GO.ID])
  GOresult$Term <- GOresult$Term_full
  
  GOresult$weight01Fisher <- as.numeric(GOresult$weight01Fisher)
  
  # Sometimes the p-value is NA when it is lower than 1e-30. We fix it with this by fixing the value in NA cases
  GOresult$weight01Fisher[is.na(GOresult$weight01Fisher)] <- 1e-30
  
  # FDR adjustment and filtering
  GOresult$adj_p_value <- p.adjust(GOresult$weight01Fisher, method = "BH")
  GOresult_filtered <- GOresult[GOresult$adj_p_value < 0.05, ]
  Save_GOresult_filtered <- as.data.frame(GOresult_filtered) %>% unique()
  
  write.csv(Save_GOresult_filtered,
            file = paste0(plots_output_path,"Enrichment/",ct,"_Enriched_topGo.csv"))
  
  write.xlsx(Save_GOresult_filtered,
             file = paste0(plots_output_path,"Excels/",ct,"_Enriched_topGo.xlsx"))
  
  # GeneRatio calculation (Significant / Total Annotated in the Gene Universe)
  GOresult_filtered$GeneRatio <- as.numeric(GOresult_filtered$Significant) / as.numeric(GOresult_filtered$Annotated)
  
  # Function to obtain GO - GENE pairs
  significantGO <- function(result_test, alpha = 0.05) {
    if (!inherits(result_test, "topGOresult")) {
      stop("result_test must be a topGOresult object.")
    }
    significant_terms <- names(score(result_test)[score(result_test) < alpha])
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
    mutate(genes = sapply(genes, function(x) paste(x, collapse = ", "))) %>% 
    unique()
  
  write_csv(significant_df,
            file = paste0(plots_output_path,"Enrichment_with_Genes/",ct,"_With_Genes_Enriched_topGo.csv"))
  
  separated_significant_df <- significant_df %>% 
    separate_rows(genes, sep = ",") %>% 
    mutate(genes = trimws(genes)) %>% 
    dplyr::rename(locusName = genes) %>% 
    dplyr::select(GO.ID,locusName,Term)
  
  # Join with Arabidopsis functional annotations
  final_info_join <- inner_join(separated_significant_df,
                                final_dicc_poplar, 
                                by = "locusName", 
                                relationship = "many-to-many") %>% 
    relocate(locusName,Arabidopsis) %>% 
    unique()
  
  write_csv(final_info_join,
            file = paste0(plots_output_path,"Enrichment_with_Genes_and_Annotations/",ct,"_With_Genes_Enriched_And_Arab_Functional_Info_topGo.csv"))
  
  # Dotplot of top 50 enriched GO terms
  GOresults_plot <- GOresult_filtered[1:50,]
  
  plot <- ggplot(GOresults_plot, aes(x = GeneRatio, y = reorder(Term, GeneRatio), 
                                     size = GeneRatio, 
                                     color = adj_p_value)) +
    geom_point() +
    scale_color_gradient(limits=c(0, 0.05), low = "purple", high = "orange", 
                         name = "p.adjusted") + 
    scale_size_continuous(name = "GeneRatio", range = c(2, 6)) +
    labs(x = "GeneRatio", y = NULL,
         title = paste0("GO Enrichment TOP 50 ",toupper(ct))) +
    theme_minimal()
  
  ggsave(plot, 
         filename = paste0(plots_output_path,"Enrichment_Plots/",ct,".svg")
  )
}
