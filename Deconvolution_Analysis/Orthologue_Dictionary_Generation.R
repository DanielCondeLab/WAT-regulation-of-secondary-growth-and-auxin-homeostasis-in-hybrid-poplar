###################################################################
# CONDA ENVIRONMENT: Bayes.yaml
# This script generates the input file from bulk data for deconvolution by orthologue mapping
# In addition, it generates the P.tremula x alba HAP2 - P.tricho dictionary based on the established criteria
# using the Phytozome Biomart file as the starting point. Rules for each relationships (1:N, N:N, N:1, 1:1) are explained below.
# File routes MUST be changed in order to use this code.
# Input files:
# - Raw bulk counts for the WATs
# - DEGs for C1 and C9 calculated separately with the additive model using edgeR
# - P.trem x alba - P.tricho ortholog file from Phytozome Biomart 
# Output files:
# - P.tremula x alba - P.trichocarpa dictionary
# - Deconvolution Input File
##############################################################

# =============================================================================
# 0. IMPORT LIBRARIES
# =============================================================================
library(tidyverse)

# =============================================================================
# 1. LOAD RAW COUNTS AND POPLAR DICTIONARY
# =============================================================================

### General path for WAT data
data_path <- "/Users/danielconde/Desktop/WATs/BayesnPrism/Data/WAT_Data/"

### Phytozome ortholog file (P.trem x alba HAP2 - P.tricho)
Data <- read_tsv("/Users/danielconde/Desktop/WATs/BayesnPrism/Data/Mapeo_Phytozome_P.tremHAP2_P.tricho/PtaHAP2vsTrico.txt")

bulk_counts_30_replicates <- read_csv("/Users/danielconde/Desktop/WATs/BayesnPrism/Data/WAT_Data/Counts_30samples.csv") %>% 
  dplyr::rename(Gene_Symbol = Gene)

# Differentially expressed genes for dictionary creation
c9 <- read_csv(paste0(data_path, "Separated_C9_FDR_005_Aditive_Model.csv"))
c9 <- c9 %>% 
  dplyr::rename(Gene = genes,
                logFC_C9  = logFC,
                logCPM_C9 = logCPM,
                PValue_C9 = PValue,
                FDR_C9 = FDR) %>% 
  dplyr::select(-"F")


c1 <- read_csv(paste0(data_path, "Separated_C1_FDR_005_Aditive_Model.csv"))

c1 <- c1 %>% 
  dplyr::rename(Gene = genes,
                logFC_C1  = logFC,
                logCPM_C1 = logCPM,
                PValue_C1 = PValue,
                FDR_C1 = FDR) %>% 
  dplyr::select(-"F")

output_path <- "/Users/danielconde/Desktop/WATs/BayesnPrism/Data/Mixture_Files/"

# ========================================================================================
# 2. CREATE DICTIONARY TO CONVERT P.TREM X ALBA HAP2 TO P.TRICHOCARPA ACCORDING TO RULES
# ========================================================================================

# Convert P.trem x alba genes to P.trichocarpa and remove genes without orthologs
# 1. 1:1 ortholog mapping
one_to_one_mapping <- Data %>% 
  dplyr::select(`Gene Name`, `Ortholog Gene Name`, Relationship) %>% 
  rename(Alba = `Gene Name`, Trichocarpa = `Ortholog Gene Name`) %>% 
  filter(Relationship == "one-to-one") %>% 
  distinct()

# 2. 1:N ortholog mapping
one_to_many_mapping <- Data %>% 
  dplyr::select(`Gene Name`, `Ortholog Gene Name`, Relationship) %>% 
  rename(Alba = `Gene Name`, Trichocarpa = `Ortholog Gene Name`) %>% 
  filter(Relationship == "one-to-many") %>% 
  distinct()

# N:N ortholog dictionary (sort alphabetically and keep the first occurrence)
many_to_many_mapping <- Data %>% 
  select(`Gene Name`, `Ortholog Gene Name`, Relationship) %>% 
  rename(Alba = `Gene Name`, Trichocarpa = `Ortholog Gene Name`) %>% 
  filter(Relationship == "many-to-many") %>% 
  arrange(Alba, Trichocarpa) %>% 
  group_by(Trichocarpa) %>% 
  slice(1) %>%        # keep the first ortholog for each Alba
  ungroup()  

# N:1 ortholog dictionary (merge and keep only those that are DEGs)
many_to_one_mapping <- Data %>% 
  dplyr::select(`Gene Name`, `Ortholog Gene Name`, Relationship) %>% 
  rename(Alba = `Gene Name`, Trichocarpa = `Ortholog Gene Name`) %>% 
  filter(Relationship == "many-to-one") 

c1_c9 <- inner_join(c1, c9, by = "Gene") %>% 
  rename("Alba" = Gene)

c1_c9_many_to_one_mapping <- inner_join(many_to_one_mapping, c1_c9, by = "Alba") %>% 
  dplyr::select("Alba", "Trichocarpa", "Relationship") %>% 
  arrange(Alba, Trichocarpa) %>% 
  distinct()

# Remove rows where chromosomes do not match
c1_c9_many_to_one_mapping <- c1_c9_many_to_one_mapping %>%
  mutate(
    num_Alba   = str_extract(Alba, "(?<=\\\\.)\\\\d+(?=G)"),
    num_Tricho = str_extract(Trichocarpa, "(?<=\\\\.)\\\\d+(?=G)")
  ) %>%
  # convert to integer to ignore leading zeros
  mutate(
    num_Alba   = as.integer(num_Alba),
    num_Tricho = as.integer(num_Tricho)
  ) %>%
  # FILTER: keep only rows where the numbers match
  filter(num_Alba == num_Tricho) %>%
  select(-num_Alba, -num_Tricho) %>% 
  rename(Gene_Symbol = Alba)

# Join with WT counts and keep the row with the highest total across WT rows
WT_Counts <- bulk_counts_30_replicates %>% 
  dplyr::select(dplyr::contains("WT"))
WT_Counts$Gene_Symbol <- bulk_counts_30_replicates$Gene_Symbol

many_to_one_c1_c9_count <- inner_join(WT_Counts, c1_c9_many_to_one_mapping, by = "Gene_Symbol") %>% 
  mutate(TOTAL = rowSums(across(where(is.numeric)))) %>% 
  dplyr::select(Gene_Symbol, Trichocarpa, TOTAL) %>% 
  group_by(Trichocarpa) %>% 
  slice_max(TOTAL, n = 1, with_ties = FALSE) %>% 
  ungroup() %>% 
  dplyr::select(-TOTAL)
  dplyr::rename(Alba = Gene_Symbol)

# Merge the 4 dictionary datasets and create final dictionary
combined_dictionary <- dplyr::bind_rows(one_to_one_mapping, 
                                        one_to_many_mapping, 
                                        many_to_many_mapping, 
                                        many_to_one_c1_c9_count) %>% 
  unique()

# Save P.tremula x alba - P.trichocarpa dictionary
write_csv(combined_dictionary,
          file = paste0("/Users/danielconde/Desktop/WATs/BayesnPrism/Data/Mapeo_Phytozome_P.tremHAP2_P.tricho/PtremxalbaHAP2_to_Ptricho.csv"))

# Variable containing the dictionary
dict_poplars <- setNames(combined_dictionary$Trichocarpa, combined_dictionary$Alba) 

# ====================================================================================
# 3. CONVERSION INTO P.TRICHO AND SAVE DECONVOLUTION INPUT FILE
# ==================================================================================== 

# Mixture file WT AND MUTANTS
Tricho_WT_Mutant_BULK_Counts <- bulk_counts_30_replicates %>%
  dplyr::mutate(
    Trichocarpa_ID = dict_poplars[Gene_Symbol]) %>% 
  dplyr::select(-Gene_Symbol) %>% 
  dplyr::rename(Gene_Symbol = Trichocarpa_ID) %>% 
  relocate(Gene_Symbol) %>% 
  drop_na()

write_tsv(
  Tricho_WT_Mutant_BULK_Counts,
  file = paste0(output_path, "Tricho_WT_Mutant_BULK_Counts.tsv")
)
