#########################################################
# CONDA ENVIRONMENT: Bayes.yaml
# This script performs the deconvolution of bulk data using a scRNA-Seq reference,
# With this, we can obtain fractions of inferred cell types
# and expression data as "pseudocounts" in inferred expression matrices for each celltype
# File routes MUST be changed in order to use this code.
# Input files:
# - RDS of P.trichocarpa scRNA-Seq (DOI: 10.1186/s13059-025-03728-x)
# - (Optional, if already computed) Result file from Seurat's FindAllMarkers of the scRNA-Seq reference
# - Input File (Obtained from Orthologue_Dictionary_Generation.R)
# Output files:
# - (Optional) Result file from Seurat's FindAllMarkers of the scRNA-Seq reference
# - Expression Matrices (pseudocounts) per celltype
# - Barplot (svg) of Cell Fractions for the mean of each genotype (WT, C1, and C9)
# Documentation on deconvolution: https://link.springer.com/article/10.1186/s41231-023-00154-8
########################################################

# =============================================================================
# 0. IMPORT DATA AND LIBRARIES
# =============================================================================

# Import Libraries
library(BayesPrism) # IMPORTANT This package MUST be installed using R: install_github("Danko-Lab/BayesPrism/BayesPrism")
library(Seurat)
library(tidyverse)
library(ggplot2)

### Import scRNA-Seq data
# RDS object of scRNA-Seq
seurat_object <- readRDS("/Users/danielconde/Desktop/WATs/BayesnPrism/Data/Integrated_dataSnRNAseqSTEMFinalclustering.rds")

# Marker Genes from Seurat's FindAllMarkers (if we already have it, we don't need to recompute it all the time)
ruta_resultado_marcadores_tejido_seurat <- "/Users/danielconde/Desktop/WATs/BayesnPrism/Results/FindAllMarkers_Seurat/Filtered_Ptricho_FindAllMarkers_Seurat.csv"
Resultado_FindAllMarkers <- read_csv(ruta_resultado_marcadores_tejido_seurat)

### Mixture File (Bulk Data / Obtained from Orthologue_Dictionary_Generation.R)
bulk_data <- read_tsv("/Users/danielconde/Desktop/WATs/BayesnPrism/Data/Mixture_Files/Tricho_WT_Mutant_BULK_Counts.tsv")

# =============================================================================
# 1. BULK DATASET PREPARATION (bk.dat in the tutorial)
# =============================================================================
# Transform bulk data to the appropriate format (bk.data in the tutorial)

# Samples (Rows) Genes (Columns)
bulk_data <- t(bulk_data)
# Add column names
colnames(bulk_data) <- bulk_data[1,]
# Transform into dataframe
bulk_data <- bulk_data %>%
  as.data.frame()
# Remove row with gene names
bulk_data <-bulk_data[-1, ]
# Transform all values into numeric
bulk_numeric_data <- as.data.frame(
  lapply(bulk_data, function(x) as.numeric(trimws(x)))
)
# Remove zeros
bulk_numeric_data <- bulk_numeric_data %>%
  select(where(~ !all(.x == 0)))
# Assign sample / replicate names
rownames(bulk_numeric_data) <- rownames(bulk_data)

# =============================================================================
# 2. SC DATA: ASSIGNMENT OF CLUSTER NUMBER - CELLTYPE MAPPING
# =============================================================================
# Depending on the granularity level, two vectors can be generated:
# One per large celltype (e.g., Xylem, Phloem...etc)
# Another per specific celltype (e.g., within phloem, distinguish Phloem Sieve Elements, Phloem Companion Cells...etc)
# It is NOT necessary to have both. If you only want to use one, when the other mapping is requested, put NULL

# Vector associating cluster number with celltype name
cell_type_rename_map <- c(
  "0" = "???",
  "1"  = "???",
  "2" = "CO",
  "3" = "MX",
  "4"  = "XMC/JX",
  "5" = "XMC/JX",
  "6"  = "ECyPC",
  "7" = "XP",
  "8" = "CR",
  "9"  = "CO",
  "10" = "CZ",
  "11" = "???",
  "12"  = "PP",
  "13" = "MX",
  "14"  = "SE",
  "15"  = "Phe",
  "16"  = "ECyPC",
  "17"  = "PMC",
  "18" = "CC",
  "19"  = "CR",
  "20" = "???",
  "21" = "CEID",
  "22" = "MX",
  "23" = "ECyPC",
  "24" = "???",
  "25" = "MX"
)

# Add "cell_type" and/or "cell_state" metadata according to the defined mapping
seurat_object@meta.data$cell_type <- cell_type_rename_map[as.character(seurat_object$seurat_clusters)]

# Set 'celltype' as the primary identity
Idents(seurat_object) <- seurat_object@meta.data$cell_type

# Exclusion of celltypes (e.g., SE in my case is a problematic artifact)
seurat_object_filtered <- subset(x = seurat_object, idents = c("SE"), invert = TRUE)

# =============================================================================
# 3. SINGLE CELL REFERENCE PREPARATION (sc.dat in the tutorial)
# =============================================================================
# Select raw counts, transform into a matrix, and transpose it
sc_data <- seurat_object_filtered[["RNA"]]$counts # raw counts (DO NOT use non-linear normalization data)
sc_data <- as.matrix(sc_data)
sc_data <- t(sc_data)

# Vector with cell - cell type association
cell_type_labels <- c(seurat_object_filtered@meta.data[["cell_type"]])
cell_type_labels <- c(seurat_object_filtered@meta.data[["cell_type"]])

# ========================================================================================
# 4. GENERATION OF CELLTYPE MARKERS (SEURAT METHOD)
# ========================================================================================
### SEURAT METHODOLOGY ###
# If the computed file already exists, load it to avoid recomputing, which takes a long time

if (exists("Resultado_FindAllMarkers") &&
    is.data.frame(Resultado_FindAllMarkers) &&
    nrow(Resultado_FindAllMarkers) > 0) {
  seurat_marker_names <- Resultado_FindAllMarkers
} else {
  # Seurat FindAllMarkers with ALL original clusters
  markers_seurat <- FindAllMarkers(seurat_object,
                                   assay = "RNA",
                                   slot = "data",
                                   test.use = "wilcox",
                                   min.pct = 0.01, # Expression in at least 1%
                                   logfc.threshold = 0.1, # Minimum expression
                                   only.pos = T,  # Only positive expression markers
                                   group.by = "cell_type")
  
  # Selection of the best markers (Excluding those associated with discarded cell types)
  seurat_marker_names <- markers_seurat %>%
    filter(p_val_adj <= 0.05,
           avg_log2FC > 0,
           cluster != "SE")
  
  write_csv(seurat_marker_names,
            file = "/Users/danielconde/Desktop/WATs/BayesnPrism/Results/FindAllMarkers_Seurat/Filtered_Ptricho_FindAllMarkers_Seurat.csv")
  
}
# Filter single cell data to only keep the matrix with selected markers
seurat_sc_data <- sc_data[, colnames(sc_data) %in% seurat_marker_names$gene]
rm(seurat_object)

# ========================================================================================
# 5. DECONVOLUTION WITH BAYESPRISM
# ========================================================================================

# Create the "prism" object
myPrism <- new.prism(
  reference=seurat_sc_data,            # Matrix marker genes - cells (the barcodes)
  mixture=bulk_numeric_data,           # Matrix genes - replicates (WT and MUTANTS)
  input.type="count.matrix",           # DO NOT modify
  cell.type.labels = cell_type_labels, # cell_type_labels 
  cell.state.labels = NULL,            # NULL since we only have defined one mapping vector 
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

# Run BayesPrism
bp.res <- run.prism(prism = myPrism, n.cores=8)
bp.res # View results

# ========================================================================================
# 6. CELL FRACTION (MEAN OF WT | C1 | C9)
# ========================================================================================

### CELLTYPE FRACTION CALCULATED IN BULK ###
# Extraction of inferred cell fractions
theta <- get.fraction (bp=bp.res,
                       which.theta="final",
                       state.or.type="type")

df <- as.data.frame(theta)

# Save rownames to a 'sample' column
df$sample <- rownames(df)

# Extract condition from the rowname (WT, C1, or C9)
df <- df %>%
  mutate(condition = str_extract(sample, "WT|C1|C9"))

# Calculate the mean per condition for all cell types
df_mean <- df %>%
  group_by(condition) %>%
  summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE))
  )

# Pivot from wide to long format (for ggplot)
df_long <- df_mean%>%
  pivot_longer(
    cols = -condition,           # all except 'sample'
    names_to  = "celltype",        # name of the new column with the celltype
    values_to = "fraction"       # numeric value
  )
# Add percentage for representation
df_long <- df_long %>% mutate(percent = fraction * 100)

# Color palette
my_dark_palette_15 <- c(
  "#4E79A7", # muted dark blue
  "#A0CBE8", # grayish blue
  "#59A14F", # dark green
  "#8CD17D", # muted green
  "#B6992D", # dark gold
  "#F28E2B", # burnt orange
  "#E15759", # muted red
  "#B07AA1", # dark purple
  "#9D7660", # grayish brown
  "#79706E", # dark gray
  "#6B4C9A", # deep violet
  "#276A7F", # dark teal
  "#5A3E36", # chocolate brown
  "#566573", # dark bluish gray
  "#7D3C98"  # muted intense purple
)


# Stacked barplot: one bar per condition with mean percentages
fraction_plot <- ggplot(df_long, aes(x = condition, y = percent, fill = celltype)) +
  geom_bar(stat = "identity") +         # use values as is
  scale_fill_manual(values = my_dark_palette_15) +
  labs(x = "Sample", y = "% of celltype" , fill = "Celltype") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

ggsave(fraction_plot,
       create.dir = T,
       width =12,
       heigh =10,
       units = "cm",
       file = paste0("/Users/danielconde/Desktop/WATs/BayesnPrism/Results/Fractions/WATs_Deconvolution_Fractions.svg"))

# ========================================================================================
# 7. EXPRESSION MATRICES PER CELLTYPE (PSEUDOCOUNTS; Z EXPRESSION MATRIX)
# ========================================================================================
# GENERATE A Z MATRIX FILE PER CELLTYPE
celltypes_labels <- unique(cell_type_labels)

for (celltype in celltypes_labels) {
  
  Z <- get.exp(bp = bp.res,
               state.or.type = "type",
               cell.name = celltype)
  
  Z_2 <- as.data.frame(Z)
  
  # Replace / with _ in the filename
  safe_celltype <- gsub("/", "_", celltype)
  
  write.csv(
    Z_2,
    file = paste0(
      "/Users/danielconde/Desktop/WATs/BayesnPrism/Results/Z_Matrix_Pseudocounts_Deconvoluted_Celltypes/",
      safe_celltype,
      ".csv"
    )
  )
}
