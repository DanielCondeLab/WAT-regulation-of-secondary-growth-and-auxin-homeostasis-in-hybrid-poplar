############################################################
# CONDA ENVIRONMENT: Bayes.yaml
# edgeR one-vs-rest by celltype using deconvoluted profiles to obtain marker genes / highly celltype-specific
# We look for marker genes for each celltype by comparing:
# target celltype expresion vs rest of celltype mean expression
# We ONLY use WT replicates because we cannot trust mutant effects
# We apply a Design: ~ 0 + celltype (the 0 removes the intercept, it serves to not base the results on any specific celltype)
# File routes MUST be changed in order to use this code.
# Input files:
# - Expression Matrices (pseudocounts) per celltype from BayesPrism_DEGs.R
# Output files:
# - .csv file with all marker genes (FDR < 0.05) from all celltypes.
############################################################

# ==============================================================================
# 0. IMPORT LIBRARIES
# ==============================================================================
library(edgeR)
library(tidyverse)

# ==============================================================================
# 1. IMPORT EXPRESSION MATRIX DATA PER CELLTYPE (PSEUDOCOUNTS)
# ==============================================================================

ruta_dir <- "/Users/danielconde/Desktop/WATs/BayesnPrism/Results/Z_Matrix_Pseudocounts_Deconvoluted_Celltypes/"
ruta_salida <- "/Users/danielconde/Desktop/WATs/BayesnPrism/Results/Z_Matrix_Pseudocounts_Deconvoluted_Celltypes/Resultados_Marcadores_One_vs_Rest/"

archivos <- list.files(
  path = ruta_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

# ==============================================================================
# 2. AUXILIARY FUNCTIONS FOR DATA PARSING
# ==============================================================================

# Extract celltype name from filename
obtener_tejido <- function(x) {
  basename(x) %>%
    str_remove("\\.csv$")
}

# Parse only WT samples
# Expected examples:
# WT1ZT23_Read_Count | WT2ZT8_Read_Count

# Anything other than WT returns NULL
parsear_muestra_wt <- function(x) {
  
  x_limpio <- str_remove(x, "_Read_Count$")
  patron_wt <- "^WT(\\d+)ZT(\\d+)$"
  
  if (!str_detect(x_limpio, patron_wt)) {
    return(NULL)
  }
  
  m <- str_match(x_limpio, patron_wt)
  
  rep <- m[, 2]
  zt  <- paste0("ZT", m[, 3])
  
  tibble(
    sample_name = x,
    genotype = "WT",
    rep = rep,
    zt = zt
  )
}

# ==============================================================================
# 3. READING FILES AND CREATING MATRICES AND METADATA
# ==============================================================================
# Each matrix will be converted to:
# genes x samples
#
# Columns will be named as:
# sampleName_celltype
#
# Example:
# WT1ZT23_Read_Count_Phe

lista_matrices <- list()
lista_metadata <- list()

for (archivo in archivos) {
  
  tejido <- obtener_tejido(archivo)
  message("Reading celltype: ", tejido)
  
  df <- read_csv(archivo, show_col_types = FALSE)
  
  # Rename the first column to sample_name
  colnames(df)[1] <- "sample_name"
  
  # Parse only WT
  meta_lista <- lapply(df$sample_name, parsear_muestra_wt)
  idx_wt <- which(!sapply(meta_lista, is.null))
  
  # If no WT, skip file
  if (length(idx_wt) == 0) {
    warning(paste("No WT samples found in", tejido, "- skipping this celltype"))
    next
  }
  
  # WT Metadata
  meta_muestras <- bind_rows(meta_lista[idx_wt])
  
  # Filter only WT rows
  df <- df[idx_wt, , drop = FALSE]
  
  # Check for duplicates
  if (anyDuplicated(meta_muestras$sample_name) > 0) {
    stop(paste("Duplicate sample_names found in", tejido))
  }
  
  # Numeric matrix: rows=samples, columns=genes
  mat_samples_x_genes <- df %>%
    select(-sample_name) %>%
    as.data.frame()
  
  rownames(mat_samples_x_genes) <- meta_muestras$sample_name
  
  # Convert to genes x samples
  mat_genes_x_samples <- t(as.matrix(mat_samples_x_genes))
  storage.mode(mat_genes_x_samples) <- "numeric"
  
  # Create unique column names
  nuevas_cols <- paste0(meta_muestras$sample_name, "__", tejido)
  colnames(mat_genes_x_samples) <- nuevas_cols
  
  # Metadata per column
  meta_cols <- meta_muestras %>%
    mutate(
      celltype = tejido,
      column_name = nuevas_cols
    )
  
  lista_matrices[[tejido]] <- mat_genes_x_samples
  lista_metadata[[tejido]] <- meta_cols
}

if (length(lista_matrices) == 0) {
  stop("No valid celltypes remained after filtering for WT.")
}


# ==============================================================================
# 4. COMBINE MATRICES AND METADATA (edgeR INPUT)
# ==============================================================================
counts <- do.call(cbind, lista_matrices)
meta <- bind_rows(lista_metadata)

# Reorder metadata to match counts
meta <- meta %>%
  slice(match(colnames(counts), column_name))

# ==============================================================================
# 5. SELECT FACTOR FOR DESIGN MATRIX AND CREATE DGELIST
# ==============================================================================
meta <- meta %>%
  mutate(
    celltype = factor(celltype)
  )

y <- DGEList(counts = counts,
             genes = rownames(counts),
             group = meta$celltype)

# ==============================================================================
# 6. FILTERING BY EXPRESSION AND NORMALIZATION (TMM)
# ==============================================================================
keep <- filterByExpr(y)

y <- y[keep, , keep.lib.sizes = FALSE]

y <- normLibSizes(y)
plotMDS(y, col= rep(1:13, each = 10))

# ==============================================================================
# 7. EXPERIMENTAL DESIGN (CELLTYPE IS WHAT WE COMPARE)
# ==============================================================================
design <- model.matrix(~ 0 + celltype, data = meta)

print(colnames(design))

# ==============================================================================
# 8. DISPERSION ESTIMATION AND FITTING
# ==============================================================================
y <- estimateDisp(y, design, robust = TRUE)

y$common.dispersion
plotBCV(y)

fit <- glmQLFit(y, design, robust = TRUE)
plotQLDisp(fit)

# ==============================================================================
# 9. CONSTRUCTION OF ONE-VS-REST CONTRASTS
# ==============================================================================
niveles_tejido <- levels(meta$celltype)

matriz_contrastes <- sapply(niveles_tejido, function(tt) {
  
  v <- rep(0, ncol(design))
  names(v) <- colnames(design)
  
  otros <- setdiff(niveles_tejido, tt)
  
  v[paste0("celltype", tt)] <- 1
  v[paste0("celltype", otros)] <- -1 / length(otros)
  
  v
})

colnames(matriz_contrastes) <- paste0(niveles_tejido, "_vs_rest")

heatmap(matriz_contrastes)

# ==============================================================================
# 10. APPLY ONE-VS-REST FOR EACH CELLTYPE
# ==============================================================================
resultados <- list()

for (i in seq_len(ncol(matriz_contrastes))) {
  
  nombre_contraste <- colnames(matriz_contrastes)[i]
  tejido_objetivo <- str_remove(nombre_contraste, "_vs_rest")
  
  message("Calculating contrast: ", nombre_contraste)
  
  qlf <- glmQLFTest(fit, contrast = matriz_contrastes[, i])
  
  tab <- topTags(qlf, n = Inf, sort.by = "PValue")$table %>%
    rownames_to_column("gene") %>%
    mutate(
      celltype = tejido_objetivo,
      contrast = nombre_contraste,
      marker_direction = ifelse(logFC > 0, "Up_in_target_celltype", "Down_in_target_celltype")
    )
  
  # Save full table
  write_csv(
    tab,
    file.path(ruta_salida, paste0(nombre_contraste, "_edgeR_full.csv"))
  )
  
  # Save positive markers for the celltype
  tab_markers <- tab %>%
    filter(logFC > 0, FDR < 0.05)
  
  write_csv(
    tab_markers,
    file.path(ruta_salida, paste0(nombre_contraste, "_markers_FDR_005.csv"))
  )
  
  resultados[[nombre_contraste]] <- tab
}

# ==============================================================================
# 11. COMBINE ALL MARKERS FROM ALL CELLTYPES INTO A SINGLE FILE
# ==============================================================================

resumen_general <- bind_rows(
  lapply(names(resultados), function(nombre_contraste) {
    resultados[[nombre_contraste]] %>%
      filter(logFC > 0, FDR < 0.05)
  })
)

write_csv(
  resumen_general,
  file.path(ruta_salida, "Resumen_General_Markers_FDR_005.csv")
)

resumen_num_genes <- resumen_general %>%
  group_by(celltype) %>%
  summarise(
    n_markers = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(n_markers))

write_csv(
  resumen_num_genes,
  file.path(ruta_salida, "Resumen_Numero_Markers_por_Tejido.csv")
)

# ==============================================================================
# 12. PERCENTAGE OF OVERLAP BETWEEN SEURAT MARKERS AND THESE
# ==============================================================================
# Load all markers already generated
archivos_markers <- list.files(
  ruta_salida,
  pattern = "_markers_FDR_005\\.csv$",
  full.names = TRUE
)

listas_markers <- lapply(archivos_markers, function(f) {
  df <- readr::read_csv(f, show_col_types = FALSE)
  unique(df$gene)
})

names(listas_markers) <- basename(archivos_markers) |>
  stringr::str_remove("_markers_FDR_005\\.csv$") |>
  stringr::str_remove("_vs_rest")

# Absolute overlap matrix
# Percentage relative to the row list
solapamiento_pct <- outer(
  names(listas_markers),
  names(listas_markers),
  Vectorize(function(a, b) {
    100 * length(intersect(listas_markers[[a]], listas_markers[[b]])) / length(listas_markers[[a]])
  })
)

rownames(solapamiento_pct) <- names(listas_markers)
colnames(solapamiento_pct) <- names(listas_markers)

solapamiento_pct %>%
  as.data.frame() %>%
  rownames_to_column("celltype1") %>%
  pivot_longer(-celltype1, names_to="celltype2", values_to="pct") %>%
  ggplot(aes(celltype1, celltype2, fill = pct)) +
  geom_tile() +
  scale_fill_viridis_c(option = "B", direction = -1) +
  coord_equal() +
  theme_minimal() +
  labs(fill = "% overlap")
