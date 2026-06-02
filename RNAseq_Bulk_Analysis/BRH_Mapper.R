############################################################
# CONDA ENVIRONMENT: Pathways_analysis.yaml
# This script identifies high-confidence Best Reciprocal Hits (BRHs)
# between the P. tremula × alba HAP2 proteome (Phytozome) and the
# P. alba proteome used by KEGG.
#
# Workflow:
# 1. Import DIAMOND blastp results from both search directions.
# 2. Retain the best hit for each query based on bitscore and e-value.
# 3. Identify reciprocal best-hit pairs.
# 4. Calculate reciprocal alignment coverage.
# 5. Filter RBHs using sequence identity and coverage thresholds.
# 6. Remove isoform suffixes from gene identifiers.
# 7. Export the final 1:1 ortholog mapping table.
#
# Input files:
# - DIAMONDs two-ways results
#  - PtXaAlbH_vs_Palba.tsv
#  - Palba_vs_PtXaAlbH.tsv
#
# Output file:
# - Best Reciprocal Hit
#   - PtXaAlbH_Palba_BRH.csv
############################################################

library(tidyverse)

# Import DIAMOND results for both search directions
AB <- read_tsv("PtXaAlbH_vs_Palba.tsv", col_names = FALSE)
BA <- read_tsv("Palba_vs_PtXaAlbH.tsv", col_names = FALSE)

# Assign informative column names
colnames(AB) <- c(
  "PtXaAlbH",
  "Palba",
  "qlen_AB",
  "slen_AB",
  "pident_AB",
  "length_AB",
  "qcovhsp_AB",
  "scovhsp_AB",
  "evalue_AB",
  "bitscore_AB"
)

colnames(BA) <- c(
  "Palba",
  "PtXaAlbH",
  "qlen_BA",
  "slen_BA",
  "pident_BA",
  "length_BA",
  "qcovhsp_BA",
  "scovhsp_BA",
  "evalue_BA",
  "bitscore_BA"
)

# Retain the best hit for each P. tremula × alba query sequence.
# Hits are ranked primarily by bitscore and secondarily by e-value.
best_AB <- AB %>%
  arrange(PtXaAlbH, desc(bitscore_AB), evalue_AB) %>%
  distinct(PtXaAlbH, .keep_all = TRUE)

# Retain the best hit for each P. alba query sequence.
best_BA <- BA %>%
  arrange(Palba, desc(bitscore_BA), evalue_BA) %>%
  distinct(Palba, .keep_all = TRUE)

# Identify Reciprocal Best Hits (RBHs) present in both searches
rbh <- inner_join(best_AB, best_BA, by = c("PtXaAlbH", "Palba"))

# Compute reciprocal coverage metrics.
# The minimum coverage observed in each direction is used to obtain
# a conservative estimate of alignment coverage.
rbh <- rbh %>%
  mutate(
    coverage_query = pmin(qcovhsp_AB, scovhsp_BA),
    coverage_subject = pmin(scovhsp_AB, qcovhsp_BA)
  )

# Retain only high-confidence ortholog candidates.
# Thresholds:
# - Amino acid identity ≥ 50%
# - Reciprocal query coverage ≥ 70%
# - Reciprocal subject coverage ≥ 70%
rbh_filtrado <- rbh %>%
  filter(
    pident_AB >= 50,
    coverage_query >= 70,
    coverage_subject >= 70
  )

# Remove transcript/isoform suffixes to obtain gene-level identifiers
rbh_final <- rbh_filtrado %>%
  mutate(
    PtXaAlbH = sub("\\.[0-9]+\\.p$", "", PtXaAlbH),
    Palba = sub("\\.[0-9]+$", "", Palba)
  )

# Export the final RBH table
write_csv(
  rbh_final,
  "PtXaAlbH_Palba_BRH.csv"
)