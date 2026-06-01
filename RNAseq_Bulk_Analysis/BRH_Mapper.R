############################################################
# CONDA ENVIRONMENT: Pathways_analysis.yaml
# This script keeps the best hit from Best Reciprocal Hit (BRH) analysis
# between P.tremula x alba HAP2 from Phytozome and P.alba from KEGG
# Input files:
# - BRH results (both ways):
#   - PtXaAlbH_vs_Palba.tsv
#   - Palba_vs_PtXaAlbH.tsv 
# Output files:
# - Consensus 1:1 BRH 
#   - PtXaAlbH_Palba_RBH.csv
############################################################

# Load BRH Data
AB <- read_tsv("/Users/danielconde/Library/CloudStorage/GoogleDrive-lab171@intranet.cbgp.upm.es/Mi unidad/Drive_Juan/Colaboracion_Agustina_2026/Data/PtXaAlbH_vs_Palba.tsv", col_names = FALSE)
BA <- read_tsv("/Users/danielconde/Library/CloudStorage/GoogleDrive-lab171@intranet.cbgp.upm.es/Mi unidad/Drive_Juan/Colaboracion_Agustina_2026/Data/Palba_vs_PtXaAlbH.tsv", col_names = FALSE)

colnames(AB) <- c("PtXaAlbH", "Palba", "pident_AB", "length_AB", "evalue_AB", "bitscore_AB")
colnames(BA) <- c("Palba", "PtXaAlbH", "pident_BA", "length_BA", "evalue_BA", "bitscore_BA")

# Arrange BRH Data
best_AB <- AB %>%
  arrange(PtXaAlbH, desc(bitscore_AB), evalue_AB) %>%
  distinct(PtXaAlbH, .keep_all = TRUE)

best_BA <- BA %>%
  arrange(Palba, desc(bitscore_BA), evalue_BA) %>%
  distinct(Palba, .keep_all = TRUE)

# Consensus
rbh <- inner_join(best_AB, best_BA, by = c("PtXaAlbH", "Palba"))

rbh <- rbh %>%
  mutate(
    PtXaAlbH = sub("\\.[0-9]+\\.p$", "", PtXaAlbH),
    Palba = sub("\\.[0-9]+$", "", Palba)
)
# Save
write_csv(rbh, "/Users/danielconde/Library/CloudStorage/GoogleDrive-lab171@intranet.cbgp.upm.es/Mi unidad/Drive_Juan/Colaboracion_Agustina_2026/Results/PtXaAlbH_Palba_RBH.csv")
