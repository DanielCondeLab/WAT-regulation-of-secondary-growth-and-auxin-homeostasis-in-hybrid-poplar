################################################
# CONDA ENVIRONMENT: orthofinder.yaml
# This script performs Best Reciprocal Hit Analysis for Proteomes of P.tremula x alba HAP2
# from Phytozome and P.alba from KEGG
# MUST be concatenated with BRH_Mapper.R script
# Input files:
# - P.tremula x alba HAP2 principal transcripts from Phytozome
# - Populus Alba Proteome from KEGG
# Output files:
# - Populus Dicctionary Mapping (P.tremula x alba HAP2 - P.Alba from KEGG)
################################################

#!/bin/bash
set -euo pipefail

source /Users/danielconde/miniforge3/etc/profile.d/conda.sh
conda activate orthofinder

diamond makedb \
  --in PtremulaxPopulusalbaHAP2_716_v5.1.protein_primaryTranscriptOnly.fa \
  -d PtXaAlbH

diamond makedb \
  --in GCF_005239225.2_ASM523922v2_protein.faa \
  -d Palba_KEGG

diamond blastp \
  -q PtremulaxPopulusalbaHAP2_716_v5.1.protein_primaryTranscriptOnly.fa \
  -d Palba_KEGG.dmnd \
  -o PtXaAlbH_vs_Palba.tsv \
  --ultra-sensitive \
  --evalue 1e-10 \
  --max-target-seqs 25 \
  --max-hsps 1 \
  --outfmt 6 qseqid sseqid qlen slen pident length qcovhsp scovhsp evalue bitscore


diamond blastp \
  -q GCF_005239225.2_ASM523922v2_protein.faa \
  -d PtXaAlbH.dmnd \
  -o Palba_vs_PtXaAlbH.tsv \
  --ultra-sensitive \
  --evalue 1e-10 \
  --max-target-seqs 25 \
  --max-hsps 1 \
  --outfmt 6 qseqid sseqid qlen slen pident length qcovhsp scovhsp evalue bitscore


