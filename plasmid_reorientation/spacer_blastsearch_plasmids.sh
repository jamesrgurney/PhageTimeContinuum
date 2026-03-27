#!/bin/bash
source /Users/sczerwinski1/miniconda3/etc/profile.d/conda.sh
conda activate bioinf-base

ORIENTED_DIR="/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/plasmid_reorientation/dnaapler_results"
SPACERS="/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/pa_dereplicated_phages/data/spacers.gt50.fna"
SPACERS2="/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/pa_dereplicated_phages/data/spacers.le50.fna"
OUT_DIR="/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/plasmid_reorientation/blast_results"
DB_DIR="${OUT_DIR}/blastdb"

mkdir -p "${OUT_DIR}" "${DB_DIR}"

# Concatenate and shorten headers to just the .part_ACCESSION
find "${ORIENTED_DIR}" -name "*_reoriented.fasta" ! -name "._*" | \
    while read f; do
        sed 's/^>\(.*\)\.part_\([^ ]*\).*/>\2/' "$f"
    done > "${DB_DIR}/plasmids_reoriented.fna"

echo "Concatenated FASTA built."

# Step 2: Make BLAST database WITHOUT -parse_seqids
makeblastdb \
    -in "${DB_DIR}/plasmids_reoriented.fna" \
    -dbtype nucl \
    -out "${DB_DIR}/orient_plasmidsSet"

echo "BLAST database built."

# Step 3: Merge spacers
cat "${SPACERS}" "${SPACERS2}" > "${OUT_DIR}/spacers_merged.fna"

# Step 4: Run BLASTn
blastn \
    -query "${OUT_DIR}/spacers_merged.fna" \
    -db "${DB_DIR}/orient_plasmidsSet" \
    -outfmt "6 qseqid sseqid sstart send sstrand pident length mismatch gapopen qlen" \
    -perc_identity 90 \
    -num_threads 8 \
    -out "${OUT_DIR}/spacers_vs_orient_plasmids.blast6.tsv"

echo "Done. Results in: ${OUT_DIR}"
