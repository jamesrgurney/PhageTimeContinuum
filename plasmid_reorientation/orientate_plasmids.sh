#!/bin/bash
source /Users/sczerwinski1/miniconda3/etc/profile.d/conda.sh
conda activate dnaapler

GENOMES_DIR="/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/25Mar26/data/individual_plasmid_fasta"
ORIENTED_DIR="/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/plasmid_reorientation/dnaapler_results"

mkdir -p "${ORIENTED_DIR}"


# Orient each genome with dnaapler

export ORIENTED_DIR

run_dnaapler() {
    fna="$1"
    base=$(basename "${fna}" .fna)
    out="${ORIENTED_DIR}/${base}"

    dnaapler plasmid \
        -i "${fna}" \
        -o "${out}" \
        -p "${base}" \
        -t 8 \
        -f \
        --autocomplete mystery \
        > "${ORIENTED_DIR}/${base}.log" 2>&1 \
    && echo "[DONE] ${base}" \
    || echo "[FAIL] ${base}"
}

export -f run_dnaapler

echo "Starting dnaapler across $(ls ${GENOMES_DIR}/*.fna | wc -l) genomes..."

find "${GENOMES_DIR}" -name "*.fna" | \
    parallel -j 10 run_dnaapler {}

echo "Done. Oriented genomes in: ${ORIENTED_DIR}"
