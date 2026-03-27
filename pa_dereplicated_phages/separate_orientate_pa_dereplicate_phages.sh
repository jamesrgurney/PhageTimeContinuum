#!/bin/bash
source /Users/sczerwinski1/miniconda3/etc/profile.d/conda.sh
conda activate dnaapler

INPUT_FA="/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/pa_dereplicated_phages/data/pseudomonas_dereplicated.fa"
SPLIT_DIR="/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/pa_dereplicated_phages/per_genome"
ORIENTED_DIR="/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/pa_dereplicated_phages/oriented"

mkdir -p "${SPLIT_DIR}" "${ORIENTED_DIR}"

# ============================================================
# Step 1: Split multi-FASTA into individual files
# ============================================================
echo "Splitting multi-FASTA..."

awk '/^>/ {
    if (outfile) close(outfile)
    # Extract accession = first word after >
    match($0, /^>([^ ]+)/, arr)
    acc = arr[1]
    outfile = "'"${SPLIT_DIR}"'/" acc ".fna"
    print > outfile
    next
}
{ print > outfile }' "${INPUT_FA}"

n_split=$(ls "${SPLIT_DIR}"/*.fna | grep -v "/\._" | wc -l)
echo "  Split into ${n_split} individual FASTA files"

# ============================================================
# Step 2: Orient each genome with dnaapler
# ============================================================
export ORIENTED_DIR

run_dnaapler() {
    fna="$1"
    base=$(basename "${fna}" .fna)
    out="${ORIENTED_DIR}/${base}"

    dnaapler phage \
        -i "${fna}" \
        -o "${out}" \
        -p "${base}" \
        -t 1 \
        -f \
        --autocomplete mystery \
        > "${ORIENTED_DIR}/${base}.log" 2>&1 \
    && echo "[DONE] ${base}" \
    || echo "[FAIL] ${base}"
}

export -f run_dnaapler

echo "Starting dnaapler on ${n_split} genomes..."
find "${SPLIT_DIR}" -name "*.fna" ! -name "._*" | \
    parallel -j 10 run_dnaapler {}

echo "Done. Oriented genomes in: ${ORIENTED_DIR}"
