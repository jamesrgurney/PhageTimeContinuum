#!/bin/bash
source /Users/sczerwinski1/miniconda3/etc/profile.d/conda.sh
conda activate defensefinder

WORK_DIR="/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/pa_dereplicated_phages"
DF_RESULTS="${WORK_DIR}/defensefinder_results"
HMM_DIR="${WORK_DIR}/../terLanalysis/hmm_profiles_terL"
OUT_DIR="${WORK_DIR}/terl_hmmsearch_results_oriented"
COMBINED_HMM="${HMM_DIR}/terl_combined.hmm"

mkdir -p "${OUT_DIR}"

cat "${HMM_DIR}"/PF03237.hmm "${HMM_DIR}"/PF03354.hmm \
    "${HMM_DIR}"/PF04466.hmm "${HMM_DIR}"/PF05876.hmm \
    > "${COMBINED_HMM}"
hmmpress -f "${COMBINED_HMM}"

export OUT_DIR COMBINED_HMM

run_one() {
    prt="$1"
    acc=$(basename "$(dirname "${prt}")")
    out_file="${OUT_DIR}/${acc}_terl_hits.tsv"
    log_file="${OUT_DIR}/${acc}.log"

    hmmsearch \
        --tblout "${out_file}" \
        --cpu 1 \
        -E 1e-5 \
        "${COMBINED_HMM}" \
        "${prt}" > "${log_file}" 2>&1 \
    && echo "[DONE] ${acc}" \
    || echo "[FAIL] ${acc}"
}

export -f run_one

echo "Starting hmmsearch on reoriented genomes..."
find "${DF_RESULTS}" -name "*_reoriented.prt" ! -name "._*" ! -name "*.idx" | \
    parallel -j 10 run_one {}

echo "Done. Results in: ${OUT_DIR}"