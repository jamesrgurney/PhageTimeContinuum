#!/bin/bash
source /Users/sczerwinski1/miniconda3/etc/profile.d/conda.sh
conda activate defensefinder
df_results="/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/plasmid_reorientation/defensefinder_results"
oriented_dir="/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/plasmid_reorientation/dnaapler_results"
mkdir -p "${df_results}"
parallel -j 10 "
    base=\$(basename \$(dirname {}));
    defense-finder run \
        -o ${df_results}/\${base} \
        -w 10 \
        -A \
        {} > ${df_results}/\${base}.log 2>&1
" ::: $(find "${oriented_dir}" -name "*_reoriented.fasta" ! -name "._*")
