#!/bin/bash

outdir="/Volumes/Extreme_SSD/CRISPR_spacer_analyses/r_PAManalysis/plasmid_reorientation/defensefinder_results"

genes_out="${outdir}/ALLplasmid_defense_finder_genes.tsv"
systems_out="${outdir}/ALLplasmid_defense_finder_systems.tsv"

rm -f "$genes_out" "$systems_out"

# combine genes
first=1
for d in "$outdir"/*/; do
    base=$(basename "$d")
    f="${d}/${base}_reoriented_defense_finder_genes.tsv"
    [ -f "$f" ] || continue

    if [ $first -eq 1 ]; then
        awk -v acc="$base" 'BEGIN{FS=OFS="\t"} NR==1{print "accession",$0} NR>1{print acc,$0}' "$f" > "$genes_out"
        first=0
    else
        awk -v acc="$base" 'BEGIN{FS=OFS="\t"} NR>1{print acc,$0}' "$f" >> "$genes_out"
    fi
done

# combine systems
first=1
for d in "$outdir"/*/; do
    base=$(basename "$d")
    f="${d}/${base}_reoriented_defense_finder_systems.tsv"
    [ -f "$f" ] || continue

    if [ $first -eq 1 ]; then
        awk -v acc="$base" 'BEGIN{FS=OFS="\t"} NR==1{print "accession",$0} NR>1{print acc,$0}' "$f" > "$systems_out"
        first=0
    else
        awk -v acc="$base" 'BEGIN{FS=OFS="\t"} NR>1{print acc,$0}' "$f" >> "$systems_out"
    fi
done
