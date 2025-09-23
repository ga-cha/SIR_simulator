#!/bin/bash

# Copies lines from input csvs to output csvs if the corresponding job has no output file
# i.e. if the job was never run or failed
# If jobs have only failed, cp_repeat_timeout.sh is preferred

in_dir="../data/gene_labels/250910_pairs"
out_dir="../data/gene_labels/250915_rerun"
chk_dir="../results/250910_rewired"
mkdir -p "$out_dir"

count=0

for csv in "$in_dir"/*.csv; do
    csvname=$(basename "$csv")
    outcsv="$out_dir/$csvname"
    > "$outcsv"  # Empty output file
    while IFS= read -r line; do
        # Expect line: GENE1,GENE2,...
        gene1=$(echo "$line" | cut -d',' -f1)
        gene2=$(echo "$line" | cut -d',' -f2)
        # Look for any file in chk_dir containing GENE1_GENE2 in its name
        found=$(find "$chk_dir" -type f -name "*${gene1}_${gene2}*.csv" | head -n1)
        if [ -z "$found" ]; then
            echo "$line" >> "$outcsv"
            ((count+=1))
        fi
    done < "$csv"
done

echo "Copied $count lines to $out_dir"