#!/bin/bash

# Copies lines from input csvs to output csvs if the corresponding job timed out

in_dir="../data/gene_labels/250908_pairs"
out_dir="../data/gene_labels/250908_rerun"
mkdir -p "$out_dir"

# Get failed job names (comma-separated) MMDDYY
failed_jobs=$(show_job --start_time 090725 | grep "TIMEOUT" | awk -F'|' '{print $3}' | sed 's/[[:space:]]*//g')

count=0
# For each failed job, grep in the csv and append to output
while IFS= read -r line; do
    # Skip lines not matching the pattern CLEAR_RISK1_SEED
    if ! [[ "$line" =~ ^[^_]+_[^_]+_[^_]+$ ]]; then
        continue
    fi
    job=$(echo "$line" | awk -F'_' '{print $1","$2}')
    seed=$(echo "$line" | awk -F'_' '{print $3}')
    # echo "Processing job: $job; seed: $seed"
    for csv in "$in_dir"/*_"$seed".csv; do
        # echo "Looking for: $csv"
        [ -e "$csv" ] || continue  # skip if no file matches
        csvname=$(basename "$csv")
        grep "$job" "$csv" >> "$out_dir/$csvname"
		((count+=1))
    done
done <<< "$failed_jobs"
echo "Successfully copied $count lines"
