#!/bin/bash

for file in ../data/gene_labels/250915_rerun/*; do
    if [[ "$file" =~ _([0-9]+)\.csv$ ]]; then
        seed="${BASH_REMATCH[1]}"
    fi
    source sbatch_w_risks.sh --input "$file" --output ../results/250910_rewired/ --null rewired --time 6:00:00 --seed "$seed"
done
