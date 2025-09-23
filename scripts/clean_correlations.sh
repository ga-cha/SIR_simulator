#!/bin/bash

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <path_to_csv_file>"
    exit 1
fi

my_path="$1"
output_path="${my_path%.csv}_clean.csv"

if [ ! -f "$my_path" ]; then
    echo "File not found: $my_path"
    exit 1
fi

# Remove empty lines and lines with missing values
header=$(head -n 1 "$my_path")
echo "$header" > "$output_path"
tail -n +2 "$my_path" | awk -F',' -v header="$header" 'BEGIN {col_count=split(header, cols, ",")} NF==col_count' >> "$output_path"

echo "Cleaned file saved to: $output_path"