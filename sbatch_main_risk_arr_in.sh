#!/bin/bash

# Check if an input CSV file is provided as an argument
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 clear_genes.csv risk_genes.csv"
  exit 1
fi

clear_csv="$1"
risk_csv="$2"
output_log="slurm_output.log"

# Check if the input CSV file exists
if [ ! -e "$clear_csv" ] || [ ! -e "$risk_csv" ]; then
  echo "One or more input files not found"
  exit 1
fi

# Read in risk genes into an array for matlab call
# This produces: a, string, array,
while IFS=',' read -r risk_gene; do
    risk_genes+="$risk_gene",
done < "$risk_csv"

# Run main with each clearance gene in a separate job
# with complete risk genes string array
while IFS=',' read -r clear_gene; do
    sbatch --output="$output_log"<<EOL
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=36:00:00
#SBATCH --job-name=random-sir-simulator
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gabriella.chan@monash.edu

# Load MATLAB module
module load matlab/r2023b

# This call should be formatted with two string arrays and two booleans; e.g.
# main("gene1", ["gene1", "gene2", "gene3"], false, false)

risk_genes_formatted=$(echo "$risk_genes" | sed 's/ *, */, /g; s/, *$/,/; s/, /", "/g; s/^/"/; s/.$/"/')
echo $risk_genes_formatted
matlab -nodesktop -nodisplay -nosplash -r "main(\"$clear_gene\", [$risk_genes_formatted], true, false); exit;"

EOL
done < "$clear_csv"
