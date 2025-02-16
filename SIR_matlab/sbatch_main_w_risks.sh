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

paste -d ',' "$clear_csv" "$risk_csv" | while IFS=',' read -r clear_gene risk_gene; do
    # Remove the quotes and extra characters from clear_gene
    clear_gene=${clear_gene:1:-1}
    clear_gene=$(echo "$clear_gene" | sed 's/""/"/g')
    
    # Remove the quotes and extra characters from risk_gene
    risk_gene=${risk_gene:1:-1}
    risk_gene=$(echo "$risk_gene" | sed 's/""/"/g')

    # Submit the job with the paired genes
    sbatch --output="$output_log" --job-name="$clear_gene-SIR" <<EOL
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gabriella.chan@monash.edu

# Load MATLAB module
module load matlab/r2023b

# Run main with each clearance gene
matlab -nodesktop -nodisplay -nosplash -r 'main([$clear_gene], risk=[$risk_gene], null="rewired", out="../../SIR_results/lme_beta/S132_corrs_rewired.csv"); exit;'
EOL
done
