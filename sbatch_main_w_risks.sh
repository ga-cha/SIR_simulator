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

while IFS=',' read -r clear_gene; do
    clear_gene=${clear_gene:1:-1}
    clear_gene=$(echo "$clear_gene" | sed 's/""/"/g')
    while IFS=',' read -r risk_gene; do
        risk_gene=${risk_gene:1:-1}
        risk_gene=$(echo "$risk_gene" | sed 's/""/"/g')
	# Remove leading/trailing whitespace
        # clear_gene=$(echo "$clear_gene" | xargs | sed 's/\r$//')  
        # risk_gene=$(echo "$risk_gene" | xargs | sed 's/\r$//')
        sbatch --output="$output_log"<<EOL
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --job-name=lyso-sir-simulator
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gabriella.chan@monash.edu

# Load MATLAB module
module load matlab/r2023b

# Run main with each clearance gene
matlab -nodesktop -nodisplay -nosplash -r 'main($clear_gene, $risk_gene); exit;'

EOL
    done < "$risk_csv"
done < "$clear_csv"
