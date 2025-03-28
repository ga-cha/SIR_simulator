#!/bin/bash

# Check if an input CSV file is provided as an argument
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 clearance_genes.csv"
  exit 1
fi

input_csv="$1"
output_log="slurm_output.log"

# Check if the input CSV file exists
if [ ! -e "$input_csv" ]; then
  echo "Input CSV file not found: $input_csv"
  exit 1
fi

# Loop through each line in the CSV file
while IFS=',' read -r clear_gene; do
  sbatch --output="$output_log" --job-name="$clear_gene-SIR"<<EOL
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=gabriella.chan@monash.edu


# Load MATLAB module
module load matlab/r2023b

# Run the MATLAB script with the arguments
matlab -nodesktop -nodisplay -nosplash -r 'main("$clear_gene", parc="S132", out="../../SIR_results/lme_beta/S132_spatial.csv", null="spatial"); exit;'

EOL
done < "$input_csv"

