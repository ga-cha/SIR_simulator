#!/bin/bash

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 genes_pairs.csv"
  exit 1
fi
input_csv="$1"
if [ ! -e "$input_csv" ]; then
  echo "Input CSV file not found: $input_csv"
  exit 1
fi

while IFS=', ' read -r clear_gene risk_gene; do

    # Submit the job with the paired genes
    sbatch --output="slurm_output.log" --job-name="$clear_gene-SIR" <<EOL
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:08:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gabriella.chan@monash.edu

# Load MATLAB module
module load matlab/r2023b

# Run the MATLAB script with the arguments
matlab -nodesktop -nodisplay -nosplash -r "main(\"$clear_gene\", risk_names=\"$risk_gene\", parc=\"S132\", null=\"spatial\", out=\"../../SIR_results/lme_beta/S132_spatial.csv\"); exit;"
EOL
done < "$input_csv"
