#!/bin/bash

# Check if an input CSV file is provided as an argument
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 clearance_genes.csv [PARC] [SEED]"
    exit 1
fi

input_csv="$1"
parc="${2:-S132}"      # Default PARC is "S132"
seed="${3:-0}"         # 0 if seed needs to be assigned
output_log="slurm_output.log"

# Check input
if [ ! -e "$input_csv" ]; then
    echo "Input CSV file not found: $input_csv"
    exit 1
fi
if ["$seed" eq "0"]; then
    case "$parc" in
        DK) seed=41 ;;
        S132) seed=51 ;;
        S332) seed=151 ;;
        *)
            echo "Invalid PARC ID. Please use DK, S132, or S332."
            exit 1
            ;;
    esac
fi

# Loop through each line in the CSV file
while IFS=',' read -r clear_gene; do
    sbatch --output="$output_log" --job-name="$clear_gene-SIR-seed$seed"<<EOL
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:10:00
#SBATCH --mail-type=all
#SBATCH --mail-user=gabriella.chan@monash.edu

# Load MATLAB module
module load matlab/r2023b

# Run the MATLAB script with the arguments
matlab -nodesktop -nodisplay -nosplash -r 'main("$clear_gene", parc="$parc", seed=$seed, out="../results/test.csv"); exit;'

EOL
    done
done < "$input_csv"