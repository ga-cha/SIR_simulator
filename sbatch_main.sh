#!/bin/bash

if [ "$#" -ne 1 ]; then
	echo "Usage: $0 input.csv"
	exit 1
fi

input_csv="$1"
output_log="slurm_output.log"

if [ ! -e "$input_csv" ]; then
	echo "Input csv not found: $input_csv"
	exit 1
fi

while IFS=',' read -r clear_gene; do
	sbatch --output="$output_log"<<EOL
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20:00:00
#SBATCH --job-name=sir-simulator
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gabriella.chan@monash.edu

# Load MATLAB module
module load matlab/r2023a

# Run main with each clearance gene
matlab -nodisplay -r "main('$clear_gene'); exit;"
EOL
done < "$input_csv"
