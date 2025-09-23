#!/bin/bash
# A job manager for SIR_SCZ simulation
# This script submits jobs to the SLURM scheduler for each gene/seed combination

function queue_jobs {
    while true; do
        n_jobs=$(squeue -u "$USER" | wc -l)
        if [ "$n_jobs" -lt 800 ]; then
            break
        fi
        sleep 10m
    done
}

function submit_job {
    clear_gene="$1"
    parc="$2"
    seed="$3"
    out_dir="$4"

#sbatch --output="${out_dir}_slurm_${seed}.log" --job-name="$clear_gene-seed$seed" --qos=shortq<<EOL
sbatch --output="${out_dir}_slurm_${seed}.log" --job-name="$clear_gene-seed$seed" <<EOL
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gabriella.chan@monash.edu

module load matlab/r2023b
matlab -nodesktop -nodisplay -nosplash -r 'main("$clear_gene", parc="$parc", seed=$seed, out="${out_dir}${parc}_${clear_gene}_${seed}.csv"); exit;'
EOL

    if [ $? -ne 0 ]; then
        echo "Failed to submit job for $clear_gene with seed $seed"
        exit 1
    fi
}

# Check arguments
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 clearance_genes.csv [PARC] [OUTPUT_DIR]"
    exit 1
fi

input_csv="$1"
parc="${2:-S132}"                       # Default "S132"
yymmdd=$(date +%y%m%d)
out_dir="${3:-../results/$yymmdd/}"     # Default "../results/YYMMDD/"

# Check input
if [ ! -e "$input_csv" ]; then
    echo "Input CSV file not found: $input_csv"
    exit 1
fi
case "$parc" in
    DK) n_rois=41 ;;
    S132) n_rois=66 ;;
    S332) n_rois=166 ;;
    *)
        echo "Invalid PARC ID. Please use DK, S132, or S332."
        exit 1
        ;;
esac
if [[ "$out_dir" != */ ]]; then
    out_dir="$out_dir/"
fi
mkdir -p "$out_dir"

# Loop through each gene/seed and submit
while IFS=',' read -r clear_gene; do
    queue_jobs
    # for seed in $(seq $n_rois); do
    for seed in $(seq 1 22); do
        submit_job "$clear_gene" "$parc" "$seed" "$out_dir"
    done
done < "$input_csv"
