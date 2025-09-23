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

# Default values
input_csv=""
yymmdd=$(date +%y%m%d)
out_dir="../results/$yymmdd/"
null_="none"
seed="51"
parc="S132"
time="0:30:00"

# Parse flags
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input) input_csv="$2"; shift 2 ;;
        --output) out_dir="$2"; shift 2 ;;
        --null) null_="$2"; shift 2 ;;
        --seed) seed="$2"; shift 2 ;;
        --parc) parc="$2"; shift 2 ;;
        --time) time="$2"; shift 2 ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 --input gene_pairs.csv [--output dir] [--null none] [--seed N] [--parc S132] [--time 1:00:00]"
            exit 1
            ;;
    esac
done

# Check required input
if [ -z "$input_csv" ]; then
    echo "Usage: $0 --input gene_pairs.csv [--output dir] [--null none] [--seed N] [--parc S132] [--time 1:00:00]"
    exit 1
fi

if [ ! -e "$input_csv" ]; then
    echo "Input CSV file not found: $input_csv"
    exit 1
fi

case "$parc" in
    DK) n_rois=41 ;;
    S132) n_rois=66 ;;
    S332) n_rois=132 ;;
    *)
        echo "Invalid PARC ID. Please use DK, S132, or S332."
        exit 1
        ;;
esac

if [[ "$out_dir" != */ ]]; then
    out_dir="$out_dir/"
fi
mkdir -p "$out_dir"

if [[ "$time" == "0:30:00" || "$time" < "0:30:00" ]]; then
    QOS="--qos=shortq"
else
    QOS=""
fi

# Read the input CSV and process each row
while IFS=',' read -r clear_gene risk_genes; do
    queue_jobs
	# Take first risk gene to differentiate between repeats
    risk_id=$(echo "$risk_genes" | cut -d',' -f1)
    # Format risk_genes into MATLAB string array syntax
    risk_names=$(echo "$risk_genes" | sed 's/,/","/g' | sed 's/^/["/' | sed 's/$/"]/')
    sbatch --output="${out_dir}_slurm_${seed}.log" --job-name="${clear_gene}_${risk_id}_${seed}" $QOS <<EOL
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=$time
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gabriella.chan@monash.edu

module load matlab/r2023b
matlab -nodesktop -nodisplay -nosplash -r 'main("$clear_gene", risk_names=$risk_names, null="$null_", pf=false, parc="$parc", seed=$seed, out="${out_dir}${parc}_${clear_gene}_${risk_id}_${null_}.csv"); exit;'
EOL

    if [ $? -ne 0 ]; then
        echo "Failed to submit job for $clear_gene with seed $seed"
        exit 1
    fi
done < "$input_csv"

