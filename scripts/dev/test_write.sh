#!/bin/bash

out=sbatch_out.log


for i in {100..120}
do
    sbatch --job-name="$i write_async" --output="$out" <<EOL
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:02:00

module load matlab/r2023b

matlab -nodesktop -nodisplay -nosplash -r 'test_write($i); exit;'

EOL
done
