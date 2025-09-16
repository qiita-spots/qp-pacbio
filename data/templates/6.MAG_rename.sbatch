#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/step-6/logs/%x-%A.out
#SBATCH -e {{output}}/step-6/logs/%x-%A.err

conda activate {{conda_environment}}
cd {{output}}

step=${SLURM_ARRAY_TASK_ID}
input=$(head -n $step {{output}}/sample_list.txt | tail -n 1)
fn=`basename ${input}`
folder=step-6/

cd 03B.DASTool_refinement/${name}/${name}_DASTool_bins

for f in *.fa; do
    k=${f##*/}
    a="${name}_${k}"
    cp "$k" "$a"
done
cp ./bc*.fa step-6/
