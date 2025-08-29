#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/step-1/logs/%x-%A.out
#SBATCH -e {{output}}/step-1/logs/%x-%A.err
#SBATCH --array {{array_params}}

conda activate {{conda-environment}}

cd {{output}}
step=${SLURM_ARRAY_TASK_ID}
input=$(head -n $step {{output}}/file_list.txt | tail -n 1)
fn=`basename ${input}`
hifiasm_meta -t 60 -o {{output}}/step-1/${fn} ${input}
