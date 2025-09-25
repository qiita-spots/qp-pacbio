#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/step-1/logs/%x-%A_%a.out
#SBATCH -e {{output}}/step-1/logs/%x-%A_%a.out
#SBATCH --array {{array_params}}

conda activate {{conda_environment}}

cd {{output}}/step-1
step=${SLURM_ARRAY_TASK_ID}
input=$(head -n $step {{output}}/sample_list.txt | tail -n 1)

sample_name=`echo $input | awk '{print $1}'`
filename=`echo $input | awk '{print $2}'`

fn=`basename ${filename}`

hifiasm_meta -t {{nprocs}} -o {{output}}/step-1/${fn} ${filename}
