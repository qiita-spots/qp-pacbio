#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/step-6/logs/%x-%A.out
#SBATCH -e {{output}}/step-6/logs/%x-%A.err
#SBATCH --array {{array_params}}

conda activate {{conda_environment}}
cd {{output}}

step=${SLURM_ARRAY_TASK_ID}
input=$(head -n $step {{output}}/sample_list.txt | tail -n 1)
sample_name=`echo $input | awk '{print $1}'`
filename=`echo $input | awk '{print $2}'`

fn=`basename ${filename}`
folder=step-6/${fn}_MAG_rename
mkdir -p ${folder}

cd {{output}}/step-5/${fn}_refinement/${fn}/${fn}_DASTool_bins

for f in *.fa; do
    k=${f##*/}
    a="${filename}_${k}"
    mv "$k" "$a"
done
cp ./${filename}*.fa step-6/${fn}_MAG_rename

#### where do we store the MAGS?
