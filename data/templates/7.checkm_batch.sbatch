#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/step-7/logs/%x-%A_%a.out
#SBATCH -e {{output}}/step-7/logs/%x-%A_%a.out
#SBATCH --array {{array_params}}


source ~/.bashrc
set -e
conda activate {{conda_environment}}
cd {{output}}

step=${SLURM_ARRAY_TASK_ID}
input=$(head -n $step {{output}}/sample_list.txt | tail -n 1)
sample_name=`echo $input | awk '{print $1}'`
filename=`echo $input | awk '{print $2}'`

fn=`basename ${filename}`

folder=step-7/checkm_out
mkdir -p ${folder}

cd step-5/${sample_name}_refinement/${sample_name}/${sample_name}_DASTool_bins

checkm lineage_wf ./ ./checkm_out -x fa -t 8 --tab_table -f ../${sample_name}_checkm_table.txt --pplacer_threads 2
cd ../
cp ${sample_name}_checkm_table.txt {{output}}/step-7/checkm_out/
