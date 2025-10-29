#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/minimap2/logs/%x-%A_%a.out
#SBATCH -e {{output}}/minimap2/logs/%x-%A_%a.err
#SBATCH --array {{array_params}}

source ~/.bashrc
set -e
{{conda_environment}}
mkdir -p {{output}}/alignments
cd {{output}}/
db=/ddn_scratch/qiita_t/working_dir/tmp/db/WoLr2.mmi

step=${SLURM_ARRAY_TASK_ID}
input=$(head -n $step {{output}}/sample_list.txt | tail -n 1)

sample_name=`echo $input | awk '{print $1}'`
filename=`echo $input | awk '{print $2}'`

fn=`basename ${filename}`

minimap2 -x map-hifi -t {{nprocs}} -a \
       --secondary=no --MD --eqx ${db} \
       ${filename} | \
   samtools sort -@ {{nprocs}} - | \
   awk 'BEGIN { FS=OFS="\t" } /^@/ { print; next } { $10="*"; $11="*" } 1' | \
   xz -1 -T1 > {{output}}/alignments/${sample_name}.sam.xz
