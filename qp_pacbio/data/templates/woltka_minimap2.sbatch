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
mkdir -p {{output}}/alignment {{output}}/filtered-alignment
cd {{output}}/
db=/scratch/qp-pacbio/minimap2/WoLr2/WoLr2.map-hifi.mmi

step=${SLURM_ARRAY_TASK_ID}
input=$(head -n $step {{output}}/sample_list.txt | tail -n 1)

sample_name=`echo $input | awk '{print $1}'`
filename=`echo $input | awk '{print $2}'`

fn=`basename ${filename}`

minimap2 -x map-hifi -t {{nprocs}} -a \
       --secondary=no --MD --eqx ${db} \
       ${filename} | \
   awk 'BEGIN { FS=OFS="\t" } /^@/ { print; next } { $10="*"; $11="*" } 1' | grep -v "^@" | sort -k 1 | \
   xz -1 -T1 > {{output}}/alignment/${sample_name}.sam.xz

xcat {{output}}/alignment/${sample_name}.sam.xz | \
    /home/mcdonadt/duckdb-miint/duckdb-2025.11.20 -f {{sql_home}}/qcov-seqident-filter.sql > {{output}}/filtered-alignment/${sample_name}.sam.gz
