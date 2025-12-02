#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/merge/logs/%x-%A.out
#SBATCH -e {{output}}/merge/logs/%x-%A.err

source ~/.bashrc
set -e
{{conda_environment}}
cd {{output}}/
