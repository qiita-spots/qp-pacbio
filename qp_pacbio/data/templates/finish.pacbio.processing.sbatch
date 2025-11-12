#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/finish/logs/finish.out
#SBATCH -e {{output}}/finish/logs/finish.err

source ~/.bashrc
set -e
{{conda_environment}}
cd {{output}}

tar zcvf results/logs.tgz step*/logs/

finish_qp_pacbio {{url}} {{qjid}} {{output}}
