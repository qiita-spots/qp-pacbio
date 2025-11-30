#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/finish/logs/%x-%A_%a.out
#SBATCH -e {{output}}/finish/logs/%x-%A_%a.err

source ~/.bashrc
set -e
{{conda_environment}}
cd {{output}}/

biom_merge_pacbio --base {{output}}/syndna --type syndna

# find {{output}}/coverages/ -iname "*.cov" > {{output}}/cov_files.txt
# micov consolidate --paths {{output}}/cov_files.txt --lengths ${len_map} --output {{output}}/coverages.tgz

# cd alignment
# tar -cvf ../alignment.tar *.sam.xz

finish_qp_pacbio {{url}} {{qjid}} {{output}}
