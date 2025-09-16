#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/step-7/logs/%x-%A.out
#SBATCH -e {{output}}/step-7/logs/%x-%A.err

conda activate {{conda_environment}}

cd {{output}}

cd 03B.DASTool_refinement/${name}/${name}_DASTool_bins
mkdir ../renamed
cp bc*.fa ../renamed/
cd ../renamed
rm -rf checkm_out
checkm lineage_wf ./ ./checkm_out -x fa -t 8 --tab_table -f ../${name}_checkm_table.txt --pplacer_threads 1
cd ../
cp ${name}_checkm_table.txt step-7
