#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/processing/logs/%x-%A_%a.out
#SBATCH -e {{output}}/processing/logs/%x-%A_%a.err
#SBATCH --array {{array_params}}

source ~/.bashrc
set -e
{{conda_environment}}
out_folder={{output}}/processing/
mkdir -p ${out_folder}/lima ${out_folder}/final

step=${SLURM_ARRAY_TASK_ID}
input=$(head -n $step {{output}}/sample_list.txt | tail -n 1)
sample_name=`echo $input | awk '{print $1}'`
filename=`echo $input | awk '{print $2}'`
fn=`basename ${filename}`
fout=${out_folder}/lima/${fn}
final=${out_folder}/final/${fn}

{{lima_cmd}}

# let's check that we actually removed some files
column=$(head -1 "${fout}.lima.counts" | tr '\t' '\n' | grep -En "Counts" | cut -d: -f1)
counts=$(cut -f "${column}" "${fout}.lima.counts" | tail -n 1)

# if not zero, then run the next command, else copy input
if [[ "$counts" -gt "0" ]]; then
    pbmarkdup "${fout}.fastq.gz" "${final}/${fn}.fastq.gz"
else
    cp "${fout}.fastq.gz" "${final}.fastq.gz"
fi

touch {{output}}/completed_${SLURM_ARRAY_TASK_ID}.log
