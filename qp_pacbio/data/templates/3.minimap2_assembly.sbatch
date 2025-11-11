#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/step-3/logs/%x-%A_%a.out
#SBATCH -e {{output}}/step-3/logs/%x-%A_%a.err
#SBATCH --array {{array_params}}

source ~/.bashrc
set -e
{{conda_environment}}
cd {{output}}

step=${SLURM_ARRAY_TASK_ID}
input=$(head -n $step {{output}}/sample_list.txt | tail -n 1)
sample_name=`echo $input | awk '{print $1}'`
filename=`echo $input | awk '{print $2}'`

fn=`basename ${filename}`

# updating the GUI when task 1 runs
if [[ "$step" == "1" ]]; then
    python -c "from qp_pacbio.util import client_connect; qclient = client_connect('{{url}}'); qclient.update_job_step('{{qjid}}', 'Running step 3: ${SLURM_ARRAY_JOB_ID}')"
fi

folder=step-3/${sample_name}_binning
mkdir -p ${folder}

if [ -f step-2/${sample_name}_noLCG.fa ]; then
    minimap2 -x map-hifi -I {{mem_in_gb}}G -t {{nprocs}} -a --MD --eqx -o ${folder}/${sample_name}.sam step-2/${sample_name}_noLCG.fa ${filename}
    samtools view -bS -@4 ${folder}/${sample_name}.sam | samtools sort -@4 -O bam -o ${folder}/${sample_name}.sorted.bam
else
    touch {{output}}/step-4/completed_${SLURM_ARRAY_TASK_ID}.log
    touch {{output}}/step-5/completed_${SLURM_ARRAY_TASK_ID}.log
    touch {{output}}/step-6/completed_${SLURM_ARRAY_TASK_ID}.log
    touch {{output}}/step-7/completed_${SLURM_ARRAY_TASK_ID}.log
fi
touch {{output}}/step-3/completed_${SLURM_ARRAY_TASK_ID}.log
