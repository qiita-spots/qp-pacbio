#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/step-6/logs/%x-%A_%a.out
#SBATCH -e {{output}}/step-6/logs/%x-%A_%a.err
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

# updating the GUI when task 1 runs
if [[ "$step" == "1" ]]; then
    python -c "from qp_pacbio.util import client_connect; qclient = client_connect('{{url}}'); qclient.update_job_step('{{qjid}}', 'Running step 6: ${SLURM_ARRAY_JOB_ID}')"
fi

folder=step-6/${sample_name}_MAG_rename
mkdir -p ${folder}

cd {{output}}/step-5/${sample_name}_refinement/${sample_name}/${sample_name}_DASTool_bins

for f in *.fa; do
    k=${f##*/}
    a="${sample_name}_${k}"
    mv "$k" "$a"
done
cp ${sample_name}*.fa {{output}}/step-6/${sample_name}_MAG_rename

mag_folder={{result_fp}}/${sample_name}/MAG/
mkdir -p ${mag_folder}
for f in `ls step-6/${sample_name}_MAG_rename/${sample_name}*.fa`; do
    sn=`basename ${f}`;
    cat $f | gzip > ${mag_folder}/${sn/.fa/.fna}.gz;
done
