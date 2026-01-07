#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/step-4/logs/%x-%A_%a.out
#SBATCH -e {{output}}/step-4/logs/%x-%A_%a.err
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
    python -c "from qp_pacbio.util import client_connect; qclient = client_connect('{{url}}'); qclient.update_job_step('{{qjid}}', 'Running step 4: ${SLURM_ARRAY_JOB_ID}')"
fi

bfolder={{output}}/step-4/${sample_name}_binning
folder=${bfolder}/work_files/
mkdir -p ${folder}
cp {{output}}/step-3/${sample_name}_binning/${sample_name}.sorted.bam ${folder}/${sample_name}.bam

ln -s ${filename} ${folder}/${sample_name}.fastq

metawrap_error_handler () {
    # there has been a couple of occasions that metawrap has failed due to
    # "MaxBin2 did not produce a single bin"; this can be real or a bug - we are not sure
    # however, during our testing we saw that when it failed unexpectedly it will _not_
    # create a bin.log file in the work_folder
    if [ -f ${folder}/maxbin2_out/bin.log ]; then
        touch {{output}}/step-4/completed_${SLURM_ARRAY_TASK_ID}.log
        touch {{output}}/step-5/completed_${SLURM_ARRAY_TASK_ID}.log
        touch {{output}}/step-6/completed_${SLURM_ARRAY_TASK_ID}.log
        touch {{output}}/step-7/completed_${SLURM_ARRAY_TASK_ID}.log
    else
        exit 1
    fi
}

trap metawrap_error_handler ERR

metawrap binning -a {{output}}/step-2/${sample_name}_noLCG.fa -o ${bfolder} \
   -t {{nprocs}} -m 100 -l 16000 --single-end --metabat2 --maxbin2 --concoct --universal ${folder}/${sample_name}.fastq

touch {{output}}/step-4/completed_${SLURM_ARRAY_TASK_ID}.log
