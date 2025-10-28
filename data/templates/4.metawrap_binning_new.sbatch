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
conda activate {{conda_environment}}
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

folder=step-4/${sample_name}_binning

mkdir -p {{output}}/${folder}/work_files/
cp {{output}}/step-3/${sample_name}_binning/${sample_name}.sorted.bam {{output}}/step-4/${sample_name}_binning/work_files/${sample_name}.bam
rm {{output}}/step-3/${sample_name}_binning/${sample_name}.sorted.bam

ln -s ${filename} {{output}}/step-4/${sample_name}_binning/work_files/${sample_name}.fastq
metawrap binning -a {{output}}/step-2/${sample_name}_noLCG.fa -o {{output}}/step-4/${sample_name}_binning \
   -t {{nprocs}} -m 100 -l 16000 --single-end --metabat2 --maxbin2 --concoct --universal {{output}}/step-4/${sample_name}_binning/work_files/${sample_name}.fastq
rm -rf {{output}}/step-4/${sample_name}_binning/work_files
