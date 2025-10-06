#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/step-4/logs/%x-%A_%a.out
#SBATCH -e {{output}}/step-4/logs/%x-%A_%a.out
#SBATCH --array {{array_params}}

conda activate {{conda_environment}}

cd {{output}}

step=${SLURM_ARRAY_TASK_ID}
input=$(head -n $step {{output}}/sample_list.txt | tail -n 1)
sample_name=`echo $input | awk '{print $1}'`
filename=`echo $input | awk '{print $2}'`

fn=`basename ${filename}`
folder=step-4/${sample_name}_binning

mkdir -p {{output}}/${folder}/work_files/
cp {{output}}/step-3/${sample_name}_binning/${sample_name}.sorted.bam {{output}}/step-4/${sample_name}_binning/work_files/${sample_name}.bam
rm {{output}}/step-3/${sample_name}_binning/${sample_name}.sorted.bam

ln -s ${filename} {{output}}/step-4/${sample_name}_binning/work_files/${sample_name}.fastq
metawrap binning -a {{output}}/step-2/${sample_name}_noLCG.fa -o {{output}}/step-4/${sample_name}_binning \
   -t {{nprocs}} -m 100 -l 16000 --single-end --metabat2 --maxbin2 --concoct --universal {{output}}/step-4/${sample_name}_binning/work_files/${sample_name}.fastq
