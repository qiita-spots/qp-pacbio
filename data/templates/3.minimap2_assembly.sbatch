#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/step-3/logs/%x-%A.out
#SBATCH -e {{output}}/step-3/logs/%x-%A.err
{% if dependency %}#SBATCH --dependency={{ dependency }}{% endif %}


conda activate {{conda_environment}}
cd {{output}}

step=${SLURM_ARRAY_TASK_ID}
input=$(head -n $step {{output}}/sample_list.txt | tail -n 1)
fn=`basename ${input}`
folder=step-3/${fn}_binning

mkdir -p ${folder}
minimap2 -x map-hifi -t {{nprocs}} -a --MD --eqx -o ${folder}/${fn}.sam ./02A_assembly_hifiasm_meta_contigs/${name}.fa  ./01B.clean_reads/${name}.fastq
samtools view -bS -@4 02A_binning/${name}_binning/${name}.sam | samtools sort -@4 -O bam -o 02A_binning/${name}_binning/${name}.sorted.bam
