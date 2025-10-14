#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/step-5/logs/%x-%A_%a.out
#SBATCH -e {{output}}/step-5/logs/%x-%A_%a.out
#SBATCH --array {{array_params}}

source ~/.bashrc

conda activate {{conda_environment}}
cd {{output}}

step=${SLURM_ARRAY_TASK_ID}
input=$(head -n $step {{output}}/sample_list.txt | tail -n 1)
sample_name=`echo $input | awk '{print $1}'`
filename=`echo $input | awk '{print $2}'`

fn=`basename ${filename}`
folder=step-5/${sample_name}_refinement

mkdir -p ${folder}


### make sure qiita have access to this DB path. Otherwise copy the folder to qiita-owned folders
DAS_db=/ddn_scratch/qiita/databases/qp-pacbio/DAS_Tool/db/
cd {{output}}/step-4/${sample_name}_binning

rm metabat2_bins/bin.unbinned.fa
rm concoct_bins/unbinned.fa
Fasta_to_Contig2Bin.sh -i ./concoct_bins -e fa > ${sample_name}.concoct.tsv
Fasta_to_Contig2Bin.sh -i ./maxbin2_bins -e fa > ${sample_name}.maxbin2.tsv
Fasta_to_Contig2Bin.sh -i ./metabat2_bins -e fa > ${sample_name}.metabat2.tsv

DAS_Tool --bins=${sample_name}.concoct.tsv,${sample_name}.maxbin2.tsv,${sample_name}.metabat2.tsv --contigs={{output}}/step-2/${sample_name}_noLCG.fa --outputbasename={{output}}/${folder}/${sample_name}/${sample_name} --labels=CONCOCT,MaxBin,MetaBAT --threads={{nprocs}} --search_engine=diamond --dbDirectory=${DAS_db} --write_bins