#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/step-5/logs/%x-%A.out
#SBATCH -e {{output}}/step-5/logs/%x-%A.err
#SBATCH --array {{array_params}}

conda activate {{conda_environment}}
cd {{output}}

step=${SLURM_ARRAY_TASK_ID}
input=$(head -n $step {{output}}/sample_list.txt | tail -n 1)
sample_name=`echo $input | awk '{print $1}'`
filename=`echo $input | awk '{print $2}'`

fn=`basename ${filename}`
folder=step-5/${fn}_refinement

mkdir -p ${folder}

DAS_db=/ddn_scratch/jiz322/THDMI_binning/DAS_Tool/db
cd {{output}}/step-4/${fn}_binning

rm metabat2_bins/bin.unbinned.fa
rm concoct_bins/unbinned.fa
Fasta_to_Contig2Bin.sh -i ./concoct_bins -e fa > ${fn}.concoct.tsv
Fasta_to_Contig2Bin.sh -i ./maxbin2_bins -e fa > ${fn}.maxbin2.tsv
Fasta_to_Contig2Bin.sh -i ./metabat2_bins -e fa > ${fn}.metabat2.tsv

DAS_Tool --bins=${fn}.concoct.tsv,${fn}.maxbin2.tsv,${fn}.metabat2.tsv --contigs={{output}}/step-2/${fn}.noLCG.fa --outputbasename={{output}}/${folder}/${fn}/${fn} --labels=CONCOCT,MaxBin,MetaBAT --threads={{nprocs}} --search_engine=diamond --dbDirectory=${DAS_db} --write_bins