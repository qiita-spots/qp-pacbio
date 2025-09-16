#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/step-5/logs/%x-%A.out
#SBATCH -e {{output}}/step-5/logs/%x-%A.err

conda activate {{conda_environment}}
cd {{output}}

step=${SLURM_ARRAY_TASK_ID}
input=$(head -n $step {{output}}/sample_list.txt | tail -n 1)
fn=`basename ${input}`
folder=step-5/${fn}_binning
DAS_db=[path to db]

mkdir -p {{output}}/${folder}/work_files/

# rm -rf 02A_binning/${name}_binning/work_files
# cd 02A_binning/${name}_binning
# rm metabat2_bins/bin.unbinned.fa
# rm concoct_bins/unbinned.fa
Fasta_to_Contig2Bin.sh -i ./concoct_bins -e fa > ${name}.concoct.tsv
Fasta_to_Contig2Bin.sh -i ./maxbin2_bins -e fa > ${name}.maxbin2.tsv
Fasta_to_Contig2Bin.sh -i ./metabat2_bins -e fa > ${name}.metabat2.tsv
mkdir -p {{output}}/step-5/${name}

DAS_Tool --bins=${name}.concoct.tsv,${name}.maxbin2.tsv,${name}.metabat2.tsv --contigs=${name}.p_ctg_asm.fa --outputbasename=${name}/${name} --labels=CONCOCT,MaxBin,MetaBAT --threads=4 --search_engine=usearch --dbDirectory=${DAS_db} --write_bins
