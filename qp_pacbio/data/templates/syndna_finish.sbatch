#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/merge/logs/%x-%A_%a.out
#SBATCH -e {{output}}/merge/logs/%x-%A_%a.err

source ~/.bashrc
set -e
{{conda_environment}}
cd {{output}}/



mkdir -p ${out_folder}/filtered/
sn_folder=${out_folder}/${sample_name}
mkdir -p ${sn_folder}
coverm contig --single $filename --reference ${db_folder}/All_synDNA_inserts.fasta --mapper minimap2-hifi \
    --min-read-percent-identity 0.95 --min-read-aligned-percent 0.0 -m mean count --threads {{nprocs}} \
    --output-file ${sn_folder}/${sample_name}_insert_counts.txt
cat ${sn_folder}/${sample_name}_insert_counts.txt | sed 's/Contig/\#OTU ID/' | \
    sed 's/ Read Count//' > ${sn_folder}/${sample_name}_insert_counts.tsv
biom convert -i ${sn_folder}/${sample_name}_insert_counts.txt -o ${sn_folder}/${sample_name}_insert_counts.biom --to-hdf5

biom_merge_pacbio --base {{output}} --type syndna

find {{output}}/coverages/ -iname "*.cov" > {{output}}/cov_files.txt
micov consolidate --paths {{output}}/cov_files.txt --lengths ${len_map} --output {{output}}/coverages.tgz

cd alignment
tar -cvf ../alignment.tar *.sam.xz

finish_qp_pacbio {{url}} {{qjid}} {{output}}
