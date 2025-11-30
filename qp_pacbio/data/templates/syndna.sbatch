#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/syndna/logs/%x-%A_%a.out
#SBATCH -e {{output}}/syndna/logs/%x-%A_%a.err
#SBATCH --array {{array_params}}

source ~/.bashrc
set -e
{{conda_environment}}
out_folder={{output}}/syndna
mkdir -p ${out_folder}
cd ${out_folder}
db_folder=/scratch/qp-pacbio/minimap2/syndna/

step=${SLURM_ARRAY_TASK_ID}
input=$(head -n $step {{output}}/sample_list.txt | tail -n 1)
sample_name=`echo $input | awk '{print $1}'`
filename=`echo $input | awk '{print $2}'`
fn=`basename ${filename}`

mkdir -p ${out_folder}/filtered/

sn_folder=${out_folder}/bioms/${sample_name}
mkdir -p ${sn_folder}

txt=${sn_folder}/${sample_name}.txt
tsv=${txt/.txt/.tsv}
coverm contig --single $filename --reference ${db_folder}/All_synDNA_inserts.fasta --mapper minimap2-hifi \
    --min-read-percent-identity 0.95 --min-read-aligned-percent 0.0 -m mean count --threads {{nprocs}} \
    --output-file ${txt}

awk 'BEGIN {FS=OFS="\t"}; {print $1,$3}' ${txt} | \
    sed 's/Contig/\#OTU ID/' | sed  's/All_synDNA_inserts.fasta\///' | \
    sed  's/ Read Count//' | sed "s/${fn}/${sample_name}/" > ${tsv}

# if counts is zero mark it as missing and stop
counts=`tail -n +2 ${tsv} | awk '{sum += $NF} END {print sum}'`
if [[ "$counts" == "0" ]]; then
    touch {{output}}/failed_${SLURM_ARRAY_TASK_ID}.log
    exit 0
fi

biom convert -i ${sn_folder}/${sample_name}.txt -o ${sn_folder}/${sample_name}.biom --to-hdf5

# removing AllsynDNA_plasmids_FASTA_ReIndexed_FINAL.fasta not coverm
minimap2 -x map-hifi -t {{nprocs}} -a --MD --eqx -o ${out_folder}/${sample_name}_plasmid.sam ${db_folder}/AllsynDNA_plasmids_FASTA_ReIndexed_FINAL.fasta $filename
samtools view -F 4 -@ {{nprocs}} ${out_folder}/${sample_name}_plasmid.sam | awk '{print $1}' | sort -u > ${out_folder}/${sample_name}_plasmid_mapped.txt
seqkit grep -v -f ${out_folder}/${sample_name}_plasmid_mapped.txt $filename > ${out_folder}/${sample_name}_no_plasmid.fastq

# removing GCF_000184185.1_ASM18418v1_genomic_chroso.fna use coverm
minimap2 -x map-hifi -t {{nprocs}} -a --MD --eqx -o ${out_folder}/${sample_name}_GCF_000184185.sam ${db_folder}/GCF_000184185.1_ASM18418v1_genomic_chroso.fna ${out_folder}/${sample_name}_no_plasmid.fastq
samtools view -bS -@ {{ nprocs/2 | int }} ${out_folder}/${sample_name}_no_plasmid.fastq | samtools sort -@ {{ nprocs/2 | int }} -O bam -o ${out_folder}/${sample_name}_GCF_000184185_sorted.bam
coverm filter --bam-files ${out_folder}/${sample_name}_GCF_000184185_sorted.bam --min-read-percent-identity 99.9 --min-read-aligned-percent 95 --threads {{nprocs}} -o ${out_folder}/${sample_name}_GCF_000184185.bam
samtools view -O SAM -o ${out_folder}/${sample_name}_no_GCF_000184185_sorted.sam ${out_folder}/${sample_name}_GCF_000184185.bam
awk '{print $1}' ${out_folder}/${sample_name}_no_GCF_000184185_sorted.sam > ${out_folder}/${sample_name}_GCF_000184185_reads_filtered.txt
seqkit grep -v -f ${out_folder}/${sample_name}_GCF_000184185_reads_filtered.txt ${out_folder}/${sample_name}_no_plasmid.fastq | gz > ${out_folder}/filtered/${fn}

touch {{output}}/completed_${SLURM_ARRAY_TASK_ID}.log
