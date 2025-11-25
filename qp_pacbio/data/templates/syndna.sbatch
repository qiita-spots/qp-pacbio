#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/minimap2/logs/%x-%A_%a.out
#SBATCH -e {{output}}/minimap2/logs/%x-%A_%a.err
#SBATCH --array {{array_params}}

source ~/.bashrc
set -e
{{conda_environment}}
out_folder={{output}}/syndna
mkdir -p
cd ${out_folder}
db_folder=/scratch/qp-pacbio/minimap2/syndna/

step=${SLURM_ARRAY_TASK_ID}
input=$(head -n $step {{output}}/sample_list.txt | tail -n 1)
sample_name=`echo $input | awk '{print $1}'`
filename=`echo $input | awk '{print $2}'`
fn=`basename ${filename}`

mkdir -p ${out_folder}/clean/

# removing AllsynDNA_plasmids_FASTA_ReIndexed_FINAL.fasta not coverm
minimap2 -x map-hifi -t {{nprocs}} -a --MD --eqx -o ${out_folder}/${sample_name}_plasmid.sam ${db_folder}/AllsynDNA_plasmids_FASTA_ReIndexed_FINAL.fasta $filename
samtools view -F 4 -@ {{nprocs}} ${out_folder}/${sample_name}_plasmid.sam | awk '{print $1}' | sort -u > ${out_folder}/${sample_name}_plasmid_mapped.txt
seqkit grep -v -f ${out_folder}/${sample_name}_plasmid_mapped.txt $filename > ${out_folder}/${sample_name}_no_plasmid.fastq

# removing All_synDNA_inserts.fasta use coverm
minimap2 -x map-hifi -t {{nprocs}} -a --MD --eqx -o ${out_folder}/${sample_name}_inserts.sam ${db_folder}/All_synDNA_inserts.fasta ${out_folder}/${sample_name}_no_plasmid.fastq
samtools view -bS -@ {{int(nprocs/2)}} ${out_folder}/${sample_name}_no_plasmid.fastq | samtools sort -@ {{int(nprocs/2)}} -O bam -o ${out_folder}/${sample_name}_inserts_sorted.sam
coverm filter --bam-files ${out_folder}/${sample_name}_inserts_sorted.sam --min-read-percent-identity 99.9 --min-read-aligned-percent 95 --threads {{nprocs}} -o ${out_folder}/${sample_name}_no_inserts.bam
samtools view -O SAM -o ${out_folder}/${sample_name}_no_inserts_sorted.sam ${out_folder}/${sample_name}_no_inserts.bam
awk '{print $1}' ${out_folder}/${sample_name}_no_inserts_sorted.sam > ${out_folder}/${sample_name}_reads_filtered.txt
seqkit grep -v -f ${out_folder}/${sample_name}_reads_filtered.txt ${out_folder}/${sample_name}_no_plasmid.fastq > ${out_folder}/${sample_name}_no_plasmid_no_inserts.fastq

# removing GCF_000184185.1_ASM18418v1_genomic_chroso.fna use coverm
minimap2 -x map-hifi -t {{nprocs}} -a --MD --eqx -o ${out_folder}/${sample_name}_GCF_000184185.sam ${db_folder}/GCF_000184185.1_ASM18418v1_genomic_chroso.fna ${out_folder}/${sample_name}_no_plasmid_no_inserts.fastq
samtools view -bS -@ {{int(nprocs/2)}} ${out_folder}/${sample_name}_no_plasmid_no_inserts.fastq | samtools sort -@ {{int(nprocs/2)}} -O bam -o ${out_folder}/${sample_name}_GCF_000184185_sorted.sam
coverm filter --bam-files ${out_folder}/${sample_name}_GCF_000184185_sorted.sam --min-read-percent-identity 99.9 --min-read-aligned-percent 95 --threads {{nprocs}} -o ${out_folder}/${sample_name}_GCF_000184185.bam
samtools view -O SAM -o ${out_folder}/${sample_name}_no_GCF_000184185_sorted.sam ${out_folder}/${sample_name}_no_inserts.bam
awk '{print $1}' ${out_folder}/${sample_name}_no_GCF_000184185_sorted.sam > ${out_folder}/${sample_name}_GCF_000184185_reads_filtered.txt
seqkit grep -v -f ${out_folder}/${sample_name}_GCF_000184185_reads_filtered.txt ${out_folder}/${sample_name}_GCF_000184185.fastq | gz > ${out_folder}/clean/${fn}
