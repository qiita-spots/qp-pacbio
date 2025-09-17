#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/step-2/logs/%x-%A.out
#SBATCH -e {{output}}/step-2/logs/%x-%A.err
#SBATCH --array {{array_params}}
{% if dependency %}#SBATCH --dependency={{ dependency }}{% endif %}

conda activate {{conda_environment}}

cd {{output}}/step-1

step=${SLURM_ARRAY_TASK_ID} ##1000_2, 1000_1

input=$(head -n $step {{output}}/sample_list.txt | tail -n 1)
fn=`basename ${input}`

cat ${fn}.p_ctg.gfa | awk '$1=="S" && ($2 ~ /.c$/) {printf ">%s\n%s\n", $2, $3} ' > ../step-2/${fn}_circ.fa
seqkit split --by-id ../step-2/${fn}_circ.fa -O ../step-2/${fn}_split

### get all contigs for each sample
cat ${fn}.p_ctg.gfa | awk '$1=="S" {printf ">%s\n%s\n", $2, $3} ' > ../step-2/${fn}_all_contigs.fa

cd ../step-2/${fn}_split

### remove small circular genomes
find . -type f -size -512k -exec rm -f {} +
### extract fasta id for all the genomes in the split folder
for f in *.fa; do
    k=${f##*/}
    n=${f%.*}
    grep -E "^>" $f >> circular_id.txt
done

sed -i 's/>//' circular_id.txt
seqkit grep -v -f circular_id.txt ../step-2/${fn}_all_contigs.fa > ../step-2/${fn}_noLCG.fa
