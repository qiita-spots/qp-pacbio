#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/step-2/logs/%x-%A.out
#SBATCH -e {{output}}/step-2/logs/%x-%A.err

cd {{output}}/step-1
for f in *.fa; do
    # k=${f##*/};
    # n=${f%.*};
    echo seqkit split --by-id $f -O ../step-2/split
done | parallel -j {{nprocs}}

### remove small circular genomes
cd split
find . -type f -size -512k -exec rm -f {} +
### extract fasta id for all the genomes in the split folder
for f in *.fa; do
    k=${f##*/}
    n=${f%.*}
    grep -E "^>" $f >> circular_id.txt
done
sed -i 's/>//' circular_id.txt
seqkit -v -f circular_id.txt all_contigs.fa > noLCG.fa
