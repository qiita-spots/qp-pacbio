#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/step-2/logs/%x-%A_%a.out
#SBATCH -e {{output}}/step-2/logs/%x-%A_%a.err
#SBATCH --array {{array_params}}

source ~/.bashrc
set -e
conda activate {{conda_environment}}
cd {{output}}/step-1

step=${SLURM_ARRAY_TASK_ID}
input=$(head -n $step {{output}}/sample_list.txt | tail -n 1)
sample_name=`echo $input | awk '{print $1}'`
filename=`echo $input | awk '{print $2}'`
fn=`basename ${filename}`

# updating the GUI when task 1 runs
if [[ "$step" == "1" ]]; then
    python -c "from qp_pacbio.util import client_connect; qclient = client_connect('{{url}}'); qclient.update_job_step('{{qjid}}', 'Running step 2: ${SLURM_ARRAY_JOB_ID}')"
fi

cat ${sample_name}.p_ctg.gfa | awk '$1=="S" && ($2 ~ /.c$/) {printf ">%s\n%s\n", $2, $3} ' > ../step-2/${sample_name}_circ.fa
seqkit split --by-id ../step-2/${sample_name}_circ.fa -O ../step-2/${sample_name}_split

### get all contigs for each sample
cat ${sample_name}.p_ctg.gfa | awk '$1=="S" {printf ">%s\n%s\n", $2, $3} ' > ../step-2/${sample_name}_all_contigs.fa

cd ../step-2/${sample_name}_split

### remove small circular genomes
find . -type f -size -512k -exec rm -f {} +
### extract fasta id for all the genomes in the split folder
for f in *.fa; do
    k=${f##*/}
    n=${f%.*}
    grep -E "^>" $f >> circular_id.txt
done

sed -i 's/>//' circular_id.txt
seqkit grep -v -f circular_id.txt ../${sample_name}_all_contigs.fa > ../${sample_name}_noLCG.fa

shopt -s nullglob

FILES=({{output}}/step-2/${sample_name}_split/*.fa)
if [ -f $FILES ]; then
    lcg_folder={{result_fp}}/${sample_name}/LCG/
    mkdir -p ${lcg_folder}
    for f in `ls {{output}}/step-2/${sample_name}_split/*.fa`; do
        sn=`basename ${f/_circ/}`;
        sn=${sn/part_/};
        cat $f | gzip > ${lcg_folder}/${sn/.fa/.fna}.gz;
    done
fi

if [ -f ${{output}}/step-2/${sample_name}_noLCG.fa ]; then
    cat {{output}}/step-2/${sample_name}_noLCG.fa | gzip > {{result_fp}}/${sample_name}/${sample_name}.noLCG.fna.gz
fi
