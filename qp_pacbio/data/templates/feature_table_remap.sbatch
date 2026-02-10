#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/remap/logs/%x-%A_%a.out
#SBATCH -e {{output}}/remap/logs/%x-%A_%a.err
#SBATCH --array {{array_params}}

source ~/.bashrc
set -e
{{conda_environment}}

step=${SLURM_ARRAY_TASK_ID}
input=$(head -n $step {{output}}/sample_list.txt | tail -n 1)
sample_name=`echo $input | awk '{print $1}'`
filename=`echo $input | awk '{print $2}'`
aid=`echo $input | awk '{print $3}'`
fn=`basename ${filename}`

out_folder={{output}}/remap
sn_folder=${out_folder}/bioms/${aid}.${sample_name}
export TMPDIR=${out_folder}/tmp
mkdir -p ${TMPDIR} ${sn_folder}

txt=${sn_folder}/${aid}.${sample_name}.txt
tsv=${txt/.txt/.tsv}
coverm genome -d {{output}}/merge/dereplicated -x fna --mapper minimap2-hifi -m covered_bases trimmed_mean mean length count \
    --trim-min 0.05 --trim-max 0.95 --min-read-percent-identity {{percent_identity}} \
    --min-read-aligned-percent 0.90 --min-covered-fraction 0 -t {{nprocs}} \
    --single ${filename} --output-file ${txt}

awk 'BEGIN {FS=OFS="\t"}; {print $1,$NF}' ${txt} | \
    sed  's/dereplicated_gtdbtk\///' | \
    sed  's/ Read Count//' | sed "s/${fn}/${aid}.${sample_name}/" > ${tsv}

# saving coverages
cov_tsv=${txt/.txt/.cov.tsv}
echo -e "mag_id\tpercent_covered" > ${cov_tsv}
tail -n +2 ${txt} | awk -F'\t' 'BEGIN {FS=OFS="\t"}; { if ($5 > 0) print $1,$2/$5*100 }' >> ${cov_tsv}

# if counts is zero do not try to convert
counts=`tail -n +2 ${tsv} | awk '{sum += $NF} END {print sum}'`
if [[ "$counts" != "0" ]]; then
    biom convert -i ${tsv} -o ${sn_folder}/counts.biom --to-hdf5
    biom convert -i ${cov_tsv} -o ${sn_folder}/coverage.biom --to-hdf5
fi

touch ${out_folder}/completed_${SLURM_ARRAY_TASK_ID}.log
