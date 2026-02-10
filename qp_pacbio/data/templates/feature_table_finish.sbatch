#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/finish/logs/%x-%A_%a.out
#SBATCH -e {{output}}/finish/logs/%x-%A_%a.err

source ~/.bashrc
set -e
{{conda_environment}}
cd {{output}}/

python -c "from qp_pacbio.util import client_connect; qclient = client_connect('{{url}}'); qclient.update_job_step('{{qjid}}', 'Merging BIOMs')"

lmap=${PWD}/length.map
ls ${PWD}/*/bioms/*/*.cov > ${PWD}/cov_files.txt
seqkit stats merge/dereplicated/*.fna | tail -n +2 | awk 'BEGIN{OFS="\t"} {print $1, $5}' | sed 's/,//g' | sed 's#merge/dereplicated/##g' | sed 's/.fna//g' > ${lmap}
micov consolidate --paths ${PWD}/cov_files.txt --lengths ${lmap} --output ${PWD}/coverages.tgz

biom_merge_pacbio --base {{output}}/remap --merge-type counts
biom_merge_pacbio --base {{output}}/remap --merge-type coverage

finish_qp_pacbio {{url}} {{qjid}} {{output}}
