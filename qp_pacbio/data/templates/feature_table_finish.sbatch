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

mkdir -p coverages
cp */bioms/*/*.cov coverages
tar czvf coverages.tgz coverages

biom_merge_pacbio --base {{output}}/remap --merge-type counts

finish_qp_pacbio {{url}} {{qjid}} {{output}}
