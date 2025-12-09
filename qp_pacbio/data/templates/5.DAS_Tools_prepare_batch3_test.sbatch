#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/step-5/logs/%x-%A_%a.out
#SBATCH -e {{output}}/step-5/logs/%x-%A_%a.err
#SBATCH --array {{array_params}}

source ~/.bashrc
set -e
{{conda_environment}}
cd {{output}}

step=${SLURM_ARRAY_TASK_ID}
input=$(head -n $step {{output}}/sample_list.txt | tail -n 1)
sample_name=`echo $input | awk '{print $1}'`
filename=`echo $input | awk '{print $2}'`

fn=`basename ${filename}`

# updating the GUI when task 1 runs
if [[ "$step" == "1" ]]; then
    python -c "from qp_pacbio.util import client_connect; qclient = client_connect('{{url}}'); qclient.update_job_step('{{qjid}}', 'Running step 5: ${SLURM_ARRAY_JOB_ID}')"
fi

folder=step-5/${sample_name}_refinement
mkdir -p ${folder}

### make sure qiita have access to this DB path. Otherwise copy the folder to qiita-owned folders
DAS_db=/scratch/qp-pacbio/DAS_Tool/db/
cd {{output}}/step-4/${sample_name}_binning

rm metabat2_bins/bin.unbinned.fa
rm metabat2_bins/bin.lowDepth.fa
rm metabat2_bins/bin.tooShort.fa
rm metabat2_bins/bin.BinInfo.txt
rm concoct_bins/unbinned.fa

Fasta_to_Contig2Bin.sh -i ./concoct_bins -e fa > ${sample_name}.concoct.tsv
Fasta_to_Contig2Bin.sh -i ./maxbin2_bins -e fa > ${sample_name}.maxbin2.tsv
Fasta_to_Contig2Bin.sh -i ./metabat2_bins -e fa | awk 'BEGIN{FS=OFS="\t"}{print $1"\t"$4}' > ${sample_name}.metabat2.tsv

mkdir {{output}}/${folder}/${sample_name}

dastool_error_handler () {
    # if we see the expected error, ignore it and complete this
    # worker, otherwise raise an error
    if grep -qF "No bins with bin-score >0.5 found" dastool.log; then
        touch {{output}}/step-5/completed_${SLURM_ARRAY_TASK_ID}.log
        touch {{output}}/step-6/completed_${SLURM_ARRAY_TASK_ID}.log
        touch {{output}}/step-7/completed_${SLURM_ARRAY_TASK_ID}.log
    else
        cat dastool.log
        exit 1
    fi
}

trap dastool_error_handler ERR

DAS_Tool --bins=${sample_name}.concoct.tsv,${sample_name}.maxbin2.tsv,${sample_name}.metabat2.tsv --contigs={{output}}/step-2/${sample_name}_noLCG.fa --outputbasename={{output}}/${folder}/${sample_name}/${sample_name} --labels=CONCOCT,MaxBin,MetaBAT --threads={{nprocs}} --search_engine=diamond --dbDirectory=${DAS_db} --write_bins > dastool.log 2>&1

touch {{output}}/step-5/completed_${SLURM_ARRAY_TASK_ID}.log
