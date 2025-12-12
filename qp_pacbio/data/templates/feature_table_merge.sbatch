#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/merge/logs/%x-%A.out
#SBATCH -e {{output}}/merge/logs/%x-%A.err

source ~/.bashrc
set -e
{{conda_environment}}
cd {{output}}/merge

python -c "from qp_pacbio.util import client_connect; qclient = client_connect('{{url}}'); qclient.update_job_step('{{qjid}}', 'Merging LCG/MAG files')"

# The expected files here are:
# {aid}_prep_info_artifact.tsv : preparation for each artifact
# {aid}_folders.tsv : LCG/MAG folders to check for eaach artifact
# {aid}_{pid}_sample_list.txt : the sample/files list of the artifact, {pid} is the parent id

# Step 1: create processing folders, note that all files need to be gunzip
mkdir -p LCG all_fna checkm

for file in $(ls ../*_folders.tsv); do
    for folder in $(cat ${file}); do
        ffn=$(basename $folder);
        aid=${ffn%%_*};
        # the pattern of the folers is [folder]/results/[LCG|MAG]/[sample-id]/
        for f in $(ls ${folder}/*/*/LCG/*fna.gz 2> /dev/null || true); do
            bn=$(basename ${f/.gz/});
            echo "gunzip -c ${f} > LCG/${aid}_${bn}"
        done
        for f in $(ls ${folder}/*/*/MAG/*fna.gz 2> /dev/null || true); do
            bn=$(basename ${f/.gz/});
            echo "gunzip -c ${f} > all_fna/${aid}_${bn}"
        done
        # the pattern for txts is [folder]/results/*checkm.txt.gz
        for f in $(ls ${folder}/*/*/*checkm.txt.gz); do
            bn=$(basename ${f/.gz/});
            echo "zcat ${f} > checkm/${aid}_${bn}";
        done
    done
done | parallel --halt now,fail=1 -j {{nprocs}}

# merging the checkm files
head -n 1 $(ls checkm/*.txt | head -n 1) > merged_checkm.txt
for f in $(ls checkm/*.txt); do
    bn=$(basename $f)
    aid=${bn%%_*}
    awk -v aid=${aid}_ 'NR > 1 { print aid $0 }' ${f} >> merged_checkm.txt
done

python -c "from qp_pacbio.util import client_connect; qclient = client_connect('{{url}}'); qclient.update_job_step('{{qjid}}', 'Running checkm')"

checkm lineage_wf LCG LCG_checkm -x fna -t {{nprocs}} --tab_table -f LCG_checkm.txt --pplacer_threads 1
awk FNR!=1 LCG_checkm.txt >> merged_checkm.txt;

mv LCG/*.fna all_fna/

python -c "from qp_pacbio.util import client_connect; qclient = client_connect('{{url}}'); qclient.update_job_step('{{qjid}}', 'Running galah')"

galah cluster --checkm-tab-table merged_checkm.txt --genome-fasta-directory all_fna \
    -x fna --min-completeness 50 --max-contamination 10 --quality-formula dRep \
    -t {{nprocs}} --cluster-method fastani \
    --output-representative-fasta-directory-copy dereplicated \
    --output-cluster-definition dereplicated.txt \
    --ani {{percent_identity}} \
    --precluster-method finch

python -c "from qp_pacbio.util import client_connect; qclient = client_connect('{{url}}'); qclient.update_job_step('{{qjid}}', 'Running gtdbtk')"

# note that we are limiting pplacer_cpus to 3 to reduce the
# memory requirements
export GTDBTK_DATA_PATH="/scratch/qp-pacbio/gtdbtk_v226_db/"
gtdbtk classify_wf --genome_dir dereplicated \
    --out_dir dereplicated_gtdbtk --cpus {{nprocs}} \
    --pplacer_cpus 3 -x fna --skip_ani_screen

python -c "from qp_pacbio.util import client_connect; qclient = client_connect('{{url}}'); qclient.update_job_step('{{qjid}}', 'Running GToTree')"

ls ${PWD}/dereplicated/* > genomes.txt
GToTree -f genomes.txt -o phylogeny -j {{nprocs}} -H Bacteria -c {{GToTree_c}} -G {{GToTree_G}}

for f in $(ls {{output}}/*_sample_list.txt); do
    cat $f >> {{output}}/sample_list.txt
    echo "" >> {{output}}/sample_list.txt
done

python -c "from qp_pacbio.util import client_connect; qclient = client_connect('{{url}}'); qclient.update_job_step('{{qjid}}', 'Remapping')"
