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
cd {{output}}/

# The expected files here are:
# {aid}_prep_info_artifact.tsv : preparation for each artifact
# {aid}_folders.tsv : LCG/MAG folders to check for eaach artifact
# {aid}_{pid}_sample_list.txt : the sample/files list of the artifact, {pid} is the parent id

# Step 1: create processing folders, note that all files need to be gunzip
mkdir -p LCG MAG checkm

for folder in $(cat *_folders.tsv); do
    ffn=$(basename $folder);
    aid=${ffn%%_*};
    # the pattern of the folers is [folder]/results/[LCG|MAG]/[sample-id]/
    for f in $(ls ${f}/*/LCG/*/*fna.gz); do
        bn=$(basename ${f/.gz/});
        echo "gunzip -c ${f} > LCG/${aid}_${bn}"
    done
    for f in $(ls ${f}/*/MAG/*/*fna.gz); do
        bn=$(basename ${f/.gz/});
        echo "gunzip -c ${f} > MAG/${aid}_${bn}"
    done
    # the pattern for txts is [folder]/results/*checkm.txt.gz
    for f in $(ls ${f}/*/*checkm.txt.gz); do
        bn=$(basename ${f/.gz/});
        echo "zcat ${f} > checkm/${aid}_${bn}";
    done
done | parallel --halt now,fail=1 -j {{nprocs}}

# merging the checkm files
checkm_files=$(ls checkm/*.txt)
head -n 1 ${checkm_files[0]} > merged_checkm.txt
for f in checkm_files; do
    bn=$(basename $f)
    aid=${bn%%_*}
    awk FNR!=1 ${f} | awk "{ print ${aid}$0 }" >> merged_checkm.txt;
done

checkm lineage_wf LCG LCG_checkm -x fna -t {{nprocs}} --tab_table -f LCG_checkm.txt --pplacer_threads 1
awk FNR!=1 LCG_checkm.txt >> merged_checkm.txt;

galah cluster --checkm-tab-table merged_checkm.txt --genome-fasta-directory all_fna \
    -x fna --min-completeness 50 --max-contamination 10 --quality-formula dRep \
    -t {{nprocs}} --cluster-method fastani \
    --output-representative-fasta-directory-copy dereplicated
    --output-cluster-definition dereplicated.txt
    --ani {{percent_identity}} \
    --precluster-method finch
