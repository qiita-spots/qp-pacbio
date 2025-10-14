#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/step-0/logs/%x-%A_%a.out
#SBATCH -e {{output}}/step-0/logs/%x-%A_%a.out

source ~/.bashrc
set -e
conda activate {{conda_environment}}
cd {{output}}/
tax=/projects/wol/qiyun/wol2/databases/minimap2/WoLr2.tax
coords=/projects/wol/qiyun/wol2/databases/minimap2/WoLr2.coords
len_map=/projects/wol/qiyun/wol2/databases/minimap2/WoLr2/length.map
functional_dir=/projects/wol/qiyun/wol2/databases/minimap2/WoLr2/

mkdir -p {{output}}/coverages/

for f in `ls alignments/*.sam.xz`; do
    sn=`basename ${f/.sam.xz/}`;
    of={{output}}/bioms/${sn};
    mkdir -p ${of};
    echo "woltka classify -i ${f} -o ${of}/none.biom --no-demux --lineage ${tax} --rank none --outcov {{output}}/coverages/";
    echo "woltka classify -i ${f} -o ${of}/per-gene.biom --no-demux -c ${coords}";
done | parallel -j {{node_count}}
wait

for f in `ls bioms/*/per-gene.biom`; do
    dn=`dirname ${f}`;
    sn=`basename ${sn}`;
    echo "woltka collapse -i ${f} -m ${functional_dir}/orf-to-ko.map.xz -o ${dn}/ko.biom; " \
        "woltka collapse -i ${dn}/ko.biom -m ${functional_dir}/ko-to-ec.map -o ${dn}/ec.biom; " \
        "woltka collapse -i ${dn}/ko.biom -m ${functional_dir}/ko-to-reaction.map -o ${dn}/reaction.biom; " \
        "woltka collapse -i ${dn}/reaction.biom -m ${functional_dir}/reaction-to-module.map -o ${dn}/module.biom; " \
        "woltka collapse -i ${dn}/module.biom -m ${functional_dir}/module-to-pathway.map -o ${dn}/pathway.biom;"
done | parallel -j {{node_count}}
wait

# MISSING:
# merge bioms!

find {{output}}/coverages/ -iname "*.cov" > {{output}}/cov_files.txt
micov consolidate --paths {{output}}/cov_files.txt --lengths ${len_map} --output {{output}}/coverages.tgz

cd alignments
tar -cvf ../alignments.tar *.sam.xz
