#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p qiita
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -o {{output}}/merge/logs/%x-%A_%a.out
#SBATCH -e {{output}}/merge/logs/%x-%A_%a.err

source ~/.bashrc
set -e
{{conda_environment}}
cd {{output}}/
tax=/scratch/qp-woltka/WoLr2/WoLr2.tax
coords=/scratch/qp-woltka/WoLr2/WoLr2.coords
len_map=/scratch/qp-woltka/WoLr2/genomes/length.map
functional_dir=/scratch/qp-woltka/WoLr2/function/kegg/

mkdir -p {{output}}/coverages/

for f in `ls alignments/*.sam.xz`; do
    sn=`basename ${f/.sam.xz/}`;
    of={{output}}/bioms/${sn};
    mkdir -p ${of};
    echo "woltka classify -i ${f} -o ${of}/none.biom --no-demux --lineage ${tax} --rank none --outcov {{output}}/coverages/";
    echo "woltka classify -i ${f} -o ${of}/per-gene.biom --no-demux -c ${coords}";
done | parallel -j {{nprocs}}
wait

for f in `ls bioms/*/per-gene.biom`; do
    dn=`dirname ${f}`;
    sn=`basename ${dn}`;
    echo "woltka collapse -i ${f} -m ${functional_dir}/orf-to-ko.map.xz -o ${dn}/ko.biom; " \
        "woltka collapse -i ${dn}/ko.biom -m ${functional_dir}/ko-to-ec.map -o ${dn}/ec.biom; " \
        "woltka collapse -i ${dn}/ko.biom -m ${functional_dir}/ko-to-reaction.map -o ${dn}/reaction.biom; " \
        "woltka collapse -i ${dn}/reaction.biom -m ${functional_dir}/reaction-to-module.map -o ${dn}/module.biom; " \
        "woltka collapse -i ${dn}/module.biom -m ${functional_dir}/module-to-pathway.map -o ${dn}/pathway.biom;"
done | parallel -j {{nprocs}}
wait

biom_merge_pacbio --base {{output}}

find {{output}}/coverages/ -iname "*.cov" > {{output}}/cov_files.txt
micov consolidate --paths {{output}}/cov_files.txt --lengths ${len_map} --output {{output}}/coverages.tgz

cd alignments
tar -cvf ../alignments.tar *.sam.xz

finish_qp_pacbio {{url}} {{qjid}} {{output}}
