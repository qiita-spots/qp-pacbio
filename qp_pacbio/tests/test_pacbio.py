# -----------------------------------------------------------------------------
# Copyright (c) 2020--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from json import dumps
from os import remove
from os.path import exists, isdir, join
from shutil import copyfile, rmtree
from tempfile import mkdtemp
from unittest import main

from qiita_client.testing import PluginTestCase

from qp_pacbio import plugin
from qp_pacbio.qp_pacbio import (
    BASEPATH,
    CONDA_ENVIRONMENT,
    generate_minimap2_processing,
    generate_sample_list,
    generate_syndna_processing,
    pacbio_generate_templates,
)

# Exact expected Step-1 script matching your template.
# Escape ${...} -> ${{...}} and awk braces -> '{{print $1}}' etc.
STEP_1_EXP = (
    "#!/bin/bash\n"
    "#SBATCH -J s1-{job_id}\n"
    f"#SBATCH -p qiita\n"
    f"#SBATCH -N 1\n"
    f"#SBATCH -n 16\n"
    f"#SBATCH --time 1-00:00:00\n"
    f"#SBATCH --mem 200G\n"
    "#SBATCH -o {out_dir}/step-1/logs/%x-%A_%a.out\n"
    "#SBATCH -e {out_dir}/step-1/logs/%x-%A_%a.err\n"
    "#SBATCH --array 1-{njobs}%16\n"
    "\n"
    "source ~/.bashrc\n"
    "set -e\n"
    f"{CONDA_ENVIRONMENT}\n"
    "cd {out_dir}/step-1\n"
    "\n"
    "step=${{SLURM_ARRAY_TASK_ID}}\n"
    "input=$(head -n $step {out_dir}/sample_list.txt | tail -n 1)\n"
    "\n"
    "sample_name=`echo $input | awk '{{print $1}}'`\n"
    "filename=`echo $input | awk '{{print $2}}'`\n"
    "\n"
    "fn=`basename ${{filename}}`\n"
    "\n"
    "# updating the GUI when task 1 runs\n"
    'if [[ "$step" == "1" ]]; then\n'
    '    python -c "from qp_pacbio.util import client_connect; '
    "qclient = client_connect('http://test.test.org'); "
    "qclient.update_job_step('my-job-id', 'Running step 1: "
    "${{SLURM_ARRAY_JOB_ID}}')\"\n"
    "fi\n"
    "\n"
    f"hifiasm_meta -t 16 -o "
    "{out_dir}/step-1/${{sample_name}} ${{filename}}\n"
    "touch {out_dir}/step-1/completed_${{SLURM_ARRAY_TASK_ID}}.log"
)


class PacBioTests(PluginTestCase):
    def setUp(self):
        plugin("https://localhost:21174", "register", "ignored")

        self._clean_up_files = []

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def _insert_data(self):
        prep_info_dict = {
            "SKB8.640193": {"run_prefix": "S22205_S104"},
            "SKD8.640184": {"run_prefix": "S22282_S102"},
        }
        data = {
            "prep_info": dumps(prep_info_dict),
            # magic #1 = testing study
            "study": 1,
            "data_type": "Metagenomic",
        }
        pid = self.qclient.post("/apitest/prep_template/", data=data)["prep"]

        # inserting artifacts
        in_dir = mkdtemp()
        self._clean_up_files.append(in_dir)

        fp1 = join(in_dir, "S22205_S104_L001_R1_001.fastq.gz")
        fp2 = join(in_dir, "S22282_S102_L001_R1_001.fastq.gz")
        fp_summary = join(in_dir, "summary.html")
        source_dir = "qp_pacbio/tests/support_files/"
        copyfile(f"{source_dir}/S22205_S104_L001_R1_001.fastq.gz", fp1)
        copyfile(f"{source_dir}/S22282_S102_L001_R1_001.fastq.gz", fp2)
        copyfile(f"{source_dir}/summary.html", fp_summary)

        data = {
            "filepaths": dumps(
                [
                    (fp1, "raw_forward_seqs"),
                    (fp2, "raw_forward_seqs"),
                    (fp_summary, "html_summary"),
                ]
            ),
            "type": "per_sample_FASTQ",
            "name": "Test artifact",
            "prep": pid,
        }
        aid = self.qclient.post("/apitest/artifact/", data=data)["artifact"]

        return aid


class PacProcessingTests(PacBioTests):
    def test_pacbio_processing(self):
        artifact_id = self._insert_data()
        job_id = "my-job-id"
        out_dir = mkdtemp()
        result_fp = join(out_dir, "result")
        self._clean_up_files.append(out_dir)

        njobs = generate_sample_list(self.qclient, artifact_id, out_dir)
        pacbio_generate_templates(
            out_dir, job_id, njobs, result_fp, "http://test.test.org"
        )

        # testing file creation, just number of lines and header
        with open(f"{out_dir}/sample_list.txt", "r") as f:
            obs_lines = f.readlines()
        self.assertEqual(2, len(obs_lines))
        njobs = len(obs_lines)

        # testing step-1
        with open(f"{out_dir}/step-1/step-1.slurm", "r") as f:
            # removing \n
            obs_lines = [ln.replace("\n", "") for ln in f.readlines()]

        self.assertCountEqual(
            STEP_1_EXP.format(
                out_dir=out_dir,
                job_id=job_id,
                njobs=njobs,
            ).split("\n"),
            obs_lines,
        )


class PacWoltkaProfilingTests(PacBioTests):
    def test_pacbio_profiling(self):
        params = {"artifact": int(self._insert_data())}
        job_id = "my-job-id"
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        url = "https://test.test.edu/"
        # this should fail cause we don't have valid data
        main_fp, merge_fp = generate_minimap2_processing(
            self.qclient, job_id, out_dir, params, url
        )
        with open(main_fp, "r") as f:
            obs_main = f.readlines()
        with open(merge_fp, "r") as f:
            obs_merge = f.readlines()

        exp_main = [
            "#!/bin/bash\n",
            "#SBATCH -J m2_my-job-id\n",
            "#SBATCH -p qiita\n",
            "#SBATCH -N 1\n",
            "#SBATCH -n 16\n",
            "#SBATCH --time 36000\n",
            "#SBATCH --mem 60G\n",
            f"#SBATCH -o {out_dir}/minimap2/logs/%x-%A_%a.out\n",
            f"#SBATCH -e {out_dir}/minimap2/logs/%x-%A_%a.err\n",
            "#SBATCH --array 1-2%16\n",
            "\n",
            "source ~/.bashrc\n",
            "set -e\n",
            f"{CONDA_ENVIRONMENT}\n",
            f"mkdir -p {out_dir}/alignment {out_dir}/filtered-alignment\n",
            f"cd {out_dir}/\n",
            "db=/scratch/qp-pacbio/minimap2/WoLr2/WoLr2.map-hifi.mmi\n",
            "\n",
            "step=${SLURM_ARRAY_TASK_ID}\n",
            f"input=$(head -n $step {out_dir}/sample_list.txt | tail -n 1)\n",
            "\n",
            "sample_name=`echo $input | awk '{print $1}'`\n",
            "filename=`echo $input | awk '{print $2}'`\n",
            "\n",
            "fn=`basename ${filename}`\n",
            "\n",
            "minimap2 -x map-hifi -t 16 -a \\\n",
            "       --secondary=no --MD --eqx ${db} \\\n",
            "       ${filename} | \\\n",
            '   awk \'BEGIN { FS=OFS="\\t" } /^@/ { print; next } '
            '{ $10="*"; $11="*" } 1\' | grep -v "^@" | sort -k 1 | \\\n',
            f"   xz -1 -T1 > {out_dir}/alignment/${{sample_name}}.sam.xz\n",
            "\n",
            f"xcat {out_dir}/alignment/${{sample_name}}.sam.xz | \\\n",
            f"    /home/mcdonadt/duckdb-miint/duckdb-2025.11.20 -f {BASEPATH}/data/sql//qcov-seqident-filter.sql > {out_dir}/filtered-alignment/${{sample_name}}.sam.gz",
        ]
        self.assertEqual(obs_main, exp_main)

        db_path = "/scratch/qp-woltka/WoLr2"
        exp_merge = [
            "#!/bin/bash\n",
            "#SBATCH -J me_my-job-id\n",
            "#SBATCH -p qiita\n",
            "#SBATCH -N 1\n",
            "#SBATCH -n 16\n",
            "#SBATCH --time 1-00:00:00\n",
            "#SBATCH --mem 120G\n",
            f"#SBATCH -o {out_dir}/merge/logs/%x-%A_%a.out\n",
            f"#SBATCH -e {out_dir}/merge/logs/%x-%A_%a.err\n",
            "\n",
            "source ~/.bashrc\n",
            "set -e\n",
            f"{CONDA_ENVIRONMENT}\n",
            f"cd {out_dir}/\n",
            f"tax={db_path}/WoLr2.tax\n",
            f"coords={db_path}/WoLr2.coords\n",
            f"len_map={db_path}/genomes/length.map\n",
            f"functional_dir={db_path}/function/kegg/\n",
            "\n",
            f"mkdir -p {out_dir}/coverages/\n",
            "\n",
            "for f in `ls filtered-alignment/*.sam.xz`; do\n",
            "    sn=`basename ${f/.sam.xz/}`;\n",
            f"    of={out_dir}/bioms/${{sn}};\n",
            "    mkdir -p ${of};\n",
            '    echo "woltka classify -i ${f} -o ${of}/none.biom '
            "--no-demux --lineage ${tax} "
            f'--rank none --outcov {out_dir}/coverages/";\n',
            '    echo "woltka classify -i ${f} -o ${of}/per-gene.biom '
            '--no-demux -c ${coords}";\n',
            "done | parallel --halt now,fail=1 -j 16\n",
            "wait\n",
            "\n",
            "for f in `ls bioms/*/per-gene.biom`; do\n",
            "    dn=`dirname ${f}`;\n",
            "    sn=`basename ${dn}`;\n",
            '    echo "woltka collapse -i ${f} -m ${functional_dir}/'
            'orf-to-ko.map.xz -o ${dn}/ko.biom; " \\\n',
            '        "woltka collapse -i ${dn}/ko.biom -m ${functional_dir}/'
            'ko-to-ec.map -o ${dn}/ec.biom; " \\\n',
            '        "woltka collapse -i ${dn}/ko.biom -m ${functional_dir}/'
            'ko-to-reaction.map -o ${dn}/reaction.biom; " \\\n',
            '        "woltka collapse -i ${dn}/reaction.biom -m '
            "${functional_dir}/reaction-to-module.map -o "
            '${dn}/module.biom; " \\\n',
            '        "woltka collapse -i ${dn}/module.biom -m '
            "${functional_dir}/module-to-pathway.map "
            '-o ${dn}/pathway.biom;"\n',
            "done | parallel --halt now,fail=1 -j 16\n",
            "wait\n",
            "\n",
            f"biom_merge_pacbio --base {out_dir} --merge-type woltka\n",
            "\n",
            (
                f'find {out_dir}/coverages/ -iname "*.cov" '
                + f"> {out_dir}/cov_files.txt\n"
            ),
            f"micov consolidate --paths {out_dir}/cov_files.txt "
            f"--lengths ${{len_map}} --output {out_dir}/coverages.tgz\n",
            "\n",
            "cd alignment\n",
            "tar -cvf ../alignment.tar *.sam.xz\n",
            "\n",
            f"finish_qp_pacbio {url} {job_id} {out_dir}",
        ]
        self.assertEqual(obs_merge, exp_merge)


class PacWoltkaSynDNATests(PacBioTests):
    def test_syndna(self):
        params = {"artifact": int(self._insert_data())}
        job_id = "my-job-id"
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        url = "https://test.test.edu/"
        # this should fail cause we don't have valid data
        main_fp, finish_fp = generate_syndna_processing(
            self.qclient, job_id, out_dir, params, url
        )
        with open(main_fp, "r") as f:
            # ignoring comments
            obs_main = [line for line in f.readlines() if not line.startswith("# ")]
        with open(finish_fp, "r") as f:
            # ignoring comments
            obs_finish = [line for line in f.readlines() if not line.startswith("# ")]

        exp_main = [
            "#!/bin/bash\n",
            "#SBATCH -J sd_my-job-id\n",
            "#SBATCH -p qiita\n",
            "#SBATCH -N 1\n",
            "#SBATCH -n 16\n",
            "#SBATCH --time 14400\n",
            "#SBATCH --mem 20G\n",
            f"#SBATCH -o {out_dir}/syndna/logs/%x-%A_%a.out\n",
            f"#SBATCH -e {out_dir}/syndna/logs/%x-%A_%a.err\n",
            "#SBATCH --array 1-2%16\n",
            "\n",
            "source ~/.bashrc\n",
            "set -e\n",
            f"{CONDA_ENVIRONMENT}\n",
            f"out_folder={out_dir}/syndna\n",
            "mkdir -p ${out_folder}\n",
            "cd ${out_folder}\n",
            "db_folder=/scratch/qp-pacbio/minimap2/syndna/\n",
            "\n",
            "step=${SLURM_ARRAY_TASK_ID}\n",
            f"input=$(head -n $step {out_dir}/sample_list.txt | tail -n 1)\n",
            "sample_name=`echo $input | awk '{print $1}'`\n",
            "filename=`echo $input | awk '{print $2}'`\n",
            "fn=`basename ${filename}`\n",
            "\n",
            "mkdir -p ${out_folder}/filtered/\n",
            "\n",
            "sn_folder=${out_folder}/bioms/${sample_name}\n",
            "mkdir -p ${sn_folder}\n",
            "\n",
            "txt=${sn_folder}/${sample_name}.txt\n",
            "tsv=${txt/.txt/.tsv}\n",
            "coverm contig --single $filename --reference ${db_folder}/All_synDNA_inserts.fasta --mapper minimap2-hifi \\\n",
            "    --min-read-percent-identity 0.95 --min-read-aligned-percent 0.0 -m mean count --threads 16 \\\n",
            "    --output-file ${txt}\n",
            "\n",
            "awk 'BEGIN {FS=OFS=\"\\t\"}; {print $1,$3}' ${txt} | \\\n",
            "    sed 's/Contig/\#OTU ID/' | sed  's/All_synDNA_inserts.fasta\///' | \\\n",
            "    sed  's/ Read Count//' | sed \"s/${fn}/${sample_name}/\" > ${tsv}\n",
            "\n",
            "counts=`tail -n +2 ${tsv} | awk '{sum += $NF} END {print sum}'`\n",
            'if [[ "$counts" == "0" ]]; then\n',
            f"    echo ${{sample_name}} > {out_dir}/failed_${{SLURM_ARRAY_TASK_ID}}.log\n",
            "    exit 0\n",
            "fi\n",
            "\n",
            "biom convert -i ${tsv} -o ${sn_folder}/syndna.biom --to-hdf5\n",
            "\n",
            "minimap2 -x map-hifi -t 16 -a --MD --eqx -o ${out_folder}/${sample_name}_plasmid.sam ${db_folder}/AllsynDNA_plasmids_FASTA_ReIndexed_FINAL.fasta $filename\n",
            "samtools view -F 4 -@ 16 ${out_folder}/${sample_name}_plasmid.sam | awk '{print $1}' | sort -u > ${out_folder}/${sample_name}_plasmid_mapped.txt\n",
            "seqkit grep -v -f ${out_folder}/${sample_name}_plasmid_mapped.txt $filename > ${out_folder}/${sample_name}_no_plasmid.fastq\n",
            "\n",
            "minimap2 -x map-hifi -t 16 -a --MD --eqx -o ${out_folder}/${sample_name}_GCF_000184185.sam ${db_folder}/GCF_000184185.1_ASM18418v1_genomic_chroso.fna ${out_folder}/${sample_name}_no_plasmid.fastq\n",
            "samtools view -bS -@ 8.0 ${out_folder}/${sample_name}_no_plasmid.fastq | samtools sort -@ 8.0 -O bam -o ${out_folder}/${sample_name}_GCF_000184185_sorted.bam\n",
            "coverm filter --bam-files ${out_folder}/${sample_name}_GCF_000184185_sorted.bam --min-read-percent-identity 99.9 --min-read-aligned-percent 95 --threads 16 -o ${out_folder}/${sample_name}_GCF_000184185.bam\n",
            "samtools view -O SAM -o ${out_folder}/${sample_name}_no_GCF_000184185_sorted.sam ${out_folder}/${sample_name}_GCF_000184185.bam\n",
            "awk '{print $1}' ${out_folder}/${sample_name}_no_GCF_000184185_sorted.sam > ${out_folder}/${sample_name}_GCF_000184185_reads_filtered.txt\n",
            "seqkit grep -v -f ${out_folder}/${sample_name}_GCF_000184185_reads_filtered.txt ${out_folder}/${sample_name}_no_plasmid.fastq | gzip > ${out_folder}/filtered/${fn}\n",
            "\n",
            f"touch {out_dir}/completed_${{SLURM_ARRAY_TASK_ID}}.log",
        ]

        self.assertEqual(obs_main, exp_main)

        exp_merge = [
            "#!/bin/bash\n",
            "#SBATCH -J me_my-job-id\n",
            "#SBATCH -p qiita\n",
            "#SBATCH -N 1\n",
            "#SBATCH -n 1\n",
            "#SBATCH --time 14400\n",
            "#SBATCH --mem 4G\n",
            f"#SBATCH -o {out_dir}/finish/logs/%x-%A_%a.out\n",
            f"#SBATCH -e {out_dir}/finish/logs/%x-%A_%a.err\n",
            "\n",
            "source ~/.bashrc\n",
            "set -e\n",
            f"{CONDA_ENVIRONMENT}\n",
            f"cd {out_dir}/\n",
            "\n",
            f"biom_merge_pacbio --base {out_dir}/syndna --merge-type syndna\n",
            "\n",
            "\n",
            "\n",
            f"finish_qp_pacbio https://test.test.edu/ my-job-id {out_dir}",
        ]
        self.assertEqual(obs_finish, exp_merge)


if __name__ == "__main__":
    main()
