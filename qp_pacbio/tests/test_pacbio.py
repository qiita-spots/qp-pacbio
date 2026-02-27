# -----------------------------------------------------------------------------
# Copyright (c) 2020--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from json import dumps
from os import makedirs, remove
from os.path import dirname, exists, isdir, join
from pathlib import Path
from shutil import copyfile, rmtree
from tempfile import mkdtemp
from unittest import main

from qiita_client.testing import PluginTestCase

from qp_pacbio import plugin
from qp_pacbio.qp_pacbio import (
    CONDA_ENVIRONMENT,
    RESOURCES,
    generate_feature_table_scripts,
    generate_minimap2_processing,
    generate_pacbio_adapter_removal,
    generate_sample_list,
    generate_syndna_processing,
    pacbio_generate_templates,
)
from qp_pacbio.util import plugin_details

# Exact expected Step-1 script matching your template.
# Escape ${...} -> ${{...}} and awk braces -> '{{print $1}}' etc.
STEP_1_EXP = (
    "#!/bin/bash\n"
    "#SBATCH -J s1-{job_id}\n"
    f"#SBATCH -p qiita\n"
    f"#SBATCH -N 1\n"
    f"#SBATCH -n 16\n"
    f"#SBATCH --time 2-00:00:00\n"
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
        params = {
            "artifact": int(self._insert_data()),
            "Database": "WoLr2",
            "percent-identity": 0.9,
            "percent-coverage": 0.9,
        }
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
        with open(f"{out_dir}/qcov-seqident-filter.sql", "r") as f:
            obs_sql = f.readlines()

        resources = RESOURCES["Woltka v0.1.7 with cov and id filter"]
        reference = resources["minimap2"]["databases"]["WoLr2"]
        reference_db = reference["reference_db"]
        reference_tax = reference["reference_tax"]
        reference_coords = reference["reference_coords"]
        reference_len_map = reference["reference_len_map"]
        reference_functional_dir = reference["reference_functional_dir"]
        miint_path = resources["minimap2"]["miint_path"]

        exp_sql = [
            "CREATE TEMP TABLE lengths AS\n",
            "    SELECT *\n",
            f'        FROM read_csv(\'{reference_len_map}\', header=false, delim=\'\\t\', columns = {{"read_id": "VARCHAR", "length": "BIGINT"}});\n',
            "\n",
            "COPY (SELECT *\n",
            "      FROM read_alignments('/dev/stdin', reference_lengths='lengths')\n",
            "          where alignment_seq_identity(cigar, tag_nm, tag_md, 'blast') > 0.9\n",
            "              and alignment_query_coverage(cigar) > 0.9\n",
            "     ) TO '/dev/stdout' (FORMAT SAM, INCLUDE_HEADER false, COMPRESSION gzip);",
        ]
        self.assertEqual(obs_sql, exp_sql)

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
            f"db={reference_db}\n",
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
            f"xzcat {out_dir}/alignment/${{sample_name}}.sam.xz | \\\n",
            f"    {miint_path} -f {out_dir}/qcov-seqident-filter.sql > {out_dir}/filtered-alignment/${{sample_name}}.sam.gz",
        ]
        self.assertEqual(obs_main, exp_main)

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
            f"tax={reference_tax}\n",
            f"coords={reference_coords}\n",
            f"len_map={reference_len_map}\n",
            f"functional_dir={reference_functional_dir}\n",
            "\n",
            f"mkdir -p {out_dir}/coverages/\n",
            "\n",
            "for f in `ls filtered-alignment/*.sam.gz`; do\n",
            '    if [[ $(gunzip -c "$f" | head -c1 | wc -c) == "0" ]]; then\n',
            "        continue\n",
            "    fi\n",
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
            "db_folder=/scratch/qp-pacbio/syndna/\n",
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
            "samtools view -bS -@ 8.0 ${out_folder}/${sample_name}_GCF_000184185.sam | samtools sort -@ 8.0 -O bam -o ${out_folder}/${sample_name}_GCF_000184185_sorted.bam\n",
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
            f"finish_qp_pacbio https://test.test.edu/ my-job-id {out_dir}",
        ]
        self.assertEqual(obs_finish, exp_merge)


class PacAdapterRmTests(PacBioTests):
    def test_pacbio_adapter_removal(self):
        params = {
            "artifact": int(self._insert_data()),
            # this adapter selection will test all options: file and direct adapter value
            "adapter_sets": "AAGCAGTGGTATCAACGCAGAGTACT,twist_adapters_231010.fasta.gz",
            "hifi-preset": "SYMMETRIC",
            "peek-guess": True,
        }
        job_id = "my-job-id"
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        url = "https://test.test.edu/"
        # this should fail cause we don't have valid data
        main_fp, finish_fp = generate_pacbio_adapter_removal(
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
            "#SBATCH -n 1\n",
            "#SBATCH --time 57600\n",
            "#SBATCH --mem 100G\n",
            f"#SBATCH -o {out_dir}/processing/logs/%x-%A_%a.out\n",
            f"#SBATCH -e {out_dir}/processing/logs/%x-%A_%a.err\n",
            "#SBATCH --array 1-2%16\n",
            "\n",
            "source ~/.bashrc\n",
            "set -e\n",
            f"{CONDA_ENVIRONMENT}\n",
            f"out_folder={out_dir}/processing/\n",
            "mkdir -p ${out_folder}/lima ${out_folder}/final\n",
            "\n",
            "step=${SLURM_ARRAY_TASK_ID}\n",
            f"input=$(head -n $step {out_dir}/sample_list.txt | tail -n 1)\n",
            "sample_name=`echo $input | awk '{print $1}'`\n",
            "filename=`echo $input | awk '{print $2}'`\n",
            "fn=`basename ${filename}`\n",
            "fout=${out_folder}/lima/${fn/.fastq.gz/}\n",
            "final=${out_folder}/final/${fn/.fastq.gz/}\n",
            "\n",
            "lima_error_handler () {\n",
            '    if grep -qF "Could not find matching barcodes" ${fout}.lima.log; then\n',
            f"        touch {out_dir}/completed_${{SLURM_ARRAY_TASK_ID}}.log\n",
            '         cp "${filename}" "${final}.fastq.gz"\n',
            "        exit 0\n",
            "    else\n",
            f'        echo "${{sample_name}}" > {out_dir}/failed.log\n',
            "        exit 0\n",
            "    fi\n",
            "}\n",
            "\n",
            "trap lima_error_handler ERR\n",
            "\n",
            f'lima "${{filename}}" {out_dir}/adapter.fasta "${{fout}}.fastq.gz" --hifi-preset SYMMETRIC --peek-guess > ${{fout}}.lima.log 2>&1\n',
            "\n",
            "column=$(head -1 \"${fout}.lima.counts\" | tr '\\t' '\\n' | grep -En \"Counts\" | cut -d: -f1)\n",
            'counts=$(cut -f "${column}" "${fout}.lima.counts" | tail -n 1)\n',
            "\n",
            'if [[ "$counts" -gt "0" ]]; then\n',
            '    pbmarkdup --rmdup "${fout}.fastq.gz" "${final}.fastq"\n',
            "fi\n",
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
            "#SBATCH --time 3600\n",
            "#SBATCH --mem 4G\n",
            f"#SBATCH -o {out_dir}/finish/logs/%x-%A_%a.out\n",
            f"#SBATCH -e {out_dir}/finish/logs/%x-%A_%a.err\n",
            "\n",
            "source ~/.bashrc\n",
            "set -e\n",
            f"{CONDA_ENVIRONMENT}\n",
            f"cd {out_dir}/\n",
            "\n",
            f"finish_qp_pacbio https://test.test.edu/ my-job-id {out_dir}",
        ]
        self.assertEqual(obs_finish, exp_merge)


class PacBioFeatureTableTests(PacBioTests):
    def _helper_woltka_bowtie(self):
        aids = []
        in_dir = mkdtemp()
        self._clean_up_files.append(in_dir)
        for i in range(2):
            parent_id = self._insert_data()

            # now that we have a parent we can create a job
            # to simulate that we process the parent to create
            # the inputs for this tests

            data = {
                "user": "demo@microbio.me",
                "command": dumps(
                    [
                        plugin_details["name"],
                        plugin_details["version"],
                        "PacBio processing",
                    ]
                ),
                "status": "running",
                "parameters": dumps({"artifact": parent_id}),
            }
            job_id = self.qclient.post("/apitest/processing_job/", data=data)["job"]

            # creating folder structures to then insert artifacts
            self._clean_up_files.append(in_dir)
            results = join(in_dir, "tmp", f"results_{i}")
            for folder in ["small_LCG", "MAG", "LCG"]:
                for sample in ["1.SKB8.640193", "1.SKD8.640184"]:
                    makedirs(f"{results}/{sample}/{folder}")
            Path(f"{results}/1.SKB8.640193/1.SKD8.640184.checkm.txt.gz").touch()
            Path(f"{results}/1.SKD8.640184/1.SKD8.640184.checkm.txt.gz").touch()
            Path(f"{results}/1.SKB8.640193/1.SKD8.640184.noLCG.fna.gz").touch()
            Path(f"{results}/1.SKD8.640184/1.SKD8.640184.noLCG.fna.gz").touch()

            Path(
                f"{results}/1.SKB8.640193/LCG/1.SKD8.640184_MaxBin_bin.1.fna.gz"
            ).touch()
            Path(
                f"{results}/1.SKB8.640193/LCG/1.SKD8.640184_MaxBin_bin.0.fna.gz"
            ).touch()
            Path(
                f"{results}/1.SKD8.640184/LCG/1.SKD8.640184_MaxBin_bin.1.fna.gz"
            ).touch()
            Path(
                f"{results}/1.SKD8.640184/LCG/1.SKD8.640184_MaxBin_bin.0.fna.gz"
            ).touch()

            Path(
                f"{results}/1.SKB8.640193/MAG/1.SKB8.640193.s37.ctg000038c.fna.gz"
            ).touch()
            Path(
                f"{results}/1.SKB8.640193/MAG/1.SKB8.640193.s40.ctg000042c.fna.gz"
            ).touch()
            Path(
                f"{results}/1.SKD8.640184/MAG/1.SKD8.640184.s37.ctg000038c.fna.gz"
            ).touch()
            Path(
                f"{results}/1.SKD8.640184/MAG/1.SKD8.640184.s40.ctg000042c.fna.gz"
            ).touch()

            data = {
                "filepaths": dumps(
                    [
                        (results, "directory"),
                    ]
                ),
                "type": "job-output-folder",
                "name": "Test artifact",
                "parents": dumps([parent_id]),
                "job_id": job_id,
            }
            aid = self.qclient.post("/apitest/artifact/", data=data)["artifact"]
            aids.append(aid)
            # retriving final location of files to delete them and avoid future errors
            fps, _ = self.qclient.artifact_and_preparation_files(aid)
            self._clean_up_files.append(dirname(fps["directory"][0]))

        data = {
            "user": "demo@microbio.me",
            "command": dumps(
                [
                    plugin_details["name"],
                    plugin_details["version"],
                    "Feature Table from LCG/MAG",
                ]
            ),
            "status": "running",
            # 2 is an existing analysis in the Qiita development database
            "parameters": dumps({"analysis": 2, "artifacts": aids}),
        }
        job_id = self.qclient.post("/apitest/processing_job/", data=data)["job"]

        return job_id, aids

    def test_pacbio_feature_table(self):
        job_id, aids = self._helper_woltka_bowtie()
        job_info = self.qclient.get_job_info(job_id)
        parameters = job_info["parameters"]

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        url = "https://test.test.edu/"
        merge_fp, remap_fp, finish_fp = generate_feature_table_scripts(
            self.qclient, job_id, out_dir, parameters, url
        )

        with open(merge_fp, "r") as f:
            obs_merge = [line for line in f.readlines() if not line.startswith("# ")]
        with open(remap_fp, "r") as f:
            obs_remap = [line for line in f.readlines() if not line.startswith("# ")]
        with open(finish_fp, "r") as f:
            obs_finish = [line for line in f.readlines() if not line.startswith("# ")]

        exp_merge = [
            "#!/bin/bash\n",
            f"#SBATCH -J m_{job_id}\n",
            "#SBATCH -p qiita\n",
            "#SBATCH -N 1\n",
            "#SBATCH -n 32\n",
            "#SBATCH --time 54000\n",
            "#SBATCH --mem 180G\n",
            f"#SBATCH -o {out_dir}/merge/logs/%x-%A.out\n",
            f"#SBATCH -e {out_dir}/merge/logs/%x-%A.err\n",
            "\n",
            "source ~/.bashrc\n",
            "set -e\n",
            f"{CONDA_ENVIRONMENT}\n",
            f"cd {out_dir}/merge\n",
            "\n",
            f"python -c \"from qp_pacbio.util import client_connect; qclient = client_connect('{url}'); qclient.update_job_step('{job_id}', 'Merging LCG/MAG files')\"\n",
            "\n",
            "\n",
            "mkdir -p LCG all_fna checkm\n",
            "\n",
            "for file in $(ls ../*_folders.tsv); do\n",
            "    for folder in $(cat ${file}); do\n",
            "        ffn=$(basename $folder);\n",
            "        aid=${ffn%%_*};\n",
            "        # the pattern of the folers is [folder]/results/[LCG|MAG]/[sample-id]/\n",
            "        for f in $(ls ${folder}/*/*/LCG/*fna.gz 2> /dev/null || true); do\n",
            "            bn=$(basename ${f/.gz/});\n",
            '            echo "gunzip -c ${f} > LCG/${aid}_${bn}"\n',
            "        done\n",
            "        for f in $(ls ${folder}/*/*/MAG/*fna.gz 2> /dev/null || true); do\n",
            "            bn=$(basename ${f/.gz/});\n",
            '            echo "gunzip -c ${f} > all_fna/${aid}_${bn}"\n',
            "        done\n",
            "        # the pattern for txts is [folder]/results/*checkm.txt.gz\n",
            "        for f in $(ls ${folder}/*/*/*checkm.txt.gz); do\n",
            "            bn=$(basename ${f/.gz/});\n",
            '            echo "zcat ${f} > checkm/${aid}_${bn}";\n',
            "        done\n",
            "    done\n",
            "done | parallel --halt now,fail=1 -j 32\n",
            "\n",
            "head -n 1 $(ls checkm/*.txt | head -n 1) > merged_checkm.txt\n",
            "for f in $(ls checkm/*.txt); do\n",
            "    bn=$(basename $f)\n",
            "    aid=${bn%%_*}\n",
            "    awk -v aid=${aid}_ 'NR > 1 { print aid $0 }' ${f} >> merged_checkm.txt\n",
            "done\n",
            "\n",
            f"python -c \"from qp_pacbio.util import client_connect; qclient = client_connect('{url}'); qclient.update_job_step('{job_id}', 'Running checkm')\"\n",
            "\n",
            'if [ -n "$(ls -A "LCG")" ]; then\n',
            "    checkm lineage_wf LCG LCG_checkm -x fna -t 32 --tab_table -f LCG_checkm.txt --pplacer_threads 1\n",
            "    awk FNR!=1 LCG_checkm.txt >> merged_checkm.txt;\n",
            "    mv LCG/*.fna all_fna/\n",
            "fi\n",
            "\n",
            f"python -c \"from qp_pacbio.util import client_connect; qclient = client_connect('{url}'); qclient.update_job_step('{job_id}', 'Running galah')\"\n",
            "\n",
            "galah cluster --checkm-tab-table merged_checkm.txt --genome-fasta-directory all_fna \\\n",
            "    -x fna --min-completeness 50 --max-contamination 10 --quality-formula dRep \\\n",
            "    -t 32 --cluster-method fastani \\\n",
            "    --output-representative-fasta-directory-copy dereplicated \\\n",
            "    --output-cluster-definition dereplicated.txt \\\n",
            "    --ani 0.995 \\\n",
            "    --precluster-method finch\n",
            "\n",
            f"python -c \"from qp_pacbio.util import client_connect; qclient = client_connect('{url}'); qclient.update_job_step('{job_id}', 'Running gtdbtk')\"\n",
            "\n",
            'export GTDBTK_DATA_PATH="/scratch/qp-pacbio/gtdbtk_v226_db/"\n',
            "gtdbtk classify_wf --genome_dir dereplicated \\\n",
            "    --out_dir dereplicated_gtdbtk --cpus 32 \\\n",
            "    --pplacer_cpus 3 -x fna --skip_ani_screen\n",
            "\n",
            f"python -c \"from qp_pacbio.util import client_connect; qclient = client_connect('{url}'); qclient.update_job_step('{job_id}', 'Running GToTree')\"\n",
            "\n",
            "ls ${PWD}/dereplicated/* > genomes.txt\n",
            "GToTree -f genomes.txt -o phylogeny -j 32 -H Bacteria -c 0.4 -G 0.4\n",
            "\n",
            f"for f in $(ls {out_dir}/*_sample_list.txt); do\n",
            "    fn=$(basename $f)\n",
            "    aid=${fn%%_*}\n",
            f'    sed "s/$/\\t${{aid}}/" ${{f}} >> {out_dir}/sample_list.txt\n',
            f'    echo "" >> {out_dir}/sample_list.txt\n',
            "done\n",
            "\n",
            f"python -c \"from qp_pacbio.util import client_connect; qclient = client_connect('{url}'); qclient.update_job_step('{job_id}', 'Remapping')\"",
        ]

        self.assertEqual(obs_merge, exp_merge)

        # checking the expected files and number of lines
        for aid in aids:
            with open(f"{out_dir}/{aid}_prep_info_artifact.tsv", "r") as fp:
                self.assertEqual(3, len(fp.readlines()))
            with open(f"{out_dir}/{aid}_folders.tsv", "r") as fp:
                self.assertEqual(1, len(fp.readlines()))
            # aid - 1, is a safe assumption as the child is created just after the parent
            with open(f"{out_dir}/{aid}_{aid - 1}_sample_list.txt", "r") as fp:
                self.assertEqual(2, len(fp.readlines()))

        exp_remap = [
            "#!/bin/bash\n",
            f"#SBATCH -J re_{job_id}\n",
            "#SBATCH -p qiita\n",
            "#SBATCH -N 1\n",
            "#SBATCH -n 12\n",
            "#SBATCH --time 7200\n",
            "#SBATCH --mem 50G\n",
            f"#SBATCH -o {out_dir}/remap/logs/%x-%A_%a.out\n",
            f"#SBATCH -e {out_dir}/remap/logs/%x-%A_%a.err\n",
            "#SBATCH --array 1-4%8\n",
            "\n",
            "source ~/.bashrc\n",
            "set -e\n",
            f"{CONDA_ENVIRONMENT}\n",
            "\n",
            "step=${SLURM_ARRAY_TASK_ID}\n",
            f"input=$(head -n $step {out_dir}/sample_list.txt | tail -n 1)\n",
            "sample_name=`echo $input | awk '{print $1}'`\n",
            "filename=`echo $input | awk '{print $2}'`\n",
            "aid=`echo $input | awk '{print $3}'`\n",
            "fn=`basename ${filename}`\n",
            "\n",
            f"out_folder={out_dir}/remap\n",
            "sn_folder=${out_folder}/bioms/${aid}.${sample_name}\n",
            "export TMPDIR=${out_folder}/tmp\n",
            "mkdir -p ${TMPDIR} ${sn_folder}\n",
            "\n",
            "txt=${sn_folder}/${aid}.${sample_name}.txt\n",
            "tsv=${txt/.txt/.tsv}\n",
            f"coverm genome -d {out_dir}/merge/dereplicated -x fna --mapper minimap2-hifi -m covered_bases trimmed_mean mean length count \\\n",
            "    --trim-min 0.05 --trim-max 0.95 --min-read-percent-identity 0.995 \\\n",
            "    --min-read-aligned-percent 0.90 --min-covered-fraction 0 -t 12 \\\n",
            "    --single ${filename} --output-file ${txt}\n",
            "\n",
            "awk 'BEGIN {FS=OFS=\"\\t\"}; {print $1,$NF}' ${txt} | \\\n",
            "    sed  's/dereplicated_gtdbtk\\///' | \\\n",
            "    sed  's/ Read Count//' | sed \"s/${fn}/${aid}.${sample_name}/\" > ${tsv}\n",
            "\n",
            "cov_tsv=${txt/.txt/.cov}\n",
            'echo -e "mag_id\\tcovered\\tlength\\tpercent_covered" > ${cov_tsv}\n',
            "tail -n +2 ${txt} | awk -F'\\t' 'BEGIN {FS=OFS=\"\\t\"}; { if ($5 > 0) print $1,$2,$5,$2/$5 }' >> ${cov_tsv}\n",
            "\n",
            "counts=`tail -n +2 ${tsv} | awk '{sum += $NF} END {print sum}'`\n",
            'if [[ "$counts" != "0" ]]; then\n',
            "    biom convert -i ${tsv} -o ${sn_folder}/counts.biom --to-hdf5\n",
            "fi\n",
            "\n",
            "touch ${out_folder}/completed_${SLURM_ARRAY_TASK_ID}.log",
        ]
        self.assertEqual(obs_remap, exp_remap)

        exp_finish = [
            "#!/bin/bash\n",
            f"#SBATCH -J f_{job_id}\n",
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
            f"python -c \"from qp_pacbio.util import client_connect; qclient = client_connect('{url}'); qclient.update_job_step('{job_id}', 'Merging BIOMs')\"\n",
            "\n",
            "lmap=${PWD}/length.map\n",
            "ls ${PWD}/*/bioms/*/*.cov > ${PWD}/cov_files.txt\n",
            "seqkit stats merge/dereplicated/*.fna | tail -n +2 | awk 'BEGIN{OFS=\"\\t\"} {print $1, $5}' | sed 's/,//g' | sed 's#merge/dereplicated/##g' | sed 's/.fna//g' > ${lmap}\n",
            "micov consolidate --paths ${PWD}/cov_files.txt --lengths ${lmap} --output ${PWD}/coverages.tgz\n",
            "\n",
            f"biom_merge_pacbio --base {out_dir}/remap --merge-type counts\n",
            "\n",
            f"finish_qp_pacbio {url} {job_id} {out_dir}",
        ]
        self.assertEqual(obs_finish, exp_finish)


if __name__ == "__main__":
    main()
