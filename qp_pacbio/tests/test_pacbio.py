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

from qiita_client import ArtifactInfo
from qiita_client.testing import PluginTestCase

from qp_pacbio import plugin
from qp_pacbio.qp_pacbio import generate_minimap2_processing, pacbio_processing

# Keep these in sync with your generator defaults
CONDA_ENV = "qp_pacbio_2025.9"
STEP1_NPROCS = 16
STEP1_WALL = 1000
STEP1_MEM_GB = 300  # generator uses 300G for step-1
NODE_COUNT = 1
PARTITION = "qiita"

# Exact expected Step-1 script matching your template.
# Escape ${...} -> ${{...}} and awk braces -> '{{print $1}}' etc.
STEP_1_EXP = (
    "#!/bin/bash\n"
    "#SBATCH -J s1-{job_id}\n"
    f"#SBATCH -p {PARTITION}\n"
    f"#SBATCH -N {NODE_COUNT}\n"
    f"#SBATCH -n {STEP1_NPROCS}\n"
    f"#SBATCH --time {STEP1_WALL}\n"
    f"#SBATCH --mem {STEP1_MEM_GB}G\n"
    "#SBATCH -o {out_dir}/step-1/logs/%x-%A_%a.out\n"
    "#SBATCH -e {out_dir}/step-1/logs/%x-%A_%a.out\n"
    "#SBATCH --array 1-{njobs}%16\n"
    "\n"
    "source ~/.bashrc\n"
    "set -e\n"
    f"conda activate {CONDA_ENV}\n"
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
    f"hifiasm_meta -t {STEP1_NPROCS} -o "
    "{out_dir}/step-1/${{sample_name}} ${{filename}}"
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
        params = {"artifact": self._insert_data()}
        job_id = "my-job-id"
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # this should fail cause we don't have valid data
        result = pacbio_processing(self.qclient, job_id, params, out_dir)
        success, ainfo, msg = result

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

        self.assertTrue(success)
        exp = [
            ArtifactInfo(
                "output",
                "job-output-folder",
                [(f"{out_dir}/results/", "directory")],
            )
        ]
        self.assertCountEqual(ainfo, exp)


class PacWoltkaProfilingTests(PacBioTests):
    def test_pacbio_processing(self):
        params = {"artifact": int(self._insert_data())}
        job_id = "my-job-id"
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # this should fail cause we don't have valid data
        main_fp, merge_fp = generate_minimap2_processing(
            self.qclient, job_id, out_dir, params
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
            "#SBATCH --time 1000\n",
            "#SBATCH --mem 120G\n",
            f"#SBATCH -o {out_dir}/minimap2/logs/%x-%A_%a.out\n",
            f"#SBATCH -e {out_dir}/minimap2/logs/%x-%A_%a.out\n",
            "#SBATCH --array 1-2%16\n",
            "\n",
            "source ~/.bashrc\n",
            "set -e\n",
            "conda activate qp_pacbio_2025.9\n",
            f"mkdir -p {out_dir}/alignments\n",
            f"cd {out_dir}/\n",
            "db=/ddn_scratch/qiita_t/working_dir/tmp/db/WoLr2.mmi\n",
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
            "   samtools sort -@ 16 - | \\\n",
            '   awk \'BEGIN { FS=OFS="\\t" } /^@/ { print; next } '
            '{ $10="*"; $11="*" } 1\' | \\\n',
            f"   xz -1 -T1 > {out_dir}/alignments/${{sample_name}}.sam.xz",
        ]
        self.assertEqual(obs_main, exp_main)

        db_path = "/projects/wol/qiyun/wol2/databases/minimap2"
        exp_merge = [
            "#!/bin/bash\n",
            "#SBATCH -J me_my-job-id\n",
            "#SBATCH -p qiita\n",
            "#SBATCH -N 1\n",
            "#SBATCH -n 16\n",
            "#SBATCH --time 1000\n",
            "#SBATCH --mem 120G\n",
            f"#SBATCH -o {out_dir}/merge/logs/%x-%A_%a.out\n",
            f"#SBATCH -e {out_dir}/merge/logs/%x-%A_%a.out\n",
            "\n",
            "source ~/.bashrc\n",
            "set -e\n",
            "conda activate qp_pacbio_2025.9\n",
            f"cd {out_dir}/\n",
            f"tax={db_path}/WoLr2.tax\n",
            f"coords={db_path}/WoLr2.coords\n",
            f"len_map={db_path}/WoLr2/length.map\n",
            f"functional_dir={db_path}/WoLr2/\n",
            "\n",
            f"mkdir -p {out_dir}/coverages/\n",
            "\n",
            "for f in `ls alignments/*.sam.xz`; do\n",
            "    sn=`basename ${f/.sam.xz/}`;\n",
            f"    of={out_dir}/bioms/${{sn}};\n",
            "    mkdir -p ${of};\n",
            '    echo "woltka classify -i ${f} -o ${of}/none.biom '
            "--no-demux --lineage ${tax} "
            f'--rank none --outcov {out_dir}/coverages/";\n',
            '    echo "woltka classify -i ${f} -o ${of}/per-gene.biom '
            '--no-demux -c ${coords}";\n',
            "done | parallel -j 1\n",
            "wait\n",
            "\n",
            "for f in `ls bioms/*/per-gene.biom`; do\n",
            "    dn=`dirname ${f}`;\n",
            "    sn=`basename ${sn}`;\n",
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
            "done | parallel -j 1\n",
            "wait\n",
            "\n",
            "# MISSING:\n",
            "# merge bioms!\n",
            "\n",
            (
                f'find {out_dir}/coverages/ -iname "*.cov" '
                + f"> {out_dir}/cov_files.txt\n"
            ),
            f"micov consolidate --paths {out_dir}/cov_files.txt "
            f"--lengths ${{len_map}} --output {out_dir}/coverages.tgz\n",
            "\n",
            "cd alignments\n",
            "tar -cvf ../alignments.tar *.sam.xz",
        ]
        self.assertEqual(obs_merge, exp_merge)


if __name__ == "__main__":
    main()
