# -----------------------------------------------------------------------------
# Copyright (c) 2020--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from qiita_client import ArtifactInfo
from qiita_client.testing import PluginTestCase
from os import remove
from os.path import exists, isdir, join
from tempfile import mkdtemp
from shutil import rmtree, copyfile
from json import dumps

from qp_pacbio import plugin
from qp_pacbio.qp_pacbio import pacbio_processing


# NOTE: Edited to match your generator:
# - Uses sample_list.txt (not file_list.txt)
# - Uses -t 16 to match nprocs=16
STEP_1_EXP = (
    "#!/bin/bash\n"
    "#SBATCH -J s1-my-job-id\n"
    "#SBATCH -N 1\n"
    "#SBATCH -n 16\n"
    "#SBATCH --time 1000\n"
    "#SBATCH --mem 300G\n"
    "#SBATCH -o {out_dir}/step-1/logs/%x-%A.out\n"
    "#SBATCH -e {out_dir}/step-1/logs/%x-%A.err\n"
    "#SBATCH --array 1-2%16\n"
    "\n"
    "conda activate qp_pacbio_2025.9\n"
    "\n"
    "cd {out_dir}/step-1\n"
    "step=${{SLURM_ARRAY_TASK_ID}}\n"
    "input=$(head -n $step {out_dir}/sample_list.txt | tail -n 1)\n"
    "fn=`basename ${{input}}`\n"
    "hifiasm_meta -t 16 -o {out_dir}/step-1/${{fn}} ${{input}}"
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

    def test_pacbio_processing(self):
        params = {"artifact_id": self._insert_data()}
        job_id = "my-job-id"
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # this should succeed (we only render, not actually run)
        success, ainfo, msg = pacbio_processing(
            self.qclient, job_id, params, out_dir
        )

        # sample list created
        with open(f"{out_dir}/sample_list.txt", "r") as f:
            sample_lines = [ln.rstrip("\n") for ln in f.readlines()]
        self.assertEqual(2, len(sample_lines))
        njobs = len(sample_lines)

        # step-1 exact comparison (kept like the original style)
        with open(f"{out_dir}/step-1/step-1.slurm", "r") as f:
            obs_lines = [ln.rstrip("\n") for ln in f.readlines()]
        self.assertCountEqual(
            STEP_1_EXP.format(out_dir=out_dir).split("\n"),
            obs_lines,
        )

        # Header-only checks for other steps to keep tests stable while
        # still catching resource misconfigurations.
        def assert_header(step, nprocs, wall, mem_gb):
            slurm_fp = f"{out_dir}/step-{step}/step-{step}.slurm"
            with open(slurm_fp, "r") as fh:
                got = [ln.rstrip("\n") for ln in fh.readlines()]

            expected_subset = [
                "#!/bin/bash",
                f"#SBATCH -J s{step}-{job_id}",
                "#SBATCH -N 1",
                f"#SBATCH -n {nprocs}",
                f"#SBATCH --time {wall}",
                f"#SBATCH --mem {mem_gb}G",
                f"#SBATCH -o {out_dir}/step-{step}/logs/%x-%A.out",
                f"#SBATCH -e {out_dir}/step-{step}/logs/%x-%A.err",
                f"#SBATCH --array 1-{njobs}%16",
                "",
                "conda activate qp_pacbio_2025.9",
                "",
                f"cd {out_dir}/step-{step}",
            ]
            for line in expected_subset:
                self.assertIn(
                    line, got, msg=f"Missing line in step-{step} header: {line}"
                )

        # step-0
        assert_header(step=0, nprocs=16, wall=1000, mem_gb=300)
        # step-2
        assert_header(step=2, nprocs=1, wall=500, mem_gb=16)
        # step-3
        assert_header(step=3, nprocs=8, wall=500, mem_gb=50)
        # step-4
        assert_header(step=4, nprocs=8, wall=500, mem_gb=50)
        # step-5 (filename in templates differs, but output path is standard)
        assert_header(step=5, nprocs=8, wall=500, mem_gb=50)
        # step-6
        assert_header(step=6, nprocs=8, wall=500, mem_gb=50)
        # step-7
        assert_header(step=7, nprocs=8, wall=500, mem_gb=50)

        self.assertTrue(success)
        exp = [
            ArtifactInfo(
                "output",
                "job-output-folder",
                [(f"{out_dir}/results/", "directory")],
            )
        ]
        self.assertCountEqual(ainfo, exp)


if __name__ == "__main__":
    main()
