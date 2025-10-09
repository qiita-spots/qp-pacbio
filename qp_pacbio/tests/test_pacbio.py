# -----------------------------------------------------------------------------
# Copyright (c) 2020--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

import difflib
import os
import re
from json import dumps
from os import remove
from os.path import exists, isdir, join
from shutil import copyfile, rmtree
from tempfile import mkdtemp
from unittest import main

from qiita_client import ArtifactInfo
from qiita_client.testing import PluginTestCase

from qp_pacbio import plugin
from qp_pacbio.qp_pacbio import pacbio_processing

# Keep these in sync with your generator defaults
CONDA_ENV = "qp_pacbio_2025.9"
STEP1_NPROCS = 16
STEP1_WALL = 1000
STEP1_MEM_GB = 300
NODE_COUNT = 1
PARTITION = "qiita"

# Order-sensitive expected template for step-1, matching your current Jinja
# (see the step-1 slurm you posted). IMPORTANT: we must escape any literal
# ${...} shell variable by doubling braces so str.format leaves them intact.
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
    "source ~/.bashrc\n"
    f"conda activate {CONDA_ENV}\n"
    "\n"
    "cd {out_dir}/step-1\n"
    "step=${{SLURM_ARRAY_TASK_ID}}\n"
    "input=$(head -n $step {out_dir}/sample_list.txt | tail -n 1)\n"
    "\n"
    "sample_name=`echo $input | awk '{print $1}'`\n"
    "filename=`echo $input | awk '{print $2}'`\n"
    "\n"
    "fn=`basename ${{filename}}`\n"
    f"hifiasm_meta -t {STEP1_NPROCS} -o "
    "{{out_dir}}/step-1/${{sample_name}} ${{filename}}"
)


class PacBioTests(PluginTestCase):
    # ---------- Test Harness Utilities ----------
    def setUp(self):
        # register plugin test server & keep full diffs visible
        plugin("https://localhost:21174", "register", "ignored")
        self._clean_up_files = []
        self.maxDiff = None

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    try:
                        rmtree(fp)
                    except Exception:
                        # best-effort cleanup
                        for root, dirs, files in os.walk(fp, topdown=False):
                            for name in files:
                                try:
                                    os.remove(os.path.join(root, name))
                                except Exception:
                                    pass
                            for name in dirs:
                                try:
                                    os.rmdir(os.path.join(root, name))
                                except Exception:
                                    pass
                        try:
                            os.rmdir(fp)
                        except Exception:
                            pass
                else:
                    try:
                        remove(fp)
                    except Exception:
                        pass

    # ---------- Fixtures ----------
    def _insert_data(self):
        prep_info_dict = {
            "SKB8.640193": {"run_prefix": "S22205_S104"},
            "SKD8.640184": {"run_prefix": "S22282_S102"},
        }
        data = {
            "prep_info": dumps(prep_info_dict),
            "study": 1,  # magic #1 = testing study
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

    # ---------- Helpers for robust comparisons ----------
    @staticmethod
    def _normalize_lines(lines, out_dir):
        """Normalize lines to remove run-to-run variance.

        - Strip trailing newlines
        - Collapse runs of spaces/tabs
        - Replace absolute out_dir with <OUT_DIR>
        - Normalize datetimes and UUIDs if present
        """
        norm = []
        for ln in lines:
            ln = ln.rstrip("\n")
            # collapse repeated whitespace (keep a single space)
            ln = re.sub(r"[ \t]+", " ", ln).rstrip()
            if out_dir:
                ln = ln.replace(out_dir, "<OUT_DIR>")
            # normalize common datetime patterns
            ln = re.sub(
                r"\b\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2}\b",
                "<DATETIME>",
                ln,
            )
            # normalize UUIDs
            ln = re.sub(
                r"\b[a-f0-9]{8}(-[a-f0-9]{4}){3}-[a-f0-9]{12}\b",
                "<UUID>",
                ln,
                flags=re.IGNORECASE,
            )
            norm.append(ln)
        return norm

    def _assert_equal_with_diff(self, expected_lines, observed_lines, out_dir):
        """Assert equality and, on failure, print a unified diff (normalized)."""
        expN = self._normalize_lines(expected_lines, out_dir)
        obsN = self._normalize_lines(observed_lines, out_dir)
        try:
            self.assertEqual(expN, obsN)
        except AssertionError:
            diff = "\n".join(
                difflib.unified_diff(
                    expN,
                    obsN,
                    fromfile="expected",
                    tofile="observed",
                    lineterm="",
                )
            )
            print("\n==== Unified diff (normalized) ====\n" + diff)
            raise

    # ---------- The Actual Test ----------
    def test_pacbio_processing(self):
        params = {"artifact_id": self._insert_data()}
        job_id = "my-job-id"
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # This renders the SLURM scripts; it should not actually run jobs.
        success, ainfo, msg = pacbio_processing(
            self.qclient, job_id, params, out_dir
        )
        self.assertTrue(success, msg=f"pacbio_processing failed: {msg}")

        # sample list created
        with open(f"{out_dir}/sample_list.txt", "r") as f:
            sample_lines = [ln.rstrip("\n") for ln in f.readlines()]
        self.assertEqual(2, len(sample_lines))
        njobs = len(sample_lines)
        self.assertGreater(njobs, 0, "No jobs derived from sample list.")

        # ----- Step 1: exact, order-sensitive comparison -----
        step1_path = f"{out_dir}/step-1/step-1.slurm"
        with open(step1_path, "r") as f:
            obs_lines = [ln.rstrip("\n") for ln in f.readlines()]

        exp_text = STEP_1_EXP.format(
            job_id=job_id,
            out_dir=out_dir,
            njobs=njobs,
        )
        exp_lines = exp_text.split("\n")

        self._assert_equal_with_diff(exp_lines, obs_lines, out_dir)

        # ----- Other steps: header sanity checks only -----
        def assert_header(step, nprocs, wall, mem_gb):
            slurm_fp = f"{out_dir}/step-{step}/step-{step}.slurm"
            with open(slurm_fp, "r") as fh:
                got = [ln.rstrip("\n") for ln in fh.readlines()]

            expected_subset = [
                "#!/bin/bash",
                f"#SBATCH -J s{step}-{job_id}",
                f"#SBATCH -p {PARTITION}",
                f"#SBATCH -N {NODE_COUNT}",
                f"#SBATCH -n {nprocs}",
                f"#SBATCH --time {wall}",
                f"#SBATCH --mem {mem_gb}G",
                (
                    "#SBATCH -o "
                    f"{out_dir}/step-{step}/logs/%x-%A_%a.out"
                ),
                (
                    "#SBATCH -e "
                    f"{out_dir}/step-{step}/logs/%x-%A_%a.out"
                ),
                f"#SBATCH --array 1-{njobs}%16",
                "source ~/.bashrc",
                f"conda activate {CONDA_ENV}",
                f"cd {out_dir}/step-{step}",
            ]
            for line in expected_subset:
                with self.subTest(step=step, line=line):
                    self.assertIn(
                        line,
                        got,
                        msg=(f"Missing line in step-{step} header: {line}"),
                    )

        # step-0
        assert_header(step=0, nprocs=16, wall=1000, mem_gb=300)
        # step-2
        assert_header(step=2, nprocs=1, wall=500, mem_gb=16)
        # step-3
        assert_header(step=3, nprocs=8, wall=500, mem_gb=50)
        # step-4
        assert_header(step=4, nprocs=8, wall=500, mem_gb=50)
        # step-5
        assert_header(step=5, nprocs=8, wall=500, mem_gb=50)
        # step-6
        assert_header(step=6, nprocs=8, wall=500, mem_gb=50)
        # step-7
        assert_header(step=7, nprocs=8, wall=500, mem_gb=50)

        # ainfo / results folder
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
