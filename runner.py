#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Copyright (c) 2025--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from qp_pacbio.qp_pacbio import generate_sample_list, pacbio_generate_templates
from qp_pacbio.util import client_connect

from os import makedirs
from os.path import join
import click
from subprocess import run

# Using woltka because we don't have a qp_pacbio command
plugin_details = {
    "name": "qp-woltka",
    "version": "2025.11",
    "description": "Woltka",
}

QIITA_URL = "https://qiita.ucsd.edu"


@click.command()
@click.option("--artifact_id", help="artifact id", required=True)
@click.option("--out_dir", help="output directory", required=True)
@click.option("--job_id", help="job id", required=True)
def runner(artifact_id, out_dir, job_id):
    def sbatch(args):
        res = run(args, check=True, capture_output=True, text=True)
        return res.stdout.strip()

    qclient = client_connect(QIITA_URL)
    njobs = generate_sample_list(qclient, artifact_id, out_dir)
    result_fp = join(out_dir, "results")
    makedirs(result_fp, exist_ok=True)
    pacbio_generate_templates(out_dir, job_id, njobs, QIITA_URL)
    print(qclient.artifact_and_preparation_files(artifact_id))

    jid1 = sbatch(["sbatch", "--parsable", f"{out_dir}/step-1/step-1.slurm"])
    jid2 = sbatch(
        [
            "sbatch",
            "--parsable",
            "--dependency",
            f"aftercorr:{jid1}",
            f"{out_dir}/step-2/step-2.slurm",
        ]
    )
    jid3 = sbatch(
        [
            "sbatch",
            "--parsable",
            "--dependency",
            f"aftercorr:{jid2}",
            f"{out_dir}/step-3/step-3.slurm",
        ]
    )
    jid4 = sbatch(
        [
            "sbatch",
            "--parsable",
            "--dependency",
            f"aftercorr:{jid3}",
            f"{out_dir}/step-4/step-4.slurm",
        ]
    )
    jid5 = sbatch(
        [
            "sbatch",
            "--parsable",
            "--dependency",
            f"aftercorr:{jid4}",
            f"{out_dir}/step-5/step-5.slurm",
        ]
    )
    jid6 = sbatch(
        [
            "sbatch",
            "--parsable",
            "--dependency",
            f"aftercorr:{jid5}",
            f"{out_dir}/step-6/step-6.slurm",
        ]
    )
    jid7 = sbatch(
        [
            "sbatch",
            "--parsable",
            "--dependency",
            f"aftercorr:{jid6}",
            f"{out_dir}/step-7/step-7.slurm",
        ]
    )

    submitted = [jid1, jid2, jid3, jid4, jid5, jid6, jid7]
    print("Submitted jobs: " + ", ".join(submitted))


if __name__ == "__main__":
    runner()
