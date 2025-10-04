#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Copyright (c) 2025--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from qp_pacbio.qp_pacbio import generate_sample_list, generate_templates
from qiita_client import QiitaClient
from configparser import ConfigParser
from os import environ
from os.path import join, expanduser
import click
from subprocess import run, PIPE

# Using woltka because we don't have a qp_pacbio command
plugin_details = {"name": "qp-woltka", "version": "2024.09", "description": "Woltka"}


def client_connect(url):
    name = plugin_details["name"]
    version = plugin_details["version"]

    config = ConfigParser()
    conf_dir = environ.get("QIITA_PLUGINS_DIR", join(expanduser("~"), ".qiita_plugins"))
    conf_fp = join(conf_dir, f"{name}_{version}.conf")

    with open(conf_fp, "r") as conf_file:
        config.readfp(conf_file)
    qclient = QiitaClient(
        url,
        config.get("oauth2", "CLIENT_ID"),
        config.get("oauth2", "CLIENT_SECRET"),
        config.get("oauth2", "SERVER_CERT"),
    )

    return qclient


@click.command()
@click.option("--artifact_id", help="artifact id", required=True)
@click.option("--out_dir", help="output directory", required=True)
@click.option("--job_id", help="job id", required=True)
def runner(artifact_id, out_dir, job_id):
    def sbatch(args):
        r = run(args, check=True, capture_output=True, text=True)
        return r.stdout.strip()

    qclient = client_connect("https://qiita.ucsd.edu")
    njobs = generate_sample_list(qclient, artifact_id, out_dir)
    generate_templates(out_dir, job_id, njobs)
    print(qclient.artifact_and_preparation_files(artifact_id))

    jid0 = sbatch(["sbatch", "--parsable", f"{out_dir}/step-0/step-0.slurm"])
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

    print(
        f"Submitted jobs: {jid0}, {jid1}, {jid2}, {jid3}, {jid4}, {jid5}, {jid6}, {jid7}"
    )


if __name__ == "__main__":
    runner()
