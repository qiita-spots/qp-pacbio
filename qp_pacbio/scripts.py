#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Copyright (c) 2025--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from glob import glob
from os import makedirs
from os.path import join
from subprocess import PIPE, run

import biom as _biom
import click
import h5py

from qp_pacbio import plugin
from qp_pacbio.qp_pacbio import (
    PACBIO_PROCESSING_STEPS,
    generate_minimap2_processing,
    generate_pacbio_adapter_removal,
    generate_sample_list,
    generate_syndna_processing,
    pacbio_generate_templates,
)
from qp_pacbio.util import client_connect


@click.command()
@click.option("--env-script", prompt="Environment script", required=True)
@click.option(
    "--ca-cert",
    prompt="Server certificate",
    required=False,
    default="None",
    show_default=True,
)
def config(env_script, ca_cert):
    """Generates the Qiita configuration files"""
    if ca_cert == "None":
        ca_cert = None
    plugin.generate_config(env_script, "start_qp_pacbio", ca_cert)


@click.command()
@click.argument("url", required=True)
@click.argument("job_id", required=True)
@click.argument("output_dir", required=True)
def execute(url, job_id, output_dir):
    """Executes the task given by job_id and puts the output in output_dir"""
    if "register" in job_id:
        plugin(url, job_id, output_dir)
    else:
        qclient = client_connect(url)
        job_info = qclient.get_job_info(job_id)
        parameters = job_info["parameters"]
        command = job_info["command"]
        artifact_id = parameters["artifact"]

        regular_commands = {
            "Woltka v0.1.7, minimap2": generate_minimap2_processing,
            "Remove SynDNA plasmid, insert, & GCF_000184185 reads (minimap2)": generate_syndna_processing,
            "PacBio adapter removal via lima/pbmarkdup": generate_pacbio_adapter_removal,
        }

        if command in regular_commands.keys():
            first_fp, second_fp = regular_commands[command](
                qclient, job_id, output_dir, parameters, url
            )

            # Submitting jobs and returning id
            first_job = run(["sbatch", first_fp], stdout=PIPE)
            first_job_id = first_job.stdout.decode("utf8").split()[-1]
            cmd = ["sbatch", "-d", f"afterok:{first_job_id}", second_fp]
            second_job = run(cmd, stdout=PIPE)
            second_job_id = second_job.stdout.decode("utf8").split()[-1]
            print(f"{first_job_id}, {second_job_id}")
            qclient.update_job_step(job_id, "Step 2 of 4: Aligning sequences")
        elif command == "PacBio processing":
            frp = join(output_dir, "results")
            makedirs(frp, exist_ok=True)

            qclient.update_job_step(job_id, "Generating commands.")

            njobs = generate_sample_list(qclient, artifact_id, output_dir)
            pacbio_generate_templates(output_dir, job_id, njobs, frp, url)

            all_jids = []
            for step_id in range(1, PACBIO_PROCESSING_STEPS + 1, 1):
                cmd = [
                    "sbatch",
                    "--parsable",
                    f"{output_dir}/step-{step_id}/step-{step_id}.slurm",
                ]
                if step_id != 1:
                    cmd.insert(2, "--dependency")
                    cmd.insert(3, f"aftercorr:{all_jids[-1]}")
                task = run(cmd, stdout=PIPE)
                jid = task.stdout.decode("utf8").split()[-1]
                all_jids.append(jid)
            # adding finish command
            cmd = [
                "sbatch",
                "--parsable",
                "--dependency",
                f"afterany:{all_jids[-1]}",
                f"{output_dir}/finish/finish.slurm",
            ]
            task = run(cmd, stdout=PIPE)
            jid = task.stdout.decode("utf8").split()[-1]
            all_jids.append(jid)
            print(", ".join(all_jids))
        else:
            # this should never happen but rather have it
            raise ValueError(f"{command} not implemented!")


@click.command()
@click.argument("url", required=True)
@click.argument("job_id", required=True)
@click.argument("output_dir", required=True)
def finish_qp_pacbio(url, job_id, output_dir):
    """Executes the task given by job_id and puts the output in output_dir"""
    plugin(url, job_id, output_dir)


def _biom_merge(tables):
    chunk_size = 30
    full = None
    for block in range(0, len(tables), chunk_size):
        lblock = block + chunk_size
        chunk = tables[block:lblock]
        loaded = []
        for c in chunk:
            skip = True
            if _biom.util.is_hdf5_file(c):
                skip = False
            else:
                with open(c) as fh:
                    for i, _l in enumerate(fh):
                        if i >= 1 and _l:
                            skip = False
                            break
            if not skip:
                temp = _biom.load_table(c)
                if temp.shape != (0, 0):
                    loaded.append(temp)

        if full is None:
            if len(loaded) == 1:
                full = loaded[0]
            else:
                full = loaded[0].concat(loaded[1:])
        else:
            full = full.concat(loaded)

    return full


@click.command()
@click.option("--base", type=click.Path(exists=True), required=True)
@click.option(
    "--merge-type", type=click.Choice(["syndna", "woltka"], case_sensitive=False)
)
def biom_merge(base, merge_type):
    """Merges all PacBio biom tables"""
    merge_type = merge_type.lower()
    if merge_type == "syndna":
        ranks = ["syndna"]
    elif merge_type == "woltka":
        ranks = ["none", "per-gene", "ko", "ec", "pathway"]
    else:
        raise ValueError(f"Type '{merge_type}' not supported")

    for rank in ranks:
        rank = rank + ".biom"
        tables = glob(f"{base}/bioms/*/{rank}")

        if not tables:
            continue

        full = _biom_merge(tables)
        with h5py.File(f"{base}/{rank}", "w") as out:
            full.to_hdf5(out, "fast-merge")
