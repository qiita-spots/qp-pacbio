#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Copyright (c) 2025--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
import click
from os import makedirs
from os.path import join
import biom as _biom
from glob import glob
import h5py


from qp_pacbio import plugin
from qp_pacbio.qp_pacbio import (
    generate_minimap2_processing,
    generate_sample_list,
    pacbio_generate_templates,
    PACBIO_PROCESSING_STEPS,
)
from qp_pacbio.util import client_connect
from subprocess import run, PIPE


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

        if command == "Woltka v0.1.7, minimap2":
            main_fp, merge_fp = generate_minimap2_processing(
                qclient, job_id, output_dir, parameters
            )

            # Submitting jobs and returning id
            main_job = run(["sbatch", main_fp], stdout=PIPE)
            main_job_id = main_job.stdout.decode("utf8").split()[-1]
            cmd = ["sbatch", "-d", f"afterok:{main_job_id}", merge_fp]
            merge_job = run(cmd, stdout=PIPE)
            merge_job_id = merge_job.stdout.decode("utf8").split()[-1]
            print(f"{main_job_id}, {merge_job_id}")
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


@click.command()
@click.option("--base", type=click.Path(exists=True), required=True)
def biom_merge(base):
    """Merges all PacBio biom tables"""
    chunk_size = 30
    for rank in ("none", "per-gene", "ko", "ec", "pathway"):
        rank = rank + ".biom"
        tables = glob(f"{base}/bioms/*/{rank}")

        if not tables:
            continue

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
                        for i, l in enumerate(fh):
                            if i >= 1 and l:
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

            with h5py.File(f"{base}/{rank}", "w") as out:
                full.to_hdf5(out, "fast-merge")
