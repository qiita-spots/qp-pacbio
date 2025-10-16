#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Copyright (c) 2025--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
import click
from qp_pacbio import plugin
from qp_pacbio.qp_pacbio import generate_minimap2_processing
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
    if 'register' in job_id:
        plugin(url, job_id, output_dir)
    else:
        qclient = client_connect(url)
        job_info = qclient.get_job_info(job_id)
        parameters = job_info['parameters']

        main_fp, merge_fp = generate_minimap2_processing(
                qclient, job_id, output_dir, parameters)

        # Submitting jobs and returning id
        main_job = run(['sbatch', main_fp], stdout=PIPE)
        main_job_id = main_job.stdout.decode('utf8').split()[-1]
        merge_job = run(['sbatch', '-d', f'afterok:{main_job_id}',
                         merge_fp], stdout=PIPE)
        merge_job_id = merge_job.stdout.decode('utf8').split()[-1]
        print(f'{main_job_id}, {merge_job_id}')

        qclient.update_job_step(
            job_id, "Step 2 of 4: Aligning sequences")
