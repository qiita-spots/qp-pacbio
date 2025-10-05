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
    plugin(url, job_id, output_dir)
