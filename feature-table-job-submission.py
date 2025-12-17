#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Copyright (c) 2025--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
import click
from qiita_db.analysis import Analysis
from qiita_db.artifact import Artifact
from qiita_db.processing_job import ProcessingJob
from qiita_db.software import Command, Parameters
from qiita_db.user import User


@click.command()
@click.option("--aid", help="artifact id", required=True, multiple=True)
@click.option("--user-email", help="user email", required=True)
@click.option(
    "--command-id",
    help="the feature table generation command id",
    required=False,
    default=850,
)
def main(aid, user_email, command_id):
    user = User(user_email)
    artifacts = [Artifact(a) for a in aid]

    command = Command(command_id)

    analysis = Analysis.create(user, command.name, "")
    samples = {a.id: list(a.prep_templates[0]) for a in artifacts}
    analysis.add_samples(samples)
    analysis._build_mapping_file(samples, rename_dup_samples=False)

    params = {
        "analysis": analysis.id,
        "artifacts": [a for a in aid],
    }
    job_params = Parameters.load(command, values_dict=params)
    job = ProcessingJob.create(user, job_params, True)
    job.submit()

    print(f"{job.id} created for Analysis {analysis.id}.")


if __name__ == "__main__":
    main()
