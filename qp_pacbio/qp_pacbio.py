# -----------------------------------------------------------------------------
# Copyright (c) 2025--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from qiita_client import ArtifactInfo
from os import makedirs


def pacbio_processing(qclient, job_id, parameters, out_dir):
    """Sequence Processing Pipeline command

    Parameters
    ----------
    qclient : tgp.qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values for this job
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    bool, list, str
        The results of the job
    """
    results_fp = f'{out_dir}/results'
    makedirs(results_fp)

    paths = [(f'{results_fp}/', 'directory')]
    return (True, [ArtifactInfo('output', 'job-output-folder', paths)],
            f'{job_id} completed.')
