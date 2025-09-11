# -----------------------------------------------------------------------------
# Copyright (c) 2025--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from qiita_client import ArtifactInfo
from os import makedirs
from os.path import basename


def search_by_filename(fname, lookup):
    if fname in lookup:
        return lookup[fname]

    original = fname
    while '_' in fname:
        fname = fname.rsplit('_', 1)[0]
        if fname in lookup:
            return lookup[fname]

    fname = original
    while '.' in fname:
        fname = fname.rsplit('.', 1)[0]
        if fname in lookup:
            return lookup[fname]

    for rp in lookup:
        if original.startswith(rp):
            return lookup[rp]

    raise KeyError("Cannot determine run_prefix for %s" % original)


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

    qclient.update_job_step(
            job_id, "Step 1 of 3: Collecting info and generating submission")
    artifact_id = parameters['artifact_id']
    files, prep = qclient.artifact_and_preparation_files(artifact_id)
    lookup = prep.set_index('run_prefix')['sample_name'].to_dict()
    lines = ['sample_name\tfilename']
    for _, (fwd, _) in files.items():
        sname = search_by_filename(basename(fwd['filepath']), lookup)
        lines.append(f'{sname}\t{fwd}')

    with open(f'{out_dir}/sample_list.txt') as f:
        f.write('\n'.append(lines))

    qclient.update_job_step(
            job_id, "Step 2 of 3: Creating submission templates")
    # TODO 1. add code to create templates
    # TODO 2. submit templates

    qclient.update_job_step(
            job_id, "Step 3 of 3: Running commands")
    # TODO 3. wait for all jobs to finish

    qclient.update_job_step(
            job_id, "Commands finished")

    paths = [(f'{results_fp}/', 'directory')]
    return (True, [ArtifactInfo('output', 'job-output-folder', paths)],
            f'{job_id} completed.')
