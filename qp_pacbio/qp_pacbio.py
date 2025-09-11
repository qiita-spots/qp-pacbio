# -----------------------------------------------------------------------------
# Copyright (c) 2025--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from qiita_client import ArtifactInfo
from jinja2 import Environment
from jinja2 import BaseLoader, TemplateNotFound
from os.path import basename, join, exists, getmtime
from os import makedirs
import pathlib


# taken from qp-woltka
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


# taken from https://jinja.palletsprojects.com/en/3.0.x/api/#jinja2.BaseLoader
class KISSLoader(BaseLoader):
    def __init__(self, path):
        # pin the path for loader to the location sequence_processing_pipeline
        # (the location of this file), along w/the relative path to the
        # templates directory.
        self.path = join(pathlib.Path(__file__).parent.resolve(), path)

    def get_source(self, environment, template):
        path = join(self.path, template)
        if not exists(path):
            raise TemplateNotFound(template)
        mtime = getmtime(path)
        with open(path) as f:
            source = f.read()
        return source, path, lambda: mtime == getmtime(path)


def generate_templates(out_dir, job_id, njobs):
    """Generates the templates

    Parameters
    ----------
    out_dir : str
        The path to the job's output directory
    job_id : str
        The job id
    njobs : int
        The number of total jobs
    """

    jinja_env = Environment(loader=KISSLoader('../data/templates'))
    template = jinja_env.get_template("1.hifiasm-meta_new.sbatch")
    cdir = f'{out_dir}/step-1'
    makedirs(cdir)
    with open(f'{cdir}/step-1.slurm', mode="w", encoding="utf-8") as f:
        f.write(template.render(
            conda_environment='qp_pacbio_2025.9',
            output=f'{out_dir}',
            job_name=f's1-{job_id}',
            node_count=1,
            nprocs=16,
            wall_time_limit=1000,
            mem_in_gb=16,
            array_params=f'1:{njobs}%16'
        ))


def generate_sample_list(qclient, artifact_id, out_dir):
    """Generates the sample_list.txt file

    Parameters
    ----------
    qclient : tgp.qiita_client.QiitaClient
        The Qiita server client
    artifact_id : int
        The artifact id
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    int
        The number of jobs/lines without header in sample_list.txt
    """
    files, prep = qclient.artifact_and_preparation_files(artifact_id)

    lookup = prep.set_index('run_prefix')['sample_name'].to_dict()
    lines = []
    for _, (fwd, _) in files.items():
        fwd = fwd['filepath']
        sname = search_by_filename(basename(fwd), lookup)
        lines.append(f'{sname}\t{fwd}')

    with open(f'{out_dir}/sample_list.txt', 'w') as f:
        f.write('\n'.join(lines))

    return len(lines)


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

    njobs = generate_sample_list(qclient, artifact_id, out_dir)

    qclient.update_job_step(
            job_id, "Step 2 of 3: Creating submission templates")

    njobs = generate_sample_list(qclient, artifact_id, out_dir)
    generate_templates(out_dir, job_id, njobs)

    # TODO 2. submit jobs

    qclient.update_job_step(
            job_id, "Step 3 of 3: Running commands")

    # TODO 3. wait for all jobs to finish

    qclient.update_job_step(
            job_id, "Commands finished")

    paths = [(f'{results_fp}/', 'directory')]
    return (True, [ArtifactInfo('output', 'job-output-folder', paths)],
            f'{job_id} completed.')
