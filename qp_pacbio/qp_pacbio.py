# -----------------------------------------------------------------------------
# Copyright (c) 2025--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from os import makedirs, mkdir
from os.path import basename, join, exists, getmtime
import pathlib
from shutil import copy2


from jinja2 import Environment, BaseLoader, TemplateNotFound
from qiita_client import ArtifactInfo
from subprocess import run


CONDA_ENV = "qp_pacbio_2025.9"
MAX_WALL_1000 = 1000
MAX_WALL_500 = 500
T5_NAME = "5.DAS_Tools_prepare_batch3_test.sbatch"


def sbatch(args):
    res = run(args, check=True, capture_output=True, text=True)
    return res.stdout.strip()


# taken from qp-woltka
def search_by_filename(fname, lookup):
    if fname in lookup:
        return lookup[fname]

    original = fname
    while "_" in fname:
        fname = fname.rsplit("_", 1)[0]
        if fname in lookup:
            return lookup[fname]

    fname = original
    while "." in fname:
        fname = fname.rsplit(".", 1)[0]
        if fname in lookup:
            return lookup[fname]

    for rp in lookup:
        if original.startswith(rp):
            return lookup[rp]

    raise KeyError("Cannot determine run_prefix for %s" % original)


# taken from the Jinja docs (BaseLoader API):
# https://jinja.palletsprojects.com/en/3.0.x/api/
class KISSLoader(BaseLoader):
    def __init__(self, path):
        base = pathlib.Path(__file__).parent.resolve()
        self.path = join(base, path)

    def get_source(self, environment, template):
        path = join(self.path, template)
        if not exists(path):
            raise TemplateNotFound(template)
        mtime = getmtime(path)
        with open(path, encoding="utf-8") as f:
            source = f.read()
        return source, path, lambda: mtime == getmtime(path)


def _write_slurm(path, template, **ctx):
    makedirs(path, exist_ok=True)
    makedirs(join(path, "logs"), exist_ok=True)
    out_fp = join(path, f"{path.rsplit('/', 1)[-1]}.slurm")
    rendered = template.render(**ctx)
    with open(out_fp, mode="w", encoding="utf-8") as f:
        f.write(rendered)

    return out_fp


def pacbio_generate_templates(out_dir, job_id, njobs):
    """Generate Slurm submission templates for PacBio processing.

    Parameters
    ----------
    out_dir : str
        Path to the job's output directory.
    job_id : str
        Qiita job id.
    njobs : int
        Number of array tasks/jobs.
    """
    jinja_env = Environment(loader=KISSLoader("../data/templates"))

    # Step 1
    template1 = jinja_env.get_template("1.hifiasm-meta_new.sbatch")
    _write_slurm(
        join(out_dir, "step-1"),
        template1,
        conda_environment=CONDA_ENV,
        output=out_dir,
        job_name=f"s1-{job_id}",
        node_count=1,
        nprocs=16,
        wall_time_limit=MAX_WALL_1000,
        mem_in_gb=300,
        array_params=f"1-{njobs}%16",
    )

    # Step 2
    template2 = jinja_env.get_template("2.get-circular-genomes.sbatch")
    _write_slurm(
        join(out_dir, "step-2"),
        template2,
        conda_environment=CONDA_ENV,
        output=out_dir,
        job_name=f"s2-{job_id}",
        node_count=1,
        nprocs=1,
        wall_time_limit=MAX_WALL_500,
        mem_in_gb=16,
        array_params=f"1-{njobs}%16",
    )

    # Step 3
    template3 = jinja_env.get_template("3.minimap2_assembly.sbatch")
    _write_slurm(
        join(out_dir, "step-3"),
        template3,
        conda_environment=CONDA_ENV,
        output=out_dir,
        job_name=f"s3-{job_id}",
        node_count=1,
        nprocs=8,
        wall_time_limit=MAX_WALL_500,
        mem_in_gb=50,
        array_params=f"1-{njobs}%16",
    )

    # Step 4
    template4 = jinja_env.get_template("4.metawrap_binning_new.sbatch")
    _write_slurm(
        join(out_dir, "step-4"),
        template4,
        conda_environment=CONDA_ENV,
        output=out_dir,
        job_name=f"s4-{job_id}",
        node_count=1,
        nprocs=8,
        wall_time_limit=MAX_WALL_500,
        mem_in_gb=50,
        array_params=f"1-{njobs}%16",
    )

    # Step 5 (long filename isolated in T5_NAME)
    template5 = jinja_env.get_template(T5_NAME)
    _write_slurm(
        join(out_dir, "step-5"),
        template5,
        conda_environment=CONDA_ENV,
        output=out_dir,
        job_name=f"s5-{job_id}",
        node_count=1,
        nprocs=8,
        wall_time_limit=MAX_WALL_500,
        mem_in_gb=50,
        array_params=f"1-{njobs}%16",
    )

    # Step 6
    template6 = jinja_env.get_template("6.MAG_rename.sbatch")
    _write_slurm(
        join(out_dir, "step-6"),
        template6,
        conda_environment=CONDA_ENV,
        output=out_dir,
        job_name=f"s6-{job_id}",
        node_count=1,
        nprocs=8,
        wall_time_limit=MAX_WALL_500,
        mem_in_gb=50,
        array_params=f"1-{njobs}%16",
    )

    # Step 7
    template7 = jinja_env.get_template("7.checkm_batch.sbatch")
    _write_slurm(
        join(out_dir, "step-7"),
        template7,
        conda_environment=CONDA_ENV,
        output=out_dir,
        job_name=f"s7-{job_id}",
        node_count=1,
        nprocs=8,
        wall_time_limit=MAX_WALL_500,
        mem_in_gb=50,
        array_params=f"1-{njobs}%16",
    )


def generate_sample_list(qclient, artifact_id, out_dir):
    """Create sample_list.txt of sample names/files for slurm arrays.

    Parameters
    ----------
    qclient : tgp.qiita_client.QiitaClient
        Qiita server client.
    artifact_id : int
        Artifact id.
    out_dir : str
        Output directory.

    Returns
    -------
    int
        Number of non-header lines written.
    """
    files, prep = qclient.artifact_and_preparation_files(artifact_id)

    lookup = prep.set_index("run_prefix")["sample_name"].to_dict()
    lines = []
    for _, (fwd, _) in files.items():
        fwd_fp = fwd["filepath"]
        sname = search_by_filename(basename(fwd_fp), lookup)
        lines.append(f"{sname}\t{fwd_fp}")

    out_fp = join(out_dir, "sample_list.txt")
    with open(out_fp, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    return len(lines)


def pacbio_processing(qclient, job_id, parameters, out_dir):
    """Sequence Processing Pipeline command.

    Parameters
    ----------
    qclient : tgp.qiita_client.QiitaClient
        Qiita server client.
    job_id : str
        Job id.
    parameters : dict
        Parameters for this job.
    out_dir : str
        Output directory.

    Returns
    -------
    bool, list, str
        Results tuple for Qiita.
    """
    results_fp = join(out_dir, "results")
    makedirs(results_fp, exist_ok=True)

    qclient.update_job_step(
        job_id,
        "Step 1 of 3: Collecting info and generating submission",
    )
    artifact_id = parameters["artifact_id"]

    njobs = generate_sample_list(qclient, artifact_id, out_dir)

    qclient.update_job_step(
        job_id,
        "Step 2 of 3: Creating submission templates",
    )
    pacbio_generate_templates(out_dir, job_id, njobs)

    # If/when you enable submission, capture Slurm job IDs and thread them:
    # jid0 = qclient.submit_job(f"{out_dir}/step-0/step-0.slurm")
    # jid1 = qclient.submit_job(f"{out_dir}/step-1/step-1.slurm")
    # For now, keep None to avoid F821.
    jid1 = None

    # Re-render Step 2 with dependency on the matching task of Step 1.
    jinja_env = Environment(loader=KISSLoader("../data/templates"))
    template2 = jinja_env.get_template("2.get-circular-genomes.sbatch")

    cdir2 = join(out_dir, "step-2")
    makedirs(cdir2, exist_ok=True)

    dep = None
    if jid1:
        # Example: chain with Slurm's aftercorr
        dep = f"aftercorr:{jid1}"

    rendered = template2.render(
        conda_environment=CONDA_ENV,
        output=out_dir,
        job_name=f"s2-{job_id}",
        node_count=1,
        nprocs=1,
        wall_time_limit=MAX_WALL_500,
        mem_in_gb=16,
        # Use "1:{njobs}%16" if template expects a colon separator.
        array_params=f"1-{njobs}%16",
        dependency=dep,
    )
    with open(join(cdir2, "step-2.slurm"), "w", encoding="utf-8") as f:
        f.write(rendered)

    qclient.update_job_step(job_id, "Step 3 of 3: Running commands")

    # TODO: wait for jobs to finish (future)
    qclient.update_job_step(job_id, "Commands finished")

    paths = [(f"{results_fp}/", "directory")]
    return (
        True,
        [ArtifactInfo("output", "job-output-folder", paths)],
        f"{job_id} completed.",
    )


def minimap2_processing(qclient, job_id, parameters, out_dir):
    """generates output for minimap2/woltka processing.

    Parameters
    ----------
    qclient : tgp.qiita_client.QiitaClient
        Qiita server client.
    job_id : str
        Job id.
    parameters : dict
        Parameters for this job.
    out_dir : str
        Output directory.

    Returns
    -------
    bool, list, str
        Results tuple for Qiita.
    """
    qclient.update_job_step(job_id, "Commands finished")

    def _coverage_copy(dest):
        fp_coverages = join(out_dir, 'coverages.tgz')
        mkdir(dest)
        dest = join(dest, 'coverages.tgz')
        copy2(fp_coverages, dest)

        return dest

    errors = []
    ainfo = []

    fp_biom = f'{out_dir}/none.biom'
    fp_alng = f'{out_dir}/alignment.tar'
    if exists(fp_biom) and exists(fp_alng):
        ainfo.append(ArtifactInfo('Per genome Predictions', 'BIOM', [
            (fp_biom, 'biom'), (fp_alng, 'log'),
            (_coverage_copy(f'{out_dir}/none/'), 'plain_text')]))
    else:
        errors.append('Table none/per-genome was not created, please contact '
                      'qiita.help@gmail.com for more information')

    bioms = [
        (f'{out_dir}/per-gene.biom', 'per_gene' 'Per gene Predictions'),
        (f'{out_dir}/ko.biom', 'ko' 'KEGG Ontology (KO)'),
        (f'{out_dir}/ec.biom', 'ec' 'KEGG Enzyme (EC)'),
        (f'{out_dir}/pathway.biom', 'pathway' 'KEGG Pathway'),
    ]

    for fb, fn, bn in bioms:
        if exists(fb):
            ainfo.append(ArtifactInfo(bn, 'BIOM', [
                (fb, 'biom'),
                (_coverage_copy(f'{out_dir}/{fn}/'), 'plain_text')]))
        else:
            errors.append(f'Table "{bn}" was not created, please contact '
                          'qiita.help@gmail.com for more information')

    if errors:
        return False, ainfo, '\n'.join(errors)
    else:

        return True, ainfo, ""


def generate_minimap2_processing(qclient, job_id, out_dir, parameters):
    """generates slurm scripts for minimap2/woltka processing.

    Parameters
    ----------
    qclient : tgp.qiita_client.QiitaClient
        Qiita server client.
    job_id : str
        Job id.
    out_dir : str
        Output directory.
    parameters : dict
        Parameters for this job.

    Returns
    -------
    str, str
        Returns the two filepaths of the slurm scripts
    """
    qclient.update_job_step(
        job_id, "Step 1 of 4: Collecting info and generating submission")

    artifact_id = parameters["artifact_id"]

    njobs = generate_sample_list(qclient, artifact_id, out_dir)

    qclient.update_job_step(
        job_id,
        "Step 2 of 4: Creating submission templates",
    )

    jinja_env = Environment(loader=KISSLoader("../data/templates"))
    minimap2_template = jinja_env.get_template("woltka_minimap2.sbatch")
    minimap2_script = _write_slurm(
        f'{out_dir}/minimap2',
        minimap2_template,
        conda_environment=CONDA_ENV,
        output=out_dir,
        job_name=f"m2_{job_id}",
        node_count=1,
        nprocs=16,
        wall_time_limit=MAX_WALL_1000,
        mem_in_gb=120,
        array_params=f"1-{njobs}%16",
    )
    minimap2_merge_template = jinja_env.get_template(
        "woltka_minimap2_merge.sbatch")
    minimap2_merge_script = _write_slurm(
        f'{out_dir}/merge',
        minimap2_merge_template,
        conda_environment=CONDA_ENV,
        output=out_dir,
        job_name=f"me_{job_id}",
        node_count=1,
        nprocs=16,
        wall_time_limit=MAX_WALL_1000,
        mem_in_gb=120
    )

    return minimap2_script, minimap2_merge_script
