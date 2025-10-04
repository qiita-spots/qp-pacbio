# -----------------------------------------------------------------------------
# Copyright (c) 2025--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from qiita_client import ArtifactInfo
from jinja2 import Environment, BaseLoader, TemplateNotFound
from os.path import basename, join, exists, getmtime
from os import makedirs
import pathlib


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


# taken from https://jinja.palletsprojects.com/en/3.0.x/api/#jinja2.BaseLoader
class KISSLoader(BaseLoader):
    def __init__(self, path):
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
    """Generates the templates to create slurm submission jobs for
    pacbio processing steps.

    Parameters
    ----------
    out_dir : str
        The path to the job's output directory
    job_id : str
        The job id
    njobs : int
        The number of total jobs
    """

    jinja_env = Environment(loader=KISSLoader("../data/templates"))

    template0 = jinja_env.get_template("0.mapping_minimap2_db.sbatch")
    cdir0 = f"{out_dir}/step-0"
    makedirs(cdir0)
    makedirs(f"{cdir0}/logs", exist_ok=True)
    with open(f"{cdir0}/step-0.slurm", mode="w", encoding="utf-8") as f:
        f.write(
            template0.render(
                conda_environment="qp_pacbio_2025.9",
                output=f"{out_dir}",
                job_name=f"s0-{job_id}",
                node_count=1,
                nprocs=16,
                wall_time_limit=1000,
                mem_in_gb=300,
                array_params=f"1-{njobs}%16",
            )
        )

    template = jinja_env.get_template("1.hifiasm-meta_new.sbatch")
    cdir = f"{out_dir}/step-1"
    makedirs(cdir)
    makedirs(f"{cdir}/logs", exist_ok=True)
    with open(f"{cdir}/step-1.slurm", mode="w", encoding="utf-8") as f:
        f.write(
            template.render(
                conda_environment="qp_pacbio_2025.9",
                output=f"{out_dir}",
                job_name=f"s1-{job_id}",
                node_count=1,
                nprocs=16,
                wall_time_limit=1000,
                mem_in_gb=300,
                array_params=f"1-{njobs}%16",
            )
        )

    template2 = jinja_env.get_template("2.get-circular-genomes.sbatch")
    cdir2 = f"{out_dir}/step-2"
    makedirs(cdir2)
    makedirs(f"{cdir2}/logs", exist_ok=True)
    with open(f"{cdir2}/step-2.slurm", mode="w", encoding="utf-8") as f:
        f.write(
            template2.render(
                conda_environment="qp_pacbio_2025.9",
                output=f"{out_dir}",
                job_name=f"s2-{job_id}",
                node_count=1,
                nprocs=1,
                wall_time_limit=500,
                mem_in_gb=16,
                array_params=f"1-{njobs}%16",
            )
        )

    template3 = jinja_env.get_template("3.minimap2_assembly.sbatch")
    cdir3 = f"{out_dir}/step-3"
    makedirs(cdir3)
    makedirs(f"{cdir3}/logs", exist_ok=True)
    with open(f"{cdir3}/step-3.slurm", mode="w", encoding="utf-8") as f:
        f.write(
            template3.render(
                conda_environment="qp_pacbio_2025.9",
                output=f"{out_dir}",
                job_name=f"s3-{job_id}",
                node_count=1,
                nprocs=8,
                wall_time_limit=500,
                mem_in_gb=50,
                array_params=f"1-{njobs}%16",
            )
        )

    template4 = jinja_env.get_template("4.metawrap_binning_new.sbatch")
    cdir4 = f"{out_dir}/step-4"
    makedirs(cdir4)
    makedirs(f"{cdir4}/logs", exist_ok=True)
    with open(f"{cdir4}/step-4.slurm", mode="w", encoding="utf-8") as f:
        f.write(
            template4.render(
                conda_environment="qp_pacbio_2025.9",
                output=f"{out_dir}",
                job_name=f"s4-{job_id}",
                node_count=1,
                nprocs=8,
                wall_time_limit=500,
                mem_in_gb=50,
                array_params=f"1-{njobs}%16",
            )
        )

    template5 = jinja_env.get_template("5.DAS_Tools_prepare_batch3_test.sbatch")
    cdir5 = f"{out_dir}/step-5"
    makedirs(cdir5)
    makedirs(f"{cdir5}/logs", exist_ok=True)
    with open(f"{cdir5}/step-5.slurm", mode="w", encoding="utf-8") as f:
        f.write(
            template5.render(
                conda_environment="qp_pacbio_2025.9",
                output=f"{out_dir}",
                job_name=f"s5-{job_id}",
                node_count=1,
                nprocs=8,
                wall_time_limit=500,
                mem_in_gb=50,
                array_params=f"1-{njobs}%16",
            )
        )

    template6 = jinja_env.get_template("6.MAG_rename.sbatch")
    cdir6 = f"{out_dir}/step-6"
    makedirs(cdir6)
    makedirs(f"{cdir6}/logs", exist_ok=True)
    with open(f"{cdir6}/step-6.slurm", mode="w", encoding="utf-8") as f:
        f.write(
            template6.render(
                conda_environment="qp_pacbio_2025.9",
                output=f"{out_dir}",
                job_name=f"s6-{job_id}",
                node_count=1,
                nprocs=8,
                wall_time_limit=500,
                mem_in_gb=50,
                array_params=f"1-{njobs}%16",
            )
        )

    template7 = jinja_env.get_template("7.checkm_batch.sbatch")
    cdir7 = f"{out_dir}/step-7"
    makedirs(cdir7)
    makedirs(f"{cdir7}/logs", exist_ok=True)
    with open(f"{cdir7}/step-7.slurm", mode="w", encoding="utf-8") as f:
        f.write(
            template7.render(
                conda_environment="qp_pacbio_2025.9",
                output=f"{out_dir}",
                job_name=f"s7-{job_id}",
                node_count=1,
                nprocs=8,
                wall_time_limit=500,
                mem_in_gb=50,
                array_params=f"1-{njobs}%16",
            )
        )


def generate_sample_list(qclient, artifact_id, out_dir):
    """Generates the sample_list.txt file of sample names/files to be
    processed with slurm templates.

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

    lookup = prep.set_index("run_prefix")["sample_name"].to_dict()
    lines = []
    for _, (fwd, _) in files.items():
        fwd = fwd["filepath"]
        sname = search_by_filename(basename(fwd), lookup)
        lines.append(f"{sname}\t{fwd}")

    with open(f"{out_dir}/sample_list.txt", "w") as f:
        f.write("\n".join(lines))

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
    results_fp = f"{out_dir}/results"
    makedirs(results_fp)

    qclient.update_job_step(
        job_id, "Step 1 of 3: Collecting info and generating submission"
    )
    artifact_id = parameters["artifact_id"]

    njobs = generate_sample_list(qclient, artifact_id, out_dir)

    qclient.update_job_step(job_id, "Step 2 of 3: Creating submission templates")

    generate_templates(out_dir, job_id, njobs)

    # Submit Step 0 (independent)
    # jid0 = qclient.submit_job(f'{out_dir}/step-0/step-0.slurm')

    # Submit Step 1 (independent) and capture its Slurm jobid
    # jid1 = qclient.submit_job(f'{out_dir}/step-1/step-1.slurm')

    # Re-render Step 2 with dependency on the corresponding task of Step 1
    # (Overwrites step-2/step-2.slurm to inject --dependency)
    jinja_env = Environment(loader=KISSLoader("../data/templates"))
    template2 = jinja_env.get_template("2.get-circular-genomes.sbatch")
    cdir2 = f"{out_dir}/step-2"
    makedirs(cdir2, exist_ok=True)
    with open(f"{cdir2}/step-2.slurm", mode="w", encoding="utf-8") as f:
        f.write(
            template2.render(
                conda_environment="qp_pacbio_2025.9",
                output=f"{out_dir}",
                job_name=f"s2-{job_id}",
                node_count=1,
                nprocs=1,
                wall_time_limit=500,
                mem_in_gb=16,
                array_params=f"1:{njobs}%16",
                dependency=(f"aftercorr:{jid1}" if jid1 else None),
            )
        )

    # Submit Step 2 (now chained to Step 1 with aftercorr)
    # jid2 = qclient.submit_job(f'{out_dir}/step-2/step-2.slurm')

    qclient.update_job_step(job_id, "Step 3 of 3: Running commands")

    # TODO 3. wait for all jobs to finish

    qclient.update_job_step(job_id, "Commands finished")

    paths = [(f"{results_fp}/", "directory")]
    return (
        True,
        [ArtifactInfo("output", "job-output-folder", paths)],
        f"{job_id} completed.",
    )
