# -----------------------------------------------------------------------------
# Copyright (c) 2025--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from glob import glob
from os import environ, makedirs, mkdir
from os.path import basename, exists, join
from shutil import copy2
from subprocess import run

import pandas as pd
import yaml
from biom import load_table
from jinja2 import Environment
from pysyndna import (
    fit_linear_regression_models_for_qiita,
)
from qiita_client import ArtifactInfo

from .util import KISSLoader, find_base_path

JENV = Environment(loader=KISSLoader("data/templates"))
JGT = JENV.get_template
RESOURCES = yaml.safe_load(open(join(find_base_path(), "data/resources.yaml")))
CONDA_ENVIRONMENT = environ["ENVIRONMENT"]
PACBIO_PROCESSING_STEPS = 7


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


def _write_slurm(path, template, **ctx):
    makedirs(path, exist_ok=True)
    makedirs(join(path, "logs"), exist_ok=True)
    out_fp = join(path, f"{path.rsplit('/', 1)[-1]}.slurm")
    rendered = template.render(**ctx)
    with open(out_fp, mode="w", encoding="utf-8") as f:
        f.write(rendered)

    return out_fp


def pacbio_generate_templates(out_dir, job_id, njobs, result_fp, url):
    """Generate Slurm submission templates for PacBio processing.

    Parameters
    ----------
    out_dir : str
        Path to the job's output directory.
    job_id : str
        Qiita job id.
    njobs : int
        Number of array tasks/jobs.
    result_fp : str, filepath
        Folder where the final results will be stored.
    url : str
        URL to update the status of the jobs
    """
    resources = RESOURCES["PacBio processing"]

    main_parameters = {
        "conda_environment": CONDA_ENVIRONMENT,
        "output": out_dir,
        "url": url,
        "qjid": job_id,
        "result_fp": result_fp,
    }

    # Step 1
    template1 = JGT("1.hifiasm-meta_new.sbatch")
    step_resources = resources["step-1"]
    params = main_parameters | {
        "job_name": f"s1-{job_id}",
        "node_count": step_resources["node_count"],
        "nprocs": step_resources["nprocs"],
        "wall_time_limit": step_resources["wall_time_limit"],
        "mem_in_gb": step_resources["mem_in_gb"],
        "array_params": f"1-{njobs}%{step_resources['max_tasks']}",
    }
    _write_slurm(join(out_dir, "step-1"), template1, **params)

    # Step 2
    template2 = JGT("2.get-circular-genomes.sbatch")
    step_resources = resources["step-2"]
    params = main_parameters | {
        "job_name": f"s2-{job_id}",
        "node_count": step_resources["node_count"],
        "nprocs": step_resources["nprocs"],
        "wall_time_limit": step_resources["wall_time_limit"],
        "mem_in_gb": step_resources["mem_in_gb"],
        "array_params": f"1-{njobs}%{step_resources['max_tasks']}",
    }
    _write_slurm(join(out_dir, "step-2"), template2, **params)

    # Step 3
    template3 = JGT("3.minimap2_assembly.sbatch")
    step_resources = resources["step-3"]
    params = main_parameters | {
        "job_name": f"s3-{job_id}",
        "node_count": step_resources["node_count"],
        "nprocs": step_resources["nprocs"],
        "wall_time_limit": step_resources["wall_time_limit"],
        "mem_in_gb": step_resources["mem_in_gb"],
        "array_params": f"1-{njobs}%{step_resources['max_tasks']}",
    }
    _write_slurm(join(out_dir, "step-3"), template3, **params)

    # Step 4
    template4 = JGT("4.metawrap_binning_new.sbatch")
    step_resources = resources["step-4"]
    params = main_parameters | {
        "job_name": f"s4-{job_id}",
        "node_count": step_resources["node_count"],
        "nprocs": step_resources["nprocs"],
        "wall_time_limit": step_resources["wall_time_limit"],
        "mem_in_gb": step_resources["mem_in_gb"],
        "array_params": f"1-{njobs}%{step_resources['max_tasks']}",
    }
    _write_slurm(join(out_dir, "step-4"), template4, **params)

    # Step 5 (long filename isolated in T5_NAME)
    template5 = JGT("5.DAS_Tools_prepare_batch3_test.sbatch")
    step_resources = resources["step-5"]
    params = main_parameters | {
        "job_name": f"s5-{job_id}",
        "node_count": step_resources["node_count"],
        "nprocs": step_resources["nprocs"],
        "wall_time_limit": step_resources["wall_time_limit"],
        "mem_in_gb": step_resources["mem_in_gb"],
        "array_params": f"1-{njobs}%{step_resources['max_tasks']}",
    }
    _write_slurm(join(out_dir, "step-5"), template5, **params)

    # Step 6
    template6 = JGT("6.MAG_rename.sbatch")
    step_resources = resources["step-6"]
    params = main_parameters | {
        "job_name": f"s6-{job_id}",
        "node_count": step_resources["node_count"],
        "nprocs": step_resources["nprocs"],
        "wall_time_limit": step_resources["wall_time_limit"],
        "mem_in_gb": step_resources["mem_in_gb"],
        "array_params": f"1-{njobs}%{step_resources['max_tasks']}",
    }
    _write_slurm(join(out_dir, "step-6"), template6, **params)

    # Step 7
    template7 = JGT("7.checkm_batch.sbatch")
    step_resources = resources["step-7"]
    params = main_parameters | {
        "job_name": f"s7-{job_id}",
        "node_count": step_resources["node_count"],
        "nprocs": step_resources["nprocs"],
        "wall_time_limit": step_resources["wall_time_limit"],
        "mem_in_gb": step_resources["mem_in_gb"],
        "array_params": f"1-{njobs}%{step_resources['max_tasks']}",
    }
    _write_slurm(join(out_dir, "step-7"), template7, **params)

    # finish command - letting qiita know that we are done
    finish = JGT("finish.pacbio.processing.sbatch")
    step_resources = resources["finish"]
    params = main_parameters | {
        "job_name": f"f-{job_id}",
        "node_count": step_resources["node_count"],
        "nprocs": step_resources["nprocs"],
        "wall_time_limit": step_resources["wall_time_limit"],
        "mem_in_gb": step_resources["mem_in_gb"],
    }
    _write_slurm(join(out_dir, "finish"), finish, **params)


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

    preparation_information = join(out_dir, "prep_info.tsv")
    prep.set_index("sample_name").to_csv(preparation_information, sep="\t")

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
    result_fp = join(out_dir, "results")
    makedirs(result_fp, exist_ok=True)

    qclient.update_job_step(job_id, "Commands finished")
    # validating that all steps are complete
    with open(join(out_dir, "sample_list.txt"), "r") as fp:
        expected = len(fp.readlines())

    failed_steps = []
    for step_id in range(1, PACBIO_PROCESSING_STEPS + 1, 1):
        obs = glob(join(out_dir, f"step-{step_id}", "completed_*.log"))
        if len(obs) != expected:
            failed_steps.append(str(step_id))
    if failed_steps:
        failed_steps = ", ".join(failed_steps)
        message = (
            "These steps have less than expected logs: "
            f" {failed_steps}; please email the admin"
        )
        return (False, [], message)

    paths = [(f"{result_fp}/", "directory")]
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
        fp_coverages = join(out_dir, "coverages.tgz")
        mkdir(dest)
        dest = join(dest, "coverages.tgz")
        copy2(fp_coverages, dest)

        return dest

    errors = []
    ainfo = []

    fp_biom = f"{out_dir}/none.biom"
    fp_alng = f"{out_dir}/alignment.tar"
    if exists(fp_biom) and exists(fp_alng):
        ainfo.append(
            ArtifactInfo(
                "Per genome Predictions",
                "BIOM",
                [
                    (fp_biom, "biom"),
                    (fp_alng, "log"),
                    (_coverage_copy(f"{out_dir}/none/"), "plain_text"),
                ],
            )
        )
    else:
        errors.append(
            "Table none/per-genome was not created, please contact "
            "qiita.help@gmail.com for more information"
        )

    bioms = [
        (f"{out_dir}/per-gene.biom", "per_gene", "Per gene Predictions"),
        (f"{out_dir}/ko.biom", "ko", "KEGG Ontology (KO)"),
        (f"{out_dir}/ec.biom", "ec", "KEGG Enzyme (EC)"),
        (f"{out_dir}/pathway.biom", "pathway", "KEGG Pathway"),
    ]

    for fb, fn, bn in bioms:
        if exists(fb):
            files = [(fb, "biom")]
            files.append((_coverage_copy(f"{out_dir}/{fn}/"), "plain_text"))
            ainfo.append(
                ArtifactInfo(
                    bn,
                    "BIOM",
                    files,
                )
            )
        else:
            errors.append(
                f'Table "{bn}" was not created, please contact '
                "qiita.help@gmail.com for more information"
            )

    if errors:
        return False, ainfo, "\n".join(errors)
    else:
        return True, ainfo, ""


def generate_minimap2_processing(qclient, job_id, out_dir, parameters, url):
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
    url : str
        URL to send the respose, finish the job.

    Returns
    -------
    str, str
        Returns the two filepaths of the slurm scripts
    """
    resources = RESOURCES["Woltka v0.1.7, minimap2"]
    main_parameters = {
        "conda_environment": CONDA_ENVIRONMENT,
        "output": out_dir,
        "qjid": job_id,
    }

    qclient.update_job_step(
        job_id, "Step 1 of 4: Collecting info and generating submission"
    )

    artifact_id = parameters["artifact"]

    njobs = generate_sample_list(qclient, artifact_id, out_dir)

    qclient.update_job_step(
        job_id,
        "Step 2 of 4: Creating submission templates",
    )

    m2t = JGT("woltka_minimap2.sbatch")
    step_resources = resources["minimap2"]
    params = main_parameters | {
        "job_name": f"m2_{job_id}",
        "node_count": step_resources["node_count"],
        "nprocs": step_resources["nprocs"],
        "wall_time_limit": step_resources["wall_time_limit"],
        "mem_in_gb": step_resources["mem_in_gb"],
        "array_params": f"1-{njobs}%{step_resources['max_tasks']}",
    }
    minimap2_script = _write_slurm(f"{out_dir}/minimap2", m2t, **params)

    m2mt = JGT("woltka_minimap2_merge.sbatch")
    step_resources = resources["merge"]
    params = main_parameters | {
        "job_name": f"me_{job_id}",
        "node_count": step_resources["node_count"],
        "nprocs": step_resources["nprocs"],
        "wall_time_limit": step_resources["wall_time_limit"],
        "mem_in_gb": step_resources["mem_in_gb"],
        "url": url,
    }
    step_resources = resources["merge"]
    minimap2_merge_script = _write_slurm(f"{out_dir}/merge", m2mt, **params)

    return minimap2_script, minimap2_merge_script


def syndna_processing(qclient, job_id, parameters, out_dir):
    """generates output for syndna processing.

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

    errors = []
    ainfo = []

    failures = glob(f"{out_dir}/failed_*.log")
    if failures:
        errors.append("Samples failed: ")
        for f in failures:
            with open(f, "r") as fp:
                errors.append(fp.read())
        return False, ainfo, "\n".join(errors)

    completed = len(glob(f"{out_dir}/completed_*.log"))
    with open(f"{out_dir}/sample_list.txt") as fp:
        samples = len(fp.readlines())

    if completed != samples:
        errors.append(f"There are {samples - completed} missing samples.")

    fp_biom = f"{out_dir}/syndna/syndna.biom"
    # do we need to stor alignments?
    # fp_alng = f'{out_dir}/sams/final/alignment.tar'
    if not errors and exists(fp_biom):  # and exists(fp_alng):
        # if we got to this point a preparation file should exist in
        # the output folder
        prep = pd.read_csv(
            f"{out_dir}/prep_info.tsv", index_col=None, sep="\t", dtype=str
        )
        output = fit_linear_regression_models_for_qiita(
            prep, load_table(fp_biom), int(parameters["min_sample_counts"])
        )
        # saving results to disk
        lin_regress_results_fp = f"{out_dir}/lin_regress_by_sample_id.yaml"
        fit_syndna_models_log_fp = f"{out_dir}/fit_syndna_models_log.txt"
        with open(lin_regress_results_fp, "w") as fp:
            fp.write(output["lin_regress_by_sample_id"])
        with open(fit_syndna_models_log_fp, "w") as fp:
            fp.write(output["fit_syndna_models_log"])
        ainfo = [
            ArtifactInfo(
                "SynDNA hits",
                "BIOM",
                [
                    (fp_biom, "biom"),
                    # rm if fp_alng is not needed
                    # (fp_alng, "log"),
                    (lin_regress_results_fp, "log"),
                    (fit_syndna_models_log_fp, "log"),
                ],
            )
        ]
    else:
        ainfo = []
        errors.append(
            'Missing files from the "SynDNA hits"; please '
            "contact qiita.help@gmail.com for more information"
        )

    fp_seqs = f"{out_dir}/syndna/filtered"
    reads = []
    for f in glob(f"{fp_seqs}/*.fastq.gz"):
        reads.append((f, "raw_forward_seqs"))

    if not errors:
        ainfo.append(ArtifactInfo("reads without SynDNA", "per_sample_FASTQ", reads))
    else:
        return False, ainfo, "\n".join(errors)

    return True, ainfo, ""


def generate_syndna_processing(qclient, job_id, out_dir, parameters, url):
    """generates slurm scripts for syndna processing.

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
    url : str
        URL to send the respose, finish the job.

    Returns
    -------
    str, str
        Returns the two filepaths of the slurm scripts
    """
    resources = RESOURCES[
        "Remove SynDNA plasmid, insert, & GCF_000184185 reads (minimap2)"
    ]
    main_parameters = {
        "conda_environment": CONDA_ENVIRONMENT,
        "output": out_dir,
        "qjid": job_id,
    }

    qclient.update_job_step(
        job_id, "Step 1 of 4: Collecting info and generating submission"
    )

    artifact_id = parameters["artifact"]

    njobs = generate_sample_list(qclient, artifact_id, out_dir)

    qclient.update_job_step(
        job_id,
        "Step 2 of 4: Creating submission templates",
    )

    m2t = JGT("syndna.sbatch")
    step_resources = resources["syndna"]
    params = main_parameters | {
        "job_name": f"sd_{job_id}",
        "node_count": step_resources["node_count"],
        "nprocs": step_resources["nprocs"],
        "wall_time_limit": step_resources["wall_time_limit"],
        "mem_in_gb": step_resources["mem_in_gb"],
        "array_params": f"1-{njobs}%{step_resources['max_tasks']}",
    }
    minimap2_script = _write_slurm(f"{out_dir}/syndna", m2t, **params)

    m2mt = JGT("syndna_finish.sbatch")
    step_resources = resources["finish"]
    params = main_parameters | {
        "job_name": f"me_{job_id}",
        "node_count": step_resources["node_count"],
        "nprocs": step_resources["nprocs"],
        "wall_time_limit": step_resources["wall_time_limit"],
        "mem_in_gb": step_resources["mem_in_gb"],
        "url": url,
    }
    minimap2_finish_script = _write_slurm(f"{out_dir}/finish", m2mt, **params)

    return minimap2_script, minimap2_finish_script


def pacbio_apdater_removal(qclient, job_id, parameters, out_dir):
    """generates output for syndna processing.

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

    errors = []
    ainfo = []

    failures = glob(f"{out_dir}/failed_*.log")
    if failures:
        errors.append("Samples failed: ")
        for f in failures:
            with open(f, "r") as fp:
                errors.append(fp.read())
        return False, ainfo, "\n".join(errors)

    completed = len(glob(f"{out_dir}/completed_*.log"))
    with open(f"{out_dir}/sample_list.txt") as fp:
        samples = len(fp.readlines())

    if completed != samples:
        errors.append(f"There are {samples - completed} missing samples.")

    fp_seqs = f"{out_dir}/processing/filtered"
    reads = []
    for f in glob(f"{fp_seqs}/*.fastq.gz"):
        reads.append((f, "raw_forward_seqs"))

    if not errors:
        ainfo.append(ArtifactInfo("no adapter reads", "per_sample_FASTQ", reads))
    else:
        return False, ainfo, "\n".join(errors)

    return True, ainfo, ""


def generate_pacbio_apdater_removal(qclient, job_id, out_dir, parameters, url):
    """generates slurm scripts for pacbio adapter removal.

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
    url : str
        URL to send the respose, finish the job.

    Returns
    -------
    str, str
        Returns the two filepaths of the slurm scripts
    """
    resources = RESOURCES["PacBio adapter removal via lima/pbmarkdup"]
    main_parameters = {
        "conda_environment": CONDA_ENVIRONMENT,
        "output": out_dir,
        "qjid": job_id,
    }

    qclient.update_job_step(
        job_id, "Step 1 of 4: Collecting info and generating submission"
    )

    artifact_id = parameters["artifact"]

    njobs = generate_sample_list(qclient, artifact_id, out_dir)

    qclient.update_job_step(
        job_id,
        "Step 2 of 4: Creating submission templates",
    )

    lima_cmd = f'lima "${{filename}}" {out_dir}/adapter.fasta "${{fout}}.fastq.gz" --hifi-preset SYMMETRIC --peek-guess'
    if parameters["css"] and parameters["css"] != "False":
        lima_cmd += " --css"
    for k, v in parameters.items():
        if k in {
            "min-score",
            "min-end-score",
            "min-ref-span",
            "min-scoring-regions",
            "min-score-lead",
            "min-length",
            "window-size-multi",
        }:
            lima_cmd += f" --{k} {v}"

    template = JGT("pacbio_adapter.sbatch")
    step_resources = resources["processing"]
    params = main_parameters | {
        "job_name": f"sd_{job_id}",
        "node_count": step_resources["node_count"],
        "nprocs": step_resources["nprocs"],
        "wall_time_limit": step_resources["wall_time_limit"],
        "mem_in_gb": step_resources["mem_in_gb"],
        "array_params": f"1-{njobs}%{step_resources['max_tasks']}",
        "lima_cmd": lima_cmd,
    }
    processing_script = _write_slurm(f"{out_dir}/processing", template, **params)

    template = JGT("pacbio_adapter_finish.sbatch")
    step_resources = resources["finish"]
    params = main_parameters | {
        "job_name": f"me_{job_id}",
        "node_count": step_resources["node_count"],
        "nprocs": step_resources["nprocs"],
        "wall_time_limit": step_resources["wall_time_limit"],
        "mem_in_gb": step_resources["mem_in_gb"],
        "url": url,
    }
    finish_script = _write_slurm(f"{out_dir}/finish", template, **params)

    with open(f"{out_dir}/adapter.fasta", "w") as fp:
        for i, adapter in enumerate(parameters["adapter"].split(",")):
            i += 1
            fp.write(f">adapter_{i}\n")
            fp.write(f"{adapter}\n")

    return processing_script, finish_script
