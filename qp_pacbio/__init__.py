# -----------------------------------------------------------------------------
# Copyright (c) 2025--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from qiita_client import QiitaCommand, QiitaPlugin

from .qp_pacbio import (
    feature_table_generation,
    minimap2_processing,
    pacbio_adapter_removal,
    pacbio_processing,
    syndna_processing,
)
from .util import plugin_details

plugin = QiitaPlugin(**plugin_details)

#
# minimap2 command
#

req_params = {"artifact": ("artifact", ["per_sample_FASTQ"])}
opt_params = {
    # ToDo, fix in the next opportunity: https://github.com/qiita-spots/qp-pacbio/issues/25
    "percent-identity": ("float", "0.9"),
    "percent-coverage": ("float", "0.9"),
    "Database": ['choice:["WoLr2"]', "WoLr2"],
}
outputs = {
    # taxonomic
    "Per genome Predictions": "BIOM",
    "Per gene Predictions": "BIOM",
    # functional
    "KEGG Ontology (KO)": "BIOM",
    "KEGG Enzyme (EC)": "BIOM",
    "KEGG Pathway": "BIOM",
}
dflt_param_set = {
    "WoLr2": {
        "Database": "WoLr2",
        "percent-identity": 0.9,
        "percent-coverage": 0.9,
    },
}
minimap2_cmd = QiitaCommand(
    "Woltka v0.1.7 with cov and id filter",
    "Functional and Taxonomic Predictions",
    minimap2_processing,
    req_params,
    opt_params,
    outputs,
    dflt_param_set,
)
plugin.register_command(minimap2_cmd)

#
# pacbio processing pipeline command
#

req_params = {
    "artifact": ("artifact", ["per_sample_FASTQ"]),
}
opt_params = {"Processing": ['choice:["default"]', "default"]}
outputs = {"output": "job-output-folder"}
dflt_param_set = {"default": {"Processing": "default"}}
pacbio_processing_cmd = QiitaCommand(
    "PacBio processing",
    "Default PacBio processing for Metagenomic Data",
    pacbio_processing,
    req_params,
    opt_params,
    outputs,
    dflt_param_set,
)

plugin.register_command(pacbio_processing_cmd)

#
# syndna filtering
#

req_params = {
    "artifact": ("artifact", ["per_sample_FASTQ"]),
}
opt_params = {"min_sample_counts": ("integer", "1")}
outputs = {
    "SynDNA hits": "BIOM",
    "reads without SynDNA": "per_sample_FASTQ",
}
dflt_param_set = {
    "SynDNA": {"min_sample_counts": 1},
}
pacbio_processing_cmd = QiitaCommand(
    "Remove SynDNA plasmid, insert, & GCF_000184185 reads (minimap2)",
    "Remove SynDNA reads using minimap2 and woltka",
    syndna_processing,
    req_params,
    opt_params,
    outputs,
    dflt_param_set,
)
plugin.register_command(pacbio_processing_cmd)

#
# pacbio adapter filtering
#

req_params = {
    "artifact": ("artifact", ["per_sample_FASTQ"]),
}
opt_params = {
    # NOTE ABOUT adapter_sets: this is a comma-separated string.  Its contents can be either (or both):
    # - a comma-separated list of sequences to filter
    # - a comma-separated list of filenames to use; these should live in qp_pacbio/data/adapters/
    # the code will "merge" all these options
    "adapter_sets": ("string", ""),
    "neighbors": ("boolean", "False"),
    "peek-guess": ("boolean", "False"),
    "hifi-preset": ('choice: ["SYMMETRIC", "ASYMMETRIC"]', "ASYMMETRIC"),
}
outputs = {
    "no adapter reads": "per_sample_FASTQ",
}
dflt_param_set = {
    "PacBio adapter": {
        "adapter_sets": "AAGCAGTGGTATCAACGCAGAGTACT",
        "neighbors": False,
        "peek-guess": True,
        "hifi-preset": "SYMMETRIC",
    },
    "PacBio twist adapter": {
        "adapter_sets": "twist_adapters_231010.fasta.gz",
        "neighbors": True,
        "peek-guess": True,
        "hifi-preset": "ASYMMETRIC",
    },
}

pacbio_adapter_removal_cmd = QiitaCommand(
    "Adapter removal via lima/pbmarkdup v2.13",
    "Remove adapter reads using lima/pbmarkdup",
    pacbio_adapter_removal,
    req_params,
    opt_params,
    outputs,
    dflt_param_set,
)
plugin.register_command(pacbio_adapter_removal_cmd)

#
# feature table generation
#

req_params = {
    "analysis": ("integer", "None"),
}
opt_params = {
    "artifacts": ("mchoice:[]", "[]"),
    # ToDo, fix in the next opportunity: https://github.com/qiita-spots/qp-pacbio/issues/25
    "percent-identity": ("float", "0.995"),
    "GToTree-c": ("float", "0.4"),
    "GToTree-G": ("float", "0.4"),
}
outputs = {
    "Merged LCG/MAG feature table": "BIOM",
}
dflt_param_set = {
    "Default": {
        "percent-identity": 0.995,
        "GToTree-c": 0.4,
        "GToTree-G": 0.4,
    },
}
ft_cmd = QiitaCommand(
    "Feature Table from LCG/MAG",
    "Feature Table Generation from LCG/MAG",
    feature_table_generation,
    req_params,
    opt_params,
    outputs,
    dflt_param_set,
)
plugin.register_command(ft_cmd)
