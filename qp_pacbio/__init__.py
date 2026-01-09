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
from .util import get_local_adapter_files, plugin_details

plugin = QiitaPlugin(**plugin_details)

#
# minimap2 command
#

req_params = {"artifact": ("artifact", ["per_sample_FASTQ"])}
opt_params = {
    "Database": ['choice:["WoLr2"]', "WoLr2"],
    "identity": (float, 0.9),
    "coverage": (float, 0.9),
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
    "WoLr2": {"Database": "WoLr2"},
    "identity": "0.9",
    "coverage": "0.9",
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
    "adapter_sets": ("string", "AAGCAGTGGTATCAACGCAGAGTACT"),
    "css": ("boolean", "False"),
    "min-score": ("integer", "0"),
    "min-end-score": ("integer", "0"),
    "min-ref-span": ("float", "0.5"),
    "min-scoring-regions": ("integer", "1"),
    "min-score-lead": ("integer", "10"),
    "min-length": ("integer", "50"),
    "window-size-multi": ("float", "3"),
}
outputs = {
    "no adapter reads": "per_sample_FASTQ",
}
dflt_param_set = {
    "PacBio adapter": {
        "adapter_sets": "AAGCAGTGGTATCAACGCAGAGTACT",
        "css": False,
        "min-score": 0,
        "min-end-score": 0,
        "min-ref-span": 0.5,
        "min-scoring-regions": 1,
        "min-score-lead": 10,
        "min-length": 50,
        "window-size-multi": 3,
    }
}
for name in get_local_adapter_files().keys():
    dflt_param_set[name] = {
        "adapter_sets": name,
        "css": False,
        "min-score": 0,
        "min-end-score": 0,
        "min-ref-span": 0.5,
        "min-scoring-regions": 1,
        "min-score-lead": 10,
        "min-length": 50,
        "window-size-multi": 3,
    }


pacbio_adapter_removal_cmd = QiitaCommand(
    "PacBio adapter removal via lima/pbmarkdup",
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
