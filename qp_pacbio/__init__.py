# -----------------------------------------------------------------------------
# Copyright (c) 2025--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from qiita_client import QiitaCommand, QiitaPlugin

from .qp_pacbio import (
    minimap2_processing,
    pacbio_processing,
    syndna_processing,
)
from .util import plugin_details

plugin = QiitaPlugin(**plugin_details)

#
# minimap2 command
#

req_params = {"artifact": ("artifact", ["per_sample_FASTQ"])}
opt_params = {"Database": ['choice:["WoLr2"]', "WoLr2"]}
outputs = {
    # taxonomic
    "Per genome Predictions": "BIOM",
    "Per gene Predictions": "BIOM",
    # functional
    "KEGG Ontology (KO)": "BIOM",
    "KEGG Enzyme (EC)": "BIOM",
    "KEGG Pathway": "BIOM",
}
dflt_param_set = {"WoLr2": {"Database": "WoLr2"}}
minimap2_cmd = QiitaCommand(
    "Woltka v0.1.7, minimap2",
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
