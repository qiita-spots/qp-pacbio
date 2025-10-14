# -----------------------------------------------------------------------------
# Copyright (c) 2025--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from qiita_client import QiitaPlugin, QiitaCommand
from .qp_pacbio import pacbio_processing, minimap2_processing
from .util import plugin_details

__version__ = "2025.9"

plugin = QiitaPlugin(**plugin_details)

#
# minimap2 command
#

req_params = {'input': ('artifact_id', ['per_sample_FASTQ'])}
opt_params = dict()
outputs = {
    # taxonomic
    'Per genome Predictions': 'BIOM',
    'Per gene Predictions': 'BIOM',
    # functional
    'KEGG Ontology (KO)': 'BIOM',
    'KEGG Enzyme (EC)': 'BIOM',
    'KEGG Pathway': 'BIOM',
    }
dflt_param_set = dict()

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
    "artifact_id": ("integer", [None]),
}
opt_params = dict()
outputs = {"output": "job-output-folder"}
dflt_param_set = dict()

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
