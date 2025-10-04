# -----------------------------------------------------------------------------
# Copyright (c) 2025--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from qiita_client import QiitaPlugin, QiitaCommand
from .qp_pacbio import pacbio_processing

__version__ = "2025.9"

plugin = QiitaPlugin("qp-pacbio", __version__, "PacBio Processing")

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
