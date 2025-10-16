# -----------------------------------------------------------------------------
# Copyright (c) 2025--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from os import environ
from os.path import join, expanduser
from configparser import ConfigParser

from qiita_client import QiitaClient


plugin_details = {'name': 'qp-pacbio',
                  'version': '2025.9',
                  'description': 'PacBio processing'}


def client_connect(url):
    name = plugin_details["name"]
    version = plugin_details["version"]

    config = ConfigParser()
    conf_dir = environ.get(
        "QIITA_PLUGINS_DIR", join(expanduser("~"), ".qiita_plugins")
    )
    conf_fp = join(conf_dir, f"{name}_{version}.conf")

    with open(conf_fp, "r", encoding="utf-8") as conf_file:
        config.read_file(conf_file)

    qclient = QiitaClient(
        url,
        config.get("oauth2", "CLIENT_ID"),
        config.get("oauth2", "CLIENT_SECRET"),
        config.get("oauth2", "SERVER_CERT"),
    )
    return qclient
