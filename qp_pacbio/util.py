# -----------------------------------------------------------------------------
# Copyright (c) 2025--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from os import environ
from os.path import join, expanduser, getmtime, exists
from configparser import ConfigParser
import pathlib
from jinja2 import BaseLoader, TemplateNotFound

from qiita_client import QiitaClient


plugin_details = {
    "name": "qp-pacbio",
    "version": "2025.09",
    "description": "PacBio processing",
}


def client_connect(url):
    name = plugin_details["name"]
    version = plugin_details["version"]

    config = ConfigParser()
    qp = join(expanduser("~"), ".qiita_plugins")
    conf_dir = environ.get("QIITA_PLUGINS_DIR", qp)
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
