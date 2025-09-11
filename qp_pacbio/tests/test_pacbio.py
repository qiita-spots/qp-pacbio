# -----------------------------------------------------------------------------
# Copyright (c) 2020--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from qiita_client import ArtifactInfo
from qiita_client.testing import PluginTestCase
from os import remove
from os.path import exists, isdir
from tempfile import mkdtemp
from shutil import rmtree

from qp_pacbio import plugin
from qp_pacbio.qp_pacbio import pacbio_processing


class WoltkaTests(PluginTestCase):
    def setUp(self):
        plugin("https://localhost:21174", 'register', 'ignored')

        self._clean_up_files = []

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_pacbio_processing(self):
        params = {'artifact_id': 5}
        job_id = 'my-job-id'
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # this should fail cause we don't have valid data
        success, ainfo, msg = pacbio_processing(
            self.qclient, job_id, params, out_dir)

        with open(f'{out_dir}/sample_list.txt', 'r') as f:
            obs_lines = f.readlines()
        self.assertCountEqual([], obs_lines)

        # TODO 1. adding tests for template processing

        self.assertTrue(success)
        exp = [ArtifactInfo('output', 'job-output-folder',
               [(f'{out_dir}/results/', 'directory')])]
        self.assertCountEqual(ainfo, exp)


if __name__ == '__main__':
    main()
