|Build Status| |Coverage Status|

Qiita plugin to process PacBio reads.

This plugin currently provides 2 commands for Qiita:

* **Woltka v0.1.7, minimap2**: which generates feature and functional profiles agains WoLr2;
  the expected output are BIOM artifacts

* **PacBio processing**: which goes from step 1 to 7 in the image below. The expected output
  is a main folder with folders per-sample and folders for each of the different outputs, as follows:

  * MAG folder: all Metagenome-Assembled Genome (MAG) generatedfor that sample
  * LCG folder: all Long-Circular Genome (LCG) generated for that sample that are over 512kb in size - approximate 515,000 bases (half a million)
  * small_LCG folder: all Long-Circular Genome (LCG) generated for that sample that are under 512kb in size
  * [sample-name].fna.gz: the no LCG reads used for MAG generation
  * [sample-name].checkm.txt.gz: MAG quanlity information from CheckM v1.2.3


.. image:: images/PacBioProcessing.png
   :alt: qp_pacbio processing steps
   :width: 90%
   :align: center


.. |Build Status| image:: https://github.com/qiita-spots/qp-pacbio/actions/workflows/qiita-plugin-ci.yml/badge.svg
   :target: https://github.com/qiita-spots/qp-pacbio/actions/workflows/qiita-plugin-ci.yml
.. |Coverage Status| image:: https://coveralls.io/repos/github/qiita-spots/qp-pacbio/badge.svg?branch=dev
   :target: https://coveralls.io/github/qiita-spots/qp-pacbio?branch=master
