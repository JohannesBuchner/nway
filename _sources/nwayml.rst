.. _install:
.. highlight:: shell

================
Machine learning
================

NWAY can be combined with machine learning :cite:p:`Salvato2022`.

For this, we reuse the ``--mag CAT:COL prior_file.txt`` option, which allows
specifying for each source, whether it is likely a target or not.
The new bit is, that we ingest prior probabilities, from
some machine learning classifier.

To use this, prepare a catalog column which gives the probability.
Let's say your catalog is called *OPT* and your column is called *prob*.
Then add to your nway command::

   ``--mag CAT:prob prob_prior_table.txt``

Here prob_prior_table.txt needs to map *prob* to probability,
which is a 1:1 mapping. The "dummy" file prob_prior_table.txt is
provided for download at: `https://raw.githubusercontent.com/JohannesBuchner/nway/refs/heads/master/doc/prob_prior_table.txt`.
