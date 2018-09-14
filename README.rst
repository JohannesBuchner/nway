nway
======================================

nway is a source cross-matching tool for arbitrarily many astronomical catalogues. 

Features
----------

* Bayesian match probabilities based on astronomical sky coordinates (RA, DEC)
* Works with arbitrarily many catalogues
* Can handle varying errors
* Can incorporate additional prior information, such as the magnitude or color distributions of the sources to match
* Works accurately and fast in small areas and all-sky catalogues

Documentation
---------------

See the manual in the doc/ folder, https://github.com/JohannesBuchner/nway/raw/master/doc/nway-manual.pdf

A full worked example with demo catalogues is included.

Installation
---------------

The easiest way to install nway is with pip::

	$ pip install nway

Add --user if you do not have admin privileges. Alternatively, 
download the source code and run "python setup.py install".

nway uses the following python packages:

* numpy, scipy, matplotlib
* astropy
* progressbar2 or progressbar or progressbar-latest
* joblib
* healpy
* pandas

nway works with both Python 3 and Python 2 and various astropy versions.

.. image:: https://codecov.io/gh/JohannesBuchner/nway/branch/master/graph/badge.svg
	:target: https://codecov.io/gh/JohannesBuchner/nway
.. image:: https://travis-ci.org/JohannesBuchner/nway.svg?branch=master
	:target: https://travis-ci.org/JohannesBuchner/nway

Reporting issues
-----------------

Please file an issue on Github if you have any problems.

Citing nway correctly
----------------------

Please cite Salvato et al (2017). https://arxiv.org/abs/1705.10711 (`Get BibTex <http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2017arXiv170510711S&data_type=BIBTEX&db_key=PRE&nocookieset=1>`_)

Licence
---------------

This program is free and open-source software, 
licensed under AGPLv3 (see LICENSE and COPYING files).
Written by Johannes Buchner (C) 2013-2017





