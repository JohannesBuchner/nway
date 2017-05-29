nway catalogue matcher
======================================

nway is a catalogue association program for arbitrarily many catalogues. 

Features
----------

* Bayesian match probabilities based on astronomical sky coordinates (RA, DEC)
* Works with arbitrarily many catalogues
* Can handle varying errors
* Can incorporate additional prior information, such as the magnitude or color distributions of the sources to match
* Works accurately and fast in small areas and all-sky catalogues

Documentation
---------------

See the manual: doc/nway-manual.pdf

A full worked example with demo catalogues is included.

Installation
---------------

The easiest way to install nway is with pip: 

	$ pip install nway

Add --user if you do not have admin privileges. Alternatively, 
download the source code and run "python setup.py install".

nway uses the following python packages:

* numpy, scipy, matplotlib
* astropy
* progressbar
* argparse
* joblib
* healpy

nway should work both with Python 3 and Python 2. 

Reporting issues
-----------------

Please file an issue on Github if you have any problems.

Citing nway correctly
----------------------

Please cite Salvato et al (in prep.).

Licence
---------------

This program is free and open-source software, 
licensed under AGPLv3 (see LICENSE and COPYING files).
Written by Johannes Buchner (C) 2013-2017





