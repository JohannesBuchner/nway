nway
====

nway is a source cross-matching tool for arbitrarily many astronomical catalogues. 

Why nway?
---------

In astrophysics, a common task is to assemble multi-wavelength
information about individual sources. This is done by taking the
detections of sources in the sky (positions, errors, and
fluxes/magnitudes) from a catalogue of one wavelength and matching it to
another from another wavelength, or multiple such catalogues. Care has
to be taken to consider all possible matches and also the possibility
that the source does not have a counterpart in a catalogue of a given
depth. For many classes of sources, the Spectral Energy Distribution
(SED) provides additional hints, which associations are likely real. For
instance, the color distribution of stars in the WISE bands is different
than that of quasars or galaxies.

nway is a generic solution to these tasks:

#. Matching of N catalogues simultaneously.

#. Consideration of all combinatorically possible matches.

#. Consideration of partial matches across catalogues, i.e. the absence
   of counterparts in some catalogues.

#. Taking into account the positional uncertainty of each source.

#. Computation of a probability for each possible match.

#. Computation of a probability that there is no match.

#. Incorporating additional prior information about the sources, such as
   magnitude, color, etc, to refine the match probabilities.

#. Accurate and fast for both small survey areas and all-sky catalogues


Contributors
------------

-  Mara Salvato – idea and leading the science case.

-  Tamás Budavári – shared basic implementation of his formulae.

-  Sotiria Fotopoulou – initial implementation for matching three
   catalogues.

-  Johannes Buchner – complete code rewrite for the general case,
   documentation, manual and adding features.

The code is the results of many discussion among colleagues and friends.
We thank in particular: Tamás Budavári, Sotiria Fotopoulou, Fabrizia
Guglielmetti, Arne Rau, Tom Dwelly, Andrea Merloni and Kirpal Nandra.


Citing nway correctly
----------------------

Please cite Salvato et al., MNRAS, 2018 (ArXiV:
https://arxiv.org/abs/1705.10711, ADS:
https://ui.adsabs.harvard.edu/abs/2018MNRAS.473.4937S/abstract).

and add a footnote with the URL to this repository::

	\footnote{\url{https://github.com/JohannesBuchner/nway/}}


Documentation
---------------

The documentation can be found at: https://johannesbuchner.github.io/nway/

The older manual is in the doc/ folder, https://github.com/JohannesBuchner/nway/raw/master/doc/nway-manual.pdf

A full worked example with demo catalogues is included.

Installation
------------

The easiest way to install nway is with pip::

	$ pip install nway

Add --user if you do not have admin privileges. Alternatively, 
download the source code and run "python setup.py install".

nway uses the following python packages:

* numpy, scipy, matplotlib
* astropy
* tqdm
* joblib
* healpy
* pandas

nway works with both Python 3 and Python 2 and various astropy versions.

.. image:: https://codecov.io/gh/JohannesBuchner/nway/branch/master/graph/badge.svg
	:target: https://codecov.io/gh/JohannesBuchner/nway
.. image:: https://github.com/JohannesBuchner/nway/actions/workflows/tests.yml/badge.svg
	:target: https://github.com/JohannesBuchner/nway/actions/workflows/tests.yml

Development, questions and issues
---------------------------------

Please let us know if you have comments about this manual or problems
running/installing nway.

For reporting bugs or requesting features, nway's issue tracker is at:
https://github.com/JohannesBuchner/nway/issues/

nway is a small but powerful tool. Most questions or apparent issues
arise to understand which information lead to a counterpart being
preferred, and therefore understanding the input data well. The manual
will guide you through the computation and the pieces of information
added together.
If you encounter issues, try simplifying.

Licence
---------------

This program is free and open-source software, 
licensed under AGPLv3 (see LICENSE and COPYING files).
Written by Johannes Buchner (C) 2013-2025
