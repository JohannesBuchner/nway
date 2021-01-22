---
title: 'NWAY - probabilistic matching for associating N catalogs'
tags:
  - Catalogues
  - Python
  - Bayesian inference
  - Combinatorics
authors:
  - name: Johannes Buchner
    orcid: 0000-0003-0872-7098
    affiliation: "1"
affiliations:
 - name: Max Planck Institute for Extraterrestrial Physics, Giessenbachstrasse, 85741 Garching, Germany. 
   index: 1

date: 22 January 2021
bibliography: paper.bib

---

# Summary

NWAY matches an arbitrary number of astrophysical catalogs,
and computes the probability of each reasonably possible linkage combination.
It provides tools to make reliable cuts to identify secure and ambiguous
associations. 
NWAY allows sources to be absent in some catalogs, and allows 
incorporating information from additional data columns that help differentiate
the target source class from contaminating field sources.

# Statement of need

In astrophysics, a common task is to assemble multi-wavelength 
information about individual sources. This is done by taking the detections of
sources in the sky (positions, errors, and fluxes/magnitudes) from a 
catalogue of one wavelength and matching it to another from another 
wavelength, or multiple such catalogues. Care has to be taken to consider all
possible matches and also the possibility that the source does not have a
counterpart in a catalogue of a given depth. For many classes of sources,
the Spectral Energy Distribution (SED) provides additional hints, which
associations are likely real. For instance, the color distribution of stars in
the WISE bands is different than that of quasars or galaxies.

NWAY is a generic and self-consistent solution to these tasks:

1. Matching of N catalogues simultaneously.
2. Consideration of all combinatorically possible matches.
3. Consideration of partial matches across catalogues, i.e. the absence
   of counterparts in some catalogues.
4. Taking into account the positional uncertainty.
5. Computation of a probability for each possible match.
6. Computation of a probability that there is no match.
7. Incorporating magnitude, color or other information about the sources
   of interest, refining the match probabilities.

NWAY is used extensively for large X-ray surveys including those of eROSITA,
Chandra, XMM-Newton and ROSAT, which sometimes feature substantial positional
uncertainties.

# Documentation

Included in NWAY is a [user manual](https://github.com/JohannesBuchner/nway/raw/master/doc/nway-manual.pdf) 
and [example catalogs](https://github.com/JohannesBuchner/nway/tree/master/doc/) to run its tutorials.

The method has been introduced and described in @Salvato2017.
This paper serves as the reference for the implementation.

Two interfaces are available: A command line tool and a Python library.

# Acknowledgements

NWAY is based on a collaboration of many people, including
Mara Salvato, Tamás Budavári, Sotiria Fotopoulou, Fabrizia
Guglielmetti, Arne Rau, Tom Dwelly, Andrea Merloni and Kirpal Nandra.

# References
