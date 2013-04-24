3way catalogue matcher
======================================

3way is a catalogue association program for arbitrarily many catalogues. 
It is useful for astronomical sky coordinates (RA, DEC).

It uses 3 methods for matching:

1) matching using exclusive search radius
  
     A rough cutting to exclude very distant matches. 
     
     *--radius* defines the radius in arcsec
  
2) computing a probability of matches based on position errors and distances.
  
     For each combination found in (1), the model that all positions belong 
     to one source is compared to the model that they are all independent.
     
     *--prior* specifies the density of false positives. E.g. nx/(1887*15e+18)
  
3) optional: weighting based on properties of close objects versus distant objects 
  
     A magnitude histogram of input catalogue matches is compared to the 
     magnitude histogram of the full input catalogue.
     
     *--mag* specifies the catalogue column to use (can be used multiple times).
     
     *--mag-radius* specifies the search radius in arcsec.

The final catalogue (*--out*) contains all input catalogues, and additional separation and probability columns.
The catalogue can be trimmed using *--min-prob*.
A flag is added for the most probable match for the primary catalogue (1), and secondary, comparably good solutions (2, *--acceptable-prob*)

Input/output formats
---------------------
The input catalogues need to specially prepared, as FITS files. 
In each file, the first table is used, and its extension name/table name is very important! 
This name is used for column specification and as a prefix in the merged catalogue.
It can be set in Table Parameters -> Name in TOPCAT.
The columns for matching are automatically searched by name.

The final catalogue is also a FITS table (table "multimatch"). Some columns deserve explaination:

* bf: Bayes factor for this association being a single source. (ln)
* bfpost: Bayesian posterior for this association being a single source.
* bias_*: individual weightings added by magnitude histogramming (ln)
* post: bfpost + biases
* match_flag: flagging of solutions: 1 for best, 2 for similarly good, 0 otherwise
*

Authors
---------
Johannes Buchner (C) 2013

References
-----------
Budavari & Szalay (2008), ApJ, 679:301-309

