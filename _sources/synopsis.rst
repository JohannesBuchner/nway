Command line arguments
======================

Program Arguments
-----------------

#. ``python ../nway.py --help`` ... display help page

.. literalinclude:: logs/help

#. ``python ../nway.py catalogue1 catalogues ...`` ... input catalogues
   as FITS files.

#. ``python ../nway.py --out`` ... Output file name (also a FITS file).

#. ``python ../nway.py --radius`` This radius (in degrees) is used to
   discard distant pairs. Always choose a value that is much larger than
   the largest positional uncertainty, then this value will not change
   the results. Smaller values make the code run faster and use less
   memory by reducing the number of combinations to explore.

#. ``python ../nway.py --prior-completeness 1`` ... set expected
   matching completeness (default: 1)

#. ``python ../nway.py --mag MAGCOLUMN MAGFILE`` ... name of
   <table>:<column> for magnitude biasing, and file name for magnitude
   histogram (use auto for auto-computation within mag-radius). Example:
   ``--mag OPT:MAG auto --mag IRAC:mag_irac1 irac_histogram.txt``

#. ``python ../nway.py --mag-radius`` ... If set, and a auto prior is
   defined, then the selected sources are taken from within this radius
   of the primary sources (in arc seconds). If not set (recommended),
   the Bayesian posterior from distance matching is used, which
   incorporates positional errors.

#. ``python ../nway.py --min-prob`` ... only retain associations in the
   output catalogue exceeding this ``p_i`` value. Recommended: 0.1.

#. ``python ../nway.py --acceptable-prob`` ... affects the flagging of
   secondary solutions (``match_flag`` column). If the secondary is
   within this difference (default: 0.005), it is marked as a secondary
   solution.

Input file specifications and units
-----------------------------------

#. Each catalogue needs to be a FITS file. The second extension should
   be the table (first extension is a header).

#. The data table needs to have a extension name.

#. The header of the data table needs the keyword SKYAREA, which
   specifies the area covered by the catalogue in **square degrees**.

#. Each catalogue needs to have a column RA and DEC **in degrees**. To
   make your life easier, Nway  tries to be a bit fuzzy and detect the
   columns named RA_something etc.

#. The primary catalogue needs to have a ID column. To make your life
   easier, Nway  tries to be a bit fuzzy and detect the columns named
   ID_something etc.

#. Positional error columns, if used, need to be **in arcseconds**.

Example catalogues are provided in the ``doc/`` directory:
COSMOS_IRAC.fits, COSMOS_OPTICAL.fits and COSMOS_XMM.fits.
