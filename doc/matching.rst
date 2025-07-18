User guide: Matching catalogs
=============================

To achieve reliable results, we recommend that you run matching with
increasing amount of information (first `distance-based matching <chap:distance-based-matching>`)
then `adding priors on source properties <chap:additional-priors>`_
and understand how each influences
the results.

Typically, the goal is a compilation of “best matches”, i.e. choosing
one reliable counterpart for each source. Whatever method used, there
are always false selection fractions and false non-selection fractions
in play, which should be characterized. To this end, we recommend to
shift the source catalogues by a distance much larger than the
positional errors to simulate the results for chance alignment. This
fake catalogue should not coincide with original positions (tool
``nway-create-fake-catalogue.py`` may help). For the matching run with
this fake catalogue, use Nway  with the same settings and from the
output, choose cut-off limits (``p_any``) that correspond to the desired
false selection fraction (tool ``nway-calibrate-cutoff.py`` may help).
The output of your first Nway  run advises how to use these tools to
characterise false selection rates and to find an appropriate ``p_any``
threshold.

Then the Nway  output on the real data can be truncated based on this
criterion (``p_any>cutoff``), and made into a best match catalogue
(``match_flag==1``). It may be worth noting ambiguous cases with
multiple solutions and to store these secondary solutions with similar
probabilities in another catalogue (``match_flag==2``).

.. _`chap:distance-based-matching`:

Simple distance-based matching
------------------------------

Before exploring the full power of Nway, we consider a simple,
illustrative case. We have three catalogues, provided as FITS files:

.. container:: float

   *Input:*

   ===================== =============== ===================
   1-1 Primary Catalogue  2nd Catalogue   3rd Catalogue
   ===================== =============== ===================
   1-1 A                  :math:`\alpha`  :math:`\mathbb{A}`
   B                      :math:`\beta`   :math:`\mathbb{B}`
   ...                    ...             ...
   ===================== =============== ===================

   *Output:*

   +--------------+----------------+--------------------+-------------+----------+
   | Primary      | 2nd            | 3rd                | Probability | Group    |
   | Catalogue    | Catalogue      | Catalogue          |             |          |
   | Entry        | Entry          | Entry              |             |          |
   +==============+================+====================+=============+==========+
   | A            | :math:`\alpha` | :math:`\mathbb{A}` | ...         | A group  |
   +--------------+----------------+--------------------+-------------+----------+
   | A            | :math:`\alpha` | :math:`\mathbb{B}` | ...         |          |
   +--------------+----------------+--------------------+-------------+----------+
   | A            | :math:`\alpha` | (none)             | ...         |          |
   +--------------+----------------+--------------------+-------------+----------+
   | A            | :math:`\beta`  | :math:`\mathbb{A}` | ...         |          |
   +--------------+----------------+--------------------+-------------+----------+
   | A            | :math:`\beta`  | :math:`\mathbb{B}` | ...         |          |
   +--------------+----------------+--------------------+-------------+----------+
   | A            | :math:`\beta`  | (none)             | ...         |          |
   +--------------+----------------+--------------------+-------------+----------+
   | A            | (none)         | :math:`\mathbb{A}` | ...         |          |
   +--------------+----------------+--------------------+-------------+----------+
   | A            | (none)         | :math:`\mathbb{B}` | ...         |          |
   +--------------+----------------+--------------------+-------------+----------+
   | A            | (none)         | (none)             | ...         |          |
   +--------------+----------------+--------------------+-------------+----------+
   | B            | ...            | ...                | ...         | B group  |
   +--------------+----------------+--------------------+-------------+----------+


In Nway, only the first catalogue (**the primary catalogue**) plays a
special role. For each entry of it, counterparts are sought from the
other catalogues.

.. container::

   .. rubric:: Example - Preparing input files
      :name: example---preparing-input-files

   Note these points about preparing a catalogue input file:

   #. | Each catalogue needs to be a FITS file. The second extension
        should be the table (first extension is a header). TOPCAT writes
        files in this way.

      .. container::

         Three example catalogues are provided for you in the ``doc/``
         directory: **COSMOS_IRAC.fits, COSMOS_OPTICAL.fits and
         COSMOS_XMM.fits**. These are the same files as in Appendix B of
         :cite:t:`2018MNRAS.473.4937S`, extracted from
         :cite:t:`Sanders2007`/:cite:t:`McCracken2007`,
         :cite:t:`Ilbert2010` and
         :cite:t:`Brusa2010` respectively.

   #. | The data table needs to have a extension name and the keyword
        SKYAREA. The extension name is used as a prefix as all columns
        are copied to the output catalogue. The SKYAREA keyword tells
        the area on the sky in **square degrees** covered by the
        catalogue. This is important for estimating the chance of random
        alignments. You can use the tool
        ``python nway-write-header.py mycat.fits mytablename myskyarea``
        to set the fits header.

      .. container::

         | For our example files we have a optical, IRAC and XMM
           catalogue covering 2 square degrees:
         | ``python nway-write-header.py COSMOS_OPTICAL.fits OPT 2``
         | ``python nway-write-header.py COSMOS_IRAC.fits IRAC 2``
         | ``python nway-write-header.py COSMOS_XMM.fits XMM 2``

   #. Each catalogue needs to have a column RA and DEC providing the
      coordinates in **degrees**. To make your life easier, Nway  tries
      to be a bit fuzzy and detect the columns named RA_something etc.
      It will print out which columns it found and used.

   #. The primary catalogue needs to have a ID column. In our example
      this is the X-ray catalogue. To make your life easier, Nway\ tries
      to be a bit fuzzy and detect the columns named ID_something etc.
      It will print out which columns it found and used.

   #. Otherwise the file can have arbitrary columns which are copied
      over to the output file.

Every possible combination of association is considered. However, in
practice you do not want an extremely large output catalogue with
extremely distant, unlikely to be physically associated. You can set the
largest distance in degrees to consider by setting ``--radius``. This
speeds up the computation. But use a value that is much larger than the
largest positional error.

.. container::

   .. rubric:: Example - Matching two catalogues
      :name: example---matching-two-catalogues

   Lets try the simplest example and match the XMM X-ray catalogue to an
   optical catalogue. The XMM catalogue has a pos_err column with the
   positional error in arcseconds. For the optical catalogue we will
   assume a fixed error of 0.1 arcseconds.

   .. container::

      Run this command in the doc/ folder:

      ``python ../nway.py COSMOS_XMM.fits :pos_err COSMOS_OPTICAL.fits 0.1 --out=example1.fits --radius 15 --prior-completeness 0.9``

   Lets understand what we put in:

   #. We passed two catalogue files: COSMOS_XMM.fits and
      COSMOS_OPTICAL.fits. For the first one, we told Nway\ to use the
      column (“:”) ``pos_err`` in that catalogue for the positional
      error (**always in arcsec**). For the second one we specified a
      fixed error of 0.1 arcsec.

   #. We specified where the output should be written (``--out``).

   #. The largest XMM error is 7.3 arcsec, so we adopt a cropping radius
      of 15 arcsec to speed up the matching (``--radius 15``). A larger
      radius produces a more complete catalogue. For dense catalogues
      larger radii can be much slower to compute, as the number of
      combinations to consider rises exponentially.

   #. The parameter ``--prior-completeness 0.9`` is mentioned below.

.. container:: float

   .. container::

      Lets understand what Nway  did:

   #. ::

         NWAY arguments:
             catalogues:  COSMOS_XMM.fits, COSMOS_OPTICAL.fits
             position errors/columns:  :pos_err, 0.1
               from catalogue "XMM" (1797), density is 3.706579e+07
               from catalogue "OPT" (560536), density is 1.156188e+10
             magnitude columns:  

      It reads the catalogues and looks at their densities.

   #. ::

         matching with 15.000000 arcsec radius
         matching: 1007283192 naive possibilities
         matching: hashing
             using RA  columns: RA, RA
             using DEC columns: DEC, DEC
         matching: healpix hashing on pixel resolution ~ 18.036304 arcsec (nside=8192)
         100% | 562333|############################################|Time: 0:00:13
         matching: collecting from 61787 buckets, creating cartesian products ...
         100%|61787|###############################################|Time: 0:00:02
         matching: 462267 unique matches from cartesian product. sorting ...
         merging in 10 columns from input catalogues ...
         100% 10|##################################################|Time: 0:00:00
             adding angular separation columns
         matching:  22435 matches after filtering by search radius

      Within 20 seconds it created a cross-match of remotely possible
      associations (1,007,283,192 in principle, 22,435 within 15
      arcseconds).

   #. It found ID, RA, DEC, and positional error columns.

   #. ::

         Computing distance-based probabilities ...
           finding position error columns ...
             Position error for "XMM": found column XMM_pos_err: Values are [0.109000..7.301000]
             Position error for "OPT": using fixed value 0.100000
           finding position columns ...
           building primary_id index ...
           computing probabilities ...
               correcting for unrelated associations ... not necessary

         Computing final probabilities ...
             grouping by column "XMM_ID" and flagging ...
         100%|  1797|###################################|Time: 0:00:00

      It computed the probability of each association.

   #. ::

         creating output FITS file ...
             writing "example1.fits" (37836 rows, 17 columns)

      It wrote the output file ``example1.fits``. This file contains all
      columns from the input catalogues and the computed probabilities
      (see below for their meaning).

So how does Nway  deal with a particular, possible association and
compute its probability?

The probability of a given association is computed by comparing the
probability of a random chance alignment of unrelated sources (prior) to
the likelihood that the source is the same. The gory mathematical
details are laid out in 'Mathematical details <math>`_, but from a user
point of view the following is important:

#. The chance of a random alignment depends on the source sky density of
   the various catalogues. **So each catalogue needs to have a FITS
   header entry** **``SKYAREA``** **which tells the area covered by the
   catalogue in square degrees.** The source density on the sky is then
   computed by the number of entries divided by that area. You can use
   the tool
   ``python nway-write-header.py mycat.fits mytablename myskyarea`` to
   set the fits header.

#. Varying depths between the catalogues and different coverage can
   further reduce the fraction of expected matches. This can be adjusted
   by setting ``--prior-completeness=0.9``, if previous experience is
   that only 90% of sources have a match with the given inputs.

The outputs catalogue then contains six important new columns along with
all columns of the input catalogues:

#. ``dist_bayesfactor``: logarithm of ratio between prior and posterior
   from distance matching

#. ``dist_post``: Distance posterior probability comparing this
   association vs. no association, as in
   :cite:t:`Budavari2008`.

#. ``p_single``: Same as ``dist_post`` unless additional information was
   added, see the `section on additional priors <chap:additional-priors>`__.

#. **``p_any``**: For each entry in the primary catalogue (e.g. A) the
   probability that one of the association is the correct one is
   computed. Because every catalogue is limited by its depth, it is
   possible that the true counterpart has not been found yet. Our
   testing suggest that the **threshold for a secure catalogue depends
   on the application.**
   The `Best practice section <sec:Best-practice-matching>`__
   explains how to calibrate a threshold.

#. | **``p_i``**: For each possible association for each entry in the
     primary catalogue (e.g. A), the relative probability is computed.
     Our testing suggest that secure, pure catalogue **should keep only
     associations where ``p_i>=0.1``. Secondary solutions down to 0.1
     may be interesting. These thresholds may depend on the application
     – please report what your testing gives.**

   .. container::

      Low ``p_any`` and ``p_i`` values by themselves do not necessarily
      mean that the counterpart is ruled out. It can also mean that
      there is not enough evidence/information to declare it a
      counterpart.

#. **``match_flag``**: The most probable match is indicated with ``1``
   for each primary catalogue entry. Secondary, almost as good solutions
   are marked with ``2``. By default, the maximum allowed ratio is at
   most 0.5, but the user can modify this threshold via the
   ``--acceptable-prob`` parameter. All other associations are marked
   with ``0``.

Use the last three columns to identify sources with one solution,
possible secondary solutions, and to build final catalogues. 
'Mathematical details <math>`_ explains how these quantities are computed. To
filter out low-probability associations (low ``p_i``) from the output
catalogue, the ``--min-prob`` parameter can be used.

.. container::

   .. rubric:: Example - Output of matching two catalogues
      :name: example---output-of-matching-two-catalogues

   Lets understand the output fits file and the associations found for a
   particular X-ray source.

   .. container::

      Open the fits file and find XMM_ID=60388. As you can see from the
      ``p_i`` column, this is a ambiguous case, where more than one
      optical counterpart is possible.

   Below is an illustration of this ambiguous case (produced with
   ``python ../nway-explain.py example1.fits 60388``).

   Two sources are at a similar distance from the X-ray source (blue,
   with error circle). Therefore their association probability (``p_i``)
   is similar. The slightly higher one is marked as match_flag=1
   (orange), the other with 2 (yellow).

   `The next section <chap:additional-priors>`_ solves this by adding more
   information (the magnitude distribution). But we can also solve this
   another way. We know AGN (the X-ray source) emit in the infrared, so
   you can also match with an IRAC catalogue.

   .. container::

      Make a three-way match like so:

      ``python ../nway.py COSMOS_XMM.fits :pos_err COSMOS_OPTICAL.fits 0.1 COSMOS_IRAC.fits 0.5 --out=example3.fits --radius 15``

   However, overall we should note that ``p_any`` is low, indicating
   that probably neither of the two candidates is the counterpart.

.. figure:: example1.fits_explain_60388.png
   :alt: Visualisation of match geometry

.. _`chap:additional-priors`:

Matching with additional information
------------------------------------

For many classes of sources, the Spectral Energy Distribution (SED)
provides additional hints, which associations are likely real. For
instance, bright X-ray sources have a different color distribution in
the WISE bands than non-X-ray emitting objects. A powerful feature of
Nway  is to take advantage of this additional information to improve the
matching. The `magnitude prior section <sec:mag-priors>`__ has the mathematical details
and a comparison to the Likelihood Ratio method.

.. container::

   .. rubric:: Example - Using magnitude information
      :name: example---using-magnitude-information

   X-ray sources (which we are looking for in our example) have a
   different optical magnitude distribution than non-X-ray emitting
   objects. Lets take advantage of this information:

   .. container::

      Run this command:

      ``python ../nway.py COSMOS_XMM.fits :pos_err COSMOS_OPTICAL.fits 0.1 --out=example2.fits --radius 15 --prior-completeness 0.9 --mag OPT:MAG auto --mag-radius 3.5``

   The last two parts are new:

   ``--mag OPT:MAG auto --mag-radius 3.5``

   We use the column ``MAG`` from the catalogue ``OPT`` (FITS table
   name), therefore ``--mag OPT:MAG``. After this follows that the
   magnitude prior histogram should be generated from the data (mode
   ``auto``), by comparing the ``MAG`` histogram of sources within 3.5
   arcsec of a X-ray source (``--mag-radius``) to that of full
   histogram.

   *(example continued below)*

There are three possible ways to specify the prior in Nway: In all cases
you specify ``--mag column-name [filename|auto]``. You can use ``--mag``
several times.

#. “File-mode”: If we know the magnitude distribution of X-ray detected
   AGN we can provide this prior distribution as a table (histogram).
   This table contains the color histogram of the sources of interest
   (X-ray detected AGN) and a histogram of other, field sources (more
   details below ).

#. “Simple auto-mode”: Specifying ``auto`` instead of a file name
   derives the two distributions from the data, as we did in our
   example: All sources inside 3.5 arcseconds (``--mag-radius``
   parameter) of a X-ray source are put into one histogram, and all
   others into another histogram.

#. “Bayesian auto-mode”: Bayesian distance probabilities (``dist_post``)
   will be used if you leave out ``--mag-radius``. This is in general
   safer and recommended. In small catalogues the histogram may not be
   sufficiently filled, in which case Nway  will give a warning (more
   details below ).

.. container::

   .. container::

      Lets look at the histograms it computed. Nway  created
      ``OPT_MAG_fit.pdf``, and also ``OPT_MAG_fit.txt`` as a histogram
      file:

      .. image:: OPT_MAG_fit.png
         :alt: histogram of optical magnitudes for target and field population

      They are clearly different: Lower magnitude (bright) sources are
      more likely associated to X-ray sources. This will help our
      matching.

   As an example, we show below the ambiguous case from before. The
   upper association has been selected because it has a better match by
   magnitude, resolving the ambiguity.

   .. image:: example2.fits_explain_60388.png
      :alt: visualisation of match geometry

Multiple priors
'''''''''''''''

You can specify as many priors as you like, taking advantage of more and
more information. Just repeat --mag.

.. container::

   The following example uses one prior from the optical catalogue and
   another prior from the IRAC catalogue. A three-way match is
   performed.
   ``python ../nway.py COSMOS_XMM.fits :pos_err COSMOS_OPTICAL.fits 0.1 COSMOS_IRAC.fits 0.5 --out=example3.fits --radius 20 --prior-completeness 0.9 --mag OPT:MAG auto --mag IRAC:mag_ch1 auto --mag-radius 3.5``

Providing a prior as a file
'''''''''''''''''''''''''''

In the paper we demonstrate the use of a WISE magnitude of X-ray
sources. If such prior information comes from previous studies, the
distributions can be passed to Nway  as a ASCII table histogram. This
table contains the histogram of the sources of interest (X-ray sources)
and a histogram of other sources (non-X-ray sources). The file
``OPT_MAG_fit.txt`` is an example of such a input file, and can be used
via ``--mag OPT:MAG OPT_MAG_fit.txt``. It contains four columns (lower
and upper bin edge, density of selected and non-selected) and looks like
this (# indicates comments):

::

   # OPT_MAG_fit.txt
   # lo         hi          selected   others
     10.76000   18.98571    0.00870    0.00183
     18.98571   20.27286    0.05562    0.01448
     ...

.. container::

   Keep in mind that a prior created from a different data set can only
   be used if it is applicable to the present data set. For example, in
   the introduction of the paper
   :cite:t:`2018MNRAS.473.4937S` we stress that a prior from
   a comparable X-ray exposure depth must be used when deriving color
   distributions.

A general approach
''''''''''''''''''

Providing priors is not limited to magnitude distributions, you can use
colors or any other information you want (e.g. morphology, variability,
etc.). The approach is very general, Nway  just looks at the
corresponding bin and reweighs the probabilities. For example, in
:cite:t:`2018MNRAS.473.4937S`, the counterparts to ROSAT
sources where found using WISE. The prior was build by using the
color-magnitude (W1-W2 vs W2) properties of ~3000 secure counterparts to
the 3XMM-Bright survey cut at the depth reached by ROSAT.

Discovering a prior from distance matching
''''''''''''''''''''''''''''''''''''''''''

If you set ``--mag OPT:MAG auto`` and do not set ``--mag-radius``, Nway 
uses the Bayesian distance matching for discovering the histogram of
``OPT:MAG``, as follows:

#. Those with ``dist_post>0.9`` are considered safe matches and are used
   for the “selected” histogram.

#. Those with ``dist_post<0.01`` are considered safe non-matches and are
   used for the “others” histogram.

#. Entries of -99 are always ignored. It is usually better to assign -99
   where the magnitude error is large, to get cleaner histograms.

This is in general more cautious, and recommended for large catalogues

.. container::

   However, if you only have a small catalogue you may build a poorly
   sampled histogram, potentially leading to biases. Nway  will warn you
   when only few sources were selected.

.. bibliography::
