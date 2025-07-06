
.. _`chap:math`:

Mathematical details and Implementation
=======================================

.. _`sec:math`:

Distance-based matching
-----------------------

Lets consider the problem of finding counterparts to a primary catalogue
(:math:`i=1`), in our example for the X-ray source position catalogue.
Let each :math:`N_{i}` denote the number of entries for the catalogues
used, and :math:`\nu_{i}=N_{i}/\Omega_{i}` denote their source surface
density on the sky.

If a counterpart is required to exist in each of the :math:`k`
catalogues, there are :math:`\prod_{i=1}^{k}N_{i}` possible
associations. If we assume that a counterpart might be missing in each
of the matching catalogues, there are
:math:`N_{1}\cdot\prod_{i=2}^{k}(N_{i}+1)` possible associations. This
minor modification, negligible for :math:`N_{i}\gg1`, is ignored in the
following for simplicity, but handled in the code.

If each catalogue covers the same area with some respective, homogeneous
source density :math:`\nu_{i}`, the probability of a chance alignment on
the sky of physically unrelated objects can then be written
(:cite:t:`Budavari2008`, eq. 25), as

.. math:: P(H)=N_{1}/\prod_{i=1}^{k}N_{i}=1/\prod_{i=2}^{k}N_{i}=1/\prod_{i=2}^{k}\nu_{i}\Omega_{i}.

Thus :math:`P(H)` is the prior probability of an association. The
posterior should strongly exceed this prior probability, to avoid false
positives.

To account for non-uniform coverage, :math:`P(H)` is modified by a
“prior completeness factor” :math:`c`, which gives the expected fraction
of sources with reliable counterpart (due to only partial coverage of
the matching catalogues :math:`\Omega_{i>1}\neq\Omega_{1}`, depth of the
catalogues and/or systematic errors in the coordinates). Our prior can
thus be written as

.. _`eq:prior`:

.. math:: P(H)=c/\prod_{i=2}^{k}\nu_{i}\Omega_{1}.\label{eq:prior}

Bayes’ theorem connects the prior probability :math:`P(H)` to the
posterior probability :math:`P(H|D)`, by incorporating information
gained from the observation data :math:`D` via

.. _`eq:bayes`:

.. math:: P(H|D)\propto P(H)\times P(D|H).\label{eq:bayes}

We now extend the approach of :cite:t:`Budavari2008`, to
allow matches where some catalogues do not participate in a match.
Comparing A12 and A14 in :cite:t:`Budavari2008`, assuming
that positions lie on the celestial sphere and adopting the expansions
developed in their Appendix B, we can write down likelihoods. For a
counterpart across :math:`k` catalogues, we obtain:

.. _`eq:nwaylikelihood`:

.. math:: P(D|H)=2^{k-1}\frac{\prod\sigma_{i}^{-2}}{\sum\sigma_{i}^{-2}}\exp\left\{ -\frac{\sum_{i<j}\psi_{ij}^{2}\sigma_{j}^{-2}\sigma_{i}^{-2}}{2\sum\sigma_{i}^{-2}}\right\} \label{eq:nwaylikelihood}

For a given association with participating members from :math:`k`
catalogues, the pairwise angular separation, :math:`\psi_{ij}`, between
catalogue :math:`i` and :math:`j` is judged with the relevant position
uncertainties :math:`\sigma`. The likelihood for the hypothesis where
some catalogues do not participate in the association has the
appropriate terms in the products and sums removed. Therefore, the
likelihood is unity for the hypothesis that there is no counterpart in
any of the catalogues.

In comparison to our method, the method of
:cite:t:`Budavari2008` only compares two hypotheses for a
association: either all sources belong to the same object
(:math:`H_{1}`), or they are coincidentally aligned (:math:`H_{0}`). In
this computation each hypothesis test is run in isolation, and relative
match probabilities for a given source are not considered. For
completeness, we also compute the posterior of this simpler model
comparison:

.. _`eq:assocPost`:

.. math::

   \begin{aligned}
   \frac{P(H_{1}|D)}{P(H_{0}|D)} & \propto & \frac{P(H_{1})}{P(H_{0})}\times\frac{P(D|H_{1})}{P(D|H_{0})}\\
   B & = & \frac{P(D|H_{1})}{P(D|H_{0})}\\
   P(H_{1}|D) & = & \left[1+\frac{1-P(H_{1})}{B\cdot P(H_{1})}\right]^{-1}\label{eq:assocPost}
   \end{aligned}

The output column ``dist_bayesfactor`` stores :math:`\log B`, while the
output column ``dist_post`` is the result of equation
:ref:`eq:assocPost`. The output column ``p_single`` gives
``dist_post`` but modified if any additional information is specified
(see Section `5.2 <sec:mag-priors>`__). As mentioned several times in
the literature, the :cite:t:`Budavari2008` approach does not
include sources absent in some of the catalogues, while the formulae we
develop below incorporate absent sources. This is similar in spirit to
:cite:p:`Pineau2016`, although the statistical approach is
different. We now go further and develop counterpart probabilities.

The first step in catalogue inference is whether the source has any
counterpart (:math:`p_{\mathrm{any}}`). The posterior probabilities
:math:`P(H|D)` are computed using Bayes theorem (eq.
`[eq:bayes] <eq:bayes>`_) with the likelihood (eq.
`[eq:nwaylikelihood] <eq:nwaylikelihood>`_) and prior (eq.
`[eq:prior] <eq:prior>`_) appropriately adopted for the number of
catalogues the particular association draws from. For each entry in the
primary catalogue, the posteriors of all possible associations are
normalised to unity, and :math:`P(H_{0}|D)`, the posterior probability
of the no-counterpart hypothesis, i.e., no catalogue participates,
computed. From this we compute:

.. _`eq:post-any`:

.. math:: p_{\mathrm{any}}=1-P(H_{0}|D)/\sum_{i}P(H_{i}|D)\label{eq:post-any}

If :math:`p_{\mathrm{any}}` is low, this indicates that there is little
evidence for any of the considered, combinatorically possible
associations, except for the no-association case. The output column
``p_any`` is the result of equation `[eq:post-any] <eq:post-any>`__.

If :math:`p_{\mathrm{any}}\approx1`, there is strong evidence for at
least one of the associations to another catalogue. To compute the
relative posterior probabilities of the options, we re-normalize with
the no-counterpart hypothesis, :math:`H_{0}`, excluded:

.. _`eq:post-assoc`:

.. math:: p_{i}=P(H_{i}|D)/\sum_{i>0}P(H_{i}|D)\label{eq:post-assoc}

If a particular association has a high :math:`p_{i}`, there is strong
evidence that it is the true one, out of all present options. The output
column ``p_i`` is the result of equation
`[eq:post-assoc] <eq:post-assoc>`__.

A “very secure” counterpart could be defined by the requirement
:math:`p_{any}>95\%` and :math:`p_{i}>95\%`, for example. However, it is
useful to run simulations to understand the rate of false positives.
Typically, much lower thresholds are acceptable.

.. _`sec:mag-priors`:

Magnitudes, Colors and other additional information
---------------------------------------------------

Astronomical objects of various classes often show distinct color and
magnitude distributions. Because most bright X-ray point-sources in deep
images are also optically bright compared to generic sources, this
information can be exploited. Previous works
(e.g., :cite:t:`Brusa2005`, ,:cite:t:`Brusa2007`) have modified the
likelihood ratio coming from the angular distance :math:`f(r)`
information (likelihood ratio method,
:cite:t:`SutherlandSaunders1992`) by a factor:

.. math:: LR=\frac{q(m)}{n(m)}\times f(r)

Here, :math:`q(m)` and :math:`n(m)` are associated with the magnitude
distributions of source (e.g. X-ray sources) and background objects
(e.g. stars, passive galaxies) respectively, but additionally contain
sky density contributions.

This idea can be put on solid footing within the Bayesian framework.
Here, two likelihoods are combined, by simply considering two
independent observations, namely one for the positions,
:math:`D_{\phi}`, and one for the magnitudes :math:`D_{m}`. The
likelihood thus becomes

.. math::

   \begin{aligned}
   P(D|H) & = & P(D_{\phi}|H)\times P(D_{m}|H)\\
    & = & P(D_{\phi}|H)\times\frac{\bar{q}(m)}{\bar{n}(m)},
   \end{aligned}

with :math:`\bar{q}(m)` and :math:`\bar{n}(m)` being the probability
that a X-ray (target) source or a generic (field) source has magnitude
:math:`m` respectively. Nway\ stores the modifying factor,
:math:`P(D_{m}|H)`, in ``bias_``\ ``*`` output columns, one for each
column giving a magnitude, color, or other distribution. This modifying
factor is however renormalized so that
:math:`P(D_{m}|H)=\frac{\bar{q}(m)}{\bar{n}(m)}/\int\frac{\bar{q}(m')}{\bar{n}(m')}\bar{n}(m')dm'`,
which makes :math:`P(D|H)=P(D_{\phi}|H)` when :math:`m` is unknown. In
that case, :math:`m` is marginalised over its distribution in the
general population, i.e. :math:`\int P(D_{m}|H)\,\bar{n}(m')\,dm`. This
has the benefit that when m is unknown, the modifying factor is unity
and the probabilities remain unmodified.

For completeness, I mention the fully generalized case. This is attained
when an arbitrary number of photometry bands are considered, each
consisting of a magnitude measurement :math:`m` and measurement
uncertainty :math:`\sigma_{m}`:

.. math:: P(D_{m}|H)=\prod\frac{\int_{m}\bar{q}(m)\,p(m|D_{m})\,dm}{\int_{m}\bar{n}(m)\,p(m|D_{m})\,dm}

Here, :math:`p(m|D_{m})` would refer to a Gaussian error distribution
with mean :math:`m` and standard deviation :math:`\sigma_{m}`. This is
convolved with the distribution properties. Alternatively,
:math:`p(m|D_{m})` can also consider upper limits. However, such options
are not yet implemented in Nway. Instead, we recommend removing
magnitude values with large uncertainties (setting them to -99).


.. _`sec:Auto-calibration`:

Auto-calibration
----------------

The probability distributions :math:`\bar{n}(m)` and :math:`\bar{q}(m)`
can be taken from other observations by computing the magnitude
histograms of the overall population and the target sub-population (e.g.
X-ray sources).

Under certain approximations and assumptions, these histograms can also
be computed during the catalogue matching procedure while also being
used for the weighting. One could perform the distance-based matching
procedure laid out above, and compute a magnitude histogram of the
secure counterparts as an approximation for :math:`\bar{q}(m)` and a
histogram of ruled out counterparts for :math:`\bar{n}(m)`. While the
weights :math:`\bar{q}(m)/\bar{n}(m)` may strongly influence the
probabilities of the associations for a single object, the bulk of the
associations will be dominated by distance-weighting. One may thus
assume that the :math:`\bar{q}(m)` and :math:`\bar{n}(m)` are computed
with and without applying the magnitude weighting are the same, which is
true in practice. When differences are noticed, they will only
strengthen :math:`\bar{q}(m)`, and the procedure may be iterated.

.. _`sec:implementation`:

Implementation
--------------

This section gives some details connecting the math above to
a concrete and efficient implementation on a computer.

The implementation for matching :math:`n` catalogues is a Python program
called nway. The input catalogues have to be in FITS format. Information
about the (shared) sky coverage has to be provided to the program as
well. The program proceeds in four steps.

First, possible associations are found. It is unfeasible and unnecessary
to consider all theoretical possibilities (complexity
:math:`O(\prod_{i=1}^{k}N_{i})`), so the sky is split first to cluster
nearby objects. For this, a hashing procedure puts each object into
HEALPix bins :cite:p:`Gorski2005`. The bin width :math:`w` is
chosen so that any association of distance :math:`w` are improbable and
negligible in practice, i.e. much larger than the largest positional
error. An object with coordinates :math:`\phi,\,\theta` is placed in the
bin corresponding to its coordinate, but also into its neighboring bins
to avoid boundary effects. This is done for each catalogue separately.
Then, in each bin, the Cartesian product across catalogues (every
possible combination of sources) is computed. All associations are
collected across the bins and filtered to be unique. The hashing
procedure adds very low effort :math:`O(\sum_{i=1}^{k}N_{i})` while the
Cartesian product is reduced drastically to
:math:`O(N_{\text{bins}}\cdot\prod_{i=1}^{k}\frac{N_{i}}{N_{bins}})`.
All primary objects that have no associations past this step have
:math:`P(\text{"any real association"}|D)=0`.

The second step is the computation of posteriors using the angular
distances between counterparts. The prior is also evaluated from the
size of the catalogue and the effective coverage, as well as the
user-supplied prior incompleteness factor. The posterior for each
association based on the distances only is calculated. These posteriors
have to be modified (“correcting for unrelated associations”), to
consider associations unrelated to primary catalogue sources (described
in the paper, :cite:t:`2018MNRAS.473.4937S`, in the appendix
section “Computing all possible matches”).

In the third step the magnitudes are considered, and the posteriors
modified. An arbitrary number of magnitude columns in the input
catalogues can be specified. It is possible to use external magnitude
histograms (e.g. for sparse matching with few objects) as well as
computing the histograms from the data itself (see Section
`[subsec:Auto-calibration] <subsec:Auto-calibration>`__). The breaks of
the histogram bins are computed adaptively based on the empirical
cumulative distribution found. Because the histogram bins are usually
larger than the magnitude measurement uncertainty, the latter is
currently not considered. The adaptive binning creates bin edges based
on the number of objects, and is thus independent of the chosen scale
(magnitudes, flux). Thus the method is not limited to magnitudes, but
can be used for virtually any other known object property (colours,
morphology, variability, etc.).

In the final step, associations are grouped by the object from the
primary catalogue (here: X-ray source catalogue). The posteriors
:math:`p_{\mathrm{any}}` and :math:`p_{i}` are computed. For the output
catalogue a cut on the posterior probability (e.g. above 80%) can be
applied, and all associations with their posterior probability are
written to the output fits catalogue file.
