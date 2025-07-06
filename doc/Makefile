# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    = 
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = autoemcee
SOURCEDIR     = .
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

reset-gh-pages:
	rm -rf build/html/
	git worktree prune
	git worktree add build/html/ gh-pages
	git worktree list

update-gh-pages:
	rm -rf build/html/.buildinfo build/html/*
	$(MAKE) html
	git -C build/html add .
	git -C build/html commit -av -m 'update doc'

%.png: %.pdf
	convert -density 200 $< -background white -alpha remove -alpha off $@

runprep:
	nway-write-header.py COSMOS_OPTICAL.fits OPT  2 | tee logs/prep-OPT
	nway-write-header.py COSMOS_IRAC.fits    IRAC 2 | tee logs/prep-IRAC
	nway-write-header.py COSMOS_XMM.fits     XMM  2 | tee logs/prep-XMM

COSMOS_XMM-fake.fits: COSMOS_XMM.fits
	nway-create-fake-catalogue.py --radius 20 COSMOS_XMM.fits COSMOS_XMM-fake.fits | tee logs/XMM-fake
COSMOS_XMM-shift.fits: COSMOS_XMM.fits
	nway-create-shifted-catalogue.py --radius 40 --shift-ra 60 COSMOS_XMM.fits COSMOS_XMM-shift.fits | tee logs/XMM-shift

example2.fits: COSMOS_XMM.fits COSMOS_OPTICAL.fits
	nway.py COSMOS_XMM.fits :pos_err COSMOS_OPTICAL.fits 0.1 --out=example2.fits --radius 20 | tee logs/match2
example2-offset.fits: COSMOS_XMM-shift.fits COSMOS_OPTICAL.fits
	nway.py COSMOS_XMM-shift.fits :pos_err COSMOS_OPTICAL.fits 0.1 --out=example2-offset.fits --radius 20 | tee logs/match2-offset
example2-mag.fits: COSMOS_XMM.fits COSMOS_OPTICAL.fits
	nway.py COSMOS_XMM.fits :pos_err COSMOS_OPTICAL.fits 0.1 --out=example2-mag.fits --radius 20 --mag OPT:MAG auto --mag-radius=4 | tee logs/match2-mag-auto
	test -e example2-mag.fits
	test -e OPT_MAG_fit.txt
example2-mag-offset.fits: COSMOS_XMM-shift.fits COSMOS_OPTICAL.fits
	nway.py COSMOS_XMM-shift.fits :pos_err COSMOS_OPTICAL.fits 0.1 --out=example2-mag-offset.fits --radius 20 --mag OPT:MAG OPT_MAG_fit.txt | tee logs/match2-mag-file

example3.fits: COSMOS_XMM.fits COSMOS_OPTICAL.fits COSMOS_IRAC.fits
	nway.py COSMOS_XMM.fits :pos_err COSMOS_OPTICAL.fits 0.1 COSMOS_IRAC.fits 0.5 --out=example3.fits --radius 20 | tee logs/match3
example3-offset.fits: COSMOS_XMM-shift.fits COSMOS_OPTICAL.fits COSMOS_IRAC.fits
	nway.py COSMOS_XMM-shift.fits :pos_err COSMOS_OPTICAL.fits 0.1 COSMOS_IRAC.fits 0.5 --out=example3-offset.fits --radius 20 | tee logs/match3-offset
example3-mag.fits: COSMOS_XMM.fits COSMOS_OPTICAL.fits COSMOS_IRAC.fits
	nway.py --radius 20 COSMOS_XMM.fits :pos_err COSMOS_OPTICAL.fits 0.1 COSMOS_IRAC.fits 0.5 --mag OPT:MAG auto --mag IRAC:mag_ch1 auto --mag-radius 4 --out=example3-mag.fits | tee logs/match3-mag-auto
	test -e example3-mag.fits
	test -e IRAC_mag_ch1_fit.pdf
	test -e IRAC_mag_ch1_fit.txt
example3-mag-offset.fits: COSMOS_XMM.fits COSMOS_OPTICAL.fits COSMOS_IRAC.fits
	nway.py COSMOS_XMM-shift.fits :pos_err COSMOS_OPTICAL.fits 0.1 COSMOS_IRAC.fits 0.5 --out=example3-mag-offset.fits --radius 20 --mag OPT:MAG OPT_MAG_fit.txt  --mag IRAC:mag_ch1 IRAC_mag_ch1_fit.txt | tee logs/match3-mag-file
logs/explain: example3.fits
	nway-explain.py example3.fits 422 | tee logs/explain
logs/cutoff2: example2.fits
	nway-calibrate-cutoff.py example2.fits example2-offset.fits | tee logs/cutoff2
	${CONVERT} example2.fits_p_any_cutoffquality.pdf 
	${CONVERT} example2.fits_p_any_p_i.pdf

logs/cutoff3: example3-mag.fits
	nway-calibrate-cutoff.py example3-mag.fits example3-mag-offset.fits | tee logs/cutoff3

random_circtest.fits: randomcatX.fits randomcatO.fits
	nway.py --radius=10.0 randomcatX.fits :pos_err randomcatO.fits 0.1 --out=random_circtest.fits --min-prob=0.01
	nway-explain.py random_circtest.fits 95
random_asymtest.fits: randomcatX.fits randomcatO.fits
	nway.py --radius=10.0 randomcatX.fits :a:b randomcatO.fits 0.1 --out=random_asymtest.fits --min-prob=0.01
	nway-explain.py random_asymtest.fits 95
random_elltest.fits: randomcatX.fits randomcatO.fits
	nway.py --radius=10.0 randomcatX.fits :a:b:phi randomcatO.fits 0.1 --out=random_elltest.fits --min-prob=0.01
	nway-explain.py random_elltest.fits 95
random3_elltest.fits: randomcatX.fits randomcatO.fits randomcatR.fits
	nway.py --radius=10.0 randomcatX.fits :a:b:phi randomcatR.fits :a:b:phi randomcatO.fits 0.1 --out=random3_elltest.fits --min-prob=0.01
	nway-explain.py random3_elltest.fits 95

runexamples: COSMOS_XMM-fake.fits COSMOS_XMM-shift.fits example2.fits example2-offset.fits example2-mag.fits example2-mag-offset.fits example3.fits example3-offset.fits example3-mag.fits example3-mag-offset.fits logs/explain logs/cutoff2 logs/cutoff3 random_circtest.fits random_asymtest.fits random_elltest.fits random3_elltest.fits
runexamples: IRAC_mag_ch1_fit.png OPT_MAG_fit.png example2.fits_p_any_cutoffquality.png example2.fits_p_any_p_i.png


cleanexamples:
	rm COSMOS_XMM-fake.fits COSMOS_XMM-shift.fits example2.fits example2-offset.fits example2-mag.fits example2-mag-offset.fits example3.fits example3-offset.fits example3-mag.fits example3-mag-offset.fits logs/explain logs/cutoff2 logs/cutoff3 random_circtest.fits random_asymtest.fits random_elltest.fits random3_elltest.fits

.PHONY: runexamples runprep update-gh-pages reset-gh-pages cleanexamples

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
html: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

