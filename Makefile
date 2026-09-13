.PHONY: test

test:
	echo "backend: Agg" > matplotlibrc
	cd doc && echo "backend: Agg" > matplotlibrc
	nway-write-header.py doc/COSMOS_OPTICAL.fits OPT  2
	nway-write-header.py doc/COSMOS_IRAC.fits    IRAC 2
	nway-write-header.py doc/COSMOS_XMM.fits     XMM  2
	cd doc && nway-create-fake-catalogue.py --radius 20 COSMOS_XMM.fits COSMOS_XMM-fake.fits
	cd doc && nway-create-shifted-catalogue.py --radius 40 --shift-ra 60 COSMOS_XMM.fits COSMOS_XMM-shift.fits
	cd doc && rm -f COSMOS_XMM-shift.fits
	cd doc && nway-create-shifted-catalogue.py --radius 40 --shift-dec 60 COSMOS_XMM.fits COSMOS_XMM-shift.fits
	cd doc && nway.py COSMOS_XMM.fits :pos_err COSMOS_OPTICAL.fits 0.1 --out=example2.fits --radius 20
	cd doc && nway.py COSMOS_XMM-shift.fits :pos_err COSMOS_OPTICAL.fits 0.1 --out=example2-offset.fits --radius 20
	cd doc && nway.py COSMOS_XMM.fits :pos_err COSMOS_OPTICAL.fits 0.1 --out=example2-mag.fits --radius 20 --mag OPT:MAG auto --mag-radius=4
	cd doc && nway.py COSMOS_XMM-shift.fits :pos_err COSMOS_OPTICAL.fits 0.1 --out=example2-mag-offset.fits --radius 20 --mag OPT:MAG OPT_MAG_fit.txt
	cd doc && nway.py COSMOS_XMM.fits :pos_err COSMOS_OPTICAL.fits 0.1 COSMOS_IRAC.fits 0.5 --out=example3.fits --radius 20
	cd doc && nway.py COSMOS_XMM-shift.fits :pos_err COSMOS_OPTICAL.fits 0.1 COSMOS_IRAC.fits 0.5 --out=example3-offset.fits --radius 20
	cd doc && nway.py --radius 20 COSMOS_XMM.fits :pos_err COSMOS_OPTICAL.fits 0.1 COSMOS_IRAC.fits 0.5 --mag OPT:MAG auto --mag IRAC:mag_ch1 auto --mag-radius 4 --out=example3-mag.fits
	cd doc && nway.py COSMOS_XMM-shift.fits :pos_err COSMOS_OPTICAL.fits 0.1 COSMOS_IRAC.fits 0.5 --out=example3-mag-offset.fits --radius 20 --mag OPT:MAG OPT_MAG_fit.txt --mag IRAC:mag_ch1 IRAC_mag_ch1_fit.txt
	cd doc && nway-explain.py example3.fits 422
	cd doc && nway-calibrate-cutoff.py example2.fits example2-offset.fits
	cd doc && nway-calibrate-cutoff.py example3-mag.fits example3-mag-offset.fits
	cd tests/elltest && nway.py --radius=10.0 randomcatX.fits :pos_err randomcatO.fits 0.1 --out=random_circtest.fits --min-prob=0.01
	cd tests/elltest && nway-explain.py random_circtest.fits 95
	cd tests/elltest && nway.py --radius=10.0 randomcatX.fits :a:b randomcatO.fits 0.1 --out=random_asymtest.fits --min-prob=0.01
	cd tests/elltest && nway-explain.py random_asymtest.fits 95
	cd tests/elltest && nway.py --radius=10.0 randomcatX.fits :a:b:phi randomcatO.fits 0.1 --out=random_elltest.fits --min-prob=0.01
	cd tests/elltest && nway-explain.py random_elltest.fits 95
	cd tests/elltest && nway.py --radius=10.0 randomcatX.fits :a:b:phi randomcatR.fits :a:b:phi randomcatO.fits 0.1 --out=random3_elltest.fits --min-prob=0.01
	cd tests/elltest && nway-explain.py random3_elltest.fits 95

release: dist ## package and upload a release
	twine upload --verbose dist/*.tar.gz

install: clean ## install the package to the active Python's site-packages
	$(PYTHON) -m pip install .
