from distutils.core import setup

setup(
	name='3way',
	version='0.7',
	author='Johannes Buchner',
	author_email='jbuchner@mpe.mpg.de',
	py_modules=['fastskymatch', 'bayesdistance', 'magnitudeweights'],
	scripts=['3way.py'],
	url='http://pypi.python.org/pypi/3way/',
	license='LICENSE',
	description='Probabilistic Cross-Identification of Astronomical Sources',
	long_description=open('README.rst').read(),
	install_requires=[
		"scipy",
		"numpy",
		"pyfits",
		"progressbar",
		"matplotlib",
		"argparse",
		"joblib"
	],
)


