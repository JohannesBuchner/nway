from distutils.core import setup

setup(
	name='nway',
	version='0.9',
	author='Johannes Buchner',
	author_email='johannes.buchner.acad@gmx.com',
	py_modules=['fastskymatch', 'bayesdistance', 'magnitudeweights'],
	scripts=['nway.py'],
	url='http://pypi.python.org/pypi/nway/',
	license='LICENSE',
	description='Probabilistic Cross-Identification of Astronomical Sources',
	long_description=open('README.rst').read(),
	install_requires=[
		"scipy",
		"numpy",
		"astropy",
		"progressbar",
		"matplotlib",
		"argparse",
		"joblib"
	],
)


