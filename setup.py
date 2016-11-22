try:
	from setuptools import setup
except ImportError:
	from distutils.core import setup

long_description = ""
with open('README.rst') as f:
	long_description = f.read()

setup(
	name='nway',
	version='1.3',
	author='Johannes Buchner',
	author_email='johannes.buchner.acad@gmx.com',
	packages=['nwaylib'],
	scripts=['nway.py', 'nway-write-header.py', 'nway-explain.py'],
	url='http://pypi.python.org/pypi/nway/',
	license='AGPLv3 (see LICENSE file)',
	description='Probabilistic Cross-Identification of Astronomical Sources',
	long_description=long_description,
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


