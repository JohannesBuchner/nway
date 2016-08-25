from distutils.core import setup

setup(
	name='nway',
	version='0.91',
	author='Johannes Buchner',
	author_email='johannes.buchner.acad@gmx.com',
	py_modules=['nwaylib'],
	scripts=['nway.py', 'nway-write-header.py', 'nway-explain.py'],
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


