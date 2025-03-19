import os
try:
	from setuptools import setup
except ImportError:
	from distutils.core import setup

long_description = ""
with open('README.rst') as f:
	long_description = f.read()

setup(
	packages=['nwaylib'],
	scripts=['nway.py', 'nway-write-header.py', 'nway-explain.py', 'nway-create-fake-catalogue.py', 'nway-create-shifted-catalogue.py', 'nway-calibrate-cutoff.py'],
	long_description=long_description,
)
