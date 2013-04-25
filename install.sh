#!/bin/bash
if hash pip 2>/dev/null; then
	echo "installing using pip"
	pip install progressbar pyfits argparse &&
	pip install -e . --user &&
	echo "all fine. 3way.py is in $HOME/.local/bin now"
else
	echo "installing using easy_install"
	export PYTHONPATH=$PYTHONPATH:$HOME/.local/lib/python2.7/site-packages/
	easy_install --user progressbar pyfits argparse &&
	python setup.py install --user &&
	echo "all fine. 3way.py is in $HOME/.local/bin now"
fi
