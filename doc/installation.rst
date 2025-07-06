.. _install:
.. highlight:: shell

============
Installation
============


Stable release
--------------

nway is a pure Python program.

To install nway, run this command in your terminal:

.. code-block:: console

    $ pip install nway

This is the preferred method to install nway, as it will always install the most recent stable release.

The command line tool ``nway.py`` is installed for you as ``~/.local/bin/nway.py``.

If that directory is in your ``$PATH``, you can run::

    $ nway.py --help

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/

Upgrading from older Nway  versions
'''''''''''''''''''''''''''''''''''

To upgrade, uninstall the older Nway  version first:

``$ pip uninstall nway``

Repeat until it says that nway is not installed. Then follow the
installation instructions above.


Download example data
---------------------

For obtaining the example data files, 
see below, installation "from sources".
The data files (*.fits) is in the doc/ folder.

From sources
------------

The sources for nway can be downloaded from the `Github repo`_.


You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/JohannesBuchner/nway

Or download the `tarball`_:

.. code-block:: console

    $ curl -OJL https://github.com/JohannesBuchner/nway/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install

You need these python packages:

``scipy, astropy, matplotlib, tqdm, argparse, joblib, healpy, pandas``.


.. _Github repo: https://github.com/JohannesBuchner/nway
.. _tarball: https://github.com/JohannesBuchner/nway/tarball/master
