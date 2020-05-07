********************
Package Installation
********************


Prerequisites
=============

This package requires both `Cython <https://cython.org/>`_ and
`Numpy <https://numpy.org/>`_ be preinstalled.

It is highly recommended to use the `conda <https://docs.conda.io/en/latest/>`_
package manager to install this and related packages.


Instructions
============

This package is installed with ``setuptools``.
Clone the repo from GitHub and ``cd`` into the directory. Then run::

	python setup.py build_ext --inplace
	python setup.py install

To do a development install (don't need to reinstall the package after making
changes to the source code), substitute ``develop`` for ``install``::

	python setup.py build_ext --inplace
	python setup.py develop

