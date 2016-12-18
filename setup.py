"""setuptools installation script for midas package"""

from setuptools import setup, find_packages
from distutils.util import convert_path
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy


# Get package version without importing it
version_ns = dict()
with open(convert_path('midas/version.py')) as fobj:
	exec(fobj.read(), version_ns)
version = version_ns['__version__']


# Dependencies
install_requires = [
	'numpy~=1.11',
	'sqlalchemy~=1.0',
	'cython~=0.24',
	'alembic~=0.8',
]


# Cython extensions
np_include = numpy.get_include()
extensions = [
	Extension('midas.cython.*', ['midas/cython/*.pyx'],
	          include_dirs=[np_include]),
]


setup(
	name='midas',
	version=version,
	description=(
		'Microbial IDentification And Surveillance through whole genome '
		'sequencing data'
	),
	author='Jared Lumpe',
	author_email='mjlumpe@gmail.com',
	packages=find_packages(),
	install_requires=install_requires,
	ext_modules=cythonize(extensions),
	include_package_data=True,
)
