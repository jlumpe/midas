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
	'numpy~=1.13',
	'sqlalchemy~=1.1',
	'biopython~=1.69',
	'alembic~=0.9',
	'pydatatypes~=0.1',
]


# Cython extensions
np_include = numpy.get_include()
extensions = [Extension(
	'midas.cython.*',
	['midas/cython/*.pyx'],
	include_dirs=[np_include],
	compiler_directives={
		'language_level': 3,
	},
	extra_compile_args=['-fopenmp'],
	extra_link_args=['-fopenmp'],
)]


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
	namespace_packages=['midas'],
	install_requires=install_requires,
	dependency_links=[
		'https://github.com/jlumpe/pydatatypes/archive/master.tar.gz#egg=pydatatypes-0.1',
	],
	ext_modules=cythonize(extensions),
	include_package_data=True,
)
