"""setuptools installation script for midas package"""

from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy


# Cython extensions
np_include = numpy.get_include()
extensions = [Extension(
	'midas._cython.*',
	['midas/_cython/*.pyx'],
	include_dirs=[np_include],
	extra_compile_args=['-fopenmp'],
	extra_link_args=['-fopenmp'],
)]


setup(ext_modules=cythonize(extensions))
