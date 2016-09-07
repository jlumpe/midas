"""setuptools installation script for midas package"""

from setuptools import setup, find_packages
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy


np_include = numpy.get_include()


extensions = [
	Extension('midas.cython.*', ['midas/cython/*.pyx'],
	          include_dirs=[np_include]),
]


setup(
	name='midas',
	version='1.0.0',
	description=('Microbial IDentification And Surveillance through whole '
	             'genome sequencing data'),
	author='Jared Lumpe',
	author_email='mjlumpe@gmail.com',
	packages=find_packages(),
	namespace_packages=['midas'],
	install_requires=[
		'numpy~=1.11',
		'sqlalchemy~=1.0',
		'cython~=0.24',
	],
	ext_modules=cythonize(extensions),
)
