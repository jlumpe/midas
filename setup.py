"""setuptools installation script for midas package"""

from setuptools import setup, find_packages


setup(
	name='midas',
	version='0.0.1',
	description=('Microbial IDentification And Surveillance through whole '
	             'genome sequencing data'),
	author='Jared Lumpe',
	author_email='mjlumpe@gmail.com',
	packages=find_packages(),
	install_requires=[
		'python >= 3.5',
		'numpy >= 1.11',
	],
)
