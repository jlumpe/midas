import os

import pytest
from py.path import local


# Path to fixtures directory
FIXTURES_DIR = local(__file__).dirpath().join('fixtures')


@pytest.fixture(scope='module')
def fixture_file():
	"""Function to get a py.path.local object for a fixture file."""

	def fixture_file_func(path):
		return FIXTURES_DIR.join(path)

	return fixture_file_func
