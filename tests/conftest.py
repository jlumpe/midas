import pytest
from py.path import local

import numpy as np


# Path to fixtures directory
FIXTURES_DIR = local(__file__).dirpath().join('fixtures')


@pytest.fixture(scope='module')
def fixture_file():
	"""Function to get a py.path.local object for a fixture file."""

	def fixture_file_func(path):
		return FIXTURES_DIR.join(path)

	return fixture_file_func


@pytest.fixture(autouse=True)
def raise_numpy_errors():
	"""Raise exceptions for all Numpy errors in all tests."""

	old_settings = np.seterr(all='raise')

	yield

	np.seterr(**old_settings)  # Not really necessary
