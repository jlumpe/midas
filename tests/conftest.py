from pathlib import Path

import numpy as np
import pytest


@pytest.fixture(scope='session')
def test_data():
	"""The directory containing test data."""
	return Path(__file__).parent / 'data'


@pytest.fixture(autouse=True)
def raise_numpy_errors():
	"""Raise exceptions for all Numpy errors in all tests."""

	old_settings = np.seterr(all='raise')

	yield

	np.seterr(**old_settings)  # Not really necessary
