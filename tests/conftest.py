from pathlib import Path

import numpy as np
import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from midas.db.models import Base as models_base


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


@pytest.fixture(scope='session')
def make_empty_db():
	"""Function which creates an empty in-memory-database with initialized schema."""
	def empty_db_factory():
		engine = create_engine('sqlite:///:memory:')
		models_base.metadata.create_all(engine)
		return engine

	return empty_db_factory


@pytest.fixture(scope='session')
def testdb_dir(test_data):
	"""Directory containing testdb_210126 data."""
	return test_data / 'testdb_210126'


@pytest.fixture(scope='session')
def testdb_session(testdb_dir):
	"""Function which creates a new session for the test database."""
	engine = create_engine('sqlite:///' + str(testdb_dir / 'testdb_210126.db'))
	return sessionmaker(engine)
