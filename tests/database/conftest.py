import pytest

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from midas.db.sqla import ReadOnlySession


@pytest.fixture(scope='session')
def testdb_dir(test_data):
	"""Directory containing testdb_210126 data."""
	return test_data / 'testdb_210126'


@pytest.fixture(scope='session')
def testdb_session(testdb_dir):
	"""Function which creates a new session for the test database."""
	engine = create_engine('sqlite:///' + str(testdb_dir / 'testdb_210126.db'))
	return sessionmaker(engine)
