"""Test the midas.db.migrate module."""

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from alembic.migration import MigrationContext
from alembic.script import ScriptDirectory

from midas.db.migrate import (current_head, current_revision, is_current_revision, init_db,
                              get_alembic_config)
from midas.db import models


# Expected current head revision
# Need to update this value every time a new revision is introduced
CURRENT_HEAD = 'd961d0698083'

# Old revision number to test. Must actually exist in the scripts directory.
OLD_REVISION = '955d795123c3'


def test_current_head():
	assert current_head() == CURRENT_HEAD


class TestCurrentRevision:
	"""Test the current_revision() and is_current_revision() functions."""

	def test_uninitialized(self):
		"""Test on uninitialized database (not stamped)."""
		engine = create_engine('sqlite:///:memory:')
		assert current_revision(engine) is None
		assert not is_current_revision(engine)

	def test_empty(self):
		"""Test on freshly initialized database."""
		engine = create_engine('sqlite:///:memory:')
		init_db(engine)
		assert current_revision(engine) == CURRENT_HEAD
		assert is_current_revision(engine)

	def test_old(self):
		"""Test on uninitialized database stamped with an old revision no."""
		engine = create_engine('sqlite:///:memory:')
		conf = get_alembic_config(engine)
		scriptdir = ScriptDirectory.from_config(conf)

		with engine.connect() as conn:
			ctx = MigrationContext.configure(conn)
			ctx.stamp(scriptdir, OLD_REVISION)

		assert current_revision(engine) == OLD_REVISION


def test_init_db():
	"""Test the init_db() function."""
	engine = create_engine('sqlite:///:memory:')
	init_db(engine)

	# Check current revision matches
	assert current_revision(engine) == current_head()

	# Check we can query all models (won't return any results, but would fail if tables didn't exist).
	session = sessionmaker(engine)()
	for model in [models.Genome, models.ReferenceGenomeSet, models.AnnotatedGenome, models.Taxon]:
		session.query(model).all()
