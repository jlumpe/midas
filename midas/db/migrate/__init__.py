"""Perform genome database migrations with Alembic.

This package also contains all Alembic data files.
"""

from typing import Optional

from alembic.config import Config
from alembic import command
from alembic.migration import MigrationContext
from alembic.script import ScriptDirectory
from pkg_resources import resource_filename
from sqlalchemy.engine import Engine


INI_PATH = resource_filename(__name__, 'alembic.ini')
SCRIPTS_DIR = resource_filename(__name__, 'alembic')


def get_alembic_config(engine: Optional[Engine] = None, **kwargs) -> Config:
	"""Get an alembic config object to perform migrations.

	Parameters
	----------
	engine
		SQLAlchemy engine specifying database connection info (optional). Assigned to ``'engine'``
		key of :attr:`alembic.config.Config.attributes`.
	\\**kwargs
		Keyword arguments to pass to :meth:`alembic.config.Config.__init__`.

	Returns
	-------
		Alembic config object.
	"""

	config = Config(INI_PATH, **kwargs)
	config.set_main_option('script_location', SCRIPTS_DIR)

	if engine is not None:
		config.attributes['engine'] = engine

	return config


def current_head() -> str:
	"""Get the current head revision number."""
	conf = get_alembic_config()
	scriptdir = ScriptDirectory.from_config(conf)
	return scriptdir.get_current_head()


def current_revision(engine: Engine) -> str:
	"""Get the current revision number of a genome database."""
	with engine.connect() as conn:
		ctx = MigrationContext.configure(conn)
		return ctx.get_current_revision()


def is_current_revision(engine: Engine):
	"""Check if the current revision of a genome database is the most recent (head) revision."""
	head = current_head()
	current = current_revision(engine)
	return current == head


def upgrade(engine: Engine, revision: str = 'head', tag=None, **kwargs):
	"""Run the alembic upgrade command.

	See :func:`alembic.command.upgrade` for more information on how this works.

	Parameters
	----------
	engine
		SQLAlchemy engine specifying genome database connection info.
	revision
		Revision to upgrade to. Passed to :func:`alembic.command.upgrade`.
	tag
		Passed to :func:`alembic.command.upgrade`.
	\\**kwargs
		Passed to :func:`.get_alembic_config`.
	"""
	config = get_alembic_config(engine, **kwargs)
	command.upgrade(config, revision, tag=tag)


def init_db(engine: Engine):
	"""
	Initialize the genome database schema by creating all tables and stamping with the latest
	Alembic revision.

	Expects a fresh database that does not already contain any tables for the :mod:`midas.db.models`
	models and has not had any migrations run on it yet.

	Parameters
	----------
	engine : sqlalchemy.engine.base.Engine
		SQLAlchemy engine specifying database connection info.

	Raises
	------
	RuntimeError
		If the database is already stamped with an Alembic revision.
	sqlalchemy.exc.OperationalError
		If any of the database tables to be created already exist.
	"""
	from midas.db.models import Base

	conf = get_alembic_config()
	script = ScriptDirectory.from_config(conf)

	with engine.connect() as conn:
		ctx = MigrationContext.configure(conn)

		# Check there is no current revision stamped
		current = ctx.get_current_revision()
		if current is not None:
			raise RuntimeError(f'Expected uninitialized database, but current alembic revision is {current}')

		# Create tables
		# Set checkfirst=false so that we get an SQL error if any tables already exist
		Base.metadata.create_all(conn, checkfirst=False)

		# Stamp latest alembic version
		ctx.stamp(script, 'head')
