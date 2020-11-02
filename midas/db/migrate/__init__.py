"""Perform database migrations with Alembic.

This package also contains all Alembic data files.
"""

from alembic.config import Config
from alembic import command
from pkg_resources import resource_filename


def get_alembic_config(engine=None, **kwargs):
	"""Get an alembic config object to perform migrations.

	Parameters
	----------
	engine : sqlalchemy.engine.base.Engine
		SQLAlchemy engine specifying database connection info (optional). Assigned to ``'engine'``
		key of :attr:`alembic.config.Config.attributes`.
	\\**kwargs
		Keyword arguments to pass to :meth:`alembic.config.Config.__init__`.

	Returns
	-------
	alembic.config.Config
		Alembic config object.
	"""
	ini_path = resource_filename(__name__, 'alembic.ini')
	script_path = resource_filename(__name__, 'alembic')

	config = Config(ini_path, **kwargs)

	config.set_main_option('script_location', script_path)

	if engine is not None:
		config.attributes['engine'] = engine

	return config


def upgrade(engine, revision='head', tag=None, **kwargs):
	"""Run the alembic upgrade command.

	See :func:`alembic.command.upgrade` for more inforamation on how this works.

	Parameters
	----------
	engine : sqlalchemy.engine.base.Engine
		SQLAlchemy engine specifying database connection info.
	revision : str
		Revision to uptrade to. Passed to :func:`alembic.command.upgrade`.
	tag
		Passed to :func:`alembic.command.upgrade`.
	\\**kwargs
		Passed to :func:`.get_alembic_config`.
	"""
	config = get_alembic_config(engine, **kwargs)
	command.upgrade(config, revision, tag=tag)
