"""Perform database migrations with Alembic.

This package also contains all Alembic data files.
"""
import os

from alembic.config import Config
from pkg_resources import resource_filename


def get_alembic_config(engine, **kwargs):
	"""Get an alembic config object to perform migrations.

	:param engine: SQLAlchemy engine specifying database connection info.
	:type engine: sqlalchemy.engine.base.Engine
	:param \\**kwargs: Keyword arguments to pass to
		:meth:`alembic.config.Config.__init__`.
	:returns: Alembic config object for the engine.
	:rtype: alembic.config.Config
	"""
	ini_path = resource_filename(__name__, 'alembic.ini')
	script_path = resource_filename(__name__, 'alembic')

	config = Config(ini_path, **kwargs)

	config.set_main_option('script_location', script_path)
	config.attributes['engine'] = engine

	return config
