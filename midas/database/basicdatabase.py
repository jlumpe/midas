"""Defines a basic database implementation that fits in a single directory"""

import os
import shutil
import gzip
import re
from contextlib import suppress

import numpy as np
from sqlalchemy import Column, String, Binary
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy import select
from sqlalchemy.orm import sessionmaker, deferred, load_only

from midas.util import kwargs_done, SubPath
from midas.kmers import KmerCoordsCollection
from . import base


# SqlAlchemy declarative base
Base = declarative_base()


class GenomeSet(Base, base.GenomeSet):
	pass


class Genome(Base, base.Genome):
	pass


class GenomeAnnotations(Base, base.GenomeAnnotations):
	pass


class Sequence(Base, base.Sequence):
	_filename = Column(String(), nullable=False, unique=True)


class KmerSetCollection(Base, base.KmerSetCollection):
	pass


class KmerSet(Base, base.KmerSet):
	_data = deferred(Column(Binary()))


class BasicDatabase(base.AbstractDatabase):
	"""Basic database that resides in a single directory.

	:param str path: Path of database directory to open.
	:param dict engine_args: Additional keyword arguments to
		:func:`sqlalchemy.create_engine`.

	.. attribute:: path

		Path to database directory.
	"""

	__root_dir_attr__ = 'path'

	_seq_dir = SubPath('sequences')

	Base = Base
	Genome = Genome
	GenomeSet = GenomeSet
	GenomeAnnotations = GenomeAnnotations
	Sequence = Sequence
	KmerSetCollection = KmerSetCollection
	KmerSet = KmerSet

	def __init__(self, path, engine_args=dict()):
		self.path = os.path.abspath(path)
		self._engine = self._make_engine(self.path, **engine_args)
		self._Session = sessionmaker(bind=self._engine)

	@classmethod
	def create(cls, path):
		"""Creates and initializes a new BasicDatabase.

		:param str path: Path to create database at.
		:rtype: .BasicDatabase
		"""
		path = os.path.abspath(path)

		# Create directory if it does not exist
		os.makedirs(path, exist_ok=True)

		# Create subdirectories
		os.mkdir(cls._seq_dir(path))

		# Create database and tables
		engine = cls._make_engine(path)
		Base.metadata.create_all(bind=engine)

		return cls(path)

	@classmethod
	def _make_engine(cls, path, **kwargs):
		"""Create the engine for a BasicDatabase located at the given path"""
		return create_engine('sqlite:///{}/db.sqlite3'.format(path), **kwargs)

	def get_session(self):
		return self._Session()

	def store_sequence(self, genome_id, src, _session=None, **kwargs):

		# Get keyword arguments
		seq_format = kwargs.pop('format', 'fasta')
		src_compression = kwargs.pop('src_compression', None)
		keep_src = kwargs.pop('keep_src', True)
		src_mode = kwargs.pop('src_mode', 't')

		kwargs_done(kwargs)

		# Check compression argument
		if src_compression not in (None, 'gzip'):
			raise ValueError("Don't know how to open src with compression {}"
			                 .format(src_compression))

		# Create session context
		with self._optional_session(_session) as session:

			# Get genome
			genome = session.query(self.Genome).get(genome_id)

			# Get destination path
			seq_fname = self._make_seq_fname(genome, ext='.gz')
			dest_path = os.path.join(self._seq_dir, seq_fname)

			# Check sequence does not already exist
			if genome.sequence is not None:
				raise RuntimeError('Genome already has a sequence')

			# Check that there is no existing sequence with the same filename
			# (this shouldn't ever happen as they should be unique per
			# Genome)
			existing = session.query(self.Sequence).\
				filter_by(_filename=seq_fname).\
				scalar()
			if existing is not None:
				raise RuntimeError('Sequence with file name already exists')

			# Remove any existing file at the destination path as it must
			# have been left behind when its Sequence instance was deleted
			if os.path.isfile(dest_path):
				os.remove(dest_path)

			# Create Sequence instance
			sequence = Sequence(genome=genome, format=seq_format,
			                    _filename=seq_fname)
			session.add(sequence)

			# Copy/move src and commit database changes together - all-or-
			# nothing
			src_moved = False
			needs_delete = False
			try:

				# If passed a string, assume file path
				if isinstance(src, str):

					# Gzip compression can copy/move at file level
					if src_compression == 'gzip':
						if keep_src:
							shutil.copyfile(src, dest_path)
						else:
							shutil.move(src, dest_path)
							src_moved = True

					# Will need to convert to gzip
					else:
						with open(src, 'rb') as src_fh:
							with gzip.open(dest_path, 'wb') as dest_fh:
								shutil.copyfileobj(src_fh, dest_fh)

						needs_delete = not keep_src

				# Otherwise assume file-like object
				else:

					# Gzip compression can write directly, otherwise compress
					if src_compression == 'gzip':
						dest_fh = open(dest_path, 'wb')
					else:
						dest_fh = gzip.open(
							dest_path,
							'wb' if 'b' in src_mode else 'wt'
						)

					with dest_fh:
						shutil.copyfileobj(src, dest_fh)

				# Finally commit the database transaction
				session.commit()

			except:

				# Need to put src file back if it was moved
				if src_moved:
					shutil.move(dest_path, src)

				# Otherwise delete the copied file
				else:
					with suppress(FileNotFoundError):
						os.remove(dest_path)

				raise

			else:
				# Delete src if needed
				if needs_delete:
					os.remove(src)

	def open_sequence(self, genome_id, _session=None):
		with self._optional_session(_session) as session:
			sequence = session.query(self.Sequence).get(genome_id)
			if sequence is None:
				raise RuntimeError(
					'Genome with ID {} has no sequence'
					.format(genome_id)
				)
			path = self._get_seq_file_path(sequence)

		return gzip.open(path, 'rt')

	def store_kset_coords(self, collection_id, genome_id, coords, _session=None):
		with self._optional_session(_session, commit=True) as session:
			kcol = session.query(self.KmerSetCollection).get(collection_id)
			coords = coords.astype(kcol.kmerspec().coords_dtype, copy=False)

			kset = KmerSet(genome_id=genome_id, collection_id=collection_id,
			               count=len(coords), _data=coords.tobytes())
			session.add(kset)

	def load_kset_coords(self, collection_id, genome_id, _session=None):
		with self._optional_session(_session) as session:
			kcol = session.query(self.KmerSetCollection).get(collection_id)
			dtype = kcol.kmerspec().coords_dtype
			return self._load_kset_coords(collection_id, genome_id, session, dtype)

	def _load_kset_coords(self, cid, gid, session, dtype):
		data, = session.query(self.KmerSet._data)\
			.filter_by(collection_id=cid, genome_id=gid)\
			.one()

		return np.frombuffer(data, dtype=dtype)

	def load_kset_collection(self, collection_id, genome_ids, callback=None, _session=None):
		"""Load stored coordinates for multiple k-mer sets in a collection.

		:param int collection_id: ID of :class:`KmerSetCollection` to load
			coordinates for.
		:param genome_ids: IDs of :class:`Genome` to load coordinates for.
		:returns: K-mer set in coordinate array format.
		:rtype: midas.kmers.KmerCoordsCollection
		"""
		with self._optional_session(_session) as session:
			kcol = session.query(self.KmerSetCollection).get(collection_id)
			dtype = kcol.kmerspec().coords_dtype

			lengths = [
				session.query(self.KmerSet.count).\
					filter_by(collection_id=collection_id, genome_id=gid)\
					.one()[0]
				for gid in genome_ids
			]
			coords = KmerCoordsCollection.empty(lengths, dtype=dtype)

			for i, gid in enumerate(genome_ids):
				coords[i] = self._load_kset_coords(collection_id, gid, session, dtype)
				if callback is not None:
					callback(gid)

			return coords

	def _get_seq_file_path(self, sequence):
		"""Gets path to sequence file"""
		seq_fname = (sequence if isinstance(sequence, str)
		             else sequence._filename)
		return os.path.join(self._seq_dir, seq_fname)

	def _make_seq_fname(self, genome, ext=None):
		"""Create a unique and descriptive file name for a Sequence's data"""
		base = '{}_{}'.format(genome.id, genome.description.upper())
		fname = re.sub('[^A-Z0-9]+', '_', base[:32])
		if ext is not None:
			fname += ext
		return fname

	def clean_seq_files(self, dry_run=False):
		"""Removes orphaned sequence files from directory

		Orphaned files are created when a Sequence is deleted from the
		SQLAlchemy database in a way that does not trigger the after_delete
		ORM event (such as from a bulk delete query). These files shouldn't
		be harmful but could waste disk space.

		:param bool dry_run: If ``True`` don't actually remove any sequences,
			just return which would have been removed.
		:returns: List of cleaned file names.
		:rtype: list
		"""

		# Get file names of existing sequences
		session = self.get_session()
		results = session.execute(select([Sequence.__table__.c._filename]))
		existing_fnames = set(fname for fname, in results.fetchall())

		# Remove others
		cleaned = []
		for fname in os.listdir(self._seq_dir):
			if fname not in existing_fnames:
				if not dry_run:
					os.unlink(self._get_seq_file_path(fname))
				cleaned.append(fname)

		# Return list of cleaned files
		return cleaned

	def get_alembic_config(self, **kwargs):
		"""Get an alembic config object to perform migrations.

		:param \\**kwargs: Keyword arguments to pass to
			:meth:`alembic.config.Config.__init__`.
		:returns: Alembic config object for the ``BasicDatabase`` class with
			this instance's database connection info.
		:rtype: alembic.config.Config
		"""

		from alembic.config import Config
		from pkg_resources import resource_filename

		migrations_dir = resource_filename('midas.database', 'migrate')
		ini_path = os.path.join(migrations_dir, 'alembic.ini')
		script_path = os.path.join(migrations_dir, 'basicdatabase')

		config = Config(ini_path, **kwargs)

		config.set_main_option('script_location', script_path)
		config.attributes['engine'] = self._engine

		return config
