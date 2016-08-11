
import os
import shutil
import gzip
import re
from contextlib import suppress

import numpy as np

from sqlalchemy import Column, String, Binary
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine, event
from sqlalchemy.orm import sessionmaker, deferred
from sqlalchemy.orm.session import Session, object_session

from . import base



class BasicDatabaseSession(Session):
	"""SQLAlchemy session that tracks the database it belongs to"""

	def __init__(self, *, midasdb, **kwargs):
		super(BasicDatabaseSession, self).__init__(**kwargs)
		self._midasdb = midasdb


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

	@classmethod
	def _after_delete(cls, mapper, connection, instance):
		"""Listener for deletion event

		Informs associated BasicDatabase that the sequence instance has been
		deleted so that its data file can be deleted as well.
		"""
		session = object_session(instance)
		if isinstance(session, BasicDatabaseSession):
			session._midasdb._after_sequence_delete(mapper, connection, instance)

	@classmethod
	def __declare_last__(cls):
		"""Register event listener"""
		event.listen(cls, 'after_delete', cls._after_delete)


class KmerSetCollection(Base, base.KmerSetCollection):
	pass


class KmerSet(Base, base.KmerSet):
	_data = deferred(Column(Binary()))


class BasicDatabase(base.AbstractDatabase):
	"""Basic database that resides in a single directory"""

	_seq_dir = 'sequences'

	def __init__(self, path):
		self.path = os.path.abspath(path)
		self._engine = self._make_engine(self.path)
		self._Session = sessionmaker(bind=self._engine,
		                             class_=BasicDatabaseSession,
		                             midasdb=self)

	@classmethod
	def create(cls, path):
		path = os.path.abspath(path)

		# Create directory if it does not exist
		os.makedirs(path, exist_ok=True)

		# Create subdirectories
		os.mkdir(os.path.join(path, cls._seq_dir))

		# Create database and tables
		engine = cls._make_engine(path)
		Base.metadata.create_all(bind=engine)

		return cls(path)

	@classmethod
	def _make_engine(cls, path):
		return create_engine('sqlite:///{}/db.sqlite3'.format(path))

	def get_session(self):
		return self._Session()

	def store_sequence(self, genome, src, **kwargs):

		# Get keyword arguments
		seq_format = kwargs.pop('format', 'fasta')
		src_compression = kwargs.pop('src_compression', None)
		keep_src = kwargs.pop('keep_src', True)
		src_mode = kwargs.pop('src_mode', 't')
		if kwargs:
			raise TypeError('Invalid keyword argument "{}"'
			                .format(next(iter(kwargs))))

		# Check compression argument
		if src_compression not in (None, 'gzip'):
			raise ValueError("Don't know how to open src with compression {}"
			                 .format(src_compression))
		
		# Get destination path
		seq_fname = self._make_seq_fname(genome) + '.gz'
		dest_path = os.path.join(self.path, self._seq_dir, seq_fname)

		# Create session context
		with self.session_context() as session:

			# Merge genome instance into current session
			merged_genome = session.merge(genome)

			# Check sequence does not already exist
			if merged_genome.sequence is not None:
				raise RuntimeError('Genome already has a sequence')

			# Create Sequence instance
			sequence = Sequence(genome=merged_genome, format=seq_format,
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
						dest_fh = gzip.open(dest_path,
						                    'wb' if 'b' in src_mode else 'wt')

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

	def open_sequence(self, sequence):
		seq_path = os.path.join(self.path, self._seq_dir,sequence._filename)
		return gzip.open(seq_path, 'rt')

	def store_kset_coords(self, genome, collection, coords):
		coords = coords.astype(collection.coords_dtype, copy=False)

		kset = KmerSet(genome_id=genome.id, collection_id=collection.id,
		               count=len(coords), _data=coords.tobytes())

		session = self.get_session()
		session.add(kset)
		session.commit()

	def load_kset_coords(self, kset):
		return np.frombuffer(kset._data, dtype=kset.collection.coords_dtype)

	def _after_sequence_delete(self, mapper, connection, target):
		"""Called by after_delete listener on Sequence"""
		os.remove(os.path.join(self.path, self._seq_dir, target._filename))

	@classmethod
	def _make_seq_fname(cls, genome):
		"""Create a unique and descriptive file name for a Sequence's data"""
		return '{}_{}'.format(
			genome.id,
			re.sub('[^A-Z0-9]+', '_', genome.description.upper())
		)

	Base = Base
	Genome = Genome
	GenomeSet = GenomeSet
	GenomeAnnotations = GenomeAnnotations
	Sequence = Sequence
	KmerSetCollection = KmerSetCollection
	KmerSet = KmerSet
