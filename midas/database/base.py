"""Defines abstract base classes for databases and SQLAlchemy models
"""

from abc import ABCMeta, abstractmethod, abstractproperty
import contextlib

from sqlalchemy import Table, ForeignKey
from sqlalchemy import Column, Integer, String, Boolean, DateTime, Enum
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy.ext.hybrid import hybrid_property

from .sqla import TrackChangesMixin, JsonType, MutableJsonDict


class AbstractDatabase(metaclass=ABCMeta):
	"""Abstract base class defining interface for all database types"""

	@abstractmethod
	def get_session(self):
		"""Create a new SQLAlchemy session"""
		pass

	def session_context(self):
		"""Creates context manager that closes session on exit"""
		return contextlib.closing(self.get_session())

	@abstractmethod
	def store_sequence(self, genome, src, **kwargs):
		"""Stores a new genome in the database"""
		pass

	@abstractmethod
	def open_sequence(self, sequence):
		"""Gets an open file handle/stream to the sequence for a stored genome"""
		pass

	@abstractmethod
	def store_kset_coords(self, genome, collection, coords):
		"""Store a k-mer set in coordinate format"""
		pass

	@abstractmethod
	def load_kset_coords(self, kset):
		"""Load stored coordinates for a k-mer set"""
		pass

	# Concrete model subclasses must be attributes on concrete database
	# subclasses
	Base = abstractproperty()
	Genome = abstractproperty()
	Sequence = abstractproperty()
	GenomeSet = abstractproperty()
	GenomeAnnotations = abstractproperty()
	KmerSetCollection = abstractproperty()
	KmerSet = abstractproperty()


class Genome(TrackChangesMixin):
	"""Base model for a reference genome queries can be run against

	Corresponds to a single assembly (one or more contigs, but at least
	partially assembled) from what should be a single sequencing run. The
	same organism may have several genome entries for it. Typically this
	will correspond directly to a record in Genbank (accession database).

	This model simply stores the metadata of the genome itself. Information
	on the actual sequence itself is stored in the sequence relationship.
	A genome has at most one sequence associated with it.
	"""

	__tablename__ = 'genomes'

	# Numeric PK
	id = Column(Integer(), primary_key=True)

	# Short, unique description
	description = Column(String(), nullable=False)

	# Whether the genome is completely assembled or has multiple contigs
	is_assembled = Column(Boolean(), nullable=False)

	# Optional metadata for genomes from genbank
	gb_db = Column(String()) # Database, e.g. "assembly"
	gb_id = Column(Integer()) # UID of genbank record
	gb_acc = Column(Integer()) # Accession number of genbank record
	gb_taxid = Column(Integer()) # Taxonomy ID

	# Optional summary (from EUtils ESummary) for record and associated
	# taxonomy entry, as mutable JSON types
	# TODO - really should be immutable
	gb_summary = Column(MutableJsonDict.as_mutable(JsonType))
	gb_tax_summary = Column(MutableJsonDict.as_mutable(JsonType))

	# Arbitrary metadata as mutable JSON dict
	meta = Column(MutableJsonDict.as_mutable(JsonType))

	# One-to-one relationship to Sequence
	@declared_attr
	def sequence(cls):
		return relationship('Sequence', uselist=False, backref='genome')

	def __repr__(self):
		return '<{}.{} id={} desc="{}">'.format(
			self.__module__,
			type(self).__name__,
			self.id,
			self.description,
		)


class Sequence:
	__tablename__ = 'sequences'

	@declared_attr
	def genome_id(cls):
		return Column(ForeignKey('genomes.id'), primary_key=True)

	# Format of sequence - e.g. fasta
	format = Column(String(), nullable=False)

	def __repr__(self):
		return '<{}.{} genome={}>'.format(
			self.__module__,
			type(self).__name__,
			self.genome,
		)


class GenomeSet:
	__tablename__ = 'genome_sets'

	id = Column(Integer(), primary_key=True)
	name = Column(String(), unique=True, nullable=False)
	description = Column(String())

	@declared_attr
	def genome_annotations(cls):
		return relationship('GenomeAnnotations', lazy=True)

	def __repr__(self):
		return '<{}.{} id={} name="{}">'.format(
			self.__module__,
			type(self).__name__,
			self.id,
			self.name,
		)


class GenomeAnnotations(TrackChangesMixin):
	__tablename__ = 'genome_annotations'

	@declared_attr
	def genome_id(cls):
		return Column(ForeignKey('genomes.id'), primary_key=True)

	@declared_attr
	def genome_set_id(cls):
		return Column(ForeignKey('genome_sets.id'), primary_key=True)

	@declared_attr
	def genome(cls):
		return relationship('Genome', backref='annotations')

	@declared_attr
	def genome_set(cls):
		return relationship('GenomeSet')

	# Taxonomy - species, subspecies, and strain. May not match original
	# Genbank annotations after curation
	tax_species = Column(String())
	tax_genus = Column(String())
	tax_strain = Column(String())

	def __repr__(self):
		return '<{}.{} genome_id={} desc="{}" set={}>'.format(
			self.__module__,
			type(self).__name__,
			self.genome_id,
			self.description,
			self.genome_set,
		)


class KmerSetCollection(TrackChangesMixin):
	"""A collection of k-mer counts/statistics for a set of genomes calculated
	with the same parameters.
	"""

	__tablename__ = 'kmer_collections'

	# Integer primary key
	id = Column(Integer(), primary_key=True)

	# Unique name
	name = Column(String(), nullable=False, unique=True)

	# Prefix - string of upper-case nucleotide codes
	prefix = Column(String(), nullable=False)

	# Number of nucleotides AFTER prefix
	k = Column(Integer(), nullable=False)

	# Additional parameters used to construct the set (if any).
	# Current reserved for future use.
	parameters = Column(MutableJsonDict.as_mutable(JsonType), nullable=False,
	                    default=dict())

	# Arbitrary metadata as mutable JSON dict
	meta = Column(MutableJsonDict.as_mutable(JsonType))

	def __repr__(self):
		return '<{}.{} id={} name="{}" k={} prefix="{}">'.format(
			self.__module__,
			type(self).__name__,
			self.id,
			self.name,
			self.k,
			self.prefix,
		)

	@property
	def coords_dtype(self):
		"""Smallest unsigned unteger numpy dtype that can store coordinates"""
		if self.k <= 4:
			return 'u1'
		elif self.k <= 8:
			return 'u2'
		elif self.k <= 12:
			return 'u4'
		elif self.k <= 16:
			return 'u8'
		else:
			return None

class KmerSet:

	__tablename__ = 'kmer_sets'

	@declared_attr
	def collection_id(cls):
		return Column(Integer(), ForeignKey('kmer_collections.id'),
		              primary_key=True)

	@declared_attr
	def genome_id(cls):
		return Column(Integer(), ForeignKey('genomes.id'),
		              primary_key=True)

	# Number of k-mers in set
	count = Column(Integer(), nullable=False)

	@declared_attr
	def collection(cls):
		return relationship(
			'KmerSetCollection',
			backref=backref('kmer_sets', lazy='dynamic',
			                cascade='all, delete-orphan')
		)

	@declared_attr
	def genome(cls):
		return relationship(
			'Genome',
			backref=backref('kmer_sets', lazy='dynamic',
			                cascade='all, delete-orphan')
		)

	def __repr__(self):
		return '<{}.{} {}, {}>'.format(
			self.__module__,
			type(self).__name__,
			repr(self.collection),
			repr(self.genome),
		)
