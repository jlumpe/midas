"""Defines abstract base classes for databases and SQLAlchemy models
"""

from abc import ABCMeta, abstractmethod, abstractproperty
import contextlib

from sqlalchemy import Table, ForeignKey, UniqueConstraint
from sqlalchemy import Column, Integer, String, Boolean, DateTime, Enum
from sqlalchemy.orm import relationship, backref, deferred
from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.ext.associationproxy import association_proxy

from .sqla import TrackChangesMixin, JsonType, MutableJsonDict, JsonableMixin
from midas.kmers import KmerSpec


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
	def open_sequence(self, genome_id):
		"""Gets an open file handle/stream to the sequence for a stored genome"""
		pass

	@abstractmethod
	def store_kset_coords(self, collection_id, genome_id, coords):
		"""Store a k-mer set in coordinate format"""
		pass

	@abstractmethod
	def load_kset_coords(self, collection_id, genome_id):
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


class KeyMixin:
	"""Mixin that defines key/key version columns

	The "key" column is intended to be a universally unique key that can be
	used to identify objects across databases on different systems. This is
	primarily intended to be used for distributing database updates.
	It can be any arbitrary string but the recommended format is a filepath-
	like structure separated by forward slashes. This results in a
	heirarchical format that supports using namespaces to avoid key
	conflicts. Example: "genbank/assembly/GCF_00000000.0", which
	corresponds to a specific genome stored in the Genbank assembly database.

	The "key_version" column is the version the keyed object's metadata,
	according to whatever source defined the key. Used to determine when the
	metadata needs to be updated. Should be in the format defined by PEP 440
	(https://www.python.org/dev/peps/pep-0440/)
	"""
	key = Column(String(), index=True)
	key_version = Column(String())

	@classmethod
	def by_key(cls, session, key, version=None):
		query = session.query(cls).filter_by(key=key)

		if version is None:
			return query.order_by(cls.key_version.desc()).first()
		else:
			return query.filter_by(key_version=version).scalar()


class Genome(KeyMixin, TrackChangesMixin, JsonableMixin):
	"""Base model for a reference genome queries can be run against

	Corresponds to a single assembly (one or more contigs, but at least
	partially assembled) from what should be a single sequencing run. The
	same organism or strain may have several genome entries for it. Typically
	this will correspond directly to a record in Genbank (assembly database).

	This model simply stores the metadata of the genome, nformation	on the
	actual sequence data itself is stored in the sequence relationship.
	A Genome entry has at most one sequence associated with it.

	The data on this model should primarily pertain to the sample and
	sequencing run itself. It would be updated if for example a better
	assembly was produced from the original raw data, however more advanced
	interpretation such as taxonomy assignments belong on an attached
	GenomeAnnotations object.
	"""

	__tablename__ = 'genomes'
	__table_args__ = (
		UniqueConstraint('key', 'key_version'),
	)

	__json_attrs__ = [
		'key',
		'key_version',
		'description',
		'is_assembled',
		'gb_db',
		'gb_id',
		'gb_acc',
		'gb_taxid',
		'gb_summary',
		'gb_tax_summary',
		'meta',
	]

	# Integer PK
	id = Column(Integer(), primary_key=True)

	# Short description. Recommended to be unique but this is not enforced.
	description = Column(String(), nullable=False)

	# Whether the genome is completely assembled or has multiple contigs
	is_assembled = Column(Boolean())

	# Optional metadata for genomes from genbank
	gb_db = Column(String()) # Database, e.g. "assembly"
	gb_id = Column(Integer(), index=True) # UID of genbank record
	gb_acc = Column(Integer(), index=True) # Accession number of genbank record
	gb_taxid = Column(Integer(), index=True) # Taxonomy ID

	# Optional summary (from EUtils ESummary) for record and associated
	# taxonomy entry, as mutable JSON types
	# TODO - really should be immutable
	@declared_attr
	def gb_summary(cls):
		return deferred(Column(MutableJsonDict.as_mutable(JsonType)))

	@declared_attr
	def gb_tax_summary(cls):
		return deferred(Column(MutableJsonDict.as_mutable(JsonType)))

	# Arbitrary metadata as mutable JSON dict
	@declared_attr
	def	meta(cls):
		return deferred(Column(MutableJsonDict.as_mutable(JsonType)))

	# One-to-one relationship to Sequence
	@declared_attr
	def sequence(cls):
		return relationship('Sequence', uselist=False, backref='genome',
		                    cascade='all, delete-orphan')

	@declared_attr
	def annotations(cls):
		return relationship('GenomeAnnotations', lazy=True,
		                    cascade='all, delete-orphan')

	def __repr__(self):
		return '<{}.{}:{} "{}">'.format(
			self.__module__,
			type(self).__name__,
			self.id,
			self.description,
		)


class Sequence:
	"""Stores metadata for genome's sequence data, if it is in the database

	Genomes may be present in the database with associated metadata as a
	Genome object, but may not have a sequence stored. This contains
	metadata for the sequence itself.
	"""
	__tablename__ = 'sequences'

	@declared_attr
	def genome_id(cls):
		return Column(ForeignKey('genomes.id'), primary_key=True)

	# Format of sequence - e.g. fasta
	format = Column(String(), nullable=False)

	def __repr__(self):
		return '<{}.{} {}>'.format(
			self.__module__,
			type(self).__name__,
			self.genome_id,
		)


class GenomeSet(KeyMixin, JsonableMixin):
	"""A collection of genomes along with additional annotations on each

	This will be used (among other things) to identify a set of genomes
	a query will be run against. Like Genomes they can be distributed via
	updates, in which case they should have a unique key and version number
	to identify them.
	"""
	__tablename__ = 'genome_sets'
	__table_args__ = (
		UniqueConstraint('key', 'key_version'),
	)

	__json_attrs__ = [
		'key',
		'key_version',
		'name',
		'description',
		'meta',
	]

	# Integer PK
	id = Column(Integer(), primary_key=True)

	# Unique name
	name = Column(String(), unique=True, nullable=False)

	# Optional text description
	description = Column(String())

	# Arbitrary metadata as mutable JSON dict
	@declared_attr
	def meta(cls):
		return deferred(Column(MutableJsonDict.as_mutable(JsonType)))

	@declared_attr
	def annotations(cls):
		return relationship('GenomeAnnotations', lazy='dynamic',
		                    cascade='all, delete-orphan')

	@declared_attr
	def genomes(cls):
		return association_proxy('annotations', 'genome')


	def __repr__(self):
		return '<{}.{}:{} "{}"">'.format(
			self.__module__,
			type(self).__name__,
			self.id,
			self.name,
		)


class GenomeAnnotations(TrackChangesMixin, JsonableMixin):
	"""Association object connecting Genomes with GenomeSets

	Can simply indicate that a Genome is contained in a GenomeSet, but can
	also carry additional annotations for the genome that are different
	between sets. Mostly holds taxonomy information as that is frequently
	a result of additional analysis on the sequence.
	"""
	__tablename__ = 'genome_annotations'

	__json_attrs__ = [
		'organism',
		'tax_species',
		'tax_genus',
		'tax_strain',
	]

	@declared_attr
	def genome_id(cls):
		return Column(ForeignKey('genomes.id'), primary_key=True)

	@declared_attr
	def genome_set_id(cls):
		return Column(ForeignKey('genome_sets.id'), primary_key=True)

	@declared_attr
	def genome(cls):
		return relationship('Genome')

	@declared_attr
	def genome_set(cls):
		return relationship('GenomeSet')

	# Single string describing the organism. May be "Genus, species[, strain]"
	# but could contain more specific information. Intended to be human-
	# readable and shouldn't have any semantic meaning for the application
	# (in contrast to the following taxonomy attributes, which are used to
	# determine a match when querying against the set).
	organism = Column(String())

	# Taxonomy - species, subspecies, and strain. May not match original
	# Genbank annotations after curation
	tax_species = Column(String(), index=True)
	tax_genus = Column(String(), index=True)
	tax_strain = Column(String(), index=True)

	def __repr__(self):
		return '<{}.{} {}:{}>'.format(
			self.__module__,
			type(self).__name__,
			self.genome_set_id,
			self.genome_id,
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
	prefix = Column(String(), nullable=False, index=True)

	# Number of nucleotides AFTER prefix
	k = Column(Integer(), nullable=False, index=True)

	# Additional parameters used to construct the set (if any).
	# Currently reserved for future use.
	@declared_attr
	def	parameters(cls):
		return deferred(Column(MutableJsonDict.as_mutable(JsonType),
		                       nullable=False, default=dict()))

	# Arbitrary metadata as mutable JSON dict
	@declared_attr
	def	meta(cls):
		return deferred(Column(MutableJsonDict.as_mutable(JsonType)))

	def __repr__(self):
		return '<{}.{}:{} "{}" {}/{}>'.format(
			self.__module__,
			type(self).__name__,
			self.id,
			self.name,
			self.prefix,
			self.k,
		)

	def kmerspec(self):
		return KmerSpec(k=self.k, prefix=self.prefix.encode('ascii'))


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
		return '<{}.{} {}:{}>'.format(
			self.__module__,
			type(self).__name__,
			self.collection_id,
			self.genome_id,
		)
