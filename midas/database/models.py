"""SQLAlchemy models for MIDAS database."""

from abc import ABCMeta

import sqlalchemy as sa
from sqlalchemy import Column, Integer, String, Boolean, Binary
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship, backref, deferred
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.ext.declarative import declarative_base, DeclarativeMeta
from sqlalchemy.ext.associationproxy import association_proxy

from midas.ncbi import SeqRecordBase
from .mixins import KeyMixin, SeqRecordMixin
from .sqla import JsonType, MutableJsonDict


__all__ = [
	'Genome',
	'GenomeSet',
	'AnnotatedGenome',
	'Taxon',
	'SignatureSet',
	'Signature',
]


class DeclarativeABCMeta(DeclarativeMeta, ABCMeta):
	"""Metaclass inheriting from both SQLAlchemy's DeclarativeMeta and ABCMeta

	This is needed because :class:`.Genome` uses a mixin which inherits from
	an abstract base class with ABCMeta as its metaclass, and all classes must
	have a metaclass that inherits from each of its base classes metaclasses.
	"""
	pass


# SqlAlchemy declarative base, using metaclass defined above
Base = declarative_base(metaclass=DeclarativeABCMeta)


class Genome(Base, SeqRecordMixin, KeyMixin):
	"""Base model for a reference genome queries can be run against.

	Corresponds to a single assembly (one or more contigs, but at least
	partially assembled) from what should be a single sequencing run. The
	same organism or strain may have several genome entries for it. Typically
	this will correspond directly to a record in Genbank (assembly database).

	The data on this model should primarily pertain to the sample and
	sequencing run itself. It would be updated if for example a better
	assembly was produced from the original raw data, however more advanced
	interpretation such as taxonomy assignments belong on an attached
	:class:`AnnotatedGenome` object.


	.. attribute:: id

		Integer primary key.

	.. attribute:: key

		String column. See :class:`midas.database.mixins.KeyMixin`.

	.. attribute:: key_version

		String column. See :class:`midas.database.mixins.KeyMixin`.

	.. attribute:: description

		String column. Short description. Recommended to be unique but this is
		not enforced.

	.. attribute:: is_assembled

		Boolean column. Whether the genome is completely assembled or has
		multiple contigs.

	.. attribute:: entrez_db

		String column. See :class:`midas.ncbi.SeqRecordBase`.

	.. attribute:: entrez_id

		String column. See :class:`midas.ncbi.SeqRecordBase`.

	.. attribute:: genbank_acc

		String column. See :class:`midas.ncbi.SeqRecordBase`.

	.. attribute:: refseq_acc

		String column. See :class:`midas.ncbi.SeqRecordBase`.

	.. attribute:: ncbi_taxid

		Integer column. Genbank taxonomy ID (as assigned in original NCBI
		record).

	.. attribute:: entrez_summary

		JSON column. Optional summary (from EUtils ESummary) for record and
		associated taxonomy entry.

	.. attribute:: extra

		JSON column. Additional arbitrary data.

	.. attribute:: annotations

		One-to-many relationship to :class:`.AnnotatedGenome`.
	"""

	__tablename__ = 'genomes'

	id = Column(Integer(), primary_key=True)
	description = Column(String(), nullable=False)
	is_assembled = Column(Boolean())
	ncbi_taxid = Column(Integer(), index=True)

	# TODO - really should be immutable
	entrez_summary = Column(MutableJsonDict.as_mutable(JsonType))
	extra = Column(MutableJsonDict.as_mutable(JsonType))

	annotations = relationship('AnnotatedGenome', lazy=True,
	                           cascade='all, delete-orphan')

	def __repr__(self):
		return '<{}.{}:{} "{}">'.format(
			self.__module__,
			type(self).__name__,
			self.id,
			self.description,
		)


class GenomeSet(Base, KeyMixin):
	"""A collection of genomes along with additional annotations on each.

	This will be used (among other things) to identify a set of genomes
	a query will be run against. Like Genomes they can be distributed via
	updates, in which case they should have a unique key and version number
	to identify them.

	.. attribute:: id

		Integer primary key.

	.. attribute:: key

		String column. See :class:`midas.database.mixins.KeyMixin`.

	.. attribute:: key_version

		String column, required. See :class:`midas.database.mixins.KeyMixin`.

	.. attribute:: name

		String column. Unique name.

	.. attribute:: description

		Text column. Optional description.

	.. attribute:: extra

		JSON column. Additional arbitrary data.

	.. attribute:: genomes

		Many-to-many relationship with :class:`.AnnotatedGenome`, annotated
		versions of genomes in this set.

	.. attribute:: base_genomes

		Unannotated :class:`Genome`\\ s in this set. Association proxy to the
		``genome`` relationship of members of :attr:`genome`.

	"""
	__tablename__ = 'genome_sets'

	id = Column(Integer(), primary_key=True)
	name = Column(String(), unique=True, nullable=False)
	description = Column(String())

	extra = Column(MutableJsonDict.as_mutable(JsonType))

	genomes = relationship('AnnotatedGenome', lazy='dynamic',
	                       cascade='all, delete-orphan')
	base_genomes = association_proxy('annotations', 'genome')

	def __repr__(self):
		return '<{}.{}:{} "{}"">'.format(
			self.__module__,
			type(self).__name__,
			self.id,
			self.name,
		)


# Association table between AnnoatedGenome and Taxon
annotations_tax_assoc = sa.Table(
	'annotations_tax_assoc',
	Base.metadata,
	Column('genome_id', primary_key=True),
	Column('genome_set_id', primary_key=True),
	Column(
		'taxon_id',
		ForeignKey('taxa.id', ondelete='CASCADE'),
		primary_key=True,
	),
	sa.ForeignKeyConstraint(
		['genome_id', 'genome_set_id'],
		['genome_annotations.genome_id', 'genome_annotations.genome_set_id'],
		ondelete='CASCADE',
	)
)


class AnnotatedGenome(Base, SeqRecordBase):
	"""A genome with additional annotations as part of a genome set.

	Technically is an association object connecting Genomes with GenomeSetss,
	but contains hybrid properties connecting to the Genome's attributes which
	effectively make it function and an extended Genome object.

	Additional annotations are optional, and the presence of this record in the
	database can simply be used to indicate that a Genome is contained in a
	GenomeSet. Mostly holds taxonomy information as that is frequently a result
	of additional analysis on the sequence.

	.. attribute:: genome_id

		Integer column, part of composite primary key. ID of :class:`.Genome`
		the annotations are or.

	.. attribute:: genome_set_id

		Integer column, part of composite primary key. ID of :class:`.GenomeSet`
		the annotations are under.

	.. attribute:: organism

		String column. Single string describing the organism. May be "Genus,
		species[, strain]" but could contain more specific information. Intended
		to be human- readable and shouldn't have any semantic meaning for the
		application (in contrast to the :attr:`taxa` relationship).

	.. attribute:: genome

		Many-to-one relationship to :class:`.Genome`.

	.. attribute:: genome_set

		Many-to-one relationship to :class:`.GenomeSet`.

	.. attribute:: taxa

		Many-to-many relationship to :class:`.Taxon`. The taxa this Genome has
		been assigned to under the associated GenomeSet. Should only link to the
		most specific taxa, as membership to all additional taxa in each taxon's
		linage is implied. Will typically contain only zero or one taxa.
		Additional non-proper (i.e. custom, not defined on GenBank, or
		paraphyletic) taxa may be present.

	.. attribute:: description

		Hybrid property connected to attribute on :attr:`genome`.

	.. attribute:: is_assembled

		Hybrid property connected to attribute on :attr:`genome`.

	.. attribute:: entrez_db

		Hybrid property connected to attribute on :attr:`genome`.

	.. attribute:: entrez_id

		Hybrid property connected to attribute on :attr:`genome`.

	.. attribute:: genbank_acc

		Hybrid property connected to attribute on :attr:`genome`.

	.. attribute:: refseq_acc

		Hybrid property connected to attribute on :attr:`genome`.

	.. attribute:: ncbi_taxid

		Hybrid property connected to attribute on :attr:`genome`.
	"""
	__tablename__ = 'genome_annotations'

	genome_id = Column(ForeignKey('genomes.id'), primary_key=True)
	genome_set_id = Column(ForeignKey('genome_sets.id'), primary_key=True)

	organism = Column(String())

	genome = relationship('Genome')
	genome_set = relationship('GenomeSet')
	taxa = relationship(
		'Taxon',
		secondary=annotations_tax_assoc,
		backref='genomes',
	)

	description = hybrid_property(lambda self: self.genome.description)
	is_assembled = hybrid_property(lambda self: self.genome.is_assembled)
	ncbi_taxid = hybrid_property(lambda self: self.genome.ncbi_taxid)
	entrez_db = hybrid_property(lambda self: self.genome.entrez_db)
	entrez_id = hybrid_property(lambda self: self.genome.entrez_id)
	genbank_acc = hybrid_property(lambda self: self.genome.genbank_acc)
	refseq_acc = hybrid_property(lambda self: self.genome.refseq_acc)

	def __repr__(self):
		return '<{}.{} {}:{}>'.format(
			self.__module__,
			type(self).__name__,
			self.genome_set_id,
			self.genome_id,
		)


class Taxon(Base):
	"""A taxon used for classifying AnnotatedGenomes.

	The taxa form a hierarchy, with each having a parent and zero or more
	children.

	.. attribute:: id

		Integer primary key.

	.. attribute:: name

		String column, required. The Taxon's scientific name (if any), e.g.
		"Escherichia coli", or otherwise any other unique descriptive name.

	.. attribute:: rank

		String column. Taxonomic rank, if any. Species, genus, family, etc.

	.. attribute:: description

		String column. Optional description of taxon. Probably blank unless it
		is a custom taxon not present in GenBank.

	.. attribute:: parent_id

		Integer column. ID of Taxon that is the direct parent of this one.

	.. attribute:: ncbi_id

		Integer column. Genbank taxonomy ID, if any.

	.. attribute:: entrez_summary

		JSON column. Entrez ESummary result of the corresponding entry in the
		NCBI taxonomy database.
	"""

	__tablename__ = 'taxa'

	id = Column(Integer(), primary_key=True)

	name = Column(String(), nullable=False, index=True)
	rank = Column(String(), index=True)
	description = Column(String())

	parent_id = Column(ForeignKey('taxa.id'))

	ncbi_id = Column(Integer(), index=True)
	entrez_summary = Column(MutableJsonDict.as_mutable(JsonType))

	parent = relationship('Taxon', remote_side=[id])
	children = relationship('Taxon')

	def lineage(self):
		"""Get sorted list of the Taxon's ancestors, starting with itself.

		:rtype: list
		"""
		l = list()
		taxon = self

		while taxon:
			l.insert(0, taxon)
			taxon = taxon.parent

		return l

	@classmethod
	def from_esummary(cls, summary):
		"""Create from Entrez ESummary result in JSON format.

		:param dict summary: JSON summary parsed to dict.
		:rtype: .Taxon
		"""
		return Taxon(
			name=summary['scientificname'],
			rank=summary['rank'],
			ncbi_id=summary['taxid'],
			entrez_summary=summary,
		)


class SignatureSet(Base):
	"""A set of k-mer signatures calculated with the same parameters.

	.. attribute:: id

		Integer primary key.

	.. attribute:: name

		String column, required. Unique name.

	.. attribute:: prefix

		String column, required. K-mer prefix to search for, string of upper-
		case nucleotide codes.

	.. attribute:: k

		Integer column, required. Number of nucleotides *after* prefix.

	.. attribute:: parameters

		JSON column. Additional parameters used to construct the set (if any).
		Currently reserved for future use.

	.. attribute:: extra

		JSON column. Additional arbitrary data.
	"""

	__tablename__ = 'signature_sets'

	id = Column(Integer(), primary_key=True)
	name = Column(String(), nullable=False, unique=True)
	prefix = Column(String(), nullable=False, index=True)
	k = Column(Integer(), nullable=False, index=True)

	parameters = Column(
		MutableJsonDict.as_mutable(JsonType),
		nullable=False,
		default=dict(),
	)

	extra = Column(MutableJsonDict.as_mutable(JsonType))

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
		"""Get the ``KmerSpec`` for these parameters.

		:rtype: :class:`midas.kmers.KmerSpec`
		"""
		from midas.kmers import KmerSpec
		return KmerSpec(k=self.k, prefix=self.prefix.encode('ascii'))

	@staticmethod
	def _coords_from_row(row):
		import numpy as np
		return np.frombuffer(row.data, dtype=row.dtype_str, count=row.count)

	def get_signature(self, genome_id, connection=None):
		"""Get the signature for a genome.

		:param int genome_id: ID of genome to get signature for.
		:returns: Signature in k-mer coordinate format, or None if no signature
			exists for the given genome.
		:rtype: numpy.ndarray
		"""

		if connection is None:
			connection = sa.inspect(self).session

		stmt = sa.select([Signature.count, Signature.dtype_str, Signature.data])\
			.where(Signature.set_id == self.id)\
			.where(Signature.genome_id == genome_id)

		row = connection.execute(stmt).first()

		if row is None:
			return None

		else:
			return self._coords_from_row(row)

	def store_signature(self, genome_id, coords, connection=None):
		"""Store the signature for a genome.

		:param int genome_id: ID of genome to store signature for.
		:param coords: k-mer signature in coordinate format.
		:type coords: numpy.ndarray
		"""

		if connection is None:
			connection = sa.inspect(self).session

		stmt = Signature.__table__.insert().values(
			set_id=self.id,
			genome_id=genome_id,
			count=coords.size,
			dtype_str=coords.dtype.str,
			data=coords.tobytes(),
		)

		connection.execute(stmt)

	def load_signatures(self, ids=None, connection=None):

		import numpy as np

		if connection is None:
			connection = sa.inspect(self).session


class Signature(Base):
	"""A stored k-mer signature for a genome.

	It's probably overkill to actually use the ORM to access this table, it
	should really be queried directly or through methods on
	:class:`.SignatureSet`.

	.. attribute:: set_id

		Integer column, part of composite primary key. ID of
		:class:`.SignatureSet` this signature belongs to.

	.. attribute:: genome_id

		Integer column, part of composite primary key. ID of :class:`.Genome`
		this signature is for.

	.. attribute:: count

		Integer column, required. Number of k-mers in signature.

	.. attribute:: dtype_str

		String(3) column, required. Numpy dtype of data in three-character
		string format specifying byte order, type, and precision. e.g. "<u4" for
		little-endian unsigned 32-bit integers.

	.. attribute:: data

		Binary column, required. Raw binary data for an array of the appropriate
		data type and length, containing the k-mer signature in coordinate
		format.

	.. attribute:: signature_set

		Many-to-one relationship to :class:`.SignatureSet`.

	.. attribute:: genome

		Many-to-one relationship to :class:`.Genome`.
	"""

	__tablename__ = 'signatures'

	set_id = Column(Integer(), ForeignKey('signature_sets.id'), primary_key=True)
	genome_id = Column(Integer(), ForeignKey('genomes.id'), primary_key=True)

	count = Column(Integer(), nullable=False)
	dtype_str = Column(String(2), nullable=False)
	data = deferred(Column(Binary(), nullable=False))

	signature_set = relationship(
		'SignatureSet',
		backref=backref(
			'signatures',
			lazy='dynamic',
			cascade='all, delete-orphan',
		),
	)

	genome = relationship(
		'Genome',
		backref=backref(
			'signatures',
			lazy='dynamic',
			cascade='all, delete-orphan',
		),
	)

	def __repr__(self):
		return '<{}.{} {}:{}>'.format(
			self.__module__,
			type(self).__name__,
			self.set_id,
			self.genome_id,
		)
