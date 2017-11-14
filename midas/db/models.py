"""SQLAlchemy models for MIDAS database."""

from abc import ABCMeta

import sqlalchemy as sa
from sqlalchemy import Column, Integer, String, Boolean
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.ext.declarative import declarative_base, DeclarativeMeta

from midas.ncbi import SeqRecordBase
from .mixins import KeyMixin, SeqRecordMixin
from .sqla import JsonType, MutableJsonDict


__all__ = [
	'Genome',
	'ReferenceGenomeSet',
	'AnnotatedGenome',
	'Taxon',
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

		String column. See :class:`midas.db.mixins.KeyMixin`.

	.. attribute:: version

		String column. See :class:`midas.db.mixins.KeyMixin`.

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
		return '<{}:{} "{}">'.format(
			type(self).__name__,
			self.id,
			self.description,
		)


class ReferenceGenomeSet(Base, KeyMixin):
	"""
	A collection of genomes along with additional annotations that enable
	queries to be run against them.

	Additional annotations include taxonomy (individual :class:`.Taxon`
	entries as well as assignments to genomes) and thresholds for determining
	membership of queries within those taxa.

	Like Genomes they can be distributed via updates, in which case they should
	have a unique key and version number to identify them.

	.. attribute:: id

		Integer primary key.

	.. attribute:: key

		String column. See :class:`midas.db.mixins.KeyMixin`.

	.. attribute:: version

		String column, required. See :class:`midas.db.mixins.KeyMixin`.

	.. attribute:: name

		String column. Unique name.

	.. attribute:: description

		Text column. Optional description.

	.. attribute:: signatureset_key

		String column. Key of signature set that should be used to query this
		genome set.

	.. attribute:: signatureset_version

		String column. Key of signature set that should be used to query this
		genome set.

	.. attribute:: extra

		JSON column. Additional arbitrary data.

	.. attribute:: genomes

		Many-to-many relationship with :class:`.AnnotatedGenome`, annotated
		versions of genomes in this set.

	.. attribute:: base_genomes

		Unannotated :class:`Genome`\\ s in this set. Association proxy to the
		``genome`` relationship of members of :attr:`genome`.

	.. attribute:: taxa

		One-to-many relationship to :class:`.Taxon`. The taxa that form the
		classification system for this reference set.
	"""
	__tablename__ = 'reference_genome_sets'

	id = Column(Integer(), primary_key=True)
	name = Column(String(), unique=True, nullable=False)
	description = Column(String())

	signatureset_key = Column(String())
	signatureset_version = Column(String())

	extra = Column(MutableJsonDict.as_mutable(JsonType))

	genomes = relationship('AnnotatedGenome', lazy='dynamic',
	                       cascade='all, delete-orphan')
	base_genomes = relationship('Genome', secondary='genome_annotations',
	                            lazy='dynamic')

	def __repr__(self):
		return '<{}:{} {!r}>'.format(
			type(self).__name__,
			self.id,
			self.name,
		)


# Association table between AnnotatedGenome and Taxon
annotations_additional_tax_assoc = sa.Table(
	'annotations_additional_tax_assoc',
	Base.metadata,
	Column('genome_id', primary_key=True),
	Column('reference_set_id', primary_key=True),
	Column(
		'taxon_id',
		ForeignKey('taxa.id', ondelete='CASCADE'),
		primary_key=True,
	),
	sa.ForeignKeyConstraint(
		['genome_id', 'reference_set_id'],
		['genome_annotations.genome_id', 'genome_annotations.reference_set_id'],
		ondelete='CASCADE',
	)
)


class AnnotatedGenome(Base, SeqRecordBase):
	"""A genome with additional annotations as part of a genome set.

	Technically is an association object connecting Genomes with
	ReferenceGenomeSets, but contains hybrid properties connecting to the
	Genome's attributes which effectively make it function and an extended
	Genome object.

	Additional annotations are optional, and the presence of this record in the
	database can simply be used to indicate that a Genome is contained in a
	ReferenceGenomeSet. Mostly holds taxonomy information as that is frequently
	a result of additional analysis on the sequence.

	.. attribute:: genome_id

		Integer column, part of composite primary key. ID of :class:`.Genome`
		the annotations are or.

	.. attribute:: reference_set_id

		Integer column, part of composite primary key. ID of
		:class:`.ReferenceGenomeSet` the annotations are under.

	.. attribute:: organism

		String column. Single string describing the organism. May be "Genus,
		species[, strain]" but could contain more specific information. Intended
		to be human- readable and shouldn't have any semantic meaning for the
		application (in contrast to the :attr:`taxa` relationship).

	.. attribute:: primary_taxon_id

		Integer column. ID of primary :class:`Taxon` this genome is classified
		as.

	.. attribute:: genome

		Many-to-one relationship to :class:`.Genome`.

	.. attribute:: reference_set

		Many-to-one relationship to :class:`.ReferenceGenomeSet`.

	.. attribute:: primary_taxon

		Many-to-one relationship to :class:`.Taxon`. The primary taxon this
		genome is classified as under the associated ReferenceGenomeSet. Should
		be the most specific and "regular" (ideally defined on NCBI) taxon this
		genome belongs to. There may be alternates, which go in the
		:attr:`alternate_taxa` attribute.

	.. attribute:: alternate_taxa

		Many-to-many relationship to :class:`.Taxon`. The additional taxa this
		Genome has been assigned to under the associated ReferenceGenomeSet.
		Should only link to the most specific taxa, as membership to all
		additional taxa in each taxon's linage is implied. Will typically be
		empty, but may contain non-proper (i.e. custom, not defined on GenBank,
		or paraphyletic) taxa.

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

	genome_id = Column(
		ForeignKey('genomes.id', ondelete='CASCADE'),
		primary_key=True
	)
	reference_set_id = Column(
		ForeignKey('reference_genome_sets.id', ondelete='CASCADE'),
		primary_key=True
	)

	primary_taxon_id = Column(
		ForeignKey('taxa.id', ondelete='SET NULL'),
		index=True
	)
	organism = Column(String())

	genome = relationship('Genome')
	reference_set = relationship('ReferenceGenomeSet')
	primary_taxon = relationship(
		'Taxon',
		backref=backref('genomes_primary', lazy='dynamic')
	)
	additional_taxa = relationship(
		'Taxon',
		secondary=annotations_additional_tax_assoc,
		backref=backref('genomes_additional', lazy='dynamic'),
	)

	description = hybrid_property(lambda self: self.genome.description)
	is_assembled = hybrid_property(lambda self: self.genome.is_assembled)
	ncbi_taxid = hybrid_property(lambda self: self.genome.ncbi_taxid)
	entrez_db = hybrid_property(lambda self: self.genome.entrez_db)
	entrez_id = hybrid_property(lambda self: self.genome.entrez_id)
	genbank_acc = hybrid_property(lambda self: self.genome.genbank_acc)
	refseq_acc = hybrid_property(lambda self: self.genome.refseq_acc)

	def __repr__(self):
		return '<{}:{}:{} {!r}/{!r}>'.format(
			type(self).__name__,
			self.reference_set_id,
			self.genome_id,
			self.reference_set.key,
			self.organism,
		)


class Taxon(Base):
	"""A taxon used for classifying AnnotatedGenomes.

	Taxa are specific to a :class:`.ReferenceGenomeSet` and form a hierarchy,
	with each having a parent and zero or more children.

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

	.. attribute:: distance_threshold

		Float column. Maximum distance from a query genome to a reference genome
		in this taxon for the query to be classified within the taxon.

	.. attribute:: report

		Boolean column. Whether to report this taxon directly as a match when
		producing a human-readable query result. Some custom taxa might need to
		be "hidden" from the user, in which case the value should be false. The
		application should then ascend the taxon's lineage and choose the first
		ancestor where this field is true. Defaults to true.

	.. attribute:: extra

		JSON column. Additional arbitrary data.

	.. attribute:: reference_set_id

		Integer column. ID of :class:`.ReferenceGenomeSet` the taxon belongs to.

	.. attribute:: parent_id

		Integer column. ID of Taxon that is the direct parent of this one.

	.. attribute:: ncbi_id

		Integer column. Genbank taxonomy ID, if any.

	.. attribute:: entrez_data

		JSON column. Entrez EFetch result of the corresponding entry in the
		NCBI taxonomy database, converted to JSON.

	.. atribute:: parent

		Many-to-one relationship with :class:`.Taxon`, the parent of this taxon
		(if any).

	.. attribute:: children

		One-to-many relationship with :class:`.Taxon`, the children of this
		taxon.

	.. attribute:: reference_set

		Many-to-one relationship to :class:`.ReferenceGenomeSet`.

	.. attribute:: genomes_primary

		One-to-many relationship with :class:`.AnnotatedGenome`, genomes which
		havet his taxon as their primary taxon.

	.. attribute:: genomes_additional

		Many-to-many relationship with :class:`.AnnotatedGenome`, genomes which
		havet his taxon in their "additional" taxa.
	"""

	__tablename__ = 'taxa'

	id = Column(Integer(), primary_key=True)

	name = Column(String(), nullable=False, index=True)
	rank = Column(String(), index=True)
	description = Column(String())
	distance_threshold = Column(sa.Float())
	report = Column(Boolean(), nullable=False, default=True, server_default=sa.true())

	reference_set_id = Column(
		ForeignKey('reference_genome_sets.id', ondelete='CASCADE'),
		nullable=False,
		index=True,
	)
	parent_id = Column(ForeignKey('taxa.id', ondelete='SET NULL'), index=True)

	ncbi_id = Column(Integer(), index=True)
	entrez_data = Column(MutableJsonDict.as_mutable(JsonType))
	extra = Column(MutableJsonDict.as_mutable(JsonType))

	reference_set = relationship(
		'ReferenceGenomeSet',
		backref=backref('taxa', lazy='dynamic', cascade='all, delete-orphan')
	)
	parent = relationship('Taxon', remote_side=[id])
	children = relationship('Taxon')

	def lineage(self):
		"""Get sorted list of the Taxon's ancestors, including itself.

		:returns: List of taxa ordered from parent to child.
		:rtype: list
		"""
		l = list()
		taxon = self

		while taxon:
			l.insert(0, taxon)
			taxon = taxon.parent

		return l

	def root(self):
		"""Get the root taxon of this taxon's tree.

		The set of taxa in a :class:`.ReferenceGenomeSet` will generally form
		a forest instead of a single tree, so there can be multiple root taxa.

		Returns self if the taxon has no parent.

		:rtype: .Taxon
		"""
		if self.parent is None:
			return self
		else:
			return self.parent.root()

	def print_tree(self, indent='  ', *, _depth=0):
		"""Print the taxon's subtree for debugging.

		:param str indent: String used to indent each level of descendants.
		"""
		print(indent * _depth + self.name)
		for child in sorted(self.children, key=lambda c: c.name):
			child.print_tree(indent=indent, _depth=_depth + 1)

	def __repr__(self):
		return '<{}:{} {!r}>'.format(
			type(self).__name__,
			self.id,
			self.name,
		)