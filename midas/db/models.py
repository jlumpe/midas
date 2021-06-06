"""SQLAlchemy models for MIDAS database."""

import sqlalchemy as sa
from sqlalchemy import Column, Integer, String, Boolean
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.ext.declarative import declarative_base

from midas.ncbi import SeqRecordBase
from .mixins import KeyMixin, SeqRecordMixin
from .sqla import JsonType, MutableJsonDict


__all__ = [
	'Genome',
	'ReferenceGenomeSet',
	'AnnotatedGenome',
	'Taxon',
]


# SqlAlchemy declarative base
Base = declarative_base()


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


	Attributes
	----------
	id : int
		Integer primary key.
	key : str
		String column. See :class:`midas.db.mixins.KeyMixin`.
	version : str
		String column. See :class:`midas.db.mixins.KeyMixin`.
	description : Optional[str]
		String column. Short description. Recommended to be unique but this is
		not enforced.
	entrez_db : Optional[str]
		String column. See :class:`midas.ncbi.SeqRecordBase`.
	entrez_id : Optional[int]
		String column. See :class:`midas.ncbi.SeqRecordBase`.
	genbank_acc : Optional[str]
		String column. See :class:`midas.ncbi.SeqRecordBase`.
	refseq_acc : Optional[str]
		String column. See :class:`midas.ncbi.SeqRecordBase`.
	extra : Optional[dict]
		JSON column. Additional arbitrary metadata.
	annotations : Collection[.AnnotatedGenome]
		One-to-many relationship to :class:`.AnnotatedGenome`.
	"""

	__tablename__ = 'genomes'

	id = Column(Integer(), primary_key=True)
	description = Column(String(), nullable=False)

	# TODO - really should be immutable
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

	Attributes
	----------
	id : int
		Integer primary key.
	key : str
		String column. See :class:`midas.db.mixins.KeyMixin`.
	version : str
		String column, required. See :class:`midas.db.mixins.KeyMixin`.
	name : str
		String column. Unique name.
	description : Optional[str]
		Text column. Optional description.
	extra : Optional[dict]
		JSON column. Additional arbitrary data.
	genomes : Collection[.AnnotatedGenome]
		Many-to-many relationship with :class:`.AnnotatedGenome`, annotated
		versions of genomes in this set.
	base_genomes : Collection[.Genome]
		Unannotated :class:`Genome`\\ s in this set. Association proxy to the
		``genome`` relationship of members of :attr:`genome`.
	taxa : Collection[.Taxon]
		One-to-many relationship to :class:`.Taxon`. The taxa that form the
		classification system for this reference set.
	"""
	__tablename__ = 'reference_genome_sets'

	id = Column(Integer(), primary_key=True)
	name = Column(String(), unique=True, nullable=False)
	description = Column(String())

	extra = Column(MutableJsonDict.as_mutable(JsonType))

	genomes = relationship('AnnotatedGenome', lazy='dynamic',
	                       cascade='all, delete-orphan')
	base_genomes = relationship('Genome', secondary='genome_annotations',
	                            lazy='dynamic', viewonly=True)

	def __repr__(self):
		return '<{}:{} {!r}>'.format(
			type(self).__name__,
			self.id,
			self.name,
		)

	def root_taxa(self):
		"""Query for root taxa belonging to the set.

		Returns
		-------
		sqlalchemy.orm.query.Query
		"""
		return self.taxa.filter_by(parent=None)


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

	Attributes
	----------
	genome_id : int
		Integer column, part of composite primary key. ID of :class:`.Genome`
		the annotations are or.
	reference_set_id : int
		Integer column, part of composite primary key. ID of
		:class:`.ReferenceGenomeSet` the annotations are under.
	organism : str
		String column. Single string describing the organism. May be "Genus,
		species[, strain]" but could contain more specific information. Intended
		to be human- readable and shouldn't have any semantic meaning for the
		application (in contrast to the :attr:`taxa` relationship).
	primary_taxon_id : int
		Integer column. ID of primary :class:`Taxon` this genome is classified
		as.
	genome : .Genome
		Many-to-one relationship to :class:`.Genome`.
	reference_set : .ReferenceGenomeSet
		Many-to-one relationship to :class:`.ReferenceGenomeSet`.
	primary_taxon : .Taxon
		Many-to-one relationship to :class:`.Taxon`. The primary taxon this
		genome is classified as under the associated ReferenceGenomeSet. Should
		be the most specific and "regular" (ideally defined on NCBI) taxon this
		genome belongs to. There may be alternates, which go in the
		:attr:`alternate_taxa` attribute.
	alternate_taxa
		Many-to-many relationship to :class:`.Taxon`. The additional taxa this
		Genome has been assigned to under the associated ReferenceGenomeSet.
		Should only link to the most specific taxa, as membership to all
		additional taxa in each taxon's linage is implied. Will typically be
		empty, but may contain non-proper (i.e. custom, not defined on GenBank,
		or paraphyletic) taxa.
	key : str
		Hybrid property connected to attribute on :attr:`genome`.
	version : str
		Hybrid property connected to attribute on :attr:`genome`.
	description : Optional[str]
		Hybrid property connected to attribute on :attr:`genome`.
	entrez_db : Optional[str]
		Hybrid property connected to attribute on :attr:`genome`.
	entrez_id : Optional[int]
		Hybrid property connected to attribute on :attr:`genome`.
	genbank_acc : Optional[str]
		Hybrid property connected to attribute on :attr:`genome`.
	refseq_acc : Optional[str]
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

	genome = relationship('Genome', back_populates='annotations')
	reference_set = relationship('ReferenceGenomeSet', back_populates='genomes')
	primary_taxon = relationship(
		'Taxon',
		backref=backref('genomes_primary', lazy='dynamic')
	)
	additional_taxa = relationship(
		'Taxon',
		secondary=annotations_additional_tax_assoc,
		backref=backref('genomes_additional', lazy='dynamic'),
	)

	key = hybrid_property(lambda self: self.genome.key)
	version = hybrid_property(lambda self: self.genome.version)
	description = hybrid_property(lambda self: self.genome.description)
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

	Attributes
	----------
	id : int
		Integer primary key.
	name : str
		String column, required. The Taxon's scientific name (if any), e.g.
		"Escherichia coli", or otherwise any other unique descriptive name.
	rank : Optional[str]
		String column. Taxonomic rank, if any. Species, genus, family, etc.
	description : Optional[str]
		String column. Optional description of taxon. Probably blank unless it
		is a custom taxon not present in GenBank.
	distance_threshold : Optional[float]
		Float column. Maximum distance from a query genome to a reference genome
		in this taxon for the query to be classified within the taxon.
	report : Bool
		Boolean column. Whether to report this taxon directly as a match when
		producing a human-readable query result. Some custom taxa might need to
		be "hidden" from the user, in which case the value should be false. The
		application should then ascend the taxon's lineage and choose the first
		ancestor where this field is true. Defaults to true.
	extra : Optional[dict]
		JSON column. Additional arbitrary data.
	reference_set_id : int
		Integer column. ID of :class:`.ReferenceGenomeSet` the taxon belongs to.
	parent_id : int
		Integer column. ID of Taxon that is the direct parent of this one.
	ncbi_id : Optional[int]
		Integer column. Genbank taxonomy ID, if any.
	parent : Optional[.Taxon]
		Many-to-one relationship with :class:`.Taxon`, the parent of this taxon
		(if any).
	children : Collection[.Taxon]
		One-to-many relationship with :class:`.Taxon`, the children of this
		taxon.
	reference_set : .ReferenceGenomeSet
		Many-to-one relationship to :class:`.ReferenceGenomeSet`.
	genomes_primary : Collection[.AnnotatedGenome]
		One-to-many relationship with :class:`.AnnotatedGenome`, genomes which
		have this taxon as their primary taxon.
	genomes_additional : Collection[.AnnotatedGenome]
		Many-to-many relationship with :class:`.AnnotatedGenome`, genomes which
		have this taxon in their "additional" taxa.
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
	extra = Column(MutableJsonDict.as_mutable(JsonType))

	reference_set = relationship(
		'ReferenceGenomeSet',
		backref=backref('taxa', lazy='dynamic', cascade='all, delete-orphan')
	)
	parent = relationship('Taxon', remote_side=[id], backref=backref('children', lazy=True))

	def ancestors(self, incself=False):
		"""Iterate through the taxon's ancestors from bottom to top.

		Parameters
		----------
		incself : bool
			If True start with self, otherwise start with parent.

		Returns
		-------
		Iterable[.Taxon]
		"""
		taxon = self if incself else self.parent
		while taxon is not None:
			yield taxon
			taxon = taxon.parent

	def lineage(self):
		"""Get a sorted list of the taxon's ancestors from top to bottom, including itself.

		Returns
		-------
		list[.Taxon]
		"""
		l = list(self.ancestors(incself=True))
		l.reverse()
		return l

	def root(self):
		"""Get the root taxon of this taxon's tree.

		The set of taxa in a :class:`.ReferenceGenomeSet` will generally form
		a forest instead of a single tree, so there can be multiple root taxa.

		Returns self if the taxon has no parent.

		Returns
		-------
		.Taxon
		"""
		if self.parent is None:
			return self
		else:
			return self.parent.root()

	def isleaf(self):
		"""Check if the taxon is a leaf (has no children).

		Returns
		-------
		bool
		"""
		return not self.children

	def descendants(self, incself=False):
		"""Iterate through taxa all of the taxon's descendants (pre-order depth-first).

		Parameters
		----------
		incself : bool
			Yield self first.

		Returns
		-------
		Iterable[.Taxon]
		"""
		if incself:
			yield self
		for child in self.children:
			yield from child.descendants(incself=True)

	def leaves(self):
		"""Iterate through all leaves in the taxon's subtree.

		For leaf taxa this will just yield the taxon itself.

		Returns
		-------
		Iterable[.Taxon]
		"""
		if self.isleaf():
			yield self
		else:
			for child in self.children:
				yield from child.leaves()

	def print_tree(self, indent='  ', *, _depth=0):
		"""Print the taxon's subtree for debugging.

		Parameters
		---------
		indent : str
			String used to indent each level of descendants.
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
