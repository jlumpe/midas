"""SQLAlchemy models for storing reference genomes and taxonomy information."""

from typing import Sequence, Union, Dict, List, Any, Optional, Tuple

import sqlalchemy as sa
from sqlalchemy import Column, Integer, String, Boolean, Float
from sqlalchemy import Table, ForeignKey, UniqueConstraint, ForeignKeyConstraint
from sqlalchemy.orm import relationship, backref, deferred
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.ext.declarative import declarative_base, declared_attr
from sqlalchemy.orm.attributes import InstrumentedAttribute

from .sqla import JsonString


__all__ = [
	'Genome',
	'ReferenceGenomeSet',
	'AnnotatedGenome',
	'Taxon',
]


# Naming convention for constraints and indices, used by SQLAlchemy when creating schema.
# Important for Alembic migration scripts, see https://alembic.sqlalchemy.org/en/latest/naming.html
NAMING_CONVENTION = {
  "ix": "ix_%(column_0_label)s",
  "uq": "uq_%(table_name)s_%(column_0_name)s",
  "ck": "ck_%(table_name)s_%(constraint_name)s",
  "fk": "fk_%(table_name)s_%(column_0_name)s_%(referred_table_name)s",
  "pk": "pk_%(table_name)s",
}

# SqlAlchemy metadata object and declarative base
metadata = sa.MetaData(naming_convention=NAMING_CONVENTION)
Base = declarative_base(metadata=metadata)


class KeyMixin:
	"""Mixin that defines key/version columns.

	Attributes
	----------
	key : str
		Intended to be a universally unique key that can be used to identify objects across
		databases on different systems. This is primarily intended to be used for distributing
		database updates. It can be any arbitrary string but the recommended format is a
		filepath-like structure separated by forward slashes. This results in a hierarchical
		format that supports using namespaces to avoid key conflicts. Example:
		``'ncbi/assembly/GCF_00000000.0'``, which corresponds to a specific genome stored in the
		NCBI assembly database.
	version : str
		Version of the keyed object according to whatever source defined the key. Used to
		determine when the metadata needs to be updated. Should be in the format defined by `PEP
		440 <https://www.python.org/dev/peps/pep-0440/>`_.
	"""
	key = Column(String(), index=True)
	version = Column(String())

	@declared_attr
	def __table_args__(cls):
		return (
			UniqueConstraint('key', 'version'),
		)

	@classmethod
	def by_key(cls, session, key, version=None):
		query = session.query(cls).filter_by(key=key)

		if version is None:
			return query.order_by(cls.version.desc()).first()
		else:
			return query.filter_by(version=version).scalar()


class Genome(Base, KeyMixin):
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
		String column. Short one-line description. Recommended to be unique but this is not enforced.
	ncbi_db : Optional[str]
		String column (optional). If the genome corresponds to a record downloaded from an NCBI
		database this column should be the database name (e.g. ``'assembly'``) and ``ncbi_id``
		should be the entry's UID. Unique along with ``ncbi_id``.
	ncbi_id : Optional[int]
		Integer column (optional). See previous.
	genbank_acc : Optional[str]
		String column. GenBank accession number for this genome, if any.
	refseq_acc : Optional[str]
		String column. RefSeq accession number for this genome, if any.
	extra : Optional[dict]
		JSON column. Additional arbitrary metadata.
	annotations : Collection[.AnnotatedGenome]
		One-to-many relationship to :class:`.AnnotatedGenome`.
	"""

	__tablename__ = 'genomes'

	@declared_attr
	def __table_args__(cls):
		return (
			UniqueConstraint('ncbi_db', 'ncbi_id'),
		)

	id = Column(Integer(), primary_key=True)
	description = Column(String(), nullable=False)
	ncbi_db = Column(String())
	ncbi_id = Column(Integer())
	genbank_acc = Column(String(), unique=True)
	refseq_acc = Column(String(), unique=True)
	extra = deferred(Column(JsonString()))

	annotations = relationship('AnnotatedGenome', lazy=True, cascade='all, delete-orphan')

	def __repr__(self):
		return f'<{type(self).__name__}:{self.id} {self.description!r}>'


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
	extra = Column(JsonString())

	genomes = relationship('AnnotatedGenome', lazy='dynamic', cascade='all, delete-orphan')
	base_genomes = relationship('Genome', secondary='genome_annotations', lazy='dynamic', viewonly=True)

	def __repr__(self):
		return f'<{type(self).__name__}:{self.id} {self.name!r}>'

	def root_taxa(self):
		"""Query for root taxa belonging to the set.

		Returns
		-------
		sqlalchemy.orm.query.Query
		"""
		return self.taxa.filter_by(parent=None)


# Association table between AnnotatedGenome and Taxon
annotations_additional_tax_assoc = Table(
	'annotations_additional_tax_assoc',
	Base.metadata,
	Column('genome_id', primary_key=True),
	Column('reference_set_id', primary_key=True),
	Column(
		'taxon_id',
		ForeignKey('taxa.id', ondelete='CASCADE'),
		primary_key=True,
	),
	ForeignKeyConstraint(
		['genome_id', 'reference_set_id'],
		['genome_annotations.genome_id', 'genome_annotations.reference_set_id'],
		ondelete='CASCADE',
	)
)


class AnnotatedGenome(Base):
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
	ncbi_db : Optional[str]
		Hybrid property connected to attribute on :attr:`genome`.
	ncbi_id : Optional[int]
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
	ncbi_db = hybrid_property(lambda self: self.genome.ncbi_db)
	ncbi_id = hybrid_property(lambda self: self.genome.ncbi_id)
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
	distance_threshold = Column(Float())
	report = Column(Boolean(), nullable=False, default=True, server_default=sa.true())

	reference_set_id = Column(
		ForeignKey('reference_genome_sets.id', ondelete='CASCADE'),
		nullable=False,
		index=True,
	)
	parent_id = Column(ForeignKey('taxa.id', ondelete='SET NULL'), index=True)

	ncbi_id = Column(Integer(), index=True)
	extra = deferred(Column(JsonString()))

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
		return f'<{type(self).__name__}:{self.id} {self.name!r}>'


#: Attributes of :class:`midas.db.models.Genome` which serve as unique IDs.
GENOME_ID_ATTRS = {
	'key': Genome.key,
	'genbank_acc': Genome.genbank_acc,
	'refseq_acc': Genome.refseq_acc,
	'ncbi_id': Genome.ncbi_id,
}


# Type alias for argument specifying genome id attribute
GenomeAttr = Union[str, InstrumentedAttribute]


def _check_genome_id_attr(attr: GenomeAttr) -> InstrumentedAttribute:
	"""Check that Genome ID attribute is valid, and convert from string argument.
	"""
	if isinstance(attr, str):
		try:
			return GENOME_ID_ATTRS[attr]
		except KeyError:
			pass

	elif isinstance(attr, InstrumentedAttribute):
		for allowed in GENOME_ID_ATTRS.values():
			if attr is allowed:
				return attr

	raise ValueError('Genome ID attribute must be one of the following: ' + ', '.join(GENOME_ID_ATTRS))


def _get_genome_id(genome: Union[Genome, AnnotatedGenome], attr: InstrumentedAttribute):
	"""Get value of ID attribute for genome."""
	if isinstance(genome, AnnotatedGenome):
		genome = genome.genome
	return attr.__get__(genome, Genome)


def _check_genomes_have_ids(genomeset: ReferenceGenomeSet, id_attr: InstrumentedAttribute):
	"""Check all genomes in ReferenceGenomeSet have values for the given ID attribute or raise a ``RuntimeError``."""
	c = genomeset.genomes \
		.join(AnnotatedGenome.genome) \
		.filter(id_attr == None) \
		.count()

	if c > 0:
		raise RuntimeError(f'{c} genomes missing value for ID attribute {id_attr.key}')


def _map_ids_to_genomes(genomeset: ReferenceGenomeSet, id_attr: Union[str, InstrumentedAttribute]) -> Dict[AnnotatedGenome, Any]:
	"""Get dict mapping ID values to AnnotatedGenome."""
	q = genomeset.genomes.join(AnnotatedGenome.genome).add_columns(id_attr)
	return {id_: g for g, id_ in q}


def genomes_by_id(genomeset: ReferenceGenomeSet, id_attr: GenomeAttr, ids: Sequence, strict: bool = True) -> List[Optional[AnnotatedGenome]]:
	"""Match a :class:`ReferenceGenomeSet`'s genomes to a set of ID values.

	This is primarily used to match genomes to signatures based on the ID values stored in a
	signature file. It is expected that the signature file may contain signatures for more genomes
	than are present in the genome set, see also :func:`.genomes_by_id_subset` for that condition.

	Parameters
	----------
	genomeset
	id_attr
		ID attribute of :class:`midas.db.models.Genome` to use for lookup. Can be used as the
		attribute itself (e.g. ``Genome.refseq_acc``) or just the name (``'refsec_acc'``).
		See :data:`.GENOME_IDS` for the set of allowed values.
	ids
		Sequence of ID values (strings or integers, matching type of attribute).
	strict
		Raise an exception if a matching genome cannot be found for any ID value.

	Returns
	-------
	List[Optional[AnnotatedGenome]]
		List of genomes of same length as ``ids``. If ``strict=False`` and a genome cannot be found
		for a given ID the list will contain ``None`` at the corresponding position.

	Raises
	------
	KeyError
		If ``strict=True`` and any ID value cannot be found.
	"""
	id_attr = _check_genome_id_attr(id_attr)
	_check_genomes_have_ids(genomeset, id_attr)
	d = _map_ids_to_genomes(genomeset, id_attr)
	if strict:
		return [d[id_] for id_ in ids]
	else:
		return [d.get(id_) for id_ in ids]


def genomes_by_id_subset(genomeset: ReferenceGenomeSet,
                         id_attr: GenomeAttr,
                         ids: Sequence,
                         ) -> Tuple[List[AnnotatedGenome], List[int]]:
	"""Match a :class:`ReferenceGenomeSet`'s genomes to a set of ID values, allowing missing genomes.

	This calls :func:`.genomes_by_id` with ``strict=False`` and filters any ``None`` values from the
	output. The filtered list is returned along with the indices of all values in ``ids`` which were
	not filtered out. The indices can be used to load only those signatures which have a matched
	genome from a signature file.

	Parameters
	----------
	genomeset
	id_attr
		ID attribute of :class:`midas.db.models.Genome` to use for lookup. Can be used as the
		attribute itself (e.g. ``Genome.refseq_acc``) or just the name (``'refsec_acc'``).
		See :data:`.GENOME_IDS` for the set of allowed values.
	ids
		Sequence of ID values (strings or integers, matching type of attribute).

	Returns
	-------
	Tuple[List[AnnotatedGenome], List[int]]
	"""
	genomes = genomes_by_id(genomeset, id_attr, ids, strict=False)
	genomes_out = []
	idxs_out = []

	for i, g in enumerate(genomes):
		if g is not None:
			genomes_out.append(g)
			idxs_out.append(i)

	return genomes_out, idxs_out
