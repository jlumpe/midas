"""Test midas.db.models.

Uses the included testdb_210126 database.
"""

import random

import pytest
from sqlalchemy.orm import sessionmaker

from midas.db import models
from midas.db.models import Genome, ReferenceGenomeSet, AnnotatedGenome, Taxon


# Some arbitrary JSON data
JSON_DATA = {
	'int': 1,
	'float': 3.14,
	'string': 'foo',
	'bool': True,
	'null': None,
	'array': [1, 3.14, 'foo', True, None, [1, 2, 3, 4], {'foo': 'bar', 'baz': 1}],
	'object': {
		'foo': 'bar',
		'baz': 1,
		'sub_list': [1, 2, 3],
		'sub_object': {'blah': 123}
	},
}

JSON_DATA2 = JSON_DATA['object']


@pytest.fixture()
def empty_db_session(make_empty_db):
	"""Session factory for empty in-memory database."""
	engine = make_empty_db()
	return sessionmaker(engine)


@pytest.fixture()
def id_lookup_session(make_empty_db):
	"""
	Session for an in-memory database containing a set of genomes which have values for all ID
	values, to test lookup/matching funcs.
	"""
	engine = make_empty_db()
	Session = sessionmaker(engine)
	session = Session()

	gset = ReferenceGenomeSet(
		key='test_gset',
		version='1.0',
		name='Test genome set',
	)
	session.add(gset)

	roottaxon = Taxon(
		name='root',
		genome_set=gset,
	)
	session.add(roottaxon)

	for i in range(20):
		g = Genome(
			key=f'test/genome_{i}',
			version='1.0',
			description=f'Test genome {i}',
			ncbi_db='assembly',
			ncbi_id=i,
			genbank_acc=f'GCA_{i:09d}.1',
			refseq_acc=f'GCF_{i:09d}.1',
		)
		ag = AnnotatedGenome(
			genome_set=gset,
			genome=g,
			taxon=roottaxon,
		)
		session.add(ag)

	session.commit()
	return session


class TestGenome:

	def test_key_lookup(self, id_lookup_session):
		session = id_lookup_session
		key = 'test/genome_0'
		genome = session.query(Genome).filter_by(key=key).one()
		assert Genome.by_key(session, genome.key) == genome
		assert Genome.by_key(session, genome.key, '1.0') == genome
		assert Genome.by_key(session, genome.key, '1.1') is None

	def test_ncbi_id_lookup(self, id_lookup_session):
		pass  # TODO

	def test_extra_json(self, empty_db_session):
		"""Test storing JSON data in the 'extra' column."""
		session = empty_db_session()

		# Save genome with JSON data
		genome = Genome(
			key='foo',
			version='1.0',
			description='test genome',
			extra=JSON_DATA,
		)
		session.add(genome)
		session.commit()

		# Reload in fresh session and check value
		session = empty_db_session()
		genome = session.query(Genome).one()
		assert genome.extra == JSON_DATA

		# Assign different data, save
		genome.extra = JSON_DATA2
		session.commit()

		# Check new value
		session = empty_db_session()
		genome = session.query(Genome).one()
		assert genome.extra == JSON_DATA2

		# Assign NULL, save
		genome.extra = None
		session.commit()

		# Check new value
		session = empty_db_session()
		genome = session.query(Genome).one()
		assert genome.extra is None


class TestReferenceGenomeSet:

	def test_key_lookup(self, testdb_session):
		session = testdb_session()
		key = 'midas/test/testdb_210126'
		gset = session.query(ReferenceGenomeSet).filter_by(key=key).one()
		assert ReferenceGenomeSet.by_key(session, gset.key) == gset
		assert ReferenceGenomeSet.by_key(session, gset.key, '1.0') == gset
		assert ReferenceGenomeSet.by_key(session, gset.key, '1.1') is None

	def test_root_taxa(self, testdb_session):
		session = testdb_session()
		gset = session.query(ReferenceGenomeSet).one()

		root = gset.taxa.filter_by(name='root').one()
		assert gset.root_taxa().all() == [root]

		# This is a read-only session specific to this test, we are free to make modifications.
		new_roots = set(root.children)
		session.delete(root)
		assert set(gset.root_taxa()) == new_roots

	def test_extra_json(self, empty_db_session):
		"""Test storing JSON data in the 'extra' column."""
		session = empty_db_session()

		# Save genome set with JSON data
		gset = ReferenceGenomeSet(
			key='foo',
			version='1.0',
			name='test',
			extra=JSON_DATA,
		)
		session.add(gset)
		session.commit()

		# Reload in fresh session and check value
		session = empty_db_session()
		gset = session.query(ReferenceGenomeSet).one()
		assert gset.extra == JSON_DATA

		# Assign different data, save
		gset.extra = JSON_DATA2
		session.commit()

		# Check new value
		session = empty_db_session()
		gset = session.query(ReferenceGenomeSet).one()
		assert gset.extra == JSON_DATA2

		# Assign NULL, save
		gset.extra = None
		session.commit()

		# Check new value
		session = empty_db_session()
		gset = session.query(ReferenceGenomeSet).one()
		assert gset.extra is None


class TestAnnotatedGenome:

	def test_hybrid_props(self, id_lookup_session):
		session = id_lookup_session

		hybrid_attrs = [
			'key',
			'version',
			'description',
			'ncbi_db',
			'ncbi_id',
			'genbank_acc',
			'refseq_acc',
		]

		for annotated in session.query(AnnotatedGenome):
			for attr in hybrid_attrs:
				assert getattr(annotated, attr) == getattr(annotated.genome, attr)


class TestTaxon:

	def test_tree(self, testdb_session):
		"""Test tree structure."""
		session = testdb_session()

		taxa = session.query(Taxon)
		root = taxa.filter_by(name='root').one()

		for taxon in taxa:
			assert taxon.root() == root
			assert taxon.isleaf() == (len(taxon.children) == 0)

			# Test parent/child relationships match
			for child in taxon.children:
				assert child.parent == taxon

			# Test ancestors() and lineage() methods
			ancestors = list(taxon.ancestors(incself=True))
			assert ancestors[0] is taxon
			assert ancestors[-1] is root
			assert list(taxon.ancestors()) == ancestors[1:]
			assert list(reversed(taxon.lineage())) == ancestors

			for i in range(len(ancestors) - 1):
				assert ancestors[i].parent is ancestors[i + 1]

			# Test descendants() and leaves() methods
			descendants = {t for t in taxa if taxon in t.ancestors(incself=True)}
			assert set(taxon.descendants()) == descendants - {taxon}
			assert set(taxon.descendants(incself=True)) == descendants
			assert set(taxon.leaves()) == {d for d in descendants if d.isleaf()}

	def test_extra_json(self, empty_db_session):
		"""Test storing JSON data in the 'extra' column."""
		session = empty_db_session()

		# Save taxon with JSON data
		gset = ReferenceGenomeSet(
			key='test',
			version='1.0',
			name='test genome set',
		)
		taxon = Taxon(
			genome_set=gset,
			name='test taxon',
			extra=JSON_DATA,
		)
		session.add(gset)
		session.add(taxon)
		session.commit()

		# Reload in fresh session and check value
		session = empty_db_session()
		taxon = session.query(Taxon).one()
		assert taxon.extra == JSON_DATA

		# Assign different data, save
		taxon.extra = JSON_DATA2
		session.commit()

		# Check new value
		session = empty_db_session()
		taxon = session.query(Taxon).one()
		assert taxon.extra == JSON_DATA2

		# Assign NULL, save
		taxon.extra = None
		session.commit()

		# Check new value
		session = empty_db_session()
		taxon = session.query(Taxon).one()
		assert taxon.extra is None


class TestGenomeIDMapping:
	"""Test mapping genomes to ID values."""

	def test__genome_id_attr(self):
		"""Test _check_genome_id_attr() function."""

		for key, attr in models.GENOME_ID_ATTRS.items():
			assert models._check_genome_id_attr(key) is attr
			assert models._check_genome_id_attr(attr) is attr

		for arg in ['description', Genome.description, AnnotatedGenome.key]:
			with pytest.raises(ValueError):
				models._check_genome_id_attr(arg)

	def test__get_genome_id(self, id_lookup_session):
		"""Test _get_genome_id() function."""
		session = id_lookup_session

		for genome in session.query(Genome):
			for key, attr in models.GENOME_ID_ATTRS.items():
				assert models._get_genome_id(genome, attr) == getattr(genome, key)

	def test_genomes_by_id(self, id_lookup_session):
		"""Test genomes_by_id() and genomes_by_id_subset() functions."""
		session = id_lookup_session
		random.seed(0)

		gset = session.query(ReferenceGenomeSet).one()
		genomes = list(gset.genomes)
		random.shuffle(genomes)

		genomes_missing = list(genomes)
		genomes_missing.insert(3, None)
		genomes_missing.insert(8, None)

		for key, attr in models.GENOME_ID_ATTRS.items():
			# Check with all valid IDs
			ids = [models._get_genome_id(g, attr) for g in genomes]
			assert models.genomes_by_id(gset, attr, ids) == genomes
			assert models.genomes_by_id(gset, key, ids) == genomes

			# with missing values
			ids_missing = [
				models._get_genome_id(g, attr) if g is not None else '!INVALID!'
				for g in genomes_missing
			]

			with pytest.raises(KeyError):
				models.genomes_by_id(gset, attr, ids_missing)

			assert models.genomes_by_id(gset, attr, ids_missing, strict=False) == genomes_missing

			genomes_sub, idxs_sub = models.genomes_by_id_subset(gset, attr, ids_missing)
			assert genomes_sub == genomes
			assert [genomes_missing[i] for i in idxs_sub] == genomes
