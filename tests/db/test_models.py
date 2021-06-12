"""Test midas.db.models.

Uses the included testdb_210126 database.
"""

import random

import pytest
from sqlalchemy.orm import sessionmaker

from midas.db import models
from midas.db.models import Genome, ReferenceGenomeSet, AnnotatedGenome, Taxon


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
		reference_set=gset,
	)
	session.add(roottaxon)

	for i in range(20):
		g = Genome(
			key=f'test/genome_{i}',
			version='1.0',
			description=f'Test genome {i}',
			entrez_db='assembly',
			entrez_id=i,
			genbank_acc=f'GCA_{i:09d}.1',
			refseq_acc=f'GCF_{i:09d}.1',
		)
		ag = AnnotatedGenome(
			reference_set=gset,
			genome=g,
			primary_taxon=roottaxon,
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


class TestAnnotatedGenome:

	def test_hybrid_props(self, id_lookup_session):
		session = id_lookup_session

		hybrid_attrs = [
			'key',
			'version',
			'description',
			'entrez_db',
			'entrez_id',
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
