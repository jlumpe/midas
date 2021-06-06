"""Test midas.db.models.

Uses the included testdb_210126 database.
"""

import pytest

from midas.db.models import Genome, ReferenceGenomeSet, AnnotatedGenome, Taxon


@pytest.fixture()
def session(testdb_session):
	return testdb_session()


class TestGenome:

	def test_key_lookup(self, session):
		key = 'midas/testdb_210126/A1_B1_C1_1'
		genome = session.query(Genome).filter_by(key=key).one()
		assert Genome.by_key(session, genome.key) == genome
		assert Genome.by_key(session, genome.key, '1.0') == genome
		assert Genome.by_key(session, genome.key, '1.1') is None

	def test_ncbi_id_lookup(self, session):
		pass  # TODO


class TestReferenceGenomeSet:

	def test_key_lookup(self, session):
		key = 'midas/test/testdb_210126'
		gset = session.query(ReferenceGenomeSet).filter_by(key=key).one()
		assert ReferenceGenomeSet.by_key(session, gset.key) == gset
		assert ReferenceGenomeSet.by_key(session, gset.key, '1.0') == gset
		assert ReferenceGenomeSet.by_key(session, gset.key, '1.1') is None


class TestAnnotatedGenome:

	def test_hybrid_props(self, session):
		hybrid_attrs = dict(
			key='test_genome',
			version='1.0',
			description='Test genome',
			entrez_db='assembly',
			entrez_id='123456',
			genbank_acc='GCA_000000.0',
			refseq_acc='GCF_000000.0',
		)

		genome = Genome(**hybrid_attrs)
		annotated = AnnotatedGenome(genome=genome)

		for attr, value in hybrid_attrs.items():
			assert getattr(annotated, attr) == value


class TestTaxon:

	def test_tree(self, session):
		"""Test tree structure."""

		root = session.query(Taxon).filter_by(name='root').one()

		for taxon in session.query(Taxon):
			assert taxon.root() == root

			# Test parent/child relationships
			for child in taxon.children:
				assert child.parent == taxon

			# Test lineage() method
			lineage = taxon.lineage()
			assert lineage[0] == root
			for i in range(len(lineage) - 1):
				assert lineage[i] == lineage[i + 1].parent
