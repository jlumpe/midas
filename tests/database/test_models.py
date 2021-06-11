"""Test midas.db.models.

Uses the included testdb_210126 database.
"""

import pytest

from sqlalchemy.orm import sessionmaker

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
