"""Test midas.database.models."""

import pytest
import numpy as np

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from midas.database import models


@pytest.fixture()
def engine():
	engine = create_engine('sqlite:///:memory:')

	models.Base.metadata.create_all(engine)

	return engine


@pytest.fixture()
def Session(engine):
	return sessionmaker(bind=engine)


def test_annotatedgenome_hybrid_props():
	"""Test hybrid properties on AnnotatedGenome."""

	hybrid_attrs = dict(
		description='Test genome',
		is_assembled=False,
		ncbi_taxid=123,
		entrez_db='assembly',
		entrez_id='123456',
		genbank_acc='GCA_000000.0',
		refseq_acc='GCF_000000.0',
	)

	genome = models.Genome(**hybrid_attrs)

	annotated = models.AnnotatedGenome(genome=genome)

	for attrname in hybrid_attrs:
		assert getattr(genome, attrname) == getattr(annotated, attrname)


class TestSignatureStorage:

	# Number of genomes to test
	N = 100

	@pytest.fixture()
	def session(self, Session):
		return Session()

	@pytest.fixture()
	def sigset(self, session):

		sigset = models.SignatureSet(
			name='test set',
			k=11,
			prefix='ATGAC',
		)

		session.add(sigset)
		session.commit()

		return sigset

	@pytest.fixture()
	def genomes(self, session):
		genomes = []

		for i in range(self.N):
			genome = models.Genome(description='test {}'.format(i + 1))

			genomes.append(genome)
			session.add(genome)

		session.commit()

		return genomes

	@pytest.fixture()
	def signatures(self, sigset, genomes):

		random = np.random.RandomState(0)

		# Use a variety of numpy integer dtypes
		dtypes = ['u2', 'u4', 'u8', 'i2', 'i4', 'i8']

		signatures = []

		for i, genome in enumerate(genomes):

			signature = np.arange(
				random.randint(5000, 15000),
				dtype=dtypes[i % len(dtypes)]
			)
			signature += i

			signatures.append(signature)

			sigset.store_signature(genome.id, signature)

		return signatures

	def test_get_one(self, sigset, genomes, signatures):
		"""Test getting one signature at a time through get_signature()."""

		for genome, signature in zip(genomes, signatures):
			assert np.array_equal(sigset.get_signature(genome.id), signature)

	def test_get_missing(self, sigset, genomes):
		"""Test trying to get a missing signature."""
		top_id = max(genome.id for genome in genomes)
		assert sigset.get_signature(top_id + 1) is None
