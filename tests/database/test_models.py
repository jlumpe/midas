"""Test midas.db.models."""

import pytest
import numpy as np

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from midas.db import models


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
		entrez_db='assembly',
		entrez_id='123456',
		genbank_acc='GCA_000000.0',
		refseq_acc='GCF_000000.0',
	)

	genome = models.Genome(**hybrid_attrs)

	annotated = models.AnnotatedGenome(genome=genome)

	for attrname in hybrid_attrs:
		assert getattr(genome, attrname) == getattr(annotated, attrname)
