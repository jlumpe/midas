"""
Run a full set of queries using the testdb_210126 database.

Tests in this file will only be run when the --midas-test-full-db option is passed to the pycharm
command.

Database files and query sequences are located in tests/data/testdb_210126, but only the genome
database file is included in version control. Other files need to be obtained separately. See the
Readme.md file in that directory for more information.
"""

from csv import DictReader

import pytest

from midas.io.seq import SequenceFile, find_kmers_in_files
from midas.signatures.hdf5 import HDF5Signatures
from midas.db.midasdb import MIDASDatabase
from midas.db.models import ReferenceGenomeSet
from midas.query import runquery


@pytest.fixture(autouse=True, scope='module')
def testdb_files(request, testdb_dir):
	"""Paths to testdb_210126 files.

	Skips all dependent tests if the --midas-test-full-db command line option is not passed.

	Checks that the directory and required files/subdirectories exist and fails tests immediately
	if they do not.
	"""
	if not request.config.getoption('midas_test_full_db'):
		pytest.skip('--midas-test-full-db option not given')

	files = dict(
		root=testdb_dir,
		db=testdb_dir / 'testdb_210126-genomes.db',
		signatures=testdb_dir / 'testdb_210126-signatures.h5',
		queries=testdb_dir / 'query-seqs/queries.csv',
	)

	for k, v in files.items():
		# Raises FileNotFoundError if does not exist
		files[k] = v.resolve()

	return files


@pytest.fixture(scope='module')
def signatures(testdb_files):
	"""K-mer signatures for test genomes."""
	return HDF5Signatures.open(testdb_files['signatures'])


@pytest.fixture(scope='module')
def testdb(testdb_session, signatures):
	"""Full MIDASDatabase object for test db."""

	session = testdb_session()
	gset = session.query(ReferenceGenomeSet).one()

	return MIDASDatabase(gset, signatures)


@pytest.fixture(scope='module')
def query_data(testdb_files):
	"""Query files and their expected taxa."""
	table_path = testdb_files['queries']
	seqs_dir = table_path.parent

	files = []
	expected_taxa = []

	with open(table_path, newline='') as f:
		for row in DictReader(f):
			file = SequenceFile(
				path=seqs_dir / (row['name'] + '.fa'),
				format='fasta',
			)
			files.append(file)
			expected_taxa.append(row['expected_taxon'])

	return files, expected_taxa


def test_query_python(testdb, query_data):
	"""Run a full query using the Python API."""

	query_files, expected_taxa = query_data
	query_sigs = find_kmers_in_files(testdb.kmerspec, query_files)

	results = runquery(testdb, query_sigs, query_files)

	assert len(results.items) == len(query_files)

	for item, file, expected in zip(results.items, query_files, expected_taxa):
		assert item.input.file == file
		assert item.success

		if expected == '':
			assert item.predicted_taxon is None
			assert item.report_taxon is None
		else:
			assert item.predicted_taxon is not None
			assert item.predicted_taxon.name == expected
			assert item.report_taxon == item.predicted_taxon

		assert not item.warnings
		assert item.error is None


def test_query_cli(testdb_files, query_data):
	"""Run a full query using the command line interface."""
	# TODO
