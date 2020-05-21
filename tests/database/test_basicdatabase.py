
import pytest

import os
from io import StringIO
import gzip
import numpy as np

from midas.database.basicdatabase import BasicDatabase
from midas.kmers import KmerCoordsCollection



@pytest.fixture()
def db(tmpdir):
	"""Creates a new database in a temporary directory"""
	return BasicDatabase.create(tmpdir.join('db').strpath)


class SeqStorageTester:
	"""Helper class to test sequence storage"""

	def __init__(self, db, seq_str):

		self.db = db
		self.session = db.get_session()

		self.seq_str = seq_str

		# Create test genome
		self.genome = self.db.Genome(description='test', is_assembled=True)
		self.session.add(self.genome)
		self.session.commit()

	def store(self, src, **kwargs):
		"""Store sequence data and verify it is readable"""

		self.db.store_sequence(self.genome.id, src, **kwargs)

		# Refresh genome as the sequence was created in a different session
		self.session.refresh(self.genome)
		assert self.genome.sequence is not None

		seq_str = self.db.open_sequence(self.genome.id).read()
		assert seq_str == self.seq_str


class TestSeqStorage:
	"""Tests sequence storage in database"""

	@pytest.fixture()
	def seq_str(self):
		return 'ATGCTAGCTAGC' * 100

	@pytest.fixture()
	def storage_tester(self, db, seq_str):
		return SeqStorageTester(db, seq_str=seq_str)

	@pytest.fixture()
	def uncomp_seq_file(seq, tmpdir, seq_str):
		"""Path for uncompressed text file containing sequence"""
		seq_path = tmpdir.join('test_seq.fasta')
		with seq_path.open('w') as fh:
			fh.write(seq_str)

		return seq_path

	@pytest.fixture()
	def gzip_seq_file(seq, tmpdir, seq_str):
		"""Path for gzip-compressed text file containing sequence"""
		seq_path = tmpdir.join('test_seq.fasta')
		with gzip.open(seq_path.strpath, 'wt') as fh:
			fh.write(seq_str)

		return seq_path

	def test_text_buffer(self, storage_tester):
		"""Add from StringIO buffer"""
		seq_buf = StringIO(storage_tester.seq_str)
		storage_tester.store(seq_buf)

	def test_uncomp_fname(self, uncomp_seq_file, storage_tester):
		"""Add from name of uncompressed text file"""
		storage_tester.store(uncomp_seq_file.strpath)

		# Check file still exists
		assert uncomp_seq_file.exists()

	def test_uncomp_fname_nokeep(self, uncomp_seq_file, storage_tester):
		"""Add from name of uncompressed text file, don't keep source"""
		storage_tester.store(uncomp_seq_file.strpath, keep_src=False)

		# Check file no longer exists at location
		assert not uncomp_seq_file.exists()

	def test_gzip_fname(self, gzip_seq_file, storage_tester):
		"""Add from name of gzip file"""
		storage_tester.store(gzip_seq_file.strpath, src_compression='gzip')

		# Check file still exists
		assert gzip_seq_file.exists()

	def test_gzip_fname_nokeep(self, gzip_seq_file, storage_tester):
		"""Add from name of gzip file, don't keep source"""
		storage_tester.store(gzip_seq_file.strpath, src_compression='gzip',
		                     keep_src=False)

		# Check file no longer exists at location
		assert not gzip_seq_file.exists()

	def test_uncomp_fileobj_text(self, uncomp_seq_file, storage_tester):
		"""Add from uncompressed text file object, text mode"""

		# Text mode
		with open(uncomp_seq_file.strpath) as fh:
			storage_tester.store(fh)

	def test_uncomp_fileobj_binary(self, uncomp_seq_file, storage_tester):
		"""Add from uncompressed text file object, binary mode"""

		# Binary mode
		with open(uncomp_seq_file.strpath, 'rb') as fh:
			storage_tester.store(fh, src_mode='b')

	def test_gzip_fileobj(self, gzip_seq_file, storage_tester):
		"""Add gzip file object"""
		with open(gzip_seq_file.strpath, 'rb') as fh:
			storage_tester.store(fh, src_compression='gzip')

	def test_clean(self, db, seq_str):
		"""Test cleaning of orphaned sequence files"""

		session = db.get_session()

		# Add some genomes
		genomes = []
		for i in range(10):
			genome = db.Genome(description='test{}'.format(i+1),
			                   is_assembled=False)
			genomes.append(genome)
			session.add(genome)
		session.commit()

		# Store real sequences
		seq_buf = StringIO(seq_str)
		for genome in genomes:
			seq_buf.seek(0)
			db.store_sequence(genome.id, seq_buf)

		# Add some orphaned files
		orphaned_files = []
		for i in range(10):
			fname = 'orphan{}.seq'.format(i+1)
			orphaned_files.append(fname)
			with open(db._get_seq_file_path(fname), 'w') as fh:
				fh.write(seq_str)

		# Clean
		cleaned_files = set(db.clean_seq_files())

		# Check orphaned files are gone and in returned list
		for fname in orphaned_files:
			assert fname in cleaned_files
			assert not os.path.exists(db._get_seq_file_path(fname))

		# Check all real files are still readable
		for genome in genomes:
			assert db.open_sequence(genome.id).read() == seq_str


class TestCoordsStorage:
	"""Tests for k-mer coordinate storage."""

	# Number of genomes to store coords for
	N = 10

	@pytest.fixture(params=['i1', 'u1', 'i2', 'u2', 'i4', 'u4', 'i8', 'u8'])
	def coords_dtype(self, request):
		return np.dtype(request.param)

	@pytest.fixture()
	def coords_db(self, db, coords_dtype):
		"""Database with test genomes and (empty) KmerSetCollection."""

		session = db.get_session()

		k = 4 * coords_dtype.itemsize
		kcol = db.KmerSetCollection(name='test', prefix='ATGAC', k=k)
		session.add(kcol)

		for i in range(self.N):
			session.add(db.Genome(description='test%d' % i))

		session.commit()

		return db

	@pytest.fixture
	def coords_col(self, coords_dtype):
		"""
		Several coordinate arrays spanning approximately the full positive range
		of the given dtype.
		"""
		max_val = np.iinfo(coords_dtype).max - self.N

		# At least n apart, if possible enough for 10000 elements up to max
		step = max(self.N, max_val / 10000)
		coords0 = np.arange(0, max_val, step)
		coords_all = [coords0 + i for i in range(self.N)]

		return KmerCoordsCollection.from_coords_seq(coords_all, dtype=coords_dtype)

	def test_coords_storage(self, coords_db, coords_col):
		"""Test that kmer coordindate arrays can be stored and retrieved"""

		db = coords_db
		session = db.get_session()

		kcol = session.query(db.KmerSetCollection).one()
		genomes = session.query(db.Genome).all()

		# Store coords
		for genome, coords in zip(genomes, coords_col):
			db.store_kset_coords(kcol.id, genome.id, coords)

		# Load back again
		for genome, coords in zip(genomes, coords_col):
			loaded = db.load_kset_coords(kcol.id, genome.id)
			assert np.array_equal(coords, loaded)
