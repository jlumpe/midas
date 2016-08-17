
import pytest

import os
from io import StringIO
import gzip
import numpy as np

from midas.database.basicdatabase import BasicDatabase



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

		self.db.store_sequence(self.genome, src, **kwargs)

		# Refresh genome as the sequence was created in a different session
		self.session.refresh(self.genome)
		assert self.genome.sequence is not None

		seq_str = self.db.open_sequence(self.genome.sequence).read()
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

	def test_deletion(self, storage_tester):
		"""Test removeal of data file after sequence instance deleted"""
		seq_buf = StringIO(storage_tester.seq_str)
		storage_tester.store(seq_buf)

		db = storage_tester.db
		session = storage_tester.session

		seq_file_name = storage_tester.genome.sequence._filename
		seq_file_path = os.path.join(db.path, db._seq_dir, seq_file_name)

		assert os.path.exists(seq_file_path)

		storage_tester.genome.sequence = None
		session.commit()

		assert not os.path.exists(seq_file_path)

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
			db.store_sequence(genome, seq_buf)

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
			assert db.open_sequence(genome.sequence).read() == seq_str


class TestCoordsStorage:
	"""Tests for k-mer coordinate storage"""

	@pytest.fixture()
	def dbobjs(self, db):
		"""Create database objects"""

		session = db.get_session()

		# Create test genome
		genome = db.Genome(description='test', is_assembled=True)
		session.add(genome)
		session.commit()

		# Create test KmerSetCollection
		kcol = db.KmerSetCollection(name='test', prefix='ATGAC', k=11)
		session.add(kcol)
		session.commit()

		return genome, kcol


int_dtypes = [np.uint8, np.int8, np.uint16, np.int16, np.uint32, np.int32,
              np.uint64, np.int64]

@pytest.fixture(params=int_dtypes)
def coords(request):
	"""
	Creates coordinate arrays for all integral dtypes spanning approximately
	their full positive range
	"""
	dtype = request.param

	max_val = np.iinfo(dtype).max

	return np.arange(0, max_val, max(max_val / 10000, 1), dtype=dtype)


def test_coords_storage(db, coords):
	"""Test that kmer coordindate arrays can be stored and retrieved"""

	session = db.get_session()

	# Create test genome
	genome = db.Genome(description='test', is_assembled=True)
	session.add(genome)
	session.commit()

	# Create test KmerSetCollection
	kcol = db.KmerSetCollection(name='test', prefix='ATGAC',
	                            k=4 * coords.dtype.itemsize)
	session.add(kcol)
	session.commit()

	# Store coords
	db.store_kset_coords(genome, kcol, coords)

	# Get kset instance
	kset = session\
		.query(db.KmerSet)\
		.filter_by(genome=genome, collection=kcol)\
		.one()

	# Load back again
	loaded = db.load_kset_coords(kset)

	# Check equal
	assert np.array_equal(coords, loaded)
