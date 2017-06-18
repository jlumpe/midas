"""Test midas.filestore.FileSequenceStore."""

import pytest

import os
from io import StringIO
import gzip

from midas.seqstore.base import SequenceNotFound, SequenceIDConflict
from midas.seqstore.filestore import FileSequenceStore
from midas import ncbi


TEST_IDS_1 = dict(
	genbank_acc='XX123456.1',
	refseq_acc='XX_123456.1',
	entrez_db='nuccore',
	entrez_id=1,
)

TEST_IDS_2 = dict(
	genbank_acc='YY123456.1',
	refseq_acc='YY_123456.1',
	entrez_db='nuccore',
	entrez_id=2,
)



@pytest.fixture()
def filestore(tmpdir):
	"""Creates a new database in a temporary directory."""
	return FileSequenceStore.create(tmpdir.join('seqstore').strpath)


@pytest.fixture(scope='module')
def seq_str():
	return '>Some FASTA file\n' + ('ATGCTAGCTAGC' * 100)


@pytest.fixture()
def uncomp_seq_file(tmpdir, seq_str):
	"""Path for uncompressed text file containing sequence."""
	seq_path = tmpdir.join('test_seq.fasta')
	with seq_path.open('w') as fh:
		fh.write(seq_str)

	return seq_path


@pytest.fixture()
def gzip_seq_file(tmpdir, seq_str):
	"""Path for gzip-compressed text file containing sequence."""
	seq_path = tmpdir.join('test_seq.fasta')
	with gzip.open(seq_path.strpath, 'wt') as fh:
		fh.write(seq_str)

	return seq_path


@pytest.fixture()
def store_test_func(filestore, seq_str):
	"""Function that stores sequence data and verifies it is readable."""

	def _test_store(src, **kwargs):

		ids = dict(genbank_acc=TEST_IDS_1['genbank_acc'])

		record = filestore.store(src, ids, **kwargs)

		assert filestore.has(**ids)

		read_str = filestore.open(record).read()
		assert read_str == seq_str

	return _test_store


def test_text_buffer(seq_str, store_test_func):
	"""Add from StringIO buffer."""
	seq_buf = StringIO(seq_str)
	store_test_func(seq_buf)


def test_uncomp_fname(uncomp_seq_file, store_test_func):
	"""Add from name of uncompressed text file."""
	store_test_func(uncomp_seq_file.strpath)

	# Check file still exists
	assert uncomp_seq_file.exists()


def test_uncomp_fname_nokeep(uncomp_seq_file, store_test_func):
	"""Add from name of uncompressed text file, don't keep source."""
	store_test_func(uncomp_seq_file.strpath, keep_src=False)

	# Check file no longer exists at location
	assert not uncomp_seq_file.exists()


def test_gzip_fname(gzip_seq_file, store_test_func):
	"""Add from name of gzip file."""
	store_test_func(gzip_seq_file.strpath, src_compression='gzip')

	# Check file still exists
	assert gzip_seq_file.exists()


def test_gzip_fname_nokeep(gzip_seq_file, store_test_func):
	"""Add from name of gzip file, don't keep source."""
	store_test_func(gzip_seq_file.strpath, src_compression='gzip',
	                     keep_src=False)

	# Check file no longer exists at location
	assert not gzip_seq_file.exists()


def test_uncomp_fileobj_text(uncomp_seq_file, store_test_func):
	"""Add from uncompressed text file object, text mode."""

	# Text mode
	with open(uncomp_seq_file.strpath) as fh:
		store_test_func(fh)


def test_uncomp_fileobj_binary(uncomp_seq_file, store_test_func):
	"""Add from uncompressed text file object, binary mode."""

	# Binary mode
	with open(uncomp_seq_file.strpath, 'rb') as fh:
		store_test_func(fh, src_mode='b')


def test_gzip_fileobj(gzip_seq_file, store_test_func):
	"""Add gzip file object."""
	with open(gzip_seq_file.strpath, 'rb') as fh:
		store_test_func(fh, src_compression='gzip')


def test_get_version(filestore):
	"""Test getting version of existing file store."""

	version = FileSequenceStore.get_version(filestore.path)
	assert version == FileSequenceStore.VERSION


def test_store_duplicate_id(filestore, uncomp_seq_file):
	"""Test storing a sequence with the same ID as an existing one."""

	ids1 = TEST_IDS_1
	filestore.store(uncomp_seq_file.strpath, ids1)


	# One full duplicated ID, others unique
	ids2 = dict(TEST_IDS_2)
	ids2['genbank_acc'] = ids1['genbank_acc']

	with pytest.raises(SequenceIDConflict):
		filestore.store(uncomp_seq_file.strpath, ids2)

	# Duplicate for only one attribute of two in index, should be fine
	ids3 = dict(TEST_IDS_2)
	ids3['entrez_db'] = ids1['entrez_db']
	ids3['entrez_id'] = ids1['entrez_id'] + 1

	filestore.store(uncomp_seq_file.strpath, ids3)


def test_get_record(filestore, uncomp_seq_file):
	"""Test getting the record for a stored sequence."""

	ids = dict(genbank_acc=TEST_IDS_1['genbank_acc'])

	# Add the sequence
	record = filestore.store(uncomp_seq_file.strpath, ids)

	# Test getting by ncbi ids
	assert filestore.get_record(**ids) == record

	# Test getting by store id
	assert filestore.get_record(record.store_id) == record

	# Check nonexistant IDs yield None
	fake_id = -999999
	assert fake_id != record.store_id
	assert filestore.get_record(fake_id) is None

	assert filestore.get_record(genbank_acc=TEST_IDS_2['genbank_acc']) is None


def test_has_ids(filestore, uncomp_seq_file):
	"""Test checking if the store has sequences with the given ids."""

	ids1 = TEST_IDS_1
	filestore.store(uncomp_seq_file.strpath, ids1)

	# Test has and has_any with correct IDs
	assert filestore.has(**ids1)
	assert filestore.has_any(**ids1)
	for index in ncbi.SEQ_ID_INDICES:
		index_ids = {key: ids1[key] for key in index}
		assert filestore.has(**index_ids)
		assert filestore.has_any(**index_ids)

	# Test with some different - should be false for has() and true for has_any()
	for index in ncbi.SEQ_ID_INDICES:

		ids2 = dict(ids1)
		ids2.update({key: TEST_IDS_2[key] for key in index})

		assert not filestore.has(**ids2)
		assert filestore.has_any(**ids2)


@pytest.mark.parametrize('use_id', [False, True])
def test_remove(filestore, uncomp_seq_file, use_id):
	"""Test removing a stored sequence."""

	ids = dict(genbank_acc=TEST_IDS_1['genbank_acc'])

	# Add the sequence
	record = filestore.store(uncomp_seq_file.strpath, ids)
	assert filestore.has(**ids)

	# Get the path of the stored sequence file (should be only file in directory)
	file_path, = filestore.seq_dir.iterdir()
	assert file_path.exists()

	# Remove it and check it no longer exists in the database or directory
	filestore.remove(record.store_id if use_id else record)

	assert not filestore.has(**ids)
	assert not file_path.exists()

	# Removing again should raise an exception as it does not exist anymore
	with pytest.raises(SequenceNotFound):
		filestore.remove(record)


def test_store_error(filestore, uncomp_seq_file):
	"""
	Test proper cleanup is performed after an error when adding a sequence.
	"""

	ids = dict(genbank_acc=TEST_IDS_1['genbank_acc'])

	# Get a closed file object
	with open(uncomp_seq_file.strpath) as fobj:
		pass

	# Try adding - IO/operation on closed file should raise ValueError
	with pytest.raises(ValueError):
		filestore.store(fobj, ids)

	# Check the record has been removed
	assert not filestore.has(**ids)
