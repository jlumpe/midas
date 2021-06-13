"""Test midas.io.seq."""

from io import StringIO
from pathlib import Path

import pytest
import numpy as np
from Bio import Seq, SeqIO

from midas.io.seq import SeqFileInfo, find_kmers_parse, FileSignatureCalculator
import midas.io.util as ioutil
from midas.kmers import dense_to_sparse
from midas.test import make_kmer_seq, random_seq


def create_sequence_records(nseqs, k, prefix_len, seq_len=10000):
	"""Create a set of random DNA sequences with known combined k-mer signature.

	Parameters
	----------
	nseqs : int
		Number of sequences to create.
	k : int
		Length of k-mers (without prefix)
	prefix_len : Int
		Length of prefix.
	seq_len : int
		Length of sequences to create.

	Returns
	-------
	tuple
		(kmerspec, records, kmer_vec) tuple.
	"""
	vec = np.zeros(4 ** k, dtype=bool)
	kspec = None

	records = []

	for i in range(nseqs):

		seq, _kspec, rec_vec = make_kmer_seq(
			seq_len,
			k=k,
			prefix_len=prefix_len,
			kmer_interval=50,
			n_interval=10,
			seed=i,
		)

		# Combine vectors of all sequences
		vec |= rec_vec

		# Just check that our specs are all the same...
		if kspec is None:
			kspec = _kspec
		else:
			assert kspec == _kspec

		# Convert every other sequence to lower case, just to switch things up...
		if i % 2:
			seq = seq.lower()

		# Create the BioPython sequence record object
		records.append(SeqIO.SeqRecord(
			seq=Seq.Seq(seq.decode('ascii')),
			id='SEQ{}'.format(i + 1),
			description='sequence {}'.format(i + 1),
		))

	return kspec, records, vec


@pytest.mark.parametrize('sparse', [False, True])
def test_find_kmers_parse(sparse):
	"""Test the find_kmers_parse function."""

	kspec, records, vec = create_sequence_records(10, k=11, prefix_len=5)

	# Write records to string buffer in FASTA format
	buf = StringIO()
	SeqIO.write(records, buf, 'fasta')
	buf.seek(0)

	# Parse from buffer
	kmers = find_kmers_parse(kspec, buf, 'fasta', sparse=sparse)

	if sparse:
		assert np.array_equal(kmers, dense_to_sparse(vec))
	else:
		assert np.array_equal(kmers, vec)


class TestSeqFileInfo:
	"""Test the SeqFileInfo class."""

	@pytest.fixture(params=['fasta'], scope='class')
	def fmt(self, request):
		"""SeqFileInfo.fmt attribute."""
		return request.param

	@pytest.fixture(params=list(ioutil.COMPRESSED_OPENERS), scope='class')
	def compression(self, request):
		"""SeqFileInfo.compression attribute."""
		return request.param

	@pytest.fixture()
	def info(self, tmpdir, fmt, compression):
		"""A SeqFileInfo instance pointing to a file in a test temporary directory.

		File does not yet exist.
		"""
		path = tmpdir.join('test.' + fmt).strpath
		return SeqFileInfo(path, fmt, compression)

	@pytest.fixture(scope='class')
	def seqrecords(self):
		"""A collection of random Bio.SeqIO.SeqRecord's."""
		records = []

		for i in range(20):
			seq = Seq.Seq(random_seq(1000))
			id_ = 'seq{}'.format(i + 1)
			descr = 'Test sequence {}'.format(i + 1)
			records.append(SeqIO.SeqRecord(seq, id=id_, description=descr))

		return tuple(records)

	@pytest.fixture
	def file_contents(self, fmt, seqrecords):
		"""String contents of a file containing the sequence records."""
		buf = StringIO()
		SeqIO.write(seqrecords, buf, fmt)
		return buf.getvalue()

	@pytest.fixture
	def info_exists(self, info, seqrecords):
		"""Copy of "info" fixture, but with "seqrecords" written to the file."""

		with info.open('rt') as fobj:
			SeqIO.write(seqrecords, fobj, info.fmt)

	def test_constructor(self):
		"""Test constructor."""

		info = SeqFileInfo('foo.fasta', 'fasta')
		assert info == SeqFileInfo('foo.fasta', 'fasta', None)
		assert info.path == Path('foo.fasta')

		# Argument fails validation
		with pytest.raises(TypeError):
			SeqFileInfo(b'foo.fasta', 'fasta')
		with pytest.raises(TypeError):
			SeqFileInfo('foo.fasta', b'fasta')
		with pytest.raises(TypeError):
			SeqFileInfo('foo.fasta', 'fasta', b'gzip')

	def test_eq(self):
		"""Test equality checking of instances."""
		infos = [
			SeqFileInfo(p, fmt, comp)
			for p in ['foo', 'bar']
			for fmt in ['fasta', 'genbank']
			for comp in [None, 'gzip']
		]

		for i, info1 in enumerate(infos):
			for j, info2 in enumerate(infos):
				if i == j:
					# Try with different instance
					assert info1 == SeqFileInfo(info1.path, info1.fmt,
												info1.compression)
				else:
					assert info1 != info2

	@pytest.mark.parametrize('binary', [False, True])
	def test_open(self, info, file_contents, binary):
		"""Test sequence file is readable and writable."""

		to_write = file_contents.encode() if binary else file_contents

		# Write data to file
		with info.open('wb' if binary else 'wt') as fobj:
			fobj.write(to_write)

		# Read it back and make sure it's the same
		with info.open('rb' if binary else 'rt') as fobj:
			read = fobj.read()

		assert read == to_write

	def test_parse(self, info, seqrecords, file_contents):
		"""Test the parse() method, ensure we get the right records back."""

		# Write pre-formatted contents to file
		with info.open('w') as fobj:
			fobj.write(file_contents)

		# Parse the sequences from it
		parsed = list(info.parse())

		# Check they match
		assert len(parsed) == len(seqrecords)

		for parsed_req, orig_req in zip(parsed, seqrecords):
			assert isinstance(parsed_req, SeqIO.SeqRecord)
			assert parsed_req.seq == orig_req.seq
			assert parsed_req.id == orig_req.id

			# This is something stupid BioPython does - when writing a SeqRecord
			# as FASTA it writes the .id attributed followed by a space and then
			# the .description attribute on the description line. When reading,
			# the entire line is used as the description attribute and so
			# includes the ID
			assert parsed_req.description == orig_req.id + ' ' + orig_req.description

	def test_path_arg(self):
		"""Test the "path" argument to the constructor."""

		path = Path('foo/bar.fasta')

		info1 = SeqFileInfo(path, 'fasta')
		assert isinstance(info1, SeqFileInfo) and info1.path == path

		info2 = SeqFileInfo(str(path), 'fasta')
		assert isinstance(info2, SeqFileInfo) and info2.path == path

	def test_absolute(self):
		"""Test the absolute() method."""

		relinfo = SeqFileInfo('foo/bar.fasta', 'fasta')
		assert not relinfo.path.is_absolute()

		absinfo = relinfo.absolute()
		assert absinfo.path.is_absolute()
		assert absinfo.path == relinfo.path.absolute()

		absinfo2 = absinfo.absolute()
		assert absinfo2 == absinfo

	def test_from_paths(self, fmt, compression):
		"""Test the from_paths() class method."""

		# List of unique path strings
		paths = ['foo/bar{}.{}'.format(i, fmt) for i in range(20)]

		infos = SeqFileInfo.from_paths(paths, fmt, compression)

		assert len(paths) == len(infos)

		for path, info in zip(paths, infos):
			assert isinstance(info, SeqFileInfo)
			assert str(info.path) == path
			assert info.fmt == fmt
			assert info.compression == compression


@pytest.mark.parametrize('fmt', ['fasta'])
@pytest.mark.parametrize('compression', list(ioutil.COMPRESSED_OPENERS))
@pytest.mark.parametrize('use_threads', [False, True])
@pytest.mark.parametrize('pass_paths', [False, True])
@pytest.mark.parametrize('ordered', [False, True])
def test_FileSignatureCalculator(tmpdir, fmt, compression, use_threads,
                                 pass_paths, ordered):
	"""Test the FileSignatureCalculator class."""

	k = 11
	prefix_len = 5

	files = []
	signatures = []

	kspec = None

	# Create some random sequence records with their signatures and write them
	# to files
	for i in range(15):

		_kspec, records, sig_vec = create_sequence_records(
			nseqs=20,
			k=k,
			prefix_len=prefix_len,
			seq_len=1000,
		)

		# create_sequence_records() creates its own KmerSpec, make sure that
		# they're all the same
		if kspec is None:
			kspec = _kspec
		else:
			assert _kspec == kspec

		# SeqFileInfo for file to eventually read
		path = tmpdir.join('seq{}.{}'.format(i, fmt)).strpath
		file = SeqFileInfo(path, fmt, compression)

		# Write sequences to file
		with file.open('wt') as fobj:
			SeqIO.write(records, fobj, fmt)

		files.append(file)
		signatures.append(dense_to_sparse(sig_vec))

	# Create the calculator
	calculator = FileSignatureCalculator(use_threads=use_threads)

	# Run jobs in context
	with calculator as enterval:

		# Check __enter__ returns self
		assert enterval is calculator

		# Start calculating
		if pass_paths:
			results = calculator.calc_signatures(
				kspec,
				[f.path for f in files],
				fmt=fmt,
				compression=compression,
				ordered=ordered,
			)
		else:
			results = calculator.calc_signatures(
				kspec,
				files,
				ordered=ordered,
			)

		last_i = -1

		# Iterate over results
		for i, sig in results:

			# Check in order
			if ordered:
				assert i == last_i + 1
				last_i = i

			# Check signature matches
			assert np.array_equal(sig, signatures[i])
