"""Tests for midas.kmers module."""

import itertools

import pytest
import numpy as np

from midas import kmers
from midas.cython.seqs import reverse_complement


# Complements to nucleotide ASCII codes
NUC_COMPLEMENTS = {
	65: 84,
	84: 65,
	71: 67,
	67: 71,
	97: 116,
	116: 97,
	103: 99,
	99: 103,
}


def make_kmerspec(k):
	"""Create a KmerSpec with some arbitrary prefix."""
	return kmers.KmerSpec(k, 'ATGC')


def test_nucleotide_order():
	"""Check k-mer indices correspond to defined nucleotide order."""

	for i, nuc in enumerate(kmers.NUCLEOTIDES):
		assert kmers.kmer_to_index(bytes([nuc])) == i


def test_index_conversion():
	"""Test converting k-mers to and from their indices."""

	# Test for k in range 0-10
	for k in range(11):

		# Test all indices to max of 1000
		for index in range(min(4 ** k, 1000)):

			kmer = kmers.index_to_kmer(index, k)

			# Check k-mer is of correct length
			assert len(kmer) == k

			# Check converting back, both cases
			assert kmers.kmer_to_index(kmer.upper()) == index
			assert kmers.kmer_to_index(kmer.lower()) == index


class TestKmerSpec:
	"""Test midas.kmers.KmerSpec."""

	def test_attributes(self):
		"""Test basic attributes."""

		# Try k from 0 to 32 (all have dtypes)
		for k in range(33):

			spec = make_kmerspec(k)

			# Check length attributes
			assert spec.prefix_len == len(spec.prefix)
			assert spec.prefix_len + spec.k == spec.total_len

		# Check prefix is bytes
		assert isinstance(kmers.KmerSpec(11, 'ATGAC').prefix, bytes)

	def test_dtype(self):
		"""Test coords_dtype attribute"""

		# Try k from 0 to 32 (all have dtypes)
		for k in range(33):

			spec = make_kmerspec(k)

			# Check dtype can store the largest index
			top_idx = spec.idx_len - 1
			assert spec.coords_dtype.type(top_idx) == top_idx

		# k > 32 should have no dtype
		assert make_kmerspec(33).coords_dtype is None

	def test_eq(self):
		"""Test equality testing."""

		kspec = kmers.KmerSpec(11, 'ATGAC')

		assert kspec == kmers.KmerSpec(11, 'ATGAC')
		assert kspec == kmers.KmerSpec(11, b'ATGAC')

		assert kspec != kmers.KmerSpec(11, 'ATGAA')
		assert kspec != kmers.KmerSpec(12, 'ATGAC')


def test_vec_coords_conversion():
	"""Test conversion between k-mer vector and coordinates."""

	for k in range(10):

		spec = make_kmerspec(k)

		# Create vector with every 3rd k-mer
		vec = np.zeros(spec.idx_len, dtype=bool)
		vec[np.arange(vec.size) % 3 == 0] = True

		# Convert to coords
		coords = kmers.vec_to_coords(vec)

		# Check coords
		assert len(coords) == vec.sum()
		for index in coords:
			assert vec[index]

		# Check coords ascending
		assert np.all(np.diff(coords) > 0)

		# Check converting back
		assert np.array_equal(vec, spec.coords_to_vec(coords))


def check_reverse_complement(seq, rc):
	"""Assert the reverse complement of a seqeuence is correct.

	:param bytes seq: Byte sequence
	:param bytes rc: Reverse complement of seq
	"""
	l = len(seq)
	for i in range(l):
		assert rc[l - i - 1] == NUC_COMPLEMENTS.get(seq[i], seq[i])


def test_revcomp():
	"""Test midas.cython.seqs.reverse_complement."""

	# Check empty
	assert reverse_complement(b'') == b''

	# Check one-nucleotide values
	for nuc1, nuc2 in NUC_COMPLEMENTS.items():
		b1, b2 = [bytes([n]) for n in [nuc1, nuc2]]
		assert reverse_complement(b1) == b2
		assert reverse_complement(b1.lower()) == b2.lower()

	# Check single invalid code
	assert reverse_complement(b'N') == b'N'
	assert reverse_complement(b'n') == b'n'

	# Check all 6-mers
	k = 6
	for i in range(4 ** k):
		kmer = kmers.index_to_kmer(i, k)

		rc = reverse_complement(kmer)

		check_reverse_complement(rc, kmer)
		check_reverse_complement(rc.lower(), kmer.lower())

		assert reverse_complement(rc) == kmer
		assert reverse_complement(rc.lower()) == kmer.lower()

	# Check longer seqs with invalid nucleotides
	seq = bytearray(b'ATGCatgc')

	for i in range(len(seq)):

		array = bytearray(seq)
		array[i] = ord(b'N')
		seq2 = bytes(array)

		rc = reverse_complement(seq2)

		check_reverse_complement(rc, seq2)
		assert reverse_complement(rc) == seq2


def fill_bytearray(pattern, length):
	"""Create a bytearray with a repeating pattern.

	:param bytes pattern: Pattern to repeat in array.
	:param int length: Length of array to create.
	:returns: Filled array
	:rtype: bytearray
	"""

	array = bytearray(length)

	n = len(pattern)

	for i in range(0, length, n):
		n2 = min(n, length - i)
		array[i:i + n2] = pattern[:n2]

	return array


def make_kmer_seq(seqlen, k, prefix_len, kmer_interval, n_interval=None, seed=0):
	"""Create a test sequence with a known set of k-mers present.

	:param int seqlen: Length of sequence.
	:param int k: Length of k-mers to find (not including prefix).
	:param int prefix_len: Prefix preceding k-mers to match.
	:param int kmer_interval: Number of nucleotides between each k-mer added.
	:param int n_interval: Every this many k-mers, add an N to the k-mer
		sequence to create a k-mer that should not be matched.

	:returns: Tuple of (seq, kmer_spec, kmer_vector).
	"""

	if kmer_interval < k + prefix_len:
		raise ValueError()

	# Seed RNG
	random = np.random.RandomState(seed)

	# Sequence of all ATN's
	seq_array = fill_bytearray(b'ATN', seqlen)

	# Choose prefix with nucleotides not found in sequence "background"
	prefix = bytes(fill_bytearray(b'CGG', prefix_len))
	kspec = kmers.KmerSpec(k, prefix)

	assert kmers.find_kmers(kspec, bytes(seq_array)).sum() == 0

	# Keep track of which kmers have been added
	vec = np.zeros(kspec.idx_len, dtype=bool)

	# Plant matches
	for i in itertools.count():

		p = i * kmer_interval

		if p + kspec.total_len >= seqlen:
			break

		# Use a k-mer without nucleotides in prefix so that we won't create
		# another accidental match
		kmer_array = bytearray(kspec.k)
		kmer_num = random.randint(2 ** k)
		for j in range(kspec.k):
			kmer_array[j] = b'AT'[(kmer_num >> j) & 1]


		# Every so often add an N just to throw things off
		invalid_kmer = n_interval is not None and i % n_interval == 0
		if invalid_kmer:
			kmer_array[0] = ord(b'N')

		kmer = bytes(kmer_array)

		if not invalid_kmer:
			vec[kmers.kmer_to_index(kmer)] = True

		match = kspec.prefix + kmer

		# Reverse every other match
		if i % 2 == 1:
			match = reverse_complement(match)

		seq_array[p:p + kspec.total_len] = match

	return bytes(seq_array), kspec, vec


class TestFindKmers:
	"""Test k-mer finding."""


	def test_basic(self):
		"""Test general k-mer finding."""

		seq, kspec, vec = make_kmer_seq(
			100000, k=11,
			prefix_len=5,
			kmer_interval=50,
			n_interval=10
		)

		# Test normal
		assert np.array_equal(kmers.find_kmers(kspec, seq), vec)

		# Test reverse complement
		assert np.array_equal(kmers.find_kmers(kspec, reverse_complement(seq)), vec)

		# Test lower case
		assert np.array_equal(kmers.find_kmers(kspec, seq.lower()), vec)

		# Test string argument
		assert np.array_equal(kmers.find_kmers(kspec, seq.decode('ascii')), vec)


	def test_bounds(self):
		"""
		Test k-mer finding at beginning and end of sequence to catch errors with
		search bounds.
		"""

		# Sequence of all ATN's
		seqlen = 100000
		seq_array = fill_bytearray(b'ATN', seqlen)

		# Choose prefix with nucleotides not found in sequence "background"
		kspec = kmers.KmerSpec(11, b'CCGGG')

		# Add at beginning
		seq_array[0:kspec.prefix_len] = kspec.prefix
		seq_array[kspec.prefix_len:kspec.total_len] = kmers.index_to_kmer(0, kspec.k)

		# Add at end
		seq_array[-kspec.total_len:-kspec.k] = kspec.prefix
		seq_array[-kspec.k:] = kmers.index_to_kmer(1, kspec.k)

		seq = bytes(seq_array)
		found_coords = kmers.vec_to_coords(kmers.find_kmers(kspec, seq))

		assert np.array_equal(found_coords, [0, 1])


	def test_overlapping(self):
		"""Test k-mer finding when k-mers overlap with each other.

		The test sequence is manually designed to have a variety of overlapping
		forwards and backwards matches
		"""

		kspec = kmers.KmerSpec(11, b'GCCGG')

		seq = b'ATATGCCGGCCGGATTATATAGCCGGCATTACATCCGATAGGATCCGGCAATAA'
		#      |    |>>>>...........
		#      |        |>>>>........... (forward match which overlaps prefix)
		#      |                     |>>>>........... (another overlapping forward match)
		#      |....<<<<| (backward match for prefix, but too close to end)
		#      |           ...........<<<<|
		#      |                                 ...........<<<<|

		expected = {
			b'CCGGATTATAT',
			b'ATTATATAGCC',
			b'CATTACATCCG',
			reverse_complement(b'GGATTATATAG'),
			reverse_complement(b'TCCGATAGGAT'),
		}

		for s in [seq, reverse_complement(seq)]:
			vec = kmers.find_kmers(kspec, s)

			found = [
				kmers.index_to_kmer(i, kspec.k)
				for i in kmers.vec_to_coords(vec)
			]

			assert len(found) == len(expected)
			assert all(kmer in expected for kmer in found)


	def test_parse(self):
		"""Test k-mer finding in FASTA files."""

		from io import StringIO

		from Bio.Seq import Seq
		from Bio.SeqIO import write as write_seqs, SeqRecord

		k = 11

		vec = np.zeros(4 ** k, dtype=bool)
		kspec = None

		records = []

		# Create FASTA file with 10 records
		for i in range(10):

			seq, _kspec, rec_vec = make_kmer_seq(
				100000,
				k=k,
				prefix_len=5,
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
			records.append(SeqRecord(
				seq=Seq(seq.decode('ascii')),
				id='SEQ{}'.format(i + 1),
				description='sequence {}'.format(i + 1),
			))

		# Write records to string buffer in FASTA format
		buf = StringIO()
		write_seqs(records, buf, 'fasta')
		buf.seek(0)

		# Parse from buffer
		parsed_vec = kmers.find_kmers_parse(kspec, buf, 'fasta')

		assert np.array_equal(vec, parsed_vec)


class TestKmerCoordsCollection:
	"""Test KmerCoordsCollection."""

	@pytest.fixture()
	def coords_list(self):
		"""List of k-mer coordinate arrays."""

		random = np.random.RandomState(seed=0)

		ncoords = 100
		nkmers_range = (5000, 15000)

		cl = [
			np.arange(random.randint(*nkmers_range)) * ncoords + i
			for i in range(ncoords)
		]

		# Add a zero-length set in there just to throw everything off
		cl.append(np.empty(0, dtype=cl[0].dtype))

		return cl

	@pytest.fixture()
	def kcol(self, coords_list):
		"""Create from coords_list using from_coords_seq()."""
		return kmers.KmerCoordsCollection.from_coords_seq(coords_list)

	def test_basic(self, kcol, coords_list):
		"""Test item access and calling methods without side effects."""

		# Check correct length
		assert len(kcol) == len(coords_list)

		# Check each coordinate set matches
		for i, coords in enumerate(coords_list):

			assert np.array_equal(kcol[i], coords)
			assert kcol.size_of(i) == len(coords)

	def test_invalid_index(self, kcol):
		"""Test passing invalid indices."""

		# Check negative integer index
		with pytest.raises(IndexError):
			kcol[-1]
		with pytest.raises(IndexError):
			kcol[-1] = 0

		# Test integer index too large
		with pytest.raises(IndexError):
			kcol[len(kcol)]
		with pytest.raises(IndexError):
			kcol[len(kcol)] = 0

		# Test slice
		with pytest.raises(TypeError):
			kcol[:]
		with pytest.raises(TypeError):
			kcol[:] = 0

	def test_from_empty(self, coords_list):
		"""Test creating from empty array."""

		lengths = list(map(len, coords_list))

		kcol = kmers.KmerCoordsCollection.empty(lengths)

		assert len(kcol) == len(coords_list)

		for i in range(len(kcol)):
			assert kcol.size_of(i) == len(coords_list[i])

	def test_assignment(self, kcol):
		"""Test item assignment."""

		for i in range(len(kcol)):

			# Assign to single number
			n = i * 11 % 23
			kcol[i] = n
			assert np.all(kcol[i] == n)

			# Assign to array
			a = (np.arange(kcol.size_of(i)) * 11 + i) % 23
			kcol[i] = a
			assert np.array_equal(kcol[i], a)

	def test_empty(self):
		"""Really an edge case, but test it anyways."""

		kcol = kmers.KmerCoordsCollection.from_coords_seq([])

		assert len(kcol) == 0

		with pytest.raises(IndexError):
			kcol[0]
