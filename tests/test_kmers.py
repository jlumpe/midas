"""Tests for midas.kmers module."""

import pytest
import numpy as np

from midas import kmers
from midas._cython.kmers import reverse_complement
import midas.io.json as mjson
from midas.test import fill_bytearray, make_kmer_seq


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


def test_kmer_coords_dtype():
	"""Test coords_dtype function."""

	# Try k from 0 to 32 (all have dtypes)
	for k in range(33):
		# Check dtype can store the largest index
		top_idx = 4**k - 1
		assert kmers.coords_dtype(k).type(top_idx) == top_idx

	# k > 32 should have no dtype
	assert kmers.coords_dtype(33) is None


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

	# Check invalid raises error
	with pytest.raises(ValueError):
		kmers.kmer_to_index(b'ATGNC')


class TestKmerSpec:
	"""Test midas.kmers.KmerSpec."""

	def test_constructor(self):
		# Prefix conversion
		assert kmers.KmerSpec(11, b'ATGAC').prefix == b'ATGAC'
		assert kmers.KmerSpec(11, 'ATGAC').prefix == b'ATGAC'
		assert kmers.KmerSpec(11, 'atgac').prefix == b'ATGAC'

		# Invalid prefix
		with pytest.raises(ValueError):
			kmers.KmerSpec(11, b'ATGAX')
			kmers.KmerSpec(11, 'ATGAX')
			kmers.KmerSpec(11, b'ATGAc')

	def test_attributes(self):
		"""Test basic attributes."""

		# Try k from 1 to 32 (all have dtypes)
		for k in range(1, 33):

			spec = make_kmerspec(k)

			# Check length attributes
			assert spec.prefix_len == len(spec.prefix)
			assert spec.prefix_len + spec.k == spec.total_len

		# Check prefix is bytes
		assert isinstance(kmers.KmerSpec(11, 'ATGAC').prefix, bytes)

	def test_eq(self):
		"""Test equality testing."""

		kspec = kmers.KmerSpec(11, 'ATGAC')
		assert kspec == kmers.KmerSpec(11, 'ATGAC')
		assert kspec != kmers.KmerSpec(11, 'ATGAA')
		assert kspec != kmers.KmerSpec(12, 'ATGAC')

	def test_pickle(self):

		import pickle

		kspec = kmers.KmerSpec(11, 'ATGAC')

		assert kspec == pickle.loads(pickle.dumps(kspec))

	def test_json(self):
		"""Test conversion to/from JSON."""

		kspec = kmers.KmerSpec(11, 'ATGAC')
		data = mjson.to_json(kspec)

		assert data == dict(
			k=kspec.k,
			prefix=kspec.prefix.decode('ascii'),
		)

		assert mjson.from_json(data, kmers.KmerSpec) == kspec


def test_vec_coords_conversion():
	"""Test conversion between k-mer vector and coordinates."""

	for k in range(1, 10):

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
	"""Assert the reverse complement of a sequence is correct.

	:param bytes seq: Byte sequence
	:param bytes rc: Reverse complement of seq
	"""
	l = len(seq)
	for i in range(l):
		assert rc[l - i - 1] == NUC_COMPLEMENTS.get(seq[i], seq[i])


def test_revcomp():
	"""Test midas._cython.kmers.reverse_complement."""

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


class TestFindKmers:
	"""Test k-mer finding."""

	def test_basic(self):
		"""Test general k-mer finding."""

		seq, kspec, vec = make_kmer_seq(
			100000,
			k=11,
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


class TestSignatureArray:
	"""Test midas.kmers.SignatureArray."""

	@pytest.fixture()
	def sigs_list(self):
		"""List of k-mer signatures as coordinate arrays."""

		random = np.random.RandomState(seed=0)

		nsigs = 100
		nkmers_range = (5000, 15000)

		cl = [
			np.arange(random.randint(*nkmers_range)) * nsigs + i
			for i in range(nsigs)
		]

		# Add a zero-length set in there just to throw everything off
		cl.append(np.empty(0, dtype=cl[0].dtype))

		return cl

	@pytest.fixture()
	def kcol(self, sigs_list):
		"""Create from sigs_list using from_signatures()."""
		return kmers.SignatureArray.from_signatures(sigs_list)

	def test_basic(self, kcol, sigs_list):
		"""Test item access and calling methods without side effects."""

		# Check correct length
		assert len(kcol) == len(sigs_list)

		# Check each coordinate set matches
		for i, sig in enumerate(sigs_list):

			# Test using Python and numpy integer types as index
			for index in [i, np.int64(i)]:

				assert np.array_equal(kcol[index], sig)
				assert kcol.sizeof(index) == len(sig)

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

	def test_from_empty(self, sigs_list):
		"""Test creating from empty array."""

		lengths = list(map(len, sigs_list))

		kcol = kmers.SignatureArray.empty(lengths)

		assert len(kcol) == len(sigs_list)

		for i in range(len(kcol)):
			assert kcol.sizeof(i) == len(sigs_list[i])

	def test_assignment(self, kcol):
		"""Test item assignment."""

		for i in range(len(kcol)):

			# Assign to single number
			n = i * 11 % 23
			kcol[i] = n
			assert np.all(kcol[i] == n)

			# Assign to array
			a = (np.arange(kcol.sizeof(i)) * 11 + i) % 23
			kcol[i] = a
			assert np.array_equal(kcol[i], a)

	def test_empty(self):
		"""Really an edge case, but test it anyways."""

		kcol = kmers.SignatureArray.from_signatures([])

		assert len(kcol) == 0

		with pytest.raises(IndexError):
			kcol[0]
