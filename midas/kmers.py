"""Core functions for searching for and working with k-mers.

Note that all code in this module operates on DNA sequences as sequences of
bytes containing ascii-encoded nucleotide codes.

.. data:: NUCLEOTIDES

	``bytes`` corresponding to the four DNA nucleotides. Ascii-encoded upper
	case letters ``ACGT``. Note that the order, while arbitrary, is important
	in this variable as it defines how unique indices are assigned to k-mer
	sequences.
"""

import collections

import numpy as np
from Bio import SeqIO

from .cython import kmers as ckmers
from .cython.kmers import kmer_to_index


# Byte representations of the four nucleotide codes in the order used for
# indexing k-mer sequences
NUCLEOTIDES = b'ACGT'


def reverse_complement(seq):
	"""Reverse complement of a (short) upper-case byte sequence"""
	complement = []
	for nuc in seq:
		if nuc == 65:  # A
			complement.append(84)
		elif nuc == 84:  # T
			complement.append(65)
		elif nuc == 67:  # C
			complement.append(71)
		elif nuc == 71:  # G
			complement.append(67)
	return bytes(reversed(complement))


class KmerSpec(object):
	"""Specifications for a k-mer search operation.

	:param int k: Value of :attr:`k` attribute.
	:param prefix: Value of :attr:`prefix` attribute. If ``str`` and not
		``bytes`` will be encoded as ascii.

	.. attribute:: prefix

		Constant prefix of k-mers to search for, upper-case nucleotide codes
		as ascii-encoded ``bytes``.

	.. attribute:: k

		Number of nucleotides in k-mer *after* prefix.

	.. attribute:: prefix_len

		Number of nucleotides in prefix.

	.. attribute:: total_len

		Sum of ``prefix_len`` and ``k``.

	.. attribute:: idx_len

		Maximum value (plus one) of integer needed to index one of the
		found k-mers. Also the number of possible k-mers fitting the spec.
		Equal to ``4 ** k``.

	.. attribute:: coords_dtype

		Smallest unsigned integer ``numpy.dtype`` that can store k-mer
		indices.
	"""

	def __init__(self, k, prefix):
		"""
		Args:
			k: Length of k-mers to find (excluding prefix).
			prefix: str. Find k-mers beginning with this subsequence.
		"""
		self.k = k

		if isinstance(prefix, str):
			self.prefix = prefix.upper().encode('ascii')
		else:
			self.prefix = bytes(prefix)

		for nuc in self.prefix:
			if nuc not in NUCLEOTIDES:
				raise ValueError('Invalid byte in prefix: {}'.format(nuc))

		self.prefix_len = len(self.prefix)
		self.total_len = self.k + self.prefix_len
		self.idx_len = 4 ** self.k

	def coords_to_vec(self, coords):
		return coords_to_vec(coords, self.idx_len)

	@property
	def coords_dtype(self):
		if self.k <= 4:
			strtype = 'u1'
		elif self.k <= 8:
			strtype = 'u2'
		elif self.k <= 12:
			strtype = 'u4'
		elif self.k <= 16:
			strtype = 'u8'
		else:
			return None

		return np.dtype(strtype)

	def __repr__(self):
		return '<{} k={} prefix="{}">'.format(
			self.__class__.__name__,
			self.k,
			self.prefix.decode('ascii')
		)


def find_kmers(kspec, seq, out=None):
	"""Find k-mers in a sequence and output a coordinate array.

	Searches sequence both backwards and forwards (reverse complement). The
	sequence may contain invalid characters (not one of the four nucleotide
	codes) which will simply not be matched.

	:param kspec: K-mer spec to use for search.
	:type kspec: .KmerSpec
	:param seq: Sequence to search within as ``bytes`` or ``str``. If ``str``
		will be encoded as ASCII.
	:param out: Existing numpy array to write output to. Should be of length
		``kspec.idx_len``. If given the same array will be returned.
	:type out: numpy.ndarray
	:returns: Array of length ``kspec.idx_len`` containing ones in the
		index of each k-mer found and zeros elsewhere. If ``out`` is not
		given the array will be of type ``bool``.
	:rtype: numpy.ndarray
	"""
	if out is None:
		out = np.zeros(kspec.idx_len, dtype=bool)

	rcprefix = reverse_complement(kspec.prefix)

	# Seq as bytes
	if isinstance(seq, bytes):
		prefix = kspec.prefix
		rcprefix = rcprefix
	else:
		prefix = kspec.prefix.decode('ascii')
		rcprefix = rcprefix.decode('ascii')

	# Search forward
	start = 0
	while True:
		loc = seq.find(prefix, start, -kspec.k)
		if loc < 0:
			break

		kmer = seq[loc + kspec.prefix_len:loc + kspec.total_len]
		if not isinstance(kmer, bytes):
			kmer = str(kmer).encode('ascii')

		try:
			out[kmer_to_index(kmer)] = 1
		except ValueError:
			pass

		start = loc + 1

	# Search backward
	start = kspec.k
	while True:
		loc = seq.find(rcprefix, start)
		if loc < 0:
			break

		rckmer = seq[loc - kspec.k:loc]
		if not isinstance(rckmer, bytes):
			rckmer = str(rckmer).encode('ascii')
		kmer = reverse_complement(rckmer)

		try:
			out[kmer_to_index(kmer)] = 1
		except ValueError:
			pass

		start = loc + 1

	return out


def find_kmers_parse(kspec, data, format, out=None):
	"""Parses sequence data with Bio.Seq.parse() and finds k-mers"""

	if out is None:
		out = np.zeros(kspec.idx_len, dtype=bool)

	for record in SeqIO.parse(data, format):
		find_kmers(kspec, record.seq, out=out)

	return out


def kmer_at_index(index, k):
	nucs_reversed = []
	for i in range(k):
		nucs_reversed.append(NUCLEOTIDES[index % 4])
		index >>= 2
	return bytes(reversed(nucs_reversed))


def vec_to_coords(vec):
	"""Convert to compressed coordinate representation"""
	coords, = np.nonzero(vec)
	return coords


def coords_to_vec(coords, idx_len):
	"""Convert from coordinate representation back to vector"""
	vec = np.zeros(idx_len, dtype=np.bool)
	vec[coords] = 1
	return vec


class KmerCoordsCollection(collections.Sequence):
	"""Stores a collection of k-mer sets in coordinate format in a single array"""

	def __init__(self, coords_array, bounds):
		self.coords_array = coords_array
		self.bounds = bounds

	def __len__(self):
		return len(self.bounds) - 1

	def __getitem__(self, index):
		return self.coords_array[self.bounds[index]:self.bounds[index + 1]]

	def __setitem__(self, index, value):
		self.coords_array[self.bounds[index]:self.bounds[index + 1]] = value

	def size_of(self, index):
		return self.bounds[index + 1] - self.bounds[index]

	@classmethod
	def from_coords_seq(cls, coord_seq, dtype=np.uint32):
		coords_col = cls.empty(list(map(len, coord_seq)), dtype=dtype)

		for i, coords in enumerate(coord_seq):
			coords_col[i] = coords

		return coords_col

	@classmethod
	def empty(cls, lengths, coords_array=None, dtype=np.uint32):
		bounds = np.zeros(len(lengths) + 1, dtype=dtype)
		bounds[1:] = np.cumsum(lengths)

		if coords_array is None:
			coords_array = np.empty(bounds[-1], dtype=dtype)

		return cls(coords_array, bounds)
