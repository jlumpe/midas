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

from .cython import seqs as cseqs
from .cython.seqs import kmer_to_index, index_to_kmer


# Byte representations of the four nucleotide codes in the order used for
# indexing k-mer sequences
NUCLEOTIDES = b'ACGT'


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

	def __eq__(self, other):
		return isinstance(other, KmerSpec) and\
			self.k == other.k and\
			self.prefix == other.prefix

	def coords_to_vec(self, coords):
		return coords_to_vec(coords, self.idx_len)

	@property
	def coords_dtype(self):
		if self.k <= 4:
			strtype = 'u1'
		elif self.k <= 8:
			strtype = 'u2'
		elif self.k <= 16:
			strtype = 'u4'
		elif self.k <= 32:
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
		will be encoded as ASCII. Lower-case characters are OK ans will be
		matched as upper-case.
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

	# Reverse complement of prefix
	rcprefix = cseqs.reverse_complement(kspec.prefix)

	# Convert sequence to bytes
	if not isinstance(seq, bytes):
		if not isinstance(seq, str):
			seq = str(seq)

		seq = seq.encode('ascii')

	# Convert to upper-case only if needed
	nucs_lower = set(NUCLEOTIDES.lower())
	for char in seq:
		if char in nucs_lower:
			seq = seq.upper()
			break

	# Search forward
	start = 0
	while True:
		loc = seq.find(kspec.prefix, start, -kspec.k)
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
		kmer = cseqs.reverse_complement(rckmer)

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


def vec_to_coords(vec):
	"""Convert k-mer vector to compressed coordinate representation.

	:param vec: Boolean vector indicating which k-mers are present.
	:type vec: numpy.ndarray
	:returns: Sorted array of coordinates of k-mers present in vector.
	:rtype: numpy.ndarray
	"""
	return np.flatnonzero(vec)


def coords_to_vec(coords, idx_len):
	"""Convert from coordinate representation back to k-mer vector.

	:param coords: Coordinate array.
	:type coords: numpy.ndarray
	:param int idx_len: Value of ``idx_len`` property of corresponding
		:class:`.KmerSpec`.
	:returns: Boolean k-mer vector.
	:rtype: numpy.ndarray
	"""
	vec = np.zeros(idx_len, dtype=np.bool)
	vec[coords] = 1
	return vec


class KmerCoordsCollection(collections.Sequence):
	"""Stores a collection of k-mer sets in coordinate format in a single array.

	Shouldn't create from constructor directly, use :meth:`from_coords_seq` or
	:meth:`empty` instead.
	"""

	def __init__(self, coords_array, bounds):
		self.coords_array = coords_array
		self.bounds = bounds

	def __len__(self):
		return len(self.bounds) - 1

	def _check_index(self, index):
		"""Check an index passed as an argument is valid."""
		if not isinstance(index, (int, np.integer)):
			raise TypeError('Index must be single integer')
		elif not 0 <= index < len(self):
			raise IndexError('Index out of bounds: {}'.format(index))

	def __getitem__(self, index):
		self._check_index(index)
		return self.coords_array[self.bounds[index]:self.bounds[index + 1]]

	def __setitem__(self, index, value):
		self._check_index(index)
		self.coords_array[self.bounds[index]:self.bounds[index + 1]] = value

	def size_of(self, index):
		"""Get the size of the coordinate set at the given index.

		Should be the case that

		    kcol.size_of(i) == len(kcol[i])

		:param int index: Index of k-mer set in collection.
		"""
		return self.bounds[index + 1] - self.bounds[index]

	@classmethod
	def from_coords_seq(cls, coord_seq, dtype=np.uint32):
		"""Create from sequence of individual coordinate arrays.

		:param coord_seq: Sequence of coordinate arrays.
		:param dtype: Numpy dtype of shared coordinates array.
		:type dtype: numpy.dtype

		:rtype: KmerCoordsCollection
		"""

		coords_col = cls.empty(list(map(len, coord_seq)), dtype=dtype)

		for i, coords in enumerate(coord_seq):
			coords_col[i] = coords

		return coords_col

	@classmethod
	def empty(cls, lengths, coords_array=None, dtype=np.uint32):
		"""Create with an empty array.

		:param lengths: Sequence of lengths for each sub-array.
		:param coords_array: Initial shared coordinates array.
		:type coords_array: numpy.ndarray
		:param dtype: Numpy dtype of shared coordinates array.
		:type dtype: numpy.dtype

		:rtype: KmerCoordsCollection
		"""

		bounds = np.zeros(len(lengths) + 1, dtype=np.intp)
		bounds[1:] = np.cumsum(lengths)

		if coords_array is None:
			coords_array = np.empty(bounds[-1], dtype=dtype)

		return cls(coords_array, bounds)
