"""Core functions for searching for and working with k-mers.

Note that all code in this module operates on DNA sequences as sequences of
bytes containing ascii-encoded nucleotide codes.

.. data:: NUCLEOTIDES

	``bytes`` corresponding to the four DNA nucleotides. Ascii-encoded upper
	case letters ``ACGT``. Note that the order, while arbitrary, is important
	in this variable as it defines how unique indices are assigned to k-mer
	sequences.
"""

from collections.abc import Sequence

import numpy as np

from .cython.seqs import kmer_to_index, index_to_kmer, reverse_complement
from midas.util.attr import attrs, attrib
from midas.io.json import Jsonable


# Byte representations of the four nucleotide codes in the order used for
# indexing k-mer sequences
NUCLEOTIDES = b'ACGT'


def validate_dna_seq_bytes(seq):
	"""Check that a sequence contains only valid nucleotide codes.

	Parameters
	----------
	seq : bytes
		ASCII-encoded nucleotide sequence.

	Raises
	------
	ValueError
		If the sequence contains an invalid nucleotide.
	"""
	for i, nuc in enumerate(seq):
		if nuc not in NUCLEOTIDES:
			raise ValueError(f'Invalid byte at position {i}: {nuc}')


def coords_dtype(k):
	"""Get the smallest unsigned integer dtype that can store k-mer indices for the given ``k``.

	Parameters
	----------
	k : int

	Returns
	-------
	numpy.dtype
	"""
	if k <= 4:
		return np.dtype('u1')
	elif k <= 8:
		return np.dtype('u2')
	elif k <= 16:
		return np.dtype('u4')
	elif k <= 32:
		return np.dtype('u8')
	else:
		return None


@attrs(frozen=True, repr=False)
class KmerSpec(Jsonable):
	"""Specifications for a k-mer search operation.

	Parameters
	----------
	k : int
		Value of :attr:`k` attribute.
	prefix : str or bytes
		Value of :attr:`prefix` attribute. If ``str`` and not ``bytes`` will be encoded as ascii.

	Attributes
	----------
	prefix : bytes
		Constant prefix of k-mers to search for, upper-case nucleotide codes
		as ascii-encoded ``bytes``.
	k : int
		Number of nucleotides in k-mer *after* prefix.
	prefix_len : int
		Number of nucleotides in prefix.
	total_len : int
		Sum of ``prefix_len`` and ``k``.
	idx_len : int
		Maximum value (plus one) of integer needed to index one of the
		found k-mers. Also the number of possible k-mers fitting the spec.
		Equal to ``4 ** k``.
	coords_dtype : numpy.dtype
		Smallest unsigned integer dtype that can store k-mer indices.
	"""
	k: int = attrib()
	prefix: bytes = attrib(
		converter=lambda v: v.upper().encode('ascii') if isinstance(v, str) else v,
	)
	prefix_len: int
	total_len: int
	idx_len: int
	coords_dtype: np.dtype

	@k.validator
	def _validate_k(self, attribute, value):
		if value < 1:
			raise ValueError('k must be positive')

	@prefix.validator
	def _validate_prefix(self, attribute, value):
		validate_dna_seq_bytes(value)

	def __attrs_post_init__(self):
		object.__setattr__(self, 'prefix_len', len(self.prefix))
		object.__setattr__(self, 'total_len', self.k + self.prefix_len)
		object.__setattr__(self, 'idx_len', 4 ** self.k)
		object.__setattr__(self, 'coords_dtype', coords_dtype(self.k))

	def __get_newargs__(self):
		return self.k, self.prefix

	def __eq__(self, other):
		return isinstance(other, KmerSpec) and\
			self.k == other.k and\
			self.prefix == other.prefix

	def coords_to_vec(self, coords):
		return coords_to_vec(coords, self.idx_len)

	def __repr__(self):
		return '{}({}, {!r})'.format(
			self.__class__.__name__,
			self.k,
			self.prefix.decode('ascii')
		)

	def __to_json__(self):
		return dict(k=self.k, prefix=self.prefix.decode('ascii'))

	@classmethod
	def __from_json__(cls, data):
		return cls(data['k'], data['prefix'])


def find_kmers(kspec, seq, out=None):
	"""Find k-mers in a sequence and output a coordinate array.

	Searches sequence both backwards and forwards (reverse complement). The
	sequence may contain invalid characters (not one of the four nucleotide
	codes) which will simply not be matched.

	Parameters
	----------
	kspec : .KmerSpec
		K-mer spec to use for search.
	seq
		Sequence to search within as ``bytes`` or ``str``. If ``str``
		will be encoded as ASCII. Lower-case characters are OK and will be
		matched as upper-case.
	out : numpy.ndarray
		Existing numpy array to write output to. Should be of length
		``kspec.idx_len``. If given the same array will be returned.

	Returns
	-------
	numpy.ndarray
		Array of length ``kspec.idx_len`` containing ones in the
		index of each k-mer found and zeros elsewhere. If ``out`` is not
		given the array will be of type ``bool``.
	"""
	if out is None:
		out = np.zeros(kspec.idx_len, dtype=bool)

	# Reverse complement of prefix
	rcprefix = reverse_complement(kspec.prefix)

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
		kmer = reverse_complement(rckmer)

		try:
			out[kmer_to_index(kmer)] = 1
		except ValueError:
			pass

		start = loc + 1

	return out


def vec_to_coords(vec):
	"""Convert boolean k-mer vector to sparse coordinate representation.

	Parameters
	----------
	vec : numpy.ndarray
		Boolean vector indicating which k-mers are present.

	Returns
	-------
	numpy.ndarray
		Sorted array of coordinates of k-mers present in vector. Data type will be ``intp``.
	"""
	return np.flatnonzero(vec)


def coords_to_vec(coords, idx_len):
	"""Convert from sparse coordinate representation back to boolean k-mer vector.

	Parameters
	----------
	coords : numpy.ndarray
		Coordinate array.
	idx_len : int
		Value of ``idx_len`` property of corresponding :class:`.KmerSpec`.

	Returns
	-------
	numpy.ndarray
		Boolean k-mer vector.
	"""
	vec = np.zeros(idx_len, dtype=np.bool_)
	vec[coords] = 1
	return vec


class SignatureArray(Sequence):
	"""
	Stores a collection of k-mer signatures (k-mer sets in sparse coordinate format) in
	a single Numpy array.

	Acts as a non-mutable sequence of :class:`numpy.ndarray` objects (however
	the ndarray's themselves may be writable), with the exception that you may
	assign a value to an index which has the same effect as ``array[:] = value``.
	Only integer indices are supported, not slices.

	Shouldn't create from constructor directly, use :meth:`from_signatures` or
	:meth:`empty` instead.
	"""

	def __init__(self, values, bounds):
		self.values = values
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
		return self.values[self.bounds[index]:self.bounds[index + 1]]

	def __setitem__(self, index, value):
		self._check_index(index)
		self.values[self.bounds[index]:self.bounds[index + 1]] = value

	def sizeof(self, index):
		"""Get the size of the signature array at the given index.

		Should be the case that

		    kcol.size_of(i) == len(kcol[i])

		Parameters
		----------
		index : int
			Index of k-mer set in collection.

		Returns
		-------
		int
		"""
		return self.bounds[index + 1] - self.bounds[index]

	@classmethod
	def from_signatures(cls, signatures, dtype=np.uint32):
		"""Create from sequence of individual signature (k-mer coordinate) arrays.

		Parameters
		----------
		signatures
			Sequence of coordinate arrays.
		dtype : numpy.dtype
			Numpy dtype of shared coordinates array.

		Returns
		-------
		.SignatureArray
		"""
		array = cls.empty(list(map(len, signatures)), dtype=dtype)

		for i, sig in enumerate(signatures):
			array[i] = sig

		return array

	@classmethod
	def empty(cls, lengths, values=None, dtype=np.uint32):
		"""Create with an empty array.

		Parameters
		----------
		lengths
			Sequence of lengths for each sub-array/signature.
		values : numpy.ndarray
			Shared data array to use (optional).
		dtype : numpy.dtype
			Numpy dtype of shared coordinates array.

		Returns
		-------
		.SignatureArray
		"""
		from midas.cython.metrics import BOUNDS_DTYPE

		bounds = np.zeros(len(lengths) + 1, dtype=BOUNDS_DTYPE)
		bounds[1:] = np.cumsum(lengths)

		if values is None:
			values = np.empty(bounds[-1], dtype=dtype)

		return cls(values, bounds)
