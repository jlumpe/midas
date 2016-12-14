"""Module for core k-mers code"""

import collections

import numpy as np
from Bio import SeqIO

from .cython import kmers as ckmers


kmer_to_index = ckmers.kmer_to_index


nucleotides = b'ACGT'


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
	"""Specifications for a k-mer search operation"""

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

		self.prefix_len = len(self.prefix)
		self.total_len = self.k + self.prefix_len
		self.idx_len = 4 ** self.k

	def coords_to_vec(self, coords):
		return coords_to_vec(coords, self.idx_len)

	@property
	def coords_dtype(self):
		"""Smallest unsigned unteger numpy dtype that can store coordinates"""
		if self.k <= 4:
			return 'u1'
		elif self.k <= 8:
			return 'u2'
		elif self.k <= 12:
			return 'u4'
		elif self.k <= 16:
			return 'u8'
		else:
			return None


def find_kmers(kspec, seq, out=None):
	"""Find k-mers in a str- or bytes-like sequence, output a coordinate array"""
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
		nucs_reversed.append(nucleotides[index % 4])
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
