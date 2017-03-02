"""Cython module for working with DNA sequences and k-mers."""

import numpy as np

from libc.stdlib cimport malloc, free


cdef class CKmerSpec:
	"""Describes a k-mer pattern to match"""

	def __cinit__(self, int k, bytes prefix):
		self.k = k
		self.prefix = prefix
		self.c_prefix = self.prefix
		self.prefix_len = len(self.prefix)
		self.idx_len = 1 << (2 * self.k)

	def __repr__(self):
		return '<{} k={} prefix="{}"">'.format(
			self.__class__.__name__,
			self.k,
			self.prefix
		)


def kmer_to_index(bytes kmer):
	"""kmer_to_index(bytes kmer)

	Convert k-mer byte string into index.

	:param bytes kmer: K-mer as bytes string.
	:returns: K-mer index as appropriate numpy integer type.
	"""
	return c_kmer_to_index(<char*>kmer, len(kmer))


def index_to_kmer(coords_t index, int k):
	"""index_to_kmer(coords_t index, int k)

	Convert k-mer index to sequence.

	:param index: K-mer index.
	:param k: Length of k-mer.

	:returns: K-mer as byte string.
	:rtype: bytes
	"""

	cdef char* buf = <char*>malloc((k + 1) * sizeof(char))

	try:

		c_index_to_kmer(index, k, buf)
		buf[k] = 0  # Null-terminate it
		return <bytes>buf

	finally:
		free(buf)


def reverse_complement(bytes seq):
	"""reverse_complement(bytes seq)

	Get the reverse complement of a nucleotide sequence.

	:param bytes seq: ASCII-encoded nucleotide sequence. Case does not matter.
	:returns: Reverse complement sequence. All characters in the input which are
		not valid nucleotide codes will appear unchanged in the cooresponding
		reverse position.
	:rtype: bytes
	"""

	cdef:
		int l = len(seq)
		char* buf = <char*>malloc((l + 1) * sizeof(char))

	try:
		c_reverse_complement(seq, l, buf)
		buf[l] = 0  # Null-terminate it
		return <bytes>buf

	finally:
		free(buf)


cdef coords_t c_kmer_to_index(const char *kmer, int k) except? 0:
	"""Convert k-mer byte string into index.

	:param kmer: Pointer to k-mer string.
	:param k: Length of k-mer string.
	:returns: Index of k-mer

	:raises ValueError: If an invalid nucleotide code is encountered.
	"""

	cdef:
		coords_t idx = 0
		char c

	for i in range(k):
		idx <<= 2

		c = kmer[i] & 0b11011111  # To upper case
		if c == <char>'A':
			idx += 0
		elif c == <char>'C':
			idx += 1
		elif c == <char>'G':
			idx += 2
		elif c == <char>'T':
			idx += 3
		else:
			raise ValueError(kmer[i])

	return idx


cdef void c_index_to_kmer(coords_t index, int k, char* out) nogil:
	"""Convert k-mer index to sequence.

	Output will be in uppper-case characters.

	:param index: Index of k-mer
	:param k: Length of k-mer
	:param out: Destination buffer of length k sequence will be written to.
	"""
	cdef:
		int i, nuc_index
		char nuc

	for i in range(k):

		nuc_index = index % 4

		if nuc_index == 0:
			nuc = 'A'
		elif nuc_index == 1:
			nuc = 'C'
		elif nuc_index == 2:
			nuc = 'G'
		else:
			nuc = 'T'

		out[k - i - 1] = nuc

		index = index / 4


cdef inline char nuc_complement(char nuc) nogil:
	"""Get the complement of a nucleotide.

	If the input is a valid nucleotide code the output will be the code of the
	complement nucleotide in the same case. If the input is not a valid
	nucleotide code it will be returned unchanged.
	"""
	if nuc == 'A':
		return 'T'
	elif nuc == 'a':
		return 't'
	elif nuc == 'T':
		return 'A'
	elif nuc == 't':
		return 'a'
	elif nuc == 'G':
		return 'C'
	elif nuc == 'g':
		return 'c'
	elif nuc == 'C':
		return 'G'
	elif nuc == 'c':
		return 'g'
	else:
		return nuc


cdef void c_reverse_complement(const char* seq, int l, char* out) nogil:
	"""Get the reverse complement of a nucleotide sequence.

	:param seq: Pointer to start of sequence
	:param l: Length of sequence
	:out: Pointer to buffer of same length as seq to store the output.
	"""

	cdef int i

	for i in range(l):
		out[l - i - 1] = nuc_complement(seq[i])