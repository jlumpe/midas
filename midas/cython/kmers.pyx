"""Cython module for core k-mers code"""

import numpy as np



cdef class CKmerSpec:
	"""Describes a k-mer pattern to match"""

	def __cinit__(self, int k, bytes prefix):
		self.k = k
		self.prefix = prefix
		self.c_prefix = self.prefix
		self.prefix_len = len(self.prefix)
		self.idx_len = 1 << (2 * self.k)

	def __repr__(self):
		return '<{} k={} prefix="{}"">'.format(self.__class__.__name__, self.k,
		                                       self.prefix)


def kmer_to_index(bytes kmer):
	"""Convert k-mer byte string into index"""
	return c_kmer_to_index(<char*>kmer, len(kmer))


cdef coords_t c_kmer_to_index(char *kmer, int l):
	"""Convert k-mer byte string into index"""

	cdef:
		coords_t idx = 0
		char c

	for i in range(l):
		idx <<= 2

		c = kmer[i] & 0b11011111 # To upper case
		if c == <char>'A':
			idx += 0
		elif c == <char>'C':
			idx += 1
		elif c == <char>'G':
			idx += 2
		elif c == <char>'T':
			idx += 3
		else:
			pass # No good way to indicate error?

	return idx
