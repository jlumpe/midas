"""Cython module for working with DNA sequences and k-mers."""

cimport numpy as np


# Type for storing k-mer coordinates - will work up to k=16 but
# should find a more generic way to do this
ctypedef np.uint32_t coords_t


cdef class CKmerSpec:

	cdef readonly:
		int k

		bytes prefix
		char* c_prefix
		int prefix_len

		np.intp_t idx_len


cdef coords_t c_kmer_to_index(const char*, int) except? 0
cdef void c_index_to_kmer(coords_t, int, char*) nogil
cdef inline char nuc_complement(char) nogil
cdef void c_reverse_complement(const char*, int, char*) nogil
