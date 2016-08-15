"""Cython module for core k-mers code"""

cimport numpy as np


# Type for storing k-mer coordinates - will work up to k=16 but
# should find a more generic way to do this
ctypedef np.uint32_t coords_t
coords_dtype = np.uint32

# Type for similarity scores
ctypedef np.float32_t score_t
score_dtype = np.float32


cdef class CKmerSpec:

	cdef readonly:
		int k

		bytes prefix
		char* c_prefix
		int prefix_len

		np.intp_t idx_len


cdef coords_t c_kmer_to_index(char*, int)
