"""Cython module for working with DNA sequences and k-mers."""

cimport numpy as np


# Fused type for storing k-mer coordinates
ctypedef fused COORDS_T:
	np.int16_t
	np.uint16_t
	np.int32_t
	np.uint32_t
	np.int64_t
	np.uint64_t

# Copy of COORDS_T, used when two arguments are of this type but not the same
ctypedef fused COORDS_T_2:
	np.int16_t
	np.uint16_t
	np.int32_t
	np.uint32_t
	np.int64_t
	np.uint64_t


cdef np.uint32_t c_kmer_to_index32(const char*, int) except? 0
cdef np.uint64_t c_kmer_to_index64(const char*, int) except? 0
cdef void c_index_to_kmer(COORDS_T, int, char*) nogil
cdef inline char nuc_complement(char) nogil
cdef void c_reverse_complement(const char*, int, char*) nogil
