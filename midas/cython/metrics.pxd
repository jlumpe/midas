"""Cython functions for calculating k-mer similarity metrics"""

from .kmers cimport coords_t


cdef float c_jaccard_coords(coords_t[:] coords1,
                            coords_t[:] coords2) nogil

cdef void c_jaccard_coords_col(coords_t[:] query, coords_t[:] ref_coords,
                               coords_t[:] ref_bounds, float[:] out)
