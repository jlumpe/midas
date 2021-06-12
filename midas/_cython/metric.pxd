"""Cython functions for calculating k-mer similarity metrics."""

from .types cimport SCORE_T, BOUNDS_T, COORDS_T, COORDS_T_2


cdef SCORE_T c_jaccard_coords(COORDS_T[:] coords1,
                              COORDS_T_2[:] coords2) nogil

cdef void c_jaccard_coords_col(COORDS_T[:] query, COORDS_T_2[:] ref_coords,
                               BOUNDS_T[:] ref_bounds, SCORE_T[:] out) nogil
