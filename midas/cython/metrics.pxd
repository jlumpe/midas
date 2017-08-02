"""Cython functions for calculating k-mer similarity metrics."""

import numpy as np
cimport numpy as np

from .seqs cimport COORDS_T, COORDS_T_2


# Type for similarity scores
ctypedef np.float32_t SCORE_T

# Type for bounds on c_jaccard_coords_col
ctypedef np.intp_t BOUNDS_T


cdef SCORE_T c_jaccard_coords(COORDS_T[:] coords1,
                              COORDS_T_2[:] coords2) nogil

cdef void c_jaccard_coords_col(COORDS_T[:] query, COORDS_T_2[:] ref_coords,
                               BOUNDS_T[:] ref_bounds, SCORE_T[:] out) nogil
