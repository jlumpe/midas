"""Cython functions for calculating k-mer similarity metrics."""

import numpy as np
cimport numpy as np

from .seqs cimport coords_t


# Type for similarity scores
ctypedef np.float32_t score_t


cdef score_t c_jaccard_coords(coords_t[:] coords1,
                              coords_t[:] coords2) nogil

cdef void c_jaccard_coords_col(coords_t[:] query, coords_t[:] ref_coords,
                               np.intp_t[:] ref_bounds, score_t[:] out)
