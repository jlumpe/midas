"""Cython functions for calculating k-mer similarity metrics"""

cimport cython
cimport numpy as np

import numpy as np
from cython.parallel import prange, parallel


# Numpy dtype equivalent to score_t
score_dtype = np.dtype(np.float32)


def jaccard_coords(coords_t[:] coords1, coords_t[:] coords2):
	"""Python wrapper for c_jaccard_coords"""
	return c_jaccard_coords(coords1, coords2)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef score_t c_jaccard_coords(coords_t[:] coords1,
                              coords_t[:] coords2) nogil:
	"""Compute the jaccard index between two k-mer sets in ordered coordinate
	format

	Declared with nogil so it can be run in parallel
	"""

	cdef:
		np.intp_t N = coords1.shape[0]
		np.intp_t M = coords2.shape[0]

		np.intp_t i = 0, j = 0
		coords_t a, b

		np.intp_t u = 0

	while i < N and j < M:
		a = coords1[i]
		b = coords2[j]

		u += 1

		if a <= b:
			i += 1

		if b <= a:
			j += 1

	u += N - i
	u += M - j

	return <score_t>(N + M - u) / u


def jaccard_coords_col(coords_t[:] query,
                       coords_t[:] ref_coords,
                       np.intp_t[:] ref_bounds):
	"""Python wrapper for c_jaccard_coords_col"""

	cdef np.ndarray[score_t, ndim=1] out
	out = np.ndarray(len(ref_bounds) - 1, dtype=score_dtype)

	c_jaccard_coords_col(query, ref_coords, ref_bounds, out)

	return out


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void c_jaccard_coords_col(coords_t[:] query,
                               coords_t[:] ref_coords,
                               np.intp_t[:] ref_bounds,
                               score_t[:] out):
	"""Calculate jaccard scores between one k-mer set and a collection

	Runs in parallel without GIL
	"""

	cdef np.intp_t N = ref_bounds.shape[0] - 1
	cdef np.intp_t begin, end
	cdef int i

	with nogil, parallel():
		for i in prange(N, schedule='dynamic'):
			begin = ref_bounds[i]
			end = ref_bounds[i+1]
			out[i] = c_jaccard_coords(query, ref_coords[begin:end])
