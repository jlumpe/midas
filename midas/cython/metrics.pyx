"""Cython functions for calculating k-mer similarity metrics"""

cimport cython
cimport numpy as np

import numpy as np
from cython.parallel import prange, parallel


# Numpy dtypes equivalent to SCORE_T and BOUNDS_T
SCORE_DTYPE = np.dtype(np.float32)
BOUNDS_DTYPE = np.dtype(np.intp)


def jaccard_coords(COORDS_T[:] coords1, COORDS_T_2[:] coords2):
	"""
	Compute the jaccard index between two k-mer sets in ordered coordinate
	format.

	Data types of array arguments may be 16, 32, or 64-bit signed or unsigned
	integers, but must match.

	:param coords1: k-mer set in ordered coordinate format.
	:type coords1: numpy.ndarray
	:param coords2: k-mer set in ordered coordinate format.
	:type coords2: numpy.ndarray

	:returns: Jaccard index for the two sets.
	:rtype: numpy.float32
	"""
	return c_jaccard_coords(coords1, coords2)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef SCORE_T c_jaccard_coords(COORDS_T[:] coords1,
                              COORDS_T_2[:] coords2) nogil:
	"""
	Compute the Jaccard index between two k-mer sets in ordered coordinate
	format.

	Declared with nogil so it can be run in parallel.
	"""

	cdef:
		# Lengths of the two arrays
		np.intp_t N = coords1.shape[0]
		np.intp_t M = coords2.shape[0]

		# Index and value of items in each array as we are iterating
		np.intp_t i = 0, j = 0
		COORDS_T a
		COORDS_T_2 b

		np.intp_t u = 0  # Size of union

	# Iterate through both arrays simultaneously, advance index for the array
	# with the smaller value. Advance both if they are equal. Increment the
	# union count each loop.
	while i < N and j < M:
		a = coords1[i]
		b = coords2[j]

		u += 1

		if a <= b:
			i += 1

		if b <= a:
			j += 1

	# In most cases we won't have i == N and j == M at the end of the loop,
	# account for the items that we didn't get through
	u += N - i
	u += M - j

	# Avoid divide by zero, define score between empty sets to be zero
	if u == 0:
		return 0

	# |A intersection B| = |A| + |B| - |A union B|
	return <SCORE_T>(N + M - u) / u


def jaccard_coords_col(COORDS_T[:] query,
                       COORDS_T_2[:] ref_coords,
                       BOUNDS_T[:] ref_bounds):
	"""
	Calculate Jaccard scores between a query k-mer set and a reference
	collection.

	Data types of k-mer coordinate arrays may be 16, 32, or 64-bit signed or
	unsigned integers, but must match.

	Internally, releases the GIL in the main loop and calculates scores
	concurrently.

	:param query: Query k-mer set in ordered coordinate format.
	:type query: numpy.ndarray
	:param ref_coords: Reference k-mer sets in ordered coordinate format,
		concatenated into a single array.
	:type ref_coords: numpy.ndarray
	:param ref_bounds: Bounds of individual k-mer sets within the ``ref_coords``
		array. The ``n``\ th k-mer set is the slice of ``ref_coords`` between
		``ref_bounds[n]`` and ``ref_bounds[n + 1]``. Length must be one greater
		than that of ``ref_coords``.
	:type ref_bounds: numpy.ndarray

	:returns: Numpy array of floats, the Jaccard score between the query and
		each reference k-mer set.
	:rtype: numpy.ndarray
	"""

	cdef np.ndarray[SCORE_T, ndim=1] out_array
	cdef SCORE_T[:] out
	out = out_array = np.ndarray(len(ref_bounds) - 1, dtype=SCORE_DTYPE)

	c_jaccard_coords_col(query, ref_coords, ref_bounds, out)

	return out_array


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void c_jaccard_coords_col(COORDS_T[:] query,
                               COORDS_T_2[:] ref_coords,
                               BOUNDS_T[:] ref_bounds,
                               SCORE_T[:] out) nogil:
	"""
	Calculate Jaccard scores between a query k-mer set and a reference
	collection.

	Loop runs in parallel without GIL.
	"""

	cdef np.intp_t N = ref_bounds.shape[0] - 1
	cdef BOUNDS_T begin, end
	cdef int i

	with parallel():
		for i in prange(N, schedule='dynamic'):
			begin = ref_bounds[i]
			end = ref_bounds[i+1]
			out[i] = c_jaccard_coords(query, ref_coords[begin:end])
