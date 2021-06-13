"""Run k-nearest-neighbor queries on a set of reference signatures."""

import numpy as np

from midas.metric import SCORE_DTYPE, jaccard_sparse_array


def nn_search(query_sig, refarray, k=None, distance=False):
	"""Find the closest reference signatures to a query signature.

	Parameters
	----------
	query_sig : numpy.ndarray
		Single query signature in coordinate format (:class:`numpy.ndarray` of increasing integer values).
	refarray : midas.kmers.SignatureArray
		Array of reference signatures to calculate scores against.
	k : int
		Number of reference signatures to find for each query. If None will only find the closest.
	distance : bool
		Report Jaccard distances instead of scores.

	Returns
	-------
	tuple
		``(index, score)`` tuple giving the index of the closest reference
		signature in the array and the score between the query and the reference.
		If ``k`` is None both will be scalars, otherwise each will be an array
		of length ``k`` and correspond to reference matches in order of
		decreasing similarity.
	"""

	# Check k
	if k is not None and not (0 < k <= len(refarray)):
		raise ValueError(
			'k must be > 0 and <= the number of reference signatures'
		)

	scores = jaccard_sparse_array(query_sig, refarray)

	if k is None:
		indices = np.argmax(scores)

	else:
		indices = np.argsort(scores)[:-k - 1:-1]

	closest_scores = scores[indices]

	return indices, 1 - closest_scores if distance else closest_scores


def nn_search_multi(query_sigs, refarray, k=None, distance=False,
                    query_size=None, has_index=False, progress=None):
	"""Find the closest reference signatures to a collection of query signatures.

	Parameters
	----------
	query
		Iterable of query signatures (e.g. :class:`midas.kmers.SignatureArray` or list).
	refarray : midas.kmers.SignatureArray
		Array of reference signatures to calculate scores against.
	k : int
		Number of reference signatures to find for each query. If None will only find the closest.
	distance : bool
		Report Jaccard distances instead of scores.
	query_size : int
		Size of query if it does not implement the ``len()`` protocol.
	has_index : bool
		If ``query_sigs`` actually contains ``(index, signature)`` pairs indicating which index in
		the output arrays the results of each signature belongs in.
	progress : callable
		Optional callable to report progress. Called with the next index after each query signature
		processed.

	Returns
	-------
	tuple
		``(index, score)`` tuple giving the index of the closest reference
		signature(s) in the array to each query signature and the distance
		between the query signature and the corresponding reference signature(s).
		Both are scalars or arrays of the same shape. If ``query`` is a sequence
		If ``k`` is not None then the last axis corresponds to reference matches
		in order of decreasing similarity.
	"""

	if query_size is None:
		query_size = len(query_sigs)

	# Create output arrays
	out_shape = (query_size,) if k is None else (query_size, k)
	indices = np.zeros(out_shape, dtype=int)
	scores = np.zeros(out_shape, dtype=SCORE_DTYPE)

	p = 0

	# Process individual signatures
	for i, signature in (query_sigs if has_index else enumerate(query_sigs)):

		indices[i], scores[i] = nn_search(
			signature,
			refarray,
			k=k,
			distance=distance
		)

		if progress is not None:
			p += 1
			progress(p)

	return indices, scores
