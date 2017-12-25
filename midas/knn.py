"""Run k-nearest-neighbor queries on a set of reference signatures."""

import numpy as np

from midas.cython import metrics


def sigarray_scores(signature, sigarray, distance=False):
	"""Calculate Jaccard scores between one signature and an array of signatures.

	This internally uses Cython code that runs in parallel over all signatures
	in ``sigarray``.

	:param signature: K-mer signature in coordinate format, increasing sequence
		of integer values.
	:type signature: numpy.ndarray
	:param sigarray: Signature array to calculate scores against.
	:type sigarray: midas.kmers.SignatureArray
	:param bool distance: Return Jaccard distances instead of scores.
	"""

	values = sigarray.values
	bounds = sigarray.bounds.astype(metrics.BOUNDS_DTYPE)

	scores = metrics.jaccard_coords_col(signature, values, bounds)

	return 1 - scores if distance else scores


def find_closest_signatures(query, refarray, *, k=None, distance=False):
	"""Find the closest reference signatures to a query.

	:param query: Single query signature in coordinate format
		(:class:`numpy.ndarray` of increasing integer values) or sequence of
		query signatures (e.g. :class:`midas.kmers.SignatureArray` or list).
	:type query: numpy.ndarray
	:param refarray: Array of reference signatures to calculate scores against.
	:type refarray: midas.kmers.SignatureArray
	:param int k: Number of reference signatures to find for each query. If
		None will only find the closest.
	:param bool distance: Report Jaccard distances instead of scores.

	:returns: ``(index, score)`` tuple giving the index of the closest reference
		signature in the array and the score between the query and the reference.
		Both are scalars or arrays of the same shape. If ``query`` is a sequence
		of signatures then the first axis corresponds to the signatures in the
		query. If ``k`` is not None then the last axis corresponds to reference
		matches in order of decreasing similarity.
	:rtype: tuple
	"""

	# Check k
	if k is not None and not (0 < k <= len(refarray)):
		raise ValueError(
			'k must be > 0 and <= the number of reference signatures'
		)

	if isinstance(query, np.ndarray):
		# Single query

		scores = sigarray_scores(query, refarray)

		if k is None:
			indices = np.argmax(scores)

		else:
			indices = np.argsort(scores)[:-k - 1:-1]

		closest_scores = scores[indices]

		return indices, 1 - closest_scores if distance else closest_scores

	else:
		# Assume sequence of queries - e.g. SignatureArray or list

		out_shape = (len(query),) if k is None else (len(query), k)
		indices = np.zeros(out_shape, dtype=int)
		scores = np.zeros(out_shape, dtype=metrics.SCORE_DTYPE)

		# Call self with individual signatures
		for i, signature in enumerate(query):
			indices[i], scores[i] = find_closest_signatures(
				signature,
				refarray,
				k=k,
				distance=distance
			)

		return indices, scores
