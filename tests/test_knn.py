"""Test midas.query.

Note: when comparing Jaccard scores with distances you need to be careful
with exact equality testing - for floating point numbers it is not always true
that 1 - (1 - x) == x. Compare distances to 1 - scores, not the other way
around.
"""

import pytest
import numpy as np

from midas.test import make_signatures
from midas import knn
from midas.kmers import SignatureArray
from midas.metric import BOUNDS_DTYPE, jaccard_sparse, jaccard_sparse_array


def todist(value, convert):
	"""Conditionally convert Jaccard scores to distances.

	Helper function for comparing output of knn functions to expected values,
	which may return Jaccard distances or scores depending on the value of a
	parameter. Converts first argument from score to distance if second arg
	is true, otherwise returns it unchanged.
	"""
	return 1 - value if convert else value


@pytest.fixture(scope='module')
def ref_sigs():
	return make_signatures(8, 40, dtype='i4')


@pytest.fixture(scope='module')
def query_sigs():
	return make_signatures(8, 10, dtype='i4')


@pytest.mark.parametrize('distance', [False, True])
@pytest.mark.parametrize('alt_bounds_dtype', [False, True])
def test_sigarray_scores(query_sigs, ref_sigs, distance, alt_bounds_dtype):
	"""Test midas.knn.sigarray_scores."""

	# The Cython function takes a specific type for the bounds array.
	# Try with this type and a different type
	if alt_bounds_dtype:
		ref_sigs = SignatureArray(ref_sigs.values, ref_sigs.bounds.astype('i4'))
		assert ref_sigs.bounds.dtype != BOUNDS_DTYPE

	else:
		assert ref_sigs.bounds.dtype == BOUNDS_DTYPE

	# Use first signature as query
	querysig = query_sigs[0]

	# Calculate scores
	scores = jaccard_sparse_array(querysig, ref_sigs, distance=distance)

	# Check shape
	assert scores.shape == (len(ref_sigs),)

	# Check scores one at a time
	for refsig, score in zip(ref_sigs, scores):

		expected = jaccard_sparse(querysig, refsig)
		assert score == todist(expected, distance)


@pytest.mark.parametrize('k', [None, 5])
@pytest.mark.parametrize('distance', [False, True])
def test_nn_search(query_sigs, ref_sigs, k, distance):

	query = query_sigs[0]

	indices, scores = knn.nn_search(query, ref_sigs, k=k, distance=distance)

	# Expected shape
	if k is None:
		assert np.isscalar(indices)
		assert np.isscalar(scores)
	else:
		assert indices.shape == (k,)
		assert scores.shape == (k,)

	# Check indices and scores match
	full_scores = jaccard_sparse_array(query, ref_sigs)
	assert np.array_equal(todist(full_scores[indices], distance), scores)

	# Check they are ordered closest to furthest
	if k is not None:
		if distance:
			assert np.all(np.diff(scores) >= 0)
		else:
			assert np.all(np.diff(scores) <= 0)

	# Check they are the k closest
	if k is None:
		assert scores == todist(np.max(full_scores), distance)

	else:
		partition = full_scores[np.argpartition(full_scores, -k)[-k]]

		if distance:
			assert np.all(scores <= 1 - partition)
		else:
			assert np.all(scores >= partition)


@pytest.mark.parametrize('k', [None, 5])
@pytest.mark.parametrize('distance', [False, True])
def test_nn_search_multi(query_sigs, ref_sigs, k, distance):

	indices, scores = knn.nn_search_multi(query_sigs, ref_sigs, k=k, distance=distance)

	# Expected shape
	expected_shape = (len(query_sigs),)
	if k is not None:
		expected_shape += (k,)

	assert indices.shape == expected_shape
	assert scores.shape == expected_shape

	# Check scalar return values aren't just zero-dimensional arrays
	if expected_shape == ():
		assert np.isscalar(indices)
		assert np.isscalar(scores)

	# Expand to 2D to make comparison easier
	if k is None:
		indices = indices[..., None]
		scores = scores[..., None]

	# Check results for each individual query sequence
	for i, qsig in enumerate(query_sigs):

		# All scores from query to ref array (not distances)
		full_scores = jaccard_sparse_array(qsig, ref_sigs)

		# Check indices and scores match
		assert np.array_equal(
			todist(full_scores[indices[i]], distance),
			scores[i]
		)

		# Check they are ordered closest to furthest
		if k is not None:
			if distance:
				assert np.all(np.diff(scores[i]) >= 0)
			else:
				assert np.all(np.diff(scores[i]) <= 0)

		# Check they are the k closest
		if k is None:
			assert scores[i, 0] == todist(np.max(full_scores), distance)

		else:
			partition = full_scores[np.argpartition(full_scores, -k)[-k]]

			if distance:
				assert np.all(scores[i] <= 1 - partition)
			else:
				assert np.all(scores[i] >= partition)
