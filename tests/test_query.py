"""Test midas.query.

Note: when comparing Jaccard scores with distances you need to be careful
with exact equality testing - for floating point numbers it is not always true
that 1 - (1 - x) == x. Compare distances to 1 - scores, not the other way
around.
"""

import pytest
import numpy as np

from midas.test import make_signatures
from midas import query
from midas.kmers import SignatureArray
from midas.cython import metrics


def todist(value, convert):
	"""Conditionally convert Jaccard scores to distances.

	Helper function for comparing output of query functions to expected values,
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
	"""Test midas.query.sigarray_scores."""

	# The Cython function takes a specific type for the bounds array.
	# Try with this type and a different type
	if alt_bounds_dtype:
		ref_sigs = SignatureArray(ref_sigs.values, ref_sigs.bounds.astype('i4'))
		assert ref_sigs.bounds.dtype != metrics.BOUNDS_DTYPE

	else:
		assert ref_sigs.bounds.dtype == metrics.BOUNDS_DTYPE

	# Use first signature as query
	querysig = query_sigs[0]

	# Calculate scores
	scores = query.sigarray_scores(querysig, ref_sigs, distance=distance)

	# Check shape
	assert scores.shape == (len(ref_sigs),)

	# Check scores one at a time
	for refsig, score in zip(ref_sigs, scores):

		expected = metrics.jaccard_coords(querysig, refsig)
		assert score == todist(expected, distance)


@pytest.mark.parametrize('single_query', [True, False])
@pytest.mark.parametrize('k', [None, 5])
@pytest.mark.parametrize('distance', [False, True])
def test_find_closest_signatures(query_sigs, ref_sigs, single_query, k, distance):

	query_arg = query_sigs[0] if single_query else query_sigs

	indices, scores = query.find_closest_signatures(query_arg, ref_sigs, k=k,
	                                                distance=distance)

	# Expected shape
	expected_shape = ()
	if not single_query:
		expected_shape += (len(query_sigs),)
	if k is not None:
		expected_shape += (k,)

	assert indices.shape == expected_shape
	assert scores.shape == expected_shape

	# Check scalar return values aren't just zero-dimensional arrays
	if expected_shape == ():
		assert np.isscalar(indices)
		assert np.isscalar(scores)

	# Expand to 2D to make comparison easier
	if single_query:
		indices = indices[None, ...]
		scores = scores[None, ...]
	if k is None:
		indices = indices[..., None]
		scores = scores[..., None]

	# Check results for each individual query sequence
	for i, qsig in enumerate([query_sigs[0]] if single_query else query_sigs):

		# All scores from query to ref array (not distances)
		full_scores = query.sigarray_scores(qsig, ref_sigs)

		# Check indices and scores match
		assert np.array_equal(
			todist(full_scores[indices[i]], distance),
			scores[i]
		)

		# Check they are ordered closest to furthest
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
