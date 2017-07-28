"""Test distance metric calcuations."""

import pickle

import pytest
import numpy as np

from midas.cython.metrics import jaccard_coords, jaccard_coords_col
from midas.kmers import SignatureArray, vec_to_coords, coords_to_vec


def slow_jaccard_coords(coords1, coords2, idx_len):
	"""Get jaccard score of two k-mer coordiate arrays the slow but correct way."""

	vec1 = coords_to_vec(coords1, idx_len)
	vec2 = coords_to_vec(coords2, idx_len)

	intersection = (vec1 & vec2).sum()
	union = (vec1 | vec2).sum()

	return intersection / union


def make_signatures(k, nsets, dtype):
	"""Make artificial k-mer signatures.

	:rtype: SignatureArray
	"""

	random = np.random.RandomState(seed=0)

	idx_len = 4 ** k

	signatures_list = []

	# Add empty and full sets as edge cases
	signatures_list.append(np.arange(0))
	signatures_list.append(np.arange(idx_len))

	# Use a core set of k-mers so that we get some overlap
	core_prob = max(0.01, 20 / idx_len)
	core_vec = random.rand(idx_len) < core_prob

	for i in range(nsets - 2):

		keep_core = random.rand(idx_len) < (random.rand() ** 2)

		new_vec = random.rand(idx_len) < random.uniform(core_prob * .1, core_prob)

		vec = (keep_core & core_vec) | new_vec
		signatures_list.append(vec_to_coords(vec))

	return SignatureArray.from_signatures(signatures_list, dtype=dtype)


@pytest.fixture(scope='module')
def load_test_coords_col(fixture_file):
	"""Function to load k-mer coordinates from file."""

	def load_test_coords_col_func():

		with fixture_file('kmer_coords/coords.pickle').open('rb') as fobj:
			signatures_list = pickle.load(fobj)

		return SignatureArray.from_signatures(signatures_list)

	return load_test_coords_col_func


@pytest.fixture(params=[
	(4, 'u2'),
	(4, 'i2'),
	(7, 'u2'),
	(7, 'i2'),
	(9, 'u4'),
	(9, 'i4'),
	(9, 'u8'),
	(9, 'i8'),
	None,
])
def coords_params(request, load_test_coords_col):
	"""Tuple of (k, SignatureArray) to test on."""

	if request.param is None:
		# Load coords from file
		k = 11
		sigs = load_test_coords_col()

	else:
		# Create coords
		k, dtype = request.param
		dtype = np.dtype(dtype)

		sigs = make_signatures(k, 25, dtype)

	return k, sigs


def test_jaccard_single(coords_params):
	"""Test calculating single scores at a time."""

	k, sigs = coords_params
	idx_len = 4 ** k

	# Iterate over all pairs
	for i, coords1 in enumerate(sigs):

		for j in range(i, len(sigs)):

			coords2 = sigs[j]

			# Skip empty set vs itself because that is undefined
			if len(coords1) == 0 and len(coords2) == 0:
				continue

			# Calc score
			score = jaccard_coords(coords1, coords2)

			# Check range
			assert 0 <= score <= 1

			# Check vs slow version
			assert np.isclose(score, slow_jaccard_coords(coords1, coords2, idx_len))

			# Check with arguments swapped
			assert score == jaccard_coords(coords2, coords1)

			# Check score vs. self is one
			if i == j:
				assert score == 1


def test_jaccard_col(coords_params):
	"""Test calculating one set against a collection of sets."""

	k, sigs = coords_params

	for i, coords1 in enumerate(sigs):

		# Get scores against all others
		scores = jaccard_coords_col(coords1, sigs.values, sigs.bounds)

		# Check against single coords
		for j, coords2 in enumerate(sigs):

			if len(coords1) == 0 and len(coords2) == 0:
				continue

			assert scores[j] == jaccard_coords(coords1, coords2)
