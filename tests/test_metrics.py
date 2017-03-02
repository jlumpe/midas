"""Test distance metric calcuations."""

import pickle

import pytest
import numpy as np

from midas.cython.metrics import jaccard_coords, jaccard_coords_col
from midas.kmers import KmerCoordsCollection, vec_to_coords, coords_to_vec


def slow_jaccard_coords(coords1, coords2, idx_len):
	"""Get jaccard score of two k-mer coordiate arrays the slow but correct way."""

	vec1 = coords_to_vec(coords1, idx_len)
	vec2 = coords_to_vec(coords2, idx_len)

	intersection = (vec1 & vec2).sum()
	union = (vec1 | vec2).sum()

	return intersection / union


def make_coords_col(k, nsets):
	"""Make artificial k-mer coord sets.

	:rtype: KmerCoordsCollection
	"""

	random = np.random.RandomState(seed=0)

	idx_len = 4 ** k

	coords_list = []

	# Add empty and full sets as edge cases
	coords_list.append(np.arange(0))
	coords_list.append(np.arange(idx_len))

	# Use a core set of k-mers so that we get some overlap
	core_prob = max(0.01, 20 / idx_len)
	core_vec = random.rand(idx_len) < core_prob

	for i in range(nsets - 2):

		keep_core = random.rand(idx_len) < (random.rand() ** 2)

		new_vec = random.rand(idx_len) < random.uniform(core_prob * .1, core_prob)

		vec = (keep_core & core_vec) | new_vec
		coords_list.append(vec_to_coords(vec))

	return KmerCoordsCollection.from_coords_seq(coords_list)


@pytest.fixture(scope='module')
def load_test_coords_col(fixture_file):
	"""Function to load k-mer coordinates from file."""

	def load_test_coords_col_func():

		with fixture_file('kmer_coords/coords.pickle').open('rb') as fobj:
			coords_list = pickle.load(fobj)

		return KmerCoordsCollection.from_coords_seq(coords_list)

	return load_test_coords_col_func


@pytest.fixture(params=[4, 8, 9, None])
def coords_params(request, load_test_coords_col):
	"""Tuple of (k, kmercoordscollection) to test on."""

	if request.param is None:
		# Load coords from file
		k = 11
		kcol = load_test_coords_col()

	else:
		# Create coords
		k = request.param
		kcol = make_coords_col(k, 25)

	return k, kcol


def test_jaccard_single(coords_params):
	"""Test calculating single scores at a time."""

	k, kcol = coords_params
	idx_len = 4 ** k

	# Iterate over all pairs
	for i, coords1 in enumerate(kcol):

		for j in range(i, len(kcol)):

			coords2 = kcol[j]

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

	k, kcol = coords_params
	idx_len = 4 ** k

	for i, coords1 in enumerate(kcol):

		# Get scores against all others
		scores = jaccard_coords_col(coords1, kcol.coords_array, kcol.bounds)

		# Check against single coords
		for j, coords2 in enumerate(kcol):

			if len(coords1) == 0 and len(coords2) == 0:
				continue

			assert scores[j] == jaccard_coords(coords1, coords2)
