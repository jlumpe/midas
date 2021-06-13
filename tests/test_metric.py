"""Test distance metric calculations."""

import pickle

import pytest
import numpy as np

from midas.metric import jaccard_sparse, jaccarddist_sparse, jaccard_bits, \
	jaccard_generic, jaccard_sparse_array
from midas.kmers import SignatureArray, sparse_to_dense
from midas.test import make_signatures


@pytest.fixture(scope='module')
def load_test_coords_col(test_data):
	"""Function to load k-mer coordinates from file."""

	def load_test_coords_col_func():
		with open(test_data / 'kmer_coords/coords.pickle', 'rb') as fobj:
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

	# Iterate over all pairs
	for i, coords1 in enumerate(sigs):

		for j in range(i, len(sigs)):

			coords2 = sigs[j]

			# Calc score
			score = jaccard_sparse(coords1, coords2)

			# Check range
			assert 0 <= score <= 1

			# Check vs slow version
			assert np.isclose(score, jaccard_generic(coords1, coords2))

			# Check with arguments swapped
			assert score == jaccard_sparse(coords2, coords1)

			# Check score vs. self is one (unless empty)
			if i == j and len(coords1) > 0:
				assert score == 1


def test_different_dtypes():

	dtypes = ['i2', 'u2', 'i4', 'u4', 'i8', 'u8']
	sigs_by_k = {}

	# Test all pairs of dtypes
	for i, dt1 in enumerate(dtypes):
		for j in range(i + 1, len(dtypes)):
			dt2 = dtypes[j]

			# Bit length of largest positive integer both dtypes can store
			nbits = min(np.iinfo(dt).max.bit_length() for dt in [dt1, dt2])

			# Try for k == 8, but not larger than will fit in dtype
			k = min(nbits // 4, 8)

			# Get signatures for this value of k, or create them
			try:
				sigs = sigs_by_k[k]
			except KeyError:
				# Create using largest dtype
				sigs = sigs_by_k[k] = make_signatures(k, 5, 'u8')

			# Convert to each dtype, making sure there is no overflow
			try:
				# Tell numpy to raise error on overflow
				old_err = np.seterr(over='raise')

				sigs1 = SignatureArray.from_signatures(sigs, dtype=dt1)
				sigs2 = SignatureArray.from_signatures(sigs, dtype=dt2)

			finally:
				np.seterr(**old_err)

			# Iterate over all pairs of signatures to test individually
			for k in range(len(sigs)):
				for l in range(k + 1, len(sigs)):

					expected = jaccard_sparse(sigs[k], sigs[l])

					assert jaccard_sparse(sigs1[k], sigs2[l]) == expected
					assert jaccard_sparse(sigs2[k], sigs1[l]) == expected

			# Test first of each against all of the other
			expected_all = jaccard_sparse_array(sigs[0], sigs)

			assert np.array_equal(
				jaccard_sparse_array(sigs1[0], sigs2),
				expected_all
			)
			assert np.array_equal(
				jaccard_sparse_array(sigs2[0], sigs1),
				expected_all
			)


def test_jaccard_col(coords_params):
	"""Test calculating one set against a collection of sets."""

	k, sigs = coords_params

	for i, coords1 in enumerate(sigs):

		# Get scores against all others
		scores = jaccard_sparse_array(coords1, sigs)

		# Check against single coords
		for j, coords2 in enumerate(sigs):

			if len(coords1) == 0 and len(coords2) == 0:
				continue

			assert scores[j] == jaccard_sparse(coords1, coords2)
