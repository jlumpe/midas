"""Test midas.signatures.SignatureArray."""

import pytest
import numpy as np

from midas.signatures import SignatureArray
from midas.test import make_signatures


@pytest.fixture(params=[None, 'i8', 'u4'])
def sigarray(request):
	np.random.seed(0)
	return make_signatures(8, 100, request.param)


@pytest.fixture()
def sigs(sigarray):
	"""Numpy array equivalent to `sigarray`."""
	return np.asarray(sigarray, dtype=object)


def sigarray_eq(a1, a2):
	"""Check two SignatureArrays or other sequences of signatures are equal."""
	return len(a1) == len(a2) and all(map(np.array_equal, a1, a2))


def test_basic(sigarray, sigs):
	"""Test basic functionality outside of __getitem__()."""

	# Check len
	assert len(sigarray) == len(sigs)

	# Check sizeof() method
	for i, sig in enumerate(sigs):
		assert sigarray.sizeof(i) == len(sig)
		assert sigarray.sizeof(np.int64(i)) == len(sig)

	with pytest.raises(IndexError):
		sigarray.sizeof(len(sigarray))

	# Assignment not supported
	with pytest.raises(TypeError):
		sigarray[0] = 0


def test_getitem_single(sigarray, sigs):
	n = len(sigarray)

	for i in range(n):
		sig = sigarray[i]
		assert np.array_equal(sig, sigs[i])
		assert sig.base is sigarray.values  # Check is view


def check_subseq(sigarray, sigs, index):
	"""Check result of indexing which results in a subsequence."""
	result = sigarray[index]
	assert isinstance(result, SignatureArray)
	assert result.values.dtype == sigarray.values.dtype
	assert sigarray_eq(result, sigs[index])
	return result


def test_getitem_slice(sigarray, sigs):
	slices = [
		slice(10, 30),
		slice(None),
		slice(None, None, 1),
		slice(None, None, 2),
	]

	for s in slices:
		result = check_subseq(sigarray, sigs, s)

		# Check values array is view for contiguous slices
		if s.step is None or s.step == 1:
			assert result.values.base is sigarray.values


def test_getitem_int_array(sigarray, sigs):
	indices = [
		[0, 10, 20, 30, 20, -10],
		[],
	]
	for index in indices:
		check_subseq(sigarray, sigs, index)


def test_getitem_bool_array(sigarray, sigs):
	index = np.arange(len(sigarray)) % 3 == 0
	check_subseq(sigarray, sigs, index)


def test_uninitialized(sigs):
	"""Test creating with uninitialized() class method."""

	lengths = list(map(len, sigs))
	sigarray = SignatureArray.uninitialized(lengths)
	assert len(sigarray) == len(sigs)

	for i in range(len(sigarray)):
		assert sigarray.sizeof(i) == len(sigs[i])


def test_construct_from_signaturearray(sigarray):
	"""Test construction from another SignatureArray."""
	sa2 = SignatureArray(sigarray)
	assert sa2.values.dtype == sigarray.values.dtype
	assert sigarray_eq(sa2, sigarray)

	for dtype in ['i8', 'u4']:
		sa2 = SignatureArray(sigarray, dtype=dtype)
		assert sa2.values.dtype == np.dtype(dtype)
		assert sigarray_eq(sa2, sigarray)


def test_empty():
	"""Really an edge case, but test it anyways."""

	sigarray = SignatureArray([])
	assert len(sigarray) == 0

	with pytest.raises(IndexError):
		sigarray[0]
