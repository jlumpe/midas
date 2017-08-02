"""Test midas.io.signaturefile"""

import io

import pytest
import numpy as np

from midas.io.signaturefile import SignatureFile
from midas.kmers import SignatureArray


@pytest.fixture(scope='module')
def sigarray():
	"""Array of random signatures."""

	random = np.random.RandomState(0)

	lengths = random.randint(2000, 10000, size=20)
	choices = np.arange(4 ** 11)

	return SignatureArray.from_signatures([
		np.sort(random.choice(choices, replace=False, size=len_))
		for len_ in lengths
	])


@pytest.fixture(params=[None, 'strings', 'ints'])
def ids(sigarray, request):
	"""None or list of strings or ints."""
	if request.param == 'strings':
		return ['id' + str(i) for i in range(len(sigarray))]

	elif request.param == 'ints':
		return list(range(len(sigarray)))


@pytest.fixture(params=[False, True])
def metadata(request):
	"""None or a dict of fake metadata."""
	if request.param:
		return dict(foo=1, bar='z', baz=[5, 4, 1], bleh=dict(a=1, b=None))
	else:
		return None


@pytest.fixture
def signature_buffer(sigarray, ids, metadata):
	"""Binary buffer with signature file written to it."""

	buf = io.BytesIO()
	SignatureFile.write(buf, sigarray, ids=ids, metadata=metadata)
	buf.seek(0)

	return buf


@pytest.fixture
def sigfile(signature_buffer):
	"""SignatureFile read from buffer."""
	return SignatureFile(signature_buffer)


def test_attrs(sigarray, sigfile):
	"""Test basic attributes set when file is first opened."""

	assert sigfile.count == len(sigarray)
	assert np.array_equal(sigfile.lengths, list(map(len, sigarray)))
	assert sigfile.nelems == sum(map(len, sigarray))
	assert sigfile.dtype == sigarray.values.dtype


def test_ids(sigarray, sigfile, ids):
	"""Test the ids attribute matches what was written."""

	if ids is None:
		assert sigfile.ids is None

	else:
		assert np.array_equal(ids, sigfile.ids)


def test_metadata(sigarray, sigfile, metadata):
	"""Test the file metadata matches what was written."""

	if metadata is None:
		assert not sigfile.has_metadata
		assert sigfile.get_metadata() is None
	else:
		assert sigfile.has_metadata
		assert sigfile.get_metadata() == metadata


def test_read_all(sigarray, sigfile):
	"""Test reading all signatures into an array."""

	sigarray2 = sigfile.get_array()
	assert np.array_equal(sigarray.values, sigarray2.values)
	assert np.array_equal(sigarray.bounds, sigarray2.bounds)


def test_read_subset(sigarray, sigfile):
	"""Test reading a subset of signatures into an array."""

	n = len(sigarray)
	random = np.random.RandomState()
	indices = random.choice(n, n // 2)

	subarray = sigfile.get_array(indices)

	assert len(subarray) == len(indices)

	for idx, sig in zip(indices, subarray):
		assert np.array_equal(sig, sigarray[idx])


def test_iter(sigarray, sigfile):
	"""Test iterating over signatures."""

	for sig1, sig2 in zip(sigarray, sigfile.iter_signatures()):
		assert np.array_equal(sig1, sig2)
