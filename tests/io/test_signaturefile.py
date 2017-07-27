"""Test midas.signaturefile"""

import io

import pytest
import numpy as np

from midas.io.signaturefile import SignatureFile
from midas.kmers import SignatureArray


@pytest.fixture
def signatures():
	"""List of random signatures."""

	random = np.random.RandomState()

	lengths = random.randint(2000, 10000, size=10)
	choices = np.arange(4 ** 11)

	return [
		np.sort(random.choice(choices, replace=False, size=len_))
		for len_ in lengths
	]


@pytest.fixture(params=[None, 'strings', 'ints'])
def ids(signatures, request):
	if request.param == 'strings':
		return ['id' + str(i) for i in range(len(signatures))]

	elif request.param == 'ints':
		return list(range(len(signatures)))


@pytest.fixture(params=[False, True])
def metadata(request):
	if request.param:
		return dict(foo=1, bar='z', baz=[5, 4, 1], bleh=dict(a=1, b=None))
	else:
		return None


def test_readwrite(signatures, ids, metadata):

	sigarray = SignatureArray.from_signatures(signatures)

	buf = io.BytesIO()

	# Write
	SignatureFile.write_sigarray(buf, sigarray, ids=ids, metadata=metadata)
	buf.seek(0)

	# Read
	sigfile = SignatureFile(buf)

	assert sigfile.count == len(signatures)
	assert np.array_equal(sigfile.lengths, list(map(len, signatures)))

	# Test reading whole array
	sigarray2 = sigfile.get_array()
	assert np.array_equal(sigarray.values, sigarray2.values)
	assert np.array_equal(sigarray.bounds, sigarray2.bounds)

	# Test iterating over signatures
	for sig1, sig2 in zip(sigarray, sigfile.iter_signatures()):
		assert np.array_equal(sig1, sig2)

	# Check IDs match
	if ids is None:
		assert sigfile.ids is None

	else:
		assert np.array_equal(ids, sigfile.ids)

	# Check metadata
	assert sigfile.get_metadata() == metadata
