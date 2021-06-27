"""Test midas.signatures.signaturefile."""

import io

import pytest
import numpy as np

from midas.signatures.signaturefile import SignatureFile, read_npy, write_npy, NamedStruct
from midas.signatures import SignatureArray


@pytest.mark.parametrize('shape', [(), (1,), (100,), (5, 3, 4)])
@pytest.mark.parametrize('dtype_str', ['i2', 'u8', 'f4', 'S10', '(2,3)i4, S4'])
def test_read_write_npy(shape, dtype_str):
	"""Test read_npy and write_npy."""

	dtype = np.dtype(dtype_str)
	count = np.prod(shape, dtype=int)

	# Array of random data
	random = np.random.RandomState(0)
	data = random.bytes(count * dtype.itemsize)

	array = np.frombuffer(data, dtype=dtype, count=count).reshape(shape)

	# Buffer to read to / write from
	buf = io.BytesIO()

	# Write
	write_npy(buf, array)
	assert buf.tell() == len(data)

	# Read
	buf.seek(0)
	array2 = read_npy(buf, dtype_str, shape)

	assert np.array_equal(array, array2)


class TestNamedStruct:
	"""Test NamedStruct class."""

	# Data type for struct to test - stolen from SignatureFile
	STRUCT_DTYPE = np.dtype([
		('magic_number', 'S4'),
		('version', 'S4'),
		('count', 'u8'),
		('dtype', 'S2'),
		('offsets', [
			('lengths', '2u8'),
			('metadata', '2u8'),
			('ids', '2u8'),
			('data', '2u8'),
		]),
		('multid_array', '(3,4,5)i4'),
	])

	@pytest.fixture
	def struct_data(self):
		"""Random bytes to create struct with."""
		random = np.random.RandomState(0)
		return random.bytes(self.STRUCT_DTYPE.itemsize)

	@pytest.fixture
	def struct(self, struct_data):
		"""Struct from data."""
		return NamedStruct._frombuffer(self.STRUCT_DTYPE, struct_data)

	def field_eq(self, val1, val2):
		"""Check if two NamedStruct field values are equal.

		Values may be scalars, arrays, or NamedStructs.
		"""

		if isinstance(val1, NamedStruct):
			val1 = val1._astuple()

		if isinstance(val2, NamedStruct):
			val2 = val2._astuple()

		if isinstance(val1, tuple):

			if isinstance(val2, tuple):
				return val1 == val2

			else:
				return False

		else:
			return np.array_equiv(val1, val2)

	def test_basic(self, struct, struct_data):
		"""Test basic attributes."""

		assert struct._dtype == self.STRUCT_DTYPE
		assert struct._fields == self.STRUCT_DTYPE.fields
		assert struct._names == self.STRUCT_DTYPE.names
		assert struct._size == self.STRUCT_DTYPE.itemsize

		assert struct._tobytes() == struct_data
		assert struct._astuple() == struct._data[()]

	def test_eq(self, struct):
		"""Test equality."""
		struct_data = struct._tobytes()

		assert struct == NamedStruct._frombuffer(self.STRUCT_DTYPE, struct_data)

		# Same dtype, zero data
		assert struct != NamedStruct(self.STRUCT_DTYPE)

		# Same data, different dtype
		assert struct != NamedStruct._frombuffer(
			[('foo', 'b', struct._size)],
			struct_data,
		)

	def test_read_fields(self, struct):
		"""Test getting field values."""

		# Check each existing field
		for name in struct._names:

			field_dtype, _ = struct._fields[name]
			value = struct[name]

			# Check attribute access matches item access
			assert np.array_equiv(getattr(struct, name), value)

			# Check value matches value of field pulled from array
			array_value = struct._data[name]

			if field_dtype.fields is not None:
				# Sub-struct
				assert isinstance(value, NamedStruct)

			assert self.field_eq(array_value, value)

		# Check nonexistant fields
		with pytest.raises(KeyError):
			struct['notafield']

		with pytest.raises(KeyError):
			struct.notafield

	def test_write_fields(self, struct):
		"""Test setting field values."""

		# Test assignment to each field
		for name in struct._names:

			prev_val = struct[name]

			if isinstance(prev_val, NamedStruct):
				# Sub-struct - try first assigning as NamedStruct

				prev_substruct = prev_val._copy()
				prev_val = prev_substruct._astuple()

				# Zero version of sub-struct (shouldn't already be equal)
				zero_substruct = NamedStruct(prev_substruct._dtype)
				assert zero_substruct != prev_substruct

				# Set to zero using item assignment syntax
				struct[name] = zero_substruct
				assert struct[name] == zero_substruct

				# Set back to original using attribute assignment syntax
				setattr(struct, name, prev_substruct)
				assert struct[name] == prev_substruct

			# Zero version of value (shouldn't already be equal)
			zeros = np.zeros_like(prev_val)[()]
			assert not self.field_eq(zeros, prev_val)

			# Set to zero using item assignment syntax
			struct[name] = zeros
			assert self.field_eq(struct._data[name][()], zeros)

			# Set back to original using attribute assignment syntax
			setattr(struct, name, prev_val)
			assert self.field_eq(struct._data[name][()], prev_val)

		# Check nonexistant fields
		with pytest.raises(KeyError):
			struct['notafield'] = 0

		with pytest.raises(KeyError):
			setattr(struct, 'notafield', 0)

	def test_sequence(self, struct):
		"""Test sequence protocol."""

		assert len(struct) == len(self.STRUCT_DTYPE.fields)

		for i, value in enumerate(struct):
			assert self.field_eq(value, struct[i])

		assert i == len(struct) - 1


class ProgressChecker:
	"""Callback for get_array method that checks it is called correctly."""

	def __init__(self, total, chunksize=None):
		self.total = total
		self.chunksize = 1 if chunksize is None else chunksize
		self.last_value = 0

	def __call__(self, completed, total):
		"""The callback function."""
		assert total == self.total

		assert completed == min(self.last_value + self.chunksize, self.total)

		self.last_value = completed

	def check_called(self):
		"""To be called at end, check we went all the way through."""
		assert self.last_value == self.total


class TestSignatureFile:

	@pytest.fixture(scope='class')
	def sigarray(self):
		"""Array of random signatures."""

		random = np.random.RandomState(0)

		lengths = random.randint(2000, 10000, size=20)

		return SignatureArray([
			np.sort(random.choice(4 ** 11, replace=False, size=len_))
			for len_ in lengths
		])

	@pytest.fixture(params=[None, 'strings', 'ints'])
	def ids(self, sigarray, request):
		"""None or list of strings or ints."""
		if request.param == 'strings':
			return ['id' + str(i) for i in range(len(sigarray))]

		elif request.param == 'ints':
			return list(range(len(sigarray)))


	@pytest.fixture(params=[False, True])
	def metadata(self, request):
		"""None or a dict of fake metadata."""
		if request.param:
			return dict(foo=1, bar='z', baz=[5, 4, 1], bleh=dict(a=1, b=None))
		else:
			return None

	@pytest.fixture
	def signature_buffer(self, sigarray, ids, metadata):
		"""Binary buffer with signature file written to it."""

		buf = io.BytesIO()
		SignatureFile.write(buf, sigarray, ids=ids, metadata=metadata)
		buf.seek(0)

		return buf

	@pytest.fixture
	def sigfile(self, signature_buffer):
		"""SignatureFile read from buffer."""
		return SignatureFile(signature_buffer)

	def test_attrs(self, sigarray, sigfile):
		"""Test basic attributes set when file is first opened."""

		assert sigfile.count == len(sigarray)
		assert np.array_equal(sigfile.lengths, list(map(len, sigarray)))
		assert sigfile.nelems == sum(map(len, sigarray))
		assert sigfile.dtype == sigarray.values.dtype

	def test_ids(self, sigarray, sigfile, ids):
		"""Test the ids attribute matches what was written."""

		if ids is None:
			assert sigfile.ids is None

		else:
			assert np.array_equal(ids, sigfile.ids)

	def test_metadata(self, sigarray, sigfile, metadata):
		"""Test the file metadata matches what was written."""

		if metadata is None:
			assert not sigfile.has_metadata
			assert sigfile.get_metadata() is None
		else:
			assert sigfile.has_metadata
			assert sigfile.get_metadata() == metadata

	@pytest.mark.parametrize('chunksize', [None, 1, 3, 100])
	@pytest.mark.parametrize('progress', [False, True])
	def test_read_all(self, sigarray, sigfile, chunksize, progress):
		"""Test reading all signatures into an array."""

		# Progress callback
		if progress:
			callback = ProgressChecker(
				len(sigarray),
				len(sigarray) if chunksize is None else chunksize
			)
		else:
			callback = None

		# Read signatures
		sigarray2 = sigfile.get_array(chunksize=chunksize, progress=callback)

		assert np.array_equal(sigarray.values, sigarray2.values)
		assert np.array_equal(sigarray.bounds, sigarray2.bounds)

		# Check progress callback was called the correct number of times
		if progress:
			callback.check_called()

	@pytest.mark.parametrize('progress', [False, True])
	def test_read_subset(self, sigarray, sigfile, progress):
		"""Test reading a subset of signatures into an array."""

		# Random set of indices
		n = len(sigarray)
		random = np.random.RandomState()
		indices = random.choice(n, n // 2)

		# Progress callback
		callback = ProgressChecker(len(indices)) if progress else None

		# Read subset of signatures
		subarray = sigfile.get_array(indices, progress=callback)

		assert len(subarray) == len(indices)

		for idx, sig in zip(indices, subarray):
			assert np.array_equal(sig, sigarray[idx])

		# Check progress callback was called the correct number of times
		if progress:
			callback.check_called()

	def test_iter(self, sigarray, sigfile):
		"""Test iterating over signatures."""

		for sig1, sig2 in zip(sigarray, sigfile.iter_signatures()):
			assert np.array_equal(sig1, sig2)

	def test_context(self, sigfile):
		"""Test context manager methods."""

		with sigfile as obj:
			assert obj is sigfile
			assert not sigfile.fobj.closed

		assert sigfile.fobj.closed
