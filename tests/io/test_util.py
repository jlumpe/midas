"""Test midas.io.util"""

import io

import pytest
import numpy as np

from midas.io import util


class TestOpenCompressed:
	"""Test midas.io.util.open_compressed."""

	@pytest.fixture(scope='class')
	def text_data(self):
		"""Random printable characters encoded as ASCII."""
		random = np.random.RandomState()
		return random.randint(32, 128, size=1000, dtype='b').tobytes()

	@pytest.fixture(scope='class', params=list(util.COMPRESSED_OPENERS))
	def compression(self, request):
		"""Compression method string."""
		return request.param

	@pytest.fixture(params=[False, True])
	def text_file(self, request, text_data, compression, tmpdir):

		file = tmpdir.join('chars.txt')

		if request.param:
			mode = 'wt'
			to_write = text_data.decode('ascii')

		else:
			mode = 'wb'
			to_write = text_data

		with util.open_compressed(compression, file.strpath, mode) as fobj:
			fobj.write(to_write)

		return file

	@pytest.mark.parametrize('text_mode', [False, True])
	def test_read(self, text_mode, text_data, compression, text_file):

		if text_mode:
			mode = 'rt'
			expected = text_data.decode('ascii')
		else:
			mode = 'rb'
			expected = text_data

		with util.open_compressed(compression, text_file.strpath, mode) as fobj:
			contents = fobj.read()

		assert contents == expected


class TestClosingIterator:
	"""Test the ClosingIterator class."""

	# Number of lines in text buffer
	NLINES = 100

	@pytest.fixture
	def fobj(self):
		"""Text buffer with a number in each line."""

		from io import StringIO

		# Write some lines to a buffer
		buf = StringIO()

		for i in range(self.NLINES):
			buf.write('{}\n'.format(i))

		buf.seek(0)

		return buf

	@pytest.fixture
	def iterator(self, fobj):
		"""Instance of ClosingIterator using fobj."""

		# Iterable which reads from fobj on every iteration."""
		iterable = (int(line.strip()) for line in fobj)

		# Create ClosingIterator object
		return util.ClosingIterator(iterable, fobj)

	def test_close_on_finish(self, iterator, fobj):
		"""Check that the stream gets closed when the iterator runs out."""

		assert not fobj.closed and not iterator.closed

		for val in iterator:
			assert not fobj.closed and not iterator.closed

		assert fobj.closed and iterator.closed

	def test_context(self, iterator, fobj):

		with iterator as rval:

			# Check __exit__() returns itself
			assert rval is iterator

			assert not fobj.closed and not iterator.closed

		# Should close after context exited
		assert fobj.closed and iterator.closed

	def test_close_method(self, iterator, fobj):
		"""Test the close() method."""

		assert not fobj.closed and not iterator.closed

		iterator.close()

		assert fobj.closed and iterator.closed


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
	util.write_npy(buf, array)
	assert buf.tell() == len(data)

	# Read
	buf.seek(0)
	array2 = util.read_npy(buf, dtype_str, shape)

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
		return util.NamedStruct._frombuffer(self.STRUCT_DTYPE, struct_data)

	def field_eq(self, val1, val2):
		"""Check if two NamedStruct field values are equal.

		Values may be scalars, arrays, or NamedStructs.
		"""

		if isinstance(val1, util.NamedStruct):
			val1 = val1._astuple()

		if isinstance(val2, util.NamedStruct):
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

		assert struct == util.NamedStruct._frombuffer(self.STRUCT_DTYPE, struct_data)

		# Same dtype, zero data
		assert struct != util.NamedStruct(self.STRUCT_DTYPE)

		# Same data, different dtype
		assert struct != util.NamedStruct._frombuffer(
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
				assert isinstance(value, util.NamedStruct)

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

			if isinstance(prev_val, util.NamedStruct):
				# Sub-struct - try first assigning as NamedStruct

				prev_substruct = prev_val._copy()
				prev_val = prev_substruct._astuple()

				# Zero version of sub-struct (shouldn't already be equal)
				zero_substruct = util.NamedStruct(prev_substruct._dtype)
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
