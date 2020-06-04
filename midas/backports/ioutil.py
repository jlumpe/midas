"""Utility code for reading/writing data files.

Backported from midas.io.util."""

import numpy as np


def read_npy(fobj, dtype, shape):
	"""Read a numpy array from raw data in a stream.

	Data expected to be in C-order.

	:param fobj: Readable file-like object in binary mode.
	:type dtype: numpy.dtype
	:param shape: Shape of array to read.

	:rtype: numpy.ndarray
	"""
	array = np.zeros(shape, dtype=dtype)
	fobj.readinto(array.data)
	return array


def write_npy(fobj, array):
	"""Write raw numpy array data to a stream.

	Warning! Data will be written in same order as the array, but
	:func:`read_npy` expects the data to be in c-order.

	:param fobj: Writeable file-like object in binary mode.
	:param array: Array to write.
	:type array: numpy.ndarray
	"""

	if array.flags.c_contiguous:
		data = array.data
	else:
		data = array.tobytes()

	written = fobj.write(data)

	if written != array.nbytes:
		raise IOError('Could not write all bytes of array')


class NamedStruct:
	"""Stores packed binary data with named fields, similar to a C struct.

	Fields are defined by a numpy structured data type. Provides read and write
	access to field values using attribute and item access syntax. May be
	recursive - fields can themselves be structs and accessing their values will
	return :class:`.NamedStruct` s as well.

	Also acts as a sequence of field values, and item access is supported by
	numeric index.

	Attributes which are not fields are prefixed by an underscore, but are
	still part of the public API.

	.. attribute:: _dtype

		:class:`numpy.dtype` defining structure.

	.. attribute:: _data

		Binary data. This is a zero-dimensional :class:`numpy.ndarray` with data
		type :attr:`_dtype`, because Numpy does not support creating scalars of
		structured data types.

	.. attribute:: _fields

		Mapping from field names to ``(dtype, offset)`` tuples.

	.. attribute:: _names

		Field names in same order as physical layout in struct.

	.. attribute:: _size

		Struct size in bytes.
	"""

	__slots__ = ['_dtype', '_data']

	def __init__(self, dtype, data=None):
		self._dtype = np.dtype(dtype)

		if self._dtype.fields is None:
			raise ValueError('Not a structured dtype')

		for name in self._dtype.fields:
			if name.startswith('_'):
				raise ValueError('Field names cannot start with underscore')

		if data is None:
			self._data = np.zeros((), dtype=dtype)
		else:
			self._data = np.asarray(data, dtype=dtype).reshape(())

	@property
	def _fields(self):
		return self._dtype.fields

	@property
	def _names(self):
		return self._dtype.names

	@property
	def _size(self):
		return self._dtype.itemsize

	def __eq__(self, other):
		return isinstance(other, NamedStruct) and \
			self._dtype == other._dtype and \
			np.array_equal(self._data, other._data)

	def __len__(self):
		return len(self._dtype)

	def __iter__(self):
		for name in self._names:
			yield self[name]

	def __getitem__(self, index):

		if isinstance(index, int):
			field_dtype = self._dtype[index]
			field_array = self._data[self._names[index]]

		elif isinstance(index, str):
			try:
				field_dtype, _ = self._dtype.fields[index]
			except KeyError:
				raise KeyError('No field with name {!r}'.format(index)) from None

			field_array = self._data[index]

		else:
			raise TypeError(type(index))

		if field_dtype.fields is not None:
			return NamedStruct(field_dtype, field_array)
		else:
			return field_array[()]

	def __setitem__(self, index, value):

		# Get name of field to set
		if isinstance(index, str):
			name = index

			if name not in self._dtype.fields:
				raise KeyError('No field with name {!r}'.format(index))

		else:
			name = self._names[index]

		if isinstance(value, NamedStruct):
			value = value._astuple()

		# Set
		self._data[name][()] = value

	def __delitem__(self, index):
		raise TypeError('Field deletion not supported')

	def __getattr__(self, name):
		if not name.startswith('_'):
			# Try getting field value by name
			return self[name]

		return object.__getattribute__(self, name)

	def __setattr__(self, name, value):
		if name.startswith('_'):
			# If prefixed by underscore, set as regular attribute
			object.__setattr__(self, name, value)

		else:
			# No underscore means setting a field value
			self[name] = value

	def __delattr__(self, name):
		if name.startswith('_'):
			# If prefixed by underscore, set as regular attribute
			object.__delattr__(self, name)

		else:
			raise TypeError('Field deletion not supported')

	def _copy(self):
		"""Get a copy of the struct with its own data.

		:rtype: .NamedStruct
		"""
		return NamedStruct(self._dtype, self._data.copy())

	def _tobytes(self):
		"""Get the byte representation of the structure's data.

		:rtype: bytes
		"""
		return self._data.tobytes()

	def _astuple(self):
		"""Get field values as tuples.

		Sub-structs are converted recursively as well.

		:rtype: tuple
		"""
		return self._data[()]

	@classmethod
	def _fromfile(cls, dtype, fobj):
		"""Read struct from a stream.

		:param dtype: Numpy data type for struct.
		:type dtype: numpy.dtype
		:param fobj: Readable file-like object in binary mode.
		:rtype: .NamedStruct
		"""
		struct = NamedStruct(dtype)
		fobj.readinto(struct._data.data)
		return struct

	@classmethod
	def _frombuffer(cls, dtype, buffer_, offset=0, copy=True):
		"""Create struct from raw data in a binary buffer.

		:param dtype: Numpy data type for struct.
		:type dtype: numpy.dtype
		:param bytes_: Object implementing buffer protocol.
		:param int offset: Offset of beginning of struct data.
		:param bool copy: If True, use a copy of the array's data.

		:rtype: .NamedStruct
		"""
		array = np.frombuffer(buffer_, dtype=dtype, count=1, offset=offset)

		if copy:
			array = array.copy()

		return NamedStruct(dtype, array)

	def __dir__(self):
		attrs = []

		attrs.extend(self.__slots__)
		attrs.extend(type(self).__dict__)
		attrs.extend(self._dtype.names)

		return attrs

	def __repr__(self):
		return '<{} {}>'.format(
			type(self).__name__,
			' '.join('{}={!r}'.format(name, self[name]) for name in self._names)
		)
