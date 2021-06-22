"""Utility code for reading/writing data files."""

import os

import numpy as np


COMPRESSED_OPENERS = {None: open}


def _compressed_opener(compression):
	"""Decorator to register opener functions for compression types."""
	def decorator(func):
		COMPRESSED_OPENERS[compression] = func
		return func
	return decorator


@_compressed_opener('gzip')
def _open_gzip(path, mode, **kwargs):
	"""Opener for gzip-compressed files."""
	import gzip

	if mode is None:
		mode = 'rt'

	# gzip defaults to binary mode, change to text instead of not specified
	if mode[-1] not in 'tb':
		mode += 't'

	return gzip.open(path, mode=mode, **kwargs)


def open_compressed(compression, path, mode=None, **kwargs):
	"""Open a file with compression method specified by a string.

	Parameters
	----------
	compression : str
		Compression method. None is no compression. Keys of :data:`COMPRESSED_OPENERS` are the
		allowed values.
	path
		Path of file to open. May be string or path-like object.
	mode : str
		Mode to open file in - same as in :func:`open`.
	\\**kwargs
		Additional text-specific keyword arguments identical to the following :func:`open`
		arguments: ``encoding``, ``errors``, and ``newlines``.

	Returns
	-------
	io.BufferedIOBase
		Open file object.
	"""

	try:
		opener = COMPRESSED_OPENERS[compression]

	except KeyError:
		raise ValueError(f'Unknown compression type {compression!r}') from None

	return opener(os.fsdecode(path), mode=mode, **kwargs)


class ClosingIterator:
	"""Wraps an iterator which reads from a stream, closes the stream when finished.

	Used to wrap return values from functions which do some sort of lazy IO
	operation (specifically :func:`Bio.SeqIO.parse`) and return an iterator
	which reads from a stream every time ``next()`` is called on it. The object
	is an iterator itself, but will close the stream automatically when it
	finishes. May also be used as a context manager which closes the stream
	on exit.

	Attributes
	----------
	fobj

		The underlying file-like object or stream which the instance is
		responsible for closing

	iterator

		The iterator which the instance wraps.

	closed

		Read-only boolean property, mirrors the same attribute of :attr:`fobj`.

	Parameters
	----------
	iterable
		Iterable to iterate over. The :attr:`iterator` attribute will be obtained from calling
		:func:`iter` on this.
	fobj
		File-like object to close when the iterator finishes, context is exited or the :meth:`close`
		method is called.
	"""

	def __init__(self, iterable, fobj):
		self.iterator = iter(iterable)
		self.fobj = fobj

	def __iter__(self):
		return self

	def __next__(self):
		try:
			return next(self.iterator)

		except StopIteration:
			# Close when iterator runs out
			self.close()
			raise

	def close(self):
		"""Close the stream.

		Just calls the ``close`` method on :attr:`fobj`.
		"""
		self.fobj.close()

	@property
	def closed(self):
		return self.fobj.closed

	def __enter__(self):
		return self

	def __exit__(self, *args):
		self.close()


def read_npy(fobj, dtype, shape):
	"""Read a numpy array from raw data in a stream.

	Data expected to be in C-order.

	Parameters
	----------
	fobj
		Readable file-like object in binary mode.
	dtype : numpy.dtype
		Data type to read.
	shape
		Shape of array to read.

	Returns
	-------
	numpy.ndarray
	"""
	array = np.zeros(shape, dtype=dtype)
	fobj.readinto(array.data)
	return array


def write_npy(fobj, array):
	"""Write raw numpy array data to a stream.

	Warning! Data will be written in same order as the array, but
	:func:`read_npy` expects the data to be in c-order.

	Parameters
	----------
	fobj
		Writeable file-like object in binary mode.
	array : numpy.ndarray
		Array to write.
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

	Attributes
	----------
	_dtype : numpy.dtype
		Numpy data type defining structure.
	_data : numpy.ndarray
		Binary data. This is a zero-dimensional array` with data
		type :attr:`_dtype`, because Numpy does not support creating scalars of
		structured data types.
	_fields
		Mapping from field names to ``(dtype, offset)`` tuples.
	_names
		Field names in same order as physical layout in struct.
	_size : int
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
				raise KeyError(f'No field with name {index!r}') from None

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
				raise KeyError(f'No field with name {index!r}')

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

		Returns
		-------
		.NamedStruct
		"""
		return NamedStruct(self._dtype, self._data.copy())

	def _tobytes(self):
		"""Get the byte representation of the structure's data.

		Returns
		-------
		bytes
		"""
		return self._data.tobytes()

	def _astuple(self):
		"""Get field values as tuples.

		Sub-structs are converted recursively as well.

		Returns
		-------
		tuple
		"""
		return self._data[()]

	@classmethod
	def _fromfile(cls, dtype, fobj):
		"""Read struct from a stream.

		Parameters
		----------
		dtype : numpy.dtype
			Numpy data type for struct.
		fobj
			Readable file-like object in binary mode.

		Returns
		-------
		.NamedStruct
		"""
		struct = NamedStruct(dtype)
		fobj.readinto(struct._data.data)
		return struct

	@classmethod
	def _frombuffer(cls, dtype, buffer_, offset=0, copy=True):
		"""Create struct from raw data in a binary buffer.

		Parameters
		----------
		dtype : numpy.dtype
			Numpy data type for struct.
		buffer_
			Object implementing buffer protocol.
		offset : int
			Offset of beginning of struct data.
		copy : bool
			If True, use a copy of the array's data.

		Returns
		-------
		.NamedStruct
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
			' '.join(f'{name}={self[name]!r}' for name in self._names)
		)
