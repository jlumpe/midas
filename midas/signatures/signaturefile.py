"""Store sets of K-mer signatures on disk in binary format."""

import io
import json

import numpy as np

from midas.signatures import SignatureArray


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


class SignatureFile:
	"""Binary file for storing k-mer signatures.

	Constructor creates an object to read data from a stream. Use the
	:meth:`write` method to write data to a file.

	Acts as a context manager which closes its stream on exit.

	Parameters
	----------
	fobj
		Readable file-like object in binary mode.

	Attributes
	----------
	fobj
		Stream file is being read from.
	count : int
		Number of signatures in the file.
	dtype : numpy.dtype
		 Data type of signatures in the file.
	lengths : numpy.ndarray
		Length of each signature in the file.
	nelems : int
		Total number of elements of all signatures in the file.
	ids : numpy.ndarray
		IDs for signatures array of integers or string objects, or None if the file doesn't have IDs.
	has_metadata : bool
		Whether the file has metadata stored.
	closed : bool
		Whether :attr:`fobj` is closed.
	"""

	# Four-byte magic number to put at beginning of file
	_MAGIC_NUMBER = b'MSF\xFF'

	# Four bytes specifying current file format version
	_VERSION = b'1.00'

	# Numpy structured layout defining binary layout of offsets
	_OFFSETS_DTYPE = np.dtype([
		('lengths', '2i8'),
		('metadata', '2i8'),
		('ids', '2i8'),
		('data', '2i8'),
	])

	# Numpy structured layout defining binary layout of header
	_HEADER_DTYPE = np.dtype([
		('magic_number', ('S', len(_MAGIC_NUMBER))),
		('version', 'S4'),
		('count', 'i8'),
		('dtype', 'S2'),
		('offsets', _OFFSETS_DTYPE),
	])

	# Data type of lengths segment
	_LENGTHS_DTYPE = np.dtype('i4')

	# Default file extension
	DEFAULT_EXT = '.midas-signatures'

	def __init__(self, fobj):

		# Read and validate header
		self.fobj = fobj
		self._header = NamedStruct._fromfile(self._HEADER_DTYPE, fobj)

		self._validate_header(self._header)

		self.version = self._header.version.decode('ascii')
		self.count = self._header.count
		self.dtype = np.dtype(self._header.dtype)

		# Read IDs
		if self._header.offsets.ids[0] > 0:
			self.ids = self._read_ids()
		else:
			self.ids = None

		# Check if metadata present
		self.has_metadata = self._header.offsets.metadata[0] > 0

		# Read lengths
		fobj.seek(self._header.offsets.lengths[0])
		self.lengths = read_npy(fobj, self._LENGTHS_DTYPE, shape=self.count)
		self.lengths.flags.writeable = False

		self.nelems = self.lengths.sum()

	def __enter__(self):
		return self

	def __exit__(self, *args):
		self.fobj.close()

	def close(self):
		"""Close the underlying file stream."""
		self.fobj.close()

	@property
	def closed(self):
		return self.fobj.closed

	@classmethod
	def _validate_header(cls, header):
		"""Validate header data."""

		if header.magic_number != cls._MAGIC_NUMBER:
			raise OSError('File does not appear to be in the correct format')

		if header.version != cls._VERSION:
			raise OSError('Unexpected version identifier')

	def get_array(self, indices=None, progress=None, chunksize=None):
		"""Read signatures from file as a SignatureArray.

		Parameters
		----------
		indices
			List or array of indices of a subset of signatures to get. Array signatures will then
			correspond to the ids `signaturefile.ids[indices]``. If None will get all signatures.
		progress : function
			Progress callback that will be called periodically with the number of signatures read.
			Will be called with positional arguments ``(ncompleted, total)`` where ``ncompleted`` is
			the total number of signatures read so far and ``total`` is the total number to read.
			Will be called after every signature has been read if ``indices`` is not None or after
			every ``chunksize`` signatures otherwise.
		chunksize : int
			Number of signatures to read at a time, if ``indices`` is None (``indices`` and
			``chunksize`` should not both be given). Only useful with the ``progress`` argument. If
			None (default) will read all signatures in one shot.

		Returns
		-------
		midas.signatures.SignatureArray
		"""

		# File position at start of data
		data_start = self._header.offsets['data'][0]

		if indices is None:
			# Reading whole array in order

			if chunksize is None:
				chunksize = self.count

			# Array to read into
			array = SignatureArray.uninitialized(self.lengths, dtype=self.dtype)

			# Read in chunks
			self.fobj.seek(data_start)

			for begin in range(0, self.count, chunksize):
				end = min(begin + chunksize, self.count)

				nitems = self.lengths[begin:end].sum()
				values = read_npy(self.fobj, self.dtype, shape=nitems)

				array.values[array.bounds[begin]:array.bounds[end]] = values

				# Progress callback
				if progress is not None:
					progress(end, self.count)

			return array

		else:
			# Reading subset or permutation

			if chunksize is not None:
				raise TypeError('Cannot specify both chunksize and indices')

			# Use this dtype as it is what the Cython metric functions expect
			from midas._cython.metric import BOUNDS_DTYPE

			# Calculate file sub-array bounds from lengths
			bounds = np.zeros(self.count + 1, dtype=BOUNDS_DTYPE)
			bounds[1:] = np.cumsum(self.lengths)

			# Calculate bounds for indices only
			array = SignatureArray.uninitialized(self.lengths[indices], dtype=self.dtype)

			# Read one at a time, in order of position in file
			index_pairs_sorted = sorted((j, i) for i, j in enumerate(indices))

			for i, (file_idx, out_idx) in enumerate(index_pairs_sorted):

				self.fobj.seek(data_start + bounds[file_idx] * self.dtype.itemsize)
				signature = read_npy(self.fobj, self.dtype, self.lengths[file_idx])

				np.copyto(array[out_idx], signature)

				# Progress callback
				if progress is not None:
					progress(i + 1, len(indices))

			return array

	def iter_signatures(self):
		"""Iterate over signatures in the file.

		Returns
		-------
			Generator yielding :class:`numpy.ndarray`.
		"""

		self.fobj.seek(self._header.offsets['data'][0])

		for length in self.lengths:
			yield read_npy(self.fobj, self.dtype, shape=length)

	@classmethod
	def write(cls, fobj, signatures, *, dtype=None, ids=None, metadata=None):
		"""Write signatures to file.

		Parameters
		----------
		fobj
			Writeable file-like object in binary mode.
		signatures
			Sequence of signatures, as numpy arrays.
		dtype : numpy.dtype
			Numpy data type of signatures to write (should be unsigned integer). If None will be
			determined automatically by signatures.
		ids
			Sequence of IDs for the signatures, for example accession numbers of the genomes the
			signatures are from. Must be the same length as ``signatures`` and items must be all
			strings or integers.
		metadata
			Arbitrary metadata that will be included with file. Must be JSONable.
		"""

		count = len(signatures)

		if ids is not None and len(ids) != count:
			raise ValueError('Number of IDs does not match number of signatures')

		# Get data type for signature arrays
		if dtype is not None:
			dtype = np.dtype(dtype)

		elif isinstance(signatures, SignatureArray):
			dtype = signatures.values.dtype

		else:
			dtype = np.find_common_type(list(set(a.dtype for a in signatures)), [])

		# Create header data
		header = NamedStruct(cls._HEADER_DTYPE)

		header.magic_number = cls._MAGIC_NUMBER
		header.version = cls._VERSION
		header.count = count
		header.dtype = dtype.str[1:]

		# Write header data (will come back and fill in offsets later)
		fobj.write(header._tobytes())

		# Write lengths
		lengths = np.fromiter(
			map(len, signatures),
			dtype=cls._LENGTHS_DTYPE,
			count=count,
		)

		header.offsets.lengths[0] = fobj.tell()
		write_npy(fobj, lengths)
		header.offsets.lengths[1] = fobj.tell() - 1

		# Write metadata
		if metadata is not None:
			header.offsets.metadata = cls._write_metadata(fobj, metadata)

		# Write IDs
		if ids is not None:
			header.offsets.ids = cls._write_ids(fobj, ids)

		# Write data
		header.offsets.data[0] = fobj.tell()
		for signature in signatures:
			write_npy(fobj, signature.astype(dtype))
		header.offsets.data[1] = fobj.tell() - 1

		# Go back and write offsets
		_, offsets_offset = cls._HEADER_DTYPE.fields['offsets']
		fobj.seek(offsets_offset)
		fobj.write(header.offsets._tobytes())

	@classmethod
	def _write_metadata(cls, fobj, metadata):
		"""Write metadata to file.

		Currently only supports JSON format.

		Returns
		-------
		tuple
			Tuple of begin and end file offsets for the metadata section.
		"""
		begin = fobj.tell()

		# Support JSON format only for now
		fobj.write(b'j')

		wrapper = io.TextIOWrapper(fobj)
		json.dump(metadata, wrapper)
		wrapper.detach()

		end = fobj.tell() - 1
		return begin, end

	@classmethod
	def _write_ids(cls, fobj, ids):
		"""Write IDs to file.

		Returns
		-------
		tuple
			Tuple of begin and end file offsets for the IDs section.
		"""
		begin = fobj.tell()

		is_strings = all(isinstance(item, str) for item in ids)

		if is_strings:
			# All strings, write directly to file terminated by null chars

			fobj.write(b's')

			for string in ids:
				fobj.write(string.encode())
				fobj.write(b'\0')

		else:
			# Should be integer array
			ids = np.asarray(ids)

			if ids.dtype.kind not in 'iu':
				raise ValueError('IDs must be strings or integers')

			# Write format character, dtype string, then raw data
			fobj.write(b'i')
			fobj.write(ids.dtype.str[1:].encode('ascii'))
			write_npy(fobj, ids)

		end = fobj.tell() - 1
		return begin, end

	def get_metadata(self):
		"""Read metadata from the file, if it has any.

		Returns
		-------
			File's metadata object (typically dict), or None if the file has no metadata.
		"""

		if not self.has_metadata:
			return None

		begin, end = map(int, self._header.offsets.metadata)

		self.fobj.seek(begin)

		# Read format character
		fmt = self.fobj.read(1)

		if fmt == b'j':
			# JSON format

			data = self.fobj.read(end - begin)
			return json.loads(data.decode())

		else:
			raise ValueError(f'Unknown metadata format character {fmt.decode()!r}')

	def _read_ids(self):
		"""Read the IDs from the file if it has any.

		Returns
		-------
			IDs stored in file.
		"""

		begin, end = map(int, self._header.offsets.ids)

		self.fobj.seek(begin)

		fmt = self.fobj.read(1)

		if fmt == b'i':
			# Raw integer array

			# Numpy format is next two bytes
			dtype = np.dtype(self.fobj.read(2).decode())

			ids = read_npy(self.fobj, dtype, self.count)

		elif fmt == b's':
			# Null-terminated strings

			data = self.fobj.read(end - begin)
			if data[-1] != 0:
				raise ValueError('Error reading signature file IDs')

			ids = np.zeros(self.count, dtype='O')

			# Read strings
			start = 0
			for i in range(self.count):
				finish = data.find(b'\0', start)

				if finish < 0:
					raise ValueError('Error reading signature file IDs')

				ids[i] = data[start:finish].decode()

				start = finish + 1
				i += 1

			if finish != len(data) - 1:
				raise ValueError('Error reading signature file IDs')

		else:
			raise ValueError(f'Unknown ID format character {fmt.decode()!r}')

		# Make read-only to prevent accidental modification
		ids.flags.writeable = False
		return ids
