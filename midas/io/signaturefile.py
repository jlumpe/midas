"""Store sets of K-mer signatures on disk in binary format."""

import io
import json

import numpy as np

from midas.kmers import KmerSpec, SignatureArray
from .util import read_npy, write_npy, NamedStruct


class SignatureFile:
	"""Binary file for storing k-mer signatures.

	Constructor creates an object to read data from a stream. Use the
	:meth:`write` method to write data to a file.

	Acts as a context manager which closes its stream on exit.

	.. attribute:: fobj

		Stream file is being read from.

	.. attribute:: count

		Number of signatures in the file.

	.. attribute:: dtype

		:class:`numpy.dtype` of signatures in the file.

	.. attribute:: lengths

		:class:`numpy.ndarray` of lengths of each signature in the file.

	.. attribute:: nelems

		Total number of elements of all signatures in the file.

	.. attribute:: ids

		IDs for signatures as :class:`numpy.ndarray` of integers or string
		objects, or None if the file doesn't have IDs.

	.. attribute:: has_metadata

		True if file has metadata stored, False otherwise.

	.. attribute:: closed

		True if :attr:`fobj` is closed, False otherwise.

	:param fobj: Readable file-like object in binary mode.
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

		:param indices: List or array of indices of a subset of signatures to
			get. Array signatures will then correspond to the ids
			`signaturefile.ids[indices]``. If None will get all signatures.
		:param progress: Progress callback that will be called periodically
			with the number of signatures read. Will be called with positional
			arguments ``(ncompleted, total)`` where ``ncompleted`` is the total
			number of signatures read so far and ``total`` is the total number
			to read. Will be called after every signature has been read if
			``indices`` is not None or after every ``chunksize`` signatures
			otherwise.
		:type progress: function
		:param int chunksize: Number of signatures to read at a time, if
			``indices`` is None (``indices`` and ``chunksize`` should not both
			be given). Only useful with the ``progress`` argument. If None
			(default) will read all signatures in one shot.
		:rtype: midas.kmers.SignatureArray
		"""

		# File position at start of data
		data_start = self._header.offsets['data'][0]

		if indices is None:
			# Reading whole array in order

			if chunksize is None:
				chunksize = self.count

			# Array to read into
			array = SignatureArray.empty(self.lengths, dtype=self.dtype)

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
				raise TypeError('chunksize and indices cannot both be None')

			# Use this dtype as it is what the Cython metric functions expect
			from midas.cython.metrics import BOUNDS_DTYPE

			# Calculate file sub-array bounds from lengths
			bounds = np.zeros(self.count + 1, dtype=BOUNDS_DTYPE)
			bounds[1:] = np.cumsum(self.lengths)

			# Calculate bounds for indices only
			array = SignatureArray.empty(self.lengths[indices], dtype=self.dtype)

			# Read one at a time, in order of position in file
			index_pairs_sorted = sorted((j, i) for i, j in enumerate(indices))

			for i, (file_idx, out_idx) in enumerate(index_pairs_sorted):

				self.fobj.seek(data_start + bounds[file_idx] * self.dtype.itemsize)
				signature = read_npy(self.fobj, self.dtype, self.lengths[file_idx])

				array[out_idx] = signature

				# Progress callback
				if progress is not None:
					progress(i + 1, len(indices))

			return array

	def iter_signatures(self):
		"""Iterate over signatures in the file.

		:returns: Generator yielding :class:`numpy.ndarray`.
		"""

		self.fobj.seek(self._header.offsets['data'][0])

		for length in self.lengths:
			yield read_npy(self.fobj, self.dtype, shape=length)

	@classmethod
	def write(cls, fobj, signatures, *, dtype=None, ids=None, metadata=None):
		"""Write signatures to file.

		:param fobj: Writeable file-like object in binary mode.
		:param signatures: Sequence of signatures, as numpy arrays.
		:param dtype: Numpy data type of signatures to write (should be unsigned
			integer). If None will be determined automatically by signatures.
		:type: numpy.dtype
		:param ids: Sequence of IDs for the signatures, for example accession
			numbers of the genomes the signatures are from. Must be the same
			length as ``signatures`` and items must be all strings or integers.
		:param metadata: Arbitrary metadata that will be included with file.
			Must be JSONable.
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

		:returns: Tuple of begin and end file offsets for the metadata section.
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

		:returns: Tuple of begin and end file offsets for the IDs section.
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

		:returns: File's metadata object (typically dict), or None if the file
			has no metadata.
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
			raise ValueError(
				'Unknown metadata format character {!r}'
				.format(fmt.decode())
			)

	def _read_ids(self):
		"""Read the IDs from the file if it has any.

		:returns: IDs stored in file.
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
			raise ValueError(
				'Unknown ID format character {!r}'
				.format(fmt.decode())
			)

		# Make read-only to prevent accidental modification
		ids.flags.writeable = False
		return ids
