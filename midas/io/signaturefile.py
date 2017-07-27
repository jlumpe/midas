"""Store sets of K-mer signatures on disk in binary format."""

import io
import json

import numpy as np

from midas.kmers import KmerSpec, SignatureArray
from .util import read_npy, write_npy, NamedStruct


class SignatureFile:
	"""Binary file for storing k-mer signatures.

	"""

	_MAGIC_NUMBER = b'MSF\xFF'

	_VERSION = b'1.00'

	_OFFSETS_DTYPE = np.dtype([
		('lengths', '2u8'),
		('metadata', '2u8'),
		('ids', '2u8'),
		('data', '2u8'),
	])

	_HEADER_DTYPE = np.dtype([
		('magic_number', ('S', len(_MAGIC_NUMBER))),
		('version', 'S4'),
		('count', 'u8'),
		('dtype', 'S2'),
		('offsets', _OFFSETS_DTYPE),
	])

	_LENGTHS_DTYPE = 'u4'

	DEFAULT_EXT = '.midas-signatures'

	def __init__(self, fobj):

		self.fobj = fobj
		self._header = NamedStruct._fromfile(self._HEADER_DTYPE, fobj)

		if self._header.magic_number != self._MAGIC_NUMBER:
			raise ValueError('File does not appear to be in the correct format')

		if self._header.version != self._VERSION:
			raise ValueError('Unexpected version identifier')

		self.version = self._header.version.decode('ascii')
		self.count = self._header.count
		self.dtype = np.dtype(self._header.dtype)

		if self._header.offsets.ids[0] > 0:
			self.ids = self._read_ids()
		else:
			self.ids = None

		self.has_metadata = self._header.offsets.metadata[0] > 0

		fobj.seek(self._header.offsets.lengths[0])
		self.lengths = read_npy(fobj, self._LENGTHS_DTYPE, shape=self.count)
		self.lengths.flags.writeable = False

		self.nelems = self.lengths.sum()

	def get_array(self):

		self.fobj.seek(self._header.offsets['data'][0])
		values = read_npy(self.fobj, self.dtype, shape=self.nelems)

		bounds = np.zeros(self.count + 1, dtype='u8')
		bounds[1:] = np.cumsum(self.lengths)

		return SignatureArray(values, bounds)

	def iter_signatures(self):

		self.fobj.seek(self._header.offsets['data'][0])

		for length in self.lengths:
			yield read_npy(self.fobj, self.dtype, shape=length)

	@classmethod
	def write_sigarray(cls, fobj, sigarray, *, dtype=None, ids=None, metadata=None):

		count = len(sigarray)

		if dtype is None:
			dtype = sigarray.values.dtype
		else:
			dtype = np.dtype(dtype)

		# Create header data
		header = NamedStruct(cls._HEADER_DTYPE)

		header.magic_number = cls._MAGIC_NUMBER
		header.version = cls._VERSION
		header.count = count
		header.dtype = dtype.str[1:]

		# Write header data (will come back and fill in offsets later)
		fobj.write(header._tobytes())

		# Write lengths
		lengths = np.diff(sigarray.bounds).astype(cls._LENGTHS_DTYPE)

		header.offsets.lengths[0] = fobj.tell()
		write_npy(fobj, lengths)
		header.offsets.lengths[1] = fobj.tell() - 1

		# Write metadata
		if metadata is not None:
			header.offsets.metadata[0] = fobj.tell()
			cls._write_metadata(fobj, metadata)
			header.offsets.metadata[1] = fobj.tell() - 1

		# Write IDs
		if ids is not None:
			if len(ids) != count:
				raise ValueError('Number of IDs does not match number of signatures')

			header.offsets.ids[0] = fobj.tell()
			cls._write_ids(fobj, ids)
			header.offsets.ids[1] = fobj.tell() - 1

		# Write data
		header.offsets.data[0] = fobj.tell()
		write_npy(fobj, sigarray.values)
		header.offsets.data[1] = fobj.tell() - 1

		# Go back and write offsets
		_, offsets_offset = cls._HEADER_DTYPE.fields['offsets']
		fobj.seek(offsets_offset)
		fobj.write(header.offsets._tobytes())

	@classmethod
	def _write_metadata(cls, fobj, metadata):

		# Support JSON format only for now

		fobj.write(b'j')

		wrapper = io.TextIOWrapper(fobj)
		json.dump(metadata, wrapper)
		wrapper.detach()

	def get_metadata(self):
		if not self.has_metadata:
			return None

		begin, end = map(int, self._header.offsets.metadata)

		self.fobj.seek(begin)

		# Read format character
		fmt = self.fobj.read(1)

		if fmt == b'j':
			# JSON format

			data = self.fobj.read(end - begin)
			return json.loads(data)

		else:
			raise ValueError(
				'Unknown metadata format character {!r}'
				.format(fmt.decode())
			)

	@classmethod
	def _write_ids(cls, fobj, ids):

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


	def _read_ids(self):

		begin, end = map(int, self._header.offsets.ids)

		self.fobj.seek(begin)

		fmt = self.fobj.read(1)

		if fmt == b'i':
			# Raw integer array

			# Numpy format is next two bytes
			dtype = np.dtype(self.fobj.read(2).decode())

			return read_npy(self.fobj, dtype, self.count)

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

			return ids

		else:
			raise ValueError(
				'Unknown ID format character {!r}'
				.format(fmt.decode())
			)
