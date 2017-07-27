"""Utility code for reading/writing data files."""

import numpy as np


def read_npy(fobj, dtype, shape=1):
	array = np.zeros(shape, dtype=dtype)
	fobj.readinto(array.data)
	return array


def write_npy(fobj, array):
	if array.ndim != 1:
		raise ValueError('Array must be one-dimensional')

	if array.flags.c_contiguous:
		data = array.data
	else:
		data = array.tobytes()

	written = fobj.write(data)

	if written != array.nbytes:
		raise IOError('Could not write all bytes of array')


class NamedStruct:

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
		return self._dtype.size

	def __getitem__(self, which):
		name = self._dtype.names[which] if isinstance(which, int) else which

		field_dtype, field_offset = self._dtype.fields[name]

		if field_dtype.fields is not None:
			return NamedStruct(field_dtype, self._data[name])
		else:
			return self._data[name][()]

	def __setitem__(self, name, value):
		self._data[name] = value

	def __getattr__(self, name):
		if not name.startswith('_'):
			try:
				return self[name]
			except KeyError:
				pass

		return object.__getattribute__(self, name)

	def __setattr__(self, name, value):
		if name.startswith('_'):
			object.__setattr__(self, name, value)

		else:
			self._data[name] = value

	def _tobytes(self):
		return self._data.tobytes()

	@classmethod
	def _fromfile(cls, dtype, fobj):
		struct = NamedStruct(dtype)
		fobj.readinto(struct._data.data)
		return struct

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
