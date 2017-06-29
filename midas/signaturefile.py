"""Store sets of K-mer signatures on disk in HDF5 format."""

import os

import numpy as np
import h5py

from midas.kmers import KmerSpec


class SignatureFile:

	FORMAT_VERSION = '1.0'

	def __init__(self, file, mode='a', **kwargs):

		if isinstance(file, h5py.File):
			if mode is not None or kwargs:
				raise TypeError(
					'Cannot specify mode or additional keyword arguments when '
					'using existing file object'
				)

			self.h5file = file

		else:
			if mode == 'r':
				file_exists = True
				self.writable = False

			elif mode == 'r+':
				file_exists = True
				self.writable = True

			elif mode in ['w', 'w-', 'x']:
				file_exists = False
				self.writable = True

			elif mode == 'a':
				self.writable = True

				if os.path.isfile(file):
					file_exists = True
					mode = 'r+'
				else:
					file_exists = False
					mode = 'x'

			else:
				raise ValueError('Invalid mode')

			if file_exists:
				self.h5file = self._open_existing(file, mode, kwargs)
			else:
				self.h5file = self._create_file(file, mode, kwargs)

	@classmethod
	def _open_existing(cls, filename, mode, kwargs):

		file = h5py.File(filename, mode, **kwargs)

		try:
			file_ver = file.attrs['MIDAS_SIGNATUREFILE_VERSION']

		except KeyError:
			raise OSError('File is not a MIDAS signature file')

		if file_ver != cls.FORMAT_VERSION:
			raise OSError(
				'File is version {}, current supported version is {}'
				.format(file_ver, cls.FORMAT_VERSION)
			)

		return file

	@classmethod
	def _create_file(cls, filename, mode, kwargs):

		file = h5py.File(filename, mode, **kwargs)

		file.attrs['MIDAS_SIGNATUREFILE_VERSION'] = cls.FORMAT_VERSION
		file.create_group('SETS')

		return file

	def _sets_group(self):
		return self.h5file['SETS']

	def list_sets(self):
		for group in self._sets_group():
			yield group.attrs['k'], group.attrs['prefix']

	@staticmethod
	def _get_set_group_name(k, prefix):

		if isinstance(prefix, bytes):
			prefix = prefix.decode('ascii')

		return '{}-{}'.format(int(k), prefix)

	@staticmethod
	def _format_signature_params(k, prefix):

		k = int(k)
		if k <= 0:
			raise ValueError('k must be a positive integer')

		if isinstance(prefix, bytes):
			prefix = prefix.decode('ascii')
		prefix = prefix.upper()
		if not all(c in 'ATGC' for c in prefix):
			raise ValueError('Prefix must contain characters ACGT only')

		return k, prefix

	def get_set(self, k, prefix):

		k, prefix = self._format_signature_params(k, prefix)

		sets_group = self._sets_group()
		group = sets_group[self._get_set_group_name(k, prefix)]

		if group is not None:
			return SignatureSet(self, group, k, prefix)
		else:
			raise KeyError((k, prefix))

	def create_set(self, k, prefix, exist_ok=False):

		k, prefix = self._format_signature_params(k, prefix)

		sets_group = self._sets_group()

		group_name = self._get_set_group_name(k, prefix)

		group = sets_group[group_name]

		if group is not None:
			if not exist_ok:
				raise ValueError('Group exists')

		else:
			group = sets_group.create_group(group_name)
			group.attrs['k'] = k
			group.attrs['prefix'] = prefix

		return SignatureSet(self, group, k, prefix)

	def remove_set(self, k, prefix):
		raise NotImplementedError()

	def __repr__(self):
		return '{}({!r}, {!r})'.format(
			type(self).__name__,
			self.h5file.filename,
			self.h5file.mode
		)


class StoredSignatureSet:

	def __init__(self, file, group, k, prefix):
		self._file = file
		self._group = group
		self._k = k
		self._prefix = prefix
		self._kspec = KmerSpec(k, prefix)

	@property
	def k(self):
		return self._k

	@property
	def prefix(self):
		return self._prefix

	def kspec(self):
		return self._kspec

	@staticmethod
	def _check_key(key):
		key = key.strip('/')

		if not all(key.split('/')):
			raise ValueError('Key may not contain consecutive forward slashes')

		if '#' in key:
			raise ValueError('Key may not contain # character')

		return key

	@staticmethod
	def _check_version(version):
		if not version:
			raise ValueError('Version string must be non-empty')

		if not all(c in '0123456789.' for c in version):
			raise ValueError('Version string must contain digits and dots only')

	def list_signatures(self):
		pass

	def get_signature(self, key, version=None):
		pass

	def store_signature(self, signature, key, version=None):
		pass

	def remove_signature(self, key, version=None):
		pass

	def __repr__(self):
		return '<{} {}/{} {!r})>'.format(
			type(self).__name__,
			self.k,
			self.prefix.decode('ascii'),
			self._file,
		)
