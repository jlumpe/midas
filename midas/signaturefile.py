"""Store sets of K-mer signatures on disk in HDF5 format."""

import os

import numpy as np
import h5py

from midas.kmers import KmerSpec


class SignatureFile:

	FORMAT_VERSION = '1.0'

	def __init__(self, file, mode=None, **kwargs):

		if isinstance(mode, h5py.File):
			if mode is not None or kwargs:
				raise TypeError(
					'Cannot specify mode or additional keyword arguemnts when '
					'using existing file object'
				)

			self.h5file = file

		else:
			if mode is None:
				mode = 'a'

			if mode in ['r', 'r+']:
				file_exists = True

			elif mode in ['w', 'w-', 'x']:
				file_exists = False

			elif mode == 'a':
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

	def list_sets(self):
		pass

	def set_exists(self, key, version):
		pass

	def get_set(self, key, version=None):

		key = self._check_key(key)
		self._check_version(version)

		sets_group = self.h5file['SETS']

		if version is not None:
			group = sets_group.get(key + '/#' + version)

			if group is not None:
				return SignatureSet(group, key, version)
			else:
				return None

		else:
			raise NotImplementedError()

	def create_set(self, key, version, k, prefix):

		key = self._check_key(key)
		self._check_version(version)

		k = int(k)
		if k <= 0:
			raise ValueError('k must be a positive integer')

		if isinstance(prefix, bytes):
			prefix = prefix.decode('ascii')
		prefix = prefix.upper()
		if not all(c in 'ATGC' for c in prefix):
			raise ValueError('Prefix must contain characters ACGT only')

		sets_group = self.h5file['SETS']

		key_group = sets_group.create_group(key)
		group = key_group.create_group('#' + version)

		group.attrs['k'] = k
		group.attrs['prefix'] = prefix

		return SignatureSet(group, key, version)

	def remove_set(self, key, version):
		pass

	def __repr__(self):
		return '{}({!r}, {!r})'.format(
			type(self).__name__,
			self.h5file.filename,
			self.h5file.mode
		)


class SignatureSet:

	def __init__(self, group, key, version):
		self._key = key
		self._version = version
		self._group = group

	@property
	def key(self):
		return self._key

	@property
	def version(self):
		return self._version

	@property
	def k(self):
		return self.group.attrs['k']

	@property
	def prefix(self):
		return self.group.attrs['prefix'].encode('ascii')

	def kspec(self):
		return KmerSpec(self.k, self.prefix)

	def list_signatures(self):
		pass

	def get_signature(self, key):
		pass

	def store_signature(self, key, signature):
		pass

	def remove_signature(self, key):
		pass

	def __repr__(self):
		return '<{} {!r} (version={!r})>'.format(
			type(self).__name__,
			self.key,
			self.version
		)
