"""Miscellaneous utility code."""

import os
import re
import collections


def kwargs_done(kwargs):
	"""Raises an error when an unknown keyword argument is encountered.

	Meant to be used with the pattern where extra keyword arguments to a
	function are caught with ``**kwargs`` and popped off as they are
	interpreted. At the end this function can be used to assert that the
	dictionary is empty.

	:param dict kwargs: Dictionary of keyword arguments, of which all should
		have been removed.
	:raises KeyError: If the dictionary is not empty, in a similar format
		to the built-in error for an unknown keyword argument.
	"""
	if kwargs:
		raise KeyError(
			'Unknown keyword argument {:r}'
			.format(next(iter(kwargs)))
		)


def sanitize_filename(name, replace='_'):
	"""Replace invalid characters in a file name.

	:param str name: File name to sanitize.
	:param str replace: String to replace invalid charcters with.
	:returns: Name with non-alphanumeric characters other than ``_.-`` replaced.
	:rtype: str
	"""
	return re.sub(r'[^A-Za-z0-9_.-]+', replace, name)


class SubPath:
	"""Data descriptor for getting absolute path from relative subpath.

	Intended for classes that manage a single root directory (e.g.,
	``BasicDatabase``) that have subpaths they need to get the absolute path
	for. To avoid many os.path.join()'s all over the place, the descriptor
	returns the absolute path on instances and a function giving the absolute
	path given the root path on classes.

	:param str path: Path relative to parent object's root directory.
	"""

	def __init__(self, path):
		self.path = path

	def get_abs_path(self, root_dir):
		return os.path.join(root_dir, self.path)

	@classmethod
	def _get_root_dir(cls, obj, type_):
		"""Get the root directory of an object."""
		root_dir_attr = getattr(type_, '__root_dir_attr__', 'root_dir')
		return getattr(obj, root_dir_attr)

	def __get__(self, obj, type_):
		if obj is None:
			return self.get_abs_path

		else:
			return self.get_abs_path(self._get_root_dir(obj, type_))

	def __set__(self, obj, value):
		raise AttributeError("Can't set attribute")
