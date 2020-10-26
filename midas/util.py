"""Miscellaneous utility code."""

from pathlib import Path
import re


def kwargs_done(kwargs):
	"""Raises an error when an unknown keyword argument is encountered.

	Meant to be used with the pattern where extra keyword arguments to a
	function are caught with ``**kwargs`` and popped off as they are
	interpreted. At the end this function can be used to assert that the
	dictionary is empty.

	Parameters
	----------
	kwargs : dict
		Dictionary of keyword arguments, of which all should have been removed.

	Raises
	------
	KeyError
		If the dictionary is not empty, in a similar format to the built-in error for an unknown
		keyword argument.
	"""
	if kwargs:
		raise KeyError(
			'Unknown keyword argument {!r}'
			.format(next(iter(kwargs)))
		)


def sanitize_filename(name, replace='_'):
	"""Replace invalid characters in a file name.

	Parameters
	----------
	name : str
		File name to sanitize.
	replace : str
		String to replace invalid characters with.

	Returns
	-------
	str
		Name with non-alphanumeric characters other than ``_.-`` replaced.
	"""
	return re.sub(r'[^A-Za-z0-9_.-]+', replace, name)


class SubPath:
	"""Data descriptor for getting absolute path from relative subpath.

	Intended for classes that manage a single root directory (e.g.,
	``BasicDatabase``) that have subpaths they need to get the absolute path
	for. To avoid many ``os.path.join()``\ 's all over the place, the descriptor
	returns the absolute path on instances and a function giving the absolute
	path given the root path on classes.

	Using the division operator with a SubPath instance on the left is syntactic
	sugar for the :meth:`join` method.

	Parameters
	----------
	path
		Path relative to parent object's root directory.
	"""

	def __init__(self, path):
		self.path = Path(path)
		if self.path.is_absolute():
			raise ValueError('Path must be relative')

	def get_abs_path(self, root_dir):
		return root_dir / self.path

	@classmethod
	def _get_root_dir(cls, obj, type_):
		"""Get the root directory of an object."""
		root_dir_attr = getattr(type_, '__path_attr__', 'path')
		return getattr(obj, root_dir_attr)

	def __get__(self, obj, type_):
		if obj is None:
			return self.get_abs_path

		else:
			return self.get_abs_path(self._get_root_dir(obj, type_))

	def __set__(self, obj, value):
		raise AttributeError("Can't set attribute")

	def join(self, *args):
		"""Create a new SubPath by joining this path with a sub(-sub?)-path.

		Parameters
		----------
		\\*args
			See :meth:`pathlib.Path.joinpath`.

		Returns
		-------
		.SubPath
		"""
		return SubPath(self.path.joinpath(*args))

	def __truediv__(self, subpath):
		return self.join(subpath)


def path_str(path):
	"""Get a path string from a path-like object.

	On Python 3.6 this uses the path-like interface (``__fspath__`` method).
	On 3.5 only recognizes :class:`pathlib.Path` objects. Strings are passed
	through. If the object is :class:`bytes` or the ``__fspath`__`` method
	returns bytes they are decoded as UTF-8.

	Parameters
	----------
	path
		String, bytes, or path-like object representing a filesystem path.

	Returns
	-------
	str
		Path string.
	"""

	if hasattr(path, '__fspath__'):
		path = path.__fspath__()

	elif isinstance(path, Path):
		path = str(path)

	elif not isinstance(path, (str, bytes)):
		raise TypeError('{!r} is not path-like'.format(type(path).__name__))

	if isinstance(path, bytes):
		return path.decode()
	else:
		return path
