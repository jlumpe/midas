"""Miscellaneous utility code"""

import os
from functools import wraps



def _get_root_dir(obj, type):
	"""Helper function for SubPath and DecoratedSubPath"""
	root_dir_attr = getattr(type, '__root_dir_attr__', 'root_dir')
	return getattr(obj, root_dir_attr)


class SubPath:
	"""Data descriptor for getting absolute path from relative subpath

	Intended for classes that manage a single root directory (e.g.,
	BasicDatabase) that have subpaths they need to get the absolute path for.
	To avoid many os.path.join()'s all over the place, the descriptor returns
	the absolute path on instances and a function giving the absolute path
	given the root path on classes.
	"""

	def __init__(self, path):
		self.path = path

	def get_abs_path(self, root_dir):
		return os.path.join(root_dir, self.path)

	def __get__(self, obj, type):
		if obj is None:
			return self.get_abs_path

		else:
			return self.get_abs_path(_get_root_dir(obj, type))

	def __set__(self, obj):
		raise AttributeError("Can't set attribute")


class SubPathMethod:

	def __init__(self, func):
		self.func = func

	def get_subpath(self, root_dir, *args, **kwargs):
		return os.path.join(root_dir, self.func(*args, **kwargs))

	def __get__(self, obj, type):
		if obj is None:
			return self.get_subpath

		else:
			@wraps(self.func)
			def get_subpath(*args, **kwargs):
				return os.path.join(_get_root_dir(obj, type),
				                    self.func(*args, **kwargs))

			return get_subpath


def subpath(arg):
	"""Creates a SubPath/SubPathMethod (as decorator) from string/callable"""
	if isinstance(arg, str):
		return SubPath(arg)
	elif callable(arg):
		return SubPathMethod(arg)
	else:
		raise TypeError(type(arg))
