"""Convert data to/from JSON format.

The :func:`.dump` and  :func:`.load` functions in this module work just like their built-in
equivalents in the :mod:`json` module, but support additional types such as ``attrs``-defined
classes.
"""

import json
from typing import Any

import cattr
import numpy as np


converter = cattr.Converter()


def to_json(obj):
	"""Convert object to JSON-writable data (anything that can be passed to :func:`json.dump`).

	Parameters
	----------
	obj
		Object to convert.
	"""
	return converter.unstructure(obj)


def from_json(data, cls=Any):
	"""Load object from parsed JSON data.

	Parameters
	----------
	data
		Data parsed from JSON format.
	cls
		Type to load.

	Returns
	-------
	Instance of ``cls``
	"""
	return converter.structure(data, cls)


def dump(obj, f):
	"""Write the JSON representation of an object to a file.

	Parameters
	----------
	obj
		Object to write.
	f
		Writeable file object in text mode.
	"""
	data = to_json(obj)
	json.dump(data, f)


def load(f, cls=Any):
	"""Load an object from a JSON file.

	Parameters
	----------
	f
		Readable file object in text mode.
	cls
		Type to load.

	Returns
	-------
	Instance of ``cls``
	"""
	data = json.load(f)
	return from_json(data, cls)


def dumps(obj):
	"""Get the JSON representation of an object as a string.

	Parameters
	----------
	obj
		Object to write.

	Returns
	-------
	str
	"""
	return json.dumps(to_json(obj))


def loads(s, cls=Any):
	"""Load an object from a JSON string.

	Parameters
	----------
	s : str
		String containing JSON-encoded data.
	cls
		Type to load.

	Returns
	-------
	Instance of ``cls``
	"""
	data = json.loads(s)
	return from_json(data, cls)


# Numpy scalars
converter.register_unstructure_hook(np.integer, int)
converter.register_unstructure_hook(np.floating, float)


class Jsonable:
	"""Mixin class that provides custom JSON conversion methods.

	Either of the special methods ``__to_json__`` and ``__from_json__`` may be set to ``None`` to
	indicate that the default conversion process should be used.

	.. py:method:: __to_json__()
		:abstractmethod:

		Convert the instance to JSON-writable data (anything that can be passed to :func:`json.dump`).


	.. py:classmethod:: __from_json__(data)
		:abstractmethod:

		Create an instance from parsed JSON data.
	"""
	__to_json__ = None
	__from_json__ = None

converter.register_structure_hook_func(
	lambda cls: issubclass(cls, Jsonable) and cls.__from_json__ is not None,
	lambda data, cls: cls.__from_json__(data),
)

converter.register_unstructure_hook_func(
	lambda cls: issubclass(cls, Jsonable) and cls.__to_json__ is not None,
	lambda obj: obj.__to_json__(),
)
