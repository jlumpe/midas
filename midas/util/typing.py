"""

"""

import sys
import typing
from typing import Union, Any


class TypeCheckError(Exception):
	"""Raised when attempting to check a value against an unsupported type annotation."""


def _is_union(T):
	"""Check if a type annotation is a *parameterized* :class:`typing.Union`.

	Parameters
	----------
	T
		Type annotation

	Returns
	-------
	bool
	"""
	if sys.version_info[1] >= 7:
		return isinstance(T, typing._GenericAlias) and T.__origin__ is typing.Union
	else:
		return isinstance(T, type(Union)) and T is not Union


def _union_types(T):
	"""Get the types from a parameterized :class:`typing.Union`.

	Parameters
	----------
	T
		Result of ``Union[A, B, ...]``.

	Returns
	-------
	tuple
	"""
	return T.__args__


def _is_generic_type(T):
	"""Check if a type annotation value corresponds to a generic type.

	Parameters
	----------
	T
		Type annotation

	Returns
	-------
	bool
	"""
	if sys.version_info[1] >= 7:
		return isinstance(T, typing._GenericAlias)
	else:
		return isinstance(T, typing.GenericMeta)


def _is_generic_parameterized(T):
	"""Check if a generic type annotation has (any) parameters specified.

	Parameters
	----------
	T
		Generic type annotation

	Returns
	-------
	bool
	"""
	if sys.version_info[1] >= 7:
		return any(not isinstance(arg, typing.TypeVar) for arg in T.__args__)
	else:
		return T.__args__ is not None


def _generic_base(T):
	"""Get the base (non-generic) type of a generic type annotation.

	This should return e.g. :class:`tuple` for ``typing.Tuple`` and
	:class:`collections.abe.Sequence` for ``typing.Sequence``.

	Parameters
	----------
	T
		Generic type annotation
	"""
	if sys.version_info[1] >= 7:
		return T.__origin__
	else:
		return T.__extra__


def _type_check_generic(value, T, ignore_params=False):
	"""Check value against generic type annotation T.

	Parameters
	----------
	value
		Value to check
	T
		Generic type annotation
	ignore_params : bool
		Ignore parameters of generic types instead of raising a :exc:`TypeCheckError`.

	Returns
	-------
	bool
	"""
	if not ignore_params and _is_generic_parameterized(T):
		raise TypeCheckError('Type checking against parameterized generic types not supported')

	base = _generic_base(T)
	if not isinstance(base, type):
		raise TypeCheckError(f'Type checking not supported for generic type {T} with base type {base}')

	return isinstance(value, base)


def type_check(value, T, ignore_params=False):
	"""Check a value against a type annotation at run time.

	Supports the following values for ``T``:

	- Any actual ``type``, in which case the builtin :func:`isinstance` is used.
	- :class:`typing.Any`
	- :class:`typing.Union`
	- Any *unparameterized* generic type
	- ``None``, which is used in type annotations to indicate :class:`NoneType`

	Parameters
	----------
	value
		Value to check
	T
		Type annotation
	ignore_params : bool
		Ignore parameters of generic types instead of raising a :exc:`TypeCheckError`.

	Returns
	-------
	bool

	Raises
	------
	TypeCheckError
		If ``T`` is not a supported type to check against.
	"""
	# None
	if T is None:
		return value is None

	# typing.Any
	if T is Any:
		return True

	# typing.Union
	if _is_union(T):
		return any(type_check(value, T2) for T2 in _union_types(T))

	# Generic
	if _is_generic_type(T):
		return _type_check_generic(value, T, ignore_params=ignore_params)

	# Actual type, use builtin isinstance()
	if isinstance(T, type):
		return isinstance(value, T)

	# Unsupported
	raise TypeCheckError(f'Type checking not supported for type {T}')
