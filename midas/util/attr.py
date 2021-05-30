"""
Extensions to functionality of :mod:`attr` package.
"""

import attr

from .typing import type_check


def validate_field_type(instance, attribute, value):
	"""Validator function which ensures a value matches an attribute's type.

	Uses :func:`midas.util.typing.type_check` which supports a subset of type objects from the
	builtin :mod:`typing` module. Ignores parameters on generic types like ``Sequence[int]``, which
	are otherwise unsupported by this function (see the ``ignore_params`` argument).

	Parameters
	----------
	instance
		Object being validated.
	attribute : attr.Attribute
		Attribute being validated.
	value
		Value to check.

	Raises
	------
	TypeError
		If ``value`` is not an instance of ``attribute.type``.
	"""
	if attribute.type is not None and not type_check(value, attribute.type, ignore_params=True):
		raise TypeError(f'{attribute.name} must be of type {attribute.type}, got {repr(value)}')


def attrs(**kw):
	"""
	Class decorator which works like :func:`attr.s` with some minor changes.

	Parameters
	----------
	\\**kw
		Keyword arguments to :func:`attr.s`.


	Returns
	-------
	Callable
		Class decorator function
	"""
	# Have order default to False
	if kw.get('order', None) is not True:
		kw['order'] = False
	return attr.s(**kw)


def attrib(optional=False, validate_type=False, **kw):
	"""
	Extended version of :func:`attr.attrib`.

	Parameters
	----------
	optional : bool
		If True the field will also accept None as well as any other valid values. In this case if
		the field is initialized with a value of None all other validators and converters will be
		bypassed. The ``default`` argument will also be set to None, if it has not been set
		already.
	validate_type : bool
		If True will add an additional validator that checks that the value matches the attribute's
		type (see :func:`.validate_field_type`).
	\\**kw
		Keyword arguments to :func:`attr.attrib`.

	Returns
	-------
		The intermediate object returned by :func:`attr.attrib`. This will be turned into a proper
		:class:`attr.Attribute` after the class is initialized.
	"""
	if validate_type:
		if kw.get('validator', None) is not None:
			kw['validator'] = attr.validators.and_(kw['validator'], validate_field_type)
		else:
			kw['validator'] = validate_field_type

	if optional:
		kw.setdefault('default', None)

		if kw.get('validator', None) is not None:
			kw['validator'] = attr.validators.optional(kw['validator'])

		if kw.get('converter', None) is not None:
			kw['converter'] = attr.converters.optional(kw['converter'])

	return attr.ib(**kw)
