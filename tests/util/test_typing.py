"""Test the midas.util.typing submodule."""

import typing
from typing import Union, Any, Optional
from collections import abc

import pytest

import midas.util.typing as mtyping
from midas.util.typing import type_check, TypeCheckError


class TestHelpers:
	"""Test private helper funcs."""

	def test_is_union(self):
		"""Test _is_union function."""
		assert mtyping._is_union(Union[int, str])
		assert mtyping._is_union(Optional[int])
		assert not mtyping._is_union(Union)
		assert not mtyping._is_union(Optional)
		assert not mtyping._is_union(int)
		assert not mtyping._is_union(None)
		assert not mtyping._is_union(typing.List)
		assert not mtyping._is_union(typing.Any)

	def test_union_types(self):
		"""Test _union_types function."""
		assert mtyping._union_types(Union[int, str]) == (int, str)
		assert mtyping._union_types(Optional[int]) == (int, type(None))

	def test_is_generic_type(self):
		"""Test _is_generic_type function."""
		assert mtyping._is_generic_type(typing.List)
		assert mtyping._is_generic_type(typing.List[int])
		assert mtyping._is_generic_type(typing.Tuple)
		assert mtyping._is_generic_type(typing.Tuple[int, str])
		assert mtyping._is_generic_type(typing.Callable)
		assert mtyping._is_generic_type(typing.Callable[[int], str])
		assert not mtyping._is_generic_type(int)
		assert not mtyping._is_generic_type(abc.Callable)
		# assert not mtyping._is_generic_type(Union[str, int])

	def test_is_generic_parameterized(self):
		"""Test _is_generic_parameterized function."""
		assert not mtyping._is_generic_parameterized(typing.List)
		assert mtyping._is_generic_parameterized(typing.List[int])
		assert not mtyping._is_generic_parameterized(typing.Tuple)
		assert mtyping._is_generic_parameterized(typing.Tuple[int, str])
		assert not mtyping._is_generic_parameterized(typing.Callable)
		assert mtyping._is_generic_parameterized(typing.Callable[[int], str])

	def test_generic_base(self):
		"""Test _generic_base function."""
		assert mtyping._generic_base(typing.List) is list
		assert mtyping._generic_base(typing.List[int]) is list
		assert mtyping._generic_base(typing.Tuple) is tuple
		assert mtyping._generic_base(typing.Tuple[int, str]) is tuple
		assert mtyping._generic_base(typing.Callable) is abc.Callable
		assert mtyping._generic_base(typing.Callable[[int], str]) is abc.Callable


class TestTypeCheck:
	"""Test type_check function."""

	def test_type(self):
		# 2nd argument is an actual type
		assert type_check(4, int)
		assert type_check(4, object)
		assert type_check('foo', str)
		assert type_check([], abc.Sequence)
		assert not type_check(4, str)
		assert not type_check([], abc.Mapping)

	def test_any(self):
		assert type_check(1, Any)
		assert type_check('foo', Any)
		assert type_check(None, Any)

	def test_union(self):
		assert type_check(1, Union[int, str])
		assert type_check('foo', Union[int, str])
		assert not type_check([], Union[int, str])
		assert type_check(1, Optional[int])
		assert type_check(None, Optional[int])
		assert not type_check('foo', Optional[int])

	def test_none(self):
		assert type_check(None, None)
		assert not type_check(int, None)

	def test_generic(self):
		assert type_check([], typing.Sequence)
		assert not type_check([], typing.Mapping)
		assert not type_check(dict(), typing.Sequence)
		assert type_check(dict(), typing.Mapping)
		assert type_check((1, 'foo'), typing.Tuple)
		assert type_check(lambda x: x**2, typing.Callable)

	def test_generic_parameterized(self):
		assert type_check([1, 2, 3], typing.Sequence[int], ignore_params=True)
		assert type_check([1, 2, 3], typing.Sequence[str], ignore_params=True)
		assert not type_check({'a': 1, 'b': 2}, typing.Sequence[int], ignore_params=True)
		assert type_check((1, 'foo'), typing.Tuple[int, str], ignore_params=True)
		assert type_check((1, 'foo'), typing.Tuple[str, int], ignore_params=True)
		assert not type_check([1, 'foo'], typing.Tuple[int, str], ignore_params=True)
		assert type_check(lambda x: x**2, typing.Callable[[int], int], ignore_params=True)

		with pytest.raises(TypeCheckError):
			type_check([1, 2, 3], typing.Sequence[int])
		with pytest.raises(TypeCheckError):
			type_check([1, 2, 3], typing.Sequence[str])
		with pytest.raises(TypeCheckError):
			type_check({'a': 1, 'b': 2}, typing.Sequence[int])
		with pytest.raises(TypeCheckError):
			type_check((1, 'foo'), typing.Tuple[int, str])
		with pytest.raises(TypeCheckError):
			type_check((1, 'foo'), typing.Tuple[str, int])
		with pytest.raises(TypeCheckError):
			type_check([1, 'foo'], typing.Tuple[int, str])
		with pytest.raises(TypeCheckError):
			type_check(lambda x: x**2, typing.Callable[[int], int])

	def test_invalid(self):
		with pytest.raises(TypeCheckError):
			type_check(0, Union)
