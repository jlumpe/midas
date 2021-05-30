"""Test the midas.util.attr_ext submodule."""

import pytest
from midas.util.attr import attrs, attrib


@pytest.mark.parametrize('order', [False, True])
def test_attrs_order_default(order):
	"""Test that the `order` argument of `attrs` defaults to False."""

	class TestCls:
		x: int = attrib()
		y: int = attrib()

	if order:
		TestCls = attrs(order=True)(TestCls)
	else:
		TestCls = attrs()(TestCls)

	# Check equality always works
	assert TestCls(1, 2) == TestCls(1, 2)
	assert TestCls(1, 2) != TestCls(3, 4)
	assert TestCls(1, 2) != 'foo'

	# Ordered comparisons should work
	if order:
		assert TestCls(1, 2) < TestCls(3, 4)

	# Methods for ordered comparisons have been removed
	else:
		with pytest.raises(TypeError):
			TestCls(1, 2) < TestCls(3, 4)


def test_attrib_optional():
	"""Test `optional` argument to `attrib`."""

	@attrs()
	class TestCls:
		x: str = attrib()
		y: int = attrib(optional=True)

	assert TestCls('foo').y is None
	assert TestCls('foo', 1).y == 1

	# Test with validator and converter
	def _validate_y(self, attribute, value):
		if value < 3:
			raise ValueError('y must be at least 3')

	@attrs()
	class TestCls:
		x: str = attrib()
		y: int = attrib(optional=True, converter=int, validator=_validate_y)

	assert TestCls('foo').y is None
	assert TestCls('foo', 3).y == 3
	assert TestCls('foo', '3').y == 3

	with pytest.raises(ValueError):
		TestCls('foo', 1)


def test_attrib_validate_type():
	"""Test `validate_type` argument to `attrib`."""

	@attrs()
	class TestCls:
		x: int = attrib(validate_type=False)
		y: int = attrib(validate_type=True)

	TestCls(1, 2)
	TestCls('1', 2)

	with pytest.raises(TypeError):
		TestCls(1, '2')

	# Test with attribute that has no type
	@attrs()
	class TestCls2:
		x = attrib(validate_type=True)

	TestCls2(1)
	TestCls2('foo')
	TestCls2(None)
