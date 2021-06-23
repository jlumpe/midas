"""Test the midas.util.attr_ext submodule."""

import pytest
from attr import Factory, NOTHING
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


@pytest.mark.parametrize('optional', [False, True])
@pytest.mark.parametrize('validator', [False, True])
@pytest.mark.parametrize('converter', [False, True])
@pytest.mark.parametrize('default', [False, True])
@pytest.mark.parametrize('validate_type', [False, True])
@pytest.mark.parametrize('use_decorators', [False, True])
def test_attrib_opts(optional, validator, converter, validate_type, default, use_decorators):
	"""Test custom arguments to attrib() and other arguments they may interact with."""

	default_value = 10

	def _validate_x(self, attribute, value):
		if value < 3:
			raise ValueError('x must be at least 3')

	if use_decorators:
		@attrs()
		class TestCls:
			x: int = attrib(
				optional=optional,
				validator=_validate_x if validator else None,
				converter=int if converter else None,
				default=default_value if default else NOTHING,
				validate_type=validate_type,
			)

	else:
		@attrs()
		class TestCls:
			x: int = attrib(
				optional=True,
				converter=int if converter else None,
				validate_type=validate_type,
			)

			if validator:
				x.validator(_validate_x)

			if default:
				@x.default
				def _default_x(self):
					return default_value

	# Explicit value that passes all checks
	assert TestCls(3).x == 3

	# Check default value, if any
	if default:
		assert TestCls().x == default_value
	elif optional:
		assert TestCls().x is None

	# Check explicit None value if optional, should bypass validators and converters
	if optional:
		assert TestCls(None).x is None

	# Check conversion
	if converter:
		for arg in ['3', 3.0]:
			x = TestCls(arg).x
			assert isinstance(x, int) and x == 3

	# Check validator
	if validator:
		with pytest.raises(ValueError):
			TestCls(1)

		# Validate converted value
		if converter:
			with pytest.raises(ValueError):
				TestCls('1')

	else:
		assert TestCls(1).x == 1

	# Check type validation (only makes sense if not using converter)
	if not converter:
		# Value not of correct type but passing validator
		if validate_type:
			with pytest.raises(TypeError):
				TestCls(3.0)

		else:
			x = TestCls(3.0).x
			assert isinstance(x, float) and x == 3
