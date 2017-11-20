"""Test SQLAlchemy utility code in midas.db.sqla."""

import pytest
import sqlalchemy as sa

from midas.db import sqla


@pytest.mark.parametrize('nullable', [False, True])
def test_keyvaluetable(tmpdir, nullable):
	"""Test KeyValueTable."""

	# Create SQLAlchemy table
	metadata = sa.MetaData()
	table = sqla.KeyValueTable.make_table('keyvalue', metadata, nullable=nullable)

	# Create test database
	engine = sa.create_engine('sqlite:///' + str(tmpdir / 'kv.db'))
	metadata.create_all(engine)

	# Create KeyValueTable instance
	kv = sqla.KeyValueTable(engine, table)

	# Check initially empty
	assert len(kv) == 0
	assert list(iter(kv)) == []
	assert 'foo' not in kv

	# Check using nonexistant key
	with pytest.raises(KeyError):
		kv['foo']

	with pytest.raises(KeyError):
		del kv['foo']

	# Add a key/value pair
	kv['foo'] = '123'
	assert len(kv) == 1
	assert 'foo' in kv
	assert kv['foo'] == '123'

	# Try deleting key
	del kv['foo']
	assert 'foo' not in kv
	assert len(kv) == 0

	# Test clearing
	kv.update(foo='123', bar='456', baz='789')
	assert len(kv) == 3

	kv.clear()

	assert len(kv) == 0
	for key in ['foo', 'bar', 'baz']:
		assert key not in kv

	# Test non-string key and value
	with pytest.raises(TypeError):
		kv[0] = 'valid_value'
	with pytest.raises(TypeError):
		kv['valid_key'] = 0

	# Check null value
	if nullable:
		kv['foo'] = None
		assert kv['foo'] is None

	else:
		with pytest.raises(TypeError):
			kv['foo'] = None
