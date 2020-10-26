"""Custom types and other utilities for SQLAlchemy."""

import json
import collections
import weakref

import sqlalchemy as sa
from sqlalchemy.types import TypeDecorator
from sqlalchemy.ext.mutable import Mutable
from sqlalchemy.orm.session import Session


class ReadOnlySession(Session):
	"""Session class that doesn't allow flushing/committing"""

	def flush(self, *args, **kwargs):
		# Make flush a no-op
		pass

	def commit(self):
		raise TypeError('Session is read-only')


# Python types corresponding to non-collection types storable in JSON
JSONABLE_SCALARS = (int, float, str, bool, type(None))


class MutableJsonCollection(Mutable):
	"""ABC for SQLAlchemy JSON collection type, supporting mutation tracking

	When nested under another collection, keeps a weak reference to its
	parent so that mutation notifications can be propagated back up to the
	root collection.
	"""

	def __init__(self, parent=None):
		if parent is not None:
			self._parent = weakref.proxy(parent)
		else:
			self._parent = None

	def __getstate__(self):
		"""For pickling"""
		return self.as_builtin()

	def __setstate__(self, state):
		"""For unpickling"""
		self.__init__(state)

	def changed(self):
		super(MutableJsonCollection, self).changed()

		# Call changes() method on parent (if exists)
		if self._parent is not None:
			try:
				self._parent.changed()
			except ReferenceError:
				self._parent = None

	def as_builtin(self):
		"""Return version of collection as builtin Python type

		(for JSON serialization)
		"""
		raise NotImplementedError()

	def _transform_element(self, elem):
		"""Transforms python types into MutableJsonCollection where possible"""

		if isinstance(elem, MutableJsonCollection):
			return elem

		elif isinstance(elem, JSONABLE_SCALARS):
			return elem

		elif isinstance(elem, collections.Mapping):
			return MutableJsonDict(elem, parent=self)

		elif isinstance(elem, collections.Sequence):
			return MutableJsonList(elem, parent=self)

		else:
			raise TypeError('{} is not a JSONable type'.format(type(elem)))

	@classmethod
	def _element_as_builtin(cls, elem):
		"""Converts an element to builtin type"""
		if isinstance(elem, MutableJsonCollection):
			return elem.as_builtin()

		else:
			return elem


class MutableJsonList(MutableJsonCollection, collections.MutableSequence):
	"""List-like object corresponding to JSON array in SQLAlchemy.

	Mutations will be tracked in this object and all nested collections.
	"""

	def __init__(self, sequence, parent=None):
		# Recursively convert any nested collections to MutableJsonCollection
		self._list = list(map(self._transform_element, sequence))

		MutableJsonCollection.__init__(self, parent)

	def __getnewargs__(self):
		"""Picklable tuple of args to __new__ for unpickling."""
		return (self.as_builtin(),)

	def __getitem__(self, index):
		return self._list[index]

	def __setitem__(self, index, value):
		self._list[index] = self._transform_element(value)
		self.changed()

	def __delitem__(self, index):
		del self._list[index]
		self.changed()

	def __iter__(self):
		return iter(self._list)

	def __len__(self):
		return len(self._list)

	def __repr__(self):
		return repr(self._list)

	def insert(self, index, value):
		self._list.insert(index, self._transform_element(value))
		self.changed()

	def as_builtin(self):
		return list(map(self._element_as_builtin, self._list))

	@classmethod
	def coerce(cls, key, value):
		if not isinstance(value, MutableJsonList):
			if isinstance(value, collections.Sequence):
				return MutableJsonList(value)
			else:
				return Mutable.coerce(key, value)
		else:
			return value


class MutableJsonDict(MutableJsonCollection, collections.MutableMapping):
	"""Dict-like object corresponding to JSON object in SQLAlchemy.

	Mutations will be tracked in this object and all nested collections.
	"""

	def __init__(self, mapping, parent=None):
		# Recursively convert any nested collections to MutableJsonCollection
		self._dict = {
			k: self._transform_element(v)
			for k, v in dict(mapping).items()
		}

		MutableJsonCollection.__init__(self, parent)

	def __getnewargs__(self):
		"""Picklable tuple of args to __new__ for unpickling."""
		return (self.as_builtin(),)

	def __getitem__(self, key):
		return self._dict[key]

	def __setitem__(self, key, value):
		if not isinstance(key, str):
			raise TypeError('Key must be string')

		self._dict[key] = self._transform_element(value)
		self.changed()

	def __delitem__(self, key):
		del self._dict[key]
		self.changed()

	def __iter__(self):
		return iter(self._dict)

	def __len__(self):
		return len(self._dict)

	def __repr__(self):
		return repr(self._dict)

	def as_builtin(self):
		return {
			k: self._element_as_builtin(v)
			for k, v in self._dict.items()
		}

	@classmethod
	def coerce(cls, key, value):
		if not isinstance(value, MutableJsonDict):
			if isinstance(value, collections.Mapping):
				return MutableJsonDict(value)
			else:
				return Mutable.coerce(key, value)
		else:
			return value


class MutableJsonCollectionEncoder(json.JSONEncoder):
	"""JSON encoder capable of serializing MutableJsonCollection objects"""

	def default(self, obj):
		if isinstance(obj, MutableJsonCollection):
			return obj.as_builtin()
		else:
			return super(MutableJsonCollectionEncoder, self).default(obj)


class JsonType(TypeDecorator):
	"""SQLA column type for JSON data"""

	impl = sa.String

	def process_bind_param(self, value, dialect):
		if value is not None:
			return json.dumps(value, separators=(',', ':'),
			                  cls=MutableJsonCollectionEncoder)
		else:
			return None

	def process_result_value(self, value, dialect):
		if value is not None:
			json_val = json.loads(value)

			if isinstance(json_val, dict):
				return MutableJsonDict(json_val)
			elif isinstance(json_val, list):
				return MutableJsonList(json_val)
			else:
				return json_val

		else:
			return None


class KeyValueTable(collections.MutableMapping):
	"""Mutable key-value store in a database table.

	Fully implements :class:`collections.MutableMapping`, with all methods
	translating into database queries.

	Parameters
	----------
	engine : sqlalchemy.engine.Engine
		Engine connected to database containing table.
	table : sqlalchemy.Table
		Table storing key-value pairs. Should be created with :meth:`make_table`.

	Attributes
	----------
	engine
		Engine connected to database containing table.
	table
		SQLAlchemy table object containing key-value pairs.
	nullable : bool
		Whether the table allows null values.
	"""

	def __init__(self, engine, table):
		self.engine = engine
		self.table = table
		self.nullable = self.table.c.value.nullable

	@classmethod
	def make_table(cls, name, metadata, *args, nullable=False, **kwargs):
		"""Create an SQLAlchemy table object to use with this class.

		Parameters
		----------
		name : str
			Name of table.
		metadata : sqlalchemy.MetaData
			SQLAlchemy metadata for table.
		\\*args
			Additional positional arguments to :class:`sqlalchemy.Table`.
		nullable : bool
			Whether to make the value column nullable.
		\\**kwargs
			Additional keyword arguments to :class:`sqlalchemy.Table`.

		Returns
		-------
		sqlalchemy.Table
			Table with string primary key column "key" and additional string column "value".
		"""
		return sa.Table(
			name,
			metadata,
			sa.Column('key', sa.String(), primary_key=True),
			sa.Column('value', sa.String(), nullable=nullable),
			*args,
			**kwargs
		)

	@staticmethod
	def _check_key(key):
		"""Check a key argument is valid. Just needs to be a string."""
		if not isinstance(key, str):
			raise TypeError('Key must be string')

	def _check_value(self, value):
		"""Check a value argument is valid."""
		if value is None:
			if not self.nullable:
				raise TypeError('Table does not accept null values')

		elif not isinstance(value, str):
			raise TypeError('Value must be string')

	def __iter__(self):
		rows = self.engine.execute(sa.select([self.table.c.key]))

		for key, in rows:
			yield key

	def __len__(self):
		stmt = sa.select([sa.func.count()]).select_from(self.table)
		return self.engine.execute(stmt).scalar()

	def __contains__(self, key):
		self._check_key(key)

		stmt = sa.select([sa.func.count()])\
			.select_from(self.table)\
			.where(self.table.c.key == key)

		return self.engine.execute(stmt).scalar() == 1

	def __getitem__(self, key):
		self._check_key(key)

		stmt = sa.select([self.table.c.value]).where(self.table.c.key == key)

		row = self.engine.execute(stmt).first()

		if row is None:
			raise KeyError(key)

		else:
			return row[0]

	def __setitem__(self, key, value):
		self._check_key(key)
		self._check_value(value)

		# Try updating an existing row
		update_stmt = self.table.update()\
			.where(self.table.c.key == key)\
			.values(value=value)

		result = self.engine.execute(update_stmt)

		# If does not already exist, insert
		if result.rowcount == 0:
			self.engine.execute(self.table.insert(), key=key, value=value)

	def __delitem__(self, key):
		self._check_key(key)

		stmt = self.table.delete().where(self.table.c.key == key)
		result = self.engine.execute(stmt)

		if result.rowcount == 0:
			raise KeyError(key)

	def clear(self):
		self.engine.execute(self.table.delete())
