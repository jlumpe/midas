"""Custom types and other utilities for SQLAlchemy"""

import datetime
import json
import collections
import weakref

from sqlalchemy import Column, Integer, String, DateTime
from sqlalchemy.event import listen
from sqlalchemy.types import TypeDecorator
from sqlalchemy.ext.mutable import Mutable
from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy.orm.session import Session


class ReadOnlySession(Session):
	"""Session class that doesn't allow flushing/committing"""

	def _flush(self, *args, **kwargs):
		raise RuntimeError('This sessison is read-only')


class VersionedMixin:
	"""Mixin for models that implement a version counter"""

	_version_id = Column('version_id', Integer(), nullable=False)

	@declared_attr
	def __mapper_args__(cls):
		return dict(version_id_col=cls._version_id)


class TrackChangesMixin:
	"""Mixin for SQLAlchemy models that tracks when updates are made"""

	created_at = Column(DateTime())
	updated_at = Column(DateTime())

	@classmethod
	def _insert_time_callback(cls, mapper, connection, instance):
		now = datetime.datetime.utcnow()
		instance.created_at = now
		instance.updated_at = now

	@classmethod
	def _update_time_callback(cls, mapper, connection, instance):
		now = datetime.datetime.utcnow()
		instance.updated_at = now

	@classmethod
	def __declare_last__(cls):
		"""Called after mapper configured, register listeners"""
		listen(cls, 'before_insert', cls._insert_time_callback)
		listen(cls, 'before_update', cls._update_time_callback)


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

	def changed(self):
		super(MutableJsonCollection, self).changed()

		# Call changes() method on parent (if exists)
		if self._parent is not None:
			try:
				self._parent.changed()
			except weakref.ReferenceError:
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
			k: self._transform_element(v) for k, v
			in dict(mapping).items()
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
		return {k: self._element_as_builtin(v) for k, v
		        in self._dict.items()}

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

	impl = String

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


class JsonableMixin:
	"""Mixin that allows model instances to be converted to/from JSON"""

	def to_json(self):
		"""Converts to JSONable dict"""

		json_data = dict()

		for name in self.__json_attrs__:
			value = getattr(self, name)
			if isinstance(value, MutableJsonCollection):
				value = value.as_builtin()

			json_data[name] = value

		return json_data

	@classmethod
	def from_json(cls, json_data):
		"""Create instance from parsed JSON dict.

		Calls the class constructor with all key/value pairs in ``json_data``
		as keyword arguments, where the key is in :attr:`__json_attrs__`.

		:param dict json_data: JSON-like dict as returned from
			:func:`json.load`.
		:returns: Instance of class.
		"""
		return cls(**{
			name: value for name, value in json_data.items()
			if name in cls.__json_attrs__
		})

	def update_from_json(self, json_data):
		"""Update attributes from parsed JSON dict.

		:param dict json_data: JSON-like dict as returned from
			:func:`json.load`.
		"""
		for name, value in json_data.items():
			if isinstance(value, MutableJsonCollection):
				value = value.as_builtin()

			if name in self.__json_attrs__:
				setattr(self, name, value)
