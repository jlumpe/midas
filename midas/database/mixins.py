"""Mixin classes for SQLA models."""

import datetime

from sqlalchemy import Column, Integer, String, DateTime
from sqlalchemy import UniqueConstraint
from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy import event

from .sqla import MutableJsonCollection
from midas.ncbi import SeqRecordBase



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
	def pass_update_time_callback(cls, mapper, connection, instance):
		now = datetime.datetime.utcnow()
		instance.updated_at = now

	@classmethod
	def __declare_last__(cls):
		"""Called after mapper configured, register listeners"""
		event.listen(cls, 'before_insert', cls._insert_time_callback)
		event.listen(cls, 'before_update', cls._update_time_callback)


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
		"""Creates from parsed JSON dict"""
		return cls(**{
			name: value for name, value in json_data.items()
			if name in cls.__json_attrs__
		})

	def update_from_json(self, json_data):
		"""Updates attributes from parsed JSON dict"""
		for name, value in json_data.items():
			if isinstance(value, MutableJsonCollection):
				value = value.as_builtin()

			if name in self.__json_attrs__:
				setattr(self, name, value)


class SeqRecordMixin(SeqRecordBase):
	"""
	Mixin for models which describe a specific sequence record in an NCBI
	database.
	"""

	__table_args__ = (
		UniqueConstraint('entrez_db', 'entrez_id'),
	)

	entrez_db = Column(String())
	entrez_id = Column(Integer())
	genbank_acc = Column(String(), unique=True)
	refseq_acc = Column(String(), unique=True)
