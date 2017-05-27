"""Mixin classes for SQLA models."""

from sqlalchemy import Column, Integer, String
from sqlalchemy import UniqueConstraint
from sqlalchemy.ext.declarative import declared_attr

from .sqla import MutableJsonCollection
from midas.ncbi import SeqRecordBase


class KeyMixin:
	"""Mixin that defines key/version columns.

	.. attribute:: key

		Intended to be a universally unique key that can be used to identify
		objects across databases on different systems. This is primarily
		intended to be used for distributing database updates. It can be any
		arbitrary string but the recommended format is a filepath-like
		structure separated by forward slashes. This results in a hierarchical
		format that supports using namespaces to avoid key conflicts. Example:
		``'ncbi/assembly/GCF_00000000.0'``, which corresponds to a specific
		genome stored in the Genbank assembly database.

	.. attribute:: version

		Version of the keyed object according to whatever source defined the key.
		Used to determine when the	metadata needs to be updated. Should be in
		the format defined by
		`PEP 440 <https://www.python.org/dev/peps/pep-0440/>`_.
	"""
	key = Column(String(), index=True)
	version = Column(String())

	@declared_attr
	def __table_args__(cls):
		return (
			UniqueConstraint('key', 'version'),
		)

	@classmethod
	def by_key(cls, session, key, version=None):
		query = session.query(cls).filter_by(key=key)

		if version is None:
			return query.order_by(cls.version.desc()).first()
		else:
			return query.filter_by(version=version).scalar()


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
