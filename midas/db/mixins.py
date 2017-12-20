"""Mixin classes for SQLA models."""

from sqlalchemy import Column, Integer, String
from sqlalchemy import UniqueConstraint
from sqlalchemy.ext.declarative import declared_attr

from .sqla import MutableJsonCollection
from pydatatypes import Jsonable, JsonConstructible
from midas import ncbi


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


@Jsonable.register
@JsonConstructible.register
class JsonableMixin:
	"""Mixin that allows model instances to be converted to/from JSON.

	Subclasses must have a class attribute ``__json_attrs__`` which is
	a sequence of field names which are converted to JSON.
	"""

	def to_json(self):
		"""Convert to a value serializable as JSON.

		:returns: Dictionary which can be passed to :func:`json.dump`.
		:rtype: dict
		"""

		data = dict()

		for name in self.__json_attrs__:
			value = getattr(self, name)
			if isinstance(value, MutableJsonCollection):
				value = value.as_builtin()

			data[name] = value

		return data

	@classmethod
	def from_json(cls, data):
		"""Create an instance of the class from JSON data.

		:param dict data: JSON object data as returned by :func:`json.load`.
		:returns: Model instance.
		"""
		return cls(**{
			name: value for name, value in data.items()
			if name in cls.__json_attrs__
		})

	def update_from_json(self, data):
		"""Updates attributes from parsed JSON dict.

		:param dict data: JSON object data as returned by :func:`json.load`.
		"""
		for name, value in data.items():
			if name in self.__json_attrs__:

				if isinstance(value, MutableJsonCollection):
					value = value.as_builtin()

				setattr(self, name, value)


class SeqRecordMixin(ncbi.SeqRecordBase):
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

	@classmethod
	def get_ncbi_id_filter(cls, *args, **kwargs):
		"""Get an expression to filter model isntances based on NCBI sequence IDs.

		If values for multiple sequence IDs are given, their sub-expressions
		will combined via logical AND.

		:param \\*args: Name of NCBI sequence ID followed by its attribute
			values (see :data:`midas.ncbi.SEQ_IDS`). Mututally exclusive with
			``**kwargs``.
		:param \\**kwargs: Valid set of NCBI sequence ID attribute values as
			keyword arguments (see :func:`midas.ncbi.get_seq_ids`).
			Mututally exclusive with``*args``.

		:returns: Expression indicating that this model's NCBI ID attributes
			match the given values.
		:rtype: sqlalchemy.sql.elements.ClauseElement
		"""

		ids = ncbi.parse_seq_id_args(args, kwargs, multiple=True, empty_ok=False)

		# Create the expression
		exp = True  # Identity for logical AND

		for id_name, id_vals in ids.items():

			# Create sub-expression for this ID
			subexp = True  # Identity for logical AND

			for attrname, val in zip(ncbi.SEQ_IDS[id_name], id_vals):
				subexp = (getattr(cls, attrname) == val) & subexp

			# AND with full expression
			exp = subexp & exp

		return exp

	@classmethod
	def by_ncbi_id(cls, session, *args, **kwargs):
		"""Get a query object on the class filtering by NCBI sequence IDs.

		If values for multiple sequence IDs are given, their sub-expressions
		will combined via logical AND (i.e., query results must match all of
		them).

		:param session: SQLAlchemy session object to create query with.
		:type session: sqlalchemy.orm.session.Session
		:param \\*args: Name of NCBI sequence ID followed by its attribute
			values (see :data:`midas.ncbi.SEQ_IDS`). Mututally exclusive with
			``**kwargs``.
		:param \\**kwargs: Valid set of NCBI sequence ID attribute values as
			keyword arguments (see :func:`midas.ncbi.get_seq_ids`).
			Mututally exclusive with``*args``.

		:returns: SQLAlchemy query on model with ID filter applied.
		:rtype: sqlalchemy.orm.query.Query
		"""
		exp = cls.get_ncbi_id_filter(*args, **kwargs)
		return session.query(cls).filter(exp)
