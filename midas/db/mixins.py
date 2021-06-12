"""Mixin classes for SQLA models."""

from sqlalchemy import Column, Integer, String
from sqlalchemy import UniqueConstraint
from sqlalchemy.ext.declarative import declared_attr

from .sqla import MutableJsonCollection
from midas import ncbi


class KeyMixin:
	"""Mixin that defines key/version columns.

	Attributes
	----------
	key : str
		Intended to be a universally unique key that can be used to identify objects across
		databases on different systems. This is primarily intended to be used for distributing
		database updates. It can be any arbitrary string but the recommended format is a
		filepath-like structure separated by forward slashes. This results in a hierarchical
		format that supports using namespaces to avoid key conflicts. Example:
		``'ncbi/assembly/GCF_00000000.0'``, which corresponds to a specific genome stored in the
		NCBI assembly database.
	version : str
		Version of the keyed object according to whatever source defined the key. Used to
		determine when the metadata needs to be updated. Should be in the format defined by `PEP
		440 <https://www.python.org/dev/peps/pep-0440/>`_.
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


class SeqRecordMixin(ncbi.SeqRecordBase):
	"""
	Mixin for models which describe a specific sequence record in an NCBI
	database.
	"""

	@declared_attr
	def __table_args__(cls):
		return (
			UniqueConstraint('entrez_db', 'entrez_id'),
		)

	entrez_db = Column(String())
	entrez_id = Column(Integer())
	genbank_acc = Column(String(), unique=True)
	refseq_acc = Column(String(), unique=True)

	@classmethod
	def get_ncbi_id_filter(cls, *args, **kwargs):
		"""Get an expression to filter model instances based on NCBI sequence IDs.

		If values for multiple sequence IDs are given, their sub-expressions
		will combined via logical AND.

		Parameters
		----------
		\\*args
			Name of NCBI sequence ID followed by its attribute values (see :data:`midas.ncbi.SEQ_IDS`).
			Mutually exclusive with``**kwargs``.
		\\**kwargs
			Valid set of NCBI sequence ID attribute values as keyword arguments (see :func:`midas.ncbi.get_seq_ids`).
			Mututally exclusive with ``*args``.

		Returns
		-------
		sqlalchemy.sql.elements.ClauseElement
			Expression indicating whether this model's NCBI ID attributes match the given values.
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

		Parameters
		----------
		session : sqlalchemy.orm.session.Session
			SQLAlchemy session object to create query with.
		\\*args
			Name of NCBI sequence ID followed by its attribute values
			(see :data:`midas.ncbi.SEQ_IDS`). Mutually exclusive with ``**kwargs``.
		\\**kwargs
			Valid set of NCBI sequence ID attribute values as keyword arguments
			(see :func:`midas.ncbi.get_seq_ids`). Mutually exclusive with``*args``.

		Returns
		-------
		sqlalchemy.orm.query.Query
			SQLAlchemy query on model with ID filter applied.
		"""
		exp = cls.get_ncbi_id_filter(*args, **kwargs)
		return session.query(cls).filter(exp)
