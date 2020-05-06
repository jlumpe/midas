"""Sequence stores that work with databases.

.. data:: default_metadata

	:class:`sqlalchemy.MetaData` used for :data:`.default_seq_table` and
	:data:`.default_meta_table`.

.. data:: default_seq_table

	Default :class:`sqlalchemy.Table` to use for storing sequence records.

.. data:: default_meta_table

	Default :class:`sqlalchemy.Table` to use for storing metadata.
"""

from abc import abstractmethod
from functools import reduce
import operator

import sqlalchemy as sa

from midas.db.mixins import SeqRecordMixin
from midas.db.sqla import KeyValueTable
from midas import ncbi
from midas.util import kwargs_done
from . import base


_seq_record_cols = {
	name: getattr(SeqRecordMixin, name)
	for name in ncbi.SEQ_ID_ATTRS
}


def make_seq_table(name, metadata, *args, **kwargs):
	"""Make an SQLAlchemy table object for storing sequence records.

	:param str name: Name of table.
	:param metadata: SQLAlchemy metadata for table.
	:type metadata: sqlalchemy.MetaData
	:param \\*args: Additional positional arguments to :class:`sqlalchemy.Table:`,
		for example for specifying additional columns.
	:param \\**kwargs: Additional keyword arguments to :class:`sqlalchemy.Table:`.

	:rtype: sqlalchemy.Table
	"""

	return sa.Table(
		name,
		metadata,
		sa.Column('store_id', sa.Integer(), primary_key=True),
		sa.Column('format', sa.String()),
		sa.Column('compression', sa.String()),
		*(
			sa.Column(name, col.type, unique=col.unique)
			for name, col in _seq_record_cols.items()
		),
		*SeqRecordMixin.__table_args__,
		*args,
		**kwargs
	)


default_metadata = sa.MetaData()
default_seq_table = make_seq_table('sequences', default_metadata)
default_meta_table = KeyValueTable.make_table('metadata', default_metadata)


class DbIndexedSequenceStore(base.SequenceStore):
	"""
	Base class/mixin for SequenceStores which store an index of sequences in a
	database.

	Handles all index/metadata management so child classes just need to handle
	storage of the actual sequence data.

	Constructor forwards all additional arguments to super constructor so this
	shoudl be able to be used with multiple inheritance.

	:Keyword Arguments:

		* **engine** --
		  SQLAlchemy engine connected to database.

		* **seq_table** --
		  SQLAlchemy table containing the sequence records. If None
		  :data:`.default_seq_table` will be used.

		* **meta_table** --
		  SQLAlchemy table containing the metadata key/value pairs. If None
		  :data:`.default_meta_table` will be used. Irrelevant if ``use_meta``
		  is False.

		* **use_meta** (*bool*) --
		  Whether to use a :class:`midas.db.sqla.KeyValueTable` to store
		  metadata in the database.

	.. attribute:: VERSION

		Class attribute equal to ``None``. If a child class overrides this
		attribute with a string, when the constructor is called it will check
		that the metadata contains this value under the "version" key and raise
		an error if it does not. This string should be changed each time the
		internal database schema changes in a way that is incompatible with the
		previous version.

	.. attribute:: _engine

		Private attribute containing SQLAlchemy engine connected to database.

	.. attribute:: _seq_table

		Private attribute containing SQLAlchemy table which stores sequence
		records.

	.. attribute:: _meta

		Private attribute containing metadata in a
		:class:`midas.db.sqla.KeyValueTable` if ``use_meta=True`` was
		passed to the constructor.
	"""

	VERSION = None

	def __init__(self, *args, **kwargs):

		self._engine = kwargs.pop('engine')
		self._seq_table = kwargs.pop('seq_table', default_seq_table)
		meta_table = kwargs.pop('meta_table', default_meta_table)
		use_meta = kwargs.pop('use_meta', False)

		if use_meta:
			self._meta = KeyValueTable(self._engine, meta_table)
		else:
			self._meta = None

		# Check version
		if self._meta is not None and self.VERSION is not None:
			version = self._meta['version']
			if version != self.VERSION:
				raise RuntimeError(
					'Attempting to create {} with version {}, current version '
					'is {}'
					.format(type(self).__name__, version, self.VERSION)
				)

		super().__init__(*args, **kwargs)

	def _make_record(self, store_id, ids, format=None):
		"""Create a SequenceStoreRecord instance."""
		return base.SequenceStoreRecord(
			store_id=store_id,
			format=format,
			**{attr: ids.get(attr, None) for attr in ncbi.SEQ_ID_ATTRS}
		)

	def _row_to_record(self, row):
		"""
		Create a SequenceStoreRecord instance from a row of the sequences table.
		"""
		return self._make_record(row['store_id'], dict(row), format=row['format'])

	def _get_row_by_id(self, store_id):
		"""Retrieve a row from the sequences table by its store id."""
		stmt = self._seq_table.select()\
			.where(self._seq_table.c.store_id == store_id)
		return self._engine.execute(stmt).first()

	def _make_where_stmt(self, ids):
		"""Create arg to an SQLA where statement from sequence IDs."""
		return reduce(
			operator.and_,
			(self._seq_table.c[attr] == val for attr, val in ids.items())
		)

	def _delete_record_from_db(self, store_id):
		"""Delete a record from the database.

		:param int store_id: Store ID of the record.
		:raises midas.seqstore.base.SequenceNotFound: If no record with the
			store ID exists.
		"""

		stmt = self._seq_table.delete()\
			.where(self._seq_table.c.store_id == store_id)
		result = self._engine.execute(stmt)

		if not result.rowcount:
			raise base.SequenceNotFound(store_id)

	def get_record(self, store_id=None, **ids):

		if store_id is not None:
			if ids:
				raise TypeError("Can't pass store_id and additional IDs")

			row = self._get_row_by_id(store_id)

		else:
			ncbi.get_seq_ids(ids, empty_ok=False)

			stmt = self._seq_table.select().where(self._make_where_stmt(ids))
			row = self._engine.execute(stmt).first()

		return self._row_to_record(row) if row is not None else None

	def has(self, **ids):
		ncbi.get_seq_ids(ids, empty_ok=False)

		stmt = self._seq_table.select().where(self._make_where_stmt(ids))
		return self._engine.execute(stmt).scalar() is not None

	def has_any(self, **ids):
		ncbi.get_seq_ids(ids, empty_ok=False)

		for index_keys in ncbi.SEQ_IDS.values():
			try:
				index_ids = {key: ids[key] for key in index_keys}
			except KeyError:
				continue

			if self.has(**index_ids):
				return True

		return False

	def update_record(self, which, **attrs):

		# Split ID attrs from others
		ids = attrs
		extra_attrs = dict()

		try:
			extra_attrs['format'] = ids.pop('format')
		except KeyError:
			pass

		record = self.get_record(self._which_arg_id(which))

		# Validate new set of IDs
		new_ids = {**record.ncbi_ids(flat=True), **ids}
		ncbi.get_seq_ids(new_ids, empty_ok=False, null_ok=True)

		# Update
		stmt = self._seq_table.update()\
			.where(self._seq_table.c.store_id == record.store_id)
		self._engine.execute(stmt, **new_ids, **extra_attrs)

	def store(self, src, ids, **kwargs):

		# Get keyword arguments
		seq_format = kwargs.pop('format', 'fasta')
		src_compression = kwargs.pop('src_compression', None)
		keep_src = kwargs.pop('keep_src', True)
		src_mode = kwargs.pop('src_mode', 't')

		kwargs_done(kwargs)

		# Check arguments
		if src_compression not in (None, 'gzip'):
			raise ValueError('Unknown src compression {}'.format(src_compression))

		# Check ids
		ncbi.get_seq_ids(ids, empty_ok=False)

		# Check duplicate
		if self.has_any(**ids):
			raise base.SequenceIDConflict(ids)

		# Insert the row
		stmt = self._seq_table.insert().values(
			format=seq_format,
			compression='gzip',
			**ids,
		)
		store_id, = self._engine.execute(stmt).inserted_primary_key

		record = self._make_record(store_id, ids, format=seq_format)

		# Try to add the file to the store
		try:
			record_updates = self._store_seq_data(
				record,
				src,
				src_compression=src_compression,
				keep_src=keep_src,
				src_mode=src_mode,
			)

		except:
			# Expect possible OSError, but regardless of type we want to remove
			# the inserted record
			self._delete_record_from_db(store_id)

			# Allow the error to propagate
			raise

		# Update record if needed
		if record_updates:
			stmt = self._seq_table.update()\
				.where(self._seq_table.c.store_id == store_id)\
				.values(**record_updates)
			self._engine.execute(stmt)

		# Return the record
		return self._make_record(store_id, ids, format=seq_format)

	def open(self, which):

		store_id = self._which_arg_id(which)
		row = self._get_row_by_id(store_id)

		if row is None:
			raise base.SequenceNotFound(store_id)

		return self._open_seq_data(row)

	def remove(self, which):
		store_id = self._which_arg_id(which)

		row = self._get_row_by_id(store_id)

		self._delete_record_from_db(store_id)

		try:
			self._remove_seq_data(row)
		except FileNotFoundError:
			pass

	@abstractmethod
	def _store_seq_data(self, record, src, src_compression, keep_src, src_mode):
		"""Store sequence data after its metadata has been added to the index.

		Abstract method to be implemented by subclass.

		Should an exception be raised from this method the sequence metadata
		will be removed from the database index during cleanup before the
		exception is propagated.

		:param src: Source data as file-like object or file name. See ``src``
			argument to :meth:`store`.
		:param record: Record for the added sequence, containing store ID and
			NCBI IDs, as well as sequence format.
		:type record: midas.seqstore.base.BaseSequenceStoreRecord
		:param str src_compression: See :meth:`store`.
		:param bool keep_src: See :meth:`store`.
		:param str src_mode: See :meth:`store`.

		:returns: ``None``, or optionally a mapping of updates to be applied to
			the sequence's database row keyed by column name. Used to store
			additional data about sequence specific to subclass implementation
			(e.g., name of file the sequence was stored as). Should only update
			additional columns defined in the subclass, no columns in
			:data:`default_seq_table` should be changed.
		"""

	@abstractmethod
	def _open_seq_data(self, row):
		"""Get open stream to stored sequence data.

		Abstract method to be implemented by subclass.

		:param row: Contents of database row for sequence.
		:returns: Open file-like object for reading sequence data.
		"""

	@abstractmethod
	def _remove_seq_data(self, row):
		"""Remove stored sequence data.

		Abstract method to be implemented by subclass.

		:param row: Contents of database row for sequence.
		"""
