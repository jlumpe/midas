"""Store sequences in the file system."""

import os
import shutil
import gzip

import sqlalchemy as sa

from midas.util import SubPath, sanitize_filename
from midas.database.sqla import KeyValueTable
from .dbstore import DbIndexedSequenceStore, make_seq_table


metadata = sa.MetaData()

file_table = make_seq_table(
	'sequences',
	metadata,
	sa.Column('filename', sa.String(), unique=True),
)

meta_table = KeyValueTable.make_table('metadata', metadata)


def make_filename(record, ext=''):
	"""Create a unique and descriptive name for a sequence file.

	:param record: Record for sequence containing store_id and NCBI IDs.
	:type record: midas.seqstore.base.BaseSequenceStoreRecord
	:param str ext: File extension, including dot.

	:rtype: str
	"""

	fname = str(record.store_id)

	if record.genbank_acc is not None:
		fname += '-GB_ACC-' + sanitize_filename(record.genbank_acc)

	elif record.refseq_acc is not None:
		fname += '-RS_ACC-' + sanitize_filename(record.refseq_acc)

	elif record.entrez_db is not None:
		fname += '-ENTREZ-{}_{}'.format(record.entrez_db, record.entrez_id)

	fname += ext

	return fname


class FileSequenceStore(DbIndexedSequenceStore):
	"""An indexed collection of sequences stored in the file system.

	All sequence files as well as index data are stored in a single directory.

	:param str path: Path to existing ``FileSequenceStore`` directory.

	.. attribute:: root_dir

		Path to root directory of sequence store.
	"""

	VERSION = '1.0'

	seq_dir = SubPath('sequences')
	db_path = SubPath('db')

	def __init__(self, path):

		super().__init__(
			engine=self._make_engine(path),
			seq_table=file_table,
			meta_table=meta_table,
			use_meta=True,
		)

		self.root_dir = os.path.abspath(path)

	@classmethod
	def _make_engine(cls, root_path):
		"""Create the SQLA engine for a file sequence store from its path."""
		return sa.create_engine('sqlite:///' + cls.db_path(root_path))

	@classmethod
	def create(cls, path):
		"""Create a new :class:`.FileSequenceStore`.

		:param str path: Path to sequence store directory. Must not already
			exist.
		:rtype: .FileSequenceStore
		"""

		# Create directories
		os.makedirs(path)
		os.mkdir(cls.seq_dir(path))

		# Create engine and initialize tables
		engine = cls._make_engine(path)
		metadata.create_all(engine)

		# Add the current version to the metadata table
		meta = KeyValueTable(engine, meta_table)
		meta['version'] = cls.VERSION

		return cls(path)

	@classmethod
	def get_version(cls, root_path):
		"""Get the version string of an existing FileSequenceStore.

		:param str root_path: Path of store's root directory.
		:rtype: str
		"""
		engine = cls._make_engine(root_path)
		meta = KeyValueTable(engine, meta_table)
		return meta['version']

	def _store_seq_data(self, record, src, src_compression, keep_src, src_mode):

		fname = make_filename(record, ext='.fasta.gz')

		# Get destination path
		dest_path = os.path.join(self.seq_dir, fname)

		# Remove any existing file at the destination path. This shouldn't
		# happen but it's possible the store was corrupted when an error or
		# version upgrade occurred.
		if os.path.isfile(dest_path):
			os.remove(dest_path)

		if isinstance(src, str):
			# If passed a string, assume file path

			if src_compression == 'gzip':
				# Gzip compression can copy/move at file level
				if keep_src:
					shutil.copyfile(src, dest_path)
				else:
					shutil.move(src, dest_path)

			else:
				# Will need to convert to gzip
				with open(src, 'rb') as src_fh:
					with gzip.open(dest_path, 'wb') as dest_fh:
						shutil.copyfileobj(src_fh, dest_fh)

				# Remove original
				if not keep_src:
					os.remove(src)

		else:
			# Otherwise assume file-like object

			# Gzip compression can write directly, otherwise compress
			if src_compression == 'gzip':
				dest_fh = open(dest_path, 'wb')
			else:
				dest_fh = gzip.open(
					dest_path,
					'wb' if 'b' in src_mode else 'wt'
				)

			with dest_fh:
				shutil.copyfileobj(src, dest_fh)

		# Return updates to row
		return dict(filename=fname)

	def _open_seq_data(self, row):
		fpath = os.path.join(self.seq_dir, row['filename'])
		return gzip.open(fpath, 'rt')

	def _remove_seq_data(self, row):
		os.remove(os.path.join(self.seq_dir, row['filename']))
