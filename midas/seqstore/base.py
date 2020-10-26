"""Base classes for sequence file storage."""

from abc import ABCMeta, abstractmethod, abstractproperty
from collections import namedtuple

from midas import ncbi


class BaseSequenceStoreRecord(ncbi.SeqRecordBase):
	"""ABC for record describing a sequence stored in a :class:`.SequenceStore`.

	Inherits all attributes from :class:`midas.ncbi.SeqRecordBase`.

	Attributes
	----------
	store_id
		Unique ID of the sequence within the :class:`.SequenceStore` it is
		stored in. Does not have any meaning outside of this.
	format : str
		String describing the file format the sequence is stored in, e.g.
		"fasta".
	"""

	store_id = abstractproperty()
	format = abstractproperty()


# Basic namedtuple implementation of BaseSequenceStoreRecord
_SequenceStoreRecord = namedtuple(
	'_SequenceStoreRecord',
	('store_id', 'format') + ncbi.SEQ_ID_ATTRS
)
class SequenceStoreRecord(_SequenceStoreRecord, BaseSequenceStoreRecord):
	"""Record describing a sequence stored in a :class:`.SequenceStore`."""


class SequenceIDError(ValueError):
	"""Base class for exceptions raised for bad/conflicting sequence IDs."""


class SequenceNotFound(SequenceIDError):
	"""Raised when a stored sequence with a given ID does not exist."""


class SequenceIDConflict(SequenceIDError):
	"""
	Raised when attempting to store a sequence with the same ID as an existing
	one.
	"""


class SequenceStore(metaclass=ABCMeta):
	"""ABC for an indexed collection of stored sequences."""

	@abstractmethod
	def store(self, src, ids, **kwargs):
		"""Add a new genome sequence to the store.

		Parameters
		----------
		src
			Open file-like object containing the sequence data, or alternatively a string containing
			the path to such a file.
		ids : dict
			NCBI IDs for the sequence. At least one is required. See :func:`midas.ncbi.get_seq_ids`
			for other requirements.
		seq_format : str
			  File format of the sequence. Defaults to "fasta".
		src_compression : str
			  Compression of the source data. Allowable values are ``None`` and
			  "gzip".
		keep_src : bool
			  If a file name is given for ``src``, whether to copy (True) or
			  move (False) the source file. Defaults to True.
		src_mode : str
			  Mode of ``src`` if it is a file-like object - "t" or "b". Only
			  applicable if compression is ``None`` (otherwise assumes binary).

		Returns
		-------
		.SeqStoreRecord
			Record for the added sequence.
		"""
		pass

	@abstractmethod
	def open(self, which):
		"""Get an open file handle/stream to the sequence for a stored genome.

		Parameters
		----------
		which
			Either the store ID of the sequence or its :class:`.BaseSequenceStoreRecord`. See
			:meth:`get_record` for getting records of stored sequences by their NCBI IDs.

		Returns
		-------
			Open file-like object in text mode from which the sequence data can be read.

		Raises
		------
		SequenceNotFoundError
			If no sequence corresponding to ``which`` exists.
		"""
		pass

	@abstractmethod
	def remove(self, which):
		"""Remove a sequence from the store.

		Parameters
		----------
		which
			Either the store ID of the sequence or its :class:`.BaseSequenceStoreRecord`. See
			:meth:`get_record` for getting records of stored sequences by their NCBI IDs.


		Raises
		------
		SequenceNotFoundError
			If no sequence corresponding to ``which`` exists.
		"""
		pass

	@abstractmethod
	def get_record(self, store_id=None, **ids):
		"""Get the record for a stored sequence by store ID or NCBI IDs.

		Parameters
		----------
		store_id
			ID of record specific to this :class:`.SequenceStore`. Can be found in the ``store_id``
			attribute of an existing record. This parameter is mutually exclusive with ``**ids``.
		\\**ids
			NCBI IDs of sequence. Must be valid set of IDs as per :func:`midas.ncbi.get_seq_ids`. If
			multiple IDs are given the returned sequence must match all of them.

		Returns
		-------
		.BaseSequenceStoreRecord
			Record of stored sequence if it exists, otherwise ``None``.
		"""
		pass

	@abstractmethod
	def has(self, **ids):
		"""Check if a sequence with all the given NCBI IDs exists in the store.

		This should be equivalent to ``seqstore.get_record(**ids) is not None``.

		Parameters
		----------
		\\**ids
			NCBI IDs of sequence. Must be valid set of IDs as per :func:`midas.ncbi.get_seq_ids`. If
			multiple IDs are given a sequence must match all of them.

		Returns
		-------
		bool
		"""
		pass

	@abstractmethod
	def has_any(self, **ids):
		"""Check if a sequence with any of the given NCBI IDs exists in the store.

		Parameters
		----------
		\\**ids
			NCBI IDs of sequence. Must be valid set of IDs as per
			:func:`midas.ncbi.get_seq_ids`. Unlike :meth:`has` if multiple IDs
			are given a sequence only needs to match the attributes in a single
			index (see :data:`midas.ncbi.SEQ_IDS`).

		Returns
		-------
		bool
		"""
		pass

	@abstractmethod
	def update_record(self, which, **attrs):
		"""Update the attributes of a sequence record.

		Parameters
		----------
		which
			Either the store ID of the sequence or its :class:`.BaseSequenceStoreRecord`. See
			:meth:`get_record` for getting records of stored sequences by their NCBI IDs.
		\\**attrs
			Attributes of record to modify. Keywords may be ``format`` or any NCBI ID attributes.
			The NCBI IDs of the record after the updates are applied must be valid. Values may be
			``None``.

		Raises
		------
		SequenceNotFoundError
			If no sequence corresponding to ``which`` exists.
		"""
		pass

	@staticmethod
	def _which_arg_id(which):
		"""
		Get the store ID from the "which" argument to many methods, which may be
		the ID itself or a SequenceStoreRecord.
		"""
		return which if isinstance(which, int) else which.store_id
