"""Base and core classes for storing signature sets."""

from abc import abstractmethod
from typing import Optional, Sequence, Mapping, Any

import numpy as np
from attr import Factory

from midas.util.attr import attrs, attrib


@attrs()
class SignaturesMeta:
	"""Metadata describing a set of k-mer signatures.

	All attributes are optional.

	Attributes
	----------
	id
		Any kind of string ID that can be used to uniquely identify the signature set.
	version
		Version string (ideally PEP 440-compliant).
	name
		Short human-readable name.
	id_attr
		Name of ``Genome`` attribute the IDs correspond to (see :data:`midas.db.models.GENOME_ID_ATTRS`).
		Optional, but signature set cannot be used as a reference for queries without it.
	description
		Human-readable description.
	extra
		Extra arbitrary metadata. Should be a ``dict`` or other mapping which can be converted to JSON.
	"""

	id : Optional[str] = attrib(optional=True, kw_only=True)
	name : Optional[str] = attrib(optional=True, kw_only=True)
	version : Optional[str] = attrib(optional=True, kw_only=True)
	id_attr : Optional[str] = attrib(optional=True, kw_only=True)
	description : Optional[str] = attrib(optional=True, kw_only=True, repr=False)
	extra : Mapping[str, Any] = attrib(default=Factory(dict), kw_only=True, repr=False)


class AbstractSignatureArray(Sequence[np.ndarray]):
	"""
	Abstract base class for types which behave as a (non-mutable) sequence of k-mer signatures
	(k-mer sets in sparse coordinate format).

	Elements should be Numpy arrays with integer data type. Should implement numpy-style advanced
	indexing, see :class:`midas.util.indexing.AdvancedIndexingMixin`. Slicing and advanced indexing
	should return another instance of ``AbstractSignatureArray``.

	Attributes
	----------
	dtype
		Numpy data type of signatures.
	"""
	dtype: np.dtype

	@abstractmethod
	def sizeof(self, index: int) -> int:
		"""Get the size/length of the signature at the given index.

		Should be the case that

		    sigarray.size_of(i) == len(sigarray[i])

		Parameters
		----------
		index
			Index of signature in array.
		"""

	def sizes(self) -> Sequence[int]:
		"""Get the sizes of all signatures in the array."""
		return np.fromiter(map(self.sizeof, range(len(self))))
