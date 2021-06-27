from abc import abstractmethod
from typing import Sequence

import numpy as np


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
