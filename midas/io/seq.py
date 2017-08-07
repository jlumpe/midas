"""Read and parse sequence files and calculate their k-mer signatures."""

import numpy as np
from Bio import SeqIO

from midas.kmers import find_kmers, vec_to_coords


def find_kmers_parse(kspec, data, format, out=None, coords=False):
	"""Parse sequence data with Bio.Seq.parse() and find k-mers.

	:param kspec: Spec for k-mer search.
	:type kspec: .KmerSpec
	:param data: Stream with sequence data. Readable file-like object in text
		mode.
	:param str format: Squence file format, as interpreted by
		:func:`Bio.SeqIO.parse`.
	:param out: Existing numpy array to write output to. Should be of length
		``kspec.idx_len``. If given the same array will be returned.
	:type out: numpy.ndarray
	:param bool coords: If True return k-mers in coordinate rather than vector
		format.

	:returns: If coords is False, returns boolean K-mer vector (same array as
		``out`` if it was given). If coords is True returns k-mers in coordinate
		format (dtype will match :func:`midas.kmers.vec_to_coords`).
	:rtype: numpy.ndarray
	"""

	if out is None:
		out = np.zeros(kspec.idx_len, dtype=bool)

	for record in SeqIO.parse(data, format):
		find_kmers(kspec, record.seq, out=out)

	if coords:
		return vec_to_coords(out)
	else:
		return out
