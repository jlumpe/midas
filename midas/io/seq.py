"""Read and parse sequence files and calculate their k-mer signatures."""

from pathlib import Path

import numpy as np
from Bio import SeqIO

from pydatatypes import dataclass, field
from midas.kmers import find_kmers, vec_to_coords
from .util import open_compressed, ClosingIterator


@dataclass(frozen=True, slots=True)
class SeqFileInfo:
	"""A reference to a DNA sequence file stored in the file system.

	Contains all the information needed to read and parse the file.

	.. attribute:: path

		Path to the file, as :class:`pathlib.Path`.

	.. attribute:: fmt

		String describing the file format, as interpreted by
		:func:`Bio.SeqIO.parse`. E.g. ``'fasta'``.

	.. attribute:: compression

		String describing compression method of the file, e.g. ``'gzip'``. None
		means no compression. See :func:`midas.io.util.open_compressed`.

	:param path: Value of :attr:`path` attribute. May be string or path-like
		object.
	:param str fmt: Value of :attr:`fmt` attribute.
	:param str compression: Value of :attr:`compression` attribute.
	"""

	path = field(Path, convert=Path)
	fmt = field(str)
	compression = field(str, optional=True)

	def open(self, mode='r', **kwargs):
		"""
		Open a stream to the file, with compression/decompression applied
		transparently.

		:param str mode: Same as equivalent argument to the built-in :func:open`.
			Some modes may not be supported by all compression types.
		:param \\**kwargs: Additional text mode specific keyword arguments to
			pass to opener. Equivlent to the following arguments of the built-in
			:func:`open`: ``encoding``, ``errors``, and ``newlines``. May not be
			supported by all compression types.
		:returns: Stream to file in given mode.
		"""
		return open_compressed(self.compression, self.path, mode, **kwargs)

	def parse(self, **kwargs):
		"""Open the file and laziy parse its contents.

		Returns iterator over sequence data in file. File is parsed lazily,
		and so must be kept open. The returned iterator is of type
		:class:`midas.io.util.ClosingIterator` so it will close the file stream
		automatically when it finishes. It may also be used as a context manager
		that closes the stream on exit. You may also close the stream explicitly
		using the iterator's ``close`` method.

		:param \\**kwargs: Keyword arguments to :meth:`open`.
		:returns: Iterator yielding :class:`Bio.SeqIO.SeqRecord` instances for
			each sequence in the file.
		:rtype: midas.io.util.ClosingIterator
		"""

		fobj = self.open('rt', **kwargs)

		try:
			records = SeqIO.parse(fobj, self.fmt)
			return ClosingIterator(records, fobj)

		except:
			fobj.close()
			raise

	def absolute(self):
		"""Make a copy of the instance with an absolute path.

		:rtype: .SeqFileInfo
		"""
		if self.path.is_absolute():
			return self
		else:
			return SeqFileInfo(self.path.absolute(), self.fmt, self.compression)

	@classmethod
	def from_paths(cls, paths, fmt, compression=None):
		"""
		Create many instances at once from a collection of paths and a single
		format and compression type.

		:param paths: Collection of paths as strings or path-like objects.d
		:param str format: Sequence file format of files.
		:param str compression: Compression method of files.

		:rtype: list[.SeqFileInfo]
		"""
		return [cls(path, fmt, compression) for path in paths]


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


class FileSignatureCalculator:
	"""
	Wrapper around a process pool that can parse sequence files and calculate
	their k-mer signatures concurrently.

	Class can be used as a context manager which will shut down the pool
	(calling the :meth:`multiprocessing.Pool.terminate` method) upon exit.

	.. attribute:: pool

		The process pool in use.

	:param int processes: Number of processes to use in the pool. If None will
		use machine's CPU count.
	:param bool use_threads: If True will parse the files in separate threads
		within the same process, using :class:`multiprocessing.dummy.Pool`
		instead of :class:`multiprocessing.Pool`.
	"""

	def __init__(self, processes=None, use_threads=False):

		if use_threads:
			from multiprocessing.dummy import Pool
		else:
			from multiprocessing import Pool

		self.pool = Pool(processes)

	def calc_signatures(self, kmerspec, files, fmt=None, compression=None, *,
	                    ordered=False):
		"""
		Parse a set of sequence files and calculate their signatures in parallel.

		:param kmerspec: K-mer spec to use for calculating signatures.
		:type kmerspec: midas.kmers.kmerspec
		:param files: Sequence or collection of files to parse. Items may be
			:class:`.SeqFileInfo` or simply file paths (as strings or path-like
			objects), in which case ``fmt`` should be specified.
		:param str fmt: Format of sequence files, if ``files`` contains
			only file paths instead of :class:`.SeqFileInfo` instances.
			Passed to :func:`Bio.SeqIO.parse`.
		:param str compression: Compresssion format of sequence files if
			`files`` contains only file paths instead of :class:`.SeqFileInfo`
			instances. Passed to :func:`midas.io.util.open_compressed`. None
			means files are uncompressed.
		:param bool ordered: If True the returned iterator will iterate over
			the results of each file in the order they were in the ``files``
			argument. This could end up being a bit slower.

		:returns: Iterator yielding  ``(i, signature)`` tuples where ``i`` is
			the index of the file in ``files`` and ``signature`` is the file's
			signature in coordinate format.
		:rtype: tuple[int, numpy.ndarray]
		"""

		# Format files to SeqFileInfo
		files = files[:]

		for i, val in enumerate(files):
			if not isinstance(val, SeqFileInfo):
				if fmt is None:
					raise TypeError(
						'Must specifiy format in "fmt" if "files" conatins '
						'paths instead of SeqFileInfo instances'
					)

				files[i] = SeqFileInfo(val, fmt, compression)

		# Create jobs and sent to process pool
		jobs = [(i, file, kmerspec) for i, file in enumerate(files)]

		if ordered:
			return self.pool.imap(self.worker, jobs)
		else:
			return self.pool.imap_unordered(self.worker, jobs)

	def __enter__(self):
		return self

	def __exit__(self, *args):
		self.pool.terminate()

	@staticmethod
	def worker(args):
		"""Thread worker function to calculate signatures from sequence files.

		:param args: Tuple of ``(i, info kmerspec)``. ``i`` is file index,
			``file`` is an instance of :class:`SeqFileInfo`, and ``kmerspec`` is
			the :class:`midas.kmers.KmerSpec` for k-mer finding.
		:returns: ``(i, signature)`` tuple
		"""
		i, file, kmerspec = args

		with file.open() as fobj:
			vec = find_kmers_parse(kmerspec, fobj, file.fmt)
			return i, vec_to_coords(vec)
