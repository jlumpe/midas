"""Run queries against a reference database."""

import numpy as np

from midas.cython import metrics
from midas.kmers import find_kmers_parse, vec_to_coords
from midas.util import kwargs_done, path_str


def sigarray_scores(signature, sigarray, distance=False):
	"""Calculate Jaccard scores between one signature and an array of signatures.

	This internally uses Cython code that runs in parallel over all signatures
	in ``sigarray``.

	:param signature: K-mer signature in coordinate format, increasing sequence
		of integer values.
	:type signature: numpy.ndarray
	:param sigarray: Signature array to calculate scores against.
	:type sigarray: midas.kmers.SignatureArray
	:param bool distance: Return Jaccard distances instead of scores.
	"""

	values = sigarray.values
	bounds = sigarray.bounds.astype(metrics.BOUNDS_DTYPE)

	scores = metrics.jaccard_coords_col(signature, values, bounds)

	return 1 - scores if distance else scores


def find_closest_signatures(query, refarray, *, k=None, distance=False):
	"""Find the closest reference signatures to a query.

	:param query: Single query signature in coordinate format
		(:class:`numpy.ndarray` of increasing integer values) or sequence of
		query signatures (e.g. :class:`midas.kmers.SignatureArray or list).
	:type query: numpy.ndarray
	:param refarray: Array of reference signatures to calculate scores against.
	:type refarray: midas.kmers.SignatureArray
	:param int k: Number of reference signatures to find for each query. If
		None will only find the closest.
	:param bool distance: Report Jaccard distances instead of scores.

	:returns: ``(index, score)`` tuple giving the index of the closest reference
		signature in the array and the score between the query and the reference.
		Both are scalars or arrays of the same shape. If ``query`` is a sequence
		of signatures then the first axis corresponds to the signatures in the
		query. If ``k`` is not None then the last axis corresponds to reference
		matches in order of decreasing similarity.
	:rtype: tuple
	"""

	# Check k
	if k is not None and not (0 < k <= len(refarray)):
		raise ValueError(
			'k must be > 0 and <= the number of reference signatures'
		)

	if isinstance(query, np.ndarray):
		# Single query

		scores = sigarray_scores(query, refarray)

		if k is None:
			indices = np.argmax(scores)

		else:
			indices = np.argsort(scores)[:-k - 1:-1]

		closest_scores = scores[indices]

		return indices, 1 - closest_scores if distance else closest_scores

	else:
		# Assume sequence of queries - e.g. SignatureArray or list

		out_shape = (len(query),) if k is None else (len(query), k)
		indices = np.zeros(out_shape, dtype=int)
		scores = np.zeros(out_shape, dtype=metrics.SCORE_DTYPE)

		# Call self with individual signatures
		for i, signature in enumerate(query):
			indices[i], scores[i] = find_closest_signatures(
				signature,
				refarray,
				k=k,
				distance=distance
			)

		return indices, scores


def _query_files_worker(args):
	"""Thread worker function to calculate signatures from sequence files.

	:param args: Tuple of ``(i, file, format_ kmerspec)``. ``i`` is file index,
		``file`` is file name, ``format_`` is file format for parsing, and
		``kmerspec`` is the :class:`midas.kmers.KmerSpec` for k-mer finding.
	:returns: ``(i, signature)`` tuple
	"""
	i, file, format_, kmerspec = args

	with open(path_str(file)) as fobj:
		vec = find_kmers_parse(kmerspec, fobj, format_)
		return i, vec_to_coords(vec)


def query_files_parallel(query_files, format_, kmerspec, refarray, *,
                         k=None, nthreads=None, **kwargs):
	"""Find closest reference genomes for a set of query sequence files.

	A process pool is set up to parse input files and calculate their
	signatures. Querying is done in the main thread as signatures become
	available, but the Cython function implementation is parallelized as well
	and releases the GIL.
	"""

	from multiprocessing import Pool

	# Get keyword arguments
	distance = kwargs.pop('distance', False)
	progress = kwargs.pop('progress', None)
	threads = kwargs.pop('threads', None)
	kwargs_done(kwargs)

	# Format strings to list
	if isinstance(format_, str):
		format_ = [format_] * len(query_files)

	elif len(format_) != len(query_files):
		raise ValueError('format_ must have same length as files argument')

	# Output arrays
	out_shape = len(query_files) if k is None else (len(query_files), k)
	closest_indices = np.zeros(out_shape, dtype=int)
	closest_scores = np.zeros(out_shape, dtype=metrics.SCORE_DTYPE)

	# Jobs for worker threads
	jobs = [
		(i, file, fmt, kmerspec)
		for i, (file, fmt) in enumerate(zip(query_files, format_))
	]

	with Pool(threads) as pool:

		# Iterate over signatures as they finish
		for i, signature in pool.imap_unordered(_query_files_worker, jobs):

			# Run query with signature
			indices, scores = find_closest_signatures(signature, refarray, k=k,
			                                          distance=distance)

			# Fill in results
			closest_indices[i] = indices
			closest_scores[i] = scores

			# Call progress function
			if progress:
				progress(
					i=i,
					file=query_files[i],
					indices=indices,
					scores=scores,
					signature=signature,
				)

	return closest_indices, closest_scores


def get_genome_by_attr(session, attrname, attrval, *, ref_set=None, force=False):

	from midas.database import models

	if ref_set is None:
		query = session.query(models.Genome).filter_by({attrname: attrval})

	else:
		ref_set_id = ref_set if isinstance(ref_set, int) else ref_set.id
		query = session.query(models.AnnotatedGenome)\
			.filter(models.AnnotatedGenome.reference_set_id == ref_set_id)\
			.filter_by({attrname: attrval})

	if force:
		return query.one()
	else:
		return query.scalar()


def genomes_from_ids(ids, session, *, id_attr='key', ref_set=None, force=False):
	return [
		get_genome_by_attr(session, id_, id_attr, ref_set=ref_set, force=force)
		for id_ in ids
	]
