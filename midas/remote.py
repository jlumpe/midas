"""Utilities for downloading remote data, e.g. from GenBank"""

import io
import shutil
import contextlib
from urllib.request import urlopen
from urllib.error import URLError
from threading import Lock

from .parallel import AbortableThreadPoolExecutor


class URLRetryError(Exception):
	"""Raised by urlopen_retry when all attempts fail.

	:attr attempt_errors:
		List of errors encountered on all attempts.
	"""

	def __init__(self, message, attempt_errors):
		attempt_errors = tuple(attempt_errors)
		super().__init__(message, attempt_errors)
		self.attempt_errors = attempt_errors


def urlopen_retry(url, attempts=3, **kwargs):
	"""Calls ``urlopen()`` but retries a given number of times after failure.

	:param str url:
	:param int attempts: Number of attempts to make.
	:param kwargs: Passed to :func:`urllib.request.urlopen`.

	:returns: Readable file-like object.

	:raises URLRetryError:
		If maximum number of attempts exceeded. Context is the error encountered
		on the most recent attempt.
	"""

	errors = []
	attempts_made = 0

	while attempts_made < attempts:

		try:
			return urlopen(url, **kwargs)

		except URLError as exc:
			errors.append(exc)

		attempts_made += 1

	msg = 'Aborted after {} attempts'.format(attempts_made)
	raise URLRetryError(msg, errors) from errors[-1]


def download_sequence(db, genome_id, url, db_lock=None, store_opts=dict(),
                      aborted=None, **kwargs):
	"""Download sequence for genome and store it in a database.

	Used as thread worker for :func:`.download_genome_sequences`.
	"""

	# Download as bytes
	# TODO - throttle after error 530?
	seq_data = urlopen_retry(url, **kwargs)

	# Place in buffer
	buf = io.BytesIO()
	shutil.copyfileobj(seq_data, buf)
	buf.seek(0)

	# If in gzip format, check it is readable
	if store_opts.get('src_compression', None) == 'gzip':
		try:
			while buf.read(2 ** 16):
				pass
		except (EOFError, IOError) as exc:
			exc2 = RuntimeError('Gzipped data appears corrupt: {}'.format(exc))
			raise exc2 from exc

		# Be kind, rewind
		buf.seek(0)

	# Store
	with contextlib.ExitStack() as exitstack:
		if db_lock:
			exitstack.enter_context(db_lock)

		if aborted is None or not aborted():
			db.store_sequence(genome_id, buf, src_mode='b', **store_opts)
			return True

		return False


def download_sequences_parallel(db, items, **kwargs):
	"""Download a large set of sequences in parallel and store in a database.

	This function will regularly be downloading tens of thousands of sequences
	at a time over several hours, so it needs to be able to exit gracefully.
	This is difficult with threads. It will attempt to abort all pending jobs
	and wait until completion if any errors are encountered in the main thread.
	This includes the first KeyboardInterrupt, which will print a warning
	message, but subsequent ones will exit immediately.

	:type db: midas.database.base.AbstractDatabase
	:param items: Collection of ``(genome, url)`` or ``(genome, url, store_opts)``
		tuples representing sequences to download.
	:param int max_workers: Maximum number of worker threads to use
	:param dict store_opts: Keyword arguments to
		:meth:`midas.database.base.AbstractDatabase.store_sequence` to use for
		all items.
	:param callable callback:
		Function with signature ``(genome, exc)`` that is called (in the main
		thread) for each item that completes. If the download failed the exception
		which caused the failure will be passed as ``exc``, otherwise this
		argument will be None. This function may return ``False`` to signal that
		the worker pool should abort all future jobs and shut down.

	Additional keyword arguments are passed to :func:`urlopen_retry`.

	:returns: False if the process was aborted before it could finish, True otherwise.
	:rtype: bool
	"""

	max_workers = kwargs.pop('max_workers', None)
	common_store_opts = kwargs.pop('store_opts', dict())
	callback = kwargs.pop('callback', None)

	with AbortableThreadPoolExecutor(max_workers) as executor:

		# Dictionary to get genomes by their future when finished
		genomes_by_future = dict()

		# Lock around database usage
		db_lock = Lock()

		# Submit the tasks
		for item in items:

			if len(item) == 2:
				genome, url = item
				store_opts = common_store_opts
			elif len(item) == 3:
				genome, url, store_opts = item
			else:
				raise ValueError('Download items shoud be 2- or 3-tuples')

			future = executor.submit(
				download_sequence,
				db,
				genome.id,
				url,
				db_lock=db_lock,
				store_opts=store_opts,
				aborted=executor.is_aborted,
				**kwargs
			)
			genomes_by_future[future] = genome

		# Iterator over futures as they complete
		# We only actually care about the result if we need to send it to
		# the callback, but we still want to iterate over all of them to ensure
		# the function does not return until all futures are completed
		# (successfully or unsuccessfully).
		for future in executor.as_completed():

			genome = genomes_by_future[future]

			# Get result for callback
			if callback is not None and not future.cancelled():

				# Getting the result will raise any exception that
				# occured in the worker thread
				try:
					genome_added = future.result()

				# Call callback with exception
				except Exception as exc:
					callback_return = callback(genome, exc)

				else:
					# Worker returns True if not aborted
					if genome_added:
						# Worker finished, call callback
						callback_return = callback(genome, None)
					else:
						# Abort has already happened, skip
						callback_return = None

				# If False returned from callback, abort
				if callback_return is False:
					executor.abort()

	return not executor.is_aborted()
