"""Utilities for downloading remote data, e.g. from GenBank"""

import sys
import io
import shutil
import contextlib
from urllib.request import urlopen
from urllib.error import URLError

from concurrent.futures import ThreadPoolExecutor, as_completed
import threading


class URLRetryError(Exception):
	"""Raised by urlopen_retry when all attempts fail"""

	def __init__(self, message, attempt_errors):
		attempt_errors = tuple(attempt_errors)
		super().__init__(message, attempt_errors)
		self.attempt_errors = attempt_errors


def urlopen_retry(url, attempts=3, **kwargs):
	"""Calls urlopen() but retries a given number of times after failure"""

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


def download_sequence(db, genome, url, db_lock=None, store_opts=dict(),
                      aborted=None, **kwargs):
	"""Download sequence for genome and store it

	Used as thread worker for download_genome_sequences.
	"""

	# Download as bytes
	# TODO - throttle after error 530?
	seq_data = urlopen_retry(url, **kwargs)

	# Place in buffer
	buf = io.BytesIO()
	shutil.copyfileobj(seq_data, buf)
	buf.seek(0)

	# Store
	with contextlib.ExitStack() as exitstack:
		if db_lock:
			exitstack.enter_context(db_lock)

		if aborted is None or not aborted():
			db.store_sequence(genome, buf, src_mode='b', **store_opts)
			return True
		else:
			return False


def download_sequences_parallel(db, items, **kwargs):
	"""Download and store set of sequences in parallel

	This function will regularly be downloading tens of thousands of sequences
	at a time over several hours, so it needs to be able to exit gracefully.
	This is difficult with threads. It will attempt to abort all pending jobs
	and wait unti completion if any errors are encountered in the main thread.
	This includes the first KeyboardInterrupt, which will print a warning
	message, but subsequent ones will exit immediately.
	"""

	max_workers = kwargs.pop('max_workers', None)
	common_store_opts = kwargs.pop('store_opts', dict())
	callback = kwargs.pop('callback', None)

	with ThreadPoolExecutor(max_workers) as executor:

		# Dictionary to get genomes by their future when finished
		genomes_by_future = dict()

		# Lock around database usage
		db_lock = threading.Lock()

		# Flag to inform worker threads that abort has happened
		is_aborted = False
		check_aborted = lambda: is_aborted

		# Cancels all pending futures and sets abort flag
		def abort():
			nonlocal is_aborted
			is_aborted = True
			for future in genomes_by_future.keys():
				future.cancel()

		try:

			# Submit the tasks
			for item in items:

				if len(item) == 2:
					genome, url = item
					store_opts = common_store_opts
				elif len(item) == 3:
					genome, url, store_opts = item
				else:
					raise ValueError('Download items shoud be 2- or 3-tuples')

				future = executor.submit(download_sequence, db, genome, url,
				                         db_lock=db_lock, store_opts=store_opts,
				                         aborted=check_aborted, **kwargs)
				genomes_by_future[future] = genome

		# If any errors in submitting jobs
		except BaseException as exc:
			abort()
			raise

		# Iterator over futures as they complete
		completed = iter(as_completed(genomes_by_future.keys()))

		# Don't use for loop as we want to wrap every call to
		# completed.__next__ in a try-catch
		while True:

			# Catch exceptions anywhere in loop
			try:

				# Get next completed future
				try:
					future = next(completed)
				except StopIteration:
					break

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
						abort()
						break

			# Cancel pending futures on first KeyboardInterrupt but keep
			# looping
			except KeyboardInterrupt:
				if not is_aborted:
					print('Keyboard interrupt caught, attempting to shut down '
					      'worker pool gracefully...',
					      file=sys.stderr)
					abort()

				else:
					raise

			# All other exceptions abort and then propagate
			except Exception:
				abort()
				raise

		return not is_aborted
