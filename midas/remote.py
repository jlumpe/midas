"""Utilities for downloading remote data, e.g. from GenBank"""

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
                      **kwargs):
	"""Download sequence for genome and store it

	Used as thread worker for download_genome_sequences.
	"""

	# Download as bytes
	seq_data = urlopen_retry(url, **kwargs)

	# Place in buffer
	buf = io.BytesIO()
	shutil.copyfileobj(seq_data, buf)
	buf.seek(0)

	# Store
	with contextlib.ExitStack() as exitstack:
		if db_lock:
			exitstack.enter_context(db_lock)

		db.store_sequence(genome, buf, src_mode='b', **store_opts)


def download_sequences_parallel(db, items, overwrite=False, callback=None,
                                **kwargs):
	"""Download and store set of sequences in parallel"""

	max_workers = kwargs.pop('max_workers', None)

	with ThreadPoolExecutor(max_workers) as executor:

		# Lock around database usage
		db_lock = threading.Lock()

		# Submit the tasks
		genomes_by_future = dict()
		for genome, url in items:
			future = executor.submit(download_sequence, db, genome, url,
			                         db_lock=db_lock, **kwargs)
			genomes_by_future[future] = genome

		# Iterate over completed tasks
		for future in as_completed(genomes_by_future.keys()):
			genome = genomes_by_future[future]

			# Get exception for callback
			try:
				result = future.result()

			except Exception as exc:
				callback(genome, exc)

			else:
				callback(genome, None)
