"""Utilities for parallel processing."""

import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
import signal


class AbortableThreadPoolExecutor(ThreadPoolExecutor):
	"""ThreadPoolExecutor with ability to abort all futures cleanly.

	Used as standard :class:`concurrent.futures.ThreadPoolExecutor`, with the
	addition of an ``abort()`` method which attempts to cleanly shut down all
	current jobs.

	The *first* time a SIGINT / ``KeyboardInterrupt`` is caught in its context,
	``abort()`` will be called and the exception will not propagate (future
	ones will, however).

	You can check whether the pool was aborted for any reason with
	:meth:`is_aborted` after the context has exited.
	"""

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self._aborted = False
		self._futures = set()

	def is_aborted(self):
		"""Check if abort() has been called."""
		return self._aborted

	def abort(self):
		"""Cancel all pending futures.

		Should only be called from main thread.
		"""
		if not self._aborted:
			self._aborted = True
			for future in self._futures:
				future.cancel()

	def _handle_sigint(self, signum, frame):
		"""Handler for SIGINT when in context"""

		# Abort and continue if not already aborted
		if not self._aborted:
			print('Keyboard interrupt caught, attempting to shut down '
			      'worker pool gracefully...',
			      file=sys.stderr)
			self.abort()

		else:
			self._original_sigint_handler(signum, frame)

	def __enter__(self):
		retval = super().__enter__()

		# Set SIGINT handler (keep track of old one)
		self._original_sigint_handler = signal.getsignal(signal.SIGINT)
		signal.signal(signal.SIGINT, self._handle_sigint)

		return retval

	def __exit__(self, exc_type, exc_value, traceback):
		# Restore default SIGINT handler
		signal.signal(signal.SIGINT, self._original_sigint_handler)

		# Cancel pending futures if exiting with an exception
		if exc_type is not None:
			self.abort()

		super().__exit__(exc_type, exc_value, traceback)

	def submit(self, *args, **kwargs):
		if self._aborted:
			return None

		future = super().submit(*args, **kwargs)
		self._futures.add(future)

		return future

	def as_completed(self, skip_cancelled=False, **kwargs):
		"""Iterate over completed futures.

		Convenience method as all futures are already stored anyways.

		Parameters
		----------
		skip_cancelled : bool
			Don't yield canceled futures.
		kwargs :
		    Passed to :func:`concurrent.futures.as_completed`.

		Returns
		-------
			Iterator over completed :class:`concurrent.futures.Future`.
		"""
		for future in as_completed(self._futures, **kwargs):
			if not (skip_cancelled and future.cancelled()):
				yield future
