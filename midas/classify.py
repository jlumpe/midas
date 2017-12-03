"""Alternative taxonomy classification methods."""

import functools

import numpy as np


def signatures_to_features(signatures, kspec, kmers):
	"""Convert a collection of signatures to a feature matrix.

	:param signatures: Sequence of signatures in coordinate format, as
		:class:`midas.kmers.SignatureArray`, list of numpy arrays, or similar.
	:param kspec: k-mer spec.
	:type ksped: midas.kmers.KmerSpec
	:param kmers: Indices of k-mers to use for the columns in the feature
		matrix.
	:type kmers: numpy.ndarray

	:returns: A 2D boolean array, where rows correspond to elementw of
		``signatures`` and columns to elements of ``kmers``. Each entry is true
		the corresponding signature contains the corresponding kmer.
	"""

	tmp = np.zeros(kspec.idx_len, dtype=bool)
	out = np.zeros((len(signatures), len(kmers)), dtype=bool)

	for i, coords in enumerate(signatures):
		tmp[:] = 0
		tmp[coords] = True
		out[i] = tmp[kmers]

	return out


class Classifier:
	"""ABC for a model which classifies genomes based on their k-mer signatures.

	Should describe the model only, not any metadata such as IDs or references
	to taxonomy nodes, etc.

	:param int n_classes: Number of classes the classifier has been trained on.
	:param kspec: K-mer spec used to calculate signatures.
	:type kspec: midas.kmers.KmerSpec
	:param kmers: Indices of k-mers which the classifier uses as features.
	:type kmers: numpy.ndarray
	"""

	def __init__(self, n_classes, kspec, kmers):
		self.n_classes = n_classes
		self.kspec = kspec

		self.kmers = np.asarray(kmers, dtype=kspec.coords_dtype)
		np.sort(self.kmers)

		self._params = dict()
		self._fit_params = dict()

	def __getstate__(self):
		return (
			self.n_classes,
			self.kspec,
			self.kmers,
			self.params,
			self.fit_params
		)

	def __setstate__(self, state):
		*init_args, params, fit_params = state
		Classifier.__init__(*init_args)
		self._params = dict(params)
		self._fit_params = dict(fit_params)

	def predict(self, x):
		"""Predict class labels for a collection of signatures or feature vectors

		:param x: Single signature in corrdinate format, sequence of signatures,
			or a feature vector or matrix.

		:returns: Predicted class labels for x. Will be a scalar integer if
			argument was a single signature or feature vector, or a 1D array
			if argument was a sequence of signatures or a feature matrix.
		"""
		raise NotImplementedError()

	def signatures_to_features(self, signatures):
		"""Convert signatures to a binary feature matrix for the model.

		:param signatures: Sequence of k-mer signatures in coordinate format.
		:returns: Matrix of boolean type where rows are feature vectors of
			elements of ``signatures``.
		:rtype: np.ndarray
		"""
		return signatures_to_features(signatures, self.kspec, self.kmers)

	@staticmethod
	def _process_features(method):
		"""Wrap a method to convert first argument to feature matrix.

		Wrapper first argument to a feature matrix, calls wrapped method,
		and drops first axis of return value if necessary.

		:param method: *Unbound* method that takes a feature matrix as its
			first argument and returns an array of values for each sample.

		:returns: Wrapper around argument.
		"""

		@functools.wraps(method)
		def wrapper(self, x, *args, **kwargs):

			# Convert argument to 2d feature matrix
			if isinstance(x, np.ndarray):
				if x.ndim == 2:
					# Feature matrix
					single_signature = False

				elif x.ndim == 1:
					if x.dtype.kind == 'b':
						# Feature vector
						single_signature = True
						x = x.reshape(1, -1)

					elif x.dtype.kind in 'ui':
						# Single signature
						single_signature = True
						x = signatures_to_features([x], self.kspec, self.kmers)

					else:
						raise ValueError('Array argument must have bool or int dtype')

				else:
					raise ValueError('Array argument must be 1 or 2 dimensional')

			else:
				# Signature sequence
				x = signatures_to_features(x, self.kspec, self.kmers)
				single_signature = False

			# Call function with feature matrix
			y = method(self, x, *args, **kwargs)

			# Drop first axis if necessary
			if single_signature:
				return y[0]
			else:
				return y

		return wrapper

	class Param:

		def __init__(self, key):
			self.key = key

		def __get__(self, obj, cls):
			if obj is None:
				return self
			else:
				return obj._params.get(self.key, None)

		def __set__(self, obj, value):
			obj._params[self.key] = value

	class FitParam:

		def __init__(self, key):
			self.key = key

		def __get__(self, obj, cls):
			if obj is None:
				return self
			else:
				return obj._fit_params.get(self.key, None)

		def __set__(self, obj, value):
			obj._fit_params[self.key] = value


def kmer_class_freqs(x, labels, k=None, sample_weights=None, smoothing=None):
	"""Find k-mer frequencies within classes.

	:param x: K-mer feature matrix (see :func:`.signatures_to_features`).
	:type x: np.ndarray
	:param labels: Nonnegative integer labels matching first axis of ``x``.
	:type labels: np.ndarray
	:param int k: Number of classes. If None will use ``1 + max(labels)``.
	:param sample_weights: Weights for rows of ``x``, to perform a weighted
		average.
	:param float smoothing: Additive smoothing to use for frequency calculations.
		Adds this amount to the count of both True and False values for each
		k-mer.
	"""

	if k is None:
		k = np.max(labels) + 1

	out = np.zeros((k, x.shape[1]))

	for i in range(k):
		in_class = labels == i

		if sample_weights is None:
			counts = np.sum(x[in_class, :], axis=0)
			total = in_class.sum()

		else:
			class_weights = sample_weights[in_class]
			counts = np.sum(x[in_class, :] * class_weights[:, None], axis=0)
			total = class_weights.sum()

		if smoothing is None:
			out[i] = counts / total
		else:
			out[i] = (counts + smoothing) / (total + 2 * smoothing)

	return out


class NaiveBayesClassifier(Classifier):
	"""Naive Bayes classifier model.

	:param float alpha: Additive smoothing parameter to use when traiing.
	"""

	alpha = Classifier.Param('alpha')
	log_p = Classifier.FitParam('log_p')
	log_1mp = Classifier.FitParam('log_1mp')

	def __init__(self, n_classes, kspec, kmers, alpha=1.0):
		Classifier.__init__(self, n_classes, kspec, kmers)
		self.alpha = alpha

	def fit(self, x, y, sample_weights=None):
		"""Fit the model.

		:param x: K-mer feature matrix.
		:type x: np.ndarray
		:param y: Array of integer label values corresonding to rows of ``x``.
		:type y: np.ndarray
		:param sample_weights: Weights of samples.
		:type sample_weights: np.ndarray
		"""

		class_freqs = kmer_class_freqs(
			x,
			y,
			self.n_classes,
			sample_weights=sample_weights,
			smoothing=self.alpha,
		)

		# Just don't support zero probabilities for now
		if np.any(class_freqs == 0):
			assert self.alpha > 0
			raise ValueError(
				'One or more class k-mer frequencies is zero, try again with '
				'alpha > 0'
			)

		# Log of probability and 1 - probability
		self.log_p = np.log(class_freqs)
		self.log_1mp = np.log(1 - class_freqs)

	@Classifier._process_features
	def log_prob(self, x, normalize=False):
		"""Get the log-likelihood/probability of each signature/class pair.

		The probabilities will be the likelihoods of elements of ``x`` given
		each class label if ``normalize`` is False. If ``normalize`` is True the
		rows will be normalized so that the probabilities sum to one, making
		them the probability of class labels given each element of ``x``.

		:param x: One or more signatures or feature vectors. See argument to
			:meth:`.Classifier.predict`.
		:param bool normalize: Normalize rows of return values so that their
			exponentials sum to one.
		:returns: Matrix of log probabilities where rows correspond to samples
			and columns to classes. Will be a vector if a single signature or
			vector is given.
		:rtype: np.ndarray
		"""

		probs = np.dot(x, self.log_p.T) + np.dot(~x, self.log_1mp.T)

		if normalize:
			# Scale rows to have max of 0 first, to avoid all elements
			# underflowing when we calculate exponential
			probs -= np.max(probs, axis=-1).reshape(-1, 1)

			# Expect underflow and divide by zero errors to occur
			with np.errstate(under='ignore', divide='ignore'):
				# Sums of exponentiated rows
				rowsums = np.sum(np.exp(probs), -1)

				# Subtract from result
				probs -= np.log(rowsums).reshape(-1, 1)

		return probs

	def prob(self, x, normalize=False):
		"""Get the likelihood/probability of each signature/class pair.

		See documentation for :meth:`log_prob`.

		:param x: One or more signatures or feature vectors.
		:param bool normalize: Normalize rows of return values so that they sum
			to one.
		:returns: Matrix of likelihoods/probabilities.
		:rtype: np.ndarray
		"""
		probs = self.log_prob(x, normalize=False)

		# Take exponential, expect underflow to occur
		with np.errstate(under='ignore'):
			probs = np.exp(probs)

		if normalize:
			probs *= np.sum(probs, -1).reshape(-1, 1)

		return probs

	def predict(self, x):
		return np.argmax(self.log_prob(x), -1)
