"""Alternative taxonomy classification methods."""

import functools
import json
import pickle
from typing import List

import numpy as np

from pydatatypes import dataclass, field
from .kmers import KmerSpec


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

	Subclasses must implement the :meth:`predict` method.

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
			self._params,
			self._fit_params
		)

	def __setstate__(self, state):
		*init_args, params, fit_params = state
		Classifier.__init__(self, *init_args)
		self._params = dict(params)
		self._fit_params = dict(fit_params)

	def predict(self, x):
		"""Predict class labels for a collection of signatures or feature vectors

		:param x: Single signature in corrdinate format, sequence of signatures,
			or a feature vector or matrix.

		:returns: Predicted class labels for x. Will be a scalar integer if
			argument was a single signature or feature vector, or a 1D array
			if argument was a sequence of signatures or a feature matrix. Array
			elements are integers from 0 to ``n_classes - 1`` which correspond
			to class assignments, or negative values which indicate no class/
			unsure.
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


@dataclass(json=True)
class ClassifierInfo:
	"""
	Data object which describes how a classifier links to MIDAS database objects.

	.. attribute:: id

		String ID intended to uniquely identify the classifier for the purposes
		of distributing updates, etc.

	.. attribute:: version

		Version string indicating the current revision of the ID. Should be
		digits separated by commas, e.g. ``'1.0'``.

	.. attribute:: parent_taxon

		Name of the :class:`midas.db.models.taxon` this classifier works within.

	.. attribute:: class_taxa

		Ordered list of taxa names that corresond to the classes of the model.

	.. attribute:: kspec

		:class:`midas.kmers.KmerSpec` used to calculate features used as input
		to the model (see :attr:`.Classifier.kspec)`.

	.. attribute:: description

		Optional string with longer description of the classifier.

	.. attribute:: metadata

		Optionable ``dict`` containing arbitrary extra metadata for the
		classifier. Should be convertible to JSON.
	"""

	id = field(str)
	version = field(str)
	parent_taxon = field(str)
	class_taxa = field(List[str])
	kspec = field(KmerSpec)
	description = field(str, optional=True)
	metadata = field(dict, optional=True)


def _check_classifier_info(classifier, info):
	"""Check for consistency between ClassifierInfo and Classifier instance.

	:param classifier: Classifier instance.
	:type classifier: .Classifier
	:param info: ClassifierInfo describing the classifier
	:type info: .ClassifierInfo

	:raises ValueError: If there is an inconsistency between the classifier and
		the info object.
	"""
	if classifier.n_classes != len(info.class_taxa):
		raise ValueError('classifier.n_classes does not match length of info.class_taxa')

	if classifier.kspec != info.kspec:
		raise ValueError('classifier.kspec does not match info.kspec')


def dump_classifier(stream, info, classifier, check=True):
	"""Save Classifier and ClassifierInfo to a binary stream.

	:param stream: Writeable stream in binary mode.
	:param info: Classifier info.
	:type info: .ClassifierInfo
	:param classifier: Classifier instance.
	:type classifier: .Classifier
	:param bool check: If True check that the classifier and info are consistent
		with one another.

	:raises ValueError: If ``check`` is True and the classifier and info objects
		are inconsistent.
	"""

	from io import TextIOWrapper

	if check:
		_check_classifier_info(classifier, info)

	# Dump info as JSON
	# New lines in strings should be encoded as "\n" so there shouldn't be any
	# newline characters until the separator.
	text = TextIOWrapper(stream)
	json.dump(info.to_json(), text)
	text.detach()

	# Separate by newline
	stream.write(b'\n')

	# Dump classifier using pickle
	pickle.dump(classifier, stream)


def load_classifier(stream, info_only=False, check=True):
	"""Read Classifier and ClassifierInfo from a binary stream.

	:param stream: Readable stream in binary mode.
	:param bool info_only: Return only the :class:`.ClassifierInfo` object, not
		the :class:`.Classifier`.
	:param bool check: If True check that the classifier and info are consistent
		with one another.

	:returns: Classifier info only if ``info_only`` is True, othewise
		``(info, classifier)`` pair.

	:raises ValueError: If ``check`` is True and the classifier and info objects
		are inconsistent.
	"""

	# Read info
	info_bytes = stream.readline()  # Reads up to newline
	info = ClassifierInfo.from_json(json.loads(info_bytes.decode()))

	if info_only:
		return info

	# Remaining data should be pickled classifier instance
	classifier = pickle.load(stream)

	if check:
		_check_classifier_info(classifier, info)

	return info, classifier


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


class GenerativeClassifier(Classifier):
	"""
	ABC for classifiers which are generative models and fit a probability
	distribution to the inputs/outputs and classify according to the likelihood/
	probability of the class labels.

	Subclasses must implement :meth:`_log_prob`.

	:param float margin: Minimum log-ratio between 1st and 2nd largest class
		probabilities for a call to be made.
	"""

	margin = Classifier.Param('margin')

	def __init__(self, n_classes, kspec, kmers, margin=0):
		Classifier.__init__(self, n_classes, kspec, kmers)
		self.margin = margin

	@Classifier._process_features
	def log_prob(self, x, normalize=False):
		"""Get the log-probability of each signature/class pair.

		The values will be the probabilities of elements of ``x`` given
		each class label if ``normalize`` is False. If ``normalize`` is True the
		rows will be normalized so that the values sum to one.

		:param x: One or more signatures or feature vectors. See argument to
			:meth:`.Classifier.predict`.
		:param bool normalize: Normalize rows of return values so that their
			exponentials sum to one.
		:returns: Matrix of log probabilities where rows correspond to samples
			and columns to classes. Will be a vector if a single signature or
			vector is given.
		:rtype: np.ndarray
		"""
		log_probs = self._log_prob(x)

		if normalize:
			return self._normalize(log_probs)

		return log_probs

	def _log_prob(self, x):
		"""Actual implementation of :meth:`log_prob`."""
		raise NotImplementedError()

	@staticmethod
	def _normalize(log_probs, return_log=True):
		"""Normalize a matix of log-probabilities so that probabilities sum to one.

		WARNING: this modifies the first argument in-place.

		:param log_probs: 2D matrix of log-probabilities in sample x class format.
		:type log_probs: np.ndarray
		:param bool return_log: Return log-probabilities (True) or regular
			probabilities (False).
		:returns: Normalized probabilities.
		:rtype: np.ndarray
		"""

		# Scale rows to have max of 0 first, to avoid all elements
		# underflowing when we calculate exponential
		log_probs -= np.max(log_probs, axis=-1).reshape(-1, 1)

		# Expect underflow errors to occur
		with np.errstate(under='ignore'):
			probs = np.exp(log_probs)

		rowsums = np.sum(probs, -1, keepdims=True)

		if return_log:
			# Subtract row sums
			log_probs -= np.log(rowsums)
			return log_probs

		else:
			# Divide by row sums
			probs /= rowsums
			return probs

	@Classifier._process_features
	def prob(self, x, normalize=False):
		"""Get the likelihood/probability of each signature/class pair.

		See documentation for :meth:`log_prob`.

		:param x: One or more signatures or feature vectors. See argument to
			:meth:`.Classifier.predict`.
		:param bool normalize: Normalize rows of return values so that they sum
			to one.
		:returns: Matrix of probabilities where rows correspond to samples
			and columns to classes. Will be a vector if a single signature or
			vector is given.
		:rtype: np.ndarray
		"""
		# Derive from log probs
		log_probs = self._log_prob(x)

		if normalize:
			return self._normalize(log_probs, return_log=False)

		else:
			# Take exponential, expect underflow to occur
			with np.errstate(under='ignore'):
				return np.exp(log_probs)

	def predict(self, x):
		log_probs = self.log_prob(x)

		# No margin, just take argmax across classes
		if self.margin == 0:
			return np.argmax(log_probs, -1)

		# More complicated...
		else:
			sort = np.argsort(log_probs, -1)
			i = np.arange(x.shape[0])

			# Difference in log p between 1st and 2nd most likely
			diffs = log_probs[i, sort[:, -1]] - log_probs[i, sort[:, -2]]

			# Replace predictions with -1 where difference less than margin
			predictions = sort[:, -1].copy()
			predictions[diffs < self.margin] = -1

			return predictions


class NaiveBayesClassifier(GenerativeClassifier):
	"""Naive Bayes classifier model.

	:param float alpha: Additive smoothing parameter to use when training.
	"""

	alpha = Classifier.Param('alpha')
	log_p_ = Classifier.FitParam('log_p')
	log_1mp_ = Classifier.FitParam('log_1mp')

	def __init__(self, n_classes, kspec, kmers, alpha=1.0, margin=0):
		GenerativeClassifier.__init__(self, n_classes, kspec, kmers, margin=margin)
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
		self.log_p_ = np.log(class_freqs)
		self.log_1mp_ = np.log(1 - class_freqs)

	def _log_prob(self, x):
		return np.dot(x, self.log_p_.T) + np.dot(~x, self.log_1mp_.T)
