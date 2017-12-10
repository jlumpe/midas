"""Test midas.classify module."""


import pytest
import numpy as np

from midas import classify
from midas.kmers import KmerSpec


def test_signatures_to_features():
	"""Test signatures_to_features function."""

	random = np.random.RandomState(0)

	# Random set of k-mers to use
	kspec = KmerSpec(8, 'AAA')
	n_kmers = int(kspec.idx_len * .01)
	kmers = random.choice(kspec.idx_len, n_kmers, replace=False)

	# Make random signatures
	signatures = [
		np.flatnonzero(random.rand(kspec.idx_len) < .01)
		for _ in range(1000)
	]

	# Convert signatures to features
	features = classify.signatures_to_features(signatures, kspec, kmers)

	# Check correct size and type
	assert features.shape == (len(signatures), n_kmers)
	assert features.dtype.kind == 'b'

	# Check individually
	for i in range(len(signatures)):
		assert np.array_equal(features[i], np.in1d(kmers, signatures[i]))


class TestNaiveBayesClassifier:
	"""Test NaiveBayesClassifer class."""

	N_CLASSES = 4
	N_KMERS = 10000
	N_TRAIN_SAMPLES = 1000
	KSPEC = KmerSpec(11, 'ATGAC')

	@pytest.fixture(scope='class')
	def kmers(self):
		"""K-mers to use."""
		random = np.random.RandomState(0)

		kmers = random.choice(self.KSPEC.idx_len, self.N_KMERS, replace=False)
		kmers.sort()

		return kmers

	@pytest.fixture(scope='class')
	def train_labels(self):
		"""Class labels for training samples."""

		random = np.random.RandomState(0)

		# Get some slightly uneven class frequencies
		class_freqs = random.dirichlet([3] * self.N_CLASSES)

		# Now just sample from a categorical variable for each training sample
		labels = random.choice(self.N_CLASSES, self.N_TRAIN_SAMPLES, p=class_freqs)

		# Make sure none of them are reeeaaallly small...
		assert np.all(np.bincount(labels) >= .05 * self.N_CLASSES)

		return labels

	@pytest.fixture(scope='class')
	def train_features(self, train_labels):
		"""Feature matrix for training samples."""

		random = np.random.RandomState(0)

		# "Prototype" signatures for each class
		prototypes = random.rand(self.N_CLASSES, self.N_KMERS) < .01

		# Flip a small fraction of bits in prototypes to get training samples
		flip = random.rand(self.N_TRAIN_SAMPLES, self.N_KMERS) < .001
		return prototypes[train_labels] ^ flip

	@pytest.fixture(scope='class', params=[.1, 1, 10])
	def model(self, request, train_labels, train_features, kmers):

		alpha = request.param

		model = classify.NaiveBayesClassifier(
			n_classes=self.N_CLASSES,
			kspec=KmerSpec(11, 'ATGAC'),
			kmers=kmers,
			alpha=alpha
		)

		# Fit, ignoring numpy divide-by-zero errors which may arise when alpha=0
		with np.errstate(divide='ignore'):
			model.fit(train_features, train_labels)

		return model

	def test_model_fit(self, model):
		"""Check fit parameters of model."""

		# Probabilities have correct range and are not zero
		# (2nd requirement is a current constraint on the model).
		assert np.all(np.isfinite(model.log_p))
		assert np.all(np.isfinite(model.log_1mp))
		assert np.all(model.log_p <= 0)
		assert np.all(model.log_1mp <= 0)

		# Probabilities sum to one
		assert np.allclose(np.exp(model.log_p) + np.exp(model.log_1mp), 1)

	@pytest.fixture(scope='class')
	def x_signatures(self):
		"""Test samples, im signature (coordinate) format.

		Don't care about modeling actual classes, sample from same distribution
		as class prototypes.
		"""

		random = np.random.RandomState(0)

		signatures = [
			np.flatnonzero(random.rand(self.KSPEC.idx_len) < .01)
			for _ in range(100)
		]

		return signatures

	def test_signatures_to_features(self, model, x_signatures):
		"""Test the signatures_to_features method."""
		rval = model.signatures_to_features(x_signatures)
		expected = classify.signatures_to_features(x_signatures, model.kspec, model.kmers)
		assert np.array_equal(rval, expected)

	def check_x_arg(self, func, kmers, x_signatures, *args, **kwargs):
		"""Check method returns consistent results for various forms of x.

		:param func: Function or method that takes signatures or feature vectors
			as its argument and returns an array of corresponding shape.
		:param kmers: K-mers used to construct feature matrix.
		:type kmers: np.ndarray
		:param x_signatures: Sequence of k-mer signatures.
		:param \\*args: Additional positional arguments to ``func``.
		:param \\**kwargs: Additional keyword arguments to ``func``.

		:returns: Return value of ``func`` on ``x_signatures``.
		:rtype: np.ndarray
		"""

		# Call with signatures
		result = func(x_signatures, *args, **kwargs)
		assert result.shape[0] == len(x_signatures)

		# Call with feature matrix, check same output
		x_features = classify.signatures_to_features(x_signatures, self.KSPEC, kmers)
		result2 = func(x_features, *args, **kwargs)
		assert np.array_equal(result2, result)

		# Check individual signatures
		for i, sig in enumerate(x_signatures):
			row_result = func(sig, *args, **kwargs)
			assert np.allclose(row_result, result[i])

		# Check individual feature vectors
		for i, vec in enumerate(x_features):
			row_result = func(vec, *args, **kwargs)
			assert np.allclose(row_result, result[i])

		return result

	def test_predict(self, model, x_signatures):
		"""Check the predict() method."""

		result = self.check_x_arg(model.predict, model.kmers, x_signatures)

		# Check result is integer vector
		assert result.ndim == 1
		assert result.dtype.kind in 'iu'

	@pytest.mark.parametrize('normalize', [False, True])
	def test_log_prob(self, model, x_signatures, normalize):
		"""Check the log_prob() method."""

		result = self.check_x_arg(model.log_prob, model.kmers, x_signatures, normalize=normalize)

		# Check result shape and dtype
		assert result.shape == (len(x_signatures), self.N_CLASSES)
		assert result.dtype.kind == 'f'

		# Check value ranges
		assert np.all(result <= 0)

		# Check probabilities sum to one
		if normalize:
			assert np.allclose(np.exp(result).sum(axis=1), 1)

	def test_compare_sklearn(self, model, train_labels, train_features, x_signatures):
		"""Compare model performance to equivalent Scikit-learn implementation."""

		# Try importing sklearn, skip test if not available
		pytest.importorskip('sklearn')
		from sklearn.naive_bayes import BernoulliNB

		# Create and train sklearn model
		skl_model = BernoulliNB(alpha=model.alpha, fit_prior=False)
		skl_model.fit(train_features, train_labels)

		# Get feature matrix
		x = model.signatures_to_features(x_signatures)

		# Compare calculated probabilities
		probs1 = model.log_prob(x, normalize=True)
		probs2 = skl_model.predict_log_proba(x)
		assert np.allclose(probs1, probs2)

		# Compare predictions
		predict1 = model.predict(x)
		predict2 = skl_model.predict(x)
		assert np.array_equal(predict1, predict2)