"""Test midas.classify module."""


import pytest
import numpy as np

from midas import classify
from midas.classify import Classifier, ClassifierInfo
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


class TestClassifierSerialization:
	"""Test serialization of classifiers."""

	class DummyClassifier(Classifier):
		"""Simple classifier for testing."""
		a = Classifier.Param('a')
		b = Classifier.Param('b')
		c = Classifier.FitParam('c')

		def __init__(self, n_classes, kspec, kmers, a, b, c):
			Classifier.__init__(self, n_classes, kspec, kmers)
			self.a = a
			self.b = b
			self.c = c

	@pytest.fixture(scope='class')
	def classifier(self):
		"""Test classifier instance."""
		kspec = KmerSpec(11, 'ATGAC')
		kmers = np.arange(10000) * 100
		return self.DummyClassifier(10, kspec, kmers, a=1, b=2, c=3)

	@pytest.fixture(scope='class')
	def info(self, classifier):
		"""Test ClassifierInfo instance."""
		return ClassifierInfo(
			id='foo',
			version='1.0',
			parent_taxon='E coli',
			class_taxa=['E coli {}'.format(i + 1) for i in range(classifier.n_classes)],
			kspec=classifier.kspec,
		)

	@pytest.fixture()
	def stream(self, info, classifier):
		"""Readable binary data stream with classifier written to it."""
		from io import BytesIO

		stream = BytesIO()
		classify.dump_classifier(stream, info, classifier)
		stream.seek(0)

		return stream

	def test_load_both(self, stream, info, classifier):
		info2, classifier2 = classify.load_classifier(stream)

		# Check info
		assert info2 == info

		# Check classifier
		assert classifier.n_classes == classifier2.n_classes
		assert classifier.kspec == classifier2.kspec
		assert np.array_equal(classifier.kmers, classifier2.kmers)

		assert classifier.a == classifier2.a
		assert classifier.b == classifier2.b
		assert classifier.c == classifier2.c

	def test_info_mismatch(self, classifier, info):

		from io import BytesIO
		from attr import asdict

		info_dict = asdict(info)

		bad_info = []

		# Number of classes doesn't match
		bad_info.append(ClassifierInfo(**{
			**info_dict,
			'class_taxa': info.class_taxa[:-1]
		}))
		bad_info.append(ClassifierInfo(**{
			**info_dict,
			'kspec': KmerSpec(info.kspec.k + 1, info.kspec.prefix)
		}))

		for info2 in bad_info:

			stream = BytesIO()

			# Check writing fails
			with pytest.raises(ValueError):
				classify.dump_classifier(stream, info2, classifier)

			# Write without checking, should succeed
			classify.dump_classifier(stream, info2, classifier, check=False)

			# Read should fail
			stream.seek(0)
			with pytest.raises(ValueError):
				classify.load_classifier(stream)

			# Should succeed if checking disabled
			stream.seek(0)
			classify.load_classifier(stream, check=False)


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
		assert np.all(np.isfinite(model.log_p_))
		assert np.all(np.isfinite(model.log_1mp_))
		assert np.all(model.log_p_ <= 0)
		assert np.all(model.log_1mp_ <= 0)

		# Probabilities sum to one
		assert np.allclose(np.exp(model.log_p_) + np.exp(model.log_1mp_), 1)

	@pytest.fixture(scope='class')
	def x_signatures(self):
		"""Test samples, in signature (coordinate) format.

		Don't care about modeling actual classes, sample from same distribution
		as class prototypes.
		"""

		random = np.random.RandomState(0)

		signatures = [
			np.flatnonzero(random.rand(self.KSPEC.idx_len) < .01)
			for _ in range(100)
		]

		return signatures

	@pytest.fixture(scope='class')
	def margin_model(self, model, x_signatures):

		# Hacky way to copy since I haven't written a method for it yet
		model2 = object.__new__(classify.NaiveBayesClassifier)
		model2.__setstate__(model.__getstate__())

		log_probs = model2.log_prob(x_signatures)

		probs_argsort = np.argsort(log_probs, -1)

		diffs = log_probs[:, probs_argsort[:, -1]] - log_probs[:, probs_argsort[:, -2]]
		assert np.all(diffs >= 0)

		# Put 10% of test signatures within margin
		model2.margin = np.percentile(diffs, 10)

		return model2


	def test_signatures_to_features(self, model, x_signatures):
		"""Test the signatures_to_features method."""
		rval = model.signatures_to_features(x_signatures)
		expected = classify.signatures_to_features(x_signatures, model.kspec, model.kmers)
		assert np.array_equal(rval, expected)

	def check_x_arg(self, func, kmers, x_signatures, *args, **kwargs):
		"""Check method returns consistent results for various forms of x.

		Parameters
		----------
		func : callable
			Function or method that takes signatures or feature vectors as its argument and returns
			an array of corresponding shape.
		kmers : np.ndarray
			K-mers used to construct feature matrix.
		x_signatures
			Sequence of k-mer signatures.
		\\*args
			Additional positional arguments to ``func``.
		\\**kwargs
			Additional keyword arguments to ``func``.

		Returns
		-------
		np.ndarray
			Return value of ``func`` on ``x_signatures``.
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

	def test_predict(self, model, margin_model, x_signatures):
		"""Check the predict() method."""

		result = self.check_x_arg(model.predict, model.kmers, x_signatures)
		margin_result = self.check_x_arg(margin_model.predict, margin_model.kmers, x_signatures)

		# Check result is integer vector
		for r in [result, margin_result]:
			assert r.shape == (len(x_signatures),)
			assert r.dtype.kind in 'iu'

		# Check margin-less result has margins within correct range
		assert np.all((result >= 0) & (result < model.n_classes))

		# Check margin results match values where non-negative
		assert np.any(margin_result < 0)
		assert np.all((margin_result == result) | (margin_result < 0))

		# TODO - check values more thoroughly

	@pytest.mark.parametrize('log', [False, True])
	def test_prob(self, model, x_signatures, log):
		"""Check the prob() and log_prob() methods."""

		method = model.log_prob if log else model.prob
		result = self.check_x_arg(method, model.kmers, x_signatures, normalize=False)
		norm_result = self.check_x_arg(method, model.kmers, x_signatures, normalize=True)

		# Check result shape and dtype
		assert result.shape == (len(x_signatures), self.N_CLASSES)
		assert result.dtype.kind == 'f'
		assert norm_result.shape == (len(x_signatures), self.N_CLASSES)
		assert norm_result.dtype.kind == 'f'

		# Check value ranges
		if log:
			assert np.all(result <= 0)
			assert np.all(norm_result <= 0)
		else:
			assert np.all((result >= 0) & (result <= 1))
			assert np.all((norm_result >= 0) & (norm_result <= 1))

		# Check normalized probabilities sum to one
		if log:
			assert np.allclose(np.exp(norm_result).sum(axis=1), 1)
		else:
			assert np.allclose(norm_result.sum(axis=1), 1)

		# Check normalized matches unnormalized (off by constant per row)
		diffs = result - norm_result if log else result / norm_result  # No zeros so OK
		assert np.allclose(diffs, diffs[:, [0]])

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
