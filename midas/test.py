"""Helper functions for tests."""


def make_signatures(k, nsets, dtype):
	"""Make artificial k-mer signatures.

	:rtype: SignatureArray
	"""

	import numpy as np
	from midas.kmers import vec_to_coords, SignatureArray

	random = np.random.RandomState(seed=0)

	idx_len = 4 ** k

	signatures_list = []

	# Add empty and full sets as edge cases
	signatures_list.append(np.arange(0))
	signatures_list.append(np.arange(idx_len))

	# Use a core set of k-mers so that we get some overlap
	core_prob = max(0.01, 20 / idx_len)
	core_vec = random.rand(idx_len) < core_prob

	for i in range(nsets - 2):

		keep_core = random.rand(idx_len) < (random.rand() ** 2)

		new_vec = random.rand(idx_len) < random.uniform(core_prob * .1, core_prob)

		vec = (keep_core & core_vec) | new_vec
		signatures_list.append(vec_to_coords(vec))

	return SignatureArray.from_signatures(signatures_list, dtype=dtype)


def fill_bytearray(pattern, length):
	"""Create a bytearray with a repeating pattern.

	:param bytes pattern: Pattern to repeat in array.
	:param int length: Length of array to create.
	:returns: Filled array
	:rtype: bytearray
	"""

	array = bytearray(length)

	n = len(pattern)

	for i in range(0, length, n):
		n2 = min(n, length - i)
		array[i:i + n2] = pattern[:n2]

	return array


def make_kmer_seq(seqlen, k, prefix_len, kmer_interval, n_interval=None, seed=0):
	"""Create a DNA sequence with a known k-mer signature.

	:param int seqlen: Length of sequence.
	:param int k: Length of k-mers to find (not including prefix).
	:param int prefix_len: Prefix preceding k-mers to match.
	:param int kmer_interval: Number of nucleotides between each k-mer added.
	:param int n_interval: Every this many k-mers, add an N to the k-mer
		sequence to create a k-mer that should not be matched.
	:param int seed: Seed for PRNG.

	:returns: Tuple of (seq, kmer_vector).
	"""

	import itertools
	import numpy as np
	from midas import kmers
	from midas.cython.seqs import reverse_complement

	if kmer_interval < k + prefix_len:
		raise ValueError()

	# Seed RNG
	random = np.random.RandomState(seed)

	# Sequence of all ATN's
	seq_array = fill_bytearray(b'ATN', seqlen)

	# Choose prefix with nucleotides not found in sequence "background"
	prefix = bytes(fill_bytearray(b'CGG', prefix_len))
	kspec = kmers.KmerSpec(k, prefix)

	assert kmers.find_kmers(kspec, bytes(seq_array)).sum() == 0

	# Keep track of which kmers have been added
	vec = np.zeros(kspec.idx_len, dtype=bool)

	# Plant matches
	for i in itertools.count():

		p = i * kmer_interval

		if p + kspec.total_len >= seqlen:
			break

		# Use a k-mer without nucleotides in prefix so that we won't create
		# another accidental match
		kmer_array = bytearray(kspec.k)
		kmer_num = random.randint(2 ** k)
		for j in range(kspec.k):
			kmer_array[j] = b'AT'[(kmer_num >> j) & 1]

		# Every so often add an N just to throw things off
		invalid_kmer = n_interval is not None and i % n_interval == 0
		if invalid_kmer:
			kmer_array[0] = ord(b'N')

		kmer = bytes(kmer_array)

		if not invalid_kmer:
			vec[kmers.kmer_to_index(kmer)] = True

		match = kspec.prefix + kmer

		# Reverse every other match
		if i % 2 == 1:
			match = reverse_complement(match)

		seq_array[p:p + kspec.total_len] = match

	return bytes(seq_array), kspec, vec
