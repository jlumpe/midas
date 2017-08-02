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
