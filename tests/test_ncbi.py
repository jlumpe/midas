"""Test midas.ncbi."""

import pytest

from midas import ncbi


def test_check_seq_ids():
	"""Test midas.ncbi.check_seq_ids."""

	# Valid non-empty
	ncbi.check_seq_ids(dict(genbank_acc='XX_123456.1'))
	ncbi.check_seq_ids(dict(entrez_db='nuccore', entrez_id=12345))

	# Empty
	ids = dict()

	with pytest.raises(TypeError):
		ncbi.check_seq_ids(dict(), empty_ok=False)

	ncbi.check_seq_ids(dict(), empty_ok=True)

	# Extra keys
	ids = dict(genbank_acc='XX_123456.1', foo='bar')

	with pytest.raises(KeyError):
		ncbi.check_seq_ids(ids)

	ncbi.check_seq_ids(ids, extra_ok=True)

	# Empty + extra keys
	ids = dict(foo='bar')

	with pytest.raises(TypeError):
		ncbi.check_seq_ids(ids, extra_ok=True, empty_ok=False)

	ncbi.check_seq_ids(ids, extra_ok=True, empty_ok=True)

	# Incomplete index
	with pytest.raises(KeyError):
		ncbi.check_seq_ids(dict(entrez_db='nuccore'))

	# Invalid types
	with pytest.raises(TypeError):
		ncbi.check_seq_ids(dict(genbank_acc=0))
	with pytest.raises(TypeError):
		ncbi.check_seq_ids(dict(entrez_db='nuccore', entrez_id='foo'))

	# Null values
	ids = dict(genbank_acc='XX_123456.1', refseq_acc=None)

	with pytest.raises(TypeError):
		ncbi.check_seq_ids(ids)

	ncbi.check_seq_ids(ids, null_ok=True)

	# Null values and empty
	ids = dict(genbank_acc=None)

	with pytest.raises(TypeError):
		ncbi.check_seq_ids(ids, null_ok=True, empty_ok=False)

	ncbi.check_seq_ids(ids, null_ok=True, empty_ok=True)
