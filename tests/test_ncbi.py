"""Test midas.ncbi."""

import pytest

from midas import ncbi

SEQ_IDS_GB = dict(genbank_acc='XX_123456.1')
SEQ_IDS_ENTREZ = dict(entrez_db='nuccore', entrez_id=12345)


def test_get_seq_ids_multiple():
	"""Test midas.ncbi.get_seq_ids with multiple=True."""

	from midas.ncbi import get_seq_ids

	# Valid non-empty
	assert get_seq_ids(SEQ_IDS_GB) == \
		{'genbank_acc': (SEQ_IDS_GB['genbank_acc'],)}
	assert get_seq_ids(SEQ_IDS_ENTREZ) == \
		{'entrez': tuple(SEQ_IDS_ENTREZ[attr] for attr in ncbi.SEQ_IDS['entrez'])}

	# Empty
	with pytest.raises(TypeError):
		get_seq_ids({}, empty_ok=False)

	assert get_seq_ids({}, empty_ok=True) == {}

	# Extra keys
	ids = dict(**SEQ_IDS_GB, foo='bar')

	with pytest.raises(KeyError):
		get_seq_ids(ids)

	assert get_seq_ids(ids, extra_ok=True) == get_seq_ids(SEQ_IDS_GB)

	# Empty + extra keys
	ids = dict(foo='bar')

	with pytest.raises(TypeError):
		get_seq_ids(ids, extra_ok=True, empty_ok=False)

	assert get_seq_ids(ids, extra_ok=True, empty_ok=True) == {}

	# Partial index, missing keys
	ids = dict(**SEQ_IDS_GB, entrez_db='nuccore')

	with pytest.raises(TypeError):
		get_seq_ids(ids)

	assert get_seq_ids(ids, ignore_partial=True) == get_seq_ids(SEQ_IDS_GB)

	# Invalid types - single attribute
	ids = dict(genbank_acc=0)
	with pytest.raises(TypeError):
		get_seq_ids(ids)

	get_seq_ids(ids, check_types=False)

	ids = dict(entrez_db='nuccore', entrez_id='foo')

	# Invalid types - multiple attribute
	with pytest.raises(TypeError):
		get_seq_ids(ids)

	get_seq_ids(ids, check_types=False)

	# Null values
	ids = dict(**SEQ_IDS_GB, refseq_acc=None)

	with pytest.raises(TypeError):
		get_seq_ids(ids)

	assert get_seq_ids(ids, null_ok=True) == get_seq_ids(SEQ_IDS_GB)

	# Null values and empty
	ids = dict(genbank_acc=None)

	with pytest.raises(TypeError):
		get_seq_ids(ids, null_ok=True, empty_ok=False)

	assert get_seq_ids(ids, null_ok=True, empty_ok=True) == {}

	# Partial index, some values None
	ids = dict(**SEQ_IDS_GB, entrez_db='nuccore', entrez_id=None)

	with pytest.raises(TypeError):
		get_seq_ids(ids, null_ok=True)

	assert get_seq_ids(ids, null_ok=True, ignore_partial=True) == \
		get_seq_ids(SEQ_IDS_GB)


def test_get_seq_ids_single():
	"""Test midas.ncbi.get_seq_ids with multiple=False."""

	from functools import partial

	get_ids = partial(ncbi.get_seq_ids, single=True)

	# Valid non-empty
	assert get_ids(SEQ_IDS_GB) == ('genbank_acc', (SEQ_IDS_GB['genbank_acc'],))
	assert get_ids(SEQ_IDS_ENTREZ) == (
		'entrez',
		tuple(SEQ_IDS_ENTREZ[attr] for attr in ncbi.SEQ_IDS['entrez'])
	)

	# Empty
	with pytest.raises(TypeError):
		get_ids({}, empty_ok=False)

	with pytest.raises(TypeError):
		get_ids({}, empty_ok=True)

	# Extra keys
	ids = dict(**SEQ_IDS_GB, foo='bar')

	with pytest.raises(KeyError):
		get_ids(ids)

	assert get_ids(ids, extra_ok=True) == get_ids(SEQ_IDS_GB)

	# Empty + extra keys
	ids = dict(foo='bar')

	with pytest.raises(TypeError):
		get_ids(ids, extra_ok=True, empty_ok=False)

	with pytest.raises(TypeError):
		get_ids(ids, extra_ok=True, empty_ok=True)

	# Partial index, missing keys
	ids = dict(**SEQ_IDS_GB, entrez_db='nuccore')

	with pytest.raises(TypeError):
		get_ids(ids)

	assert get_ids(ids, ignore_partial=True) == get_ids(SEQ_IDS_GB)

	# Invalid types
	ids = dict(genbank_acc=0)
	with pytest.raises(TypeError):
		get_ids(ids)

	get_ids(ids, check_types=False)

	ids = dict(entrez_db='nuccore', entrez_id='foo')

	# Null values
	ids = dict(**SEQ_IDS_GB, refseq_acc=None)

	with pytest.raises(TypeError):
		get_ids(ids)

	assert get_ids(ids, null_ok=True) == get_ids(SEQ_IDS_GB)

	# Null values and empty
	ids = dict(genbank_acc=None)

	with pytest.raises(TypeError):
		get_ids(ids, null_ok=True, empty_ok=False)

	with pytest.raises(TypeError):
		get_ids(ids, null_ok=True, empty_ok=False)

	# Partial index, some values None
	ids = dict(**SEQ_IDS_GB, entrez_db='nuccore', entrez_id=None)

	with pytest.raises(TypeError):
		get_ids(ids, null_ok=True)

	assert get_ids(ids, null_ok=True, ignore_partial=True) == get_ids(SEQ_IDS_GB)
