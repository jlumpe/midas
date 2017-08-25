"""Tools for working with NCBI databases.


.. data:: NCBI_BASE_URL

	Base URL for NCBI web site.

.. data:: SEQ_ID_ATTRS

	Attributes which represent IDs for sequences in NCBI databases. List of
	strings.

.. data:: SEQ_ID_ATTR_TYPES

	Mapping from attribute names in :data:`.SEQ_ID_ATTRS` to their types.

.. data:: SEQ_ID_INDICES

	A partition of the attribute names in :data:`.SEQ_ID_ATTRS` into sets which
	form unique indices - a combination of values for the attributes in an index
	matches	at most one NCBI sequence. List of sets.
"""

from abc import ABCMeta, abstractproperty
from types import MappingProxyType


NCBI_BASE_URL = 'https://www.ncbi.nlm.nih.gov/'


SEQ_ID_ATTRS = [
	'entrez_db',
	'entrez_id',
	'genbank_acc',
	'refseq_acc',
]


SEQ_ID_ATTR_TYPES = MappingProxyType(dict(
	entrez_db=str,
	entrez_id=int,
	genbank_acc=str,
	refseq_acc=str,
))


SEQ_ID_INDICES = [
	frozenset({'entrez_db', 'entrez_id'}),
	frozenset({'genbank_acc'}),
	frozenset({'refseq_acc'}),
]


def check_seq_ids(ids, extra_ok=False, empty_ok=True, null_ok=False):
	"""Raise an exception if a dictionary of NCBI sequence IDs is not valid.

	Checks all ID values are of the correct type and that no incomplete indices
	are present (keys are present for some but not all attributes in an item of
	:data:`.SEQ_ID_INDICES`\\ ).

	:param dict ids: Dictionary of ID values, with keys matching the names in
		:data:`.SEQ_ID_ATTRS`.
	:param bool extra_ok: If additional keys are allowed in ``ids``.
	:param bool empty_ok: If False and ``ids`` is empty or all values are
		``None`` a :exc:`KeyError` will be raised.
	:param bool null_ok: If values are allowed to be ``None``.

	:raises KeyError: If extra keys are found and ``extra_ok`` is False or if an
		incomplete set of keys for a unique index is specified.
	:raises TypeError: If any ID values are of an incorect type (including
		``None`` if ``null_ok`` is False) or  if ``ids`` is empty (or all Nones)
		and ``empty_ok`` is False,
	"""

	# Check empty / values all None
	if not empty_ok and\
			all(ids.get(key, None) is None for key in SEQ_ID_ATTRS):
		raise TypeError('Must give at least one ID value')

	# Check unique indices
	for index_keys in SEQ_ID_INDICES:
		in_ids = set(filter(lambda key: ids.get(key, None) is not None, index_keys))

		if 0 < len(in_ids) < len(index_keys):
			raise KeyError(
				'Partial set of keys given for unique index {{{}}}'
				.format(', '.join(map(repr, index_keys)))
			)

	# Check extra keys
	extra_keys = set(ids).difference(SEQ_ID_ATTRS)
	if extra_keys and not extra_ok:
		raise KeyError('Invalid ID key: {!r}'.format(extra_keys.pop()))

	# Check types
	for key in SEQ_ID_ATTRS:

		try:
			value = ids[key]
		except KeyError:
			continue

		type_ = SEQ_ID_ATTR_TYPES[key]

		if value is None:
			if not null_ok:
				raise TypeError('Null value not allowed for key {!r}'.format(key))

		elif not isinstance(value, type_):
			raise TypeError(
				'Key {!r} should be {!r}, got {!r}'
				.format(key, type_, type(value))
			)


class SeqRecordBase(metaclass=ABCMeta):
	"""ABC for a record describing a sequence from an NCBI database.

	Should include all attributes from :data:`.SEQ_ID_ATTRS`.
	"""

	entrez_db = abstractproperty()
	entrez_id = abstractproperty()
	genbank_acc = abstractproperty()
	refseq_acc = abstractproperty()

	def ncbi_ids(self):
		"""Get a dictionary of all NCBI sequence ID attribute values.

		:rtype: dict
		"""
		return {key: getattr(self, key) for key in SEQ_ID_ATTRS}


def ncbi_sequence_url(**ncbi_ids):
	"""Attempt to guess the URL for a sequence on NCBI.

	Not too smart about it right now, only operates off the entrez IDs.

	:param \\**ncbi_ids: Valid set of NCBI sequence IDs as keyword arguments.
		See :func:`check_seq_ids`.
	:returns: Guess for URL, or None no guess could be made.
	:rtype: str
	"""
	from urllib.parse import urljoin, quote_plus

	ncbi_ids = {k: v for k, v in ncbi_ids.items() if v is not None}
	check_seq_ids(ncbi_ids)

	try:
		entrez_db = ncbi_ids['entrez_db']
		entrez_id = ncbi_ids['entrez_id']
	except KeyError:
		return None

	path = [quote_plus(str(p)) for p in [entrez_db, entrez_id]]
	return urljoin(NCBI_BASE_URL, '/'.join(path))


def ncbi_search_url(term):
	"""Get the URL for the NCBI search results page for a given term.

	:param str term: Search term.
	:returns; URL for search page.
	:rtype str:
	"""
	from urllib.parse import urljoin, urlencode
	return urljoin(NCBI_BASE_URL, 'gquery/?' + urlencode({'term': term}))
