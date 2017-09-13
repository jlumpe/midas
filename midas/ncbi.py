"""Tools for working with NCBI databases.


.. data:: NCBI_BASE_URL

	Base URL for NCBI web site.

.. data:: SEQ_IDS

	Indices/IDs which can uniquely identify an NCBI sequence. Mapping from index
	names to tuples of attribute names. Simple indices have a single attribute,
	compound indices have more than one. The current indices are:

	* ``entrez (entrez_db, entrez_id)`` -- Compound index of Entrez database
	  name (``str``) and ID within that database (``int``).
	* ``genbank_acc (genbank_acc)`` -- Genbank accession number.
	* ``refseq_acc (refseq_acc)`` -- RefSeq accession number.

.. data:: SEQ_ID_ATTRS

	Names of attributes which represent IDs for sequences in NCBI databases.
	Sequence of strings, corresponds to concatenation of values in
	:data:`.SEQ_IDS`.

.. data:: SEQ_ID_ATTR_TYPES

	Mapping from attribute names in :data:`.SEQ_ID_ATTRS` to their types.
"""

from abc import ABCMeta, abstractproperty
from collections import OrderedDict
from types import MappingProxyType


NCBI_BASE_URL = 'https://www.ncbi.nlm.nih.gov/'


SEQ_IDS = MappingProxyType(OrderedDict([
	('entrez', ('entrez_db', 'entrez_id',)),
	('genbank_acc', ('genbank_acc',)),
	('refseq_acc', ('refseq_acc',)),
]))


SEQ_ID_ATTRS = tuple(attr for attrs in SEQ_IDS.values() for attr in attrs)


SEQ_ID_ATTR_TYPES = MappingProxyType({
	'entrez_db': str,
	'entrez_id': int,
	'genbank_acc': str,
	'refseq_acc': str,
})


def check_seq_id_value(id_name, id_vals):
	"""Check the values for the attributes of an NCBI sequence ID.

	:param str id_name: Name of ID (key of :data:`.SEQ_IDS`).
	:param id_vals: Tuple of values corresponding to the attributes of the ID.

	:raises ValueError: If ``id_vals`` is not the correct length.
	:raises TypeError: If any elements of ``id_vals`` are not of the correct
		type (see :data:`SEQ_ID_ATTR_TYPES`).
	"""

	attr_names = SEQ_IDS[id_name]

	if len(attr_names) != len(id_vals):
		raise ValueError(
			'NCBI sequence ID {!r} has {} elements but got {}'
			.format(id_name, len(attr_names), len(id_vals))
		)

	for name, val in zip(attr_names, id_vals):
		type_ = SEQ_ID_ATTR_TYPES[name]
		if not isinstance(val, type_):
			raise TypeError(
				'NCBI sequence ID attribute {!r} should be of type {}, not {}'
				.format(name, type_.__name__, type(val).__name__)
			)


def get_seq_ids(mapping, single=False, *, extra_ok=False, empty_ok=True,
                null_ok=False, ignore_partial=False, check_types=True):
	"""Get NCBI sequence IDs from a mapping of attribute values.

	Checks all ID values are of the correct type and that no incomplete indices
	are present (keys are present for some but not all attributes in an item of
	:data:`.SEQ_IDS`\\ ).

	Can also be used to check if a set of attribute values are valid, in which
	case the return value can be ignored.

	:param mapping: Mapping from NCBI sequence ID attribute names to their
		values. See :data:`.SEQ_ID_ATTRS`.
	:param bool single: If True, expect attribute values for a single ID to be
		present in the mapping and return a pair of ID name and value. If False
		(default), return a mapping of a variable number ID names to values.
	:param bool extra_ok: If True, ignore keys in ``mapping`` that are not in
		:data:`.SEQ_ID_ATTRS`. Default False.
	:param bool empty_ok: If False and no IDs found, raise :exc:TypeError:.
		Default True.
	:param bool null_ok: If True allow values of None for attributes. These will
		be treated as if the key for the attribute was not present. If False
		(default) raise :exc:`TypeError` when value is None.
	:param bool ignore_partial: If True ignore cases where some but not all
		attributes of a compound ID are present (or not None). If False (default)
		will raise :exc:`TypeError`.
	:param bool check_types: If True check types of non-null attribute values
		and raise :exc:`TypeError`  if they do not match
		:data:`.SEQ_ID_ATTR_TYPES`.

	:raises KeyError: If extra keys are found and ``extra_ok`` is False.
	:raises TypeError: If ``single`` is True and more or less than one ID value
		found, if ``empty_ok`` is False and no ID values found, if ``null_ok``
		is False any any ID attribute values are None, if ``ignore_partial`` is
		False and any IDs have only a partial set of attributes present, if
		``check_types`` is True and any attribute values are not None or of the
		correct type.
	"""

	mapping = dict(mapping)

	# Mapping from ID names to values
	ids = dict()

	# Look for each NCBI ID
	for id_name, id_attrs in SEQ_IDS.items():

		attr_vals = []

		# Try to get the value for each attribute from the mapping
		for attrname in id_attrs:
			try:
				attrval = mapping.pop(attrname)

			except KeyError:
				attrval = None

			else:
				# Key there but value is None - raise exception if not allowed
				if not null_ok and attrval is None:
					raise TypeError(
						'Value of NCBI sequence ID attribute {!r} cannot be None'
						.format(attrname)
					)

			attr_vals.append(attrval)

		if None in attr_vals:
			# Either at least one key was missing or had a null value

			# All None - Nothing at all found, just skip
			if all(val is None for val in attr_vals):
				continue

			# We have a partial index, skip or raise exception
			if ignore_partial:
				continue

			# Format exception message and raise
			null_attr = None
			not_null_attr = None
			for attrname, val in zip(id_attrs, attr_vals):
				if null_attr is None and val is None:
					null_attr = attrname
				if not_null_attr is None and val is not None:
					not_null_attr = attrname

			raise TypeError(
				'Incomplete value for NCBI sequence ID {!r}: got value for '
				'{!r} but not {!r}'
				.format(id_name, not_null_attr, null_attr)
			)

		# We're good, record it
		ids[id_name] = tuple(attr_vals)

	# Check for leftover keys
	if mapping and not extra_ok:
		raise KeyError(
			'{!r} is not an NCBI sequence ID attribute'
			.format(next(iter(mapping)))
		)

	# Check empty
	if not ids and (not empty_ok or single):
		raise TypeError('No NCBI sequence ID attributes present')

	# Check types
	if check_types:
		for id_name, id_vals in ids.items():
			check_seq_id_value(id_name, id_vals)

	if single:
		# Return exactly one

		# Got more than one ID
		if len(ids) > 1:
			raise TypeError(
				'Got values for more than one NCBI sequence ID: {!r}'
				.format(list(ids))
			)

		# Return the only item
		return next(iter(ids.items()))

	else:
		# Multiple requested, return full mapping
		return ids


def parse_seq_id_args(args, kwargs, multiple=False, **more):
	"""Parse NCBI sequence IDS from positional and keyword arguments.

	Meant to be used as a convenience in other functions where the signature
	allows you to specify a single NCBI sequence ID in positional arguments
	as ``(id_name1, *id_attr_vals...)`` or possibly multiple IDs by giving their
	attribute values in keyword arguments like
	``(... attr1=val1, attr2=val2, ...)``.

	:param tuple args: Zero or one sequence IDs encoded in positional arguments.
	:param dict kwargs: Any number of sequence ID attribute values in keyword
		arguments.
	:param bool multiple: If True get values for multiple sequence IDs. If
		False (default) expect
	:param \\**more: Additional keyword arguments to forward to
		:func:`get_seq_ids`.
	:returns: If ``multiple`` is False returns a``(id_name, id_values)`` pair,
		where ``id_name`` is a key of :data:`.SEQ_IDS` and ``id_values`` is a
		tuple of values for the corresponding attributes. If ``multiple`` is
		True returns a dictionary where elements of ``items()`` are formatted
		in the same way.
	:rtype: tuple or dict

	:raises TypeError: If both or neither of ``args`` and ``kwargs`` are empty.
	"""
	if args:
		if kwargs:
			raise TypeError(
				'Must give NCBI sequence IDS in either positional or keyword '
				'arguments, not both'
			)

		# Single ID name followed by values
		id_name, *id_vals = args

		# Check that number and type of values are OK
		check_seq_id_value(id_name, id_vals)

		# Return in correct format
		if multiple:
			return {id_name: tuple(id_vals)}
		else:
			return id_name, tuple(id_vals)

	elif kwargs:
		# Attribute values passed as keyword arguments
		return get_seq_ids(kwargs, single=not multiple, **more)

	else:
		raise TypeError('Must pass either positional or keyword arguments')


class SeqRecordBase(metaclass=ABCMeta):
	"""ABC for a record describing a sequence from an NCBI database.

	Should include all attributes from :data:`.SEQ_ID_ATTRS`.
	"""

	entrez_db = abstractproperty()
	entrez_id = abstractproperty()
	genbank_acc = abstractproperty()
	refseq_acc = abstractproperty()

	def ncbi_ids(self, flat=False):
		"""Get a dictionary of all NCBI sequence ID attribute values.

		:param bool flat: If True, get a flat mapping of attribute names to
			values. False (default) gets mapping from ID names to tuples of
			attribute values (see :data:`.SEQ_IDS`).

		:rtype: dict
		"""
		if flat:
			return {name: getattr(self, name) for name in SEQ_ID_ATTRS}
		else:
			return {
				k: tuple(getattr(self, name) for name in names)
				for k, names in SEQ_IDS.items()
			}


def entrez_url(entrez_db, entrez_id):
	"""Get the URL for an entry in Entrez, more or less.

	:param str entrez_db: Name of the Entrez database - e.g. "nuccore".
	:param int entrez_id: ID of the database entry.
	:rtype: str:
	"""
	from urllib.parse import urljoin, quote_plus

	path = [quote_plus(str(p)) for p in [entrez_db, entrez_id]]
	return urljoin(NCBI_BASE_URL, '/'.join(path))


def ncbi_sequence_url(*args, **kwargs):
	"""Attempt to guess the URL for a sequence on get_seq_ids.

	Not too smart about it right now, only operates off the entrez IDs.

	:param \\*args: Name of NCBI sequence ID followed by its attribute values
		(see :data:`.SEQ_IDS`). Mututally exclusive with ``**kwargs``.
	:param \\**kwargs: Valid set of NCBI sequence ID attribute values as keyword
		arguments (see :func:`.get_seq_ids`). Mututally exclusive with ``*args``.
	:returns: Guess for URL, or None no guess could be made.
	:rtype: str
	"""
	ids = parse_seq_id_args(args, kwargs, multiple=True, null_ok=True)

	try:
		entrez_id_vals = ids['entrez']

	except KeyError:
		return None

	return entrez_url(*entrez_id_vals)


def ncbi_search_url(term):
	"""Get the URL for the NCBI search results page for a given term.

	:param str term: Search term.
	:returns; URL for search page.
	:rtype str:
	"""
	from urllib.parse import urljoin, urlencode
	return urljoin(NCBI_BASE_URL, 'gquery/?' + urlencode({'term': term}))
