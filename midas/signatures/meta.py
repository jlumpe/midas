from typing import Optional, Mapping, Any

from attr import Factory

from midas.util.attr import attrs, attrib


@attrs()
class SignaturesMeta:
	"""Metadata describing a set of k-mer signatures.

	All attributes are optional.

	Attributes
	----------
	id
		Any kind of string ID that can be used to uniquely identify the signature set.
	version
		Version string (ideally PEP 440-compliant).
	name
		Short human-readable name.
	id_attr
		Name of ``Genome`` attribute the IDs correspond to (see :data:`midas.db.models.GENOME_ID_ATTRS`).
		Optional, but signature set cannot be used as a reference for queries without it.
	description
		Human-readable description.
	extra
		Extra arbitrary metadata. Should be a ``dict`` or other mapping which can be converted to JSON.
	"""

	id : Optional[str] = attrib(optional=True, kw_only=True)
	name : Optional[str] = attrib(optional=True, kw_only=True)
	version : Optional[str] = attrib(optional=True, kw_only=True)
	id_attr : Optional[str] = attrib(optional=True, kw_only=True)
	description : Optional[str] = attrib(optional=True, kw_only=True, repr=False)
	extra : Mapping[str, Any] = attrib(default=Factory(dict), kw_only=True, repr=False)
