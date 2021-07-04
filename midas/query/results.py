"""Data classes to represent query results."""

from typing import Optional, List, Union
from datetime import datetime

from attr import attrs, attrib

from midas.io.seq import SequenceFile
from midas.db.models import ReferenceGenomeSet, Taxon
from midas.signatures import SignaturesMeta


@attrs()
class QueryInput:
	"""Information on a query genome.

	Attributes
	----------
	label
		Some unique label for the input, probably the file name.
	file
		Source file (optional).
	"""
	label: str = attrib()
	file: Optional[SequenceFile] = attrib(default=None, repr=False)

	@classmethod
	def convert(cls, x: Union['QueryInput', SequenceFile, str]) -> 'QueryInput':
		"""Convenience function to convert flexible argument types into QueryInput.

		Accepts single string label, ``SequenceFile`` (uses file name for label), or existing
		``QueryInput`` instance (returned unchanged).
		"""
		if isinstance(x, QueryInput):
			return x
		if isinstance(x, str):
			return QueryInput(x)
		if isinstance(x, SequenceFile):
			return QueryInput(x.path.name, x)
		raise TypeError(f'Cannot convert {type(x)} instance to QueryInput')


@attrs()
class QueryResultItem:
	"""Result for a single query sequence.

	Attributes
	----------
	input
		Information on input genome.
	success
		If the query ran successfully with no fatal errors. It is still possible no match was found.
	predicted_taxon
		Predicted taxon for query genome, if any.
	report_taxon
		Predicted taxon or its first ancestor with ``report=True``.
	warnings
		List of non-fatal warning messages to report.
	error
		Message describing a fatal error which occurred, if any.
	"""
	input: QueryInput = attrib()
	success: bool = attrib()
	predicted_taxon: Optional[Taxon] = attrib(default=None)
	report_taxon: Optional[Taxon] = attrib(default=None)
	warnings: List[str] = attrib(factory=list, repr=False)
	error: Optional[str] = attrib(default=None, repr=False)


@attrs(repr=False)
class QueryResults:
	"""Results for a set of queries, as well as information on database and parameters used.

	Attributes
	----------
	items
		Results for each query sequence.
	genomeset
		Genome set used.
	signaturesmeta
		Metadata for signatures set used.
	timestamp
		Time query was completed.
	"""
	items: List[QueryResultItem] = attrib()
	genomeset: Optional[ReferenceGenomeSet] = attrib(default=None)
	signaturesmeta: Optional[SignaturesMeta] = attrib(default=None)
	timestamp: datetime = attrib(factory=datetime.now)
