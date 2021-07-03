"""Export results to JSON."""

import json
from functools import singledispatchmethod

from attr import attrs, attrib, asdict

from .base import AbstractResultsExporter
from midas.query.results import QueryResults, QueryResultItem
from midas.db.models import ReferenceGenomeSet, Taxon
import midas.io.json as mjson


@attrs()
class JSONResultsExporter(AbstractResultsExporter):
	"""Exports query results in JSON format.

	Attributes
	----------
	dense
		Write with no whitespace to cut down on file size. Disable to produce more human-friendly
		output. Defaults to True.
	"""
	dense: bool = attrib(default=True)

	@staticmethod
	def _todict(obj, attrs):
		return {a: getattr(obj, a) for a in attrs}

	@singledispatchmethod
	def _to_json(self, obj):
		# Base case
		return mjson.to_json(obj)

	@_to_json.register
	def _results_to_json(self, results: QueryResults):
		return asdict(results)

	@_to_json.register
	def _item_to_json(self, item: QueryResultItem):
		return asdict(item)

	@_to_json.register
	def _genomeset_to_json(self, gset: ReferenceGenomeSet):
		return self._todict(gset, ['id', 'key', 'version', 'name', 'description'])

	@_to_json.register
	def _taxon_to_json(self, taxon: Taxon):
		return self._todict(taxon, ['id', 'name', 'ncbi_id', 'distance_threshold'])

	def export(self, f, results: QueryResults):
		json.dump(results, f, default=self._to_json)
