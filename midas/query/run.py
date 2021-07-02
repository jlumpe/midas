"""Query a MIDAS database."""

from typing import Collection, Sequence, Optional, Union

import numpy as np

from midas.db.midasdb import MIDASDatabase
from midas.kmers import KmerSignature
from midas.io.seq import SequenceFile
from midas.metric import jaccard_sparse_array
from .classify import find_matches, consensus_taxon, reportable_taxon
from .results import QueryInput, QueryResultItem, QueryResults


def _taxon_repr(taxon):
	"""Get a short string representation of a Taxon for logging and warning/error messages."""
	return f'{taxon.id}:{taxon.name}'


def runquery(db: MIDASDatabase,
             queries: Sequence[KmerSignature],
             inputs: Optional[Sequence[Union[QueryInput, SequenceFile, str]]],
             ) -> QueryResults:
	"""Predict the taxonomy of one or more query genomes using a given MIDAS reference database.

	Parameters
	----------
	db
		Database to query.
	queries
		Sequence of k-mer signatures.
	inputs
		Description for each input, converted to :class:`midas.query.result.QueryInput` in results
		object. Only used for reporting, does not any other aspect of results. Items can be
		``QueryInput``, ``SequenceFile`` or ``str``.
	"""
	queries = list(queries)

	if len(queries) == 0:
		raise ValueError('Must supply at least one query.')

	if inputs is not None:
		inputs = list(map(QueryInput.convert, inputs))
		if len(inputs) != len(queries):
			raise ValueError('Number of inputs does not match number of queries.')
	else:
		inputs = [QueryInput(str(i + 1)) for i in range(len(queries))]

	items = [_query_single(db, query, input) for query, input in zip(queries, inputs)]
	return QueryResults(items=items, genomeset=db.genomeset, signaturesmeta=db.signatures_meta)


def _query_single(db: MIDASDatabase, sig: np.ndarray, input: QueryInput):
	dists = jaccard_sparse_array(sig, db.genome_signatures, distance=True)
	matches = find_matches(zip(db.genomes, dists))
	consensus, others = consensus_taxon(matches.keys())

	# No matches found
	if not matches:
		return QueryResultItem(
			input=input,
			success=True,
			predicted_taxon=None,
			report_taxon=None,
		)

	item = QueryResultItem(
		input=input,
		success=True,
		predicted_taxon=consensus,
		report_taxon=None if consensus is None else reportable_taxon(consensus),
	)

	# Warn of inconsistent matches
	if others:
		msg = f'Query matched {len(others)} inconsistent taxa: '
		msg += ', '.join(map(_taxon_repr, others))
		msg += '. Reporting lowest common ancestor of this set.'
		item.warnings.append(msg)

	# No consensus found - matches do not have common ancestor
	if consensus is None:
		item.success = False
		item.error = 'Matched taxa have no common ancestor.'

	# Consensus has no reportable ancestor
	elif item.report_taxon is None:
		item.success = False
		item.error = (
			f'Matched taxon {_taxon_repr(consensus)} has no reportable ancestor. '
			'This indicates a problem with the database.'
		)

	return item
