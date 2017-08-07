"""Test midas.io.seq."""

from io import StringIO

import pytest
import numpy as np
from Bio import Seq, SeqIO

import midas.io.seq as mseqio
from midas.kmers import vec_to_coords
from midas.test import make_kmer_seq


@pytest.mark.parametrize('as_coords', [False, True])
def test_find_kmers_parse(as_coords):
	"""Test the find_kmers_parse function."""

	k = 11

	vec = np.zeros(4 ** k, dtype=bool)
	kspec = None

	records = []

	# Create FASTA file with 10 records
	for i in range(10):

		seq, _kspec, rec_vec = make_kmer_seq(
			100000,
			k=k,
			prefix_len=5,
			kmer_interval=50,
			n_interval=10,
			seed=i,
		)

		# Combine vectors of all sequences
		vec |= rec_vec

		# Just check that our specs are all the same...
		if kspec is None:
			kspec = _kspec
		else:
			assert kspec == _kspec

		# Convert every other sequence to lower case, just to switch things up...
		if i % 2:
			seq = seq.lower()

		# Create the BioPython sequence record object
		records.append(SeqIO.SeqRecord(
			seq=Seq.Seq(seq.decode('ascii')),
			id='SEQ{}'.format(i + 1),
			description='sequence {}'.format(i + 1),
		))

	# Write records to string buffer in FASTA format
	buf = StringIO()
	SeqIO.write(records, buf, 'fasta')
	buf.seek(0)

	# Parse from buffer
	kmers = mseqio.find_kmers_parse(kspec, buf, 'fasta', coords=as_coords)

	if as_coords:
		assert np.array_equal(vec_to_coords(vec), kmers)
	else:
		assert np.array_equal(vec, kmers)
