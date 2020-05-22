"""Tests for midas.archive module"""

import pytest

import midas.database.basicdatabase as basicdb
from midas.archive import DatabaseArchive
from midas.database.io import extract_archive, extract_genomes, extract_genome_set


genome_attrs = [
	('key', False),
	('key_version', False),
	('description', False),
	('is_assembled', False),
	('gb_db', False),
	('gb_id', False),
	('gb_acc', False),
	('gb_taxid', False),
	('gb_summary', True),
	('gb_tax_summary', True),
	('meta', True),
]

genome_set_attrs = [
	('key', False),
	('key_version', False),
	('name', False),
	('description', False),
	('meta', True),
]

annotations_attrs = [
	('organism', False),
	('tax_species', False),
	('tax_genus', False),
	('tax_strain', False),
]


def json_equal(d1, d2):
	"""Check parsed JSON values for equality."""
	for data in (d1, d2):
		if not isinstance(data, (str, int, float, list, dict)):
			raise TypeError('{} is not a JSONable type'.format(type(data)))

	if isinstance(d1, dict) and isinstance(d2, dict):
		for k, v in d1.items():
			if k not in d2 or not json_equal(v, d2[k]):
				return False

		return True

	elif type(d1) != type(d2):
		return False

	else:
		return d1 == d2


def obj_equal(o1, o2, attrs):
	"""Check if two ORM object instances ar equal."""
	for name, is_sqla_json in attrs:
		val1 = getattr(o1, name)
		val2 = getattr(o2, name)

		if is_sqla_json:
			if not json_equal(val1.as_builtin(), val2.as_builtin()):
				return False
		elif val1 != val2:
			return False

	return True


def obj_json_equal(obj, jsondata, attrs):
	"""Check if ORM object instance is equivalent to JSON data."""
	for name, is_sqla_json in attrs:
		val1 = getattr(obj, name)
		val2 = jsondata[name]

		if is_sqla_json:
			if not json_equal(val1.as_builtin(), val2):
				return False
		elif val1 != val2:
			return False

	return True


@pytest.fixture()
def genomes():
	"""Genomes to store and retrieve from the archive"""
	genomes = []
	for i in range(20):
		genome = basicdb.Genome(
			description='test genome {}'.format(i + 1),
			key='test/{}'.format(i + 1),
			key_version='1.0',
			is_assembled=i % 2 == 0,
			gb_db='assembly',
			gb_id=i + 1,
			gb_acc='GCF_{:08d}.1'.format(i + 1),
			gb_taxid=i + 1,
		)

		genome.gb_summary = dict(
			database=genome.gb_db,
			uid=genome.gb_id,
			accession=genome.gb_acc,
		)

		genome.gb_tax_summary = dict(
			uid=genome.gb_taxid,
		)

		genome.meta = dict(
			key=genome.key,
		)

		genomes.append(genome)

	return genomes


@pytest.fixture()
def genome_set(genomes):
	genome_set = basicdb.GenomeSet(
		key='test/test_set',
		key_version='1.0',
		name='test set',
		description='some description text',
		meta=dict(
			some_key='some_value',
		),
	)

	for i, genome in enumerate(genomes):
		annotations = basicdb.GenomeAnnotations(
			genome=genome,
			organism='organism {}'.format(i + 1),
			tax_species='species {}'.format(i + 1),
			tax_genus='genus {}'.format(i + 1),
			tax_strain='strain {}'.format(i + 1),
		)
		genome_set.annotations.append(annotations)

	return genome_set


@pytest.fixture(params=[False, True])
def archive(request, tmpdir, genomes, genome_set):
	use_store_all = request.param

	archive_path = tmpdir.\
	               join('archive{}.zip'.format(2 if use_store_all else 1)).\
	               strpath

	with DatabaseArchive.create(archive_path) as archive:

		if use_store_all:
			for genome in genomes:
				archive.store_genome(genome)

			archive.store_genome_set(genome_set)

		else:
			archive.store_genome_set(genome_set, store_genomes=True)

	return DatabaseArchive(archive_path)


class TestGenomeRetrieval:

	def test_list(self, archive, genomes):
		keys = set(g.key for g in genomes)
		stored_keys = set(archive.list_genomes())
		assert keys == stored_keys

	def test_has_genome(self, archive, genomes):

		for genome in genomes:

			assert archive.has_genome(genome.key)

		assert not archive.has_genome('test/notinarchive')

	def test_get_json(self, archive, genomes):

		for genome in genomes:
			json_data = archive.get_genome(genome.key)
			assert obj_json_equal(genome, json_data, genome_attrs)

	def test_get_obj(self, archive, genomes):

		for genome in genomes:
			loaded_genome = archive.get_genome(genome.key, basicdb.Genome)
			assert obj_equal(genome, loaded_genome, genome_attrs)

	# TODO - test overwrite old version


class TestGenomeSetRetrieval:

	def test_list(self, archive, genome_set):
		stored_keys = archive.list_genome_sets()
		assert stored_keys == [genome_set.key]

	def test_has_genome_set(self, archive, genome_set):

		assert archive.has_genome_set(genome_set.key)
		assert not archive.has_genome_set('test/notinarchive')

	def test_get_json(self, archive, genome_set):

		set_json, annotations_json = archive.get_genome_set(genome_set.key)

		assert obj_json_equal(genome_set, set_json, genome_set_attrs)

		for annotations in genome_set.annotations:

			json_data = annotations_json[annotations.genome.key]
			assert obj_json_equal(annotations, json_data, annotations_attrs)

	def test_get_obj(self, archive, genome_set):

		gset_classes = (basicdb.GenomeSet, basicdb.GenomeAnnotations)
		loaded_data = archive.get_genome_set(genome_set.key, gset_classes)
		loaded_set, loaded_annotations_dict = loaded_data

		assert obj_equal(genome_set, loaded_set, genome_set_attrs)

		for annotations in genome_set.annotations:

			loaded_annotations = loaded_annotations_dict[annotations.genome.key]
			assert obj_equal(annotations, loaded_annotations, annotations_attrs)


class TestExtract:

	@pytest.fixture
	def db(self, tmpdir):
		return basicdb.BasicDatabase.create(tmpdir.join('testdb').strpath)


	def test_extract_genomes(self, archive, db, genomes, genome_set):
		# TODO...
		extract_genomes(archive, db)


	def test_extract_genome_set(self, archive, db, genomes, genome_set):
		# TODO...
		extract_genomes(archive, db)
		extract_genome_set(archive, db, genome_set.key)


	def test_extract_all(self, archive, db, genomes, genome_set):
		# TODO...
		extract_archive(archive, db)


	def test_extract_method(self, archive, db, genomes, genome_set):
		# TODO...
		archive.extract(db)
