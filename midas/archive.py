import io
from zipfile import ZipFile
import json
import warnings

from midas.database.io import extract_archive


class DatabaseArchive:
	"""Zipped file for storing and distributing MIDAS database data.

	Data is stored internally in JSON format, and may be retrieved as parsed
	JSON (e.g. ``dict``\\ s) or as ORM model instances. Storing data requires the
	ORM model format.
	"""

	_CURRENT_FORMAT_VERSION = '1.0'

	def __init__(self, file):
		self._zipfile = ZipFile(file, 'r')

		with self._open_text('info') as fobj:
			info = json.load(fobj)
		if info['archive_version'] != self._CURRENT_FORMAT_VERSION:
			raise IOError('Archive format is not the current version')

	@property
	def writable(self):
		return self._zipfile.mode != 'r'

	def list_genomes(self):
		"""List keys of all genomes stored in archive.

		:rtype: list[str]
		"""
		return [
			f.split('/', 1)[1] for f in self._zipfile.namelist()
			if f.startswith('genomes/')
		]

	def has_genome(self, key):
		"""Check if genome with ``key`` exists in archive.

		:type key: str
		:rtype: bool
		"""
		try:
			self._zipfile.getinfo(self._genome_path(key))
		except KeyError:
			return False
		else:
			return True

	def store_genome(self, genome):
		"""Store a genome in the archive.

		:type genome: midas.database.base.Genome

		:raises KeyError: If a genome with the same key already exists in the archive.
		"""
		if genome.key is None or genome.key_version is None:
			raise TypeError('Genome key and version must be present')

		path = self._genome_path(genome.key)
		if self.has_genome(genome.key):
			raise KeyError(
				'Genome already exists in archive with key "{}"'
				.format(genome.key)
			)

		data = json.dumps(genome.to_json())
		self._zipfile.writestr(path, data)

	def get_genome(self, key, class_=None):
		"""Get the genome in the archive with a given key.

		:param str key:
		:param type class_: ORM model class to return instance of, if any (subclass of
			:class:`midas.database.base.Genome`).

		:returns: ORM model instance if ``class_=True``, otherwise its JSON representation.
		:raises KeyError: If no genome with the given key exists.
		"""
		with self._open_text(self._genome_path(key)) as fobj:
			data = json.load(fobj)

		if class_ is None:
			return data
		else:
			return class_.from_json(data)

	def all_genomes(self, class_=None):
		"""Iterate over all genomes stored in the archive.

		:param type class_: ORM model class to yield instances of, if any (subclass of
			:class:`midas.database.base.Genome`).

		:returns: Iterator over genome data.

		See also: :meth:`.get_genome`
		"""
		for key in self.list_genomes():
			yield self.get_genome(key, class_=class_)

	def list_genome_sets(self):
		"""List keys of all genome sets stored in the archive.

		:rtype: list[str]
		"""
		return [
			f.split('/', 1)[1] for f in self._zipfile.namelist()
			if f.startswith('genome_sets/')
		]

	def has_genome_set(self, key):
		"""Check if a genome set with the given key is stored in the archive.

		:param str key:
		:rtype: bool
		"""
		try:
			self._zipfile.getinfo(self._genome_set_path(key))
		except KeyError:
			return False
		else:
			return True

	def store_genome_set(self, genome_set, store_genomes=False):
		"""Store a genome set in the archive.

		:type genome_set: midas.database.base.GenomeSet
		:param bool store_genomes: Whether to also store corresponding genome objects for all
			of the genome set's annotations.

		:raises KeyError: If a genome set with the given key already exists in the archive.
		"""
		if genome_set.key is None or genome_set.key_version is None:
			raise TypeError('Genome set key and version must be present')

		path = self._genome_set_path(genome_set.key)
		if self.has_genome_set(genome_set.key):
			raise KeyError(
				'Genome set already exists in archive with key "{}"'
				.format(genome_set.key)
			)

		json_data = genome_set.to_json()
		annotations_dict = json_data['annotations'] = dict()

		for annotations in genome_set.annotations:
			annotations_dict[annotations.genome.key] = annotations.to_json()

		self._zipfile.writestr(path, json.dumps(json_data))

		if store_genomes:
			for genome in genome_set.genomes:
				self.store_genome(genome)

	def get_genome_set(self, key, classes=None):
		"""Get a stored genome set by key.

		:type key: str
		:param classes: Optional, 2-tuple of ORM classes for genome set and genome
			annotations, e.g.``(GenomeSet, GenomeAnnotations)``.

		:returns: ``(genome_set, annotations)`` pair. If ``classes`` is given,
			the first item is a ``GenomeSet`` instance and the second is a
			mapping from genome keys to ``GenomeAnnotations`` instances. If
			``classes`` is omitted JSON data is substituted for the ORM model
			instances.

		:raises KeyError: If no genome with the given key exists.
		"""
		with self._open_text(self._genome_set_path(key)) as fobj:
			data = json.load(fobj)

		annotations_dict = data.pop('annotations')

		if classes is None:
			return data, annotations_dict
		else:
			gset_class, annotations_class = classes
			gset = gset_class.from_json(data)
			annotations = {
				k: annotations_class.from_json(v)
				for k, v in annotations_dict.items()
			}
			return gset, annotations

	def all_genome_sets(self, classes=None):
		"""Iterate through all stored genome sets.

		:param classes: Optional, 2-tuple of ORM classes for genome set and genome
			annotations, e.g.``(GenomeSet, GenomeAnnotations)``.

		:returns: Iterator over ``(metadata, annotations)`` pairs.

		See also: :meth:`.get_genome_set`
		"""
		for key in self.list_genome_sets():
			yield self.get_genome_set(key, classes=classes)

	def extract(self, db, session=None):
		"""Extract all contents of archive into MIDAS database.

		This method is deprecated, use :func:`midas.database.io.extract_archive` instead.

		:type archive: .DatabaseArchive
		:type db: midas.database.base.AbstractDatabase
		:param session: SQLAlchemy session to use.
		:type session: sqlalchemy.orm.session.Session
		"""
		warnings.warn(
			'DatabaseArchive.extract() is deprecated, use midas.database.io.extract_archive() instead.',
			DeprecationWarning
		)
		extract_archive(self, db, session)

	def _open_text(self, name):
		"""Open archive file in text mode"""
		return io.TextIOWrapper(self._zipfile.open(name))

	@classmethod
	def _genome_path(cls, key):
		return 'genomes/' + key

	@classmethod
	def _genome_set_path(cls, key):
		return 'genome_sets/' + key

	@classmethod
	def create(cls, file, overwrite=False):
		"""Create a new empty archive file.

		:param file: Path to file to create.
		:type file: str
		:param overwrite: If False, raise an exception if file already exists.
		:type overwrite: bool

		:rtype: .DatabaseArchive

		:raises FileExistsError: If ``file`` exists and ``overwrite`` is false.
		"""
		zf = ZipFile(file, 'w' if overwrite else 'x')

		info = dict(
			archive_version=cls._CURRENT_FORMAT_VERSION
		)
		zf.writestr('info', json.dumps(info))

		archive = DatabaseArchive.__new__(cls)
		archive._zipfile = zf
		return archive

	def __enter__(self):
		self._zipfile.__enter__()
		return self

	def __exit__(self, *args):
		self._zipfile.__exit__(*args)
