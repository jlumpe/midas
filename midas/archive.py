import io
from zipfile import ZipFile
import json
from distutils.version import LooseVersion


class DatabaseArchive:
	"""Zipped file for storing and distributing MIDAS database data.

	Data is stored internally in JSON format, and may be retrieved as parsed
	JSON (e.g. ``dict``s) or as ORM model instances. Storing data requires the
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

		Returns
		-------
		list of str
		"""
		return [
			f.split('/', 1)[1] for f in self._zipfile.namelist()
			if f.startswith('genomes/')
		]

	def has_genome(self, key):
		"""Check if genome with ``key`` exists in archive.

		Parameters
		----------
		key : str

		Returns
		-------
		bool
		"""
		try:
			self._zipfile.getinfo(self._genome_path(key))
		except KeyError:
			return False
		else:
			return True

	def store_genome(self, genome):
		"""Store a genome in the archive.

		Parameters
		----------
		genome : midas.database.base.Genome

		Raises
		------
		KeyError
			If a genome with the same key already exists in the archive.
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

		Parameters
		----------
		key :
		class_ : type
			ORM model class to return instance of, if any (subclass of
			:type:`midas.database.base.Genome`).

		Returns
		-------
		ORM model instance if ``class_=True``, otherwise its JSON representation.

		Raises
		------
		KeyError
			If no genome with the given key exists.
		"""
		with self._open_text(self._genome_path(key)) as fobj:
			data = json.load(fobj)

		if class_ is None:
			return data
		else:
			return class_.from_json(data)

	def all_genomes(self, class_=None):
		"""Iterate over all genomes stored in the archive.

		Parameters
		----------
		class_ : type
			ORM model class to yield instances of, if any (subclass of
			:type:`midas.database.base.Genome`).

		Returns
		-------
		Iterator over genome data.

		See Also
		--------
		:method:`.get_genome`
		"""
		for key in self.list_genomes():
			yield self.get_genome(key, class_=class_)

	def list_genome_sets(self):
		"""List keys of all genome sets stored in the archive.

		Returns
		-------
		list of str
		"""
		return [
			f.split('/', 1)[1] for f in self._zipfile.namelist()
			if f.startswith('genome_sets/')
		]

	def has_genome_set(self, key):
		"""Check if a genome set with the given key is stored in the archive.

		Parameters
		----------
		key : str

		Returns
		-------
		bool
		"""
		try:
			self._zipfile.getinfo(self._genome_set_path(key))
		except KeyError:
			return False
		else:
			return True

	def store_genome_set(self, genome_set, store_genomes=False):
		"""Store a genome set in the archive.

		Parameters
		----------
		genome_set : midas.database.base.GenomeSet
		store_genomes : bool
			Whether to also store corresponding genome objects for all
			of the genome set's annotations.

		Raises
		------
		KeyError
			If a genome set with the given key already exists in the archive.
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

		Parameters
		----------
		key : str
		classes :
			Optional, 2-tuple of ORM classes for genome set and genome
			annotations, e.g.``(GenomeSet, GenomeAnnotations)``.

		Returns
		-------
			``(genome_set, annotations)`` pair. If ``classes`` is given,
			the first item is a ``GenomeSet`` instance and the second is a
			mapping from genome keys to ``GenomeAnnotations`` instances. If
			``classes`` is omitted JSON data is substituted for the ORM model
			instances.

		Raises
		------
		KeyError
			If no genome with the given key exists.
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

		Parameters
		----------
		classes :
			Optional, 2-tuple of ORM classes for genome set and genome
			annotations, e.g.``(GenomeSet, GenomeAnnotations)``.

		Returns
		-------
			Iterator over ``(metadata, annotations)`` pairs.

		See Also
		--------
		:meth:`.get_genome_set`
		"""
		for key in self.list_genome_sets():
			yield self.get_genome_set(key, classes=classes)

	def extract(self, db, session=None):
		"""Extract archive into MIDAS database.

		Parameters
		----------
		db : midas.database.base.BasicDatabase
		session : sqlalchemy.orm.session.Session
			SQLAlchemy session to use.
		"""

		# Create session if needed (remember to commit after)
		if session is None:
			commit_after = True
			session = db.get_session()
		else:
			commit_after = False

		# Extract genomes
		for genome_data in self.all_genomes():
			existing = session.query(db.Genome)\
				.filter_by(key=genome_data['key'])\
				.scalar()

			if existing is None:
				# Create new
				genome = db.Genome.from_json(genome_data)
				session.add(genome)

			else:
				# Update existing
				current_version = LooseVersion(existing.key_version)
				new_version = LooseVersion(genome_data['key_version'])

				if current_version < new_version:
					existing.update_from_json(genome_data)
				else:
					pass  # TODO - warn?

		# Extract genome sets
		for gset_data, annotations_data_dicts in self.all_genome_sets():

			existing = session.query(db.GenomeSet)\
				.filter_by(key=gset_data['key'])\
				.scalar()

			if existing is None:
				# Create new
				gset = db.GenomeSet.from_json(gset_data)
				session.add(gset)

			else:
				# Update existing
				current_version = LooseVersion(existing.key_version)
				new_version = LooseVersion(gset_data['key_version'])

				if current_version >= new_version:
					continue  # TODO - warn?

				existing.update_from_json(gset_data)
				gset = existing

			# Reset exising annotations (if any) and recreate
			gset.annotations = []
			for genome_key, annotations_data in annotations_data_dicts.items():

				genome = session.query(db.Genome)\
					.filter_by(key=genome_key)\
					.scalar()

				if genome is None:
					raise RuntimeError(
						'Genome key {} not found in database'
						.format(genome_key)
					)

				annotations = db.GenomeAnnotations.from_json(annotations_data)
				annotations.genome = genome

				gset.annotations.append(annotations)

		if commit_after:
			session.commit()

	def _open_text(self, name):
		"""Open archive file in test mode"""
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

		Parameters
		----------
		file : str
			Path to file to create.
		overwrite : bool
			If False, raise an exception if file already exists.

		Returns
		-------
		.DatabaseArchive

		Raises
		------
		FileExistsError
			If ``file`` exists and ``overwrite`` is false.
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
