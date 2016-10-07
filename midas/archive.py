import io
from zipfile import ZipFile
import json
from distutils.version import LooseVersion


class DatabaseArchive:
	"""Zipped file for storing and distributing MIDAS database data"""

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
		return [ f.split('/', 1)[1] for f in self._zipfile.namelist()
		         if f.startswith('genomes/') ]

	def has_genome(self, key):
		try:
			self._zipfile.getinfo(self._genome_path(key))
		except KeyError:
			return False
		else:
			return True

	def store_genome(self, genome):
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
		with self._open_text(self._genome_path(key)) as fobj:
			data = json.load(fobj)

		if class_ is None:
			return data
		else:
			return class_.from_json(data)

	def all_genomes(self, class_=None):
		for key in self.list_genomes():
			yield self.get_genome(key, class_=class_)

	def list_genome_sets(self):
		return [ f.split('/', 1)[1] for f in self._zipfile.namelist()
		         if f.startswith('genome_sets/') ]

	def has_genome_set(self, key):
		try:
			self._zipfile.getinfo(self._genome_set_path(key))
		except KeyError:
			return False
		else:
			return True

	def store_genome_set(self, genome_set, store_genomes=False):
		if genome_set.key is None or genome_set.key_version is None:
			raise TypeError('Genome set key and version must be present')

		path = self._genome_set_path(genome_set.key)
		if self.has_genome(genome_set.key):
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
		with self._open_text(self._genome_set_path(key)) as fobj:
			data = json.load(fobj)

		annotations_dict = data.pop('annotations')

		if classes is None:
			return data, annotations_dict
		else:
			gset_class, annotations_class = classes
			gset = gset_class.from_json(data)
			annotations = { k: annotations_class.from_json(v)
			                for k, v in annotations_dict.items() }
			return gset, annotations

	def all_genome_sets(self, classes=None):
		for key in self.list_genome_sets():
			yield self.get_genome_set(key, classes=classes)

	def extract(self, db, session=None):
		"""Extract archive into database"""

		# Create session if needed (remember to commit after)
		if session is None:
			commit_after = True
			session = db.get_session()
		else:
			commit_after = False

		# Extract genomes
		for genome_data in self.all_genomes():
			existing = session.query(db.Genome).\
			           filter_by(key=genome_data['key']).\
			           scalar()

			if existing is None:
				# Create new
				genome = db.Genome.from_json(genome_data)
				session.add(genome)

			else:
				# Update existing
				current_version = LooseVersion(existing.key_version)
				new_version = LooseVersion(genome_data['key_version'])

				if current_version < new_version:
					existing.update_from_json(json_data)
				else:
					pass # TODO - warn?

		# Extract genome sets
		for gset_data, annotations_data_dicts in self.all_genome_sets():

			existing = session.query(db.GenomeSet).\
			           filter_by(key=gset_data['key']).\
			           scalar()

			if existing is None:
				# Create new
				gset = db.GenomeSet.from_json(gset_data)
				session.add(gset)

			else:
				# Update existing
				current_version = LooseVersion(existing.key_version)
				new_version = LooseVersion(gset_data['key_version'])

				if current_version >= new_version:
					continue # TODO - warn?

				existing.update_from_json(gset_data)
				gset = existing

			# Reset exising annotations (if any) and recreate
			gset.annotations = []
			for genome_key, annotations_data in annotations_data_dicts.items():

				genome = session.query(db.Genome).\
				         filter_by(key=genome_key).\
				         scalar()

				if genome is None:
					raise RuntimeError('Genome key {} not found in database'
					                   .format(genome_key))

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
		zf = ZipFile(file, 'x' if overwrite else 'w')

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
