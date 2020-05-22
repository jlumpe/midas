"""Import information into database from other formats."""

from distutils.version import LooseVersion


_ARCHIVE_UPDATE_VALUES = ('raise', 'newer', 'same', 'always', 'skip')


def _extract_genome(data, db, session, update):
	key = data['key']

	existing = session.query(db.Genome).filter_by(key=key).scalar()

	if existing is None:
		# Create new
		genome = db.Genome.from_json(data)
		session.add(genome)
		return True

	if update is False:
		raise KeyError('Genome with key %r already exists in database' % key)

	if update == 'skip':
		return False

	current_version = LooseVersion(existing.key_version)
	new_version = LooseVersion(data['key_version'])

	if update is True and (current_version <= new_version):
		return False
	if update is 'overwrite' and (current_version < new_version):
		return False

	existing.update_from_json(data)
	return True

def extract_genomes(archive, db, update='raise', _session=None):
	"""Extract genomes from an archive file into a MIDAS database.

	:type archive: .DatabaseArchive
	:type db: midas.database.base.AbstractDatabase
	:param update: Action to take when attempting to extract a genome that is
	    already present in the database. ``'raise'`` raises a ``KeyError`` when
	    this situation occurs. ``'skip'`` skips extracting that genome.
	    ``'newer'`` and ``'same'`` update the existing genome's data with the
	    archive data if the archive data's version is greater than or greater
	    than or equal to the existing version, respectively. ``'always'``
	    always updates existing genomes even if it is with an older version.
	:type update: str

	:returns: Dict mapping genome keys to bools indicating which genomes were
		extracted. Those with ``False`` values were skipped due to the update
		policy.
	:rtype: dict[str, bool]
	"""
	if update not in _ARCHIVE_UPDATE_VALUES:
		raise ValueError('update must be one of %s' % ', '.join(map(repr, _ARCHIVE_UPDATE_VALUES)))

	extracted = {}

	with db._optional_session(_session, commit=True) as session:

		for genome_data in archive.all_genomes():
			key = genome_data['key']
			extracted[key] = _extract_genome(genome_data, db, session, update)

	return extracted


def extract_genome_set(archive, db, key, _session=None):
	"""Extract a genome set from an archive file into a MIDAS database.

	:type archive: .DatabaseArchive
	:type db: midas.database.base.AbstractDatabase
	:param key: Key of genome set to extract
	:type key: str
	:param _session:
	:type _session:
	"""

	gset_data, annotations_data_dicts = archive.get_genome_set(key)

	with db._optional_session(_session, commit=True) as session:

		kv = gset_data['key_version']
		existing = session.query(db.GenomeSet) \
			.filter_by(key=key, key_version=kv) \
			.scalar()
		if existing is not None:
			raise KeyError(
				'A genome set with key %r and version %r already exists in the database'
				% (key, kv))

		# Create
		gset = db.GenomeSet.from_json(gset_data)
		session.add(gset)

		# Add annotations
		for gkey, annotations_data in annotations_data_dicts.items():

			genome = session.query(db.Genome).filter_by(key=gkey).scalar()

			if genome is None:
				raise RuntimeError('Genome key %r not found in database' % gkey)

			annotations = db.GenomeAnnotations.from_json(annotations_data)
			annotations.genome = genome

			gset.annotations.append(annotations)


def extract_archive(archive, db, session=None):
	"""Extract all contents of archive into MIDAS database.

	:type archive: .DatabaseArchive
	:type db: midas.database.base.AbstractDatabase
	:param session: SQLAlchemy session to use.
	:type session: sqlalchemy.orm.session.Session
	"""
	with db._optional_session(session, commit=True) as session:
		extract_genomes(archive, db, _session=session)

		for key in archive.list_genome_sets():
			extract_genome_set(archive, db, key, _session=session)
