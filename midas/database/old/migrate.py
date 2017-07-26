"""Migrate old database (1.X) to new (2.X)."""

from types import MappingProxyType

from sqlalchemy.orm import sessionmaker

from midas.database import models


# Map from attributes of new Genome to old
GENOME_ATTR_MAP = MappingProxyType(dict(
	key='key',
	version='key_version',
	description='description',
	is_assembled='is_assembled',
	ncbi_taxid='gb_taxid',
	entrez_summary='gb_summary',
	extra='meta',
	entrez_db='gb_db',
	entrez_id='gb_id',
	genbank_acc='gb_acc',
))

# Map from attributes of new GenomeSet to old
GENOMESET_ATTR_MAP = MappingProxyType(dict(
	key='key',
	version='key_version',
	name='name',
	description='description',
	extra='meta',
))


def copy_attrs(src, dest, attr_map):
	"""Copy attributes from one object to another.

	:param src: Source object.
	:param dest: Dest object.
	:param dict attr_map: Map from names of dest attributes to src.
	"""
	for dest_attr, src_attr in attr_map.items():
		setattr(dest, dest_attr, getattr(src, src_attr))


def migrate_instances(src_objs, dest_model, dest_session, attr_map, progress=False):
	"""Migrate instances to a different database session.

	:param src_objs: Collection of source object instances.
	:param dest_model: SQLAlchemy model class to migrate to.
	:param dest_session: Destination SQLAlchemy session to add migrated
		instances to.
	:param attr_map: Map from names of dest attributes to src.
	:returns: Dict mapping src IDs to dest instances.
	"""

	id_map = dict()

	if progress:
		from tqdm import tqdm
		objs_iter = tqdm(src_objs)
	else:
		objs_iter = src_objs

	for src in objs_iter:

		dest = dest_model()
		copy_attrs(src, dest, attr_map)

		dest_session.add(dest)
		dest_session.flush()

		id_map[src.id] = dest

	return id_map


def migrate_genomes(old_session, old_db, new_session):
	return migrate_instances(
		old_session.query(old_db.Genome),
		models.Genome,
		new_session,
		GENOME_ATTR_MAP,
	)


def migrate_db(old_db, new_engine):

	new_session = sessionmaker(bind=new_engine)()
	old_session = old_db.get_session()

	# Migrate genomes
	genome_map = migrate_genomes(old_session, old_db, new_session)

	# Migrate genome sets
	gset_map = migrate_instances(
		old_session.query(old_db.GenomeSet),
		models.GenomeSet,
		new_session,
		GENOMESET_ATTR_MAP,
	)

	new_session.commit()
