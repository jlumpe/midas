"""Rename ReferenceGenomeSet table and relationships

Revision ID: b1c356705db2
Revises: 1c060eb1fc83
Create Date: 2021-07-06 00:18:45.416430

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'b1c356705db2'
down_revision = '1c060eb1fc83'
branch_labels = None
depends_on = None


def upgrade():
	op.rename_table('reference_genome_sets', 'genome_sets')

	with op.batch_alter_table('genome_sets') as batch_op:
		batch_op.drop_index('ix_reference_genome_sets_key')
		batch_op.create_index(None, ['key'], unique=False)

	with op.batch_alter_table('genome_annotations') as batch_op:
		batch_op.alter_column('reference_set_id', new_column_name='genome_set_id')

	with op.batch_alter_table('taxa') as batch_op:
		batch_op.drop_index('ix_taxa_reference_set_id')
		batch_op.alter_column('reference_set_id', new_column_name='genome_set_id')

	op.create_index(None, 'taxa', ['genome_set_id'], unique=False)

	# op.drop_constraint(None, 'genome_annotations', type_='foreignkey')
	# op.create_foreign_key(None, 'genome_annotations', 'genome_sets', ['genome_set_id'], ['id'], ondelete='CASCADE')
	# op.create_index(op.f('ix_taxa_genome_set_id'), 'taxa', ['genome_set_id'], unique=False)
	# op.drop_constraint(None, 'taxa', type_='foreignkey')
	# op.create_foreign_key(None, 'taxa', 'genome_sets', ['genome_set_id'], ['id'], ondelete='CASCADE')
	# ### end Alembic commands ###


def downgrade():
	raise NotImplementedError()
