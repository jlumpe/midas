"""Update key version cols

Revision ID: c43540b80d50
Revises: b1c356705db2
Create Date: 2021-07-06 20:39:52.289431

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'c43540b80d50'
down_revision = 'b1c356705db2'
branch_labels = None
depends_on = None


def upgrade():
	with op.batch_alter_table('genomes') as batch_op:
		batch_op.drop_index('ix_genomes_key')
		batch_op.drop_column('version')
		batch_op.alter_column('key', nullable=False)

	with op.batch_alter_table('genome_sets') as batch_op:
		batch_op.alter_column('key', nullable=False)

	with op.batch_alter_table('taxa') as batch_op:
		# Add key column but don't make non-nullable yet
		batch_op.add_column(sa.Column('key', sa.String()))
		batch_op.create_unique_constraint(None, ['key'])

	# Populate taxon.key with values of taxon.name
	conn = op.get_bind()
	if conn is not None:
		conn.execute('UPDATE taxa SET key = name;')

	# Now make taxon.key non-nullable
	with op.batch_alter_table('taxa') as batch_op:
		batch_op.alter_column('key', nullable=False)


def downgrade():
	raise NotImplementedError()
