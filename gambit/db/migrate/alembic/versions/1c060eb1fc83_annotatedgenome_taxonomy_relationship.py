"""AnnotatedGenome taxonomy relationship

Revision ID: 1c060eb1fc83
Revises: 7c1a8837b74d
Create Date: 2021-07-06 00:08:06.845391

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '1c060eb1fc83'
down_revision = '7c1a8837b74d'
branch_labels = None
depends_on = None


def upgrade():
    op.drop_table('annotations_additional_tax_assoc')

    with op.batch_alter_table('genome_annotations') as batch_op:
        batch_op.drop_index('ix_genome_annotations_primary_taxon_id')

        batch_op.alter_column('primary_taxon_id', new_column_name='taxon_id')

    with op.batch_alter_table('genome_annotations') as batch_op:
        batch_op.create_index(None, ['taxon_id'], unique=False)


def downgrade():
    raise NotImplementedError()
