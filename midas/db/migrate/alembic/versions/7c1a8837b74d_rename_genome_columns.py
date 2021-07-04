"""Rename Genome columns

Revision ID: 7c1a8837b74d
Revises: d961d0698083
Create Date: 2021-07-06 00:01:30.975726

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '7c1a8837b74d'
down_revision = 'd961d0698083'
branch_labels = None
depends_on = None


def upgrade():
    with op.batch_alter_table('genomes') as batch_op:
        batch_op.alter_column('entrez_db', new_column_name='ncbi_db')
        batch_op.alter_column('entrez_id', new_column_name='ncbi_id')


def downgrade():
    with op.batch_alter_table('genomes') as batch_op:
        batch_op.alter_column('ncbi_db', new_column_name='entrez_db')
        batch_op.alter_column('ncbi_id', new_column_name='entrez_id')
