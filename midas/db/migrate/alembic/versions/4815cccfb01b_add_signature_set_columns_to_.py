"""Add signature set columns to ReferenceGenomeSet table

Revision ID: 4815cccfb01b
Revises: 8e4286c25e33
Create Date: 2017-11-05 11:28:35.579800

"""
from alembic import op
import sqlalchemy as sa
import midas.db.sqla


# revision identifiers, used by Alembic.
revision = '4815cccfb01b'
down_revision = '8e4286c25e33'
branch_labels = None
depends_on = None


def upgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    with op.batch_alter_table('reference_genome_sets', schema=None) as batch_op:
        batch_op.add_column(sa.Column('signatureset_key', sa.String(), nullable=True))
        batch_op.add_column(sa.Column('signatureset_version', sa.String(), nullable=True))

    # ### end Alembic commands ###


def downgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    with op.batch_alter_table('reference_genome_sets', schema=None) as batch_op:
        batch_op.drop_column('signatureset_version')
        batch_op.drop_column('signatureset_key')

    # ### end Alembic commands ###
