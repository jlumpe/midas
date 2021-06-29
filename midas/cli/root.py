"""Define the root CLI command group."""

import click

from midas import __version__ as MIDAS_VERSION
from .context import CLIContext


# Top-level cli group
@click.group()
@click.option(
	'-d', '--db', 'db_path',
	type=click.Path(exists=True, file_okay=False),
	envvar='MIDAS_DB_PATH',
	help='Directory containing MIDAS database files.'
)
@click.version_option(MIDAS_VERSION, prog_name='midas')
@click.pass_context
def cli(ctx: click.Context, db_path):
	"""Tool for rapid taxonomic identification of microbial pathogens from genomic data."""
	ctx.obj = CLIContext(db_path)
