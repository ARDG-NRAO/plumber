#!/usr/bin/env python3

import click
from casatools import table, msmetadata
from casatools import ms as MS

tb = table()
msmd = msmetadata()
ms = MS()

import logging

# Add colours for warnings and errors
logging.addLevelName(logging.WARNING, "\033[1;31m%s\033[1;0m" % logging.getLevelName(logging.WARNING))
logging.addLevelName(logging.ERROR, "\033[1;41m%s\033[1;0m" % logging.getLevelName(logging.ERROR))

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)-20s %(levelname)-8s %(message)s',
    handlers=[
        logging.FileHandler("plumber.log"),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger()

ctx = dict(help_option_names=['-h', '--help'])
@click.command(context_settings = ctx)
@click.argument('MS', type=click.Path(exists=True))
@click.option('--field', type=str, default=None,
              help='Name of field to consider, must match exactly the name in the MS. '
              'If this is not specified, will use the first field with the '
              'TARGET intent.')
def main(ms, field):
    """
    Helper script to determine the range of parallactic angles in an input MS.

    Part of the plumber package to generate full Stokes beam models.
    """

    msmd.open(ms)
    telescope_name = msmd.observatorynames()
    telescope_pos = msmd.observatoryposition()

    field_names = msmd.fieldnames()
    field_dirs = msmd.sourcedirs()
    target_fields = msmd.fieldsforintent('TARGET', asnames=True)

    msmd.close()

    print(telescope_name, telescope_pos, field_names, field_dirs, target_fields)
