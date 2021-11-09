
#!/usr/bin/env python3

import os
import logging

import click
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter

from plumber.sky import ParallacticAngle

# Add colours for warnings and errors
logging.addLevelName(
    logging.WARNING, "\033[1;31m%s\033[1;0m" % logging.getLevelName(logging.WARNING))
logging.addLevelName(
    logging.ERROR, "\033[1;41m%s\033[1;0m" % logging.getLevelName(logging.ERROR))

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

@click.command(context_settings=ctx)
@click.argument('MS', type=click.Path(exists=True))
@click.option('--field', type=str, default=None,
              help='Name of field to consider, must match exactly the name in the MS. '
              'If this is not specified, will use the first field with the '
              'TARGET intent.')
@click.option('--use-astropy', is_flag=True,
              help='Use astropy convention rather than CASA convention to '
              '  calculate the parallactic angle')
def main(ms, field, use_astropy):
    """
    Helper script to determine the range of parallactic angles in an input MS.

    Part of the plumber package to generate full Stokes beam models.
    """

    parang = ParallacticAngle(ms=ms, use_astropy=use_astropy)

    logger.info(f"First parallactic angle is : {parang.parang_start:.2f}")
    logger.info(f"Middle parallactic angle is : {parang.parang_mid:.2f}")
    logger.info(f"Last parallactic angle is : {parang.parang_end:.2f}")

    plt.style.use('seaborn-poster')

    fig, ax = plt.subplots(1, 1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.plot(parang.source_times.to_datetime(), parang.parangs, 'o')
    fmt = DateFormatter('%Y-%m-%d %H:%M')
    ax.xaxis.set_major_formatter(fmt)
    ax.xaxis.set_tick_params(rotation=45)

    ax.set_xlabel('Date')
    ax.set_ylabel('Parallactic Angle (deg)')

    ms_basename = os.path.basename(ms).strip('.mms')
    ms_basename = os.path.basename(ms).strip('.ms')

    outname = f"{ms_basename}_parallactic_angles.png"
    plt.savefig(outname, bbox_inches='tight')
