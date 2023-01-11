
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
@click.option('--calcweights', '-c', is_flag=True,
              help='Print a list of parallactic angle and associated weights')
@click.option('--binwidth', '-w', type=int, default=5,
              help='Bin width in degrees')
def main(ms, field, use_astropy, calcweights, binwidth):
    """
    Helper script to determine the range of parallactic angles in an input MS.

    Part of the plumber package to generate full Stokes beam models.
    """

    parang = ParallacticAngle(ms=ms, use_astropy=use_astropy)

    if calcweights:
        weights, bins = parang.parang_weights(bin_width=binwidth)
        # Create fractional weights - this doesn't change anything in the averaging
        weights = weights/np.sum(weights)

        with open('parang_weights.txt', 'w') as fptr:
            fptr.write("#parang  weight\n")
            for ww, bb in zip(weights, bins):
                fptr.write("%.2f   % .2f\n" % (bb, ww))

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
