#!/usr/bin/env python3

import os
import logging

import click
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord

from astroplan import Observer
from astroplan.constraints import observability_table

from casatools import msmetadata, measures, quanta
msmd = msmetadata()
qa = quanta()
me = measures()

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
def main(ms, field):
    """
    Helper script to determine the range of parallactic angles in an input MS.

    Part of the plumber package to generate full Stokes beam models.
    """

    msmd.open(ms)
    telescope_name = msmd.observatorynames()
    telescope_pos = msmd.observatoryposition()

    # In coordinates of rad, rad, m
    telescope_pos = [np.rad2deg(telescope_pos['m0']['value']),
                     np.rad2deg(telescope_pos['m1']['value']),
                     telescope_pos['m2']['value']]

    logger.info(f"Telescope {telescope_name[0]} is at co-ordinates "
                f"{telescope_pos[0]:.2f} deg, "
                f"{telescope_pos[1]:.2f} deg, "
                f"{telescope_pos[2]:.2f} m")

    field_names = msmd.fieldnames()
    target_fields = msmd.fieldsforintent('TARGET', asnames=True)

    if len(field_names) == 0:
        raise ValueError(
            "No fields found in input MS. Please check your inputs.")

    if field is None:
        if len(target_fields) == 0:
            logger.warning(
                f'No fields found with intent TARGET. Using field 0 {field_names[0]}')
            logger.warning(
                'If this is not desired, pass in --field on the command line.')
            field = field_names[0]

        elif len(target_fields) > 1:
            logger.warning(
                f'Multiple fields found with intent TARGET. Using the first one {target_fields[0]}')
            field = target_fields[0]

        elif len(target_fields) == 1:
            logger.info(f'Using field {target_fields[0]}')
            field = target_fields[0]

    else:
        if field not in field_names:
            outmsg = f"Input field {field} not found in MS. Available fields are {','.join(field_names)}"
            raise ValueError(outmsg)

        logger.info(f"Using field {field}")

    # CASA likes to store the time as mjd seconds
    times = msmd.timesforfield(msmd.fieldsforname(field)[0])
    times = Time(times/(3600.*24), format='mjd')

    ph_centre = msmd.phasecenter(msmd.fieldsforname(field)[0])
    ph_centre = SkyCoord(
        ph_centre['m0']['value']*u.rad, ph_centre['m1']['value']*u.rad, unit='radian')
    logger.info(f"Co-ordinates of source are {ph_centre.to_string('hmsdms')}")

    obs = Observer(longitude=telescope_pos[0]*u.deg, latitude=telescope_pos[1]*u.deg,
                   elevation=telescope_pos[2]*u.m, name=telescope_name[0], timezone='Etc/GMT0')

    parangs = np.rad2deg(obs.parallactic_angle(times, ph_centre))

    parang_beg = parangs[0]
    parang_mid = parangs[parangs.size//2]
    parang_end = parangs[-1]

    logger.info(f"First parallactic angle is : {parang_beg:.2f}")
    logger.info(f"Middle parallactic angle is : {parang_mid:.2f}")
    logger.info(f"Last parallactic angle is : {parang_end:.2f}")

    plt.style.use('seaborn-poster')

    fig, ax = plt.subplots(1, 1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.plot(times.to_datetime(), parangs, 'o')
    fmt = DateFormatter('%Y-%m-%d %H:%M')
    ax.xaxis.set_major_formatter(fmt)
    ax.xaxis.set_tick_params(rotation=45)

    ax.set_xlabel('Date')
    ax.set_ylabel('Parallactic Angle (deg)')

    ms_basename = os.path.basename(ms).strip('.mms')
    outname = f"{ms_basename}_parallactic_angles.png"
    plt.savefig(outname, bbox_inches='tight')
