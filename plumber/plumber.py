#!/usr/bin/env python3

import click
from plumber.image import parse_image
from plumber.zernike import get_zcoeffs, zernikeBeam
#from plumber.parang_finder import parallctic_angle

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
@click.argument('imagename', type=click.Path(exists=True))
@click.argument('CSV', type=click.File())
@click.option('-a', '--padding', type=int, default=8, help='Padding factor for aperture, affects smoothness of output beam', show_default=True)
@click.option('-d', '--dish_dia', type=float, default=None, help='Diameter of the antenna dish. If not one of VLA, ALMA, MeerKAT or GMRT, must be specified.')
@click.option('-l', '--islinear', is_flag=True, help='Specifies if the telescope has linear feeds. If not one of VLA, ALMA, MeerKAT or GMRT, must be specified')
@click.option('-I', '--stokesI', is_flag=True, help='Only generate the Stokes I beam, not the full Stokes beams')
@click.option('-P', '--parallel', is_flag=True, help='Use parallel processing (no MPI) to speed things up')
@click.option('-p', '--parang', type=float, multiple=True, help='Beginning (and optionally end) parallactic angle for the PB. Pass --parang twice to specify beg and end.', show_default=True)
def main(imagename, csv, padding, dish_dia, islinear, stokesi, parallel, parang):
    """
    Given the input image and the coefficient CSV file, generate the full Stokes
    primary beam at the image centre frequency. If the input is a cube, a
    corresponding output cube of PBs will be generated, for each of the
    Stokes-I, -Q, -U and -V beams.

    To pass in a beginning and end parallactic angle, pass in --parang twice,
    like

    plumber <image> <CSV> --parang 20 --parang 110

    The CSV file must be of the format

    #stokes,freq,ind,real,imag,eta

    The eta column in the CSV is optional.
    """

    parang = sorted(parang)
    if len(parang) > 2:
        raise ValueError(f"Either pass in a single PA or two values of PA. Currently set to {parang}")

    imsize, imfreq, is_cube = parse_image(imagename)
    zdf, zfreq, nstokes = get_zcoeffs(csv, imfreq)

    logger.info(f"Image is at {imfreq:.2f} MHz. Model PB will be generated at {zfreq:.2f} MHz")

    zb = zernikeBeam()
    zb.initialize(zdf, imagename, padfac=padding, dish_dia=dish_dia,
                        islinear=islinear, stokesi=stokesi, parang=parang, parallel=parallel)

    jones_beams = zb.gen_jones_beams()
    stokes_beams = zb.jones_to_mueller(jones_beams)

    zb.regrid_to_template(stokes_beams, imagename)

    #if parang is not None:
        #stokes_beams = zb.rotate_beam(stokes_beams, parang)
