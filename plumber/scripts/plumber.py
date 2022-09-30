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
    level=logging.DEBUG,
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
@click.option('-d', '--dish-dia', type=float, default=None, help='Diameter of the antenna dish. If not one of VLA, ALMA, MeerKAT or GMRT, must be specified.')
@click.option('-l', '--linear', is_flag=True, help='Specifies if the telescope has linear feeds. If not one of VLA, ALMA, MeerKAT or GMRT, must be specified')
@click.option('-c', '--circular', is_flag=True, help='Specifies if the telescope has circular feeds. If not one of VLA, ALMA, MeerKAT or GMRT, must be specified')
@click.option('-I', '--stokesI', is_flag=True, help='Only generate the Stokes I beam, not the full Stokes beams')
@click.option('-P', '--parallel', is_flag=True, help='Use parallel processing (no MPI) to speed things up')
@click.option('-p', '--parang', type=float, default=0, help='Parallactic angle at which to generate the PB', show_default=True)
@click.option('--parang-file', type=click.Path(exists=True), help='Pass a file containing a list of parallactic angles and weights')
@click.option('--scale', nargs=2, type=float, help='X and Y scaling factors for number of pixels', default=[None, None], show_default=True)
@click.option('--scale-cdelt', nargs=2, type=float, help='X and Y scaling factors for cdelt', default=[None, None], show_default=True)
def main(imagename, csv, padding, dish_dia, linear, circular, stokesi, parallel, parang, parang_file, scale, scale_cdelt):
    """
    Given the input image and the coefficient CSV file, generate the full Stokes
    primary beam at the image centre frequency. If the input is a cube, a
    corresponding output cube of PBs will be generated, for each of the
    Stokes-I, -Q, -U and -V beams.

    Specifying a --parang-file turns on parallactic angle averaging. The file
    should contain two columns - the parallactic angle and the weight. This list
    will be used to generate a weighted average beam over the listed parallactic
    angles. For a simple arithmetic mean over parallactic angle, setting all the
    weights to one is sufficient.

    The CSV file must be of the format

    #stokes,freq,ind,real,imag,eta

    The eta column in the CSV is optional.
    """

    islinear = None
    if linear is True:
        islinear = True

    #parang = sorted(parang)
    #if len(parang) > 2:
    #    raise ValueError(f"Either pass in a single PA or two values of PA. Currently set to {parang}")

    imsize, imfreq, is_stokes_cube = parse_image(imagename)
    zdflist, zfreqlist, nstokes = get_zcoeffs(csv, imfreq)

    logger.info(f"Image is at {imfreq[0].value/1e6:.2f} MHz. Model PB will be generated at {zfreqlist[0]:.2f} MHz")
    logger.warn(f"The above frequency is the first channel frequency if the input image is a spectral cube")

    zb = zernikeBeam()

    for zdf in zdflist:
        zb.initialize(zdf, imagename, padfac=padding, dish_dia=dish_dia,
                            islinear=islinear, stokesi=stokesi, parang=parang, parang_file=parang_file, parallel=parallel,
                            scale=scale, scale_cdelt=scale_cdelt)

        jones_beams = zb.gen_jones_beams()
        stokes_beams = zb.jones_to_mueller(jones_beams)

        zb.regrid_to_template(stokes_beams, imagename)

    #if parang is not None:
        #stokes_beams = zb.rotate_beam(stokes_beams, parang)
