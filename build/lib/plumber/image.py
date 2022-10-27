#!/usr/bin/env python3

import spectral_cube
print(spectral_cube.__path__)

from spectral_cube import SpectralCube, StokesSpectralCube
from casatools import image
ia = image()

from astropy.io.registry import IORegistryError

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)-15s %(levelname)s: %(message)s",
                        level=logging.INFO)

def parse_image(imagename):
    """
    Parse metadata from the image, such as imsize, central frequency, whether
    it is a cube etc.

    Inputs:
    imagename           Name of the input image, string

    Returns:
    imsize              Size of the direction coordinates in pixels, array
    reffreq             Reference frequency of the image in MHz, float
    is_cube             Flag to indicate if the image is a cube, bool
    """

    try:
        cube = StokesSpectralCube.read(imagename)
    except IORegistryError:
        # Probably CASA image
        cube = StokesSpectralCube.read(imagename, format='casa_image')

    if len(cube.shape) > 3 and cube.shape[-1] > 1:
        is_stokes_cube = True
    else:
        is_stokes_cube = False

    shape = cube.shape
    imsize = [shape[0], shape[1]]

    # This is always a list
    imfreqs = cube.stokes_data['I'].spectral_axis

    return imsize, imfreqs, is_stokes_cube
