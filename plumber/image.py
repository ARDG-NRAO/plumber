#!/usr/bin/env python3

from casatools import image
ia = image()


import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)-15s %(levelname)s: %(message)s",
                        level=logging.INFO)

def check_if_cube(summary):
    """
    Check if the input shape is a cube

    Inputs:
    summary     Dict, output of ia.summary()

    Returns:
    is_cube     Flag to indicate if the image is a cube, bool
    """

    shape = summary['shape']
    freqidx = [ff for ff, unit in enumerate(summary['axisunits']) if 'hz' in unit.lower()][0]
    stokesidx = [ff for ff, unit in enumerate(summary['axisunits']) if len(unit) == 0][0]

    if shape[freqidx] > 1 or shape[stokesidx] > 1:
        return True
    else:
        return False


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

    ia.open(imagename)
    shape = ia.shape()
    csys = ia.coordsys().torecord()
    summary = ia.summary(list=False, verbose=False)
    ia.close()

    is_cube = check_if_cube(summary)
    imsize = [shape[0], shape[1]]

    freqidx = [ff for ff, unit in enumerate(summary['axisunits']) if 'hz' in unit.lower()][0]
    reffreq = summary['refval'][freqidx]

    if summary['axisunits'][freqidx] == 'Hz':
        factor = 1e6
    elif summary['axisunits'][freqidx] == 'MHz':
        factor = 1
    elif summary['axisunits'][freqidx] == 'GHz':
        factor = 1e-3

    reffreq /= factor

    return imsize, reffreq, is_cube
