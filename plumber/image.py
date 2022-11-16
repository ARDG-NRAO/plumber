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


class MuellerExpression():
    """
    Store the `immath` expressions required to compute the Mueller beams
    (Stokes basis) from the Jones beams (feed basis).

    The expressions must be a list of strings, each element of the list (i.e.,
    each string) must contain the complete `immath` equation for that particular
    Mueller term.

    Given the feed basis and telescope name as input, a pre-defined set of
    Mueller expressions will be defined. However these may be over-written by
    user input.
    """

    def __init__(self, islinear:bool, telescope:str) -> None:
        """
        Initialize the MuellerExpression() class.

        Inputs:
        islinear    Is the feed basis linear (islinear=True) or circular (islinear=False), bool
        telescope   Name of the instrument/facility, str
        """
        self._expr = []
        self.islinear = islinear
        self._telescope = telescope


    @property
    def telescope(self) -> str:
        """
        Set the name of the telescope for which the Mueller terms are needed.
        """

        return self._telescope


    @telescope.setter
    def telescope(self, name:str) -> None:
        """
        Setter method for telescope
        """

        self._telescope = name
        if telescope.lower() not in ['vla', 'alma', 'meerkat']:
            if self.islinear is False:
                warnings.warn(f'The input telescope {name} does not have pre-defined Mueller expressions. Using the equations for VLA.')
            else:
                warnings.warn(f'The input telescope {name} does not have pre-defined Mueller expressions. Using the equations for ALMA.')


    @telescope.deleter
    def telescope(self):
        """
        Deleter method for telescope
        """

        del self._telescope

    @property
    def expr(self):
        """
        The immath expressions to convert from the input Jones beams in feed
        basis to the output Mueller beams in Stokes basis.

        Depending on the input telescope and polarization basis will automatically select the
        expressions. However it can be over-written by a user supplied list of expressions.
        """

        if len(self._expr) == 0:
            if self.islinear:
                # This is probably true for ASKAP as well
                if 'meerkat' in self.telescope.lower():
                    self._expr = [
                        'real(IM0*CONJ(IM0) + IM1*CONJ(IM1) + IM2*CONJ(IM2) + IM3*CONJ(IM3))',
                        'real(-IM0*CONJ(IM0) - IM1*CONJ(IM1) + IM2*CONJ(IM2) + IM3*CONJ(IM3))',
                        'real(-IM0*CONJ(IM2) - IM1*CONJ(IM3) - IM2*CONJ(IM0) - IM3*CONJ(IM1))',
                        'real(1i*(IM0*CONJ(IM2) + IM1*CONJ(IM3) - IM2*CONJ(IM0) - IM3*CONJ(IM1)))'
                    ]
                else: # This works for ALMA
                    self._expr = [
                        'real(CONJ(IM0)*IM0 + CONJ(IM1)*IM1 + CONJ(IM2)*IM2 + CONJ(IM3)*IM3)',
                        'real(CONJ(IM0)*IM0 - CONJ(IM1)*IM1 + CONJ(IM2)*IM2 - CONJ(IM3)*IM3)',
                        'real(CONJ(IM0)*IM1 + CONJ(IM1)*IM0 + CONJ(IM2)*IM3 + CONJ(IM3)*IM2)',
                        'real(1i*(CONJ(IM0)*IM1 - CONJ(IM1)*IM0 + CONJ(IM2)*IM3 - CONJ(IM3)*IM2))'
                    ]
            else:
                self._expr = [
                    'real( CONJ(IM0)*IM0 + CONJ(IM1)*IM1 + CONJ(IM2)*IM2 + CONJ(IM3)*IM3 )/2.',
                    'real( CONJ(IM0)*IM1 + CONJ(IM1)*IM0 + CONJ(IM2)*IM3 + CONJ(IM3)*IM2 )/2.',
                    'real( 1i*(CONJ(IM0)*IM1 + CONJ(IM2)*IM3) - 1i*(CONJ(IM1)*IM0 + CONJ(IM3)*IM2) )/2.',
                    'real( CONJ(IM0)*IM0 - CONJ(IM1)*IM1 + CONJ(IM2)*IM2 - CONJ(IM3)*IM3 )/2.'
                ]

        return self._expr


    @expr.setter
    def expr(self, Sdag_M_S):
        """
        Setter function for expr
        """

        self._expr = Sdag_M_S

    @expr.deleter
    def expr(self):
        """
        Deleter function for expr
        """

        del self._expr



