#!/usr/bin/env python3

from __future__ import annotations

import spectral_cube
print(spectral_cube.__path__)

from spectral_cube import SpectralCube, StokesSpectralCube
from casatools import image
ia = image()

from casatasks import imsubimage, immath

from plumber.misc import wipe_file, make_unique
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
        self._jones_to_mueller_expr = []
        self._leakage_correct_expr = []
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
    def jones_to_mueller_expr(self):
        """
        The immath expressions to convert from the input Jones beams in feed
        basis to the output Mueller beams in Stokes basis.

        Depending on the input telescope and polarization basis will automatically select the
        expressions. However it can be over-written by a user supplied list of expressions.
        """

        if len(self._jones_to_mueller_expr) == 0:
            if self.islinear:
                # This is probably true for ASKAP as well
                if 'meerkat' in self.telescope.lower():
                    self._jones_to_mueller_expr = [
                        'real(IM0*CONJ(IM0) + IM1*CONJ(IM1) + IM2*CONJ(IM2) + IM3*CONJ(IM3))',
                        'real(IM0*CONJ(IM0) + IM1*CONJ(IM1) - IM2*CONJ(IM2) - IM3*CONJ(IM3))',
                        'real(IM0*CONJ(IM2) + IM1*CONJ(IM3) + IM2*CONJ(IM0) + IM3*CONJ(IM1))',
                        'real(1i*(-IM0*CONJ(IM2) - IM1*CONJ(IM3) + IM2*CONJ(IM0) + IM3*CONJ(IM1)))',
                    ]
                    #    'real(IM0*CONJ(IM0) - IM1*CONJ(IM1) + IM2*CONJ(IM2) - IM3*CONJ(IM3))',
                    #    'real(IM0*CONJ(IM0) - IM1*CONJ(IM1) - IM2*CONJ(IM2) + IM3*CONJ(IM3))',
                    #    'real(IM0*CONJ(IM2) - IM1*CONJ(IM3) + IM2*CONJ(IM0) - IM3*CONJ(IM1))',
                    #    'real(1i*(-IM0*CONJ(IM2) + IM1*CONJ(IM3) + IM2*CONJ(IM0) - IM3*CONJ(IM1)))',
                    #    'real(IM0*CONJ(IM1) + IM1*CONJ(IM0) + IM2*CONJ(IM3) + IM3*CONJ(IM2))',
                    #    'real(IM0*CONJ(IM1) + IM1*CONJ(IM0) - IM2*CONJ(IM3) - IM3*CONJ(IM2))',
                    #    'real(IM0*CONJ(IM3) + IM1*CONJ(IM2) + IM2*CONJ(IM1) + IM3*CONJ(IM0))',
                    #    'real(1i*(-IM0*CONJ(IM3) - IM1*CONJ(IM2) + IM2*CONJ(IM1) + IM3*CONJ(IM0)))',
                    #    'real(1i*(IM0*CONJ(IM1) - IM1*CONJ(IM0) + IM2*CONJ(IM3) - IM3*CONJ(IM2)))',
                    #    'real(1i*(IM0*CONJ(IM1) - IM1*CONJ(IM0) - IM2*CONJ(IM3) + IM3*CONJ(IM2)))',
                    #    'real(1i*(IM0*CONJ(IM3) - IM1*CONJ(IM2) + IM2*CONJ(IM1) - IM3*CONJ(IM0)))',
                    #    'real(IM0*CONJ(IM3) - IM1*CONJ(IM2) - IM2*CONJ(IM1) + IM3*CONJ(IM0))',
                    #]
                else: # This works for ALMA
                    self._jones_to_mueller_expr = [
                        'real(CONJ(IM0)*IM0 + CONJ(IM1)*IM1 + CONJ(IM2)*IM2 + CONJ(IM3)*IM3)',
                        'real(CONJ(IM0)*IM0 - CONJ(IM1)*IM1 + CONJ(IM2)*IM2 - CONJ(IM3)*IM3)',
                        'real(CONJ(IM0)*IM1 + CONJ(IM1)*IM0 + CONJ(IM2)*IM3 + CONJ(IM3)*IM2)',
                        'real(1i*(CONJ(IM0)*IM1 - CONJ(IM1)*IM0 + CONJ(IM2)*IM3 - CONJ(IM3)*IM2))'
                    ]
            else:
                self._jones_to_mueller_expr = [
                    'real( CONJ(IM0)*IM0 + CONJ(IM1)*IM1 + CONJ(IM2)*IM2 + CONJ(IM3)*IM3 )/2.',
                    'real( CONJ(IM0)*IM1 + CONJ(IM1)*IM0 + CONJ(IM2)*IM3 + CONJ(IM3)*IM2 )/2.',
                    'real( 1i*(CONJ(IM0)*IM1 + CONJ(IM2)*IM3) - 1i*(CONJ(IM1)*IM0 + CONJ(IM3)*IM2) )/2.',
                    'real( CONJ(IM0)*IM0 - CONJ(IM1)*IM1 + CONJ(IM2)*IM2 - CONJ(IM3)*IM3 )/2.'
                ]

        return self._jones_to_mueller_expr


    @jones_to_mueller_expr.setter
    def jones_to_mueller_expr(self, Sdag_M_S):
        """
        Setter function for expr
        """

        self._jones_to_mueller_expr = Sdag_M_S

    @jones_to_mueller_expr.deleter
    def jones_to_mueller_expr(self):
        """
        Deleter function for expr
        """

        del self._jones_to_mueller_expr


    @property
    def leakage_correct_expr(self):
        """
        The immath expressions to produce the leakage corrected images.

        Each Stokes must be a single image, Stokes cubes will not work.

        Depending on the input telescope and polarization basis will automatically select the
        expressions. However it can be over-written by a user supplied list of expressions.

        The convention is :
        IM0, IM1, IM2, IM3 --> Input Jones beams (J_p, J_pq, J_qp, J_q)
        IM4, IM5, IM6, IM7 --> The original Stokes images (I, Q, U, V)
        """

        if len(self._leakage_correct_expr) == 0:
            if 'meerkat' in self.telescope.lower():
                self._leakage_correct_expr = [
                    'IM0*IM16 + IM1*IM17 + IM18*IM2 + IM19*IM3',
                    'IM16*IM4 + IM17*IM5 + IM18*IM6 + IM19*IM7',
                    'IM10*IM18 + IM11*IM19 + IM16*IM8 + IM17*IM9',
                    'IM12*IM16 + IM13*IM17 + IM14*IM18 + IM15*IM19',
                    ]
            else:
                raise NotImplementedError(f'Only MeerKAT is supported for image plane leakage at the moment.')

        return self._leakage_correct_expr


    @leakage_correct_expr.setter
    def leakage_correct_expr(self, Sdag_M_S):
        """
        Setter function for expr
        """

        self._leakage_correct_expr = Sdag_M_S

    @leakage_correct_expr.deleter
    def leakage_correct_expr(self):
        """
        Deleter function for expr
        """

        del self._leakage_correct_expr



class ImagePlaneCorrect():
    """
    Class to hold and perform all the image plane leakage correction math and
    procedures.
    """

    def __init__(self, islinear, telescope):
        self.parallel = False
        self.islinear = islinear
        self.telescope = telescope
        self.mueller_expr = MuellerExpression(self.islinear, self.telescope)


    def leakage_correct(self, mueller_beams:list[str], templateim:list[str]) -> list[str]:
        """
        Perform image plane leakage correction on the input image.
        At present, the output will be four images with each Stokes plane split out independently.

        Inputs:
        mueller_beams       List of names of the Mueller beam images, list[str]
        templateim          Name(s) of the image to perform the correction on, str

        Returns:
        corrected_ims       List of corrected Stokes images.
        """

        imstokes = ['I', 'Q', 'U', 'V']
        stokes_templates = []

        if len(templateim) == 1:
            # Split out each Stokes image
            for stok in imstokes:
                outname = templateim.replace('.im', f'.Stokes{stok}.im')
                imsubimage(templateim, outfile=outname, stokes=stok)

                stokes_templates.append(outname)
        else:
            stokes_templates = list(templateim)


        stokes_corrected = [stok.replace('.im', '_corrected.im') for stok in stokes_templates]

        print(mueller_beams + stokes_templates)
        print(self.mueller_expr.leakage_correct_expr)
        for idx, stok in enumerate(imstokes):
            wipe_file(stokes_corrected[idx])
            print(stokes_corrected[idx])
            print(self.mueller_expr.leakage_correct_expr[idx])

            immath(imagename = mueller_beams + stokes_templates, outfile=stokes_corrected[idx], expr=self.mueller_expr.leakage_correct_expr[idx])
