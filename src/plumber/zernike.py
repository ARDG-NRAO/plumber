#!/usr/bin/env python3

"""
Module to hold all the function that deal with the Zernike generation or CSV
parsing
"""

from __future__ import annotations

import shutil
import pandas as pd
import numpy as np

from scipy.ndimage import rotate
import astropy.units as u

from copy import deepcopy
from typing import Union

import multiprocessing
from functools import partial, wraps

import logging
logger = logging.getLogger(__name__)
logger.setLevel('INFO')

from plumber.misc import wipe_file, make_unique
from plumber.image import parse_image

from casatasks import immath, imregrid
from casatools import image
ia = image()

# Use up to 4 concunrrent processes
NCPU = min(multiprocessing.cpu_count(), 4)


def get_zcoeffs(csv: str, imfreq: float) -> Union[pd.Dataframe, float, int]:
    """
    Given the input frequency of the image, returns the Pandas dataframe with
    the coefficients corresponding to the input frequency. The frequencies in
    the input CSV file are expected to be in MHz.

    Inputs:
    csv         Input CSV filename, string
    imfreq      Image frequency in MHz, float

    Returns:
    zdf         Dataframe containing the subset of the CSV file which is closest
                to imfreq.
    zfreq       Frequency of the coefficients in MHz, float
    nstokes     Number of stokes in the CSV, integer
    """

    zdflist = []
    freqlist = []

    df = pd.read_csv(csv, skipinitialspace=True)
    nstokes = df['#stokes'].unique().size
    freqs = df['freq'].unique()

    #imfreq is an array astropy.units, so convert to right unit and grab
    #numerical value
    for ifreq in imfreq:
        idx = np.argmin(np.abs(freqs - ifreq.to(u.MHz).value))
        zfreq = freqs[idx]

        zdf = df[df['freq'] == zfreq]

        zdflist.append(zdf)
        freqlist.append(zfreq)

    return zdflist, freqlist, nstokes



class zernikeBeam():
    """
    Generate the model beam on a 2D pixel grid.
    """

    def __init__(self):
        self.freq = None
        self.padfac = None
        self.dish_dia = []
        self.islinear = None
        self.parallel = False
        self.telescope = None

        # Everything with ft_ is in aperture domain
        self.ft_cdelt = [None, None]
        self.ft_cdelt_os  = [None, None]
        self.ftcoords = None

        self.ft_npix = 0
        self.ft_npix_os = 0

        self.lambd = 0
        self.scale = [None, None]
        self.eta = [None, None]

        self.oversamp = 20


    def initialize(self, df, templateim, padfac=8, dish_dia=[], islinear=None, stokesi=False, parang=None, parang_file=None, parallel=False, scale=None):
        """
        Initialize the class with an input DataFrame, and optionally padding
        factor for the FFT.

        Inputs:
        df              Input DataFrame containing only the frequency of interest
        templateim      Name of the template image, string
        padfac          Padding factor for FFT, integer  (default: 8)
        islinear        Determines whether telescope is in circular or linear basis, bool
                        If one of MeerKAT, VLA, ALMA or GMRT automatically determined. (default: None)

        Returns:
        None
        """

        self.df = df
        self.telescope = self.get_telescope(templateim)
        self.eta = [1., 1.]

        self.dish_dia = dish_dia
        self.get_dish_diameter()
        self.islinear = islinear
        self.get_feed_basis()

        self.padfac = int(padfac)
        self.parallel = parallel
        self.parang = parang
        self.parang_file = parang_file

        # In MHz
        self.freq = df['freq'].unique()[0]
        self.do_stokesi_only = stokesi

        # List of floats, to scale in X and Y respectively
        self.scale = scale
        logger.debug(f"Scale is {self.scale}")

        # Creates ftcoords and sets ft_cdelt and ft_cdelt_os
        # Sets ft_npix and ft_npix_os
        self.get_npix_aperture(templateim)


    def get_telescope(self, templateim: str) -> str:
        """
        Get the telescope name from the image header

        Input:
        templateim      Name of the input template image, string

        Returns:
        telescope       Name of the input telescope
        """

        ia.open(templateim)
        csys = ia.coordsys().torecord()
        ia.close()

        telescope = csys['telescope']
        logger.debug(f"Telescope is {telescope}.")

        return csys['telescope']


    def get_feed_basis(self) -> None:
        """
        Get the feed basis, from the telescope name.
        """

        if self.islinear is not None:
            outstr = 'circular' if self.islinear == False else 'linear'
            logger.info(f"Overriding automatic determination of feed basis, using user supplied value of {outstr}")

        if "vla" in self.telescope.lower():
            self.islinear = False
        elif "meerkat" in self.telescope.lower():
            self.islinear = True
        elif "alma" in self.telescope.lower():
            self.islinear = True
        elif "gmrt" in self.telescope.lower() and self.freq < 9e2:
            self.islinear = False
        elif "gmrt" in self.telescope.lower() and self.freq >= 9e2:
            self.islinear = True
        else:
            raise ValueError("Unable to determine feed basis, unknown telescope. "
            "Please pass in the feed basis via the islinear paramter to "
            ".initialize()")


    def get_dish_diameter(self) -> None:
        """
        Get the dish diameter for different instruments.

        Returns:
        None, sets self.dish_dia in place
        """

        if self.dish_dia is not None:
            logger.info(f"Overriding automatic determination of telescope dish "
                        "diameter. Using the input of {self.dish_dia} m")
            return

        if 'vla' in self.telescope.lower():
            self.dish_dia = [25, 25]
        elif 'meerkat' in self.telescope.lower():
            #self.dish_dia = 13.5
            self.dish_dia = [13.5, 13.5]
        elif 'alma' in self.telescope.lower():
            self.dish_dia = [12, 12]
        elif 'gmrt' in self.telescope.lower():
            self.dish_dia = [45, 45]
        else:
            raise ValueError('Unknown telescope type. Please initialize the '
                             'class with the dish_diameter value in metres.')

        logger.info(f"Using dish diameter of ({self.dish_dia[0]}m, {self.dish_dia[1]}m) for telescope {self.telescope}.")



    def set_eta_values(self) -> None:
        """
        Set the eta scaling values for the aperture at the selected frequency.

        It can be a single value (for a circular aperture), a separate value for
        each dimension (like for an offset Gregorian), or none.

        If the `--scale` option is passed in on the command line, over-ride all
        pre-existing eta values.
        """

        # Set eta values
        if all(self.scale):
            self.eta = deepcopy(self.scale)
        else:
            if 'eta' in self.df.columns:
                eta = self.df['eta'].unique()[0]
                self.eta = [eta, eta]
                logger.info(f"Using eta value of {self.eta}.")
            elif 'etax' in self.df.columns and 'etay' in self.df.columns:
                self.eta = [self.df['etax'].unique()[0], self.df['etay'].unique()[0]]
                logger.info(f"Using eta value of {self.eta}.")
            else:
                logger.info("Eta column(s) do not exist in the input CSV file, not "
                            "applying frequency dependent scaling.")
                self.eta = [1.0, 1.0]



    def get_npix_aperture(self, templateim: str) -> Union[int, float]:
        """
        Calculate the number of pixels across the aperture from the template
        image.

        Inputs:
        templateim      Name of the input template image, str

        Returns:
        npix            The number of pixels across the aperture, int
        cdelt           The pixel delta in lambda, float
        """

        imsize, imfreq, is_stokes_cube = parse_image(templateim)
        freq = imfreq[0].value

        ia.open(templateim)
        template_csys = ia.coordsys().torecord()
        self.template_csys = deepcopy(template_csys)
        ia.close()

        lambd = 299792458/freq
        self.lambd = lambd

        logger.debug(f"lambda {self.lambd}")
        logger.debug(f"freq {self.freq}")

        outname = make_unique('fft.im')
        wipe_file(outname)

        ia.open(templateim)
        ia.fft(complex=outname)
        ia.close()

        ia.open(outname)
        csys = ia.coordsys().torecord()
        cdelt = deepcopy(csys['linear0']['cdelt'])
        cdelt = np.abs(cdelt[0])

        self.ftcoords = csys
        self.ft_cdelt = [cdelt, cdelt]
        self.ft_cdelt_os = [cdelt/self.oversamp, cdelt/self.oversamp]

        logger.debug(f"self.ft_cdelt {self.ft_cdelt}")
        logger.debug(f"self.ft_cdelt_os {self.ft_cdelt_os}")
        logger.debug(f"self.oversamp {self.oversamp}")

        # Figure out the eta values (if any) from the input CSV
        # It can be a single value (like for VLA)
        # Or an (etax, etay) tuple (like for MeerKAT)
        # or non-existent.
        self.set_eta_values()

        # Number of pixels across the aperture
        npix = [int(np.round(dia*eta/lambd)) for dia, eta in zip(self.dish_dia, self.eta)]
        self.ft_npix = [int(np.floor(nn/ft_cdelt)) for nn, ft_cdelt in zip(npix, self.ft_cdelt)]
        self.ft_npix_os = [int(np.floor(nn/ft_cdelt_os)) for nn, ft_cdelt_os in zip(npix, self.ft_cdelt_os)]

        # So aperture is centred on a single pixel, to avoid offsets in the image
        #self.ft_npix_os = [nn + 1 if nn % 2 == 0 else nn for nn in self.ft_npix_os]
        #self.ft_npix = [nn + 1 if nn % 2 == 0 else nn for nn in self.ft_npix]

        # Scale it by this amount
        #if all(self.scale):
        #    self.ft_npix_os[0] *= self.scale[0]
        #    self.ft_npix_os[1] *= self.scale[1]

        #    self.ft_npix_os = [int(np.round(nn)) for nn in self.ft_npix_os]

        logger.debug(f"Scale is {self.scale}")
        logger.debug(f"npix is {npix}")
        logger.debug(f"self.ft_npix {self.ft_npix}")
        logger.debug(f"self.ft_npix_os {self.ft_npix_os}")
        logger.debug(f"self.ft_cdelt {self.ft_cdelt}")
        logger.debug(f"self.ft_cdelt_os {self.ft_cdelt_os}")

        wipe_file(outname)



    def powl(self, base: float, exp: float) -> float:
        """
        Raise base to an exponent quickly.
        Algorithm taken from https://stackoverflow.com/questions/2198138/calculating-powers-e-g-211-quickly

        Inputs:
        base        Input base, float
        exp         Exponent to raise to, float

        Returns:
        base^exp    Base raised to the exponent, float
        """
        if exp == 0:
            return 1
        elif exp == 1:
            return base
        elif (exp & 1) != 0:
            return base * self.powl(base * base, exp // 2)
        else:
            return self.powl(base * base, exp // 2)


    def gen_zernike_surface(self, coeffs: np.ndarray, x: np.ndarray, y:np.ndarray) -> np.ndarray:

        '''
        Use polynomial approximations to generate Zernike surface

        Inputs:
        coeffs      List of coefficients, np.ndarray
        x           2D array of X coordinates, from np.meshgrid
        y           2D array of Y coordinates, from np.meshgrid

        Returns:
        ZW          The Zernike surface, np.ndarray
        '''

        #Setting coefficients array
        Z = np.zeros(67)
        if len(coeffs) < len(Z):
            c = Z.copy()
            c[:len(coeffs)] += coeffs
        else:
            c = Z.copy()
            c[:len(coeffs)] += coeffs

        #Setting the equations for the Zernike polynomials
        #r = np.sqrt(powl(x,2) + powl(y,2))
        Z1  =  c[0]  * 1 # m = 0    n = 0
        Z2  =  c[1]  * x # m = -1   n = 1
        Z3  =  c[2]  * y # m = 1    n = 1
        Z4  =  c[3]  * 2*x*y # m = -2   n = 2
        Z5  =  c[4]  * (2*self.powl(x,2) + 2*self.powl(y,2) -1)# m = 0  n = 2
        Z6  =  c[5]  * (-1*self.powl(x,2) + self.powl(y,2))# m = 2  n = 2
        Z7  =  c[6]  * (-1*self.powl(x,3) + 3*x*self.powl(y,2)) #m = -3     n = 3
        Z8  =  c[7]  * (-2*x + 3*(self.powl(x,3)) + 3*x*(self.powl(y,2)))# m = -1   n = 3
        Z9  =  c[8]  * (-2*y + 3*self.powl(y,3) + 3*(self.powl(x,2)) * y)# m = 1    n = 3
        Z10 =  c[9]  * (self.powl(y,3)-3*(self.powl(x,2))*y)# m = 3 n =3
        Z11 =  c[10] * (-4*(self.powl(x,3))*y + 4*x*(self.powl(y,3)))#m = -4    n = 4
        Z12 =  c[11] * (-6*x*y + 8*(self.powl(x,3))*y + 8*x*(self.powl(y,3)))# m = -2   n = 4
        Z13 =  c[12] * (1-6*self.powl(x,2) - 6*self.powl(y,2) + 6*self.powl(x,4) + 12*(self.powl(x,2))*(self.powl(y,2)) + 6*self.powl(y,4))# m = 0  n = 4
        Z14 =  c[13] * (3*self.powl(x,2) - 3*self.powl(y,2) - 4*self.powl(x,4) + 4*self.powl(y,4))#m = 2    n = 4
        Z15 =  c[14] * (self.powl(x,4) - 6*(self.powl(x,2))*(self.powl(y,2)) + self.powl(y,4))# m = 4   n = 4
        Z16 =  c[15] * (self.powl(x,5)-10*(self.powl(x,3))*self.powl(y,2) + 5*x*(self.powl(y,4)))# m = -5   n = 5
        Z17 =  c[16] * (4*self.powl(x,3) - 12*x*(self.powl(y,2)) -5*self.powl(x,5) + 10*(self.powl(x,3))*(self.powl(y,2)) + 15*x*self.powl(y,4))# m =-3     n = 5
        Z18 =  c[17] * (3*x - 12*self.powl(x,3) - 12*x*(self.powl(y,2)) + 10*self.powl(x,5) + 20*(self.powl(x,3))*(self.powl(y,2)) + 10*x*(self.powl(y,4)))# m= -1  n = 5
        Z19 =  c[18] * (3*y - 12*self.powl(y,3) - 12*y*(self.powl(x,2)) + 10*self.powl(y,5) + 20*(self.powl(y,3))*(self.powl(x,2)) + 10*y*(self.powl(x,4)))# m = 1  n = 5
        Z20 =  c[19] * (-4*self.powl(y,3) + 12*y*(self.powl(x,2)) + 5*self.powl(y,5) - 10*(self.powl(y,3))*(self.powl(x,2)) - 15*y*self.powl(x,4))# m = 3   n = 5
        Z21 =  c[20] * (self.powl(y,5)-10*(self.powl(y,3))*self.powl(x,2) + 5*y*(self.powl(x,4)))#m = 5 n = 5
        Z22 =  c[21] * (6*(self.powl(x,5))*y - 20*(self.powl(x,3))*(self.powl(y,3)) + 6*x*(self.powl(y,5)))# m = -6 n = 6
        Z23 =  c[22] * (20*(self.powl(x,3))*y - 20*x*(self.powl(y,3)) - 24*(self.powl(x,5))*y + 24*x*(self.powl(y,5)))#m = -4   n = 6
        Z24 =  c[23] * (12*x*y + 40*(self.powl(x,3))*y - 40*x*(self.powl(y,3)) + 30*(self.powl(x,5))*y + 60*(self.powl(x,3))*(self.powl(y,3)) - 30*x*(self.powl(y,5)))#m = -2   n = 6
        Z25 =  c[24] * (-1 + 12*(self.powl(x,2)) + 12*(self.powl(y,2)) - 30*(self.powl(x,4)) - 60*(self.powl(x,2))*(self.powl(y,2)) - 30*(self.powl(y,4)) + 20*(self.powl(x,6)) + 60*(self.powl(x,4))*self.powl(y,2) + 60 *(self.powl(x,2))*(self.powl(y,4)) + 20*(self.powl(y,6)))#m = 0   n = 6
        Z26 =  c[25] * (-6*(self.powl(x,2)) + 6*(self.powl(y,2)) + 20*(self.powl(x,4)) - 20*(self.powl(y,4)) - 15*(self.powl(x,6)) - 15*(self.powl(x,4))*(self.powl(y,2)) + 15*(self.powl(x,2))*(self.powl(y,4)) + 15*(self.powl(y,6)))#m = 2   n = 6
        Z27 =  c[26] * (-5*(self.powl(x,4)) + 30*(self.powl(x,2))*(self.powl(y,2)) - 5*(self.powl(y,4)) + 6*(self.powl(x,6)) - 30*(self.powl(x,4))*self.powl(y,2) - 30*(self.powl(x,2))*(self.powl(y,4)) + 6*(self.powl(y,6)))#m = 4    n = 6
        Z28 =  c[27] * (-1*(self.powl(x,6)) + 15*(self.powl(x,4))*(self.powl(y,2)) - 15*(self.powl(x,2))*(self.powl(y,4)) + self.powl(y,6))#m = 6   n = 6
        Z29 =  c[28] * (-1*(self.powl(x,7)) + 21*(self.powl(x,5))*(self.powl(y,2)) - 35*(self.powl(x,3))*(self.powl(y,4)) + 7*x*(self.powl(y,6)))#m = -7    n = 7
        Z30 =  c[29] * (-6*(self.powl(x,5)) + 60*(self.powl(x,3))*(self.powl(y,2)) - 30*x*(self.powl(y,4)) + 7*self.powl(x,7) - 63*(self.powl(x,5))*(self.powl(y,2)) - 35*(self.powl(x,3))*(self.powl(y,4)) + 35*x*(self.powl(y,6))) #m = -5    n = 7
        Z31 =  c[30] * (-10*(self.powl(x,3)) + 30*x*(self.powl(y,2)) + 30*self.powl(x,5) - 60*(self.powl(x,3))*(self.powl(y,2)) - 90*x*(self.powl(y,4)) - 21*self.powl(x,7) + 21*(self.powl(x,5))*(self.powl(y,2)) + 105*(self.powl(x,3))*(self.powl(y,4)) + 63*x*(self.powl(y,6)))#m =-3       n = 7
        Z32 =  c[31] * (-4*x + 30*self.powl(x,3) + 30*x*(self.powl(y,2)) - 60*(self.powl(x,5)) - 120*(self.powl(x,3))*(self.powl(y,2)) - 60*x*(self.powl(y,4)) + 35*self.powl(x,7) + 105*(self.powl(x,5))*(self.powl(y,2)) + 105*(self.powl(x,3))*(self.powl(y,4)) + 35*x*(self.powl(y,6)))#m = -1  n = 7
        Z33 =  c[32] * (-4*y + 30*self.powl(y,3) + 30*y*(self.powl(x,2)) - 60*(self.powl(y,5)) - 120*(self.powl(y,3))*(self.powl(x,2)) - 60*y*(self.powl(x,4)) + 35*self.powl(y,7) + 105*(self.powl(y,5))*(self.powl(x,2)) + 105*(self.powl(y,3))*(self.powl(x,4)) + 35*y*(self.powl(x,6)))#m = 1   n = 7
        Z34 =  c[33] * (10*(self.powl(y,3)) - 30*y*(self.powl(x,2)) - 30*self.powl(y,5) + 60*(self.powl(y,3))*(self.powl(x,2)) + 90*y*(self.powl(x,4)) + 21*self.powl(y,7) - 21*(self.powl(y,5))*(self.powl(x,2)) - 105*(self.powl(y,3))*(self.powl(x,4)) - 63*y*(self.powl(x,6)))#m =3     n = 7
        Z35 =  c[34] * (-6*(self.powl(y,5)) + 60*(self.powl(y,3))*(self.powl(x,2)) - 30*y*(self.powl(x,4)) + 7*self.powl(y,7) - 63*(self.powl(y,5))*(self.powl(x,2)) - 35*(self.powl(y,3))*(self.powl(x,4)) + 35*y*(self.powl(x,6)))#m = 5  n = 7
        Z36 =  c[35] * (self.powl(y,7) - 21*(self.powl(y,5))*(self.powl(x,2)) + 35*(self.powl(y,3))*(self.powl(x,4)) - 7*y*(self.powl(x,6)))#m = 7  n = 7
        Z37 =  c[36] * (-8*(self.powl(x,7))*y + 56*(self.powl(x,5))*(self.powl(y,3)) - 56*(self.powl(x,3))*(self.powl(y,5)) + 8*x*(self.powl(y,7)))#m = -8  n = 8
        Z38 =  c[37] * (-42*(self.powl(x,5))*y + 140*(self.powl(x,3))*(self.powl(y,3)) - 42*x*(self.powl(y,5)) + 48*(self.powl(x,7))*y - 112*(self.powl(x,5))*(self.powl(y,3)) - 112*(self.powl(x,3))*(self.powl(y,5)) + 48*x*(self.powl(y,7)))#m = -6  n = 8
        Z39 =  c[38] * (-60*(self.powl(x,3))*y + 60*x*(self.powl(y,3)) + 168*(self.powl(x,5))*y -168*x*(self.powl(y,5)) - 112*(self.powl(x,7))*y - 112*(self.powl(x,5))*(self.powl(y,3)) + 112*(self.powl(x,3))*(self.powl(y,5)) + 112*x*(self.powl(y,7)))#m = -4   n = 8
        Z40 =  c[39] * (-20*x*y + 120*(self.powl(x,3))*y + 120*x*(self.powl(y,3)) - 210*(self.powl(x,5))*y - 420*(self.powl(x,3))*(self.powl(y,3)) - 210*x*(self.powl(y,5)) - 112*(self.powl(x,7))*y + 336*(self.powl(x,5))*(self.powl(y,3)) + 336*(self.powl(x,3))*(self.powl(y,5)) + 112*x*(self.powl(y,7)))#m = -2   n = 8
        Z41 =  c[40] * (1 - 20*self.powl(x,2) - 20*self.powl(y,2) + 90*self.powl(x,4) + 180*(self.powl(x,2))*(self.powl(y,2)) + 90*self.powl(y,4) - 140*self.powl(x,6) - 420*(self.powl(x,4))*(self.powl(y,2)) - 420*(self.powl(x,2))*(self.powl(y,4)) - 140*(self.powl(y,6)) + 70*self.powl(x,8) + 280*(self.powl(x,6))*(self.powl(y,2)) + 420*(self.powl(x,4))*(self.powl(y,4)) + 280*(self.powl(x,2))*(self.powl(y,6)) + 70*self.powl(y,8))#m = 0    n = 8
        Z42 =  c[41] * (10*self.powl(x,2) - 10*self.powl(y,2) - 60*self.powl(x,4) + 105*(self.powl(x,4))*(self.powl(y,2)) - 105*(self.powl(x,2))*(self.powl(y,4)) + 60*self.powl(y,4) + 105*self.powl(x,6) - 105*self.powl(y,6) - 56*self.powl(x,8) - 112*(self.powl(x,6))*(self.powl(y,2)) + 112*(self.powl(x,2))*(self.powl(y,6)) + 56*self.powl(y,8))#m = 2  n = 8
        Z43 =  c[42] * (15*self.powl(x,4) - 90*(self.powl(x,2))*(self.powl(y,2)) + 15*self.powl(y,4) - 42*self.powl(x,6) + 210*(self.powl(x,4))*(self.powl(y,2)) + 210*(self.powl(x,2))*(self.powl(y,4)) - 42*self.powl(y,6) + 28*self.powl(x,8) - 112*(self.powl(x,6))*(self.powl(y,2)) - 280*(self.powl(x,4))*(self.powl(y,4)) - 112*(self.powl(x,2))*(self.powl(y,6)) + 28*self.powl(y,8))#m = 4     n = 8
        Z44 =  c[43] * (7*self.powl(x,6) - 105*(self.powl(x,4))*(self.powl(y,2)) + 105*(self.powl(x,2))*(self.powl(y,4)) - 7*self.powl(y,6) - 8*self.powl(x,8) + 112*(self.powl(x,6))*(self.powl(y,2)) - 112*(self.powl(x,2))*(self.powl(y,6)) + 8*self.powl(y,8))#m = 6    n = 8
        Z45 =  c[44] * (self.powl(x,8) - 28*(self.powl(x,6))*(self.powl(y,2)) + 70*(self.powl(x,4))*(self.powl(y,4)) - 28*(self.powl(x,2))*(self.powl(y,6)) + self.powl(y,8))#m = 8     n = 9
        Z46 =  c[45] * (self.powl(x,9) - 36*(self.powl(x,7))*(self.powl(y,2)) + 126*(self.powl(x,5))*(self.powl(y,4)) - 84*(self.powl(x,3))*(self.powl(y,6)) + 9*x*(self.powl(y,8)))#m = -9     n = 9
        Z47 =  c[46] * (8*self.powl(x,7) - 168*(self.powl(x,5))*(self.powl(y,2)) + 280*(self.powl(x,3))*(self.powl(y,4)) - 56 *x*(self.powl(y,6)) - 9*self.powl(x,9) + 180*(self.powl(x,7))*(self.powl(y,2)) - 126*(self.powl(x,5))*(self.powl(y,4)) - 252*(self.powl(x,3))*(self.powl(y,6)) + 63*x*(self.powl(y,8)))#m = -7    n = 9
        Z48 =  c[47] * (21*self.powl(x,5) - 210*(self.powl(x,3))*(self.powl(y,2)) + 105*x*(self.powl(y,4)) - 56*self.powl(x,7) + 504*(self.powl(x,5))*(self.powl(y,2)) + 280*(self.powl(x,3))*(self.powl(y,4)) - 280*x*(self.powl(y,6)) + 36*self.powl(x,9) - 288*(self.powl(x,7))*(self.powl(y,2)) - 504*(self.powl(x,5))*(self.powl(y,4)) + 180*x*(self.powl(y,8)))#m = -5    n = 9
        Z49 =  c[48] * (20*self.powl(x,3) - 60*x*(self.powl(y,2)) - 105*self.powl(x,5) + 210*(self.powl(x,3))*(self.powl(y,2)) + 315*x*(self.powl(y,4)) + 168*self.powl(x,7) - 168*(self.powl(x,5))*(self.powl(y,2)) - 840*(self.powl(x,3))*(self.powl(y,4)) - 504*x*(self.powl(y,6)) - 84*self.powl(x,9) + 504*(self.powl(x,5))*(self.powl(y,4)) + 672*(self.powl(x,3))*(self.powl(y,6)) + 252*x*(self.powl(y,8)))#m = -3  n = 9
        Z50 =  c[49] * (5*x - 60*self.powl(x,3) - 60*x*(self.powl(y,2)) + 210*self.powl(x,5) + 420*(self.powl(x,3))*(self.powl(y,2)) + 210*x*(self.powl(y,4)) - 280*self.powl(x,7) - 840*(self.powl(x,5))*(self.powl(y,2)) - 840*(self.powl(x,3))*(self.powl(y,4)) - 280*x*(self.powl(y,6)) + 126*self.powl(x,9) + 504*(self.powl(x,7))*(self.powl(y,2)) + 756*(self.powl(x,5))*(self.powl(y,4)) + 504*(self.powl(x,3))*(self.powl(y,6)) + 126*x*(self.powl(y,8)))#m = -1   n = 9
        Z51 =  c[50] * (5*y - 60*self.powl(y,3) - 60*y*(self.powl(x,2)) + 210*self.powl(y,5) + 420*(self.powl(y,3))*(self.powl(x,2)) + 210*y*(self.powl(x,4)) - 280*self.powl(y,7) - 840*(self.powl(y,5))*(self.powl(x,2)) - 840*(self.powl(y,3))*(self.powl(x,4)) - 280*y*(self.powl(x,6)) + 126*self.powl(y,9) + 504*(self.powl(y,7))*(self.powl(x,2)) + 756*(self.powl(y,5))*(self.powl(x,4)) + 504*(self.powl(y,3))*(self.powl(x,6)) + 126*y*(self.powl(x,8)))#m = -1   n = 9
        Z52 =  c[51] * (-20*self.powl(y,3) + 60*y*(self.powl(x,2)) + 105*self.powl(y,5) - 210*(self.powl(y,3))*(self.powl(x,2)) - 315*y*(self.powl(x,4)) - 168*self.powl(y,7) + 168*(self.powl(y,5))*(self.powl(x,2)) + 840*(self.powl(y,3))*(self.powl(x,4)) + 504*y*(self.powl(x,6)) + 84*self.powl(y,9) - 504*(self.powl(y,5))*(self.powl(x,4)) - 672*(self.powl(y,3))*(self.powl(x,6)) - 252*y*(self.powl(x,8)))#m = 3  n = 9
        Z53 =  c[52] * (21*self.powl(y,5) - 210*(self.powl(y,3))*(self.powl(x,2)) + 105*y*(self.powl(x,4)) - 56*self.powl(y,7) + 504*(self.powl(y,5))*(self.powl(x,2)) + 280*(self.powl(y,3))*(self.powl(x,4)) - 280*y*(self.powl(x,6)) + 36*self.powl(y,9) - 288*(self.powl(y,7))*(self.powl(x,2)) - 504*(self.powl(y,5))*(self.powl(x,4)) + 180*y*(self.powl(x,8)))#m = 5     n = 9
        Z54 =  c[53] * (-8*self.powl(y,7) + 168*(self.powl(y,5))*(self.powl(x,2)) - 280*(self.powl(y,3))*(self.powl(x,4)) + 56 *y*(self.powl(x,6)) + 9*self.powl(y,9) - 180*(self.powl(y,7))*(self.powl(x,2)) + 126*(self.powl(y,5))*(self.powl(x,4)) - 252*(self.powl(y,3))*(self.powl(x,6)) - 63*y*(self.powl(x,8)))#m = 7     n = 9
        Z55 =  c[54] * (self.powl(y,9) - 36*(self.powl(y,7))*(self.powl(x,2)) + 126*(self.powl(y,5))*(self.powl(x,4)) - 84*(self.powl(y,3))*(self.powl(x,6)) + 9*y*(self.powl(x,8)))#m = 9       n = 9
        Z56 =  c[55] * (10*(self.powl(x,9))*y - 120*(self.powl(x,7))*(self.powl(y,3)) + 252*(self.powl(x,5))*(self.powl(y,5)) - 120*(self.powl(x,3))*(self.powl(y,7)) + 10*x*(self.powl(y,9)))#m = -10   n = 10
        Z57 =  c[56] * (72*(self.powl(x,7))*y - 504*(self.powl(x,5))*(self.powl(y,3)) + 504*(self.powl(x,3))*(self.powl(y,5)) - 72*x*(self.powl(y,7)) - 80*(self.powl(x,9))*y + 480*(self.powl(x,7))*(self.powl(y,3)) - 480*(self.powl(x,3))*(self.powl(y,7)) + 80*x*(self.powl(y,9)))#m = -8    n = 10
        Z58 =  c[57] * (270*(self.powl(x,9))*y - 360*(self.powl(x,7))*(self.powl(y,3)) - 1260*(self.powl(x,5))*(self.powl(y,5)) - 360*(self.powl(x,3))*(self.powl(y,7)) + 270*x*(self.powl(y,9)) - 432*(self.powl(x,7))*y + 1008*(self.powl(x,5))*(self.powl(y,3)) + 1008*(self.powl(x,3))*(self.powl(y,5)) - 432*x*(self.powl(y,7)) + 168*(self.powl(x,5))*y - 560*(self.powl(x,3))*(self.powl(y,3)) + 168*x*(self.powl(y,5)))#m = -6   n = 10
        Z59 =  c[58] * (140*(self.powl(x,3))*y - 140*x*(self.powl(y,3)) - 672*(self.powl(x,5))*y + 672*x*(self.powl(y,5)) + 1008*(self.powl(x,7))*y + 1008*(self.powl(x,5))*(self.powl(y,3)) - 1008*(self.powl(x,3))*(self.powl(y,5)) - 1008*x*(self.powl(y,7)) - 480*(self.powl(x,9))*y - 960*(self.powl(x,7))*(self.powl(y,3)) + 960*(self.powl(x,3))*(self.powl(y,7)) + 480 *x*(self.powl(y,9)))#m = -4   n = 10
        Z60 =  c[59] * (30*x*y - 280*(self.powl(x,3))*y - 280*x*(self.powl(y,3)) + 840*(self.powl(x,5))*y + 1680*(self.powl(x,3))*(self.powl(y,3)) +840*x*(self.powl(y,5)) - 1008*(self.powl(x,7))*y - 3024*(self.powl(x,5))*(self.powl(y,3)) - 3024*(self.powl(x,3))*(self.powl(y,5)) - 1008*x*(self.powl(y,7)) + 420*(self.powl(x,9))*y + 1680*(self.powl(x,7))*(self.powl(y,3)) + 2520*(self.powl(x,5))*(self.powl(y,5)) + 1680*(self.powl(x,3))*(self.powl(y,7)) + 420*x*(self.powl(y,9)) )#m = -2   n = 10
        Z61 =  c[60] * (-1 + 30*self.powl(x,2) + 30*self.powl(y,2) - 210*self.powl(x,4) - 420*(self.powl(x,2))*(self.powl(y,2)) - 210*self.powl(y,4) + 560*self.powl(x,6) + 1680*(self.powl(x,4))*(self.powl(y,2)) + 1680*(self.powl(x,2))*(self.powl(y,4)) + 560*self.powl(y,6) - 630*self.powl(x,8) - 2520*(self.powl(x,6))*(self.powl(y,2)) - 3780*(self.powl(x,4))*(self.powl(y,4)) - 2520*(self.powl(x,2))*(self.powl(y,6)) - 630*self.powl(y,8) + 252*self.powl(x,10) + 1260*(self.powl(x,8))*(self.powl(y,2)) + 2520*(self.powl(x,6))*(self.powl(y,4)) + 2520*(self.powl(x,4))*(self.powl(y,6)) + 1260*(self.powl(x,2))*(self.powl(y,8)) + 252*self.powl(y,10))#m = 0    n = 10
        Z62 =  c[61] * (-15*self.powl(x,2) + 15*self.powl(y,2) + 140*self.powl(x,4) - 140*self.powl(y,4) - 420*self.powl(x,6) - 420*(self.powl(x,4))*(self.powl(y,2)) + 420*(self.powl(x,2))*(self.powl(y,4)) + 420*self.powl(y,6) + 504*self.powl(x,8) + 1008*(self.powl(x,6))*(self.powl(y,2)) - 1008*(self.powl(x,2))*(self.powl(y,6)) - 504*self.powl(y,8) - 210*self.powl(x,10) - 630*(self.powl(x,8))*(self.powl(y,2)) - 420*(self.powl(x,6))*(self.powl(y,4)) + 420*(self.powl(x,4))*(self.powl(y,6)) + 630*(self.powl(x,2))*(self.powl(y,8)) + 210*self.powl(y,10))# m = 2  n = 10
        Z63 =  c[62] * (-35*self.powl(x,4) + 210*(self.powl(x,2))*(self.powl(y,2)) - 35*self.powl(y,4) + 168*self.powl(x,6) - 840*(self.powl(x,4))*(self.powl(y,2)) - 840*(self.powl(x,2))*(self.powl(y,4)) + 168*self.powl(y,6) - 252*self.powl(x,8) + 1008*(self.powl(x,6))*(self.powl(y,2)) + 2520*(self.powl(x,4))*(self.powl(y,4)) + 1008*(self.powl(x,2))*(self.powl(y,6)) - 252*(self.powl(y,8)) + 120*self.powl(x,10) - 360*(self.powl(x,8))*(self.powl(y,2)) - 1680*(self.powl(x,6))*(self.powl(y,4)) - 1680*(self.powl(x,4))*(self.powl(y,6)) - 360*(self.powl(x,2))*(self.powl(y,8)) + 120*self.powl(y,10))#m = 4     n = 10
        Z64 =  c[63] * (-28*self.powl(x,6) + 420*(self.powl(x,4))*(self.powl(y,2)) - 420*(self.powl(x,2))*(self.powl(y,4)) + 28*self.powl(y,6) + 72*self.powl(x,8) - 1008*(self.powl(x,6))*(self.powl(y,2)) + 1008*(self.powl(x,2))*(self.powl(y,6)) - 72*self.powl(y,8) - 45*self.powl(x,10) + 585*(self.powl(x,8))*(self.powl(y,2)) + 630*(self.powl(x,6))*(self.powl(y,4)) - 630*(self.powl(x,4))*(self.powl(y,6)) - 585*(self.powl(x,2))*(self.powl(y,8)) + 45*self.powl(y,10))#m = 6    n = 10
        Z65 =  c[64] * (-9*self.powl(x,8) + 252*(self.powl(x,6))*(self.powl(y,2)) - 630*(self.powl(x,4))*(self.powl(y,4)) + 252*(self.powl(x,2))*(self.powl(y,6)) - 9*self.powl(y,8) + 10*self.powl(x,10) - 270*(self.powl(x,8))*(self.powl(y,2)) + 420*(self.powl(x,6))*(self.powl(y,4)) + 420*(self.powl(x,4))*(self.powl(y,6)) - 270*(self.powl(x,2))*(self.powl(y,8)) + 10*self.powl(y,10))#m = 8    n = 10
        Z66 =  c[65] * (-1*self.powl(x,10) + 45*(self.powl(x,8))*(self.powl(y,2)) - 210*(self.powl(x,6))*(self.powl(y,4)) + 210*(self.powl(x,4))*(self.powl(y,6)) - 45*(self.powl(x,2))*(self.powl(y,8)) + self.powl(y,10))#m = 10   n = 10

        ZW =    Z1 + Z2 +  Z3+  Z4+  Z5+  Z6+  Z7+  Z8+  Z9+  Z10+ Z11+ Z12+ Z13+ Z14+ Z15+ Z16+ Z17+ Z18+ Z19+ Z20+ Z21+ Z22+ Z23+ Z24+ Z25+ Z26+ Z27+ Z28+ Z29+ Z30+ Z31+ Z32+ Z33+ Z34+ Z35+ Z36+ Z37+ Z38+ Z39+ Z40+ Z41+ Z42+ Z43+ Z44+ Z45+ Z46+ Z47+ Z48+ Z49+Z50+ Z51+ Z52+ Z53+ Z54+ Z55+ Z56+ Z57+ Z58+ Z59+Z60+ Z61+ Z62+ Z63+ Z64+ Z65+ Z66
        return ZW


    def pad_image(self, inpdat: np.ndarray) -> np.ndarray:
        """
        Pad the input image by self.padfac and return the numpy array

        Inputs:
        inpdat      Input image data, 2D array

        Returns:
        paddat      Padded image data, 2D array
        """

        shape = inpdat.shape
        logger.debug(f"Aperture shape {shape}")

        try:
            padshape = [shape[0]*self.padfac, shape[1]*self.padfac, shape[2], shape[3]]
        except IndexError:
            try:
                padshape = [shape[0]*self.padfac, shape[1]*self.padfac, shape[2], 1]
            except IndexError:
                padshape = [shape[0]*self.padfac, shape[1]*self.padfac, 1, 1]

        logger.debug(f"Padded aperture shape {shape}")

        paddat = np.zeros(padshape, dtype=complex)

        cx, cy = shape[0]//2, shape[1]//2
        pcx, pcy = padshape[0]//2, padshape[1]//2

        if shape[0] % 2:
            slicex = slice(pcx-cx, pcx+cx+1)
        else:
            slicex = slice(pcx-cx, pcx+cx)

        if shape[1] % 2:
            slicey = slice(pcy-cy, pcy+cy+1)
        else:
            slicey = slice(pcy-cy, pcy+cy)


        paddat[slicex, slicey] = inpdat

        if self.parang != 0:
            logger.info(f"Rotating to parallactic angle {self.parang}")
            paddat = self.rotate_image(paddat, self.parang)

        if self.parang_file is not None:
            logger.info(f"Calculating weighted parang average from {self.parang_file}")

            pfile_dat = np.loadtxt(self.parang_file)
            parangs = pfile_dat.T[0]
            weights = pfile_dat.T[1]

            paddat_wavg = np.zeros_like(paddat)
            for pp, ww in zip(parangs, weights):
                rotdat = self.rotate_image(paddat, pp)
                paddat_wavg += rotdat * ww
            paddat_wavg /= np.sum(ww)

            paddat = paddat_wavg

        return paddat


    def gen_jones_beams(self) -> list[str]:
        """
        Generate the Jones beams from the aperture coefficients.

        Inputs:
        None

        Returns:
        beamnames       List of filenames for Jones beams, list of str
        """

        stepsize = [(os/npix) for os, npix in zip(self.ft_cdelt_os, self.ft_npix_os)]

        logger.debug(f"ft_cdelt_os {self.ft_cdelt_os}")

        if self.ft_npix_os[0] % 2 == 0:
            x = np.linspace(-1, 1, self.ft_npix_os[0]+1)
        else:
            x = np.linspace(-1, 1, self.ft_npix_os[0])

        if self.ft_npix_os[1] % 2 == 0:
            y = np.linspace(-1, 1, self.ft_npix_os[1]+1)
        else:
            y = np.linspace(-1, 1, self.ft_npix_os[1])

        logger.debug(f"number of pixels in X (oversamp) {x.size}")
        logger.debug(f"number of pixels in Y (oversamp) {y.size}")

        xx, yy = np.meshgrid(x, y)
        maskidx = np.where(np.sqrt(xx**2 + yy**2) > 1)

        # Oversampled aperture
        jonesnames_os = ['J00_ap_os.im', 'J01_ap_os.im', 'J10_ap_os.im', 'J11_ap_os.im']
        jonesnames = ['J00_ap.im', 'J01_ap.im', 'J10_ap.im', 'J11_ap.im']

        [wipe_file(jj) for jj in jonesnames_os]
        [wipe_file(jj) for jj in jonesnames]

        jonesnames_os = [make_unique(jj) for jj in jonesnames_os]
        jonesnames = [make_unique(jj) for jj in jonesnames]

        stokes = self.df['#stokes'].unique()
        for idx, ss in enumerate(stokes):
            sdf = self.df[self.df['#stokes'] == ss]

            zreal = self.gen_zernike_surface(sdf['real'].values, xx, yy)
            zimag = self.gen_zernike_surface(sdf['imag'].values, xx, yy)

            zaperture = zreal + 1j*zimag
            zaperture[maskidx] = 0

            # Python has the array flipped relative to what CASA wants
            zaperture = np.fliplr(zaperture)

            ft_csys_os = deepcopy(self.ftcoords)
            ft_csys_os['linear0']['cdelt'] = self.ft_cdelt_os
            ft_csys_os['linear0']['crpix'] = [xx.shape[0]//2, xx.shape[1]//2]
            # The ref pixel is the input imsize//2 (ex: 4096,4096) but our
            # aperture does not span that many pixels, so adjust
            self.ftcoords['linear0']['crpix'] = [xx.shape[0]//2, xx.shape[1]//2]

            ia.fromarray(jonesnames_os[idx], zaperture[:,:,None,None], linear=True, csys=ft_csys_os)
            ia.close()

            #imregrid_coords = imregrid(jonesnames_os[idx])
            #imregrid_coords['csys'] = self.ftcoords

            ## Resample to original UV coords
            #imregrid(jonesnames_os[idx], output=jonesnames[idx], template=imregrid_coords)

        padjonesnames = ['pad_J00_ap.im', 'pad_J01_ap.im', 'pad_J10_ap.im', 'pad_J11_ap.im']
        beamnames = ['J00.im', 'J01.im', 'J10.im', 'J11.im']

        [wipe_file(jj) for jj in padjonesnames]
        [wipe_file(jj) for jj in beamnames]

        padjonesnames = [make_unique(jj) for jj in padjonesnames]
        beamnames = [make_unique(jj) for jj in beamnames]

        for jj, pp, bb in zip(jonesnames_os, padjonesnames, beamnames):
            ia.open(jj)
            dat = ia.getchunk()
            ia.close()

            paddat = self.pad_image(dat)

            ia.fromarray(pp, paddat, linear=True)
            ia.close()

            ia.open(pp)
            csys = ia.coordsys().torecord()
            padcdelt = self.ft_cdelt[0]/self.padfac
            csys['linear0']['cdelt'] = [padcdelt, padcdelt, 1, 1]
            csys['linear0']['units'] = ['lambda', 'lambda', 'Hz', '']
            csys['telescope'] = self.telescope
            ia.setcoordsys(csys)
            ia.putchunk(paddat)
            ia.close()

            ia.open(pp)
            ia.fft(complex=bb, axes=[0,1])
            ia.close()

            self.fix_jones_header(bb)

        [wipe_file(jj) for jj in jonesnames]
        [wipe_file(jj) for jj in jonesnames_os]
        [wipe_file(jj) for jj in padjonesnames]

        return beamnames


    def fix_jones_header(self, beamname : str) -> None:
        """
        Fix the header of the Jones beams generated from the apertures.

        The `ia.fft` function preserves the header information when performing a
        reverse transformation, however it does not correctly assign the RA and
        DEC axes, but rather sets them to be `1./linear` which imregrid does not
        handle correctly.

        The header is fixed in-place.

        Inputs:
        beamname    Name of the Jones beam, str

        Returns:
        None
        """

        ia.open(beamname)
        csys = ia.coordsys().torecord()
        beamdat = ia.getchunk()
        ia.close()

        wipe_file(beamname)

        ia.fromarray(beamname, beamdat)
        ia.close()

        ia.open(beamname)
        shape = ia.shape()
        beamcsys = ia.coordsys().torecord()

        beamcsys['direction0']['cdelt'][:2] = 2.0 * csys['linear0']['cdelt'][:2]
        beamcsys['direction0']['crval'] = self.template_csys['direction0']['crval']
        beamcsys['direction0']['latpole'] = self.template_csys['direction0']['latpole']
        beamcsys['direction0']['longpole'] = self.template_csys['direction0']['longpole']
        beamcsys['direction0']['units'] = ['rad', 'rad']

        # Make sure the Stokes and Spectral axes are ordered in the same way
        if 'spectral1' in self.template_csys:
            beamcsys['spectral1'] = self.template_csys['spectral1']
            beamcsys.pop('spectral2', None)
        elif 'spectral2' in self.template_csys:
            beamcsys['spectral2'] = self.template_csys['spectral2']
            beamcsys.pop('spectral1', None)

        if 'stokes1' in self.template_csys:
            beamcsys['stokes1'] = self.template_csys['stokes1']
            beamcsys.pop('stokes2', None)
        elif 'stokes2' in self.template_csys:
            beamcsys['stokes2'] = self.template_csys['stokes2']
            beamcsys.pop('stokes1', None)

        ia.setcoordsys(beamcsys)
        ia.close()


    def jones_to_mueller(self, beamnames: list[str]) -> list[str]:
        """
        Convert the input Jones beams to Mueller beams

        Inputs:
        beamnames       List of Jones beam names, list of str

        Returns:
        stokesnames     List of Stokes beam names, list of str
        """

        stokesnames = [f'I_{self.telescope}_{self.freq}MHz.im',
                    f'IQ_{self.telescope}_{self.freq}MHz.im',
                    f'IU_{self.telescope}_{self.freq}MHz.im',
                    f'IV_{self.telescope}_{self.freq}MHz.im']

        [wipe_file(ss) for ss in stokesnames]

        if self.islinear:
            # This is probably true for ASKAP as well
            if 'meerkat' in self.telescope.lower():
                Sdag_M_S = [
                    'real(IM0*CONJ(IM0) + IM1*CONJ(IM1) + IM2*CONJ(IM2) + IM3*CONJ(IM3))',
                    'real(-IM0*CONJ(IM0) - IM1*CONJ(IM1) + IM2*CONJ(IM2) + IM3*CONJ(IM3))',
                    'real(-IM0*CONJ(IM2) - IM1*CONJ(IM3) - IM2*CONJ(IM0) - IM3*CONJ(IM1))',
                    'real(1i*(IM0*CONJ(IM2) + IM1*CONJ(IM3) - IM2*CONJ(IM0) - IM3*CONJ(IM1)))'
                ]
            else: # This works for ALMA
                Sdag_M_S = [
                    'real(CONJ(IM0)*IM0 + CONJ(IM1)*IM1 + CONJ(IM2)*IM2 + CONJ(IM3)*IM3)',
                    'real(CONJ(IM0)*IM0 - CONJ(IM1)*IM1 + CONJ(IM2)*IM2 - CONJ(IM3)*IM3)',
                    'real(CONJ(IM0)*IM1 + CONJ(IM1)*IM0 + CONJ(IM2)*IM3 + CONJ(IM3)*IM2)',
                    'real(1i*(CONJ(IM0)*IM1 - CONJ(IM1)*IM0 + CONJ(IM2)*IM3 - CONJ(IM3)*IM2))'
                ]
        else:
            Sdag_M_S = [
                'real( CONJ(IM0)*IM0 + CONJ(IM1)*IM1 + CONJ(IM2)*IM2 + CONJ(IM3)*IM3 )/2.',
                'real( CONJ(IM0)*IM1 + CONJ(IM1)*IM0 + CONJ(IM2)*IM3 + CONJ(IM3)*IM2 )/2.',
                'real( 1i*(CONJ(IM0)*IM1 + CONJ(IM2)*IM3) - 1i*(CONJ(IM1)*IM0 + CONJ(IM3)*IM2) )/2.',
                'real( CONJ(IM0)*IM0 - CONJ(IM1)*IM1 + CONJ(IM2)*IM2 - CONJ(IM3)*IM3 )/2.'
            ]


        immath(imagename=beamnames, outfile=stokesnames[0], expr=Sdag_M_S[0])
        if not self.do_stokesi_only:
            immath(imagename=beamnames, outfile=stokesnames[1], expr=Sdag_M_S[1])
            immath(imagename=beamnames, outfile=stokesnames[2], expr=Sdag_M_S[2])
            immath(imagename=beamnames, outfile=stokesnames[3], expr=Sdag_M_S[3])

        ia.open(stokesnames[0])
        dat = ia.getchunk()
        maxv = np.amax(np.abs(dat))
        ia.close()

        for iss, ss in enumerate(stokesnames):
            if self.do_stokesi_only and iss > 0:
                    break

            ia.open(ss)
            dat = ia.getchunk()
            # Normalize so StokesI peaks at 1.0
            dat /= maxv
            ia.putchunk(dat)
            ia.close()

            # ia.fft() put it into linear coords, dump it back into SkyCoord
            #shutil.rmtree(ss)
            #ia.fromarray(ss, dat)
            #ia.close()

        [wipe_file(ss) for ss in beamnames]

        return stokesnames


    def _do_regrid_to_template(self, inname: str, templatecoord: dict, outcsys: dict) -> None:
        """
        Actually perform the regrid. Abstracted into it's own function for
        multi-processing purposes.

        Inputs:
        inname              Name of the input image, string
        templatecoord       Template coordinate system, dict returned from imregrid()
        outcsys             The co-ordinate system to apply to the input image prior to regrid,
                            dict returned from imregrid()

        Returns:
        None
        """

        ia.open(inname)
        dat = ia.getchunk()
        ia.setcoordsys(outcsys)
        ia.close()

        # Check if input is complex to speedup imregrid
        # Splitting into real and imag can be up to a 100X speedup within
        # imregrid
        is_complex = np.iscomplexobj(dat)

        if is_complex:
            inname_r = make_unique(f'{inname}_real.im')
            inname_i = make_unique(f'{inname}_imag.im')

            outname_r = make_unique(f'tmp_{inname}_real.im')
            outname_i = make_unique(f'tmp_{inname}_imag.im')

            wipe_file(outname_r)
            wipe_file(outname_i)

            immath(inname, outfile=inname_r, expr='REAL(IM0)')
            immath(inname, outfile=inname_i, expr='IMAG(IM0)')

            imregrid(inname_r, template=templatecoord, output=outname_r, overwrite=True, axes=[0,1], interpolation='cubic', decimate=10)
            imregrid(inname_i, template=templatecoord, output=outname_i, overwrite=True, axes=[0,1], interpolation='cubic', decimate=10)

            outname = f'tmp_{inname}.im'
            outname = make_unique(outname)
            wipe_file(outname)
            immath([outname_r, outname_i], outfile=outname, expr='complex(IM0, IM1)')

            wipe_file(outname_r)
            wipe_file(outname_i)
            wipe_file(inname_r)
            wipe_file(inname_i)

        else:
            outname = f'tmp_{inname}.im'
            outname = make_unique(outname)
            wipe_file(outname)

            imregrid(inname, template=templatecoord, output=outname, overwrite=True, axes=[0,1], interpolation='cubic', decimate=10)

        shutil.rmtree(inname)
        shutil.move(outname, inname)


    def regrid_to_template(self, stokes_beams: list[str], templateim: str) -> None:
        """
        Match co-ordinates between the template image and the PBs, and regrid
        the PBs to lie on the same co-ordinate grid.

        Inputs:
        stokes_beams    List of names of the Stokes beam models, list of str
        templateim      Name of the template image, str

        Returns:
        None
        """

        ia.open(stokes_beams[0])
        outcsys = ia.coordsys().torecord()
        ia.close()

        templatecoord = imregrid(templateim)
        # We are only generating one stokes plane at a time, so force stokes planes = 1
        templatecoord['shap'][-1] = 1

        _do_regrid_partial = partial(self._do_regrid_to_template, templatecoord=templatecoord, outcsys=outcsys)

        if self.parallel and not self.do_stokesi_only:
            pool = multiprocessing.Pool(NCPU)
            pool.map(_do_regrid_partial, stokes_beams)
        else:
            for iss, ss in enumerate(stokes_beams):
                if self.do_stokesi_only and iss > 0:
                        break

                _do_regrid_partial(ss)


    def _do_regrid_pointing_offset(self, inname: str, templatecoord: dict, offset : list[float]) -> None:
        """
        Actually perform the regrid for the pointing offset. Abstracted into
        it's own function for multi-processing purposes.

        Inputs:
        inname              Name of the input image, string
        templatecoord       Template coordinate system, dict returned from imregrid(). This
                            must have the reference pixel as the centre of the image.
        offset              The value of the offset to correct for, in pixels. list of float

        Returns:
        None
        """

        outname = f'tmp_{inname}.im'
        outname = make_unique(outname)
        wipe_file(outname)

        # Move the crpix to the offset so that it gets regridded back to the centre
        ia.open(inname)
        csys = ia.coordsys().torecord()
        crpix = csys['direction0']['crpix']
        csys['direction0']['crpix'] = [crpix[0] + offset[0], crpix[1] + offset[1]]
        ia.setcoordsys(csys)
        ia.close()

        imregrid(inname, template=templatecoord, output=outname, overwrite=True, axes=[0,1], interpolation='cubic', decimate=10)

        shutil.rmtree(inname)
        shutil.move(outname, inname)


    def fix_pointing_offset(self, stokes_beams : list[str]) -> None:
        """

        The Stokes I beam can be produced with a few pixel offset from the
        centre of the image. Force it such that the peak is centred at the
        central pixel of the image. The other Stokes images will be regridded by
        an identical amount.

        Inputs:
        stokes_beams    List of names of the Stokes beam models, list of str

        Returns:
        None
        """

        csys_template = imregrid(stokes_beams[0])
        shape = csys_template['shap']
        cx, cy = csys_template['csys']['direction0']['crpix']

        ia.open(stokes_beams[0])
        dat = ia.getchunk()
        ia.close()

        maxpos = np.unravel_index(np.argmax(dat), dat.shape)

        offset = [maxpos[0] - cx, maxpos[1] - cy]

        _do_regrid_partial = partial(self._do_regrid_pointing_offset, templatecoord=csys_template, offset=offset)

        if self.parallel and not self.do_stokesi_only:
            pool = multiprocessing.Pool(NCPU)
            pool.map(_do_regrid_partial, stokes_beams)
        else:
            for iss, ss in enumerate(stokes_beams):
                if self.do_stokesi_only and iss > 0:
                        break

                _do_regrid_partial(ss)



    def fix_imhead(self, stokes_beams: list[str]) -> None:
        """
        Put in the correct metadata into the final images. They contain the default telescope name (ALMA)
        and the frequency of the template image, not of the generated PB.

        Inputs:
        stokes_beams    List of names of the Stokes beam models, list of str

        Returns:
        None
        """

        for idx, stok_beam in enumerate(stokes_beams):
            if self.do_stokesi_only and idx > 0:
                break

            ia.open(stok_beam)
            csys = ia.coordsys()

            init_tel = csys.telescope()
            csys.settelescope(self.telescope)
            final_tel = csys.telescope()

            csys_dict = csys.torecord()
            ia.setcoordsys(csys_dict)

            logger.debug(f"Telescope changed from {init_tel} to {final_tel}")

            # Below is probably not needed, so it's commented out.
            # However in case it is, it's here.
            #spectral_crval = csys.referencevalue('spectral')['numeric'][0]
            #spectral_crpix = csys.referencepixel('spectral')['numeric'][0]
            #spectral_cdelt = csys.increment('spectral')['numeric'][0]

            #if crpix != 1:
            #    if crpix < 0:
            #        refval = crval - cdelt*crpix
            #    else:
            #        refval = crval + cdelt*crpix
            #else:
            #    refval = crval
            #
            #csys.setreferencevalue(refval)

            #logger.debug(f"Reference frequency is {refval} Hz")
            #logger.debug(f"self.freq is {self.freq} MHz")

            csys.done()
            ia.close()



    def rotate_image(self, indat: np.ndarray, parang: float) -> np.ndarray:
        """
        Rotate the input data by parang degrees in place.

        Inputs:
        indat           Input image data, 2D array
        parang          Parallactic angle in deg, float

        Returns:
        rotdat          Rotated image, 2D array
        """

        if indat.dtype == complex:
            indat_r = rotate(indat.real, angle=parang, reshape=False)
            indat_i = rotate(indat.imag, angle=parang, reshape=False)

            rotdat = indat_r + 1j*indat_i
        else:
            rotdat = rotate(indat, angle=parang, reshape=False)

        return rotdat
