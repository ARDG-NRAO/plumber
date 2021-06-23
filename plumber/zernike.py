#!/usr/bin/env python3

"""
Module to hold all the function that deal with the Zernike generation or CSV
parsing
"""

import shutil
import pandas as pd
import numpy as np

from scipy.ndimage import rotate

from copy import deepcopy

import multiprocessing
from functools import partial, wraps

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)-15s %(levelname)s: %(message)s",
                        level=logging.INFO)

from plumber.misc import wipe_file

from casatasks import immath, imregrid
from casatools import image
ia = image()


def background(f):
    @wraps(f)
    async def wrapped(*args, **kwargs):
        loop = asyncio.get_event_loop()
        # Parliament funkadelic sends their regards
        p_func = partial(f, *args, **kwargs)
        return await loop.run_in_executor(ProcessPoolExecutor(), p_func)

    return wrapped


def get_zcoeffs(csv, imfreq):
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

    df = pd.read_csv(csv, skipinitialspace=True)
    nstokes = df['#stokes'].unique().size
    freqs = df['freq'].unique()

    idx = np.argmin(np.abs(freqs - imfreq))
    zfreq = freqs[idx]

    zdf = df[df['freq'] == zfreq]

    return zdf, zfreq, nstokes



class zernikeBeam():
    """
    Generate the model beam on a 2D pixel grid.
    """

    def __init__(self):
        self.freq = None
        self.padfac = None
        self.dish_dia = None
        self.islinear = None
        self.parallel = False


    def initialize(self, df, templateim, padfac=8, dish_dia=None, islinear=None, stokesi=False, parang=None, parallel=False):
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

        self.dish_dia = dish_dia
        self.get_dish_diameter(templateim)
        self.islinear = islinear
        self.get_feed_basis()

        self.npix, self.cdelt_aperture = self.get_npix_aperture(templateim)

        self.padfac = int(padfac)
        self.parallel = parallel
        self.parang = parang

        # In MHz
        self.freq = df['freq'].unique()[0]
        self.do_stokesi_only = stokesi


    def get_feed_basis(self):
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
        elif "gmrt" in self.telescope.lower() and self.freq < 1e3:
            self.islinear = False
        elif "gmrt" in self.telescope.lower() and self.freq >= 1e3:
            self.islinear = True
        else:
            raise ValueError("Unable to determine feed bsis, unknown telescope. "
            "Please pass in the feed basis via the islinear paramter to "
            ".initialize()")



    def get_dish_diameter(self, templateim):
        """
        Get the dish diameter for different instruments.

        Input:
        templateim      Name of the input template image, string

        Returns:
        None, sets self.dish_dia in place
        """

        if self.dish_dia is not None:
            logger.info(f"Overriding automatic determination of telescope dish "
                        "diameter. Using the input of {self.dish_dia} m")
            return

        ia.open(templateim)
        csys = ia.coordsys().torecord()
        ia.close()

        self.telescope = csys['telescope']

        if 'vla' in self.telescope.lower():
            self.dish_dia = 25
        elif 'meerkat' in self.telescope.lower():
            self.dish_dia = 13.5
        elif 'alma' in self.telescope.lower():
            self.dish_dia = 12
        elif 'gmrt' in self.telescope.lower():
            self.dish_dia = 45
        else:
            raise ValueError('Unknown telescope type. Please initialize the '
                             'class with the dish_diameter value in metres.')


    def get_npix_aperture(self, templateim):
        """
        Calculate the number of pixels across the aperture from the template
        image.
        """

        ia.open(templateim)
        template_csys = ia.coordsys().torecord()
        ia.close()

        try:
            unit = template_csys['spectral2']['unit']
        except KeyError:
            unit = template_csys['spectral1']['unit']

        if unit.lower() == 'hz':
            fac = 1
        elif unit.lower() == 'mhz':
            fac = 1e-6
        elif unit.lower() == 'ghz':
            fac = 1e-9
        else:
            raise ValueError("Unknown unit in image.")

        try:
            freq = template_csys['spectral2']['restfreq']
        except KeyError:
            freq = template_csys['spectral1']['restfreq']

        lambd = 299792458/(freq*fac)

        outname = 'fft.im'
        wipe_file(outname)
        ia.open(templateim)
        ia.fft(complex=outname)
        ia.close()

        ia.open(outname)
        csys = ia.coordsys().torecord()
        cdelt = csys['linear0']['cdelt']
        cdelt = np.abs(cdelt[0])

        # Lambda scaling - accounts for changes in effective illumination across the band
        try:
            eta = self.df['eta'].unique()[0]
        except KeyError:
            logger.info("Eta column does not exist in the input CSV file, not "
                        "applying frequency dependent scaling.")
            eta = 1.

        # Number of pixels across a 25m array
        npix = self.dish_dia*eta/lambd
        npix = int(np.ceil(npix/cdelt))

        if npix < 10:
            npix *= 2
            cdelt /= 2.

        # So aperture is centred on a single pixel, to avoid offsets in the image
        if npix % 2 == 0:
            npix += 1

        wipe_file(outname)

        return int(npix), cdelt



    def powl(self, base, exp):
        """
        Algorithm taken from https://stackoverflow.com/questions/2198138/calculating-powers-e-g-211-quickly
        """
        if exp == 0:
            return 1
        elif exp == 1:
            return base
        elif (exp & 1) != 0:
            return base * self.powl(base * base, exp // 2)
        else:
            return self.powl(base * base, exp // 2)


    def gen_zernike_surface(self, coeffs, x, y):

        '''
        Use polynomial approximations to generate Zernike surface
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


    def pad_image(self, inpdat):
        """
        Pad the input image by self.padfac and return the numpy array
        """

        shape = inpdat.shape
        try:
            padshape = [shape[0]*self.padfac, shape[1]*self.padfac, shape[2], shape[3]]
        except IndexError:
            try:
                padshape = [shape[0]*self.padfac, shape[1]*self.padfac, shape[2], 1]
            except IndexError:
                padshape = [shape[0]*self.padfac, shape[1]*self.padfac, 1, 1]

        paddat = np.zeros(padshape, dtype=complex)

        cx, cy = shape[0]//2, shape[1]//2
        pcx, pcy = padshape[0]//2, padshape[1]//2

        if shape[0] % 2:
            paddat[pcx-cx:pcx+cx+1, pcy-cy:pcy+cy+1] = inpdat
        else:
            paddat[pcx-cx:pcx+cx, pcy-cy:pcy+cy] = inpdat

        if self.parang != None:
            paddat_r = rotate(paddat.real, angle=self.parang, reshape=False)
            paddat_i = rotate(paddat.imag, angle=self.parang, reshape=False)

            paddat = paddat_r + 1j*paddat_i

        return paddat


    def gen_jones_beams(self):
        """
        Generate the Jones beams from the aperture coefficients.

        Inputs:
        None

        Returns:
        beamnames       List of filenames for Jones beams, list of str
        """
        x = np.linspace(-1, 1, self.npix)
        y = np.linspace(-1, 1, self.npix)

        xx, yy = np.meshgrid(x, y)
        maskidx = np.where(np.sqrt(xx**2 + yy**2) > 1)

        jonesnames = ['J00_ap.im', 'J01_ap.im', 'J10_ap.im', 'J11_ap.im']

        [wipe_file(jj) for jj in jonesnames]

        stokes = self.df['#stokes'].unique()
        for idx, ss in enumerate(stokes):
            sdf = self.df[self.df['#stokes'] == ss]

            zreal = self.gen_zernike_surface(sdf['real'].values, xx, yy)
            zimag = self.gen_zernike_surface(sdf['imag'].values, xx, yy)


            zaperture = zreal + 1j*zimag
            zaperture[maskidx] = 0

            ia.fromarray(jonesnames[idx], zaperture[:,:,None,None], linear=True)
            ia.close()

        padjonesnames = ['pad_J00_ap.im', 'pad_J01_ap.im', 'pad_J10_ap.im', 'pad_J11_ap.im']
        beamnames = ['J00.im', 'J01.im', 'J10.im', 'J11.im']
        [wipe_file(jj) for jj in padjonesnames]
        [wipe_file(jj) for jj in beamnames]

        for jj, pp, bb in zip(jonesnames, padjonesnames, beamnames):
            ia.open(jj)
            dat = ia.getchunk()
            ia.close()

            paddat = self.pad_image(dat)

            ia.fromarray(pp, paddat, linear=True)
            ia.close()

            ia.open(pp)
            csys = ia.coordsys().torecord()
            padcdelt = self.cdelt_aperture/self.padfac
            csys['linear0']['cdelt'] = [padcdelt, padcdelt, 1, 1]
            csys['linear0']['units'] = ['lambda', 'lambda', 'Hz', '']
            csys['telescope'] = self.telescope
            ia.setcoordsys(csys)
            ia.close()


            ia.open(pp)
            ia.fft(complex=bb, axes=[0,1])
            ia.close()

        [wipe_file(jj) for jj in jonesnames]
        [wipe_file(jj) for jj in padjonesnames]

        return beamnames


    def jones_to_mueller(self, beamnames):
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
            Sdag_M_S = [
                'real( CONJ(IM0)*IM0 + CONJ(IM1)*IM1 + CONJ(IM2)*IM2 + CONJ(IM3)*IM3 )',
                'real( CONJ(IM0)*IM0 - CONJ(IM1)*IM1 + CONJ(IM2)*IM2 - CONJ(IM3)*IM3 )',
                'real( CONJ(IM0)*IM1 + CONJ(IM1)*IM0 + CONJ(IM2)*IM3 + CONJ(IM3)*IM2 )',
                'real( 1i*(CONJ(IM0)*IM1 - CONJ(IM1)*IM0 + CONJ(IM2)*IM3 - CONJ(IM3)*IM2) )'
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
            ia.close()

            # ia.fft() put it into lienar coords, dump it back into SkyCoord
            shutil.rmtree(ss)
            ia.fromarray(ss, dat)
            ia.close()

        [wipe_file(ss) for ss in beamnames]

        return stokesnames

    def _do_regrid(self, inname, templatecoord, outcsys):
        """
        Actually perform the regird. Abstracted into it's own function for
        multi-processing purposes.

        Inputs:
        templatecoord       Template coordinate system, dict returned from imregrid()
        inname              Name of the input image, string

        Returns:
        None
        """

        ia.open(inname)
        ia.setcoordsys(outcsys)
        ia.close()

        outname = f'tmp_{inname}.im'
        wipe_file(outname)

        imregrid(inname, template=templatecoord, output=outname, overwrite=True, axes=[0,1], interpolation='cubic', decimate=10)

        shutil.rmtree(inname)
        shutil.move(outname, inname)


    def regrid_to_template(self, stokes_beams, templateim):
        """
        Match co-ordinates between the template image and the PBs, and regrid
        the PBs to lie on the same co-ordinate grid.

        Inputs:
        stokes_beams    List of names of the Stokes beam models, list of str
        templateim      Name of the template image, str

        Returns:
        None
        """

        ia.open(templateim)
        templatecsys = ia.coordsys().torecord()
        ia.close()

        ia.open(stokes_beams[0])
        dat = ia.getchunk()
        shape = ia.shape()
        ia.close()

        imextent = 1./(self.cdelt_aperture)
        imcdelt = (imextent/shape[0])

        templatecoord = imregrid(templateim)

        crpix = np.unravel_index(np.argmax(dat), shape)
        outcsys = deepcopy(templatecsys)
        outcsys['direction0']['crpix'] = np.asarray(crpix).astype(int)
        outcsys['direction0']['cdelt'] = [-imcdelt, imcdelt]
        outcsys['direction0']['units'] = ['rad', 'rad']
        outcsys['direction0']['latpole'] = templatecoord['csys']['direction0']['latpole']
        outcsys['direction0']['longpole'] = templatecoord['csys']['direction0']['longpole']

        _do_regrid_partial = partial(self._do_regrid, templatecoord=templatecoord, outcsys=outcsys)

        if self.parallel and not self.do_stokesi_only:
            pool = multiprocessing.Pool(4)
            pool.map(_do_regrid_partial, stokes_beams)
        else:
            for iss, ss in enumerate(stokes_beams):
                if self.do_stokesi_only and iss > 0:
                        break

                self._do_regrid(ss, templatecoord, outcsys)


    def _do_rotate(self, imname, parang):
        """
        Rotate the input image by parang degrees in place.

        Inputs:
        imname          Input image name, string
        parang          Parallactic angle in deg, float

        Returns:
        None
        """

        tmpname = 'tmp_' + imname
        ia.open(imname)
        ia.rotate(outfile=tmpname, pa=f'{parang}deg')
        ia.close()

        shutil.move(tmpname, imname)

    def rotate_beam(self, stokes_beams, parang):
        """
        Rotate the beam by parang degrees

        Inputs:
        stokes_beams    Names of the input beams, array of strings
        parang          Parallactic angle in deg, float

        Returns:
        stokes_beams    Rotated beams (in place), array of strings
        """

        if self.parallel and not self.do_stokesi_only:
            _do_rotate_partial = partian(self._do_rotate, parang=parang)
            pool = multiprocessing.Pool(4)
            pool.map(_do_rotate_partial, stokes_beams)
        else:
            for iss, ss in enumerate(stokes_beams):
                if self.do_stokesi_only and iss > 0:
                        break

                self._do_rotate(ss, parang)
