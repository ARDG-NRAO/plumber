#!/usr/bin/env python3

"""
Module to parse the input images and CSV files.
"""

import os
import numpy as np
import numpy.typing as npt
import pandas as pd

import astropy.units as u

import logging
logger = logging.getLogger(__name__)
logger.setLevel('INFO')

from typing import Union, List, TypeVar
Quantity = TypeVar('astropy.units.quantity.Quantity')



class FileParser:
    """
    A base class to parse files of different kinds. Doesn't do much on it's own.
    """

    def __init__(self):
        self._filepath = None

    @property
    def filepath(self) -> str:
        return self._filepath


    @filepath.setter
    def filepath(self, filepath:str) -> None:
        """
        Setter function for the filepath attribute.

        Checks the existence of the filepath before setting it.
        """

        if not os.path.exists(filepath):
            err = f"File {filepath} does not exist."
            raise FileNotFoundError(err)

        self._filepath = filepath



class CSVParser(FileParser):
    """
    Parse the input CSV file(s) and return a list of pandas DataFrames.
    """

    def __init__(self):
        self._dataframe_list = None
        self._frequency_list = None
        self._nstokes_csv = None

    @property
    def dataframe_list(self) -> List[pd.DataFrame]:
        return self._dataframe_list

    @dataframe_list.setter
    def dataframe_list(self, value: List[pd.DataFrame]) -> None:
        self._dataframe_list = value

    @property
    def frequency_list(self):
        return self._frequency_list

    @frequency_list.setter
    def frequency_list(self, value: Quantity) -> None:
        self._frequency_list = value

    @property
    def nstokes_csv(self):
        return self._nstokes_csv

    @nstokes_csv.setter
    def nstokes_csv(self, value: float) -> None:
        self._nstokes_csv = value


    def csv_to_df(self, csv: str) -> pd.DataFrame:
        """
        Convert the input CSV into a Pandas DataFrame.

        Inputs:
        csv             Input CSV filename, str

        Returns:
        df              Pandas dataframe
        """

        self.filepath = csv
        df = pd.read_csv(self.filepath, skipinitialspace=True)

        return df


    def get_zcoeffs(self, df: pd.DataFrame, imfreq: Quantity) -> Union[pd.DataFrame, float, int]:
        """
        Given the input frequency of the image, returns the Pandas dataframe with
        the coefficients from the input dataframe corresponding to the input
        frequency. The frequencies in the input dataframe file are expected to be in
        MHz.

        `imfreq` should be a list of floats, that are `astropy.units` quantities.

        Inputs:
        df                  Input pandas DataFrame
        imfreq              Image frequency in MHz, astropy.units.quantity.Quantity

        Returns:
        dataframe_list      List of dataframe containing the subset of the CSV file which is
                            closest to imfreq.
        frequency_list      List of frequencies that were found in the CSV file,
                            closest matching the input frequencies.
        nstokes_csv         The number of Stokes parameters in the input CSV file.
        """

        self.dataframe_list = []
        self.frequency_list = []

        if not type(imfreq) is u.quantity.Quantity:
            raise TypeError(f'Input {imfreq} must be an astropy.units Quantity.')

        self.nstokes_csv = df['#stokes'].unique().size
        freqs = df['freq'].unique()

        # imfreq is astropy.units, so explicitly convert into MHz before compare.
        for ifreq in imfreq:
            idx = np.argmin(np.abs(freqs - ifreq.to(u.MHz).value))
            zfreq = freqs[idx]

            zdf = df[df['freq'] == zfreq]

            self.dataframe_list.append(zdf)
            self.frequency_list.append(zfreq)

        return self.dataframe_list, self.frequency_list, self.nstokes_csv


class ImageParser(FileParser):
    """
    Parse the input FITS/CASA image, and return the data and metadata.
    """

    def __init__(self):
        self.imsize = [0,0]
        self.imfreqs = [-1,]
        self.is_stokes_cube = False

    @staticmethod
    def get_telescope(templateim: str) -> str:
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
        logger.info(f"Telescope is {telescope}.")

        return csys['telescope']


    def parse_image(self, imagename : str) -> Union[npt.NDArray[np.float], np.NDArray[np.float], bool]:
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
            self.is_stokes_cube = True
        else:
            self.is_stokes_cube = False

        shape = cube.shape
        self.imsize = [shape[0], shape[1]]

        # This is always a list
        self.imfreqs = cube.stokes_data['I'].spectral_axis

        return self.imsize, self.imfreqs, self.is_stokes_cube
