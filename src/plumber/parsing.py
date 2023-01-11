#!/usr/bin/env python3

"""
Module to parse the input images and CSV files.
"""

import os
import numpy as np
import pandas as pd
from typing import Union

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
    def dataframe_list(self):
        return self._dataframe_list

    @property
    def frequency_list(self):
        return self._frequency_list

    @property
    def nstokes_csv(self):
        return self._nstokes_csv


    def get_zcoeffs(self, csv: str, imfreq: list[float]) -> Union[pd.Dataframe, float, int]:
    """
    Given the input frequency of the image, returns the Pandas dataframe with
    the coefficients from the CSV file corresponding to the input frequency. The
    frequencies in the input CSV file are expected to be in MHz.

    `imfreq` should be a list of floats, that are `astropy.units` quantities.

    Inputs:
    csv                 Input CSV filename, string
    imfreq              Image frequency in MHz, list[float]

    Returns:
    dataframe_list      List of dataframe containing the subset of the CSV file which is
                        closest to imfreq.
    frequency_list      List of frequencies that were found in the CSV file,
                        closest matching the input frequencies.
    nstokes_csv         The number of Stokes parameters in the input CSV file.
    """

    self.dataframe_list = []
    self.frequency_list = []

    df = pd.read_csv(csv, skipinitialspace=True)
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
