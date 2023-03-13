#!/usr/bin/env python3

"""
Module that holds functionality related to telescopes.  The class and functions
to return the metadata relating to the different instruments, and the class and
functions to setup the aperture illumination, and obtain the model PB.
"""

import numpy as np
import numpy.typing as npt
from typing import Union

import logging
logger = logging.getLogger(__name__)
logger.setLevel('INFO')

class TelescopeInfo():
    """
    Given the telescope name, return various parameters related to the
    instrument.
    """

    def __init__(self):
        self._telescope_name = ''
        self.known_telescopes = ['MeerKAT', 'VLA', 'JVLA', 'EVLA', 'ALMA', 'uGMRT', 'GMRT']
        self.is_known = True
        self.is_linear = False
        self.dish_diameter = [-1, -1]
        self.band = ''


    @property
    def telescope_name(self) -> str:
        return self._telescope_name

    @telescope_name.setter
    def telescope_name(self, name:str) :
        # If telescope name is unknown, print warning and set flag
        if not any([name.lower() in tname.lower() for tname in self.known_telescopes]):
            self.is_known = False
            logger.warning(f'Telescope name {name} is unknown. Using default parameters everywhere. '
                           'Please pass in specific values through the command line if necessary.')

        # XXX Set the telescope name to a standard string here, based on the input
        # XXX Makes parsing it again much easier.
        self._telescope_name = name


    def get_feed_basis(self) -> None:
        """
        Get the feed basis from the telescope name.
        """

        # XXX : Need to take into acount uGMRT linear feeds at L-Band and VLA
        # linear feeds at P Band

        if self.telescope_name == '':
            raise ValueError('Please set the telescope name to determine the band.')

        if not self.is_known:
            logger.warning(f"Using default linear basis for unknown telescope {self.telescope_name}")
            self.is_linear = False
            return

        if "vla" in self.telescope_name.lower():
            self.is_linear = False
        elif "meerkat" in self.telescope_name.lower():
            self.is_linear = True
        elif "alma" in self.telescope_name.lower():
            self.is_linear = True
        elif "gmrt" in self.telescope_name.lower():
            self.is_linear = False
        else: # This should never happen
            raise ValueError("Unable to determine feed basis, something went wrong.")

        basis = "linear" if self.is_linear else "circular"
        logger.info(f"Using {basis} basis for telescope {self.telescope_name}")
        return basis


    def get_dish_diameter(self) -> None:
        """
        Get the dish diameter from the telescope name.
        """

        if self.telescope_name == '':
            raise ValueError('Please set the telescope name to determine the band.')

        if not self.is_known:
            logger.warning(f"Using default dish diameter of 1m for unknown telescope {self.telescope_name}")
            self.dish_diameter = [1,1]
            return

        #if self.dish_diameter is not None:
        #    logger.info(f"Overriding automatic determination of telescope dish "
        #                "diameter. Using the input of {self.dish_diameter} m")
        #    return

        if 'vla' in self.telescope_name.lower():
            self.dish_diameter = [25, 25]
        elif 'meerkat' in self.telescope_name.lower():
            self.dish_diameter = [13.5, 13.5]
        elif 'alma' in self.telescope_name.lower():
            self.dish_diameter = [12, 12]
        elif 'gmrt' in self.telescope_name.lower():
            self.dish_diameter = [45, 45]
        else: # This should never happen.
            raise ValueError("Unable to determine telescope type. Something went wrong.")

        return self.dish_diameter
        logger.info(f"Using dish diameter of ({self.dish_diameter[0]}m, {self.dish_diameter[1]}m) for telescope {self.telescope_name}.")


    def get_band(self, inp_frequency: Union[float,npt.NDArray[float]]) -> None:
        """
        Given the input frequency/frequencies, determine the band of the telescope.

        This is necessary, because the PB models can be different for different bands, even at overlapping frequencies.
        This sets the TelescopeInfo().band attribute.

        Inputs:
        inp_frequency           The input frequency, astropy.quantity.Quantity

        Returns:
        None
        """

        if self.telescope_name == '':
            raise ValueError('Please set the telescope name to determine the band.')

        if self.telescope_name == 'MeerKAT':
            pass
