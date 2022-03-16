#!/usr/bin/env python3

import os
import logging

#from typing import Any
#from _typeshed import Self

import numpy as np

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord

from astroplan import Observer
from astroplan.constraints import observability_table

from typing import Union

from casatools import msmetadata, measures, quanta
msmd = msmetadata()
qa = quanta()
me = measures()

# Add colours for warnings and errors
logging.addLevelName(
    logging.WARNING, "\033[1;31m%s\033[1;0m" % logging.getLevelName(logging.WARNING))
logging.addLevelName(
    logging.ERROR, "\033[1;41m%s\033[1;0m" % logging.getLevelName(logging.ERROR))

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)-20s %(levelname)-8s %(message)s',
    handlers=[
        logging.FileHandler("plumber.log"),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger()


class ParallacticAngle():
    def __init__(self, ms: str=None, use_astropy: bool=False) -> None:
        super().__init__()

        self.parangs = []
        self.parang_slope = []

        self.ha = []

        self.telescope_name = "VLA"
        self.telescope_pos = []

        self.source_pos = []
        self.source_times = []

        self.field_names = []
        self.target_fields = []

        self.use_astropy = use_astropy

        if ms is not None:
            self.fill_ms_attributes(ms)


    def fill_ms_attributes(self, measurement_set: str, field: str=None) -> None:
        msmd.open(measurement_set)
        self.telescope_name = msmd.observatorynames()

        if len(self.telescope_name) > 1:
            logger.warning(f"Multiple observatory names found : {','.join(self.telescope_name)}. "
                           f"Using the first one, {self.telesocpe_name[0]}")

        self.telescope_name = self.telescope_name[0]
        self.telescope_pos = msmd.observatoryposition()

        # Convert from rad, rad, m to deg, deg, m
        self.telescope_pos = [
                    self.telescope_pos['m0']['value'],
                    self.telescope_pos['m1']['value'],
                    self.telescope_pos['m2']['value']
        ]

        logger.info(f"Telescope {self.telescope_name} is at co-ordinates "
                f"{np.rad2deg(self.telescope_pos[0]):.2f} deg, "
                f"{np.rad2deg(self.telescope_pos[1]):.2f} deg, "
                f"{self.telescope_pos[2]:.2f} m")

        self.field_names = msmd.fieldnames()
        self.target_fields = msmd.fieldsforintent('TARGET', asnames=True)

        if len(self.field_names) == 0:
            raise ValueError(
                "No fields found in input MS. Please check your inputs."
            )

        if field is None:
            if len(self.target_fields) == 0:
                logger.warning(
                    f'No fields found with intent TARGET. Using field 0 {self.field_names[0]}'
                )

                logger.warning(
                    'If this is not desired, pass in --field on the command line.'
                )

                field = self.field_names[0]

            elif len(self.target_fields) > 1:
                logger.warning(
                    f'Multiple fields found with intent TARGET. Using the first one {self.target_fields[0]}'
                )
                field = self.target_fields[0]

            elif len(self.target_fields) == 1:
                logger.info(f'Using field {self.target_fields[0]}')
                field = self.target_fields[0]
        else:
            if field not in self.field_names:
                outmsg = f"Input field {field} not found in MS. Available fields are {', '.join(self.field_names)}"
                raise ValueError(outmsg)

        # CASA likes to store the time as mjd seconds
        self.source_times = msmd.timesforfield(msmd.fieldsforname(field)[0])
        self.source_times = Time(self.source_times/(3600.* 24), format='mjd')

        self.source_pos = msmd.phasecenter(msmd.fieldsforname(field)[0])
        self.source_pos = SkyCoord(
                            self.source_pos['m0']['value'],
                            self.source_pos['m1']['value'], unit='radian'
                          )

        logger.info(f"Co-ordinates of source are {self.source_pos.to_string('hmsdms')}")

        if self.use_astropy:
            self.observer = Observer(longitude=self.telescope_pos[0]*u.rad, latitude=self.telescope_pos[1]*u.rad,
                    elevation=self.telescope_pos[2]*u.m, name=self.telescope_name[0], timezone='Etc/GMT0')

            self.parangs = np.rad2deg(self.observer.parallactic_angle(self.source_times, self.source_pos))
            self.parangs = self.parangs.value

        else:
            self.ha = np.asarray([self.hour_angle(self.source_pos.ra.hour, tt.mjd*24*3600.,  observatory=self.telescope_name) for tt in self.source_times])
            #self.ha *= u.rad

            # TODO: This can be very slow with a large number of timestamps/HAs to iterate over.
            # Store only the parangs, skip the slope
            self.parangs = np.asarray([self.parallactic_angle(hh, self.telescope_pos[1], self.source_pos.dec.rad)[0] for hh in self.ha])
            self.parangs = np.rad2deg(self.parangs)

            # Store only the slope
            self.parang_slope = np.asarray([self.parallactic_angle(hh, self.telescope_pos[1], self.source_pos.dec.rad)[1] for hh in self.ha])

        self.parang_start = self.parangs[0]
        self.parang_mid = self.parangs[self.parangs.size//2]
        self.parang_end = self.parangs[-1]

        logger.debug(f"First parallactic angle is : {self.parang_start:.2f}")
        logger.debug(f"Middle parallactic angle is : {self.parang_mid:.2f}")
        logger.debug(f"Last parallactic angle is : {self.parang_end:.2f}")

        msmd.close()


    def parallactic_angle(self, HA: float, lat: float, dec: float) -> Union[float,float]:
        """
        Inputs:

        HA      Hour angle, in radian
        lat     Latitude of observation, in deg
        dec     Declination of source, in radian

        Returns:
        eta     Parallactic angle of source, in radian

        Function to calculate the parallactic angle.
        Equations from:
        "A treatise on spherical astronomy" By Sir Robert Stawell Ball
        (p. 91, as viewed on Google Books)

        sin(eta)*sin(z) = cos(lat)*sin(HA)
        cos(eta)*sin(z) = sin(lat)*cos(dec) - cos(lat)*sin(dec)*cos(HA)
        Where eta is the parallactic angle, z is the zenith angle, lat is the
        observer's latitude, dec is the declination, and HA is the hour angle.

        Thus:
        tan(eta) = cos(lat)*sin(HA) / (sin(lat)*cos(dec)-cos(lat)*sin(dec)*cos(HA))

        NOTE : This function returns a parallactic angle that is consistent with
        CASA conventions. However it is 2*pi radians different from the
        parallactic angles calculated using the astropy ecosystem (i.e., via
        astroquery.Observer)
        """

        eta = np.arctan2(
                np.cos(lat)*np.sin(HA),
                (np.sin(lat) * np.cos(dec) - np.cos(lat)*np.sin(dec)*np.cos(HA))
        )

        slope = 1.

        if(eta < 0):
            slope = -1.0
            eta += 2.0*np.pi

        return eta, slope


    def hour_angle(self, ra: float, time: float, observatory:str ='VLA') -> float:
        """
        Function to compute the hourangle of a source at a given time, given an
        observatory name. This function uses the measures tool to compute the
        epoch and the corresponding position

        Inputs:
        ra              Right Ascension of object in radians, float
        time            Timestamp of observation in MJD-seconds, float
        observatory     Name of the observatory, str

        Outputs:
        ha              Hour angle in rad corresponding to the input time stamp, float
        """

        if hasattr(time, '__len__') and (not isinstance(time, str)):
            raise ValueError("Input time should be a single float, not an array/list.")

        # Check if observatory is valid
        obslist = me.obslist()

        if observatory.lower() not in [oo.lower() for oo in obslist]:
            msg = f"Observatory {observatory} is unknown. "
            msg += f"The ability to add a custom observatory will be added in a future update. "
            msg += f"Currently the following observatories are known : {','.join(obslist)}"
            raise ValueError(msg)

        me.doframe(me.observatory(observatory))

        tm = me.epoch('UTC', f'{str(time)}s')
        last = me.measure(tm, 'LAST')['m0']['value']
        last -= np.floor(last)  # days
        last *= 24.0  # hours
        ha = last - ra  # hours
        ha = ha * 2*np.pi/24.  # rad

        return ha


    def parang_weights(self, bin_width:int = 5) -> np.ndarray:
        """

        Bin the calculated parallactic angles and calculate the weight per
        parang bin.

        Inputs:
        bin_width   Width of each bin in degrees, int

        Returns:
        hist        Normalized histogram, i.e., weights per bin, array of floats
        bins        List of bin centres in degrees, array of floats
        """

        if len(self.parangs) == 0:
            raise ValueError("No valid parallactic angles found. Please run "
                                "ParallacticAngle.fill_ms_attributes() first.")

        parang_range = [np.amin(self.parangs), np.amax(self.parangs)]
        parang_bins = np.arange(parang_range[0], parang_range[1], bin_width)

        hist, bins = np.histogram(self.parangs, bins=parang_bins, density=True)

        bins = (bins[1:] + bins[:-1])/2.

        return hist, bins
