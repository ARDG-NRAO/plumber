#!/usr/bin/env python3

import shutil
import numpy as np
import os

from plumber.scripts.plumber import main as plumber_main

from casatools import image
ia = image()

def get_imdat(imname):
    ia.open(imname)
    dat =  np.squeeze(ia.getchunk())
    ia.close()

    return dat


def test_StokesI_generation():
    """
    Given a fixed CSV file, compare against a known beam and confirm the PB
    generation is working.
    """

    template_image = './plumber-data/images/I_MeerKAT_1577.832MHz_template.im'
    csv_path = './plumber-data/csv/MeerKAT_coeffs_for_testing.csv'

    try:
        # Need to pass args as a list since plumnber.main is wrapped with Click
        plumber_main([f'{template_image}', f'{csv_path}', '--stokesI'])
    except SystemExit:
        pass

    imname = 'I_MeerKAT_1577.832MHz.im'

    imdat_gen = get_imdat(imname)
    imdat_template = get_imdat(template_image)

    shutil.rmtree('I_MeerKAT_1577.832MHz.im')

    np.testing.assert_allclose(imdat_gen, imdat_template, rtol=1e-06, atol=1e-06)
