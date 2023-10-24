#!/usr/bin/env python3

from plumber.sky import ParallacticAngle
import numpy as np

def test_parang_casa(input_ms='./plumber-data/measurement_sets/J1939-6342_MeerKAT_1400MHz.ms'):
    parang_casa = ParallacticAngle(ms=input_ms, use_astropy=False)
    print(parang_casa.parangs[0])
    print(parang_casa.parangs[-1])

    assert np.allclose(parang_casa.parangs[0], 262.576893)
    assert np.allclose(parang_casa.parangs[-1], 263.558544)

def test_parang_astropy(input_ms='./plumber-data/measurement_sets/J1939-6342_MeerKAT_1400MHz.ms'):
    parang_astro = ParallacticAngle(ms=input_ms, use_astropy=True)
    print(parang_astro.parangs[0])
    print(parang_astro.parangs[-1])

    assert np.allclose(parang_astro.parangs[0], -97.4202206)
    assert np.allclose(parang_astro.parangs[-1], -96.4385573)
