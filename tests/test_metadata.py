#! /usr/bin/env python

"""
Module to test all the metadata-related functionality in plumber.

Currently tests the TelescopeInfo class, but will be expanded to support other
metadata classes.
"""

from plumber.telescope import TelescopeInfo
import numpy as np


def test_meerkat_info():
    ti = TelescopeInfo()
    ti.telescope_name = 'MeerKAT'

    feed_basis = ti.get_feed_basis()
    dish_dia = ti.get_dish_diameter()

    print(f"Feed basis {feed_basis}, Dish_dia {dish_dia}")

    np.testing.assert_string_equal(feed_basis, "linear")
    np.testing.assert_array_equal(dish_dia, [13.5,13.5])


def test_vla_info():
    ti = TelescopeInfo()
    ti.telescope_name = 'VLA'

    feed_basis = ti.get_feed_basis()
    dish_dia = ti.get_dish_diameter()

    np.testing.assert_string_equal(feed_basis, "circular")
    np.testing.assert_array_equal(dish_dia, [25,25])


def test_alma_info():
    ti = TelescopeInfo()
    ti.telescope_name = 'ALMA'

    feed_basis = ti.get_feed_basis()
    dish_dia = ti.get_dish_diameter()

    np.testing.assert_string_equal(feed_basis, "linear")
    np.testing.assert_array_equal(dish_dia, [12,12])


def test_gmrt_info():
    ti = TelescopeInfo()
    ti.telescope_name = 'GMRT'

    feed_basis = ti.get_feed_basis()
    dish_dia = ti.get_dish_diameter()

    np.testing.assert_string_equal(feed_basis, "circular")
    np.testing.assert_array_equal(dish_dia, [45,45])
