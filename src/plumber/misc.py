#!/usr/bin/env python3

import os
import shutil
import uuid

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)-15s %(levelname)s: %(message)s",
                        level=logging.INFO)

def wipe_file(filename):
    """
    Check if the input file exists, and if so, delete it.

    Input:
    filename    Input file name, string

    Returns:
    None
    """

    if os.path.exists(filename):
        shutil.rmtree(filename)


def make_unique(filename):
    """
    Given an input filename, append a random hexadecimal string to the end.
    This allows for multiple independent parallel runs to avoid naming conflicts
    of temporary files created within plumber.

    Input:
    filename    Input file name, string

    Returns:
    uniquename  Output name, with hex string appendend, string
    """

    hex = uuid.uuid4().hex
    return filename + hex
