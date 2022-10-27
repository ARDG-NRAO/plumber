#!/usr/bin/env python3

import os
import shutil

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
