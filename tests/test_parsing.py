#!/usr/bin/env python3

from plumber.parsing import CSVParser

import pdb
import os
import numpy as np
import pandas as pd

import astropy.units as u

import pytest

mock_csv = """
#stokes,freq,ind,real,imag,etax,etay
9,1000,0,1,1,1,1
10,1000,0,1,1,1,1
11,1000,0,1,1,1,1
12,1000,0,1,1,1,1
"""

def test_csv_to_df():
   """
   Test that the CSV to DataFrame parsing works correctly.
   """

   tmpfile = 'tmp_parsing.txt'
   with open(tmpfile, 'w') as fptr:
      fptr.writelines(mock_csv)

   parser = CSVParser()
   df = parser.csv_to_df(tmpfile)

   os.remove(tmpfile)

   assert type(df) is pd.DataFrame


def test_set_filepath():
   """
   Test that the filepath method does the right thing.
   """

   tmpfile = 'tmp_parsing.txt'
   with open(tmpfile, 'w') as fptr:
      fptr.writelines(mock_csv)

   parser = CSVParser()
   parser.filepath = tmpfile

   assert type(parser.filepath) is str


def test_set_filepath_exception():
   """
   Test that the filepath method throws an exception when there is no file
   passed in.
   """

   parser = CSVParser()

   with pytest.raises(FileNotFoundError):
      parser.filepath = mock_csv


def test_imfreq_units_exception_scalar():
   """
   Test that the imfreq exception raising works.
   """

   tmpfile = 'tmp_parsing.txt'
   with open(tmpfile, 'w') as fptr:
      fptr.writelines(mock_csv)

   parser = CSVParser()
   df = parser.csv_to_df(tmpfile)

   os.remove(tmpfile)

   # imfreq is not astropy.unit so should throw exception
   with pytest.raises(TypeError):
      outdf, outfreq, nstokes = parser.get_zcoeffs(df, 1000)


def test_imfreq_units_exception_vector():
   """
   Test that the imfreq exception raising works.
   """

   tmpfile = 'tmp_parsing.txt'
   with open(tmpfile, 'w') as fptr:
      fptr.writelines(mock_csv)

   parser = CSVParser()
   df = parser.csv_to_df(tmpfile)

   os.remove(tmpfile)

   # imfreq is not astropy.unit so should throw exception
   with pytest.raises(TypeError):
      outdf, outfreq, nstokes = parser.get_zcoeffs(df, [1000,2000])



def test_imfreq_units_vector():
   """
   Test that the imfreq works for vectors
   """

   tmpfile = 'tmp_parsing.txt'
   with open(tmpfile, 'w') as fptr:
      fptr.writelines(mock_csv)

   parser = CSVParser()
   df = parser.csv_to_df(tmpfile)

   os.remove(tmpfile)

   freqs = np.arange(1000, 2000, 500) * u.MHz
   outdf, outfreq, nstokes = parser.get_zcoeffs(df, freqs)
   #pdb.set_trace()

   assert outfreq == [1000, 1000]
   assert nstokes == 4
   assert len(outdf) == 2
   assert outdf[0].equals(df)
