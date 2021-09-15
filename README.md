# Plumber

[![DOI](https://zenodo.org/badge/402479973.svg)](https://zenodo.org/badge/latestdoi/402479973)

Image plane polarization leakage correction. Uses Zernike models of antenna
apertures to generate the primary beam across the relevant field of view.

Requires modular CASA 6 (casatools and casatasks) and the relevant coefficient
CSV file for the telescope and band in question.

As of June 2021 this supports VLA S Band and MeerKAT L band observations.

## Installation
Clone the repo, and then run 

```pip install -e setup.py```


## Usage:

```plumber <path_to_template_image> <path_to_coefficient_csv> [options]```


The full command line arguments :

```
Usage: plumber [OPTIONS] IMAGENAME CSV

  Given the input image and the coefficient CSV file, generate the full
  Stokes primary beam at the image centre frequency. If the input is a cube,
  a corresponding output cube of PBs will be generated, for each of the
  Stokes-I, -Q, -U and -V beams.

  To pass in a beginning and end parallactic angle, pass in --parang twice,
  like

  plumber <image> <CSV> --parang 20 --parang 110

  The CSV file must be of the format

  #stokes,freq,ind,real,imag,eta

  The eta column in the CSV is optional.

Options:
  -a, --padding INTEGER  Padding factor for aperture, affects smoothness of
                         output beam  [default: 8]

  -d, --dish_dia FLOAT   Diameter of the antenna dish. If not one of VLA,
                         ALMA, MeerKAT or GMRT, must be specified.

  -l, --islinear         Specifies if the telescope has linear feeds. If not
                         one of VLA, ALMA, MeerKAT or GMRT, must be specified

  -I, --stokesI          Only generate the Stokes I beam, not the full Stokes
                         beams

  -P, --parallel         Use parallel processing (no MPI) to speed things up
  -p, --parang FLOAT     Beginning (and optionally end) parallactic angle for
                         the PB. Pass --parang twice to specify beg and end.

  -h, --help             Show this message and exit.
```


## Determining parallactic angle

In order to rotate the beams to the right parallactic angle, the
`parang_finder` helper script is bundled with `plumber`. 

```
Usage: parang_finder [OPTIONS] MS

  Helper script to determine the range of parallactic angles in an input MS.

  Part of the plumber package to generate full Stokes beam models.

Options:
  --field TEXT  Name of field to consider, must match exactly the name in the
                MS. If this is not specified, will use the first field with
                the TARGET intent.

  -h, --help    Show this message and exit.
```

## Location of coefficients

The coefficient files are within the `plumber/data/` directory of this
repository. When the repository is cloned locally, the CSV files can be passed
in via command line.
