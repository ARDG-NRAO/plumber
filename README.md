<div align="center">
  
  # Plumber

  Image plane polarization leakage correction for radio interferometers
  
  [![arXiv](https://img.shields.io/badge/arXiv-2107.10009-b31b1b.svg)](https://arxiv.org/abs/2107.10009)
  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5484098.svg)](https://doi.org/10.5281/zenodo.5484098)

  
</div>
  

Plumber generates full Stokes primary beam models for radio interferometers using Zernike model coefficients of the antenna aperture illumination pattern. The generated PB models are scaled and matched to the input image coordinate system in order to use standard PB correction tools (such as those found in [CASA](https://casadocs.readthedocs.io/en/stable/)).

As of 2021 this supports VLA S Band and MeerKAT L band observations, although support for other telescopes and bands is simply a question of providing the coefficients via the input CSV file.

---

### Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Citation](#citation)



# Installation

> **Note**
> The `spectral-cube` dependency points to the [master branch](https://github.com/radio-astro-tools/spectral-cube.git), pending the stable release of certain bleeding-edge features that are required by `plumber`.

## Installation via pip

Clone the repo, then `cd` into the cloned directory and run

```pip install .```

## Installation via conda
First ensure you have `conda` available from [miniconda](https://docs.conda.io/en/latest/miniconda.html) or [anaconda](https://docs.anaconda.com/anaconda/index.html). [mamba](https://github.com/mamba-org/mamba) is also recommended to speed up the installation process. 

Clone the repo, and then run 

```conda env create plumber```

This will create a conda virtual environment called `plumber`. To use `plumber` simply run

```conda activate plumber```

If you encounter a casa-data error like:

```ImportError: measures data is not available, visit http://go.nrao.edu/casadata-info for more information```

You can fix this by updating your casa-data:

```python -m casatools --update-user-data```

# Usage:

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


### Determining parallactic angle

In order to rotate the beams to the right parallactic angle, the
`parang_finder` helper script is bundled with `plumber`. 

```
Usage: parang_finder [OPTIONS] MS

  Helper script to determine the range of parallactic angles in an input MS.

  Part of the plumber package to generate full Stokes beam models.

Options:
  --field TEXT            Name of field to consider, must match exactly the
                          name in the MS. If this is not specified, will use
                          the first field with the TARGET intent.
  --use-astropy           Use astropy convention rather than CASA convention
                          to   calculate the parallactic angle
  -c, --calcweights       Print a list of parallactic angle and associated
                          weights
  -w, --binwidth INTEGER  Bin width in degrees
  -h, --help              Show this message and exit.
```

### Location of coefficients

The coefficient files are within the `./plumber/data/` directory of this
repository. When the repository is cloned locally, the path of the CSV files can be passed
in via command line.


# Citation

If you have used `plumber` in your research, please cite the following paper : 

Sekhar, S., Jagannathan, P., Kirk, B., Bhatnagar, S., & Taylor, R. (2022). Direction-dependent Corrections in Polarimetric Radio Imaging. III. A-to-Z Solver—Modeling the Full Jones Antenna Aperture Illumination Pattern. The Astronomical Journal, 163(2), 87.

An example BibTeX entry is below : 
```
@article{Sekhar_2022,
doi = {10.3847/1538-3881/ac41c4},
url = {https://dx.doi.org/10.3847/1538-3881/ac41c4},
year = {2022},
month = {jan},
publisher = {The American Astronomical Society},
volume = {163},
number = {2},
pages = {87},
author = {Srikrishna Sekhar and Preshanth Jagannathan and Brian Kirk and Sanjay Bhatnagar and Russ Taylor},
title = {Direction-dependent Corrections in Polarimetric Radio Imaging. III. A-to-Z Solver—Modeling the Full Jones Antenna Aperture Illumination Pattern},
journal = {The Astronomical Journal},
}
```
