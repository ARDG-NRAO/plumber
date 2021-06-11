# Plumber

Image plane polarization leakage correction. Uses Zernike models of antenna
apertures to generate the primary beam across the relevant field of view.

Requires modular CASA 6 (casatools and casatasks) and the relevant coefficient
CSV file for the telescope and band in question.

As of June 2021 this supports VLA S Band and MeerKAT L band observations.
