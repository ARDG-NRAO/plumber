//# ZernikePolyCalc_xt2py.cc: Implementation for ZernikeCalc
//# Copyright (C) 1996,1997,1998,1999,2000,2002
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be adressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//#
//# $Id$
// #include <xtensor/xarray.hpp>
// #include <xtensor/xio.hpp>
#include "ZernikePolyCalc_xt.h"
// #include <xtensor/xmath.hpp>
// #include <xtensor/xtensor.hpp>
// #include <xtensor/xview.hpp>
// #include <xtensor/xadapt.hpp>
// #include <xtensor/xsort.hpp>
// #include <xtensor/xindex_view.hpp>
#define FORCE_IMPORT_ARRAY
#include <xtensor-python/pytensor.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace std;

xt::pytensor<double,2> ZernikePolyCalc_xt2py (xt::pytensor<double,1>& coeffs, xt::pytensor<double,2>& xgrid, xt::pytensor<double,2>& ygrid) {
  int ncoeffs = coeffs.shape()[0];
  int ngrid = xgrid.shape()[0];
  xt::pytensor<double,2> surface = xt::zeros<double>({ngrid, ngrid});
  ZernikePolyCalcXT *zpcxt;
  zpcxt = new ZernikePolyCalcXT(coeffs, xgrid, ygrid);
//   zpcxt.zernikesurface(coeffs, xgrid, ygrid);
  surface = zpcxt->getsurface(surface);
  return surface;
}
PYBIND11_MODULE(ZernikePolyCalc_xt2py, m) {
    xt::import_numpy();
    m.def("ZernikePolyCalc_xt2py", &ZernikePolyCalc_xt2py, "Calculate Zernike surface");
    }

