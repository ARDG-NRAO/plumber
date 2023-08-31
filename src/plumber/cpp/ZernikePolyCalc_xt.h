//# ZernikePolyCalc_xt.cc: Implementation for ZernikeCalc
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

#ifndef ZERNIKEPOLYCALCXT_H
#define ZERNIKEPOLYCALCXT_H

#include <xtensor/xtensor.hpp>
#include <xtensor/xio.hpp>
#include "ZernikePolyCalc_xt.h"
#include <xtensor/xmath.hpp>

// #include <xtensor/xview.hpp>
// #include <xtensor/xadapt.hpp>
// #include <xtensor/xsort.hpp>
// #include <xtensor/xindex_view.hpp>
// #include <xtensor-python/pyarray.hpp>
// #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define FORCE_IMPORT_ARRAY
// #include <xtensor-python/pytensor.hpp>
// #include <xtensor-python/pyvectorize.hpp>
// #include <xtensor-python/xtensor.hpp>

using namespace std;
using namespace xt;
  class ZernikePolyCalcXT
  {
  public:
    ZernikePolyCalcXT(){};
    ZernikePolyCalcXT(xt::xtensor<double,1> coeffs,  xt::xtensor<double,2> xgrid,  xt::xtensor<double,2> ygrid) {
      coeffs_p = coeffs;
      xgrid_p = xgrid;
      ygrid_p = ygrid;}
      // evaluate (surface_p);}
    ~ZernikePolyCalcXT(){};
    xt::xtensor<double,2> zsurface ( xt::xtensor<double,1>& coeffs,  xt::xtensor<double,2>& xgrid, xt::xtensor<double,2>& ygrid);
    void setcoeffs(xt::xtensor<double,2>& coeffs);
    void setxgrid(xt::xtensor<double,2>& xgrid);
    void setygrid(xt::xtensor<double,2>& ygrid);
    void evaluate(xt::xtensor<double,2>& surface_p) {surface_p = zsurface(coeffs_p, xgrid_p, ygrid_p);};
    xt::xtensor<double,2> getsurface(xt::xtensor<double,2> surface){ evaluate(surface); return surface;};
     // xt::xtensor<double,2> surface_p;
    xt::xtensor<double,1> coeffs_p;
    xt::xtensor<double,2> xgrid_p, ygrid_p, surface_p;

  };
#endif
