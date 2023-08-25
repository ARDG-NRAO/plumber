//# ZernikePolyCalc2py.h: Definition for Pybind11 wrapper for ZernikePolyCalc
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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "ZernikePolyCalc.h"
#include <string>
#include <iostream>

using namespace std;
using namespace pybind11;
using namespace pybind11::literals;
namespace py = pybind11;

PYBIND11_MODULE(ZernikePolyCalc2py, m)
{
  m.doc() = "pybind11-based wrapper for the ZenikePolyCalc class in C++"; // optional module docstring
  py::class_<ZernikePolyCalc>(m, "ZernikePolyCalc")
  .def(py::init<>())
  .def("powl", &ZernikePolyCalc::powl, "base"_a, "exp"_a)
  .def("zernikeSurface", &ZernikePolyCalc::zernikeSurface, "amp"_a, "xCoords"_a, "yCoords"_a, "xSize"_a, "ySize"_a, "surface"_a)
  .def("ncoeffs", &ZernikePolyCalc::ncoeffs, "Return number of coefficients")
  .def("getCoeffs", &ZernikePolyCalc::getCoeffs, "Return coefficients")
  .def("getZPolyvals", &ZernikePolyCalc::getZPolyvals, "Return Zernike polynomial values")
  .def_readwrite("Z_p", &ZernikePolyCalc::Z_p, "Zernike polynomial values")
  .def_readwrite("C_p", &ZernikePolyCalc::C_p, "Zernike polynomial coefficients")
  .def("__repr__",
    [](const ZernikePolyCalc &a) {
        return "<zernikepolycalc2py.ZernikePolyCalc>";
    }
  );

}