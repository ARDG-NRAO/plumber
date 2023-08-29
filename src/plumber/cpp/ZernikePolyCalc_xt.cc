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
#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include "ZernikePolyCalc_xt.h"
#include <xtensor/xmath.hpp>
#include <xtensor/xview.hpp>
// #include <xtensor/xadapt.hpp>
// #include <xtensor/xsort.hpp>
// #include <xtensor/xindex_view.hpp>
// #include <xtensor-python/pyarray.hpp>
// #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
// #define FORCE_IMPORT_ARRAY
// #include <xtensor-python/pytensor.hpp>
// #include <xtensor-python/pyvectorize.hpp>
// #include <xtensor-python/xtensor.hpp>

using namespace std;

    // Calculate the Zernike polynomials for a given set of coefficients based on the function gen_zernike_surface in zernike.py to xtensor c++
    // coeffs: array of Zernike coefficients
    // xgrid: x coordinates of the grid
    // ygrid: y coordinates of the grid
    // return: Zernike surface
    xt::xtensor<double,2> ZernikePolyCalcXT::zernikesurface (xt::xtensor<double,1>& coeffs, xt::xtensor<double,2>& xgrid, xt::xtensor<double,2>& ygrid)
    {
        xt::xtensor<double,1> c = xt::zeros<double>({67});
        auto coeffs_shape = coeffs.size();
        auto c_shape = c.size();
        if (coeffs_shape <= c_shape){
            c += coeffs;
        }
        xt::xtensor<double,2> x = xgrid;
        xt::xtensor<double,2> y = ygrid;
        auto x_shape = xgrid.shape();
        auto unit_vector = xt::ones<double>({x_shape[0], x_shape[1]});
        xt::xtensor<double,2> Z1  =  c(0) * unit_vector;//#m = 0    n = 0
        // xt::xtensor<double,2> Z1  =  c( 0) * xt::ones<double> ({x_shape, x_shape}) ;//#m = 0    n = 0
        xt::xtensor<double,2> Z2  =  c( 1) * x ;//#m = -1   n = 1
        xt::xtensor<double,2> Z3  =  c( 2) * y ;//#m = 1    n = 1
        xt::xtensor<double,2> Z4  =  c( 3) * 2*x*y ;//#m = -2   n = 2
        xt::xtensor<double,2> Z5  =  c( 4) * (2*xt::pow(x,2) + 2*xt::pow(y,2) -1);//#m = 0  n = 2
        xt::xtensor<double,2> Z6  =  c( 5) * (-1*xt::pow(x,2) + xt::pow(y,2));//#m = 2  n = 2
        xt::xtensor<double,2> Z7  =  c( 6) * (-1*xt::pow(x,3) + 3*x*xt::pow(y,2)) ;//#m = -3     n = 3
        xt::xtensor<double,2> Z8  =  c( 7) * (-2*x + 3*(xt::pow(x,3)) + 3*x*(xt::pow(y,2)));//#m = -1   n = 3
        xt::xtensor<double,2> Z9  =  c( 8) * (-2*y + 3*xt::pow(y,3) + 3*(xt::pow(x,2)) * y);//#m = 1    n = 3
        xt::xtensor<double,2> Z10 =  c( 9) * (xt::pow(y,3)-3*(xt::pow(x,2))*y);//#m = 3 n =3
        xt::xtensor<double,2> Z11 =  c( 10) * (-4*(xt::pow(x,3))*y + 4*x*(xt::pow(y,3)));//#m = -4    n = 4
        xt::xtensor<double,2> Z12 =  c( 11) * (-6*x*y + 8*(xt::pow(x,3))*y + 8*x*(xt::pow(y,3)));//#m = -2   n = 4
        xt::xtensor<double,2> Z13 =  c( 12) * (1-6*xt::pow(x,2) - 6*xt::pow(y,2) + 6*xt::pow(x,4) + 12*(xt::pow(x,2))*(xt::pow(y,2)) + 6*xt::pow(y,4));//#m = 0  n = 4
        xt::xtensor<double,2> Z14 =  c( 13) * (3*xt::pow(x,2) - 3*xt::pow(y,2) - 4*xt::pow(x,4) + 4*xt::pow(y,4));//#m = 2    n = 4
        xt::xtensor<double,2> Z15 =  c( 14) * (xt::pow(x,4) - 6*(xt::pow(x,2))*(xt::pow(y,2)) + xt::pow(y,4));//#m = 4   n = 4
        xt::xtensor<double,2> Z16 =  c( 15) * (xt::pow(x,5)-10*(xt::pow(x,3))*xt::pow(y,2) + 5*x*(xt::pow(y,4)));//#m = -5   n = 5
        xt::xtensor<double,2> Z17 =  c( 16) * (4*xt::pow(x,3) - 12*x*(xt::pow(y,2)) -5*xt::pow(x,5) + 10*(xt::pow(x,3))*(xt::pow(y,2)) + 15*x*xt::pow(y,4));//#m =-3     n = 5
        xt::xtensor<double,2> Z18 =  c( 17) * (3*x - 12*xt::pow(x,3) - 12*x*(xt::pow(y,2)) + 10*xt::pow(x,5) + 20*(xt::pow(x,3))*(xt::pow(y,2)) + 10*x*(xt::pow(y,4)));//#m= -1  n = 5
        xt::xtensor<double,2> Z19 =  c( 18) * (3*y - 12*xt::pow(y,3) - 12*y*(xt::pow(x,2)) + 10*xt::pow(y,5) + 20*(xt::pow(y,3))*(xt::pow(x,2)) + 10*y*(xt::pow(x,4)));//#m = 1  n = 5
        xt::xtensor<double,2> Z20 =  c( 19) * (-4*xt::pow(y,3) + 12*y*(xt::pow(x,2)) + 5*xt::pow(y,5) - 10*(xt::pow(y,3))*(xt::pow(x,2)) - 15*y*xt::pow(x,4));//#m = 3   n = 5
        xt::xtensor<double,2> Z21 =  c( 20) * (xt::pow(y,5)-10*(xt::pow(y,3))*xt::pow(x,2) + 5*y*(xt::pow(x,4)));//#m = 5 n = 5
        xt::xtensor<double,2> Z22 =  c( 21) * (6*(xt::pow(x,5))*y - 20*(xt::pow(x,3))*(xt::pow(y,3)) + 6*x*(xt::pow(y,5)));//#m = -6 n = 6
        xt::xtensor<double,2> Z23 =  c( 22) * (20*(xt::pow(x,3))*y - 20*x*(xt::pow(y,3)) - 24*(xt::pow(x,5))*y + 24*x*(xt::pow(y,5)));//#m = -4   n = 6
        xt::xtensor<double,2> Z24 =  c( 23) * (12*x*y + 40*(xt::pow(x,3))*y - 40*x*(xt::pow(y,3)) + 30*(xt::pow(x,5))*y + 60*(xt::pow(x,3))*(xt::pow(y,3)) - 30*x*(xt::pow(y,5)));//#m = -2   n = 6
        xt::xtensor<double,2> Z25 =  c( 24) * (-1 + 12*(xt::pow(x,2)) + 12*(xt::pow(y,2)) - 30*(xt::pow(x,4)) - 60*(xt::pow(x,2))*(xt::pow(y,2)) - 30*(xt::pow(y,4)) + 20*(xt::pow(x,6)) + 60*(xt::pow(x,4))*xt::pow(y,2) + 60 *(xt::pow(x,2))*(xt::pow(y,4)) + 20*(xt::pow(y,6)));//#m = 0   n = 6
        xt::xtensor<double,2> Z26 =  c( 25) * (-6*(xt::pow(x,2)) + 6*(xt::pow(y,2)) + 20*(xt::pow(x,4)) - 20*(xt::pow(y,4)) - 15*(xt::pow(x,6)) - 15*(xt::pow(x,4))*(xt::pow(y,2)) + 15*(xt::pow(x,2))*(xt::pow(y,4)) + 15*(xt::pow(y,6)));//#m = 2   n = 6
        xt::xtensor<double,2> Z27 =  c( 26) * (-5*(xt::pow(x,4)) + 30*(xt::pow(x,2))*(xt::pow(y,2)) - 5*(xt::pow(y,4)) + 6*(xt::pow(x,6)) - 30*(xt::pow(x,4))*xt::pow(y,2) - 30*(xt::pow(x,2))*(xt::pow(y,4)) + 6*(xt::pow(y,6)));//#m = 4    n = 6
        xt::xtensor<double,2> Z28 =  c( 27) * (-1*(xt::pow(x,6)) + 15*(xt::pow(x,4))*(xt::pow(y,2)) - 15*(xt::pow(x,2))*(xt::pow(y,4)) + xt::pow(y,6));//#m = 6   n = 6
        xt::xtensor<double,2> Z29 =  c( 28) * (-1*(xt::pow(x,7)) + 21*(xt::pow(x,5))*(xt::pow(y,2)) - 35*(xt::pow(x,3))*(xt::pow(y,4)) + 7*x*(xt::pow(y,6)));//#m = -7    n = 7
        xt::xtensor<double,2> Z30 =  c( 29) * (-6*(xt::pow(x,5)) + 60*(xt::pow(x,3))*(xt::pow(y,2)) - 30*x*(xt::pow(y,4)) + 7*xt::pow(x,7) - 63*(xt::pow(x,5))*(xt::pow(y,2)) - 35*(xt::pow(x,3))*(xt::pow(y,4)) + 35*x*(xt::pow(y,6))) ;//#m = -5    n = 7
        xt::xtensor<double,2> Z31 =  c( 30) * (-10*(xt::pow(x,3)) + 30*x*(xt::pow(y,2)) + 30*xt::pow(x,5) - 60*(xt::pow(x,3))*(xt::pow(y,2)) - 90*x*(xt::pow(y,4)) - 21*xt::pow(x,7) + 21*(xt::pow(x,5))*(xt::pow(y,2)) + 105*(xt::pow(x,3))*(xt::pow(y,4)) + 63*x*(xt::pow(y,6)));//#m =-3       n = 7
        xt::xtensor<double,2> Z32 =  c( 31) * (-4*x + 30*xt::pow(x,3) + 30*x*(xt::pow(y,2)) - 60*(xt::pow(x,5)) - 120*(xt::pow(x,3))*(xt::pow(y,2)) - 60*x*(xt::pow(y,4)) + 35*xt::pow(x,7) + 105*(xt::pow(x,5))*(xt::pow(y,2)) + 105*(xt::pow(x,3))*(xt::pow(y,4)) + 35*x*(xt::pow(y,6)));//#m = -1  n = 7
        xt::xtensor<double,2> Z33 =  c( 32) * (-4*y + 30*xt::pow(y,3) + 30*y*(xt::pow(x,2)) - 60*(xt::pow(y,5)) - 120*(xt::pow(y,3))*(xt::pow(x,2)) - 60*y*(xt::pow(x,4)) + 35*xt::pow(y,7) + 105*(xt::pow(y,5))*(xt::pow(x,2)) + 105*(xt::pow(y,3))*(xt::pow(x,4)) + 35*y*(xt::pow(x,6)));//#m = 1   n = 7
        xt::xtensor<double,2> Z34 =  c( 33) * (10*(xt::pow(y,3)) - 30*y*(xt::pow(x,2)) - 30*xt::pow(y,5) + 60*(xt::pow(y,3))*(xt::pow(x,2)) + 90*y*(xt::pow(x,4)) + 21*xt::pow(y,7) - 21*(xt::pow(y,5))*(xt::pow(x,2)) - 105*(xt::pow(y,3))*(xt::pow(x,4)) - 63*y*(xt::pow(x,6)));//#m =3     n = 7
        xt::xtensor<double,2> Z35 =  c( 34) * (-6*(xt::pow(y,5)) + 60*(xt::pow(y,3))*(xt::pow(x,2)) - 30*y*(xt::pow(x,4)) + 7*xt::pow(y,7) - 63*(xt::pow(y,5))*(xt::pow(x,2)) - 35*(xt::pow(y,3))*(xt::pow(x,4)) + 35*y*(xt::pow(x,6)));//#m = 5  n = 7
        xt::xtensor<double,2> Z36 =  c( 35) * (xt::pow(y,7) - 21*(xt::pow(y,5))*(xt::pow(x,2)) + 35*(xt::pow(y,3))*(xt::pow(x,4)) - 7*y*(xt::pow(x,6)));//#m = 7  n = 7
        xt::xtensor<double,2> Z37 =  c( 36) * (-8*(xt::pow(x,7))*y + 56*(xt::pow(x,5))*(xt::pow(y,3)) - 56*(xt::pow(x,3))*(xt::pow(y,5)) + 8*x*(xt::pow(y,7)));//#m = -8  n = 8
        xt::xtensor<double,2> Z38 =  c( 37) * (-42*(xt::pow(x,5))*y + 140*(xt::pow(x,3))*(xt::pow(y,3)) - 42*x*(xt::pow(y,5)) + 48*(xt::pow(x,7))*y - 112*(xt::pow(x,5))*(xt::pow(y,3)) - 112*(xt::pow(x,3))*(xt::pow(y,5)) + 48*x*(xt::pow(y,7)));//#m = -6  n = 8
        xt::xtensor<double,2> Z39 =  c( 38) * (-60*(xt::pow(x,3))*y + 60*x*(xt::pow(y,3)) + 168*(xt::pow(x,5))*y -168*x*(xt::pow(y,5)) - 112*(xt::pow(x,7))*y - 112*(xt::pow(x,5))*(xt::pow(y,3)) + 112*(xt::pow(x,3))*(xt::pow(y,5)) + 112*x*(xt::pow(y,7)));//#m = -4   n = 8
        xt::xtensor<double,2> Z40 =  c( 39) * (-20*x*y + 120*(xt::pow(x,3))*y + 120*x*(xt::pow(y,3)) - 210*(xt::pow(x,5))*y - 420*(xt::pow(x,3))*(xt::pow(y,3)) - 210*x*(xt::pow(y,5)) - 112*(xt::pow(x,7))*y + 336*(xt::pow(x,5))*(xt::pow(y,3)) + 336*(xt::pow(x,3))*(xt::pow(y,5)) + 112*x*(xt::pow(y,7)));//#m = -2   n = 8
        xt::xtensor<double,2> Z41 =  c( 40) * (1 - 20*xt::pow(x,2) - 20*xt::pow(y,2) + 90*xt::pow(x,4) + 180*(xt::pow(x,2))*(xt::pow(y,2)) + 90*xt::pow(y,4) - 140*xt::pow(x,6) - 420*(xt::pow(x,4))*(xt::pow(y,2)) - 420*(xt::pow(x,2))*(xt::pow(y,4)) - 140*(xt::pow(y,6)) + 70*xt::pow(x,8) + 280*(xt::pow(x,6))*(xt::pow(y,2)) + 420*(xt::pow(x,4))*(xt::pow(y,4)) + 280*(xt::pow(x,2))*(xt::pow(y,6)) + 70*xt::pow(y,8));//#m = 0    n = 8
        xt::xtensor<double,2> Z42 =  c( 41) * (10*xt::pow(x,2) - 10*xt::pow(y,2) - 60*xt::pow(x,4) + 105*(xt::pow(x,4))*(xt::pow(y,2)) - 105*(xt::pow(x,2))*(xt::pow(y,4)) + 60*xt::pow(y,4) + 105*xt::pow(x,6) - 105*xt::pow(y,6) - 56*xt::pow(x,8) - 112*(xt::pow(x,6))*(xt::pow(y,2)) + 112*(xt::pow(x,2))*(xt::pow(y,6)) + 56*xt::pow(y,8));//#m = 2  n = 8
        xt::xtensor<double,2> Z43 =  c( 42) * (15*xt::pow(x,4) - 90*(xt::pow(x,2))*(xt::pow(y,2)) + 15*xt::pow(y,4) - 42*xt::pow(x,6) + 210*(xt::pow(x,4))*(xt::pow(y,2)) + 210*(xt::pow(x,2))*(xt::pow(y,4)) - 42*xt::pow(y,6) + 28*xt::pow(x,8) - 112*(xt::pow(x,6))*(xt::pow(y,2)) - 280*(xt::pow(x,4))*(xt::pow(y,4)) - 112*(xt::pow(x,2))*(xt::pow(y,6)) + 28*xt::pow(y,8));//#m = 4     n = 8
        xt::xtensor<double,2> Z44 =  c( 43) * (7*xt::pow(x,6) - 105*(xt::pow(x,4))*(xt::pow(y,2)) + 105*(xt::pow(x,2))*(xt::pow(y,4)) - 7*xt::pow(y,6) - 8*xt::pow(x,8) + 112*(xt::pow(x,6))*(xt::pow(y,2)) - 112*(xt::pow(x,2))*(xt::pow(y,6)) + 8*xt::pow(y,8));//#m = 6    n = 8
        xt::xtensor<double,2> Z45 =  c( 44) * (xt::pow(x,8) - 28*(xt::pow(x,6))*(xt::pow(y,2)) + 70*(xt::pow(x,4))*(xt::pow(y,4)) - 28*(xt::pow(x,2))*(xt::pow(y,6)) + xt::pow(y,8));//#m = 8     n = 9
        xt::xtensor<double,2> Z46 =  c( 45) * (xt::pow(x,9) - 36*(xt::pow(x,7))*(xt::pow(y,2)) + 126*(xt::pow(x,5))*(xt::pow(y,4)) - 84*(xt::pow(x,3))*(xt::pow(y,6)) + 9*x*(xt::pow(y,8)));//#m = -9     n = 9
        xt::xtensor<double,2> Z47 =  c( 46) * (8*xt::pow(x,7) - 168*(xt::pow(x,5))*(xt::pow(y,2)) + 280*(xt::pow(x,3))*(xt::pow(y,4)) - 56 *x*(xt::pow(y,6)) - 9*xt::pow(x,9) + 180*(xt::pow(x,7))*(xt::pow(y,2)) - 126*(xt::pow(x,5))*(xt::pow(y,4)) - 252*(xt::pow(x,3))*(xt::pow(y,6)) + 63*x*(xt::pow(y,8)));//#m = -7    n = 9
        xt::xtensor<double,2> Z48 =  c( 47) * (21*xt::pow(x,5) - 210*(xt::pow(x,3))*(xt::pow(y,2)) + 105*x*(xt::pow(y,4)) - 56*xt::pow(x,7) + 504*(xt::pow(x,5))*(xt::pow(y,2)) + 280*(xt::pow(x,3))*(xt::pow(y,4)) - 280*x*(xt::pow(y,6)) + 36*xt::pow(x,9) - 288*(xt::pow(x,7))*(xt::pow(y,2)) - 504*(xt::pow(x,5))*(xt::pow(y,4)) + 180*x*(xt::pow(y,8)));//#m = -5    n = 9
        xt::xtensor<double,2> Z49 =  c( 48) * (20*xt::pow(x,3) - 60*x*(xt::pow(y,2)) - 105*xt::pow(x,5) + 210*(xt::pow(x,3))*(xt::pow(y,2)) + 315*x*(xt::pow(y,4)) + 168*xt::pow(x,7) - 168*(xt::pow(x,5))*(xt::pow(y,2)) - 840*(xt::pow(x,3))*(xt::pow(y,4)) - 504*x*(xt::pow(y,6)) - 84*xt::pow(x,9) + 504*(xt::pow(x,5))*(xt::pow(y,4)) + 672*(xt::pow(x,3))*(xt::pow(y,6)) + 252*x*(xt::pow(y,8)));//#m = -3  n = 9
        xt::xtensor<double,2> Z50 =  c( 49) * (5*x - 60*xt::pow(x,3) - 60*x*(xt::pow(y,2)) + 210*xt::pow(x,5) + 420*(xt::pow(x,3))*(xt::pow(y,2)) + 210*x*(xt::pow(y,4)) - 280*xt::pow(x,7) - 840*(xt::pow(x,5))*(xt::pow(y,2)) - 840*(xt::pow(x,3))*(xt::pow(y,4)) - 280*x*(xt::pow(y,6)) + 126*xt::pow(x,9) + 504*(xt::pow(x,7))*(xt::pow(y,2)) + 756*(xt::pow(x,5))*(xt::pow(y,4)) + 504*(xt::pow(x,3))*(xt::pow(y,6)) + 126*x*(xt::pow(y,8)));//#m = -1   n = 9
        xt::xtensor<double,2> Z51 =  c( 50) * (5*y - 60*xt::pow(y,3) - 60*y*(xt::pow(x,2)) + 210*xt::pow(y,5) + 420*(xt::pow(y,3))*(xt::pow(x,2)) + 210*y*(xt::pow(x,4)) - 280*xt::pow(y,7) - 840*(xt::pow(y,5))*(xt::pow(x,2)) - 840*(xt::pow(y,3))*(xt::pow(x,4)) - 280*y*(xt::pow(x,6)) + 126*xt::pow(y,9) + 504*(xt::pow(y,7))*(xt::pow(x,2)) + 756*(xt::pow(y,5))*(xt::pow(x,4)) + 504*(xt::pow(y,3))*(xt::pow(x,6)) + 126*y*(xt::pow(x,8)));//#m = -1   n = 9
        xt::xtensor<double,2> Z52 =  c( 51) * (-20*xt::pow(y,3) + 60*y*(xt::pow(x,2)) + 105*xt::pow(y,5) - 210*(xt::pow(y,3))*(xt::pow(x,2)) - 315*y*(xt::pow(x,4)) - 168*xt::pow(y,7) + 168*(xt::pow(y,5))*(xt::pow(x,2)) + 840*(xt::pow(y,3))*(xt::pow(x,4)) + 504*y*(xt::pow(x,6)) + 84*xt::pow(y,9) - 504*(xt::pow(y,5))*(xt::pow(x,4)) - 672*(xt::pow(y,3))*(xt::pow(x,6)) - 252*y*(xt::pow(x,8)));//#m = 3  n = 9
        xt::xtensor<double,2> Z53 =  c( 52) * (21*xt::pow(y,5) - 210*(xt::pow(y,3))*(xt::pow(x,2)) + 105*y*(xt::pow(x,4)) - 56*xt::pow(y,7) + 504*(xt::pow(y,5))*(xt::pow(x,2)) + 280*(xt::pow(y,3))*(xt::pow(x,4)) - 280*y*(xt::pow(x,6)) + 36*xt::pow(y,9) - 288*(xt::pow(y,7))*(xt::pow(x,2)) - 504*(xt::pow(y,5))*(xt::pow(x,4)) + 180*y*(xt::pow(x,8)));//#m = 5     n = 9
        xt::xtensor<double,2> Z54 =  c( 53) * (-8*xt::pow(y,7) + 168*(xt::pow(y,5))*(xt::pow(x,2)) - 280*(xt::pow(y,3))*(xt::pow(x,4)) + 56 *y*(xt::pow(x,6)) + 9*xt::pow(y,9) - 180*(xt::pow(y,7))*(xt::pow(x,2)) + 126*(xt::pow(y,5))*(xt::pow(x,4)) - 252*(xt::pow(y,3))*(xt::pow(x,6)) - 63*y*(xt::pow(x,8)));//#m = 7     n = 9
        xt::xtensor<double,2> Z55 =  c( 54) * (xt::pow(y,9) - 36*(xt::pow(y,7))*(xt::pow(x,2)) + 126*(xt::pow(y,5))*(xt::pow(x,4)) - 84*(xt::pow(y,3))*(xt::pow(x,6)) + 9*y*(xt::pow(x,8)));//#m = 9       n = 9
        xt::xtensor<double,2> Z56 =  c( 55) * (10*(xt::pow(x,9))*y - 120*(xt::pow(x,7))*(xt::pow(y,3)) + 252*(xt::pow(x,5))*(xt::pow(y,5)) - 120*(xt::pow(x,3))*(xt::pow(y,7)) + 10*x*(xt::pow(y,9)));//#m = -10   n = 10
        xt::xtensor<double,2> Z57 =  c( 56) * (72*(xt::pow(x,7))*y - 504*(xt::pow(x,5))*(xt::pow(y,3)) + 504*(xt::pow(x,3))*(xt::pow(y,5)) - 72*x*(xt::pow(y,7)) - 80*(xt::pow(x,9))*y + 480*(xt::pow(x,7))*(xt::pow(y,3)) - 480*(xt::pow(x,3))*(xt::pow(y,7)) + 80*x*(xt::pow(y,9)));//#m = -8    n = 10
        xt::xtensor<double,2> Z58 =  c( 57) * (270*(xt::pow(x,9))*y - 360*(xt::pow(x,7))*(xt::pow(y,3)) - 1260*(xt::pow(x,5))*(xt::pow(y,5)) - 360*(xt::pow(x,3))*(xt::pow(y,7)) + 270*x*(xt::pow(y,9)) - 432*(xt::pow(x,7))*y + 1008*(xt::pow(x,5))*(xt::pow(y,3)) + 1008*(xt::pow(x,3))*(xt::pow(y,5)) - 432*x*(xt::pow(y,7)) + 168*(xt::pow(x,5))*y - 560*(xt::pow(x,3))*(xt::pow(y,3)) + 168*x*(xt::pow(y,5)));//#m = -6   n = 10
        xt::xtensor<double,2> Z59 =  c( 58) * (140*(xt::pow(x,3))*y - 140*x*(xt::pow(y,3)) - 672*(xt::pow(x,5))*y + 672*x*(xt::pow(y,5)) + 1008*(xt::pow(x,7))*y + 1008*(xt::pow(x,5))*(xt::pow(y,3)) - 1008*(xt::pow(x,3))*(xt::pow(y,5)) - 1008*x*(xt::pow(y,7)) - 480*(xt::pow(x,9))*y - 960*(xt::pow(x,7))*(xt::pow(y,3)) + 960*(xt::pow(x,3))*(xt::pow(y,7)) + 480 *x*(xt::pow(y,9)));//#m = -4   n = 10
        xt::xtensor<double,2> Z60 =  c( 59) * (30*x*y - 280*(xt::pow(x,3))*y - 280*x*(xt::pow(y,3)) + 840*(xt::pow(x,5))*y + 1680*(xt::pow(x,3))*(xt::pow(y,3)) +840*x*(xt::pow(y,5)) - 1008*(xt::pow(x,7))*y - 3024*(xt::pow(x,5))*(xt::pow(y,3)) - 3024*(xt::pow(x,3))*(xt::pow(y,5)) - 1008*x*(xt::pow(y,7)) + 420*(xt::pow(x,9))*y + 1680*(xt::pow(x,7))*(xt::pow(y,3)) + 2520*(xt::pow(x,5))*(xt::pow(y,5)) + 1680*(xt::pow(x,3))*(xt::pow(y,7)) + 420*x*(xt::pow(y,9)) );//#m = -2   n = 10
        xt::xtensor<double,2> Z61 =  c( 60) * (-1 + 30*xt::pow(x,2) + 30*xt::pow(y,2) - 210*xt::pow(x,4) - 420*(xt::pow(x,2))*(xt::pow(y,2)) - 210*xt::pow(y,4) + 560*xt::pow(x,6) + 1680*(xt::pow(x,4))*(xt::pow(y,2)) + 1680*(xt::pow(x,2))*(xt::pow(y,4)) + 560*xt::pow(y,6) - 630*xt::pow(x,8) - 2520*(xt::pow(x,6))*(xt::pow(y,2)) - 3780*(xt::pow(x,4))*(xt::pow(y,4)) - 2520*(xt::pow(x,2))*(xt::pow(y,6)) - 630*xt::pow(y,8) + 252*xt::pow(x,10) + 1260*(xt::pow(x,8))*(xt::pow(y,2)) + 2520*(xt::pow(x,6))*(xt::pow(y,4)) + 2520*(xt::pow(x,4))*(xt::pow(y,6)) + 1260*(xt::pow(x,2))*(xt::pow(y,8)) + 252*xt::pow(y,10));//#m = 0    n = 10
        xt::xtensor<double,2> Z62 =  c( 61) * (-15*xt::pow(x,2) + 15*xt::pow(y,2) + 140*xt::pow(x,4) - 140*xt::pow(y,4) - 420*xt::pow(x,6) - 420*(xt::pow(x,4))*(xt::pow(y,2)) + 420*(xt::pow(x,2))*(xt::pow(y,4)) + 420*xt::pow(y,6) + 504*xt::pow(x,8) + 1008*(xt::pow(x,6))*(xt::pow(y,2)) - 1008*(xt::pow(x,2))*(xt::pow(y,6)) - 504*xt::pow(y,8) - 210*xt::pow(x,10) - 630*(xt::pow(x,8))*(xt::pow(y,2)) - 420*(xt::pow(x,6))*(xt::pow(y,4)) + 420*(xt::pow(x,4))*(xt::pow(y,6)) + 630*(xt::pow(x,2))*(xt::pow(y,8)) + 210*xt::pow(y,10));//#m = 2  n = 10
        xt::xtensor<double,2> Z63 =  c( 62) * (-35*xt::pow(x,4) + 210*(xt::pow(x,2))*(xt::pow(y,2)) - 35*xt::pow(y,4) + 168*xt::pow(x,6) - 840*(xt::pow(x,4))*(xt::pow(y,2)) - 840*(xt::pow(x,2))*(xt::pow(y,4)) + 168*xt::pow(y,6) - 252*xt::pow(x,8) + 1008*(xt::pow(x,6))*(xt::pow(y,2)) + 2520*(xt::pow(x,4))*(xt::pow(y,4)) + 1008*(xt::pow(x,2))*(xt::pow(y,6)) - 252*(xt::pow(y,8)) + 120*xt::pow(x,10) - 360*(xt::pow(x,8))*(xt::pow(y,2)) - 1680*(xt::pow(x,6))*(xt::pow(y,4)) - 1680*(xt::pow(x,4))*(xt::pow(y,6)) - 360*(xt::pow(x,2))*(xt::pow(y,8)) + 120*xt::pow(y,10));//#m = 4     n = 10
        xt::xtensor<double,2> Z64 =  c( 63) * (-28*xt::pow(x,6) + 420*(xt::pow(x,4))*(xt::pow(y,2)) - 420*(xt::pow(x,2))*(xt::pow(y,4)) + 28*xt::pow(y,6) + 72*xt::pow(x,8) - 1008*(xt::pow(x,6))*(xt::pow(y,2)) + 1008*(xt::pow(x,2))*(xt::pow(y,6)) - 72*xt::pow(y,8) - 45*xt::pow(x,10) + 585*(xt::pow(x,8))*(xt::pow(y,2)) + 630*(xt::pow(x,6))*(xt::pow(y,4)) - 630*(xt::pow(x,4))*(xt::pow(y,6)) - 585*(xt::pow(x,2))*(xt::pow(y,8)) + 45*xt::pow(y,10));//#m = 6    n = 10
        xt::xtensor<double,2> Z65 =  c( 64) * (-9*xt::pow(x,8) + 252*(xt::pow(x,6))*(xt::pow(y,2)) - 630*(xt::pow(x,4))*(xt::pow(y,4)) + 252*(xt::pow(x,2))*(xt::pow(y,6)) - 9*xt::pow(y,8) + 10*xt::pow(x,10) - 270*(xt::pow(x,8))*(xt::pow(y,2)) + 420*(xt::pow(x,6))*(xt::pow(y,4)) + 420*(xt::pow(x,4))*(xt::pow(y,6)) - 270*(xt::pow(x,2))*(xt::pow(y,8)) + 10*xt::pow(y,10));//#m = 8    n = 10
        xt::xtensor<double,2> Z66 =  c( 65) * (-1*xt::pow(x,10) + 45*(xt::pow(x,8))*(xt::pow(y,2)) - 210*(xt::pow(x,6))*(xt::pow(y,4)) + 210*(xt::pow(x,4))*(xt::pow(y,6)) - 45*(xt::pow(x,2))*(xt::pow(y,8)) + xt::pow(y,10));//#m = 10   n = 10
        xt::xtensor<double,2> ZW =    Z1 + Z2 +  Z3+  Z4+  Z5+  Z6+  Z7+  Z8+  Z9+  Z10+ Z11+ Z12+ Z13+ Z14+ Z15+ Z16+ Z17+ Z18+ Z19+ Z20+ Z21+ Z22+ Z23+ Z24+ Z25+ Z26+ Z27+ Z28+ Z29+ Z30+ Z31+ Z32+ Z33+ Z34+ Z35+ Z36+ Z37+ Z38+ Z39+ Z40+ Z41+ Z42+ Z43+ Z44+ Z45+ Z46+ Z47+ Z48+ Z49+Z50+ Z51+ Z52+ Z53+ Z54+ Z55+ Z56+ Z57+ Z58+ Z59+Z60+ Z61+ Z62+ Z63+ Z64+ Z65+ Z66 ;

        return ZW ;

    }
