//# ZernikeCalc.h: Definition for ZernikeCalc
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

#ifndef SYNTHESIS_ZERNIKEPOLYCALC_H
#define SYNTHESIS_ZERNIKEPOLYCALC_H

#include <stdlib.h>
#include <iostream> 
#include <fstream>
#include <vector>
using namespace std;  
  class ZernikePolyCalc
  {
  public:
    ZernikePolyCalc():Z_p(),C_p()
    {};

    ~ZernikePolyCalc(){};

    ZernikePolyCalc& operator=(const ZernikePolyCalc& other);

    double powl(double base, int exp);
    void nullAperture(vector<vector<double>>& surface, int xSize, int ySize);
    vector<vector<double>> zernikeSurface(vector<double>& amp, vector<vector<float>>& xCoords, 
						    vector<vector<float>>& yCoords, int xSize, int ySize,
						    vector<vector<double>>& surface);

    inline int ncoeffs() {return C_p.size();};
    inline vector<float> getCoeffs() {return C_p;};
    inline vector<float> getZPolyvals() {return C_p;};
    
    vector<float> Z_p, C_p; // The _p is meant to represent the private variables in this class.

  };
#endif
