//# ZernikePolyCalc.cc: Implementation for ZernikeCalc
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

#include <math.h>
#include "ZernikePolyCalc.h"
#include <sys/stat.h>

using namespace std;


  double ZernikePolyCalc::powl(double base, int exp)
  {
    /* Fast power function */

    if (exp == 0)
      return 1;
    else if (exp == 1)
      return base;
    else if ((exp & 1) != 0)
      return base * powl(base*base, floor(exp/2));
    else
      return pow(base*base, floor(exp/2));
  }

  vector< vector <double> > ZernikePolyCalc::zernikeSurface(vector<double>& amp, vector<vector<float> >& xCoords, 
				      vector<vector<float> >& yCoords,int xSize, int ySize, 
				      vector<vector<double> >& surface)
  {

    Z_p.resize(66,0);
    C_p.resize(66,0);
	
    // Copy input ampR and ampI into cR and cI. In case all 64 terms are not required.
    copy(amp.begin(), amp.end(), C_p.begin());
    
    // summon cthulhu
    for(int ii = 0; ii < xSize; ii++) {
      for(int jj = 0; jj < ySize; jj++) {
	float x = xCoords[ii][jj];
	float y = yCoords[ii][jj];
			
    if (sqrt((pow((x),2)+pow((y),2))) < 1.0){
        Z_p[0] = C_p[0];

        Z_p[1] = C_p[1] * x;

        Z_p[2] = C_p[2] * y;

        Z_p[3] = C_p[3] * 2 * x * y;

        Z_p[4] = C_p[4] * (2*powl(x,2) + 2*powl(y,2) - 1);

        Z_p[5] = C_p[5] * (-1*powl(x,2) + powl(y,2));

        Z_p[6] = C_p[6] * (-1*powl(x,3) + 3*x*powl(y,2));

        Z_p[7] = C_p[7] * (-2*x + 3*powl(x,3) + 3*x*powl(y,2));

        Z_p[8] = C_p[8] * (-2*y + 3*powl(y,3) + 3*powl(x,2)*y);

        Z_p[9] = C_p[9] * (powl(y,3) - 3*powl(x,2)*y);

        Z_p[10] = C_p[10] * (-4*powl(x,3)*y + 4*x*powl(y,3));

        Z_p[11] = C_p[11] * (-6*x*y + 8*powl(x,3)*y + 8*x*powl(y,3));

        Z_p[12] = C_p[12] * (1-6*powl(x,2) -6*powl(y,2) + 6*powl(x,4) + 12*(powl(x,2)*powl(y,2)) + 6*powl(y,4));

        Z_p[13] = C_p[13] * (3*powl(x,2) - 3*powl(y,2) - 4*powl(x,4) + 4*powl(y,4));

        Z_p[14] = C_p[14] * (powl(x,4) - 6*powl(x,2)*powl(y,2) + powl(y,4));

        Z_p[15] = C_p[15] * (powl(x,5)-10*powl(x,3)*powl(y,2) + 5*x*powl(y,4));

        Z_p[16] = C_p[16] * (4*powl(x,3) - 12*x*(powl(y,2)) -5*powl(x,5) + 10*(powl(x,3))*(powl(y,2)) + 15*x*powl(y,4));

        Z_p[17] = C_p[17] * (3*x - 12*powl(x,3) - 12*x*(powl(y,2)) + 10*powl(x,5) + 20*(powl(x,3))*(powl(y,2)) + 10*x*(powl(y,4))); //m= -1  n = 5

        Z_p[18] = C_p[18] * (3*y - 12*powl(y,3) - 12*y*(powl(x,2)) + 10*powl(y,5) + 20*(powl(y,3))*(powl(x,2)) + 10*y*(powl(x,4))); //m = 1  n = 5

        Z_p[19] = C_p[19] * (-4*powl(y,3) + 12*y*(powl(x,2)) + 5*powl(y,5) - 10*(powl(y,3))*(powl(x,2)) - 15*y*powl(x,4)); //m = 3   n = 5

        Z_p[20] =  C_p[20] * (powl(y,5)-10*(powl(y,3))*powl(x,2) + 5*y*(powl(x,4))); //m = 5 n = 5

        Z_p[21] =  C_p[21] * (6*(powl(x,5))*y - 20*(powl(x,3))*(powl(y,3)) + 6*x*(powl(y,5))); //m = -6 n = 6

        Z_p[22] =  C_p[22] * (20*(powl(x,3))*y - 20*x*(powl(y,3)) - 24*(powl(x,5))*y + 24*x*(powl(y,5))); //m = -4   n = 6

        Z_p[23] =  C_p[23] * (12*x*y + 40*(powl(x,3))*y - 40*x*(powl(y,3)) + 30*(powl(x,5))*y + 60*(powl(x,3))*(powl(y,3)) - 30*x*(powl(y,5))); //m = -2   n = 6

        Z_p[24] =  C_p[24] * (-1 + 12*(powl(x,2)) + 12*(powl(y,2)) - 30*(powl(x,4)) - 60*(powl(x,2))*(powl(y,2)) - 30*(powl(y,4)) + 20*(powl(x,6))
                      + 60*(powl(x,4))*powl(y,2) + 60 *(powl(x,2))*(powl(y,4)) + 20*(powl(y,6))); //m = 0   n = 6

        Z_p[25] =  C_p[25] * (-6*(powl(x,2)) + 6*(powl(y,2)) + 20*(powl(x,4)) - 20*(powl(y,4)) - 15*(powl(x,6)) - 15*(powl(x,4))*(powl(y,2))
                      + 15*(powl(x,2))*(powl(y,4)) + 15*(powl(y,6))); //m = 2   n = 6

        Z_p[26] =  C_p[26] * (-5*(powl(x,4)) + 30*(powl(x,2))*(powl(y,2)) - 5*(powl(y,4)) + 6*(powl(x,6)) - 30*(powl(x,4))*powl(y,2)
                      - 30*(powl(x,2))*(powl(y,4)) + 6*(powl(y,6))); //m = 4    n = 6

        Z_p[27] =  C_p[27] * (-1*(powl(x,6)) + 15*(powl(x,4))*(powl(y,2)) - 15*(powl(x,2))*(powl(y,4)) + powl(y,6)); //m = 6   n = 6

        Z_p[28] =  C_p[28] * (-1*(powl(x,7)) + 21*(powl(x,5))*(powl(y,2)) - 35*(powl(x,3))*(powl(y,4)) + 7*x*(powl(y,6))); //m = -7    n = 7

        Z_p[29] =  C_p[29] * (-6*(powl(x,5)) + 60*(powl(x,3))*(powl(y,2)) - 30*x*(powl(y,4)) + 7*powl(x,7) - 63*(powl(x,5))*(powl(y,2))
                      - 35*(powl(x,3))*(powl(y,4)) + 35*x*(powl(y,6))) ; //m = -5    n = 7

        Z_p[30] =  C_p[30] * (-10*(powl(x,3)) + 30*x*(powl(y,2)) + 30*powl(x,5) - 60*(powl(x,3))*(powl(y,2)) - 90*x*(powl(y,4)) - 21*powl(x,7)
                      + 21*(powl(x,5))*(powl(y,2)) + 105*(powl(x,3))*(powl(y,4)) + 63*x*(powl(y,6))); //m =-3       n = 7

        Z_p[31] =  C_p[31] * (-4*x + 30*powl(x,3) + 30*x*(powl(y,2)) - 60*(powl(x,5)) - 120*(powl(x,3))*(powl(y,2)) - 60*x*(powl(y,4))
                      + 35*powl(x,7) + 105*(powl(x,5))*(powl(y,2)) + 105*(powl(x,3))*(powl(y,4)) + 35*x*(powl(y,6))); //m = -1  n = 7

        Z_p[32] =  C_p[32] * (-4*y + 30*powl(y,3) + 30*y*(powl(x,2)) - 60*(powl(y,5)) - 120*(powl(y,3))*(powl(x,2)) - 60*y*(powl(x,4))
                      + 35*powl(y,7) + 105*(powl(y,5))*(powl(x,2)) + 105*(powl(y,3))*(powl(x,4)) + 35*y*(powl(x,6))); //m = 1   n = 7

        Z_p[33] =  C_p[33] * (10*(powl(y,3)) - 30*y*(powl(x,2)) - 30*powl(y,5) + 60*(powl(y,3))*(powl(x,2)) + 90*y*(powl(x,4)) + 21*powl(y,7)
                      - 21*(powl(y,5))*(powl(x,2)) - 105*(powl(y,3))*(powl(x,4)) - 63*y*(powl(x,6))); //m =3     n = 7

        Z_p[34] =  C_p[34] * (-6*(powl(y,5)) + 60*(powl(y,3))*(powl(x,2)) - 30*y*(powl(x,4)) + 7*powl(y,7) - 63*(powl(y,5))*(powl(x,2))
                      - 35*(powl(y,3))*(powl(x,4)) + 35*y*(powl(x,6))); //m = 5  n = 7

        Z_p[35] =  C_p[35] * (powl(y,7) - 21*(powl(y,5))*(powl(x,2)) + 35*(powl(y,3))*(powl(x,4)) - 7*y*(powl(x,6))); //m = 7  n = 7

        Z_p[36] =  C_p[36] * (-8*(powl(x,7))*y + 56*(powl(x,5))*(powl(y,3)) - 56*(powl(x,3))*(powl(y,5)) + 8*x*(powl(y,7))); //m = -8  n = 8

        Z_p[37] =  C_p[37] * (-42*(powl(x,5))*y + 140*(powl(x,3))*(powl(y,3)) - 42*x*(powl(y,5)) + 48*(powl(x,7))*y - 112*(powl(x,5))*(powl(y,3))
                      - 112*(powl(x,3))*(powl(y,5)) + 48*x*(powl(y,7))); //m = -6  n = 8

        Z_p[38] =  C_p[38] * (-60*(powl(x,3))*y + 60*x*(powl(y,3)) + 168*(powl(x,5))*y -168*x*(powl(y,5)) - 112*(powl(x,7))*y
                      - 112*(powl(x,5))*(powl(y,3)) + 112*(powl(x,3))*(powl(y,5)) + 112*x*(powl(y,7))); //m = -4   n = 8

        Z_p[39] =  C_p[39] * (-20*x*y + 120*(powl(x,3))*y + 120*x*(powl(y,3)) - 210*(powl(x,5))*y - 420*(powl(x,3))*(powl(y,3)) - 210*x*(powl(y,5))
                      - 112*(powl(x,7))*y + 336*(powl(x,5))*(powl(y,3)) + 336*(powl(x,3))*(powl(y,5)) + 112*x*(powl(y,7))); //m = -2   n = 8

        Z_p[40] =  C_p[40] * (1 - 20*powl(x,2) - 20*powl(y,2) + 90*powl(x,4) + 180*(powl(x,2))*(powl(y,2)) + 90*powl(y,4) - 140*powl(x,6)
                      - 420*(powl(x,4))*(powl(y,2)) - 420*(powl(x,2))*(powl(y,4)) - 140*(powl(y,6)) + 70*powl(x,8)
                      + 280*(powl(x,6))*(powl(y,2)) + 420*(powl(x,4))*(powl(y,4)) + 280*(powl(x,2))*(powl(y,6)) + 70*powl(y,8)); //m = 0    n = 8

        Z_p[41] =  C_p[41] * (10*powl(x,2) - 10*powl(y,2) - 60*powl(x,4) + 105*(powl(x,4))*(powl(y,2)) - 105*(powl(x,2))*(powl(y,4))
                      + 60*powl(y,4) + 105*powl(x,6) - 105*powl(y,6) - 56*powl(x,8) - 112*(powl(x,6))*(powl(y,2))
                      + 112*(powl(x,2))*(powl(y,6)) + 56*powl(y,8)); //m = 2  n = 8

        Z_p[42] =  C_p[42] * (15*powl(x,4) - 90*(powl(x,2))*(powl(y,2)) + 15*powl(y,4) - 42*powl(x,6) + 210*(powl(x,4))*(powl(y,2))
                      + 210*(powl(x,2))*(powl(y,4)) - 42*powl(y,6) + 28*powl(x,8) - 112*(powl(x,6))*(powl(y,2))
                      - 280*(powl(x,4))*(powl(y,4)) - 112*(powl(x,2))*(powl(y,6)) + 28*powl(y,8)); //m = 4     n = 8

        Z_p[43] =  C_p[43] * (7*powl(x,6) - 105*(powl(x,4))*(powl(y,2)) + 105*(powl(x,2))*(powl(y,4)) - 7*powl(y,6) - 8*powl(x,8)
                      + 112*(powl(x,6))*(powl(y,2)) - 112*(powl(x,2))*(powl(y,6)) + 8*powl(y,8)); //m = 6    n = 8

        Z_p[44] =  C_p[44] * (powl(x,8) - 28*(powl(x,6))*(powl(y,2)) + 70*(powl(x,4))*(powl(y,4)) - 28*(powl(x,2))*(powl(y,6)) + powl(y,8)); //m = 8     n = 9

        Z_p[45] =  C_p[45] * (powl(x,9) - 36*(powl(x,7))*(powl(y,2)) + 126*(powl(x,5))*(powl(y,4)) - 84*(powl(x,3))*(powl(y,6)) + 9*x*(powl(y,8))); //m = -9     n = 9

        Z_p[46] =  C_p[46] * (8*powl(x,7) - 168*(powl(x,5))*(powl(y,2)) + 280*(powl(x,3))*(powl(y,4)) - 56 *x*(powl(y,6)) - 9*powl(x,9)
                      + 180*(powl(x,7))*(powl(y,2)) - 126*(powl(x,5))*(powl(y,4)) - 252*(powl(x,3))*(powl(y,6)) + 63*x*(powl(y,8))); //m = -7    n = 9

        Z_p[47] =  C_p[47] * (21*powl(x,5) - 210*(powl(x,3))*(powl(y,2)) + 105*x*(powl(y,4)) - 56*powl(x,7) + 504*(powl(x,5))*(powl(y,2))
                      + 280*(powl(x,3))*(powl(y,4)) - 280*x*(powl(y,6)) + 36*powl(x,9) - 288*(powl(x,7))*(powl(y,2))
                      - 504*(powl(x,5))*(powl(y,4)) + 180*x*(powl(y,8))); //m = -5    n = 9

        Z_p[48] =  C_p[48] * (20*powl(x,3) - 60*x*(powl(y,2)) - 105*powl(x,5) + 210*(powl(x,3))*(powl(y,2)) + 315*x*(powl(y,4)) + 168*powl(x,7)
                      - 168*(powl(x,5))*(powl(y,2)) - 840*(powl(x,3))*(powl(y,4)) - 504*x*(powl(y,6)) - 84*powl(x,9)
                      + 504*(powl(x,5))*(powl(y,4)) + 672*(powl(x,3))*(powl(y,6)) + 252*x*(powl(y,8))); //m = -3  n = 9

        Z_p[49] =  C_p[49] * (5*x - 60*powl(x,3) - 60*x*(powl(y,2)) + 210*powl(x,5) + 420*(powl(x,3))*(powl(y,2)) + 210*x*(powl(y,4))
                      - 280*powl(x,7) - 840*(powl(x,5))*(powl(y,2)) - 840*(powl(x,3))*(powl(y,4)) - 280*x*(powl(y,6)) + 126*powl(x,9)
                      + 504*(powl(x,7))*(powl(y,2)) + 756*(powl(x,5))*(powl(y,4)) + 504*(powl(x,3))*(powl(y,6)) + 126*x*(powl(y,8))); //m = -1   n = 9

        Z_p[50] =  C_p[50] * (5*y - 60*powl(y,3) - 60*y*(powl(x,2)) + 210*powl(y,5) + 420*(powl(y,3))*(powl(x,2)) + 210*y*(powl(x,4))
                      - 280*powl(y,7) - 840*(powl(y,5))*(powl(x,2)) - 840*(powl(y,3))*(powl(x,4)) - 280*y*(powl(x,6)) + 126*powl(y,9)
                      + 504*(powl(y,7))*(powl(x,2)) + 756*(powl(y,5))*(powl(x,4)) + 504*(powl(y,3))*(powl(x,6)) + 126*y*(powl(x,8))); //m = -1   n = 9

        Z_p[51] =  C_p[51] * (-20*powl(y,3) + 60*y*(powl(x,2)) + 105*powl(y,5) - 210*(powl(y,3))*(powl(x,2)) - 315*y*(powl(x,4)) - 168*powl(y,7)
                      + 168*(powl(y,5))*(powl(x,2)) + 840*(powl(y,3))*(powl(x,4)) + 504*y*(powl(x,6)) + 84*powl(y,9)
                      - 504*(powl(y,5))*(powl(x,4)) - 672*(powl(y,3))*(powl(x,6)) - 252*y*(powl(x,8))); //m = 3  n = 9

        Z_p[52] =  C_p[52] * (21*powl(y,5) - 210*(powl(y,3))*(powl(x,2)) + 105*y*(powl(x,4)) - 56*powl(y,7) + 504*(powl(y,5))*(powl(x,2))
                      + 280*(powl(y,3))*(powl(x,4)) - 280*y*(powl(x,6)) + 36*powl(y,9) - 288*(powl(y,7))*(powl(x,2))
                      - 504*(powl(y,5))*(powl(x,4)) + 180*y*(powl(x,8))); //m = 5     n = 9

        Z_p[53] =  C_p[53] *(-8*powl(y,7) + 168*(powl(y,5))*(powl(x,2)) - 280*(powl(y,3))*(powl(x,4)) + 56 *y*(powl(x,6)) + 9*powl(y,9)
                     - 180*(powl(y,7))*(powl(x,2)) + 126*(powl(y,5))*(powl(x,4)) - 252*(powl(y,3))*(powl(x,6)) - 63*y*(powl(x,8))); //m = 7     n = 9

        Z_p[54] =  C_p[54] *(powl(y,9) - 36*(powl(y,7))*(powl(x,2)) + 126*(powl(y,5))*(powl(x,4)) - 84*(powl(y,3))*(powl(x,6)) + 9*y*(powl(x,8))); //m = 9       n = 9

        Z_p[55] =  C_p[55] *(10*(powl(x,9))*y - 120*(powl(x,7))*(powl(y,3)) + 252*(powl(x,5))*(powl(y,5)) - 120*(powl(x,3))*(powl(y,7)) + 10*x*(powl(y,9))); //m = -10   n = 10

        Z_p[56] =  C_p[56] *(72*(powl(x,7))*y - 504*(powl(x,5))*(powl(y,3)) + 504*(powl(x,3))*(powl(y,5)) - 72*x*(powl(y,7)) - 80*(powl(x,9))*y
                     + 480*(powl(x,7))*(powl(y,3)) - 480*(powl(x,3))*(powl(y,7)) + 80*x*(powl(y,9))); //m = -8    n = 10

        Z_p[57] =  C_p[57] *(270*(powl(x,9))*y - 360*(powl(x,7))*(powl(y,3)) - 1260*(powl(x,5))*(powl(y,5)) - 360*(powl(x,3))*(powl(y,7))
                     + 270*x*(powl(y,9)) - 432*(powl(x,7))*y + 1008*(powl(x,5))*(powl(y,3)) + 1008*(powl(x,3))*(powl(y,5))
                     - 432*x*(powl(y,7)) + 168*(powl(x,5))*y - 560*(powl(x,3))*(powl(y,3)) + 168*x*(powl(y,5))); //m = -6   n = 10

        Z_p[58] =  C_p[58] *(140*(powl(x,3))*y - 140*x*(powl(y,3)) - 672*(powl(x,5))*y + 672*x*(powl(y,5)) + 1008*(powl(x,7))*y
                     + 1008*(powl(x,5))*(powl(y,3)) - 1008*(powl(x,3))*(powl(y,5)) - 1008*x*(powl(y,7)) - 480*(powl(x,9))*y
                     - 960*(powl(x,7))*(powl(y,3)) + 960*(powl(x,3))*(powl(y,7)) + 480 *x*(powl(y,9))); //m = -4   n = 10

        Z_p[59] =  C_p[59] *(30*x*y - 280*(powl(x,3))*y - 280*x*(powl(y,3)) + 840*(powl(x,5))*y + 1680*(powl(x,3))*(powl(y,3))
                     +840*x*(powl(y,5)) - 1008*(powl(x,7))*y - 3024*(powl(x,5))*(powl(y,3)) - 3024*(powl(x,3))*(powl(y,5))
                     - 1008*x*(powl(y,7)) + 420*(powl(x,9))*y + 1680*(powl(x,7))*(powl(y,3)) + 2520*(powl(x,5))*(powl(y,5))
                     + 1680*(powl(x,3))*(powl(y,7)) + 420*x*(powl(y,9)) ); //m = -2   n = 10

        Z_p[60] =  C_p[60] * (-1 + 30*powl(x,2) + 30*powl(y,2) - 210*powl(x,4) - 420*(powl(x,2))*(powl(y,2)) - 210*powl(y,4)
                      + 560*powl(x,6) + 1680*(powl(x,4))*(powl(y,2)) + 1680*(powl(x,2))*(powl(y,4)) + 560*powl(y,6)
                      - 630*powl(x,8) - 2520*(powl(x,6))*(powl(y,2)) - 3780*(powl(x,4))*(powl(y,4)) - 2520*(powl(x,2))*(powl(y,6))
                      - 630*powl(y,8) + 252*powl(x,10) + 1260*(powl(x,8))*(powl(y,2)) + 2520*(powl(x,6))*(powl(y,4))
                      + 2520*(powl(x,4))*(powl(y,6)) + 1260*(powl(x,2))*(powl(y,8)) + 252*powl(y,10)); //m = 0    n = 10

        Z_p[61] =  C_p[61] * (-15*powl(x,2) + 15*powl(y,2) + 140*powl(x,4) - 140*powl(y,4) - 420*powl(x,6) - 420*(powl(x,4))*(powl(y,2))
                      + 420*(powl(x,2))*(powl(y,4)) + 420*powl(y,6) + 504*powl(x,8) + 1008*(powl(x,6))*(powl(y,2))
                      - 1008*(powl(x,2))*(powl(y,6)) - 504*powl(y,8) - 210*powl(x,10) - 630*(powl(x,8))*(powl(y,2))
                      - 420*(powl(x,6))*(powl(y,4)) + 420*(powl(x,4))*(powl(y,6)) + 630*(powl(x,2))*(powl(y,8)) + 210*powl(y,10)); // m = 2  n = 10

        Z_p[62] =  C_p[62] *(-35*powl(x,4) + 210*(powl(x,2))*(powl(y,2)) - 35*powl(y,4) + 168*powl(x,6) - 840*(powl(x,4))*(powl(y,2))
                     - 840*(powl(x,2))*(powl(y,4)) + 168*powl(y,6) - 252*powl(x,8) + 1008*(powl(x,6))*(powl(y,2))
                     + 2520*(powl(x,4))*(powl(y,4)) + 1008*(powl(x,2))*(powl(y,6)) - 252*(powl(y,8)) + 120*powl(x,10)
                     - 360*(powl(x,8))*(powl(y,2)) - 1680*(powl(x,6))*(powl(y,4)) - 1680*(powl(x,4))*(powl(y,6))
                     - 360*(powl(x,2))*(powl(y,8)) + 120*powl(y,10)); //m = 4     n = 10

        Z_p[63] =  C_p[63] *(-28*powl(x,6) + 420*(powl(x,4))*(powl(y,2)) - 420*(powl(x,2))*(powl(y,4)) + 28*powl(y,6) + 72*powl(x,8)
                     - 1008*(powl(x,6))*(powl(y,2)) + 1008*(powl(x,2))*(powl(y,6)) - 72*powl(y,8) - 45*powl(x,10)
                     + 585*(powl(x,8))*(powl(y,2)) + 630*(powl(x,6))*(powl(y,4)) - 630*(powl(x,4))*(powl(y,6))
                     - 585*(powl(x,2))*(powl(y,8)) + 45*powl(y,10)); //m = 6    n = 10

        Z_p[64] =  C_p[64] *(-9*powl(x,8) + 252*(powl(x,6))*(powl(y,2)) - 630*(powl(x,4))*(powl(y,4)) + 252*(powl(x,2))*(powl(y,6))
                     - 9*powl(y,8) + 10*powl(x,10) - 270*(powl(x,8))*(powl(y,2)) + 420*(powl(x,6))*(powl(y,4))
                     + 420*(powl(x,4))*(powl(y,6)) - 270*(powl(x,2))*(powl(y,8)) + 10*powl(y,10)); //m = 8    n = 10

        Z_p[65] =  C_p[65] *(-1*powl(x,10) + 45*(powl(x,8))*(powl(y,2)) - 210*(powl(x,6))*(powl(y,4)) + 210*(powl(x,4))*(powl(y,6))
                     - 45*(powl(x,2))*(powl(y,8)) + powl(y,10)); //m = 10   n = 10


         for(unsigned int kk = 0; kk < Z_p.size(); kk++) {
              surface[ii][jj] += Z_p[kk]; //accumulate(Z_p.begin(),Z_p.end(),0);
             }
        }
      }
    }
		

    return surface;
  }

  
