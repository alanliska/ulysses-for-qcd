/*ULYSSES, a semi-empirical package
    Copyright (C) 2023- Filipe Menezes (filipe.menezes@helmholtz-munich.de)
                        Grzegorz Popowicz (grzegorz.popowicz@helmholtz-munich.de)

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software 
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA*/

#ifndef _Solvation_Grid_
#define _Solvation_Grid_
#include <vector>
#include <string>
#include "../UnitConversion.hpp"
#include "LebedevLaikovGridData.hpp"

//descrition:
//parameters to auxiliate the generation of grids for the calculation of charges of atoms in molecules

std::vector<int> LebedevGrid() {
  //function returning Lebedev grid
  std::vector<int> lebedev(32,0);
  lebedev[0] = 6;
  lebedev[1] = 14;
  lebedev[2] = 26;
  lebedev[3] = 38;
  lebedev[4] = 50;
  lebedev[5] = 74;
  lebedev[6] = 86;
  lebedev[7] = 110;
  lebedev[8] = 146;
  lebedev[9] = 170;
  lebedev[10] = 194;
  lebedev[11] = 230;
  lebedev[12] = 266;
  lebedev[13] = 302;
  lebedev[14] = 350;
  lebedev[15] = 434;
  lebedev[16] = 590;
  lebedev[17] = 770;
  lebedev[18] = 974;
  lebedev[19] = 1202;
  lebedev[20] = 1454;
  lebedev[21] = 1730;
  lebedev[22] = 2030;
  lebedev[23] = 2354;
  lebedev[24] = 2702;
  lebedev[25] = 3074;
  lebedev[26] = 3470;
  lebedev[27] = 3890;
  lebedev[28] = 4334;
  lebedev[29] = 4802;
  lebedev[30] = 5294;
  lebedev[31] = 5810;
  return lebedev;
}
void getLebedevLaikovGrid(int pos, std::vector<double> & Weights, matrixE & Grid) {
  //function to get the Lebedev Laikov grid to integrate Born radii
  int gridsz = Weights.size();
  for (size_t idx = 0; idx < gridsz; ++idx) {
    Weights[idx] = 0.0;
    Grid(1,idx + 1) = 0.0;
    Grid(2,idx + 1) = 0.0;
    Grid(3,idx + 1) = 0.0;
  }
  switch (pos) {
    case 1:
      ld0006(Weights,Grid);
      break;
    case 2:
      ld0014(Weights,Grid);
      break;
    case 3:
      ld0026(Weights,Grid);
      break;
    case 4:
      ld0038(Weights,Grid);
      break;
    case 5:
      ld0050(Weights,Grid);
      break;
    case 6:
      ld0074(Weights,Grid);
      break;
    case 7:
      ld0086(Weights,Grid);
      break;
    case 8:
      ld0110(Weights,Grid);
      break;
    case 9:
      ld0146(Weights,Grid);
      break;
    case 10:
      ld0170(Weights,Grid);
      break;
    case 11:
      ld0194(Weights,Grid);
      break;
    case 12:
      ld0230(Weights,Grid);
      break;
    case 13:
      ld0266(Weights,Grid);
      break;
    case 14:
      ld0302(Weights,Grid);
      break;
    case 15:
      ld0350(Weights,Grid);
      break;
    case 16:
      ld0434(Weights,Grid);
      break;
    case 17:
      ld0590(Weights,Grid);
      break;
    case 18:
      ld0770(Weights,Grid);
      break;
    case 19:
      ld0974(Weights,Grid);
      break;
    case 20:
      ld1202(Weights,Grid);
      break;
    case 21:
      ld1454(Weights,Grid);
      break;
    case 22:
      ld1730(Weights,Grid);
      break;
    case 23:
      ld2030(Weights,Grid);
      break;
    case 24:
      ld2354(Weights,Grid);
      break;
    case 25:
      ld2702(Weights,Grid);
      break;
    case 26:
      ld3074(Weights,Grid);
      break;
    case 27:
      ld3470(Weights,Grid);
      break;
    case 28:
      ld3890(Weights,Grid);
      break;
    case 29:
      ld4334(Weights,Grid);
      break;
    case 30:
      ld4802(Weights,Grid);
      break;
    case 31:
      ld5294(Weights,Grid);
      break;
    case 32:
      ld5810(Weights,Grid);
      break;
  }
}


#endif //_Solvation_Grid_
