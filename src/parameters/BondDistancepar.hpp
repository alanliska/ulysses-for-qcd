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

#ifndef _Bond_Distance_Parameters_
#define _Bond_Distance_Parameters_

//description:
//experimental bond distances

double ExperimentalBondDistance(int atom1, int atom2, double bondorder) {
  //collection of experimental bond distances from CCCBDB, as collected from the comparison with B3LYP calculations
  //data averaged
  int atmA = std::min(atom1,atom2);
  int atmB = std::max(atom1,atom2);
  double rdist = 0.0;
  switch (atmA) {
    case 1:
      if ((atmB == 3)&&(bondorder == 1.0)) {rdist = 1.595;}
      else if ((atmB == 4)&&(bondorder == 1.0)) {rdist = 1.335;}
      else if ((atmB == 5)&&(bondorder == 1.0)) {rdist = 1.238;}
      else if ((atmB == 6)&&(bondorder == 1.0)) {rdist = 1.092;}
      else if ((atmB == 7)&&(bondorder == 1.0)) {rdist = 1.022;}
      else if ((atmB == 8)&&(bondorder == 1.0)) {rdist = 0.964;}
      else if ((atmB == 11)&&(bondorder == 1.0)) {rdist = 1.887;}
      else if ((atmB == 12)&&(bondorder == 1.0)) {rdist = 1.694;}
      else if ((atmB == 13)&&(bondorder == 1.0)) {rdist = 1.645;}
      else if ((atmB == 14)&&(bondorder == 1.0)) {rdist = 1.482;}
      else if ((atmB == 15)&&(bondorder == 1.0)) {rdist = 1.418;}
      else if ((atmB == 16)&&(bondorder == 1.0)) {rdist = 1.351;}
      else if ((atmB == 19)&&(bondorder == 1.0)) {rdist = 2.243;}
      else if ((atmB == 20)&&(bondorder == 1.0)) {rdist = 2.003;}
      else if ((atmB == 21)&&(bondorder == 1.0)) {rdist = 1.775;}
      else if ((atmB == 24)&&(bondorder == 1.0)) {rdist = 1.655;}
      else if ((atmB == 29)&&(bondorder == 1.0)) {rdist = 1.463;}
      else if ((atmB == 30)&&(bondorder == 1.0)) {rdist = 1.515;}
      else if ((atmB == 32)&&(bondorder == 1.0)) {rdist = 1.588;}
      else if ((atmB == 33)&&(bondorder == 1.0)) {rdist = 1.524;}
      else if ((atmB == 34)&&(bondorder == 1.0)) {rdist = 1.469;}
      break;
    case 3:
      if ((atmB == 3)&&(bondorder == 1.0)) {rdist = 2.673;}
      else if ((atmB == 7)&&(bondorder == 1.0)) {rdist = 1.736;}
      else if ((atmB == 9)&&(bondorder == 1.0)) {rdist = 1.564;}
      else if ((atmB == 11)&&(bondorder == 1.0)) {rdist = 2.889;}
      else if ((atmB == 17)&&(bondorder == 1.0)) {rdist = 2.029;}
      else if ((atmB == 35)&&(bondorder == 1.0)) {rdist = 2.170;}
      else if ((atmB == 53)&&(bondorder == 1.0)) {rdist = 2.392;}
      break;
    case 4:
      if ((atmB == 8)&&(bondorder == 2.0)) {rdist = 1.331;}
      else if ((atmB == 9)&&(bondorder == 1.0)) {rdist = 1.361;}
      else if ((atmB == 16)&&(bondorder == 2.0)) {rdist = 1.742;}
      else if ((atmB == 17)&&(bondorder == 1.0)) {rdist = 1.797;}
      break;
    case 5:
      if ((atmB == 5)&&(bondorder == 1.0)) {rdist = 1.590;}
      else if ((atmB == 7)&&(bondorder == 2.0)) {rdist = 1.325;}
      else if ((atmB == 8)&&(bondorder == 2.0)) {rdist = 1.205;}
      else if ((atmB == 9)&&(bondorder == 1.0)) {rdist = 1.307;}
      else if ((atmB == 16)&&(bondorder == 1.0)) {rdist = 1.609;}
      else if ((atmB == 17)&&(bondorder == 1.0)) {rdist = 1.719;}
      break;
    case 6:
      if ((atmB == 6)&&(bondorder == 1.0)) {rdist = 1.514;}
      else if ((atmB == 6)&&(bondorder == 1.5)) {rdist = 1.400;}
      else if ((atmB == 6)&&(bondorder == 2.0)) {rdist = 1.306;}
      else if ((atmB == 6)&&(bondorder == 3.0)) {rdist = 1.222;}
      else if ((atmB == 7)&&(bondorder == 1.0)) {rdist = 1.410;}
      else if ((atmB == 7)&&(bondorder == 1.5)) {rdist = 1.339;}
      else if ((atmB == 7)&&(bondorder == 2.0)) {rdist = 1.282;}
      else if ((atmB == 7)&&(bondorder == 3.0)) {rdist = 1.169;}
      else if ((atmB == 8)&&(bondorder == 1.0)) {rdist = 1.399;}
      else if ((atmB == 8)&&(bondorder == 1.5)) {rdist = 1.362;}
      else if ((atmB == 8)&&(bondorder == 2.0)) {rdist = 1.200;}
      else if ((atmB == 8)&&(bondorder == 3.0)) {rdist = 1.128;}
      else if ((atmB == 9)&&(bondorder == 1.0)) {rdist = 1.339;}
      else if ((atmB == 13)&&(bondorder == 2.0)) {rdist = 1.955;}
      else if ((atmB == 14)&&(bondorder == 1.0)) {rdist = 1.799;}
      else if ((atmB == 15)&&(bondorder == 3.0)) {rdist = 1.562;}
      else if ((atmB == 16)&&(bondorder == 1.0)) {rdist = 1.794;}
      else if ((atmB == 16)&&(bondorder == 1.5)) {rdist = 1.785;}
      else if ((atmB == 16)&&(bondorder == 2.0)) {rdist = 1.603;}
      else if ((atmB == 16)&&(bondorder == 3.0)) {rdist = 1.535;}
      else if ((atmB == 17)&&(bondorder == 1.0)) {rdist = 1.716;}
      else if ((atmB == 20)&&(bondorder == 3.0)) {rdist = 2.302;}
      else if ((atmB == 34)&&(bondorder == 1.0)) {rdist = 1.980;}
      else if ((atmB == 34)&&(bondorder == 2.0)) {rdist = 1.694;}
      else if ((atmB == 35)&&(bondorder == 1.0)) {rdist = 1.932;}
      else if ((atmB == 53)&&(bondorder == 1.0)) {rdist = 2.089;}
      break;
    case 7:
      if ((atmB == 7)&&(bondorder == 1.0)) {rdist = 1.181;}
      else if ((atmB == 7)&&(bondorder == 2.0)) {rdist = 1.139;}
      else if ((atmB == 7)&&(bondorder == 3.0)) {rdist = 1.098;}
      else if ((atmB == 8)&&(bondorder == 1.0)) {rdist = 1.395;}
      else if ((atmB == 8)&&(bondorder == 2.0)) {rdist = 1.175;}
      else if ((atmB == 9)&&(bondorder == 1.0)) {rdist = 1.403;}
      else if ((atmB == 14)&&(bondorder == 1.0)) {rdist = 1.575;}
      else if ((atmB == 15)&&(bondorder == 3.0)) {rdist = 1.491;}
      else if ((atmB == 16)&&(bondorder == 1.0)) {rdist = 1.470;}
      else if ((atmB == 17)&&(bondorder == 1.0)) {rdist = 1.611;}
      else if ((atmB == 34)&&(bondorder == 2.0)) {rdist = 1.652;}
      break;
    case 8:
      if ((atmB == 8)&&(bondorder == 1.0)) {rdist = 1.313;}
      else if ((atmB == 9)&&(bondorder == 1.0)) {rdist = 1.377;}
      else if ((atmB == 12)&&(bondorder == 1.0)) {rdist = 1.758;}
      else if ((atmB == 13)&&(bondorder == 1.0)) {rdist = 1.618;}
      else if ((atmB == 14)&&(bondorder == 1.0)) {rdist = 1.510;}
      else if ((atmB == 15)&&(bondorder == 1.0)) {rdist = 1.540;}
      else if ((atmB == 15)&&(bondorder == 2.0)) {rdist = 1.469;}
      else if ((atmB == 16)&&(bondorder == 1.0)) {rdist = 1.574;}
      else if ((atmB == 16)&&(bondorder == 2.0)) {rdist = 1.436;}
      else if ((atmB == 17)&&(bondorder == 1.0)) {rdist = 1.642;}
      else if ((atmB == 17)&&(bondorder == 2.0)) {rdist = 1.404;}
      else if ((atmB == 20)&&(bondorder == 1.0)) {rdist = 1.880;}
      else if ((atmB == 26)&&(bondorder == 1.0)) {rdist = 1.950;}
      else if ((atmB == 26)&&(bondorder == 2.0)) {rdist = 1.650;}
      else if ((atmB == 29)&&(bondorder == 1.0)) {rdist = 1.724;}
      else if ((atmB == 34)&&(bondorder == 1.0)) {rdist = 1.623;}
      else if ((atmB == 35)&&(bondorder == 1.0)) {rdist = 1.776;}
      else if ((atmB == 53)&&(bondorder == 1.0)) {rdist = 1.868;}
      break;
    case 9:
      if ((atmB == 11)&&(bondorder == 1.0)) {rdist = 1.926;}
      else if ((atmB == 12)&&(bondorder == 1.0)) {rdist = 1.750;}
      else if ((atmB == 13)&&(bondorder == 1.0)) {rdist = 1.654;}
      else if ((atmB == 14)&&(bondorder == 1.0)) {rdist = 1.588;}
      else if ((atmB == 15)&&(bondorder == 1.0)) {rdist = 1.556;}
      else if ((atmB == 16)&&(bondorder == 1.0)) {rdist = 1.585;}
      else if ((atmB == 19)&&(bondorder == 1.0)) {rdist = 2.171;}
      else if ((atmB == 20)&&(bondorder == 1.0)) {rdist = 1.967;}
      else if ((atmB == 29)&&(bondorder == 1.0)) {rdist = 1.745;}
      else if ((atmB == 31)&&(bondorder == 1.0)) {rdist = 1.774;}
      else if ((atmB == 32)&&(bondorder == 1.0)) {rdist = 1.705;}
      else if ((atmB == 33)&&(bondorder == 1.0)) {rdist = 1.736;}
      break;
    case 11:
      if ((atmB == 17)&&(bondorder == 1.0)) {rdist = 2.361;}
      else if ((atmB == 35)&&(bondorder == 1.0)) {rdist = 2.502;}
      break;
    case 12:
      if ((atmB == 16)&&(bondorder == 1.0)) {rdist = 2.143;}
      else if ((atmB == 17)&&(bondorder == 1.0)) {rdist = 2.199;}
      break;
    case 13:
      if ((atmB == 16)&&(bondorder == 1.0)) {rdist = 2.029;}
      else if ((atmB == 17)&&(bondorder == 1.0)) {rdist = 2.088;}
      else if ((atmB == 35)&&(bondorder == 1.0)) {rdist = 2.295;}
      break;
    case 14:
      if ((atmB == 14)&&(bondorder == 1.0)) {rdist = 2.246;}
      else if ((atmB == 15)&&(bondorder == 1.0)) {rdist = 2.078;}
      else if ((atmB == 16)&&(bondorder == 1.0)) {rdist = 1.929;}
      else if ((atmB == 17)&&(bondorder == 1.0)) {rdist = 2.040;}
      else if ((atmB == 35)&&(bondorder == 1.0)) {rdist = 2.095;}
      break;
    case 15:
      if ((atmB == 15)&&(bondorder == 1.0)) {rdist = 1.986;}
      else if ((atmB == 16)&&(bondorder == 2.0)) {rdist = 1.900;}
      else if ((atmB == 17)&&(bondorder == 1.0)) {rdist = 2.057;}
      break;
    case 16:
      if ((atmB == 16)&&(bondorder == 1.0)) {rdist = 2.054;}
      else if ((atmB == 17)&&(bondorder == 1.0)) {rdist = 1.993;}
      else if ((atmB == 34)&&(bondorder == 2.0)) {rdist = 2.037;}
      break;
    case 17:
      if ((atmB == 19)&&(bondorder == 1.0)) {rdist = 2.667;}
      else if ((atmB == 20)&&(bondorder == 1.0)) {rdist = 2.437;}
      else if ((atmB == 22)&&(bondorder == 1.0)) {rdist = 2.265;}
      else if ((atmB == 26)&&(bondorder == 1.0)) {rdist = 2.500;}
      else if ((atmB == 29)&&(bondorder == 1.0)) {rdist = 2.051;}
      else if ((atmB == 32)&&(bondorder == 1.0)) {rdist = 2.164;}
      break;
    case 19:
      if ((atmB == 35)&&(bondorder == 1.0)) {rdist = 2.821;}
      break;
    case 20:
      if ((atmB == 35)&&(bondorder == 1.0)) {rdist = 2.594;}
      break;
    case 29:
      if ((atmB == 29)&&(bondorder == 1.0)) {rdist = 2.220;}
      break;
    case 31:
      if ((atmB == 33)&&(bondorder == 3.0)) {rdist = 2.530;}
      break;
    case 33:
      if ((atmB == 33)&&(bondorder == 1.0)) {rdist = 2.435;}
      break;
    case 34:
      if ((atmB == 34)&&(bondorder == 2.0)) {rdist = 2.166;}
      break;
  }
  return rdist;
}

#endif //_Bond_Distance_Parameters_
