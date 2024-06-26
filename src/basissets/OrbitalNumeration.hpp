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

#ifndef _Orbital_Numeration_
#define _Orbital_Numeration_
#include <vector>
#include <string>

//description:
//A bunch of functions related to orbital numeration and their offsets

int OrbitalOffset(int Lorb) {
  //function that determines the position of a given orbital; this is strictly valid and usable with GTOs
  int oo = 0;
  if (Lorb == 1) {oo = 1;}
  else if (Lorb == 2) {oo = 4;}
  else if (Lorb == 3) {oo = 10;}
  return oo;
}
//numeration for the following 3 functions represents the orbital ordering used in GTO bases
//s,px,py,pz,dx2,dy2,dz2,dxy,dxz,dyz,...
//0,1, 2, 3, 4,5,6,  7,  8,  9,  10, ...
int lx(int AOnumber) {
  int LX = 0;
  switch(AOnumber) {
    case 2:
      LX = 1;
      break;
    case 3:
      LX = 0;
      break;
    case 4:
      LX = 0;
      break;
    case 5:
      LX = 2;
      break;
    case 6:
      LX = 0;
      break;
    case 7:
      LX = 0;
      break;
    case 8:
      LX = 1;
      break;
    case 9:
      LX = 1;
      break;
    case 10: 
      LX = 0;
      break;
    case 11: 
      LX = 3;
      break;
    case 12: 
      LX = 0;
      break;
    case 13: 
      LX = 0;
      break;
    case 14: 
      LX = 2;
      break;
    case 15: 
      LX = 2;
      break;
    case 16: 
      LX = 1;
      break;
    case 17: 
      LX = 0;
      break;
    case 18: 
      LX = 1;
      break;
    case 19: 
      LX = 0;
      break;
    case 20: 
      LX = 1;
      break;
    case 21: 
      LX = 4;
      break;
    case 22: 
      LX = 0;
      break;
    case 23: 
      LX = 0;
      break;
    case 24: 
      LX = 3;
      break;
    case 25: 
      LX = 3;
      break;
    case 26: 
      LX = 1;
      break;
    case 27: 
      LX = 0;
      break;
    case 28: 
      LX = 1;
      break;
    case 29: 
      LX = 0;
      break;
    case 30: 
      LX = 2;
      break;
    case 31: 
      LX = 2;
      break;
    case 32: 
      LX = 0;
      break;
    case 33: 
      LX = 2;
      break;
    case 34: 
      LX = 1;
      break;
    case 35: 
      LX = 1;
      break;
    case 36: 
      LX = 5;
      break;
    case 37: 
      LX = 0;
      break;
    case 38: 
      LX = 0;
      break;
    case 39: 
      LX = 3;
      break;
    case 40: 
      LX = 3;
      break;
    case 41: 
      LX = 2;
      break;
    case 42: 
      LX = 2;
      break;
    case 43: 
      LX = 0;
      break;
    case 44: 
      LX = 0;
      break;
    case 45: 
      LX = 4;
      break;
    case 46: 
      LX = 4;
      break;
    case 47: 
      LX = 1;
      break;
    case 48: 
      LX = 0;
      break;
    case 49: 
      LX = 0;
      break;
    case 50: 
      LX = 1;
      break;
    case 51: 
      LX = 1;
      break;
    case 52: 
      LX = 3;
      break;
    case 53: 
      LX = 1;
      break;
    case 54: 
      LX = 2;
      break;
    case 55: 
      LX = 2;
      break;
    case 56: 
      LX = 1;
      break;
    case 57: 
      LX = 6;
      break;
    case 58: 
      LX = 0;
      break;
    case 59: 
      LX = 0;
      break;
    case 60: 
      LX = 3;
      break;
    case 61: 
      LX = 3;
      break;
    case 62: 
      LX = 0;
      break;
    case 63: 
      LX = 5;
      break;
    case 64: 
      LX = 5;
      break;
    case 65: 
      LX = 1;
      break;
    case 66: 
      LX = 0;
      break;
    case 67: 
      LX = 0;
      break;
    case 68: 
      LX = 1;
      break;
    case 69: 
      LX = 4;
      break;
    case 70: 
      LX = 4;
      break;
    case 71: 
      LX = 2;
      break;
    case 72: 
      LX = 0;
      break;
    case 73: 
      LX = 2;
      break;
    case 74: 
      LX = 0;
      break;
    case 75: 
      LX = 3;
      break;
    case 76: 
      LX = 3;
      break;
    case 77: 
      LX = 1;
      break;
    case 78: 
      LX = 2;
      break;
    case 79: 
      LX = 2;
      break;
    case 80: 
      LX = 1;
      break;
    case 81: 
      LX = 4;
      break;
    case 82: 
      LX = 1;
      break;
    case 83: 
      LX = 1;
      break;
    case 84: 
      LX = 2;
      break;
  }
  return LX;
}
int ly(int AOnumber) {
  int LY = 0;
  switch (AOnumber) {
    case 2:
      LY = 0;
      break;
    case 3:
      LY = 1;
      break;
    case 4:
      LY = 0;
      break;
    case 5:
      LY = 0;
      break;
    case 6:
      LY = 2;
      break;
    case 7:
      LY = 0;
      break;
    case 8:
      LY = 1;
      break;
    case 9:
      LY = 0;
      break;
    case 10: 
      LY = 1;
      break;
    case 11: 
      LY = 0;
      break;
    case 12: 
      LY = 3;
      break;
    case 13: 
      LY = 0;
      break;
    case 14: 
      LY = 1;
      break;
    case 15: 
      LY = 0;
      break;
    case 16: 
      LY = 2;
      break;
    case 17: 
      LY = 2;
      break;
    case 18: 
      LY = 0;
      break;
    case 19: 
      LY = 1;
      break;
    case 20: 
      LY = 1;
      break;
    case 21: 
      LY = 0;
      break;
    case 22: 
      LY = 4;
      break;
    case 23: 
      LY = 0;
      break;
    case 24: 
      LY = 1;
      break;
    case 25: 
      LY = 0;
      break;
    case 26: 
      LY = 3;
      break;
    case 27: 
      LY = 3;
      break;
    case 28: 
      LY = 0;
      break;
    case 29: 
      LY = 1;
      break;
    case 30: 
      LY = 2;
      break;
    case 31: 
      LY = 0;
      break;
    case 32: 
      LY = 2;
      break;
    case 33: 
      LY = 1;
      break;
    case 34: 
      LY = 2;
      break;
    case 35: 
      LY = 1;
      break;
    case 36: 
      LY = 0;
      break;
    case 37: 
      LY = 5;
      break;
    case 38: 
      LY = 0;
      break;
    case 39: 
      LY = 2;
      break;
    case 40: 
      LY = 0;
      break;
    case 41: 
      LY = 3;
      break;
    case 42: 
      LY = 0;
      break;
    case 43: 
      LY = 3;
      break;
    case 44: 
      LY = 2;
      break;
    case 45: 
      LY = 1;
      break;
    case 46: 
      LY = 0;
      break;
    case 47: 
      LY = 4;
      break;
    case 48: 
      LY = 4;
      break;
    case 49: 
      LY = 1;
      break;
    case 50: 
      LY = 0;
      break;
    case 51: 
      LY = 1;
      break;
    case 52: 
      LY = 1;
      break;
    case 53: 
      LY = 3;
      break;
    case 54: 
      LY = 2;
      break;
    case 55: 
      LY = 1;
      break;
    case 56: 
      LY = 2;
      break;
    case 57: 
      LY = 0;
      break;
    case 58: 
      LY = 6;
      break;
    case 59: 
      LY = 0;
      break;
    case 60: 
      LY = 3;
      break;
    case 61: 
      LY = 0;
      break;
    case 62: 
      LY = 3;
      break;
    case 63: 
      LY = 1;
      break;
    case 64: 
      LY = 0;
      break;
    case 65: 
      LY = 0;
      break;
    case 66: 
      LY = 1;
      break;
    case 67: 
      LY = 5;
      break;
    case 68: 
      LY = 5;
      break;
    case 69: 
      LY = 2;
      break;
    case 70: 
      LY = 0;
      break;
    case 71: 
      LY = 0;
      break;
    case 72: 
      LY = 2;
      break;
    case 73: 
      LY = 4;
      break;
    case 74: 
      LY = 4;
      break;
    case 75: 
      LY = 2;
      break;
    case 76: 
      LY = 1;
      break;
    case 77: 
      LY = 3;
      break;
    case 78: 
      LY = 1;
      break;
    case 79: 
      LY = 3;
      break;
    case 80: 
      LY = 2;
      break;
    case 81: 
      LY = 1;
      break;
    case 82: 
      LY = 4;
      break;
    case 83: 
      LY = 1;
      break;
    case 84: 
      LY = 2;
      break;
  }
  return LY;
}
int lz(int AOnumber) {
  int LZ = 0;
  switch (AOnumber) {
    case 2:
      LZ = 0;
      break;
    case 3:
      LZ = 0;
      break;
    case 4:
      LZ = 1;
      break;
    case 5:
      LZ = 0;
      break;
    case 6:
      LZ = 0;
      break;
    case 7:
      LZ = 2;
      break;
    case 8:
      LZ = 0;
      break;
    case 9:
      LZ = 1;
      break;
    case 10: 
      LZ = 1;
      break;
    case 11: 
      LZ = 0;
      break;
    case 12: 
      LZ = 0;
      break;
    case 13: 
      LZ = 3;
      break;
    case 14: 
      LZ = 0;
      break;
    case 15: 
      LZ = 1;
      break;
    case 16: 
      LZ = 0;
      break;
    case 17: 
      LZ = 1;
      break;
    case 18: 
      LZ = 2;
      break;
    case 19: 
      LZ = 2;
      break;
    case 20: 
      LZ = 1;
      break;
    case 21: 
      LZ = 0;
      break;
    case 22: 
      LZ = 0;
      break;
    case 23: 
      LZ = 4;
      break;
    case 24: 
      LZ = 0;
      break;
    case 25: 
      LZ = 1;
      break;
    case 26: 
      LZ = 0;
      break;
    case 27: 
      LZ = 1;
      break;
    case 28: 
      LZ = 3;
      break;
    case 29: 
      LZ = 3;
      break;
    case 30: 
      LZ = 0;
      break;
    case 31: 
      LZ = 2;
      break;
    case 32: 
      LZ = 2;
      break;
    case 33: 
      LZ = 1;
      break;
    case 34: 
      LZ = 1;
      break;
    case 35: 
      LZ = 2;
      break;
    case 36: 
      LZ = 0;
      break;
    case 37: 
      LZ = 0;
      break;
    case 38: 
      LZ = 5;
      break;
    case 39: 
      LZ = 0;
      break;
    case 40: 
      LZ = 2;
      break;
    case 41: 
      LZ = 0;
      break;
    case 42: 
      LZ = 3;
      break;
    case 43: 
      LZ = 2;
      break;
    case 44: 
      LZ = 3;
      break;
    case 45: 
      LZ = 0;
      break;
    case 46: 
      LZ = 1;
      break;
    case 47: 
      LZ = 0;
      break;
    case 48: 
      LZ = 1;
      break;
    case 49: 
      LZ = 4;
      break;
    case 50: 
      LZ = 4;
      break;
    case 51: 
      LZ = 3;
      break;
    case 52: 
      LZ = 1;
      break;
    case 53: 
      LZ = 1;
      break;
    case 54: 
      LZ = 1;
      break;
    case 55: 
      LZ = 2;
      break;
    case 56: 
      LZ = 2;
      break;
    case 57: 
      LZ = 0;
      break;
    case 58: 
      LZ = 0;
      break;
    case 59: 
      LZ = 6;
      break;
    case 60: 
      LZ = 0;
      break;
    case 61: 
      LZ = 3;
      break;
    case 62: 
      LZ = 3;
      break;
    case 63: 
      LZ = 0;
      break;
    case 64: 
      LZ = 1;
      break;
    case 65: 
      LZ = 5;
      break;
    case 66: 
      LZ = 5;
      break;
    case 67: 
      LZ = 1;
      break;
    case 68: 
      LZ = 0;
      break;
    case 69: 
      LZ = 0;
      break;
    case 70: 
      LZ = 2;
      break;
    case 71: 
      LZ = 4;
      break;
    case 72: 
      LZ = 4;
      break;
    case 73: 
      LZ = 0;
      break;
    case 74: 
      LZ = 2;
      break;
    case 75: 
      LZ = 1;
      break;
    case 76: 
      LZ = 2;
      break;
    case 77: 
      LZ = 2;
      break;
    case 78: 
      LZ = 3;
      break;
    case 79: 
      LZ = 1;
      break;
    case 80: 
      LZ = 3;
      break;
    case 81: 
      LZ = 1;
      break;
    case 82: 
      LZ = 1;
      break;
    case 83: 
      LZ = 4;
      break;
    case 84: 
      LZ = 2;
      break;
  }
  return LZ;
}

#endif //_Orbital_Numeration_
