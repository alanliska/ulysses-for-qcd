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

#ifndef _EISOL_AUX_Parameters_
#define _EISOL_AUX_Parameters_

//description:
//parameters for calculation of atomic heats of formation to be used in NDDO methods using the EISOL procedure

int NrS(int atmnr) {
  //function returning the number of s electrons for each atomic number
  int nrs = 0;
  switch (atmnr) {
    case 1: 
      nrs = 1;
      break;
    case 2: 
      nrs = 2;
      break;
    case 3: 
      nrs = 1;
      break;
    case 4: 
      nrs = 2;
      break;
    case 5: 
      nrs = 2;
      break;
    case 6: 
      nrs = 2;
      break;
    case 7: 
      nrs = 2;
      break;
    case 8: 
      nrs = 2;
      break;
    case 9: 
      nrs = 2;
      break;
    case 10: 
      nrs = 0;
      break;
    case 11: 
      nrs = 1;
      break;
    case 12: 
      nrs = 2;
      break;
    case 13: 
      nrs = 2;
      break;
    case 14: 
      nrs = 2;
      break;
    case 15: 
      nrs = 2;
      break;
    case 16: 
      nrs = 2;
      break;
    case 17: 
      nrs = 2;
      break;
    case 18: 
      nrs = 0;
      break;
    case 19: 
      nrs = 1;
      break;
    case 20: 
      nrs = 2;
      break;
    case 21: 
      nrs = 2;
      break;
    case 22: 
      nrs = 2;
      break;
    case 23: 
      nrs = 2;
      break;
    case 24: 
      nrs = 1;
      break;
    case 25: 
      nrs = 2;
      break;
    case 26: 
      nrs = 2;
      break;
    case 27: 
      nrs = 2;
      break;
    case 28: 
      nrs = 2;
      break;
    case 29: 
      nrs = 1;
      break;
    case 30: 
      nrs = 2;
      break;
    case 31: 
      nrs = 2;
      break;
    case 32: 
      nrs = 2;
      break;
    case 33: 
      nrs = 2;
      break;
    case 34: 
      nrs = 2;
      break;
    case 35: 
      nrs = 2;
      break;
    case 36: 
      nrs = 0;
      break;
    case 37: 
      nrs = 1;
      break;
    case 38: 
      nrs = 2;
      break;
    case 39: 
      nrs = 2;
      break;
    case 40: 
      nrs = 2;
      break;
    case 41: 
      nrs = 1;
      break;
    case 42: 
      nrs = 1;
      break;
    case 43: 
      nrs = 2;
      break;
    case 44: 
      nrs = 1;
      break;
    case 45: 
      nrs = 1;
      break;
    case 46: 
      nrs = 0;
      break;
    case 47: 
      nrs = 1;
      break;
    case 48: 
      nrs = 2;
      break;
    case 49: 
      nrs = 2;
      break;
    case 50: 
      nrs = 2;
      break;
    case 51: 
      nrs = 2;
      break;
    case 52: 
      nrs = 2;
      break;
    case 53: 
      nrs = 2;
      break;
    case 54: 
      nrs = 0;
      break;
    case 55: 
      nrs = 1;
      break;
    case 56: 
      nrs = 2;
      break;
    case 57: 
      nrs = 2;
      break;
    case 58: 
      nrs = 0;
      break;
    case 59: 
      nrs = 0;
      break;
    case 60: 
      nrs = 0;
      break;
    case 61: 
      nrs = 0;
      break;
    case 62: 
      nrs = 0;
      break;
    case 63: 
      nrs = 2;
      break;
    case 64: 
      nrs = 2;
      break;
    case 65: 
      nrs = 2;
      break;
    case 66: 
      nrs = 2;
      break;
    case 67: 
      nrs = 2;
      break;
    case 68: 
      nrs = 2;
      break;
    case 69: 
      nrs = 2;
      break;
    case 70: 
      nrs = 2;
      break;
    case 71: 
      nrs = 2;
      break;
    case 72: 
      nrs = 2;
      break;
    case 73: 
      nrs = 2;
      break;
    case 74: 
      nrs = 1;
      break;
    case 75: 
      nrs = 2;
      break;
    case 76: 
      nrs = 2;
      break;
    case 77: 
      nrs = 2;
      break;
    case 78: 
      nrs = 1;
      break;
    case 79: 
      nrs = 1;
      break;
    case 80: 
      nrs = 2;
      break;
    case 81: 
      nrs = 2;
      break;
    case 82: 
      nrs = 2;
      break;
    case 83: 
      nrs = 2;
      break;
    case 84: 
      nrs = 2;
      break;
    case 85: 
      nrs = 2;
      break;
    case 86: 
      nrs = 0;
      break;
  }
  return nrs;
}
int NrP(int atmnr) {
  //function returning the number of p electrons for each atomic number
  int nrp = 0;
  switch (atmnr) {
    case 1:
      nrp = 0;
      break;
    case 2:
      nrp = 0;
      break;
    case 3:
      nrp = 0;
      break;
    case 4:
      nrp = 0;
      break;
    case 5:
      nrp = 1;
      break;
    case 6:
      nrp = 2;
      break;
    case 7:
      nrp = 3;
      break;
    case 8:
      nrp = 4;
      break;
    case 9:
      nrp = 5;
      break;
    case 10:
      nrp = 6;
      break;
    case 11:
      nrp = 0;
      break;
    case 12:
      nrp = 0;
      break;
    case 13:
      nrp = 1;
      break;
    case 14:
      nrp = 2;
      break;
    case 15:
      nrp = 3;
      break;
    case 16:
      nrp = 4;
      break;
    case 17:
      nrp = 5;
      break;
    case 18:
      nrp = 6;
      break;
    case 19:
      nrp = 0;
      break;
    case 20:
      nrp = 0;
      break;
    case 21:
      nrp = 0;
      break;
    case 22:
      nrp = 0;
      break;
    case 23:
      nrp = 0;
      break;
    case 24:
      nrp = 0;
      break;
    case 25:
      nrp = 0;
      break;
    case 26:
      nrp = 0;
      break;
    case 27:
      nrp = 0;
      break;
    case 28:
      nrp = 0;
      break;
    case 29:
      nrp = 0;
      break;
    case 30:
      nrp = 0;
      break;
    case 31:
      nrp = 1;
      break;
    case 32:
      nrp = 2;
      break;
    case 33:
      nrp = 3;
      break;
    case 34:
      nrp = 4;
      break;
    case 35:
      nrp = 5;
      break;
    case 36:
      nrp = 6;
      break;
    case 37:
      nrp = 0;
      break;
    case 38:
      nrp = 0;
      break;
    case 39:
      nrp = 0;
      break;
    case 40:
      nrp = 0;
      break;
    case 41:
      nrp = 0;
      break;
    case 42:
      nrp = 0;
      break;
    case 43:
      nrp = 0;
      break;
    case 44:
      nrp = 0;
      break;
    case 45:
      nrp = 0;
      break;
    case 46:
      nrp = 0;
      break;
    case 47:
      nrp = 0;
      break;
    case 48:
      nrp = 0;
      break;
    case 49:
      nrp = 1;
      break;
    case 50:
      nrp = 2;
      break;
    case 51:
      nrp = 3;
      break;
    case 52:
      nrp = 4;
      break;
    case 53:
      nrp = 5;
      break;
    case 54:
      nrp = 6;
      break;
    case 55:
      nrp = 0;
      break;
    case 56:
      nrp = 0;
      break;
    case 57:
      nrp = 0;
      break;
    case 58:
      nrp = 0;
      break;
    case 59:
      nrp = 0;
      break;
    case 60:
      nrp = 0;
      break;
    case 61:
      nrp = 0;
      break;
    case 62:
      nrp = 0;
      break;
    case 63:
      nrp = 0;
      break;
    case 64:
      nrp = 0;
      break;
    case 65:
      nrp = 0;
      break;
    case 66:
      nrp = 0;
      break;
    case 67:
      nrp = 0;
      break;
    case 68:
      nrp = 0;
      break;
    case 69:
      nrp = 0;
      break;
    case 70:
      nrp = 0;
      break;
    case 71:
      nrp = 0;
      break;
    case 72:
      nrp = 0;
      break;
    case 73:
      nrp = 0;
      break;
    case 74:
      nrp = 0;
      break;
    case 75:
      nrp = 0;
      break;
    case 76:
      nrp = 0;
      break;
    case 77:
      nrp = 0;
      break;
    case 78:
      nrp = 0;
      break;
    case 79:
      nrp = 0;
      break;
    case 80:
      nrp = 0;
      break;
    case 81:
      nrp = 1;
      break;
    case 82:
      nrp = 2;
      break;
    case 83:
      nrp = 3;
      break;
    case 84:
      nrp = 4;
      break;
    case 85:
      nrp = 5;
      break;
    case 86:
      nrp = 6;
      break;
  }
  return nrp;
}
int NrD(int atmnr) {
  //function returning the number of d electrons for each atomic number
  int nrd = 0;
  switch (atmnr) {
    case 1:
      nrd = 0;
      break;
    case 2:
      nrd = 0;
      break;
    case 3:
      nrd = 0;
      break;
    case 4:
      nrd = 0;
      break;
    case 5:
      nrd = 0;
      break;
    case 6:
      nrd = 0;
      break;
    case 7:
      nrd = 0;
      break;
    case 8:
      nrd = 0;
      break;
    case 9:
      nrd = 0;
      break;
    case 10:
      nrd = 0;
      break;
    case 11:
      nrd = 0;
      break;
    case 12:
      nrd = 0;
      break;
    case 13:
      nrd = 0;
      break;
    case 14:
      nrd = 0;
      break;
    case 15:
      nrd = 0;
      break;
    case 16:
      nrd = 0;
      break;
    case 17:
      nrd = 0;
      break;
    case 18:
      nrd = 0;
      break;
    case 19:
      nrd = 0;
      break;
    case 20:
      nrd = 0;
      break;
    case 21:
      nrd = 1;
      break;
    case 22:
      nrd = 2;
      break;
    case 23:
      nrd = 3;
      break;
    case 24:
      nrd = 5;
      break;
    case 25:
      nrd = 5;
      break;
    case 26:
      nrd = 6;
      break;
    case 27:
      nrd = 7;
      break;
    case 28:
      nrd = 8;
      break;
    case 29:
      nrd = 10;
      break;
    case 30:
      nrd = 0;
      break;
    case 31:
      nrd = 0;
      break;
    case 32:
      nrd = 0;
      break;
    case 33:
      nrd = 0;
      break;
    case 34:
      nrd = 0;
      break;
    case 35:
      nrd = 0;
      break;
    case 36:
      nrd = 0;
      break;
    case 37:
      nrd = 0;
      break;
    case 38:
      nrd = 0;
      break;
    case 39:
      nrd = 1;
      break;
    case 40:
      nrd = 2;
      break;
    case 41:
      nrd = 4;
      break;
    case 42:
      nrd = 5;
      break;
    case 43:
      nrd = 5;
      break;
    case 44:
      nrd = 7;
      break;
    case 45:
      nrd = 8;
      break;
    case 46:
      nrd = 10;
      break;
    case 47:
      nrd = 10;
      break;
    case 48:
      nrd = 0;
      break;
    case 49:
      nrd = 0;
      break;
    case 50:
      nrd = 0;
      break;
    case 51:
      nrd = 0;
      break;
    case 52:
      nrd = 0;
      break;
    case 53:
      nrd = 0;
      break;
    case 54:
      nrd = 0;
      break;
    case 55:
      nrd = 0;
      break;
    case 56:
      nrd = 0;
      break;
    case 57:
      nrd = 1;
      break;
    case 58:
      nrd = 0;
      break;
    case 59:
      nrd = 0;
      break;
    case 60:
      nrd = 0;
      break;
    case 61:
      nrd = 0;
      break;
    case 62:
      nrd = 0;
      break;
    case 63:
      nrd = 0;
      break;
    case 64:
      nrd = 0;
      break;
    case 65:
      nrd = 0;
      break;
    case 66:
      nrd = 0;
      break;
    case 67:
      nrd = 0;
      break;
    case 68:
      nrd = 0;
      break;
    case 69:
      nrd = 0;
      break;
    case 70:
      nrd = 0;
      break;
    case 71:
      nrd = 1;
      break;
    case 72:
      nrd = 2;
      break;
    case 73:
      nrd = 3;
      break;
    case 74:
      nrd = 5;
      break;
    case 75:
      nrd = 5;
      break;
    case 76:
      nrd = 6;
      break;
    case 77:
      nrd = 7;
      break;
    case 78:
      nrd = 9;
      break;
    case 79:
      nrd = 10;
      break;
    case 80:
      nrd = 0;
      break;
    case 81:
      nrd = 0;
      break;
    case 82:
      nrd = 0;
      break;
    case 83:
      nrd = 0;
      break;
    case 84:
      nrd = 0;
      break;
    case 85:
      nrd = 0;
      break;
    case 86:
      nrd = 0;
      break;
  }
  return nrd;
}
int Nr016(int atmnr) {
  //function returning the number of first terms for the <ss|dd> integrals
  int nr016 = 0;
  switch (atmnr) {
    case 1:
      nr016 = 0;
      break;
    case 2:
      nr016 = 0;
      break;
    case 3:
      nr016 = 0;
      break;
    case 4:
      nr016 = 0;
      break;
    case 5:
      nr016 = 0;
      break;
    case 6:
      nr016 = 0;
      break;
    case 7:
      nr016 = 0;
      break;
    case 8:
      nr016 = 0;
      break;
    case 9:
      nr016 = 0;
      break;
    case 10:
      nr016 = 0;
      break;
    case 11:
      nr016 = 0;
      break;
    case 12:
      nr016 = 0;
      break;
    case 13:
      nr016 = 0;
      break;
    case 14:
      nr016 = 0;
      break;
    case 15:
      nr016 = 0;
      break;
    case 16:
      nr016 = 0;
      break;
    case 17:
      nr016 = 0;
      break;
    case 18:
      nr016 = 0;
      break;
    case 19:
      nr016 = 0;
      break;
    case 20:
      nr016 = 0;
      break;
    case 21:
      nr016 = 2;
      break;
    case 22:
      nr016 = 4;
      break;
    case 23:
      nr016 = 6;
      break;
    case 24:
      nr016 = 5;
      break;
    case 25:
      nr016 = 10;
      break;
    case 26:
      nr016 = 12;
      break;
    case 27:
      nr016 = 14;
      break;
    case 28:
      nr016 = 16;
      break;
    case 29:
      nr016 = 10;
      break;
    case 30:
      nr016 = 0;
      break;
    case 31:
      nr016 = 0;
      break;
    case 32:
      nr016 = 0;
      break;
    case 33:
      nr016 = 0;
      break;
    case 34:
      nr016 = 0;
      break;
    case 35:
      nr016 = 0;
      break;
    case 36:
      nr016 = 0;
      break;
    case 37:
      nr016 = 0;
      break;
    case 38:
      nr016 = 0;
      break;
    case 39:
      nr016 = 4;
      break;
    case 40:
      nr016 = 4;
      break;
    case 41:
      nr016 = 4;
      break;
    case 42:
      nr016 = 5;
      break;
    case 43:
      nr016 = 10;
      break;
    case 44:
      nr016 = 7;
      break;
    case 45:
      nr016 = 8;
      break;
    case 46:
      nr016 = 0;
      break;
    case 47:
      nr016 = 10;
      break;
    case 48:
      nr016 = 0;
      break;
    case 49:
      nr016 = 0;
      break;
    case 50:
      nr016 = 0;
      break;
    case 51:
      nr016 = 0;
      break;
    case 52:
      nr016 = 0;
      break;
    case 53:
      nr016 = 0;
      break;
    case 54:
      nr016 = 0;
      break;
    case 55:
      nr016 = 0;
      break;
    case 56:
      nr016 = 0;
      break;
    case 57:
      nr016 = 0;
      break;
    case 58:
      nr016 = 0;
      break;
    case 59:
      nr016 = 0;
      break;
    case 60:
      nr016 = 0;
      break;
    case 61:
      nr016 = 0;
      break;
    case 62:
      nr016 = 0;
      break;
    case 63:
      nr016 = 0;
      break;
    case 64:
      nr016 = 0;
      break;
    case 65:
      nr016 = 0;
      break;
    case 66:
      nr016 = 0;
      break;
    case 67:
      nr016 = 0;
      break;
    case 68:
      nr016 = 0;
      break;
    case 69:
      nr016 = 0;
      break;
    case 70:
      nr016 = 0;
      break;
    case 71:
      nr016 = 0;
      break;
    case 72:
      nr016 = 4;
      break;
    case 73:
      nr016 = 6;
      break;
    case 74:
      nr016 = 8;
      break;
    case 75:
      nr016 = 10;
      break;
    case 76:
      nr016 = 12;
      break;
    case 77:
      nr016 = 14;
      break;
    case 78:
      nr016 = 9;
      break;
    case 79:
      nr016 = 10;
      break;
    case 80:
      nr016 = 0;
      break;
    case 81:
      nr016 = 0;
      break;
    case 82:
      nr016 = 0;
      break;
    case 83:
      nr016 = 0;
      break;
    case 84:
      nr016 = 0;
      break;
    case 85:
      nr016 = 0;
      break;
    case 86:
      nr016 = 0;
      break;
  }
  return nr016;
}
int Nr066(int atmnr) {
  //function returning the number of first terms for the <dd|dd> integrals
  int nr066 = 0;
  switch (atmnr) {
    case 1:
      nr066 = 0;
      break;
    case 2:
      nr066 = 0;
      break;
    case 3:
      nr066 = 0;
      break;
    case 4:
      nr066 = 0;
      break;
    case 5:
      nr066 = 0;
      break;
    case 6:
      nr066 = 0;
      break;
    case 7:
      nr066 = 0;
      break;
    case 8:
      nr066 = 0;
      break;
    case 9:
      nr066 = 0;
      break;
    case 10:
      nr066 = 0;
      break;
    case 11:
      nr066 = 0;
      break;
    case 12:
      nr066 = 0;
      break;
    case 13:
      nr066 = 0;
      break;
    case 14:
      nr066 = 0;
      break;
    case 15:
      nr066 = 0;
      break;
    case 16:
      nr066 = 0;
      break;
    case 17:
      nr066 = 0;
      break;
    case 18:
      nr066 = 0;
      break;
    case 19:
      nr066 = 0;
      break;
    case 20:
      nr066 = 0;
      break;
    case 21:
      nr066 = 0;
      break;
    case 22:
      nr066 = 1;
      break;
    case 23:
      nr066 = 3;
      break;
    case 24:
      nr066 = 10;
      break;
    case 25:
      nr066 = 10;
      break;
    case 26:
      nr066 = 15;
      break;
    case 27:
      nr066 = 21;
      break;
    case 28:
      nr066 = 28;
      break;
    case 29:
      nr066 = 45;
      break;
    case 30:
      nr066 = 0;
      break;
    case 31:
      nr066 = 0;
      break;
    case 32:
      nr066 = 0;
      break;
    case 33:
      nr066 = 0;
      break;
    case 34:
      nr066 = 0;
      break;
    case 35:
      nr066 = 0;
      break;
    case 36:
      nr066 = 0;
      break;
    case 37:
      nr066 = 0;
      break;
    case 38:
      nr066 = 0;
      break;
    case 39:
      nr066 = 0;
      break;
    case 40:
      nr066 = 1;
      break;
    case 41:
      nr066 = 6;
      break;
    case 42:
      nr066 = 10;
      break;
    case 43:
      nr066 = 10;
      break;
    case 44:
      nr066 = 21;
      break;
    case 45:
      nr066 = 28;
      break;
    case 46:
      nr066 = 45;
      break;
    case 47:
      nr066 = 45;
      break;
    case 48:
      nr066 = 0;
      break;
    case 49:
      nr066 = 0;
      break;
    case 50:
      nr066 = 0;
      break;
    case 51:
      nr066 = 0;
      break;
    case 52:
      nr066 = 0;
      break;
    case 53:
      nr066 = 0;
      break;
    case 54:
      nr066 = 0;
      break;
    case 55:
      nr066 = 0;
      break;
    case 56:
      nr066 = 0;
      break;
    case 57:
      nr066 = 0;
      break;
    case 58:
      nr066 = 0;
      break;
    case 59:
      nr066 = 0;
      break;
    case 60:
      nr066 = 0;
      break;
    case 61:
      nr066 = 0;
      break;
    case 62:
      nr066 = 0;
      break;
    case 63:
      nr066 = 0;
      break;
    case 64:
      nr066 = 0;
      break;
    case 65:
      nr066 = 0;
      break;
    case 66:
      nr066 = 0;
      break;
    case 67:
      nr066 = 0;
      break;
    case 68:
      nr066 = 0;
      break;
    case 69:
      nr066 = 0;
      break;
    case 70:
      nr066 = 0;
      break;
    case 71:
      nr066 = 0;
      break;
    case 72:
      nr066 = 1;
      break;
    case 73:
      nr066 = 3;
      break;
    case 74:
      nr066 = 6;
      break;
    case 75:
      nr066 = 10;
      break;
    case 76:
      nr066 = 15;
      break;
    case 77:
      nr066 = 21;
      break;
    case 78:
      nr066 = 36;
      break;
    case 79:
      nr066 = 45;
      break;
    case 80:
      nr066 = 0;
      break;
    case 81:
      nr066 = 0;
      break;
    case 82:
      nr066 = 0;
      break;
    case 83:
      nr066 = 0;
      break;
    case 84:
      nr066 = 0;
      break;
    case 85:
      nr066 = 0;
      break;
    case 86:
      nr066 = 0;
      break;
  }
  return nr066;
}
int Nr244(int atmnr) {
  //function returning the number of <sd|sd> integrals
  int nr244 = 0;
  switch (atmnr) {
    case 1:
      nr244 = 0;
      break;
    case 2:
      nr244 = 0;
      break;
    case 3:
      nr244 = 0;
      break;
    case 4:
      nr244 = 0;
      break;
    case 5:
      nr244 = 0;
      break;
    case 6:
      nr244 = 0;
      break;
    case 7:
      nr244 = 0;
      break;
    case 8:
      nr244 = 0;
      break;
    case 9:
      nr244 = 0;
      break;
    case 10:
      nr244 = 0;
      break;
    case 11:
      nr244 = 0;
      break;
    case 12:
      nr244 = 0;
      break;
    case 13:
      nr244 = 0;
      break;
    case 14:
      nr244 = 0;
      break;
    case 15:
      nr244 = 0;
      break;
    case 16:
      nr244 = 0;
      break;
    case 17:
      nr244 = 0;
      break;
    case 18:
      nr244 = 0;
      break;
    case 19:
      nr244 = 0;
      break;
    case 20:
      nr244 = 0;
      break;
    case 21:
      nr244 = 1;
      break;
    case 22:
      nr244 = 2;
      break;
    case 23:
      nr244 = 3;
      break;
    case 24:
      nr244 = 5;
      break;
    case 25:
      nr244 = 5;
      break;
    case 26:
      nr244 = 6;
      break;
    case 27:
      nr244 = 7;
      break;
    case 28:
      nr244 = 8;
      break;
    case 29:
      nr244 = 5;
      break;
    case 30:
      nr244 = 0;
      break;
    case 31:
      nr244 = 0;
      break;
    case 32:
      nr244 = 0;
      break;
    case 33:
      nr244 = 0;
      break;
    case 34:
      nr244 = 0;
      break;
    case 35:
      nr244 = 0;
      break;
    case 36:
      nr244 = 0;
      break;
    case 37:
      nr244 = 0;
      break;
    case 38:
      nr244 = 0;
      break;
    case 39:
      nr244 = 1;
      break;
    case 40:
      nr244 = 2;
      break;
    case 41:
      nr244 = 4;
      break;
    case 42:
      nr244 = 5;
      break;
    case 43:
      nr244 = 5;
      break;
    case 44:
      nr244 = 5;
      break;
    case 45:
      nr244 = 5;
      break;
    case 46:
      nr244 = 0;
      break;
    case 47:
      nr244 = 5;
      break;
    case 48:
      nr244 = 0;
      break;
    case 49:
      nr244 = 0;
      break;
    case 50:
      nr244 = 0;
      break;
    case 51:
      nr244 = 0;
      break;
    case 52:
      nr244 = 0;
      break;
    case 53:
      nr244 = 0;
      break;
    case 54:
      nr244 = 0;
      break;
    case 55:
      nr244 = 0;
      break;
    case 56:
      nr244 = 0;
      break;
    case 57:
      nr244 = 0;
      break;
    case 58:
      nr244 = 0;
      break;
    case 59:
      nr244 = 0;
      break;
    case 60:
      nr244 = 0;
      break;
    case 61:
      nr244 = 0;
      break;
    case 62:
      nr244 = 0;
      break;
    case 63:
      nr244 = 0;
      break;
    case 64:
      nr244 = 0;
      break;
    case 65:
      nr244 = 0;
      break;
    case 66:
      nr244 = 0;
      break;
    case 67:
      nr244 = 0;
      break;
    case 68:
      nr244 = 0;
      break;
    case 69:
      nr244 = 0;
      break;
    case 70:
      nr244 = 0;
      break;
    case 71:
      nr244 = 0;
      break;
    case 72:
      nr244 = 2;
      break;
    case 73:
      nr244 = 3;
      break;
    case 74:
      nr244 = 4;
      break;
    case 75:
      nr244 = 5;
      break;
    case 76:
      nr244 = 6;
      break;
    case 77:
      nr244 = 7;
      break;
    case 78:
      nr244 = 5;
      break;
    case 79:
      nr244 = 5;
      break;
    case 80:
      nr244 = 0;
      break;
    case 81:
      nr244 = 0;
      break;
    case 82:
      nr244 = 0;
      break;
    case 83:
      nr244 = 0;
      break;
    case 84:
      nr244 = 0;
      break;
    case 85:
      nr244 = 0;
      break;
    case 86:
      nr244 = 0;
      break;
  }
  return nr244;
}
int Nr266(int atmnr) {
  //function returning the number of second terms for the <dd|dd> integrals
  int nr266 = 0;
  switch (atmnr) {
    case 1: 
      nr266 = 0;
      break;
    case 2: 
      nr266 = 0;
      break;
    case 3: 
      nr266 = 0;
      break;
    case 4: 
      nr266 = 0;
      break;
    case 5: 
      nr266 = 0;
      break;
    case 6: 
      nr266 = 0;
      break;
    case 7: 
      nr266 = 0;
      break;
    case 8: 
      nr266 = 0;
      break;
    case 9: 
      nr266 = 0;
      break;
    case 10: 
      nr266 = 0;
      break;
    case 11: 
      nr266 = 0;
      break;
    case 12: 
      nr266 = 0;
      break;
    case 13: 
      nr266 = 0;
      break;
    case 14: 
      nr266 = 0;
      break;
    case 15: 
      nr266 = 0;
      break;
    case 16: 
      nr266 = 0;
      break;
    case 17: 
      nr266 = 0;
      break;
    case 18: 
      nr266 = 0;
      break;
    case 19: 
      nr266 = 0;
      break;
    case 20: 
      nr266 = 0;
      break;
    case 21: 
      nr266 = 0;
      break;
    case 22: 
      nr266 = 8;
      break;
    case 23: 
      nr266 = 15;
      break;
    case 24: 
      nr266 = 35;
      break;
    case 25: 
      nr266 = 35;
      break;
    case 26: 
      nr266 = 35;
      break;
    case 27: 
      nr266 = 43;
      break;
    case 28: 
      nr266 = 50;
      break;
    case 29: 
      nr266 = 70;
      break;
    case 30: 
      nr266 = 0;
      break;
    case 31: 
      nr266 = 0;
      break;
    case 32: 
      nr266 = 0;
      break;
    case 33: 
      nr266 = 0;
      break;
    case 34: 
      nr266 = 0;
      break;
    case 35: 
      nr266 = 0;
      break;
    case 36: 
      nr266 = 0;
      break;
    case 37: 
      nr266 = 0;
      break;
    case 38: 
      nr266 = 0;
      break;
    case 39: 
      nr266 = 0;
      break;
    case 40: 
      nr266 = 8;
      break;
    case 41: 
      nr266 = 21;
      break;
    case 42: 
      nr266 = 35;
      break;
    case 43: 
      nr266 = 35;
      break;
    case 44: 
      nr266 = 43;
      break;
    case 45: 
      nr266 = 50;
      break;
    case 46: 
      nr266 = 70;
      break;
    case 47: 
      nr266 = 70;
      break;
    case 48: 
      nr266 = 0;
      break;
    case 49: 
      nr266 = 0;
      break;
    case 50: 
      nr266 = 0;
      break;
    case 51: 
      nr266 = 0;
      break;
    case 52: 
      nr266 = 0;
      break;
    case 53: 
      nr266 = 0;
      break;
    case 54: 
      nr266 = 0;
      break;
    case 55: 
      nr266 = 0;
      break;
    case 56: 
      nr266 = 0;
      break;
    case 57: 
      nr266 = 0;
      break;
    case 58: 
      nr266 = 0;
      break;
    case 59: 
      nr266 = 0;
      break;
    case 60: 
      nr266 = 0;
      break;
    case 61: 
      nr266 = 0;
      break;
    case 62: 
      nr266 = 0;
      break;
    case 63: 
      nr266 = 0;
      break;
    case 64: 
      nr266 = 0;
      break;
    case 65: 
      nr266 = 0;
      break;
    case 66: 
      nr266 = 0;
      break;
    case 67: 
      nr266 = 0;
      break;
    case 68: 
      nr266 = 0;
      break;
    case 69: 
      nr266 = 0;
      break;
    case 70: 
      nr266 = 0;
      break;
    case 71: 
      nr266 = 0;
      break;
    case 72: 
      nr266 = 8;
      break;
    case 73: 
      nr266 = 15;
      break;
    case 74: 
      nr266 = 21;
      break;
    case 75: 
      nr266 = 35;
      break;
    case 76: 
      nr266 = 35;
      break;
    case 77: 
      nr266 = 43;
      break;
    case 78: 
      nr266 = 56;
      break;
    case 79: 
      nr266 = 70;
      break;
    case 80: 
      nr266 = 0;
      break;
    case 81: 
      nr266 = 0;
      break;
    case 82: 
      nr266 = 0;
      break;
    case 83: 
      nr266 = 0;
      break;
    case 84: 
      nr266 = 0;
      break;
    case 85: 
      nr266 = 0;
      break;
    case 86: 
      nr266 = 0;
      break;
  }
  return nr266;
}
int Nr466(int atmnr) {
  //function returning the number of third terms for the <dd|dd> integrals
  int nr466 = 0;
  switch (atmnr) {
    case 1:
      nr466 = 0;
      break;
    case 2:
      nr466 = 0;
      break;
    case 3:
      nr466 = 0;
      break;
    case 4:
      nr466 = 0;
      break;
    case 5:
      nr466 = 0;
      break;
    case 6:
      nr466 = 0;
      break;
    case 7:
      nr466 = 0;
      break;
    case 8:
      nr466 = 0;
      break;
    case 9:
      nr466 = 0;
      break;
    case 10:
      nr466 = 0;
      break;
    case 11:
      nr466 = 0;
      break;
    case 12:
      nr466 = 0;
      break;
    case 13:
      nr466 = 0;
      break;
    case 14:
      nr466 = 0;
      break;
    case 15:
      nr466 = 0;
      break;
    case 16:
      nr466 = 0;
      break;
    case 17:
      nr466 = 0;
      break;
    case 18:
      nr466 = 0;
      break;
    case 19:
      nr466 = 0;
      break;
    case 20:
      nr466 = 0;
      break;
    case 21:
      nr466 = 0;
      break;
    case 22:
      nr466 = 1;
      break;
    case 23:
      nr466 = 8;
      break;
    case 24:
      nr466 = 35;
      break;
    case 25:
      nr466 = 35;
      break;
    case 26:
      nr466 = 35;
      break;
    case 27:
      nr466 = 36;
      break;
    case 28:
      nr466 = 43;
      break;
    case 29:
      nr466 = 70;
      break;
    case 30:
      nr466 = 0;
      break;
    case 31:
      nr466 = 0;
      break;
    case 32:
      nr466 = 0;
      break;
    case 33:
      nr466 = 0;
      break;
    case 34:
      nr466 = 0;
      break;
    case 35:
      nr466 = 0;
      break;
    case 36:
      nr466 = 0;
      break;
    case 37:
      nr466 = 0;
      break;
    case 38:
      nr466 = 0;
      break;
    case 39:
      nr466 = 0;
      break;
    case 40:
      nr466 = 1;
      break;
    case 41:
      nr466 = 21;
      break;
    case 42:
      nr466 = 35;
      break;
    case 43:
      nr466 = 35;
      break;
    case 44:
      nr466 = 36;
      break;
    case 45:
      nr466 = 43;
      break;
    case 46:
      nr466 = 70;
      break;
    case 47:
      nr466 = 70;
      break;
    case 48:
      nr466 = 0;
      break;
    case 49:
      nr466 = 0;
      break;
    case 50:
      nr466 = 0;
      break;
    case 51:
      nr466 = 0;
      break;
    case 52:
      nr466 = 0;
      break;
    case 53:
      nr466 = 0;
      break;
    case 54:
      nr466 = 0;
      break;
    case 55:
      nr466 = 0;
      break;
    case 56:
      nr466 = 0;
      break;
    case 57:
      nr466 = 0;
      break;
    case 58:
      nr466 = 0;
      break;
    case 59:
      nr466 = 0;
      break;
    case 60:
      nr466 = 0;
      break;
    case 61:
      nr466 = 0;
      break;
    case 62:
      nr466 = 0;
      break;
    case 63:
      nr466 = 0;
      break;
    case 64:
      nr466 = 0;
      break;
    case 65:
      nr466 = 0;
      break;
    case 66:
      nr466 = 0;
      break;
    case 67:
      nr466 = 0;
      break;
    case 68:
      nr466 = 0;
      break;
    case 69:
      nr466 = 0;
      break;
    case 70:
      nr466 = 0;
      break;
    case 71:
      nr466 = 0;
      break;
    case 72:
      nr466 = 1;
      break;
    case 73:
      nr466 = 8;
      break;
    case 74:
      nr466 = 21;
      break;
    case 75:
      nr466 = 35;
      break;
    case 76:
      nr466 = 35;
      break;
    case 77:
      nr466 = 36;
      break;
    case 78:
      nr466 = 56;
      break;
    case 79:
      nr466 = 70;
      break;
    case 80:
      nr466 = 0;
      break;
    case 81:
      nr466 = 0;
      break;
    case 82:
      nr466 = 0;
      break;
    case 83:
      nr466 = 0;
      break;
    case 84:
      nr466 = 0;
      break;
    case 85:
      nr466 = 0;
      break;
    case 86:
      nr466 = 0;
      break;
  }
  return nr466;
}

#endif //_EISOL_AUX_Parameters_
