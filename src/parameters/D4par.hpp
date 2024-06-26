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

#ifndef _Dispersion4_Parameters_
#define _Dispersion4_Parameters_
#include "D4C6par.hpp"

//description:
//parameters for D4 correction term

double s6D4(const std::string & method) {
  //function returning the scaling factor s6
  double ssix = 1.0;
  //GFN2-xTB
  if ((method == "GFN2")||(method == "gfn2")) {ssix = 1.0;}
  return ssix;
}
double s8D4(const std::string & method) {
  //function returning the scaling factor s8; no parameter yet optimized
  double seight = 1.0;
  //GFN2-xTB
  if ((method == "GFN2")||(method == "gfn2")) {seight = 2.7;}
  return seight;
}
double s9D4(const std::string & method) {
  //function returning the scaling factor s9; no parameter yet optimized
  double snine = 0.0;
  //GFN2-xTB
  if ((method == "GFN2")||(method == "gfn2")) {snine = 5.0;}
  return snine;
}
double s10D4(const std::string & method) {
  //function returning the scaling factor s9; no parameter yet optimized
  double sten = 0.0;
  if ((method == "GFN2")||(method == "gfn2")) {sten = 0.0;}
  return sten;
}
double A1D4(const std::string & method) {
  //function returning the a1 parameter in Becke-Johnson damping
  double a1 = 0.40;           //this is just a general value obtained by visually averaging what Grimme has in his program
  //a1 = 0.4289;
  //GFN2-xTB
  if ((method == "GFN2")||(method == "gfn2")) {a1 = 0.52;}
  return a1;
}
double A2D4(const std::string & method) {
  //function returning the a2 parameter in Becke-Johnson damping
  double a2 = 4.0;           //this is just a general value obtained by visually averaging what Grimme has in his program
  //a2 = 4.4407;
  //GFN2-xTB
  if ((method == "GFN2")||(method == "gfn2")) {a2 = 5.0;}
  return a2;
}
double AlphaD4(const std::string & method, int order) {
  //function returning the exponent alpha
  double alpha = 0.0;
  if (order == 9) {
    //alpha = 16.0;
    //GFN2
    if ((method == "GFN2")||(method == "gfn2")) {alpha = 16.0;}
  }
  return alpha;
}
double k4D4() {           //k4 parameter
  return 4.10451;
}
double k5D4() {           //k5 parameter
  return 19.08857;
}
double k6D4() {           //k6 parameter
  return 254.55531485519995; //2.0*11.28174*11.28174; note that this is the inverse of Grimme's parameter
}
double knD4() {           //kn parameter
  return 7.50;
}
double Zeff(int atmnr) {
  //effective atomic numbers
  double zeff = 0.0;
  switch (atmnr) {
    case 1:
      zeff = 1.0;
      break;
    case 2:
      zeff = 2.0;
      break;
    case 3:
      zeff = 3.0;
      break;
    case 4:
      zeff = 4.0;
      break;
    case 5:
      zeff = 5.0;
      break;
    case 6:
      zeff = 6.0;
      break;
    case 7:
      zeff = 7.0;
      break;
    case 8:
      zeff = 8.0;
      break;
    case 9:
      zeff = 9.0;
      break;
    case 10:
      zeff = 10.0;
      break;
    case 11:
      zeff = 11.0;
      break;
    case 12:
      zeff = 12.0;
      break;
    case 13:
      zeff = 13.0;
      break;
    case 14:
      zeff = 14.0;
      break;
    case 15:
      zeff = 15.0;
      break;
    case 16:
      zeff = 16.0;
      break;
    case 17:
      zeff = 17.0;
      break;
    case 18:
      zeff = 18.0;
      break;
    case 19:
      zeff = 19.0;
      break;
    case 20:
      zeff = 20.0;
      break;
    case 21:
      zeff = 21.0;
      break;
    case 22:
      zeff = 22.0;
      break;
    case 23:
      zeff = 23.0;
      break;
    case 24:
      zeff = 24.0;
      break;
    case 25:
      zeff = 25.0;
      break;
    case 26:
      zeff = 26.0;
      break;
    case 27:
      zeff = 27.0;
      break;
    case 28:
      zeff = 28.0;
      break;
    case 29:
      zeff = 29.0;
      break;
    case 30:
      zeff = 30.0;
      break;
    case 31:
      zeff = 31.0;
      break;
    case 32:
      zeff = 32.0;
      break;
    case 33:
      zeff = 33.0;
      break;
    case 34:
      zeff = 34.0;
      break;
    case 35:
      zeff = 35.0;
      break;
    case 36:
      zeff = 36.0;
      break;
    case 37:
      zeff = 9.0;
      break;
    case 38:
      zeff = 10.0;
      break;
    case 39:
      zeff = 11.0;
      break;
    case 40:
      zeff = 12.0;
      break;
    case 41:
      zeff = 13.0;
      break;
    case 42:
      zeff = 14.0;
      break;
    case 43:
      zeff = 15.0;
      break;
    case 44:
      zeff = 16.0;
      break;
    case 45:
      zeff = 17.0;
      break;
    case 46:
      zeff = 18.0;
      break;
    case 47:
      zeff = 19.0;
      break;
    case 48:
      zeff = 20.0;
      break;
    case 49:
      zeff = 21.0;
      break;
    case 50:
      zeff = 22.0;
      break;
    case 51:
      zeff = 23.0;
      break;
    case 52:
      zeff = 24.0;
      break;
    case 53:
      zeff = 25.0;
      break;
    case 54:
      zeff = 26.0;
      break;
    case 55:
      zeff = 9.0;
      break;
    case 56:
      zeff = 10.0;
      break;
    case 57:
      zeff = 11.0;
      break;
    case 58:
      zeff = 30.0;
      break;
    case 59:
      zeff = 31.0;
      break;
    case 60:
      zeff = 32.0;
      break;
    case 61:
      zeff = 33.0;
      break;
    case 62:
      zeff = 34.0;
      break;
    case 63:
      zeff = 35.0;
      break;
    case 64:
      zeff = 36.0;
      break;
    case 65:
      zeff = 37.0;
      break;
    case 66:
      zeff = 38.0;
      break;
    case 67:
      zeff = 39.0;
      break;
    case 68:
      zeff = 40.0;
      break;
    case 69:
      zeff = 41.0;
      break;
    case 70:
      zeff = 42.0;
      break;
    case 71:
      zeff = 43.0;
      break;
    case 72:
      zeff = 12.0;
      break;
    case 73:
      zeff = 13.0;
      break;
    case 74:
      zeff = 14.0;
      break;
    case 75:
      zeff = 15.0;
      break;
    case 76:
      zeff = 16.0;
      break;
    case 77:
      zeff = 17.0;
      break;
    case 78:
      zeff = 18.0;
      break;
    case 79:
      zeff = 19.0;
      break;
    case 80:
      zeff = 20.0;
      break;
    case 81:
      zeff = 21.0;
      break;
    case 82:
      zeff = 22.0;
      break;
    case 83:
      zeff = 23.0;
      break;
    case 84:
      zeff = 24.0;
      break;
    case 85:
      zeff = 25.0;
      break;
    case 86:
      zeff = 26.0;
      break;
    case 87:
      zeff = 9.0;
      break;
    case 88:
      zeff = 10.0;
      break;
    case 89:
      zeff = 11.0;
      break;
    case 90:
      zeff = 30.0;
      break;
    case 91:
      zeff = 31.0;
      break;
    case 92:
      zeff = 32.0;
      break;
    case 93:
      zeff = 33.0;
      break;
    case 94:
      zeff = 34.0;
      break;
    case 95:
      zeff = 35.0;
      break;
    case 96:
      zeff = 36.0;
      break;
    case 97:
      zeff = 37.0;
      break;
    case 98:
      zeff = 38.0;
      break;
    case 99:
      zeff = 39.0;
      break;
    case 100:
      zeff = 40.0;
      break;
    case 101:
      zeff = 41.0;
      break;
    case 102:
      zeff = 42.0;
      break;
    case 103:
      zeff = 43.0;
      break;
    case 104:
      zeff = 12.0;
      break;
    case 105:
      zeff = 13.0;
      break;
    case 106:
      zeff = 14.0;
      break;
    case 107:
      zeff = 15.0;
      break;
    case 108:
      zeff = 16.0;
      break;
    case 109:
      zeff = 17.0;
      break;
    case 110:
      zeff = 18.0;
      break;
    case 111:
      zeff = 19.0;
      break;
    case 112:
      zeff = 20.0;
      break;
    case 113:
      zeff = 21.0;
      break;
    case 114:
      zeff = 22.0;
      break;
    case 115:
      zeff = 23.0;
      break;
    case 116:
      zeff = 24.0;
      break;
    case 117:
      zeff = 25.0;
      break;
    case 118:
      zeff = 26.0;
      break;
  }
  return zeff;
}
int RefN(int atmnr) {
  //size of refq and refh vectors used in D4
  int refn = 2;
  if ((atmnr > 86)&&(atmnr < 112)) {refn = 0;}
  else {
    switch (atmnr) {
      case 1:
        break;
      case 2:
        refn = 1;
        break;
      case 3:
        refn = 3;
        break;
      case 4:
        refn = 4;
        break;
      case 5:
        refn = 5;
        break;
      case 6:
        refn = 7;
        break;
      case 7:
        refn = 5;
        break;
      case 8:
        refn = 4;
        break;
      case 9:
        break;
      case 10:
        refn = 1;
        break;
      case 11:
        refn = 3;
        break;
      case 12:
        refn = 4;
        break;
      case 13:
        refn = 4;
        break;
      case 14:
        refn = 5;
        break;
      case 15:
        refn = 4;
        break;
      case 16:
        refn = 3;
        break;
      case 17:
        break;
      case 18:
        refn = 1;
        break;
      case 19:
        refn = 3;
        break;
      case 20:
        refn = 4;
        break;
      case 21:
        refn = 4;
        break;
      case 22:
        refn = 4;
        break;
      case 23:
        refn = 4;
        break;
      case 24:
        refn = 4;
        break;
      case 25:
        refn = 3;
        break;
      case 26:
        refn = 3;
        break;
      case 27:
        refn = 4;
        break;
      case 28:
        refn = 4;
        break;
      case 29:
        break;
      case 30:
        break;
      case 31:
        refn = 3;
        break;
      case 32:
        refn = 5;
        break;
      case 33:
        refn = 4;
        break;
      case 34:
        refn = 3;
        break;
      case 35:
        break;
      case 36:
        refn = 1;
        break;
      case 37:
        refn = 3;
        break;
      case 38:
        refn = 4;
        break;
      case 39:
        refn = 3;
        break;
      case 40:
        refn = 4;
        break;
      case 41:
        refn = 4;
        break;
      case 42:
        refn = 4;
        break;
      case 43:
        refn = 3;
        break;
      case 44:
        refn = 3;
        break;
      case 45:
        refn = 4;
        break;
      case 46:
        refn = 3;
        break;
      case 47:
        break;
      case 48:
        break;
      case 49:
        refn = 4;
        break;
      case 50:
        refn = 5;
        break;
      case 51:
        refn = 4;
        break;
      case 52:
        refn = 3;
        break;
      case 53:
        break;
      case 54:
        refn = 1;
        break;
      case 55:
        refn = 3;
        break;
      case 56:
        refn = 4;
        break;
      case 57:
        refn = 3;
        break;
      case 58:
        refn = 1;
        break;
      case 59:
        break;
      case 60:
        break;
      case 61:
        break;
      case 62:
        break;
      case 63:
        break;
      case 64:
        break;
      case 65:
        break;
      case 66:
        break;
      case 67:
        break;
      case 68:
        break;
      case 69:
        break;
      case 70:
        break;
      case 71:
        break;
      case 72:
        refn = 4;
        break;
      case 73:
        refn = 4;
        break;
      case 74:
        refn = 3;
        break;
      case 75:
        refn = 3;
        break;
      case 76:
        refn = 3;
        break;
      case 77:
        refn = 5;
        break;
      case 78:
        refn = 3;
        break;
      case 79:
        break;
      case 80:
        break;
      case 81:
        refn = 4;
        break;
      case 82:
        refn = 5;
        break;
      case 83:
        refn = 4;
        break;
      case 84:
        refn = 3;
        break;
      case 85:
        break;
      case 86:
        refn = 1;
        break;
      case 112:
        refn = 1;
        break;
      case 113:
        refn = 1;
        break;
      case 114:
        refn = 1;
        break;
      case 115:
        refn = 1;
        break;
      case 116:
        refn = 1;
        break;
      case 117:
        refn = 1;
        break;
      case 118:
        refn = 1;
        break;
    }
  }
  return refn;
}
void RefQ(int atmnr, std::vector<double> & refq, size_t reference = 0) {
  //refq values for D4
  //two references available
  // -> D4, reference = 0
  // -> GFN2, reference = 1
  //vector not zeroed because we use RefN to access the vector's formal size
  if (reference == 0) {
    switch (atmnr) {
      case 1:
        refq[0] = 0.0;
        refq[1] = 0.0;
        break;
      case 2:
        refq[0] = 0.0;
        break;
      case 3:
        refq[0] = 0.0;
        refq[1] = 0.37061732228445;
        refq[2] = 0.444;
        break;
      case 4:
        refq[0] = 0.0;
        refq[1] = 0.17355187691088;
        refq[2] = 0.31121011558254;
        refq[3] = 0.685;
        break;
      case 5:
        refq[0] = 0.0;
        refq[1] = -0.10464538187607;
        refq[2] = -0.20153708141;
        refq[3] = -0.28222637549158;
        refq[4] = -0.34038820700608;
        break;
      case 6:
        refq[0] = 0.0;
        refq[1] = -0.17123705791311;
        refq[2] = -0.14210003554128;
        refq[3] = -0.22250854670718;
        refq[4] = -0.25653918216622;
        refq[5] = -0.10441353918078;
        refq[6] = 0.21030312699251;
        break;
      case 7:
        refq[0] = 0.0;
        refq[1] = -0.35298592555229;
        refq[2] = -0.30125788296606;
        refq[3] = -0.82429723124422;
        refq[4] = 0.0;
        break;
      case 8:
        refq[0] = 0.0;
        refq[1] = -0.35471470167212;
        refq[2] = -0.5988152606046;
        refq[3] = 0.0;
        break;
      case 9:
        refq[0] = 0.0;
        refq[1] = -0.24773915491465;
        break;
      case 10:
        refq[0] = 0.0;
        break;
      case 11:
        refq[0] = 0.0;
        refq[1] = 0.37895090683878;
        refq[2] = 0.482;
        break;
      case 12:
        refq[0] = 0.0;
        refq[1] = 0.32933884303798;
        refq[2] = 0.51074920852765;
        refq[3] = 0.721;
        break;
      case 13:
        refq[0] = 0.0;
        refq[1] = 0.03709152257935;
        refq[2] = 0.05079767549293;
        refq[3] = 0.05291396113772;
        break;
      case 14:
        refq[0] = 0.0;
        refq[1] = 0.01649127941162;
        refq[2] = -0.00640083798505;
        refq[3] = -0.03084420203156;
        refq[4] = -0.06585975737745;
        break;
      case 15:
        refq[0] = 0.0;
        refq[1] = -0.17987135979307;
        refq[2] = -0.1267112413984;
        refq[3] = -0.17467927422329;
        break;
      case 16:
        refq[0] = 0.0;
        refq[1] = -0.16267234588041;
        refq[2] = -0.23775171448878;
        break;
      case 17:
        refq[0] = 0.0;
        refq[1] = -0.15045057451823;
        break;
      case 18:
        refq[0] = 0.0;
        break;
      case 19:
        refq[0] = 0.0;
        refq[1] = 0.4437624278749;
        refq[2] = 0.62148041473613;
        break;
      case 20:
        refq[0] = 0.0;
        refq[1] = 0.31814452979636;
        refq[2] = 0.50078886624431;
        refq[3] = 0.708;
        break;
      case 21:
        refq[0] = 0.0;
        refq[1] = 0.50144727255861;
        refq[2] = 0.6248537329654;
        refq[3] = 1.198;
        break;
      case 22:
        refq[0] = 0.0;
        refq[1] = 0.32259167829296;
        refq[2] = 0.53505953733387;
        refq[3] = 1.4528449023815;
        break;
      case 23:
        refq[0] = 0.0;
        refq[1] = 0.30538713298844;
        refq[2] = 0.3974593135865;
        refq[3] = 1.29423363276103;
        break;
      case 24:
        refq[0] = 0.0;
        refq[1] = 0.17721817570715;
        refq[2] = 0.16754955881404;
        refq[3] = 0.61467776168125;
        break;
      case 25:
        refq[0] = 0.0;
        refq[1] = 0.26888428905691;
        refq[2] = -0.01438794807264;
        break;
      case 26:
        refq[0] = 0.0;
        refq[1] = 0.1610833014447;
        refq[2] = 0.05493242513468;
        break;
      case 27:
        refq[0] = 0.0;
        refq[1] = 0.2156815158038;
        refq[2] = 0.2059004977811;
        refq[3] = -0.0157789481861;
        break;
      case 28:
        refq[0] = 0.0;
        refq[1] = 0.09850839901243;
        refq[2] = 0.02948298358915;
        refq[3] = 0.02439689318964;
        break;
      case 29:
        refq[0] = 0.0;
        refq[1] = 0.05450628884215;
        break;
      case 30:
        refq[0] = 0.0;
        refq[1] = 0.25796918840353;
        break;
      case 31:
        refq[0] = 0.0;
        refq[1] = 0.02870304944818;
        refq[2] = 0.04511416923234;
        break;
      case 32:
        refq[0] = 0.0;
        refq[1] = 0.00938786608492;
        refq[2] = -0.00592798448224;
        refq[3] = -0.0316158546311;
        refq[4] = -0.06445558548675;
        break;
      case 33:
        refq[0] = 0.0;
        refq[1] = -0.09954248086481;
        refq[2] = -0.11896736577051;
        refq[3] = -0.10225384717106;
        break;
      case 34:
        refq[0] = 0.0;
        refq[1] = -0.12327157100582;
        refq[2] = -0.17535803722366;
        break;
      case 35:
        refq[0] = 0.0;
        refq[1] = -0.14490782122377;
        break;
      case 36:
        refq[0] = 0.0;
        break;
      case 37:
        refq[0] = 0.0;
        refq[1] = 0.44852683448077;
        refq[2] = 0.566;
        break;
      case 38:
        refq[0] = 0.0;
        refq[1] = 0.42643340166796;
        refq[2] = 0.6490281427178;
        refq[3] = 0.843;
        break;
      case 39:
        refq[0] = 0.0;
        refq[1] = 0.5761854611428;
        refq[2] = 0.72312847632101;
        break;
      case 40:
        refq[0] = 0.0;
        refq[1] = 0.34661262612707;
        refq[2] = 0.56547785565182;
        refq[3] = 1.59446204953052;
        break;
      case 41:
        refq[0] = 0.0;
        refq[1] = 0.42830519499403;
        refq[2] = 0.5605094860094;
        refq[3] = 1.638;
        break;
      case 42:
        refq[0] = 0.0;
        refq[1] = 0.46738200215089;
        refq[2] = 0.473386703959;
        refq[3] = 0.59427677788562;
        break;
      case 43:
        refq[0] = 0.0;
        refq[1] = 0.24709379453256;
        refq[2] = 0.18198173917164;
        break;
      case 44:
        refq[0] = 0.0;
        refq[1] = 0.32415456703039;
        refq[2] = 0.36375674498582;
        break;
      case 45:
        refq[0] = 0.0;
        refq[1] = 0.14509956942986;
        refq[2] = 0.18388211183689;
        refq[3] = 0.24549798804335;
        break;
      case 46:
        refq[0] = 0.0;
        refq[1] = 0.22956514930613;
        refq[2] = 0.55444971429388;
        break;
      case 47:
        refq[0] = 0.0;
        refq[1] = 0.10670361771604;
        break;
      case 48:
        refq[0] = 0.0;
        refq[1] = 0.2189854054517;
        break;
      case 49:
        refq[0] = 0.0;
        refq[1] = 0.06391398485994;
        refq[2] = 0.10806425732282;
        refq[3] = 0.14171698994803;
        break;
      case 50:
        refq[0] = 0.0;
        refq[1] = 0.03599066704066;
        refq[2] = 0.04368851386536;
        refq[3] = 0.03463534758499;
        refq[4] = 0.01418221277464;
        break;
      case 51:
        refq[0] = 0.0;
        refq[1] = -0.01439446368141;
        refq[2] = -0.015432015813;
        refq[3] = -0.00985357197761;
        break;
      case 52:
        refq[0] = 0.0;
        refq[1] = -0.09115200126999;
        refq[2] = -0.09678003345023;
        break;
      case 53:
        refq[0] = 0.0;
        refq[1] = -0.13740767685282;
        break;
      case 54:
        refq[0] = 0.0;
        break;
      case 55:
        refq[0] = 0.0;
        refq[1] = 0.46524059775602;
        refq[2] = 0.605;
        break;
      case 56:
        refq[0] = 0.0;
        refq[1] = 0.42198236936188;
        refq[2] = 0.63397117832958;
        refq[3] = 0.82;
        break;
      case 57:
        refq[0] = 0.0;
        refq[1] = 0.65599616057853;
        refq[2] = 0.8308486053033;
        break;
      case 58:
        refq[0] = 0.78733393312177;
        break;
      case 59:
        refq[0] = 0.0;
        refq[1] = 0.78255214121928;
        break;
      case 60:
        refq[0] = 0.0;
        refq[1] = 0.78608430866652;
        break;
      case 61:
        refq[0] = 0.0;
        refq[1] = 0.79124653412356;
        break;
      case 62:
        refq[0] = 0.0;
        refq[1] = 0.78972506109909;
        break;
      case 63:
        refq[0] = 0.0;
        refq[1] = 0.65326281414848;
        break;
      case 64:
        refq[0] = 0.0;
        refq[1] = 0.78994446367069;
        break;
      case 65:
        refq[0] = 0.0;
        refq[1] = 0.7946243928496;
        break;
      case 66:
        refq[0] = 0.0;
        refq[1] = 0.79518107342432;
        break;
      case 67:
        refq[0] = 0.0;
        refq[1] = 0.80930480242405;
        break;
      case 68:
        refq[0] = 0.0;
        refq[1] = 0.81306730354597;
        break;
      case 69:
        refq[0] = 0.0;
        refq[1] = 0.81408370489997;
        break;
      case 70:
        refq[0] = 0.0;
        refq[1] = 0.7841584373539;
        break;
      case 71:
        refq[0] = 0.0;
        refq[1] = 0.8226318727305;
        break;
      case 72:
        refq[0] = 0.0;
        refq[1] = 0.44266476936835;
        refq[2] = 0.68407375970491;
        refq[3] = 1.59380889275784;
        break;
      case 73:
        refq[0] = 0.0;
        refq[1] = 0.23523422154086;
        refq[2] = 0.30980440015125;
        refq[3] = 1.51762325430455;
        break;
      case 74:
        refq[0] = 0.0;
        refq[1] = 0.31319918056428;
        refq[2] = 0.23020705290823;
        break;
      case 75:
        refq[0] = 0.0;
        refq[1] = 0.10191716059118;
        refq[2] = -0.03451303623343;
        break;
      case 76:
        refq[0] = 0.0;
        refq[1] = 0.1986439176555;
        refq[2] = -0.06474020857185;
        break;
      case 77:
        refq[0] = 0.0;
        refq[1] = -0.01399684403519;
        refq[2] = -0.04786532465344;
        refq[3] = -0.07756443085155;
        refq[4] = -0.12042733196092;
        break;
      case 78:
        refq[0] = 0.0;
        refq[1] = 0.00248540621951;
        refq[2] = 0.10111888298865;
        break;
      case 79:
        refq[0] = 0.0;
        refq[1] = 0.02560952126341;
        break;
      case 80:
        refq[0] = 0.0;
        refq[1] = 0.1199114691992;
        break;
      case 81:
        refq[0] = 0.0;
        refq[1] = 0.05732791283824;
        refq[2] = 0.10187349213041;
        refq[3] = 0.14012325640618;
        break;
      case 82:
        refq[0] = 0.0;
        refq[1] = 0.11449566751294;
        refq[2] = 0.1675845932319;
        refq[3] = 0.18265663815354;
        refq[4] = 0.17825662719618;
        break;
      case 83:
        refq[0] = 0.0;
        refq[1] = 0.04276095694715;
        refq[2] = 0.0654518397337;
        refq[3] = 0.07662886467999;
        break;
      case 84:
        refq[0] = 0.0;
        refq[1] = -0.04479748091835;
        refq[2] = -0.02201995669726;
        break;
      case 85:
        refq[0] = 0.0;
        refq[1] = -0.07596582020292;
        break;
      case 86:
        refq[0] = 0.0;
        break;
      case 112:
        refq[0] = 0.0;
        break;
      case 113:
        refq[0] = 0.10288904092396;
        break;
      case 114:
        refq[0] = 0.0;
        break;
      case 116:
        refq[0] = 0.04481663412472;
        break;
      case 115:
        refq[0] = 0.10991065534975;
        break;
      case 117:
        refq[0] = 0.07044735360796;
        break;
      case 118:
        refq[0] = 0.0;
        break;
    }
  }
  else if (reference == 1) {
    switch (atmnr) {
      case 1:
        refq[0] = 0.0;
        refq[1] = 0.0;
        break;
      case 2:
        refq[0] = 0.0;
        break;
      case 3:
        refq[0] = 0.0;
        refq[1] = 0.35924680736483;
        refq[2] = 0.60553105026862;
        break;
      case 4:
        refq[0] = 0.0;
        refq[1] = 0.35318171705506;
        refq[2] = 0.47868871444650;
        refq[3] = -0.65605016262773;
        break;
      case 5:
        refq[0] = -0.00000000000001;
        refq[1] = 0.15557155274674;
        refq[2] = 0.19433603553394;
        refq[3] = 0.24475443834742;
        refq[4] = 0.14920702115395;
        break;
      case 6:
        refq[0] = -0.00000000000001;
        refq[1] = -0.12994792301999;
        refq[2] = -0.07885200894330;
        refq[3] = -0.07676292658733;
        refq[4] = -0.09463126916388;
        refq[5] = -0.02858894127595;
        refq[6] = 0.05061953301341;
        break;
      case 7:
        refq[0] = 0.0;
        refq[1] = -0.17975369721743;
        refq[2] = -0.15258555866677;
        refq[3] = -0.42898107659879;
        refq[4] = 0.0;
        break;
      case 8:
        refq[0] = 0.00000000000001;
        refq[1] = -0.33069420009468;
        refq[2] = -0.56376716743110;
        refq[3] = -0.00000000000002;
        break;
      case 9:
        refq[0] = 0.00000000000001;
        refq[1] = -0.39039862127573;
        break;
      case 10:
        refq[0] = 0.0;
        break;
      case 11:
        refq[0] = 0.0;
        refq[1] = 0.44612149004129;
        refq[2] = 0.69570459023927;
        break;
      case 12:
        refq[0] = 0.0;
        refq[1] = 0.31858696714581;
        refq[2] = 0.41366273285585;
        refq[3] = 0.50571212038828;
        break;
      case 13:
        refq[0] = -0.00000000000001;
        refq[1] = 0.24477359917277;
        refq[2] = 0.36087089581165;
        refq[3] = 0.43660315534647;
        break;
      case 14:
        refq[0] = 0.0;
        refq[1] = 0.16996090538399;
        refq[2] = 0.24464651050076;
        refq[3] = 0.20893161783021;
        refq[4] = 0.36245469258615;
        break;
      case 15:
        refq[0] = -0.00000000000002;
        refq[1] = -0.08781933942405;
        refq[2] = -0.03728416546523;
        refq[3] = -0.11027568692625;
        break;
      case 16:
        refq[0] = 0.0;
        refq[1] = -0.17222533230765;
        refq[2] = -0.34208440200689;
        break;
      case 17:
        refq[0] = 0.0;
        refq[1] = -0.33096614692963;
        break;
      case 18:
        refq[0] = 0.0;
        break;
      case 19:
        refq[0] = 0.0;
        refq[1] = 0.44919564669363;
        refq[2] = 0.70700012466991;
        break;
      case 20:
        refq[0] = 0.0;
        refq[1] = 0.42382301144653;
        refq[2] = 0.53704626474621;
        refq[3] = 0.63912284405165;
        break;
      case 21:
        refq[0] = -0.00000000016423;
        refq[1] = 0.19081329495662;
        refq[2] = 0.57390951001388;
        refq[3] = 1.02910987985129;
        break;
      case 22:
        refq[0] = -0.00000000032847;
        refq[1] = 0.38774830140171;
        refq[2] = 0.97207545160459;
        refq[3] = 1.45284490238150;
        break;
      case 23:
        refq[0] = -0.00000000016528;
        refq[1] = 0.30616936989879;
        refq[2] = 0.39411696843052;
        refq[3] = 0.69749634235669;
        break;
      case 24:
        refq[0] = -0.00000000000208;
        refq[1] = 0.53605366789797;
        refq[2] = 0.32326600457791;
        refq[3] = 0.17434704715481;
        break;
      case 25:
        refq[0] = 0.00000000016527;
        refq[1] = 0.72906362172268;
        refq[2] = -1.38551399356320;
        break;
      case 26:
        refq[0] = 0.00000000032847;
        refq[1] = 0.89971355662264;
        refq[2] = -0.19729173677108;
        break;
      case 27:
        refq[0] = 0.00000000016425;
        refq[1] = 0.33885434982022;
        refq[2] = 0.06673726947568;
        refq[3] = -0.99922236530237;
        break;
      case 28:
        refq[0] = 0.0;
        refq[1] = -0.02893652143275;
        refq[2] = -1.41046704877362;
        refq[3] = -1.49873637353726;
        break;
      case 29:
        refq[0] = 0.0;
        refq[1] = 0.17543548739457;
        break;
      case 30:
        refq[0] = 0.0;
        refq[1] = 0.22218813020625;
        break;
      case 31:
        refq[0] = -0.00000000000001;
        refq[1] = 0.30692560016815;
        refq[2] = 0.25583231781832;
        break;
      case 32:
        refq[0] = -0.00000000000001;
        refq[1] = 0.15565415588149;
        refq[2] = 0.19489066026694;
        refq[3] = 0.07796365954686;
        refq[4] = 0.08330185602026;
        break;
      case 33:
        refq[0] = 0.0;
        refq[1] = 0.05341421485122;
        refq[2] = 0.08123466402805;
        refq[3] = 0.09121599779474;
        break;
      case 34:
        refq[0] = 0.00000000000001;
        refq[1] = -0.14680263633976;
        refq[2] = -0.30550786355706;
        break;
      case 35:
        refq[0] = 0.0;
        refq[1] = -0.26466311607527;
        break;
      case 36:
        refq[0] = 0.0;
        break;
      case 37:
        refq[0] = 0.0;
        refq[1] = 0.38168923853998;
        refq[2] = 0.61945426066521;
        break;
      case 38:
        refq[0] = 0.0;
        refq[1] = 0.35645143561998;
        refq[2] = 0.49180081018397;
        refq[3] = 0.58944749671323;
        break;
      case 39:
        refq[0] = -0.00000000016424;
        refq[1] = 0.24428414680300;
        refq[2] = 0.57366289721017;
        break;
      case 40:
        refq[0] = -0.00000000032846;
        refq[1] = 0.23434724050528;
        refq[2] = 1.20638669472971;
        refq[3] = 1.59446204953052;
        break;
      case 41:
        refq[0] = -0.00000000016530;
        refq[1] = 0.07227863443605;
        refq[2] = 0.29185185542919;
        refq[3] = 1.63800000000000;
        break;
      case 42:
        refq[0] = 0.00000000000215;
        refq[1] = 0.41711608738091;
        refq[2] = -4.49557215843238;
        refq[3] = -0.06533097052885;
        break;
      case 43:
        refq[0] = 0.00000000016529;
        refq[1] = 0.64591594055856;
        refq[2] = -0.23793642312608;
        break;
      case 44:
        refq[0] = 0.00000000032850;
        refq[1] = 0.33761249164822;
        refq[2] = -0.70670795998300;
        break;
      case 45:
        refq[0] = 0.00000000016424;
        refq[1] = 0.13533447502854;
        refq[2] = -0.06808099954071;
        refq[3] = -0.00946705792926;
        break;
      case 46:
        refq[0] = 0.0;
        refq[1] = 0.13729986752066;
        refq[2] = -1.40283988548386;
        break;
      case 47:
        refq[0] = 0.00000000000088;
        refq[1] = 0.07743633228111;
        break;
      case 48:
        refq[0] = 0.0;
        refq[1] = 0.33653584549362;
        break;
      case 49:
        refq[0] = 0.0;
        refq[1] = 0.34258129843160;
        refq[2] = 0.26652372419494;
        refq[3] = 0.20507380029383;
        break;
      case 50:
        refq[0] = -0.00000000000001;
        refq[1] = 0.11709105905770;
        refq[2] = 0.20024277847179;
        refq[3] = 0.07115823557122;
        refq[4] = 0.02040246592599;
        break;
      case 51:
        refq[0] = 0.0;
        refq[1] = 0.11698834506861;
        refq[2] = 0.18480068371662;
        refq[3] = 0.23470973824995;
        break;
      case 52:
        refq[0] = 0.00000000000001;
        refq[1] = -0.00520338227711;
        refq[2] = -0.00151877920060;
        break;
      case 53:
        refq[0] = 0.00000000000001;
        refq[1] = -0.18609272966361;
        break;
      case 54:
        refq[0] = 0.0;
        break;
      case 55:
        refq[0] = 0.0;
        refq[1] = 0.44073166763575;
        refq[2] = 0.68068249931321;
        break;
      case 56:
        refq[0] = 0.0;
        refq[1] = 0.39462085293344;
        refq[2] = 0.56717887047949;
        refq[3] = 0.62956728884112;
        break;
      case 57:
        refq[0] = -0.00000000016424;
        refq[1] = 0.14182729573672;
        refq[2] = 0.61948745805358;
        break;
      case 58:
        refq[0] = 0.54768224095698;
        break;
      case 59:
        refq[0] = -0.00000000021208;
        refq[1] = 0.56091390202022;
        break;
      case 60:
        refq[0] = -0.00000000019384;
        refq[1] = 0.57175663319949;
        break;
      case 61:
        refq[0] = -0.00000000047809;
        refq[1] = 0.57848494573230;
        break;
      case 62:
        refq[0] = -0.00000000067069;
        refq[1] = 0.58363743729144;
        break;
      case 63:
        refq[0] = -0.00000000029832;
        refq[1] = 0.14255010862838;
        break;
      case 64:
        refq[0] = -0.00000000020011;
        refq[1] = 0.58503142163331;
        break;
      case 65:
        refq[0] = -0.00000000017383;
        refq[1] = 0.60523146786190;
        break;
      case 66:
        refq[0] = -0.00000000016680;
        refq[1] = 0.61019655094139;
        break;
      case 67:
        refq[0] = -0.00000000016492;
        refq[1] = 0.62162856785073;
        break;
      case 68:
        refq[0] = -0.00000000016442;
        refq[1] = 0.62470266440724;
        break;
      case 69:
        refq[0] = -0.00000000016430;
        refq[1] = 0.63183032821235;
        break;
      case 70:
        refq[0] = -0.00000000016425;
        refq[1] = 0.60550144239070;
        break;
      case 71:
        refq[0] = -0.00000000016423;
        refq[1] = 0.64112998199707;
        break;
      case 72:
        refq[0] = -0.00000000032846;
        refq[1] = 0.11182987101245;
        refq[2] = 1.03176947330435;
        refq[3] = 1.35275862077626;
        break;
      case 73:
        refq[0] = -0.00000000016531;
        refq[1] = 0.16442900959316;
        refq[2] = 0.50275891608761;
        refq[3] = 1.84791460855243;
        break;
      case 74:
        refq[0] = -0.00000000000215;
        refq[1] = -0.33545451190268;
        refq[2] = -0.88683226154612;
        break;
      case 75:
        refq[0] = 0.00000000000001;
        refq[1] = 0.24763704610429;
        refq[2] = -1.18136092555004;
        break;
      case 76:
        refq[0] = 0.00000000032846;
        refq[1] = 0.39810831827982;
        refq[2] = -0.74972152576721;
        break;
      case 77:
        refq[0] = 0.00000000016428;
        refq[1] = 0.25293750608113;
        refq[2] = 0.07099664544311;
        refq[3] = 0.08089441537446;
        refq[4] = 0.07168988850722;
        break;
      case 78:
        refq[0] = 0.00000000000057;
        refq[1] = 0.06077075971517;
        refq[2] = -1.31315249908318;
        break;
      case 79:
        refq[0] = 0.0;
        refq[1] = 0.16039939518296;
        break;
      case 80:
        refq[0] = 0.0;
        refq[1] = -0.04321502684659;
        break;
      case 81:
        refq[0] = -0.00000000000001;
        refq[1] = 0.30915981534272;
        refq[2] = 0.08292819294371;
        refq[3] = -0.09717467138010;
        break;
      case 82:
        refq[0] = -0.00000000000001;
        refq[1] = 0.25127061741389;
        refq[2] = 0.36893958472478;
        refq[3] = 0.20035516415189;
        refq[4] = 0.01348187443840;
        break;
      case 83:
        refq[0] = 0.0;
        refq[1] = 0.16648528849973;
        refq[2] = 0.27582968178310;
        refq[3] = 0.35078718856604;
        break;
      case 84:
        refq[0] = 0.00000000000001;
        refq[1] = 0.13640050127609;
        refq[2] = 0.28152411931918;
        break;
      case 85:
        refq[0] = 0.00000000000001;
        refq[1] = 0.03800713219800;
        break;
      case 86:
        refq[0] = 0.0;
        break;
      case 112:
        refq[0] = 0.0;
        break;
      case 113:
        refq[0] = 0.10288904092396;
        break;
      case 114:
        refq[0] = 0.0;
        break;
      case 116:
        refq[0] = 0.04481663412472;
        break;
      case 115:
        refq[0] = 0.10991065534975;
        break;
      case 117:
        refq[0] = 0.07044735360796;
        break;
      case 118:
        refq[0] = 0.0;
        break;
    }
  }
}
void RefH(int atmnr, std::vector<double> & refh) {
  //refh values for D4
  //vector not zeroed because we use RefN to access the vector's formal size
  switch (atmnr) {
    case 1:
      refh[0] = 0.0;
      refh[1] = 0.0;
      break;
    case 2:
      refh[0] = 0.0;
      break;
    case 3:
      refh[0] = 0.0;
      refh[1] = -0.35924680736483;
      refh[2] = -0.60553105026862;
      break;
    case 4:
      refh[0] = 0.0;
      refh[1] = -0.35318171705506;
      refh[2] = -0.23934435722325;
      refh[3] = 0.0;
      break;
    case 5:
      refh[0] = 0.0;
      refh[1] = -0.15557155274761;
      refh[2] = -0.09716801776697;
      refh[3] = -0.08158481278247;
      refh[4] = -0.04973567371798;
      break;
    case 6:
      refh[0] = 0.0;
      refh[1] = 0.12994792301999;
      refh[2] = 0.07885200894330;
      refh[3] = 0.03838146329367;
      refh[4] = 0.03154375638796;
      refh[5] = 0.02858894127595;
      refh[6] = -0.05061953301341;
      break;
    case 7:
      refh[0] = 0.0;
      refh[1] = 0.17975369722031;
      refh[2] = 0.15258555866677;
      refh[3] = 0.14299369219960;
      refh[4] = 0.0;
      break;
    case 8:
      refh[0] = 0.0;
      refh[1] = 0.33069420009468;
      refh[2] = 0.28188358371555;
      refh[3] = 0.0;
      break;
    case 9:
      refh[0] = 0.0;
      refh[1] = 0.39039862127574;
      break;
    case 10:
      refh[0] = 0.0;
      break;
    case 11:
      refh[0] = 0.0;
      refh[1] = -0.44612149004129;
      refh[2] = -0.69570459023927;
      break;
    case 12:
      refh[0] = 0.0;
      refh[1] = -0.31858696714581;
      refh[2] = -0.20683136642793;
      refh[3] = 0.0;
      break;
    case 13:
      refh[0] = 0.0;
      refh[1] = -0.24477359917277;
      refh[2] = -0.18043544790582;
      refh[3] = -0.14553438511549;
      break;
    case 14:
      refh[0] = 0.0;
      refh[1] = -0.16996090538399;
      refh[2] = -0.12232325525038;
      refh[3] = -0.10446580891510;
      refh[4] = -0.12081823086205;
      break;
    case 15:
      refh[0] = 0.0;
      refh[1] = 0.08781933942403;
      refh[2] = 0.03728416546523;
      refh[3] = 0.03675856230875;
      break;
    case 16:
      refh[0] = 0.0;
      refh[1] = 0.17222533230765;
      refh[2] = 0.17104220100345;
      break;
    case 17:
      refh[0] = 0.0;
      refh[1] = 0.33096614692963;
      break;
    case 18:
      refh[0] = 0.0;
      break;
    case 19:
      refh[0] = 0.0;
      refh[1] = -0.44919564669363;
      refh[2] = -0.70700012466991;
      break;
    case 20:
      refh[0] = 0.0;
      refh[1] = -0.42382301144653;
      refh[2] = -0.26852313237311;
      refh[3] = -0.63912284405165;
      break;
    case 21:
      refh[0] = 0.0;
      refh[1] = -0.09540664747831;
      refh[2] = -0.19130317000463;
      refh[3] = -0.34303662000000;
      break;
    case 22:
      refh[0] = 0.0;
      refh[1] = -0.19387415070086;
      refh[2] = -0.24301886290115;
      refh[3] = -0.36600000000000;
      break;
    case 23:
      refh[0] = 0.0;
      refh[1] = -0.15308468494977;
      refh[2] = -0.13137232283301;
      refh[3] = -0.13960000000000;
      break;
    case 24:
      refh[0] = 0.0;
      refh[1] = -0.26802683372543;
      refh[2] = -0.02693883371483;
      refh[3] = -0.12895656370357;
      break;
    case 25:
      refh[0] = 0.0;
      refh[1] = -0.36453181087692;
      refh[2] = 0.12595581759665;
      break;
    case 26:
      refh[0] = 0.0;
      refh[1] = -0.44985677831131;
      refh[2] = 0.01972917367711;
      break;
    case 27:
      refh[0] = 0.0;
      refh[1] = -0.16942717490891;
      refh[2] = -0.02224575649166;
      refh[3] = 0.11102470725582;
      break;
    case 28:
      refh[0] = 0.0;
      refh[1] = 0.01446826083399;
      refh[2] = 0.17630838109670;
      refh[3] = 0.18734204669216;
      break;
    case 29:
      refh[0] = 0.0;
      refh[1] = -0.17543548739456;
      break;
    case 30:
      refh[0] = 0.0;
      refh[1] = -0.11109406510313;
      break;
    case 31:
      refh[0] = 0.0;
      refh[1] = -0.30692560016815;
      refh[2] = -0.12791615890916;
      break;
    case 32:
      refh[0] = 0.0;
      refh[1] = -0.15565415588148;
      refh[2] = -0.09744533013347;
      refh[3] = -0.02598788651562;
      refh[4] = -0.02082546400506;
      break;
    case 33:
      refh[0] = 0.0;
      refh[1] = -0.05341421485122;
      refh[2] = -0.04061733201402;
      refh[3] = -0.03040533259825;
      break;
    case 34:
      refh[0] = 0.0;
      refh[1] = 0.14680263633976;
      refh[2] = 0.15275393177853;
      break;
    case 35:
      refh[0] = 0.0;
      refh[1] = 0.26466311607527;
      break;
    case 36:
      refh[0] = 0.0;
      break;
    case 37:
      refh[0] = 0.0;
      refh[1] = -0.38168923853998;
      refh[2] = -0.61945426066521;
      break;
    case 38:
      refh[0] = 0.0;
      refh[1] = -0.35645143561998;
      refh[2] = -0.24590040509199;
      refh[3] = -0.58944749671323;
      break;
    case 39:
      refh[0] = 0.0;
      refh[1] = -0.12214207340151;
      refh[2] = -0.19122096573671;
      break;
    case 40:
      refh[0] = 0.0;
      refh[1] = -0.11717362032611;
      refh[2] = -0.30159667368243;
      refh[3] = -0.39900000000000;
      break;
    case 41:
      refh[0] = 0.0;
      refh[1] = -0.03613931721803;
      refh[2] = -0.09728395180973;
      refh[3] = -0.32760000000000;
      break;
    case 42:
      refh[0] = 0.0;
      refh[1] = -0.20855804369167;
      refh[2] = 0.37463101320270;
      refh[3] = 0.01088849508812;
      break;
    case 43:
      refh[0] = 0.0;
      refh[1] = -0.32295797028009;
      refh[2] = 0.02163058392055;
      break;
    case 44:
      refh[0] = 0.0;
      refh[1] = -0.16880624582411;
      refh[2] = 0.07067079599830;
      break;
    case 45:
      refh[0] = 0.0;
      refh[1] = -0.06766723747867;
      refh[2] = 0.02269366651357;
      refh[3] = 0.00189341152691;
      break;
    case 46:
      refh[0] = 0.0;
      refh[1] = -0.06864993376014;
      refh[2] = 0.17535498568548;
      break;
    case 47:
      refh[0] = 0.0;
      refh[1] = -0.07743633228111;
      break;
    case 48:
      refh[0] = 0.0;
      refh[1] = -0.16826792274681;
      break;
    case 49:
      refh[0] = 0.0;
      refh[1] = -0.34258129843160;
      refh[2] = -0.13326186209747;
      refh[3] = -0.06835793343128;
      break;
    case 50:
      refh[0] = 0.0;
      refh[1] = -0.11709105905770;
      refh[2] = -0.10012138923590;
      refh[3] = -0.02371941185707;
      refh[4] = -0.00510061648150;
      break;
    case 51:
      refh[0] = 0.0;
      refh[1] = -0.11698834506861;
      refh[2] = -0.09240034185831;
      refh[3] = -0.07823657941665;
      break;
    case 52:
      refh[0] = 0.0;
      refh[1] = 0.00520338227710;
      refh[2] = 0.00075938960030;
      break;
    case 53:
      refh[0] = 0.0;
      refh[1] = 0.18609272966361;
      break;
    case 54:
      refh[0] = 0.0;
      break;
    case 55:
      refh[0] = 0.0;
      refh[1] = -0.44073166763575;
      refh[2] = -0.68068249931321;
      break;
    case 56:
      refh[0] = 0.0;
      refh[1] = -0.39462085293344;
      refh[2] = -0.28358943523975;
      refh[3] = -0.62956728884112;
      break;
    case 57:
      refh[0] = 0.0;
      refh[1] = -0.07091364828303;
      refh[2] = -0.20649581935119;
      break;
    case 58:
      refh[0] = -0.18256074698566;
      break;
    case 59:
      refh[0] = 0.0;
      refh[1] = -0.18697130067341;
      break;
    case 60:
      refh[0] = 0.0;
      refh[1] = -0.19058554439983;
      break;
    case 61:
      refh[0] = 0.0;
      refh[1] = -0.19282831524410;
      break;
    case 62:
      refh[0] = 0.0;
      refh[1] = -0.19454581243048;
      break;
    case 63:
      refh[0] = 0.0;
      refh[1] = -0.07127505431420;
      break;
    case 64:
      refh[0] = 0.0;
      refh[1] = -0.19501047387777;
      break;
    case 65:
      refh[0] = 0.0;
      refh[1] = -0.20174382262063;
      break;
    case 66:
      refh[0] = 0.0;
      refh[1] = -0.20339885031380;
      break;
    case 67:
      refh[0] = 0.0;
      refh[1] = -0.20720952261691;
      break;
    case 68:
      refh[0] = 0.0;
      refh[1] = -0.20823422146908;
      break;
    case 69:
      refh[0] = 0.0;
      refh[1] = -0.21061010940412;
      break;
    case 70:
      refh[0] = 0.0;
      refh[1] = -0.20183381413024;
      break;
    case 71:
      refh[0] = 0.0;
      refh[1] = -0.21370999399902;
      break;
    case 72:
      refh[0] = 0.0;
      refh[1] = -0.05591493550625;
      refh[2] = -0.25794236832609;
      refh[3] = -0.33818965519407;
      break;
    case 73:
      refh[0] = 0.0;
      refh[1] = -0.08221450482839;
      refh[2] = -0.16758630536254;
      refh[3] = -0.36958292171049;
      break;
    case 74:
      refh[0] = 0.0;
      refh[1] = 0.16772725595657;
      refh[2] = 0.07390268846218;
      break;
    case 75:
      refh[0] = 0.0;
      refh[1] = -0.12381852305011;
      refh[2] = 0.10739644777728;
      break;
    case 76:
      refh[0] = 0.0;
      refh[1] = -0.19905415913991;
      refh[2] = 0.07497215257672;
      break;
    case 77:
      refh[0] = 0.0;
      refh[1] = -0.12646875304056;
      refh[2] = -0.02366554848104;
      refh[3] = -0.02022360384360;
      refh[4] = -0.01433797770974;
      break;
    case 78:
      refh[0] = 0.0;
      refh[1] = -0.03038537970592;
      refh[2] = 0.16414406238540;
      break;
    case 79:
      refh[0] = 0.0;
      refh[1] = -0.16039939518296;
      break;
    case 80:
      refh[0] = 0.0;
      refh[1] = 0.02160751342330;
      break;
    case 81:
      refh[0] = 0.0;
      refh[1] = -0.30915981534272;
      refh[2] = -0.04146409647185;
      refh[3] = 0.03239155712670;
      break;
    case 82:
      refh[0] = 0.0;
      refh[1] = -0.25127061741389;
      refh[2] = -0.18446979236239;
      refh[3] = -0.06678505471730;
      refh[4] = -0.00337046860960;
      break;
    case 83:
      refh[0] = 0.0;
      refh[1] = -0.16648528849972;
      refh[2] = -0.13791484089155;
      refh[3] = -0.11692906285535;
      break;
    case 84:
      refh[0] = 0.0;
      refh[1] = -0.13640050127609;
      refh[2] = -0.14076205965959;
      break;
    case 85:
      refh[0] = 0.0;
      refh[1] = -0.03800713219800;
      break;
    case 86:
      refh[0] = 0.0;
      break;
    case 112:
      refh[0] = 0.0;
      break;
    case 113:
      refh[0] = -0.10288904092396;
      break;
    case 114:
      refh[0] = 0.0;
      break;
    case 116:
      refh[0] = -0.02240831706236;
      break;
    case 115:
      refh[0] = -0.10991065534975;
      break;
    case 117:
      refh[0] = -0.07044735360796;
      break;
    case 118:
      refh[0] = 0.0;
      break;
  }
}
void RefSystem(int atmnr, std::vector<int> & refsys) {
  //reference system values for D4
  //vector not zeroed because we use RefN to access the vector's formal size
  switch (atmnr) {
    case 1:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 2:
      refsys[0] = 1;
      break;
    case 3:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 17;
      break;
    case 4:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 9;
      break;
    case 5:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 1;
      refsys[4] = 1;
      break;
    case 6:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 1;
      refsys[4] = 1;
      refsys[5] = 1;
      refsys[6] = 2;
      break;
    case 7:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 1;
      refsys[4] = 1;
      break;
    case 8:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 1;
      break;
    case 9:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 10:
      refsys[0] = 1;
      break;
    case 11:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 17;
      break;
    case 12:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 9;
      break;
    case 13:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 1;
      break;
    case 14:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 1;
      refsys[4] = 1;
      break;
    case 15:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 1;
      break;
    case 16:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      break;
    case 17:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 18:
      refsys[0] = 1;
      break;
    case 19:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 17;
      break;
    case 20:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 9;
      break;
    case 21:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 11;
      break;
    case 22:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 11;
      break;
    case 23:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 11;
      break;
    case 24:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 8;
      break;
    case 25:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      break;
    case 26:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      break;
    case 27:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 1;
      break;
    case 28:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 1;
      break;
    case 29:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 30:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 31:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      break;
    case 32:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 1;
      refsys[4] = 1;
      break;
    case 33:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 1;
      break;
    case 34:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      break;
    case 35:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 36:
      refsys[0] = 1;
      break;
    case 37:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 17;
      break;
    case 38:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 9;
      break;
    case 39:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      break;
    case 40:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 11;
      break;
    case 41:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 11;
      break;
    case 42:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 1;
      break;
    case 43:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      break;
    case 44:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      break;
    case 45:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 1;
      break;
    case 46:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      break;
    case 47:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 48:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 49:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 1;
      break;
    case 50:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 1;
      refsys[4] = 1;
      break;
    case 51:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 1;
      break;
    case 52:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      break;
    case 53:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 54:
      refsys[0] = 1;
      break;
    case 55:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 17;
      break;
    case 56:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 9;
      break;
    case 57:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      break;
    case 58:
      refsys[0] = 1;
      break;
    case 59:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 60:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 61:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 62:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 63:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 64:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 65:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 66:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 67:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 68:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 69:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 70:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 71:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 72:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 11;
      break;
    case 73:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 11;
      break;
    case 74:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      break;
    case 75:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      break;
    case 76:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      break;
    case 77:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 1;
      refsys[4] = 1;
      break;
    case 78:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      break;
    case 79:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 80:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 81:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 1;
      break;
    case 82:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 1;
      refsys[4] = 1;
      break;
    case 83:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      refsys[3] = 1;
      break;
    case 84:
      refsys[0] = 1;
      refsys[1] = 1;
      refsys[2] = 1;
      break;
    case 85:
      refsys[0] = 1;
      refsys[1] = 1;
      break;
    case 86:
      refsys[0] = 1;
      break;
    case 112:
      refsys[0] = 1;
      break;
    case 113:
      refsys[0] = 1;
      break;
    case 114:
      refsys[0] = 1;
      break;
    case 116:
      refsys[0] = 1;
      break;
    case 115:
      refsys[0] = 1;
      break;
    case 117:
      refsys[0] = 1;
      break;
    case 118:
      refsys[0] = 1;
      break;
  }
}
double SScale(int refsys) {
  //function returning the parameter sscale
  double sscale = 0.5;
  if (refsys == 2) {sscale = 0.0;}
  else if (refsys == 6) {sscale = 0.1666667;}
  else if (refsys == 7) {sscale = 1.0;}
  else if (refsys == 8) {sscale = 1.0;}
  return sscale;
}
double SecQ(int refsys) {
  //function returning SecQ
  return 0.0;
}
void SecaWeights(int refsys, std::vector<double> & secaiw) {
  //secaiw values for D4
  if ((refsys == 1)||(refsys == 2)) {
    secaiw[0]  = 5.4415160;
    secaiw[1]  = 5.3912720;
    secaiw[2]  = 5.2466780;
    secaiw[3]  = 4.7462570;
    secaiw[4]  = 4.1122050;
    secaiw[5]  = 3.4827990;
    secaiw[6]  = 2.9256260;
    secaiw[7]  = 2.4586020;
    secaiw[8]  = 2.0763900;
    secaiw[9]  = 1.7660350;
    secaiw[10] = 1.5138980;
    secaiw[11] = 1.3080740;
    secaiw[12] = 0.9987770;
    secaiw[13] = 0.7833600;
    secaiw[14] = 0.6286810;
    secaiw[15] = 0.5145050;
    secaiw[16] = 0.4281480;
    secaiw[17] = 0.2867670;
    secaiw[18] = 0.2047270;
    secaiw[19] = 0.1187560;
    secaiw[20] = 0.0772270;
    secaiw[21] = 0.0349350;
    secaiw[22] = 0.0197880;
  }
  else if (refsys == 6) {
    secaiw[0]  = 68.5832590;
    secaiw[1]  = 67.5115260;
    secaiw[2]  = 64.6123080;
    secaiw[3]  = 56.1286650;
    secaiw[4]  = 47.4318310;
    secaiw[5]  = 39.9459190;
    secaiw[6]  = 33.7814890;
    secaiw[7]  = 28.7553020;
    secaiw[8]  = 24.6561470;
    secaiw[9]  = 21.2992860;
    secaiw[10] = 18.5340330;
    secaiw[11] = 16.2406480;
    secaiw[12] = 12.7133690;
    secaiw[13] = 10.1832050;
    secaiw[14] = 8.3194640;
    secaiw[15] = 6.9133790;
    secaiw[16] = 5.8298100;
    secaiw[17] = 4.0106600;
    secaiw[18] = 2.9230920;
    secaiw[19] = 1.7494800;
    secaiw[20] = 1.1654830;
    secaiw[21] = 0.5523060;
    secaiw[22] = 0.3242020;
  }
  else if (refsys == 7) {
    secaiw[0]  = 13.8928580;
    secaiw[1]  = 13.7335660;
    secaiw[2]  = 13.2948950;
    secaiw[3]  = 11.9342710;
    secaiw[4]  = 10.4022050;
    secaiw[5]  = 8.9706190;
    secaiw[6]  = 7.7218140;
    secaiw[7]  = 6.6635680;
    secaiw[8]  = 5.7772340;
    secaiw[9]  = 5.0371340;
    secaiw[10] = 4.4181730;
    secaiw[11] = 3.8984410;
    secaiw[12] = 3.0872240;
    secaiw[13] = 2.4956330;
    secaiw[14] = 2.0539790;
    secaiw[15] = 1.7170460;
    secaiw[16] = 1.4549570;
    secaiw[17] = 1.0095450;
    secaiw[18] = 0.7395630;
    secaiw[19] = 0.4445600;
    secaiw[20] = 0.2961500;
    secaiw[21] = 0.1392520;
    secaiw[22] = 0.0809340;
  }
  else if (refsys == 8) {
    secaiw[0]  = 12.9390420;
    secaiw[1]  = 12.8215110;
    secaiw[2]  = 12.4887870;
    secaiw[3]  = 11.3861750;
    secaiw[4]  = 10.0547130;
    secaiw[5]  = 8.7593060;
    secaiw[6]  = 7.6055220;
    secaiw[7]  = 6.6166590;
    secaiw[8]  = 5.7823860;
    secaiw[9]  = 5.0816990;
    secaiw[10] = 4.4924840;
    secaiw[11] = 3.9949630;
    secaiw[12] = 3.2117190;
    secaiw[13] = 2.6333830;
    secaiw[14] = 2.1960400;
    secaiw[15] = 1.8580780;
    secaiw[16] = 1.5918740;
    secaiw[17] = 1.1306230;
    secaiw[18] = 0.8437390;
    secaiw[19] = 0.5212300;
    secaiw[20] = 0.3539970;
    secaiw[21] = 0.1716310;
    secaiw[22] = 0.1015580;
  }
  else if (refsys == 9) {
    secaiw[0]  = 10.3975870;
    secaiw[1]  = 10.3144230;
    secaiw[2]  = 10.0802310;
    secaiw[3]  = 9.3150480;
    secaiw[4]  = 8.4027020;
    secaiw[5]  = 7.5095120;
    secaiw[6]  = 6.6939110;
    secaiw[7]  = 5.9700690;
    secaiw[8]  = 5.3361820;
    secaiw[9]  = 4.7845300;
    secaiw[10] = 4.3054570;
    secaiw[11] = 3.8891980;
    secaiw[12] = 3.2100590;
    secaiw[13] = 2.6878840;
    secaiw[14] = 2.2798670;
    secaiw[15] = 1.9559530;
    secaiw[16] = 1.6950280;
    secaiw[17] = 1.2295900;
    secaiw[18] = 0.9306730;
    secaiw[19] = 0.5848620;
    secaiw[20] = 0.4008950;
    secaiw[21] = 0.1959850;
    secaiw[22] = 0.1160240;
  }
  else if (refsys == 10) {
    secaiw[0]  = 11.6125430;
    secaiw[1]  = 11.5410500;
    secaiw[2]  = 11.3332410;
    secaiw[3]  = 10.5895270;
    secaiw[4]  = 9.5886300;
    secaiw[5]  = 8.5228410;
    secaiw[6]  = 7.5106350;
    secaiw[7]  = 6.6043140;
    secaiw[8]  = 5.8164280;
    secaiw[9]  = 5.1406550;
    secaiw[10] = 4.5636850;
    secaiw[11] = 4.0709240;
    secaiw[12] = 3.2858130;
    secaiw[13] = 2.6994530;
    secaiw[14] = 2.2527560;
    secaiw[15] = 1.9059510;
    secaiw[16] = 1.6320270;
    secaiw[17] = 1.1565830;
    secaiw[18] = 0.8609730;
    secaiw[19] = 0.5294970;
    secaiw[20] = 0.3582750;
    secaiw[21] = 0.1725840;
    secaiw[22] = 0.1018070;
  }
  else if (refsys == 11) {
    secaiw[0]  = 8.1455493;
    secaiw[1]  = 8.1058081;
    secaiw[2]  = 7.9903083;
    secaiw[3]  = 7.5756433;
    secaiw[4]  = 7.0109396;
    secaiw[5]  = 6.3972855;
    secaiw[6]  = 5.7985381;
    secaiw[7]  = 5.2452249;
    secaiw[8]  = 4.7476246;
    secaiw[9]  = 4.3059510;
    secaiw[10] = 3.9161167;
    secaiw[11] = 3.5725710;
    secaiw[12] = 3.0018411;
    secaiw[13] = 2.5535671;
    secaiw[14] = 2.1967930;
    secaiw[15] = 1.9088609;
    secaiw[16] = 1.6733812;
    secaiw[17] = 1.2435728;
    secaiw[18] = 0.9590400;
    secaiw[19] = 0.6186562;
    secaiw[20] = 0.4313276;
    secaiw[21] = 0.2159017;
    secaiw[22] = 0.1291912;
  }
  else if (refsys == 17) {
    secaiw[0]  = 30.5968154;
    secaiw[1]  = 30.2907908;
    secaiw[2]  = 29.4281025;
    secaiw[3]  = 26.5588368;
    secaiw[4]  = 23.0571447;
    secaiw[5]  = 19.6470649;
    secaiw[6]  = 16.6420247;
    secaiw[7]  = 14.1130996;
    secaiw[8]  = 12.0267253;
    secaiw[9]  = 10.3165485;
    secaiw[10] = 8.9140109;
    secaiw[11] = 7.7590140;
    secaiw[12] = 6.0042243;
    secaiw[13] = 4.7671234;
    secaiw[14] = 3.8705622;
    secaiw[15] = 3.2039005;
    secaiw[16] = 2.6965878;
    secaiw[17] = 1.8588040;
    secaiw[18] = 1.3665825;
    secaiw[19] = 0.8411737;
    secaiw[20] = 0.5795824;
    secaiw[21] = 0.2981883;
    secaiw[22] = 0.1870572;
  }
}
void RefCovCN(int atmnr, std::vector<double> & refcovcn) {
  //reference covalent CN values for D4
  //vector not zeroed because we use RefN to access the vector's formal size
  switch (atmnr) {
    case 1:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.89422436546040;
      break;
    case 2:
      refcovcn[0] = 0.0;
      break;
    case 3:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.80196797025252;
      refcovcn[2] = 4.55700000000000;
      break;
    case 4:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.87483634146801;
      refcovcn[2] = 1.75692293756846;
      refcovcn[3] = 8.44900000000000;
      break;
    case 5:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.93017535174568;
      refcovcn[2] = 1.86318070483743;
      refcovcn[3] = 2.79150308430417;
      refcovcn[4] = 4.41327468359181;
      break;
    case 6:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.91894761153698;
      refcovcn[2] = 1.90786062124811;
      refcovcn[3] = 2.83122393221341;
      refcovcn[4] = 3.74870030209262;
      refcovcn[5] = 2.91832701804388;
      refcovcn[6] = 0.85560990866661;
      break;
    case 7:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.85827791612782;
      refcovcn[2] = 1.83985112339643;
      refcovcn[3] = 2.58039902538908;
      refcovcn[4] = 0.98085535554904;
      break;
    case 8:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.80416789582615;
      refcovcn[2] = 1.61123940672270;
      refcovcn[3] = 0.97984400975835;
      break;
    case 9:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.74076991274917;
      break;
    case 10:
      refcovcn[0] = 0.0;
      break;
    case 11:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.78077255822701;
      refcovcn[2] = 6.00100000000000;
      break;
    case 12:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.82437946770105;
      refcovcn[2] = 1.66972660811953;
      refcovcn[3] = 13.95600000000000;
      break;
    case 13:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.86567622356503;
      refcovcn[2] = 1.73262223552501;
      refcovcn[3] = 2.61541321571983;
      break;
    case 14:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.89153537011814;
      refcovcn[2] = 1.82343155163634;
      refcovcn[3] = 2.80029114356397;
      refcovcn[4] = 3.66675219599495;
      break;
    case 15:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.97498969761258;
      refcovcn[2] = 1.95559430261502;
      refcovcn[3] = 2.92666412421238;
      break;
    case 16:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.92198495968462;
      refcovcn[2] = 1.84453813661707;
      break;
    case 17:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.84439969796621;
      break;
    case 18:
      refcovcn[0] = 0.0;
      break;
    case 19:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.77389074718844;
      refcovcn[2] = 11.45000000000000;
      break;
    case 20:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.80173786145238;
      refcovcn[2] = 1.57758711825521;
      refcovcn[3] = 15.14000000000000;
      break;
    case 21:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.60644512750973;
      refcovcn[2] = 2.50274820582596;
      refcovcn[3] = 3.83500000000000;
      break;
    case 22:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.62278979228613;
      refcovcn[2] = 3.43281247525220;
      refcovcn[3] = 4.43849000000000;
      break;
    case 23:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.72211010302514;
      refcovcn[2] = 2.62026342999598;
      refcovcn[3] = 3.38800000000000;
      break;
    case 24:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.65041632145300;
      refcovcn[2] = 9.59333666980809;
      refcovcn[3] = 5.09941611449080;
      break;
    case 25:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.45861721080984;
      refcovcn[2] = 8.78012041793345;
      break;
    case 26:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.52985788829933;
      refcovcn[2] = 8.47495927159456;
      break;
    case 27:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.60281178943204;
      refcovcn[2] = 2.73732238386544;
      refcovcn[3] = 7.27710111777824;
      break;
    case 28:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.68368353980441;
      refcovcn[2] = 6.15501721386632;
      refcovcn[3] = 5.91752241519496;
      break;
    case 29:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.89816244646957;
      break;
    case 30:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.75336283145604;
      break;
    case 31:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.88820287155353;
      refcovcn[2] = 1.78721228422408;
      break;
    case 32:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.89941149965103;
      refcovcn[2] = 1.85525243419195;
      refcovcn[3] = 2.78448725734282;
      refcovcn[4] = 3.71823460680579;
      break;
    case 33:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.96810324138508;
      refcovcn[2] = 1.93762986034873;
      refcovcn[3] = 2.90829444959119;
      break;
    case 34:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.92195274896170;
      refcovcn[2] = 1.84493749809528;
      break;
    case 35:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.87141367228098;
      break;
    case 36:
      refcovcn[0] = 0.0;
      break;
    case 37:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.77155167960451;
      refcovcn[2] = 12.39000000000000;
      break;
    case 38:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.79294768050705;
      refcovcn[2] = 1.54805734413866;
      refcovcn[3] = 15.12100000000000;
      break;
    case 39:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.61670701598223;
      refcovcn[2] = 2.43996146811665;
      break;
    case 40:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.66263906536582;
      refcovcn[2] = 3.35991209283567;
      refcovcn[3] = 6.09300000000000;
      break;
    case 41:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.75136276843804;
      refcovcn[2] = 2.61865286219105;
      refcovcn[3] = 4.02900000000000;
      break;
    case 42:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.89510707051732;
      refcovcn[2] = 10.77850977361032;
      refcovcn[3] = 5.59646145748005;
      break;
    case 43:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.56460355314490;
      refcovcn[2] = 8.94426823783675;
      break;
    case 44:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.82293327778554;
      refcovcn[2] = 8.72336154530085;
      break;
    case 45:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.84241892648043;
      refcovcn[2] = 2.88090049147235;
      refcovcn[3] = 4.69798180984017;
      break;
    case 46:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.82751376419277;
      refcovcn[2] = 5.62156480548985;
      break;
    case 47:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.91215552676727;
      break;
    case 48:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.77531555102521;
      break;
    case 49:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.88719377220275;
      refcovcn[2] = 1.78495637384615;
      refcovcn[3] = 2.70420959057537;
      break;
    case 50:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.90042437628478;
      refcovcn[2] = 1.84687551160093;
      refcovcn[3] = 2.77047146659601;
      refcovcn[4] = 3.70466711294240;
      break;
    case 51:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.93580947053485;
      refcovcn[2] = 1.87397783659382;
      refcovcn[3] = 2.81393944504634;
      break;
    case 52:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.94900106848812;
      refcovcn[2] = 1.89954395376061;
      break;
    case 53:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.91239969640123;
      break;
    case 54:
      refcovcn[0] = 0.0;
      break;
    case 55:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.78131762864252;
      refcovcn[2] = 13.94000000000000;
      break;
    case 56:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.78501458609867;
      refcovcn[2] = 1.47819220471197;
      refcovcn[3] = 14.67700000000000;
      break;
    case 57:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.60493415835961;
      refcovcn[2] = 2.41715875994685;
      break;
    case 58:
      refcovcn[0] = 2.32446531201143;
      break;
    case 59:
      refcovcn[0] = 0.0;
      refcovcn[1] = 2.44976242184767;
      break;
    case 60:
      refcovcn[0] = 0.0;
      refcovcn[1] = 2.45621204848065;
      break;
    case 61:
      refcovcn[0] = 0.0;
      refcovcn[1] = 2.45653273553474;
      break;
    case 62:
      refcovcn[0] = 0.0;
      refcovcn[1] = 2.45472309204146;
      break;
    case 63:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.53272635693059;
      break;
    case 64:
      refcovcn[0] = 0.0;
      refcovcn[1] = 2.41754019361739;
      break;
    case 65:
      refcovcn[0] = 0.0;
      refcovcn[1] = 2.45190848719608;
      break;
    case 66:
      refcovcn[0] = 0.0;
      refcovcn[1] = 2.44535641509585;
      break;
    case 67:
      refcovcn[0] = 0.0;
      refcovcn[1] = 2.47312888776754;
      break;
    case 68:
      refcovcn[0] = 0.0;
      refcovcn[1] = 2.48043640078965;
      break;
    case 69:
      refcovcn[0] = 0.0;
      refcovcn[1] = 2.48120800603580;
      break;
    case 70:
      refcovcn[0] = 0.0;
      refcovcn[1] = 2.41884356027485;
      break;
    case 71:
      refcovcn[0] = 0.0;
      refcovcn[1] = 2.48673110207243;
      break;
    case 72:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.64596625526870;
      refcovcn[2] = 3.32095096313483;
      refcovcn[3] = 4.62300000000000;
      break;
    case 73:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.71671961147963;
      refcovcn[2] = 2.55699041057149;
      refcovcn[3] = 3.96300000000000;
      break;
    case 74:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.84852525804219;
      refcovcn[2] = 10.44998431099169;
      break;
    case 75:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.69673556648333;
      refcovcn[2] = 9.19053013143453;
      break;
    case 76:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.85014246784508;
      refcovcn[2] = 8.97974115543403;
      break;
    case 77:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.82924504495576;
      refcovcn[2] = 2.88910972828218;
      refcovcn[3] = 3.78801363268039;
      refcovcn[4] = 4.66121308138933;
      break;
    case 78:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.86072268763509;
      refcovcn[2] = 6.46855936294889;
      break;
    case 79:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.91830179564855;
      break;
    case 80:
      refcovcn[0] = 0.0;
      refcovcn[1] = 1.89515764632510;
      break;
    case 81:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.83218042800594;
      refcovcn[2] = 1.73451304607259;
      refcovcn[3] = 2.64455909521743;
      break;
    case 82:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.90267202371440;
      refcovcn[2] = 1.86326667138745;
      refcovcn[3] = 2.79951773292937;
      refcovcn[4] = 3.76444392234350;
      break;
    case 83:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.93840398774753;
      refcovcn[2] = 1.87829610127455;
      refcovcn[3] = 2.81907131341312;
      break;
    case 84:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.93512916979697;
      refcovcn[2] = 1.87110482058901;
      break;
    case 85:
      refcovcn[0] = 0.0;
      refcovcn[1] = 0.97711087656722;
      break;
    case 86:
      refcovcn[0] = 0.0;
      break;
    case 112:
      refcovcn[0] = 0.0;
      break;
    case 113:
      refcovcn[0] = 0.83118288370659;
      break;
    case 114:
      refcovcn[0] = 0.0;
      break;
    case 116:
      refcovcn[0] = 1.74812370366573;
      break;
    case 115:
      refcovcn[0] = 0.80392173226442;
      break;
    case 117:
      refcovcn[0] = 0.85819693463106;
      break;
    case 118:
      refcovcn[0] = 0.0;
      break;
  }
}
double HardnessRefSys(int refsys) {
  //hardness values for D4's reference systems
  double gamma = 0.0;
  if (refsys == 1) {gamma = 0.47259288;}
  else if (refsys == 2) {gamma = 0.92203391;}
  //else if (refsys == 6) {gamma = ;}
  //else if (refsys == 7) {gamma = ;}
  else if (refsys == 8) {gamma = 0.58691863;}
  else if (refsys == 9) {gamma = 0.66931351;}
  //else if (refsys == 10) {gamma = ;}
  else if (refsys == 11) {gamma = 0.17964105;}
  else if (refsys == 17) {gamma = 0.4311567;}
  return gamma;
}

double HardnessD4(int atmnr) {
  //hardness values for D4's atoms as in D. C. Ghosh, N. Islam, Int. J. Quant. Chem., 110(6), 1206, 2010
  //values in the paper multiplied by two because (IP - EA) = \frac{{d}^{2}E}{d{N}^{2}} but the hardness is \frac{1}{2}\frac{{d}^{2}E}{d{N}^{2}} (in Eh)
  double gamma = 0.0;
  switch (atmnr) {
    case 1: 
      gamma = 0.47259288;
      break;
    case 2: 
      gamma = 0.92203391;
      break;
    case 3: 
      gamma = 0.17452888;
      break;
    case 4: 
      gamma = 0.25700733;
      break;
    case 5: 
      gamma = 0.33949086;
      break;
    case 6: 
      gamma = 0.42195412;
      break;
    case 7: 
      gamma = 0.50438193;
      break;
    case 8: 
      gamma = 0.58691863;
      break;
    case 9: 
      gamma = 0.66931351;
      break;
    case 10: 
      gamma = 0.75191607;
      break;
    case 11: 
      gamma = 0.17964105;
      break;
    case 12: 
      gamma = 0.22157276;
      break;
    case 13: 
      gamma = 0.26348578;
      break;
    case 14: 
      gamma = 0.30539645;
      break;
    case 15: 
      gamma = 0.34734014;
      break;
    case 16: 
      gamma = 0.38924725;
      break;
    case 17: 
      gamma = 0.43115670;
      break;
    case 18: 
      gamma = 0.47308269;
      break;
    case 19: 
      gamma = 0.17105469;
      break;
    case 20: 
      gamma = 0.20276244;
      break;
    case 21: 
      gamma = 0.21007322;
      break;
    case 22: 
      gamma = 0.21739647;
      break;
    case 23: 
      gamma = 0.22471039;
      break;
    case 24: 
      gamma = 0.23201501;
      break;
    case 25: 
      gamma = 0.23933969;
      break;
    case 26: 
      gamma = 0.24665638;
      break;
    case 27: 
      gamma = 0.25398255;
      break;
    case 28: 
      gamma = 0.26128863;
      break;
    case 29: 
      gamma = 0.26859476;
      break;
    case 30: 
      gamma = 0.27592565;
      break;
    case 31: 
      gamma = 0.30762999;
      break;
    case 32: 
      gamma = 0.33931580;
      break;
    case 33: 
      gamma = 0.37235985;
      break;
    case 34: 
      gamma = 0.40273549;
      break;
    case 35: 
      gamma = 0.43445776;
      break;
    case 36: 
      gamma = 0.46611708;
      break;
    case 37: 
      gamma = 0.15585079;
      break;
    case 38: 
      gamma = 0.18649324;
      break;
    case 39: 
      gamma = 0.19356210;
      break;
    case 40: 
      gamma = 0.20063311;
      break;
    case 41: 
      gamma = 0.20770522;
      break;
    case 42: 
      gamma = 0.21477254;
      break;
    case 43: 
      gamma = 0.22184614;
      break;
    case 44: 
      gamma = 0.22891872;
      break;
    case 45: 
      gamma = 0.23598621;
      break;
    case 46: 
      gamma = 0.24305612;
      break;
    case 47: 
      gamma = 0.25013018;
      break;
    case 48: 
      gamma = 0.25719937;
      break;
    case 49: 
      gamma = 0.28784780;
      break;
    case 50: 
      gamma = 0.31848673;
      break;
    case 51: 
      gamma = 0.34912431;
      break;
    case 52: 
      gamma = 0.37976593;
      break;
    case 53: 
      gamma = 0.41040808;
      break;
    case 54: 
      gamma = 0.44105777;
      break;
    case 55: 
      gamma = 0.05019332;
      break;
    case 56: 
      gamma = 0.06762570;
      break;
    case 57: 
      gamma = 0.08504445;
      break;
    case 58: 
      gamma = 0.10247736;
      break;
    case 59: 
      gamma = 0.11991105;
      break;
    case 60: 
      gamma = 0.13732772;
      break;
    case 61: 
      gamma = 0.15476297;
      break;
    case 62: 
      gamma = 0.17218265;
      break;
    case 63: 
      gamma = 0.18961288;
      break;
    case 64: 
      gamma = 0.20704760;
      break;
    case 65: 
      gamma = 0.22446752;
      break;
    case 66: 
      gamma = 0.24189645;
      break;
    case 67: 
      gamma = 0.25932503;
      break;
    case 68: 
      gamma = 0.27676094;
      break;
    case 69: 
      gamma = 0.29418231;
      break;
    case 70: 
      gamma = 0.31159587;
      break;
    case 71: 
      gamma = 0.32902274;
      break;
    case 72: 
      gamma = 0.34592298;
      break;
    case 73: 
      gamma = 0.36388048;
      break;
    case 74: 
      gamma = 0.38130586;
      break;
    case 75: 
      gamma = 0.39877476;
      break;
    case 76: 
      gamma = 0.41614298;
      break;
    case 77: 
      gamma = 0.43364510;
      break;
    case 78: 
      gamma = 0.45104014;
      break;
    case 79: 
      gamma = 0.46848986;
      break;
    case 80: 
      gamma = 0.48584550;
      break;
    case 81: 
      gamma = 0.12526730;
      break;
    case 82: 
      gamma = 0.14268677;
      break;
    case 83: 
      gamma = 0.16011615;
      break;
    case 84: 
      gamma = 0.17755889;
      break;
    case 85: 
      gamma = 0.19497557;
      break;
    case 86: 
      gamma = 0.21240778;
      break;
    case 87: 
      gamma = 0.07263525;
      break;
    case 88: 
      gamma = 0.09422158;
      break;
    case 89: 
      gamma = 0.09920295;
      break;
    case 90: 
      gamma = 0.10418621;
      break;
    case 91: 
      gamma = 0.14235633;
      break;
    case 92: 
      gamma = 0.16394294;
      break;
    case 93: 
      gamma = 0.18551941;
      break;
    case 94: 
      gamma = 0.22370139;
      break;
  }
  return gamma;
}
double GA() {return 3.0;}
double GC() {return 2.0;} 
void RefSystem(int atmnr, int refsys, std::vector<double> & alphaiw) {
  //reference system values for D4
  //vector not zeroed because we use RefN to access the vector's formal size
  switch (atmnr) {
    case 1:
      if (refsys == 1) {
        alphaiw[0] = 5.0540160;
        alphaiw[1] = 4.9668210;
        alphaiw[2] = 4.7244390;
        alphaiw[3] = 3.9707860;
        alphaiw[4] = 3.1655030;
        alphaiw[5] = 2.4886460;
        alphaiw[6] = 1.9670460;
        alphaiw[7] = 1.5750840;
        alphaiw[8] = 1.2804290;
        alphaiw[9] = 1.0565330;
        alphaiw[10] = 0.8839250;
        alphaiw[11] = 0.7488080;
        alphaiw[12] = 0.5550060;
        alphaiw[13] = 0.4262920;
        alphaiw[14] = 0.3369390;
        alphaiw[15] = 0.2726050;
        alphaiw[16] = 0.2248580;
        alphaiw[17] = 0.1483870;
        alphaiw[18] = 0.1049920;
        alphaiw[19] = 0.0603100;
        alphaiw[20] = 0.0390250;
        alphaiw[21] = 0.0175590;
        alphaiw[22] = 0.0099250;
      }
      else if (refsys == 2) {
        alphaiw[0] = 5.4415160;
        alphaiw[1] = 5.3912720;
        alphaiw[2] = 5.2466780;
        alphaiw[3] = 4.7462570;
        alphaiw[4] = 4.1122050;
        alphaiw[5] = 3.4827990;
        alphaiw[6] = 2.9256260;
        alphaiw[7] = 2.4586020;
        alphaiw[8] = 2.0763900;
        alphaiw[9] = 1.7660350;
        alphaiw[10] = 1.5138980;
        alphaiw[11] = 1.3080740;
        alphaiw[12] = 0.9987770;
        alphaiw[13] = 0.7833600;
        alphaiw[14] = 0.6286810;
        alphaiw[15] = 0.5145050;
        alphaiw[16] = 0.4281480;
        alphaiw[17] = 0.2867670;
        alphaiw[18] = 0.2047270;
        alphaiw[19] = 0.1187560;
        alphaiw[20] = 0.0772270;
        alphaiw[21] = 0.0349350;
        alphaiw[22] = 0.0197880;
      }
      break;
    case 2:
      if (refsys == 1) {
        alphaiw[0] = 1.4546360;
        alphaiw[1] = 1.4502130;
        alphaiw[2] = 1.4371400;
        alphaiw[3] = 1.3875900;
        alphaiw[4] = 1.3135450;
        alphaiw[5] = 1.2243600;
        alphaiw[6] = 1.1286130;
        alphaiw[7] = 1.0327800;
        alphaiw[8] = 0.9410580;
        alphaiw[9] = 0.8557720;
        alphaiw[10] = 0.7779290;
        alphaiw[11] = 0.7077040;
        alphaiw[12] = 0.5886400;
        alphaiw[13] = 0.4940120;
        alphaiw[14] = 0.4187090;
        alphaiw[15] = 0.3583530;
        alphaiw[16] = 0.3095160;
        alphaiw[17] = 0.2223000;
        alphaiw[18] = 0.1665090;
        alphaiw[19] = 0.1026240;
        alphaiw[20] = 0.0691900;
        alphaiw[21] = 0.0327930;
        alphaiw[22] = 0.0189800;
      }
      break;
    case 3:
      if (refsys == 1) {
        alphaiw[0] = 149.6846930;
        alphaiw[1] = 101.0937750;
        alphaiw[2] = 51.6672830;
        alphaiw[3] = 18.0181080;
        alphaiw[4] = 8.9093460;
        alphaiw[5] = 5.3309030;
        alphaiw[6] = 3.5697190;
        alphaiw[7] = 2.5730730;
        alphaiw[8] = 1.9537540;
        alphaiw[9] = 1.5421740;
        alphaiw[10] = 1.2543870;
        alphaiw[11] = 1.0449370;
        alphaiw[12] = 0.7658330;
        alphaiw[13] = 0.5922860;
        alphaiw[14] = 0.4760820;
        alphaiw[15] = 0.3938140;
        alphaiw[16] = 0.3330090;
        alphaiw[17] = 0.2345470;
        alphaiw[18] = 0.1766850;
        alphaiw[19] = 0.1133040;
        alphaiw[20] = 0.0802740;
        alphaiw[21] = 0.0424220;
        alphaiw[22] = 0.0264160;
      }
      else if (refsys == 2) {
        alphaiw[0] = 27.7934790;
        alphaiw[1] = 25.9220070;
        alphaiw[2] = 21.8214140;
        alphaiw[3] = 14.1717810;
        alphaiw[4] = 9.4824880;
        alphaiw[5] = 6.7129090;
        alphaiw[6] = 4.9857250;
        alphaiw[7] = 3.8465020;
        alphaiw[8] = 3.0589760;
        alphaiw[9] = 2.4932250;
        alphaiw[10] = 2.0736610;
        alphaiw[11] = 1.7541010;
        alphaiw[12] = 1.3073430;
        alphaiw[13] = 1.0165020;
        alphaiw[14] = 0.8161580;
        alphaiw[15] = 0.6719100;
        alphaiw[16] = 0.5643080;
        alphaiw[17] = 0.3893270;
        alphaiw[18] = 0.2870170;
        alphaiw[19] = 0.1768550;
        alphaiw[20] = 0.1211260;
        alphaiw[21] = 0.0602840;
        alphaiw[22] = 0.0362580;
      }
      else if (refsys == 3) {
        alphaiw[0] = 84.6459010;
        alphaiw[1] = 83.6841238;
        alphaiw[2] = 80.9655079;
        alphaiw[3] = 72.0284169;
        alphaiw[4] = 61.4818132;
        alphaiw[5] = 51.5922762;
        alphaiw[6] = 43.1612588;
        alphaiw[7] = 36.2527416;
        alphaiw[8] = 30.6710866;
        alphaiw[9] = 26.1700904;
        alphaiw[10] = 22.5260442;
        alphaiw[11] = 19.5557200;
        alphaiw[12] = 15.0928030;
        alphaiw[13] = 11.9795535;
        alphaiw[14] = 9.7376006;
        alphaiw[15] = 8.0762262;
        alphaiw[16] = 6.8136584;
        alphaiw[17] = 4.7267579;
        alphaiw[18] = 3.4948520;
        alphaiw[19] = 2.1681444;
        alphaiw[20] = 1.4997611;
        alphaiw[21] = 0.7718181;
        alphaiw[22] = 0.4815951;
      }
      break;
    case 4:
      if (refsys == 1) {
        alphaiw[0] = 42.8790460;
        alphaiw[1] = 39.9295910;
        alphaiw[2] = 33.1319060;
        alphaiw[3] = 19.8584460;
        alphaiw[4] = 12.0562070;
        alphaiw[5] = 7.8731240;
        alphaiw[6] = 5.5003700;
        alphaiw[7] = 4.0522140;
        alphaiw[8] = 3.1097180;
        alphaiw[9] = 2.4633930;
        alphaiw[10] = 2.0011630;
        alphaiw[11] = 1.6591860;
        alphaiw[12] = 1.1965720;
        alphaiw[13] = 0.9061710;
        alphaiw[14] = 0.7119670;
        alphaiw[15] = 0.5756820;
        alphaiw[16] = 0.4763260;
        alphaiw[17] = 0.3199380;
        alphaiw[18] = 0.2321350;
        alphaiw[19] = 0.1411470;
        alphaiw[20] = 0.0965570;
        alphaiw[21] = 0.0488560;
        alphaiw[22] = 0.0302000;
      }
      else if (refsys == 2) {
        alphaiw[0] = 33.9191750;
        alphaiw[1] = 30.6229690;
        alphaiw[2] = 25.1809150;
        alphaiw[3] = 17.1197860;
        alphaiw[4] = 12.0494810;
        alphaiw[5] = 8.7823960;
        alphaiw[6] = 6.6233290;
        alphaiw[7] = 5.1483500;
        alphaiw[8] = 4.1059930;
        alphaiw[9] = 3.3462160;
        alphaiw[10] = 2.7771820;
        alphaiw[11] = 2.3408840;
        alphaiw[12] = 1.7276510;
        alphaiw[13] = 1.3276250;
        alphaiw[14] = 1.0528460;
        alphaiw[15] = 0.8562260;
        alphaiw[16] = 0.7107890;
        alphaiw[17] = 0.4781520;
        alphaiw[18] = 0.3456310;
        alphaiw[19] = 0.2073370;
        alphaiw[20] = 0.1396630;
        alphaiw[21] = 0.0681800;
        alphaiw[22] = 0.0409330;
      }
      else if (refsys == 3) {
        alphaiw[0] = 20.3763510;
        alphaiw[1] = 19.9881660;
        alphaiw[2] = 18.9207290;
        alphaiw[3] = 15.6997430;
        alphaiw[4] = 12.3941440;
        alphaiw[5] = 9.6994520;
        alphaiw[6] = 7.6613340;
        alphaiw[7] = 6.1448630;
        alphaiw[8] = 5.0098820;
        alphaiw[9] = 4.1485350;
        alphaiw[10] = 3.4841730;
        alphaiw[11] = 2.9633890;
        alphaiw[12] = 2.2143260;
        alphaiw[13] = 1.7146750;
        alphaiw[14] = 1.3662970;
        alphaiw[15] = 1.1143890;
        alphaiw[16] = 0.9266450;
        alphaiw[17] = 0.6239030;
        alphaiw[18] = 0.4502460;
        alphaiw[19] = 0.2684610;
        alphaiw[20] = 0.1795900;
        alphaiw[21] = 0.0862930;
        alphaiw[22] = 0.0511440;
      }
      else if (refsys == 4) {
        alphaiw[0] = 50.8062040;
        alphaiw[1] = 50.0592716;
        alphaiw[2] = 48.0905154;
        alphaiw[3] = 42.5028663;
        alphaiw[4] = 36.7056665;
        alphaiw[5] = 31.5232793;
        alphaiw[6] = 27.1066190;
        alphaiw[7] = 23.4118521;
        alphaiw[8] = 20.3403941;
        alphaiw[9] = 17.7870106;
        alphaiw[10] = 15.6568911;
        alphaiw[11] = 13.8703586;
        alphaiw[12] = 11.0817107;
        alphaiw[13] = 9.0435022;
        alphaiw[14] = 7.5156407;
        alphaiw[15] = 6.3438199;
        alphaiw[16] = 5.4266858;
        alphaiw[17] = 3.8505033;
        alphaiw[18] = 2.8781608;
        alphaiw[19] = 1.9058183;
        alphaiw[20] = 1.2248615;
        alphaiw[21] = 0.6029055;
        alphaiw[22] = 0.3602096;
      }
      break;
    case 5:
      if (refsys == 1) {
        alphaiw[0] = 22.1068880;
        alphaiw[1] = 21.3442260;
        alphaiw[2] = 19.4317910;
        alphaiw[3] = 14.6694210;
        alphaiw[4] = 10.6975470;
        alphaiw[5] = 7.9023200;
        alphaiw[6] = 5.9922100;
        alphaiw[7] = 4.6684490;
        alphaiw[8] = 3.7269880;
        alphaiw[9] = 3.0389440;
        alphaiw[10] = 2.5231990;
        alphaiw[11] = 2.1277390;
        alphaiw[12] = 1.5720620;
        alphaiw[13] = 1.2095680;
        alphaiw[14] = 0.9603270;
        alphaiw[15] = 0.7816570;
        alphaiw[16] = 0.6491900;
        alphaiw[17] = 0.4364560;
        alphaiw[18] = 0.3146920;
        alphaiw[19] = 0.1873730;
        alphaiw[20] = 0.1252700;
        alphaiw[21] = 0.0602010;
        alphaiw[22] = 0.0356010;
      }
      else if (refsys == 2) {
        alphaiw[0] = 20.5521070;
        alphaiw[1] = 19.9101690;
        alphaiw[2] = 18.3552210;
        alphaiw[3] = 14.5802680;
        alphaiw[4] = 11.3056180;
        alphaiw[5] = 8.8165140;
        alphaiw[6] = 6.9794390;
        alphaiw[7] = 5.6206180;
        alphaiw[8] = 4.6022690;
        alphaiw[9] = 3.8263850;
        alphaiw[10] = 3.2251530;
        alphaiw[11] = 2.7516720;
        alphaiw[12] = 2.0664510;
        alphaiw[13] = 1.6059710;
        alphaiw[14] = 1.2828710;
        alphaiw[15] = 1.0479810;
        alphaiw[16] = 0.8721210;
        alphaiw[17] = 0.5869300;
        alphaiw[18] = 0.4224710;
        alphaiw[19] = 0.2500750;
        alphaiw[20] = 0.1661340;
        alphaiw[21] = 0.0788140;
        alphaiw[22] = 0.0462310;
      }
      else if (refsys == 3) {
        alphaiw[0] = 20.2806350;
        alphaiw[1] = 19.1635070;
        alphaiw[2] = 17.8009360;
        alphaiw[3] = 14.8840980;
        alphaiw[4] = 12.0983160;
        alphaiw[5] = 9.7836060;
        alphaiw[6] = 7.9569310;
        alphaiw[7] = 6.5370310;
        alphaiw[8] = 5.4328530;
        alphaiw[9] = 4.5678260;
        alphaiw[10] = 3.8830910;
        alphaiw[11] = 3.3348630;
        alphaiw[12] = 2.5273350;
        alphaiw[13] = 1.9752950;
        alphaiw[14] = 1.5835380;
        alphaiw[15] = 1.2965230;
        alphaiw[16] = 1.0804740;
        alphaiw[17] = 0.7281530;
        alphaiw[18] = 0.5240250;
        alphaiw[19] = 0.3095240;
        alphaiw[20] = 0.2050450;
        alphaiw[21] = 0.0966460;
        alphaiw[22] = 0.0564440;
      }
      else if (refsys == 4) {
        alphaiw[0] = 17.3839000;
        alphaiw[1] = 17.1992340;
        alphaiw[2] = 16.6725580;
        alphaiw[3] = 14.8986510;
        alphaiw[4] = 12.7407520;
        alphaiw[5] = 10.6754850;
        alphaiw[6] = 8.8969080;
        alphaiw[7] = 7.4349140;
        alphaiw[8] = 6.2546920;
        alphaiw[9] = 5.3056760;
        alphaiw[10] = 4.5401840;
        alphaiw[11] = 3.9186610;
        alphaiw[12] = 2.9899200;
        alphaiw[13] = 2.3464510;
        alphaiw[14] = 1.8858730;
        alphaiw[15] = 1.5464990;
        alphaiw[16] = 1.2900300;
        alphaiw[17] = 0.8701190;
        alphaiw[18] = 0.6260260;
        alphaiw[19] = 0.3691130;
        alphaiw[20] = 0.2439620;
        alphaiw[21] = 0.1143850;
        alphaiw[22] = 0.0665550;
      }
      else if (refsys == 5) {
        alphaiw[0] = 33.3389310;
        alphaiw[1] = 33.0011270;
        alphaiw[2] = 32.0357710;
        alphaiw[3] = 28.7639060;
        alphaiw[4] = 24.7422160;
        alphaiw[5] = 20.8503600;
        alphaiw[6] = 17.4641170;
        alphaiw[7] = 14.6555110;
        alphaiw[8] = 12.3709610;
        alphaiw[9] = 10.5223950;
        alphaiw[10] = 9.0236300;
        alphaiw[11] = 7.8016400;
        alphaiw[12] = 5.9670990;
        alphaiw[13] = 4.6901690;
        alphaiw[14] = 3.7733460;
        alphaiw[15] = 3.0963780;
        alphaiw[16] = 2.5840490;
        alphaiw[17] = 1.7440470;
        alphaiw[18] = 1.2552220;
        alphaiw[19] = 0.7404570;
        alphaiw[20] = 0.4896440;
        alphaiw[21] = 0.2299100;
        alphaiw[22] = 0.1339700;
      }
      break;
    case 6:
      if (refsys == 1) {
        alphaiw[0] = 12.3555340;
        alphaiw[1] = 12.1386640;
        alphaiw[2] = 11.5492610;
        alphaiw[3] = 9.7999490;
        alphaiw[4] = 7.9916820;
        alphaiw[5] = 6.4691120;
        alphaiw[6] = 5.2715550;
        alphaiw[7] = 4.3463700;
        alphaiw[8] = 3.6298580;
        alphaiw[9] = 3.0691570;
        alphaiw[10] = 2.6246820;
        alphaiw[11] = 2.2676780;
        alphaiw[12] = 1.7384300;
        alphaiw[13] = 1.3729830;
        alphaiw[14] = 1.1110210;
        alphaiw[15] = 0.9172590;
        alphaiw[16] = 0.7701000;
        alphaiw[17] = 0.5269110;
        alphaiw[18] = 0.3835030;
        alphaiw[19] = 0.2298730;
        alphaiw[20] = 0.1535810;
        alphaiw[21] = 0.0732580;
        alphaiw[22] = 0.0431920;
      }
      else if (refsys == 2) {
        alphaiw[0] = 14.5482810;
        alphaiw[1] = 14.2754640;
        alphaiw[2] = 13.3786710;
        alphaiw[3] = 11.3492310;
        alphaiw[4] = 9.4185310;
        alphaiw[5] = 7.7719120;
        alphaiw[6] = 6.4361850;
        alphaiw[7] = 5.3728750;
        alphaiw[8] = 4.5287720;
        alphaiw[9] = 3.8553340;
        alphaiw[10] = 3.3134960;
        alphaiw[11] = 2.8732790;
        alphaiw[12] = 2.2127670;
        alphaiw[13] = 1.7514880;
        alphaiw[14] = 1.4184590;
        alphaiw[15] = 1.1710210;
        alphaiw[16] = 0.9825780;
        alphaiw[17] = 0.6705550;
        alphaiw[18] = 0.4865370;
        alphaiw[19] = 0.2899040;
        alphaiw[20] = 0.1927160;
        alphaiw[21] = 0.0910580;
        alphaiw[22] = 0.0533330;
      }
      else if (refsys == 3) {
        alphaiw[0] = 23.0116230;
        alphaiw[1] = 22.7312670;
        alphaiw[2] = 21.9447640;
        alphaiw[3] = 19.4160980;
        alphaiw[4] = 16.5257670;
        alphaiw[5] = 13.8835800;
        alphaiw[6] = 11.6644030;
        alphaiw[7] = 9.8566400;
        alphaiw[8] = 8.3954220;
        alphaiw[9] = 7.2118160;
        alphaiw[10] = 6.2469570;
        alphaiw[11] = 5.4540580;
        alphaiw[12] = 4.2475560;
        alphaiw[13] = 3.3914730;
        alphaiw[14] = 2.7655050;
        alphaiw[15] = 2.2955930;
        alphaiw[16] = 1.9346700;
        alphaiw[17] = 1.3305920;
        alphaiw[18] = 0.9701360;
        alphaiw[19] = 0.5811850;
        alphaiw[20] = 0.3874320;
        alphaiw[21] = 0.1836940;
        alphaiw[22] = 0.1078100;
      }
      else if (refsys == 4) {
        alphaiw[0] = 27.4537640;
        alphaiw[1] = 27.0921900;
        alphaiw[2] = 26.0934140;
        alphaiw[3] = 23.0050050;
        alphaiw[4] = 19.6170750;
        alphaiw[5] = 16.5717410;
        alphaiw[6] = 14.0107840;
        alphaiw[7] = 11.9057630;
        alphaiw[8] = 10.1859920;
        alphaiw[9] = 8.7791490;
        alphaiw[10] = 7.6227700;
        alphaiw[11] = 6.6661200;
        alphaiw[12] = 5.2000230;
        alphaiw[13] = 4.1529630;
        alphaiw[14] = 3.3845150;
        alphaiw[15] = 2.8065480;
        alphaiw[16] = 2.3623130;
        alphaiw[17] = 1.6190990;
        alphaiw[18] = 1.1765830;
        alphaiw[19] = 0.7009620;
        alphaiw[20] = 0.4652560;
        alphaiw[21] = 0.2188710;
        alphaiw[22] = 0.1277540;
      }
      else if (refsys == 5) {
        alphaiw[0] = 28.3352760;
        alphaiw[1] = 28.0932470;
        alphaiw[2] = 27.3990360;
        alphaiw[3] = 25.0151020;
        alphaiw[4] = 22.0076950;
        alphaiw[5] = 18.9989890;
        alphaiw[6] = 16.2867160;
        alphaiw[7] = 13.9588710;
        alphaiw[8] = 12.0051730;
        alphaiw[9] = 10.3795410;
        alphaiw[10] = 9.0286370;
        alphaiw[11] = 7.9030880;
        alphaiw[12] = 6.1676970;
        alphaiw[13] = 4.9231280;
        alphaiw[14] = 4.0083160;
        alphaiw[15] = 3.3200870;
        alphaiw[16] = 2.7913060;
        alphaiw[17] = 1.9077980;
        alphaiw[18] = 1.3830030;
        alphaiw[19] = 0.8207370;
        alphaiw[20] = 0.5431290;
        alphaiw[21] = 0.2541200;
        alphaiw[22] = 0.1477570;
      }
      else if (refsys == 6) {
        alphaiw[0] = 68.5832590;
        alphaiw[1] = 67.5115260;
        alphaiw[2] = 64.6123080;
        alphaiw[3] = 56.1286650;
        alphaiw[4] = 47.4318310;
        alphaiw[5] = 39.9459190;
        alphaiw[6] = 33.7814890;
        alphaiw[7] = 28.7553020;
        alphaiw[8] = 24.6561470;
        alphaiw[9] = 21.2992860;
        alphaiw[10] = 18.5340330;
        alphaiw[11] = 16.2406480;
        alphaiw[12] = 12.7133690;
        alphaiw[13] = 10.1832050;
        alphaiw[14] = 8.3194640;
        alphaiw[15] = 6.9133790;
        alphaiw[16] = 5.8298100;
        alphaiw[17] = 4.0106600;
        alphaiw[18] = 2.9230920;
        alphaiw[19] = 1.7494800;
        alphaiw[20] = 1.1654830;
        alphaiw[21] = 0.5523060;
        alphaiw[22] = 0.3242020;
      }
      else if (refsys == 7) {
        alphaiw[0] = 7.7402500;
        alphaiw[1] = 7.6643000;
        alphaiw[2] = 7.4486700;
        alphaiw[3] = 6.7286500;
        alphaiw[4] = 5.8533600;
        alphaiw[5] = 5.0045500;
        alphaiw[6] = 4.2585700;
        alphaiw[7] = 3.6316200;
        alphaiw[8] = 3.1142900;
        alphaiw[9] = 2.6894300;
        alphaiw[10] = 2.3397600;
        alphaiw[11] = 2.0503600;
        alphaiw[12] = 1.6066900;
        alphaiw[13] = 1.2894400;
        alphaiw[14] = 1.0561100;
        alphaiw[15] = 0.8801010;
        alphaiw[16] = 0.7443600;
        alphaiw[17] = 0.5158280;
        alphaiw[18] = 0.3784030;
        alphaiw[19] = 0.2287990;
        alphaiw[20] = 0.1535500;
        alphaiw[21] = 0.0736385;
        alphaiw[22] = 0.0435460;
      }
      break;
    case 7:
      if (refsys == 1) {
        alphaiw[0] = 7.5385090;
        alphaiw[1] = 7.4657840;
        alphaiw[2] = 7.2592230;
        alphaiw[3] = 6.5708420;
        alphaiw[4] = 5.7412090;
        alphaiw[5] = 4.9445040;
        alphaiw[6] = 4.2474110;
        alphaiw[7] = 3.6605820;
        alphaiw[8] = 3.1734720;
        alphaiw[9] = 2.7700380;
        alphaiw[10] = 2.4347170;
        alphaiw[11] = 2.1542670;
        alphaiw[12] = 1.7174610;
        alphaiw[13] = 1.3984350;
        alphaiw[14] = 1.1591760;
        alphaiw[15] = 0.9755810;
        alphaiw[16] = 0.8318860;
        alphaiw[17] = 0.5852590;
        alphaiw[18] = 0.4337270;
        alphaiw[19] = 0.2654650;
        alphaiw[20] = 0.1792040;
        alphaiw[21] = 0.0860960;
        alphaiw[22] = 0.0507150;
      }
      else if (refsys == 2) {
        alphaiw[0] = 9.8340890;
        alphaiw[1] = 9.6980440;
        alphaiw[2] = 9.3521500;
        alphaiw[3] = 8.3811310;
        alphaiw[4] = 7.3184230;
        alphaiw[5] = 6.3183850;
        alphaiw[6] = 5.4407240;
        alphaiw[7] = 4.6956550;
        alphaiw[8] = 4.0720980;
        alphaiw[9] = 3.5523310;
        alphaiw[10] = 3.1184230;
        alphaiw[11] = 2.7545880;
        alphaiw[12] = 2.1872670;
        alphaiw[13] = 1.7734490;
        alphaiw[14] = 1.4640490;
        alphaiw[15] = 1.2275070;
        alphaiw[16] = 1.0430710;
        alphaiw[17] = 0.7283800;
        alphaiw[18] = 0.5365400;
        alphaiw[19] = 0.3254160;
        alphaiw[20] = 0.2182610;
        alphaiw[21] = 0.1038080;
        alphaiw[22] = 0.0607780;
      }
      else if (refsys == 3) {
        alphaiw[0] = 18.5281530;
        alphaiw[1] = 18.3270930;
        alphaiw[2] = 17.7675350;
        alphaiw[3] = 15.9940890;
        alphaiw[4] = 13.9674690;
        alphaiw[5] = 12.0720960;
        alphaiw[6] = 10.4237930;
        alphaiw[7] = 9.0298300;
        alphaiw[8] = 7.8623050;
        alphaiw[9] = 6.8857300;
        alphaiw[10] = 6.0666250;
        alphaiw[11] = 5.3763210;
        alphaiw[12] = 4.2923510;
        alphaiw[13] = 3.4950320;
        alphaiw[14] = 2.8948440;
        alphaiw[15] = 2.4334700;
        alphaiw[16] = 2.0721210;
        alphaiw[17] = 1.4520950;
        alphaiw[18] = 1.0718490;
        alphaiw[19] = 0.6513790;
        alphaiw[20] = 0.4372510;
        alphaiw[21] = 0.2081480;
        alphaiw[22] = 0.1219450;
      }
      else if (refsys == 4) {
        alphaiw[0] = 13.8928580;
        alphaiw[1] = 13.7335660;
        alphaiw[2] = 13.2948950;
        alphaiw[3] = 11.9342710;
        alphaiw[4] = 10.4022050;
        alphaiw[5] = 8.9706190;
        alphaiw[6] = 7.7218140;
        alphaiw[7] = 6.6635680;
        alphaiw[8] = 5.7772340;
        alphaiw[9] = 5.0371340;
        alphaiw[10] = 4.4181730;
        alphaiw[11] = 3.8984410;
        alphaiw[12] = 3.0872240;
        alphaiw[13] = 2.4956330;
        alphaiw[14] = 2.0539790;
        alphaiw[15] = 1.7170460;
        alphaiw[16] = 1.4549570;
        alphaiw[17] = 1.0095450;
        alphaiw[18] = 0.7395630;
        alphaiw[19] = 0.4445600;
        alphaiw[20] = 0.2961500;
        alphaiw[21] = 0.1392520;
        alphaiw[22] = 0.0809340;
      }
      else if (refsys == 5) {
        alphaiw[0] = 11.6125430;
        alphaiw[1] = 11.5410500;
        alphaiw[2] = 11.3332410;
        alphaiw[3] = 10.5895270;
        alphaiw[4] = 9.5886300;
        alphaiw[5] = 8.5228410;
        alphaiw[6] = 7.5106350;
        alphaiw[7] = 6.6043140;
        alphaiw[8] = 5.8164280;
        alphaiw[9] = 5.1406550;
        alphaiw[10] = 4.5636850;
        alphaiw[11] = 4.0709240;
        alphaiw[12] = 3.2858130;
        alphaiw[13] = 2.6994530;
        alphaiw[14] = 2.2527560;
        alphaiw[15] = 1.9059510;
        alphaiw[16] = 1.6320270;
        alphaiw[17] = 1.1565830;
        alphaiw[18] = 0.8609730;
        alphaiw[19] = 0.5294970;
        alphaiw[20] = 0.3582750;
        alphaiw[21] = 0.1725840;
        alphaiw[22] = 0.1018070;
      }
      break;
    case 8:
      if (refsys == 1) {
        alphaiw[0] = 5.1967090;
        alphaiw[1] = 5.1644300;
        alphaiw[2] = 5.0715690;
        alphaiw[3] = 4.7484500;
        alphaiw[4] = 4.3273530;
        alphaiw[5] = 3.8854260;
        alphaiw[6] = 3.4648010;
        alphaiw[7] = 3.0838640;
        alphaiw[8] = 2.7478290;
        alphaiw[9] = 2.4553310;
        alphaiw[10] = 2.2021600;
        alphaiw[11] = 1.9832410;
        alphaiw[12] = 1.6287590;
        alphaiw[13] = 1.3587170;
        alphaiw[14] = 1.1492470;
        alphaiw[15] = 0.9838790;
        alphaiw[16] = 0.8512370;
        alphaiw[17] = 0.6158310;
        alphaiw[18] = 0.4654260;
        alphaiw[19] = 0.2921270;
        alphaiw[20] = 0.2001820;
        alphaiw[21] = 0.0978590;
        alphaiw[22] = 0.0579230;
      }
      else if (refsys == 2) {
        alphaiw[0] = 7.2941570;
        alphaiw[1] = 7.2356790;
        alphaiw[2] = 7.0728280;
        alphaiw[3] = 6.5406890;
        alphaiw[4] = 5.8872890;
        alphaiw[5] = 5.2288430;
        alphaiw[6] = 4.6198180;
        alphaiw[7] = 4.0794700;
        alphaiw[8] = 3.6098730;
        alphaiw[9] = 3.2056430;
        alphaiw[10] = 2.8587920;
        alphaiw[11] = 2.5610180;
        alphaiw[12] = 2.0831650;
        alphaiw[13] = 1.7231230;
        alphaiw[14] = 1.4466600;
        alphaiw[15] = 1.2304460;
        alphaiw[16] = 1.0585030;
        alphaiw[17] = 0.7570550;
        alphaiw[18] = 0.5672680;
        alphaiw[19] = 0.3517810;
        alphaiw[20] = 0.2391180;
        alphaiw[21] = 0.1155270;
        alphaiw[22] = 0.0679560;
      }
      else if (refsys == 3) {
        alphaiw[0] = 9.3715840;
        alphaiw[1] = 9.2917670;
        alphaiw[2] = 9.0663960;
        alphaiw[3] = 8.3227770;
        alphaiw[4] = 7.4235050;
        alphaiw[5] = 6.5392070;
        alphaiw[6] = 5.7386200;
        alphaiw[7] = 5.0398020;
        alphaiw[8] = 4.4396540;
        alphaiw[9] = 3.9274640;
        alphaiw[10] = 3.4907530;
        alphaiw[11] = 3.1176580;
        alphaiw[12] = 2.5222470;
        alphaiw[13] = 2.0764190;
        alphaiw[14] = 1.7359450;
        alphaiw[15] = 1.4709930;
        alphaiw[16] = 1.2612520;
        alphaiw[17] = 0.8959550;
        alphaiw[18] = 0.6678320;
        alphaiw[19] = 0.4109820;
        alphaiw[20] = 0.2778660;
        alphaiw[21] = 0.1331600;
        alphaiw[22] = 0.0779780;
      }
      else if (refsys == 4) {
        alphaiw[0] = 10.3975870;
        alphaiw[1] = 10.3144230;
        alphaiw[2] = 10.0802310;
        alphaiw[3] = 9.3150480;
        alphaiw[4] = 8.4027020;
        alphaiw[5] = 7.5095120;
        alphaiw[6] = 6.6939110;
        alphaiw[7] = 5.9700690;
        alphaiw[8] = 5.3361820;
        alphaiw[9] = 4.7845300;
        alphaiw[10] = 4.3054570;
        alphaiw[11] = 3.8891980;
        alphaiw[12] = 3.2100590;
        alphaiw[13] = 2.6878840;
        alphaiw[14] = 2.2798670;
        alphaiw[15] = 1.9559530;
        alphaiw[16] = 1.6950280;
        alphaiw[17] = 1.2295900;
        alphaiw[18] = 0.9306730;
        alphaiw[19] = 0.5848620;
        alphaiw[20] = 0.4008950;
        alphaiw[21] = 0.1959850;
        alphaiw[22] = 0.1160240;
      }
      break;
    case 9:
      if (refsys == 1) {
        alphaiw[0] = 3.6338390;
        alphaiw[1] = 3.6196960;
        alphaiw[2] = 3.5783940;
        alphaiw[3] = 3.4279380;
        alphaiw[4] = 3.2173720;
        alphaiw[5] = 2.9799290;
        alphaiw[6] = 2.7383490;
        alphaiw[7] = 2.5057910;
        alphaiw[8] = 2.2889990;
        alphaiw[9] = 2.0908140;
        alphaiw[10] = 1.9117910;
        alphaiw[11] = 1.7512000;
        alphaiw[12] = 1.4795260;
        alphaiw[13] = 1.2627740;
        alphaiw[14] = 1.0886880;
        alphaiw[15] = 0.9474060;
        alphaiw[16] = 0.8314380;
        alphaiw[17] = 0.6189740;
        alphaiw[18] = 0.4778550;
        alphaiw[19] = 0.3086410;
        alphaiw[20] = 0.2153540;
        alphaiw[21] = 0.1079130;
        alphaiw[22] = 0.0645990;
      }
      else if (refsys == 2) {
        alphaiw[0] = 5.4170520;
        alphaiw[1] = 5.3885820;
        alphaiw[2] = 5.3062630;
        alphaiw[3] = 5.0155720;
        alphaiw[4] = 4.6285580;
        alphaiw[5] = 4.2137360;
        alphaiw[6] = 3.8104020;
        alphaiw[7] = 3.4368990;
        alphaiw[8] = 3.0998560;
        alphaiw[9] = 2.7999440;
        alphaiw[10] = 2.5350070;
        alphaiw[11] = 2.3017150;
        alphaiw[12] = 1.9155660;
        alphaiw[13] = 1.6145960;
        alphaiw[14] = 1.3772910;
        alphaiw[15] = 1.1876550;
        alphaiw[16] = 1.0340620;
        alphaiw[17] = 0.7578460;
        alphaiw[18] = 0.5783980;
        alphaiw[19] = 0.3678490;
        alphaiw[20] = 0.2541230;
        alphaiw[21] = 0.1255490;
        alphaiw[22] = 0.0746120;
      }
      break;
    case 10:
      if (refsys == 1) {
        alphaiw[0] = 2.6227940;
        alphaiw[1] = 2.6160810;
        alphaiw[2] = 2.5963010;
        alphaiw[3] = 2.5221150;
        alphaiw[4] = 2.4130070;
        alphaiw[5] = 2.2831850;
        alphaiw[6] = 2.1442240;
        alphaiw[7] = 2.0040890;
        alphaiw[8] = 1.8677460;
        alphaiw[9] = 1.7380750;
        alphaiw[10] = 1.6165850;
        alphaiw[11] = 1.5039070;
        alphaiw[12] = 1.3049180;
        alphaiw[13] = 1.1383170;
        alphaiw[14] = 0.9993630;
        alphaiw[15] = 0.8832020;
        alphaiw[16] = 0.7855560;
        alphaiw[17] = 0.6010260;
        alphaiw[18] = 0.4739350;
        alphaiw[19] = 0.3155110;
        alphaiw[20] = 0.2245720;
        alphaiw[21] = 0.1157090;
        alphaiw[22] = 0.0702160;
      }
      break;
    case 11:
      if (refsys == 1) {
        alphaiw[0] = 170.9285310;
        alphaiw[1] = 119.1946230;
        alphaiw[2] = 62.8166710;
        alphaiw[3] = 22.2250350;
        alphaiw[4] = 11.0609650;
        alphaiw[5] = 6.7225680;
        alphaiw[6] = 4.6167630;
        alphaiw[7] = 3.4360280;
        alphaiw[8] = 2.7038020;
        alphaiw[9] = 2.2145400;
        alphaiw[10] = 1.8682650;
        alphaiw[11] = 1.6117710;
        alphaiw[12] = 1.2586520;
        alphaiw[13] = 1.0272610;
        alphaiw[14] = 0.8636930;
        alphaiw[15] = 0.7418060;
        alphaiw[16] = 0.6474800;
        alphaiw[17] = 0.4849210;
        alphaiw[18] = 0.3823130;
        alphaiw[19] = 0.2620130;
        alphaiw[20] = 0.1947130;
        alphaiw[21] = 0.1108050;
        alphaiw[22] = 0.0717200;
      }
      else if (refsys == 2) {
        alphaiw[0] = 40.8166710;
        alphaiw[1] = 37.3885700;
        alphaiw[2] = 30.2959980;
        alphaiw[3] = 18.3845350;
        alphaiw[4] = 11.8437670;
        alphaiw[5] = 8.2528940;
        alphaiw[6] = 6.1199710;
        alphaiw[7] = 4.7572100;
        alphaiw[8] = 3.8333160;
        alphaiw[9] = 3.1762600;
        alphaiw[10] = 2.6903860;
        alphaiw[11] = 2.3193520;
        alphaiw[12] = 1.7948680;
        alphaiw[13] = 1.4454310;
        alphaiw[14] = 1.1979870;
        alphaiw[15] = 1.0146980;
        alphaiw[16] = 0.8742120;
        alphaiw[17] = 0.6365370;
        alphaiw[18] = 0.4905090;
        alphaiw[19] = 0.3247020;
        alphaiw[20] = 0.2353690;
        alphaiw[21] = 0.1290210;
        alphaiw[22] = 0.0819540;
      }
      else if (refsys == 3) {
        alphaiw[0] = 92.6411774;
        alphaiw[1] = 91.4243783;
        alphaiw[2] = 88.0445743;
        alphaiw[3] = 77.4364073;
        alphaiw[4] = 65.6163384;
        alphaiw[5] = 54.9685877;
        alphaiw[6] = 46.0979628;
        alphaiw[7] = 38.9126798;
        alphaiw[8] = 33.1338025;
        alphaiw[9] = 28.4756146;
        alphaiw[10] = 24.6963717;
        alphaiw[11] = 21.6046368;
        alphaiw[12] = 16.9269705;
        alphaiw[13] = 13.6278912;
        alphaiw[14] = 11.2244983;
        alphaiw[15] = 9.4230737;
        alphaiw[16] = 8.0391231;
        alphaiw[17] = 5.7141662;
        alphaiw[18] = 4.3121698;
        alphaiw[19] = 2.9652901;
        alphaiw[20] = 1.9616546;
        alphaiw[21] = 1.0496861;
        alphaiw[22] = 0.6661557;
      }
      break;
    case 12:
      if (refsys == 1) {
        alphaiw[0] = 75.6666810;
        alphaiw[1] = 68.6955210;
        alphaiw[2] = 53.9153050;
        alphaiw[3] = 29.2692840;
        alphaiw[4] = 16.8382330;
        alphaiw[5] = 10.6997010;
        alphaiw[6] = 7.3793160;
        alphaiw[7] = 5.4168470;
        alphaiw[8] = 4.1702110;
        alphaiw[9] = 3.3314920;
        alphaiw[10] = 2.7405460;
        alphaiw[11] = 2.3081520;
        alphaiw[12] = 1.7285400;
        alphaiw[13] = 1.3653680;
        alphaiw[14] = 1.1201210;
        alphaiw[15] = 0.9447260;
        alphaiw[16] = 0.8135180;
        alphaiw[17] = 0.5959160;
        alphaiw[18] = 0.4626290;
        alphaiw[19] = 0.3079420;
        alphaiw[20] = 0.2219990;
        alphaiw[21] = 0.1188580;
        alphaiw[22] = 0.0745220;
      }
      else if (refsys == 2) {
        alphaiw[0] = 65.0041490;
        alphaiw[1] = 57.3670330;
        alphaiw[2] = 44.6926220;
        alphaiw[3] = 26.8427890;
        alphaiw[4] = 17.2249920;
        alphaiw[5] = 11.8363720;
        alphaiw[6] = 8.6109420;
        alphaiw[7] = 6.5554130;
        alphaiw[8] = 5.1738650;
        alphaiw[9] = 4.2034870;
        alphaiw[10] = 3.4965440;
        alphaiw[11] = 2.9655110;
        alphaiw[12] = 2.2334680;
        alphaiw[13] = 1.7623610;
        alphaiw[14] = 1.4391360;
        alphaiw[15] = 1.2060360;
        alphaiw[16] = 1.0311360;
        alphaiw[17] = 0.7419490;
        alphaiw[18] = 0.5670970;
        alphaiw[19] = 0.3688080;
        alphaiw[20] = 0.2617440;
        alphaiw[21] = 0.1369550;
        alphaiw[22] = 0.0847940;
      }
      else if (refsys == 3) {
        alphaiw[0] = 34.5442750;
        alphaiw[1] = 33.5344340;
        alphaiw[2] = 30.8656630;
        alphaiw[3] = 23.6561780;
        alphaiw[4] = 17.3397030;
        alphaiw[5] = 12.8432530;
        alphaiw[6] = 9.7723270;
        alphaiw[7] = 7.6521380;
        alphaiw[8] = 6.1504810;
        alphaiw[9] = 5.0567630;
        alphaiw[10] = 4.2388110;
        alphaiw[11] = 3.6123090;
        alphaiw[12] = 2.7314370;
        alphaiw[13] = 2.1542960;
        alphaiw[14] = 1.7542370;
        alphaiw[15] = 1.4641950;
        alphaiw[16] = 1.2461470;
        alphaiw[17] = 0.8862330;
        alphaiw[18] = 0.6703080;
        alphaiw[19] = 0.4289380;
        alphaiw[20] = 0.3010120;
        alphaiw[21] = 0.1548430;
        alphaiw[22] = 0.0949500;
      }
      else if (refsys == 4) {
        alphaiw[0] = 74.6298821;
        alphaiw[1] = 72.6677188;
        alphaiw[2] = 67.9314704;
        alphaiw[3] = 56.5001225;
        alphaiw[4] = 46.4504007;
        alphaiw[5] = 38.4956042;
        alphaiw[6] = 32.3037179;
        alphaiw[7] = 27.4626661;
        alphaiw[8] = 23.6358080;
        alphaiw[9] = 20.5710128;
        alphaiw[10] = 18.0838250;
        alphaiw[11] = 16.0395958;
        alphaiw[12] = 12.9096858;
        alphaiw[13] = 10.6544671;
        alphaiw[14] = 8.9710814;
        alphaiw[15] = 7.6771288;
        alphaiw[16] = 6.6579370;
        alphaiw[17] = 4.8778230;
        alphaiw[18] = 3.8586312;
        alphaiw[19] = 2.8394394;
        alphaiw[20] = 1.7107782;
        alphaiw[21] = 0.8784693;
        alphaiw[22] = 0.5365232;
      }
      break;
    case 13:
      if (refsys == 1) {
        alphaiw[0] = 60.8923430;
        alphaiw[1] = 56.3100140;
        alphaiw[2] = 46.7886820;
        alphaiw[3] = 29.5825800;
        alphaiw[4] = 18.9999900;
        alphaiw[5] = 12.8554770;
        alphaiw[6] = 9.1653720;
        alphaiw[7] = 6.8351170;
        alphaiw[8] = 5.2901970;
        alphaiw[9] = 4.2211880;
        alphaiw[10] = 3.4539770;
        alphaiw[11] = 2.8860210;
        alphaiw[12] = 2.1186190;
        alphaiw[13] = 1.6376640;
        alphaiw[14] = 1.3157820;
        alphaiw[15] = 1.0890010;
        alphaiw[16] = 0.9224950;
        alphaiw[17] = 0.6556710;
        alphaiw[18] = 0.5002950;
        alphaiw[19] = 0.3292210;
        alphaiw[20] = 0.2377920;
        alphaiw[21] = 0.1286420;
        alphaiw[22] = 0.0807650;
      }
      else if (refsys == 2) {
        alphaiw[0] = 56.2109890;
        alphaiw[1] = 52.4017100;
        alphaiw[2] = 44.3132520;
        alphaiw[3] = 29.4084250;
        alphaiw[4] = 19.9025080;
        alphaiw[5] = 14.0717810;
        alphaiw[6] = 10.3773750;
        alphaiw[7] = 7.9372390;
        alphaiw[8] = 6.2599590;
        alphaiw[9] = 5.0653090;
        alphaiw[10] = 4.1877030;
        alphaiw[11] = 3.5255730;
        alphaiw[12] = 2.6117110;
        alphaiw[13] = 2.0262590;
        alphaiw[14] = 1.6284670;
        alphaiw[15] = 1.3453070;
        alphaiw[16] = 1.1360060;
        alphaiw[17] = 0.7989100;
        alphaiw[18] = 0.6026620;
        alphaiw[19] = 0.3886840;
        alphaiw[20] = 0.2764850;
        alphaiw[21] = 0.1461420;
        alphaiw[22] = 0.0906750;
      }
      else if (refsys == 3) {
        alphaiw[0] = 43.0078740;
        alphaiw[1] = 40.3049500;
        alphaiw[2] = 35.7272810;
        alphaiw[3] = 26.8099970;
        alphaiw[4] = 19.8081180;
        alphaiw[5] = 14.8306090;
        alphaiw[6] = 11.3606550;
        alphaiw[7] = 8.9171350;
        alphaiw[8] = 7.1607780;
        alphaiw[9] = 5.8687510;
        alphaiw[10] = 4.8965000;
        alphaiw[11] = 4.1493490;
        alphaiw[12] = 3.0979680;
        alphaiw[13] = 2.4116390;
        alphaiw[14] = 1.9395410;
        alphaiw[15] = 1.6007580;
        alphaiw[16] = 1.3490360;
        alphaiw[17] = 0.9419730;
        alphaiw[18] = 0.7048990;
        alphaiw[19] = 0.4480100;
        alphaiw[20] = 0.3150480;
        alphaiw[21] = 0.1635630;
        alphaiw[22] = 0.1005420;
      }
      else if (refsys == 4) {
        alphaiw[0] = 32.2831180;
        alphaiw[1] = 31.7115930;
        alphaiw[2] = 30.1255060;
        alphaiw[3] = 25.2170820;
        alphaiw[4] = 20.0120360;
        alphaiw[5] = 15.6749350;
        alphaiw[6] = 12.3626340;
        alphaiw[7] = 9.8945680;
        alphaiw[8] = 8.0530620;
        alphaiw[9] = 6.6626970;
        alphaiw[10] = 5.5965060;
        alphaiw[11] = 4.7654850;
        alphaiw[12] = 3.5787890;
        alphaiw[13] = 2.7932030;
        alphaiw[14] = 2.2479080;
        alphaiw[15] = 1.8542580;
        alphaiw[16] = 1.5606340;
        alphaiw[17] = 1.0843300;
        alphaiw[18] = 0.8067590;
        alphaiw[19] = 0.5072070;
        alphaiw[20] = 0.3535600;
        alphaiw[21] = 0.1809750;
        alphaiw[22] = 0.1104040;
      }
      break;
    case 14:
      if (refsys == 1) {
        alphaiw[0] = 38.9490650;
        alphaiw[1] = 37.5127500;
        alphaiw[2] = 33.9244170;
        alphaiw[3] = 25.1761170;
        alphaiw[4] = 18.0845250;
        alphaiw[5] = 13.1765080;
        alphaiw[6] = 9.8616700;
        alphaiw[7] = 7.5917060;
        alphaiw[8] = 5.9982420;
        alphaiw[9] = 4.8493400;
        alphaiw[10] = 3.9995310;
        alphaiw[11] = 3.3561160;
        alphaiw[12] = 2.4668840;
        alphaiw[13] = 1.8982190;
        alphaiw[14] = 1.5134620;
        alphaiw[15] = 1.2411020;
        alphaiw[16] = 1.0410670;
        alphaiw[17] = 0.7226000;
        alphaiw[18] = 0.5405560;
        alphaiw[19] = 0.3465530;
        alphaiw[20] = 0.2471760;
        alphaiw[21] = 0.1332610;
        alphaiw[22] = 0.0845050;
      }
      else if (refsys == 2) {
        alphaiw[0] = 38.0097750;
        alphaiw[1] = 36.6441270;
        alphaiw[2] = 33.3701090;
        alphaiw[3] = 25.5897470;
        alphaiw[4] = 19.0941400;
        alphaiw[5] = 14.3691430;
        alphaiw[6] = 11.0297080;
        alphaiw[7] = 8.6569570;
        alphaiw[8] = 6.9417830;
        alphaiw[9] = 5.6758970;
        alphaiw[10] = 4.7217950;
        alphaiw[11] = 3.9882980;
        alphaiw[12] = 2.9570580;
        alphaiw[13] = 2.2858030;
        alphaiw[14] = 1.8259530;
        alphaiw[15] = 1.4975500;
        alphaiw[16] = 1.2548440;
        alphaiw[17] = 0.8660940;
        alphaiw[18] = 0.6430810;
        alphaiw[19] = 0.4060490;
        alphaiw[20] = 0.2858610;
        alphaiw[21] = 0.1507420;
        alphaiw[22] = 0.0943940;
      }
      else if (refsys == 3) {
        alphaiw[0] = 37.2763200;
        alphaiw[1] = 36.2054630;
        alphaiw[2] = 33.4711240;
        alphaiw[3] = 26.4376500;
        alphaiw[4] = 20.2438260;
        alphaiw[5] = 15.5716430;
        alphaiw[6] = 12.1648030;
        alphaiw[7] = 9.6798740;
        alphaiw[8] = 7.8448670;
        alphaiw[9] = 6.4670510;
        alphaiw[10] = 5.4140820;
        alphaiw[11] = 4.5954000;
        alphaiw[12] = 3.4297600;
        alphaiw[13] = 2.6610590;
        alphaiw[14] = 2.1295730;
        alphaiw[15] = 1.7474840;
        alphaiw[16] = 1.4637430;
        alphaiw[17] = 1.0070410;
        alphaiw[18] = 0.7441640;
        alphaiw[19] = 0.4650030;
        alphaiw[20] = 0.3243180;
        alphaiw[21] = 0.1681920;
        alphaiw[22] = 0.1042810;
      }
      else if (refsys == 4) {
        alphaiw[0] = 75.2662880;
        alphaiw[1] = 72.3612970;
        alphaiw[2] = 65.6345910;
        alphaiw[3] = 50.8100850;
        alphaiw[4] = 38.9033790;
        alphaiw[5] = 30.0787730;
        alphaiw[6] = 23.6268450;
        alphaiw[7] = 18.8896320;
        alphaiw[8] = 15.3681050;
        alphaiw[9] = 12.7085630;
        alphaiw[10] = 10.6661350;
        alphaiw[11] = 9.0717390;
        alphaiw[12] = 6.7912830;
        alphaiw[13] = 5.2803260;
        alphaiw[14] = 4.2321590;
        alphaiw[15] = 3.4767820;
        alphaiw[16] = 2.9147970;
        alphaiw[17] = 2.0082810;
        alphaiw[18] = 1.4853250;
        alphaiw[19] = 0.9289700;
        alphaiw[20] = 0.6481770;
        alphaiw[21] = 0.3362410;
        alphaiw[22] = 0.2084750;
      }
      else if (refsys == 5) {
        alphaiw[0] = 61.7861480;
        alphaiw[1] = 60.8116450;
        alphaiw[2] = 58.0975290;
        alphaiw[3] = 49.5798710;
        alphaiw[4] = 40.2758690;
        alphaiw[5] = 32.2331470;
        alphaiw[6] = 25.8672360;
        alphaiw[7] = 20.9751250;
        alphaiw[8] = 17.2319240;
        alphaiw[9] = 14.3488200;
        alphaiw[10] = 12.1032080;
        alphaiw[11] = 10.3316220;
        alphaiw[12] = 7.7698820;
        alphaiw[13] = 6.0547620;
        alphaiw[14] = 4.8568220;
        alphaiw[15] = 3.9895410;
        alphaiw[16] = 3.3422950;
        alphaiw[17] = 2.2952740;
        alphaiw[18] = 1.6903860;
        alphaiw[19] = 1.0479590;
        alphaiw[20] = 0.7255320;
        alphaiw[21] = 0.3711770;
        alphaiw[22] = 0.2282300;
      }
      break;
    case 15:
      if (refsys == 1) {
        alphaiw[0] = 25.9470070;
        alphaiw[1] = 25.4290340;
        alphaiw[2] = 24.0216260;
        alphaiw[3] = 19.8723550;
        alphaiw[4] = 15.6749860;
        alphaiw[5] = 12.2514680;
        alphaiw[6] = 9.6516650;
        alphaiw[7] = 7.7129690;
        alphaiw[8] = 6.2628530;
        alphaiw[9] = 5.1656700;
        alphaiw[10] = 4.3233980;
        alphaiw[11] = 3.6669330;
        alphaiw[12] = 2.7311570;
        alphaiw[13] = 2.1144260;
        alphaiw[14] = 1.6888290;
        alphaiw[15] = 1.3835900;
        alphaiw[16] = 1.1574800;
        alphaiw[17] = 0.7949870;
        alphaiw[18] = 0.5875340;
        alphaiw[19] = 0.3687610;
        alphaiw[20] = 0.2592740;
        alphaiw[21] = 0.1378720;
        alphaiw[22] = 0.0875080;
      }
      else if (refsys == 2) {
        alphaiw[0] = 27.8753970;
        alphaiw[1] = 27.3122270;
        alphaiw[2] = 25.8232190;
        alphaiw[3] = 21.5634300;
        alphaiw[4] = 17.2591270;
        alphaiw[5] = 13.6872210;
        alphaiw[6] = 10.9190600;
        alphaiw[7] = 8.8162050;
        alphaiw[8] = 7.2183460;
        alphaiw[9] = 5.9933920;
        alphaiw[10] = 5.0426620;
        alphaiw[11] = 4.2948080;
        alphaiw[12] = 3.2172770;
        alphaiw[13] = 2.4989620;
        alphaiw[14] = 1.9991280;
        alphaiw[15] = 1.6384630;
        alphaiw[16] = 1.3701110;
        alphaiw[17] = 0.9379270;
        alphaiw[18] = 0.6897680;
        alphaiw[19] = 0.4281650;
        alphaiw[20] = 0.2979330;
        alphaiw[21] = 0.1553700;
        alphaiw[22] = 0.0974170;
      }
      else if (refsys == 3) {
        alphaiw[0] = 55.7853670;
        alphaiw[1] = 54.5314410;
        alphaiw[2] = 51.2683850;
        alphaiw[3] = 42.3425800;
        alphaiw[4] = 33.7346300;
        alphaiw[5] = 26.7484370;
        alphaiw[6] = 21.3755400;
        alphaiw[7] = 17.2987970;
        alphaiw[8] = 14.1966540;
        alphaiw[9] = 11.8128320;
        alphaiw[10] = 9.9578290;
        alphaiw[11] = 8.4949500;
        alphaiw[12] = 6.3802250;
        alphaiw[13] = 4.9651180;
        alphaiw[14] = 3.9775280;
        alphaiw[15] = 3.2633040;
        alphaiw[16] = 2.7309600;
        alphaiw[17] = 1.8718860;
        alphaiw[18] = 1.3776140;
        alphaiw[19] = 0.8557670;
        alphaiw[20] = 0.5956760;
        alphaiw[21] = 0.3107310;
        alphaiw[22] = 0.1948320;
      }
      else if (refsys == 4) {
        alphaiw[0] = 30.4586930;
        alphaiw[1] = 29.9832590;
        alphaiw[2] = 28.6718190;
        alphaiw[3] = 24.6351360;
        alphaiw[4] = 20.2683100;
        alphaiw[5] = 16.4561420;
        alphaiw[6] = 13.3810200;
        alphaiw[7] = 10.9693570;
        alphaiw[8] = 9.0896850;
        alphaiw[9] = 7.6189190;
        alphaiw[10] = 6.4582620;
        alphaiw[11] = 5.5327160;
        alphaiw[12] = 4.1782340;
        alphaiw[13] = 3.2605350;
        alphaiw[14] = 2.6145330;
        alphaiw[15] = 2.1444890;
        alphaiw[16] = 1.7926280;
        alphaiw[17] = 1.2223890;
        alphaiw[18] = 0.8934200;
        alphaiw[19] = 0.5466340;
        alphaiw[20] = 0.3750780;
        alphaiw[21] = 0.1903100;
        alphaiw[22] = 0.1172040;
      }
      break;
    case 16:
      if (refsys == 1) {
        alphaiw[0] = 19.6135110;
        alphaiw[1] = 19.3396950;
        alphaiw[2] = 18.5777950;
        alphaiw[3] = 16.1740460;
        alphaiw[4] = 13.4828010;
        alphaiw[5] = 11.0662650;
        alphaiw[6] = 9.0764420;
        alphaiw[7] = 7.4918500;
        alphaiw[8] = 6.2419920;
        alphaiw[9] = 5.2546810;
        alphaiw[10] = 4.4695850;
        alphaiw[11] = 3.8397190;
        alphaiw[12] = 2.9119750;
        alphaiw[13] = 2.2796060;
        alphaiw[14] = 1.8328050;
        alphaiw[15] = 1.5069310;
        alphaiw[16] = 1.2625940;
        alphaiw[17] = 0.8658660;
        alphaiw[18] = 0.6364800;
        alphaiw[19] = 0.3940160;
        alphaiw[20] = 0.2736140;
        alphaiw[21] = 0.1429400;
        alphaiw[22] = 0.0902880;
      }
      else if (refsys == 2) {
        alphaiw[0] = 22.2041430;
        alphaiw[1] = 21.8891090;
        alphaiw[2] = 21.0218000;
        alphaiw[3] = 18.3279640;
        alphaiw[4] = 15.3337320;
        alphaiw[5] = 12.6407210;
        alphaiw[6] = 10.4119420;
        alphaiw[7] = 8.6265960;
        alphaiw[8] = 7.2105080;
        alphaiw[9] = 6.0862970;
        alphaiw[10] = 5.1884110;
        alphaiw[11] = 4.4652670;
        alphaiw[12] = 3.3950220;
        alphaiw[13] = 2.6615450;
        alphaiw[14] = 2.1411100;
        alphaiw[15] = 1.7603030;
        alphaiw[16] = 1.4740840;
        alphaiw[17] = 1.0081850;
        alphaiw[18] = 0.7383390;
        alphaiw[19] = 0.4532570;
        alphaiw[20] = 0.3121980;
        alphaiw[21] = 0.1604390;
        alphaiw[22] = 0.1002160;
      }
      else if (refsys == 3) {
        alphaiw[0] = 24.5118610;
        alphaiw[1] = 24.1774360;
        alphaiw[2] = 23.2503900;
        alphaiw[3] = 20.3440110;
        alphaiw[4] = 17.0929990;
        alphaiw[5] = 14.1530500;
        alphaiw[6] = 11.7050130;
        alphaiw[7] = 9.7320240;
        alphaiw[8] = 8.1583820;
        alphaiw[9] = 6.9030290;
        alphaiw[10] = 5.8962250;
        alphaiw[11] = 5.0824590;
        alphaiw[12] = 3.8729230;
        alphaiw[13] = 3.0400830;
        alphaiw[14] = 2.4470470;
        alphaiw[15] = 2.0119580;
        alphaiw[16] = 1.6842890;
        alphaiw[17] = 1.1498080;
        alphaiw[18] = 0.8397730;
        alphaiw[19] = 0.5122980;
        alphaiw[20] = 0.3506710;
        alphaiw[21] = 0.1779010;
        alphaiw[22] = 0.1101270;
      }
      break;
    case 17:
      if (refsys == 1) {
        alphaiw[0] = 14.6510610;
        alphaiw[1] = 14.5166850;
        alphaiw[2] = 14.1329320;
        alphaiw[3] = 12.8305660;
        alphaiw[4] = 11.2118050;
        alphaiw[5] = 9.6110460;
        alphaiw[6] = 8.1817260;
        alphaiw[7] = 6.9658720;
        alphaiw[8] = 5.9541940;
        alphaiw[9] = 5.1193500;
        alphaiw[10] = 4.4310730;
        alphaiw[11] = 3.8619390;
        alphaiw[12] = 2.9935150;
        alphaiw[13] = 2.3788950;
        alphaiw[14] = 1.9324470;
        alphaiw[15] = 1.6000440;
        alphaiw[16] = 1.3469040;
        alphaiw[17] = 0.9286160;
        alphaiw[18] = 0.6827690;
        alphaiw[19] = 0.4202900;
        alphaiw[20] = 0.2895830;
        alphaiw[21] = 0.1489740;
        alphaiw[22] = 0.0934510;
      }
      else if (refsys == 2) {
        alphaiw[0] = 17.2630540;
        alphaiw[1] = 17.0983080;
        alphaiw[2] = 16.6294740;
        alphaiw[3] = 15.0531340;
        alphaiw[4] = 13.1168290;
        alphaiw[5] = 11.2193720;
        alphaiw[6] = 9.5356020;
        alphaiw[7] = 8.1091300;
        alphaiw[8] = 6.9254200;
        alphaiw[9] = 5.9504380;
        alphaiw[10] = 5.1476860;
        alphaiw[11] = 4.4845120;
        alphaiw[12] = 3.4734290;
        alphaiw[13] = 2.7582120;
        alphaiw[14] = 2.2387270;
        alphaiw[15] = 1.8518950;
        alphaiw[16] = 1.5572580;
        alphaiw[17] = 1.0703600;
        alphaiw[18] = 0.7842990;
        alphaiw[19] = 0.4793850;
        alphaiw[20] = 0.3280880;
        alphaiw[21] = 0.1664540;
        alphaiw[22] = 0.1033810;
      }
      break;
    case 18:
      if (refsys == 1) {
        alphaiw[0] = 11.1080370;
        alphaiw[1] = 11.0394140;
        alphaiw[2] = 10.8402330;
        alphaiw[3] = 10.1295640;
        alphaiw[4] = 9.1733700;
        alphaiw[5] = 8.1487470;
        alphaiw[6] = 7.1654820;
        alphaiw[7] = 6.2757630;
        alphaiw[8] = 5.4960060;
        alphaiw[9] = 4.8240080;
        alphaiw[10] = 4.2494460;
        alphaiw[11] = 3.7594860;
        alphaiw[12] = 2.9839560;
        alphaiw[13] = 2.4128240;
        alphaiw[14] = 1.9852170;
        alphaiw[15] = 1.6593100;
        alphaiw[16] = 1.4065440;
        alphaiw[17] = 0.9797740;
        alphaiw[18] = 0.7235180;
        alphaiw[19] = 0.4456840;
        alphaiw[20] = 0.3060220;
        alphaiw[21] = 0.1557480;
        alphaiw[22] = 0.0970350;
      }
      break;
    case 19:
      if (refsys == 1) {
        alphaiw[0] = 366.6173270;
        alphaiw[1] = 195.2407260;
        alphaiw[2] = 83.8654450;
        alphaiw[3] = 28.8051330;
        alphaiw[4] = 15.9216120;
        alphaiw[5] = 10.9496410;
        alphaiw[6] = 8.4035800;
        alphaiw[7] = 6.8431600;
        alphaiw[8] = 5.7666680;
        alphaiw[9] = 4.9640380;
        alphaiw[10] = 4.3344440;
        alphaiw[11] = 3.8237620;
        alphaiw[12] = 3.0427910;
        alphaiw[13] = 2.4756550;
        alphaiw[14] = 2.0494570;
        alphaiw[15] = 1.7214130;
        alphaiw[16] = 1.4641050;
        alphaiw[17] = 1.0222030;
        alphaiw[18] = 0.7516860;
        alphaiw[19] = 0.4544310;
        alphaiw[20] = 0.3045490;
        alphaiw[21] = 0.1458010;
        alphaiw[22] = 0.0862690;
      }
      else if (refsys == 2) {
        alphaiw[0] = 53.4389170;
        alphaiw[1] = 47.8463210;
        alphaiw[2] = 37.7399140;
        alphaiw[3] = 23.4699070;
        alphaiw[4] = 16.3225240;
        alphaiw[5] = 12.3576000;
        alphaiw[6] = 9.8698810;
        alphaiw[7] = 8.1579680;
        alphaiw[8] = 6.9001080;
        alphaiw[9] = 5.9321790;
        alphaiw[10] = 5.1625350;
        alphaiw[11] = 4.5358880;
        alphaiw[12] = 3.5804740;
        alphaiw[13] = 2.8930340;
        alphaiw[14] = 2.3815560;
        alphaiw[15] = 1.9913210;
        alphaiw[16] = 1.6874660;
        alphaiw[17] = 1.1703840;
        alphaiw[18] = 0.8568800;
        alphaiw[19] = 0.5151140;
        alphaiw[20] = 0.3439260;
        alphaiw[21] = 0.1635940;
        alphaiw[22] = 0.0963540;
      }
      else if (refsys == 3) {
        alphaiw[0] = 112.8769963;
        alphaiw[1] = 111.5382559;
        alphaiw[2] = 107.8201226;
        alphaiw[3] = 96.1458563;
        alphaiw[4] = 83.0747922;
        alphaiw[5] = 71.1307898;
        alphaiw[6] = 60.9421945;
        alphaiw[7] = 52.4432551;
        alphaiw[8] = 45.3906629;
        alphaiw[9] = 39.5298970;
        alphaiw[10] = 34.6393638;
        alphaiw[11] = 30.5368708;
        alphaiw[12] = 24.1382573;
        alphaiw[13] = 19.4781015;
        alphaiw[14] = 16.0067586;
        alphaiw[15] = 13.3663894;
        alphaiw[16] = 11.3195232;
        alphaiw[17] = 7.8610363;
        alphaiw[18] = 5.7811058;
        alphaiw[19] = 3.5249696;
        alphaiw[20] = 2.3940353;
        alphaiw[21] = 1.1879380;
        alphaiw[22] = 0.7245068;
      }
      break;
    case 20:
      if (refsys == 1) {
        alphaiw[0] = 167.3331200;
        alphaiw[1] = 138.7350820;
        alphaiw[2] = 92.3660100;
        alphaiw[3] = 41.0298690;
        alphaiw[4] = 22.5025050;
        alphaiw[5] = 14.5451620;
        alphaiw[6] = 10.4632190;
        alphaiw[7] = 8.0723760;
        alphaiw[8] = 6.5267080;
        alphaiw[9] = 5.4500180;
        alphaiw[10] = 4.6561940;
        alphaiw[11] = 4.0450480;
        alphaiw[12] = 3.1618950;
        alphaiw[13] = 2.5524970;
        alphaiw[14] = 2.1073550;
        alphaiw[15] = 1.7698140;
        alphaiw[16] = 1.5069710;
        alphaiw[17] = 1.0569390;
        alphaiw[18] = 0.7807210;
        alphaiw[19] = 0.4750180;
        alphaiw[20] = 0.3195680;
        alphaiw[21] = 0.1538130;
        alphaiw[22] = 0.0914010;
      }
      else if (refsys == 2) {
        alphaiw[0] = 158.1508830;
        alphaiw[1] = 111.6461600;
        alphaiw[2] = 67.9213530;
        alphaiw[3] = 33.9702430;
        alphaiw[4] = 21.3038760;
        alphaiw[5] = 15.0163380;
        alphaiw[6] = 11.3779670;
        alphaiw[7] = 9.0482850;
        alphaiw[8] = 7.4431400;
        alphaiw[9] = 6.2742820;
        alphaiw[10] = 5.3861210;
        alphaiw[11] = 4.6885960;
        alphaiw[12] = 3.6641970;
        alphaiw[13] = 2.9506000;
        alphaiw[14] = 2.4285760;
        alphaiw[15] = 2.0334750;
        alphaiw[16] = 1.7267590;
        alphaiw[17] = 1.2044810;
        alphaiw[18] = 0.8862100;
        alphaiw[19] = 0.5363960;
        alphaiw[20] = 0.3596090;
        alphaiw[21] = 0.1720510;
        alphaiw[22] = 0.1017830;
      }
      else if (refsys == 3) {
        alphaiw[0] = 47.8971840;
        alphaiw[1] = 45.8022250;
        alphaiw[2] = 40.7933190;
        alphaiw[3] = 29.6764790;
        alphaiw[4] = 21.5460190;
        alphaiw[5] = 16.2269040;
        alphaiw[6] = 12.6953380;
        alphaiw[7] = 10.2550730;
        alphaiw[8] = 8.4973390;
        alphaiw[9] = 7.1832110;
        alphaiw[10] = 6.1691650;
        alphaiw[11] = 5.3658830;
        alphaiw[12] = 4.1801140;
        alphaiw[13] = 3.3536760;
        alphaiw[14] = 2.7508620;
        alphaiw[15] = 2.2964010;
        alphaiw[16] = 1.9450020;
        alphaiw[17] = 1.3500690;
        alphaiw[18] = 0.9899710;
        alphaiw[19] = 0.5965960;
        alphaiw[20] = 0.3988370;
        alphaiw[21] = 0.1899030;
        alphaiw[22] = 0.1119470;
      }
      else if (refsys == 4) {
        alphaiw[0] = 92.8563023;
        alphaiw[1] = 91.1707612;
        alphaiw[2] = 86.7243129;
        alphaiw[3] = 74.2666422;
        alphaiw[4] = 62.0928688;
        alphaiw[5] = 52.1273713;
        alphaiw[6] = 44.3017181;
        alphaiw[7] = 38.1421847;
        alphaiw[8] = 33.2177260;
        alphaiw[9] = 29.2098267;
        alphaiw[10] = 25.8942119;
        alphaiw[11] = 23.1129142;
        alphaiw[12] = 18.7302605;
        alphaiw[13] = 15.4634027;
        alphaiw[14] = 12.9627169;
        alphaiw[15] = 11.0083621;
        alphaiw[16] = 9.4545882;
        alphaiw[17] = 6.7335695;
        alphaiw[18] = 5.0260366;
        alphaiw[19] = 3.0995061;
        alphaiw[20] = 2.1015110;
        alphaiw[21] = 1.0184280;
        alphaiw[22] = 0.6045719;
      }
      break;
    case 21:
      if (refsys == 1) {
        alphaiw[0] = 125.6799050;
        alphaiw[1] = 108.3144390;
        alphaiw[2] = 77.2929700;
        alphaiw[3] = 37.1320260;
        alphaiw[4] = 20.7670580;
        alphaiw[5] = 13.4119080;
        alphaiw[6] = 9.5811350;
        alphaiw[7] = 7.3353090;
        alphaiw[8] = 5.8931080;
        alphaiw[9] = 4.8989950;
        alphaiw[10] = 4.1745440;
        alphaiw[11] = 3.6229180;
        alphaiw[12] = 2.8357880;
        alphaiw[13] = 2.2984530;
        alphaiw[14] = 1.9074130;
        alphaiw[15] = 1.6106400;
        alphaiw[16] = 1.3787470;
        alphaiw[17] = 0.9784770;
        alphaiw[18] = 0.7296280;
        alphaiw[19] = 0.4499330;
        alphaiw[20] = 0.3052790;
        alphaiw[21] = 0.1484840;
        alphaiw[22] = 0.0886520;
      }
      else if (refsys == 2) {
        alphaiw[0] = 43.3830910;
        alphaiw[1] = 45.0959100;
        alphaiw[2] = 39.5276970;
        alphaiw[3] = 28.9987430;
        alphaiw[4] = 21.3879050;
        alphaiw[5] = 16.2717510;
        alphaiw[6] = 12.7943430;
        alphaiw[7] = 10.3540150;
        alphaiw[8] = 8.5813500;
        alphaiw[9] = 7.2517530;
        alphaiw[10] = 6.2258690;
        alphaiw[11] = 5.4147740;
        alphaiw[12] = 4.2221260;
        alphaiw[13] = 3.3947970;
        alphaiw[14] = 2.7929060;
        alphaiw[15] = 2.3393880;
        alphaiw[16] = 1.9883760;
        alphaiw[17] = 1.3919190;
        alphaiw[18] = 1.0282690;
        alphaiw[19] = 0.6267560;
        alphaiw[20] = 0.4221330;
        alphaiw[21] = 0.2026710;
        alphaiw[22] = 0.1197200;
      }
      else if (refsys == 3) {
        alphaiw[0] = 45.7861990;
        alphaiw[1] = 44.3266490;
        alphaiw[2] = 40.7084650;
        alphaiw[3] = 31.7823730;
        alphaiw[4] = 24.2343690;
        alphaiw[5] = 18.7383620;
        alphaiw[6] = 14.8308240;
        alphaiw[7] = 12.0175160;
        alphaiw[8] = 9.9446720;
        alphaiw[9] = 8.3787230;
        alphaiw[10] = 7.1672090;
        alphaiw[11] = 6.2095210;
        alphaiw[12] = 4.8060710;
        alphaiw[13] = 3.8394390;
        alphaiw[14] = 3.1416970;
        alphaiw[15] = 2.6197890;
        alphaiw[16] = 2.2184260;
        alphaiw[17] = 1.5421450;
        alphaiw[18] = 1.1337000;
        alphaiw[19] = 0.6865730;
        alphaiw[20] = 0.4605190;
        alphaiw[21] = 0.2197990;
        alphaiw[22] = 0.1294090;
      }
      else if (refsys == 4) {
        alphaiw[0] = 25.7885393;
        alphaiw[1] = 25.5284322;
        alphaiw[2] = 24.8018061;
        alphaiw[3] = 22.4847968;
        alphaiw[4] = 19.8432092;
        alphaiw[5] = 17.4046425;
        alphaiw[6] = 15.3140809;
        alphaiw[7] = 13.5610822;
        alphaiw[8] = 12.0932363;
        alphaiw[9] = 10.8561131;
        alphaiw[10] = 9.8039073;
        alphaiw[11] = 8.9004931;
        alphaiw[12] = 7.4345994;
        alphaiw[13] = 6.3020196;
        alphaiw[14] = 5.4061411;
        alphaiw[15] = 4.6844025;
        alphaiw[16] = 4.0943751;
        alphaiw[17] = 3.0185882;
        alphaiw[18] = 2.3095039;
        alphaiw[19] = 1.4697240;
        alphaiw[20] = 1.0145692;
        alphaiw[21] = 0.5013176;
        alphaiw[22] = 0.2990179;
      }
      break;
    case 22:
      if (refsys == 1) {
        alphaiw[0] = 114.7243640;
        alphaiw[1] = 99.7786920;
        alphaiw[2] = 73.3046090;
        alphaiw[3] = 37.5872700;
        alphaiw[4] = 22.0065650;
        alphaiw[5] = 14.6330060;
        alphaiw[6] = 10.6399730;
        alphaiw[7] = 8.2278720;
        alphaiw[8] = 6.6441090;
        alphaiw[9] = 5.5356020;
        alphaiw[10] = 4.7202080;
        alphaiw[11] = 4.0964710;
        alphaiw[12] = 3.2058500;
        alphaiw[13] = 2.6004600;
        alphaiw[14] = 2.1623770;
        alphaiw[15] = 1.8313070;
        alphaiw[16] = 1.5731530;
        alphaiw[17] = 1.1271400;
        alphaiw[18] = 0.8480460;
        alphaiw[19] = 0.5303920;
        alphaiw[20] = 0.3633330;
        alphaiw[21] = 0.1787040;
        alphaiw[22] = 0.1068890;
      }
      else if (refsys == 2) {
        alphaiw[0] = 46.0919720;
        alphaiw[1] = 42.4892580;
        alphaiw[2] = 37.6213730;
        alphaiw[3] = 28.0117340;
        alphaiw[4] = 20.8457830;
        alphaiw[5] = 15.9626110;
        alphaiw[6] = 12.6095600;
        alphaiw[7] = 10.2373540;
        alphaiw[8] = 8.5034920;
        alphaiw[9] = 7.1973670;
        alphaiw[10] = 6.1869270;
        alphaiw[11] = 5.3870300;
        alphaiw[12] = 4.2107080;
        alphaiw[13] = 3.3956060;
        alphaiw[14] = 2.8031890;
        alphaiw[15] = 2.3567820;
        alphaiw[16] = 2.0108410;
        alphaiw[17] = 1.4205320;
        alphaiw[18] = 1.0576160;
        alphaiw[19] = 0.6522060;
        alphaiw[20] = 0.4427770;
        alphaiw[21] = 0.2149140;
        alphaiw[22] = 0.1275310;
      }
      else if (refsys == 3) {
        alphaiw[0] = 48.7114030;
        alphaiw[1] = 47.2968810;
        alphaiw[2] = 43.6774280;
        alphaiw[3] = 34.4392690;
        alphaiw[4] = 26.5310300;
        alphaiw[5] = 20.7154850;
        alphaiw[6] = 16.5179260;
        alphaiw[7] = 13.4481710;
        alphaiw[8] = 11.1565010;
        alphaiw[9] = 9.4084180;
        alphaiw[10] = 8.0473770;
        alphaiw[11] = 6.9676210;
        alphaiw[12] = 5.3829550;
        alphaiw[13] = 4.2933920;
        alphaiw[14] = 3.5097360;
        alphaiw[15] = 2.9257120;
        alphaiw[16] = 2.4779210;
        alphaiw[17] = 1.7256070;
        alphaiw[18] = 1.2717960;
        alphaiw[19] = 0.7740440;
        alphaiw[20] = 0.5212220;
        alphaiw[21] = 0.2500910;
        alphaiw[22] = 0.1474550;
      }
      else if (refsys == 4) {
        alphaiw[0] = 31.8693190;
        alphaiw[1] = 31.5954595;
        alphaiw[2] = 30.8154120;
        alphaiw[3] = 28.1936583;
        alphaiw[4] = 24.9990546;
        alphaiw[5] = 21.9123454;
        alphaiw[6] = 19.2042984;
        alphaiw[7] = 16.9162029;
        alphaiw[8] = 15.0030479;
        alphaiw[9] = 13.4005228;
        alphaiw[10] = 12.0488331;
        alphaiw[11] = 10.8986824;
        alphaiw[12] = 9.0558940;
        alphaiw[13] = 7.6529726;
        alphaiw[14] = 6.5554549;
        alphaiw[15] = 5.6780283;
        alphaiw[16] = 4.9642843;
        alphaiw[17] = 3.6680250;
        alphaiw[18] = 2.8144887;
        alphaiw[19] = 1.8007388;
        alphaiw[20] = 1.2482416;
        alphaiw[21] = 0.6204756;
        alphaiw[22] = 0.3709690;
      }
      break;
    case 23:
      if (refsys == 1) {
        alphaiw[0] = 100.9910610;
        alphaiw[1] = 88.6597800;
        alphaiw[2] = 66.6115730;
        alphaiw[3] = 35.5163660;
        alphaiw[4] = 21.2160720;
        alphaiw[5] = 14.2442310;
        alphaiw[6] = 10.4067670;
        alphaiw[7] = 8.0678270;
        alphaiw[8] = 6.5248700;
        alphaiw[9] = 5.4425580;
        alphaiw[10] = 4.6458490;
        alphaiw[11] = 4.0364450;
        alphaiw[12] = 3.1669610;
        alphaiw[13] = 2.5767360;
        alphaiw[14] = 2.1500770;
        alphaiw[15] = 1.8277610;
        alphaiw[16] = 1.5763140;
        alphaiw[17] = 1.1406790;
        alphaiw[18] = 0.8660900;
        alphaiw[19] = 0.5494430;
        alphaiw[20] = 0.3799770;
        alphaiw[21] = 0.1891620;
        alphaiw[22] = 0.1136620;
      }
      else if (refsys == 2) {
        alphaiw[0] = 44.7007430;
        alphaiw[1] = 42.0873530;
        alphaiw[2] = 37.1472140;
        alphaiw[3] = 27.4044260;
        alphaiw[4] = 20.3770970;
        alphaiw[5] = 15.6299000;
        alphaiw[6] = 12.3710770;
        alphaiw[7] = 10.0613500;
        alphaiw[8] = 8.3699050;
        alphaiw[9] = 7.0937900;
        alphaiw[10] = 6.1054940;
        alphaiw[11] = 5.3225590;
        alphaiw[12] = 4.1705140;
        alphaiw[13] = 3.3719800;
        alphaiw[14] = 2.7915250;
        alphaiw[15] = 2.3540000;
        alphaiw[16] = 2.0147100;
        alphaiw[17] = 1.4344300;
        alphaiw[18] = 1.0757290;
        alphaiw[19] = 0.6710700;
        alphaiw[20] = 0.4591900;
        alphaiw[21] = 0.2251960;
        alphaiw[22] = 0.1341810;
      }
      else if (refsys == 3) {
        alphaiw[0] = 43.8589710;
        alphaiw[1] = 41.9185700;
        alphaiw[2] = 37.9612250;
        alphaiw[3] = 29.4769420;
        alphaiw[4] = 22.6653000;
        alphaiw[5] = 17.7162500;
        alphaiw[6] = 14.1578990;
        alphaiw[7] = 11.5615470;
        alphaiw[8] = 9.6256600;
        alphaiw[9] = 8.1492470;
        alphaiw[10] = 6.9989060;
        alphaiw[11] = 6.0850100;
        alphaiw[12] = 4.7397660;
        alphaiw[13] = 3.8101490;
        alphaiw[14] = 3.1377770;
        alphaiw[15] = 2.6337640;
        alphaiw[16] = 2.2450370;
        alphaiw[17] = 1.5855880;
        alphaiw[18] = 1.1821140;
        alphaiw[19] = 0.7316980;
        alphaiw[20] = 0.4982510;
        alphaiw[21] = 0.2427270;
        alphaiw[22] = 0.1441100;
      }
      else if (refsys == 4) {
        alphaiw[0] = 39.4919722;
        alphaiw[1] = 39.0027923;
        alphaiw[2] = 37.6740133;
        alphaiw[3] = 33.6693555;
        alphaiw[4] = 29.3313719;
        alphaiw[5] = 25.4515734;
        alphaiw[6] = 22.2023666;
        alphaiw[7] = 19.5294838;
        alphaiw[8] = 17.3268203;
        alphaiw[9] = 15.4949224;
        alphaiw[10] = 13.9540832;
        alphaiw[11] = 12.6434247;
        alphaiw[12] = 10.5397592;
        alphaiw[13] = 8.9322131;
        alphaiw[14] = 7.6697405;
        alphaiw[15] = 6.6569547;
        alphaiw[16] = 5.8306559;
        alphaiw[17] = 4.3239152;
        alphaiw[18] = 3.3270116;
        alphaiw[19] = 2.1370028;
        alphaiw[20] = 1.4850445;
        alphaiw[21] = 0.7404650;
        alphaiw[22] = 0.4431346;
      }
      break;
    case 24:
      if (refsys == 1) {
        alphaiw[0] = 81.4204420;
        alphaiw[1] = 70.3861320;
        alphaiw[2] = 51.0334860;
        alphaiw[3] = 26.6035920;
        alphaiw[4] = 16.4401450;
        alphaiw[5] = 11.5586160;
        alphaiw[6] = 8.7956620;
        alphaiw[7] = 7.0404390;
        alphaiw[8] = 5.8344750;
        alphaiw[9] = 4.9584520;
        alphaiw[10] = 4.2948990;
        alphaiw[11] = 3.7754900;
        alphaiw[12] = 3.0150310;
        alphaiw[13] = 2.4847490;
        alphaiw[14] = 2.0937350;
        alphaiw[15] = 1.7938190;
        alphaiw[16] = 1.5570520;
        alphaiw[17] = 1.1409650;
        alphaiw[18] = 0.8747790;
        alphaiw[19] = 0.5634010;
        alphaiw[20] = 0.3939510;
        alphaiw[21] = 0.1991850;
        alphaiw[22] = 0.1204340;
      }
      else if (refsys == 2) {
        alphaiw[0] = 45.3318020;
        alphaiw[1] = 42.4171170;
        alphaiw[2] = 37.0714340;
        alphaiw[3] = 26.9188730;
        alphaiw[4] = 19.8754030;
        alphaiw[5] = 15.2233940;
        alphaiw[6] = 12.0588400;
        alphaiw[7] = 9.8211890;
        alphaiw[8] = 8.1820250;
        alphaiw[9] = 6.9441660;
        alphaiw[10] = 5.9847000;
        alphaiw[11] = 5.2241990;
        alphaiw[12] = 4.1047330;
        alphaiw[13] = 3.3284220;
        alphaiw[14] = 2.7636460;
        alphaiw[15] = 2.3373710;
        alphaiw[16] = 2.0062470;
        alphaiw[17] = 1.4381260;
        alphaiw[18] = 1.0851790;
        alphaiw[19] = 0.6841550;
        alphaiw[20] = 0.4719940;
        alphaiw[21] = 0.2343080;
        alphaiw[22] = 0.1403490;
      }
      else if (refsys == 3) {
        alphaiw[0] = 56.4761320;
        alphaiw[1] = 55.4404440;
        alphaiw[2] = 52.7176090;
        alphaiw[3] = 45.1454930;
        alphaiw[4] = 37.6466880;
        alphaiw[5] = 31.2740070;
        alphaiw[6] = 26.0862050;
        alphaiw[7] = 21.9197000;
        alphaiw[8] = 18.5802620;
        alphaiw[9] = 15.8937900;
        alphaiw[10] = 13.7180050;
        alphaiw[11] = 11.9412500;
        alphaiw[12] = 9.2607000;
        alphaiw[13] = 7.3778240;
        alphaiw[14] = 6.0116860;
        alphaiw[15] = 4.9917560;
        alphaiw[16] = 4.2112130;
        alphaiw[17] = 2.9083100;
        alphaiw[18] = 2.1307140;
        alphaiw[19] = 1.2874000;
        alphaiw[20] = 0.8631570;
        alphaiw[21] = 0.4110260;
        alphaiw[22] = 0.2406390;
      }
      else if (refsys == 4) {
        alphaiw[0] = 119.6067900;
        alphaiw[1] = 117.3611800;
        alphaiw[2] = 111.4540000;
        alphaiw[3] = 95.2456700;
        alphaiw[4] = 79.7390500;
        alphaiw[5] = 66.9517800;
        alphaiw[6] = 56.7007000;
        alphaiw[7] = 48.4921200;
        alphaiw[8] = 41.8757900;
        alphaiw[9] = 36.4928900;
        alphaiw[10] = 32.0685600;
        alphaiw[11] = 28.3948300;
        alphaiw[12] = 22.7080700;
        alphaiw[13] = 18.5713200;
        alphaiw[14] = 15.4707700;
        alphaiw[15] = 13.0875400;
        alphaiw[16] = 11.2163200;
        alphaiw[17] = 7.9816700;
        alphaiw[18] = 5.9708100;
        alphaiw[19] = 3.7056100;
        alphaiw[20] = 2.5262500;
        alphaiw[21] = 1.2328500;
        alphaiw[22] = 0.7325600;
      }
      break;
    case 25:
      if (refsys == 1) {
        alphaiw[0] = 82.1201120;
        alphaiw[1] = 73.0648270;
        alphaiw[2] = 56.5578050;
        alphaiw[3] = 31.8457630;
        alphaiw[4] = 19.5848590;
        alphaiw[5] = 13.3420540;
        alphaiw[6] = 9.8197330;
        alphaiw[7] = 7.6416120;
        alphaiw[8] = 6.1931620;
        alphaiw[9] = 5.1733980;
        alphaiw[10] = 4.4221950;
        alphaiw[11] = 3.8482630;
        alphaiw[12] = 3.0320160;
        alphaiw[13] = 2.4803430;
        alphaiw[14] = 2.0824750;
        alphaiw[15] = 1.7818650;
        alphaiw[16] = 1.5468520;
        alphaiw[17] = 1.1371080;
        alphaiw[18] = 0.8757490;
        alphaiw[19] = 0.5688430;
        alphaiw[20] = 0.4003450;
        alphaiw[21] = 0.2041230;
        alphaiw[22] = 0.1236650;
      }
      else if (refsys == 2) {
        alphaiw[0] = 37.5942320;
        alphaiw[1] = 36.5102540;
        alphaiw[2] = 33.7065530;
        alphaiw[3] = 26.3778710;
        alphaiw[4] = 20.0194940;
        alphaiw[5] = 15.4093840;
        alphaiw[6] = 12.1689910;
        alphaiw[7] = 9.8609100;
        alphaiw[8] = 8.1746400;
        alphaiw[9] = 6.9089940;
        alphaiw[10] = 5.9348110;
        alphaiw[11] = 5.1678620;
        alphaiw[12] = 4.0486360;
        alphaiw[13] = 3.2798880;
        alphaiw[14] = 2.7243850;
        alphaiw[15] = 2.3069400;
        alphaiw[16] = 1.9834830;
        alphaiw[17] = 1.4290960;
        alphaiw[18] = 1.0840060;
        alphaiw[19] = 0.6895570;
        alphaiw[20] = 0.4788650;
        alphaiw[21] = 0.2397600;
        alphaiw[22] = 0.1439400;
      }
      else if (refsys == 3) {
        alphaiw[0] = 48.8093290;
        alphaiw[1] = 48.0980790;
        alphaiw[2] = 46.1748290;
        alphaiw[3] = 40.4031130;
        alphaiw[4] = 34.1760690;
        alphaiw[5] = 28.6218290;
        alphaiw[6] = 23.9924040;
        alphaiw[7] = 20.2296980;
        alphaiw[8] = 17.1934700;
        alphaiw[9] = 14.7402990;
        alphaiw[10] = 12.7472650;
        alphaiw[11] = 11.1157910;
        alphaiw[12] = 8.6475930;
        alphaiw[13] = 6.9086080;
        alphaiw[14] = 5.6437710;
        alphaiw[15] = 4.6974660;
        alphaiw[16] = 3.9718750;
        alphaiw[17] = 2.7571320;
        alphaiw[18] = 2.0290510;
        alphaiw[19] = 1.2349290;
        alphaiw[20] = 0.8324270;
        alphaiw[21] = 0.3994510;
        alphaiw[22] = 0.2345940;
      }
      break;
    case 26:
      if (refsys == 1) {
        alphaiw[0] = 65.1578920;
        alphaiw[1] = 57.4961000;
        alphaiw[2] = 43.3419800;
        alphaiw[3] = 23.8760970;
        alphaiw[4] = 15.1156130;
        alphaiw[5] = 10.7673250;
        alphaiw[6] = 8.2713540;
        alphaiw[7] = 6.6686050;
        alphaiw[8] = 5.5550430;
        alphaiw[9] = 4.7375840;
        alphaiw[10] = 4.1131870;
        alphaiw[11] = 3.6217190;
        alphaiw[12] = 2.8999580;
        alphaiw[13] = 2.3973090;
        alphaiw[14] = 2.0280480;
        alphaiw[15] = 1.7456620;
        alphaiw[16] = 1.5229730;
        alphaiw[17] = 1.1305990;
        alphaiw[18] = 0.8773050;
        alphaiw[19] = 0.5763700;
        alphaiw[20] = 0.4092270;
        alphaiw[21] = 0.2118590;
        alphaiw[22] = 0.1294110;
      }
      else if (refsys == 2) {
        alphaiw[0] = 34.5401790;
        alphaiw[1] = 33.5667420;
        alphaiw[2] = 31.1228640;
        alphaiw[3] = 24.7676090;
        alphaiw[4] = 19.0811590;
        alphaiw[5] = 14.8373870;
        alphaiw[6] = 11.7962700;
        alphaiw[7] = 9.6020150;
        alphaiw[8] = 7.9844480;
        alphaiw[9] = 6.7625510;
        alphaiw[10] = 5.8177240;
        alphaiw[11] = 5.0715270;
        alphaiw[12] = 3.9797330;
        alphaiw[13] = 3.2288080;
        alphaiw[14] = 2.6862250;
        alphaiw[15] = 2.2786720;
        alphaiw[16] = 1.9629540;
        alphaiw[17] = 1.4214720;
        alphaiw[18] = 1.0834860;
        alphaiw[19] = 0.6950870;
        alphaiw[20] = 0.4860810;
        alphaiw[21] = 0.2464070;
        alphaiw[22] = 0.1489560;
      }
      else if (refsys == 3) {
        alphaiw[0] = 45.9056800;
        alphaiw[1] = 45.2093150;
        alphaiw[2] = 43.3396340;
        alphaiw[3] = 37.7811330;
        alphaiw[4] = 31.8440890;
        alphaiw[5] = 26.6022340;
        alphaiw[6] = 22.2699010;
        alphaiw[7] = 18.7694260;
        alphaiw[8] = 15.9553430;
        alphaiw[9] = 13.6865850;
        alphaiw[10] = 11.8454300;
        alphaiw[11] = 10.3389430;
        alphaiw[12] = 8.0596090;
        alphaiw[13] = 6.4524730;
        alphaiw[14] = 5.2822830;
        alphaiw[15] = 4.4057560;
        alphaiw[16] = 3.7328380;
        alphaiw[17] = 2.6038570;
        alphaiw[18] = 1.9248540;
        alphaiw[19] = 1.1805840;
        alphaiw[20] = 0.8007680;
        alphaiw[21] = 0.3884140;
        alphaiw[22] = 0.2295320;
      }
      break;
    case 27:
      if (refsys == 1) {
        alphaiw[0] = 60.7910330;
        alphaiw[1] = 56.2385650;
        alphaiw[2] = 46.3237860;
        alphaiw[3] = 28.3000530;
        alphaiw[4] = 18.0784830;
        alphaiw[5] = 12.5506200;
        alphaiw[6] = 9.3329520;
        alphaiw[7] = 7.3060380;
        alphaiw[8] = 5.9416390;
        alphaiw[9] = 4.9729450;
        alphaiw[10] = 4.2552790;
        alphaiw[11] = 3.7050500;
        alphaiw[12] = 2.9212400;
        alphaiw[13] = 2.3925770;
        alphaiw[14] = 2.0131210;
        alphaiw[15] = 1.7278630;
        alphaiw[16] = 1.5057010;
        alphaiw[17] = 1.1191290;
        alphaiw[18] = 0.8715490;
        alphaiw[19] = 0.5771900;
        alphaiw[20] = 0.4126650;
        alphaiw[21] = 0.2164470;
        alphaiw[22] = 0.1333190;
      }
      else if (refsys == 2) {
        alphaiw[0] = 33.6766540;
        alphaiw[1] = 32.2795850;
        alphaiw[2] = 29.8110510;
        alphaiw[3] = 23.6860370;
        alphaiw[4] = 18.3057380;
        alphaiw[5] = 14.2980140;
        alphaiw[6] = 11.4150610;
        alphaiw[7] = 9.3244580;
        alphaiw[8] = 7.7757940;
        alphaiw[9] = 6.6008390;
        alphaiw[10] = 5.6888800;
        alphaiw[11] = 4.9663640;
        alphaiw[12] = 3.9055760;
        alphaiw[13] = 3.1737960;
        alphaiw[14] = 2.6443420;
        alphaiw[15] = 2.2465040;
        alphaiw[16] = 1.9383190;
        alphaiw[17] = 1.4096790;
        alphaiw[18] = 1.0792130;
        alphaiw[19] = 0.6978430;
        alphaiw[20] = 0.4912380;
        alphaiw[21] = 0.2521090;
        alphaiw[22] = 0.1535920;
      }
      else if (refsys == 3) {
        alphaiw[0] = 33.0681610;
        alphaiw[1] = 31.8814040;
        alphaiw[2] = 29.3779300;
        alphaiw[3] = 23.6057900;
        alphaiw[4] = 18.6463550;
        alphaiw[5] = 14.9145580;
        alphaiw[6] = 12.1627890;
        alphaiw[7] = 10.1071300;
        alphaiw[8] = 8.5395170;
        alphaiw[9] = 7.3190970;
        alphaiw[10] = 6.3509310;
        alphaiw[11] = 5.5699740;
        alphaiw[12] = 4.4005400;
        alphaiw[13] = 3.5786790;
        alphaiw[14] = 2.9775460;
        alphaiw[15] = 2.5232770;
        alphaiw[16] = 2.1705940;
        alphaiw[17] = 1.5662690;
        alphaiw[18] = 1.1907560;
        alphaiw[19] = 0.7619390;
        alphaiw[20] = 0.5325880;
        alphaiw[21] = 0.2706440;
        alphaiw[22] = 0.1640920;
      }
      else if (refsys == 4) {
        alphaiw[0] = 44.7198600;
        alphaiw[1] = 43.9726200;
        alphaiw[2] = 41.9679160;
        alphaiw[3] = 36.1202890;
        alphaiw[4] = 30.0864670;
        alphaiw[5] = 24.9191710;
        alphaiw[6] = 20.7437800;
        alphaiw[7] = 17.4226860;
        alphaiw[8] = 14.7811930;
        alphaiw[9] = 12.6668970;
        alphaiw[10] = 10.9594050;
        alphaiw[11] = 9.5668070;
        alphaiw[12] = 7.4655900;
        alphaiw[13] = 5.9866860;
        alphaiw[14] = 4.9103160;
        alphaiw[15] = 4.1038230;
        alphaiw[16] = 3.4842290;
        alphaiw[17] = 2.4429010;
        alphaiw[18] = 1.8145610;
        alphaiw[19] = 1.1222760;
        alphaiw[20] = 0.7664210;
        alphaiw[21] = 0.3762730;
        alphaiw[22] = 0.2239990;
      }
      break;
    case 28:
      if (refsys == 1) {
        alphaiw[0] = 65.1471670;
        alphaiw[1] = 59.7497700;
        alphaiw[2] = 48.4286590;
        alphaiw[3] = 28.5849020;
        alphaiw[4] = 17.8269450;
        alphaiw[5] = 12.1955640;
        alphaiw[6] = 8.9866060;
        alphaiw[7] = 6.9942690;
        alphaiw[8] = 5.6675310;
        alphaiw[9] = 4.7335770;
        alphaiw[10] = 4.0463770;
        alphaiw[11] = 3.5223700;
        alphaiw[12] = 2.7801490;
        alphaiw[13] = 2.2820840;
        alphaiw[14] = 1.9256220;
        alphaiw[15] = 1.6580700;
        alphaiw[16] = 1.4498280;
        alphaiw[17] = 1.0871500;
        alphaiw[18] = 0.8537840;
        alphaiw[19] = 0.5733100;
        alphaiw[20] = 0.4140840;
        alphaiw[21] = 0.2208130;
        alphaiw[22] = 0.1373770;
      }
      else if (refsys == 2) {
        alphaiw[0] = 33.2010460;
        alphaiw[1] = 31.9092910;
        alphaiw[2] = 29.3670970;
        alphaiw[3] = 23.1912300;
        alphaiw[4] = 17.8686070;
        alphaiw[5] = 13.9453180;
        alphaiw[6] = 11.1368280;
        alphaiw[7] = 9.1041990;
        alphaiw[8] = 7.5992180;
        alphaiw[9] = 6.4570670;
        alphaiw[10] = 5.5699260;
        alphaiw[11] = 4.8664170;
        alphaiw[12] = 3.8320980;
        alphaiw[13] = 3.1174810;
        alphaiw[14] = 2.6000250;
        alphaiw[15] = 2.2111480;
        alphaiw[16] = 1.9099840;
        alphaiw[17] = 1.3937010;
        alphaiw[18] = 1.0709720;
        alphaiw[19] = 0.6975870;
        alphaiw[20] = 0.4941560;
        alphaiw[21] = 0.2566390;
        alphaiw[22] = 0.1576070;
      }
      else if (refsys == 3) {
        alphaiw[0] = 43.8557940;
        alphaiw[1] = 42.9928330;
        alphaiw[2] = 40.7502500;
        alphaiw[3] = 34.4926310;
        alphaiw[4] = 28.3244810;
        alphaiw[5] = 23.2253070;
        alphaiw[6] = 19.2045670;
        alphaiw[7] = 16.0592920;
        alphaiw[8] = 13.5864050;
        alphaiw[9] = 11.6233040;
        alphaiw[10] = 10.0473580;
        alphaiw[11] = 8.7676580;
        alphaiw[12] = 6.8450210;
        alphaiw[13] = 5.4965520;
        alphaiw[14] = 4.5167510;
        alphaiw[15] = 3.7830110;
        alphaiw[16] = 3.2192030;
        alphaiw[17] = 2.2703070;
        alphaiw[18] = 1.6957880;
        alphaiw[19] = 1.0589930;
        alphaiw[20] = 0.7288230;
        alphaiw[21] = 0.3627090;
        alphaiw[22] = 0.2177520;
      }
      else if (refsys == 4) {
        alphaiw[0] = 46.9971430;
        alphaiw[1] = 45.7683980;
        alphaiw[2] = 42.6546550;
        alphaiw[3] = 34.8257600;
        alphaiw[4] = 28.0599290;
        alphaiw[5] = 22.8623880;
        alphaiw[6] = 18.8881430;
        alphaiw[7] = 15.8119240;
        alphaiw[8] = 13.3988200;
        alphaiw[9] = 11.4815020;
        alphaiw[10] = 9.9393620;
        alphaiw[11] = 8.6844580;
        alphaiw[12] = 6.7937410;
        alphaiw[13] = 5.4634420;
        alphaiw[14] = 4.4945350;
        alphaiw[15] = 3.7676420;
        alphaiw[16] = 3.2083130;
        alphaiw[17] = 2.2653470;
        alphaiw[18] = 1.6933670;
        alphaiw[19] = 1.0583620;
        alphaiw[20] = 0.7286770;
        alphaiw[21] = 0.3627700;
        alphaiw[22] = 0.2178100;
      }
      break;
    case 29:
      if (refsys == 1) {
        alphaiw[0] = 51.7659850;
        alphaiw[1] = 46.4291690;
        alphaiw[2] = 36.0204980;
        alphaiw[3] = 20.5952700;
        alphaiw[4] = 13.2155350;
        alphaiw[5] = 9.4590740;
        alphaiw[6] = 7.2960120;
        alphaiw[7] = 5.9145750;
        alphaiw[8] = 4.9593020;
        alphaiw[9] = 4.2583130;
        alphaiw[10] = 3.7208220;
        alphaiw[11] = 3.2950100;
        alphaiw[12] = 2.6629290;
        alphaiw[13] = 2.2171860;
        alphaiw[14] = 1.8871890;
        alphaiw[15] = 1.6338740;
        alphaiw[16] = 1.4337810;
        alphaiw[17] = 1.0804860;
        alphaiw[18] = 0.8509580;
        alphaiw[19] = 0.5739100;
        alphaiw[20] = 0.4162210;
        alphaiw[21] = 0.2239360;
        alphaiw[22] = 0.1401690;
      }
      else if (refsys == 2) {
        alphaiw[0] = 32.2640590;
        alphaiw[1] = 30.7643190;
        alphaiw[2] = 27.3346460;
        alphaiw[3] = 19.9931860;
        alphaiw[4] = 14.6605450;
        alphaiw[5] = 11.1745880;
        alphaiw[6] = 8.8578670;
        alphaiw[7] = 7.2504740;
        alphaiw[8] = 6.0856020;
        alphaiw[9] = 5.2091410;
        alphaiw[10] = 4.5290430;
        alphaiw[11] = 3.9879440;
        alphaiw[12] = 3.1856560;
        alphaiw[13] = 2.6238350;
        alphaiw[14] = 2.2117070;
        alphaiw[15] = 1.8983580;
        alphaiw[16] = 1.6531670;
        alphaiw[17] = 1.2266060;
        alphaiw[18] = 0.9549290;
        alphaiw[19] = 0.6340400;
        alphaiw[20] = 0.4552650;
        alphaiw[21] = 0.2415160;
        alphaiw[22] = 0.1500870;
      }
      break;
    case 30:
      if (refsys == 1) {
        alphaiw[0] = 44.8764280;
        alphaiw[1] = 42.4015380;
        alphaiw[2] = 36.5011020;
        alphaiw[3] = 24.0580460;
        alphaiw[4] = 15.9961690;
        alphaiw[5] = 11.3275720;
        alphaiw[6] = 8.5125980;
        alphaiw[7] = 6.7050110;
        alphaiw[8] = 5.4750510;
        alphaiw[9] = 4.5962800;
        alphaiw[10] = 3.9426060;
        alphaiw[11] = 3.4399690;
        alphaiw[12] = 2.7216070;
        alphaiw[13] = 2.2355020;
        alphaiw[14] = 1.8860970;
        alphaiw[15] = 1.6235520;
        alphaiw[16] = 1.4194560;
        alphaiw[17] = 1.0657210;
        alphaiw[18] = 0.8398960;
        alphaiw[19] = 0.5698780;
        alphaiw[20] = 0.4161790;
        alphaiw[21] = 0.2271770;
        alphaiw[22] = 0.1438000;
      }
      else if (refsys == 2) {
        alphaiw[0] = 29.2304490;
        alphaiw[1] = 28.6434220;
        alphaiw[2] = 27.0389880;
        alphaiw[3] = 22.2870810;
        alphaiw[4] = 17.5546960;
        alphaiw[5] = 13.8059590;
        alphaiw[6] = 11.0322270;
        alphaiw[7] = 8.9976380;
        alphaiw[8] = 7.4858160;
        alphaiw[9] = 6.3399040;
        alphaiw[10] = 5.4528720;
        alphaiw[11] = 4.7523760;
        alphaiw[12] = 3.7286630;
        alphaiw[13] = 3.0264940;
        alphaiw[14] = 2.5210460;
        alphaiw[15] = 2.1430610;
        alphaiw[16] = 1.8515770;
        alphaiw[17] = 1.3547510;
        alphaiw[18] = 1.0459910;
        alphaiw[19] = 0.6893360;
        alphaiw[20] = 0.4939420;
        alphaiw[21] = 0.2624940;
        alphaiw[22] = 0.1638640;
      }
      break;
    case 31:
      if (refsys == 1) {
        alphaiw[0] = 56.2798430;
        alphaiw[1] = 51.9167870;
        alphaiw[2] = 43.2310870;
        alphaiw[3] = 28.1416200;
        alphaiw[4] = 18.8462820;
        alphaiw[5] = 13.2994650;
        alphaiw[6] = 9.8600450;
        alphaiw[7] = 7.6193980;
        alphaiw[8] = 6.0894620;
        alphaiw[9] = 5.0007640;
        alphaiw[10] = 4.1981820;
        alphaiw[11] = 3.5884760;
        alphaiw[12] = 2.7347440;
        alphaiw[13] = 2.1738840;
        alphaiw[14] = 1.7818110;
        alphaiw[15] = 1.4944090;
        alphaiw[16] = 1.2758120;
        alphaiw[17] = 0.9086190;
        alphaiw[18] = 0.6839730;
        alphaiw[19] = 0.4300690;
        alphaiw[20] = 0.2959040;
        alphaiw[21] = 0.1458990;
        alphaiw[22] = 0.0867970;
      }
      else if (refsys == 2) {
        alphaiw[0] = 53.3046110;
        alphaiw[1] = 49.9067250;
        alphaiw[2] = 42.6926500;
        alphaiw[3] = 29.2253470;
        alphaiw[4] = 20.4243950;
        alphaiw[5] = 14.8882280;
        alphaiw[6] = 11.2916180;
        alphaiw[7] = 8.8592270;
        alphaiw[8] = 7.1500550;
        alphaiw[9] = 5.9073300;
        alphaiw[10] = 4.9763880;
        alphaiw[11] = 4.2607350;
        alphaiw[12] = 3.2471250;
        alphaiw[13] = 2.5751450;
        alphaiw[14] = 2.1036730;
        alphaiw[15] = 1.7579590;
        alphaiw[16] = 1.4954520;
        alphaiw[17] = 1.0567390;
        alphaiw[18] = 0.7907280;
        alphaiw[19] = 0.4933050;
        alphaiw[20] = 0.3378040;
        alphaiw[21] = 0.1654770;
        alphaiw[22] = 0.0980990;
      }
      else if (refsys == 3) {
        alphaiw[0] = 40.6953470;
        alphaiw[1] = 38.6598060;
        alphaiw[2] = 34.9271010;
        alphaiw[3] = 27.0292080;
        alphaiw[4] = 20.5104160;
        alphaiw[5] = 15.7311360;
        alphaiw[6] = 12.3171900;
        alphaiw[7] = 9.8626390;
        alphaiw[8] = 8.0655210;
        alphaiw[9] = 6.7210170;
        alphaiw[10] = 5.6931530;
        alphaiw[11] = 4.8912630;
        alphaiw[12] = 3.7391050;
        alphaiw[13] = 2.9659310;
        alphaiw[14] = 2.4200560;
        alphaiw[15] = 2.0186870;
        alphaiw[16] = 1.7137410;
        alphaiw[17] = 1.2049770;
        alphaiw[18] = 0.8979620;
        alphaiw[19] = 0.5570650;
        alphaiw[20] = 0.3801550;
        alphaiw[21] = 0.1853670;
        alphaiw[22] = 0.1096290;
      }
      break;
    case 32:
      if (refsys == 1) {
        alphaiw[0] = 41.6184280;
        alphaiw[1] = 39.9833160;
        alphaiw[2] = 35.9872960;
        alphaiw[3] = 26.6103970;
        alphaiw[4] = 19.2286830;
        alphaiw[5] = 14.1743150;
        alphaiw[6] = 10.7638470;
        alphaiw[7] = 8.4179500;
        alphaiw[8] = 6.7587960;
        alphaiw[9] = 5.5514820;
        alphaiw[10] = 4.6493180;
        alphaiw[11] = 3.9588270;
        alphaiw[12] = 2.9887560;
        alphaiw[13] = 2.3535080;
        alphaiw[14] = 1.9131390;
        alphaiw[15] = 1.5937360;
        alphaiw[16] = 1.3534810;
        alphaiw[17] = 0.9566170;
        alphaiw[18] = 0.7184900;
        alphaiw[19] = 0.4529630;
        alphaiw[20] = 0.3132500;
        alphaiw[21] = 0.1561230;
        alphaiw[22] = 0.0935260;
      }
      else if (refsys == 2) {
        alphaiw[0] = 40.4520440;
        alphaiw[1] = 39.0264160;
        alphaiw[2] = 35.5809840;
        alphaiw[3] = 27.3666530;
        alphaiw[4] = 20.5342330;
        alphaiw[5] = 15.5791490;
        alphaiw[6] = 12.0769830;
        alphaiw[7] = 9.5821540;
        alphaiw[8] = 7.7708270;
        alphaiw[9] = 6.4264770;
        alphaiw[10] = 5.4066710;
        alphaiw[11] = 4.6170600;
        alphaiw[12] = 3.4943900;
        alphaiw[13] = 2.7511970;
        alphaiw[14] = 2.2328710;
        alphaiw[15] = 1.8558220;
        alphaiw[16] = 1.5719620;
        alphaiw[17] = 1.1037790;
        alphaiw[18] = 0.8243160;
        alphaiw[19] = 0.5154070;
        alphaiw[20] = 0.3545620;
        alphaiw[21] = 0.1754620;
        alphaiw[22] = 0.1047250;
      }
      else if (refsys == 3) {
        alphaiw[0] = 40.4418220;
        alphaiw[1] = 39.2545360;
        alphaiw[2] = 36.2478860;
        alphaiw[3] = 28.6222220;
        alphaiw[4] = 21.9792320;
        alphaiw[5] = 16.9925560;
        alphaiw[6] = 13.3626540;
        alphaiw[7] = 10.7138950;
        alphaiw[8] = 8.7539710;
        alphaiw[9] = 7.2776830;
        alphaiw[10] = 6.1448790;
        alphaiw[11] = 5.2599430;
        alphaiw[12] = 3.9900180;
        alphaiw[13] = 3.1421830;
        alphaiw[14] = 2.5480020;
        alphaiw[15] = 2.1146930;
        alphaiw[16] = 1.7881720;
        alphaiw[17] = 1.2499950;
        alphaiw[18] = 0.9298180;
        alphaiw[19] = 0.5780120;
        alphaiw[20] = 0.3961690;
        alphaiw[21] = 0.1950880;
        alphaiw[22] = 0.1161420;
      }
      else if (refsys == 4) {
        alphaiw[0] = 37.6716300;
        alphaiw[1] = 36.8525860;
        alphaiw[2] = 34.6926590;
        alphaiw[3] = 28.6826580;
        alphaiw[4] = 22.8485090;
        alphaiw[5] = 18.1179570;
        alphaiw[6] = 14.4936950;
        alphaiw[7] = 11.7568160;
        alphaiw[8] = 9.6831810;
        alphaiw[9] = 8.0946770;
        alphaiw[10] = 6.8606380;
        alphaiw[11] = 5.8877070;
        alphaiw[12] = 4.4784020;
        alphaiw[13] = 3.5295510;
        alphaiw[14] = 2.8613170;
        alphaiw[15] = 2.3726930;
        alphaiw[16] = 2.0040290;
        alphaiw[17] = 1.3963730;
        alphaiw[18] = 1.0356150;
        alphaiw[19] = 0.6409590;
        alphaiw[20] = 0.4381220;
        alphaiw[21] = 0.2150380;
        alphaiw[22] = 0.1278320;
      }
      else if (refsys == 5) {
        alphaiw[0] = 34.9443880;
        alphaiw[1] = 34.4589550;
        alphaiw[2] = 33.0932300;
        alphaiw[3] = 28.6849650;
        alphaiw[4] = 23.6815050;
        alphaiw[5] = 19.2193090;
        alphaiw[6] = 15.6074480;
        alphaiw[7] = 12.7869910;
        alphaiw[8] = 10.6029230;
        alphaiw[9] = 8.9045950;
        alphaiw[10] = 7.5710690;
        alphaiw[11] = 6.5114240;
        alphaiw[12] = 4.9643950;
        alphaiw[13] = 3.9154690;
        alphaiw[14] = 3.1737400;
        alphaiw[15] = 2.6301480;
        alphaiw[16] = 2.2195640;
        alphaiw[17] = 1.5427250;
        alphaiw[18] = 1.1415220;
        alphaiw[19] = 0.7041370;
        alphaiw[20] = 0.4803530;
        alphaiw[21] = 0.2352710;
        alphaiw[22] = 0.1397600;
      }
      break;
    case 33:
      if (refsys == 1) {
        alphaiw[0] = 30.9622720;
        alphaiw[1] = 30.2682480;
        alphaiw[2] = 28.4117290;
        alphaiw[3] = 23.1419060;
        alphaiw[4] = 18.0459760;
        alphaiw[5] = 14.0238200;
        alphaiw[6] = 11.0354630;
        alphaiw[7] = 8.8374630;
        alphaiw[8] = 7.2065030;
        alphaiw[9] = 5.9774600;
        alphaiw[10] = 5.0352070;
        alphaiw[11] = 4.3004040;
        alphaiw[12] = 3.2496190;
        alphaiw[13] = 2.5521200;
        alphaiw[14] = 2.0662400;
        alphaiw[15] = 1.7139300;
        alphaiw[16] = 1.4497910;
        alphaiw[17] = 1.0172000;
        alphaiw[18] = 0.7612360;
        alphaiw[19] = 0.4798100;
        alphaiw[20] = 0.3329600;
        alphaiw[21] = 0.1673800;
        alphaiw[22] = 0.1007870;
      }
      else if (refsys == 2) {
        alphaiw[0] = 32.9857230;
        alphaiw[1] = 32.2738720;
        alphaiw[2] = 30.3903080;
        alphaiw[3] = 25.0725470;
        alphaiw[4] = 19.8545550;
        alphaiw[5] = 15.6403290;
        alphaiw[6] = 12.4395230;
        alphaiw[7] = 10.0416420;
        alphaiw[8] = 8.2363030;
        alphaiw[9] = 6.8602630;
        alphaiw[10] = 5.7958010;
        alphaiw[11] = 4.9597520;
        alphaiw[12] = 3.7548970;
        alphaiw[13] = 2.9491190;
        alphaiw[14] = 2.3851990;
        alphaiw[15] = 1.9752150;
        alphaiw[16] = 1.6674560;
        alphaiw[17] = 1.1635310;
        alphaiw[18] = 0.8662610;
        alphaiw[19] = 0.5416150;
        alphaiw[20] = 0.3738430;
        alphaiw[21] = 0.1866530;
        alphaiw[22] = 0.1120540;
      }
      else if (refsys == 3) {
        alphaiw[0] = 34.5403260;
        alphaiw[1] = 33.8506380;
        alphaiw[2] = 32.0374310;
        alphaiw[3] = 26.8393630;
        alphaiw[4] = 21.5730010;
        alphaiw[5] = 17.1996430;
        alphaiw[6] = 13.8049410;
        alphaiw[7] = 11.2190870;
        alphaiw[8] = 9.2473150;
        alphaiw[9] = 7.7296440;
        alphaiw[10] = 6.5466390;
        alphaiw[11] = 5.6118890;
        alphaiw[12] = 4.2560500;
        alphaiw[13] = 3.3436390;
        alphaiw[14] = 2.7026150;
        alphaiw[15] = 2.2355180;
        alphaiw[16] = 1.8844910;
        alphaiw[17] = 1.3096800;
        alphaiw[18] = 0.9712870;
        alphaiw[19] = 0.6035540;
        alphaiw[20] = 0.4148920;
        alphaiw[21] = 0.2060850;
        alphaiw[22] = 0.1234480;
      }
      else if (refsys == 4) {
        alphaiw[0] = 35.5626970;
        alphaiw[1] = 34.9738250;
        alphaiw[2] = 33.3553240;
        alphaiw[3] = 28.4265670;
        alphaiw[4] = 23.1867600;
        alphaiw[5] = 18.6947170;
        alphaiw[6] = 15.1290090;
        alphaiw[7] = 12.3689250;
        alphaiw[8] = 10.2393320;
        alphaiw[9] = 8.5855950;
        alphaiw[10] = 7.2877370;
        alphaiw[11] = 6.2567940;
        alphaiw[12] = 4.7529590;
        alphaiw[13] = 3.7354930;
        alphaiw[14] = 3.0182630;
        alphaiw[15] = 2.4945970;
        alphaiw[16] = 2.1006490;
        alphaiw[17] = 1.4553950;
        alphaiw[18] = 1.0760650;
        alphaiw[19] = 0.6653780;
        alphaiw[20] = 0.4558740;
        alphaiw[21] = 0.2254990;
        alphaiw[22] = 0.1348450;
      }
      break;
    case 34:
      if (refsys == 1) {
        alphaiw[0] = 25.7234860;
        alphaiw[1] = 25.2967560;
        alphaiw[2] = 24.1266960;
        alphaiw[3] = 20.5754680;
        alphaiw[4] = 16.8021780;
        alphaiw[5] = 13.5652820;
        alphaiw[6] = 10.9952140;
        alphaiw[7] = 9.0046390;
        alphaiw[8] = 7.4667080;
        alphaiw[9] = 6.2701390;
        alphaiw[10] = 5.3290860;
        alphaiw[11] = 4.5800740;
        alphaiw[12] = 3.4849480;
        alphaiw[13] = 2.7424650;
        alphaiw[14] = 2.2186000;
        alphaiw[15] = 1.8360660;
        alphaiw[16] = 1.5483750;
        alphaiw[17] = 1.0775180;
        alphaiw[18] = 0.8008720;
        alphaiw[19] = 0.5007840;
        alphaiw[20] = 0.3466260;
        alphaiw[21] = 0.1746160;
        alphaiw[22] = 0.1055700;
      }
      else if (refsys == 2) {
        alphaiw[0] = 28.5144970;
        alphaiw[1] = 28.0475240;
        alphaiw[2] = 26.7740720;
        alphaiw[3] = 22.9293920;
        alphaiw[4] = 18.8329680;
        alphaiw[5] = 15.2899050;
        alphaiw[6] = 12.4503800;
        alphaiw[7] = 10.2324530;
        alphaiw[8] = 8.5069500;
        alphaiw[9] = 7.1570530;
        alphaiw[10] = 6.0907950;
        alphaiw[11] = 5.2391750;
        alphaiw[12] = 3.9892260;
        alphaiw[13] = 3.1385090;
        alphaiw[14] = 2.5367660;
        alphaiw[15] = 2.0966900;
        alphaiw[16] = 1.7654690;
        alphaiw[17] = 1.2233860;
        alphaiw[18] = 0.9054680;
        alphaiw[19] = 0.5621910;
        alphaiw[20] = 0.3871530;
        alphaiw[21] = 0.1936670;
        alphaiw[22] = 0.1167080;
      }
      else if (refsys == 3) {
        alphaiw[0] = 30.9890680;
        alphaiw[1] = 30.5032140;
        alphaiw[2] = 29.1715810;
        alphaiw[3] = 25.1156940;
        alphaiw[4] = 20.7515920;
        alphaiw[5] = 16.9397500;
        alphaiw[6] = 13.8556780;
        alphaiw[7] = 11.4267060;
        alphaiw[8] = 9.5242290;
        alphaiw[9] = 8.0279190;
        alphaiw[10] = 6.8410310;
        alphaiw[11] = 5.8898910;
        alphaiw[12] = 4.4887500;
        alphaiw[13] = 3.5316780;
        alphaiw[14] = 2.8531070;
        alphaiw[15] = 2.3561140;
        alphaiw[16] = 1.9817540;
        alphaiw[17] = 1.3689360;
        alphaiw[18] = 1.0099510;
        alphaiw[19] = 0.6236450;
        alphaiw[20] = 0.4277790;
        alphaiw[21] = 0.2128400;
        alphaiw[22] = 0.1279540;
      }
      break;
    case 35:
      if (refsys == 1) {
        alphaiw[0] = 20.8160780;
        alphaiw[1] = 20.5760200;
        alphaiw[2] = 19.8987130;
        alphaiw[3] = 17.6813900;
        alphaiw[4] = 15.0736970;
        alphaiw[5] = 12.6318500;
        alphaiw[6] = 10.5528970;
        alphaiw[7] = 8.8524810;
        alphaiw[8] = 7.4812900;
        alphaiw[9] = 6.3774160;
        alphaiw[10] = 5.4848940;
        alphaiw[11] = 4.7581530;
        alphaiw[12] = 3.6676360;
        alphaiw[13] = 2.9082370;
        alphaiw[14] = 2.3624100;
        alphaiw[15] = 1.9587510;
        alphaiw[16] = 1.6526030;
        alphaiw[17] = 1.1479460;
        alphaiw[18] = 0.8506710;
        alphaiw[19] = 0.5297310;
        alphaiw[20] = 0.3663790;
        alphaiw[21] = 0.1854320;
        alphaiw[22] = 0.1127450;
      }
      else if (refsys == 2) {
        alphaiw[0] = 23.7818040;
        alphaiw[1] = 23.5043600;
        alphaiw[2] = 22.7232600;
        alphaiw[3] = 20.1781530;
        alphaiw[4] = 17.1975750;
        alphaiw[5] = 14.4112140;
        alphaiw[6] = 12.0390110;
        alphaiw[7] = 10.0975330;
        alphaiw[8] = 8.5309690;
        alphaiw[9] = 7.2693050;
        alphaiw[10] = 6.2490260;
        alphaiw[11] = 5.4182350;
        alphaiw[12] = 4.1716820;
        alphaiw[13] = 3.3037540;
        alphaiw[14] = 2.6800190;
        alphaiw[15] = 2.2188460;
        alphaiw[16] = 1.8691920;
        alphaiw[17] = 1.2933190;
        alphaiw[18] = 0.9547540;
        alphaiw[19] = 0.5906060;
        alphaiw[20] = 0.4063970;
        alphaiw[21] = 0.2041110;
        alphaiw[22] = 0.1236350;
      }
      break;
    case 36:
      if (refsys == 1) {
        alphaiw[0] = 16.8771250;
        alphaiw[1] = 16.7404080;
        alphaiw[2] = 16.3473440;
        alphaiw[3] = 14.9868570;
        alphaiw[4] = 13.2465280;
        alphaiw[5] = 11.4802080;
        alphaiw[6] = 9.8697990;
        alphaiw[7] = 8.4773240;
        alphaiw[8] = 7.3035050;
        alphaiw[9] = 6.3243740;
        alphaiw[10] = 5.5096050;
        alphaiw[11] = 4.8302660;
        alphaiw[12] = 3.7825130;
        alphaiw[13] = 3.0313660;
        alphaiw[14] = 2.4797800;
        alphaiw[15] = 2.0653280;
        alphaiw[16] = 1.7472760;
        alphaiw[17] = 1.2164410;
        alphaiw[18] = 0.9007200;
        alphaiw[19] = 0.5591630;
        alphaiw[20] = 0.3861250;
        alphaiw[21] = 0.1959080;
        alphaiw[22] = 0.1197080;
      }
      break;
    case 37:
      if (refsys == 1) {
        alphaiw[0] = 410.1570650;
        alphaiw[1] = 214.3629420;
        alphaiw[2] = 92.2189500;
        alphaiw[3] = 33.0762760;
        alphaiw[4] = 19.1702550;
        alphaiw[5] = 13.6453490;
        alphaiw[6] = 10.6953320;
        alphaiw[7] = 8.8115010;
        alphaiw[8] = 7.4699530;
        alphaiw[9] = 6.4487380;
        alphaiw[10] = 5.6381240;
        alphaiw[11] = 4.9766630;
        alphaiw[12] = 3.9619960;
        alphaiw[13] = 3.2244460;
        alphaiw[14] = 2.6697510;
        alphaiw[15] = 2.2421050;
        alphaiw[16] = 1.9059170;
        alphaiw[17] = 1.3262900;
        alphaiw[18] = 0.9698860;
        alphaiw[19] = 0.5775830;
        alphaiw[20] = 0.3804570;
        alphaiw[21] = 0.1743570;
        alphaiw[22] = 0.0992380;
      }
      else if (refsys == 2) {
        alphaiw[0] = 62.1046880;
        alphaiw[1] = 55.3737420;
        alphaiw[2] = 43.5682870;
        alphaiw[3] = 27.4983580;
        alphaiw[4] = 19.5497690;
        alphaiw[5] = 15.0663500;
        alphaiw[6] = 12.1785660;
        alphaiw[7] = 10.1413920;
        alphaiw[8] = 8.6158260;
        alphaiw[9] = 7.4268580;
        alphaiw[10] = 6.4740940;
        alphaiw[11] = 5.6949060;
        alphaiw[12] = 4.5031460;
        alphaiw[13] = 3.6435490;
        alphaiw[14] = 3.0024770;
        alphaiw[15] = 2.5119750;
        alphaiw[16] = 2.1288570;
        alphaiw[17] = 1.4737190;
        alphaiw[18] = 1.0743370;
        alphaiw[19] = 0.6377080;
        alphaiw[20] = 0.4194260;
        alphaiw[21] = 0.1919380;
        alphaiw[22] = 0.1091900;
      }
      else if (refsys == 3) {
        alphaiw[0] = 126.6310029;
        alphaiw[1] = 125.1559058;
        alphaiw[2] = 121.0742788;
        alphaiw[3] = 108.3482901;
        alphaiw[4] = 94.1247762;
        alphaiw[5] = 81.0337905;
        alphaiw[6] = 69.7439741;
        alphaiw[7] = 60.2269348;
        alphaiw[8] = 52.2628997;
        alphaiw[9] = 45.6040488;
        alphaiw[10] = 40.0235464;
        alphaiw[11] = 35.3278801;
        alphaiw[12] = 27.9808202;
        alphaiw[13] = 22.6114560;
        alphaiw[14] = 18.5992009;
        alphaiw[15] = 15.5379934;
        alphaiw[16] = 13.1579386;
        alphaiw[17] = 9.1200087;
        alphaiw[18] = 6.6812408;
        alphaiw[19] = 4.0309298;
        alphaiw[20] = 2.7053916;
        alphaiw[21] = 1.3052810;
        alphaiw[22] = 0.7781749;
      }
      break;
    case 38:
      if (refsys == 1) {
        alphaiw[0] = 207.0959830;
        alphaiw[1] = 167.4062640;
        alphaiw[2] = 107.5395190;
        alphaiw[3] = 46.8864590;
        alphaiw[4] = 26.2027860;
        alphaiw[5] = 17.4151700;
        alphaiw[6] = 12.8603380;
        alphaiw[7] = 10.1317340;
        alphaiw[8] = 8.3169810;
        alphaiw[9] = 7.0156000;
        alphaiw[10] = 6.0305450;
        alphaiw[11] = 5.2554110;
        alphaiw[12] = 4.1095420;
        alphaiw[13] = 3.3041900;
        alphaiw[14] = 2.7115560;
        alphaiw[15] = 2.2617260;
        alphaiw[16] = 1.9123040;
        alphaiw[17] = 1.3185610;
        alphaiw[18] = 0.9589950;
        alphaiw[19] = 0.5679990;
        alphaiw[20] = 0.3733060;
        alphaiw[21] = 0.1708320;
        alphaiw[22] = 0.0972390;
      }
      else if (refsys == 2) {
        alphaiw[0] = 203.7569620;
        alphaiw[1] = 136.2520380;
        alphaiw[2] = 79.3342830;
        alphaiw[3] = 39.0709150;
        alphaiw[4] = 24.8495540;
        alphaiw[5] = 17.8598820;
        alphaiw[6] = 13.7760240;
        alphaiw[7] = 11.1115300;
        alphaiw[8] = 9.2347870;
        alphaiw[9] = 7.8382880;
        alphaiw[10] = 6.7567000;
        alphaiw[11] = 5.8938450;
        alphaiw[12] = 4.6058940;
        alphaiw[13] = 3.6968050;
        alphaiw[14] = 3.0281860;
        alphaiw[15] = 2.5217330;
        alphaiw[16] = 2.1292640;
        alphaiw[17] = 1.4647060;
        alphaiw[18] = 1.0638290;
        alphaiw[19] = 0.6292740;
        alphaiw[20] = 0.4133770;
        alphaiw[21] = 0.1891240;
        alphaiw[22] = 0.1076560;
      }
      else if (refsys == 3) {
        alphaiw[0] = 56.7366690;
        alphaiw[1] = 53.8312150;
        alphaiw[2] = 47.2695200;
        alphaiw[3] = 33.9461020;
        alphaiw[4] = 24.8182410;
        alphaiw[5] = 18.9521460;
        alphaiw[6] = 15.0374310;
        alphaiw[7] = 12.2916810;
        alphaiw[8] = 10.2772390;
        alphaiw[9] = 8.7433090;
        alphaiw[10] = 7.5399850;
        alphaiw[11] = 6.5735100;
        alphaiw[12] = 5.1255840;
        alphaiw[13] = 4.1035260;
        alphaiw[14] = 3.3535750;
        alphaiw[15] = 2.7871740;
        alphaiw[16] = 2.3494990;
        alphaiw[17] = 1.6113550;
        alphaiw[18] = 1.1681150;
        alphaiw[19] = 0.6895060;
        alphaiw[20] = 0.4524650;
        alphaiw[21] = 0.2067690;
        alphaiw[22] = 0.1176430;
      }
      else if (refsys == 4) {
        alphaiw[0] = 116.4505849;
        alphaiw[1] = 114.1968844;
        alphaiw[2] = 108.2963555;
        alphaiw[3] = 92.1203844;
        alphaiw[4] = 76.7209726;
        alphaiw[5] = 64.3033634;
        alphaiw[6] = 54.5977463;
        alphaiw[7] = 46.9458079;
        alphaiw[8] = 40.8022639;
        alphaiw[9] = 35.7806868;
        alphaiw[10] = 31.6136150;
        alphaiw[11] = 28.1128216;
        alphaiw[12] = 22.5997820;
        alphaiw[13] = 18.5077567;
        alphaiw[14] = 15.3957648;
        alphaiw[15] = 12.9813495;
        alphaiw[16] = 11.0757782;
        alphaiw[17] = 7.7757352;
        alphaiw[18] = 5.7345816;
        alphaiw[19] = 3.4676638;
        alphaiw[20] = 2.3132969;
        alphaiw[21] = 1.0842547;
        alphaiw[22] = 0.6262415;
      }
      break;
    case 39:
      if (refsys == 1) {
        alphaiw[0] = 156.5816530;
        alphaiw[1] = 133.4096570;
        alphaiw[2] = 94.6810390;
        alphaiw[3] = 47.1887390;
        alphaiw[4] = 27.8621270;
        alphaiw[5] = 18.8348770;
        alphaiw[6] = 13.9092620;
        alphaiw[7] = 10.8900460;
        alphaiw[8] = 8.8730220;
        alphaiw[9] = 7.4352670;
        alphaiw[10] = 6.3585500;
        alphaiw[11] = 5.5211390;
        alphaiw[12] = 4.3014730;
        alphaiw[13] = 3.4566600;
        alphaiw[14] = 2.8399260;
        alphaiw[15] = 2.3734080;
        alphaiw[16] = 2.0112610;
        alphaiw[17] = 1.3945770;
        alphaiw[18] = 1.0190510;
        alphaiw[19] = 0.6075350;
        alphaiw[20] = 0.4008910;
        alphaiw[21] = 0.1843540;
        alphaiw[22] = 0.1051420;
      }
      else if (refsys == 2) {
        alphaiw[0] = 126.8052320;
        alphaiw[1] = 88.7747170;
        alphaiw[2] = 64.1464690;
        alphaiw[3] = 40.2998810;
        alphaiw[4] = 28.2499780;
        alphaiw[5] = 21.0643660;
        alphaiw[6] = 16.4186290;
        alphaiw[7] = 13.2334170;
        alphaiw[8] = 10.9441490;
        alphaiw[9] = 9.2339520;
        alphaiw[10] = 7.9151340;
        alphaiw[11] = 6.8713270;
        alphaiw[12] = 5.3330370;
        alphaiw[13] = 4.2636390;
        alphaiw[14] = 3.4855460;
        alphaiw[15] = 2.9002550;
        alphaiw[16] = 2.4485960;
        alphaiw[17] = 1.6861640;
        alphaiw[18] = 1.2266060;
        alphaiw[19] = 0.7275190;
        alphaiw[20] = 0.4787640;
        alphaiw[21] = 0.2194880;
        alphaiw[22] = 0.1250190;
      }
      else if (refsys == 3) {
        alphaiw[0] = 51.9455750;
        alphaiw[1] = 50.4162010;
        alphaiw[2] = 46.5373240;
        alphaiw[3] = 36.6790590;
        alphaiw[4] = 28.1986330;
        alphaiw[5] = 21.9785510;
        alphaiw[6] = 17.5295770;
        alphaiw[7] = 14.3039920;
        alphaiw[8] = 11.9077650;
        alphaiw[9] = 10.0810740;
        alphaiw[10] = 8.6547010;
        alphaiw[11] = 7.5170970;
        alphaiw[12] = 5.8310100;
        alphaiw[13] = 4.6555680;
        alphaiw[14] = 3.8004240;
        alphaiw[15] = 3.1579800;
        alphaiw[16] = 2.6630160;
        alphaiw[17] = 1.8296590;
        alphaiw[18] = 1.3289880;
        alphaiw[19] = 0.7868860;
        alphaiw[20] = 0.5173750;
        alphaiw[21] = 0.2369610;
        alphaiw[22] = 0.1349160;
      }
      break;
    case 40:
      if (refsys == 1) {
        alphaiw[0] = 129.1649250;
        alphaiw[1] = 112.9115260;
        alphaiw[2] = 84.2744150;
        alphaiw[3] = 45.2123570;
        alphaiw[4] = 27.6899320;
        alphaiw[5] = 19.0856380;
        alphaiw[6] = 14.2390790;
        alphaiw[7] = 11.2024000;
        alphaiw[8] = 9.1447760;
        alphaiw[9] = 7.6666660;
        alphaiw[10] = 6.5563690;
        alphaiw[11] = 5.6928180;
        alphaiw[12] = 4.4382110;
        alphaiw[13] = 3.5726300;
        alphaiw[14] = 2.9422010;
        alphaiw[15] = 2.4654780;
        alphaiw[16] = 2.0949780;
        alphaiw[17] = 1.4616900;
        alphaiw[18] = 1.0734540;
        alphaiw[19] = 0.6444360;
        alphaiw[20] = 0.4270850;
        alphaiw[21] = 0.1974770;
        alphaiw[22] = 0.1128890;
      }
      else if (refsys == 2) {
        alphaiw[0] = 59.5329360;
        alphaiw[1] = 59.0731440;
        alphaiw[2] = 49.3004260;
        alphaiw[3] = 35.6901590;
        alphaiw[4] = 26.5930640;
        alphaiw[5] = 20.4777520;
        alphaiw[6] = 16.2621680;
        alphaiw[7] = 13.2552040;
        alphaiw[8] = 11.0380980;
        alphaiw[9] = 9.3543390;
        alphaiw[10] = 8.0423850;
        alphaiw[11] = 6.9973630;
        alphaiw[12] = 5.4495930;
        alphaiw[13] = 4.3698600;
        alphaiw[14] = 3.5825450;
        alphaiw[15] = 2.9890610;
        alphaiw[16] = 2.5300090;
        alphaiw[17] = 1.7518420;
        alphaiw[18] = 1.2798730;
        alphaiw[19] = 0.7635970;
        alphaiw[20] = 0.5043350;
        alphaiw[21] = 0.2322730;
        alphaiw[22] = 0.1325630;
      }
      else if (refsys == 3) {
        alphaiw[0] = 53.7033680;
        alphaiw[1] = 52.3694770;
        alphaiw[2] = 48.8771380;
        alphaiw[3] = 39.4788990;
        alphaiw[4] = 30.9424830;
        alphaiw[5] = 24.4286920;
        alphaiw[6] = 19.6266710;
        alphaiw[7] = 16.0687640;
        alphaiw[8] = 13.3881100;
        alphaiw[9] = 11.3281120;
        alphaiw[10] = 9.7137520;
        alphaiw[11] = 8.4253900;
        alphaiw[12] = 6.5198670;
        alphaiw[13] = 5.1977820;
        alphaiw[14] = 4.2403830;
        alphaiw[15] = 3.5235900;
        alphaiw[16] = 2.9725540;
        alphaiw[17] = 2.0460840;
        alphaiw[18] = 1.4892320;
        alphaiw[19] = 0.8847240;
        alphaiw[20] = 0.5830200;
        alphaiw[21] = 0.2678100;
        alphaiw[22] = 0.1526730;
      }
      else if (refsys == 4) {
        alphaiw[0] = 33.6465243;
        alphaiw[1] = 33.4038388;
        alphaiw[2] = 32.7056564;
        alphaiw[3] = 30.2880677;
        alphaiw[4] = 27.2060499;
        alphaiw[5] = 24.0984692;
        alphaiw[6] = 21.2765610;
        alphaiw[7] = 18.8296599;
        alphaiw[8] = 16.7447048;
        alphaiw[9] = 14.9743210;
        alphaiw[10] = 13.4663529;
        alphaiw[11] = 12.1742122;
        alphaiw[12] = 10.0904600;
        alphaiw[13] = 8.4969686;
        alphaiw[14] = 7.2491379;
        alphaiw[15] = 6.2526961;
        alphaiw[16] = 5.4441107;
        alphaiw[17] = 3.9836870;
        alphaiw[18] = 3.0306974;
        alphaiw[19] = 1.9121250;
        alphaiw[20] = 1.3109813;
        alphaiw[21] = 0.6388831;
        alphaiw[22] = 0.3766424;
      }
      break;
    case 41:
      if (refsys == 1) {
        alphaiw[0] = 112.5868220;
        alphaiw[1] = 98.9946210;
        alphaiw[2] = 75.5219100;
        alphaiw[3] = 42.5430350;
        alphaiw[4] = 26.8101620;
        alphaiw[5] = 18.7904950;
        alphaiw[6] = 14.1713590;
        alphaiw[7] = 11.2319370;
        alphaiw[8] = 9.2165940;
        alphaiw[9] = 7.7556700;
        alphaiw[10] = 6.6508620;
        alphaiw[11] = 5.7874800;
        alphaiw[12] = 4.5279270;
        alphaiw[13] = 3.6564020;
        alphaiw[14] = 3.0206270;
        alphaiw[15] = 2.5390860;
        alphaiw[16] = 2.1640690;
        alphaiw[17] = 1.5202440;
        alphaiw[18] = 1.1226240;
        alphaiw[19] = 0.6790940;
        alphaiw[20] = 0.4521310;
        alphaiw[21] = 0.2102320;
        alphaiw[22] = 0.1204560;
      }
      else if (refsys == 2) {
        alphaiw[0] = 55.9986780;
        alphaiw[1] = 60.8749560;
        alphaiw[2] = 50.6953400;
        alphaiw[3] = 35.6643390;
        alphaiw[4] = 26.4531740;
        alphaiw[5] = 20.4359360;
        alphaiw[6] = 16.2988900;
        alphaiw[7] = 13.3351670;
        alphaiw[8] = 11.1377280;
        alphaiw[9] = 9.4604950;
        alphaiw[10] = 8.1484630;
        alphaiw[11] = 7.1004180;
        alphaiw[12] = 5.5443520;
        alphaiw[13] = 4.4569430;
        alphaiw[14] = 3.6632070;
        alphaiw[15] = 3.0641820;
        alphaiw[16] = 2.6001160;
        alphaiw[17] = 1.8107380;
        alphaiw[18] = 1.3291220;
        alphaiw[19] = 0.7982160;
        alphaiw[20] = 0.5293550;
        alphaiw[21] = 0.2450550;
        alphaiw[22] = 0.1401740;
      }
      else if (refsys == 3) {
        alphaiw[0] = 51.3329100;
        alphaiw[1] = 49.0284950;
        alphaiw[2] = 44.6433950;
        alphaiw[3] = 35.3677840;
        alphaiw[4] = 27.6830280;
        alphaiw[5] = 21.9413960;
        alphaiw[6] = 17.7220800;
        alphaiw[7] = 14.5880150;
        alphaiw[8] = 12.2159700;
        alphaiw[9] = 10.3838810;
        alphaiw[10] = 8.9410440;
        alphaiw[11] = 7.7843510;
        alphaiw[12] = 6.0637690;
        alphaiw[13] = 4.8619510;
        alphaiw[14] = 3.9865780;
        alphaiw[15] = 3.3277110;
        alphaiw[16] = 2.8186730;
        alphaiw[17] = 1.9562960;
        alphaiw[18] = 1.4327150;
        alphaiw[19] = 0.8581570;
        alphaiw[20] = 0.5683180;
        alphaiw[21] = 0.2626950;
        alphaiw[22] = 0.1501750;
      }
      else if (refsys == 4) {
        alphaiw[0] = 39.2699458;
        alphaiw[1] = 38.9559085;
        alphaiw[2] = 38.0619736;
        alphaiw[3] = 35.0545164;
        alphaiw[4] = 31.3566941;
        alphaiw[5] = 27.7212247;
        alphaiw[6] = 24.4670081;
        alphaiw[7] = 21.6658391;
        alphaiw[8] = 19.2870190;
        alphaiw[9] = 17.2695766;
        alphaiw[10] = 15.5513454;
        alphaiw[11] = 14.0783331;
        alphaiw[12] = 11.7000325;
        alphaiw[13] = 9.8776054;
        alphaiw[14] = 8.4472870;
        alphaiw[15] = 7.3024041;
        alphaiw[16] = 6.3710983;
        alphaiw[17] = 4.6822974;
        alphaiw[18] = 3.5742534;
        alphaiw[19] = 2.2657569;
        alphaiw[20] = 1.5581961;
        alphaiw[21] = 0.7625966;
        alphaiw[22] = 0.4505494;
      }
      break;
    case 42:
      if (refsys == 1) {
        alphaiw[0] = 80.1166370;
        alphaiw[1] = 72.0033490;
        alphaiw[2] = 56.3464510;
        alphaiw[3] = 33.3622730;
        alphaiw[4] = 22.2521730;
        alphaiw[5] = 16.3666450;
        alphaiw[6] = 12.8045350;
        alphaiw[7] = 10.4296810;
        alphaiw[8] = 8.7352680;
        alphaiw[9] = 7.4656260;
        alphaiw[10] = 6.4789340;
        alphaiw[11] = 5.6904860;
        alphaiw[12] = 4.5115160;
        alphaiw[13] = 3.6758650;
        alphaiw[14] = 3.0566210;
        alphaiw[15] = 2.5826100;
        alphaiw[16] = 2.2106410;
        alphaiw[17] = 1.5662560;
        alphaiw[18] = 1.1641470;
        alphaiw[19] = 0.7106370;
        alphaiw[20] = 0.4758290;
        alphaiw[21] = 0.2228080;
        alphaiw[22] = 0.1280230;
      }
      else if (refsys == 2) {
        alphaiw[0] = 57.3696220;
        alphaiw[1] = 53.6619600;
        alphaiw[2] = 46.4519190;
        alphaiw[3] = 33.7721100;
        alphaiw[4] = 25.3796990;
        alphaiw[5] = 19.8008640;
        alphaiw[6] = 15.9255030;
        alphaiw[7] = 13.1223040;
        alphaiw[8] = 11.0249000;
        alphaiw[9] = 9.4107880;
        alphaiw[10] = 8.1389620;
        alphaiw[11] = 7.1166930;
        alphaiw[12] = 5.5878340;
        alphaiw[13] = 4.5113880;
        alphaiw[14] = 3.7214650;
        alphaiw[15] = 3.1229210;
        alphaiw[16] = 2.6576930;
        alphaiw[17] = 1.8625700;
        alphaiw[18] = 1.3742440;
        alphaiw[19] = 0.8315500;
        alphaiw[20] = 0.5541640;
        alphaiw[21] = 0.2581730;
        alphaiw[22] = 0.1480840;
      }
      else if (refsys == 3) {
        alphaiw[0] = 67.2624950;
        alphaiw[1] = 66.0303950;
        alphaiw[2] = 62.7575620;
        alphaiw[3] = 53.5269770;
        alphaiw[4] = 44.3751560;
        alphaiw[5] = 36.6776440;
        alphaiw[6] = 30.4829790;
        alphaiw[7] = 25.5527440;
        alphaiw[8] = 21.6261800;
        alphaiw[9] = 18.4804960;
        alphaiw[10] = 15.9392930;
        alphaiw[11] = 13.8670690;
        alphaiw[12] = 10.7428020;
        alphaiw[13] = 8.5472490;
        alphaiw[14] = 6.9523710;
        alphaiw[15] = 5.7601520;
        alphaiw[16] = 4.8468120;
        alphaiw[17] = 3.3211190;
        alphaiw[18] = 2.4114970;
        alphaiw[19] = 1.4303190;
        alphaiw[20] = 0.9425060;
        alphaiw[21] = 0.4333420;
        alphaiw[22] = 0.2472210;
      }
      else if (refsys == 4) {
        alphaiw[0] = 57.1122440;
        alphaiw[1] = 54.5507140;
        alphaiw[2] = 49.8607540;
        alphaiw[3] = 40.5094150;
        alphaiw[4] = 32.7396620;
        alphaiw[5] = 26.6350470;
        alphaiw[6] = 21.9139730;
        alphaiw[7] = 18.2593580;
        alphaiw[8] = 15.4061580;
        alphaiw[9] = 13.1523240;
        alphaiw[10] = 11.3489830;
        alphaiw[11] = 9.8875440;
        alphaiw[12] = 7.6939570;
        alphaiw[13] = 6.1540760;
        alphaiw[14] = 5.0326460;
        alphaiw[15] = 4.1907210;
        alphaiw[16] = 3.5424610;
        alphaiw[17] = 2.4502640;
        alphaiw[18] = 1.7914150;
        alphaiw[19] = 1.0718980;
        alphaiw[20] = 0.7099120;
        alphaiw[21] = 0.3284040;
        alphaiw[22] = 0.1878520;
      }
      break;
    case 43:
      if (refsys == 1) {
        alphaiw[0] = 92.9387160;
        alphaiw[1] = 81.9452450;
        alphaiw[2] = 63.3698580;
        alphaiw[3] = 37.3296730;
        alphaiw[4] = 24.3564850;
        alphaiw[5] = 17.4779390;
        alphaiw[6] = 13.4110250;
        alphaiw[7] = 10.7755340;
        alphaiw[8] = 8.9429790;
        alphaiw[9] = 7.5986600;
        alphaiw[10] = 6.5712030;
        alphaiw[11] = 5.7605070;
        alphaiw[12] = 4.5631700;
        alphaiw[13] = 3.7227110;
        alphaiw[14] = 3.1024520;
        alphaiw[15] = 2.6281370;
        alphaiw[16] = 2.2556760;
        alphaiw[17] = 1.6085080;
        alphaiw[18] = 1.2022890;
        alphaiw[19] = 0.7402780;
        alphaiw[20] = 0.4986150;
        alphaiw[21] = 0.2353880;
        alphaiw[22] = 0.1357590;
      }
      else if (refsys == 2) {
        alphaiw[0] = 75.8730320;
        alphaiw[1] = 65.9104100;
        alphaiw[2] = 51.9040420;
        alphaiw[3] = 34.7083320;
        alphaiw[4] = 25.4784390;
        alphaiw[5] = 19.7433960;
        alphaiw[6] = 15.8501840;
        alphaiw[7] = 13.0603640;
        alphaiw[8] = 10.9826550;
        alphaiw[9] = 9.3875760;
        alphaiw[10] = 8.1320480;
        alphaiw[11] = 7.1229420;
        alphaiw[12] = 5.6121410;
        alphaiw[13] = 4.5457470;
        alphaiw[14] = 3.7609510;
        alphaiw[15] = 3.1646320;
        alphaiw[16] = 2.6999290;
        alphaiw[17] = 1.9026050;
        alphaiw[18] = 1.4102760;
        alphaiw[19] = 0.8594030;
        alphaiw[20] = 0.5755190;
        alphaiw[21] = 0.2699220;
        alphaiw[22] = 0.1552940;
      }
      else if (refsys == 3) {
        alphaiw[0] = 58.5751980;
        alphaiw[1] = 57.7047570;
        alphaiw[2] = 55.3283600;
        alphaiw[3] = 48.1684600;
        alphaiw[4] = 40.5340630;
        alphaiw[5] = 33.8274610;
        alphaiw[6] = 28.3028200;
        alphaiw[7] = 23.8467240;
        alphaiw[8] = 20.2669430;
        alphaiw[9] = 17.3809560;
        alphaiw[10] = 15.0378770;
        alphaiw[11] = 13.1192270;
        alphaiw[12] = 10.2118960;
        alphaiw[13] = 8.1570940;
        alphaiw[14] = 6.6575030;
        alphaiw[15] = 5.5321270;
        alphaiw[16] = 4.6670770;
        alphaiw[17] = 3.2152040;
        alphaiw[18] = 2.3443700;
        alphaiw[19] = 1.3987970;
        alphaiw[20] = 0.9253150;
        alphaiw[21] = 0.4276690;
        alphaiw[22] = 0.2445830;
      }
      break;
    case 44:
      if (refsys == 1) {
        alphaiw[0] = 64.2873650;
        alphaiw[1] = 58.2984130;
        alphaiw[2] = 46.6241190;
        alphaiw[3] = 28.9081460;
        alphaiw[4] = 19.8860070;
        alphaiw[5] = 14.9414970;
        alphaiw[6] = 11.8820490;
        alphaiw[7] = 9.8067370;
        alphaiw[8] = 8.3050910;
        alphaiw[9] = 7.1669030;
        alphaiw[10] = 6.2736860;
        alphaiw[11] = 5.5536060;
        alphaiw[12] = 4.4638560;
        alphaiw[13] = 3.6792510;
        alphaiw[14] = 3.0894480;
        alphaiw[15] = 2.6322260;
        alphaiw[16] = 2.2694820;
        alphaiw[17] = 1.6317090;
        alphaiw[18] = 1.2266730;
        alphaiw[19] = 0.7613680;
        alphaiw[20] = 0.5156110;
        alphaiw[21] = 0.2453210;
        alphaiw[22] = 0.1420300;
      }
      else if (refsys == 2) {
        alphaiw[0] = 45.7015170;
        alphaiw[1] = 44.1640550;
        alphaiw[2] = 40.0828620;
        alphaiw[3] = 30.9952860;
        alphaiw[4] = 23.8618350;
        alphaiw[5] = 18.8259100;
        alphaiw[6] = 15.2521470;
        alphaiw[7] = 12.6439180;
        alphaiw[8] = 10.6826650;
        alphaiw[9] = 9.1676420;
        alphaiw[10] = 7.9697090;
        alphaiw[11] = 7.0033880;
        alphaiw[12] = 5.5501890;
        alphaiw[13] = 4.5187320;
        alphaiw[14] = 3.7556240;
        alphaiw[15] = 3.1728410;
        alphaiw[16] = 2.7165260;
        alphaiw[17] = 1.9281300;
        alphaiw[18] = 1.4370560;
        alphaiw[19] = 0.8825770;
        alphaiw[20] = 0.5941370;
        alphaiw[21] = 0.2807420;
        alphaiw[22] = 0.1621180;
      }
      else if (refsys == 3) {
        alphaiw[0] = 53.6449610;
        alphaiw[1] = 52.9021650;
        alphaiw[2] = 50.8530030;
        alphaiw[3] = 44.5291150;
        alphaiw[4] = 37.6245730;
        alphaiw[5] = 31.4884170;
        alphaiw[6] = 26.4097500;
        alphaiw[7] = 22.3040630;
        alphaiw[8] = 19.0002040;
        alphaiw[9] = 16.3321480;
        alphaiw[10] = 14.1620830;
        alphaiw[11] = 12.3817320;
        alphaiw[12] = 9.6763600;
        alphaiw[13] = 7.7570110;
        alphaiw[14] = 6.3512380;
        alphaiw[15] = 5.2927530;
        alphaiw[16] = 4.4766230;
        alphaiw[17] = 3.1006590;
        alphaiw[18] = 2.2704610;
        alphaiw[19] = 1.3630470;
        alphaiw[20] = 0.9053930;
        alphaiw[21] = 0.4209110;
        alphaiw[22] = 0.2414160;
      }
      break;
    case 45:
      if (refsys == 1) {
        alphaiw[0] = 68.1542450;
        alphaiw[1] = 62.9960050;
        alphaiw[2] = 52.0807400;
        alphaiw[3] = 32.9144330;
        alphaiw[4] = 22.1117060;
        alphaiw[5] = 16.1257920;
        alphaiw[6] = 12.5115460;
        alphaiw[7] = 10.1406800;
        alphaiw[8] = 8.4788070;
        alphaiw[9] = 7.2525150;
        alphaiw[10] = 6.3108490;
        alphaiw[11] = 5.5647590;
        alphaiw[12] = 4.4563970;
        alphaiw[13] = 3.6717400;
        alphaiw[14] = 3.0873140;
        alphaiw[15] = 2.6361620;
        alphaiw[16] = 2.2786150;
        alphaiw[17] = 1.6486870;
        alphaiw[18] = 1.2462650;
        alphaiw[19] = 0.7799260;
        alphaiw[20] = 0.5311430;
        alphaiw[21] = 0.2547010;
        alphaiw[22] = 0.1480120;
      }
      else if (refsys == 2) {
        alphaiw[0] = 43.3417530;
        alphaiw[1] = 41.4841810;
        alphaiw[2] = 37.5198890;
        alphaiw[3] = 29.0651560;
        alphaiw[4] = 22.5178410;
        alphaiw[5] = 17.8848750;
        alphaiw[6] = 14.5797210;
        alphaiw[7] = 12.1532300;
        alphaiw[8] = 10.3179030;
        alphaiw[9] = 8.8923860;
        alphaiw[10] = 7.7596700;
        alphaiw[11] = 6.8419520;
        alphaiw[12] = 5.4542350;
        alphaiw[13] = 4.4628030;
        alphaiw[14] = 3.7251260;
        alphaiw[15] = 3.1588890;
        alphaiw[16] = 2.7134840;
        alphaiw[17] = 1.9388590;
        alphaiw[18] = 1.4524300;
        alphaiw[19] = 0.8985520;
        alphaiw[20] = 0.6078940;
        alphaiw[21] = 0.2892870;
        alphaiw[22] = 0.1676290;
      }
      else if (refsys == 3) {
        alphaiw[0] = 37.7442490;
        alphaiw[1] = 36.6313080;
        alphaiw[2] = 34.0210470;
        alphaiw[3] = 27.9047180;
        alphaiw[4] = 22.6219420;
        alphaiw[5] = 18.5364410;
        alphaiw[6] = 15.4279090;
        alphaiw[7] = 13.0366770;
        alphaiw[8] = 11.1653460;
        alphaiw[9] = 9.6752930;
        alphaiw[10] = 8.4697300;
        alphaiw[11] = 7.4802140;
        alphaiw[12] = 5.9656450;
        alphaiw[13] = 4.8738020;
        alphaiw[14] = 4.0587720;
        alphaiw[15] = 3.4331020;
        alphaiw[16] = 2.9417130;
        alphaiw[17] = 2.0906290;
        alphaiw[18] = 1.5597110;
        alphaiw[19] = 0.9597480;
        alphaiw[20] = 0.6472660;
        alphaiw[21] = 0.3069220;
        alphaiw[22] = 0.1776240;
      }
      else if (refsys == 4) {
        alphaiw[0] = 50.7362600;
        alphaiw[1] = 49.4436100;
        alphaiw[2] = 46.2123700;
        alphaiw[3] = 37.7336800;
        alphaiw[4] = 30.0265600;
        alphaiw[5] = 24.1196600;
        alphaiw[6] = 19.7259100;
        alphaiw[7] = 16.4228700;
        alphaiw[8] = 13.8898500;
        alphaiw[9] = 11.9075300;
        alphaiw[10] = 10.3272700;
        alphaiw[11] = 9.0468100;
        alphaiw[12] = 7.1181400;
        alphaiw[13] = 5.7535600;
        alphaiw[14] = 4.7506500;
        alphaiw[15] = 3.9906900;
        alphaiw[16] = 3.4003100;
        alphaiw[17] = 2.3920200;
        alphaiw[18] = 1.7726100;
        alphaiw[19] = 1.0819100;
        alphaiw[20] = 0.7262700;
        alphaiw[21] = 0.3424400;
        alphaiw[22] = 0.1977000;
      }
      break;
    case 46:
      if (refsys == 1) {
        alphaiw[0] = 60.8864630;
        alphaiw[1] = 56.8478940;
        alphaiw[2] = 47.8540600;
        alphaiw[3] = 30.9841900;
        alphaiw[4] = 21.0472500;
        alphaiw[5] = 15.4485530;
        alphaiw[6] = 12.0425860;
        alphaiw[7] = 9.7984690;
        alphaiw[8] = 8.2204320;
        alphaiw[9] = 7.0529890;
        alphaiw[10] = 6.1545310;
        alphaiw[11] = 5.4412960;
        alphaiw[12] = 4.3791150;
        alphaiw[13] = 3.6247120;
        alphaiw[14] = 3.0609020;
        alphaiw[15] = 2.6240310;
        alphaiw[16] = 2.2764150;
        alphaiw[17] = 1.6598380;
        alphaiw[18] = 1.2622490;
        alphaiw[19] = 0.7967230;
        alphaiw[20] = 0.5457380;
        alphaiw[21] = 0.2639640;
        alphaiw[22] = 0.1541040;
      }
      else if (refsys == 2) {
        alphaiw[0] = 40.2110290;
        alphaiw[1] = 38.5072500;
        alphaiw[2] = 34.8932240;
        alphaiw[3] = 27.2138020;
        alphaiw[4] = 21.2371640;
        alphaiw[5] = 16.9805470;
        alphaiw[6] = 13.9260910;
        alphaiw[7] = 11.6707350;
        alphaiw[8] = 9.9551980;
        alphaiw[9] = 8.6156180;
        alphaiw[10] = 7.5460330;
        alphaiw[11] = 6.6757400;
        alphaiw[12] = 5.3527110;
        alphaiw[13] = 4.4016410;
        alphaiw[14] = 3.6902360;
        alphaiw[15] = 3.1415360;
        alphaiw[16] = 2.7079930;
        alphaiw[17] = 1.9489220;
        alphaiw[18] = 1.4681030;
        alphaiw[19] = 0.9154940;
        alphaiw[20] = 0.6227080;
        alphaiw[21] = 0.2987130;
        alphaiw[22] = 0.1738250;
      }
      else if (refsys == 3) {
        alphaiw[0] = 52.9449710;
        alphaiw[1] = 51.8101820;
        alphaiw[2] = 48.8650090;
        alphaiw[3] = 40.9527250;
        alphaiw[4] = 33.5198470;
        alphaiw[5] = 27.5254650;
        alphaiw[6] = 22.8486080;
        alphaiw[7] = 19.2026620;
        alphaiw[8] = 16.3345560;
        alphaiw[9] = 14.0509400;
        alphaiw[10] = 12.2094700;
        alphaiw[11] = 10.7060500;
        alphaiw[12] = 8.4272930;
        alphaiw[13] = 6.8085780;
        alphaiw[14] = 5.6176560;
        alphaiw[15] = 4.7154850;
        alphaiw[16] = 4.0151810;
        alphaiw[17] = 2.8211020;
        alphaiw[18] = 2.0891520;
        alphaiw[19] = 1.2745700;
        alphaiw[20] = 0.8557330;
        alphaiw[21] = 0.4038840;
        alphaiw[22] = 0.2333960;
      }
      break;
    case 47:
      if (refsys == 1) {
        alphaiw[0] = 54.1989880;
        alphaiw[1] = 49.0489660;
        alphaiw[2] = 38.9829400;
        alphaiw[3] = 23.9821320;
        alphaiw[4] = 16.6050290;
        alphaiw[5] = 12.6373760;
        alphaiw[6] = 10.1995360;
        alphaiw[7] = 8.5437370;
        alphaiw[8] = 7.3366070;
        alphaiw[9] = 6.4116550;
        alphaiw[10] = 5.6772120;
        alphaiw[11] = 5.0784690;
        alphaiw[12] = 4.1592690;
        alphaiw[13] = 3.4860530;
        alphaiw[14] = 2.9720170;
        alphaiw[15] = 2.5673730;
        alphaiw[16] = 2.2414010;
        alphaiw[17] = 1.6540770;
        alphaiw[18] = 1.2685970;
        alphaiw[19] = 0.8095530;
        alphaiw[20] = 0.5581500;
        alphaiw[21] = 0.2721620;
        alphaiw[22] = 0.1594640;
      }
      else if (refsys == 2) {
        alphaiw[0] = 36.6000980;
        alphaiw[1] = 35.1208740;
        alphaiw[2] = 31.6461790;
        alphaiw[3] = 24.0228290;
        alphaiw[4] = 18.3425350;
        alphaiw[5] = 14.4988310;
        alphaiw[6] = 11.8463360;
        alphaiw[7] = 9.9378870;
        alphaiw[8] = 8.5081750;
        alphaiw[9] = 7.4002440;
        alphaiw[10] = 6.5179040;
        alphaiw[11] = 5.7995690;
        alphaiw[12] = 4.7030760;
        alphaiw[13] = 3.9082530;
        alphaiw[14] = 3.3079420;
        alphaiw[15] = 2.8402610;
        alphaiw[16] = 2.4670520;
        alphaiw[17] = 1.8034180;
        alphaiw[18] = 1.3743830;
        alphaiw[19] = 0.8703890;
        alphaiw[20] = 0.5975480;
        alphaiw[21] = 0.2899250;
        alphaiw[22] = 0.1695250;
      }
      break;
    case 48:
      if (refsys == 1) {
        alphaiw[0] = 51.2430300;
        alphaiw[1] = 48.3502990;
        alphaiw[2] = 41.5746300;
        alphaiw[3] = 27.7890000;
        alphaiw[4] = 19.1321820;
        alphaiw[5] = 14.1382230;
        alphaiw[6] = 11.0784560;
        alphaiw[7] = 9.0587800;
        alphaiw[8] = 7.6378310;
        alphaiw[9] = 6.5859030;
        alphaiw[10] = 5.7753370;
        alphaiw[11] = 5.1307120;
        alphaiw[12] = 4.1675160;
        alphaiw[13] = 3.4798590;
        alphaiw[14] = 2.9630260;
        alphaiw[15] = 2.5600920;
        alphaiw[16] = 2.2373290;
        alphaiw[17] = 1.6578750;
        alphaiw[18] = 1.2772960;
        alphaiw[19] = 0.8215180;
        alphaiw[20] = 0.5698100;
        alphaiw[21] = 0.2806600;
        alphaiw[22] = 0.1654590;
      }
      else if (refsys == 2) {
        alphaiw[0] = 35.9140330;
        alphaiw[1] = 35.1566180;
        alphaiw[2] = 33.1068480;
        alphaiw[3] = 27.1969530;
        alphaiw[4] = 21.5145310;
        alphaiw[5] = 17.1191360;
        alphaiw[6] = 13.8992970;
        alphaiw[7] = 11.5355770;
        alphaiw[8] = 9.7648350;
        alphaiw[9] = 8.4054680;
        alphaiw[10] = 7.3368940;
        alphaiw[11] = 6.4788070;
        alphaiw[12] = 5.1930790;
        alphaiw[13] = 4.2809850;
        alphaiw[14] = 3.6036050;
        alphaiw[15] = 3.0826690;
        alphaiw[16] = 2.6710330;
        alphaiw[17] = 1.9469920;
        alphaiw[18] = 1.4831570;
        alphaiw[19] = 0.9406720;
        alphaiw[20] = 0.6472420;
        alphaiw[21] = 0.3156970;
        alphaiw[22] = 0.1853400;
      }
      break;
    case 49:
      if (refsys == 1) {
        alphaiw[0] = 72.9150960;
        alphaiw[1] = 66.0753600;
        alphaiw[2] = 53.3578280;
        alphaiw[3] = 33.5644110;
        alphaiw[4] = 22.4687340;
        alphaiw[5] = 16.1472430;
        alphaiw[6] = 12.3105260;
        alphaiw[7] = 9.8245790;
        alphaiw[8] = 8.1173940;
        alphaiw[9] = 6.8860560;
        alphaiw[10] = 5.9610840;
        alphaiw[11] = 5.2425400;
        alphaiw[12] = 4.1995990;
        alphaiw[13] = 3.4778910;
        alphaiw[14] = 2.9472540;
        alphaiw[15] = 2.5398030;
        alphaiw[16] = 2.2168390;
        alphaiw[17] = 1.6430710;
        alphaiw[18] = 1.2690370;
        alphaiw[19] = 0.8214960;
        alphaiw[20] = 0.5729700;
        alphaiw[21] = 0.2844110;
        alphaiw[22] = 0.1681820;
      }
      else if (refsys == 2) {
        alphaiw[0] = 68.5448320;
        alphaiw[1] = 63.3499550;
        alphaiw[2] = 52.9158640;
        alphaiw[3] = 35.1620400;
        alphaiw[4] = 24.4638730;
        alphaiw[5] = 18.0151280;
        alphaiw[6] = 13.9204010;
        alphaiw[7] = 11.1780650;
        alphaiw[8] = 9.2513440;
        alphaiw[9] = 7.8406570;
        alphaiw[10] = 6.7710770;
        alphaiw[11] = 5.9359020;
        alphaiw[12] = 4.7209280;
        alphaiw[13] = 3.8821540;
        alphaiw[14] = 3.2690480;
        alphaiw[15] = 2.8016450;
        alphaiw[16] = 2.4338840;
        alphaiw[17] = 1.7877840;
        alphaiw[18] = 1.3722790;
        alphaiw[19] = 0.8814810;
        alphaiw[20] = 0.6120400;
        alphaiw[21] = 0.3021170;
        alphaiw[22] = 0.1782190;
      }
      else if (refsys == 3) {
        alphaiw[0] = 50.7622140;
        alphaiw[1] = 48.0385490;
        alphaiw[2] = 43.0567370;
        alphaiw[3] = 32.8278830;
        alphaiw[4] = 24.7634540;
        alphaiw[5] = 19.0651230;
        alphaiw[6] = 15.0969960;
        alphaiw[7] = 12.2851460;
        alphaiw[8] = 10.2375520;
        alphaiw[9] = 8.7031490;
        alphaiw[10] = 7.5220080;
        alphaiw[11] = 6.5906520;
        alphaiw[12] = 5.2252970;
        alphaiw[13] = 4.2790840;
        alphaiw[14] = 3.5880180;
        alphaiw[15] = 3.0628180;
        alphaiw[16] = 2.6512680;
        alphaiw[17] = 1.9334490;
        alphaiw[18] = 1.4763250;
        alphaiw[19] = 0.9418740;
        alphaiw[20] = 0.6513060;
        alphaiw[21] = 0.3198580;
        alphaiw[22] = 0.1882630;
      }
      else if (refsys == 4) {
        alphaiw[0] = 41.1871070;
        alphaiw[1] = 40.4130170;
        alphaiw[2] = 38.2892500;
        alphaiw[3] = 31.9210850;
        alphaiw[4] = 25.4481130;
        alphaiw[5] = 20.2225040;
        alphaiw[6] = 16.3010840;
        alphaiw[7] = 13.3954500;
        alphaiw[8] = 11.2196820;
        alphaiw[9] = 9.5597380;
        alphaiw[10] = 8.2670060;
        alphaiw[11] = 7.2399920;
        alphaiw[12] = 5.7255130;
        alphaiw[13] = 4.6728770;
        alphaiw[14] = 3.9045910;
        alphaiw[15] = 3.3221260;
        alphaiw[16] = 2.8671790;
        alphaiw[17] = 2.0782430;
        alphaiw[18] = 1.5798240;
        alphaiw[19] = 1.0020180;
        alphaiw[20] = 0.6904450;
        alphaiw[21] = 0.3375630;
        alphaiw[22] = 0.1982910;
      }
      break;
    case 50:
      if (refsys == 1) {
        alphaiw[0] = 57.3040270;
        alphaiw[1] = 54.4476260;
        alphaiw[2] = 47.7957070;
        alphaiw[3] = 33.6841460;
        alphaiw[4] = 23.7376670;
        alphaiw[5] = 17.3960520;
        alphaiw[6] = 13.3024310;
        alphaiw[7] = 10.5620400;
        alphaiw[8] = 8.6525610;
        alphaiw[9] = 7.2708170;
        alphaiw[10] = 6.2366020;
        alphaiw[11] = 5.4393040;
        alphaiw[12] = 4.2988630;
        alphaiw[13] = 3.5266460;
        alphaiw[14] = 2.9701000;
        alphaiw[15] = 2.5497820;
        alphaiw[16] = 2.2208960;
        alphaiw[17] = 1.6446800;
        alphaiw[18] = 1.2728700;
        alphaiw[19] = 0.8288590;
        alphaiw[20] = 0.5812410;
        alphaiw[21] = 0.2912600;
        alphaiw[22] = 0.1732200;
      }
      else if (refsys == 2) {
        alphaiw[0] = 53.7298320;
        alphaiw[1] = 51.4329790;
        alphaiw[2] = 46.0716700;
        alphaiw[3] = 34.2024090;
        alphaiw[4] = 25.1302790;
        alphaiw[5] = 18.9316540;
        alphaiw[6] = 14.7246270;
        alphaiw[7] = 11.8066270;
        alphaiw[8] = 9.7217990;
        alphaiw[9] = 8.1863100;
        alphaiw[10] = 7.0228060;
        alphaiw[11] = 6.1183050;
        alphaiw[12] = 4.8155710;
        alphaiw[13] = 3.9303330;
        alphaiw[14] = 3.2929680;
        alphaiw[15] = 2.8132790;
        alphaiw[16] = 2.4396930;
        alphaiw[17] = 1.7907330;
        alphaiw[18] = 1.3769680;
        alphaiw[19] = 0.8891430;
        alphaiw[20] = 0.6203820;
        alphaiw[21] = 0.3088870;
        alphaiw[22] = 0.1831710;
      }
      else if (refsys == 3) {
        alphaiw[0] = 53.3912050;
        alphaiw[1] = 51.4589370;
        alphaiw[2] = 46.7536920;
        alphaiw[3] = 35.7231120;
        alphaiw[4] = 26.8682930;
        alphaiw[5] = 20.5832790;
        alphaiw[6] = 16.1845180;
        alphaiw[7] = 13.0610500;
        alphaiw[8] = 10.7905050;
        alphaiw[9] = 9.0973330;
        alphaiw[10] = 7.8031170;
        alphaiw[11] = 6.7910320;
        alphaiw[12] = 5.3263460;
        alphaiw[13] = 4.3288220;
        alphaiw[14] = 3.6113910;
        alphaiw[15] = 3.0730140;
        alphaiw[16] = 2.6553170;
        alphaiw[17] = 1.9347190;
        alphaiw[18] = 1.4797110;
        alphaiw[19] = 0.9488200;
        alphaiw[20] = 0.6592330;
        alphaiw[21] = 0.3264620;
        alphaiw[22] = 0.1931150;
      }
      else if (refsys == 4) {
        alphaiw[0] = 48.1721000;
        alphaiw[1] = 46.9654470;
        alphaiw[2] = 43.8406630;
        alphaiw[3] = 35.4940690;
        alphaiw[4] = 27.7965480;
        alphaiw[5] = 21.8303470;
        alphaiw[6] = 17.4239690;
        alphaiw[7] = 14.1866160;
        alphaiw[8] = 11.7798410;
        alphaiw[9] = 9.9575870;
        alphaiw[10] = 8.5500830;
        alphaiw[11] = 7.4414740;
        alphaiw[12] = 5.8269560;
        alphaiw[13] = 4.7227970;
        alphaiw[14] = 3.9280720;
        alphaiw[15] = 3.3323920;
        alphaiw[16] = 2.8712630;
        alphaiw[17] = 2.0794750;
        alphaiw[18] = 1.5831290;
        alphaiw[19] = 1.0088640;
        alphaiw[20] = 0.6982730;
        alphaiw[21] = 0.3440760;
        alphaiw[22] = 0.2030650;
      }
      else if (refsys == 5) {
        alphaiw[0] = 44.3809730;
        alphaiw[1] = 43.6657010;
        alphaiw[2] = 41.6751380;
        alphaiw[3] = 35.4583760;
        alphaiw[4] = 28.7575380;
        alphaiw[5] = 23.0706960;
        alphaiw[6] = 18.6492860;
        alphaiw[7] = 15.2988620;
        alphaiw[8] = 12.7583260;
        alphaiw[9] = 10.8093490;
        alphaiw[10] = 9.2904440;
        alphaiw[11] = 8.0867420;
        alphaiw[12] = 6.3242650;
        alphaiw[13] = 5.1145440;
        alphaiw[14] = 4.2431700;
        alphaiw[15] = 3.5905930;
        alphaiw[16] = 3.0863030;
        alphaiw[17] = 2.2237000;
        alphaiw[18] = 1.6862000;
        alphaiw[19] = 1.0687380;
        alphaiw[20] = 0.7372220;
        alphaiw[21] = 0.3616620;
        alphaiw[22] = 0.2130010;
      }
      break;
    case 51:
      if (refsys == 1) {
        alphaiw[0] = 44.8358450;
        alphaiw[1] = 43.5130560;
        alphaiw[2] = 40.0926450;
        alphaiw[3] = 31.1774780;
        alphaiw[4] = 23.4511520;
        alphaiw[5] = 17.8471750;
        alphaiw[6] = 13.9239600;
        alphaiw[7] = 11.1566290;
        alphaiw[8] = 9.1627500;
        alphaiw[9] = 7.6898730;
        alphaiw[10] = 6.5744680;
        alphaiw[11] = 5.7098380;
        alphaiw[12] = 4.4725460;
        alphaiw[13] = 3.6404470;
        alphaiw[14] = 3.0472370;
        alphaiw[15] = 2.6044190;
        alphaiw[16] = 2.2616350;
        alphaiw[17] = 1.6693460;
        alphaiw[18] = 1.2920700;
        alphaiw[19] = 0.8442200;
        alphaiw[20] = 0.5942970;
        alphaiw[21] = 0.2999110;
        alphaiw[22] = 0.1791200;
      }
      else if (refsys == 2) {
        alphaiw[0] = 45.9916410;
        alphaiw[1] = 44.7450290;
        alphaiw[2] = 41.5208600;
        alphaiw[3] = 32.9820310;
        alphaiw[4] = 25.3160810;
        alphaiw[5] = 19.5590470;
        alphaiw[6] = 15.4150740;
        alphaiw[7] = 12.4297230;
        alphaiw[8] = 10.2444980;
        alphaiw[9] = 8.6111660;
        alphaiw[10] = 7.3634730;
        alphaiw[11] = 6.3902070;
        alphaiw[12] = 4.9895020;
        alphaiw[13] = 4.0440290;
        alphaiw[14] = 3.3698710;
        alphaiw[15] = 2.8676280;
        alphaiw[16] = 2.4801170;
        alphaiw[17] = 1.8150710;
        alphaiw[18] = 1.3958700;
        alphaiw[19] = 0.9043070;
        alphaiw[20] = 0.6333220;
        alphaiw[21] = 0.3175040;
        alphaiw[22] = 0.1890550;
      }
      else if (refsys == 3) {
        alphaiw[0] = 46.6496880;
        alphaiw[1] = 45.5197710;
        alphaiw[2] = 42.5927060;
        alphaiw[3] = 34.6039430;
        alphaiw[4] = 27.0775160;
        alphaiw[5] = 21.2059700;
        alphaiw[6] = 16.8636310;
        alphaiw[7] = 13.6742900;
        alphaiw[8] = 11.3067400;
        alphaiw[9] = 9.5188430;
        alphaiw[10] = 8.1427680;
        alphaiw[11] = 7.0635060;
        alphaiw[12] = 5.5024820;
        alphaiw[13] = 4.4452150;
        alphaiw[14] = 3.6909720;
        alphaiw[15] = 3.1298080;
        alphaiw[16] = 2.6978790;
        alphaiw[17] = 1.9604550;
        alphaiw[18] = 1.4994800;
        alphaiw[19] = 0.9643160;
        alphaiw[20] = 0.6723070;
        alphaiw[21] = 0.3350840;
        alphaiw[22] = 0.1989840;
      }
      else if (refsys == 4) {
        alphaiw[0] = 46.7552750;
        alphaiw[1] = 45.8323020;
        alphaiw[2] = 43.3346180;
        alphaiw[3] = 36.0394380;
        alphaiw[4] = 28.7272660;
        alphaiw[5] = 22.7845200;
        alphaiw[6] = 18.2688140;
        alphaiw[7] = 14.8903330;
        alphaiw[8] = 12.3496050;
        alphaiw[9] = 10.4129670;
        alphaiw[10] = 8.9123290;
        alphaiw[11] = 7.7296420;
        alphaiw[12] = 6.0113260;
        alphaiw[13] = 4.8438310;
        alphaiw[14] = 4.0103740;
        alphaiw[15] = 3.3908020;
        alphaiw[16] = 2.9147740;
        alphaiw[17] = 2.1053690;
        alphaiw[18] = 1.6027900;
        alphaiw[19] = 1.0241710;
        alphaiw[20] = 0.7112070;
        alphaiw[21] = 0.3526470;
        alphaiw[22] = 0.2089120;
      }
      break;
    case 52:
      if (refsys == 1) {
        alphaiw[0] = 38.7558670;
        alphaiw[1] = 37.8891050;
        alphaiw[2] = 35.5765770;
        alphaiw[3] = 29.0541470;
        alphaiw[4] = 22.7985480;
        alphaiw[5] = 17.8842260;
        alphaiw[6] = 14.2354170;
        alphaiw[7] = 11.5477160;
        alphaiw[8] = 9.5486370;
        alphaiw[9] = 8.0375580;
        alphaiw[10] = 6.8745420;
        alphaiw[11] = 5.9631130;
        alphaiw[12] = 4.6477480;
        alphaiw[13] = 3.7603970;
        alphaiw[14] = 3.1299390;
        alphaiw[15] = 2.6624490;
        alphaiw[16] = 2.3034360;
        alphaiw[17] = 1.6910060;
        alphaiw[18] = 1.3067460;
        alphaiw[19] = 0.8554460;
        alphaiw[20] = 0.6045040;
        alphaiw[21] = 0.3075880;
        alphaiw[22] = 0.1846460;
      }
      else if (refsys == 2) {
        alphaiw[0] = 41.2800660;
        alphaiw[1] = 40.4036870;
        alphaiw[2] = 38.0648030;
        alphaiw[3] = 31.4147320;
        alphaiw[4] = 24.9165410;
        alphaiw[5] = 19.7095440;
        alphaiw[6] = 15.7789570;
        alphaiw[7] = 12.8463230;
        alphaiw[8] = 10.6436490;
        alphaiw[9] = 8.9662810;
        alphaiw[10] = 7.6680330;
        alphaiw[11] = 6.6463830;
        alphaiw[12] = 5.1661270;
        alphaiw[13] = 4.1647840;
        alphaiw[14] = 3.4530650;
        alphaiw[15] = 2.9259730;
        alphaiw[16] = 2.5221220;
        alphaiw[17] = 1.8367870;
        alphaiw[18] = 1.4105310;
        alphaiw[19] = 0.9154730;
        alphaiw[20] = 0.6434690;
        alphaiw[21] = 0.3251430;
        alphaiw[22] = 0.1945560;
      }
      else if (refsys == 3) {
        alphaiw[0] = 43.3767540;
        alphaiw[1] = 42.5220370;
        alphaiw[2] = 40.2233370;
        alphaiw[3] = 33.5679210;
        alphaiw[4] = 26.9051580;
        alphaiw[5] = 21.4529440;
        alphaiw[6] = 17.2695530;
        alphaiw[7] = 14.1098870;
        alphaiw[8] = 11.7148730;
        alphaiw[9] = 9.8784450;
        alphaiw[10] = 8.4497140;
        alphaiw[11] = 7.3210330;
        alphaiw[12] = 5.6796140;
        alphaiw[13] = 4.5661840;
        alphaiw[14] = 3.7742540;
        alphaiw[15] = 3.1881760;
        alphaiw[16] = 2.7398680;
        alphaiw[17] = 1.9821020;
        alphaiw[18] = 1.5140480;
        alphaiw[19] = 0.9753820;
        alphaiw[20] = 0.6823680;
        alphaiw[21] = 0.3426700;
        alphaiw[22] = 0.2044490;
      }
      break;
    case 53:
      if (refsys == 1) {
        alphaiw[0] = 32.5364480;
        alphaiw[1] = 32.0135530;
        alphaiw[2] = 30.5701660;
        alphaiw[3] = 26.1371840;
        alphaiw[4] = 21.4024280;
        alphaiw[5] = 17.3505230;
        alphaiw[6] = 14.1439250;
        alphaiw[7] = 11.6664030;
        alphaiw[8] = 9.7557500;
        alphaiw[9] = 8.2712240;
        alphaiw[10] = 7.1046010;
        alphaiw[11] = 6.1759850;
        alphaiw[12] = 4.8153900;
        alphaiw[13] = 3.8867100;
        alphaiw[14] = 3.2241390;
        alphaiw[15] = 2.7330620;
        alphaiw[16] = 2.3571070;
        alphaiw[17] = 1.7206990;
        alphaiw[18] = 1.3260350;
        alphaiw[19] = 0.8674990;
        alphaiw[20] = 0.6140980;
        alphaiw[21] = 0.3141040;
        alphaiw[22] = 0.1892860;
      }
      else if (refsys == 2) {
        alphaiw[0] = 35.5939760;
        alphaiw[1] = 35.0391840;
        alphaiw[2] = 33.5058920;
        alphaiw[3] = 28.7725970;
        alphaiw[4] = 23.6664920;
        alphaiw[5] = 19.2510560;
        alphaiw[6] = 15.7268020;
        alphaiw[7] = 12.9864250;
        alphaiw[8] = 10.8631350;
        alphaiw[9] = 9.2076450;
        alphaiw[10] = 7.9032670;
        alphaiw[11] = 6.8629800;
        alphaiw[12] = 5.3359890;
        alphaiw[13] = 4.2925500;
        alphaiw[14] = 3.5482480;
        alphaiw[15] = 2.9972560;
        alphaiw[16] = 2.5762480;
        alphaiw[17] = 1.8666410;
        alphaiw[18] = 1.4298640;
        alphaiw[19] = 0.9275130;
        alphaiw[20] = 0.6530570;
        alphaiw[21] = 0.3316710;
        alphaiw[22] = 0.1992050;
      }
      break;
    case 54:
      if (refsys == 1) {
        alphaiw[0] = 27.3277050;
        alphaiw[1] = 27.0110760;
        alphaiw[2] = 26.1157860;
        alphaiw[3] = 23.1734500;
        alphaiw[4] = 19.7140260;
        alphaiw[5] = 16.4952160;
        alphaiw[6] = 13.7782840;
        alphaiw[7] = 11.5745650;
        alphaiw[8] = 9.8106950;
        alphaiw[9] = 8.3998660;
        alphaiw[10] = 7.2654260;
        alphaiw[11] = 6.3458240;
        alphaiw[12] = 4.9718440;
        alphaiw[13] = 4.0169210;
        alphaiw[14] = 3.3286250;
        alphaiw[15] = 2.8159860;
        alphaiw[16] = 2.4230090;
        alphaiw[17] = 1.7595470;
        alphaiw[18] = 1.3513270;
        alphaiw[19] = 0.8819580;
        alphaiw[20] = 0.6246870;
        alphaiw[21] = 0.3206970;
        alphaiw[22] = 0.1938380;
      }
      break;
    case 55:
      if (refsys == 1) {
        alphaiw[0] = 537.0897540;
        alphaiw[1] = 251.2031120;
        alphaiw[2] = 104.6716850;
        alphaiw[3] = 40.1908390;
        alphaiw[4] = 24.9381150;
        alphaiw[5] = 18.3967650;
        alphaiw[6] = 14.5779860;
        alphaiw[7] = 11.9651660;
        alphaiw[8] = 10.0264350;
        alphaiw[9] = 8.5234170;
        alphaiw[10] = 7.3270290;
        alphaiw[11] = 6.3574360;
        alphaiw[12] = 4.8990430;
        alphaiw[13] = 3.8742590;
        alphaiw[14] = 3.1301840;
        alphaiw[15] = 2.5751510;
        alphaiw[16] = 2.1515440;
        alphaiw[17] = 1.4501300;
        alphaiw[18] = 1.0384390;
        alphaiw[19] = 0.6036570;
        alphaiw[20] = 0.3927120;
        alphaiw[21] = 0.1776630;
        alphaiw[22] = 0.1006850;
      }
      else if (refsys == 2) {
        alphaiw[0] = 65.7018290;
        alphaiw[1] = 59.8918430;
        alphaiw[2] = 49.3703130;
        alphaiw[3] = 33.8726930;
        alphaiw[4] = 25.2509890;
        alphaiw[5] = 19.8576040;
        alphaiw[6] = 16.1061580;
        alphaiw[7] = 13.3278610;
        alphaiw[8] = 11.1931740;
        alphaiw[9] = 9.5136360;
        alphaiw[10] = 8.1692610;
        alphaiw[11] = 7.0782580;
        alphaiw[12] = 5.4392600;
        alphaiw[13] = 4.2914390;
        alphaiw[14] = 3.4609450;
        alphaiw[15] = 2.8433200;
        alphaiw[16] = 2.3731020;
        alphaiw[17] = 1.5968080;
        alphaiw[18] = 1.1424900;
        alphaiw[19] = 0.6636650;
        alphaiw[20] = 0.4316520;
        alphaiw[21] = 0.1952610;
        alphaiw[22] = 0.1106580;
      }
      else if (refsys == 3) {
        alphaiw[0] = 159.0053270;
        alphaiw[1] = 157.1436212;
        alphaiw[2] = 151.9731787;
        alphaiw[3] = 135.7037659;
        alphaiw[4] = 117.3153187;
        alphaiw[5] = 100.2482652;
        alphaiw[6] = 85.4662885;
        alphaiw[7] = 73.0126871;
        alphaiw[8] = 62.6425221;
        alphaiw[9] = 54.0415555;
        alphaiw[10] = 46.9052446;
        alphaiw[11] = 40.9663710;
        alphaiw[12] = 31.8279948;
        alphaiw[13] = 25.2977404;
        alphaiw[14] = 20.5185492;
        alphaiw[15] = 16.9398730;
        alphaiw[16] = 14.2032476;
        alphaiw[17] = 9.6642916;
        alphaiw[18] = 6.9942505;
        alphaiw[19] = 4.1611863;
        alphaiw[20] = 2.7731001;
        alphaiw[21] = 1.3285730;
        alphaiw[22] = 0.7903369;
      }
      break;
    case 56:
      if (refsys == 1) {
        alphaiw[0] = 288.2408830;
        alphaiw[1] = 220.1240810;
        alphaiw[2] = 131.4902760;
        alphaiw[3] = 55.2586600;
        alphaiw[4] = 31.5893280;
        alphaiw[5] = 21.6257080;
        alphaiw[6] = 16.3339130;
        alphaiw[7] = 13.0455460;
        alphaiw[8] = 10.7783530;
        alphaiw[9] = 9.1050960;
        alphaiw[10] = 7.8133370;
        alphaiw[11] = 6.7850620;
        alphaiw[12] = 5.2566000;
        alphaiw[13] = 4.1863800;
        alphaiw[14] = 3.4065330;
        alphaiw[15] = 2.8211540;
        alphaiw[16] = 2.3711900;
        alphaiw[17] = 1.6174560;
        alphaiw[18] = 1.1682080;
        alphaiw[19] = 0.6862840;
        alphaiw[20] = 0.4490360;
        alphaiw[21] = 0.2044320;
        alphaiw[22] = 0.1161510;
      }
      else if (refsys == 2) {
        alphaiw[0] = 309.5716070;
        alphaiw[1] = 175.3545830;
        alphaiw[2] = 94.6066250;
        alphaiw[3] = 46.0830780;
        alphaiw[4] = 30.0147260;
        alphaiw[5] = 22.0684410;
        alphaiw[6] = 17.2886740;
        alphaiw[7] = 14.0639130;
        alphaiw[8] = 11.7253160;
        alphaiw[9] = 9.9472460;
        alphaiw[10] = 8.5510400;
        alphaiw[11] = 7.4290820;
        alphaiw[12] = 5.7511370;
        alphaiw[13] = 4.5735640;
        alphaiw[14] = 3.7161310;
        alphaiw[15] = 3.0735940;
        alphaiw[16] = 2.5806020;
        alphaiw[17] = 1.7570230;
        alphaiw[18] = 1.2676610;
        alphaiw[19] = 0.7440160;
        alphaiw[20] = 0.4866820;
        alphaiw[21] = 0.2215860;
        alphaiw[22] = 0.1259200;
      }
      else if (refsys == 3) {
        alphaiw[0] = 66.9605790;
        alphaiw[1] = 63.0034800;
        alphaiw[2] = 54.8221920;
        alphaiw[3] = 39.7272000;
        alphaiw[4] = 29.6372990;
        alphaiw[5] = 23.0260300;
        alphaiw[6] = 18.4860590;
        alphaiw[7] = 15.2106190;
        alphaiw[8] = 12.7502140;
        alphaiw[9] = 10.8439290;
        alphaiw[10] = 9.3316440;
        alphaiw[11] = 8.1096160;
        alphaiw[12] = 6.2753070;
        alphaiw[13] = 4.9860280;
        alphaiw[14] = 4.0475020;
        alphaiw[15] = 3.3448150;
        alphaiw[16] = 2.8062400;
        alphaiw[17] = 1.9080120;
        alphaiw[18] = 1.3753850;
        alphaiw[19] = 0.8064740;
        alphaiw[20] = 0.5272920;
        alphaiw[21] = 0.2399460;
        alphaiw[22] = 0.1363100;
      }
      else if (refsys == 4) {
        alphaiw[0] = 155.5094396;
        alphaiw[1] = 152.2957198;
        alphaiw[2] = 143.9110133;
        alphaiw[3] = 121.2148688;
        alphaiw[4] = 100.0521440;
        alphaiw[5] = 83.2245940;
        alphaiw[6] = 70.1280996;
        alphaiw[7] = 60.7626607;
        alphaiw[8] = 51.5155269;
        alphaiw[9] = 44.7636011;
        alphaiw[10] = 39.1938022;
        alphaiw[11] = 34.5522205;
        alphaiw[12] = 27.3444515;
        alphaiw[13] = 22.1012480;
        alphaiw[14] = 18.1884416;
        alphaiw[15] = 15.2025580;
        alphaiw[16] = 12.8788334;
        alphaiw[17] = 8.9266992;
        alphaiw[18] = 6.5294153;
        alphaiw[19] = 3.9103194;
        alphaiw[20] = 2.5946296;
        alphaiw[21] = 1.2081345;
        alphaiw[22] = 0.6957280;
      }
      break;
    case 57:
      if (refsys == 1) {
        alphaiw[0] = 215.8140220;
        alphaiw[1] = 176.8449270;
        alphaiw[2] = 118.3158760;
        alphaiw[3] = 56.5506640;
        alphaiw[4] = 33.7025040;
        alphaiw[5] = 23.2152200;
        alphaiw[6] = 17.4335570;
        alphaiw[7] = 13.8111350;
        alphaiw[8] = 11.3289880;
        alphaiw[9] = 9.5173300;
        alphaiw[10] = 8.1343090;
        alphaiw[11] = 7.0437020;
        alphaiw[12] = 5.4382140;
        alphaiw[13] = 4.3231050;
        alphaiw[14] = 3.5139470;
        alphaiw[15] = 2.9079640;
        alphaiw[16] = 2.4428210;
        alphaiw[17] = 1.6648090;
        alphaiw[18] = 1.2017740;
        alphaiw[19] = 0.7056550;
        alphaiw[20] = 0.4616280;
        alphaiw[21] = 0.2101180;
        alphaiw[22] = 0.1193250;
      }
      else if (refsys == 2) {
        alphaiw[0] = 174.1703670;
        alphaiw[1] = 109.8104790;
        alphaiw[2] = 78.4184010;
        alphaiw[3] = 48.3320410;
        alphaiw[4] = 33.8987960;
        alphaiw[5] = 25.4468050;
        alphaiw[6] = 19.9742230;
        alphaiw[7] = 16.1822840;
        alphaiw[8] = 13.4194490;
        alphaiw[9] = 11.3279880;
        alphaiw[10] = 9.6974550;
        alphaiw[11] = 8.3967050;
        alphaiw[12] = 6.4683870;
        alphaiw[13] = 5.1269150;
        alphaiw[14] = 4.1557730;
        alphaiw[15] = 3.4309420;
        alphaiw[16] = 2.8764740;
        alphaiw[17] = 1.9535170;
        alphaiw[18] = 1.4072110;
        alphaiw[19] = 0.8245040;
        alphaiw[20] = 0.5388600;
        alphaiw[21] = 0.2450480;
        alphaiw[22] = 0.1391110;
      }
      else if (refsys == 3) {
        alphaiw[0] = 66.6599540;
        alphaiw[1] = 64.3784460;
        alphaiw[2] = 58.7168170;
        alphaiw[3] = 45.1581050;
        alphaiw[4] = 34.3061210;
        alphaiw[5] = 26.6622110;
        alphaiw[6] = 21.2769290;
        alphaiw[7] = 17.3778510;
        alphaiw[8] = 14.4673110;
        alphaiw[9] = 12.2338560;
        alphaiw[10] = 10.4793230;
        alphaiw[11] = 9.0738430;
        alphaiw[12] = 6.9849920;
        alphaiw[13] = 5.5308960;
        alphaiw[14] = 4.4789930;
        alphaiw[15] = 3.6947640;
        alphaiw[16] = 3.0955470;
        alphaiw[17] = 2.0997350;
        alphaiw[18] = 1.5114130;
        alphaiw[19] = 0.8848790;
        alphaiw[20] = 0.5781220;
        alphaiw[21] = 0.2628230;
        alphaiw[22] = 0.1491850;
      }
      break;
    case 58:
      if (refsys == 1) {
        alphaiw[0] = 63.1803350;
        alphaiw[1] = 60.8655150;
        alphaiw[2] = 55.6077000;
        alphaiw[3] = 42.9730420;
        alphaiw[4] = 32.8160560;
        alphaiw[5] = 25.6410650;
        alphaiw[6] = 20.5717660;
        alphaiw[7] = 16.8914020;
        alphaiw[8] = 14.1378050;
        alphaiw[9] = 12.0210120;
        alphaiw[10] = 10.3558470;
        alphaiw[11] = 9.0204850;
        alphaiw[12] = 7.0330480;
        alphaiw[13] = 5.6466750;
        alphaiw[14] = 4.6410990;
        alphaiw[15] = 3.8887730;
        alphaiw[16] = 3.3113370;
        alphaiw[17] = 2.3420650;
        alphaiw[18] = 1.7576480;
        alphaiw[19] = 1.1111050;
        alphaiw[20] = 0.7737220;
        alphaiw[21] = 0.3918750;
        alphaiw[22] = 0.2365100;
      }
      break;
    case 59:
      if (refsys == 1) {
        alphaiw[0] = 253.6075830;
        alphaiw[1] = 185.3478080;
        alphaiw[2] = 117.3131200;
        alphaiw[3] = 51.4986990;
        alphaiw[4] = 29.5194860;
        alphaiw[5] = 20.1469030;
        alphaiw[6] = 15.2121030;
        alphaiw[7] = 12.1925180;
        alphaiw[8] = 10.1406100;
        alphaiw[9] = 8.6417550;
        alphaiw[10] = 7.4910100;
        alphaiw[11] = 6.5762820;
        alphaiw[12] = 5.2119620;
        alphaiw[13] = 4.2469660;
        alphaiw[14] = 3.5345800;
        alphaiw[15] = 2.9924220;
        alphaiw[16] = 2.5699090;
        alphaiw[17] = 1.8461480;
        alphaiw[18] = 1.3998220;
        alphaiw[19] = 0.8963550;
        alphaiw[20] = 0.6292150;
        alphaiw[21] = 0.3225380;
        alphaiw[22] = 0.1960510;
      }
      else if (refsys == 2) {
        alphaiw[0] = 61.2373490;
        alphaiw[1] = 58.9335430;
        alphaiw[2] = 53.9743400;
        alphaiw[3] = 41.9343430;
        alphaiw[4] = 32.1281250;
        alphaiw[5] = 25.1481370;
        alphaiw[6] = 20.2000000;
        alphaiw[7] = 16.6031610;
        alphaiw[8] = 13.9109010;
        alphaiw[9] = 11.8406070;
        alphaiw[10] = 10.2112030;
        alphaiw[11] = 8.9035100;
        alphaiw[12] = 6.9541880;
        alphaiw[13] = 5.5908190;
        alphaiw[14] = 4.5991840;
        alphaiw[15] = 3.8553720;
        alphaiw[16] = 3.2832050;
        alphaiw[17] = 2.3202580;
        alphaiw[18] = 1.7385060;
        alphaiw[19] = 1.0952390;
        alphaiw[20] = 0.7607110;
        alphaiw[21] = 0.3842190;
        alphaiw[22] = 0.2318250;
      }
      break;
    case 60:
      if (refsys == 1) {
        alphaiw[0] = 224.5762180;
        alphaiw[1] = 179.2684400;
        alphaiw[2] = 114.5509860;
        alphaiw[3] = 50.7711900;
        alphaiw[4] = 29.0959830;
        alphaiw[5] = 19.8086410;
        alphaiw[6] = 14.9251210;
        alphaiw[7] = 11.9483050;
        alphaiw[8] = 9.9338740;
        alphaiw[9] = 8.4673600;
        alphaiw[10] = 7.3439400;
        alphaiw[11] = 6.4519280;
        alphaiw[12] = 5.1215350;
        alphaiw[13] = 4.1790220;
        alphaiw[14] = 3.4814770;
        alphaiw[15] = 2.9491840;
        alphaiw[16] = 2.5333290;
        alphaiw[17] = 1.8187150;
        alphaiw[18] = 1.3768350;
        alphaiw[19] = 0.8783290;
        alphaiw[20] = 0.6146620;
        alphaiw[21] = 0.3138300;
        alphaiw[22] = 0.1905480;
      }
      else if (refsys == 2) {
        alphaiw[0] = 59.8908920;
        alphaiw[1] = 57.8663720;
        alphaiw[2] = 53.0214580;
        alphaiw[3] = 41.2308960;
        alphaiw[4] = 31.6026170;
        alphaiw[5] = 24.7396280;
        alphaiw[6] = 19.8738380;
        alphaiw[7] = 16.3389450;
        alphaiw[8] = 13.6950980;
        alphaiw[9] = 11.6632880;
        alphaiw[10] = 10.0646540;
        alphaiw[11] = 8.7815920;
        alphaiw[12] = 6.8677690;
        alphaiw[13] = 5.5271730;
        alphaiw[14] = 4.5502380;
        alphaiw[15] = 3.8160260;
        alphaiw[16] = 3.2502340;
        alphaiw[17] = 2.2958240;
        alphaiw[18] = 1.7180490;
        alphaiw[19] = 1.0790660;
        alphaiw[20] = 0.7475530;
        alphaiw[21] = 0.3762800;
        alphaiw[22] = 0.2268030;
      }
      break;
    case 61:
      if (refsys == 1) {
        alphaiw[0] = 216.5540280;
        alphaiw[1] = 174.0685600;
        alphaiw[2] = 112.0657200;
        alphaiw[3] = 49.9521540;
        alphaiw[4] = 28.5964740;
        alphaiw[5] = 19.4194800;
        alphaiw[6] = 14.6027420;
        alphaiw[7] = 11.6781260;
        alphaiw[8] = 9.7072940;
        alphaiw[9] = 8.2773840;
        alphaiw[10] = 7.1844350;
        alphaiw[11] = 6.3175570;
        alphaiw[12] = 5.0245550;
        alphaiw[13] = 4.1068450;
        alphaiw[14] = 3.4257530;
        alphaiw[15] = 2.9044660;
        alphaiw[16] = 2.4960780;
        alphaiw[17] = 1.7917290;
        alphaiw[18] = 1.3546830;
        alphaiw[19] = 0.8610840;
        alphaiw[20] = 0.6005800;
        alphaiw[21] = 0.3050960;
        alphaiw[22] = 0.1848870;
      }
      else if (refsys == 2) {
        alphaiw[0] = 57.9857080;
        alphaiw[1] = 56.0520210;
        alphaiw[2] = 51.4036550;
        alphaiw[3] = 40.1715870;
        alphaiw[4] = 30.9043600;
        alphaiw[5] = 24.2360810;
        alphaiw[6] = 19.4864830;
        alphaiw[7] = 16.0303390;
        alphaiw[8] = 13.4446280;
        alphaiw[9] = 11.4577250;
        alphaiw[10] = 9.8945360;
        alphaiw[11] = 8.6397290;
        alphaiw[12] = 6.7667310;
        alphaiw[13] = 5.4526280;
        alphaiw[14] = 4.4930820;
        alphaiw[15] = 3.7704490;
        alphaiw[16] = 3.2125050;
        alphaiw[17] = 2.2688980;
        alphaiw[18] = 1.6962120;
        alphaiw[19] = 1.0622980;
        alphaiw[20] = 0.7339400;
        alphaiw[21] = 0.3678790;
        alphaiw[22] = 0.2213700;
      }
      break;
    case 62:
      if (refsys == 1) {
        alphaiw[0] = 209.4821630;
        alphaiw[1] = 169.3264270;
        alphaiw[2] = 109.7496270;
        alphaiw[3] = 49.1677750;
        alphaiw[4] = 28.1221600;
        alphaiw[5] = 19.0508730;
        alphaiw[6] = 14.2954700;
        alphaiw[7] = 11.4178440;
        alphaiw[8] = 9.4862320;
        alphaiw[9] = 8.0895200;
        alphaiw[10] = 7.0245330;
        alphaiw[11] = 6.1810260;
        alphaiw[12] = 4.9234940;
        alphaiw[13] = 4.0299890;
        alphaiw[14] = 3.3654240;
        alphaiw[15] = 2.8555240;
        alphaiw[16] = 2.4550980;
        alphaiw[17] = 1.7622060;
        alphaiw[18] = 1.3308620;
        alphaiw[19] = 0.8430750;
        alphaiw[20] = 0.5860700;
        alphaiw[21] = 0.2961160;
        alphaiw[22] = 0.1790130;
      }
      else if (refsys == 2) {
        alphaiw[0] = 60.4975550;
        alphaiw[1] = 56.8914260;
        alphaiw[2] = 51.7310830;
        alphaiw[3] = 40.1100700;
        alphaiw[4] = 30.7211290;
        alphaiw[5] = 24.0196790;
        alphaiw[6] = 19.2720050;
        alphaiw[7] = 15.8320160;
        alphaiw[8] = 13.2674100;
        alphaiw[9] = 11.3021990;
        alphaiw[10] = 9.7592770;
        alphaiw[11] = 8.5225050;
        alphaiw[12] = 6.6783810;
        alphaiw[13] = 5.3848030;
        alphaiw[14] = 4.4395640;
        alphaiw[15] = 3.7268880;
        alphaiw[16] = 3.1759310;
        alphaiw[17] = 2.2423610;
        alphaiw[18] = 1.6746190;
        alphaiw[19] = 1.0457120;
        alphaiw[20] = 0.7204260;
        alphaiw[21] = 0.3593990;
        alphaiw[22] = 0.2157950;
      }
      break;
    case 63:
      if (refsys == 1) {
        alphaiw[0] = 203.2111900;
        alphaiw[1] = 165.0774090;
        alphaiw[2] = 107.6941350;
        alphaiw[3] = 48.4835280;
        alphaiw[4] = 27.7173710;
        alphaiw[5] = 18.7379040;
        alphaiw[6] = 14.0329820;
        alphaiw[7] = 11.1932480;
        alphaiw[8] = 9.2933470;
        alphaiw[9] = 7.9237680;
        alphaiw[10] = 6.8819260;
        alphaiw[11] = 6.0580060;
        alphaiw[12] = 4.8307140;
        alphaiw[13] = 3.9582910;
        alphaiw[14] = 3.3084100;
        alphaiw[15] = 2.8088210;
        alphaiw[16] = 2.4157350;
        alphaiw[17] = 1.7336990;
        alphaiw[18] = 1.3079430;
        alphaiw[19] = 0.8259450;
        alphaiw[20] = 0.5723790;
        alphaiw[21] = 0.2877150;
        alphaiw[22] = 0.1735300;
      }
      else if (refsys == 2) {
        alphaiw[0] = 58.3846820;
        alphaiw[1] = 55.5074830;
        alphaiw[2] = 48.9596600;
        alphaiw[3] = 35.5239180;
        alphaiw[4] = 26.2422590;
        alphaiw[5] = 20.2406920;
        alphaiw[6] = 16.2113610;
        alphaiw[7] = 13.3691670;
        alphaiw[8] = 11.2734740;
        alphaiw[9] = 9.6705860;
        alphaiw[10] = 8.4080470;
        alphaiw[11] = 7.3900200;
        alphaiw[12] = 5.8559040;
        alphaiw[13] = 4.7635490;
        alphaiw[14] = 3.9543760;
        alphaiw[15] = 3.3370760;
        alphaiw[16] = 2.8551160;
        alphaiw[17] = 2.0281870;
        alphaiw[18] = 1.5187450;
        alphaiw[19] = 0.9491330;
        alphaiw[20] = 0.6530540;
        alphaiw[21] = 0.3246600;
        alphaiw[22] = 0.1946040;
      }
      break;
    case 64:
      if (refsys == 1) {
        alphaiw[0] = 141.8420640;
        alphaiw[1] = 116.8571040;
        alphaiw[2] = 82.4531150;
        alphaiw[3] = 42.8638860;
        alphaiw[4] = 26.4230410;
        alphaiw[5] = 18.5963810;
        alphaiw[6] = 14.2465310;
        alphaiw[7] = 11.5139220;
        alphaiw[8] = 9.6332520;
        alphaiw[9] = 8.2500980;
        alphaiw[10] = 7.1830900;
        alphaiw[11] = 6.3311290;
        alphaiw[12] = 5.0514610;
        alphaiw[13] = 4.1367300;
        alphaiw[14] = 3.4541820;
        alphaiw[15] = 2.9295050;
        alphaiw[16] = 2.5169890;
        alphaiw[17] = 1.8023770;
        alphaiw[18] = 1.3572570;
        alphaiw[19] = 0.8545860;
        alphaiw[20] = 0.5909780;
        alphaiw[21] = 0.2962180;
        alphaiw[22] = 0.1784510;
      }
      else if (refsys == 2) {
        alphaiw[0] = 59.9521870;
        alphaiw[1] = 57.6299620;
        alphaiw[2] = 52.3452910;
        alphaiw[3] = 40.4588930;
        alphaiw[4] = 30.8901980;
        alphaiw[5] = 24.0499010;
        alphaiw[6] = 19.2156040;
        alphaiw[7] = 15.7302290;
        alphaiw[8] = 13.1464310;
        alphaiw[9] = 11.1770130;
        alphaiw[10] = 9.6378050;
        alphaiw[11] = 8.4085080;
        alphaiw[12] = 6.5823240;
        alphaiw[13] = 5.3050770;
        alphaiw[14] = 4.3727770;
        alphaiw[15] = 3.6698530;
        alphaiw[16] = 3.1261280;
        alphaiw[17] = 2.2035920;
        alphaiw[18] = 1.6416690;
        alphaiw[19] = 1.0191840;
        alphaiw[20] = 0.6982840;
        alphaiw[21] = 0.3451090;
        alphaiw[22] = 0.2062760;
      }
      break;
    case 65:
      if (refsys == 1) {
        alphaiw[0] = 229.3015960;
        alphaiw[1] = 148.7572060;
        alphaiw[2] = 83.8992590;
        alphaiw[3] = 37.7525280;
        alphaiw[4] = 23.0411680;
        alphaiw[5] = 16.4030610;
        alphaiw[6] = 12.7361960;
        alphaiw[7] = 10.4183080;
        alphaiw[8] = 8.8072510;
        alphaiw[9] = 7.6095540;
        alphaiw[10] = 6.6757310;
        alphaiw[11] = 5.9225500;
        alphaiw[12] = 4.7760310;
        alphaiw[13] = 3.9431820;
        alphaiw[14] = 3.3133000;
        alphaiw[15] = 2.8236490;
        alphaiw[16] = 2.4350780;
        alphaiw[17] = 1.7540070;
        alphaiw[18] = 1.3245690;
        alphaiw[19] = 0.8351270;
        alphaiw[20] = 0.5770300;
        alphaiw[21] = 0.2882710;
        alphaiw[22] = 0.1732830;
      }
      else if (refsys == 2) {
        alphaiw[0] = 55.7777620;
        alphaiw[1] = 53.9540020;
        alphaiw[2] = 49.4627260;
        alphaiw[3] = 38.5259830;
        alphaiw[4] = 29.5462910;
        alphaiw[5] = 23.1049190;
        alphaiw[6] = 18.5294200;
        alphaiw[7] = 15.2134150;
        alphaiw[8] = 12.7449150;
        alphaiw[9] = 10.8576310;
        alphaiw[10] = 9.3793010;
        alphaiw[11] = 8.1965330;
        alphaiw[12] = 6.4356150;
        alphaiw[13] = 5.2004630;
        alphaiw[14] = 4.2962800;
        alphaiw[15] = 3.6126080;
        alphaiw[16] = 3.0823350;
        alphaiw[17] = 2.1790430;
        alphaiw[18] = 1.6261760;
        alphaiw[19] = 1.0110080;
        alphaiw[20] = 0.6927930;
        alphaiw[21] = 0.3420770;
        alphaiw[22] = 0.2042870;
      }
      break;
    case 66:
      if (refsys == 1) {
        alphaiw[0] = 232.4656730;
        alphaiw[1] = 135.1836000;
        alphaiw[2] = 80.1286950;
        alphaiw[3] = 37.7486970;
        alphaiw[4] = 23.2434510;
        alphaiw[5] = 16.5330840;
        alphaiw[6] = 12.7969630;
        alphaiw[7] = 10.4341740;
        alphaiw[8] = 8.7966470;
        alphaiw[9] = 7.5842480;
        alphaiw[10] = 6.6428660;
        alphaiw[11] = 5.8863450;
        alphaiw[12] = 4.7393580;
        alphaiw[13] = 3.9091410;
        alphaiw[14] = 3.2824410;
        alphaiw[15] = 2.7957150;
        alphaiw[16] = 2.4096230;
        alphaiw[17] = 1.7330220;
        alphaiw[18] = 1.3064720;
        alphaiw[19] = 0.8207680;
        alphaiw[20] = 0.5653110;
        alphaiw[21] = 0.2809900;
        alphaiw[22] = 0.1685370;
      }
      else if (refsys == 2) {
        alphaiw[0] = 60.4090430;
        alphaiw[1] = 55.9872550;
        alphaiw[2] = 50.3819660;
        alphaiw[3] = 38.5219880;
        alphaiw[4] = 29.3573150;
        alphaiw[5] = 22.8882210;
        alphaiw[6] = 18.3247830;
        alphaiw[7] = 15.0311820;
        alphaiw[8] = 12.5861420;
        alphaiw[9] = 10.7203330;
        alphaiw[10] = 9.2606660;
        alphaiw[11] = 8.0937380;
        alphaiw[12] = 6.3571910;
        alphaiw[13] = 5.1389330;
        alphaiw[14] = 4.2464670;
        alphaiw[15] = 3.5710020;
        alphaiw[16] = 3.0465710;
        alphaiw[17] = 2.1519680;
        alphaiw[18] = 1.6036830;
        alphaiw[19] = 0.9935670;
        alphaiw[20] = 0.6786190;
        alphaiw[21] = 0.3332680;
        alphaiw[22] = 0.1985350;
      }
      break;
    case 67:
      if (refsys == 1) {
        alphaiw[0] = 174.3579710;
        alphaiw[1] = 145.3475730;
        alphaiw[2] = 98.3386140;
        alphaiw[3] = 45.6653390;
        alphaiw[4] = 26.1535730;
        alphaiw[5] = 17.5697070;
        alphaiw[6] = 13.0686310;
        alphaiw[7] = 10.3734350;
        alphaiw[8] = 8.5909150;
        alphaiw[9] = 7.3206360;
        alphaiw[10] = 6.3633450;
        alphaiw[11] = 5.6111750;
        alphaiw[12] = 4.4954000;
        alphaiw[13] = 3.7014620;
        alphaiw[14] = 3.1066380;
        alphaiw[15] = 2.6458480;
        alphaiw[16] = 2.2803750;
        alphaiw[17] = 1.6387150;
        alphaiw[18] = 1.2329370;
        alphaiw[19] = 0.7701030;
        alphaiw[20] = 0.5271880;
        alphaiw[21] = 0.2591200;
        alphaiw[22] = 0.1545480;
      }
      else if (refsys == 2) {
        alphaiw[0] = 52.4790850;
        alphaiw[1] = 50.9455980;
        alphaiw[2] = 47.0541390;
        alphaiw[3] = 37.1641820;
        alphaiw[4] = 28.6862320;
        alphaiw[5] = 22.4764090;
        alphaiw[6] = 18.0300540;
        alphaiw[7] = 14.8000680;
        alphaiw[8] = 12.3955370;
        alphaiw[9] = 10.5586820;
        alphaiw[10] = 9.1212810;
        alphaiw[11] = 7.9722170;
        alphaiw[12] = 6.2625620;
        alphaiw[13] = 5.0632080;
        alphaiw[14] = 4.1843110;
        alphaiw[15] = 3.5187310;
        alphaiw[16] = 3.0016100;
        alphaiw[17] = 2.1184470;
        alphaiw[18] = 1.5764470;
        alphaiw[19] = 0.9730470;
        alphaiw[20] = 0.6620460;
        alphaiw[21] = 0.3227820;
        alphaiw[22] = 0.1915480;
      }
      break;
    case 68:
      if (refsys == 1) {
        alphaiw[0] = 169.5086700;
        alphaiw[1] = 141.8311300;
        alphaiw[2] = 96.4854630;
        alphaiw[3] = 44.9898390;
        alphaiw[4] = 25.7372660;
        alphaiw[5] = 17.2457070;
        alphaiw[6] = 12.7960440;
        alphaiw[7] = 10.1383860;
        alphaiw[8] = 8.3867110;
        alphaiw[9] = 7.1428030;
        alphaiw[10] = 6.2082610;
        alphaiw[11] = 5.4756830;
        alphaiw[12] = 4.3910800;
        alphaiw[13] = 3.6198210;
        alphaiw[14] = 3.0414510;
        alphaiw[15] = 2.5926580;
        alphaiw[16] = 2.2360160;
        alphaiw[17] = 1.6079710;
        alphaiw[18] = 1.2093860;
        alphaiw[19] = 0.7535860;
        alphaiw[20] = 0.5142640;
        alphaiw[21] = 0.2510880;
        alphaiw[22] = 0.1491740;
      }
      else if (refsys == 2) {
        alphaiw[0] = 50.9243700;
        alphaiw[1] = 49.4869020;
        alphaiw[2] = 45.8178160;
        alphaiw[3] = 36.3744420;
        alphaiw[4] = 28.1606380;
        alphaiw[5] = 22.0904510;
        alphaiw[6] = 17.7252050;
        alphaiw[7] = 14.5487700;
        alphaiw[8] = 12.1834110;
        alphaiw[9] = 10.3771780;
        alphaiw[10] = 8.9646510;
        alphaiw[11] = 7.8362150;
        alphaiw[12] = 6.1584380;
        alphaiw[13] = 4.9818960;
        alphaiw[14] = 4.1194840;
        alphaiw[15] = 3.4659310;
        alphaiw[16] = 2.9576820;
        alphaiw[17] = 2.0882320;
        alphaiw[18] = 1.5534830;
        alphaiw[19] = 0.9571210;
        alphaiw[20] = 0.6496440;
        alphaiw[21] = 0.3150810;
        alphaiw[22] = 0.1863840;
      }
      break;
    case 69:
      if (refsys == 1) {
        alphaiw[0] = 164.8744790;
        alphaiw[1] = 138.5034170;
        alphaiw[2] = 94.7797640;
        alphaiw[3] = 44.4385580;
        alphaiw[4] = 25.4312120;
        alphaiw[5] = 17.0158970;
        alphaiw[6] = 12.6028120;
        alphaiw[7] = 9.9697610;
        alphaiw[8] = 8.2378580;
        alphaiw[9] = 7.0109450;
        alphaiw[10] = 6.0912830;
        alphaiw[11] = 5.3717500;
        alphaiw[12] = 4.3084370;
        alphaiw[13] = 3.5531670;
        alphaiw[14] = 2.9867270;
        alphaiw[15] = 2.5468500;
        alphaiw[16] = 2.1969270;
        alphaiw[17] = 1.5796200;
        alphaiw[18] = 1.1870370;
        alphaiw[19] = 0.7375970;
        alphaiw[20] = 0.5018220;
        alphaiw[21] = 0.2436440;
        alphaiw[22] = 0.1443670;
      }
      else if (refsys == 2) {
        alphaiw[0] = 50.9258390;
        alphaiw[1] = 49.3627740;
        alphaiw[2] = 45.6624910;
        alphaiw[3] = 36.1857370;
        alphaiw[4] = 27.9704110;
        alphaiw[5] = 21.9122080;
        alphaiw[6] = 17.5625150;
        alphaiw[7] = 14.4018520;
        alphaiw[8] = 12.0514680;
        alphaiw[9] = 10.2590760;
        alphaiw[10] = 8.8591280;
        alphaiw[11] = 7.7419770;
        alphaiw[12] = 6.0830330;
        alphaiw[13] = 4.9209200;
        alphaiw[14] = 4.0693970;
        alphaiw[15] = 3.4240270;
        alphaiw[16] = 2.9219510;
        alphaiw[17] = 2.0623180;
        alphaiw[18] = 1.5329870;
        alphaiw[19] = 0.9422880;
        alphaiw[20] = 0.6379740;
        alphaiw[21] = 0.3079830;
        alphaiw[22] = 0.1817680;
      }
      break;
    case 70:
      if (refsys == 1) {
        alphaiw[0] = 163.7379070;
        alphaiw[1] = 137.4016000;
        alphaiw[2] = 93.8094750;
        alphaiw[3] = 43.7828040;
        alphaiw[4] = 24.9508110;
        alphaiw[5] = 16.6301360;
        alphaiw[6] = 12.2781120;
        alphaiw[7] = 9.6905770;
        alphaiw[8] = 7.9955170;
        alphaiw[9] = 6.7996180;
        alphaiw[10] = 5.9064770;
        alphaiw[11] = 5.2097330;
        alphaiw[12] = 4.1828480;
        alphaiw[13] = 3.4544290;
        alphaiw[14] = 2.9077830;
        alphaiw[15] = 2.4825740;
        alphaiw[16] = 2.1436170;
        alphaiw[17] = 1.5435770;
        alphaiw[18] = 1.1603320;
        alphaiw[19] = 0.7200560;
        alphaiw[20] = 0.4887560;
        alphaiw[21] = 0.2360610;
        alphaiw[22] = 0.1394570;
      }
      else if (refsys == 2) {
        alphaiw[0] = 54.1811670;
        alphaiw[1] = 52.3584220;
        alphaiw[2] = 48.1909980;
        alphaiw[3] = 37.7829980;
        alphaiw[4] = 28.8680310;
        alphaiw[5] = 22.3694590;
        alphaiw[6] = 17.7640250;
        alphaiw[7] = 14.4601420;
        alphaiw[8] = 12.0316690;
        alphaiw[9] = 10.1983910;
        alphaiw[10] = 8.7787560;
        alphaiw[11] = 7.6539450;
        alphaiw[12] = 5.9967880;
        alphaiw[13] = 4.8445850;
        alphaiw[14] = 4.0040030;
        alphaiw[15] = 3.3683890;
        alphaiw[16] = 2.8743980;
        alphaiw[17] = 2.0286710;
        alphaiw[18] = 1.5073600;
        alphaiw[19] = 0.9248310;
        alphaiw[20] = 0.6246590;
        alphaiw[21] = 0.3000110;
        alphaiw[22] = 0.1765320;
      }
      break;
    case 71:
      if (refsys == 1) {
        alphaiw[0] = 136.8530670;
        alphaiw[1] = 104.2975620;
        alphaiw[2] = 76.5562240;
        alphaiw[3] = 40.2834860;
        alphaiw[4] = 24.1515190;
        alphaiw[5] = 16.4764390;
        alphaiw[6] = 12.3100420;
        alphaiw[7] = 9.7790490;
        alphaiw[8] = 8.0987470;
        alphaiw[9] = 6.9030530;
        alphaiw[10] = 6.0050390;
        alphaiw[11] = 5.3018750;
        alphaiw[12] = 4.2621270;
        alphaiw[13] = 3.5227910;
        alphaiw[14] = 2.9673480;
        alphaiw[15] = 2.5350740;
        alphaiw[16] = 2.1903720;
        alphaiw[17] = 1.5797660;
        alphaiw[18] = 1.1891950;
        alphaiw[19] = 0.7393030;
        alphaiw[20] = 0.5021560;
        alphaiw[21] = 0.2423920;
        alphaiw[22] = 0.1430370;
      }
      else if (refsys == 2) {
        alphaiw[0] = 48.9452970;
        alphaiw[1] = 47.6317230;
        alphaiw[2] = 44.2461460;
        alphaiw[3] = 35.3576300;
        alphaiw[4] = 27.4507560;
        alphaiw[5] = 21.5311680;
        alphaiw[6] = 17.2508500;
        alphaiw[7] = 14.1326890;
        alphaiw[8] = 11.8136070;
        alphaiw[9] = 10.0470280;
        alphaiw[10] = 8.6694960;
        alphaiw[11] = 7.5721550;
        alphaiw[12] = 5.9462740;
        alphaiw[13] = 4.8097520;
        alphaiw[14] = 3.9777670;
        alphaiw[15] = 3.3472210;
        alphaiw[16] = 2.8563980;
        alphaiw[17] = 2.0147310;
        alphaiw[18] = 1.4952460;
        alphaiw[19] = 0.9145890;
        alphaiw[20] = 0.6157550;
        alphaiw[21] = 0.2938290;
        alphaiw[22] = 0.1722380;
      }
      break;
    case 72:
      if (refsys == 1) {
        alphaiw[0] = 106.2816210;
        alphaiw[1] = 95.4891650;
        alphaiw[2] = 75.3034870;
        alphaiw[3] = 44.1393100;
        alphaiw[4] = 28.1355730;
        alphaiw[5] = 19.6843670;
        alphaiw[6] = 14.7522820;
        alphaiw[7] = 11.6120870;
        alphaiw[8] = 9.4717210;
        alphaiw[9] = 7.9330520;
        alphaiw[10] = 6.7792700;
        alphaiw[11] = 5.8843260;
        alphaiw[12] = 4.5893730;
        alphaiw[13] = 3.6999000;
        alphaiw[14] = 3.0536720;
        alphaiw[15] = 2.5653680;
        alphaiw[16] = 2.1856620;
        alphaiw[17] = 1.5348190;
        alphaiw[18] = 1.1334250;
        alphaiw[19] = 0.6859220;
        alphaiw[20] = 0.4568140;
        alphaiw[21] = 0.2123490;
        alphaiw[22] = 0.1215650;
      }
      else if (refsys == 2) {
        alphaiw[0] = 82.4536570;
        alphaiw[1] = 70.5184470;
        alphaiw[2] = 56.7403060;
        alphaiw[3] = 38.8993780;
        alphaiw[4] = 28.3733700;
        alphaiw[5] = 21.5865770;
        alphaiw[6] = 16.9965970;
        alphaiw[7] = 13.7660430;
        alphaiw[8] = 11.4103040;
        alphaiw[9] = 9.6381650;
        alphaiw[10] = 8.2684630;
        alphaiw[11] = 7.1847890;
        alphaiw[12] = 5.5917520;
        alphaiw[13] = 4.4881220;
        alphaiw[14] = 3.6863220;
        alphaiw[15] = 3.0827340;
        alphaiw[16] = 2.6157700;
        alphaiw[17] = 1.8222580;
        alphaiw[18] = 1.3383190;
        alphaiw[19] = 0.8045950;
        alphaiw[20] = 0.5339520;
        alphaiw[21] = 0.2472490;
        alphaiw[22] = 0.1413400;
      }
      else if (refsys == 3) {
        alphaiw[0] = 50.2657100;
        alphaiw[1] = 49.2005530;
        alphaiw[2] = 46.3604080;
        alphaiw[3] = 38.3345580;
        alphaiw[4] = 30.5672240;
        alphaiw[5] = 24.3724470;
        alphaiw[6] = 19.6879820;
        alphaiw[7] = 16.1687550;
        alphaiw[8] = 13.4975460;
        alphaiw[9] = 11.4368740;
        alphaiw[10] = 9.8188790;
        alphaiw[11] = 8.5264490;
        alphaiw[12] = 6.6139230;
        alphaiw[13] = 5.2863370;
        alphaiw[14] = 4.3241350;
        alphaiw[15] = 3.6027470;
        alphaiw[16] = 3.0471630;
        alphaiw[17] = 2.1096360;
        alphaiw[18] = 1.5428050;
        alphaiw[19] = 0.9227930;
        alphaiw[20] = 0.6106910;
        alphaiw[21] = 0.2819190;
        alphaiw[22] = 0.1609740;
      }
      else if (refsys == 4) {
        alphaiw[0] = 32.7095574;
        alphaiw[1] = 32.4889582;
        alphaiw[2] = 31.8527540;
        alphaiw[3] = 29.6320422;
        alphaiw[4] = 26.7617762;
        alphaiw[5] = 23.8245539;
        alphaiw[6] = 21.1219459;
        alphaiw[7] = 18.7537050;
        alphaiw[8] = 16.7198153;
        alphaiw[9] = 14.9828615;
        alphaiw[10] = 13.4971673;
        alphaiw[11] = 12.2201143;
        alphaiw[12] = 10.1536814;
        alphaiw[13] = 8.5678059;
        alphaiw[14] = 7.3223689;
        alphaiw[15] = 6.3253685;
        alphaiw[16] = 5.5145527;
        alphaiw[17] = 4.0455469;
        alphaiw[18] = 3.0833508;
        alphaiw[19] = 1.9496957;
        alphaiw[20] = 1.3382618;
        alphaiw[21] = 0.6527061;
        alphaiw[22] = 0.3847355;
      }
      break;
    case 73:
      if (refsys == 1) {
        alphaiw[0] = 88.9256120;
        alphaiw[1] = 81.9318660;
        alphaiw[2] = 67.3277280;
        alphaiw[3] = 41.9543150;
        alphaiw[4] = 27.7130410;
        alphaiw[5] = 19.8298060;
        alphaiw[6] = 15.0833980;
        alphaiw[7] = 11.9897820;
        alphaiw[8] = 9.8427490;
        alphaiw[9] = 8.2781960;
        alphaiw[10] = 7.0935190;
        alphaiw[11] = 6.1684980;
        alphaiw[12] = 4.8227300;
        alphaiw[13] = 3.8952850;
        alphaiw[14] = 3.2207380;
        alphaiw[15] = 2.7107210;
        alphaiw[16] = 2.3138110;
        alphaiw[17] = 1.6320700;
        alphaiw[18] = 1.2099340;
        alphaiw[19] = 0.7365720;
        alphaiw[20] = 0.4925290;
        alphaiw[21] = 0.2302270;
        alphaiw[22] = 0.1321400;
      }
      else if (refsys == 2) {
        alphaiw[0] = 69.3801270;
        alphaiw[1] = 62.6525490;
        alphaiw[2] = 52.3321990;
        alphaiw[3] = 37.2824880;
        alphaiw[4] = 27.8804010;
        alphaiw[5] = 21.6007000;
        alphaiw[6] = 17.2276620;
        alphaiw[7] = 14.0778720;
        alphaiw[8] = 11.7401200;
        alphaiw[9] = 9.9582810;
        alphaiw[10] = 8.5679630;
        alphaiw[11] = 7.4606290;
        alphaiw[12] = 5.8231300;
        alphaiw[13] = 4.6837700;
        alphaiw[14] = 3.8543220;
        alphaiw[15] = 3.2291230;
        alphaiw[16] = 2.7448770;
        alphaiw[17] = 1.9201830;
        alphaiw[18] = 1.4153120;
        alphaiw[19] = 0.8555430;
        alphaiw[20] = 0.5698730;
        alphaiw[21] = 0.2652310;
        alphaiw[22] = 0.1519760;
      }
      else if (refsys == 3) {
        alphaiw[0] = 53.1343510;
        alphaiw[1] = 50.1385310;
        alphaiw[2] = 45.3135370;
        alphaiw[3] = 36.1734470;
        alphaiw[4] = 28.5997780;
        alphaiw[5] = 22.7928480;
        alphaiw[6] = 18.4525250;
        alphaiw[7] = 15.2001580;
        alphaiw[8] = 12.7291810;
        alphaiw[9] = 10.8187030;
        alphaiw[10] = 9.3147510;
        alphaiw[11] = 8.1103590;
        alphaiw[12] = 6.3220940;
        alphaiw[13] = 5.0755740;
        alphaiw[14] = 4.1686360;
        alphaiw[15] = 3.4861320;
        alphaiw[16] = 2.9585540;
        alphaiw[17] = 2.0630420;
        alphaiw[18] = 1.5171890;
        alphaiw[19] = 0.9145820;
        alphaiw[20] = 0.6082630;
        alphaiw[21] = 0.2826130;
        alphaiw[22] = 0.1618310;
      }
      else if (refsys == 4) {
        alphaiw[0] = 38.0167314;
        alphaiw[1] = 37.7505450;
        alphaiw[2] = 36.9868016;
        alphaiw[3] = 34.3581519;
        alphaiw[4] = 31.0182858;
        alphaiw[5] = 27.6354658;
        alphaiw[6] = 24.5337826;
        alphaiw[7] = 21.8145582;
        alphaiw[8] = 19.4741943;
        alphaiw[9] = 17.4704983;
        alphaiw[10] = 15.7528316;
        alphaiw[11] = 14.2738200;
        alphaiw[12] = 11.8766483;
        alphaiw[13] = 10.0346794;
        alphaiw[14] = 8.5871211;
        alphaiw[15] = 7.4275650;
        alphaiw[16] = 6.4837949;
        alphaiw[17] = 4.7709367;
        alphaiw[18] = 3.6456491;
        alphaiw[19] = 2.3144054;
        alphaiw[20] = 1.5929899;
        alphaiw[21] = 0.7800977;
        alphaiw[22] = 0.4607918;
      }
      break;
    case 74:
      if (refsys == 1) {
        alphaiw[0] = 72.6835210;
        alphaiw[1] = 66.5725390;
        alphaiw[2] = 54.5218660;
        alphaiw[3] = 35.1042250;
        alphaiw[4] = 24.3743050;
        alphaiw[5] = 18.1885020;
        alphaiw[6] = 14.2758210;
        alphaiw[7] = 11.6114380;
        alphaiw[8] = 9.6940500;
        alphaiw[9] = 8.2552810;
        alphaiw[10] = 7.1400870;
        alphaiw[11] = 6.2531290;
        alphaiw[12] = 4.9373670;
        alphaiw[13] = 4.0140520;
        alphaiw[14] = 3.3348740;
        alphaiw[15] = 2.8174560;
        alphaiw[16] = 2.4125540;
        alphaiw[17] = 1.7123560;
        alphaiw[18] = 1.2753760;
        alphaiw[19] = 0.7814210;
        alphaiw[20] = 0.5246740;
        alphaiw[21] = 0.2466040;
        alphaiw[22] = 0.1418950;
      }
      else if (refsys == 2) {
        alphaiw[0] = 57.9535230;
        alphaiw[1] = 54.3254100;
        alphaiw[2] = 47.3716060;
        alphaiw[3] = 35.2077750;
        alphaiw[4] = 26.9133880;
        alphaiw[5] = 21.1900880;
        alphaiw[6] = 17.1120840;
        alphaiw[7] = 14.1184440;
        alphaiw[8] = 11.8610010;
        alphaiw[9] = 10.1176970;
        alphaiw[10] = 8.7429580;
        alphaiw[11] = 7.6387710;
        alphaiw[12] = 5.9912880;
        alphaiw[13] = 4.8355390;
        alphaiw[14] = 3.9899060;
        alphaiw[15] = 3.3503580;
        alphaiw[16] = 2.8537360;
        alphaiw[17] = 2.0050600;
        alphaiw[18] = 1.4831160;
        alphaiw[19] = 0.9012110;
        alphaiw[20] = 0.6023820;
        alphaiw[21] = 0.2816950;
        alphaiw[22] = 0.1617650;
      }
      else if (refsys == 3) {
        alphaiw[0] = 69.9811980;
        alphaiw[1] = 68.7138070;
        alphaiw[2] = 65.3381200;
        alphaiw[3] = 55.7596720;
        alphaiw[4] = 46.2137800;
        alphaiw[5] = 38.1741940;
        alphaiw[6] = 31.7064340;
        alphaiw[7] = 26.5625450;
        alphaiw[8] = 22.4690330;
        alphaiw[9] = 19.1923040;
        alphaiw[10] = 16.5475190;
        alphaiw[11] = 14.3927290;
        alphaiw[12] = 11.1480040;
        alphaiw[13] = 8.8712130;
        alphaiw[14] = 7.2191730;
        alphaiw[15] = 5.9850930;
        alphaiw[16] = 5.0399870;
        alphaiw[17] = 3.4609450;
        alphaiw[18] = 2.5183910;
        alphaiw[19] = 1.4990770;
        alphaiw[20] = 0.9903650;
        alphaiw[21] = 0.4569140;
        alphaiw[22] = 0.2609890;
      }
      break;
    case 75:
      if (refsys == 1) {
        alphaiw[0] = 70.9067150;
        alphaiw[1] = 65.7037140;
        alphaiw[2] = 55.2189840;
        alphaiw[3] = 36.8309030;
        alphaiw[4] = 25.6365300;
        alphaiw[5] = 19.0080410;
        alphaiw[6] = 14.8342930;
        alphaiw[7] = 12.0225100;
        alphaiw[8] = 10.0181560;
        alphaiw[9] = 8.5243540;
        alphaiw[10] = 7.3716210;
        alphaiw[11] = 6.4572650;
        alphaiw[12] = 5.1034730;
        alphaiw[13] = 4.1544350;
        alphaiw[14] = 3.4564620;
        alphaiw[15] = 2.9246120;
        alphaiw[16] = 2.5081950;
        alphaiw[17] = 1.7869670;
        alphaiw[18] = 1.3354090;
        alphaiw[19] = 0.8224190;
        alphaiw[20] = 0.5541320;
        alphaiw[21] = 0.2617050;
        alphaiw[22] = 0.1509270;
      }
      else if (refsys == 2) {
        alphaiw[0] = 63.2047780;
        alphaiw[1] = 58.0228720;
        alphaiw[2] = 48.8925050;
        alphaiw[3] = 34.6221050;
        alphaiw[4] = 25.8457800;
        alphaiw[5] = 20.1387260;
        alphaiw[6] = 16.1920480;
        alphaiw[7] = 13.3396770;
        alphaiw[8] = 11.2072140;
        alphaiw[9] = 9.5682700;
        alphaiw[10] = 8.2790220;
        alphaiw[11] = 7.2445840;
        alphaiw[12] = 5.7011200;
        alphaiw[13] = 4.6167910;
        alphaiw[14] = 3.8216470;
        alphaiw[15] = 3.2187300;
        alphaiw[16] = 2.7492670;
        alphaiw[17] = 1.9432630;
        alphaiw[18] = 1.4442130;
        alphaiw[19] = 0.8833880;
        alphaiw[20] = 0.5929380;
        alphaiw[21] = 0.2787900;
        alphaiw[22] = 0.1604910;
      }
      else if (refsys == 3) {
        alphaiw[0] = 61.9391170;
        alphaiw[1] = 61.0033150;
        alphaiw[2] = 58.4522600;
        alphaiw[3] = 50.8001250;
        alphaiw[4] = 42.6897350;
        alphaiw[5] = 35.5953060;
        alphaiw[6] = 29.7653350;
        alphaiw[7] = 25.0690300;
        alphaiw[8] = 21.2990290;
        alphaiw[9] = 18.2612110;
        alphaiw[10] = 15.7959520;
        alphaiw[11] = 13.7781660;
        alphaiw[12] = 10.7227960;
        alphaiw[13] = 8.5655510;
        alphaiw[14] = 6.9926320;
        alphaiw[15] = 5.8130460;
        alphaiw[16] = 4.9067170;
        alphaiw[17] = 3.3857530;
        alphaiw[18] = 2.4728560;
        alphaiw[19] = 1.4797000;
        alphaiw[20] = 0.9808630;
        alphaiw[21] = 0.4545450;
        alphaiw[22] = 0.2601570;
      }
      break;
    case 76:
      if (refsys == 1) {
        alphaiw[0] = 59.6083520;
        alphaiw[1] = 56.2806020;
        alphaiw[2] = 48.8648650;
        alphaiw[3] = 34.1600550;
        alphaiw[4] = 24.4374970;
        alphaiw[5] = 18.4318460;
        alphaiw[6] = 14.5524170;
        alphaiw[7] = 11.8947470;
        alphaiw[8] = 9.9775490;
        alphaiw[9] = 8.5355150;
        alphaiw[10] = 7.4142660;
        alphaiw[11] = 6.5190680;
        alphaiw[12] = 5.1830340;
        alphaiw[13] = 4.2381040;
        alphaiw[14] = 3.5384920;
        alphaiw[15] = 3.0026950;
        alphaiw[16] = 2.5815390;
        alphaiw[17] = 1.8484760;
        alphaiw[18] = 1.3867840;
        alphaiw[19] = 0.8588690;
        alphaiw[20] = 0.5808420;
        alphaiw[21] = 0.2757360;
        alphaiw[22] = 0.1594300;
      }
      else if (refsys == 2) {
        alphaiw[0] = 46.6731420;
        alphaiw[1] = 44.2599930;
        alphaiw[2] = 40.9119690;
        alphaiw[3] = 32.7985670;
        alphaiw[4] = 25.8416100;
        alphaiw[5] = 20.6679230;
        alphaiw[6] = 16.8773970;
        alphaiw[7] = 14.0534960;
        alphaiw[8] = 11.9009790;
        alphaiw[9] = 10.2230790;
        alphaiw[10] = 8.8884890;
        alphaiw[11] = 7.8079950;
        alphaiw[12] = 6.1792880;
        alphaiw[13] = 5.0231940;
        alphaiw[14] = 4.1694910;
        alphaiw[15] = 3.5191180;
        alphaiw[16] = 3.0110530;
        alphaiw[17] = 2.1356650;
        alphaiw[18] = 1.5915890;
        alphaiw[19] = 0.9776240;
        alphaiw[20] = 0.6581110;
        alphaiw[21] = 0.3107410;
        alphaiw[22] = 0.1792660;
      }
      else if (refsys == 3) {
        alphaiw[0] = 55.8892400;
        alphaiw[1] = 55.1763930;
        alphaiw[2] = 53.1894280;
        alphaiw[3] = 46.9039320;
        alphaiw[4] = 39.8479720;
        alphaiw[5] = 33.4618590;
        alphaiw[6] = 28.1191240;
        alphaiw[7] = 23.7725720;
        alphaiw[8] = 20.2616430;
        alphaiw[9] = 17.4199010;
        alphaiw[10] = 15.1054600;
        alphaiw[11] = 13.2052620;
        alphaiw[12] = 10.3167870;
        alphaiw[13] = 8.2680590;
        alphaiw[14] = 6.7685540;
        alphaiw[15] = 5.6403850;
        alphaiw[16] = 4.7711490;
        alphaiw[17] = 3.3068100;
        alphaiw[18] = 2.4236090;
        alphaiw[19] = 1.4575360;
        alphaiw[20] = 0.9693850;
        alphaiw[21] = 0.4512540;
        alphaiw[22] = 0.2588140;
      }
      break;
    case 77:
      if (refsys == 1) {
        alphaiw[0] = 51.5493740;
        alphaiw[1] = 49.2234270;
        alphaiw[2] = 43.7106450;
        alphaiw[3] = 31.7225470;
        alphaiw[4] = 23.2089290;
        alphaiw[5] = 17.7566590;
        alphaiw[6] = 14.1614650;
        alphaiw[7] = 11.6653050;
        alphaiw[8] = 9.8468960;
        alphaiw[9] = 8.4684340;
        alphaiw[10] = 7.3894730;
        alphaiw[11] = 6.5229570;
        alphaiw[12] = 5.2199510;
        alphaiw[13] = 4.2901010;
        alphaiw[14] = 3.5965860;
        alphaiw[15] = 3.0623010;
        alphaiw[16] = 2.6403160;
        alphaiw[17] = 1.9013040;
        alphaiw[18] = 1.4325800;
        alphaiw[19] = 0.8926780;
        alphaiw[20] = 0.6061140;
        alphaiw[21] = 0.2892620;
        alphaiw[22] = 0.1676580;
      }
      else if (refsys == 2) {
        alphaiw[0] = 42.7200080;
        alphaiw[1] = 41.3931980;
        alphaiw[2] = 38.3655160;
        alphaiw[3] = 31.0080660;
        alphaiw[4] = 24.6653650;
        alphaiw[5] = 19.9066380;
        alphaiw[6] = 16.3877970;
        alphaiw[7] = 13.7425380;
        alphaiw[8] = 11.7090270;
        alphaiw[9] = 10.1115020;
        alphaiw[10] = 8.8318110;
        alphaiw[11] = 7.7890930;
        alphaiw[12] = 6.2044270;
        alphaiw[13] = 5.0689240;
        alphaiw[14] = 4.2241180;
        alphaiw[15] = 3.5767130;
        alphaiw[16] = 3.0685910;
        alphaiw[17] = 2.1879850;
        alphaiw[18] = 1.6370730;
        alphaiw[19] = 1.0112500;
        alphaiw[20] = 0.6832680;
        alphaiw[21] = 0.3242250;
        alphaiw[22] = 0.1874750;
      }
      else if (refsys == 3) {
        alphaiw[0] = 40.1998410;
        alphaiw[1] = 39.0125180;
        alphaiw[2] = 36.3126410;
        alphaiw[3] = 30.1829760;
        alphaiw[4] = 24.8286850;
        alphaiw[5] = 20.5587890;
        alphaiw[6] = 17.2259500;
        alphaiw[7] = 14.6163560;
        alphaiw[8] = 12.5495230;
        alphaiw[9] = 10.8901200;
        alphaiw[10] = 9.5396530;
        alphaiw[11] = 8.4265790;
        alphaiw[12] = 6.7165300;
        alphaiw[13] = 5.4810030;
        alphaiw[14] = 4.5587850;
        alphaiw[15] = 3.8517520;
        alphaiw[16] = 3.2974500;
        alphaiw[17] = 2.3400930;
        alphaiw[18] = 1.7446250;
        alphaiw[19] = 1.0727640;
        alphaiw[20] = 0.7229650;
        alphaiw[21] = 0.3420560;
        alphaiw[22] = 0.1975510;
      }
      else if (refsys == 4) {
        alphaiw[0] = 43.2638615;
        alphaiw[1] = 41.8727089;
        alphaiw[2] = 38.9775692;
        alphaiw[3] = 32.5354682;
        alphaiw[4] = 26.8687902;
        alphaiw[5] = 22.3369089;
        alphaiw[6] = 18.7845502;
        alphaiw[7] = 15.9853930;
        alphaiw[8] = 13.7519014;
        alphaiw[9] = 11.9455547;
        alphaiw[10] = 10.4662014;
        alphaiw[11] = 9.2409010;
        alphaiw[12] = 7.3500768;
        alphaiw[13] = 5.9807691;
        alphaiw[14] = 4.9594258;
        alphaiw[15] = 4.1782069;
        alphaiw[16] = 3.5675888;
        alphaiw[17] = 2.5185219;
        alphaiw[18] = 1.8706839;
        alphaiw[19] = 1.1451485;
        alphaiw[20] = 0.7699740;
        alphaiw[21] = 0.3635000;
        alphaiw[22] = 0.2098079;
      }
      else if (refsys == 5) {
        alphaiw[0] = 49.3567100;
        alphaiw[1] = 48.1701200;
        alphaiw[2] = 45.3609500;
        alphaiw[3] = 38.0546700;
        alphaiw[4] = 31.0840500;
        alphaiw[5] = 25.4572400;
        alphaiw[6] = 21.1008600;
        alphaiw[7] = 17.7296600;
        alphaiw[8] = 15.0901300;
        alphaiw[9] = 12.9932400;
        alphaiw[10] = 11.3032400;
        alphaiw[11] = 9.9227200;
        alphaiw[12] = 7.8265900;
        alphaiw[13] = 6.3334100;
        alphaiw[14] = 5.2320700;
        alphaiw[15] = 4.3961000;
        alphaiw[16] = 3.7462100;
        alphaiw[17] = 2.6360300;
        alphaiw[18] = 1.9540400;
        alphaiw[19] = 1.1932200;
        alphaiw[20] = 0.8010500;
        alphaiw[21] = 0.3773100;
        alphaiw[22] = 0.2175100;
      }
      break;
    case 78:
      if (refsys == 1) {
        alphaiw[0] = 41.1424820;
        alphaiw[1] = 39.3180510;
        alphaiw[2] = 35.0148090;
        alphaiw[3] = 25.8710280;
        alphaiw[4] = 19.5216660;
        alphaiw[5] = 15.4256100;
        alphaiw[6] = 12.6507620;
        alphaiw[7] = 10.6590770;
        alphaiw[8] = 9.1609870;
        alphaiw[9] = 7.9929620;
        alphaiw[10] = 7.0565850;
        alphaiw[11] = 6.2892190;
        alphaiw[12] = 5.1074490;
        alphaiw[13] = 4.2422610;
        alphaiw[14] = 3.5845540;
        alphaiw[15] = 3.0705210;
        alphaiw[16] = 2.6600560;
        alphaiw[17] = 1.9319420;
        alphaiw[18] = 1.4641410;
        alphaiw[19] = 0.9193200;
        alphaiw[20] = 0.6271910;
        alphaiw[21] = 0.3012380;
        alphaiw[22] = 0.1751400;
      }
      else if (refsys == 2) {
        alphaiw[0] = 40.5572290;
        alphaiw[1] = 39.3599010;
        alphaiw[2] = 36.5944450;
        alphaiw[3] = 29.7359090;
        alphaiw[4] = 23.7699580;
        alphaiw[5] = 19.2758940;
        alphaiw[6] = 15.9444510;
        alphaiw[7] = 13.4329020;
        alphaiw[8] = 11.4953860;
        alphaiw[9] = 9.9672740;
        alphaiw[10] = 8.7381500;
        alphaiw[11] = 7.7325090;
        alphaiw[12] = 6.1953270;
        alphaiw[13] = 5.0856920;
        alphaiw[14] = 4.2547590;
        alphaiw[15] = 3.6144440;
        alphaiw[16] = 3.1095350;
        alphaiw[17] = 2.2291660;
        alphaiw[18] = 1.6746000;
        alphaiw[19] = 1.0403320;
        alphaiw[20] = 0.7055790;
        alphaiw[21] = 0.3365800;
        alphaiw[22] = 0.1951260;
      }
      else if (refsys == 3) {
        alphaiw[0] = 54.9491410;
        alphaiw[1] = 53.9578470;
        alphaiw[2] = 51.3088730;
        alphaiw[3] = 43.7409950;
        alphaiw[4] = 36.2053630;
        alphaiw[5] = 29.9365010;
        alphaiw[6] = 24.9624590;
        alphaiw[7] = 21.0441630;
        alphaiw[8] = 17.9395030;
        alphaiw[9] = 15.4543750;
        alphaiw[10] = 13.4423040;
        alphaiw[11] = 11.7944940;
        alphaiw[12] = 9.2888530;
        alphaiw[13] = 7.5040460;
        alphaiw[14] = 6.1891800;
        alphaiw[15] = 5.1927300;
        alphaiw[16] = 4.4193840;
        alphaiw[17] = 3.1019380;
        alphaiw[18] = 2.2955520;
        alphaiw[19] = 1.3992410;
        alphaiw[20] = 0.9386560;
        alphaiw[21] = 0.4419430;
        alphaiw[22] = 0.2547950;
      }
      break;
    case 79:
      if (refsys == 1) {
        alphaiw[0] = 37.7290460;
        alphaiw[1] = 36.1091840;
        alphaiw[2] = 32.2722340;
        alphaiw[3] = 24.0716250;
        alphaiw[4] = 18.3385860;
        alphaiw[5] = 14.6143630;
        alphaiw[6] = 12.0774260;
        alphaiw[7] = 10.2474000;
        alphaiw[8] = 8.8635320;
        alphaiw[9] = 7.7783310;
        alphaiw[10] = 6.9031680;
        alphaiw[11] = 6.1817050;
        alphaiw[12] = 5.0613040;
        alphaiw[13] = 4.2321700;
        alphaiw[14] = 3.5957550;
        alphaiw[15] = 3.0941580;
        alphaiw[16] = 2.6907340;
        alphaiw[17] = 1.9683670;
        alphaiw[18] = 1.4994510;
        alphaiw[19] = 0.9481020;
        alphaiw[20] = 0.6497440;
        alphaiw[21] = 0.3139430;
        alphaiw[22] = 0.1830170;
      }
      else if (refsys == 2) {
        alphaiw[0] = 36.6153320;
        alphaiw[1] = 35.4187520;
        alphaiw[2] = 32.5691470;
        alphaiw[3] = 25.9574740;
        alphaiw[4] = 20.5780280;
        alphaiw[5] = 16.6720520;
        alphaiw[6] = 13.8374270;
        alphaiw[7] = 11.7244450;
        alphaiw[8] = 10.1016410;
        alphaiw[9] = 8.8215770;
        alphaiw[10] = 7.7887820;
        alphaiw[11] = 6.9395760;
        alphaiw[12] = 5.6295900;
        alphaiw[13] = 4.6707850;
        alphaiw[14] = 3.9429090;
        alphaiw[15] = 3.3749190;
        alphaiw[16] = 2.9220340;
        alphaiw[17] = 2.1203720;
        alphaiw[18] = 1.6065850;
        alphaiw[19] = 1.0093190;
        alphaiw[20] = 0.6892520;
        alphaiw[21] = 0.3316860;
        alphaiw[22] = 0.1930410;
      }
      break;
    case 80:
      if (refsys == 1) {
        alphaiw[0] = 37.1668150;
        alphaiw[1] = 35.9903170;
        alphaiw[2] = 32.9876470;
        alphaiw[3] = 25.5237630;
        alphaiw[4] = 19.5407590;
        alphaiw[5] = 15.4358360;
        alphaiw[6] = 12.6124980;
        alphaiw[7] = 10.5949600;
        alphaiw[8] = 9.0927810;
        alphaiw[9] = 7.9334090;
        alphaiw[10] = 7.0116660;
        alphaiw[11] = 6.2609290;
        alphaiw[12] = 5.1108470;
        alphaiw[13] = 4.2707020;
        alphaiw[14] = 3.6306640;
        alphaiw[15] = 3.1280940;
        alphaiw[16] = 2.7244100;
        alphaiw[17] = 2.0009140;
        alphaiw[18] = 1.5295470;
        alphaiw[19] = 0.9722000;
        alphaiw[20] = 0.6686340;
        alphaiw[21] = 0.3246570;
        alphaiw[22] = 0.1896630;
      }
      else if (refsys == 2) {
        alphaiw[0] = 35.5596890;
        alphaiw[1] = 34.9453750;
        alphaiw[2] = 33.2600440;
        alphaiw[3] = 28.2046230;
        alphaiw[4] = 23.0451020;
        alphaiw[5] = 18.8299000;
        alphaiw[6] = 15.6024300;
        alphaiw[7] = 13.1484490;
        alphaiw[8] = 11.2575930;
        alphaiw[9] = 9.7724990;
        alphaiw[10] = 8.5831260;
        alphaiw[11] = 7.6133220;
        alphaiw[12] = 6.1348920;
        alphaiw[13] = 5.0678290;
        alphaiw[14] = 4.2663670;
        alphaiw[15] = 3.6457370;
        alphaiw[16] = 3.1535120;
        alphaiw[17] = 2.2866750;
        alphaiw[18] = 1.7331020;
        alphaiw[19] = 1.0902810;
        alphaiw[20] = 0.7455680;
        alphaiw[21] = 0.3596280;
        alphaiw[22] = 0.2095320;
      }
      break;
    case 81:
      if (refsys == 1) {
        alphaiw[0] = 75.4787320;
        alphaiw[1] = 66.9797680;
        alphaiw[2] = 52.7293310;
        alphaiw[3] = 33.2800020;
        alphaiw[4] = 23.1647510;
        alphaiw[5] = 17.3598020;
        alphaiw[6] = 13.7123950;
        alphaiw[7] = 11.2540670;
        alphaiw[8] = 9.5015210;
        alphaiw[9] = 8.1942640;
        alphaiw[10] = 7.1828260;
        alphaiw[11] = 6.3767320;
        alphaiw[12] = 5.1706100;
        alphaiw[13] = 4.3089140;
        alphaiw[14] = 3.6612600;
        alphaiw[15] = 3.1566110;
        alphaiw[16] = 2.7528300;
        alphaiw[17] = 2.0300880;
        alphaiw[18] = 1.5579660;
        alphaiw[19] = 0.9964200;
        alphaiw[20] = 0.6883510;
        alphaiw[21] = 0.3365720;
        alphaiw[22] = 0.1974540;
      }
      else if (refsys == 2) {
        alphaiw[0] = 73.7909520;
        alphaiw[1] = 67.2077030;
        alphaiw[2] = 55.1633950;
        alphaiw[3] = 36.6818870;
        alphaiw[4] = 26.1091350;
        alphaiw[5] = 19.7299200;
        alphaiw[6] = 15.6034900;
        alphaiw[7] = 12.7753400;
        alphaiw[8] = 10.7414800;
        alphaiw[9] = 9.2191530;
        alphaiw[10] = 8.0412980;
        alphaiw[11] = 7.1046070;
        alphaiw[12] = 5.7107080;
        alphaiw[13] = 4.7240350;
        alphaiw[14] = 3.9895440;
        alphaiw[15] = 3.4223730;
        alphaiw[16] = 2.9722160;
        alphaiw[17] = 2.1752900;
        alphaiw[18] = 1.6611010;
        alphaiw[19] = 1.0561290;
        alphaiw[20] = 0.7272340;
        alphaiw[21] = 0.3542350;
        alphaiw[22] = 0.2074800;
      }
      else if (refsys == 3) {
        alphaiw[0] = 51.2159910;
        alphaiw[1] = 49.0618590;
        alphaiw[2] = 44.5281440;
        alphaiw[3] = 34.5954290;
        alphaiw[4] = 26.6152220;
        alphaiw[5] = 20.8897950;
        alphaiw[6] = 16.8343740;
        alphaiw[7] = 13.9090070;
        alphaiw[8] = 11.7406120;
        alphaiw[9] = 10.0876440;
        alphaiw[10] = 8.7945760;
        alphaiw[11] = 7.7596980;
        alphaiw[12] = 6.2135910;
        alphaiw[13] = 5.1188320;
        alphaiw[14] = 4.3061810;
        alphaiw[15] = 3.6812160;
        alphaiw[16] = 3.1873790;
        alphaiw[17] = 2.3191740;
        alphaiw[18] = 1.7638260;
        alphaiw[19] = 1.1158640;
        alphaiw[20] = 0.7661990;
        alphaiw[21] = 0.3719620;
        alphaiw[22] = 0.2175450;
      }
      else if (refsys == 4) {
        alphaiw[0] = 44.4232940;
        alphaiw[1] = 43.5695950;
        alphaiw[2] = 41.2468990;
        alphaiw[3] = 34.4223640;
        alphaiw[4] = 27.6300950;
        alphaiw[5] = 22.1919640;
        alphaiw[6] = 18.1026400;
        alphaiw[7] = 15.0481770;
        alphaiw[8] = 12.7356250;
        alphaiw[9] = 10.9495660;
        alphaiw[10] = 9.5411710;
        alphaiw[11] = 8.4086990;
        alphaiw[12] = 6.7118360;
        alphaiw[13] = 5.5101390;
        alphaiw[14] = 4.6201450;
        alphaiw[15] = 3.9379650;
        alphaiw[16] = 3.4008700;
        alphaiw[17] = 2.4620390;
        alphaiw[18] = 1.8658870;
        alphaiw[19] = 1.1752790;
        alphaiw[20] = 0.8049870;
        alphaiw[21] = 0.3896270;
        alphaiw[22] = 0.2275810;
      }
      break;
    case 82:
      if (refsys == 1) {
        alphaiw[0] = 63.2905510;
        alphaiw[1] = 59.5286430;
        alphaiw[2] = 51.2607860;
        alphaiw[3] = 35.4025850;
        alphaiw[4] = 25.1394190;
        alphaiw[5] = 18.7947900;
        alphaiw[6] = 14.7040930;
        alphaiw[7] = 11.9322080;
        alphaiw[8] = 9.9660590;
        alphaiw[9] = 8.5145470;
        alphaiw[10] = 7.4056250;
        alphaiw[11] = 6.5334430;
        alphaiw[12] = 5.2518040;
        alphaiw[13] = 4.3552580;
        alphaiw[14] = 3.6919830;
        alphaiw[15] = 3.1808980;
        alphaiw[16] = 2.7750050;
        alphaiw[17] = 2.0530120;
        alphaiw[18] = 1.5823440;
        alphaiw[19] = 1.0201520;
        alphaiw[20] = 0.7090500;
        alphaiw[21] = 0.3497720;
        alphaiw[22] = 0.2060120;
      }
      else if (refsys == 2) {
        alphaiw[0] = 57.4725650;
        alphaiw[1] = 54.8466790;
        alphaiw[2] = 48.8575900;
        alphaiw[3] = 36.1749940;
        alphaiw[4] = 26.8405620;
        alphaiw[5] = 20.5425950;
        alphaiw[6] = 16.2594910;
        alphaiw[7] = 13.2605030;
        alphaiw[8] = 11.0893760;
        alphaiw[9] = 9.4661890;
        alphaiw[10] = 8.2168000;
        alphaiw[11] = 7.2301920;
        alphaiw[12] = 5.7780700;
        alphaiw[13] = 4.7643580;
        alphaiw[14] = 4.0179560;
        alphaiw[15] = 3.4461300;
        alphaiw[16] = 2.9946910;
        alphaiw[17] = 2.1989640;
        alphaiw[18] = 1.6860470;
        alphaiw[19] = 1.0800310;
        alphaiw[20] = 0.7479170;
        alphaiw[21] = 0.3673300;
        alphaiw[22] = 0.2159550;
      }
      else if (refsys == 3) {
        alphaiw[0] = 53.0414450;
        alphaiw[1] = 51.0751360;
        alphaiw[2] = 46.8366890;
        alphaiw[3] = 36.7990540;
        alphaiw[4] = 28.3107680;
        alphaiw[5] = 22.1045460;
        alphaiw[6] = 17.6826350;
        alphaiw[7] = 14.4967630;
        alphaiw[8] = 12.1481600;
        alphaiw[9] = 10.3717920;
        alphaiw[10] = 8.9944640;
        alphaiw[11] = 7.9020540;
        alphaiw[12] = 6.2898170;
        alphaiw[13] = 5.1644070;
        alphaiw[14] = 4.3379580;
        alphaiw[15] = 3.7072310;
        alphaiw[16] = 3.2113990;
        alphaiw[17] = 2.3434100;
        alphaiw[18] = 1.7888720;
        alphaiw[19] = 1.1395320;
        alphaiw[20] = 0.7865880;
        alphaiw[21] = 0.3848350;
        alphaiw[22] = 0.2258800;
      }
      else if (refsys == 4) {
        alphaiw[0] = 53.1115280;
        alphaiw[1] = 51.7198390;
        alphaiw[2] = 48.1733620;
        alphaiw[3] = 38.8959490;
        alphaiw[4] = 30.4802980;
        alphaiw[5] = 24.0395110;
        alphaiw[6] = 19.3203190;
        alphaiw[7] = 15.8625500;
        alphaiw[8] = 13.2876830;
        alphaiw[9] = 11.3287250;
        alphaiw[10] = 9.8050600;
        alphaiw[11] = 8.5950040;
        alphaiw[12] = 6.8096590;
        alphaiw[13] = 5.5666540;
        alphaiw[14] = 4.6574470;
        alphaiw[15] = 3.9666200;
        alphaiw[16] = 3.4259370;
        alphaiw[17] = 2.4857500;
        alphaiw[18] = 1.8900620;
        alphaiw[19] = 1.1981450;
        alphaiw[20] = 0.8247670;
        alphaiw[21] = 0.4021920;
        alphaiw[22] = 0.2357430;
      }
      else if (refsys == 5) {
        alphaiw[0] = 50.0879780;
        alphaiw[1] = 49.2042630;
        alphaiw[2] = 46.7719060;
        alphaiw[3] = 39.3982430;
        alphaiw[4] = 31.7542100;
        alphaiw[5] = 25.4539200;
        alphaiw[6] = 20.6426050;
        alphaiw[7] = 17.0297710;
        alphaiw[8] = 14.2980180;
        alphaiw[9] = 12.1992850;
        alphaiw[10] = 10.5566110;
        alphaiw[11] = 9.2468670;
        alphaiw[12] = 7.3088260;
        alphaiw[13] = 5.9581810;
        alphaiw[14] = 4.9713460;
        alphaiw[15] = 4.2231630;
        alphaiw[16] = 3.6391300;
        alphaiw[17] = 2.6281690;
        alphaiw[18] = 1.9916070;
        alphaiw[19] = 1.2570530;
        alphaiw[20] = 0.8631290;
        alphaiw[21] = 0.4196110;
        alphaiw[22] = 0.2456340;
      }
      break;
    case 83:
      if (refsys == 1) {
        alphaiw[0] = 51.1916750;
        alphaiw[1] = 49.4069930;
        alphaiw[2] = 44.9458450;
        alphaiw[3] = 34.1505430;
        alphaiw[4] = 25.5045870;
        alphaiw[5] = 19.5073500;
        alphaiw[6] = 15.3898030;
        alphaiw[7] = 12.4998210;
        alphaiw[8] = 10.4110930;
        alphaiw[9] = 8.8561170;
        alphaiw[10] = 7.6660440;
        alphaiw[11] = 6.7322080;
        alphaiw[12] = 5.3700990;
        alphaiw[13] = 4.4292980;
        alphaiw[14] = 3.7418080;
        alphaiw[15] = 3.2175040;
        alphaiw[16] = 2.8044220;
        alphaiw[17] = 2.0757420;
        alphaiw[18] = 1.6033620;
        alphaiw[19] = 1.0390450;
        alphaiw[20] = 0.7253430;
        alphaiw[21] = 0.3603050;
        alphaiw[22] = 0.2129720;
      }
      else if (refsys == 2) {
        alphaiw[0] = 51.6147800;
        alphaiw[1] = 50.0629560;
        alphaiw[2] = 46.1204110;
        alphaiw[3] = 36.1339090;
        alphaiw[4] = 27.6118910;
        alphaiw[5] = 21.4077340;
        alphaiw[6] = 17.0101490;
        alphaiw[7] = 13.8592580;
        alphaiw[8] = 11.5510200;
        alphaiw[9] = 9.8174290;
        alphaiw[10] = 8.4832410;
        alphaiw[11] = 7.4329030;
        alphaiw[12] = 5.8982610;
        alphaiw[13] = 4.8394220;
        alphaiw[14] = 4.0683880;
        alphaiw[15] = 3.4831190;
        alphaiw[16] = 3.0243540;
        alphaiw[17] = 2.2217570;
        alphaiw[18] = 1.7070420;
        alphaiw[19] = 1.0988490;
        alphaiw[20] = 0.7641350;
        alphaiw[21] = 0.3778150;
        alphaiw[22] = 0.2228850;
      }
      else if (refsys == 3) {
        alphaiw[0] = 51.7479250;
        alphaiw[1] = 50.4406230;
        alphaiw[2] = 47.0464190;
        alphaiw[3] = 37.9441110;
        alphaiw[4] = 29.6097990;
        alphaiw[5] = 23.2406220;
        alphaiw[6] = 18.5881070;
        alphaiw[7] = 15.1912320;
        alphaiw[8] = 12.6725630;
        alphaiw[9] = 10.7660510;
        alphaiw[10] = 9.2914280;
        alphaiw[11] = 8.1270340;
        alphaiw[12] = 6.4226970;
        alphaiw[13] = 5.2472620;
        alphaiw[14] = 4.3934730;
        alphaiw[15] = 3.7477030;
        alphaiw[16] = 3.2435420;
        alphaiw[17] = 2.3673900;
        alphaiw[18] = 1.8104910;
        alphaiw[19] = 1.1585430;
        alphaiw[20] = 0.8028630;
        alphaiw[21] = 0.3953020;
        alphaiw[22] = 0.2327880;
      }
      else if (refsys == 4) {
        alphaiw[0] = 51.5178550;
        alphaiw[1] = 50.4932930;
        alphaiw[2] = 47.7156380;
        alphaiw[3] = 39.5948620;
        alphaiw[4] = 31.5084490;
        alphaiw[5] = 25.0123180;
        alphaiw[6] = 20.1276560;
        alphaiw[7] = 16.4983230;
        alphaiw[8] = 13.7774350;
        alphaiw[9] = 11.7031450;
        alphaiw[10] = 10.0914150;
        alphaiw[11] = 8.8151820;
        alphaiw[12] = 6.9437270;
        alphaiw[13] = 5.6530110;
        alphaiw[14] = 4.7171930;
        alphaiw[15] = 4.0113480;
        alphaiw[16] = 3.4620540;
        alphaiw[17] = 2.5126730;
        alphaiw[18] = 1.9137260;
        alphaiw[19] = 1.2181320;
        alphaiw[20] = 0.8415320;
        alphaiw[21] = 0.4127710;
        alphaiw[22] = 0.2426840;
      }
      break;
    case 84:
      if (refsys == 1) {
        alphaiw[0] = 45.7354800;
        alphaiw[1] = 44.5177550;
        alphaiw[2] = 41.3554830;
        alphaiw[3] = 32.9793560;
        alphaiw[4] = 25.5092680;
        alphaiw[5] = 19.9242900;
        alphaiw[6] = 15.8977680;
        alphaiw[7] = 12.9790260;
        alphaiw[8] = 10.8242560;
        alphaiw[9] = 9.1984520;
        alphaiw[10] = 7.9445240;
        alphaiw[11] = 6.9570500;
        alphaiw[12] = 5.5167450;
        alphaiw[13] = 4.5271540;
        alphaiw[14] = 3.8097850;
        alphaiw[15] = 3.2672340;
        alphaiw[16] = 2.8429600;
        alphaiw[17] = 2.1014180;
        alphaiw[18] = 1.6245240;
        alphaiw[19] = 1.0563850;
        alphaiw[20] = 0.7399150;
        alphaiw[21] = 0.3695680;
        alphaiw[22] = 0.2190670;
      }
      else if (refsys == 2) {
        alphaiw[0] = 47.9831520;
        alphaiw[1] = 46.8085170;
        alphaiw[2] = 43.7444240;
        alphaiw[3] = 35.4613210;
        alphaiw[4] = 27.8088110;
        alphaiw[5] = 21.9082420;
        alphaiw[6] = 17.5604940;
        alphaiw[7] = 14.3629380;
        alphaiw[8] = 11.9797160;
        alphaiw[9] = 10.1702660;
        alphaiw[10] = 8.7691240;
        alphaiw[11] = 7.6631150;
        alphaiw[12] = 6.0479100;
        alphaiw[13] = 4.9390680;
        alphaiw[14] = 4.1374920;
        alphaiw[15] = 3.5335950;
        alphaiw[16] = 3.0634060;
        alphaiw[17] = 2.2476620;
        alphaiw[18] = 1.7283170;
        alphaiw[19] = 1.1162240;
        alphaiw[20] = 0.7787240;
        alphaiw[21] = 0.3870890;
        alphaiw[22] = 0.2289900;
      }
      else if (refsys == 3) {
        alphaiw[0] = 49.9757880;
        alphaiw[1] = 48.8536140;
        alphaiw[2] = 45.9097740;
        alphaiw[3] = 37.7785020;
        alphaiw[4] = 29.9993350;
        alphaiw[5] = 23.8231110;
        alphaiw[6] = 19.1794330;
        alphaiw[7] = 15.7185060;
        alphaiw[8] = 13.1162960;
        alphaiw[9] = 11.1291240;
        alphaiw[10] = 9.5845740;
        alphaiw[11] = 8.3625420;
        alphaiw[12] = 6.5753180;
        alphaiw[13] = 5.3486700;
        alphaiw[14] = 4.4636720;
        alphaiw[15] = 3.7988870;
        alphaiw[16] = 3.2830670;
        alphaiw[17] = 2.3934820;
        alphaiw[18] = 1.8318420;
        alphaiw[19] = 1.1759270;
        alphaiw[20] = 0.8174500;
        alphaiw[21] = 0.4045760;
        alphaiw[22] = 0.2388960;
      }
      break;
    case 85:
      if (refsys == 1) {
        alphaiw[0] = 39.1638670;
        alphaiw[1] = 38.4250150;
        alphaiw[2] = 36.4213580;
        alphaiw[3] = 30.5411780;
        alphaiw[4] = 24.6134160;
        alphaiw[5] = 19.7657890;
        alphaiw[6] = 16.0524090;
        alphaiw[7] = 13.2471470;
        alphaiw[8] = 11.1156700;
        alphaiw[9] = 9.4745710;
        alphaiw[10] = 8.1909470;
        alphaiw[11] = 7.1705270;
        alphaiw[12] = 5.6711490;
        alphaiw[13] = 4.6378720;
        alphaiw[14] = 3.8904420;
        alphaiw[15] = 3.3277200;
        alphaiw[16] = 2.8899970;
        alphaiw[17] = 2.1310810;
        alphaiw[18] = 1.6471780;
        alphaiw[19] = 1.0734420;
        alphaiw[20] = 0.7539160;
        alphaiw[21] = 0.3784910;
        alphaiw[22] = 0.2250150;
      }
      else if (refsys == 2) {
        alphaiw[0] = 42.3026280;
        alphaiw[1] = 41.5366110;
        alphaiw[2] = 39.4596880;
        alphaiw[3] = 33.3254640;
        alphaiw[4] = 27.0380510;
        alphaiw[5] = 21.8068890;
        alphaiw[6] = 17.7463530;
        alphaiw[7] = 14.6509070;
        alphaiw[8] = 12.2850030;
        alphaiw[9] = 10.4566080;
        alphaiw[10] = 9.0232990;
        alphaiw[11] = 7.8825660;
        alphaiw[12] = 6.2059600;
        alphaiw[13] = 5.0520660;
        alphaiw[14] = 4.2196050;
        alphaiw[15] = 3.5950280;
        alphaiw[16] = 3.1110670;
        alphaiw[17] = 2.2775460;
        alphaiw[18] = 1.7510350;
        alphaiw[19] = 1.1332590;
        alphaiw[20] = 0.7926970;
        alphaiw[21] = 0.3959980;
        alphaiw[22] = 0.2349320;
      }
      break;
    case 86:
      if (refsys == 1) {
        alphaiw[0] = 33.6514050;
        alphaiw[1] = 33.1852420;
        alphaiw[2] = 31.8852830;
        alphaiw[3] = 27.7786920;
        alphaiw[4] = 23.2154680;
        alphaiw[5] = 19.1736050;
        alphaiw[6] = 15.8877090;
        alphaiw[7] = 13.2945880;
        alphaiw[8] = 11.2595530;
        alphaiw[9] = 9.6544200;
        alphaiw[10] = 8.3760900;
        alphaiw[11] = 7.3462510;
        alphaiw[12] = 5.8137240;
        alphaiw[13] = 4.7476360;
        alphaiw[14] = 3.9742240;
        alphaiw[15] = 3.3924070;
        alphaiw[16] = 2.9410550;
        alphaiw[17] = 2.1630540;
        alphaiw[18] = 1.6707900;
        alphaiw[19] = 1.0904110;
        alphaiw[20] = 0.7676850;
        alphaiw[21] = 0.3873300;
        alphaiw[22] = 0.2309680;
      }
      break;
    case 112:
      if (refsys == 1) {
        alphaiw[0] = 28.9658506;
        alphaiw[1] = 28.5229043;
        alphaiw[2] = 27.3194441;
        alphaiw[3] = 23.7497041;
        alphaiw[4] = 20.0600679;
        alphaiw[5] = 16.9382662;
        alphaiw[6] = 14.4456995;
        alphaiw[7] = 12.4709213;
        alphaiw[8] = 10.8915816;
        alphaiw[9] = 9.6105637;
        alphaiw[10] = 8.5564757;
        alphaiw[11] = 7.6775158;
        alphaiw[12] = 6.3032570;
        alphaiw[13] = 5.2860251;
        alphaiw[14] = 4.5087413;
        alphaiw[15] = 3.8995856;
        alphaiw[16] = 3.4123111;
        alphaiw[17] = 2.5460363;
        alphaiw[18] = 1.9872976;
        alphaiw[19] = 1.3293939;
        alphaiw[20] = 0.9663593;
        alphaiw[21] = 0.5321666;
        alphaiw[22] = 0.3423450;
      }
      break;
    case 113:
      if (refsys == 1) {
        alphaiw[0] = 36.8515202;
        alphaiw[1] = 35.9314075;
        alphaiw[2] = 33.7296086;
        alphaiw[3] = 28.2311084;
        alphaiw[4] = 23.2301521;
        alphaiw[5] = 19.2714758;
        alphaiw[6] = 16.2269059;
        alphaiw[7] = 13.8738534;
        alphaiw[8] = 12.0270465;
        alphaiw[9] = 10.5514261;
        alphaiw[10] = 9.3517015;
        alphaiw[11] = 8.3607091;
        alphaiw[12] = 6.8266410;
        alphaiw[13] = 5.7014944;
        alphaiw[14] = 4.8468778;
        alphaiw[15] = 4.1800460;
        alphaiw[16] = 3.6485283;
        alphaiw[17] = 2.7081324;
        alphaiw[18] = 2.1050561;
        alphaiw[19] = 1.3991594;
        alphaiw[20] = 1.0122747;
        alphaiw[21] = 0.5535296;
        alphaiw[22] = 0.3548537;
      }
      break;
    case 114:
      if (refsys == 1) {
        alphaiw[0] = 31.0177261;
        alphaiw[1] = 30.3409187;
        alphaiw[2] = 28.5968922;
        alphaiw[3] = 24.0264841;
        alphaiw[4] = 19.8927260;
        alphaiw[5] = 16.6470264;
        alphaiw[6] = 14.1460801;
        alphaiw[7] = 12.2005462;
        alphaiw[8] = 10.6619593;
        alphaiw[9] = 9.4233104;
        alphaiw[10] = 8.4088155;
        alphaiw[11] = 7.5647344;
        alphaiw[12] = 6.2444458;
        alphaiw[13] = 5.2627485;
        alphaiw[14] = 4.5076687;
        alphaiw[15] = 3.9118848;
        alphaiw[16] = 3.4323326;
        alphaiw[17] = 2.5724899;
        alphaiw[18] = 2.0125562;
        alphaiw[19] = 1.3479270;
        alphaiw[20] = 0.9794339;
        alphaiw[21] = 0.5389514;
        alphaiw[22] = 0.3469517;
      }
      break;
    case 116:
      if (refsys == 1) {
        alphaiw[0] = 62.5965347;
        alphaiw[1] = 59.9793469;
        alphaiw[2] = 54.2098284;
        alphaiw[3] = 41.6476762;
        alphaiw[4] = 31.8504866;
        alphaiw[5] = 24.9639388;
        alphaiw[6] = 20.1291251;
        alphaiw[7] = 16.6500840;
        alphaiw[8] = 14.0715476;
        alphaiw[9] = 12.1054557;
        alphaiw[10] = 10.5676443;
        alphaiw[11] = 9.3376437;
        alphaiw[12] = 7.5027755;
        alphaiw[13] = 6.2067008;
        alphaiw[14] = 5.2468060;
        alphaiw[15] = 4.5100150;
        alphaiw[16] = 3.9287041;
        alphaiw[17] = 2.9084035;
        alphaiw[18] = 2.2561760;
        alphaiw[19] = 1.4921393;
        alphaiw[20] = 1.0737727;
        alphaiw[21] = 0.5817618;
        alphaiw[22] = 0.3715195;
      }
      break;
    case 115:
      if (refsys == 1) {
        alphaiw[0] = 50.5879876;
        alphaiw[1] = 48.5224886;
        alphaiw[2] = 43.8104441;
        alphaiw[3] = 33.6795593;
        alphaiw[4] = 26.0680976;
        alphaiw[5] = 20.7866748;
        alphaiw[6] = 17.0594839;
        alphaiw[7] = 14.3406936;
        alphaiw[8] = 12.2916388;
        alphaiw[9] = 10.7020689;
        alphaiw[10] = 9.4379170;
        alphaiw[11] = 8.4110146;
        alphaiw[12] = 6.8482249;
        alphaiw[13] = 5.7186423;
        alphaiw[14] = 4.8669543;
        alphaiw[15] = 4.2042872;
        alphaiw[16] = 3.6761256;
        alphaiw[17] = 2.7387147;
        alphaiw[18] = 2.1335072;
        alphaiw[19] = 1.4195987;
        alphaiw[20] = 1.0262708;
        alphaiw[21] = 0.5601697;
        alphaiw[22] = 0.3590984;
      }
      break;
    case 117:
      if (refsys == 1) {
        alphaiw[0] = 68.0864127;
        alphaiw[1] = 64.1264282;
        alphaiw[2] = 56.6470176;
        alphaiw[3] = 42.0637390;
        alphaiw[4] = 31.4475652;
        alphaiw[5] = 24.2894637;
        alphaiw[6] = 19.3977301;
        alphaiw[7] = 15.9410604;
        alphaiw[8] = 13.4120039;
        alphaiw[9] = 11.5024670;
        alphaiw[10] = 10.0205855;
        alphaiw[11] = 8.8430108;
        alphaiw[12] = 7.0995079;
        alphaiw[13] = 5.8772976;
        alphaiw[14] = 4.9763464;
        alphaiw[15] = 4.2864060;
        alphaiw[16] = 3.7423450;
        alphaiw[17] = 2.7856521;
        alphaiw[18] = 2.1708261;
        alphaiw[19] = 1.4447260;
        alphaiw[20] = 1.0436089;
        alphaiw[21] = 0.5684047;
        alphaiw[22] = 0.3641795;
      }
      break;
    case 118:
      if (refsys == 1) {
        alphaiw[0] = 60.1915989;
        alphaiw[1] = 57.8350021;
        alphaiw[2] = 52.2126738;
        alphaiw[3] = 39.5649630;
        alphaiw[4] = 29.8330516;
        alphaiw[5] = 23.1060271;
        alphaiw[6] = 18.4549609;
        alphaiw[7] = 15.1539330;
        alphaiw[8] = 12.7366393;
        alphaiw[9] = 10.9128902;
        alphaiw[10] = 9.4998323;
        alphaiw[11] = 8.3792180;
        alphaiw[12] = 6.7254432;
        alphaiw[13] = 5.5711554;
        alphaiw[14] = 4.7232296;
        alphaiw[15] = 4.0753582;
        alphaiw[16] = 3.5649854;
        alphaiw[17] = 2.6668767;
        alphaiw[18] = 2.0874183;
        alphaiw[19] = 1.3980570;
        alphaiw[20] = 1.0139007;
        alphaiw[21] = 0.5552450;
        alphaiw[22] = 0.3569580;
      }
      break;
  }
}
void Hcount(int atmnr, std::vector<double> & hcount) {
  //hcount values for D4
  //vector not zeroed because we use RefN to access the vector's formal size
  switch (atmnr) {
    case 1:
      hcount[0] = 0.0;
      hcount[1] = 0.0;
      break;
    case 2:
      hcount[0] = 0.0;
      break;
    case 3:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 4.0;
      break;
    case 4:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      hcount[3] = 4.0;
      break;
    case 5:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      hcount[3] = 3.0;
      hcount[4] = 6.0;
      break;
    case 6:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      hcount[3] = 4.0;
      hcount[4] = 6.0;
      hcount[5] = 6.0;
      hcount[6] = 1.0;
      break;
    case 7:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      hcount[3] = 3.0;
      hcount[4] = 0.0;
      break;
    case 8:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      hcount[3] = 0.0;
      break;
    case 9:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      break;
    case 10:
      hcount[0] = 0.0;
      break;
    case 11:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 4.0;
      break;
    case 12:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      hcount[3] = 4.0;
      break;
    case 13:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      hcount[3] = 3.0;
      break;
    case 14:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      hcount[3] = 4.0;
      hcount[4] = 6.0;
      break;
    case 15:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      hcount[3] = 3.0;
      break;
    case 16:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      break;
    case 17:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      break;
    case 18:
      hcount[0] = 0.0;
      break;
    case 19:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 4.0;
      break;
    case 20:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      hcount[3] = 4.0;
      break;
    case 21:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 3.0;
      hcount[3] = 1.5;
      break;
    case 22:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 4.0;
      hcount[3] = 4.0;
      break;
    case 23:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 3.0;
      hcount[3] = 5.0;
      break;
    case 24:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 12.0;
      hcount[3] = 6.0;
      break;
    case 25:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 11.0;
      break;
    case 26:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 10.0;
      break;
    case 27:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 3.0;
      hcount[3] = 9.0;
      break;
    case 28:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 8.0;
      hcount[3] = 8.0;
      break;
    case 29:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      break;
    case 30:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      break;
    case 31:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      break;
    case 32:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      hcount[3] = 3.0;
      hcount[4] = 4.0;
      break;
    case 33:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      hcount[3] = 3.0;
      break;
    case 34:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      break;
    case 35:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      break;
    case 36:
      hcount[0] = 0.0;
      break;
    case 37:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 4.0;
      break;
    case 38:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      hcount[3] = 4.0;
      break;
    case 39:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 3.0;
      break;
    case 40:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 4.0;
      hcount[3] = 4.0;
      break;
    case 41:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 3.0;
      hcount[3] = 5.0;
      break;
    case 42:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 12.0;
      hcount[3] = 6.0;
      break;
    case 43:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 11.0;
      break;
    case 44:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 10.0;
      break;
    case 45:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 3.0;
      hcount[3] = 5.0;
      break;
    case 46:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 8.0;
      break;
    case 47:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      break;
    case 48:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      break;
    case 49:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      hcount[3] = 3.0;
      break;
    case 50:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      hcount[3] = 3.0;
      hcount[4] = 4.0;
      break;
    case 51:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      hcount[3] = 3.0;
      break;
    case 52:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      break;
    case 53:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      break;
    case 54:
      hcount[0] = 0.0;
      break;
    case 55:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 4.0;
      break;
    case 56:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      hcount[3] = 4.0;
      break;
    case 57:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 3.0;
      break;
    case 58:
      hcount[0] = 3.0;
      break;
    case 59:
      hcount[0] = 0.0;
      hcount[1] = 3.0;
      break;
    case 60:
      hcount[0] = 0.0;
      hcount[1] = 3.0;
      break;
    case 61:
      hcount[0] = 0.0;
      hcount[1] = 3.0;
      break;
    case 62:
      hcount[0] = 0.0;
      hcount[1] = 3.0;
      break;
    case 63:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      break;
    case 64:
      hcount[0] = 0.0;
      hcount[1] = 3.0;
      break;
    case 65:
      hcount[0] = 0.0;
      hcount[1] = 3.0;
      break;
    case 66:
      hcount[0] = 0.0;
      hcount[1] = 3.0;
      break;
    case 67:
      hcount[0] = 0.0;
      hcount[1] = 3.0;
      break;
    case 68:
      hcount[0] = 0.0;
      hcount[1] = 3.0;
      break;
    case 69:
      hcount[0] = 0.0;
      hcount[1] = 3.0;
      break;
    case 70:
      hcount[0] = 0.0;
      hcount[1] = 3.0;
      break;
    case 71:
      hcount[0] = 0.0;
      hcount[1] = 3.0;
      break;
    case 72:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 4.0;
      hcount[3] = 4.0;
      break;
    case 73:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 3.0;
      hcount[3] = 5.0;
      break;
    case 74:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 12.0;
      break;
    case 75:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 11.0;
      break;
    case 76:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 10.0;
      break;
    case 77:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 3.0;
      hcount[3] = 4.0;
      hcount[4] = 5.0;
      break;
    case 78:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      hcount[2] = 8.0;
      break;
    case 79:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      break;
    case 80:
      hcount[0] = 0.0;
      hcount[1] = 2.0;
      break;
    case 81:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      hcount[3] = 3.0;
      break;
    case 82:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      hcount[3] = 3.0;
      hcount[4] = 4.0;
      break;
    case 83:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      hcount[3] = 3.0;
      break;
    case 84:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      hcount[2] = 2.0;
      break;
    case 85:
      hcount[0] = 0.0;
      hcount[1] = 1.0;
      break;
    case 86:
      hcount[0] = 0.0;
      break;
    case 112:
      hcount[0] = 0.0;
      break;
    case 113:
      hcount[0] = 1.0;
      break;
    case 114:
      hcount[0] = 0.0;
      break;
    case 116:
      hcount[0] = 2.0;
      break;
    case 115:
      hcount[0] = 1.0;
      break;
    case 117:
      hcount[0] = 1.0;
      break;
    case 118:
      hcount[0] = 0.0;
      break;
  }
}
double AScale(int atmnr, int refsys) {
  //the ascale parameter
  double ascal = 1.0;
  if ((atmnr == 1)&&(refsys == 2)) {ascal = 0.5;}
  else if ((atmnr == 6)&&(refsys == 6)) {ascal = 0.16666666666667;}
  else if (((atmnr == 3)||(atmnr == 11)||(atmnr == 19)||(atmnr == 37)||(atmnr == 55))&&(refsys == 3)) {ascal = 0.25;}
  else if (((atmnr == 4)||(atmnr == 12)||(atmnr == 20)||(atmnr == 38)||(atmnr == 56))&&(refsys == 4)) {ascal = 0.25;}
  else if (((atmnr == 6)||(atmnr == 7)||(atmnr == 15))&&(refsys == 3)) {ascal = 0.5;}
  else if (((atmnr == 6)||(atmnr == 8)||(atmnr == 14))&&(refsys == 4)) {ascal = 0.5;}
  else if (((atmnr == 5)||(atmnr == 6)||(atmnr == 7)||(atmnr == 14))&&(refsys == 5)) {ascal = 0.5;}
  return ascal;
}
void RefCN(int atmnr, std::vector<double> & refcn) {
  //reference covalent CN values for D4
  //vector not zeroed because we use RefN to access the vector's formal size
  refcn[0] = 0.00000000000000;
  switch (atmnr) {
    case 1:
      refcn[1] = 0.91179218222140;
      break;
    case 2:
      break;
    case 3:
      refcn[1] = 0.98647651737145;
      refcn[2] = 6.93900000000000;
      break;
    case 4:
      refcn[1] = 0.98079323297698;
      refcn[2] = 1.96965813264154;
      refcn[3] = 10.26200000000000;
      break;
    case 5:
      refcn[1] = 0.97061615724771;
      refcn[2] = 1.94412533738189;
      refcn[3] = 2.91284702438220;
      refcn[4] = 4.58564454717390;
      break;
    case 6:
      refcn[1] = 0.98678926403687;
      refcn[2] = 1.99847658843687;
      refcn[3] = 2.99872598643336;
      refcn[4] = 3.98439894415100;
      refcn[5] = 3.14177128194115;
      refcn[6] = 0.99995005170393;
      break;
    case 7:
      refcn[1] = 0.99440725902444;
      refcn[2] = 2.01427212990570;
      refcn[3] = 2.99026427338514;
      refcn[4] = 0.99998952959314;
      break;
    case 8:
      refcn[1] = 0.99245941779060;
      refcn[2] = 1.98869274482158;
      refcn[3] = 0.99856725131399;
      break;
    case 9:
      refcn[1] = 0.99820363068818;
      break;
    case 10:
      break;
    case 11:
      refcn[1] = 0.96835675605614;
      refcn[2] = 8.18700000000000;
      break;
    case 12:
      refcn[1] = 0.96281367993504;
      refcn[2] = 1.94960242607378;
      refcn[3] = 16.20300000000000;
      break;
    case 13:
      refcn[1] = 0.96483397822083;
      refcn[2] = 1.93105228668158;
      refcn[3] = 2.91455515763177;
      break;
    case 14:
      refcn[1] = 0.95070518838146;
      refcn[2] = 1.94352476341216;
      refcn[3] = 2.94071103313097;
      refcn[4] = 3.86773979105472;
      break;
    case 15:
      refcn[1] = 0.99467058048882;
      refcn[2] = 2.01017301871106;
      refcn[3] = 2.98587564194422;
      break;
    case 16:
      refcn[1] = 0.99484136522295;
      refcn[2] = 1.99034385492992;
      break;
    case 17:
      refcn[1] = 0.99718907461633;
      break;
    case 18:
      break;
    case 19:
      refcn[1] = 0.97673829344237;
      refcn[2] = 13.53400000000000;
      break;
    case 20:
      refcn[1] = 0.98305940539640;
      refcn[2] = 1.93490105201102;
      refcn[3] = 18.07900000000000;
      break;
    case 21:
      refcn[1] = 1.86271675054418;
      refcn[2] = 2.89994293581383;
      refcn[3] = 6.69400000000000;
      break;
    case 22:
      refcn[1] = 1.82992898343369;
      refcn[2] = 3.86754432133472;
      refcn[3] = 6.63750000000000;
      break;
    case 23:
      refcn[1] = 1.91379548691600;
      refcn[2] = 2.91102328524604;
      refcn[3] = 5.15100000000000;
      break;
    case 24:
      refcn[1] = 1.82687247869787;
      refcn[2] = 10.61913477887437;
      refcn[3] = 6.07975116715200;
      break;
    case 25:
      refcn[1] = 1.64060352754295;
      refcn[2] = 9.88490304290262;
      break;
    case 26:
      refcn[1] = 1.64827170778593;
      refcn[2] = 9.13759391858428;
      break;
    case 27:
      refcn[1] = 1.71485433521559;
      refcn[2] = 2.92629831246542;
      refcn[3] = 7.77854420648529;
      break;
    case 28:
      refcn[1] = 1.79373527663037;
      refcn[2] = 6.54584816193237;
      refcn[3] = 6.29178843690998;
      break;
    case 29:
      refcn[1] = 0.95762492307496;
      break;
    case 30:
      refcn[1] = 1.94189578186459;
      break;
    case 31:
      refcn[1] = 0.96005556460526;
      refcn[2] = 1.93153571165990;
      break;
    case 32:
      refcn[1] = 0.94335166609313;
      refcn[2] = 1.94465170362984;
      refcn[3] = 2.91863046426809;
      refcn[4] = 3.89724939850564;
      break;
    case 33:
      refcn[1] = 0.98892125897221;
      refcn[2] = 1.97932099000612;
      refcn[3] = 2.97091083178345;
      break;
    case 34:
      refcn[1] = 0.99005644237283;
      refcn[2] = 1.98124719453101;
      break;
    case 35:
      refcn[1] = 0.99735391054501;
      break;
    case 36:
      break;
    case 37:
      refcn[1] = 0.97383869672670;
      refcn[2] = 14.53800000000000;
      break;
    case 38:
      refcn[1] = 0.98009983712321;
      refcn[2] = 1.91434965086468;
      refcn[3] = 18.12500000000000;
      break;
    case 39:
      refcn[1] = 1.91534612543569;
      refcn[2] = 2.89028721234797;
      break;
    case 40:
      refcn[1] = 1.93553367207877;
      refcn[2] = 3.91060485457820;
      refcn[3] = 9.43900000000000;
      break;
    case 41:
      refcn[1] = 1.95448265742483;
      refcn[2] = 2.92252286735086;
      refcn[3] = 6.19100000000000;
      break;
    case 42:
      refcn[1] = 1.94202903408340;
      refcn[2] = 11.05564046968760;
      refcn[3] = 5.73689802594352;
      break;
    case 43:
      refcn[1] = 1.66816908014472;
      refcn[2] = 9.54020187746512;
      break;
    case 44:
      refcn[1] = 1.85843916069775;
      refcn[2] = 8.88947582399905;
      break;
    case 45:
      refcn[1] = 1.90032831305941;
      refcn[2] = 2.96957012334100;
      refcn[3] = 4.84416063571149;
      break;
    case 46:
      refcn[1] = 1.86304976118092;
      refcn[2] = 5.70949602385229;
      break;
    case 47:
      refcn[1] = 0.96788403280472;
      break;
    case 48:
      refcn[1] = 1.95389378202525;
      break;
    case 49:
      refcn[1] = 0.96330831417334;
      refcn[2] = 1.93784277625627;
      refcn[3] = 2.93528671852710;
      break;
    case 50:
      refcn[1] = 0.95144638754491;
      refcn[2] = 1.95049412276785;
      refcn[3] = 2.92590527066264;
      refcn[4] = 3.91232527967725;
      break;
    case 51:
      refcn[1] = 0.97493375701249;
      refcn[2] = 1.95228160320121;
      refcn[3] = 2.93146774469466;
      break;
    case 52:
      refcn[1] = 0.98114838477859;
      refcn[2] = 1.96387521222899;
      break;
    case 53:
      refcn[1] = 0.99682171802295;
      break;
    case 54:
      break;
    case 55:
      refcn[1] = 0.99085899461464;
      refcn[2] = 16.33600000000000;
      break;
    case 56:
      refcn[1] = 0.97966156728304;
      refcn[2] = 1.84668170996467;
      refcn[3] = 17.41200000000000;
      break;
    case 57:
      refcn[1] = 1.93733589109731;
      refcn[2] = 2.91753935596640;
      break;
    case 58:
      refcn[0] = 2.79905712137645;
      break;
    case 59:
      refcn[1] = 2.94245250980475;
      break;
    case 60:
      refcn[1] = 2.94550013327159;
      break;
    case 61:
      refcn[1] = 2.94127110730758;
      break;
    case 62:
      refcn[1] = 2.92998749109076;
      break;
    case 63:
      refcn[1] = 1.82859934010382;
      break;
    case 64:
      refcn[1] = 2.87315743546354;
      break;
    case 65:
      refcn[1] = 2.90861742800663;
      break;
    case 66:
      refcn[1] = 2.89653794321150;
      break;
    case 67:
      refcn[1] = 2.92423628373451;
      break;
    case 68:
      refcn[1] = 2.92818783679168;
      break;
    case 69:
      refcn[1] = 2.92455515909581;
      break;
    case 70:
      refcn[1] = 2.84817486472051;
      break;
    case 71:
      refcn[1] = 2.92190497244653;
      break;
    case 72:
      refcn[1] = 1.92539163150370;
      refcn[2] = 3.88399563447207;
      refcn[3] = 6.33700000000000;
      break;
    case 73:
      refcn[1] = 1.94592039225300;
      refcn[2] = 2.89880253611560;
      refcn[3] = 6.20500000000000;
      break;
    case 74:
      refcn[1] = 1.92915611172061;
      refcn[2] = 10.91528020745492;
      break;
    case 75:
      refcn[1] = 1.81037867017969;
      refcn[2] = 9.80535616498394;
      break;
    case 76:
      refcn[1] = 1.88576348924552;
      refcn[2] = 9.15269257898964;
      break;
    case 77:
      refcn[1] = 1.86479126833318;
      refcn[2] = 2.94241057251882;
      refcn[3] = 3.85906061728294;
      refcn[4] = 4.74946213886664;
      break;
    case 78:
      refcn[1] = 1.91881173481654;
      refcn[2] = 6.66687770578862;
      break;
    case 79:
      refcn[1] = 0.98458933772055;
      break;
    case 80:
      refcn[1] = 1.98963020776929;
      break;
    case 81:
      refcn[1] = 0.92674492705574;
      refcn[2] = 1.93019298965987;
      refcn[3] = 2.94201332191403;
      break;
    case 82:
      refcn[1] = 0.93830467603587;
      refcn[2] = 1.93560176827475;
      refcn[3] = 2.90809757363721;
      refcn[4] = 3.90982050474815;
      break;
    case 83:
      refcn[1] = 0.98198108329616;
      refcn[2] = 1.96550673996239;
      refcn[3] = 2.94995150877077;
      break;
    case 84:
      refcn[1] = 0.98152738589119;
      refcn[2] = 1.96393490432407;
      break;
    case 85:
      refcn[1] = 0.99539313449513;
      break;
    case 86:
      break;
    case 112:
      break;
    case 113:
      refcn[0] = 0.94277714762458;
      break;
    case 114:
      break;
    case 116:
      refcn[0] = 1.98131399263162;
      break;
    case 115:
      refcn[0] = 0.91220333076895;
      break;
    case 117:
      refcn[0] = 0.97277918745352;
      break;
    case 118:
      break;
  }
}
void getNcount(int atmnr, std::vector<int> & nc) {
  switch (atmnr) {
    case 1:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 2:
      nc[0] = 3;
      break;
    case 3:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      break;
    case 4:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 5:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      nc[4] = 1;
      break;
    case 6:
      nc[0] = 3;
      nc[1] = 3;
      nc[2] = 1;
      nc[3] = 3;
      nc[4] = 1;
      nc[5] = 3;
      nc[6] = 3;
      break;
    case 7:
      nc[0] = 3;
      nc[1] = 3;
      nc[2] = 1;
      nc[3] = 1;
      nc[4] = 3;
      break;
    case 8:
      nc[0] = 3;
      nc[1] = 3;
      nc[2] = 1;
      nc[3] = 3;
      break;
    case 9:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 10:
      nc[0] = 3;
      break;
    case 11:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      break;
    case 12:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 13:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 14:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      nc[4] = 1;
      break;
    case 15:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 16:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      break;
    case 17:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 18:
      nc[0] = 3;
      break;
    case 19:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      break;
    case 20:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 21:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 22:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 23:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 24:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 25:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      break;
    case 26:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      break;
    case 27:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 28:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 29:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 30:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 31:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      break;
    case 32:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      nc[4] = 1;
      break;
    case 33:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 34:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      break;
    case 35:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 36:
      nc[0] = 3;
      break;
    case 37:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      break;
    case 38:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 39:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      break;
    case 40:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 41:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 42:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 43:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      break;
    case 44:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      break;
    case 45:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 46:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      break;
    case 47:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 48:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 49:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 50:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      nc[4] = 1;
      break;
    case 51:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 52:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      break;
    case 53:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 54:
      nc[0] = 3;
      break;
    case 55:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      break;
    case 56:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 57:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      break;
    case 58:
      nc[0] = 1;
      break;
    case 59:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 60:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 61:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 62:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 63:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 64:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 65:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 66:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 67:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 68:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 69:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 70:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 71:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 72:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 73:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 74:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      break;
    case 75:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      break;
    case 76:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      break;
    case 77:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      nc[4] = 1;
      break;
    case 78:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      break;
    case 79:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 80:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 81:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 82:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      nc[4] = 1;
      break;
    case 83:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      nc[3] = 1;
      break;
    case 84:
      nc[0] = 3;
      nc[1] = 1;
      nc[2] = 1;
      break;
    case 85:
      nc[0] = 3;
      nc[1] = 1;
      break;
    case 86:
      nc[0] = 3;
      break;
    case 87:
      nc[0] = 0;
      break;
    case 88:
      nc[0] = 0;
      break;
    case 89:
      nc[0] = 0;
      break;
    case 90:
      nc[0] = 0;
      break;
    case 91:
      nc[0] = 0;
      break;
    case 92:
      nc[0] = 0;
      break;
    case 93:
      nc[0] = 0;
      break;
    case 94:
      nc[0] = 0;
      break;
    case 95:
      nc[0] = 0;
      break;
    case 96:
      nc[0] = 0;
      break;
    case 97:
      nc[0] = 0;
      break;
    case 98:
      nc[0] = 0;
      break;
    case 99:
      nc[0] = 0;
      break;
    case 100:
      nc[0] = 0;
      break;
    case 101:
      nc[0] = 0;
      break;
    case 102:
      nc[0] = 0;
      break;
    case 103:
      nc[0] = 0;
      break;
    case 104:
      nc[0] = 0;
      break;
    case 105:
      nc[0] = 0;
      break;
    case 106:
      nc[0] = 0;
      break;
    case 107:
      nc[0] = 0;
      break;
    case 108:
      nc[0] = 0;
      break;
    case 109:
      nc[0] = 0;
      break;
    case 110:
      nc[0] = 0;
      break;
    case 111:
      nc[0] = 0;
      break;
    case 112:
      nc[0] = 3;
      break;
    case 113:
      nc[0] = 1;
      break;
    case 114:
      nc[0] = 3;
      break;
    case 115:
      nc[0] = 1;
      break;
    case 116:
      nc[0] = 1;
      break;
    case 117:
      nc[0] = 1;
      break;
    case 118:
      nc[0] = 3;
      break;
  }
}
double WF() {return 6.0;}

#endif //_Dispersion4_Parameters_
