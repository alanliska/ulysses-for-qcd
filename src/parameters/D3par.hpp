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

#ifndef _Dispersion3_Parameters_
#define _Dispersion3_Parameters_
#include "D3C6par.hpp"

//description:
//parameters for D3 correction term
//some functions are recycled to be used in D4 as well

double s6D3(const std::string & method) {
  //function returning the scaling factor s6
  double ssix = 1.0;
  //D3H+
  if ((method == "PM6D3H+")||(method == "pm6D3H+")) {ssix = 1.0;}
  //D3H4X
  else if ((method == "PM6D3H4X")||(method == "pm6D3H4X")) {ssix = 0.88;}
  else if ((method == "MNDOD3H4X")||(method == "mndoD3H4X")) {ssix = 0.0;}
  else if ((method == "AM1D3H4X")||(method == "am1D3H4X")) {ssix = 1.0;}
  else if ((method == "RM1D3H4X")||(method == "rm1D3H4X")) {ssix = 1.0;}
  else if ((method == "PM3D3H4X")||(method == "pm3D3H4X")) {ssix = 1.0;}
  return ssix;
}
double s8D3(const std::string & method) {
  //function returning the scaling factor s8
  double seight = 1.0;
  //D3H+
  if ((method == "PM6D3H+")||(method == "pm6D3H+")) {seight = 1.009;}
  //D3H4X
  else if ((method == "PM6D3H4X")||(method == "pm6D3H4X")) {seight = 0.0;}
  else if ((method == "MNDOD3H4X")||(method == "mndoD3H4X")) {seight = 0.0;}
  else if ((method == "AM1D3H4X")||(method == "am1D3H4X")) {seight = 0.0;}
  else if ((method == "RM1D3H4X")||(method == "rm1D3H4X")) {seight = 0.0;}
  else if ((method == "PM3D3H4X")||(method == "pm3D3H4X")) {seight = 0.0;}
  return seight;
}
double s9D3(const std::string & method) {
  //function returning the scaling factor s9
  double snine = 0.0;
  return snine;
}
double s10D3(const std::string & method) {
  //function returning the scaling factor s10
  double sten = 0.0;
  return sten;
}
double sR(const std::string & method, int order) {
  //function returning the scaling factor sR
  double sr = 0.0;
  switch(order) {
    case 6:
      //SCC-DFTB
      if ((method == "SCC-DFTB")||(method == "scc-dftb")) {sr = 1.215;}
      //D3H+
      else if ((method == "PM6D3H+")||(method == "pm6D3H+")) {sr = 1.560;}
      //D3H4X
      else if ((method == "PM6D3H4X")||(method == "pm6D3H4X")) {sr = 1.18;}
      else if ((method == "RM1D3H4X")||(method == "rm1D3H4X")) {sr = 1.00;}
      else if ((method == "AM1D3H4X")||(method == "am1D3H4X")) {sr = 0.90;}
      else if ((method == "PM3D3H4X")||(method == "pm3D3H4X")) {sr = 0.90;}
      break;
    case 8:
      sr = 1.0;
      break;
    case 9:
      sr = 4.0/3.0;
      break;
  }
  return sr;
}
double AlphaD3(const std::string & method, int order) {
  //function returning the exponent alpha
  double alpha = 0.0;
  if (order == 6) {
    //SCC-DFTB
    if ((method == "SCC-DFTB")||(method == "scc-dftb")) {alpha = 30.0;}
    //D3H+
    else if ((method == "PM6D3H+")||(method == "pm6D3H+")) {alpha = 14.0;}
    //D3H4X
    else if ((method == "PM6D3H4X")||(method == "pm6D3H4X")) {alpha = 22.0;}
    else if ((method == "RM1D3H4X")||(method == "rm1D3H4X")) {alpha = 16.0;}
    else if ((method == "AM1D3H4X")||(method == "am1D3H4X")) {alpha = 15.0;}
    else if ((method == "PM3D3H4X")||(method == "pm3D3H4X")) {alpha = 22.0;}
  }
  else if (order == 8) {
    //SCC-DFTB
    if ((method == "SCC-DFTB")||(method == "scc-dftb")) {alpha = 32.0;}
    //D3H+
    else if ((method == "PM6D3H+")||(method == "pm6D3H+")) {alpha = 16.0;}
    //D3H4X
    else if ((method == "PM6D3H4X")||(method == "pm6D3H4X")) {alpha = 24.0;}
    else if ((method == "RM1D3H4X")||(method == "rm1D3H4X")) {alpha = 18.0;}
    else if ((method == "AM1D3H4X")||(method == "am1D3H4X")) {alpha = 17.0;}
    else if ((method == "PM3D3H4X")||(method == "pm3D3H4X")) {alpha = 24.0;}
  }
  else if (order == 9) {
    alpha = 16.0;
  }
  return alpha;
}
double A1D3(const std::string & method) {
  //function returning the a1 parameter in Becke-Johnson damping
  double a1 = 0.40;           //this is just a general value obtained by visually averaging what Grimme has in his program
  //a1 = 0.4289;
  return a1;
}
double A2D3(const std::string & method) {
  //function returning the a2 parameter in Becke-Johnson damping
  double a2 = 4.0;           //this is just a general value obtained by visually averaging what Grimme has in his program
  //a2 = 4.4407;
  return a2;
}
double covRadius(size_t atmnr) {
  //function returning covalent radii of atoms; values returned in Angstrom and scaled by 4/3
  //parameters taken from P. Pyykkoe, M. Atsumi, Chem. Eur. J., 15, 188, 2009
  //values for metals decreased by 10 %
  double crad = 0.0;
    switch (atmnr) {
      case 1:
        crad = 0.32;
        break;
      case 2:
        crad = 0.46;
        break;
      case 3:
        crad = 1.20;
        break;
      case 4:
        crad = 0.94;
        break;
      case 5:
        crad = 0.77;
        break;
      case 6:
        crad = 0.75;
        break;
      case 7:
        crad = 0.71;
        break;
      case 8:
        crad = 0.63;
        break;
      case 9:
        crad = 0.64;
        break;
      case 10:
        crad = 0.67;
        break;
      case 11:
        crad = 1.40;
        break;
      case 12:
        crad = 1.25;
        break;
      case 13:
        crad = 1.13;
        break;
      case 14:
        crad = 1.04;
        break;
      case 15:
        crad = 1.10;
        break;
      case 16:
        crad = 1.02;
        break;
      case 17:
        crad = 0.99;
        break;
      case 18:
        crad = 0.96;
        break;
      case 19:
        crad = 1.76;
        break;
      case 20:
        crad = 1.54;
        break;
      case 21:
        crad = 1.33;
        break;
      case 22:
        crad = 1.22;
        break;
      case 23:
        crad = 1.21;
        break;
      case 24:
        crad = 1.10;
        break;
      case 25:
        crad = 1.07;
        break;
      case 26:
        crad = 1.04;
        break;
      case 27:
        crad = 1.00;
        break;
      case 28:
        crad = 0.99;
        break;
      case 29:
        crad = 1.01;
        break;
      case 30:
        crad = 1.09;
        break;
      case 31:
        crad = 1.12;
        break;
      case 32:
        crad = 1.09;
        break;
      case 33:
        crad = 1.15;
        break;
      case 34:
        crad = 1.10;
        break;
      case 35:
        crad = 1.14;
        break;
      case 36:
        crad = 1.17;
        break;
      case 37:
        crad = 1.89;
        break;
      case 38:
        crad = 1.67;
        break;
      case 39:
        crad = 1.47;
        break;
      case 40:
        crad = 1.39;
        break;
      case 41:
        crad = 1.32;
        break;
      case 42:
        crad = 1.24;
        break;
      case 43:
        crad = 1.15;
        break;
      case 44:
        crad = 1.13;
        break;
      case 45:
        crad = 1.13;
        break;
      case 46:
        crad = 1.08;
        break;
      case 47:
        crad = 1.15;
        break;
      case 48:
        crad = 1.23;
        break;
      case 49:
        crad = 1.28;
        break;
      case 50:
        crad = 1.26;
        break;
      case 51:
        crad = 1.26;
        break;
      case 52:
        crad = 1.23;
        break;
      case 53:
        crad = 1.32;
        break;
      case 54:
        crad = 1.31;
        break;
      case 55:
        crad = 2.09;
        break;
      case 56:
        crad = 1.76;
        break;
      case 57:
        crad = 1.62;
        break;
      case 58:
        crad = 1.47;
        break;
      case 59:
        crad = 1.58;
        break;
      case 60:
        crad = 1.57;
        break;
      case 61:
        crad = 1.56;
        break;
      case 62:
        crad = 1.55;
        break;
      case 63:
        crad = 1.51;
        break;
      case 64:
        crad = 1.52;
        break;
      case 65:
        crad = 1.51;
        break;
      case 66:
        crad = 1.50;
        break;
      case 67:
        crad = 1.49;
        break;
      case 68:
        crad = 1.49;
        break;
      case 69:
        crad = 1.48;
        break;
      case 70:
        crad = 1.53;
        break;
      case 71:
        crad = 1.46;
        break;
      case 72:
        crad = 1.37;
        break;
      case 73:
        crad = 1.31;
        break;
      case 74:
        crad = 1.23;
        break;
      case 75:
        crad = 1.18;
        break;
      case 76:
        crad = 1.16;
        break;
      case 77:
        crad = 1.11;
        break;
      case 78:
        crad = 1.12;
        break;
      case 79:
        crad = 1.13;
        break;
      case 80:
        crad = 1.32;
        break;
      case 81:
        crad = 1.30;
        break;
      case 82:
        crad = 1.30;
        break;
      case 83:
        crad = 1.36;
        break;
      case 84:
        crad = 1.31;
        break;
      case 85:
        crad = 1.38;
        break;
      case 86:
        crad = 1.42;
        break;
      case 87:
        crad = 2.01;
        break;
      case 88:
        crad = 1.81;
        break;
      case 89:
        crad = 1.67;
        break;
      case 90:
        crad = 1.58;
        break;
      case 91:
        crad = 1.52;
        break;
      case 92:
        crad = 1.53;
        break;
      case 93:
        crad = 1.54;
        break;
      case 94:
        crad = 1.55;
        break;
      case 95:
        crad = 1.49;
        break;
      case 96:
        crad = 1.49;
        break;
      case 97:
        crad = 1.51;
        break;
      case 98:
        crad = 1.51;
        break;
      case 99:
        crad = 1.48;
        break;
      case 100:
        crad = 1.50;
        break;
      case 101:
        crad = 1.56;
        break;
      case 102:
        crad = 1.58;
        break;
      case 103:
        crad = 1.45;
        break;
      case 104:
        crad = 1.41;
        break;
      case 105:
        crad = 1.34;
        break;
      case 106:
        crad = 1.29;
        break;
      case 107:
        crad = 1.27;
        break;
      case 108:
        crad = 1.21;
        break;
      case 109:
        crad = 1.16;
        break;
      case 110:
        crad = 1.15;
        break;
      case 111:
        crad = 1.09;
        break;
      case 112:
        crad = 1.22;
        break;
      case 113:
        crad = 1.36;
        break;
      case 114:
        crad = 1.43;
        break;
      case 115:
        crad = 1.46;
        break;
      case 116:
        crad = 1.58;
        break;
      case 117:
        crad = 1.48;
        break;
      case 118:
        crad = 1.57;
        break;
  }
  return 4.0*crad/3.0;
}
double R0AB(size_t atm1, size_t atm2) {
  //function returning cut-off radii; values stored in Angstrom, returned in au
  //note that these values differ slightly from Grimme's, since we use different conversion factors
  double r0ab = 0.0;
  size_t atmA = std::min(atm1,atm2);
  size_t atmB = std::max(atm1,atm2);
  switch (atmB) {
    case 1:
      switch (atmA) {
        case 1:
          r0ab = 2.1823000907897949;
          break;
      }
      break;
    case 2:
      switch (atmA) {
        case 1:
          r0ab = 1.8546999692916870;
          break;
        case 2:
          r0ab = 1.7346999645233154;
          break;
      }
      break;
    case 3:
      switch (atmA) {
        case 1:
          r0ab = 2.9086000919342041;
          break;
        case 2:
          r0ab = 2.5731999874114990;
          break;
        case 3:
          r0ab = 3.4955999851226807;
          break;
      }
      break;
    case 4:
      switch (atmA) {
        case 1:
          r0ab = 2.3550000190734863;
          break;
        case 2:
          r0ab = 2.5095000267028809;
          break;
        case 3:
          r0ab = 2.9802000522613525;
          break;
        case 4:
          r0ab = 3.0982000827789307;
          break;
      }
      break;
    case 5:
      switch (atmA) {
        case 1:
          r0ab = 2.5141000747680664;
          break;
        case 2:
          r0ab = 2.3917000293731689;
          break;
        case 3:
          r0ab = 2.9976999759674072;
          break;
        case 4:
          r0ab = 2.9484000205993652;
          break;
        case 5:
          r0ab = 3.2160000801086426;
          break;
      }
      break;
    case 6:
      switch (atmA) {
        case 1:
          r0ab = 2.4491999149322510;
          break;
        case 2:
          r0ab = 2.2527000904083252;
          break;
        case 3:
          r0ab = 3.1933000087738037;
          break;
        case 4:
          r0ab = 3.0213999748229980;
          break;
        case 5:
          r0ab = 2.9530999660491943;
          break;
        case 6:
          r0ab = 2.9103000164031982;
          break;
      }
      break;
    case 7:
      switch (atmA) {
        case 1:
          r0ab = 2.3666999340057373;
          break;
        case 2:
          r0ab = 2.1328001022338867;
          break;
        case 3:
          r0ab = 2.8784000873565674;
          break;
        case 4:
          r0ab = 2.7660000324249268;
          break;
        case 5:
          r0ab = 2.7776000499725342;
          break;
        case 6:
          r0ab = 2.7063000202178955;
          break;
        case 7:
          r0ab = 2.6224999427795410;
          break;
      }
      break;
    case 8:
      switch (atmA) {
        case 1:
          r0ab = 2.1768000125885010;
          break;
        case 2:
          r0ab = 2.0625000000000000;
          break;
        case 3:
          r0ab = 2.6394999027252197;
          break;
        case 4:
          r0ab = 2.6647999286651611;
          break;
        case 5:
          r0ab = 2.6482000350952148;
          break;
        case 6:
          r0ab = 2.5697000026702881;
          break;
        case 7:
          r0ab = 2.4846000671386719;
          break;
        case 8:
          r0ab = 2.4816999435424805;
          break;
      }
      break;
    case 9:
      switch (atmA) {
        case 1:
          r0ab = 2.0645999908447266;
          break;
        case 2:
          r0ab = 1.9890999794006348;
          break;
        case 3:
          r0ab = 2.5085999965667725;
          break;
        case 4:
          r0ab = 2.6907999515533447;
          break;
        case 5:
          r0ab = 2.6233000755310059;
          break;
        case 6:
          r0ab = 2.4769999980926514;
          break;
        case 7:
          r0ab = 2.3884999752044678;
          break;
        case 8:
          r0ab = 2.3510999679565430;
          break;
        case 9:
          r0ab = 2.2995998859405518;
          break;
      }
      break;
    case 10:
      switch (atmA) {
        case 1:
          r0ab = 1.9891999959945679;
          break;
        case 2:
          r0ab = 1.9250999689102173;
          break;
        case 3:
          r0ab = 2.4189999103546143;
          break;
        case 4:
          r0ab = 2.5473001003265381;
          break;
        case 5:
          r0ab = 2.4993999004364014;
          break;
        case 6:
          r0ab = 2.4091000556945801;
          break;
        case 7:
          r0ab = 2.3176000118255615;
          break;
        case 8:
          r0ab = 2.2571001052856445;
          break;
        case 9:
          r0ab = 2.1946001052856445;
          break;
        case 10:
          r0ab = 2.1373999118804932;
          break;
      }
      break;
    case 11:
      switch (atmA) {
        case 1:
          r0ab = 2.9897999763488770;
          break;
        case 2:
          r0ab = 2.6396999359130859;
          break;
        case 3:
          r0ab = 3.6031000614166260;
          break;
        case 4:
          r0ab = 3.1219000816345215;
          break;
        case 5:
          r0ab = 3.7620000839233398;
          break;
        case 6:
          r0ab = 3.2485001087188721;
          break;
        case 7:
          r0ab = 2.9356999397277832;
          break;
        case 8:
          r0ab = 2.7093000411987305;
          break;
        case 9:
          r0ab = 2.5780999660491943;
          break;
        case 10:
          r0ab = 2.4839000701904297;
          break;
        case 11:
          r0ab = 3.7081999778747559;
          break;
      }
      break;
    case 12:
      switch (atmA) {
        case 1:
          r0ab = 2.5129001140594482;
          break;
        case 2:
          r0ab = 2.7321000099182129;
          break;
        case 3:
          r0ab = 3.1052000522613525;
          break;
        case 4:
          r0ab = 3.2962000370025635;
          break;
        case 5:
          r0ab = 3.1331000328063965;
          break;
        case 6:
          r0ab = 3.2000000476837158;
          break;
        case 7:
          r0ab = 2.9586000442504883;
          break;
        case 8:
          r0ab = 3.0822000503540039;
          break;
        case 9:
          r0ab = 2.8582000732421875;
          break;
        case 10:
          r0ab = 2.7119998931884766;
          break;
        case 11:
          r0ab = 3.2569999694824219;
          break;
        case 12:
          r0ab = 3.4839000701904297;
          break;
      }
      break;
    case 13:
      switch (atmA) {
        case 1:
          r0ab = 2.8766000270843506;
          break;
        case 2:
          r0ab = 2.7427000999450684;
          break;
        case 3:
          r0ab = 3.2776000499725342;
          break;
        case 4:
          r0ab = 3.2362999916076660;
          break;
        case 5:
          r0ab = 3.5929000377655029;
          break;
        case 6:
          r0ab = 3.2825999259948730;
          break;
        case 7:
          r0ab = 3.0910999774932861;
          break;
        case 8:
          r0ab = 2.9368999004364014;
          break;
        case 9:
          r0ab = 2.9030001163482666;
          break;
        case 10:
          r0ab = 2.7788999080657959;
          break;
        case 11:
          r0ab = 3.3921000957489014;
          break;
        case 12:
          r0ab = 3.3970000743865967;
          break;
        case 13:
          r0ab = 4.0106000900268555;
          break;
      }
      break;
    case 14:
      switch (atmA) {
        case 1:
          r0ab = 2.8884000778198242;
          break;
        case 2:
          r0ab = 2.6605000495910645;
          break;
        case 3:
          r0ab = 3.7513000965118408;
          break;
        case 4:
          r0ab = 3.1612999439239502;
          break;
        case 5:
          r0ab = 3.3605000972747803;
          break;
        case 6:
          r0ab = 3.3324999809265137;
          break;
        case 7:
          r0ab = 3.0991001129150391;
          break;
        case 8:
          r0ab = 2.9296998977661133;
          break;
        case 9:
          r0ab = 2.8673999309539795;
          break;
        case 10:
          r0ab = 2.7571001052856445;
          break;
        case 11:
          r0ab = 3.8129000663757324;
          break;
        case 12:
          r0ab = 3.3266000747680664;
          break;
        case 13:
          r0ab = 3.7105000019073486;
          break;
        case 14:
          r0ab = 3.7916998863220215;
          break;
      }
      break;
    case 15:
      switch (atmA) {
        case 1:
          r0ab = 2.8303999900817871;
          break;
        case 2:
          r0ab = 2.5538001060485840;
          break;
        case 3:
          r0ab = 3.3931999206542969;
          break;
        case 4:
          r0ab = 3.1192998886108398;
          break;
        case 5:
          r0ab = 3.1865999698638916;
          break;
        case 6:
          r0ab = 3.1245000362396240;
          break;
        case 7:
          r0ab = 3.0464999675750732;
          break;
        case 8:
          r0ab = 2.8726999759674072;
          break;
        case 9:
          r0ab = 2.7664000988006592;
          break;
        case 10:
          r0ab = 2.6926000118255615;
          break;
        case 11:
          r0ab = 3.4607999324798584;
          break;
        case 12:
          r0ab = 3.2983999252319336;
          break;
        case 13:
          r0ab = 3.5141999721527100;
          break;
        case 14:
          r0ab = 3.5418000221252441;
          break;
        case 15:
          r0ab = 3.5016999244689941;
          break;
      }
      break;
    case 16:
      switch (atmA) {
        case 1:
          r0ab = 2.6189999580383301;
          break;
        case 2:
          r0ab = 2.4797000885009766;
          break;
        case 3:
          r0ab = 3.1331000328063965;
          break;
        case 4:
          r0ab = 3.0539999008178711;
          break;
        case 5:
          r0ab = 3.0650999546051025;
          break;
        case 6:
          r0ab = 2.9879000186920166;
          break;
        case 7:
          r0ab = 2.9054000377655029;
          break;
        case 8:
          r0ab = 2.8805000782012939;
          break;
        case 9:
          r0ab = 2.7330000400543213;
          break;
        case 10:
          r0ab = 2.6331000328063965;
          break;
        case 11:
          r0ab = 3.2095999717712402;
          break;
        case 12:
          r0ab = 3.5668001174926758;
          break;
        case 13:
          r0ab = 3.3684000968933105;
          break;
        case 14:
          r0ab = 3.3685998916625977;
          break;
        case 15:
          r0ab = 3.3180000782012939;
          break;
        case 16:
          r0ab = 3.3106999397277832;
          break;
      }
      break;
    case 17:
      switch (atmA) {
        case 1:
          r0ab = 2.4756999015808105;
          break;
        case 2:
          r0ab = 2.4019000530242920;
          break;
        case 3:
          r0ab = 2.9788999557495117;
          break;
        case 4:
          r0ab = 3.1468000411987305;
          break;
        case 5:
          r0ab = 2.9767999649047852;
          break;
        case 6:
          r0ab = 2.8847999572753906;
          break;
        case 7:
          r0ab = 2.7952001094818115;
          break;
        case 8:
          r0ab = 2.7456998825073242;
          break;
        case 9:
          r0ab = 2.6881000995635986;
          break;
        case 10:
          r0ab = 2.5727999210357666;
          break;
        case 11:
          r0ab = 3.0573999881744385;
          break;
        case 12:
          r0ab = 3.3264000415802002;
          break;
        case 13:
          r0ab = 3.3561999797821045;
          break;
        case 14:
          r0ab = 3.2528998851776123;
          break;
        case 15:
          r0ab = 3.1916000843048096;
          break;
        case 16:
          r0ab = 3.1522998809814453;
          break;
        case 17:
          r0ab = 3.1045999526977539;
          break;
      }
      break;
    case 18:
      switch (atmA) {
        case 1:
          r0ab = 2.3724999427795410;
          break;
        case 2:
          r0ab = 2.3289000988006592;
          break;
        case 3:
          r0ab = 2.8759999275207520;
          break;
        case 4:
          r0ab = 2.9804000854492188;
          break;
        case 5:
          r0ab = 2.9093000888824463;
          break;
        case 6:
          r0ab = 2.8039999008178711;
          break;
        case 7:
          r0ab = 2.7070999145507812;
          break;
        case 8:
          r0ab = 2.6386001110076904;
          break;
        case 9:
          r0ab = 2.5720000267028809;
          break;
        case 10:
          r0ab = 2.5139000415802002;
          break;
        case 11:
          r0ab = 2.9516999721527100;
          break;
        case 12:
          r0ab = 3.1605999469757080;
          break;
        case 13:
          r0ab = 3.2084999084472656;
          break;
        case 14:
          r0ab = 3.1691999435424805;
          break;
        case 15:
          r0ab = 3.0982000827789307;
          break;
        case 16:
          r0ab = 3.0352001190185547;
          break;
        case 17:
          r0ab = 2.9730000495910645;
          break;
        case 18:
          r0ab = 2.9147999286651611;
          break;
      }
      break;
    case 19:
      switch (atmA) {
        case 1:
          r0ab = 3.2146999835968018;
          break;
        case 2:
          r0ab = 2.8315000534057617;
          break;
        case 3:
          r0ab = 3.8724000453948975;
          break;
        case 4:
          r0ab = 3.4621000289916992;
          break;
        case 5:
          r0ab = 3.8822999000549316;
          break;
        case 6:
          r0ab = 3.3759999275207520;
          break;
        case 7:
          r0ab = 3.0745999813079834;
          break;
        case 8:
          r0ab = 2.8817000389099121;
          break;
        case 9:
          r0ab = 2.7551999092102051;
          break;
        case 10:
          r0ab = 2.6605000495910645;
          break;
        case 11:
          r0ab = 3.9739999771118164;
          break;
        case 12:
          r0ab = 3.6191999912261963;
          break;
        case 13:
          r0ab = 3.6568999290466309;
          break;
        case 14:
          r0ab = 3.9586000442504883;
          break;
        case 15:
          r0ab = 3.6187999248504639;
          break;
        case 16:
          r0ab = 3.3917000293731689;
          break;
        case 17:
          r0ab = 3.2479000091552734;
          break;
        case 18:
          r0ab = 3.1433999538421631;
          break;
        case 19:
          r0ab = 4.2410998344421387;
          break;
      }
      break;
    case 20:
      switch (atmA) {
        case 1:
          r0ab = 2.7597000598907471;
          break;
        case 2:
          r0ab = 3.0587999820709229;
          break;
        case 3:
          r0ab = 3.3473999500274658;
          break;
        case 4:
          r0ab = 3.6214001178741455;
          break;
        case 5:
          r0ab = 3.4353001117706299;
          break;
        case 6:
          r0ab = 3.4728999137878418;
          break;
        case 7:
          r0ab = 3.2486999034881592;
          break;
        case 8:
          r0ab = 3.3199999332427979;
          break;
        case 9:
          r0ab = 3.0913999080657959;
          break;
        case 10:
          r0ab = 2.9402999877929688;
          break;
        case 11:
          r0ab = 3.4972000122070312;
          break;
        case 12:
          r0ab = 3.7992999553680420;
          break;
        case 13:
          r0ab = 3.6772999763488770;
          break;
        case 14:
          r0ab = 3.8677999973297119;
          break;
        case 15:
          r0ab = 3.5808000564575195;
          break;
        case 16:
          r0ab = 3.8243000507354736;
          break;
        case 17:
          r0ab = 3.5826001167297363;
          break;
        case 18:
          r0ab = 3.4156000614166260;
          break;
        case 19:
          r0ab = 3.8764998912811279;
          break;
        case 20:
          r0ab = 4.1034998893737793;
          break;
      }
      break;
    case 21:
      switch (atmA) {
        case 1:
          r0ab = 2.7360999584197998;
          break;
        case 2:
          r0ab = 2.9765000343322754;
          break;
        case 3:
          r0ab = 3.2474999427795410;
          break;
        case 4:
          r0ab = 3.5004000663757324;
          break;
        case 5:
          r0ab = 3.4184999465942383;
          break;
        case 6:
          r0ab = 3.4377999305725098;
          break;
        case 7:
          r0ab = 3.2084000110626221;
          break;
        case 8:
          r0ab = 3.2787001132965088;
          break;
        case 9:
          r0ab = 3.0604000091552734;
          break;
        case 10:
          r0ab = 2.9186999797821045;
          break;
        case 11:
          r0ab = 3.4037001132965088;
          break;
        case 12:
          r0ab = 3.6758999824523926;
          break;
        case 13:
          r0ab = 3.6586000919342041;
          break;
        case 14:
          r0ab = 3.8327000141143799;
          break;
        case 15:
          r0ab = 3.5371999740600586;
          break;
        case 16:
          r0ab = 3.7664999961853027;
          break;
        case 17:
          r0ab = 3.5309998989105225;
          break;
        case 18:
          r0ab = 3.3699998855590820;
          break;
        case 19:
          r0ab = 3.7788000106811523;
          break;
        case 20:
          r0ab = 3.9804000854492188;
          break;
        case 21:
          r0ab = 3.8903000354766846;
          break;
      }
      break;
    case 22:
      switch (atmA) {
        case 1:
          r0ab = 2.6831998825073242;
          break;
        case 2:
          r0ab = 2.9059998989105225;
          break;
        case 3:
          r0ab = 3.2613000869750977;
          break;
        case 4:
          r0ab = 3.4358999729156494;
          break;
        case 5:
          r0ab = 3.3538000583648682;
          break;
        case 6:
          r0ab = 3.3859999179840088;
          break;
        case 7:
          r0ab = 3.1549999713897705;
          break;
        case 8:
          r0ab = 3.2300000190734863;
          break;
        case 9:
          r0ab = 3.0132999420166016;
          break;
        case 10:
          r0ab = 2.8736000061035156;
          break;
        case 11:
          r0ab = 3.4024000167846680;
          break;
        case 12:
          r0ab = 3.6142001152038574;
          break;
        case 13:
          r0ab = 3.5978999137878418;
          break;
        case 14:
          r0ab = 3.5295000076293945;
          break;
        case 15:
          r0ab = 3.4834001064300537;
          break;
        case 16:
          r0ab = 3.7139999866485596;
          break;
        case 17:
          r0ab = 3.4781999588012695;
          break;
        case 18:
          r0ab = 3.3169999122619629;
          break;
        case 19:
          r0ab = 3.7434000968933105;
          break;
        case 20:
          r0ab = 3.9623000621795654;
          break;
        case 21:
          r0ab = 3.8180999755859375;
          break;
        case 22:
          r0ab = 3.7641999721527100;
          break;
      }
      break;
    case 23:
      switch (atmA) {
        case 1:
          r0ab = 2.6379001140594482;
          break;
        case 2:
          r0ab = 2.8494000434875488;
          break;
        case 3:
          r0ab = 3.1840000152587891;
          break;
        case 4:
          r0ab = 3.4224998950958252;
          break;
        case 5:
          r0ab = 3.2771000862121582;
          break;
        case 6:
          r0ab = 3.3401000499725342;
          break;
        case 7:
          r0ab = 3.1071999073028564;
          break;
        case 8:
          r0ab = 3.1884999275207520;
          break;
        case 9:
          r0ab = 2.9714000225067139;
          break;
        case 10:
          r0ab = 2.8318998813629150;
          break;
        case 11:
          r0ab = 3.3315000534057617;
          break;
        case 12:
          r0ab = 3.5978999137878418;
          break;
        case 13:
          r0ab = 3.5255999565124512;
          break;
        case 14:
          r0ab = 3.4979999065399170;
          break;
        case 15:
          r0ab = 3.4375998973846436;
          break;
        case 16:
          r0ab = 3.6714000701904297;
          break;
        case 17:
          r0ab = 3.4346001148223877;
          break;
        case 18:
          r0ab = 3.2723000049591064;
          break;
        case 19:
          r0ab = 3.6858999729156494;
          break;
        case 20:
          r0ab = 3.8984999656677246;
          break;
        case 21:
          r0ab = 3.7918000221252441;
          break;
        case 22:
          r0ab = 3.7372000217437744;
          break;
        case 23:
          r0ab = 3.7211000919342041;
          break;
      }
      break;
    case 24:
      switch (atmA) {
        case 1:
          r0ab = 2.9230000972747803;
          break;
        case 2:
          r0ab = 2.6222999095916748;
          break;
        case 3:
          r0ab = 3.4161000251770020;
          break;
        case 4:
          r0ab = 2.8998999595642090;
          break;
        case 5:
          r0ab = 3.0557000637054443;
          break;
        case 6:
          r0ab = 3.3308000564575195;
          break;
        case 7:
          r0ab = 3.0555000305175781;
          break;
        case 8:
          r0ab = 2.8508000373840332;
          break;
        case 9:
          r0ab = 2.7385001182556152;
          break;
        case 10:
          r0ab = 2.6640000343322754;
          break;
        case 11:
          r0ab = 3.5262999534606934;
          break;
        case 12:
          r0ab = 3.0276999473571777;
          break;
        case 13:
          r0ab = 3.2990000247955322;
          break;
        case 14:
          r0ab = 3.7720999717712402;
          break;
        case 15:
          r0ab = 3.5016999244689941;
          break;
        case 16:
          r0ab = 3.2750999927520752;
          break;
        case 17:
          r0ab = 3.1368000507354736;
          break;
        case 18:
          r0ab = 3.0434999465942383;
          break;
        case 19:
          r0ab = 3.7873001098632812;
          break;
        case 20:
          r0ab = 3.2857999801635742;
          break;
        case 21:
          r0ab = 3.2139999866485596;
          break;
        case 22:
          r0ab = 3.1726999282836914;
          break;
        case 23:
          r0ab = 3.2177999019622803;
          break;
        case 24:
          r0ab = 3.4414000511169434;
          break;
      }
      break;
    case 25:
      switch (atmA) {
        case 1:
          r0ab = 2.5490000247955322;
          break;
        case 2:
          r0ab = 2.7623000144958496;
          break;
        case 3:
          r0ab = 3.0991001129150391;
          break;
        case 4:
          r0ab = 3.3252000808715820;
          break;
        case 5:
          r0ab = 3.1835999488830566;
          break;
        case 6:
          r0ab = 3.2427999973297119;
          break;
        case 7:
          r0ab = 3.0258998870849609;
          break;
        case 8:
          r0ab = 3.1224999427795410;
          break;
        case 9:
          r0ab = 2.9031999111175537;
          break;
        case 10:
          r0ab = 2.7620999813079834;
          break;
        case 11:
          r0ab = 3.2490000724792480;
          break;
        case 12:
          r0ab = 3.5109999179840088;
          break;
        case 13:
          r0ab = 3.4428999423980713;
          break;
        case 14:
          r0ab = 3.3845000267028809;
          break;
        case 15:
          r0ab = 3.3573999404907227;
          break;
        case 16:
          r0ab = 3.6045000553131104;
          break;
        case 17:
          r0ab = 3.3657999038696289;
          break;
        case 18:
          r0ab = 3.2012999057769775;
          break;
        case 19:
          r0ab = 3.6110000610351562;
          break;
        case 20:
          r0ab = 3.8241000175476074;
          break;
        case 21:
          r0ab = 3.7090001106262207;
          break;
        case 22:
          r0ab = 3.6496000289916992;
          break;
        case 23:
          r0ab = 3.6333000659942627;
          break;
        case 24:
          r0ab = 3.0896000862121582;
          break;
        case 25:
          r0ab = 3.5462000370025635;
          break;
      }
      break;
    case 26:
      switch (atmA) {
        case 1:
          r0ab = 2.4925999641418457;
          break;
        case 2:
          r0ab = 2.7135999202728271;
          break;
        case 3:
          r0ab = 3.0692999362945557;
          break;
        case 4:
          r0ab = 3.2699000835418701;
          break;
        case 5:
          r0ab = 3.1271998882293701;
          break;
        case 6:
          r0ab = 3.1893000602722168;
          break;
        case 7:
          r0ab = 2.9658000469207764;
          break;
        case 8:
          r0ab = 3.0971999168395996;
          break;
        case 9:
          r0ab = 2.8777999877929688;
          break;
        case 10:
          r0ab = 2.7358000278472900;
          break;
        case 11:
          r0ab = 3.2205998897552490;
          break;
        case 12:
          r0ab = 3.4565999507904053;
          break;
        case 13:
          r0ab = 3.3896000385284424;
          break;
        case 14:
          r0ab = 3.3257000446319580;
          break;
        case 15:
          r0ab = 3.2946000099182129;
          break;
        case 16:
          r0ab = 3.5692999362945557;
          break;
        case 17:
          r0ab = 3.3311998844146729;
          break;
        case 18:
          r0ab = 3.1670000553131104;
          break;
        case 19:
          r0ab = 3.5804998874664307;
          break;
        case 20:
          r0ab = 3.7711000442504883;
          break;
        case 21:
          r0ab = 3.6535999774932861;
          break;
        case 22:
          r0ab = 3.5927000045776367;
          break;
        case 23:
          r0ab = 3.5775001049041748;
          break;
        case 24:
          r0ab = 3.0411000251770020;
          break;
        case 25:
          r0ab = 3.4885001182556152;
          break;
        case 26:
          r0ab = 3.4421000480651855;
          break;
      }
      break;
    case 27:
      switch (atmA) {
        case 1:
          r0ab = 2.4667000770568848;
          break;
        case 2:
          r0ab = 2.6709001064300537;
          break;
        case 3:
          r0ab = 3.0574998855590820;
          break;
        case 4:
          r0ab = 3.2356998920440674;
          break;
        case 5:
          r0ab = 3.0908000469207764;
          break;
        case 6:
          r0ab = 3.1537001132965088;
          break;
        case 7:
          r0ab = 2.9235000610351562;
          break;
        case 8:
          r0ab = 3.0669000148773193;
          break;
        case 9:
          r0ab = 2.8475999832153320;
          break;
        case 10:
          r0ab = 2.7053999900817871;
          break;
        case 11:
          r0ab = 3.2063999176025391;
          break;
        case 12:
          r0ab = 3.4519000053405762;
          break;
        case 13:
          r0ab = 3.3592998981475830;
          break;
        case 14:
          r0ab = 3.2920999526977539;
          break;
        case 15:
          r0ab = 3.2576999664306641;
          break;
        case 16:
          r0ab = 3.2160999774932861;
          break;
        case 17:
          r0ab = 3.2981998920440674;
          break;
        case 18:
          r0ab = 3.1338999271392822;
          break;
        case 19:
          r0ab = 3.5606000423431396;
          break;
        case 20:
          r0ab = 3.7581999301910400;
          break;
        case 21:
          r0ab = 3.6431999206542969;
          break;
        case 22:
          r0ab = 3.5833001136779785;
          break;
        case 23:
          r0ab = 3.5690999031066895;
          break;
        case 24:
          r0ab = 3.0160999298095703;
          break;
        case 25:
          r0ab = 3.4811999797821045;
          break;
        case 26:
          r0ab = 3.4339001178741455;
          break;
        case 27:
          r0ab = 3.4326999187469482;
          break;
      }
      break;
    case 28:
      switch (atmA) {
        case 1:
          r0ab = 2.4514999389648438;
          break;
        case 2:
          r0ab = 2.6338000297546387;
          break;
        case 3:
          r0ab = 3.0511000156402588;
          break;
        case 4:
          r0ab = 3.2228999137878418;
          break;
        case 5:
          r0ab = 3.0629999637603760;
          break;
        case 6:
          r0ab = 3.1264998912811279;
          break;
        case 7:
          r0ab = 2.8908998966217041;
          break;
        case 8:
          r0ab = 3.0253000259399414;
          break;
        case 9:
          r0ab = 2.8183999061584473;
          break;
        case 10:
          r0ab = 2.6763999462127686;
          break;
        case 11:
          r0ab = 3.1967999935150146;
          break;
        case 12:
          r0ab = 3.4114000797271729;
          break;
        case 13:
          r0ab = 3.3492000102996826;
          break;
        case 14:
          r0ab = 3.2690999507904053;
          break;
        case 15:
          r0ab = 3.2320001125335693;
          break;
        case 16:
          r0ab = 3.1786000728607178;
          break;
        case 17:
          r0ab = 3.2679998874664307;
          break;
        case 18:
          r0ab = 3.1036000251770020;
          break;
        case 19:
          r0ab = 3.5453000068664551;
          break;
        case 20:
          r0ab = 3.7258999347686768;
          break;
        case 21:
          r0ab = 3.6089999675750732;
          break;
        case 22:
          r0ab = 3.5473001003265381;
          break;
        case 23:
          r0ab = 3.5327000617980957;
          break;
        case 24:
          r0ab = 3.0018000602722168;
          break;
        case 25:
          r0ab = 3.4412999153137207;
          break;
        case 26:
          r0ab = 3.3907001018524170;
          break;
        case 27:
          r0ab = 3.3592998981475830;
          break;
        case 28:
          r0ab = 3.3461999893188477;
          break;
      }
      break;
    case 29:
      switch (atmA) {
        case 1:
          r0ab = 2.4412999153137207;
          break;
        case 2:
          r0ab = 2.6006000041961670;
          break;
        case 3:
          r0ab = 3.0539999008178711;
          break;
        case 4:
          r0ab = 3.1986999511718750;
          break;
        case 5:
          r0ab = 3.0490000247955322;
          break;
        case 6:
          r0ab = 3.1057999134063721;
          break;
        case 7:
          r0ab = 2.8643000125885010;
          break;
        case 8:
          r0ab = 2.9948000907897949;
          break;
        case 9:
          r0ab = 2.7908000946044922;
          break;
        case 10:
          r0ab = 2.6491000652313232;
          break;
        case 11:
          r0ab = 3.1949999332427979;
          break;
        case 12:
          r0ab = 3.3921999931335449;
          break;
        case 13:
          r0ab = 3.3315999507904053;
          break;
        case 14:
          r0ab = 3.2585000991821289;
          break;
        case 15:
          r0ab = 3.2135999202728271;
          break;
        case 16:
          r0ab = 3.1515998840332031;
          break;
        case 17:
          r0ab = 3.2363998889923096;
          break;
        case 18:
          r0ab = 3.0752000808715820;
          break;
        case 19:
          r0ab = 3.5367999076843262;
          break;
        case 20:
          r0ab = 3.7116999626159668;
          break;
        case 21:
          r0ab = 3.5940999984741211;
          break;
        case 22:
          r0ab = 3.5313000679016113;
          break;
        case 23:
          r0ab = 3.5164000988006592;
          break;
        case 24:
          r0ab = 2.9962000846862793;
          break;
        case 25:
          r0ab = 3.4224998950958252;
          break;
        case 26:
          r0ab = 3.3698999881744385;
          break;
        case 27:
          r0ab = 3.3369998931884766;
          break;
        case 28:
          r0ab = 3.3234000205993652;
          break;
        case 29:
          r0ab = 3.3008000850677490;
          break;
      }
      break;
    case 30:
      switch (atmA) {
        case 1:
          r0ab = 2.4317998886108398;
          break;
        case 2:
          r0ab = 2.5729000568389893;
          break;
        case 3:
          r0ab = 3.0415999889373779;
          break;
        case 4:
          r0ab = 3.1638998985290527;
          break;
        case 5:
          r0ab = 3.0195999145507812;
          break;
        case 6:
          r0ab = 3.0843000411987305;
          break;
        case 7:
          r0ab = 2.8413000106811523;
          break;
        case 8:
          r0ab = 2.7435998916625977;
          break;
        case 9:
          r0ab = 2.7607998847961426;
          break;
        case 10:
          r0ab = 2.6270999908447266;
          break;
        case 11:
          r0ab = 3.1810998916625977;
          break;
        case 12:
          r0ab = 3.3591001033782959;
          break;
        case 13:
          r0ab = 3.3045001029968262;
          break;
        case 14:
          r0ab = 3.2348999977111816;
          break;
        case 15:
          r0ab = 3.1942000389099121;
          break;
        case 16:
          r0ab = 3.1291000843048096;
          break;
        case 17:
          r0ab = 3.2111001014709473;
          break;
        case 18:
          r0ab = 3.0534000396728516;
          break;
        case 19:
          r0ab = 3.5188999176025391;
          break;
        case 20:
          r0ab = 3.6809000968933105;
          break;
        case 21:
          r0ab = 3.5634999275207520;
          break;
        case 22:
          r0ab = 3.5000998973846436;
          break;
        case 23:
          r0ab = 3.4853999614715576;
          break;
        case 24:
          r0ab = 2.9856998920440674;
          break;
        case 25:
          r0ab = 3.3896999359130859;
          break;
        case 26:
          r0ab = 3.3362998962402344;
          break;
        case 27:
          r0ab = 3.3027000427246094;
          break;
        case 28:
          r0ab = 3.2890000343322754;
          break;
        case 29:
          r0ab = 3.2655000686645508;
          break;
        case 30:
          r0ab = 3.2309000492095947;
          break;
      }
      break;
    case 31:
      switch (atmA) {
        case 1:
          r0ab = 2.8501999378204346;
          break;
        case 2:
          r0ab = 2.6933999061584473;
          break;
        case 3:
          r0ab = 3.2467000484466553;
          break;
        case 4:
          r0ab = 3.1921000480651855;
          break;
        case 5:
          r0ab = 3.5662999153137207;
          break;
        case 6:
          r0ab = 3.2541000843048096;
          break;
        case 7:
          r0ab = 3.0571000576019287;
          break;
        case 8:
          r0ab = 2.9047999382019043;
          break;
        case 9:
          r0ab = 2.8657000064849854;
          break;
        case 10:
          r0ab = 2.7437999248504639;
          break;
        case 11:
          r0ab = 3.3547000885009766;
          break;
        case 12:
          r0ab = 3.3510000705718994;
          break;
        case 13:
          r0ab = 3.9837000370025635;
          break;
        case 14:
          r0ab = 3.6870999336242676;
          break;
        case 15:
          r0ab = 3.4862000942230225;
          break;
        case 16:
          r0ab = 3.3389000892639160;
          break;
        case 17:
          r0ab = 3.2413001060485840;
          break;
        case 18:
          r0ab = 3.1707999706268311;
          break;
        case 19:
          r0ab = 3.6096000671386719;
          break;
        case 20:
          r0ab = 3.6280000209808350;
          break;
        case 21:
          r0ab = 3.6860001087188721;
          break;
        case 22:
          r0ab = 3.5567998886108398;
          break;
        case 23:
          r0ab = 3.4835999011993408;
          break;
        case 24:
          r0ab = 3.2867999076843262;
          break;
        case 25:
          r0ab = 3.3993999958038330;
          break;
        case 26:
          r0ab = 3.3475999832153320;
          break;
        case 27:
          r0ab = 3.3169999122619629;
          break;
        case 28:
          r0ab = 3.2950000762939453;
          break;
        case 29:
          r0ab = 3.2874000072479248;
          break;
        case 30:
          r0ab = 3.2606000900268555;
          break;
        case 31:
          r0ab = 3.9579000473022461;
          break;
      }
      break;
    case 32:
      switch (atmA) {
        case 1:
          r0ab = 2.9226000308990479;
          break;
        case 2:
          r0ab = 2.6837999820709229;
          break;
        case 3:
          r0ab = 3.7867000102996826;
          break;
        case 4:
          r0ab = 3.1731998920440674;
          break;
        case 5:
          r0ab = 3.3872001171112061;
          break;
        case 6:
          r0ab = 3.3643000125885010;
          break;
        case 7:
          r0ab = 3.1266999244689941;
          break;
        case 8:
          r0ab = 2.9540998935699463;
          break;
        case 9:
          r0ab = 2.8505001068115234;
          break;
        case 10:
          r0ab = 2.7781000137329102;
          break;
        case 11:
          r0ab = 3.8475000858306885;
          break;
        case 12:
          r0ab = 3.3336000442504883;
          break;
        case 13:
          r0ab = 3.7358999252319336;
          break;
        case 14:
          r0ab = 3.8266000747680664;
          break;
        case 15:
          r0ab = 3.5732998847961426;
          break;
        case 16:
          r0ab = 3.3959000110626221;
          break;
        case 17:
          r0ab = 3.2774999141693115;
          break;
        case 18:
          r0ab = 3.1914999485015869;
          break;
        case 19:
          r0ab = 3.9877998828887939;
          break;
        case 20:
          r0ab = 3.8815999031066895;
          break;
        case 21:
          r0ab = 3.5810000896453857;
          break;
        case 22:
          r0ab = 3.5364000797271729;
          break;
        case 23:
          r0ab = 3.5060000419616699;
          break;
        case 24:
          r0ab = 3.8097000122070312;
          break;
        case 25:
          r0ab = 3.3924999237060547;
          break;
        case 26:
          r0ab = 3.3348000049591064;
          break;
        case 27:
          r0ab = 3.3018999099731445;
          break;
        case 28:
          r0ab = 3.2795999050140381;
          break;
        case 29:
          r0ab = 3.2662000656127930;
          break;
        case 30:
          r0ab = 3.2464001178741455;
          break;
        case 31:
          r0ab = 3.7135999202728271;
          break;
        case 32:
          r0ab = 3.8619000911712646;
          break;
      }
      break;
    case 33:
      switch (atmA) {
        case 1:
          r0ab = 2.9140000343322754;
          break;
        case 2:
          r0ab = 2.6270999908447266;
          break;
        case 3:
          r0ab = 3.4770998954772949;
          break;
        case 4:
          r0ab = 3.1774001121520996;
          break;
        case 5:
          r0ab = 3.2560000419616699;
          break;
        case 6:
          r0ab = 3.1970000267028809;
          break;
        case 7:
          r0ab = 3.1206998825073242;
          break;
        case 8:
          r0ab = 2.9405999183654785;
          break;
        case 9:
          r0ab = 2.8322000503540039;
          break;
        case 10:
          r0ab = 2.7571001052856445;
          break;
        case 11:
          r0ab = 3.5455000400543213;
          break;
        case 12:
          r0ab = 3.3513998985290527;
          break;
        case 13:
          r0ab = 3.5836999416351318;
          break;
        case 14:
          r0ab = 3.6177000999450684;
          break;
        case 15:
          r0ab = 3.5815999507904053;
          break;
        case 16:
          r0ab = 3.3901998996734619;
          break;
        case 17:
          r0ab = 3.2604000568389893;
          break;
        case 18:
          r0ab = 3.1651999950408936;
          break;
        case 19:
          r0ab = 3.7037000656127930;
          break;
        case 20:
          r0ab = 3.6282999515533447;
          break;
        case 21:
          r0ab = 3.5857999324798584;
          break;
        case 22:
          r0ab = 3.5329999923706055;
          break;
        case 23:
          r0ab = 3.4883999824523926;
          break;
        case 24:
          r0ab = 3.5789000988006592;
          break;
        case 25:
          r0ab = 3.4093999862670898;
          break;
        case 26:
          r0ab = 3.3473000526428223;
          break;
        case 27:
          r0ab = 3.3118000030517578;
          break;
        case 28:
          r0ab = 3.2876000404357910;
          break;
        case 29:
          r0ab = 3.2706999778747559;
          break;
        case 30:
          r0ab = 3.2520999908447266;
          break;
        case 31:
          r0ab = 3.5569999217987061;
          break;
        case 32:
          r0ab = 3.6496000289916992;
          break;
        case 33:
          r0ab = 3.6624999046325684;
          break;
      }
      break;
    case 34:
      switch (atmA) {
        case 1:
          r0ab = 2.7300000190734863;
          break;
        case 2:
          r0ab = 2.5869998931884766;
          break;
        case 3:
          r0ab = 3.2471001148223877;
          break;
        case 4:
          r0ab = 3.1486999988555908;
          break;
        case 5:
          r0ab = 3.1666998863220215;
          break;
        case 6:
          r0ab = 3.0913999080657959;
          break;
        case 7:
          r0ab = 3.0106999874114990;
          break;
        case 8:
          r0ab = 2.9811999797821045;
          break;
        case 9:
          r0ab = 2.8299999237060547;
          break;
        case 10:
          r0ab = 2.7283999919891357;
          break;
        case 11:
          r0ab = 3.3259000778198242;
          break;
        case 12:
          r0ab = 3.3182001113891602;
          break;
        case 13:
          r0ab = 3.4707000255584717;
          break;
        case 14:
          r0ab = 3.4748001098632812;
          break;
        case 15:
          r0ab = 3.4279000759124756;
          break;
        case 16:
          r0ab = 3.4182000160217285;
          break;
        case 17:
          r0ab = 3.2546999454498291;
          break;
        case 18:
          r0ab = 3.1352999210357666;
          break;
        case 19:
          r0ab = 3.5116000175476074;
          break;
        case 20:
          r0ab = 3.9432001113891602;
          break;
        case 21:
          r0ab = 3.8828001022338867;
          break;
        case 22:
          r0ab = 3.8303000926971436;
          break;
        case 23:
          r0ab = 3.7880001068115234;
          break;
        case 24:
          r0ab = 3.3759999275207520;
          break;
        case 25:
          r0ab = 3.7218000888824463;
          break;
        case 26:
          r0ab = 3.3408000469207764;
          break;
        case 27:
          r0ab = 3.3059000968933105;
          break;
        case 28:
          r0ab = 3.2697999477386475;
          break;
        case 29:
          r0ab = 3.2446000576019287;
          break;
        case 30:
          r0ab = 3.2228999137878418;
          break;
        case 31:
          r0ab = 3.4421999454498291;
          break;
        case 32:
          r0ab = 3.5023000240325928;
          break;
        case 33:
          r0ab = 3.5009000301361084;
          break;
        case 34:
          r0ab = 3.5267999172210693;
          break;
      }
      break;
    case 35:
      switch (atmA) {
        case 1:
          r0ab = 2.6026000976562500;
          break;
        case 2:
          r0ab = 2.5355000495910645;
          break;
        case 3:
          r0ab = 3.1129000186920166;
          break;
        case 4:
          r0ab = 3.2862999439239502;
          break;
        case 5:
          r0ab = 3.1029000282287598;
          break;
        case 6:
          r0ab = 3.0107998847961426;
          break;
        case 7:
          r0ab = 2.9226999282836914;
          break;
        case 8:
          r0ab = 2.8694000244140625;
          break;
        case 9:
          r0ab = 2.8108999729156494;
          break;
        case 10:
          r0ab = 2.6928999423980713;
          break;
        case 11:
          r0ab = 3.1958000659942627;
          break;
        case 12:
          r0ab = 3.4670000076293945;
          break;
        case 13:
          r0ab = 3.4017999172210693;
          break;
        case 14:
          r0ab = 3.3805000782012939;
          break;
        case 15:
          r0ab = 3.3217999935150146;
          break;
        case 16:
          r0ab = 3.2815001010894775;
          break;
        case 17:
          r0ab = 3.2346000671386719;
          break;
        case 18:
          r0ab = 3.0994000434875488;
          break;
        case 19:
          r0ab = 3.3936998844146729;
          break;
        case 20:
          r0ab = 3.7265999317169189;
          break;
        case 21:
          r0ab = 3.6696999073028564;
          break;
        case 22:
          r0ab = 3.6164000034332275;
          break;
        case 23:
          r0ab = 3.5729999542236328;
          break;
        case 24:
          r0ab = 3.2521998882293701;
          break;
        case 25:
          r0ab = 3.5051000118255615;
          break;
        case 26:
          r0ab = 3.4686000347137451;
          break;
        case 27:
          r0ab = 3.4354999065399170;
          break;
        case 28:
          r0ab = 3.4084000587463379;
          break;
        case 29:
          r0ab = 3.3747999668121338;
          break;
        case 30:
          r0ab = 3.3496000766754150;
          break;
        case 31:
          r0ab = 3.3691999912261963;
          break;
        case 32:
          r0ab = 3.4052000045776367;
          break;
        case 33:
          r0ab = 3.3910000324249268;
          break;
        case 34:
          r0ab = 3.3849000930786133;
          break;
        case 35:
          r0ab = 3.3661999702453613;
          break;
      }
      break;
    case 36:
      switch (atmA) {
        case 1:
          r0ab = 2.5086998939514160;
          break;
        case 2:
          r0ab = 2.4814000129699707;
          break;
        case 3:
          r0ab = 3.0239000320434570;
          break;
        case 4:
          r0ab = 3.1312000751495361;
          break;
        case 5:
          r0ab = 3.0534999370574951;
          break;
        case 6:
          r0ab = 2.9456999301910400;
          break;
        case 7:
          r0ab = 2.8496000766754150;
          break;
        case 8:
          r0ab = 2.7780001163482666;
          break;
        case 9:
          r0ab = 2.7827999591827393;
          break;
        case 10:
          r0ab = 2.6531999111175537;
          break;
        case 11:
          r0ab = 3.1063001155853271;
          break;
        case 12:
          r0ab = 3.3143000602722168;
          break;
        case 13:
          r0ab = 3.3548998832702637;
          break;
        case 14:
          r0ab = 3.3120000362396240;
          break;
        case 15:
          r0ab = 3.2421000003814697;
          break;
        case 16:
          r0ab = 3.1786999702453613;
          break;
        case 17:
          r0ab = 3.1175999641418457;
          break;
        case 18:
          r0ab = 3.0613000392913818;
          break;
        case 19:
          r0ab = 3.3081998825073242;
          break;
        case 20:
          r0ab = 3.5755000114440918;
          break;
        case 21:
          r0ab = 3.5222001075744629;
          break;
        case 22:
          r0ab = 3.4677999019622803;
          break;
        case 23:
          r0ab = 3.4230999946594238;
          break;
        case 24:
          r0ab = 3.1684000492095947;
          break;
        case 25:
          r0ab = 3.3527998924255371;
          break;
        case 26:
          r0ab = 3.3162000179290771;
          break;
        case 27:
          r0ab = 3.2827000617980957;
          break;
        case 28:
          r0ab = 3.2527000904083252;
          break;
        case 29:
          r0ab = 3.2307999134063721;
          break;
        case 30:
          r0ab = 3.2028999328613281;
          break;
        case 31:
          r0ab = 3.3173000812530518;
          break;
        case 32:
          r0ab = 3.3343000411987305;
          break;
        case 33:
          r0ab = 3.3092000484466553;
          break;
        case 34:
          r0ab = 3.2795000076293945;
          break;
        case 35:
          r0ab = 3.2451999187469482;
          break;
        case 36:
          r0ab = 3.2095999717712402;
          break;
      }
      break;
    case 37:
      switch (atmA) {
        case 1:
          r0ab = 3.2892999649047852;
          break;
        case 2:
          r0ab = 2.8991000652313232;
          break;
        case 3:
          r0ab = 4.0387997627258301;
          break;
        case 4:
          r0ab = 3.6099998950958252;
          break;
        case 5:
          r0ab = 3.9388000965118408;
          break;
        case 6:
          r0ab = 3.4474999904632568;
          break;
        case 7:
          r0ab = 3.1589999198913574;
          break;
        case 8:
          r0ab = 2.9811999797821045;
          break;
        case 9:
          r0ab = 2.8585999011993408;
          break;
        case 10:
          r0ab = 2.7683000564575195;
          break;
        case 11:
          r0ab = 4.1427998542785645;
          break;
        case 12:
          r0ab = 3.7911000251770020;
          break;
        case 13:
          r0ab = 3.8224999904632568;
          break;
        case 14:
          r0ab = 4.0371999740600586;
          break;
        case 15:
          r0ab = 3.7058999538421631;
          break;
        case 16:
          r0ab = 3.4934999942779541;
          break;
        case 17:
          r0ab = 3.3529000282287598;
          break;
        case 18:
          r0ab = 3.2492001056671143;
          break;
        case 19:
          r0ab = 4.4352002143859863;
          break;
        case 20:
          r0ab = 4.0826001167297363;
          break;
        case 21:
          r0ab = 3.9732999801635742;
          break;
        case 22:
          r0ab = 3.9254000186920166;
          break;
        case 23:
          r0ab = 3.8645999431610107;
          break;
        case 24:
          r0ab = 3.9314999580383301;
          break;
        case 25:
          r0ab = 3.7836999893188477;
          break;
        case 26:
          r0ab = 3.7465000152587891;
          break;
        case 27:
          r0ab = 3.7211000919342041;
          break;
        case 28:
          r0ab = 3.7012000083923340;
          break;
        case 29:
          r0ab = 3.6893000602722168;
          break;
        case 30:
          r0ab = 3.6675999164581299;
          break;
        case 31:
          r0ab = 3.7736001014709473;
          break;
        case 32:
          r0ab = 4.0659999847412109;
          break;
        case 33:
          r0ab = 3.7925999164581299;
          break;
        case 34:
          r0ab = 3.6157999038696289;
          break;
        case 35:
          r0ab = 3.5016999244689941;
          break;
        case 36:
          r0ab = 3.4165999889373779;
          break;
        case 37:
          r0ab = 4.6175999641418457;
          break;
      }
      break;
    case 38:
      switch (atmA) {
        case 1:
          r0ab = 2.8785998821258545;
          break;
        case 2:
          r0ab = 3.1658000946044922;
          break;
        case 3:
          r0ab = 3.5822999477386475;
          break;
        case 4:
          r0ab = 3.7688999176025391;
          break;
        case 5:
          r0ab = 3.5762000083923340;
          break;
        case 6:
          r0ab = 3.5789000988006592;
          break;
        case 7:
          r0ab = 3.3552000522613525;
          break;
        case 8:
          r0ab = 3.4003999233245850;
          break;
        case 9:
          r0ab = 3.1721999645233154;
          break;
        case 10:
          r0ab = 3.0211999416351318;
          break;
        case 11:
          r0ab = 3.7241001129150391;
          break;
        case 12:
          r0ab = 3.9604001045227051;
          break;
        case 13:
          r0ab = 3.8499999046325684;
          break;
        case 14:
          r0ab = 3.9844000339508057;
          break;
        case 15:
          r0ab = 3.7035000324249268;
          break;
        case 16:
          r0ab = 3.9161000251770020;
          break;
        case 17:
          r0ab = 3.6751000881195068;
          break;
        case 18:
          r0ab = 3.5074999332427979;
          break;
        case 19:
          r0ab = 4.1150999069213867;
          break;
        case 20:
          r0ab = 4.2877001762390137;
          break;
        case 21:
          r0ab = 4.1578998565673828;
          break;
        case 22:
          r0ab = 4.1247000694274902;
          break;
        case 23:
          r0ab = 4.0616998672485352;
          break;
        case 24:
          r0ab = 3.4874000549316406;
          break;
        case 25:
          r0ab = 3.9848001003265381;
          break;
        case 26:
          r0ab = 3.9279999732971191;
          break;
        case 27:
          r0ab = 3.9079000949859619;
          break;
        case 28:
          r0ab = 3.8750998973846436;
          break;
        case 29:
          r0ab = 3.8603999614715576;
          break;
        case 30:
          r0ab = 3.8276998996734619;
          break;
        case 31:
          r0ab = 3.8001999855041504;
          break;
        case 32:
          r0ab = 3.9981000423431396;
          break;
        case 33:
          r0ab = 3.7544000148773193;
          break;
        case 34:
          r0ab = 4.0370998382568359;
          break;
        case 35:
          r0ab = 3.8224999904632568;
          break;
        case 36:
          r0ab = 3.6717998981475830;
          break;
        case 37:
          r0ab = 4.3091998100280762;
          break;
        case 38:
          r0ab = 4.4763998985290527;
          break;
      }
      break;
    case 39:
      switch (atmA) {
        case 1:
          r0ab = 2.8996999263763428;
          break;
        case 2:
          r0ab = 3.0952999591827393;
          break;
        case 3:
          r0ab = 3.4523999691009521;
          break;
        case 4:
          r0ab = 3.6106998920440674;
          break;
        case 5:
          r0ab = 3.6061999797821045;
          break;
        case 6:
          r0ab = 3.5782999992370605;
          break;
        case 7:
          r0ab = 3.3462998867034912;
          break;
        case 8:
          r0ab = 3.3854999542236328;
          break;
        case 9:
          r0ab = 3.1745998859405518;
          break;
        case 10:
          r0ab = 3.0381000041961670;
          break;
        case 11:
          r0ab = 3.6019001007080078;
          break;
        case 12:
          r0ab = 3.7938001155853271;
          break;
        case 13:
          r0ab = 3.8696999549865723;
          break;
        case 14:
          r0ab = 3.9781000614166260;
          break;
        case 15:
          r0ab = 3.6877000331878662;
          break;
        case 16:
          r0ab = 3.8736000061035156;
          break;
        case 17:
          r0ab = 3.6451001167297363;
          break;
        case 18:
          r0ab = 3.4890000820159912;
          break;
        case 19:
          r0ab = 3.9858000278472900;
          break;
        case 20:
          r0ab = 4.1178998947143555;
          break;
        case 21:
          r0ab = 4.0430002212524414;
          break;
        case 22:
          r0ab = 3.9563000202178955;
          break;
        case 23:
          r0ab = 3.9182000160217285;
          break;
        case 24:
          r0ab = 3.4001998901367188;
          break;
        case 25:
          r0ab = 3.8310000896453857;
          break;
        case 26:
          r0ab = 3.7716000080108643;
          break;
        case 27:
          r0ab = 3.7543001174926758;
          break;
        case 28:
          r0ab = 3.7202999591827393;
          break;
        case 29:
          r0ab = 3.7053000926971436;
          break;
        case 30:
          r0ab = 3.6742000579833984;
          break;
        case 31:
          r0ab = 3.8317999839782715;
          break;
        case 32:
          r0ab = 3.7630999088287354;
          break;
        case 33:
          r0ab = 3.7392001152038574;
          break;
        case 34:
          r0ab = 3.9892001152038574;
          break;
        case 35:
          r0ab = 3.7832000255584717;
          break;
        case 36:
          r0ab = 3.6405999660491943;
          break;
        case 37:
          r0ab = 4.1701002120971680;
          break;
        case 38:
          r0ab = 4.3015999794006348;
          break;
        case 39:
          r0ab = 4.2196002006530762;
          break;
      }
      break;
    case 40:
      switch (atmA) {
        case 1:
          r0ab = 2.8534998893737793;
          break;
        case 2:
          r0ab = 3.0167000293731689;
          break;
        case 3:
          r0ab = 3.3977999687194824;
          break;
        case 4:
          r0ab = 3.5362999439239502;
          break;
        case 5:
          r0ab = 3.5392999649047852;
          break;
        case 6:
          r0ab = 3.5301001071929932;
          break;
        case 7:
          r0ab = 3.2960000038146973;
          break;
        case 8:
          r0ab = 3.3352000713348389;
          break;
        case 9:
          r0ab = 3.1287000179290771;
          break;
        case 10:
          r0ab = 2.9967000484466553;
          break;
        case 11:
          r0ab = 3.6658999919891357;
          break;
        case 12:
          r0ab = 3.7239000797271729;
          break;
        case 13:
          r0ab = 3.8069999217987061;
          break;
        case 14:
          r0ab = 3.7165000438690186;
          break;
        case 15:
          r0ab = 3.6368000507354736;
          break;
        case 16:
          r0ab = 3.8162000179290771;
          break;
        case 17:
          r0ab = 3.5885000228881836;
          break;
        case 18:
          r0ab = 3.4335999488830566;
          break;
        case 19:
          r0ab = 3.9828999042510986;
          break;
        case 20:
          r0ab = 4.0528998374938965;
          break;
        case 21:
          r0ab = 3.9584000110626221;
          break;
        case 22:
          r0ab = 3.9024999141693115;
          break;
        case 23:
          r0ab = 3.8606998920440674;
          break;
        case 24:
          r0ab = 3.3673000335693359;
          break;
        case 25:
          r0ab = 3.7657999992370605;
          break;
        case 26:
          r0ab = 3.7035000324249268;
          break;
        case 27:
          r0ab = 3.6865999698638916;
          break;
        case 28:
          r0ab = 3.6503999233245850;
          break;
        case 29:
          r0ab = 3.6338999271392822;
          break;
        case 30:
          r0ab = 3.6024000644683838;
          break;
        case 31:
          r0ab = 3.7708001136779785;
          break;
        case 32:
          r0ab = 3.7283000946044922;
          break;
        case 33:
          r0ab = 3.6895999908447266;
          break;
        case 34:
          r0ab = 3.9314999580383301;
          break;
        case 35:
          r0ab = 3.7249999046325684;
          break;
        case 36:
          r0ab = 3.5818998813629150;
          break;
        case 37:
          r0ab = 4.1456999778747559;
          break;
        case 38:
          r0ab = 4.2280001640319824;
          break;
        case 39:
          r0ab = 4.1129999160766602;
          break;
        case 40:
          r0ab = 4.0597000122070312;
          break;
      }
      break;
    case 41:
      switch (atmA) {
        case 1:
          r0ab = 3.0905001163482666;
          break;
        case 2:
          r0ab = 2.7997999191284180;
          break;
        case 3:
          r0ab = 3.6447999477386475;
          break;
        case 4:
          r0ab = 3.0738999843597412;
          break;
        case 5:
          r0ab = 3.2995998859405518;
          break;
        case 6:
          r0ab = 3.5262000560760498;
          break;
        case 7:
          r0ab = 3.2558999061584473;
          break;
        case 8:
          r0ab = 3.0518000125885010;
          break;
        case 9:
          r0ab = 2.9393999576568604;
          break;
        case 10:
          r0ab = 2.8657999038696289;
          break;
        case 11:
          r0ab = 3.7513999938964844;
          break;
        case 12:
          r0ab = 3.2295000553131104;
          break;
        case 13:
          r0ab = 3.5643000602722168;
          break;
        case 14:
          r0ab = 3.7808001041412354;
          break;
        case 15:
          r0ab = 3.6930999755859375;
          break;
        case 16:
          r0ab = 3.4723000526428223;
          break;
        case 17:
          r0ab = 3.3357000350952148;
          break;
        case 18:
          r0ab = 3.2428998947143555;
          break;
        case 19:
          r0ab = 4.0279998779296875;
          break;
        case 20:
          r0ab = 3.5589001178741455;
          break;
        case 21:
          r0ab = 3.4635999202728271;
          break;
        case 22:
          r0ab = 3.4993999004364014;
          break;
        case 23:
          r0ab = 3.4309000968933105;
          break;
        case 24:
          r0ab = 3.6177000999450684;
          break;
        case 25:
          r0ab = 3.2946000099182129;
          break;
        case 26:
          r0ab = 3.2376000881195068;
          break;
        case 27:
          r0ab = 3.2049999237060547;
          break;
        case 28:
          r0ab = 3.1847000122070312;
          break;
        case 29:
          r0ab = 3.1714999675750732;
          break;
        case 30:
          r0ab = 3.1598999500274658;
          break;
        case 31:
          r0ab = 3.5555000305175781;
          break;
        case 32:
          r0ab = 3.8111000061035156;
          break;
        case 33:
          r0ab = 3.7692999839782715;
          break;
        case 34:
          r0ab = 3.5717999935150146;
          break;
        case 35:
          r0ab = 3.4498000144958496;
          break;
        case 36:
          r0ab = 3.3661999702453613;
          break;
        case 37:
          r0ab = 4.1607999801635742;
          break;
        case 38:
          r0ab = 3.7416999340057373;
          break;
        case 39:
          r0ab = 3.6535999774932861;
          break;
        case 40:
          r0ab = 3.6154000759124756;
          break;
        case 41:
          r0ab = 3.8596000671386719;
          break;
      }
      break;
    case 42:
      switch (atmA) {
        case 1:
          r0ab = 3.0301001071929932;
          break;
        case 2:
          r0ab = 2.7311999797821045;
          break;
        case 3:
          r0ab = 3.5820999145507812;
          break;
        case 4:
          r0ab = 3.0473001003265381;
          break;
        case 5:
          r0ab = 3.2137000560760498;
          break;
        case 6:
          r0ab = 3.4679000377655029;
          break;
        case 7:
          r0ab = 3.1974999904632568;
          break;
        case 8:
          r0ab = 2.9969000816345215;
          break;
        case 9:
          r0ab = 2.8847000598907471;
          break;
        case 10:
          r0ab = 2.8110001087188721;
          break;
        case 11:
          r0ab = 3.6930999755859375;
          break;
        case 12:
          r0ab = 3.2076001167297363;
          break;
        case 13:
          r0ab = 3.4942998886108398;
          break;
        case 14:
          r0ab = 3.5955998897552490;
          break;
        case 15:
          r0ab = 3.6379001140594482;
          break;
        case 16:
          r0ab = 3.4189999103546143;
          break;
        case 17:
          r0ab = 3.2808001041412354;
          break;
        case 18:
          r0ab = 3.1860001087188721;
          break;
        case 19:
          r0ab = 3.9849998950958252;
          break;
        case 20:
          r0ab = 3.5104999542236328;
          break;
        case 21:
          r0ab = 3.4330000877380371;
          break;
        case 22:
          r0ab = 3.3796999454498291;
          break;
        case 23:
          r0ab = 3.4154999256134033;
          break;
        case 24:
          r0ab = 3.6033000946044922;
          break;
        case 25:
          r0ab = 3.2736999988555908;
          break;
        case 26:
          r0ab = 3.2144999504089355;
          break;
        case 27:
          r0ab = 3.1807000637054443;
          break;
        case 28:
          r0ab = 3.1596000194549561;
          break;
        case 29:
          r0ab = 3.1461000442504883;
          break;
        case 30:
          r0ab = 3.1336998939514160;
          break;
        case 31:
          r0ab = 3.4811999797821045;
          break;
        case 32:
          r0ab = 3.6250998973846436;
          break;
        case 33:
          r0ab = 3.7151999473571777;
          break;
        case 34:
          r0ab = 3.5201001167297363;
          break;
        case 35:
          r0ab = 3.3966000080108643;
          break;
        case 36:
          r0ab = 3.3106999397277832;
          break;
        case 37:
          r0ab = 4.1128001213073730;
          break;
        case 38:
          r0ab = 3.6898999214172363;
          break;
        case 39:
          r0ab = 3.6082000732421875;
          break;
        case 40:
          r0ab = 3.5604000091552734;
          break;
        case 41:
          r0ab = 3.7834000587463379;
          break;
        case 42:
          r0ab = 3.7543001174926758;
          break;
      }
      break;
    case 43:
      switch (atmA) {
        case 1:
          r0ab = 2.9189000129699707;
          break;
        case 2:
          r0ab = 2.6777000427246094;
          break;
        case 3:
          r0ab = 3.4925000667572021;
          break;
        case 4:
          r0ab = 2.9647998809814453;
          break;
        case 5:
          r0ab = 3.1215999126434326;
          break;
        case 6:
          r0ab = 3.2939999103546143;
          break;
        case 7:
          r0ab = 3.0975000858306885;
          break;
        case 8:
          r0ab = 2.9756999015808105;
          break;
        case 9:
          r0ab = 2.8492999076843262;
          break;
        case 10:
          r0ab = 2.7637999057769775;
          break;
        case 11:
          r0ab = 3.6085000038146973;
          break;
        case 12:
          r0ab = 3.1214001178741455;
          break;
        case 13:
          r0ab = 3.4005999565124512;
          break;
        case 14:
          r0ab = 3.4793000221252441;
          break;
        case 15:
          r0ab = 3.5146999359130859;
          break;
        case 16:
          r0ab = 3.3805999755859375;
          break;
        case 17:
          r0ab = 3.2355999946594238;
          break;
        case 18:
          r0ab = 3.1335000991821289;
          break;
        case 19:
          r0ab = 3.9144001007080078;
          break;
        case 20:
          r0ab = 3.4182999134063721;
          break;
        case 21:
          r0ab = 3.3368999958038330;
          break;
        case 22:
          r0ab = 3.2802999019622803;
          break;
        case 23:
          r0ab = 3.2678999900817871;
          break;
        case 24:
          r0ab = 3.4870998859405518;
          break;
        case 25:
          r0ab = 3.1714000701904297;
          break;
        case 26:
          r0ab = 3.1521000862121582;
          break;
        case 27:
          r0ab = 3.1101000308990479;
          break;
        case 28:
          r0ab = 3.0843000411987305;
          break;
        case 29:
          r0ab = 3.0669999122619629;
          break;
        case 30:
          r0ab = 3.0539000034332275;
          break;
        case 31:
          r0ab = 3.3889999389648438;
          break;
        case 32:
          r0ab = 3.5085999965667725;
          break;
        case 33:
          r0ab = 3.5894999504089355;
          break;
        case 34:
          r0ab = 3.4783000946044922;
          break;
        case 35:
          r0ab = 3.3484001159667969;
          break;
        case 36:
          r0ab = 3.2558999061584473;
          break;
        case 37:
          r0ab = 4.0422000885009766;
          break;
        case 38:
          r0ab = 3.5966999530792236;
          break;
        case 39:
          r0ab = 3.5113000869750977;
          break;
        case 40:
          r0ab = 3.4576001167297363;
          break;
        case 41:
          r0ab = 3.6593999862670898;
          break;
        case 42:
          r0ab = 3.6312999725341797;
          break;
        case 43:
          r0ab = 3.5690000057220459;
          break;
      }
      break;
    case 44:
      switch (atmA) {
        case 1:
          r0ab = 2.8578000068664551;
          break;
        case 2:
          r0ab = 2.6333999633789062;
          break;
        case 3:
          r0ab = 3.4672999382019043;
          break;
        case 4:
          r0ab = 2.9244999885559082;
          break;
        case 5:
          r0ab = 3.0731999874114990;
          break;
        case 6:
          r0ab = 3.2434999942779541;
          break;
        case 7:
          r0ab = 3.0337998867034912;
          break;
        case 8:
          r0ab = 2.9461998939514160;
          break;
        case 9:
          r0ab = 2.8143000602722168;
          break;
        case 10:
          r0ab = 2.7239999771118164;
          break;
        case 11:
          r0ab = 3.5831999778747559;
          break;
        case 12:
          r0ab = 3.0789000988006592;
          break;
        case 13:
          r0ab = 3.3617000579833984;
          break;
        case 14:
          r0ab = 3.4245998859405518;
          break;
        case 15:
          r0ab = 3.4505000114440918;
          break;
        case 16:
          r0ab = 3.3443000316619873;
          break;
        case 17:
          r0ab = 3.1963999271392822;
          break;
        case 18:
          r0ab = 3.0913000106811523;
          break;
        case 19:
          r0ab = 3.8921000957489014;
          break;
        case 20:
          r0ab = 3.3712999820709229;
          break;
        case 21:
          r0ab = 3.2873001098632812;
          break;
        case 22:
          r0ab = 3.2281000614166260;
          break;
        case 23:
          r0ab = 3.2165000438690186;
          break;
        case 24:
          r0ab = 3.4386000633239746;
          break;
        case 25:
          r0ab = 3.1164000034332275;
          break;
        case 26:
          r0ab = 3.1219999790191650;
          break;
        case 27:
          r0ab = 3.0761001110076904;
          break;
        case 28:
          r0ab = 3.0480000972747803;
          break;
        case 29:
          r0ab = 3.0295000076293945;
          break;
        case 30:
          r0ab = 3.0155000686645508;
          break;
        case 31:
          r0ab = 3.3494999408721924;
          break;
        case 32:
          r0ab = 3.4542999267578125;
          break;
        case 33:
          r0ab = 3.5260000228881836;
          break;
        case 34:
          r0ab = 3.4412999153137207;
          break;
        case 35:
          r0ab = 3.3085000514984131;
          break;
        case 36:
          r0ab = 3.2133998870849609;
          break;
        case 37:
          r0ab = 4.0170001983642578;
          break;
        case 38:
          r0ab = 3.5464000701904297;
          break;
        case 39:
          r0ab = 3.4586999416351318;
          break;
        case 40:
          r0ab = 3.4005999565124512;
          break;
        case 41:
          r0ab = 3.6026999950408936;
          break;
        case 42:
          r0ab = 3.5729999542236328;
          break;
        case 43:
          r0ab = 3.4944999217987061;
          break;
        case 44:
          r0ab = 3.4623000621795654;
          break;
      }
      break;
    case 45:
      switch (atmA) {
        case 1:
          r0ab = 2.8239998817443848;
          break;
        case 2:
          r0ab = 2.5959999561309814;
          break;
        case 3:
          r0ab = 3.4635000228881836;
          break;
        case 4:
          r0ab = 2.9031999111175537;
          break;
        case 5:
          r0ab = 3.0431001186370850;
          break;
        case 6:
          r0ab = 3.2114999294281006;
          break;
        case 7:
          r0ab = 2.9892001152038574;
          break;
        case 8:
          r0ab = 2.9147999286651611;
          break;
        case 9:
          r0ab = 2.7801001071929932;
          break;
        case 10:
          r0ab = 2.6872999668121338;
          break;
        case 11:
          r0ab = 3.5776000022888184;
          break;
        case 12:
          r0ab = 3.0567998886108398;
          break;
        case 13:
          r0ab = 3.3433001041412354;
          break;
        case 14:
          r0ab = 3.3949000835418701;
          break;
        case 15:
          r0ab = 3.4131999015808105;
          break;
        case 16:
          r0ab = 3.3115999698638916;
          break;
        case 17:
          r0ab = 3.1616001129150391;
          break;
        case 18:
          r0ab = 3.0548000335693359;
          break;
        case 19:
          r0ab = 3.8859000205993652;
          break;
        case 20:
          r0ab = 3.3719000816345215;
          break;
        case 21:
          r0ab = 3.2916998863220215;
          break;
        case 22:
          r0ab = 3.2344999313354492;
          break;
        case 23:
          r0ab = 3.2274000644683838;
          break;
        case 24:
          r0ab = 3.4170999526977539;
          break;
        case 25:
          r0ab = 3.1293001174926758;
          break;
        case 26:
          r0ab = 3.0566999912261963;
          break;
        case 27:
          r0ab = 3.0564999580383301;
          break;
        case 28:
          r0ab = 3.0274000167846680;
          break;
        case 29:
          r0ab = 3.0086998939514160;
          break;
        case 30:
          r0ab = 2.9939000606536865;
          break;
        case 31:
          r0ab = 3.3292999267578125;
          break;
        case 32:
          r0ab = 3.4249000549316406;
          break;
        case 33:
          r0ab = 3.4902000427246094;
          break;
        case 34:
          r0ab = 3.4091000556945801;
          break;
        case 35:
          r0ab = 3.2743999958038330;
          break;
        case 36:
          r0ab = 3.1775999069213867;
          break;
        case 37:
          r0ab = 4.0078001022338867;
          break;
        case 38:
          r0ab = 3.5374000072479248;
          break;
        case 39:
          r0ab = 3.4537000656127930;
          break;
        case 40:
          r0ab = 3.3956000804901123;
          break;
        case 41:
          r0ab = 3.5747001171112061;
          break;
        case 42:
          r0ab = 3.5429999828338623;
          break;
        case 43:
          r0ab = 3.4521999359130859;
          break;
        case 44:
          r0ab = 3.4159998893737793;
          break;
        case 45:
          r0ab = 3.3975000381469727;
          break;
      }
      break;
    case 46:
      switch (atmA) {
        case 1:
          r0ab = 2.8004000186920166;
          break;
        case 2:
          r0ab = 2.5620999336242676;
          break;
        case 3:
          r0ab = 3.4616999626159668;
          break;
        case 4:
          r0ab = 2.9154000282287598;
          break;
        case 5:
          r0ab = 3.0202999114990234;
          break;
        case 6:
          r0ab = 3.1875000000000000;
          break;
        case 7:
          r0ab = 2.9547998905181885;
          break;
        case 8:
          r0ab = 2.8038001060485840;
          break;
        case 9:
          r0ab = 2.7472000122070312;
          break;
        case 10:
          r0ab = 2.6530001163482666;
          break;
        case 11:
          r0ab = 3.5736000537872314;
          break;
        case 12:
          r0ab = 3.0583999156951904;
          break;
        case 13:
          r0ab = 3.3303999900817871;
          break;
        case 14:
          r0ab = 3.3747999668121338;
          break;
        case 15:
          r0ab = 3.3870999813079834;
          break;
        case 16:
          r0ab = 3.2028000354766846;
          break;
        case 17:
          r0ab = 3.1296000480651855;
          break;
        case 18:
          r0ab = 3.0213999748229980;
          break;
        case 19:
          r0ab = 3.8796000480651855;
          break;
        case 20:
          r0ab = 3.3336999416351318;
          break;
        case 21:
          r0ab = 3.2492001056671143;
          break;
        case 22:
          r0ab = 3.1882998943328857;
          break;
        case 23:
          r0ab = 3.1802000999450684;
          break;
        case 24:
          r0ab = 3.4049999713897705;
          break;
        case 25:
          r0ab = 3.0755999088287354;
          break;
        case 26:
          r0ab = 3.0478000640869141;
          break;
        case 27:
          r0ab = 3.0322000980377197;
          break;
        case 28:
          r0ab = 3.0322999954223633;
          break;
        case 29:
          r0ab = 3.0162999629974365;
          break;
        case 30:
          r0ab = 3.0018999576568604;
          break;
        case 31:
          r0ab = 3.3145000934600830;
          break;
        case 32:
          r0ab = 3.4049999713897705;
          break;
        case 33:
          r0ab = 3.4656000137329102;
          break;
        case 34:
          r0ab = 3.3020999431610107;
          break;
        case 35:
          r0ab = 3.2432999610900879;
          break;
        case 36:
          r0ab = 3.1452999114990234;
          break;
        case 37:
          r0ab = 3.9990999698638916;
          break;
        case 38:
          r0ab = 3.5016999244689941;
          break;
        case 39:
          r0ab = 3.4140999317169189;
          break;
        case 40:
          r0ab = 3.3519999980926514;
          break;
        case 41:
          r0ab = 3.5583000183105469;
          break;
        case 42:
          r0ab = 3.5250999927520752;
          break;
        case 43:
          r0ab = 3.4242999553680420;
          break;
        case 44:
          r0ab = 3.3850998878479004;
          break;
        case 45:
          r0ab = 3.3661999702453613;
          break;
        case 46:
          r0ab = 3.3524999618530273;
          break;
      }
      break;
    case 47:
      switch (atmA) {
        case 1:
          r0ab = 2.7846000194549561;
          break;
        case 2:
          r0ab = 2.5323998928070068;
          break;
        case 3:
          r0ab = 3.4651999473571777;
          break;
        case 4:
          r0ab = 2.8759000301361084;
          break;
        case 5:
          r0ab = 3.0051000118255615;
          break;
        case 6:
          r0ab = 3.1691999435424805;
          break;
        case 7:
          r0ab = 2.9272999763488770;
          break;
        case 8:
          r0ab = 2.7614998817443848;
          break;
        case 9:
          r0ab = 2.7163999080657959;
          break;
        case 10:
          r0ab = 2.6212000846862793;
          break;
        case 11:
          r0ab = 3.5743999481201172;
          break;
        case 12:
          r0ab = 3.0274999141693115;
          break;
        case 13:
          r0ab = 3.3248999118804932;
          break;
        case 14:
          r0ab = 3.3626999855041504;
          break;
        case 15:
          r0ab = 3.3685998916625977;
          break;
        case 16:
          r0ab = 3.1668999195098877;
          break;
        case 17:
          r0ab = 3.0583999156951904;
          break;
        case 18:
          r0ab = 2.9914999008178711;
          break;
        case 19:
          r0ab = 3.8773000240325928;
          break;
        case 20:
          r0ab = 3.3099000453948975;
          break;
        case 21:
          r0ab = 3.2230999469757080;
          break;
        case 22:
          r0ab = 3.1600000858306885;
          break;
        case 23:
          r0ab = 3.1519999504089355;
          break;
        case 24:
          r0ab = 3.4022998809814453;
          break;
        case 25:
          r0ab = 3.0425999164581299;
          break;
        case 26:
          r0ab = 3.0099000930786133;
          break;
        case 27:
          r0ab = 2.9920001029968262;
          break;
        case 28:
          r0ab = 2.9809000492095947;
          break;
        case 29:
          r0ab = 2.9800000190734863;
          break;
        case 30:
          r0ab = 2.9646000862121582;
          break;
        case 31:
          r0ab = 3.3067998886108398;
          break;
        case 32:
          r0ab = 3.3929998874664307;
          break;
        case 33:
          r0ab = 3.4486000537872314;
          break;
        case 34:
          r0ab = 3.2681999206542969;
          break;
        case 35:
          r0ab = 3.1728999614715576;
          break;
        case 36:
          r0ab = 3.1168000698089600;
          break;
        case 37:
          r0ab = 3.9951999187469482;
          break;
        case 38:
          r0ab = 3.4795999526977539;
          break;
        case 39:
          r0ab = 3.3901000022888184;
          break;
        case 40:
          r0ab = 3.3255000114440918;
          break;
        case 41:
          r0ab = 3.5529999732971191;
          break;
        case 42:
          r0ab = 3.5183000564575195;
          break;
        case 43:
          r0ab = 3.4096999168395996;
          break;
        case 44:
          r0ab = 3.3682999610900879;
          break;
        case 45:
          r0ab = 3.3492000102996826;
          break;
        case 46:
          r0ab = 3.3359999656677246;
          break;
        case 47:
          r0ab = 3.3308000564575195;
          break;
      }
      break;
    case 48:
      switch (atmA) {
        case 1:
          r0ab = 2.5423998832702637;
          break;
        case 2:
          r0ab = 2.6600999832153320;
          break;
        case 3:
          r0ab = 3.2555000782012939;
          break;
        case 4:
          r0ab = 3.2806999683380127;
          break;
        case 5:
          r0ab = 3.1384000778198242;
          break;
        case 6:
          r0ab = 3.1737000942230225;
          break;
        case 7:
          r0ab = 2.9396998882293701;
          break;
        case 8:
          r0ab = 2.8429000377655029;
          break;
        case 9:
          r0ab = 2.8492000102996826;
          break;
        case 10:
          r0ab = 2.7225000858306885;
          break;
        case 11:
          r0ab = 3.3875000476837158;
          break;
        case 12:
          r0ab = 3.4909999370574951;
          break;
        case 13:
          r0ab = 3.4519999027252197;
          break;
        case 14:
          r0ab = 3.3608000278472900;
          break;
        case 15:
          r0ab = 3.3036000728607178;
          break;
        case 16:
          r0ab = 3.2344999313354492;
          break;
        case 17:
          r0ab = 3.2999000549316406;
          break;
        case 18:
          r0ab = 3.1486999988555908;
          break;
        case 19:
          r0ab = 3.7409000396728516;
          break;
        case 20:
          r0ab = 3.8392000198364258;
          break;
        case 21:
          r0ab = 3.7147998809814453;
          break;
        case 22:
          r0ab = 3.6438999176025391;
          break;
        case 23:
          r0ab = 3.6182000637054443;
          break;
        case 24:
          r0ab = 3.1752998828887939;
          break;
        case 25:
          r0ab = 3.5209999084472656;
          break;
        case 26:
          r0ab = 3.4639000892639160;
          break;
        case 27:
          r0ab = 3.4265000820159912;
          break;
        case 28:
          r0ab = 3.4075000286102295;
          break;
        case 29:
          r0ab = 3.3828001022338867;
          break;
        case 30:
          r0ab = 3.3473999500274658;
          break;
        case 31:
          r0ab = 3.4070999622344971;
          break;
        case 32:
          r0ab = 3.3754000663757324;
          break;
        case 33:
          r0ab = 3.3645999431610107;
          break;
        case 34:
          r0ab = 3.3308000564575195;
          break;
        case 35:
          r0ab = 3.4393000602722168;
          break;
        case 36:
          r0ab = 3.2992999553680420;
          break;
        case 37:
          r0ab = 3.8768000602722168;
          break;
        case 38:
          r0ab = 3.9890999794006348;
          break;
        case 39:
          r0ab = 3.8310000896453857;
          break;
        case 40:
          r0ab = 3.7483000755310059;
          break;
        case 41:
          r0ab = 3.3417000770568848;
          break;
        case 42:
          r0ab = 3.3018999099731445;
          break;
        case 43:
          r0ab = 3.2249999046325684;
          break;
        case 44:
          r0ab = 3.1831998825073242;
          break;
        case 45:
          r0ab = 3.1577999591827393;
          break;
        case 46:
          r0ab = 3.1563999652862549;
          break;
        case 47:
          r0ab = 3.1224000453948975;
          break;
        case 48:
          r0ab = 3.4619998931884766;
          break;
      }
      break;
    case 49:
      switch (atmA) {
        case 1:
          r0ab = 2.9742999076843262;
          break;
        case 2:
          r0ab = 2.8057999610900879;
          break;
        case 3:
          r0ab = 3.4830000400543213;
          break;
        case 4:
          r0ab = 3.3473999500274658;
          break;
        case 5:
          r0ab = 3.6863000392913818;
          break;
        case 6:
          r0ab = 3.3617000579833984;
          break;
        case 7:
          r0ab = 3.1607999801635742;
          break;
        case 8:
          r0ab = 3.0069000720977783;
          break;
        case 9:
          r0ab = 2.9639999866485596;
          break;
        case 10:
          r0ab = 2.8427000045776367;
          break;
        case 11:
          r0ab = 3.5885000228881836;
          break;
        case 12:
          r0ab = 3.5218999385833740;
          break;
        case 13:
          r0ab = 4.1314001083374023;
          break;
        case 14:
          r0ab = 3.8120000362396240;
          break;
        case 15:
          r0ab = 3.6015000343322754;
          break;
        case 16:
          r0ab = 3.4502000808715820;
          break;
        case 17:
          r0ab = 3.3498001098632812;
          break;
        case 18:
          r0ab = 3.2776999473571777;
          break;
        case 19:
          r0ab = 3.8635001182556152;
          break;
        case 20:
          r0ab = 3.8231999874114990;
          break;
        case 21:
          r0ab = 3.8485999107360840;
          break;
        case 22:
          r0ab = 3.7214999198913574;
          break;
        case 23:
          r0ab = 3.6486999988555908;
          break;
        case 24:
          r0ab = 3.4723999500274658;
          break;
        case 25:
          r0ab = 3.5627000331878662;
          break;
        case 26:
          r0ab = 3.5086998939514160;
          break;
        case 27:
          r0ab = 3.4756999015808105;
          break;
        case 28:
          r0ab = 3.4516999721527100;
          break;
        case 29:
          r0ab = 3.4423000812530518;
          break;
        case 30:
          r0ab = 3.4138998985290527;
          break;
        case 31:
          r0ab = 4.1027998924255371;
          break;
        case 32:
          r0ab = 3.8387999534606934;
          break;
        case 33:
          r0ab = 3.6744999885559082;
          break;
        case 34:
          r0ab = 3.5562000274658203;
          break;
        case 35:
          r0ab = 3.4806001186370850;
          break;
        case 36:
          r0ab = 3.4272000789642334;
          break;
        case 37:
          r0ab = 4.0181999206542969;
          break;
        case 38:
          r0ab = 3.9990999698638916;
          break;
        case 39:
          r0ab = 4.0006999969482422;
          break;
        case 40:
          r0ab = 3.9282000064849854;
          break;
        case 41:
          r0ab = 3.7237999439239502;
          break;
        case 42:
          r0ab = 3.6498000621795654;
          break;
        case 43:
          r0ab = 3.5604999065399170;
          break;
        case 44:
          r0ab = 3.5211000442504883;
          break;
        case 45:
          r0ab = 3.5009000301361084;
          break;
        case 46:
          r0ab = 3.4858999252319336;
          break;
        case 47:
          r0ab = 3.4784998893737793;
          break;
        case 48:
          r0ab = 3.5620999336242676;
          break;
        case 49:
          r0ab = 4.2623000144958496;
          break;
      }
      break;
    case 50:
      switch (atmA) {
        case 1:
          r0ab = 3.0775001049041748;
          break;
        case 2:
          r0ab = 2.8275001049041748;
          break;
        case 3:
          r0ab = 4.0180997848510742;
          break;
        case 4:
          r0ab = 3.3385000228881836;
          break;
        case 5:
          r0ab = 3.5378999710083008;
          break;
        case 6:
          r0ab = 3.5035998821258545;
          break;
        case 7:
          r0ab = 3.2588999271392822;
          break;
        case 8:
          r0ab = 3.0803999900817871;
          break;
        case 9:
          r0ab = 3.0093998908996582;
          break;
        case 10:
          r0ab = 2.9003000259399414;
          break;
        case 11:
          r0ab = 4.0869002342224121;
          break;
        case 12:
          r0ab = 3.5088000297546387;
          break;
        case 13:
          r0ab = 3.9105000495910645;
          break;
        case 14:
          r0ab = 3.9832999706268311;
          break;
        case 15:
          r0ab = 3.7176001071929932;
          break;
        case 16:
          r0ab = 3.5322999954223633;
          break;
        case 17:
          r0ab = 3.4102001190185547;
          break;
        case 18:
          r0ab = 3.3227000236511230;
          break;
        case 19:
          r0ab = 4.2701997756958008;
          break;
        case 20:
          r0ab = 4.0887999534606934;
          break;
        case 21:
          r0ab = 3.7560000419616699;
          break;
        case 22:
          r0ab = 3.7686998844146729;
          break;
        case 23:
          r0ab = 3.6681001186370850;
          break;
        case 24:
          r0ab = 3.6405000686645508;
          break;
        case 25:
          r0ab = 3.5569000244140625;
          break;
        case 26:
          r0ab = 3.4990000724792480;
          break;
        case 27:
          r0ab = 3.4658999443054199;
          break;
        case 28:
          r0ab = 3.4433000087738037;
          break;
        case 29:
          r0ab = 3.4330000877380371;
          break;
        case 30:
          r0ab = 3.4091999530792236;
          break;
        case 31:
          r0ab = 3.8866999149322510;
          break;
        case 32:
          r0ab = 4.0190000534057617;
          break;
        case 33:
          r0ab = 3.7960999011993408;
          break;
        case 34:
          r0ab = 3.6412000656127930;
          break;
        case 35:
          r0ab = 3.5404999256134033;
          break;
        case 36:
          r0ab = 3.4681000709533691;
          break;
        case 37:
          r0ab = 4.3537998199462891;
          break;
        case 38:
          r0ab = 4.2136001586914062;
          break;
        case 39:
          r0ab = 3.9381000995635986;
          break;
        case 40:
          r0ab = 3.8912000656127930;
          break;
        case 41:
          r0ab = 3.9681000709533691;
          break;
        case 42:
          r0ab = 3.7908999919891357;
          break;
        case 43:
          r0ab = 3.6774001121520996;
          break;
        case 44:
          r0ab = 3.6261999607086182;
          break;
        case 45:
          r0ab = 3.5999000072479248;
          break;
        case 46:
          r0ab = 3.5822999477386475;
          break;
        case 47:
          r0ab = 3.5727000236511230;
          break;
        case 48:
          r0ab = 3.5418999195098877;
          break;
        case 49:
          r0ab = 4.0244998931884766;
          break;
        case 50:
          r0ab = 4.1873998641967773;
          break;
      }
      break;
    case 51:
      switch (atmA) {
        case 1:
          r0ab = 3.0892999172210693;
          break;
        case 2:
          r0ab = 2.7916998863220215;
          break;
        case 3:
          r0ab = 3.7262001037597656;
          break;
        case 4:
          r0ab = 3.3517999649047852;
          break;
        case 5:
          r0ab = 3.4240999221801758;
          break;
        case 6:
          r0ab = 3.5432999134063721;
          break;
        case 7:
          r0ab = 3.2772998809814453;
          break;
        case 8:
          r0ab = 3.0889999866485596;
          break;
        case 9:
          r0ab = 2.9774999618530273;
          break;
        case 10:
          r0ab = 2.9010000228881836;
          break;
        case 11:
          r0ab = 3.8048000335693359;
          break;
        case 12:
          r0ab = 3.5362000465393066;
          break;
        case 13:
          r0ab = 3.7746000289916992;
          break;
        case 14:
          r0ab = 3.7911000251770020;
          break;
        case 15:
          r0ab = 3.7511000633239746;
          break;
        case 16:
          r0ab = 3.5494999885559082;
          break;
        case 17:
          r0ab = 3.4149000644683838;
          break;
        case 18:
          r0ab = 3.3176999092102051;
          break;
        case 19:
          r0ab = 4.0128998756408691;
          break;
        case 20:
          r0ab = 3.8369998931884766;
          break;
        case 21:
          r0ab = 3.7739000320434570;
          break;
        case 22:
          r0ab = 3.7125000953674316;
          break;
        case 23:
          r0ab = 3.7151999473571777;
          break;
        case 24:
          r0ab = 3.7701001167297363;
          break;
        case 25:
          r0ab = 3.5813000202178955;
          break;
        case 26:
          r0ab = 3.5186998844146729;
          break;
        case 27:
          r0ab = 3.4835000038146973;
          break;
        case 28:
          r0ab = 3.4595000743865967;
          break;
        case 29:
          r0ab = 3.4439001083374023;
          break;
        case 30:
          r0ab = 3.4242000579833984;
          break;
        case 31:
          r0ab = 3.7476000785827637;
          break;
        case 32:
          r0ab = 3.8238999843597412;
          break;
        case 33:
          r0ab = 3.8345999717712402;
          break;
        case 34:
          r0ab = 3.6626999378204346;
          break;
        case 35:
          r0ab = 3.5478999614715576;
          break;
        case 36:
          r0ab = 3.4639000892639160;
          break;
        case 37:
          r0ab = 4.1026000976562500;
          break;
        case 38:
          r0ab = 3.9732999801635742;
          break;
        case 39:
          r0ab = 3.9291999340057373;
          break;
        case 40:
          r0ab = 3.8666999340057373;
          break;
        case 41:
          r0ab = 3.9512999057769775;
          break;
        case 42:
          r0ab = 3.8959000110626221;
          break;
        case 43:
          r0ab = 3.7697999477386475;
          break;
        case 44:
          r0ab = 3.7088999748229980;
          break;
        case 45:
          r0ab = 3.6765000820159912;
          break;
        case 46:
          r0ab = 3.6547999382019043;
          break;
        case 47:
          r0ab = 3.6408998966217041;
          break;
        case 48:
          r0ab = 3.5397999286651611;
          break;
        case 49:
          r0ab = 3.8759000301361084;
          break;
        case 50:
          r0ab = 3.9804000854492188;
          break;
        case 51:
          r0ab = 4.0149998664855957;
          break;
      }
      break;
    case 52:
      switch (atmA) {
        case 1:
          r0ab = 2.9091000556945801;
          break;
        case 2:
          r0ab = 2.7637999057769775;
          break;
        case 3:
          r0ab = 3.5065999031066895;
          break;
        case 4:
          r0ab = 3.3376998901367188;
          break;
        case 5:
          r0ab = 3.3480999469757080;
          break;
        case 6:
          r0ab = 3.2632999420166016;
          break;
        case 7:
          r0ab = 3.1809999942779541;
          break;
        case 8:
          r0ab = 3.1428000926971436;
          break;
        case 9:
          r0ab = 2.9872000217437744;
          break;
        case 10:
          r0ab = 2.8836998939514160;
          break;
        case 11:
          r0ab = 3.5929000377655029;
          break;
        case 12:
          r0ab = 3.5183000564575195;
          break;
        case 13:
          r0ab = 3.6728999614715576;
          break;
        case 14:
          r0ab = 3.6596000194549561;
          break;
        case 15:
          r0ab = 3.6082000732421875;
          break;
        case 16:
          r0ab = 3.5927000045776367;
          break;
        case 17:
          r0ab = 3.4223999977111816;
          break;
        case 18:
          r0ab = 3.2997000217437744;
          break;
        case 19:
          r0ab = 3.8190000057220459;
          break;
        case 20:
          r0ab = 4.1865000724792480;
          break;
        case 21:
          r0ab = 4.1114001274108887;
          break;
        case 22:
          r0ab = 4.0539999008178711;
          break;
        case 23:
          r0ab = 3.6324999332427979;
          break;
        case 24:
          r0ab = 3.5697000026702881;
          break;
        case 25:
          r0ab = 3.5560998916625977;
          break;
        case 26:
          r0ab = 3.5258998870849609;
          break;
        case 27:
          r0ab = 3.4900999069213867;
          break;
        case 28:
          r0ab = 3.4551999568939209;
          break;
        case 29:
          r0ab = 3.4314999580383301;
          break;
        case 30:
          r0ab = 3.4091000556945801;
          break;
        case 31:
          r0ab = 3.6438000202178955;
          break;
        case 32:
          r0ab = 3.6879000663757324;
          break;
        case 33:
          r0ab = 3.6831998825073242;
          break;
        case 34:
          r0ab = 3.7042999267578125;
          break;
        case 35:
          r0ab = 3.5557000637054443;
          break;
        case 36:
          r0ab = 3.4465999603271484;
          break;
        case 37:
          r0ab = 3.9203000068664551;
          break;
        case 38:
          r0ab = 4.2919001579284668;
          break;
        case 39:
          r0ab = 4.2196002006530762;
          break;
        case 40:
          r0ab = 4.1542000770568848;
          break;
        case 41:
          r0ab = 3.7572999000549316;
          break;
        case 42:
          r0ab = 3.7039000988006592;
          break;
        case 43:
          r0ab = 3.6545999050140381;
          break;
        case 44:
          r0ab = 3.6150999069213867;
          break;
        case 45:
          r0ab = 3.5292999744415283;
          break;
        case 46:
          r0ab = 3.4848999977111816;
          break;
        case 47:
          r0ab = 3.4551999568939209;
          break;
        case 48:
          r0ab = 3.5192000865936279;
          break;
        case 49:
          r0ab = 3.7672998905181885;
          break;
        case 50:
          r0ab = 3.8359000682830811;
          break;
        case 51:
          r0ab = 3.8524999618530273;
          break;
        case 52:
          r0ab = 3.8901000022888184;
          break;
      }
      break;
    case 53:
      switch (atmA) {
        case 1:
          r0ab = 2.7806000709533691;
          break;
        case 2:
          r0ab = 2.7209000587463379;
          break;
        case 3:
          r0ab = 3.3812000751495361;
          break;
        case 4:
          r0ab = 3.4958000183105469;
          break;
        case 5:
          r0ab = 3.2913000583648682;
          break;
        case 6:
          r0ab = 3.1888000965118408;
          break;
        case 7:
          r0ab = 3.0989999771118164;
          break;
        case 8:
          r0ab = 3.0394001007080078;
          break;
        case 9:
          r0ab = 2.9788999557495117;
          break;
        case 10:
          r0ab = 2.8582000732421875;
          break;
        case 11:
          r0ab = 3.4716000556945801;
          break;
        case 12:
          r0ab = 3.6882998943328857;
          break;
        case 13:
          r0ab = 3.6105000972747803;
          break;
        case 14:
          r0ab = 3.5703999996185303;
          break;
        case 15:
          r0ab = 3.5058999061584473;
          break;
        case 16:
          r0ab = 3.4618999958038330;
          break;
        case 17:
          r0ab = 3.4138000011444092;
          break;
        case 18:
          r0ab = 3.2741999626159668;
          break;
        case 19:
          r0ab = 3.7079999446868896;
          break;
        case 20:
          r0ab = 3.9772999286651611;
          break;
        case 21:
          r0ab = 3.9010000228881836;
          break;
        case 22:
          r0ab = 3.8408999443054199;
          break;
        case 23:
          r0ab = 3.7943999767303467;
          break;
        case 24:
          r0ab = 3.4465000629425049;
          break;
        case 25:
          r0ab = 3.7235000133514404;
          break;
        case 26:
          r0ab = 3.6807999610900879;
          break;
        case 27:
          r0ab = 3.6452999114990234;
          break;
        case 28:
          r0ab = 3.6168000698089600;
          break;
        case 29:
          r0ab = 3.5843999385833740;
          break;
        case 30:
          r0ab = 3.5576000213623047;
          break;
        case 31:
          r0ab = 3.5771999359130859;
          break;
        case 32:
          r0ab = 3.5959000587463379;
          break;
        case 33:
          r0ab = 3.5768001079559326;
          break;
        case 34:
          r0ab = 3.5678000450134277;
          break;
        case 35:
          r0ab = 3.5485999584197998;
          break;
        case 36:
          r0ab = 3.4228000640869141;
          break;
        case 37:
          r0ab = 3.8106999397277832;
          break;
        case 38:
          r0ab = 4.0865998268127441;
          break;
        case 39:
          r0ab = 4.0169000625610352;
          break;
        case 40:
          r0ab = 3.9475998878479004;
          break;
        case 41:
          r0ab = 3.6357998847961426;
          break;
        case 42:
          r0ab = 3.5799999237060547;
          break;
        case 43:
          r0ab = 3.5260000228881836;
          break;
        case 44:
          r0ab = 3.4837999343872070;
          break;
        case 45:
          r0ab = 3.4500999450683594;
          break;
        case 46:
          r0ab = 3.4203999042510986;
          break;
        case 47:
          r0ab = 3.3552999496459961;
          break;
        case 48:
          r0ab = 3.6486999988555908;
          break;
        case 49:
          r0ab = 3.6972999572753906;
          break;
        case 50:
          r0ab = 3.7397999763488770;
          break;
        case 51:
          r0ab = 3.7404999732971191;
          break;
        case 52:
          r0ab = 3.7458999156951904;
          break;
        case 53:
          r0ab = 3.7379999160766602;
          break;
      }
      break;
    case 54:
      switch (atmA) {
        case 1:
          r0ab = 2.6847999095916748;
          break;
        case 2:
          r0ab = 2.6740000247955322;
          break;
        case 3:
          r0ab = 3.2925000190734863;
          break;
        case 4:
          r0ab = 3.3385999202728271;
          break;
        case 5:
          r0ab = 3.2472999095916748;
          break;
        case 6:
          r0ab = 3.1284000873565674;
          break;
        case 7:
          r0ab = 3.0301001071929932;
          break;
        case 8:
          r0ab = 2.9530999660491943;
          break;
        case 9:
          r0ab = 2.9602000713348389;
          break;
        case 10:
          r0ab = 2.8271999359130859;
          break;
        case 11:
          r0ab = 3.3829998970031738;
          break;
        case 12:
          r0ab = 3.5357999801635742;
          break;
        case 13:
          r0ab = 3.5671999454498291;
          break;
        case 14:
          r0ab = 3.5048999786376953;
          break;
        case 15:
          r0ab = 3.4284000396728516;
          break;
        case 16:
          r0ab = 3.3620998859405518;
          break;
        case 17:
          r0ab = 3.3001000881195068;
          break;
        case 18:
          r0ab = 3.2451000213623047;
          break;
        case 19:
          r0ab = 3.6208999156951904;
          break;
        case 20:
          r0ab = 3.8299000263214111;
          break;
        case 21:
          r0ab = 3.7543001174926758;
          break;
        case 22:
          r0ab = 3.6919999122619629;
          break;
        case 23:
          r0ab = 3.6435999870300293;
          break;
        case 24:
          r0ab = 3.3598001003265381;
          break;
        case 25:
          r0ab = 3.5701000690460205;
          break;
        case 26:
          r0ab = 3.5265998840332031;
          break;
        case 27:
          r0ab = 3.4904000759124756;
          break;
        case 28:
          r0ab = 3.4590001106262207;
          break;
        case 29:
          r0ab = 3.4363999366760254;
          break;
        case 30:
          r0ab = 3.4077000617980957;
          break;
        case 31:
          r0ab = 3.5287001132965088;
          break;
        case 32:
          r0ab = 3.5280001163482666;
          break;
        case 33:
          r0ab = 3.4969000816345215;
          break;
        case 34:
          r0ab = 3.4649999141693115;
          break;
        case 35:
          r0ab = 3.4303998947143555;
          break;
        case 36:
          r0ab = 3.3963000774383545;
          break;
        case 37:
          r0ab = 3.7228999137878418;
          break;
        case 38:
          r0ab = 3.9402000904083252;
          break;
        case 39:
          r0ab = 3.8752999305725098;
          break;
        case 40:
          r0ab = 3.8034999370574951;
          break;
        case 41:
          r0ab = 3.5499000549316406;
          break;
        case 42:
          r0ab = 3.4913001060485840;
          break;
        case 43:
          r0ab = 3.4319000244140625;
          break;
        case 44:
          r0ab = 3.3873000144958496;
          break;
        case 45:
          r0ab = 3.3519999980926514;
          break;
        case 46:
          r0ab = 3.3208999633789062;
          break;
        case 47:
          r0ab = 3.2948000431060791;
          break;
        case 48:
          r0ab = 3.5051999092102051;
          break;
        case 49:
          r0ab = 3.6465001106262207;
          break;
        case 50:
          r0ab = 3.6696000099182129;
          break;
        case 51:
          r0ab = 3.6577000617980957;
          break;
        case 52:
          r0ab = 3.6387999057769775;
          break;
        case 53:
          r0ab = 3.6142001152038574;
          break;
        case 54:
          r0ab = 3.5889000892639160;
          break;
      }
      break;
    case 55:
      switch (atmA) {
        case 1:
          r0ab = 3.3968000411987305;
          break;
        case 2:
          r0ab = 3.0122001171112061;
          break;
        case 3:
          r0ab = 4.2241001129150391;
          break;
        case 4:
          r0ab = 3.7887001037597656;
          break;
        case 5:
          r0ab = 4.0048999786376953;
          break;
        case 6:
          r0ab = 3.5383999347686768;
          break;
        case 7:
          r0ab = 3.2697999477386475;
          break;
        case 8:
          r0ab = 3.1082999706268311;
          break;
        case 9:
          r0ab = 2.9916999340057373;
          break;
        case 10:
          r0ab = 2.9056999683380127;
          break;
        case 11:
          r0ab = 4.3340001106262207;
          break;
        case 12:
          r0ab = 3.9900000095367432;
          break;
        case 13:
          r0ab = 4.6588001251220703;
          break;
        case 14:
          r0ab = 4.1277999877929688;
          break;
        case 15:
          r0ab = 3.8125000000000000;
          break;
        case 16:
          r0ab = 3.6189000606536865;
          break;
        case 17:
          r0ab = 3.4851000308990479;
          break;
        case 18:
          r0ab = 3.3859000205993652;
          break;
        case 19:
          r0ab = 4.6531000137329102;
          break;
        case 20:
          r0ab = 4.3133997917175293;
          break;
        case 21:
          r0ab = 4.2258000373840332;
          break;
        case 22:
          r0ab = 4.1308999061584473;
          break;
        case 23:
          r0ab = 4.0692000389099121;
          break;
        case 24:
          r0ab = 4.0943999290466309;
          break;
        case 25:
          r0ab = 3.9849998950958252;
          break;
        case 26:
          r0ab = 3.9416000843048096;
          break;
        case 27:
          r0ab = 3.9112000465393066;
          break;
        case 28:
          r0ab = 3.8873000144958496;
          break;
        case 29:
          r0ab = 3.8736000061035156;
          break;
        case 30:
          r0ab = 3.8473000526428223;
          break;
        case 31:
          r0ab = 4.6027002334594727;
          break;
        case 32:
          r0ab = 4.1538000106811523;
          break;
        case 33:
          r0ab = 3.8993999958038330;
          break;
        case 34:
          r0ab = 3.7418999671936035;
          break;
        case 35:
          r0ab = 3.6356000900268555;
          break;
        case 36:
          r0ab = 3.5548000335693359;
          break;
        case 37:
          r0ab = 4.8352999687194824;
          break;
        case 38:
          r0ab = 4.5412998199462891;
          break;
        case 39:
          r0ab = 4.3891000747680664;
          break;
        case 40:
          r0ab = 4.3415999412536621;
          break;
        case 41:
          r0ab = 4.3242998123168945;
          break;
        case 42:
          r0ab = 4.2753000259399414;
          break;
        case 43:
          r0ab = 4.2052998542785645;
          break;
        case 44:
          r0ab = 4.1789999008178711;
          break;
        case 45:
          r0ab = 4.1684999465942383;
          break;
        case 46:
          r0ab = 4.1585001945495605;
          break;
        case 47:
          r0ab = 4.1536002159118652;
          break;
        case 48:
          r0ab = 4.0578999519348145;
          break;
        case 49:
          r0ab = 4.1979999542236328;
          break;
        case 50:
          r0ab = 4.4563999176025391;
          break;
        case 51:
          r0ab = 4.2192001342773438;
          break;
        case 52:
          r0ab = 4.0528001785278320;
          break;
        case 53:
          r0ab = 3.9488999843597412;
          break;
        case 54:
          r0ab = 3.8642001152038574;
          break;
        case 55:
          r0ab = 5.0567002296447754;
          break;
      }
      break;
    case 56:
      switch (atmA) {
        case 1:
          r0ab = 3.0629999637603760;
          break;
        case 2:
          r0ab = 3.3271000385284424;
          break;
        case 3:
          r0ab = 4.0432000160217285;
          break;
        case 4:
          r0ab = 4.0046000480651855;
          break;
        case 5:
          r0ab = 4.1554999351501465;
          break;
        case 6:
          r0ab = 3.7425999641418457;
          break;
        case 7:
          r0ab = 3.5130000114440918;
          break;
        case 8:
          r0ab = 3.5174000263214111;
          break;
        case 9:
          r0ab = 3.2883999347686768;
          break;
        case 10:
          r0ab = 3.1377999782562256;
          break;
        case 11:
          r0ab = 4.1894001960754395;
          break;
        case 12:
          r0ab = 4.2321000099182129;
          break;
        case 13:
          r0ab = 4.1725001335144043;
          break;
        case 14:
          r0ab = 4.1833000183105469;
          break;
        case 15:
          r0ab = 3.8928999900817871;
          break;
        case 16:
          r0ab = 4.0543999671936035;
          break;
        case 17:
          r0ab = 3.8118000030517578;
          break;
        case 18:
          r0ab = 3.6414000988006592;
          break;
        case 19:
          r0ab = 4.6373000144958496;
          break;
        case 20:
          r0ab = 4.6268000602722168;
          break;
        case 21:
          r0ab = 4.4749999046325684;
          break;
        case 22:
          r0ab = 4.4134001731872559;
          break;
        case 23:
          r0ab = 4.3457999229431152;
          break;
        case 24:
          r0ab = 3.8582000732421875;
          break;
        case 25:
          r0ab = 4.2582998275756836;
          break;
        case 26:
          r0ab = 4.1897997856140137;
          break;
        case 27:
          r0ab = 4.1561999320983887;
          break;
        case 28:
          r0ab = 4.1191000938415527;
          break;
        case 29:
          r0ab = 4.1069002151489258;
          break;
        case 30:
          r0ab = 4.0638999938964844;
          break;
        case 31:
          r0ab = 4.1256999969482422;
          break;
        case 32:
          r0ab = 4.1974000930786133;
          break;
        case 33:
          r0ab = 3.9532001018524170;
          break;
        case 34:
          r0ab = 4.1793999671936035;
          break;
        case 35:
          r0ab = 3.9660000801086426;
          break;
        case 36:
          r0ab = 3.8129999637603760;
          break;
        case 37:
          r0ab = 4.8159999847412109;
          break;
        case 38:
          r0ab = 4.8271999359130859;
          break;
        case 39:
          r0ab = 4.6293997764587402;
          break;
        case 40:
          r0ab = 4.5840001106262207;
          break;
        case 41:
          r0ab = 4.0770001411437988;
          break;
        case 42:
          r0ab = 4.0088000297546387;
          break;
        case 43:
          r0ab = 3.9103000164031982;
          break;
        case 44:
          r0ab = 3.8536000251770020;
          break;
        case 45:
          r0ab = 3.8324000835418701;
          break;
        case 46:
          r0ab = 3.7994999885559082;
          break;
        case 47:
          r0ab = 3.7825999259948730;
          break;
        case 48:
          r0ab = 4.2294001579284668;
          break;
        case 49:
          r0ab = 4.3379998207092285;
          break;
        case 50:
          r0ab = 4.4352002143859863;
          break;
        case 51:
          r0ab = 4.1932997703552246;
          break;
        case 52:
          r0ab = 4.4580001831054688;
          break;
        case 53:
          r0ab = 4.2554001808166504;
          break;
        case 54:
          r0ab = 4.1072001457214355;
          break;
        case 55:
          r0ab = 5.0454001426696777;
          break;
        case 56:
          r0ab = 5.1813998222351074;
          break;
      }
      break;
    case 57:
      switch (atmA) {
        case 1:
          r0ab = 3.0631999969482422;
          break;
        case 2:
          r0ab = 3.2662000656127930;
          break;
        case 3:
          r0ab = 3.6431999206542969;
          break;
        case 4:
          r0ab = 3.8087999820709229;
          break;
        case 5:
          r0ab = 3.7909998893737793;
          break;
        case 6:
          r0ab = 3.7381000518798828;
          break;
        case 7:
          r0ab = 3.5092999935150146;
          break;
        case 8:
          r0ab = 3.5155000686645508;
          break;
        case 9:
          r0ab = 3.3046998977661133;
          break;
        case 10:
          r0ab = 3.1681001186370850;
          break;
        case 11:
          r0ab = 3.7871000766754150;
          break;
        case 12:
          r0ab = 3.9923999309539795;
          break;
        case 13:
          r0ab = 4.0637001991271973;
          break;
        case 14:
          r0ab = 4.1381998062133789;
          break;
        case 15:
          r0ab = 3.8591001033782959;
          break;
        case 16:
          r0ab = 4.0163998603820801;
          break;
        case 17:
          r0ab = 3.7878000736236572;
          break;
        case 18:
          r0ab = 3.6315999031066895;
          break;
        case 19:
          r0ab = 4.1740999221801758;
          break;
        case 20:
          r0ab = 4.3165998458862305;
          break;
        case 21:
          r0ab = 4.2395000457763672;
          break;
        case 22:
          r0ab = 4.1831002235412598;
          break;
        case 23:
          r0ab = 4.1107001304626465;
          break;
        case 24:
          r0ab = 3.5857000350952148;
          break;
        case 25:
          r0ab = 4.0269999504089355;
          break;
        case 26:
          r0ab = 3.9676001071929932;
          break;
        case 27:
          r0ab = 3.9463000297546387;
          break;
        case 28:
          r0ab = 3.9149999618530273;
          break;
        case 29:
          r0ab = 3.9021000862121582;
          break;
        case 30:
          r0ab = 3.8708000183105469;
          break;
        case 31:
          r0ab = 4.0240001678466797;
          break;
        case 32:
          r0ab = 4.1550998687744141;
          break;
        case 33:
          r0ab = 3.9107999801635742;
          break;
        case 34:
          r0ab = 4.1336998939514160;
          break;
        case 35:
          r0ab = 3.9289000034332275;
          break;
        case 36:
          r0ab = 3.7873001098632812;
          break;
        case 37:
          r0ab = 4.3666000366210938;
          break;
        case 38:
          r0ab = 4.5079998970031738;
          break;
        case 39:
          r0ab = 4.4232001304626465;
          break;
        case 40:
          r0ab = 4.3154997825622559;
          break;
        case 41:
          r0ab = 3.8461000919342041;
          break;
        case 42:
          r0ab = 3.8006999492645264;
          break;
        case 43:
          r0ab = 3.6991000175476074;
          break;
        case 44:
          r0ab = 3.6447000503540039;
          break;
        case 45:
          r0ab = 3.6308000087738037;
          break;
        case 46:
          r0ab = 3.5959000587463379;
          break;
        case 47:
          r0ab = 3.5748999118804932;
          break;
        case 48:
          r0ab = 4.0359001159667969;
          break;
        case 49:
          r0ab = 4.3123998641967773;
          break;
        case 50:
          r0ab = 4.3538999557495117;
          break;
        case 51:
          r0ab = 4.1121997833251953;
          break;
        case 52:
          r0ab = 4.3772001266479492;
          break;
        case 53:
          r0ab = 4.1785001754760742;
          break;
        case 54:
          r0ab = 4.0385999679565430;
          break;
        case 55:
          r0ab = 4.7003998756408691;
          break;
        case 56:
          r0ab = 4.8604001998901367;
          break;
        case 57:
          r0ab = 4.6261000633239746;
          break;
      }
      break;
    case 58:
      switch (atmA) {
        case 1:
          r0ab = 2.9454998970031738;
          break;
        case 2:
          r0ab = 3.2469999790191650;
          break;
        case 3:
          r0ab = 3.6108000278472900;
          break;
        case 4:
          r0ab = 3.8522000312805176;
          break;
        case 5:
          r0ab = 3.6624999046325684;
          break;
        case 6:
          r0ab = 3.6598000526428223;
          break;
        case 7:
          r0ab = 3.4410998821258545;
          break;
        case 8:
          r0ab = 3.4660000801086426;
          break;
        case 9:
          r0ab = 3.2414999008178711;
          break;
        case 10:
          r0ab = 3.0943999290466309;
          break;
        case 11:
          r0ab = 3.7513999938964844;
          break;
        case 12:
          r0ab = 4.0397000312805176;
          break;
        case 13:
          r0ab = 3.9230999946594238;
          break;
        case 14:
          r0ab = 4.0560998916625977;
          break;
        case 15:
          r0ab = 3.7860000133514404;
          break;
        case 16:
          r0ab = 3.9844999313354492;
          break;
        case 17:
          r0ab = 3.7453999519348145;
          break;
        case 18:
          r0ab = 3.5801999568939209;
          break;
        case 19:
          r0ab = 4.1366000175476074;
          break;
        case 20:
          r0ab = 4.3580999374389648;
          break;
        case 21:
          r0ab = 4.2350997924804688;
          break;
        case 22:
          r0ab = 4.2010998725891113;
          break;
        case 23:
          r0ab = 4.1402001380920410;
          break;
        case 24:
          r0ab = 3.5381000041961670;
          break;
        case 25:
          r0ab = 4.0652999877929688;
          break;
        case 26:
          r0ab = 4.0093002319335938;
          break;
        case 27:
          r0ab = 3.9883000850677490;
          break;
        case 28:
          r0ab = 3.9570000171661377;
          break;
        case 29:
          r0ab = 3.9428999423980713;
          break;
        case 30:
          r0ab = 3.9112000465393066;
          break;
        case 31:
          r0ab = 3.8728001117706299;
          break;
        case 32:
          r0ab = 4.0682001113891602;
          break;
        case 33:
          r0ab = 3.8350999355316162;
          break;
        case 34:
          r0ab = 4.1054000854492188;
          break;
        case 35:
          r0ab = 3.8928000926971436;
          break;
        case 36:
          r0ab = 3.7444999217987061;
          break;
        case 37:
          r0ab = 4.3414998054504395;
          break;
        case 38:
          r0ab = 4.5496997833251953;
          break;
        case 39:
          r0ab = 4.3832998275756836;
          break;
        case 40:
          r0ab = 4.3122000694274902;
          break;
        case 41:
          r0ab = 3.8050999641418457;
          break;
        case 42:
          r0ab = 3.7583000659942627;
          break;
        case 43:
          r0ab = 3.6621999740600586;
          break;
        case 44:
          r0ab = 3.6108000278472900;
          break;
        case 45:
          r0ab = 3.5971000194549561;
          break;
        case 46:
          r0ab = 3.5627999305725098;
          break;
        case 47:
          r0ab = 3.5408000946044922;
          break;
        case 48:
          r0ab = 4.0780000686645508;
          break;
        case 49:
          r0ab = 4.0727000236511230;
          break;
        case 50:
          r0ab = 4.2835998535156250;
          break;
        case 51:
          r0ab = 4.0553002357482910;
          break;
        case 52:
          r0ab = 4.3646998405456543;
          break;
        case 53:
          r0ab = 4.1621999740600586;
          break;
        case 54:
          r0ab = 4.0177998542785645;
          break;
        case 55:
          r0ab = 4.5802001953125000;
          break;
        case 56:
          r0ab = 4.9124999046325684;
          break;
        case 57:
          r0ab = 4.5861001014709473;
          break;
        case 58:
          r0ab = 4.6201000213623047;
          break;
      }
      break;
    case 59:
      switch (atmA) {
        case 1:
          r0ab = 2.9244000911712646;
          break;
        case 2:
          r0ab = 3.2241001129150391;
          break;
        case 3:
          r0ab = 3.5848000049591064;
          break;
        case 4:
          r0ab = 3.8292999267578125;
          break;
        case 5:
          r0ab = 3.6394999027252197;
          break;
        case 6:
          r0ab = 3.6400001049041748;
          break;
        case 7:
          r0ab = 3.4203999042510986;
          break;
        case 8:
          r0ab = 3.4498999118804932;
          break;
        case 9:
          r0ab = 3.2253000736236572;
          break;
        case 10:
          r0ab = 3.0778999328613281;
          break;
        case 11:
          r0ab = 3.7256999015808105;
          break;
        case 12:
          r0ab = 4.0170001983642578;
          break;
        case 13:
          r0ab = 3.9003000259399414;
          break;
        case 14:
          r0ab = 4.0371999740600586;
          break;
        case 15:
          r0ab = 3.7653000354766846;
          break;
        case 16:
          r0ab = 3.9672000408172607;
          break;
        case 17:
          r0ab = 3.7283000946044922;
          break;
        case 18:
          r0ab = 3.5629999637603760;
          break;
        case 19:
          r0ab = 4.1092000007629395;
          break;
        case 20:
          r0ab = 4.3347001075744629;
          break;
        case 21:
          r0ab = 4.2116999626159668;
          break;
        case 22:
          r0ab = 4.1792998313903809;
          break;
        case 23:
          r0ab = 4.1178998947143555;
          break;
        case 24:
          r0ab = 3.5139000415802002;
          break;
        case 25:
          r0ab = 4.0426001548767090;
          break;
        case 26:
          r0ab = 3.9867000579833984;
          break;
        case 27:
          r0ab = 3.9660999774932861;
          break;
        case 28:
          r0ab = 3.9344999790191650;
          break;
        case 29:
          r0ab = 3.9200000762939453;
          break;
        case 30:
          r0ab = 3.8882999420166016;
          break;
        case 31:
          r0ab = 3.8498001098632812;
          break;
        case 32:
          r0ab = 4.0496001243591309;
          break;
        case 33:
          r0ab = 3.8145000934600830;
          break;
        case 34:
          r0ab = 4.0880999565124512;
          break;
        case 35:
          r0ab = 3.8756000995635986;
          break;
        case 36:
          r0ab = 3.7270998954772949;
          break;
        case 37:
          r0ab = 4.3127999305725098;
          break;
        case 38:
          r0ab = 4.5241999626159668;
          break;
        case 39:
          r0ab = 4.3578000068664551;
          break;
        case 40:
          r0ab = 4.2870001792907715;
          break;
        case 41:
          r0ab = 3.7795999050140381;
          break;
        case 42:
          r0ab = 3.7318000793457031;
          break;
        case 43:
          r0ab = 3.6363999843597412;
          break;
        case 44:
          r0ab = 3.5854001045227051;
          break;
        case 45:
          r0ab = 3.5725998878479004;
          break;
        case 46:
          r0ab = 3.5378000736236572;
          break;
        case 47:
          r0ab = 3.5155000686645508;
          break;
        case 48:
          r0ab = 4.0527000427246094;
          break;
        case 49:
          r0ab = 4.0478000640869141;
          break;
        case 50:
          r0ab = 4.2630000114440918;
          break;
        case 51:
          r0ab = 4.0321998596191406;
          break;
        case 52:
          r0ab = 4.3449001312255859;
          break;
        case 53:
          r0ab = 4.1420998573303223;
          break;
        case 54:
          r0ab = 3.9974999427795410;
          break;
        case 55:
          r0ab = 4.5499000549316406;
          break;
        case 56:
          r0ab = 4.8825001716613770;
          break;
        case 57:
          r0ab = 4.5601000785827637;
          break;
        case 58:
          r0ab = 4.5949997901916504;
          break;
        case 59:
          r0ab = 4.5701999664306641;
          break;
      }
      break;
    case 60:
      switch (atmA) {
        case 1:
          r0ab = 2.9045999050140381;
          break;
        case 2:
          r0ab = 3.2044000625610352;
          break;
        case 3:
          r0ab = 3.5620999336242676;
          break;
        case 4:
          r0ab = 3.8078000545501709;
          break;
        case 5:
          r0ab = 3.6184999942779541;
          break;
        case 6:
          r0ab = 3.6219999790191650;
          break;
        case 7:
          r0ab = 3.4019000530242920;
          break;
        case 8:
          r0ab = 3.4358999729156494;
          break;
        case 9:
          r0ab = 3.2109999656677246;
          break;
        case 10:
          r0ab = 3.0634999275207520;
          break;
        case 11:
          r0ab = 3.7037000656127930;
          break;
        case 12:
          r0ab = 3.9958000183105469;
          break;
        case 13:
          r0ab = 3.8791999816894531;
          break;
        case 14:
          r0ab = 4.0194001197814941;
          break;
        case 15:
          r0ab = 3.7460000514984131;
          break;
        case 16:
          r0ab = 3.9516999721527100;
          break;
        case 17:
          r0ab = 3.7128000259399414;
          break;
        case 18:
          r0ab = 3.5473999977111816;
          break;
        case 19:
          r0ab = 4.0872001647949219;
          break;
        case 20:
          r0ab = 4.3137998580932617;
          break;
        case 21:
          r0ab = 4.1905999183654785;
          break;
        case 22:
          r0ab = 4.1592998504638672;
          break;
        case 23:
          r0ab = 4.0973000526428223;
          break;
        case 24:
          r0ab = 3.4918999671936035;
          break;
        case 25:
          r0ab = 4.0215997695922852;
          break;
        case 26:
          r0ab = 3.9656999111175537;
          break;
        case 27:
          r0ab = 3.9453999996185303;
          break;
        case 28:
          r0ab = 3.9133999347686768;
          break;
        case 29:
          r0ab = 3.8986001014709473;
          break;
        case 30:
          r0ab = 3.8668999671936035;
          break;
        case 31:
          r0ab = 3.8289000988006592;
          break;
        case 32:
          r0ab = 4.0322999954223633;
          break;
        case 33:
          r0ab = 3.7953999042510986;
          break;
        case 34:
          r0ab = 4.0725002288818359;
          break;
        case 35:
          r0ab = 3.8598001003265381;
          break;
        case 36:
          r0ab = 3.7112998962402344;
          break;
        case 37:
          r0ab = 4.2895998954772949;
          break;
        case 38:
          r0ab = 4.5020999908447266;
          break;
        case 39:
          r0ab = 4.3324999809265137;
          break;
        case 40:
          r0ab = 4.2645001411437988;
          break;
        case 41:
          r0ab = 3.7571001052856445;
          break;
        case 42:
          r0ab = 3.7083001136779785;
          break;
        case 43:
          r0ab = 3.6136000156402588;
          break;
        case 44:
          r0ab = 3.5627999305725098;
          break;
        case 45:
          r0ab = 3.5506999492645264;
          break;
        case 46:
          r0ab = 3.5155000686645508;
          break;
        case 47:
          r0ab = 3.4928998947143555;
          break;
        case 48:
          r0ab = 4.0296998023986816;
          break;
        case 49:
          r0ab = 4.0233998298645020;
          break;
        case 50:
          r0ab = 4.2442002296447754;
          break;
        case 51:
          r0ab = 4.0111999511718750;
          break;
        case 52:
          r0ab = 4.3274002075195312;
          break;
        case 53:
          r0ab = 4.1240000724792480;
          break;
        case 54:
          r0ab = 3.9793000221252441;
          break;
        case 55:
          r0ab = 4.5257000923156738;
          break;
        case 56:
          r0ab = 4.8568000793457031;
          break;
        case 57:
          r0ab = 4.5352997779846191;
          break;
        case 58:
          r0ab = 4.5732998847961426;
          break;
        case 59:
          r0ab = 4.5485000610351562;
          break;
        case 60:
          r0ab = 4.5271000862121582;
          break;
      }
      break;
    case 61:
      switch (atmA) {
        case 1:
          r0ab = 2.8877999782562256;
          break;
        case 2:
          r0ab = 3.1889998912811279;
          break;
        case 3:
          r0ab = 3.5411999225616455;
          break;
        case 4:
          r0ab = 3.7908000946044922;
          break;
        case 5:
          r0ab = 3.5973999500274658;
          break;
        case 6:
          r0ab = 3.6078000068664551;
          break;
        case 7:
          r0ab = 3.3870999813079834;
          break;
        case 8:
          r0ab = 3.4242999553680420;
          break;
        case 9:
          r0ab = 3.1991999149322510;
          break;
        case 10:
          r0ab = 3.0513000488281250;
          break;
        case 11:
          r0ab = 3.6830999851226807;
          break;
        case 12:
          r0ab = 3.9783999919891357;
          break;
        case 13:
          r0ab = 3.8578999042510986;
          break;
        case 14:
          r0ab = 4.0048999786376953;
          break;
        case 15:
          r0ab = 3.7304000854492188;
          break;
        case 16:
          r0ab = 3.9391999244689941;
          break;
        case 17:
          r0ab = 3.7002000808715820;
          break;
        case 18:
          r0ab = 3.5346999168395996;
          break;
        case 19:
          r0ab = 4.0657000541687012;
          break;
        case 20:
          r0ab = 4.2954998016357422;
          break;
        case 21:
          r0ab = 4.1704998016357422;
          break;
        case 22:
          r0ab = 4.1423997879028320;
          break;
        case 23:
          r0ab = 4.0799999237060547;
          break;
        case 24:
          r0ab = 3.4716999530792236;
          break;
        case 25:
          r0ab = 4.0043001174926758;
          break;
        case 26:
          r0ab = 3.9484999179840088;
          break;
        case 27:
          r0ab = 3.9286000728607178;
          break;
        case 28:
          r0ab = 3.8965001106262207;
          break;
        case 29:
          r0ab = 3.8815000057220459;
          break;
        case 30:
          r0ab = 3.8499999046325684;
          break;
        case 31:
          r0ab = 3.8073000907897949;
          break;
        case 32:
          r0ab = 4.0180001258850098;
          break;
        case 33:
          r0ab = 3.7795999050140381;
          break;
        case 34:
          r0ab = 4.0598001480102539;
          break;
        case 35:
          r0ab = 3.8469998836517334;
          break;
        case 36:
          r0ab = 3.6982998847961426;
          break;
        case 37:
          r0ab = 4.2677998542785645;
          break;
        case 38:
          r0ab = 4.4829998016357422;
          break;
        case 39:
          r0ab = 4.3131999969482422;
          break;
        case 40:
          r0ab = 4.2444000244140625;
          break;
        case 41:
          r0ab = 3.7369999885559082;
          break;
        case 42:
          r0ab = 3.6875998973846436;
          break;
        case 43:
          r0ab = 3.5934998989105225;
          break;
        case 44:
          r0ab = 3.5427999496459961;
          break;
        case 45:
          r0ab = 3.5313999652862549;
          break;
        case 46:
          r0ab = 3.4958000183105469;
          break;
        case 47:
          r0ab = 3.4730000495910645;
          break;
        case 48:
          r0ab = 4.0117001533508301;
          break;
        case 49:
          r0ab = 4.0043001174926758;
          break;
        case 50:
          r0ab = 4.2287001609802246;
          break;
        case 51:
          r0ab = 3.9939000606536865;
          break;
        case 52:
          r0ab = 4.3133997917175293;
          break;
        case 53:
          r0ab = 4.1096000671386719;
          break;
        case 54:
          r0ab = 3.9646000862121582;
          break;
        case 55:
          r0ab = 4.5032000541687012;
          break;
        case 56:
          r0ab = 4.8355998992919922;
          break;
        case 57:
          r0ab = 4.5156002044677734;
          break;
        case 58:
          r0ab = 4.5543999671936035;
          break;
        case 59:
          r0ab = 4.5296998023986816;
          break;
        case 60:
          r0ab = 4.5082998275756836;
          break;
        case 61:
          r0ab = 4.4896001815795898;
          break;
      }
      break;
    case 62:
      switch (atmA) {
        case 1:
          r0ab = 2.8708999156951904;
          break;
        case 2:
          r0ab = 3.1737000942230225;
          break;
        case 3:
          r0ab = 3.5199000835418701;
          break;
        case 4:
          r0ab = 3.7734000682830811;
          break;
        case 5:
          r0ab = 3.5801999568939209;
          break;
        case 6:
          r0ab = 3.5934000015258789;
          break;
        case 7:
          r0ab = 3.3724000453948975;
          break;
        case 8:
          r0ab = 3.4128000736236572;
          break;
        case 9:
          r0ab = 3.1877000331878662;
          break;
        case 10:
          r0ab = 3.0395998954772949;
          break;
        case 11:
          r0ab = 3.6624000072479248;
          break;
        case 12:
          r0ab = 3.9607999324798584;
          break;
        case 13:
          r0ab = 3.8396999835968018;
          break;
        case 14:
          r0ab = 3.9893000125885010;
          break;
        case 15:
          r0ab = 3.7144999504089355;
          break;
        case 16:
          r0ab = 3.9265999794006348;
          break;
        case 17:
          r0ab = 3.6877000331878662;
          break;
        case 18:
          r0ab = 3.5222001075744629;
          break;
        case 19:
          r0ab = 4.0447998046875000;
          break;
        case 20:
          r0ab = 4.2771000862121582;
          break;
        case 21:
          r0ab = 4.1522998809814453;
          break;
        case 22:
          r0ab = 4.1247000694274902;
          break;
        case 23:
          r0ab = 4.0626001358032227;
          break;
        case 24:
          r0ab = 3.4530000686645508;
          break;
        case 25:
          r0ab = 3.9865999221801758;
          break;
        case 26:
          r0ab = 3.9309999942779541;
          break;
        case 27:
          r0ab = 3.9114999771118164;
          break;
        case 28:
          r0ab = 3.8791999816894531;
          break;
        case 29:
          r0ab = 3.8640999794006348;
          break;
        case 30:
          r0ab = 3.8326001167297363;
          break;
        case 31:
          r0ab = 3.7892000675201416;
          break;
        case 32:
          r0ab = 4.0025000572204590;
          break;
        case 33:
          r0ab = 3.7636001110076904;
          break;
        case 34:
          r0ab = 4.0471000671386719;
          break;
        case 35:
          r0ab = 3.8343000411987305;
          break;
        case 36:
          r0ab = 3.6854000091552734;
          break;
        case 37:
          r0ab = 4.2463998794555664;
          break;
        case 38:
          r0ab = 4.4635000228881836;
          break;
        case 39:
          r0ab = 4.2939000129699707;
          break;
        case 40:
          r0ab = 4.2252001762390137;
          break;
        case 41:
          r0ab = 3.7169001102447510;
          break;
        case 42:
          r0ab = 3.6675000190734863;
          break;
        case 43:
          r0ab = 3.5738999843597412;
          break;
        case 44:
          r0ab = 3.5234999656677246;
          break;
        case 45:
          r0ab = 3.5125999450683594;
          break;
        case 46:
          r0ab = 3.4767999649047852;
          break;
        case 47:
          r0ab = 3.4537000656127930;
          break;
        case 48:
          r0ab = 3.9932000637054443;
          break;
        case 49:
          r0ab = 3.9853999614715576;
          break;
        case 50:
          r0ab = 4.2122998237609863;
          break;
        case 51:
          r0ab = 3.9765000343322754;
          break;
        case 52:
          r0ab = 4.2992000579833984;
          break;
        case 53:
          r0ab = 4.0950999259948730;
          break;
        case 54:
          r0ab = 3.9500000476837158;
          break;
        case 55:
          r0ab = 4.4811000823974609;
          break;
        case 56:
          r0ab = 4.8134999275207520;
          break;
        case 57:
          r0ab = 4.4959001541137695;
          break;
        case 58:
          r0ab = 4.5350999832153320;
          break;
        case 59:
          r0ab = 4.5104999542236328;
          break;
        case 60:
          r0ab = 4.4890999794006348;
          break;
        case 61:
          r0ab = 4.4704999923706055;
          break;
        case 62:
          r0ab = 4.4514999389648438;
          break;
      }
      break;
    case 63:
      switch (atmA) {
        case 1:
          r0ab = 2.8568000793457031;
          break;
        case 2:
          r0ab = 3.1607999801635742;
          break;
        case 3:
          r0ab = 3.5050001144409180;
          break;
        case 4:
          r0ab = 3.7597999572753906;
          break;
        case 5:
          r0ab = 3.5664999485015869;
          break;
        case 6:
          r0ab = 3.5803000926971436;
          break;
        case 7:
          r0ab = 3.3601000308990479;
          break;
        case 8:
          r0ab = 3.4031000137329102;
          break;
        case 9:
          r0ab = 3.1779000759124756;
          break;
        case 10:
          r0ab = 3.0295999050140381;
          break;
        case 11:
          r0ab = 3.6479001045227051;
          break;
        case 12:
          r0ab = 3.9470999240875244;
          break;
        case 13:
          r0ab = 3.8262000083923340;
          break;
        case 14:
          r0ab = 3.9772999286651611;
          break;
        case 15:
          r0ab = 3.7014999389648438;
          break;
        case 16:
          r0ab = 3.9161999225616455;
          break;
        case 17:
          r0ab = 3.6770999431610107;
          break;
        case 18:
          r0ab = 3.5114998817443848;
          break;
        case 19:
          r0ab = 4.0306000709533691;
          break;
        case 20:
          r0ab = 4.2634000778198242;
          break;
        case 21:
          r0ab = 4.1385002136230469;
          break;
        case 22:
          r0ab = 4.1115999221801758;
          break;
        case 23:
          r0ab = 4.0489001274108887;
          break;
        case 24:
          r0ab = 3.4365999698638916;
          break;
        case 25:
          r0ab = 3.9732000827789307;
          break;
        case 26:
          r0ab = 3.9175999164581299;
          break;
        case 27:
          r0ab = 3.8982999324798584;
          break;
        case 28:
          r0ab = 3.8659000396728516;
          break;
        case 29:
          r0ab = 3.8506999015808105;
          break;
        case 30:
          r0ab = 3.8190999031066895;
          break;
        case 31:
          r0ab = 3.7757000923156738;
          break;
        case 32:
          r0ab = 3.9907000064849854;
          break;
        case 33:
          r0ab = 3.7506000995635986;
          break;
        case 34:
          r0ab = 4.0364999771118164;
          break;
        case 35:
          r0ab = 3.8234999179840088;
          break;
        case 36:
          r0ab = 3.6744999885559082;
          break;
        case 37:
          r0ab = 4.2314000129699707;
          break;
        case 38:
          r0ab = 4.4489998817443848;
          break;
        case 39:
          r0ab = 4.2792000770568848;
          break;
        case 40:
          r0ab = 4.2104997634887695;
          break;
        case 41:
          r0ab = 3.7002999782562256;
          break;
        case 42:
          r0ab = 3.6510000228881836;
          break;
        case 43:
          r0ab = 3.5578000545501709;
          break;
        case 44:
          r0ab = 3.5074999332427979;
          break;
        case 45:
          r0ab = 3.4971001148223877;
          break;
        case 46:
          r0ab = 3.4609000682830811;
          break;
        case 47:
          r0ab = 3.4377000331878662;
          break;
        case 48:
          r0ab = 3.9788000583648682;
          break;
        case 49:
          r0ab = 3.9711999893188477;
          break;
        case 50:
          r0ab = 4.1996998786926270;
          break;
        case 51:
          r0ab = 3.9623999595642090;
          break;
        case 52:
          r0ab = 4.2877001762390137;
          break;
        case 53:
          r0ab = 4.0830998420715332;
          break;
        case 54:
          r0ab = 3.9377999305725098;
          break;
        case 55:
          r0ab = 4.4654998779296875;
          break;
        case 56:
          r0ab = 4.7973999977111816;
          break;
        case 57:
          r0ab = 4.4812998771667480;
          break;
        case 58:
          r0ab = 4.5208997726440430;
          break;
        case 59:
          r0ab = 4.4963998794555664;
          break;
        case 60:
          r0ab = 4.4749999046325684;
          break;
        case 61:
          r0ab = 4.4565000534057617;
          break;
        case 62:
          r0ab = 4.4375000000000000;
          break;
        case 63:
          r0ab = 4.4233999252319336;
          break;
      }
      break;
    case 64:
      switch (atmA) {
        case 1:
          r0ab = 2.6798000335693359;
          break;
        case 2:
          r0ab = 3.0151000022888184;
          break;
        case 3:
          r0ab = 3.2585999965667725;
          break;
        case 4:
          r0ab = 3.5292000770568848;
          break;
        case 5:
          r0ab = 3.5390999317169189;
          break;
        case 6:
          r0ab = 3.4902000427246094;
          break;
        case 7:
          r0ab = 3.2887001037597656;
          break;
        case 8:
          r0ab = 3.3322000503540039;
          break;
        case 9:
          r0ab = 3.1228001117706299;
          break;
        case 10:
          r0ab = 2.9888000488281250;
          break;
        case 11:
          r0ab = 3.4012000560760498;
          break;
        case 12:
          r0ab = 3.7144999504089355;
          break;
        case 13:
          r0ab = 3.7829999923706055;
          break;
        case 14:
          r0ab = 3.6665000915527344;
          break;
        case 15:
          r0ab = 3.5897998809814453;
          break;
        case 16:
          r0ab = 3.8076999187469482;
          break;
        case 17:
          r0ab = 3.5810000896453857;
          break;
        case 18:
          r0ab = 3.4265000820159912;
          break;
        case 19:
          r0ab = 3.7725999355316162;
          break;
        case 20:
          r0ab = 4.0307002067565918;
          break;
        case 21:
          r0ab = 3.9763000011444092;
          break;
        case 22:
          r0ab = 3.8889999389648438;
          break;
        case 23:
          r0ab = 3.8489000797271729;
          break;
        case 24:
          r0ab = 3.2706000804901123;
          break;
        case 25:
          r0ab = 3.7595000267028809;
          break;
        case 26:
          r0ab = 3.6984000205993652;
          break;
        case 27:
          r0ab = 3.6772000789642334;
          break;
        case 28:
          r0ab = 3.6428000926971436;
          break;
        case 29:
          r0ab = 3.6243000030517578;
          break;
        case 30:
          r0ab = 3.5950999259948730;
          break;
        case 31:
          r0ab = 3.7497000694274902;
          break;
        case 32:
          r0ab = 3.6775000095367432;
          break;
        case 33:
          r0ab = 3.6363999843597412;
          break;
        case 34:
          r0ab = 3.9203000068664551;
          break;
        case 35:
          r0ab = 3.7156999111175537;
          break;
        case 36:
          r0ab = 3.5745999813079834;
          break;
        case 37:
          r0ab = 3.9493999481201172;
          break;
        case 38:
          r0ab = 4.2076001167297363;
          break;
        case 39:
          r0ab = 4.1563000679016113;
          break;
        case 40:
          r0ab = 4.0507998466491699;
          break;
        case 41:
          r0ab = 3.5329000949859619;
          break;
        case 42:
          r0ab = 3.4779999256134033;
          break;
        case 43:
          r0ab = 3.3731000423431396;
          break;
        case 44:
          r0ab = 3.3125998973846436;
          break;
        case 45:
          r0ab = 3.2846000194549561;
          break;
        case 46:
          r0ab = 3.2425999641418457;
          break;
        case 47:
          r0ab = 3.2135000228881836;
          break;
        case 48:
          r0ab = 3.7490999698638916;
          break;
        case 49:
          r0ab = 3.9005999565124512;
          break;
        case 50:
          r0ab = 3.8331999778747559;
          break;
        case 51:
          r0ab = 3.8029000759124756;
          break;
        case 52:
          r0ab = 4.1435999870300293;
          break;
        case 53:
          r0ab = 3.9407000541687012;
          break;
        case 54:
          r0ab = 3.7997999191284180;
          break;
        case 55:
          r0ab = 4.1662998199462891;
          break;
        case 56:
          r0ab = 4.5309000015258789;
          break;
        case 57:
          r0ab = 4.3481001853942871;
          break;
        case 58:
          r0ab = 4.2911000251770020;
          break;
        case 59:
          r0ab = 4.2670998573303223;
          break;
        case 60:
          r0ab = 4.2414999008178711;
          break;
        case 61:
          r0ab = 4.2230000495910645;
          break;
        case 62:
          r0ab = 4.2046999931335449;
          break;
        case 63:
          r0ab = 4.1908001899719238;
          break;
        case 64:
          r0ab = 4.1243000030517578;
          break;
      }
      break;
    case 65:
      switch (atmA) {
        case 1:
          r0ab = 2.5188999176025391;
          break;
        case 2:
          r0ab = 2.9702999591827393;
          break;
        case 3:
          r0ab = 3.3062999248504639;
          break;
        case 4:
          r0ab = 3.6235001087188721;
          break;
        case 5:
          r0ab = 3.4516999721527100;
          break;
        case 6:
          r0ab = 3.3989000320434570;
          break;
        case 7:
          r0ab = 3.2107000350952148;
          break;
        case 8:
          r0ab = 3.2434000968933105;
          break;
        case 9:
          r0ab = 3.0093998908996582;
          break;
        case 10:
          r0ab = 2.8580000400543213;
          break;
        case 11:
          r0ab = 3.4252998828887939;
          break;
        case 12:
          r0ab = 3.8157000541687012;
          break;
        case 13:
          r0ab = 3.7258000373840332;
          break;
        case 14:
          r0ab = 3.6131999492645264;
          break;
        case 15:
          r0ab = 3.5297000408172607;
          break;
        case 16:
          r0ab = 3.7565999031066895;
          break;
        case 17:
          r0ab = 3.5095000267028809;
          break;
        case 18:
          r0ab = 3.3368000984191895;
          break;
        case 19:
          r0ab = 3.7890000343322754;
          break;
        case 20:
          r0ab = 4.1297998428344727;
          break;
        case 21:
          r0ab = 4.0190000534057617;
          break;
        case 22:
          r0ab = 3.9572999477386475;
          break;
        case 23:
          r0ab = 3.9237000942230225;
          break;
        case 24:
          r0ab = 3.2676999568939209;
          break;
        case 25:
          r0ab = 3.8480000495910645;
          break;
        case 26:
          r0ab = 3.8157000541687012;
          break;
        case 27:
          r0ab = 3.7655999660491943;
          break;
        case 28:
          r0ab = 3.7316999435424805;
          break;
        case 29:
          r0ab = 3.7125999927520752;
          break;
        case 30:
          r0ab = 3.6814000606536865;
          break;
        case 31:
          r0ab = 3.6793000698089600;
          break;
        case 32:
          r0ab = 3.6217999458312988;
          break;
        case 33:
          r0ab = 3.5787999629974365;
          break;
        case 34:
          r0ab = 3.8763000965118408;
          break;
        case 35:
          r0ab = 3.6572000980377197;
          break;
        case 36:
          r0ab = 3.5021998882293701;
          break;
        case 37:
          r0ab = 3.9737000465393066;
          break;
        case 38:
          r0ab = 4.3255000114440918;
          break;
        case 39:
          r0ab = 4.1827998161315918;
          break;
        case 40:
          r0ab = 4.1157999038696289;
          break;
        case 41:
          r0ab = 3.5078001022338867;
          break;
        case 42:
          r0ab = 3.4595000743865967;
          break;
        case 43:
          r0ab = 3.3599998950958252;
          break;
        case 44:
          r0ab = 3.3087999820709229;
          break;
        case 45:
          r0ab = 3.2574999332427979;
          break;
        case 46:
          r0ab = 3.2163999080657959;
          break;
        case 47:
          r0ab = 3.1856000423431396;
          break;
        case 48:
          r0ab = 3.8522000312805176;
          break;
        case 49:
          r0ab = 3.8664999008178711;
          break;
        case 50:
          r0ab = 3.8074998855590820;
          break;
        case 51:
          r0ab = 3.7771999835968018;
          break;
        case 52:
          r0ab = 4.1391000747680664;
          break;
        case 53:
          r0ab = 3.9296000003814697;
          break;
        case 54:
          r0ab = 3.7771999835968018;
          break;
        case 55:
          r0ab = 4.2133998870849609;
          break;
        case 56:
          r0ab = 4.7308001518249512;
          break;
        case 57:
          r0ab = 4.3786997795104980;
          break;
        case 58:
          r0ab = 4.3894000053405762;
          break;
        case 59:
          r0ab = 4.3649001121520996;
          break;
        case 60:
          r0ab = 4.3440999984741211;
          break;
        case 61:
          r0ab = 4.3256998062133789;
          break;
        case 62:
          r0ab = 4.3073000907897949;
          break;
        case 63:
          r0ab = 4.2940998077392578;
          break;
        case 64:
          r0ab = 4.1251997947692871;
          break;
        case 65:
          r0ab = 4.2427000999450684;
          break;
      }
      break;
    case 66:
      switch (atmA) {
        case 1:
          r0ab = 3.0480999946594238;
          break;
        case 2:
          r0ab = 2.9584000110626221;
          break;
        case 3:
          r0ab = 3.6919000148773193;
          break;
        case 4:
          r0ab = 3.5989999771118164;
          break;
        case 5:
          r0ab = 3.8880999088287354;
          break;
        case 6:
          r0ab = 3.4209001064300537;
          break;
        case 7:
          r0ab = 3.1605999469757080;
          break;
        case 8:
          r0ab = 3.1937999725341797;
          break;
        case 9:
          r0ab = 2.9974999427795410;
          break;
        case 10:
          r0ab = 2.8645999431610107;
          break;
        case 11:
          r0ab = 3.8138000965118408;
          break;
        case 12:
          r0ab = 3.7934999465942383;
          break;
        case 13:
          r0ab = 3.7081000804901123;
          break;
        case 14:
          r0ab = 3.9154999256134033;
          break;
        case 15:
          r0ab = 3.5910000801086426;
          break;
        case 16:
          r0ab = 3.4807999134063721;
          break;
        case 17:
          r0ab = 3.4886000156402588;
          break;
        case 18:
          r0ab = 3.3396999835968018;
          break;
        case 19:
          r0ab = 4.1336002349853516;
          break;
        case 20:
          r0ab = 4.1121997833251953;
          break;
        case 21:
          r0ab = 3.9888000488281250;
          break;
        case 22:
          r0ab = 3.9542999267578125;
          break;
        case 23:
          r0ab = 3.8917000293731689;
          break;
        case 24:
          r0ab = 3.5894000530242920;
          break;
        case 25:
          r0ab = 3.8131000995635986;
          break;
        case 26:
          r0ab = 3.7634999752044678;
          break;
        case 27:
          r0ab = 3.7418999671936035;
          break;
        case 28:
          r0ab = 3.7070999145507812;
          break;
        case 29:
          r0ab = 3.6879999637603760;
          break;
        case 30:
          r0ab = 3.6573998928070068;
          break;
        case 31:
          r0ab = 3.6545999050140381;
          break;
        case 32:
          r0ab = 3.9375000000000000;
          break;
        case 33:
          r0ab = 3.6579000949859619;
          break;
        case 34:
          r0ab = 3.5869998931884766;
          break;
        case 35:
          r0ab = 3.6361000537872314;
          break;
        case 36:
          r0ab = 3.5039000511169434;
          break;
        case 37:
          r0ab = 4.3148999214172363;
          break;
        case 38:
          r0ab = 4.2978000640869141;
          break;
        case 39:
          r0ab = 4.1321001052856445;
          break;
        case 40:
          r0ab = 4.1297998428344727;
          break;
        case 41:
          r0ab = 3.8164000511169434;
          break;
        case 42:
          r0ab = 3.7679998874664307;
          break;
        case 43:
          r0ab = 3.7153999805450439;
          break;
        case 44:
          r0ab = 3.6858000755310059;
          break;
        case 45:
          r0ab = 3.6709001064300537;
          break;
        case 46:
          r0ab = 3.6665999889373779;
          break;
        case 47:
          r0ab = 3.6517000198364258;
          break;
        case 48:
          r0ab = 3.8173999786376953;
          break;
        case 49:
          r0ab = 3.8608000278472900;
          break;
        case 50:
          r0ab = 4.1805000305175781;
          break;
        case 51:
          r0ab = 3.9102001190185547;
          break;
        case 52:
          r0ab = 3.8394000530242920;
          break;
        case 53:
          r0ab = 3.8968000411987305;
          break;
        case 54:
          r0ab = 3.7672998905181885;
          break;
        case 55:
          r0ab = 4.5274000167846680;
          break;
        case 56:
          r0ab = 4.6682000160217285;
          break;
        case 57:
          r0ab = 4.3344001770019531;
          break;
        case 58:
          r0ab = 4.3639001846313477;
          break;
        case 59:
          r0ab = 4.3383998870849609;
          break;
        case 60:
          r0ab = 4.3161997795104980;
          break;
        case 61:
          r0ab = 4.2972002029418945;
          break;
        case 62:
          r0ab = 4.2779002189636230;
          break;
        case 63:
          r0ab = 4.2635998725891113;
          break;
        case 64:
          r0ab = 4.0253000259399414;
          break;
        case 65:
          r0ab = 4.1167998313903809;
          break;
        case 66:
          r0ab = 4.1540999412536621;
          break;
      }
      break;
    case 67:
      switch (atmA) {
        case 1:
          r0ab = 2.8136000633239746;
          break;
        case 2:
          r0ab = 3.0950999259948730;
          break;
        case 3:
          r0ab = 3.4635000228881836;
          break;
        case 4:
          r0ab = 3.6875000000000000;
          break;
        case 5:
          r0ab = 3.4986999034881592;
          break;
        case 6:
          r0ab = 3.5183000564575195;
          break;
        case 7:
          r0ab = 3.2936999797821045;
          break;
        case 8:
          r0ab = 3.3580000400543213;
          break;
        case 9:
          r0ab = 3.1324999332427979;
          break;
        case 10:
          r0ab = 2.9832000732421875;
          break;
        case 11:
          r0ab = 3.6078000068664551;
          break;
        case 12:
          r0ab = 3.8756999969482422;
          break;
        case 13:
          r0ab = 3.7616000175476074;
          break;
        case 14:
          r0ab = 3.9221999645233154;
          break;
        case 15:
          r0ab = 3.6370000839233398;
          break;
        case 16:
          r0ab = 3.8647000789642334;
          break;
        case 17:
          r0ab = 3.6256000995635986;
          break;
        case 18:
          r0ab = 3.4595000743865967;
          break;
        case 19:
          r0ab = 3.9874000549316406;
          break;
        case 20:
          r0ab = 4.1937999725341797;
          break;
        case 21:
          r0ab = 4.0679001808166504;
          break;
        case 22:
          r0ab = 4.0430002212524414;
          break;
        case 23:
          r0ab = 3.9781000614166260;
          break;
        case 24:
          r0ab = 3.3886001110076904;
          break;
        case 25:
          r0ab = 3.9007999897003174;
          break;
        case 26:
          r0ab = 3.8462998867034912;
          break;
        case 27:
          r0ab = 3.8287999629974365;
          break;
        case 28:
          r0ab = 3.7950000762939453;
          break;
        case 29:
          r0ab = 3.7790000438690186;
          break;
        case 30:
          r0ab = 3.7472000122070312;
          break;
        case 31:
          r0ab = 3.7116999626159668;
          break;
        case 32:
          r0ab = 3.9370999336242676;
          break;
        case 33:
          r0ab = 3.6872999668121338;
          break;
        case 34:
          r0ab = 3.9846000671386719;
          break;
        case 35:
          r0ab = 3.7709000110626221;
          break;
        case 36:
          r0ab = 3.6210000514984131;
          break;
        case 37:
          r0ab = 4.1812000274658203;
          break;
        case 38:
          r0ab = 4.3750000000000000;
          break;
        case 39:
          r0ab = 4.2044000625610352;
          break;
        case 40:
          r0ab = 4.1339998245239258;
          break;
        case 41:
          r0ab = 3.6459000110626221;
          break;
        case 42:
          r0ab = 3.5929000377655029;
          break;
        case 43:
          r0ab = 3.5035998821258545;
          break;
        case 44:
          r0ab = 3.4577000141143799;
          break;
        case 45:
          r0ab = 3.4528000354766846;
          break;
        case 46:
          r0ab = 3.4145998954772949;
          break;
        case 47:
          r0ab = 3.3903999328613281;
          break;
        case 48:
          r0ab = 3.9014000892639160;
          break;
        case 49:
          r0ab = 3.9031000137329102;
          break;
        case 50:
          r0ab = 4.1442999839782715;
          break;
        case 51:
          r0ab = 3.8961000442504883;
          break;
        case 52:
          r0ab = 4.2294998168945312;
          break;
        case 53:
          r0ab = 4.0226998329162598;
          break;
        case 54:
          r0ab = 3.8763000965118408;
          break;
        case 55:
          r0ab = 4.4085998535156250;
          break;
        case 56:
          r0ab = 4.7097001075744629;
          break;
        case 57:
          r0ab = 4.4064002037048340;
          break;
        case 58:
          r0ab = 4.4488000869750977;
          break;
        case 59:
          r0ab = 4.4243001937866211;
          break;
        case 60:
          r0ab = 4.4029002189636230;
          break;
        case 61:
          r0ab = 4.3842000961303711;
          break;
        case 62:
          r0ab = 4.3654999732971191;
          break;
        case 63:
          r0ab = 4.3513998985290527;
          break;
        case 64:
          r0ab = 4.1161999702453613;
          break;
        case 65:
          r0ab = 4.2204999923706055;
          break;
        case 66:
          r0ab = 4.1953001022338867;
          break;
        case 67:
          r0ab = 4.2793998718261719;
          break;
      }
      break;
    case 68:
      switch (atmA) {
        case 1:
          r0ab = 2.8032000064849854;
          break;
        case 2:
          r0ab = 3.0804998874664307;
          break;
        case 3:
          r0ab = 3.4519000053405762;
          break;
        case 4:
          r0ab = 3.6700000762939453;
          break;
        case 5:
          r0ab = 3.4827001094818115;
          break;
        case 6:
          r0ab = 3.5050001144409180;
          break;
        case 7:
          r0ab = 3.2799000740051270;
          break;
        case 8:
          r0ab = 3.3482000827789307;
          break;
        case 9:
          r0ab = 3.1233000755310059;
          break;
        case 10:
          r0ab = 2.9746999740600586;
          break;
        case 11:
          r0ab = 3.5971000194549561;
          break;
        case 12:
          r0ab = 3.8585999011993408;
          break;
        case 13:
          r0ab = 3.7460999488830566;
          break;
        case 14:
          r0ab = 3.9100000858306885;
          break;
        case 15:
          r0ab = 3.6228001117706299;
          break;
        case 16:
          r0ab = 3.8534998893737793;
          break;
        case 17:
          r0ab = 3.6147000789642334;
          break;
        case 18:
          r0ab = 3.4489998817443848;
          break;
        case 19:
          r0ab = 3.9763998985290527;
          break;
        case 20:
          r0ab = 4.1772999763488770;
          break;
        case 21:
          r0ab = 4.0510997772216797;
          break;
        case 22:
          r0ab = 4.0269999504089355;
          break;
        case 23:
          r0ab = 3.9614000320434570;
          break;
        case 24:
          r0ab = 3.3754000663757324;
          break;
        case 25:
          r0ab = 3.8835999965667725;
          break;
        case 26:
          r0ab = 3.8290998935699463;
          break;
        case 27:
          r0ab = 3.8120999336242676;
          break;
        case 28:
          r0ab = 3.7780001163482666;
          break;
        case 29:
          r0ab = 3.7618999481201172;
          break;
        case 30:
          r0ab = 3.7300000190734863;
          break;
        case 31:
          r0ab = 3.6965000629425049;
          break;
        case 32:
          r0ab = 3.9252998828887939;
          break;
        case 33:
          r0ab = 3.6733999252319336;
          break;
        case 34:
          r0ab = 3.9732999801635742;
          break;
        case 35:
          r0ab = 3.7597000598907471;
          break;
        case 36:
          r0ab = 3.6098999977111816;
          break;
        case 37:
          r0ab = 4.1683001518249512;
          break;
        case 38:
          r0ab = 4.3572001457214355;
          break;
        case 39:
          r0ab = 4.1862001419067383;
          break;
        case 40:
          r0ab = 4.1153001785278320;
          break;
        case 41:
          r0ab = 3.6312000751495361;
          break;
        case 42:
          r0ab = 3.5771999359130859;
          break;
        case 43:
          r0ab = 3.4881000518798828;
          break;
        case 44:
          r0ab = 3.4428999423980713;
          break;
        case 45:
          r0ab = 3.4395000934600830;
          break;
        case 46:
          r0ab = 3.4008998870849609;
          break;
        case 47:
          r0ab = 3.3766000270843506;
          break;
        case 48:
          r0ab = 3.8826999664306641;
          break;
        case 49:
          r0ab = 3.8868000507354736;
          break;
        case 50:
          r0ab = 4.1315999031066895;
          break;
        case 51:
          r0ab = 3.8807001113891602;
          break;
        case 52:
          r0ab = 4.2164001464843750;
          break;
        case 53:
          r0ab = 4.0092000961303711;
          break;
        case 54:
          r0ab = 3.8626999855041504;
          break;
        case 55:
          r0ab = 4.3935999870300293;
          break;
        case 56:
          r0ab = 4.6870999336242676;
          break;
        case 57:
          r0ab = 4.3881998062133789;
          break;
        case 58:
          r0ab = 4.4316000938415527;
          break;
        case 59:
          r0ab = 4.4072999954223633;
          break;
        case 60:
          r0ab = 4.3857998847961426;
          break;
        case 61:
          r0ab = 4.3671998977661133;
          break;
        case 62:
          r0ab = 4.3484997749328613;
          break;
        case 63:
          r0ab = 4.3344001770019531;
          break;
        case 64:
          r0ab = 4.0984001159667969;
          break;
        case 65:
          r0ab = 4.2035999298095703;
          break;
        case 66:
          r0ab = 4.1791000366210938;
          break;
        case 67:
          r0ab = 4.2621998786926270;
          break;
        case 68:
          r0ab = 4.2449998855590820;
          break;
      }
      break;
    case 69:
      switch (atmA) {
        case 1:
          r0ab = 2.7967000007629395;
          break;
        case 2:
          r0ab = 3.0689001083374023;
          break;
        case 3:
          r0ab = 3.4444999694824219;
          break;
        case 4:
          r0ab = 3.6580998897552490;
          break;
        case 5:
          r0ab = 3.4716999530792236;
          break;
        case 6:
          r0ab = 3.4951000213623047;
          break;
        case 7:
          r0ab = 3.2693998813629150;
          break;
        case 8:
          r0ab = 3.3396999835968018;
          break;
        case 9:
          r0ab = 3.1147000789642334;
          break;
        case 10:
          r0ab = 2.9660999774932861;
          break;
        case 11:
          r0ab = 3.5897998809814453;
          break;
        case 12:
          r0ab = 3.8468000888824463;
          break;
        case 13:
          r0ab = 3.7358000278472900;
          break;
        case 14:
          r0ab = 3.9014000892639160;
          break;
        case 15:
          r0ab = 3.6129000186920166;
          break;
        case 16:
          r0ab = 3.8443000316619873;
          break;
        case 17:
          r0ab = 3.6054000854492188;
          break;
        case 18:
          r0ab = 3.4395999908447266;
          break;
        case 19:
          r0ab = 3.9683001041412354;
          break;
        case 20:
          r0ab = 4.1655998229980469;
          break;
        case 21:
          r0ab = 4.0394001007080078;
          break;
        case 22:
          r0ab = 4.0157999992370605;
          break;
        case 23:
          r0ab = 3.9498000144958496;
          break;
        case 24:
          r0ab = 3.3677000999450684;
          break;
        case 25:
          r0ab = 3.8717999458312988;
          break;
        case 26:
          r0ab = 3.8164000511169434;
          break;
        case 27:
          r0ab = 3.8004999160766602;
          break;
        case 28:
          r0ab = 3.7662000656127930;
          break;
        case 29:
          r0ab = 3.7500000000000000;
          break;
        case 30:
          r0ab = 3.7181000709533691;
          break;
        case 31:
          r0ab = 3.6863000392913818;
          break;
        case 32:
          r0ab = 3.9170000553131104;
          break;
        case 33:
          r0ab = 3.6637001037597656;
          break;
        case 34:
          r0ab = 3.9640998840332031;
          break;
        case 35:
          r0ab = 3.7502999305725098;
          break;
        case 36:
          r0ab = 3.6003999710083008;
          break;
        case 37:
          r0ab = 4.1589999198913574;
          break;
        case 38:
          r0ab = 4.3447999954223633;
          break;
        case 39:
          r0ab = 4.1739001274108887;
          break;
        case 40:
          r0ab = 4.1029000282287598;
          break;
        case 41:
          r0ab = 3.6224000453948975;
          break;
        case 42:
          r0ab = 3.5676999092102051;
          break;
        case 43:
          r0ab = 3.4784998893737793;
          break;
        case 44:
          r0ab = 3.4314000606536865;
          break;
        case 45:
          r0ab = 3.4312999248504639;
          break;
        case 46:
          r0ab = 3.3922998905181885;
          break;
        case 47:
          r0ab = 3.3680000305175781;
          break;
        case 48:
          r0ab = 3.8698000907897949;
          break;
        case 49:
          r0ab = 3.8757998943328857;
          break;
        case 50:
          r0ab = 4.1229000091552734;
          break;
        case 51:
          r0ab = 3.8703999519348145;
          break;
        case 52:
          r0ab = 4.2062997817993164;
          break;
        case 53:
          r0ab = 3.9986999034881592;
          break;
        case 54:
          r0ab = 3.8519001007080078;
          break;
        case 55:
          r0ab = 4.3832001686096191;
          break;
        case 56:
          r0ab = 4.6728000640869141;
          break;
        case 57:
          r0ab = 4.3758997917175293;
          break;
        case 58:
          r0ab = 4.4194998741149902;
          break;
        case 59:
          r0ab = 4.3951997756958008;
          break;
        case 60:
          r0ab = 4.3737001419067383;
          break;
        case 61:
          r0ab = 4.3551001548767090;
          break;
        case 62:
          r0ab = 4.3364000320434570;
          break;
        case 63:
          r0ab = 4.3222999572753906;
          break;
        case 64:
          r0ab = 4.0861001014709473;
          break;
        case 65:
          r0ab = 4.1911001205444336;
          break;
        case 66:
          r0ab = 4.1676001548767090;
          break;
        case 67:
          r0ab = 4.2501001358032227;
          break;
        case 68:
          r0ab = 4.2329001426696777;
          break;
        case 69:
          r0ab = 4.2207999229431152;
          break;
      }
      break;
    case 70:
      switch (atmA) {
        case 1:
          r0ab = 2.7897000312805176;
          break;
        case 2:
          r0ab = 3.0636000633239746;
          break;
        case 3:
          r0ab = 3.4344000816345215;
          break;
        case 4:
          r0ab = 3.6480000019073486;
          break;
        case 5:
          r0ab = 3.4625999927520752;
          break;
        case 6:
          r0ab = 3.4892001152038574;
          break;
        case 7:
          r0ab = 3.2625999450683594;
          break;
        case 8:
          r0ab = 3.3343999385833740;
          break;
        case 9:
          r0ab = 3.1087999343872070;
          break;
        case 10:
          r0ab = 2.9597001075744629;
          break;
        case 11:
          r0ab = 3.5803999900817871;
          break;
        case 12:
          r0ab = 3.8359000682830811;
          break;
        case 13:
          r0ab = 3.7251000404357910;
          break;
        case 14:
          r0ab = 3.8940000534057617;
          break;
        case 15:
          r0ab = 3.6047000885009766;
          break;
        case 16:
          r0ab = 3.8375000953674316;
          break;
        case 17:
          r0ab = 3.5989999771118164;
          break;
        case 18:
          r0ab = 3.4328999519348145;
          break;
        case 19:
          r0ab = 3.9597001075744629;
          break;
        case 20:
          r0ab = 4.1542000770568848;
          break;
        case 21:
          r0ab = 4.0278000831604004;
          break;
        case 22:
          r0ab = 4.0047998428344727;
          break;
        case 23:
          r0ab = 3.9389998912811279;
          break;
        case 24:
          r0ab = 3.3571000099182129;
          break;
        case 25:
          r0ab = 3.8608000278472900;
          break;
        case 26:
          r0ab = 3.8055999279022217;
          break;
        case 27:
          r0ab = 3.7899000644683838;
          break;
        case 28:
          r0ab = 3.7560000419616699;
          break;
        case 29:
          r0ab = 3.7400000095367432;
          break;
        case 30:
          r0ab = 3.7081000804901123;
          break;
        case 31:
          r0ab = 3.6758000850677490;
          break;
        case 32:
          r0ab = 3.9094998836517334;
          break;
        case 33:
          r0ab = 3.6552000045776367;
          break;
        case 34:
          r0ab = 3.9572000503540039;
          break;
        case 35:
          r0ab = 3.7435998916625977;
          break;
        case 36:
          r0ab = 3.5933001041412354;
          break;
        case 37:
          r0ab = 4.1508002281188965;
          break;
        case 38:
          r0ab = 4.3337001800537109;
          break;
        case 39:
          r0ab = 4.1623997688293457;
          break;
        case 40:
          r0ab = 4.0915999412536621;
          break;
        case 41:
          r0ab = 3.6126000881195068;
          break;
        case 42:
          r0ab = 3.5581998825073242;
          break;
        case 43:
          r0ab = 3.4684000015258789;
          break;
        case 44:
          r0ab = 3.4212000370025635;
          break;
        case 45:
          r0ab = 3.4207000732421875;
          break;
        case 46:
          r0ab = 3.3828999996185303;
          break;
        case 47:
          r0ab = 3.3585999011993408;
          break;
        case 48:
          r0ab = 3.8603999614715576;
          break;
        case 49:
          r0ab = 3.8657999038696289;
          break;
        case 50:
          r0ab = 4.1156001091003418;
          break;
        case 51:
          r0ab = 3.8619999885559082;
          break;
        case 52:
          r0ab = 4.1993999481201172;
          break;
        case 53:
          r0ab = 3.9916999340057373;
          break;
        case 54:
          r0ab = 3.8445999622344971;
          break;
        case 55:
          r0ab = 4.3750000000000000;
          break;
        case 56:
          r0ab = 4.6616997718811035;
          break;
        case 57:
          r0ab = 4.3643999099731445;
          break;
        case 58:
          r0ab = 4.4082999229431152;
          break;
        case 59:
          r0ab = 4.3839998245239258;
          break;
        case 60:
          r0ab = 4.3625001907348633;
          break;
        case 61:
          r0ab = 4.3439002037048340;
          break;
        case 62:
          r0ab = 4.3253002166748047;
          break;
        case 63:
          r0ab = 4.3112001419067383;
          break;
        case 64:
          r0ab = 4.0745000839233398;
          break;
        case 65:
          r0ab = 4.1806998252868652;
          break;
        case 66:
          r0ab = 4.1578001976013184;
          break;
        case 67:
          r0ab = 4.2389998435974121;
          break;
        case 68:
          r0ab = 4.2217998504638672;
          break;
        case 69:
          r0ab = 4.2097001075744629;
          break;
        case 70:
          r0ab = 4.1985998153686523;
          break;
      }
      break;
    case 71:
      switch (atmA) {
        case 1:
          r0ab = 2.8394999504089355;
          break;
        case 2:
          r0ab = 3.0081000328063965;
          break;
        case 3:
          r0ab = 3.3171000480651855;
          break;
        case 4:
          r0ab = 3.4877998828887939;
          break;
        case 5:
          r0ab = 3.5360000133514404;
          break;
        case 6:
          r0ab = 3.5144999027252197;
          break;
        case 7:
          r0ab = 3.2809000015258789;
          break;
        case 8:
          r0ab = 3.3306999206542969;
          break;
        case 9:
          r0ab = 3.1259999275207520;
          break;
        case 10:
          r0ab = 2.9939999580383301;
          break;
        case 11:
          r0ab = 3.4741001129150391;
          break;
        case 12:
          r0ab = 3.6675000190734863;
          break;
        case 13:
          r0ab = 3.7832000255584717;
          break;
        case 14:
          r0ab = 3.6786999702453613;
          break;
        case 15:
          r0ab = 3.6156001091003418;
          break;
        case 16:
          r0ab = 3.8041000366210938;
          break;
        case 17:
          r0ab = 3.5813000202178955;
          break;
        case 18:
          r0ab = 3.4300999641418457;
          break;
        case 19:
          r0ab = 3.8480000495910645;
          break;
        case 20:
          r0ab = 3.9848999977111816;
          break;
        case 21:
          r0ab = 3.9314000606536865;
          break;
        case 22:
          r0ab = 3.8405001163482666;
          break;
        case 23:
          r0ab = 3.8029000759124756;
          break;
        case 24:
          r0ab = 3.2962000370025635;
          break;
        case 25:
          r0ab = 3.7104001045227051;
          break;
        case 26:
          r0ab = 3.6514999866485596;
          break;
        case 27:
          r0ab = 3.6377999782562256;
          break;
        case 28:
          r0ab = 3.6019999980926514;
          break;
        case 29:
          r0ab = 3.5848999023437500;
          break;
        case 30:
          r0ab = 3.5550000667572021;
          break;
        case 31:
          r0ab = 3.7493999004364014;
          break;
        case 32:
          r0ab = 3.6893000602722168;
          break;
        case 33:
          r0ab = 3.6665999889373779;
          break;
        case 34:
          r0ab = 3.9170000553131104;
          break;
        case 35:
          r0ab = 3.7149999141693115;
          break;
        case 36:
          r0ab = 3.5759999752044678;
          break;
        case 37:
          r0ab = 4.0268001556396484;
          break;
        case 38:
          r0ab = 4.1595997810363770;
          break;
        case 39:
          r0ab = 4.1107001304626465;
          break;
        case 40:
          r0ab = 3.9995000362396240;
          break;
        case 41:
          r0ab = 3.5573999881744385;
          break;
        case 42:
          r0ab = 3.5102999210357666;
          break;
        case 43:
          r0ab = 3.4163000583648682;
          break;
        case 44:
          r0ab = 3.3654999732971191;
          break;
        case 45:
          r0ab = 3.3677000999450684;
          break;
        case 46:
          r0ab = 3.3243000507354736;
          break;
        case 47:
          r0ab = 3.2974998950958252;
          break;
        case 48:
          r0ab = 3.7070999145507812;
          break;
        case 49:
          r0ab = 3.9047000408172607;
          break;
        case 50:
          r0ab = 3.8513998985290527;
          break;
        case 51:
          r0ab = 3.8422000408172607;
          break;
        case 52:
          r0ab = 3.8022000789642334;
          break;
        case 53:
          r0ab = 3.9323000907897949;
          break;
        case 54:
          r0ab = 3.7932000160217285;
          break;
        case 55:
          r0ab = 4.2343001365661621;
          break;
        case 56:
          r0ab = 4.4583001136779785;
          break;
        case 57:
          r0ab = 4.3115000724792480;
          break;
        case 58:
          r0ab = 4.2456998825073242;
          break;
        case 59:
          r0ab = 4.2213001251220703;
          break;
        case 60:
          r0ab = 4.1944999694824219;
          break;
        case 61:
          r0ab = 4.1756000518798828;
          break;
        case 62:
          r0ab = 4.1568999290466309;
          break;
        case 63:
          r0ab = 4.1423997879028320;
          break;
        case 64:
          r0ab = 4.0619997978210449;
          break;
        case 65:
          r0ab = 4.0493998527526855;
          break;
        case 66:
          r0ab = 3.9953000545501709;
          break;
        case 67:
          r0ab = 4.0693998336791992;
          break;
        case 68:
          r0ab = 4.0515999794006348;
          break;
        case 69:
          r0ab = 4.0395998954772949;
          break;
        case 70:
          r0ab = 4.0279998779296875;
          break;
        case 71:
          r0ab = 4.0130000114440918;
          break;
      }
      break;
    case 72:
      switch (atmA) {
        case 1:
          r0ab = 2.9007000923156738;
          break;
        case 2:
          r0ab = 2.9674000740051270;
          break;
        case 3:
          r0ab = 3.8173999786376953;
          break;
        case 4:
          r0ab = 3.5855998992919922;
          break;
        case 5:
          r0ab = 3.6486001014709473;
          break;
        case 6:
          r0ab = 3.5339000225067139;
          break;
        case 7:
          r0ab = 3.2832000255584717;
          break;
        case 8:
          r0ab = 3.3153998851776123;
          break;
        case 9:
          r0ab = 3.1143999099731445;
          break;
        case 10:
          r0ab = 2.9865999221801758;
          break;
        case 11:
          r0ab = 3.9618000984191895;
          break;
        case 12:
          r0ab = 3.8429999351501465;
          break;
        case 13:
          r0ab = 3.9979999065399170;
          break;
        case 14:
          r0ab = 3.8134000301361084;
          break;
        case 15:
          r0ab = 3.6651999950408936;
          break;
        case 16:
          r0ab = 3.7985000610351562;
          break;
        case 17:
          r0ab = 3.5755999088287354;
          break;
        case 18:
          r0ab = 3.4207000732421875;
          break;
        case 19:
          r0ab = 4.4060997962951660;
          break;
        case 20:
          r0ab = 4.2817001342773438;
          break;
        case 21:
          r0ab = 4.1476998329162598;
          break;
        case 22:
          r0ab = 4.0616002082824707;
          break;
        case 23:
          r0ab = 3.9979000091552734;
          break;
        case 24:
          r0ab = 3.6491999626159668;
          break;
        case 25:
          r0ab = 3.8833000659942627;
          break;
        case 26:
          r0ab = 3.8027000427246094;
          break;
        case 27:
          r0ab = 3.7660000324249268;
          break;
        case 28:
          r0ab = 3.7183001041412354;
          break;
        case 29:
          r0ab = 3.6953999996185303;
          break;
        case 30:
          r0ab = 3.6524999141693115;
          break;
        case 31:
          r0ab = 3.9669001102447510;
          break;
        case 32:
          r0ab = 3.8371000289916992;
          break;
        case 33:
          r0ab = 3.7325000762939453;
          break;
        case 34:
          r0ab = 3.9159998893737793;
          break;
        case 35:
          r0ab = 3.7156000137329102;
          break;
        case 36:
          r0ab = 3.5713999271392822;
          break;
        case 37:
          r0ab = 4.6036000251770020;
          break;
        case 38:
          r0ab = 4.4619998931884766;
          break;
        case 39:
          r0ab = 4.3091998100280762;
          break;
        case 40:
          r0ab = 4.2122001647949219;
          break;
        case 41:
          r0ab = 3.8478000164031982;
          break;
        case 42:
          r0ab = 3.7572000026702881;
          break;
        case 43:
          r0ab = 3.6596999168395996;
          break;
        case 44:
          r0ab = 3.5968999862670898;
          break;
        case 45:
          r0ab = 3.5574998855590820;
          break;
        case 46:
          r0ab = 3.5385999679565430;
          break;
        case 47:
          r0ab = 3.5153000354766846;
          break;
        case 48:
          r0ab = 3.7818000316619873;
          break;
        case 49:
          r0ab = 4.1335000991821289;
          break;
        case 50:
          r0ab = 4.0152997970581055;
          break;
        case 51:
          r0ab = 3.9177000522613525;
          break;
        case 52:
          r0ab = 3.8603000640869141;
          break;
        case 53:
          r0ab = 3.9365000724792480;
          break;
        case 54:
          r0ab = 3.7906000614166260;
          break;
        case 55:
          r0ab = 4.7936000823974609;
          break;
        case 56:
          r0ab = 4.7410001754760742;
          break;
        case 57:
          r0ab = 4.5461001396179199;
          break;
        case 58:
          r0ab = 4.5661997795104980;
          break;
        case 59:
          r0ab = 4.5339999198913574;
          break;
        case 60:
          r0ab = 4.5058999061584473;
          break;
        case 61:
          r0ab = 4.4832000732421875;
          break;
        case 62:
          r0ab = 4.4604001045227051;
          break;
        case 63:
          r0ab = 4.4429001808166504;
          break;
        case 64:
          r0ab = 4.2346000671386719;
          break;
        case 65:
          r0ab = 4.4204001426696777;
          break;
        case 66:
          r0ab = 4.3119001388549805;
          break;
        case 67:
          r0ab = 4.3449997901916504;
          break;
        case 68:
          r0ab = 4.3193001747131348;
          break;
        case 69:
          r0ab = 4.3035001754760742;
          break;
        case 70:
          r0ab = 4.2933001518249512;
          break;
        case 71:
          r0ab = 4.1581997871398926;
          break;
        case 72:
          r0ab = 4.2449998855590820;
          break;
      }
      break;
    case 73:
      switch (atmA) {
        case 1:
          r0ab = 2.8559000492095947;
          break;
        case 2:
          r0ab = 2.9049999713897705;
          break;
        case 3:
          r0ab = 3.8324999809265137;
          break;
        case 4:
          r0ab = 3.5441999435424805;
          break;
        case 5:
          r0ab = 3.5076999664306641;
          break;
        case 6:
          r0ab = 3.4904999732971191;
          break;
        case 7:
          r0ab = 3.2395999431610107;
          break;
        case 8:
          r0ab = 3.2720000743865967;
          break;
        case 9:
          r0ab = 3.0725998878479004;
          break;
        case 10:
          r0ab = 2.9467000961303711;
          break;
        case 11:
          r0ab = 3.9644000530242920;
          break;
        case 12:
          r0ab = 3.8050000667572021;
          break;
        case 13:
          r0ab = 3.8980998992919922;
          break;
        case 14:
          r0ab = 3.7762000560760498;
          break;
        case 15:
          r0ab = 3.6215999126434326;
          break;
        case 16:
          r0ab = 3.7530999183654785;
          break;
        case 17:
          r0ab = 3.5297000408172607;
          break;
        case 18:
          r0ab = 3.3742001056671143;
          break;
        case 19:
          r0ab = 4.3814001083374023;
          break;
        case 20:
          r0ab = 4.2817997932434082;
          break;
        case 21:
          r0ab = 4.1026000976562500;
          break;
        case 22:
          r0ab = 4.0293998718261719;
          break;
        case 23:
          r0ab = 3.9639999866485596;
          break;
        case 24:
          r0ab = 3.6208000183105469;
          break;
        case 25:
          r0ab = 3.8464000225067139;
          break;
        case 26:
          r0ab = 3.7648000717163086;
          break;
        case 27:
          r0ab = 3.7281000614166260;
          break;
        case 28:
          r0ab = 3.6789999008178711;
          break;
        case 29:
          r0ab = 3.6542000770568848;
          break;
        case 30:
          r0ab = 3.6117000579833984;
          break;
        case 31:
          r0ab = 3.8650000095367432;
          break;
        case 32:
          r0ab = 3.8010001182556152;
          break;
        case 33:
          r0ab = 3.6893999576568604;
          break;
        case 34:
          r0ab = 3.8712999820709229;
          break;
        case 35:
          r0ab = 3.6698999404907227;
          break;
        case 36:
          r0ab = 3.5243999958038330;
          break;
        case 37:
          r0ab = 4.5151000022888184;
          break;
        case 38:
          r0ab = 4.4517002105712891;
          break;
        case 39:
          r0ab = 4.2537999153137207;
          break;
        case 40:
          r0ab = 4.1483001708984375;
          break;
        case 41:
          r0ab = 3.8640999794006348;
          break;
        case 42:
          r0ab = 3.7244000434875488;
          break;
        case 43:
          r0ab = 3.6243000030517578;
          break;
        case 44:
          r0ab = 3.5589001178741455;
          break;
        case 45:
          r0ab = 3.5171999931335449;
          break;
        case 46:
          r0ab = 3.4972999095916748;
          break;
        case 47:
          r0ab = 3.4714999198913574;
          break;
        case 48:
          r0ab = 3.7339999675750732;
          break;
        case 49:
          r0ab = 4.0315999984741211;
          break;
        case 50:
          r0ab = 3.9958000183105469;
          break;
        case 51:
          r0ab = 3.8687000274658203;
          break;
        case 52:
          r0ab = 3.8115000724792480;
          break;
        case 53:
          r0ab = 3.8861999511718750;
          break;
        case 54:
          r0ab = 3.7379000186920166;
          break;
        case 55:
          r0ab = 4.7090997695922852;
          break;
        case 56:
          r0ab = 4.7156000137329102;
          break;
        case 57:
          r0ab = 4.5198998451232910;
          break;
        case 58:
          r0ab = 4.5542001724243164;
          break;
        case 59:
          r0ab = 4.5229997634887695;
          break;
        case 60:
          r0ab = 4.4959001541137695;
          break;
        case 61:
          r0ab = 4.4749999046325684;
          break;
        case 62:
          r0ab = 4.4528999328613281;
          break;
        case 63:
          r0ab = 4.4361000061035156;
          break;
        case 64:
          r0ab = 4.1774001121520996;
          break;
        case 65:
          r0ab = 4.3773999214172363;
          break;
        case 66:
          r0ab = 4.2962999343872070;
          break;
        case 67:
          r0ab = 4.3406000137329102;
          break;
        case 68:
          r0ab = 4.3158998489379883;
          break;
        case 69:
          r0ab = 4.3006000518798828;
          break;
        case 70:
          r0ab = 4.2909998893737793;
          break;
        case 71:
          r0ab = 4.1008000373840332;
          break;
        case 72:
          r0ab = 4.1567997932434082;
          break;
        case 73:
          r0ab = 4.0980000495910645;
          break;
      }
      break;
    case 74:
      switch (atmA) {
        case 1:
          r0ab = 2.8110001087188721;
          break;
        case 2:
          r0ab = 2.8519999980926514;
          break;
        case 3:
          r0ab = 3.7479999065399170;
          break;
        case 4:
          r0ab = 3.5104999542236328;
          break;
        case 5:
          r0ab = 3.4346001148223877;
          break;
        case 6:
          r0ab = 3.3461000919342041;
          break;
        case 7:
          r0ab = 3.1970999240875244;
          break;
        case 8:
          r0ab = 3.2325999736785889;
          break;
        case 9:
          r0ab = 3.0329000949859619;
          break;
        case 10:
          r0ab = 2.9070000648498535;
          break;
        case 11:
          r0ab = 3.8822999000549316;
          break;
        case 12:
          r0ab = 3.7927999496459961;
          break;
        case 13:
          r0ab = 3.8264000415802002;
          break;
        case 14:
          r0ab = 3.7005999088287354;
          break;
        case 15:
          r0ab = 3.5796999931335449;
          break;
        case 16:
          r0ab = 3.7140998840332031;
          break;
        case 17:
          r0ab = 3.4893999099731445;
          break;
        case 18:
          r0ab = 3.3326001167297363;
          break;
        case 19:
          r0ab = 4.3048000335693359;
          break;
        case 20:
          r0ab = 4.2217001914978027;
          break;
        case 21:
          r0ab = 4.0785999298095703;
          break;
        case 22:
          r0ab = 3.9900000095367432;
          break;
        case 23:
          r0ab = 3.9356999397277832;
          break;
        case 24:
          r0ab = 3.6331000328063965;
          break;
        case 25:
          r0ab = 3.8333001136779785;
          break;
        case 26:
          r0ab = 3.7316999435424805;
          break;
        case 27:
          r0ab = 3.6956999301910400;
          break;
        case 28:
          r0ab = 3.6459999084472656;
          break;
        case 29:
          r0ab = 3.6196999549865723;
          break;
        case 30:
          r0ab = 3.5778999328613281;
          break;
        case 31:
          r0ab = 3.7908999919891357;
          break;
        case 32:
          r0ab = 3.7256999015808105;
          break;
        case 33:
          r0ab = 3.6475999355316162;
          break;
        case 34:
          r0ab = 3.5729000568389893;
          break;
        case 35:
          r0ab = 3.6303999423980713;
          break;
        case 36:
          r0ab = 3.4834001064300537;
          break;
        case 37:
          r0ab = 4.4368000030517578;
          break;
        case 38:
          r0ab = 4.3920998573303223;
          break;
        case 39:
          r0ab = 4.2206997871398926;
          break;
        case 40:
          r0ab = 4.1132998466491699;
          break;
        case 41:
          r0ab = 3.8066999912261963;
          break;
        case 42:
          r0ab = 3.7421000003814697;
          break;
        case 43:
          r0ab = 3.6140000820159912;
          break;
        case 44:
          r0ab = 3.5490999221801758;
          break;
        case 45:
          r0ab = 3.5076999664306641;
          break;
        case 46:
          r0ab = 3.4886999130249023;
          break;
        case 47:
          r0ab = 3.4623000621795654;
          break;
        case 48:
          r0ab = 3.6956000328063965;
          break;
        case 49:
          r0ab = 3.9567999839782715;
          break;
        case 50:
          r0ab = 3.8975999355316162;
          break;
        case 51:
          r0ab = 3.8239998817443848;
          break;
        case 52:
          r0ab = 3.7683999538421631;
          break;
        case 53:
          r0ab = 3.8450999259948730;
          break;
        case 54:
          r0ab = 3.6949000358581543;
          break;
        case 55:
          r0ab = 4.6318001747131348;
          break;
        case 56:
          r0ab = 4.6559000015258789;
          break;
        case 57:
          r0ab = 4.4532999992370605;
          break;
        case 58:
          r0ab = 4.4956002235412598;
          break;
        case 59:
          r0ab = 4.4640998840332031;
          break;
        case 60:
          r0ab = 4.4366002082824707;
          break;
        case 61:
          r0ab = 4.4155001640319824;
          break;
        case 62:
          r0ab = 4.3935999870300293;
          break;
        case 63:
          r0ab = 4.3763999938964844;
          break;
        case 64:
          r0ab = 4.1301999092102051;
          break;
        case 65:
          r0ab = 4.3397998809814453;
          break;
        case 66:
          r0ab = 4.2283000946044922;
          break;
        case 67:
          r0ab = 4.2796001434326172;
          break;
        case 68:
          r0ab = 4.2547001838684082;
          break;
        case 69:
          r0ab = 4.2390999794006348;
          break;
        case 70:
          r0ab = 4.2295999526977539;
          break;
        case 71:
          r0ab = 4.0699000358581543;
          break;
        case 72:
          r0ab = 4.1083002090454102;
          break;
        case 73:
          r0ab = 4.0318999290466309;
          break;
        case 74:
          r0ab = 3.9855000972747803;
          break;
      }
      break;
    case 75:
      switch (atmA) {
        case 1:
          r0ab = 2.7676000595092773;
          break;
        case 2:
          r0ab = 2.8078000545501709;
          break;
        case 3:
          r0ab = 3.6724998950958252;
          break;
        case 4:
          r0ab = 3.4804000854492188;
          break;
        case 5:
          r0ab = 3.3775000572204590;
          break;
        case 6:
          r0ab = 3.2411000728607178;
          break;
        case 7:
          r0ab = 3.1580998897552490;
          break;
        case 8:
          r0ab = 3.1982998847961426;
          break;
        case 9:
          r0ab = 2.9972999095916748;
          break;
        case 10:
          r0ab = 2.8705000877380371;
          break;
        case 11:
          r0ab = 3.8069999217987061;
          break;
        case 12:
          r0ab = 3.7392001152038574;
          break;
        case 13:
          r0ab = 3.7667999267578125;
          break;
        case 14:
          r0ab = 3.6263000965118408;
          break;
        case 15:
          r0ab = 3.5401999950408936;
          break;
        case 16:
          r0ab = 3.6807000637054443;
          break;
        case 17:
          r0ab = 3.4544999599456787;
          break;
        case 18:
          r0ab = 3.2962000370025635;
          break;
        case 19:
          r0ab = 4.2283000946044922;
          break;
        case 20:
          r0ab = 4.1697998046875000;
          break;
        case 21:
          r0ab = 4.0240001678466797;
          break;
        case 22:
          r0ab = 3.9340999126434326;
          break;
        case 23:
          r0ab = 3.8710999488830566;
          break;
        case 24:
          r0ab = 3.5488998889923096;
          break;
        case 25:
          r0ab = 3.7797999382019043;
          break;
        case 26:
          r0ab = 3.7000000476837158;
          break;
        case 27:
          r0ab = 3.6654000282287598;
          break;
        case 28:
          r0ab = 3.6154000759124756;
          break;
        case 29:
          r0ab = 3.5882000923156738;
          break;
        case 30:
          r0ab = 3.5471999645233154;
          break;
        case 31:
          r0ab = 3.7288999557495117;
          break;
        case 32:
          r0ab = 3.6510000228881836;
          break;
        case 33:
          r0ab = 3.6078000068664551;
          break;
        case 34:
          r0ab = 3.5355000495910645;
          break;
        case 35:
          r0ab = 3.5962998867034912;
          break;
        case 36:
          r0ab = 3.4479999542236328;
          break;
        case 37:
          r0ab = 4.3586997985839844;
          break;
        case 38:
          r0ab = 4.3390002250671387;
          break;
        case 39:
          r0ab = 4.1634998321533203;
          break;
        case 40:
          r0ab = 4.0535998344421387;
          break;
        case 41:
          r0ab = 3.7193000316619873;
          break;
        case 42:
          r0ab = 3.6528999805450439;
          break;
        case 43:
          r0ab = 3.5511999130249023;
          break;
        case 44:
          r0ab = 3.4837000370025635;
          break;
        case 45:
          r0ab = 3.4400000572204590;
          break;
        case 46:
          r0ab = 3.4191000461578369;
          break;
        case 47:
          r0ab = 3.3891000747680664;
          break;
        case 48:
          r0ab = 3.6621999740600586;
          break;
        case 49:
          r0ab = 3.8933999538421631;
          break;
        case 50:
          r0ab = 3.8234999179840088;
          break;
        case 51:
          r0ab = 3.7822999954223633;
          break;
        case 52:
          r0ab = 3.7291998863220215;
          break;
        case 53:
          r0ab = 3.8106000423431396;
          break;
        case 54:
          r0ab = 3.6589000225067139;
          break;
        case 55:
          r0ab = 4.5535001754760742;
          break;
        case 56:
          r0ab = 4.6012997627258301;
          break;
        case 57:
          r0ab = 4.3961000442504883;
          break;
        case 58:
          r0ab = 4.4422998428344727;
          break;
        case 59:
          r0ab = 4.4109001159667969;
          break;
        case 60:
          r0ab = 4.3835000991821289;
          break;
        case 61:
          r0ab = 4.3625001907348633;
          break;
        case 62:
          r0ab = 4.3407001495361328;
          break;
        case 63:
          r0ab = 4.3236999511718750;
          break;
        case 64:
          r0ab = 4.0862998962402344;
          break;
        case 65:
          r0ab = 4.2835001945495605;
          break;
        case 66:
          r0ab = 4.1675000190734863;
          break;
        case 67:
          r0ab = 4.2272000312805176;
          break;
        case 68:
          r0ab = 4.2024998664855957;
          break;
        case 69:
          r0ab = 4.1869001388549805;
          break;
        case 70:
          r0ab = 4.1774001121520996;
          break;
        case 71:
          r0ab = 4.0125999450683594;
          break;
        case 72:
          r0ab = 4.0460000038146973;
          break;
        case 73:
          r0ab = 3.9814999103546143;
          break;
        case 74:
          r0ab = 3.9340000152587891;
          break;
        case 75:
          r0ab = 3.8954999446868896;
          break;
      }
      break;
    case 76:
      switch (atmA) {
        case 1:
          r0ab = 2.6912000179290771;
          break;
        case 2:
          r0ab = 2.7604000568389893;
          break;
        case 3:
          r0ab = 3.6036999225616455;
          break;
        case 4:
          r0ab = 3.4193999767303467;
          break;
        case 5:
          r0ab = 3.3094000816345215;
          break;
        case 6:
          r0ab = 3.1710000038146973;
          break;
        case 7:
          r0ab = 3.0861999988555908;
          break;
        case 8:
          r0ab = 3.1789000034332275;
          break;
        case 9:
          r0ab = 2.9737999439239502;
          break;
        case 10:
          r0ab = 2.8427000045776367;
          break;
        case 11:
          r0ab = 3.7377998828887939;
          break;
        case 12:
          r0ab = 3.6742000579833984;
          break;
        case 13:
          r0ab = 3.6928000450134277;
          break;
        case 14:
          r0ab = 3.5511999130249023;
          break;
        case 15:
          r0ab = 3.4614000320434570;
          break;
        case 16:
          r0ab = 3.4086999893188477;
          break;
        case 17:
          r0ab = 3.4200999736785889;
          break;
        case 18:
          r0ab = 3.2606999874114990;
          break;
        case 19:
          r0ab = 4.1526999473571777;
          break;
        case 20:
          r0ab = 4.0977001190185547;
          break;
        case 21:
          r0ab = 3.9523000717163086;
          break;
        case 22:
          r0ab = 3.8627998828887939;
          break;
        case 23:
          r0ab = 3.8001999855041504;
          break;
        case 24:
          r0ab = 3.4758999347686768;
          break;
        case 25:
          r0ab = 3.7102000713348389;
          break;
        case 26:
          r0ab = 3.6466000080108643;
          break;
        case 27:
          r0ab = 3.6105999946594238;
          break;
        case 28:
          r0ab = 3.5580000877380371;
          break;
        case 29:
          r0ab = 3.5281999111175537;
          break;
        case 30:
          r0ab = 3.4877998828887939;
          break;
        case 31:
          r0ab = 3.6547000408172607;
          break;
        case 32:
          r0ab = 3.5762999057769775;
          break;
        case 33:
          r0ab = 3.5288999080657959;
          break;
        case 34:
          r0ab = 3.5085999965667725;
          break;
        case 35:
          r0ab = 3.5592999458312988;
          break;
        case 36:
          r0ab = 3.4098999500274658;
          break;
        case 37:
          r0ab = 4.2788000106811523;
          break;
        case 38:
          r0ab = 4.2624001502990723;
          break;
        case 39:
          r0ab = 4.0872998237609863;
          break;
        case 40:
          r0ab = 3.9769999980926514;
          break;
        case 41:
          r0ab = 3.6407001018524170;
          break;
        case 42:
          r0ab = 3.5743000507354736;
          break;
        case 43:
          r0ab = 3.5178000926971436;
          break;
        case 44:
          r0ab = 3.4753000736236572;
          break;
        case 45:
          r0ab = 3.3931000232696533;
          break;
        case 46:
          r0ab = 3.3694000244140625;
          break;
        case 47:
          r0ab = 3.3338999748229980;
          break;
        case 48:
          r0ab = 3.6001999378204346;
          break;
        case 49:
          r0ab = 3.8164000511169434;
          break;
        case 50:
          r0ab = 3.7478001117706299;
          break;
        case 51:
          r0ab = 3.7028000354766846;
          break;
        case 52:
          r0ab = 3.6951999664306641;
          break;
        case 53:
          r0ab = 3.7669000625610352;
          break;
        case 54:
          r0ab = 3.6136999130249023;
          break;
        case 55:
          r0ab = 4.4697999954223633;
          break;
        case 56:
          r0ab = 4.5487999916076660;
          break;
        case 57:
          r0ab = 4.3168001174926758;
          break;
        case 58:
          r0ab = 4.3646001815795898;
          break;
        case 59:
          r0ab = 4.3337998390197754;
          break;
        case 60:
          r0ab = 4.3067002296447754;
          break;
        case 61:
          r0ab = 4.2859997749328613;
          break;
        case 62:
          r0ab = 4.2645001411437988;
          break;
        case 63:
          r0ab = 4.2477998733520508;
          break;
        case 64:
          r0ab = 4.0067000389099121;
          break;
        case 65:
          r0ab = 4.2348999977111816;
          break;
        case 66:
          r0ab = 4.0957999229431152;
          break;
        case 67:
          r0ab = 4.1543002128601074;
          break;
        case 68:
          r0ab = 4.1301999092102051;
          break;
        case 69:
          r0ab = 4.1140999794006348;
          break;
        case 70:
          r0ab = 4.1048002243041992;
          break;
        case 71:
          r0ab = 3.9409999847412109;
          break;
        case 72:
          r0ab = 3.9595000743865967;
          break;
        case 73:
          r0ab = 3.8940999507904053;
          break;
        case 74:
          r0ab = 3.8464999198913574;
          break;
        case 75:
          r0ab = 3.8089001178741455;
          break;
        case 76:
          r0ab = 3.7490000724792480;
          break;
      }
      break;
    case 77:
      switch (atmA) {
        case 1:
          r0ab = 2.7894999980926514;
          break;
        case 2:
          r0ab = 2.5848999023437500;
          break;
        case 3:
          r0ab = 3.6484000682830811;
          break;
        case 4:
          r0ab = 3.0162000656127930;
          break;
        case 5:
          r0ab = 3.1266999244689941;
          break;
        case 6:
          r0ab = 3.2125000953674316;
          break;
        case 7:
          r0ab = 3.0043001174926758;
          break;
        case 8:
          r0ab = 2.9572000503540039;
          break;
        case 9:
          r0ab = 2.8197000026702881;
          break;
        case 10:
          r0ab = 2.7260999679565430;
          break;
        case 11:
          r0ab = 3.7701001167297363;
          break;
        case 12:
          r0ab = 3.2446000576019287;
          break;
        case 13:
          r0ab = 3.5239000320434570;
          break;
        case 14:
          r0ab = 3.4695999622344971;
          break;
        case 15:
          r0ab = 3.4261000156402588;
          break;
        case 16:
          r0ab = 3.3508000373840332;
          break;
        case 17:
          r0ab = 3.1967999935150146;
          break;
        case 18:
          r0ab = 3.0848000049591064;
          break;
        case 19:
          r0ab = 4.1496000289916992;
          break;
        case 20:
          r0ab = 3.6598000526428223;
          break;
        case 21:
          r0ab = 3.5111000537872314;
          break;
        case 22:
          r0ab = 3.4198999404907227;
          break;
        case 23:
          r0ab = 3.3808999061584473;
          break;
        case 24:
          r0ab = 3.5381999015808105;
          break;
        case 25:
          r0ab = 3.2572000026702881;
          break;
        case 26:
          r0ab = 3.2100000381469727;
          break;
        case 27:
          r0ab = 3.1916999816894531;
          break;
        case 28:
          r0ab = 3.1519000530242920;
          break;
        case 29:
          r0ab = 3.1198000907897949;
          break;
        case 30:
          r0ab = 3.1005001068115234;
          break;
        case 31:
          r0ab = 3.5071001052856445;
          break;
        case 32:
          r0ab = 3.5085999965667725;
          break;
        case 33:
          r0ab = 3.5072999000549316;
          break;
        case 34:
          r0ab = 3.4509000778198242;
          break;
        case 35:
          r0ab = 3.3120000362396240;
          break;
        case 36:
          r0ab = 3.2081999778747559;
          break;
        case 37:
          r0ab = 4.2610998153686523;
          break;
        case 38:
          r0ab = 3.8117001056671143;
          break;
        case 39:
          r0ab = 3.6988000869750977;
          break;
        case 40:
          r0ab = 3.5645999908447266;
          break;
        case 41:
          r0ab = 3.6925001144409180;
          break;
        case 42:
          r0ab = 3.6294999122619629;
          break;
        case 43:
          r0ab = 3.5383000373840332;
          break;
        case 44:
          r0ab = 3.4909999370574951;
          break;
        case 45:
          r0ab = 3.4625000953674316;
          break;
        case 46:
          r0ab = 3.4233000278472900;
          break;
        case 47:
          r0ab = 3.4007000923156738;
          break;
        case 48:
          r0ab = 3.2328999042510986;
          break;
        case 49:
          r0ab = 3.6723001003265381;
          break;
        case 50:
          r0ab = 3.6844999790191650;
          break;
        case 51:
          r0ab = 3.6875998973846436;
          break;
        case 52:
          r0ab = 3.6196999549865723;
          break;
        case 53:
          r0ab = 3.4798998832702637;
          break;
        case 54:
          r0ab = 3.3736999034881592;
          break;
        case 55:
          r0ab = 4.4341001510620117;
          break;
        case 56:
          r0ab = 4.0524997711181641;
          break;
        case 57:
          r0ab = 3.9010999202728271;
          break;
        case 58:
          r0ab = 3.8945000171661377;
          break;
        case 59:
          r0ab = 3.8635001182556152;
          break;
        case 60:
          r0ab = 3.8368000984191895;
          break;
        case 61:
          r0ab = 3.8152999877929688;
          break;
        case 62:
          r0ab = 3.7936000823974609;
          break;
        case 63:
          r0ab = 3.7757999897003174;
          break;
        case 64:
          r0ab = 3.4944000244140625;
          break;
        case 65:
          r0ab = 3.4872999191284180;
          break;
        case 66:
          r0ab = 3.9040000438690186;
          break;
        case 67:
          r0ab = 3.7109999656677246;
          break;
        case 68:
          r0ab = 3.6921999454498291;
          break;
        case 69:
          r0ab = 3.6798999309539795;
          break;
        case 70:
          r0ab = 3.6723999977111816;
          break;
        case 71:
          r0ab = 3.5622000694274902;
          break;
        case 72:
          r0ab = 3.6080999374389648;
          break;
        case 73:
          r0ab = 3.5425999164581299;
          break;
        case 74:
          r0ab = 3.4921998977661133;
          break;
        case 75:
          r0ab = 3.4498000144958496;
          break;
        case 76:
          r0ab = 3.3984000682830811;
          break;
        case 77:
          r0ab = 3.4456000328063965;
          break;
      }
      break;
    case 78:
      switch (atmA) {
        case 1:
          r0ab = 2.7521998882293701;
          break;
        case 2:
          r0ab = 2.5524001121520996;
          break;
        case 3:
          r0ab = 3.5741999149322510;
          break;
        case 4:
          r0ab = 2.9507999420166016;
          break;
        case 5:
          r0ab = 3.0750999450683594;
          break;
        case 6:
          r0ab = 3.0157999992370605;
          break;
        case 7:
          r0ab = 2.9644000530242920;
          break;
        case 8:
          r0ab = 2.8338000774383545;
          break;
        case 9:
          r0ab = 2.7890999317169189;
          break;
        case 10:
          r0ab = 2.6933000087738037;
          break;
        case 11:
          r0ab = 3.6926000118255615;
          break;
        case 12:
          r0ab = 3.1814000606536865;
          break;
        case 13:
          r0ab = 3.4528000354766846;
          break;
        case 14:
          r0ab = 3.4186000823974609;
          break;
        case 15:
          r0ab = 3.3835999965667725;
          break;
        case 16:
          r0ab = 3.2212998867034912;
          break;
        case 17:
          r0ab = 3.1626000404357910;
          break;
        case 18:
          r0ab = 3.0506999492645264;
          break;
        case 19:
          r0ab = 4.0548000335693359;
          break;
        case 20:
          r0ab = 3.5311999320983887;
          break;
        case 21:
          r0ab = 3.4244000911712646;
          break;
        case 22:
          r0ab = 3.3408999443054199;
          break;
        case 23:
          r0ab = 3.2809998989105225;
          break;
        case 24:
          r0ab = 3.4781999588012695;
          break;
        case 25:
          r0ab = 3.1905000209808350;
          break;
        case 26:
          r0ab = 3.1493999958038330;
          break;
        case 27:
          r0ab = 3.1221001148223877;
          break;
        case 28:
          r0ab = 3.1127998828887939;
          break;
        case 29:
          r0ab = 3.0852999687194824;
          break;
        case 30:
          r0ab = 3.0383999347686768;
          break;
        case 31:
          r0ab = 3.4365999698638916;
          break;
        case 32:
          r0ab = 3.4561998844146729;
          break;
        case 33:
          r0ab = 3.4637999534606934;
          break;
        case 34:
          r0ab = 3.3210999965667725;
          break;
        case 35:
          r0ab = 3.2762000560760498;
          break;
        case 36:
          r0ab = 3.1730000972747803;
          break;
        case 37:
          r0ab = 4.1631999015808105;
          break;
        case 38:
          r0ab = 3.6824998855590820;
          break;
        case 39:
          r0ab = 3.5822000503540039;
          break;
        case 40:
          r0ab = 3.4869999885559082;
          break;
        case 41:
          r0ab = 3.6324999332427979;
          break;
        case 42:
          r0ab = 3.5739998817443848;
          break;
        case 43:
          r0ab = 3.4732999801635742;
          break;
        case 44:
          r0ab = 3.4247000217437744;
          break;
        case 45:
          r0ab = 3.3968999385833740;
          break;
        case 46:
          r0ab = 3.3763999938964844;
          break;
        case 47:
          r0ab = 3.3524999618530273;
          break;
        case 48:
          r0ab = 3.1984000205993652;
          break;
        case 49:
          r0ab = 3.5989000797271729;
          break;
        case 50:
          r0ab = 3.6298999786376953;
          break;
        case 51:
          r0ab = 3.6433000564575195;
          break;
        case 52:
          r0ab = 3.4937000274658203;
          break;
        case 53:
          r0ab = 3.4416999816894531;
          break;
        case 54:
          r0ab = 3.3364999294281006;
          break;
        case 55:
          r0ab = 4.3303999900817871;
          break;
        case 56:
          r0ab = 3.9242000579833984;
          break;
        case 57:
          r0ab = 3.7792999744415283;
          break;
        case 58:
          r0ab = 3.7623000144958496;
          break;
        case 59:
          r0ab = 3.7327001094818115;
          break;
        case 60:
          r0ab = 3.7070999145507812;
          break;
        case 61:
          r0ab = 3.6860001087188721;
          break;
        case 62:
          r0ab = 3.6649999618530273;
          break;
        case 63:
          r0ab = 3.6475999355316162;
          break;
        case 64:
          r0ab = 3.3849000930786133;
          break;
        case 65:
          r0ab = 3.3533999919891357;
          break;
        case 66:
          r0ab = 3.8215999603271484;
          break;
        case 67:
          r0ab = 3.5869998931884766;
          break;
        case 68:
          r0ab = 3.5694999694824219;
          break;
        case 69:
          r0ab = 3.5583999156951904;
          break;
        case 70:
          r0ab = 3.5508000850677490;
          break;
        case 71:
          r0ab = 3.4855999946594238;
          break;
        case 72:
          r0ab = 3.5522999763488770;
          break;
        case 73:
          r0ab = 3.4934000968933105;
          break;
        case 74:
          r0ab = 3.4463999271392822;
          break;
        case 75:
          r0ab = 3.4054999351501465;
          break;
        case 76:
          r0ab = 3.3550999164581299;
          break;
        case 77:
          r0ab = 3.3887999057769775;
          break;
        case 78:
          r0ab = 3.3524999618530273;
          break;
      }
      break;
    case 79:
      switch (atmA) {
        case 1:
          r0ab = 2.7202000617980957;
          break;
        case 2:
          r0ab = 2.5183000564575195;
          break;
        case 3:
          r0ab = 3.4946999549865723;
          break;
        case 4:
          r0ab = 2.8731000423431396;
          break;
        case 5:
          r0ab = 3.0197999477386475;
          break;
        case 6:
          r0ab = 3.1456999778747559;
          break;
        case 7:
          r0ab = 2.9275999069213867;
          break;
        case 8:
          r0ab = 2.7825999259948730;
          break;
        case 9:
          r0ab = 2.7574000358581543;
          break;
        case 10:
          r0ab = 2.6605999469757080;
          break;
        case 11:
          r0ab = 3.6089999675750732;
          break;
        case 12:
          r0ab = 3.0580999851226807;
          break;
        case 13:
          r0ab = 3.3747000694274902;
          break;
        case 14:
          r0ab = 3.3677000999450684;
          break;
        case 15:
          r0ab = 3.3450000286102295;
          break;
        case 16:
          r0ab = 3.1651000976562500;
          break;
        case 17:
          r0ab = 3.1259000301361084;
          break;
        case 18:
          r0ab = 3.0146999359130859;
          break;
        case 19:
          r0ab = 3.9498000144958496;
          break;
        case 20:
          r0ab = 3.3856999874114990;
          break;
        case 21:
          r0ab = 3.2916998863220215;
          break;
        case 22:
          r0ab = 3.2153999805450439;
          break;
        case 23:
          r0ab = 3.1603999137878418;
          break;
        case 24:
          r0ab = 3.4173998832702637;
          break;
        case 25:
          r0ab = 3.0734999179840088;
          break;
        case 26:
          r0ab = 3.0341999530792236;
          break;
        case 27:
          r0ab = 3.0095999240875244;
          break;
        case 28:
          r0ab = 3.0136001110076904;
          break;
        case 29:
          r0ab = 2.9855000972747803;
          break;
        case 30:
          r0ab = 2.9679999351501465;
          break;
        case 31:
          r0ab = 3.3603999614715576;
          break;
        case 32:
          r0ab = 3.4037001132965088;
          break;
        case 33:
          r0ab = 3.4242999553680420;
          break;
        case 34:
          r0ab = 3.2632999420166016;
          break;
        case 35:
          r0ab = 3.1809999942779541;
          break;
        case 36:
          r0ab = 3.1350998878479004;
          break;
        case 37:
          r0ab = 4.0556998252868652;
          break;
        case 38:
          r0ab = 3.5367999076843262;
          break;
        case 39:
          r0ab = 3.4526000022888184;
          break;
        case 40:
          r0ab = 3.3698999881744385;
          break;
        case 41:
          r0ab = 3.5706999301910400;
          break;
        case 42:
          r0ab = 3.5183999538421631;
          break;
        case 43:
          r0ab = 3.4084999561309814;
          break;
        case 44:
          r0ab = 3.3594999313354492;
          break;
        case 45:
          r0ab = 3.3333001136779785;
          break;
        case 46:
          r0ab = 3.3143000602722168;
          break;
        case 47:
          r0ab = 3.3041000366210938;
          break;
        case 48:
          r0ab = 3.1094000339508057;
          break;
        case 49:
          r0ab = 3.5192999839782715;
          break;
        case 50:
          r0ab = 3.5745000839233398;
          break;
        case 51:
          r0ab = 3.6024999618530273;
          break;
        case 52:
          r0ab = 3.4337999820709229;
          break;
        case 53:
          r0ab = 3.3447999954223633;
          break;
        case 54:
          r0ab = 3.2952001094818115;
          break;
        case 55:
          r0ab = 4.2157998085021973;
          break;
        case 56:
          r0ab = 3.7802000045776367;
          break;
        case 57:
          r0ab = 3.6431000232696533;
          break;
        case 58:
          r0ab = 3.6129000186920166;
          break;
        case 59:
          r0ab = 3.5852999687194824;
          break;
        case 60:
          r0ab = 3.5610001087188721;
          break;
        case 61:
          r0ab = 3.5406000614166260;
          break;
        case 62:
          r0ab = 3.5204000473022461;
          break;
        case 63:
          r0ab = 3.5035998821258545;
          break;
        case 64:
          r0ab = 3.2678999900817871;
          break;
        case 65:
          r0ab = 3.2162001132965088;
          break;
        case 66:
          r0ab = 3.7067999839782715;
          break;
        case 67:
          r0ab = 3.4482998847961426;
          break;
        case 68:
          r0ab = 3.4323000907897949;
          break;
        case 69:
          r0ab = 3.4221000671386719;
          break;
        case 70:
          r0ab = 3.4138000011444092;
          break;
        case 71:
          r0ab = 3.3652000427246094;
          break;
        case 72:
          r0ab = 3.4576001167297363;
          break;
        case 73:
          r0ab = 3.4052999019622803;
          break;
        case 74:
          r0ab = 3.3617999553680420;
          break;
        case 75:
          r0ab = 3.3224000930786133;
          break;
        case 76:
          r0ab = 3.2711000442504883;
          break;
        case 77:
          r0ab = 3.3326001167297363;
          break;
        case 78:
          r0ab = 3.2950000762939453;
          break;
        case 79:
          r0ab = 3.2564001083374023;
          break;
      }
      break;
    case 80:
      switch (atmA) {
        case 1:
          r0ab = 2.5315001010894775;
          break;
        case 2:
          r0ab = 2.6103999614715576;
          break;
        case 3:
          r0ab = 3.2734000682830811;
          break;
        case 4:
          r0ab = 3.2298998832702637;
          break;
        case 5:
          r0ab = 3.1089999675750732;
          break;
        case 6:
          r0ab = 2.9941999912261963;
          break;
        case 7:
          r0ab = 2.9158999919891357;
          break;
        case 8:
          r0ab = 2.8324000835418701;
          break;
        case 9:
          r0ab = 2.8350000381469727;
          break;
        case 10:
          r0ab = 2.7216000556945801;
          break;
        case 11:
          r0ab = 3.3993999958038330;
          break;
        case 12:
          r0ab = 3.4474999904632568;
          break;
        case 13:
          r0ab = 3.4354000091552734;
          break;
        case 14:
          r0ab = 3.3438000679016113;
          break;
        case 15:
          r0ab = 3.2806999683380127;
          break;
        case 16:
          r0ab = 3.2169001102447510;
          break;
        case 17:
          r0ab = 3.2676999568939209;
          break;
        case 18:
          r0ab = 3.1296000480651855;
          break;
        case 19:
          r0ab = 3.7493000030517578;
          break;
        case 20:
          r0ab = 3.8074998855590820;
          break;
        case 21:
          r0ab = 3.6846001148223877;
          break;
        case 22:
          r0ab = 3.6103999614715576;
          break;
        case 23:
          r0ab = 3.5576999187469482;
          break;
        case 24:
          r0ab = 3.2051999568939209;
          break;
        case 25:
          r0ab = 3.4802999496459961;
          break;
        case 26:
          r0ab = 3.4235999584197998;
          break;
        case 27:
          r0ab = 3.3845000267028809;
          break;
        case 28:
          r0ab = 3.3640000820159912;
          break;
        case 29:
          r0ab = 3.3364999294281006;
          break;
        case 30:
          r0ab = 3.3010001182556152;
          break;
        case 31:
          r0ab = 3.3938000202178955;
          break;
        case 32:
          r0ab = 3.3624000549316406;
          break;
        case 33:
          r0ab = 3.3440001010894775;
          break;
        case 34:
          r0ab = 3.3131999969482422;
          break;
        case 35:
          r0ab = 3.4035000801086426;
          break;
        case 36:
          r0ab = 3.2753999233245850;
          break;
        case 37:
          r0ab = 3.8701000213623047;
          break;
        case 38:
          r0ab = 3.9523000717163086;
          break;
        case 39:
          r0ab = 3.8018000125885010;
          break;
        case 40:
          r0ab = 3.7149000167846680;
          break;
        case 41:
          r0ab = 3.3673000335693359;
          break;
        case 42:
          r0ab = 3.3199000358581543;
          break;
        case 43:
          r0ab = 3.2483000755310059;
          break;
        case 44:
          r0ab = 3.2068998813629150;
          break;
        case 45:
          r0ab = 3.1793000698089600;
          break;
        case 46:
          r0ab = 3.1558001041412354;
          break;
        case 47:
          r0ab = 3.1394999027252197;
          break;
        case 48:
          r0ab = 3.4096999168395996;
          break;
        case 49:
          r0ab = 3.5409998893737793;
          break;
        case 50:
          r0ab = 3.5227999687194824;
          break;
        case 51:
          r0ab = 3.5116000175476074;
          break;
        case 52:
          r0ab = 3.4921000003814697;
          break;
        case 53:
          r0ab = 3.4781000614166260;
          break;
        case 54:
          r0ab = 3.4690001010894775;
          break;
        case 55:
          r0ab = 4.0419998168945312;
          break;
        case 56:
          r0ab = 4.1758999824523926;
          break;
        case 57:
          r0ab = 4.0078001022338867;
          break;
        case 58:
          r0ab = 4.0450000762939453;
          break;
        case 59:
          r0ab = 4.0188999176025391;
          break;
        case 60:
          r0ab = 3.9951999187469482;
          break;
        case 61:
          r0ab = 3.9769999980926514;
          break;
        case 62:
          r0ab = 3.9583001136779785;
          break;
        case 63:
          r0ab = 3.9433999061584473;
          break;
        case 64:
          r0ab = 3.7216999530792236;
          break;
        case 65:
          r0ab = 3.8227999210357666;
          break;
        case 66:
          r0ab = 3.7825999259948730;
          break;
        case 67:
          r0ab = 3.8640000820159912;
          break;
        case 68:
          r0ab = 3.8445999622344971;
          break;
        case 69:
          r0ab = 3.8313999176025391;
          break;
        case 70:
          r0ab = 3.8224999904632568;
          break;
        case 71:
          r0ab = 3.6816999912261963;
          break;
        case 72:
          r0ab = 3.7067999839782715;
          break;
        case 73:
          r0ab = 3.6554999351501465;
          break;
        case 74:
          r0ab = 3.6159000396728516;
          break;
        case 75:
          r0ab = 3.5831000804901123;
          break;
        case 76:
          r0ab = 3.5257000923156738;
          break;
        case 77:
          r0ab = 3.2132999897003174;
          break;
        case 78:
          r0ab = 3.1689000129699707;
          break;
        case 79:
          r0ab = 3.1196000576019287;
          break;
        case 80:
          r0ab = 3.3598999977111816;
          break;
      }
      break;
    case 81:
      switch (atmA) {
        case 1:
          r0ab = 2.9851999282836914;
          break;
        case 2:
          r0ab = 2.7881000041961670;
          break;
        case 3:
          r0ab = 3.5283999443054199;
          break;
        case 4:
          r0ab = 3.3492999076843262;
          break;
        case 5:
          r0ab = 3.6958000659942627;
          break;
        case 6:
          r0ab = 3.3642001152038574;
          break;
        case 7:
          r0ab = 3.1568000316619873;
          break;
        case 8:
          r0ab = 3.0055000782012939;
          break;
        case 9:
          r0ab = 2.9558000564575195;
          break;
        case 10:
          r0ab = 2.8392999172210693;
          break;
        case 11:
          r0ab = 3.6287000179290771;
          break;
        case 12:
          r0ab = 3.5283000469207764;
          break;
        case 13:
          r0ab = 4.1511001586914062;
          break;
        case 14:
          r0ab = 3.8259000778198242;
          break;
        case 15:
          r0ab = 3.6066000461578369;
          break;
        case 16:
          r0ab = 3.4526998996734619;
          break;
        case 17:
          r0ab = 3.3480000495910645;
          break;
        case 18:
          r0ab = 3.2713000774383545;
          break;
        case 19:
          r0ab = 3.9037001132965088;
          break;
        case 20:
          r0ab = 3.8361001014709473;
          break;
        case 21:
          r0ab = 3.8578999042510986;
          break;
        case 22:
          r0ab = 3.7311000823974609;
          break;
        case 23:
          r0ab = 3.6575000286102295;
          break;
        case 24:
          r0ab = 3.5176000595092773;
          break;
        case 25:
          r0ab = 3.5692999362945557;
          break;
        case 26:
          r0ab = 3.5157001018524170;
          break;
        case 27:
          r0ab = 3.4814000129699707;
          break;
        case 28:
          r0ab = 3.4558999538421631;
          break;
        case 29:
          r0ab = 3.4444999694824219;
          break;
        case 30:
          r0ab = 3.4159998893737793;
          break;
        case 31:
          r0ab = 4.1230998039245605;
          break;
        case 32:
          r0ab = 3.8543000221252441;
          break;
        case 33:
          r0ab = 3.6816000938415527;
          break;
        case 34:
          r0ab = 3.5601999759674072;
          break;
        case 35:
          r0ab = 3.4797999858856201;
          break;
        case 36:
          r0ab = 3.4207999706268311;
          break;
        case 37:
          r0ab = 4.0542001724243164;
          break;
        case 38:
          r0ab = 4.0138998031616211;
          break;
        case 39:
          r0ab = 4.0164999961853027;
          break;
        case 40:
          r0ab = 3.9412000179290771;
          break;
        case 41:
          r0ab = 3.7697999477386475;
          break;
        case 42:
          r0ab = 3.6914999485015869;
          break;
        case 43:
          r0ab = 3.6043000221252441;
          break;
        case 44:
          r0ab = 3.5638999938964844;
          break;
        case 45:
          r0ab = 3.5415999889373779;
          break;
        case 46:
          r0ab = 3.5246999263763428;
          break;
        case 47:
          r0ab = 3.5153000354766846;
          break;
        case 48:
          r0ab = 3.5653998851776123;
          break;
        case 49:
          r0ab = 4.2862000465393066;
          break;
        case 50:
          r0ab = 4.0437002182006836;
          break;
        case 51:
          r0ab = 3.8870999813079834;
          break;
        case 52:
          r0ab = 3.7741000652313232;
          break;
        case 53:
          r0ab = 3.6984999179840088;
          break;
        case 54:
          r0ab = 3.6412999629974365;
          break;
        case 55:
          r0ab = 4.2344999313354492;
          break;
        case 56:
          r0ab = 4.3663001060485840;
          break;
        case 57:
          r0ab = 4.3256998062133789;
          break;
        case 58:
          r0ab = 4.0869002342224121;
          break;
        case 59:
          r0ab = 4.0612001419067383;
          break;
        case 60:
          r0ab = 4.0363998413085938;
          break;
        case 61:
          r0ab = 4.0170001983642578;
          break;
        case 62:
          r0ab = 3.9978001117706299;
          break;
        case 63:
          r0ab = 3.9834001064300537;
          break;
        case 64:
          r0ab = 3.9137001037597656;
          break;
        case 65:
          r0ab = 3.8824999332427979;
          break;
        case 66:
          r0ab = 3.8757998943328857;
          break;
        case 67:
          r0ab = 3.9142999649047852;
          break;
        case 68:
          r0ab = 3.8975999355316162;
          break;
        case 69:
          r0ab = 3.8863999843597412;
          break;
        case 70:
          r0ab = 3.8768000602722168;
          break;
        case 71:
          r0ab = 3.9189999103546143;
          break;
        case 72:
          r0ab = 4.1613001823425293;
          break;
        case 73:
          r0ab = 4.0566000938415527;
          break;
        case 74:
          r0ab = 3.9783999919891357;
          break;
        case 75:
          r0ab = 3.9116001129150391;
          break;
        case 76:
          r0ab = 3.8326001167297363;
          break;
        case 77:
          r0ab = 3.7121999263763428;
          break;
        case 78:
          r0ab = 3.6377999782562256;
          break;
        case 79:
          r0ab = 3.5576000213623047;
          break;
        case 80:
          r0ab = 3.5457000732421875;
          break;
        case 81:
          r0ab = 4.3126997947692871;
          break;
      }
      break;
    case 82:
      switch (atmA) {
        case 1:
          r0ab = 3.1159999370574951;
          break;
        case 2:
          r0ab = 2.8482000827789307;
          break;
        case 3:
          r0ab = 4.0739002227783203;
          break;
        case 4:
          r0ab = 3.3598999977111816;
          break;
        case 5:
          r0ab = 3.5697999000549316;
          break;
        case 6:
          r0ab = 3.5366001129150391;
          break;
        case 7:
          r0ab = 3.2853999137878418;
          break;
        case 8:
          r0ab = 3.1038999557495117;
          break;
        case 9:
          r0ab = 2.9953000545501709;
          break;
        case 10:
          r0ab = 2.9191999435424805;
          break;
        case 11:
          r0ab = 4.1431999206542969;
          break;
        case 12:
          r0ab = 3.5320000648498535;
          break;
        case 13:
          r0ab = 3.9477999210357666;
          break;
        case 14:
          r0ab = 4.0230998992919922;
          break;
        case 15:
          r0ab = 3.7509000301361084;
          break;
        case 16:
          r0ab = 3.5604000091552734;
          break;
        case 17:
          r0ab = 3.4340000152587891;
          break;
        case 18:
          r0ab = 3.3426001071929932;
          break;
        case 19:
          r0ab = 4.3327999114990234;
          break;
        case 20:
          r0ab = 3.8287999629974365;
          break;
        case 21:
          r0ab = 3.7822000980377197;
          break;
        case 22:
          r0ab = 3.7908999919891357;
          break;
        case 23:
          r0ab = 3.6907000541687012;
          break;
        case 24:
          r0ab = 3.6863999366760254;
          break;
        case 25:
          r0ab = 3.5792999267578125;
          break;
        case 26:
          r0ab = 3.5220999717712402;
          break;
        case 27:
          r0ab = 3.4883000850677490;
          break;
        case 28:
          r0ab = 3.4649000167846680;
          break;
        case 29:
          r0ab = 3.4514000415802002;
          break;
        case 30:
          r0ab = 3.4300999641418457;
          break;
        case 31:
          r0ab = 3.9256000518798828;
          break;
        case 32:
          r0ab = 4.0595998764038086;
          break;
        case 33:
          r0ab = 3.8306999206542969;
          break;
        case 34:
          r0ab = 3.6702001094818115;
          break;
        case 35:
          r0ab = 3.5650999546051025;
          break;
        case 36:
          r0ab = 3.4883999824523926;
          break;
        case 37:
          r0ab = 4.4182000160217285;
          break;
        case 38:
          r0ab = 4.2515997886657715;
          break;
        case 39:
          r0ab = 3.9686999320983887;
          break;
        case 40:
          r0ab = 3.9186000823974609;
          break;
        case 41:
          r0ab = 3.9484999179840088;
          break;
        case 42:
          r0ab = 3.8369998931884766;
          break;
        case 43:
          r0ab = 3.7255001068115234;
          break;
        case 44:
          r0ab = 3.6744000911712646;
          break;
        case 45:
          r0ab = 3.6475999355316162;
          break;
        case 46:
          r0ab = 3.6294999122619629;
          break;
        case 47:
          r0ab = 3.6192998886108398;
          break;
        case 48:
          r0ab = 3.5659000873565674;
          break;
        case 49:
          r0ab = 4.0662999153137207;
          break;
        case 50:
          r0ab = 4.2308998107910156;
          break;
        case 51:
          r0ab = 4.0183000564575195;
          break;
        case 52:
          r0ab = 3.8680000305175781;
          break;
        case 53:
          r0ab = 3.7671999931335449;
          break;
        case 54:
          r0ab = 3.6923000812530518;
          break;
        case 55:
          r0ab = 4.5240001678466797;
          break;
        case 56:
          r0ab = 4.4833998680114746;
          break;
        case 57:
          r0ab = 4.1570000648498535;
          break;
        case 58:
          r0ab = 4.3204002380371094;
          break;
        case 59:
          r0ab = 4.2993001937866211;
          break;
        case 60:
          r0ab = 4.2803997993469238;
          break;
        case 61:
          r0ab = 4.2646999359130859;
          break;
        case 62:
          r0ab = 4.2480998039245605;
          break;
        case 63:
          r0ab = 4.2354001998901367;
          break;
        case 64:
          r0ab = 3.8626000881195068;
          break;
        case 65:
          r0ab = 3.8447999954223633;
          break;
        case 66:
          r0ab = 4.2266998291015625;
          break;
        case 67:
          r0ab = 4.1799001693725586;
          break;
        case 68:
          r0ab = 4.1669998168945312;
          break;
        case 69:
          r0ab = 3.8738000392913818;
          break;
        case 70:
          r0ab = 3.8643000125885010;
          break;
        case 71:
          r0ab = 3.8796000480651855;
          break;
        case 72:
          r0ab = 4.0574998855590820;
          break;
        case 73:
          r0ab = 4.0353999137878418;
          break;
        case 74:
          r0ab = 3.9365000724792480;
          break;
        case 75:
          r0ab = 3.8610999584197998;
          break;
        case 76:
          r0ab = 3.7846999168395996;
          break;
        case 77:
          r0ab = 3.7388000488281250;
          break;
        case 78:
          r0ab = 3.6826000213623047;
          break;
        case 79:
          r0ab = 3.6250998973846436;
          break;
        case 80:
          r0ab = 3.5492000579833984;
          break;
        case 81:
          r0ab = 4.0889000892639160;
          break;
        case 82:
          r0ab = 4.2764000892639160;
          break;
      }
      break;
    case 83:
      switch (atmA) {
        case 1:
          r0ab = 3.1415998935699463;
          break;
        case 2:
          r0ab = 2.8324999809265137;
          break;
        case 3:
          r0ab = 3.7734999656677246;
          break;
        case 4:
          r0ab = 3.3787000179290771;
          break;
        case 5:
          r0ab = 3.4632000923156738;
          break;
        case 6:
          r0ab = 3.5922999382019043;
          break;
        case 7:
          r0ab = 3.3213999271392822;
          break;
        case 8:
          r0ab = 3.1284999847412109;
          break;
        case 9:
          r0ab = 3.0146999359130859;
          break;
        case 10:
          r0ab = 2.9365999698638916;
          break;
        case 11:
          r0ab = 3.8526999950408936;
          break;
        case 12:
          r0ab = 3.5601999759674072;
          break;
        case 13:
          r0ab = 3.8131000995635986;
          break;
        case 14:
          r0ab = 3.8348999023437500;
          break;
        case 15:
          r0ab = 3.7994999885559082;
          break;
        case 16:
          r0ab = 3.5919001102447510;
          break;
        case 17:
          r0ab = 3.4539000988006592;
          break;
        case 18:
          r0ab = 3.3540000915527344;
          break;
        case 19:
          r0ab = 4.0654001235961914;
          break;
        case 20:
          r0ab = 3.8603000640869141;
          break;
        case 21:
          r0ab = 3.7971999645233154;
          break;
        case 22:
          r0ab = 3.7358000278472900;
          break;
        case 23:
          r0ab = 3.7392001152038574;
          break;
        case 24:
          r0ab = 3.8157000541687012;
          break;
        case 25:
          r0ab = 3.6054999828338623;
          break;
        case 26:
          r0ab = 3.5438001155853271;
          break;
        case 27:
          r0ab = 3.5088999271392822;
          break;
        case 28:
          r0ab = 3.4853000640869141;
          break;
        case 29:
          r0ab = 3.4697999954223633;
          break;
        case 30:
          r0ab = 3.4507999420166016;
          break;
        case 31:
          r0ab = 3.7881999015808105;
          break;
        case 32:
          r0ab = 3.8682000637054443;
          break;
        case 33:
          r0ab = 3.8836998939514160;
          break;
        case 34:
          r0ab = 3.7054998874664307;
          break;
        case 35:
          r0ab = 3.5869998931884766;
          break;
        case 36:
          r0ab = 3.5000000000000000;
          break;
        case 37:
          r0ab = 4.1572999954223633;
          break;
        case 38:
          r0ab = 4.0005002021789551;
          break;
        case 39:
          r0ab = 3.9567999839782715;
          break;
        case 40:
          r0ab = 3.8935999870300293;
          break;
        case 41:
          r0ab = 3.9990000724792480;
          break;
        case 42:
          r0ab = 3.9433000087738037;
          break;
        case 43:
          r0ab = 3.8171999454498291;
          break;
        case 44:
          r0ab = 3.7565999031066895;
          break;
        case 45:
          r0ab = 3.7246000766754150;
          break;
        case 46:
          r0ab = 3.7032999992370605;
          break;
        case 47:
          r0ab = 3.6900000572204590;
          break;
        case 48:
          r0ab = 3.5697000026702881;
          break;
        case 49:
          r0ab = 3.9182999134063721;
          break;
        case 50:
          r0ab = 4.0261998176574707;
          break;
        case 51:
          r0ab = 4.0658998489379883;
          break;
        case 52:
          r0ab = 3.8968999385833740;
          break;
        case 53:
          r0ab = 3.7809000015258789;
          break;
        case 54:
          r0ab = 3.6949000358581543;
          break;
        case 55:
          r0ab = 4.2765002250671387;
          break;
        case 56:
          r0ab = 4.2312002182006836;
          break;
        case 57:
          r0ab = 4.1401000022888184;
          break;
        case 58:
          r0ab = 4.0815000534057617;
          break;
        case 59:
          r0ab = 4.0580000877380371;
          break;
        case 60:
          r0ab = 4.0369000434875488;
          break;
        case 61:
          r0ab = 4.0194001197814941;
          break;
        case 62:
          r0ab = 4.0016999244689941;
          break;
        case 63:
          r0ab = 3.9874000549316406;
          break;
        case 64:
          r0ab = 3.8311998844146729;
          break;
        case 65:
          r0ab = 3.8120000362396240;
          break;
        case 66:
          r0ab = 3.9453999996185303;
          break;
        case 67:
          r0ab = 3.9210000038146973;
          break;
        case 68:
          r0ab = 3.9054999351501465;
          break;
        case 69:
          r0ab = 3.8951001167297363;
          break;
        case 70:
          r0ab = 3.8866000175476074;
          break;
        case 71:
          r0ab = 3.8689000606536865;
          break;
        case 72:
          r0ab = 3.9602999687194824;
          break;
        case 73:
          r0ab = 3.9109001159667969;
          break;
        case 74:
          r0ab = 3.9121999740600586;
          break;
        case 75:
          r0ab = 3.8232998847961426;
          break;
        case 76:
          r0ab = 3.7437999248504639;
          break;
        case 77:
          r0ab = 3.7435998916625977;
          break;
        case 78:
          r0ab = 3.6981000900268555;
          break;
        case 79:
          r0ab = 3.6554999351501465;
          break;
        case 80:
          r0ab = 3.5452001094818115;
          break;
        case 81:
          r0ab = 3.9326999187469482;
          break;
        case 82:
          r0ab = 4.0658001899719238;
          break;
        case 83:
          r0ab = 4.1174998283386230;
          break;
      }
      break;
    case 84:
      switch (atmA) {
        case 1:
          r0ab = 2.9663999080657959;
          break;
        case 2:
          r0ab = 2.8208999633789062;
          break;
        case 3:
          r0ab = 3.5546998977661133;
          break;
        case 4:
          r0ab = 3.3796000480651855;
          break;
        case 5:
          r0ab = 3.3984999656677246;
          break;
        case 6:
          r0ab = 3.3164000511169434;
          break;
        case 7:
          r0ab = 3.2363998889923096;
          break;
        case 8:
          r0ab = 3.1956000328063965;
          break;
        case 9:
          r0ab = 3.0369999408721924;
          break;
        case 10:
          r0ab = 2.9312999248504639;
          break;
        case 11:
          r0ab = 3.6424999237060547;
          break;
        case 12:
          r0ab = 3.5564999580383301;
          break;
        case 13:
          r0ab = 3.7209000587463379;
          break;
        case 14:
          r0ab = 3.7107999324798584;
          break;
        case 15:
          r0ab = 3.6638998985290527;
          break;
        case 16:
          r0ab = 3.6484000682830811;
          break;
        case 17:
          r0ab = 3.4744999408721924;
          break;
        case 18:
          r0ab = 3.3492000102996826;
          break;
        case 19:
          r0ab = 3.8754999637603760;
          break;
        case 20:
          r0ab = 4.2456998825073242;
          break;
        case 21:
          r0ab = 3.7757999897003174;
          break;
        case 22:
          r0ab = 3.7160999774932861;
          break;
        case 23:
          r0ab = 3.6693000793457031;
          break;
        case 24:
          r0ab = 3.6154999732971191;
          break;
        case 25:
          r0ab = 3.5940999984741211;
          break;
        case 26:
          r0ab = 3.5643000602722168;
          break;
        case 27:
          r0ab = 3.5292000770568848;
          break;
        case 28:
          r0ab = 3.4949998855590820;
          break;
        case 29:
          r0ab = 3.4719998836517334;
          break;
        case 30:
          r0ab = 3.4502999782562256;
          break;
        case 31:
          r0ab = 3.6935999393463135;
          break;
        case 32:
          r0ab = 3.7392001152038574;
          break;
        case 33:
          r0ab = 3.7388000488281250;
          break;
        case 34:
          r0ab = 3.7602000236511230;
          break;
        case 35:
          r0ab = 3.6078000068664551;
          break;
        case 36:
          r0ab = 3.4960000514984131;
          break;
        case 37:
          r0ab = 3.9800000190734863;
          break;
        case 38:
          r0ab = 4.3517999649047852;
          break;
        case 39:
          r0ab = 4.2802000045776367;
          break;
        case 40:
          r0ab = 3.8580000400543213;
          break;
        case 41:
          r0ab = 3.8055999279022217;
          break;
        case 42:
          r0ab = 3.7527000904083252;
          break;
        case 43:
          r0ab = 3.7019000053405762;
          break;
        case 44:
          r0ab = 3.6614999771118164;
          break;
        case 45:
          r0ab = 3.5768001079559326;
          break;
        case 46:
          r0ab = 3.5329999923706055;
          break;
        case 47:
          r0ab = 3.5037999153137207;
          break;
        case 48:
          r0ab = 3.5638999938964844;
          break;
        case 49:
          r0ab = 3.8192000389099121;
          break;
        case 50:
          r0ab = 3.8882999420166016;
          break;
        case 51:
          r0ab = 3.9091999530792236;
          break;
        case 52:
          r0ab = 3.9477999210357666;
          break;
        case 53:
          r0ab = 3.7994999885559082;
          break;
        case 54:
          r0ab = 3.6895999908447266;
          break;
        case 55:
          r0ab = 4.1164999008178711;
          break;
        case 56:
          r0ab = 4.5232000350952148;
          break;
        case 57:
          r0ab = 4.4356999397277832;
          break;
        case 58:
          r0ab = 4.4225997924804688;
          break;
        case 59:
          r0ab = 4.4031000137329102;
          break;
        case 60:
          r0ab = 4.3860001564025879;
          break;
        case 61:
          r0ab = 4.3720998764038086;
          break;
        case 62:
          r0ab = 4.3579998016357422;
          break;
        case 63:
          r0ab = 4.3466000556945801;
          break;
        case 64:
          r0ab = 4.2035999298095703;
          break;
        case 65:
          r0ab = 4.2037000656127930;
          break;
        case 66:
          r0ab = 3.8866999149322510;
          break;
        case 67:
          r0ab = 4.2895002365112305;
          break;
        case 68:
          r0ab = 4.2765998840332031;
          break;
        case 69:
          r0ab = 4.2662000656127930;
          break;
        case 70:
          r0ab = 4.2597999572753906;
          break;
        case 71:
          r0ab = 3.8408000469207764;
          break;
        case 72:
          r0ab = 3.9168999195098877;
          break;
        case 73:
          r0ab = 3.8680999279022217;
          break;
        case 74:
          r0ab = 3.8250000476837158;
          break;
        case 75:
          r0ab = 3.7855000495910645;
          break;
        case 76:
          r0ab = 3.7500998973846436;
          break;
        case 77:
          r0ab = 3.6752998828887939;
          break;
        case 78:
          r0ab = 3.5499000549316406;
          break;
        case 79:
          r0ab = 3.4872000217437744;
          break;
        case 80:
          r0ab = 3.5401000976562500;
          break;
        case 81:
          r0ab = 3.8287999629974365;
          break;
        case 82:
          r0ab = 3.9217000007629395;
          break;
        case 83:
          r0ab = 3.9537999629974365;
          break;
        case 84:
          r0ab = 4.0054001808166504;
          break;
      }
      break;
    case 85:
      switch (atmA) {
        case 1:
          r0ab = 2.8387999534606934;
          break;
        case 2:
          r0ab = 2.7890000343322754;
          break;
        case 3:
          r0ab = 3.4328999519348145;
          break;
        case 4:
          r0ab = 3.5592999458312988;
          break;
        case 5:
          r0ab = 3.3487999439239502;
          break;
        case 6:
          r0ab = 3.2486000061035156;
          break;
        case 7:
          r0ab = 3.1614999771118164;
          break;
        case 8:
          r0ab = 3.0999999046325684;
          break;
        case 9:
          r0ab = 3.0394001007080078;
          break;
        case 10:
          r0ab = 2.9165000915527344;
          break;
        case 11:
          r0ab = 3.5267000198364258;
          break;
        case 12:
          r0ab = 3.7479000091552734;
          break;
        case 13:
          r0ab = 3.6649999618530273;
          break;
        case 14:
          r0ab = 3.6263000965118408;
          break;
        case 15:
          r0ab = 3.5657999515533447;
          break;
        case 16:
          r0ab = 3.5223999023437500;
          break;
        case 17:
          r0ab = 3.4762001037597656;
          break;
        case 18:
          r0ab = 3.3341999053955078;
          break;
        case 19:
          r0ab = 3.7737998962402344;
          break;
        case 20:
          r0ab = 4.0332999229431152;
          break;
        case 21:
          r0ab = 3.9567999839782715;
          break;
        case 22:
          r0ab = 3.8975000381469727;
          break;
        case 23:
          r0ab = 3.8520998954772949;
          break;
        case 24:
          r0ab = 3.4928998947143555;
          break;
        case 25:
          r0ab = 3.7829999923706055;
          break;
        case 26:
          r0ab = 3.7409000396728516;
          break;
        case 27:
          r0ab = 3.7061998844146729;
          break;
        case 28:
          r0ab = 3.6786000728607178;
          break;
        case 29:
          r0ab = 3.6470999717712402;
          break;
        case 30:
          r0ab = 3.6208000183105469;
          break;
        case 31:
          r0ab = 3.6336998939514160;
          break;
        case 32:
          r0ab = 3.6519000530242920;
          break;
        case 33:
          r0ab = 3.6363000869750977;
          break;
        case 34:
          r0ab = 3.6277999877929688;
          break;
        case 35:
          r0ab = 3.6110000610351562;
          break;
        case 36:
          r0ab = 3.4825000762939453;
          break;
        case 37:
          r0ab = 3.8794999122619629;
          break;
        case 38:
          r0ab = 4.1448001861572266;
          break;
        case 39:
          r0ab = 4.0735998153686523;
          break;
        case 40:
          r0ab = 4.0044999122619629;
          break;
        case 41:
          r0ab = 3.6842999458312988;
          break;
        case 42:
          r0ab = 3.6291000843048096;
          break;
        case 43:
          r0ab = 3.5741000175476074;
          break;
        case 44:
          r0ab = 3.5311999320983887;
          break;
        case 45:
          r0ab = 3.4974000453948975;
          break;
        case 46:
          r0ab = 3.4472000598907471;
          break;
        case 47:
          r0ab = 3.4033999443054199;
          break;
        case 48:
          r0ab = 3.7130999565124512;
          break;
        case 49:
          r0ab = 3.7557001113891602;
          break;
        case 50:
          r0ab = 3.7966001033782959;
          break;
        case 51:
          r0ab = 3.8004999160766602;
          break;
        case 52:
          r0ab = 3.8067998886108398;
          break;
        case 53:
          r0ab = 3.8015000820159912;
          break;
        case 54:
          r0ab = 3.6747000217437744;
          break;
        case 55:
          r0ab = 4.0222001075744629;
          break;
        case 56:
          r0ab = 4.3207001686096191;
          break;
        case 57:
          r0ab = 4.2347002029418945;
          break;
        case 58:
          r0ab = 4.2190999984741211;
          break;
        case 59:
          r0ab = 4.1989998817443848;
          break;
        case 60:
          r0ab = 4.1810998916625977;
          break;
        case 61:
          r0ab = 4.1666002273559570;
          break;
        case 62:
          r0ab = 4.1521000862121582;
          break;
        case 63:
          r0ab = 4.1401000022888184;
          break;
        case 64:
          r0ab = 3.9969999790191650;
          break;
        case 65:
          r0ab = 3.9942998886108398;
          break;
        case 66:
          r0ab = 3.9591999053955078;
          break;
        case 67:
          r0ab = 4.0799999237060547;
          break;
        case 68:
          r0ab = 4.0664000511169434;
          break;
        case 69:
          r0ab = 4.0559000968933105;
          break;
        case 70:
          r0ab = 4.0487999916076660;
          break;
        case 71:
          r0ab = 3.9881999492645264;
          break;
        case 72:
          r0ab = 4.0034999847412109;
          break;
        case 73:
          r0ab = 3.9539000988006592;
          break;
        case 74:
          r0ab = 3.9138000011444092;
          break;
        case 75:
          r0ab = 3.8798000812530518;
          break;
        case 76:
          r0ab = 3.8355000019073486;
          break;
        case 77:
          r0ab = 3.5359001159667969;
          break;
        case 78:
          r0ab = 3.4953999519348145;
          break;
        case 79:
          r0ab = 3.3961999416351318;
          break;
        case 80:
          r0ab = 3.5339000225067139;
          break;
        case 81:
          r0ab = 3.7595000267028809;
          break;
        case 82:
          r0ab = 3.8250000476837158;
          break;
        case 83:
          r0ab = 3.8408000469207764;
          break;
        case 84:
          r0ab = 3.8599998950958252;
          break;
        case 85:
          r0ab = 3.8643999099731445;
          break;
      }
      break;
    case 86:
      switch (atmA) {
        case 1:
          r0ab = 2.7411999702453613;
          break;
        case 2:
          r0ab = 2.7488999366760254;
          break;
        case 3:
          r0ab = 3.3373999595642090;
          break;
        case 4:
          r0ab = 3.3949999809265137;
          break;
        case 5:
          r0ab = 3.3076000213623047;
          break;
        case 6:
          r0ab = 3.1909999847412109;
          break;
        case 7:
          r0ab = 3.0961000919342041;
          break;
        case 8:
          r0ab = 3.0174999237060547;
          break;
        case 9:
          r0ab = 3.0280001163482666;
          break;
        case 10:
          r0ab = 2.8928999900817871;
          break;
        case 11:
          r0ab = 3.4328000545501709;
          break;
        case 12:
          r0ab = 3.5882999897003174;
          break;
        case 13:
          r0ab = 3.6226999759674072;
          break;
        case 14:
          r0ab = 3.5615999698638916;
          break;
        case 15:
          r0ab = 3.4893999099731445;
          break;
        case 16:
          r0ab = 3.4240999221801758;
          break;
        case 17:
          r0ab = 3.3640999794006348;
          break;
        case 18:
          r0ab = 3.3120000362396240;
          break;
        case 19:
          r0ab = 3.6814999580383301;
          break;
        case 20:
          r0ab = 3.8789000511169434;
          break;
        case 21:
          r0ab = 3.8031001091003418;
          break;
        case 22:
          r0ab = 3.7413001060485840;
          break;
        case 23:
          r0ab = 3.6939001083374023;
          break;
        case 24:
          r0ab = 3.4010000228881836;
          break;
        case 25:
          r0ab = 3.6224999427795410;
          break;
        case 26:
          r0ab = 3.5796999931335449;
          break;
        case 27:
          r0ab = 3.5443000793457031;
          break;
        case 28:
          r0ab = 3.5139000415802002;
          break;
        case 29:
          r0ab = 3.4923000335693359;
          break;
        case 30:
          r0ab = 3.4642000198364258;
          break;
        case 31:
          r0ab = 3.5859999656677246;
          break;
        case 32:
          r0ab = 3.5848999023437500;
          break;
        case 33:
          r0ab = 3.5569999217987061;
          break;
        case 34:
          r0ab = 3.5257000923156738;
          break;
        case 35:
          r0ab = 3.4935998916625977;
          break;
        case 36:
          r0ab = 3.4628000259399414;
          break;
        case 37:
          r0ab = 3.7874000072479248;
          break;
        case 38:
          r0ab = 3.9916000366210938;
          break;
        case 39:
          r0ab = 3.9249000549316406;
          break;
        case 40:
          r0ab = 3.8529999256134033;
          break;
        case 41:
          r0ab = 3.5931999683380127;
          break;
        case 42:
          r0ab = 3.5355000495910645;
          break;
        case 43:
          r0ab = 3.4756999015808105;
          break;
        case 44:
          r0ab = 3.4305999279022217;
          break;
        case 45:
          r0ab = 3.3952999114990234;
          break;
        case 46:
          r0ab = 3.3645999431610107;
          break;
        case 47:
          r0ab = 3.3389999866485596;
          break;
        case 48:
          r0ab = 3.5636999607086182;
          break;
        case 49:
          r0ab = 3.7053000926971436;
          break;
        case 50:
          r0ab = 3.7265999317169189;
          break;
        case 51:
          r0ab = 3.7177000045776367;
          break;
        case 52:
          r0ab = 3.6995999813079834;
          break;
        case 53:
          r0ab = 3.6775000095367432;
          break;
        case 54:
          r0ab = 3.6558001041412354;
          break;
        case 55:
          r0ab = 3.9330999851226807;
          break;
        case 56:
          r0ab = 4.1655001640319824;
          break;
        case 57:
          r0ab = 4.0879001617431641;
          break;
        case 58:
          r0ab = 4.0680999755859375;
          break;
        case 59:
          r0ab = 4.0479001998901367;
          break;
        case 60:
          r0ab = 4.0299000740051270;
          break;
        case 61:
          r0ab = 4.0152001380920410;
          break;
        case 62:
          r0ab = 4.0005998611450195;
          break;
        case 63:
          r0ab = 3.9883000850677490;
          break;
        case 64:
          r0ab = 3.8499999046325684;
          break;
        case 65:
          r0ab = 3.8359000682830811;
          break;
        case 66:
          r0ab = 3.8248999118804932;
          break;
        case 67:
          r0ab = 3.9268999099731445;
          break;
        case 68:
          r0ab = 3.9133000373840332;
          break;
        case 69:
          r0ab = 3.9024999141693115;
          break;
        case 70:
          r0ab = 3.8947999477386475;
          break;
        case 71:
          r0ab = 3.8422000408172607;
          break;
        case 72:
          r0ab = 3.8508999347686768;
          break;
        case 73:
          r0ab = 3.7990000247955322;
          break;
        case 74:
          r0ab = 3.7569999694824219;
          break;
        case 75:
          r0ab = 3.7218999862670898;
          break;
        case 76:
          r0ab = 3.6761999130249023;
          break;
        case 77:
          r0ab = 3.4260001182556152;
          break;
        case 78:
          r0ab = 3.3866000175476074;
          break;
        case 79:
          r0ab = 3.3424999713897705;
          break;
        case 80:
          r0ab = 3.5294001102447510;
          break;
        case 81:
          r0ab = 3.7021999359130859;
          break;
        case 82:
          r0ab = 3.7497000694274902;
          break;
        case 83:
          r0ab = 3.7541999816894531;
          break;
        case 84:
          r0ab = 3.7493999004364014;
          break;
        case 85:
          r0ab = 3.7369999885559082;
          break;
        case 86:
          r0ab = 3.7216000556945801;
          break;
      }
      break;
    case 87:
      switch (atmA) {
        case 1:
          r0ab = 3.4154999256134033;
          break;
        case 2:
          r0ab = 3.0522000789642334;
          break;
        case 3:
          r0ab = 4.2540998458862305;
          break;
        case 4:
          r0ab = 3.8217999935150146;
          break;
        case 5:
          r0ab = 4.0437998771667480;
          break;
        case 6:
          r0ab = 3.5875000953674316;
          break;
        case 7:
          r0ab = 3.3285999298095703;
          break;
        case 8:
          r0ab = 3.1682000160217285;
          break;
        case 9:
          r0ab = 3.0566000938415527;
          break;
        case 10:
          r0ab = 2.9746000766754150;
          break;
        case 11:
          r0ab = 4.3626999855041504;
          break;
        case 12:
          r0ab = 4.0248999595642090;
          break;
        case 13:
          r0ab = 4.6946997642517090;
          break;
        case 14:
          r0ab = 4.1718001365661621;
          break;
        case 15:
          r0ab = 3.8638999462127686;
          break;
        case 16:
          r0ab = 3.6735000610351562;
          break;
        case 17:
          r0ab = 3.5434999465942383;
          break;
        case 18:
          r0ab = 3.4479000568389893;
          break;
        case 19:
          r0ab = 4.6806001663208008;
          break;
        case 20:
          r0ab = 4.3484997749328613;
          break;
        case 21:
          r0ab = 4.2667999267578125;
          break;
        case 22:
          r0ab = 4.1690001487731934;
          break;
        case 23:
          r0ab = 4.1061000823974609;
          break;
        case 24:
          r0ab = 4.1244997978210449;
          break;
        case 25:
          r0ab = 4.0205998420715332;
          break;
        case 26:
          r0ab = 3.9765000343322754;
          break;
        case 27:
          r0ab = 3.9458000659942627;
          break;
        case 28:
          r0ab = 3.9217000007629395;
          break;
        case 29:
          r0ab = 3.9075000286102295;
          break;
        case 30:
          r0ab = 3.8812999725341797;
          break;
        case 31:
          r0ab = 3.9946999549865723;
          break;
        case 32:
          r0ab = 4.1989002227783203;
          break;
        case 33:
          r0ab = 3.9507000446319580;
          break;
        case 34:
          r0ab = 3.7960000038146973;
          break;
        case 35:
          r0ab = 3.6925001144409180;
          break;
        case 36:
          r0ab = 3.6150000095367432;
          break;
        case 37:
          r0ab = 4.8534998893737793;
          break;
        case 38:
          r0ab = 4.5641999244689941;
          break;
        case 39:
          r0ab = 4.4134001731872559;
          break;
        case 40:
          r0ab = 4.3688001632690430;
          break;
        case 41:
          r0ab = 4.3396000862121582;
          break;
        case 42:
          r0ab = 4.2878999710083008;
          break;
        case 43:
          r0ab = 4.2165999412536621;
          break;
        case 44:
          r0ab = 4.1887998580932617;
          break;
        case 45:
          r0ab = 4.1767997741699219;
          break;
        case 46:
          r0ab = 4.1659998893737793;
          break;
        case 47:
          r0ab = 4.1607999801635742;
          break;
        case 48:
          r0ab = 4.0745000839233398;
          break;
        case 49:
          r0ab = 4.2288999557495117;
          break;
        case 50:
          r0ab = 4.4862999916076660;
          break;
        case 51:
          r0ab = 4.2512998580932617;
          break;
        case 52:
          r0ab = 4.0897002220153809;
          break;
        case 53:
          r0ab = 3.9876000881195068;
          break;
        case 54:
          r0ab = 3.9061000347137451;
          break;
        case 55:
          r0ab = 5.0689997673034668;
          break;
        case 56:
          r0ab = 5.0446000099182129;
          break;
        case 57:
          r0ab = 4.6185998916625977;
          break;
        case 58:
          r0ab = 4.6078000068664551;
          break;
        case 59:
          r0ab = 4.5780000686645508;
          break;
        case 60:
          r0ab = 4.5538001060485840;
          break;
        case 61:
          r0ab = 4.5318999290466309;
          break;
        case 62:
          r0ab = 4.5100998878479004;
          break;
        case 63:
          r0ab = 4.4945001602172852;
          break;
        case 64:
          r0ab = 4.1911997795104980;
          break;
        case 65:
          r0ab = 4.2315001487731934;
          break;
        case 66:
          r0ab = 4.5534000396728516;
          break;
        case 67:
          r0ab = 4.4373002052307129;
          break;
        case 68:
          r0ab = 4.4223999977111816;
          break;
        case 69:
          r0ab = 4.4120001792907715;
          break;
        case 70:
          r0ab = 4.4039998054504395;
          break;
        case 71:
          r0ab = 4.2634000778198242;
          break;
        case 72:
          r0ab = 4.7769999504089355;
          break;
        case 73:
          r0ab = 4.6890001296997070;
          break;
        case 74:
          r0ab = 4.6107001304626465;
          break;
        case 75:
          r0ab = 4.5331001281738281;
          break;
        case 76:
          r0ab = 4.4496002197265625;
          break;
        case 77:
          r0ab = 4.4081997871398926;
          break;
        case 78:
          r0ab = 4.3095002174377441;
          break;
        case 79:
          r0ab = 4.2023000717163086;
          break;
        case 80:
          r0ab = 4.0500998497009277;
          break;
        case 81:
          r0ab = 4.2595000267028809;
          break;
        case 82:
          r0ab = 4.5496997833251953;
          break;
        case 83:
          r0ab = 4.3056001663208008;
          break;
        case 84:
          r0ab = 4.1505999565124512;
          break;
        case 85:
          r0ab = 4.0574002265930176;
          break;
        case 86:
          r0ab = 3.9725000858306885;
          break;
        case 87:
          r0ab = 5.0795998573303223;
          break;
      }
      break;
    case 88:
      switch (atmA) {
        case 1:
          r0ab = 3.0548000335693359;
          break;
        case 2:
          r0ab = 3.3206000328063965;
          break;
        case 3:
          r0ab = 3.8131999969482422;
          break;
        case 4:
          r0ab = 3.9719998836517334;
          break;
        case 5:
          r0ab = 3.7674999237060547;
          break;
        case 6:
          r0ab = 3.7351000308990479;
          break;
        case 7:
          r0ab = 3.5167000293731689;
          break;
        case 8:
          r0ab = 3.5274000167846680;
          break;
        case 9:
          r0ab = 3.3085000514984131;
          break;
        case 10:
          r0ab = 3.1652998924255371;
          break;
        case 11:
          r0ab = 3.9500000476837158;
          break;
        case 12:
          r0ab = 4.1729998588562012;
          break;
        case 13:
          r0ab = 4.0612998008728027;
          break;
        case 14:
          r0ab = 4.1493000984191895;
          break;
        case 15:
          r0ab = 3.8822999000549316;
          break;
        case 16:
          r0ab = 4.0536999702453613;
          break;
        case 17:
          r0ab = 3.8199999332427979;
          break;
        case 18:
          r0ab = 3.6582000255584717;
          break;
        case 19:
          r0ab = 4.3421998023986816;
          break;
        case 20:
          r0ab = 4.5110998153686523;
          break;
        case 21:
          r0ab = 4.3794999122619629;
          break;
        case 22:
          r0ab = 4.3362002372741699;
          break;
        case 23:
          r0ab = 4.2751002311706543;
          break;
        case 24:
          r0ab = 3.7102999687194824;
          break;
        case 25:
          r0ab = 4.1972999572753906;
          break;
        case 26:
          r0ab = 4.1385002136230469;
          break;
        case 27:
          r0ab = 4.1128997802734375;
          break;
        case 28:
          r0ab = 4.0799999237060547;
          break;
        case 29:
          r0ab = 4.0647001266479492;
          break;
        case 30:
          r0ab = 4.0307998657226562;
          break;
        case 31:
          r0ab = 4.0096001625061035;
          break;
        case 32:
          r0ab = 4.1619000434875488;
          break;
        case 33:
          r0ab = 3.9360001087188721;
          break;
        case 34:
          r0ab = 4.1765999794006348;
          break;
        case 35:
          r0ab = 3.9704999923706055;
          break;
        case 36:
          r0ab = 3.8262000083923340;
          break;
        case 37:
          r0ab = 4.5348000526428223;
          break;
        case 38:
          r0ab = 4.7024998664855957;
          break;
        case 39:
          r0ab = 4.5268001556396484;
          break;
        case 40:
          r0ab = 4.5075998306274414;
          break;
        case 41:
          r0ab = 3.9561998844146729;
          break;
        case 42:
          r0ab = 3.9065001010894775;
          break;
        case 43:
          r0ab = 3.8118999004364014;
          break;
        case 44:
          r0ab = 3.7604999542236328;
          break;
        case 45:
          r0ab = 3.7446999549865723;
          break;
        case 46:
          r0ab = 3.7118999958038330;
          break;
        case 47:
          r0ab = 3.6916000843048096;
          break;
        case 48:
          r0ab = 4.1950001716613770;
          break;
        case 49:
          r0ab = 4.2109999656677246;
          break;
        case 50:
          r0ab = 4.3843002319335938;
          break;
        case 51:
          r0ab = 4.1630997657775879;
          break;
        case 52:
          r0ab = 4.4426999092102051;
          break;
        case 53:
          r0ab = 4.2463002204895020;
          break;
        case 54:
          r0ab = 4.1054000854492188;
          break;
        case 55:
          r0ab = 4.7692999839782715;
          break;
        case 56:
          r0ab = 5.0648999214172363;
          break;
        case 57:
          r0ab = 4.7364997863769531;
          break;
        case 58:
          r0ab = 4.7761001586914062;
          break;
        case 59:
          r0ab = 4.7498002052307129;
          break;
        case 60:
          r0ab = 4.7272000312805176;
          break;
        case 61:
          r0ab = 4.7076001167297363;
          break;
        case 62:
          r0ab = 4.6876997947692871;
          break;
        case 63:
          r0ab = 4.6729998588562012;
          break;
        case 64:
          r0ab = 4.4274001121520996;
          break;
        case 65:
          r0ab = 4.5472998619079590;
          break;
        case 66:
          r0ab = 4.5169000625610352;
          break;
        case 67:
          r0ab = 4.5974998474121094;
          break;
        case 68:
          r0ab = 4.5792999267578125;
          break;
        case 69:
          r0ab = 4.5666999816894531;
          break;
        case 70:
          r0ab = 4.5559000968933105;
          break;
        case 71:
          r0ab = 4.3804001808166504;
          break;
        case 72:
          r0ab = 4.6919999122619629;
          break;
        case 73:
          r0ab = 4.6730999946594238;
          break;
        case 74:
          r0ab = 4.6142001152038574;
          break;
        case 75:
          r0ab = 4.5599999427795410;
          break;
        case 76:
          r0ab = 4.4801001548767090;
          break;
        case 77:
          r0ab = 4.0149002075195312;
          break;
        case 78:
          r0ab = 3.8856000900268555;
          break;
        case 79:
          r0ab = 3.7407000064849854;
          break;
        case 80:
          r0ab = 4.1545000076293945;
          break;
        case 81:
          r0ab = 4.2252998352050781;
          break;
        case 82:
          r0ab = 4.4229001998901367;
          break;
        case 83:
          r0ab = 4.1922998428344727;
          break;
        case 84:
          r0ab = 4.5022001266479492;
          break;
        case 85:
          r0ab = 4.3059000968933105;
          break;
        case 86:
          r0ab = 4.1591000556945801;
          break;
        case 87:
          r0ab = 4.7883000373840332;
          break;
        case 88:
          r0ab = 4.9293999671936035;
          break;
      }
      break;
    case 89:
      switch (atmA) {
        case 1:
          r0ab = 3.3849999904632568;
          break;
        case 2:
          r0ab = 3.4207999706268311;
          break;
        case 3:
          r0ab = 3.7004001140594482;
          break;
        case 4:
          r0ab = 3.8800001144409180;
          break;
        case 5:
          r0ab = 3.9886000156402588;
          break;
        case 6:
          r0ab = 3.9040000438690186;
          break;
        case 7:
          r0ab = 3.6719000339508057;
          break;
        case 8:
          r0ab = 3.6547000408172607;
          break;
        case 9:
          r0ab = 3.4625000953674316;
          break;
        case 10:
          r0ab = 3.3369998931884766;
          break;
        case 11:
          r0ab = 3.8394000530242920;
          break;
        case 12:
          r0ab = 4.0335001945495605;
          break;
        case 13:
          r0ab = 4.2372999191284180;
          break;
        case 14:
          r0ab = 4.3022999763488770;
          break;
        case 15:
          r0ab = 4.0306000709533691;
          break;
        case 16:
          r0ab = 4.1407999992370605;
          break;
        case 17:
          r0ab = 3.9296998977661133;
          break;
        case 18:
          r0ab = 3.7857000827789307;
          break;
        case 19:
          r0ab = 4.1907000541687012;
          break;
        case 20:
          r0ab = 4.3229999542236328;
          break;
        case 21:
          r0ab = 4.2663998603820801;
          break;
        case 22:
          r0ab = 4.2172999382019043;
          break;
        case 23:
          r0ab = 4.1482000350952148;
          break;
        case 24:
          r0ab = 3.6823000907897949;
          break;
        case 25:
          r0ab = 4.0711002349853516;
          break;
        case 26:
          r0ab = 4.0180001258850098;
          break;
        case 27:
          r0ab = 4.0016999244689941;
          break;
        case 28:
          r0ab = 3.9746999740600586;
          break;
        case 29:
          r0ab = 3.9633998870849609;
          break;
        case 30:
          r0ab = 3.9382998943328857;
          break;
        case 31:
          r0ab = 4.1992998123168945;
          break;
        case 32:
          r0ab = 4.3204998970031738;
          break;
        case 33:
          r0ab = 4.0820999145507812;
          break;
        case 34:
          r0ab = 4.2547001838684082;
          break;
        case 35:
          r0ab = 4.0658998489379883;
          break;
        case 36:
          r0ab = 3.9358999729156494;
          break;
        case 37:
          r0ab = 4.3951997756958008;
          break;
        case 38:
          r0ab = 4.5176000595092773;
          break;
        case 39:
          r0ab = 4.3888001441955566;
          break;
        case 40:
          r0ab = 4.3607001304626465;
          break;
        case 41:
          r0ab = 3.9583001136779785;
          break;
        case 42:
          r0ab = 3.9279999732971191;
          break;
        case 43:
          r0ab = 3.8389999866485596;
          break;
        case 44:
          r0ab = 3.7971000671386719;
          break;
        case 45:
          r0ab = 3.7955000400543213;
          break;
        case 46:
          r0ab = 3.7674000263214111;
          break;
        case 47:
          r0ab = 3.7520999908447266;
          break;
        case 48:
          r0ab = 4.1062002182006836;
          break;
        case 49:
          r0ab = 4.3632998466491699;
          break;
        case 50:
          r0ab = 4.2990999221801758;
          break;
        case 51:
          r0ab = 4.2767000198364258;
          break;
        case 52:
          r0ab = 4.4857001304626465;
          break;
        case 53:
          r0ab = 4.3038997650146484;
          break;
        case 54:
          r0ab = 4.1761999130249023;
          break;
        case 55:
          r0ab = 4.6196999549865723;
          break;
        case 56:
          r0ab = 4.8653998374938965;
          break;
        case 57:
          r0ab = 4.6633000373840332;
          break;
        case 58:
          r0ab = 4.5878000259399414;
          break;
        case 59:
          r0ab = 4.5640001296997070;
          break;
        case 60:
          r0ab = 4.5422000885009766;
          break;
        case 61:
          r0ab = 4.5230998992919922;
          break;
        case 62:
          r0ab = 4.5041999816894531;
          break;
        case 63:
          r0ab = 4.4900999069213867;
          break;
        case 64:
          r0ab = 4.3281998634338379;
          break;
        case 65:
          r0ab = 4.3977999687194824;
          break;
        case 66:
          r0ab = 4.3482999801635742;
          break;
        case 67:
          r0ab = 4.4201998710632324;
          break;
        case 68:
          r0ab = 4.4039001464843750;
          break;
        case 69:
          r0ab = 4.3926000595092773;
          break;
        case 70:
          r0ab = 4.3807001113891602;
          break;
        case 71:
          r0ab = 4.2649002075195312;
          break;
        case 72:
          r0ab = 4.6135001182556152;
          break;
        case 73:
          r0ab = 4.5605001449584961;
          break;
        case 74:
          r0ab = 4.5232000350952148;
          break;
        case 75:
          r0ab = 4.4675998687744141;
          break;
        case 76:
          r0ab = 4.3948001861572266;
          break;
        case 77:
          r0ab = 4.0988998413085938;
          break;
        case 78:
          r0ab = 3.9863998889923096;
          break;
        case 79:
          r0ab = 3.8596000671386719;
          break;
        case 80:
          r0ab = 4.0942001342773438;
          break;
        case 81:
          r0ab = 4.2719998359680176;
          break;
        case 82:
          r0ab = 4.3270001411437988;
          break;
        case 83:
          r0ab = 4.3021998405456543;
          break;
        case 84:
          r0ab = 4.5409998893737793;
          break;
        case 85:
          r0ab = 4.3576002120971680;
          break;
        case 86:
          r0ab = 4.2234997749328613;
          break;
        case 87:
          r0ab = 4.6545000076293945;
          break;
        case 88:
          r0ab = 4.7446999549865723;
          break;
        case 89:
          r0ab = 4.7042999267578125;
          break;
      }
      break;
    case 90:
      switch (atmA) {
        case 1:
          r0ab = 3.0941998958587646;
          break;
        case 2:
          r0ab = 3.2074999809265137;
          break;
        case 3:
          r0ab = 3.5151998996734619;
          break;
        case 4:
          r0ab = 3.6658999919891357;
          break;
        case 5:
          r0ab = 3.8289000988006592;
          break;
        case 6:
          r0ab = 3.7458999156951904;
          break;
        case 7:
          r0ab = 3.5155999660491943;
          break;
        case 8:
          r0ab = 3.5197000503540039;
          break;
        case 9:
          r0ab = 3.3289999961853027;
          break;
        case 10:
          r0ab = 3.2068998813629150;
          break;
        case 11:
          r0ab = 3.6702001094818115;
          break;
        case 12:
          r0ab = 3.8447999954223633;
          break;
        case 13:
          r0ab = 4.0339999198913574;
          break;
        case 14:
          r0ab = 3.9509000778198242;
          break;
        case 15:
          r0ab = 3.8585000038146973;
          break;
        case 16:
          r0ab = 3.9893999099731445;
          break;
        case 17:
          r0ab = 3.7787001132965088;
          break;
        case 18:
          r0ab = 3.6364998817443848;
          break;
        case 19:
          r0ab = 4.1424999237060547;
          break;
        case 20:
          r0ab = 4.1617999076843262;
          break;
        case 21:
          r0ab = 4.0939998626708984;
          break;
        case 22:
          r0ab = 4.0465998649597168;
          break;
        case 23:
          r0ab = 3.9941000938415527;
          break;
        case 24:
          r0ab = 3.5425999164581299;
          break;
        case 25:
          r0ab = 3.8952000141143799;
          break;
        case 26:
          r0ab = 3.8327000141143799;
          break;
        case 27:
          r0ab = 3.8125998973846436;
          break;
        case 28:
          r0ab = 3.7795999050140381;
          break;
        case 29:
          r0ab = 3.7634999752044678;
          break;
        case 30:
          r0ab = 3.7355999946594238;
          break;
        case 31:
          r0ab = 4.0047001838684082;
          break;
        case 32:
          r0ab = 3.9655001163482666;
          break;
        case 33:
          r0ab = 3.9116001129150391;
          break;
        case 34:
          r0ab = 4.1009998321533203;
          break;
        case 35:
          r0ab = 3.9102001190185547;
          break;
        case 36:
          r0ab = 3.7799999713897705;
          break;
        case 37:
          r0ab = 4.2964000701904297;
          break;
        case 38:
          r0ab = 4.3330001831054688;
          break;
        case 39:
          r0ab = 4.2621998786926270;
          break;
        case 40:
          r0ab = 4.2253999710083008;
          break;
        case 41:
          r0ab = 3.8194999694824219;
          break;
        case 42:
          r0ab = 3.7560000419616699;
          break;
        case 43:
          r0ab = 3.6512999534606934;
          break;
        case 44:
          r0ab = 3.5940999984741211;
          break;
        case 45:
          r0ab = 3.5810000896453857;
          break;
        case 46:
          r0ab = 3.5420000553131104;
          break;
        case 47:
          r0ab = 3.5178000926971436;
          break;
        case 48:
          r0ab = 3.8861000537872314;
          break;
        case 49:
          r0ab = 4.1458997726440430;
          break;
        case 50:
          r0ab = 4.1146998405456543;
          break;
        case 51:
          r0ab = 4.0771999359130859;
          break;
        case 52:
          r0ab = 4.3119997978210449;
          break;
        case 53:
          r0ab = 4.1206998825073242;
          break;
        case 54:
          r0ab = 3.9900000095367432;
          break;
        case 55:
          r0ab = 4.4732999801635742;
          break;
        case 56:
          r0ab = 4.6156997680664062;
          break;
        case 57:
          r0ab = 4.4580001831054688;
          break;
        case 58:
          r0ab = 4.4194002151489258;
          break;
        case 59:
          r0ab = 4.3954000473022461;
          break;
        case 60:
          r0ab = 4.3738999366760254;
          break;
        case 61:
          r0ab = 4.3530998229980469;
          break;
        case 62:
          r0ab = 4.3343000411987305;
          break;
        case 63:
          r0ab = 4.3196001052856445;
          break;
        case 64:
          r0ab = 4.2140002250671387;
          break;
        case 65:
          r0ab = 4.2339000701904297;
          break;
        case 66:
          r0ab = 4.1737999916076660;
          break;
        case 67:
          r0ab = 4.2458000183105469;
          break;
        case 68:
          r0ab = 4.2277998924255371;
          break;
        case 69:
          r0ab = 4.2157998085021973;
          break;
        case 70:
          r0ab = 4.2038998603820801;
          break;
        case 71:
          r0ab = 4.1658000946044922;
          break;
        case 72:
          r0ab = 4.3594999313354492;
          break;
        case 73:
          r0ab = 4.2856998443603516;
          break;
        case 74:
          r0ab = 4.2444000244140625;
          break;
        case 75:
          r0ab = 4.1855001449584961;
          break;
        case 76:
          r0ab = 4.1121997833251953;
          break;
        case 77:
          r0ab = 3.7839000225067139;
          break;
        case 78:
          r0ab = 3.6879000663757324;
          break;
        case 79:
          r0ab = 3.5815999507904053;
          break;
        case 80:
          r0ab = 3.8633000850677490;
          break;
        case 81:
          r0ab = 4.1585001945495605;
          break;
        case 82:
          r0ab = 4.1402001380920410;
          break;
        case 83:
          r0ab = 4.1036000251770020;
          break;
        case 84:
          r0ab = 4.3694000244140625;
          break;
        case 85:
          r0ab = 4.1735000610351562;
          break;
        case 86:
          r0ab = 4.0367999076843262;
          break;
        case 87:
          r0ab = 4.5095000267028809;
          break;
        case 88:
          r0ab = 4.5538001060485840;
          break;
        case 89:
          r0ab = 4.5240001678466797;
          break;
        case 90:
          r0ab = 4.4251999855041504;
          break;
      }
      break;
    case 91:
      switch (atmA) {
        case 1:
          r0ab = 3.0186998844146729;
          break;
        case 2:
          r0ab = 3.1918001174926758;
          break;
        case 3:
          r0ab = 3.5127000808715820;
          break;
        case 4:
          r0ab = 3.6875000000000000;
          break;
        case 5:
          r0ab = 3.7404000759124756;
          break;
        case 6:
          r0ab = 3.6942999362945557;
          break;
        case 7:
          r0ab = 3.4702000617980957;
          break;
        case 8:
          r0ab = 3.4888000488281250;
          break;
        case 9:
          r0ab = 3.2913999557495117;
          break;
        case 10:
          r0ab = 3.1642999649047852;
          break;
        case 11:
          r0ab = 3.6668999195098877;
          break;
        case 12:
          r0ab = 3.8724000453948975;
          break;
        case 13:
          r0ab = 3.9939999580383301;
          break;
        case 14:
          r0ab = 4.0816001892089844;
          break;
        case 15:
          r0ab = 3.8053998947143555;
          break;
        case 16:
          r0ab = 3.9660999774932861;
          break;
        case 17:
          r0ab = 3.7492001056671143;
          break;
        case 18:
          r0ab = 3.6024000644683838;
          break;
        case 19:
          r0ab = 4.0427999496459961;
          break;
        case 20:
          r0ab = 4.1950998306274414;
          break;
        case 21:
          r0ab = 4.1465997695922852;
          break;
        case 22:
          r0ab = 4.0514998435974121;
          break;
        case 23:
          r0ab = 4.0075001716613770;
          break;
        case 24:
          r0ab = 3.5020000934600830;
          break;
        case 25:
          r0ab = 3.9158000946044922;
          break;
        case 26:
          r0ab = 3.8545999526977539;
          break;
        case 27:
          r0ab = 3.8341999053955078;
          break;
        case 28:
          r0ab = 3.8008000850677490;
          break;
        case 29:
          r0ab = 3.7844998836517334;
          break;
        case 30:
          r0ab = 3.7548999786376953;
          break;
        case 31:
          r0ab = 3.9602000713348389;
          break;
        case 32:
          r0ab = 3.8872001171112061;
          break;
        case 33:
          r0ab = 3.8564000129699707;
          break;
        case 34:
          r0ab = 4.0792999267578125;
          break;
        case 35:
          r0ab = 3.8835000991821289;
          break;
        case 36:
          r0ab = 3.7495000362396240;
          break;
        case 37:
          r0ab = 4.2213001251220703;
          break;
        case 38:
          r0ab = 4.3703999519348145;
          break;
        case 39:
          r0ab = 4.3299999237060547;
          break;
        case 40:
          r0ab = 4.2121000289916992;
          break;
        case 41:
          r0ab = 3.7643001079559326;
          break;
        case 42:
          r0ab = 3.7130000591278076;
          break;
        case 43:
          r0ab = 3.6143999099731445;
          break;
        case 44:
          r0ab = 3.5599000453948975;
          break;
        case 45:
          r0ab = 3.5473999977111816;
          break;
        case 46:
          r0ab = 3.5092999935150146;
          break;
        case 47:
          r0ab = 3.4853000640869141;
          break;
        case 48:
          r0ab = 3.9075000286102295;
          break;
        case 49:
          r0ab = 4.1114997863769531;
          break;
        case 50:
          r0ab = 4.0472998619079590;
          break;
        case 51:
          r0ab = 4.0317997932434082;
          break;
        case 52:
          r0ab = 4.2999000549316406;
          break;
        case 53:
          r0ab = 4.1050000190734863;
          break;
        case 54:
          r0ab = 3.9709999561309814;
          break;
        case 55:
          r0ab = 4.4320001602172852;
          break;
        case 56:
          r0ab = 4.6705999374389648;
          break;
        case 57:
          r0ab = 4.5272998809814453;
          break;
        case 58:
          r0ab = 4.4580998420715332;
          break;
        case 59:
          r0ab = 4.4331998825073242;
          break;
        case 60:
          r0ab = 4.4064002037048340;
          break;
        case 61:
          r0ab = 4.3873000144958496;
          break;
        case 62:
          r0ab = 4.3684000968933105;
          break;
        case 63:
          r0ab = 4.3537001609802246;
          break;
        case 64:
          r0ab = 4.2727999687194824;
          break;
        case 65:
          r0ab = 4.2548999786376953;
          break;
        case 66:
          r0ab = 4.2031998634338379;
          break;
        case 67:
          r0ab = 4.2793998718261719;
          break;
        case 68:
          r0ab = 4.2613000869750977;
          break;
        case 69:
          r0ab = 4.2491002082824707;
          break;
        case 70:
          r0ab = 4.2375001907348633;
          break;
        case 71:
          r0ab = 4.2322001457214355;
          break;
        case 72:
          r0ab = 4.3664999008178711;
          break;
        case 73:
          r0ab = 4.3060998916625977;
          break;
        case 74:
          r0ab = 4.2713999748229980;
          break;
        case 75:
          r0ab = 4.2154998779296875;
          break;
        case 76:
          r0ab = 4.1416001319885254;
          break;
        case 77:
          r0ab = 3.7660000324249268;
          break;
        case 78:
          r0ab = 3.6628000736236572;
          break;
        case 79:
          r0ab = 3.5476000308990479;
          break;
        case 80:
          r0ab = 3.8789999485015869;
          break;
        case 81:
          r0ab = 4.1233000755310059;
          break;
        case 82:
          r0ab = 4.0738000869750977;
          break;
        case 83:
          r0ab = 4.0574998855590820;
          break;
        case 84:
          r0ab = 4.3575000762939453;
          break;
        case 85:
          r0ab = 4.1585998535156250;
          break;
        case 86:
          r0ab = 4.0183000564575195;
          break;
        case 87:
          r0ab = 4.4593000411987305;
          break;
        case 88:
          r0ab = 4.5927000045776367;
          break;
        case 89:
          r0ab = 4.4864997863769531;
          break;
        case 90:
          r0ab = 4.3812999725341797;
          break;
        case 91:
          r0ab = 4.4594001770019531;
          break;
      }
      break;
    case 92:
      switch (atmA) {
        case 1:
          r0ab = 2.9874999523162842;
          break;
        case 2:
          r0ab = 3.1673998832702637;
          break;
        case 3:
          r0ab = 3.4971001148223877;
          break;
        case 4:
          r0ab = 3.6714999675750732;
          break;
        case 5:
          r0ab = 3.7114000320434570;
          break;
        case 6:
          r0ab = 3.6691999435424805;
          break;
        case 7:
          r0ab = 3.4446001052856445;
          break;
        case 8:
          r0ab = 3.4676001071929932;
          break;
        case 9:
          r0ab = 3.2685000896453857;
          break;
        case 10:
          r0ab = 3.1405000686645508;
          break;
        case 11:
          r0ab = 3.6545999050140381;
          break;
        case 12:
          r0ab = 3.8578999042510986;
          break;
        case 13:
          r0ab = 3.9637000560760498;
          break;
        case 14:
          r0ab = 4.0581002235412598;
          break;
        case 15:
          r0ab = 3.7795999050140381;
          break;
        case 16:
          r0ab = 3.9463000297546387;
          break;
        case 17:
          r0ab = 3.7274999618530273;
          break;
        case 18:
          r0ab = 3.5792000293731689;
          break;
        case 19:
          r0ab = 4.0295000076293945;
          break;
        case 20:
          r0ab = 4.1824002265930176;
          break;
        case 21:
          r0ab = 4.1247000694274902;
          break;
        case 22:
          r0ab = 4.0356998443603516;
          break;
        case 23:
          r0ab = 3.9925999641418457;
          break;
        case 24:
          r0ab = 3.4827001094818115;
          break;
        case 25:
          r0ab = 3.9007000923156738;
          break;
        case 26:
          r0ab = 3.8392000198364258;
          break;
        case 27:
          r0ab = 3.8190999031066895;
          break;
        case 28:
          r0ab = 3.7850999832153320;
          break;
        case 29:
          r0ab = 3.7686998844146729;
          break;
        case 30:
          r0ab = 3.7386999130249023;
          break;
        case 31:
          r0ab = 3.9289999008178711;
          break;
        case 32:
          r0ab = 3.8605999946594238;
          break;
        case 33:
          r0ab = 3.8306000232696533;
          break;
        case 34:
          r0ab = 4.0601000785827637;
          break;
        case 35:
          r0ab = 3.8624999523162842;
          break;
        case 36:
          r0ab = 3.7269001007080078;
          break;
        case 37:
          r0ab = 4.2062001228332520;
          break;
        case 38:
          r0ab = 4.3565998077392578;
          break;
        case 39:
          r0ab = 4.3021998405456543;
          break;
        case 40:
          r0ab = 4.1929001808166504;
          break;
        case 41:
          r0ab = 3.7400999069213867;
          break;
        case 42:
          r0ab = 3.6888000965118408;
          break;
        case 43:
          r0ab = 3.5899999141693115;
          break;
        case 44:
          r0ab = 3.5350000858306885;
          break;
        case 45:
          r0ab = 3.5225999355316162;
          break;
        case 46:
          r0ab = 3.4837999343872070;
          break;
        case 47:
          r0ab = 3.4593999385833740;
          break;
        case 48:
          r0ab = 3.8887999057769775;
          break;
        case 49:
          r0ab = 4.0812997817993164;
          break;
        case 50:
          r0ab = 4.0208997726440430;
          break;
        case 51:
          r0ab = 4.0058999061584473;
          break;
        case 52:
          r0ab = 4.2810001373291016;
          break;
        case 53:
          r0ab = 4.0843000411987305;
          break;
        case 54:
          r0ab = 3.9486000537872314;
          break;
        case 55:
          r0ab = 4.4162001609802246;
          break;
        case 56:
          r0ab = 4.6542000770568848;
          break;
        case 57:
          r0ab = 4.5005002021789551;
          break;
        case 58:
          r0ab = 4.4443998336791992;
          break;
        case 59:
          r0ab = 4.4196000099182129;
          break;
        case 60:
          r0ab = 4.3933000564575195;
          break;
        case 61:
          r0ab = 4.3741002082824707;
          break;
        case 62:
          r0ab = 4.3551998138427734;
          break;
        case 63:
          r0ab = 4.3406000137329102;
          break;
        case 64:
          r0ab = 4.2484002113342285;
          break;
        case 65:
          r0ab = 4.2413001060485840;
          break;
        case 66:
          r0ab = 4.1907000541687012;
          break;
        case 67:
          r0ab = 4.2656002044677734;
          break;
        case 68:
          r0ab = 4.2473998069763184;
          break;
        case 69:
          r0ab = 4.2351999282836914;
          break;
        case 70:
          r0ab = 4.2235999107360840;
          break;
        case 71:
          r0ab = 4.2067999839782715;
          break;
        case 72:
          r0ab = 4.3410000801086426;
          break;
        case 73:
          r0ab = 4.2817001342773438;
          break;
        case 74:
          r0ab = 4.2479000091552734;
          break;
        case 75:
          r0ab = 4.1921000480651855;
          break;
        case 76:
          r0ab = 4.1181998252868652;
          break;
        case 77:
          r0ab = 3.7346000671386719;
          break;
        case 78:
          r0ab = 3.6314001083374023;
          break;
        case 79:
          r0ab = 3.5167999267578125;
          break;
        case 80:
          r0ab = 3.8582000732421875;
          break;
        case 81:
          r0ab = 4.0927000045776367;
          break;
        case 82:
          r0ab = 4.0468997955322266;
          break;
        case 83:
          r0ab = 4.0313000679016113;
          break;
        case 84:
          r0ab = 4.3390998840332031;
          break;
        case 85:
          r0ab = 4.1381001472473145;
          break;
        case 86:
          r0ab = 3.9962000846862793;
          break;
        case 87:
          r0ab = 4.4429001808166504;
          break;
        case 88:
          r0ab = 4.5787000656127930;
          break;
        case 89:
          r0ab = 4.4731001853942871;
          break;
        case 90:
          r0ab = 4.3587999343872070;
          break;
        case 91:
          r0ab = 4.4270000457763672;
          break;
        case 92:
          r0ab = 4.3956999778747559;
          break;
      }
      break;
    case 93:
      switch (atmA) {
        case 1:
          r0ab = 2.9658999443054199;
          break;
        case 2:
          r0ab = 3.1442000865936279;
          break;
        case 3:
          r0ab = 3.4795000553131104;
          break;
        case 4:
          r0ab = 3.6503000259399414;
          break;
        case 5:
          r0ab = 3.6814000606536865;
          break;
        case 6:
          r0ab = 3.6475999355316162;
          break;
        case 7:
          r0ab = 3.4221999645233154;
          break;
        case 8:
          r0ab = 3.4491000175476074;
          break;
        case 9:
          r0ab = 3.2493999004364014;
          break;
        case 10:
          r0ab = 3.1208999156951904;
          break;
        case 11:
          r0ab = 3.6324000358581543;
          break;
        case 12:
          r0ab = 3.8375000953674316;
          break;
        case 13:
          r0ab = 3.9396998882293701;
          break;
        case 14:
          r0ab = 3.8310999870300293;
          break;
        case 15:
          r0ab = 3.7581000328063965;
          break;
        case 16:
          r0ab = 3.9274001121520996;
          break;
        case 17:
          r0ab = 3.7084999084472656;
          break;
        case 18:
          r0ab = 3.5597999095916748;
          break;
        case 19:
          r0ab = 4.0079998970031738;
          break;
        case 20:
          r0ab = 4.1641001701354980;
          break;
        case 21:
          r0ab = 4.1057000160217285;
          break;
        case 22:
          r0ab = 4.0157999992370605;
          break;
        case 23:
          r0ab = 3.9725999832153320;
          break;
        case 24:
          r0ab = 3.4667000770568848;
          break;
        case 25:
          r0ab = 3.8801999092102051;
          break;
        case 26:
          r0ab = 3.8187999725341797;
          break;
        case 27:
          r0ab = 3.7988998889923096;
          break;
        case 28:
          r0ab = 3.7644000053405762;
          break;
        case 29:
          r0ab = 3.7474000453948975;
          break;
        case 30:
          r0ab = 3.7172999382019043;
          break;
        case 31:
          r0ab = 3.9049000740051270;
          break;
        case 32:
          r0ab = 3.8424000740051270;
          break;
        case 33:
          r0ab = 3.8094999790191650;
          break;
        case 34:
          r0ab = 4.0412001609802246;
          break;
        case 35:
          r0ab = 3.8436000347137451;
          break;
        case 36:
          r0ab = 3.7077000141143799;
          break;
        case 37:
          r0ab = 4.1837000846862793;
          break;
        case 38:
          r0ab = 4.3365998268127441;
          break;
        case 39:
          r0ab = 4.2815999984741211;
          break;
        case 40:
          r0ab = 4.1686000823974609;
          break;
        case 41:
          r0ab = 3.7293000221252441;
          break;
        case 42:
          r0ab = 3.6709001064300537;
          break;
        case 43:
          r0ab = 3.5699999332427979;
          break;
        case 44:
          r0ab = 3.5153000354766846;
          break;
        case 45:
          r0ab = 3.5039000511169434;
          break;
        case 46:
          r0ab = 3.4684000015258789;
          break;
        case 47:
          r0ab = 3.4437000751495361;
          break;
        case 48:
          r0ab = 3.8663001060485840;
          break;
        case 49:
          r0ab = 4.0574998855590820;
          break;
        case 50:
          r0ab = 4.0019998550415039;
          break;
        case 51:
          r0ab = 3.9842000007629395;
          break;
        case 52:
          r0ab = 4.2611999511718750;
          break;
        case 53:
          r0ab = 4.0643000602722168;
          break;
        case 54:
          r0ab = 3.9284999370574951;
          break;
        case 55:
          r0ab = 4.3927998542785645;
          break;
        case 56:
          r0ab = 4.6307997703552246;
          break;
        case 57:
          r0ab = 4.4798998832702637;
          break;
        case 58:
          r0ab = 4.4243998527526855;
          break;
        case 59:
          r0ab = 4.3996000289916992;
          break;
        case 60:
          r0ab = 4.3737001419067383;
          break;
        case 61:
          r0ab = 4.3547000885009766;
          break;
        case 62:
          r0ab = 4.3358001708984375;
          break;
        case 63:
          r0ab = 4.3211998939514160;
          break;
        case 64:
          r0ab = 4.2274999618530273;
          break;
        case 65:
          r0ab = 4.2216000556945801;
          break;
        case 66:
          r0ab = 4.1676001548767090;
          break;
        case 67:
          r0ab = 4.2465000152587891;
          break;
        case 68:
          r0ab = 4.2283000946044922;
          break;
        case 69:
          r0ab = 4.2161002159118652;
          break;
        case 70:
          r0ab = 4.2045001983642578;
          break;
        case 71:
          r0ab = 4.1841001510620117;
          break;
        case 72:
          r0ab = 4.3134999275207520;
          break;
        case 73:
          r0ab = 4.2561998367309570;
          break;
        case 74:
          r0ab = 4.2225999832153320;
          break;
        case 75:
          r0ab = 4.1666998863220215;
          break;
        case 76:
          r0ab = 4.0932002067565918;
          break;
        case 77:
          r0ab = 3.7133998870849609;
          break;
        case 78:
          r0ab = 3.6108999252319336;
          break;
        case 79:
          r0ab = 3.4962000846862793;
          break;
        case 80:
          r0ab = 3.8352000713348389;
          break;
        case 81:
          r0ab = 4.0687999725341797;
          break;
        case 82:
          r0ab = 4.0281000137329102;
          break;
        case 83:
          r0ab = 4.0099000930786133;
          break;
        case 84:
          r0ab = 4.3199000358581543;
          break;
        case 85:
          r0ab = 4.1188001632690430;
          break;
        case 86:
          r0ab = 3.9767999649047852;
          break;
        case 87:
          r0ab = 4.4191999435424805;
          break;
        case 88:
          r0ab = 4.5577001571655273;
          break;
        case 89:
          r0ab = 4.4516000747680664;
          break;
        case 90:
          r0ab = 4.3365001678466797;
          break;
        case 91:
          r0ab = 4.4057998657226562;
          break;
        case 92:
          r0ab = 4.3744997978210449;
          break;
        case 93:
          r0ab = 4.3538999557495117;
          break;
      }
      break;
    case 94:
      switch (atmA) {
        case 1:
          r0ab = 2.8763000965118408;
          break;
        case 2:
          r0ab = 3.1294000148773193;
          break;
        case 3:
          r0ab = 3.5597999095916748;
          break;
        case 4:
          r0ab = 3.7465000152587891;
          break;
        case 5:
          r0ab = 3.5659000873565674;
          break;
        case 6:
          r0ab = 3.5815999507904053;
          break;
        case 7:
          r0ab = 3.3598999977111816;
          break;
        case 8:
          r0ab = 3.4024000167846680;
          break;
        case 9:
          r0ab = 3.1877000331878662;
          break;
        case 10:
          r0ab = 3.0483999252319336;
          break;
        case 11:
          r0ab = 3.7009000778198242;
          break;
        case 12:
          r0ab = 3.9451000690460205;
          break;
        case 13:
          r0ab = 3.8464999198913574;
          break;
        case 14:
          r0ab = 3.9872999191284180;
          break;
        case 15:
          r0ab = 3.7079000473022461;
          break;
        case 16:
          r0ab = 3.9082999229431152;
          break;
        case 17:
          r0ab = 3.6756000518798828;
          break;
        case 18:
          r0ab = 3.5150001049041748;
          break;
        case 19:
          r0ab = 4.0829000473022461;
          break;
        case 20:
          r0ab = 4.2779998779296875;
          break;
        case 21:
          r0ab = 4.1511001586914062;
          break;
        case 22:
          r0ab = 4.1259999275207520;
          break;
        case 23:
          r0ab = 4.0570998191833496;
          break;
        case 24:
          r0ab = 3.4865000247955322;
          break;
        case 25:
          r0ab = 3.9744000434875488;
          break;
        case 26:
          r0ab = 3.9149999618530273;
          break;
        case 27:
          r0ab = 3.8929998874664307;
          break;
        case 28:
          r0ab = 3.8578000068664551;
          break;
        case 29:
          r0ab = 3.8401999473571777;
          break;
        case 30:
          r0ab = 3.8073000907897949;
          break;
        case 31:
          r0ab = 3.7976999282836914;
          break;
        case 32:
          r0ab = 4.0036001205444336;
          break;
        case 33:
          r0ab = 3.7604000568389893;
          break;
        case 34:
          r0ab = 4.0288000106811523;
          break;
        case 35:
          r0ab = 3.8210000991821289;
          break;
        case 36:
          r0ab = 3.6756999492645264;
          break;
        case 37:
          r0ab = 4.2645998001098633;
          break;
        case 38:
          r0ab = 4.4558000564575195;
          break;
        case 39:
          r0ab = 4.2862000465393066;
          break;
        case 40:
          r0ab = 4.2122001647949219;
          break;
        case 41:
          r0ab = 3.7088000774383545;
          break;
        case 42:
          r0ab = 3.6728999614715576;
          break;
        case 43:
          r0ab = 3.5799999237060547;
          break;
        case 44:
          r0ab = 3.5276000499725342;
          break;
        case 45:
          r0ab = 3.5164999961853027;
          break;
        case 46:
          r0ab = 3.4783000946044922;
          break;
        case 47:
          r0ab = 3.4539000988006592;
          break;
        case 48:
          r0ab = 3.9553000926971436;
          break;
        case 49:
          r0ab = 3.9818000793457031;
          break;
        case 50:
          r0ab = 4.2039999961853027;
          break;
        case 51:
          r0ab = 3.9604001045227051;
          break;
        case 52:
          r0ab = 4.2718000411987305;
          break;
        case 53:
          r0ab = 4.0689001083374023;
          break;
        case 54:
          r0ab = 3.9252998828887939;
          break;
        case 55:
          r0ab = 4.4868998527526855;
          break;
        case 56:
          r0ab = 4.7792000770568848;
          break;
        case 57:
          r0ab = 4.4917998313903809;
          break;
        case 58:
          r0ab = 4.5342001914978027;
          break;
        case 59:
          r0ab = 4.5089998245239258;
          break;
        case 60:
          r0ab = 4.4868001937866211;
          break;
        case 61:
          r0ab = 4.4679999351501465;
          break;
        case 62:
          r0ab = 4.4485998153686523;
          break;
        case 63:
          r0ab = 4.4341001510620117;
          break;
        case 64:
          r0ab = 4.2023000717163086;
          break;
        case 65:
          r0ab = 4.3122000694274902;
          break;
        case 66:
          r0ab = 4.2709999084472656;
          break;
        case 67:
          r0ab = 4.3586997985839844;
          break;
        case 68:
          r0ab = 4.3407001495361328;
          break;
        case 69:
          r0ab = 4.3281002044677734;
          break;
        case 70:
          r0ab = 4.3173999786376953;
          break;
        case 71:
          r0ab = 4.1498999595642090;
          break;
        case 72:
          r0ab = 4.3940000534057617;
          break;
        case 73:
          r0ab = 4.3895001411437988;
          break;
        case 74:
          r0ab = 4.3260002136230469;
          break;
        case 75:
          r0ab = 4.2725000381469727;
          break;
        case 76:
          r0ab = 4.1961002349853516;
          break;
        case 77:
          r0ab = 3.7360999584197998;
          break;
        case 78:
          r0ab = 3.6192998886108398;
          break;
        case 79:
          r0ab = 3.4916000366210938;
          break;
        case 80:
          r0ab = 3.9114999771118164;
          break;
        case 81:
          r0ab = 3.9914000034332275;
          break;
        case 82:
          r0ab = 3.9809000492095947;
          break;
        case 83:
          r0ab = 3.9865999221801758;
          break;
        case 84:
          r0ab = 4.3329000473022461;
          break;
        case 85:
          r0ab = 4.1276001930236816;
          break;
        case 86:
          r0ab = 3.9781999588012695;
          break;
        case 87:
          r0ab = 4.5096998214721680;
          break;
        case 88:
          r0ab = 4.6768999099731445;
          break;
        case 89:
          r0ab = 4.5157999992370605;
          break;
        case 90:
          r0ab = 4.3291001319885254;
          break;
        case 91:
          r0ab = 4.3608999252319336;
          break;
        case 92:
          r0ab = 4.3461999893188477;
          break;
        case 93:
          r0ab = 4.3264999389648438;
          break;
        case 94:
          r0ab = 4.4341001510620117;
          break;
      }
    break;
  }
  return r0ab*dist_Angstrom2aum1;
}
double R2R4(int atmnr) {
  //the r2r4 parameters of D3 and D4, defined as \sqrt{\frac{1}{2} r2r4(A) {\sqrt{{Z}_{A}}}}
  double r2r4 = 0.0;
  switch (atmnr) {
    case 1:
      r2r4 = 2.00734899805689;
      break;
    case 2:
      r2r4 = 1.5663713191197939;
      break;
    case 3:
      r2r4 = 5.019869279580627;
      break;
    case 4:
      r2r4 = 3.8537903419880015;
      break;
    case 5:
      r2r4 = 3.6444659393318353;
      break;
    case 6:
      r2r4 = 3.104928220612981;
      break;
    case 7:
      r2r4 = 2.7117524212163846;
      break;
    case 8:
      r2r4 = 2.5936168242020377;
      break;
    case 9:
      r2r4 = 2.3882524992135985;
      break;
    case 10:
      r2r4 = 2.2152251517733568;
      break;
    case 11:
      r2r4 = 6.585855476576431;
      break;
    case 12:
      r2r4 = 5.462959658896886;
      break;
    case 13:
      r2r4 = 5.65216658249283;
      break;
    case 14:
      r2r4 = 4.882849069371511;
      break;
    case 15:
      r2r4 = 4.297275676970734;
      break;
    case 16:
      r2r4 = 4.041088962148693;
      break;
    case 17:
      r2r4 = 3.7293234964336643;
      break;
    case 18:
      r2r4 = 3.4467727589836006;
      break;
    case 19:
      r2r4 = 7.977627461536416;
      break;
    case 20:
      r2r4 = 7.076239442482412;
      break;
    case 21:
      r2r4 = 6.608440672001199;
      break;
    case 22:
      r2r4 = 6.287913775445048;
      break;
    case 23:
      r2r4 = 6.077287053921326;
      break;
    case 24:
      r2r4 = 5.546430990343218;
      break;
    case 25:
      r2r4 = 5.8049117133682575;
      break;
    case 26:
      r2r4 = 5.584156058385277;
      break;
    case 27:
      r2r4 = 5.413745277777264;
      break;
    case 28:
      r2r4 = 5.284972281457849;
      break;
    case 29:
      r2r4 = 5.225928128371258;
      break;
    case 30:
      r2r4 = 5.098171474418196;
      break;
    case 31:
      r2r4 = 6.121496915882168;
      break;
    case 32:
      r2r4 = 5.540837426133801;
      break;
    case 33:
      r2r4 = 5.066968900257076;
      break;
    case 34:
      r2r4 = 4.870051009207528;
      break;
    case 35:
      r2r4 = 4.590896430032106;
      break;
    case 36:
      r2r4 = 4.311762980498812;
      break;
    case 37:
      r2r4 = 9.554617121107926;
      break;
    case 38:
      r2r4 = 8.673960892714286;
      break;
    case 39:
      r2r4 = 7.972101848320851;
      break;
    case 40:
      r2r4 = 7.434399030481411;
      break;
    case 41:
      r2r4 = 6.587118572072311;
      break;
    case 42:
      r2r4 = 6.195362041346371;
      break;
    case 43:
      r2r4 = 6.0151730164200075;
      break;
    case 44:
      r2r4 = 5.8162339941054615;
      break;
    case 45:
      r2r4 = 5.657104250477866;
      break;
    case 46:
      r2r4 = 5.5264067015115925;
      break;
    case 47:
      r2r4 = 5.442633110925525;
      break;
    case 48:
      r2r4 = 5.582853596701417;
      break;
    case 49:
      r2r4 = 7.020819040539359;
      break;
    case 50:
      r2r4 = 6.4681553296959295;
      break;
    case 51:
      r2r4 = 5.980891064644942;
      break;
    case 52:
      r2r4 = 5.816866457936288;
      break;
    case 53:
      r2r4 = 5.5332180636017965;
      break;
    case 54:
      r2r4 = 5.254770020890183;
      break;
    case 55:
      r2r4 = 11.022045592219385;
      break;
    case 56:
      r2r4 = 10.156795138195266;
      break;
    case 57:
      r2r4 = 9.351678174139682;
      break;
    case 58:
      r2r4 = 9.069260824857563;
      break;
    case 59:
      r2r4 = 8.972411509564969;
      break;
    case 60:
      r2r4 = 8.900928149731811;
      break;
    case 61:
      r2r4 = 8.859848439318803;
      break;
    case 62:
      r2r4 = 8.817368385523206;
      break;
    case 63:
      r2r4 = 8.793177211043448;
      break;
    case 64:
      r2r4 = 7.8996961966900985;
      break;
    case 65:
      r2r4 = 8.805884467156623;
      break;
    case 66:
      r2r4 = 8.424392019436123;
      break;
    case 67:
      r2r4 = 8.542892635498433;
      break;
    case 68:
      r2r4 = 8.475833616174544;
      break;
    case 69:
      r2r4 = 8.450908826050055;
      break;
    case 70:
      r2r4 = 8.473393564979439;
      break;
    case 71:
      r2r4 = 7.8352563466301675;
      break;
    case 72:
      r2r4 = 8.207028391602428;
      break;
    case 73:
      r2r4 = 7.705590621613434;
      break;
    case 74:
      r2r4 = 7.327559868011995;
      break;
    case 75:
      r2r4 = 7.038873941420576;
      break;
    case 76:
      r2r4 = 6.68978709251845;
      break;
    case 77:
      r2r4 = 6.05450039062772;
      break;
    case 78:
      r2r4 = 5.887520280063885;
      break;
    case 79:
      r2r4 = 5.706615067235794;
      break;
    case 80:
      r2r4 = 5.784506998466409;
      break;
    case 81:
      r2r4 = 7.79780738413049;
      break;
    case 82:
      r2r4 = 7.264438585983326;
      break;
    case 83:
      r2r4 = 6.78151997607037;
      break;
    case 84:
      r2r4 = 6.678831620478252;
      break;
    case 85:
      r2r4 = 6.390243055646761;
      break;
    case 86:
      r2r4 = 6.095279657044429;
      break;
    case 87:
      r2r4 = 11.791560867825417;
      break;
    case 88:
      r2r4 = 11.109976328092937;
      break;
    case 89:
      r2r4 = 9.513778086345827;
      break;
    case 90:
      r2r4 = 8.671970578810924;
      break;
    case 91:
      r2r4 = 8.7714070389123;
      break;
    case 92:
      r2r4 = 8.654027077101373;
      break;
    case 93:
      r2r4 = 8.539235123774027;
      break;
    case 94:
      r2r4 = 8.850247005829276;
      break;
    case 95:
      break;
    case 96:
      break;
    case 97:
      break;
    case 98:
      break;
    case 99:
      break;
    case 100:
      break;
    case 101:
      break;
    case 102:
      break;
    case 103:
      break;
    case 104:
      break;
    case 105:
      break;
    case 106:
      break;
    case 107:
      break;
    case 108:
      break;
    case 109:
      break;
    case 110:
      break;
    case 111:
      break;
    case 112:
      r2r4 = 5.3912609613237406;
      break;
    case 113:
      r2r4 = 5.980217350379766;
      break;
    case 114:
      r2r4 = 5.897234884462072;
      break;
    case 115:
      r2r4 = 7.650840150676422;
      break;
    case 116:
      r2r4 = 7.4692909571065345;
      break;
    case 117:
      r2r4 = 7.157464287012987;
      break;
    case 118:
      r2r4 = 6.859891269320811;
      break;
  }
  return r2r4;
}

#endif //_Dispersion3_Parameters_
