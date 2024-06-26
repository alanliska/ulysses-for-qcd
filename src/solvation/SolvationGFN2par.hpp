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

#ifndef _Solvation_Parameters_for_GFN2_
#define _Solvation_Parameters_for_GFN2_
#include <vector>
#include <string>
#include "../UnitConversion.hpp"

//descrition:
//parameters to auxiliate the calculation of charges of atoms in molecules for GFN2

double FreeEnergyShiftGFN2(std::string solvent) {
  //function that returns the free energy shifts for generalized Born model of solvation for GFN2 method
  //data taken from xTB, stored in kcal/mol
  double shift = 0.0;
  if (solvent == "acetone") {shift = 3.71745477;}
  else if (solvent == "acetonitrile") {shift = 3.26220511;}
  else if (solvent == "aniline") {shift = 1.86522103;}
  else if (solvent == "benzaldehyde") {shift = 1.41868733;}
  else if (solvent == "benzene") {shift = 2.93326080;}
  else if (solvent == "dichloromethane") {shift = 2.60334690;}
  else if (solvent == "chloroform") {shift = 2.44165144;}
  else if (solvent == "carbon disulfide") {shift = 3.12529087;}
  else if (solvent == "dioxane") {shift = 1.73889651;}
  else if (solvent == "dmf") {shift = 3.77045120;}
  else if (solvent == "dmso") {shift = 4.51473868;}
  else if (solvent == "ethanol") {shift = 2.14060025;}
  else if (solvent == "diethyl ether") {shift = 3.10832901;}
  else if (solvent == "ethyl acetate") {shift = 1.51781785;}
  else if (solvent == "furane") {shift = 1.24661341;}
  else if (solvent == "hexadecane") {shift = 2.75911340;}
  else if (solvent == "hexane") {shift = 3.02042638;}
  else if (solvent == "methanol") {shift = 2.47547456;}
  else if (solvent == "nitromethane") {shift = 1.51623302;}
  else if (solvent == "octanol") {shift = 1.46816470;}
  else if (solvent == "phenol") {shift = 1.52339711;}
  else if (solvent == "thf") {shift = 3.34137214;}
  else if (solvent == "toluene") {shift = 1.38563074;}
  else if (solvent == "water") {shift = 0.67818695;}
  else if (solvent == "octanol wet") {shift = 1.57636999;}
  return shift/au2kcalmol;
}
double BornRadiusScaleGFN2(std::string solvent) {
  //function that returns the scaling of Born radii for generalized Born model of solvation for GFN2 method
  //data taken from xTB
  double scale = 0.0;
  if (solvent == "acetone") {scale = 1.82744517;}
  else if (solvent == "acetonitrile") {scale = 1.72136166;}
  else if (solvent == "aniline") {scale = 1.48549561;}
  else if (solvent == "benzaldehyde") {scale = 1.81873738;}
  else if (solvent == "benzene") {scale = 1.69138605;}
  else if (solvent == "dichloromethane") {scale = 1.56911675;}
  else if (solvent == "chloroform") {scale = 1.32742806;}
  else if (solvent == "carbon disulfide") {scale = 1.35477845;}
  else if (solvent == "dioxane") {scale = 1.08594401;}
  else if (solvent == "dmf") {scale = 1.85402904;}
  else if (solvent == "dmso") {scale = 1.79272757;}
  else if (solvent == "ethanol") {scale = 1.62446646;}
  else if (solvent == "diethyl ether") {scale = 1.58465468;}
  else if (solvent == "ethyl acetate") {scale = 1.64733795;}
  else if (solvent == "furane") {scale = 1.28872589;}
  else if (solvent == "hexadecane") {scale = 1.25929425;}
  else if (solvent == "hexane") {scale = 1.10400200;}
  else if (solvent == "methanol") {scale = 1.66003361;}
  else if (solvent == "nitromethane") {scale = 1.70029849;}
  else if (solvent == "octanol") {scale = 1.58983390;}
  else if (solvent == "phenol") {scale = 1.29875185;}
  else if (solvent == "thf") {scale = 1.73108801;}
  else if (solvent == "toluene") {scale = 1.71719241;}
  else if (solvent == "water") {scale = 1.47438678;}
  else if (solvent == "octanol wet") {scale = 1.61503543;}
  return scale;
}
double ProbeRadiusGFN2(std::string solvent) {
  //function that returns the probe radii for generalized Born model of solvation for GFN2 method
  //data taken from xTB, stored and returned in Angstroem
  double proberad = 0.0;
  if (solvent == "acetone") {proberad = 1.16056371;}
  else if (solvent == "acetonitrile") {proberad = 1.02173881;}
  else if (solvent == "aniline") {proberad = 1.07420745;}
  else if (solvent == "benzaldehyde") {proberad = 0.27702266;}
  else if (solvent == "benzene") {proberad = 0.51651114;}
  else if (solvent == "dichloromethane") {proberad = 0.45821664;}
  else if (solvent == "chloroform") {proberad = 0.75620958;}
  else if (solvent == "carbon disulfide") {proberad = 0.57997496;}
  else if (solvent == "dioxane") {proberad = 1.00107567;}
  else if (solvent == "dmf") {proberad = 1.28494293;}
  else if (solvent == "dmso") {proberad = 1.49215217;}
  else if (solvent == "ethanol") {proberad = 0.91279578;}
  else if (solvent == "diethyl ether") {proberad = 1.14930964;}
  else if (solvent == "ethyl acetate") {proberad = 0.50241618;}
  else if (solvent == "furane") {proberad = 0.48900006;}
  else if (solvent == "hexadecane") {proberad = 0.55336799;}
  else if (solvent == "hexane") {proberad = 0.66693748;}
  else if (solvent == "methanol") {proberad = 1.08058328;}
  else if (solvent == "nitromethane") {proberad = 0.25325726;}
  else if (solvent == "octanol") {proberad = 0.49307217;}
  else if (solvent == "phenol") {proberad = 1.18848348;}
  else if (solvent == "thf") {proberad = 1.22908546;}
  else if (solvent == "toluene") {proberad = 0.09081163;}
  else if (solvent == "water") {proberad = 1.13409020;}
  else if (solvent == "octanol wet") {proberad = 0.48474691;}
  return proberad;
}
double GammaScaleGFN2(std::string solvent, size_t atmnr) {
  //function that returns the gamma scale for generalized Born model of solvation for GFN2 method
  //data taken from xTB; missing parameters are zero, but since they are not on relevant atoms for common applications,
  //I decided to keep these out
  double gammascal = 0.0;
  if (solvent == "acetone") {
    switch (atmnr) {
      case 1: 
        gammascal = -1.58616844;
        break;
      case 2: 
        gammascal = -0.30709560;
        break;
      case 3: 
        gammascal = -5.94928113;
        break;
      case 4: 
        gammascal = -13.65649967;
        break;
      case 5: 
        gammascal = -2.50699930;
        break;
      case 6: 
        gammascal = -2.56310412;
        break;
      case 7: 
        gammascal = -3.15675044;
        break;
      case 8: 
        gammascal = -2.44808825;
        break;
      case 9: 
        gammascal = -0.89945455;
        break;
      case 10: 
        gammascal = -0.33634476;
        break;
      case 11: 
        gammascal = -3.87007975;
        break;
      case 12: 
        gammascal = -2.72929146;
        break;
      case 13: 
        gammascal = -2.04849151;
        break;
      case 14: 
        gammascal = 0.17586358;
        break;
      case 15: 
        gammascal = -2.20812285;
        break;
      case 16: 
        gammascal = -2.31423806;
        break;
      case 17: 
        gammascal = -1.73119355;
        break;
      case 18: 
        gammascal = -1.28889221;
        break;
      case 19: 
        gammascal = -3.21905703;
        break;
      case 20: 
        gammascal = -6.03598657;
        break;
      case 21: 
        gammascal = -5.57005246;
        break;
      case 22: 
        gammascal = -4.25136024;
        break;
      case 23: 
        gammascal = -4.31991206;
        break;
      case 24: 
        gammascal = -2.99519395;
        break;
      case 25: 
        gammascal = -2.02492855;
        break;
      case 26: 
        gammascal = 1.71668551;
        break;
      case 27: 
        gammascal = 11.68225676;
        break;
      case 28: 
        gammascal = -1.82459822;
        break;
      case 29: 
        gammascal = -2.05246870;
        break;
      case 30: 
        gammascal = -3.66452610;
        break;
      case 31: 
        gammascal = -2.77606678;
        break;
      case 32: 
        gammascal = -1.88897991;
        break;
      case 33: 
        gammascal = -2.29491449;
        break;
      case 34: 
        gammascal = -1.82087329;
        break;
      case 35: 
        gammascal = -2.34857165;
        break;
      case 36: 
        gammascal = -1.31367125;
        break;
      case 37: 
        gammascal = -3.22967523;
        break;
      case 38: 
        gammascal = -6.92349673;
        break;
      case 39: 
        gammascal = -8.03352261;
        break;
      case 40: 
        gammascal = -0.95352196;
        break;
      case 41: 
        gammascal = -4.07110792;
        break;
      case 42: 
        gammascal = -4.02516701;
        break;
      case 43: 
        gammascal = 2.24676734;
        break;
      case 44: 
        gammascal = -2.78928205;
        break;
      case 45: 
        gammascal = -1.53836954;
        break;
      case 46: 
        gammascal = -2.15615111;
        break;
      case 47: 
        gammascal = -2.34319672;
        break;
      case 48: 
        gammascal = -2.51689036;
        break;
      case 49: 
        gammascal = -3.12829238;
        break;
      case 50: 
        gammascal = -1.59731993;
        break;
      case 51: 
        gammascal = -1.86791451;
        break;
      case 52: 
        gammascal = -2.16601206;
        break;
      case 53: 
        gammascal = -2.94891199;
        break;
      case 54: 
        gammascal = -1.12781227;
        break;
      case 55: 
        gammascal = -3.18229731;
        break;
      case 56: 
        gammascal = -5.79886288;
        break;
      case 72: 
        gammascal = -3.02379488;
        break;
      case 73: 
        gammascal = -3.21255866;
        break;
      case 74: 
        gammascal = 2.57275593;
        break;
      case 75: 
        gammascal = 0.69395644;
        break;
      case 76: 
        gammascal = -0.79978668;
        break;
      case 77: 
        gammascal = -2.02076917;
        break;
      case 78: 
        gammascal = -0.95964178;
        break;
      case 79: 
        gammascal = -1.31659237;
        break;
      case 80: 
        gammascal = -2.71478938;
        break;
      case 81: 
        gammascal = -2.12346632;
        break;
      case 82: 
        gammascal = -0.90197761;
        break;
      case 83: 
        gammascal = -1.40814465;
        break;
    }
  }
  else if (solvent == "acetonitrile") {
    switch (atmnr) {
      case 1: 
        gammascal = -1.30829234;
        break;
      case 2: 
        gammascal = 0.07869756;
        break;
      case 3: 
        gammascal = -6.34782905;
        break;
      case 4: 
        gammascal = -15.52350747;
        break;
      case 5: 
        gammascal = -2.39455571;
        break;
      case 6: 
        gammascal = -2.52921042;
        break;
      case 7: 
        gammascal = -3.51706691;
        break;
      case 8: 
        gammascal = -2.39172770;
        break;
      case 9: 
        gammascal = -0.60311200;
        break;
      case 10: 
        gammascal = 0.04022038;
        break;
      case 11: 
        gammascal = -4.14669337;
        break;
      case 12: 
        gammascal = -3.26116104;
        break;
      case 13: 
        gammascal = -1.95553354;
        break;
      case 14: 
        gammascal = 0.42302196;
        break;
      case 15: 
        gammascal = -2.07953398;
        break;
      case 16: 
        gammascal = -2.20329733;
        break;
      case 17: 
        gammascal = -1.52271664;
        break;
      case 18: 
        gammascal = -0.87279161;
        break;
      case 19: 
        gammascal = -3.33652372;
        break;
      case 20: 
        gammascal = -7.13358295;
        break;
      case 21: 
        gammascal = -7.53566194;
        break;
      case 22: 
        gammascal = -4.56573371;
        break;
      case 23: 
        gammascal = -5.13232712;
        break;
      case 24: 
        gammascal = -2.99559920;
        break;
      case 25: 
        gammascal = -2.17079646;
        break;
      case 26: 
        gammascal = 1.03547747;
        break;
      case 27: 
        gammascal = 15.00866716;
        break;
      case 28: 
        gammascal = -1.96685553;
        break;
      case 29: 
        gammascal = -1.84653361;
        break;
      case 30: 
        gammascal = -3.90820543;
        break;
      case 31: 
        gammascal = -2.96022711;
        break;
      case 32: 
        gammascal = -1.45693350;
        break;
      case 33: 
        gammascal = -2.19372337;
        break;
      case 34: 
        gammascal = -1.60615814;
        break;
      case 35: 
        gammascal = -1.86875461;
        break;
      case 36: 
        gammascal = -0.90835488;
        break;
      case 37: 
        gammascal = -3.37828045;
        break;
      case 38: 
        gammascal = -8.14805044;
        break;
      case 39: 
        gammascal = -9.32105410;
        break;
      case 40: 
        gammascal = -1.06999861;
        break;
      case 41: 
        gammascal = -4.97269666;
        break;
      case 42: 
        gammascal = -4.28496846;
        break;
      case 43: 
        gammascal = 0.65739714;
        break;
      case 44: 
        gammascal = -3.31577377;
        break;
      case 45: 
        gammascal = -1.26229117;
        break;
      case 46: 
        gammascal = -2.37622966;
        break;
      case 47: 
        gammascal = -2.16541127;
        break;
      case 48: 
        gammascal = -2.89715920;
        break;
      case 49: 
        gammascal = -3.27231922;
        break;
      case 50: 
        gammascal = -1.28366406;
        break;
      case 51: 
        gammascal = -1.77091204;
        break;
      case 52: 
        gammascal = -2.01166126;
        break;
      case 53: 
        gammascal = -2.50433473;
        break;
      case 54: 
        gammascal = -0.76859785;
        break;
      case 55: 
        gammascal = -3.15818239;
        break;
      case 56: 
        gammascal = -6.65964080;
        break;
      case 72: 
        gammascal = -3.42834747;
        break;
      case 73: 
        gammascal = -3.84978423;
        break;
      case 74: 
        gammascal = 0.90801348;
        break;
      case 75: 
        gammascal = 0.58973148;
        break;
      case 76: 
        gammascal = -1.08176895;
        break;
      case 77: 
        gammascal = -1.80527023;
        break;
      case 78: 
        gammascal = -0.90475874;
        break;
      case 79: 
        gammascal = -1.06282789;
        break;
      case 80: 
        gammascal = -2.86573574;
        break;
      case 81: 
        gammascal = -2.06342409;
        break;
      case 82: 
        gammascal = -0.58497698;
        break;
      case 83: 
        gammascal = -1.14790964;
        break;
    }
  }
  else if (solvent == "aniline") {
    switch (atmnr) {
      case 1:
        gammascal = -1.11127792;
        break;
      case 2:
        gammascal = -0.00755122;
        break;
      case 3:
        gammascal = -4.78401252;
        break;
      case 4:
        gammascal = -14.18157696;
        break;
      case 5:
        gammascal = -2.14093643;
        break;
      case 6:
        gammascal = -2.26934003;
        break;
      case 7:
        gammascal = -2.37133105;
        break;
      case 8:
        gammascal = -0.25020795;
        break;
      case 9:
        gammascal = 0.29186478;
        break;
      case 10:
        gammascal = -0.06627081;
        break;
      case 11:
        gammascal = -3.37712198;
        break;
      case 12:
        gammascal = -2.29418433;
        break;
      case 13:
        gammascal = -1.55411577;
        break;
      case 14:
        gammascal = 0.52667190;
        break;
      case 15:
        gammascal = -1.89402854;
        break;
      case 16:
        gammascal = -1.48875065;
        break;
      case 17:
        gammascal = -1.28068651;
        break;
      case 18:
        gammascal = -0.86024457;
        break;
      case 19:
        gammascal = -2.83699421;
        break;
      case 20:
        gammascal = -6.84000195;
        break;
      case 21:
        gammascal = -6.27688128;
        break;
      case 22:
        gammascal = -3.13882726;
        break;
      case 23:
        gammascal = -5.38467709;
        break;
      case 24:
        gammascal = -2.18311181;
        break;
      case 25:
        gammascal = -2.80545653;
        break;
      case 26:
        gammascal = -0.19204474;
        break;
      case 27:
        gammascal = 1.00000000;
        break;
      case 28:
        gammascal = -1.02971367;
        break;
      case 29:
        gammascal = -1.76075659;
        break;
      case 30:
        gammascal = -2.99671708;
        break;
      case 31:
        gammascal = -2.76891595;
        break;
      case 32:
        gammascal = -1.46369860;
        break;
      case 33:
        gammascal = -1.98160058;
        break;
      case 34:
        gammascal = -1.16424874;
        break;
      case 35:
        gammascal = -1.50974562;
        break;
      case 36:
        gammascal = -0.92107172;
        break;
      case 37:
        gammascal = -3.17418901;
        break;
      case 38:
        gammascal = -7.91243730;
        break;
      case 39:
        gammascal = -7.95108977;
        break;
      case 40:
        gammascal = 0.29924239;
        break;
      case 41:
        gammascal = -5.65543732;
        break;
      case 42:
        gammascal = -2.41710525;
        break;
      case 43:
        gammascal = -2.86025348;
        break;
      case 44:
        gammascal = -3.44412426;
        break;
      case 45:
        gammascal = -1.90499988;
        break;
      case 46:
        gammascal = -1.96078529;
        break;
      case 47:
        gammascal = -1.95466975;
        break;
      case 48:
        gammascal = -2.04487534;
        break;
      case 49:
        gammascal = -3.06482026;
        break;
      case 50:
        gammascal = -1.40793726;
        break;
      case 51:
        gammascal = -1.68527092;
        break;
      case 52:
        gammascal = -1.84822443;
        break;
      case 53:
        gammascal = -2.05159327;
        break;
      case 54:
        gammascal = -0.87797886;
        break;
      case 55:
        gammascal = -3.15507821;
        break;
      case 56:
        gammascal = -6.36518430;
        break;
      case 72:
        gammascal = -2.10301188;
        break;
      case 73:
        gammascal = -4.14190735;
        break;
      case 74:
        gammascal = 3.56821358;
        break;
      case 75:
        gammascal = -3.95574473;
        break;
      case 76:
        gammascal = -1.76679550;
        break;
      case 77:
        gammascal = -1.89776298;
        break;
      case 78:
        gammascal = -0.74504857;
        break;
      case 79:
        gammascal = -0.84934565;
        break;
      case 80:
        gammascal = -2.65870028;
        break;
      case 81:
        gammascal = -1.60578659;
        break;
      case 82:
        gammascal = -0.71401912;
        break;
      case 83:
        gammascal = -0.44774034;
        break;
    }
  }
  else if (solvent == "benzaldehyde") {
    switch (atmnr) {
      case 1:
        gammascal = -1.46702255;
        break;
      case 2:
        gammascal = -0.06148252;
        break;
      case 3:
        gammascal = -7.66687950;
        break;
      case 4:
        gammascal = -21.04489640;
        break;
      case 5:
        gammascal = -3.19790444;
        break;
      case 6:
        gammascal = -2.87023244;
        break;
      case 7:
        gammascal = -4.68567000;
        break;
      case 8:
        gammascal = -3.54578103;
        break;
      case 9:
        gammascal = -0.36961486;
        break;
      case 10:
        gammascal = -0.09580315;
        break;
      case 11:
        gammascal = -5.15518940;
        break;
      case 12:
        gammascal = -3.73204762;
        break;
      case 13:
        gammascal = -1.88800489;
        break;
      case 14:
        gammascal = 0.25660828;
        break;
      case 15:
        gammascal = -2.54580556;
        break;
      case 16:
        gammascal = -2.78674456;
        break;
      case 17:
        gammascal = -1.79878573;
        break;
      case 18:
        gammascal = -1.10277999;
        break;
      case 19:
        gammascal = -3.40792693;
        break;
      case 20:
        gammascal = -7.23574094;
        break;
      case 21:
        gammascal = -3.65302081;
        break;
      case 22:
        gammascal = -2.92425016;
        break;
      case 23:
        gammascal = -5.00985749;
        break;
      case 24:
        gammascal = 1.31782112;
        break;
      case 25:
        gammascal = 4.57672665;
        break;
      case 26:
        gammascal = 3.74920531;
        break;
      case 27:
        gammascal = 14.83823511;
        break;
      case 28:
        gammascal = -0.91784711;
        break;
      case 29:
        gammascal = -2.00332067;
        break;
      case 30:
        gammascal = -4.88556958;
        break;
      case 31:
        gammascal = -3.45892351;
        break;
      case 32:
        gammascal = -1.55737795;
        break;
      case 33:
        gammascal = -2.74248296;
        break;
      case 34:
        gammascal = -2.00878407;
        break;
      case 35:
        gammascal = -2.04723148;
        break;
      case 36:
        gammascal = -1.14305823;
        break;
      case 37:
        gammascal = -3.20791240;
        break;
      case 38:
        gammascal = -9.22235201;
        break;
      case 39:
        gammascal = -10.42734295;
        break;
      case 40:
        gammascal = -3.07438550;
        break;
      case 41:
        gammascal = -5.74242280;
        break;
      case 42:
        gammascal = -2.71645716;
        break;
      case 43:
        gammascal = 5.08899389;
        break;
      case 44:
        gammascal = -1.02644353;
        break;
      case 45:
        gammascal = 1.19774424;
        break;
      case 46:
        gammascal = -2.32241807;
        break;
      case 47:
        gammascal = -2.58800777;
        break;
      case 48:
        gammascal = -3.22784738;
        break;
      case 49:
        gammascal = -3.67099805;
        break;
      case 50:
        gammascal = -1.38846468;
        break;
      case 51:
        gammascal = -2.44562704;
        break;
      case 52:
        gammascal = -2.47775416;
        break;
      case 53:
        gammascal = -2.52758619;
        break;
      case 54:
        gammascal = -0.71329061;
        break;
      case 55:
        gammascal = -1.17609485;
        break;
      case 56:
        gammascal = -7.05060535;
        break;
      case 72:
        gammascal = -4.37875898;
        break;
      case 73:
        gammascal = -5.32093886;
        break;
      case 74:
        gammascal = 3.97697453;
        break;
      case 75:
        gammascal = 5.04415442;
        break;
      case 76:
        gammascal = -0.60943167;
        break;
      case 77:
        gammascal = 0.57381406;
        break;
      case 78:
        gammascal = -1.26494355;
        break;
      case 79:
        gammascal = -1.05946557;
        break;
      case 80:
        gammascal = -2.59415591;
        break;
      case 81:
        gammascal = -2.60250239;
        break;
      case 82:
        gammascal = -0.89855659;
        break;
      case 83:
        gammascal = -1.97691485;
        break;
    }
  }
  else if (solvent == "benzene") {
    switch (atmnr) {
      case 1:
        gammascal = -1.82491712;
        break;
      case 2:
        gammascal = -0.22190385;
        break;
      case 3:
        gammascal = -5.64774682;
        break;
      case 4:
        gammascal = -15.56185215;
        break;
      case 5:
        gammascal = -3.31382345;
        break;
      case 6:
        gammascal = -3.02635919;
        break;
      case 7:
        gammascal = -3.80560008;
        break;
      case 8:
        gammascal = -2.95045534;
        break;
      case 9:
        gammascal = -0.88599329;
        break;
      case 10:
        gammascal = -0.23739093;
        break;
      case 11:
        gammascal = -4.07769540;
        break;
      case 12:
        gammascal = -3.29476667;
        break;
      case 13:
        gammascal = -2.47947364;
        break;
      case 14:
        gammascal = -0.03922808;
        break;
      case 15:
        gammascal = -3.07265214;
        break;
      case 16:
        gammascal = -2.97425165;
        break;
      case 17:
        gammascal = -2.13703079;
        break;
      case 18:
        gammascal = -1.29380105;
        break;
      case 19:
        gammascal = -3.16966887;
        break;
      case 20:
        gammascal = -5.29716763;
        break;
      case 21:
        gammascal = -5.01022748;
        break;
      case 22:
        gammascal = -3.01979201;
        break;
      case 23:
        gammascal = -4.04930089;
        break;
      case 24:
        gammascal = -0.64154702;
        break;
      case 25:
        gammascal = -1.06418916;
        break;
      case 26:
        gammascal = -0.57559699;
        break;
      case 27:
        gammascal = 7.32318069;
        break;
      case 28:
        gammascal = -2.47299450;
        break;
      case 29:
        gammascal = -2.35486094;
        break;
      case 30:
        gammascal = -4.03524732;
        break;
      case 31:
        gammascal = -3.38965614;
        break;
      case 32:
        gammascal = -2.18350982;
        break;
      case 33:
        gammascal = -2.85131617;
        break;
      case 34:
        gammascal = -2.23764222;
        break;
      case 35:
        gammascal = -2.29085994;
        break;
      case 36:
        gammascal = -1.36448221;
        break;
      case 37:
        gammascal = -3.19544832;
        break;
      case 38:
        gammascal = -6.32582126;
        break;
      case 39:
        gammascal = -8.29616111;
        break;
      case 40:
        gammascal = -2.23752420;
        break;
      case 41:
        gammascal = -4.59687421;
        break;
      case 42:
        gammascal = -2.78786191;
        break;
      case 43:
        gammascal = 1.88755867;
        break;
      case 44:
        gammascal = -2.26504777;
        break;
      case 45:
        gammascal = -1.08358365;
        break;
      case 46:
        gammascal = -3.00792883;
        break;
      case 47:
        gammascal = -2.86224027;
        break;
      case 48:
        gammascal = -3.28733192;
        break;
      case 49:
        gammascal = -3.45415205;
        break;
      case 50:
        gammascal = -1.73360669;
        break;
      case 51:
        gammascal = -2.45173756;
        break;
      case 52:
        gammascal = -2.60607353;
        break;
      case 53:
        gammascal = -2.87046623;
        break;
      case 54:
        gammascal = -0.82010829;
        break;
      case 55:
        gammascal = -2.70187431;
        break;
      case 56:
        gammascal = -5.18419436;
        break;
      case 72:
        gammascal = -3.92108430;
        break;
      case 73:
        gammascal = -4.05877218;
        break;
      case 74:
        gammascal = 2.82019057;
        break;
      case 75:
        gammascal = 1.44160582;
        break;
      case 76:
        gammascal = -0.53266645;
        break;
      case 77:
        gammascal = -2.05711862;
        break;
      case 78:
        gammascal = -1.54216144;
        break;
      case 79:
        gammascal = -1.59857149;
        break;
      case 80:
        gammascal = -2.35979380;
        break;
      case 81:
        gammascal = -2.31965110;
        break;
      case 82:
        gammascal = -1.10448115;
        break;
      case 83:
        gammascal = -1.67224873;
        break;
    }
  }
  else if (solvent == "dichloromethane") {
    switch (atmnr) {
      case 1:
        gammascal = -1.87729929;
        break;
      case 2:
        gammascal = -0.25023386;
        break;
      case 3:
        gammascal = -6.20317370;
        break;
      case 4:
        gammascal = -17.83076437;
        break;
      case 5:
        gammascal = -3.49513745;
        break;
      case 6:
        gammascal = -2.69888049;
        break;
      case 7:
        gammascal = -4.62381563;
        break;
      case 8:
        gammascal = -2.59784469;
        break;
      case 9:
        gammascal = -0.44703707;
        break;
      case 10:
        gammascal = -0.26155304;
        break;
      case 11:
        gammascal = -4.65721100;
        break;
      case 12:
        gammascal = -3.31383144;
        break;
      case 13:
        gammascal = -2.57206719;
        break;
      case 14:
        gammascal = -0.04078322;
        break;
      case 15:
        gammascal = -3.17353079;
        break;
      case 16:
        gammascal = -3.04933155;
        break;
      case 17:
        gammascal = -2.08989058;
        break;
      case 18:
        gammascal = -1.35119955;
        break;
      case 19:
        gammascal = -3.57321863;
        break;
      case 20:
        gammascal = -7.44660287;
        break;
      case 21:
        gammascal = -5.39298648;
        break;
      case 22:
        gammascal = -3.39519575;
        break;
      case 23:
        gammascal = -4.72991591;
        break;
      case 24:
        gammascal = -0.63470925;
        break;
      case 25:
        gammascal = -1.92305209;
        break;
      case 26:
        gammascal = -0.70471468;
        break;
      case 27:
        gammascal = 10.17280212;
        break;
      case 28:
        gammascal = -2.15635298;
        break;
      case 29:
        gammascal = -2.47307994;
        break;
      case 30:
        gammascal = -4.04988034;
        break;
      case 31:
        gammascal = -3.83986237;
        break;
      case 32:
        gammascal = -2.21072317;
        break;
      case 33:
        gammascal = -2.80776694;
        break;
      case 34:
        gammascal = -2.40219736;
        break;
      case 35:
        gammascal = -2.22755851;
        break;
      case 36:
        gammascal = -1.42068154;
        break;
      case 37:
        gammascal = -3.70147701;
        break;
      case 38:
        gammascal = -8.42548707;
        break;
      case 39:
        gammascal = -9.17205775;
        break;
      case 40:
        gammascal = -0.99611531;
        break;
      case 41:
        gammascal = -6.12266771;
        break;
      case 42:
        gammascal = -3.03434318;
        break;
      case 43:
        gammascal = -2.58724432;
        break;
      case 44:
        gammascal = -2.94152428;
        break;
      case 45:
        gammascal = -2.08406140;
        break;
      case 46:
        gammascal = -3.30886859;
        break;
      case 47:
        gammascal = -2.97219928;
        break;
      case 48:
        gammascal = -3.06469315;
        break;
      case 49:
        gammascal = -3.96435075;
        break;
      case 50:
        gammascal = -1.75982006;
        break;
      case 51:
        gammascal = -2.73007709;
        break;
      case 52:
        gammascal = -2.77639176;
        break;
      case 53:
        gammascal = -3.12670692;
        break;
      case 54:
        gammascal = -0.83353793;
        break;
      case 55:
        gammascal = -3.47549649;
        break;
      case 56:
        gammascal = -5.88276127;
        break;
      case 72:
        gammascal = -3.03729539;
        break;
      case 73:
        gammascal = -4.97135726;
        break;
      case 74:
        gammascal = 4.16790013;
        break;
      case 75:
        gammascal = 0.30630646;
        break;
      case 76:
        gammascal = -1.33068824;
        break;
      case 77:
        gammascal = -2.62686662;
        break;
      case 78:
        gammascal = -1.54642144;
        break;
      case 79:
        gammascal = -1.57554917;
        break;
      case 80:
        gammascal = -2.44018189;
        break;
      case 81:
        gammascal = -2.40510708;
        break;
      case 82:
        gammascal = -1.14149533;
        break;
      case 83:
        gammascal = -1.54226592;
        break;
    }
  }
  else if (solvent == "chloroform") {
    switch (atmnr) {
      case 1:
        gammascal = -1.45124803;
        break;
      case 2:
        gammascal = -0.34288175;
        break;
      case 3:
        gammascal = -4.90550620;
        break;
      case 4:
        gammascal = -12.24883340;
        break;
      case 5:
        gammascal = -2.97051057;
        break;
      case 6:
        gammascal = -2.46046595;
        break;
      case 7:
        gammascal = -3.20499136;
        break;
      case 8:
        gammascal = -0.79186459;
        break;
      case 9:
        gammascal = 0.18544629;
        break;
      case 10:
        gammascal = -0.38955820;
        break;
      case 11:
        gammascal = -3.76356128;
        break;
      case 12:
        gammascal = -2.68971807;
        break;
      case 13:
        gammascal = -2.31501600;
        break;
      case 14:
        gammascal = 0.13149009;
        break;
      case 15:
        gammascal = -2.76183333;
        break;
      case 16:
        gammascal = -2.49678227;
        break;
      case 17:
        gammascal = -1.77429756;
        break;
      case 18:
        gammascal = -1.43838426;
        break;
      case 19:
        gammascal = -2.82408642;
        break;
      case 20:
        gammascal = -6.95175503;
        break;
      case 21:
        gammascal = -5.70468363;
        break;
      case 22:
        gammascal = -2.84574964;
        break;
      case 23:
        gammascal = -5.28512286;
        break;
      case 24:
        gammascal = -1.26239123;
        break;
      case 25:
        gammascal = -2.26516813;
        break;
      case 26:
        gammascal = -1.48532205;
        break;
      case 27:
        gammascal = 4.85394490;
        break;
      case 28:
        gammascal = -1.19463251;
        break;
      case 29:
        gammascal = -2.42708314;
        break;
      case 30:
        gammascal = -3.41088828;
        break;
      case 31:
        gammascal = -3.26901849;
        break;
      case 32:
        gammascal = -2.44120314;
        break;
      case 33:
        gammascal = -2.63812244;
        break;
      case 34:
        gammascal = -1.74626228;
        break;
      case 35:
        gammascal = -2.13180988;
        break;
      case 36:
        gammascal = -1.49256064;
        break;
      case 37:
        gammascal = -3.38516728;
        break;
      case 38:
        gammascal = -8.81329234;
        break;
      case 39:
        gammascal = -7.93257249;
        break;
      case 40:
        gammascal = 1.70913687;
        break;
      case 41:
        gammascal = -7.18032181;
        break;
      case 42:
        gammascal = -1.31323641;
        break;
      case 43:
        gammascal = -2.28930032;
        break;
      case 44:
        gammascal = -3.58478095;
        break;
      case 45:
        gammascal = -2.49391434;
        break;
      case 46:
        gammascal = -1.99369134;
        break;
      case 47:
        gammascal = -2.58226479;
        break;
      case 48:
        gammascal = 1.00000000;
        break;
      case 49:
        gammascal = -3.46355219;
        break;
      case 50:
        gammascal = -1.99680141;
        break;
      case 51:
        gammascal = -2.30696505;
        break;
      case 52:
        gammascal = -2.44521966;
        break;
      case 53:
        gammascal = -2.98470401;
        break;
      case 54:
        gammascal = -1.09297370;
        break;
      case 55:
        gammascal = -4.50182510;
        break;
      case 56:
        gammascal = -6.75767347;
        break;
      case 72:
        gammascal = -1.76894995;
        break;
      case 73:
        gammascal = -4.84348395;
        break;
      case 74:
        gammascal = 4.48381165;
        break;
      case 75:
        gammascal = -2.31673205;
        break;
      case 76:
        gammascal = -2.11643302;
        break;
      case 77:
        gammascal = -3.30136147;
        break;
      case 78:
        gammascal = -1.21261417;
        break;
      case 79:
        gammascal = -1.37980608;
        break;
      case 80:
        gammascal = -2.39587135;
        break;
      case 81:
        gammascal = -1.82204248;
        break;
      case 82:
        gammascal = -1.19017725;
        break;
      case 83:
        gammascal = -1.08346722;
        break;
    }
  }
  else if (solvent == "carbon disulfide") {
    switch (atmnr) {
      case 1:
        gammascal = -1.67630824;
        break;
      case 2:
        gammascal = -0.18644409;
        break;
      case 3:
        gammascal = -4.08142241;
        break;
      case 4:
        gammascal = -11.58201744;
        break;
      case 5:
        gammascal = -2.90130414;
        break;
      case 6:
        gammascal = -3.36892881;
        break;
      case 7:
        gammascal = -2.78632521;
        break;
      case 8:
        gammascal = -1.63610693;
        break;
      case 9:
        gammascal = -1.49297157;
        break;
      case 10:
        gammascal = -0.24699289;
        break;
      case 11:
        gammascal = -2.71615343;
        break;
      case 12:
        gammascal = -2.42746739;
        break;
      case 13:
        gammascal = -2.04324889;
        break;
      case 14:
        gammascal = 0.65718705;
        break;
      case 15:
        gammascal = -2.75409442;
        break;
      case 16:
        gammascal = -2.71058465;
        break;
      case 17:
        gammascal = -2.61131358;
        break;
      case 18:
        gammascal = -1.16644278;
        break;
      case 19:
        gammascal = -2.13314591;
        break;
      case 20:
        gammascal = -2.63632255;
        break;
      case 21:
        gammascal = -3.16877355;
        break;
      case 22:
        gammascal = -1.39796916;
        break;
      case 23:
        gammascal = -2.75518562;
        break;
      case 24:
        gammascal = -0.44225905;
        break;
      case 25:
        gammascal = -0.69011451;
        break;
      case 26:
        gammascal = -0.10776830;
        break;
      case 27:
        gammascal = -4.16810114;
        break;
      case 28:
        gammascal = -1.91968646;
        break;
      case 29:
        gammascal = -1.96185727;
        break;
      case 30:
        gammascal = -2.88161963;
        break;
      case 31:
        gammascal = -2.71223036;
        break;
      case 32:
        gammascal = -2.11892755;
        break;
      case 33:
        gammascal = -2.42274714;
        break;
      case 34:
        gammascal = -1.82161043;
        break;
      case 35:
        gammascal = -3.22845209;
        break;
      case 36:
        gammascal = -1.22681052;
        break;
      case 37:
        gammascal = -2.34829452;
        break;
      case 38:
        gammascal = -3.04588406;
        break;
      case 39:
        gammascal = -5.57003905;
        break;
      case 40:
        gammascal = 1.40101560;
        break;
      case 41:
        gammascal = -2.40146173;
        break;
      case 42:
        gammascal = -1.00476201;
        break;
      case 43:
        gammascal = 1.88708763;
        break;
      case 44:
        gammascal = -1.51182430;
        break;
      case 45:
        gammascal = -1.75652430;
        break;
      case 46:
        gammascal = -1.54489737;
        break;
      case 47:
        gammascal = -2.34681469;
        break;
      case 48:
        gammascal = -2.26824929;
        break;
      case 49:
        gammascal = -2.62671852;
        break;
      case 50:
        gammascal = -1.58167271;
        break;
      case 51:
        gammascal = -2.08854809;
        break;
      case 52:
        gammascal = -2.25711606;
        break;
      case 53:
        gammascal = -3.27762061;
        break;
      case 54:
        gammascal = -0.81986677;
        break;
      case 55:
        gammascal = -2.34291712;
        break;
      case 56:
        gammascal = -2.98284235;
        break;
      case 72:
        gammascal = -1.80591619;
        break;
      case 73:
        gammascal = -1.84572202;
        break;
      case 74:
        gammascal = 2.39890896;
        break;
      case 75:
        gammascal = -0.17262619;
        break;
      case 76:
        gammascal = -0.91545732;
        break;
      case 77:
        gammascal = -2.33520652;
        break;
      case 78:
        gammascal = -1.13222036;
        break;
      case 79:
        gammascal = -1.29972175;
        break;
      case 80:
        gammascal = -1.31761075;
        break;
      case 81:
        gammascal = -1.53817638;
        break;
      case 82:
        gammascal = -0.89984085;
        break;
      case 83:
        gammascal = -0.85591904;
        break;
    }
  }
  else if (solvent == "dioxane") {
    switch (atmnr) {
      case 1:
        gammascal = -1.30290774;
        break;
      case 2:
        gammascal = -0.11747889;
        break;
      case 3:
        gammascal = -5.21763288;
        break;
      case 4:
        gammascal = -12.06254155;
        break;
      case 5:
        gammascal = -2.26649580;
        break;
      case 6:
        gammascal = -2.32061003;
        break;
      case 7:
        gammascal = -2.10913167;
        break;
      case 8:
        gammascal = 0.41984122;
        break;
      case 9:
        gammascal = 0.11160374;
        break;
      case 10:
        gammascal = -0.16720505;
        break;
      case 11:
        gammascal = -3.15925017;
        break;
      case 12:
        gammascal = -2.35890206;
        break;
      case 13:
        gammascal = -1.78234075;
        break;
      case 14:
        gammascal = 1.00019486;
        break;
      case 15:
        gammascal = -1.50992486;
        break;
      case 16:
        gammascal = -1.55783018;
        break;
      case 17:
        gammascal = -1.41889343;
        break;
      case 18:
        gammascal = -0.90873964;
        break;
      case 19:
        gammascal = -2.45741294;
        break;
      case 20:
        gammascal = -5.14009725;
        break;
      case 21:
        gammascal = -3.82643718;
        break;
      case 22:
        gammascal = -1.84498043;
        break;
      case 23:
        gammascal = -5.65375333;
        break;
      case 24:
        gammascal = -1.93935987;
        break;
      case 25:
        gammascal = -1.34357575;
        break;
      case 26:
        gammascal = 1.26875485;
        break;
      case 27:
        gammascal = -10.07963841;
        break;
      case 28:
        gammascal = -1.37816150;
        break;
      case 29:
        gammascal = -1.64650261;
        break;
      case 30:
        gammascal = -3.08213950;
        break;
      case 31:
        gammascal = -2.70966327;
        break;
      case 32:
        gammascal = -1.47471128;
        break;
      case 33:
        gammascal = -1.95292540;
        break;
      case 34:
        gammascal = -1.11231973;
        break;
      case 35:
        gammascal = -1.78670355;
        break;
      case 36:
        gammascal = -0.91293395;
        break;
      case 37:
        gammascal = -2.65076547;
        break;
      case 38:
        gammascal = -5.98063844;
        break;
      case 39:
        gammascal = -7.79249444;
        break;
      case 40:
        gammascal = 3.08009398;
        break;
      case 41:
        gammascal = -4.64823076;
        break;
      case 42:
        gammascal = -0.37879616;
        break;
      case 43:
        gammascal = 0.04905680;
        break;
      case 44:
        gammascal = -3.37456181;
        break;
      case 45:
        gammascal = -2.38573910;
        break;
      case 46:
        gammascal = -1.97354440;
        break;
      case 47:
        gammascal = -2.08001446;
        break;
      case 48:
        gammascal = -1.63445987;
        break;
      case 49:
        gammascal = -2.78036246;
        break;
      case 50:
        gammascal = -1.37470577;
        break;
      case 51:
        gammascal = -1.44947390;
        break;
      case 52:
        gammascal = -1.90241170;
        break;
      case 53:
        gammascal = -2.14102405;
        break;
      case 54:
        gammascal = -0.95822732;
        break;
      case 55:
        gammascal = -3.09207921;
        break;
      case 56:
        gammascal = -5.17995903;
        break;
      case 72:
        gammascal = 0.40509807;
        break;
      case 73:
        gammascal = -2.74316309;
        break;
      case 74:
        gammascal = 5.47707232;
        break;
      case 75:
        gammascal = -3.93096056;
        break;
      case 76:
        gammascal = -0.96010077;
        break;
      case 77:
        gammascal = -1.99348798;
        break;
      case 78:
        gammascal = -0.65434419;
        break;
      case 79:
        gammascal = -0.92143731;
        break;
      case 80:
        gammascal = -2.71423264;
        break;
      case 81:
        gammascal = -1.27668848;
        break;
      case 82:
        gammascal = -0.82650501;
        break;
      case 83:
        gammascal = -0.31800653;
        break;
    }
  }
  else if (solvent == "dmf") {
    switch (atmnr) {
      case 1:
        gammascal = -1.48208914;
        break;
      case 2:
        gammascal = -0.28924933;
        break;
      case 3:
        gammascal = -6.08890005;
        break;
      case 4:
        gammascal = -13.56322009;
        break;
      case 5:
        gammascal = -2.36024155;
        break;
      case 6:
        gammascal = -2.51386998;
        break;
      case 7:
        gammascal = -3.06520025;
        break;
      case 8:
        gammascal = -2.46380756;
        break;
      case 9:
        gammascal = -0.91151625;
        break;
      case 10:
        gammascal = -0.32645958;
        break;
      case 11:
        gammascal = -3.89741810;
        break;
      case 12:
        gammascal = -2.69211933;
        break;
      case 13:
        gammascal = -1.97835397;
        break;
      case 14:
        gammascal = 0.20458102;
        break;
      case 15:
        gammascal = -2.03561547;
        break;
      case 16:
        gammascal = -2.18873123;
        break;
      case 17:
        gammascal = -1.62204391;
        break;
      case 18:
        gammascal = -1.25679804;
        break;
      case 19:
        gammascal = -3.49876406;
        break;
      case 20:
        gammascal = -6.29423117;
        break;
      case 21:
        gammascal = -5.81709631;
        break;
      case 22:
        gammascal = -4.47956803;
        break;
      case 23:
        gammascal = -4.64094542;
        break;
      case 24:
        gammascal = -3.65813236;
        break;
      case 25:
        gammascal = -1.89789079;
        break;
      case 26:
        gammascal = 2.41669184;
        break;
      case 27:
        gammascal = 13.35689270;
        break;
      case 28:
        gammascal = -1.79029475;
        break;
      case 29:
        gammascal = -1.92742742;
        break;
      case 30:
        gammascal = -3.66031801;
        break;
      case 31:
        gammascal = -2.72610721;
        break;
      case 32:
        gammascal = -1.79015528;
        break;
      case 33:
        gammascal = -2.21054855;
        break;
      case 34:
        gammascal = -1.73401472;
        break;
      case 35:
        gammascal = -2.27602615;
        break;
      case 36:
        gammascal = -1.26418738;
        break;
      case 37:
        gammascal = -3.21304702;
        break;
      case 38:
        gammascal = -7.22048375;
        break;
      case 39:
        gammascal = -8.15841847;
        break;
      case 40:
        gammascal = -0.87578501;
        break;
      case 41:
        gammascal = -4.16517801;
        break;
      case 42:
        gammascal = -4.36451780;
        break;
      case 43:
        gammascal = 2.69328761;
        break;
      case 44:
        gammascal = -2.85560707;
        break;
      case 45:
        gammascal = -1.26368814;
        break;
      case 46:
        gammascal = -2.04573691;
        break;
      case 47:
        gammascal = -2.23767811;
        break;
      case 48:
        gammascal = -3.06188154;
        break;
      case 49:
        gammascal = -3.10002805;
        break;
      case 50:
        gammascal = -1.54588935;
        break;
      case 51:
        gammascal = -1.85028119;
        break;
      case 52:
        gammascal = -2.08450522;
        break;
      case 53:
        gammascal = -2.87721546;
        break;
      case 54:
        gammascal = -1.18600787;
        break;
      case 55:
        gammascal = -3.30211598;
        break;
      case 56:
        gammascal = -6.03310028;
        break;
      case 72:
        gammascal = -2.90152918;
        break;
      case 73:
        gammascal = -3.32546395;
        break;
      case 74:
        gammascal = 2.93771818;
        break;
      case 75:
        gammascal = 0.85423436;
        break;
      case 76:
        gammascal = -0.69768791;
        break;
      case 77:
        gammascal = -1.77911712;
        break;
      case 78:
        gammascal = -0.92456925;
        break;
      case 79:
        gammascal = -1.24218444;
        break;
      case 80:
        gammascal = -2.73624618;
        break;
      case 81:
        gammascal = -2.11592322;
        break;
      case 82:
        gammascal = -1.30882885;
        break;
      case 83:
        gammascal = -1.40473043;
        break;
    }
  }
  else if (solvent == "dmso") {
    switch (atmnr) {
      case 1:
        gammascal = -1.38664253;
        break;
      case 2:
        gammascal = -0.23793625;
        break;
      case 3:
        gammascal = -6.03843971;
        break;
      case 4:
        gammascal = -13.22919070;
        break;
      case 5:
        gammascal = -2.22701589;
        break;
      case 6:
        gammascal = -2.38738228;
        break;
      case 7:
        gammascal = -2.90820463;
        break;
      case 8:
        gammascal = -2.31874520;
        break;
      case 9:
        gammascal = -0.99639634;
        break;
      case 10:
        gammascal = -0.26161245;
        break;
      case 11:
        gammascal = -3.82454484;
        break;
      case 12:
        gammascal = -2.64362668;
        break;
      case 13:
        gammascal = -1.86876630;
        break;
      case 14:
        gammascal = 0.19339626;
        break;
      case 15:
        gammascal = -1.87968108;
        break;
      case 16:
        gammascal = -2.02769642;
        break;
      case 17:
        gammascal = -1.62247140;
        break;
      case 18:
        gammascal = -1.24077419;
        break;
      case 19:
        gammascal = -3.03583284;
        break;
      case 20:
        gammascal = -6.49296876;
        break;
      case 21:
        gammascal = -6.09881747;
        break;
      case 22:
        gammascal = -4.58817578;
        break;
      case 23:
        gammascal = -4.82700240;
        break;
      case 24:
        gammascal = -4.20180767;
        break;
      case 25:
        gammascal = -1.98834227;
        break;
      case 26:
        gammascal = 2.87126034;
        break;
      case 27:
        gammascal = 16.64033536;
        break;
      case 28:
        gammascal = -1.97960055;
        break;
      case 29:
        gammascal = -1.81744327;
        break;
      case 30:
        gammascal = -3.58357716;
        break;
      case 31:
        gammascal = -2.65452713;
        break;
      case 32:
        gammascal = -1.73348397;
        break;
      case 33:
        gammascal = -2.12848796;
        break;
      case 34:
        gammascal = -1.67987360;
        break;
      case 35:
        gammascal = -2.33227071;
        break;
      case 36:
        gammascal = -1.19822635;
        break;
      case 37:
        gammascal = -3.15766157;
        break;
      case 38:
        gammascal = -7.38729953;
        break;
      case 39:
        gammascal = -8.04591955;
        break;
      case 40:
        gammascal = -0.68974293;
        break;
      case 41:
        gammascal = -4.34602557;
        break;
      case 42:
        gammascal = -4.57828188;
        break;
      case 43:
        gammascal = 2.79315217;
        break;
      case 44:
        gammascal = -3.00760331;
        break;
      case 45:
        gammascal = -1.44232566;
        break;
      case 46:
        gammascal = -1.91915984;
        break;
      case 47:
        gammascal = -2.11409085;
        break;
      case 48:
        gammascal = -2.80634213;
        break;
      case 49:
        gammascal = -3.03358640;
        break;
      case 50:
        gammascal = -1.49292716;
        break;
      case 51:
        gammascal = -1.73235161;
        break;
      case 52:
        gammascal = -1.97756737;
        break;
      case 53:
        gammascal = -2.73254858;
        break;
      case 54:
        gammascal = -1.15534479;
        break;
      case 55:
        gammascal = -2.97559221;
        break;
      case 56:
        gammascal = -6.39036337;
        break;
      case 72:
        gammascal = -2.96470778;
        break;
      case 73:
        gammascal = -3.31841657;
        break;
      case 74:
        gammascal = 3.53259054;
        break;
      case 75:
        gammascal = 0.43466852;
        break;
      case 76:
        gammascal = -0.63758586;
        break;
      case 77:
        gammascal = -1.93676931;
        break;
      case 78:
        gammascal = -0.88498933;
        break;
      case 79:
        gammascal = -1.17566028;
        break;
      case 80:
        gammascal = -2.52788484;
        break;
      case 81:
        gammascal = -2.02962235;
        break;
      case 82:
        gammascal = -0.74677606;
        break;
      case 83:
        gammascal = -1.23774498;
        break;
    }
  }
  else if (solvent == "ethanol") {
    switch (atmnr) {
      case 1:
        gammascal = -1.25218220;
        break;
      case 2:
        gammascal = -0.13897627;
        break;
      case 3:
        gammascal = -0.86743826;
        break;
      case 4:
        gammascal = -1.81320472;
        break;
      case 5:
        gammascal = -2.27731021;
        break;
      case 6:
        gammascal = -2.39031618;
        break;
      case 7:
        gammascal = -2.47198050;
        break;
      case 8:
        gammascal = -0.08505090;
        break;
      case 9:
        gammascal = 0.08097631;
        break;
      case 10:
        gammascal = -0.17916016;
        break;
      case 11:
        gammascal = -7.72113272;
        break;
      case 12:
        gammascal = -1.03737823;
        break;
      case 13:
        gammascal = -1.88138258;
        break;
      case 14:
        gammascal = 0.51149221;
        break;
      case 15:
        gammascal = -1.84167649;
        break;
      case 16:
        gammascal = -1.69336780;
        break;
      case 17:
        gammascal = -1.26036493;
        break;
      case 18:
        gammascal = -0.65422551;
        break;
      case 19:
        gammascal = -1.58921466;
        break;
      case 20:
        gammascal = -9.11059797;
        break;
      case 21:
        gammascal = -7.22982390;
        break;
      case 22:
        gammascal = -4.61868820;
        break;
      case 23:
        gammascal = -4.57645783;
        break;
      case 24:
        gammascal = -2.67385778;
        break;
      case 25:
        gammascal = -2.93711233;
        break;
      case 26:
        gammascal = 2.01257589;
        break;
      case 27:
        gammascal = 9.74833738;
        break;
      case 28:
        gammascal = -1.24821122;
        break;
      case 29:
        gammascal = -1.74669202;
        break;
      case 30:
        gammascal = -0.33760131;
        break;
      case 31:
        gammascal = -2.78116024;
        break;
      case 32:
        gammascal = -1.74297664;
        break;
      case 33:
        gammascal = -2.01121412;
        break;
      case 34:
        gammascal = -1.39066829;
        break;
      case 35:
        gammascal = -1.61519239;
        break;
      case 36:
        gammascal = -1.26500230;
        break;
      case 37:
        gammascal = -4.62202053;
        break;
      case 38:
        gammascal = 11.17395333;
        break;
      case 39:
        gammascal = -8.32584947;
        break;
      case 40:
        gammascal = 0.87582141;
        break;
      case 41:
        gammascal = -5.16854115;
        break;
      case 42:
        gammascal = -2.71195928;
        break;
      case 43:
        gammascal = -2.90742248;
        break;
      case 44:
        gammascal = -3.41807262;
        break;
      case 45:
        gammascal = -1.34857246;
        break;
      case 46:
        gammascal = -1.00979441;
        break;
      case 47:
        gammascal = -2.06679863;
        break;
      case 48:
        gammascal = -1.97683363;
        break;
      case 49:
        gammascal = -3.18881568;
        break;
      case 50:
        gammascal = -1.54608577;
        break;
      case 51:
        gammascal = -1.87562748;
        break;
      case 52:
        gammascal = -2.09978128;
        break;
      case 53:
        gammascal = -2.25403543;
        break;
      case 54:
        gammascal = -0.98014403;
        break;
      case 55:
        gammascal = -1.99502357;
        break;
      case 56:
        gammascal = -5.98891556;
        break;
      case 72:
        gammascal = -2.77351612;
        break;
      case 73:
        gammascal = -4.04830880;
        break;
      case 74:
        gammascal = 4.11271145;
        break;
      case 75:
        gammascal = -3.49989324;
        break;
      case 76:
        gammascal = -1.43077966;
        break;
      case 77:
        gammascal = -2.13182301;
        break;
      case 78:
        gammascal = -0.72044902;
        break;
      case 79:
        gammascal = -0.93399195;
        break;
      case 80:
        gammascal = -2.51417568;
        break;
      case 81:
        gammascal = -1.84648806;
        break;
      case 82:
        gammascal = -0.95977914;
        break;
      case 83:
        gammascal = -0.55418129;
        break;
    }
  }
  else if (solvent == "diethyl ether") {
    switch (atmnr) {
      case 1:
        gammascal = -1.37103872;
        break;
      case 2:
        gammascal = -0.41977745;
        break;
      case 3:
        gammascal = -3.69586174;
        break;
      case 4:
        gammascal = -11.61278609;
        break;
      case 5:
        gammascal = -2.40421376;
        break;
      case 6:
        gammascal = -2.42446747;
        break;
      case 7:
        gammascal = -2.00586580;
        break;
      case 8:
        gammascal = -0.87351482;
        break;
      case 9:
        gammascal = -0.42023641;
        break;
      case 10:
        gammascal = -0.46118972;
        break;
      case 11:
        gammascal = -2.98655641;
        break;
      case 12:
        gammascal = -2.44766523;
        break;
      case 13:
        gammascal = -1.72431597;
        break;
      case 14:
        gammascal = 0.22255920;
        break;
      case 15:
        gammascal = -2.24449155;
        break;
      case 16:
        gammascal = -2.06128996;
        break;
      case 17:
        gammascal = -1.61171641;
        break;
      case 18:
        gammascal = -1.41599498;
        break;
      case 19:
        gammascal = -3.24984604;
        break;
      case 20:
        gammascal = -4.21734547;
        break;
      case 21:
        gammascal = -4.59486463;
        break;
      case 22:
        gammascal = -3.61337416;
        break;
      case 23:
        gammascal = -4.20015559;
        break;
      case 24:
        gammascal = -3.44971035;
        break;
      case 25:
        gammascal = -2.53825202;
        break;
      case 26:
        gammascal = 2.54721106;
        break;
      case 27:
        gammascal = -1.74362570;
        break;
      case 28:
        gammascal = -1.53793503;
        break;
      case 29:
        gammascal = -2.00296184;
        break;
      case 30:
        gammascal = -3.02669579;
        break;
      case 31:
        gammascal = -2.56020451;
        break;
      case 32:
        gammascal = -2.12420203;
        break;
      case 33:
        gammascal = -2.20339162;
        break;
      case 34:
        gammascal = -1.67875918;
        break;
      case 35:
        gammascal = -2.12984787;
        break;
      case 36:
        gammascal = -1.40156757;
        break;
      case 37:
        gammascal = -2.89780651;
        break;
      case 38:
        gammascal = -5.31327909;
        break;
      case 39:
        gammascal = -6.57325469;
        break;
      case 40:
        gammascal = 0.33434000;
        break;
      case 41:
        gammascal = -3.84490382;
        break;
      case 42:
        gammascal = -2.74384260;
        break;
      case 43:
        gammascal = -2.79468289;
        break;
      case 44:
        gammascal = -2.81191075;
        break;
      case 45:
        gammascal = -1.66841518;
        break;
      case 46:
        gammascal = -1.90612170;
        break;
      case 47:
        gammascal = -2.19547580;
        break;
      case 48:
        gammascal = -2.05732619;
        break;
      case 49:
        gammascal = -2.72176429;
        break;
      case 50:
        gammascal = -1.80202693;
        break;
      case 51:
        gammascal = -1.92324306;
        break;
      case 52:
        gammascal = -2.09793839;
        break;
      case 53:
        gammascal = -2.63669739;
        break;
      case 54:
        gammascal = -1.23868974;
        break;
      case 55:
        gammascal = -2.31036608;
        break;
      case 56:
        gammascal = -4.36607834;
        break;
      case 72:
        gammascal = -2.03725494;
        break;
      case 73:
        gammascal = -2.67608710;
        break;
      case 74:
        gammascal = 3.20814838;
        break;
      case 75:
        gammascal = -2.65692482;
        break;
      case 76:
        gammascal = -1.22214666;
        break;
      case 77:
        gammascal = -2.15626456;
        break;
      case 78:
        gammascal = -1.05683034;
        break;
      case 79:
        gammascal = -1.23143100;
        break;
      case 80:
        gammascal = -2.32551963;
        break;
      case 81:
        gammascal = -1.66013720;
        break;
      case 82:
        gammascal = -1.04271426;
        break;
      case 83:
        gammascal = -1.14802293;
        break;
    }
  }
  else if (solvent == "ethyl acetate") {
    switch (atmnr) {
      case 1:
        gammascal = -1.47464227;
        break;
      case 2:
        gammascal = -0.10752409;
        break;
      case 3:
        gammascal = -6.76065469;
        break;
      case 4:
        gammascal = -18.22167040;
        break;
      case 5:
        gammascal = -2.92965152;
        break;
      case 6:
        gammascal = -2.56430588;
        break;
      case 7:
        gammascal = -3.40983221;
        break;
      case 8:
        gammascal = -1.91105293;
        break;
      case 9:
        gammascal = -0.01270725;
        break;
      case 10:
        gammascal = -0.13755626;
        break;
      case 11:
        gammascal = -4.92259034;
        break;
      case 12:
        gammascal = -3.60583392;
        break;
      case 13:
        gammascal = -2.01072681;
        break;
      case 14:
        gammascal = 0.23949435;
        break;
      case 15:
        gammascal = -2.35667752;
        break;
      case 16:
        gammascal = -2.40564101;
        break;
      case 17:
        gammascal = -1.92060622;
        break;
      case 18:
        gammascal = -1.01258184;
        break;
      case 19:
        gammascal = -3.26559937;
        break;
      case 20:
        gammascal = -6.45262589;
        break;
      case 21:
        gammascal = -5.52115691;
        break;
      case 22:
        gammascal = -3.59313979;
        break;
      case 23:
        gammascal = -6.59576807;
        break;
      case 24:
        gammascal = -3.13763691;
        break;
      case 25:
        gammascal = -0.88670591;
        break;
      case 26:
        gammascal = 1.74780192;
        break;
      case 27:
        gammascal = 9.88340512;
        break;
      case 28:
        gammascal = -0.88297167;
        break;
      case 29:
        gammascal = -2.44026373;
        break;
      case 30:
        gammascal = -4.57961109;
        break;
      case 31:
        gammascal = -3.24252108;
        break;
      case 32:
        gammascal = -1.88416436;
        break;
      case 33:
        gammascal = -2.67771231;
        break;
      case 34:
        gammascal = -1.81379923;
        break;
      case 35:
        gammascal = -1.96620377;
        break;
      case 36:
        gammascal = -1.07379349;
        break;
      case 37:
        gammascal = -3.22324039;
        break;
      case 38:
        gammascal = -8.16078288;
        break;
      case 39:
        gammascal = -9.67527276;
        break;
      case 40:
        gammascal = -1.97565069;
        break;
      case 41:
        gammascal = -7.30870212;
        break;
      case 42:
        gammascal = -2.99011992;
        break;
      case 43:
        gammascal = -0.98908968;
        break;
      case 44:
        gammascal = -1.86634080;
        break;
      case 45:
        gammascal = -1.51427958;
        break;
      case 46:
        gammascal = -2.72024442;
        break;
      case 47:
        gammascal = -2.63735961;
        break;
      case 48:
        gammascal = -3.20722817;
        break;
      case 49:
        gammascal = -3.50287298;
        break;
      case 50:
        gammascal = -1.61913490;
        break;
      case 51:
        gammascal = -2.39438994;
        break;
      case 52:
        gammascal = -2.35173908;
        break;
      case 53:
        gammascal = -2.43955645;
        break;
      case 54:
        gammascal = -0.73066180;
        break;
      case 55:
        gammascal = -3.41083482;
        break;
      case 56:
        gammascal = -6.53195023;
        break;
      case 72:
        gammascal = -1.29112956;
        break;
      case 73:
        gammascal = -5.25328160;
        break;
      case 74:
        gammascal = 2.71032461;
        break;
      case 75:
        gammascal = 2.41078006;
        break;
      case 76:
        gammascal = -0.69666985;
        break;
      case 77:
        gammascal = -2.53087586;
        break;
      case 78:
        gammascal = -1.11903365;
        break;
      case 79:
        gammascal = -1.30417534;
        break;
      case 80:
        gammascal = -2.65071281;
        break;
      case 81:
        gammascal = -2.26744852;
        break;
      case 82:
        gammascal = -1.01854238;
        break;
      case 83:
        gammascal = -1.52083385;
        break;
    }
  }
  else if (solvent == "furane") {
    switch (atmnr) {
      case 1:
        gammascal = -1.34331419;
        break;
      case 2:
        gammascal = -0.10528890;
        break;
      case 3:
        gammascal = -5.44529786;
        break;
      case 4:
        gammascal = -15.64566091;
        break;
      case 5:
        gammascal = -2.89065815;
        break;
      case 6:
        gammascal = -2.76752119;
        break;
      case 7:
        gammascal = -2.95774936;
        break;
      case 8:
        gammascal = -0.18194139;
        break;
      case 9:
        gammascal = 0.38260590;
        break;
      case 10:
        gammascal = -0.13731173;
        break;
      case 11:
        gammascal = -4.40895349;
        break;
      case 12:
        gammascal = -2.96664573;
        break;
      case 13:
        gammascal = -1.81939300;
        break;
      case 14:
        gammascal = 0.53703682;
        break;
      case 15:
        gammascal = -2.46510797;
        break;
      case 16:
        gammascal = -2.34661217;
        break;
      case 17:
        gammascal = -1.68263052;
        break;
      case 18:
        gammascal = -0.94180070;
        break;
      case 19:
        gammascal = -2.67597579;
        break;
      case 20:
        gammascal = -5.96743843;
        break;
      case 21:
        gammascal = -4.53958647;
        break;
      case 22:
        gammascal = -2.37519332;
        break;
      case 23:
        gammascal = -6.48695622;
        break;
      case 24:
        gammascal = -0.95637608;
        break;
      case 25:
        gammascal = -2.33324934;
        break;
      case 26:
        gammascal = -2.11005655;
        break;
      case 27:
        gammascal = -2.47431927;
        break;
      case 28:
        gammascal = -0.70122124;
        break;
      case 29:
        gammascal = -2.32818118;
        break;
      case 30:
        gammascal = -3.63711001;
        break;
      case 31:
        gammascal = -3.18443664;
        break;
      case 32:
        gammascal = -1.82164033;
        break;
      case 33:
        gammascal = -2.55821661;
        break;
      case 34:
        gammascal = -1.49033789;
        break;
      case 35:
        gammascal = -2.04893421;
        break;
      case 36:
        gammascal = -1.01036070;
        break;
      case 37:
        gammascal = -3.03695448;
        break;
      case 38:
        gammascal = -7.25156713;
        break;
      case 39:
        gammascal = -8.65094836;
        break;
      case 40:
        gammascal = 0.90539857;
        break;
      case 41:
        gammascal = -6.82658746;
        break;
      case 42:
        gammascal = -1.04614201;
        break;
      case 43:
        gammascal = -2.12328069;
        break;
      case 44:
        gammascal = -3.12626442;
        break;
      case 45:
        gammascal = -3.89828505;
        break;
      case 46:
        gammascal = -1.83026210;
        break;
      case 47:
        gammascal = -2.52136489;
        break;
      case 48:
        gammascal = -2.35134609;
        break;
      case 49:
        gammascal = -3.50011757;
        break;
      case 50:
        gammascal = -1.57550664;
        break;
      case 51:
        gammascal = -2.38268715;
        break;
      case 52:
        gammascal = -2.34460162;
        break;
      case 53:
        gammascal = -2.47225415;
        break;
      case 54:
        gammascal = -0.73527495;
        break;
      case 55:
        gammascal = -3.02936007;
        break;
      case 56:
        gammascal = -5.70415926;
        break;
      case 72:
        gammascal = -0.41590025;
        break;
      case 73:
        gammascal = -4.60920960;
        break;
      case 74:
        gammascal = 1.37287726;
        break;
      case 75:
        gammascal = -3.03161593;
        break;
      case 76:
        gammascal = -2.77466931;
        break;
      case 77:
        gammascal = -3.98630896;
        break;
      case 78:
        gammascal = -1.12315108;
        break;
      case 79:
        gammascal = -1.14369672;
        break;
      case 80:
        gammascal = -2.43237503;
        break;
      case 81:
        gammascal = -1.70071275;
        break;
      case 82:
        gammascal = -1.00202900;
        break;
      case 83:
        gammascal = -0.58455271;
        break;
    }
  }
  else if (solvent == "hexadecane") {
    switch (atmnr) {
      case 1:
        gammascal = -1.41302846;
        break;
      case 2:
        gammascal = -0.21817594;
        break;
      case 3:
        gammascal = -3.96046413;
        break;
      case 4:
        gammascal = -6.37565919;
        break;
      case 5:
        gammascal = -4.99798651;
        break;
      case 6:
        gammascal = -2.95180128;
        break;
      case 7:
        gammascal = -2.19535198;
        break;
      case 8:
        gammascal = -0.40316839;
        break;
      case 9:
        gammascal = -0.34216736;
        break;
      case 10:
        gammascal = -0.22904045;
        break;
      case 11:
        gammascal = -4.26653719;
        break;
      case 12:
        gammascal = 1.00000000;
        break;
      case 13:
        gammascal = -1.53810072;
        break;
      case 14:
        gammascal = 0.53791645;
        break;
      case 15:
        gammascal = 0.31766295;
        break;
      case 16:
        gammascal = -2.65615380;
        break;
      case 17:
        gammascal = -1.99636321;
        break;
      case 18:
        gammascal = -1.32626136;
        break;
      case 19:
        gammascal = -1.59213275;
        break;
      case 20:
        gammascal = -1.78922916;
        break;
      case 21:
        gammascal = -2.07782586;
        break;
      case 22:
        gammascal = -0.54299058;
        break;
      case 23:
        gammascal = -4.09491593;
        break;
      case 24:
        gammascal = -4.41626389;
        break;
      case 25:
        gammascal = 1.91772631;
        break;
      case 26:
        gammascal = -1.79675624;
        break;
      case 27:
        gammascal = 1.61743962;
        break;
      case 28:
        gammascal = 5.22134905;
        break;
      case 29:
        gammascal = -2.00488443;
        break;
      case 30:
        gammascal = -3.03173384;
        break;
      case 31:
        gammascal = -2.82289972;
        break;
      case 32:
        gammascal = -2.31282947;
        break;
      case 33:
        gammascal = -2.02622711;
        break;
      case 34:
        gammascal = -1.87815147;
        break;
      case 35:
        gammascal = -1.92922352;
        break;
      case 36:
        gammascal = -1.36093524;
        break;
      case 37:
        gammascal = -2.27884213;
        break;
      case 38:
        gammascal = -0.92152770;
        break;
      case 39:
        gammascal = -7.04205120;
        break;
      case 40:
        gammascal = -1.19076792;
        break;
      case 41:
        gammascal = -7.85868436;
        break;
      case 42:
        gammascal = -3.15491305;
        break;
      case 43:
        gammascal = 9.53387271;
        break;
      case 44:
        gammascal = 21.52049825;
        break;
      case 45:
        gammascal = -6.08919614;
        break;
      case 46:
        gammascal = 1.90117938;
        break;
      case 47:
        gammascal = -2.65656836;
        break;
      case 48:
        gammascal = -2.47211880;
        break;
      case 49:
        gammascal = -2.92293593;
        break;
      case 50:
        gammascal = -1.82469302;
        break;
      case 51:
        gammascal = -2.20543823;
        break;
      case 52:
        gammascal = -2.44141473;
        break;
      case 53:
        gammascal = -2.39889034;
        break;
      case 54:
        gammascal = -0.80676041;
        break;
      case 55:
        gammascal = -1.86448685;
        break;
      case 56:
        gammascal = -1.20279093;
        break;
      case 72:
        gammascal = -2.22589358;
        break;
      case 73:
        gammascal = -2.87188948;
        break;
      case 74:
        gammascal = -6.51431543;
        break;
      case 75:
        gammascal = -2.24985596;
        break;
      case 76:
        gammascal = 3.27655073;
        break;
      case 77:
        gammascal = -2.47752678;
        break;
      case 78:
        gammascal = -0.12881989;
        break;
      case 79:
        gammascal = -1.95683437;
        break;
      case 80:
        gammascal = -1.64658014;
        break;
      case 81:
        gammascal = -1.74249857;
        break;
      case 82:
        gammascal = -1.04977272;
        break;
      case 83:
        gammascal = -0.85285748;
        break;
    }
  }
  else if (solvent == "hexane") {
    switch (atmnr) {
      case 1:
        gammascal = -1.52102878;
        break;
      case 2:
        gammascal = -0.29572459;
        break;
      case 3:
        gammascal = -3.22402901;
        break;
      case 4:
        gammascal = -9.28422093;
        break;
      case 5:
        gammascal = -2.92454598;
        break;
      case 6:
        gammascal = -2.61044437;
        break;
      case 7:
        gammascal = -1.70103214;
        break;
      case 8:
        gammascal = -0.07571818;
        break;
      case 9:
        gammascal = -0.47601293;
        break;
      case 10:
        gammascal = -0.32203048;
        break;
      case 11:
        gammascal = -2.30502588;
        break;
      case 12:
        gammascal = -2.23312483;
        break;
      case 13:
        gammascal = -1.90407247;
        break;
      case 14:
        gammascal = 0.33954238;
        break;
      case 15:
        gammascal = -2.84448645;
        break;
      case 16:
        gammascal = -2.16182302;
        break;
      case 17:
        gammascal = -2.01000174;
        break;
      case 18:
        gammascal = -1.31362932;
        break;
      case 19:
        gammascal = -2.14681519;
        break;
      case 20:
        gammascal = -1.80787190;
        break;
      case 21:
        gammascal = -2.19811180;
        break;
      case 22:
        gammascal = -1.03665074;
        break;
      case 23:
        gammascal = -4.22573924;
        break;
      case 24:
        gammascal = -0.66982776;
        break;
      case 25:
        gammascal = -1.66668578;
        break;
      case 26:
        gammascal = -0.33282416;
        break;
      case 27:
        gammascal = -17.21674320;
        break;
      case 28:
        gammascal = -2.08441258;
        break;
      case 29:
        gammascal = -2.14713230;
        break;
      case 30:
        gammascal = -2.71400958;
        break;
      case 31:
        gammascal = -2.74969834;
        break;
      case 32:
        gammascal = -2.44700645;
        break;
      case 33:
        gammascal = -2.49105589;
        break;
      case 34:
        gammascal = -1.67123353;
        break;
      case 35:
        gammascal = -1.92922352;
        break;
      case 36:
        gammascal = -1.36228565;
        break;
      case 37:
        gammascal = -2.11661217;
        break;
      case 38:
        gammascal = -1.85566928;
        break;
      case 39:
        gammascal = -5.47115227;
        break;
      case 40:
        gammascal = 3.31129723;
        break;
      case 41:
        gammascal = -3.22204731;
        break;
      case 42:
        gammascal = 0.37929174;
        break;
      case 43:
        gammascal = -0.12147642;
        break;
      case 44:
        gammascal = -2.95656362;
        break;
      case 45:
        gammascal = -3.76018325;
        break;
      case 46:
        gammascal = -2.22228769;
        break;
      case 47:
        gammascal = -2.51103975;
        break;
      case 48:
        gammascal = -2.20345787;
        break;
      case 49:
        gammascal = -2.53170602;
        break;
      case 50:
        gammascal = -1.97666701;
        break;
      case 51:
        gammascal = -2.12863417;
        break;
      case 52:
        gammascal = -2.33015886;
        break;
      case 53:
        gammascal = -2.39889034;
        break;
      case 54:
        gammascal = -0.96420681;
        break;
      case 55:
        gammascal = -1.77296385;
        break;
      case 56:
        gammascal = -1.96187875;
        break;
      case 72:
        gammascal = 1.40973485;
        break;
      case 73:
        gammascal = -1.99055268;
        break;
      case 74:
        gammascal = 1.95951904;
        break;
      case 75:
        gammascal = -3.91096637;
        break;
      case 76:
        gammascal = -1.69948253;
        break;
      case 77:
        gammascal = -3.33070532;
        break;
      case 78:
        gammascal = -1.22150714;
        break;
      case 79:
        gammascal = -1.41239427;
        break;
      case 80:
        gammascal = -1.83343380;
        break;
      case 81:
        gammascal = -1.40576541;
        break;
      case 82:
        gammascal = -1.12763505;
        break;
      case 83:
        gammascal = -0.75529202;
        break;
    }
  }
  else if (solvent == "methanol") {
    switch (atmnr) {
      case 1:
        gammascal = -1.15883854;
        break;
      case 2:
        gammascal = -0.08489676;
        break;
      case 3:
        gammascal = -6.19150578;
        break;
      case 4:
        gammascal = -17.34818113;
        break;
      case 5:
        gammascal = -2.16450184;
        break;
      case 6:
        gammascal = -2.41928336;
        break;
      case 7:
        gammascal = -2.09250314;
        break;
      case 8:
        gammascal = 0.08009064;
        break;
      case 9:
        gammascal = 0.09231303;
        break;
      case 10:
        gammascal = -0.10961343;
        break;
      case 11:
        gammascal = -3.71341215;
        break;
      case 12:
        gammascal = -1.03737823;
        break;
      case 13:
        gammascal = -1.84333839;
        break;
      case 14:
        gammascal = 0.32165527;
        break;
      case 15:
        gammascal = -1.39958467;
        break;
      case 16:
        gammascal = -1.59264452;
        break;
      case 17:
        gammascal = -1.16811046;
        break;
      case 18:
        gammascal = -0.99317649;
        break;
      case 19:
        gammascal = -3.31190335;
        break;
      case 20:
        gammascal = -9.07931008;
        break;
      case 21:
        gammascal = -7.46009442;
        break;
      case 22:
        gammascal = -4.68809582;
        break;
      case 23:
        gammascal = -5.33324810;
        break;
      case 24:
        gammascal = -3.46147417;
        break;
      case 25:
        gammascal = -3.32816812;
        break;
      case 26:
        gammascal = 1.08658130;
        break;
      case 27:
        gammascal = 10.00917344;
        break;
      case 28:
        gammascal = -1.63541603;
        break;
      case 29:
        gammascal = -1.74976633;
        break;
      case 30:
        gammascal = -3.33574171;
        break;
      case 31:
        gammascal = -2.76257240;
        break;
      case 32:
        gammascal = -1.59630173;
        break;
      case 33:
        gammascal = -1.89641826;
        break;
      case 34:
        gammascal = -1.19638371;
        break;
      case 35:
        gammascal = -1.55115682;
        break;
      case 36:
        gammascal = -1.07629721;
        break;
      case 37:
        gammascal = -3.36260929;
        break;
      case 38:
        gammascal = 10.88581368;
        break;
      case 39:
        gammascal = -8.37014419;
        break;
      case 40:
        gammascal = -0.08461763;
        break;
      case 41:
        gammascal = -5.86923780;
        break;
      case 42:
        gammascal = -3.61454676;
        break;
      case 43:
        gammascal = -3.53195625;
        break;
      case 44:
        gammascal = -3.66518218;
        break;
      case 45:
        gammascal = -1.71753606;
        break;
      case 46:
        gammascal = -1.44785600;
        break;
      case 47:
        gammascal = -1.98342667;
        break;
      case 48:
        gammascal = -1.95568677;
        break;
      case 49:
        gammascal = -3.29192109;
        break;
      case 50:
        gammascal = -1.57241888;
        break;
      case 51:
        gammascal = -1.77179655;
        break;
      case 52:
        gammascal = -1.86523585;
        break;
      case 53:
        gammascal = -2.15610833;
        break;
      case 54:
        gammascal = -0.98122873;
        break;
      case 55:
        gammascal = -3.57188311;
        break;
      case 56:
        gammascal = -6.16452009;
        break;
      case 72:
        gammascal = -3.04478959;
        break;
      case 73:
        gammascal = -4.54489752;
        break;
      case 74:
        gammascal = 4.29067433;
        break;
      case 75:
        gammascal = -4.17192593;
        break;
      case 76:
        gammascal = -1.64076232;
        break;
      case 77:
        gammascal = -2.32526625;
        break;
      case 78:
        gammascal = -0.75662991;
        break;
      case 79:
        gammascal = -0.90385530;
        break;
      case 80:
        gammascal = -2.51763493;
        break;
      case 81:
        gammascal = -1.85485514;
        break;
      case 82:
        gammascal = -0.96014493;
        break;
      case 83:
        gammascal = -0.59400861;
        break;
    }
  }
  else if (solvent == "nitromethane") {
    switch (atmnr) {
      case 1:
        gammascal = -1.47954453;
        break;
      case 2:
        gammascal = 0.38952459;
        break;
      case 3:
        gammascal = -7.99375919;
        break;
      case 4:
        gammascal = -22.35078225;
        break;
      case 5:
        gammascal = -2.91089286;
        break;
      case 6:
        gammascal = -2.31214284;
        break;
      case 7:
        gammascal = -5.72766908;
        break;
      case 8:
        gammascal = -2.56534156;
        break;
      case 9:
        gammascal = 0.12654011;
        break;
      case 10:
        gammascal = 0.29685423;
        break;
      case 11:
        gammascal = -5.15086929;
        break;
      case 12:
        gammascal = -3.73148185;
        break;
      case 13:
        gammascal = -1.47444504;
        break;
      case 14:
        gammascal = 0.88472679;
        break;
      case 15:
        gammascal = -1.70086816;
        break;
      case 16:
        gammascal = -2.23408026;
        break;
      case 17:
        gammascal = -1.39779591;
        break;
      case 18:
        gammascal = -0.41947036;
        break;
      case 19:
        gammascal = -3.67595497;
        break;
      case 20:
        gammascal = -8.84845671;
        break;
      case 21:
        gammascal = -3.21242227;
        break;
      case 22:
        gammascal = -3.08630839;
        break;
      case 23:
        gammascal = -4.49485623;
        break;
      case 24:
        gammascal = 2.50659683;
        break;
      case 25:
        gammascal = -4.35823038;
        break;
      case 26:
        gammascal = 3.72961758;
        break;
      case 27:
        gammascal = 16.30898893;
        break;
      case 28:
        gammascal = -0.61183584;
        break;
      case 29:
        gammascal = -1.48392732;
        break;
      case 30:
        gammascal = -4.82145932;
        break;
      case 31:
        gammascal = -3.59054893;
        break;
      case 32:
        gammascal = -0.64794863;
        break;
      case 33:
        gammascal = -2.48330309;
        break;
      case 34:
        gammascal = -1.48852864;
        break;
      case 35:
        gammascal = -1.54446697;
        break;
      case 36:
        gammascal = -0.49591744;
        break;
      case 37:
        gammascal = -3.48949272;
        break;
      case 38:
        gammascal = 11.16347212;
        break;
      case 39:
        gammascal = -10.98193865;
        break;
      case 40:
        gammascal = -2.30355673;
        break;
      case 41:
        gammascal = -5.52666306;
        break;
      case 42:
        gammascal = -2.94829269;
        break;
      case 43:
        gammascal = 2.64731036;
        break;
      case 44:
        gammascal = -1.63047574;
        break;
      case 45:
        gammascal = 0.40465711;
        break;
      case 46:
        gammascal = -2.80106522;
        break;
      case 47:
        gammascal = -2.23096480;
        break;
      case 48:
        gammascal = -3.04992350;
        break;
      case 49:
        gammascal = -3.84516620;
        break;
      case 50:
        gammascal = -0.78652326;
        break;
      case 51:
        gammascal = -2.15576776;
        break;
      case 52:
        gammascal = -2.16532595;
        break;
      case 53:
        gammascal = -2.12183124;
        break;
      case 54:
        gammascal = -0.41186652;
        break;
      case 55:
        gammascal = -5.46492125;
        break;
      case 56:
        gammascal = -8.41525424;
        break;
      case 72:
        gammascal = -3.89300514;
        break;
      case 73:
        gammascal = -4.77046824;
        break;
      case 74:
        gammascal = 5.60494472;
        break;
      case 75:
        gammascal = 5.59934195;
        break;
      case 76:
        gammascal = -0.06192965;
        break;
      case 77:
        gammascal = -0.19372051;
        break;
      case 78:
        gammascal = -0.76424414;
        break;
      case 79:
        gammascal = -0.62242422;
        break;
      case 80:
        gammascal = -2.71154153;
        break;
      case 81:
        gammascal = -2.37541089;
        break;
      case 82:
        gammascal = -0.40184269;
        break;
      case 83:
        gammascal = -1.44436799;
        break;
    }
  }
  else if (solvent == "octanol") {
    switch (atmnr) {
      case 1:
        gammascal = -0.99751210;
        break;
      case 2:
        gammascal = 0.14126340;
        break;
      case 3:
        gammascal = -5.94420799;
        break;
      case 4:
        gammascal = -6.21933351;
        break;
      case 5:
        gammascal = -1.91336597;
        break;
      case 6:
        gammascal = -2.58084256;
        break;
      case 7:
        gammascal = -1.39167721;
        break;
      case 8:
        gammascal = 1.18610532;
        break;
      case 9:
        gammascal = 0.66287534;
        break;
      case 10:
        gammascal = -0.02352835;
        break;
      case 11:
        gammascal = -11.28613227;
        break;
      case 12:
        gammascal = -3.04087824;
        break;
      case 13:
        gammascal = -2.34203056;
        break;
      case 14:
        gammascal = -0.13651019;
        break;
      case 15:
        gammascal = 3.00190223;
        break;
      case 16:
        gammascal = -1.14716215;
        break;
      case 17:
        gammascal = -1.38358737;
        break;
      case 18:
        gammascal = -0.75196120;
        break;
      case 19:
        gammascal = -0.82964249;
        break;
      case 20:
        gammascal = -2.61116785;
        break;
      case 21:
        gammascal = -4.35500206;
        break;
      case 22:
        gammascal = -1.54577684;
        break;
      case 23:
        gammascal = -4.58452949;
        break;
      case 24:
        gammascal = -8.81484099;
        break;
      case 25:
        gammascal = -4.73458013;
        break;
      case 26:
        gammascal = -6.00637700;
        break;
      case 27:
        gammascal = -4.65653689;
        break;
      case 28:
        gammascal = -1.74048544;
        break;
      case 29:
        gammascal = -1.42017360;
        break;
      case 30:
        gammascal = -3.92383157;
        break;
      case 31:
        gammascal = -3.04313605;
        break;
      case 32:
        gammascal = -1.72539307;
        break;
      case 33:
        gammascal = -2.93434489;
        break;
      case 34:
        gammascal = -1.39984928;
        break;
      case 35:
        gammascal = -1.71850451;
        break;
      case 36:
        gammascal = -0.87425568;
        break;
      case 37:
        gammascal = -1.93902191;
        break;
      case 38:
        gammascal = -3.86160814;
        break;
      case 39:
        gammascal = -4.52377738;
        break;
      case 40:
        gammascal = 3.65228470;
        break;
      case 41:
        gammascal = -3.56700546;
        break;
      case 42:
        gammascal = -3.33888912;
        break;
      case 43:
        gammascal = -7.82820091;
        break;
      case 44:
        gammascal = -6.18641633;
        break;
      case 45:
        gammascal = -1.51355814;
        break;
      case 46:
        gammascal = -1.89776883;
        break;
      case 47:
        gammascal = -2.53294592;
        break;
      case 48:
        gammascal = -2.52832288;
        break;
      case 49:
        gammascal = -3.22643150;
        break;
      case 50:
        gammascal = -1.41870327;
        break;
      case 51:
        gammascal = -2.75152789;
        break;
      case 52:
        gammascal = -2.57548249;
        break;
      case 53:
        gammascal = -2.25886984;
        break;
      case 54:
        gammascal = -0.78149273;
        break;
      case 55:
        gammascal = -1.84324509;
        break;
      case 56:
        gammascal = -4.72184513;
        break;
      case 72:
        gammascal = 1.24832154;
        break;
      case 73:
        gammascal = -2.29708010;
        break;
      case 74:
        gammascal = -0.59154657;
        break;
      case 75:
        gammascal = -7.19400147;
        break;
      case 76:
        gammascal = -5.92283414;
        break;
      case 77:
        gammascal = -2.62344386;
        break;
      case 78:
        gammascal = -2.44415941;
        break;
      case 79:
        gammascal = -2.04678953;
        break;
      case 80:
        gammascal = -2.25621066;
        break;
      case 81:
        gammascal = -2.59752403;
        break;
      case 82:
        gammascal = -2.26874279;
        break;
      case 83:
        gammascal = -1.73826293;
        break;
    }
  }
  else if (solvent == "phenol") {
    switch (atmnr) {
      case 1:
        gammascal = -0.96839923;
        break;
      case 2:
        gammascal = -0.00900322;
        break;
      case 3:
        gammascal = -5.58352520;
        break;
      case 4:
        gammascal = -17.67562525;
        break;
      case 5:
        gammascal = -1.91227300;
        break;
      case 6:
        gammascal = -2.05126808;
        break;
      case 7:
        gammascal = -2.03602629;
        break;
      case 8:
        gammascal = 0.14465721;
        break;
      case 9:
        gammascal = 0.82724006;
        break;
      case 10:
        gammascal = -0.10703607;
        break;
      case 11:
        gammascal = -3.55307037;
        break;
      case 12:
        gammascal = -2.06307394;
        break;
      case 13:
        gammascal = -1.09876276;
        break;
      case 14:
        gammascal = 1.17369466;
        break;
      case 15:
        gammascal = -1.70191289;
        break;
      case 16:
        gammascal = -1.33447061;
        break;
      case 17:
        gammascal = -0.57641270;
        break;
      case 18:
        gammascal = -1.01075068;
        break;
      case 19:
        gammascal = -2.88120550;
        break;
      case 20:
        gammascal = -10.49749479;
        break;
      case 21:
        gammascal = -7.16133016;
        break;
      case 22:
        gammascal = -1.81547354;
        break;
      case 23:
        gammascal = -4.97224150;
        break;
      case 24:
        gammascal = -2.63002803;
        break;
      case 25:
        gammascal = -1.29258922;
        break;
      case 26:
        gammascal = 0.57687093;
        break;
      case 27:
        gammascal = 0.00000000;
        break;
      case 28:
        gammascal = -0.39543776;
        break;
      case 29:
        gammascal = -0.95163672;
        break;
      case 30:
        gammascal = -2.74863366;
        break;
      case 31:
        gammascal = -3.09862387;
        break;
      case 32:
        gammascal = -1.38820677;
        break;
      case 33:
        gammascal = -1.83557031;
        break;
      case 34:
        gammascal = -0.85345969;
        break;
      case 35:
        gammascal = -1.43511796;
        break;
      case 36:
        gammascal = -1.04909364;
        break;
      case 37:
        gammascal = -3.62341517;
        break;
      case 38:
        gammascal = -11.01027657;
        break;
      case 39:
        gammascal = -7.53075769;
        break;
      case 40:
        gammascal = 1.83364257;
        break;
      case 41:
        gammascal = -6.66701430;
        break;
      case 42:
        gammascal = -1.80598377;
        break;
      case 43:
        gammascal = -0.02824780;
        break;
      case 44:
        gammascal = -4.04091704;
        break;
      case 45:
        gammascal = -0.84747254;
        break;
      case 46:
        gammascal = -1.16957244;
        break;
      case 47:
        gammascal = -1.78343067;
        break;
      case 48:
        gammascal = -2.71656682;
        break;
      case 49:
        gammascal = -3.82893776;
        break;
      case 50:
        gammascal = -1.47390709;
        break;
      case 51:
        gammascal = -1.65873077;
        break;
      case 52:
        gammascal = -1.72434474;
        break;
      case 53:
        gammascal = -2.02651720;
        break;
      case 54:
        gammascal = -0.95822200;
        break;
      case 55:
        gammascal = -3.97555323;
        break;
      case 56:
        gammascal = -9.27732995;
        break;
      case 72:
        gammascal = -0.13556388;
        break;
      case 73:
        gammascal = -5.02103977;
        break;
      case 74:
        gammascal = 10.45055312;
        break;
      case 75:
        gammascal = 0.90327146;
        break;
      case 76:
        gammascal = -1.96126586;
        break;
      case 77:
        gammascal = -1.97086444;
        break;
      case 78:
        gammascal = -0.42030415;
        break;
      case 79:
        gammascal = -0.66030833;
        break;
      case 80:
        gammascal = -2.41319691;
        break;
      case 81:
        gammascal = -1.53856791;
        break;
      case 82:
        gammascal = -0.65689063;
        break;
      case 83:
        gammascal = -0.67200346;
        break;
    }
  }
  else if (solvent == "thf") {
    switch (atmnr) {
      case 1:
        gammascal = -1.39382608;
        break;
      case 2:
        gammascal = -0.41673512;
        break;
      case 3:
        gammascal = -5.24166548;
        break;
      case 4:
        gammascal = -12.28002293;
        break;
      case 5:
        gammascal = -2.39879944;
        break;
      case 6:
        gammascal = -2.37361122;
        break;
      case 7:
        gammascal = -2.28412966;
        break;
      case 8:
        gammascal = -1.29340598;
        break;
      case 9:
        gammascal = -0.59946682;
        break;
      case 10:
        gammascal = -0.45617923;
        break;
      case 11:
        gammascal = -3.28531560;
        break;
      case 12:
        gammascal = -2.60390263;
        break;
      case 13:
        gammascal = -2.01369319;
        break;
      case 14:
        gammascal = 0.16334916;
        break;
      case 15:
        gammascal = -2.19322298;
        break;
      case 16:
        gammascal = -2.07853613;
        break;
      case 17:
        gammascal = -1.65297751;
        break;
      case 18:
        gammascal = -1.34212056;
        break;
      case 19:
        gammascal = -2.83157992;
        break;
      case 20:
        gammascal = -5.21371290;
        break;
      case 21:
        gammascal = -5.30650410;
        break;
      case 22:
        gammascal = -4.42395931;
        break;
      case 23:
        gammascal = -4.20357161;
        break;
      case 24:
        gammascal = -4.46225699;
        break;
      case 25:
        gammascal = -2.76811704;
        break;
      case 26:
        gammascal = 3.22683787;
        break;
      case 27:
        gammascal = 4.92263961;
        break;
      case 28:
        gammascal = -2.07981512;
        break;
      case 29:
        gammascal = -2.08023751;
        break;
      case 30:
        gammascal = -3.35541496;
        break;
      case 31:
        gammascal = -2.66621144;
        break;
      case 32:
        gammascal = -2.00734620;
        break;
      case 33:
        gammascal = -2.22988638;
        break;
      case 34:
        gammascal = -1.67860150;
        break;
      case 35:
        gammascal = -1.66179813;
        break;
      case 36:
        gammascal = -1.36804717;
        break;
      case 37:
        gammascal = -2.90001002;
        break;
      case 38:
        gammascal = -5.88833988;
        break;
      case 39:
        gammascal = -7.16970264;
        break;
      case 40:
        gammascal = -0.21731622;
        break;
      case 41:
        gammascal = -3.60458626;
        break;
      case 42:
        gammascal = -3.77610601;
        break;
      case 43:
        gammascal = -1.40220558;
        break;
      case 44:
        gammascal = -3.05574179;
        break;
      case 45:
        gammascal = -1.34006342;
        break;
      case 46:
        gammascal = -2.23355539;
        break;
      case 47:
        gammascal = -2.28942451;
        break;
      case 48:
        gammascal = -2.90334663;
        break;
      case 49:
        gammascal = -2.87078852;
        break;
      case 50:
        gammascal = -1.79385735;
        break;
      case 51:
        gammascal = -1.86935973;
        break;
      case 52:
        gammascal = -2.13111698;
        break;
      case 53:
        gammascal = -2.21044812;
        break;
      case 54:
        gammascal = -1.28526914;
        break;
      case 55:
        gammascal = -3.06455933;
        break;
      case 56:
        gammascal = -5.03972111;
        break;
      case 72:
        gammascal = -1.58111347;
        break;
      case 73:
        gammascal = -2.63438997;
        break;
      case 74:
        gammascal = 1.92177355;
        break;
      case 75:
        gammascal = -1.48283989;
        break;
      case 76:
        gammascal = -1.33964787;
        break;
      case 77:
        gammascal = -2.00840956;
        break;
      case 78:
        gammascal = -1.11872636;
        break;
      case 79:
        gammascal = -1.29675924;
        break;
      case 80:
        gammascal = -2.74774148;
        break;
      case 81:
        gammascal = -2.00042498;
        break;
      case 82:
        gammascal = -1.04613567;
        break;
      case 83:
        gammascal = -1.29528026;
        break;
    }
  }
  else if (solvent == "toluene") {
    switch (atmnr) {
      case 1:
        gammascal = -1.57713161;
        break;
      case 2:
        gammascal = -0.34551336;
        break;
      case 3:
        gammascal = -6.30332548;
        break;
      case 4:
        gammascal = -18.98649997;
        break;
      case 5:
        gammascal = -3.66060739;
        break;
      case 6:
        gammascal = -3.00030249;
        break;
      case 7:
        gammascal = -4.17462021;
        break;
      case 8:
        gammascal = -3.12456571;
        break;
      case 9:
        gammascal = -0.34762131;
        break;
      case 10:
        gammascal = -0.34387770;
        break;
      case 11:
        gammascal = -3.79758210;
        break;
      case 12:
        gammascal = -3.48313415;
        break;
      case 13:
        gammascal = -2.64930731;
        break;
      case 14:
        gammascal = 0.06877995;
        break;
      case 15:
        gammascal = -3.21112656;
        break;
      case 16:
        gammascal = -3.06348464;
        break;
      case 17:
        gammascal = -2.09906469;
        break;
      case 18:
        gammascal = -1.54101094;
        break;
      case 19:
        gammascal = -2.83099083;
        break;
      case 20:
        gammascal = -4.98641655;
        break;
      case 21:
        gammascal = -4.80452663;
        break;
      case 22:
        gammascal = -2.01073769;
        break;
      case 23:
        gammascal = -3.03998118;
        break;
      case 24:
        gammascal = 2.18902787;
        break;
      case 25:
        gammascal = 8.35600917;
        break;
      case 26:
        gammascal = 2.90721631;
        break;
      case 27:
        gammascal = 8.20948262;
        break;
      case 28:
        gammascal = -1.66062338;
        break;
      case 29:
        gammascal = -1.98564808;
        break;
      case 30:
        gammascal = -3.87337313;
        break;
      case 31:
        gammascal = -3.62089746;
        break;
      case 32:
        gammascal = -1.89074832;
        break;
      case 33:
        gammascal = -3.02387557;
        break;
      case 34:
        gammascal = -2.31730325;
        break;
      case 35:
        gammascal = -2.83207161;
        break;
      case 36:
        gammascal = -1.66381358;
        break;
      case 37:
        gammascal = -2.90724772;
        break;
      case 38:
        gammascal = -6.48619927;
        break;
      case 39:
        gammascal = -9.04057593;
        break;
      case 40:
        gammascal = -2.91113008;
        break;
      case 41:
        gammascal = -3.29996433;
        break;
      case 42:
        gammascal = -1.51966800;
        break;
      case 43:
        gammascal = 4.05033533;
        break;
      case 44:
        gammascal = -0.57755575;
        break;
      case 45:
        gammascal = 2.41836271;
        break;
      case 46:
        gammascal = -2.89821150;
        break;
      case 47:
        gammascal = -2.83024392;
        break;
      case 48:
        gammascal = -1.82881890;
        break;
      case 49:
        gammascal = -3.59606977;
        break;
      case 50:
        gammascal = -1.61888881;
        break;
      case 51:
        gammascal = -2.67083583;
        break;
      case 52:
        gammascal = -2.76398086;
        break;
      case 53:
        gammascal = -3.30677278;
        break;
      case 54:
        gammascal = -1.02925461;
        break;
      case 55:
        gammascal = -2.53155988;
        break;
      case 56:
        gammascal = -5.26980224;
        break;
      case 72:
        gammascal = -4.47869281;
        break;
      case 73:
        gammascal = -3.86512339;
        break;
      case 74:
        gammascal = 1.18683963;
        break;
      case 75:
        gammascal = 3.22740259;
        break;
      case 76:
        gammascal = 0.18806102;
        break;
      case 77:
        gammascal = 1.55354823;
        break;
      case 78:
        gammascal = -1.99571944;
        break;
      case 79:
        gammascal = -1.29611030;
        break;
      case 80:
        gammascal = -2.16555007;
        break;
      case 81:
        gammascal = -2.53244490;
        break;
      case 82:
        gammascal = -1.13144107;
        break;
      case 83:
        gammascal = -1.80817292;
        break;
    }
  }
  else if (solvent == "water") {
    switch (atmnr) {
      case 1:
        gammascal = -0.08533368;
        break;
      case 2:
        gammascal = 1.13711777;
        break;
      case 3:
        gammascal = -5.83236188;
        break;
      case 4:
        gammascal = -17.01350889;
        break;
      case 5:
        gammascal = -0.82891231;
        break;
      case 6:
        gammascal = -0.51385188;
        break;
      case 7:
        gammascal = -3.24932941;
        break;
      case 8:
        gammascal = 2.03601297;
        break;
      case 9:
        gammascal = 1.26697909;
        break;
      case 10:
        gammascal = 1.00726547;
        break;
      case 11:
        gammascal = -3.05732760;
        break;
      case 12:
        gammascal = -1.61595951;
        break;
      case 13:
        gammascal = -0.10749537;
        break;
      case 14:
        gammascal = 1.73814581;
        break;
      case 15:
        gammascal = 1.01416364;
        break;
      case 16:
        gammascal = 1.20810520;
        break;
      case 17:
        gammascal = -0.06835879;
        break;
      case 18:
        gammascal = 0.43292409;
        break;
      case 19:
        gammascal = -2.69545448;
        break;
      case 20:
        gammascal = -9.68055177;
        break;
      case 21:
        gammascal = -6.86171080;
        break;
      case 22:
        gammascal = -2.32183024;
        break;
      case 23:
        gammascal = -2.10041202;
        break;
      case 24:
        gammascal = -4.15076452;
        break;
      case 25:
        gammascal = -2.32406972;
        break;
      case 26:
        gammascal = 9.24779587;
        break;
      case 27:
        gammascal = 4.48206277;
        break;
      case 28:
        gammascal = -1.15972411;
        break;
      case 29:
        gammascal = -0.21281688;
        break;
      case 30:
        gammascal = -2.64500225;
        break;
      case 31:
        gammascal = -2.36086956;
        break;
      case 32:
        gammascal = -0.46303904;
        break;
      case 33:
        gammascal = -0.88149455;
        break;
      case 34:
        gammascal = 0.23523157;
        break;
      case 35:
        gammascal = -0.18620262;
        break;
      case 36:
        gammascal = 0.35105321;
        break;
      case 37:
        gammascal = -2.88851792;
        break;
      case 38:
        gammascal = -11.04377179;
        break;
      case 39:
        gammascal = -7.95128133;
        break;
      case 40:
        gammascal = 3.90876499;
        break;
      case 41:
        gammascal = -2.48432528;
        break;
      case 42:
        gammascal = -3.66936332;
        break;
      case 43:
        gammascal = -4.43332314;
        break;
      case 44:
        gammascal = -2.94937418;
        break;
      case 45:
        gammascal = 0.00028790;
        break;
      case 46:
        gammascal = -0.93137790;
        break;
      case 47:
        gammascal = -0.79778296;
        break;
      case 48:
        gammascal = -0.92747581;
        break;
      case 49:
        gammascal = -2.70394304;
        break;
      case 50:
        gammascal = -0.43878679;
        break;
      case 51:
        gammascal = -0.70393148;
        break;
      case 52:
        gammascal = -0.77480977;
        break;
      case 53:
        gammascal = -0.76873446;
        break;
      case 54:
        gammascal = -0.06431749;
        break;
      case 55:
        gammascal = -3.15995511;
        break;
      case 56:
        gammascal = -5.92642054;
        break;
      case 72:
        gammascal = -0.91732017;
        break;
      case 73:
        gammascal = -0.76038638;
        break;
      case 74:
        gammascal = 4.55364802;
        break;
      case 75:
        gammascal = -5.19397805;
        break;
      case 76:
        gammascal = -0.97455175;
        break;
      case 77:
        gammascal = -0.19095469;
        break;
      case 78:
        gammascal = 0.37263783;
        break;
      case 79:
        gammascal = 0.41221465;
        break;
      case 80:
        gammascal = -1.77134941;
        break;
      case 81:
        gammascal = -0.89456867;
        break;
      case 82:
        gammascal = 0.24690462;
        break;
      case 83:
        gammascal = 0.62621722;
        break;
    }
  }
  else if (solvent == "octanol wet") {
    switch (atmnr) {
      case 1:
        gammascal = -1.16531035;
        break;
      case 2:
        gammascal = 0.15400251;
        break;
      case 3:
        gammascal = -7.29505708;
        break;
      case 4:
        gammascal = -6.92108479;
        break;
      case 5:
        gammascal = -1.90066124;
        break;
      case 6:
        gammascal = -2.39668191;
        break;
      case 7:
        gammascal = -1.93160159;
        break;
      case 8:
        gammascal = 1.09209376;
        break;
      case 9:
        gammascal = 0.32771424;
        break;
      case 10:
        gammascal = -0.00692900;
        break;
      case 11:
        gammascal = -10.04402154;
        break;
      case 12:
        gammascal = -3.29244006;
        break;
      case 13:
        gammascal = -2.11370570;
        break;
      case 14:
        gammascal = -0.17668798;
        break;
      case 15:
        gammascal = -2.22054477;
        break;
      case 16:
        gammascal = -2.08917815;
        break;
      case 17:
        gammascal = -1.67594007;
        break;
      case 18:
        gammascal = -0.72693210;
        break;
      case 19:
        gammascal = -1.41680106;
        break;
      case 20:
        gammascal = -3.25078235;
        break;
      case 21:
        gammascal = -4.31771202;
        break;
      case 22:
        gammascal = -3.56728566;
        break;
      case 23:
        gammascal = -2.79213857;
        break;
      case 24:
        gammascal = -7.17461986;
        break;
      case 25:
        gammascal = -5.03691775;
        break;
      case 26:
        gammascal = -4.12095913;
        break;
      case 27:
        gammascal = -3.11311579;
        break;
      case 28:
        gammascal = -1.67264344;
        break;
      case 29:
        gammascal = -1.53159485;
        break;
      case 30:
        gammascal = -4.10993799;
        break;
      case 31:
        gammascal = -2.92842594;
        break;
      case 32:
        gammascal = -2.13168832;
        break;
      case 33:
        gammascal = -2.85015528;
        break;
      case 34:
        gammascal = -1.49735630;
        break;
      case 35:
        gammascal = -2.02695519;
        break;
      case 36:
        gammascal = -0.82672883;
        break;
      case 37:
        gammascal = -2.16257109;
        break;
      case 38:
        gammascal = -4.30867113;
        break;
      case 39:
        gammascal = -5.48068119;
        break;
      case 40:
        gammascal = 2.96187570;
        break;
      case 41:
        gammascal = -3.38544740;
        break;
      case 42:
        gammascal = -1.72989709;
        break;
      case 43:
        gammascal = -7.03494832;
        break;
      case 44:
        gammascal = -4.81575078;
        break;
      case 45:
        gammascal = -1.12140213;
        break;
      case 46:
        gammascal = -2.37375979;
        break;
      case 47:
        gammascal = -2.55847796;
        break;
      case 48:
        gammascal = -1.88753081;
        break;
      case 49:
        gammascal = -3.05082878;
        break;
      case 50:
        gammascal = -1.29182261;
        break;
      case 51:
        gammascal = -2.66606103;
        break;
      case 52:
        gammascal = -2.49069454;
        break;
      case 53:
        gammascal = -2.19407569;
        break;
      case 54:
        gammascal = -0.76000453;
        break;
      case 55:
        gammascal = -2.14509302;
        break;
      case 56:
        gammascal = -4.94583201;
        break;
      case 72:
        gammascal = -0.57355318;
        break;
      case 73:
        gammascal = -2.14357748;
        break;
      case 74:
        gammascal = -3.46058590;
        break;
      case 75:
        gammascal = -7.44068284;
        break;
      case 76:
        gammascal = -5.63337544;
        break;
      case 77:
        gammascal = -2.44278043;
        break;
      case 78:
        gammascal = -2.20054033;
        break;
      case 79:
        gammascal = -2.03319381;
        break;
      case 80:
        gammascal = -2.11730852;
        break;
      case 81:
        gammascal = -2.40760023;
        break;
      case 82:
        gammascal = -2.12279947;
        break;
      case 83:
        gammascal = -1.71809404;
        break;
    }
  }
  return gammascal;
}
double HbondStrengthGFN2(std::string solvent, size_t atmnr) {
  //function returning the hydrogen-bonding strength of each atom for generalized Born model of solvation for GFN2 method
  //data taken from xTB
  double hbondstrength = 0.0;
  if (solvent == "acetone") {
    if (atmnr == 1) {hbondstrength = 6.76435360;}
  }
  else if (solvent == "acetonitrile") {
    if (atmnr == 1) {hbondstrength = 6.01652192;}
  }
  else if (solvent == "aniline") {
    if (atmnr == 1) {hbondstrength = 5.07468979;}
    else if (atmnr == 6) {hbondstrength = 0.16700779;}
    else if (atmnr == 7) {hbondstrength = 0.19747334;}
    else if (atmnr == 8) {hbondstrength = 0.73961102;}
    else if (atmnr == 9) {hbondstrength = 0.92008779;}
    else if (atmnr == 15) {hbondstrength = 1.16300541;}
    else if (atmnr == 16) {hbondstrength = 2.66688398;}
    else if (atmnr == 35) {hbondstrength = 2.45236657;}
    else if (atmnr == 53) {hbondstrength = 2.01252622;}
  }
  else if (solvent == "benzaldehyde") {
    if (atmnr == 1) {hbondstrength = 3.07480768;}
  }
  else if (solvent == "benzene") {hbondstrength = 0.0;}
  else if (solvent == "dichloromethane") {hbondstrength = 0.0;}
  else if (solvent == "chloroform") {hbondstrength = 0.0;}
  else if (solvent == "carbon disulfide") {hbondstrength = 0.0;}
  else if (solvent == "dioxane") {
    if (atmnr == 1) {hbondstrength = 5.79628587;}
  }
  else if (solvent == "dmf") {
    if (atmnr == 1) {hbondstrength = 7.61302606;}
  }
  else if (solvent == "dmso") {
    if (atmnr == 1) {hbondstrength = 8.48193404;}
  }
  else if (solvent == "ethanol") {
    if (atmnr == 1) {hbondstrength = 6.42795353;}
    else if (atmnr == 6) {hbondstrength = 2.17349835;}
    else if (atmnr == 7) {hbondstrength = 2.61707358;}
    else if (atmnr == 8) {hbondstrength = 0.46312477;}
    else if (atmnr == 9) {hbondstrength = 0.01355631;}
    else if (atmnr == 15) {hbondstrength = 2.78593962;}
    else if (atmnr == 16) {hbondstrength = 0.10540762;}
    else if (atmnr == 17) {hbondstrength = 0.06513631;}
  }
  else if (solvent == "diethyl ether") {
    if (atmnr == 1) {hbondstrength = 7.09710513;}
  }
  else if (solvent == "ethyl acetate") {
    if (atmnr == 1) {hbondstrength = 4.31040661;}
  }
  else if (solvent == "furane") {
    if (atmnr == 1) {hbondstrength = 0.33594644;}
  }
  else if (solvent == "hexadecane") {hbondstrength = 0.0;}
  else if (solvent == "hexane") {hbondstrength = 0.0;}
  else if (solvent == "methanol") {
    if (atmnr == 1) {hbondstrength = 6.45931259;}
    else if (atmnr == 6) {hbondstrength = 0.00797465;}
    else if (atmnr == 7) {hbondstrength = 2.98775273;}
    else if (atmnr == 8) {hbondstrength = 1.51391089;}
    else if (atmnr == 9) {hbondstrength = 0.51134039;}
    else if (atmnr == 15) {hbondstrength = 2.79185940;}
    else if (atmnr == 16) {hbondstrength = 0.09199181;}
    else if (atmnr == 17) {hbondstrength = 1.93913555;}
    else if (atmnr == 34) {hbondstrength = 2.10235367;}
    else if (atmnr == 35) {hbondstrength = 2.03381074;}
    else if (atmnr == 53) {hbondstrength = 2.40567409;}
  }
  else if (solvent == "nitromethane") {
    if (atmnr == 1) {hbondstrength = 2.65877401;}
  }
  else if (solvent == "octanol") {
    if (atmnr == 1) {hbondstrength = 6.95344958;}
    else if (atmnr == 6) {hbondstrength = 4.25383208;}
    else if (atmnr == 7) {hbondstrength = 3.30143242;}
    else if (atmnr == 8) {hbondstrength = 1.63706074;}
    else if (atmnr == 9) {hbondstrength = 1.30661877;}
    else if (atmnr == 15) {hbondstrength = 3.99388076;}
    else if (atmnr == 16) {hbondstrength = 1.99734857;}
    else if (atmnr == 17) {hbondstrength = 0.39402401;}
    else if (atmnr == 35) {hbondstrength = 1.97758816;}
    else if (atmnr == 53) {hbondstrength = 1.55106831;}
  }
  else if (solvent == "phenol") {
    if (atmnr == 1) {hbondstrength = 3.25508042;}
    else if (atmnr == 6) {hbondstrength = -0.00337953;}
    else if (atmnr == 7) {hbondstrength = 2.93395996;}
    else if (atmnr == 8) {hbondstrength = 1.50994779;}
    else if (atmnr == 9) {hbondstrength = 0.51818401;}
    else if (atmnr == 15) {hbondstrength = 1.98710351;}
    else if (atmnr == 16) {hbondstrength = 2.04284860;}
    else if (atmnr == 35) {hbondstrength = 2.89274707;}
    else if (atmnr == 53) {hbondstrength = 2.29819357;}
  }
  else if (solvent == "thf") {
    if (atmnr == 1) {hbondstrength = 8.03952870;}
  }
  else if (solvent == "toluene") {hbondstrength = 0.0;}
  else if (solvent == "water") {
    if (atmnr == 1) {hbondstrength = 8.09964704;}
    else if (atmnr == 6) {hbondstrength = 2.23420044;}
    else if (atmnr == 7) {hbondstrength = 1.80633279;}
    else if (atmnr == 8) {hbondstrength = 2.22319193;}
    else if (atmnr == 9) {hbondstrength = 3.56359195;}
    else if (atmnr == 15) {hbondstrength = 3.62926820;}
    else if (atmnr == 16) {hbondstrength = 0.26212102;}
    else if (atmnr == 17) {hbondstrength = 0.15334756;}
    else if (atmnr == 34) {hbondstrength = 2.37162766;}
    else if (atmnr == 35) {hbondstrength = 2.02275702;}
    else if (atmnr == 53) {hbondstrength = 2.61716906;}
  }
  else if (solvent == "octanol wet") {
    if (atmnr == 1) {hbondstrength = 6.78861979;}
    else if (atmnr == 6) {hbondstrength = 2.92814324;}
    else if (atmnr == 7) {hbondstrength = 3.68750347;}
    else if (atmnr == 8) {hbondstrength = 1.48885028;}
    else if (atmnr == 9) {hbondstrength = 0.01079927;}
    else if (atmnr == 15) {hbondstrength = 2.02936796;}
    else if (atmnr == 16) {hbondstrength = 0.08907567;}
    else if (atmnr == 17) {hbondstrength = 0.08617914;}
    else if (atmnr == 35) {hbondstrength = 0.01057765;}
    else if (atmnr == 53) {hbondstrength = 3.18651808;}
  }
  return hbondstrength;
}
double DeScreeningGFN2(std::string solvent, size_t atmnr) {
  //descreening constants for generalized Born models with GFN2
  double descreen = 0.0;
  if (solvent == "acetone") {
    switch (atmnr) {
      case 1:
        descreen = 0.48130656;
        break;
      case 2:
        descreen = 1.32818914;
        break;
      case 3:
        descreen = 0.30000000;
        break;
      case 4:
        descreen = 0.82280446;
        break;
      case 5:
        descreen = 0.79342871;
        break;
      case 6:
        descreen = 0.84069542;
        break;
      case 7:
        descreen = 0.30000000;
        break;
      case 8:
        descreen = 0.30000000;
        break;
      case 9:
        descreen = 0.76391376;
        break;
      case 10:
        descreen = 0.30000000;
        break;
      case 11:
        descreen = 0.15633089;
        break;
      case 12:
        descreen = 1.09358429;
        break;
      case 13:
        descreen = 0.69841172;
        break;
      case 14:
        descreen = 0.87924324;
        break;
      case 15:
        descreen = 0.76121319;
        break;
      case 16:
        descreen = 0.94722356;
        break;
      case 17:
        descreen = 0.74062885;
        break;
      case 18:
        descreen = 0.30000000;
        break;
      case 19:
        descreen = 0.30000000;
        break;
      case 20:
        descreen = 0.30000000;
        break;
      case 21:
        descreen = 1.02523016;
        break;
      case 22:
        descreen = 1.06672104;
        break;
      case 23:
        descreen = 1.19745079;
        break;
      case 24:
        descreen = 1.37842025;
        break;
      case 25:
        descreen = 1.13846943;
        break;
      case 26:
        descreen = 1.42662550;
        break;
      case 27:
        descreen = 0.84079786;
        break;
      case 28:
        descreen = 0.89696274;
        break;
      case 29:
        descreen = 0.49020444;
        break;
      case 30:
        descreen = 0.30000000;
        break;
      case 31:
        descreen = 1.11091916;
        break;
      case 32:
        descreen = 0.85944054;
        break;
      case 33:
        descreen = 1.00706733;
        break;
      case 34:
        descreen = 1.21342516;
        break;
      case 35:
        descreen = 0.59847219;
        break;
      case 36:
        descreen = 1.20576503;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 0.30000000;
        break;
      case 39:
        descreen = 1.20873605;
        break;
      case 40:
        descreen = 0.71782314;
        break;
      case 41:
        descreen = 1.23226470;
        break;
      case 42:
        descreen = 1.32705365;
        break;
      case 43:
        descreen = 1.30515572;
        break;
      case 44:
        descreen = 1.35680638;
        break;
      case 45:
        descreen = 1.26591080;
        break;
      case 46:
        descreen = 0.67728437;
        break;
      case 47:
        descreen = 0.98164468;
        break;
      case 48:
        descreen = 1.14661537;
        break;
      case 49:
        descreen = 0.86519215;
        break;
      case 50:
        descreen = 0.64492595;
        break;
      case 51:
        descreen = 1.08495415;
        break;
      case 52:
        descreen = 0.30000000;
        break;
      case 53:
        descreen = 0.65422829;
        break;
      case 54:
        descreen = 1.14641416;
        break;
      case 55:
        descreen = 0.61548674;
        break;
      case 56:
        descreen = 0.97521406;
        break;
      case 57:
        descreen = 0.80000000;
        break;
      case 58:
        descreen = 0.80000000;
        break;
      case 59:
        descreen = 0.80000000;
        break;
      case 60:
        descreen = 0.80000000;
        break;
      case 61:
        descreen = 0.80000000;
        break;
      case 62:
        descreen = 0.80000000;
        break;
      case 63:
        descreen = 0.80000000;
        break;
      case 64:
        descreen = 0.80000000;
        break;
      case 65:
        descreen = 0.80000000;
        break;
      case 66:
        descreen = 0.80000000;
        break;
      case 67:
        descreen = 0.80000000;
        break;
      case 68:
        descreen = 0.80000000;
        break;
      case 69:
        descreen = 0.80000000;
        break;
      case 70:
        descreen = 0.80000000;
        break;
      case 71:
        descreen = 0.80000000;
        break;
      case 72:
        descreen = 1.27818607;
        break;
      case 73:
        descreen = 1.08627283;
        break;
      case 74:
        descreen = 1.12592398;
        break;
      case 75:
        descreen = 1.04507244;
        break;
      case 76:
        descreen = 1.29039929;
        break;
      case 77:
        descreen = 1.23173428;
        break;
      case 78:
        descreen = 0.73926867;
        break;
      case 79:
        descreen = 1.02083793;
        break;
      case 80:
        descreen = 0.30000000;
        break;
      case 81:
        descreen = 0.84958526;
        break;
      case 82:
        descreen = 0.98861174;
        break;
      case 83:
        descreen = 0.53599756;
        break;
      case 84:
        descreen = 0.80000000;
        break;
      case 85:
        descreen = 0.80000000;
        break;
      case 86:
        descreen = 0.80000000;
        break;
      case 87:
        descreen = 0.80000000;
        break;
      case 88:
        descreen = 0.80000000;
        break;
      case 89:
        descreen = 0.80000000;
        break;
      case 90:
        descreen = 0.80000000;
        break;
      case 91:
        descreen = 0.80000000;
        break;
      case 92:
        descreen = 0.80000000;
        break;
      case 93:
        descreen = 0.80000000;
        break;
      case 94:
        descreen = 0.80000000;
        break;
    }
  }
  else if (solvent == "acetonitrile") {
    switch (atmnr) {
      case 1:
        descreen = 0.57074107;
        break;
      case 2:
        descreen = 1.68244035;
        break;
      case 3:
        descreen = 0.05496586;
        break;
      case 4:
        descreen = 0.80546054;
        break;
      case 5:
        descreen = 0.79513406;
        break;
      case 6:
        descreen = 0.83570107;
        break;
      case 7:
        descreen = 0.70688560;
        break;
      case 8:
        descreen = 0.30000000;
        break;
      case 9:
        descreen = 0.73686633;
        break;
      case 10:
        descreen = 0.30000000;
        break;
      case 11:
        descreen = 0.30000000;
        break;
      case 12:
        descreen = 0.80058936;
        break;
      case 13:
        descreen = 0.71089231;
        break;
      case 14:
        descreen = 0.83507313;
        break;
      case 15:
        descreen = 0.75643432;
        break;
      case 16:
        descreen = 0.93261460;
        break;
      case 17:
        descreen = 0.70992920;
        break;
      case 18:
        descreen = 0.30000000;
        break;
      case 19:
        descreen = 0.30000000;
        break;
      case 20:
        descreen = 0.30000000;
        break;
      case 21:
        descreen = 1.00871604;
        break;
      case 22:
        descreen = 1.05156360;
        break;
      case 23:
        descreen = 1.15737083;
        break;
      case 24:
        descreen = 1.38505998;
        break;
      case 25:
        descreen = 1.11092644;
        break;
      case 26:
        descreen = 1.45293995;
        break;
      case 27:
        descreen = 0.71760627;
        break;
      case 28:
        descreen = 0.85265458;
        break;
      case 29:
        descreen = 0.53117479;
        break;
      case 30:
        descreen = 0.30000000;
        break;
      case 31:
        descreen = 0.96713214;
        break;
      case 32:
        descreen = 0.82433940;
        break;
      case 33:
        descreen = 1.02048126;
        break;
      case 34:
        descreen = 1.05813104;
        break;
      case 35:
        descreen = 0.76000744;
        break;
      case 36:
        descreen = 1.12053043;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 0.30000000;
        break;
      case 39:
        descreen = 1.20076271;
        break;
      case 40:
        descreen = 0.42555962;
        break;
      case 41:
        descreen = 1.18726757;
        break;
      case 42:
        descreen = 1.34130360;
        break;
      case 43:
        descreen = 1.24256188;
        break;
      case 44:
        descreen = 1.32903149;
        break;
      case 45:
        descreen = 1.15531771;
        break;
      case 46:
        descreen = 0.88346638;
        break;
      case 47:
        descreen = 0.86584805;
        break;
      case 48:
        descreen = 0.30000000;
        break;
      case 49:
        descreen = 0.84764033;
        break;
      case 50:
        descreen = 0.51389128;
        break;
      case 51:
        descreen = 0.73022974;
        break;
      case 52:
        descreen = 0.55572302;
        break;
      case 53:
        descreen = 0.65793507;
        break;
      case 54:
        descreen = 1.05764484;
        break;
      case 55:
        descreen = 0.30000000;
        break;
      case 56:
        descreen = 0.95722319;
        break;
      case 57:
        descreen = 0.80000000;
        break;
      case 58:
        descreen = 0.80000000;
        break;
      case 59:
        descreen = 0.80000000;
        break;
      case 60:
        descreen = 0.80000000;
        break;
      case 61:
        descreen = 0.80000000;
        break;
      case 62:
        descreen = 0.80000000;
        break;
      case 63:
        descreen = 0.80000000;
        break;
      case 64:
        descreen = 0.80000000;
        break;
      case 65:
        descreen = 0.80000000;
        break;
      case 66:
        descreen = 0.80000000;
        break;
      case 67:
        descreen = 0.80000000;
        break;
      case 68:
        descreen = 0.80000000;
        break;
      case 69:
        descreen = 0.80000000;
        break;
      case 70:
        descreen = 0.80000000;
        break;
      case 71:
        descreen = 0.80000000;
        break;
      case 72:
        descreen = 1.22867793;
        break;
      case 73:
        descreen = 1.03450363;
        break;
      case 74:
        descreen = 1.14226044;
        break;
      case 75:
        descreen = 0.94780050;
        break;
      case 76:
        descreen = 1.30072500;
        break;
      case 77:
        descreen = 1.22069757;
        break;
      case 78:
        descreen = 1.00206526;
        break;
      case 79:
        descreen = 0.97064876;
        break;
      case 80:
        descreen = 0.80000000;
        break;
      case 81:
        descreen = 0.76551161;
        break;
      case 82:
        descreen = 0.99298353;
        break;
      case 83:
        descreen = 0.41893473;
        break;
      case 84:
        descreen = 0.80000000;
        break;
      case 85:
        descreen = 0.80000000;
        break;
      case 86:
        descreen = 0.80000000;
        break;
      case 87:
        descreen = 0.80000000;
        break;
      case 88:
        descreen = 0.80000000;
        break;
      case 89:
        descreen = 0.80000000;
        break;
      case 90:
        descreen = 0.80000000;
        break;
      case 91:
        descreen = 0.80000000;
        break;
      case 92:
        descreen = 0.80000000;
        break;
      case 93:
        descreen = 0.80000000;
        break;
      case 94:
        descreen = 0.80000000;
        break;
    }
  }
  else if (solvent == "aniline") {
    switch (atmnr) {
      case 1:
        descreen = 0.63606235;
        break;
      case 2:
        descreen = 1.90461008;
        break;
      case 3:
        descreen = 0.30000000;
        break;
      case 4:
        descreen = 0.66432605;
        break;
      case 5:
        descreen = 0.86557022;
        break;
      case 6:
        descreen = 0.33655338;
        break;
      case 7:
        descreen = 0.81296684;
        break;
      case 8:
        descreen = 0.04519102;
        break;
      case 9:
        descreen = 0.80446803;
        break;
      case 10:
        descreen = 0.30000000;
        break;
      case 11:
        descreen = 0.09880557;
        break;
      case 12:
        descreen = 1.03423058;
        break;
      case 13:
        descreen = 1.03328171;
        break;
      case 14:
        descreen = 0.88613511;
        break;
      case 15:
        descreen = 0.36689316;
        break;
      case 16:
        descreen = 0.83271482;
        break;
      case 17:
        descreen = 1.12754750;
        break;
      case 18:
        descreen = 0.70866652;
        break;
      case 19:
        descreen = 0.30000000;
        break;
      case 20:
        descreen = 0.30000000;
        break;
      case 21:
        descreen = 1.07299680;
        break;
      case 22:
        descreen = 1.06506047;
        break;
      case 23:
        descreen = 1.18796933;
        break;
      case 24:
        descreen = 1.10493342;
        break;
      case 25:
        descreen = 1.23555482;
        break;
      case 26:
        descreen = 1.34459642;
        break;
      case 27:
        descreen = 0.30000000;
        break;
      case 28:
        descreen = 0.86958534;
        break;
      case 29:
        descreen = 0.95178949;
        break;
      case 30:
        descreen = 0.69345353;
        break;
      case 31:
        descreen = 1.14302096;
        break;
      case 32:
        descreen = 0.74926530;
        break;
      case 33:
        descreen = 1.12796278;
        break;
      case 34:
        descreen = 0.87516562;
        break;
      case 35:
        descreen = 0.66418114;
        break;
      case 36:
        descreen = 1.19372266;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 0.48326820;
        break;
      case 39:
        descreen = 1.05595120;
        break;
      case 40:
        descreen = 0.65625684;
        break;
      case 41:
        descreen = 1.22915663;
        break;
      case 42:
        descreen = 0.81502057;
        break;
      case 43:
        descreen = 1.17973220;
        break;
      case 44:
        descreen = 1.49439343;
        break;
      case 45:
        descreen = 1.26260832;
        break;
      case 46:
        descreen = 1.13168322;
        break;
      case 47:
        descreen = 0.98389441;
        break;
      case 48:
        descreen = 1.09193368;
        break;
      case 49:
        descreen = 0.82405083;
        break;
      case 50:
        descreen = 0.30000000;
        break;
      case 51:
        descreen = 1.12708752;
        break;
      case 52:
        descreen = 0.66720759;
        break;
      case 53:
        descreen = 0.74844322;
        break;
      case 54:
        descreen = 1.00653653;
        break;
      case 55:
        descreen = 0.26744866;
        break;
      case 56:
        descreen = 0.91697745;
        break;
      case 57:
        descreen = 1.00000000;
        break;
      case 58:
        descreen = 1.00000000;
        break;
      case 59:
        descreen = 1.00000000;
        break;
      case 60:
        descreen = 1.00000000;
        break;
      case 61:
        descreen = 1.00000000;
        break;
      case 62:
        descreen = 1.00000000;
        break;
      case 63:
        descreen = 1.00000000;
        break;
      case 64:
        descreen = 1.00000000;
        break;
      case 65:
        descreen = 1.00000000;
        break;
      case 66:
        descreen = 1.00000000;
        break;
      case 67:
        descreen = 1.00000000;
        break;
      case 68:
        descreen = 1.00000000;
        break;
      case 69:
        descreen = 1.00000000;
        break;
      case 70:
        descreen = 1.00000000;
        break;
      case 71:
        descreen = 1.00000000;
        break;
      case 72:
        descreen = 1.26533804;
        break;
      case 73:
        descreen = 1.09633120;
        break;
      case 74:
        descreen = 0.78905670;
        break;
      case 75:
        descreen = 1.29186761;
        break;
      case 76:
        descreen = 1.22830941;
        break;
      case 77:
        descreen = 1.13900129;
        break;
      case 78:
        descreen = 1.17881856;
        break;
      case 79:
        descreen = 1.02488217;
        break;
      case 80:
        descreen = 0.30000000;
        break;
      case 81:
        descreen = 0.57825236;
        break;
      case 82:
        descreen = 0.99096833;
        break;
      case 83:
        descreen = 0.54120328;
        break;
      case 84:
        descreen = 1.00000000;
        break;
      case 85:
        descreen = 1.00000000;
        break;
      case 86:
        descreen = 1.00000000;
        break;
      case 87:
        descreen = 1.00000000;
        break;
      case 88:
        descreen = 1.00000000;
        break;
      case 89:
        descreen = 1.00000000;
        break;
      case 90:
        descreen = 1.00000000;
        break;
      case 91:
        descreen = 1.00000000;
        break;
      case 92:
        descreen = 1.00000000;
        break;
      case 93:
        descreen = 1.00000000;
        break;
      case 94:
        descreen = 1.00000000;
        break;
    }
  }
  else if (solvent == "benzaldehyde") {
    switch (atmnr) {
      case 1:
        descreen = 0.30000000;
        break;
      case 2:
        descreen = 1.61346257;
        break;
      case 3:
        descreen = 0.30000000;
        break;
      case 4:
        descreen = 0.72329908;
        break;
      case 5:
        descreen = 0.90415406;
        break;
      case 6:
        descreen = 0.88650074;
        break;
      case 7:
        descreen = 0.30000000;
        break;
      case 8:
        descreen = 0.30000000;
        break;
      case 9:
        descreen = 0.62508407;
        break;
      case 10:
        descreen = 0.30000000;
        break;
      case 11:
        descreen = 1.39661904;
        break;
      case 12:
        descreen = 0.40000272;
        break;
      case 13:
        descreen = 1.01280579;
        break;
      case 14:
        descreen = 0.87840431;
        break;
      case 15:
        descreen = 0.72565292;
        break;
      case 16:
        descreen = 0.99589149;
        break;
      case 17:
        descreen = 0.69936160;
        break;
      case 18:
        descreen = 0.30000000;
        break;
      case 19:
        descreen = 0.30000000;
        break;
      case 20:
        descreen = 0.30000000;
        break;
      case 21:
        descreen = 1.01819363;
        break;
      case 22:
        descreen = 1.08326840;
        break;
      case 23:
        descreen = 1.12282316;
        break;
      case 24:
        descreen = 1.31387356;
        break;
      case 25:
        descreen = 1.03269268;
        break;
      case 26:
        descreen = 1.39788491;
        break;
      case 27:
        descreen = 0.30000000;
        break;
      case 28:
        descreen = 0.30000000;
        break;
      case 29:
        descreen = 0.30000000;
        break;
      case 30:
        descreen = 0.51166801;
        break;
      case 31:
        descreen = 1.04118715;
        break;
      case 32:
        descreen = 0.66053269;
        break;
      case 33:
        descreen = 0.68451764;
        break;
      case 34:
        descreen = 1.18173725;
        break;
      case 35:
        descreen = 0.86780744;
        break;
      case 36:
        descreen = 0.91189787;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 0.30000000;
        break;
      case 39:
        descreen = 1.08094157;
        break;
      case 40:
        descreen = 0.65297835;
        break;
      case 41:
        descreen = 1.14126726;
        break;
      case 42:
        descreen = 1.12987591;
        break;
      case 43:
        descreen = 1.29974033;
        break;
      case 44:
        descreen = 1.30361238;
        break;
      case 45:
        descreen = 1.35870404;
        break;
      case 46:
        descreen = 0.30000000;
        break;
      case 47:
        descreen = 0.83581563;
        break;
      case 48:
        descreen = 0.72159784;
        break;
      case 49:
        descreen = 0.86859167;
        break;
      case 50:
        descreen = 0.30000000;
        break;
      case 51:
        descreen = 0.79322763;
        break;
      case 52:
        descreen = 0.95456497;
        break;
      case 53:
        descreen = 1.01152475;
        break;
      case 54:
        descreen = 1.01695596;
        break;
      case 55:
        descreen = 1.00001606;
        break;
      case 56:
        descreen = 0.94824319;
        break;
      case 57:
        descreen = 1.00000000;
        break;
      case 58:
        descreen = 1.00000000;
        break;
      case 59:
        descreen = 1.00000000;
        break;
      case 60:
        descreen = 1.00000000;
        break;
      case 61:
        descreen = 1.00000000;
        break;
      case 62:
        descreen = 1.00000000;
        break;
      case 63:
        descreen = 1.00000000;
        break;
      case 64:
        descreen = 1.00000000;
        break;
      case 65:
        descreen = 1.00000000;
        break;
      case 66:
        descreen = 1.00000000;
        break;
      case 67:
        descreen = 1.00000000;
        break;
      case 68:
        descreen = 1.00000000;
        break;
      case 69:
        descreen = 1.00000000;
        break;
      case 70:
        descreen = 1.00000000;
        break;
      case 71:
        descreen = 1.00000000;
        break;
      case 72:
        descreen = 1.20600000;
        break;
      case 73:
        descreen = 1.03370566;
        break;
      case 74:
        descreen = 0.09474026;
        break;
      case 75:
        descreen = 1.03678522;
        break;
      case 76:
        descreen = 1.28943937;
        break;
      case 77:
        descreen = 1.19776428;
        break;
      case 78:
        descreen = 0.96819305;
        break;
      case 79:
        descreen = 0.79640895;
        break;
      case 80:
        descreen = 1.03355029;
        break;
      case 81:
        descreen = 0.87396085;
        break;
      case 82:
        descreen = 0.98726110;
        break;
      case 83:
        descreen = 0.87780299;
        break;
      case 84:
        descreen = 1.00000000;
        break;
      case 85:
        descreen = 1.00000000;
        break;
      case 86:
        descreen = 1.00000000;
        break;
      case 87:
        descreen = 1.00000000;
        break;
      case 88:
        descreen = 1.00000000;
        break;
      case 89:
        descreen = 1.00000000;
        break;
      case 90:
        descreen = 1.00000000;
        break;
      case 91:
        descreen = 1.00000000;
        break;
      case 92:
        descreen = 1.00000000;
        break;
      case 93:
        descreen = 1.00000000;
        break;
      case 94:
        descreen = 1.00000000;
        break;
    }
  }
  else if (solvent == "benzene") {
    switch (atmnr) {
      case 1:
        descreen = 0.30000000;
        break;
      case 2:
        descreen = 1.83908723;
        break;
      case 3:
        descreen = 0.08101223;
        break;
      case 4:
        descreen = 0.85405447;
        break;
      case 5:
        descreen = 0.79885844;
        break;
      case 6:
        descreen = 0.89956631;
        break;
      case 7:
        descreen = 0.69373209;
        break;
      case 8:
        descreen = 0.97484683;
        break;
      case 9:
        descreen = 0.67756560;
        break;
      case 10:
        descreen = 1.60107651;
        break;
      case 11:
        descreen = 0.30000000;
        break;
      case 12:
        descreen = 1.02908641;
        break;
      case 13:
        descreen = 0.76635684;
        break;
      case 14:
        descreen = 0.91624669;
        break;
      case 15:
        descreen = 0.83857854;
        break;
      case 16:
        descreen = 0.92415015;
        break;
      case 17:
        descreen = 0.83876766;
        break;
      case 18:
        descreen = 0.81381344;
        break;
      case 19:
        descreen = 0.30000000;
        break;
      case 20:
        descreen = 0.30000000;
        break;
      case 21:
        descreen = 0.97871615;
        break;
      case 22:
        descreen = 1.04968870;
        break;
      case 23:
        descreen = 1.11852537;
        break;
      case 24:
        descreen = 1.17326074;
        break;
      case 25:
        descreen = 1.08808540;
        break;
      case 26:
        descreen = 1.24646400;
        break;
      case 27:
        descreen = 0.94870717;
        break;
      case 28:
        descreen = 0.82222625;
        break;
      case 29:
        descreen = 0.30000000;
        break;
      case 30:
        descreen = 1.13515462;
        break;
      case 31:
        descreen = 1.11046532;
        break;
      case 32:
        descreen = 0.78852512;
        break;
      case 33:
        descreen = 1.03237961;
        break;
      case 34:
        descreen = 1.21390807;
        break;
      case 35:
        descreen = 0.30000000;
        break;
      case 36:
        descreen = 1.05153339;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 0.30000000;
        break;
      case 39:
        descreen = 0.95554419;
        break;
      case 40:
        descreen = 0.77380944;
        break;
      case 41:
        descreen = 1.15543062;
        break;
      case 42:
        descreen = 1.12734500;
        break;
      case 43:
        descreen = 1.35413983;
        break;
      case 44:
        descreen = 1.25561286;
        break;
      case 45:
        descreen = 1.26194500;
        break;
      case 46:
        descreen = 1.16174826;
        break;
      case 47:
        descreen = 1.01518058;
        break;
      case 48:
        descreen = 0.84597912;
        break;
      case 49:
        descreen = 0.88104467;
        break;
      case 50:
        descreen = 0.30000000;
        break;
      case 51:
        descreen = 1.06112470;
        break;
      case 52:
        descreen = 0.76212761;
        break;
      case 53:
        descreen = 0.59834033;
        break;
      case 54:
        descreen = 0.96053668;
        break;
      case 55:
        descreen = 0.00158803;
        break;
      case 56:
        descreen = 1.01900329;
        break;
      case 57:
        descreen = 0.80000000;
        break;
      case 58:
        descreen = 0.80000000;
        break;
      case 59:
        descreen = 0.80000000;
        break;
      case 60:
        descreen = 0.80000000;
        break;
      case 61:
        descreen = 0.80000000;
        break;
      case 62:
        descreen = 0.80000000;
        break;
      case 63:
        descreen = 0.80000000;
        break;
      case 64:
        descreen = 0.80000000;
        break;
      case 65:
        descreen = 0.80000000;
        break;
      case 66:
        descreen = 0.80000000;
        break;
      case 67:
        descreen = 0.80000000;
        break;
      case 68:
        descreen = 0.80000000;
        break;
      case 69:
        descreen = 0.80000000;
        break;
      case 70:
        descreen = 0.80000000;
        break;
      case 71:
        descreen = 0.80000000;
        break;
      case 72:
        descreen = 1.11068729;
        break;
      case 73:
        descreen = 1.01340753;
        break;
      case 74:
        descreen = 0.05417569;
        break;
      case 75:
        descreen = 1.03834101;
        break;
      case 76:
        descreen = 1.21151920;
        break;
      case 77:
        descreen = 1.14185967;
        break;
      case 78:
        descreen = 1.07046574;
        break;
      case 79:
        descreen = 0.97928478;
        break;
      case 80:
        descreen = 1.26108955;
        break;
      case 81:
        descreen = 0.85241365;
        break;
      case 82:
        descreen = 0.98379651;
        break;
      case 83:
        descreen = 0.69795610;
        break;
      case 84:
        descreen = 0.80000000;
        break;
      case 85:
        descreen = 0.80000000;
        break;
      case 86:
        descreen = 0.80000000;
        break;
      case 87:
        descreen = 0.80000000;
        break;
      case 88:
        descreen = 0.80000000;
        break;
      case 89:
        descreen = 0.80000000;
        break;
      case 90:
        descreen = 0.80000000;
        break;
      case 91:
        descreen = 0.80000000;
        break;
      case 92:
        descreen = 0.80000000;
        break;
      case 93:
        descreen = 0.80000000;
        break;
      case 94:
        descreen = 0.80000000;
        break;
    }
  }
  else if (solvent == "dichloromethane") {
    switch (atmnr) {
      case 1:
        descreen = 0.28926333;
        break;
      case 2:
        descreen = 1.99261093;
        break;
      case 3:
        descreen = 0.30000000;
        break;
      case 4:
        descreen = 0.79793587;
        break;
      case 5:
        descreen = 0.81135795;
        break;
      case 6:
        descreen = 0.83616685;
        break;
      case 7:
        descreen = 0.79057315;
        break;
      case 8:
        descreen = 0.86271642;
        break;
      case 9:
        descreen = 0.70573060;
        break;
      case 10:
        descreen = 1.38322206;
        break;
      case 11:
        descreen = 0.30000000;
        break;
      case 12:
        descreen = 0.30000000;
        break;
      case 13:
        descreen = 0.70041023;
        break;
      case 14:
        descreen = 0.90673216;
        break;
      case 15:
        descreen = 0.76272228;
        break;
      case 16:
        descreen = 0.91021623;
        break;
      case 17:
        descreen = 0.81532007;
        break;
      case 18:
        descreen = 0.53937688;
        break;
      case 19:
        descreen = 0.30000000;
        break;
      case 20:
        descreen = 0.30000000;
        break;
      case 21:
        descreen = 1.01472709;
        break;
      case 22:
        descreen = 1.04010262;
        break;
      case 23:
        descreen = 1.10949252;
        break;
      case 24:
        descreen = 1.33303258;
        break;
      case 25:
        descreen = 1.12658111;
        break;
      case 26:
        descreen = 1.42195776;
        break;
      case 27:
        descreen = 0.80273663;
        break;
      case 28:
        descreen = 0.81746015;
        break;
      case 29:
        descreen = 0.65703514;
        break;
      case 30:
        descreen = 1.18208681;
        break;
      case 31:
        descreen = 1.01508856;
        break;
      case 32:
        descreen = 0.80119384;
        break;
      case 33:
        descreen = 1.05604075;
        break;
      case 34:
        descreen = 1.20897255;
        break;
      case 35:
        descreen = 0.39287452;
        break;
      case 36:
        descreen = 0.96296618;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 0.30000000;
        break;
      case 39:
        descreen = 1.00370492;
        break;
      case 40:
        descreen = 0.09377212;
        break;
      case 41:
        descreen = 1.16652332;
        break;
      case 42:
        descreen = 1.29719519;
        break;
      case 43:
        descreen = 1.27712571;
        break;
      case 44:
        descreen = 1.28018285;
        break;
      case 45:
        descreen = 1.28811425;
        break;
      case 46:
        descreen = 1.10052024;
        break;
      case 47:
        descreen = 0.95578477;
        break;
      case 48:
        descreen = 1.09164864;
        break;
      case 49:
        descreen = 0.84497183;
        break;
      case 50:
        descreen = 0.30000000;
        break;
      case 51:
        descreen = 0.83981245;
        break;
      case 52:
        descreen = 0.88552660;
        break;
      case 53:
        descreen = 0.39084153;
        break;
      case 54:
        descreen = 0.81863824;
        break;
      case 55:
        descreen = 0.30000000;
        break;
      case 56:
        descreen = 1.04090665;
        break;
      case 57:
        descreen = 0.80000000;
        break;
      case 58:
        descreen = 0.80000000;
        break;
      case 59:
        descreen = 0.80000000;
        break;
      case 60:
        descreen = 0.80000000;
        break;
      case 61:
        descreen = 0.80000000;
        break;
      case 62:
        descreen = 0.80000000;
        break;
      case 63:
        descreen = 0.80000000;
        break;
      case 64:
        descreen = 0.80000000;
        break;
      case 65:
        descreen = 0.80000000;
        break;
      case 66:
        descreen = 0.80000000;
        break;
      case 67:
        descreen = 0.80000000;
        break;
      case 68:
        descreen = 0.80000000;
        break;
      case 69:
        descreen = 0.80000000;
        break;
      case 70:
        descreen = 0.80000000;
        break;
      case 71:
        descreen = 0.80000000;
        break;
      case 72:
        descreen = 1.20679413;
        break;
      case 73:
        descreen = 1.01528947;
        break;
      case 74:
        descreen = 0.27374694;
        break;
      case 75:
        descreen = 1.13101163;
        break;
      case 76:
        descreen = 1.24496642;
        break;
      case 77:
        descreen = 1.22853764;
        break;
      case 78:
        descreen = 1.17172086;
        break;
      case 79:
        descreen = 0.96700342;
        break;
      case 80:
        descreen = 0.82371044;
        break;
      case 81:
        descreen = 0.64251626;
        break;
      case 82:
        descreen = 0.99215602;
        break;
      case 83:
        descreen = 0.64849560;
        break;
      case 84:
        descreen = 0.80000000;
        break;
      case 85:
        descreen = 0.80000000;
        break;
      case 86:
        descreen = 0.80000000;
        break;
      case 87:
        descreen = 0.80000000;
        break;
      case 88:
        descreen = 0.80000000;
        break;
      case 89:
        descreen = 0.80000000;
        break;
      case 90:
        descreen = 0.80000000;
        break;
      case 91:
        descreen = 0.80000000;
        break;
      case 92:
        descreen = 0.80000000;
        break;
      case 93:
        descreen = 0.80000000;
        break;
      case 94:
        descreen = 0.80000000;
        break;
    }
  }
  else if (solvent == "chloroform") {
    switch (atmnr) {
      case 1:
        descreen = 0.54425653;
        break;
      case 2:
        descreen = 0.30000000;
        break;
      case 3:
        descreen = 0.68352647;
        break;
      case 4:
        descreen = 0.30024907;
        break;
      case 5:
        descreen = 0.87131168;
        break;
      case 6:
        descreen = 0.44274116;
        break;
      case 7:
        descreen = 0.91118733;
        break;
      case 8:
        descreen = 1.02990206;
        break;
      case 9:
        descreen = 0.87224972;
        break;
      case 10:
        descreen = 1.36064412;
        break;
      case 11:
        descreen = 0.04523267;
        break;
      case 12:
        descreen = 1.01826727;
        break;
      case 13:
        descreen = 0.80291179;
        break;
      case 14:
        descreen = 0.91994934;
        break;
      case 15:
        descreen = 0.27703877;
        break;
      case 16:
        descreen = 0.79671366;
        break;
      case 17:
        descreen = 1.12191178;
        break;
      case 18:
        descreen = 1.06327945;
        break;
      case 19:
        descreen = 0.30000000;
        break;
      case 20:
        descreen = 0.30000000;
        break;
      case 21:
        descreen = 1.05134616;
        break;
      case 22:
        descreen = 1.02524534;
        break;
      case 23:
        descreen = 1.15488259;
        break;
      case 24:
        descreen = 1.15659192;
        break;
      case 25:
        descreen = 1.15582331;
        break;
      case 26:
        descreen = 1.35956403;
        break;
      case 27:
        descreen = 1.02994976;
        break;
      case 28:
        descreen = 0.89578087;
        break;
      case 29:
        descreen = 1.07980658;
        break;
      case 30:
        descreen = 1.18018780;
        break;
      case 31:
        descreen = 1.15034660;
        break;
      case 32:
        descreen = 0.80126086;
        break;
      case 33:
        descreen = 1.06719436;
        break;
      case 34:
        descreen = 0.97305695;
        break;
      case 35:
        descreen = 0.30000000;
        break;
      case 36:
        descreen = 1.17576974;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 1.26859525;
        break;
      case 39:
        descreen = 0.91060526;
        break;
      case 40:
        descreen = 0.30000000;
        break;
      case 41:
        descreen = 1.26737715;
        break;
      case 42:
        descreen = 0.88598699;
        break;
      case 43:
        descreen = 1.28128367;
        break;
      case 44:
        descreen = 1.41019625;
        break;
      case 45:
        descreen = 1.24068510;
        break;
      case 46:
        descreen = 1.22207835;
        break;
      case 47:
        descreen = 1.05317124;
        break;
      case 48:
        descreen = 0.30000000;
        break;
      case 49:
        descreen = 0.78929937;
        break;
      case 50:
        descreen = 0.30000000;
        break;
      case 51:
        descreen = 0.54876707;
        break;
      case 52:
        descreen = 0.85496180;
        break;
      case 53:
        descreen = 0.30000000;
        break;
      case 54:
        descreen = 0.90436962;
        break;
      case 55:
        descreen = 1.11376639;
        break;
      case 56:
        descreen = 1.08812332;
        break;
      case 57:
        descreen = 0.80000000;
        break;
      case 58:
        descreen = 0.80000000;
        break;
      case 59:
        descreen = 0.80000000;
        break;
      case 60:
        descreen = 0.80000000;
        break;
      case 61:
        descreen = 0.80000000;
        break;
      case 62:
        descreen = 0.80000000;
        break;
      case 63:
        descreen = 0.80000000;
        break;
      case 64:
        descreen = 0.80000000;
        break;
      case 65:
        descreen = 0.80000000;
        break;
      case 66:
        descreen = 0.80000000;
        break;
      case 67:
        descreen = 0.80000000;
        break;
      case 68:
        descreen = 0.80000000;
        break;
      case 69:
        descreen = 0.80000000;
        break;
      case 70:
        descreen = 0.80000000;
        break;
      case 71:
        descreen = 0.80000000;
        break;
      case 72:
        descreen = 1.22592678;
        break;
      case 73:
        descreen = 1.08931843;
        break;
      case 74:
        descreen = 0.62355860;
        break;
      case 75:
        descreen = 1.25551232;
        break;
      case 76:
        descreen = 1.28912417;
        break;
      case 77:
        descreen = 1.37088359;
        break;
      case 78:
        descreen = 1.21255461;
        break;
      case 79:
        descreen = 1.05341027;
        break;
      case 80:
        descreen = 0.74644294;
        break;
      case 81:
        descreen = 0.54998313;
        break;
      case 82:
        descreen = 0.76413109;
        break;
      case 83:
        descreen = 1.12560117;
        break;
      case 84:
        descreen = 0.80000000;
        break;
      case 85:
        descreen = 0.80000000;
        break;
      case 86:
        descreen = 0.80000000;
        break;
      case 87:
        descreen = 0.80000000;
        break;
      case 88:
        descreen = 0.80000000;
        break;
      case 89:
        descreen = 0.80000000;
        break;
      case 90:
        descreen = 0.80000000;
        break;
      case 91:
        descreen = 0.80000000;
        break;
      case 92:
        descreen = 0.80000000;
        break;
      case 93:
        descreen = 0.80000000;
        break;
      case 94:
        descreen = 0.80000000;
        break;
    }
  }
  else if (solvent == "carbon disulfide") {
    switch (atmnr) {
      case 1:
        descreen = 0.67101392;
        break;
      case 2:
        descreen = 1.83162475;
        break;
      case 3:
        descreen = 0.74845704;
        break;
      case 4:
        descreen = 0.81714736;
        break;
      case 5:
        descreen = 0.91247423;
        break;
      case 6:
        descreen = 0.89039802;
        break;
      case 7:
        descreen = 0.34641431;
        break;
      case 8:
        descreen = 1.00584787;
        break;
      case 9:
        descreen = 0.81852179;
        break;
      case 10:
        descreen = 0.30000000;
        break;
      case 11:
        descreen = 0.30000000;
        break;
      case 12:
        descreen = 0.98259438;
        break;
      case 13:
        descreen = 0.76203165;
        break;
      case 14:
        descreen = 0.95238120;
        break;
      case 15:
        descreen = 0.64489838;
        break;
      case 16:
        descreen = 0.67618683;
        break;
      case 17:
        descreen = 0.91744526;
        break;
      case 18:
        descreen = 1.31005168;
        break;
      case 19:
        descreen = 0.59274765;
        break;
      case 20:
        descreen = 0.30000000;
        break;
      case 21:
        descreen = 0.93866362;
        break;
      case 22:
        descreen = 1.04942980;
        break;
      case 23:
        descreen = 1.15552558;
        break;
      case 24:
        descreen = 1.19633812;
        break;
      case 25:
        descreen = 1.09018770;
        break;
      case 26:
        descreen = 1.33937956;
        break;
      case 27:
        descreen = 1.29141196;
        break;
      case 28:
        descreen = 0.97055753;
        break;
      case 29:
        descreen = 0.88131074;
        break;
      case 30:
        descreen = 1.16378129;
        break;
      case 31:
        descreen = 0.95344329;
        break;
      case 32:
        descreen = 0.86373088;
        break;
      case 33:
        descreen = 1.07656363;
        break;
      case 34:
        descreen = 1.25642469;
        break;
      case 35:
        descreen = 0.30000000;
        break;
      case 36:
        descreen = 1.22913065;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 0.30000000;
        break;
      case 39:
        descreen = 0.90701663;
        break;
      case 40:
        descreen = 0.04955879;
        break;
      case 41:
        descreen = 1.08947296;
        break;
      case 42:
        descreen = 0.98903611;
        break;
      case 43:
        descreen = 1.52634108;
        break;
      case 44:
        descreen = 1.33495934;
        break;
      case 45:
        descreen = 1.26075308;
        break;
      case 46:
        descreen = 1.25921485;
        break;
      case 47:
        descreen = 1.28794815;
        break;
      case 48:
        descreen = 1.08985190;
        break;
      case 49:
        descreen = 0.90115381;
        break;
      case 50:
        descreen = 0.67061995;
        break;
      case 51:
        descreen = 1.08462735;
        break;
      case 52:
        descreen = 0.74799904;
        break;
      case 53:
        descreen = 0.30000000;
        break;
      case 54:
        descreen = 0.84506316;
        break;
      case 55:
        descreen = 0.68945810;
        break;
      case 56:
        descreen = 1.08066679;
        break;
      case 57:
        descreen = 0.80000000;
        break;
      case 58:
        descreen = 0.80000000;
        break;
      case 59:
        descreen = 0.80000000;
        break;
      case 60:
        descreen = 0.80000000;
        break;
      case 61:
        descreen = 0.80000000;
        break;
      case 62:
        descreen = 0.80000000;
        break;
      case 63:
        descreen = 0.80000000;
        break;
      case 64:
        descreen = 0.80000000;
        break;
      case 65:
        descreen = 0.80000000;
        break;
      case 66:
        descreen = 0.80000000;
        break;
      case 67:
        descreen = 0.80000000;
        break;
      case 68:
        descreen = 0.80000000;
        break;
      case 69:
        descreen = 0.80000000;
        break;
      case 70:
        descreen = 0.80000000;
        break;
      case 71:
        descreen = 0.80000000;
        break;
      case 72:
        descreen = 1.10139934;
        break;
      case 73:
        descreen = 0.91376504;
        break;
      case 74:
        descreen = 0.82801025;
        break;
      case 75:
        descreen = 1.22830578;
        break;
      case 76:
        descreen = 1.24664000;
        break;
      case 77:
        descreen = 1.23171074;
        break;
      case 78:
        descreen = 1.19796332;
        break;
      case 79:
        descreen = 1.20890185;
        break;
      case 80:
        descreen = 1.29655080;
        break;
      case 81:
        descreen = 0.83370444;
        break;
      case 82:
        descreen = 1.00609132;
        break;
      case 83:
        descreen = 1.11767725;
        break;
      case 84:
        descreen = 0.80000000;
        break;
      case 85:
        descreen = 0.80000000;
        break;
      case 86:
        descreen = 0.80000000;
        break;
      case 87:
        descreen = 0.80000000;
        break;
      case 88:
        descreen = 0.80000000;
        break;
      case 89:
        descreen = 0.80000000;
        break;
      case 90:
        descreen = 0.80000000;
        break;
      case 91:
        descreen = 0.80000000;
        break;
      case 92:
        descreen = 0.80000000;
        break;
      case 93:
        descreen = 0.80000000;
        break;
      case 94:
        descreen = 0.80000000;
        break;
    }
  }
  else if (solvent == "dioxane") {
    switch (atmnr) {
      case 1:
        descreen = 0.99501947;
        break;
      case 2:
        descreen = 0.54939611;
        break;
      case 3:
        descreen = 0.84749233;
        break;
      case 4:
        descreen = 0.08192417;
        break;
      case 5:
        descreen = 0.92829185;
        break;
      case 6:
        descreen = 0.30000000;
        break;
      case 7:
        descreen = 0.30000000;
        break;
      case 8:
        descreen = 0.30000000;
        break;
      case 9:
        descreen = 1.02157731;
        break;
      case 10:
        descreen = 0.30000000;
        break;
      case 11:
        descreen = 0.50727977;
        break;
      case 12:
        descreen = 0.66635680;
        break;
      case 13:
        descreen = 1.06317195;
        break;
      case 14:
        descreen = 0.90780637;
        break;
      case 15:
        descreen = 0.30000000;
        break;
      case 16:
        descreen = 0.86175540;
        break;
      case 17:
        descreen = 0.93435708;
        break;
      case 18:
        descreen = 1.39628743;
        break;
      case 19:
        descreen = 0.75104383;
        break;
      case 20:
        descreen = 0.30000000;
        break;
      case 21:
        descreen = 1.05628751;
        break;
      case 22:
        descreen = 1.08591184;
        break;
      case 23:
        descreen = 1.13615426;
        break;
      case 24:
        descreen = 0.94969656;
        break;
      case 25:
        descreen = 1.19054921;
        break;
      case 26:
        descreen = 1.38854191;
        break;
      case 27:
        descreen = 0.95877392;
        break;
      case 28:
        descreen = 1.16632036;
        break;
      case 29:
        descreen = 1.04862692;
        break;
      case 30:
        descreen = 0.30000000;
        break;
      case 31:
        descreen = 0.98417023;
        break;
      case 32:
        descreen = 0.78173613;
        break;
      case 33:
        descreen = 0.59487457;
        break;
      case 34:
        descreen = 0.80927909;
        break;
      case 35:
        descreen = 0.84402819;
        break;
      case 36:
        descreen = 1.43968875;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 0.30000000;
        break;
      case 39:
        descreen = 1.24646806;
        break;
      case 40:
        descreen = 0.91395622;
        break;
      case 41:
        descreen = 1.10166853;
        break;
      case 42:
        descreen = 0.53876725;
        break;
      case 43:
        descreen = 0.63549151;
        break;
      case 44:
        descreen = 1.50683885;
        break;
      case 45:
        descreen = 1.22017940;
        break;
      case 46:
        descreen = 1.10791748;
        break;
      case 47:
        descreen = 1.13222922;
        break;
      case 48:
        descreen = 0.61488695;
        break;
      case 49:
        descreen = 0.83615148;
        break;
      case 50:
        descreen = 0.48577020;
        break;
      case 51:
        descreen = 0.53599960;
        break;
      case 52:
        descreen = 0.61714433;
        break;
      case 53:
        descreen = 0.62668099;
        break;
      case 54:
        descreen = 1.15581020;
        break;
      case 55:
        descreen = 0.99978295;
        break;
      case 56:
        descreen = 0.89248506;
        break;
      case 57:
        descreen = 1.00000000;
        break;
      case 58:
        descreen = 1.00000000;
        break;
      case 59:
        descreen = 1.00000000;
        break;
      case 60:
        descreen = 1.00000000;
        break;
      case 61:
        descreen = 1.00000000;
        break;
      case 62:
        descreen = 1.00000000;
        break;
      case 63:
        descreen = 1.00000000;
        break;
      case 64:
        descreen = 1.00000000;
        break;
      case 65:
        descreen = 1.00000000;
        break;
      case 66:
        descreen = 1.00000000;
        break;
      case 67:
        descreen = 1.00000000;
        break;
      case 68:
        descreen = 1.00000000;
        break;
      case 69:
        descreen = 1.00000000;
        break;
      case 70:
        descreen = 1.00000000;
        break;
      case 71:
        descreen = 1.00000000;
        break;
      case 72:
        descreen = 0.30000000;
        break;
      case 73:
        descreen = 1.09239868;
        break;
      case 74:
        descreen = 0.85801984;
        break;
      case 75:
        descreen = 1.22312847;
        break;
      case 76:
        descreen = 0.85940608;
        break;
      case 77:
        descreen = 0.93453487;
        break;
      case 78:
        descreen = 0.58858953;
        break;
      case 79:
        descreen = 1.12905516;
        break;
      case 80:
        descreen = 0.81905399;
        break;
      case 81:
        descreen = 0.62548456;
        break;
      case 82:
        descreen = 1.03207940;
        break;
      case 83:
        descreen = 0.86738337;
        break;
      case 84:
        descreen = 1.00000000;
        break;
      case 85:
        descreen = 1.00000000;
        break;
      case 86:
        descreen = 1.00000000;
        break;
      case 87:
        descreen = 1.00000000;
        break;
      case 88:
        descreen = 1.00000000;
        break;
      case 89:
        descreen = 1.00000000;
        break;
      case 90:
        descreen = 1.00000000;
        break;
      case 91:
        descreen = 1.00000000;
        break;
      case 92:
        descreen = 1.00000000;
        break;
      case 93:
        descreen = 1.00000000;
        break;
      case 94:
        descreen = 1.00000000;
        break;
    }
  }
  else if (solvent == "dmf") {
    switch (atmnr) {
      case 1:
        descreen = 0.30000000;
        break;
      case 2:
        descreen = 1.11740044;
        break;
      case 3:
        descreen = 0.30000000;
        break;
      case 4:
        descreen = 0.79377214;
        break;
      case 5:
        descreen = 0.79128875;
        break;
      case 6:
        descreen = 0.84746424;
        break;
      case 7:
        descreen = 0.30000000;
        break;
      case 8:
        descreen = 0.30000000;
        break;
      case 9:
        descreen = 0.73272246;
        break;
      case 10:
        descreen = 0.30000000;
        break;
      case 11:
        descreen = 0.30000000;
        break;
      case 12:
        descreen = 1.09781017;
        break;
      case 13:
        descreen = 0.68104256;
        break;
      case 14:
        descreen = 0.87613172;
        break;
      case 15:
        descreen = 0.71094960;
        break;
      case 16:
        descreen = 0.97250194;
        break;
      case 17:
        descreen = 0.69149551;
        break;
      case 18:
        descreen = 0.30000000;
        break;
      case 19:
        descreen = 0.81267529;
        break;
      case 20:
        descreen = 0.30000000;
        break;
      case 21:
        descreen = 1.02960081;
        break;
      case 22:
        descreen = 1.06755215;
        break;
      case 23:
        descreen = 1.19817888;
        break;
      case 24:
        descreen = 1.37244147;
        break;
      case 25:
        descreen = 1.13028964;
        break;
      case 26:
        descreen = 1.42931164;
        break;
      case 27:
        descreen = 0.74089201;
        break;
      case 28:
        descreen = 0.87324519;
        break;
      case 29:
        descreen = 0.48264953;
        break;
      case 30:
        descreen = 0.30000000;
        break;
      case 31:
        descreen = 1.10948488;
        break;
      case 32:
        descreen = 0.84166902;
        break;
      case 33:
        descreen = 0.94321119;
        break;
      case 34:
        descreen = 1.20887165;
        break;
      case 35:
        descreen = 0.62809355;
        break;
      case 36:
        descreen = 1.20375255;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 0.30000000;
        break;
      case 39:
        descreen = 1.24805370;
        break;
      case 40:
        descreen = 0.78521635;
        break;
      case 41:
        descreen = 1.21838310;
        break;
      case 42:
        descreen = 1.31736493;
        break;
      case 43:
        descreen = 1.29076324;
        break;
      case 44:
        descreen = 1.34483110;
        break;
      case 45:
        descreen = 1.17454078;
        break;
      case 46:
        descreen = 0.44931579;
        break;
      case 47:
        descreen = 0.95950423;
        break;
      case 48:
        descreen = 0.00103496;
        break;
      case 49:
        descreen = 0.86421402;
        break;
      case 50:
        descreen = 0.62028336;
        break;
      case 51:
        descreen = 0.96766772;
        break;
      case 52:
        descreen = 0.30000000;
        break;
      case 53:
        descreen = 0.65792997;
        break;
      case 54:
        descreen = 1.15170978;
        break;
      case 55:
        descreen = 0.75961623;
        break;
      case 56:
        descreen = 0.96980166;
        break;
      case 57:
        descreen = 0.80000000;
        break;
      case 58:
        descreen = 0.80000000;
        break;
      case 59:
        descreen = 0.80000000;
        break;
      case 60:
        descreen = 0.80000000;
        break;
      case 61:
        descreen = 0.80000000;
        break;
      case 62:
        descreen = 0.80000000;
        break;
      case 63:
        descreen = 0.80000000;
        break;
      case 64:
        descreen = 0.80000000;
        break;
      case 65:
        descreen = 0.80000000;
        break;
      case 66:
        descreen = 0.80000000;
        break;
      case 67:
        descreen = 0.80000000;
        break;
      case 68:
        descreen = 0.80000000;
        break;
      case 69:
        descreen = 0.80000000;
        break;
      case 70:
        descreen = 0.80000000;
        break;
      case 71:
        descreen = 0.80000000;
        break;
      case 72:
        descreen = 1.30020823;
        break;
      case 73:
        descreen = 1.06635336;
        break;
      case 74:
        descreen = 1.12536532;
        break;
      case 75:
        descreen = 1.04133374;
        break;
      case 76:
        descreen = 1.28582930;
        break;
      case 77:
        descreen = 1.22444454;
        break;
      case 78:
        descreen = 0.59366301;
        break;
      case 79:
        descreen = 1.01406884;
        break;
      case 80:
        descreen = 0.80000000;
        break;
      case 81:
        descreen = 0.85197539;
        break;
      case 82:
        descreen = 0.86646859;
        break;
      case 83:
        descreen = 0.55870388;
        break;
      case 84:
        descreen = 0.80000000;
        break;
      case 85:
        descreen = 0.80000000;
        break;
      case 86:
        descreen = 0.80000000;
        break;
      case 87:
        descreen = 0.80000000;
        break;
      case 88:
        descreen = 0.80000000;
        break;
      case 89:
        descreen = 0.80000000;
        break;
      case 90:
        descreen = 0.80000000;
        break;
      case 91:
        descreen = 0.80000000;
        break;
      case 92:
        descreen = 0.80000000;
        break;
      case 93:
        descreen = 0.80000000;
        break;
      case 94:
        descreen = 0.80000000;
        break;
    }
  }
  else if (solvent == "dmso") {
    switch (atmnr) {
      case 1:
        descreen = 0.30000000;
        break;
      case 2:
        descreen = 1.15182062;
        break;
      case 3:
        descreen = 0.36346925;
        break;
      case 4:
        descreen = 0.80891240;
        break;
      case 5:
        descreen = 0.79863748;
        break;
      case 6:
        descreen = 0.84683070;
        break;
      case 7:
        descreen = 0.30000000;
        break;
      case 8:
        descreen = 0.30000000;
        break;
      case 9:
        descreen = 0.79523149;
        break;
      case 10:
        descreen = 0.30000000;
        break;
      case 11:
        descreen = 0.30000000;
        break;
      case 12:
        descreen = 1.11733780;
        break;
      case 13:
        descreen = 0.68855764;
        break;
      case 14:
        descreen = 0.85275203;
        break;
      case 15:
        descreen = 0.64999526;
        break;
      case 16:
        descreen = 0.94531218;
        break;
      case 17:
        descreen = 0.71461583;
        break;
      case 18:
        descreen = 0.06447625;
        break;
      case 19:
        descreen = 0.30000000;
        break;
      case 20:
        descreen = 0.30000000;
        break;
      case 21:
        descreen = 1.04052720;
        break;
      case 22:
        descreen = 1.06355128;
        break;
      case 23:
        descreen = 1.19997332;
        break;
      case 24:
        descreen = 1.37553047;
        break;
      case 25:
        descreen = 1.12749336;
        break;
      case 26:
        descreen = 1.44031207;
        break;
      case 27:
        descreen = 0.73244215;
        break;
      case 28:
        descreen = 0.91336675;
        break;
      case 29:
        descreen = 0.55665998;
        break;
      case 30:
        descreen = 0.30000000;
        break;
      case 31:
        descreen = 1.00566608;
        break;
      case 32:
        descreen = 0.86585707;
        break;
      case 33:
        descreen = 0.92570177;
        break;
      case 34:
        descreen = 1.19682566;
        break;
      case 35:
        descreen = 0.70301713;
        break;
      case 36:
        descreen = 1.22450147;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 0.30000000;
        break;
      case 39:
        descreen = 1.27558039;
        break;
      case 40:
        descreen = 0.86739237;
        break;
      case 41:
        descreen = 1.21827824;
        break;
      case 42:
        descreen = 1.31824236;
        break;
      case 43:
        descreen = 1.28151078;
        break;
      case 44:
        descreen = 1.34103395;
        break;
      case 45:
        descreen = 1.15206036;
        break;
      case 46:
        descreen = 0.53212940;
        break;
      case 47:
        descreen = 0.95490440;
        break;
      case 48:
        descreen = 0.30000000;
        break;
      case 49:
        descreen = 0.86177123;
        break;
      case 50:
        descreen = 0.64518231;
        break;
      case 51:
        descreen = 0.87758629;
        break;
      case 52:
        descreen = 0.30000000;
        break;
      case 53:
        descreen = 0.58562896;
        break;
      case 54:
        descreen = 1.15408952;
        break;
      case 55:
        descreen = 0.30000000;
        break;
      case 56:
        descreen = 0.97793465;
        break;
      case 57:
        descreen = 0.80000000;
        break;
      case 58:
        descreen = 0.80000000;
        break;
      case 59:
        descreen = 0.80000000;
        break;
      case 60:
        descreen = 0.80000000;
        break;
      case 61:
        descreen = 0.80000000;
        break;
      case 62:
        descreen = 0.80000000;
        break;
      case 63:
        descreen = 0.80000000;
        break;
      case 64:
        descreen = 0.80000000;
        break;
      case 65:
        descreen = 0.80000000;
        break;
      case 66:
        descreen = 0.80000000;
        break;
      case 67:
        descreen = 0.80000000;
        break;
      case 68:
        descreen = 0.80000000;
        break;
      case 69:
        descreen = 0.80000000;
        break;
      case 70:
        descreen = 0.80000000;
        break;
      case 71:
        descreen = 0.80000000;
        break;
      case 72:
        descreen = 1.29747497;
        break;
      case 73:
        descreen = 1.05674546;
        break;
      case 74:
        descreen = 1.13719614;
        break;
      case 75:
        descreen = 1.07832899;
        break;
      case 76:
        descreen = 1.28825126;
        break;
      case 77:
        descreen = 1.22762234;
        break;
      case 78:
        descreen = 0.57649023;
        break;
      case 79:
        descreen = 1.01284711;
        break;
      case 80:
        descreen = 0.80000000;
        break;
      case 81:
        descreen = 0.81705969;
        break;
      case 82:
        descreen = 0.98935435;
        break;
      case 83:
        descreen = 0.47102566;
        break;
      case 84:
        descreen = 0.80000000;
        break;
      case 85:
        descreen = 0.80000000;
        break;
      case 86:
        descreen = 0.80000000;
        break;
      case 87:
        descreen = 0.80000000;
        break;
      case 88:
        descreen = 0.80000000;
        break;
      case 89:
        descreen = 0.80000000;
        break;
      case 90:
        descreen = 0.80000000;
        break;
      case 91:
        descreen = 0.80000000;
        break;
      case 92:
        descreen = 0.80000000;
        break;
      case 93:
        descreen = 0.80000000;
        break;
      case 94:
        descreen = 0.80000000;
        break;
    }
  }
  else if (solvent == "ethanol") {
    switch (atmnr) {
      case 1:
        descreen = 0.67300950;
        break;
      case 2:
        descreen = 1.48082514;
        break;
      case 3:
        descreen = 0.30000000;
        break;
      case 4:
        descreen = 0.30000000;
        break;
      case 5:
        descreen = 0.93890156;
        break;
      case 6:
        descreen = 0.66599452;
        break;
      case 7:
        descreen = 0.64141786;
        break;
      case 8:
        descreen = 0.30000000;
        break;
      case 9:
        descreen = 0.58410696;
        break;
      case 10:
        descreen = 0.10806216;
        break;
      case 11:
        descreen = 0.30000000;
        break;
      case 12:
        descreen = 0.30000000;
        break;
      case 13:
        descreen = 0.69868374;
        break;
      case 14:
        descreen = 0.88202513;
        break;
      case 15:
        descreen = 0.56538613;
        break;
      case 16:
        descreen = 0.30000000;
        break;
      case 17:
        descreen = 0.67076347;
        break;
      case 18:
        descreen = 0.30000000;
        break;
      case 19:
        descreen = 0.30000000;
        break;
      case 20:
        descreen = 1.06221632;
        break;
      case 21:
        descreen = 1.03845566;
        break;
      case 22:
        descreen = 1.03686698;
        break;
      case 23:
        descreen = 1.18601149;
        break;
      case 24:
        descreen = 1.14819850;
        break;
      case 25:
        descreen = 1.19382065;
        break;
      case 26:
        descreen = 1.27807253;
        break;
      case 27:
        descreen = 1.10136260;
        break;
      case 28:
        descreen = 1.10327617;
        break;
      case 29:
        descreen = 1.00338308;
        break;
      case 30:
        descreen = 0.30000000;
        break;
      case 31:
        descreen = 1.13425281;
        break;
      case 32:
        descreen = 0.75272852;
        break;
      case 33:
        descreen = 1.12004575;
        break;
      case 34:
        descreen = 1.00444769;
        break;
      case 35:
        descreen = 0.61819443;
        break;
      case 36:
        descreen = 1.02837404;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 1.01897022;
        break;
      case 39:
        descreen = 1.40172230;
        break;
      case 40:
        descreen = 0.50606322;
        break;
      case 41:
        descreen = 1.23472719;
        break;
      case 42:
        descreen = 0.01428596;
        break;
      case 43:
        descreen = 1.30913755;
        break;
      case 44:
        descreen = 1.40802475;
        break;
      case 45:
        descreen = 1.24326656;
        break;
      case 46:
        descreen = 1.22477940;
        break;
      case 47:
        descreen = 1.27877712;
        break;
      case 48:
        descreen = 0.67969962;
        break;
      case 49:
        descreen = 0.81985787;
        break;
      case 50:
        descreen = 0.30000000;
        break;
      case 51:
        descreen = 0.52373132;
        break;
      case 52:
        descreen = 0.02799162;
        break;
      case 53:
        descreen = 1.47551118;
        break;
      case 54:
        descreen = 0.91807593;
        break;
      case 55:
        descreen = 0.30000000;
        break;
      case 56:
        descreen = 0.03606160;
        break;
      case 57:
        descreen = 0.80000000;
        break;
      case 58:
        descreen = 0.80000000;
        break;
      case 59:
        descreen = 0.80000000;
        break;
      case 60:
        descreen = 0.80000000;
        break;
      case 61:
        descreen = 0.80000000;
        break;
      case 62:
        descreen = 0.80000000;
        break;
      case 63:
        descreen = 0.80000000;
        break;
      case 64:
        descreen = 0.80000000;
        break;
      case 65:
        descreen = 0.80000000;
        break;
      case 66:
        descreen = 0.80000000;
        break;
      case 67:
        descreen = 0.80000000;
        break;
      case 68:
        descreen = 0.80000000;
        break;
      case 69:
        descreen = 0.80000000;
        break;
      case 70:
        descreen = 0.80000000;
        break;
      case 71:
        descreen = 0.80000000;
        break;
      case 72:
        descreen = 1.25539355;
        break;
      case 73:
        descreen = 1.09526296;
        break;
      case 74:
        descreen = 0.79240070;
        break;
      case 75:
        descreen = 1.27346989;
        break;
      case 76:
        descreen = 1.22445252;
        break;
      case 77:
        descreen = 1.41515774;
        break;
      case 78:
        descreen = 1.17588727;
        break;
      case 79:
        descreen = 1.08983202;
        break;
      case 80:
        descreen = 0.82623981;
        break;
      case 81:
        descreen = 0.56362162;
        break;
      case 82:
        descreen = 0.77451211;
        break;
      case 83:
        descreen = 0.56368506;
        break;
      case 84:
        descreen = 0.80000000;
        break;
      case 85:
        descreen = 0.80000000;
        break;
      case 86:
        descreen = 0.80000000;
        break;
      case 87:
        descreen = 0.80000000;
        break;
      case 88:
        descreen = 0.80000000;
        break;
      case 89:
        descreen = 0.80000000;
        break;
      case 90:
        descreen = 0.80000000;
        break;
      case 91:
        descreen = 0.80000000;
        break;
      case 92:
        descreen = 0.80000000;
        break;
      case 93:
        descreen = 0.80000000;
        break;
      case 94:
        descreen = 0.80000000;
        break;
    }
  }
  else if (solvent == "diethyl ether") {
    switch (atmnr) {
      case 1:
        descreen = 0.85545181;
        break;
      case 2:
        descreen = 1.26078863;
        break;
      case 3:
        descreen = 1.44092061;
        break;
      case 4:
        descreen = 0.46761613;
        break;
      case 5:
        descreen = 0.91673162;
        break;
      case 6:
        descreen = 0.02464722;
        break;
      case 7:
        descreen = 0.34310332;
        break;
      case 8:
        descreen = 0.30000000;
        break;
      case 9:
        descreen = 0.85470383;
        break;
      case 10:
        descreen = 0.30000000;
        break;
      case 11:
        descreen = 0.30000000;
        break;
      case 12:
        descreen = 1.00408544;
        break;
      case 13:
        descreen = 1.02412550;
        break;
      case 14:
        descreen = 0.89056740;
        break;
      case 15:
        descreen = 0.44358004;
        break;
      case 16:
        descreen = 0.85825998;
        break;
      case 17:
        descreen = 0.84280882;
        break;
      case 18:
        descreen = 1.03427028;
        break;
      case 19:
        descreen = 1.27691271;
        break;
      case 20:
        descreen = 1.34442198;
        break;
      case 21:
        descreen = 1.03159986;
        break;
      case 22:
        descreen = 1.07670987;
        break;
      case 23:
        descreen = 1.35114132;
        break;
      case 24:
        descreen = 1.21955226;
        break;
      case 25:
        descreen = 1.27263052;
        break;
      case 26:
        descreen = 1.54570751;
        break;
      case 27:
        descreen = 1.07348117;
        break;
      case 28:
        descreen = 1.04116392;
        break;
      case 29:
        descreen = 0.97593637;
        break;
      case 30:
        descreen = 0.07266388;
        break;
      case 31:
        descreen = 1.13433161;
        break;
      case 32:
        descreen = 0.75819466;
        break;
      case 33:
        descreen = 1.12659303;
        break;
      case 34:
        descreen = 1.08552585;
        break;
      case 35:
        descreen = 0.39180649;
        break;
      case 36:
        descreen = 1.33367354;
        break;
      case 37:
        descreen = 1.23210461;
        break;
      case 38:
        descreen = 1.25112664;
        break;
      case 39:
        descreen = 1.36363262;
        break;
      case 40:
        descreen = 0.87729504;
        break;
      case 41:
        descreen = 1.37232120;
        break;
      case 42:
        descreen = 0.56546146;
        break;
      case 43:
        descreen = 1.43209178;
        break;
      case 44:
        descreen = 1.56196204;
        break;
      case 45:
        descreen = 1.29955995;
        break;
      case 46:
        descreen = 1.14954036;
        break;
      case 47:
        descreen = 1.27260356;
        break;
      case 48:
        descreen = 1.51423549;
        break;
      case 49:
        descreen = 0.83495527;
        break;
      case 50:
        descreen = 0.61288562;
        break;
      case 51:
        descreen = 0.64530264;
        break;
      case 52:
        descreen = 0.24151854;
        break;
      case 53:
        descreen = 0.30000000;
        break;
      case 54:
        descreen = 1.14775060;
        break;
      case 55:
        descreen = 1.14084866;
        break;
      case 56:
        descreen = 0.93117727;
        break;
      case 57:
        descreen = 0.80000000;
        break;
      case 58:
        descreen = 0.80000000;
        break;
      case 59:
        descreen = 0.80000000;
        break;
      case 60:
        descreen = 0.80000000;
        break;
      case 61:
        descreen = 0.80000000;
        break;
      case 62:
        descreen = 0.80000000;
        break;
      case 63:
        descreen = 0.80000000;
        break;
      case 64:
        descreen = 0.80000000;
        break;
      case 65:
        descreen = 0.80000000;
        break;
      case 66:
        descreen = 0.80000000;
        break;
      case 67:
        descreen = 0.80000000;
        break;
      case 68:
        descreen = 0.80000000;
        break;
      case 69:
        descreen = 0.80000000;
        break;
      case 70:
        descreen = 0.80000000;
        break;
      case 71:
        descreen = 0.80000000;
        break;
      case 72:
        descreen = 1.33413719;
        break;
      case 73:
        descreen = 1.17474693;
        break;
      case 74:
        descreen = 0.89109996;
        break;
      case 75:
        descreen = 1.14121649;
        break;
      case 76:
        descreen = 1.45473344;
        break;
      case 77:
        descreen = 1.42630812;
        break;
      case 78:
        descreen = 1.10052032;
        break;
      case 79:
        descreen = 1.08765113;
        break;
      case 80:
        descreen = 0.80000000;
        break;
      case 81:
        descreen = 1.16676059;
        break;
      case 82:
        descreen = 0.99157401;
        break;
      case 83:
        descreen = 1.13438161;
        break;
      case 84:
        descreen = 0.80000000;
        break;
      case 85:
        descreen = 0.80000000;
        break;
      case 86:
        descreen = 0.80000000;
        break;
      case 87:
        descreen = 0.80000000;
        break;
      case 88:
        descreen = 0.80000000;
        break;
      case 89:
        descreen = 0.80000000;
        break;
      case 90:
        descreen = 0.80000000;
        break;
      case 91:
        descreen = 0.80000000;
        break;
      case 92:
        descreen = 0.80000000;
        break;
      case 93:
        descreen = 0.80000000;
        break;
      case 94:
        descreen = 0.80000000;
        break;
    }
  }
  else if (solvent == "ethyl acetate") {
    switch (atmnr) {
      case 1:
        descreen = 0.75327067;
        break;
      case 2:
        descreen = 1.48107773;
        break;
      case 3:
        descreen = 0.30000000;
        break;
      case 4:
        descreen = 0.72056458;
        break;
      case 5:
        descreen = 0.84852207;
        break;
      case 6:
        descreen = 0.48790277;
        break;
      case 7:
        descreen = 0.30000000;
        break;
      case 8:
        descreen = 0.30000000;
        break;
      case 9:
        descreen = 0.68876350;
        break;
      case 10:
        descreen = 0.30000000;
        break;
      case 11:
        descreen = 1.33884934;
        break;
      case 12:
        descreen = 1.06391564;
        break;
      case 13:
        descreen = 1.03031025;
        break;
      case 14:
        descreen = 0.87150586;
        break;
      case 15:
        descreen = 0.52556850;
        break;
      case 16:
        descreen = 0.91708861;
        break;
      case 17:
        descreen = 1.08370970;
        break;
      case 18:
        descreen = 0.30000000;
        break;
      case 19:
        descreen = 0.30000000;
        break;
      case 20:
        descreen = 1.38117372;
        break;
      case 21:
        descreen = 1.03118497;
        break;
      case 22:
        descreen = 1.08835081;
        break;
      case 23:
        descreen = 1.29293858;
        break;
      case 24:
        descreen = 1.45814614;
        break;
      case 25:
        descreen = 1.20013819;
        break;
      case 26:
        descreen = 1.57654044;
        break;
      case 27:
        descreen = 0.66816646;
        break;
      case 28:
        descreen = 0.49846436;
        break;
      case 29:
        descreen = 0.89653473;
        break;
      case 30:
        descreen = 0.30000000;
        break;
      case 31:
        descreen = 1.11103491;
        break;
      case 32:
        descreen = 0.68994820;
        break;
      case 33:
        descreen = 1.03098261;
        break;
      case 34:
        descreen = 1.00425454;
        break;
      case 35:
        descreen = 0.84978396;
        break;
      case 36:
        descreen = 1.10682104;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 0.30000000;
        break;
      case 39:
        descreen = 0.95097273;
        break;
      case 40:
        descreen = 0.30000000;
        break;
      case 41:
        descreen = 1.39571242;
        break;
      case 42:
        descreen = 0.63676822;
        break;
      case 43:
        descreen = 1.43256971;
        break;
      case 44:
        descreen = 1.46774994;
        break;
      case 45:
        descreen = 1.30373035;
        break;
      case 46:
        descreen = 1.08475178;
        break;
      case 47:
        descreen = 1.27165741;
        break;
      case 48:
        descreen = 1.09117135;
        break;
      case 49:
        descreen = 0.84805784;
        break;
      case 50:
        descreen = 0.30000000;
        break;
      case 51:
        descreen = 0.93668467;
        break;
      case 52:
        descreen = 0.70968282;
        break;
      case 53:
        descreen = 0.96861202;
        break;
      case 54:
        descreen = 1.02972030;
        break;
      case 55:
        descreen = 1.00011722;
        break;
      case 56:
        descreen = 0.95132992;
        break;
      case 57:
        descreen = 1.00000000;
        break;
      case 58:
        descreen = 1.00000000;
        break;
      case 59:
        descreen = 1.00000000;
        break;
      case 60:
        descreen = 1.00000000;
        break;
      case 61:
        descreen = 1.00000000;
        break;
      case 62:
        descreen = 1.00000000;
        break;
      case 63:
        descreen = 1.00000000;
        break;
      case 64:
        descreen = 1.00000000;
        break;
      case 65:
        descreen = 1.00000000;
        break;
      case 66:
        descreen = 1.00000000;
        break;
      case 67:
        descreen = 1.00000000;
        break;
      case 68:
        descreen = 1.00000000;
        break;
      case 69:
        descreen = 1.00000000;
        break;
      case 70:
        descreen = 1.00000000;
        break;
      case 71:
        descreen = 1.00000000;
        break;
      case 72:
        descreen = 1.44373140;
        break;
      case 73:
        descreen = 1.20752219;
        break;
      case 74:
        descreen = 0.23928831;
        break;
      case 75:
        descreen = 0.99697589;
        break;
      case 76:
        descreen = 1.37637072;
        break;
      case 77:
        descreen = 1.32145117;
        break;
      case 78:
        descreen = 0.95526657;
        break;
      case 79:
        descreen = 1.03585397;
        break;
      case 80:
        descreen = 0.97599615;
        break;
      case 81:
        descreen = 0.76100333;
        break;
      case 82:
        descreen = 0.98966698;
        break;
      case 83:
        descreen = 0.62279966;
        break;
      case 84:
        descreen = 1.00000000;
        break;
      case 85:
        descreen = 1.00000000;
        break;
      case 86:
        descreen = 1.00000000;
        break;
      case 87:
        descreen = 1.00000000;
        break;
      case 88:
        descreen = 1.00000000;
        break;
      case 89:
        descreen = 1.00000000;
        break;
      case 90:
        descreen = 1.00000000;
        break;
      case 91:
        descreen = 1.00000000;
        break;
      case 92:
        descreen = 1.00000000;
        break;
      case 93:
        descreen = 1.00000000;
        break;
      case 94:
        descreen = 1.00000000;
        break;
    }
  }
  else if (solvent == "furane") {
    switch (atmnr) {
      case 1:
        descreen = 0.58301411;
        break;
      case 2:
        descreen = 1.81350872;
        break;
      case 3:
        descreen = 0.75206698;
        break;
      case 4:
        descreen = 0.61619223;
        break;
      case 5:
        descreen = 0.86799156;
        break;
      case 6:
        descreen = 0.23981977;
        break;
      case 7:
        descreen = 0.30000000;
        break;
      case 8:
        descreen = 0.86919153;
        break;
      case 9:
        descreen = 0.87791119;
        break;
      case 10:
        descreen = 1.26285484;
        break;
      case 11:
        descreen = 1.42990107;
        break;
      case 12:
        descreen = 1.00546613;
        break;
      case 13:
        descreen = 0.92998604;
        break;
      case 14:
        descreen = 0.95474113;
        break;
      case 15:
        descreen = 0.30000000;
        break;
      case 16:
        descreen = 0.66201697;
        break;
      case 17:
        descreen = 1.13833490;
        break;
      case 18:
        descreen = 1.05845403;
        break;
      case 19:
        descreen = 0.30000000;
        break;
      case 20:
        descreen = 0.30000000;
        break;
      case 21:
        descreen = 1.07129882;
        break;
      case 22:
        descreen = 1.07814047;
        break;
      case 23:
        descreen = 1.24208512;
        break;
      case 24:
        descreen = 1.01436759;
        break;
      case 25:
        descreen = 1.20260572;
        break;
      case 26:
        descreen = 1.56541659;
        break;
      case 27:
        descreen = 1.06281194;
        break;
      case 28:
        descreen = 0.76863040;
        break;
      case 29:
        descreen = 1.06281844;
        break;
      case 30:
        descreen = 0.63675905;
        break;
      case 31:
        descreen = 1.15958898;
        break;
      case 32:
        descreen = 0.75122278;
        break;
      case 33:
        descreen = 1.11516041;
        break;
      case 34:
        descreen = 0.81956508;
        break;
      case 35:
        descreen = 0.84910390;
        break;
      case 36:
        descreen = 1.19791276;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 0.30000000;
        break;
      case 39:
        descreen = 0.84196403;
        break;
      case 40:
        descreen = 0.30000000;
        break;
      case 41:
        descreen = 1.23531583;
        break;
      case 42:
        descreen = 0.46887795;
        break;
      case 43:
        descreen = 1.14991646;
        break;
      case 44:
        descreen = 1.47823465;
        break;
      case 45:
        descreen = 1.30663546;
        break;
      case 46:
        descreen = 1.24692545;
        break;
      case 47:
        descreen = 1.03605678;
        break;
      case 48:
        descreen = 0.81294699;
        break;
      case 49:
        descreen = 0.89306108;
        break;
      case 50:
        descreen = 0.30000000;
        break;
      case 51:
        descreen = 1.03252233;
        break;
      case 52:
        descreen = 0.86195584;
        break;
      case 53:
        descreen = 0.73022054;
        break;
      case 54:
        descreen = 0.91170779;
        break;
      case 55:
        descreen = 0.99759987;
        break;
      case 56:
        descreen = 0.87098751;
        break;
      case 57:
        descreen = 1.00000000;
        break;
      case 58:
        descreen = 1.00000000;
        break;
      case 59:
        descreen = 1.00000000;
        break;
      case 60:
        descreen = 1.00000000;
        break;
      case 61:
        descreen = 1.00000000;
        break;
      case 62:
        descreen = 1.00000000;
        break;
      case 63:
        descreen = 1.00000000;
        break;
      case 64:
        descreen = 1.00000000;
        break;
      case 65:
        descreen = 1.00000000;
        break;
      case 66:
        descreen = 1.00000000;
        break;
      case 67:
        descreen = 1.00000000;
        break;
      case 68:
        descreen = 1.00000000;
        break;
      case 69:
        descreen = 1.00000000;
        break;
      case 70:
        descreen = 1.00000000;
        break;
      case 71:
        descreen = 1.00000000;
        break;
      case 72:
        descreen = 1.49296628;
        break;
      case 73:
        descreen = 1.07474522;
        break;
      case 74:
        descreen = 0.51977584;
        break;
      case 75:
        descreen = 1.27241266;
        break;
      case 76:
        descreen = 1.41683942;
        break;
      case 77:
        descreen = 1.19581174;
        break;
      case 78:
        descreen = 1.14400185;
        break;
      case 79:
        descreen = 1.02974405;
        break;
      case 80:
        descreen = 0.81932870;
        break;
      case 81:
        descreen = 0.61068641;
        break;
      case 82:
        descreen = 1.01652783;
        break;
      case 83:
        descreen = 0.73329767;
        break;
      case 84:
        descreen = 1.00000000;
        break;
      case 85:
        descreen = 1.00000000;
        break;
      case 86:
        descreen = 1.00000000;
        break;
      case 87:
        descreen = 1.00000000;
        break;
      case 88:
        descreen = 1.00000000;
        break;
      case 89:
        descreen = 1.00000000;
        break;
      case 90:
        descreen = 1.00000000;
        break;
      case 91:
        descreen = 1.00000000;
        break;
      case 92:
        descreen = 1.00000000;
        break;
      case 93:
        descreen = 1.00000000;
        break;
      case 94:
        descreen = 1.00000000;
        break;
    }
  }
  else if (solvent == "hexadecane") {
    switch (atmnr) {
      case 1:
        descreen = 0.76479783;
        break;
      case 2:
        descreen = 0.30000000;
        break;
      case 3:
        descreen = 1.02085775;
        break;
      case 4:
        descreen = 1.03333496;
        break;
      case 5:
        descreen = 0.79334429;
        break;
      case 6:
        descreen = 0.86233306;
        break;
      case 7:
        descreen = 0.45709544;
        break;
      case 8:
        descreen = 1.07579741;
        break;
      case 9:
        descreen = 0.93167572;
        break;
      case 10:
        descreen = 0.30000000;
        break;
      case 11:
        descreen = 0.71084540;
        break;
      case 12:
        descreen = 0.30000000;
        break;
      case 13:
        descreen = 0.24092595;
        break;
      case 14:
        descreen = 0.94156054;
        break;
      case 15:
        descreen = 0.91957660;
        break;
      case 16:
        descreen = 0.01642782;
        break;
      case 17:
        descreen = 0.80716343;
        break;
      case 18:
        descreen = 1.53628596;
        break;
      case 19:
        descreen = 0.31332033;
        break;
      case 20:
        descreen = 0.30000000;
        break;
      case 21:
        descreen = 0.98257164;
        break;
      case 22:
        descreen = 1.04846793;
        break;
      case 23:
        descreen = 1.21052162;
        break;
      case 24:
        descreen = 0.99524912;
        break;
      case 25:
        descreen = 0.29999966;
        break;
      case 26:
        descreen = 1.64151921;
        break;
      case 27:
        descreen = 0.84683933;
        break;
      case 28:
        descreen = 0.81594847;
        break;
      case 29:
        descreen = 1.07781728;
        break;
      case 30:
        descreen = 0.87088973;
        break;
      case 31:
        descreen = 0.78581903;
        break;
      case 32:
        descreen = 0.80872723;
        break;
      case 33:
        descreen = 1.06348294;
        break;
      case 34:
        descreen = 0.99342894;
        break;
      case 35:
        descreen = 0.40057871;
        break;
      case 36:
        descreen = 1.18569987;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 0.22796659;
        break;
      case 39:
        descreen = 0.74131939;
        break;
      case 40:
        descreen = 0.70936521;
        break;
      case 41:
        descreen = 1.06800972;
        break;
      case 42:
        descreen = 0.29999915;
        break;
      case 43:
        descreen = 1.36418361;
        break;
      case 44:
        descreen = 0.99019169;
        break;
      case 45:
        descreen = 1.18248069;
        break;
      case 46:
        descreen = 1.20939488;
        break;
      case 47:
        descreen = 1.10302503;
        break;
      case 48:
        descreen = 0.90747278;
        break;
      case 49:
        descreen = 0.82296070;
        break;
      case 50:
        descreen = 0.67425979;
        break;
      case 51:
        descreen = 1.06945825;
        break;
      case 52:
        descreen = 0.30000000;
        break;
      case 53:
        descreen = 0.30000000;
        break;
      case 54:
        descreen = 0.45726383;
        break;
      case 55:
        descreen = 0.72521971;
        break;
      case 56:
        descreen = 0.68105377;
        break;
      case 57:
        descreen = 0.80000000;
        break;
      case 58:
        descreen = 0.80000000;
        break;
      case 59:
        descreen = 0.80000000;
        break;
      case 60:
        descreen = 0.80000000;
        break;
      case 61:
        descreen = 0.80000000;
        break;
      case 62:
        descreen = 0.80000000;
        break;
      case 63:
        descreen = 0.80000000;
        break;
      case 64:
        descreen = 0.80000000;
        break;
      case 65:
        descreen = 0.80000000;
        break;
      case 66:
        descreen = 0.80000000;
        break;
      case 67:
        descreen = 0.80000000;
        break;
      case 68:
        descreen = 0.80000000;
        break;
      case 69:
        descreen = 0.80000000;
        break;
      case 70:
        descreen = 0.80000000;
        break;
      case 71:
        descreen = 0.80000000;
        break;
      case 72:
        descreen = 1.01859677;
        break;
      case 73:
        descreen = 0.83518238;
        break;
      case 74:
        descreen = 0.80000000;
        break;
      case 75:
        descreen = 0.78019293;
        break;
      case 76:
        descreen = 1.15244777;
        break;
      case 77:
        descreen = 1.09219744;
        break;
      case 78:
        descreen = 1.08008536;
        break;
      case 79:
        descreen = 1.21039755;
        break;
      case 80:
        descreen = 1.26722417;
        break;
      case 81:
        descreen = 0.80209242;
        break;
      case 82:
        descreen = 0.99702187;
        break;
      case 83:
        descreen = 0.89177874;
        break;
      case 84:
        descreen = 0.80000000;
        break;
      case 85:
        descreen = 0.80000000;
        break;
      case 86:
        descreen = 0.80000000;
        break;
      case 87:
        descreen = 0.80000000;
        break;
      case 88:
        descreen = 0.80000000;
        break;
      case 89:
        descreen = 0.80000000;
        break;
      case 90:
        descreen = 0.80000000;
        break;
      case 91:
        descreen = 0.80000000;
        break;
      case 92:
        descreen = 0.80000000;
        break;
      case 93:
        descreen = 0.80000000;
        break;
      case 94:
        descreen = 0.80000000;
        break;
    }
  }
  else if (solvent == "hexane") {
    switch (atmnr) {
      case 1:
        descreen = 0.69492201;
        break;
      case 2:
        descreen = 1.46480233;
        break;
      case 3:
        descreen = 0.88266262;
        break;
      case 4:
        descreen = 0.00144423;
        break;
      case 5:
        descreen = 0.90227475;
        break;
      case 6:
        descreen = 0.33834686;
        break;
      case 7:
        descreen = 0.38580647;
        break;
      case 8:
        descreen = 1.07977497;
        break;
      case 9:
        descreen = 1.01946398;
        break;
      case 10:
        descreen = 0.30000000;
        break;
      case 11:
        descreen = 0.57684540;
        break;
      case 12:
        descreen = 0.90781813;
        break;
      case 13:
        descreen = 0.67929421;
        break;
      case 14:
        descreen = 0.97897332;
        break;
      case 15:
        descreen = 0.30000000;
        break;
      case 16:
        descreen = 0.50372430;
        break;
      case 17:
        descreen = 1.03507144;
        break;
      case 18:
        descreen = 1.57196365;
        break;
      case 19:
        descreen = 0.90206611;
        break;
      case 20:
        descreen = 0.30000000;
        break;
      case 21:
        descreen = 1.05593017;
        break;
      case 22:
        descreen = 1.06821630;
        break;
      case 23:
        descreen = 1.30187340;
        break;
      case 24:
        descreen = 0.83631747;
        break;
      case 25:
        descreen = 1.17331634;
        break;
      case 26:
        descreen = 1.69903371;
        break;
      case 27:
        descreen = 1.49816563;
        break;
      case 28:
        descreen = 1.24499888;
        break;
      case 29:
        descreen = 1.12692444;
        break;
      case 30:
        descreen = 1.16527439;
        break;
      case 31:
        descreen = 0.96261700;
        break;
      case 32:
        descreen = 0.82563991;
        break;
      case 33:
        descreen = 1.10554100;
        break;
      case 34:
        descreen = 0.95042189;
        break;
      case 35:
        descreen = 0.40057871;
        break;
      case 36:
        descreen = 1.40007868;
        break;
      case 37:
        descreen = 0.44393709;
        break;
      case 38:
        descreen = 0.30000000;
        break;
      case 39:
        descreen = 1.23057528;
        break;
      case 40:
        descreen = 0.30000000;
        break;
      case 41:
        descreen = 1.10084797;
        break;
      case 42:
        descreen = 0.30000000;
        break;
      case 43:
        descreen = 0.73933959;
        break;
      case 44:
        descreen = 1.57771337;
        break;
      case 45:
        descreen = 1.40159338;
        break;
      case 46:
        descreen = 1.66948169;
        break;
      case 47:
        descreen = 1.29849830;
        break;
      case 48:
        descreen = 0.96226919;
        break;
      case 49:
        descreen = 0.82060847;
        break;
      case 50:
        descreen = 0.71404762;
        break;
      case 51:
        descreen = 0.62364944;
        break;
      case 52:
        descreen = 0.77643816;
        break;
      case 53:
        descreen = 0.30000000;
        break;
      case 54:
        descreen = 1.06922363;
        break;
      case 55:
        descreen = 0.78943942;
        break;
      case 56:
        descreen = 0.77872529;
        break;
      case 57:
        descreen = 0.80000000;
        break;
      case 58:
        descreen = 0.80000000;
        break;
      case 59:
        descreen = 0.80000000;
        break;
      case 60:
        descreen = 0.80000000;
        break;
      case 61:
        descreen = 0.80000000;
        break;
      case 62:
        descreen = 0.80000000;
        break;
      case 63:
        descreen = 0.80000000;
        break;
      case 64:
        descreen = 0.80000000;
        break;
      case 65:
        descreen = 0.80000000;
        break;
      case 66:
        descreen = 0.80000000;
        break;
      case 67:
        descreen = 0.80000000;
        break;
      case 68:
        descreen = 0.80000000;
        break;
      case 69:
        descreen = 0.80000000;
        break;
      case 70:
        descreen = 0.80000000;
        break;
      case 71:
        descreen = 0.80000000;
        break;
      case 72:
        descreen = 1.52041452;
        break;
      case 73:
        descreen = 1.03401131;
        break;
      case 74:
        descreen = 0.60970029;
        break;
      case 75:
        descreen = 1.29146134;
        break;
      case 76:
        descreen = 1.10022245;
        break;
      case 77:
        descreen = 1.07142522;
        break;
      case 78:
        descreen = 1.24845561;
        break;
      case 79:
        descreen = 1.16615750;
        break;
      case 80:
        descreen = 1.31600834;
        break;
      case 81:
        descreen = 0.81779417;
        break;
      case 82:
        descreen = 1.03567786;
        break;
      case 83:
        descreen = 1.13705067;
        break;
      case 84:
        descreen = 0.80000000;
        break;
      case 85:
        descreen = 0.80000000;
        break;
      case 86:
        descreen = 0.80000000;
        break;
      case 87:
        descreen = 0.80000000;
        break;
      case 88:
        descreen = 0.80000000;
        break;
      case 89:
        descreen = 0.80000000;
        break;
      case 90:
        descreen = 0.80000000;
        break;
      case 91:
        descreen = 0.80000000;
        break;
      case 92:
        descreen = 0.80000000;
        break;
      case 93:
        descreen = 0.80000000;
        break;
      case 94:
        descreen = 0.80000000;
        break;
    }
  }
  else if (solvent == "methanol") {
    switch (atmnr) {
      case 1:
        descreen = 0.30000000;
        break;
      case 2:
        descreen = 1.66734117;
        break;
      case 3:
        descreen = 0.30000000;
        break;
      case 4:
        descreen = 0.30000000;
        break;
      case 5:
        descreen = 0.92476227;
        break;
      case 6:
        descreen = 0.57913331;
        break;
      case 7:
        descreen = 0.70841493;
        break;
      case 8:
        descreen = 0.30000000;
        break;
      case 9:
        descreen = 0.62979186;
        break;
      case 10:
        descreen = 0.30000000;
        break;
      case 11:
        descreen = 0.30000000;
        break;
      case 12:
        descreen = 0.30000000;
        break;
      case 13:
        descreen = 0.68727478;
        break;
      case 14:
        descreen = 0.89791803;
        break;
      case 15:
        descreen = 0.30000000;
        break;
      case 16:
        descreen = 0.30000000;
        break;
      case 17:
        descreen = 0.72026418;
        break;
      case 18:
        descreen = 0.11959209;
        break;
      case 19:
        descreen = 0.30000000;
        break;
      case 20:
        descreen = 1.05739468;
        break;
      case 21:
        descreen = 1.08100951;
        break;
      case 22:
        descreen = 1.05255421;
        break;
      case 23:
        descreen = 1.19979404;
        break;
      case 24:
        descreen = 1.10394208;
        break;
      case 25:
        descreen = 1.22055934;
        break;
      case 26:
        descreen = 1.33316412;
        break;
      case 27:
        descreen = 1.06388051;
        break;
      case 28:
        descreen = 1.13562868;
        break;
      case 29:
        descreen = 1.00404400;
        break;
      case 30:
        descreen = 0.42063679;
        break;
      case 31:
        descreen = 1.13731433;
        break;
      case 32:
        descreen = 0.71952599;
        break;
      case 33:
        descreen = 1.12219028;
        break;
      case 34:
        descreen = 1.02909322;
        break;
      case 35:
        descreen = 0.56687555;
        break;
      case 36:
        descreen = 1.12191212;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 1.00073021;
        break;
      case 39:
        descreen = 1.40232008;
        break;
      case 40:
        descreen = 0.62162675;
        break;
      case 41:
        descreen = 1.25434250;
        break;
      case 42:
        descreen = 0.30000000;
        break;
      case 43:
        descreen = 1.32082304;
        break;
      case 44:
        descreen = 1.44966217;
        break;
      case 45:
        descreen = 1.27758140;
        break;
      case 46:
        descreen = 1.21596161;
        break;
      case 47:
        descreen = 1.28170823;
        break;
      case 48:
        descreen = 0.60166517;
        break;
      case 49:
        descreen = 0.83425879;
        break;
      case 50:
        descreen = 0.30000000;
        break;
      case 51:
        descreen = 0.49853397;
        break;
      case 52:
        descreen = 0.70709437;
        break;
      case 53:
        descreen = 1.45485315;
        break;
      case 54:
        descreen = 1.01437992;
        break;
      case 55:
        descreen = 0.30000000;
        break;
      case 56:
        descreen = 0.82700576;
        break;
      case 57:
        descreen = 0.80000000;
        break;
      case 58:
        descreen = 0.80000000;
        break;
      case 59:
        descreen = 0.80000000;
        break;
      case 60:
        descreen = 0.80000000;
        break;
      case 61:
        descreen = 0.80000000;
        break;
      case 62:
        descreen = 0.80000000;
        break;
      case 63:
        descreen = 0.80000000;
        break;
      case 64:
        descreen = 0.80000000;
        break;
      case 65:
        descreen = 0.80000000;
        break;
      case 66:
        descreen = 0.80000000;
        break;
      case 67:
        descreen = 0.80000000;
        break;
      case 68:
        descreen = 0.80000000;
        break;
      case 69:
        descreen = 0.80000000;
        break;
      case 70:
        descreen = 0.80000000;
        break;
      case 71:
        descreen = 0.80000000;
        break;
      case 72:
        descreen = 1.29507907;
        break;
      case 73:
        descreen = 1.13915241;
        break;
      case 74:
        descreen = 0.74105172;
        break;
      case 75:
        descreen = 1.31265698;
        break;
      case 76:
        descreen = 1.29112147;
        break;
      case 77:
        descreen = 1.43939838;
        break;
      case 78:
        descreen = 1.16855363;
        break;
      case 79:
        descreen = 1.07319850;
        break;
      case 80:
        descreen = 0.37117984;
        break;
      case 81:
        descreen = 0.50699717;
        break;
      case 82:
        descreen = 0.74663209;
        break;
      case 83:
        descreen = 0.38505688;
        break;
      case 84:
        descreen = 0.80000000;
        break;
      case 85:
        descreen = 0.80000000;
        break;
      case 86:
        descreen = 0.80000000;
        break;
      case 87:
        descreen = 0.80000000;
        break;
      case 88:
        descreen = 0.80000000;
        break;
      case 89:
        descreen = 0.80000000;
        break;
      case 90:
        descreen = 0.80000000;
        break;
      case 91:
        descreen = 0.80000000;
        break;
      case 92:
        descreen = 0.80000000;
        break;
      case 93:
        descreen = 0.80000000;
        break;
      case 94:
        descreen = 0.80000000;
        break;
    }
  }
  else if (solvent == "nitromethane") {
    switch (atmnr) {
      case 1:
        descreen = 0.54545897;
        break;
      case 2:
        descreen = 1.81867594;
        break;
      case 3:
        descreen = 0.30000000;
        break;
      case 4:
        descreen = 0.80637882;
        break;
      case 5:
        descreen = 0.90063169;
        break;
      case 6:
        descreen = 0.84417780;
        break;
      case 7:
        descreen = 0.55326728;
        break;
      case 8:
        descreen = 0.30000000;
        break;
      case 9:
        descreen = 0.62755862;
        break;
      case 10:
        descreen = 0.30090875;
        break;
      case 11:
        descreen = 0.30000000;
        break;
      case 12:
        descreen = 1.08266495;
        break;
      case 13:
        descreen = 1.01141441;
        break;
      case 14:
        descreen = 0.87362989;
        break;
      case 15:
        descreen = 0.78614730;
        break;
      case 16:
        descreen = 0.94616088;
        break;
      case 17:
        descreen = 0.76931269;
        break;
      case 18:
        descreen = 0.30000000;
        break;
      case 19:
        descreen = 0.30000000;
        break;
      case 20:
        descreen = 0.30000000;
        break;
      case 21:
        descreen = 1.00975751;
        break;
      case 22:
        descreen = 1.07599883;
        break;
      case 23:
        descreen = 1.09618488;
        break;
      case 24:
        descreen = 1.39578535;
        break;
      case 25:
        descreen = 1.25769739;
        break;
      case 26:
        descreen = 1.47835942;
        break;
      case 27:
        descreen = 0.14034901;
        break;
      case 28:
        descreen = 0.05435367;
        break;
      case 29:
        descreen = 0.30000000;
        break;
      case 30:
        descreen = 0.61039860;
        break;
      case 31:
        descreen = 0.98374894;
        break;
      case 32:
        descreen = 0.70465720;
        break;
      case 33:
        descreen = 0.91158866;
        break;
      case 34:
        descreen = 0.90365475;
        break;
      case 35:
        descreen = 0.86995096;
        break;
      case 36:
        descreen = 0.85040817;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 0.30000000;
        break;
      case 39:
        descreen = 1.05596541;
        break;
      case 40:
        descreen = 0.30000000;
        break;
      case 41:
        descreen = 1.14094110;
        break;
      case 42:
        descreen = 1.30074682;
        break;
      case 43:
        descreen = 1.30773420;
        break;
      case 44:
        descreen = 1.28373226;
        break;
      case 45:
        descreen = 1.32597636;
        break;
      case 46:
        descreen = 0.88357203;
        break;
      case 47:
        descreen = 0.77252305;
        break;
      case 48:
        descreen = 1.08651889;
        break;
      case 49:
        descreen = 0.85311149;
        break;
      case 50:
        descreen = 0.30000000;
        break;
      case 51:
        descreen = 0.93716128;
        break;
      case 52:
        descreen = 0.95357198;
        break;
      case 53:
        descreen = 1.01760593;
        break;
      case 54:
        descreen = 0.91464786;
        break;
      case 55:
        descreen = 0.99892740;
        break;
      case 56:
        descreen = 0.95075944;
        break;
      case 57:
        descreen = 1.00000000;
        break;
      case 58:
        descreen = 1.00000000;
        break;
      case 59:
        descreen = 1.00000000;
        break;
      case 60:
        descreen = 1.00000000;
        break;
      case 61:
        descreen = 1.00000000;
        break;
      case 62:
        descreen = 1.00000000;
        break;
      case 63:
        descreen = 1.00000000;
        break;
      case 64:
        descreen = 1.00000000;
        break;
      case 65:
        descreen = 1.00000000;
        break;
      case 66:
        descreen = 1.00000000;
        break;
      case 67:
        descreen = 1.00000000;
        break;
      case 68:
        descreen = 1.00000000;
        break;
      case 69:
        descreen = 1.00000000;
        break;
      case 70:
        descreen = 1.00000000;
        break;
      case 71:
        descreen = 1.00000000;
        break;
      case 72:
        descreen = 1.23243197;
        break;
      case 73:
        descreen = 1.00100748;
        break;
      case 74:
        descreen = 0.30014100;
        break;
      case 75:
        descreen = 0.97529690;
        break;
      case 76:
        descreen = 1.30486029;
        break;
      case 77:
        descreen = 1.22948266;
        break;
      case 78:
        descreen = 0.97712571;
        break;
      case 79:
        descreen = 0.85423836;
        break;
      case 80:
        descreen = 0.30000000;
        break;
      case 81:
        descreen = 0.70334223;
        break;
      case 82:
        descreen = 0.99403972;
        break;
      case 83:
        descreen = 0.59349783;
        break;
      case 84:
        descreen = 1.00000000;
        break;
      case 85:
        descreen = 1.00000000;
        break;
      case 86:
        descreen = 1.00000000;
        break;
      case 87:
        descreen = 1.00000000;
        break;
      case 88:
        descreen = 1.00000000;
        break;
      case 89:
        descreen = 1.00000000;
        break;
      case 90:
        descreen = 1.00000000;
        break;
      case 91:
        descreen = 1.00000000;
        break;
      case 92:
        descreen = 1.00000000;
        break;
      case 93:
        descreen = 1.00000000;
        break;
      case 94:
        descreen = 1.00000000;
        break;
    }
  }
  else if (solvent == "octanol") {
    switch (atmnr) {
      case 1:
        descreen = 0.60595605;
        break;
      case 2:
        descreen = 1.50092262;
        break;
      case 3:
        descreen = 0.30000000;
        break;
      case 4:
        descreen = 0.88335269;
        break;
      case 5:
        descreen = 0.95929612;
        break;
      case 6:
        descreen = 0.71841955;
        break;
      case 7:
        descreen = 0.30000000;
        break;
      case 8:
        descreen = 0.97405294;
        break;
      case 9:
        descreen = 0.68654417;
        break;
      case 10:
        descreen = 1.67610880;
        break;
      case 11:
        descreen = 1.26953387;
        break;
      case 12:
        descreen = 1.02342209;
        break;
      case 13:
        descreen = 0.92816603;
        break;
      case 14:
        descreen = 0.85331623;
        break;
      case 15:
        descreen = 1.08775462;
        break;
      case 16:
        descreen = 1.17022711;
        break;
      case 17:
        descreen = 0.03352880;
        break;
      case 18:
        descreen = 1.11273480;
        break;
      case 19:
        descreen = 0.30000000;
        break;
      case 20:
        descreen = 1.26136129;
        break;
      case 21:
        descreen = 0.99810251;
        break;
      case 22:
        descreen = 1.05324997;
        break;
      case 23:
        descreen = 1.19891417;
        break;
      case 24:
        descreen = 1.15274461;
        break;
      case 25:
        descreen = 1.14151588;
        break;
      case 26:
        descreen = 1.12131571;
        break;
      case 27:
        descreen = 1.45576586;
        break;
      case 28:
        descreen = 1.26892081;
        break;
      case 29:
        descreen = 0.99905738;
        break;
      case 30:
        descreen = 1.19671061;
        break;
      case 31:
        descreen = 1.12420265;
        break;
      case 32:
        descreen = 0.68616739;
        break;
      case 33:
        descreen = 1.09094258;
        break;
      case 34:
        descreen = 1.25179022;
        break;
      case 35:
        descreen = 0.89916614;
        break;
      case 36:
        descreen = 1.08812054;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 0.68020607;
        break;
      case 39:
        descreen = 1.54032209;
        break;
      case 40:
        descreen = 0.98448430;
        break;
      case 41:
        descreen = 1.18532041;
        break;
      case 42:
        descreen = 1.56616088;
        break;
      case 43:
        descreen = 1.03875963;
        break;
      case 44:
        descreen = 1.29493511;
        break;
      case 45:
        descreen = 1.19846948;
        break;
      case 46:
        descreen = 1.16870153;
        break;
      case 47:
        descreen = 1.26440981;
        break;
      case 48:
        descreen = 1.18601017;
        break;
      case 49:
        descreen = 0.91262881;
        break;
      case 50:
        descreen = 0.84986655;
        break;
      case 51:
        descreen = 1.07665288;
        break;
      case 52:
        descreen = 0.30000000;
        break;
      case 53:
        descreen = 0.74665903;
        break;
      case 54:
        descreen = 0.92072223;
        break;
      case 55:
        descreen = 0.30000000;
        break;
      case 56:
        descreen = 1.02176174;
        break;
      case 57:
        descreen = 1.00000000;
        break;
      case 58:
        descreen = 1.00000000;
        break;
      case 59:
        descreen = 1.00000000;
        break;
      case 60:
        descreen = 1.00000000;
        break;
      case 61:
        descreen = 1.00000000;
        break;
      case 62:
        descreen = 1.00000000;
        break;
      case 63:
        descreen = 1.00000000;
        break;
      case 64:
        descreen = 1.00000000;
        break;
      case 65:
        descreen = 1.00000000;
        break;
      case 66:
        descreen = 1.00000000;
        break;
      case 67:
        descreen = 1.00000000;
        break;
      case 68:
        descreen = 1.00000000;
        break;
      case 69:
        descreen = 1.00000000;
        break;
      case 70:
        descreen = 1.00000000;
        break;
      case 71:
        descreen = 1.00000000;
        break;
      case 72:
        descreen = 1.25958770;
        break;
      case 73:
        descreen = 1.00757395;
        break;
      case 74:
        descreen = 1.36792732;
        break;
      case 75:
        descreen = 1.19710104;
        break;
      case 76:
        descreen = 1.09977357;
        break;
      case 77:
        descreen = 1.37108016;
        break;
      case 78:
        descreen = 1.18284767;
        break;
      case 79:
        descreen = 1.09762567;
        break;
      case 80:
        descreen = 1.32310198;
        break;
      case 81:
        descreen = 0.81664106;
        break;
      case 82:
        descreen = 0.87125841;
        break;
      case 83:
        descreen = 1.10653215;
        break;
      case 84:
        descreen = 1.00000000;
        break;
      case 85:
        descreen = 1.00000000;
        break;
      case 86:
        descreen = 1.00000000;
        break;
      case 87:
        descreen = 1.00000000;
        break;
      case 88:
        descreen = 1.00000000;
        break;
      case 89:
        descreen = 1.00000000;
        break;
      case 90:
        descreen = 1.00000000;
        break;
      case 91:
        descreen = 1.00000000;
        break;
      case 92:
        descreen = 1.00000000;
        break;
      case 93:
        descreen = 1.00000000;
        break;
      case 94:
        descreen = 1.00000000;
        break;
    }
  }
  else if (solvent == "phenol") {
    switch (atmnr) {
      case 1:
        descreen = 0.59187702;
        break;
      case 2:
        descreen = 0.30000000;
        break;
      case 3:
        descreen = 0.76775826;
        break;
      case 4:
        descreen = 0.79943133;
        break;
      case 5:
        descreen = 0.76401539;
        break;
      case 6:
        descreen = 0.29587189;
        break;
      case 7:
        descreen = 0.72550165;
        break;
      case 8:
        descreen = 1.05554190;
        break;
      case 9:
        descreen = 0.89719313;
        break;
      case 10:
        descreen = 0.30000000;
        break;
      case 11:
        descreen = 0.30000000;
        break;
      case 12:
        descreen = 1.09268922;
        break;
      case 13:
        descreen = 0.98648472;
        break;
      case 14:
        descreen = 0.70542229;
        break;
      case 15:
        descreen = 0.37763231;
        break;
      case 16:
        descreen = 1.06429459;
        break;
      case 17:
        descreen = 1.51943692;
        break;
      case 18:
        descreen = 1.15277819;
        break;
      case 19:
        descreen = 0.30000000;
        break;
      case 20:
        descreen = 1.08425669;
        break;
      case 21:
        descreen = 0.73486167;
        break;
      case 22:
        descreen = 1.04153584;
        break;
      case 23:
        descreen = 1.12843007;
        break;
      case 24:
        descreen = 1.49573736;
        break;
      case 25:
        descreen = 1.14841065;
        break;
      case 26:
        descreen = 1.28556943;
        break;
      case 27:
        descreen = 1.00000000;
        break;
      case 28:
        descreen = 0.85503425;
        break;
      case 29:
        descreen = 0.81217971;
        break;
      case 30:
        descreen = 1.20751137;
        break;
      case 31:
        descreen = 1.03918370;
        break;
      case 32:
        descreen = 0.76807984;
        break;
      case 33:
        descreen = 1.05362080;
        break;
      case 34:
        descreen = 0.95541971;
        break;
      case 35:
        descreen = 0.07651436;
        break;
      case 36:
        descreen = 1.38094206;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 0.30000000;
        break;
      case 39:
        descreen = 0.30000000;
        break;
      case 40:
        descreen = 0.76974940;
        break;
      case 41:
        descreen = 1.20584949;
        break;
      case 42:
        descreen = 0.99868109;
        break;
      case 43:
        descreen = 1.22358411;
        break;
      case 44:
        descreen = 1.38842844;
        break;
      case 45:
        descreen = 1.15301205;
        break;
      case 46:
        descreen = 1.21724323;
        break;
      case 47:
        descreen = 1.04033395;
        break;
      case 48:
        descreen = 1.36333247;
        break;
      case 49:
        descreen = 0.87900885;
        break;
      case 50:
        descreen = 0.30000000;
        break;
      case 51:
        descreen = 1.08369143;
        break;
      case 52:
        descreen = 0.85168772;
        break;
      case 53:
        descreen = 0.20788426;
        break;
      case 54:
        descreen = 1.26408986;
        break;
      case 55:
        descreen = 0.30000000;
        break;
      case 56:
        descreen = 1.01146385;
        break;
      case 57:
        descreen = 1.00000000;
        break;
      case 58:
        descreen = 1.00000000;
        break;
      case 59:
        descreen = 1.00000000;
        break;
      case 60:
        descreen = 1.00000000;
        break;
      case 61:
        descreen = 1.00000000;
        break;
      case 62:
        descreen = 1.00000000;
        break;
      case 63:
        descreen = 1.00000000;
        break;
      case 64:
        descreen = 1.00000000;
        break;
      case 65:
        descreen = 1.00000000;
        break;
      case 66:
        descreen = 1.00000000;
        break;
      case 67:
        descreen = 1.00000000;
        break;
      case 68:
        descreen = 1.00000000;
        break;
      case 69:
        descreen = 1.00000000;
        break;
      case 70:
        descreen = 1.00000000;
        break;
      case 71:
        descreen = 1.00000000;
        break;
      case 72:
        descreen = 1.24585931;
        break;
      case 73:
        descreen = 1.05975023;
        break;
      case 74:
        descreen = 0.61491811;
        break;
      case 75:
        descreen = 0.76996490;
        break;
      case 76:
        descreen = 1.16306253;
        break;
      case 77:
        descreen = 1.38575810;
        break;
      case 78:
        descreen = 1.22563783;
        break;
      case 79:
        descreen = 1.08568538;
        break;
      case 80:
        descreen = 0.12982227;
        break;
      case 81:
        descreen = 0.40091242;
        break;
      case 82:
        descreen = 1.00917281;
        break;
      case 83:
        descreen = 0.79104123;
        break;
      case 84:
        descreen = 1.00000000;
        break;
      case 85:
        descreen = 1.00000000;
        break;
      case 86:
        descreen = 1.00000000;
        break;
      case 87:
        descreen = 1.00000000;
        break;
      case 88:
        descreen = 1.00000000;
        break;
      case 89:
        descreen = 1.00000000;
        break;
      case 90:
        descreen = 1.00000000;
        break;
      case 91:
        descreen = 1.00000000;
        break;
      case 92:
        descreen = 1.00000000;
        break;
      case 93:
        descreen = 1.00000000;
        break;
      case 94:
        descreen = 1.00000000;
        break;
    }
  }
  else if (solvent == "thf") {
    switch (atmnr) {
      case 1:
        descreen = 0.87250470;
        break;
      case 2:
        descreen = 1.10416477;
        break;
      case 3:
        descreen = 0.31945941;
        break;
      case 4:
        descreen = 0.65872837;
        break;
      case 5:
        descreen = 0.79162434;
        break;
      case 6:
        descreen = 0.30000000;
        break;
      case 7:
        descreen = 0.30000000;
        break;
      case 8:
        descreen = 0.30000000;
        break;
      case 9:
        descreen = 0.79792828;
        break;
      case 10:
        descreen = 0.30000000;
        break;
      case 11:
        descreen = 0.30000000;
        break;
      case 12:
        descreen = 1.02374372;
        break;
      case 13:
        descreen = 0.70685905;
        break;
      case 14:
        descreen = 0.86421406;
        break;
      case 15:
        descreen = 0.58100025;
        break;
      case 16:
        descreen = 0.95186765;
        break;
      case 17:
        descreen = 0.76927451;
        break;
      case 18:
        descreen = 0.84904260;
        break;
      case 19:
        descreen = 0.30000000;
        break;
      case 20:
        descreen = 0.30000000;
        break;
      case 21:
        descreen = 1.01186506;
        break;
      case 22:
        descreen = 1.07013676;
        break;
      case 23:
        descreen = 1.31744070;
        break;
      case 24:
        descreen = 1.35550797;
        break;
      case 25:
        descreen = 1.25692141;
        break;
      case 26:
        descreen = 1.53278820;
        break;
      case 27:
        descreen = 0.86935765;
        break;
      case 28:
        descreen = 1.17794267;
        break;
      case 29:
        descreen = 0.92855498;
        break;
      case 30:
        descreen = 0.30000000;
        break;
      case 31:
        descreen = 1.11570278;
        break;
      case 32:
        descreen = 0.72884747;
        break;
      case 33:
        descreen = 1.11051328;
        break;
      case 34:
        descreen = 1.05780570;
        break;
      case 35:
        descreen = 0.61925510;
        break;
      case 36:
        descreen = 1.31374189;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 0.30000000;
        break;
      case 39:
        descreen = 1.37177771;
        break;
      case 40:
        descreen = 0.98992867;
        break;
      case 41:
        descreen = 1.35774541;
        break;
      case 42:
        descreen = 0.77673697;
        break;
      case 43:
        descreen = 1.32720260;
        break;
      case 44:
        descreen = 1.52379000;
        break;
      case 45:
        descreen = 1.26437762;
        break;
      case 46:
        descreen = 1.01132798;
        break;
      case 47:
        descreen = 1.25960274;
        break;
      case 48:
        descreen = 0.30000000;
        break;
      case 49:
        descreen = 0.82622509;
        break;
      case 50:
        descreen = 0.67512388;
        break;
      case 51:
        descreen = 1.10710921;
        break;
      case 52:
        descreen = 0.30000000;
        break;
      case 53:
        descreen = 0.57462767;
        break;
      case 54:
        descreen = 1.18921728;
        break;
      case 55:
        descreen = 0.66080140;
        break;
      case 56:
        descreen = 0.95565365;
        break;
      case 57:
        descreen = 0.80000000;
        break;
      case 58:
        descreen = 0.80000000;
        break;
      case 59:
        descreen = 0.80000000;
        break;
      case 60:
        descreen = 0.80000000;
        break;
      case 61:
        descreen = 0.80000000;
        break;
      case 62:
        descreen = 0.80000000;
        break;
      case 63:
        descreen = 0.80000000;
        break;
      case 64:
        descreen = 0.80000000;
        break;
      case 65:
        descreen = 0.80000000;
        break;
      case 66:
        descreen = 0.80000000;
        break;
      case 67:
        descreen = 0.80000000;
        break;
      case 68:
        descreen = 0.80000000;
        break;
      case 69:
        descreen = 0.80000000;
        break;
      case 70:
        descreen = 0.80000000;
        break;
      case 71:
        descreen = 0.80000000;
        break;
      case 72:
        descreen = 1.38799438;
        break;
      case 73:
        descreen = 1.17053682;
        break;
      case 74:
        descreen = 0.92249293;
        break;
      case 75:
        descreen = 0.98656359;
        break;
      case 76:
        descreen = 1.40300471;
        break;
      case 77:
        descreen = 1.35177042;
        break;
      case 78:
        descreen = 1.02852347;
        break;
      case 79:
        descreen = 1.07324743;
        break;
      case 80:
        descreen = 0.30000000;
        break;
      case 81:
        descreen = 0.83275069;
        break;
      case 82:
        descreen = 0.99004346;
        break;
      case 83:
        descreen = 0.51381674;
        break;
      case 84:
        descreen = 0.80000000;
        break;
      case 85:
        descreen = 0.80000000;
        break;
      case 86:
        descreen = 0.80000000;
        break;
      case 87:
        descreen = 0.80000000;
        break;
      case 88:
        descreen = 0.80000000;
        break;
      case 89:
        descreen = 0.80000000;
        break;
      case 90:
        descreen = 0.80000000;
        break;
      case 91:
        descreen = 0.80000000;
        break;
      case 92:
        descreen = 0.80000000;
        break;
      case 93:
        descreen = 0.80000000;
        break;
      case 94:
        descreen = 0.80000000;
        break;
    }
  }
  else if (solvent == "toluene") {
    switch (atmnr) {
      case 1:
        descreen = 0.52616268;
        break;
      case 2:
        descreen = 1.86915069;
        break;
      case 3:
        descreen = 0.30000000;
        break;
      case 4:
        descreen = 0.78775869;
        break;
      case 5:
        descreen = 0.89614831;
        break;
      case 6:
        descreen = 0.91253854;
        break;
      case 7:
        descreen = 0.45660103;
        break;
      case 8:
        descreen = 0.96435408;
        break;
      case 9:
        descreen = 0.59657922;
        break;
      case 10:
        descreen = 1.76204982;
        break;
      case 11:
        descreen = 0.08541653;
        break;
      case 12:
        descreen = 1.03228913;
        break;
      case 13:
        descreen = 0.73804453;
        break;
      case 14:
        descreen = 0.91311889;
        break;
      case 15:
        descreen = 0.80348462;
        break;
      case 16:
        descreen = 0.89182317;
        break;
      case 17:
        descreen = 0.84925544;
        break;
      case 18:
        descreen = 0.67075474;
        break;
      case 19:
        descreen = 0.03719049;
        break;
      case 20:
        descreen = 0.30000000;
        break;
      case 21:
        descreen = 0.93960870;
        break;
      case 22:
        descreen = 1.06511743;
        break;
      case 23:
        descreen = 1.05875564;
        break;
      case 24:
        descreen = 1.09639660;
        break;
      case 25:
        descreen = 0.92755875;
        break;
      case 26:
        descreen = 1.26905678;
        break;
      case 27:
        descreen = 0.57409504;
        break;
      case 28:
        descreen = 0.29086316;
        break;
      case 29:
        descreen = 0.00684059;
        break;
      case 30:
        descreen = 1.19792739;
        break;
      case 31:
        descreen = 0.95251007;
        break;
      case 32:
        descreen = 0.74733021;
        break;
      case 33:
        descreen = 1.03881529;
        break;
      case 34:
        descreen = 1.20383771;
        break;
      case 35:
        descreen = 0.97501677;
        break;
      case 36:
        descreen = 0.95791995;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 0.30000000;
        break;
      case 39:
        descreen = 0.94902740;
        break;
      case 40:
        descreen = 1.02136500;
        break;
      case 41:
        descreen = 1.02262099;
        break;
      case 42:
        descreen = 0.94022824;
        break;
      case 43:
        descreen = 1.32204147;
        break;
      case 44:
        descreen = 1.21056822;
        break;
      case 45:
        descreen = 1.23200391;
        break;
      case 46:
        descreen = 1.17583968;
        break;
      case 47:
        descreen = 0.97207478;
        break;
      case 48:
        descreen = 1.53729651;
        break;
      case 49:
        descreen = 0.87757472;
        break;
      case 50:
        descreen = 0.30000000;
        break;
      case 51:
        descreen = 1.05664166;
        break;
      case 52:
        descreen = 0.99310287;
        break;
      case 53:
        descreen = 0.87799759;
        break;
      case 54:
        descreen = 0.94698088;
        break;
      case 55:
        descreen = 0.30000000;
        break;
      case 56:
        descreen = 0.89414162;
        break;
      case 57:
        descreen = 0.80000000;
        break;
      case 58:
        descreen = 0.80000000;
        break;
      case 59:
        descreen = 0.80000000;
        break;
      case 60:
        descreen = 0.80000000;
        break;
      case 61:
        descreen = 0.80000000;
        break;
      case 62:
        descreen = 0.80000000;
        break;
      case 63:
        descreen = 0.80000000;
        break;
      case 64:
        descreen = 0.80000000;
        break;
      case 65:
        descreen = 0.80000000;
        break;
      case 66:
        descreen = 0.80000000;
        break;
      case 67:
        descreen = 0.80000000;
        break;
      case 68:
        descreen = 0.80000000;
        break;
      case 69:
        descreen = 0.80000000;
        break;
      case 70:
        descreen = 0.80000000;
        break;
      case 71:
        descreen = 0.80000000;
        break;
      case 72:
        descreen = 1.04546757;
        break;
      case 73:
        descreen = 0.90634017;
        break;
      case 74:
        descreen = 0.30000000;
        break;
      case 75:
        descreen = 1.18211304;
        break;
      case 76:
        descreen = 1.19646933;
        break;
      case 77:
        descreen = 1.11828221;
        break;
      case 78:
        descreen = 1.82095370;
        break;
      case 79:
        descreen = 0.86706889;
        break;
      case 80:
        descreen = 1.26606799;
        break;
      case 81:
        descreen = 0.89753087;
        break;
      case 82:
        descreen = 0.95602208;
        break;
      case 83:
        descreen = 0.80335480;
        break;
      case 84:
        descreen = 0.80000000;
        break;
      case 85:
        descreen = 0.80000000;
        break;
      case 86:
        descreen = 0.80000000;
        break;
      case 87:
        descreen = 0.80000000;
        break;
      case 88:
        descreen = 0.80000000;
        break;
      case 89:
        descreen = 0.80000000;
        break;
      case 90:
        descreen = 0.80000000;
        break;
      case 91:
        descreen = 0.80000000;
        break;
      case 92:
        descreen = 0.80000000;
        break;
      case 93:
        descreen = 0.80000000;
        break;
      case 94:
        descreen = 0.80000000;
        break;
    }
  }
  else if (solvent == "water") {
    switch (atmnr) {
      case 1:
        descreen = 0.18678116;
        break;
      case 2:
        descreen = 1.99854836;
        break;
      case 3:
        descreen = 0.50934487;
        break;
      case 4:
        descreen = 0.30000000;
        break;
      case 5:
        descreen = 0.93372749;
        break;
      case 6:
        descreen = 0.73948749;
        break;
      case 7:
        descreen = 0.77003311;
        break;
      case 8:
        descreen = 0.30000000;
        break;
      case 9:
        descreen = 0.62227524;
        break;
      case 10:
        descreen = 1.22892076;
        break;
      case 11:
        descreen = 0.65065466;
        break;
      case 12:
        descreen = 1.15301804;
        break;
      case 13:
        descreen = 1.00517744;
        break;
      case 14:
        descreen = 0.93204996;
        break;
      case 15:
        descreen = 0.80000000;
        break;
      case 16:
        descreen = 0.30000000;
        break;
      case 17:
        descreen = 0.80000000;
        break;
      case 18:
        descreen = 0.90930774;
        break;
      case 19:
        descreen = 0.30000000;
        break;
      case 20:
        descreen = 1.13114016;
        break;
      case 21:
        descreen = 1.06981655;
        break;
      case 22:
        descreen = 1.42134411;
        break;
      case 23:
        descreen = 1.09146204;
        break;
      case 24:
        descreen = 1.14487481;
        break;
      case 25:
        descreen = 0.93238901;
        break;
      case 26:
        descreen = 1.66030803;
        break;
      case 27:
        descreen = 1.03734066;
        break;
      case 28:
        descreen = 1.21328994;
        break;
      case 29:
        descreen = 0.96459170;
        break;
      case 30:
        descreen = 0.70218231;
        break;
      case 31:
        descreen = 0.94458397;
        break;
      case 32:
        descreen = 0.87849331;
        break;
      case 33:
        descreen = 1.04104293;
        break;
      case 34:
        descreen = 1.01890919;
        break;
      case 35:
        descreen = 1.42471942;
        break;
      case 36:
        descreen = 1.31471665;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 0.99506730;
        break;
      case 39:
        descreen = 1.31750068;
        break;
      case 40:
        descreen = 0.91790577;
        break;
      case 41:
        descreen = 1.13352069;
        break;
      case 42:
        descreen = 0.90365194;
        break;
      case 43:
        descreen = 0.97192416;
        break;
      case 44:
        descreen = 1.13693166;
        break;
      case 45:
        descreen = 1.13908614;
        break;
      case 46:
        descreen = 1.03723586;
        break;
      case 47:
        descreen = 0.99373444;
        break;
      case 48:
        descreen = 0.57959099;
        break;
      case 49:
        descreen = 0.80138068;
        break;
      case 50:
        descreen = 0.45705909;
        break;
      case 51:
        descreen = 1.12709787;
        break;
      case 52:
        descreen = 0.80213876;
        break;
      case 53:
        descreen = 1.37998347;
        break;
      case 54:
        descreen = 1.10387200;
        break;
      case 55:
        descreen = 0.30000000;
        break;
      case 56:
        descreen = 0.30000000;
        break;
      case 57:
        descreen = 0.80000000;
        break;
      case 58:
        descreen = 0.80000000;
        break;
      case 59:
        descreen = 0.80000000;
        break;
      case 60:
        descreen = 0.80000000;
        break;
      case 61:
        descreen = 0.80000000;
        break;
      case 62:
        descreen = 0.80000000;
        break;
      case 63:
        descreen = 0.80000000;
        break;
      case 64:
        descreen = 0.80000000;
        break;
      case 65:
        descreen = 0.80000000;
        break;
      case 66:
        descreen = 0.80000000;
        break;
      case 67:
        descreen = 0.80000000;
        break;
      case 68:
        descreen = 0.80000000;
        break;
      case 69:
        descreen = 0.80000000;
        break;
      case 70:
        descreen = 0.80000000;
        break;
      case 71:
        descreen = 0.80000000;
        break;
      case 72:
        descreen = 1.11415978;
        break;
      case 73:
        descreen = 0.95102453;
        break;
      case 74:
        descreen = 1.31024463;
        break;
      case 75:
        descreen = 1.21914682;
        break;
      case 76:
        descreen = 0.91559959;
        break;
      case 77:
        descreen = 0.99579666;
        break;
      case 78:
        descreen = 0.90140483;
        break;
      case 79:
        descreen = 1.09420509;
        break;
      case 80:
        descreen = 0.02609846;
        break;
      case 81:
        descreen = 0.49510247;
        break;
      case 82:
        descreen = 0.69445869;
        break;
      case 83:
        descreen = 0.54304437;
        break;
      case 84:
        descreen = 0.80000000;
        break;
      case 85:
        descreen = 0.80000000;
        break;
      case 86:
        descreen = 0.80000000;
        break;
      case 87:
        descreen = 0.80000000;
        break;
      case 88:
        descreen = 0.80000000;
        break;
      case 89:
        descreen = 0.80000000;
        break;
      case 90:
        descreen = 0.80000000;
        break;
      case 91:
        descreen = 0.80000000;
        break;
      case 92:
        descreen = 0.80000000;
        break;
      case 93:
        descreen = 0.80000000;
        break;
      case 94:
        descreen = 0.80000000;
        break;
    }
  }
  else if (solvent == "octanol wet") {
    switch (atmnr) {
      case 1:
        descreen = 0.79740323;
        break;
      case 2:
        descreen = 1.42704421;
        break;
      case 3:
        descreen = 0.30000000;
        break;
      case 4:
        descreen = 0.84930829;
        break;
      case 5:
        descreen = 0.94053705;
        break;
      case 6:
        descreen = 0.71294102;
        break;
      case 7:
        descreen = 0.41315557;
        break;
      case 8:
        descreen = 0.10903752;
        break;
      case 9:
        descreen = 0.54083051;
        break;
      case 10:
        descreen = 1.49011806;
        break;
      case 11:
        descreen = 0.30000000;
        break;
      case 12:
        descreen = 1.13700779;
        break;
      case 13:
        descreen = 0.79561367;
        break;
      case 14:
        descreen = 0.81464854;
        break;
      case 15:
        descreen = 0.78522026;
        break;
      case 16:
        descreen = 0.94325024;
        break;
      case 17:
        descreen = 0.59762370;
        break;
      case 18:
        descreen = 0.87618008;
        break;
      case 19:
        descreen = 0.30000000;
        break;
      case 20:
        descreen = 1.25627726;
        break;
      case 21:
        descreen = 0.99040000;
        break;
      case 22:
        descreen = 1.05129906;
        break;
      case 23:
        descreen = 1.16850823;
        break;
      case 24:
        descreen = 1.21731707;
        break;
      case 25:
        descreen = 1.17468219;
        break;
      case 26:
        descreen = 1.11410026;
        break;
      case 27:
        descreen = 1.46827948;
        break;
      case 28:
        descreen = 1.12202311;
        break;
      case 29:
        descreen = 1.00062525;
        break;
      case 30:
        descreen = 1.18824199;
        break;
      case 31:
        descreen = 1.13090107;
        break;
      case 32:
        descreen = 0.88700962;
        break;
      case 33:
        descreen = 1.13804632;
        break;
      case 34:
        descreen = 1.23577058;
        break;
      case 35:
        descreen = 1.16618260;
        break;
      case 36:
        descreen = 0.95747171;
        break;
      case 37:
        descreen = 0.30000000;
        break;
      case 38:
        descreen = 0.59614255;
        break;
      case 39:
        descreen = 1.49375672;
        break;
      case 40:
        descreen = 0.95206613;
        break;
      case 41:
        descreen = 1.21021349;
        break;
      case 42:
        descreen = 1.56301734;
        break;
      case 43:
        descreen = 1.06430351;
        break;
      case 44:
        descreen = 1.38889074;
        break;
      case 45:
        descreen = 1.27424891;
        break;
      case 46:
        descreen = 1.13131865;
        break;
      case 47:
        descreen = 1.26716290;
        break;
      case 48:
        descreen = 1.18569309;
        break;
      case 49:
        descreen = 0.78154893;
        break;
      case 50:
        descreen = 0.83488723;
        break;
      case 51:
        descreen = 1.11082387;
        break;
      case 52:
        descreen = 0.07128553;
        break;
      case 53:
        descreen = 1.33670711;
        break;
      case 54:
        descreen = 0.86747398;
        break;
      case 55:
        descreen = 0.30000000;
        break;
      case 56:
        descreen = 0.99793477;
        break;
      case 57:
        descreen = 1.00000000;
        break;
      case 58:
        descreen = 1.00000000;
        break;
      case 59:
        descreen = 1.00000000;
        break;
      case 60:
        descreen = 1.00000000;
        break;
      case 61:
        descreen = 1.00000000;
        break;
      case 62:
        descreen = 1.00000000;
        break;
      case 63:
        descreen = 1.00000000;
        break;
      case 64:
        descreen = 1.00000000;
        break;
      case 65:
        descreen = 1.00000000;
        break;
      case 66:
        descreen = 1.00000000;
        break;
      case 67:
        descreen = 1.00000000;
        break;
      case 68:
        descreen = 1.00000000;
        break;
      case 69:
        descreen = 1.00000000;
        break;
      case 70:
        descreen = 1.00000000;
        break;
      case 71:
        descreen = 1.00000000;
        break;
      case 72:
        descreen = 1.28484313;
        break;
      case 73:
        descreen = 1.04263261;
        break;
      case 74:
        descreen = 1.31898601;
        break;
      case 75:
        descreen = 1.22064579;
        break;
      case 76:
        descreen = 1.07815523;
        break;
      case 77:
        descreen = 1.44458556;
        break;
      case 78:
        descreen = 1.16104008;
        break;
      case 79:
        descreen = 1.10531249;
        break;
      case 80:
        descreen = 1.32222214;
        break;
      case 81:
        descreen = 0.70207407;
        break;
      case 82:
        descreen = 0.98315389;
        break;
      case 83:
        descreen = 1.13356860;
        break;
      case 84:
        descreen = 1.00000000;
        break;
      case 85:
        descreen = 1.00000000;
        break;
      case 86:
        descreen = 1.00000000;
        break;
      case 87:
        descreen = 1.00000000;
        break;
      case 88:
        descreen = 1.00000000;
        break;
      case 89:
        descreen = 1.00000000;
        break;
      case 90:
        descreen = 1.00000000;
        break;
      case 91:
        descreen = 1.00000000;
        break;
      case 92:
        descreen = 1.00000000;
        break;
      case 93:
        descreen = 1.00000000;
        break;
      case 94:
        descreen = 1.00000000;
        break;
    }
  }
  return descreen;
}

#endif //_Solvation_Parameters_for_GFN2_
