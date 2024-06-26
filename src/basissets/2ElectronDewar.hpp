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

#ifndef _2_Electron_Functions_Dewar_
#define _2_Electron_Functions_Dewar_
#include <math.h>
#include "OvFunctions.hpp"
#include "../math/IntegralPackage.hpp"

//description:
//collection of functions for the calculation of two-center electron repulsion integrals for MNDO methods

std::vector<std::vector<int> > NonZeroMultipoleMom_sp(int LA, int MA, int LB, int MB) {
  //function that returns pairs of non-zero charges as adapted from 
  //Thiel, Voityuk, Theor. Chim. Acta, 81, 391, 1992
  //this is an auxiliary function for the calculation of two-eri in NDDO methods
  if ((MA > LA)||(MB > LB)) {throw("2Electron.hpp - NonZeroMultipoleMom_sp: wrong orbital definition");}
  //atom 1 should contain smallest L; if L1 == L2 then 1 contains lowest M
  int L1 = LA;
  int M1 = MA;
  int L2 = LB;
  int M2 = MB;
  if (L1 > L2) {
    L1 = LB;
    M1 = MB;
    L2 = LA;
    M2 = MA;
  }
  if (L1 == L2) {
    if (MB < MA) {
      M1 = MB;
      M2 = MA;
    }
  }
  std::vector<std::vector<int> > nonzero;
  std::vector<int> aux(2,0);
  if ((L1 == 0)&&(L2 == 0)) {
    //ss case
    nonzero.push_back(aux);
  }
  else if ((L1 == 0)&&(L2 == 1)) {
    //sp case
    aux[0] = 1;
    if (M2 == 0) {nonzero.push_back(aux);}
    else {
      aux[1] = 1*M2;               //multiply by M2 so that it takes signal
      nonzero.push_back(aux);
    }
  }
  else if ((L1 == 1)&&(L2 == 1)) {
    //pp case
    if ((M1 == 0)&&(M2 == 0)) {
      //sigma-sigma
      nonzero.push_back(aux);
      aux[0] = 2;
      nonzero.push_back(aux);
    }
    else if (((M1 == 0)&&(M2 == 1))||((M1 == 1)&&(M2 == 0))) {
      //sigma-pi, pi-sigma
      aux[0] = 2;
      aux[1] = 1;
      nonzero.push_back(aux);
    }
    else if (((M1 == 0)&&(M2 == -1))||((M1 == -1)&&(M2 == 0))) {
      //sigma-pib, pib-sigma
      aux[0] = 2;
      aux[1] = -1;
      nonzero.push_back(aux);
    }
    else if (((M1 == 1)&&(M2 == 1))||((M1 == -1)&&(M2 == -1))) {
      //pi-pi, pib-pib
      nonzero.push_back(aux);
      aux[0] = 2;
      aux[1] = 2*M2;               //multiply by M2 so that it takes signal
      nonzero.push_back(aux);
    }
    else if (((M1 == 1)&&(M2 == -1))||((M1 == -1)&&(M2 == 1))) {
      //pi-pib, pib-pi
      aux[0] = 2;
      aux[1] = 3;
      nonzero.push_back(aux);
    }
  }
  return nonzero;
}
std::vector<std::vector<int> > NonZeroMultipoleMom_d(int LA, int MA, int LB, int MB) {
  //function that returns pairs of non-zero charges as adapted from 
  //Thiel, Voityuk, Theor. Chim. Acta, 81, 391, 1992
  //this is an auxiliary function for the calculation of two-eri in NDDO methods
  if ((MA > LA)||(MB > LB)) {throw("2Electron.hpp - NonZeroMultipoleMom_d: wrong orbital definition");}
  //atom 1 should contain smallest L; if L1 == L2 then 1 contains lowest M
  int L1 = LA;
  int M1 = MA;
  int L2 = LB;
  int M2 = MB;
  if (L1 > L2) {
    L1 = LB;
    M1 = MB;
    L2 = LA;
    M2 = MA;
  }
  if (L1 == L2) {
    if (abs(MB) < abs(MA)) {
      M1 = MB;
      M2 = MA;
    }
  }
  std::vector<std::vector<int> > nonzero;
  std::vector<int> aux(2,0);
  if ((L1 == 0)&&(L2 == 0)) {
    //ss case
    nonzero.push_back(aux);
  }
  else if ((L1 == 0)&&(L2 == 1)) {
    //sp case
    aux[0] = 1;
    if (M2 == 0) {nonzero.push_back(aux);}
    else {
      aux[1] = M2;                 //multiply by M2 so that it takes signal
      nonzero.push_back(aux);
    }
  }
  else if ((L1 == 1)&&(L2 == 1)) {
    //pp case
    if ((M1 == 0)&&(M2 == 0)) {
      //sigma-sigma
      nonzero.push_back(aux);
      aux[0] = 2;
      nonzero.push_back(aux);
    }
    else if (((M1 == 0)&&(M2 == 1))||((M1 == 1)&&(M2 == 0))) {
      //sigma-pi, pi-sigma
      aux[0] = 2;
      aux[1] = 1;
      nonzero.push_back(aux);
    }
    else if (((M1 == 0)&&(M2 == -1))||((M1 == -1)&&(M2 == 0))) {
      //sigma-pib, pib-sigma
      aux[0] = 2;
      aux[1] = -1;
      nonzero.push_back(aux);
    }
    else if (((M1 == 1)&&(M2 == 1))||((M1 == -1)&&(M2 == -1))) {
      //pi-pi, pib-pib
      nonzero.push_back(aux);
      aux[0] = 2;
      nonzero.push_back(aux);
      aux[1] = 2;
      nonzero.push_back(aux);
    }
    else if (((M1 == 1)&&(M2 == -1))||((M1 == -1)&&(M2 == 1))) {
      //pi-pib, pib-pi
      aux[0] = 2;
      aux[1] = 2;
      nonzero.push_back(aux);
    }
  }
  else if ((L1 == 0)&&(L2 == 2)) {
    //sd case
    aux[0] = 2;
    if (M2 == 0) {nonzero.push_back(aux);}
    else {
      aux[1] = M2;                 //multiply by M2 so that it takes signal
      nonzero.push_back(aux);
    }
  }
  else if ((L1 == 1)&&(L2 == 2)) {
    //pd case
    aux[0] = 1;
    if ((M1 == 0)&&(M2 == 0)) {nonzero.push_back(aux);}     //sigma-sigma
    else if ((M1 == 0)&&(abs(M2) == 1)) {                   //sigma-pi
      aux[1] = M2;                 //multiply by M2 so that it takes signal
      nonzero.push_back(aux);
    }
    else if ((M1 == 0)&&(abs(M2) == 2)) {}                  //sigma-delta    -> do nothing
    else if ((abs(M1) == 1)&&(M2 == 0)) {                   //pi-sigma
      aux[1] = M1;                 //multiply by M1 so that it takes signal
      nonzero.push_back(aux);
    }
    else if ((abs(M1) == 1)&&(abs(M2) == 1)&&(M1 == M2)) {                   //pi-pi
      nonzero.push_back(aux);
    }
    else if ((abs(M1) == 1)&&(abs(M2) == 1)&&(M1 != M2)) {}                  //pi-pib
    else if ((M1 == 1)&&(M2 == 2)) {                                         //pi-delta
      aux[1] = 1;
      nonzero.push_back(aux);
    }
    else if ((M1 == 1)&&(M2 == -2)) {                                         //pi-deltab
      aux[1] = -1;
      nonzero.push_back(aux);
    }
    else if ((M1 == -1)&&(M2 == 2)) {                                         //pib-delta
      aux[1] = -1;
      nonzero.push_back(aux);
    }
    else if ((M1 == -1)&&(M2 == -2)) {                                         //pib-deltab
      aux[1] = 1;
      nonzero.push_back(aux);
    }
  }
  else if ((L1 == 2)&&(L2 == 2)) {
    //dd case
    if ((M1 == 0)&&(M2 == 0)) {                             //sigma-sigma
      nonzero.push_back(aux);
      aux[0] = 2;
      nonzero.push_back(aux);
    }
    else if ((M1 == 0)&&(abs(M2) == 1)) {                   //sigma-pi
      aux[0] = 2;
      aux[1] = M2;                 //multiply by M2 so that it takes signal
      nonzero.push_back(aux);
    }
    else if ((M1 == 0)&&(abs(M2) == 2)) {                   //sigma-delta
      aux[0] = 2;
      aux[1] = M2;                 //multiply by M2 so that it takes signal
      nonzero.push_back(aux);
    }
    else if ((abs(M1) == 1)&&(abs(M2) == 1)&&(M1 == M2)) {                   //pi-pi
      nonzero.push_back(aux);
      aux[0] = 2;
      nonzero.push_back(aux);
      aux[1] = 2;
      nonzero.push_back(aux);
    }
    else if ((abs(M1) == 1)&&(abs(M2) == 1)&&(M1 != M2)) {                   //pi-pib
      aux[0] = 2;
      aux[1] = 2;
      nonzero.push_back(aux);
    }
    else if ((M1 == 1)&&(abs(M2) == 2)) {                                    //pi-delta
      aux[0] = 2;
      aux[1] = M2/2;             //multiply by M2 so that it takes signal
      nonzero.push_back(aux);
    }
    else if ((M1 == -1)&&(abs(M2) == 2)) {                                   //pib-delta
      aux[0] = 2;
      aux[1] = -M2/2;
      nonzero.push_back(aux);
    }
    else if ((abs(M1) == 2)&&(abs(M2) == 2)&&(M1 == M2)) {                   //delta-delta
      nonzero.push_back(aux);
      aux[0] = 2;
      nonzero.push_back(aux);
    }
    else if ((abs(M1) == 2)&&(abs(M2) == 2)&&(M1 != M2)) {}                  //delta-deltab
  }
  return nonzero;
}
double se_multipole(int l1, int m1, int l2, int m2, double RAB, double aLL, std::vector<double> D) {
  //ERI auxiliary function in the calculation of Dewar's electron-repulsion integrals
  //Dewar, Thiel, Theor. Chim. Acta (Berl.), 46, 89, 1977
  //fuction giving back [l1m1,l2m2], where limi (i=1,2) is given by
  //00 -> q
  //10 -> {mu}_{z}
  //11 -> {mu}_{x}
  //20 -> {Q}_{zz}
  //21 -> {Q}_{xz}
  //22 -> {Q}_{xx}
  //2-2 -> {Q}_{yy}
  //23 -> {Q}_{xy}
  double DA1 = D[0];
  double DA2 = D[1];
  double DB1 = D[2];
  double DB2 = D[3];
  double multipole = 0.0;
  if ((l1 == 0)&&(l2 == 0)) {multipole = 1.0/sqrt(RAB*RAB + aLL*aLL);} 
  else if ((l1 == 0)&&(l2 == 1)&&(m2 == 0)) {multipole = -0.5*(1.0/sqrt((RAB + DB1)*(RAB + DB1) + aLL*aLL) - 1.0/sqrt((RAB - DB1)*(RAB - DB1) + aLL*aLL));}
  else if ((l1 == 1)&&(m1 == 0)&&(l2 == 0)) {
    DA1 = D[2];
    DA2 = D[3];
    DB1 = D[0];
    DB2 = D[1];
    multipole = 0.5*(1.0/sqrt((RAB + DB1)*(RAB + DB1) + aLL*aLL) - 1.0/sqrt((RAB - DB1)*(RAB - DB1) + aLL*aLL));
  }
  else if ((l1 == 0)&&(l2 == 2)&&(abs(m2) == 2)) {multipole = 0.5*(1.0/sqrt(RAB*RAB + 4.0*DB2*DB2 + aLL*aLL) - 1.0/sqrt(RAB*RAB + aLL*aLL));}
  else if ((l1 == 2)&&(abs(m1) == 2)&&(l2 == 0)) {
    DA1 = D[2];
    DA2 = D[3];
    DB1 = D[0];
    DB2 = D[1];
    multipole = 0.5*(1.0/sqrt(RAB*RAB + 4.0*DB2*DB2 + aLL*aLL) - 1.0/sqrt(RAB*RAB + aLL*aLL));
  }
  else if ((l1 == 0)&&(l2 == 2)&&(m2 == 0)) {multipole = 0.25*(1.0/sqrt((RAB + 2.0*DB2)*(RAB + 2.0*DB2) + aLL*aLL) - 2.0/sqrt(RAB*RAB + aLL*aLL) + 1.0/sqrt((RAB - 2.0*DB2)*(RAB - 2.0*DB2) + aLL*aLL));}
  else if ((l1 == 2)&&(m1 == 0)&&(l2 == 0)) {
    DA1 = D[2];
    DA2 = D[3];
    DB1 = D[0];
    DB2 = D[1];
    multipole = 0.25*(1.0/sqrt((RAB + 2.0*DB2)*(RAB + 2.0*DB2) + aLL*aLL) - 2.0/sqrt(RAB*RAB + aLL*aLL) + 1.0/sqrt((RAB - 2.0*DB2)*(RAB - 2.0*DB2) + aLL*aLL));
  }
  else if ((l1 == 1)&&(abs(m1) == 1)&&(l2 == 1)&&(abs(m2) == 1)) {multipole = 0.5*(1.0/sqrt(RAB*RAB + (DA1 - DB1)*(DA1 - DB1) + aLL*aLL) - 1.0/sqrt(RAB*RAB + (DA1 + DB1)*(DA1 + DB1) + aLL*aLL));}
  else if ((l1 == 1)&&(m1 == 0)&&(l2 == 1)&&(m2 == 0)) {multipole = 0.25*(1.0/sqrt((RAB + DA1 - DB1)*(RAB + DA1 - DB1) + aLL*aLL) - 1.0/sqrt((RAB + DA1 + DB1)*(RAB + DA1 + DB1) + aLL*aLL) - 1.0/sqrt((RAB - DA1 - DB1)*(RAB - DA1 - DB1) + aLL*aLL) + 1.0/sqrt((RAB - DA1 + DB1)*(RAB - DA1 + DB1) + aLL*aLL));}
  else if ((l1 == 1)&&(abs(m1) == 1)&&(l2 == 2)&&(m2 == 1)) {multipole = -0.25*(-1.0/sqrt((RAB - DB2)*(RAB - DB2) + (DA1 - DB2)*(DA1 - DB2) + aLL*aLL) + 1.0/sqrt((RAB - DB2)*(RAB - DB2) + (DA1 + DB2)*(DA1 + DB2) + aLL*aLL) + 1.0/sqrt((RAB + DB2)*(RAB + DB2) + (DA1 - DB2)*(DA1 - DB2) + aLL*aLL) - 1.0/sqrt((RAB + DB2)*(RAB + DB2) + (DA1 + DB2)*(DA1 + DB2) + aLL*aLL));}
  else if ((l1 == 2)&&(m1 == 1)&&(l2 == 1)&&(abs(m2) == 1)) {
    DA1 = D[2];
    DA2 = D[3];
    DB1 = D[0];
    DB2 = D[1];
    multipole = 0.25*(-1.0/sqrt((RAB - DB2)*(RAB - DB2) + (DA1 - DB2)*(DA1 - DB2) + aLL*aLL) + 1.0/sqrt((RAB - DB2)*(RAB - DB2) + (DA1 + DB2)*(DA1 + DB2) + aLL*aLL) + 1.0/sqrt((RAB + DB2)*(RAB + DB2) + (DA1 - DB2)*(DA1 - DB2) + aLL*aLL) - 1.0/sqrt((RAB + DB2)*(RAB + DB2) + (DA1 + DB2)*(DA1 + DB2) + aLL*aLL));
  }
  else if ((l1 == 1)&&(m1 == 0)&&(l2 == 2)&&(abs(m2) == 2)) {multipole = -0.25*(-1.0/sqrt((RAB + DA1)*(RAB + DA1) + 4.0*DB2*DB2 + aLL*aLL) + 1.0/sqrt((RAB - DA1)*(RAB - DA1) + 4.0*DB2*DB2 + aLL*aLL) + 1.0/sqrt((RAB + DA1)*(RAB + DA1) + aLL*aLL) - 1.0/sqrt((RAB - DA1)*(RAB - DA1) + aLL*aLL));}
  else if ((l1 == 2)&&(abs(m1) == 2)&&(l2 == 1)&&(m2 == 0)) {
    DA1 = D[2];
    DA2 = D[3];
    DB1 = D[0];
    DB2 = D[1];
    multipole = 0.25*(-1.0/sqrt((RAB + DA1)*(RAB + DA1) + 4.0*DB2*DB2 + aLL*aLL) + 1.0/sqrt((RAB - DA1)*(RAB - DA1) + 4.0*DB2*DB2 + aLL*aLL) + 1.0/sqrt((RAB + DA1)*(RAB + DA1) + aLL*aLL) - 1.0/sqrt((RAB - DA1)*(RAB - DA1) + aLL*aLL));
  }
  else if ((l1 == 1)&&(m1 == 0)&&(l2 == 2)&&(m2 == 0)) {multipole = -0.125*(-1.0/sqrt((RAB + DA1 - 2.0*DB2)*(RAB + DA1 - 2.0*DB2) + aLL*aLL) + 1.0/sqrt((RAB - DA1 - 2.0*DB2)*(RAB - DA1 - 2.0*DB2) + aLL*aLL) - 1.0/sqrt((RAB + DA1 + 2.0*DB2)*(RAB + DA1 + 2.0*DB2) + aLL*aLL) + 1.0/sqrt((RAB - DA1 + 2.0*DB2)*(RAB - DA1 + 2.0*DB2) + aLL*aLL) + 2.0/sqrt((RAB + DA1)*(RAB + DA1) + aLL*aLL) - 2.0/sqrt((RAB - DA1)*(RAB - DA1) + aLL*aLL));}
  else if ((l1 == 2)&&(m1 == 0)&&(l2 == 1)&&(m2 == 0)) {
    DA1 = D[2];
    DA2 = D[3];
    DB1 = D[0];
    DB2 = D[1];
    multipole = 0.125*(-1.0/sqrt((RAB + DA1 - 2.0*DB2)*(RAB + DA1 - 2.0*DB2) + aLL*aLL) + 1.0/sqrt((RAB - DA1 - 2.0*DB2)*(RAB - DA1 - 2.0*DB2) + aLL*aLL) - 1.0/sqrt((RAB + DA1 + 2.0*DB2)*(RAB + DA1 + 2.0*DB2) + aLL*aLL) + 1.0/sqrt((RAB - DA1 + 2.0*DB2)*(RAB - DA1 + 2.0*DB2) + aLL*aLL) + 2.0/sqrt((RAB + DA1)*(RAB + DA1) + aLL*aLL) - 2.0/sqrt((RAB - DA1)*(RAB - DA1) + aLL*aLL));
  }
  else if ((l1 == 2)&&(abs(m1) == 2)&&(l2 == 2)&&(abs(m2) == 2)&&(m1 == m2)) {multipole = 0.125*(1.0/sqrt(RAB*RAB + 4.0*(DA2 - DB2)*(DA2 - DB2) + aLL*aLL) + 1.0/sqrt(RAB*RAB + 4.0*(DA2 + DB2)*(DA2 + DB2) + aLL*aLL) - 2.0/sqrt(RAB*RAB + 4.0*DA2*DA2 + aLL*aLL) - 2.0/sqrt(RAB*RAB + 4.0*DB2*DB2 + aLL*aLL) + 2.0/sqrt(RAB*RAB + aLL*aLL));}
  else if ((l1 == 2)&&(m1 == 2)&&(l2 == 2)&&(m2 == -2)) {multipole = 0.25*(1.0/sqrt(RAB*RAB + 4.0*DA2*DA2 + 4.0*DB2*DB2 + aLL*aLL) - 1.0/sqrt(RAB*RAB + 4.0*DA2*DA2 + aLL*aLL) - 1.0/sqrt(RAB*RAB + 4.0*DB2*DB2 + aLL*aLL) + 1.0/sqrt(RAB*RAB + aLL*aLL));}
  else if ((l1 == 2)&&(m1 == -2)&&(l2 == 2)&&(m2 == 2)) {multipole = 0.25*(1.0/sqrt(RAB*RAB + 4.0*DA2*DA2 + 4.0*DB2*DB2 + aLL*aLL) - 1.0/sqrt(RAB*RAB + 4.0*DA2*DA2 + aLL*aLL) - 1.0/sqrt(RAB*RAB + 4.0*DB2*DB2 + aLL*aLL) + 1.0/sqrt(RAB*RAB + aLL*aLL));}
  else if ((l1 == 2)&&(abs(m1) == 2)&&(l2 == 2)&&(m2 == 0)) {multipole = 0.125*(1.0/sqrt((RAB - 2.0*DB2)*(RAB - 2.0*DB2) + 4.0*DA2*DA2 + aLL*aLL) + 1.0/sqrt((RAB + 2.0*DB2)*(RAB + 2.0*DB2) + 4.0*DA2*DA2 + aLL*aLL) - 1.0/sqrt((RAB - 2.0*DB2)*(RAB - 2.0*DB2) + aLL*aLL) - 1.0/sqrt((RAB + 2.0*DB2)*(RAB + 2.0*DB2) + aLL*aLL) - 2.0/sqrt(RAB*RAB + 4.0*DA2*DA2 + aLL*aLL) + 2.0/sqrt(RAB*RAB + aLL*aLL));}
  else if ((l1 == 2)&&(m1 == 0)&&(l2 == 2)&&(abs(m2) == 2)) {
    DA1 = D[2];
    DA2 = D[3];
    DB1 = D[0];
    DB2 = D[1];
    multipole = 0.125*(1.0/sqrt((RAB - 2.0*DB2)*(RAB - 2.0*DB2) + 4.0*DA2*DA2 + aLL*aLL) + 1.0/sqrt((RAB + 2.0*DB2)*(RAB + 2.0*DB2) + 4.0*DA2*DA2 + aLL*aLL) - 1.0/sqrt((RAB - 2.0*DB2)*(RAB - 2.0*DB2) + aLL*aLL) - 1.0/sqrt((RAB + 2.0*DB2)*(RAB + 2.0*DB2) + aLL*aLL) - 2.0/sqrt(RAB*RAB + 4.0*DA2*DA2 + aLL*aLL) + 2.0/sqrt(RAB*RAB + aLL*aLL));
  }
  else if ((l1 == 2)&&(m1 == 0)&&(l2 == 2)&&(m2 == 0)) {multipole = 0.0625*(1.0/sqrt((RAB + 2.0*DA2 - 2.0*DB2)*(RAB + 2.0*DA2 - 2.0*DB2) + aLL*aLL) + 1.0/sqrt((RAB + 2.0*DA2 + 2.0*DB2)*(RAB + 2.0*DA2 + 2.0*DB2) + aLL*aLL) + 1.0/sqrt((RAB - 2.0*DA2 - 2.0*DB2)*(RAB - 2.0*DA2 - 2.0*DB2) + aLL*aLL) + 1.0/sqrt((RAB - 2.0*DA2 + 2.0*DB2)*(RAB - 2.0*DA2 + 2.0*DB2) + aLL*aLL) - 2.0/sqrt((RAB + 2.0*DA2)*(RAB + 2.0*DA2) + aLL*aLL) - 2.0/sqrt((RAB - 2.0*DA2)*(RAB - 2.0*DA2) + aLL*aLL) - 2.0/sqrt((RAB + 2.0*DB2)*(RAB + 2.0*DB2) + aLL*aLL) - 2.0/sqrt((RAB - 2.0*DB2)*(RAB - 2.0*DB2) + aLL*aLL) + 4.0/sqrt(RAB*RAB + aLL*aLL));}
  else if ((l1 == 2)&&(m1 == 1)&&(l2 == 2)&&(m2 == 1)) {multipole = 0.125*(1.0/sqrt((RAB + DA2 - DB2)*(RAB + DA2 - DB2) + (DA2 - DB2)*(DA2 - DB2) + aLL*aLL) - 1.0/sqrt((RAB + DA2 - DB2)*(RAB + DA2 - DB2) + (DA2 + DB2)*(DA2 + DB2) + aLL*aLL) - 1.0/sqrt((RAB + DA2 + DB2)*(RAB + DA2 + DB2) + (DA2 - DB2)*(DA2 - DB2) + aLL*aLL) + 1.0/sqrt((RAB + DA2 + DB2)*(RAB + DA2 + DB2) + (DA2 + DB2)*(DA2 + DB2) + aLL*aLL) - 1.0/sqrt((RAB - DA2 - DB2)*(RAB - DA2 - DB2) + (DA2 - DB2)*(DA2 - DB2) + aLL*aLL) + 1.0/sqrt((RAB - DA2 - DB2)*(RAB - DA2 - DB2) + (DA2 + DB2)*(DA2 + DB2) + aLL*aLL) + 1.0/sqrt((RAB - DA2 + DB2)*(RAB - DA2 + DB2) + (DA2 - DB2)*(DA2 - DB2) + aLL*aLL) - 1.0/sqrt((RAB - DA2 + DB2)*(RAB - DA2 + DB2) + (DA2 + DB2)*(DA2 + DB2) + aLL*aLL));}
  else if ((l1 == 2)&&(m1 == 3)&&(l2 == 2)&&(m2 == 3)) {multipole = 0.25*(1.0/sqrt(RAB*RAB + 2.0*(DA2 - DB2)*(DA2 - DB2) + aLL*aLL) + 1.0/sqrt(RAB*RAB + 2.0*(DA2 + DB2)*(DA2 + DB2) + aLL*aLL) - 2.0/sqrt(RAB*RAB + 2.0*DA2*DA2 + 2.0*DB2*DB2 + aLL*aLL));}
  return multipole;
}
double se_multipole_dR(int l1, int m1, int l2, int m2, double RAB, double aLL, std::vector<double> D) {
  //auxiliary function in the calculation of the derivative of ERIs for Dewar's electron-repulsion integrals 
  //with respect to internuclear distance
  //fuction giving back (d/dR)[l1m1,l2m2], where limi (i=1,2) is given by
  //00 -> q
  //10 -> {mu}_{z}
  //11 -> {mu}_{x}
  //20 -> {Q}_{zz}
  //21 -> {Q}_{xz}
  //22 -> {Q}_{xx}
  //2-2 -> {Q}_{yy}
  //23 -> {Q}_{xy}
  double DA1 = D[0];
  double DA2 = D[1];
  double DB1 = D[2];
  double DB2 = D[3];
  double gmultipole = 0.0;
  double aux;
  if ((l1 == 0)&&(l2 == 0)) {
    aux = sqrt(RAB*RAB + aLL*aLL);
    gmultipole = -RAB/(aux*aux*aux);
  } 
  else if ((l1 == 0)&&(l2 == 1)&&(m2 == 0)) {
    aux = sqrt((RAB + DB1)*(RAB + DB1) + aLL*aLL);
    gmultipole = 0.5*(RAB + DB1)/(aux*aux*aux);
    aux = sqrt((RAB - DB1)*(RAB - DB1) + aLL*aLL);
    gmultipole -= 0.5*(RAB - DB1)/(aux*aux*aux);
  }
  else if ((l1 == 1)&&(m1 == 0)&&(l2 == 0)) {
    DA1 = D[2];
    DA2 = D[3];
    DB1 = D[0];
    DB2 = D[1];
    aux = sqrt((RAB + DB1)*(RAB + DB1) + aLL*aLL);
    gmultipole = -0.5*(RAB + DB1)/(aux*aux*aux);
    aux = sqrt((RAB - DB1)*(RAB - DB1) + aLL*aLL);
    gmultipole += 0.5*(RAB - DB1)/(aux*aux*aux);
  }
  else if ((l1 == 0)&&(l2 == 2)&&(abs(m2) == 2)) {
    aux = sqrt(RAB*RAB + 4.0*DB2*DB2 + aLL*aLL);
    gmultipole = -0.5*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + aLL*aLL);
    gmultipole += 0.5*RAB/(aux*aux*aux);
  }
  else if ((l1 == 2)&&(abs(m1) == 2)&&(l2 == 0)) {
    DA1 = D[2];
    DA2 = D[3];
    DB1 = D[0];
    DB2 = D[1];
    aux = sqrt(RAB*RAB + 4.0*DB2*DB2 + aLL*aLL);
    gmultipole = -0.5*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + aLL*aLL);
    gmultipole += 0.5*RAB/(aux*aux*aux);
  }
  else if ((l1 == 0)&&(l2 == 2)&&(m2 == 0)) {
    aux = sqrt((RAB + 2.0*DB2)*(RAB + 2.0*DB2) + aLL*aLL);
    gmultipole = -0.25*(RAB + 2.0*DB2)/(aux*aux*aux);
    aux = sqrt(RAB*RAB + aLL*aLL);
    gmultipole += 0.5*RAB/(aux*aux*aux);
    aux = sqrt((RAB - 2.0*DB2)*(RAB - 2.0*DB2) + aLL*aLL);
    gmultipole -= 0.25*(RAB - 2.0*DB2)/(aux*aux*aux);
  }
  else if ((l1 == 2)&&(m1 == 0)&&(l2 == 0)) {
    DA1 = D[2];
    DA2 = D[3];
    DB1 = D[0];
    DB2 = D[1];
    aux = sqrt((RAB + 2.0*DB2)*(RAB + 2.0*DB2) + aLL*aLL);
    gmultipole = -0.25*(RAB + 2.0*DB2)/(aux*aux*aux);
    aux = sqrt(RAB*RAB + aLL*aLL);
    gmultipole += 0.5*RAB/(aux*aux*aux);
    aux = sqrt((RAB - 2.0*DB2)*(RAB - 2.0*DB2) + aLL*aLL);
    gmultipole -= 0.25*(RAB - 2.0*DB2)/(aux*aux*aux);
  }
  else if ((l1 == 1)&&(abs(m1) == 1)&&(l2 == 1)&&(abs(m2) == 1)) {
    aux = sqrt(RAB*RAB + (DA1 - DB1)*(DA1 - DB1) + aLL*aLL);
    gmultipole = -0.5*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + (DA1 + DB1)*(DA1 + DB1) + aLL*aLL);
    gmultipole += 0.5*RAB/(aux*aux*aux);
  }
  else if ((l1 == 1)&&(m1 == 0)&&(l2 == 1)&&(m2 == 0)) {
    aux = sqrt((RAB + DA1 - DB1)*(RAB + DA1 - DB1) + aLL*aLL);
    gmultipole = -0.25*(RAB + DA1 - DB1)/(aux*aux*aux);
    aux = sqrt((RAB + DA1 + DB1)*(RAB + DA1 + DB1) + aLL*aLL);
    gmultipole += 0.25*(RAB + DA1 + DB1)/(aux*aux*aux);
    aux = sqrt((RAB - DA1 - DB1)*(RAB - DA1 - DB1) + aLL*aLL);
    gmultipole += 0.25*(RAB - DA1 - DB1)/(aux*aux*aux);
    aux = sqrt((RAB - DA1 + DB1)*(RAB - DA1 + DB1) + aLL*aLL);
    gmultipole -= 0.25*(RAB - DA1 + DB1)/(aux*aux*aux);
  }
  else if ((l1 == 1)&&(abs(m1) == 1)&&(l2 == 2)&&(m2 == 1)) {
    aux = sqrt((RAB - DB2)*(RAB - DB2) + (DA1 - DB2)*(DA1 - DB2) + aLL*aLL);
    gmultipole = -0.25*(RAB - DB2)/(aux*aux*aux);
    aux = sqrt((RAB - DB2)*(RAB - DB2) + (DA1 + DB2)*(DA1 + DB2) + aLL*aLL);
    gmultipole += 0.25*(RAB - DB2)/(aux*aux*aux);
    aux = sqrt((RAB + DB2)*(RAB + DB2) + (DA1 - DB2)*(DA1 - DB2) + aLL*aLL);
    gmultipole += 0.25*(RAB + DB2)/(aux*aux*aux);
    aux = sqrt((RAB + DB2)*(RAB + DB2) + (DA1 + DB2)*(DA1 + DB2) + aLL*aLL);
    gmultipole -= 0.25*(RAB + DB2)/(aux*aux*aux);
  }
  else if ((l1 == 2)&&(m1 == 1)&&(l2 == 1)&&(abs(m2) == 1)) {
    DA1 = D[2];
    DA2 = D[3];
    DB1 = D[0];
    DB2 = D[1];
    aux = sqrt((RAB - DB2)*(RAB - DB2) + (DA1 - DB2)*(DA1 - DB2) + aLL*aLL);
    gmultipole = 0.25*(RAB - DB2)/(aux*aux*aux);
    aux = sqrt((RAB - DB2)*(RAB - DB2) + (DA1 + DB2)*(DA1 + DB2) + aLL*aLL);
    gmultipole -= 0.25*(RAB - DB2)/(aux*aux*aux);
    aux = sqrt((RAB + DB2)*(RAB + DB2) + (DA1 - DB2)*(DA1 - DB2) + aLL*aLL);
    gmultipole -= 0.25*(RAB + DB2)/(aux*aux*aux);
    aux = sqrt((RAB + DB2)*(RAB + DB2) + (DA1 + DB2)*(DA1 + DB2) + aLL*aLL);
    gmultipole += 0.25*(RAB + DB2)/(aux*aux*aux);
  }
  else if ((l1 == 1)&&(m1 == 0)&&(l2 == 2)&&(abs(m2) == 2)) {
    aux = sqrt((RAB + DA1)*(RAB + DA1) + 4.0*DB2*DB2 + aLL*aLL);
    gmultipole = -0.25*(RAB + DA1)/(aux*aux*aux);
    aux = sqrt((RAB - DA1)*(RAB - DA1) + 4.0*DB2*DB2 + aLL*aLL);
    gmultipole += 0.25*(RAB - DA1)/(aux*aux*aux);
    aux = sqrt((RAB + DA1)*(RAB + DA1) + aLL*aLL);
    gmultipole += 0.25*(RAB + DA1)/(aux*aux*aux);
    aux = sqrt((RAB - DA1)*(RAB - DA1) + aLL*aLL);
    gmultipole -= 0.25*(RAB - DA1)/(aux*aux*aux);
  }
  else if ((l1 == 2)&&(abs(m1) == 2)&&(l2 == 1)&&(m2 == 0)) {
    DA1 = D[2];
    DA2 = D[3];
    DB1 = D[0];
    DB2 = D[1];
    aux = sqrt((RAB + DA1)*(RAB + DA1) + 4.0*DB2*DB2 + aLL*aLL);
    gmultipole = 0.25*(RAB + DA1)/(aux*aux*aux);
    aux = sqrt((RAB - DA1)*(RAB - DA1) + 4.0*DB2*DB2 + aLL*aLL);
    gmultipole -= 0.25*(RAB - DA1)/(aux*aux*aux);
    aux = sqrt((RAB + DA1)*(RAB + DA1) + aLL*aLL);
    gmultipole -= 0.25*(RAB + DA1)/(aux*aux*aux);
    aux = sqrt((RAB - DA1)*(RAB - DA1) + aLL*aLL);
    gmultipole += 0.25*(RAB - DA1)/(aux*aux*aux);
  }
  else if ((l1 == 1)&&(m1 == 0)&&(l2 == 2)&&(m2 == 0)) {
    aux = sqrt((RAB + DA1 - 2.0*DB2)*(RAB + DA1 - 2.0*DB2) + aLL*aLL);
    gmultipole = -0.125*(RAB + DA1 - 2.0*DB2)/(aux*aux*aux);
    aux = sqrt((RAB - DA1 - 2.0*DB2)*(RAB - DA1 - 2.0*DB2) + aLL*aLL);
    gmultipole += 0.125*(RAB - DA1 - 2.0*DB2)/(aux*aux*aux);
    aux = sqrt((RAB + DA1 + 2.0*DB2)*(RAB + DA1 + 2.0*DB2) + aLL*aLL);
    gmultipole -= 0.125*(RAB + DA1 + 2.0*DB2)/(aux*aux*aux);
    aux = sqrt((RAB - DA1 + 2.0*DB2)*(RAB - DA1 + 2.0*DB2) + aLL*aLL);
    gmultipole += 0.125*(RAB - DA1 + 2.0*DB2)/(aux*aux*aux);
    aux = sqrt((RAB + DA1)*(RAB + DA1) + aLL*aLL);
    gmultipole += 0.25*(RAB + DA1)/(aux*aux*aux);
    aux = sqrt((RAB - DA1)*(RAB - DA1) + aLL*aLL);
    gmultipole -= 0.25*(RAB - DA1)/(aux*aux*aux);
  }
  else if ((l1 == 2)&&(m1 == 0)&&(l2 == 1)&&(m2 == 0)) {
    DA1 = D[2];
    DA2 = D[3];
    DB1 = D[0];
    DB2 = D[1];
    aux = sqrt((RAB + DA1 - 2.0*DB2)*(RAB + DA1 - 2.0*DB2) + aLL*aLL);
    gmultipole = 0.125*(RAB + DA1 - 2.0*DB2)/(aux*aux*aux);
    aux = sqrt((RAB - DA1 - 2.0*DB2)*(RAB - DA1 - 2.0*DB2) + aLL*aLL);
    gmultipole -= 0.125*(RAB - DA1 - 2.0*DB2)/(aux*aux*aux);
    aux = sqrt((RAB + DA1 + 2.0*DB2)*(RAB + DA1 + 2.0*DB2) + aLL*aLL);
    gmultipole += 0.125*(RAB + DA1 + 2.0*DB2)/(aux*aux*aux);
    aux = sqrt((RAB - DA1 + 2.0*DB2)*(RAB - DA1 + 2.0*DB2) + aLL*aLL);
    gmultipole -= 0.125*(RAB - DA1 + 2.0*DB2)/(aux*aux*aux);
    aux = sqrt((RAB + DA1)*(RAB + DA1) + aLL*aLL);
    gmultipole -= 0.25*(RAB + DA1)/(aux*aux*aux);
    aux = sqrt((RAB - DA1)*(RAB - DA1) + aLL*aLL);
    gmultipole += 0.25*(RAB - DA1)/(aux*aux*aux);
  }
  else if ((l1 == 2)&&(abs(m1) == 2)&&(l2 == 2)&&(abs(m2) == 2)&&(m1 == m2)) {
    aux = sqrt(RAB*RAB + 4.0*(DA2 - DB2)*(DA2 - DB2) + aLL*aLL);
    gmultipole = -0.125*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + 4.0*(DA2 + DB2)*(DA2 + DB2) + aLL*aLL);
    gmultipole -= 0.125*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + 4.0*DA2*DA2 + aLL*aLL);
    gmultipole += 0.25*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + 4.0*DB2*DB2 + aLL*aLL);
    gmultipole += 0.25*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + aLL*aLL);
    gmultipole -= 0.25*RAB/(aux*aux*aux);
  }
  else if ((l1 == 2)&&(m1 == 2)&&(l2 == 2)&&(m2 == -2)) {
    aux = sqrt(RAB*RAB + 4.0*DA2*DA2 + 4.0*DB2*DB2 + aLL*aLL);
    gmultipole = -0.25*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + 4.0*DA2*DA2 + aLL*aLL);
    gmultipole += 0.25*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + 4.0*DB2*DB2 + aLL*aLL);
    gmultipole += 0.25*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + aLL*aLL);
    gmultipole -= 0.25*RAB/(aux*aux*aux);
  }
  else if ((l1 == 2)&&(m1 == -2)&&(l2 == 2)&&(m2 == 2)) {
    aux = sqrt(RAB*RAB + 4.0*DA2*DA2 + 4.0*DB2*DB2 + aLL*aLL);
    gmultipole = -0.25*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + 4.0*DA2*DA2 + aLL*aLL);
    gmultipole += 0.25*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + 4.0*DB2*DB2 + aLL*aLL);
    gmultipole += 0.25*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + aLL*aLL);
    gmultipole -= 0.25*RAB/(aux*aux*aux);
  }
  else if ((l1 == 2)&&(abs(m1) == 2)&&(l2 == 2)&&(m2 == 0)) {
    aux = sqrt((RAB - 2.0*DB2)*(RAB - 2.0*DB2) + 4.0*DA2*DA2 + aLL*aLL);
    gmultipole = -0.125*(RAB - 2.0*DB2)/(aux*aux*aux);
    aux = sqrt((RAB + 2.0*DB2)*(RAB + 2.0*DB2) + 4.0*DA2*DA2 + aLL*aLL);
    gmultipole -= 0.125*(RAB + 2.0*DB2)/(aux*aux*aux);
    aux = sqrt((RAB - 2.0*DB2)*(RAB - 2.0*DB2) + aLL*aLL);
    gmultipole += 0.125*(RAB - 2.0*DB2)/(aux*aux*aux);
    aux = sqrt((RAB + 2.0*DB2)*(RAB + 2.0*DB2) + aLL*aLL);
    gmultipole += 0.125*(RAB + 2.0*DB2)/(aux*aux*aux);
    aux = sqrt(RAB*RAB + 4.0*DA2*DA2 + aLL*aLL);
    gmultipole += 0.25*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + aLL*aLL);
    gmultipole -= 0.25*RAB/(aux*aux*aux);
  }
  else if ((l1 == 2)&&(m1 == 0)&&(l2 == 2)&&(abs(m2) == 2)) {
    DA1 = D[2];
    DA2 = D[3];
    DB1 = D[0];
    DB2 = D[1];
    aux = sqrt((RAB - 2.0*DB2)*(RAB - 2.0*DB2) + 4.0*DA2*DA2 + aLL*aLL);
    gmultipole = -0.125*(RAB - 2.0*DB2)/(aux*aux*aux);
    aux = sqrt((RAB + 2.0*DB2)*(RAB + 2.0*DB2) + 4.0*DA2*DA2 + aLL*aLL);
    gmultipole -= 0.125*(RAB + 2.0*DB2)/(aux*aux*aux);
    aux = sqrt((RAB - 2.0*DB2)*(RAB - 2.0*DB2) + aLL*aLL);
    gmultipole += 0.125*(RAB - 2.0*DB2)/(aux*aux*aux);
    aux = sqrt((RAB + 2.0*DB2)*(RAB + 2.0*DB2) + aLL*aLL);
    gmultipole += 0.125*(RAB + 2.0*DB2)/(aux*aux*aux);
    aux = sqrt(RAB*RAB + 4.0*DA2*DA2 + aLL*aLL);
    gmultipole += 0.25*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + aLL*aLL);
    gmultipole -= 0.25*RAB/(aux*aux*aux);
  }
  else if ((l1 == 2)&&(m1 == 0)&&(l2 == 2)&&(m2 == 0)) {
    aux = sqrt((RAB + 2.0*DA2 - 2.0*DB2)*(RAB + 2.0*DA2 - 2.0*DB2) + aLL*aLL);
    gmultipole = -0.0625*(RAB + 2.0*DA2 - 2.0*DB2)/(aux*aux*aux);
    aux = sqrt((RAB + 2.0*DA2 + 2.0*DB2)*(RAB + 2.0*DA2 + 2.0*DB2) + aLL*aLL);
    gmultipole -= 0.0625*(RAB + 2.0*DA2 + 2.0*DB2)/(aux*aux*aux);
    aux = sqrt((RAB - 2.0*DA2 - 2.0*DB2)*(RAB - 2.0*DA2 - 2.0*DB2) + aLL*aLL);
    gmultipole -= 0.0625*(RAB - 2.0*DA2 - 2.0*DB2)/(aux*aux*aux);
    aux = sqrt((RAB - 2.0*DA2 + 2.0*DB2)*(RAB - 2.0*DA2 + 2.0*DB2) + aLL*aLL);
    gmultipole -= 0.0625*(RAB - 2.0*DA2 + 2.0*DB2)/(aux*aux*aux);
    aux = sqrt((RAB + 2.0*DA2)*(RAB + 2.0*DA2) + aLL*aLL);
    gmultipole += 0.125*(RAB + 2.0*DA2)/(aux*aux*aux);
    aux = sqrt((RAB - 2.0*DA2)*(RAB - 2.0*DA2) + aLL*aLL);
    gmultipole += 0.125*(RAB - 2.0*DA2)/(aux*aux*aux);
    aux = sqrt((RAB + 2.0*DB2)*(RAB + 2.0*DB2) + aLL*aLL);
    gmultipole += 0.125*(RAB + 2.0*DB2)/(aux*aux*aux);
    aux = sqrt((RAB - 2.0*DB2)*(RAB - 2.0*DB2) + aLL*aLL);
    gmultipole += 0.125*(RAB - 2.0*DB2)/(aux*aux*aux);
    aux = sqrt(RAB*RAB + aLL*aLL);
    gmultipole -= 0.25*RAB/(aux*aux*aux);
  }
  else if ((l1 == 2)&&(m1 == 1)&&(l2 == 2)&&(m2 == 1)) {
    aux = sqrt((RAB + DA2 - DB2)*(RAB + DA2 - DB2) + (DA2 - DB2)*(DA2 - DB2) + aLL*aLL);
    gmultipole = -0.125*(RAB + DA2 - DB2)/(aux*aux*aux);
    aux = sqrt((RAB + DA2 - DB2)*(RAB + DA2 - DB2) + (DA2 + DB2)*(DA2 + DB2) + aLL*aLL);
    gmultipole += 0.125*(RAB + DA2 - DB2)/(aux*aux*aux);
    aux = sqrt((RAB + DA2 + DB2)*(RAB + DA2 + DB2) + (DA2 - DB2)*(DA2 - DB2) + aLL*aLL);
    gmultipole += 0.125*(RAB + DA2 + DB2)/(aux*aux*aux);
    aux = sqrt((RAB + DA2 + DB2)*(RAB + DA2 + DB2) + (DA2 + DB2)*(DA2 + DB2) + aLL*aLL);
    gmultipole -= 0.125*(RAB + DA2 + DB2)/(aux*aux*aux);
    aux = sqrt((RAB - DA2 - DB2)*(RAB - DA2 - DB2) + (DA2 - DB2)*(DA2 - DB2) + aLL*aLL);
    gmultipole += 0.125*(RAB - DA2 - DB2)/(aux*aux*aux);
    aux = sqrt((RAB - DA2 - DB2)*(RAB - DA2 - DB2) + (DA2 + DB2)*(DA2 + DB2) + aLL*aLL);
    gmultipole -= 0.125*(RAB - DA2 - DB2)/(aux*aux*aux);
    aux = sqrt((RAB - DA2 + DB2)*(RAB - DA2 + DB2) + (DA2 - DB2)*(DA2 - DB2) + aLL*aLL);
    gmultipole -= 0.125*(RAB - DA2 + DB2)/(aux*aux*aux);
    aux = sqrt((RAB - DA2 + DB2)*(RAB - DA2 + DB2) + (DA2 + DB2)*(DA2 + DB2) + aLL*aLL);
    gmultipole += 0.125*(RAB - DA2 + DB2)/(aux*aux*aux);
  }
  else if ((l1 == 2)&&(m1 == 3)&&(l2 == 2)&&(m2 == 3)) {
    aux = sqrt(RAB*RAB + 2.0*(DA2 - DB2)*(DA2 - DB2) + aLL*aLL);
    gmultipole = -0.25*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + 2.0*(DA2 + DB2)*(DA2 + DB2) + aLL*aLL);
    gmultipole -= 0.25*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + 2.0*DA2*DA2 + 2.0*DB2*DB2 + aLL*aLL);
    gmultipole += 0.5*RAB/(aux*aux*aux);
  }
  return gmultipole;
}

double se_multipole_dR2(int l1, int m1, int l2, int m2, double RAB, double aLL, std::vector<double> D) {
  //auxiliary function in the calculation of the second-derivatives of ERIs for Dewar's electron-repulsion integrals 
  //with respect to internuclear distance
  //fuction giving back (d2/dR2)[l1m1,l2m2], where limi (i=1,2) is given by
  //00 -> q
  //10 -> {mu}_{z}
  //11 -> {mu}_{x}
  //20 -> {Q}_{zz}
  //21 -> {Q}_{xz}
  //22 -> {Q}_{xx}
  //2-2 -> {Q}_{yy}
  //23 -> {Q}_{xy}
  double DA1 = D[0];
  double DA2 = D[1];
  double DB1 = D[2];
  double DB2 = D[3];
  double hmultipole = 0.0;
  double aux;
  if ((l1 == 0)&&(l2 == 0)) {
    aux = RAB*RAB + aLL*aLL;
    hmultipole = (3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
  } 
  else if ((l1 == 0)&&(l2 == 1)&&(m2 == 0)) {
    aux = (RAB + DB1)*(RAB + DB1) + aLL*aLL;
    hmultipole = -0.5*(3.0*(RAB + DB1)*(RAB + DB1) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DB1)*(RAB - DB1) + aLL*aLL;
    hmultipole += 0.5*(3.0*(RAB - DB1)*(RAB - DB1) - aux)/(aux*aux*sqrt(aux));
  }
  else if ((l1 == 1)&&(m1 == 0)&&(l2 == 0)) {
    DA1 = D[2];
    DA2 = D[3];
    DB1 = D[0];
    DB2 = D[1];
    aux = (RAB + DB1)*(RAB + DB1) + aLL*aLL;
    hmultipole = 0.5*(3.0*(RAB + DB1)*(RAB + DB1) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DB1)*(RAB - DB1) + aLL*aLL;
    hmultipole -= 0.5*(3.0*(RAB - DB1)*(RAB - DB1) - aux)/(aux*aux*sqrt(aux));
  }
  else if ((l1 == 0)&&(l2 == 2)&&(abs(m2) == 2)) {
    aux = RAB*RAB + 4.0*DB2*DB2 + aLL*aLL;
    hmultipole = 0.5*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + aLL*aLL;
    hmultipole -= 0.5*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
  }
  else if ((l1 == 2)&&(abs(m1) == 2)&&(l2 == 0)) {
    DA1 = D[2];
    DA2 = D[3];
    DB1 = D[0];
    DB2 = D[1];
    aux = RAB*RAB + 4.0*DB2*DB2 + aLL*aLL;
    hmultipole = 0.5*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + aLL*aLL;
    hmultipole -= 0.5*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
  }
  else if ((l1 == 0)&&(l2 == 2)&&(m2 == 0)) {
    aux = (RAB + 2.0*DB2)*(RAB + 2.0*DB2) + aLL*aLL;
    hmultipole = 0.25*(3.0*(RAB + 2.0*DB2)*(RAB + 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + aLL*aLL;
    hmultipole -= 0.5*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - 2.0*DB2)*(RAB - 2.0*DB2) + aLL*aLL;
    hmultipole += 0.25*(3.0*(RAB - 2.0*DB2)*(RAB - 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
  }
  else if ((l1 == 2)&&(m1 == 0)&&(l2 == 0)) {
    DA1 = D[2];
    DA2 = D[3];
    DB1 = D[0];
    DB2 = D[1];
    aux = (RAB + 2.0*DB2)*(RAB + 2.0*DB2) + aLL*aLL;
    hmultipole = 0.25*(3.0*(RAB + 2.0*DB2)*(RAB + 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + aLL*aLL;
    hmultipole -= 0.5*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - 2.0*DB2)*(RAB - 2.0*DB2) + aLL*aLL;
    hmultipole += 0.25*(3.0*(RAB - 2.0*DB2)*(RAB - 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
  }
  else if ((l1 == 1)&&(abs(m1) == 1)&&(l2 == 1)&&(abs(m2) == 1)) {
    aux = RAB*RAB + (DA1 - DB1)*(DA1 - DB1) + aLL*aLL;
    hmultipole = 0.5*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + (DA1 + DB1)*(DA1 + DB1) + aLL*aLL;
    hmultipole -= 0.5*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
  }
  else if ((l1 == 1)&&(m1 == 0)&&(l2 == 1)&&(m2 == 0)) {
    aux = (RAB + DA1 - DB1)*(RAB + DA1 - DB1) + aLL*aLL;
    hmultipole = 0.25*(3.0*(RAB + DA1 - DB1)*(RAB + DA1 - DB1) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DA1 + DB1)*(RAB + DA1 + DB1) + aLL*aLL;
    hmultipole -= 0.25*(3.0*(RAB + DA1 + DB1)*(RAB + DA1 + DB1) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA1 - DB1)*(RAB - DA1 - DB1) + aLL*aLL;
    hmultipole -= 0.25*(3.0*(RAB - DA1 - DB1)*(RAB - DA1 - DB1) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA1 + DB1)*(RAB - DA1 + DB1) + aLL*aLL;
    hmultipole += 0.25*(3.0*(RAB - DA1 + DB1)*(RAB - DA1 + DB1) - aux)/(aux*aux*sqrt(aux));
  }
  else if ((l1 == 1)&&(abs(m1) == 1)&&(l2 == 2)&&(m2 == 1)) {
    aux = (RAB - DB2)*(RAB - DB2) + (DA1 - DB2)*(DA1 - DB2) + aLL*aLL;
    hmultipole = 0.25*(3.0*(RAB - DB2)*(RAB - DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DB2)*(RAB - DB2) + (DA1 + DB2)*(DA1 + DB2) + aLL*aLL;
    hmultipole -= 0.25*(3.0*(RAB - DB2)*(RAB - DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DB2)*(RAB + DB2) + (DA1 - DB2)*(DA1 - DB2) + aLL*aLL;
    hmultipole -= 0.25*(3.0*(RAB + DB2)*(RAB + DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DB2)*(RAB + DB2) + (DA1 + DB2)*(DA1 + DB2) + aLL*aLL;
    hmultipole += 0.25*(3.0*(RAB + DB2)*(RAB + DB2) - aux)/(aux*aux*sqrt(aux));
  }
  else if ((l1 == 2)&&(m1 == 1)&&(l2 == 1)&&(abs(m2) == 1)) {
    DA1 = D[2];
    DA2 = D[3];
    DB1 = D[0];
    DB2 = D[1];
    aux = (RAB - DB2)*(RAB - DB2) + (DA1 - DB2)*(DA1 - DB2) + aLL*aLL;
    hmultipole = -0.25*(3.0*(RAB - DB2)*(RAB - DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DB2)*(RAB - DB2) + (DA1 + DB2)*(DA1 + DB2) + aLL*aLL;
    hmultipole += 0.25*(3.0*(RAB - DB2)*(RAB - DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DB2)*(RAB + DB2) + (DA1 - DB2)*(DA1 - DB2) + aLL*aLL;
    hmultipole += 0.25*(3.0*(RAB + DB2)*(RAB + DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DB2)*(RAB + DB2) + (DA1 + DB2)*(DA1 + DB2) + aLL*aLL;
    hmultipole -= 0.25*(3.0*(RAB + DB2)*(RAB + DB2) - aux)/(aux*aux*sqrt(aux));
  }
  else if ((l1 == 1)&&(m1 == 0)&&(l2 == 2)&&(abs(m2) == 2)) {
    aux = (RAB + DA1)*(RAB + DA1) + 4.0*DB2*DB2 + aLL*aLL;
    hmultipole = 0.25*(3.0*(RAB + DA1)*(RAB + DA1) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA1)*(RAB - DA1) + 4.0*DB2*DB2 + aLL*aLL;
    hmultipole -= 0.25*(3.0*(RAB - DA1)*(RAB - DA1) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DA1)*(RAB + DA1) + aLL*aLL;
    hmultipole -= 0.25*(3.0*(RAB + DA1)*(RAB + DA1) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA1)*(RAB - DA1) + aLL*aLL;
    hmultipole += 0.25*(3.0*(RAB - DA1)*(RAB - DA1) - aux)/(aux*aux*sqrt(aux));
  }
  else if ((l1 == 2)&&(abs(m1) == 2)&&(l2 == 1)&&(m2 == 0)) {
    DA1 = D[2];
    DA2 = D[3];
    DB1 = D[0];
    DB2 = D[1];
    aux = (RAB + DA1)*(RAB + DA1) + 4.0*DB2*DB2 + aLL*aLL;
    hmultipole = -0.25*(3.0*(RAB + DA1)*(RAB + DA1) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA1)*(RAB - DA1) + 4.0*DB2*DB2 + aLL*aLL;
    hmultipole += 0.25*(3.0*(RAB - DA1)*(RAB - DA1) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DA1)*(RAB + DA1) + aLL*aLL;
    hmultipole += 0.25*(3.0*(RAB + DA1)*(RAB + DA1) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA1)*(RAB - DA1) + aLL*aLL;
    hmultipole -= 0.25*(3.0*(RAB - DA1)*(RAB - DA1) - aux)/(aux*aux*sqrt(aux));
  }
  else if ((l1 == 1)&&(m1 == 0)&&(l2 == 2)&&(m2 == 0)) {
    aux = (RAB + DA1 - 2.0*DB2)*(RAB + DA1 - 2.0*DB2) + aLL*aLL;
    hmultipole = 0.125*(3.0*(RAB + DA1 - 2.0*DB2)*(RAB + DA1 - 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA1 - 2.0*DB2)*(RAB - DA1 - 2.0*DB2) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB - DA1 - 2.0*DB2)*(RAB - DA1 - 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DA1 + 2.0*DB2)*(RAB + DA1 + 2.0*DB2) + aLL*aLL;
    hmultipole += 0.125*(3.0*(RAB + DA1 + 2.0*DB2)*(RAB + DA1 + 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA1 + 2.0*DB2)*(RAB - DA1 + 2.0*DB2) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB - DA1 + 2.0*DB2)*(RAB - DA1 + 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DA1)*(RAB + DA1) + aLL*aLL;
    hmultipole -= 0.25*(3.0*(RAB + DA1)*(RAB + DA1) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA1)*(RAB - DA1) + aLL*aLL;
    hmultipole += 0.25*(3.0*(RAB - DA1)*(RAB - DA1) - aux)/(aux*aux*sqrt(aux));
  }
  else if ((l1 == 2)&&(m1 == 0)&&(l2 == 1)&&(m2 == 0)) {
    DA1 = D[2];
    DA2 = D[3];
    DB1 = D[0];
    DB2 = D[1];
    aux = (RAB + DA1 - 2.0*DB2)*(RAB + DA1 - 2.0*DB2) + aLL*aLL;
    hmultipole = -0.125*(3.0*(RAB + DA1 - 2.0*DB2)*(RAB + DA1 - 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA1 - 2.0*DB2)*(RAB - DA1 - 2.0*DB2) + aLL*aLL;
    hmultipole += 0.125*(3.0*(RAB - DA1 - 2.0*DB2)*(RAB - DA1 - 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DA1 + 2.0*DB2)*(RAB + DA1 + 2.0*DB2) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB + DA1 + 2.0*DB2)*(RAB + DA1 + 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA1 + 2.0*DB2)*(RAB - DA1 + 2.0*DB2) + aLL*aLL;
    hmultipole += 0.125*(3.0*(RAB - DA1 + 2.0*DB2)*(RAB - DA1 + 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DA1)*(RAB + DA1) + aLL*aLL;
    hmultipole += 0.25*(3.0*(RAB + DA1)*(RAB + DA1) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA1)*(RAB - DA1) + aLL*aLL;
    hmultipole -= 0.25*(3.0*(RAB - DA1)*(RAB - DA1) - aux)/(aux*aux*sqrt(aux));
  }
  else if ((l1 == 2)&&(abs(m1) == 2)&&(l2 == 2)&&(abs(m2) == 2)&&(m1 == m2)) {
    aux = RAB*RAB + 4.0*(DA2 - DB2)*(DA2 - DB2) + aLL*aLL;
    hmultipole = 0.125*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + 4.0*(DA2 + DB2)*(DA2 + DB2) + aLL*aLL;
    hmultipole += 0.125*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + 4.0*DA2*DA2 + aLL*aLL;
    hmultipole -= 0.25*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + 4.0*DB2*DB2 + aLL*aLL;
    hmultipole -= 0.25*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + aLL*aLL;
    hmultipole += 0.25*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
  }
  else if ((l1 == 2)&&(m1 == 2)&&(l2 == 2)&&(m2 == -2)) {
    aux = RAB*RAB + 4.0*DA2*DA2 + 4.0*DB2*DB2 + aLL*aLL;
    hmultipole = 0.25*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + 4.0*DA2*DA2 + aLL*aLL;
    hmultipole -= 0.25*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + 4.0*DB2*DB2 + aLL*aLL;
    hmultipole -= 0.25*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + aLL*aLL;
    hmultipole += 0.25*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
  }
  else if ((l1 == 2)&&(m1 == -2)&&(l2 == 2)&&(m2 == 2)) {
    aux = RAB*RAB + 4.0*DA2*DA2 + 4.0*DB2*DB2 + aLL*aLL;
    hmultipole = 0.25*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + 4.0*DA2*DA2 + aLL*aLL;
    hmultipole -= 0.25*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + 4.0*DB2*DB2 + aLL*aLL;
    hmultipole -= 0.25*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + aLL*aLL;
    hmultipole += 0.25*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
  }
  else if ((l1 == 2)&&(abs(m1) == 2)&&(l2 == 2)&&(m2 == 0)) {
    aux = (RAB - 2.0*DB2)*(RAB - 2.0*DB2) + 4.0*DA2*DA2 + aLL*aLL;
    hmultipole = 0.125*(3.0*(RAB - 2.0*DB2)*(RAB - 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + 2.0*DB2)*(RAB + 2.0*DB2) + 4.0*DA2*DA2 + aLL*aLL;
    hmultipole += 0.125*(3.0*(RAB + 2.0*DB2)*(RAB + 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - 2.0*DB2)*(RAB - 2.0*DB2) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB - 2.0*DB2)*(RAB - 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + 2.0*DB2)*(RAB + 2.0*DB2) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB + 2.0*DB2)*(RAB + 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + 4.0*DA2*DA2 + aLL*aLL;
    hmultipole -= 0.25*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + aLL*aLL;
    hmultipole += 0.25*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
  }
  else if ((l1 == 2)&&(m1 == 0)&&(l2 == 2)&&(abs(m2) == 2)) {
    DA1 = D[2];
    DA2 = D[3];
    DB1 = D[0];
    DB2 = D[1];
    aux = (RAB - 2.0*DB2)*(RAB - 2.0*DB2) + 4.0*DA2*DA2 + aLL*aLL;
    hmultipole = 0.125*(3.0*(RAB - 2.0*DB2)*(RAB - 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + 2.0*DB2)*(RAB + 2.0*DB2) + 4.0*DA2*DA2 + aLL*aLL;
    hmultipole += 0.125*(3.0*(RAB + 2.0*DB2)*(RAB + 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - 2.0*DB2)*(RAB - 2.0*DB2) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB - 2.0*DB2)*(RAB - 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + 2.0*DB2)*(RAB + 2.0*DB2) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB + 2.0*DB2)*(RAB + 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + 4.0*DA2*DA2 + aLL*aLL;
    hmultipole -= 0.25*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + aLL*aLL;
    hmultipole += 0.25*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
  }
  else if ((l1 == 2)&&(m1 == 0)&&(l2 == 2)&&(m2 == 0)) {
    aux = (RAB + 2.0*DA2 - 2.0*DB2)*(RAB + 2.0*DA2 - 2.0*DB2) + aLL*aLL;
    hmultipole = 0.0625*(3.0*(RAB + 2.0*DA2 - 2.0*DB2)*(RAB + 2.0*DA2 - 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + 2.0*DA2 + 2.0*DB2)*(RAB + 2.0*DA2 + 2.0*DB2) + aLL*aLL;
    hmultipole += 0.0625*(3.0*(RAB + 2.0*DA2 + 2.0*DB2)*(RAB + 2.0*DA2 + 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - 2.0*DA2 - 2.0*DB2)*(RAB - 2.0*DA2 - 2.0*DB2) + aLL*aLL;
    hmultipole += 0.0625*(3.0*(RAB - 2.0*DA2 - 2.0*DB2)*(RAB - 2.0*DA2 - 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - 2.0*DA2 + 2.0*DB2)*(RAB - 2.0*DA2 + 2.0*DB2) + aLL*aLL;
    hmultipole += 0.0625*(3.0*(RAB - 2.0*DA2 + 2.0*DB2)*(RAB - 2.0*DA2 + 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + 2.0*DA2)*(RAB + 2.0*DA2) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB + 2.0*DA2)*(RAB + 2.0*DA2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - 2.0*DA2)*(RAB - 2.0*DA2) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB - 2.0*DA2)*(RAB - 2.0*DA2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + 2.0*DB2)*(RAB + 2.0*DB2) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB + 2.0*DB2)*(RAB + 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - 2.0*DB2)*(RAB - 2.0*DB2) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB - 2.0*DB2)*(RAB - 2.0*DB2) - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + aLL*aLL;
    hmultipole += 0.25*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
  }
  else if ((l1 == 2)&&(m1 == 1)&&(l2 == 2)&&(m2 == 1)) {
    aux = (RAB + DA2 - DB2)*(RAB + DA2 - DB2) + (DA2 - DB2)*(DA2 - DB2) + aLL*aLL;
    hmultipole = 0.125*(3.0*(RAB + DA2 - DB2)*(RAB + DA2 - DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DA2 - DB2)*(RAB + DA2 - DB2) + (DA2 + DB2)*(DA2 + DB2) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB + DA2 - DB2)*(RAB + DA2 - DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DA2 + DB2)*(RAB + DA2 + DB2) + (DA2 - DB2)*(DA2 - DB2) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB + DA2 + DB2)*(RAB + DA2 + DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DA2 + DB2)*(RAB + DA2 + DB2) + (DA2 + DB2)*(DA2 + DB2) + aLL*aLL;
    hmultipole += 0.125*(3.0*(RAB + DA2 + DB2)*(RAB + DA2 + DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA2 - DB2)*(RAB - DA2 - DB2) + (DA2 - DB2)*(DA2 - DB2) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB - DA2 - DB2)*(RAB - DA2 - DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA2 - DB2)*(RAB - DA2 - DB2) + (DA2 + DB2)*(DA2 + DB2) + aLL*aLL;
    hmultipole += 0.125*(3.0*(RAB - DA2 - DB2)*(RAB - DA2 - DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA2 + DB2)*(RAB - DA2 + DB2) + (DA2 - DB2)*(DA2 - DB2) + aLL*aLL;
    hmultipole += 0.125*(3.0*(RAB - DA2 + DB2)*(RAB - DA2 + DB2) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA2 + DB2)*(RAB - DA2 + DB2) + (DA2 + DB2)*(DA2 + DB2) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB - DA2 + DB2)*(RAB - DA2 + DB2) - aux)/(aux*aux*sqrt(aux));
  }
  else if ((l1 == 2)&&(m1 == 3)&&(l2 == 2)&&(m2 == 3)) {
    aux = RAB*RAB + 2.0*(DA2 - DB2)*(DA2 - DB2) + aLL*aLL;
    hmultipole = 0.25*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + 2.0*(DA2 + DB2)*(DA2 + DB2) + aLL*aLL;
    hmultipole += 0.25*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + 2.0*DA2*DA2 + 2.0*DB2*DB2 + aLL*aLL;
    hmultipole -= 0.5*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
  }
  return hmultipole;
}

double d_se_multipole(int l1, int m1, int l2, int m2, double RAB, double aLL, std::vector<double> D) {
  //ERI auxiliary function in the calculation of Thiel and Voityuk's extension to electron-repulsion integrals involving d orbitals
  //Thiel, Voityuk Theor. Chim. Acta (Berl.), 81, 391, 1992
  //function giving back [l1m1,l2m2], where limi (i=1,2) is given by
  //00 -> q
  //10 -> {mu}_{z}
  //11 -> {mu}_{x}
  //20 -> {Q}_{zz}
  //21 -> {Q}_{xz}
  //22 -> {Q}_{xx}
  double DA = D[0];
  double DB = D[1];
  double multipole = 0.0;
  //[q,q]
  if ((l1 == 0)&&(l2 == 0)) {multipole = 1.0/sqrt(RAB*RAB + aLL*aLL);} 
  //[q,muz]
  else if ((l1 == 0)&&(l2 == 1)&&(m2 == 0)) {multipole = -0.5*(1.0/sqrt((RAB + DB)*(RAB + DB) + aLL*aLL) - 1.0/sqrt((RAB - DB)*(RAB - DB) + aLL*aLL));}
  //[muz,q]
  else if ((l1 == 1)&&(m1 == 0)&&(l2 == 0)) {
    DB = D[0];
    DA = D[1];
    multipole = 0.5*(1.0/sqrt((RAB + DB)*(RAB + DB) + aLL*aLL) - 1.0/sqrt((RAB - DB)*(RAB - DB) + aLL*aLL));
  }
  //[q,Qzz]
  else if ((l1 == 0)&&(l2 == 2)&&(m2 == 0)) {
    multipole = 0.25*(1.0/sqrt((RAB + DB)*(RAB + DB) + aLL*aLL) - 2.0/sqrt(RAB*RAB + DB*DB + aLL*aLL) + 1.0/sqrt((RAB - DB)*(RAB - DB) + aLL*aLL));
  }
  //[Qzz,q]
  else if ((l1 == 2)&&(m1 == 0)&&(l2 == 0)) {
    DB = D[0];
    DA = D[1];
    multipole = 0.25*(1.0/sqrt((RAB + DB)*(RAB + DB) + aLL*aLL) - 2.0/sqrt(RAB*RAB + DB*DB + aLL*aLL) + 1.0/sqrt((RAB - DB)*(RAB - DB) + aLL*aLL));
  }
  //[mux,mux]
  else if ((l1 == 1)&&(abs(m1) == 1)&&(l2 == 1)&&(abs(m2) == 1)) {
    multipole = 0.5*(1.0/sqrt(RAB*RAB + (DA - DB)*(DA - DB) + aLL*aLL) - 1.0/sqrt(RAB*RAB + (DA + DB)*(DA + DB) + aLL*aLL));
  }
  //[muz,muz]
  else if ((l1 == 1)&&(m1 == 0)&&(l2 == 1)&&(m2 == 0)) {
    multipole = 0.25*(1.0/sqrt((RAB + DA - DB)*(RAB + DA - DB) + aLL*aLL) - 1.0/sqrt((RAB + DA + DB)*(RAB + DA + DB) + aLL*aLL) - 1.0/sqrt((RAB - DA - DB)*(RAB - DA - DB) + aLL*aLL) + 1.0/sqrt((RAB - DA + DB)*(RAB - DA + DB) + aLL*aLL));
  }
  //[mux,Qxz]
  else if ((l1 == 1)&&(abs(m1) == 1)&&(l2 == 2)&&(abs(m2) == 1)&&(m1 == m2)) {
    double sq2 = 1.0/sqrt(2.0);
    multipole = -0.25*(-1.0/sqrt((RAB - sq2*DB)*(RAB - sq2*DB) + (DA - sq2*DB)*(DA - sq2*DB) + aLL*aLL) + 1.0/sqrt((RAB - sq2*DB)*(RAB - sq2*DB) + (DA + sq2*DB)*(DA + sq2*DB) + aLL*aLL) + 1.0/sqrt((RAB + sq2*DB)*(RAB + sq2*DB) + (DA - sq2*DB)*(DA - sq2*DB) + aLL*aLL) - 1.0/sqrt((RAB + sq2*DB)*(RAB + sq2*DB) + (DA + sq2*DB)*(DA + sq2*DB) + aLL*aLL));
  }
  //[Qxz,mux]
  else if ((l1 == 2)&&(abs(m1) == 1)&&(l2 == 1)&&(abs(m2) == 1)&&(m1 == m2)) {
    double sq2 = 1.0/sqrt(2.0);
    DB = D[0];
    DA = D[1];
    multipole = 0.25*(-1.0/sqrt((RAB - sq2*DB)*(RAB - sq2*DB) + (DA - sq2*DB)*(DA - sq2*DB) + aLL*aLL) + 1.0/sqrt((RAB - sq2*DB)*(RAB - sq2*DB) + (DA + sq2*DB)*(DA + sq2*DB) + aLL*aLL) + 1.0/sqrt((RAB + sq2*DB)*(RAB + sq2*DB) + (DA - sq2*DB)*(DA - sq2*DB) + aLL*aLL) - 1.0/sqrt((RAB + sq2*DB)*(RAB + sq2*DB) + (DA + sq2*DB)*(DA + sq2*DB) + aLL*aLL));
  }
  //[muz,Qzz]
  else if ((l1 == 1)&&(m1 == 0)&&(l2 == 2)&&(m2 == 0)) {
    multipole = -0.125*(-1.0/sqrt((RAB + DA - DB)*(RAB + DA - DB) + aLL*aLL) + 1.0/sqrt((RAB - DA - DB)*(RAB - DA - DB) + aLL*aLL) - 1.0/sqrt((RAB + DA + DB)*(RAB + DA + DB) + aLL*aLL) + 1.0/sqrt((RAB - DA + DB)*(RAB - DA + DB) + aLL*aLL) + 2.0/sqrt((RAB + DA)*(RAB + DA) + DB*DB + aLL*aLL) - 2.0/sqrt((RAB - DA)*(RAB - DA) + DB*DB + aLL*aLL));
  }
  //[Qzz,muz]
  else if ((l1 == 2)&&(m1 == 0)&&(l2 == 1)&&(m2 == 0)) {
    DB = D[0];
    DA = D[1];
    multipole = 0.125*(-1.0/sqrt((RAB + DA - DB)*(RAB + DA - DB) + aLL*aLL) + 1.0/sqrt((RAB - DA - DB)*(RAB - DA - DB) + aLL*aLL) - 1.0/sqrt((RAB + DA + DB)*(RAB + DA + DB) + aLL*aLL) + 1.0/sqrt((RAB - DA + DB)*(RAB - DA + DB) + aLL*aLL) + 2.0/sqrt((RAB + DA)*(RAB + DA) + DB*DB + aLL*aLL) - 2.0/sqrt((RAB - DA)*(RAB - DA) + DB*DB + aLL*aLL));
  }
  //[Qxy,Qxy]
  else if ((l1 == 2)&&(abs(m1) == 2)&&(l2 == 2)&&(abs(m2) == 2)&&(m1 == m2)) {
    multipole = 0.25*(1.0/sqrt(RAB*RAB + (DA - DB)*(DA - DB) + aLL*aLL) + 1.0/sqrt(RAB*RAB + (DA + DB)*(DA + DB) + aLL*aLL) - 2.0/sqrt(RAB*RAB + DA*DA + DB*DB + aLL*aLL));
  }
  //[Qzz,Qzz]
  else if ((l1 == 2)&&(m1 == 0)&&(l2 == 2)&&(m2 == 0)) {
    multipole = 0.0625*(1.0/sqrt((RAB + DA - DB)*(RAB + DA - DB) + aLL*aLL) + 1.0/sqrt((RAB + DA + DB)*(RAB + DA + DB) + aLL*aLL) + 1.0/sqrt((RAB - DA - DB)*(RAB - DA - DB) + aLL*aLL) + 1.0/sqrt((RAB - DA + DB)*(RAB - DA + DB) + aLL*aLL) - 2.0/sqrt((RAB + DA)*(RAB + DA) + DB*DB + aLL*aLL) - 2.0/sqrt((RAB - DA)*(RAB - DA) + DB*DB + aLL*aLL) - 2.0/sqrt((RAB + DB)*(RAB + DB) + DA*DA + aLL*aLL) - 2.0/sqrt((RAB - DB)*(RAB - DB) + DA*DA + aLL*aLL) + 2.0/sqrt(RAB*RAB + DA*DA + DB*DB + aLL*aLL));
    multipole += 0.0625*(1.0/sqrt(RAB*RAB + (DA - DB)*(DA - DB) + aLL*aLL) + 1.0/sqrt(RAB*RAB + (DA + DB)*(DA + DB) + aLL*aLL));
  }
  //[Qxz,Qxz]
  else if ((l1 == 2)&&(abs(m1) == 1)&&(l2 == 2)&&(abs(m2) == 1)&&(m1 == m2)) {
    double sq2 = 1.0/sqrt(2.0);
    multipole = 0.125*(1.0/sqrt((RAB + sq2*DA - sq2*DB)*(RAB + sq2*DA - sq2*DB) + 0.5*(DA - DB)*(DA - DB) + aLL*aLL) - 1.0/sqrt((RAB + sq2*DA - sq2*DB)*(RAB + sq2*DA - sq2*DB) + 0.5*(DA + DB)*(DA + DB) + aLL*aLL) - 1.0/sqrt((RAB + sq2*DA + sq2*DB)*(RAB + sq2*DA + sq2*DB) + 0.5*(DA - DB)*(DA - DB) + aLL*aLL) + 1.0/sqrt((RAB + sq2*DA + sq2*DB)*(RAB + sq2*DA + sq2*DB) + 0.5*(DA + DB)*(DA + DB) + aLL*aLL) - 1.0/sqrt((RAB - sq2*DA - sq2*DB)*(RAB - sq2*DA - sq2*DB) + 0.5*(DA - DB)*(DA - DB) + aLL*aLL) + 1.0/sqrt((RAB - sq2*DA - sq2*DB)*(RAB - sq2*DA - sq2*DB) + 0.5*(DA + DB)*(DA + DB) + aLL*aLL) + 1.0/sqrt((RAB - sq2*DA + sq2*DB)*(RAB - sq2*DA + sq2*DB) + 0.5*(DA - DB)*(DA - DB) + aLL*aLL) - 1.0/sqrt((RAB - sq2*DA + sq2*DB)*(RAB - sq2*DA + sq2*DB) + 0.5*(DA + DB)*(DA + DB) + aLL*aLL));
  }
  return multipole;
}
double d_se_multipole_dR(int l1, int m1, int l2, int m2, double RAB, double aLL, std::vector<double> D) {
  //auxiliary function in the calculation of the derivative of ERIs of Thiel and Voityuk's extension to electron-repulsion integrals involving d orbitals
  //function giving back [l1m1,l2m2], where limi (i=1,2) is given by
  //00 -> q
  //10 -> {mu}_{z}
  //11 -> {mu}_{x}
  //20 -> {Q}_{zz}
  //21 -> {Q}_{xz}
  //22 -> {Q}_{xx}
  double DA = D[0];
  double DB = D[1];
  double gmultipole = 0.0;
  double aux;
  //[q,q]
  if ((l1 == 0)&&(l2 == 0)) {
    aux = sqrt(RAB*RAB + aLL*aLL);
    gmultipole = -RAB/(aux*aux*aux);
  }
  //[q,muz]
  else if ((l1 == 0)&&(l2 == 1)&&(m2 == 0)) {
    aux = sqrt((RAB + DB)*(RAB + DB) + aLL*aLL);
    gmultipole = 0.5*(RAB + DB)/(aux*aux*aux);
    aux = sqrt((RAB - DB)*(RAB - DB) + aLL*aLL);
    gmultipole -= 0.5*(RAB - DB)/(aux*aux*aux);
  }
  //[muz,q]
  else if ((l1 == 1)&&(m1 == 0)&&(l2 == 0)) {
    DB = D[0];
    DA = D[1];
    aux = sqrt((RAB + DB)*(RAB + DB) + aLL*aLL);
    gmultipole = -0.5*(RAB + DB)/(aux*aux*aux);
    aux = sqrt((RAB - DB)*(RAB - DB) + aLL*aLL);
    gmultipole += 0.5*(RAB - DB)/(aux*aux*aux);
  }
  //[q,Qzz]
  else if ((l1 == 0)&&(l2 == 2)&&(m2 == 0)) {
    aux = sqrt((RAB + DB)*(RAB + DB) + aLL*aLL);
    gmultipole -= 0.25*(RAB + DB)/(aux*aux*aux);
    aux = sqrt(RAB*RAB + DB*DB + aLL*aLL);
    gmultipole += 0.5*RAB/(aux*aux*aux);
    aux = sqrt((RAB - DB)*(RAB - DB) + aLL*aLL);
    gmultipole -= 0.25*(RAB - DB)/(aux*aux*aux);
  }
  //[Qzz,q]
  else if ((l1 == 2)&&(m1 == 0)&&(l2 == 0)) {
    DB = D[0];
    DA = D[1];
    aux = sqrt((RAB + DB)*(RAB + DB) + aLL*aLL);
    gmultipole -= 0.25*(RAB + DB)/(aux*aux*aux);
    aux = sqrt(RAB*RAB + DB*DB + aLL*aLL);
    gmultipole += 0.5*RAB/(aux*aux*aux);
    aux = sqrt((RAB - DB)*(RAB - DB) + aLL*aLL);
    gmultipole -= 0.25*(RAB - DB)/(aux*aux*aux);
  }
  //[mux,mux]
  else if ((l1 == 1)&&(abs(m1) == 1)&&(l2 == 1)&&(abs(m2) == 1)) {
    aux = sqrt(RAB*RAB + (DA - DB)*(DA - DB) + aLL*aLL);
    gmultipole -= 0.5*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + (DA + DB)*(DA + DB) + aLL*aLL);
    gmultipole += 0.5*RAB/(aux*aux*aux);
  }
  //[muz,muz]
  else if ((l1 == 1)&&(m1 == 0)&&(l2 == 1)&&(m2 == 0)) {
    aux = sqrt((RAB + DA - DB)*(RAB + DA - DB) + aLL*aLL);
    gmultipole -= 0.25*(RAB + DA - DB)/(aux*aux*aux);
    aux = sqrt((RAB + DA + DB)*(RAB + DA + DB) + aLL*aLL);
    gmultipole += 0.25*(RAB + DA + DB)/(aux*aux*aux);
    aux = sqrt((RAB - DA - DB)*(RAB - DA - DB) + aLL*aLL);
    gmultipole += 0.25*(RAB - DA - DB)/(aux*aux*aux);
    aux = sqrt((RAB - DA + DB)*(RAB - DA + DB) + aLL*aLL);
    gmultipole -= 0.25*(RAB - DA + DB)/(aux*aux*aux);
  }
  //[mux,Qxz]
  else if ((l1 == 1)&&(abs(m1) == 1)&&(l2 == 2)&&(abs(m2) == 1)&&(m1 == m2)) {
    double sq2 = 1.0/sqrt(2.0);
    aux = sqrt((RAB - sq2*DB)*(RAB - sq2*DB) + (DA - sq2*DB)*(DA - sq2*DB) + aLL*aLL);
    gmultipole -= 0.25*(RAB - sq2*DB)/(aux*aux*aux);
    aux = sqrt((RAB - sq2*DB)*(RAB - sq2*DB) + (DA + sq2*DB)*(DA + sq2*DB) + aLL*aLL);
    gmultipole += 0.25*(RAB - sq2*DB)/(aux*aux*aux);
    aux = sqrt((RAB + sq2*DB)*(RAB + sq2*DB) + (DA - sq2*DB)*(DA - sq2*DB) + aLL*aLL);
    gmultipole += 0.25*(RAB + sq2*DB)/(aux*aux*aux);
    aux = sqrt((RAB + sq2*DB)*(RAB + sq2*DB) + (DA + sq2*DB)*(DA + sq2*DB) + aLL*aLL);
    gmultipole -= 0.25*(RAB + sq2*DB)/(aux*aux*aux);
  }
  //[Qxz,mux]
  else if ((l1 == 2)&&(abs(m1) == 1)&&(l2 == 1)&&(abs(m2) == 1)&&(m1 == m2)) {
    double sq2 = 1.0/sqrt(2.0);
    DB = D[0];
    DA = D[1];
    aux = sqrt((RAB - sq2*DB)*(RAB - sq2*DB) + (DA - sq2*DB)*(DA - sq2*DB) + aLL*aLL);
    gmultipole += 0.25*(RAB - sq2*DB)/(aux*aux*aux);
    aux = sqrt((RAB - sq2*DB)*(RAB - sq2*DB) + (DA + sq2*DB)*(DA + sq2*DB) + aLL*aLL);
    gmultipole -= 0.25*(RAB - sq2*DB)/(aux*aux*aux);
    aux = sqrt((RAB + sq2*DB)*(RAB + sq2*DB) + (DA - sq2*DB)*(DA - sq2*DB) + aLL*aLL);
    gmultipole -= 0.25*(RAB + sq2*DB)/(aux*aux*aux);
    aux = sqrt((RAB + sq2*DB)*(RAB + sq2*DB) + (DA + sq2*DB)*(DA + sq2*DB) + aLL*aLL);
    gmultipole += 0.25*(RAB + sq2*DB)/(aux*aux*aux);
  }
  //[muz,Qzz]
  else if ((l1 == 1)&&(m1 == 0)&&(l2 == 2)&&(m2 == 0)) {
    aux = sqrt((RAB + DA - DB)*(RAB + DA - DB) + aLL*aLL);
    gmultipole -= 0.125*(RAB + DA - DB)/(aux*aux*aux);
    aux = sqrt((RAB - DA - DB)*(RAB - DA - DB) + aLL*aLL);
    gmultipole += 0.125*(RAB - DA - DB)/(aux*aux*aux);
    aux = sqrt((RAB + DA + DB)*(RAB + DA + DB) + aLL*aLL);
    gmultipole -= 0.125*(RAB + DA + DB)/(aux*aux*aux);
    aux = sqrt((RAB - DA + DB)*(RAB - DA + DB) + aLL*aLL);
    gmultipole += 0.125*(RAB - DA + DB)/(aux*aux*aux);
    aux = sqrt((RAB + DA)*(RAB + DA) + DB*DB + aLL*aLL);
    gmultipole += 0.25*(RAB + DA)/(aux*aux*aux);
    aux = sqrt((RAB - DA)*(RAB - DA) + DB*DB + aLL*aLL);
    gmultipole -= 0.25*(RAB - DA)/(aux*aux*aux);
  }
  //[Qzz,muz]
  else if ((l1 == 2)&&(m1 == 0)&&(l2 == 1)&&(m2 == 0)) {
    DB = D[0];
    DA = D[1];
    aux = sqrt((RAB + DA - DB)*(RAB + DA - DB) + aLL*aLL);
    gmultipole += 0.125*(RAB + DA - DB)/(aux*aux*aux);
    aux = sqrt((RAB - DA - DB)*(RAB - DA - DB) + aLL*aLL);
    gmultipole -= 0.125*(RAB - DA - DB)/(aux*aux*aux);
    aux = sqrt((RAB + DA + DB)*(RAB + DA + DB) + aLL*aLL);
    gmultipole += 0.125*(RAB + DA + DB)/(aux*aux*aux);
    aux = sqrt((RAB - DA + DB)*(RAB - DA + DB) + aLL*aLL);
    gmultipole -= 0.125*(RAB - DA + DB)/(aux*aux*aux);
    aux = sqrt((RAB + DA)*(RAB + DA) + DB*DB + aLL*aLL);
    gmultipole -= 0.25*(RAB + DA)/(aux*aux*aux);
    aux = sqrt((RAB - DA)*(RAB - DA) + DB*DB + aLL*aLL);
    gmultipole += 0.25*(RAB - DA)/(aux*aux*aux);
  }
  //[Qxy,Qxy]
  else if ((l1 == 2)&&(abs(m1) == 2)&&(l2 == 2)&&(abs(m2) == 2)&&(m1 == m2)) {
    aux = sqrt(RAB*RAB + (DA - DB)*(DA - DB) + aLL*aLL);
    gmultipole -= 0.25*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + (DA + DB)*(DA + DB) + aLL*aLL);
    gmultipole -= 0.25*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + DA*DA + DB*DB + aLL*aLL);
    gmultipole += 0.5*RAB/(aux*aux*aux);
  }
  //[Qzz,Qzz]
  else if ((l1 == 2)&&(m1 == 0)&&(l2 == 2)&&(m2 == 0)) {
    aux = sqrt((RAB + DA - DB)*(RAB + DA - DB) + aLL*aLL);
    gmultipole -= 0.0625*(RAB + DA - DB)/(aux*aux*aux);
    aux = sqrt((RAB + DA + DB)*(RAB + DA + DB) + aLL*aLL);
    gmultipole -= 0.0625*(RAB + DA + DB)/(aux*aux*aux);
    aux = sqrt((RAB - DA - DB)*(RAB - DA - DB) + aLL*aLL);
    gmultipole -= 0.0625*(RAB - DA - DB)/(aux*aux*aux);
    aux = sqrt((RAB - DA + DB)*(RAB - DA + DB) + aLL*aLL);
    gmultipole -= 0.0625*(RAB - DA + DB)/(aux*aux*aux);
    aux = sqrt((RAB + DA)*(RAB + DA) + DB*DB + aLL*aLL);
    gmultipole += 0.125*(RAB + DA)/(aux*aux*aux);
    aux = sqrt((RAB - DA)*(RAB - DA) + DB*DB + aLL*aLL);
    gmultipole += 0.125*(RAB - DA)/(aux*aux*aux);
    aux = sqrt((RAB + DB)*(RAB + DB) + DA*DA + aLL*aLL);
    gmultipole += 0.125*(RAB + DB)/(aux*aux*aux);
    aux = sqrt((RAB - DB)*(RAB - DB) + DA*DA + aLL*aLL);
    gmultipole += 0.125*(RAB - DB)/(aux*aux*aux);
    aux = sqrt(RAB*RAB + DA*DA + DB*DB + aLL*aLL);
    gmultipole -= 0.125*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + (DA - DB)*(DA - DB) + aLL*aLL);
    gmultipole -= 0.0625*RAB/(aux*aux*aux);
    aux = sqrt(RAB*RAB + (DA + DB)*(DA + DB) + aLL*aLL);
    gmultipole -= 0.0625*RAB/(aux*aux*aux);
  }
  //[Qxz,Qxz]
  else if ((l1 == 2)&&(abs(m1) == 1)&&(l2 == 2)&&(abs(m2) == 1)&&(m1 == m2)) {
    double sq2 = 1.0/sqrt(2.0);
    aux = sqrt((RAB + sq2*DA - sq2*DB)*(RAB + sq2*DA - sq2*DB) + 0.5*(DA - DB)*(DA - DB) + aLL*aLL);
    gmultipole -= 0.125*(RAB + sq2*DA - sq2*DB)/(aux*aux*aux);
    aux = sqrt((RAB + sq2*DA - sq2*DB)*(RAB + sq2*DA - sq2*DB) + 0.5*(DA + DB)*(DA + DB) + aLL*aLL);
    gmultipole += 0.125*(RAB + sq2*DA - sq2*DB)/(aux*aux*aux);
    aux = sqrt((RAB + sq2*DA + sq2*DB)*(RAB + sq2*DA + sq2*DB) + 0.5*(DA - DB)*(DA - DB) + aLL*aLL);
    gmultipole += 0.125*(RAB + sq2*DA + sq2*DB)/(aux*aux*aux);
    aux = sqrt((RAB + sq2*DA + sq2*DB)*(RAB + sq2*DA + sq2*DB) + 0.5*(DA + DB)*(DA + DB) + aLL*aLL);
    gmultipole -= 0.125*(RAB + sq2*DA + sq2*DB)/(aux*aux*aux);
    aux = sqrt((RAB - sq2*DA - sq2*DB)*(RAB - sq2*DA - sq2*DB) + 0.5*(DA - DB)*(DA - DB) + aLL*aLL);
    gmultipole += 0.125*(RAB - sq2*DA - sq2*DB)/(aux*aux*aux);
    aux = sqrt((RAB - sq2*DA - sq2*DB)*(RAB - sq2*DA - sq2*DB) + 0.5*(DA + DB)*(DA + DB) + aLL*aLL);
    gmultipole -= 0.125*(RAB - sq2*DA - sq2*DB)/(aux*aux*aux);
    aux = sqrt((RAB - sq2*DA + sq2*DB)*(RAB - sq2*DA + sq2*DB) + 0.5*(DA - DB)*(DA - DB) + aLL*aLL);
    gmultipole -= 0.125*(RAB - sq2*DA + sq2*DB)/(aux*aux*aux);
    aux = sqrt((RAB - sq2*DA + sq2*DB)*(RAB - sq2*DA + sq2*DB) + 0.5*(DA + DB)*(DA + DB) + aLL*aLL);
    gmultipole += 0.125*(RAB - sq2*DA + sq2*DB)/(aux*aux*aux);
  }
  return gmultipole;
}
double d_se_multipole_dR2(int l1, int m1, int l2, int m2, double RAB, double aLL, std::vector<double> D) {
  //auxiliary function in the calculation of the second-derivatives of ERIs of Thiel and Voityuk's extension to electron-repulsion integrals involving d orbitals
  //function giving back [l1m1,l2m2], where limi (i=1,2) is given by
  //00 -> q
  //10 -> {mu}_{z}
  //11 -> {mu}_{x}
  //20 -> {Q}_{zz}
  //21 -> {Q}_{xz}
  //22 -> {Q}_{xx}
  double DA = D[0];
  double DB = D[1];
  double hmultipole = 0.0;
  double aux;
  //[q,q]
  if ((l1 == 0)&&(l2 == 0)) {
    aux = RAB*RAB + aLL*aLL;
    hmultipole = (3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
  }
  //[q,muz]
  else if ((l1 == 0)&&(l2 == 1)&&(m2 == 0)) {
    aux = (RAB + DB)*(RAB + DB) + aLL*aLL;
    hmultipole = -0.5*(3.0*(RAB + DB)*(RAB + DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DB)*(RAB - DB) + aLL*aLL;
    hmultipole += 0.5*(3.0*(RAB - DB)*(RAB - DB) - aux)/(aux*aux*sqrt(aux));
  }
  //[muz,q]
  else if ((l1 == 1)&&(m1 == 0)&&(l2 == 0)) {
    DB = D[0];
    DA = D[1];
    aux = (RAB + DB)*(RAB + DB) + aLL*aLL;
    hmultipole = 0.5*(3.0*(RAB + DB)*(RAB + DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DB)*(RAB - DB) + aLL*aLL;
    hmultipole -= 0.5*(3.0*(RAB - DB)*(RAB - DB) - aux)/(aux*aux*sqrt(aux));
  }
  //[q,Qzz]
  else if ((l1 == 0)&&(l2 == 2)&&(m2 == 0)) {
    aux = (RAB + DB)*(RAB + DB) + aLL*aLL;
    hmultipole = 0.25*(3.0*(RAB + DB)*(RAB + DB) - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + DB*DB + aLL*aLL;
    hmultipole -= 0.5*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DB)*(RAB - DB) + aLL*aLL;
    hmultipole += 0.25*(3.0*(RAB - DB)*(RAB - DB) - aux)/(aux*aux*sqrt(aux));
  }
  //[Qzz,q]
  else if ((l1 == 2)&&(m1 == 0)&&(l2 == 0)) {
    DB = D[0];
    DA = D[1];
    aux = (RAB + DB)*(RAB + DB) + aLL*aLL;
    hmultipole = 0.25*(3.0*(RAB + DB)*(RAB + DB) - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + DB*DB + aLL*aLL;
    hmultipole -= 0.5*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DB)*(RAB - DB) + aLL*aLL;
    hmultipole += 0.25*(3.0*(RAB - DB)*(RAB - DB) - aux)/(aux*aux*sqrt(aux));
  }
  //[mux,mux]
  else if ((l1 == 1)&&(abs(m1) == 1)&&(l2 == 1)&&(abs(m2) == 1)) {
    aux = RAB*RAB + (DA - DB)*(DA - DB) + aLL*aLL;
    hmultipole = 0.5*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + (DA + DB)*(DA + DB) + aLL*aLL;
    hmultipole -= 0.5*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
  }
  //[muz,muz]
  else if ((l1 == 1)&&(m1 == 0)&&(l2 == 1)&&(m2 == 0)) {
    aux = (RAB + DA - DB)*(RAB + DA - DB) + aLL*aLL;
    hmultipole = 0.25*(3.0*(RAB + DA - DB)*(RAB + DA - DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DA + DB)*(RAB + DA + DB) + aLL*aLL;
    hmultipole -= 0.25*(3.0*(RAB + DA + DB)*(RAB + DA + DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA - DB)*(RAB - DA - DB) + aLL*aLL;
    hmultipole -= 0.25*(3.0*(RAB - DA - DB)*(RAB - DA - DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA + DB)*(RAB - DA + DB) + aLL*aLL;
    hmultipole += 0.25*(3.0*(RAB - DA + DB)*(RAB - DA + DB) - aux)/(aux*aux*sqrt(aux));
  }
  //[mux,Qxz]
  else if ((l1 == 1)&&(abs(m1) == 1)&&(l2 == 2)&&(abs(m2) == 1)&&(m1 == m2)) {
    double sq2 = 1.0/sqrt(2.0);
    aux = (RAB - sq2*DB)*(RAB - sq2*DB) + (DA - sq2*DB)*(DA - sq2*DB) + aLL*aLL;
    hmultipole = 0.25*(3.0*(RAB - sq2*DB)*(RAB - sq2*DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - sq2*DB)*(RAB - sq2*DB) + (DA + sq2*DB)*(DA + sq2*DB) + aLL*aLL;
    hmultipole -= 0.25*(3.0*(RAB - sq2*DB)*(RAB - sq2*DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + sq2*DB)*(RAB + sq2*DB) + (DA - sq2*DB)*(DA - sq2*DB) + aLL*aLL;
    hmultipole -= 0.25*(3.0*(RAB + sq2*DB)*(RAB + sq2*DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + sq2*DB)*(RAB + sq2*DB) + (DA + sq2*DB)*(DA + sq2*DB) + aLL*aLL;
    hmultipole += 0.25*(3.0*(RAB + sq2*DB)*(RAB + sq2*DB) - aux)/(aux*aux*sqrt(aux));
  }
  //[Qxz,mux]
  else if ((l1 == 2)&&(abs(m1) == 1)&&(l2 == 1)&&(abs(m2) == 1)&&(m1 == m2)) {
    double sq2 = 1.0/sqrt(2.0);
    DB = D[0];
    DA = D[1];
    aux = (RAB - sq2*DB)*(RAB - sq2*DB) + (DA - sq2*DB)*(DA - sq2*DB) + aLL*aLL;
    hmultipole = -0.25*(3.0*(RAB - sq2*DB)*(RAB - sq2*DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - sq2*DB)*(RAB - sq2*DB) + (DA + sq2*DB)*(DA + sq2*DB) + aLL*aLL;
    hmultipole += 0.25*(3.0*(RAB - sq2*DB)*(RAB - sq2*DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + sq2*DB)*(RAB + sq2*DB) + (DA - sq2*DB)*(DA - sq2*DB) + aLL*aLL;
    hmultipole += 0.25*(3.0*(RAB + sq2*DB)*(RAB + sq2*DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + sq2*DB)*(RAB + sq2*DB) + (DA + sq2*DB)*(DA + sq2*DB) + aLL*aLL;
    hmultipole -= 0.25*(3.0*(RAB + sq2*DB)*(RAB + sq2*DB) - aux)/(aux*aux*sqrt(aux));
  }
  //[muz,Qzz]
  else if ((l1 == 1)&&(m1 == 0)&&(l2 == 2)&&(m2 == 0)) {
    aux = (RAB + DA - DB)*(RAB + DA - DB) + aLL*aLL;
    hmultipole = 0.125*(3.0*(RAB + DA - DB)*(RAB + DA - DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA - DB)*(RAB - DA - DB) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB - DA - DB)*(RAB - DA - DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DA + DB)*(RAB + DA + DB) + aLL*aLL;
    hmultipole += 0.125*(3.0*(RAB + DA + DB)*(RAB + DA + DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA + DB)*(RAB - DA + DB) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB - DA + DB)*(RAB - DA + DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DA)*(RAB + DA) + DB*DB + aLL*aLL;
    hmultipole -= 0.25*(3.0*(RAB + DA)*(RAB + DA) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA)*(RAB - DA) + DB*DB + aLL*aLL;
    hmultipole += 0.25*(3.0*(RAB - DA)*(RAB - DA) - aux)/(aux*aux*sqrt(aux));
  }
  //[Qzz,muz]
  else if ((l1 == 2)&&(m1 == 0)&&(l2 == 1)&&(m2 == 0)) {
    DB = D[0];
    DA = D[1];
    aux = (RAB + DA - DB)*(RAB + DA - DB) + aLL*aLL;
    hmultipole = -0.125*(3.0*(RAB + DA - DB)*(RAB + DA - DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA - DB)*(RAB - DA - DB) + aLL*aLL;
    hmultipole += 0.125*(3.0*(RAB - DA - DB)*(RAB - DA - DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DA + DB)*(RAB + DA + DB) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB + DA + DB)*(RAB + DA + DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA + DB)*(RAB - DA + DB) + aLL*aLL;
    hmultipole += 0.125*(3.0*(RAB - DA + DB)*(RAB - DA + DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DA)*(RAB + DA) + DB*DB + aLL*aLL;
    hmultipole += 0.25*(3.0*(RAB + DA)*(RAB + DA) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA)*(RAB - DA) + DB*DB + aLL*aLL;
    hmultipole -= 0.25*(3.0*(RAB - DA)*(RAB - DA) - aux)/(aux*aux*sqrt(aux));
  }
  //[Qxy,Qxy]
  else if ((l1 == 2)&&(abs(m1) == 2)&&(l2 == 2)&&(abs(m2) == 2)&&(m1 == m2)) {
    aux = RAB*RAB + (DA - DB)*(DA - DB) + aLL*aLL;
    hmultipole = 0.25*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + (DA + DB)*(DA + DB) + aLL*aLL;
    hmultipole += 0.25*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + DA*DA + DB*DB + aLL*aLL;
    hmultipole -= 0.5*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
  }
  //[Qzz,Qzz]
  else if ((l1 == 2)&&(m1 == 0)&&(l2 == 2)&&(m2 == 0)) {
    aux = (RAB + DA - DB)*(RAB + DA - DB) + aLL*aLL;
    hmultipole = 0.0625*(3.0*(RAB + DA - DB)*(RAB + DA - DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DA + DB)*(RAB + DA + DB) + aLL*aLL;
    hmultipole += 0.0625*(3.0*(RAB + DA + DB)*(RAB + DA + DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA - DB)*(RAB - DA - DB) + aLL*aLL;
    hmultipole += 0.0625*(3.0*(RAB - DA - DB)*(RAB - DA - DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA + DB)*(RAB - DA + DB) + aLL*aLL;
    hmultipole += 0.0625*(3.0*(RAB - DA + DB)*(RAB - DA + DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DA)*(RAB + DA) + DB*DB + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB + DA)*(RAB + DA) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DA)*(RAB - DA) + DB*DB + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB - DA)*(RAB - DA) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + DB)*(RAB + DB) + DA*DA + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB + DB)*(RAB + DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - DB)*(RAB - DB) + DA*DA + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB - DB)*(RAB - DB) - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + DA*DA + DB*DB + aLL*aLL;
    hmultipole += 0.125*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + (DA - DB)*(DA - DB) + aLL*aLL;
    hmultipole += 0.0625*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
    aux = RAB*RAB + (DA + DB)*(DA + DB) + aLL*aLL;
    hmultipole += 0.0625*(3.0*RAB*RAB - aux)/(aux*aux*sqrt(aux));
  }
  //[Qxz,Qxz]
  else if ((l1 == 2)&&(abs(m1) == 1)&&(l2 == 2)&&(abs(m2) == 1)&&(m1 == m2)) {
    double sq2 = 1.0/sqrt(2.0);
    aux = (RAB + sq2*DA - sq2*DB)*(RAB + sq2*DA - sq2*DB) + 0.5*(DA - DB)*(DA - DB) + aLL*aLL;
    hmultipole = 0.125*(3.0*(RAB + sq2*DA - sq2*DB)*(RAB + sq2*DA - sq2*DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + sq2*DA - sq2*DB)*(RAB + sq2*DA - sq2*DB) + 0.5*(DA + DB)*(DA + DB) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB + sq2*DA - sq2*DB)*(RAB + sq2*DA - sq2*DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + sq2*DA + sq2*DB)*(RAB + sq2*DA + sq2*DB) + 0.5*(DA - DB)*(DA - DB) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB + sq2*DA + sq2*DB)*(RAB + sq2*DA + sq2*DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB + sq2*DA + sq2*DB)*(RAB + sq2*DA + sq2*DB) + 0.5*(DA + DB)*(DA + DB) + aLL*aLL;
    hmultipole += 0.125*(3.0*(RAB + sq2*DA + sq2*DB)*(RAB + sq2*DA + sq2*DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - sq2*DA - sq2*DB)*(RAB - sq2*DA - sq2*DB) + 0.5*(DA - DB)*(DA - DB) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB - sq2*DA - sq2*DB)*(RAB - sq2*DA - sq2*DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - sq2*DA - sq2*DB)*(RAB - sq2*DA - sq2*DB) + 0.5*(DA + DB)*(DA + DB) + aLL*aLL;
    hmultipole += 0.125*(3.0*(RAB - sq2*DA - sq2*DB)*(RAB - sq2*DA - sq2*DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - sq2*DA + sq2*DB)*(RAB - sq2*DA + sq2*DB) + 0.5*(DA - DB)*(DA - DB) + aLL*aLL;
    hmultipole += 0.125*(3.0*(RAB - sq2*DA + sq2*DB)*(RAB - sq2*DA + sq2*DB) - aux)/(aux*aux*sqrt(aux));
    aux = (RAB - sq2*DA + sq2*DB)*(RAB - sq2*DA + sq2*DB) + 0.5*(DA + DB)*(DA + DB) + aLL*aLL;
    hmultipole -= 0.125*(3.0*(RAB - sq2*DA + sq2*DB)*(RAB - sq2*DA + sq2*DB) - aux)/(aux*aux*sqrt(aux));
  }
  return hmultipole;
}
double clm_munu(int l1, int m1, int l2, int m2, int l, int m) {
  //function returning the factors required for the calculation of the Thiel and Voityuk integrals
  //correction for factors from Husch, Vaucher, Reiher, Semiempirical Molecular Orbital Models based on the Neglect of Diatomic Differential Overlap Approximation
  double clm = 1.0;             //this is the default value
  if ((l == 1)&&(m == 0)) {
    if ((l1 == 1)&&(m1 == 0)&&(l2 == 2)&&(m2 == 0)) {clm = sqrt(4.0/3.0);}
    else if ((l1 == 2)&&(m1 == 0)&&(l2 == 1)&&(m2 == 0)) {clm = sqrt(4.0/3.0);}
  }
  else if ((l == 2)&&(m == 0)) {
    clm = 4.0/3.0;
    if ((l1 == l2)&&(m1 == m2)&&(abs(m1) == 1)) {clm = 2.0*pow(-1.0,l1)/3.0;}
    if ((l1 == l2)&&(m1 == m2)&&(l1 == 2)&&(abs(m1) == 2)) {clm = -4.0/3.0;}
    else if ((l1 == 0)&&(l2 == 2)&&(m2 == 0)) {clm = sqrt(4.0/3.0);}
    else if ((l2 == 0)&&(l1 == 2)&&(m1 == 0)) {clm = sqrt(4.0/3.0);}
  }
  else if ((l == 1)&&(abs(m) == 1)) {
    if ((l1 == 1)&&(abs(m1) == 1)&&(l2 == 2)&&(m2 == 0)) {clm = -sqrt(1.0/3.0);}
    else if ((l2 == 1)&&(abs(m2) == 1)&&(l1 == 2)&&(m1 == 0)) {clm = sqrt(1.0/3.0);}
    else if ((l1 == 1)&&(m1 == -1)&&(l2 == 2)&&(m2 == 2)) {clm = -1.0;}
    else if ((l2 == 1)&&(m2 == -1)&&(l1 == 2)&&(m1 == 2)) {clm = -1.0;}
  }
  else if ((l == 2)&&(abs(m) == 1)) {
    if ((l1 == l2)&&(l2 == 2)&&(m1 == 0)&&(abs(m2) == 1)) {clm = sqrt(1.0/3.0);}
    else if ((l1 == l2)&&(l2 == 2)&&(abs(m1) == 1)&&(m2 == 0)) {clm = sqrt(1.0/3.0);}
    else if ((l1 == l2)&&(m1 == -1)&&(m2 == 2)&&(l1 == 2)) {clm = -1.0;}
    else if ((l1 == l2)&&(m2 == -1)&&(m1 == 2)&&(l1 == 2)) {clm = -1.0;}
  }
  else if ((l == 2)&&(abs(m) == 2)) {
    if ((l1 == l2)&&(l2 == 2)&&(m1 == 0)&&(abs(m2) == 2)) {clm = -sqrt(4.0/3.0);}
    else if ((l1 == l2)&&(l2 == 2)&&(abs(m1) == 2)&&(m2 == 0)) {clm = -sqrt(4.0/3.0);}
    else if ((l1 == l2)&&(m1 == m2)&&(l1 == 1)&&(abs(m1) == 1)) {clm = double(m1);}
    else if ((l1 == l2)&&(m1 == m2)&&(l1 == 2)&&(m1 == -1)) {clm = -1.0;}
  }
  return clm;
}
#endif //_2_Electron_Functions_Dewar_