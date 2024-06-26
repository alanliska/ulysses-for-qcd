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

#ifndef _Two_Center_Overlaps_STO_
#define _Two_Center_Overlaps_STO_
#include <vector>
#include <string>
#include <math.h>
#include "OvRot.hpp"
#include "OvFunctions.hpp"
#include "../math/MatrixPackage.hpp"
#include "../math/IntegralPackage.hpp"

//description:
//functions for the calculation of two-center overlap functions using STOs

double OvFactor(int nbra, int lbra, int mbra, double expbra, int nket, int lket, int mket, double expket, double RAB) {
  //function that calculates the overlap factor used in the calculation of the overlap matrix
  double ovfactor = 0.0;
  double x = 1.0;
  double y = 1.0;
  double w = 1.0; 
  //factors from Lofthus paper
  if ((lbra == 0)&&(lket == 0)) {w = 2.0;}
  else if ((lbra == 1)&&(lket == 1)) {
    if ((mbra == 0)&&(mket == 0)) {
      x = 3.0;
      w = 2.0;
    }
    else if ((abs(mbra) == 1)&&(abs(mket) == 1)) {
      x = 3.0;
      w = 4.0;
    }
    else {throw("2CenterOvSTO.hpp: OvFactor(): invalid combination for overlap integral (p-p)");}
  }
  else if ((lbra == 2)&&(lket == 2)) {
    if ((mbra == 0)&&(mket == 0)) {
      x = 5.0;
      w = 8.0;
    }
    else if ((abs(mbra) == 1)&&(abs(mket) == 1)) {
      x = 15.0;
      w = 4.0;
    }
    else if ((abs(mbra) == 2)&&(abs(mket) == 2)) {
      x = 15.0;
      w = 16.0;
    }
    else {throw("2CenterOvSTO.hpp: OvFactor(): invalid combination for overlap integral (d-d)");}
  }
  else if ((lbra == 3)&&(lket == 3)) {
    if ((mbra == 0)&&(mket == 0)) {
      x = 7.0;
      w = 8.0;
    }
    else if ((abs(mbra) == 1)&&(abs(mket) == 1)) {
      x = 21.0;
      w = 32.0;
    }
    else if ((abs(mbra) == 2)&&(abs(mket) == 2)) {
      x = 105.0;
      w = 16.0;
    }
    else if ((abs(mbra) == 3)&&(abs(mket) == 3)) {
      x = 35.0;
      w = 32.0;
    }
    else {throw("2CenterOvSTO.hpp: OvFactor(): invalid combination for overlap integral (f-f)");}
  }
  else if (((lbra == 0)&&(lket == 1))||((lbra == 1)&&(lket == 0))) {
    if (mket == 0) {
      y = 3.0;
      w = 2.0;
    }
    else {throw("2CenterOvSTO.hpp: OvFactor(): invalid combination for overlap integral (s-p/p-s)");}
  }
  else if (((lbra == 0)&&(lket == 2))||((lbra == 2)&&(lket == 0))) {
    if (mket == 0) {
      y = 5.0;
      w = 4.0;
    }
    else {throw("2CenterOvSTO.hpp: OvFactor(): invalid combination for overlap integral (s-d/d-s)");}
  }
  else if (((lbra == 0)&&(lket == 3))||((lbra == 3)&&(lket == 0))) {
    if (mket == 0) {
      y = 7.0;
      w = 4.0;
    }
    else {throw("2CenterOvSTO.hpp: OvFactor(): invalid combination for overlap integral (s-f/f-s)");}
  }
  else if (((lbra == 1)&&(lket == 2))||((lbra == 2)&&(lket == 1))) {
    if ((mbra == 0)&&(mket == 0)) {
      y = 15.0;
      w = 4.0;
    }
    else if ((abs(mbra) == 1)&&(abs(mket) == 1)) {
      x = 3.0;
      y = 5.0;
      w = 4.0;
    }
    else {throw("2CenterOvSTO.hpp: OvFactor(): invalid combination for overlap integral (p-d/d-p)");}
  }
  else if (((lbra == 1)&&(lket == 3))||((lbra == 3)&&(lket == 1))) {
    if ((mbra == 0)&&(mket == 0)) {
      y = 21.0;
      w = 4.0;
    }
    else if ((abs(mbra) == 1)&&(abs(mket) == 1)) {
      x = 3.0;
      y = 14.0;
      w = 16.0;
    }
    else {throw("2CenterOvSTO.hpp: OvFactor(): invalid combination for overlap integral (p-f/f-p)");}
  }
  else if (((lbra == 2)&&(lket == 3))||((lbra == 3)&&(lket == 2))) {
    if ((mbra == 0)&&(mket == 0)) {
      y = 35.0;
      w = 8.0;
    }
    else if ((abs(mbra) == 1)&&(abs(mket) == 1)) {
      x = 3.0;
      y = 70.0;
      w = 16.0;
    }
    else if ((abs(mbra) == 2)&&(abs(mket) == 2)) {
      x = 15.0;
      y = 7.0;
      w = 16.0;
    }
    else {throw("2CenterOvSTO.hpp: OvFactor(): invalid combination for overlap integral (d-f/f-d)");}
  }
  double t = FactorialRatio(2*std::max(nbra,nket),2*std::min(nbra,nket));
  ovfactor = (x/w)*sqrt(y/t)*pow(RAB,double(nbra + nket + 1))*sqrt(pow(expbra,double(2*nbra + 1))*pow(expket,double(2*nket + 1)))/Factorial(2*std::min(nbra,nket));
  return ovfactor;
}
void TwoCenterOv(matrixE & overlaps, size_t nao1, size_t nao2, double dist, std::vector<double> & exponents1, matrix<size_t> & QN1, std::vector<double> & exponents2, matrix<size_t> & QN2, std::vector<double> & orientation, bool rotate = true, double tolerance = 1e-6) {
  //implementation overlap matrix function; only up to shell f
  //calculation of the part of the overlap between atom A - atom B (A -> 1; B -> 2)
  overlaps.resize(nao1,nao2);
  overlaps.zero();
  //auxiliary variables for integrals
  double sigma_sigma = 0.0;
  double pi_pi = 0.0;
  double delta_delta = 0.0;
  double phi_phi = 0.0;
  //angles needed to reorient the integrals according to the position of atoms in the molecule
  double cost = orientation[2];
  double sint = sqrt(1 - cost*cost);
  double cosp = 1.0;
  double sinp = 0.0;
  if (fabs(sint) > tolerance) {
    cosp = orientation[0]/sint;
    sinp = orientation[1]/sint;
  }
  //matrices with overlap rotations
  matrixE SProt(1,1);
  matrixE SDrot(1,1);
  matrixE SFrot(1,1);
  matrixE PProt(1,1);
  matrixE PDrot(1,1);
  matrixE PFrot(1,1);
  matrixE DDrot(1,1);
  matrixE DFrot(1,1);
  matrixE FFrot(1,1);
  //other auxiliaries
  std::vector<std::vector<int> > polynomial;
  double a;
  double b;
  double c;
  double argA;
  double argB;
  std::vector<double> Ak;
  std::vector<double> Bk;
  double aux;
  //determine maximum azimuthal quantum number
  size_t maxL1 = 0;
  size_t maxL2 = 0;
  for (size_t idx = 0; idx < QN1.rows(); ++idx) {
    if (QN1(idx + 1,2) > maxL1) {maxL1 = QN1(idx + 1,2);}
  }
  for (size_t idx = 0; idx < QN2.rows(); ++idx) {
    if (QN2(idx + 1,2) > maxL2) {maxL2 = QN2(idx + 1,2);}
  }
  if ((maxL1 > 0)||(maxL2 > 0)) {
    //in this case we have p orbitals
    SProt = SPtransf(cost,sint,cosp,sinp,rotate);               //the orientation vector already contains this transformation, so this would not be actually necessary; however, a 3D vector is basically for free and like this I have debug options
    PProt = PPtransf(cost,sint,cosp,sinp,rotate);
  }
  if ((maxL1 > 1)||(maxL2 > 1)) {
    //in this case we have d orbitals
    SDrot = SDtransf(cost,sint,cosp,sinp,rotate);
    PDrot = PDtransf(cost,sint,cosp,sinp,rotate);
    DDrot = DDtransf(cost,sint,cosp,sinp,rotate);
  }
  if ((maxL1 > 2)||(maxL2 > 2)) {
    //in this case we have f orbitals
    SFrot = SFtransf(cost,sint,cosp,sinp,rotate);
    PFrot = PFtransf(cost,sint,cosp,sinp,rotate);
    DFrot = DFtransf(cost,sint,cosp,sinp,rotate);
    FFrot = FFtransf(cost,sint,cosp,sinp,rotate);
  }
  size_t iNA;
  size_t iNB;
  size_t iLA;
  size_t iLB;
  size_t irw = 1;
  size_t icl = 1;
  for (size_t isbra = 0; isbra < QN1.rows(); ++isbra) {           //loop over bra shells
    a = exponents1[isbra];
    if (a == 0.0) continue;
    iNA = QN1(isbra + 1,1) - 1;
    iLA = QN1(isbra + 1,2);
    icl = 1;
    for (size_t isket = 0; isket < QN2.rows(); ++isket) {         //loop over ket shells
      b = exponents2[isket];
      if (b == 0.0) continue;
      iNB = QN2(isket + 1,1) - 1;
      iLB = QN2(isket + 1,2);
      //getting integrals
      argA = 0.5*dist*(a + b);
      argB = 0.5*dist*(a - b);
      Ak = int_1_infty_xpk_expmyx_dx(argA,iNA + iNB + 2);
      Bk = int_m1_1_xpk_expmyx_dx_Pople(argB,iNA + iNB + 2);
      aux = OvFactor(iNA + 1,iLA,0,a,iNB + 1,iLB,0,b,dist);
      //different azimuthal quantum number cases
      if ((iLA == 0)&&(iLB == 0)) {
        //ss
        polynomial = SSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          overlaps(irw,icl) += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
      }
      else if ((iLA == 1)&&(iLB == 1)) {
        //pp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = PPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 9 dimensions
        overlaps(irw + 0,icl + 0) = -sigma_sigma*PProt(1,1) + pi_pi*(PProt(1,2) + PProt(1,3));      //<x|x>
        overlaps(irw + 0,icl + 1) = -sigma_sigma*PProt(2,1) + pi_pi*(PProt(2,2) + PProt(2,3));      //<x|y>
        overlaps(irw + 0,icl + 2) = -sigma_sigma*PProt(3,1) + pi_pi*(PProt(3,2) + PProt(3,3));      //<x|z>
        overlaps(irw + 1,icl + 0) = -sigma_sigma*PProt(2,1) + pi_pi*(PProt(2,2) + PProt(2,3));      //<y|x>
        overlaps(irw + 1,icl + 1) = -sigma_sigma*PProt(4,1) + pi_pi*(PProt(4,2) + PProt(4,3));      //<y|y>
        overlaps(irw + 1,icl + 2) = -sigma_sigma*PProt(5,1) + pi_pi*(PProt(5,2) + PProt(5,3));      //<y|z>
        overlaps(irw + 2,icl + 0) = -sigma_sigma*PProt(3,1) + pi_pi*(PProt(3,2) + PProt(3,3));      //<z|x>
        overlaps(irw + 2,icl + 1) = -sigma_sigma*PProt(5,1) + pi_pi*(PProt(5,2) + PProt(5,3));      //<z|y>
        overlaps(irw + 2,icl + 2) = -sigma_sigma*PProt(6,1) + pi_pi*(PProt(6,2) + PProt(6,3));      //<z|z>
      }
      else if ((iLA == 2)&&(iLB == 2)) {
        //dd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = DDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = DDdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5x5 dimensions (xz,yz,z2,xy,x2-y2)(xz,yz,z2,xy,x2-y2)
        overlaps(irw + 0,icl + 0) = sigma_sigma*DDrot(1,1)  - pi_pi*(DDrot(1,2)  + DDrot(1,3))  + delta_delta*(DDrot(1,4)  + DDrot(1,5));             //<xz|xz>
        overlaps(irw + 0,icl + 1) = sigma_sigma*DDrot(2,1)  - pi_pi*(DDrot(2,2)  + DDrot(2,3))  + delta_delta*(DDrot(2,4)  + DDrot(2,5));             //<xz|yz>
        overlaps(irw + 0,icl + 2) = sigma_sigma*DDrot(3,1)  - pi_pi*(DDrot(3,2)  + DDrot(3,3))  + delta_delta*(DDrot(3,4)  + DDrot(3,5));             //<xz|z2>
        overlaps(irw + 0,icl + 3) = sigma_sigma*DDrot(4,1)  - pi_pi*(DDrot(4,2)  + DDrot(4,3))  + delta_delta*(DDrot(4,4)  + DDrot(4,5));             //<xz|xy>
        overlaps(irw + 0,icl + 4) = sigma_sigma*DDrot(5,1)  - pi_pi*(DDrot(5,2)  + DDrot(5,3))  + delta_delta*(DDrot(5,4)  + DDrot(5,5));             //<xz|x2-y2>
        overlaps(irw + 1,icl + 0) = sigma_sigma*DDrot(2,1)  - pi_pi*(DDrot(2,2)  + DDrot(2,3))  + delta_delta*(DDrot(2,4)  + DDrot(2,5));             //<yz|xz>
        overlaps(irw + 1,icl + 1) = sigma_sigma*DDrot(6,1)  - pi_pi*(DDrot(6,2)  + DDrot(6,3))  + delta_delta*(DDrot(6,4)  + DDrot(6,5));             //<yz|yz>
        overlaps(irw + 1,icl + 2) = sigma_sigma*DDrot(7,1)  - pi_pi*(DDrot(7,2)  + DDrot(7,3))  + delta_delta*(DDrot(7,4)  + DDrot(7,5));             //<yz|z2>
        overlaps(irw + 1,icl + 3) = sigma_sigma*DDrot(8,1)  - pi_pi*(DDrot(8,2)  + DDrot(8,3))  + delta_delta*(DDrot(8,4)  + DDrot(8,5));             //<yz|xy>
        overlaps(irw + 1,icl + 4) = sigma_sigma*DDrot(9,1)  - pi_pi*(DDrot(9,2)  + DDrot(9,3))  + delta_delta*(DDrot(9,4)  + DDrot(9,5));             //<yz|x2-y2>
        overlaps(irw + 2,icl + 0) = sigma_sigma*DDrot(3,1)  - pi_pi*(DDrot(3,2)  + DDrot(3,3))  + delta_delta*(DDrot(3,4)  + DDrot(3,5));             //<z2|xz>
        overlaps(irw + 2,icl + 1) = sigma_sigma*DDrot(7,1)  - pi_pi*(DDrot(7,2)  + DDrot(7,3))  + delta_delta*(DDrot(7,4)  + DDrot(7,5));             //<z2|yz>
        overlaps(irw + 2,icl + 2) = sigma_sigma*DDrot(10,1) - pi_pi*(DDrot(10,2) + DDrot(10,3)) + delta_delta*(DDrot(10,4) + DDrot(10,5));            //<z2|z2>
        overlaps(irw + 2,icl + 3) = sigma_sigma*DDrot(11,1) - pi_pi*(DDrot(11,2) + DDrot(11,3)) + delta_delta*(DDrot(11,4) + DDrot(11,5));            //<z2|xy>
        overlaps(irw + 2,icl + 4) = sigma_sigma*DDrot(12,1) - pi_pi*(DDrot(12,2) + DDrot(12,3)) + delta_delta*(DDrot(12,4) + DDrot(12,5));            //<z2|x2-y2>
        overlaps(irw + 3,icl + 0) = sigma_sigma*DDrot(4,1)  - pi_pi*(DDrot(4,2)  + DDrot(4,3))  + delta_delta*(DDrot(4,4)  + DDrot(4,5));             //<xy|xz>
        overlaps(irw + 3,icl + 1) = sigma_sigma*DDrot(8,1)  - pi_pi*(DDrot(8,2)  + DDrot(8,3))  + delta_delta*(DDrot(8,4)  + DDrot(8,5));             //<xy|yz>
        overlaps(irw + 3,icl + 2) = sigma_sigma*DDrot(11,1) - pi_pi*(DDrot(11,2) + DDrot(11,3)) + delta_delta*(DDrot(11,4) + DDrot(11,5));            //<xy|z2>
        overlaps(irw + 3,icl + 3) = sigma_sigma*DDrot(13,1) - pi_pi*(DDrot(13,2) + DDrot(13,3)) + delta_delta*(DDrot(13,4) + DDrot(13,5));            //<xy|xy>
        overlaps(irw + 3,icl + 4) = sigma_sigma*DDrot(14,1) - pi_pi*(DDrot(14,2) + DDrot(14,3)) + delta_delta*(DDrot(14,4) + DDrot(14,5));            //<xy|x2-y2>
        overlaps(irw + 4,icl + 0) = sigma_sigma*DDrot(5,1)  - pi_pi*(DDrot(5,2)  + DDrot(5,3))  + delta_delta*(DDrot(5,4)  + DDrot(5,5));             //<x2-y2|xz>
        overlaps(irw + 4,icl + 1) = sigma_sigma*DDrot(9,1)  - pi_pi*(DDrot(9,2)  + DDrot(9,3))  + delta_delta*(DDrot(9,4)  + DDrot(9,5));             //<x2-y2|yz>
        overlaps(irw + 4,icl + 2) = sigma_sigma*DDrot(12,1) - pi_pi*(DDrot(12,2) + DDrot(12,3)) + delta_delta*(DDrot(12,4) + DDrot(12,5));            //<x2-y2|z2>
        overlaps(irw + 4,icl + 3) = sigma_sigma*DDrot(14,1) - pi_pi*(DDrot(14,2) + DDrot(14,3)) + delta_delta*(DDrot(14,4) + DDrot(14,5));            //<x2-y2|xy>
        overlaps(irw + 4,icl + 4) = sigma_sigma*DDrot(15,1) - pi_pi*(DDrot(15,2) + DDrot(15,3)) + delta_delta*(DDrot(15,4) + DDrot(15,5));            //<x2-y2|x2-y2>
      }
      else if ((iLA == 3)&&(iLB == 3)) {
        //ff
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = FFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = FFdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //phi integral
        aux = OvFactor(iNA + 1,iLA,3,a,iNB + 1,iLB,3,b,dist);
        phi_phi = 0.0;
        polynomial = FFphiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          phi_phi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 7x7 dimensions (fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        overlaps(irw + 0,icl + 0) = -sigma_sigma*FFrot(1,1) + pi_pi*(FFrot(1,2) + FFrot(1,3)) - delta_delta*(FFrot(1,4) + FFrot(1,5)) + phi_phi*(FFrot(1,6) + FFrot(1,7));
        overlaps(irw + 0,icl + 1) = -sigma_sigma*FFrot(2,1) + pi_pi*(FFrot(2,2) + FFrot(2,3)) - delta_delta*(FFrot(2,4) + FFrot(2,5)) + phi_phi*(FFrot(2,6) + FFrot(2,7));
        overlaps(irw + 0,icl + 2) = -sigma_sigma*FFrot(3,1) + pi_pi*(FFrot(3,2) + FFrot(3,3)) - delta_delta*(FFrot(3,4) + FFrot(3,5)) + phi_phi*(FFrot(3,6) + FFrot(3,7));
        overlaps(irw + 0,icl + 3) = -sigma_sigma*FFrot(4,1) + pi_pi*(FFrot(4,2) + FFrot(4,3)) - delta_delta*(FFrot(4,4) + FFrot(4,5)) + phi_phi*(FFrot(4,6) + FFrot(4,7));
        overlaps(irw + 0,icl + 4) = -sigma_sigma*FFrot(5,1) + pi_pi*(FFrot(5,2) + FFrot(5,3)) - delta_delta*(FFrot(5,4) + FFrot(5,5)) + phi_phi*(FFrot(5,6) + FFrot(5,7));
        overlaps(irw + 0,icl + 5) = -sigma_sigma*FFrot(6,1) + pi_pi*(FFrot(6,2) + FFrot(6,3)) - delta_delta*(FFrot(6,4) + FFrot(6,5)) + phi_phi*(FFrot(6,6) + FFrot(6,7));
        overlaps(irw + 0,icl + 6) = -sigma_sigma*FFrot(7,1) + pi_pi*(FFrot(7,2) + FFrot(7,3)) - delta_delta*(FFrot(7,4) + FFrot(7,5)) + phi_phi*(FFrot(7,6) + FFrot(7,7));
        overlaps(irw + 1,icl + 0) = -sigma_sigma*FFrot(2,1) + pi_pi*(FFrot(2,2) + FFrot(2,3)) - delta_delta*(FFrot(2,4) + FFrot(2,5)) + phi_phi*(FFrot(2,6) + FFrot(2,7));
        overlaps(irw + 1,icl + 1) = -sigma_sigma*FFrot(8,1) + pi_pi*(FFrot(8,2) + FFrot(8,3)) - delta_delta*(FFrot(8,4) + FFrot(8,5)) + phi_phi*(FFrot(8,6) + FFrot(8,7));
        overlaps(irw + 1,icl + 2) = -sigma_sigma*FFrot(9,1) + pi_pi*(FFrot(9,2) + FFrot(9,3)) - delta_delta*(FFrot(9,4) + FFrot(9,5)) + phi_phi*(FFrot(9,6) + FFrot(9,7));
        overlaps(irw + 1,icl + 3) = -sigma_sigma*FFrot(10,1) + pi_pi*(FFrot(10,2) + FFrot(10,3)) - delta_delta*(FFrot(10,4) + FFrot(10,5)) + phi_phi*(FFrot(10,6) + FFrot(10,7));
        overlaps(irw + 1,icl + 4) = -sigma_sigma*FFrot(11,1) + pi_pi*(FFrot(11,2) + FFrot(11,3)) - delta_delta*(FFrot(11,4) + FFrot(11,5)) + phi_phi*(FFrot(11,6) + FFrot(11,7));
        overlaps(irw + 1,icl + 5) = -sigma_sigma*FFrot(12,1) + pi_pi*(FFrot(12,2) + FFrot(12,3)) - delta_delta*(FFrot(12,4) + FFrot(12,5)) + phi_phi*(FFrot(12,6) + FFrot(12,7));
        overlaps(irw + 1,icl + 6) = -sigma_sigma*FFrot(13,1) + pi_pi*(FFrot(13,2) + FFrot(13,3)) - delta_delta*(FFrot(13,4) + FFrot(13,5)) + phi_phi*(FFrot(13,6) + FFrot(13,7));
        overlaps(irw + 2,icl + 0) = -sigma_sigma*FFrot(3,1) +  pi_pi*(FFrot(3,2) + FFrot(3,3)) - delta_delta*(FFrot(3,4) + FFrot(3,5)) + phi_phi*(FFrot(3,6) + FFrot(3,7));
        overlaps(irw + 2,icl + 1) = -sigma_sigma*FFrot(9,1) +  pi_pi*(FFrot(9,2) + FFrot(9,3)) - delta_delta*(FFrot(9,4) + FFrot(9,5)) + phi_phi*(FFrot(9,6) + FFrot(9,7));
        overlaps(irw + 2,icl + 2) = -sigma_sigma*FFrot(14,1) + pi_pi*(FFrot(14,2) + FFrot(14,3)) - delta_delta*(FFrot(14,4) + FFrot(14,5)) + phi_phi*(FFrot(14,6) + FFrot(14,7));
        overlaps(irw + 2,icl + 3) = -sigma_sigma*FFrot(15,1) + pi_pi*(FFrot(15,2) + FFrot(15,3)) - delta_delta*(FFrot(15,4) + FFrot(15,5)) + phi_phi*(FFrot(15,6) + FFrot(15,7));
        overlaps(irw + 2,icl + 4) = -sigma_sigma*FFrot(16,1) + pi_pi*(FFrot(16,2) + FFrot(16,3)) - delta_delta*(FFrot(16,4) + FFrot(16,5)) + phi_phi*(FFrot(16,6) + FFrot(16,7));
        overlaps(irw + 2,icl + 5) = -sigma_sigma*FFrot(17,1) + pi_pi*(FFrot(17,2) + FFrot(17,3)) - delta_delta*(FFrot(17,4) + FFrot(17,5)) + phi_phi*(FFrot(17,6) + FFrot(17,7));
        overlaps(irw + 2,icl + 6) = -sigma_sigma*FFrot(18,1) + pi_pi*(FFrot(18,2) + FFrot(18,3)) - delta_delta*(FFrot(18,4) + FFrot(18,5)) + phi_phi*(FFrot(18,6) + FFrot(18,7));
        overlaps(irw + 3,icl + 0) = -sigma_sigma*FFrot(4,1) +  pi_pi*(FFrot(4,2) + FFrot(4,3)) - delta_delta*(FFrot(4,4) + FFrot(4,5)) + phi_phi*(FFrot(4,6) + FFrot(4,7));
        overlaps(irw + 3,icl + 1) = -sigma_sigma*FFrot(10,1) + pi_pi*(FFrot(10,2) + FFrot(10,3)) - delta_delta*(FFrot(10,4) + FFrot(10,5)) + phi_phi*(FFrot(10,6) + FFrot(10,7));
        overlaps(irw + 3,icl + 2) = -sigma_sigma*FFrot(15,1) + pi_pi*(FFrot(15,2) + FFrot(15,3)) - delta_delta*(FFrot(15,4) + FFrot(15,5)) + phi_phi*(FFrot(15,6) + FFrot(15,7));
        overlaps(irw + 3,icl + 3) = -sigma_sigma*FFrot(19,1) + pi_pi*(FFrot(19,2) + FFrot(19,3)) - delta_delta*(FFrot(19,4) + FFrot(19,5)) + phi_phi*(FFrot(19,6) + FFrot(19,7));
        overlaps(irw + 3,icl + 4) = -sigma_sigma*FFrot(20,1) + pi_pi*(FFrot(20,2) + FFrot(20,3)) - delta_delta*(FFrot(20,4) + FFrot(20,5)) + phi_phi*(FFrot(20,6) + FFrot(20,7));
        overlaps(irw + 3,icl + 5) = -sigma_sigma*FFrot(21,1) + pi_pi*(FFrot(21,2) + FFrot(21,3)) - delta_delta*(FFrot(21,4) + FFrot(21,5)) + phi_phi*(FFrot(21,6) + FFrot(21,7));
        overlaps(irw + 3,icl + 6) = -sigma_sigma*FFrot(22,1) + pi_pi*(FFrot(22,2) + FFrot(22,3)) - delta_delta*(FFrot(22,4) + FFrot(22,5)) + phi_phi*(FFrot(22,6) + FFrot(22,7));
        overlaps(irw + 4,icl + 0) = -sigma_sigma*FFrot(5,1) +  pi_pi*(FFrot(5,2) + FFrot(5,3)) - delta_delta*(FFrot(5,4) + FFrot(5,5)) + phi_phi*(FFrot(5,6) + FFrot(5,7));
        overlaps(irw + 4,icl + 1) = -sigma_sigma*FFrot(11,1) + pi_pi*(FFrot(11,2) + FFrot(11,3)) - delta_delta*(FFrot(11,4) + FFrot(11,5)) + phi_phi*(FFrot(11,6) + FFrot(11,7));
        overlaps(irw + 4,icl + 2) = -sigma_sigma*FFrot(16,1) + pi_pi*(FFrot(16,2) + FFrot(16,3)) - delta_delta*(FFrot(16,4) + FFrot(16,5)) + phi_phi*(FFrot(16,6) + FFrot(16,7));
        overlaps(irw + 4,icl + 3) = -sigma_sigma*FFrot(20,1) + pi_pi*(FFrot(20,2) + FFrot(20,3)) - delta_delta*(FFrot(20,4) + FFrot(20,5)) + phi_phi*(FFrot(20,6) + FFrot(20,7));
        overlaps(irw + 4,icl + 4) = -sigma_sigma*FFrot(23,1) + pi_pi*(FFrot(23,2) + FFrot(23,3)) - delta_delta*(FFrot(23,4) + FFrot(23,5)) + phi_phi*(FFrot(23,6) + FFrot(23,7));
        overlaps(irw + 4,icl + 5) = -sigma_sigma*FFrot(24,1) + pi_pi*(FFrot(24,2) + FFrot(24,3)) - delta_delta*(FFrot(24,4) + FFrot(24,5)) + phi_phi*(FFrot(24,6) + FFrot(24,7));
        overlaps(irw + 4,icl + 6) = -sigma_sigma*FFrot(25,1) + pi_pi*(FFrot(25,2) + FFrot(25,3)) - delta_delta*(FFrot(25,4) + FFrot(25,5)) + phi_phi*(FFrot(25,6) + FFrot(25,7));
        overlaps(irw + 5,icl + 0) = -sigma_sigma*FFrot(6,1) +  pi_pi*(FFrot(6,2) + FFrot(6,3)) - delta_delta*(FFrot(6,4) + FFrot(6,5)) + phi_phi*(FFrot(6,6) + FFrot(6,7));
        overlaps(irw + 5,icl + 1) = -sigma_sigma*FFrot(12,1) + pi_pi*(FFrot(12,2) + FFrot(12,3)) - delta_delta*(FFrot(12,4) + FFrot(12,5)) + phi_phi*(FFrot(12,6) + FFrot(12,7));
        overlaps(irw + 5,icl + 2) = -sigma_sigma*FFrot(17,1) + pi_pi*(FFrot(17,2) + FFrot(17,3)) - delta_delta*(FFrot(17,4) + FFrot(17,5)) + phi_phi*(FFrot(17,6) + FFrot(17,7));
        overlaps(irw + 5,icl + 3) = -sigma_sigma*FFrot(21,1) + pi_pi*(FFrot(21,2) + FFrot(21,3)) - delta_delta*(FFrot(21,4) + FFrot(21,5)) + phi_phi*(FFrot(21,6) + FFrot(21,7));
        overlaps(irw + 5,icl + 4) = -sigma_sigma*FFrot(24,1) + pi_pi*(FFrot(24,2) + FFrot(24,3)) - delta_delta*(FFrot(24,4) + FFrot(24,5)) + phi_phi*(FFrot(24,6) + FFrot(24,7));
        overlaps(irw + 5,icl + 5) = -sigma_sigma*FFrot(26,1) + pi_pi*(FFrot(26,2) + FFrot(26,3)) - delta_delta*(FFrot(26,4) + FFrot(26,5)) + phi_phi*(FFrot(26,6) + FFrot(26,7));
        overlaps(irw + 5,icl + 6) = -sigma_sigma*FFrot(27,1) + pi_pi*(FFrot(27,2) + FFrot(27,3)) - delta_delta*(FFrot(27,4) + FFrot(27,5)) + phi_phi*(FFrot(27,6) + FFrot(27,7));
        overlaps(irw + 6,icl + 0) = -sigma_sigma*FFrot(7,1) +  pi_pi*(FFrot(7,2) + FFrot(7,3)) - delta_delta*(FFrot(7,4) + FFrot(7,5)) + phi_phi*(FFrot(7,6) + FFrot(7,7));
        overlaps(irw + 6,icl + 1) = -sigma_sigma*FFrot(13,1) + pi_pi*(FFrot(13,2) + FFrot(13,3)) - delta_delta*(FFrot(13,4) + FFrot(13,5)) + phi_phi*(FFrot(13,6) + FFrot(13,7));
        overlaps(irw + 6,icl + 2) = -sigma_sigma*FFrot(18,1) + pi_pi*(FFrot(18,2) + FFrot(18,3)) - delta_delta*(FFrot(18,4) + FFrot(18,5)) + phi_phi*(FFrot(18,6) + FFrot(18,7));
        overlaps(irw + 6,icl + 3) = -sigma_sigma*FFrot(22,1) + pi_pi*(FFrot(22,2) + FFrot(22,3)) - delta_delta*(FFrot(22,4) + FFrot(22,5)) + phi_phi*(FFrot(22,6) + FFrot(22,7));
        overlaps(irw + 6,icl + 4) = -sigma_sigma*FFrot(25,1) + pi_pi*(FFrot(25,2) + FFrot(25,3)) - delta_delta*(FFrot(25,4) + FFrot(25,5)) + phi_phi*(FFrot(25,6) + FFrot(25,7));
        overlaps(irw + 6,icl + 5) = -sigma_sigma*FFrot(27,1) + pi_pi*(FFrot(27,2) + FFrot(27,3)) - delta_delta*(FFrot(27,4) + FFrot(27,5)) + phi_phi*(FFrot(27,6) + FFrot(27,7));
        overlaps(irw + 6,icl + 6) = -sigma_sigma*FFrot(28,1) + pi_pi*(FFrot(28,2) + FFrot(28,3)) - delta_delta*(FFrot(28,4) + FFrot(28,5)) + phi_phi*(FFrot(28,6) + FFrot(28,7));
      }
      else if ((iLA == 0)&&(iLB == 1)) {
        //sp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SPOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3 dimensions (px,py,pz)
        overlaps(irw,icl + 0) = -SProt(1,1)*sigma_sigma;
        overlaps(irw,icl + 1) = -SProt(2,1)*sigma_sigma;
        overlaps(irw,icl + 2) = -SProt(3,1)*sigma_sigma;
      }
      else if ((iLA == 1)&&(iLB == 0)) {
        //ps
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3 dimensions (px,py,pz)
        overlaps(irw + 0,icl) = SProt(1,1)*sigma_sigma;
        overlaps(irw + 1,icl) = SProt(2,1)*sigma_sigma;
        overlaps(irw + 2,icl) = SProt(3,1)*sigma_sigma;
      }
      else if ((iLA == 0)&&(iLB == 2)) {
        //sd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SDOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          overlaps(irw,icl + index) = SDrot(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 2)&&(iLB == 0)) {
        //ds
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          overlaps(irw + index,icl) = SDrot(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 0)&&(iLB == 3)) {
        //sf
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SFOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 7 dimensions (fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 7; ++index) {
          overlaps(irw,icl + index) = -SFrot(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 3)&&(iLB == 0)) {
        //fs
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 7; ++index) {
          overlaps(irw + index,icl) = SFrot(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 1)&&(iLB == 2)) {
        //pd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = PDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 3; ++index2) {
            overlaps(irw + index2,icl + index) = PDrot(1+index2+3*index,1)*sigma_sigma - (PDrot(1+index2+3*index,2) + PDrot(1+index2+3*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 2)&&(iLB == 1)) {
        //dp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = DPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 3; ++index2) {
            overlaps(irw + index,icl + index2) = -PDrot(1+index2+3*index,1)*sigma_sigma + (PDrot(1+index2+3*index,2) + PDrot(1+index2+3*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 1)&&(iLB == 3)) {
        //pf
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = PFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 3; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index,icl + index2) = -PFrot(1+index2+7*index,1)*sigma_sigma + (PFrot(1+index2+7*index,2) + PFrot(1+index2+7*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 3)&&(iLB == 1)) {
        //fp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = FPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 3; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index2,icl + index) = -PFrot(1+index2+7*index,1)*sigma_sigma + (PFrot(1+index2+7*index,2) + PFrot(1+index2+7*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 2)&&(iLB == 3)) {
        //df
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = DFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = DFdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (xz,yz,z2,xy,x2-y2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index,icl + index2) = -DFrot(1+index2+7*index,1)*sigma_sigma + (DFrot(1+index2+7*index,2) + DFrot(1+index2+7*index,3))*pi_pi - (DFrot(1+index2+7*index,4) + DFrot(1+index2+7*index,5))*delta_delta;
          }
        }
      }
      else if ((iLA == 3)&&(iLB == 2)) {
        //fd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = FDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = FDdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (xz,yz,z2,xy,x2-y2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index2,icl + index) = DFrot(1+index2+7*index,1)*sigma_sigma - (DFrot(1+index2+7*index,2) + DFrot(1+index2+7*index,3))*pi_pi + (DFrot(1+index2+7*index,4) + DFrot(1+index2+7*index,5))*delta_delta;
          }
        }
      }
      icl += 2*QN2(isket + 1,2) + 1;
    }
    irw += 2*QN1(isbra + 1,2) + 1;
  }
}
void TwoCenterOv_dR(matrixE & overlapdR, size_t nao1, size_t nao2, double dist, std::vector<double> & exponents1, matrix<size_t> & QN1, std::vector<double> & exponents2, matrix<size_t> & QN2, std::vector<double> & orientation, bool rotate = true, double tolerance = 1e-6) {
  //implementation of first-derivative of overlap matrix with respect to internuclear distance; only up to shell f
  //calculation of the derivative of the overlap function between atoms A and B (A -> 1; B -> 2) with respect to internuclear distance
  overlapdR.resize(nao1,nao2);
  overlapdR.zero();
  //auxiliary variables for integrals
  double sigma_sigma = 0.0;
  double pi_pi = 0.0;
  double delta_delta = 0.0;
  double phi_phi = 0.0;
  //angles needed to reorient the integrals according to the position of atoms in the molecule
  double cost = orientation[2];
  double sint = sqrt(1 - cost*cost);
  double cosp = 1.0;
  double sinp = 0.0;
  if (fabs(sint) > tolerance) {
    cosp = orientation[0]/sint;
    sinp = orientation[1]/sint;
  }
  //matrices with overlap rotations
  matrixE SProt(1,1);
  matrixE SDrot(1,1);
  matrixE SFrot(1,1);
  matrixE PProt(1,1);
  matrixE PDrot(1,1);
  matrixE PFrot(1,1);
  matrixE DDrot(1,1);
  matrixE DFrot(1,1);
  matrixE FFrot(1,1);
  //other auxiliaries
  std::vector<std::vector<int> > polynomial;
  double a;
  double b;
  double c;
  double argA;
  double argB;
  std::vector<double> Ak;
  std::vector<double> Bk;
  double aux;
  //determine maximum azimuthal quantum number
  size_t maxL1 = 0;
  size_t maxL2 = 0;
  for (size_t idx = 0; idx < QN1.rows(); ++idx) {
    if (QN1(idx + 1,2) > maxL1) {maxL1 = QN1(idx + 1,2);}
  }
  for (size_t idx = 0; idx < QN2.rows(); ++idx) {
    if (QN2(idx + 1,2) > maxL2) {maxL2 = QN2(idx + 1,2);}
  }
  if ((maxL1 > 0)||(maxL2 > 0)) {
    //in this case we have p orbitals
    SProt = SPtransf(cost,sint,cosp,sinp,rotate);               //the orientation vector already contains this transformation, so this would not be actually necessary; however, a 3D vector is basically for free and like this I have debug options
    PProt = PPtransf(cost,sint,cosp,sinp,rotate);
  }
  if ((maxL1 > 1)||(maxL2 > 1)) {
    //in this case we have d orbitals
    SDrot = SDtransf(cost,sint,cosp,sinp,rotate);
    PDrot = PDtransf(cost,sint,cosp,sinp,rotate);
    DDrot = DDtransf(cost,sint,cosp,sinp,rotate);
  }
  if ((maxL1 > 2)||(maxL2 > 2)) {
    //in this case we have f orbitals
    SFrot = SFtransf(cost,sint,cosp,sinp,rotate);
    PFrot = PFtransf(cost,sint,cosp,sinp,rotate);
    DFrot = DFtransf(cost,sint,cosp,sinp,rotate);
    FFrot = FFtransf(cost,sint,cosp,sinp,rotate);
  }
  size_t iNA;
  size_t iNB;
  size_t iLA;
  size_t iLB;
  size_t irw = 1;
  size_t icl = 1;
  for (size_t isbra = 0; isbra < QN1.rows(); ++isbra) {           //loop over bra shells
    a = exponents1[isbra];
    if (a == 0.0) continue;
    iNA = QN1(isbra + 1,1) - 1;             //note that iNA is subtracted -1 !
    iLA = QN1(isbra + 1,2);
    icl = 1;
    for (size_t isket = 0; isket < QN2.rows(); ++isket) {         //loop over ket shells
      b = exponents2[isket];
      if (b == 0.0) continue;
      iNB = QN2(isket + 1,1) - 1;           //note that iNB is subtracted -1 !
      iLB = QN2(isket + 1,2);
      //getting integrals
      argA = 0.5*dist*(a + b);
      argB = 0.5*dist*(a - b);
      Ak = int_1_infty_xpk_expmyx_dx(argA,iNA + iNB + 3);
      Bk = int_m1_1_xpk_expmyx_dx_Pople(argB,iNA + iNB + 3);
      aux = OvFactor(iNA + 1,iLA,0,a,iNB + 1,iLB,0,b,dist)/dist;
      //different azimuthal quantum number cases
      if ((iLA == 0)&&(iLB == 0)) {
        //ss
        polynomial = SSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          overlapdR(irw,icl) += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          overlapdR(irw,icl) -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          overlapdR(irw,icl) -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
      }
      else if ((iLA == 1)&&(iLB == 1)) {
        //pp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = PPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 9 dimensions
        overlapdR(irw + 0,icl + 0) = -sigma_sigma*PProt(1,1) + pi_pi*(PProt(1,2) + PProt(1,3));      //<x|x>
        overlapdR(irw + 0,icl + 1) = -sigma_sigma*PProt(2,1) + pi_pi*(PProt(2,2) + PProt(2,3));      //<x|y>
        overlapdR(irw + 0,icl + 2) = -sigma_sigma*PProt(3,1) + pi_pi*(PProt(3,2) + PProt(3,3));      //<x|z>
        overlapdR(irw + 1,icl + 0) = -sigma_sigma*PProt(2,1) + pi_pi*(PProt(2,2) + PProt(2,3));      //<y|x>
        overlapdR(irw + 1,icl + 1) = -sigma_sigma*PProt(4,1) + pi_pi*(PProt(4,2) + PProt(4,3));      //<y|y>
        overlapdR(irw + 1,icl + 2) = -sigma_sigma*PProt(5,1) + pi_pi*(PProt(5,2) + PProt(5,3));      //<y|z>
        overlapdR(irw + 2,icl + 0) = -sigma_sigma*PProt(3,1) + pi_pi*(PProt(3,2) + PProt(3,3));      //<z|x>
        overlapdR(irw + 2,icl + 1) = -sigma_sigma*PProt(5,1) + pi_pi*(PProt(5,2) + PProt(5,3));      //<z|y>
        overlapdR(irw + 2,icl + 2) = -sigma_sigma*PProt(6,1) + pi_pi*(PProt(6,2) + PProt(6,3));      //<z|z>
      }
      else if ((iLA == 2)&&(iLB == 2)) {
        //dd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = DDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist)/dist;
        delta_delta = 0.0;
        polynomial = DDdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 5x5 dimensions (xz,yz,z2,xy,x2-y2)(xz,yz,z2,xy,x2-y2)
        overlapdR(irw + 0,icl + 0) = sigma_sigma*DDrot(1,1)  - pi_pi*(DDrot(1,2)  + DDrot(1,3))  + delta_delta*(DDrot(1,4)  + DDrot(1,5));             //<xz|xz>
        overlapdR(irw + 0,icl + 1) = sigma_sigma*DDrot(2,1)  - pi_pi*(DDrot(2,2)  + DDrot(2,3))  + delta_delta*(DDrot(2,4)  + DDrot(2,5));             //<xz|yz>
        overlapdR(irw + 0,icl + 2) = sigma_sigma*DDrot(3,1)  - pi_pi*(DDrot(3,2)  + DDrot(3,3))  + delta_delta*(DDrot(3,4)  + DDrot(3,5));             //<xz|z2>
        overlapdR(irw + 0,icl + 3) = sigma_sigma*DDrot(4,1)  - pi_pi*(DDrot(4,2)  + DDrot(4,3))  + delta_delta*(DDrot(4,4)  + DDrot(4,5));             //<xz|xy>
        overlapdR(irw + 0,icl + 4) = sigma_sigma*DDrot(5,1)  - pi_pi*(DDrot(5,2)  + DDrot(5,3))  + delta_delta*(DDrot(5,4)  + DDrot(5,5));             //<xz|x2-y2>
        overlapdR(irw + 1,icl + 0) = sigma_sigma*DDrot(2,1)  - pi_pi*(DDrot(2,2)  + DDrot(2,3))  + delta_delta*(DDrot(2,4)  + DDrot(2,5));             //<yz|xz>
        overlapdR(irw + 1,icl + 1) = sigma_sigma*DDrot(6,1)  - pi_pi*(DDrot(6,2)  + DDrot(6,3))  + delta_delta*(DDrot(6,4)  + DDrot(6,5));             //<yz|yz>
        overlapdR(irw + 1,icl + 2) = sigma_sigma*DDrot(7,1)  - pi_pi*(DDrot(7,2)  + DDrot(7,3))  + delta_delta*(DDrot(7,4)  + DDrot(7,5));             //<yz|z2>
        overlapdR(irw + 1,icl + 3) = sigma_sigma*DDrot(8,1)  - pi_pi*(DDrot(8,2)  + DDrot(8,3))  + delta_delta*(DDrot(8,4)  + DDrot(8,5));             //<yz|xy>
        overlapdR(irw + 1,icl + 4) = sigma_sigma*DDrot(9,1)  - pi_pi*(DDrot(9,2)  + DDrot(9,3))  + delta_delta*(DDrot(9,4)  + DDrot(9,5));             //<yz|x2-y2>
        overlapdR(irw + 2,icl + 0) = sigma_sigma*DDrot(3,1)  - pi_pi*(DDrot(3,2)  + DDrot(3,3))  + delta_delta*(DDrot(3,4)  + DDrot(3,5));             //<z2|xz>
        overlapdR(irw + 2,icl + 1) = sigma_sigma*DDrot(7,1)  - pi_pi*(DDrot(7,2)  + DDrot(7,3))  + delta_delta*(DDrot(7,4)  + DDrot(7,5));             //<z2|yz>
        overlapdR(irw + 2,icl + 2) = sigma_sigma*DDrot(10,1) - pi_pi*(DDrot(10,2) + DDrot(10,3)) + delta_delta*(DDrot(10,4) + DDrot(10,5));            //<z2|z2>
        overlapdR(irw + 2,icl + 3) = sigma_sigma*DDrot(11,1) - pi_pi*(DDrot(11,2) + DDrot(11,3)) + delta_delta*(DDrot(11,4) + DDrot(11,5));            //<z2|xy>
        overlapdR(irw + 2,icl + 4) = sigma_sigma*DDrot(12,1) - pi_pi*(DDrot(12,2) + DDrot(12,3)) + delta_delta*(DDrot(12,4) + DDrot(12,5));            //<z2|x2-y2>
        overlapdR(irw + 3,icl + 0) = sigma_sigma*DDrot(4,1)  - pi_pi*(DDrot(4,2)  + DDrot(4,3))  + delta_delta*(DDrot(4,4)  + DDrot(4,5));             //<xy|xz>
        overlapdR(irw + 3,icl + 1) = sigma_sigma*DDrot(8,1)  - pi_pi*(DDrot(8,2)  + DDrot(8,3))  + delta_delta*(DDrot(8,4)  + DDrot(8,5));             //<xy|yz>
        overlapdR(irw + 3,icl + 2) = sigma_sigma*DDrot(11,1) - pi_pi*(DDrot(11,2) + DDrot(11,3)) + delta_delta*(DDrot(11,4) + DDrot(11,5));            //<xy|z2>
        overlapdR(irw + 3,icl + 3) = sigma_sigma*DDrot(13,1) - pi_pi*(DDrot(13,2) + DDrot(13,3)) + delta_delta*(DDrot(13,4) + DDrot(13,5));            //<xy|xy>
        overlapdR(irw + 3,icl + 4) = sigma_sigma*DDrot(14,1) - pi_pi*(DDrot(14,2) + DDrot(14,3)) + delta_delta*(DDrot(14,4) + DDrot(14,5));            //<xy|x2-y2>
        overlapdR(irw + 4,icl + 0) = sigma_sigma*DDrot(5,1)  - pi_pi*(DDrot(5,2)  + DDrot(5,3))  + delta_delta*(DDrot(5,4)  + DDrot(5,5));             //<x2-y2|xz>
        overlapdR(irw + 4,icl + 1) = sigma_sigma*DDrot(9,1)  - pi_pi*(DDrot(9,2)  + DDrot(9,3))  + delta_delta*(DDrot(9,4)  + DDrot(9,5));             //<x2-y2|yz>
        overlapdR(irw + 4,icl + 2) = sigma_sigma*DDrot(12,1) - pi_pi*(DDrot(12,2) + DDrot(12,3)) + delta_delta*(DDrot(12,4) + DDrot(12,5));            //<x2-y2|z2>
        overlapdR(irw + 4,icl + 3) = sigma_sigma*DDrot(14,1) - pi_pi*(DDrot(14,2) + DDrot(14,3)) + delta_delta*(DDrot(14,4) + DDrot(14,5));            //<x2-y2|xy>
        overlapdR(irw + 4,icl + 4) = sigma_sigma*DDrot(15,1) - pi_pi*(DDrot(15,2) + DDrot(15,3)) + delta_delta*(DDrot(15,4) + DDrot(15,5));            //<x2-y2|x2-y2>
      }
      else if ((iLA == 3)&&(iLB == 3)) {
        //ff
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = FFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist)/dist;
        delta_delta = 0.0;
        polynomial = FFdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //phi integral
        aux = OvFactor(iNA + 1,iLA,3,a,iNB + 1,iLB,3,b,dist)/dist;
        phi_phi = 0.0;
        polynomial = FFphiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          phi_phi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          phi_phi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          phi_phi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 7x7 dimensions (fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        overlapdR(irw + 0,icl + 0) = -sigma_sigma*FFrot(1,1) + pi_pi*(FFrot(1,2) + FFrot(1,3)) - delta_delta*(FFrot(1,4) + FFrot(1,5)) + phi_phi*(FFrot(1,6) + FFrot(1,7));
        overlapdR(irw + 0,icl + 1) = -sigma_sigma*FFrot(2,1) + pi_pi*(FFrot(2,2) + FFrot(2,3)) - delta_delta*(FFrot(2,4) + FFrot(2,5)) + phi_phi*(FFrot(2,6) + FFrot(2,7));
        overlapdR(irw + 0,icl + 2) = -sigma_sigma*FFrot(3,1) + pi_pi*(FFrot(3,2) + FFrot(3,3)) - delta_delta*(FFrot(3,4) + FFrot(3,5)) + phi_phi*(FFrot(3,6) + FFrot(3,7));
        overlapdR(irw + 0,icl + 3) = -sigma_sigma*FFrot(4,1) + pi_pi*(FFrot(4,2) + FFrot(4,3)) - delta_delta*(FFrot(4,4) + FFrot(4,5)) + phi_phi*(FFrot(4,6) + FFrot(4,7));
        overlapdR(irw + 0,icl + 4) = -sigma_sigma*FFrot(5,1) + pi_pi*(FFrot(5,2) + FFrot(5,3)) - delta_delta*(FFrot(5,4) + FFrot(5,5)) + phi_phi*(FFrot(5,6) + FFrot(5,7));
        overlapdR(irw + 0,icl + 5) = -sigma_sigma*FFrot(6,1) + pi_pi*(FFrot(6,2) + FFrot(6,3)) - delta_delta*(FFrot(6,4) + FFrot(6,5)) + phi_phi*(FFrot(6,6) + FFrot(6,7));
        overlapdR(irw + 0,icl + 6) = -sigma_sigma*FFrot(7,1) + pi_pi*(FFrot(7,2) + FFrot(7,3)) - delta_delta*(FFrot(7,4) + FFrot(7,5)) + phi_phi*(FFrot(7,6) + FFrot(7,7));
        overlapdR(irw + 1,icl + 0) = -sigma_sigma*FFrot(2,1) + pi_pi*(FFrot(2,2) + FFrot(2,3)) - delta_delta*(FFrot(2,4) + FFrot(2,5)) + phi_phi*(FFrot(2,6) + FFrot(2,7));
        overlapdR(irw + 1,icl + 1) = -sigma_sigma*FFrot(8,1) + pi_pi*(FFrot(8,2) + FFrot(8,3)) - delta_delta*(FFrot(8,4) + FFrot(8,5)) + phi_phi*(FFrot(8,6) + FFrot(8,7));
        overlapdR(irw + 1,icl + 2) = -sigma_sigma*FFrot(9,1) + pi_pi*(FFrot(9,2) + FFrot(9,3)) - delta_delta*(FFrot(9,4) + FFrot(9,5)) + phi_phi*(FFrot(9,6) + FFrot(9,7));
        overlapdR(irw + 1,icl + 3) = -sigma_sigma*FFrot(10,1) + pi_pi*(FFrot(10,2) + FFrot(10,3)) - delta_delta*(FFrot(10,4) + FFrot(10,5)) + phi_phi*(FFrot(10,6) + FFrot(10,7));
        overlapdR(irw + 1,icl + 4) = -sigma_sigma*FFrot(11,1) + pi_pi*(FFrot(11,2) + FFrot(11,3)) - delta_delta*(FFrot(11,4) + FFrot(11,5)) + phi_phi*(FFrot(11,6) + FFrot(11,7));
        overlapdR(irw + 1,icl + 5) = -sigma_sigma*FFrot(12,1) + pi_pi*(FFrot(12,2) + FFrot(12,3)) - delta_delta*(FFrot(12,4) + FFrot(12,5)) + phi_phi*(FFrot(12,6) + FFrot(12,7));
        overlapdR(irw + 1,icl + 6) = -sigma_sigma*FFrot(13,1) + pi_pi*(FFrot(13,2) + FFrot(13,3)) - delta_delta*(FFrot(13,4) + FFrot(13,5)) + phi_phi*(FFrot(13,6) + FFrot(13,7));
        overlapdR(irw + 2,icl + 0) = -sigma_sigma*FFrot(3,1) +  pi_pi*(FFrot(3,2) + FFrot(3,3)) - delta_delta*(FFrot(3,4) + FFrot(3,5)) + phi_phi*(FFrot(3,6) + FFrot(3,7));
        overlapdR(irw + 2,icl + 1) = -sigma_sigma*FFrot(9,1) +  pi_pi*(FFrot(9,2) + FFrot(9,3)) - delta_delta*(FFrot(9,4) + FFrot(9,5)) + phi_phi*(FFrot(9,6) + FFrot(9,7));
        overlapdR(irw + 2,icl + 2) = -sigma_sigma*FFrot(14,1) + pi_pi*(FFrot(14,2) + FFrot(14,3)) - delta_delta*(FFrot(14,4) + FFrot(14,5)) + phi_phi*(FFrot(14,6) + FFrot(14,7));
        overlapdR(irw + 2,icl + 3) = -sigma_sigma*FFrot(15,1) + pi_pi*(FFrot(15,2) + FFrot(15,3)) - delta_delta*(FFrot(15,4) + FFrot(15,5)) + phi_phi*(FFrot(15,6) + FFrot(15,7));
        overlapdR(irw + 2,icl + 4) = -sigma_sigma*FFrot(16,1) + pi_pi*(FFrot(16,2) + FFrot(16,3)) - delta_delta*(FFrot(16,4) + FFrot(16,5)) + phi_phi*(FFrot(16,6) + FFrot(16,7));
        overlapdR(irw + 2,icl + 5) = -sigma_sigma*FFrot(17,1) + pi_pi*(FFrot(17,2) + FFrot(17,3)) - delta_delta*(FFrot(17,4) + FFrot(17,5)) + phi_phi*(FFrot(17,6) + FFrot(17,7));
        overlapdR(irw + 2,icl + 6) = -sigma_sigma*FFrot(18,1) + pi_pi*(FFrot(18,2) + FFrot(18,3)) - delta_delta*(FFrot(18,4) + FFrot(18,5)) + phi_phi*(FFrot(18,6) + FFrot(18,7));
        overlapdR(irw + 3,icl + 0) = -sigma_sigma*FFrot(4,1) +  pi_pi*(FFrot(4,2) + FFrot(4,3)) - delta_delta*(FFrot(4,4) + FFrot(4,5)) + phi_phi*(FFrot(4,6) + FFrot(4,7));
        overlapdR(irw + 3,icl + 1) = -sigma_sigma*FFrot(10,1) + pi_pi*(FFrot(10,2) + FFrot(10,3)) - delta_delta*(FFrot(10,4) + FFrot(10,5)) + phi_phi*(FFrot(10,6) + FFrot(10,7));
        overlapdR(irw + 3,icl + 2) = -sigma_sigma*FFrot(15,1) + pi_pi*(FFrot(15,2) + FFrot(15,3)) - delta_delta*(FFrot(15,4) + FFrot(15,5)) + phi_phi*(FFrot(15,6) + FFrot(15,7));
        overlapdR(irw + 3,icl + 3) = -sigma_sigma*FFrot(19,1) + pi_pi*(FFrot(19,2) + FFrot(19,3)) - delta_delta*(FFrot(19,4) + FFrot(19,5)) + phi_phi*(FFrot(19,6) + FFrot(19,7));
        overlapdR(irw + 3,icl + 4) = -sigma_sigma*FFrot(20,1) + pi_pi*(FFrot(20,2) + FFrot(20,3)) - delta_delta*(FFrot(20,4) + FFrot(20,5)) + phi_phi*(FFrot(20,6) + FFrot(20,7));
        overlapdR(irw + 3,icl + 5) = -sigma_sigma*FFrot(21,1) + pi_pi*(FFrot(21,2) + FFrot(21,3)) - delta_delta*(FFrot(21,4) + FFrot(21,5)) + phi_phi*(FFrot(21,6) + FFrot(21,7));
        overlapdR(irw + 3,icl + 6) = -sigma_sigma*FFrot(22,1) + pi_pi*(FFrot(22,2) + FFrot(22,3)) - delta_delta*(FFrot(22,4) + FFrot(22,5)) + phi_phi*(FFrot(22,6) + FFrot(22,7));
        overlapdR(irw + 4,icl + 0) = -sigma_sigma*FFrot(5,1) +  pi_pi*(FFrot(5,2) + FFrot(5,3)) - delta_delta*(FFrot(5,4) + FFrot(5,5)) + phi_phi*(FFrot(5,6) + FFrot(5,7));
        overlapdR(irw + 4,icl + 1) = -sigma_sigma*FFrot(11,1) + pi_pi*(FFrot(11,2) + FFrot(11,3)) - delta_delta*(FFrot(11,4) + FFrot(11,5)) + phi_phi*(FFrot(11,6) + FFrot(11,7));
        overlapdR(irw + 4,icl + 2) = -sigma_sigma*FFrot(16,1) + pi_pi*(FFrot(16,2) + FFrot(16,3)) - delta_delta*(FFrot(16,4) + FFrot(16,5)) + phi_phi*(FFrot(16,6) + FFrot(16,7));
        overlapdR(irw + 4,icl + 3) = -sigma_sigma*FFrot(20,1) + pi_pi*(FFrot(20,2) + FFrot(20,3)) - delta_delta*(FFrot(20,4) + FFrot(20,5)) + phi_phi*(FFrot(20,6) + FFrot(20,7));
        overlapdR(irw + 4,icl + 4) = -sigma_sigma*FFrot(23,1) + pi_pi*(FFrot(23,2) + FFrot(23,3)) - delta_delta*(FFrot(23,4) + FFrot(23,5)) + phi_phi*(FFrot(23,6) + FFrot(23,7));
        overlapdR(irw + 4,icl + 5) = -sigma_sigma*FFrot(24,1) + pi_pi*(FFrot(24,2) + FFrot(24,3)) - delta_delta*(FFrot(24,4) + FFrot(24,5)) + phi_phi*(FFrot(24,6) + FFrot(24,7));
        overlapdR(irw + 4,icl + 6) = -sigma_sigma*FFrot(25,1) + pi_pi*(FFrot(25,2) + FFrot(25,3)) - delta_delta*(FFrot(25,4) + FFrot(25,5)) + phi_phi*(FFrot(25,6) + FFrot(25,7));
        overlapdR(irw + 5,icl + 0) = -sigma_sigma*FFrot(6,1) +  pi_pi*(FFrot(6,2) + FFrot(6,3)) - delta_delta*(FFrot(6,4) + FFrot(6,5)) + phi_phi*(FFrot(6,6) + FFrot(6,7));
        overlapdR(irw + 5,icl + 1) = -sigma_sigma*FFrot(12,1) + pi_pi*(FFrot(12,2) + FFrot(12,3)) - delta_delta*(FFrot(12,4) + FFrot(12,5)) + phi_phi*(FFrot(12,6) + FFrot(12,7));
        overlapdR(irw + 5,icl + 2) = -sigma_sigma*FFrot(17,1) + pi_pi*(FFrot(17,2) + FFrot(17,3)) - delta_delta*(FFrot(17,4) + FFrot(17,5)) + phi_phi*(FFrot(17,6) + FFrot(17,7));
        overlapdR(irw + 5,icl + 3) = -sigma_sigma*FFrot(21,1) + pi_pi*(FFrot(21,2) + FFrot(21,3)) - delta_delta*(FFrot(21,4) + FFrot(21,5)) + phi_phi*(FFrot(21,6) + FFrot(21,7));
        overlapdR(irw + 5,icl + 4) = -sigma_sigma*FFrot(24,1) + pi_pi*(FFrot(24,2) + FFrot(24,3)) - delta_delta*(FFrot(24,4) + FFrot(24,5)) + phi_phi*(FFrot(24,6) + FFrot(24,7));
        overlapdR(irw + 5,icl + 5) = -sigma_sigma*FFrot(26,1) + pi_pi*(FFrot(26,2) + FFrot(26,3)) - delta_delta*(FFrot(26,4) + FFrot(26,5)) + phi_phi*(FFrot(26,6) + FFrot(26,7));
        overlapdR(irw + 5,icl + 6) = -sigma_sigma*FFrot(27,1) + pi_pi*(FFrot(27,2) + FFrot(27,3)) - delta_delta*(FFrot(27,4) + FFrot(27,5)) + phi_phi*(FFrot(27,6) + FFrot(27,7));
        overlapdR(irw + 6,icl + 0) = -sigma_sigma*FFrot(7,1) +  pi_pi*(FFrot(7,2) + FFrot(7,3)) - delta_delta*(FFrot(7,4) + FFrot(7,5)) + phi_phi*(FFrot(7,6) + FFrot(7,7));
        overlapdR(irw + 6,icl + 1) = -sigma_sigma*FFrot(13,1) + pi_pi*(FFrot(13,2) + FFrot(13,3)) - delta_delta*(FFrot(13,4) + FFrot(13,5)) + phi_phi*(FFrot(13,6) + FFrot(13,7));
        overlapdR(irw + 6,icl + 2) = -sigma_sigma*FFrot(18,1) + pi_pi*(FFrot(18,2) + FFrot(18,3)) - delta_delta*(FFrot(18,4) + FFrot(18,5)) + phi_phi*(FFrot(18,6) + FFrot(18,7));
        overlapdR(irw + 6,icl + 3) = -sigma_sigma*FFrot(22,1) + pi_pi*(FFrot(22,2) + FFrot(22,3)) - delta_delta*(FFrot(22,4) + FFrot(22,5)) + phi_phi*(FFrot(22,6) + FFrot(22,7));
        overlapdR(irw + 6,icl + 4) = -sigma_sigma*FFrot(25,1) + pi_pi*(FFrot(25,2) + FFrot(25,3)) - delta_delta*(FFrot(25,4) + FFrot(25,5)) + phi_phi*(FFrot(25,6) + FFrot(25,7));
        overlapdR(irw + 6,icl + 5) = -sigma_sigma*FFrot(27,1) + pi_pi*(FFrot(27,2) + FFrot(27,3)) - delta_delta*(FFrot(27,4) + FFrot(27,5)) + phi_phi*(FFrot(27,6) + FFrot(27,7));
        overlapdR(irw + 6,icl + 6) = -sigma_sigma*FFrot(28,1) + pi_pi*(FFrot(28,2) + FFrot(28,3)) - delta_delta*(FFrot(28,4) + FFrot(28,5)) + phi_phi*(FFrot(28,6) + FFrot(28,7));
      }
      else if ((iLA == 0)&&(iLB == 1)) {
        //sp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SPOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3 dimensions (px,py,pz)
        overlapdR(irw,icl + 0) = -SProt(1,1)*sigma_sigma;
        overlapdR(irw,icl + 1) = -SProt(2,1)*sigma_sigma;
        overlapdR(irw,icl + 2) = -SProt(3,1)*sigma_sigma;
      }
      else if ((iLA == 1)&&(iLB == 0)) {
        //ps
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3 dimensions (px,py,pz)
        overlapdR(irw + 0,icl) = SProt(1,1)*sigma_sigma;
        overlapdR(irw + 1,icl) = SProt(2,1)*sigma_sigma;
        overlapdR(irw + 2,icl) = SProt(3,1)*sigma_sigma;
      }
      else if ((iLA == 0)&&(iLB == 2)) {
        //sd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SDOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          overlapdR(irw,icl + index) = SDrot(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 2)&&(iLB == 0)) {
        //ds
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          overlapdR(irw + index,icl) = SDrot(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 0)&&(iLB == 3)) {
        //sf
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SFOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 7 dimensions (fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 7; ++index) {
          overlapdR(irw,icl + index) = -SFrot(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 3)&&(iLB == 0)) {
        //fs
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 7; ++index) {
          overlapdR(irw + index,icl) = SFrot(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 1)&&(iLB == 2)) {
        //pd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = PDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3x5 dimensions (px,py,pz)(xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 3; ++index2) {
            overlapdR(irw + index2,icl + index) = PDrot(1+index2+3*index,1)*sigma_sigma - (PDrot(1+index2+3*index,2) + PDrot(1+index2+3*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 2)&&(iLB == 1)) {
        //dp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = DPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3x5 dimensions (px,py,pz)(xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 3; ++index2) {
            overlapdR(irw + index,icl + index2) = -PDrot(1+index2+3*index,1)*sigma_sigma + (PDrot(1+index2+3*index,2) + PDrot(1+index2+3*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 1)&&(iLB == 3)) {
        //pf
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = PFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3x5 dimensions (px,py,pz)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 3; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlapdR(irw + index,icl + index2) = -PFrot(1+index2+7*index,1)*sigma_sigma + (PFrot(1+index2+7*index,2) + PFrot(1+index2+7*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 3)&&(iLB == 1)) {
        //fp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = FPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3x5 dimensions (px,py,pz)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 3; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlapdR(irw + index2,icl + index) = -PFrot(1+index2+7*index,1)*sigma_sigma + (PFrot(1+index2+7*index,2) + PFrot(1+index2+7*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 2)&&(iLB == 3)) {
        //df
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = DFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist)/dist;
        delta_delta = 0.0;
        polynomial = DFdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3x5 dimensions (xz,yz,z2,xy,x2-y2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlapdR(irw + index,icl + index2) = -DFrot(1+index2+7*index,1)*sigma_sigma + (DFrot(1+index2+7*index,2) + DFrot(1+index2+7*index,3))*pi_pi - (DFrot(1+index2+7*index,4) + DFrot(1+index2+7*index,5))*delta_delta;
          }
        }
      }
      else if ((iLA == 3)&&(iLB == 2)) {
        //fd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = FDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist)/dist;
        delta_delta = 0.0;
        polynomial = FDdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3x5 dimensions (xz,yz,z2,xy,x2-y2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlapdR(irw + index2,icl + index) = DFrot(1+index2+7*index,1)*sigma_sigma - (DFrot(1+index2+7*index,2) + DFrot(1+index2+7*index,3))*pi_pi + (DFrot(1+index2+7*index,4) + DFrot(1+index2+7*index,5))*delta_delta;
          }
        }
      }
      icl += 2*QN2(isket + 1,2) + 1;
    }
    irw += 2*QN1(isbra + 1,2) + 1;
  }
}
void TwoCenterOv_dt(matrixE & overlaps, size_t nao1, size_t nao2, double dist, std::vector<double> & exponents1, matrix<size_t> & QN1, std::vector<double> & exponents2, matrix<size_t> & QN2, std::vector<double> & orientation, bool rotate = true, double tolerance = 1e-6) {
  //implementation of first-derivative of overlap matrix function with respect to theta; only up to shell f
  //calculation of the derivative of the overlap function between atoms A and B (A -> 1; B -> 2) with respect to angle theta
  overlaps.resize(nao1,nao2);
  overlaps.zero();
  //auxiliary variables for integrals
  double sigma_sigma = 0.0;
  double pi_pi = 0.0;
  double delta_delta = 0.0;
  double phi_phi = 0.0;
  //angles needed to reorient the integrals according to the position of atoms in the molecule
  double cost = orientation[2];
  double sint = sqrt(1 - cost*cost);
  double cosp = 1.0;
  double sinp = 0.0;
  if (fabs(sint) > tolerance) {
    cosp = orientation[0]/sint;
    sinp = orientation[1]/sint;
  }
  //matrices with derivatives of overlap rotations with respect to theta
  matrixE SProt_dt(1,1);
  matrixE SDrot_dt(1,1);
  matrixE SFrot_dt(1,1);
  matrixE PProt_dt(1,1);
  matrixE PDrot_dt(1,1);
  matrixE PFrot_dt(1,1);
  matrixE DDrot_dt(1,1);
  matrixE DFrot_dt(1,1);
  matrixE FFrot_dt(1,1);
  //other auxiliaries
  std::vector<std::vector<int> > polynomial;
  double a;
  double b;
  double c;
  double argA;
  double argB;
  std::vector<double> Ak;
  std::vector<double> Bk;
  double aux;
  //determine maximum azimuthal quantum number
  size_t maxL1 = 0;
  size_t maxL2 = 0;
  for (size_t idx = 0; idx < QN1.rows(); ++idx) {
    if (QN1(idx + 1,2) > maxL1) {maxL1 = QN1(idx + 1,2);}
  }
  for (size_t idx = 0; idx < QN2.rows(); ++idx) {
    if (QN2(idx + 1,2) > maxL2) {maxL2 = QN2(idx + 1,2);}
  }
  if ((maxL1 > 0)||(maxL2 > 0)) {
    //in this case we have p orbitals
    SProt_dt = SPtransf_dt(cost,sint,cosp,sinp,rotate);
    PProt_dt = PPtransf_dt(cost,sint,cosp,sinp,rotate);
  }
  if ((maxL1 > 1)||(maxL2 > 1)) {
    //in this case we have d orbitals
    SDrot_dt = SDtransf_dt(cost,sint,cosp,sinp,rotate);
    PDrot_dt = PDtransf_dt(cost,sint,cosp,sinp,rotate);
    DDrot_dt = DDtransf_dt(cost,sint,cosp,sinp,rotate);
  }
  if ((maxL1 > 2)||(maxL2 > 2)) {
    //in this case we have f orbitals
    SFrot_dt = SFtransf_dt(cost,sint,cosp,sinp,rotate);
    PFrot_dt = PFtransf_dt(cost,sint,cosp,sinp,rotate);
    DFrot_dt = DFtransf_dt(cost,sint,cosp,sinp,rotate);
    FFrot_dt = FFtransf_dt(cost,sint,cosp,sinp,rotate);
  }
  size_t iNA;
  size_t iNB;
  size_t iLA;
  size_t iLB;
  size_t irw = 1;
  size_t icl = 1;
  for (size_t isbra = 0; isbra < QN1.rows(); ++isbra) {           //loop over bra shells
    a = exponents1[isbra];
    if (a == 0.0) continue;
    iNA = QN1(isbra + 1,1) - 1;
    iLA = QN1(isbra + 1,2);
    icl = 1;
    for (size_t isket = 0; isket < QN2.rows(); ++isket) {         //loop over ket shells
      b = exponents2[isket];
      if (b == 0.0) continue;
      iNB = QN2(isket + 1,1) - 1;
      iLB = QN2(isket + 1,2);
      //getting integrals
      argA = 0.5*dist*(a + b);
      argB = 0.5*dist*(a - b);
      Ak = int_1_infty_xpk_expmyx_dx(argA,iNA + iNB + 2);
      Bk = int_m1_1_xpk_expmyx_dx_Pople(argB,iNA + iNB + 2);
      aux = OvFactor(iNA + 1,iLA,0,a,iNB + 1,iLB,0,b,dist);
      //different azimuthal quantum number cases
      if ((iLA == 0)&&(iLB == 0)) {
        //ss
        overlaps(irw,icl) = 0.0;
      }
      //ss = 0
      else if ((iLA == 1)&&(iLB == 1)) {
        //pp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = PPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 9 dimensions
        overlaps(irw + 0,icl + 0) = -sigma_sigma*PProt_dt(1,1) + pi_pi*(PProt_dt(1,2) + PProt_dt(1,3));      //<x|x>
        overlaps(irw + 0,icl + 1) = -sigma_sigma*PProt_dt(2,1) + pi_pi*(PProt_dt(2,2) + PProt_dt(2,3));      //<x|y>
        overlaps(irw + 0,icl + 2) = -sigma_sigma*PProt_dt(3,1) + pi_pi*(PProt_dt(3,2) + PProt_dt(3,3));      //<x|z>
        overlaps(irw + 1,icl + 0) = -sigma_sigma*PProt_dt(2,1) + pi_pi*(PProt_dt(2,2) + PProt_dt(2,3));      //<y|x>
        overlaps(irw + 1,icl + 1) = -sigma_sigma*PProt_dt(4,1) + pi_pi*(PProt_dt(4,2) + PProt_dt(4,3));      //<y|y>
        overlaps(irw + 1,icl + 2) = -sigma_sigma*PProt_dt(5,1) + pi_pi*(PProt_dt(5,2) + PProt_dt(5,3));      //<y|z>
        overlaps(irw + 2,icl + 0) = -sigma_sigma*PProt_dt(3,1) + pi_pi*(PProt_dt(3,2) + PProt_dt(3,3));      //<z|x>
        overlaps(irw + 2,icl + 1) = -sigma_sigma*PProt_dt(5,1) + pi_pi*(PProt_dt(5,2) + PProt_dt(5,3));      //<z|y>
        overlaps(irw + 2,icl + 2) = -sigma_sigma*PProt_dt(6,1) + pi_pi*(PProt_dt(6,2) + PProt_dt(6,3));      //<z|z>
      }
      else if ((iLA == 2)&&(iLB == 2)) {
        //dd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = DDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = DDdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5x5 dimensions (xz,yz,z2,xy,x2-y2)(xz,yz,z2,xy,x2-y2)
        overlaps(irw + 0,icl + 0) = sigma_sigma*DDrot_dt(1,1)  - pi_pi*(DDrot_dt(1,2)  + DDrot_dt(1,3))  + delta_delta*(DDrot_dt(1,4)  + DDrot_dt(1,5));             //<xz|xz>
        overlaps(irw + 0,icl + 1) = sigma_sigma*DDrot_dt(2,1)  - pi_pi*(DDrot_dt(2,2)  + DDrot_dt(2,3))  + delta_delta*(DDrot_dt(2,4)  + DDrot_dt(2,5));             //<xz|yz>
        overlaps(irw + 0,icl + 2) = sigma_sigma*DDrot_dt(3,1)  - pi_pi*(DDrot_dt(3,2)  + DDrot_dt(3,3))  + delta_delta*(DDrot_dt(3,4)  + DDrot_dt(3,5));             //<xz|z2>
        overlaps(irw + 0,icl + 3) = sigma_sigma*DDrot_dt(4,1)  - pi_pi*(DDrot_dt(4,2)  + DDrot_dt(4,3))  + delta_delta*(DDrot_dt(4,4)  + DDrot_dt(4,5));             //<xz|xy>
        overlaps(irw + 0,icl + 4) = sigma_sigma*DDrot_dt(5,1)  - pi_pi*(DDrot_dt(5,2)  + DDrot_dt(5,3))  + delta_delta*(DDrot_dt(5,4)  + DDrot_dt(5,5));             //<xz|x2-y2>
        overlaps(irw + 1,icl + 0) = sigma_sigma*DDrot_dt(2,1)  - pi_pi*(DDrot_dt(2,2)  + DDrot_dt(2,3))  + delta_delta*(DDrot_dt(2,4)  + DDrot_dt(2,5));             //<yz|xz>
        overlaps(irw + 1,icl + 1) = sigma_sigma*DDrot_dt(6,1)  - pi_pi*(DDrot_dt(6,2)  + DDrot_dt(6,3))  + delta_delta*(DDrot_dt(6,4)  + DDrot_dt(6,5));             //<yz|yz>
        overlaps(irw + 1,icl + 2) = sigma_sigma*DDrot_dt(7,1)  - pi_pi*(DDrot_dt(7,2)  + DDrot_dt(7,3))  + delta_delta*(DDrot_dt(7,4)  + DDrot_dt(7,5));             //<yz|z2>
        overlaps(irw + 1,icl + 3) = sigma_sigma*DDrot_dt(8,1)  - pi_pi*(DDrot_dt(8,2)  + DDrot_dt(8,3))  + delta_delta*(DDrot_dt(8,4)  + DDrot_dt(8,5));             //<yz|xy>
        overlaps(irw + 1,icl + 4) = sigma_sigma*DDrot_dt(9,1)  - pi_pi*(DDrot_dt(9,2)  + DDrot_dt(9,3))  + delta_delta*(DDrot_dt(9,4)  + DDrot_dt(9,5));             //<yz|x2-y2>
        overlaps(irw + 2,icl + 0) = sigma_sigma*DDrot_dt(3,1)  - pi_pi*(DDrot_dt(3,2)  + DDrot_dt(3,3))  + delta_delta*(DDrot_dt(3,4)  + DDrot_dt(3,5));             //<z2|xz>
        overlaps(irw + 2,icl + 1) = sigma_sigma*DDrot_dt(7,1)  - pi_pi*(DDrot_dt(7,2)  + DDrot_dt(7,3))  + delta_delta*(DDrot_dt(7,4)  + DDrot_dt(7,5));             //<z2|yz>
        overlaps(irw + 2,icl + 2) = sigma_sigma*DDrot_dt(10,1) - pi_pi*(DDrot_dt(10,2) + DDrot_dt(10,3)) + delta_delta*(DDrot_dt(10,4) + DDrot_dt(10,5));            //<z2|z2>
        overlaps(irw + 2,icl + 3) = sigma_sigma*DDrot_dt(11,1) - pi_pi*(DDrot_dt(11,2) + DDrot_dt(11,3)) + delta_delta*(DDrot_dt(11,4) + DDrot_dt(11,5));            //<z2|xy>
        overlaps(irw + 2,icl + 4) = sigma_sigma*DDrot_dt(12,1) - pi_pi*(DDrot_dt(12,2) + DDrot_dt(12,3)) + delta_delta*(DDrot_dt(12,4) + DDrot_dt(12,5));            //<z2|x2-y2>
        overlaps(irw + 3,icl + 0) = sigma_sigma*DDrot_dt(4,1)  - pi_pi*(DDrot_dt(4,2)  + DDrot_dt(4,3))  + delta_delta*(DDrot_dt(4,4)  + DDrot_dt(4,5));             //<xy|xz>
        overlaps(irw + 3,icl + 1) = sigma_sigma*DDrot_dt(8,1)  - pi_pi*(DDrot_dt(8,2)  + DDrot_dt(8,3))  + delta_delta*(DDrot_dt(8,4)  + DDrot_dt(8,5));             //<xy|yz>
        overlaps(irw + 3,icl + 2) = sigma_sigma*DDrot_dt(11,1) - pi_pi*(DDrot_dt(11,2) + DDrot_dt(11,3)) + delta_delta*(DDrot_dt(11,4) + DDrot_dt(11,5));            //<xy|z2>
        overlaps(irw + 3,icl + 3) = sigma_sigma*DDrot_dt(13,1) - pi_pi*(DDrot_dt(13,2) + DDrot_dt(13,3)) + delta_delta*(DDrot_dt(13,4) + DDrot_dt(13,5));            //<xy|xy>
        overlaps(irw + 3,icl + 4) = sigma_sigma*DDrot_dt(14,1) - pi_pi*(DDrot_dt(14,2) + DDrot_dt(14,3)) + delta_delta*(DDrot_dt(14,4) + DDrot_dt(14,5));            //<xy|x2-y2>
        overlaps(irw + 4,icl + 0) = sigma_sigma*DDrot_dt(5,1)  - pi_pi*(DDrot_dt(5,2)  + DDrot_dt(5,3))  + delta_delta*(DDrot_dt(5,4)  + DDrot_dt(5,5));             //<x2-y2|xz>
        overlaps(irw + 4,icl + 1) = sigma_sigma*DDrot_dt(9,1)  - pi_pi*(DDrot_dt(9,2)  + DDrot_dt(9,3))  + delta_delta*(DDrot_dt(9,4)  + DDrot_dt(9,5));             //<x2-y2|yz>
        overlaps(irw + 4,icl + 2) = sigma_sigma*DDrot_dt(12,1) - pi_pi*(DDrot_dt(12,2) + DDrot_dt(12,3)) + delta_delta*(DDrot_dt(12,4) + DDrot_dt(12,5));            //<x2-y2|z2>
        overlaps(irw + 4,icl + 3) = sigma_sigma*DDrot_dt(14,1) - pi_pi*(DDrot_dt(14,2) + DDrot_dt(14,3)) + delta_delta*(DDrot_dt(14,4) + DDrot_dt(14,5));            //<x2-y2|xy>
        overlaps(irw + 4,icl + 4) = sigma_sigma*DDrot_dt(15,1) - pi_pi*(DDrot_dt(15,2) + DDrot_dt(15,3)) + delta_delta*(DDrot_dt(15,4) + DDrot_dt(15,5));            //<x2-y2|x2-y2>
      }
      else if ((iLA == 3)&&(iLB == 3)) {
        //ff
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = FFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = FFdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //phi integral
        aux = OvFactor(iNA + 1,iLA,3,a,iNB + 1,iLB,3,b,dist);
        phi_phi = 0.0;
        polynomial = FFphiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          phi_phi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 7x7 dimensions (fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        overlaps(irw + 0,icl + 0) = -sigma_sigma*FFrot_dt(1,1)  + pi_pi*(FFrot_dt(1,2)  + FFrot_dt(1,3))  - delta_delta*(FFrot_dt(1,4)  + FFrot_dt(1,5))  + phi_phi*(FFrot_dt(1,6)  + FFrot_dt(1,7));
        overlaps(irw + 0,icl + 1) = -sigma_sigma*FFrot_dt(2,1)  + pi_pi*(FFrot_dt(2,2)  + FFrot_dt(2,3))  - delta_delta*(FFrot_dt(2,4)  + FFrot_dt(2,5))  + phi_phi*(FFrot_dt(2,6)  + FFrot_dt(2,7));
        overlaps(irw + 0,icl + 2) = -sigma_sigma*FFrot_dt(3,1)  + pi_pi*(FFrot_dt(3,2)  + FFrot_dt(3,3))  - delta_delta*(FFrot_dt(3,4)  + FFrot_dt(3,5))  + phi_phi*(FFrot_dt(3,6)  + FFrot_dt(3,7));
        overlaps(irw + 0,icl + 3) = -sigma_sigma*FFrot_dt(4,1)  + pi_pi*(FFrot_dt(4,2)  + FFrot_dt(4,3))  - delta_delta*(FFrot_dt(4,4)  + FFrot_dt(4,5))  + phi_phi*(FFrot_dt(4,6)  + FFrot_dt(4,7));
        overlaps(irw + 0,icl + 4) = -sigma_sigma*FFrot_dt(5,1)  + pi_pi*(FFrot_dt(5,2)  + FFrot_dt(5,3))  - delta_delta*(FFrot_dt(5,4)  + FFrot_dt(5,5))  + phi_phi*(FFrot_dt(5,6)  + FFrot_dt(5,7));
        overlaps(irw + 0,icl + 5) = -sigma_sigma*FFrot_dt(6,1)  + pi_pi*(FFrot_dt(6,2)  + FFrot_dt(6,3))  - delta_delta*(FFrot_dt(6,4)  + FFrot_dt(6,5))  + phi_phi*(FFrot_dt(6,6)  + FFrot_dt(6,7));
        overlaps(irw + 0,icl + 6) = -sigma_sigma*FFrot_dt(7,1)  + pi_pi*(FFrot_dt(7,2)  + FFrot_dt(7,3))  - delta_delta*(FFrot_dt(7,4)  + FFrot_dt(7,5))  + phi_phi*(FFrot_dt(7,6)  + FFrot_dt(7,7));
        overlaps(irw + 1,icl + 0) = -sigma_sigma*FFrot_dt(2,1)  + pi_pi*(FFrot_dt(2,2)  + FFrot_dt(2,3))  - delta_delta*(FFrot_dt(2,4)  + FFrot_dt(2,5))  + phi_phi*(FFrot_dt(2,6)  + FFrot_dt(2,7));
        overlaps(irw + 1,icl + 1) = -sigma_sigma*FFrot_dt(8,1)  + pi_pi*(FFrot_dt(8,2)  + FFrot_dt(8,3))  - delta_delta*(FFrot_dt(8,4)  + FFrot_dt(8,5))  + phi_phi*(FFrot_dt(8,6)  + FFrot_dt(8,7));
        overlaps(irw + 1,icl + 2) = -sigma_sigma*FFrot_dt(9,1)  + pi_pi*(FFrot_dt(9,2)  + FFrot_dt(9,3))  - delta_delta*(FFrot_dt(9,4)  + FFrot_dt(9,5))  + phi_phi*(FFrot_dt(9,6)  + FFrot_dt(9,7));
        overlaps(irw + 1,icl + 3) = -sigma_sigma*FFrot_dt(10,1) + pi_pi*(FFrot_dt(10,2) + FFrot_dt(10,3)) - delta_delta*(FFrot_dt(10,4) + FFrot_dt(10,5)) + phi_phi*(FFrot_dt(10,6) + FFrot_dt(10,7));
        overlaps(irw + 1,icl + 4) = -sigma_sigma*FFrot_dt(11,1) + pi_pi*(FFrot_dt(11,2) + FFrot_dt(11,3)) - delta_delta*(FFrot_dt(11,4) + FFrot_dt(11,5)) + phi_phi*(FFrot_dt(11,6) + FFrot_dt(11,7));
        overlaps(irw + 1,icl + 5) = -sigma_sigma*FFrot_dt(12,1) + pi_pi*(FFrot_dt(12,2) + FFrot_dt(12,3)) - delta_delta*(FFrot_dt(12,4) + FFrot_dt(12,5)) + phi_phi*(FFrot_dt(12,6) + FFrot_dt(12,7));
        overlaps(irw + 1,icl + 6) = -sigma_sigma*FFrot_dt(13,1) + pi_pi*(FFrot_dt(13,2) + FFrot_dt(13,3)) - delta_delta*(FFrot_dt(13,4) + FFrot_dt(13,5)) + phi_phi*(FFrot_dt(13,6) + FFrot_dt(13,7));
        overlaps(irw + 2,icl + 0) = -sigma_sigma*FFrot_dt(3,1) +  pi_pi*(FFrot_dt(3,2)  + FFrot_dt(3,3))  - delta_delta*(FFrot_dt(3,4)  + FFrot_dt(3,5))  + phi_phi*(FFrot_dt(3,6)  + FFrot_dt(3,7));
        overlaps(irw + 2,icl + 1) = -sigma_sigma*FFrot_dt(9,1) +  pi_pi*(FFrot_dt(9,2)  + FFrot_dt(9,3))  - delta_delta*(FFrot_dt(9,4)  + FFrot_dt(9,5))  + phi_phi*(FFrot_dt(9,6)  + FFrot_dt(9,7));
        overlaps(irw + 2,icl + 2) = -sigma_sigma*FFrot_dt(14,1) + pi_pi*(FFrot_dt(14,2) + FFrot_dt(14,3)) - delta_delta*(FFrot_dt(14,4) + FFrot_dt(14,5)) + phi_phi*(FFrot_dt(14,6) + FFrot_dt(14,7));
        overlaps(irw + 2,icl + 3) = -sigma_sigma*FFrot_dt(15,1) + pi_pi*(FFrot_dt(15,2) + FFrot_dt(15,3)) - delta_delta*(FFrot_dt(15,4) + FFrot_dt(15,5)) + phi_phi*(FFrot_dt(15,6) + FFrot_dt(15,7));
        overlaps(irw + 2,icl + 4) = -sigma_sigma*FFrot_dt(16,1) + pi_pi*(FFrot_dt(16,2) + FFrot_dt(16,3)) - delta_delta*(FFrot_dt(16,4) + FFrot_dt(16,5)) + phi_phi*(FFrot_dt(16,6) + FFrot_dt(16,7));
        overlaps(irw + 2,icl + 5) = -sigma_sigma*FFrot_dt(17,1) + pi_pi*(FFrot_dt(17,2) + FFrot_dt(17,3)) - delta_delta*(FFrot_dt(17,4) + FFrot_dt(17,5)) + phi_phi*(FFrot_dt(17,6) + FFrot_dt(17,7));
        overlaps(irw + 2,icl + 6) = -sigma_sigma*FFrot_dt(18,1) + pi_pi*(FFrot_dt(18,2) + FFrot_dt(18,3)) - delta_delta*(FFrot_dt(18,4) + FFrot_dt(18,5)) + phi_phi*(FFrot_dt(18,6) + FFrot_dt(18,7));
        overlaps(irw + 3,icl + 0) = -sigma_sigma*FFrot_dt(4,1) +  pi_pi*(FFrot_dt(4,2)  + FFrot_dt(4,3))  - delta_delta*(FFrot_dt(4,4)  + FFrot_dt(4,5))  + phi_phi*(FFrot_dt(4,6)  + FFrot_dt(4,7));
        overlaps(irw + 3,icl + 1) = -sigma_sigma*FFrot_dt(10,1) + pi_pi*(FFrot_dt(10,2) + FFrot_dt(10,3)) - delta_delta*(FFrot_dt(10,4) + FFrot_dt(10,5)) + phi_phi*(FFrot_dt(10,6) + FFrot_dt(10,7));
        overlaps(irw + 3,icl + 2) = -sigma_sigma*FFrot_dt(15,1) + pi_pi*(FFrot_dt(15,2) + FFrot_dt(15,3)) - delta_delta*(FFrot_dt(15,4) + FFrot_dt(15,5)) + phi_phi*(FFrot_dt(15,6) + FFrot_dt(15,7));
        overlaps(irw + 3,icl + 3) = -sigma_sigma*FFrot_dt(19,1) + pi_pi*(FFrot_dt(19,2) + FFrot_dt(19,3)) - delta_delta*(FFrot_dt(19,4) + FFrot_dt(19,5)) + phi_phi*(FFrot_dt(19,6) + FFrot_dt(19,7));
        overlaps(irw + 3,icl + 4) = -sigma_sigma*FFrot_dt(20,1) + pi_pi*(FFrot_dt(20,2) + FFrot_dt(20,3)) - delta_delta*(FFrot_dt(20,4) + FFrot_dt(20,5)) + phi_phi*(FFrot_dt(20,6) + FFrot_dt(20,7));
        overlaps(irw + 3,icl + 5) = -sigma_sigma*FFrot_dt(21,1) + pi_pi*(FFrot_dt(21,2) + FFrot_dt(21,3)) - delta_delta*(FFrot_dt(21,4) + FFrot_dt(21,5)) + phi_phi*(FFrot_dt(21,6) + FFrot_dt(21,7));
        overlaps(irw + 3,icl + 6) = -sigma_sigma*FFrot_dt(22,1) + pi_pi*(FFrot_dt(22,2) + FFrot_dt(22,3)) - delta_delta*(FFrot_dt(22,4) + FFrot_dt(22,5)) + phi_phi*(FFrot_dt(22,6) + FFrot_dt(22,7));
        overlaps(irw + 4,icl + 0) = -sigma_sigma*FFrot_dt(5,1) +  pi_pi*(FFrot_dt(5,2)  + FFrot_dt(5,3))  - delta_delta*(FFrot_dt(5,4)  + FFrot_dt(5,5))  + phi_phi*(FFrot_dt(5,6)  + FFrot_dt(5,7));
        overlaps(irw + 4,icl + 1) = -sigma_sigma*FFrot_dt(11,1) + pi_pi*(FFrot_dt(11,2) + FFrot_dt(11,3)) - delta_delta*(FFrot_dt(11,4) + FFrot_dt(11,5)) + phi_phi*(FFrot_dt(11,6) + FFrot_dt(11,7));
        overlaps(irw + 4,icl + 2) = -sigma_sigma*FFrot_dt(16,1) + pi_pi*(FFrot_dt(16,2) + FFrot_dt(16,3)) - delta_delta*(FFrot_dt(16,4) + FFrot_dt(16,5)) + phi_phi*(FFrot_dt(16,6) + FFrot_dt(16,7));
        overlaps(irw + 4,icl + 3) = -sigma_sigma*FFrot_dt(20,1) + pi_pi*(FFrot_dt(20,2) + FFrot_dt(20,3)) - delta_delta*(FFrot_dt(20,4) + FFrot_dt(20,5)) + phi_phi*(FFrot_dt(20,6) + FFrot_dt(20,7));
        overlaps(irw + 4,icl + 4) = -sigma_sigma*FFrot_dt(23,1) + pi_pi*(FFrot_dt(23,2) + FFrot_dt(23,3)) - delta_delta*(FFrot_dt(23,4) + FFrot_dt(23,5)) + phi_phi*(FFrot_dt(23,6) + FFrot_dt(23,7));
        overlaps(irw + 4,icl + 5) = -sigma_sigma*FFrot_dt(24,1) + pi_pi*(FFrot_dt(24,2) + FFrot_dt(24,3)) - delta_delta*(FFrot_dt(24,4) + FFrot_dt(24,5)) + phi_phi*(FFrot_dt(24,6) + FFrot_dt(24,7));
        overlaps(irw + 4,icl + 6) = -sigma_sigma*FFrot_dt(25,1) + pi_pi*(FFrot_dt(25,2) + FFrot_dt(25,3)) - delta_delta*(FFrot_dt(25,4) + FFrot_dt(25,5)) + phi_phi*(FFrot_dt(25,6) + FFrot_dt(25,7));
        overlaps(irw + 5,icl + 0) = -sigma_sigma*FFrot_dt(6,1) +  pi_pi*(FFrot_dt(6,2)  + FFrot_dt(6,3))  - delta_delta*(FFrot_dt(6,4)  + FFrot_dt(6,5))  + phi_phi*(FFrot_dt(6,6)  + FFrot_dt(6,7));
        overlaps(irw + 5,icl + 1) = -sigma_sigma*FFrot_dt(12,1) + pi_pi*(FFrot_dt(12,2) + FFrot_dt(12,3)) - delta_delta*(FFrot_dt(12,4) + FFrot_dt(12,5)) + phi_phi*(FFrot_dt(12,6) + FFrot_dt(12,7));
        overlaps(irw + 5,icl + 2) = -sigma_sigma*FFrot_dt(17,1) + pi_pi*(FFrot_dt(17,2) + FFrot_dt(17,3)) - delta_delta*(FFrot_dt(17,4) + FFrot_dt(17,5)) + phi_phi*(FFrot_dt(17,6) + FFrot_dt(17,7));
        overlaps(irw + 5,icl + 3) = -sigma_sigma*FFrot_dt(21,1) + pi_pi*(FFrot_dt(21,2) + FFrot_dt(21,3)) - delta_delta*(FFrot_dt(21,4) + FFrot_dt(21,5)) + phi_phi*(FFrot_dt(21,6) + FFrot_dt(21,7));
        overlaps(irw + 5,icl + 4) = -sigma_sigma*FFrot_dt(24,1) + pi_pi*(FFrot_dt(24,2) + FFrot_dt(24,3)) - delta_delta*(FFrot_dt(24,4) + FFrot_dt(24,5)) + phi_phi*(FFrot_dt(24,6) + FFrot_dt(24,7));
        overlaps(irw + 5,icl + 5) = -sigma_sigma*FFrot_dt(26,1) + pi_pi*(FFrot_dt(26,2) + FFrot_dt(26,3)) - delta_delta*(FFrot_dt(26,4) + FFrot_dt(26,5)) + phi_phi*(FFrot_dt(26,6) + FFrot_dt(26,7));
        overlaps(irw + 5,icl + 6) = -sigma_sigma*FFrot_dt(27,1) + pi_pi*(FFrot_dt(27,2) + FFrot_dt(27,3)) - delta_delta*(FFrot_dt(27,4) + FFrot_dt(27,5)) + phi_phi*(FFrot_dt(27,6) + FFrot_dt(27,7));
        overlaps(irw + 6,icl + 0) = -sigma_sigma*FFrot_dt(7,1) +  pi_pi*(FFrot_dt(7,2)  + FFrot_dt(7,3))  - delta_delta*(FFrot_dt(7,4)  + FFrot_dt(7,5))  + phi_phi*(FFrot_dt(7,6)  + FFrot_dt(7,7));
        overlaps(irw + 6,icl + 1) = -sigma_sigma*FFrot_dt(13,1) + pi_pi*(FFrot_dt(13,2) + FFrot_dt(13,3)) - delta_delta*(FFrot_dt(13,4) + FFrot_dt(13,5)) + phi_phi*(FFrot_dt(13,6) + FFrot_dt(13,7));
        overlaps(irw + 6,icl + 2) = -sigma_sigma*FFrot_dt(18,1) + pi_pi*(FFrot_dt(18,2) + FFrot_dt(18,3)) - delta_delta*(FFrot_dt(18,4) + FFrot_dt(18,5)) + phi_phi*(FFrot_dt(18,6) + FFrot_dt(18,7));
        overlaps(irw + 6,icl + 3) = -sigma_sigma*FFrot_dt(22,1) + pi_pi*(FFrot_dt(22,2) + FFrot_dt(22,3)) - delta_delta*(FFrot_dt(22,4) + FFrot_dt(22,5)) + phi_phi*(FFrot_dt(22,6) + FFrot_dt(22,7));
        overlaps(irw + 6,icl + 4) = -sigma_sigma*FFrot_dt(25,1) + pi_pi*(FFrot_dt(25,2) + FFrot_dt(25,3)) - delta_delta*(FFrot_dt(25,4) + FFrot_dt(25,5)) + phi_phi*(FFrot_dt(25,6) + FFrot_dt(25,7));
        overlaps(irw + 6,icl + 5) = -sigma_sigma*FFrot_dt(27,1) + pi_pi*(FFrot_dt(27,2) + FFrot_dt(27,3)) - delta_delta*(FFrot_dt(27,4) + FFrot_dt(27,5)) + phi_phi*(FFrot_dt(27,6) + FFrot_dt(27,7));
        overlaps(irw + 6,icl + 6) = -sigma_sigma*FFrot_dt(28,1) + pi_pi*(FFrot_dt(28,2) + FFrot_dt(28,3)) - delta_delta*(FFrot_dt(28,4) + FFrot_dt(28,5)) + phi_phi*(FFrot_dt(28,6) + FFrot_dt(28,7));
      }
      else if ((iLA == 0)&&(iLB == 1)) {
        //sp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SPOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3 dimensions (px,py,pz)
        overlaps(irw,icl + 0) = -SProt_dt(1,1)*sigma_sigma;
        overlaps(irw,icl + 1) = -SProt_dt(2,1)*sigma_sigma;
        overlaps(irw,icl + 2) = -SProt_dt(3,1)*sigma_sigma;
      }
      else if ((iLA == 1)&&(iLB == 0)) {
        //ps
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3 dimensions (px,py,pz)
        overlaps(irw + 0,icl) = SProt_dt(1,1)*sigma_sigma;
        overlaps(irw + 1,icl) = SProt_dt(2,1)*sigma_sigma;
        overlaps(irw + 2,icl) = SProt_dt(3,1)*sigma_sigma;
      }
      else if ((iLA == 0)&&(iLB == 2)) {
        //sd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SDOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          overlaps(irw,icl + index) = SDrot_dt(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 2)&&(iLB == 0)) {
        //ds
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          overlaps(irw + index,icl) = SDrot_dt(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 0)&&(iLB == 3)) {
        //sf
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SFOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 7 dimensions (fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 7; ++index) {
          overlaps(irw,icl + index) = -SFrot_dt(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 3)&&(iLB == 0)) {
        //fs
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 7; ++index) {
          overlaps(irw + index,icl) = SFrot_dt(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 1)&&(iLB == 2)) {
        //pd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = PDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 3; ++index2) {
            overlaps(irw + index2,icl + index) = PDrot_dt(1+index2+3*index,1)*sigma_sigma - (PDrot_dt(1+index2+3*index,2) + PDrot_dt(1+index2+3*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 2)&&(iLB == 1)) {
        //dp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = DPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 3; ++index2) {
            overlaps(irw + index,icl + index2) = -PDrot_dt(1+index2+3*index,1)*sigma_sigma + (PDrot_dt(1+index2+3*index,2) + PDrot_dt(1+index2+3*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 1)&&(iLB == 3)) {
        //pf
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = PFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 3; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index,icl + index2) = -PFrot_dt(1+index2+7*index,1)*sigma_sigma + (PFrot_dt(1+index2+7*index,2) + PFrot_dt(1+index2+7*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 3)&&(iLB == 1)) {
        //fp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = FPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 3; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index2,icl + index) = -PFrot_dt(1+index2+7*index,1)*sigma_sigma + (PFrot_dt(1+index2+7*index,2) + PFrot_dt(1+index2+7*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 2)&&(iLB == 3)) {
        //df
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = DFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = DFdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (xz,yz,z2,xy,x2-y2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index,icl + index2) = -DFrot_dt(1+index2+7*index,1)*sigma_sigma + (DFrot_dt(1+index2+7*index,2) + DFrot_dt(1+index2+7*index,3))*pi_pi - (DFrot_dt(1+index2+7*index,4) + DFrot_dt(1+index2+7*index,5))*delta_delta;
          }
        }
      }
      else if ((iLA == 3)&&(iLB == 2)) {
        //fd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = FDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = FDdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (xz,yz,z2,xy,x2-y2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index2,icl + index) = DFrot_dt(1+index2+7*index,1)*sigma_sigma - (DFrot_dt(1+index2+7*index,2) + DFrot_dt(1+index2+7*index,3))*pi_pi + (DFrot_dt(1+index2+7*index,4) + DFrot_dt(1+index2+7*index,5))*delta_delta;
          }
        }
      }
      icl += 2*QN2(isket + 1,2) + 1;
    }
    irw += 2*QN1(isbra + 1,2) + 1;
  }
}
void TwoCenterOv_dp(matrixE & overlaps, size_t nao1, size_t nao2, double dist, std::vector<double> & exponents1, matrix<size_t> & QN1, std::vector<double> & exponents2, matrix<size_t> & QN2, std::vector<double> & orientation, bool rotate = true, double tolerance = 1e-6) {
  //implementation of first-derivative of overlap matrix function with respect to phi; only up to shell f
  //calculation of the derivative of the overlap function between atoms A and B (A -> 1; B -> 2) with respect to angle phi
  overlaps.resize(nao1,nao2);
  overlaps.zero();
  //auxiliary variables for integrals
  double sigma_sigma = 0.0;
  double pi_pi = 0.0;
  double delta_delta = 0.0;
  double phi_phi = 0.0;
  //angles needed to reorient the integrals according to the position of atoms in the molecule
  double cost = orientation[2];
  double sint = sqrt(1 - cost*cost);
  double cosp = 1.0;
  double sinp = 0.0;
  if (fabs(sint) > tolerance) {
    cosp = orientation[0]/sint;
    sinp = orientation[1]/sint;
  }
  //matrices with derivatives of overlap rotations with respect to phi
  matrixE SProt_dp(1,1);
  matrixE SDrot_dp(1,1);
  matrixE SFrot_dp(1,1);
  matrixE PProt_dp(1,1);
  matrixE PDrot_dp(1,1);
  matrixE PFrot_dp(1,1);
  matrixE DDrot_dp(1,1);
  matrixE DFrot_dp(1,1);
  matrixE FFrot_dp(1,1);
  //other auxiliaries
  std::vector<std::vector<int> > polynomial;
  double a;
  double b;
  double c;
  double argA;
  double argB;
  std::vector<double> Ak;
  std::vector<double> Bk;
  double aux;
  //determine maximum azimuthal quantum number
  size_t maxL1 = 0;
  size_t maxL2 = 0;
  for (size_t idx = 0; idx < QN1.rows(); ++idx) {
    if (QN1(idx + 1,2) > maxL1) {maxL1 = QN1(idx + 1,2);}
  }
  for (size_t idx = 0; idx < QN2.rows(); ++idx) {
    if (QN2(idx + 1,2) > maxL2) {maxL2 = QN2(idx + 1,2);}
  }
  if ((maxL1 > 0)||(maxL2 > 0)) {
    //in this case we have p orbitals
    SProt_dp = SPtransf_dp(cost,sint,cosp,sinp,rotate);
    PProt_dp = PPtransf_dp(cost,sint,cosp,sinp,rotate);
  }
  if ((maxL1 > 1)||(maxL2 > 1)) {
    //in this case we have d orbitals
    SDrot_dp = SDtransf_dp(cost,sint,cosp,sinp,rotate);
    PDrot_dp = PDtransf_dp(cost,sint,cosp,sinp,rotate);
    DDrot_dp = DDtransf_dp(cost,sint,cosp,sinp,rotate);
  }
  if ((maxL1 > 2)||(maxL2 > 2)) {
    //in this case we have f orbitals
    SFrot_dp = SFtransf_dp(cost,sint,cosp,sinp,rotate);
    PFrot_dp = PFtransf_dp(cost,sint,cosp,sinp,rotate);
    DFrot_dp = DFtransf_dp(cost,sint,cosp,sinp,rotate);
    FFrot_dp = FFtransf_dp(cost,sint,cosp,sinp,rotate);
  }
  size_t iNA;
  size_t iNB;
  size_t iLA;
  size_t iLB;
  size_t irw = 1;
  size_t icl = 1;
  for (size_t isbra = 0; isbra < QN1.rows(); ++isbra) {           //loop over bra shells
    a = exponents1[isbra];
    if (a == 0.0) continue;
    iNA = QN1(isbra + 1,1) - 1;
    iLA = QN1(isbra + 1,2);
    icl = 1;
    for (size_t isket = 0; isket < QN2.rows(); ++isket) {         //loop over ket shells
      b = exponents2[isket];
      if (b == 0.0) continue;
      iNB = QN2(isket + 1,1) - 1;
      iLB = QN2(isket + 1,2);
      //getting integrals
      argA = 0.5*dist*(a + b);
      argB = 0.5*dist*(a - b);
      Ak = int_1_infty_xpk_expmyx_dx(argA,iNA + iNB + 2);
      Bk = int_m1_1_xpk_expmyx_dx_Pople(argB,iNA + iNB + 2);
      aux = OvFactor(iNA + 1,iLA,0,a,iNB + 1,iLB,0,b,dist);
      //different azimuthal quantum number cases
      //ss = 0
      if ((iLA == 1)&&(iLB == 1)) {
        //pp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = PPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 9 dimensions
        overlaps(irw + 0,icl + 0) = -sigma_sigma*PProt_dp(1,1) + pi_pi*(PProt_dp(1,2) + PProt_dp(1,3));      //<x|x>
        overlaps(irw + 0,icl + 1) = -sigma_sigma*PProt_dp(2,1) + pi_pi*(PProt_dp(2,2) + PProt_dp(2,3));      //<x|y>
        overlaps(irw + 0,icl + 2) = -sigma_sigma*PProt_dp(3,1) + pi_pi*(PProt_dp(3,2) + PProt_dp(3,3));      //<x|z>
        overlaps(irw + 1,icl + 0) = -sigma_sigma*PProt_dp(2,1) + pi_pi*(PProt_dp(2,2) + PProt_dp(2,3));      //<y|x>
        overlaps(irw + 1,icl + 1) = -sigma_sigma*PProt_dp(4,1) + pi_pi*(PProt_dp(4,2) + PProt_dp(4,3));      //<y|y>
        overlaps(irw + 1,icl + 2) = -sigma_sigma*PProt_dp(5,1) + pi_pi*(PProt_dp(5,2) + PProt_dp(5,3));      //<y|z>
        overlaps(irw + 2,icl + 0) = -sigma_sigma*PProt_dp(3,1) + pi_pi*(PProt_dp(3,2) + PProt_dp(3,3));      //<z|x>
        overlaps(irw + 2,icl + 1) = -sigma_sigma*PProt_dp(5,1) + pi_pi*(PProt_dp(5,2) + PProt_dp(5,3));      //<z|y>
        overlaps(irw + 2,icl + 2) = -sigma_sigma*PProt_dp(6,1) + pi_pi*(PProt_dp(6,2) + PProt_dp(6,3));      //<z|z>
      }
      else if ((iLA == 2)&&(iLB == 2)) {
        //dd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = DDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = DDdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5x5 dimensions (xz,yz,z2,xy,x2-y2)(xz,yz,z2,xy,x2-y2)
        overlaps(irw + 0,icl + 0) = sigma_sigma*DDrot_dp(1,1)  - pi_pi*(DDrot_dp(1,2)  + DDrot_dp(1,3))  + delta_delta*(DDrot_dp(1,4)  + DDrot_dp(1,5));             //<xz|xz>
        overlaps(irw + 0,icl + 1) = sigma_sigma*DDrot_dp(2,1)  - pi_pi*(DDrot_dp(2,2)  + DDrot_dp(2,3))  + delta_delta*(DDrot_dp(2,4)  + DDrot_dp(2,5));             //<xz|yz>
        overlaps(irw + 0,icl + 2) = sigma_sigma*DDrot_dp(3,1)  - pi_pi*(DDrot_dp(3,2)  + DDrot_dp(3,3))  + delta_delta*(DDrot_dp(3,4)  + DDrot_dp(3,5));             //<xz|z2>
        overlaps(irw + 0,icl + 3) = sigma_sigma*DDrot_dp(4,1)  - pi_pi*(DDrot_dp(4,2)  + DDrot_dp(4,3))  + delta_delta*(DDrot_dp(4,4)  + DDrot_dp(4,5));             //<xz|xy>
        overlaps(irw + 0,icl + 4) = sigma_sigma*DDrot_dp(5,1)  - pi_pi*(DDrot_dp(5,2)  + DDrot_dp(5,3))  + delta_delta*(DDrot_dp(5,4)  + DDrot_dp(5,5));             //<xz|x2-y2>
        overlaps(irw + 1,icl + 0) = sigma_sigma*DDrot_dp(2,1)  - pi_pi*(DDrot_dp(2,2)  + DDrot_dp(2,3))  + delta_delta*(DDrot_dp(2,4)  + DDrot_dp(2,5));             //<yz|xz>
        overlaps(irw + 1,icl + 1) = sigma_sigma*DDrot_dp(6,1)  - pi_pi*(DDrot_dp(6,2)  + DDrot_dp(6,3))  + delta_delta*(DDrot_dp(6,4)  + DDrot_dp(6,5));             //<yz|yz>
        overlaps(irw + 1,icl + 2) = sigma_sigma*DDrot_dp(7,1)  - pi_pi*(DDrot_dp(7,2)  + DDrot_dp(7,3))  + delta_delta*(DDrot_dp(7,4)  + DDrot_dp(7,5));             //<yz|z2>
        overlaps(irw + 1,icl + 3) = sigma_sigma*DDrot_dp(8,1)  - pi_pi*(DDrot_dp(8,2)  + DDrot_dp(8,3))  + delta_delta*(DDrot_dp(8,4)  + DDrot_dp(8,5));             //<yz|xy>
        overlaps(irw + 1,icl + 4) = sigma_sigma*DDrot_dp(9,1)  - pi_pi*(DDrot_dp(9,2)  + DDrot_dp(9,3))  + delta_delta*(DDrot_dp(9,4)  + DDrot_dp(9,5));             //<yz|x2-y2>
        overlaps(irw + 2,icl + 0) = sigma_sigma*DDrot_dp(3,1)  - pi_pi*(DDrot_dp(3,2)  + DDrot_dp(3,3))  + delta_delta*(DDrot_dp(3,4)  + DDrot_dp(3,5));             //<z2|xz>
        overlaps(irw + 2,icl + 1) = sigma_sigma*DDrot_dp(7,1)  - pi_pi*(DDrot_dp(7,2)  + DDrot_dp(7,3))  + delta_delta*(DDrot_dp(7,4)  + DDrot_dp(7,5));             //<z2|yz>
        overlaps(irw + 2,icl + 2) = sigma_sigma*DDrot_dp(10,1) - pi_pi*(DDrot_dp(10,2) + DDrot_dp(10,3)) + delta_delta*(DDrot_dp(10,4) + DDrot_dp(10,5));            //<z2|z2>
        overlaps(irw + 2,icl + 3) = sigma_sigma*DDrot_dp(11,1) - pi_pi*(DDrot_dp(11,2) + DDrot_dp(11,3)) + delta_delta*(DDrot_dp(11,4) + DDrot_dp(11,5));            //<z2|xy>
        overlaps(irw + 2,icl + 4) = sigma_sigma*DDrot_dp(12,1) - pi_pi*(DDrot_dp(12,2) + DDrot_dp(12,3)) + delta_delta*(DDrot_dp(12,4) + DDrot_dp(12,5));            //<z2|x2-y2>
        overlaps(irw + 3,icl + 0) = sigma_sigma*DDrot_dp(4,1)  - pi_pi*(DDrot_dp(4,2)  + DDrot_dp(4,3))  + delta_delta*(DDrot_dp(4,4)  + DDrot_dp(4,5));             //<xy|xz>
        overlaps(irw + 3,icl + 1) = sigma_sigma*DDrot_dp(8,1)  - pi_pi*(DDrot_dp(8,2)  + DDrot_dp(8,3))  + delta_delta*(DDrot_dp(8,4)  + DDrot_dp(8,5));             //<xy|yz>
        overlaps(irw + 3,icl + 2) = sigma_sigma*DDrot_dp(11,1) - pi_pi*(DDrot_dp(11,2) + DDrot_dp(11,3)) + delta_delta*(DDrot_dp(11,4) + DDrot_dp(11,5));            //<xy|z2>
        overlaps(irw + 3,icl + 3) = sigma_sigma*DDrot_dp(13,1) - pi_pi*(DDrot_dp(13,2) + DDrot_dp(13,3)) + delta_delta*(DDrot_dp(13,4) + DDrot_dp(13,5));            //<xy|xy>
        overlaps(irw + 3,icl + 4) = sigma_sigma*DDrot_dp(14,1) - pi_pi*(DDrot_dp(14,2) + DDrot_dp(14,3)) + delta_delta*(DDrot_dp(14,4) + DDrot_dp(14,5));            //<xy|x2-y2>
        overlaps(irw + 4,icl + 0) = sigma_sigma*DDrot_dp(5,1)  - pi_pi*(DDrot_dp(5,2)  + DDrot_dp(5,3))  + delta_delta*(DDrot_dp(5,4)  + DDrot_dp(5,5));             //<x2-y2|xz>
        overlaps(irw + 4,icl + 1) = sigma_sigma*DDrot_dp(9,1)  - pi_pi*(DDrot_dp(9,2)  + DDrot_dp(9,3))  + delta_delta*(DDrot_dp(9,4)  + DDrot_dp(9,5));             //<x2-y2|yz>
        overlaps(irw + 4,icl + 2) = sigma_sigma*DDrot_dp(12,1) - pi_pi*(DDrot_dp(12,2) + DDrot_dp(12,3)) + delta_delta*(DDrot_dp(12,4) + DDrot_dp(12,5));            //<x2-y2|z2>
        overlaps(irw + 4,icl + 3) = sigma_sigma*DDrot_dp(14,1) - pi_pi*(DDrot_dp(14,2) + DDrot_dp(14,3)) + delta_delta*(DDrot_dp(14,4) + DDrot_dp(14,5));            //<x2-y2|xy>
        overlaps(irw + 4,icl + 4) = sigma_sigma*DDrot_dp(15,1) - pi_pi*(DDrot_dp(15,2) + DDrot_dp(15,3)) + delta_delta*(DDrot_dp(15,4) + DDrot_dp(15,5));            //<x2-y2|x2-y2>
      }
      else if ((iLA == 3)&&(iLB == 3)) {
        //ff
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = FFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = FFdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //phi integral
        aux = OvFactor(iNA + 1,iLA,3,a,iNB + 1,iLB,3,b,dist);
        phi_phi = 0.0;
        polynomial = FFphiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          phi_phi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 7x7 dimensions (fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        overlaps(irw + 0,icl + 0) = -sigma_sigma*FFrot_dp(1,1)  + pi_pi*(FFrot_dp(1,2)  + FFrot_dp(1,3))  - delta_delta*(FFrot_dp(1,4)  + FFrot_dp(1,5))  + phi_phi*(FFrot_dp(1,6)  + FFrot_dp(1,7));
        overlaps(irw + 0,icl + 1) = -sigma_sigma*FFrot_dp(2,1)  + pi_pi*(FFrot_dp(2,2)  + FFrot_dp(2,3))  - delta_delta*(FFrot_dp(2,4)  + FFrot_dp(2,5))  + phi_phi*(FFrot_dp(2,6)  + FFrot_dp(2,7));
        overlaps(irw + 0,icl + 2) = -sigma_sigma*FFrot_dp(3,1)  + pi_pi*(FFrot_dp(3,2)  + FFrot_dp(3,3))  - delta_delta*(FFrot_dp(3,4)  + FFrot_dp(3,5))  + phi_phi*(FFrot_dp(3,6)  + FFrot_dp(3,7));
        overlaps(irw + 0,icl + 3) = -sigma_sigma*FFrot_dp(4,1)  + pi_pi*(FFrot_dp(4,2)  + FFrot_dp(4,3))  - delta_delta*(FFrot_dp(4,4)  + FFrot_dp(4,5))  + phi_phi*(FFrot_dp(4,6)  + FFrot_dp(4,7));
        overlaps(irw + 0,icl + 4) = -sigma_sigma*FFrot_dp(5,1)  + pi_pi*(FFrot_dp(5,2)  + FFrot_dp(5,3))  - delta_delta*(FFrot_dp(5,4)  + FFrot_dp(5,5))  + phi_phi*(FFrot_dp(5,6)  + FFrot_dp(5,7));
        overlaps(irw + 0,icl + 5) = -sigma_sigma*FFrot_dp(6,1)  + pi_pi*(FFrot_dp(6,2)  + FFrot_dp(6,3))  - delta_delta*(FFrot_dp(6,4)  + FFrot_dp(6,5))  + phi_phi*(FFrot_dp(6,6)  + FFrot_dp(6,7));
        overlaps(irw + 0,icl + 6) = -sigma_sigma*FFrot_dp(7,1)  + pi_pi*(FFrot_dp(7,2)  + FFrot_dp(7,3))  - delta_delta*(FFrot_dp(7,4)  + FFrot_dp(7,5))  + phi_phi*(FFrot_dp(7,6)  + FFrot_dp(7,7));
        overlaps(irw + 1,icl + 0) = -sigma_sigma*FFrot_dp(2,1)  + pi_pi*(FFrot_dp(2,2)  + FFrot_dp(2,3))  - delta_delta*(FFrot_dp(2,4)  + FFrot_dp(2,5))  + phi_phi*(FFrot_dp(2,6)  + FFrot_dp(2,7));
        overlaps(irw + 1,icl + 1) = -sigma_sigma*FFrot_dp(8,1)  + pi_pi*(FFrot_dp(8,2)  + FFrot_dp(8,3))  - delta_delta*(FFrot_dp(8,4)  + FFrot_dp(8,5))  + phi_phi*(FFrot_dp(8,6)  + FFrot_dp(8,7));
        overlaps(irw + 1,icl + 2) = -sigma_sigma*FFrot_dp(9,1)  + pi_pi*(FFrot_dp(9,2)  + FFrot_dp(9,3))  - delta_delta*(FFrot_dp(9,4)  + FFrot_dp(9,5))  + phi_phi*(FFrot_dp(9,6)  + FFrot_dp(9,7));
        overlaps(irw + 1,icl + 3) = -sigma_sigma*FFrot_dp(10,1) + pi_pi*(FFrot_dp(10,2) + FFrot_dp(10,3)) - delta_delta*(FFrot_dp(10,4) + FFrot_dp(10,5)) + phi_phi*(FFrot_dp(10,6) + FFrot_dp(10,7));
        overlaps(irw + 1,icl + 4) = -sigma_sigma*FFrot_dp(11,1) + pi_pi*(FFrot_dp(11,2) + FFrot_dp(11,3)) - delta_delta*(FFrot_dp(11,4) + FFrot_dp(11,5)) + phi_phi*(FFrot_dp(11,6) + FFrot_dp(11,7));
        overlaps(irw + 1,icl + 5) = -sigma_sigma*FFrot_dp(12,1) + pi_pi*(FFrot_dp(12,2) + FFrot_dp(12,3)) - delta_delta*(FFrot_dp(12,4) + FFrot_dp(12,5)) + phi_phi*(FFrot_dp(12,6) + FFrot_dp(12,7));
        overlaps(irw + 1,icl + 6) = -sigma_sigma*FFrot_dp(13,1) + pi_pi*(FFrot_dp(13,2) + FFrot_dp(13,3)) - delta_delta*(FFrot_dp(13,4) + FFrot_dp(13,5)) + phi_phi*(FFrot_dp(13,6) + FFrot_dp(13,7));
        overlaps(irw + 2,icl + 0) = -sigma_sigma*FFrot_dp(3,1) +  pi_pi*(FFrot_dp(3,2)  + FFrot_dp(3,3))  - delta_delta*(FFrot_dp(3,4)  + FFrot_dp(3,5))  + phi_phi*(FFrot_dp(3,6)  + FFrot_dp(3,7));
        overlaps(irw + 2,icl + 1) = -sigma_sigma*FFrot_dp(9,1) +  pi_pi*(FFrot_dp(9,2)  + FFrot_dp(9,3))  - delta_delta*(FFrot_dp(9,4)  + FFrot_dp(9,5))  + phi_phi*(FFrot_dp(9,6)  + FFrot_dp(9,7));
        overlaps(irw + 2,icl + 2) = -sigma_sigma*FFrot_dp(14,1) + pi_pi*(FFrot_dp(14,2) + FFrot_dp(14,3)) - delta_delta*(FFrot_dp(14,4) + FFrot_dp(14,5)) + phi_phi*(FFrot_dp(14,6) + FFrot_dp(14,7));
        overlaps(irw + 2,icl + 3) = -sigma_sigma*FFrot_dp(15,1) + pi_pi*(FFrot_dp(15,2) + FFrot_dp(15,3)) - delta_delta*(FFrot_dp(15,4) + FFrot_dp(15,5)) + phi_phi*(FFrot_dp(15,6) + FFrot_dp(15,7));
        overlaps(irw + 2,icl + 4) = -sigma_sigma*FFrot_dp(16,1) + pi_pi*(FFrot_dp(16,2) + FFrot_dp(16,3)) - delta_delta*(FFrot_dp(16,4) + FFrot_dp(16,5)) + phi_phi*(FFrot_dp(16,6) + FFrot_dp(16,7));
        overlaps(irw + 2,icl + 5) = -sigma_sigma*FFrot_dp(17,1) + pi_pi*(FFrot_dp(17,2) + FFrot_dp(17,3)) - delta_delta*(FFrot_dp(17,4) + FFrot_dp(17,5)) + phi_phi*(FFrot_dp(17,6) + FFrot_dp(17,7));
        overlaps(irw + 2,icl + 6) = -sigma_sigma*FFrot_dp(18,1) + pi_pi*(FFrot_dp(18,2) + FFrot_dp(18,3)) - delta_delta*(FFrot_dp(18,4) + FFrot_dp(18,5)) + phi_phi*(FFrot_dp(18,6) + FFrot_dp(18,7));
        overlaps(irw + 3,icl + 0) = -sigma_sigma*FFrot_dp(4,1) +  pi_pi*(FFrot_dp(4,2)  + FFrot_dp(4,3))  - delta_delta*(FFrot_dp(4,4)  + FFrot_dp(4,5))  + phi_phi*(FFrot_dp(4,6)  + FFrot_dp(4,7));
        overlaps(irw + 3,icl + 1) = -sigma_sigma*FFrot_dp(10,1) + pi_pi*(FFrot_dp(10,2) + FFrot_dp(10,3)) - delta_delta*(FFrot_dp(10,4) + FFrot_dp(10,5)) + phi_phi*(FFrot_dp(10,6) + FFrot_dp(10,7));
        overlaps(irw + 3,icl + 2) = -sigma_sigma*FFrot_dp(15,1) + pi_pi*(FFrot_dp(15,2) + FFrot_dp(15,3)) - delta_delta*(FFrot_dp(15,4) + FFrot_dp(15,5)) + phi_phi*(FFrot_dp(15,6) + FFrot_dp(15,7));
        overlaps(irw + 3,icl + 3) = -sigma_sigma*FFrot_dp(19,1) + pi_pi*(FFrot_dp(19,2) + FFrot_dp(19,3)) - delta_delta*(FFrot_dp(19,4) + FFrot_dp(19,5)) + phi_phi*(FFrot_dp(19,6) + FFrot_dp(19,7));
        overlaps(irw + 3,icl + 4) = -sigma_sigma*FFrot_dp(20,1) + pi_pi*(FFrot_dp(20,2) + FFrot_dp(20,3)) - delta_delta*(FFrot_dp(20,4) + FFrot_dp(20,5)) + phi_phi*(FFrot_dp(20,6) + FFrot_dp(20,7));
        overlaps(irw + 3,icl + 5) = -sigma_sigma*FFrot_dp(21,1) + pi_pi*(FFrot_dp(21,2) + FFrot_dp(21,3)) - delta_delta*(FFrot_dp(21,4) + FFrot_dp(21,5)) + phi_phi*(FFrot_dp(21,6) + FFrot_dp(21,7));
        overlaps(irw + 3,icl + 6) = -sigma_sigma*FFrot_dp(22,1) + pi_pi*(FFrot_dp(22,2) + FFrot_dp(22,3)) - delta_delta*(FFrot_dp(22,4) + FFrot_dp(22,5)) + phi_phi*(FFrot_dp(22,6) + FFrot_dp(22,7));
        overlaps(irw + 4,icl + 0) = -sigma_sigma*FFrot_dp(5,1) +  pi_pi*(FFrot_dp(5,2)  + FFrot_dp(5,3))  - delta_delta*(FFrot_dp(5,4)  + FFrot_dp(5,5))  + phi_phi*(FFrot_dp(5,6)  + FFrot_dp(5,7));
        overlaps(irw + 4,icl + 1) = -sigma_sigma*FFrot_dp(11,1) + pi_pi*(FFrot_dp(11,2) + FFrot_dp(11,3)) - delta_delta*(FFrot_dp(11,4) + FFrot_dp(11,5)) + phi_phi*(FFrot_dp(11,6) + FFrot_dp(11,7));
        overlaps(irw + 4,icl + 2) = -sigma_sigma*FFrot_dp(16,1) + pi_pi*(FFrot_dp(16,2) + FFrot_dp(16,3)) - delta_delta*(FFrot_dp(16,4) + FFrot_dp(16,5)) + phi_phi*(FFrot_dp(16,6) + FFrot_dp(16,7));
        overlaps(irw + 4,icl + 3) = -sigma_sigma*FFrot_dp(20,1) + pi_pi*(FFrot_dp(20,2) + FFrot_dp(20,3)) - delta_delta*(FFrot_dp(20,4) + FFrot_dp(20,5)) + phi_phi*(FFrot_dp(20,6) + FFrot_dp(20,7));
        overlaps(irw + 4,icl + 4) = -sigma_sigma*FFrot_dp(23,1) + pi_pi*(FFrot_dp(23,2) + FFrot_dp(23,3)) - delta_delta*(FFrot_dp(23,4) + FFrot_dp(23,5)) + phi_phi*(FFrot_dp(23,6) + FFrot_dp(23,7));
        overlaps(irw + 4,icl + 5) = -sigma_sigma*FFrot_dp(24,1) + pi_pi*(FFrot_dp(24,2) + FFrot_dp(24,3)) - delta_delta*(FFrot_dp(24,4) + FFrot_dp(24,5)) + phi_phi*(FFrot_dp(24,6) + FFrot_dp(24,7));
        overlaps(irw + 4,icl + 6) = -sigma_sigma*FFrot_dp(25,1) + pi_pi*(FFrot_dp(25,2) + FFrot_dp(25,3)) - delta_delta*(FFrot_dp(25,4) + FFrot_dp(25,5)) + phi_phi*(FFrot_dp(25,6) + FFrot_dp(25,7));
        overlaps(irw + 5,icl + 0) = -sigma_sigma*FFrot_dp(6,1) +  pi_pi*(FFrot_dp(6,2)  + FFrot_dp(6,3))  - delta_delta*(FFrot_dp(6,4)  + FFrot_dp(6,5))  + phi_phi*(FFrot_dp(6,6)  + FFrot_dp(6,7));
        overlaps(irw + 5,icl + 1) = -sigma_sigma*FFrot_dp(12,1) + pi_pi*(FFrot_dp(12,2) + FFrot_dp(12,3)) - delta_delta*(FFrot_dp(12,4) + FFrot_dp(12,5)) + phi_phi*(FFrot_dp(12,6) + FFrot_dp(12,7));
        overlaps(irw + 5,icl + 2) = -sigma_sigma*FFrot_dp(17,1) + pi_pi*(FFrot_dp(17,2) + FFrot_dp(17,3)) - delta_delta*(FFrot_dp(17,4) + FFrot_dp(17,5)) + phi_phi*(FFrot_dp(17,6) + FFrot_dp(17,7));
        overlaps(irw + 5,icl + 3) = -sigma_sigma*FFrot_dp(21,1) + pi_pi*(FFrot_dp(21,2) + FFrot_dp(21,3)) - delta_delta*(FFrot_dp(21,4) + FFrot_dp(21,5)) + phi_phi*(FFrot_dp(21,6) + FFrot_dp(21,7));
        overlaps(irw + 5,icl + 4) = -sigma_sigma*FFrot_dp(24,1) + pi_pi*(FFrot_dp(24,2) + FFrot_dp(24,3)) - delta_delta*(FFrot_dp(24,4) + FFrot_dp(24,5)) + phi_phi*(FFrot_dp(24,6) + FFrot_dp(24,7));
        overlaps(irw + 5,icl + 5) = -sigma_sigma*FFrot_dp(26,1) + pi_pi*(FFrot_dp(26,2) + FFrot_dp(26,3)) - delta_delta*(FFrot_dp(26,4) + FFrot_dp(26,5)) + phi_phi*(FFrot_dp(26,6) + FFrot_dp(26,7));
        overlaps(irw + 5,icl + 6) = -sigma_sigma*FFrot_dp(27,1) + pi_pi*(FFrot_dp(27,2) + FFrot_dp(27,3)) - delta_delta*(FFrot_dp(27,4) + FFrot_dp(27,5)) + phi_phi*(FFrot_dp(27,6) + FFrot_dp(27,7));
        overlaps(irw + 6,icl + 0) = -sigma_sigma*FFrot_dp(7,1) +  pi_pi*(FFrot_dp(7,2)  + FFrot_dp(7,3))  - delta_delta*(FFrot_dp(7,4)  + FFrot_dp(7,5))  + phi_phi*(FFrot_dp(7,6)  + FFrot_dp(7,7));
        overlaps(irw + 6,icl + 1) = -sigma_sigma*FFrot_dp(13,1) + pi_pi*(FFrot_dp(13,2) + FFrot_dp(13,3)) - delta_delta*(FFrot_dp(13,4) + FFrot_dp(13,5)) + phi_phi*(FFrot_dp(13,6) + FFrot_dp(13,7));
        overlaps(irw + 6,icl + 2) = -sigma_sigma*FFrot_dp(18,1) + pi_pi*(FFrot_dp(18,2) + FFrot_dp(18,3)) - delta_delta*(FFrot_dp(18,4) + FFrot_dp(18,5)) + phi_phi*(FFrot_dp(18,6) + FFrot_dp(18,7));
        overlaps(irw + 6,icl + 3) = -sigma_sigma*FFrot_dp(22,1) + pi_pi*(FFrot_dp(22,2) + FFrot_dp(22,3)) - delta_delta*(FFrot_dp(22,4) + FFrot_dp(22,5)) + phi_phi*(FFrot_dp(22,6) + FFrot_dp(22,7));
        overlaps(irw + 6,icl + 4) = -sigma_sigma*FFrot_dp(25,1) + pi_pi*(FFrot_dp(25,2) + FFrot_dp(25,3)) - delta_delta*(FFrot_dp(25,4) + FFrot_dp(25,5)) + phi_phi*(FFrot_dp(25,6) + FFrot_dp(25,7));
        overlaps(irw + 6,icl + 5) = -sigma_sigma*FFrot_dp(27,1) + pi_pi*(FFrot_dp(27,2) + FFrot_dp(27,3)) - delta_delta*(FFrot_dp(27,4) + FFrot_dp(27,5)) + phi_phi*(FFrot_dp(27,6) + FFrot_dp(27,7));
        overlaps(irw + 6,icl + 6) = -sigma_sigma*FFrot_dp(28,1) + pi_pi*(FFrot_dp(28,2) + FFrot_dp(28,3)) - delta_delta*(FFrot_dp(28,4) + FFrot_dp(28,5)) + phi_phi*(FFrot_dp(28,6) + FFrot_dp(28,7));
      }
      else if ((iLA == 0)&&(iLB == 1)) {
        //sp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SPOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3 dimensions (px,py,pz)
        overlaps(irw,icl + 0) = -SProt_dp(1,1)*sigma_sigma;
        overlaps(irw,icl + 1) = -SProt_dp(2,1)*sigma_sigma;
        overlaps(irw,icl + 2) = -SProt_dp(3,1)*sigma_sigma;
      }
      else if ((iLA == 1)&&(iLB == 0)) {
        //ps
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3 dimensions (px,py,pz)
        overlaps(irw + 0,icl) = SProt_dp(1,1)*sigma_sigma;
        overlaps(irw + 1,icl) = SProt_dp(2,1)*sigma_sigma;
        overlaps(irw + 2,icl) = SProt_dp(3,1)*sigma_sigma;
      }
      else if ((iLA == 0)&&(iLB == 2)) {
        //sd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SDOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          overlaps(irw,icl + index) = SDrot_dp(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 2)&&(iLB == 0)) {
        //ds
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          overlaps(irw + index,icl) = SDrot_dp(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 0)&&(iLB == 3)) {
        //sf
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SFOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 7 dimensions (fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 7; ++index) {
          overlaps(irw,icl + index) = -SFrot_dp(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 3)&&(iLB == 0)) {
        //fs
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 7; ++index) {
          overlaps(irw + index,icl) = SFrot_dp(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 1)&&(iLB == 2)) {
        //pd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = PDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 3; ++index2) {
            overlaps(irw + index2,icl + index) = PDrot_dp(1+index2+3*index,1)*sigma_sigma - (PDrot_dp(1+index2+3*index,2) + PDrot_dp(1+index2+3*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 2)&&(iLB == 1)) {
        //dp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = DPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 3; ++index2) {
            overlaps(irw + index,icl + index2) = -PDrot_dp(1+index2+3*index,1)*sigma_sigma + (PDrot_dp(1+index2+3*index,2) + PDrot_dp(1+index2+3*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 1)&&(iLB == 3)) {
        //pf
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = PFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 3; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index,icl + index2) = -PFrot_dp(1+index2+7*index,1)*sigma_sigma + (PFrot_dp(1+index2+7*index,2) + PFrot_dp(1+index2+7*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 3)&&(iLB == 1)) {
        //fp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = FPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 3; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index2,icl + index) = -PFrot_dp(1+index2+7*index,1)*sigma_sigma + (PFrot_dp(1+index2+7*index,2) + PFrot_dp(1+index2+7*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 2)&&(iLB == 3)) {
        //df
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = DFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = DFdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (xz,yz,z2,xy,x2-y2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index,icl + index2) = -DFrot_dp(1+index2+7*index,1)*sigma_sigma + (DFrot_dp(1+index2+7*index,2) + DFrot_dp(1+index2+7*index,3))*pi_pi - (DFrot_dp(1+index2+7*index,4) + DFrot_dp(1+index2+7*index,5))*delta_delta;
          }
        }
      }
      else if ((iLA == 3)&&(iLB == 2)) {
        //fd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = FDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = FDdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (xz,yz,z2,xy,x2-y2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index2,icl + index) = DFrot_dp(1+index2+7*index,1)*sigma_sigma - (DFrot_dp(1+index2+7*index,2) + DFrot_dp(1+index2+7*index,3))*pi_pi + (DFrot_dp(1+index2+7*index,4) + DFrot_dp(1+index2+7*index,5))*delta_delta;
          }
        }
      }
      icl += 2*QN2(isket + 1,2) + 1;
    }
    irw += 2*QN1(isbra + 1,2) + 1;
  }
}
void TwoCenterOv_dR2(matrixE & overlapdR2, size_t nao1, size_t nao2, double dist, std::vector<double> & exponents1, matrix<size_t> & QN1, std::vector<double> & exponents2, matrix<size_t> & QN2, std::vector<double> & orientation, bool rotate = true, double tolerance = 1e-6) {
  //implementation of second-derivatives of overlap matrix with respect to internuclear distance; only up to shell f
  //calculation of the derivative of the overlap function between atoms A and B (A -> 1; B -> 2) with respect to internuclear distance
  overlapdR2.resize(nao1,nao2);
  overlapdR2.zero();
  //auxiliary variables for integrals
  double sigma_sigma = 0.0;
  double pi_pi = 0.0;
  double delta_delta = 0.0;
  double phi_phi = 0.0;
  //angles needed to reorient the integrals according to the position of atoms in the molecule
  double cost = orientation[2];
  double sint = sqrt(1 - cost*cost);
  double cosp = 1.0;
  double sinp = 0.0;
  if (fabs(sint) > tolerance) {
    cosp = orientation[0]/sint;
    sinp = orientation[1]/sint;
  }
  //matrices with overlap rotations
  matrixE SProt(1,1);
  matrixE SDrot(1,1);
  matrixE SFrot(1,1);
  matrixE PProt(1,1);
  matrixE PDrot(1,1);
  matrixE PFrot(1,1);
  matrixE DDrot(1,1);
  matrixE DFrot(1,1);
  matrixE FFrot(1,1);
  //other auxiliaries
  std::vector<std::vector<int> > polynomial;
  double a;
  double b;
  double c;
  double argA;
  double argB;
  std::vector<double> Ak;
  std::vector<double> Bk;
  double aux;
  //determine maximum azimuthal quantum number
  size_t maxL1 = 0;
  size_t maxL2 = 0;
  for (size_t idx = 0; idx < QN1.rows(); ++idx) {
    if (QN1(idx + 1,2) > maxL1) {maxL1 = QN1(idx + 1,2);}
  }
  for (size_t idx = 0; idx < QN2.rows(); ++idx) {
    if (QN2(idx + 1,2) > maxL2) {maxL2 = QN2(idx + 1,2);}
  }
  if ((maxL1 > 0)||(maxL2 > 0)) {
    //in this case we have p orbitals
    SProt = SPtransf(cost,sint,cosp,sinp,rotate);               //the orientation vector already contains this transformation, so this would not be actually necessary; however, a 3D vector is basically for free and like this I have debug options
    PProt = PPtransf(cost,sint,cosp,sinp,rotate);
  }
  if ((maxL1 > 1)||(maxL2 > 1)) {
    //in this case we have d orbitals
    SDrot = SDtransf(cost,sint,cosp,sinp,rotate);
    PDrot = PDtransf(cost,sint,cosp,sinp,rotate);
    DDrot = DDtransf(cost,sint,cosp,sinp,rotate);
  }
  if ((maxL1 > 2)||(maxL2 > 2)) {
    //in this case we have f orbitals
    SFrot = SFtransf(cost,sint,cosp,sinp,rotate);
    PFrot = PFtransf(cost,sint,cosp,sinp,rotate);
    DFrot = DFtransf(cost,sint,cosp,sinp,rotate);
    FFrot = FFtransf(cost,sint,cosp,sinp,rotate);
  }
  size_t iNA;
  size_t iNB;
  size_t iLA;
  size_t iLB;
  size_t irw = 1;
  size_t icl = 1;
  int iNANB;
  for (size_t isbra = 0; isbra < QN1.rows(); ++isbra) {           //loop over bra shells
    a = exponents1[isbra];
    if (a == 0.0) continue;
    iNA = QN1(isbra + 1,1) - 1;             //note that iNA is subtracted -1 !
    iLA = QN1(isbra + 1,2);
    icl = 1;
    for (size_t isket = 0; isket < QN2.rows(); ++isket) {         //loop over ket shells
      b = exponents2[isket];
      if (b == 0.0) continue;
      iNB = QN2(isket + 1,1) - 1;           //note that iNB is subtracted -1 !
      iLB = QN2(isket + 1,2);
      iNANB = iNA + iNB;
      //getting integrals
      argA = 0.5*dist*(a + b);
      argB = 0.5*dist*(a - b);
      Ak = int_1_infty_xpk_expmyx_dx(argA,iNA + iNB + 4);
      Bk = int_m1_1_xpk_expmyx_dx_Pople(argB,iNA + iNB + 4);
      aux = OvFactor(iNA + 1,iLA,0,a,iNB + 1,iLB,0,b,dist)/dist;
      //different azimuthal quantum number cases
      if ((iLA == 0)&&(iLB == 0)) {
        //ss
        polynomial = SSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          overlapdR2(irw,icl) += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          overlapdR2(irw,icl) -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          overlapdR2(irw,icl) -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          overlapdR2(irw,icl) += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          overlapdR2(irw,icl) += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          overlapdR2(irw,icl) += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
      }
      else if ((iLA == 1)&&(iLB == 1)) {
        //pp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = PPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          pi_pi -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          pi_pi -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          pi_pi += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          pi_pi += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          pi_pi += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //partition over 9 dimensions
        overlapdR2(irw + 0,icl + 0) = -sigma_sigma*PProt(1,1) + pi_pi*(PProt(1,2) + PProt(1,3));      //<x|x>
        overlapdR2(irw + 0,icl + 1) = -sigma_sigma*PProt(2,1) + pi_pi*(PProt(2,2) + PProt(2,3));      //<x|y>
        overlapdR2(irw + 0,icl + 2) = -sigma_sigma*PProt(3,1) + pi_pi*(PProt(3,2) + PProt(3,3));      //<x|z>
        overlapdR2(irw + 1,icl + 0) = -sigma_sigma*PProt(2,1) + pi_pi*(PProt(2,2) + PProt(2,3));      //<y|x>
        overlapdR2(irw + 1,icl + 1) = -sigma_sigma*PProt(4,1) + pi_pi*(PProt(4,2) + PProt(4,3));      //<y|y>
        overlapdR2(irw + 1,icl + 2) = -sigma_sigma*PProt(5,1) + pi_pi*(PProt(5,2) + PProt(5,3));      //<y|z>
        overlapdR2(irw + 2,icl + 0) = -sigma_sigma*PProt(3,1) + pi_pi*(PProt(3,2) + PProt(3,3));      //<z|x>
        overlapdR2(irw + 2,icl + 1) = -sigma_sigma*PProt(5,1) + pi_pi*(PProt(5,2) + PProt(5,3));      //<z|y>
        overlapdR2(irw + 2,icl + 2) = -sigma_sigma*PProt(6,1) + pi_pi*(PProt(6,2) + PProt(6,3));      //<z|z>
      }
      else if ((iLA == 2)&&(iLB == 2)) {
        //dd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = DDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          pi_pi -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          pi_pi -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          pi_pi += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          pi_pi += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          pi_pi += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist)/dist;
        delta_delta = 0.0;
        polynomial = DDdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          delta_delta -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          delta_delta -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          delta_delta += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          delta_delta += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          delta_delta += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //partition over 5x5 dimensions (xz,yz,z2,xy,x2-y2)(xz,yz,z2,xy,x2-y2)
        overlapdR2(irw + 0,icl + 0) = sigma_sigma*DDrot(1,1)  - pi_pi*(DDrot(1,2)  + DDrot(1,3))  + delta_delta*(DDrot(1,4)  + DDrot(1,5));             //<xz|xz>
        overlapdR2(irw + 0,icl + 1) = sigma_sigma*DDrot(2,1)  - pi_pi*(DDrot(2,2)  + DDrot(2,3))  + delta_delta*(DDrot(2,4)  + DDrot(2,5));             //<xz|yz>
        overlapdR2(irw + 0,icl + 2) = sigma_sigma*DDrot(3,1)  - pi_pi*(DDrot(3,2)  + DDrot(3,3))  + delta_delta*(DDrot(3,4)  + DDrot(3,5));             //<xz|z2>
        overlapdR2(irw + 0,icl + 3) = sigma_sigma*DDrot(4,1)  - pi_pi*(DDrot(4,2)  + DDrot(4,3))  + delta_delta*(DDrot(4,4)  + DDrot(4,5));             //<xz|xy>
        overlapdR2(irw + 0,icl + 4) = sigma_sigma*DDrot(5,1)  - pi_pi*(DDrot(5,2)  + DDrot(5,3))  + delta_delta*(DDrot(5,4)  + DDrot(5,5));             //<xz|x2-y2>
        overlapdR2(irw + 1,icl + 0) = sigma_sigma*DDrot(2,1)  - pi_pi*(DDrot(2,2)  + DDrot(2,3))  + delta_delta*(DDrot(2,4)  + DDrot(2,5));             //<yz|xz>
        overlapdR2(irw + 1,icl + 1) = sigma_sigma*DDrot(6,1)  - pi_pi*(DDrot(6,2)  + DDrot(6,3))  + delta_delta*(DDrot(6,4)  + DDrot(6,5));             //<yz|yz>
        overlapdR2(irw + 1,icl + 2) = sigma_sigma*DDrot(7,1)  - pi_pi*(DDrot(7,2)  + DDrot(7,3))  + delta_delta*(DDrot(7,4)  + DDrot(7,5));             //<yz|z2>
        overlapdR2(irw + 1,icl + 3) = sigma_sigma*DDrot(8,1)  - pi_pi*(DDrot(8,2)  + DDrot(8,3))  + delta_delta*(DDrot(8,4)  + DDrot(8,5));             //<yz|xy>
        overlapdR2(irw + 1,icl + 4) = sigma_sigma*DDrot(9,1)  - pi_pi*(DDrot(9,2)  + DDrot(9,3))  + delta_delta*(DDrot(9,4)  + DDrot(9,5));             //<yz|x2-y2>
        overlapdR2(irw + 2,icl + 0) = sigma_sigma*DDrot(3,1)  - pi_pi*(DDrot(3,2)  + DDrot(3,3))  + delta_delta*(DDrot(3,4)  + DDrot(3,5));             //<z2|xz>
        overlapdR2(irw + 2,icl + 1) = sigma_sigma*DDrot(7,1)  - pi_pi*(DDrot(7,2)  + DDrot(7,3))  + delta_delta*(DDrot(7,4)  + DDrot(7,5));             //<z2|yz>
        overlapdR2(irw + 2,icl + 2) = sigma_sigma*DDrot(10,1) - pi_pi*(DDrot(10,2) + DDrot(10,3)) + delta_delta*(DDrot(10,4) + DDrot(10,5));            //<z2|z2>
        overlapdR2(irw + 2,icl + 3) = sigma_sigma*DDrot(11,1) - pi_pi*(DDrot(11,2) + DDrot(11,3)) + delta_delta*(DDrot(11,4) + DDrot(11,5));            //<z2|xy>
        overlapdR2(irw + 2,icl + 4) = sigma_sigma*DDrot(12,1) - pi_pi*(DDrot(12,2) + DDrot(12,3)) + delta_delta*(DDrot(12,4) + DDrot(12,5));            //<z2|x2-y2>
        overlapdR2(irw + 3,icl + 0) = sigma_sigma*DDrot(4,1)  - pi_pi*(DDrot(4,2)  + DDrot(4,3))  + delta_delta*(DDrot(4,4)  + DDrot(4,5));             //<xy|xz>
        overlapdR2(irw + 3,icl + 1) = sigma_sigma*DDrot(8,1)  - pi_pi*(DDrot(8,2)  + DDrot(8,3))  + delta_delta*(DDrot(8,4)  + DDrot(8,5));             //<xy|yz>
        overlapdR2(irw + 3,icl + 2) = sigma_sigma*DDrot(11,1) - pi_pi*(DDrot(11,2) + DDrot(11,3)) + delta_delta*(DDrot(11,4) + DDrot(11,5));            //<xy|z2>
        overlapdR2(irw + 3,icl + 3) = sigma_sigma*DDrot(13,1) - pi_pi*(DDrot(13,2) + DDrot(13,3)) + delta_delta*(DDrot(13,4) + DDrot(13,5));            //<xy|xy>
        overlapdR2(irw + 3,icl + 4) = sigma_sigma*DDrot(14,1) - pi_pi*(DDrot(14,2) + DDrot(14,3)) + delta_delta*(DDrot(14,4) + DDrot(14,5));            //<xy|x2-y2>
        overlapdR2(irw + 4,icl + 0) = sigma_sigma*DDrot(5,1)  - pi_pi*(DDrot(5,2)  + DDrot(5,3))  + delta_delta*(DDrot(5,4)  + DDrot(5,5));             //<x2-y2|xz>
        overlapdR2(irw + 4,icl + 1) = sigma_sigma*DDrot(9,1)  - pi_pi*(DDrot(9,2)  + DDrot(9,3))  + delta_delta*(DDrot(9,4)  + DDrot(9,5));             //<x2-y2|yz>
        overlapdR2(irw + 4,icl + 2) = sigma_sigma*DDrot(12,1) - pi_pi*(DDrot(12,2) + DDrot(12,3)) + delta_delta*(DDrot(12,4) + DDrot(12,5));            //<x2-y2|z2>
        overlapdR2(irw + 4,icl + 3) = sigma_sigma*DDrot(14,1) - pi_pi*(DDrot(14,2) + DDrot(14,3)) + delta_delta*(DDrot(14,4) + DDrot(14,5));            //<x2-y2|xy>
        overlapdR2(irw + 4,icl + 4) = sigma_sigma*DDrot(15,1) - pi_pi*(DDrot(15,2) + DDrot(15,3)) + delta_delta*(DDrot(15,4) + DDrot(15,5));            //<x2-y2|x2-y2>
      }
      else if ((iLA == 3)&&(iLB == 3)) {
        //ff
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = FFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          pi_pi -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          pi_pi -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          pi_pi += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          pi_pi += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          pi_pi += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist)/dist;
        delta_delta = 0.0;
        polynomial = FFdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          delta_delta -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          delta_delta -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          delta_delta += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          delta_delta += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          delta_delta += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //phi integral
        aux = OvFactor(iNA + 1,iLA,3,a,iNB + 1,iLB,3,b,dist)/dist;
        phi_phi = 0.0;
        polynomial = FFphiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          phi_phi += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          phi_phi -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          phi_phi -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          phi_phi += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          phi_phi += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          phi_phi += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //partition over 7x7 dimensions (fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        overlapdR2(irw + 0,icl + 0) = -sigma_sigma*FFrot(1,1) + pi_pi*(FFrot(1,2) + FFrot(1,3)) - delta_delta*(FFrot(1,4) + FFrot(1,5)) + phi_phi*(FFrot(1,6) + FFrot(1,7));
        overlapdR2(irw + 0,icl + 1) = -sigma_sigma*FFrot(2,1) + pi_pi*(FFrot(2,2) + FFrot(2,3)) - delta_delta*(FFrot(2,4) + FFrot(2,5)) + phi_phi*(FFrot(2,6) + FFrot(2,7));
        overlapdR2(irw + 0,icl + 2) = -sigma_sigma*FFrot(3,1) + pi_pi*(FFrot(3,2) + FFrot(3,3)) - delta_delta*(FFrot(3,4) + FFrot(3,5)) + phi_phi*(FFrot(3,6) + FFrot(3,7));
        overlapdR2(irw + 0,icl + 3) = -sigma_sigma*FFrot(4,1) + pi_pi*(FFrot(4,2) + FFrot(4,3)) - delta_delta*(FFrot(4,4) + FFrot(4,5)) + phi_phi*(FFrot(4,6) + FFrot(4,7));
        overlapdR2(irw + 0,icl + 4) = -sigma_sigma*FFrot(5,1) + pi_pi*(FFrot(5,2) + FFrot(5,3)) - delta_delta*(FFrot(5,4) + FFrot(5,5)) + phi_phi*(FFrot(5,6) + FFrot(5,7));
        overlapdR2(irw + 0,icl + 5) = -sigma_sigma*FFrot(6,1) + pi_pi*(FFrot(6,2) + FFrot(6,3)) - delta_delta*(FFrot(6,4) + FFrot(6,5)) + phi_phi*(FFrot(6,6) + FFrot(6,7));
        overlapdR2(irw + 0,icl + 6) = -sigma_sigma*FFrot(7,1) + pi_pi*(FFrot(7,2) + FFrot(7,3)) - delta_delta*(FFrot(7,4) + FFrot(7,5)) + phi_phi*(FFrot(7,6) + FFrot(7,7));
        overlapdR2(irw + 1,icl + 0) = -sigma_sigma*FFrot(2,1) + pi_pi*(FFrot(2,2) + FFrot(2,3)) - delta_delta*(FFrot(2,4) + FFrot(2,5)) + phi_phi*(FFrot(2,6) + FFrot(2,7));
        overlapdR2(irw + 1,icl + 1) = -sigma_sigma*FFrot(8,1) + pi_pi*(FFrot(8,2) + FFrot(8,3)) - delta_delta*(FFrot(8,4) + FFrot(8,5)) + phi_phi*(FFrot(8,6) + FFrot(8,7));
        overlapdR2(irw + 1,icl + 2) = -sigma_sigma*FFrot(9,1) + pi_pi*(FFrot(9,2) + FFrot(9,3)) - delta_delta*(FFrot(9,4) + FFrot(9,5)) + phi_phi*(FFrot(9,6) + FFrot(9,7));
        overlapdR2(irw + 1,icl + 3) = -sigma_sigma*FFrot(10,1) + pi_pi*(FFrot(10,2) + FFrot(10,3)) - delta_delta*(FFrot(10,4) + FFrot(10,5)) + phi_phi*(FFrot(10,6) + FFrot(10,7));
        overlapdR2(irw + 1,icl + 4) = -sigma_sigma*FFrot(11,1) + pi_pi*(FFrot(11,2) + FFrot(11,3)) - delta_delta*(FFrot(11,4) + FFrot(11,5)) + phi_phi*(FFrot(11,6) + FFrot(11,7));
        overlapdR2(irw + 1,icl + 5) = -sigma_sigma*FFrot(12,1) + pi_pi*(FFrot(12,2) + FFrot(12,3)) - delta_delta*(FFrot(12,4) + FFrot(12,5)) + phi_phi*(FFrot(12,6) + FFrot(12,7));
        overlapdR2(irw + 1,icl + 6) = -sigma_sigma*FFrot(13,1) + pi_pi*(FFrot(13,2) + FFrot(13,3)) - delta_delta*(FFrot(13,4) + FFrot(13,5)) + phi_phi*(FFrot(13,6) + FFrot(13,7));
        overlapdR2(irw + 2,icl + 0) = -sigma_sigma*FFrot(3,1) +  pi_pi*(FFrot(3,2) + FFrot(3,3)) - delta_delta*(FFrot(3,4) + FFrot(3,5)) + phi_phi*(FFrot(3,6) + FFrot(3,7));
        overlapdR2(irw + 2,icl + 1) = -sigma_sigma*FFrot(9,1) +  pi_pi*(FFrot(9,2) + FFrot(9,3)) - delta_delta*(FFrot(9,4) + FFrot(9,5)) + phi_phi*(FFrot(9,6) + FFrot(9,7));
        overlapdR2(irw + 2,icl + 2) = -sigma_sigma*FFrot(14,1) + pi_pi*(FFrot(14,2) + FFrot(14,3)) - delta_delta*(FFrot(14,4) + FFrot(14,5)) + phi_phi*(FFrot(14,6) + FFrot(14,7));
        overlapdR2(irw + 2,icl + 3) = -sigma_sigma*FFrot(15,1) + pi_pi*(FFrot(15,2) + FFrot(15,3)) - delta_delta*(FFrot(15,4) + FFrot(15,5)) + phi_phi*(FFrot(15,6) + FFrot(15,7));
        overlapdR2(irw + 2,icl + 4) = -sigma_sigma*FFrot(16,1) + pi_pi*(FFrot(16,2) + FFrot(16,3)) - delta_delta*(FFrot(16,4) + FFrot(16,5)) + phi_phi*(FFrot(16,6) + FFrot(16,7));
        overlapdR2(irw + 2,icl + 5) = -sigma_sigma*FFrot(17,1) + pi_pi*(FFrot(17,2) + FFrot(17,3)) - delta_delta*(FFrot(17,4) + FFrot(17,5)) + phi_phi*(FFrot(17,6) + FFrot(17,7));
        overlapdR2(irw + 2,icl + 6) = -sigma_sigma*FFrot(18,1) + pi_pi*(FFrot(18,2) + FFrot(18,3)) - delta_delta*(FFrot(18,4) + FFrot(18,5)) + phi_phi*(FFrot(18,6) + FFrot(18,7));
        overlapdR2(irw + 3,icl + 0) = -sigma_sigma*FFrot(4,1) +  pi_pi*(FFrot(4,2) + FFrot(4,3)) - delta_delta*(FFrot(4,4) + FFrot(4,5)) + phi_phi*(FFrot(4,6) + FFrot(4,7));
        overlapdR2(irw + 3,icl + 1) = -sigma_sigma*FFrot(10,1) + pi_pi*(FFrot(10,2) + FFrot(10,3)) - delta_delta*(FFrot(10,4) + FFrot(10,5)) + phi_phi*(FFrot(10,6) + FFrot(10,7));
        overlapdR2(irw + 3,icl + 2) = -sigma_sigma*FFrot(15,1) + pi_pi*(FFrot(15,2) + FFrot(15,3)) - delta_delta*(FFrot(15,4) + FFrot(15,5)) + phi_phi*(FFrot(15,6) + FFrot(15,7));
        overlapdR2(irw + 3,icl + 3) = -sigma_sigma*FFrot(19,1) + pi_pi*(FFrot(19,2) + FFrot(19,3)) - delta_delta*(FFrot(19,4) + FFrot(19,5)) + phi_phi*(FFrot(19,6) + FFrot(19,7));
        overlapdR2(irw + 3,icl + 4) = -sigma_sigma*FFrot(20,1) + pi_pi*(FFrot(20,2) + FFrot(20,3)) - delta_delta*(FFrot(20,4) + FFrot(20,5)) + phi_phi*(FFrot(20,6) + FFrot(20,7));
        overlapdR2(irw + 3,icl + 5) = -sigma_sigma*FFrot(21,1) + pi_pi*(FFrot(21,2) + FFrot(21,3)) - delta_delta*(FFrot(21,4) + FFrot(21,5)) + phi_phi*(FFrot(21,6) + FFrot(21,7));
        overlapdR2(irw + 3,icl + 6) = -sigma_sigma*FFrot(22,1) + pi_pi*(FFrot(22,2) + FFrot(22,3)) - delta_delta*(FFrot(22,4) + FFrot(22,5)) + phi_phi*(FFrot(22,6) + FFrot(22,7));
        overlapdR2(irw + 4,icl + 0) = -sigma_sigma*FFrot(5,1) +  pi_pi*(FFrot(5,2) + FFrot(5,3)) - delta_delta*(FFrot(5,4) + FFrot(5,5)) + phi_phi*(FFrot(5,6) + FFrot(5,7));
        overlapdR2(irw + 4,icl + 1) = -sigma_sigma*FFrot(11,1) + pi_pi*(FFrot(11,2) + FFrot(11,3)) - delta_delta*(FFrot(11,4) + FFrot(11,5)) + phi_phi*(FFrot(11,6) + FFrot(11,7));
        overlapdR2(irw + 4,icl + 2) = -sigma_sigma*FFrot(16,1) + pi_pi*(FFrot(16,2) + FFrot(16,3)) - delta_delta*(FFrot(16,4) + FFrot(16,5)) + phi_phi*(FFrot(16,6) + FFrot(16,7));
        overlapdR2(irw + 4,icl + 3) = -sigma_sigma*FFrot(20,1) + pi_pi*(FFrot(20,2) + FFrot(20,3)) - delta_delta*(FFrot(20,4) + FFrot(20,5)) + phi_phi*(FFrot(20,6) + FFrot(20,7));
        overlapdR2(irw + 4,icl + 4) = -sigma_sigma*FFrot(23,1) + pi_pi*(FFrot(23,2) + FFrot(23,3)) - delta_delta*(FFrot(23,4) + FFrot(23,5)) + phi_phi*(FFrot(23,6) + FFrot(23,7));
        overlapdR2(irw + 4,icl + 5) = -sigma_sigma*FFrot(24,1) + pi_pi*(FFrot(24,2) + FFrot(24,3)) - delta_delta*(FFrot(24,4) + FFrot(24,5)) + phi_phi*(FFrot(24,6) + FFrot(24,7));
        overlapdR2(irw + 4,icl + 6) = -sigma_sigma*FFrot(25,1) + pi_pi*(FFrot(25,2) + FFrot(25,3)) - delta_delta*(FFrot(25,4) + FFrot(25,5)) + phi_phi*(FFrot(25,6) + FFrot(25,7));
        overlapdR2(irw + 5,icl + 0) = -sigma_sigma*FFrot(6,1) +  pi_pi*(FFrot(6,2) + FFrot(6,3)) - delta_delta*(FFrot(6,4) + FFrot(6,5)) + phi_phi*(FFrot(6,6) + FFrot(6,7));
        overlapdR2(irw + 5,icl + 1) = -sigma_sigma*FFrot(12,1) + pi_pi*(FFrot(12,2) + FFrot(12,3)) - delta_delta*(FFrot(12,4) + FFrot(12,5)) + phi_phi*(FFrot(12,6) + FFrot(12,7));
        overlapdR2(irw + 5,icl + 2) = -sigma_sigma*FFrot(17,1) + pi_pi*(FFrot(17,2) + FFrot(17,3)) - delta_delta*(FFrot(17,4) + FFrot(17,5)) + phi_phi*(FFrot(17,6) + FFrot(17,7));
        overlapdR2(irw + 5,icl + 3) = -sigma_sigma*FFrot(21,1) + pi_pi*(FFrot(21,2) + FFrot(21,3)) - delta_delta*(FFrot(21,4) + FFrot(21,5)) + phi_phi*(FFrot(21,6) + FFrot(21,7));
        overlapdR2(irw + 5,icl + 4) = -sigma_sigma*FFrot(24,1) + pi_pi*(FFrot(24,2) + FFrot(24,3)) - delta_delta*(FFrot(24,4) + FFrot(24,5)) + phi_phi*(FFrot(24,6) + FFrot(24,7));
        overlapdR2(irw + 5,icl + 5) = -sigma_sigma*FFrot(26,1) + pi_pi*(FFrot(26,2) + FFrot(26,3)) - delta_delta*(FFrot(26,4) + FFrot(26,5)) + phi_phi*(FFrot(26,6) + FFrot(26,7));
        overlapdR2(irw + 5,icl + 6) = -sigma_sigma*FFrot(27,1) + pi_pi*(FFrot(27,2) + FFrot(27,3)) - delta_delta*(FFrot(27,4) + FFrot(27,5)) + phi_phi*(FFrot(27,6) + FFrot(27,7));
        overlapdR2(irw + 6,icl + 0) = -sigma_sigma*FFrot(7,1) +  pi_pi*(FFrot(7,2) + FFrot(7,3)) - delta_delta*(FFrot(7,4) + FFrot(7,5)) + phi_phi*(FFrot(7,6) + FFrot(7,7));
        overlapdR2(irw + 6,icl + 1) = -sigma_sigma*FFrot(13,1) + pi_pi*(FFrot(13,2) + FFrot(13,3)) - delta_delta*(FFrot(13,4) + FFrot(13,5)) + phi_phi*(FFrot(13,6) + FFrot(13,7));
        overlapdR2(irw + 6,icl + 2) = -sigma_sigma*FFrot(18,1) + pi_pi*(FFrot(18,2) + FFrot(18,3)) - delta_delta*(FFrot(18,4) + FFrot(18,5)) + phi_phi*(FFrot(18,6) + FFrot(18,7));
        overlapdR2(irw + 6,icl + 3) = -sigma_sigma*FFrot(22,1) + pi_pi*(FFrot(22,2) + FFrot(22,3)) - delta_delta*(FFrot(22,4) + FFrot(22,5)) + phi_phi*(FFrot(22,6) + FFrot(22,7));
        overlapdR2(irw + 6,icl + 4) = -sigma_sigma*FFrot(25,1) + pi_pi*(FFrot(25,2) + FFrot(25,3)) - delta_delta*(FFrot(25,4) + FFrot(25,5)) + phi_phi*(FFrot(25,6) + FFrot(25,7));
        overlapdR2(irw + 6,icl + 5) = -sigma_sigma*FFrot(27,1) + pi_pi*(FFrot(27,2) + FFrot(27,3)) - delta_delta*(FFrot(27,4) + FFrot(27,5)) + phi_phi*(FFrot(27,6) + FFrot(27,7));
        overlapdR2(irw + 6,icl + 6) = -sigma_sigma*FFrot(28,1) + pi_pi*(FFrot(28,2) + FFrot(28,3)) - delta_delta*(FFrot(28,4) + FFrot(28,5)) + phi_phi*(FFrot(28,6) + FFrot(28,7));
      }
      else if ((iLA == 0)&&(iLB == 1)) {
        //sp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SPOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //partition over 3 dimensions (px,py,pz)
        overlapdR2(irw,icl + 0) = -SProt(1,1)*sigma_sigma;
        overlapdR2(irw,icl + 1) = -SProt(2,1)*sigma_sigma;
        overlapdR2(irw,icl + 2) = -SProt(3,1)*sigma_sigma;
      }
      else if ((iLA == 1)&&(iLB == 0)) {
        //ps
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //partition over 3 dimensions (px,py,pz)
        overlapdR2(irw + 0,icl) = SProt(1,1)*sigma_sigma;
        overlapdR2(irw + 1,icl) = SProt(2,1)*sigma_sigma;
        overlapdR2(irw + 2,icl) = SProt(3,1)*sigma_sigma;
      }
      else if ((iLA == 0)&&(iLB == 2)) {
        //sd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SDOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          overlapdR2(irw,icl + index) = SDrot(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 2)&&(iLB == 0)) {
        //ds
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          overlapdR2(irw + index,icl) = SDrot(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 0)&&(iLB == 3)) {
        //sf
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SFOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //partition over 7 dimensions (fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 7; ++index) {
          overlapdR2(irw,icl + index) = -SFrot(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 3)&&(iLB == 0)) {
        //fs
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 7; ++index) {
          overlapdR2(irw + index,icl) = SFrot(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 1)&&(iLB == 2)) {
        //pd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = PDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          pi_pi -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          pi_pi -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          pi_pi += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          pi_pi += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          pi_pi += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //partition over 3x5 dimensions (px,py,pz)(xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 3; ++index2) {
            overlapdR2(irw + index2,icl + index) = PDrot(1+index2+3*index,1)*sigma_sigma - (PDrot(1+index2+3*index,2) + PDrot(1+index2+3*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 2)&&(iLB == 1)) {
        //dp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = DPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          pi_pi -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          pi_pi -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          pi_pi += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          pi_pi += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          pi_pi += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //partition over 3x5 dimensions (px,py,pz)(xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 3; ++index2) {
            overlapdR2(irw + index,icl + index2) = -PDrot(1+index2+3*index,1)*sigma_sigma + (PDrot(1+index2+3*index,2) + PDrot(1+index2+3*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 1)&&(iLB == 3)) {
        //pf
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = PFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          pi_pi -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          pi_pi -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          pi_pi += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          pi_pi += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          pi_pi += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //partition over 3x5 dimensions (px,py,pz)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 3; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlapdR2(irw + index,icl + index2) = -PFrot(1+index2+7*index,1)*sigma_sigma + (PFrot(1+index2+7*index,2) + PFrot(1+index2+7*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 3)&&(iLB == 1)) {
        //fp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = FPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          pi_pi -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          pi_pi -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          pi_pi += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          pi_pi += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          pi_pi += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //partition over 3x5 dimensions (px,py,pz)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 3; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlapdR2(irw + index2,icl + index) = -PFrot(1+index2+7*index,1)*sigma_sigma + (PFrot(1+index2+7*index,2) + PFrot(1+index2+7*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 2)&&(iLB == 3)) {
        //df
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = DFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          pi_pi -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          pi_pi -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          pi_pi += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          pi_pi += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          pi_pi += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist)/dist;
        delta_delta = 0.0;
        polynomial = DFdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          delta_delta -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          delta_delta -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          delta_delta += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          delta_delta += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          delta_delta += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //partition over 3x5 dimensions (xz,yz,z2,xy,x2-y2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlapdR2(irw + index,icl + index2) = -DFrot(1+index2+7*index,1)*sigma_sigma + (DFrot(1+index2+7*index,2) + DFrot(1+index2+7*index,3))*pi_pi - (DFrot(1+index2+7*index,4) + DFrot(1+index2+7*index,5))*delta_delta;
          }
        }
      }
      else if ((iLA == 3)&&(iLB == 2)) {
        //fd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          sigma_sigma += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          sigma_sigma += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = FDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          pi_pi -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          pi_pi -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          pi_pi += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          pi_pi += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          pi_pi += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist)/dist;
        delta_delta = 0.0;
        polynomial = FDdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double((iNANB + 3)*(iNANB + 2))*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]]/dist;
          delta_delta -= 2.0*aux*argA*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]]/dist;
          delta_delta -= 2.0*aux*argB*double(iNANB + 3)*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1]/dist;
          delta_delta += 2.0*aux*argA*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx] + 1]/dist;
          delta_delta += aux*argA*argA*double(polynomial[0][idx])*Ak[polynomial[1][idx] + 2]*Bk[polynomial[2][idx]]/dist;
          delta_delta += aux*argB*argB*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 2]/dist;
        }
        //partition over 3x5 dimensions (xz,yz,z2,xy,x2-y2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlapdR2(irw + index2,icl + index) = DFrot(1+index2+7*index,1)*sigma_sigma - (DFrot(1+index2+7*index,2) + DFrot(1+index2+7*index,3))*pi_pi + (DFrot(1+index2+7*index,4) + DFrot(1+index2+7*index,5))*delta_delta;
          }
        }
      }
      icl += 2*QN2(isket + 1,2) + 1;
    }
    irw += 2*QN1(isbra + 1,2) + 1;
  }
}
void TwoCenterOv_dRdt(matrixE & overlapdRdt, size_t nao1, size_t nao2, double dist, std::vector<double> & exponents1, matrix<size_t> & QN1, std::vector<double> & exponents2, matrix<size_t> & QN2, std::vector<double> & orientation, bool rotate = true, double tolerance = 1e-6) {
  //implementation of second-derivative of overlap matrix with respect to internuclear distance and angle theta; only up to shell f
  //calculation of the derivative of the overlap function between atoms A and B (A -> 1; B -> 2) with respect to internuclear distance
  overlapdRdt.resize(nao1,nao2);
  overlapdRdt.zero();
  //auxiliary variables for integrals
  double sigma_sigma = 0.0;
  double pi_pi = 0.0;
  double delta_delta = 0.0;
  double phi_phi = 0.0;
  //angles needed to reorient the integrals according to the position of atoms in the molecule
  double cost = orientation[2];
  double sint = sqrt(1 - cost*cost);
  double cosp = 1.0;
  double sinp = 0.0;
  if (fabs(sint) > tolerance) {
    cosp = orientation[0]/sint;
    sinp = orientation[1]/sint;
  }
  //matrices with overlap rotations
  matrixE SProtdt(1,1);
  matrixE SDrotdt(1,1);
  matrixE SFrotdt(1,1);
  matrixE PProtdt(1,1);
  matrixE PDrotdt(1,1);
  matrixE PFrotdt(1,1);
  matrixE DDrotdt(1,1);
  matrixE DFrotdt(1,1);
  matrixE FFrotdt(1,1);
  //other auxiliaries
  std::vector<std::vector<int> > polynomial;
  double a;
  double b;
  double c;
  double argA;
  double argB;
  std::vector<double> Ak;
  std::vector<double> Bk;
  double aux;
  //determine maximum azimuthal quantum number
  size_t maxL1 = 0;
  size_t maxL2 = 0;
  for (size_t idx = 0; idx < QN1.rows(); ++idx) {
    if (QN1(idx + 1,2) > maxL1) {maxL1 = QN1(idx + 1,2);}
  }
  for (size_t idx = 0; idx < QN2.rows(); ++idx) {
    if (QN2(idx + 1,2) > maxL2) {maxL2 = QN2(idx + 1,2);}
  }
  if ((maxL1 > 0)||(maxL2 > 0)) {
    //in this case we have p orbitals
    SProtdt = SPtransf_dt(cost,sint,cosp,sinp,rotate);               //the orientation vector already contains this transformation, so this would not be actually necessary; however, a 3D vector is basically for free and like this I have debug options
    PProtdt = PPtransf_dt(cost,sint,cosp,sinp,rotate);
  }
  if ((maxL1 > 1)||(maxL2 > 1)) {
    //in this case we have d orbitals
    SDrotdt = SDtransf_dt(cost,sint,cosp,sinp,rotate);
    PDrotdt = PDtransf_dt(cost,sint,cosp,sinp,rotate);
    DDrotdt = DDtransf_dt(cost,sint,cosp,sinp,rotate);
  }
  if ((maxL1 > 2)||(maxL2 > 2)) {
    //in this case we have f orbitals
    SFrotdt = SFtransf_dt(cost,sint,cosp,sinp,rotate);
    PFrotdt = PFtransf_dt(cost,sint,cosp,sinp,rotate);
    DFrotdt = DFtransf_dt(cost,sint,cosp,sinp,rotate);
    FFrotdt = FFtransf_dt(cost,sint,cosp,sinp,rotate);
  }
  size_t iNA;
  size_t iNB;
  size_t iLA;
  size_t iLB;
  size_t irw = 1;
  size_t icl = 1;
  for (size_t isbra = 0; isbra < QN1.rows(); ++isbra) {           //loop over bra shells
    a = exponents1[isbra];
    if (a == 0.0) continue;
    iNA = QN1(isbra + 1,1) - 1;             //note that iNA is subtracted -1 !
    iLA = QN1(isbra + 1,2);
    icl = 1;
    for (size_t isket = 0; isket < QN2.rows(); ++isket) {         //loop over ket shells
      b = exponents2[isket];
      if (b == 0.0) continue;
      iNB = QN2(isket + 1,1) - 1;           //note that iNB is subtracted -1 !
      iLB = QN2(isket + 1,2);
      //getting integrals
      argA = 0.5*dist*(a + b);
      argB = 0.5*dist*(a - b);
      Ak = int_1_infty_xpk_expmyx_dx(argA,iNA + iNB + 3);
      Bk = int_m1_1_xpk_expmyx_dx_Pople(argB,iNA + iNB + 3);
      aux = OvFactor(iNA + 1,iLA,0,a,iNB + 1,iLB,0,b,dist)/dist;
      //different azimuthal quantum number cases
      //ss = 0
      if ((iLA == 1)&&(iLB == 1)) {
        //pp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = PPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 9 dimensions
        overlapdRdt(irw + 0,icl + 0) = -sigma_sigma*PProtdt(1,1) + pi_pi*(PProtdt(1,2) + PProtdt(1,3));      //<x|x>
        overlapdRdt(irw + 0,icl + 1) = -sigma_sigma*PProtdt(2,1) + pi_pi*(PProtdt(2,2) + PProtdt(2,3));      //<x|y>
        overlapdRdt(irw + 0,icl + 2) = -sigma_sigma*PProtdt(3,1) + pi_pi*(PProtdt(3,2) + PProtdt(3,3));      //<x|z>
        overlapdRdt(irw + 1,icl + 0) = -sigma_sigma*PProtdt(2,1) + pi_pi*(PProtdt(2,2) + PProtdt(2,3));      //<y|x>
        overlapdRdt(irw + 1,icl + 1) = -sigma_sigma*PProtdt(4,1) + pi_pi*(PProtdt(4,2) + PProtdt(4,3));      //<y|y>
        overlapdRdt(irw + 1,icl + 2) = -sigma_sigma*PProtdt(5,1) + pi_pi*(PProtdt(5,2) + PProtdt(5,3));      //<y|z>
        overlapdRdt(irw + 2,icl + 0) = -sigma_sigma*PProtdt(3,1) + pi_pi*(PProtdt(3,2) + PProtdt(3,3));      //<z|x>
        overlapdRdt(irw + 2,icl + 1) = -sigma_sigma*PProtdt(5,1) + pi_pi*(PProtdt(5,2) + PProtdt(5,3));      //<z|y>
        overlapdRdt(irw + 2,icl + 2) = -sigma_sigma*PProtdt(6,1) + pi_pi*(PProtdt(6,2) + PProtdt(6,3));      //<z|z>
      }
      else if ((iLA == 2)&&(iLB == 2)) {
        //dd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = DDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist)/dist;
        delta_delta = 0.0;
        polynomial = DDdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 5x5 dimensions (xz,yz,z2,xy,x2-y2)(xz,yz,z2,xy,x2-y2)
        overlapdRdt(irw + 0,icl + 0) = sigma_sigma*DDrotdt(1,1)  - pi_pi*(DDrotdt(1,2)  + DDrotdt(1,3))  + delta_delta*(DDrotdt(1,4)  + DDrotdt(1,5));             //<xz|xz>
        overlapdRdt(irw + 0,icl + 1) = sigma_sigma*DDrotdt(2,1)  - pi_pi*(DDrotdt(2,2)  + DDrotdt(2,3))  + delta_delta*(DDrotdt(2,4)  + DDrotdt(2,5));             //<xz|yz>
        overlapdRdt(irw + 0,icl + 2) = sigma_sigma*DDrotdt(3,1)  - pi_pi*(DDrotdt(3,2)  + DDrotdt(3,3))  + delta_delta*(DDrotdt(3,4)  + DDrotdt(3,5));             //<xz|z2>
        overlapdRdt(irw + 0,icl + 3) = sigma_sigma*DDrotdt(4,1)  - pi_pi*(DDrotdt(4,2)  + DDrotdt(4,3))  + delta_delta*(DDrotdt(4,4)  + DDrotdt(4,5));             //<xz|xy>
        overlapdRdt(irw + 0,icl + 4) = sigma_sigma*DDrotdt(5,1)  - pi_pi*(DDrotdt(5,2)  + DDrotdt(5,3))  + delta_delta*(DDrotdt(5,4)  + DDrotdt(5,5));             //<xz|x2-y2>
        overlapdRdt(irw + 1,icl + 0) = sigma_sigma*DDrotdt(2,1)  - pi_pi*(DDrotdt(2,2)  + DDrotdt(2,3))  + delta_delta*(DDrotdt(2,4)  + DDrotdt(2,5));             //<yz|xz>
        overlapdRdt(irw + 1,icl + 1) = sigma_sigma*DDrotdt(6,1)  - pi_pi*(DDrotdt(6,2)  + DDrotdt(6,3))  + delta_delta*(DDrotdt(6,4)  + DDrotdt(6,5));             //<yz|yz>
        overlapdRdt(irw + 1,icl + 2) = sigma_sigma*DDrotdt(7,1)  - pi_pi*(DDrotdt(7,2)  + DDrotdt(7,3))  + delta_delta*(DDrotdt(7,4)  + DDrotdt(7,5));             //<yz|z2>
        overlapdRdt(irw + 1,icl + 3) = sigma_sigma*DDrotdt(8,1)  - pi_pi*(DDrotdt(8,2)  + DDrotdt(8,3))  + delta_delta*(DDrotdt(8,4)  + DDrotdt(8,5));             //<yz|xy>
        overlapdRdt(irw + 1,icl + 4) = sigma_sigma*DDrotdt(9,1)  - pi_pi*(DDrotdt(9,2)  + DDrotdt(9,3))  + delta_delta*(DDrotdt(9,4)  + DDrotdt(9,5));             //<yz|x2-y2>
        overlapdRdt(irw + 2,icl + 0) = sigma_sigma*DDrotdt(3,1)  - pi_pi*(DDrotdt(3,2)  + DDrotdt(3,3))  + delta_delta*(DDrotdt(3,4)  + DDrotdt(3,5));             //<z2|xz>
        overlapdRdt(irw + 2,icl + 1) = sigma_sigma*DDrotdt(7,1)  - pi_pi*(DDrotdt(7,2)  + DDrotdt(7,3))  + delta_delta*(DDrotdt(7,4)  + DDrotdt(7,5));             //<z2|yz>
        overlapdRdt(irw + 2,icl + 2) = sigma_sigma*DDrotdt(10,1) - pi_pi*(DDrotdt(10,2) + DDrotdt(10,3)) + delta_delta*(DDrotdt(10,4) + DDrotdt(10,5));            //<z2|z2>
        overlapdRdt(irw + 2,icl + 3) = sigma_sigma*DDrotdt(11,1) - pi_pi*(DDrotdt(11,2) + DDrotdt(11,3)) + delta_delta*(DDrotdt(11,4) + DDrotdt(11,5));            //<z2|xy>
        overlapdRdt(irw + 2,icl + 4) = sigma_sigma*DDrotdt(12,1) - pi_pi*(DDrotdt(12,2) + DDrotdt(12,3)) + delta_delta*(DDrotdt(12,4) + DDrotdt(12,5));            //<z2|x2-y2>
        overlapdRdt(irw + 3,icl + 0) = sigma_sigma*DDrotdt(4,1)  - pi_pi*(DDrotdt(4,2)  + DDrotdt(4,3))  + delta_delta*(DDrotdt(4,4)  + DDrotdt(4,5));             //<xy|xz>
        overlapdRdt(irw + 3,icl + 1) = sigma_sigma*DDrotdt(8,1)  - pi_pi*(DDrotdt(8,2)  + DDrotdt(8,3))  + delta_delta*(DDrotdt(8,4)  + DDrotdt(8,5));             //<xy|yz>
        overlapdRdt(irw + 3,icl + 2) = sigma_sigma*DDrotdt(11,1) - pi_pi*(DDrotdt(11,2) + DDrotdt(11,3)) + delta_delta*(DDrotdt(11,4) + DDrotdt(11,5));            //<xy|z2>
        overlapdRdt(irw + 3,icl + 3) = sigma_sigma*DDrotdt(13,1) - pi_pi*(DDrotdt(13,2) + DDrotdt(13,3)) + delta_delta*(DDrotdt(13,4) + DDrotdt(13,5));            //<xy|xy>
        overlapdRdt(irw + 3,icl + 4) = sigma_sigma*DDrotdt(14,1) - pi_pi*(DDrotdt(14,2) + DDrotdt(14,3)) + delta_delta*(DDrotdt(14,4) + DDrotdt(14,5));            //<xy|x2-y2>
        overlapdRdt(irw + 4,icl + 0) = sigma_sigma*DDrotdt(5,1)  - pi_pi*(DDrotdt(5,2)  + DDrotdt(5,3))  + delta_delta*(DDrotdt(5,4)  + DDrotdt(5,5));             //<x2-y2|xz>
        overlapdRdt(irw + 4,icl + 1) = sigma_sigma*DDrotdt(9,1)  - pi_pi*(DDrotdt(9,2)  + DDrotdt(9,3))  + delta_delta*(DDrotdt(9,4)  + DDrotdt(9,5));             //<x2-y2|yz>
        overlapdRdt(irw + 4,icl + 2) = sigma_sigma*DDrotdt(12,1) - pi_pi*(DDrotdt(12,2) + DDrotdt(12,3)) + delta_delta*(DDrotdt(12,4) + DDrotdt(12,5));            //<x2-y2|z2>
        overlapdRdt(irw + 4,icl + 3) = sigma_sigma*DDrotdt(14,1) - pi_pi*(DDrotdt(14,2) + DDrotdt(14,3)) + delta_delta*(DDrotdt(14,4) + DDrotdt(14,5));            //<x2-y2|xy>
        overlapdRdt(irw + 4,icl + 4) = sigma_sigma*DDrotdt(15,1) - pi_pi*(DDrotdt(15,2) + DDrotdt(15,3)) + delta_delta*(DDrotdt(15,4) + DDrotdt(15,5));            //<x2-y2|x2-y2>
      }
      else if ((iLA == 3)&&(iLB == 3)) {
        //ff
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = FFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist)/dist;
        delta_delta = 0.0;
        polynomial = FFdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //phi integral
        aux = OvFactor(iNA + 1,iLA,3,a,iNB + 1,iLB,3,b,dist)/dist;
        phi_phi = 0.0;
        polynomial = FFphiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          phi_phi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          phi_phi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          phi_phi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 7x7 dimensions (fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        overlapdRdt(irw + 0,icl + 0) = -sigma_sigma*FFrotdt(1,1) + pi_pi*(FFrotdt(1,2) + FFrotdt(1,3)) - delta_delta*(FFrotdt(1,4) + FFrotdt(1,5)) + phi_phi*(FFrotdt(1,6) + FFrotdt(1,7));
        overlapdRdt(irw + 0,icl + 1) = -sigma_sigma*FFrotdt(2,1) + pi_pi*(FFrotdt(2,2) + FFrotdt(2,3)) - delta_delta*(FFrotdt(2,4) + FFrotdt(2,5)) + phi_phi*(FFrotdt(2,6) + FFrotdt(2,7));
        overlapdRdt(irw + 0,icl + 2) = -sigma_sigma*FFrotdt(3,1) + pi_pi*(FFrotdt(3,2) + FFrotdt(3,3)) - delta_delta*(FFrotdt(3,4) + FFrotdt(3,5)) + phi_phi*(FFrotdt(3,6) + FFrotdt(3,7));
        overlapdRdt(irw + 0,icl + 3) = -sigma_sigma*FFrotdt(4,1) + pi_pi*(FFrotdt(4,2) + FFrotdt(4,3)) - delta_delta*(FFrotdt(4,4) + FFrotdt(4,5)) + phi_phi*(FFrotdt(4,6) + FFrotdt(4,7));
        overlapdRdt(irw + 0,icl + 4) = -sigma_sigma*FFrotdt(5,1) + pi_pi*(FFrotdt(5,2) + FFrotdt(5,3)) - delta_delta*(FFrotdt(5,4) + FFrotdt(5,5)) + phi_phi*(FFrotdt(5,6) + FFrotdt(5,7));
        overlapdRdt(irw + 0,icl + 5) = -sigma_sigma*FFrotdt(6,1) + pi_pi*(FFrotdt(6,2) + FFrotdt(6,3)) - delta_delta*(FFrotdt(6,4) + FFrotdt(6,5)) + phi_phi*(FFrotdt(6,6) + FFrotdt(6,7));
        overlapdRdt(irw + 0,icl + 6) = -sigma_sigma*FFrotdt(7,1) + pi_pi*(FFrotdt(7,2) + FFrotdt(7,3)) - delta_delta*(FFrotdt(7,4) + FFrotdt(7,5)) + phi_phi*(FFrotdt(7,6) + FFrotdt(7,7));
        overlapdRdt(irw + 1,icl + 0) = -sigma_sigma*FFrotdt(2,1) + pi_pi*(FFrotdt(2,2) + FFrotdt(2,3)) - delta_delta*(FFrotdt(2,4) + FFrotdt(2,5)) + phi_phi*(FFrotdt(2,6) + FFrotdt(2,7));
        overlapdRdt(irw + 1,icl + 1) = -sigma_sigma*FFrotdt(8,1) + pi_pi*(FFrotdt(8,2) + FFrotdt(8,3)) - delta_delta*(FFrotdt(8,4) + FFrotdt(8,5)) + phi_phi*(FFrotdt(8,6) + FFrotdt(8,7));
        overlapdRdt(irw + 1,icl + 2) = -sigma_sigma*FFrotdt(9,1) + pi_pi*(FFrotdt(9,2) + FFrotdt(9,3)) - delta_delta*(FFrotdt(9,4) + FFrotdt(9,5)) + phi_phi*(FFrotdt(9,6) + FFrotdt(9,7));
        overlapdRdt(irw + 1,icl + 3) = -sigma_sigma*FFrotdt(10,1) + pi_pi*(FFrotdt(10,2) + FFrotdt(10,3)) - delta_delta*(FFrotdt(10,4) + FFrotdt(10,5)) + phi_phi*(FFrotdt(10,6) + FFrotdt(10,7));
        overlapdRdt(irw + 1,icl + 4) = -sigma_sigma*FFrotdt(11,1) + pi_pi*(FFrotdt(11,2) + FFrotdt(11,3)) - delta_delta*(FFrotdt(11,4) + FFrotdt(11,5)) + phi_phi*(FFrotdt(11,6) + FFrotdt(11,7));
        overlapdRdt(irw + 1,icl + 5) = -sigma_sigma*FFrotdt(12,1) + pi_pi*(FFrotdt(12,2) + FFrotdt(12,3)) - delta_delta*(FFrotdt(12,4) + FFrotdt(12,5)) + phi_phi*(FFrotdt(12,6) + FFrotdt(12,7));
        overlapdRdt(irw + 1,icl + 6) = -sigma_sigma*FFrotdt(13,1) + pi_pi*(FFrotdt(13,2) + FFrotdt(13,3)) - delta_delta*(FFrotdt(13,4) + FFrotdt(13,5)) + phi_phi*(FFrotdt(13,6) + FFrotdt(13,7));
        overlapdRdt(irw + 2,icl + 0) = -sigma_sigma*FFrotdt(3,1) +  pi_pi*(FFrotdt(3,2) + FFrotdt(3,3)) - delta_delta*(FFrotdt(3,4) + FFrotdt(3,5)) + phi_phi*(FFrotdt(3,6) + FFrotdt(3,7));
        overlapdRdt(irw + 2,icl + 1) = -sigma_sigma*FFrotdt(9,1) +  pi_pi*(FFrotdt(9,2) + FFrotdt(9,3)) - delta_delta*(FFrotdt(9,4) + FFrotdt(9,5)) + phi_phi*(FFrotdt(9,6) + FFrotdt(9,7));
        overlapdRdt(irw + 2,icl + 2) = -sigma_sigma*FFrotdt(14,1) + pi_pi*(FFrotdt(14,2) + FFrotdt(14,3)) - delta_delta*(FFrotdt(14,4) + FFrotdt(14,5)) + phi_phi*(FFrotdt(14,6) + FFrotdt(14,7));
        overlapdRdt(irw + 2,icl + 3) = -sigma_sigma*FFrotdt(15,1) + pi_pi*(FFrotdt(15,2) + FFrotdt(15,3)) - delta_delta*(FFrotdt(15,4) + FFrotdt(15,5)) + phi_phi*(FFrotdt(15,6) + FFrotdt(15,7));
        overlapdRdt(irw + 2,icl + 4) = -sigma_sigma*FFrotdt(16,1) + pi_pi*(FFrotdt(16,2) + FFrotdt(16,3)) - delta_delta*(FFrotdt(16,4) + FFrotdt(16,5)) + phi_phi*(FFrotdt(16,6) + FFrotdt(16,7));
        overlapdRdt(irw + 2,icl + 5) = -sigma_sigma*FFrotdt(17,1) + pi_pi*(FFrotdt(17,2) + FFrotdt(17,3)) - delta_delta*(FFrotdt(17,4) + FFrotdt(17,5)) + phi_phi*(FFrotdt(17,6) + FFrotdt(17,7));
        overlapdRdt(irw + 2,icl + 6) = -sigma_sigma*FFrotdt(18,1) + pi_pi*(FFrotdt(18,2) + FFrotdt(18,3)) - delta_delta*(FFrotdt(18,4) + FFrotdt(18,5)) + phi_phi*(FFrotdt(18,6) + FFrotdt(18,7));
        overlapdRdt(irw + 3,icl + 0) = -sigma_sigma*FFrotdt(4,1) +  pi_pi*(FFrotdt(4,2) + FFrotdt(4,3)) - delta_delta*(FFrotdt(4,4) + FFrotdt(4,5)) + phi_phi*(FFrotdt(4,6) + FFrotdt(4,7));
        overlapdRdt(irw + 3,icl + 1) = -sigma_sigma*FFrotdt(10,1) + pi_pi*(FFrotdt(10,2) + FFrotdt(10,3)) - delta_delta*(FFrotdt(10,4) + FFrotdt(10,5)) + phi_phi*(FFrotdt(10,6) + FFrotdt(10,7));
        overlapdRdt(irw + 3,icl + 2) = -sigma_sigma*FFrotdt(15,1) + pi_pi*(FFrotdt(15,2) + FFrotdt(15,3)) - delta_delta*(FFrotdt(15,4) + FFrotdt(15,5)) + phi_phi*(FFrotdt(15,6) + FFrotdt(15,7));
        overlapdRdt(irw + 3,icl + 3) = -sigma_sigma*FFrotdt(19,1) + pi_pi*(FFrotdt(19,2) + FFrotdt(19,3)) - delta_delta*(FFrotdt(19,4) + FFrotdt(19,5)) + phi_phi*(FFrotdt(19,6) + FFrotdt(19,7));
        overlapdRdt(irw + 3,icl + 4) = -sigma_sigma*FFrotdt(20,1) + pi_pi*(FFrotdt(20,2) + FFrotdt(20,3)) - delta_delta*(FFrotdt(20,4) + FFrotdt(20,5)) + phi_phi*(FFrotdt(20,6) + FFrotdt(20,7));
        overlapdRdt(irw + 3,icl + 5) = -sigma_sigma*FFrotdt(21,1) + pi_pi*(FFrotdt(21,2) + FFrotdt(21,3)) - delta_delta*(FFrotdt(21,4) + FFrotdt(21,5)) + phi_phi*(FFrotdt(21,6) + FFrotdt(21,7));
        overlapdRdt(irw + 3,icl + 6) = -sigma_sigma*FFrotdt(22,1) + pi_pi*(FFrotdt(22,2) + FFrotdt(22,3)) - delta_delta*(FFrotdt(22,4) + FFrotdt(22,5)) + phi_phi*(FFrotdt(22,6) + FFrotdt(22,7));
        overlapdRdt(irw + 4,icl + 0) = -sigma_sigma*FFrotdt(5,1) +  pi_pi*(FFrotdt(5,2) + FFrotdt(5,3)) - delta_delta*(FFrotdt(5,4) + FFrotdt(5,5)) + phi_phi*(FFrotdt(5,6) + FFrotdt(5,7));
        overlapdRdt(irw + 4,icl + 1) = -sigma_sigma*FFrotdt(11,1) + pi_pi*(FFrotdt(11,2) + FFrotdt(11,3)) - delta_delta*(FFrotdt(11,4) + FFrotdt(11,5)) + phi_phi*(FFrotdt(11,6) + FFrotdt(11,7));
        overlapdRdt(irw + 4,icl + 2) = -sigma_sigma*FFrotdt(16,1) + pi_pi*(FFrotdt(16,2) + FFrotdt(16,3)) - delta_delta*(FFrotdt(16,4) + FFrotdt(16,5)) + phi_phi*(FFrotdt(16,6) + FFrotdt(16,7));
        overlapdRdt(irw + 4,icl + 3) = -sigma_sigma*FFrotdt(20,1) + pi_pi*(FFrotdt(20,2) + FFrotdt(20,3)) - delta_delta*(FFrotdt(20,4) + FFrotdt(20,5)) + phi_phi*(FFrotdt(20,6) + FFrotdt(20,7));
        overlapdRdt(irw + 4,icl + 4) = -sigma_sigma*FFrotdt(23,1) + pi_pi*(FFrotdt(23,2) + FFrotdt(23,3)) - delta_delta*(FFrotdt(23,4) + FFrotdt(23,5)) + phi_phi*(FFrotdt(23,6) + FFrotdt(23,7));
        overlapdRdt(irw + 4,icl + 5) = -sigma_sigma*FFrotdt(24,1) + pi_pi*(FFrotdt(24,2) + FFrotdt(24,3)) - delta_delta*(FFrotdt(24,4) + FFrotdt(24,5)) + phi_phi*(FFrotdt(24,6) + FFrotdt(24,7));
        overlapdRdt(irw + 4,icl + 6) = -sigma_sigma*FFrotdt(25,1) + pi_pi*(FFrotdt(25,2) + FFrotdt(25,3)) - delta_delta*(FFrotdt(25,4) + FFrotdt(25,5)) + phi_phi*(FFrotdt(25,6) + FFrotdt(25,7));
        overlapdRdt(irw + 5,icl + 0) = -sigma_sigma*FFrotdt(6,1) +  pi_pi*(FFrotdt(6,2) + FFrotdt(6,3)) - delta_delta*(FFrotdt(6,4) + FFrotdt(6,5)) + phi_phi*(FFrotdt(6,6) + FFrotdt(6,7));
        overlapdRdt(irw + 5,icl + 1) = -sigma_sigma*FFrotdt(12,1) + pi_pi*(FFrotdt(12,2) + FFrotdt(12,3)) - delta_delta*(FFrotdt(12,4) + FFrotdt(12,5)) + phi_phi*(FFrotdt(12,6) + FFrotdt(12,7));
        overlapdRdt(irw + 5,icl + 2) = -sigma_sigma*FFrotdt(17,1) + pi_pi*(FFrotdt(17,2) + FFrotdt(17,3)) - delta_delta*(FFrotdt(17,4) + FFrotdt(17,5)) + phi_phi*(FFrotdt(17,6) + FFrotdt(17,7));
        overlapdRdt(irw + 5,icl + 3) = -sigma_sigma*FFrotdt(21,1) + pi_pi*(FFrotdt(21,2) + FFrotdt(21,3)) - delta_delta*(FFrotdt(21,4) + FFrotdt(21,5)) + phi_phi*(FFrotdt(21,6) + FFrotdt(21,7));
        overlapdRdt(irw + 5,icl + 4) = -sigma_sigma*FFrotdt(24,1) + pi_pi*(FFrotdt(24,2) + FFrotdt(24,3)) - delta_delta*(FFrotdt(24,4) + FFrotdt(24,5)) + phi_phi*(FFrotdt(24,6) + FFrotdt(24,7));
        overlapdRdt(irw + 5,icl + 5) = -sigma_sigma*FFrotdt(26,1) + pi_pi*(FFrotdt(26,2) + FFrotdt(26,3)) - delta_delta*(FFrotdt(26,4) + FFrotdt(26,5)) + phi_phi*(FFrotdt(26,6) + FFrotdt(26,7));
        overlapdRdt(irw + 5,icl + 6) = -sigma_sigma*FFrotdt(27,1) + pi_pi*(FFrotdt(27,2) + FFrotdt(27,3)) - delta_delta*(FFrotdt(27,4) + FFrotdt(27,5)) + phi_phi*(FFrotdt(27,6) + FFrotdt(27,7));
        overlapdRdt(irw + 6,icl + 0) = -sigma_sigma*FFrotdt(7,1) +  pi_pi*(FFrotdt(7,2) + FFrotdt(7,3)) - delta_delta*(FFrotdt(7,4) + FFrotdt(7,5)) + phi_phi*(FFrotdt(7,6) + FFrotdt(7,7));
        overlapdRdt(irw + 6,icl + 1) = -sigma_sigma*FFrotdt(13,1) + pi_pi*(FFrotdt(13,2) + FFrotdt(13,3)) - delta_delta*(FFrotdt(13,4) + FFrotdt(13,5)) + phi_phi*(FFrotdt(13,6) + FFrotdt(13,7));
        overlapdRdt(irw + 6,icl + 2) = -sigma_sigma*FFrotdt(18,1) + pi_pi*(FFrotdt(18,2) + FFrotdt(18,3)) - delta_delta*(FFrotdt(18,4) + FFrotdt(18,5)) + phi_phi*(FFrotdt(18,6) + FFrotdt(18,7));
        overlapdRdt(irw + 6,icl + 3) = -sigma_sigma*FFrotdt(22,1) + pi_pi*(FFrotdt(22,2) + FFrotdt(22,3)) - delta_delta*(FFrotdt(22,4) + FFrotdt(22,5)) + phi_phi*(FFrotdt(22,6) + FFrotdt(22,7));
        overlapdRdt(irw + 6,icl + 4) = -sigma_sigma*FFrotdt(25,1) + pi_pi*(FFrotdt(25,2) + FFrotdt(25,3)) - delta_delta*(FFrotdt(25,4) + FFrotdt(25,5)) + phi_phi*(FFrotdt(25,6) + FFrotdt(25,7));
        overlapdRdt(irw + 6,icl + 5) = -sigma_sigma*FFrotdt(27,1) + pi_pi*(FFrotdt(27,2) + FFrotdt(27,3)) - delta_delta*(FFrotdt(27,4) + FFrotdt(27,5)) + phi_phi*(FFrotdt(27,6) + FFrotdt(27,7));
        overlapdRdt(irw + 6,icl + 6) = -sigma_sigma*FFrotdt(28,1) + pi_pi*(FFrotdt(28,2) + FFrotdt(28,3)) - delta_delta*(FFrotdt(28,4) + FFrotdt(28,5)) + phi_phi*(FFrotdt(28,6) + FFrotdt(28,7));
      }
      else if ((iLA == 0)&&(iLB == 1)) {
        //sp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SPOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3 dimensions (px,py,pz)
        overlapdRdt(irw,icl + 0) = -SProtdt(1,1)*sigma_sigma;
        overlapdRdt(irw,icl + 1) = -SProtdt(2,1)*sigma_sigma;
        overlapdRdt(irw,icl + 2) = -SProtdt(3,1)*sigma_sigma;
      }
      else if ((iLA == 1)&&(iLB == 0)) {
        //ps
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3 dimensions (px,py,pz)
        overlapdRdt(irw + 0,icl) = SProtdt(1,1)*sigma_sigma;
        overlapdRdt(irw + 1,icl) = SProtdt(2,1)*sigma_sigma;
        overlapdRdt(irw + 2,icl) = SProtdt(3,1)*sigma_sigma;
      }
      else if ((iLA == 0)&&(iLB == 2)) {
        //sd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SDOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          overlapdRdt(irw,icl + index) = SDrotdt(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 2)&&(iLB == 0)) {
        //ds
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          overlapdRdt(irw + index,icl) = SDrotdt(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 0)&&(iLB == 3)) {
        //sf
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SFOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 7 dimensions (fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 7; ++index) {
          overlapdRdt(irw,icl + index) = -SFrotdt(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 3)&&(iLB == 0)) {
        //fs
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 7; ++index) {
          overlapdRdt(irw + index,icl) = SFrotdt(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 1)&&(iLB == 2)) {
        //pd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = PDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3x5 dimensions (px,py,pz)(xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 3; ++index2) {
            overlapdRdt(irw + index2,icl + index) = PDrotdt(1+index2+3*index,1)*sigma_sigma - (PDrotdt(1+index2+3*index,2) + PDrotdt(1+index2+3*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 2)&&(iLB == 1)) {
        //dp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = DPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3x5 dimensions (px,py,pz)(xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 3; ++index2) {
            overlapdRdt(irw + index,icl + index2) = -PDrotdt(1+index2+3*index,1)*sigma_sigma + (PDrotdt(1+index2+3*index,2) + PDrotdt(1+index2+3*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 1)&&(iLB == 3)) {
        //pf
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = PFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3x5 dimensions (px,py,pz)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 3; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlapdRdt(irw + index,icl + index2) = -PFrotdt(1+index2+7*index,1)*sigma_sigma + (PFrotdt(1+index2+7*index,2) + PFrotdt(1+index2+7*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 3)&&(iLB == 1)) {
        //fp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = FPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3x5 dimensions (px,py,pz)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 3; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlapdRdt(irw + index2,icl + index) = -PFrotdt(1+index2+7*index,1)*sigma_sigma + (PFrotdt(1+index2+7*index,2) + PFrotdt(1+index2+7*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 2)&&(iLB == 3)) {
        //df
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = DFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist)/dist;
        delta_delta = 0.0;
        polynomial = DFdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3x5 dimensions (xz,yz,z2,xy,x2-y2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlapdRdt(irw + index,icl + index2) = -DFrotdt(1+index2+7*index,1)*sigma_sigma + (DFrotdt(1+index2+7*index,2) + DFrotdt(1+index2+7*index,3))*pi_pi - (DFrotdt(1+index2+7*index,4) + DFrotdt(1+index2+7*index,5))*delta_delta;
          }
        }
      }
      else if ((iLA == 3)&&(iLB == 2)) {
        //fd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = FDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist)/dist;
        delta_delta = 0.0;
        polynomial = FDdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3x5 dimensions (xz,yz,z2,xy,x2-y2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlapdRdt(irw + index2,icl + index) = DFrotdt(1+index2+7*index,1)*sigma_sigma - (DFrotdt(1+index2+7*index,2) + DFrotdt(1+index2+7*index,3))*pi_pi + (DFrotdt(1+index2+7*index,4) + DFrotdt(1+index2+7*index,5))*delta_delta;
          }
        }
      }
      icl += 2*QN2(isket + 1,2) + 1;
    }
    irw += 2*QN1(isbra + 1,2) + 1;
  }
}
void TwoCenterOv_dRdp(matrixE & overlapdRdp, size_t nao1, size_t nao2, double dist, std::vector<double> & exponents1, matrix<size_t> & QN1, std::vector<double> & exponents2, matrix<size_t> & QN2, std::vector<double> & orientation, bool rotate = true, double tolerance = 1e-6) {
  //implementation of second-derivative of overlap matrix with respect to internuclear distance and angle phi; only up to shell f
  //calculation of the derivative of the overlap function between atoms A and B (A -> 1; B -> 2) with respect to internuclear distance
  overlapdRdp.resize(nao1,nao2);
  overlapdRdp.zero();
  //auxiliary variables for integrals
  double sigma_sigma = 0.0;
  double pi_pi = 0.0;
  double delta_delta = 0.0;
  double phi_phi = 0.0;
  //angles needed to reorient the integrals according to the position of atoms in the molecule
  double cost = orientation[2];
  double sint = sqrt(1 - cost*cost);
  double cosp = 1.0;
  double sinp = 0.0;
  if (fabs(sint) > tolerance) {
    cosp = orientation[0]/sint;
    sinp = orientation[1]/sint;
  }
  //matrices with overlap rotations
  matrixE SProtdp(1,1);
  matrixE SDrotdp(1,1);
  matrixE SFrotdp(1,1);
  matrixE PProtdp(1,1);
  matrixE PDrotdp(1,1);
  matrixE PFrotdp(1,1);
  matrixE DDrotdp(1,1);
  matrixE DFrotdp(1,1);
  matrixE FFrotdp(1,1);
  //other auxiliaries
  std::vector<std::vector<int> > polynomial;
  double a;
  double b;
  double c;
  double argA;
  double argB;
  std::vector<double> Ak;
  std::vector<double> Bk;
  double aux;
  //determine maximum azimuthal quantum number
  size_t maxL1 = 0;
  size_t maxL2 = 0;
  for (size_t idx = 0; idx < QN1.rows(); ++idx) {
    if (QN1(idx + 1,2) > maxL1) {maxL1 = QN1(idx + 1,2);}
  }
  for (size_t idx = 0; idx < QN2.rows(); ++idx) {
    if (QN2(idx + 1,2) > maxL2) {maxL2 = QN2(idx + 1,2);}
  }
  if ((maxL1 > 0)||(maxL2 > 0)) {
    //in this case we have p orbitals
    SProtdp = SPtransf_dp(cost,sint,cosp,sinp,rotate);               //the orientation vector already contains this transformation, so this would not be actually necessary; however, a 3D vector is basically for free and like this I have debug options
    PProtdp = PPtransf_dp(cost,sint,cosp,sinp,rotate);
  }
  if ((maxL1 > 1)||(maxL2 > 1)) {
    //in this case we have d orbitals
    SDrotdp = SDtransf_dp(cost,sint,cosp,sinp,rotate);
    PDrotdp = PDtransf_dp(cost,sint,cosp,sinp,rotate);
    DDrotdp = DDtransf_dp(cost,sint,cosp,sinp,rotate);
  }
  if ((maxL1 > 2)||(maxL2 > 2)) {
    //in this case we have f orbitals
    SFrotdp = SFtransf_dp(cost,sint,cosp,sinp,rotate);
    PFrotdp = PFtransf_dp(cost,sint,cosp,sinp,rotate);
    DFrotdp = DFtransf_dp(cost,sint,cosp,sinp,rotate);
    FFrotdp = FFtransf_dp(cost,sint,cosp,sinp,rotate);
  }
  size_t iNA;
  size_t iNB;
  size_t iLA;
  size_t iLB;
  size_t irw = 1;
  size_t icl = 1;
  for (size_t isbra = 0; isbra < QN1.rows(); ++isbra) {           //loop over bra shells
    a = exponents1[isbra];
    if (a == 0.0) continue;
    iNA = QN1(isbra + 1,1) - 1;             //note that iNA is subtracted -1 !
    iLA = QN1(isbra + 1,2);
    icl = 1;
    for (size_t isket = 0; isket < QN2.rows(); ++isket) {         //loop over ket shells
      b = exponents2[isket];
      if (b == 0.0) continue;
      iNB = QN2(isket + 1,1) - 1;           //note that iNB is subtracted -1 !
      iLB = QN2(isket + 1,2);
      //getting integrals
      argA = 0.5*dist*(a + b);
      argB = 0.5*dist*(a - b);
      Ak = int_1_infty_xpk_expmyx_dx(argA,iNA + iNB + 3);
      Bk = int_m1_1_xpk_expmyx_dx_Pople(argB,iNA + iNB + 3);
      aux = OvFactor(iNA + 1,iLA,0,a,iNB + 1,iLB,0,b,dist)/dist;
      //different azimuthal quantum number cases
      //ss = 0
      if ((iLA == 1)&&(iLB == 1)) {
        //pp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = PPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 9 dimensions
        overlapdRdp(irw + 0,icl + 0) = -sigma_sigma*PProtdp(1,1) + pi_pi*(PProtdp(1,2) + PProtdp(1,3));      //<x|x>
        overlapdRdp(irw + 0,icl + 1) = -sigma_sigma*PProtdp(2,1) + pi_pi*(PProtdp(2,2) + PProtdp(2,3));      //<x|y>
        overlapdRdp(irw + 0,icl + 2) = -sigma_sigma*PProtdp(3,1) + pi_pi*(PProtdp(3,2) + PProtdp(3,3));      //<x|z>
        overlapdRdp(irw + 1,icl + 0) = -sigma_sigma*PProtdp(2,1) + pi_pi*(PProtdp(2,2) + PProtdp(2,3));      //<y|x>
        overlapdRdp(irw + 1,icl + 1) = -sigma_sigma*PProtdp(4,1) + pi_pi*(PProtdp(4,2) + PProtdp(4,3));      //<y|y>
        overlapdRdp(irw + 1,icl + 2) = -sigma_sigma*PProtdp(5,1) + pi_pi*(PProtdp(5,2) + PProtdp(5,3));      //<y|z>
        overlapdRdp(irw + 2,icl + 0) = -sigma_sigma*PProtdp(3,1) + pi_pi*(PProtdp(3,2) + PProtdp(3,3));      //<z|x>
        overlapdRdp(irw + 2,icl + 1) = -sigma_sigma*PProtdp(5,1) + pi_pi*(PProtdp(5,2) + PProtdp(5,3));      //<z|y>
        overlapdRdp(irw + 2,icl + 2) = -sigma_sigma*PProtdp(6,1) + pi_pi*(PProtdp(6,2) + PProtdp(6,3));      //<z|z>
      }
      else if ((iLA == 2)&&(iLB == 2)) {
        //dd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = DDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist)/dist;
        delta_delta = 0.0;
        polynomial = DDdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 5x5 dimensions (xz,yz,z2,xy,x2-y2)(xz,yz,z2,xy,x2-y2)
        overlapdRdp(irw + 0,icl + 0) = sigma_sigma*DDrotdp(1,1)  - pi_pi*(DDrotdp(1,2)  + DDrotdp(1,3))  + delta_delta*(DDrotdp(1,4)  + DDrotdp(1,5));             //<xz|xz>
        overlapdRdp(irw + 0,icl + 1) = sigma_sigma*DDrotdp(2,1)  - pi_pi*(DDrotdp(2,2)  + DDrotdp(2,3))  + delta_delta*(DDrotdp(2,4)  + DDrotdp(2,5));             //<xz|yz>
        overlapdRdp(irw + 0,icl + 2) = sigma_sigma*DDrotdp(3,1)  - pi_pi*(DDrotdp(3,2)  + DDrotdp(3,3))  + delta_delta*(DDrotdp(3,4)  + DDrotdp(3,5));             //<xz|z2>
        overlapdRdp(irw + 0,icl + 3) = sigma_sigma*DDrotdp(4,1)  - pi_pi*(DDrotdp(4,2)  + DDrotdp(4,3))  + delta_delta*(DDrotdp(4,4)  + DDrotdp(4,5));             //<xz|xy>
        overlapdRdp(irw + 0,icl + 4) = sigma_sigma*DDrotdp(5,1)  - pi_pi*(DDrotdp(5,2)  + DDrotdp(5,3))  + delta_delta*(DDrotdp(5,4)  + DDrotdp(5,5));             //<xz|x2-y2>
        overlapdRdp(irw + 1,icl + 0) = sigma_sigma*DDrotdp(2,1)  - pi_pi*(DDrotdp(2,2)  + DDrotdp(2,3))  + delta_delta*(DDrotdp(2,4)  + DDrotdp(2,5));             //<yz|xz>
        overlapdRdp(irw + 1,icl + 1) = sigma_sigma*DDrotdp(6,1)  - pi_pi*(DDrotdp(6,2)  + DDrotdp(6,3))  + delta_delta*(DDrotdp(6,4)  + DDrotdp(6,5));             //<yz|yz>
        overlapdRdp(irw + 1,icl + 2) = sigma_sigma*DDrotdp(7,1)  - pi_pi*(DDrotdp(7,2)  + DDrotdp(7,3))  + delta_delta*(DDrotdp(7,4)  + DDrotdp(7,5));             //<yz|z2>
        overlapdRdp(irw + 1,icl + 3) = sigma_sigma*DDrotdp(8,1)  - pi_pi*(DDrotdp(8,2)  + DDrotdp(8,3))  + delta_delta*(DDrotdp(8,4)  + DDrotdp(8,5));             //<yz|xy>
        overlapdRdp(irw + 1,icl + 4) = sigma_sigma*DDrotdp(9,1)  - pi_pi*(DDrotdp(9,2)  + DDrotdp(9,3))  + delta_delta*(DDrotdp(9,4)  + DDrotdp(9,5));             //<yz|x2-y2>
        overlapdRdp(irw + 2,icl + 0) = sigma_sigma*DDrotdp(3,1)  - pi_pi*(DDrotdp(3,2)  + DDrotdp(3,3))  + delta_delta*(DDrotdp(3,4)  + DDrotdp(3,5));             //<z2|xz>
        overlapdRdp(irw + 2,icl + 1) = sigma_sigma*DDrotdp(7,1)  - pi_pi*(DDrotdp(7,2)  + DDrotdp(7,3))  + delta_delta*(DDrotdp(7,4)  + DDrotdp(7,5));             //<z2|yz>
        overlapdRdp(irw + 2,icl + 2) = sigma_sigma*DDrotdp(10,1) - pi_pi*(DDrotdp(10,2) + DDrotdp(10,3)) + delta_delta*(DDrotdp(10,4) + DDrotdp(10,5));            //<z2|z2>
        overlapdRdp(irw + 2,icl + 3) = sigma_sigma*DDrotdp(11,1) - pi_pi*(DDrotdp(11,2) + DDrotdp(11,3)) + delta_delta*(DDrotdp(11,4) + DDrotdp(11,5));            //<z2|xy>
        overlapdRdp(irw + 2,icl + 4) = sigma_sigma*DDrotdp(12,1) - pi_pi*(DDrotdp(12,2) + DDrotdp(12,3)) + delta_delta*(DDrotdp(12,4) + DDrotdp(12,5));            //<z2|x2-y2>
        overlapdRdp(irw + 3,icl + 0) = sigma_sigma*DDrotdp(4,1)  - pi_pi*(DDrotdp(4,2)  + DDrotdp(4,3))  + delta_delta*(DDrotdp(4,4)  + DDrotdp(4,5));             //<xy|xz>
        overlapdRdp(irw + 3,icl + 1) = sigma_sigma*DDrotdp(8,1)  - pi_pi*(DDrotdp(8,2)  + DDrotdp(8,3))  + delta_delta*(DDrotdp(8,4)  + DDrotdp(8,5));             //<xy|yz>
        overlapdRdp(irw + 3,icl + 2) = sigma_sigma*DDrotdp(11,1) - pi_pi*(DDrotdp(11,2) + DDrotdp(11,3)) + delta_delta*(DDrotdp(11,4) + DDrotdp(11,5));            //<xy|z2>
        overlapdRdp(irw + 3,icl + 3) = sigma_sigma*DDrotdp(13,1) - pi_pi*(DDrotdp(13,2) + DDrotdp(13,3)) + delta_delta*(DDrotdp(13,4) + DDrotdp(13,5));            //<xy|xy>
        overlapdRdp(irw + 3,icl + 4) = sigma_sigma*DDrotdp(14,1) - pi_pi*(DDrotdp(14,2) + DDrotdp(14,3)) + delta_delta*(DDrotdp(14,4) + DDrotdp(14,5));            //<xy|x2-y2>
        overlapdRdp(irw + 4,icl + 0) = sigma_sigma*DDrotdp(5,1)  - pi_pi*(DDrotdp(5,2)  + DDrotdp(5,3))  + delta_delta*(DDrotdp(5,4)  + DDrotdp(5,5));             //<x2-y2|xz>
        overlapdRdp(irw + 4,icl + 1) = sigma_sigma*DDrotdp(9,1)  - pi_pi*(DDrotdp(9,2)  + DDrotdp(9,3))  + delta_delta*(DDrotdp(9,4)  + DDrotdp(9,5));             //<x2-y2|yz>
        overlapdRdp(irw + 4,icl + 2) = sigma_sigma*DDrotdp(12,1) - pi_pi*(DDrotdp(12,2) + DDrotdp(12,3)) + delta_delta*(DDrotdp(12,4) + DDrotdp(12,5));            //<x2-y2|z2>
        overlapdRdp(irw + 4,icl + 3) = sigma_sigma*DDrotdp(14,1) - pi_pi*(DDrotdp(14,2) + DDrotdp(14,3)) + delta_delta*(DDrotdp(14,4) + DDrotdp(14,5));            //<x2-y2|xy>
        overlapdRdp(irw + 4,icl + 4) = sigma_sigma*DDrotdp(15,1) - pi_pi*(DDrotdp(15,2) + DDrotdp(15,3)) + delta_delta*(DDrotdp(15,4) + DDrotdp(15,5));            //<x2-y2|x2-y2>
      }
      else if ((iLA == 3)&&(iLB == 3)) {
        //ff
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = FFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist)/dist;
        delta_delta = 0.0;
        polynomial = FFdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //phi integral
        aux = OvFactor(iNA + 1,iLA,3,a,iNB + 1,iLB,3,b,dist)/dist;
        phi_phi = 0.0;
        polynomial = FFphiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          phi_phi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          phi_phi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          phi_phi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 7x7 dimensions (fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        overlapdRdp(irw + 0,icl + 0) = -sigma_sigma*FFrotdp(1,1) + pi_pi*(FFrotdp(1,2) + FFrotdp(1,3)) - delta_delta*(FFrotdp(1,4) + FFrotdp(1,5)) + phi_phi*(FFrotdp(1,6) + FFrotdp(1,7));
        overlapdRdp(irw + 0,icl + 1) = -sigma_sigma*FFrotdp(2,1) + pi_pi*(FFrotdp(2,2) + FFrotdp(2,3)) - delta_delta*(FFrotdp(2,4) + FFrotdp(2,5)) + phi_phi*(FFrotdp(2,6) + FFrotdp(2,7));
        overlapdRdp(irw + 0,icl + 2) = -sigma_sigma*FFrotdp(3,1) + pi_pi*(FFrotdp(3,2) + FFrotdp(3,3)) - delta_delta*(FFrotdp(3,4) + FFrotdp(3,5)) + phi_phi*(FFrotdp(3,6) + FFrotdp(3,7));
        overlapdRdp(irw + 0,icl + 3) = -sigma_sigma*FFrotdp(4,1) + pi_pi*(FFrotdp(4,2) + FFrotdp(4,3)) - delta_delta*(FFrotdp(4,4) + FFrotdp(4,5)) + phi_phi*(FFrotdp(4,6) + FFrotdp(4,7));
        overlapdRdp(irw + 0,icl + 4) = -sigma_sigma*FFrotdp(5,1) + pi_pi*(FFrotdp(5,2) + FFrotdp(5,3)) - delta_delta*(FFrotdp(5,4) + FFrotdp(5,5)) + phi_phi*(FFrotdp(5,6) + FFrotdp(5,7));
        overlapdRdp(irw + 0,icl + 5) = -sigma_sigma*FFrotdp(6,1) + pi_pi*(FFrotdp(6,2) + FFrotdp(6,3)) - delta_delta*(FFrotdp(6,4) + FFrotdp(6,5)) + phi_phi*(FFrotdp(6,6) + FFrotdp(6,7));
        overlapdRdp(irw + 0,icl + 6) = -sigma_sigma*FFrotdp(7,1) + pi_pi*(FFrotdp(7,2) + FFrotdp(7,3)) - delta_delta*(FFrotdp(7,4) + FFrotdp(7,5)) + phi_phi*(FFrotdp(7,6) + FFrotdp(7,7));
        overlapdRdp(irw + 1,icl + 0) = -sigma_sigma*FFrotdp(2,1) + pi_pi*(FFrotdp(2,2) + FFrotdp(2,3)) - delta_delta*(FFrotdp(2,4) + FFrotdp(2,5)) + phi_phi*(FFrotdp(2,6) + FFrotdp(2,7));
        overlapdRdp(irw + 1,icl + 1) = -sigma_sigma*FFrotdp(8,1) + pi_pi*(FFrotdp(8,2) + FFrotdp(8,3)) - delta_delta*(FFrotdp(8,4) + FFrotdp(8,5)) + phi_phi*(FFrotdp(8,6) + FFrotdp(8,7));
        overlapdRdp(irw + 1,icl + 2) = -sigma_sigma*FFrotdp(9,1) + pi_pi*(FFrotdp(9,2) + FFrotdp(9,3)) - delta_delta*(FFrotdp(9,4) + FFrotdp(9,5)) + phi_phi*(FFrotdp(9,6) + FFrotdp(9,7));
        overlapdRdp(irw + 1,icl + 3) = -sigma_sigma*FFrotdp(10,1) + pi_pi*(FFrotdp(10,2) + FFrotdp(10,3)) - delta_delta*(FFrotdp(10,4) + FFrotdp(10,5)) + phi_phi*(FFrotdp(10,6) + FFrotdp(10,7));
        overlapdRdp(irw + 1,icl + 4) = -sigma_sigma*FFrotdp(11,1) + pi_pi*(FFrotdp(11,2) + FFrotdp(11,3)) - delta_delta*(FFrotdp(11,4) + FFrotdp(11,5)) + phi_phi*(FFrotdp(11,6) + FFrotdp(11,7));
        overlapdRdp(irw + 1,icl + 5) = -sigma_sigma*FFrotdp(12,1) + pi_pi*(FFrotdp(12,2) + FFrotdp(12,3)) - delta_delta*(FFrotdp(12,4) + FFrotdp(12,5)) + phi_phi*(FFrotdp(12,6) + FFrotdp(12,7));
        overlapdRdp(irw + 1,icl + 6) = -sigma_sigma*FFrotdp(13,1) + pi_pi*(FFrotdp(13,2) + FFrotdp(13,3)) - delta_delta*(FFrotdp(13,4) + FFrotdp(13,5)) + phi_phi*(FFrotdp(13,6) + FFrotdp(13,7));
        overlapdRdp(irw + 2,icl + 0) = -sigma_sigma*FFrotdp(3,1) +  pi_pi*(FFrotdp(3,2) + FFrotdp(3,3)) - delta_delta*(FFrotdp(3,4) + FFrotdp(3,5)) + phi_phi*(FFrotdp(3,6) + FFrotdp(3,7));
        overlapdRdp(irw + 2,icl + 1) = -sigma_sigma*FFrotdp(9,1) +  pi_pi*(FFrotdp(9,2) + FFrotdp(9,3)) - delta_delta*(FFrotdp(9,4) + FFrotdp(9,5)) + phi_phi*(FFrotdp(9,6) + FFrotdp(9,7));
        overlapdRdp(irw + 2,icl + 2) = -sigma_sigma*FFrotdp(14,1) + pi_pi*(FFrotdp(14,2) + FFrotdp(14,3)) - delta_delta*(FFrotdp(14,4) + FFrotdp(14,5)) + phi_phi*(FFrotdp(14,6) + FFrotdp(14,7));
        overlapdRdp(irw + 2,icl + 3) = -sigma_sigma*FFrotdp(15,1) + pi_pi*(FFrotdp(15,2) + FFrotdp(15,3)) - delta_delta*(FFrotdp(15,4) + FFrotdp(15,5)) + phi_phi*(FFrotdp(15,6) + FFrotdp(15,7));
        overlapdRdp(irw + 2,icl + 4) = -sigma_sigma*FFrotdp(16,1) + pi_pi*(FFrotdp(16,2) + FFrotdp(16,3)) - delta_delta*(FFrotdp(16,4) + FFrotdp(16,5)) + phi_phi*(FFrotdp(16,6) + FFrotdp(16,7));
        overlapdRdp(irw + 2,icl + 5) = -sigma_sigma*FFrotdp(17,1) + pi_pi*(FFrotdp(17,2) + FFrotdp(17,3)) - delta_delta*(FFrotdp(17,4) + FFrotdp(17,5)) + phi_phi*(FFrotdp(17,6) + FFrotdp(17,7));
        overlapdRdp(irw + 2,icl + 6) = -sigma_sigma*FFrotdp(18,1) + pi_pi*(FFrotdp(18,2) + FFrotdp(18,3)) - delta_delta*(FFrotdp(18,4) + FFrotdp(18,5)) + phi_phi*(FFrotdp(18,6) + FFrotdp(18,7));
        overlapdRdp(irw + 3,icl + 0) = -sigma_sigma*FFrotdp(4,1) +  pi_pi*(FFrotdp(4,2) + FFrotdp(4,3)) - delta_delta*(FFrotdp(4,4) + FFrotdp(4,5)) + phi_phi*(FFrotdp(4,6) + FFrotdp(4,7));
        overlapdRdp(irw + 3,icl + 1) = -sigma_sigma*FFrotdp(10,1) + pi_pi*(FFrotdp(10,2) + FFrotdp(10,3)) - delta_delta*(FFrotdp(10,4) + FFrotdp(10,5)) + phi_phi*(FFrotdp(10,6) + FFrotdp(10,7));
        overlapdRdp(irw + 3,icl + 2) = -sigma_sigma*FFrotdp(15,1) + pi_pi*(FFrotdp(15,2) + FFrotdp(15,3)) - delta_delta*(FFrotdp(15,4) + FFrotdp(15,5)) + phi_phi*(FFrotdp(15,6) + FFrotdp(15,7));
        overlapdRdp(irw + 3,icl + 3) = -sigma_sigma*FFrotdp(19,1) + pi_pi*(FFrotdp(19,2) + FFrotdp(19,3)) - delta_delta*(FFrotdp(19,4) + FFrotdp(19,5)) + phi_phi*(FFrotdp(19,6) + FFrotdp(19,7));
        overlapdRdp(irw + 3,icl + 4) = -sigma_sigma*FFrotdp(20,1) + pi_pi*(FFrotdp(20,2) + FFrotdp(20,3)) - delta_delta*(FFrotdp(20,4) + FFrotdp(20,5)) + phi_phi*(FFrotdp(20,6) + FFrotdp(20,7));
        overlapdRdp(irw + 3,icl + 5) = -sigma_sigma*FFrotdp(21,1) + pi_pi*(FFrotdp(21,2) + FFrotdp(21,3)) - delta_delta*(FFrotdp(21,4) + FFrotdp(21,5)) + phi_phi*(FFrotdp(21,6) + FFrotdp(21,7));
        overlapdRdp(irw + 3,icl + 6) = -sigma_sigma*FFrotdp(22,1) + pi_pi*(FFrotdp(22,2) + FFrotdp(22,3)) - delta_delta*(FFrotdp(22,4) + FFrotdp(22,5)) + phi_phi*(FFrotdp(22,6) + FFrotdp(22,7));
        overlapdRdp(irw + 4,icl + 0) = -sigma_sigma*FFrotdp(5,1) +  pi_pi*(FFrotdp(5,2) + FFrotdp(5,3)) - delta_delta*(FFrotdp(5,4) + FFrotdp(5,5)) + phi_phi*(FFrotdp(5,6) + FFrotdp(5,7));
        overlapdRdp(irw + 4,icl + 1) = -sigma_sigma*FFrotdp(11,1) + pi_pi*(FFrotdp(11,2) + FFrotdp(11,3)) - delta_delta*(FFrotdp(11,4) + FFrotdp(11,5)) + phi_phi*(FFrotdp(11,6) + FFrotdp(11,7));
        overlapdRdp(irw + 4,icl + 2) = -sigma_sigma*FFrotdp(16,1) + pi_pi*(FFrotdp(16,2) + FFrotdp(16,3)) - delta_delta*(FFrotdp(16,4) + FFrotdp(16,5)) + phi_phi*(FFrotdp(16,6) + FFrotdp(16,7));
        overlapdRdp(irw + 4,icl + 3) = -sigma_sigma*FFrotdp(20,1) + pi_pi*(FFrotdp(20,2) + FFrotdp(20,3)) - delta_delta*(FFrotdp(20,4) + FFrotdp(20,5)) + phi_phi*(FFrotdp(20,6) + FFrotdp(20,7));
        overlapdRdp(irw + 4,icl + 4) = -sigma_sigma*FFrotdp(23,1) + pi_pi*(FFrotdp(23,2) + FFrotdp(23,3)) - delta_delta*(FFrotdp(23,4) + FFrotdp(23,5)) + phi_phi*(FFrotdp(23,6) + FFrotdp(23,7));
        overlapdRdp(irw + 4,icl + 5) = -sigma_sigma*FFrotdp(24,1) + pi_pi*(FFrotdp(24,2) + FFrotdp(24,3)) - delta_delta*(FFrotdp(24,4) + FFrotdp(24,5)) + phi_phi*(FFrotdp(24,6) + FFrotdp(24,7));
        overlapdRdp(irw + 4,icl + 6) = -sigma_sigma*FFrotdp(25,1) + pi_pi*(FFrotdp(25,2) + FFrotdp(25,3)) - delta_delta*(FFrotdp(25,4) + FFrotdp(25,5)) + phi_phi*(FFrotdp(25,6) + FFrotdp(25,7));
        overlapdRdp(irw + 5,icl + 0) = -sigma_sigma*FFrotdp(6,1) +  pi_pi*(FFrotdp(6,2) + FFrotdp(6,3)) - delta_delta*(FFrotdp(6,4) + FFrotdp(6,5)) + phi_phi*(FFrotdp(6,6) + FFrotdp(6,7));
        overlapdRdp(irw + 5,icl + 1) = -sigma_sigma*FFrotdp(12,1) + pi_pi*(FFrotdp(12,2) + FFrotdp(12,3)) - delta_delta*(FFrotdp(12,4) + FFrotdp(12,5)) + phi_phi*(FFrotdp(12,6) + FFrotdp(12,7));
        overlapdRdp(irw + 5,icl + 2) = -sigma_sigma*FFrotdp(17,1) + pi_pi*(FFrotdp(17,2) + FFrotdp(17,3)) - delta_delta*(FFrotdp(17,4) + FFrotdp(17,5)) + phi_phi*(FFrotdp(17,6) + FFrotdp(17,7));
        overlapdRdp(irw + 5,icl + 3) = -sigma_sigma*FFrotdp(21,1) + pi_pi*(FFrotdp(21,2) + FFrotdp(21,3)) - delta_delta*(FFrotdp(21,4) + FFrotdp(21,5)) + phi_phi*(FFrotdp(21,6) + FFrotdp(21,7));
        overlapdRdp(irw + 5,icl + 4) = -sigma_sigma*FFrotdp(24,1) + pi_pi*(FFrotdp(24,2) + FFrotdp(24,3)) - delta_delta*(FFrotdp(24,4) + FFrotdp(24,5)) + phi_phi*(FFrotdp(24,6) + FFrotdp(24,7));
        overlapdRdp(irw + 5,icl + 5) = -sigma_sigma*FFrotdp(26,1) + pi_pi*(FFrotdp(26,2) + FFrotdp(26,3)) - delta_delta*(FFrotdp(26,4) + FFrotdp(26,5)) + phi_phi*(FFrotdp(26,6) + FFrotdp(26,7));
        overlapdRdp(irw + 5,icl + 6) = -sigma_sigma*FFrotdp(27,1) + pi_pi*(FFrotdp(27,2) + FFrotdp(27,3)) - delta_delta*(FFrotdp(27,4) + FFrotdp(27,5)) + phi_phi*(FFrotdp(27,6) + FFrotdp(27,7));
        overlapdRdp(irw + 6,icl + 0) = -sigma_sigma*FFrotdp(7,1) +  pi_pi*(FFrotdp(7,2) + FFrotdp(7,3)) - delta_delta*(FFrotdp(7,4) + FFrotdp(7,5)) + phi_phi*(FFrotdp(7,6) + FFrotdp(7,7));
        overlapdRdp(irw + 6,icl + 1) = -sigma_sigma*FFrotdp(13,1) + pi_pi*(FFrotdp(13,2) + FFrotdp(13,3)) - delta_delta*(FFrotdp(13,4) + FFrotdp(13,5)) + phi_phi*(FFrotdp(13,6) + FFrotdp(13,7));
        overlapdRdp(irw + 6,icl + 2) = -sigma_sigma*FFrotdp(18,1) + pi_pi*(FFrotdp(18,2) + FFrotdp(18,3)) - delta_delta*(FFrotdp(18,4) + FFrotdp(18,5)) + phi_phi*(FFrotdp(18,6) + FFrotdp(18,7));
        overlapdRdp(irw + 6,icl + 3) = -sigma_sigma*FFrotdp(22,1) + pi_pi*(FFrotdp(22,2) + FFrotdp(22,3)) - delta_delta*(FFrotdp(22,4) + FFrotdp(22,5)) + phi_phi*(FFrotdp(22,6) + FFrotdp(22,7));
        overlapdRdp(irw + 6,icl + 4) = -sigma_sigma*FFrotdp(25,1) + pi_pi*(FFrotdp(25,2) + FFrotdp(25,3)) - delta_delta*(FFrotdp(25,4) + FFrotdp(25,5)) + phi_phi*(FFrotdp(25,6) + FFrotdp(25,7));
        overlapdRdp(irw + 6,icl + 5) = -sigma_sigma*FFrotdp(27,1) + pi_pi*(FFrotdp(27,2) + FFrotdp(27,3)) - delta_delta*(FFrotdp(27,4) + FFrotdp(27,5)) + phi_phi*(FFrotdp(27,6) + FFrotdp(27,7));
        overlapdRdp(irw + 6,icl + 6) = -sigma_sigma*FFrotdp(28,1) + pi_pi*(FFrotdp(28,2) + FFrotdp(28,3)) - delta_delta*(FFrotdp(28,4) + FFrotdp(28,5)) + phi_phi*(FFrotdp(28,6) + FFrotdp(28,7));
      }
      else if ((iLA == 0)&&(iLB == 1)) {
        //sp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SPOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3 dimensions (px,py,pz)
        overlapdRdp(irw,icl + 0) = -SProtdp(1,1)*sigma_sigma;
        overlapdRdp(irw,icl + 1) = -SProtdp(2,1)*sigma_sigma;
        overlapdRdp(irw,icl + 2) = -SProtdp(3,1)*sigma_sigma;
      }
      else if ((iLA == 1)&&(iLB == 0)) {
        //ps
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3 dimensions (px,py,pz)
        overlapdRdp(irw + 0,icl) = SProtdp(1,1)*sigma_sigma;
        overlapdRdp(irw + 1,icl) = SProtdp(2,1)*sigma_sigma;
        overlapdRdp(irw + 2,icl) = SProtdp(3,1)*sigma_sigma;
      }
      else if ((iLA == 0)&&(iLB == 2)) {
        //sd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SDOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          overlapdRdp(irw,icl + index) = SDrotdp(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 2)&&(iLB == 0)) {
        //ds
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          overlapdRdp(irw + index,icl) = SDrotdp(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 0)&&(iLB == 3)) {
        //sf
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SFOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 7 dimensions (fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 7; ++index) {
          overlapdRdp(irw,icl + index) = -SFrotdp(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 3)&&(iLB == 0)) {
        //fs
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 7; ++index) {
          overlapdRdp(irw + index,icl) = SFrotdp(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 1)&&(iLB == 2)) {
        //pd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = PDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3x5 dimensions (px,py,pz)(xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 3; ++index2) {
            overlapdRdp(irw + index2,icl + index) = PDrotdp(1+index2+3*index,1)*sigma_sigma - (PDrotdp(1+index2+3*index,2) + PDrotdp(1+index2+3*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 2)&&(iLB == 1)) {
        //dp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = DPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3x5 dimensions (px,py,pz)(xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 3; ++index2) {
            overlapdRdp(irw + index,icl + index2) = -PDrotdp(1+index2+3*index,1)*sigma_sigma + (PDrotdp(1+index2+3*index,2) + PDrotdp(1+index2+3*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 1)&&(iLB == 3)) {
        //pf
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = PFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3x5 dimensions (px,py,pz)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 3; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlapdRdp(irw + index,icl + index2) = -PFrotdp(1+index2+7*index,1)*sigma_sigma + (PFrotdp(1+index2+7*index,2) + PFrotdp(1+index2+7*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 3)&&(iLB == 1)) {
        //fp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = FPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3x5 dimensions (px,py,pz)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 3; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlapdRdp(irw + index2,icl + index) = -PFrotdp(1+index2+7*index,1)*sigma_sigma + (PFrotdp(1+index2+7*index,2) + PFrotdp(1+index2+7*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 2)&&(iLB == 3)) {
        //df
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = DFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist)/dist;
        delta_delta = 0.0;
        polynomial = DFdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3x5 dimensions (xz,yz,z2,xy,x2-y2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlapdRdp(irw + index,icl + index2) = -DFrotdp(1+index2+7*index,1)*sigma_sigma + (DFrotdp(1+index2+7*index,2) + DFrotdp(1+index2+7*index,3))*pi_pi - (DFrotdp(1+index2+7*index,4) + DFrotdp(1+index2+7*index,5))*delta_delta;
          }
        }
      }
      else if ((iLA == 3)&&(iLB == 2)) {
        //fd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          sigma_sigma -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist)/dist;
        pi_pi = 0.0;
        polynomial = FDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          pi_pi -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist)/dist;
        delta_delta = 0.0;
        polynomial = FDdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*double(iNA + iNB + 3)*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argA*Ak[polynomial[1][idx] + 1]*Bk[polynomial[2][idx]];
          delta_delta -= aux*double(polynomial[0][idx])*argB*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx] + 1];
        }
        //partition over 3x5 dimensions (xz,yz,z2,xy,x2-y2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlapdRdp(irw + index2,icl + index) = DFrotdp(1+index2+7*index,1)*sigma_sigma - (DFrotdp(1+index2+7*index,2) + DFrotdp(1+index2+7*index,3))*pi_pi + (DFrotdp(1+index2+7*index,4) + DFrotdp(1+index2+7*index,5))*delta_delta;
          }
        }
      }
      icl += 2*QN2(isket + 1,2) + 1;
    }
    irw += 2*QN1(isbra + 1,2) + 1;
  }
}
void TwoCenterOv_dt2(matrixE & overlaps, size_t nao1, size_t nao2, double dist, std::vector<double> & exponents1, matrix<size_t> & QN1, std::vector<double> & exponents2, matrix<size_t> & QN2, std::vector<double> & orientation, bool rotate = true, double tolerance = 1e-6) {
  //implementation of second-derivatives of overlap matrix function with respect to theta; only up to shell f
  //calculation of the derivative of the overlap function between atoms A and B (A -> 1; B -> 2) with respect to angle theta
  overlaps.resize(nao1,nao2);
  overlaps.zero();
  //auxiliary variables for integrals
  double sigma_sigma = 0.0;
  double pi_pi = 0.0;
  double delta_delta = 0.0;
  double phi_phi = 0.0;
  //angles needed to reorient the integrals according to the position of atoms in the molecule
  double cost = orientation[2];
  double sint = sqrt(1 - cost*cost);
  double cosp = 1.0;
  double sinp = 0.0;
  if (fabs(sint) > tolerance) {
    cosp = orientation[0]/sint;
    sinp = orientation[1]/sint;
  }
  //matrices with derivatives of overlap rotations with respect to theta
  matrixE SProt_dt2(1,1);
  matrixE SDrot_dt2(1,1);
  matrixE SFrot_dt2(1,1);
  matrixE PProt_dt2(1,1);
  matrixE PDrot_dt2(1,1);
  matrixE PFrot_dt2(1,1);
  matrixE DDrot_dt2(1,1);
  matrixE DFrot_dt2(1,1);
  matrixE FFrot_dt2(1,1);
  //other auxiliaries
  std::vector<std::vector<int> > polynomial;
  double a;
  double b;
  double c;
  double argA;
  double argB;
  std::vector<double> Ak;
  std::vector<double> Bk;
  double aux;
  //determine maximum azimuthal quantum number
  size_t maxL1 = 0;
  size_t maxL2 = 0;
  for (size_t idx = 0; idx < QN1.rows(); ++idx) {
    if (QN1(idx + 1,2) > maxL1) {maxL1 = QN1(idx + 1,2);}
  }
  for (size_t idx = 0; idx < QN2.rows(); ++idx) {
    if (QN2(idx + 1,2) > maxL2) {maxL2 = QN2(idx + 1,2);}
  }
  if ((maxL1 > 0)||(maxL2 > 0)) {
    //in this case we have p orbitals
    SProt_dt2 = SPtransf_dt2(cost,sint,cosp,sinp,rotate);
    PProt_dt2 = PPtransf_dt2(cost,sint,cosp,sinp,rotate);
  }
  if ((maxL1 > 1)||(maxL2 > 1)) {
    //in this case we have d orbitals
    SDrot_dt2 = SDtransf_dt2(cost,sint,cosp,sinp,rotate);
    PDrot_dt2 = PDtransf_dt2(cost,sint,cosp,sinp,rotate);
    DDrot_dt2 = DDtransf_dt2(cost,sint,cosp,sinp,rotate);
  }
  if ((maxL1 > 2)||(maxL2 > 2)) {
    //in this case we have f orbitals
    SFrot_dt2 = SFtransf_dt2(cost,sint,cosp,sinp,rotate);
    PFrot_dt2 = PFtransf_dt2(cost,sint,cosp,sinp,rotate);
    DFrot_dt2 = DFtransf_dt2(cost,sint,cosp,sinp,rotate);
    FFrot_dt2 = FFtransf_dt2(cost,sint,cosp,sinp,rotate);
  }
  size_t iNA;
  size_t iNB;
  size_t iLA;
  size_t iLB;
  size_t irw = 1;
  size_t icl = 1;
  for (size_t isbra = 0; isbra < QN1.rows(); ++isbra) {           //loop over bra shells
    a = exponents1[isbra];
    if (a == 0.0) continue;
    iNA = QN1(isbra + 1,1) - 1;
    iLA = QN1(isbra + 1,2);
    icl = 1;
    for (size_t isket = 0; isket < QN2.rows(); ++isket) {         //loop over ket shells
      b = exponents2[isket];
      if (b == 0.0) continue;
      iNB = QN2(isket + 1,1) - 1;
      iLB = QN2(isket + 1,2);
      //getting integrals
      argA = 0.5*dist*(a + b);
      argB = 0.5*dist*(a - b);
      Ak = int_1_infty_xpk_expmyx_dx(argA,iNA + iNB + 2);
      Bk = int_m1_1_xpk_expmyx_dx_Pople(argB,iNA + iNB + 2);
      aux = OvFactor(iNA + 1,iLA,0,a,iNB + 1,iLB,0,b,dist);
      //different azimuthal quantum number cases
      //ss = 0
      if ((iLA == 1)&&(iLB == 1)) {
        //pp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = PPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 9 dimensions
        overlaps(irw + 0,icl + 0) = -sigma_sigma*PProt_dt2(1,1) + pi_pi*(PProt_dt2(1,2) + PProt_dt2(1,3));      //<x|x>
        overlaps(irw + 0,icl + 1) = -sigma_sigma*PProt_dt2(2,1) + pi_pi*(PProt_dt2(2,2) + PProt_dt2(2,3));      //<x|y>
        overlaps(irw + 0,icl + 2) = -sigma_sigma*PProt_dt2(3,1) + pi_pi*(PProt_dt2(3,2) + PProt_dt2(3,3));      //<x|z>
        overlaps(irw + 1,icl + 0) = -sigma_sigma*PProt_dt2(2,1) + pi_pi*(PProt_dt2(2,2) + PProt_dt2(2,3));      //<y|x>
        overlaps(irw + 1,icl + 1) = -sigma_sigma*PProt_dt2(4,1) + pi_pi*(PProt_dt2(4,2) + PProt_dt2(4,3));      //<y|y>
        overlaps(irw + 1,icl + 2) = -sigma_sigma*PProt_dt2(5,1) + pi_pi*(PProt_dt2(5,2) + PProt_dt2(5,3));      //<y|z>
        overlaps(irw + 2,icl + 0) = -sigma_sigma*PProt_dt2(3,1) + pi_pi*(PProt_dt2(3,2) + PProt_dt2(3,3));      //<z|x>
        overlaps(irw + 2,icl + 1) = -sigma_sigma*PProt_dt2(5,1) + pi_pi*(PProt_dt2(5,2) + PProt_dt2(5,3));      //<z|y>
        overlaps(irw + 2,icl + 2) = -sigma_sigma*PProt_dt2(6,1) + pi_pi*(PProt_dt2(6,2) + PProt_dt2(6,3));      //<z|z>
      }
      else if ((iLA == 2)&&(iLB == 2)) {
        //dd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = DDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = DDdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5x5 dimensions (xz,yz,z2,xy,x2-y2)(xz,yz,z2,xy,x2-y2)
        overlaps(irw + 0,icl + 0) = sigma_sigma*DDrot_dt2(1,1)  - pi_pi*(DDrot_dt2(1,2)  + DDrot_dt2(1,3))  + delta_delta*(DDrot_dt2(1,4)  + DDrot_dt2(1,5));             //<xz|xz>
        overlaps(irw + 0,icl + 1) = sigma_sigma*DDrot_dt2(2,1)  - pi_pi*(DDrot_dt2(2,2)  + DDrot_dt2(2,3))  + delta_delta*(DDrot_dt2(2,4)  + DDrot_dt2(2,5));             //<xz|yz>
        overlaps(irw + 0,icl + 2) = sigma_sigma*DDrot_dt2(3,1)  - pi_pi*(DDrot_dt2(3,2)  + DDrot_dt2(3,3))  + delta_delta*(DDrot_dt2(3,4)  + DDrot_dt2(3,5));             //<xz|z2>
        overlaps(irw + 0,icl + 3) = sigma_sigma*DDrot_dt2(4,1)  - pi_pi*(DDrot_dt2(4,2)  + DDrot_dt2(4,3))  + delta_delta*(DDrot_dt2(4,4)  + DDrot_dt2(4,5));             //<xz|xy>
        overlaps(irw + 0,icl + 4) = sigma_sigma*DDrot_dt2(5,1)  - pi_pi*(DDrot_dt2(5,2)  + DDrot_dt2(5,3))  + delta_delta*(DDrot_dt2(5,4)  + DDrot_dt2(5,5));             //<xz|x2-y2>
        overlaps(irw + 1,icl + 0) = sigma_sigma*DDrot_dt2(2,1)  - pi_pi*(DDrot_dt2(2,2)  + DDrot_dt2(2,3))  + delta_delta*(DDrot_dt2(2,4)  + DDrot_dt2(2,5));             //<yz|xz>
        overlaps(irw + 1,icl + 1) = sigma_sigma*DDrot_dt2(6,1)  - pi_pi*(DDrot_dt2(6,2)  + DDrot_dt2(6,3))  + delta_delta*(DDrot_dt2(6,4)  + DDrot_dt2(6,5));             //<yz|yz>
        overlaps(irw + 1,icl + 2) = sigma_sigma*DDrot_dt2(7,1)  - pi_pi*(DDrot_dt2(7,2)  + DDrot_dt2(7,3))  + delta_delta*(DDrot_dt2(7,4)  + DDrot_dt2(7,5));             //<yz|z2>
        overlaps(irw + 1,icl + 3) = sigma_sigma*DDrot_dt2(8,1)  - pi_pi*(DDrot_dt2(8,2)  + DDrot_dt2(8,3))  + delta_delta*(DDrot_dt2(8,4)  + DDrot_dt2(8,5));             //<yz|xy>
        overlaps(irw + 1,icl + 4) = sigma_sigma*DDrot_dt2(9,1)  - pi_pi*(DDrot_dt2(9,2)  + DDrot_dt2(9,3))  + delta_delta*(DDrot_dt2(9,4)  + DDrot_dt2(9,5));             //<yz|x2-y2>
        overlaps(irw + 2,icl + 0) = sigma_sigma*DDrot_dt2(3,1)  - pi_pi*(DDrot_dt2(3,2)  + DDrot_dt2(3,3))  + delta_delta*(DDrot_dt2(3,4)  + DDrot_dt2(3,5));             //<z2|xz>
        overlaps(irw + 2,icl + 1) = sigma_sigma*DDrot_dt2(7,1)  - pi_pi*(DDrot_dt2(7,2)  + DDrot_dt2(7,3))  + delta_delta*(DDrot_dt2(7,4)  + DDrot_dt2(7,5));             //<z2|yz>
        overlaps(irw + 2,icl + 2) = sigma_sigma*DDrot_dt2(10,1) - pi_pi*(DDrot_dt2(10,2) + DDrot_dt2(10,3)) + delta_delta*(DDrot_dt2(10,4) + DDrot_dt2(10,5));            //<z2|z2>
        overlaps(irw + 2,icl + 3) = sigma_sigma*DDrot_dt2(11,1) - pi_pi*(DDrot_dt2(11,2) + DDrot_dt2(11,3)) + delta_delta*(DDrot_dt2(11,4) + DDrot_dt2(11,5));            //<z2|xy>
        overlaps(irw + 2,icl + 4) = sigma_sigma*DDrot_dt2(12,1) - pi_pi*(DDrot_dt2(12,2) + DDrot_dt2(12,3)) + delta_delta*(DDrot_dt2(12,4) + DDrot_dt2(12,5));            //<z2|x2-y2>
        overlaps(irw + 3,icl + 0) = sigma_sigma*DDrot_dt2(4,1)  - pi_pi*(DDrot_dt2(4,2)  + DDrot_dt2(4,3))  + delta_delta*(DDrot_dt2(4,4)  + DDrot_dt2(4,5));             //<xy|xz>
        overlaps(irw + 3,icl + 1) = sigma_sigma*DDrot_dt2(8,1)  - pi_pi*(DDrot_dt2(8,2)  + DDrot_dt2(8,3))  + delta_delta*(DDrot_dt2(8,4)  + DDrot_dt2(8,5));             //<xy|yz>
        overlaps(irw + 3,icl + 2) = sigma_sigma*DDrot_dt2(11,1) - pi_pi*(DDrot_dt2(11,2) + DDrot_dt2(11,3)) + delta_delta*(DDrot_dt2(11,4) + DDrot_dt2(11,5));            //<xy|z2>
        overlaps(irw + 3,icl + 3) = sigma_sigma*DDrot_dt2(13,1) - pi_pi*(DDrot_dt2(13,2) + DDrot_dt2(13,3)) + delta_delta*(DDrot_dt2(13,4) + DDrot_dt2(13,5));            //<xy|xy>
        overlaps(irw + 3,icl + 4) = sigma_sigma*DDrot_dt2(14,1) - pi_pi*(DDrot_dt2(14,2) + DDrot_dt2(14,3)) + delta_delta*(DDrot_dt2(14,4) + DDrot_dt2(14,5));            //<xy|x2-y2>
        overlaps(irw + 4,icl + 0) = sigma_sigma*DDrot_dt2(5,1)  - pi_pi*(DDrot_dt2(5,2)  + DDrot_dt2(5,3))  + delta_delta*(DDrot_dt2(5,4)  + DDrot_dt2(5,5));             //<x2-y2|xz>
        overlaps(irw + 4,icl + 1) = sigma_sigma*DDrot_dt2(9,1)  - pi_pi*(DDrot_dt2(9,2)  + DDrot_dt2(9,3))  + delta_delta*(DDrot_dt2(9,4)  + DDrot_dt2(9,5));             //<x2-y2|yz>
        overlaps(irw + 4,icl + 2) = sigma_sigma*DDrot_dt2(12,1) - pi_pi*(DDrot_dt2(12,2) + DDrot_dt2(12,3)) + delta_delta*(DDrot_dt2(12,4) + DDrot_dt2(12,5));            //<x2-y2|z2>
        overlaps(irw + 4,icl + 3) = sigma_sigma*DDrot_dt2(14,1) - pi_pi*(DDrot_dt2(14,2) + DDrot_dt2(14,3)) + delta_delta*(DDrot_dt2(14,4) + DDrot_dt2(14,5));            //<x2-y2|xy>
        overlaps(irw + 4,icl + 4) = sigma_sigma*DDrot_dt2(15,1) - pi_pi*(DDrot_dt2(15,2) + DDrot_dt2(15,3)) + delta_delta*(DDrot_dt2(15,4) + DDrot_dt2(15,5));            //<x2-y2|x2-y2>
      }
      else if ((iLA == 3)&&(iLB == 3)) {
        //ff
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = FFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = FFdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //phi integral
        aux = OvFactor(iNA + 1,iLA,3,a,iNB + 1,iLB,3,b,dist);
        phi_phi = 0.0;
        polynomial = FFphiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          phi_phi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 7x7 dimensions (fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        overlaps(irw + 0,icl + 0) = -sigma_sigma*FFrot_dt2(1,1)  + pi_pi*(FFrot_dt2(1,2)  + FFrot_dt2(1,3))  - delta_delta*(FFrot_dt2(1,4)  + FFrot_dt2(1,5))  + phi_phi*(FFrot_dt2(1,6)  + FFrot_dt2(1,7));
        overlaps(irw + 0,icl + 1) = -sigma_sigma*FFrot_dt2(2,1)  + pi_pi*(FFrot_dt2(2,2)  + FFrot_dt2(2,3))  - delta_delta*(FFrot_dt2(2,4)  + FFrot_dt2(2,5))  + phi_phi*(FFrot_dt2(2,6)  + FFrot_dt2(2,7));
        overlaps(irw + 0,icl + 2) = -sigma_sigma*FFrot_dt2(3,1)  + pi_pi*(FFrot_dt2(3,2)  + FFrot_dt2(3,3))  - delta_delta*(FFrot_dt2(3,4)  + FFrot_dt2(3,5))  + phi_phi*(FFrot_dt2(3,6)  + FFrot_dt2(3,7));
        overlaps(irw + 0,icl + 3) = -sigma_sigma*FFrot_dt2(4,1)  + pi_pi*(FFrot_dt2(4,2)  + FFrot_dt2(4,3))  - delta_delta*(FFrot_dt2(4,4)  + FFrot_dt2(4,5))  + phi_phi*(FFrot_dt2(4,6)  + FFrot_dt2(4,7));
        overlaps(irw + 0,icl + 4) = -sigma_sigma*FFrot_dt2(5,1)  + pi_pi*(FFrot_dt2(5,2)  + FFrot_dt2(5,3))  - delta_delta*(FFrot_dt2(5,4)  + FFrot_dt2(5,5))  + phi_phi*(FFrot_dt2(5,6)  + FFrot_dt2(5,7));
        overlaps(irw + 0,icl + 5) = -sigma_sigma*FFrot_dt2(6,1)  + pi_pi*(FFrot_dt2(6,2)  + FFrot_dt2(6,3))  - delta_delta*(FFrot_dt2(6,4)  + FFrot_dt2(6,5))  + phi_phi*(FFrot_dt2(6,6)  + FFrot_dt2(6,7));
        overlaps(irw + 0,icl + 6) = -sigma_sigma*FFrot_dt2(7,1)  + pi_pi*(FFrot_dt2(7,2)  + FFrot_dt2(7,3))  - delta_delta*(FFrot_dt2(7,4)  + FFrot_dt2(7,5))  + phi_phi*(FFrot_dt2(7,6)  + FFrot_dt2(7,7));
        overlaps(irw + 1,icl + 0) = -sigma_sigma*FFrot_dt2(2,1)  + pi_pi*(FFrot_dt2(2,2)  + FFrot_dt2(2,3))  - delta_delta*(FFrot_dt2(2,4)  + FFrot_dt2(2,5))  + phi_phi*(FFrot_dt2(2,6)  + FFrot_dt2(2,7));
        overlaps(irw + 1,icl + 1) = -sigma_sigma*FFrot_dt2(8,1)  + pi_pi*(FFrot_dt2(8,2)  + FFrot_dt2(8,3))  - delta_delta*(FFrot_dt2(8,4)  + FFrot_dt2(8,5))  + phi_phi*(FFrot_dt2(8,6)  + FFrot_dt2(8,7));
        overlaps(irw + 1,icl + 2) = -sigma_sigma*FFrot_dt2(9,1)  + pi_pi*(FFrot_dt2(9,2)  + FFrot_dt2(9,3))  - delta_delta*(FFrot_dt2(9,4)  + FFrot_dt2(9,5))  + phi_phi*(FFrot_dt2(9,6)  + FFrot_dt2(9,7));
        overlaps(irw + 1,icl + 3) = -sigma_sigma*FFrot_dt2(10,1) + pi_pi*(FFrot_dt2(10,2) + FFrot_dt2(10,3)) - delta_delta*(FFrot_dt2(10,4) + FFrot_dt2(10,5)) + phi_phi*(FFrot_dt2(10,6) + FFrot_dt2(10,7));
        overlaps(irw + 1,icl + 4) = -sigma_sigma*FFrot_dt2(11,1) + pi_pi*(FFrot_dt2(11,2) + FFrot_dt2(11,3)) - delta_delta*(FFrot_dt2(11,4) + FFrot_dt2(11,5)) + phi_phi*(FFrot_dt2(11,6) + FFrot_dt2(11,7));
        overlaps(irw + 1,icl + 5) = -sigma_sigma*FFrot_dt2(12,1) + pi_pi*(FFrot_dt2(12,2) + FFrot_dt2(12,3)) - delta_delta*(FFrot_dt2(12,4) + FFrot_dt2(12,5)) + phi_phi*(FFrot_dt2(12,6) + FFrot_dt2(12,7));
        overlaps(irw + 1,icl + 6) = -sigma_sigma*FFrot_dt2(13,1) + pi_pi*(FFrot_dt2(13,2) + FFrot_dt2(13,3)) - delta_delta*(FFrot_dt2(13,4) + FFrot_dt2(13,5)) + phi_phi*(FFrot_dt2(13,6) + FFrot_dt2(13,7));
        overlaps(irw + 2,icl + 0) = -sigma_sigma*FFrot_dt2(3,1) +  pi_pi*(FFrot_dt2(3,2)  + FFrot_dt2(3,3))  - delta_delta*(FFrot_dt2(3,4)  + FFrot_dt2(3,5))  + phi_phi*(FFrot_dt2(3,6)  + FFrot_dt2(3,7));
        overlaps(irw + 2,icl + 1) = -sigma_sigma*FFrot_dt2(9,1) +  pi_pi*(FFrot_dt2(9,2)  + FFrot_dt2(9,3))  - delta_delta*(FFrot_dt2(9,4)  + FFrot_dt2(9,5))  + phi_phi*(FFrot_dt2(9,6)  + FFrot_dt2(9,7));
        overlaps(irw + 2,icl + 2) = -sigma_sigma*FFrot_dt2(14,1) + pi_pi*(FFrot_dt2(14,2) + FFrot_dt2(14,3)) - delta_delta*(FFrot_dt2(14,4) + FFrot_dt2(14,5)) + phi_phi*(FFrot_dt2(14,6) + FFrot_dt2(14,7));
        overlaps(irw + 2,icl + 3) = -sigma_sigma*FFrot_dt2(15,1) + pi_pi*(FFrot_dt2(15,2) + FFrot_dt2(15,3)) - delta_delta*(FFrot_dt2(15,4) + FFrot_dt2(15,5)) + phi_phi*(FFrot_dt2(15,6) + FFrot_dt2(15,7));
        overlaps(irw + 2,icl + 4) = -sigma_sigma*FFrot_dt2(16,1) + pi_pi*(FFrot_dt2(16,2) + FFrot_dt2(16,3)) - delta_delta*(FFrot_dt2(16,4) + FFrot_dt2(16,5)) + phi_phi*(FFrot_dt2(16,6) + FFrot_dt2(16,7));
        overlaps(irw + 2,icl + 5) = -sigma_sigma*FFrot_dt2(17,1) + pi_pi*(FFrot_dt2(17,2) + FFrot_dt2(17,3)) - delta_delta*(FFrot_dt2(17,4) + FFrot_dt2(17,5)) + phi_phi*(FFrot_dt2(17,6) + FFrot_dt2(17,7));
        overlaps(irw + 2,icl + 6) = -sigma_sigma*FFrot_dt2(18,1) + pi_pi*(FFrot_dt2(18,2) + FFrot_dt2(18,3)) - delta_delta*(FFrot_dt2(18,4) + FFrot_dt2(18,5)) + phi_phi*(FFrot_dt2(18,6) + FFrot_dt2(18,7));
        overlaps(irw + 3,icl + 0) = -sigma_sigma*FFrot_dt2(4,1) +  pi_pi*(FFrot_dt2(4,2)  + FFrot_dt2(4,3))  - delta_delta*(FFrot_dt2(4,4)  + FFrot_dt2(4,5))  + phi_phi*(FFrot_dt2(4,6)  + FFrot_dt2(4,7));
        overlaps(irw + 3,icl + 1) = -sigma_sigma*FFrot_dt2(10,1) + pi_pi*(FFrot_dt2(10,2) + FFrot_dt2(10,3)) - delta_delta*(FFrot_dt2(10,4) + FFrot_dt2(10,5)) + phi_phi*(FFrot_dt2(10,6) + FFrot_dt2(10,7));
        overlaps(irw + 3,icl + 2) = -sigma_sigma*FFrot_dt2(15,1) + pi_pi*(FFrot_dt2(15,2) + FFrot_dt2(15,3)) - delta_delta*(FFrot_dt2(15,4) + FFrot_dt2(15,5)) + phi_phi*(FFrot_dt2(15,6) + FFrot_dt2(15,7));
        overlaps(irw + 3,icl + 3) = -sigma_sigma*FFrot_dt2(19,1) + pi_pi*(FFrot_dt2(19,2) + FFrot_dt2(19,3)) - delta_delta*(FFrot_dt2(19,4) + FFrot_dt2(19,5)) + phi_phi*(FFrot_dt2(19,6) + FFrot_dt2(19,7));
        overlaps(irw + 3,icl + 4) = -sigma_sigma*FFrot_dt2(20,1) + pi_pi*(FFrot_dt2(20,2) + FFrot_dt2(20,3)) - delta_delta*(FFrot_dt2(20,4) + FFrot_dt2(20,5)) + phi_phi*(FFrot_dt2(20,6) + FFrot_dt2(20,7));
        overlaps(irw + 3,icl + 5) = -sigma_sigma*FFrot_dt2(21,1) + pi_pi*(FFrot_dt2(21,2) + FFrot_dt2(21,3)) - delta_delta*(FFrot_dt2(21,4) + FFrot_dt2(21,5)) + phi_phi*(FFrot_dt2(21,6) + FFrot_dt2(21,7));
        overlaps(irw + 3,icl + 6) = -sigma_sigma*FFrot_dt2(22,1) + pi_pi*(FFrot_dt2(22,2) + FFrot_dt2(22,3)) - delta_delta*(FFrot_dt2(22,4) + FFrot_dt2(22,5)) + phi_phi*(FFrot_dt2(22,6) + FFrot_dt2(22,7));
        overlaps(irw + 4,icl + 0) = -sigma_sigma*FFrot_dt2(5,1) +  pi_pi*(FFrot_dt2(5,2)  + FFrot_dt2(5,3))  - delta_delta*(FFrot_dt2(5,4)  + FFrot_dt2(5,5))  + phi_phi*(FFrot_dt2(5,6)  + FFrot_dt2(5,7));
        overlaps(irw + 4,icl + 1) = -sigma_sigma*FFrot_dt2(11,1) + pi_pi*(FFrot_dt2(11,2) + FFrot_dt2(11,3)) - delta_delta*(FFrot_dt2(11,4) + FFrot_dt2(11,5)) + phi_phi*(FFrot_dt2(11,6) + FFrot_dt2(11,7));
        overlaps(irw + 4,icl + 2) = -sigma_sigma*FFrot_dt2(16,1) + pi_pi*(FFrot_dt2(16,2) + FFrot_dt2(16,3)) - delta_delta*(FFrot_dt2(16,4) + FFrot_dt2(16,5)) + phi_phi*(FFrot_dt2(16,6) + FFrot_dt2(16,7));
        overlaps(irw + 4,icl + 3) = -sigma_sigma*FFrot_dt2(20,1) + pi_pi*(FFrot_dt2(20,2) + FFrot_dt2(20,3)) - delta_delta*(FFrot_dt2(20,4) + FFrot_dt2(20,5)) + phi_phi*(FFrot_dt2(20,6) + FFrot_dt2(20,7));
        overlaps(irw + 4,icl + 4) = -sigma_sigma*FFrot_dt2(23,1) + pi_pi*(FFrot_dt2(23,2) + FFrot_dt2(23,3)) - delta_delta*(FFrot_dt2(23,4) + FFrot_dt2(23,5)) + phi_phi*(FFrot_dt2(23,6) + FFrot_dt2(23,7));
        overlaps(irw + 4,icl + 5) = -sigma_sigma*FFrot_dt2(24,1) + pi_pi*(FFrot_dt2(24,2) + FFrot_dt2(24,3)) - delta_delta*(FFrot_dt2(24,4) + FFrot_dt2(24,5)) + phi_phi*(FFrot_dt2(24,6) + FFrot_dt2(24,7));
        overlaps(irw + 4,icl + 6) = -sigma_sigma*FFrot_dt2(25,1) + pi_pi*(FFrot_dt2(25,2) + FFrot_dt2(25,3)) - delta_delta*(FFrot_dt2(25,4) + FFrot_dt2(25,5)) + phi_phi*(FFrot_dt2(25,6) + FFrot_dt2(25,7));
        overlaps(irw + 5,icl + 0) = -sigma_sigma*FFrot_dt2(6,1) +  pi_pi*(FFrot_dt2(6,2)  + FFrot_dt2(6,3))  - delta_delta*(FFrot_dt2(6,4)  + FFrot_dt2(6,5))  + phi_phi*(FFrot_dt2(6,6)  + FFrot_dt2(6,7));
        overlaps(irw + 5,icl + 1) = -sigma_sigma*FFrot_dt2(12,1) + pi_pi*(FFrot_dt2(12,2) + FFrot_dt2(12,3)) - delta_delta*(FFrot_dt2(12,4) + FFrot_dt2(12,5)) + phi_phi*(FFrot_dt2(12,6) + FFrot_dt2(12,7));
        overlaps(irw + 5,icl + 2) = -sigma_sigma*FFrot_dt2(17,1) + pi_pi*(FFrot_dt2(17,2) + FFrot_dt2(17,3)) - delta_delta*(FFrot_dt2(17,4) + FFrot_dt2(17,5)) + phi_phi*(FFrot_dt2(17,6) + FFrot_dt2(17,7));
        overlaps(irw + 5,icl + 3) = -sigma_sigma*FFrot_dt2(21,1) + pi_pi*(FFrot_dt2(21,2) + FFrot_dt2(21,3)) - delta_delta*(FFrot_dt2(21,4) + FFrot_dt2(21,5)) + phi_phi*(FFrot_dt2(21,6) + FFrot_dt2(21,7));
        overlaps(irw + 5,icl + 4) = -sigma_sigma*FFrot_dt2(24,1) + pi_pi*(FFrot_dt2(24,2) + FFrot_dt2(24,3)) - delta_delta*(FFrot_dt2(24,4) + FFrot_dt2(24,5)) + phi_phi*(FFrot_dt2(24,6) + FFrot_dt2(24,7));
        overlaps(irw + 5,icl + 5) = -sigma_sigma*FFrot_dt2(26,1) + pi_pi*(FFrot_dt2(26,2) + FFrot_dt2(26,3)) - delta_delta*(FFrot_dt2(26,4) + FFrot_dt2(26,5)) + phi_phi*(FFrot_dt2(26,6) + FFrot_dt2(26,7));
        overlaps(irw + 5,icl + 6) = -sigma_sigma*FFrot_dt2(27,1) + pi_pi*(FFrot_dt2(27,2) + FFrot_dt2(27,3)) - delta_delta*(FFrot_dt2(27,4) + FFrot_dt2(27,5)) + phi_phi*(FFrot_dt2(27,6) + FFrot_dt2(27,7));
        overlaps(irw + 6,icl + 0) = -sigma_sigma*FFrot_dt2(7,1) +  pi_pi*(FFrot_dt2(7,2)  + FFrot_dt2(7,3))  - delta_delta*(FFrot_dt2(7,4)  + FFrot_dt2(7,5))  + phi_phi*(FFrot_dt2(7,6)  + FFrot_dt2(7,7));
        overlaps(irw + 6,icl + 1) = -sigma_sigma*FFrot_dt2(13,1) + pi_pi*(FFrot_dt2(13,2) + FFrot_dt2(13,3)) - delta_delta*(FFrot_dt2(13,4) + FFrot_dt2(13,5)) + phi_phi*(FFrot_dt2(13,6) + FFrot_dt2(13,7));
        overlaps(irw + 6,icl + 2) = -sigma_sigma*FFrot_dt2(18,1) + pi_pi*(FFrot_dt2(18,2) + FFrot_dt2(18,3)) - delta_delta*(FFrot_dt2(18,4) + FFrot_dt2(18,5)) + phi_phi*(FFrot_dt2(18,6) + FFrot_dt2(18,7));
        overlaps(irw + 6,icl + 3) = -sigma_sigma*FFrot_dt2(22,1) + pi_pi*(FFrot_dt2(22,2) + FFrot_dt2(22,3)) - delta_delta*(FFrot_dt2(22,4) + FFrot_dt2(22,5)) + phi_phi*(FFrot_dt2(22,6) + FFrot_dt2(22,7));
        overlaps(irw + 6,icl + 4) = -sigma_sigma*FFrot_dt2(25,1) + pi_pi*(FFrot_dt2(25,2) + FFrot_dt2(25,3)) - delta_delta*(FFrot_dt2(25,4) + FFrot_dt2(25,5)) + phi_phi*(FFrot_dt2(25,6) + FFrot_dt2(25,7));
        overlaps(irw + 6,icl + 5) = -sigma_sigma*FFrot_dt2(27,1) + pi_pi*(FFrot_dt2(27,2) + FFrot_dt2(27,3)) - delta_delta*(FFrot_dt2(27,4) + FFrot_dt2(27,5)) + phi_phi*(FFrot_dt2(27,6) + FFrot_dt2(27,7));
        overlaps(irw + 6,icl + 6) = -sigma_sigma*FFrot_dt2(28,1) + pi_pi*(FFrot_dt2(28,2) + FFrot_dt2(28,3)) - delta_delta*(FFrot_dt2(28,4) + FFrot_dt2(28,5)) + phi_phi*(FFrot_dt2(28,6) + FFrot_dt2(28,7));
      }
      else if ((iLA == 0)&&(iLB == 1)) {
        //sp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SPOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3 dimensions (px,py,pz)
        overlaps(irw,icl + 0) = -SProt_dt2(1,1)*sigma_sigma;
        overlaps(irw,icl + 1) = -SProt_dt2(2,1)*sigma_sigma;
        overlaps(irw,icl + 2) = -SProt_dt2(3,1)*sigma_sigma;
      }
      else if ((iLA == 1)&&(iLB == 0)) {
        //ps
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3 dimensions (px,py,pz)
        overlaps(irw + 0,icl) = SProt_dt2(1,1)*sigma_sigma;
        overlaps(irw + 1,icl) = SProt_dt2(2,1)*sigma_sigma;
        overlaps(irw + 2,icl) = SProt_dt2(3,1)*sigma_sigma;
      }
      else if ((iLA == 0)&&(iLB == 2)) {
        //sd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SDOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          overlaps(irw,icl + index) = SDrot_dt2(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 2)&&(iLB == 0)) {
        //ds
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          overlaps(irw + index,icl) = SDrot_dt2(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 0)&&(iLB == 3)) {
        //sf
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SFOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 7 dimensions (fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 7; ++index) {
          overlaps(irw,icl + index) = -SFrot_dt2(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 3)&&(iLB == 0)) {
        //fs
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 7; ++index) {
          overlaps(irw + index,icl) = SFrot_dt2(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 1)&&(iLB == 2)) {
        //pd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = PDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 3; ++index2) {
            overlaps(irw + index2,icl + index) = PDrot_dt2(1+index2+3*index,1)*sigma_sigma - (PDrot_dt2(1+index2+3*index,2) + PDrot_dt2(1+index2+3*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 2)&&(iLB == 1)) {
        //dp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = DPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 3; ++index2) {
            overlaps(irw + index,icl + index2) = -PDrot_dt2(1+index2+3*index,1)*sigma_sigma + (PDrot_dt2(1+index2+3*index,2) + PDrot_dt2(1+index2+3*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 1)&&(iLB == 3)) {
        //pf
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = PFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 3; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index,icl + index2) = -PFrot_dt2(1+index2+7*index,1)*sigma_sigma + (PFrot_dt2(1+index2+7*index,2) + PFrot_dt2(1+index2+7*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 3)&&(iLB == 1)) {
        //fp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = FPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 3; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index2,icl + index) = -PFrot_dt2(1+index2+7*index,1)*sigma_sigma + (PFrot_dt2(1+index2+7*index,2) + PFrot_dt2(1+index2+7*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 2)&&(iLB == 3)) {
        //df
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = DFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = DFdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (xz,yz,z2,xy,x2-y2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index,icl + index2) = -DFrot_dt2(1+index2+7*index,1)*sigma_sigma + (DFrot_dt2(1+index2+7*index,2) + DFrot_dt2(1+index2+7*index,3))*pi_pi - (DFrot_dt2(1+index2+7*index,4) + DFrot_dt2(1+index2+7*index,5))*delta_delta;
          }
        }
      }
      else if ((iLA == 3)&&(iLB == 2)) {
        //fd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = FDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = FDdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (xz,yz,z2,xy,x2-y2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index2,icl + index) = DFrot_dt2(1+index2+7*index,1)*sigma_sigma - (DFrot_dt2(1+index2+7*index,2) + DFrot_dt2(1+index2+7*index,3))*pi_pi + (DFrot_dt2(1+index2+7*index,4) + DFrot_dt2(1+index2+7*index,5))*delta_delta;
          }
        }
      }
      icl += 2*QN2(isket + 1,2) + 1;
    }
    irw += 2*QN1(isbra + 1,2) + 1;
  }
}
void TwoCenterOv_dp2(matrixE & overlaps, size_t nao1, size_t nao2, double dist, std::vector<double> & exponents1, matrix<size_t> & QN1, std::vector<double> & exponents2, matrix<size_t> & QN2, std::vector<double> & orientation, bool rotate = true, double tolerance = 1e-6) {
  //implementation of second-derivatives of overlap matrix function with respect to phi; only up to shell f
  //calculation of the derivative of the overlap function between atoms A and B (A -> 1; B -> 2) with respect to angle phi
  overlaps.resize(nao1,nao2);
  overlaps.zero();
  //auxiliary variables for integrals
  double sigma_sigma = 0.0;
  double pi_pi = 0.0;
  double delta_delta = 0.0;
  double phi_phi = 0.0;
  //angles needed to reorient the integrals according to the position of atoms in the molecule
  double cost = orientation[2];
  double sint = sqrt(1 - cost*cost);
  double cosp = 1.0;
  double sinp = 0.0;
  if (fabs(sint) > tolerance) {
    cosp = orientation[0]/sint;
    sinp = orientation[1]/sint;
  }
  //matrices with derivatives of overlap rotations with respect to phi
  matrixE SProt_dp2(1,1);
  matrixE SDrot_dp2(1,1);
  matrixE SFrot_dp2(1,1);
  matrixE PProt_dp2(1,1);
  matrixE PDrot_dp2(1,1);
  matrixE PFrot_dp2(1,1);
  matrixE DDrot_dp2(1,1);
  matrixE DFrot_dp2(1,1);
  matrixE FFrot_dp2(1,1);
  //other auxiliaries
  std::vector<std::vector<int> > polynomial;
  double a;
  double b;
  double c;
  double argA;
  double argB;
  std::vector<double> Ak;
  std::vector<double> Bk;
  double aux;
  //determine maximum azimuthal quantum number
  size_t maxL1 = 0;
  size_t maxL2 = 0;
  for (size_t idx = 0; idx < QN1.rows(); ++idx) {
    if (QN1(idx + 1,2) > maxL1) {maxL1 = QN1(idx + 1,2);}
  }
  for (size_t idx = 0; idx < QN2.rows(); ++idx) {
    if (QN2(idx + 1,2) > maxL2) {maxL2 = QN2(idx + 1,2);}
  }
  if ((maxL1 > 0)||(maxL2 > 0)) {
    //in this case we have p orbitals
    SProt_dp2 = SPtransf_dp2(cost,sint,cosp,sinp,rotate);
    PProt_dp2 = PPtransf_dp2(cost,sint,cosp,sinp,rotate);
  }
  if ((maxL1 > 1)||(maxL2 > 1)) {
    //in this case we have d orbitals
    SDrot_dp2 = SDtransf_dp2(cost,sint,cosp,sinp,rotate);
    PDrot_dp2 = PDtransf_dp2(cost,sint,cosp,sinp,rotate);
    DDrot_dp2 = DDtransf_dp2(cost,sint,cosp,sinp,rotate);
  }
  if ((maxL1 > 2)||(maxL2 > 2)) {
    //in this case we have f orbitals
    SFrot_dp2 = SFtransf_dp2(cost,sint,cosp,sinp,rotate);
    PFrot_dp2 = PFtransf_dp2(cost,sint,cosp,sinp,rotate);
    DFrot_dp2 = DFtransf_dp2(cost,sint,cosp,sinp,rotate);
    FFrot_dp2 = FFtransf_dp2(cost,sint,cosp,sinp,rotate);
  }
  size_t iNA;
  size_t iNB;
  size_t iLA;
  size_t iLB;
  size_t irw = 1;
  size_t icl = 1;
  for (size_t isbra = 0; isbra < QN1.rows(); ++isbra) {           //loop over bra shells
    a = exponents1[isbra];
    if (a == 0.0) continue;
    iNA = QN1(isbra + 1,1) - 1;
    iLA = QN1(isbra + 1,2);
    icl = 1;
    for (size_t isket = 0; isket < QN2.rows(); ++isket) {         //loop over ket shells
      b = exponents2[isket];
      if (b == 0.0) continue;
      iNB = QN2(isket + 1,1) - 1;
      iLB = QN2(isket + 1,2);
      //getting integrals
      argA = 0.5*dist*(a + b);
      argB = 0.5*dist*(a - b);
      Ak = int_1_infty_xpk_expmyx_dx(argA,iNA + iNB + 2);
      Bk = int_m1_1_xpk_expmyx_dx_Pople(argB,iNA + iNB + 2);
      aux = OvFactor(iNA + 1,iLA,0,a,iNB + 1,iLB,0,b,dist);
      //different azimuthal quantum number cases
      //ss = 0
      if ((iLA == 1)&&(iLB == 1)) {
        //pp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = PPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 9 dimensions
        overlaps(irw + 0,icl + 0) = -sigma_sigma*PProt_dp2(1,1) + pi_pi*(PProt_dp2(1,2) + PProt_dp2(1,3));      //<x|x>
        overlaps(irw + 0,icl + 1) = -sigma_sigma*PProt_dp2(2,1) + pi_pi*(PProt_dp2(2,2) + PProt_dp2(2,3));      //<x|y>
        overlaps(irw + 0,icl + 2) = -sigma_sigma*PProt_dp2(3,1) + pi_pi*(PProt_dp2(3,2) + PProt_dp2(3,3));      //<x|z>
        overlaps(irw + 1,icl + 0) = -sigma_sigma*PProt_dp2(2,1) + pi_pi*(PProt_dp2(2,2) + PProt_dp2(2,3));      //<y|x>
        overlaps(irw + 1,icl + 1) = -sigma_sigma*PProt_dp2(4,1) + pi_pi*(PProt_dp2(4,2) + PProt_dp2(4,3));      //<y|y>
        overlaps(irw + 1,icl + 2) = -sigma_sigma*PProt_dp2(5,1) + pi_pi*(PProt_dp2(5,2) + PProt_dp2(5,3));      //<y|z>
        overlaps(irw + 2,icl + 0) = -sigma_sigma*PProt_dp2(3,1) + pi_pi*(PProt_dp2(3,2) + PProt_dp2(3,3));      //<z|x>
        overlaps(irw + 2,icl + 1) = -sigma_sigma*PProt_dp2(5,1) + pi_pi*(PProt_dp2(5,2) + PProt_dp2(5,3));      //<z|y>
        overlaps(irw + 2,icl + 2) = -sigma_sigma*PProt_dp2(6,1) + pi_pi*(PProt_dp2(6,2) + PProt_dp2(6,3));      //<z|z>
      }
      else if ((iLA == 2)&&(iLB == 2)) {
        //dd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = DDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = DDdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5x5 dimensions (xz,yz,z2,xy,x2-y2)(xz,yz,z2,xy,x2-y2)
        overlaps(irw + 0,icl + 0) = sigma_sigma*DDrot_dp2(1,1)  - pi_pi*(DDrot_dp2(1,2)  + DDrot_dp2(1,3))  + delta_delta*(DDrot_dp2(1,4)  + DDrot_dp2(1,5));             //<xz|xz>
        overlaps(irw + 0,icl + 1) = sigma_sigma*DDrot_dp2(2,1)  - pi_pi*(DDrot_dp2(2,2)  + DDrot_dp2(2,3))  + delta_delta*(DDrot_dp2(2,4)  + DDrot_dp2(2,5));             //<xz|yz>
        overlaps(irw + 0,icl + 2) = sigma_sigma*DDrot_dp2(3,1)  - pi_pi*(DDrot_dp2(3,2)  + DDrot_dp2(3,3))  + delta_delta*(DDrot_dp2(3,4)  + DDrot_dp2(3,5));             //<xz|z2>
        overlaps(irw + 0,icl + 3) = sigma_sigma*DDrot_dp2(4,1)  - pi_pi*(DDrot_dp2(4,2)  + DDrot_dp2(4,3))  + delta_delta*(DDrot_dp2(4,4)  + DDrot_dp2(4,5));             //<xz|xy>
        overlaps(irw + 0,icl + 4) = sigma_sigma*DDrot_dp2(5,1)  - pi_pi*(DDrot_dp2(5,2)  + DDrot_dp2(5,3))  + delta_delta*(DDrot_dp2(5,4)  + DDrot_dp2(5,5));             //<xz|x2-y2>
        overlaps(irw + 1,icl + 0) = sigma_sigma*DDrot_dp2(2,1)  - pi_pi*(DDrot_dp2(2,2)  + DDrot_dp2(2,3))  + delta_delta*(DDrot_dp2(2,4)  + DDrot_dp2(2,5));             //<yz|xz>
        overlaps(irw + 1,icl + 1) = sigma_sigma*DDrot_dp2(6,1)  - pi_pi*(DDrot_dp2(6,2)  + DDrot_dp2(6,3))  + delta_delta*(DDrot_dp2(6,4)  + DDrot_dp2(6,5));             //<yz|yz>
        overlaps(irw + 1,icl + 2) = sigma_sigma*DDrot_dp2(7,1)  - pi_pi*(DDrot_dp2(7,2)  + DDrot_dp2(7,3))  + delta_delta*(DDrot_dp2(7,4)  + DDrot_dp2(7,5));             //<yz|z2>
        overlaps(irw + 1,icl + 3) = sigma_sigma*DDrot_dp2(8,1)  - pi_pi*(DDrot_dp2(8,2)  + DDrot_dp2(8,3))  + delta_delta*(DDrot_dp2(8,4)  + DDrot_dp2(8,5));             //<yz|xy>
        overlaps(irw + 1,icl + 4) = sigma_sigma*DDrot_dp2(9,1)  - pi_pi*(DDrot_dp2(9,2)  + DDrot_dp2(9,3))  + delta_delta*(DDrot_dp2(9,4)  + DDrot_dp2(9,5));             //<yz|x2-y2>
        overlaps(irw + 2,icl + 0) = sigma_sigma*DDrot_dp2(3,1)  - pi_pi*(DDrot_dp2(3,2)  + DDrot_dp2(3,3))  + delta_delta*(DDrot_dp2(3,4)  + DDrot_dp2(3,5));             //<z2|xz>
        overlaps(irw + 2,icl + 1) = sigma_sigma*DDrot_dp2(7,1)  - pi_pi*(DDrot_dp2(7,2)  + DDrot_dp2(7,3))  + delta_delta*(DDrot_dp2(7,4)  + DDrot_dp2(7,5));             //<z2|yz>
        overlaps(irw + 2,icl + 2) = sigma_sigma*DDrot_dp2(10,1) - pi_pi*(DDrot_dp2(10,2) + DDrot_dp2(10,3)) + delta_delta*(DDrot_dp2(10,4) + DDrot_dp2(10,5));            //<z2|z2>
        overlaps(irw + 2,icl + 3) = sigma_sigma*DDrot_dp2(11,1) - pi_pi*(DDrot_dp2(11,2) + DDrot_dp2(11,3)) + delta_delta*(DDrot_dp2(11,4) + DDrot_dp2(11,5));            //<z2|xy>
        overlaps(irw + 2,icl + 4) = sigma_sigma*DDrot_dp2(12,1) - pi_pi*(DDrot_dp2(12,2) + DDrot_dp2(12,3)) + delta_delta*(DDrot_dp2(12,4) + DDrot_dp2(12,5));            //<z2|x2-y2>
        overlaps(irw + 3,icl + 0) = sigma_sigma*DDrot_dp2(4,1)  - pi_pi*(DDrot_dp2(4,2)  + DDrot_dp2(4,3))  + delta_delta*(DDrot_dp2(4,4)  + DDrot_dp2(4,5));             //<xy|xz>
        overlaps(irw + 3,icl + 1) = sigma_sigma*DDrot_dp2(8,1)  - pi_pi*(DDrot_dp2(8,2)  + DDrot_dp2(8,3))  + delta_delta*(DDrot_dp2(8,4)  + DDrot_dp2(8,5));             //<xy|yz>
        overlaps(irw + 3,icl + 2) = sigma_sigma*DDrot_dp2(11,1) - pi_pi*(DDrot_dp2(11,2) + DDrot_dp2(11,3)) + delta_delta*(DDrot_dp2(11,4) + DDrot_dp2(11,5));            //<xy|z2>
        overlaps(irw + 3,icl + 3) = sigma_sigma*DDrot_dp2(13,1) - pi_pi*(DDrot_dp2(13,2) + DDrot_dp2(13,3)) + delta_delta*(DDrot_dp2(13,4) + DDrot_dp2(13,5));            //<xy|xy>
        overlaps(irw + 3,icl + 4) = sigma_sigma*DDrot_dp2(14,1) - pi_pi*(DDrot_dp2(14,2) + DDrot_dp2(14,3)) + delta_delta*(DDrot_dp2(14,4) + DDrot_dp2(14,5));            //<xy|x2-y2>
        overlaps(irw + 4,icl + 0) = sigma_sigma*DDrot_dp2(5,1)  - pi_pi*(DDrot_dp2(5,2)  + DDrot_dp2(5,3))  + delta_delta*(DDrot_dp2(5,4)  + DDrot_dp2(5,5));             //<x2-y2|xz>
        overlaps(irw + 4,icl + 1) = sigma_sigma*DDrot_dp2(9,1)  - pi_pi*(DDrot_dp2(9,2)  + DDrot_dp2(9,3))  + delta_delta*(DDrot_dp2(9,4)  + DDrot_dp2(9,5));             //<x2-y2|yz>
        overlaps(irw + 4,icl + 2) = sigma_sigma*DDrot_dp2(12,1) - pi_pi*(DDrot_dp2(12,2) + DDrot_dp2(12,3)) + delta_delta*(DDrot_dp2(12,4) + DDrot_dp2(12,5));            //<x2-y2|z2>
        overlaps(irw + 4,icl + 3) = sigma_sigma*DDrot_dp2(14,1) - pi_pi*(DDrot_dp2(14,2) + DDrot_dp2(14,3)) + delta_delta*(DDrot_dp2(14,4) + DDrot_dp2(14,5));            //<x2-y2|xy>
        overlaps(irw + 4,icl + 4) = sigma_sigma*DDrot_dp2(15,1) - pi_pi*(DDrot_dp2(15,2) + DDrot_dp2(15,3)) + delta_delta*(DDrot_dp2(15,4) + DDrot_dp2(15,5));            //<x2-y2|x2-y2>
      }
      else if ((iLA == 3)&&(iLB == 3)) {
        //ff
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = FFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = FFdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //phi integral
        aux = OvFactor(iNA + 1,iLA,3,a,iNB + 1,iLB,3,b,dist);
        phi_phi = 0.0;
        polynomial = FFphiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          phi_phi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 7x7 dimensions (fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        overlaps(irw + 0,icl + 0) = -sigma_sigma*FFrot_dp2(1,1)  + pi_pi*(FFrot_dp2(1,2)  + FFrot_dp2(1,3))  - delta_delta*(FFrot_dp2(1,4)  + FFrot_dp2(1,5))  + phi_phi*(FFrot_dp2(1,6)  + FFrot_dp2(1,7));
        overlaps(irw + 0,icl + 1) = -sigma_sigma*FFrot_dp2(2,1)  + pi_pi*(FFrot_dp2(2,2)  + FFrot_dp2(2,3))  - delta_delta*(FFrot_dp2(2,4)  + FFrot_dp2(2,5))  + phi_phi*(FFrot_dp2(2,6)  + FFrot_dp2(2,7));
        overlaps(irw + 0,icl + 2) = -sigma_sigma*FFrot_dp2(3,1)  + pi_pi*(FFrot_dp2(3,2)  + FFrot_dp2(3,3))  - delta_delta*(FFrot_dp2(3,4)  + FFrot_dp2(3,5))  + phi_phi*(FFrot_dp2(3,6)  + FFrot_dp2(3,7));
        overlaps(irw + 0,icl + 3) = -sigma_sigma*FFrot_dp2(4,1)  + pi_pi*(FFrot_dp2(4,2)  + FFrot_dp2(4,3))  - delta_delta*(FFrot_dp2(4,4)  + FFrot_dp2(4,5))  + phi_phi*(FFrot_dp2(4,6)  + FFrot_dp2(4,7));
        overlaps(irw + 0,icl + 4) = -sigma_sigma*FFrot_dp2(5,1)  + pi_pi*(FFrot_dp2(5,2)  + FFrot_dp2(5,3))  - delta_delta*(FFrot_dp2(5,4)  + FFrot_dp2(5,5))  + phi_phi*(FFrot_dp2(5,6)  + FFrot_dp2(5,7));
        overlaps(irw + 0,icl + 5) = -sigma_sigma*FFrot_dp2(6,1)  + pi_pi*(FFrot_dp2(6,2)  + FFrot_dp2(6,3))  - delta_delta*(FFrot_dp2(6,4)  + FFrot_dp2(6,5))  + phi_phi*(FFrot_dp2(6,6)  + FFrot_dp2(6,7));
        overlaps(irw + 0,icl + 6) = -sigma_sigma*FFrot_dp2(7,1)  + pi_pi*(FFrot_dp2(7,2)  + FFrot_dp2(7,3))  - delta_delta*(FFrot_dp2(7,4)  + FFrot_dp2(7,5))  + phi_phi*(FFrot_dp2(7,6)  + FFrot_dp2(7,7));
        overlaps(irw + 1,icl + 0) = -sigma_sigma*FFrot_dp2(2,1)  + pi_pi*(FFrot_dp2(2,2)  + FFrot_dp2(2,3))  - delta_delta*(FFrot_dp2(2,4)  + FFrot_dp2(2,5))  + phi_phi*(FFrot_dp2(2,6)  + FFrot_dp2(2,7));
        overlaps(irw + 1,icl + 1) = -sigma_sigma*FFrot_dp2(8,1)  + pi_pi*(FFrot_dp2(8,2)  + FFrot_dp2(8,3))  - delta_delta*(FFrot_dp2(8,4)  + FFrot_dp2(8,5))  + phi_phi*(FFrot_dp2(8,6)  + FFrot_dp2(8,7));
        overlaps(irw + 1,icl + 2) = -sigma_sigma*FFrot_dp2(9,1)  + pi_pi*(FFrot_dp2(9,2)  + FFrot_dp2(9,3))  - delta_delta*(FFrot_dp2(9,4)  + FFrot_dp2(9,5))  + phi_phi*(FFrot_dp2(9,6)  + FFrot_dp2(9,7));
        overlaps(irw + 1,icl + 3) = -sigma_sigma*FFrot_dp2(10,1) + pi_pi*(FFrot_dp2(10,2) + FFrot_dp2(10,3)) - delta_delta*(FFrot_dp2(10,4) + FFrot_dp2(10,5)) + phi_phi*(FFrot_dp2(10,6) + FFrot_dp2(10,7));
        overlaps(irw + 1,icl + 4) = -sigma_sigma*FFrot_dp2(11,1) + pi_pi*(FFrot_dp2(11,2) + FFrot_dp2(11,3)) - delta_delta*(FFrot_dp2(11,4) + FFrot_dp2(11,5)) + phi_phi*(FFrot_dp2(11,6) + FFrot_dp2(11,7));
        overlaps(irw + 1,icl + 5) = -sigma_sigma*FFrot_dp2(12,1) + pi_pi*(FFrot_dp2(12,2) + FFrot_dp2(12,3)) - delta_delta*(FFrot_dp2(12,4) + FFrot_dp2(12,5)) + phi_phi*(FFrot_dp2(12,6) + FFrot_dp2(12,7));
        overlaps(irw + 1,icl + 6) = -sigma_sigma*FFrot_dp2(13,1) + pi_pi*(FFrot_dp2(13,2) + FFrot_dp2(13,3)) - delta_delta*(FFrot_dp2(13,4) + FFrot_dp2(13,5)) + phi_phi*(FFrot_dp2(13,6) + FFrot_dp2(13,7));
        overlaps(irw + 2,icl + 0) = -sigma_sigma*FFrot_dp2(3,1) +  pi_pi*(FFrot_dp2(3,2)  + FFrot_dp2(3,3))  - delta_delta*(FFrot_dp2(3,4)  + FFrot_dp2(3,5))  + phi_phi*(FFrot_dp2(3,6)  + FFrot_dp2(3,7));
        overlaps(irw + 2,icl + 1) = -sigma_sigma*FFrot_dp2(9,1) +  pi_pi*(FFrot_dp2(9,2)  + FFrot_dp2(9,3))  - delta_delta*(FFrot_dp2(9,4)  + FFrot_dp2(9,5))  + phi_phi*(FFrot_dp2(9,6)  + FFrot_dp2(9,7));
        overlaps(irw + 2,icl + 2) = -sigma_sigma*FFrot_dp2(14,1) + pi_pi*(FFrot_dp2(14,2) + FFrot_dp2(14,3)) - delta_delta*(FFrot_dp2(14,4) + FFrot_dp2(14,5)) + phi_phi*(FFrot_dp2(14,6) + FFrot_dp2(14,7));
        overlaps(irw + 2,icl + 3) = -sigma_sigma*FFrot_dp2(15,1) + pi_pi*(FFrot_dp2(15,2) + FFrot_dp2(15,3)) - delta_delta*(FFrot_dp2(15,4) + FFrot_dp2(15,5)) + phi_phi*(FFrot_dp2(15,6) + FFrot_dp2(15,7));
        overlaps(irw + 2,icl + 4) = -sigma_sigma*FFrot_dp2(16,1) + pi_pi*(FFrot_dp2(16,2) + FFrot_dp2(16,3)) - delta_delta*(FFrot_dp2(16,4) + FFrot_dp2(16,5)) + phi_phi*(FFrot_dp2(16,6) + FFrot_dp2(16,7));
        overlaps(irw + 2,icl + 5) = -sigma_sigma*FFrot_dp2(17,1) + pi_pi*(FFrot_dp2(17,2) + FFrot_dp2(17,3)) - delta_delta*(FFrot_dp2(17,4) + FFrot_dp2(17,5)) + phi_phi*(FFrot_dp2(17,6) + FFrot_dp2(17,7));
        overlaps(irw + 2,icl + 6) = -sigma_sigma*FFrot_dp2(18,1) + pi_pi*(FFrot_dp2(18,2) + FFrot_dp2(18,3)) - delta_delta*(FFrot_dp2(18,4) + FFrot_dp2(18,5)) + phi_phi*(FFrot_dp2(18,6) + FFrot_dp2(18,7));
        overlaps(irw + 3,icl + 0) = -sigma_sigma*FFrot_dp2(4,1) +  pi_pi*(FFrot_dp2(4,2)  + FFrot_dp2(4,3))  - delta_delta*(FFrot_dp2(4,4)  + FFrot_dp2(4,5))  + phi_phi*(FFrot_dp2(4,6)  + FFrot_dp2(4,7));
        overlaps(irw + 3,icl + 1) = -sigma_sigma*FFrot_dp2(10,1) + pi_pi*(FFrot_dp2(10,2) + FFrot_dp2(10,3)) - delta_delta*(FFrot_dp2(10,4) + FFrot_dp2(10,5)) + phi_phi*(FFrot_dp2(10,6) + FFrot_dp2(10,7));
        overlaps(irw + 3,icl + 2) = -sigma_sigma*FFrot_dp2(15,1) + pi_pi*(FFrot_dp2(15,2) + FFrot_dp2(15,3)) - delta_delta*(FFrot_dp2(15,4) + FFrot_dp2(15,5)) + phi_phi*(FFrot_dp2(15,6) + FFrot_dp2(15,7));
        overlaps(irw + 3,icl + 3) = -sigma_sigma*FFrot_dp2(19,1) + pi_pi*(FFrot_dp2(19,2) + FFrot_dp2(19,3)) - delta_delta*(FFrot_dp2(19,4) + FFrot_dp2(19,5)) + phi_phi*(FFrot_dp2(19,6) + FFrot_dp2(19,7));
        overlaps(irw + 3,icl + 4) = -sigma_sigma*FFrot_dp2(20,1) + pi_pi*(FFrot_dp2(20,2) + FFrot_dp2(20,3)) - delta_delta*(FFrot_dp2(20,4) + FFrot_dp2(20,5)) + phi_phi*(FFrot_dp2(20,6) + FFrot_dp2(20,7));
        overlaps(irw + 3,icl + 5) = -sigma_sigma*FFrot_dp2(21,1) + pi_pi*(FFrot_dp2(21,2) + FFrot_dp2(21,3)) - delta_delta*(FFrot_dp2(21,4) + FFrot_dp2(21,5)) + phi_phi*(FFrot_dp2(21,6) + FFrot_dp2(21,7));
        overlaps(irw + 3,icl + 6) = -sigma_sigma*FFrot_dp2(22,1) + pi_pi*(FFrot_dp2(22,2) + FFrot_dp2(22,3)) - delta_delta*(FFrot_dp2(22,4) + FFrot_dp2(22,5)) + phi_phi*(FFrot_dp2(22,6) + FFrot_dp2(22,7));
        overlaps(irw + 4,icl + 0) = -sigma_sigma*FFrot_dp2(5,1) +  pi_pi*(FFrot_dp2(5,2)  + FFrot_dp2(5,3))  - delta_delta*(FFrot_dp2(5,4)  + FFrot_dp2(5,5))  + phi_phi*(FFrot_dp2(5,6)  + FFrot_dp2(5,7));
        overlaps(irw + 4,icl + 1) = -sigma_sigma*FFrot_dp2(11,1) + pi_pi*(FFrot_dp2(11,2) + FFrot_dp2(11,3)) - delta_delta*(FFrot_dp2(11,4) + FFrot_dp2(11,5)) + phi_phi*(FFrot_dp2(11,6) + FFrot_dp2(11,7));
        overlaps(irw + 4,icl + 2) = -sigma_sigma*FFrot_dp2(16,1) + pi_pi*(FFrot_dp2(16,2) + FFrot_dp2(16,3)) - delta_delta*(FFrot_dp2(16,4) + FFrot_dp2(16,5)) + phi_phi*(FFrot_dp2(16,6) + FFrot_dp2(16,7));
        overlaps(irw + 4,icl + 3) = -sigma_sigma*FFrot_dp2(20,1) + pi_pi*(FFrot_dp2(20,2) + FFrot_dp2(20,3)) - delta_delta*(FFrot_dp2(20,4) + FFrot_dp2(20,5)) + phi_phi*(FFrot_dp2(20,6) + FFrot_dp2(20,7));
        overlaps(irw + 4,icl + 4) = -sigma_sigma*FFrot_dp2(23,1) + pi_pi*(FFrot_dp2(23,2) + FFrot_dp2(23,3)) - delta_delta*(FFrot_dp2(23,4) + FFrot_dp2(23,5)) + phi_phi*(FFrot_dp2(23,6) + FFrot_dp2(23,7));
        overlaps(irw + 4,icl + 5) = -sigma_sigma*FFrot_dp2(24,1) + pi_pi*(FFrot_dp2(24,2) + FFrot_dp2(24,3)) - delta_delta*(FFrot_dp2(24,4) + FFrot_dp2(24,5)) + phi_phi*(FFrot_dp2(24,6) + FFrot_dp2(24,7));
        overlaps(irw + 4,icl + 6) = -sigma_sigma*FFrot_dp2(25,1) + pi_pi*(FFrot_dp2(25,2) + FFrot_dp2(25,3)) - delta_delta*(FFrot_dp2(25,4) + FFrot_dp2(25,5)) + phi_phi*(FFrot_dp2(25,6) + FFrot_dp2(25,7));
        overlaps(irw + 5,icl + 0) = -sigma_sigma*FFrot_dp2(6,1) +  pi_pi*(FFrot_dp2(6,2)  + FFrot_dp2(6,3))  - delta_delta*(FFrot_dp2(6,4)  + FFrot_dp2(6,5))  + phi_phi*(FFrot_dp2(6,6)  + FFrot_dp2(6,7));
        overlaps(irw + 5,icl + 1) = -sigma_sigma*FFrot_dp2(12,1) + pi_pi*(FFrot_dp2(12,2) + FFrot_dp2(12,3)) - delta_delta*(FFrot_dp2(12,4) + FFrot_dp2(12,5)) + phi_phi*(FFrot_dp2(12,6) + FFrot_dp2(12,7));
        overlaps(irw + 5,icl + 2) = -sigma_sigma*FFrot_dp2(17,1) + pi_pi*(FFrot_dp2(17,2) + FFrot_dp2(17,3)) - delta_delta*(FFrot_dp2(17,4) + FFrot_dp2(17,5)) + phi_phi*(FFrot_dp2(17,6) + FFrot_dp2(17,7));
        overlaps(irw + 5,icl + 3) = -sigma_sigma*FFrot_dp2(21,1) + pi_pi*(FFrot_dp2(21,2) + FFrot_dp2(21,3)) - delta_delta*(FFrot_dp2(21,4) + FFrot_dp2(21,5)) + phi_phi*(FFrot_dp2(21,6) + FFrot_dp2(21,7));
        overlaps(irw + 5,icl + 4) = -sigma_sigma*FFrot_dp2(24,1) + pi_pi*(FFrot_dp2(24,2) + FFrot_dp2(24,3)) - delta_delta*(FFrot_dp2(24,4) + FFrot_dp2(24,5)) + phi_phi*(FFrot_dp2(24,6) + FFrot_dp2(24,7));
        overlaps(irw + 5,icl + 5) = -sigma_sigma*FFrot_dp2(26,1) + pi_pi*(FFrot_dp2(26,2) + FFrot_dp2(26,3)) - delta_delta*(FFrot_dp2(26,4) + FFrot_dp2(26,5)) + phi_phi*(FFrot_dp2(26,6) + FFrot_dp2(26,7));
        overlaps(irw + 5,icl + 6) = -sigma_sigma*FFrot_dp2(27,1) + pi_pi*(FFrot_dp2(27,2) + FFrot_dp2(27,3)) - delta_delta*(FFrot_dp2(27,4) + FFrot_dp2(27,5)) + phi_phi*(FFrot_dp2(27,6) + FFrot_dp2(27,7));
        overlaps(irw + 6,icl + 0) = -sigma_sigma*FFrot_dp2(7,1) +  pi_pi*(FFrot_dp2(7,2)  + FFrot_dp2(7,3))  - delta_delta*(FFrot_dp2(7,4)  + FFrot_dp2(7,5))  + phi_phi*(FFrot_dp2(7,6)  + FFrot_dp2(7,7));
        overlaps(irw + 6,icl + 1) = -sigma_sigma*FFrot_dp2(13,1) + pi_pi*(FFrot_dp2(13,2) + FFrot_dp2(13,3)) - delta_delta*(FFrot_dp2(13,4) + FFrot_dp2(13,5)) + phi_phi*(FFrot_dp2(13,6) + FFrot_dp2(13,7));
        overlaps(irw + 6,icl + 2) = -sigma_sigma*FFrot_dp2(18,1) + pi_pi*(FFrot_dp2(18,2) + FFrot_dp2(18,3)) - delta_delta*(FFrot_dp2(18,4) + FFrot_dp2(18,5)) + phi_phi*(FFrot_dp2(18,6) + FFrot_dp2(18,7));
        overlaps(irw + 6,icl + 3) = -sigma_sigma*FFrot_dp2(22,1) + pi_pi*(FFrot_dp2(22,2) + FFrot_dp2(22,3)) - delta_delta*(FFrot_dp2(22,4) + FFrot_dp2(22,5)) + phi_phi*(FFrot_dp2(22,6) + FFrot_dp2(22,7));
        overlaps(irw + 6,icl + 4) = -sigma_sigma*FFrot_dp2(25,1) + pi_pi*(FFrot_dp2(25,2) + FFrot_dp2(25,3)) - delta_delta*(FFrot_dp2(25,4) + FFrot_dp2(25,5)) + phi_phi*(FFrot_dp2(25,6) + FFrot_dp2(25,7));
        overlaps(irw + 6,icl + 5) = -sigma_sigma*FFrot_dp2(27,1) + pi_pi*(FFrot_dp2(27,2) + FFrot_dp2(27,3)) - delta_delta*(FFrot_dp2(27,4) + FFrot_dp2(27,5)) + phi_phi*(FFrot_dp2(27,6) + FFrot_dp2(27,7));
        overlaps(irw + 6,icl + 6) = -sigma_sigma*FFrot_dp2(28,1) + pi_pi*(FFrot_dp2(28,2) + FFrot_dp2(28,3)) - delta_delta*(FFrot_dp2(28,4) + FFrot_dp2(28,5)) + phi_phi*(FFrot_dp2(28,6) + FFrot_dp2(28,7));
      }
      else if ((iLA == 0)&&(iLB == 1)) {
        //sp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SPOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3 dimensions (px,py,pz)
        overlaps(irw,icl + 0) = -SProt_dp2(1,1)*sigma_sigma;
        overlaps(irw,icl + 1) = -SProt_dp2(2,1)*sigma_sigma;
        overlaps(irw,icl + 2) = -SProt_dp2(3,1)*sigma_sigma;
      }
      else if ((iLA == 1)&&(iLB == 0)) {
        //ps
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3 dimensions (px,py,pz)
        overlaps(irw + 0,icl) = SProt_dp2(1,1)*sigma_sigma;
        overlaps(irw + 1,icl) = SProt_dp2(2,1)*sigma_sigma;
        overlaps(irw + 2,icl) = SProt_dp2(3,1)*sigma_sigma;
      }
      else if ((iLA == 0)&&(iLB == 2)) {
        //sd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SDOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          overlaps(irw,icl + index) = SDrot_dp2(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 2)&&(iLB == 0)) {
        //ds
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          overlaps(irw + index,icl) = SDrot_dp2(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 0)&&(iLB == 3)) {
        //sf
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SFOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 7 dimensions (fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 7; ++index) {
          overlaps(irw,icl + index) = -SFrot_dp2(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 3)&&(iLB == 0)) {
        //fs
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 7; ++index) {
          overlaps(irw + index,icl) = SFrot_dp2(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 1)&&(iLB == 2)) {
        //pd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = PDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 3; ++index2) {
            overlaps(irw + index2,icl + index) = PDrot_dp2(1+index2+3*index,1)*sigma_sigma - (PDrot_dp2(1+index2+3*index,2) + PDrot_dp2(1+index2+3*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 2)&&(iLB == 1)) {
        //dp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = DPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 3; ++index2) {
            overlaps(irw + index,icl + index2) = -PDrot_dp2(1+index2+3*index,1)*sigma_sigma + (PDrot_dp2(1+index2+3*index,2) + PDrot_dp2(1+index2+3*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 1)&&(iLB == 3)) {
        //pf
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = PFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 3; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index,icl + index2) = -PFrot_dp2(1+index2+7*index,1)*sigma_sigma + (PFrot_dp2(1+index2+7*index,2) + PFrot_dp2(1+index2+7*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 3)&&(iLB == 1)) {
        //fp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = FPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 3; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index2,icl + index) = -PFrot_dp2(1+index2+7*index,1)*sigma_sigma + (PFrot_dp2(1+index2+7*index,2) + PFrot_dp2(1+index2+7*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 2)&&(iLB == 3)) {
        //df
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = DFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = DFdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (xz,yz,z2,xy,x2-y2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index,icl + index2) = -DFrot_dp2(1+index2+7*index,1)*sigma_sigma + (DFrot_dp2(1+index2+7*index,2) + DFrot_dp2(1+index2+7*index,3))*pi_pi - (DFrot_dp2(1+index2+7*index,4) + DFrot_dp2(1+index2+7*index,5))*delta_delta;
          }
        }
      }
      else if ((iLA == 3)&&(iLB == 2)) {
        //fd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = FDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = FDdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (xz,yz,z2,xy,x2-y2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index2,icl + index) = DFrot_dp2(1+index2+7*index,1)*sigma_sigma - (DFrot_dp2(1+index2+7*index,2) + DFrot_dp2(1+index2+7*index,3))*pi_pi + (DFrot_dp2(1+index2+7*index,4) + DFrot_dp2(1+index2+7*index,5))*delta_delta;
          }
        }
      }
      icl += 2*QN2(isket + 1,2) + 1;
    }
    irw += 2*QN1(isbra + 1,2) + 1;
  }
}
void TwoCenterOv_dtdp(matrixE & overlaps, size_t nao1, size_t nao2, double dist, std::vector<double> & exponents1, matrix<size_t> & QN1, std::vector<double> & exponents2, matrix<size_t> & QN2, std::vector<double> & orientation, bool rotate = true, double tolerance = 1e-6) {
  //implementation of second-derivatives of overlap matrix function with respect to theta and phi; only up to shell f
  //calculation of the derivative of the overlap function between atoms A and B (A -> 1; B -> 2) with respect to angle theta
  overlaps.resize(nao1,nao2);
  overlaps.zero();
  //auxiliary variables for integrals
  double sigma_sigma = 0.0;
  double pi_pi = 0.0;
  double delta_delta = 0.0;
  double phi_phi = 0.0;
  //angles needed to reorient the integrals according to the position of atoms in the molecule
  double cost = orientation[2];
  double sint = sqrt(1 - cost*cost);
  double cosp = 1.0;
  double sinp = 0.0;
  if (fabs(sint) > tolerance) {
    cosp = orientation[0]/sint;
    sinp = orientation[1]/sint;
  }
  //matrices with derivatives of overlap rotations with respect to theta
  matrixE SProt_dtdp(1,1);
  matrixE SDrot_dtdp(1,1);
  matrixE SFrot_dtdp(1,1);
  matrixE PProt_dtdp(1,1);
  matrixE PDrot_dtdp(1,1);
  matrixE PFrot_dtdp(1,1);
  matrixE DDrot_dtdp(1,1);
  matrixE DFrot_dtdp(1,1);
  matrixE FFrot_dtdp(1,1);
  //other auxiliaries
  std::vector<std::vector<int> > polynomial;
  double a;
  double b;
  double c;
  double argA;
  double argB;
  std::vector<double> Ak;
  std::vector<double> Bk;
  double aux;
  //determine maximum azimuthal quantum number
  size_t maxL1 = 0;
  size_t maxL2 = 0;
  for (size_t idx = 0; idx < QN1.rows(); ++idx) {
    if (QN1(idx + 1,2) > maxL1) {maxL1 = QN1(idx + 1,2);}
  }
  for (size_t idx = 0; idx < QN2.rows(); ++idx) {
    if (QN2(idx + 1,2) > maxL2) {maxL2 = QN2(idx + 1,2);}
  }
  if ((maxL1 > 0)||(maxL2 > 0)) {
    //in this case we have p orbitals
    SProt_dtdp = SPtransf_dtdp(cost,sint,cosp,sinp,rotate);
    PProt_dtdp = PPtransf_dtdp(cost,sint,cosp,sinp,rotate);
  }
  if ((maxL1 > 1)||(maxL2 > 1)) {
    //in this case we have d orbitals
    SDrot_dtdp = SDtransf_dtdp(cost,sint,cosp,sinp,rotate);
    PDrot_dtdp = PDtransf_dtdp(cost,sint,cosp,sinp,rotate);
    DDrot_dtdp = DDtransf_dtdp(cost,sint,cosp,sinp,rotate);
  }
  if ((maxL1 > 2)||(maxL2 > 2)) {
    //in this case we have f orbitals
    SFrot_dtdp = SFtransf_dtdp(cost,sint,cosp,sinp,rotate);
    PFrot_dtdp = PFtransf_dtdp(cost,sint,cosp,sinp,rotate);
    DFrot_dtdp = DFtransf_dtdp(cost,sint,cosp,sinp,rotate);
    FFrot_dtdp = FFtransf_dtdp(cost,sint,cosp,sinp,rotate);
  }
  size_t iNA;
  size_t iNB;
  size_t iLA;
  size_t iLB;
  size_t irw = 1;
  size_t icl = 1;
  for (size_t isbra = 0; isbra < QN1.rows(); ++isbra) {           //loop over bra shells
    a = exponents1[isbra];
    if (a == 0.0) continue;
    iNA = QN1(isbra + 1,1) - 1;
    iLA = QN1(isbra + 1,2);
    icl = 1;
    for (size_t isket = 0; isket < QN2.rows(); ++isket) {         //loop over ket shells
      b = exponents2[isket];
      if (b == 0.0) continue;
      iNB = QN2(isket + 1,1) - 1;
      iLB = QN2(isket + 1,2);
      //getting integrals
      argA = 0.5*dist*(a + b);
      argB = 0.5*dist*(a - b);
      Ak = int_1_infty_xpk_expmyx_dx(argA,iNA + iNB + 2);
      Bk = int_m1_1_xpk_expmyx_dx_Pople(argB,iNA + iNB + 2);
      aux = OvFactor(iNA + 1,iLA,0,a,iNB + 1,iLB,0,b,dist);
      //different azimuthal quantum number cases
      //ss = 0
      if ((iLA == 1)&&(iLB == 1)) {
        //pp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = PPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 9 dimensions
        overlaps(irw + 0,icl + 0) = -sigma_sigma*PProt_dtdp(1,1) + pi_pi*(PProt_dtdp(1,2) + PProt_dtdp(1,3));      //<x|x>
        overlaps(irw + 0,icl + 1) = -sigma_sigma*PProt_dtdp(2,1) + pi_pi*(PProt_dtdp(2,2) + PProt_dtdp(2,3));      //<x|y>
        overlaps(irw + 0,icl + 2) = -sigma_sigma*PProt_dtdp(3,1) + pi_pi*(PProt_dtdp(3,2) + PProt_dtdp(3,3));      //<x|z>
        overlaps(irw + 1,icl + 0) = -sigma_sigma*PProt_dtdp(2,1) + pi_pi*(PProt_dtdp(2,2) + PProt_dtdp(2,3));      //<y|x>
        overlaps(irw + 1,icl + 1) = -sigma_sigma*PProt_dtdp(4,1) + pi_pi*(PProt_dtdp(4,2) + PProt_dtdp(4,3));      //<y|y>
        overlaps(irw + 1,icl + 2) = -sigma_sigma*PProt_dtdp(5,1) + pi_pi*(PProt_dtdp(5,2) + PProt_dtdp(5,3));      //<y|z>
        overlaps(irw + 2,icl + 0) = -sigma_sigma*PProt_dtdp(3,1) + pi_pi*(PProt_dtdp(3,2) + PProt_dtdp(3,3));      //<z|x>
        overlaps(irw + 2,icl + 1) = -sigma_sigma*PProt_dtdp(5,1) + pi_pi*(PProt_dtdp(5,2) + PProt_dtdp(5,3));      //<z|y>
        overlaps(irw + 2,icl + 2) = -sigma_sigma*PProt_dtdp(6,1) + pi_pi*(PProt_dtdp(6,2) + PProt_dtdp(6,3));      //<z|z>
      }
      else if ((iLA == 2)&&(iLB == 2)) {
        //dd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = DDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = DDdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5x5 dimensions (xz,yz,z2,xy,x2-y2)(xz,yz,z2,xy,x2-y2)
        overlaps(irw + 0,icl + 0) = sigma_sigma*DDrot_dtdp(1,1)  - pi_pi*(DDrot_dtdp(1,2)  + DDrot_dtdp(1,3))  + delta_delta*(DDrot_dtdp(1,4)  + DDrot_dtdp(1,5));             //<xz|xz>
        overlaps(irw + 0,icl + 1) = sigma_sigma*DDrot_dtdp(2,1)  - pi_pi*(DDrot_dtdp(2,2)  + DDrot_dtdp(2,3))  + delta_delta*(DDrot_dtdp(2,4)  + DDrot_dtdp(2,5));             //<xz|yz>
        overlaps(irw + 0,icl + 2) = sigma_sigma*DDrot_dtdp(3,1)  - pi_pi*(DDrot_dtdp(3,2)  + DDrot_dtdp(3,3))  + delta_delta*(DDrot_dtdp(3,4)  + DDrot_dtdp(3,5));             //<xz|z2>
        overlaps(irw + 0,icl + 3) = sigma_sigma*DDrot_dtdp(4,1)  - pi_pi*(DDrot_dtdp(4,2)  + DDrot_dtdp(4,3))  + delta_delta*(DDrot_dtdp(4,4)  + DDrot_dtdp(4,5));             //<xz|xy>
        overlaps(irw + 0,icl + 4) = sigma_sigma*DDrot_dtdp(5,1)  - pi_pi*(DDrot_dtdp(5,2)  + DDrot_dtdp(5,3))  + delta_delta*(DDrot_dtdp(5,4)  + DDrot_dtdp(5,5));             //<xz|x2-y2>
        overlaps(irw + 1,icl + 0) = sigma_sigma*DDrot_dtdp(2,1)  - pi_pi*(DDrot_dtdp(2,2)  + DDrot_dtdp(2,3))  + delta_delta*(DDrot_dtdp(2,4)  + DDrot_dtdp(2,5));             //<yz|xz>
        overlaps(irw + 1,icl + 1) = sigma_sigma*DDrot_dtdp(6,1)  - pi_pi*(DDrot_dtdp(6,2)  + DDrot_dtdp(6,3))  + delta_delta*(DDrot_dtdp(6,4)  + DDrot_dtdp(6,5));             //<yz|yz>
        overlaps(irw + 1,icl + 2) = sigma_sigma*DDrot_dtdp(7,1)  - pi_pi*(DDrot_dtdp(7,2)  + DDrot_dtdp(7,3))  + delta_delta*(DDrot_dtdp(7,4)  + DDrot_dtdp(7,5));             //<yz|z2>
        overlaps(irw + 1,icl + 3) = sigma_sigma*DDrot_dtdp(8,1)  - pi_pi*(DDrot_dtdp(8,2)  + DDrot_dtdp(8,3))  + delta_delta*(DDrot_dtdp(8,4)  + DDrot_dtdp(8,5));             //<yz|xy>
        overlaps(irw + 1,icl + 4) = sigma_sigma*DDrot_dtdp(9,1)  - pi_pi*(DDrot_dtdp(9,2)  + DDrot_dtdp(9,3))  + delta_delta*(DDrot_dtdp(9,4)  + DDrot_dtdp(9,5));             //<yz|x2-y2>
        overlaps(irw + 2,icl + 0) = sigma_sigma*DDrot_dtdp(3,1)  - pi_pi*(DDrot_dtdp(3,2)  + DDrot_dtdp(3,3))  + delta_delta*(DDrot_dtdp(3,4)  + DDrot_dtdp(3,5));             //<z2|xz>
        overlaps(irw + 2,icl + 1) = sigma_sigma*DDrot_dtdp(7,1)  - pi_pi*(DDrot_dtdp(7,2)  + DDrot_dtdp(7,3))  + delta_delta*(DDrot_dtdp(7,4)  + DDrot_dtdp(7,5));             //<z2|yz>
        overlaps(irw + 2,icl + 2) = sigma_sigma*DDrot_dtdp(10,1) - pi_pi*(DDrot_dtdp(10,2) + DDrot_dtdp(10,3)) + delta_delta*(DDrot_dtdp(10,4) + DDrot_dtdp(10,5));            //<z2|z2>
        overlaps(irw + 2,icl + 3) = sigma_sigma*DDrot_dtdp(11,1) - pi_pi*(DDrot_dtdp(11,2) + DDrot_dtdp(11,3)) + delta_delta*(DDrot_dtdp(11,4) + DDrot_dtdp(11,5));            //<z2|xy>
        overlaps(irw + 2,icl + 4) = sigma_sigma*DDrot_dtdp(12,1) - pi_pi*(DDrot_dtdp(12,2) + DDrot_dtdp(12,3)) + delta_delta*(DDrot_dtdp(12,4) + DDrot_dtdp(12,5));            //<z2|x2-y2>
        overlaps(irw + 3,icl + 0) = sigma_sigma*DDrot_dtdp(4,1)  - pi_pi*(DDrot_dtdp(4,2)  + DDrot_dtdp(4,3))  + delta_delta*(DDrot_dtdp(4,4)  + DDrot_dtdp(4,5));             //<xy|xz>
        overlaps(irw + 3,icl + 1) = sigma_sigma*DDrot_dtdp(8,1)  - pi_pi*(DDrot_dtdp(8,2)  + DDrot_dtdp(8,3))  + delta_delta*(DDrot_dtdp(8,4)  + DDrot_dtdp(8,5));             //<xy|yz>
        overlaps(irw + 3,icl + 2) = sigma_sigma*DDrot_dtdp(11,1) - pi_pi*(DDrot_dtdp(11,2) + DDrot_dtdp(11,3)) + delta_delta*(DDrot_dtdp(11,4) + DDrot_dtdp(11,5));            //<xy|z2>
        overlaps(irw + 3,icl + 3) = sigma_sigma*DDrot_dtdp(13,1) - pi_pi*(DDrot_dtdp(13,2) + DDrot_dtdp(13,3)) + delta_delta*(DDrot_dtdp(13,4) + DDrot_dtdp(13,5));            //<xy|xy>
        overlaps(irw + 3,icl + 4) = sigma_sigma*DDrot_dtdp(14,1) - pi_pi*(DDrot_dtdp(14,2) + DDrot_dtdp(14,3)) + delta_delta*(DDrot_dtdp(14,4) + DDrot_dtdp(14,5));            //<xy|x2-y2>
        overlaps(irw + 4,icl + 0) = sigma_sigma*DDrot_dtdp(5,1)  - pi_pi*(DDrot_dtdp(5,2)  + DDrot_dtdp(5,3))  + delta_delta*(DDrot_dtdp(5,4)  + DDrot_dtdp(5,5));             //<x2-y2|xz>
        overlaps(irw + 4,icl + 1) = sigma_sigma*DDrot_dtdp(9,1)  - pi_pi*(DDrot_dtdp(9,2)  + DDrot_dtdp(9,3))  + delta_delta*(DDrot_dtdp(9,4)  + DDrot_dtdp(9,5));             //<x2-y2|yz>
        overlaps(irw + 4,icl + 2) = sigma_sigma*DDrot_dtdp(12,1) - pi_pi*(DDrot_dtdp(12,2) + DDrot_dtdp(12,3)) + delta_delta*(DDrot_dtdp(12,4) + DDrot_dtdp(12,5));            //<x2-y2|z2>
        overlaps(irw + 4,icl + 3) = sigma_sigma*DDrot_dtdp(14,1) - pi_pi*(DDrot_dtdp(14,2) + DDrot_dtdp(14,3)) + delta_delta*(DDrot_dtdp(14,4) + DDrot_dtdp(14,5));            //<x2-y2|xy>
        overlaps(irw + 4,icl + 4) = sigma_sigma*DDrot_dtdp(15,1) - pi_pi*(DDrot_dtdp(15,2) + DDrot_dtdp(15,3)) + delta_delta*(DDrot_dtdp(15,4) + DDrot_dtdp(15,5));            //<x2-y2|x2-y2>
      }
      else if ((iLA == 3)&&(iLB == 3)) {
        //ff
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = FFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = FFdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //phi integral
        aux = OvFactor(iNA + 1,iLA,3,a,iNB + 1,iLB,3,b,dist);
        phi_phi = 0.0;
        polynomial = FFphiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          phi_phi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 7x7 dimensions (fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        overlaps(irw + 0,icl + 0) = -sigma_sigma*FFrot_dtdp(1,1)  + pi_pi*(FFrot_dtdp(1,2)  + FFrot_dtdp(1,3))  - delta_delta*(FFrot_dtdp(1,4)  + FFrot_dtdp(1,5))  + phi_phi*(FFrot_dtdp(1,6)  + FFrot_dtdp(1,7));
        overlaps(irw + 0,icl + 1) = -sigma_sigma*FFrot_dtdp(2,1)  + pi_pi*(FFrot_dtdp(2,2)  + FFrot_dtdp(2,3))  - delta_delta*(FFrot_dtdp(2,4)  + FFrot_dtdp(2,5))  + phi_phi*(FFrot_dtdp(2,6)  + FFrot_dtdp(2,7));
        overlaps(irw + 0,icl + 2) = -sigma_sigma*FFrot_dtdp(3,1)  + pi_pi*(FFrot_dtdp(3,2)  + FFrot_dtdp(3,3))  - delta_delta*(FFrot_dtdp(3,4)  + FFrot_dtdp(3,5))  + phi_phi*(FFrot_dtdp(3,6)  + FFrot_dtdp(3,7));
        overlaps(irw + 0,icl + 3) = -sigma_sigma*FFrot_dtdp(4,1)  + pi_pi*(FFrot_dtdp(4,2)  + FFrot_dtdp(4,3))  - delta_delta*(FFrot_dtdp(4,4)  + FFrot_dtdp(4,5))  + phi_phi*(FFrot_dtdp(4,6)  + FFrot_dtdp(4,7));
        overlaps(irw + 0,icl + 4) = -sigma_sigma*FFrot_dtdp(5,1)  + pi_pi*(FFrot_dtdp(5,2)  + FFrot_dtdp(5,3))  - delta_delta*(FFrot_dtdp(5,4)  + FFrot_dtdp(5,5))  + phi_phi*(FFrot_dtdp(5,6)  + FFrot_dtdp(5,7));
        overlaps(irw + 0,icl + 5) = -sigma_sigma*FFrot_dtdp(6,1)  + pi_pi*(FFrot_dtdp(6,2)  + FFrot_dtdp(6,3))  - delta_delta*(FFrot_dtdp(6,4)  + FFrot_dtdp(6,5))  + phi_phi*(FFrot_dtdp(6,6)  + FFrot_dtdp(6,7));
        overlaps(irw + 0,icl + 6) = -sigma_sigma*FFrot_dtdp(7,1)  + pi_pi*(FFrot_dtdp(7,2)  + FFrot_dtdp(7,3))  - delta_delta*(FFrot_dtdp(7,4)  + FFrot_dtdp(7,5))  + phi_phi*(FFrot_dtdp(7,6)  + FFrot_dtdp(7,7));
        overlaps(irw + 1,icl + 0) = -sigma_sigma*FFrot_dtdp(2,1)  + pi_pi*(FFrot_dtdp(2,2)  + FFrot_dtdp(2,3))  - delta_delta*(FFrot_dtdp(2,4)  + FFrot_dtdp(2,5))  + phi_phi*(FFrot_dtdp(2,6)  + FFrot_dtdp(2,7));
        overlaps(irw + 1,icl + 1) = -sigma_sigma*FFrot_dtdp(8,1)  + pi_pi*(FFrot_dtdp(8,2)  + FFrot_dtdp(8,3))  - delta_delta*(FFrot_dtdp(8,4)  + FFrot_dtdp(8,5))  + phi_phi*(FFrot_dtdp(8,6)  + FFrot_dtdp(8,7));
        overlaps(irw + 1,icl + 2) = -sigma_sigma*FFrot_dtdp(9,1)  + pi_pi*(FFrot_dtdp(9,2)  + FFrot_dtdp(9,3))  - delta_delta*(FFrot_dtdp(9,4)  + FFrot_dtdp(9,5))  + phi_phi*(FFrot_dtdp(9,6)  + FFrot_dtdp(9,7));
        overlaps(irw + 1,icl + 3) = -sigma_sigma*FFrot_dtdp(10,1) + pi_pi*(FFrot_dtdp(10,2) + FFrot_dtdp(10,3)) - delta_delta*(FFrot_dtdp(10,4) + FFrot_dtdp(10,5)) + phi_phi*(FFrot_dtdp(10,6) + FFrot_dtdp(10,7));
        overlaps(irw + 1,icl + 4) = -sigma_sigma*FFrot_dtdp(11,1) + pi_pi*(FFrot_dtdp(11,2) + FFrot_dtdp(11,3)) - delta_delta*(FFrot_dtdp(11,4) + FFrot_dtdp(11,5)) + phi_phi*(FFrot_dtdp(11,6) + FFrot_dtdp(11,7));
        overlaps(irw + 1,icl + 5) = -sigma_sigma*FFrot_dtdp(12,1) + pi_pi*(FFrot_dtdp(12,2) + FFrot_dtdp(12,3)) - delta_delta*(FFrot_dtdp(12,4) + FFrot_dtdp(12,5)) + phi_phi*(FFrot_dtdp(12,6) + FFrot_dtdp(12,7));
        overlaps(irw + 1,icl + 6) = -sigma_sigma*FFrot_dtdp(13,1) + pi_pi*(FFrot_dtdp(13,2) + FFrot_dtdp(13,3)) - delta_delta*(FFrot_dtdp(13,4) + FFrot_dtdp(13,5)) + phi_phi*(FFrot_dtdp(13,6) + FFrot_dtdp(13,7));
        overlaps(irw + 2,icl + 0) = -sigma_sigma*FFrot_dtdp(3,1) +  pi_pi*(FFrot_dtdp(3,2)  + FFrot_dtdp(3,3))  - delta_delta*(FFrot_dtdp(3,4)  + FFrot_dtdp(3,5))  + phi_phi*(FFrot_dtdp(3,6)  + FFrot_dtdp(3,7));
        overlaps(irw + 2,icl + 1) = -sigma_sigma*FFrot_dtdp(9,1) +  pi_pi*(FFrot_dtdp(9,2)  + FFrot_dtdp(9,3))  - delta_delta*(FFrot_dtdp(9,4)  + FFrot_dtdp(9,5))  + phi_phi*(FFrot_dtdp(9,6)  + FFrot_dtdp(9,7));
        overlaps(irw + 2,icl + 2) = -sigma_sigma*FFrot_dtdp(14,1) + pi_pi*(FFrot_dtdp(14,2) + FFrot_dtdp(14,3)) - delta_delta*(FFrot_dtdp(14,4) + FFrot_dtdp(14,5)) + phi_phi*(FFrot_dtdp(14,6) + FFrot_dtdp(14,7));
        overlaps(irw + 2,icl + 3) = -sigma_sigma*FFrot_dtdp(15,1) + pi_pi*(FFrot_dtdp(15,2) + FFrot_dtdp(15,3)) - delta_delta*(FFrot_dtdp(15,4) + FFrot_dtdp(15,5)) + phi_phi*(FFrot_dtdp(15,6) + FFrot_dtdp(15,7));
        overlaps(irw + 2,icl + 4) = -sigma_sigma*FFrot_dtdp(16,1) + pi_pi*(FFrot_dtdp(16,2) + FFrot_dtdp(16,3)) - delta_delta*(FFrot_dtdp(16,4) + FFrot_dtdp(16,5)) + phi_phi*(FFrot_dtdp(16,6) + FFrot_dtdp(16,7));
        overlaps(irw + 2,icl + 5) = -sigma_sigma*FFrot_dtdp(17,1) + pi_pi*(FFrot_dtdp(17,2) + FFrot_dtdp(17,3)) - delta_delta*(FFrot_dtdp(17,4) + FFrot_dtdp(17,5)) + phi_phi*(FFrot_dtdp(17,6) + FFrot_dtdp(17,7));
        overlaps(irw + 2,icl + 6) = -sigma_sigma*FFrot_dtdp(18,1) + pi_pi*(FFrot_dtdp(18,2) + FFrot_dtdp(18,3)) - delta_delta*(FFrot_dtdp(18,4) + FFrot_dtdp(18,5)) + phi_phi*(FFrot_dtdp(18,6) + FFrot_dtdp(18,7));
        overlaps(irw + 3,icl + 0) = -sigma_sigma*FFrot_dtdp(4,1) +  pi_pi*(FFrot_dtdp(4,2)  + FFrot_dtdp(4,3))  - delta_delta*(FFrot_dtdp(4,4)  + FFrot_dtdp(4,5))  + phi_phi*(FFrot_dtdp(4,6)  + FFrot_dtdp(4,7));
        overlaps(irw + 3,icl + 1) = -sigma_sigma*FFrot_dtdp(10,1) + pi_pi*(FFrot_dtdp(10,2) + FFrot_dtdp(10,3)) - delta_delta*(FFrot_dtdp(10,4) + FFrot_dtdp(10,5)) + phi_phi*(FFrot_dtdp(10,6) + FFrot_dtdp(10,7));
        overlaps(irw + 3,icl + 2) = -sigma_sigma*FFrot_dtdp(15,1) + pi_pi*(FFrot_dtdp(15,2) + FFrot_dtdp(15,3)) - delta_delta*(FFrot_dtdp(15,4) + FFrot_dtdp(15,5)) + phi_phi*(FFrot_dtdp(15,6) + FFrot_dtdp(15,7));
        overlaps(irw + 3,icl + 3) = -sigma_sigma*FFrot_dtdp(19,1) + pi_pi*(FFrot_dtdp(19,2) + FFrot_dtdp(19,3)) - delta_delta*(FFrot_dtdp(19,4) + FFrot_dtdp(19,5)) + phi_phi*(FFrot_dtdp(19,6) + FFrot_dtdp(19,7));
        overlaps(irw + 3,icl + 4) = -sigma_sigma*FFrot_dtdp(20,1) + pi_pi*(FFrot_dtdp(20,2) + FFrot_dtdp(20,3)) - delta_delta*(FFrot_dtdp(20,4) + FFrot_dtdp(20,5)) + phi_phi*(FFrot_dtdp(20,6) + FFrot_dtdp(20,7));
        overlaps(irw + 3,icl + 5) = -sigma_sigma*FFrot_dtdp(21,1) + pi_pi*(FFrot_dtdp(21,2) + FFrot_dtdp(21,3)) - delta_delta*(FFrot_dtdp(21,4) + FFrot_dtdp(21,5)) + phi_phi*(FFrot_dtdp(21,6) + FFrot_dtdp(21,7));
        overlaps(irw + 3,icl + 6) = -sigma_sigma*FFrot_dtdp(22,1) + pi_pi*(FFrot_dtdp(22,2) + FFrot_dtdp(22,3)) - delta_delta*(FFrot_dtdp(22,4) + FFrot_dtdp(22,5)) + phi_phi*(FFrot_dtdp(22,6) + FFrot_dtdp(22,7));
        overlaps(irw + 4,icl + 0) = -sigma_sigma*FFrot_dtdp(5,1) +  pi_pi*(FFrot_dtdp(5,2)  + FFrot_dtdp(5,3))  - delta_delta*(FFrot_dtdp(5,4)  + FFrot_dtdp(5,5))  + phi_phi*(FFrot_dtdp(5,6)  + FFrot_dtdp(5,7));
        overlaps(irw + 4,icl + 1) = -sigma_sigma*FFrot_dtdp(11,1) + pi_pi*(FFrot_dtdp(11,2) + FFrot_dtdp(11,3)) - delta_delta*(FFrot_dtdp(11,4) + FFrot_dtdp(11,5)) + phi_phi*(FFrot_dtdp(11,6) + FFrot_dtdp(11,7));
        overlaps(irw + 4,icl + 2) = -sigma_sigma*FFrot_dtdp(16,1) + pi_pi*(FFrot_dtdp(16,2) + FFrot_dtdp(16,3)) - delta_delta*(FFrot_dtdp(16,4) + FFrot_dtdp(16,5)) + phi_phi*(FFrot_dtdp(16,6) + FFrot_dtdp(16,7));
        overlaps(irw + 4,icl + 3) = -sigma_sigma*FFrot_dtdp(20,1) + pi_pi*(FFrot_dtdp(20,2) + FFrot_dtdp(20,3)) - delta_delta*(FFrot_dtdp(20,4) + FFrot_dtdp(20,5)) + phi_phi*(FFrot_dtdp(20,6) + FFrot_dtdp(20,7));
        overlaps(irw + 4,icl + 4) = -sigma_sigma*FFrot_dtdp(23,1) + pi_pi*(FFrot_dtdp(23,2) + FFrot_dtdp(23,3)) - delta_delta*(FFrot_dtdp(23,4) + FFrot_dtdp(23,5)) + phi_phi*(FFrot_dtdp(23,6) + FFrot_dtdp(23,7));
        overlaps(irw + 4,icl + 5) = -sigma_sigma*FFrot_dtdp(24,1) + pi_pi*(FFrot_dtdp(24,2) + FFrot_dtdp(24,3)) - delta_delta*(FFrot_dtdp(24,4) + FFrot_dtdp(24,5)) + phi_phi*(FFrot_dtdp(24,6) + FFrot_dtdp(24,7));
        overlaps(irw + 4,icl + 6) = -sigma_sigma*FFrot_dtdp(25,1) + pi_pi*(FFrot_dtdp(25,2) + FFrot_dtdp(25,3)) - delta_delta*(FFrot_dtdp(25,4) + FFrot_dtdp(25,5)) + phi_phi*(FFrot_dtdp(25,6) + FFrot_dtdp(25,7));
        overlaps(irw + 5,icl + 0) = -sigma_sigma*FFrot_dtdp(6,1) +  pi_pi*(FFrot_dtdp(6,2)  + FFrot_dtdp(6,3))  - delta_delta*(FFrot_dtdp(6,4)  + FFrot_dtdp(6,5))  + phi_phi*(FFrot_dtdp(6,6)  + FFrot_dtdp(6,7));
        overlaps(irw + 5,icl + 1) = -sigma_sigma*FFrot_dtdp(12,1) + pi_pi*(FFrot_dtdp(12,2) + FFrot_dtdp(12,3)) - delta_delta*(FFrot_dtdp(12,4) + FFrot_dtdp(12,5)) + phi_phi*(FFrot_dtdp(12,6) + FFrot_dtdp(12,7));
        overlaps(irw + 5,icl + 2) = -sigma_sigma*FFrot_dtdp(17,1) + pi_pi*(FFrot_dtdp(17,2) + FFrot_dtdp(17,3)) - delta_delta*(FFrot_dtdp(17,4) + FFrot_dtdp(17,5)) + phi_phi*(FFrot_dtdp(17,6) + FFrot_dtdp(17,7));
        overlaps(irw + 5,icl + 3) = -sigma_sigma*FFrot_dtdp(21,1) + pi_pi*(FFrot_dtdp(21,2) + FFrot_dtdp(21,3)) - delta_delta*(FFrot_dtdp(21,4) + FFrot_dtdp(21,5)) + phi_phi*(FFrot_dtdp(21,6) + FFrot_dtdp(21,7));
        overlaps(irw + 5,icl + 4) = -sigma_sigma*FFrot_dtdp(24,1) + pi_pi*(FFrot_dtdp(24,2) + FFrot_dtdp(24,3)) - delta_delta*(FFrot_dtdp(24,4) + FFrot_dtdp(24,5)) + phi_phi*(FFrot_dtdp(24,6) + FFrot_dtdp(24,7));
        overlaps(irw + 5,icl + 5) = -sigma_sigma*FFrot_dtdp(26,1) + pi_pi*(FFrot_dtdp(26,2) + FFrot_dtdp(26,3)) - delta_delta*(FFrot_dtdp(26,4) + FFrot_dtdp(26,5)) + phi_phi*(FFrot_dtdp(26,6) + FFrot_dtdp(26,7));
        overlaps(irw + 5,icl + 6) = -sigma_sigma*FFrot_dtdp(27,1) + pi_pi*(FFrot_dtdp(27,2) + FFrot_dtdp(27,3)) - delta_delta*(FFrot_dtdp(27,4) + FFrot_dtdp(27,5)) + phi_phi*(FFrot_dtdp(27,6) + FFrot_dtdp(27,7));
        overlaps(irw + 6,icl + 0) = -sigma_sigma*FFrot_dtdp(7,1) +  pi_pi*(FFrot_dtdp(7,2)  + FFrot_dtdp(7,3))  - delta_delta*(FFrot_dtdp(7,4)  + FFrot_dtdp(7,5))  + phi_phi*(FFrot_dtdp(7,6)  + FFrot_dtdp(7,7));
        overlaps(irw + 6,icl + 1) = -sigma_sigma*FFrot_dtdp(13,1) + pi_pi*(FFrot_dtdp(13,2) + FFrot_dtdp(13,3)) - delta_delta*(FFrot_dtdp(13,4) + FFrot_dtdp(13,5)) + phi_phi*(FFrot_dtdp(13,6) + FFrot_dtdp(13,7));
        overlaps(irw + 6,icl + 2) = -sigma_sigma*FFrot_dtdp(18,1) + pi_pi*(FFrot_dtdp(18,2) + FFrot_dtdp(18,3)) - delta_delta*(FFrot_dtdp(18,4) + FFrot_dtdp(18,5)) + phi_phi*(FFrot_dtdp(18,6) + FFrot_dtdp(18,7));
        overlaps(irw + 6,icl + 3) = -sigma_sigma*FFrot_dtdp(22,1) + pi_pi*(FFrot_dtdp(22,2) + FFrot_dtdp(22,3)) - delta_delta*(FFrot_dtdp(22,4) + FFrot_dtdp(22,5)) + phi_phi*(FFrot_dtdp(22,6) + FFrot_dtdp(22,7));
        overlaps(irw + 6,icl + 4) = -sigma_sigma*FFrot_dtdp(25,1) + pi_pi*(FFrot_dtdp(25,2) + FFrot_dtdp(25,3)) - delta_delta*(FFrot_dtdp(25,4) + FFrot_dtdp(25,5)) + phi_phi*(FFrot_dtdp(25,6) + FFrot_dtdp(25,7));
        overlaps(irw + 6,icl + 5) = -sigma_sigma*FFrot_dtdp(27,1) + pi_pi*(FFrot_dtdp(27,2) + FFrot_dtdp(27,3)) - delta_delta*(FFrot_dtdp(27,4) + FFrot_dtdp(27,5)) + phi_phi*(FFrot_dtdp(27,6) + FFrot_dtdp(27,7));
        overlaps(irw + 6,icl + 6) = -sigma_sigma*FFrot_dtdp(28,1) + pi_pi*(FFrot_dtdp(28,2) + FFrot_dtdp(28,3)) - delta_delta*(FFrot_dtdp(28,4) + FFrot_dtdp(28,5)) + phi_phi*(FFrot_dtdp(28,6) + FFrot_dtdp(28,7));
      }
      else if ((iLA == 0)&&(iLB == 1)) {
        //sp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SPOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3 dimensions (px,py,pz)
        overlaps(irw,icl + 0) = -SProt_dtdp(1,1)*sigma_sigma;
        overlaps(irw,icl + 1) = -SProt_dtdp(2,1)*sigma_sigma;
        overlaps(irw,icl + 2) = -SProt_dtdp(3,1)*sigma_sigma;
      }
      else if ((iLA == 1)&&(iLB == 0)) {
        //ps
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3 dimensions (px,py,pz)
        overlaps(irw + 0,icl) = SProt_dtdp(1,1)*sigma_sigma;
        overlaps(irw + 1,icl) = SProt_dtdp(2,1)*sigma_sigma;
        overlaps(irw + 2,icl) = SProt_dtdp(3,1)*sigma_sigma;
      }
      else if ((iLA == 0)&&(iLB == 2)) {
        //sd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SDOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          overlaps(irw,icl + index) = SDrot_dtdp(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 2)&&(iLB == 0)) {
        //ds
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          overlaps(irw + index,icl) = SDrot_dtdp(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 0)&&(iLB == 3)) {
        //sf
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = SFOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 7 dimensions (fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 7; ++index) {
          overlaps(irw,icl + index) = -SFrot_dtdp(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 3)&&(iLB == 0)) {
        //fs
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FSOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 5 dimensions (xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 7; ++index) {
          overlaps(irw + index,icl) = SFrot_dtdp(index + 1,1)*sigma_sigma;
        }
      }
      else if ((iLA == 1)&&(iLB == 2)) {
        //pd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = PDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 3; ++index2) {
            overlaps(irw + index2,icl + index) = PDrot_dtdp(1+index2+3*index,1)*sigma_sigma - (PDrot_dtdp(1+index2+3*index,2) + PDrot_dtdp(1+index2+3*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 2)&&(iLB == 1)) {
        //dp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = DPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(xz,yz,z2,xy,x2-y2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 3; ++index2) {
            overlaps(irw + index,icl + index2) = -PDrot_dtdp(1+index2+3*index,1)*sigma_sigma + (PDrot_dtdp(1+index2+3*index,2) + PDrot_dtdp(1+index2+3*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 1)&&(iLB == 3)) {
        //pf
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = PFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = PFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 3; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index,icl + index2) = -PFrot_dtdp(1+index2+7*index,1)*sigma_sigma + (PFrot_dtdp(1+index2+7*index,2) + PFrot_dtdp(1+index2+7*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 3)&&(iLB == 1)) {
        //fp
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FPsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = FPpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (px,py,pz)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 3; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index2,icl + index) = -PFrot_dtdp(1+index2+7*index,1)*sigma_sigma + (PFrot_dtdp(1+index2+7*index,2) + PFrot_dtdp(1+index2+7*index,3))*pi_pi;
          }
        }
      }
      else if ((iLA == 2)&&(iLB == 3)) {
        //df
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = DFsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = DFpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = DFdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (xz,yz,z2,xy,x2-y2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index,icl + index2) = -DFrot_dtdp(1+index2+7*index,1)*sigma_sigma + (DFrot_dtdp(1+index2+7*index,2) + DFrot_dtdp(1+index2+7*index,3))*pi_pi - (DFrot_dtdp(1+index2+7*index,4) + DFrot_dtdp(1+index2+7*index,5))*delta_delta;
          }
        }
      }
      else if ((iLA == 3)&&(iLB == 2)) {
        //fd
        //sigma integral
        sigma_sigma = 0.0;
        polynomial = FDsigmaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          sigma_sigma += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //pi integral
        aux = OvFactor(iNA + 1,iLA,1,a,iNB + 1,iLB,1,b,dist);
        pi_pi = 0.0;
        polynomial = FDpiOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          pi_pi += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //delta integral
        aux = OvFactor(iNA + 1,iLA,2,a,iNB + 1,iLB,2,b,dist);
        delta_delta = 0.0;
        polynomial = FDdeltaOvIntIndex(iNA + 1,iNB + 1);
        for (size_t idx = 0; idx < polynomial[0].size(); ++idx) {
          delta_delta += aux*double(polynomial[0][idx])*Ak[polynomial[1][idx]]*Bk[polynomial[2][idx]];
        }
        //partition over 3x5 dimensions (xz,yz,z2,xy,x2-y2)(fz(x2-y2),fz(x2-y2),fy(3x2-y2),fxyz,fz3,fxz2,fyz2)
        for (size_t index = 0; index < 5; ++index) {
          for (size_t index2 = 0; index2 < 7; ++index2) {
            overlaps(irw + index2,icl + index) = DFrot_dtdp(1+index2+7*index,1)*sigma_sigma - (DFrot_dtdp(1+index2+7*index,2) + DFrot_dtdp(1+index2+7*index,3))*pi_pi + (DFrot_dtdp(1+index2+7*index,4) + DFrot_dtdp(1+index2+7*index,5))*delta_delta;
          }
        }
      }
      icl += 2*QN2(isket + 1,2) + 1;
    }
    irw += 2*QN1(isbra + 1,2) + 1;
  }
}
#endif //_Two_Center_Overlaps_STO_
