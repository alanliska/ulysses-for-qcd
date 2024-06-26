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

#ifndef _Gaussian_Package_
#define _Gaussian_Package_
#include <math.h>
#include <vector>
#include <string>
#include "../ConstantsPackage.hpp"

//description:
//stuff related to gaussian functions

std::vector<double> GPTheorem(double alpha, std::vector<double> & RA, double beta, std::vector<double> & RB) {
  //function thst applies gaussian product theorem
  //this is a C++ translation of gpcenter in xTB
  size_t ndim = RA.size();
  if (ndim != RB.size()) {throw("ERROR: GaussianPackage.hpp: GPTheorem(): Gaussians with different dimensions");}
  std::vector<double> RC(ndim,0.0);
  double ac = alpha/(alpha + beta);
  double bc = beta/(alpha + beta);
  for (size_t idx = 0; idx < ndim; ++idx) {
    RC[idx] = RA[idx]*ac + RB[idx]*bc;
  }
  return RC;
}
void BuildKAB(double alpha, std::vector<double> & RA, double beta, std::vector<double> & RB, double & gamma, double & kab) {
  //function to compute the center, exponent, and multiplying factor of a single gaussian that replaces the product of two gaussian centers RA and RB, and exponents alpha and beta
  //this is a C++ translation of build_kab in xTB
  size_t ndim = RA.size();
  if (ndim != RB.size()) {throw("ERROR: GaussianPackage.hpp: BuildKAB(): Gaussians with different dimensions");}
  gamma = alpha + beta;
  double gm1 = 1.0/gamma;
  double est;
  double RAB2 = 0.0;
  for (size_t idcoord = 0; idcoord < ndim; ++idcoord) {
    est = RA[idcoord] - RB[idcoord];
    RAB2 += est*est;
  }
  est = RAB2*alpha*beta*gm1;
  RAB2 = sqrt(pi*gm1);
  kab = exp(-est)*RAB2*RAB2*RAB2;
}
double Gauss(double alpha, double xx) {
  return exp(-alpha*xx*xx);
}

#endif //_Gaussian_Package_
