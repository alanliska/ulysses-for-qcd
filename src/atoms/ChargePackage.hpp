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

#ifndef _Charge_Package_
#define _Charge_Package_
#include <vector>
#include <string>
#include "../parameters/Chargepar.hpp"

//descrition:
//functions that may help calculate charges of atoms in molecules

double GBEffectiveCoulombStill(double RAB2, double RA, double RB) {
  //the effective Coulomb operator of Still etal as in
  //W. C. Still, A. Tempczyk, R. C. Hawley, T. Hendrickson, J. Am. Chem. Soc., 112, 6127, 1990
  //the default value for the constant c (cc) is 4.0, as recommended by Still etal
  //other values are 
  //1.64  -> B. Jayaram, Y. Liu, D. L. Beveridge, J. Chem. Phys., 109, 1465, 1998
  //3.7   -> A. V. Marenich, C. J. Cramer, D. G. Truhlar, J. Chem. Theo. Comput., 5, 2447, 2009
  //5.249 -> J. A. Grant, B. T. Pickup, M. J. Sykes, C. A. Kitchen, A. Nicholls, Phys. Chem. Chem. Phys., 9, 4913, 2007
  //8.0   -> M. S. Lee, M. Feig, J. F. R. Salsbury, C. L. Brooks, J. Comput. Chem., 24, 1348, 2003
  //RAB2 is the square of distance between atoms A and B, whereas RA and RB are the Born radii
  double cc = 4.0;
  double RARB = RA*RB;
  double Bij = RAB2/(cc*RARB);
  return sqrt(RAB2 + RARB*exp(-Bij));
}
void GBEffectiveCoulombStill_dX(double RAB2, double RA, double RB, double & dRAB, double & dRI, double keps) {
//void GBEffectiveCoulombStill_dX(double RAB2, double RA, double RB, double & dRAB, double & dRI, double & ionder, double keps, double kappa, double ionicterm) {
  //the derivatives of the effective Coulomb operator of Still etal as in
  //W. C. Still, A. Tempczyk, R. C. Hawley, T. Hendrickson, J. Am. Chem. Soc., 112, 6127, 1990
  //the default value for the constant c (cc) is 4.0, as recommended by Still etal
  //other values are 
  //1.64  -> B. Jayaram, Y. Liu, D. L. Beveridge, J. Chem. Phys., 109, 1465, 1998
  //3.7   -> A. V. Marenich, C. J. Cramer, D. G. Truhlar, J. Chem. Theo. Comput., 5, 2447, 2009
  //5.249 -> J. A. Grant, B. T. Pickup, M. J. Sykes, C. A. Kitchen, A. Nicholls, Phys. Chem. Chem. Phys., 9, 4913, 2007
  //8.0   -> M. S. Lee, M. Feig, J. F. R. Salsbury, C. L. Brooks, J. Comput. Chem., 24, 1348, 2003
  //RAB2 is the square of distance between atoms A and B, whereas RA and RB are the Born radii
  double ccd1 = 0.25;
  double RARB = RA*RB;
  double exponent = ccd1*RAB2/RARB;
  double expd = exp(-exponent);
  double dfgb2 = 1.0/(RAB2 + RARB*expd);
  double dfgb3 = dfgb2*sqrt(dfgb2)*keps;
  dRAB = (1.0 - ccd1*expd)*dfgb3;
  dRI = -0.5*expd*(1.0 + exponent)*dfgb3;
}
double GBEffectiveCoulombP16(double RAB2, double RA, double RB) {
  //the effective Coulomb operator of Lange and Herbert as in equations 4.11 and 4.12
  //A. W. Lange, J. M. Herbert, J. Chem. Theo. Comput., 8, 1999, 2012
  //RAB2 is the square of distance between atoms A and B, whereas RA and RB are the Born radii
  double zeta = 1.028;
  double sqrtRARB = sqrt(RA*RB);
  double RAB = sqrt(RAB2);
  double aux = sqrtRARB/(sqrtRARB + 0.0625*zeta*RAB);    //x in equation 5.7
  double res = aux*aux;                                  //x^2
  aux = res*res;                                         //x^4
  res = aux*aux;                                         //x^8
  aux = res*res;                                         //x^16 = {\Omega}_{AB}^{p16}, equation 4.12
  return RAB + aux*sqrtRARB;                             //equation 4.11
}
void GBEffectiveCoulombP16_dX(double RAB2, double RA, double RB, double & dRAB, double & dRI, double keps) {
//void GBEffectiveCoulombP16_dX(double RAB2, double RA, double RB, double & dRAB, double & dRI, double & ionder, double keps, double kappa, double ionicterm) {
  //the derivatives of the effective Coulomb operator of Lange and Herbert as in equations 4.11 and 4.12
  //A. W. Lange, J. M. Herbert, J. Chem. Theo. Comput., 8, 1999, 2012
  //RAB2 is the square of distance between atoms A and B, whereas RA and RB are the Born radii
  //dRAB -> for gradient contribution from atomic distance; still needs to be multiplied by distance vector and product of charges
  //dRI -> for gradient contribution from the Born radii; still needs multiplication by Born radii
  double zeta = 1.028;
  double sqrtRARB = sqrt(RA*RB);
  double RAB = sqrt(RAB2);
  double aux = sqrtRARB/(sqrtRARB + 0.0625*zeta*RAB);    //x in equation 5.7
  double singleprod = aux;
  double res = aux*aux;                                  //x^2
  aux = res*res;                                         //x^4
  res = aux*aux;                                         //x^8
  aux = res*res;                                         //x^16 = {\Omega}_{AB}^{p16}, equation 4.12
  res = 1.0/(RAB + aux*sqrtRARB);
  dRI = -0.5*keps*(RAB*zeta*singleprod/sqrtRARB + 1.0)*aux*res*res/sqrtRARB;
  dRAB = (1.0 - zeta*singleprod*aux)*keps*res*res/RAB;
}

#endif //_Charge_Package_
