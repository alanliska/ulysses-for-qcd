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

#ifndef _Kinetics_Package_
#define _Kinetics_Package_
#include <vector>
#include "ConstantsPackage.hpp"

//description:
//functions for the calculation of kinetic constants

double EyringRateT(double T, double H, double S = 0, double c = 1) {
  //function that calculates the TST reaction rate as defined by Eyring in thermodynamic formulation
  //Eyring, J. Chem. Phys., 3, 107, 1935
  //     T = temperature in Kelvin
  //     H = DH≠  (activation enthalpy change) -> {H}_{TS} - {\sum}_{i}{{H}_{i}} in J/mol
  //     S = DS≠  (activation entropy change) -> {S}_{TS} - {\sum}_{i}{{S}_{i}} in J/[K.mol]
  //     c is Eyring correction factor to account for activated complexes that actually do not react (c ~ 1)
  //if S = 0 then it is assumed that H is the activation Gibbs energy change
  double kr = (c*KB*T/h)*exp(S/R)*exp(-H/(R*T));
  return kr;
}
double WignerRateT(double T, double omega, double H, double S = 0) {
  //function that calculates the TST reaction rate as defined by Eyring in thermodynamic formulation
  //and corrects for quantal effects using Wigner's transmission coefficient
  //Garrett, Truhlar, J. Phys. Chem., 83(1), 200, 1979
  //     omega is the imaginary vibrational frequency (which is here represented by a real number)
  double om = fabs(100.0*c*omega);
  double wigner = 1.0 + (h*om)*(h*om)/24.0;
  //hbar or simply h?
  double kr = EyringRateT(T,H,S,wigner);
  return kr;
}

#endif //_Kinetics_Package_
