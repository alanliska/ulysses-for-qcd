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

#ifndef _Model_Potentials_
#define _Model_Potentials_
#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include "BSet.hpp"
#include "Molecule.hpp"
#include "ConstantsPackage.hpp"
#include "UnitConversion.hpp"
#include "math/FunctionPackage.hpp"

//description:
//classes to calculate energies according to a given model potential

class EllipsoidalConstraint {
  //class that evaluates constraint potential of ellipsoidal form
  //only parameters evaluated
  double radii[3];                        //ellipsoid's radii in Angstroem
  double maxrad;
  double anisotropy[3];
  double exponent;                        //exponent in the case of power potential
  double beta;                            //beta for log-Fermi
  double temperature;                     //temperature scaling parameter
  double center[3];                       //origin of the potential
  double KBEh;                            //Boltzmann constant in Hartree
  int potentialtype;                      //type of potential to be used: 0 -> log-Fermi (default)
                                          //                              1 -> power
public:
  //radii must be positive
  EllipsoidalConstraint(double r1 = 1.0, double r2 = -1.0, double r3 = -1.0) {
    radii[0] = r1;
    if (r2 < 0.0) {radii[1] = r1;}
    else {radii[1] = r2;}
    if (r3 < 0.0) {radii[2] = r1;}
    else {radii[2] = r3;}
    maxrad = r1;
    if (r2 > maxrad) {maxrad = r2;}
    if (r3 > maxrad) {maxrad = r3;}
    for (size_t idx = 0; idx < 3; ++idx) {
      anisotropy[idx] = maxrad/radii[idx];
    }
    maxrad *= dist_Angstrom2aum1;
    exponent = 30.0;
    beta = 6.0;
    temperature = 300.0;
    center[0] = 0.0;
    center[1] = 0.0;
    center[2] = 0.0;
    potentialtype = 0;
    KBEh = KB/au2J;
  }
  ~EllipsoidalConstraint() {}
  //getters
  std::vector<double> Radii() {
    std::vector<double> rad(3,radii[0]);
    rad[1] = radii[1];
    rad[2] = radii[2];
    return rad;
  }
  double Exponent() {return exponent;}
  double Beta() {return beta;}
  double Temperature() {return temperature;}
  std::vector<double> Center() {
    std::vector<double> centre(3,center[0]);
    centre[1] = center[1];
    centre[2] = center[2];
    return centre;
  }
  int TypeOfPotential() {return potentialtype;}
  //setters
  void setRadii(double r1, double r2 = -1.0, double r3 = -1.0) {
    radii[0] = r1;
    if (r2 < 0.0) {radii[1] = r1;}
    else {radii[1] = r2;}
    if (r3 < 0.0) {radii[2] = r1;}
    else {radii[2] = r3;}
    maxrad = r1;
    if (r2 > maxrad) {maxrad = r2;}
    if (r3 > maxrad) {maxrad = r3;}
    for (size_t idx = 0; idx < 3; ++idx) {
      anisotropy[idx] = maxrad/radii[idx];
    }
    maxrad *= dist_Angstrom2aum1;
  }
  void setExponent(double newexp) {exponent = newexp;}
  void setBeta(double newb) {beta = newb;}
  void setTemperature(double Temp) {temperature = Temp;}
  void setCenter(double c1, double c2, double c3) {
    center[0] = c1;
    center[1] = c2;
    center[2] = c3;
  }
  void setTypeOfPotential(int ptype) {potentialtype = ptype;}
  //other functions
  void LogFermiPotential(double & result, const matrixE & geometry, matrixE & grad) {
    //the log-Fermi potential
    double rAO[3];
    double RAO;
    double expterm;
    double RAOANIS[3];
    double Fermi;
    size_t Natoms = geometry.rows();
    result = 0.0;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      expterm = 0.0;
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        rAO[idcoord] = anisotropy[idcoord]*(geometry(idAtm + 1,idcoord + 1) - center[idcoord])*dist_Angstrom2aum1;
        RAOANIS[idcoord] = rAO[idcoord]*anisotropy[idcoord];
        expterm += rAO[idcoord]*rAO[idcoord];
      }
      RAO = sqrt(expterm);
      expterm = exp(beta*(RAO - maxrad));
      Fermi = dist_Angstrom2aum1*KBEh*temperature/(1.0 + expterm);
      result += log(1.0 + expterm);
      grad(3*idAtm + 1,1) += beta*expterm*Fermi*RAOANIS[0]/(RAO + 1.0e-14);
      grad(3*idAtm + 2,1) += beta*expterm*Fermi*RAOANIS[1]/(RAO + 1.0e-14);
      grad(3*idAtm + 3,1) += beta*expterm*Fermi*RAOANIS[2]/(RAO + 1.0e-14);
    }
    result *= KBEh*temperature;
  }
  void PowerPotential(double & result, const matrixE & geometry, matrixE & grad) {
    //the power potential function
    double RAO2;
    double polyterm;
    double RAOANIS[3];
    size_t Natoms = geometry.rows();
    result = 0.0;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      RAO2 = 0.0;
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        polyterm = anisotropy[idcoord]*(geometry(idAtm + 1,idcoord + 1) - center[idcoord])*dist_Angstrom2aum1;
        RAOANIS[idcoord] = polyterm*anisotropy[idcoord];
        RAO2 += polyterm*polyterm;
      }
      polyterm = pow(RAO2/(maxrad*maxrad),0.5*exponent);
      result += polyterm;
      grad(3*idAtm + 1,1) += exponent*polyterm*RAOANIS[0]/(RAO2 + 1.0e-14);
      grad(3*idAtm + 2,1) += exponent*polyterm*RAOANIS[1]/(RAO2 + 1.0e-14);
      grad(3*idAtm + 3,1) += exponent*polyterm*RAOANIS[2]/(RAO2 + 1.0e-14);
    }
  }
  double Calculate(const matrixE & geometry, matrixE & grad) {
    //function that evaluates potential according to setup
    double result;
    if (potentialtype == 0) {LogFermiPotential(result,geometry,grad);}
    else if (potentialtype == 1) {PowerPotential(result,geometry,grad);}
    return result;
  }
};

#endif //_Model_Potentials_
