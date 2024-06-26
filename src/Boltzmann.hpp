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

#ifndef _Boltzmann_
#define _Boltzmann_
#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include "ConstantsPackage.hpp"
#include "UnitConversion.hpp"

//description:
//basis class to calculate energy levels using, e.g., quantum mechanics

class Boltzmann {
  //for the full model one needs the whole set of energies and a temperature
  double Temp;                           //temperature
  std::vector<double> energy;            //energy level
  std::vector<double> gE;                //degeneracy of each energy level
  size_t nlevels;                        //total number of levels to be considered
  double eps0;
public:
  Boltzmann() {
    Temp = 1.0;
    nlevels = 1;
    energy.resize(nlevels);
    gE.resize(nlevels);
    eps0 = 0.0;
  }
  ~Boltzmann() {
    energy.erase(energy.begin(),energy.end());
    gE.erase(gE.begin(),gE.end());
  }
  double E0() {return eps0;}
  const double & T() {return Temp;}
  virtual double Volume() {return 0.0;}
  void T(double Tnew) {
    if (Tnew <= 0.0) {throw std::string("ERROR: Boltzmann.hpp: Boltzmann: T(): temperature must be positive");}
    Temp = Tnew;
  }
  const size_t & nLevels() {return nlevels;}
  void nLevels(size_t nlevelsnew) {
    nlevels = nlevelsnew;
    energy.resize(nlevels);
    gE.resize(nlevels);
  }
  const std::vector<double> & Energy() {return energy;}
  const std::vector<double> & GE() {return gE;}
  const double & Energy(size_t pos) {return energy[pos-1];}
  void setEnergy(const std::vector<double> & _en, std::string unit = "SI") {
    energy.clear();
    double convert = 1.0;
    if (nlevels != _en.size()) {nlevels = _en.size();}
    if (unit == "au") {convert = au2J;}
    eps0 = _en[0]*convert;
    for (size_t idx = 0; idx < _en.size(); ++idx) {
      energy.push_back(_en[idx]*convert - eps0);
    }
  }
  void setgE(const std::vector<double> & _ge) {
    if (nlevels != _ge.size()) {nlevels = _ge.size();}
    gE = _ge;
  }
  virtual double epsj() {return 1.0;}
  virtual double Eshift() {return eps0;}
  virtual double qX() {
    //partition function
    double qx = 0.0;
    double beta = 1.0/(KB*Temp);
    for (size_t idx = 0; idx < nlevels; ++idx) {
      qx += gE[idx]*exp(-energy[idx]*beta);
    }
    return qx;
  }
  virtual double dqXdT() {
    //first derivative of the partition function with respect to temperature
    double dqxdt = 0.0;
    double tempm1 = 1.0/Temp;
    double beta = tempm1/KB;
    for (size_t idx = 0; idx < nlevels; ++idx) {
      dqxdt += gE[idx]*exp(-energy[idx]*beta)*energy[idx];
    }
    return tempm1*beta*dqxdt;
  }
  virtual double d2qXdT2() {
    //first derivative of the partition function with respect to temperature
    double d2qxdt2 = 0.0;
    double kBT = KB*Temp;
    double beta = 1.0/(kBT);
    for (size_t idx = 0; idx < nlevels; ++idx) {
      d2qxdt2 += gE[idx]*exp(-energy[idx]*beta)*energy[idx]*(energy[idx] - 2.0*kBT);
    }
    return beta*beta*d2qxdt2/(Temp*Temp);
  }
  virtual double logqX() {
    //logarithm of partition function
    return log(this->qX());
  }
  virtual double dlogqXdT() {
    //first temperature derivative of logarithm of partition function
    double aux0 = 0.0;
    double aux1 = 0.0;
    double tempm1 = 1.0/Temp;
    double beta = tempm1/KB;
    for (size_t idx = 0; idx < nlevels; ++idx) {
      aux0 += gE[idx]*exp(-energy[idx]*beta);
      aux1 += energy[idx]*gE[idx]*exp(-energy[idx]*beta);
    }
    double dlogqxdt = aux1*beta*tempm1/aux0;
    return dlogqxdt;
  }
  virtual double d2logqXdT2() {
    //second temperature derivative of logarithm of partition function
    double E0 = energy[0];
    double der0 = 0.0;
    double der1 = 0.0;
    double der2 = 0.0;
    double kBT = KB*Temp;
    double tempm1 = 1.0/Temp;
    double beta = tempm1/KB;
    for (size_t idx = 0; idx < nlevels; ++idx) {
      der0 += gE[idx]*exp(-energy[idx]*beta);
      der1 += energy[idx]*gE[idx]*exp(-energy[idx]*beta)*beta*tempm1;
      der2 += energy[idx]*(energy[idx] - 2.0*kBT)*gE[idx]*exp(-energy[idx]*beta)*beta*beta*tempm1*tempm1;
    }
    return der2/der0 - der1*der1/(der0*der0);
  }
  virtual double dlogqXdlogT() {return 0.0;}
  virtual double dlogqXdV() {return 0.0;}
  virtual double dlogqXdlogV() {return 0.0;}
  virtual double U() {return KB*this->T()*this->T()*this->dlogqXdT();}
  virtual double H() {return KB*this->T()*(this->T()*this->dlogqXdT() + this->Volume()*this->dlogqXdV());}
  virtual double S() {return KB*(this->T()*this->dlogqXdT() + this->logqX());}
  virtual double G() {return this->H() - this->T()*this->S();}
  virtual double A() {return -KB*this->T()*this->logqX();}
  virtual double P() {return KB*this->T()*this->dlogqXdV();}
  virtual double CV() {return KB*this->T()*(2.0*this->dlogqXdT() + this->T()*this->d2logqXdT2());}
};

#endif //_Boltzmann_
