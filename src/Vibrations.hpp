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

#ifndef _Vibrations_
#define _Vibrations_
#include <vector>
#include <string>
#include <math.h>
#include "ConstantsPackage.hpp"
#include "Boltzmann.hpp"

//description:
//models for vibrational partition functions

class HarmonicOscillator: public Boltzmann {
  //quantum mechanical description of harmonic oscillator
  //if vibrational frequency is negative or imaginary, then skip, in accordance to Transition State Theory
  std::vector<double> vibrations;      //vibrational frequencies in cm-1
  std::vector<double> ThetaVibr;       //characteristic vibrational temperatures in K
  size_t tvsize;
public:
  HarmonicOscillator(double _T, std::vector<double> _vibr) {
    T(_T);
    vibrations = _vibr;
    //convert the vibrational energies to Kelvin
    CalcTheta();
  }
  HarmonicOscillator(double _T) {
    T(_T);
  }
  ~HarmonicOscillator() {vibrations.erase(vibrations.begin(),vibrations.end());}
//----------------------------------------------------------------------------------------------------------------
  void setVibrations(std::vector<double> _vibr) {
    vibrations = _vibr;
    CalcTheta();
  }
  std::vector<double> getVibrations() {return vibrations;}
  std::vector<double> getThetas() {return ThetaVibr;}
  void CalcTheta() {
    //function that calculates characteristic vibrational temperatures
    //ThetaVibr.resize(vibrations.size());
    ThetaVibr.clear();
    for (size_t idx = 0; idx < vibrations.size(); ++idx) {
      //ThetaVibr[idx] = 29979245800.0*vibrations[idx]*h/KB;
      if (vibrations[idx] < 0) {
        std::cout << "skipping imaginary frequency " << vibrations[idx] << std::endl;
        ThetaVibr.push_back(0.0);
        continue;
      }
      ThetaVibr.push_back(29979245800.0*vibrations[idx]*h/KB);
    }
    tvsize = ThetaVibr.size();
  }
  double TSfrequency() {
    //function returning the first negative frequency showing up in the frequency vector
    double freq = 0.0;
    for (size_t idx = 0; idx < vibrations.size(); ++idx) {
      if (vibrations[idx] < 0) {
        freq = vibrations[idx];
        break;
      }
    }
    return freq;
  }
//----------------------------------------------------------------------------------------------------------------
  double qHOsingle(double theta, double tm1) {
    //evaluation of the harmonic oscillator's partition function for a single particle
    return 1.0/(1.0 - exp(-theta*tm1));
  }
  double dqHOsingledT(double theta, double tm1) {
    //first derivative with respect to temperature of a single element of the harmonic oscillator partition function
    double expterm = exp(theta*tm1);
    return theta*tm1*tm1*expterm/((expterm - 1)*(expterm - 1));
  }
  double dqHOdT(double tm1, size_t idj) {
    //evaluation of the first derivative with respect to temperature for the harmonic oscillator's partition function
    double otherthetai = 1.0;
    for (size_t idy = 0; idy < tvsize; ++idy) {
      if (idy == idj) {continue;}
      otherthetai *= qHOsingle(ThetaVibr[idy],tm1);
    }
    return otherthetai*dqHOsingledT(ThetaVibr[idj],tm1);
  }
  double d2qHOsingledT2(double theta, double temp, double tm1) {
    //second derivative with respect to temperature of a single element of the harmonic oscillator partition function
    double expterm = exp(theta*tm1);
    return theta*tm1*tm1*expterm/((expterm - 1)*(expterm - 1))*(2.0*temp + theta - (2.0*temp - theta)*expterm)*tm1*tm1/(expterm - 1.0);
  }
  double d2qHOdT2(double temp, double tm1, size_t idj) {
    double othertheta1 = 1.0;
    double othertheta2 = 0.0;
    double othertheta_aux;
    for (size_t idy = 0; idy < tvsize; ++idy) {
      if (idy == idj) {continue;}
      othertheta1 *= qHOsingle(ThetaVibr[idy],tm1);
      othertheta_aux = 1.0;
      for (size_t idz = 0; idz < tvsize; ++idz) {
        if ((idz == idj)||(idz == idy)) {continue;}
        othertheta_aux *= qHOsingle(ThetaVibr[idz],tm1);
      }
      othertheta2 += othertheta_aux*dqHOsingledT(ThetaVibr[idy],tm1);
    }
    othertheta1 *= d2qHOsingledT2(ThetaVibr[idj],temp,tm1);
    othertheta2 *= dqHOsingledT(ThetaVibr[idj],tm1);
    return othertheta1 + othertheta2;
  }
  double qFRsingle(double vibr, double temp, double Bavg) {
    //evaluation of the free rotor's partition function for a single particle
    double beta = au2J/(KB*temp);
    double omega = vibr*cmm12au;
    double mu = 0.5/(omega + 1.0e-14);
    double mup = mu*Bavg/(mu + Bavg);
    return sqrt(2.0*pi*mup/beta);
  }
  double dqFRsingledT(double vibr, double temp, double Bavg) {
    //evaluation of the free rotor's partition function for a single particle
    return 0.5*qFRsingle(vibr,temp,Bavg)/temp;
  }
  double dqFRdT(double temp, double Bavg, size_t idj) {
    //evaluation of the first derivative with respect to temperature for the harmonic oscillator's partition function
    double otherthetai = 1.0;
    for (size_t idy = 0; idy < tvsize; ++idy) {
      if (idy == idj) {continue;}
      otherthetai *= qFRsingle(vibrations[idy],temp,Bavg);
    }
    return otherthetai*dqFRsingledT(vibrations[idj],temp,Bavg);
  }
  double d2qFRsingledT2(double vibr, double temp, double Bavg) {
    //evaluation of the free rotor's partition function for a single particle
    return -0.25*qFRsingle(vibr,temp,Bavg)/(temp*temp);
  }
  double d2qFRdT2(double temp, double Bavg, size_t idj) {
    double othertheta1 = 1.0;
    double othertheta2 = 0.0;
    double othertheta_aux;
    for (size_t idy = 0; idy < tvsize; ++idy) {
      if (idy == idj) {continue;}
      othertheta1 *= qFRsingle(vibrations[idy],temp,Bavg);
      othertheta_aux = 1.0;
      for (size_t idz = 0; idz < tvsize; ++idz) {
        if ((idz == idj)||(idz == idy)) {continue;}
        othertheta_aux *= qFRsingle(vibrations[idz],temp,Bavg);
      }
      othertheta2 += othertheta_aux*dqFRsingledT(vibrations[idy],temp,Bavg);
    }
    othertheta1 *= d2qFRsingledT2(vibrations[idj],temp,Bavg);
    othertheta2 *= dqFRsingledT(vibrations[idj],temp,Bavg);
    return othertheta1 + othertheta2;
  }
  double qMIXsingle(double weight, double theta, double vibr, double temp, double tm1, double Bavg) {
    return pow(qHOsingle(theta,tm1),weight)*pow(qFRsingle(vibr,temp,Bavg),1.0 - weight);
  }
  double dqMIXdTsingle(double weight, double theta, double vibr, double temp, double tm1, double Bavg) {
    double qho = qHOsingle(theta,tm1);
    double qfr = qFRsingle(vibr,temp,Bavg);
    double qhop = dqHOsingledT(theta,tm1);
    double qfrp = dqFRsingledT(vibr,temp,Bavg);
    double qhowm1 = pow(qho,weight - 1.0);
    double qfrmw = pow(qfr,-weight);
    return qhowm1*qfrmw*((1.0 - weight)*qfrp*qho + weight*qfr*qhop);
  }
  double d2qMIXdT2single(double weight, double theta, double vibr, double temp, double tm1, double Bavg) {
    double qho = qHOsingle(theta,tm1);              //Q(HO)
    double qfr = qFRsingle(vibr,temp,Bavg);         //Q(FR)
    double qhop = dqHOsingledT(theta,tm1);          //dQ(HO)/dT
    double qfrp = dqFRsingledT(vibr,temp,Bavg);     //dQ(FR)/dT
    double qhopp = d2qHOsingledT2(theta,temp,tm1);  //d2Q(HO)/dT2
    double qfrpp = d2qFRsingledT2(vibr,temp,Bavg);  //d2Q(FR)/dT2
    //auxiliary variables
    double gamma = qfrp*qho;
    double delta = qfr*qhop;
    double aux = qfr*qho*((1.0 - weight)*qfrpp*qho + weight*qfr*qhopp) + weight*(weight - 1.0)*(gamma - delta)*(gamma - delta);
    return pow(qho,weight - 2.0)*pow(qfr,-weight - 1.0)*aux;
  }
  double dqMIXdT(size_t idj, double weight, double temp, double tm1, double Bavg, double omegazero, double aparam, bool grimmecorr) {
    //evaluation of the first derivative with respect to temperature for the mixed partition function
    double otherthetai = 1.0;
    double omega;
    double weight2;
    double gfactor = double(grimmecorr);
    for (size_t idy = 0; idy < tvsize; ++idy) {
      if (idy == idj) {continue;}
      omega = vibrations[idy]*cmm12au;
      if (omega < 0.0) {continue;}
      weight2 = gfactor*weightGrimme(omega,omegazero,aparam) + 1.0 - gfactor;
      otherthetai *= qMIXsingle(weight2,ThetaVibr[idy],vibrations[idy],temp,tm1,Bavg);
    }
    return otherthetai*dqMIXdTsingle(weight,ThetaVibr[idj],vibrations[idj],temp,tm1,Bavg);
  }
  double d2qMIXdT2(size_t idj, double weight, double temp, double tm1, double Bavg, double omegazero, double aparam, bool grimmecorr) {
    //evaluation of the second derivative with respect to temperature for the mixed partition function
    double otherthetai = 1.0;
    double otherthetaj = 0.0;
    double othertheta_aux;
    double omega;
    double weight2;
    double gfactor = double(grimmecorr);
    for (size_t idy = 0; idy < tvsize; ++idy) {
      if (idy == idj) {continue;}
      othertheta_aux = 1.0;
      for (size_t idz = 0; idz < tvsize; ++idz) {
        if ((idz == idj)||(idz == idy)) {continue;}
        omega = vibrations[idz]*cmm12au;
        if (omega < 0.0) {continue;}
        weight2 = gfactor*weightGrimme(omega,omegazero,aparam) + 1.0 - gfactor;
        othertheta_aux *= qMIXsingle(weight2,ThetaVibr[idz],vibrations[idz],temp,tm1,Bavg);
      }
      omega = vibrations[idy]*cmm12au;
      if (omega < 0.0) {continue;}
      weight2 = gfactor*weightGrimme(omega,omegazero,aparam) + 1.0 - gfactor;
      otherthetai *= qMIXsingle(weight2,ThetaVibr[idy],vibrations[idy],temp,tm1,Bavg);
      otherthetaj += othertheta_aux*dqMIXdTsingle(weight2,ThetaVibr[idy],vibrations[idy],temp,tm1,Bavg);
    }
    otherthetai *= d2qMIXdT2single(weight,ThetaVibr[idj],vibrations[idj],temp,tm1,Bavg);
    otherthetaj *= dqMIXdTsingle(weight,ThetaVibr[idj],vibrations[idj],temp,tm1,Bavg);
    return otherthetai + otherthetaj;
  }
  double qX(bool grimmecorr, double omegazero = 75.0, double aparam = 4.0, double Bavg = 1.0e-44) {
    //partition function
    //weight = 1 -> full harmonic oscillator
    //weight = 0 -> full free rotor
    double qvib = 1.0;
    double temperature = T();
    double temperaturem1 = 1.0/temperature;
    double omega;
    double weight;
    double gfactor = double(grimmecorr);
    Bavg *= kg2melec*dist_Angstrom2aum1*dist_Angstrom2aum1*1.0e20;      //unit conversion
    omegazero *= cmm12au;                                               //unit conversion
    for (size_t idx = 0; idx < tvsize; ++idx) {
      omega = vibrations[idx]*cmm12au;
      if (omega < 0.0) {continue;}
      weight = gfactor*weightGrimme(omega,omegazero,aparam) + 1.0 - gfactor;
      qvib *= qMIXsingle(weight,ThetaVibr[idx],vibrations[idx],temperature,temperaturem1,Bavg);
    }
    return qvib;
  }
  double dqXdV() {return 0.0;}
  double dqXdT(bool grimmecorr, double omegazero = 75.0, double aparam = 4.0, double Bavg = 1.0e-44) {
    //first derivative of the partition function with respect to temperature
    //weight = 1 -> full harmonic oscillator
    //weight = 0 -> full free rotor
    double dqvibdT = 0.0;
    double temperature = T();
    double temperaturem1 = 1.0/temperature;
    double omega;
    double weight;
    double gfactor = double(grimmecorr);
    Bavg *= kg2melec*dist_Angstrom2aum1*dist_Angstrom2aum1*1.0e20;      //unit conversion
    omegazero *= cmm12au;                                               //unit conversion
    for (size_t idx = 0; idx < tvsize; ++idx) {
      omega = vibrations[idx]*cmm12au;
      if (omega < 0.0) {continue;}
      weight = gfactor*weightGrimme(omega,omegazero,aparam) + 1.0 - gfactor;
      dqvibdT += dqMIXdT(idx,weight,temperature,temperaturem1,Bavg,omegazero,aparam,grimmecorr);
    }
    return dqvibdT;
  }
  double d2qXdT2(bool grimmecorr, double omegazero = 75.0, double aparam = 4.0, double Bavg = 1.0e-44) {
    //second derivative of the partition function with respect to temperature
    double d2qvibdT2 = 0.0;
    double temperature = T();
    double temperaturem1 = 1.0/temperature;
    double omega;
    double weight;
    double gfactor = double(grimmecorr);
    Bavg *= kg2melec*dist_Angstrom2aum1*dist_Angstrom2aum1*1.0e20;      //unit conversion
    omegazero *= cmm12au;                                               //unit conversion
    for (size_t idx = 0; idx < tvsize; ++idx) {
      omega = vibrations[idx]*cmm12au;
      if (omega < 0.0) {continue;}
      weight = gfactor*weightGrimme(omega,omegazero,aparam) + 1.0 - gfactor;
      d2qvibdT2 += d2qMIXdT2(idx,weight,temperature,temperaturem1,Bavg,omegazero,aparam,grimmecorr);
    }
    return d2qvibdT2;
  }
  void getDerivatives(int level, std::vector<double> & derivs, bool grimmecorr, double omegazero = 75.0, double aparam = 4.0, double Bavg = 1.0e-44) {
    //function that gets all the derivatives up to a certain degree, to speed up calculation
    //weight = 1 -> full harmonic oscillator
    //weight = 0 -> full free rotor
    derivs.resize(level + 1);
    for (size_t idx = 0; idx < level; ++idx) {
      derivs[idx] = 0.0;
    }
    double qvib = 1.0;
    double dqvibdT = 0.0;
    double d2qvibdT2 = 0.0;
    double temperature = T();
    double temperaturem1 = 1.0/temperature;
    double omega;
    double weight;
    double gfactor = double(grimmecorr);
    double otherthetai = 1.0;
    double otherthetaj = 0.0;
    double othertheta_aux = 0.0;
    Bavg *= kg2melec*dist_Angstrom2aum1*dist_Angstrom2aum1*1.0e20;      //unit conversion
    omegazero *= cmm12au;                                               //unit conversion
    //precalculation of all factors needed
    std::vector<double> Qsingle(tvsize,0.0);
    std::vector<double> dQsingledT(tvsize,0.0);
    std::vector<double> d2QsingledT2(tvsize,0.0);
    for (size_t idx = 0; idx < tvsize; ++idx) {
      omega = vibrations[idx]*cmm12au;
      if (omega < 0.0) {
        Qsingle[idx] = 0.0;
        dQsingledT[idx] = 0.0;
        continue;
      }
      weight = gfactor*weightGrimme(omega,omegazero,aparam) + 1.0 - gfactor;
      Qsingle[idx] = qMIXsingle(weight,ThetaVibr[idx],vibrations[idx],temperature,temperaturem1,Bavg);
      dQsingledT[idx] = dqMIXdTsingle(weight,ThetaVibr[idx],vibrations[idx],temperature,temperaturem1,Bavg);
    }
    if (level == 2) {
      for (size_t idx = 0; idx < tvsize; ++idx) {
        omega = vibrations[idx]*cmm12au;
        if (omega < 0.0) {
          d2QsingledT2[idx] = 0.0;
          continue;
        }
        weight = gfactor*weightGrimme(omega,omegazero,aparam) + 1.0 - gfactor;
        d2QsingledT2[idx] = d2qMIXdT2single(weight,ThetaVibr[idx],vibrations[idx],temperature,temperaturem1,Bavg);
      }
    }
    //calculation
    for (size_t idx = 0; idx < tvsize; ++idx) {
      otherthetai = 1.0;
      otherthetaj = 0.0;
      for (size_t idy = 0; idy < tvsize; ++idy) {
        if (idy == idx) {continue;}
        if (level == 2) {
          othertheta_aux = 1.0;
          for (size_t idz = 0; idz < tvsize; ++idz) {
            if ((idz == idx)||(idz == idy)) {continue;}
            othertheta_aux *= Qsingle[idz];
          }
        }
        otherthetai *= Qsingle[idy];
        otherthetaj += othertheta_aux*dQsingledT[idy];
      }
      qvib *= Qsingle[idx];
      dqvibdT += otherthetai*dQsingledT[idx];
      d2qvibdT2 += otherthetai*d2QsingledT2[idx] + otherthetaj*dQsingledT[idx];
    }
    derivs[0] = qvib;
    derivs[1] = dqvibdT;
    if (level == 2) {derivs[2] = d2qvibdT2;}
  }
  double logqX(bool grimmecorr, double omegazero = 75.0, double aparam = 4.0, double Bavg = 1.0e-44) {
    //logarithm of partition function
    double logqvib = 0.0;
    double temperature = T();
    double temperaturem1 = 1.0/temperature;
    double omega;
    double weight;
    double gfactor = double(grimmecorr);
    Bavg *= kg2melec*dist_Angstrom2aum1*dist_Angstrom2aum1*1.0e20;      //unit conversion
    omegazero *= cmm12au;                                               //unit conversion
    for (size_t idx = 0; idx < tvsize; ++idx) {
      omega = vibrations[idx]*cmm12au;
      if (omega < 0.0) {continue;}
      weight = gfactor*weightGrimme(omega,omegazero,aparam) + 1.0 - gfactor;
      logqvib += (1.0 - weight)*log(qFRsingle(vibrations[idx],temperature,Bavg)) + weight*log(qHOsingle(ThetaVibr[idx],temperaturem1));
    }
    return logqvib;
  }
  double dlogqXdT(bool grimmecorr, double omegazero = 75.0, double aparam = 4.0, double Bavg = 1.0e-44) {
    //first derivative to temperature of logarithm of the partition function
    double dlogqxdt = 0.0;
    double temperature = T();
    double temperaturem1 = 1.0/temperature;
    double omega;
    double weight;
    double gfactor = double(grimmecorr);
    Bavg *= kg2melec*dist_Angstrom2aum1*dist_Angstrom2aum1*1.0e20;      //unit conversion
    omegazero *= cmm12au;                                               //unit conversion
    for (size_t idx = 0; idx < tvsize; ++idx) {
      omega = vibrations[idx]*cmm12au;
      if (omega < 0.0) {continue;}
      weight = gfactor*weightGrimme(omega,omegazero,aparam) + 1.0 - gfactor;
      dlogqxdt += 0.5*(1.0 - weight)*temperaturem1 + weight*dqHOsingledT(ThetaVibr[idx],temperaturem1)/qHOsingle(ThetaVibr[idx],temperaturem1);
    }
    return dlogqxdt;
  }
  double d2logqXdT2(bool grimmecorr, double omegazero = 75.0, double aparam = 4.0, double Bavg = 1.0e-44) {
    //second derivative to temperature of logarithm of the partition function
    double d2logqxdt2 = 0.0;
    double temperature = T();
    double temperaturem1 = 1.0/temperature;
    double temperaturem2 = temperaturem1*temperaturem1;
    double omega;
    double weight;
    double qhom1;
    double qhop;
    double gfactor = double(grimmecorr);
    Bavg *= kg2melec*dist_Angstrom2aum1*dist_Angstrom2aum1*1.0e20;      //unit conversion
    omegazero *= cmm12au;                                               //unit conversion
    for (size_t idx = 0; idx < tvsize; ++idx) {
      omega = vibrations[idx]*cmm12au;
      if (omega < 0.0) {continue;}
      weight = gfactor*weightGrimme(omega,omegazero,aparam) + 1.0 - gfactor;
      qhom1 = 1.0/qHOsingle(ThetaVibr[idx],temperaturem1);
      qhop = dqHOsingledT(ThetaVibr[idx],temperaturem1);
      d2logqxdt2 += 0.5*(weight - 1.0)*temperaturem2 + weight*qhom1*(d2qHOsingledT2(ThetaVibr[idx],temperature,temperaturem1) - qhop*qhop*qhom1);
    }
    return d2logqxdt2;
  }
  //derivatives to volume
  double dlogqXdV() {return 0.0;}
  double dlogqXdlogV() {return 0.0;}
  double weightGrimme(double omega, double omegazero = 100.0, double aparam = 4.0) {
    double freqfactor = pow(omegazero/omega,aparam);
    return 1.0/(1.0 + freqfactor);
  }
  double ZPVE() {
    //calculation of the zero point vibrational energy contribution to the enthalpy in units of Kelvin
    double Ezpve = 0.0;
    for (size_t idx = 0; idx < tvsize; ++idx) {
      Ezpve += ThetaVibr[idx];
    }
    return 0.5*KB*Ezpve;
  }
  double U(bool grimmecorr = true, double Bavg = 1.0e-44, double omegazero = 75.0, double aparam = 4.0) {
    //Grimme correction as described in S. Grimme, Chem. Eur. J., 18(32), 9955, 2012
    //used by default; this is controlled by the boolean variable grimmecorr
    //even though Grimme suggested using this correction only for entropies, we apply for other thermodynamic quantities
    //Bavg is the average molecular moment of inertia in kg/m2
    //omegazero is the reference frequency for the weights in cm-1
    //aparam is the exponent a in the weight function
    double Uvibr = this->T()*this->T()*this->dlogqXdT(grimmecorr,omegazero,aparam,Bavg);
    return KB*Uvibr;
  }
  double S(bool grimmecorr = true, double Bavg = 1.0e-44, double omegazero = 75.0, double aparam = 4.0) {
    //Grimme correction as described in S. Grimme, Chem. Eur. J., 18(32), 9955, 2012
    //used by default; this is controlled by the boolean variable grimmecorr
    //Bavg is the average molecular moment of inertia in kg/m2
    //omegazero is the reference frequency for the weights in cm-1
    //aparam is the exponent a in the weight function
    double Svibr = this->T()*this->dlogqXdT(grimmecorr,omegazero,aparam,Bavg) + this->logqX(grimmecorr,omegazero,aparam,Bavg);
    return KB*Svibr;
  }
  double H(bool grimmecorr = true, double Bavg = 1.0e-44, double omegazero = 75.0, double aparam = 4.0) {
    //Grimme correction as described in S. Grimme, Chem. Eur. J., 18(32), 9955, 2012
    //used by default; this is controlled by the boolean variable grimmecorr
    //even though Grimme suggested using this correction only for entropies, we apply for other thermodynamic quantities
    //Bavg is the average molecular moment of inertia in kg/m2
    //omegazero is the reference frequency for the weights in cm-1
    //aparam is the exponent a in the weight function
    double Hvibr = this->T()*this->T()*this->dlogqXdT(grimmecorr,omegazero,aparam,Bavg);
    return KB*Hvibr;
  }
  double G(bool grimmecorr = true, double Bavg = 1.0e-44, double omegazero = 75.0, double aparam = 4.0) {return this->H(grimmecorr,Bavg,omegazero,aparam) - this->T()*this->S(grimmecorr,Bavg,omegazero,aparam);}
  virtual double A(bool grimmecorr = true, double Bavg = 1.0e-44, double omegazero = 75.0, double aparam = 4.0) {
    //Grimme correction as described in S. Grimme, Chem. Eur. J., 18(32), 9955, 2012
    //used by default; this is controlled by the boolean variable grimmecorr
    //Bavg is the average molecular moment of inertia in kg/m2
    //omegazero is the reference frequency for the weights in cm-1
    //aparam is the exponent a in the weight function
    double Avibr = -this->T()*this->logqX(grimmecorr,omegazero,aparam,Bavg);
    return KB*Avibr;
  }
  double CV(bool grimmecorr = true, double Bavg = 1.0e-44, double omegazero = 75.0, double aparam = 4.0) {
    //Grimme correction as described in S. Grimme, Chem. Eur. J., 18(32), 9955, 2012
    //used by default; this is controlled by the boolean variable grimmecorr
    //even though Grimme suggested using this correction only for entropies, we apply for other thermodynamic quantities
    //Bavg is the average molecular moment of inertia in kg/m2
    //omegazero is the reference frequency for the weights in cm-1
    //aparam is the exponent a in the weight function
    double CVvibr = this->T()*(2.0*this->dlogqXdT(grimmecorr,omegazero,aparam,Bavg) + this->T()*this->d2logqXdT2(grimmecorr,omegazero,aparam,Bavg));
    return KB*CVvibr;
  }
};
double ScalingFactor(std::string method, std::string basis) {
  //function to store vibrational scaling factors
  //strings take only capital letters
  double sfac = 0.0;
  if (method == "AM1") {sfac = 0.9535;}             //https://cccbdb.nist.gov/vibscalejustx.asp
  else if (method == "PM3") {sfac = 0.9761;}        //A. P. Scott, L. Radom, J. Phys. Chem., 100, 16502, 1996
  else if (method == "PM6") {sfac = 1.062;}         //https://cccbdb.nist.gov/vibscalejustx.asp
  else if (method == "BLYP") {
    if (basis == "6-31Gd") {sfac = 0.9945;}         //A. P. Scott, L. Radom, J. Phys. Chem., 100, 16502, 1996
    else if (basis == "6-311Gdfp") {sfac = 0.9986;} //A. P. Scott, L. Radom, J. Phys. Chem., 100, 16502, 1996
  }
  else if (method == "BP86") {
    if (basis == "6-31Gd") {sfac = 0.9914;}         //A. P. Scott, L. Radom, J. Phys. Chem., 100, 16502, 1996
  }
  else if (method == "B3LYP") {
    if (basis == "6-31Gd") {sfac = 0.9614;}         //A. P. Scott, L. Radom, J. Phys. Chem., 100, 16502, 1996
  }
  else if (method == "B3P86") {
    if (basis == "6-31Gd") {sfac = 0.9558;}         //A. P. Scott, L. Radom, J. Phys. Chem., 100, 16502, 1996
  }
  else if (method == "B3PW91") {
    if (basis == "6-31Gd") {sfac = 0.9573;}         //A. P. Scott, L. Radom, J. Phys. Chem., 100, 16502, 1996
  }
  return sfac;
}

#endif //_Vibrations_
