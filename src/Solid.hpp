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

#ifndef _Solids_
#define _Solids_
#include <string>
#include <vector>
#include "ConstantsPackage.hpp"
#include "math/IntegralPackage.hpp"
#include "rapidjson/json_interface.hpp"

//description:
//methods for the calculation of properties of solids

class Debye {
  //class that calculates thermodynamical properties of the Debye crystal
  double temperature;        //K
  double density;            //g/m3
  double speedsound;         //
  double MW;                 //g/mol
  double phi0;               //J/mol
  double nu_m;
  double x_m;
  double thetaD;
public:
  Debye(double _temp, double _rho = 0.0, double _sos = 0.0, double _mw = 0.0, double _phi0 = 0.0) {
    temperature = _temp;
    density = _rho;
    speedsound = _sos;
    MW = _mw;
    phi0 = _phi0;
    calc_nu_m();
    calc_xm();
    calc_thetaD();
  }
  ~Debye() {}
  double getMW() {return MW;}
  double getPhi0() {return phi0;}
  double getThetaD() {return thetaD;}
  double getDensity() {return density;}
  double getSpeedSound() {return speedsound;}
  void calc_nu_m() {
    //fuction that calculates the threshold frequency nu_m in the determination of radial distribution function
    //this should be called only once for each solid
    double aux = (3.0*NA*density)/(4.0*pi*MW);
    nu_m = pow(aux,1.0/3.0)*speedsound;
  }
  void calc_xm() {
    //fuction that calculates the threshold frequency x_m in the determination of radial distribution function
    //this should be called everytime the temperature changes
    x_m = (h*nu_m)/(KB*temperature);
  }
  void calc_thetaD() {
    //fuction that calculates the Debye characteristic temperature theta_D
    //this should be called only once for each solid
    thetaD = h*nu_m/KB;
  }
  void getDB(std::string compound) {
    //function that gets data from database for calculation
    std::vector<std::string> rstring;
    std::vector<std::string> * prstring = &rstring;
    std::vector<double> rdouble;
    std::vector<double> * prdouble = &rdouble;
    DebyeJSON(compound,prstring,prdouble);
    density = rdouble[2];
    speedsound = rdouble[3];
    MW = rdouble[1];
    phi0 = rdouble[0];
    thetaD = rdouble[4];
  }
  void changeT(double newT) {
    temperature = newT;
    calc_xm();
  }
  void setThetaD(double _newtheta) {
    //function that simply sets thetaD to overcome defining all the solid's parameters
    thetaD = _newtheta;
  }
  double CV() {
    //function calculating the heat capacity at constant volume for a Debye crystal
    double cv = (9.0*R*temperature*temperature*temperature*debye_type1(thetaD/temperature,3))/(thetaD*thetaD*thetaD);
    return cv;
  }
  double S() {
    //function calculating the entropy for a Debye crystal
    double entropy = 0.0;
    double Debye2integral;
    if (fabs(thetaD/temperature) < 1.0) {Debye2integral = debye_type2(thetaD/temperature,3);}
    else {Debye2integral = pi*pi*pi*pi/15.0;}
    entropy += 9.0*R*(temperature/thetaD)*(temperature/thetaD)*(temperature/thetaD)*(Debye2integral - int_0_a_xpmlog1mexpmx(thetaD/temperature,2));
    return entropy;
  }
  double U() {
    //function calculating the internal energy of a Debye crystal
    if (phi0 == 0.0) {std::cout << "WARNING: Solid.hpp: Debye: U(): the calculation of internal energies might be in serious error since phi0 = 0";}
    double intenergy = 0.5*NA*phi0;
    double thetadT = thetaD/temperature;
    double Debye2integral;
    if (fabs(thetadT) < 1.0) {Debye2integral = debye_type2(thetadT,3);}
    else {Debye2integral = pi*pi*pi*pi/15.0;}
    intenergy += 9.0*R*temperature*(Debye2integral + 2.0*(0.5*thetadT)*(0.5*thetadT)*(0.5*thetadT)*(0.5*thetadT))/(thetadT*thetadT*thetadT);
    return intenergy;
  }
  double A() {
    //function calculating the Helmoltz energy of a Debye crystal
    if (phi0 == 0.0) {std::cout << "WARNING: Solid.hpp: Debye: A(): the calculation of Helmoltz energies might be in serious error since phi0 = 0";}
    double helmoltz = 0.5*NA*phi0 + 9.0*R*thetaD/8.0;
    helmoltz += (9.0*R*temperature*temperature*temperature*temperature*int_0_a_xpmlog1mexpmx(thetaD/temperature,2))/(thetaD*thetaD*thetaD);
    return helmoltz;
  }
  double mu() {
    //function that calculates the chemical potential of the Debye solid
    if (phi0 == 0.0) {std::cout << "WARNING: Solid.hpp: Debye: mu(): the calculation of chemical potential might be in serious error since phi0 = 0";}
    double mu = 0.5*phi0 + 9.0*KB*thetaD/8.0;
    mu += (9.0*KB*temperature*temperature*temperature*temperature*int_0_a_xpmlog1mexpmx(thetaD/temperature,2))/(thetaD*thetaD*thetaD);
    return mu;
  }
  double logPvap(std::vector<double> vibrations, double sigma, std::vector<double> inertia) {
    //function that calculates the vapour pressure of the Debye crystal given an ideal gas
    //vibrations is the vector containing all the vibrational frequencies of a single molecule in cm-1
    //sigma is the symmetry factor used in rotations
    //inertia is the vector of eigenvalues of the inertia matrix
    if (phi0 == 0.0) {std::cout << "WARNING: Solid.hpp: Debye: logPvap(): the calculation of vapour pressures might be in serious error since phi0 = 0";}
    double logpvap = this->mu()/(KB*temperature) + 4*log(KB*temperature) + 0.5*log(pi) - log(sigma);
    //get single particle mass
    double mass = 1000.0*MW/NA;
    logpvap += 1.5*log(2.0*pi*mass/(h*h)) + 1.5*log(8*pi*pi/(h*h));
    //get vibrational frequencies in K
    double theta;
    for (size_t idx = 0; idx < vibrations.size(); ++idx) {
      theta = 100.0*c*h*vibrations[idx]/KB;
      logpvap -= (0.5*theta/temperature + log(1.0 - exp(-theta/temperature)));
    }
    logpvap += 0.5*log(inertia[0]*inertia[1]*inertia[2]);
    return logpvap;
  }
  double Hsub(std::vector<double> vibrations, size_t Tdependence = 0) {
    //function that calculates the heat of sublimation for the Debye crystal given an ideal gas
    //vibrations is the vector containing all the vibrational frequencies of a single molecule in cm-1
    //only if Tdependence == 0 we have the most simple expression. Otherwise use temperature dependence.
    if (phi0 == 0.0) {std::cout << "WARNING: Solid.hpp: Debye: Hsub(): the calculation of heats of sublimation might be in serious error since phi0 = 0";}
    double hsub = -9.0*R*thetaD/8.0 - 0.5*NA*phi0;
    if (Tdependence != 0) {
      //addition of temperature dependent terms
      hsub += 4.0*R*temperature;
      hsub -= (9.0*R*temperature*temperature*temperature*temperature*debye_type2(thetaD/temperature,3))/(thetaD*thetaD*thetaD);
    }
    //get vibrational frequencies in K
    double theta;
    for (size_t idx = 0; idx < vibrations.size(); ++idx) {
      theta = 100.0*c*h*vibrations[idx]/KB;
      hsub += 0.5*theta*R;
      if (Tdependence != 0) {hsub += R*theta/(exp(theta/temperature) - 1.0);}
    }
    return hsub;
  }
};

#endif //_Solids_
