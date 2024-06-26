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

#ifndef _Rotations_
#define _Rotations_
#include <vector>
#include <string>
#include <math.h>
#include "ConstantsPackage.hpp"
#include "Boltzmann.hpp"

//description:
//methods for the statistial mechanical treatment of rotations

class RigidRotor: public Boltzmann {
  //quantum mechanical description of rigid rotor
  std::vector<double> Inertia;         //eigenvalues of inertia matrix. They are supposed to be ordered from lower to greater.
  int symmgroup;                       //0 atom; 1 diatomic; 2 spherical top; 3 prolate top; 4 oblate top; 5 asymmetric top
  unsigned int j;                      //main quantum number
public:
  RigidRotor(double _T) {
    T(_T);
    Inertia.resize(3);
    for (size_t ind = 0; ind < 3; ++ind) {
      Inertia[ind] = 1.0;
    }
    SymmGroup();
    j = 1;
  }
  RigidRotor(double _T, std::vector<double> _Inertia) {
    T(_T);
    Inertia.resize(3);
    Inertia = _Inertia;
    SymmGroup();
    j = 1;
  }
  ~RigidRotor() {Inertia.erase(Inertia.begin(),Inertia.end());}
//----------------------------------------------------------------------------------------------------------------
  void SymmGroup(int _sym) {symmgroup = _sym;}
  void SymmGroup(std::string _sym) {
    if (_sym == "atom") {symmgroup = 0;}
    else if (_sym == "linear") {symmgroup = 1;}
    else if (_sym == "spherical top") {symmgroup = 2;}
    else if (_sym == "prolate top") {symmgroup = 3;}
    else if (_sym == "oblate top") {symmgroup = 4;}
    else if (_sym == "asymmetric") {symmgroup = 5;}
  }
  unsigned int getSymmGroup() {return symmgroup;}
  void SymmGroup() {
    //function to determine the symmetry group. Atom and diatomic do not come from eigenvalue calculations!
    //inertia = (-1,-1,-1) -> atom
    //inertia = (I,-1,-1)  -> diatomic
    if ((Inertia[1] == 0)&&(Inertia[2] == 0)) {
      if (Inertia[0] == 0) {symmgroup = 0;}            //atom
      else {symmgroup = 1;}                            //diatomic
    }
    //polyatomic
    else {
      if (Inertia[0] == Inertia[1]) {
        if (Inertia[1] == Inertia[2]) {symmgroup = 2;} //spherical top
        else {symmgroup = 4;}                          //oblate top
      }
      else {
        if (Inertia[1] == Inertia[2]) {symmgroup = 3;} //prolate top
        else {symmgroup = 5;}                          //asymmetric top
      }
    }
  }
  void setInertia(std::vector<double> _inertia) {
    if (_inertia.size() != 3) {throw std::string("ERROR: Rotations.hpp: RigidRotor: setInertia(): the inertia matrix has only 3 eigenvalues");}
    Inertia = _inertia;
    SymmGroup();
  }
  std::vector<double> getInertia() {return Inertia;}
//----------------------------------------------------------------------------------------------------------------
  unsigned int getj() {return j;}
  void ppj() {++j;}
  void mmj() {--j;}
  void setj(unsigned int _j) {j = _j;}
//----------------------------------------------------------------------------------------------------------------
  double A(unsigned int k = 0) {
    double _A = 0.0;           //atom case
    if ((symmgroup == 1)||((symmgroup == 2))) {
      _A = (hbar*hbar*double(j)*double(j+1))/(2.0*KB*Inertia[0]);
    }
    else if ((symmgroup == 3)||((symmgroup == 4))) {
      _A = (hbar*hbar*double(j)*double(j+1))/(2.0*KB*Inertia[1]);
      double _Ix = Inertia[0];
      if (symmgroup == 4) {_Ix = Inertia[3];}
      _A += (hbar*hbar*double(k)*double(k)*(Inertia[1]-_Ix))/(2.0*KB*Inertia[1]*_Ix);
    }
    else if (symmgroup == 5) {
      throw std::string("WARNING: Rotations.hpp: RigidRotor: A(): asymmetric tops have no analytical solution for the rotational energy");
    }
    return _A;
  }
//----------------------------------------------------------------------------------------------------------------
  double qX(double EnMax = 1.5*KB*300) {
    //saving original state of j
    unsigned int savej = j;
    //zeroing j
    j = 1;
    double _qX = 0.0;
    while (KB*A() <= EnMax) {
      _qX += double(2*j+1)*exp((-A())/(T()));
      if ((symmgroup == 3)||((symmgroup == 4))) {
        for (unsigned int k = 1; k < j; ++k) {
          _qX += 2.0*double(2*j+1)*exp((-A(k))/(T()));
          if (KB*A(k) > EnMax) {break;}
        }
      }
      ppj();
    }
    j = savej;
    return _qX;
  }
  virtual double logqX(double EnMax = 1.5*KB*300) {return (log(qX(EnMax)));}
  double EnExpVal(double EnMax = 1.5*KB*300) {
    //saving original state of j
    unsigned int savej = j;
    //zeroing j
    j = 1;
    double _expval = 0.0;
    while (KB*A() <= EnMax) {
      _expval += double(2*j+1)*A()*exp((-A())/(T()));
      if ((symmgroup == 3)||((symmgroup == 4))) {
        for (unsigned int k = 1; k < j; ++k) {
          _expval += 2.0*double(2*j+1)*A(k)*exp((-A(k))/(T()));
          if (KB*A(k) > EnMax) {break;}
        }
      }
      ppj();
    }
    j = savej;
    return (_expval/qX(EnMax));
  }
  double EnExpVal2(double EnMax = 1.5*KB*300) {
    //saving original state of j
    unsigned int savej = j;
    //zeroing j
    j = 1;
    double _expval2 = 0.0;
    while (KB*A() <= EnMax) {
      _expval2 += double(2*j+1)*A()*A()*exp((-A())/(T()));
      if ((symmgroup == 3)||((symmgroup == 4))) {
        for (unsigned int k = 1; k < j; ++k) {
          _expval2 += 2.0*double(2*j+1)*A(k)*A(k)*exp((-A(k))/(T()));
          if (KB*A(k) > EnMax) {break;}
        }
      }
      ppj();
    }
    j = savej;
    return (_expval2/qX(EnMax));
  }
  virtual double dlogqXdT(double EnMax = 1.5*KB*300) {return (EnExpVal(EnMax)/(T()*T()));}
  virtual double d2logqXdT2(double EnMax = 1.5*KB*300) {
    double der2 = -2.0*EnExpVal(EnMax)/(T()*T()*T());
    der2 += EnExpVal2(EnMax)/(T()*T()*T()*T());
    der2 -= EnExpVal(EnMax)*EnExpVal(EnMax)/(T()*T()*T()*T());
    return der2;
  }
  virtual double dlogqXdlogT(double EnMax = 1.5*KB*300) {return (EnExpVal(EnMax)/(T()));}
  virtual double dlogqXdV(double EnMax = 1.5*KB*300) {return 0.0;}
  virtual double dlogqXdlogV(double EnMax = 1.5*KB*300) {return 0.0;}
};
class RigidRotorlim: public Boltzmann {
  //quantum mechanical description of rigid rotor in the classical limit
  unsigned int sigma;
  std::vector<double> Inertia;         //eigenvalues of inertia matrix. They are supposed to be ordered from lower to greater.
  unsigned int symmgroup;              //0 atom; 1 diatomic; 2 spherical top; 3 prolate top; 4 oblate top; 5 asymmetric top
public:
  RigidRotorlim(double _T, unsigned int _sigma = 1) {
    this->T(_T);
    sigma = _sigma;
    Inertia.resize(3);
    for (size_t ind = 0; ind < 3; ++ind) {
      Inertia[ind] = 1.0;
    }
    this->SymmGroup();
  }
  RigidRotorlim(double _T, std::vector<double> _Inertia, unsigned int _sigma = 1) {
    this->T(_T);
    sigma = _sigma;
    Inertia.resize(3);
    Inertia = _Inertia;
    this->SymmGroup();
  }
  ~RigidRotorlim() {Inertia.erase(Inertia.begin(),Inertia.end());}
//----------------------------------------------------------------------------------------------------------------
  void SymmGroup(unsigned int _sym) {symmgroup = _sym;}
  void SymmGroup(std::string _sym) {
    if (_sym == "atom") {symmgroup = 0;}
    else if (_sym == "linear") {symmgroup = 1;}
    else if (_sym == "spherical top") {symmgroup = 2;}
    else if (_sym == "prolate top") {symmgroup = 3;}
    else if (_sym == "oblate top") {symmgroup = 4;}
    else if (_sym == "asymmetric") {symmgroup = 5;}
  }
  void SymmGroup() {
    //function to determine the symmetry group. Atom and diatomic do not come from eigenvalue calculations!
    //inertia = (-1,-1,-1) -> atom
    //inertia = (I,-1,-1)  -> diatomic
    if ((Inertia[1] == 0)&&(Inertia[2] == 0)) {
      if (Inertia[0] == 0) {symmgroup = 0;}            //atom
      else {symmgroup = 1;}                            //diatomic
    }
    //polyatomic
    else {
      if (Inertia[0] == Inertia[1]) {
        if (Inertia[1] == Inertia[2]) {symmgroup = 2;} //spherical top
        else {symmgroup = 4;}                          //oblate top
      }
      else {
        if (Inertia[1] == Inertia[2]) {symmgroup = 3;} //prolate top
        else {symmgroup = 5;}                          //asymmetric top
      }
    }
  }
  unsigned int getSymmGroup() {return symmgroup;}
  void setInertia(std::vector<double> _inertia) {
    if (_inertia.size() != 3) {throw std::string("ERROR: Rotations.hpp: RigidRotorlim: setInertia(): the inertia matrix has only 3 eigenvalues");}
    Inertia = _inertia;
    this->SymmGroup();
  }
  std::vector<double> getInertia() {return Inertia;}
  unsigned int getSigma() {return sigma;}
  void setSigma(unsigned int _sigma) {sigma = _sigma;}
//----------------------------------------------------------------------------------------------------------------
  double qX() {
    double _qX = 1.0;                                //atom
    if ((symmgroup == 1)||(symmgroup == 2)) {        //diatomic and spherical top
      double B = (hbar)/(4*pi*c*Inertia[0]);
      _qX = (KB*this->T())/(sigma*h*c*B);
    }
    else if ((symmgroup == 3)||(symmgroup == 4)) {   //oblate and prolate tops
      double cube = sqrt((2.0*KB*this->T()*Inertia[1])/(hbar*hbar));
      double Ix = Inertia[0]*(symmgroup == 3) + Inertia[2]*(symmgroup == 4);
      _qX = cube*cube*cube*sqrt(pi*Ix/fabs(Inertia[1] - Ix));
    }
    else if (symmgroup == 5) {                       //asymmetric top
      double A = (hbar)/(4*pi*c*Inertia[0]);
      double B = (hbar)/(4*pi*c*Inertia[1]);
      double C = (hbar)/(4*pi*c*Inertia[2]);
      double rotA = sqrt((KB*this->T())/(h*c*A));
      double rotB = sqrt((KB*this->T())/(h*c*B));
      double rotC = sqrt((KB*this->T())/(h*c*C));
      _qX = (sqrt(pi)/sigma)*rotA*rotB*rotC;
    }
    return _qX;
  }
  double dqXdV() {return 0.0;}
  double dqXdT() {
    double _dqXdT = 0.0;                             //atom
    if ((symmgroup == 1)||(symmgroup == 2)) {        //diatomic and spherical top
      double B = (hbar)/(4*pi*c*Inertia[0]);
      _dqXdT = (KB)/(sigma*h*c*B);
    }
    else if ((symmgroup == 3)||(symmgroup == 4)) {   //oblate and prolate tops
      double cube = sqrt((2.0*KB*this->T()*Inertia[1])/(hbar*hbar));
      double Ix = Inertia[0]*(symmgroup == 3) + Inertia[2]*(symmgroup == 4);
      _dqXdT = 1.5*cube*cube*cube*sqrt(pi*Ix/fabs(Inertia[1] - Ix))/this->T();
    }
    else if (symmgroup == 5) {                       //asymmetric top
      double A = (hbar)/(4*pi*c*Inertia[0]);
      double B = (hbar)/(4*pi*c*Inertia[1]);
      double C = (hbar)/(4*pi*c*Inertia[2]);
      double rotA = sqrt((KB*this->T())/(h*c*A));
      double rotB = sqrt((KB*this->T())/(h*c*B));
      double rotC = sqrt((KB*this->T())/(h*c*C));
      _dqXdT = 1.5*(sqrt(pi)/sigma)*rotA*rotB*rotC/this->T();
    }
    return _dqXdT;
  }
  double d2qXdT2() {
    double _d2qXdT2 = 0.0;                           //atom
    if ((symmgroup == 1)||(symmgroup == 2)) {        //diatomic and spherical top
      _d2qXdT2 = 0.0;
    }
    else if ((symmgroup == 3)||(symmgroup == 4)) {   //oblate and prolate tops
      double cube = sqrt((2.0*KB*this->T()*Inertia[1])/(hbar*hbar));
      double Ix = Inertia[0]*(symmgroup == 3) + Inertia[2]*(symmgroup == 4);
      _d2qXdT2 = 0.75*cube*cube*cube*sqrt(pi*Ix/fabs(Inertia[1] - Ix))/(this->T()*this->T());
    }
    else if (symmgroup == 5) {                       //asymmetric top
      double A = (hbar)/(4*pi*c*Inertia[0]);
      double B = (hbar)/(4*pi*c*Inertia[1]);
      double C = (hbar)/(4*pi*c*Inertia[2]);
      double rotA = sqrt((KB*this->T())/(h*c*A));
      double rotB = sqrt((KB*this->T())/(h*c*B));
      double rotC = sqrt((KB*this->T())/(h*c*C));
      _d2qXdT2 = 0.75*(sqrt(pi)/sigma)*rotA*rotB*rotC/(this->T()*this->T());
    }
    return _d2qXdT2;
  }
  double logqX() {return (log(this->qX()));}
  double dlogqXdT() {
    double _dlogqXdT = 0;                                            //atom
    if ((symmgroup == 1)||(symmgroup == 2)) {                        //diatomic and spherical top
      _dlogqXdT = 1.0/this->T();
    }
    else if ((symmgroup == 3)||(symmgroup == 4)||(symmgroup == 5)) { //oblate, prolate and asymmetric tops
      _dlogqXdT = 3.0/(2.0*this->T());
    }
    return _dlogqXdT;
  }
  double d2logqXdT2() {
    double d2logqXdT2 = 0;                                            //atom
    if ((symmgroup == 1)||(symmgroup == 2)) {                         //diatomic and spherical top
      d2logqXdT2 = -1.0/(this->T()*this->T());
    }
    else if ((symmgroup == 3)||(symmgroup == 4)||(symmgroup == 5)) {  //oblate, prolate and asymmetric tops
      d2logqXdT2 = -3.0/(2.0*this->T()*this->T());
    }

    return d2logqXdT2;
  }
  double dlogqXdlogT() {
    double dlogqXdlogT = 0;                                           //atom
    if ((symmgroup == 1)||(symmgroup == 2)) {                         //diatomic and spherical top
      dlogqXdlogT = 1.0;
    }
    else if ((symmgroup == 3)||(symmgroup == 4)||(symmgroup == 5)) {  //oblate, prolate and asymmetric tops
      dlogqXdlogT = 1.5;
    }

    return dlogqXdlogT;
  }
  double dlogqXdV() {return 0.0;}
  double dlogqXdlogV() {return 0.0;}
};

#endif //_Rotations_
