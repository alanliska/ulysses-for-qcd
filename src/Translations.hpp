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

#ifndef _Translations_
#define _Translations_
#include <vector>
#include <string>
#include <iostream>
#include "ConstantsPackage.hpp"
#include "Boltzmann.hpp"

//description:
//models for translational motion of molecules

class PartBox: public Boltzmann {
  //quantum mechanical description of particle in a box model, ndim dimensional
  double mass;
  size_t ndim;
  std::vector<double> L;
  std::vector<double> k;
  std::vector<unsigned int> ni;
public:
  PartBox(double _T, double _mass, size_t _nlevels, size_t _ndim = 1) {
    T(_T);
    nLevels(_nlevels);
    ndim = _ndim;
    mass = _mass;
    std::vector<double> en_(nLevels(),0.0);
    setEnergy(en_);
    L.resize(ndim);
    ni.resize(ndim);
    for (size_t ind = 0; ind < ndim; ++ind) {
      L[ind] = 1.0;
      ni[ind] = 1;
    }
    k.resize(ndim);
    calck();
  }
  ~PartBox() {
    L.erase(L.begin(),L.end());
    ni.erase(ni.begin(),ni.end());
    k.erase(k.begin(),k.end());
  }
//----------------------------------------------------------------------------------------------------------------
  size_t nDimensions() {return ndim;}
  void nDimensions(size_t _ndim) {
    ndim = _ndim;
    L.resize(ndim);
    ni.resize(ndim);
    for (size_t ind = 0; ind < ndim; ++ind) {
      L[ind] = 1.0;
      ni[ind] = 1;
    }
    k.resize(ndim);
    calck();
  }
  void setL(std::vector<double> dimbox) {
    if (dimbox.size() == ndim) {
      L = dimbox;
      calck();
    }
    else {throw std::string("ERROR: Translations.hpp: PartBox: setL(): the vector introduced as argument does not have the correct dimension");}
  }
  const std::vector<double> & Lprint() {return L;}
  const std::vector<double> & kprint() {return k;}
  const std::vector<unsigned int> & niprint() {return ni;}
  void setni(std::vector<unsigned int> levels) {
    if (levels.size() == ndim) {
      ni = levels;
      calck();
    }
    else {throw std::string("ERROR: Translations.hpp: PartBox: setni(): the vector introduced as argument does not have the correct dimension");}
  }
  void ppni(size_t dim) {
    ++ni[dim-1];
    calck();
  }
  void mmni(size_t dim) {
    --ni[dim-1];
    calck();
  }
  void changeni(size_t dim, size_t chg, char sign) {
    if (sign == '+') {
      ni[dim-1] += chg;
      calck();
    }
    else if (sign == '-') {
      if (ni[dim-1] > chg) {
        ni[dim-1] -= chg;
        calck();
      }
      else {throw std::string("ERROR: Translations.hpp: PartBox: changeni(): change too large and the quantum number will zero or become negative");}
    }
    else {throw std::string("ERROR: Translations.hpp: PartBox: changeni(): no sign identified");}
  }
  void calck() {
    for (size_t ind = 0; ind < ndim; ++ind) {
      if (L[ind] != 0.0) {k[ind] = double(ni[ind])*pi/L[ind];}
      else {k[ind] = 0.0;}                                     //If L[ind] == 0.0 then consider the dimension inexistent
    }
  }
  double ksquare() {
    double _ksquare = 0.0;
    for (size_t ind = 0; ind < ndim; ++ind) {
      _ksquare += k[ind]*k[ind];
    }
    return _ksquare;
  }
//----------------------------------------------------------------------------------------------------------------
  double Volume() {
    double vol = 1.0;
    for (size_t ind = 0; ind < Lprint().size(); ++ind) {
      if (Lprint()[ind] != 0.0) {vol *= Lprint()[ind];}
    }
    return vol;
  }
  double Mass() {return mass;}
//----------------------------------------------------------------------------------------------------------------
    virtual double epsj() {return ((hbar*hbar*ksquare())/(2.0*mass));}
//----------------------------------------------------------------------------------------------------------------
  double A() {
    double naux = 0.0;
    double Laux = 1.0;
    for (size_t ind = 0; ind < nDimensions(); ++ind) {
      Laux = 1.0;
      for (size_t ind2 = 0; ind2 < nDimensions(); ++ind2) {
        if (ind2 != ind) {Laux *= L[ind2];}
      }
      naux += ni[ind]*ni[ind]*Laux*Laux;
    }
    return (hbar*hbar*pi*pi*naux)/(2.0*mass*KB);
  }
  double qX(double EnMax = 1.5*KB*300) {
    //saving original state of ni
    std::vector<unsigned int> saveni = ni;
    //zeroing ni
    for (size_t ind = 0; ind < nDimensions(); ++ind) {
      ni[ind] = 1;
    }
    double _qX = 0.0;
    //implementation only for 1, 2, and 3 dimensions. Anything else is not reasonable
    if (nDimensions() == 1) {
      size_t nmax = size_t((L[0]*sqrt(2.0*mass*EnMax))/(hbar*pi));
      //std::cout << "nmax = " << nmax << std::endl;
      for (size_t ind = 0; ind < nmax; ++ind) {
        _qX += exp((-A())/(T()*Volume()*Volume()));
        ppni(1);
      }
    }
    else if (nDimensions() == 2) {
      size_t nxmax = size_t(sqrt(((2.0*mass*EnMax*L[0]*L[0])/(hbar*pi*hbar*pi))-((L[0]*L[0])/(L[1]*L[1]))));
      size_t nymax = size_t(sqrt(((2.0*mass*EnMax*L[1]*L[1])/(hbar*hbar*pi*pi))-((L[1]*L[1])/(L[0]*L[0]))));
      for (size_t indy = 0; indy < nymax; ++indy) {
        ni[0] = 1;
        calck();
        for (size_t indx = 0; indx < nxmax; ++indx) {
          if (epsj() > EnMax) {break;}
          _qX += exp((-A())/(T()*Volume()*Volume()));
          ppni(1);
        }
        ppni(2);
      }
    }
    else if (nDimensions() == 3) {
      size_t nxmax = size_t(sqrt(((2.0*mass*EnMax*L[0]*L[0])/(hbar*hbar*pi*pi))-((L[0]*L[0])/(L[1]*L[1]))-((L[0]*L[0])/(L[2]*L[2]))));
      size_t nymax = size_t(sqrt(((2.0*mass*EnMax*L[1]*L[1])/(hbar*hbar*pi*pi))-((L[1]*L[1])/(L[0]*L[0]))-((L[1]*L[1])/(L[2]*L[2]))));
      size_t nzmax = size_t(sqrt(((2.0*mass*EnMax*L[2]*L[2])/(hbar*hbar*pi*pi))-((L[2]*L[2])/(L[0]*L[0]))-((L[2]*L[2])/(L[1]*L[1]))));
      for (size_t indz = 0; indz < nzmax; ++indz) {
        ni[1] = 1;
        calck();
        for (size_t indy = 0; indy < nymax; ++indy) {
          ni[0] = 1;
          calck();
          for (size_t indx = 0; indx < nxmax; ++indx) {
            if (epsj() > EnMax) {break;}
            _qX += exp((-A())/(T()*Volume()*Volume()));
            ppni(1);
          }
          ppni(2);
        }
        ppni(3);
      }
    }
    else {throw std::string("ERROR: Translations.hpp: PartBox: qX(): number of dimensions not implemented");}
    setni(saveni);
    return _qX;
  }
  double EnExpVal(double EnMax = 1.5*KB*300) {
    //saving original state of ni
    std::vector<unsigned int> saveni = ni;
    //zeroing ni
    for (size_t ind = 0; ind < nDimensions(); ++ind) {
      ni[ind] = 1;
    }
    double _expval = 0.0;
    //implementation only for 1, 2, and 3 dimensions. Anything else is not reasonable
    if (nDimensions() == 1) {
      size_t nmax = size_t((L[0]*sqrt(2.0*mass*EnMax))/(hbar*pi));
      for (size_t ind = 0; ind < nmax; ++ind) {
        _expval += A()*exp((-A())/(T()*Volume()*Volume()));
        ppni(1);
      }
    }
    else if (nDimensions() == 2) {
      size_t nxmax = size_t(sqrt(((2.0*mass*EnMax*L[0]*L[0])/(hbar*pi*hbar*pi))-((L[0]*L[0])/(L[1]*L[1]))));
      size_t nymax = size_t(sqrt(((2.0*mass*EnMax*L[1]*L[1])/(hbar*hbar*pi*pi))-((L[1]*L[1])/(L[0]*L[0]))));
      for (size_t indy = 0; indy < nymax; ++indy) {
        ni[0] = 1;
        calck();
        for (size_t indx = 0; indx < nxmax; ++indx) {
          if (epsj() > EnMax) {break;}
          _expval += A()*exp((-A())/(T()*Volume()*Volume()));
          ppni(1);
        }
        ppni(2);
      }
    }
    else if (nDimensions() == 3) {
      size_t nxmax = size_t(sqrt(((2.0*mass*EnMax*L[0]*L[0])/(hbar*hbar*pi*pi))-((L[0]*L[0])/(L[1]*L[1]))-((L[0]*L[0])/(L[2]*L[2]))));
      size_t nymax = size_t(sqrt(((2.0*mass*EnMax*L[1]*L[1])/(hbar*hbar*pi*pi))-((L[1]*L[1])/(L[0]*L[0]))-((L[1]*L[1])/(L[2]*L[2]))));
      size_t nzmax = size_t(sqrt(((2.0*mass*EnMax*L[2]*L[2])/(hbar*hbar*pi*pi))-((L[2]*L[2])/(L[0]*L[0]))-((L[2]*L[2])/(L[1]*L[1]))));
      for (size_t indz = 0; indz < nzmax; ++indz) {
        ni[1] = 1;
        calck();
        for (size_t indy = 0; indy < nymax; ++indy) {
          ni[0] = 1;
          calck();
          for (size_t indx = 0; indx < nxmax; ++indx) {
            if (epsj() > EnMax) {break;}
            _expval += A()*exp((-A())/(T()*Volume()*Volume()));
            ppni(1);
          }
          ppni(2);
        }
        ppni(3);
      }
    }
    else {throw std::string("ERROR: Translations.hpp: PartBox: EnExpVal(): number of dimensions not implemented.");}
    setni(saveni);
    return (_expval/qX(EnMax));
  }
  double EnExpVal2(double EnMax = 1.5*KB*300) {
    //saving original state of ni
    std::vector<unsigned int> saveni = ni;
    //zeroing ni
    for (size_t ind = 0; ind < nDimensions(); ++ind) {
      ni[ind] = 1;
    }
    double _expval2 = 0.0;
    //implementation only for 1, 2, and 3 dimensions. Anything else is not reasonable
    if (nDimensions() == 1) {
      size_t nmax = size_t((L[0]*sqrt(2.0*mass*EnMax))/(hbar*pi));
      for (size_t ind = 0; ind < nmax; ++ind) {
        _expval2 += A()*A()*exp((-A())/(T()*Volume()*Volume()));
        ppni(1);
      }
    }
    else if (nDimensions() == 2) {
      size_t nxmax = size_t(sqrt(((2.0*mass*EnMax*L[0]*L[0])/(hbar*pi*hbar*pi))-((L[0]*L[0])/(L[1]*L[1]))));
      size_t nymax = size_t(sqrt(((2.0*mass*EnMax*L[1]*L[1])/(hbar*hbar*pi*pi))-((L[1]*L[1])/(L[0]*L[0]))));
      for (size_t indy = 0; indy < nymax; ++indy) {
        ni[0] = 1;
        calck();
        for (size_t indx = 0; indx < nxmax; ++indx) {
          if (epsj() > EnMax) {break;}
          _expval2 += A()*A()*exp((-A())/(T()*Volume()*Volume()));
          ppni(1);
        }
        ppni(2);
      }
    }
    else if (nDimensions() == 3) {
      size_t nxmax = size_t(sqrt(((2.0*mass*EnMax*L[0]*L[0])/(hbar*hbar*pi*pi))-((L[0]*L[0])/(L[1]*L[1]))-((L[0]*L[0])/(L[2]*L[2]))));
      size_t nymax = size_t(sqrt(((2.0*mass*EnMax*L[1]*L[1])/(hbar*hbar*pi*pi))-((L[1]*L[1])/(L[0]*L[0]))-((L[1]*L[1])/(L[2]*L[2]))));
      size_t nzmax = size_t(sqrt(((2.0*mass*EnMax*L[2]*L[2])/(hbar*hbar*pi*pi))-((L[2]*L[2])/(L[0]*L[0]))-((L[2]*L[2])/(L[1]*L[1]))));
      for (size_t indz = 0; indz < nzmax; ++indz) {
        ni[1] = 1;
        calck();
        for (size_t indy = 0; indy < nymax; ++indy) {
          ni[0] = 1;
          calck();
          for (size_t indx = 0; indx < nxmax; ++indx) {
            if (epsj() > EnMax) {break;}
            _expval2 += A()*A()*exp((-A())/(T()*Volume()*Volume()));
            ppni(1);
          }
          ppni(2);
        }
        ppni(3);
      }
    }
    else {throw std::string("ERROR: Translations.hpp: PartBox: EnExpVal2(): number of dimensions not implemented.");}
    setni(saveni);
    return (_expval2/qX(EnMax));
  }
  virtual double logqX(double EnMax = 1.5*KB*300) {return log(qX(EnMax));}
  virtual double dlogqXdT(double EnMax = 1.5*KB*300) {return (EnExpVal(EnMax)/(T()*T()*Volume()*Volume()));}
  virtual double d2logqXdT2(double EnMax = 1.5*KB*300) {
    double der2 = -2.0*EnExpVal(EnMax)/(T()*T()*T()*Volume()*Volume());
    der2 += EnExpVal2(EnMax)/(T()*T()*T()*T()*Volume()*Volume()*Volume()*Volume());
    der2 -= EnExpVal(EnMax)*EnExpVal(EnMax)/(T()*T()*T()*T()*Volume()*Volume()*Volume()*Volume());
    return der2;
  }
  virtual double dlogqXdlogT(double EnMax = 1.5*KB*300) {return (EnExpVal(EnMax)/(T()*Volume()*Volume()));}
  virtual double dlogqXdV(double EnMax = 1.5*KB*300) {return (2.0*EnExpVal(EnMax)/(T()*Volume()*Volume()*Volume()));}
  virtual double dlogqXdlogV(double EnMax = 1.5*KB*300) {return (2.0*EnExpVal(EnMax)/(T()*Volume()*Volume()));}
//----------------------------------------------------------------------------------------------------------------
  std::vector<double> Psiv(std::vector<double> pos) {
    if (pos.size() != ndim) {throw std::string("ERROR: Translations.hpp: PartBox: Psiv(): the vector introduced as argument does not have the correct dimension");}
    bool goodpos = true;
    for (size_t ind = 0; ind < ndim; ++ind) {
      if ((pos[ind] < 0.0)||(pos[ind] > L[ind])) {
        goodpos = false;
        break;
      }
    }
    std::vector<double> _psiv(ndim,0.0);
    if (goodpos) {
      for (size_t ind = 0; ind < ndim; ++ind) {
        if (L[ind] != 0.0) {_psiv[ind] = sqrt(2.0/L[ind])*sin(double(ni[ind])*pi*pos[ind]/L[ind]);}
        else {_psiv[ind] = 0.0;}                                 //If L[ind] == 0.0 then consider the dimension inexistent
      }
    }
    return _psiv;
  }
  double Psi(std::vector<double> pos) {
    std::vector<double> psivec = Psiv(pos);
    double _psi = 1.0;
    for (size_t ind = 0; ind < ndim; ++ind) {
      _psi *= psivec[ind];
    }
    return _psi;
  }
//----------------------------------------------------------------------------------------------------------------
  std::vector<double> Densv(std::vector<double> pos) {
    //this is valid if the vectorial basis is orthogonal
    if (pos.size() != ndim) {throw std::string("ERROR: Translations.hpp: PartBox: Densv(): the vector introduced as argument does not have the correct dimension");}
    bool goodpos = true;
    for (size_t ind = 0; ind < ndim; ++ind) {
      if ((pos[ind] < 0.0)||(pos[ind] > L[ind])) {
        goodpos = false;
        break;
      }
    }
    std::vector<double> _probdens(ndim,0.0);
    if (goodpos) {
      double sine = 0.0;
      for (size_t ind = 0; ind < ndim; ++ind) {
        if (L[ind] != 0.0) {
          sine = sin(double(ni[ind])*pi*pos[ind]/L[ind]);
          _probdens[ind] = (2.0/L[ind])*sine*sine;
        }
        else {_probdens[ind] = 0.0;}                             //If L[ind] == 0.0 then consider the dimension inexistent
      }
    }
    return _probdens;
  }
  double Dens(std::vector<double> pos) {
    std::vector<double> probdens = Densv(pos);
    double _probdens = 1.0;
    for (size_t ind = 0; ind < ndim; ++ind) {
      _probdens *= probdens[ind];
    }
    return _probdens;
  }
//----------------------------------------------------------------------------------------------------------------
  std::vector<double> Probv(std::vector<double> ipos, std::vector<double> epos) {
    //this is valid if the vectorial basis is orthogonal.
    if (ipos.size() != ndim) {throw std::string("ERROR: Translations.hpp: PartBox: Probv(): the vector with the initial position does not have the correct dimension");}
    if (epos.size() != ndim) {throw std::string("ERROR: Translations.hpp: PartBox: Probv(): the vector with the end position does not have the correct dimension");}
    bool goodposi = true;
    for (size_t ind = 0; ind < ndim; ++ind) {
      if ((ipos[ind] < 0.0)||(ipos[ind] > L[ind])) {
        goodposi = false;
        break;
      }
    }
    bool goodpose = true;
    for (size_t ind = 0; ind < ndim; ++ind) {
      if ((epos[ind] < 0.0)||(epos[ind] > L[ind])) {
        goodpose = false;
        break;
      }
    }
    std::vector<double> _dens(ndim,0.0);
    if ((goodposi)&&(goodpose)) {
      double sineb = 0.0;
      double sinea = 0.0;
      for (size_t ind = 0; ind < ndim; ++ind) {
        if (L[ind] != 0.0) {
          sineb = sin((4.0*double(ni[ind])*pi*epos[ind])/L[ind]);
          sineb = sin((4.0*double(ni[ind])*pi*ipos[ind])/L[ind]);
          _dens[ind] = ((epos[ind]-ipos[ind])/L[ind]) - ((sineb-sinea)/(4.0*ni[ind])*pi);
        }
        else {_dens[ind] = 0.0;}                             //If L[ind] == 0.0 then consider the dimension inexistent
      }
    }
    return _dens;
  }
  double Prob(std::vector<double> ipos, std::vector<double> epos) {
    std::vector<double> dens = Probv(ipos,epos);
    double _dens = 1.0;
    for (size_t ind = 0; ind < ndim; ++ind) {
      _dens *= dens[ind];
    }
    return _dens;
  }
};
class PartBoxlim: public Boltzmann {
  //quantum m echanical description of particle in a box model, ndim dimensional
  double mass;
  size_t ndim;
  std::vector<double> L;
public:
  PartBoxlim(double _T, double _mass, size_t _ndim = 1) {
    this->T(_T);
    ndim = _ndim;
    mass = _mass;
    L.resize(ndim);
    for (size_t ind = 0; ind < ndim; ++ind) {
      L[ind] = 1.0;
    }
  }
  ~PartBoxlim() {L.erase(L.begin(),L.end());}
//----------------------------------------------------------------------------------------------------------------
  size_t nDimensions() {return ndim;}
  void nDimensions(size_t _ndim) {
    ndim = _ndim;
    L.resize(ndim);
    for (size_t ind = 0; ind < ndim; ++ind) {
      L[ind] = 1.0;
    }
  }
  void setL(std::vector<double> dimbox) {
    if (dimbox.size() == ndim) {L = dimbox;}
    else {throw std::string("ERROR: Translations.hpp: PartBoxlim: ProsetLbv(): the vector introduced as argument does not have the correct dimension.");}
  }
  const std::vector<double> & Lprint() {return L;}
  double Mass() {return mass;}
  void setMass(double _mass) {mass = _mass;}
//----------------------------------------------------------------------------------------------------------------
  double Volume() {
    double vol = 1.0;
    size_t Lpdim = L.size();
    for (size_t ind = 0; ind < Lpdim; ++ind) {
      if (L[ind] != 0.0) {vol *= L[ind];}
    }
    return vol;
  }
  void Volume(double vol) {
    size_t ndimen = this->nDimensions();
    double _l = pow(vol,1.0/double(ndimen));
    setL(std::vector<double> (ndimen,_l));
  }
//----------------------------------------------------------------------------------------------------------------
  std::vector<double> qXv() {
    double _qXv = sqrt((2.0*pi*mass*KB*T())/(h*h));
    size_t ndimen = this->nDimensions();
    std::vector<double> qXvector(ndimen,0.0);
    for (size_t ind = 0; ind < ndimen; ++ind) {
      qXvector[ind] = L[ind]*_qXv;
    }
    return qXvector;
  }
  double qX() {
    double _qX = 1.0;
    size_t ndimen = this->nDimensions();
    std::vector<double> vec = this->qXv();
    for (size_t ind = 0; ind < ndimen; ++ind) {
      _qX *= vec[ind];
    }
    return _qX;
  }
  double logqX() {return (log(this->qX()));}
  double dlogqXdT() {return (1.5/this->T());}
  double d2logqXdT2() {return (-1.5/(this->T()*this->T()));}
  double dlogqXdlogT() {return (1.5);}
  double dlogqXdV() {return (1.0/this->Volume());}
  double dlogqXdlogV() {return (1.0);}
};

#endif //_Translations_
