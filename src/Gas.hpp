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

#ifndef _Gas_
#define _Gas_
#include <string>
#include "rapidjson/json_interface.hpp"
#include "math/MatrixPackage.hpp"
#include "ConstantsPackage.hpp"
#include "Translations.hpp"
#include "Rotations.hpp"
#include "Vibrations.hpp"
#include "Molecule.hpp"

//description:
//statistical mechanical gas models

class PBlRRlHOEbase {
  //base class for gas model composed by:
  //  - Translations: particle in box in classical limit (PBl)
  //  - Rotations: rigid rotor in classical limit (RRl)
  //  - Vibrations: harmonic oscilator (HO)
  //  - Electrons: normal (E)
  //this is just the template model
protected:
  PartBoxlim translations;
  std::vector<RigidRotorlim> rotations;
  std::vector<HarmonicOscillator> vibrations;
  std::vector<Boltzmann> electrons;
  std::vector<Molecule> molecules;
  std::vector<double> Bavg;
  double N;
  double V;
  double T;
  double logNfact;
  double agrimme;                  //exponent in Grimme's weight factor
  double omega0;                   //reference frequency in Grimme's weight factor
  double DformH;                   //enthalpy of formation at 298.15 K
  bool grimmecorr;
  bool isatom;
  size_t nconf;
  size_t iconfepsmin;
public:
  PBlRRlHOEbase(double temperature, bool _Grimmecorr = true, double _N = NA, double volume = 0.0224, double mass = 0.0): translations(50.0,1.0,3) {
    //standard conditions is one mol in the volume of ideal gas at standard temperature and pressure
    //the constructor sets only basic information
    N = _N;
    Stirling();
    grimmecorr = _Grimmecorr;
    agrimme = 4.0;
    omega0 = 75.0;
    iconfepsmin = 0;
    V = volume;
    T = temperature;
    nconf = 0;
    isatom = false;
    setTranslations(mass,3);
  }
  ~PBlRRlHOEbase(){}
  //getters
  bool useGrimmeModel() {return grimmecorr;}
  size_t NumberConformers() {return nconf;}
  double aGrimme() {return agrimme;}
  double Omega0() {return omega0;}
  double Hformation() {return DformH;}
  size_t MinEnergyConformer() {return iconfepsmin;}
  //setters
  void setGrimmeModel(bool useGrimme) {grimmecorr = useGrimme;}
  void setaGrimme(double newag) {agrimme = newag;}
  void setOmega0(double newomega) {omega0 = newomega;}
  void setHformation(double hform) {DformH = hform;}
  void setNumberConformers(size_t nconf_new) {nconf = nconf_new;}
  void setMinEnergyConformer(size_t newmin) {iconfepsmin = newmin;}
  void setTranslations(double _mass, size_t _dim = 3) {
    //function that sets all required variables for translations class, except temperature and volume
    translations.setMass(_mass);
    translations.nDimensions(_dim);
    translations.T(T);
    translations.Volume(V);
  }
  //other functions
  void getNumberConformers(bool usemax = true) {
    size_t nrot = rotations.size();
    size_t nvibr = vibrations.size();
    size_t nelec = electrons.size();
    nconf = nrot;
    if (usemax) {       //search for largest of arrays
      if (nvibr > nconf) {nconf = nvibr;}
      if (nelec > nconf) {nconf = nelec;}
    }
    else {              //search for smallest of arrays
      if (nvibr < nconf) {nconf = nvibr;}
      if (nelec < nconf) {nconf = nelec;}
    }
  }
  void Stirling() {
    //calculates the Stirling correction factor log(N!)
    logNfact = N*(log(N) - 1.0);
  }
  void changeT(double _T) {
    //function that changes the temperature
    T = _T;
    translations.T(_T);
    for (size_t iconf = 0; iconf < nconf; ++iconf) {
      rotations[iconf].T(_T);
      vibrations[iconf].T(_T);
      electrons[iconf].T(_T);
    }
  }
  void changeV(double _V) {
    //function that changes the volume
    V = _V;
    translations.Volume(_V);
  }
  void changeN(double _N) {
    //function that changes the number of particles
    N = _N;
    Stirling();
  }
  void setNatoms(int _natom) {
    if (_natom == 1) {isatom = true;}
    else {isatom = false;}
  }
  void SearchMinEnergyConformer() {
    //function that searches for the conformer with minimum energy
    double epsmin = electrons[0].Eshift() + vibrations[0].ZPVE();
    iconfepsmin = 0;
    double eps;
    for (size_t iconf = 1; iconf < nconf; ++iconf) {
      eps = electrons[iconf].Eshift() + vibrations[iconf].ZPVE();
      if (eps < epsmin) {
        epsmin = eps;
        iconfepsmin = iconf;
      }
    }
  }
  void PrintEnergy(int iconf = -1, std::string unit = "Eh") {
    double convert = 1.0/au2J;
    if (unit == "J/mol") {convert = NA;}
    else if (unit == "kJ/mol") {convert = 0.001*NA;}
    else if (unit == "kcal/mol") {convert = 0.001*NA/cal2J;}
    if (iconf > -1) {std::cout << "Energy(" << iconf << ") = " << (electrons[iconf].Eshift() + vibrations[iconf].ZPVE())*convert << std::endl;}
    else {
      for (size_t icnf = 0; icnf < nconf; ++icnf) {
        std::cout << "Energy(" << icnf << ") = " << (electrons[icnf].Eshift() + vibrations[icnf].ZPVE())*convert << std::endl;
      }
    }
  }
  void addRotationModel(const std::vector<double> & _inertia, Molecule & gas, std::string _rotsymm = "0", bool increasenconf = false) {
    //addition of a rotation model
    RigidRotorlim single_rotation(T);
    if (increasenconf) {++nconf;}
    //get inertia's eigenvectors
    double sumI = 0.0;
    std::vector<double> IVec;
    for (size_t ival = 0; ival < _inertia.size(); ++ival) {
      sumI += _inertia[ival];
    }
    if ((sumI == 0)||(_inertia.size() != 3)) {
      //if there is no inertia read from db, then obtain it
      matrixE IMat = gas.Inertia();
      IVec = MatDiag(IMat);
      if ((IMat.rows() != 3)||(IMat.cols() != 3)) {
        IVec.push_back(IMat(1,1));
        IVec.push_back(IMat(1,1));
      }
    }
    else {IVec = _inertia;}
    Bavg[nconf - 1] = pow(IVec[0]*IVec[1]*IVec[2],1.0/3.0);
    if (Bavg[nconf - 1] < 1.0e-44) {Bavg[nconf - 1] = (IVec[0] + IVec[1] + IVec[2])/3.0;}
    single_rotation.setInertia(IVec);
    //get symmetry related info
    single_rotation.setSigma(PointGroup2Sigma("C1"));
    if (_rotsymm != "0") {single_rotation.SymmGroup(_rotsymm);}
    else {single_rotation.SymmGroup();}
    //printout
    std::cout << "rigid rotor: " << std::endl;
    std::cout << "   symmetry group    " << single_rotation.getSymmGroup() << std::endl;
    std::cout << "   sigma             " << single_rotation.getSigma() << std::endl;
    std::cout << "   inertia           ";
    std::cout << std::scientific;
    for (size_t idx = 0; idx < single_rotation.getInertia().size(); ++idx) {
      std::cout << single_rotation.getInertia()[idx] << " ";
    }
    std::cout << std::fixed;
    std::cout << std::endl;
    rotations.push_back(single_rotation);
  }
  void addVibrationModel(const std::vector<double> & _vibrations) {
    //addition of a vibration model
    HarmonicOscillator single_vibration(T);
    //set vibrations and get the thetas
    single_vibration.setVibrations(_vibrations);
    vibrations.push_back(single_vibration);
    //printout
    std::cout << "harmonic oscillator: " << std::endl;
    std::cout << "   vibrations        ";
    std::vector<double> vibraux = single_vibration.getVibrations();
    size_t nvibrations = vibraux.size();
    for (size_t idx = 0; idx < nvibrations; ++idx) {
      std::cout << vibraux[idx] << " ";
    }
    std::cout << std::endl;
    std::cout << "   ThetaVibr         ";
    vibraux = single_vibration.getThetas();
    for (size_t idx = 0; idx < nvibrations; ++idx) {
      std::cout << vibraux[idx] << " ";
    }
    std::cout << std::endl;
  }
  void addElectronModel(const std::vector<double> & _Eel, const std::vector<double> & _gel) {
    //addition of an electronic model
    Boltzmann single_electron;
    //get electrons
    single_electron.setEnergy(_Eel,"au");
    single_electron.setgE(_gel);
    electrons.push_back(single_electron);
    std::cout << "Electronic energies: " << std::endl;
    std::vector<double> ge = single_electron.GE();
    std::vector<double> ee = single_electron.Energy();
    size_t nelectronic = ge.size();
    double eps0 = single_electron.E0();
    for (size_t idx = 0; idx < nelectronic; ++idx) {
      std::cout << "(" << ge[idx] << ") " << (eps0 + ee[idx])/au2J << "  ";
    }
    std::cout << std::endl;
  }
  void addConformer(std::string gfile, const std::vector<double> & _inertia, const std::vector<double> & _vibrations, const std::vector<double> & _Eel, const std::vector<double> & _gel, int _chg = 0, size_t _mult = 1, std::string _symm = "0", std::string _rotsymm = "0", double thresholdzero = 1.0e-10) {
    //function that adds a conformer in full
    //allocate temporary particle objects
    Molecule gas(gfile,_chg,_mult,_symm);
    molecules.push_back(gas);
    Bavg.push_back(1.0e-44);
    ++nconf;
    double massconf = translations.Mass();
    if (massconf < thresholdzero) {translations.setMass(gas.m());}    //set new mass
    else {                                                            //compare masses, these must match
      if (fabs(massconf - gas.m()) > thresholdzero) {
        std::cout << "WARNING: Gas.hpp: PBlRRlHOEbase: addConformer(): added conformation has mass inconsistent with the rest of the ensemble\n";
      }
    }
    //check whether we have an atom
    size_t NAtoms = 0;
    if (gfile != "") {NAtoms = gas.Natoms();}
    this->setNatoms(NAtoms);
    std::cout << "print output:" << std::endl;
    std::cout << "gas constant:     " << NA*KB << std::endl;
    std::cout << "number molecules: " << N << std::endl;
    std::cout << "particle in box: " << std::endl;
    std::cout << "   n dimensions      " << translations.nDimensions() << std::endl;
    std::cout << std::scientific;
    std::cout << "   mass              " << translations.Mass() << std::endl;
    std::cout << std::fixed;
    std::cout << "   volume            " << translations.Volume() << std::endl;
    std::cout << std::endl;
    //add rest of models
    this->addRotationModel(_inertia,gas,_rotsymm,false);
    this->addVibrationModel(_vibrations);
    this->addElectronModel(_Eel,_gel);
    //check for minimum in energy
    massconf = electrons[iconfepsmin].Eshift() + vibrations[iconfepsmin].ZPVE();
    double eps = electrons[nconf - 1].Eshift() + vibrations[nconf - 1].ZPVE();
    if (eps < massconf) {iconfepsmin = nconf - 1;}
  }
  void removeConformer(int pos = -1) {
    //function that removes a certain conformer from the ensemble
    if (pos < 0) {pos = nconf - 1;}
    rotations.erase(rotations.begin() + pos);
    vibrations.erase(vibrations.begin() + pos);
    electrons.erase(electrons.begin() + pos);
    molecules.erase(molecules.begin() + pos);
    Bavg.erase(Bavg.begin() + pos);
    --nconf;
  }
  void ReadDB(std::vector<std::string> method, int _chg = 0, size_t _mult = 1, std::string db = "main", double thresholdzero = 1.0e-10) {
    //auxiliary variables and arrays
    std::vector<std::string> outstrng;
    std::vector<std::string> * poutstrng = & outstrng;
    std::vector<double> outdouble;
    std::vector<double> * poutdouble = & outdouble;
    std::vector<std::vector<double> > outvector;
    std::vector<std::vector<double> > * poutvector = & outvector;
    int nconformations;
    //getting information about molecules
    PBlRRlHOEJSON(db,method,poutstrng,poutdouble,poutvector,nconformations);
    std::string symm = outstrng[3];
    std::string rotsymm = outstrng[4];
    //generating molecules
    //outstrng[3*idx - 1] -> geometries
    //outvector[4*idx - 4] -> inertia
    //outvector[4*idx - 3] -> vibrations
    //outvector[4*idx - 2] -> Eelectronic
    //outvector[4*idx - 1] -> gE
    DformH = outdouble[0];
    for (size_t iconf = 1; iconf < nconformations + 1; ++iconf) {
      this->addConformer("../db/geom/" + outstrng[3*iconf - 1],outvector[4*iconf - 4],outvector[4*iconf - 3],outvector[4*iconf - 2],outvector[4*iconf - 1],_chg,_mult,symm,rotsymm,thresholdzero);
    }
  }
  virtual double qX() {
    //get the partition function for a multistate "single particle"
    return translations.qX();
  }
  virtual double logqX() {
    //get the logarithm of the partition function for a multistate "single particle"
    return translations.logqX();
  }
  virtual double dlogqXdT() {
    //get the first derivative of the logarithm of the partition function for a multistate "single particle" with respect to temperature
    return translations.dlogqXdT();
  }
  virtual double dlogqXdV() {
    //get the first derivative of the logarithm of the partition function for a multistate "single particle" with respect to volume
    return translations.dlogqXdV();
  }
  virtual double d2logqXdT2() {
    //get the second derivative of the logarithm of the partition function for a multistate "single particle" with respect to temperature
    return translations.d2logqXdT2();
  }
  double Usingle(size_t iconf) {
    double u = N*(translations.U() + electrons[iconf].U() + electrons[iconf].Eshift());
    if (!isatom) {u += N*(rotations[iconf].U() + vibrations[iconf].U(grimmecorr,Bavg[iconf],omega0,agrimme) + vibrations[iconf].ZPVE());}
    return u;
  }
  double Hsingle(size_t iconf) {
    double enth = N*(translations.H() + electrons[iconf].H() + electrons[iconf].Eshift());
    if (!isatom) {enth += N*(rotations[iconf].H() + vibrations[iconf].H(grimmecorr,Bavg[iconf],omega0,agrimme) + vibrations[iconf].ZPVE());}
    return enth;
  }
  double Ssingle(size_t iconf) {
    double s = N*(translations.S() + electrons[iconf].S()) - KB*logNfact;
    if (!isatom) {s += N*(rotations[iconf].S() + vibrations[iconf].S(grimmecorr,Bavg[iconf],omega0,agrimme));}
    return s;
  }
  double Gsingle(size_t iconf) {
    double g = N*(translations.G() + electrons[iconf].G() + electrons[iconf].Eshift()) + KB*T*logNfact;
    if (!isatom) {g += N*(rotations[iconf].G() + vibrations[iconf].G(grimmecorr,Bavg[iconf],omega0,agrimme) + vibrations[iconf].ZPVE());}
    return g;
  }
  double Asingle(size_t iconf) {
    double a = KB*T*logNfact + N*(translations.A() + electrons[iconf].A() + electrons[iconf].Eshift());
    if (!isatom) {a += N*(rotations[iconf].A() + vibrations[iconf].A(grimmecorr,Bavg[iconf],omega0,agrimme) + vibrations[iconf].ZPVE());}
    return a;
  }
  double CPsingle(size_t iconf) {
    double cp = N*(translations.CV() + KB + electrons[iconf].CV());
    if (!isatom) {cp += N*(rotations[iconf].CV() + vibrations[iconf].CV(grimmecorr,Bavg[iconf],omega0,agrimme));}
    return cp;
  }
  double CVsingle(size_t iconf) {
    double cv = N*(translations.CV() + electrons[iconf].CV());
    if (!isatom) {cv += N*(rotations[iconf].CV() + vibrations[iconf].CV(grimmecorr,Bavg[iconf],omega0,agrimme));}
    return cv;
  }
  double Hshift() {
    //enthalpy in J/mol with a shift so that H(298.15) = DformH
    double oldT = T;
    changeT(298.15);
    double shift = DformH - this->H();
    changeT(T);
    return shift;
  }  //these functions are temperature specific, therefore I do not save the partition function between calls, since that is also temperature dependent
  virtual double U() {
    //internal energy in J/mol
    return Usingle(0);
  }
  virtual double H() {
    //enthalpy in J/mol
    return Hsingle(0);
  }
  virtual double S() {
    //entropy in J/(K.mol)
    return Ssingle(0);
  }
  virtual double G() {
    //Gibbs energy in J/mol
    return Gsingle(0);
  }
  virtual double A() {
    //Helmoltz energy in J/mol
    return Asingle(0);
  }
  virtual double CP() {
    //CP in J/(K.mol)
    return CPsingle(0);
  }
  virtual double CV() {
    //CV in J/(K.mol)
    return CVsingle(0);
  }
  virtual double P(int iconf = -1) {
    //pressure in Pa
    double p = N*translations.P();
    return p;
  }
};
class PBlRRlHOE: public PBlRRlHOEbase {
  //this is a single structure model
public:
  PBlRRlHOE(double temperature = 298.15, bool _Grimmecorr = true, double _N = NA, double volume = 0.0224): PBlRRlHOEbase(temperature,_Grimmecorr,_N,volume,0.0) {}
  PBlRRlHOE(double temperature, std::string _file, std::vector<double> _inertia, std::vector<double> _vibrations, std::vector<double> _Eel, std::vector<double> _gel, int _chg = 0, size_t _mult = 1, std::string _symm = "0", std::string _rotsymm = "0", bool _Grimmecorr = true, double _N = NA, double volume = 0.0224, double thresholdzero = 1.0e-8):
  PBlRRlHOEbase(temperature,_Grimmecorr,_N,volume,0.0) {
    addConformer(_file,_inertia,_vibrations,_Eel,_gel,_chg,_mult,_symm,_rotsymm,thresholdzero);
  }
  ~PBlRRlHOE(){}
  void setRotations(std::string _symm, std::vector<double> IVec) {
    //function that sets all required variables for rotations class, except temperature
    rotations[0].setSigma(PointGroup2Sigma(_symm));
    rotations[0].setInertia(IVec);
  }
  void setVibrations(std::vector<double> _vibrations) {
    //function that sets all required variables for vibrations class, except temperature
    vibrations[0].setVibrations(_vibrations);
  }
  void setElectrons(std::vector<double> _Eel, std::vector<double> _gel) {
    //function that sets all required variables for electron class, except temperature
    electrons[0].setEnergy(_Eel,"au");
    electrons[0].setgE(_gel);
  }
  double TSfrequency() {
    //function returning the negative frequency (first showing up in the frequency vector)
    return vibrations[0].TSfrequency();
  }
};
class MSPBlRRlHOE: public PBlRRlHOEbase {
  //this is the multistructure model
public:
  MSPBlRRlHOE(double temperature, bool _Grimmecorr = true, double _N = NA, double volume = 0.0224, double mass = 0.0): PBlRRlHOEbase(temperature,_Grimmecorr,_N,volume,mass) {}
  ~MSPBlRRlHOE(){}
  double TSfrequency(int iconf = -1) {
    //function returning the negative frequency (first showing up in the frequency vector)
    double ifreq = 0.0;
    if (iconf > -1) {ifreq = vibrations[iconf].TSfrequency();}
    else {}  //to be implemented
    return ifreq;
  }
  //these functions are temperature specific, therefore I do not save the partition function between calls, since that is also temperature dependent
  double Eshift() {
    //function that averages the energy shifts for all conformers, to be added to all energies
    double eweight = 0.0;
    double weight = 0.0;
    double beta = 1.0/(KB*T);
    double eps;
    double epsmin = electrons[iconfepsmin].Eshift() + vibrations[iconfepsmin].ZPVE();
    for (size_t iconf = 0; iconf < nconf; ++iconf) {
      eps = electrons[iconf].Eshift() + vibrations[iconf].ZPVE();
      eweight += eps*exp(-(eps - epsmin)*beta);
      weight += exp(-(eps - epsmin)*beta);
    }
    return eweight/weight;
  }
  double qRVE() {
    //get the partition function coupled for rotations, vibrations and electrons
    double qrve = 0.0;
    double eps;
    double epsmin = electrons[iconfepsmin].Eshift() + vibrations[iconfepsmin].ZPVE();
    double beta = 1.0/(KB*T);
    double weight = 0.0;
    for (size_t iconf = 0; iconf < nconf; ++iconf) {
      eps = electrons[iconf].Eshift() + vibrations[iconf].ZPVE();
      weight = exp(-(eps - epsmin)*beta);
      qrve += weight*rotations[iconf].qX()*vibrations[iconf].qX(grimmecorr,omega0,agrimme,Bavg[iconf])*electrons[iconf].qX();
    }
    return qrve;
  }
  double logqRVE() {
    //get the logarithm of the partition function coupled for rotations, vibrations and electrons
    return log(this->qRVE());
  }
  double dqRVEdT() {
    //get the first derivative of the partition function coupled for rotations, vibrations and electrons with respect to the temperature
    double dqrvedt = 0.0;
    double eps;
    double epsmin = electrons[iconfepsmin].Eshift() + vibrations[iconfepsmin].ZPVE();
    double depsbeta;
    double tempm1 = 1.0/T;
    double beta = 1.0/(KB*T);
    double weight = 0.0;
    double qr;
    double qv;
    double qe;
    std::vector<double> qvder(2,0.0);
    for (size_t iconf = 0; iconf < nconf; ++iconf) {
      eps = electrons[iconf].Eshift() + vibrations[iconf].ZPVE();
      depsbeta = (eps - epsmin)*beta;
      weight = exp(-depsbeta);
      qr = rotations[iconf].qX();
      vibrations[iconf].getDerivatives(1,qvder,grimmecorr,omega0,agrimme,Bavg[iconf]);
      qe = electrons[iconf].qX();
      dqrvedt += weight*(rotations[iconf].dqXdT()*qvder[0]*qe + qr*qvder[1]*qe + qr*qvder[0]*electrons[iconf].dqXdT() + depsbeta*tempm1*qr*qvder[0]*qe);
    }
    return dqrvedt;
  }
  double d2qRVEdT2() {
    //get the second derivative of the partition function coupled for rotations, vibrations and electrons with respect to the temperature
    double d2qrvedt2 = 0.0;
    double eps;
    double epsmin = electrons[iconfepsmin].Eshift() + vibrations[iconfepsmin].ZPVE();
    double depsbeta;
    double tempm1 = 1.0/T;
    double KBT = KB*T;
    double beta = 1.0/KBT;
    double KBT3m1 = beta*tempm1*tempm1;
    double weight = 0.0;
    double qr;
    double qv;
    double qe;
    double qrp;
    double qep;
    std::vector<double> qvder(3,0.0);
    for (size_t iconf = 0; iconf < nconf; ++iconf) {
      eps = electrons[iconf].Eshift() + vibrations[iconf].ZPVE();
      depsbeta = (eps - epsmin)*beta;
      weight = exp(-depsbeta);
      qr = rotations[iconf].qX();
      qrp = rotations[iconf].dqXdT();
      vibrations[iconf].getDerivatives(2,qvder,grimmecorr,omega0,agrimme,Bavg[iconf]);
      qe = electrons[iconf].qX();
      qep = electrons[iconf].dqXdT();
      d2qrvedt2 += weight*(rotations[iconf].d2qXdT2()*qvder[0]*qe + qr*qvder[2]*qe + qr*qvder[0]*electrons[iconf].d2qXdT2() + qr*qvder[0]*qe*depsbeta*KBT3m1*(eps - epsmin - 2.0*KBT));
      d2qrvedt2 += 2.0*weight*(qrp*qvder[1]*qe + qr*qvder[1]*qep + qrp*qvder[0]*qep + depsbeta*tempm1*(qrp*qvder[0]*qe + qr*qvder[1]*qe + qr*qvder[0]*qep));
    }
    return d2qrvedt2;
  }
  double qX() {
    //get the partition function for a multistate "single particle"
    return this->qRVE()*translations.qX();
  }
  double logqX() {
    //get the logarithm of the partition function for a multistate "single particle"
    return translations.logqX() + this->logqRVE();
  }
  double dlogqXdT() {
    //get the first derivative of the logarithm of the partition function for a multistate "single particle" with respect to temperature
    double qrve = 0.0;
    double dqrvedt = 0.0;
    double qr;
    double qe;
    double eps;
    double epsmin = electrons[iconfepsmin].Eshift() + vibrations[iconfepsmin].ZPVE();
    double depsbeta;
    double tempm1 = 1.0/T;
    double beta = 1.0/(KB*T);
    double weight = 0.0;
    std::vector<double> qvder(2,0.0);
    for (size_t iconf = 0; iconf < nconf; ++iconf) {
      eps = electrons[iconf].Eshift() + vibrations[iconf].ZPVE();
      depsbeta = (eps - epsmin)*beta;
      weight = exp(-depsbeta);
      qr = rotations[iconf].qX();
      vibrations[iconf].getDerivatives(1,qvder,grimmecorr,omega0,agrimme,Bavg[iconf]);
      qe = electrons[iconf].qX();
      qrve += weight*qr*qvder[0]*qe;
      dqrvedt += weight*(rotations[iconf].dqXdT()*qvder[0]*qe + qr*qvder[1]*qe + qr*qvder[0]*electrons[iconf].dqXdT() + depsbeta*tempm1*qr*qvder[0]*qe);
    }
    return translations.dlogqXdT() + dqrvedt/qrve;
  }
  double dlogqXdV() {
    //get the first derivative of the logarithm of the partition function for a multistate "single particle" with respect to volume
    return translations.dlogqXdV();
  }
  double d2logqXdT2() {
    //get the second derivative of the logarithm of the partition function for a multistate "single particle" with respect to temperature
    double qrve = 0.0;
    double dqrvedt = 0.0;
    double d2qrvedt2 = 0.0;
    double qr;
    double qe;
    double qrp;
    double qep;
    double eps;
    double epsmin = electrons[iconfepsmin].Eshift() + vibrations[iconfepsmin].ZPVE();
    double depsbeta;
    double tempm1 = 1.0/T;
    double KBT = KB*T;
    double beta = 1.0/KBT;
    double KBT3m1 = beta*tempm1*tempm1;
    double weight = 0.0;
    std::vector<double> qvder(3,0.0);
    for (size_t iconf = 0; iconf < nconf; ++iconf) {
      eps = electrons[iconf].Eshift() + vibrations[iconf].ZPVE();
      depsbeta = (eps - epsmin)*beta;
      weight = exp(-depsbeta);
      qr = rotations[iconf].qX();
      qe = electrons[iconf].qX();
      qrp = rotations[iconf].dqXdT();
      qep = electrons[iconf].dqXdT();
      vibrations[iconf].getDerivatives(2,qvder,grimmecorr,omega0,agrimme,Bavg[iconf]);
      qrve += weight*qr*qvder[0]*qe;
      dqrvedt += weight*(qrp*qvder[0]*qe + qr*qvder[1]*qe + qr*qvder[0]*qep + depsbeta*tempm1*qr*qvder[0]*qe);
      d2qrvedt2 += weight*(rotations[iconf].d2qXdT2()*qvder[0]*qe + qr*qvder[2]*qe + qr*qvder[0]*electrons[iconf].d2qXdT2() + qr*qvder[0]*qe*depsbeta*KBT3m1*(eps - epsmin - 2.0*KBT));
      d2qrvedt2 += 2.0*weight*(qrp*qvder[1]*qe + qr*qvder[1]*qep + qrp*qvder[0]*qep + depsbeta*tempm1*(qrp*qvder[0]*qe + qr*qvder[1]*qe + qr*qvder[0]*qep));
    }
    return translations.d2logqXdT2() + d2qrvedt2/qrve - dqrvedt*dqrvedt/(qrve*qrve);
  }
  double U(int iconf = -1) {
    //internal energy in J/mol
    double u = 0.0;
    if (iconf > -1) {
      u = N*(translations.U() + electrons[iconf].U() + electrons[iconf].Eshift());
      if (!isatom) {u += N*(rotations[iconf].U() + vibrations[iconf].U(grimmecorr,Bavg[iconf],omega0,agrimme) + vibrations[iconf].ZPVE());}
    }
    else {u = N*(KB*T*T*this->dlogqXdT() + this->Eshift());}
    return u;
  }
  double H(int iconf = -1) {
    //enthalpy in J/mol
    double enth = 0.0;
    if (iconf > -1) {
      enth = N*(translations.H() + electrons[iconf].H() + electrons[iconf].Eshift());
      if (!isatom) {enth += N*(rotations[iconf].H() + vibrations[iconf].H(grimmecorr,Bavg[iconf],omega0,agrimme) + vibrations[iconf].ZPVE());}
    }
    else {enth = N*(KB*T*(T*this->dlogqXdT() + V*this->dlogqXdV()) + this->Eshift());}
    return enth;
  }
  double S(int iconf = -1) {
    //entropy in J/(K.mol)
    double s = 0.0;
    if (iconf > -1) {
      s = N*(translations.S() + electrons[iconf].S()) - KB*logNfact;
      if (!isatom) {s += N*(rotations[iconf].S() + vibrations[iconf].S(grimmecorr,Bavg[iconf],omega0,agrimme));}
    }
    else {s = N*KB*(T*this->dlogqXdT() + this->logqX()) - KB*logNfact;}
    return s;
  }
  double G(int iconf = -1) {
    //Gibbs energy in J/mol
    double g = 0.0;
    if (iconf > -1) {
      g = N*(translations.G() + electrons[iconf].G() + electrons[iconf].Eshift()) + KB*T*logNfact;
      if (!isatom) {g += N*(rotations[iconf].G() + vibrations[iconf].G(grimmecorr,Bavg[iconf],omega0,agrimme) + vibrations[iconf].ZPVE());}
    }
    else {g = N*(KB*T*(V*this->dlogqXdV() - this->logqX()) + this->Eshift()) + KB*T*logNfact;}
    return g;
  }
  double A(int iconf = -1) {
    //Helmoltz energy in J/mol
    double a = 0.0;
    if (iconf > -1) {
      a = KB*T*logNfact + N*(translations.A() + electrons[iconf].A() + electrons[iconf].Eshift());
      if (!isatom) {a += N*(rotations[iconf].A() + vibrations[iconf].A(grimmecorr,Bavg[iconf],omega0,agrimme) + vibrations[iconf].ZPVE());}
    }
    else {a = KB*T*logNfact + N*(this->Eshift() - KB*T*this->logqX());}
    return a;
  }
  double CP(int iconf = -1) {
    //CP in J/(K.mol)
    double cp = 0.0;
    if (iconf > -1) {
      cp = N*(translations.CV() + KB + electrons[iconf].CV());
      if (!isatom) {cp += N*(rotations[iconf].CV() + vibrations[iconf].CV(grimmecorr,Bavg[iconf],omega0,agrimme));}
    }
    else {cp = N*KB*(1.0 + 2.0*T*this->dlogqXdT() + T*T*this->d2logqXdT2());}
    return cp;
  }
  double CV(int iconf = -1) {
    //CV in J/(K.mol)
    double cv = 0.0;
    if (iconf > -1) {
      cv = N*(translations.CV() + electrons[iconf].CV());
      if (!isatom) {cv += N*(rotations[iconf].CV() + vibrations[iconf].CV(grimmecorr,Bavg[iconf],omega0,agrimme));}
    }
    else {cv = N*KB*T*(2.0*this->dlogqXdT() + T*this->d2logqXdT2());}
    return cv;
  }
};
class MSavgPBlRRlHOE: public PBlRRlHOEbase {
  //this is a multistructure model that weight averages the thermodynamic properties
  std::vector<PBlRRlHOE> conformers;            //contains series of molecules that must all be conformers
  std::vector<double> cweights;                 //contains the Boltzmann weights for each conformer
  double Told;
  double Vold;
  bool gibbsweight;
public:
  MSavgPBlRRlHOE(double temperature, bool _Grimmecorr = true, double _N = NA, double volume = 0.0224, double mass = 0.0): PBlRRlHOEbase(temperature,_Grimmecorr,_N,volume,mass) {
    Told = temperature - 1.0;
    Vold = volume - 1.0;
    gibbsweight = true;
  }
  bool GibbsWeight() {return gibbsweight;}
  void setGibbsWeight(bool _gibbsweight) {gibbsweight = _gibbsweight;}
  double TSfrequency() {
    //function returning the negative frequency of the first conformer (assume there is oly one conformer)
    return conformers[0].TSfrequency();
  }
  void setWeights() {
    //interface for setting weights
    if (gibbsweight) {this->GWeights();}
    else {this->AWeights();}
  }
  void GWeights() {
    //function calculating the weight of conformer
    //weight is given from Gibbs averaging!
    cweights.resize(nconf);
    double g0 = Gsingle(0);
    double sumgs = 0.0;
    double expx;
    for (size_t icnf = 0; icnf < nconf; ++icnf) {
      expx = exp(-(Gsingle(icnf) - g0)/(N*KB*T));
      sumgs += expx;
      cweights[icnf] = expx;
    }
    //renormalizing
    for (size_t icnf = 0; icnf < nconf; ++icnf) {
      cweights[icnf] /= sumgs;
    }
  }
  void AWeights() {
    //function calculating the weight of conformer
    //weight is given from Helmoltz averaging!
    cweights.resize(nconf);
    double a0 = Asingle(0);
    double sumgs = 0.0;
    double expx;
    for (size_t icnf = 0; icnf < nconf; ++icnf) {
      expx = exp(-(Asingle(icnf) - a0)/(N*KB*T));
      sumgs += expx;
      cweights[icnf] = expx;
    }
    //renormalizing
    for (size_t icnf = 0; icnf < nconf; ++icnf) {
      cweights[icnf] /= sumgs;
    }
  }
  double U(int iconf = -1) {
    //internal energy in J/mol
    if ((Told != T)||(Vold != V)) {this->setWeights();}
    double u = 0.0;
    if (iconf > -1) {u = Usingle(iconf);}
    else {
      for (size_t icnf = 0; icnf < nconf; ++icnf) {
        u += Usingle(icnf)*cweights[icnf];
      }
    }
    return u;
  }
  double H(int iconf = -1) {
    //enthalpy in J/mol
    if ((Told != T)||(Vold != V)) {this->setWeights();}
    double enth = 0.0;
    if (iconf > -1) {enth = Hsingle(iconf);}
    else {
      for (size_t icnf = 0; icnf < nconf; ++icnf) {
        enth += Hsingle(icnf)*cweights[icnf];
      }
    }
    return enth;
  }
  double S(int iconf = -1) {
    //entropy in J/(K.mol)
    if ((Told != T)||(Vold != V)) {this->setWeights();}
    double s = 0.0;
    if (iconf > -1) {s = Ssingle(iconf);}
    else {
      for (size_t icnf = 0; icnf < nconf; ++icnf) {
        s += Ssingle(icnf)*cweights[icnf] - N*KB*cweights[icnf]*log(cweights[icnf]);
      }
    }
    return s;
  }
  double G(int iconf = -1) {
    //Gibbs energy in J/mol
    if ((Told != T)||(Vold != V)) {this->setWeights();}
    double g = 0.0;
    if (iconf > -1) {g = Gsingle(iconf);}
    else {
      for (size_t icnf = 0; icnf < nconf; ++icnf) {
        g += Gsingle(icnf)*cweights[icnf] + N*KB*T*cweights[icnf]*log(cweights[icnf]);
      }
    }
    return g;
  }
  double A(int iconf = -1) {
    //Helmoltz energy in J/mol
    if ((Told != T)||(Vold != V)) {this->setWeights();}
    double a = 0.0;
    if (iconf > -1) {a = Asingle(iconf);}
    else {
      for (size_t icnf = 0; icnf < nconf; ++icnf) {
        a += Asingle(icnf)*cweights[icnf] + N*KB*T*cweights[icnf]*log(cweights[icnf]);
      }
    }
    return a;
  }
  double CP(int iconf = -1) {
    //CP in J/(K.mol)
    if ((Told != T)||(Vold != V)) {this->setWeights();}
    double cp = 0.0;
    if (iconf > -1) {cp = CPsingle(iconf);}
    else {
      for (size_t icnf = 0; icnf < nconf; ++icnf) {
        cp += CPsingle(icnf)*cweights[icnf];
      }
    }
    return cp;
  }
  double CV(int iconf = -1) {
    //CV in J/(K.mol)
    if ((Told != T)||(Vold != V)) {this->setWeights();}
    double cv = 0.0;
    if (iconf > -1) {cv = CVsingle(iconf);}
    else {
      for (size_t icnf = 0; icnf < nconf; ++icnf) {
        cv += CVsingle(icnf)*cweights[icnf];
      }
    }
    return cv;
  }
};

#endif //_Gas_
