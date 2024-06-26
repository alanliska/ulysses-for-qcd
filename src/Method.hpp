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

#ifndef _METHOD_METHODS_
#define _METHOD_METHODS_
#include <string>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include "Hessian.hpp"
#include "Molecule.hpp"
#include "Intermolecular.hpp"
#include "UnitConversion.hpp"
#include "atoms/ChargePackage.hpp"
#include "atoms/AtomPackage.hpp"
#include "math/MatrixPackage.hpp"
#include "math/VectorPackage.hpp"
#include "math/SearchPackage.hpp"
#include "math/SolverPackage.hpp"

//descrition:
//base class for any method

class Method {
protected:
  std::vector<size_t> atoms;                          //vector with atomic numbers
  std::vector<double> QAtoms;                         //atomic charges
  std::vector<double> rAB;
  matrixE BWilson;                                    //Wilson's B matrix to change to internal coordinates
  matrixE auxMAT;
  matrixE auxMAT3;
  Molecule mol;                                       //molecule object
  int Natoms;                                         //total number of atoms
  int print;                                          //determine whether to print
  int restart;                                        //control whether to restart the calculation is accepted
  double etot;                                        //total energy
  double corecorrection;                              //empirical correction to quantum chemistry, e.g., from dispersion or hydrogen bonding
public:
  double thresholdzero;                               //whenever needed, this is the numerical zero
  Method(Molecule & _mol): mol(_mol) {
    etot = 0.0;
    corecorrection = 0.0;
    atoms = mol.Atoms();
    Natoms = atoms.size();
    thresholdzero = 1.0e-10;
    print = 1;
    restart = 0;
    rAB.resize(3);
  }
  ~Method() {}
  //getters
  virtual std::string Type() {return "Method";}
  virtual std::string MechanicsType() {return "none";}
  matrixE WilsonBmatrix() {return BWilson;}
  std::vector<double> getQAtoms() {return QAtoms;}
  matrixE Geometry() {return mol.Geometry();}
  std::vector<size_t> Atoms() {return atoms;}
  size_t NAtoms() {return Natoms;}
  Molecule & Component() {return mol;}
  double getEnergy(bool complete = 0) {
    double enrg = etot;
    if (complete) {enrg += corecorrection;}
    return enrg;
  }
  virtual void calcEnergy() {}
  virtual void calcEnergy(matrixE & geom) {}
  double getThresholdZero() {return thresholdzero;}
  int Restart() {return restart;}
  //setters
  virtual void setCheck4Aromatics(bool newvar) {}
  virtual void ReSet(Molecule & _newmol, matrix<int> & newc, matrixE & newbo) {}
  virtual void setRingMatrix(matrix<int> & newrmat) {}
  void setQAtoms(std::vector<double> & _QAtoms) {QAtoms = _QAtoms;}
  virtual void setMolecule(Molecule & _component) {
    mol = _component;
    atoms = mol.Atoms();
    Natoms = atoms.size();
  }
  virtual void setMoleculeStructure(matrixE & newgeom, std::vector<size_t> & newatoms) {
    mol.setAtoms(newatoms);
    mol.setGeometry(newgeom);
    mol.calcNatoms();
    mol.setCharge(0);
    mol.CountElectrons();
    mol.setMultiplicity(1);
    mol.masses();
    mol.setPGroup("C1");
  }
  void updGeometry(matrixE & displacement, bool project = 0) {
    //function that adds a change in geometry to the already existing one
    //changes assuming that only atomic positions are modified
    if ((project)&&(BWilson.rows()*BWilson.cols() > 1)) {
      //if project and we have Bmatrix, then use it to "fix" the geometry
      auxMAT = displacement;
      displacement = BWilson*auxMAT;
    }
    auxMAT = mol.Geometry();
    int nrows = auxMAT.rows();
    int ncols = auxMAT.cols();
    int count = 1;
    for (size_t idatom = 0; idatom < nrows; ++idatom) {
      for (size_t idcoord = 0; idcoord < ncols; ++idcoord, ++count) {
        auxMAT(idatom + 1,idcoord + 1) -= displacement(count,1);
      }
    }
    this->setGeometry(auxMAT);
  }
  virtual void setGeometry(matrixE & newgeom) {}
  void setThresholdZero(double _thresholdzero) {thresholdzero = _thresholdzero;}
  void setRestart(int _restart) {
    //when changing molecule or initializing, restart = 0 (or -1)
    //if restart set to -1, then never to be considered
    //if restart == 0, by solving SCF restart set to 1
    //if restart == 1, then density reused
    restart = _restart;
  }
  void setAtomList(std::vector<size_t> & newatm) {
    atoms = newatm;
    Natoms = atoms.size();
  }
  //needed for MD
  virtual void setElectronTemp(double _etemp) {}
  virtual void setThresholdDensity(double new_thresh) {}
  virtual void setThresholdEnergy(double new_thresh) {}
  virtual void MayerBondOrder(matrixE & MBO, double distance_cutoff = 100.0) {}
  virtual void AtomicPolarizabilities(std::vector<double> & polarizabilities, const std::vector<double> & AtmCharge) {}
  virtual void TotalDispersion(double & mc6, double & mc8, const std::vector<double> & AtmCharge) {}
  virtual void TotalPolarizability(double & mpol, const std::vector<double> & AtmCharge) {}
  virtual double getHeatFormation() {return 0.0;}
  virtual double IonizationPotential(bool Koopman = true) {return 0.0;}
  virtual double ElectronAffinity() {return 0.0;}
  virtual void HSABdata(double & chi, double & eta) {}
  virtual double ReactivityIndices(matrixE & RxData, bool simple = false) {return 0.0;}
  //Calculation
  virtual void Calculate(int _print = 1, size_t maxiter = 200, bool _DIIS = true, bool _RCA = true, double d2threshold = 1600.0) {}
  //derivatives
  virtual void AnalyticalGrad(matrixE & gen, double tolerance = 1.0e-8) {}
  void NumericalGrad(matrixE & gen, double step = 1.0e-4) {
    //numerical gradients with central differences
    matrixE initialGeom = mol.Geometry();            //saving the initial geometry
    matrixE stepGeom = mol.Geometry();               //matrix containing the actual step
    double Eplus;
    double Eminus;
    double invstep = 0.5/step;
    size_t counter = 1;
    gen.resize(3*Natoms,1);
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      for (size_t idcoord = 0; idcoord < 3; ++idcoord,++counter) {
        //do geometry + step
        stepGeom(idatm + 1,idcoord + 1) = initialGeom(idatm + 1,idcoord + 1) + step;
        this->setGeometry(stepGeom);
        if (restart != -1) {restart = 0;}
        this->Calculate(0);
        Eplus = this->getEnergy(1);
        //do geometry - step
        stepGeom(idatm + 1,idcoord + 1) = initialGeom(idatm + 1,idcoord + 1) - step;
        this->setGeometry(stepGeom);
        if (restart != -1) {restart = 0;}
        this->Calculate(0);
        Eminus = this->getEnergy(1);
        //reset geometry
        stepGeom(idatm + 1,idcoord + 1) = initialGeom(idatm + 1,idcoord + 1);
        //now gradient
        gen(counter,1) = (Eplus - Eminus)*invstep;
      }
    }
  }
  void gEnergy(matrixE & gen, int type = 1, bool project = 0, double tolerance = 1.0e-8) {
    //function getting the gradients of the energy
    //all final derivatives come for nuclear coordinates in Angstrom!
    if (type == 0) {NumericalGrad(gen,tolerance);}                  //numerical derivatives
    else if (type == 1) {AnalyticalGrad(gen,tolerance);}
    if ((project)&&(BWilson.rows()*BWilson.cols() > 1)) {
      //if project and we have Bmatrix, then use it
      auxMAT = gen;
      gen = BWilson.trans()*auxMAT;
    }
  }
  void ProjectGradients(matrixE & gen) {
    //function that applies the projection to internal degrees of freedom only to some gradients available from somewhere else
    if (BWilson.rows()*BWilson.cols() > 1) {
      //if we have Bmatrix, then use it
      auxMAT = gen;
      gen = BWilson.trans()*auxMAT;
    }
  }
  virtual void AnalyticalHess(matrixE & hessian, int algorithm = 2, double tolerance = 1.0e-7, int maxiter = 50) {}
  void NumericalHessian(matrixE & hessian, double step = 1.0e-3, int gradtype = 1, bool setrestart1 = false) {
    //numerical Hessian with central differences
    matrixE initialGeom = mol.Geometry();            //saving the initial geometry
    matrixE stepGeom = mol.Geometry();               //matrix containing the actual step
    matrixE gplus;
    matrixE gminus;
    double invstep = 0.5/step;
    size_t icounter = 1;
    size_t jcounter = 1;
    hessian.resize(3*Natoms,3*Natoms);
    hessian.zero();
    //initial calculation to set restart to 1
    this->Calculate(0);
    if (setrestart1) {restart = 1;}
    else {restart = -1;}
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      for (size_t idcoord = 0; idcoord < 3; ++idcoord,++icounter) {
        //do geometry + step
        stepGeom(idatm + 1,idcoord + 1) = initialGeom(idatm + 1,idcoord + 1) + step;
        this->setGeometry(stepGeom);
        this->Calculate(0);
        this->gEnergy(gplus,gradtype,0);
        //do geometry - step
        stepGeom(idatm + 1,idcoord + 1) = initialGeom(idatm + 1,idcoord + 1) - step;
        this->setGeometry(stepGeom);
        this->Calculate(0);
        this->gEnergy(gminus,gradtype,0);
        //reset geometry
        stepGeom(idatm + 1,idcoord + 1) = initialGeom(idatm + 1,idcoord + 1);
        //now Hessian
        jcounter = 1;
        for (size_t idbtm = 0; idbtm < Natoms; ++idbtm) {
          for (size_t jdcoord = 0; jdcoord < 3; ++jdcoord,++jcounter) {
            hessian(icounter,jcounter) += (gplus(jcounter,1) - gminus(jcounter,1))*invstep;
          }
        }
      }
    }
    for (size_t idatm = 0; idatm < 3*Natoms; ++idatm) {
      for (size_t idbtm = 0; idbtm < idatm; ++idbtm) {
        invstep = 0.5*(hessian(idatm + 1,idbtm + 1) + hessian(idbtm + 1,idatm + 1));
        hessian(idatm + 1,idbtm + 1) = invstep;
        hessian(idbtm + 1,idatm + 1) = invstep;
      }
    }
  }
  virtual int AvailableGradients() {return 0;}
  void hEnergy(matrixE & hessian, int algorithm = 0, bool project = 0, bool setrestart1 = false, double damp = 0.01, double tolerance = 1.0e-7, int maxiter = 50) {
    //function calculating the hessian of the energy
    //all final derivatives come for nuclear coordinates in Angstrom!
    //algorithm defines the algorithm used: -1 -> identity
    //                                      0  -> numerical
    //                                      1  -> Head-Zerner V1 approximation
    //                                      2  -> full analytical
    //                                      3  -> Schlegel's approximate force field
    //                                      4  -> Lindh's approximate force field
    //                                      5  -> estimate Hessian using BFGS algorithm
    if (algorithm == 0) {NumericalHessian(hessian,1.0e-3,AvailableGradients(),setrestart1);}
    else if (algorithm == -1) {hessian = Identity(3*Natoms);}
    //THESE BELOW ARE BUGGY, DO NOT USE
    //else if ((algorithm == 1)||(algorithm == 2)) {AnalyticalHess(hessian,algorithm,tolerance,maxiter);}
    else if (algorithm == 3) {SchlegelApproxHessian(hessian,mol);}
    else if (algorithm == 4) {LindhApproxHessian(hessian,mol,2,1.0e-14,70.0,1.0e-10);}
    if (project) {
      auxMAT = mol.Geometry();
      ProjectForceConstants(hessian,auxMAT);
      auxMAT = mol.Geometry();
      WilsonBmatrix(hessian,auxMAT,damp);              //this modifies the Hessian as well -> diagonal matrix
    }
  }
  virtual void checkAtoms() {}
  bool TransitionMetal(size_t atomicnr) {
    //function to determine whether the atom is transition metal
    return (((atomicnr > 20)&&(atomicnr < 30))||((atomicnr > 38)&&(atomicnr < 48))||((atomicnr > 56)&&(atomicnr < 80)));
  }
  bool Metal(size_t atomicnr) {
    //is the element a metal?
    return MetalHead(atomicnr);
  }
  //charges
  virtual std::vector<double> CM1charge() {
    //just placeholder
    std::vector<double> cm1;
    return cm1;
  }
  virtual std::vector<double> CM2charge() {
    //just placeholder
    std::vector<double> cm2;
    return cm2;
  }
  virtual std::vector<double> CM3charge() {
    //just placeholder
    std::vector<double> cm3;
    return cm3;
  }
  virtual std::vector<double> CM5charge(int type = 0) {
    //just placeholder
    std::vector<double> cm5;
    return cm5;
  }
  void GoedeckerCharges(std::vector<double> & charges, std::vector<size_t> & atoms, matrixE & Geometry, double totalcharge) {
    //Ghasemi, Hofstetter, Saha, Goedecker, Phys. Rev. B, 92, 045131, 2015
    double temp;
    double RAB;
    double gammaAB;
    double sqrt2dpi = sqrt(2.0/pi);
    charges.resize(Natoms + 1);
    auxMAT3.resize(Natoms + 1,Natoms + 1);
    //get coordination numbers
    std::vector<double> ncoord;
    NCoord(ncoord,atoms,Geometry,ERFcount,false,1600.0,true);
    //get required quantities
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      temp = CNFACGFN0(atoms[idatm])/(sqrt(ncoord[idatm]) + 1.0e-14);
      charges[idatm] = -XIGFN0(atoms[idatm]) + temp*ncoord[idatm];                      //{X}_{i} = -{EN}_{i} + {k}_{i}\sqrt{{CN}_{i}}
      temp = ALPGFN0(atoms[idatm]);                                                     //alphas
      ncoord[idatm] = temp*temp;
      //    {J}_{ii}                + \frac{2{\gamma}_{ii}}{\sqrt{\pi}}
      RAB = GAMMAGFN0(atoms[idatm]) + sqrt2dpi/temp;
      auxMAT3(idatm + 1,idatm + 1) = RAB;                                          //def. Aii below equation 3
      rAB[0] = Geometry(idatm + 1,1);
      rAB[1] = Geometry(idatm + 1,2);
      rAB[2] = Geometry(idatm + 1,3);
      for (size_t idbtm = 0; idbtm < idatm; ++idbtm) {
        temp = 0.0;
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          RAB = rAB[idcoord] - Geometry(idbtm + 1,idcoord + 1);
          temp += RAB*RAB;
        }
        RAB = sqrt(temp)*dist_Angstrom2aum1;
        gammaAB = 1.0/sqrt(ncoord[idatm] + ncoord[idbtm]);                              //{\gamma}_{AB} = \frac{1}{\sqrt{{\alpha}^{2}_{A} + {\alpha}^{2}_{B}}}, as defined in text
        temp = erf(gammaAB*RAB)/RAB;                                                    //def. Aij below equation 3
        auxMAT3(idatm + 1,idbtm + 1) = temp;
        auxMAT3(idbtm + 1,idatm + 1) = temp;
      }
      auxMAT3(Natoms + 1,idatm + 1) = 1.0;                                         //equation 5
      auxMAT3(idatm + 1,Natoms + 1) = 1.0;                                         //equation 5
    }
    auxMAT3(Natoms + 1,Natoms + 1) = 0.0;                                          //equation 5
    charges[Natoms] = totalcharge;
    //solve to get partial charges
    Solve_Ax_eq_b(auxMAT3,charges);                                                        //equation 4
    charges.resize(Natoms);
  }
  std::vector<double> getCharges(std::string mode = "0", int type = 0) {
    //this is a more complete alternative to getQAtoms()
    //different modes for distinct results
    if ((mode == "Mulliken")||(mode == "mulliken")||(mode == "MULLIKEN")) {MullikenCharges(QAtoms);}
    else if ((mode == "Loewdin")||(mode == "loewdin")||(mode == "LOEWDIN")) {LowdinCharges();}
    else if ((mode == "CM1")||(mode == "cm1")) {QAtoms = CM1charge();}
    else if ((mode == "CM2")||(mode == "cm2")) {QAtoms = CM2charge();}
    else if ((mode == "CM3")||(mode == "cm3")) {QAtoms = CM3charge();}
    else if ((mode == "CM5")||(mode == "cm5")) {QAtoms = CM5charge(type);}
    else if ((mode == "Goedecker")||(mode == "goedecker")||(mode == "GOEDECKER")) {
      matrixE geom = mol.Geometry();
      GoedeckerCharges(QAtoms,atoms,geom,mol.Charge());
    }
    return QAtoms;
  }
  double DipoleMoment(std::string mode = "0", int type = 0) {
    std::vector<double> qk;
    if ((mode == "CM1")||(mode == "cm1")) {qk = CM1charge();}
    else if ((mode == "CM2")||(mode == "cm2")) {qk = CM2charge();}
    else if ((mode == "CM3")||(mode == "cm3")) {qk = CM3charge();}
    else if ((mode == "CM5")||(mode == "cm5")) {qk = CM5charge(type);}
    else {throw("ERROR: Method.hpp: Method: DipoleMoment(): charge model not identified");}
    matrixE geom = mol.Geometry();
    double xkqk = 0.0;
    double ykqk = 0.0;
    double zkqk = 0.0;
    for (size_t iatm = 0; iatm < qk.size(); ++iatm) {
      xkqk += qk[iatm]*geom(iatm + 1,1);
      ykqk += qk[iatm]*geom(iatm + 1,2);
      zkqk += qk[iatm]*geom(iatm + 1,3);
    }
    return sqrt(xkqk*xkqk + ykqk*ykqk + zkqk*zkqk)/D2eA;
  }
  //population analysis
  virtual void LowdinCharges() {}
  virtual void LMcharges(std::vector<double> & charges, matrixE & ds, std::vector<size_t> & AOs, std::vector<double> & Ccharge) {}
  virtual void Mulliken(matrixE & mulliken) {}
  void getMullikenCharges() {
    //set internal charges to Mulliken charges
    MullikenCharges(QAtoms);
  }
  virtual void MullikenCharges(std::vector<double> & atomPop) {}
  //module for vibrational frequencies
  std::vector<double> CalcVibrFrequencies(bool restartcalcs = true) {
    //wrapper to get vibrational frequencies without headaches
    matrixE hessian;
    this->hEnergy(hessian,0,0,restartcalcs,0.0,1.0e-7,50);
    std::vector<double> Vfreq = this->Vibrations(hessian,1,1,true);
    return Vfreq;
  }
  matrixE PerturbGeometry(double fraction, int nevec, bool restartcalcs = false) {
    //function that perturbs the geometry using 100*fraction % of eigenvector number nevec
    //if nevec < 0 then looks for TS frequency
    matrixE hessian;
    this->hEnergy(hessian,0,0,restartcalcs,0.0,1.0e-7,50);
    std::vector<double> Vfreq = this->Vibrations(hessian,1,1,true);
    int eigenvector_index = nevec - 1;
    size_t nmodes = Vfreq.size();
    if (nevec < 0) {
      for (size_t ideig = 0; ideig < nmodes; ++ideig) {
        if (Vfreq[ideig] < 0.0) {
          eigenvector_index = ideig;
          break;
        }
      }
    }
    if (eigenvector_index < 0) {eigenvector_index = 0;}
    matrixE newgeom = mol.Geometry();
    nmodes = 0;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      for (size_t idcoord = 0; idcoord < 3; ++idcoord, ++nmodes) {
        newgeom(idAtm + 1,idcoord + 1) += fraction*hessian(nmodes + 1,eigenvector_index + 1);
      }
    }
    return newgeom;
  }
  std::vector<double> Vibrations(matrixE & hessian, int unitout = 1, int unitin = 1, bool zerofirstvibr = true, double negativefreqthresh = 0.1) {
    //function that calculates vibrational frequencies and the vibrational modes
    //unitout determines which units the result comes out with: 0 -> atomic units;
    //                                                          1 -> cm-1
    //unitin defines the units in which the hessian comes in:   0 -> Bohr
    //                                                          1 -> Angstroem
    //zerovibr determines whether to zero the first f vibrational modes
    //the vibrational modes are returned in place of the Hessian, meaning the latter is overwritten
    //get the inverse square root of the atomic masses
    matrixE geom = mol.Geometry();
    std::vector<double> sqrtmim1(Natoms,0.0);
    double unitconv = 1.0/(au2kg*NA);
    double factor;
    double factor2 = 1.0;
    if (unitin == 1) {factor2 = dist_Angstrom2au*dist_Angstrom2au;}
    int AtmA;
    int AtmB;
    size_t zerovibr = 6;
    size_t zerovibstart = 0;
    if (mol.Linear()) {zerovibr = 5;}
    if (!zerofirstvibr) {zerovibr = 0;}
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      sqrtmim1[idAtm] = 1.0/sqrt(Weight(atoms[idAtm])*unitconv);
    }
    //project out the translations and rotations
    if (!mol.Linear()) {ProjectForceConstants(hessian,geom);}
    //mass-weight the Hessian
    for (size_t idrow = 0; idrow < 3*Natoms; ++idrow) {
      AtmA = idrow/3;
      hessian(idrow + 1,idrow + 1) *= sqrtmim1[AtmA]*sqrtmim1[AtmA]*factor2;
      for (size_t idcol = 0; idcol < idrow; ++idcol) {
        AtmB = idcol/3;
        factor = sqrtmim1[AtmA]*sqrtmim1[AtmB]*factor2;
        hessian(idrow + 1,idcol + 1) *= factor;
        hessian(idcol + 1,idrow + 1) *= factor;
      }
    }
    std::vector<double> VibrFreq = MatDiag(hessian);
    //get vibrational frequencies
    unitconv = (au2cmm1 - 1)*(unitout == 1) + 1;             //unit conversion to cm-1
    for (size_t idfreq = 0; idfreq < 3*Natoms; ++idfreq) {
      factor = sqrt(fabs(VibrFreq[idfreq]));
      if (VibrFreq[idfreq] < 0.0) {factor *= -1.0;}
      VibrFreq[idfreq] = factor*unitconv;
      if (VibrFreq[idfreq] < -negativefreqthresh) {++zerovibstart;}
    }
    zerovibr += zerovibstart;
    for (size_t idx = zerovibstart; idx < zerovibr; ++idx) {
      VibrFreq[idx] = 0.0;
    }
    return VibrFreq;
  }
  void WilsonBmatrix(matrixE & hessian, matrixE & geom, double damp = 0.01) {
    //function that calculates a Wilson B matrix used in geometry optimization
    int count = 1;
    int counteaux = 0;
    int counttrartomodes = 0;           //counter to see how many translational and rotational modes we have found
    double Rnormal;
    double Rdistorted;
    double Rtotal;
    double aux;
    BWilson = hessian*dist_Angstrom2au*dist_Angstrom2au;
    matrixE distgeom(Natoms,3);
    std::vector<double> eigval = MatDiag(BWilson);
    std::vector<int> eindex(3*Natoms,0);
    std::vector<double> eaux(3*Natoms,-1.0);
    //check for translations and rotations
    for (size_t ideig = 0; ideig < 3*Natoms; ++ideig) {
      if (fabs(eigval[ideig]) > 0.05) {continue;}
      //distort molecule along the vibrational modi
      distgeom = geom*dist_Angstrom2aum1;
      count = 1;
      Rtotal = 0.0;              //this is used to compare the two geometries
      for (size_t idatm = 1; idatm < Natoms + 1; ++idatm) {
        for (size_t idcoord = 1; idcoord < 4; ++idcoord, ++count) {
          distgeom(idatm,idcoord) += BWilson(count,ideig + 1);
        }
        for (size_t idatm2 = 1; idatm2 < idatm; ++idatm2) {
          Rnormal = 0.0;
          Rdistorted = 0.0;
          for (size_t idc = 1; idc < 4; ++idc) {
            aux = geom(idatm,idc) - geom(idatm2,idc);
            Rnormal += aux*aux;
            aux = distgeom(idatm,idc) - distgeom(idatm2,idc);
            Rdistorted += aux*aux;
          }
          aux = sqrt(Rnormal)*dist_Angstrom2aum1 - sqrt(Rdistorted);
          Rtotal += aux*aux;
        }
      }
      eaux[counteaux] = sqrt(Rtotal/Natoms)*fabs(eigval[ideig]);
      eindex[counteaux] = ideig;
      ++counteaux;
      ++counttrartomodes;
    }
    if ((counttrartomodes != 5)&&(counttrartomodes != 6)) {
      //reorder
      for (size_t idstep = 1; idstep < counttrartomodes; ++idstep) {
        for (size_t idswap = 0; idswap < counttrartomodes - idstep; ++idswap) {
          if (eaux[idswap] > eaux[idswap + 1]) {
            Rnormal = eaux[idswap + 1];
            eaux[idswap + 1] = eaux[idswap];
            eaux[idswap] = Rnormal;
            counteaux = eindex[idswap];
            eindex[idswap] = eindex[idswap + 1];
            eindex[idswap + 1] = counteaux;
          }
        }
      }
    }
    for (size_t idx = 0; idx < 5; ++idx) {
      eigval[eindex[idx]] = 0.0;
    }
    aux = eigval[eindex[5]];                   //minimum eigenvalue
    if (!mol.Linear()) {
      eigval[eindex[5]] = 0.0;
      aux = eigval[eindex[6]];
    }
    if (fabs(aux) < 1.0e-10) {aux = 1.0e-10;}
    //shift so that smallest eigenvalue has at least a minimum value
    Rnormal = fmax(damp - aux,0.0);
    for (size_t ideig = 0; ideig < 3*Natoms; ++ideig) {
      aux = (fabs(eigval[ideig]) > 1.0e-11)*(eigval[ideig] + Rnormal);
      eigval[ideig] = aux;
    }
    //build the B matrix and the new Hessian; in this case the Hessian is diagonal
    counttrartomodes = 3*Natoms - mol.Linear()*5 - (!mol.Linear())*6;
    distgeom = BWilson;
    BWilson.resize(3*Natoms,counttrartomodes);
    hessian.resize(counttrartomodes,counttrartomodes);
    hessian.zero();
    count = 0;
    for (size_t idr = 0; idr < 3*Natoms; ++idr) {
      if (fabs(eigval[idr]) < 1.0e-11) {continue;}
      ++count;
      hessian(count,count) = fmin(fmax(eigval[idr],damp),5.0)*dist_Angstrom2aum1*dist_Angstrom2aum1;       //damp the Hessian on both ends to ensure a clean optimization
      for (size_t idx = 1; idx < 3*Natoms + 1; ++idx) {
        BWilson(idx,count) = distgeom(idx,idr + 1);
      }
    }
  }
};

#endif //_METHOD_METHODS_
