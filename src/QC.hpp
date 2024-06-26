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

#ifndef _QC_METHODS_
#define _QC_METHODS_
#include "BSet.hpp"
#include "Density.hpp"
#include "Method.hpp"
#include "solvation/DiElectricConsts.hpp"
#include "solvation/SolvationGrid.hpp"

//descrition:
//base class for quantum chemical methods
//quantities without prefix "b" are always used in the closed-shell case; for open-shell systems, these become containers for alpha electrons
//quantities with prefix "b" refer to beta electrons, and are used strictly in open-shell calculations

class QCbasis: public Method {
protected:
  std::vector<int> CoreCharge;
  std::vector<double> occupancy;                      //occupancy vector for orbitals
  std::vector<double> boccupancy;
  std::vector<double> EMOs;                           //MO energy vectors
  std::vector<double> bEMOs;
  std::vector<double> solvForce;
  std::vector<size_t> AOs;
  matrixE CMO;                                        //orbital coefficient matrix
  matrixE bCMO;
  matrixE Hcore;                                      //core Hamiltonian
  matrixE Fock;                                       //this is a hard one to guess...
  matrixE bFock;
  matrixE VAB;
  matrixE sao;
  matrixE dens;
  matrixE bdens;
  matrixE auxMAT2;
  matrixE DensAUX;
  matrixE BornMat;                                    //Born matrix, as required by generalized Born models
  std::string openclosed;
  std::string shell;                                  //control string to determine whether system is open- or closed-shell
  BSet basis;                                         //basis set object
  std::string solvationmodel;                         //the type of solvation model used
                                                      //ALPB -> default
                                                      //GB
  std::string solvent;
  size_t NAOs;                                        //total number of atomic orbitals used
  int CoulombSolv;                                    //the effective Coulomb operator in generalized Born solvation
                                                      //0 -> Still
                                                      //1 -> P16 (default)
  int ngrid;                                          //number of grid points for integration in solvation model
  int startpseudodiag;                                //iteration number when to start using pseudo diagonalization
  bool convEn;                                        //energy convergence control
  bool convDen;                                       //density convergence control
  bool usesolvation;                                  //determines whether to use solvation
  bool uselevelshift;                                 //determines whether to use level-shift in the SCF
  double enuc;                                        //nuclear contribution to the energy
  double esolv;
  double epsilon_in;                                  //{\varepsilon}_{in} parameter required by generalized Born models
  double saltmolarity;                                //if solvent has salt, its molarity (mol/L)
  double saltradius;                                  //if solvent has salt, its radius (Bohr)
  double temperature;                                 //for solvation models
  double SASA;                                        //energy contribution from solvent available surface area
  double Telectron;                                   //electronic tempeperature for fractional occupation numbers in DFT(B) methods
  double epsilonS;                                    //the minimum value for HOMO-LUMO gap in the level-shift approach
  double thresh_dens;
  double thresh_en;
public:
  double thresholdzero;                               //whenever needed, this is the numerical zero
  QCbasis(BSet & _bset, Molecule & _mol): Method(_mol), basis(_bset) {
    enuc = 0.0;
    NAOs = basis.NAO();
    CoreCharge.resize(Natoms);
    usesolvation = false;
    solvationmodel = "ALPB";
    temperature = 298.15;
    ngrid = 230;
    epsilon_in = 1.0;
    saltmolarity = 0.0;
    saltradius = 0.0;
    SASA = 0.0;
    solvent = "0";
    esolv = 0.0;
    solvForce.resize(Natoms);
    CoulombSolv = 1;
    Telectron = 0.0;
    uselevelshift = true;
    epsilonS = 0.1;
    thresh_dens = 1e-7;
    thresh_en = 1e-7;
    startpseudodiag = 1;
  }
  ~QCbasis() {}
  //getters
  virtual std::string Type() {return "QCbasis";}
  virtual std::string TypeOverlap() {return "0";}
  std::string ShellType() {return openclosed;}
  std::string MechanicsType() {return "quantum";}
  matrixE MOs() {return CMO;}
  matrixE bMOs() {return bCMO;}
  matrixE getVAB() {return VAB;}
  std::vector<double> EnMOs() {return EMOs;}
  std::vector<double> bEnMOs() {return bEMOs;}
  std::vector<double> occ() {return occupancy;}
  std::vector<double> bocc() {return boccupancy;}
  size_t NAO() {return NAOs;}
  std::vector<size_t> getAtomicOrbitals() {
    AOs = basis.AtomNAOs(atoms);
    return AOs;
  }
  std::string Shell() {return shell;}
  std::vector<int> CoreCharges() {return CoreCharge;}
  double ThresholdDensity() {return thresh_dens;}
  double ThresholdEnergy() {return thresh_en;}
  double getEnuclear() {return enuc;}
  virtual double Energy() {return 0.0;}
  virtual double getHeatFormation() {return 0.0;}
  matrixE getSAO() {return sao;}
  matrixE getHcore() {return Hcore;}
  matrixE getFock() {return Fock;}
  matrixE getFockb() {return bFock;}
  void getDens(matrixE & dens) {
    if ((dens.rows()*dens.cols() != this->NAOs*this->NAOs)||(dens.rows() != dens.cols())) {dens.resize(this->NAOs,this->NAOs);}
    dens.zero();
    Dens(dens,this->occupancy,this->CMO,DensAUX);
  }
  void getbDens(matrixE & dens) {
    if ((dens.rows()*dens.cols() != this->NAOs*this->NAOs)||(dens.rows() != dens.cols())) {dens.resize(this->NAOs,this->NAOs);}
    dens.zero();
    Dens(dens,this->boccupancy,this->bCMO,DensAUX);
  }
  void fetchDens(matrixE & dens_) {dens_ = dens;}
  void fetchbDens(matrixE & bdens_) {bdens_ = bdens;}
  int StartPseudoDiagonalization() {return startpseudodiag;}
  bool UseSolvation() {return usesolvation;}
  bool UseLevelShift() {return uselevelshift;}
  double EpsilonS() {return epsilonS;}
  std::string SolvationModel() {return solvationmodel;}
  std::string Solvent() {return solvent;}
  double SaltMolarity() {return saltmolarity;}
  double SaltRadius() {return saltradius;}
  double Temperature() {return temperature;}
  int NGridPoints() {return ngrid;}
  double EpsilonIn() {return epsilon_in;}
  double ElectronTemp() {return Telectron;}
  matrixE BornMatrix() {return BornMat;}
  double SolventAccessibleSurfaceArea() {return SASA;}
  std::vector<double> SolvationForce() {return solvForce;}
  int EffectiveCoulombOperator() {return CoulombSolv;}
  int HOMOposition() {
    //function to find the HOMO
    int szocc = occupancy.size();
    int homo = 0;
    for (int idx = 0; idx < szocc; ++idx) {
      if (int(occupancy[idx]) == 0) {
        homo = idx - 1;
        break;
      }
    }
    if ((homo == 0)&&(szocc != 0)) {homo = szocc - 1;}
    return homo;
  }
  double getOrbital(int orbnum, std::vector<double> & ORB) {
    //function to extract an orbital
    ORB.resize(NAOs);
    for (size_t idx = 0; idx < NAOs; ++idx) {
      ORB[idx] = CMO(idx + 1,orbnum);
    }
    return EMOs[orbnum - 1];
  }
  double getHOMO(std::vector<double> & HOMO) {
    //function to extract the HOMO
    int homo = HOMOposition();
    HOMO.resize(NAOs);
    for (size_t idx = 0; idx < NAOs; ++idx) {
      HOMO[idx] = CMO(idx + 1,homo + 1);
    }
    return EMOs[homo];
  }
  double getLUMO(std::vector<double> & LUMO) {
    //function to extract the LUMO
    int lumo = HOMOposition() + 1;
    LUMO.resize(NAOs);
    for (size_t idx = 0; idx < NAOs; ++idx) {
      LUMO[idx] = CMO(idx + 1,lumo + 1);
    }
    return EMOs[lumo];
  }
  virtual matrixE getQShells() {return Hcore;}
  virtual matrixE getmuAtoms() {return Hcore;}
  virtual matrixE getthetaAtoms() {return Hcore;}
  virtual std::vector<int> AtomicShells() {return CoreCharge;}
  //setters
  void setMolecule(Molecule & _component) {
    mol = _component;
    atoms = mol.Atoms();
    Natoms = atoms.size();
    solvForce.resize(Natoms);
    if (restart != -1) {restart = 0;}
    basis.setMolecule(mol);
    NAOs = basis.NAO();
  }
  void setGeometry(matrixE & newgeom) {
    mol.setGeometry(newgeom);
    basis.setMolecule(mol);
  }
  void setNAOs(size_t newNAOs) {NAOs = newNAOs;}
  void setAOBasis(const BSet & _bset) {
    basis = _bset;
    NAOs = basis.NAO();
    if (restart != -1) {restart = 0;}
  }
  void setMOs(const matrixE & newMO) {CMO = newMO;}
  void setbMOs(const matrixE & newMO) {bCMO = newMO;}
  void setOcc(const std::vector<double> & newocc) {occupancy = newocc;}
  void setbOcc(const std::vector<double> & newocc) {boccupancy = newocc;}
  void setElectronTemp(double _etemp) {Telectron = _etemp;}
  void setThresholdDensity(double new_thresh) {thresh_dens = new_thresh;}
  void setThresholdEnergy(double new_thresh) {thresh_en = new_thresh;}
  void setThresholdZero(double _thresholdzero) {thresholdzero = _thresholdzero;}
  void setEMOs(std::vector<double> _EMOs) {EMOs = _EMOs;}
  void setbEMOs(std::vector<double> _bEMOs) {bEMOs = _bEMOs;}
  void setSAO(matrixE & _sao) {sao = _sao;}
  void setDens(matrixE & _dens) {dens = _dens;}
  void setbDens(matrixE & _bdens) {bdens = _bdens;}
  void setHcore(matrixE & hcore) {Hcore = hcore;}
  void setStartPseudoDiagonalization(int newiterval) {startpseudodiag = newiterval - 1;}        //C++ has to be special and start on zero
  void setUseSolvation(bool usesolvation_) {usesolvation = usesolvation_;}
  void setUseLevelShift(bool levelshift) {uselevelshift = levelshift;}
  void setEpsilonS(double eps) {epsilonS = eps;}
  void setSolvationModel(std::string solvationmodel_) {solvationmodel = solvationmodel_;}
  void setSolvent(std::string solvent_, double molarity_ = -1.0, double radius_ = 0.0) {
    //even though the radius of salt ions is in Angstroem, it is converted to Bohr
    solvent = solvent_;
    usesolvation = true;
    if (molarity_ > 0.0) {
      saltmolarity = molarity_;
      saltradius = radius_*dist_Angstrom2aum1;
    }
  }
  void setSaltMolarity(double molarity_) {if (molarity_ > 0.0) {saltmolarity = molarity_;}}
  void setSaltRadius(double saltradius_) {
    //even though the radius of salt ions is in Angstroem, it is converted to Bohr
    saltradius = saltradius_*dist_Angstrom2aum1;
  }
  void setTemperature(double temperature_) {temperature = temperature_;}
  void setNGridPoints(int ngr) {ngrid = ngr;}
  void setEpsilonIn(double epsilon_in_) {epsilon_in = epsilon_in_;}
  void setEffectiveCoulombOperator(int CoulombSolv_) {CoulombSolv = CoulombSolv_;}
  virtual void setQShells(matrixE & _QShells) {}
  virtual void setmuAtoms(matrixE & _muAtoms) {}
  virtual void setthetaAtoms(matrixE & _thetaAtoms) {}
  //Calculation
  virtual void calcFock(matrixE & SAO) {}
  virtual void transformFock(matrixE & XAO) {}
  virtual void transformCMO(matrixE & XAO) {}
  //calculation of orbitals
  void calcMOs() {
    //function that calculates the MOs
    CMO = Fock;
    EMOs = MatDiag(CMO);
    if ((shell == "open")&&(bFock.rows()*bFock.cols() > 1)) {
      bCMO = bFock;
      bEMOs = MatDiag(bCMO);
    }
  }
  void ElectronicDensity(std::string density_filename = "density") {
    //function that calculates the electronic density of the molecule
    matrixE geometry = mol.Geometry();
    matrixE NQN(Natoms,3);                                            //principal quantum numbers
    matrixE LQN(Natoms,3);                                            //azimuthal quantum numbers
    basis.getQuantumNumbers(NQN,LQN,atoms);
    std::vector<std::vector<double> > AOexp = basis.GetExponents();   //AO exponents
    std::vector<double> step(3,0.1);
    ElectronicDensityCalc(density_filename,dens,geometry,atoms,AOexp,NQN,LQN,step,0.001,1.0e-4,3.0);
  }
  void PseudoDiagonalization(matrixE & MOcoeff, matrixE & FockM, std::vector<double> & MOen, int nocc, double thresholdzero = 1.0e-8) {
    //function applying the pseudo-diagonalization algorithm with the mopac7 modifications
    //J. J. P. Stewart, P. Csaszar, P. Pulay, J. Comp. Chem., 3(2), 227, 1982
    //get occupied block of MO coefficients
    int ncolfetch = int(NAOs) - nocc;
    MOcoeff.getBlock(auxMAT,1,1,NAOs,nocc);
    auxMAT2 = auxMAT.trans()*FockM;
    //now the virtual block
    MOcoeff.getBlock(auxMAT,1,nocc + 1,NAOs,ncolfetch);
    auxMAT3 = auxMAT2*auxMAT;
    double fia;
    double ei;
    double ea;
    double dia;
    double ss;
    double cc;
    double tt;
    double factor;
    double vi;
    double va;
    //now rotate to get the new MO coefficients
    if ((auxMAT.rows() != NAOs)||(auxMAT.cols() != NAOs)) {auxMAT.resize(NAOs,NAOs);}
    auxMAT = MOcoeff;
    for (int idocc = 0; idocc < nocc; ++idocc) {
      for (int idvirt = nocc; idvirt < NAOs; ++idvirt) {
        fia = auxMAT3(idocc + 1,idvirt - nocc + 1);
        if (fabs(fia) < thresholdzero) {continue;}
        ei = MOen[idocc];
        ea = MOen[idvirt];
        dia = ei - ea;
        factor = double(dia > 0) - double(dia < 0);
        cc = factor*sqrt(4.0*fia*fia + dia*dia);
        ss = sqrt(0.5*(1.0 + dia/cc));
        factor = double(fia > 0) - double(fia < 0);
        tt = -factor*sqrt(1.0 - ss*ss);
        for (int idrow = 0; idrow < NAOs; ++idrow) {
          vi = auxMAT(idrow + 1,idocc + 1);
          va = auxMAT(idrow + 1,idvirt + 1);
          MOcoeff(idrow + 1,idocc + 1) = ss*vi + tt*va;
          MOcoeff(idrow + 1,idvirt + 1) = ss*va - tt*vi;
        }
      }
    }
  }
  void ShiftFock(bool openshell = false, double thresholdzero = 1.0e-6) {
    //function applying level-shift to the Fock matrix
    if (uselevelshift) {
      int LUMO = occupancy.size();
      int HOMO = LUMO - 1;
      double factor = 0.5;
      double shift = 2.0;
      double HLgap = EMOs[LUMO] - EMOs[HOMO];
      if (openshell) {
        factor = 1.0;
        shift = 1.0;
      }
      ApplyLevelShift(HLgap,factor,shift,dens,Fock);
      if (openshell) {
        LUMO = boccupancy.size();
        HOMO = LUMO - 1;
        HLgap = bEMOs[LUMO] - bEMOs[HOMO];
        ApplyLevelShift(HLgap,factor,shift,bdens,bFock);
      }
    }
  }
  void ApplyLevelShift(double HLgap, double factor, double shift, matrixE & densM, matrixE & FockM, double thresholdzero = 1.0e-6) {
    //the actual level-shift
    double alphaS = fmax(factor*(epsilonS - HLgap),0.0);
    if (alphaS > thresholdzero) {
      auxMAT = densM*(-1.0);
      for (size_t idao = 0; idao < NAOs; ++idao) {
        auxMAT(idao + 1,idao + 1) += shift;
      }
      //now shift the Fock matrix
      FockM += alphaS*auxMAT;
    }
  }
  //occupancy
  virtual void Occ() {}
  //core Hamiltonian
  virtual void calcHcore(matrixE & sao) {}
  //other, somehow energy related
  virtual void PrepareEnergyCalculation(matrixE & geometry) {}
  void CalcOverlap(matrixE & PbarAlpha, matrixE & PbarBeta, bool openshell, int endmonomer, bool orthogonalbasis) {
    //this is the function that calculates the overlap between two monomer wavefunctions
    //orthogonalbasis determines whether to use NDDO (true) or something else (false)
    //this is because the AO basis in these methods is orthogonal for the monomers
    matrixE SAOsupermol = basis.SAO();
    std::vector<double> EvalMonomer(occupancy.size(),1.0);
    std::vector<double> EvalSupermol(boccupancy.size(),1.0);
    matrixE auxSAOmonomer = Identity(NAOs);
    matrixE auxSAOsupermol = Identity(NAOs);
    if (!orthogonalbasis) {
      matrixE SAOmonomer = SAOsupermol;
      matrixE auxcalc;
      for (size_t idorb = 0; idorb < endmonomer; ++idorb) {
        for (size_t jdorb = endmonomer; jdorb < NAOs; ++jdorb) {
          SAOmonomer(jdorb + 1,idorb + 1) = 0.0;
          SAOmonomer(idorb + 1,jdorb + 1) = 0.0;
        }
      }
      //get the transformation matrices
      EvalMonomer = MatDiag(SAOmonomer);
      EvalSupermol = MatDiag(SAOsupermol);
      for (size_t idao = 0; idao < NAOs; ++idao) {
        auxSAOmonomer(idao + 1,idao + 1) = sqrt(EvalMonomer[idao]);
        auxSAOsupermol(idao + 1,idao + 1) = 1.0/sqrt(EvalSupermol[idao]);
      }
      //get the matrix S1/2 for the monomer basis
      auxcalc = SAOmonomer*auxSAOmonomer;
      auxcalc *= SAOmonomer.trans();
      SAOmonomer = auxcalc;
      //get the matrix S-1/2 for the supermolecule basis
      auxcalc = SAOsupermol*auxSAOsupermol;
      auxcalc *= SAOsupermol.trans();
      SAOsupermol = auxcalc;
      //convert MOs to orthogonal basis
      auxcalc = SAOmonomer*CMO;
      //convert MOs to the supermolecule basis
      SAOmonomer = SAOsupermol*auxcalc;
      SAOsupermol = SAOmonomer.trans();
      generalizedDens(PbarAlpha,occupancy,SAOmonomer,SAOsupermol,1.0e-10);
    }
    else if (orthogonalbasis) {
      //modify AO basis for monomers
      for (size_t idAO = 0; idAO < endmonomer; ++idAO) {
        for (size_t idBO = 0; idBO < idAO; ++idBO) {
          SAOsupermol(idAO + 1,idBO + 1) = 0.0;
          SAOsupermol(idBO + 1,idAO + 1) = 0.0;
        }
        SAOsupermol(idAO + 1,idAO + 1) = 1.0;
      }
      for (size_t idAO = endmonomer; idAO < NAOs; ++idAO) {
        for (size_t idBO = endmonomer; idBO < idAO; ++idBO) {
          SAOsupermol(idAO + 1,idBO + 1) = 0.0;
          SAOsupermol(idBO + 1,idAO + 1) = 0.0;
        }
        SAOsupermol(idAO + 1,idAO + 1) = 1.0;
      }
      basis.setSAO(SAOsupermol);
      SAOsupermol = basis.S12();
      generalizedDens(auxSAOmonomer,EvalMonomer,SAOsupermol,CMO,1.0e-10);
      auxSAOsupermol = auxSAOmonomer.trans();
      generalizedDens(PbarAlpha,occupancy,auxSAOmonomer,auxSAOsupermol,1.0e-10);
      dens = PbarAlpha;
      if (openshell) {
        //now the same for beta spin
        generalizedDens(auxSAOmonomer,EvalSupermol,SAOsupermol,bCMO,1.0e-10);
        auxSAOsupermol = auxSAOmonomer.trans();
        generalizedDens(PbarBeta,occupancy,auxSAOmonomer,auxSAOsupermol,1.0e-10);
        bdens = PbarBeta;
      }
    }
  }
  void CleanMatrix(matrixE & matrix, int endmonomer) {
    for (size_t idAO = 0; idAO < endmonomer; ++idAO) {
      for (size_t idBO = endmonomer; idBO < NAOs; ++idBO) {
        matrix(idAO + 1,idBO + 1) = 0.0;
        matrix(idBO + 1,idAO + 1) = 0.0;
      }
    }
  }
  virtual void GFermi(double thresh = 1.0e-9, double expthresh = 50.0, double maxiter = 200) {}
  //other auxiliary functions
  virtual matrixE EnuclearMat() {
    matrixE enuc(1,1);
    return enuc;
  }
  virtual void Enuclear(matrixE & geometry, matrixE & Emat, double d2threshold = 1600.0) {}
  virtual void calcExchangeFock(matrixE & xFockA, matrixE & xFockB, matrixE & SAO) {
    xFockA.resize(NAOs,NAOs);
    xFockB.resize(NAOs,NAOs);
    for (size_t idAO = 0; idAO < NAOs; ++idAO) {
      for (size_t jdAO = 0; jdAO < idAO + 1; ++jdAO) {
        xFockA(idAO + 1,jdAO + 1) = 0.0;
        xFockA(jdAO + 1,idAO + 1) = 0.0;
        xFockB(idAO + 1,jdAO + 1) = 0.0;
        xFockB(jdAO + 1,idAO + 1) = 0.0;
      }
    }
  }
  virtual void EnergyDispersion(matrixE & geometry, matrixE & Edisp, double d2threshold = 1600.0) {}
  virtual void EnergyMatrix(matrixE & geometry, matrixE & Emat, bool redoCAMM = false, int cleanmat = -1) {}
  virtual void ConstrainedSCF(int endmonomerA, double NelecA, double NelecB, bool fromscratch, bool popconstrain, int _print = 1, size_t maxiter = 200, double thresholdVk = 1.0e-6) {}
  //population analysis
  void LowdinCharges() {
    //wrapper around LMcharges
    AOs = basis.AtomNAOs(atoms);
    std::vector<double> ccharge(Natoms);
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      ccharge[idatm] = double(CoreCharge[idatm]);
    }
    this->getDens(dens);
    if (this->Shell() == "open") {this->getbDens(bdens);}
    dens += bdens;
    LMcharges(QAtoms,dens,AOs,ccharge);
  }
  void LMcharges(std::vector<double> & charges, matrixE & ds, std::vector<size_t> & AOs, std::vector<double> & Ccharge) {
    //function returning the Loewdin/Mulliken charges for each atom; note that at NDDO level these are the same because S = I
    int icnt = 1;
    charges.resize(Natoms);
    for (size_t iatom = 0; iatom < Natoms; ++iatom) {
      charges[iatom] = double(atoms[iatom]) - Ccharge[iatom];            //this assumes that Z - Ncore > 0, as it should, for Z the nuclear charge and Ncore the number of pairs of core electrons
      for (size_t iaos = 0; iaos < AOs[iatom]; ++iaos, ++icnt) {
        charges[iatom] -= ds(icnt,icnt);
      }
    }
  }
  void Mulliken(matrixE & mulliken) {
    //function that gets the Mulliken matrix
    mulliken = dens;
    if ((shell == "open")&&(bdens.rows()*bdens.cols() > 1)) {mulliken += bdens;}
    //calculate the Mulliken matrix
    if (TypeOverlap() != "Orthogonal") {
      for (size_t irw = 0; irw < NAOs; ++irw) {
        for (size_t icl = irw; icl < NAOs; ++icl) {
          mulliken(irw + 1,icl + 1) *= sao(irw + 1,icl + 1);
          mulliken(icl + 1,irw + 1) *= sao(icl + 1,irw + 1);
        }
      }
    }
    else {
      for (size_t irw = 0; irw < NAOs; ++irw) {
        for (size_t icl = irw + 1; icl < NAOs; ++icl) {
          mulliken(irw + 1,icl + 1) = 0.0;
          mulliken(icl + 1,irw + 1) = 0.0;
        }
      }
    }
  }
  void MullikenCharges(std::vector<double> & atomPop) {
    //function returning a vector with atomic partial charges
    matrixE mulliken(1,1);
    Mulliken(mulliken);
    atomPop.resize(Natoms);
    size_t icnt = 0;
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      atomPop[idatm] = double(atoms[idatm]) - CoreCharge[idatm];
      for (size_t iAOrb = 0; iAOrb < AOs[idatm]; ++iAOrb) {
        for (size_t iao = 0; iao < NAOs; ++iao) {
          atomPop[idatm] -= mulliken(icnt + iAOrb + 1,iao + 1);
        }
      }
      icnt += AOs[idatm];
    }
  }
  //bond-order functions
  void MullikenBondOrder(matrixE & MBO, double distance_cutoff = 100.0*dist_Angstrom2au*dist_Angstrom2au) {
    //function with the Mulliken bond orders
    matrixE mulliken(1,1);
    Mulliken(mulliken);
    MBO.resize(Natoms,Natoms);
    size_t icntA = 0;
    size_t icntB = 0;
    double bondorder;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      icntB = 0;
      MBO(idAtm + 1,idAtm + 1) = 0.0;
      for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
        bondorder = 0.0;
        for (size_t idAO = 0; idAO < AOs[idAtm]; ++idAO) {
          for (size_t idBO = 0; idBO < AOs[idBtm]; ++idBO) {
            bondorder += mulliken(icntA + idAO + 1,icntB + idBO + 1)*mulliken(icntA + idAO + 1,icntB + idBO + 1);
          }
        }
        MBO(idAtm + 1,idBtm + 1) = bondorder;
        MBO(idBtm + 1,idAtm + 1) = bondorder;
        icntB += AOs[idBtm];
      }
      icntA += AOs[idAtm];
    }
  }
  void MayerBondOrder(matrixE & MBO, double distance_cutoff = 100.0*dist_Angstrom2au*dist_Angstrom2au) {
    //function returning Mayer bond orders
    //I. Mayer, Chem. Phys. Lett., 97, 270, 1983
    int icnt = 1;
    int icns = 1;
    int naoA;
    int naoB;
    double aux;
    double rAB2;
    AOs = basis.AtomNAOs(atoms);
    MBO.resize(Natoms,Natoms);
    MBO.zero();
    matrixE geom = mol.Geometry();
    matrixE PS = dens*sao;
    if ((shell == "open")&&(bdens.rows()*bdens.cols() > 1)) {PS += bdens*sao;}
    for (size_t irow = 0; irow < Natoms; ++irow) {
      icns = 1;
      naoA = AOs[irow];
      MBO(irow + 1,irow + 1) = 0.0;
      for (size_t icol = 0; icol < irow; ++icol) {
        naoB = AOs[icol];
        rAB2 = 0.0;
        for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
          aux = geom(irow + 1,idcoord) - geom(icol + 1,idcoord);
          rAB2 += aux*aux;
        }
        aux = 0.0;
        if (rAB2 < distance_cutoff) {
          for (size_t imu = 0; imu < naoA; ++imu) {
            for (size_t inu = 0; inu < naoB; ++inu) {
              aux += PS(icnt + imu,icns + inu)*PS(icns + inu,icnt + imu);
            }
          }
        }
        MBO(irow + 1,icol + 1) = aux;
        MBO(icol + 1,irow + 1) = aux;
        icns += naoB;
      }
      icnt += naoA;
    }
  }
  void ArmstrongBondOrder(matrixE & ABO, double distance_cutoff = 100.0*dist_Angstrom2au*dist_Angstrom2au) {
    //function returning Armstrong bond orders
    //D. R. Armstrong, R. G. Perkins, J. J. P. Stewart, J. Chem. Soc. Dalton Trans., 838, 1973
    int icnt = 1;
    int icns = 1;
    int naoA;
    int naoB;
    double aux;
    double rAB2;
    AOs = basis.AtomNAOs(atoms);
    ABO.resize(Natoms,Natoms);
    ABO.zero();
    matrixE geom = mol.Geometry();
    matrixE Pmat = dens;
    if ((shell == "open")&&(bdens.rows()*bdens.cols() > 1)) {Pmat += bdens;}
    for (size_t irow = 0; irow < Natoms; ++irow) {                //I could simply use here the number of rows for the matrix; by using atomvec I check at least that pBKKp has not too less atoms.
      icns = 1;
      naoA = AOs[irow];
      ABO(irow + 1,irow + 1) = 0.0;
      for (size_t icol = 0; icol < irow; ++icol) {
        naoB = AOs[icol];
        rAB2 = 0.0;
        for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
          aux = geom(irow + 1,idcoord) - geom(icol + 1,idcoord);
          rAB2 += aux*aux;
        }
        aux = 0.0;
        if (rAB2 < distance_cutoff) {
          for (size_t imu = 0; imu < naoA; ++imu) {
            for (size_t inu = 0; inu < naoB; ++inu) {
              aux += Pmat(icnt + imu,icns + inu)*Pmat(icnt + imu,icns + inu);
            }
          }
        }
        ABO(irow + 1,icol + 1) = aux;
        ABO(icol + 1,irow + 1) = aux;
        icns += naoB;
      }
      icnt += naoA;
    }
  }
  //reactivity functions
  double IonizationPotential(bool Koopman = true) {
    //function determining the ionization of potential for last calculated structure
    //if Koopman then use Koopman's theorem, otherwise use definition
    //Koopman only for closed-shells
    //ionization potential in Hartree
    double ip = 0.0;
    if (Koopman) {
      int homo = HOMOposition();
      ip = -EMOs[homo];
    }
    else {
      int rest = restart;
      restart = -1;
      Calculate(0);
      double E0 = getEnergy(1);
      int charge = mol.Charge();
      int multiplicity = mol.Multiplicity();
      int multiplicity_new = multiplicity - 1;
      if (multiplicity_new == 0) {multiplicity_new = 2;}     //if it was singlet, now should be duplet
      mol.setCharge(charge + 1);
      mol.setMultiplicity(multiplicity_new);
      Calculate(0);
      double Eplus = getEnergy(1);
      mol.setCharge(charge);
      mol.setMultiplicity(multiplicity);
      ip = Eplus - E0;
      restart = rest;
    }
    return ip;
  }
  double ElectronAffinity() {
    //function determining electron affinities for last calculated structure using the definition
    int rest = restart;
    restart = -1;
    Calculate(0);
    double E0 = getEnergy(1);
    int charge = mol.Charge();
    int multiplicity = mol.Multiplicity();
    int multiplicity_new = multiplicity - 1;
    if (multiplicity_new == 0) {multiplicity_new = 2;}       //if it was singlet, now should be duplet
    mol.setCharge(charge - 1);
    mol.setMultiplicity(multiplicity_new);
    Calculate(0);
    double Eminus = getEnergy(1);
    mol.setCharge(charge);
    mol.setMultiplicity(multiplicity);
    double ea = E0 - Eminus;
    restart = rest;
    return ea;
  }
  void HSABdata(double & chi, double & eta) {
    //function returning hardness and electronegativies for last calculated structure
    //chi -> electronegativity
    //eta -> hardness
    int rest = restart;
    restart = -1;
    Calculate(0);
    double E0 = getEnergy(1);
    int charge = mol.Charge();
    int multiplicity = mol.Multiplicity();
    int multiplicity_new = multiplicity - 1;
    if (multiplicity_new == 0) {multiplicity_new = 2;}       //if it was singlet, now should be duplet
    mol.setCharge(charge + 1);
    mol.setMultiplicity(multiplicity_new);
    Calculate(0);
    double Eplus = getEnergy(1);
    mol.setCharge(charge - 1);
    Calculate(0);
    double Eminus = getEnergy(1);
    mol.setCharge(charge);
    double ip = Eplus - E0;
    double ea = E0 - Eminus;
    restart = rest;
    chi = 0.5*(ip + ea);
    eta = 0.5*(ip - ea);
  }
  void FukuiIndices(std::vector<double> & fplus, std::vector<double> & fminus, std::vector<double> & f0, bool simple = false) {
    //function returning Fukui indices/functions for last calculated structure
    //f+ -> describes electrophilic centers
    //f- -> describes nucleophilic centers
    //f0 -> describes radical centers
    //negative values of f* functions indicate points in which density increases when adding/removing electrons
    //simple determines how to calculate the indices; if true, it goes via FMOs, meaning faster; otherwise goes via electronic populations
    fplus.resize(Natoms);
    fminus.resize(Natoms);
    f0.resize(Natoms);
    if (simple) {
      //book of Jensen
      int homo = HOMOposition();
      int lumo = homo + 1;
      double CHOMO;
      double CLUMO;
      std::vector<size_t> AOs = basis.AtomNAOs(atoms);
      size_t icnt = 1;
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        fplus[idAtm] = 0.0;
        fminus[idAtm] = 0.0;
        for (size_t idao = 0; idao < AOs[idAtm]; ++idao, ++icnt) {
          CHOMO = fabs(CMO(icnt,homo + 1));
          CLUMO = fabs(CMO(icnt,lumo + 1));
          if (fplus[idAtm] < CLUMO) {fplus[idAtm] = CLUMO;}
          if (fminus[idAtm] < CHOMO) {fminus[idAtm] = CHOMO;}
        }
        f0[idAtm] = 0.5*(fplus[idAtm] + fminus[idAtm]);
      }
    }
    else {
      std::vector<double> chg_plus;
      std::vector<double> chg_minus;
      std::vector<double> chg_0;
      int rest = restart;
      restart = -1;
      Calculate(0);
      if (Type() != "GFN2") {getMullikenCharges();}
      chg_0 = getQAtoms();
      int charge = mol.Charge();
      int multiplicity = mol.Multiplicity();
      int multiplicity_new = multiplicity - 1;
      if (multiplicity_new == 0) {multiplicity_new = 2;}       //if it was singlet, now should be duplet
      mol.setCharge(charge + 1);
      mol.setMultiplicity(multiplicity_new);
      Calculate(0);
      if (Type() != "GFN2") {getMullikenCharges();}
      chg_plus = getQAtoms();
      mol.setCharge(charge - 1);
      Calculate(0);
      if (Type() != "GFN2") {getMullikenCharges();}
      chg_minus = getQAtoms();
      mol.setCharge(charge);
      mol.setMultiplicity(multiplicity);
      restart = rest;
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        fplus[idAtm] = chg_0[idAtm] - chg_minus[idAtm];
        fminus[idAtm] = chg_plus[idAtm] - chg_0[idAtm];
        f0[idAtm] = 0.5*(fplus[idAtm] + fminus[idAtm]);
      }
    }
  }
  double ReactivityIndices(matrixE & RxData, bool simple = false) {
    //function determining the reactivity indices for a system
    //matrix with reactivity data: f+, f-, f0, softness+, softness-, softness0
    RxData.resize(Natoms,6);
    double softness = 0.0;
    Calculate(0);
    int rest = restart;
    restart = -1;
    double E0 = getEnergy(1);
    int charge = mol.Charge();
    int multiplicity = mol.Multiplicity();
    int multiplicity_new = multiplicity - 1;
    if (multiplicity_new == 0) {multiplicity_new = 2;}       //if it was singlet, now should be duplet
    if (simple) {
      int homo = HOMOposition();
      int lumo = homo + 1;
      double CHOMO;
      double CLUMO;
      std::vector<size_t> AOs = basis.AtomNAOs(atoms);
      size_t icnt = 1;
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        RxData(idAtm + 1,1) = 0.0;
        RxData(idAtm + 1,2) = 0.0;
        for (size_t idao = 0; idao < AOs[idAtm]; ++idao, ++icnt) {
          CHOMO = fabs(CMO(icnt,homo + 1));
          CLUMO = fabs(CMO(icnt,lumo + 1));
          if (RxData(idAtm + 1,1) < CLUMO) {RxData(idAtm + 1,1) = CLUMO;}
          if (RxData(idAtm + 1,2) < CHOMO) {RxData(idAtm + 1,2) = CHOMO;}
        }
        RxData(idAtm + 1,3) = 0.5*(RxData(idAtm + 1,1) + RxData(idAtm + 1,2));
      }
      mol.setCharge(charge + 1);
      mol.setMultiplicity(multiplicity_new);
      Calculate(0);
      double Eplus = getEnergy(1);
      mol.setCharge(charge - 1);
      Calculate(0);
      double Eminus = getEnergy(1);
      mol.setCharge(charge);
      mol.setMultiplicity(multiplicity);
      double ip = Eplus - E0;
      double ea = E0 - Eminus;
      restart = rest;
      softness = 2.0/(ip - ea);
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        RxData(idAtm + 1,4) = softness*RxData(idAtm + 1,1);
        RxData(idAtm + 1,5) = softness*RxData(idAtm + 1,2);
        RxData(idAtm + 1,6) = softness*RxData(idAtm + 1,3);
      }
    }
    else {
      std::vector<double> chg_plus;
      std::vector<double> chg_minus;
      std::vector<double> chg_0;
      if (Type() != "GFN2") {getMullikenCharges();}
      chg_0 = getQAtoms();
      mol.setCharge(charge + 1);
      mol.setMultiplicity(multiplicity_new);
      Calculate(0);
      double Eplus = getEnergy(1);
      if (Type() != "GFN2") {getMullikenCharges();}
      chg_plus = getQAtoms();
      mol.setCharge(charge - 1);
      Calculate(0);
      double Eminus = getEnergy(1);
      if (Type() != "GFN2") {getMullikenCharges();}
      chg_minus = getQAtoms();
      mol.setCharge(charge);
      mol.setMultiplicity(multiplicity);
      restart = rest;
      double ip = Eplus - E0;
      double ea = E0 - Eminus;
      softness = 2.0/(ip - ea);
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        RxData(idAtm + 1,1) = chg_0[idAtm] - chg_minus[idAtm];
        RxData(idAtm + 1,2) = chg_plus[idAtm] - chg_0[idAtm];
        RxData(idAtm + 1,3) = 0.5*(RxData(idAtm + 1,1) + RxData(idAtm + 1,2));
        RxData(idAtm + 1,4) = softness*RxData(idAtm + 1,1);
        RxData(idAtm + 1,5) = softness*RxData(idAtm + 1,2);
        RxData(idAtm + 1,6) = softness*RxData(idAtm + 1,3);
      }
    }
    return softness;
  }
  //solvation model stuff
  virtual double DeScreening(std::string solvent, size_t atmnr) {
    //descreening constant used to calculate Born-radii
    return 1.0;
  }
  virtual double BornOffset() {return 0.0;}
  virtual double FreeEnergyShift(std::string solvent_) {
    //free energy shift for ALPB model
    return 0.0;
  }
  virtual double BornRadiusScale(std::string solvent_) {
    //scaling parameter for Born radii
    return 0.0;
  }
  virtual double ProbeRadius(std::string solvent_) {
    //probe radii
    return 0.0;
  }
  virtual double GammaScale(std::string solvent, size_t atmnr) {
    //function returning the gamma's for surface tension
    return 0.0;
  }
  virtual double HbondStrength(std::string solvent, size_t atmnr) {
    //function returning the hydrogen-bonding strength of each atom
    return 0.0;
  }
  void IntegralNonOverlap(const double & RAB, const double & RABm1, const double & rho, double & I, double & dI, double & Ukk, double & Lkk, double & rhAB, double & logAB) {
    //function calculating the integration for Born radii in case atoms do not overlap
    //RAB and rho in Bohr; RABm1 in inverse Bohr; four last arguments are just to avoid constant alloc/delloc
    //based on equation 13 of 
    //A. V. Onufriev, D. A. Case, Annu. Rev. Biophys., 48, 275, 2019
    Ukk = RAB + rho;
    Lkk = RAB - rho;
    rhAB = rho/(Ukk*Lkk);
    logAB = 0.5*log(Lkk/Ukk)*RABm1;
    I = rhAB + logAB;
    dI = -2.0*rhAB/(Ukk*Lkk) + (rhAB - logAB)*RABm1*RABm1;
  }
  void IntegralOverlap(const double & RAB, const double & RABm1, const double & rho, double & I, double & dI, double & Ukk, double & Lkk, double & rhAB, double & logAB, double & RABm2, double & vdwm1, double & UkkvdWm1) {
    //function calculating the integration for Born radii in case atoms overlap
    //RAB and rho in Bohr; RABm1 in inverse Bohr
    //what is the reference for this one?
    RABm2 = RABm1*RABm1;
    Ukk = RAB + rho;
    Lkk = RAB - rho;
    rhAB = 1.0/Ukk;
    UkkvdWm1 = Ukk*vdwm1;
    logAB = log(UkkvdWm1);
    I = vdwm1 - rhAB + 0.5*RABm1*(0.5*Lkk*(rhAB - vdwm1*UkkvdWm1) - logAB);
    dI = rhAB*rhAB*(1.0 - 0.25*Lkk*RABm1*(1.0 + UkkvdWm1*UkkvdWm1))*RABm1 + 0.25*rho*RABm2*RABm1*(rhAB - vdwm1*UkkvdWm1) + 0.5*RABm2*(RABm1*logAB - rhAB);
  }
  void initGBVariables(double & shortrangecutoff) {
    //function that gets the short range cutoff distance
    double proberad = ProbeRadius(solvent);
    double ww = 0.3;
    double vdWradA;
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      vdWradA = vdWaalsRadii(atoms[idatm]) + proberad + ww;
      //update shortrangecutoff?
      if (vdWradA > shortrangecutoff) {shortrangecutoff = vdWradA;}
    }
    shortrangecutoff *= 2.0;
    shortrangecutoff += 2.0;
  }
  void calcBornRadii(std::vector<double> & BornRadii, matrixE & geom, matrix<int> & ListShortRange, double shortrangecutoff, double longrangecutoff = 35.0) {
    //function that calculates Born radii in Angstroem following
    //a) A. Onufriev, D. Bashford, D. A. Case, proteins, 55, 383, 2004
    if (BornRadii.size() != Natoms) {BornRadii.resize(Natoms);}
    double boffset = BornOffset();
    double bornradscal = BornRadiusScale(solvent);
    double rhoA;
    double rhoB;
    double vdWradA;
    double vdWradB;
    double aux;
    double RAB;
    double RAB2;
    double rhAB;
    double aprh1;
    double rh1;
    double ap;
    double am;
    double logAB;
    double grA;
    double grB;
    double dgrA;
    double dgrB;
    //parameters for GBOBC(II)
    double alphaLR = 1.0;           //equation 8, reference a)
    double betaLR = 0.8;            //equation 8, reference a)
    double gammaLR = 4.85;          //equation 8, reference a)
    int atomA;
    int atomB;
    bool ABoverlap;
    bool BAoverlap;
    std::vector<int> counter(Natoms,0);
    //get the integral I (actually 2I) of reference a)
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      BornRadii[idatm] = 0.0;
      atomA = atoms[idatm];
      vdWradA = vdWaalsRadii(atomA);
      rhoA = vdWradA*DeScreening(solvent,atomA);
      for (size_t idbtm = 0; idbtm < idatm; ++idbtm) {
        RAB2 = 0.0;
        for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
          rAB[idcoord - 1] = geom(idatm + 1,idcoord) - geom(idbtm + 1,idcoord);
          RAB2 += rAB[idcoord - 1]*rAB[idcoord - 1];
        }
        if (RAB2 > longrangecutoff*longrangecutoff) {continue;}
        if (RAB2 <= shortrangecutoff*shortrangecutoff) {
          ++counter[idatm];
          ++counter[idbtm];
          ListShortRange(idatm + 1,counter[idatm]) = idbtm + 1;
          ListShortRange(idbtm + 1,counter[idbtm]) = idatm + 1;
        }
        atomB = atoms[idbtm];
        RAB = sqrt(RAB2);
        vdWradB = vdWaalsRadii(atomB);
        rhoB = vdWradB*DeScreening(solvent,atomB);
        ABoverlap = (RAB < vdWradA + rhoB);
        BAoverlap = (RAB < vdWradB + rhoA);
        RAB2 = 1.0/RAB;
        grA = 0.0;
        grB = 0.0;
        if ((!ABoverlap)&&(!BAoverlap)) {             //AB and BA do not overlap
          IntegralNonOverlap(RAB,RAB2,rhoB,grA,dgrA,ap,am,rhAB,logAB);
          if (fabs(rhoA - rhoB) < 1.0e-8) {grB = grA;}
          else {IntegralNonOverlap(RAB,RAB2,rhoA,grB,dgrB,ap,am,rhAB,logAB);}
        }
        else if ((!ABoverlap)&&(BAoverlap)) {         //AB does not overlap but BA does
          IntegralNonOverlap(RAB,RAB2,rhoB,grA,dgrA,ap,am,rhAB,logAB);
          if (RAB + rhoA > vdWradB) {
            rh1 = 1.0/vdWradB;
            IntegralOverlap(RAB,RAB2,rhoA,grB,dgrB,ap,am,rhAB,logAB,aux,rh1,aprh1);
          }
        }
        else if ((ABoverlap)&&(!BAoverlap)) {         //BA does not overlap but AB does
          if (RAB + rhoB > vdWradA) {
            rh1 = 1.0/vdWradA;
            IntegralOverlap(RAB,RAB2,rhoB,grA,dgrA,ap,am,rhAB,logAB,aux,rh1,aprh1);
          }
          IntegralNonOverlap(RAB,RAB2,rhoA,grB,dgrB,ap,am,rhAB,logAB);
        }
        else if ((ABoverlap)&&(BAoverlap)) {          //AB and BA overlap
          if (RAB + rhoB > vdWradA) {
            rh1 = 1.0/vdWradA;
            IntegralOverlap(RAB,RAB2,rhoB,grA,dgrA,ap,am,rhAB,logAB,aux,rh1,aprh1);
          }
          if (RAB + rhoA > vdWradB) {
            rh1 = 1.0/vdWradB;
            IntegralOverlap(RAB,RAB2,rhoA,grB,dgrB,ap,am,rhAB,logAB,aux,rh1,aprh1);
          }
        }
        //accumulate radii
        BornRadii[idatm] += grA;
        BornRadii[idbtm] += grB;
      }
    }
    //now get the Born radii
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      atomA = atoms[idatm];
      vdWradA = vdWaalsRadii(atomA);
      vdWradB = vdWradA - boffset;                              //{\tilde{\rho}}_{i} in reference a)
      ap = 1.0/vdWradA;
      am = 1.0/vdWradB;
      aux = 0.5*BornRadii[idatm]*vdWradB;                       //{\Psi} in reference a)
      aprh1 = aux*(alphaLR + aux*(gammaLR*aux - betaLR));       //argument of tanh in equation 6, reference a)
      rh1 = tanh(aprh1);
      aux = 1.0/(am - ap*rh1);                                  //equation 6 in reference a)
      BornRadii[idatm] = aux*bornradscal*dist_Angstrom2aum1;
    }
  }
  double calcBornRadii_dX(std::vector<double> & BornRadii, std::vector<matrixE> & dBornRadii, std::vector<double> & comass, matrixE & geom, matrix<int> & ListShortRange, double shortrangecutoff, double longrangecutoff = 35.0) {
    //function that calculates the first derivatives of Born radii and it also gets the Born radii (these are for free) in units of Angstroem
    //a) A. Onufriev, D. Bashford, D. A. Case, proteins, 55, 383, 2004
    if (dBornRadii.size() != 3) {
      if ((BornMat.rows() != Natoms)||(BornMat.cols() != Natoms)) {BornMat.resize(Natoms,Natoms);}
      if (dBornRadii.size() > 0) {dBornRadii.clear();}
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        dBornRadii.push_back(BornMat);           //Born matrix is supposedly allocated by now
      }
    }
    else {
      bool check;
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        check = ((dBornRadii[idcoord].rows() == Natoms)&&(dBornRadii[idcoord].cols() == Natoms));
        if (!check) {dBornRadii[idcoord].resize(Natoms,Natoms);}
      }
    }
    if (BornRadii.size() != Natoms) {BornRadii.resize(Natoms);}
    double boffset = BornOffset();
    double bornradscal = BornRadiusScale(solvent);
    double rhoA;
    double rhoB;
    double vdWradA;
    double vdWradB;
    double aux;
    double RAB;
    double RAB2;
    double rhAB;
    double aprh1;
    double rh1;
    double ap;
    double am;
    double logAB;
    double grA;
    double grB;
    double dgrA;
    double dgrB;
    //parameters for GBOBC(II)
    double alphaLR = 1.0;           //equation 8, reference a)
    double betaLR = 0.8;            //equation 8, reference a)
    double gammaLR = 4.85;          //equation 8, reference a)
    double vdWr3;
    double trad3;
    int atomA;
    int atomB;
    bool ABoverlap;
    bool BAoverlap;
    std::vector<int> counter(Natoms,0);
    //get the integral I (actually 2I) of reference a) and its derivatives
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      BornRadii[idatm] = 0.0;
      for (size_t idx = 0; idx < 3; ++idx) {
        dBornRadii[idx](idatm + 1,idatm + 1) = 0.0;
        for (size_t idbtm = 0; idbtm < idatm; ++idbtm) {
          dBornRadii[idx](idatm + 1,idbtm + 1) = 0.0;
          dBornRadii[idx](idbtm + 1,idatm + 1) = 0.0;
        }
      }
      atomA = atoms[idatm];
      vdWradA = vdWaalsRadii(atomA);
      rhoA = vdWradA*DeScreening(solvent,atomA);
      for (size_t idbtm = 0; idbtm < idatm; ++idbtm) {
        RAB2 = 0.0;
        for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
          rAB[idcoord - 1] = geom(idatm + 1,idcoord) - geom(idbtm + 1,idcoord);
          RAB2 += rAB[idcoord - 1]*rAB[idcoord - 1];
        }
        if (RAB2 > longrangecutoff*longrangecutoff) {continue;}
        if (RAB2 <= shortrangecutoff*shortrangecutoff) {
          ++counter[idatm];
          ++counter[idbtm];
          ListShortRange(idatm + 1,counter[idatm]) = idbtm + 1;
          ListShortRange(idbtm + 1,counter[idbtm]) = idatm + 1;
        }
        atomB = atoms[idbtm];
        RAB = sqrt(RAB2);
        vdWradB = vdWaalsRadii(atomB);
        rhoB = vdWradB*DeScreening(solvent,atomB);
        ABoverlap = (RAB < vdWradA + rhoB);
        BAoverlap = (RAB < vdWradB + rhoA);
        RAB2 = 1.0/RAB;
        grA = 0.0;
        dgrA = 0.0;
        grB = 0.0;
        dgrB = 0.0;
        if ((!ABoverlap)&&(!BAoverlap)) {             //AB and BA do not overlap
          IntegralNonOverlap(RAB,RAB2,rhoB,grA,dgrA,ap,am,rhAB,logAB);
          if (fabs(rhoA - rhoB) < 1.0e-8) {
            grB = grA;
            dgrB = dgrA;
          }
          else {IntegralNonOverlap(RAB,RAB2,rhoA,grB,dgrB,ap,am,rhAB,logAB);}
        }
        else if ((!ABoverlap)&&(BAoverlap)) {         //AB does not overlap but BA does
          IntegralNonOverlap(RAB,RAB2,rhoB,grA,dgrA,ap,am,rhAB,logAB);
          if (RAB + rhoA > vdWradB) {
            rh1 = 1.0/vdWradB;
            IntegralOverlap(RAB,RAB2,rhoA,grB,dgrB,ap,am,rhAB,logAB,aux,rh1,aprh1);
          }
        }
        else if ((ABoverlap)&&(!BAoverlap)) {         //BA does not overlap but AB does
          if (RAB + rhoB > vdWradA) {
            rh1 = 1.0/vdWradA;
            IntegralOverlap(RAB,RAB2,rhoB,grA,dgrA,ap,am,rhAB,logAB,aux,rh1,aprh1);
          }
          IntegralNonOverlap(RAB,RAB2,rhoA,grB,dgrB,ap,am,rhAB,logAB);
        }
        else if ((ABoverlap)&&(BAoverlap)) {          //AB and BA overlap
          if (RAB + rhoB > vdWradA) {
            rh1 = 1.0/vdWradA;
            IntegralOverlap(RAB,RAB2,rhoB,grA,dgrA,ap,am,rhAB,logAB,aux,rh1,aprh1);
          }
          if (RAB + rhoA > vdWradB) {
            rh1 = 1.0/vdWradB;
            IntegralOverlap(RAB,RAB2,rhoA,grB,dgrB,ap,am,rhAB,logAB,aux,rh1,aprh1);
          }
        }
        //accumulate radii
        BornRadii[idatm] += grA;
        BornRadii[idbtm] += grB;
        //accumulate gradient
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          dBornRadii[idcoord](idatm + 1,idatm + 1) += dgrA*rAB[idcoord];
          dBornRadii[idcoord](idatm + 1,idbtm + 1) += dgrB*rAB[idcoord];
          dBornRadii[idcoord](idbtm + 1,idatm + 1) -= dgrA*rAB[idcoord];
          dBornRadii[idcoord](idbtm + 1,idbtm + 1) -= dgrB*rAB[idcoord];
        }
      }
    }
    //now get the Born radii and their derivatives
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      atomA = atoms[idatm];
      vdWradA = vdWaalsRadii(atomA);
      vdWradB = vdWradA - boffset;                              //{\tilde{\rho}}_{i} in reference a)
      ap = 1.0/vdWradA;
      am = 1.0/vdWradB;
      aux = 0.5*BornRadii[idatm]*vdWradB;                       //{\Psi} in reference a)
      rhAB = aux*(gammaLR*aux - betaLR);
      aprh1 = aux*(alphaLR + rhAB);                             //argument of tanh in equation 6, reference a)
      rhAB *= 2.0;
      rhAB += alphaLR + gammaLR*aux*aux;
      rh1 = tanh(aprh1);
      logAB = cosh(aprh1);
      aux = 1.0/(am - ap*rh1);                                  //equation 6 in reference a)
      BornRadii[idatm] = aux*bornradscal*dist_Angstrom2aum1;
      rhoA = logAB*(am - ap*rh1);
      rhoB = 0.5*vdWradB*ap*rhAB/(rhoA*rhoA);
      rhoA = bornradscal*rhoB;
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        for (size_t idbtm = 0; idbtm < Natoms; ++idbtm) {
          dBornRadii[idcoord](idbtm + 1,idatm + 1) *= rhoA;
        }
      }
      vdWradA *= dist_Angstrom2aum1;
      vdWr3 = vdWradA*vdWradA*vdWradA;
      trad3 += vdWr3;
      comass[0] += geom(idatm + 1,1)*dist_Angstrom2aum1*vdWr3;
      comass[1] += geom(idatm + 1,2)*dist_Angstrom2aum1*vdWr3;
      comass[2] += geom(idatm + 1,3)*dist_Angstrom2aum1*vdWr3;
    }
    return trad3;
  }
  double calcSolventAccessibleSurfaceArea(std::vector<double> & hbond, matrixE & geom, matrix<int> & ListShortRange, double shortrangecutoff, double threshzero = 1.0e-6) {
    //W. Im, M. S. Lee, C. L. Brooks, J. Comput. Chem. 24(14), 1691, 2003
    double sasa = 0.0;
    double sasaA;
    double vdWradsa;
    double proberad = ProbeRadius(solvent);
    double wrp;
    double aux;
    double sasapoint;
    double tB2;
    double vdWsapw;
    double vdWsamw;
    double vdWsa;
    double unitconv = 1.0/au2kcalmol;
    //parameters to determine radial weights for numerical integration
    double ww = 0.3*dist_Angstrom2aum1;
    double w3 = ww*ww*ww;
    double ah0 = 0.5;
    double ah1 = 0.75/ww;
    double ah3 = -0.25/w3;
    //surface tension parameter
    double stparam = 1.0e-5;    //Eh
    int atomA;
    int position;
    int idbtm;
    //get Grid
    std::vector<int> lebgrid = LebedevGrid();
    BisectionSearch(position,ngrid,lebgrid);
    int ngridpoints = lebgrid[position - 1];
    std::vector<double> GridWeights(ngridpoints);
    double rgp[3];                                  //position of the grid point
    double tB[3];
    matrixE Grid(3,ngridpoints);
    getLebedevLaikovGrid(position,GridWeights,Grid);
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      atomA = atoms[idatm];
      sasaA = 0.0;
      vdWradsa = (vdWaalsRadii(atomA) + proberad)*dist_Angstrom2aum1;
      rAB[0] = geom(idatm + 1,1)*dist_Angstrom2aum1;
      rAB[1] = geom(idatm + 1,2)*dist_Angstrom2aum1;
      rAB[2] = geom(idatm + 1,3)*dist_Angstrom2aum1;
      //get radial weights for the numerical integration
      aux = vdWradsa + ww;
      wrp = (0.25/ww + 3.0*ah3*(0.2*aux*aux - 0.5*aux*vdWradsa + vdWradsa*vdWradsa/3.0))*aux*aux*aux;
      aux = vdWradsa - ww;
      wrp -= (0.25/ww + 3.0*ah3*(0.2*aux*aux - 0.5*aux*vdWradsa + vdWradsa*vdWradsa/3.0))*aux*aux*aux;
      //now loop over grid points
      for (size_t igrpnt = 0; igrpnt < ngridpoints; ++igrpnt) {
        //get the grid point
        rgp[0] = rAB[0] + vdWradsa*Grid(1,igrpnt + 1);
        rgp[1] = rAB[1] + vdWradsa*Grid(2,igrpnt + 1);
        rgp[2] = rAB[2] + vdWradsa*Grid(3,igrpnt + 1);
        sasapoint = 1.0;
        //get atomic surface at grid point
        for (size_t ibtm = 0; ibtm < Natoms; ++ibtm) {
          idbtm = ListShortRange(idatm + 1,ibtm + 1) - 1;
          if (idbtm < 0) {break;}
          tB2 = 0.0;
          for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
            tB[idcoord] = rgp[idcoord] - geom(idbtm + 1,idcoord + 1)*dist_Angstrom2aum1;
            tB2 += tB[idcoord]*tB[idcoord];
          }
          vdWsa = (vdWaalsRadii(atoms[idbtm]) + proberad)*dist_Angstrom2aum1;
          vdWsapw = vdWsa + ww;
          vdWsamw = vdWsa - ww;
          if (tB2 < vdWsapw*vdWsapw) {
            if (tB2 <= vdWsamw*vdWsamw) {
              sasapoint = 0.0;
              break;
            }
            else {
              vdWsapw = sqrt(tB2);
              aux = vdWsapw - vdWsa;
              vdWsamw = ah3*aux*aux;
              tB2 = ah0 + (ah1 + vdWsamw)*aux;
              sasapoint *= tB2;                   //accumulate molecular surface
            }
          }
        }
        if(sasapoint > threshzero) {
          aux = GridWeights[igrpnt]*wrp*sasapoint;                    //weight of numerical quadrature
          sasaA += aux;                                               //accumulate to surface area
        }
      }
      aux = HbondStrength(solvent,atomA);
      hbond[idatm] = -sasaA*aux*aux*unitconv/(vdWradsa*vdWradsa);
      sasa += sasaA*GammaScale(solvent,atomA);
    }
    return 4.0*pi*stparam*sasa;
  }
  double calcSolventAccessibleSurfaceArea_dX(matrixE & dsasa, std::vector<matrixE> & dSASAdrt, std::vector<double> & hbond, std::vector<double> & dhbond, matrixE & geom, matrix<int> & ListShortRange, double shortrangecutoff, double threshzero = 1.0e-6) {
    if (dsasa.rows() != 3*Natoms) {dsasa.resize(3*Natoms,1);}
    dsasa.zero();
    dSASAdrt[0].zero();
    dSASAdrt[1] = dSASAdrt[0];
    dSASAdrt[2] = dSASAdrt[0];
    double sasa = 0.0;
    double sasaA;
    double vdWradsa;
    double proberad = ProbeRadius(solvent);
    double wrp;
    double aux;
    double sasapoint;
    double tB2;
    double vdWsapw;
    double vdWsamw;
    double vdWsa;
    double gammascaA;
    double unitconv = 1.0/au2kcalmol;
    //parameters to determine radial weights for numerical integration
    double ww = 0.3*dist_Angstrom2aum1;
    double w3 = ww*ww*ww;
    double ah0 = 0.5;
    double ah1 = 0.75/ww;
    double ah3 = -0.25/w3;
    //surface tension parameter
    double stparam = 1.0e-5;    //Eh
    int atomA;
    int position;
    int idBtm;
    int npoints;
    //get Grid
    std::vector<int> lebgrid = LebedevGrid();
    BisectionSearch(position,ngrid,lebgrid);
    int ngridpoints = lebgrid[position - 1];
    std::vector<double> GridWeights(ngridpoints);
    double rgp[3];                                  //position of the grid point
    double tB[3];
    std::vector<int> grdi(Natoms);
    matrixE Grid(3,ngridpoints);
    matrixE grds(Natoms,3);
    getLebedevLaikovGrid(position,GridWeights,Grid);
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      atomA = atoms[idAtm];
      sasaA = 0.0;
      vdWradsa = (vdWaalsRadii(atomA) + proberad)*dist_Angstrom2aum1;
      rAB[0] = geom(idAtm + 1,1)*dist_Angstrom2aum1;
      rAB[1] = geom(idAtm + 1,2)*dist_Angstrom2aum1;
      rAB[2] = geom(idAtm + 1,3)*dist_Angstrom2aum1;
      gammascaA = GammaScale(solvent,atomA);
      //get radial weights for the numerical integration
      aux = vdWradsa + ww;
      wrp = (0.25/ww + 3.0*ah3*(0.2*aux*aux - 0.5*aux*vdWradsa + vdWradsa*vdWradsa/3.0))*aux*aux*aux;
      aux = vdWradsa - ww;
      wrp -= (0.25/ww + 3.0*ah3*(0.2*aux*aux - 0.5*aux*vdWradsa + vdWradsa*vdWradsa/3.0))*aux*aux*aux;
      //now loop over grid points
      for (size_t igrpnt = 0; igrpnt < ngridpoints; ++igrpnt) {
        //get the grid point
        rgp[0] = rAB[0] + vdWradsa*Grid(1,igrpnt + 1);
        rgp[1] = rAB[1] + vdWradsa*Grid(2,igrpnt + 1);
        rgp[2] = rAB[2] + vdWradsa*Grid(3,igrpnt + 1);
        npoints = 0;
        sasapoint = 1.0;
        //get atomic surface at grid point
        for (size_t ibtm = 0; ibtm < Natoms; ++ibtm) {
          idBtm = ListShortRange(idAtm + 1,ibtm + 1) - 1;
          if (idBtm < 0) {break;}
          tB2 = 0.0;
          for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
            tB[idcoord] = rgp[idcoord] - geom(idBtm + 1,idcoord + 1)*dist_Angstrom2aum1;
            tB2 += tB[idcoord]*tB[idcoord];
          }
          vdWsa = (vdWaalsRadii(atoms[idBtm]) + proberad)*dist_Angstrom2aum1;
          vdWsapw = vdWsa + ww;
          vdWsamw = vdWsa - ww;
          if (tB2 < vdWsapw*vdWsapw) {
            if (tB2 <= vdWsamw*vdWsamw) {
              sasapoint = 0.0;
              break;
            }
            else {
              vdWsapw = sqrt(tB2);
              aux = vdWsapw - vdWsa;
              vdWsamw = ah3*aux*aux;
              tB2 = ah0 + (ah1 + vdWsamw)*aux;
              sasapoint *= tB2;                   //accumulate molecular surface
              aux = ah1 + 3.0*vdWsamw;
              aux /= tB2*vdWsapw;                 //gradients now
              grdi[npoints] = idBtm;
              ++npoints;
              grds(npoints,1) = aux*tB[0];
              grds(npoints,2) = aux*tB[1];
              grds(npoints,3) = aux*tB[2];
            }
          }
        }
        if(sasapoint > threshzero) {
          aux = GridWeights[igrpnt]*wrp*sasapoint;                    //weight of numerical quadrature
          sasaA += aux;                                               //accumulate to surface area
          for (size_t ipoint = 0; ipoint < npoints; ++ipoint) {       //accumulate to surface gradient
            idBtm = grdi[ipoint];
            for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
              dSASAdrt[idcoord](idAtm + 1,idAtm + 1) += aux*grds(ipoint + 1,idcoord + 1);
              dSASAdrt[idcoord](idBtm + 1,idAtm + 1) -= aux*grds(ipoint + 1,idcoord + 1);
              dsasa(3*idAtm + idcoord + 1,1) += 4.0*pi*stparam*aux*grds(ipoint + 1,idcoord + 1)*gammascaA;
              dsasa(3*idBtm + idcoord + 1,1) -= 4.0*pi*stparam*aux*grds(ipoint + 1,idcoord + 1)*gammascaA;
            }
          }
        }
      }
      aux = HbondStrength(solvent,atomA);
      hbond[idAtm] = -sasaA*aux*aux*unitconv/(vdWradsa*vdWradsa);
      dhbond[idAtm] = -aux*aux*unitconv/(vdWradsa*vdWradsa);
      sasa += sasaA*gammascaA;
    }
    return 4.0*pi*stparam*sasa;
  }
  double calcAdeterminant(matrixE & geom) {
    //function calculating Adeterminant as defined in equation 23 of
    //G. Sigalov, A. Fenley, A. Onufriev, J. Chem. Phys., 124, 124902, 2006
    double vdWr;
    double vdWr2;
    double vdWr3;
    double sumvdWr3 = 0.0;
    double comass[3] = {0.0,0.0,0.0};
    matrixE inertia(3,3);
    //get van der Waals based center of mass
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      vdWr = vdWaalsRadii(atoms[idatm])*dist_Angstrom2aum1;
      vdWr3 = vdWr*vdWr*vdWr;
      sumvdWr3 += vdWr3;
      comass[0] += geom(idatm + 1,1)*dist_Angstrom2aum1*vdWr3;
      comass[1] += geom(idatm + 1,2)*dist_Angstrom2aum1*vdWr3;
      comass[2] += geom(idatm + 1,3)*dist_Angstrom2aum1*vdWr3;
    }
    comass[0] /= sumvdWr3;
    comass[1] /= sumvdWr3;
    comass[2] /= sumvdWr3;
    //get inertia of ellipsoid
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      vdWr = vdWaalsRadii(atoms[idatm])*dist_Angstrom2aum1;
      vdWr2 = vdWr*vdWr;
      vdWr3 = vdWr2*vdWr;
      rAB[0] = geom(idatm + 1,1)*dist_Angstrom2aum1 - comass[0];
      rAB[1] = geom(idatm + 1,2)*dist_Angstrom2aum1 - comass[1];
      rAB[2] = geom(idatm + 1,3)*dist_Angstrom2aum1 - comass[2];
      for (size_t idrow = 0; idrow < 3; ++idrow) {
        for (size_t idcol = 0; idcol < 3; ++idcol) {
          if (idrow != idcol) {
            inertia(idrow + 1,idrow + 1) += vdWr3*(rAB[idcol]*rAB[idcol] + 0.2*vdWr2);
            inertia(idrow + 1,idcol + 1) -= vdWr3*rAB[idrow]*rAB[idcol];
          }
        }
      }
    }
    comass[0] = 0.5*(inertia(1,2) + inertia(2,1));
    comass[1] = 0.5*(inertia(1,3) + inertia(3,1));
    comass[2] = 0.5*(inertia(3,2) + inertia(2,3));
    vdWr = inertia(1,1)*inertia(2,2)*inertia(3,3) + 2.0*comass[0]*comass[1]*comass[2];
    vdWr -= inertia(1,1)*comass[2]*comass[2] + inertia(2,2)*comass[1]*comass[1] + inertia(3,3)*comass[0]*comass[0];
    return sqrt(2.5*cbrt(vdWr)/sumvdWr3);
  }
  double DebyeHueckel(std::vector<double> & ionexcl, double betaGB,std::vector<double> & BornRadii) {
    //ion-exclusion term from Debye-Hueckel
    double kappaDH = sqrt(saltmolarity*betaGB/(0.7897e-3*temperature))*dist_Angstrom2au;
    double aux;
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      aux = kappaDH*(BornRadii[idatm] + saltradius);
      ionexcl[idatm] = 0.5*exp(aux)/(1.0 + aux);
    }
    return kappaDH;
  }
  double dDebyeHueckel(std::vector<double> & ionexcl, std::vector<double> & dionexcl, double betaGB,std::vector<double> & BornRadii) {
    //ion-exclusion term from Debye-Hueckel and its derivatives
    double kappaDH = sqrt(saltmolarity*betaGB/(0.7897e-3*temperature))*dist_Angstrom2au;
    double aux;
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      aux = kappaDH*(BornRadii[idatm] + saltradius);
      ionexcl[idatm] = 0.5*exp(aux)/(1.0 + aux);
      dionexcl[idatm] = ionexcl[idatm]*kappaDH*aux/(1.0 + aux);
    }
    return kappaDH;
  }
  void calcBornMatrix(matrixE & geom, double (*EffCoulModel)(double, double, double), double zerosaltmolarity = 1.0e-10) {
    //function that calculates the Born matrix in generalized Born models
    //based on the generalization from ALPB method, equation 25 of
    //G. Sigalov, A. Fenley, A. Onufriev, J. Chem. Phys., 124, 124902, 2006
    //as implemented in xTB
    //EffCoulModel is a function defining the effective Coulomb operator (c.f. ChargePackage)
    BornMat.resize(Natoms,Natoms);
    matrix<int> ListShortRange(Natoms,Natoms);      //matrix containing pairs of atoms not further than "shortrange" distance
    std::vector<double> BornRadii(Natoms);
    std::vector<double> hbond(Natoms);
    std::vector<double> ionexcl(Natoms,0.5);
    double shortrangecutoff = 0.0;
    double gAB;
    double fAB;
    double alphaGB = 0.571412*(solvationmodel == "ALPB");
    double betaGB = 1.0/SolventDiElectricConstant(solvent,this->Type());
    double kappaDH = 0.0;
    double gamma;
    double delta;
    double Adet = calcAdeterminant(geom);
    initGBVariables(shortrangecutoff);
    calcBornRadii(BornRadii,geom,ListShortRange,shortrangecutoff);
    SASA = calcSolventAccessibleSurfaceArea(hbond,geom,ListShortRange,shortrangecutoff);
    //Debye-Hueckel
    if (saltmolarity > zerosaltmolarity) {
      //salt molarity must be positive and non-zero
      kappaDH = DebyeHueckel(ionexcl,betaGB,BornRadii);
    }
    betaGB *= epsilon_in;
    gamma = 1.0/(epsilon_in*(1.0 + alphaGB*betaGB));
    delta = alphaGB*betaGB/Adet;
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      fAB = BornRadii[idatm];
      alphaGB = 1.0 - 2.0*ionexcl[idatm]*betaGB*exp(-kappaDH*fAB);     //Debye-Hueckel screening term
      gAB = gamma*alphaGB*(1.0/fAB + delta);
      BornMat(idatm + 1,idatm + 1) = gAB - 2.0*hbond[idatm];
      for (size_t idbtm = 0; idbtm < idatm; ++idbtm) {
        gAB = 0.0;
        for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
          fAB = (geom(idatm + 1,idcoord) - geom(idbtm + 1,idcoord))*dist_Angstrom2aum1;
          gAB += fAB*fAB;
        }
        fAB = (*EffCoulModel)(gAB,BornRadii[idatm],BornRadii[idbtm]);
        alphaGB = 1.0 - (ionexcl[idatm] + ionexcl[idbtm])*betaGB*exp(-kappaDH*fAB);     //Debye-Hueckel screening term
        gAB = gamma*alphaGB*(1.0/fAB + delta);
        BornMat(idatm + 1,idbtm + 1) = gAB;
        BornMat(idbtm + 1,idatm + 1) = gAB;
      }
    }
  }
  void calcBornMatrix_dX(matrixE & gradient, std::vector<matrixE> & dBornRadii, matrixE & aux1, matrixE & aux2, matrixE & aux3, matrixE & geom, void (*gEffCoulModel)(double, double, double, double &, double &, double), double zerosaltmolarity = 1.0e-10, double thresholdzero = 1.0e-8) {
    //function that calculates the first derivative of the Born matrix in generalized Born models
    //based on the generalization from ALPB method, equation 25 of
    //G. Sigalov, A. Fenley, A. Onufriev, J. Chem. Phys., 124, 124902, 2006
    //as implemented in xTB
    //EffCoulModel is a function defining the effective Coulomb operator (c.f. ChargePackage)
    dBornRadii[0] = BornMat;
    dBornRadii[1] = BornMat;
    dBornRadii[2] = BornMat;
    matrixE dSASA(3*Natoms,1);
    matrixE dEdBr(Natoms,1);
    matrixE inertia(3,3);
    inertia.zero();
    matrix<int> ListShortRange(Natoms,Natoms);      //matrix containing pairs of atoms not further than "shortrange" distance
    std::vector<double> BornRadii(Natoms);
    std::vector<double> hbond(Natoms);
    std::vector<double> dhbond(Natoms);
    std::vector<double> comass(3,0.0);
    double shortrangecutoff = 0.0;
    double rAB2;
    double drAB;
    double dRI;
    double kqAqB;
    double alphaGB = 0.571412*(solvationmodel == "ALPB");
    double alphap;
    double betaGB = epsilon_in/SolventDiElectricConstant(solvent,this->Type());
    double gamma;
    double charge = 0.0;
    double vdWradA;
    double vdWr2;
    double vdWr3;
    double Adet;
    double coeff[6];
    initGBVariables(shortrangecutoff);
    double trad3 = calcBornRadii_dX(BornRadii,dBornRadii,comass,geom,ListShortRange,shortrangecutoff);
    comass[0] /= trad3;
    comass[1] /= trad3;
    comass[2] /= trad3;
    //Debye-Hueckel not implemented in gradients
    gamma = 1.0/(epsilon_in*(1.0 + alphaGB*betaGB));
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      alphap = 1.0 - betaGB;
      rAB2 = 1.0/BornRadii[idAtm];
      kqAqB = QAtoms[idAtm]*rAB2;
      dEdBr(idAtm + 1,1) = 0.5*gamma*alphap*kqAqB*rAB2*QAtoms[idAtm];
      for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
        rAB2 = 0.0;
        for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
          rAB[idcoord - 1] = (geom(idAtm + 1,idcoord) - geom(idBtm + 1,idcoord))*dist_Angstrom2aum1;
          rAB2 += rAB[idcoord - 1]*rAB[idcoord - 1];
        }
        (*gEffCoulModel)(rAB2,BornRadii[idAtm],BornRadii[idBtm],drAB,dRI,gamma*alphap);
        kqAqB = QAtoms[idAtm]*QAtoms[idBtm];
        drAB *= kqAqB;
        dRI *= kqAqB;
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          gradient(3*idAtm + idcoord + 1,1) += rAB[idcoord]*drAB;
          gradient(3*idBtm + idcoord + 1,1) -= rAB[idcoord]*drAB;
        }
        dEdBr(idAtm + 1,1) -= BornRadii[idBtm]*dRI;
        dEdBr(idBtm + 1,1) -= BornRadii[idAtm]*dRI;
      }
      //for derivatives of Adet, the inertia of the ellipsoid
      vdWradA = vdWaalsRadii(atoms[idAtm])*dist_Angstrom2aum1;
      vdWr2 = vdWradA*vdWradA;
      vdWr3 = vdWr2*vdWradA;
      rAB[0] = geom(idAtm + 1,1)*dist_Angstrom2aum1 - comass[0];
      rAB[1] = geom(idAtm + 1,2)*dist_Angstrom2aum1 - comass[1];
      rAB[2] = geom(idAtm + 1,3)*dist_Angstrom2aum1 - comass[2];
      for (size_t idrow = 0; idrow < 3; ++idrow) {
        for (size_t idcol = 0; idcol < 3; ++idcol) {
          if (idrow != idcol) {
            inertia(idrow + 1,idrow + 1) += vdWr3*(rAB[idcol]*rAB[idcol] + 0.2*vdWr2);
            inertia(idrow + 1,idcol + 1) -= vdWr3*rAB[idrow]*rAB[idcol];
          }
        }
      }
      charge += QAtoms[idAtm];
    }
    aux1 = dBornRadii[0]*dEdBr;
    aux2 = dBornRadii[1]*dEdBr;
    aux3 = dBornRadii[2]*dEdBr;
    SASA = calcSolventAccessibleSurfaceArea_dX(dSASA,dBornRadii,hbond,dhbond,geom,ListShortRange,shortrangecutoff);
    //last terms for derivatives of Adet
    rAB[0] = 0.5*(inertia(1,2) + inertia(2,1));
    rAB[1] = 0.5*(inertia(1,3) + inertia(3,1));
    rAB[2] = 0.5*(inertia(3,2) + inertia(2,3));
    vdWradA = inertia(1,1)*inertia(2,2)*inertia(3,3) + 2.0*rAB[0]*rAB[1]*rAB[2];
    vdWradA -= inertia(1,1)*rAB[2]*rAB[2] + inertia(2,2)*rAB[1]*rAB[1] + inertia(3,3)*rAB[0]*rAB[0];
    Adet = sqrt(2.5*cbrt(vdWradA)/trad3);
    vdWr2 = Adet*Adet;
    vdWradA = (250.0/(48.0*trad3*trad3*trad3*vdWr2*vdWr2*Adet))*(0.5*gamma*(1.0 - betaGB)*alphaGB*betaGB*charge*charge/vdWr2);
    coeff[0] = vdWradA*(inertia(1,1)*(inertia(2,2) + inertia(3,3)) - inertia(1,2)*inertia(1,2) - inertia(1,3)*inertia(1,3));
    coeff[1] = vdWradA*(inertia(1,2)*inertia(3,3) - inertia(1,3)*inertia(2,3));
    coeff[2] = vdWradA*(inertia(2,2)*(inertia(1,1) + inertia(3,3)) - inertia(1,2)*inertia(1,2) - inertia(2,3)*inertia(2,3));
    coeff[3] = vdWradA*(inertia(1,3)*inertia(2,2) - inertia(1,2)*inertia(3,2));
    coeff[4] = vdWradA*(inertia(1,1)*inertia(2,3) - inertia(1,2)*inertia(1,3));
    coeff[5] = vdWradA*(inertia(3,3)*(inertia(1,1) + inertia(2,2)) - inertia(1,3)*inertia(1,3) - inertia(2,3)*inertia(2,3));
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      vdWradA = vdWaalsRadii(atoms[idAtm])*dist_Angstrom2aum1;
      vdWr2 = vdWradA*vdWradA;
      vdWr3 = vdWr2*vdWradA;
      rAB[0] = geom(idAtm + 1,1)*dist_Angstrom2aum1 - comass[0];
      rAB[1] = geom(idAtm + 1,2)*dist_Angstrom2aum1 - comass[1];
      rAB[2] = geom(idAtm + 1,3)*dist_Angstrom2aum1 - comass[2];
      gradient(3*idAtm + 1,1) += aux1(idAtm + 1,1) + dSASA(3*idAtm + 1,1) + vdWr3*(coeff[0]*rAB[0] + coeff[1]*rAB[1] + coeff[3]*rAB[2]);
      gradient(3*idAtm + 2,1) += aux2(idAtm + 1,1) + dSASA(3*idAtm + 2,1) + vdWr3*(coeff[1]*rAB[0] + coeff[2]*rAB[1] + coeff[4]*rAB[2]);
      gradient(3*idAtm + 3,1) += aux3(idAtm + 1,1) + dSASA(3*idAtm + 3,1) + vdWr3*(coeff[3]*rAB[0] + coeff[4]*rAB[1] + coeff[5]*rAB[2]);
      //for derivatives of Hbond terms
      if (fabs(dhbond[idAtm]) < thresholdzero) {continue;}
      kqAqB = QAtoms[idAtm]*QAtoms[idAtm]*dhbond[idAtm];
      for (size_t idBtm = 0; idBtm < Natoms; ++idBtm) {
        gradient(3*idBtm + 1,1) += kqAqB*dBornRadii[0](idBtm + 1,idAtm + 1);
        gradient(3*idBtm + 2,1) += kqAqB*dBornRadii[1](idBtm + 1,idAtm + 1);
        gradient(3*idBtm + 3,1) += kqAqB*dBornRadii[2](idBtm + 1,idAtm + 1);
      }
    }
  }
  void calcSolvationForces(std::vector<double> & force, std::vector<double> & charges) {
    //function that calculates the force on atoms due to solvent effects
    if (force.size() != Natoms) {force.resize(Natoms);}
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      force[idatm] = 0.0;
      for (size_t idbtm = 0; idbtm < Natoms; ++idbtm) {
        force[idatm] -= BornMat(idatm + 1,idbtm + 1)*charges[idbtm];
      }
    }
  }
  double SolvationEnergy(std::vector<double> & charges) {
    //function that calculates the solvation energy on a molecule
    double Esolv = SASA + FreeEnergyShift(solvent);
    double aux;
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      aux = 0.0;
      for (size_t idbtm = 0; idbtm < Natoms; ++idbtm) {
        aux -= BornMat(idatm + 1,idbtm + 1)*charges[idbtm];
      }
      Esolv += 0.5*charges[idatm]*aux;
    }
    return Esolv;
  }
  matrixE CurrentSolvationEnergyMatrix(matrixE & geom) {
    //function returning the current solvation energy matrix, i.e., with the last set of charges
    matrix<int> ListShortRange(Natoms,Natoms);      //matrix containing pairs of atoms not further than "shortrange" distance
    double shortrangecutoff;
    double RAB2;
    double aux;
    std::vector<int> counter(Natoms,0);
    initGBVariables(aux);
    shortrangecutoff = aux*aux;
    //get the list of close neighbours; note that this calculation is embedded in the loops in the calculation of the Born-Matrix
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      for (size_t idbtm = 0; idbtm < idatm; ++idbtm) {
        RAB2 = 0.0;
        for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
          aux = geom(idatm + 1,idcoord) - geom(idbtm + 1,idcoord);
          RAB2 += aux*aux;
        }
        if (RAB2 <= shortrangecutoff) {
          ++counter[idatm];
          ++counter[idbtm];
          ListShortRange(idatm + 1,counter[idatm]) = idbtm + 1;
          ListShortRange(idbtm + 1,counter[idbtm]) = idatm + 1;
        }
      }
    }
    return SolvationEnergyMatrix(QAtoms,geom,ListShortRange);
  }
  matrixE SolvationEnergyMatrix(std::vector<double> & charges, matrixE & geom, matrix<int> & ListShortRange, double threshzero = 1.0e-6) {
    //function that calculates the solvation energy on a molecule as a matrix
    matrixE smat(Natoms,Natoms);
    double free_energy = FreeEnergyShift(solvent)/Natoms;
    double stparam = 1.0e-5;    //surface tension parameter in Eh
    double sasaA;
    double aux;
    double sasapoint;
    double wrp;
    double tB2;
    double ww = 0.3*dist_Angstrom2aum1;
    double w3 = ww*ww*ww;
    double ah0 = 0.5;
    double ah1 = 0.75/ww;
    double ah3 = -0.25/w3;
    double vdWsamw;
    double vdWradsa;
    double proberad = ProbeRadius(solvent);
    double vdWsapw;
    double vdWsa;
    int atomA;
    int position;
    int idbtm;
    //get Grid
    std::vector<int> lebgrid = LebedevGrid();
    BisectionSearch(position,ngrid,lebgrid);
    int ngridpoints = lebgrid[position - 1];
    std::vector<double> GridWeights(ngridpoints);
    double rgp[3];                                  //position of the grid point
    double tB[3];
    matrixE Grid(3,ngridpoints);
    getLebedevLaikovGrid(position,GridWeights,Grid);
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      //Born model contributions as well as free energy
      smat(idatm + 1,idatm + 1) = -0.5*BornMat(idatm + 1,idatm + 1)*charges[idatm]*charges[idatm] + free_energy;
      //SASA terms
      atomA = atoms[idatm];
      sasaA = 0.0;
      vdWradsa = (vdWaalsRadii(atomA) + proberad)*dist_Angstrom2aum1;
      rAB[0] = geom(idatm + 1,1)*dist_Angstrom2aum1;
      rAB[1] = geom(idatm + 1,2)*dist_Angstrom2aum1;
      rAB[2] = geom(idatm + 1,3)*dist_Angstrom2aum1;
      //get radial weights for the numerical integration
      aux = vdWradsa + ww;
      wrp = (0.25/ww + 3.0*ah3*(0.2*aux*aux - 0.5*aux*vdWradsa + vdWradsa*vdWradsa/3.0))*aux*aux*aux;
      aux = vdWradsa - ww;
      wrp -= (0.25/ww + 3.0*ah3*(0.2*aux*aux - 0.5*aux*vdWradsa + vdWradsa*vdWradsa/3.0))*aux*aux*aux;
      //now loop over grid points
      for (size_t igrpnt = 0; igrpnt < ngridpoints; ++igrpnt) {
        //get the grid point
        rgp[0] = rAB[0] + vdWradsa*Grid(1,igrpnt + 1);
        rgp[1] = rAB[1] + vdWradsa*Grid(2,igrpnt + 1);
        rgp[2] = rAB[2] + vdWradsa*Grid(3,igrpnt + 1);
        sasapoint = 1.0;
        //get atomic surface at grid point
        for (size_t ibtm = 0; ibtm < Natoms; ++ibtm) {
          idbtm = ListShortRange(idatm + 1,ibtm + 1) - 1;
          if (idbtm < 0) {break;}
          tB2 = 0.0;
          for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
            tB[idcoord] = rgp[idcoord] - geom(idbtm + 1,idcoord + 1)*dist_Angstrom2aum1;
            tB2 += tB[idcoord]*tB[idcoord];
          }
          vdWsa = (vdWaalsRadii(atoms[idbtm]) + proberad)*dist_Angstrom2aum1;
          vdWsapw = vdWsa + ww;
          vdWsamw = vdWsa - ww;
          if (tB2 < vdWsapw*vdWsapw) {
            if (tB2 <= vdWsamw*vdWsamw) {
              sasapoint = 0.0;
              break;
            }
            else {
              vdWsapw = sqrt(tB2);
              aux = vdWsapw - vdWsa;
              vdWsamw = ah3*aux*aux;
              tB2 = ah0 + (ah1 + vdWsamw)*aux;
              sasapoint *= tB2;                   //accumulate molecular surface
            }
          }
        }
        if(sasapoint > threshzero) {
          aux = GridWeights[igrpnt]*wrp*sasapoint;                    //weight of numerical quadrature
          sasaA += aux;                                               //accumulate to surface area
        }
      }
      smat(idatm + 1,idatm + 1) += 4.0*pi*stparam*sasaA*GammaScale(solvent,atomA);
      for (size_t idbtm = 0; idbtm < idatm; ++idbtm) {
        smat(idatm + 1,idbtm + 1) = -0.5*charges[idatm]*(BornMat(idatm + 1,idbtm + 1) + BornMat(idbtm + 1,idatm + 1))*charges[idbtm];
        smat(idbtm + 1,idatm + 1) = smat(idatm + 1,idbtm + 1);
      }
    }
    return smat;
  }
  //dispersion related
  virtual void AtomicPolarizabilities(std::vector<double> & polarizabilities, const std::vector<double> & AtmCharge) {}
  virtual void TotalDispersion(double & mc6, double & mc8, const std::vector<double> & AtmCharge) {}
  virtual void TotalPolarizability(double & mpol, const std::vector<double> & AtmCharge) {}
};

#endif //_QC_METHODS_
