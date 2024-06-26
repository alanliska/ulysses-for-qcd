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

#ifndef _MNDO_
#define _MNDO_
#include "QC.hpp"
#include "Intermolecular.hpp"
#include "atoms/AtomPackage.hpp"
#include "math/ConvergenceAccelaration.hpp"
#include "math/SolverPackage.hpp"
#include "math/FunctionPackage.hpp"
#include "parameters/EISOLpar.hpp"
#include "parameters/ZeroOvpar.hpp"
#include "basissets/2ElectronDewar.hpp"

//descrition:
//NDDO methods with sp bases

class MNDO: public QCbasis {
  //this is the implementation of Dewar's MNDO
  //M. J. S. Dewar, W. Thiel, J. Am. Chem. Soc., 99(15), 4899, 1977
  //M. J. S. Dewar, W. Thiel, Theor. Chim. Acta (Berl.), 46, 89, 1977
  //bX quantities are "barred" tensors, which are used only for the open-shell case
protected:
  matrixE gammaSE;
  matrixE SProt;
  matrixE PProt;
  matrixE SDrot;
  matrixE PDrot;
  matrixE DDrot;
  matrixE olddens;
  matrixE bolddens;
  matrixE Bij;
  matrixE AUXmat;
  matrixE integrals;
  matrixE iBlockdR;
  matrixE iBlockdt;
  matrixE iBlockdp;
  matrixE iBlockdRVC;
  matrixE iBlockdtVC;
  matrixE iBlockdpVC;
  matrixE iBlockdRVD;
  matrixE iBlockdtVD;
  matrixE iBlockdpVD;
  std::vector<std::vector<int> > chg1;
  std::vector<std::vector<int> > chg2;
  std::vector<double> AUXvec;
  std::vector<double> AUXvec2;
  std::vector<double> D;
  std::vector<double> Dd;
  std::vector<double> KA;
  std::vector<double> KB;
  std::vector<double> LA;
  std::vector<double> LB;
  std::vector<double> MA;
  std::vector<double> MB;
  std::vector<double> intn;
  std::vector<double> enuc_dR;
  std::vector<double> enuc_dR2;
  std::vector<double> chg;
  std::vector<double> p1;
  std::vector<double> bp1;
  std::vector<int> pos;
  std::vector<int> positions_pp;
  std::vector<size_t> nintegrals;
  std::vector<size_t> atomslocal;
  std::vector<size_t> ipos;
  std::vector<size_t> AOs_local;
  std::string corecorr;                 //this is a control variable to tell which empirical correction should be used for core-core terms
  size_t DIISstatus;
  bool doDIIS;
  bool doRCA;
  bool RCAstatus;
  bool RCAdone;
  bool doPDDG;
  bool oscillationdamp;
  int sizepos_pp;
  int print;
  int ncoreelectrons;
  double oscillationdampthresh;
  double natocc;
  double factorA;
  double factorB;
  double factor;
  double sumA;
  double sumB;
  double deltaA;
  double dval;
  double odval;
  double dampfactor;
public:
  MNDO(BSet _bset, Molecule _mol, const std::string & _openclosed = "0", std::string _corecorrection = "0"): QCbasis(_bset,_mol) {
    openclosed = _openclosed;
    corecorr = _corecorrection;
    pos.resize(7);
    D.resize(4);
    intn.resize(10);
    sizepos_pp = 3;
    positions_pp.resize(sizepos_pp);
    positions_pp[0] = 5;
    positions_pp[1] = 8;
    positions_pp[2] = 10;
    nintegrals.resize(2);
    atomslocal.resize(2);
    AOs_local.resize(2);
    ipos.resize(2);
    chg.resize(2);
    doPDDG = false;
    oscillationdamp = false;
    oscillationdampthresh = 0.05;
  }
  ~MNDO() {}
  //getters
  virtual std::string Type() {return "MNDO";}
  virtual std::string TypeOverlap() {return "Orthogonal";}
  double getHeatFormation() {return EnthForm(etot);}
  virtual bool AtomWithDOrbitals(size_t atomicnr) {
    //function that determines whether atom has d orbitals
    return false;
  }
  bool OscillationDamping() {return oscillationdamp;}
  double ThresholdOscillationDamping() {return oscillationdampthresh;}
  //setters
  void setOscillationDamping(bool newoscdamp) {oscillationdamp = newoscdamp;}
  void setThresholdOscillationDamping(double newthresh) {oscillationdampthresh = newthresh;}
  //Calculation
  void Calculate(int _print = 1, size_t maxiter = 200, bool _DIIS = true, bool _RCA = true, double d2threshold = 1600.0) {
    //initialize necessary quantities
    matrixE geometry = mol.Geometry();
    print = _print;
    DIISstatus = 0;
    doDIIS = _DIIS;
    doRCA = _RCA;
    RCAdone = false;
    RCAstatus = false;
    double energy = 0.0;
    double oldenergy;
    double denconv;
    double aux;
    //required for DIIS; the dimensions are variable, so it is not worth to save these globally for the object
    Bij.resize(1,1);
    Bij(1,1) = 0.0;
    std::vector<matrixE> vError;        //if not declared here, these must be resized, which is the same
    std::vector<matrixE> vFock;         //if not declared here, these must be resized, which is the same
    std::vector<matrixE> vbFock;        //if not declared here, these must be resized, which is the same
    std::vector<matrixE> vDens;         //if not declared here, these must be resized, which is the same
    std::vector<matrixE> vbDens;        //if not declared here, these must be resized, which is the same
    //check whether calculation is possible
    checkAtoms();
    //determine whether system is closed- or open-shell
    AOs = basis.AtomNAOs(atoms);
    OpenClosed();
    if (openclosed != "0") {
      //in this case we want to override what was previously determined by OpenClosed()
      if ((openclosed == "open")||(openclosed == "UHF")||(openclosed == "uhf")) {shell = "open";}
      else if ((openclosed == "closed")||(openclosed == "RHF")||(openclosed == "rhf")) {shell = "closed";}
      else {throw("ERROR: MNDOd.hpp: MNDOd: Calculate(): calculation type not possible to determine");}
    }
    //get occupation vector
    Occ();
    if (print > 0) {
      std::cout << shell << "-shell calculation" << std::endl;
      std::cout << "occupation vector " << occupancy.size() << "/" << NAOs << std::endl;
      for (size_t id = 0; id < occupancy.size(); ++id) {
        std::cout << occupancy[id] << " ";
      }
      std::cout << std::endl;
      if (shell == "open") {
        std::cout << "occupation vector (bar) " << boccupancy.size() << "/" << NAOs << std::endl;
        for (size_t id = 0; id < boccupancy.size(); ++id) {
          std::cout << boccupancy[id] << " ";
        }
        std::cout << std::endl;
      }
    }
    se_integrals(geometry);
    sao = basis.SAO();
    //get Hcore
    calcHcore(sao,geometry);
    //initiliaze Fock Matrix with the density initialization
    if (restart < 1) {
      dens.resize(NAOs,NAOs);
      dens.zero();
      olddens = dens;
      if (shell == "open") {
        bdens = dens;
        bolddens = bdens;
      }
      initDens(geometry);
      calcFock(sao);
      calcMOs();
      olddens = dens;
      Dens(dens,occupancy,CMO,DensAUX);
      if (shell == "open") {
        bolddens = bdens;
        Dens(bdens,boccupancy,bCMO,DensAUX);
      }
      if (restart == 0) {restart = 1;}
    }
    oldenergy = Energy();
    if (print > 0) {std::cout << "starting energy: " << oldenergy << std::endl;}
    vDens.push_back(dens);
    if (shell == "open") {vbDens.push_back(bdens);}
    //SCF
    for (size_t iter = 1; iter < maxiter; ++iter) {
      if (print  > 0) {std::cout << "iteration " << iter << std::endl;}
      //get new MOs
      calcFock(sao);
      ShiftFock(shell == "open",1.0e-6);
      if ((!doDIIS)||(DIISstatus == 0)) {
        if (iter > startpseudodiag) {
          PseudoDiagonalization(CMO,Fock,EMOs,occupancy.size());
          if (shell == "open") {PseudoDiagonalization(bCMO,bFock,bEMOs,boccupancy.size());}
        }
        else {calcMOs();}
      }
      //save old density
      olddens = dens;
      if (shell == "open") {bolddens = bdens;}
      //get new energy
      energy = Energy();
      if (print > 0) {
        std::cout << "CMO " << std::endl;
        CMO.Print();
        if (shell == "open") {
          std::cout << "CMO (bar) " << std::endl;
          bCMO.Print();
        }
        if (print > 1) {
          std::cout << "density " << std::endl;
          dens.Print();
          if (shell == "open") {
            std::cout << "density (bar) " << std::endl;
            bdens.Print();
          }
        }
        std::cout << "Energy = " << energy << std::endl;
      }
      //DIIS
      if (doDIIS) {
        //preparing DIIS
        vFock.push_back(Fock);
        if (shell == "open") {vbFock.push_back(bFock);}
        //DIIS step
        DIISse(vError,vFock,vbFock,dens,bdens,Bij,occupancy,boccupancy,DIISstatus,AUXmat,AUXvec,AUXvec2,DensAUX);
        if (DIISstatus == 2) {
          if (print > 0) {std::cout << "DIIS convergence after " << iter << " iterations" << std::endl;}
          break;
        }
        if ((RCAstatus)&&(DIISstatus == 0)) {
          //If DIIS still does not start converging, then RCA
          RCA(vFock,vbFock,vDens,vbDens,occupancy,boccupancy,energy - enuc,oldenergy - enuc,AUXmat,AUXvec,DensAUX,1.e-7);
          dens = vDens[1]; 
          vDens.clear();
          vDens.push_back(dens);
          if (shell == "open") {
            bdens = vbDens[1]; 
            vbDens.clear();
            vbDens.push_back(bdens);
          }
          RCAdone = true;
        }
      }
      //assume that in open-shell cases RCA is on or off simultaneously for both alpha and beta cases
      if (((!doDIIS)||(DIISstatus == 0))&&(!(RCAstatus))) {
        //get new density
        Dens(dens,occupancy,CMO,DensAUX);
        if (doDIIS) {vDens.push_back(dens);}
        if ((vDens.size() > 1)&&(doRCA)) {RCAstatus = true;}
        if (shell == "open") {
          Dens(bdens,boccupancy,bCMO,DensAUX);
          if (doDIIS) {vbDens.push_back(bdens);}
          if ((vbDens.size() > 1)&&(doRCA)&&(!(RCAstatus))) {RCAstatus = true;}
        }
        if (print > 0) {
          std::cout << "density" << std::endl;
          dens.Print();
          std::cout << "Energy = " << energy << std::endl;
          std::cout << "orbital energies: ";
          for (size_t id = 0; id < EMOs.size(); ++id) {
            std::cout << EMOs[id] << " ";
          }
          std::cout << std::endl;
          if (shell == "open") {
            std::cout << "density (bar)" << std::endl;
            bdens.Print();
            std::cout << "orbital energies (bar): ";
            for (size_t id = 0; id < bEMOs.size(); ++id) {
              std::cout << bEMOs[id] << " ";
            }
            std::cout << std::endl;
          }
        }
        //check for convergence on density
        if (shell == "open") {
          denconv = 0.0;
          for (size_t idc = 0; idc < NAOs; ++idc) {
            for (size_t idr = 0; idr < NAOs; ++idr) {
              aux = dens(idc + 1, idr + 1) + bdens(idc + 1, idr + 1) - olddens(idc + 1, idr + 1) - bolddens(idc + 1, idr + 1);
              denconv += aux*aux;
            }
          }
        }
        else {
          denconv = 0.0;
          for (size_t idc = 0; idc < NAOs; ++idc) {
            for (size_t idr = 0; idr < NAOs; ++idr) {
              aux = dens(idc + 1, idr + 1) - olddens(idc + 1, idr + 1);
              denconv += aux*aux;
            }
          }
        }
        if (denconv < thresh_dens*thresh_dens) {
          if (print > 0) {std::cout << "convergence in density after " << iter << " iterations" << std::endl;}
          break;
        }
      }
      if ((oscillationdamp)&&((iter + 1)%3 == 0)) {
        if (p1.size() != NAOs) {p1.resize(NAOs);}
        DampOscillations(iter,dens,olddens,p1);
        if (shell == "open") {
          if (bp1.size() != NAOs) {bp1.resize(NAOs);}
          DampOscillations(iter,bdens,bolddens,bp1);
        }
      }
      //check for convergence on energy
      if (fabs(oldenergy - energy) < thresh_en) {
        if (print > 0) {std::cout << "convergence in energy after " << iter << " iterations" << std::endl;}
        break;
      }
      oldenergy = energy;
      if ((RCAdone)&&(RCAstatus)) {
        RCAdone = false;
        RCAstatus = false;
      }
    }
    //almost all the techniques to avoid directly solving the SCF require calculating the MOs at the end;
    //only the "normal" way does not require this, but that should have ampered convergence as well, 
    //so one last matrix diagonalization does not affect the calculation times by much
    calcMOs();
    dens = olddens;
    if (shell == "open") {bdens = bolddens;}
    if (corecorr != "0"){
      if ((corecorr == "D3H4X")||(corecorr == "d3h4x")) {corecorrection = D3H4X(atoms,geometry,this->Type());}
      else if ((corecorr == "D3H+")||(corecorr == "d3h+")) {corecorrection = D3Hplus(atoms,geometry,this->Type(),true);}
    }
    if (print > 0) {
      std::cout << "Optimized MOs" << std::endl;
      CMO.Print();
      std::cout << "Optimized orbital energies: ";
      for (size_t id = 0; id < EMOs.size(); ++id) {
        std::cout << EMOs[id] << " ";
      }
      std::cout << std::endl;
      if (shell == "open") {
        std::cout << "Optimized MOs (bar)" << std::endl;
        bCMO.Print();
        std::cout << "Optimized orbital energies (bar): ";
        for (size_t id = 0; id < bEMOs.size(); ++id) {
          std::cout << bEMOs[id] << " ";
        }
        std::cout << std::endl;
      }
      std::cout << "core correction: " << corecorrection << std::endl;
      std::cout << "Energy with correction = " << energy + corecorrection << std::endl;
      std::cout << "Energy = " << energy << std::endl;
      std::cout << "Nuclear Repulsion = " << enuc << std::endl;
      std::cout << "Electronic Energy = " << energy - enuc << std::endl;
      std::cout << "Enthalpy Formation = " << EnthForm(energy) << std::endl;
    }
    etot = energy;
  }
  //calculate MO occupancy
  void Occ() {
    //function to calculate the occupancy vector for a system
    occupancy.clear();
    boccupancy.clear();
    size_t Nelec = size_t(mol.Nelectrons());
    size_t Nelec_aux = mol.Multiplicity() - 1;             //total number of unpaired electrons
    if (print > 0) {std::cout << "multiplicity " << mol.Multiplicity() << std::endl;}
    Nelec -= Nelec_aux + ncoreelectrons;
    Nelec /= 2;
    for (size_t idx = 0; idx < Nelec; ++idx) {
      occupancy.push_back(2.0);
    }
    for (size_t idx = 0; idx < Nelec_aux; ++idx) {
      occupancy.push_back(1.0);
    }
    if (shell == "open") {
      for (size_t idx = 0; idx < occupancy.size(); ++idx) {
        if (int(occupancy[idx]) == 2) {
          occupancy[idx] -= 1.0;
          boccupancy.push_back(1.0);
        }
      }
    }
  }
  //calculate energy
  double Energy() {
    double energy = enuc;
    //other terms
    iBlockdR = Hcore + Fock;
    energy += 0.5*trace(iBlockdR,dens);
    if (shell == "open") {
      iBlockdR = Hcore + bFock;
      energy += 0.5*trace(iBlockdR,bdens);
    }
    return energy;
  }
  void calcEnergy() {etot = this->Energy();}
  double EnthForm(double energy) {
    //function calculating heats of formation
    double enth = energy;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      enth += ElementFormExp(atoms[iatm]) - ElementFormTheo(atoms[iatm]);
    }
    return enth*au2kcalmol;
  }
  //functions required to build the Hamiltonian
  int ShiftAO(int iorb, bool transmetal) {
    //just a shift function to simplify code
    int im = iorb;
    if (transmetal) {
      if (iorb > 4) {im = iorb - 5;}
      else if (iorb < 5) {im = iorb + 4;}
    }
    return im;
  }
  void calcHcore(matrixE & sao, matrixE & geom) {
    //function that calculates the core Hamiltonian
    Hcore.resize(NAOs,NAOs);
    Hcore.zero();
    int counter = 0;
    int counter2 = 0;
    int counter3 = 0;
    int pgamma;
    int L = 0;
    int Lp = 0;
    int im;
    int in;
    bool transitionmetalA;
    bool transitionmetalB;
    double RAAmax;
    double Raux;
    double RAB;
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      transitionmetalA = TransitionMetal(atoms[idatm]);
      L = 2*int(transitionmetalA);
      rAB[0] = geom(idatm + 1,1);
      rAB[1] = geom(idatm + 1,2);
      rAB[2] = geom(idatm + 1,3);
      RAAmax = ZeroOverlap(atoms[idatm]);
      for (size_t imu = 0; imu < AOs[idatm]; ++imu) {
        Hcore(counter + imu + 1, counter + imu + 1) = UlX(atoms[idatm],L);
        im = ShiftAO(imu,transitionmetalA);
        Hcore(counter + imu + 1, counter + imu + 1) += VAB(counter2 + posgamma(im,im) + 1,1);
        for (size_t inu = imu + 1; inu < AOs[idatm]; ++inu) {
          //determine the position in Vmunu required
          in = ShiftAO(inu,transitionmetalA);
          pgamma = posgamma(im,in);
          Hcore(counter + imu + 1, counter + inu + 1) += VAB(counter2 + pgamma + 1,1);
          Hcore(counter + inu + 1, counter + imu + 1) += VAB(counter2 + pgamma + 1,1);
        }
        counter3 = 0;
        for (size_t idbtm = 0; idbtm < idatm; ++idbtm) {
          RAB = 0.0;
          for (size_t idcoord = 0; idcoord < 3; ++idcoord){
            Raux = rAB[idcoord] - geom(idbtm + 1,idcoord + 1);
            RAB += Raux*Raux;
          }
          Raux = fmax(RAAmax,ZeroOverlap(atoms[idbtm]));         //this here is an approximation, but works conservatively well
          if (RAB < Raux*Raux) {
            transitionmetalB = TransitionMetal(atoms[idbtm]);
            Lp = 2*int(transitionmetalB);
            for (size_t jorb = 0; jorb < AOs[idbtm]; ++jorb) {
              Hcore(counter + imu + 1, counter3 + jorb + 1) = 0.5*(betaA0(atoms[idatm],L) + betaA0(atoms[idbtm],Lp))*sao(counter + imu + 1, counter3 + jorb + 1);
              Hcore(counter3 + jorb + 1,counter + imu + 1) = Hcore(counter + imu + 1, counter3 + jorb + 1);
              if (transitionmetalB) {
                if (jorb == 4) {Lp = 0;}
                if (jorb == 5) {++Lp;}
              }
              else {
                if ((jorb == 0)||(jorb == 3)) {++Lp;}
              }
            }
          }
          counter3 += AOs[idbtm];
        }
        if (transitionmetalA) {
          if (imu == 4) {L = 0;}
          if (imu == 5) {++L;}
        }
        else {
          if ((imu == 0)||(imu == 3)) {++L;}
        }
      }
      counter += AOs[idatm];
      counter2 += AOs[idatm]*(AOs[idatm] + 1)/2;
    }
    if (print > 0) {
      std::cout << "Hcore matrix" << std::endl;
      Hcore.Print();
    }
  }
  void HcoreCD_dX(matrixE & SAO_dX, size_t atmnrC, size_t atmnrD, size_t naoC, size_t naoD) {
    //function that calculates the off-diagonal block of the core Hamiltonian or any of its derivatives, involving the pair of atoms C and D
    //input is a matrix containing already the overlap matrix SAO (or its derivatives), which is modified to get the respective Hcore block
    //note that this applies only to off-diagonal blocks in terms of atoms, i.e. C != D
    bool transitionmetalC = TransitionMetal(atmnrC);
    bool transitionmetalD = TransitionMetal(atmnrD);
    int L = 2*int(transitionmetalC);
    int Lp = 0;
    for (size_t iorb = 0; iorb < naoC; ++iorb) {
      Lp = 2*int(transitionmetalD);
      for (size_t jorb = 0; jorb < naoD; ++jorb) {
        SAO_dX(iorb + 1,jorb + 1) *= 0.5*(betaA0(atmnrC,L) + betaA0(atmnrD,Lp));
        if (transitionmetalD) {
          if (jorb == 4) {Lp = 0;}
          if (jorb == 5) {++Lp;}
        }
        else {
          if ((jorb == 0)||(jorb == 3)) {++Lp;}
        }
      }
      if (transitionmetalC) {
        if (iorb == 4) {L = 0;}
        if (iorb == 5) {++L;}
      }
      else {
        if ((iorb == 0)||(iorb == 3)) {++L;}
      }
    }
  }
  void initDens(matrixE & geometry) {
    //function that calculates starting values for the density matrix according to what is done in MOPAC
    //only diagonal elements changed
    double dens_el;
    double bdens_el;
    double factorSym = 1.0;
    double factorShell = 1.0;
    bool transitionmetal;
    int cntA = 0;
    //generate asymmetry and open-shell conditions
    if (shell == "open") {
      factorShell = double(occupancy.size())/double(occupancy.size() + boccupancy.size());
      if (occupancy.size() != boccupancy.size()) {factorSym = 1.1;}
    }
    GoedeckerCharges(QAtoms,atoms,geometry,mol.Charge());
    for (size_t idA = 0; idA < Natoms; ++idA) {
      QAtoms[idA] *= -1.0;
      QAtoms[idA] += atoms[idA] - CoreCharge[idA];
      transitionmetal = (((atoms[idA] > 20)&&(atoms[idA] < 30))||((atoms[idA] > 38)&&(atoms[idA] < 48))||((atoms[idA] > 56)&&(atoms[idA] < 80)));
      for (size_t imu = 0; imu < AOs[idA]; ++imu) {
        if (!transitionmetal) {
          if (imu > 3) {break;}
          dens_el = factorShell*factorSym*QAtoms[idA]/4.0;
          bdens_el = (1.0 - factorShell)*QAtoms[idA]/(4.0*factorSym);
        }
        else {
          dens_el = factorShell*factorSym*QAtoms[idA]/9.0;
          bdens_el = (1.0 - factorShell)*QAtoms[idA]/(9.0*factorSym);
        }
        if (CoreCharge[idA] == 0) {
          dens_el *= 4.0;
          bdens_el *= 4.0;
        }
        dens(cntA + imu + 1,cntA + imu + 1) = dens_el;
        if (shell == "open") {bdens(cntA + imu + 1,cntA + imu + 1) = bdens_el;}
      }
      cntA += AOs[idA];
    }
  }
  void calcFock(matrixE & SAO) {
    //function that calculates the Fock matrix
    Fock = Hcore;
    if (shell == "open") {bFock = Hcore;}
    //atom position control for Fock and density
    int cntA = 0;
    int cntB = 0;
    //atom position control for gamma matrix
    int gposA = 0;
    int gposB = 0;
    int im;
    int in;
    int il;
    int is;
    bool transitionmetalA;
    bool transitionmetalB;
    for (size_t idA = 0; idA < Natoms; ++idA) {
      transitionmetalA = (((atoms[idA] > 20)&&(atoms[idA] < 30))||((atoms[idA] > 38)&&(atoms[idA] < 48))||((atoms[idA] > 56)&&(atoms[idA] < 80)));
      for (size_t imu = 0; imu < AOs[idA]; ++imu) {
        //get position in gamma matrix
        im = ShiftAO(imu,transitionmetalA);
        pos[0] = posgamma(im,im);
        for (size_t inu = 0; inu < AOs[idA]; ++inu) {
          in = ShiftAO(inu,transitionmetalA);
          pos[1] = posgamma(im,in);
          for (size_t ilambda = 0; ilambda < AOs[idA]; ++ ilambda) {
            il = ShiftAO(ilambda,transitionmetalA);
            pos[2] = posgamma(in,il);
            pos[3] = posgamma(im,il);
            Fock(cntA + imu + 1,cntA + imu + 1) += dens(cntA + inu + 1,cntA + ilambda + 1)*(gammaSE(gposA + pos[0] + 1,gposA + pos[2] + 1) - 0.5*gammaSE(gposA + pos[1] + 1,gposA + pos[3] + 1));
            if (shell == "open") {
              Fock(cntA + imu + 1,cntA + imu + 1) += bdens(cntA + inu + 1,cntA + ilambda + 1)*gammaSE(gposA + pos[0] + 1,gposA + pos[2] + 1) - 0.5*dens(cntA + inu + 1,cntA + ilambda + 1)*gammaSE(gposA + pos[1] + 1,gposA + pos[3] + 1);
              bFock(cntA + imu + 1,cntA + imu + 1) += bdens(cntA + inu + 1,cntA + ilambda + 1)*(gammaSE(gposA + pos[0] + 1,gposA + pos[2] + 1) - gammaSE(gposA + pos[1] + 1,gposA + pos[3] + 1));
              bFock(cntA + imu + 1,cntA + imu + 1) += dens(cntA + inu + 1,cntA + ilambda + 1)*gammaSE(gposA + pos[0] + 1,gposA + pos[2] + 1);
            }
            if (imu != inu) {
              for (size_t isigma = 0; isigma < AOs[idA]; ++ isigma) {
                is = ShiftAO(isigma,transitionmetalA);
                pos[2] = posgamma(il,is);
                pos[4] = posgamma(is,in);
                Fock(cntA + imu + 1,cntA + inu + 1) += dens(cntA + ilambda + 1,cntA + isigma + 1)*(gammaSE(gposA + pos[1] + 1,gposA + pos[2] + 1) - 0.5*gammaSE(gposA + pos[3] + 1,gposA + pos[4] + 1));
                if (shell == "open") {
                  Fock(cntA + imu + 1,cntA + inu + 1) += bdens(cntA + ilambda + 1,cntA + isigma + 1)*gammaSE(gposA + pos[1] + 1,gposA + pos[2] + 1);
                  Fock(cntA + imu + 1,cntA + inu + 1) -= 0.5*dens(cntA + ilambda + 1,cntA + isigma + 1)*gammaSE(gposA + pos[3] + 1,gposA + pos[4] + 1);
                  bFock(cntA + imu + 1,cntA + inu + 1) += (bdens(cntA + ilambda + 1,cntA + isigma + 1) + dens(cntA + ilambda + 1,cntA + isigma + 1))*gammaSE(gposA + pos[1] + 1,gposA + pos[2] + 1);
                  bFock(cntA + imu + 1,cntA + inu + 1) -= bdens(cntA + ilambda + 1,cntA + isigma + 1)*gammaSE(gposA + pos[3] + 1,gposA + pos[4] + 1);
                }
              }
            }
          }
          cntB = 0;
          gposB = 0;
          for (size_t idB = 0; idB < idA; ++idB) {
            transitionmetalB = (((atoms[idB] > 20)&&(atoms[idB] < 30))||((atoms[idB] > 38)&&(atoms[idB] < 48))||((atoms[idB] > 56)&&(atoms[idB] < 80)));
            for (size_t ilambda = 0; ilambda < AOs[idB]; ++ilambda) {
              il = ShiftAO(ilambda,transitionmetalB);
              for (size_t isigma = 0; isigma < AOs[idB]; ++isigma) {
                is = ShiftAO(isigma,transitionmetalB);
                pos[2] = posgamma(il,is);
                Fock(cntA + imu + 1,cntB + ilambda + 1) -= 0.5*dens(cntA + inu + 1,cntB + isigma + 1)*gammaSE(gposA + pos[1] + 1,gposB + pos[2] + 1);
                if (shell == "open") {
                  Fock(cntA + imu + 1,cntB + ilambda + 1) -= 0.5*dens(cntA + inu + 1,cntB + isigma + 1)*gammaSE(gposA + pos[1] + 1,gposB + pos[2] + 1);
                  bFock(cntA + imu + 1,cntB + ilambda + 1) -= bdens(cntA + inu + 1,cntB + isigma + 1)*gammaSE(gposA + pos[1] + 1,gposB + pos[2] + 1);
                }
                Fock(cntB + ilambda + 1,cntA + imu + 1) -= 0.5*dens(cntB + isigma + 1,cntA + inu + 1)*gammaSE(gposB + pos[2] + 1,gposA + pos[1] + 1);
                if (shell == "open") {
                  Fock(cntB + ilambda + 1,cntA + imu + 1) -= 0.5*dens(cntB + isigma + 1,cntA + inu + 1)*gammaSE(gposB + pos[2] + 1,gposA + pos[1] + 1);
                  bFock(cntB + ilambda + 1,cntA + imu + 1) -= bdens(cntB + isigma + 1,cntA + inu + 1)*gammaSE(gposB + pos[2] + 1,gposA + pos[1] + 1);
                }
                Fock(cntA + imu + 1,cntA + inu + 1) += dens(cntB + ilambda + 1,cntB + isigma + 1)*gammaSE(gposA + pos[1] + 1,gposB + pos[2] + 1);
                if (shell == "open") {
                  Fock(cntA + imu + 1,cntA + inu + 1) += bdens(cntB + ilambda + 1,cntB + isigma + 1)*gammaSE(gposA + pos[1] + 1,gposB + pos[2] + 1);
                  bFock(cntA + imu + 1,cntA + inu + 1) += (dens(cntB + ilambda + 1,cntB + isigma + 1) + bdens(cntB + ilambda + 1,cntB + isigma + 1))*gammaSE(gposA + pos[1] + 1,gposB + pos[2] + 1);
                }
                Fock(cntB + ilambda + 1,cntB + isigma + 1) += dens(cntA + imu + 1,cntA + inu + 1)*gammaSE(gposB + pos[2] + 1,gposA + pos[1] + 1);
                if (shell == "open") {
                  Fock(cntB + ilambda + 1,cntB + isigma + 1) += bdens(cntA + imu + 1,cntA + inu + 1)*gammaSE(gposB + pos[2] + 1,gposA + pos[1] + 1);
                  bFock(cntB + ilambda + 1,cntB + isigma + 1) += (dens(cntA + imu + 1,cntA + inu + 1) + bdens(cntA + imu + 1,cntA + inu + 1))*gammaSE(gposB + pos[2] + 1,gposA + pos[1] + 1);
                }
              }
            }
            cntB += AOs[idB];
            gposB += AOs[idB]*(AOs[idB] + 1)/2;
          }
        }
      }
      cntA += AOs[idA];
      gposA += AOs[idA]*(AOs[idA] + 1)/2;
    }
    if (print > 0) {
      std::cout << "Fock matrix" << std::endl;
      Fock.Print();
      if (shell == "open") {
        std::cout << "Fock matrix (bar)" << std::endl;
        bFock.Print();
      }
    }
  }
  void AnalyticalGrad(matrixE & gen, double tolerance = 1.0e-6) {
    //function calculating the gradients of the energy for MNDO methods
    //all final derivatives come for nuclear coordinates in atomic units!
    this->getDens(olddens);                //direct use of density matrix can lead to problems with open-shells
    if (shell == "open") {
      this->getbDens(bolddens);
      AUXmat = olddens - bolddens;
      olddens += bolddens;
      bolddens = AUXmat;
    }
    AUXmat = mol.Geometry();
    AOs = basis.AtomNAOs(atoms);
    Enuclear_dR(enuc_dR);
    int icntC = 0;
    int icntD = 0;
    int icnt = 0;
    int jcnt = 0;
    int kcnt = 0;
    int naoC;
    int naoD;
    int nintC;
    int nintD;
    int mu;
    int nu;
    int lambda;
    int sigma;
    bool transitionmetalC;
    bool transitionmetalD;
    gen.resize(3*Natoms,1);
    gen.zero();
    double enel_dR = 0.0;                            //derivatives of energy with respect to {R}_{CD}
    double enel_dt = 0.0;                            //derivatives of energy with respect to {\theta}_{CD}
    double enel_dp = 0.0;                            //derivatives of energy with respect to {\phi}_{CD}
    double RCD;
    double cost;
    double sint;
    double cosp;
    double sinp;
    double auxD;
    double auxV;
    double dRCDdxD;
    double dRCDdyD;
    double dRCDdzD;
    double dTCDdxD;
    double dTCDdyD;
    double dTCDdzD;
    double dPCDdxD;
    double dPCDdyD;
    double dPCDdzD;
    double chgC;
    double chgD;
    for (size_t idC = 0; idC < Natoms; ++idC) {
      chgC = double(int(atoms[idC]) - int(CoreCharge[idC]));
      naoC = AOs[idC];
      nintC = naoC*(naoC + 1)/2;
      icntD = icntC + naoC;
      transitionmetalC = TransitionMetal(atoms[idC]);
      for (size_t idD = idC + 1; idD < Natoms; ++idD) {
        chgD = double(int(atoms[idD]) - int(CoreCharge[idD]));
        naoD = AOs[idD];
        nintD = naoD*(naoD + 1)/2;
        transitionmetalD = TransitionMetal(atoms[idD]);
        //getting the orientation vector
        rAB[0] = (AUXmat(idD + 1,1) - AUXmat(idC + 1,1))*dist_Angstrom2aum1;                  //Delta x
        rAB[1] = (AUXmat(idD + 1,2) - AUXmat(idC + 1,2))*dist_Angstrom2aum1;                  //Delta y
        rAB[2] = (AUXmat(idD + 1,3) - AUXmat(idC + 1,3))*dist_Angstrom2aum1;                  //Delta z
        //normalizing it
        RCD = sqrt(rAB[0]*rAB[0] + rAB[1]*rAB[1] + rAB[2]*rAB[2]);
        //getting chain rule terms
        dRCDdxD = rAB[0]/RCD;
        dRCDdyD = rAB[1]/RCD;
        dRCDdzD = rAB[2]/RCD;
        dTCDdxD = 0.0;
        dTCDdyD = 0.0;
        dTCDdzD = -sqrt(rAB[0]*rAB[0] + rAB[1]*rAB[1])/(RCD*RCD);
        dPCDdxD = 0.0;
        dPCDdyD = 0.0;
        dPCDdzD = 0.0;
        if (fabs(rAB[0] + rAB[1]) > tolerance) {
          dTCDdyD = rAB[1]*rAB[2]/(sqrt(rAB[0]*rAB[0] + rAB[1]*rAB[1])*RCD*RCD);
          dTCDdxD = rAB[0]*rAB[2]/(sqrt(rAB[0]*rAB[0] + rAB[1]*rAB[1])*RCD*RCD);
          dPCDdyD = rAB[0]/(rAB[0]*rAB[0] + rAB[1]*rAB[1]);
          dPCDdxD = -rAB[1]/(rAB[0]*rAB[0] + rAB[1]*rAB[1]);
        }
        rAB[0] /= RCD;
        rAB[1] /= RCD;
        rAB[2] /= RCD;
        //getting trigonometric functions for rotations
        cost = rAB[2];
        sint = sqrt(1.0 - cost*cost);
        cosp = 1.0;
        sinp = 0.0;
        if (fabs(sint) > tolerance) {
          cosp = rAB[0]/sint;
          sinp = rAB[1]/sint;
        }
        //integrals
        iBlockdR.resize(nintC,nintD);              //resize the blocks containing the integral first-derivatives with respect to internuclear distance
        iBlockdt.resize(nintC,nintD);              //resize the blocks containing the integral first-derivatives with respect to angle theta
        iBlockdp.resize(nintC,nintD);              //resize the blocks containing the integral first-derivatives with respect to angle phi
        iBlockdRVC.resize(nintC,1);                //resize the blocks containing the integral first-derivatives with respect to internuclear distance
        iBlockdtVC.resize(nintC,1);                //resize the blocks containing the integral first-derivatives with respect to angle theta
        iBlockdpVC.resize(nintC,1);                //resize the blocks containing the integral first-derivatives with respect to angle phi
        iBlockdRVD.resize(nintD,1);                //resize the blocks containing the integral first-derivatives with respect to internuclear distance
        iBlockdtVD.resize(nintD,1);                //resize the blocks containing the integral first-derivatives with respect to angle theta
        iBlockdpVD.resize(nintD,1);                //resize the blocks containing the integral first-derivatives with respect to angle phi
        IntegralBlock2C_dR(0,iBlockdR,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);
        IntegralBlock2C_dA(0,iBlockdt,2,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);
        IntegralBlock2C_dA(0,iBlockdp,3,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);
        IntegralBlock2C_dR(1,iBlockdRVC,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);
        IntegralBlock2C_dA(1,iBlockdtVC,2,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);
        IntegralBlock2C_dA(1,iBlockdpVC,3,atoms[idC],atoms[idD],RCD,cost,sint,cosp,sinp);
        IntegralBlock2C_dR(1,iBlockdRVD,atoms[idD],atoms[idC],RCD,-cost,sint,-cosp,-sinp);
        IntegralBlock2C_dA(1,iBlockdtVD,2,atoms[idD],atoms[idC],RCD,-cost,sint,-cosp,-sinp);
        IntegralBlock2C_dA(1,iBlockdpVD,3,atoms[idD],atoms[idC],RCD,-cost,sint,-cosp,-sinp);
        enel_dR = 0.0;
        enel_dt = 0.0;
        enel_dp = 0.0;
        for (size_t idmu = 0; idmu < naoC; ++idmu) {
          mu = ShiftAO(idmu,transitionmetalC);
          for (size_t idnu = idmu; idnu < naoC; ++idnu) {
            nu = ShiftAO(idnu,transitionmetalC);
            jcnt = posgamma(mu,nu);
            auxV = olddens(icntC + idmu + 1,icntC + idnu + 1);
            if (idmu != idnu) {auxV *= 2.0;}
            enel_dR -= auxV*chgD*iBlockdRVC(jcnt + 1,1);
            enel_dt -= auxV*chgD*iBlockdtVC(jcnt + 1,1);
            enel_dp -= auxV*chgD*iBlockdpVC(jcnt + 1,1);
            for (size_t idlambda = 0; idlambda < naoD; ++idlambda) {
              lambda = ShiftAO(idlambda,transitionmetalD);
              for (size_t idsigma = idlambda; idsigma < naoD; ++idsigma) {
                sigma = ShiftAO(idsigma,transitionmetalD);
                kcnt = posgamma(lambda,sigma);
                if ((mu == 0)&&(nu == 0)) {
                  auxV = olddens(icntD + idlambda + 1,icntD + idsigma + 1);
                  if (idlambda != idsigma) {auxV *= 2.0;}
                  enel_dR -= auxV*chgC*iBlockdRVD(kcnt + 1,1);
                  enel_dt += auxV*chgC*iBlockdtVD(kcnt + 1,1);
                  enel_dp -= auxV*chgC*iBlockdpVD(kcnt + 1,1);
                }
                auxD = olddens(icntC + idmu + 1,icntC + idnu + 1)*olddens(icntD + idlambda + 1,icntD + idsigma + 1) - 0.5*olddens(icntC + idmu + 1,icntD + idlambda + 1)*olddens(icntC + idnu + 1,icntD + idsigma + 1);
                if (shell == "open") {auxD -= 0.5*bolddens(icntC + idmu + 1,icntD + idlambda + 1)*bolddens(icntC + idnu + 1,icntD + idsigma + 1);}
                if (idlambda != idsigma) {
                  auxD += olddens(icntC + idmu + 1,icntC + idnu + 1)*olddens(icntD + idsigma + 1,icntD + idlambda + 1) - 0.5*olddens(icntC + idmu + 1,icntD + idsigma + 1)*olddens(icntC + idnu + 1,icntD + idlambda + 1);
                  if (shell == "open") {auxD -= 0.5*bolddens(icntC + idmu + 1,icntD + idsigma + 1)*bolddens(icntC + idnu + 1,icntD + idlambda + 1);}
                }
                if (idmu != idnu) {
                  auxD += olddens(icntC + idnu + 1,icntC + idmu + 1)*olddens(icntD + idlambda + 1,icntD + idsigma + 1) - 0.5*olddens(icntC + idnu + 1,icntD + idlambda + 1)*olddens(icntC + idmu + 1,icntD + idsigma + 1);
                  if (shell == "open") {auxD -= 0.5*bolddens(icntC + idnu + 1,icntD + idlambda + 1)*bolddens(icntC + idmu + 1,icntD + idsigma + 1);}
                }
                if ((idlambda != idsigma)&&(idmu != idnu)) {
                  auxD += olddens(icntC + idnu + 1,icntC + idmu + 1)*olddens(icntD + idsigma + 1,icntD + idlambda + 1) - 0.5*olddens(icntC + idnu + 1,icntD + idsigma + 1)*olddens(icntC + idmu + 1,icntD + idlambda + 1);
                  if (shell == "open") {auxD -= 0.5*bolddens(icntC + idnu + 1,icntD + idsigma + 1)*bolddens(icntC + idmu + 1,icntD + idlambda + 1);}
                }
                enel_dR += auxD*iBlockdR(jcnt + 1,kcnt + 1);
                enel_dt += auxD*iBlockdt(jcnt + 1,kcnt + 1);
                enel_dp += auxD*iBlockdp(jcnt + 1,kcnt + 1);
              }
            }
          }
        }
        //Overlap terms
        basis.Overlap_dXCD(iBlockdR,"R",atoms[idC],atoms[idD],RCD,rAB);
        basis.Overlap_dXCD(iBlockdt,"t",atoms[idC],atoms[idD],RCD,rAB);
        basis.Overlap_dXCD(iBlockdp,"p",atoms[idC],atoms[idD],RCD,rAB);
        //transform to Hcore
        HcoreCD_dX(iBlockdR,atoms[idC],atoms[idD],naoC,naoD);
        HcoreCD_dX(iBlockdt,atoms[idC],atoms[idD],naoC,naoD);
        HcoreCD_dX(iBlockdp,atoms[idC],atoms[idD],naoC,naoD);
        //even though I could solve this with trace, I merge into a single loop
        for (size_t idmu = 0; idmu < naoC; ++idmu) {
          for (size_t idlambda = 0; idlambda < naoD; ++idlambda) {
            enel_dR += 2.0*olddens(icntC + idmu + 1,icntD + idlambda + 1)*iBlockdR(idmu + 1,idlambda + 1);
            enel_dt += 2.0*olddens(icntC + idmu + 1,icntD + idlambda + 1)*iBlockdt(idmu + 1,idlambda + 1);
            enel_dp += 2.0*olddens(icntC + idmu + 1,icntD + idlambda + 1)*iBlockdp(idmu + 1,idlambda + 1);
          }
        }
        gen(3*idC + 1,1) -= dist_Angstrom2aum1*((enel_dR + enuc_dR[icnt])*dRCDdxD + enel_dt*dTCDdxD + enel_dp*dPCDdxD);        //xC
        gen(3*idC + 2,1) -= dist_Angstrom2aum1*((enel_dR + enuc_dR[icnt])*dRCDdyD + enel_dt*dTCDdyD + enel_dp*dPCDdyD);        //yC
        gen(3*idC + 3,1) -= dist_Angstrom2aum1*((enel_dR + enuc_dR[icnt])*dRCDdzD + enel_dt*dTCDdzD + enel_dp*dPCDdzD);        //zC
        gen(3*idD + 1,1) += dist_Angstrom2aum1*((enel_dR + enuc_dR[icnt])*dRCDdxD + enel_dt*dTCDdxD + enel_dp*dPCDdxD);        //xD
        gen(3*idD + 2,1) += dist_Angstrom2aum1*((enel_dR + enuc_dR[icnt])*dRCDdyD + enel_dt*dTCDdyD + enel_dp*dPCDdyD);        //yD
        gen(3*idD + 3,1) += dist_Angstrom2aum1*((enel_dR + enuc_dR[icnt])*dRCDdzD + enel_dt*dTCDdzD + enel_dp*dPCDdzD);        //zD
        ++icnt;
        icntD += naoD;
      }
      icntC += naoC;
    }
    //add the empirical corrections
    if (corecorr != "0"){
      std::vector<double> dispersioncorrection;
      if ((corecorr == "D3H4X")||(corecorr == "d3h4x")) {gD3H4X(dispersioncorrection,atoms,AUXmat,this->Type());}
      else if ((corecorr == "D3H+")||(corecorr == "d3h+")) {gD3Hplus(dispersioncorrection,atoms,AUXmat,this->Type(),true,9000.0,10.5,1.4);}
      for (size_t idx = 0; idx < 3*Natoms; ++idx) {
        gen(idx + 1,1) += dispersioncorrection[idx];
      }
    }
  }
  int AvailableGradients() {return 1;}
  size_t derindex(size_t ic1, size_t ic2) {
    //function returning the index of the derivatives for a pair of internal coordinates ic1,ic2
    size_t index = 0;
    if (ic1 == 0) {
      if (ic2 == 0) {index = 3;}
      else if (ic2 == 1) {index = 6;}
      else if (ic2 == 2) {index = 7;}
    }
    else if (ic1 == 1) {
      if (ic2 == 0) {index = 9;}
      else if (ic2 == 1) {index = 4;}
      else if (ic2 == 2) {index = 8;}
    }
    else if (ic1 == 2) {
      if (ic2 == 0) {index = 10;}
      else if (ic2 == 1) {index = 11;}
      else if (ic2 == 2) {index = 5;}
    }
    return index;
  }
  void IntegralBlock1C(matrixE & block, int atom) {
    //function calculating a block of one-center integrals
    block.zero();
    block(1,1) = eri1Center(atom,0,0);
    size_t nrows = block.rows();
    if (nrows > 1) {
      intn[0] = eri1Center(atom,1,1);
      for (size_t idsp = 2; idsp < 5; ++idsp) {
        block(idsp,idsp) = intn[0];                                //<spx|spx>,<spy|spy>,<spz|spz>
      }
      intn[0] = eri1Center(atom,2,2);
      intn[1] = eri1Center(atom,0,2);
      intn[2] = eri1Center(atom,2,-2);
      for (size_t idx = 0; idx < sizepos_pp; ++idx) {
        block(positions_pp[idx],positions_pp[idx]) = intn[0];      //<pxpx|pxpx>,<pypy|pypy>,<pzpz|pzpz>
        block(1,positions_pp[idx]) = intn[1];                      //<ss|pxpx>,<ss|pypy>,<ss|pzpz>
        for (size_t idy = idx + 1; idy < sizepos_pp; ++idy) {
          block(positions_pp[idx],positions_pp[idy]) = intn[2];    //<pxpx|pypy>,<pxpx|pzpz>,<pypy|pzpz>
        }
      }
      intn[0] = 0.5*(eri1Center(atom,2,2) - eri1Center(atom,2,-2));
      block(6,6) = intn[0];                                        //<pxpy|pxpy>
      block(7,7) = intn[0];                                        //<pxpz|pxpz>
      block(9,9) = intn[0];                                        //<pypz|pypz>
      if (nrows > 10) {                                    //then d orbitals
        //basic quantities
        intn[0] = SlaterCondonRadialIntegral(atom,0,1,6);        //F0sd
        intn[1] = SlaterCondonRadialIntegral(atom,0,3,6);        //F0pd
        intn[2] = SlaterCondonRadialIntegral(atom,0,6,6);        //F0dd
        intn[3] = SlaterCondonRadialIntegral(atom,1,5,5);        //G1pd
        intn[4] = SlaterCondonRadialIntegral(atom,1,2,5);        //r1sppd
        intn[5] = SlaterCondonRadialIntegral(atom,2,4,4);        //G2sd
        intn[6] = SlaterCondonRadialIntegral(atom,2,3,6);        //F2pd
        intn[7] = SlaterCondonRadialIntegral(atom,2,6,6);        //F2dd
        intn[8] = SlaterCondonRadialIntegral(atom,2,3,4);        //r2sdpp
        intn[9] = SlaterCondonRadialIntegral(atom,2,4,6);        //r2sddd
        intn[10] = SlaterCondonRadialIntegral(atom,3,5,5);       //G3pd
        intn[11] = SlaterCondonRadialIntegral(atom,4,6,6);       //F4dd
        for (size_t idp1 = 1; idp1 < 4; ++idp1) {
          for (size_t idp2 = 1; idp2 < 4; ++idp2) {
            if (idp1 == idp2) {
              if (idp1 < 3) {
                block(1 + idp1,5*(2 + idp2) + 3) = -sqrt(1.0/45.0)*intn[4];                                           //(spx|pxdz2),(spy|pydz2)
                block(1 + idp1,5*(2 + idp2) + 5) = pow(-1.0,idp1 - 1)*sqrt(3.0/45.0)*intn[4];                         //(spx|pxdx2-y2),(spy|pydx2-y2)
              }
            }
            else if ((idp1 < 3)&&(idp2 == 3)) {
              block(1 + idp2,5*(2 + idp1) + idp1) = sqrt(3.0/45.0)*intn[4];                                           //(spz|pxdxz),(spz|pydyz)
              block(1 + idp1,25 + idp1) = sqrt(3.0/45.0)*intn[4];                                                     //(spx|pzdxz),(spy|pydxz)
            }
          }
        }
        int cnt1 = 1;
        int cnt2 = 1;
        for (size_t ippb1 = 1; ippb1 < 4; ++ippb1) {
          for (size_t ippb2 = ippb1; ippb2 < 4; ++ippb2, ++cnt1) {
            cnt2 = 1;
            for (size_t iddk1 = 1; iddk1 < 6; ++iddk1) {
              if ((ippb1 == ippb2)&&(iddk1 == 3)) {
                if (ippb1 < 3) {block(4 + cnt1,10 + iddk1) = -intn[8]/sqrt(125.0);}                                   //(pxpx|sdz2),(pypy|sdz2)
              }
              else if ((ippb1 == ippb2)&&(ippb1 < 3)&&(iddk1 == 5)) {block(4 + cnt1,10 + iddk1) = pow(-1.0,ippb1 - 1)*sqrt(3.0/125.0)*intn[8];} //(pxpx|sdx2-y2),(pypy|sdx2-y2)
              else if ((ippb1 == iddk1)&&(ippb2 == 3)) {block(4 + cnt1,10 + iddk1) = sqrt(3.0/125.0)*intn[8];}        //(pxpz|sdxz),(pypz|sdyz)
              for (size_t iddk2 = iddk1; iddk2 < 6; ++iddk2, ++cnt2) {
                if ((ippb1 == ippb2)&&(iddk1 == iddk2)) {
                  if (ippb1 == 3) {
                    if (iddk1 < 3) {intn[12] = 2.0;}
                    else if (iddk1 == 3) {intn[12] = 4.0;}
                    else if (iddk1 > 3) {intn[12] = -4.0;}
                  }
                  else if (ippb1 < 3) {
                    if ((iddk1 > 3)||(ippb1 == iddk1)) {intn[12] = 2.0;}
                    else if ((iddk1 < 3)&&(ippb1 != iddk1)) {intn[12] = -4.0;}
                    else if (iddk1 == 3) {intn[12] = -2.0;}
                  }
                  block(4 + cnt1,30 + cnt2) = intn[1] + intn[12]*intn[6]/35.0;                                         //(pp|dd)
                }
                else if ((ippb1 != ippb2)&&(iddk1 != iddk2)) {
                  if (ippb2 == 3) {
                    if ((ippb1 == iddk1)&&(iddk2 == 5)) {block(4 + cnt1,30 + cnt2) = pow(-1.0,iddk1 - 1)*3.0*intn[6]/35.0;} //(pxpz|dxzdx2-y2),(pypz|dyzdx2-y2)
                    else if ((iddk2 == 4)&&(ippb1 + iddk1 == 3)) {block(4 + cnt1,30 + cnt2) = 3.0*intn[6]/35.0;}      //(pxpz|dyzdxy),(pypz|dxzdxy)
                    else if ((iddk2 == 3)&&(ippb1 == iddk1)) {block(4 + cnt1,30 + cnt2) = sqrt(3.0)*intn[6]/35.0;}    //(pxpz|dxzdz2),(pypz|dyzdz2)
                  }
                }
              }
            }
          }
        }
        cnt1 = 1;
        for (int idddb1 = 1; idddb1 < 6; ++idddb1) {
          block(10 + idddb1,10 + idddb1) = 0.2*intn[5];                                                               //(sdxz|sdxz),(sdyz|sdyz),(sdz2|sdz2),(sdxy|sdxy),(sdx2-y2|sdx2-y2)
          if (idddb1 < 3) {block(10 + idddb1,37 + idddb1) = pow(-1.0,idddb1 - 1)*sqrt(3.0/245.0)*intn[9];}            //(sdxz|dyzdxy),(sdyz|dyzdx2-y2)
          if (idddb1 != 3) {block(10 + idddb1,36 - idddb1) = sqrt(3.0/245.0)*intn[9];}                                //(sdxz|dxzdx2-y2),(sdyz|dxzdxy),(sdxy|dxzdyz),(sdx2-y2|dxzdxz)
          for (int idddb2 = idddb1; idddb2 < 6; ++idddb2, ++cnt1) {
            //bra and ket alike
            if (idddb1 == idddb2) {
              block(30 + cnt1,30 + cnt1) = intn[2] + 4.0*intn[7]/49.0 + 36.0*intn[11]/441.0;                          //(dxzdxz|dxzdxz),(dyzdyz|dyzdyz),(dz2dz2|dz2dz2),(dxydxy|dxydxy),(dx2-y2dx2-y2|dx2-y2dx2-y2)
              block(1,30 + cnt1) = intn[0];                                                                           //(ss|dxzdxz),(ss|dyzdyz),(ss|dz2dxz2),(ss|dxydxy),(ss|dx2-y2dx2-y2)
              intn[12] = 1.0;
              if (idddb1 == 3) {intn[12] = 2.0;}
              else if (idddb1 > 3) {intn[12] = -2.0;}
              block(13,30 + cnt1) = intn[12]*intn[9]/sqrt(245.0);                                                      //(sdz2|dd)
            }
            else {
              if ((idddb1 < 3)&&(idddb2 == 3)) {
                block(30 + cnt1,30 + cnt1) = intn[7]/49.0 + 30.0*intn[11]/441.0;                                      //(dxzdz2|dxzdz2),(dyzdz2|dyzdz2)
                block(10 + idddb1,30 + cnt1) = intn[9]/sqrt(245.0);                                                   //(sdxz|dz2dxz),(sdyz|dz2dyz)
              }
              else if (idddb1 == 3) {
                block(30 + cnt1,30 + cnt1) = 4.0*intn[7]/49.0 + 15.0*intn[11]/441.0;                                  //(dz2dxy|dz2dxy),(dz2dx2-y2|dz2dx2-y2)
                block(10 + idddb2,30 + cnt1) = -2.0*intn[9]/sqrt(245.0);                                              //(sdxy|dz2dxy),(sdx2-y2|dz2dx2-y2)
              }
              else if ((idddb1 < 3)&&(idddb2 != 3)) {block(30 + cnt1,30 + cnt1) = 3.0*intn[7]/49.0 + 20.0*intn[11]/441.0;}//(dxzdxy|dxzdxy),(dyzdxy|dyzdxy),(dxzdyz|dxzdyz),(dxzdx2-y2|dxzdx2-y2),(dyzdx2-y2|dyzdx2-y2)
            }
            cnt2 = 1;
            for (int idddk1 = 1; idddk1 < 6; ++idddk1) {
              for (int idddk2 = idddk1; idddk2 < 6; ++idddk2, ++cnt2) {
                if (cnt2 < cnt1) continue;
                if ((idddb1 == idddb2)&&(idddk1 == idddk2)&&(idddb1 == 3)&&((idddk1 > 3))) {block(30 + cnt1,30 + cnt2) = intn[2] - 4.0*intn[7]/49.0 + 6.0*intn[11]/441.0;}       //(dz2dz2|dxydxy),(dz2dz2|dx2-y2dx2-y2)
                else if ((idddb1 == idddb2)&&(idddk1 == idddk2)&&(idddk1 == 3)&&((idddb1 < 3))) {block(30 + cnt1,30 + cnt2) = intn[2] + 2.0*intn[7]/49.0 - 24.0*intn[11]/441.0;} //(dxzdxz|dz2dz2),(dyzdyz|dz2dz2)
                else if ((idddb1 == idddb2)&&(idddk1 == idddk2)&&(idddb1 != idddk1)&&(idddb1 != 3)&&((idddk1 != 3))) {
                  if (((idddb1 == 4)&&(idddk1 == 5))||((idddb1 == 5)&&(idddk1 == 4))) {block(30 + cnt1,30 + cnt2) = intn[2] + 4.0*intn[7]/49.0 - 34.0*intn[11]/441.0;}           //(dxydxy|dx2-y2dx2-y2)
                  else {block(30 + cnt1,30 + cnt2) = intn[2] - 2.0*intn[7]/49.0 - 4.0*intn[11]/441.0;}                                                                           //(dxzdxz|dyzdyz),(dxzdxz|dxydxy),(dxzdxz|dx2-y2dx2-y2),(dyzdyz|dxydxy),(dyzdyz|dx2-y2dx2-y2),(dyzdyz|dxzdxz),(dxydxy|dxzdxz),(dx2-y2dx2-y2|dxzdxz),(dxydxy|dyzdyz),(dx2-y2dx2-y2|dyzdyz)
                }
                else if ((idddk1 == 3)&&(idddk2 == 5)&&(idddb1 < 3)&&(idddb1 == idddb2)) {block(30 + cnt1,30 + cnt2) = pow(-1,idddb1+1)*(-sqrt(12.0)*intn[7]/49.0 + sqrt(300.0)*intn[11]/441.0);} //(dxzdxz|dz2dx2-y2),(dyzdyz|dz2dx2-y2)
                else if (((idddb1 < 3)&&(idddb2 == 3)&&(idddk1 < 3)&&(idddk2 > 3))||((idddb1 < 3)&&(idddb2 > 3)&&(idddk1 < 3)&&(idddk2 == 3))) {
                  if ((idddb1 + idddb2 + idddk1 + idddk2)%2 != 0) {continue;}
                  block(30 + cnt1,30 + cnt2) = pow(-1.0,DeltaDirac(idddk1+idddk2,7))*((sqrt(3.0)/49.0)*intn[7] - (sqrt(75.0)/441.0)*intn[11]);                                   //(dxzdz2|dxzdx2-y2),(dxzdz2|dyzdxy),(dyzdz2|dyzdx2-y2),(dxzdxy|dyzdz2)
                }
                else if (((idddb1 < 3)&&(idddb2 == 5)&&(idddk1 < 3)&&(idddk2 == 4))||((idddb1 < 3)&&(idddb2 == 4)&&(idddk1 < 3)&&(idddk2 == 5))) {
                  if (idddb1 == idddk1) {continue;}
                  block(30 + cnt1,30 + cnt2) = pow(-1.0,idddb1 + idddb2)*((3.0/49.0)*intn[7] - (15.0/441.0)*intn[11]);                                                           //(dxzdx2-y2|dyzdxy),(dxzdxy|dyzdx2-y2)
                }
              }
            }
          }
        }
        for (size_t idx = 0; idx < 12; ++idx) {
          if ((idx == 1)||(idx == 2)||(idx == 5)||(idx == 7)) {continue;}
          block(16 + idx,16 + idx) = 0.2*intn[3] + 24.0*intn[10]/245.0;                                               //(pxdxy|pxdxy),(pxdx2-y2|pxdx2-y2),(pydyz|pydyz),(pydxy|pydxy),(pydx2-y2|pydx2-y2),(pzdxz|pzdxz),(pzdyz|pzdyz)
        }
        cnt1 = 0;
        for (size_t idx = 4; idx < 13; cnt1 += idx, idx += 4) {
          cnt2 = cnt1;
          for (size_t idy = idx; idy < 13; cnt2 += idy, idy += 4) {
            block(17 + cnt1,17 + cnt2) = 3.0*intn[10]/49.0;                                                           //(pxdyz|pzdxy),(pydxz|pzdxy),(pxdyz|pydxz),(pxdyz|pxdyz),(pydxz|pydxz),(pzdxy|pzdxy)
          }
        }
        block(6,32) = 3.0*intn[6]/35.0;                                                                               //(pxpy|dxzdyz)
        block(5,42) = -sqrt(12.0)*intn[6]/35.0;                                                                       //(pxpx|dz2dx2-y2)
        block(6,41) = -sqrt(12.0)*intn[6]/35.0;                                                                       //(pxpy|dz2dxy)
        block(8,42) =  sqrt(12.0)*intn[6]/35.0;                                                                       //(pypy|dz2dx2-y2)
        block(15,36) = -sqrt(3.0/245.0)*intn[9];                                                                      //(sdx2-y2|dyzdyz)
        block(6,14) =  sqrt(3.0/125.0)*intn[8];                                                                       //(pxpy|sdxy)
        block(2,24) =  sqrt(3.0/45.0)*intn[4];                                                                        //(spx|pydxy)
        block(3,19) =  sqrt(3.0/45.0)*intn[4];                                                                        //(spy|pxdxy)
        block(4,28) =  2.0*intn[4]/sqrt(45.0);                                                                        //(spz|pzdz2)
        block(10,13) = 2.0*intn[8]/sqrt(125.0);                                                                       //(pzpz|sdz2)
        block(44,44) = 35.0*intn[11]/441.0;                                                                           //(dxydx2-y2|dxydx2-y2)
        block(32,41) = -sqrt(12.0)*intn[7]/49.0 + sqrt(300.0)*intn[11]/441.0;                                         //(dxzdyz|dz2dxy)
        block(18,18) = intn[3]/15.0 + 18.0*intn[10]/245.0;                                                            //(pxdz2|pxdz2)
        block(23,23) = intn[3]/15.0 + 18.0*intn[10]/245.0;                                                            //(pydz2|pydz2)
        block(28,28) = 4.0*intn[3]/15.0 + 27.0*intn[10]/245.0;                                                        //(pzdz2|pzdz2)
        block(18,20) = -sqrt(3.0)*intn[3]/15.0 - sqrt(27.0)*intn[10]/245.0;                                           //(pxdz2|pxdx2-y2)
        block(18,24) = -sqrt(3.0)*intn[3]/15.0 - sqrt(27.0)*intn[10]/245.0;                                           //(pxdz2|pydxy)
        block(19,23) = -sqrt(3.0)*intn[3]/15.0 - sqrt(27.0)*intn[10]/245.0;                                           //(pxdxy|pydz2)
        block(23,25) =  sqrt(3.0)*intn[3]/15.0 + sqrt(27.0)*intn[10]/245.0;                                           //(pydz2|pydx2-y2)
        block(16,30) =  3.0*intn[10]/49.0;                                                                            //(pxdxz|pzdx2-y2)
        block(22,30) = -3.0*intn[10]/49.0;                                                                            //(pydyz|pzdx2-y2)
        block(30,30) =  3.0*intn[10]/49.0;                                                                            //(pzdx2-y2|pzdx2-y2)
        block(20,24) =  0.2*intn[3] - 21.0*intn[10]/245.0;                                                            //(pxdx2-y2|pydxy)
        block(19,25) = -0.2*intn[3] + 21.0*intn[10]/245.0;                                                            //(pxdxy|pydx2-y2)
        block(16,28) = sqrt(12.0)*intn[3]/15.0 - sqrt(243.0)*intn[10]/245.0;                                          //(pxdxz|pzdz2)
        block(22,28) = sqrt(12.0)*intn[3]/15.0 - sqrt(243.0)*intn[10]/245.0;                                          //(pydyz|pzdz2)
        block(16,22) =  0.2*intn[3] - 6.0*intn[10]/245.0;                                                             //(pxdxz|pydyz)
        block(20,26) =  0.2*intn[3] - 6.0*intn[10]/245.0;                                                             //(pxdx2-y2|pzdxz)
        block(24,26) =  0.2*intn[3] - 6.0*intn[10]/245.0;                                                             //(pydxy|pzdxz)
        block(19,27) =  0.2*intn[3] - 6.0*intn[10]/245.0;                                                             //(pxdxy|pzdyz)
        block(25,27) = -0.2*intn[3] + 6.0*intn[10]/245.0;                                                             //(pydx2-y2|pzdyz)
        block(23,27) = -sqrt(3.0)*intn[3]/15.0 + sqrt(432.0)*intn[10]/245.0;                                          //(pydz2|pzdyz)
        block(18,26) = -sqrt(3.0)*intn[3]/15.0 + sqrt(432.0)*intn[10]/245.0;                                          //(pxdz2|pzdxz)
      }
      //symmetrize
      for (size_t idr = 1; idr < nrows + 1; ++idr) {
        for (size_t idc = idr + 1; idc < nrows + 1; ++idc) {
          block(idc,idr) = block(idr,idc);
        }
      }
    }
  }
  void IntegralBlock2C(int type, matrixE & block, int atmC, int atmD, double RCD, double cost, double sint, double cosp, double sinp) {
    //function calculating a block of two-center integrals
    int cnt1;
    int cnt2;
    size_t nrows = block.rows();
    size_t ncols = block.cols();
    if ((nrows > 1)||(ncols > 1)) {
      SProt = SPtransf(cost,sint,cosp,sinp);
      PProt = PPtransf(cost,sint,cosp,sinp);
      if ((nrows > 10)||(ncols > 10)) {
        SDrot = SDtransf(cost,sint,cosp,sinp);
        PDrot = PDtransf(cost,sint,cosp,sinp);
        DDrot = DDtransf(cost,sint,cosp,sinp);
      }
    }
    D[0] = Dvalue(atmC,1);
    D[1] = Dvalue(atmC,2);
    D[2] = Dvalue(atmD,1);
    D[3] = Dvalue(atmD,2);
    //integral calculation
    block(1,1) = eri2Center(0,0,0,0,0,0,0,0,RCD,D,atmC,atmD,type);       //(ss|ss)
    if (nrows > 1) {
      intn[0] = eri2Center(0,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center(1,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center(1,1,1,1,0,0,0,0,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(idx + 2,1) = -intn[0]*SProt(idx + 1,1);}                                   //(sp|ss) integrals
        block(5 + idx,1) = intn[1]*PProt(idx + 1,1) + intn[2]*(PProt(idx + 1,2) + PProt(idx + 1,3));   //(pp|ss) integrals
      }
      if (nrows > 10) {
        //(sd|ss) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(0,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(10 + idsd,1) = intn[0]*SDrot(idsd,1);
        }
        //(pd|ss) integrals
        Dd[0] = Dvalue(atmC,5);
        intn[0] = d_eri2Center(1,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center(1,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(5*(idpd2+2) + idpd1,1) = -(intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3)));
          }
        }
        //(dd|ss) integrals
        Dd[0] = Dvalue(atmC,6);
        intn[0] = d_eri2Center(2,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(30 + iddd,1) = intn[0]*DDrot(iddd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3)) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5));
        }
      }
    }
    if ((ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center(0,0,0,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center(0,0,0,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center(0,0,0,0,1,1,1,1,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(1,idx + 2) = -intn[0]*SProt(idx + 1,1);}                                   //(ss|sp) integrals
        block(1,5 + idx) = intn[1]*PProt(idx + 1,1) + intn[2]*(PProt(idx + 1,2) + PProt(idx + 1,3));   //(ss|pp) integrals
      }
      if (ncols > 10) {
        //(ss|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(0,0,0,0,0,0,2,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(1,10 + idsd) = intn[0]*SDrot(idsd,1);
        }
        //(ss|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(0,0,0,0,1,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center(0,0,0,0,1,1,2,1,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(1,5*(idpd2 + 2) + idpd1) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3));
          }
        }
        //(ss|dd) integrals
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(0,0,0,0,2,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center(0,0,0,0,2,1,2,1,RCD,Dd,atmC,atmD,type);
        intn[2] = d_eri2Center(0,0,0,0,2,2,2,2,RCD,Dd,atmC,atmD,type);
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(1,30 + iddd) = intn[0]*DDrot(iddd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3)) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5));
        }
      }
    }
    if ((nrows > 1)&&(ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center(0,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center(0,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center(0,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[3] = eri2Center(0,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);
      intn[4] = eri2Center(0,0,1,1,1,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idpb = 1; idpb < 4; ++idpb) {
        for (size_t idpk = 1; idpk < 4; ++idpk) {           //(sp|sp)
          block(idpb + 1,idpk + 1) = intn[0]*SProt(idpb,1)*SProt(idpk,1) + intn[1]*(SProt(idpb,2)*SProt(idpk,2) + SProt(idpb,3)*SProt(idpk,3));
        }
        for (size_t idc = 1; idc < 7; ++idc) {              //(sp|pp)
          block(idpb + 1,idc + 4) = -SProt(idpb,1)*(intn[2]*PProt(idc,1) + intn[3]*(PProt(idc,2) + PProt(idc,3))) - intn[4]*(SProt(idpb,2)*PProt(idc,4) + SProt(idpb,3)*PProt(idc,5));
        }
      }
      intn[0] = eri2Center(1,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);                          //(ps,ps|ps,ps)
      intn[1] = eri2Center(1,1,1,1,1,0,1,0,RCD,D,atmC,atmD,type);                          //(pp,pp|ps,ps)
      intn[2] = eri2Center(1,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);                          //(ps,ps|pp,pp)
      intn[3] = eri2Center(1,1,1,1,1,1,1,1,RCD,D,atmC,atmD,type);                          //(pp,pp|pp,pp)
      intn[4] = eri2Center(1,1,1,0,1,1,1,0,RCD,D,atmC,atmD,type);                          //(pp,ps|pp,ps)
      intn[5] = eri2Center(1,1,1,1,1,-1,1,-1,RCD,D,atmC,atmD,type);                        //(pp,pp|pp*,pp*)
      intn[6] = eri2Center(1,-1,1,1,1,-1,1,1,RCD,D,atmC,atmD,type);                        //(pp*,pp|pp*,pp)
      intn[7] = eri2Center(1,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[8] = eri2Center(1,1,1,1,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[9] = eri2Center(1,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idr = 1; idr < 7; ++idr) {
        for (size_t idc = 1; idc < 4; ++idc) {                //(pp|sp)
          block(idr + 4,idc + 1) = -SProt(idc,1)*(intn[7]*PProt(idr,1) + intn[8]*(PProt(idr,2) + PProt(idr,3))) - intn[9]*(SProt(idc,2)*PProt(idr,4) + SProt(idc,3)*PProt(idr,5));
        }
        for (size_t idc = 1; idc < 7; ++idc) {                //(pp|pp)
          block(idr + 4,idc + 4) = intn[0]*PProt(idr,1)*PProt(idc,1) + intn[1]*(PProt(idr,2) + PProt(idr,3))*PProt(idc,1) + intn[2]*PProt(idr,1)*(PProt(idc,2) + PProt(idc,3)) + intn[3]*(PProt(idr,2)*PProt(idc,2) + PProt(idr,3)*PProt(idc,3)) + intn[4]*(PProt(idr,4)*PProt(idc,4) + PProt(idr,5)*PProt(idc,5)) + intn[5]*(PProt(idr,2)*PProt(idc,3) + PProt(idr,3)*PProt(idc,2)) + intn[6]*PProt(idr,6)*PProt(idc,6);
        }
      }
      if (nrows > 10) {
        //(sd|sp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(0,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sd|sp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idsp = 1; idsp < 4; ++idsp) {
            block(10 + idsd,1 + idsp) = -intn[0]*SDrot(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot(idsd,2)*SProt(idsp,2) + SDrot(idsd,3)*SProt(idsp,3));
          }
        }
        //(sd|pp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center(0,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(0,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(sd|pp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(10 + idsd,4 + idpp) = intn[0]*SDrot(idsd,1)*PProt(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot(idsd,2)*PProt(idpp,4) + SDrot(idsd,3)*PProt(idpp,5));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot(idsd,4)*PProt(idpp,2) - SDrot(idsd,4)*PProt(idpp,3) + SDrot(idsd,5)*PProt(idpp,6));
          }
        }
        //(pd|sp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(1,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(1,1,2,0,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        intn[4] = d_eri2Center(1,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,1 + idsp) = intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
            }
          }
        }
        //(dd|sp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(2,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center(2,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(2,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(30 + iddd,1 + idsp) = -intn[0]*DDrot(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt(idsp,1);
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot(iddd,6)*SProt(idsp,2) + DDrot(iddd,7)*SProt(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt(idsp,3));
          }
        }
        //(pd|pp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center(1,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center(1,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center(1,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center(1,1,2,0,1,1,1,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center(1,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|pp) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,4 + idpp) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
            }
          }
        }
        //(dd|pp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center(2,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center(2,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[4] = d_eri2Center(2,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center(2,2,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center(2,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(2,-1,2,-1,1,1,1,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center(2,1,2,-1,1,1,1,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center(2,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=(sg,-dl|pi,-pi)=-(sg,dl|-pi,-pi)
        intn[10] = d_eri2Center(2,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|pp) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(30 + iddd,4 + idpp) = intn[0]*DDrot(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt(idpp,1);
            block(30 + iddd,4 + idpp) += intn[3]*DDrot(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt(idpp,2) + DDrot(iddd,3)*PProt(idpp,3));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt(idpp,4) + DDrot(iddd,7)*PProt(idpp,5));
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot(iddd,2)*PProt(idpp,3) + DDrot(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt(idpp,6);
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot(iddd,9)*PProt(idpp,6));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt(idpp,5));
          }
        }
      }
      if (ncols > 10) {
        //(sp|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(0,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sp|sd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(1 + idsp,10 + idsd) = -intn[0]*SDrot(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot(idsd,2)*SProt(idsp,2) + SDrot(idsd,3)*SProt(idsp,3));
          }
        }
        //(pp|sd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(1,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,1,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(1,1,1,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        //(pp|sd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,10 + idsd) = intn[0]*SDrot(idsd,1)*PProt(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot(idsd,2)*PProt(idpp,4) + SDrot(idsd,3)*PProt(idpp,5));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot(idsd,4)*PProt(idpp,2) - SDrot(idsd,4)*PProt(idpp,3) + SDrot(idsd,5)*PProt(idpp,6));
          }
        }
        //(sp|pd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(0,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(0,0,1,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[4] = d_eri2Center(0,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|pd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(1 + idsp,5*(idpd2+2) + idpd1) = intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
            }
          }
        }
        //(sp|dd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(0,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(0,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(0,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|dd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(1 + idsp,30 + iddd) = -intn[0]*DDrot(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot(iddd,6)*SProt(idsp,2) + DDrot(iddd,7)*SProt(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt(idsp,3));
          }
        }
        //(pp|pd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(1,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,1,1,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center(1,1,1,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center(1,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center(1,1,1,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center(1,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|pd) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(4 + idpp,5*(idpd2+2) + idpd1) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
            }
          }
        }
        //(pp|dd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(1,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(1,1,1,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center(1,1,1,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center(1,1,1,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)=(-pi,-pi|-dl,-dl)
        intn[6] = d_eri2Center(1,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(1,-1,1,-1,2,1,2,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center(1,1,1,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center(1,1,1,1,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=(pi,-pi|sg,-dl)=-(-pi,-pi|sg,dl)
        intn[10] = d_eri2Center(1,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,30 + iddd) = intn[0]*DDrot(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt(idpp,1);
            block(4 + idpp,30 + iddd) += intn[3]*DDrot(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt(idpp,2) + DDrot(iddd,3)*PProt(idpp,3));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt(idpp,4) + DDrot(iddd,7)*PProt(idpp,5));
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot(iddd,2)*PProt(idpp,3) + DDrot(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt(idpp,6);
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot(iddd,9)*PProt(idpp,6));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt(idpp,5));
          }
        }
      }
      if ((nrows > 10)&&(ncols > 10)) {
        //(sd|sd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(0,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[2] = d_eri2Center(0,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        //(sd|sd) rotate
        for (size_t idsdb = 1; idsdb < 6; ++idsdb) {
          for (size_t idsdk = 1; idsdk < 6; ++idsdk) {
            block(10 + idsdb,10 + idsdk) = intn[0]*SDrot(idsdb,1)*SDrot(idsdk,1) + intn[1]*(SDrot(idsdb,2)*SDrot(idsdk,2) + SDrot(idsdb,3)*SDrot(idsdk,3)) + intn[2]*(SDrot(idsdb,4)*SDrot(idsdk,4) + SDrot(idsdb,5)*SDrot(idsdk,5));
          }
        }
        //(sd|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(0,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(0,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(sg,pi|-pi,-dl)
        intn[4] = d_eri2Center(0,0,2,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        //(sd|pd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(10 + idsd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
            }
          }
        }
        //(pd|sd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(1,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(1,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[4] = d_eri2Center(1,1,2,0,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        //(pd|sd) rotate
        for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
            for (size_t idsd = 1; idsd < 6; ++idsd) {
              block(5*(idpd2+2) + idpd1,10 + idsd) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
            }
          }
        }
        //(dd|sd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(2,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center(2,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(2,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center(2,1,2,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        intn[6] = d_eri2Center(2,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)
        //(dd|sd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(30 + iddd,10 + idsd) = intn[0]*DDrot(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot(idsd,1);
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot(iddd,6)*SDrot(idsd,2) + DDrot(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot(idsd,4) + DDrot(iddd,9)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot(idsd,4) + DDrot(iddd,10)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot(idsd,3));
          }
        }
        //(sd|dd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(0,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(0,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(0,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center(0,0,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        intn[6] = d_eri2Center(0,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        //(sd|dd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(10 + idsd,30 + iddd) = intn[0]*DDrot(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot(idsd,1);
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot(iddd,6)*SDrot(idsd,2) + DDrot(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot(idsd,4) + DDrot(iddd,9)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot(idsd,4) + DDrot(iddd,10)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot(idsd,3));
          }
        }
        //(pd|pd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(1,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)=(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|-pi,-pi)
        intn[3] = d_eri2Center(1,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(1,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[5] = d_eri2Center(1,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(pi,-dl|pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)
        intn[6] = d_eri2Center(1,1,2,0,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)=(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[7] = d_eri2Center(1,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[8] = d_eri2Center(1,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)=(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        //(pd|pd) rotate
        for (size_t idpdb1 = 1; idpdb1 < 6; ++idpdb1) {            //d orbital on bra
          for (size_t idpdb2 = 1; idpdb2 < 4; ++idpdb2) {          //p orbital on bra
            for (size_t idpdk1 = 1; idpdk1 < 6; ++idpdk1) {        //d orbital on ket
              for (size_t idpdk2 = 1; idpdk2 < 4; ++idpdk2) {      //p orbital on ket
                cnt1 = (idpdb1 - 1)*3 + idpdb2;
                cnt2 = (idpdk1 - 1)*3 + idpdk2;
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) = intn[0]*PDrot(cnt1,1)*PDrot(cnt2,1) + intn[1]*((PDrot(cnt1,2) + PDrot(cnt1,3))*PDrot(cnt2,1) + PDrot(cnt1,1)*(PDrot(cnt2,2) + PDrot(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(PDrot(cnt2,2) + PDrot(cnt2,3)) + intn[3]*(PDrot(cnt1,4)*PDrot(cnt2,4) + PDrot(cnt1,5)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot(cnt1,8)*PDrot(cnt2,8) + PDrot(cnt1,12)*PDrot(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(PDrot(cnt2,14) - PDrot(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot(cnt1,4)*PDrot(cnt2,8) + PDrot(cnt1,8)*PDrot(cnt2,4) + PDrot(cnt1,5)*PDrot(cnt2,12) + PDrot(cnt1,12)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot(cnt2,4) + PDrot(cnt1,4)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot(cnt2,5) + PDrot(cnt1,5)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot(cnt2,8) + PDrot(cnt1,8)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot(cnt2,12) + PDrot(cnt1,12)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
              }
            }
          }
        }
        //(pd|dd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(1,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center(1,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center(1,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[5] = d_eri2Center(1,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|-dl,-dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)
        intn[6] = d_eri2Center(1,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(1,1,2,0,2,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center(1,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=(pi,-dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center(1,1,2,0,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)
        intn[10] = d_eri2Center(1,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        intn[11] = d_eri2Center(1,1,2,2,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(5*(idpd2+2) + idpd1,30 + iddd) = -intn[0]*PDrot(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot(cnt1,4)*DDrot(iddd,6) + PDrot(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot(iddd,6) + PDrot(cnt1,12)*DDrot(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot(iddd,7));
            }
          }
        }
        //(dd|pd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(2,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center(2,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center(2,2,2,2,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center(2,2,2,2,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|-pi,-pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)
        intn[6] = d_eri2Center(2,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(2,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center(2,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center(2,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        intn[10] = d_eri2Center(2,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[11] = d_eri2Center(2,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(dd|pd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(30 + iddd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot(cnt1,4)*DDrot(iddd,6) + PDrot(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot(iddd,6) + PDrot(cnt1,12)*DDrot(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot(iddd,7));
            }
          }
        }
        //(dd|dd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(2,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(2,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(2,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center(2,2,2,2,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center(2,2,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center(2,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|-dl,-dl)
        intn[7] = d_eri2Center(2,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[8] = d_eri2Center(2,2,2,2,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(dl,dl|dl,dl)=(-dl,-dl|-dl,-dl)
        intn[9] = d_eri2Center(2,1,2,1,2,-1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[10] = d_eri2Center(2,-2,2,-2,2,2,2,2,RCD,Dd,atmC,atmD,type); //(-dl,-dl|dl,dl)=(dl,dl|-dl,-dl)
        intn[11] = d_eri2Center(2,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[12] = d_eri2Center(2,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[13] = d_eri2Center(2,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)
        intn[14] = d_eri2Center(2,1,2,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type); //(pi,-pi|pi,-pi)
        intn[15] = d_eri2Center(2,1,2,-2,2,0,2,-1,RCD,Dd,atmC,atmD,type); //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)=(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        intn[16] = d_eri2Center(2,1,2,1,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)=(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(dd|dd) rotate
        for (size_t iddd1 = 1; iddd1 < 16; ++iddd1) {
          for (size_t iddd2 = 1; iddd2 < 16; ++iddd2) {
            block(30 + iddd1,30 + iddd2) = intn[0]*DDrot(iddd1,1)*DDrot(iddd2,1) + intn[1]*DDrot(iddd1,1)*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[2]*DDrot(iddd1,1)*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot(iddd1,2) + DDrot(iddd1,3))*DDrot(iddd2,1) + intn[4]*(DDrot(iddd1,4) + DDrot(iddd1,5))*DDrot(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot(iddd1,4) + DDrot(iddd1,5))*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[6]*(DDrot(iddd1,2) + DDrot(iddd1,3))*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot(iddd1,2)*DDrot(iddd2,2) + DDrot(iddd1,3)*DDrot(iddd2,3)) + intn[8]*(DDrot(iddd1,4)*DDrot(iddd2,4) + DDrot(iddd1,5)*DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot(iddd1,2)*DDrot(iddd2,3) + DDrot(iddd1,3)*DDrot(iddd2,2)) + intn[10]*(DDrot(iddd1,4)*DDrot(iddd2,5) + DDrot(iddd1,5)*DDrot(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot(iddd1,6)*DDrot(iddd2,6) + DDrot(iddd1,7)*DDrot(iddd2,7)) + intn[12]*(DDrot(iddd1,8)*DDrot(iddd2,8) + DDrot(iddd1,9)*DDrot(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot(iddd1,11) + DDrot(iddd1,14))*(DDrot(iddd2,11) + DDrot(iddd2,14)) + (DDrot(iddd1,12) - DDrot(iddd1,13))*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[14]*DDrot(iddd1,10)*DDrot(iddd2,10);
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot(iddd1,11) + DDrot(iddd1,14))*DDrot(iddd2,6) + (DDrot(iddd1,12) - DDrot(iddd1,13))*DDrot(iddd2,7) + DDrot(iddd1,6)*(DDrot(iddd2,11) + DDrot(iddd2,14)) + DDrot(iddd1,7)*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot(iddd1,2) - DDrot(iddd1,3))*DDrot(iddd2,8) + DDrot(iddd1,10)*DDrot(iddd2,9) + DDrot(iddd1,8)*(DDrot(iddd2,2) - DDrot(iddd2,3)) + DDrot(iddd1,9)*DDrot(iddd2,10));
          }
        }
      }
    }
  }
  void IntegralBlock2C_dR(int type, matrixE & block, int atmC, int atmD, double RCD, double cost, double sint, double cosp, double sinp) {
    //function calculating a block of first-derivatives of two-center integrals with respect to the internuclear distance
    int cnt1;
    int cnt2;
    size_t nrows = block.rows();
    size_t ncols = block.cols();
    if ((nrows > 1)||(ncols > 1)) {
      SProt = SPtransf(cost,sint,cosp,sinp);
      PProt = PPtransf(cost,sint,cosp,sinp);
      if ((nrows > 10)||(ncols > 10)) {
        SDrot = SDtransf(cost,sint,cosp,sinp);
        PDrot = PDtransf(cost,sint,cosp,sinp);
        DDrot = DDtransf(cost,sint,cosp,sinp);
      }
    }
    D[0] = Dvalue(atmC,1);
    D[1] = Dvalue(atmC,2);
    D[2] = Dvalue(atmD,1);
    D[3] = Dvalue(atmD,2);
    //integral calculation
    block(1,1) = eri2Center_dR(0,0,0,0,0,0,0,0,RCD,D,atmC,atmD,type);       //(ss|ss)
    if (nrows > 1) {
      intn[0] = eri2Center_dR(0,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center_dR(1,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center_dR(1,1,1,1,0,0,0,0,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(idx + 2,1) = -intn[0]*SProt(idx + 1,1);}                                                   //(sp|ss) integrals
        block(5 + idx,1) = intn[1]*PProt(idx + 1,1) + intn[2]*(PProt(idx + 1,2) + PProt(idx + 1,3));   //(pp|ss) integrals
      }
      if (nrows > 10) {
        //(sd|ss) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR(0,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(10 + idsd,1) = intn[0]*SDrot(idsd,1);
        }
        //(pd|ss) integrals
        Dd[0] = Dvalue(atmC,5);
        intn[0] = d_eri2Center_dR(1,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center_dR(1,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(5*(idpd2+2) + idpd1,1) = -(intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3)));
          }
        }
        //(dd|ss) integrals
        Dd[0] = Dvalue(atmC,6);
        intn[0] = d_eri2Center_dR(2,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(2,2,2,2,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(30 + iddd,1) = intn[0]*DDrot(iddd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3)) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5));
        }
      }
    }
    if ((ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center_dR(0,0,0,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center_dR(0,0,0,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center_dR(0,0,0,0,1,1,1,1,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(1,idx + 2) = -intn[0]*SProt(idx + 1,1);}                                                   //(ss|sp) integrals
        block(1,5 + idx) = intn[1]*PProt(idx + 1,1) + intn[2]*(PProt(idx + 1,2) + PProt(idx + 1,3));   //(ss|pp) integrals
      }
      if (ncols > 10) {
        //(ss|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(0,0,0,0,0,0,2,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(1,10 + idsd) = intn[0]*SDrot(idsd,1);
        }
        //(ss|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(0,0,0,0,1,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center_dR(0,0,0,0,1,1,2,1,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(1,5*(idpd2 + 2) + idpd1) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3));
          }
        }
        //(ss|dd) integrals
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(0,0,0,0,2,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center_dR(0,0,0,0,2,1,2,1,RCD,Dd,atmC,atmD,type);
        intn[2] = d_eri2Center_dR(0,0,0,0,2,2,2,2,RCD,Dd,atmC,atmD,type);
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(1,30 + iddd) = intn[0]*DDrot(iddd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3)) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5));
        }
      }
    }
    if ((nrows > 1)&&(ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center_dR(0,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center_dR(0,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center_dR(0,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[3] = eri2Center_dR(0,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);
      intn[4] = eri2Center_dR(0,0,1,1,1,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idpb = 1; idpb < 4; ++idpb) {
        for (size_t idpk = 1; idpk < 4; ++idpk) {           //(sp|sp)
          block(idpb + 1,idpk + 1) = intn[0]*SProt(idpb,1)*SProt(idpk,1) + intn[1]*(SProt(idpb,2)*SProt(idpk,2) + SProt(idpb,3)*SProt(idpk,3));
        }
        for (size_t idc = 1; idc < 7; ++idc) {              //(sp|pp)
          block(idpb + 1,idc + 4) = -SProt(idpb,1)*(intn[2]*PProt(idc,1) + intn[3]*(PProt(idc,2) + PProt(idc,3))) - intn[4]*(SProt(idpb,2)*PProt(idc,4) + SProt(idpb,3)*PProt(idc,5));
        }
      }
      intn[0] = eri2Center_dR(1,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);                          //(ps,ps|ps,ps)
      intn[1] = eri2Center_dR(1,1,1,1,1,0,1,0,RCD,D,atmC,atmD,type);                          //(pp,pp|ps,ps)
      intn[2] = eri2Center_dR(1,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);                          //(ps,ps|pp,pp)
      intn[3] = eri2Center_dR(1,1,1,1,1,1,1,1,RCD,D,atmC,atmD,type);                          //(pp,pp|pp,pp)
      intn[4] = eri2Center_dR(1,1,1,0,1,1,1,0,RCD,D,atmC,atmD,type);                          //(pp,ps|pp,ps)
      intn[5] = eri2Center_dR(1,1,1,1,1,-1,1,-1,RCD,D,atmC,atmD,type);                        //(pp,pp|pp*,pp*)
      intn[6] = eri2Center_dR(1,-1,1,1,1,-1,1,1,RCD,D,atmC,atmD,type);                        //(pp*,pp|pp*,pp)
      intn[7] = eri2Center_dR(1,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[8] = eri2Center_dR(1,1,1,1,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[9] = eri2Center_dR(1,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idr = 1; idr < 7; ++idr) {
        for (size_t idc = 1; idc < 4; ++idc) {                //(pp|sp)
          block(idr + 4,idc + 1) = -SProt(idc,1)*(intn[7]*PProt(idr,1) + intn[8]*(PProt(idr,2) + PProt(idr,3))) - intn[9]*(SProt(idc,2)*PProt(idr,4) + SProt(idc,3)*PProt(idr,5));
        }
        for (size_t idc = 1; idc < 7; ++idc) {                //(pp|pp)
          block(idr + 4,idc + 4) = intn[0]*PProt(idr,1)*PProt(idc,1) + intn[1]*(PProt(idr,2) + PProt(idr,3))*PProt(idc,1) + intn[2]*PProt(idr,1)*(PProt(idc,2) + PProt(idc,3)) + intn[3]*(PProt(idr,2)*PProt(idc,2) + PProt(idr,3)*PProt(idc,3)) + intn[4]*(PProt(idr,4)*PProt(idc,4) + PProt(idr,5)*PProt(idc,5)) + intn[5]*(PProt(idr,2)*PProt(idc,3) + PProt(idr,3)*PProt(idc,2)) + intn[6]*PProt(idr,6)*PProt(idc,6);
        }
      }
      if (nrows > 10) {
        //(sd|sp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR(0,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sd|sp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idsp = 1; idsp < 4; ++idsp) {
            block(10 + idsd,1 + idsp) = -intn[0]*SDrot(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot(idsd,2)*SProt(idsp,2) + SDrot(idsd,3)*SProt(idsp,3));
          }
        }
        //(sd|pp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center_dR(0,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(0,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(0,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(sd|pp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(10 + idsd,4 + idpp) = intn[0]*SDrot(idsd,1)*PProt(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot(idsd,2)*PProt(idpp,4) + SDrot(idsd,3)*PProt(idpp,5));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot(idsd,4)*PProt(idpp,2) - SDrot(idsd,4)*PProt(idpp,3) + SDrot(idsd,5)*PProt(idpp,6));
          }
        }
        //(pd|sp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR(1,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(1,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(1,1,2,0,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        intn[4] = d_eri2Center_dR(1,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,1 + idsp) = intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
            }
          }
        }
        //(dd|sp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR(2,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(2,2,2,2,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center_dR(2,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR(2,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(30 + iddd,1 + idsp) = -intn[0]*DDrot(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt(idsp,1);
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot(iddd,6)*SProt(idsp,2) + DDrot(iddd,7)*SProt(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt(idsp,3));
          }
        }
        //(pd|pp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center_dR(1,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(1,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center_dR(1,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center_dR(1,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center_dR(1,1,2,0,1,1,1,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center_dR(1,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|pp) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,4 + idpp) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
            }
          }
        }
        //(dd|pp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center_dR(2,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(2,2,2,2,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center_dR(2,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[4] = d_eri2Center_dR(2,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center_dR(2,2,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center_dR(2,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR(2,-1,2,-1,1,1,1,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center_dR(2,1,2,-1,1,1,1,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center_dR(2,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=(sg,-dl|pi,-pi)=-(sg,dl|-pi,-pi)
        intn[10] = d_eri2Center_dR(2,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|pp) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(30 + iddd,4 + idpp) = intn[0]*DDrot(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt(idpp,1);
            block(30 + iddd,4 + idpp) += intn[3]*DDrot(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt(idpp,2) + DDrot(iddd,3)*PProt(idpp,3));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt(idpp,4) + DDrot(iddd,7)*PProt(idpp,5));
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot(iddd,2)*PProt(idpp,3) + DDrot(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt(idpp,6);
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot(iddd,9)*PProt(idpp,6));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt(idpp,5));
          }
        }
      }
      if (ncols > 10) {
        //(sp|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(0,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sp|sd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(1 + idsp,10 + idsd) = -intn[0]*SDrot(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot(idsd,2)*SProt(idsp,2) + SDrot(idsd,3)*SProt(idsp,3));
          }
        }
        //(pp|sd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(1,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,1,1,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(1,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(1,1,1,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        //(pp|sd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,10 + idsd) = intn[0]*SDrot(idsd,1)*PProt(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot(idsd,2)*PProt(idpp,4) + SDrot(idsd,3)*PProt(idpp,5));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot(idsd,4)*PProt(idpp,2) - SDrot(idsd,4)*PProt(idpp,3) + SDrot(idsd,5)*PProt(idpp,6));
          }
        }
        //(sp|pd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(0,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(0,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(0,0,1,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[4] = d_eri2Center_dR(0,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|pd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(1 + idsp,5*(idpd2+2) + idpd1) = intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
            }
          }
        }
        //(sp|dd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(0,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(0,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR(0,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR(0,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|dd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(1 + idsp,30 + iddd) = -intn[0]*DDrot(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot(iddd,6)*SProt(idsp,2) + DDrot(iddd,7)*SProt(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt(idsp,3));
          }
        }
        //(pp|pd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(1,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(1,1,1,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center_dR(1,1,1,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center_dR(1,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center_dR(1,1,1,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center_dR(1,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|pd) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(4 + idpp,5*(idpd2+2) + idpd1) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
            }
          }
        }
        //(pp|dd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(1,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(1,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR(1,1,1,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center_dR(1,1,1,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center_dR(1,1,1,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)=(-pi,-pi|-dl,-dl)
        intn[6] = d_eri2Center_dR(1,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR(1,-1,1,-1,2,1,2,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center_dR(1,1,1,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center_dR(1,1,1,1,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=(pi,-pi|sg,-dl)=-(-pi,-pi|sg,dl)
        intn[10] = d_eri2Center_dR(1,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,30 + iddd) = intn[0]*DDrot(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt(idpp,1);
            block(4 + idpp,30 + iddd) += intn[3]*DDrot(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt(idpp,2) + DDrot(iddd,3)*PProt(idpp,3));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt(idpp,4) + DDrot(iddd,7)*PProt(idpp,5));
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot(iddd,2)*PProt(idpp,3) + DDrot(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt(idpp,6);
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot(iddd,9)*PProt(idpp,6));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt(idpp,5));
          }
        }
      }
      if ((nrows > 10)&&(ncols > 10)) {
        //(sd|sd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(0,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[2] = d_eri2Center_dR(0,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        //(sd|sd) rotate
        for (size_t idsdb = 1; idsdb < 6; ++idsdb) {
          for (size_t idsdk = 1; idsdk < 6; ++idsdk) {
            block(10 + idsdb,10 + idsdk) = intn[0]*SDrot(idsdb,1)*SDrot(idsdk,1) + intn[1]*(SDrot(idsdb,2)*SDrot(idsdk,2) + SDrot(idsdb,3)*SDrot(idsdk,3)) + intn[2]*(SDrot(idsdb,4)*SDrot(idsdk,4) + SDrot(idsdb,5)*SDrot(idsdk,5));
          }
        }
        //(sd|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(0,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(0,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(0,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(sg,pi|-pi,-dl)
        intn[4] = d_eri2Center_dR(0,0,2,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        //(sd|pd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(10 + idsd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
            }
          }
        }
        //(pd|sd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(1,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(1,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(1,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[4] = d_eri2Center_dR(1,1,2,0,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        //(pd|sd) rotate
        for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
            for (size_t idsd = 1; idsd < 6; ++idsd) {
              block(5*(idpd2+2) + idpd1,10 + idsd) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
            }
          }
        }
        //(dd|sd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(2,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(2,2,2,2,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center_dR(2,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR(2,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center_dR(2,1,2,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        intn[6] = d_eri2Center_dR(2,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)
        //(dd|sd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(30 + iddd,10 + idsd) = intn[0]*DDrot(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot(idsd,1);
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot(iddd,6)*SDrot(idsd,2) + DDrot(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot(idsd,4) + DDrot(iddd,9)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot(idsd,4) + DDrot(iddd,10)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot(idsd,3));
          }
        }
        //(sd|dd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(0,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(0,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR(0,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR(0,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center_dR(0,0,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        intn[6] = d_eri2Center_dR(0,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        //(sd|dd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(10 + idsd,30 + iddd) = intn[0]*DDrot(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot(idsd,1);
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot(iddd,6)*SDrot(idsd,2) + DDrot(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot(idsd,4) + DDrot(iddd,9)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot(idsd,4) + DDrot(iddd,10)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot(idsd,3));
          }
        }
        //(pd|pd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(1,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)=(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(1,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|-pi,-pi)
        intn[3] = d_eri2Center_dR(1,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR(1,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[5] = d_eri2Center_dR(1,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(pi,-dl|pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)
        intn[6] = d_eri2Center_dR(1,1,2,0,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)=(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[7] = d_eri2Center_dR(1,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[8] = d_eri2Center_dR(1,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)=(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        //(pd|pd) rotate
        for (size_t idpdb1 = 1; idpdb1 < 6; ++idpdb1) {            //d orbital on bra
          for (size_t idpdb2 = 1; idpdb2 < 4; ++idpdb2) {          //p orbital on bra
            for (size_t idpdk1 = 1; idpdk1 < 6; ++idpdk1) {        //d orbital on ket
              for (size_t idpdk2 = 1; idpdk2 < 4; ++idpdk2) {      //p orbital on ket
                cnt1 = (idpdb1 - 1)*3 + idpdb2;
                cnt2 = (idpdk1 - 1)*3 + idpdk2;
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) = intn[0]*PDrot(cnt1,1)*PDrot(cnt2,1) + intn[1]*((PDrot(cnt1,2) + PDrot(cnt1,3))*PDrot(cnt2,1) + PDrot(cnt1,1)*(PDrot(cnt2,2) + PDrot(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(PDrot(cnt2,2) + PDrot(cnt2,3)) + intn[3]*(PDrot(cnt1,4)*PDrot(cnt2,4) + PDrot(cnt1,5)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot(cnt1,8)*PDrot(cnt2,8) + PDrot(cnt1,12)*PDrot(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(PDrot(cnt2,14) - PDrot(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot(cnt1,4)*PDrot(cnt2,8) + PDrot(cnt1,8)*PDrot(cnt2,4) + PDrot(cnt1,5)*PDrot(cnt2,12) + PDrot(cnt1,12)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot(cnt2,4) + PDrot(cnt1,4)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot(cnt2,5) + PDrot(cnt1,5)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot(cnt2,8) + PDrot(cnt1,8)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot(cnt2,12) + PDrot(cnt1,12)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
              }
            }
          }
        }
        //(pd|dd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(1,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(1,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center_dR(1,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center_dR(1,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[5] = d_eri2Center_dR(1,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|-dl,-dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)
        intn[6] = d_eri2Center_dR(1,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR(1,1,2,0,2,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center_dR(1,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=(pi,-dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center_dR(1,1,2,0,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)
        intn[10] = d_eri2Center_dR(1,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        intn[11] = d_eri2Center_dR(1,1,2,2,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(5*(idpd2+2) + idpd1,30 + iddd) = -intn[0]*PDrot(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot(cnt1,4)*DDrot(iddd,6) + PDrot(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot(iddd,6) + PDrot(cnt1,12)*DDrot(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot(iddd,7));
            }
          }
        }
        //(dd|pd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(2,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(2,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center_dR(2,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center_dR(2,2,2,2,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center_dR(2,2,2,2,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|-pi,-pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)
        intn[6] = d_eri2Center_dR(2,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR(2,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center_dR(2,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center_dR(2,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        intn[10] = d_eri2Center_dR(2,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[11] = d_eri2Center_dR(2,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(dd|pd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(30 + iddd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot(cnt1,4)*DDrot(iddd,6) + PDrot(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot(iddd,6) + PDrot(cnt1,12)*DDrot(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot(iddd,7));
            }
          }
        }
        //(dd|dd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(2,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(2,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR(2,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center_dR(2,2,2,2,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center_dR(2,2,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center_dR(2,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|-dl,-dl)
        intn[7] = d_eri2Center_dR(2,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[8] = d_eri2Center_dR(2,2,2,2,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(dl,dl|dl,dl)=(-dl,-dl|-dl,-dl)
        intn[9] = d_eri2Center_dR(2,1,2,1,2,-1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[10] = d_eri2Center_dR(2,-2,2,-2,2,2,2,2,RCD,Dd,atmC,atmD,type); //(-dl,-dl|dl,dl)=(dl,dl|-dl,-dl)
        intn[11] = d_eri2Center_dR(2,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[12] = d_eri2Center_dR(2,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[13] = d_eri2Center_dR(2,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)
        intn[14] = d_eri2Center_dR(2,1,2,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type); //(pi,-pi|pi,-pi)
        intn[15] = d_eri2Center_dR(2,1,2,-2,2,0,2,-1,RCD,Dd,atmC,atmD,type); //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)=(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        intn[16] = d_eri2Center_dR(2,1,2,1,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)=(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(dd|dd) rotate
        for (size_t iddd1 = 1; iddd1 < 16; ++iddd1) {
          for (size_t iddd2 = 1; iddd2 < 16; ++iddd2) {
            block(30 + iddd1,30 + iddd2) = intn[0]*DDrot(iddd1,1)*DDrot(iddd2,1) + intn[1]*DDrot(iddd1,1)*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[2]*DDrot(iddd1,1)*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot(iddd1,2) + DDrot(iddd1,3))*DDrot(iddd2,1) + intn[4]*(DDrot(iddd1,4) + DDrot(iddd1,5))*DDrot(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot(iddd1,4) + DDrot(iddd1,5))*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[6]*(DDrot(iddd1,2) + DDrot(iddd1,3))*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot(iddd1,2)*DDrot(iddd2,2) + DDrot(iddd1,3)*DDrot(iddd2,3)) + intn[8]*(DDrot(iddd1,4)*DDrot(iddd2,4) + DDrot(iddd1,5)*DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot(iddd1,2)*DDrot(iddd2,3) + DDrot(iddd1,3)*DDrot(iddd2,2)) + intn[10]*(DDrot(iddd1,4)*DDrot(iddd2,5) + DDrot(iddd1,5)*DDrot(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot(iddd1,6)*DDrot(iddd2,6) + DDrot(iddd1,7)*DDrot(iddd2,7)) + intn[12]*(DDrot(iddd1,8)*DDrot(iddd2,8) + DDrot(iddd1,9)*DDrot(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot(iddd1,11) + DDrot(iddd1,14))*(DDrot(iddd2,11) + DDrot(iddd2,14)) + (DDrot(iddd1,12) - DDrot(iddd1,13))*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[14]*DDrot(iddd1,10)*DDrot(iddd2,10);
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot(iddd1,11) + DDrot(iddd1,14))*DDrot(iddd2,6) + (DDrot(iddd1,12) - DDrot(iddd1,13))*DDrot(iddd2,7) + DDrot(iddd1,6)*(DDrot(iddd2,11) + DDrot(iddd2,14)) + DDrot(iddd1,7)*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot(iddd1,2) - DDrot(iddd1,3))*DDrot(iddd2,8) + DDrot(iddd1,10)*DDrot(iddd2,9) + DDrot(iddd1,8)*(DDrot(iddd2,2) - DDrot(iddd2,3)) + DDrot(iddd1,9)*DDrot(iddd2,10));
          }
        }
      }
    }
  }
  void IntegralBlock2C_dA(int type, matrixE & block, int der, int atmC, int atmD, double RCD, double cost, double sint, double cosp, double sinp) {
    //function calculating a block of first-derivatives of two-center integrals with respect to angles
    //2 -> theta; 3 -> phi
    matrixE SProt_dA(1,1);
    matrixE PProt_dA(1,1);
    matrixE SDrot_dA(1,1);
    matrixE PDrot_dA(1,1);
    matrixE DDrot_dA(1,1);
    int cnt1;
    int cnt2;
    size_t nrows = block.rows();
    size_t ncols = block.cols();
    if ((nrows > 1)||(ncols > 1)) {
      SProt = SPtransf(cost,sint,cosp,sinp);
      PProt = PPtransf(cost,sint,cosp,sinp);
      if (der == 2) {
        SProt_dA = SPtransf_dt(cost,sint,cosp,sinp);
        PProt_dA = PPtransf_dt(cost,sint,cosp,sinp);
      }
      else if (der == 3) {
        SProt_dA = SPtransf_dp(cost,sint,cosp,sinp);
        PProt_dA = PPtransf_dp(cost,sint,cosp,sinp);
      }
      if ((nrows > 10)||(ncols > 10)) {
        SDrot = SDtransf(cost,sint,cosp,sinp);
        PDrot = PDtransf(cost,sint,cosp,sinp);
        DDrot = DDtransf(cost,sint,cosp,sinp);
        if (der == 2) {
          SDrot_dA = SDtransf_dt(cost,sint,cosp,sinp);
          PDrot_dA = PDtransf_dt(cost,sint,cosp,sinp);
          DDrot_dA = DDtransf_dt(cost,sint,cosp,sinp);
        }
        else if (der == 3) {
          SDrot_dA = SDtransf_dp(cost,sint,cosp,sinp);
          PDrot_dA = PDtransf_dp(cost,sint,cosp,sinp);
          DDrot_dA = DDtransf_dp(cost,sint,cosp,sinp);
        }
      }
    }
    D[0] = Dvalue(atmC,1);
    D[1] = Dvalue(atmC,2);
    D[2] = Dvalue(atmD,1);
    D[3] = Dvalue(atmD,2);
    //integral calculation
    block(1,1) = 0.0;                                                                                                                   //(ss|ss)
    if (nrows > 1) {
      intn[0] = eri2Center(0,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center(1,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center(1,1,1,1,0,0,0,0,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(idx + 2,1) = -intn[0]*SProt_dA(idx + 1,1);}                                                         //(sp|ss) integrals
        block(5 + idx,1) = intn[1]*PProt_dA(idx + 1,1) + intn[2]*(PProt_dA(idx + 1,2) + PProt_dA(idx + 1,3));   //(pp|ss) integrals
      }
      if (nrows > 10) {
        //(sd|ss) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(0,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(10 + idsd,1) = intn[0]*SDrot_dA(idsd,1);
        }
        //(pd|ss) integrals
        Dd[0] = Dvalue(atmC,5);
        intn[0] = d_eri2Center(1,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center(1,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(5*(idpd2+2) + idpd1,1) = -(intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3)));
          }
        }
        //(dd|ss) integrals
        Dd[0] = Dvalue(atmC,6);
        intn[0] = d_eri2Center(2,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(30 + iddd,1) = intn[0]*DDrot_dA(iddd,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3)) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5));
        }
      }
    }
    if ((ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center(0,0,0,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center(0,0,0,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center(0,0,0,0,1,1,1,1,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(1,idx + 2) = -intn[0]*SProt_dA(idx + 1,1);}                                                         //(ss|sp) integrals
        block(1,5 + idx) = intn[1]*PProt_dA(idx + 1,1) + intn[2]*(PProt_dA(idx + 1,2) + PProt_dA(idx + 1,3));   //(ss|pp) integrals
      }
      if (ncols > 10) {
        //(ss|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(0,0,0,0,0,0,2,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(1,10 + idsd) = intn[0]*SDrot_dA(idsd,1);
        }
        //(ss|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(0,0,0,0,1,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center(0,0,0,0,1,1,2,1,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(1,5*(idpd2 + 2) + idpd1) = -intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1) - intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3));
          }
        }
        //(ss|dd) integrals
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(0,0,0,0,2,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center(0,0,0,0,2,1,2,1,RCD,Dd,atmC,atmD,type);
        intn[2] = d_eri2Center(0,0,0,0,2,2,2,2,RCD,Dd,atmC,atmD,type);
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(1,30 + iddd) = intn[0]*DDrot_dA(iddd,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3)) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5));
        }
      }
    }
    if ((nrows > 1)&&(ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center(0,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center(0,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center(0,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[3] = eri2Center(0,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);
      intn[4] = eri2Center(0,0,1,1,1,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idpb = 1; idpb < 4; ++idpb) {
        for (size_t idpk = 1; idpk < 4; ++idpk) {           //(sp|sp)
          block(idpb + 1,idpk + 1) = intn[0]*SProt_dA(idpb,1)*SProt(idpk,1) + intn[1]*(SProt_dA(idpb,2)*SProt(idpk,2) + SProt_dA(idpb,3)*SProt(idpk,3));
          block(idpb + 1,idpk + 1) += intn[0]*SProt(idpb,1)*SProt_dA(idpk,1) + intn[1]*(SProt(idpb,2)*SProt_dA(idpk,2) + SProt(idpb,3)*SProt_dA(idpk,3));
        }
        for (size_t idc = 1; idc < 7; ++idc) {              //(sp|pp)
          block(idpb + 1,idc + 4) = -SProt_dA(idpb,1)*(intn[2]*PProt(idc,1) + intn[3]*(PProt(idc,2) + PProt(idc,3))) - intn[4]*(SProt_dA(idpb,2)*PProt(idc,4) + SProt_dA(idpb,3)*PProt(idc,5));
          block(idpb + 1,idc + 4) += -SProt(idpb,1)*(intn[2]*PProt_dA(idc,1) + intn[3]*(PProt_dA(idc,2) + PProt_dA(idc,3))) - intn[4]*(SProt(idpb,2)*PProt_dA(idc,4) + SProt(idpb,3)*PProt_dA(idc,5));
        }
      }
      intn[0] = eri2Center(1,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);                          //(ps,ps|ps,ps)
      intn[1] = eri2Center(1,1,1,1,1,0,1,0,RCD,D,atmC,atmD,type);                          //(pp,pp|ps,ps)
      intn[2] = eri2Center(1,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);                          //(ps,ps|pp,pp)
      intn[3] = eri2Center(1,1,1,1,1,1,1,1,RCD,D,atmC,atmD,type);                          //(pp,pp|pp,pp)
      intn[4] = eri2Center(1,1,1,0,1,1,1,0,RCD,D,atmC,atmD,type);                          //(pp,ps|pp,ps)
      intn[5] = eri2Center(1,1,1,1,1,-1,1,-1,RCD,D,atmC,atmD,type);                        //(pp,pp|pp*,pp*)
      intn[6] = eri2Center(1,-1,1,1,1,-1,1,1,RCD,D,atmC,atmD,type);                        //(pp*,pp|pp*,pp)
      intn[7] = eri2Center(1,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[8] = eri2Center(1,1,1,1,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[9] = eri2Center(1,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idr = 1; idr < 7; ++idr) {
        for (size_t idc = 1; idc < 4; ++idc) {                //(pp|sp)
          block(idr + 4,idc + 1) = -SProt_dA(idc,1)*(intn[7]*PProt(idr,1) + intn[8]*(PProt(idr,2) + PProt(idr,3))) - intn[9]*(SProt_dA(idc,2)*PProt(idr,4) + SProt_dA(idc,3)*PProt(idr,5));
          block(idr + 4,idc + 1) += -SProt(idc,1)*(intn[7]*PProt_dA(idr,1) + intn[8]*(PProt_dA(idr,2) + PProt_dA(idr,3))) - intn[9]*(SProt(idc,2)*PProt_dA(idr,4) + SProt(idc,3)*PProt_dA(idr,5));
        }
        for (size_t idc = 1; idc < 7; ++idc) {                //(pp|pp)
          block(idr + 4,idc + 4) = intn[0]*PProt_dA(idr,1)*PProt(idc,1) + intn[1]*(PProt_dA(idr,2) + PProt_dA(idr,3))*PProt(idc,1) + intn[2]*PProt_dA(idr,1)*(PProt(idc,2) + PProt(idc,3)) + intn[3]*(PProt_dA(idr,2)*PProt(idc,2) + PProt_dA(idr,3)*PProt(idc,3)) + intn[4]*(PProt_dA(idr,4)*PProt(idc,4) + PProt_dA(idr,5)*PProt(idc,5)) + intn[5]*(PProt_dA(idr,2)*PProt(idc,3) + PProt_dA(idr,3)*PProt(idc,2)) + intn[6]*PProt_dA(idr,6)*PProt(idc,6);
          block(idr + 4,idc + 4) += intn[0]*PProt(idr,1)*PProt_dA(idc,1) + intn[1]*(PProt(idr,2) + PProt(idr,3))*PProt_dA(idc,1) + intn[2]*PProt(idr,1)*(PProt_dA(idc,2) + PProt_dA(idc,3)) + intn[3]*(PProt(idr,2)*PProt_dA(idc,2) + PProt(idr,3)*PProt_dA(idc,3)) + intn[4]*(PProt(idr,4)*PProt_dA(idc,4) + PProt(idr,5)*PProt_dA(idc,5)) + intn[5]*(PProt(idr,2)*PProt_dA(idc,3) + PProt(idr,3)*PProt_dA(idc,2)) + intn[6]*PProt(idr,6)*PProt_dA(idc,6);
        }
      }
      if (nrows > 10) {
        //(sd|sp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(0,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sd|sp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idsp = 1; idsp < 4; ++idsp) {
            block(10 + idsd,1 + idsp) = -intn[0]*SDrot_dA(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot_dA(idsd,2)*SProt(idsp,2) + SDrot_dA(idsd,3)*SProt(idsp,3));
            block(10 + idsd,1 + idsp) += -intn[0]*SDrot(idsd,1)*SProt_dA(idsp,1) - intn[1]*(SDrot(idsd,2)*SProt_dA(idsp,2) + SDrot(idsd,3)*SProt_dA(idsp,3));
          }
        }
        //(sd|pp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center(0,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(0,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(sd|pp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(10 + idsd,4 + idpp) = intn[0]*SDrot_dA(idsd,1)*PProt(idpp,1) + intn[1]*SDrot_dA(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(10 + idsd,4 + idpp) += intn[0]*SDrot(idsd,1)*PProt_dA(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot_dA(idsd,2)*PProt(idpp,4) + SDrot_dA(idsd,3)*PProt(idpp,5));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot(idsd,2)*PProt_dA(idpp,4) + SDrot(idsd,3)*PProt_dA(idpp,5));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot_dA(idsd,4)*PProt(idpp,2) - SDrot_dA(idsd,4)*PProt(idpp,3) + SDrot_dA(idsd,5)*PProt(idpp,6));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot(idsd,4)*PProt_dA(idpp,2) - SDrot(idsd,4)*PProt_dA(idpp,3) + SDrot(idsd,5)*PProt_dA(idpp,6));
          }
        }
        //(pd|sp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(1,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(1,1,2,0,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        intn[4] = d_eri2Center(1,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,1 + idsp) = intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt_dA(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt_dA(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt_dA(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt_dA(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt_dA(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt_dA(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt_dA(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt_dA(idsp,3));
            }
          }
        }
        //(dd|sp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(2,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center(2,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(2,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(30 + iddd,1 + idsp) = -intn[0]*DDrot_dA(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SProt(idsp,1);
            block(30 + iddd,1 + idsp) += -intn[0]*DDrot(iddd,1)*SProt_dA(idsp,1) - intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt_dA(idsp,1) - intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt_dA(idsp,1);
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot_dA(iddd,6)*SProt(idsp,2) + DDrot_dA(iddd,7)*SProt(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot(iddd,6)*SProt_dA(idsp,2) + DDrot(iddd,7)*SProt_dA(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SProt(idsp,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SProt(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt_dA(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt_dA(idsp,3));
          }
        }
        //(pd|pp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center(1,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center(1,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center(1,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center(1,1,2,0,1,1,1,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center(1,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|pp) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,4 + idpp) = -intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) += -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt_dA(idpp,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt_dA(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt_dA(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt_dA(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt_dA(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt_dA(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt_dA(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt_dA(idpp,5));
            }
          }
        }
        //(dd|pp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center(2,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center(2,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[4] = d_eri2Center(2,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center(2,2,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center(2,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(2,-1,2,-1,1,1,1,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center(2,1,2,-1,1,1,1,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center(2,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=(sg,-dl|pi,-pi)=-(sg,dl|-pi,-pi)
        intn[10] = d_eri2Center(2,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|pp) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(30 + iddd,4 + idpp) = intn[0]*DDrot_dA(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*PProt(idpp,1);
            block(30 + iddd,4 + idpp) += intn[0]*DDrot(iddd,1)*PProt_dA(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt_dA(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt_dA(idpp,1);
            block(30 + iddd,4 + idpp) += intn[3]*DDrot_dA(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot_dA(iddd,2)*PProt(idpp,2) + DDrot_dA(iddd,3)*PProt(idpp,3));
            block(30 + iddd,4 + idpp) += intn[3]*DDrot(iddd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt_dA(idpp,2) + DDrot(iddd,3)*PProt_dA(idpp,3));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot_dA(iddd,6)*PProt(idpp,4) + DDrot_dA(iddd,7)*PProt(idpp,5));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt_dA(idpp,4) + DDrot(iddd,7)*PProt_dA(idpp,5));
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot_dA(iddd,2)*PProt(idpp,3) + DDrot_dA(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot_dA(iddd,10)*PProt(idpp,6);
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot(iddd,2)*PProt_dA(idpp,3) + DDrot(iddd,3)*PProt_dA(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt_dA(idpp,6);
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot_dA(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot_dA(iddd,9)*PProt(idpp,6));
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot(iddd,8)*(PProt_dA(idpp,2) - PProt_dA(idpp,3)) + DDrot(iddd,9)*PProt_dA(idpp,6));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*PProt(idpp,4) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*PProt(idpp,5));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt_dA(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt_dA(idpp,5));
          }
        }
      }
      if (ncols > 10) {
        //(sp|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(0,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);   //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sp|sd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(1 + idsp,10 + idsd) = -intn[0]*SDrot_dA(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot_dA(idsd,2)*SProt(idsp,2) + SDrot_dA(idsd,3)*SProt(idsp,3));
            block(1 + idsp,10 + idsd) -= intn[0]*SDrot(idsd,1)*SProt_dA(idsp,1) + intn[1]*(SDrot(idsd,2)*SProt_dA(idsp,2) + SDrot(idsd,3)*SProt_dA(idsp,3));
          }
        }
        //(pp|sd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(1,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,1,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(1,1,1,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        //(pp|sd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,10 + idsd) = intn[0]*SDrot_dA(idsd,1)*PProt(idpp,1) + intn[1]*SDrot_dA(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(4 + idpp,10 + idsd) += intn[0]*SDrot(idsd,1)*PProt_dA(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot_dA(idsd,2)*PProt(idpp,4) + SDrot_dA(idsd,3)*PProt(idpp,5));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot(idsd,2)*PProt_dA(idpp,4) + SDrot(idsd,3)*PProt_dA(idpp,5));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot_dA(idsd,4)*PProt(idpp,2) - SDrot_dA(idsd,4)*PProt(idpp,3) + SDrot_dA(idsd,5)*PProt(idpp,6));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot(idsd,4)*PProt_dA(idpp,2) - SDrot(idsd,4)*PProt_dA(idpp,3) + SDrot(idsd,5)*PProt_dA(idpp,6));
          }
        }
        //(sp|pd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(0,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(0,0,1,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[4] = d_eri2Center(0,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|pd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(1 + idsp,5*(idpd2+2) + idpd1) = intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt_dA(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt_dA(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt_dA(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt_dA(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt_dA(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt_dA(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt_dA(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt_dA(idsp,3));
            }
          }
        }
        //(sp|dd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(0,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(0,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(0,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|dd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(1 + idsp,30 + iddd) = -intn[0]*DDrot_dA(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SProt(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[0]*DDrot(iddd,1)*SProt_dA(idsp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt_dA(idsp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt_dA(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot_dA(iddd,6)*SProt(idsp,2) + DDrot_dA(iddd,7)*SProt(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot(iddd,6)*SProt_dA(idsp,2) + DDrot(iddd,7)*SProt_dA(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SProt(idsp,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SProt(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt_dA(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt_dA(idsp,3));
          }
        }
        //(pp|pd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(1,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,1,1,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center(1,1,1,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center(1,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center(1,1,1,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center(1,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|pd) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(4 + idpp,5*(idpd2+2) + idpd1) = -intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt_dA(idpp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt_dA(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt_dA(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt_dA(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt_dA(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt_dA(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt_dA(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt_dA(idpp,5));
            }
          }
        }
        //(pp|dd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(1,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(1,1,1,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center(1,1,1,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center(1,1,1,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)=(-pi,-pi|-dl,-dl)
        intn[6] = d_eri2Center(1,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(1,-1,1,-1,2,1,2,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center(1,1,1,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center(1,1,1,1,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=(pi,-pi|sg,-dl)=-(-pi,-pi|sg,dl)
        intn[10] = d_eri2Center(1,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,30 + iddd) = intn[0]*DDrot_dA(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*PProt(idpp,1);
            block(4 + idpp,30 + iddd) += intn[0]*DDrot(iddd,1)*PProt_dA(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt_dA(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt_dA(idpp,1);
            block(4 + idpp,30 + iddd) += intn[3]*DDrot_dA(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot_dA(iddd,2)*PProt(idpp,2) + DDrot_dA(iddd,3)*PProt(idpp,3));
            block(4 + idpp,30 + iddd) += intn[3]*DDrot(iddd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt_dA(idpp,2) + DDrot(iddd,3)*PProt_dA(idpp,3));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot_dA(iddd,6)*PProt(idpp,4) + DDrot_dA(iddd,7)*PProt(idpp,5));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt_dA(idpp,4) + DDrot(iddd,7)*PProt_dA(idpp,5));
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot_dA(iddd,2)*PProt(idpp,3) + DDrot_dA(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot_dA(iddd,10)*PProt(idpp,6);
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot(iddd,2)*PProt_dA(idpp,3) + DDrot(iddd,3)*PProt_dA(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt_dA(idpp,6);
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot_dA(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot_dA(iddd,9)*PProt(idpp,6));
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot(iddd,8)*(PProt_dA(idpp,2) - PProt_dA(idpp,3)) + DDrot(iddd,9)*PProt_dA(idpp,6));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*PProt(idpp,4) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*PProt(idpp,5));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt_dA(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt_dA(idpp,5));
          }
        }
      }
      if ((nrows > 10)&&(ncols > 10)) {
        //(sd|sd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(0,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[2] = d_eri2Center(0,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        //(sd|sd) rotate
        for (size_t idsdb = 1; idsdb < 6; ++idsdb) {
          for (size_t idsdk = 1; idsdk < 6; ++idsdk) {
            block(10 + idsdb,10 + idsdk) = intn[0]*SDrot_dA(idsdb,1)*SDrot(idsdk,1) + intn[1]*(SDrot_dA(idsdb,2)*SDrot(idsdk,2) + SDrot_dA(idsdb,3)*SDrot(idsdk,3)) + intn[2]*(SDrot_dA(idsdb,4)*SDrot(idsdk,4) + SDrot_dA(idsdb,5)*SDrot(idsdk,5));
            block(10 + idsdb,10 + idsdk) += intn[0]*SDrot(idsdb,1)*SDrot_dA(idsdk,1) + intn[1]*(SDrot(idsdb,2)*SDrot_dA(idsdk,2) + SDrot(idsdb,3)*SDrot_dA(idsdk,3)) + intn[2]*(SDrot(idsdb,4)*SDrot_dA(idsdk,4) + SDrot(idsdb,5)*SDrot_dA(idsdk,5));
          }
        }
        //(sd|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(0,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(0,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(sg,pi|-pi,-dl)
        intn[4] = d_eri2Center(0,0,2,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        //(sd|pd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(10 + idsd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot_dA(idsd,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot_dA(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot_dA(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot_dA(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot_dA(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot_dA(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot_dA(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot_dA(idsd,3));
            }
          }
        }
        //(pd|sd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(1,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(1,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[4] = d_eri2Center(1,1,2,0,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        //(pd|sd) rotate
        for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
            for (size_t idsd = 1; idsd < 6; ++idsd) {
              block(5*(idpd2+2) + idpd1,10 + idsd) = -intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot_dA(idsd,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot_dA(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot_dA(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot_dA(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot_dA(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot_dA(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot_dA(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot_dA(idsd,3));
            }
          }
        }
        //(dd|sd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(2,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center(2,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(2,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center(2,1,2,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        intn[6] = d_eri2Center(2,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)
        //(dd|sd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(30 + iddd,10 + idsd) = intn[0]*DDrot_dA(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot_dA(idsd,1);
            block(30 + iddd,10 + idsd) += intn[0]*DDrot(iddd,1)*SDrot_dA(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot_dA(idsd,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SDrot(idsd,1);
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot_dA(iddd,6)*SDrot(idsd,2) + DDrot_dA(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot_dA(iddd,8)*SDrot(idsd,4) + DDrot_dA(iddd,9)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot(iddd,6)*SDrot_dA(idsd,2) + DDrot(iddd,7)*SDrot_dA(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot_dA(idsd,4) + DDrot(iddd,9)*SDrot_dA(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot_dA(iddd,2) - DDrot_dA(iddd,3))*SDrot(idsd,4) + DDrot_dA(iddd,10)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot_dA(idsd,4) + DDrot(iddd,10)*SDrot_dA(idsd,5));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SDrot(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot_dA(idsd,3));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot_dA(idsd,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SDrot(idsd,3));
          }
        }
        //(sd|dd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(0,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(0,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(0,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center(0,0,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        intn[6] = d_eri2Center(0,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        //(sd|dd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(10 + idsd,30 + iddd) = intn[0]*DDrot_dA(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot_dA(idsd,1);
            block(10 + idsd,30 + iddd) += intn[0]*DDrot(iddd,1)*SDrot_dA(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot_dA(idsd,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SDrot(idsd,1);
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot_dA(iddd,6)*SDrot(idsd,2) + DDrot_dA(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot_dA(iddd,8)*SDrot(idsd,4) + DDrot_dA(iddd,9)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot(iddd,6)*SDrot_dA(idsd,2) + DDrot(iddd,7)*SDrot_dA(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot_dA(idsd,4) + DDrot(iddd,9)*SDrot_dA(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot_dA(iddd,2) - DDrot_dA(iddd,3))*SDrot(idsd,4) + DDrot_dA(iddd,10)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot_dA(idsd,4) + DDrot(iddd,10)*SDrot_dA(idsd,5));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SDrot(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot_dA(idsd,3));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot_dA(idsd,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SDrot(idsd,3));
          }
        }
        //(pd|pd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(1,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)=(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|-pi,-pi)
        intn[3] = d_eri2Center(1,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(1,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[5] = d_eri2Center(1,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(pi,-dl|pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)
        intn[6] = d_eri2Center(1,1,2,0,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)=(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[7] = d_eri2Center(1,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[8] = d_eri2Center(1,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)=(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        //(pd|pd) rotate
        for (size_t idpdb1 = 1; idpdb1 < 6; ++idpdb1) {            //d orbital on bra
          for (size_t idpdb2 = 1; idpdb2 < 4; ++idpdb2) {          //p orbital on bra
            for (size_t idpdk1 = 1; idpdk1 < 6; ++idpdk1) {        //d orbital on ket
              for (size_t idpdk2 = 1; idpdk2 < 4; ++idpdk2) {      //p orbital on ket
                cnt1 = (idpdb1 - 1)*3 + idpdb2;
                cnt2 = (idpdk1 - 1)*3 + idpdk2;
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) = intn[0]*PDrot_dA(cnt1,1)*PDrot(cnt2,1) + intn[1]*((PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*PDrot(cnt2,1) + PDrot_dA(cnt1,1)*(PDrot(cnt2,2) + PDrot(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[0]*PDrot(cnt1,1)*PDrot_dA(cnt2,1) + intn[1]*((PDrot(cnt1,2) + PDrot(cnt1,3))*PDrot_dA(cnt2,1) + PDrot(cnt1,1)*(PDrot_dA(cnt2,2) + PDrot_dA(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(PDrot(cnt2,2) + PDrot(cnt2,3)) + intn[3]*(PDrot_dA(cnt1,4)*PDrot(cnt2,4) + PDrot_dA(cnt1,5)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(PDrot_dA(cnt2,2) + PDrot_dA(cnt2,3)) + intn[3]*(PDrot(cnt1,4)*PDrot_dA(cnt2,4) + PDrot(cnt1,5)*PDrot_dA(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot_dA(cnt1,8)*PDrot(cnt2,8) + PDrot_dA(cnt1,12)*PDrot(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot(cnt1,8)*PDrot_dA(cnt2,8) + PDrot(cnt1,12)*PDrot_dA(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot_dA(cnt1,14) - PDrot_dA(cnt1,11))*(PDrot(cnt2,14) - PDrot(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(PDrot_dA(cnt2,10) + PDrot_dA(cnt2,15)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(PDrot_dA(cnt2,14) - PDrot_dA(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot_dA(cnt1,4)*PDrot(cnt2,8) + PDrot_dA(cnt1,8)*PDrot(cnt2,4) + PDrot_dA(cnt1,5)*PDrot(cnt2,12) + PDrot_dA(cnt1,12)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot(cnt1,4)*PDrot_dA(cnt2,8) + PDrot(cnt1,8)*PDrot_dA(cnt2,4) + PDrot(cnt1,5)*PDrot_dA(cnt2,12) + PDrot(cnt1,12)*PDrot_dA(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*PDrot(cnt2,4) + PDrot_dA(cnt1,4)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*PDrot(cnt2,5) + PDrot_dA(cnt1,5)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot_dA(cnt2,4) + PDrot(cnt1,4)*(PDrot_dA(cnt2,10) + PDrot_dA(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot_dA(cnt2,5) + PDrot(cnt1,5)*(PDrot_dA(cnt2,11) - PDrot_dA(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*PDrot(cnt2,8) + PDrot_dA(cnt1,8)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*PDrot(cnt2,12) + PDrot_dA(cnt1,12)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot_dA(cnt2,8) + PDrot(cnt1,8)*(PDrot_dA(cnt2,10) + PDrot_dA(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot_dA(cnt2,12) + PDrot(cnt1,12)*(PDrot_dA(cnt2,11) - PDrot_dA(cnt2,14)));
              }
            }
          }
        }
        //(pd|dd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(1,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center(1,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center(1,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[5] = d_eri2Center(1,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|-dl,-dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)
        intn[6] = d_eri2Center(1,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(1,1,2,0,2,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center(1,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=(pi,-dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center(1,1,2,0,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)
        intn[10] = d_eri2Center(1,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        intn[11] = d_eri2Center(1,1,2,2,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(5*(idpd2+2) + idpd1,30 + iddd) = -intn[0]*PDrot_dA(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot_dA(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[0]*PDrot(cnt1,1)*DDrot_dA(iddd,1) + intn[1]*PDrot(cnt1,1)*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot_dA(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot_dA(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot(cnt1,1)*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot_dA(cnt1,4)*DDrot(iddd,6) + PDrot_dA(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot_dA(cnt1,8)*DDrot(iddd,6) + PDrot_dA(cnt1,12)*DDrot(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot(cnt1,4)*DDrot_dA(iddd,6) + PDrot(cnt1,5)*DDrot_dA(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot_dA(iddd,6) + PDrot(cnt1,12)*DDrot_dA(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot_dA(cnt1,14) - PDrot_dA(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot_dA(iddd,13) - DDrot_dA(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot_dA(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot_dA(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot(cnt1,8)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot(cnt1,12)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot_dA(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,5)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot(cnt1,4)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot_dA(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*DDrot(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot_dA(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot_dA(iddd,6) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*DDrot(iddd,7));
            }
          }
        }
        //(dd|pd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(2,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center(2,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center(2,2,2,2,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center(2,2,2,2,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|-pi,-pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)
        intn[6] = d_eri2Center(2,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(2,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center(2,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center(2,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        intn[10] = d_eri2Center(2,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[11] = d_eri2Center(2,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(dd|pd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(30 + iddd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot_dA(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot_dA(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[0]*PDrot(cnt1,1)*DDrot_dA(iddd,1) + intn[1]*PDrot(cnt1,1)*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot_dA(iddd,1) + intn[3]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot_dA(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot(cnt1,1)*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot_dA(cnt1,4)*DDrot(iddd,6) + PDrot_dA(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot_dA(cnt1,8)*DDrot(iddd,6) + PDrot_dA(cnt1,12)*DDrot(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot(cnt1,4)*DDrot_dA(iddd,6) + PDrot(cnt1,5)*DDrot_dA(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot_dA(iddd,6) + PDrot(cnt1,12)*DDrot_dA(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot_dA(cnt1,14) - PDrot_dA(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot_dA(iddd,13) - DDrot_dA(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot_dA(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,12)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot(cnt1,8)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot_dA(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot_dA(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,5)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot(cnt1,4)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot_dA(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*DDrot(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot_dA(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot_dA(iddd,6) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*DDrot(iddd,7));
            }
          }
        }
        //(dd|dd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(2,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(2,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(2,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center(2,2,2,2,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center(2,2,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center(2,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|-dl,-dl)
        intn[7] = d_eri2Center(2,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[8] = d_eri2Center(2,2,2,2,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(dl,dl|dl,dl)=(-dl,-dl|-dl,-dl)
        intn[9] = d_eri2Center(2,1,2,1,2,-1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[10] = d_eri2Center(2,-2,2,-2,2,2,2,2,RCD,Dd,atmC,atmD,type); //(-dl,-dl|dl,dl)=(dl,dl|-dl,-dl)
        intn[11] = d_eri2Center(2,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[12] = d_eri2Center(2,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[13] = d_eri2Center(2,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)
        intn[14] = d_eri2Center(2,1,2,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type); //(pi,-pi|pi,-pi)
        intn[15] = d_eri2Center(2,1,2,-2,2,0,2,-1,RCD,Dd,atmC,atmD,type); //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)=(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        intn[16] = d_eri2Center(2,1,2,1,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)=(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(dd|dd) rotate
        for (size_t iddd1 = 1; iddd1 < 16; ++iddd1) {
          for (size_t iddd2 = 1; iddd2 < 16; ++iddd2) {
            block(30 + iddd1,30 + iddd2) = intn[0]*DDrot_dA(iddd1,1)*DDrot(iddd2,1) + intn[1]*DDrot_dA(iddd1,1)*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[2]*DDrot(iddd1,1)*(DDrot_dA(iddd2,4) + DDrot_dA(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[0]*DDrot(iddd1,1)*DDrot_dA(iddd2,1) + intn[1]*DDrot(iddd1,1)*(DDrot_dA(iddd2,2) + DDrot_dA(iddd2,3)) + intn[2]*DDrot_dA(iddd1,1)*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot_dA(iddd1,2) + DDrot_dA(iddd1,3))*DDrot(iddd2,1) + intn[4]*(DDrot(iddd1,4) + DDrot(iddd1,5))*DDrot_dA(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot(iddd1,2) + DDrot(iddd1,3))*DDrot_dA(iddd2,1) + intn[4]*(DDrot_dA(iddd1,4) + DDrot_dA(iddd1,5))*DDrot(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot_dA(iddd1,4) + DDrot_dA(iddd1,5))*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[6]*(DDrot_dA(iddd1,2) + DDrot_dA(iddd1,3))*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot(iddd1,4) + DDrot(iddd1,5))*(DDrot_dA(iddd2,2) + DDrot_dA(iddd2,3)) + intn[6]*(DDrot(iddd1,2) + DDrot(iddd1,3))*(DDrot_dA(iddd2,4) + DDrot_dA(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot_dA(iddd1,2)*DDrot(iddd2,2) + DDrot_dA(iddd1,3)*DDrot(iddd2,3)) + intn[8]*(DDrot_dA(iddd1,4)*DDrot(iddd2,4) + DDrot_dA(iddd1,5)*DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot(iddd1,2)*DDrot_dA(iddd2,2) + DDrot(iddd1,3)*DDrot_dA(iddd2,3)) + intn[8]*(DDrot(iddd1,4)*DDrot_dA(iddd2,4) + DDrot(iddd1,5)*DDrot_dA(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot_dA(iddd1,2)*DDrot(iddd2,3) + DDrot_dA(iddd1,3)*DDrot(iddd2,2)) + intn[10]*(DDrot_dA(iddd1,4)*DDrot(iddd2,5) + DDrot_dA(iddd1,5)*DDrot(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot(iddd1,2)*DDrot_dA(iddd2,3) + DDrot(iddd1,3)*DDrot_dA(iddd2,2)) + intn[10]*(DDrot(iddd1,4)*DDrot_dA(iddd2,5) + DDrot(iddd1,5)*DDrot_dA(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot_dA(iddd1,6)*DDrot(iddd2,6) + DDrot_dA(iddd1,7)*DDrot(iddd2,7)) + intn[12]*(DDrot_dA(iddd1,8)*DDrot(iddd2,8) + DDrot_dA(iddd1,9)*DDrot(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot(iddd1,6)*DDrot_dA(iddd2,6) + DDrot(iddd1,7)*DDrot_dA(iddd2,7)) + intn[12]*(DDrot(iddd1,8)*DDrot_dA(iddd2,8) + DDrot(iddd1,9)*DDrot_dA(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot_dA(iddd1,11) + DDrot_dA(iddd1,14))*(DDrot(iddd2,11) + DDrot(iddd2,14)) + (DDrot_dA(iddd1,12) - DDrot_dA(iddd1,13))*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot(iddd1,11) + DDrot(iddd1,14))*(DDrot_dA(iddd2,11) + DDrot_dA(iddd2,14)) + (DDrot(iddd1,12) - DDrot(iddd1,13))*(DDrot_dA(iddd2,12) - DDrot_dA(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[14]*DDrot_dA(iddd1,10)*DDrot(iddd2,10) + intn[14]*DDrot(iddd1,10)*DDrot_dA(iddd2,10);
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot_dA(iddd1,11) + DDrot_dA(iddd1,14))*DDrot(iddd2,6) + (DDrot(iddd1,12) - DDrot(iddd1,13))*DDrot_dA(iddd2,7) + DDrot_dA(iddd1,6)*(DDrot(iddd2,11) + DDrot(iddd2,14)) + DDrot(iddd1,7)*(DDrot_dA(iddd2,12) - DDrot_dA(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot(iddd1,11) + DDrot(iddd1,14))*DDrot_dA(iddd2,6) + (DDrot_dA(iddd1,12) - DDrot_dA(iddd1,13))*DDrot(iddd2,7) + DDrot(iddd1,6)*(DDrot_dA(iddd2,11) + DDrot_dA(iddd2,14)) + DDrot_dA(iddd1,7)*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot_dA(iddd1,2) - DDrot_dA(iddd1,3))*DDrot(iddd2,8) + DDrot_dA(iddd1,10)*DDrot(iddd2,9) + DDrot_dA(iddd1,8)*(DDrot(iddd2,2) - DDrot(iddd2,3)) + DDrot_dA(iddd1,9)*DDrot(iddd2,10));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot(iddd1,2) - DDrot(iddd1,3))*DDrot_dA(iddd2,8) + DDrot(iddd1,10)*DDrot_dA(iddd2,9) + DDrot(iddd1,8)*(DDrot_dA(iddd2,2) - DDrot_dA(iddd2,3)) + DDrot(iddd1,9)*DDrot_dA(iddd2,10));
          }
        }
      }
    }
  }
  void IntegralBlock2C_dR2(int type, matrixE & block, int atmC, int atmD, double RCD, double cost, double sint, double cosp, double sinp) {
    //function calculating a block of second-derivatives of two-center integrals with respect to the internuclear distance
    int cnt1;
    int cnt2;
    size_t nrows = block.rows();
    size_t ncols = block.cols();
    if ((nrows > 1)||(ncols > 1)) {
      SProt = SPtransf(cost,sint,cosp,sinp);
      PProt = PPtransf(cost,sint,cosp,sinp);
      if ((nrows > 10)||(ncols > 10)) {
        SDrot = SDtransf(cost,sint,cosp,sinp);
        PDrot = PDtransf(cost,sint,cosp,sinp);
        DDrot = DDtransf(cost,sint,cosp,sinp);
      }
    }
    D[0] = Dvalue(atmC,1);
    D[1] = Dvalue(atmC,2);
    D[2] = Dvalue(atmD,1);
    D[3] = Dvalue(atmD,2);
    //integral calculation
    block(1,1) = eri2Center_dR2(0,0,0,0,0,0,0,0,RCD,D,atmC,atmD,type);       //(ss|ss)
    if (nrows > 1) {
      intn[0] = eri2Center_dR2(0,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center_dR2(1,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center_dR2(1,1,1,1,0,0,0,0,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(idx + 2,1) = -intn[0]*SProt(idx + 1,1);}                                                   //(sp|ss) integrals
        block(5 + idx,1) = intn[1]*PProt(idx + 1,1) + intn[2]*(PProt(idx + 1,2) + PProt(idx + 1,3));   //(pp|ss) integrals
      }
      if (nrows > 10) {
        //(sd|ss) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR2(0,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(10 + idsd,1) = intn[0]*SDrot(idsd,1);
        }
        //(pd|ss) integrals
        Dd[0] = Dvalue(atmC,5);
        intn[0] = d_eri2Center_dR2(1,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center_dR2(1,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(5*(idpd2+2) + idpd1,1) = -(intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3)));
          }
        }
        //(dd|ss) integrals
        Dd[0] = Dvalue(atmC,6);
        intn[0] = d_eri2Center_dR2(2,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(2,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR2(2,2,2,2,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(30 + iddd,1) = intn[0]*DDrot(iddd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3)) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5));
        }
      }
    }
    if ((ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center_dR2(0,0,0,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center_dR2(0,0,0,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center_dR2(0,0,0,0,1,1,1,1,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(1,idx + 2) = -intn[0]*SProt(idx + 1,1);}                                                   //(ss|sp) integrals
        block(1,5 + idx) = intn[1]*PProt(idx + 1,1) + intn[2]*(PProt(idx + 1,2) + PProt(idx + 1,3));   //(ss|pp) integrals
      }
      if (ncols > 10) {
        //(ss|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR2(0,0,0,0,0,0,2,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(1,10 + idsd) = intn[0]*SDrot(idsd,1);
        }
        //(ss|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR2(0,0,0,0,1,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center_dR2(0,0,0,0,1,1,2,1,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(1,5*(idpd2 + 2) + idpd1) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3));
          }
        }
        //(ss|dd) integrals
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR2(0,0,0,0,2,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center_dR2(0,0,0,0,2,1,2,1,RCD,Dd,atmC,atmD,type);
        intn[2] = d_eri2Center_dR2(0,0,0,0,2,2,2,2,RCD,Dd,atmC,atmD,type);
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(1,30 + iddd) = intn[0]*DDrot(iddd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3)) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5));
        }
      }
    }
    if ((nrows > 1)&&(ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center_dR2(0,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center_dR2(0,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center_dR2(0,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[3] = eri2Center_dR2(0,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);
      intn[4] = eri2Center_dR2(0,0,1,1,1,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idpb = 1; idpb < 4; ++idpb) {
        for (size_t idpk = 1; idpk < 4; ++idpk) {           //(sp|sp)
          block(idpb + 1,idpk + 1) = intn[0]*SProt(idpb,1)*SProt(idpk,1) + intn[1]*(SProt(idpb,2)*SProt(idpk,2) + SProt(idpb,3)*SProt(idpk,3));
        }
        for (size_t idc = 1; idc < 7; ++idc) {              //(sp|pp)
          block(idpb + 1,idc + 4) = -SProt(idpb,1)*(intn[2]*PProt(idc,1) + intn[3]*(PProt(idc,2) + PProt(idc,3))) - intn[4]*(SProt(idpb,2)*PProt(idc,4) + SProt(idpb,3)*PProt(idc,5));
        }
      }
      intn[0] = eri2Center_dR2(1,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);                          //(ps,ps|ps,ps)
      intn[1] = eri2Center_dR2(1,1,1,1,1,0,1,0,RCD,D,atmC,atmD,type);                          //(pp,pp|ps,ps)
      intn[2] = eri2Center_dR2(1,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);                          //(ps,ps|pp,pp)
      intn[3] = eri2Center_dR2(1,1,1,1,1,1,1,1,RCD,D,atmC,atmD,type);                          //(pp,pp|pp,pp)
      intn[4] = eri2Center_dR2(1,1,1,0,1,1,1,0,RCD,D,atmC,atmD,type);                          //(pp,ps|pp,ps)
      intn[5] = eri2Center_dR2(1,1,1,1,1,-1,1,-1,RCD,D,atmC,atmD,type);                        //(pp,pp|pp*,pp*)
      intn[6] = eri2Center_dR2(1,-1,1,1,1,-1,1,1,RCD,D,atmC,atmD,type);                        //(pp*,pp|pp*,pp)
      intn[7] = eri2Center_dR2(1,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[8] = eri2Center_dR2(1,1,1,1,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[9] = eri2Center_dR2(1,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idr = 1; idr < 7; ++idr) {
        for (size_t idc = 1; idc < 4; ++idc) {                //(pp|sp)
          block(idr + 4,idc + 1) = -SProt(idc,1)*(intn[7]*PProt(idr,1) + intn[8]*(PProt(idr,2) + PProt(idr,3))) - intn[9]*(SProt(idc,2)*PProt(idr,4) + SProt(idc,3)*PProt(idr,5));
        }
        for (size_t idc = 1; idc < 7; ++idc) {                //(pp|pp)
          block(idr + 4,idc + 4) = intn[0]*PProt(idr,1)*PProt(idc,1) + intn[1]*(PProt(idr,2) + PProt(idr,3))*PProt(idc,1) + intn[2]*PProt(idr,1)*(PProt(idc,2) + PProt(idc,3)) + intn[3]*(PProt(idr,2)*PProt(idc,2) + PProt(idr,3)*PProt(idc,3)) + intn[4]*(PProt(idr,4)*PProt(idc,4) + PProt(idr,5)*PProt(idc,5)) + intn[5]*(PProt(idr,2)*PProt(idc,3) + PProt(idr,3)*PProt(idc,2)) + intn[6]*PProt(idr,6)*PProt(idc,6);
        }
      }
      if (nrows > 10) {
        //(sd|sp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR2(0,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(0,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sd|sp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idsp = 1; idsp < 4; ++idsp) {
            block(10 + idsd,1 + idsp) = -intn[0]*SDrot(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot(idsd,2)*SProt(idsp,2) + SDrot(idsd,3)*SProt(idsp,3));
          }
        }
        //(sd|pp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center_dR2(0,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(0,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR2(0,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR2(0,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(sd|pp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(10 + idsd,4 + idpp) = intn[0]*SDrot(idsd,1)*PProt(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot(idsd,2)*PProt(idpp,4) + SDrot(idsd,3)*PProt(idpp,5));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot(idsd,4)*PProt(idpp,2) - SDrot(idsd,4)*PProt(idpp,3) + SDrot(idsd,5)*PProt(idpp,6));
          }
        }
        //(pd|sp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR2(1,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(1,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR2(1,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR2(1,1,2,0,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        intn[4] = d_eri2Center_dR2(1,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,1 + idsp) = intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
            }
          }
        }
        //(dd|sp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR2(2,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(2,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR2(2,2,2,2,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center_dR2(2,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR2(2,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(30 + iddd,1 + idsp) = -intn[0]*DDrot(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt(idsp,1);
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot(iddd,6)*SProt(idsp,2) + DDrot(iddd,7)*SProt(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt(idsp,3));
          }
        }
        //(pd|pp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center_dR2(1,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(1,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR2(1,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center_dR2(1,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center_dR2(1,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center_dR2(1,1,2,0,1,1,1,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center_dR2(1,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|pp) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,4 + idpp) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
            }
          }
        }
        //(dd|pp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center_dR2(2,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(2,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR2(2,2,2,2,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center_dR2(2,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[4] = d_eri2Center_dR2(2,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center_dR2(2,2,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center_dR2(2,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR2(2,-1,2,-1,1,1,1,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center_dR2(2,1,2,-1,1,1,1,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center_dR2(2,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=(sg,-dl|pi,-pi)=-(sg,dl|-pi,-pi)
        intn[10] = d_eri2Center_dR2(2,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|pp) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(30 + iddd,4 + idpp) = intn[0]*DDrot(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt(idpp,1);
            block(30 + iddd,4 + idpp) += intn[3]*DDrot(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt(idpp,2) + DDrot(iddd,3)*PProt(idpp,3));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt(idpp,4) + DDrot(iddd,7)*PProt(idpp,5));
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot(iddd,2)*PProt(idpp,3) + DDrot(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt(idpp,6);
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot(iddd,9)*PProt(idpp,6));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt(idpp,5));
          }
        }
      }
      if (ncols > 10) {
        //(sp|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR2(0,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(0,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sp|sd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(1 + idsp,10 + idsd) = -intn[0]*SDrot(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot(idsd,2)*SProt(idsp,2) + SDrot(idsd,3)*SProt(idsp,3));
          }
        }
        //(pp|sd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR2(1,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(1,1,1,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR2(1,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR2(1,1,1,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        //(pp|sd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,10 + idsd) = intn[0]*SDrot(idsd,1)*PProt(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot(idsd,2)*PProt(idpp,4) + SDrot(idsd,3)*PProt(idpp,5));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot(idsd,4)*PProt(idpp,2) - SDrot(idsd,4)*PProt(idpp,3) + SDrot(idsd,5)*PProt(idpp,6));
          }
        }
        //(sp|pd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR2(0,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(0,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR2(0,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR2(0,0,1,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[4] = d_eri2Center_dR2(0,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|pd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(1 + idsp,5*(idpd2+2) + idpd1) = intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
            }
          }
        }
        //(sp|dd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR2(0,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(0,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR2(0,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR2(0,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR2(0,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|dd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(1 + idsp,30 + iddd) = -intn[0]*DDrot(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot(iddd,6)*SProt(idsp,2) + DDrot(iddd,7)*SProt(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt(idsp,3));
          }
        }
        //(pp|pd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR2(1,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(1,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR2(1,1,1,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center_dR2(1,1,1,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center_dR2(1,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center_dR2(1,1,1,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center_dR2(1,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|pd) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(4 + idpp,5*(idpd2+2) + idpd1) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
            }
          }
        }
        //(pp|dd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR2(1,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(1,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR2(1,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR2(1,1,1,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center_dR2(1,1,1,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center_dR2(1,1,1,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)=(-pi,-pi|-dl,-dl)
        intn[6] = d_eri2Center_dR2(1,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR2(1,-1,1,-1,2,1,2,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center_dR2(1,1,1,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center_dR2(1,1,1,1,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=(pi,-pi|sg,-dl)=-(-pi,-pi|sg,dl)
        intn[10] = d_eri2Center_dR2(1,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,30 + iddd) = intn[0]*DDrot(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt(idpp,1);
            block(4 + idpp,30 + iddd) += intn[3]*DDrot(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt(idpp,2) + DDrot(iddd,3)*PProt(idpp,3));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt(idpp,4) + DDrot(iddd,7)*PProt(idpp,5));
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot(iddd,2)*PProt(idpp,3) + DDrot(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt(idpp,6);
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot(iddd,9)*PProt(idpp,6));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt(idpp,5));
          }
        }
      }
      if ((nrows > 10)&&(ncols > 10)) {
        //(sd|sd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR2(0,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(0,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[2] = d_eri2Center_dR2(0,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        //(sd|sd) rotate
        for (size_t idsdb = 1; idsdb < 6; ++idsdb) {
          for (size_t idsdk = 1; idsdk < 6; ++idsdk) {
            block(10 + idsdb,10 + idsdk) = intn[0]*SDrot(idsdb,1)*SDrot(idsdk,1) + intn[1]*(SDrot(idsdb,2)*SDrot(idsdk,2) + SDrot(idsdb,3)*SDrot(idsdk,3)) + intn[2]*(SDrot(idsdb,4)*SDrot(idsdk,4) + SDrot(idsdb,5)*SDrot(idsdk,5));
          }
        }
        //(sd|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR2(0,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(0,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR2(0,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR2(0,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(sg,pi|-pi,-dl)
        intn[4] = d_eri2Center_dR2(0,0,2,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        //(sd|pd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(10 + idsd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
            }
          }
        }
        //(pd|sd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR2(1,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(1,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR2(1,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR2(1,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[4] = d_eri2Center_dR2(1,1,2,0,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        //(pd|sd) rotate
        for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
            for (size_t idsd = 1; idsd < 6; ++idsd) {
              block(5*(idpd2+2) + idpd1,10 + idsd) = -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
            }
          }
        }
        //(dd|sd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR2(2,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(2,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR2(2,2,2,2,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center_dR2(2,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR2(2,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center_dR2(2,1,2,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        intn[6] = d_eri2Center_dR2(2,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)
        //(dd|sd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(30 + iddd,10 + idsd) = intn[0]*DDrot(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot(idsd,1);
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot(iddd,6)*SDrot(idsd,2) + DDrot(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot(idsd,4) + DDrot(iddd,9)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot(idsd,4) + DDrot(iddd,10)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot(idsd,3));
          }
        }
        //(sd|dd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR2(0,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(0,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR2(0,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR2(0,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR2(0,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center_dR2(0,0,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        intn[6] = d_eri2Center_dR2(0,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        //(sd|dd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(10 + idsd,30 + iddd) = intn[0]*DDrot(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot(idsd,1);
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot(iddd,6)*SDrot(idsd,2) + DDrot(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot(idsd,4) + DDrot(iddd,9)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot(idsd,4) + DDrot(iddd,10)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot(idsd,3));
          }
        }
        //(pd|pd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR2(1,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(1,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)=(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR2(1,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|-pi,-pi)
        intn[3] = d_eri2Center_dR2(1,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR2(1,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[5] = d_eri2Center_dR2(1,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(pi,-dl|pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)
        intn[6] = d_eri2Center_dR2(1,1,2,0,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)=(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[7] = d_eri2Center_dR2(1,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[8] = d_eri2Center_dR2(1,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)=(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        //(pd|pd) rotate
        for (size_t idpdb1 = 1; idpdb1 < 6; ++idpdb1) {            //d orbital on bra
          for (size_t idpdb2 = 1; idpdb2 < 4; ++idpdb2) {          //p orbital on bra
            for (size_t idpdk1 = 1; idpdk1 < 6; ++idpdk1) {        //d orbital on ket
              for (size_t idpdk2 = 1; idpdk2 < 4; ++idpdk2) {      //p orbital on ket
                cnt1 = (idpdb1 - 1)*3 + idpdb2;
                cnt2 = (idpdk1 - 1)*3 + idpdk2;
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) = intn[0]*PDrot(cnt1,1)*PDrot(cnt2,1) + intn[1]*((PDrot(cnt1,2) + PDrot(cnt1,3))*PDrot(cnt2,1) + PDrot(cnt1,1)*(PDrot(cnt2,2) + PDrot(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(PDrot(cnt2,2) + PDrot(cnt2,3)) + intn[3]*(PDrot(cnt1,4)*PDrot(cnt2,4) + PDrot(cnt1,5)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot(cnt1,8)*PDrot(cnt2,8) + PDrot(cnt1,12)*PDrot(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(PDrot(cnt2,14) - PDrot(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot(cnt1,4)*PDrot(cnt2,8) + PDrot(cnt1,8)*PDrot(cnt2,4) + PDrot(cnt1,5)*PDrot(cnt2,12) + PDrot(cnt1,12)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot(cnt2,4) + PDrot(cnt1,4)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot(cnt2,5) + PDrot(cnt1,5)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot(cnt2,8) + PDrot(cnt1,8)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot(cnt2,12) + PDrot(cnt1,12)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
              }
            }
          }
        }
        //(pd|dd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR2(1,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(1,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR2(1,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center_dR2(1,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center_dR2(1,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[5] = d_eri2Center_dR2(1,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|-dl,-dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)
        intn[6] = d_eri2Center_dR2(1,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR2(1,1,2,0,2,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center_dR2(1,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=(pi,-dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center_dR2(1,1,2,0,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)
        intn[10] = d_eri2Center_dR2(1,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        intn[11] = d_eri2Center_dR2(1,1,2,2,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(5*(idpd2+2) + idpd1,30 + iddd) = -intn[0]*PDrot(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot(cnt1,4)*DDrot(iddd,6) + PDrot(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot(iddd,6) + PDrot(cnt1,12)*DDrot(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot(iddd,7));
            }
          }
        }
        //(dd|pd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR2(2,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(2,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR2(2,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center_dR2(2,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center_dR2(2,2,2,2,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center_dR2(2,2,2,2,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|-pi,-pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)
        intn[6] = d_eri2Center_dR2(2,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR2(2,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center_dR2(2,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center_dR2(2,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        intn[10] = d_eri2Center_dR2(2,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[11] = d_eri2Center_dR2(2,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(dd|pd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(30 + iddd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot(cnt1,4)*DDrot(iddd,6) + PDrot(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot(iddd,6) + PDrot(cnt1,12)*DDrot(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot(iddd,7));
            }
          }
        }
        //(dd|dd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR2(2,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR2(2,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR2(2,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR2(2,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center_dR2(2,2,2,2,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center_dR2(2,2,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center_dR2(2,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|-dl,-dl)
        intn[7] = d_eri2Center_dR2(2,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[8] = d_eri2Center_dR2(2,2,2,2,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(dl,dl|dl,dl)=(-dl,-dl|-dl,-dl)
        intn[9] = d_eri2Center_dR2(2,1,2,1,2,-1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[10] = d_eri2Center_dR2(2,-2,2,-2,2,2,2,2,RCD,Dd,atmC,atmD,type); //(-dl,-dl|dl,dl)=(dl,dl|-dl,-dl)
        intn[11] = d_eri2Center_dR2(2,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[12] = d_eri2Center_dR2(2,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[13] = d_eri2Center_dR2(2,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)
        intn[14] = d_eri2Center_dR2(2,1,2,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type); //(pi,-pi|pi,-pi)
        intn[15] = d_eri2Center_dR2(2,1,2,-2,2,0,2,-1,RCD,Dd,atmC,atmD,type); //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)=(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        intn[16] = d_eri2Center_dR2(2,1,2,1,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)=(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(dd|dd) rotate
        for (size_t iddd1 = 1; iddd1 < 16; ++iddd1) {
          for (size_t iddd2 = 1; iddd2 < 16; ++iddd2) {
            block(30 + iddd1,30 + iddd2) = intn[0]*DDrot(iddd1,1)*DDrot(iddd2,1) + intn[1]*DDrot(iddd1,1)*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[2]*DDrot(iddd1,1)*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot(iddd1,2) + DDrot(iddd1,3))*DDrot(iddd2,1) + intn[4]*(DDrot(iddd1,4) + DDrot(iddd1,5))*DDrot(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot(iddd1,4) + DDrot(iddd1,5))*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[6]*(DDrot(iddd1,2) + DDrot(iddd1,3))*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot(iddd1,2)*DDrot(iddd2,2) + DDrot(iddd1,3)*DDrot(iddd2,3)) + intn[8]*(DDrot(iddd1,4)*DDrot(iddd2,4) + DDrot(iddd1,5)*DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot(iddd1,2)*DDrot(iddd2,3) + DDrot(iddd1,3)*DDrot(iddd2,2)) + intn[10]*(DDrot(iddd1,4)*DDrot(iddd2,5) + DDrot(iddd1,5)*DDrot(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot(iddd1,6)*DDrot(iddd2,6) + DDrot(iddd1,7)*DDrot(iddd2,7)) + intn[12]*(DDrot(iddd1,8)*DDrot(iddd2,8) + DDrot(iddd1,9)*DDrot(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot(iddd1,11) + DDrot(iddd1,14))*(DDrot(iddd2,11) + DDrot(iddd2,14)) + (DDrot(iddd1,12) - DDrot(iddd1,13))*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[14]*DDrot(iddd1,10)*DDrot(iddd2,10);
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot(iddd1,11) + DDrot(iddd1,14))*DDrot(iddd2,6) + (DDrot(iddd1,12) - DDrot(iddd1,13))*DDrot(iddd2,7) + DDrot(iddd1,6)*(DDrot(iddd2,11) + DDrot(iddd2,14)) + DDrot(iddd1,7)*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot(iddd1,2) - DDrot(iddd1,3))*DDrot(iddd2,8) + DDrot(iddd1,10)*DDrot(iddd2,9) + DDrot(iddd1,8)*(DDrot(iddd2,2) - DDrot(iddd2,3)) + DDrot(iddd1,9)*DDrot(iddd2,10));
          }
        }
      }
    }
  }
  void IntegralBlock2C_dAB(int type, matrixE & block, int derA, int derB, int atmC, int atmD, double RCD, double cost, double sint, double cosp, double sinp) {
    //function calculating a block of second-derivatives with respect to pair of angular variables of two-center integrals
    //2 -> theta; 3 -> phi
    matrixE SProt_dA(1,1);
    matrixE PProt_dA(1,1);
    matrixE SDrot_dA(1,1);
    matrixE PDrot_dA(1,1);
    matrixE DDrot_dA(1,1);
    matrixE SProt_dB(1,1);
    matrixE PProt_dB(1,1);
    matrixE SDrot_dB(1,1);
    matrixE PDrot_dB(1,1);
    matrixE DDrot_dB(1,1);
    matrixE SProt_dAB(1,1);
    matrixE PProt_dAB(1,1);
    matrixE SDrot_dAB(1,1);
    matrixE PDrot_dAB(1,1);
    matrixE DDrot_dAB(1,1);
    int cnt1;
    int cnt2;
    size_t nrows = block.rows();
    size_t ncols = block.cols();
    if ((nrows > 1)||(ncols > 1)) {
      SProt = SPtransf(cost,sint,cosp,sinp);
      PProt = PPtransf(cost,sint,cosp,sinp);
      if ((nrows > 10)||(ncols > 10)) {
        SDrot = SDtransf(cost,sint,cosp,sinp);
        PDrot = PDtransf(cost,sint,cosp,sinp);
        DDrot = DDtransf(cost,sint,cosp,sinp);
      }
      if (derA == 2) {
        SProt_dA = SPtransf_dt(cost,sint,cosp,sinp);
        PProt_dA = PPtransf_dt(cost,sint,cosp,sinp);
        if ((nrows > 10)||(ncols > 10)) {
          SDrot_dA = SDtransf_dt(cost,sint,cosp,sinp);
          PDrot_dA = PDtransf_dt(cost,sint,cosp,sinp);
          DDrot_dA = DDtransf_dt(cost,sint,cosp,sinp);
        }
      }
      else if (derA == 3) {
        SProt_dA = SPtransf_dp(cost,sint,cosp,sinp);
        PProt_dA = PPtransf_dp(cost,sint,cosp,sinp);
        if ((nrows > 10)||(ncols > 10)) {
          SDrot_dA = SDtransf_dp(cost,sint,cosp,sinp);
          PDrot_dA = PDtransf_dp(cost,sint,cosp,sinp);
          DDrot_dA = DDtransf_dp(cost,sint,cosp,sinp);
        }
      }
      if (derB == 2) {
        SProt_dB = SPtransf_dt(cost,sint,cosp,sinp);
        PProt_dB = PPtransf_dt(cost,sint,cosp,sinp);
        if ((nrows > 10)||(ncols > 10)) {
          SDrot_dB = SDtransf_dt(cost,sint,cosp,sinp);
          PDrot_dB = PDtransf_dt(cost,sint,cosp,sinp);
          DDrot_dB = DDtransf_dt(cost,sint,cosp,sinp);
        }
      }
      else if (derB == 3) {
        SProt_dB = SPtransf_dp(cost,sint,cosp,sinp);
        PProt_dB = PPtransf_dp(cost,sint,cosp,sinp);
        if ((nrows > 10)||(ncols > 10)) {
          SDrot_dB = SDtransf_dp(cost,sint,cosp,sinp);
          PDrot_dB = PDtransf_dp(cost,sint,cosp,sinp);
          DDrot_dB = DDtransf_dp(cost,sint,cosp,sinp);
        }
      }
      if (derA == derB) {
        if (derA == 2) {
          SProt_dAB = SPtransf_dt2(cost,sint,cosp,sinp);
          PProt_dAB = PPtransf_dt2(cost,sint,cosp,sinp);
          if ((nrows > 10)||(ncols > 10)) {
            SDrot_dAB = SDtransf_dt2(cost,sint,cosp,sinp);
            PDrot_dAB = PDtransf_dt2(cost,sint,cosp,sinp);
            DDrot_dAB = DDtransf_dt2(cost,sint,cosp,sinp);
          }
        }
        else if (derA == 3) {
          SProt_dAB = SPtransf_dp2(cost,sint,cosp,sinp);
          PProt_dAB = PPtransf_dp2(cost,sint,cosp,sinp);
          if ((nrows > 10)||(ncols > 10)) {
            SDrot_dAB = SDtransf_dp2(cost,sint,cosp,sinp);
            PDrot_dAB = PDtransf_dp2(cost,sint,cosp,sinp);
            DDrot_dAB = DDtransf_dp2(cost,sint,cosp,sinp);
          }
        }
      }
      else {
        SProt_dAB = SPtransf_dtdp(cost,sint,cosp,sinp);
        PProt_dAB = PPtransf_dtdp(cost,sint,cosp,sinp);
        if ((nrows > 10)||(ncols > 10)) {
          SDrot_dAB = SDtransf_dtdp(cost,sint,cosp,sinp);
          PDrot_dAB = PDtransf_dtdp(cost,sint,cosp,sinp);
          DDrot_dAB = DDtransf_dtdp(cost,sint,cosp,sinp);
        }
      }
    }
    D[0] = Dvalue(atmC,1);
    D[1] = Dvalue(atmC,2);
    D[2] = Dvalue(atmD,1);
    D[3] = Dvalue(atmD,2);
    //integral calculation
    block(1,1) = 0.0;                                                                                                                      //(ss|ss)
    if (nrows > 1) {
      intn[0] = eri2Center(0,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center(1,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center(1,1,1,1,0,0,0,0,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(idx + 2,1) = -intn[0]*SProt_dAB(idx + 1,1);}                                                           //(sp|ss) integrals
        block(5 + idx,1) = intn[1]*PProt_dAB(idx + 1,1) + intn[2]*(PProt_dAB(idx + 1,2) + PProt_dAB(idx + 1,3));   //(pp|ss) integrals
      }
      if (nrows > 10) {
        //(sd|ss) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(0,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(10 + idsd,1) = intn[0]*SDrot_dAB(idsd,1);
        }
        //(pd|ss) integrals
        Dd[0] = Dvalue(atmC,5);
        intn[0] = d_eri2Center(1,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center(1,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(5*(idpd2+2) + idpd1,1) = -(intn[0]*PDrot_dAB((idpd1 - 1)*3 + idpd2,1) + intn[1]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,3)));
          }
        }
        //(dd|ss) integrals
        Dd[0] = Dvalue(atmC,6);
        intn[0] = d_eri2Center(2,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(30 + iddd,1) = intn[0]*DDrot_dAB(iddd,1) + intn[1]*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3)) + intn[2]*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5));
        }
      }
    }
    if ((ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center(0,0,0,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center(0,0,0,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center(0,0,0,0,1,1,1,1,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(1,idx + 2) = -intn[0]*SProt_dAB(idx + 1,1);}                                                           //(ss|sp) integrals
        block(1,5 + idx) = intn[1]*PProt_dAB(idx + 1,1) + intn[2]*(PProt_dAB(idx + 1,2) + PProt_dAB(idx + 1,3));   //(ss|pp) integrals
      }
      if (ncols > 10) {
        //(ss|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(0,0,0,0,0,0,2,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(1,10 + idsd) = intn[0]*SDrot_dAB(idsd,1);
        }
        //(ss|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(0,0,0,0,1,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center(0,0,0,0,1,1,2,1,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(1,5*(idpd2 + 2) + idpd1) = -intn[0]*PDrot_dAB((idpd1 - 1)*3 + idpd2,1) - intn[1]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,3));
          }
        }
        //(ss|dd) integrals
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(0,0,0,0,2,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center(0,0,0,0,2,1,2,1,RCD,Dd,atmC,atmD,type);
        intn[2] = d_eri2Center(0,0,0,0,2,2,2,2,RCD,Dd,atmC,atmD,type);
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(1,30 + iddd) = intn[0]*DDrot_dAB(iddd,1) + intn[1]*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3)) + intn[2]*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5));
        }
      }
    }
    if ((nrows > 1)&&(ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center(0,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center(0,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center(0,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[3] = eri2Center(0,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);
      intn[4] = eri2Center(0,0,1,1,1,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idpb = 1; idpb < 4; ++idpb) {
        for (size_t idpk = 1; idpk < 4; ++idpk) {           //(sp|sp)
          block(idpb + 1,idpk + 1) = intn[0]*SProt_dAB(idpb,1)*SProt(idpk,1) + intn[1]*(SProt_dAB(idpb,2)*SProt(idpk,2) + SProt_dAB(idpb,3)*SProt(idpk,3));
          block(idpb + 1,idpk + 1) += intn[0]*SProt_dA(idpb,1)*SProt_dB(idpk,1) + intn[1]*(SProt_dA(idpb,2)*SProt_dB(idpk,2) + SProt_dA(idpb,3)*SProt_dB(idpk,3));
          block(idpb + 1,idpk + 1) += intn[0]*SProt_dB(idpb,1)*SProt_dA(idpk,1) + intn[1]*(SProt_dB(idpb,2)*SProt_dA(idpk,2) + SProt_dB(idpb,3)*SProt_dA(idpk,3));
          block(idpb + 1,idpk + 1) += intn[0]*SProt(idpb,1)*SProt_dAB(idpk,1) + intn[1]*(SProt(idpb,2)*SProt_dAB(idpk,2) + SProt(idpb,3)*SProt_dAB(idpk,3));
        }
        for (size_t idc = 1; idc < 7; ++idc) {              //(sp|pp)
          block(idpb + 1,idc + 4) = -SProt_dAB(idpb,1)*(intn[2]*PProt(idc,1) + intn[3]*(PProt(idc,2) + PProt(idc,3))) - intn[4]*(SProt_dAB(idpb,2)*PProt(idc,4) + SProt_dAB(idpb,3)*PProt(idc,5));
          block(idpb + 1,idc + 4) -= SProt_dA(idpb,1)*(intn[2]*PProt_dB(idc,1) + intn[3]*(PProt_dB(idc,2) + PProt_dB(idc,3))) + intn[4]*(SProt_dA(idpb,2)*PProt_dB(idc,4) + SProt_dA(idpb,3)*PProt_dB(idc,5));
          block(idpb + 1,idc + 4) -= SProt_dB(idpb,1)*(intn[2]*PProt_dA(idc,1) + intn[3]*(PProt_dA(idc,2) + PProt_dA(idc,3))) + intn[4]*(SProt_dB(idpb,2)*PProt_dA(idc,4) + SProt_dB(idpb,3)*PProt_dA(idc,5));
          block(idpb + 1,idc + 4) -= SProt(idpb,1)*(intn[2]*PProt_dAB(idc,1) + intn[3]*(PProt_dAB(idc,2) + PProt_dAB(idc,3))) + intn[4]*(SProt(idpb,2)*PProt_dAB(idc,4) + SProt(idpb,3)*PProt_dAB(idc,5));
        }
      }
      intn[0] = eri2Center(1,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);                          //(ps,ps|ps,ps)
      intn[1] = eri2Center(1,1,1,1,1,0,1,0,RCD,D,atmC,atmD,type);                          //(pp,pp|ps,ps)
      intn[2] = eri2Center(1,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);                          //(ps,ps|pp,pp)
      intn[3] = eri2Center(1,1,1,1,1,1,1,1,RCD,D,atmC,atmD,type);                          //(pp,pp|pp,pp)
      intn[4] = eri2Center(1,1,1,0,1,1,1,0,RCD,D,atmC,atmD,type);                          //(pp,ps|pp,ps)
      intn[5] = eri2Center(1,1,1,1,1,-1,1,-1,RCD,D,atmC,atmD,type);                        //(pp,pp|pp*,pp*)
      intn[6] = eri2Center(1,-1,1,1,1,-1,1,1,RCD,D,atmC,atmD,type);                        //(pp*,pp|pp*,pp)
      intn[7] = eri2Center(1,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[8] = eri2Center(1,1,1,1,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[9] = eri2Center(1,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idr = 1; idr < 7; ++idr) {
        for (size_t idc = 1; idc < 4; ++idc) {                //(pp|sp)
          block(idr + 4,idc + 1) = -SProt_dAB(idc,1)*(intn[7]*PProt(idr,1) + intn[8]*(PProt(idr,2) + PProt(idr,3))) - intn[9]*(SProt_dAB(idc,2)*PProt(idr,4) + SProt_dAB(idc,3)*PProt(idr,5));
          block(idr + 4,idc + 1) -= SProt_dA(idc,1)*(intn[7]*PProt_dB(idr,1) + intn[8]*(PProt_dB(idr,2) + PProt_dB(idr,3))) + intn[9]*(SProt_dA(idc,2)*PProt_dB(idr,4) + SProt_dA(idc,3)*PProt_dB(idr,5));
          block(idr + 4,idc + 1) -= SProt_dB(idc,1)*(intn[7]*PProt_dA(idr,1) + intn[8]*(PProt_dA(idr,2) + PProt_dA(idr,3))) + intn[9]*(SProt_dB(idc,2)*PProt_dA(idr,4) + SProt_dB(idc,3)*PProt_dA(idr,5));
          block(idr + 4,idc + 1) -= SProt(idc,1)*(intn[7]*PProt_dAB(idr,1) + intn[8]*(PProt_dAB(idr,2) + PProt_dAB(idr,3))) + intn[9]*(SProt(idc,2)*PProt_dAB(idr,4) + SProt(idc,3)*PProt_dAB(idr,5));
        }
        for (size_t idc = 1; idc < 7; ++idc) {                //(pp|pp)
          block(idr + 4,idc + 4) = intn[0]*PProt_dAB(idr,1)*PProt(idc,1) + intn[1]*(PProt_dAB(idr,2) + PProt_dAB(idr,3))*PProt(idc,1) + intn[2]*PProt_dAB(idr,1)*(PProt(idc,2) + PProt(idc,3)) + intn[3]*(PProt_dAB(idr,2)*PProt(idc,2) + PProt_dAB(idr,3)*PProt(idc,3)) + intn[4]*(PProt_dAB(idr,4)*PProt(idc,4) + PProt_dAB(idr,5)*PProt(idc,5)) + intn[5]*(PProt_dAB(idr,2)*PProt(idc,3) + PProt_dAB(idr,3)*PProt(idc,2)) + intn[6]*PProt_dAB(idr,6)*PProt(idc,6);
          block(idr + 4,idc + 4) += intn[0]*PProt_dA(idr,1)*PProt_dB(idc,1) + intn[1]*(PProt_dA(idr,2) + PProt_dA(idr,3))*PProt_dB(idc,1) + intn[2]*PProt_dA(idr,1)*(PProt_dB(idc,2) + PProt_dB(idc,3)) + intn[3]*(PProt_dA(idr,2)*PProt_dB(idc,2) + PProt_dA(idr,3)*PProt_dB(idc,3)) + intn[4]*(PProt_dA(idr,4)*PProt_dB(idc,4) + PProt_dA(idr,5)*PProt_dB(idc,5)) + intn[5]*(PProt_dA(idr,2)*PProt_dB(idc,3) + PProt_dA(idr,3)*PProt_dB(idc,2)) + intn[6]*PProt_dA(idr,6)*PProt_dB(idc,6);
          block(idr + 4,idc + 4) += intn[0]*PProt_dB(idr,1)*PProt_dA(idc,1) + intn[1]*(PProt_dB(idr,2) + PProt_dB(idr,3))*PProt_dA(idc,1) + intn[2]*PProt_dB(idr,1)*(PProt_dA(idc,2) + PProt_dA(idc,3)) + intn[3]*(PProt_dB(idr,2)*PProt_dA(idc,2) + PProt_dB(idr,3)*PProt_dA(idc,3)) + intn[4]*(PProt_dB(idr,4)*PProt_dA(idc,4) + PProt_dB(idr,5)*PProt_dA(idc,5)) + intn[5]*(PProt_dB(idr,2)*PProt_dA(idc,3) + PProt_dB(idr,3)*PProt_dA(idc,2)) + intn[6]*PProt_dB(idr,6)*PProt_dA(idc,6);
          block(idr + 4,idc + 4) += intn[0]*PProt(idr,1)*PProt_dAB(idc,1) + intn[1]*(PProt(idr,2) + PProt(idr,3))*PProt_dAB(idc,1) + intn[2]*PProt(idr,1)*(PProt_dAB(idc,2) + PProt_dAB(idc,3)) + intn[3]*(PProt(idr,2)*PProt_dAB(idc,2) + PProt(idr,3)*PProt_dAB(idc,3)) + intn[4]*(PProt(idr,4)*PProt_dAB(idc,4) + PProt(idr,5)*PProt_dAB(idc,5)) + intn[5]*(PProt(idr,2)*PProt_dAB(idc,3) + PProt(idr,3)*PProt_dAB(idc,2)) + intn[6]*PProt(idr,6)*PProt_dAB(idc,6);
        }
      }
      if (nrows > 10) {
        //(sd|sp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(0,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sd|sp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idsp = 1; idsp < 4; ++idsp) {
            block(10 + idsd,1 + idsp) = -intn[0]*SDrot_dAB(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot_dAB(idsd,2)*SProt(idsp,2) + SDrot_dAB(idsd,3)*SProt(idsp,3));
            block(10 + idsd,1 + idsp) -= intn[0]*SDrot_dA(idsd,1)*SProt_dB(idsp,1) + intn[1]*(SDrot_dA(idsd,2)*SProt_dB(idsp,2) + SDrot_dA(idsd,3)*SProt_dB(idsp,3));
            block(10 + idsd,1 + idsp) -= intn[0]*SDrot_dB(idsd,1)*SProt_dA(idsp,1) + intn[1]*(SDrot_dB(idsd,2)*SProt_dA(idsp,2) + SDrot_dB(idsd,3)*SProt_dA(idsp,3));
            block(10 + idsd,1 + idsp) -= intn[0]*SDrot(idsd,1)*SProt_dAB(idsp,1) + intn[1]*(SDrot(idsd,2)*SProt_dAB(idsp,2) + SDrot(idsd,3)*SProt_dAB(idsp,3));
          }
        }
        //(sd|pp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center(0,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(0,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(sd|pp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(10 + idsd,4 + idpp) = intn[0]*SDrot_dAB(idsd,1)*PProt(idpp,1) + intn[1]*SDrot_dAB(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(10 + idsd,4 + idpp) += intn[0]*SDrot_dA(idsd,1)*PProt_dB(idpp,1) + intn[1]*SDrot_dA(idsd,1)*(PProt_dB(idpp,2) + PProt_dB(idpp,3));
            block(10 + idsd,4 + idpp) += intn[0]*SDrot_dB(idsd,1)*PProt_dA(idpp,1) + intn[1]*SDrot_dB(idsd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
            block(10 + idsd,4 + idpp) += intn[0]*SDrot(idsd,1)*PProt_dAB(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt_dAB(idpp,2) + PProt_dAB(idpp,3));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot_dAB(idsd,2)*PProt(idpp,4) + SDrot_dAB(idsd,3)*PProt(idpp,5));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot_dA(idsd,2)*PProt_dB(idpp,4) + SDrot_dA(idsd,3)*PProt_dB(idpp,5));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot_dB(idsd,2)*PProt_dA(idpp,4) + SDrot_dB(idsd,3)*PProt_dA(idpp,5));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot(idsd,2)*PProt_dAB(idpp,4) + SDrot(idsd,3)*PProt_dAB(idpp,5));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot_dAB(idsd,4)*PProt(idpp,2) - SDrot_dAB(idsd,4)*PProt(idpp,3) + SDrot_dAB(idsd,5)*PProt(idpp,6));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot_dA(idsd,4)*PProt_dB(idpp,2) - SDrot_dA(idsd,4)*PProt_dB(idpp,3) + SDrot_dA(idsd,5)*PProt_dB(idpp,6));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot_dB(idsd,4)*PProt_dA(idpp,2) - SDrot_dB(idsd,4)*PProt_dA(idpp,3) + SDrot_dB(idsd,5)*PProt_dA(idpp,6));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot(idsd,4)*PProt_dAB(idpp,2) - SDrot(idsd,4)*PProt_dAB(idpp,3) + SDrot(idsd,5)*PProt_dAB(idpp,6));
          }
        }
        //(pd|sp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(1,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(1,1,2,0,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        intn[4] = d_eri2Center(1,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,1 + idsp) = intn[0]*PDrot_dAB((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SProt_dB(idsp,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SProt_dB(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[0]*PDrot_dB((idpd1 - 1)*3 + idpd2,1)*SProt_dA(idsp,1) + intn[1]*(PDrot_dB((idpd1 - 1)*3 + idpd2,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,3))*SProt_dA(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt_dAB(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt_dAB(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SProt_dB(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SProt_dB(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot_dB((idpd1 - 1)*3 + idpd2,4)*SProt_dA(idsp,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,5)*SProt_dA(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt_dAB(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt_dAB(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SProt_dB(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SProt_dB(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot_dB((idpd1 - 1)*3 + idpd2,8)*SProt_dA(idsp,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,12)*SProt_dA(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt_dAB(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt_dAB(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot_dAB((idpd1 - 1)*3 + idpd2,10) + PDrot_dAB((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot_dAB((idpd1 - 1)*3 + idpd2,11) - PDrot_dAB((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SProt_dB(idsp,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SProt_dB(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot_dB((idpd1 - 1)*3 + idpd2,10) + PDrot_dB((idpd1 - 1)*3 + idpd2,15))*SProt_dA(idsp,2) + (PDrot_dB((idpd1 - 1)*3 + idpd2,11) - PDrot_dB((idpd1 - 1)*3 + idpd2,14))*SProt_dA(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt_dAB(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt_dAB(idsp,3));
            }
          }
        }
        //(dd|sp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center(2,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center(2,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(2,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(30 + iddd,1 + idsp) = -intn[0]*DDrot_dAB(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5))*SProt(idsp,1);
            block(30 + iddd,1 + idsp) -= intn[0]*DDrot_dA(iddd,1)*SProt_dB(idsp,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SProt_dB(idsp,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SProt_dB(idsp,1);
            block(30 + iddd,1 + idsp) -= intn[0]*DDrot_dB(iddd,1)*SProt_dA(idsp,1) + intn[1]*(DDrot_dB(iddd,2) + DDrot_dB(iddd,3))*SProt_dA(idsp,1) + intn[2]*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5))*SProt_dA(idsp,1);
            block(30 + iddd,1 + idsp) -= intn[0]*DDrot(iddd,1)*SProt_dAB(idsp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt_dAB(idsp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt_dAB(idsp,1);
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot_dAB(iddd,6)*SProt(idsp,2) + DDrot_dAB(iddd,7)*SProt(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot_dA(iddd,6)*SProt_dB(idsp,2) + DDrot_dA(iddd,7)*SProt_dB(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot_dB(iddd,6)*SProt_dA(idsp,2) + DDrot_dB(iddd,7)*SProt_dA(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot(iddd,6)*SProt_dAB(idsp,2) + DDrot(iddd,7)*SProt_dAB(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14))*SProt(idsp,2) + (DDrot_dAB(iddd,12) - DDrot_dAB(iddd,13))*SProt(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SProt_dB(idsp,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SProt_dB(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot_dB(iddd,11) + DDrot_dB(iddd,14))*SProt_dA(idsp,2) + (DDrot_dB(iddd,12) - DDrot_dB(iddd,13))*SProt_dA(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt_dAB(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt_dAB(idsp,3));
          }
        }
        //(pd|pp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center(1,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center(1,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center(1,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center(1,1,2,0,1,1,1,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center(1,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|pp) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,4 + idpp) = -intn[0]*PDrot_dAB((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*PProt_dB(idpp,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*PProt_dB(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[0]*PDrot_dB((idpd1 - 1)*3 + idpd2,1)*PProt_dA(idpp,1) + intn[1]*(PDrot_dB((idpd1 - 1)*3 + idpd2,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,3))*PProt_dA(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt_dAB(idpp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt_dAB(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot_dAB((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*(PProt_dB(idpp,2) + PProt_dB(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot_dB((idpd1 - 1)*3 + idpd2,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt_dAB(idpp,2) + PProt_dAB(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*(PProt_dB(idpp,2) + PProt_dB(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot_dB((idpd1 - 1)*3 + idpd2,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,3))*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt_dAB(idpp,2) + PProt_dAB(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot_dAB((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*PProt_dB(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*PProt_dB(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot_dB((idpd1 - 1)*3 + idpd2,4)*PProt_dA(idpp,4) + PDrot_dB((idpd1 - 1)*3 + idpd2,5)*PProt_dA(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt_dAB(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt_dAB(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot_dAB((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*PProt_dB(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*PProt_dB(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot_dB((idpd1 - 1)*3 + idpd2,8)*PProt_dA(idpp,4) + PDrot_dB((idpd1 - 1)*3 + idpd2,12)*PProt_dA(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt_dAB(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt_dAB(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot_dAB((idpd1 - 1)*3 + idpd2,10) + PDrot_dAB((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot_dAB((idpd1 - 1)*3 + idpd2,11) - PDrot_dAB((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*PProt_dB(idpp,4) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*PProt_dB(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot_dB((idpd1 - 1)*3 + idpd2,10) + PDrot_dB((idpd1 - 1)*3 + idpd2,15))*PProt_dA(idpp,4) + (PDrot_dB((idpd1 - 1)*3 + idpd2,11) - PDrot_dB((idpd1 - 1)*3 + idpd2,14))*PProt_dA(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt_dAB(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt_dAB(idpp,5));
            }
          }
        }
        //(dd|pp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center(2,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center(2,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[4] = d_eri2Center(2,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center(2,2,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center(2,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(2,-1,2,-1,1,1,1,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center(2,1,2,-1,1,1,1,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center(2,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=(sg,-dl|pi,-pi)=-(sg,dl|-pi,-pi)
        intn[10] = d_eri2Center(2,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|pp) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(30 + iddd,4 + idpp) = intn[0]*DDrot_dAB(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5))*PProt(idpp,1);
            block(30 + iddd,4 + idpp) += intn[0]*DDrot_dA(iddd,1)*PProt_dB(idpp,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*PProt_dB(idpp,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*PProt_dB(idpp,1);
            block(30 + iddd,4 + idpp) += intn[0]*DDrot_dB(iddd,1)*PProt_dA(idpp,1) + intn[1]*(DDrot_dB(iddd,2) + DDrot_dB(iddd,3))*PProt_dA(idpp,1) + intn[2]*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5))*PProt_dA(idpp,1);
            block(30 + iddd,4 + idpp) += intn[0]*DDrot(iddd,1)*PProt_dAB(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt_dAB(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt_dAB(idpp,1);
            block(30 + iddd,4 + idpp) += intn[3]*DDrot_dAB(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot_dAB(iddd,2)*PProt(idpp,2) + DDrot_dAB(iddd,3)*PProt(idpp,3));
            block(30 + iddd,4 + idpp) += intn[3]*DDrot_dA(iddd,1)*(PProt_dB(idpp,2) + PProt_dB(idpp,3)) + intn[4]*(DDrot_dA(iddd,2)*PProt_dB(idpp,2) + DDrot_dA(iddd,3)*PProt_dB(idpp,3));
            block(30 + iddd,4 + idpp) += intn[3]*DDrot_dB(iddd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[4]*(DDrot_dB(iddd,2)*PProt_dA(idpp,2) + DDrot_dB(iddd,3)*PProt_dA(idpp,3));
            block(30 + iddd,4 + idpp) += intn[3]*DDrot(iddd,1)*(PProt_dAB(idpp,2) + PProt_dAB(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt_dAB(idpp,2) + DDrot(iddd,3)*PProt_dAB(idpp,3));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot_dAB(iddd,6)*PProt(idpp,4) + DDrot_dAB(iddd,7)*PProt(idpp,5));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*(PProt_dB(idpp,2) + PProt_dB(idpp,3)) + intn[6]*(DDrot_dA(iddd,6)*PProt_dB(idpp,4) + DDrot_dA(iddd,7)*PProt_dB(idpp,5));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5))*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[6]*(DDrot_dB(iddd,6)*PProt_dA(idpp,4) + DDrot_dB(iddd,7)*PProt_dA(idpp,5));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt_dAB(idpp,2) + PProt_dAB(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt_dAB(idpp,4) + DDrot(iddd,7)*PProt_dAB(idpp,5));
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot_dAB(iddd,2)*PProt(idpp,3) + DDrot_dAB(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot_dAB(iddd,10)*PProt(idpp,6);
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot_dA(iddd,2)*PProt_dB(idpp,3) + DDrot_dA(iddd,3)*PProt_dB(idpp,2)) + intn[8]*DDrot_dA(iddd,10)*PProt_dB(idpp,6);
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot_dB(iddd,2)*PProt_dA(idpp,3) + DDrot_dB(iddd,3)*PProt_dA(idpp,2)) + intn[8]*DDrot_dB(iddd,10)*PProt_dA(idpp,6);
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot(iddd,2)*PProt_dAB(idpp,3) + DDrot(iddd,3)*PProt_dAB(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt_dAB(idpp,6);
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot_dAB(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot_dAB(iddd,9)*PProt(idpp,6));
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot_dA(iddd,8)*(PProt_dB(idpp,2) - PProt_dB(idpp,3)) + DDrot_dA(iddd,9)*PProt_dB(idpp,6));
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot_dB(iddd,8)*(PProt_dA(idpp,2) - PProt_dA(idpp,3)) + DDrot_dB(iddd,9)*PProt_dA(idpp,6));
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot(iddd,8)*(PProt_dAB(idpp,2) - PProt_dAB(idpp,3)) + DDrot(iddd,9)*PProt_dAB(idpp,6));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14))*PProt(idpp,4) + (DDrot_dAB(iddd,12) - DDrot_dAB(iddd,13))*PProt(idpp,5));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*PProt_dB(idpp,4) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*PProt_dB(idpp,5));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot_dB(iddd,11) + DDrot_dB(iddd,14))*PProt_dA(idpp,4) + (DDrot_dB(iddd,12) - DDrot_dB(iddd,13))*PProt_dA(idpp,5));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt_dAB(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt_dAB(idpp,5));
          }
        }
      }
      if (ncols > 10) {
        //(sp|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(0,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sp|sd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(1 + idsp,10 + idsd) = -intn[0]*SDrot_dAB(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot_dAB(idsd,2)*SProt(idsp,2) + SDrot_dAB(idsd,3)*SProt(idsp,3));
            block(1 + idsp,10 + idsd) -= intn[0]*SDrot_dA(idsd,1)*SProt_dB(idsp,1) + intn[1]*(SDrot_dA(idsd,2)*SProt_dB(idsp,2) + SDrot_dA(idsd,3)*SProt_dB(idsp,3));
            block(1 + idsp,10 + idsd) -= intn[0]*SDrot_dB(idsd,1)*SProt_dA(idsp,1) + intn[1]*(SDrot_dB(idsd,2)*SProt_dA(idsp,2) + SDrot_dB(idsd,3)*SProt_dA(idsp,3));
            block(1 + idsp,10 + idsd) -= intn[0]*SDrot(idsd,1)*SProt_dAB(idsp,1) + intn[1]*(SDrot(idsd,2)*SProt_dAB(idsp,2) + SDrot(idsd,3)*SProt_dAB(idsp,3));
          }
        }
        //(pp|sd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(1,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,1,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(1,1,1,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        //(pp|sd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,10 + idsd) = intn[0]*SDrot_dAB(idsd,1)*PProt(idpp,1) + intn[1]*SDrot_dAB(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(4 + idpp,10 + idsd) += intn[0]*SDrot_dA(idsd,1)*PProt_dB(idpp,1) + intn[1]*SDrot_dA(idsd,1)*(PProt_dB(idpp,2) + PProt_dB(idpp,3));
            block(4 + idpp,10 + idsd) += intn[0]*SDrot_dB(idsd,1)*PProt_dA(idpp,1) + intn[1]*SDrot_dB(idsd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
            block(4 + idpp,10 + idsd) += intn[0]*SDrot(idsd,1)*PProt_dAB(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt_dAB(idpp,2) + PProt_dAB(idpp,3));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot_dAB(idsd,2)*PProt(idpp,4) + SDrot_dAB(idsd,3)*PProt(idpp,5));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot_dA(idsd,2)*PProt_dB(idpp,4) + SDrot_dA(idsd,3)*PProt_dB(idpp,5));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot_dB(idsd,2)*PProt_dA(idpp,4) + SDrot_dB(idsd,3)*PProt_dA(idpp,5));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot(idsd,2)*PProt_dAB(idpp,4) + SDrot(idsd,3)*PProt_dAB(idpp,5));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot_dAB(idsd,4)*PProt(idpp,2) - SDrot_dAB(idsd,4)*PProt(idpp,3) + SDrot_dAB(idsd,5)*PProt(idpp,6));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot_dA(idsd,4)*PProt_dB(idpp,2) - SDrot_dA(idsd,4)*PProt_dB(idpp,3) + SDrot_dA(idsd,5)*PProt_dB(idpp,6));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot_dB(idsd,4)*PProt_dA(idpp,2) - SDrot_dB(idsd,4)*PProt_dA(idpp,3) + SDrot_dB(idsd,5)*PProt_dA(idpp,6));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot(idsd,4)*PProt_dAB(idpp,2) - SDrot(idsd,4)*PProt_dAB(idpp,3) + SDrot(idsd,5)*PProt_dAB(idpp,6));
          }
        }
        //(sp|pd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(0,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(0,0,1,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[4] = d_eri2Center(0,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|pd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(1 + idsp,5*(idpd2+2) + idpd1) = intn[0]*PDrot_dAB((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SProt_dB(idsp,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SProt_dB(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[0]*PDrot_dB((idpd1 - 1)*3 + idpd2,1)*SProt_dA(idsp,1) + intn[1]*(PDrot_dB((idpd1 - 1)*3 + idpd2,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,3))*SProt_dA(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt_dAB(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt_dAB(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SProt_dB(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SProt_dB(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot_dB((idpd1 - 1)*3 + idpd2,4)*SProt_dA(idsp,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,5)*SProt_dA(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt_dAB(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt_dAB(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SProt_dB(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SProt_dB(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot_dB((idpd1 - 1)*3 + idpd2,8)*SProt_dA(idsp,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,12)*SProt_dA(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt_dAB(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt_dAB(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot_dAB((idpd1 - 1)*3 + idpd2,10) + PDrot_dAB((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot_dAB((idpd1 - 1)*3 + idpd2,11) - PDrot_dAB((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SProt_dB(idsp,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SProt_dB(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot_dB((idpd1 - 1)*3 + idpd2,10) + PDrot_dB((idpd1 - 1)*3 + idpd2,15))*SProt_dA(idsp,2) + (PDrot_dB((idpd1 - 1)*3 + idpd2,11) - PDrot_dB((idpd1 - 1)*3 + idpd2,14))*SProt_dA(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt_dAB(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt_dAB(idsp,3));
            }
          }
        }
        //(sp|dd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(0,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(0,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(0,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|dd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(1 + idsp,30 + iddd) = -intn[0]*DDrot_dAB(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5))*SProt(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[0]*DDrot_dA(iddd,1)*SProt_dB(idsp,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SProt_dB(idsp,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SProt_dB(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[0]*DDrot_dB(iddd,1)*SProt_dA(idsp,1) + intn[1]*(DDrot_dB(iddd,2) + DDrot_dB(iddd,3))*SProt_dA(idsp,1) + intn[2]*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5))*SProt_dA(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[0]*DDrot(iddd,1)*SProt_dAB(idsp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt_dAB(idsp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt_dAB(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot_dAB(iddd,6)*SProt(idsp,2) + DDrot_dAB(iddd,7)*SProt(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot_dA(iddd,6)*SProt_dB(idsp,2) + DDrot_dA(iddd,7)*SProt_dB(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot_dB(iddd,6)*SProt_dA(idsp,2) + DDrot_dB(iddd,7)*SProt_dA(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot(iddd,6)*SProt_dAB(idsp,2) + DDrot(iddd,7)*SProt_dAB(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14))*SProt(idsp,2) + (DDrot_dAB(iddd,12) - DDrot_dAB(iddd,13))*SProt(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SProt_dB(idsp,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SProt_dB(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot_dB(iddd,11) + DDrot_dB(iddd,14))*SProt_dA(idsp,2) + (DDrot_dB(iddd,12) - DDrot_dB(iddd,13))*SProt_dA(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt_dAB(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt_dAB(idsp,3));
          }
        }
        //(pp|pd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(1,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,1,1,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center(1,1,1,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center(1,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center(1,1,1,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center(1,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|pd) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(4 + idpp,5*(idpd2+2) + idpd1) = -intn[0]*PDrot_dAB((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*PProt_dB(idpp,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*PProt_dB(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[0]*PDrot_dB((idpd1 - 1)*3 + idpd2,1)*PProt_dA(idpp,1) + intn[1]*(PDrot_dB((idpd1 - 1)*3 + idpd2,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,3))*PProt_dA(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt_dAB(idpp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt_dAB(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot_dAB((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*(PProt_dB(idpp,2) + PProt_dB(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot_dB((idpd1 - 1)*3 + idpd2,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt_dAB(idpp,2) + PProt_dAB(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*(PProt_dB(idpp,2) + PProt_dB(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot_dB((idpd1 - 1)*3 + idpd2,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,3))*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt_dAB(idpp,2) + PProt_dAB(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot_dAB((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*PProt_dB(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*PProt_dB(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot_dB((idpd1 - 1)*3 + idpd2,4)*PProt_dA(idpp,4) + PDrot_dB((idpd1 - 1)*3 + idpd2,5)*PProt_dA(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt_dAB(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt_dAB(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot_dAB((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*PProt_dB(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*PProt_dB(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot_dB((idpd1 - 1)*3 + idpd2,8)*PProt_dA(idpp,4) + PDrot_dB((idpd1 - 1)*3 + idpd2,12)*PProt_dA(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt_dAB(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt_dAB(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot_dAB((idpd1 - 1)*3 + idpd2,10) + PDrot_dAB((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot_dAB((idpd1 - 1)*3 + idpd2,11) - PDrot_dAB((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*PProt_dB(idpp,4) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*PProt_dB(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot_dB((idpd1 - 1)*3 + idpd2,10) + PDrot_dB((idpd1 - 1)*3 + idpd2,15))*PProt_dA(idpp,4) + (PDrot_dB((idpd1 - 1)*3 + idpd2,11) - PDrot_dB((idpd1 - 1)*3 + idpd2,14))*PProt_dA(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt_dAB(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt_dAB(idpp,5));
            }
          }
        }
        //(pp|dd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(1,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(1,1,1,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center(1,1,1,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center(1,1,1,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)=(-pi,-pi|-dl,-dl)
        intn[6] = d_eri2Center(1,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(1,-1,1,-1,2,1,2,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center(1,1,1,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center(1,1,1,1,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=(pi,-pi|sg,-dl)=-(-pi,-pi|sg,dl)
        intn[10] = d_eri2Center(1,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,30 + iddd) = intn[0]*DDrot_dAB(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5))*PProt(idpp,1);
            block(4 + idpp,30 + iddd) += intn[0]*DDrot_dA(iddd,1)*PProt_dB(idpp,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*PProt_dB(idpp,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*PProt_dB(idpp,1);
            block(4 + idpp,30 + iddd) += intn[0]*DDrot_dB(iddd,1)*PProt_dA(idpp,1) + intn[1]*(DDrot_dB(iddd,2) + DDrot_dB(iddd,3))*PProt_dA(idpp,1) + intn[2]*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5))*PProt_dA(idpp,1);
            block(4 + idpp,30 + iddd) += intn[0]*DDrot(iddd,1)*PProt_dAB(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt_dAB(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt_dAB(idpp,1);
            block(4 + idpp,30 + iddd) += intn[3]*DDrot_dAB(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot_dAB(iddd,2)*PProt(idpp,2) + DDrot_dAB(iddd,3)*PProt(idpp,3));
            block(4 + idpp,30 + iddd) += intn[3]*DDrot_dA(iddd,1)*(PProt_dB(idpp,2) + PProt_dB(idpp,3)) + intn[4]*(DDrot_dA(iddd,2)*PProt_dB(idpp,2) + DDrot_dA(iddd,3)*PProt_dB(idpp,3));
            block(4 + idpp,30 + iddd) += intn[3]*DDrot_dB(iddd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[4]*(DDrot_dB(iddd,2)*PProt_dA(idpp,2) + DDrot_dB(iddd,3)*PProt_dA(idpp,3));
            block(4 + idpp,30 + iddd) += intn[3]*DDrot(iddd,1)*(PProt_dAB(idpp,2) + PProt_dAB(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt_dAB(idpp,2) + DDrot(iddd,3)*PProt_dAB(idpp,3));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot_dAB(iddd,6)*PProt(idpp,4) + DDrot_dAB(iddd,7)*PProt(idpp,5));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*(PProt_dB(idpp,2) + PProt_dB(idpp,3)) + intn[6]*(DDrot_dA(iddd,6)*PProt_dB(idpp,4) + DDrot_dA(iddd,7)*PProt_dB(idpp,5));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5))*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[6]*(DDrot_dB(iddd,6)*PProt_dA(idpp,4) + DDrot_dB(iddd,7)*PProt_dA(idpp,5));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt_dAB(idpp,2) + PProt_dAB(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt_dAB(idpp,4) + DDrot(iddd,7)*PProt_dAB(idpp,5));
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot_dAB(iddd,2)*PProt(idpp,3) + DDrot_dAB(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot_dAB(iddd,10)*PProt(idpp,6);
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot_dA(iddd,2)*PProt_dB(idpp,3) + DDrot_dA(iddd,3)*PProt_dB(idpp,2)) + intn[8]*DDrot_dA(iddd,10)*PProt_dB(idpp,6);
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot_dB(iddd,2)*PProt_dA(idpp,3) + DDrot_dB(iddd,3)*PProt_dA(idpp,2)) + intn[8]*DDrot_dB(iddd,10)*PProt_dA(idpp,6);
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot(iddd,2)*PProt_dAB(idpp,3) + DDrot(iddd,3)*PProt_dAB(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt_dAB(idpp,6);
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot_dAB(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot_dAB(iddd,9)*PProt(idpp,6));
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot_dA(iddd,8)*(PProt_dB(idpp,2) - PProt_dB(idpp,3)) + DDrot_dA(iddd,9)*PProt_dB(idpp,6));
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot_dB(iddd,8)*(PProt_dA(idpp,2) - PProt_dA(idpp,3)) + DDrot_dB(iddd,9)*PProt_dA(idpp,6));
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot(iddd,8)*(PProt_dAB(idpp,2) - PProt_dAB(idpp,3)) + DDrot(iddd,9)*PProt_dAB(idpp,6));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14))*PProt(idpp,4) + (DDrot_dAB(iddd,12) - DDrot_dAB(iddd,13))*PProt(idpp,5));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*PProt_dB(idpp,4) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*PProt_dB(idpp,5));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot_dB(iddd,11) + DDrot_dB(iddd,14))*PProt_dA(idpp,4) + (DDrot_dB(iddd,12) - DDrot_dB(iddd,13))*PProt_dA(idpp,5));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt_dAB(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt_dAB(idpp,5));
          }
        }
      }
      if ((nrows > 10)&&(ncols > 10)) {
        //(sd|sd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(0,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[2] = d_eri2Center(0,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        //(sd|sd) rotate
        for (size_t idsdb = 1; idsdb < 6; ++idsdb) {
          for (size_t idsdk = 1; idsdk < 6; ++idsdk) {
            block(10 + idsdb,10 + idsdk) = intn[0]*SDrot_dAB(idsdb,1)*SDrot(idsdk,1) + intn[1]*(SDrot_dAB(idsdb,2)*SDrot(idsdk,2) + SDrot_dAB(idsdb,3)*SDrot(idsdk,3)) + intn[2]*(SDrot_dAB(idsdb,4)*SDrot(idsdk,4) + SDrot_dAB(idsdb,5)*SDrot(idsdk,5));
            block(10 + idsdb,10 + idsdk) += intn[0]*SDrot_dA(idsdb,1)*SDrot_dB(idsdk,1) + intn[1]*(SDrot_dA(idsdb,2)*SDrot_dB(idsdk,2) + SDrot_dA(idsdb,3)*SDrot_dB(idsdk,3)) + intn[2]*(SDrot_dA(idsdb,4)*SDrot_dB(idsdk,4) + SDrot_dA(idsdb,5)*SDrot_dB(idsdk,5));
            block(10 + idsdb,10 + idsdk) += intn[0]*SDrot_dB(idsdb,1)*SDrot_dA(idsdk,1) + intn[1]*(SDrot_dB(idsdb,2)*SDrot_dA(idsdk,2) + SDrot_dB(idsdb,3)*SDrot_dA(idsdk,3)) + intn[2]*(SDrot_dB(idsdb,4)*SDrot_dA(idsdk,4) + SDrot_dB(idsdb,5)*SDrot_dA(idsdk,5));
            block(10 + idsdb,10 + idsdk) += intn[0]*SDrot(idsdb,1)*SDrot_dAB(idsdk,1) + intn[1]*(SDrot(idsdb,2)*SDrot_dAB(idsdk,2) + SDrot(idsdb,3)*SDrot_dAB(idsdk,3)) + intn[2]*(SDrot(idsdb,4)*SDrot_dAB(idsdk,4) + SDrot(idsdb,5)*SDrot_dAB(idsdk,5));
          }
        }
        //(sd|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(0,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(0,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(sg,pi|-pi,-dl)
        intn[4] = d_eri2Center(0,0,2,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        //(sd|pd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(10 + idsd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot_dAB((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SDrot_dB(idsd,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SDrot_dB(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[0]*PDrot_dB((idpd1 - 1)*3 + idpd2,1)*SDrot_dA(idsd,1) + intn[1]*(PDrot_dB((idpd1 - 1)*3 + idpd2,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,3))*SDrot_dA(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot_dAB(idsd,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot_dAB(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SDrot_dB(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SDrot_dB(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot_dB((idpd1 - 1)*3 + idpd2,4)*SDrot_dA(idsd,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,5)*SDrot_dA(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot_dAB(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot_dAB(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot_dAB((idpd1 - 1)*3 + idpd2,10) + PDrot_dAB((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot_dAB((idpd1 - 1)*3 + idpd2,11) - PDrot_dAB((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SDrot_dB(idsd,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SDrot_dB(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot_dB((idpd1 - 1)*3 + idpd2,10) + PDrot_dB((idpd1 - 1)*3 + idpd2,15))*SDrot_dA(idsd,2) + (PDrot_dB((idpd1 - 1)*3 + idpd2,11) - PDrot_dB((idpd1 - 1)*3 + idpd2,14))*SDrot_dA(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot_dAB(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot_dAB(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SDrot_dB(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SDrot_dB(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot_dB((idpd1 - 1)*3 + idpd2,8)*SDrot_dA(idsd,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,12)*SDrot_dA(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot_dAB(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot_dAB(idsd,3));
            }
          }
        }
        //(pd|sd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(1,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(1,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center(1,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[4] = d_eri2Center(1,1,2,0,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        //(pd|sd) rotate
        for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
            for (size_t idsd = 1; idsd < 6; ++idsd) {
              block(5*(idpd2+2) + idpd1,10 + idsd) = -intn[0]*PDrot_dAB((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SDrot_dB(idsd,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SDrot_dB(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[0]*PDrot_dB((idpd1 - 1)*3 + idpd2,1)*SDrot_dA(idsd,1) + intn[1]*(PDrot_dB((idpd1 - 1)*3 + idpd2,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,3))*SDrot_dA(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot_dAB(idsd,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot_dAB(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SDrot_dB(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SDrot_dB(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot_dB((idpd1 - 1)*3 + idpd2,4)*SDrot_dA(idsd,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,5)*SDrot_dA(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot_dAB(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot_dAB(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot_dAB((idpd1 - 1)*3 + idpd2,10) + PDrot_dAB((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot_dAB((idpd1 - 1)*3 + idpd2,11) - PDrot_dAB((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SDrot_dB(idsd,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SDrot_dB(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot_dB((idpd1 - 1)*3 + idpd2,10) + PDrot_dB((idpd1 - 1)*3 + idpd2,15))*SDrot_dA(idsd,2) + (PDrot_dB((idpd1 - 1)*3 + idpd2,11) - PDrot_dB((idpd1 - 1)*3 + idpd2,14))*SDrot_dA(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot_dAB(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot_dAB(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot_dAB((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot_dAB((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SDrot_dB(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SDrot_dB(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot_dB((idpd1 - 1)*3 + idpd2,8)*SDrot_dA(idsd,2) + PDrot_dB((idpd1 - 1)*3 + idpd2,12)*SDrot_dA(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot_dAB(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot_dAB(idsd,3));
            }
          }
        }
        //(dd|sd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center(2,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,2,2,2,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center(2,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(2,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center(2,1,2,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        intn[6] = d_eri2Center(2,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)
        //(dd|sd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(30 + iddd,10 + idsd) = intn[0]*DDrot_dAB(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5))*SDrot_dA(idsd,1);
            block(30 + iddd,10 + idsd) += intn[0]*DDrot_dA(iddd,1)*SDrot_dB(idsd,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SDrot_dB(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot_dAB(idsd,1);
            block(30 + iddd,10 + idsd) += intn[0]*DDrot_dB(iddd,1)*SDrot_dA(idsd,1) + intn[1]*(DDrot_dB(iddd,2) + DDrot_dB(iddd,3))*SDrot_dA(idsd,1) + intn[2]*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5))*SDrot(idsd,1);
            block(30 + iddd,10 + idsd) += intn[0]*DDrot(iddd,1)*SDrot_dAB(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot_dAB(idsd,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SDrot_dB(idsd,1);
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot_dAB(iddd,6)*SDrot(idsd,2) + DDrot_dAB(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot_dAB(iddd,8)*SDrot(idsd,4) + DDrot_dAB(iddd,9)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot_dA(iddd,6)*SDrot_dB(idsd,2) + DDrot_dA(iddd,7)*SDrot_dB(idsd,3)) + intn[4]*(DDrot_dA(iddd,8)*SDrot_dB(idsd,4) + DDrot_dA(iddd,9)*SDrot_dB(idsd,5));
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot_dB(iddd,6)*SDrot_dA(idsd,2) + DDrot_dB(iddd,7)*SDrot_dA(idsd,3)) + intn[4]*(DDrot_dB(iddd,8)*SDrot_dA(idsd,4) + DDrot_dB(iddd,9)*SDrot_dA(idsd,5));
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot(iddd,6)*SDrot_dAB(idsd,2) + DDrot(iddd,7)*SDrot_dAB(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot_dAB(idsd,4) + DDrot(iddd,9)*SDrot_dAB(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot_dAB(iddd,2) - DDrot_dAB(iddd,3))*SDrot(idsd,4) + DDrot_dAB(iddd,10)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot_dA(iddd,2) - DDrot_dA(iddd,3))*SDrot_dB(idsd,4) + DDrot_dA(iddd,10)*SDrot_dB(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot_dB(iddd,2) - DDrot_dB(iddd,3))*SDrot_dA(idsd,4) + DDrot_dB(iddd,10)*SDrot_dA(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot_dAB(idsd,4) + DDrot(iddd,10)*SDrot_dAB(idsd,5));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14))*SDrot(idsd,2) + (DDrot_dB(iddd,12) - DDrot_dB(iddd,13))*SDrot_dA(idsd,3));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SDrot_dB(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot_dAB(idsd,3));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot_dB(iddd,11) + DDrot_dB(iddd,14))*SDrot_dA(idsd,2) + (DDrot_dAB(iddd,12) - DDrot_dAB(iddd,13))*SDrot(idsd,3));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot_dAB(idsd,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SDrot_dB(idsd,3));
          }
        }
        //(sd|dd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(0,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(0,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(0,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(0,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(0,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center(0,0,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        intn[6] = d_eri2Center(0,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        //(sd|dd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(10 + idsd,30 + iddd) = intn[0]*DDrot_dAB(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5))*SDrot_dA(idsd,1);
            block(10 + idsd,30 + iddd) += intn[0]*DDrot_dA(iddd,1)*SDrot_dB(idsd,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SDrot_dB(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot_dAB(idsd,1);
            block(10 + idsd,30 + iddd) += intn[0]*DDrot_dB(iddd,1)*SDrot_dA(idsd,1) + intn[1]*(DDrot_dB(iddd,2) + DDrot_dB(iddd,3))*SDrot_dA(idsd,1) + intn[2]*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5))*SDrot(idsd,1);
            block(10 + idsd,30 + iddd) += intn[0]*DDrot(iddd,1)*SDrot_dAB(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot_dAB(idsd,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SDrot_dB(idsd,1);
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot_dAB(iddd,6)*SDrot(idsd,2) + DDrot_dAB(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot_dAB(iddd,8)*SDrot(idsd,4) + DDrot_dAB(iddd,9)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot_dA(iddd,6)*SDrot_dB(idsd,2) + DDrot_dA(iddd,7)*SDrot_dB(idsd,3)) + intn[4]*(DDrot_dA(iddd,8)*SDrot_dB(idsd,4) + DDrot_dA(iddd,9)*SDrot_dB(idsd,5));
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot_dB(iddd,6)*SDrot_dA(idsd,2) + DDrot_dB(iddd,7)*SDrot_dA(idsd,3)) + intn[4]*(DDrot_dB(iddd,8)*SDrot_dA(idsd,4) + DDrot_dB(iddd,9)*SDrot_dA(idsd,5));
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot(iddd,6)*SDrot_dAB(idsd,2) + DDrot(iddd,7)*SDrot_dAB(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot_dAB(idsd,4) + DDrot(iddd,9)*SDrot_dAB(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot_dAB(iddd,2) - DDrot_dAB(iddd,3))*SDrot(idsd,4) + DDrot_dAB(iddd,10)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot_dA(iddd,2) - DDrot_dA(iddd,3))*SDrot_dB(idsd,4) + DDrot_dA(iddd,10)*SDrot_dB(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot_dB(iddd,2) - DDrot_dB(iddd,3))*SDrot_dA(idsd,4) + DDrot_dB(iddd,10)*SDrot_dA(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot_dAB(idsd,4) + DDrot(iddd,10)*SDrot_dAB(idsd,5));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14))*SDrot(idsd,2) + (DDrot_dB(iddd,12) - DDrot_dB(iddd,13))*SDrot_dA(idsd,3));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SDrot_dB(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot_dAB(idsd,3));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot_dB(iddd,11) + DDrot_dB(iddd,14))*SDrot_dA(idsd,2) + (DDrot_dAB(iddd,12) - DDrot_dAB(iddd,13))*SDrot(idsd,3));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot_dAB(idsd,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SDrot_dB(idsd,3));
          }
        }
        //(pd|pd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(1,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)=(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|-pi,-pi)
        intn[3] = d_eri2Center(1,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center(1,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[5] = d_eri2Center(1,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(pi,-dl|pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)
        intn[6] = d_eri2Center(1,1,2,0,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)=(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[7] = d_eri2Center(1,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[8] = d_eri2Center(1,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)=(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        //(pd|pd) rotate
        for (size_t idpdb1 = 1; idpdb1 < 6; ++idpdb1) {            //d orbital on bra
          for (size_t idpdb2 = 1; idpdb2 < 4; ++idpdb2) {          //p orbital on bra
            for (size_t idpdk1 = 1; idpdk1 < 6; ++idpdk1) {        //d orbital on ket
              for (size_t idpdk2 = 1; idpdk2 < 4; ++idpdk2) {      //p orbital on ket
                cnt1 = (idpdb1 - 1)*3 + idpdb2;
                cnt2 = (idpdk1 - 1)*3 + idpdk2;
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) = intn[0]*PDrot_dAB(cnt1,1)*PDrot(cnt2,1) + intn[1]*((PDrot_dAB(cnt1,2) + PDrot_dAB(cnt1,3))*PDrot(cnt2,1) + PDrot_dAB(cnt1,1)*(PDrot(cnt2,2) + PDrot(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[0]*PDrot_dA(cnt1,1)*PDrot_dB(cnt2,1) + intn[1]*((PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*PDrot_dB(cnt2,1) + PDrot_dA(cnt1,1)*(PDrot_dB(cnt2,2) + PDrot_dB(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[0]*PDrot_dB(cnt1,1)*PDrot_dA(cnt2,1) + intn[1]*((PDrot_dB(cnt1,2) + PDrot_dB(cnt1,3))*PDrot_dA(cnt2,1) + PDrot_dB(cnt1,1)*(PDrot_dA(cnt2,2) + PDrot_dA(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[0]*PDrot(cnt1,1)*PDrot_dAB(cnt2,1) + intn[1]*((PDrot(cnt1,2) + PDrot(cnt1,3))*PDrot_dAB(cnt2,1) + PDrot(cnt1,1)*(PDrot_dAB(cnt2,2) + PDrot_dAB(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot_dAB(cnt1,2) + PDrot_dAB(cnt1,3))*(PDrot(cnt2,2) + PDrot(cnt2,3)) + intn[3]*(PDrot_dAB(cnt1,4)*PDrot(cnt2,4) + PDrot_dAB(cnt1,5)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(PDrot_dB(cnt2,2) + PDrot_dB(cnt2,3)) + intn[3]*(PDrot_dA(cnt1,4)*PDrot_dB(cnt2,4) + PDrot_dA(cnt1,5)*PDrot_dB(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot_dB(cnt1,2) + PDrot_dB(cnt1,3))*(PDrot_dA(cnt2,2) + PDrot_dA(cnt2,3)) + intn[3]*(PDrot_dB(cnt1,4)*PDrot_dA(cnt2,4) + PDrot_dB(cnt1,5)*PDrot_dA(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(PDrot_dAB(cnt2,2) + PDrot_dAB(cnt2,3)) + intn[3]*(PDrot(cnt1,4)*PDrot_dAB(cnt2,4) + PDrot(cnt1,5)*PDrot_dAB(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot_dAB(cnt1,8)*PDrot(cnt2,8) + PDrot_dAB(cnt1,12)*PDrot(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot_dA(cnt1,8)*PDrot_dB(cnt2,8) + PDrot_dA(cnt1,12)*PDrot_dB(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot_dB(cnt1,8)*PDrot_dA(cnt2,8) + PDrot_dB(cnt1,12)*PDrot_dA(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot(cnt1,8)*PDrot_dAB(cnt2,8) + PDrot(cnt1,12)*PDrot_dAB(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot_dAB(cnt1,10) + PDrot_dAB(cnt1,15))*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot_dAB(cnt1,14) - PDrot_dAB(cnt1,11))*(PDrot(cnt2,14) - PDrot(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*(PDrot_dB(cnt2,10) + PDrot_dB(cnt2,15)) + (PDrot_dA(cnt1,14) - PDrot_dA(cnt1,11))*(PDrot_dB(cnt2,14) - PDrot_dB(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot_dB(cnt1,10) + PDrot_dB(cnt1,15))*(PDrot_dA(cnt2,10) + PDrot_dA(cnt2,15)) + (PDrot_dB(cnt1,14) - PDrot_dB(cnt1,11))*(PDrot_dA(cnt2,14) - PDrot_dA(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(PDrot_dAB(cnt2,10) + PDrot_dAB(cnt2,15)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(PDrot_dAB(cnt2,14) - PDrot_dAB(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot_dAB(cnt1,4)*PDrot(cnt2,8) + PDrot_dAB(cnt1,8)*PDrot(cnt2,4) + PDrot_dAB(cnt1,5)*PDrot(cnt2,12) + PDrot_dAB(cnt1,12)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot_dA(cnt1,4)*PDrot_dB(cnt2,8) + PDrot_dA(cnt1,8)*PDrot_dB(cnt2,4) + PDrot_dA(cnt1,5)*PDrot_dB(cnt2,12) + PDrot_dA(cnt1,12)*PDrot_dB(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot_dB(cnt1,4)*PDrot_dA(cnt2,8) + PDrot_dB(cnt1,8)*PDrot_dA(cnt2,4) + PDrot_dB(cnt1,5)*PDrot_dA(cnt2,12) + PDrot_dB(cnt1,12)*PDrot_dA(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot(cnt1,4)*PDrot_dAB(cnt2,8) + PDrot(cnt1,8)*PDrot_dAB(cnt2,4) + PDrot(cnt1,5)*PDrot_dAB(cnt2,12) + PDrot(cnt1,12)*PDrot_dAB(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot_dAB(cnt1,10) + PDrot_dAB(cnt1,15))*PDrot(cnt2,4) + PDrot_dAB(cnt1,4)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot_dAB(cnt1,11) - PDrot_dAB(cnt1,14))*PDrot(cnt2,5) + PDrot_dAB(cnt1,5)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*PDrot_dB(cnt2,4) + PDrot_dA(cnt1,4)*(PDrot_dB(cnt2,10) + PDrot_dB(cnt2,15)) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*PDrot_dB(cnt2,5) + PDrot_dA(cnt1,5)*(PDrot_dB(cnt2,11) - PDrot_dB(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot_dB(cnt1,10) + PDrot_dB(cnt1,15))*PDrot_dA(cnt2,4) + PDrot_dB(cnt1,4)*(PDrot_dA(cnt2,10) + PDrot_dA(cnt2,15)) + (PDrot_dB(cnt1,11) - PDrot_dB(cnt1,14))*PDrot_dA(cnt2,5) + PDrot_dB(cnt1,5)*(PDrot_dA(cnt2,11) - PDrot_dA(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot_dAB(cnt2,4) + PDrot(cnt1,4)*(PDrot_dAB(cnt2,10) + PDrot_dAB(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot_dAB(cnt2,5) + PDrot(cnt1,5)*(PDrot_dAB(cnt2,11) - PDrot_dAB(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot_dAB(cnt1,10) + PDrot_dAB(cnt1,15))*PDrot(cnt2,8) + PDrot_dAB(cnt1,8)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot_dAB(cnt1,11) - PDrot_dAB(cnt1,14))*PDrot(cnt2,12) + PDrot_dAB(cnt1,12)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*PDrot_dB(cnt2,8) + PDrot_dA(cnt1,8)*(PDrot_dB(cnt2,10) + PDrot_dB(cnt2,15)) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*PDrot_dB(cnt2,12) + PDrot_dA(cnt1,12)*(PDrot_dB(cnt2,11) - PDrot_dB(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot_dB(cnt1,10) + PDrot_dB(cnt1,15))*PDrot_dA(cnt2,8) + PDrot_dB(cnt1,8)*(PDrot_dA(cnt2,10) + PDrot_dA(cnt2,15)) + (PDrot_dB(cnt1,11) - PDrot_dB(cnt1,14))*PDrot_dA(cnt2,12) + PDrot_dB(cnt1,12)*(PDrot_dA(cnt2,11) - PDrot_dA(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot_dAB(cnt2,8) + PDrot(cnt1,8)*(PDrot_dAB(cnt2,10) + PDrot_dAB(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot_dAB(cnt2,12) + PDrot(cnt1,12)*(PDrot_dAB(cnt2,11) - PDrot_dAB(cnt2,14)));
              }
            }
          }
        }
        //(pd|dd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(1,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(1,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(1,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center(1,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center(1,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[5] = d_eri2Center(1,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|-dl,-dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)
        intn[6] = d_eri2Center(1,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(1,1,2,0,2,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center(1,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=(pi,-dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center(1,1,2,0,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)
        intn[10] = d_eri2Center(1,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        intn[11] = d_eri2Center(1,1,2,2,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(5*(idpd2+2) + idpd1,30 + iddd) = -intn[0]*PDrot_dAB(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot_dAB(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[0]*PDrot_dA(cnt1,1)*DDrot_dB(iddd,1) + intn[1]*PDrot_dA(cnt1,1)*(DDrot_dB(iddd,2) + DDrot_dB(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[0]*PDrot_dB(cnt1,1)*DDrot_dA(iddd,1) + intn[1]*PDrot_dB(cnt1,1)*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[0]*PDrot(cnt1,1)*DDrot_dAB(iddd,1) + intn[1]*PDrot(cnt1,1)*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot_dAB(cnt1,2) + PDrot_dAB(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot_dAB(cnt1,2) + PDrot_dAB(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*DDrot_dB(iddd,1) + intn[3]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot_dB(iddd,2) + DDrot_dB(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot_dB(cnt1,2) + PDrot_dB(cnt1,3))*DDrot_dA(iddd,1) + intn[3]*(PDrot_dB(cnt1,2) + PDrot_dB(cnt1,3))*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot_dAB(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot_dAB(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot_dAB(cnt1,2) + PDrot_dAB(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot_dA(cnt1,1)*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5)) + intn[5]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot_dB(cnt1,1)*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5)) + intn[5]*(PDrot_dB(cnt1,2) + PDrot_dB(cnt1,3))*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot(cnt1,1)*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot_dAB(cnt1,4)*DDrot(iddd,6) + PDrot_dAB(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot_dAB(cnt1,8)*DDrot(iddd,6) + PDrot_dAB(cnt1,12)*DDrot(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot_dA(cnt1,4)*DDrot_dB(iddd,6) + PDrot_dA(cnt1,5)*DDrot_dB(iddd,7)) + intn[7]*(PDrot_dA(cnt1,8)*DDrot_dB(iddd,6) + PDrot_dA(cnt1,12)*DDrot_dB(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot_dB(cnt1,4)*DDrot_dA(iddd,6) + PDrot_dB(cnt1,5)*DDrot_dA(iddd,7)) + intn[7]*(PDrot_dB(cnt1,8)*DDrot_dA(iddd,6) + PDrot_dB(cnt1,12)*DDrot_dA(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot(cnt1,4)*DDrot_dAB(iddd,6) + PDrot(cnt1,5)*DDrot_dAB(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot_dAB(iddd,6) + PDrot(cnt1,12)*DDrot_dAB(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot_dAB(cnt1,10) + PDrot_dAB(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot_dAB(cnt1,14) - PDrot_dAB(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*(DDrot_dB(iddd,11) + DDrot_dB(iddd,14)) + (PDrot_dA(cnt1,14) - PDrot_dA(cnt1,11))*(DDrot_dB(iddd,13) - DDrot_dB(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot_dB(cnt1,10) + PDrot_dB(cnt1,15))*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + (PDrot_dB(cnt1,14) - PDrot_dB(cnt1,11))*(DDrot_dA(iddd,13) - DDrot_dA(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot_dAB(iddd,13) - DDrot_dAB(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot_dAB(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot_dAB(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot_dA(cnt1,8)*(DDrot_dB(iddd,11) + DDrot_dB(iddd,14)) + PDrot_dA(cnt1,12)*(DDrot_dB(iddd,12) - DDrot_dB(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot_dB(cnt1,8)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot_dB(cnt1,12)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot(cnt1,8)*(DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14)) + PDrot(cnt1,12)*(DDrot_dAB(iddd,12) - DDrot_dAB(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot_dAB(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot_dB(cnt1,5)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot_dA(cnt1,4)*(DDrot_dB(iddd,11) + DDrot_dB(iddd,14)) + PDrot(cnt1,5)*(DDrot_dAB(iddd,12) - DDrot_dAB(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot_dB(cnt1,4)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot_dAB(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot(cnt1,4)*(DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14)) + PDrot_dA(cnt1,5)*(DDrot_dB(iddd,12) - DDrot_dB(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot_dAB(cnt1,10) + PDrot_dAB(cnt1,15))*DDrot(iddd,6) + (PDrot_dB(cnt1,11) - PDrot_dB(cnt1,14))*DDrot_dA(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*DDrot_dB(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot_dAB(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot_dB(cnt1,10) + PDrot_dB(cnt1,15))*DDrot_dA(iddd,6) + (PDrot_dAB(cnt1,11) - PDrot_dAB(cnt1,14))*DDrot(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot_dAB(iddd,6) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*DDrot_dB(iddd,7));
            }
          }
        }
        //(dd|pd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center(2,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center(2,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center(2,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center(2,2,2,2,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center(2,2,2,2,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|-pi,-pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)
        intn[6] = d_eri2Center(2,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center(2,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center(2,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center(2,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        intn[10] = d_eri2Center(2,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[11] = d_eri2Center(2,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(dd|pd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(30 + iddd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot_dAB(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot_dAB(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[0]*PDrot_dA(cnt1,1)*DDrot_dB(iddd,1) + intn[1]*PDrot_dA(cnt1,1)*(DDrot_dB(iddd,2) + DDrot_dB(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[0]*PDrot_dB(cnt1,1)*DDrot_dA(iddd,1) + intn[1]*PDrot_dB(cnt1,1)*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[0]*PDrot(cnt1,1)*DDrot_dAB(iddd,1) + intn[1]*PDrot(cnt1,1)*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot_dAB(cnt1,2) + PDrot_dAB(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot_dB(cnt1,2) + PDrot_dB(cnt1,3))*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*DDrot_dB(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dAB(iddd,2) + DDrot_dAB(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot_dB(cnt1,2) + PDrot_dB(cnt1,3))*DDrot_dA(iddd,1) + intn[3]*(PDrot_dAB(cnt1,2) + PDrot_dAB(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot_dAB(iddd,1) + intn[3]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot_dB(iddd,2) + DDrot_dB(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot_dAB(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot_dAB(cnt1,2) + PDrot_dAB(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot_dA(cnt1,1)*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5)) + intn[5]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot_dB(iddd,4) + DDrot_dB(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot_dB(cnt1,1)*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5)) + intn[5]*(PDrot_dB(cnt1,2) + PDrot_dB(cnt1,3))*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot(cnt1,1)*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dAB(iddd,4) + DDrot_dAB(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot_dAB(cnt1,4)*DDrot(iddd,6) + PDrot_dAB(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot_dAB(cnt1,8)*DDrot(iddd,6) + PDrot_dAB(cnt1,12)*DDrot(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot_dA(cnt1,4)*DDrot_dB(iddd,6) + PDrot_dA(cnt1,5)*DDrot_dB(iddd,7)) + intn[7]*(PDrot_dA(cnt1,8)*DDrot_dB(iddd,6) + PDrot_dA(cnt1,12)*DDrot_dB(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot_dB(cnt1,4)*DDrot_dA(iddd,6) + PDrot_dB(cnt1,5)*DDrot_dA(iddd,7)) + intn[7]*(PDrot_dB(cnt1,8)*DDrot_dA(iddd,6) + PDrot_dB(cnt1,12)*DDrot_dA(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot(cnt1,4)*DDrot_dAB(iddd,6) + PDrot(cnt1,5)*DDrot_dAB(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot_dAB(iddd,6) + PDrot(cnt1,12)*DDrot_dAB(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot_dAB(cnt1,10) + PDrot_dAB(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot_dAB(cnt1,14) - PDrot_dAB(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*(DDrot_dB(iddd,11) + DDrot_dB(iddd,14)) + (PDrot_dA(cnt1,14) - PDrot_dA(cnt1,11))*(DDrot_dB(iddd,13) - DDrot_dB(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot_dB(cnt1,10) + PDrot_dB(cnt1,15))*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + (PDrot_dB(cnt1,14) - PDrot_dB(cnt1,11))*(DDrot_dA(iddd,13) - DDrot_dA(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot_dAB(iddd,13) - DDrot_dAB(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot_dAB(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot_dB(cnt1,12)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot_dA(cnt1,8)*(DDrot_dB(iddd,11) + DDrot_dB(iddd,14)) + PDrot(cnt1,12)*(DDrot_dAB(iddd,12) - DDrot_dAB(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot_dB(cnt1,8)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot_dAB(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot(cnt1,8)*(DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14)) + PDrot_dA(cnt1,12)*(DDrot_dB(iddd,12) - DDrot_dB(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot_dAB(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot_dB(cnt1,5)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot_dA(cnt1,4)*(DDrot_dB(iddd,11) + DDrot_dB(iddd,14)) + PDrot(cnt1,5)*(DDrot_dAB(iddd,12) - DDrot_dAB(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot_dB(cnt1,4)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot_dAB(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot(cnt1,4)*(DDrot_dAB(iddd,11) + DDrot_dAB(iddd,14)) + PDrot_dA(cnt1,5)*(DDrot_dB(iddd,12) - DDrot_dB(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot_dAB(cnt1,10) + PDrot_dAB(cnt1,15))*DDrot(iddd,6) + (PDrot_dB(cnt1,11) - PDrot_dB(cnt1,14))*DDrot_dA(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*DDrot_dB(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot_dAB(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot_dB(cnt1,10) + PDrot_dB(cnt1,15))*DDrot_dA(iddd,6) + (PDrot_dAB(cnt1,11) - PDrot_dAB(cnt1,14))*DDrot(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot_dAB(iddd,6) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*DDrot_dB(iddd,7));
            }
          }
        }
        //(dd|dd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center(2,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center(2,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center(2,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center(2,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center(2,2,2,2,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center(2,2,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center(2,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|-dl,-dl)
        intn[7] = d_eri2Center(2,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[8] = d_eri2Center(2,2,2,2,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(dl,dl|dl,dl)=(-dl,-dl|-dl,-dl)
        intn[9] = d_eri2Center(2,1,2,1,2,-1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[10] = d_eri2Center(2,-2,2,-2,2,2,2,2,RCD,Dd,atmC,atmD,type); //(-dl,-dl|dl,dl)=(dl,dl|-dl,-dl)
        intn[11] = d_eri2Center(2,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[12] = d_eri2Center(2,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[13] = d_eri2Center(2,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)
        intn[14] = d_eri2Center(2,1,2,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type); //(pi,-pi|pi,-pi)
        intn[15] = d_eri2Center(2,1,2,-2,2,0,2,-1,RCD,Dd,atmC,atmD,type); //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)=(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        intn[16] = d_eri2Center(2,1,2,1,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)=(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(dd|dd) rotate
        for (size_t iddd1 = 1; iddd1 < 16; ++iddd1) {
          for (size_t iddd2 = 1; iddd2 < 16; ++iddd2) {
            block(30 + iddd1,30 + iddd2) = intn[0]*DDrot_dAB(iddd1,1)*DDrot(iddd2,1) + intn[1]*DDrot_dAB(iddd1,1)*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[2]*DDrot_dB(iddd1,1)*(DDrot_dA(iddd2,4) + DDrot_dA(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[0]*DDrot_dA(iddd1,1)*DDrot_dB(iddd2,1) + intn[1]*DDrot_dA(iddd1,1)*(DDrot_dB(iddd2,2) + DDrot_dB(iddd2,3)) + intn[2]*DDrot(iddd1,1)*(DDrot_dAB(iddd2,4) + DDrot_dAB(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[0]*DDrot_dB(iddd1,1)*DDrot_dA(iddd2,1) + intn[1]*DDrot_dB(iddd1,1)*(DDrot_dA(iddd2,2) + DDrot_dA(iddd2,3)) + intn[2]*DDrot_dAB(iddd1,1)*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[0]*DDrot(iddd1,1)*DDrot_dAB(iddd2,1) + intn[1]*DDrot(iddd1,1)*(DDrot_dAB(iddd2,2) + DDrot_dAB(iddd2,3)) + intn[2]*DDrot_dA(iddd1,1)*(DDrot_dB(iddd2,4) + DDrot_dB(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot_dAB(iddd1,2) + DDrot_dAB(iddd1,3))*DDrot(iddd2,1) + intn[4]*(DDrot_dB(iddd1,4) + DDrot_dB(iddd1,5))*DDrot_dA(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot_dA(iddd1,2) + DDrot_dA(iddd1,3))*DDrot_dB(iddd2,1) + intn[4]*(DDrot(iddd1,4) + DDrot(iddd1,5))*DDrot_dAB(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot_dB(iddd1,2) + DDrot_dB(iddd1,3))*DDrot_dA(iddd2,1) + intn[4]*(DDrot_dAB(iddd1,4) + DDrot_dAB(iddd1,5))*DDrot(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot(iddd1,2) + DDrot(iddd1,3))*DDrot_dAB(iddd2,1) + intn[4]*(DDrot_dA(iddd1,4) + DDrot_dA(iddd1,5))*DDrot_dB(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot_dAB(iddd1,4) + DDrot_dAB(iddd1,5))*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[6]*(DDrot_dAB(iddd1,2) + DDrot_dAB(iddd1,3))*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot_dA(iddd1,4) + DDrot_dA(iddd1,5))*(DDrot_dB(iddd2,2) + DDrot_dB(iddd2,3)) + intn[6]*(DDrot_dA(iddd1,2) + DDrot_dA(iddd1,3))*(DDrot_dB(iddd2,4) + DDrot_dB(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot_dB(iddd1,4) + DDrot_dB(iddd1,5))*(DDrot_dA(iddd2,2) + DDrot_dA(iddd2,3)) + intn[6]*(DDrot_dB(iddd1,2) + DDrot_dB(iddd1,3))*(DDrot_dA(iddd2,4) + DDrot_dA(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot(iddd1,4) + DDrot(iddd1,5))*(DDrot_dAB(iddd2,2) + DDrot_dAB(iddd2,3)) + intn[6]*(DDrot(iddd1,2) + DDrot(iddd1,3))*(DDrot_dAB(iddd2,4) + DDrot_dAB(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot_dAB(iddd1,2)*DDrot(iddd2,2) + DDrot_dAB(iddd1,3)*DDrot(iddd2,3)) + intn[8]*(DDrot_dAB(iddd1,4)*DDrot(iddd2,4) + DDrot_dAB(iddd1,5)*DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot_dA(iddd1,2)*DDrot_dB(iddd2,2) + DDrot_dA(iddd1,3)*DDrot_dB(iddd2,3)) + intn[8]*(DDrot_dA(iddd1,4)*DDrot_dB(iddd2,4) + DDrot_dA(iddd1,5)*DDrot_dB(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot_dB(iddd1,2)*DDrot_dA(iddd2,2) + DDrot_dB(iddd1,3)*DDrot_dA(iddd2,3)) + intn[8]*(DDrot_dB(iddd1,4)*DDrot_dA(iddd2,4) + DDrot_dB(iddd1,5)*DDrot_dA(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot(iddd1,2)*DDrot_dAB(iddd2,2) + DDrot(iddd1,3)*DDrot_dAB(iddd2,3)) + intn[8]*(DDrot(iddd1,4)*DDrot_dAB(iddd2,4) + DDrot(iddd1,5)*DDrot_dAB(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot_dAB(iddd1,2)*DDrot(iddd2,3) + DDrot_dAB(iddd1,3)*DDrot(iddd2,2)) + intn[10]*(DDrot_dAB(iddd1,4)*DDrot(iddd2,5) + DDrot_dAB(iddd1,5)*DDrot(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot_dA(iddd1,2)*DDrot_dB(iddd2,3) + DDrot_dA(iddd1,3)*DDrot_dB(iddd2,2)) + intn[10]*(DDrot_dA(iddd1,4)*DDrot_dB(iddd2,5) + DDrot_dA(iddd1,5)*DDrot_dB(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot_dB(iddd1,2)*DDrot_dA(iddd2,3) + DDrot_dB(iddd1,3)*DDrot_dA(iddd2,2)) + intn[10]*(DDrot_dB(iddd1,4)*DDrot_dA(iddd2,5) + DDrot_dB(iddd1,5)*DDrot_dA(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot(iddd1,2)*DDrot_dAB(iddd2,3) + DDrot(iddd1,3)*DDrot_dAB(iddd2,2)) + intn[10]*(DDrot(iddd1,4)*DDrot_dAB(iddd2,5) + DDrot(iddd1,5)*DDrot_dAB(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot_dAB(iddd1,6)*DDrot(iddd2,6) + DDrot_dAB(iddd1,7)*DDrot(iddd2,7)) + intn[12]*(DDrot_dAB(iddd1,8)*DDrot(iddd2,8) + DDrot_dAB(iddd1,9)*DDrot(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot_dA(iddd1,6)*DDrot_dB(iddd2,6) + DDrot_dA(iddd1,7)*DDrot_dB(iddd2,7)) + intn[12]*(DDrot_dA(iddd1,8)*DDrot_dB(iddd2,8) + DDrot_dA(iddd1,9)*DDrot_dB(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot_dB(iddd1,6)*DDrot_dA(iddd2,6) + DDrot_dB(iddd1,7)*DDrot_dA(iddd2,7)) + intn[12]*(DDrot_dB(iddd1,8)*DDrot_dA(iddd2,8) + DDrot_dB(iddd1,9)*DDrot_dA(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot(iddd1,6)*DDrot_dAB(iddd2,6) + DDrot(iddd1,7)*DDrot_dAB(iddd2,7)) + intn[12]*(DDrot(iddd1,8)*DDrot_dAB(iddd2,8) + DDrot(iddd1,9)*DDrot_dAB(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot_dAB(iddd1,11) + DDrot_dAB(iddd1,14))*(DDrot(iddd2,11) + DDrot(iddd2,14)) + (DDrot_dAB(iddd1,12) - DDrot_dAB(iddd1,13))*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot_dA(iddd1,11) + DDrot_dA(iddd1,14))*(DDrot_dB(iddd2,11) + DDrot_dB(iddd2,14)) + (DDrot_dA(iddd1,12) - DDrot_dA(iddd1,13))*(DDrot_dB(iddd2,12) - DDrot_dB(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot_dB(iddd1,11) + DDrot_dB(iddd1,14))*(DDrot_dA(iddd2,11) + DDrot_dA(iddd2,14)) + (DDrot_dB(iddd1,12) - DDrot_dB(iddd1,13))*(DDrot_dA(iddd2,12) - DDrot_dA(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot(iddd1,11) + DDrot(iddd1,14))*(DDrot_dAB(iddd2,11) + DDrot_dAB(iddd2,14)) + (DDrot(iddd1,12) - DDrot(iddd1,13))*(DDrot_dAB(iddd2,12) - DDrot_dAB(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[14]*DDrot_dAB(iddd1,10)*DDrot(iddd2,10) + intn[14]*DDrot_dB(iddd1,10)*DDrot_dA(iddd2,10);
            block(30 + iddd1,30 + iddd2) += intn[14]*DDrot_dA(iddd1,10)*DDrot_dB(iddd2,10) + intn[14]*DDrot(iddd1,10)*DDrot_dAB(iddd2,10);
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot_dAB(iddd1,11) + DDrot_dAB(iddd1,14))*DDrot(iddd2,6) + (DDrot_dB(iddd1,12) - DDrot_dB(iddd1,13))*DDrot_dA(iddd2,7) + DDrot_dAB(iddd1,6)*(DDrot(iddd2,11) + DDrot(iddd2,14)) + DDrot_dB(iddd1,7)*(DDrot_dA(iddd2,12) - DDrot_dA(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot_dA(iddd1,11) + DDrot_dA(iddd1,14))*DDrot_dB(iddd2,6) + (DDrot(iddd1,12) - DDrot(iddd1,13))*DDrot_dAB(iddd2,7) + DDrot_dA(iddd1,6)*(DDrot_dB(iddd2,11) + DDrot_dB(iddd2,14)) + DDrot(iddd1,7)*(DDrot_dAB(iddd2,12) - DDrot_dAB(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot_dB(iddd1,11) + DDrot_dB(iddd1,14))*DDrot_dA(iddd2,6) + (DDrot_dAB(iddd1,12) - DDrot_dAB(iddd1,13))*DDrot(iddd2,7) + DDrot_dB(iddd1,6)*(DDrot_dA(iddd2,11) + DDrot_dA(iddd2,14)) + DDrot_dAB(iddd1,7)*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot(iddd1,11) + DDrot(iddd1,14))*DDrot_dAB(iddd2,6) + (DDrot_dA(iddd1,12) - DDrot_dA(iddd1,13))*DDrot_dB(iddd2,7) + DDrot(iddd1,6)*(DDrot_dAB(iddd2,11) + DDrot_dAB(iddd2,14)) + DDrot_dA(iddd1,7)*(DDrot_dB(iddd2,12) - DDrot_dB(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot_dAB(iddd1,2) - DDrot_dAB(iddd1,3))*DDrot(iddd2,8) + DDrot_dAB(iddd1,10)*DDrot(iddd2,9) + DDrot_dAB(iddd1,8)*(DDrot(iddd2,2) - DDrot(iddd2,3)) + DDrot_dAB(iddd1,9)*DDrot(iddd2,10));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot_dA(iddd1,2) - DDrot_dA(iddd1,3))*DDrot_dB(iddd2,8) + DDrot_dA(iddd1,10)*DDrot_dB(iddd2,9) + DDrot_dA(iddd1,8)*(DDrot_dB(iddd2,2) - DDrot_dB(iddd2,3)) + DDrot_dA(iddd1,9)*DDrot_dB(iddd2,10));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot_dB(iddd1,2) - DDrot_dB(iddd1,3))*DDrot_dA(iddd2,8) + DDrot_dB(iddd1,10)*DDrot_dA(iddd2,9) + DDrot_dB(iddd1,8)*(DDrot_dA(iddd2,2) - DDrot_dA(iddd2,3)) + DDrot_dB(iddd1,9)*DDrot_dA(iddd2,10));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot(iddd1,2) - DDrot(iddd1,3))*DDrot_dAB(iddd2,8) + DDrot(iddd1,10)*DDrot_dAB(iddd2,9) + DDrot(iddd1,8)*(DDrot_dAB(iddd2,2) - DDrot_dAB(iddd2,3)) + DDrot(iddd1,9)*DDrot_dAB(iddd2,10));
          }
        }
      }
    }
  }
  void IntegralBlock2C_dRA(int type, matrixE & block, int der, int atmC, int atmD, double RCD, double cost, double sint, double cosp, double sinp) {
    //function calculating a block of second-derivatives of two-center integrals with respect to an angle and the internuclear distance
    //2 -> theta; 3 -> phi
    matrixE SProt_dA(1,1);
    matrixE PProt_dA(1,1);
    matrixE SDrot_dA(1,1);
    matrixE PDrot_dA(1,1);
    matrixE DDrot_dA(1,1);
    int cnt1;
    int cnt2;
    size_t nrows = block.rows();
    size_t ncols = block.cols();
    if ((nrows > 1)||(ncols > 1)) {
      SProt = SPtransf(cost,sint,cosp,sinp);
      PProt = PPtransf(cost,sint,cosp,sinp);
      if (der == 2) {
        SProt_dA = SPtransf_dt(cost,sint,cosp,sinp);
        PProt_dA = PPtransf_dt(cost,sint,cosp,sinp);
      }
      else if (der == 3) {
        SProt_dA = SPtransf_dp(cost,sint,cosp,sinp);
        PProt_dA = PPtransf_dp(cost,sint,cosp,sinp);
      }
      if ((nrows > 10)||(ncols > 10)) {
        SDrot = SDtransf(cost,sint,cosp,sinp);
        PDrot = PDtransf(cost,sint,cosp,sinp);
        DDrot = DDtransf(cost,sint,cosp,sinp);
        if (der == 2) {
          SDrot_dA = SDtransf_dt(cost,sint,cosp,sinp);
          PDrot_dA = PDtransf_dt(cost,sint,cosp,sinp);
          DDrot_dA = DDtransf_dt(cost,sint,cosp,sinp);
        }
        else if (der == 3) {
          SDrot_dA = SDtransf_dp(cost,sint,cosp,sinp);
          PDrot_dA = PDtransf_dp(cost,sint,cosp,sinp);
          DDrot_dA = DDtransf_dp(cost,sint,cosp,sinp);
        }
      }
    }
    D[0] = Dvalue(atmC,1);
    D[1] = Dvalue(atmC,2);
    D[2] = Dvalue(atmD,1);
    D[3] = Dvalue(atmD,2);
    //integral calculation
    block(1,1) = 0.0;                                                                                                                   //(ss|ss)
    if (nrows > 1) {
      intn[0] = eri2Center_dR(0,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center_dR(1,0,1,0,0,0,0,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center_dR(1,1,1,1,0,0,0,0,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(idx + 2,1) = -intn[0]*SProt_dA(idx + 1,1);}                                                         //(sp|ss) integrals
        block(5 + idx,1) = intn[1]*PProt_dA(idx + 1,1) + intn[2]*(PProt_dA(idx + 1,2) + PProt_dA(idx + 1,3));   //(pp|ss) integrals
      }
      if (nrows > 10) {
        //(sd|ss) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR(0,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(10 + idsd,1) = intn[0]*SDrot_dA(idsd,1);
        }
        //(pd|ss) integrals
        Dd[0] = Dvalue(atmC,5);
        intn[0] = d_eri2Center_dR(1,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center_dR(1,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(5*(idpd2+2) + idpd1,1) = -(intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3)));
          }
        }
        //(dd|ss) integrals
        Dd[0] = Dvalue(atmC,6);
        intn[0] = d_eri2Center_dR(2,0,2,0,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,1,2,1,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(2,2,2,2,0,0,0,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(30 + iddd,1) = intn[0]*DDrot_dA(iddd,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3)) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5));
        }
      }
    }
    if ((ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center_dR(0,0,0,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center_dR(0,0,0,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center_dR(0,0,0,0,1,1,1,1,RCD,D,atmC,atmD,type);
      for (size_t idx = 0; idx < 6; ++idx) {
        if (idx < 3) {block(1,idx + 2) = -intn[0]*SProt_dA(idx + 1,1);}                                                         //(ss|sp) integrals
        block(1,5 + idx) = intn[1]*PProt_dA(idx + 1,1) + intn[2]*(PProt_dA(idx + 1,2) + PProt_dA(idx + 1,3));   //(ss|pp) integrals
      }
      if (ncols > 10) {
        //(ss|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(0,0,0,0,0,0,2,0,RCD,Dd,atmC,atmD,type);
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          block(1,10 + idsd) = intn[0]*SDrot_dA(idsd,1);
        }
        //(ss|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(0,0,0,0,1,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center_dR(0,0,0,0,1,1,2,1,RCD,Dd,atmC,atmD,type);
        for (int idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (int idpd2 = 1; idpd2 < 4; ++idpd2) {
            block(1,5*(idpd2 + 2) + idpd1) = -intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1) - intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3));
          }
        }
        //(ss|dd) integrals
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(0,0,0,0,2,0,2,0,RCD,Dd,atmC,atmD,type);
        intn[1] = d_eri2Center_dR(0,0,0,0,2,1,2,1,RCD,Dd,atmC,atmD,type);
        intn[2] = d_eri2Center_dR(0,0,0,0,2,2,2,2,RCD,Dd,atmC,atmD,type);
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          block(1,30 + iddd) = intn[0]*DDrot_dA(iddd,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3)) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5));
        }
      }
    }
    if ((nrows > 1)&&(ncols > 1)&&(type == 0)) {
      intn[0] = eri2Center_dR(0,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[1] = eri2Center_dR(0,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      intn[2] = eri2Center_dR(0,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);
      intn[3] = eri2Center_dR(0,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);
      intn[4] = eri2Center_dR(0,0,1,1,1,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idpb = 1; idpb < 4; ++idpb) {
        for (size_t idpk = 1; idpk < 4; ++idpk) {           //(sp|sp)
          block(idpb + 1,idpk + 1) = intn[0]*SProt_dA(idpb,1)*SProt(idpk,1) + intn[1]*(SProt_dA(idpb,2)*SProt(idpk,2) + SProt_dA(idpb,3)*SProt(idpk,3));
          block(idpb + 1,idpk + 1) += intn[0]*SProt(idpb,1)*SProt_dA(idpk,1) + intn[1]*(SProt(idpb,2)*SProt_dA(idpk,2) + SProt(idpb,3)*SProt_dA(idpk,3));
        }
        for (size_t idc = 1; idc < 7; ++idc) {              //(sp|pp)
          block(idpb + 1,idc + 4) = -SProt_dA(idpb,1)*(intn[2]*PProt(idc,1) + intn[3]*(PProt(idc,2) + PProt(idc,3))) - intn[4]*(SProt_dA(idpb,2)*PProt(idc,4) + SProt_dA(idpb,3)*PProt(idc,5));
          block(idpb + 1,idc + 4) += -SProt(idpb,1)*(intn[2]*PProt_dA(idc,1) + intn[3]*(PProt_dA(idc,2) + PProt_dA(idc,3))) - intn[4]*(SProt(idpb,2)*PProt_dA(idc,4) + SProt(idpb,3)*PProt_dA(idc,5));
        }
      }
      intn[0] = eri2Center_dR(1,0,1,0,1,0,1,0,RCD,D,atmC,atmD,type);                          //(ps,ps|ps,ps)
      intn[1] = eri2Center_dR(1,1,1,1,1,0,1,0,RCD,D,atmC,atmD,type);                          //(pp,pp|ps,ps)
      intn[2] = eri2Center_dR(1,0,1,0,1,1,1,1,RCD,D,atmC,atmD,type);                          //(ps,ps|pp,pp)
      intn[3] = eri2Center_dR(1,1,1,1,1,1,1,1,RCD,D,atmC,atmD,type);                          //(pp,pp|pp,pp)
      intn[4] = eri2Center_dR(1,1,1,0,1,1,1,0,RCD,D,atmC,atmD,type);                          //(pp,ps|pp,ps)
      intn[5] = eri2Center_dR(1,1,1,1,1,-1,1,-1,RCD,D,atmC,atmD,type);                        //(pp,pp|pp*,pp*)
      intn[6] = eri2Center_dR(1,-1,1,1,1,-1,1,1,RCD,D,atmC,atmD,type);                        //(pp*,pp|pp*,pp)
      intn[7] = eri2Center_dR(1,0,1,0,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[8] = eri2Center_dR(1,1,1,1,0,0,1,0,RCD,D,atmC,atmD,type);
      intn[9] = eri2Center_dR(1,0,1,1,0,0,1,1,RCD,D,atmC,atmD,type);
      for (size_t idr = 1; idr < 7; ++idr) {
        for (size_t idc = 1; idc < 4; ++idc) {                //(pp|sp)
          block(idr + 4,idc + 1) = -SProt_dA(idc,1)*(intn[7]*PProt(idr,1) + intn[8]*(PProt(idr,2) + PProt(idr,3))) - intn[9]*(SProt_dA(idc,2)*PProt(idr,4) + SProt_dA(idc,3)*PProt(idr,5));
          block(idr + 4,idc + 1) += -SProt(idc,1)*(intn[7]*PProt_dA(idr,1) + intn[8]*(PProt_dA(idr,2) + PProt_dA(idr,3))) - intn[9]*(SProt(idc,2)*PProt_dA(idr,4) + SProt(idc,3)*PProt_dA(idr,5));
        }
        for (size_t idc = 1; idc < 7; ++idc) {                //(pp|pp)
          block(idr + 4,idc + 4) = intn[0]*PProt_dA(idr,1)*PProt(idc,1) + intn[1]*(PProt_dA(idr,2) + PProt_dA(idr,3))*PProt(idc,1) + intn[2]*PProt_dA(idr,1)*(PProt(idc,2) + PProt(idc,3)) + intn[3]*(PProt_dA(idr,2)*PProt(idc,2) + PProt_dA(idr,3)*PProt(idc,3)) + intn[4]*(PProt_dA(idr,4)*PProt(idc,4) + PProt_dA(idr,5)*PProt(idc,5)) + intn[5]*(PProt_dA(idr,2)*PProt(idc,3) + PProt_dA(idr,3)*PProt(idc,2)) + intn[6]*PProt_dA(idr,6)*PProt(idc,6);
          block(idr + 4,idc + 4) += intn[0]*PProt(idr,1)*PProt_dA(idc,1) + intn[1]*(PProt(idr,2) + PProt(idr,3))*PProt_dA(idc,1) + intn[2]*PProt(idr,1)*(PProt_dA(idc,2) + PProt_dA(idc,3)) + intn[3]*(PProt(idr,2)*PProt_dA(idc,2) + PProt(idr,3)*PProt_dA(idc,3)) + intn[4]*(PProt(idr,4)*PProt_dA(idc,4) + PProt(idr,5)*PProt_dA(idc,5)) + intn[5]*(PProt(idr,2)*PProt_dA(idc,3) + PProt(idr,3)*PProt_dA(idc,2)) + intn[6]*PProt(idr,6)*PProt_dA(idc,6);
        }
      }
      if (nrows > 10) {
        //(sd|sp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR(0,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sd|sp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idsp = 1; idsp < 4; ++idsp) {
            block(10 + idsd,1 + idsp) = -intn[0]*SDrot_dA(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot_dA(idsd,2)*SProt(idsp,2) + SDrot_dA(idsd,3)*SProt(idsp,3));
            block(10 + idsd,1 + idsp) += -intn[0]*SDrot(idsd,1)*SProt_dA(idsp,1) - intn[1]*(SDrot(idsd,2)*SProt_dA(idsp,2) + SDrot(idsd,3)*SProt_dA(idsp,3));
          }
        }
        //(sd|pp) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center_dR(0,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(0,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(0,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(sd|pp) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(10 + idsd,4 + idpp) = intn[0]*SDrot_dA(idsd,1)*PProt(idpp,1) + intn[1]*SDrot_dA(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(10 + idsd,4 + idpp) += intn[0]*SDrot(idsd,1)*PProt_dA(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot_dA(idsd,2)*PProt(idpp,4) + SDrot_dA(idsd,3)*PProt(idpp,5));
            block(10 + idsd,4 + idpp) += intn[2]*(SDrot(idsd,2)*PProt_dA(idpp,4) + SDrot(idsd,3)*PProt_dA(idpp,5));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot_dA(idsd,4)*PProt(idpp,2) - SDrot_dA(idsd,4)*PProt(idpp,3) + SDrot_dA(idsd,5)*PProt(idpp,6));
            block(10 + idsd,4 + idpp) += intn[3]*(SDrot(idsd,4)*PProt_dA(idpp,2) - SDrot(idsd,4)*PProt_dA(idpp,3) + SDrot(idsd,5)*PProt_dA(idpp,6));
          }
        }
        //(pd|sp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR(1,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(1,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(1,1,2,0,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        intn[4] = d_eri2Center_dR(1,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,1 + idsp) = intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt_dA(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt_dA(idsp,1);
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt_dA(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt_dA(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt_dA(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt_dA(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
              block(5*(idpd2+2) + idpd1,1 + idsp) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt_dA(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt_dA(idsp,3));
            }
          }
        }
        //(dd|sp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,1);
        intn[0] = d_eri2Center_dR(2,0,2,0,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,1,2,1,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(2,2,2,2,0,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center_dR(2,0,2,1,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR(2,1,2,2,0,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|sp) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(30 + iddd,1 + idsp) = -intn[0]*DDrot_dA(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SProt(idsp,1);
            block(30 + iddd,1 + idsp) -= intn[0]*DDrot(iddd,1)*SProt_dA(idsp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt_dA(idsp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt_dA(idsp,1);
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot_dA(iddd,6)*SProt(idsp,2) + DDrot_dA(iddd,7)*SProt(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[3]*(DDrot(iddd,6)*SProt_dA(idsp,2) + DDrot(iddd,7)*SProt_dA(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SProt(idsp,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SProt(idsp,3));
            block(30 + iddd,1 + idsp) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt_dA(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt_dA(idsp,3));
          }
        }
        //(pd|pp) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center_dR(1,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(1,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center_dR(1,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center_dR(1,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center_dR(1,1,2,0,1,1,1,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center_dR(1,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|pp) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(5*(idpd2+2) + idpd1,4 + idpp) = -intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) += -intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt_dA(idpp,1) - intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt_dA(idpp,1);
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt_dA(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt_dA(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt_dA(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt_dA(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
              block(5*(idpd2+2) + idpd1,4 + idpp) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt_dA(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt_dA(idpp,5));
            }
          }
        }
        //(dd|pp) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,3);
        intn[0] = d_eri2Center_dR(2,0,2,0,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,1,2,1,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(2,2,2,2,1,0,1,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center_dR(2,0,2,0,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[4] = d_eri2Center_dR(2,1,2,1,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center_dR(2,2,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center_dR(2,0,2,1,1,0,1,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR(2,-1,2,-1,1,1,1,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center_dR(2,1,2,-1,1,1,1,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center_dR(2,0,2,2,1,1,1,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=(sg,-dl|pi,-pi)=-(sg,dl|-pi,-pi)
        intn[10] = d_eri2Center_dR(2,1,2,2,1,0,1,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(dd|pp) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(30 + iddd,4 + idpp) = intn[0]*DDrot_dA(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*PProt(idpp,1);
            block(30 + iddd,4 + idpp) += intn[0]*DDrot(iddd,1)*PProt_dA(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt_dA(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt_dA(idpp,1);
            block(30 + iddd,4 + idpp) += intn[3]*DDrot_dA(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot_dA(iddd,2)*PProt(idpp,2) + DDrot_dA(iddd,3)*PProt(idpp,3));
            block(30 + iddd,4 + idpp) += intn[3]*DDrot(iddd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt_dA(idpp,2) + DDrot(iddd,3)*PProt_dA(idpp,3));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot_dA(iddd,6)*PProt(idpp,4) + DDrot_dA(iddd,7)*PProt(idpp,5));
            block(30 + iddd,4 + idpp) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt_dA(idpp,4) + DDrot(iddd,7)*PProt_dA(idpp,5));
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot_dA(iddd,2)*PProt(idpp,3) + DDrot_dA(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot_dA(iddd,10)*PProt(idpp,6);
            block(30 + iddd,4 + idpp) += intn[7]*(DDrot(iddd,2)*PProt_dA(idpp,3) + DDrot(iddd,3)*PProt_dA(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt_dA(idpp,6);
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot_dA(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot_dA(iddd,9)*PProt(idpp,6));
            block(30 + iddd,4 + idpp) += intn[9]*(DDrot(iddd,8)*(PProt_dA(idpp,2) - PProt_dA(idpp,3)) + DDrot(iddd,9)*PProt_dA(idpp,6));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*PProt(idpp,4) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*PProt(idpp,5));
            block(30 + iddd,4 + idpp) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt_dA(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt_dA(idpp,5));
          }
        }
      }
      if (ncols > 10) {
        //(sp|sd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(0,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        //(sp|sd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(1 + idsp,10 + idsd) = -intn[0]*SDrot_dA(idsd,1)*SProt(idsp,1) - intn[1]*(SDrot_dA(idsd,2)*SProt(idsp,2) + SDrot_dA(idsd,3)*SProt(idsp,3));
            block(1 + idsp,10 + idsd) -= intn[0]*SDrot(idsd,1)*SProt_dA(idsp,1) + intn[1]*(SDrot(idsd,2)*SProt_dA(idsp,2) + SDrot(idsd,3)*SProt_dA(idsp,3));
          }
        }
        //(pp|sd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(1,0,1,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,1,1,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(1,0,1,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(1,1,1,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        //(pp|sd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,10 + idsd) = intn[0]*SDrot_dA(idsd,1)*PProt(idpp,1) + intn[1]*SDrot_dA(idsd,1)*(PProt(idpp,2) + PProt(idpp,3));
            block(4 + idpp,10 + idsd) += intn[0]*SDrot(idsd,1)*PProt_dA(idpp,1) + intn[1]*SDrot(idsd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot_dA(idsd,2)*PProt(idpp,4) + SDrot_dA(idsd,3)*PProt(idpp,5));
            block(4 + idpp,10 + idsd) += intn[2]*(SDrot(idsd,2)*PProt_dA(idpp,4) + SDrot(idsd,3)*PProt_dA(idpp,5));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot_dA(idsd,4)*PProt(idpp,2) - SDrot_dA(idsd,4)*PProt(idpp,3) + SDrot_dA(idsd,5)*PProt(idpp,6));
            block(4 + idpp,10 + idsd) += intn[3]*(SDrot(idsd,4)*PProt_dA(idpp,2) - SDrot(idsd,4)*PProt_dA(idpp,3) + SDrot(idsd,5)*PProt_dA(idpp,6));
          }
        }
        //(sp|pd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(0,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(0,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(0,0,1,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[4] = d_eri2Center_dR(0,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|pd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(1 + idsp,5*(idpd2+2) + idpd1) = intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SProt(idsp,1) + intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SProt(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SProt_dA(idsp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SProt_dA(idsp,1);
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SProt(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SProt_dA(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SProt_dA(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SProt(idsp,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SProt_dA(idsp,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SProt_dA(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SProt(idsp,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SProt(idsp,3));
              block(1 + idsp,5*(idpd2+2) + idpd1) += intn[4]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SProt_dA(idsp,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SProt_dA(idsp,3));
            }
          }
        }
        //(sp|dd) integrals
        Dd[0] = Dvalue(atmC,1);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(0,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(0,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR(0,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR(0,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(sp|dd) rotate
        for (size_t idsp = 1; idsp < 4; ++idsp) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(1 + idsp,30 + iddd) = -intn[0]*DDrot_dA(iddd,1)*SProt(idsp,1) - intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SProt(idsp,1) - intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SProt(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[0]*DDrot(iddd,1)*SProt_dA(idsp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SProt_dA(idsp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SProt_dA(idsp,1);
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot_dA(iddd,6)*SProt(idsp,2) + DDrot_dA(iddd,7)*SProt(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[3]*(DDrot(iddd,6)*SProt_dA(idsp,2) + DDrot(iddd,7)*SProt_dA(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SProt(idsp,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SProt(idsp,3));
            block(1 + idsp,30 + iddd) -= intn[4]*((DDrot(iddd,11) + DDrot(iddd,14))*SProt_dA(idsp,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SProt_dA(idsp,3));
          }
        }
        //(pp|pd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(1,0,1,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,0,1,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(1,1,1,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center_dR(1,1,1,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[4] = d_eri2Center_dR(1,0,1,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[5] = d_eri2Center_dR(1,1,1,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[6] = d_eri2Center_dR(1,0,1,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|pd) rotate
        for (size_t idpp = 1; idpp < 7; ++idpp) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(4 + idpp,5*(idpd2+2) + idpd1) = -intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*PProt(idpp,1) - intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*PProt(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*PProt_dA(idpp,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*PProt_dA(idpp,1);
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[2]*PDrot((idpd1 - 1)*3 + idpd2,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*(PProt(idpp,2) + PProt(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[3]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*(PProt_dA(idpp,2) + PProt_dA(idpp,3));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*PProt(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,4)*PProt_dA(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,5)*PProt_dA(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*PProt(idpp,4) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[5]*(PDrot((idpd1 - 1)*3 + idpd2,8)*PProt_dA(idpp,4) + PDrot((idpd1 - 1)*3 + idpd2,12)*PProt_dA(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*PProt(idpp,4) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*PProt(idpp,5));
              block(4 + idpp,5*(idpd2+2) + idpd1) -= intn[6]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*PProt_dA(idpp,4) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*PProt_dA(idpp,5));
            }
          }
        }
        //(pp|dd) integrals
        Dd[0] = Dvalue(atmC,3);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(1,0,1,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,0,1,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(1,0,1,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR(1,1,1,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center_dR(1,1,1,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[5] = d_eri2Center_dR(1,1,1,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)=(-pi,-pi|-dl,-dl)
        intn[6] = d_eri2Center_dR(1,0,1,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR(1,-1,1,-1,2,1,2,1,RCD,Dd,atmC,atmD,type);  //(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)
        intn[8] = d_eri2Center_dR(1,1,1,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,-pi|pi,-pi)
        intn[9] = d_eri2Center_dR(1,1,1,1,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=(pi,-pi|sg,-dl)=-(-pi,-pi|sg,dl)
        intn[10] = d_eri2Center_dR(1,0,1,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(pp|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpp = 1; idpp < 7; ++idpp) {
            block(4 + idpp,30 + iddd) = intn[0]*DDrot_dA(iddd,1)*PProt(idpp,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*PProt(idpp,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*PProt(idpp,1);
            block(4 + idpp,30 + iddd) += intn[0]*DDrot(iddd,1)*PProt_dA(idpp,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*PProt_dA(idpp,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*PProt_dA(idpp,1);
            block(4 + idpp,30 + iddd) += intn[3]*DDrot_dA(iddd,1)*(PProt(idpp,2) + PProt(idpp,3)) + intn[4]*(DDrot_dA(iddd,2)*PProt(idpp,2) + DDrot_dA(iddd,3)*PProt(idpp,3));
            block(4 + idpp,30 + iddd) += intn[3]*DDrot(iddd,1)*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[4]*(DDrot(iddd,2)*PProt_dA(idpp,2) + DDrot(iddd,3)*PProt_dA(idpp,3));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*(PProt(idpp,2) + PProt(idpp,3)) + intn[6]*(DDrot_dA(iddd,6)*PProt(idpp,4) + DDrot_dA(iddd,7)*PProt(idpp,5));
            block(4 + idpp,30 + iddd) += intn[5]*(DDrot(iddd,4) + DDrot(iddd,5))*(PProt_dA(idpp,2) + PProt_dA(idpp,3)) + intn[6]*(DDrot(iddd,6)*PProt_dA(idpp,4) + DDrot(iddd,7)*PProt_dA(idpp,5));
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot_dA(iddd,2)*PProt(idpp,3) + DDrot_dA(iddd,3)*PProt(idpp,2)) + intn[8]*DDrot_dA(iddd,10)*PProt(idpp,6);
            block(4 + idpp,30 + iddd) += intn[7]*(DDrot(iddd,2)*PProt_dA(idpp,3) + DDrot(iddd,3)*PProt_dA(idpp,2)) + intn[8]*DDrot(iddd,10)*PProt_dA(idpp,6);
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot_dA(iddd,8)*(PProt(idpp,2) - PProt(idpp,3)) + DDrot_dA(iddd,9)*PProt(idpp,6));
            block(4 + idpp,30 + iddd) += intn[9]*(DDrot(iddd,8)*(PProt_dA(idpp,2) - PProt_dA(idpp,3)) + DDrot(iddd,9)*PProt_dA(idpp,6));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*PProt(idpp,4) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*PProt(idpp,5));
            block(4 + idpp,30 + iddd) += intn[10]*((DDrot(iddd,11) + DDrot(iddd,14))*PProt_dA(idpp,4) + (DDrot(iddd,12) - DDrot(iddd,13))*PProt_dA(idpp,5));
          }
        }
      }
      if ((nrows > 10)&&(ncols > 10)) {
        //(sd|sd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(0,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[2] = d_eri2Center_dR(0,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        //(sd|sd) rotate
        for (size_t idsdb = 1; idsdb < 6; ++idsdb) {
          for (size_t idsdk = 1; idsdk < 6; ++idsdk) {
            block(10 + idsdb,10 + idsdk) = intn[0]*SDrot_dA(idsdb,1)*SDrot(idsdk,1) + intn[1]*(SDrot_dA(idsdb,2)*SDrot(idsdk,2) + SDrot_dA(idsdb,3)*SDrot(idsdk,3)) + intn[2]*(SDrot_dA(idsdb,4)*SDrot(idsdk,4) + SDrot_dA(idsdb,5)*SDrot(idsdk,5));
            block(10 + idsdb,10 + idsdk) += intn[0]*SDrot(idsdb,1)*SDrot_dA(idsdk,1) + intn[1]*(SDrot(idsdb,2)*SDrot_dA(idsdk,2) + SDrot(idsdb,3)*SDrot_dA(idsdk,3)) + intn[2]*(SDrot(idsdb,4)*SDrot_dA(idsdk,4) + SDrot(idsdb,5)*SDrot_dA(idsdk,5));
          }
        }
        //(sd|pd) integrals
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(0,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(0,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(0,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(sg,pi|-pi,-dl)
        intn[4] = d_eri2Center_dR(0,0,2,1,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        //(sd|pd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              block(10 + idsd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot_dA(idsd,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot_dA(idsd,1);
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot_dA(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot_dA(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot_dA(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot_dA(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
              block(10 + idsd,5*(idpd2+2) + idpd1) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot_dA(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot_dA(idsd,3));
            }
          }
        }
        //(pd|sd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(1,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(1,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[3] = d_eri2Center_dR(1,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[4] = d_eri2Center_dR(1,1,2,0,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)
        //(pd|sd) rotate
        for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
          for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
            for (size_t idsd = 1; idsd < 6; ++idsd) {
              block(5*(idpd2+2) + idpd1,10 + idsd) = -intn[0]*PDrot_dA((idpd1 - 1)*3 + idpd2,1)*SDrot(idsd,1) - intn[1]*(PDrot_dA((idpd1 - 1)*3 + idpd2,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,3))*SDrot(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[0]*PDrot((idpd1 - 1)*3 + idpd2,1)*SDrot_dA(idsd,1) + intn[1]*(PDrot((idpd1 - 1)*3 + idpd2,2) + PDrot((idpd1 - 1)*3 + idpd2,3))*SDrot_dA(idsd,1);
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot_dA((idpd1 - 1)*3 + idpd2,4)*SDrot(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,5)*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[2]*(PDrot((idpd1 - 1)*3 + idpd2,4)*SDrot_dA(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,5)*SDrot_dA(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot_dA((idpd1 - 1)*3 + idpd2,10) + PDrot_dA((idpd1 - 1)*3 + idpd2,15))*SDrot(idsd,2) + (PDrot_dA((idpd1 - 1)*3 + idpd2,11) - PDrot_dA((idpd1 - 1)*3 + idpd2,14))*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[3]*((PDrot((idpd1 - 1)*3 + idpd2,10) + PDrot((idpd1 - 1)*3 + idpd2,15))*SDrot_dA(idsd,2) + (PDrot((idpd1 - 1)*3 + idpd2,11) - PDrot((idpd1 - 1)*3 + idpd2,14))*SDrot_dA(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot_dA((idpd1 - 1)*3 + idpd2,8)*SDrot(idsd,2) + PDrot_dA((idpd1 - 1)*3 + idpd2,12)*SDrot(idsd,3));
              block(5*(idpd2+2) + idpd1,10 + idsd) -= intn[4]*(PDrot((idpd1 - 1)*3 + idpd2,8)*SDrot_dA(idsd,2) + PDrot((idpd1 - 1)*3 + idpd2,12)*SDrot_dA(idsd,3));
            }
          }
        }
        //(dd|sd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,4);
        intn[0] = d_eri2Center_dR(2,0,2,0,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,1,2,1,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(2,2,2,2,0,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[3] = d_eri2Center_dR(2,0,2,1,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR(2,0,2,2,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center_dR(2,1,2,1,0,0,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)
        intn[6] = d_eri2Center_dR(2,1,2,2,0,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)
        //(dd|sd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idsd = 1; idsd < 6; ++idsd) {
            block(30 + iddd,10 + idsd) = intn[0]*DDrot_dA(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot_dA(idsd,1);
            block(30 + iddd,10 + idsd) += intn[0]*DDrot(iddd,1)*SDrot_dA(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot_dA(idsd,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SDrot(idsd,1);
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot_dA(iddd,6)*SDrot(idsd,2) + DDrot_dA(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot_dA(iddd,8)*SDrot(idsd,4) + DDrot_dA(iddd,9)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[3]*(DDrot(iddd,6)*SDrot_dA(idsd,2) + DDrot(iddd,7)*SDrot_dA(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot_dA(idsd,4) + DDrot(iddd,9)*SDrot_dA(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot_dA(iddd,2) - DDrot_dA(iddd,3))*SDrot(idsd,4) + DDrot_dA(iddd,10)*SDrot(idsd,5));
            block(30 + iddd,10 + idsd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot_dA(idsd,4) + DDrot(iddd,10)*SDrot_dA(idsd,5));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SDrot(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot_dA(idsd,3));
            block(30 + iddd,10 + idsd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot_dA(idsd,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SDrot(idsd,3));
          }
        }
        //(sd|dd) integrals
        Dd[0] = Dvalue(atmC,4);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(0,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(0,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(0,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR(0,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR(0,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);    //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[5] = d_eri2Center_dR(0,0,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        intn[6] = d_eri2Center_dR(0,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        //(sd|dd) rotate
        for (size_t idsd = 1; idsd < 6; ++idsd) {
          for (size_t iddd = 1; iddd < 16; ++iddd) {
            block(10 + idsd,30 + iddd) = intn[0]*DDrot_dA(iddd,1)*SDrot(idsd,1) + intn[1]*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3))*SDrot(idsd,1) + intn[2]*(DDrot(iddd,4) + DDrot(iddd,5))*SDrot_dA(idsd,1);
            block(10 + idsd,30 + iddd) += intn[0]*DDrot(iddd,1)*SDrot_dA(idsd,1) + intn[1]*(DDrot(iddd,2) + DDrot(iddd,3))*SDrot_dA(idsd,1) + intn[2]*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5))*SDrot(idsd,1);
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot_dA(iddd,6)*SDrot(idsd,2) + DDrot_dA(iddd,7)*SDrot(idsd,3)) + intn[4]*(DDrot_dA(iddd,8)*SDrot(idsd,4) + DDrot_dA(iddd,9)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[3]*(DDrot(iddd,6)*SDrot_dA(idsd,2) + DDrot(iddd,7)*SDrot_dA(idsd,3)) + intn[4]*(DDrot(iddd,8)*SDrot_dA(idsd,4) + DDrot(iddd,9)*SDrot_dA(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot_dA(iddd,2) - DDrot_dA(iddd,3))*SDrot(idsd,4) + DDrot_dA(iddd,10)*SDrot(idsd,5));
            block(10 + idsd,30 + iddd) += intn[5]*((DDrot(iddd,2) - DDrot(iddd,3))*SDrot_dA(idsd,4) + DDrot(iddd,10)*SDrot_dA(idsd,5));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot_dA(iddd,11) + DDrot_dA(iddd,14))*SDrot(idsd,2) + (DDrot(iddd,12) - DDrot(iddd,13))*SDrot_dA(idsd,3));
            block(10 + idsd,30 + iddd) += intn[6]*((DDrot(iddd,11) + DDrot(iddd,14))*SDrot_dA(idsd,2) + (DDrot_dA(iddd,12) - DDrot_dA(iddd,13))*SDrot(idsd,3));
          }
        }
        //(pd|pd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(1,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)=(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(1,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|pi,pi)=(pi,pi|-pi,-pi)=(-pi,-pi|-pi,-pi)
        intn[3] = d_eri2Center_dR(1,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[4] = d_eri2Center_dR(1,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[5] = d_eri2Center_dR(1,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(pi,-dl|pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)
        intn[6] = d_eri2Center_dR(1,1,2,0,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(pi,sg|sg,pi)=(-pi,sg|sg,-pi)=(sg,pi|pi,sg)=(sg,-pi|-pi,sg)
        intn[7] = d_eri2Center_dR(1,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)=(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[8] = d_eri2Center_dR(1,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)=(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        //(pd|pd) rotate
        for (size_t idpdb1 = 1; idpdb1 < 6; ++idpdb1) {            //d orbital on bra
          for (size_t idpdb2 = 1; idpdb2 < 4; ++idpdb2) {          //p orbital on bra
            for (size_t idpdk1 = 1; idpdk1 < 6; ++idpdk1) {        //d orbital on ket
              for (size_t idpdk2 = 1; idpdk2 < 4; ++idpdk2) {      //p orbital on ket
                cnt1 = (idpdb1 - 1)*3 + idpdb2;
                cnt2 = (idpdk1 - 1)*3 + idpdk2;
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) = intn[0]*PDrot_dA(cnt1,1)*PDrot(cnt2,1) + intn[1]*((PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*PDrot(cnt2,1) + PDrot_dA(cnt1,1)*(PDrot(cnt2,2) + PDrot(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[0]*PDrot(cnt1,1)*PDrot_dA(cnt2,1) + intn[1]*((PDrot(cnt1,2) + PDrot(cnt1,3))*PDrot_dA(cnt2,1) + PDrot(cnt1,1)*(PDrot_dA(cnt2,2) + PDrot_dA(cnt2,3)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(PDrot(cnt2,2) + PDrot(cnt2,3)) + intn[3]*(PDrot_dA(cnt1,4)*PDrot(cnt2,4) + PDrot_dA(cnt1,5)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(PDrot_dA(cnt2,2) + PDrot_dA(cnt2,3)) + intn[3]*(PDrot(cnt1,4)*PDrot_dA(cnt2,4) + PDrot(cnt1,5)*PDrot_dA(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot_dA(cnt1,8)*PDrot(cnt2,8) + PDrot_dA(cnt1,12)*PDrot(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[4]*(PDrot(cnt1,8)*PDrot_dA(cnt2,8) + PDrot(cnt1,12)*PDrot_dA(cnt2,12));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot_dA(cnt1,14) - PDrot_dA(cnt1,11))*(PDrot(cnt2,14) - PDrot(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[5]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(PDrot_dA(cnt2,10) + PDrot_dA(cnt2,15)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(PDrot_dA(cnt2,14) - PDrot_dA(cnt2,11)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot_dA(cnt1,4)*PDrot(cnt2,8) + PDrot_dA(cnt1,8)*PDrot(cnt2,4) + PDrot_dA(cnt1,5)*PDrot(cnt2,12) + PDrot_dA(cnt1,12)*PDrot(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[6]*(PDrot(cnt1,4)*PDrot_dA(cnt2,8) + PDrot(cnt1,8)*PDrot_dA(cnt2,4) + PDrot(cnt1,5)*PDrot_dA(cnt2,12) + PDrot(cnt1,12)*PDrot_dA(cnt2,5));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*PDrot(cnt2,4) + PDrot_dA(cnt1,4)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*PDrot(cnt2,5) + PDrot_dA(cnt1,5)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[7]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot_dA(cnt2,4) + PDrot(cnt1,4)*(PDrot_dA(cnt2,10) + PDrot_dA(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot_dA(cnt2,5) + PDrot(cnt1,5)*(PDrot_dA(cnt2,11) - PDrot_dA(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*PDrot(cnt2,8) + PDrot_dA(cnt1,8)*(PDrot(cnt2,10) + PDrot(cnt2,15)) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*PDrot(cnt2,12) + PDrot_dA(cnt1,12)*(PDrot(cnt2,11) - PDrot(cnt2,14)));
                block(5*(idpdb2+2) + idpdb1,5*(idpdk2+2) + idpdk1) += intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*PDrot_dA(cnt2,8) + PDrot(cnt1,8)*(PDrot_dA(cnt2,10) + PDrot_dA(cnt2,15)) + (PDrot(cnt1,11) - PDrot(cnt1,14))*PDrot_dA(cnt2,12) + PDrot(cnt1,12)*(PDrot_dA(cnt2,11) - PDrot_dA(cnt2,14)));
              }
            }
          }
        }
        //(pd|dd) integrals
        Dd[0] = Dvalue(atmC,5);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(1,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(1,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(1,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[3] = d_eri2Center_dR(1,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center_dR(1,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[5] = d_eri2Center_dR(1,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|-dl,-dl)=(pi,pi|-dl,-dl)=(-pi,-pi|dl,dl)
        intn[6] = d_eri2Center_dR(1,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR(1,1,2,0,2,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center_dR(1,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,dl|-pi,dl)=-(-pi,dl|pi,-dl)=(pi,-dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center_dR(1,1,2,0,2,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,dl)=(pi,sg|-pi,-dl)=-(-pi,sg|-pi,dl)=(-pi,sg|pi,-dl)
        intn[10] = d_eri2Center_dR(1,0,2,1,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        intn[11] = d_eri2Center_dR(1,1,2,2,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)
        //(pd|dd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(5*(idpd2+2) + idpd1,30 + iddd) = -intn[0]*PDrot_dA(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot_dA(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[0]*PDrot(cnt1,1)*DDrot_dA(iddd,1) + intn[1]*PDrot(cnt1,1)*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot_dA(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot_dA(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[4]*PDrot(cnt1,1)*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot_dA(cnt1,4)*DDrot(iddd,6) + PDrot_dA(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot_dA(cnt1,8)*DDrot(iddd,6) + PDrot_dA(cnt1,12)*DDrot(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[6]*(PDrot(cnt1,4)*DDrot_dA(iddd,6) + PDrot(cnt1,5)*DDrot_dA(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot_dA(iddd,6) + PDrot(cnt1,12)*DDrot_dA(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot_dA(cnt1,14) - PDrot_dA(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot_dA(iddd,13) - DDrot_dA(iddd,12)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot_dA(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot_dA(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[9]*(PDrot(cnt1,8)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot(cnt1,12)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot_dA(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,5)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[10]*(PDrot(cnt1,4)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot_dA(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*DDrot(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot_dA(iddd,7));
              block(5*(idpd2+2) + idpd1,30 + iddd) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot_dA(iddd,6) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*DDrot(iddd,7));
            }
          }
        }
        //(dd|pd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,5);
        intn[0] = d_eri2Center_dR(2,0,2,0,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,1,2,1,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[2] = d_eri2Center_dR(2,0,2,0,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[3] = d_eri2Center_dR(2,1,2,1,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)=(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[4] = d_eri2Center_dR(2,2,2,2,1,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center_dR(2,2,2,2,1,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|-pi,-pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)
        intn[6] = d_eri2Center_dR(2,0,2,1,1,0,2,1,RCD,Dd,atmC,atmD,type);    //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[7] = d_eri2Center_dR(2,1,2,0,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,sg|pi,sg)=(-pi,sg|-pi,sg)
        intn[8] = d_eri2Center_dR(2,1,2,2,1,1,2,2,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,dl)=(-pi,-dl|-pi,-dl)=(-pi,-dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)
        intn[9] = d_eri2Center_dR(2,1,2,2,1,1,2,0,RCD,Dd,atmC,atmD,type);    //(pi,dl|pi,sg)=-(-pi,dl|-pi,sg)=(pi,-dl|-pi,sg)=(-pi,-dl|pi,sg)
        intn[10] = d_eri2Center_dR(2,1,2,2,1,0,2,1,RCD,Dd,atmC,atmD,type);   //(pi,dl|sg,pi)=-(-pi,dl|sg,-pi)=(pi,-dl|sg,-pi)=(-pi,-dl|sg,pi)
        intn[11] = d_eri2Center_dR(2,0,2,1,1,1,2,2,RCD,Dd,atmC,atmD,type);   //(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=-(sg,-pi|-pi,dl)=(sg,-pi|pi,-dl)
        //(dd|pd) rotate
        for (size_t iddd = 1; iddd < 16; ++iddd) {
          for (size_t idpd1 = 1; idpd1 < 6; ++idpd1) {
            for (size_t idpd2 = 1; idpd2 < 4; ++idpd2) {
              cnt1 = (idpd1 - 1)*3 + idpd2;
              block(30 + iddd,5*(idpd2+2) + idpd1) = -intn[0]*PDrot_dA(cnt1,1)*DDrot(iddd,1) - intn[1]*PDrot_dA(cnt1,1)*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[0]*PDrot(cnt1,1)*DDrot_dA(iddd,1) + intn[1]*PDrot(cnt1,1)*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*DDrot(iddd,1) + intn[3]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dA(iddd,2) + DDrot_dA(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[2]*(PDrot(cnt1,2) + PDrot(cnt1,3))*DDrot_dA(iddd,1) + intn[3]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot(iddd,2) + DDrot(iddd,3));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot_dA(cnt1,1)*(DDrot(iddd,4) + DDrot(iddd,5)) + intn[5]*(PDrot_dA(cnt1,2) + PDrot_dA(cnt1,3))*(DDrot(iddd,4) + DDrot(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[4]*PDrot(cnt1,1)*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5)) + intn[5]*(PDrot(cnt1,2) + PDrot(cnt1,3))*(DDrot_dA(iddd,4) + DDrot_dA(iddd,5));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot_dA(cnt1,4)*DDrot(iddd,6) + PDrot_dA(cnt1,5)*DDrot(iddd,7)) + intn[7]*(PDrot_dA(cnt1,8)*DDrot(iddd,6) + PDrot_dA(cnt1,12)*DDrot(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[6]*(PDrot(cnt1,4)*DDrot_dA(iddd,6) + PDrot(cnt1,5)*DDrot_dA(iddd,7)) + intn[7]*(PDrot(cnt1,8)*DDrot_dA(iddd,6) + PDrot(cnt1,12)*DDrot_dA(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*(DDrot(iddd,11) + DDrot(iddd,14)) + (PDrot_dA(cnt1,14) - PDrot_dA(cnt1,11))*(DDrot(iddd,13) - DDrot(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[8]*((PDrot(cnt1,10) + PDrot(cnt1,15))*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + (PDrot(cnt1,14) - PDrot(cnt1,11))*(DDrot_dA(iddd,13) - DDrot_dA(iddd,12)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot_dA(cnt1,8)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,12)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[9]*(PDrot(cnt1,8)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot_dA(cnt1,12)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot_dA(cnt1,4)*(DDrot(iddd,11) + DDrot(iddd,14)) + PDrot(cnt1,5)*(DDrot_dA(iddd,12) - DDrot_dA(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[10]*(PDrot(cnt1,4)*(DDrot_dA(iddd,11) + DDrot_dA(iddd,14)) + PDrot_dA(cnt1,5)*(DDrot(iddd,12) - DDrot(iddd,13)));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot_dA(cnt1,10) + PDrot_dA(cnt1,15))*DDrot(iddd,6) + (PDrot(cnt1,11) - PDrot(cnt1,14))*DDrot_dA(iddd,7));
              block(30 + iddd,5*(idpd2+2) + idpd1) -= intn[11]*((PDrot(cnt1,10) + PDrot(cnt1,15))*DDrot_dA(iddd,6) + (PDrot_dA(cnt1,11) - PDrot_dA(cnt1,14))*DDrot(iddd,7));
            }
          }
        }
        //(dd|dd) integrals
        Dd[0] = Dvalue(atmC,6);
        Dd[1] = Dvalue(atmD,6);
        intn[0] = d_eri2Center_dR(2,0,2,0,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(sg,sg|sg,sg)
        intn[1] = d_eri2Center_dR(2,0,2,0,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(sg,sg|pi,pi)=(sg,sg|-pi,-pi)
        intn[2] = d_eri2Center_dR(2,0,2,0,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(sg,sg|dl,dl)=(sg,sg|-dl,-dl)
        intn[3] = d_eri2Center_dR(2,1,2,1,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(pi,pi|sg,sg)=(-pi,-pi|sg,sg)
        intn[4] = d_eri2Center_dR(2,2,2,2,2,0,2,0,RCD,Dd,atmC,atmD,type);    //(dl,dl|sg,sg)=(-dl,-dl|sg,sg)
        intn[5] = d_eri2Center_dR(2,2,2,2,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(dl,dl|pi,pi)=(-dl,-dl|pi,pi)=(dl,dl|-pi,-pi)=(-dl,-dl|-pi,-pi)
        intn[6] = d_eri2Center_dR(2,1,2,1,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(pi,pi|dl,dl)=(-pi,-pi|dl,dl)=(pi,pi|-dl,-dl)=(-pi,-pi|-dl,-dl)
        intn[7] = d_eri2Center_dR(2,1,2,1,2,1,2,1,RCD,Dd,atmC,atmD,type);    //(pi,pi|pi,pi)=(-pi,-pi|-pi,-pi)
        intn[8] = d_eri2Center_dR(2,2,2,2,2,2,2,2,RCD,Dd,atmC,atmD,type);    //(dl,dl|dl,dl)=(-dl,-dl|-dl,-dl)
        intn[9] = d_eri2Center_dR(2,1,2,1,2,-1,2,-1,RCD,Dd,atmC,atmD,type);  //(pi,pi|-pi,-pi)=(-pi,-pi|pi,pi)
        intn[10] = d_eri2Center_dR(2,-2,2,-2,2,2,2,2,RCD,Dd,atmC,atmD,type); //(-dl,-dl|dl,dl)=(dl,dl|-dl,-dl)
        intn[11] = d_eri2Center_dR(2,0,2,1,2,0,2,1,RCD,Dd,atmC,atmD,type);   //(sg,pi|sg,pi)=(sg,-pi|sg,-pi)
        intn[12] = d_eri2Center_dR(2,0,2,2,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(sg,dl|sg,dl)=(sg,-dl|sg,-dl)
        intn[13] = d_eri2Center_dR(2,1,2,2,2,1,2,2,RCD,Dd,atmC,atmD,type);   //(pi,dl|pi,dl)=(pi,dl|-pi,-dl)=(-pi,-dl|pi,dl)=(-pi,-dl|-pi,-dl)=-(-pi,dl|pi,-dl)=-(pi,-dl|-pi,dl)=(-pi,dl|-pi,dl)=(pi,-dl|pi,-dl)
        intn[14] = d_eri2Center_dR(2,1,2,-1,2,1,2,-1,RCD,Dd,atmC,atmD,type); //(pi,-pi|pi,-pi)
        intn[15] = d_eri2Center_dR(2,1,2,-2,2,0,2,-1,RCD,Dd,atmC,atmD,type); //(pi,dl|sg,pi)=(-pi,-dl|sg,pi)=(pi,-dl|sg,-pi)=-(-pi,dl|sg,-pi)=(sg,pi|pi,dl)=(sg,pi|-pi,-dl)=(sg,-pi|pi,-dl)=-(sg,-pi|-pi,dl)
        intn[16] = d_eri2Center_dR(2,1,2,1,2,0,2,2,RCD,Dd,atmC,atmD,type);   //(pi,pi|sg,dl)=-(-pi,-pi|sg,dl)=(pi,-pi|sg,-dl)=(sg,dl|pi,pi)=-(sg,dl|-pi,-pi)=(sg,-dl|pi,-pi)
        //(dd|dd) rotate
        for (size_t iddd1 = 1; iddd1 < 16; ++iddd1) {
          for (size_t iddd2 = 1; iddd2 < 16; ++iddd2) {
            block(30 + iddd1,30 + iddd2) = intn[0]*DDrot_dA(iddd1,1)*DDrot(iddd2,1) + intn[1]*DDrot_dA(iddd1,1)*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[2]*DDrot(iddd1,1)*(DDrot_dA(iddd2,4) + DDrot_dA(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[0]*DDrot(iddd1,1)*DDrot_dA(iddd2,1) + intn[1]*DDrot(iddd1,1)*(DDrot_dA(iddd2,2) + DDrot_dA(iddd2,3)) + intn[2]*DDrot_dA(iddd1,1)*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot_dA(iddd1,2) + DDrot_dA(iddd1,3))*DDrot(iddd2,1) + intn[4]*(DDrot(iddd1,4) + DDrot(iddd1,5))*DDrot_dA(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[3]*(DDrot(iddd1,2) + DDrot(iddd1,3))*DDrot_dA(iddd2,1) + intn[4]*(DDrot_dA(iddd1,4) + DDrot_dA(iddd1,5))*DDrot(iddd2,1);
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot_dA(iddd1,4) + DDrot_dA(iddd1,5))*(DDrot(iddd2,2) + DDrot(iddd2,3)) + intn[6]*(DDrot_dA(iddd1,2) + DDrot_dA(iddd1,3))*(DDrot(iddd2,4) + DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[5]*(DDrot(iddd1,4) + DDrot(iddd1,5))*(DDrot_dA(iddd2,2) + DDrot_dA(iddd2,3)) + intn[6]*(DDrot(iddd1,2) + DDrot(iddd1,3))*(DDrot_dA(iddd2,4) + DDrot_dA(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot_dA(iddd1,2)*DDrot(iddd2,2) + DDrot_dA(iddd1,3)*DDrot(iddd2,3)) + intn[8]*(DDrot_dA(iddd1,4)*DDrot(iddd2,4) + DDrot_dA(iddd1,5)*DDrot(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[7]*(DDrot(iddd1,2)*DDrot_dA(iddd2,2) + DDrot(iddd1,3)*DDrot_dA(iddd2,3)) + intn[8]*(DDrot(iddd1,4)*DDrot_dA(iddd2,4) + DDrot(iddd1,5)*DDrot_dA(iddd2,5));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot_dA(iddd1,2)*DDrot(iddd2,3) + DDrot_dA(iddd1,3)*DDrot(iddd2,2)) + intn[10]*(DDrot_dA(iddd1,4)*DDrot(iddd2,5) + DDrot_dA(iddd1,5)*DDrot(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[9]*(DDrot(iddd1,2)*DDrot_dA(iddd2,3) + DDrot(iddd1,3)*DDrot_dA(iddd2,2)) + intn[10]*(DDrot(iddd1,4)*DDrot_dA(iddd2,5) + DDrot(iddd1,5)*DDrot_dA(iddd2,4));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot_dA(iddd1,6)*DDrot(iddd2,6) + DDrot_dA(iddd1,7)*DDrot(iddd2,7)) + intn[12]*(DDrot_dA(iddd1,8)*DDrot(iddd2,8) + DDrot_dA(iddd1,9)*DDrot(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[11]*(DDrot(iddd1,6)*DDrot_dA(iddd2,6) + DDrot(iddd1,7)*DDrot_dA(iddd2,7)) + intn[12]*(DDrot(iddd1,8)*DDrot_dA(iddd2,8) + DDrot(iddd1,9)*DDrot_dA(iddd2,9));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot_dA(iddd1,11) + DDrot_dA(iddd1,14))*(DDrot(iddd2,11) + DDrot(iddd2,14)) + (DDrot_dA(iddd1,12) - DDrot_dA(iddd1,13))*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[13]*((DDrot(iddd1,11) + DDrot(iddd1,14))*(DDrot_dA(iddd2,11) + DDrot_dA(iddd2,14)) + (DDrot(iddd1,12) - DDrot(iddd1,13))*(DDrot_dA(iddd2,12) - DDrot_dA(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[14]*DDrot_dA(iddd1,10)*DDrot(iddd2,10) + intn[14]*DDrot(iddd1,10)*DDrot_dA(iddd2,10);
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot_dA(iddd1,11) + DDrot_dA(iddd1,14))*DDrot(iddd2,6) + (DDrot(iddd1,12) - DDrot(iddd1,13))*DDrot_dA(iddd2,7) + DDrot_dA(iddd1,6)*(DDrot(iddd2,11) + DDrot(iddd2,14)) + DDrot(iddd1,7)*(DDrot_dA(iddd2,12) - DDrot_dA(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[15]*((DDrot(iddd1,11) + DDrot(iddd1,14))*DDrot_dA(iddd2,6) + (DDrot_dA(iddd1,12) - DDrot_dA(iddd1,13))*DDrot(iddd2,7) + DDrot(iddd1,6)*(DDrot_dA(iddd2,11) + DDrot_dA(iddd2,14)) + DDrot_dA(iddd1,7)*(DDrot(iddd2,12) - DDrot(iddd2,13)));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot_dA(iddd1,2) - DDrot_dA(iddd1,3))*DDrot(iddd2,8) + DDrot_dA(iddd1,10)*DDrot(iddd2,9) + DDrot_dA(iddd1,8)*(DDrot(iddd2,2) - DDrot(iddd2,3)) + DDrot_dA(iddd1,9)*DDrot(iddd2,10));
            block(30 + iddd1,30 + iddd2) += intn[16]*((DDrot(iddd1,2) - DDrot(iddd1,3))*DDrot_dA(iddd2,8) + DDrot(iddd1,10)*DDrot_dA(iddd2,9) + DDrot(iddd1,8)*(DDrot_dA(iddd2,2) - DDrot_dA(iddd2,3)) + DDrot(iddd1,9)*DDrot_dA(iddd2,10));
          }
        }
      }
    }
  }
  void dFockdX(int der, matrixE & FockA, matrixE & FockB, matrixE & dens, matrixE & qdens, matrixE & hcoreCD, double RCD, double cost, double sint, double cosp, double sinp, int iC, int iD, int atmC, int atmD, int posC, int posD) {
    //function that calculates the first-static-derivatives of the Fock matrix; note that this function cannot be used for the full calculation of the Fock matrix because it only loops over two atoms
    //the derivative variable controls which derivatives are calculated, 1 -> with respect to RCD, 2 -> with respect to thetaCD, 3 -> with respect to phiCD
    //results in atomic units
    size_t natoms = 2;
    int cntA = 0;
    int cntB = 0;
    int im;
    int in;
    int il;
    int is;
    int nAOsC = 1;
    int nAOsD = 1;
    if (atmC > 1) {
      nAOsC = 4;
      if (AtomWithDOrbitals(atmC)) {nAOsC = 9;}
    }
    if (atmD > 1) {
      nAOsD = 4;
      if (AtomWithDOrbitals(atmD)) {nAOsD = 9;}
    }
    chg[0] = double(atmD) - double(CoreCharge[iD]);
    chg[1] = double(atmC) - double(CoreCharge[iC]);
    nintegrals[0] = nAOsC*(nAOsC + 1)/2;
    nintegrals[1] = nAOsD*(nAOsD + 1)/2;
    atomslocal[0] = atmC;
    atomslocal[1] = atmD;
    ipos[0] = posC;
    ipos[1] = posD;
    AOs_local[0] = nAOsC;
    AOs_local[1] = nAOsD;
    std::vector<matrixE> integrals;
    matrixE block(nintegrals[0],1);
    integrals.push_back(block);
    block.resize(nintegrals[1],1);
    integrals.push_back(block);
    block.resize(nintegrals[0],nintegrals[1]);
    bool transitionmetalA;
    bool transitionmetalB;
    double factor = 1.0;
    if (der == 1) {
      IntegralBlock2C_dR(0,block,atmC,atmD,RCD,cost,sint,cosp,sinp);
      IntegralBlock2C_dR(1,integrals[0],atmC,atmD,RCD,cost,sint,cosp,sinp);
      IntegralBlock2C_dR(1,integrals[1],atmD,atmC,RCD,-cost,sint,-cosp,-sinp);
    }
    else if (der > 1) {
      IntegralBlock2C_dA(0,block,der,atmC,atmD,RCD,cost,sint,cosp,sinp);
      IntegralBlock2C_dA(1,integrals[0],der,atmC,atmD,RCD,cost,sint,cosp,sinp);
      IntegralBlock2C_dA(1,integrals[1],der,atmD,atmC,RCD,-cost,sint,-cosp,-sinp);
    }
    //zero Fock matrix and resize
    FockA.resize(nAOsC + nAOsD,nAOsC + nAOsD);
    FockA.zero();
    if (shell == "open") {
      FockB.resize(nAOsC + nAOsD,nAOsC + nAOsD);
      FockB.zero();
    }
    for (size_t idA = 0; idA < natoms; ++idA) {
      transitionmetalA = TransitionMetal(atomslocal[idA]);
      factor = 1.0;
      if ((idA == 1)&&(der == 2)) {factor = -1.0;}
      for (size_t imu = 0; imu < AOs_local[idA]; ++imu) {
        //get position in gamma matrix
        im = imu;
        if (transitionmetalA) {
          if (imu > 4) {im = imu - 5;}
          else if (imu < 5) {im = imu + 4;}
        }
        pos[0] = posgamma(im,im);
        for (size_t inu = 0; inu < AOs_local[idA]; ++inu) {
          in = inu;
          if (transitionmetalA) {
            if (inu > 4) {in = inu - 5;}
            else if (inu < 5) {in = inu + 4;}
          }
          pos[1] = posgamma(im,in);
          FockA(cntA + imu + 1,cntA + inu + 1) -= factor*chg[idA]*integrals[idA](pos[1] + 1,1);
          if (shell == "open") {FockB(cntA + imu + 1,cntA + inu + 1) -= factor*chg[idA]*integrals[idA](pos[1] + 1,1);}
          cntB = 0;
          for (size_t idB = 0; idB < natoms; ++idB) {
            transitionmetalB = TransitionMetal(atomslocal[idB]);
            if (idB != idA) {
              for (size_t ilambda = 0; ilambda < AOs_local[idB]; ++ilambda) {
                il = ilambda;
                if (transitionmetalB) {
                  if (ilambda > 4) {il = ilambda - 5;}
                  else if (ilambda < 5) {il = ilambda + 4;}
                }
                if (imu == inu) {
                  FockA(cntA + imu + 1,cntB + ilambda + 1) += hcoreCD(imu*(idA == 0) + ilambda*(idA == 1) + 1,imu*(idA == 1) + ilambda*(idA == 0) + 1);
                  if (shell == "open") {FockB(cntA + imu + 1,cntB + ilambda + 1) += hcoreCD(imu*(idA == 0) + ilambda*(idA == 1) + 1,imu*(idA == 1) + ilambda*(idA == 0) + 1);}
                }
                for (size_t isigma = 0; isigma < AOs_local[idB]; ++isigma) {
                  is = isigma;
                  if (transitionmetalB) {
                    if (isigma > 4) {is = isigma - 5;}
                    else if (isigma < 5) {is = isigma + 4;}
                  }
                  pos[2] = posgamma(il,is);
                  FockA(cntA + imu + 1,cntB + ilambda + 1) -= 0.5*dens(ipos[idA] + inu + 1,ipos[idB] + isigma + 1)*block(pos[1]*(idA == 0) + pos[2]*(idA == 1) + 1,pos[2]*(idA == 0) + pos[1]*(idA == 1) + 1);
                  if (shell == "open") {
                    FockA(cntA + imu + 1,cntB + ilambda + 1) -= 0.5*qdens(ipos[idA] + inu + 1,ipos[idB] + isigma + 1)*block(pos[1]*(idA == 0) + pos[2]*(idA == 1) + 1,pos[2]*(idA == 0) + pos[1]*(idA == 1) + 1);
                    FockB(cntA + imu + 1,cntB + ilambda + 1) -= 0.5*(dens(ipos[idA] + inu + 1,ipos[idB] + isigma + 1) - qdens(ipos[idA] + inu + 1,ipos[idB] + isigma + 1))*block(pos[1]*(idA == 0) + pos[2]*(idA == 1) + 1,pos[2]*(idA == 0) + pos[1]*(idA == 1) + 1);
                  }
                  FockA(cntA + imu + 1,cntA + inu + 1) += dens(ipos[idB] + ilambda + 1,ipos[idB] + isigma + 1)*block(pos[1]*(idA == 0) + pos[2]*(idA == 1) + 1,pos[2]*(idA == 0) + pos[1]*(idA == 1) + 1);
                  if (shell == "open") {
                    FockB(cntA + imu + 1,cntA + inu + 1) += dens(ipos[idB] + ilambda + 1,ipos[idB] + isigma + 1)*block(pos[1]*(idA == 0) + pos[2]*(idA == 1) + 1,pos[2]*(idA == 0) + pos[1]*(idA == 1) + 1);
                  }
                }
              }
            }
            cntB += AOs_local[idB];
          }
        }
      }
      cntA += AOs_local[idA];
    }
  }
  void dFockdXTEST(int der, matrixE & FockA, matrixE & FockB, matrixE & dens, matrixE & qdens, matrixE & hcoreCD, double RCD, double cost, double sint, double cosp, double sinp, int iC, int iD, int atmC, int atmD, int posC, int posD) {
    size_t natoms = 2;
    int cntA = 0;
    int cntB = 0;
    int im;
    int in;
    int il;
    int is;
    int nAOsC = 1;
    int nAOsD = 1;
    if (atmC > 1) {
      nAOsC = 4;
      if (AtomWithDOrbitals(atmC)) {nAOsC = 9;}
    }
    if (atmD > 1) {
      nAOsD = 4;
      if (AtomWithDOrbitals(atmD)) {nAOsD = 9;}
    }
    chg[0] = double(atmD) - double(CoreCharge[iD]);
    chg[1] = double(atmC) - double(CoreCharge[iC]);
    nintegrals[0] = nAOsC*(nAOsC + 1)/2;
    nintegrals[1] = nAOsD*(nAOsD + 1)/2;
    atomslocal[0] = atmC;
    atomslocal[1] = atmD;
    ipos[0] = posC;
    ipos[1] = posD;
    AOs_local[0] = nAOsC;
    AOs_local[1] = nAOsD;
    std::vector<matrixE> integrals;
    matrixE block(nintegrals[0],1);
    integrals.push_back(block);
    block.resize(nintegrals[1],1);
    integrals.push_back(block);
    block.resize(nintegrals[0],nintegrals[1]);
    bool transitionmetalA;
    bool transitionmetalB;
    double factor = 1.0;
    IntegralBlock2C(0,block,atmC,atmD,RCD,cost,sint,cosp,sinp);
    IntegralBlock2C(1,integrals[0],atmC,atmD,RCD,cost,sint,cosp,sinp);
    IntegralBlock2C(1,integrals[1],atmD,atmC,RCD,-cost,sint,-cosp,-sinp);
    //zero Fock matrix and resize
    FockA.resize(nAOsC + nAOsD,nAOsC + nAOsD);
    FockA.zero();
    if (shell == "open") {
      FockB.resize(nAOsC + nAOsD,nAOsC + nAOsD);
      FockB.zero();
    }
    for (size_t idA = 0; idA < natoms; ++idA) {
      transitionmetalA = TransitionMetal(atomslocal[idA]);
      factor = 1.0;
      if ((idA == 1)&&(der == 2)) {factor = -1.0;}
      for (size_t imu = 0; imu < AOs_local[idA]; ++imu) {
        //get position in gamma matrix
        im = imu;
        if (transitionmetalA) {
          if (imu > 4) {im = imu - 5;}
          else if (imu < 5) {im = imu + 4;}
        }
        pos[0] = posgamma(im,im);
        for (size_t inu = 0; inu < AOs_local[idA]; ++inu) {
          in = inu;
          if (transitionmetalA) {
            if (inu > 4) {in = inu - 5;}
            else if (inu < 5) {in = inu + 4;}
          }
          pos[1] = posgamma(im,in);
          FockA(cntA + imu + 1,cntA + inu + 1) -= factor*chg[idA]*integrals[idA](pos[1] + 1,1);
          if (shell == "open") {FockB(cntA + imu + 1,cntA + inu + 1) -= factor*chg[idA]*integrals[idA](pos[1] + 1,1);}
          cntB = 0;
          for (size_t idB = 0; idB < natoms; ++idB) {
            transitionmetalB = TransitionMetal(atomslocal[idB]);
            if (idB != idA) {
              for (size_t ilambda = 0; ilambda < AOs_local[idB]; ++ilambda) {
                il = ilambda;
                if (transitionmetalB) {
                  if (ilambda > 4) {il = ilambda - 5;}
                  else if (ilambda < 5) {il = ilambda + 4;}
                }
                if (imu == inu) {
                  FockA(cntA + imu + 1,cntB + ilambda + 1) += hcoreCD(imu*(idA == 0) + ilambda*(idA == 1) + 1,imu*(idA == 1) + ilambda*(idA == 0) + 1);
                  if (shell == "open") {FockB(cntA + imu + 1,cntB + ilambda + 1) += hcoreCD(imu*(idA == 0) + ilambda*(idA == 1) + 1,imu*(idA == 1) + ilambda*(idA == 0) + 1);}
                }
                for (size_t isigma = 0; isigma < AOs_local[idB]; ++isigma) {
                  is = isigma;
                  if (transitionmetalB) {
                    if (isigma > 4) {is = isigma - 5;}
                    else if (isigma < 5) {is = isigma + 4;}
                  }
                  pos[2] = posgamma(il,is);
                  FockA(cntA + imu + 1,cntB + ilambda + 1) -= 0.5*dens(ipos[idA] + inu + 1,ipos[idB] + isigma + 1)*block(pos[1]*(idA == 0) + pos[2]*(idA == 1) + 1,pos[2]*(idA == 0) + pos[1]*(idA == 1) + 1);
                  if (shell == "open") {
                    FockA(cntA + imu + 1,cntB + ilambda + 1) -= 0.5*qdens(ipos[idA] + inu + 1,ipos[idB] + isigma + 1)*block(pos[1]*(idA == 0) + pos[2]*(idA == 1) + 1,pos[2]*(idA == 0) + pos[1]*(idA == 1) + 1);
                    FockB(cntA + imu + 1,cntB + ilambda + 1) -= 0.5*(dens(ipos[idA] + inu + 1,ipos[idB] + isigma + 1) - qdens(ipos[idA] + inu + 1,ipos[idB] + isigma + 1))*block(pos[1]*(idA == 0) + pos[2]*(idA == 1) + 1,pos[2]*(idA == 0) + pos[1]*(idA == 1) + 1);
                  }
                  FockA(cntA + imu + 1,cntA + inu + 1) += dens(ipos[idB] + ilambda + 1,ipos[idB] + isigma + 1)*block(pos[1]*(idA == 0) + pos[2]*(idA == 1) + 1,pos[2]*(idA == 0) + pos[1]*(idA == 1) + 1);
                  if (shell == "open") {
                    FockB(cntA + imu + 1,cntA + inu + 1) += dens(ipos[idB] + ilambda + 1,ipos[idB] + isigma + 1)*block(pos[1]*(idA == 0) + pos[2]*(idA == 1) + 1,pos[2]*(idA == 0) + pos[1]*(idA == 1) + 1);
                  }
                }
              }
            }
            cntB += AOs_local[idB];
          }
        }
      }
      cntA += AOs_local[idA];
    }
  }
  int posgamma(int mu1, int nu1) {
    //function returning the gamma matrix element for the orbital pair (mu,nu)
    int mu = std::min(mu1,nu1);
    int nu = std::max(mu1,nu1);
    int pos = 0;
    if (mu == 0) {
      if (nu == 1) {pos = 1;}       //s,px
      else if (nu == 2) {pos = 2;}  //s,py
      else if (nu == 3) {pos = 3;}  //s,pz
      else if (nu == 4) {pos = 10;} //s,dxz
      else if (nu == 5) {pos = 11;} //s,dyz
      else if (nu == 6) {pos = 12;} //s,dz2
      else if (nu == 7) {pos = 13;} //s,dxy
      else if (nu == 8) {pos = 14;} //s,dx2-y2
    }
    else if (mu == 1) {
      if (nu == 1) {pos = 4;}       //px,px
      else if (nu == 2) {pos = 5;}  //px,py
      else if (nu == 3) {pos = 6;}  //px,pz
      else if (nu == 4) {pos = 15;} //px,dxz
      else if (nu == 5) {pos = 16;} //px,dyz
      else if (nu == 6) {pos = 17;} //px,dz2
      else if (nu == 7) {pos = 18;} //px,dxy
      else if (nu == 8) {pos = 19;} //px,dx2-y2
    }
    else if (mu == 2) {
      if (nu == 2) {pos = 7;}       //py,py
      else if (nu == 3) {pos = 8;}  //py,pz
      else if (nu == 4) {pos = 20;} //py,dxz
      else if (nu == 5) {pos = 21;} //py,dyz
      else if (nu == 6) {pos = 22;} //py,dz2
      else if (nu == 7) {pos = 23;} //py,dxy
      else if (nu == 8) {pos = 24;} //py,dx2-y2
    }
    else if (mu == 3) {
      if (nu == 3) {pos = 9;}       //pz,pz
      else if (nu == 4) {pos = 25;} //pz,dxz
      else if (nu == 5) {pos = 26;} //pz,dyz
      else if (nu == 6) {pos = 27;} //pz,dz2
      else if (nu == 7) {pos = 28;} //pz,dxy
      else if (nu == 8) {pos = 29;} //pz,dx2-y2
    }
    else if (mu == 4) {
      if (nu == 4) {pos = 30;}      //dxz,dxz
      else if (nu == 5) {pos = 31;} //dxz,dyz
      else if (nu == 6) {pos = 32;} //dxz,dz2
      else if (nu == 7) {pos = 33;} //dxz,dxy
      else if (nu == 8) {pos = 34;} //dxz,dx2-y2
    }
    else if (mu == 5) {
      if (nu == 5) {pos = 35;}      //dyz,dyz
      else if (nu == 6) {pos = 36;} //dyz,dz2
      else if (nu == 7) {pos = 37;} //dyz,dxy
      else if (nu == 8) {pos = 38;} //dyz,dx2-y2
    }
    else if (mu == 6) {
      if (nu == 6) {pos = 39;}      //dz2,dz2
      else if (nu == 7) {pos = 40;} //dz2,dxy
      else if (nu == 8) {pos = 41;} //dz2,dx2-y2
    }
    else if (mu == 7) {
      if (nu == 7) {pos = 42;}      //dxy,dxy
      else if (nu == 8) {pos = 43;} //dxy,dx2-y2
    }
    else if (mu == 8) {
      if (nu == 8) {pos = 44;}      //dx2-y2,dx2-y2
    }
    return pos;
  }
  virtual double gfactor(size_t iatm1, size_t iatm2, double RAB) {
    //function determining special factors in the nuclear repulsion terms of MNDO
    double gfac = 1.0;
    if (iatm2 == 1) {
      if ((iatm1 == 7)||(iatm1 == 8)) {gfac = RAB*dist_Angstrom2au;}
    }
    return gfac;
  }
  virtual double gfactor_dR(size_t iatm1, size_t iatm2) {
    //function determining the derivative of the special factors in the nuclear repulsion terms of MNDO
    //with respect to internuclear distance
    double gfac = 0.0;
    if (iatm2 == 1) {
      if ((iatm1 == 7)||(iatm1 == 8)) {gfac = dist_Angstrom2au;}
    }
    return gfac;
  }
  virtual double AM1factor(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {return 0.0;}
  virtual double AM1factor_dR(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {return 0.0;}
  virtual double AM1factor_dR2(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {return 0.0;}
  double EPDDG(int atmA, int atmB, double RAB) {
    //function calculating the PDDG contribution to the nuclear energy
    double nA = ValenceElectrons(atmA);
    double nB = ValenceElectrons(atmB);
    double PAi;
    double DAi;
    double epddg = 0.0;
    double aux;
    for (size_t idi = 1; idi < 3; ++idi) {
      PAi = PAPDDG(atmA,idi);
      DAi = DAPDDG(atmA,idi);
      for (size_t idj = 1; idj < 3; ++idj) {
        aux = RAB - DAi - DAPDDG(atmB,idj);
        epddg += (PAi*nA + PAPDDG(atmB,idj)*nB)*exp(-10.0*aux*aux);
      }
    }
    return epddg/(nA + nB);
  }
  double dEPDDG(int atmA, int atmB, double RAB) {
    //function calculating the PDDG contribution to the nuclear energy
    double nA = ValenceElectrons(atmA);
    double nB = ValenceElectrons(atmB);
    double PAi;
    double DAi;
    double depddg = 0.0;
    double aux;
    for (size_t idi = 1; idi < 3; ++idi) {
      PAi = PAPDDG(atmA,idi);
      DAi = DAPDDG(atmA,idi);
      for (size_t idj = 1; idj < 3; ++idj) {
        aux = RAB - DAi - DAPDDG(atmB,idj);
        depddg += (PAi*nA + PAPDDG(atmB,idj)*nB)*exp(-10.0*aux*aux)*aux;
      }
    }
    return -20.0*depddg/(nA + nB);
  }
  double d2EPDDG(int atmA, int atmB, double RAB) {
    //function calculating the PDDG contribution to the nuclear energy
    double nA = ValenceElectrons(atmA);
    double nB = ValenceElectrons(atmB);
    double PAi;
    double DAi;
    double d2epddg = 0.0;
    double aux;
    for (size_t idi = 1; idi < 3; ++idi) {
      PAi = PAPDDG(atmA,idi);
      DAi = DAPDDG(atmA,idi);
      for (size_t idj = 1; idj < 3; ++idj) {
        aux = RAB - DAi - DAPDDG(atmB,idj);
        d2epddg += (PAi*nA + PAPDDG(atmB,idj)*nB)*exp(-10.0*aux*aux)*(1.0 - 20.0*aux*aux);
      }
    }
    return -20.0*d2epddg/(nA + nB);
  }
  matrixE EnuclearMat() {
    //function that calculates the nuclear energy
    double XAB;
    double alphaA;
    double chgA;
    double alphaB;
    double chgB;
    double RAB;
    double integral;
    double aux;
    double enucaux = 0.0;
    double gA;
    double gB;
    double factorA;
    double factorB;
    double factorC;
    double unitconv = 1.0/(au2eV*dist_Angstrom2au);
    matrixE Enuclear(Natoms,Natoms);
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      //charge and nuclear energy variables
      KA = AM1K(atoms[iatm]);
      LA = AM1L(atoms[iatm]);
      MA = AM1M(atoms[iatm]);
      chgA = double(int(atoms[iatm]) - int(CoreCharge[iatm]));
      Enuclear(iatm + 1,iatm + 1) = 0.0;
      for (size_t ibtm = iatm + 1; ibtm < Natoms; ++ibtm) {
        alphaA = alpha(atoms[iatm],atoms[ibtm]);
        XAB = xAB(atoms[iatm],atoms[ibtm]);
        //charge and nuclear energy variables
        KB = AM1K(atoms[ibtm],atoms[iatm]);
        LB = AM1L(atoms[ibtm],atoms[iatm]);
        MB = AM1M(atoms[ibtm],atoms[iatm]);
        if (atoms[iatm] == 5) {
          KA = AM1K(atoms[iatm],atoms[ibtm]);
          LA = AM1L(atoms[iatm],atoms[ibtm]);
          MA = AM1M(atoms[iatm],atoms[ibtm]);
        }
        chgB = double(int(atoms[ibtm]) - int(CoreCharge[ibtm]));
        //interatomic distance
        RAB = mol.AUdistance(iatm + 1,ibtm + 1);
        //integral
        integral = eri2Center(0,0,0,0,0,0,0,0,RAB,D,atoms[iatm],atoms[ibtm],2);
        factorA = 1.0;
        factorB = 0.0;
        factorC = 0.0;
        if (XAB > 0.0) {           //PM6
          //assume here that these parameters take values considerably different from zero
          if (((atoms[iatm] == 1)&&(atoms[ibtm] == 6))||((atoms[iatm] == 6)&&(atoms[ibtm] == 1))||((atoms[iatm] == 1)&&(atoms[ibtm] == 7))||((atoms[iatm] == 7)&&(atoms[ibtm] == 1))||((atoms[iatm] == 1)&&(atoms[ibtm] == 8))||((atoms[iatm] == 8)&&(atoms[ibtm] == 1))) {
            factorA += 2.0*XAB*exp(-alphaA*RAB*RAB*dist_Angstrom2au);
          }
          else {
            factorA += 2.0*XAB*exp(-alphaA*RAB*(1.0 + 0.0003*RAB*RAB*RAB*RAB*RAB*dist_Angstrom2au*dist_Angstrom2au*dist_Angstrom2au*dist_Angstrom2au*dist_Angstrom2au));
          }
          if ((atoms[iatm] == 6)&&(atoms[ibtm] == 6)) {factorA += 9.28*exp(-5.98*RAB*dist_Angstrom2au);}
          else if (((atoms[iatm] == 8)&&(atoms[ibtm] == 14))||((atoms[iatm] == 14)&&(atoms[ibtm] == 8))) {
            factorA -= 0.0007*exp(-(RAB*dist_Angstrom2au - 2.9)*(RAB*dist_Angstrom2au - 2.9));
          }
          factorB = 1.0e-8*pow((pow(chgA,1.0/3.0) + pow(chgB,1.0/3.0))/(RAB*dist_Angstrom2au),12.0)/au2eV;
        }
        gA = gfactor(atoms[iatm],atoms[ibtm],RAB);
        gB = gfactor(atoms[ibtm],atoms[iatm],RAB);
        alphaB = alpha(atoms[ibtm],atoms[iatm]);
        factorA += gA*exp(-alphaA*RAB*dist_Angstrom2au) + gB*exp(-alphaB*RAB*dist_Angstrom2au);
        if (fabs(KA[0]) > 1.0e-5) {factorC += AM1factor(RAB,KA,LA,MA);}
        if (fabs(KB[0]) > 1.0e-5) {factorC += AM1factor(RAB,KB,LB,MB);}
        aux = chgA*chgB*integral*factorA + factorB + chgA*chgB*factorC;
        if (doPDDG) {aux += EPDDG(atoms[iatm],atoms[ibtm],RAB*dist_Angstrom2au);}
        enucaux += aux;
        Enuclear(iatm + 1,ibtm + 1) = aux;
        Enuclear(ibtm + 1,iatm + 1) = aux;
      }
    }
    return Enuclear;
  }
  void Enuclear_dR(std::vector<double> & enuclear_dr) {
    //function that calculates the derivative of the nuclear energy with respect to internuclear distance
    size_t index = 0;
    enuclear_dr.resize(Natoms*(Natoms - 1)/2);
    for (size_t idx = 0; idx < Natoms*(Natoms - 1)/2; ++idx) {
      enuclear_dr[idx] = 0.0;
    }
    double alphaA;
    double chgA;
    double alphaB;
    double chgB;
    double XAB;
    double RAB;
    double RABangstroem;
    double integral;
    double integral_dR;
    double aux;
    double gA;
    double gB;
    double gAdR;
    double gBdR;
    double expA;
    double expB;
    double fac1;
    double fac2;
    double fac3;
    double fac4;
    double A2au5 = dist_Angstrom2au*dist_Angstrom2au*dist_Angstrom2au*dist_Angstrom2au*dist_Angstrom2au;
    double RAB5;
    double unitconv = 1.0/(au2eV*dist_Angstrom2au);
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      //charge and nuclear energy variables
      KA = AM1K(atoms[iatm]);
      LA = AM1L(atoms[iatm]);
      MA = AM1M(atoms[iatm]);
      chgA = double(int(atoms[iatm]) - int(CoreCharge[iatm]));
      for (size_t ibtm = iatm + 1; ibtm < Natoms; ++ibtm) {
        alphaA = alpha(atoms[iatm],atoms[ibtm]);
        alphaB = alpha(atoms[ibtm],atoms[iatm]);
        //charge and nuclear energy variables
        XAB = xAB(atoms[iatm],atoms[ibtm]);
        KB = AM1K(atoms[ibtm],atoms[iatm]);
        LB = AM1L(atoms[ibtm],atoms[iatm]);
        MB = AM1M(atoms[ibtm],atoms[iatm]);
        if (atoms[iatm] == 5) {
          KA = AM1K(atoms[iatm],atoms[ibtm]);
          LA = AM1L(atoms[iatm],atoms[ibtm]);
          MA = AM1M(atoms[iatm],atoms[ibtm]);
        }
        chgB = double(int(atoms[ibtm]) - int(CoreCharge[ibtm]));
        //interatomic distance
        RAB = mol.AUdistance(iatm + 1,ibtm + 1);
        RABangstroem = RAB*dist_Angstrom2au;
        //integrals required
        integral = eri2Center(0,0,0,0,0,0,0,0,RAB,D,atoms[iatm],atoms[ibtm],2);
        integral_dR = eri2Center_dR(0,0,0,0,0,0,0,0,RAB,D,atoms[iatm],atoms[ibtm],2);
        //auxiliaries
        fac1 = 1.0;
        fac2 = 0.0;
        fac3 = 0.0;
        fac4 = 0.0;
        if (XAB > 0.0) {           //PM6
          if (((atoms[iatm] == 1)&&(atoms[ibtm] == 6))||((atoms[iatm] == 6)&&(atoms[ibtm] == 1))||((atoms[iatm] == 1)&&(atoms[ibtm] == 7))||((atoms[iatm] == 7)&&(atoms[ibtm] == 1))||((atoms[iatm] == 1)&&(atoms[ibtm] == 8))||((atoms[iatm] == 8)&&(atoms[ibtm] == 1))) {
            expA = 2.0*exp(-alphaA*RAB*RABangstroem);
            fac1 += XAB*expA;
            fac2 += 2.0*XAB*alphaA*expA*RABangstroem;
          }
          else {
            RAB5 = RAB*RAB*RAB*RAB*RAB*A2au5;
            expA = 2.0*exp(-alphaA*RAB*(1.0 + 0.0003*RAB5));
            fac1 += XAB*expA;
            fac2 += XAB*alphaA*expA*(1.0 + 0.0018*RAB5);
          }
          if ((atoms[iatm] == 6)&&(atoms[ibtm] == 6)) {
            expA = exp(-5.98*RABangstroem);
            fac1 += 9.28*expA;
            fac2 += 55.4944*expA*dist_Angstrom2au;
          }
          else if (((atoms[iatm] == 8)&&(atoms[ibtm] == 14))||((atoms[iatm] == 14)&&(atoms[ibtm] == 8))) {
            aux = RABangstroem - 2.9;
            expA = exp(-aux*aux);
            fac1 -= 0.0007*expA;
            fac2 -= 0.0014*aux*expA*dist_Angstrom2au;
          }
          fac3 = 1.2e-7*pow((pow(chgA,1.0/3.0) + pow(chgB,1.0/3.0))/(RABangstroem),12.0)/(au2eV*RAB);
        }
        gA = gfactor(atoms[iatm],atoms[ibtm],RAB);
        gB = gfactor(atoms[ibtm],atoms[iatm],RAB);
        if ((gA > 0.0)||(gB > 0.0)) {
          expA = exp(-alphaA*RABangstroem);
          expB = exp(-alphaB*RABangstroem);
          gAdR = gfactor_dR(atoms[iatm],atoms[ibtm]);
          gBdR = gfactor_dR(atoms[ibtm],atoms[iatm]);
          fac1 += gA*expA + gB*expB;
          fac2 -= gAdR*expA - gA*alphaA*dist_Angstrom2au*expA + gBdR*expB - gB*alphaB*dist_Angstrom2au*expB;
        }
        if (fabs(KA[0]) > 1.0e-5) {fac4 += AM1factor_dR(RAB,KA,LA,MA);}
        if (fabs(KB[0]) > 1.0e-5) {fac4 += AM1factor_dR(RAB,KB,LB,MB);}
        aux = chgA*chgB*(integral_dR*fac1 - integral*fac2) - fac3 - chgA*chgB*fac4;
        if (doPDDG) {aux += dEPDDG(atoms[iatm],atoms[ibtm],RABangstroem);}
        enuclear_dr[index] = aux;
        ++index;
      }
    }
  }
  virtual void Enuclear_dR2(std::vector<double> & enuclear_dr2) {
    //function that calculates the second-derivative of the nuclear energy with respect to internuclear distance
    size_t index = 0;
    enuclear_dr2.resize(Natoms*(Natoms - 1)/2);
    for (size_t idx = 0; idx < Natoms*(Natoms - 1)/2; ++idx) {
      enuclear_dr2[idx] = 0.0;
    }
    double alphaA;
    double chgA;
    double alphaB;
    double chgB;
    double RAB;
    double XAB;
    double integral;
    double integral_dR;
    double integral_dR2;
    double aux;
    double gA;
    double gB;
    double gAdR;
    double gBdR;
    double expA;
    double expB;
    double fac1;
    double fac2;
    double fac3;
    double fac4;
    double fac5;
    double RAB5;
    double RR;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      //charge and nuclear energy variables
      KA = AM1K(atoms[iatm]);
      LA = AM1L(atoms[iatm]);
      MA = AM1M(atoms[iatm]);
      chgA = double(int(atoms[iatm]) - int(CoreCharge[iatm]));
      for (size_t ibtm = iatm + 1; ibtm < Natoms; ++ibtm) {
        alphaA = alpha(atoms[iatm],atoms[ibtm]);
        //charge and nuclear energy variables
        XAB = xAB(atoms[iatm],atoms[ibtm]);
        alphaB = alpha(atoms[ibtm],atoms[iatm]);
        KB = AM1K(atoms[ibtm],atoms[iatm]);
        LB = AM1L(atoms[ibtm],atoms[iatm]);
        MB = AM1M(atoms[ibtm],atoms[iatm]);
        if (atoms[iatm] == 5) {
          KA = AM1K(atoms[iatm],atoms[ibtm]);
          LA = AM1L(atoms[iatm],atoms[ibtm]);
          MA = AM1M(atoms[iatm],atoms[ibtm]);
        }
        chgB = double(int(atoms[ibtm]) - int(CoreCharge[ibtm]));
        //interatomic distance
        RAB = mol.AUdistance(iatm + 1,ibtm + 1);
        RR = RAB*dist_Angstrom2au;
        //integrals required
        integral = eri2Center(0,0,0,0,0,0,0,0,RAB,D,atoms[iatm],atoms[ibtm],2);
        integral_dR = eri2Center_dR(0,0,0,0,0,0,0,0,RAB,D,atoms[iatm],atoms[ibtm],2);
        integral_dR2 = eri2Center_dR2(0,0,0,0,0,0,0,0,RAB,D,atoms[iatm],atoms[ibtm],2);
        //auxiliaries
        fac1 = 1.0;                                     //factor assiciated to (mn|sl) dR2
        fac2 = 0.0;                                     //factor assiciated to (mn|sl) dR
        fac3 = 0.0;                                     //factor assiciated to (mn|sl)
        fac4 = 0.0;
        fac5 = 0.0;
        if (XAB > 0.0) {           //PM6
          if (((atoms[iatm] == 1)&&(atoms[ibtm] == 6))||((atoms[iatm] == 6)&&(atoms[ibtm] == 1))||((atoms[iatm] == 1)&&(atoms[ibtm] == 7))||((atoms[iatm] == 7)&&(atoms[ibtm] == 1))||((atoms[iatm] == 1)&&(atoms[ibtm] == 8))||((atoms[iatm] == 8)&&(atoms[ibtm] == 1))) {
            expA = 2.0*exp(-alphaA*RAB*RR);
            fac1 += XAB*expA;
            fac2 += 4.0*XAB*alphaA*expA*RR;
            fac3 += 2.0*XAB*alphaA*expA*(dist_Angstrom2au - 2.0*alphaA*RR*RR);
          }
          else {
            RAB5 = RR*RR*RR*RR*RR;
            expA = 2.0*exp(-alphaA*RAB*(1.0 + 0.0003*RAB5));
            fac1 += XAB*expA;
            fac2 += 2.0*XAB*alphaA*expA*(1.0 + 0.0018*RAB5);
            fac3 += XAB*alphaA*expA*(0.009*RAB5/RAB - alphaA*(1.0 + 0.0018*RAB5)*(1.0 + 0.0018*RAB5));
          }
          if ((atoms[iatm] == 6)&&(atoms[ibtm] == 6)) {
            expA = exp(-5.98*RR);
            fac1 += 9.28*expA;
            fac2 += 110.9888*expA*dist_Angstrom2au;
            fac3 -= 331.856512*expA*dist_Angstrom2au*dist_Angstrom2au;
          }
          else if (((atoms[iatm] == 8)&&(atoms[ibtm] == 14))||((atoms[iatm] == 14)&&(atoms[ibtm] == 8))) {
            aux = RR - 2.9;
            expA = exp(-aux*aux);
            fac1 -= 0.0007*expA;
            fac2 -= 0.0028*aux*expA*dist_Angstrom2au;
            fac3 -= 0.0014*expA*(1.0 - 2.0*aux*aux)*dist_Angstrom2au*dist_Angstrom2au;
          }
          fac4 = 1.56e-6*pow((pow(chgA,1.0/3.0) + pow(chgB,1.0/3.0))/RR,12.0)/(au2eV*RAB*RAB);
        }
        gA = gfactor(atoms[iatm],atoms[ibtm],RAB);
        gB = gfactor(atoms[ibtm],atoms[iatm],RAB);
        if ((gA > 0.0)||(gB > 0.0)) {
          gAdR = gfactor_dR(atoms[iatm],atoms[ibtm]);
          gBdR = gfactor_dR(atoms[ibtm],atoms[iatm]);
          expA = exp(-alphaA*RR);
          expB = exp(-alphaB*RR);
          fac1 += gA*expA + gB*expB;
          fac2 -= 2.0*((gAdR - alphaA*gA*dist_Angstrom2au)*expA + (gBdR - alphaB*gB*dist_Angstrom2au)*expB);
          fac3 -= (alphaA*(alphaA*gA*dist_Angstrom2au - 2.0*gAdR)*expA + alphaB*(alphaB*gB*dist_Angstrom2au - 2.0*gBdR)*expB)*dist_Angstrom2au;
        }
        if (fabs(KA[0]) > 1.0e-5) {fac5 += AM1factor_dR2(RAB,KA,LA,MA);}
        if (fabs(KB[0]) > 1.0e-5) {fac5 += AM1factor_dR2(RAB,KB,LB,MB);}
        aux = chgA*chgB*(integral_dR2*fac1 - integral_dR*fac2 - integral*fac3) + fac4 - chgA*chgB*fac5;
        if (doPDDG) {aux += d2EPDDG(atoms[iatm],atoms[ibtm],RR);}
        enuclear_dr2[index] = aux;
        ++index;
      }
    }
  }
  void se_integrals(matrixE & geometry, double tolerance = 1e-6) {
    //function calculating ALL the required integrals according to Dewar's semiempirical formulas
    //get total number of orbital pairs in system
    int irow = 0;
    int icol = 0;
    int naosA;
    int naosB;
    int nintA;
    int nintB;
    int AOpair = 0;
    for (size_t iatm = 0; iatm < Natoms; ++ iatm) {
      AOpair += AOs[iatm]*(AOs[iatm] + 1)/2;                    //saving solely unique pairs
    }
    gammaSE.resize(AOpair,AOpair);
    VAB.resize(AOpair,1);
    enuc = 0.0;
    //auxiliaries
    double intn;
    double RAB;
    double RR;
    double cost;
    double sint;
    double cosp;
    double sinp;
    double chgA;
    double chgB;
    double alphaA;
    double alphaB;
    double XAB;
    double gA;
    double gB;
    double factorA;
    double factorB;
    double factorC;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      //charge and nuclear energy variables
      KA = AM1K(atoms[iatm]);
      LA = AM1L(atoms[iatm]);
      MA = AM1M(atoms[iatm]);
      chgA = double(int(atoms[iatm]) - int(CoreCharge[iatm]));
      naosA = AOs[iatm];
      nintA = naosA*(naosA + 1)/2;
      icol = 0;
      for (size_t idr = 1; idr < nintA + 1; ++idr) {
        VAB(irow + idr,1) = 0.0;
      }
      for (size_t ibtm = 0; ibtm < Natoms; ++ibtm) {
        naosB = AOs[ibtm];
        nintB = naosB*(naosB + 1)/2;
        iBlockdR.resize(nintA,nintB);
        if (iatm == ibtm) {
          IntegralBlock1C(iBlockdR,atoms[iatm]);
          for (size_t idr = 1; idr < nintA + 1; ++idr) {
            gammaSE(irow + idr,icol + idr) = iBlockdR(idr,idr);
            for (size_t idc = idr + 1; idc < nintA + 1; ++idc) {
              gammaSE(irow + idr,icol + idc) = iBlockdR(idr,idc);
              gammaSE(icol + idc,irow + idr) = iBlockdR(idr,idc);
            }
          }
        }
        else {
          chgB = double(int(atoms[ibtm]) - int(CoreCharge[ibtm]));
          //getting the orientation vector
          rAB[0] = (geometry(ibtm + 1,1) - geometry(iatm + 1,1));                  //Delta x
          rAB[1] = (geometry(ibtm + 1,2) - geometry(iatm + 1,2));                  //Delta y
          rAB[2] = (geometry(ibtm + 1,3) - geometry(iatm + 1,3));                  //Delta z
          //normalizing it
          RAB = sqrt(rAB[0]*rAB[0] + rAB[1]*rAB[1] + rAB[2]*rAB[2]);
          rAB[0] /= RAB;
          rAB[1] /= RAB;
          rAB[2] /= RAB;
          //getting trigonometric functions for rotations
          cost = rAB[2];
          sint = sqrt(1.0 - cost*cost);
          cosp = 1.0;
          sinp = 0.0;
          if (fabs(sint) > tolerance) {
            cosp = rAB[0]/sint;
            sinp = rAB[1]/sint;
          }
          RR = RAB;
          RAB *= dist_Angstrom2aum1;
          iBlockdRVC.resize(nintA,1);
          IntegralBlock2C(1,iBlockdRVC,atoms[iatm],atoms[ibtm],RAB,cost,sint,cosp,sinp);
          for (size_t idr = 1; idr < nintA + 1; ++idr) {
            VAB(irow + idr,1) -= chgB*iBlockdRVC(idr,1);
          }
          if (ibtm > iatm) {
            //nuclear energy terms
            alphaA = alpha(atoms[iatm],atoms[ibtm]);
            XAB = xAB(atoms[iatm],atoms[ibtm]);
            KB = AM1K(atoms[ibtm],atoms[iatm]);
            LB = AM1L(atoms[ibtm],atoms[iatm]);
            MB = AM1M(atoms[ibtm],atoms[iatm]);
            if (atoms[iatm] == 5) {
              KA = AM1K(atoms[iatm],atoms[ibtm]);
              LA = AM1L(atoms[iatm],atoms[ibtm]);
              MA = AM1M(atoms[iatm],atoms[ibtm]);
            }
            factorA = 1.0;
            factorB = 0.0;
            factorC = 0.0;
            if (XAB > 0.0) {
              //assume here that these parameters take values considerably different from zero
              if (((atoms[iatm] == 1)&&(atoms[ibtm] == 6))||((atoms[iatm] == 6)&&(atoms[ibtm] == 1))||((atoms[iatm] == 1)&&(atoms[ibtm] == 7))||((atoms[iatm] == 7)&&(atoms[ibtm] == 1))||((atoms[iatm] == 1)&&(atoms[ibtm] == 8))||((atoms[iatm] == 8)&&(atoms[ibtm] == 1))) {
                factorA += 2.0*XAB*exp(-alphaA*RAB*RR);
              }
              else {
                factorA += 2.0*XAB*exp(-alphaA*RAB*(1.0 + 0.0003*RR*RR*RR*RR*RR));
              }
              if ((atoms[iatm] == 6)&&(atoms[ibtm] == 6)) {factorA += 9.28*exp(-5.98*RR);}
              else if (((atoms[iatm] == 8)&&(atoms[ibtm] == 14))||((atoms[iatm] == 14)&&(atoms[ibtm] == 8))) {
                factorA -= 0.0007*exp(-(RR - 2.9)*(RR - 2.9));
              }
              factorB = 1.0e-8*pow((pow(chgA,1.0/3.0) + pow(chgB,1.0/3.0))/(RR),12.0)/au2eV;
            }
            gA = gfactor(atoms[iatm],atoms[ibtm],RAB);
            gB = gfactor(atoms[ibtm],atoms[iatm],RAB);
            alphaB = alpha(atoms[ibtm],atoms[iatm]);
            factorA += gA*exp(-alphaA*RR) + gB*exp(-alphaB*RR);
            if (fabs(KA[0]) > 1.0e-5) {factorC += AM1factor(RAB,KA,LA,MA);}
            if (fabs(KB[0]) > 1.0e-5) {factorC += AM1factor(RAB,KB,LB,MB);}
            intn = eri2Center(0,0,0,0,0,0,0,0,RAB,D,atoms[iatm],atoms[ibtm],2);
            enuc += chgA*chgB*intn*factorA + factorB + chgA*chgB*factorC;
            if (doPDDG) {enuc += EPDDG(atoms[iatm],atoms[ibtm],RR);}
            //2-electron integrals
            IntegralBlock2C(0,iBlockdR,atoms[iatm],atoms[ibtm],RAB,cost,sint,cosp,sinp);
            for (size_t idr = 1; idr < nintA + 1; ++idr) {
              for (size_t idc = 1; idc < nintB + 1; ++idc) {
                gammaSE(irow + idr,icol + idc) = iBlockdR(idr,idc);
                gammaSE(icol + idc,irow + idr) = iBlockdR(idr,idc);
              }
            }
          }
        }
        icol += nintB;
      }
      irow += nintA;
    }
    if (print > 0) {
      std::cout << "Matrix Semiempirical Gammas:" << std::endl;
      gammaSE.Print();
      std::cout << "VAB:" << std::endl;
      VAB.Print();
    }
  }
  virtual double SlaterCondonRadialIntegral(int atmnr, int L, int II, int JJ, double threshold = 1e-7) {return 0.0;}
  virtual std::vector<double> AM1K(int atomicnr, int atm2 = 0) {
    std::vector<double> am1k(1,0.0);
    return am1k;
  }
  virtual std::vector<double> AM1L(int atomicnr, int atm2 = 0) {
    std::vector<double> am1l;
    return am1l;
  }
  virtual std::vector<double> AM1M(int atomicnr, int atm2 = 0) {
    std::vector<double> am1m;
    return am1m;
  }
  double eri2Center(int L1, int M1, int L2, int M2, int L3, int M3, int L4, int M4, double RAB, std::vector<double> D,size_t atmnrA,size_t atmnrB, int core = 0) {
    //function that calculates 2 electron repulsion integrals as in the NDDO approximation as applied in MNDO
    //AOs contains AO information on the 4 orbitals involved in the integral, meaning that this is only for s,p basis sets/atoms
    //RAB is distance between centers
    //D is a vector with all necessary D parameters: {D1A, D2A, D1B, D2B}
    //M. J. S. Dewar, W. Thiel, Theoret. Chim. Acta (Berl.), 46, 89, 1977
    double eri = 0.0;
    double rhoA;
    double rhoB;
    int lA;
    int lB;
    int mA;
    int mB;
    chg1 = NonZeroMultipoleMom_sp(L1,M1,L2,M2);
    chg2 = NonZeroMultipoleMom_sp(L3,M3,L4,M4);
    for (size_t mp1 = 0; mp1 < chg1.size(); ++mp1) {
      lA = chg1[mp1][0];
      mA = chg1[mp1][1];
      for (size_t mp2 = 0; mp2 < chg2.size(); ++mp2) {
        lB = chg2[mp2][0];
        mB = chg2[mp2][1];
        if (((abs(mA) % 2 == 0)&&(abs(mB) % 2 == 0))||((abs(mA) % 2 != 0)&&(abs(mB) % 2 != 0))) {
          //even-even or odd-odd only
          if (abs(mA) != abs(mB)) {
            mA = abs(mA);
            mB = abs(mB);
          }
          rhoA = rho(atmnrA,lA);
          rhoB = rho(atmnrB,lB);
          if (core > 0) {rhoB = rhocore(atmnrB);}                        //in this case use the correct core rho for the ss atom
          if ((core > 1)||(core < 0)) {rhoA = rhocore(atmnrA);}          //in this case use the correct core rho for the ss atom
          if (mA != 3) {eri += se_multipole(lA,mA,lB,mB,RAB,rhoA + rhoB,D);}
          else {eri += 0.5*(se_multipole(2,2,2,2,RAB,rhoA + rhoB,D) - se_multipole(2,2,2,-2,RAB,rhoA + rhoB,D));}
        }
      }
    }
    return eri;
  }
  double eri2Center_dR(int L1, int M1, int L2, int M2, int L3, int M3, int L4, int M4, double RAB, std::vector<double> D,size_t atmnrA,size_t atmnrB, int core = 0) {
    //function that calculates the first-derivatives of 2 electron repulsion integrals as in the NDDO approximation as applied in MNDO
    //AOs contains AO information on the 4 orbitals involved in the integral, meaning that this is only for s,p basis sets/atoms
    //RAB is distance between centers
    //D is a vector with all necessary D parameters: {D1A, D2A, D1B, D2B}
    double eri = 0.0;
    double rhoA;
    double rhoB;
    int lA;
    int lB;
    int mA;
    int mB;
    chg1 = NonZeroMultipoleMom_sp(L1,M1,L2,M2);
    chg2 = NonZeroMultipoleMom_sp(L3,M3,L4,M4);
    for (size_t mp1 = 0; mp1 < chg1.size(); ++mp1) {
      lA = chg1[mp1][0];
      mA = chg1[mp1][1];
      for (size_t mp2 = 0; mp2 < chg2.size(); ++mp2) {
        lB = chg2[mp2][0];
        mB = chg2[mp2][1];
        if (((abs(mA) % 2 == 0)&&(abs(mB) % 2 == 0))||((abs(mA) % 2 != 0)&&(abs(mB) % 2 != 0))) {
          //even-even or odd-odd only
          if (abs(mA) != abs(mB)) {
            mA = abs(mA);
            mB = abs(mB);
          }
          rhoA = rho(atmnrA,lA);
          rhoB = rho(atmnrB,lB);
          if (core > 0) {rhoB = rhocore(atmnrB);}                        //in this case use the correct core rho for the ss atom
          if ((core > 1)||(core < 0)) {rhoA = rhocore(atmnrA);}          //in this case use the correct core rho for the ss atom
          if (mA != 3) {eri += se_multipole_dR(lA,mA,lB,mB,RAB,rhoA + rhoB,D);}
          else {eri += 0.5*(se_multipole_dR(2,2,2,2,RAB,rhoA + rhoB,D) - se_multipole_dR(2,2,2,-2,RAB,rhoA + rhoB,D));}
        }
      }
    }
    return eri;
  }
  double eri2Center_dR2(int L1, int M1, int L2, int M2, int L3, int M3, int L4, int M4, double RAB, std::vector<double> D,size_t atmnrA,size_t atmnrB, int core = 0) {
    //function that calculates the second-derivatives of 2 electron repulsion integrals as in the NDDO approximation as applied in MNDO
    //AOs contains AO information on the 4 orbitals involved in the integral, meaning that this is only for s,p basis sets/atoms
    //RAB is distance between centers
    //D is a vector with all necessary D parameters: {D1A, D2A, D1B, D2B}
    double eri = 0.0;
    double rhoA;
    double rhoB;
    int lA;
    int lB;
    int mA;
    int mB;
    chg1 = NonZeroMultipoleMom_sp(L1,M1,L2,M2);
    chg2 = NonZeroMultipoleMom_sp(L3,M3,L4,M4);
    for (size_t mp1 = 0; mp1 < chg1.size(); ++mp1) {
      lA = chg1[mp1][0];
      mA = chg1[mp1][1];
      for (size_t mp2 = 0; mp2 < chg2.size(); ++mp2) {
        lB = chg2[mp2][0];
        mB = chg2[mp2][1];
        if (((abs(mA) % 2 == 0)&&(abs(mB) % 2 == 0))||((abs(mA) % 2 != 0)&&(abs(mB) % 2 != 0))) {
          //even-even or odd-odd only
          if (abs(mA) != abs(mB)) {
            mA = abs(mA);
            mB = abs(mB);
          }
          rhoA = rho(atmnrA,lA);
          rhoB = rho(atmnrB,lB);
          if (core > 0) {rhoB = rhocore(atmnrB);}                        //in this case use the correct core rho for the ss atom
          if ((core > 1)||(core < 0)) {rhoA = rhocore(atmnrA);}          //in this case use the correct core rho for the ss atom
          if (mA != 3) {eri += se_multipole_dR2(lA,mA,lB,mB,RAB,rhoA + rhoB,D);}
          else {eri += 0.5*(se_multipole_dR2(2,2,2,2,RAB,rhoA + rhoB,D) - se_multipole_dR2(2,2,2,-2,RAB,rhoA + rhoB,D));}
        }
      }
    }
    return eri;
  }
  virtual double d_eri2Center(int L1, int M1, int L2, int M2, int L3, int M3, int L4, int M4, double RAB, std::vector<double> D, size_t atmnrA, size_t atmnrB, int core = 0) {
    return 0.0;
  }
  virtual double d_eri2Center_dR(int L1, int M1, int L2, int M2, int L3, int M3, int L4, int M4, double RAB, std::vector<double> D, size_t atmnrA, size_t atmnrB, int core = 0) {
    return 0.0;
  }
  virtual double d_eri2Center_dR2(int L1, int M1, int L2, int M2, int L3, int M3, int L4, int M4, double RAB, std::vector<double> D, size_t atmnrA, size_t atmnrB, int core = 0) {
    return 0.0;
  }
  //other auxiliary functions
  virtual void checkAtoms() {
    //function that checks whether the atoms in the molecule are fully defined in the respective MNDO theories
    bool def;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      def = false;
      if (atoms[iatm] == 1) {def = true;}                                                       //H
      else if ((atoms[iatm] > 2)&&(atoms[iatm] < 10)) {def = true;}                             //Li,Be,B,C,N,O,F
      else if ((atoms[iatm] > 12)&&(atoms[iatm] < 18)) {def = true;}                            //Al,Si,P,S,Cl
      else if ((atoms[iatm] == 30)||(atoms[iatm] == 32)) {def = true;}                          //Zn,Ge
      else if ((atoms[iatm] == 35)||(atoms[iatm] == 53)) {def = true;}                          //Br,I
      else if ((atoms[iatm] == 50)||(atoms[iatm] == 80)||(atoms[iatm] == 82)) {def = true;}     //Sn,Hg,Pb
      if (!def) {throw("ERROR: MNDO.hpp: MNDO: checkAtoms(): atom not fully specified for MNDO-theory");}
    }
  }
  void OpenClosed() {
    //function to determine whether system is open- or closed-shell; also gets the number of core electrons per atom
    shell = "closed";
    ncoreelectrons = 0;
    if (mol.Multiplicity() > 1) {shell = "open";}
    else if (mol.Nelectrons()%2 != 0) {shell = "open";}
    if (CoreCharge.size() != Natoms) {CoreCharge.resize(Natoms);}
    for (size_t idx = 0; idx < Natoms; ++idx) {
      if (atoms[idx] < 3) {CoreCharge[idx] = 0;}                                //H,He
      else if ((atoms[idx] > 2)&&(atoms[idx] < 11))  {CoreCharge[idx] = 2;}     //Li-Ne
      else if ((atoms[idx] > 10)&&(atoms[idx] < 19)) {CoreCharge[idx] = 10;}    //Na-Ar
      else if ((atoms[idx] > 18)&&(atoms[idx] < 30)) {CoreCharge[idx] = 18;}    //K-Cu
      else if ((atoms[idx] > 29)&&(atoms[idx] < 37)) {CoreCharge[idx] = 28;}    //Zn-Kr
      else if ((atoms[idx] > 36)&&(atoms[idx] < 48)) {CoreCharge[idx] = 36;}    //Rb-Ag
      else if ((atoms[idx] > 47)&&(atoms[idx] < 55)) {CoreCharge[idx] = 46;}    //Cd-Xe
      else if ((atoms[idx] > 54)&&(atoms[idx] < 58)) {CoreCharge[idx] = 54;}    //Cs-La
      else if ((atoms[idx] > 70)&&(atoms[idx] < 80)) {CoreCharge[idx] = 68;}    //Lu-Au
      else if ((atoms[idx] > 79)&&(atoms[idx] < 84)) {CoreCharge[idx] = 78;}    //Hg-Bi
      if ((atoms[idx] == 10)||(atoms[idx] == 18)||(atoms[idx] == 36)||(atoms[idx] == 54)) {CoreCharge[idx] += 2;}
      ncoreelectrons += CoreCharge[idx];
    }
    if (print > 0) {
      std::cout << "core charge" << std::endl;
      for (size_t idx = 0; idx < CoreCharge.size(); ++idx) {
        std::cout << CoreCharge[idx] << " ";
      }
      std::cout << std::endl;
    }
  }
  void DampOscillations(int iteration, matrixE & DM, matrixE & oDM, std::vector<double> & pvec) {
    //apply mopac's oscillation damping, i.e., check whether the density matrix is changing by more than a certain value between iterations
    natocc = 1.0 + double(shell != "open");
    factorA = 0.0;
    factorB = 0.0;
    factor = 0.0;
    dampfactor = 10.0;
    if (iteration > 3) {dampfactor = oscillationdampthresh;}
    oscillationdamp = ((iteration + 1)%3 != 0);
    sumA = 0.0;
    sumB = 0.0;
    for (size_t iorb = 0; iorb < NAOs; ++iorb) {
      dval = DM(iorb + 1,iorb + 1);
      sumA += dval;
      deltaA = fabs(dval - oDM(iorb + 1,iorb + 1));
      if (!oscillationdamp) {
        factorA += deltaA*deltaA;
        deltaA = (dval - 2.0*oDM(iorb + 1,iorb + 1) + pvec[iorb]);
        factorB += deltaA*deltaA;
      }
      pvec[iorb] = oDM(iorb + 1,iorb + 1);
      oDM(iorb + 1,iorb + 1) = dval;
    }
    if (factorB >= 0.0) {
      if (factorA < 100.0*factorB) {factor = sqrt(factorA/factorB);}
    }
    for (size_t iorb = 0; iorb < NAOs; ++iorb) {
      for (size_t jorb = 0; jorb < iorb; ++jorb) {
        dval = DM(iorb + 1,jorb + 1);
        odval = oDM(iorb + 1,jorb + 1);
        oDM(iorb + 1,jorb + 1) = dval + factor*(dval - odval);
        DM(iorb + 1,jorb + 1) = oDM(iorb + 1,jorb + 1);
        DM(jorb + 1,iorb + 1) = oDM(iorb + 1,jorb + 1);
      }
      if (fabs(oDM(iorb + 1,iorb + 1) - pvec[iorb]) > dampfactor) {
        odval = oDM(iorb + 1,iorb + 1) - pvec[iorb];
        if (odval != 0.0) {
          dval = fabs(odval);
          odval /= dval;
        }
        oDM(iorb + 1,iorb + 1) = pvec[iorb] + odval*dampfactor;
      }
      else {oDM(iorb + 1,iorb + 1) += factor*(oDM(iorb + 1,iorb + 1) - pvec[iorb]);}
      odval = oDM(iorb + 1,iorb + 1);
      oDM(iorb + 1,iorb + 1) = fmin(natocc,fmax(0.0,odval));
      sumB += oDM(iorb + 1,iorb + 1);
      DM(iorb + 1,iorb + 1) = oDM(iorb + 1,iorb + 1);
    }
    deltaA = sumA;
label:
    dval = 0.0;
    if (sumB > 1.0e-3) {dval = sumA/sumB;}
    sumA = deltaA;
    if ((sumB > 1.0e-3)&&(fabs(dval - 1.0) > 1.0e-5)) {
      sumB = 0.0;
      for (size_t iorb = 0; iorb < NAOs; ++iorb) {
        oDM(iorb + 1,iorb + 1) *= dval;
        odval = oDM(iorb + 1,iorb + 1) + 1.0e-20;
        oDM(iorb + 1,iorb + 1) = fmax(odval,0.0);
        if (oDM(iorb + 1,iorb + 1) > natocc) {
          oDM(iorb + 1,iorb + 1) = natocc;
          sumA -= natocc;
        }
        else {sumB += oDM(iorb + 1,iorb + 1);}
        DM(iorb + 1,iorb + 1) = oDM(iorb + 1,iorb + 1);
      }
      goto label;
    }
    oscillationdamp = true;
  }
  //charge model stuff
  std::vector<double> CM1charge() {
    //function returning the CM1 partial charges
    //J. W. Storer, D. J. Giesen, C. J. Cramer, D. G. Truhlar, J. Comp.-Aid. Mol. Des., 9, 87, 1995
    //get Mulliken charges; note that actual population analysis not required here
    AOs = basis.AtomNAOs(atoms);
    this->getDens(dens);
    if (this->Shell() == "open") {this->getbDens(bdens);}
    std::vector<double> q0k;
    std::vector<double> Dqk(Natoms,0.0);
    std::vector<double> qk(Natoms,0.0);
    std::vector<double> Bkkp(Natoms,0.0);
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      qk[idatm] = double(CoreCharge[idatm]);
    }
    LMcharges(q0k,dens,AOs,qk);
    matrixE BKKp(Natoms,Natoms);
    ArmstrongBondOrder(BKKp);
    for (size_t irow = 0; irow < Natoms; ++irow) {
      qk[irow] = q0k[irow];
      Dqk[irow] += CM1ck(atoms[irow])*q0k[irow] + CM1dk(atoms[irow]);
      for (size_t icol = 0; icol < Natoms; ++icol) {
        if (irow == icol) {continue;}
        Bkkp[irow] += BKKp(irow + 1,icol + 1);
        if (atoms[irow] == 1) {Dqk[irow] += BKKp(irow + 1,icol + 1)*CM1dkkp(atoms[irow],atoms[icol]);}
        else if (atoms[irow] == 7) {
          if (atoms[icol] == 6) {Dqk[irow] += 0.5*(CM1ckkp(atoms[irow],atoms[icol])*q0k[irow] + CM1dkkp(atoms[irow],atoms[icol]))*(tanh(10.0*(BKKp(irow + 1,icol + 1) - 2.3)) + 1.0);}
          else if (atoms[icol] == 8) {Dqk[irow] += BKKp(irow + 1,icol + 1)*CM1dkkp(atoms[irow],atoms[icol]);}
        }
        else if (atoms[irow] == 8) {
          if (atoms[icol] == 16) {Dqk[irow] += BKKp(irow + 1,icol + 1)*CM1dkkp(atoms[irow],atoms[icol]);}
        }
      }
    }
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      qk[iatm] += Bkkp[iatm]*Dqk[iatm];
      for (size_t ibtm = 0; ibtm < Natoms; ++ibtm) {
        if (iatm == ibtm) {continue;}
        qk[iatm] -= BKKp(iatm + 1,ibtm + 1)*Dqk[ibtm];
      }
    }
    return qk;
  }
  std::vector<double> CM2charge() {
    //function returning the CM2 partial charges
    //J. Li, T. Zhu, C. J. Cramer, D. G. Truhlar, J. Phys. Chem. A, 102, 1820, 1998
    AOs = basis.AtomNAOs(atoms);
    this->getDens(dens);
    if (this->Shell() == "open") {this->getbDens(bdens);}
    sao = Identity(NAOs);                                   //for NDDO methods the overlap is supposed to be the identity matrix
    std::vector<double> q0k(Natoms,0.0);
    std::vector<double> qk(Natoms,0.0);
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      q0k[idatm] = double(CoreCharge[idatm]);
    }
    LMcharges(qk,dens,AOs,q0k);
    matrixE BKKp(Natoms,Natoms);
    MayerBondOrder(BKKp);
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      for (size_t ibtm = 0; ibtm < Natoms; ++ibtm) {
        if (iatm == ibtm) {continue;}
        qk[iatm] += BKKp(iatm + 1,ibtm + 1)*(CM2dkkp(atoms[iatm],atoms[ibtm]) + CM2ckkp(atoms[iatm],atoms[ibtm])*BKKp(iatm + 1,ibtm + 1));
      }
    }
    return qk;
  }
  std::vector<double> CM3charge() {
    //function returning the CM3 partial charges
    //P. Winget, J. D. Thompson, J. D. Xidos, C. J. Cramer, D. G. Truhlar, J. Phys. Chem. A, 106(44), 10707 2002
    //J. D. Thompson, C. J. Cramer, D. G. Truhlar, J. Comp. Chem., 24(11), 1291, 2003
    int other;
    AOs = basis.AtomNAOs(atoms);
    this->getDens(dens);
    if (this->Shell() == "open") {this->getbDens(bdens);}
    sao = Identity(NAOs);                                   //for NDDO methods the overlap is supposed to be the identity matrix
    std::vector<double> q0k(Natoms,0.0);
    std::vector<double> qk(Natoms,0.0);
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      q0k[idatm] = double(CoreCharge[idatm]);
    }
    LMcharges(qk,dens,AOs,q0k);
    matrixE BKKp(Natoms,Natoms);
    MayerBondOrder(BKKp);
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      if ((atoms[iatm] != 7)&&(atoms[iatm] != 8)) {
        for (size_t ibtm = 0; ibtm < Natoms; ++ibtm) {
          if (iatm == ibtm) {continue;}
          qk[iatm] += BKKp(iatm + 1,ibtm + 1)*(CM3dkkp(atoms[iatm],atoms[ibtm]) + CM3ckkp(atoms[iatm],atoms[ibtm])*BKKp(iatm + 1,ibtm + 1));
        }
      }
      else {
        if (atoms[iatm] == 8) {other = 7;}
        else {other = 8;}
        for (size_t ibtm = 0; ibtm < Natoms; ++ibtm) {
          if (iatm == ibtm) {continue;}
          if ((atoms[iatm] != 7)||(atoms[iatm] != 8)) {
            qk[iatm] += BKKp(iatm + 1,ibtm + 1)*(CM3dkkp(atoms[iatm],atoms[ibtm]) + CM3ckkp(atoms[iatm],atoms[ibtm])*BKKp(iatm + 1,ibtm + 1));}
          if (atoms[iatm] == other) {
            qk[iatm] += BKKp(iatm + 1,ibtm + 1)*(CM3dkkp(atoms[iatm],atoms[ibtm]) + CM3ckkp(atoms[iatm],atoms[ibtm])*exp(-(BKKp(iatm + 1,ibtm + 1)*BKKp(iatm + 1,ibtm + 1))/(B0()*B0())));}
        }
      }
    }
    return qk;
  }
  //estimate parameters for heat of formation
  virtual double calcEISOL(int atmnr) {
    //function that calculates the theoretical heats of formation for atoms
    //this function is not used directly by the methods, its purpose is the generation of the parameters that are then tabulated
    //integrals
    double Uss = UlX(atmnr,0);
    double Upp = UlX(atmnr,1);
    double Udd = UlX(atmnr,2);
    double Gss = eri1Center(atmnr,0,0);
    double Gpp = eri1Center(atmnr,2,2);
    double Gsp = eri1Center(atmnr,0,2);
    double Gp2 = eri1Center(atmnr,2,-2);
    double Hsp = eri1Center(atmnr,1,1);
    if (Hsp < 1.0e-7) {Hsp = 1.0e-7;}
    //Take into account constraints on the values of the integrals
    double Hpp = 0.5*(Gpp - Gp2);
    if (Hpp < 0.1) {Hpp = 0.1;}
    //counting how often integrals show up
    int nrs = NrS(atmnr);
    int nrp = NrP(atmnr);
    int nrd = NrD(atmnr);
    int ngss = (nrs - 1)*int(nrs > 1);                     //number of integrals of type <ss|ss>
    int ngsp = nrs*nrp;                                    //number of integrals of type <ss|pp>
    int iaux = std::min(nrp,6 - nrp);
    int npp = -iaux*(iaux - 1)/4;
    int ngp2 = nrp*(nrp - 1)/2 - npp;                      //number of integrals of type <pp|pp>; note that hpp is replaced by 0.5*(gpp - gp2) to insure rotational invariance
    int nhsp = -nrp;                                       //number of integrals of type <sp|sp>; note that if nrp != 0 then nrs = 2 due to Aufbau principle
    double eisol = Uss*nrs + Upp*nrp + Udd*nrd + Gss*ngss + Gpp*npp + Gsp*ngsp + Gp2*ngp2 + Hsp*nhsp;
    return eisol;
  }
  //parameters
  virtual double B0() {return 1.0;}
  virtual double ZeroOverlap(size_t atm) {return ZeroOverlapMNDO(atm);}
  virtual double ElementFormTheo(size_t atomicnr) {
    //function returning theoretical heats of formation for elements; values directly in a.u.
    double enth = 0.0;
    switch (atomicnr) {
      case 1:            //H
        enth = -0.437547617217;
        break;
      case 3:            //Li
        enth = -0.188450543318;
        break;
      case 4:            //Be
        enth = -0.889508467058;
        break;
      case 5:            //B
        enth = -2.3635678084;
        break;
      case 6:            //C
        enth = -4.40553146138;
        break;
      case 7:            //N
        enth = -7.41863055941;
        break;
      case 8:            //O
        enth = -11.6531563457;
        break;
      case 9:            //F
        enth = -17.5178076287;
        break;
      case 13:           //Al
        enth = -1.63475966018;
        break;
      case 14:           //Si
        enth = -3.03014442762;
        break;
      case 15:           //P
        enth = -5.60353574631;
        break;
      case 16:           //S
        enth = -8.2861420927;
        break;
      case 17:           //Cl
        enth = -12.9768362327;
        break;
      case 30:           //Zn
        enth = -1.09804898487;
        break;
      case 32:           //Ge
        enth = -2.78739756374;
        break;
      case 35:           //Br
        enth = -12.7403019342;
        break;
      case 50:           //Sn
        enth = -3.3781487843;
        break;
      case 53:           //I
        enth = -12.5167626458;
        break;
      case 80:           //Hg
        enth = -1.05908426259;
        break;
      case 82:           //Pb
        enth = -3.87464695071;
        break;
    }
    return enth;
  }
  virtual double betaA0(size_t atomicnr, int L) {
    //function returning the bonding parameters; values stored in eV, returned however in a.u.
    double betaa0 = 0.0;
    switch (atomicnr) {
      case 1:                   //H
        betaa0 = -6.989064;  
        break;
      case 3:                   //Li
        betaa0 = -1.3500400;
        break;
      case 4:                   //Be
        betaa0 = -4.0170960;
        break;
      case 5:                   //B
        betaa0 = -8.252054;
        break;
      case 6:                   //C
        if (L == 0) {betaa0 = -18.985044;}
        else if (L == 1) {betaa0 = -7.934122;}
        break;
      case 7:                   //N
        betaa0 = -20.495758;
        break;
      case 8:                   //O
        betaa0 = -32.688082;
        break;
      case 9:                   //F
        if (L == 0) {betaa0 = -48.290460;}
        else if (L == 1) {betaa0 = -36.508540;}
        break;
      case 13:                  //Al
        betaa0 = -2.6702840;
        break;
      case 14:                  //Si
        if (L == 0) {betaa0 = -9.0868040;}
        else if (L == 1) {betaa0 = -1.0758270;}
        break;
      case 15:                  //P
        betaa0 = -6.7916000;
        break;
      case 16:                  //S
        if (L == 0) {betaa0 = -10.7616700;}
        else if (L == 1) {betaa0 = -10.1084330;}
        break;
      case 17:                  //Cl
        betaa0 = -14.2623200;
        break;
      case 30:                  //Zn
        if (L == 0) {betaa0 = -1.000;}
        else if (L == 1) {betaa0 = -2.000;}
        break;
      case 32:                  //Ge
        if (L == 0) {betaa0 = -4.5164790;}
        else if (L == 1) {betaa0 = -1.7555170;}
        break;
      case 35:                  //Br
        if (L == 0) {betaa0 = -8.91710680;}
        else if (L == 1) {betaa0 = -9.9437398;}
        break;
      case 50:                  //Sn
        if (L == 0) {betaa0 = -3.2351470;}
        else if (L == 1) {betaa0 = -4.2904160;}
        break;
      case 53:                  //I
        if (L == 0) {betaa0 = -7.4144514;}
        else if (L == 1) {betaa0 = -6.1967812;}
        break;
      case 80:                  //Hg
        if (L == 0) {betaa0 = -0.4045250;}
        else if (L == 1) {betaa0 = -6.2066830;}
        break;
      case 82:                  //Pb
        if (L == 0) {betaa0 = -8.0423870;}
        else if (L == 1) {betaa0 = -3.00000;}
        break;
    }
    return betaa0/au2eV;
  }
  virtual double UlX(size_t atomicnr, int L) {
    //function that returns the atomic potential U
    double ulx = 0.0;
    switch (atomicnr) {
      case 1:                 //H
        ulx = -11.906276;
        break;
      case 3:                 //Li
        if (L == 0) {ulx = -5.1280000;}
        else if (L == 1) {ulx = -2.7212000;}
        break;
      case 4:                 //Be
        if (L == 0) {ulx = -16.6023780;}
        else if (L == 1) {ulx = -10.7037710;}
        break;
      case 5:                 //B
        if (L == 0) {ulx = -34.547130;}
        else if (L == 1) {ulx = -23.121690;}
        break;
      case 6:                 //C
        if (L == 0) {ulx = -52.279745;}
        else if (L == 1) {ulx = -39.205558;}
        break;
      case 7:                 //N
        if (L == 0) {ulx = -71.932122;}
        else if (L == 1) {ulx = -57.172319;}
        break;
      case 8:                 //O
        if (L == 0) {ulx = -99.644309;}
        else if (L == 1) {ulx = -77.797472;}
        break;
      case 9:                 //F
        if (L == 0) {ulx = -131.071548;}
        else if (L == 1) {ulx = -105.782137;}
        break;
      case 13:                //Al
        if (L == 0) {ulx = -23.8070970;}
        else if (L == 1) {ulx = -17.5198780;}
        break;
      case 14:                //Si
        if (L == 0) {ulx = -37.0375330;}
        else if (L == 1) {ulx = -27.76967800;}
        break;
      case 15:                //P
        if (L == 0) {ulx = -56.1433600;}
        else if (L == 1) {ulx = -42.8510800;}
        break;
      case 16:                //S
        if (L == 0) {ulx = -72.2422810;}
        else if (L == 1) {ulx = -56.9732070;}
        break;
      case 17:                //Cl
        if (L == 0) {ulx = -100.2271660;}
        else if (L == 1) {ulx = -77.3786670;}
        break;
      case 30:                //Zn
        if (L == 0) {ulx = -20.8397160;}
        else if (L == 1) {ulx = -19.6252240;}
        break;
      case 32:                //Ge
        if (L == 0) {ulx = -33.9493670;}
        else if (L == 1) {ulx = -27.4251050;}
        break;
      case 35:                //Br
        if (L == 0) {ulx = -99.98644054;}
        else if (L == 1) {ulx = -75.67130754;}
        break;
      case 50:                //Sn
        if (L == 0) {ulx = -40.8518020;}
        else if (L == 1) {ulx = -28.5602490;}
        break;
      case 53:                //I
        if (L == 0) {ulx = -100.0030538;}
        else if (L == 1) {ulx = -74.61146919;}
        break;
      case 80:                //Hg
        if (L == 0) {ulx = -19.8095740;}
        else if (L == 1) {ulx = -13.1025300;}
        break;
      case 82:                //Pb
        if (L == 0) {ulx = -47.3196920;}
        else if (L == 1) {ulx = -28.8475600;}
        break;
    }
    return ulx/au2eV;
  }
  virtual double alpha(int atomicnr, int atm2) {
    //function that returns alpha values for MNDO
    //alphas stored and returned in 1/Angstrom
    double alpha = 0.0;
    switch (atomicnr) {
      case 1:     //H
        alpha = 2.544134;
        break;
      case 3:     //Li
        alpha = 1.250140;
        break;
      case 4:     //Be
        alpha = 1.669434;
        break;
      case 5:     //B
        alpha = 2.134993;
        break;
      case 6:     //C
        alpha = 2.546380;
        break;
      case 7:     //N
        alpha = 2.861342;
        break;
      case 8:     //O
        alpha = 3.160604;
        break;
      case 9:     //F
        alpha = 3.419661;
        break;
      case 13:    //Al
        alpha = 1.868839;
        break;
      case 14:    //Si
        alpha = 2.205316;
        break;
      case 15:    //P
        alpha = 2.415280;
        break;
      case 16:    //S
        alpha = 2.478026;
        break;
      case 17:    //Cl
        alpha = 2.542201;
        break;
      case 30:    //Zn
        alpha = 1.506457;
        break;
      case 32:    //Ge
        alpha = 1.978498;
        break;
      case 35:    //Br
        alpha = 2.44570512;
        break;
      case 50:    //Sn
        alpha = 1.800814;
        break;
      case 53:    //I
        alpha = 2.20732001;
        break;
      case 80:    //Hg
        alpha = 1.335641;
        break;
      case 82:    //Pb
        alpha = 1.728333;
        break;
    }
    return alpha;
  }
  virtual double Dvalue(size_t atmnr, size_t idx) {
    //function returning the D values needed to calculate eris; values stored in Angstrom, but returned in a.u.
    //note that this function returns both D1 and D2; idx is then either 1 or 2
    double dval = 0.0;
    switch (atmnr) {
      case 1:                    //H
        dval = 0.0;
        break;
      case 3:                    //Li
        if (idx == 1) {dval = 1.0874267;}
        else if (idx == 2) {dval = 0.9227121;}
        break;
      case 4:                    //Be
        if (idx == 1) {dval = 0.7605847;}
        else if (idx == 2) {dval = 0.6453775;}
        break;
      case 5:                    //B
        if (idx == 1) {dval = 0.5068929;}
        else if (idx == 2) {dval = 0.4301129;}
        break;
      case 6:                    //C
        if (idx == 1) {dval = 0.4272845;}
        else if (idx == 2) {dval = 0.3625629;}
        break;
      case 7:                    //N
        if (idx == 1) {dval = 0.3386159;}
        else if (idx == 2) {dval = 0.2873251;}
        break;
      case 8:                    //O
        if (idx == 1) {dval = 0.2828939;}
        else if (idx == 2) {dval = 0.2400435;}
        break;
      case 9:                    //F
        if (idx == 1) {dval = 0.2681377;}
        else if (idx == 2) {dval = 0.2275224;}
        break;
      case 13:                   //Al
        if (idx == 1) {dval = 0.7404309;}
        else if (idx == 2) {dval = 0.6131351;}
        break;
      case 14:                   //Si
        if (idx == 1) {dval = 0.6657106;}
        else if (idx == 2) {dval = 0.5178335;}
        break;
      case 15:                   //P
        if (idx == 1) {dval = 0.5360302;}
        else if (idx == 2) {dval = 0.4958342;}
        break;
      case 16:                   //S
        if (idx == 1) {dval = 0.4863010;}
        else if (idx == 2) {dval = 0.4407175;}
        break;
      case 17:                   //Cl
        if (idx == 1) {dval = 0.2638887;}
        else if (idx == 2) {dval = 0.4348484;}
        break;
      case 30:                   //Zn
        if (idx == 1) {dval = 0.6899187;}
        else if (idx == 2) {dval = 0.7683602;}
        break;
      case 32:                   //Ge
        if (idx == 1) {dval = 0.6644269;}
        else if (idx == 2) {dval = 0.5555542;}
        break;
      case 35:                   //Br
        if (idx == 1) {dval = 0.3202029;}
        else if (idx == 2) {dval = 0.5104278;}
        break;
      case 50:                   //Sn
        if (idx == 1) {dval = 0.8306740;}
        else if (idx == 2) {dval = 0.7017967;}
        break;
      case 53:                   //I
        if (idx == 1) {dval = 0.7542341;}
        else if (idx == 2) {dval = 0.6266241;}
        break;
      case 80:                   //Hg
        if (idx == 1) {dval = 0.9195890;}
        else if (idx == 2) {dval = 0.7730105;}
        break;
      case 82:                   //Pb
        if (idx == 1) {dval = 0.8216177;}
        else if (idx == 2) {dval = 0.7666867;}
        break;
    }
    return dval*dist_Angstrom2aum1;
  }
  virtual double rhocore(size_t atmnr) {
    //function returning rho values used to calculate VAB integrals; values already in atomic units
    double rho_ = rho(atmnr,0);
    return rho_;
  }
  virtual double rho(size_t atmnr, size_t l) {
    //function returning the rho values needed to calculate eris; values stored in Angstrom but returned in atomic units
    double rho = 0.0;
    switch (atmnr) {
      case 1:                    //H
        if (l == 0) {rho = 0.560345403;}
        else if (l == 1) {rho = 0.0;}
        else if (l == 2) {rho = 0.0;}
        break;
      case 3:                    //Li
        if (l == 0) {rho = 0.986207884;}
        else if (l == 1) {rho = 1.165672376;}
        else if (l == 2) {rho = 1.011953732;}
        break;
      case 4:                    //Be
        if (l == 0) {rho = 0.799924235;}
        else if (l == 1) {rho = 0.788356095;}
        else if (l == 2) {rho = 0.687877905;}
        break;
      case 5:                    //B
        if (l == 0) {rho = 0.679822279;}
        else if (l == 1) {rho = 0.539445596;}
        else if (l == 2) {rho = 0.476128306;}
        break;
      case 6:                    //C
        if (l == 0) {rho = 0.588660438;}
        else if (l == 1) {rho = 0.430253872;}
        else if (l == 2) {rho = 0.395733736;}
        break;
      case 7:                    //N
        if (l == 0) {rho = 0.529751104;}
        else if (l == 1) {rho = 0.337322211;}
        else if (l == 2) {rho = 0.325583327;}
        break;
      case 8:                    //O
        if (l == 0) {rho = 0.466881794;}
        else if (l == 1) {rho = 0.275821517;}
        else if (l == 2) {rho = 0.2786282;}
        break;
      case 9:                    //F
        if (l == 0) {rho = 0.425491557;}
        else if (l == 1) {rho = 0.243848996;}
        else if (l == 2) {rho = 0.255793341;}
        break;
      case 13:                   //Al
        if (l == 0) {rho = 0.889903107;}
        else if (l == 1) {rho = 1.00389327;}
        else if (l == 2) {rho = 0.720237317;}
        break;
      case 14:                   //Si
        if (l == 0) {rho = 0.733128067;}
        else if (l == 1) {rho = 0.722068454;}
        else if (l == 2) {rho = 0.587084012;}
        break;
      case 15:                   //P
        if (l == 0) {rho = 0.622778301;}
        else if (l == 1) {rho = 0.541910569;}
        else if (l == 2) {rho = 0.531355547;}
        break;
      case 16:                   //S
        if (l == 0) {rho = 0.558953167;}
        else if (l == 1) {rho = 0.477199756;}
        else if (l == 2) {rho = 0.473718785;}
        break;
      case 17:                   //Cl
        if (l == 0) {rho = 0.478996434;}
        else if (l == 1) {rho = 0.328217689;}
        else if (l == 2) {rho = 0.435988;}
        break;
      case 30:                   //Zn
        if (l == 0) {rho = 0.610111605;}
        else if (l == 1) {rho = 1.113608164;}
        else if (l == 2) {rho = 0.966035844;}
        break;
      case 32:                   //Ge
        if (l == 0) {rho = 0.734624198;}
        else if (l == 1) {rho = 0.726135254;}
        else if (l == 2) {rho = 0.608610513;}
        break;
      case 35:                   //Br
        if (l == 0) {rho = 0.478791611;}
        else if (l == 1) {rho = 0.364523933;}
        else if (l == 2) {rho = 0.474624228;}
        break;
      case 50:                   //Sn
        if (l == 0) {rho = 0.734624198;}
        else if (l == 1) {rho = 0.821688399;}
        else if (l == 2) {rho = 0.712428177;}
        break;
      case 53:                   //I
        if (l == 0) {rho = 0.478664021;}
        else if (l == 1) {rho = 0.576001573;}
        else if (l == 2) {rho = 0.577015931;}
        break;
      case 80:                   //Hg
        if (l == 0) {rho = 0.666603429;}
        else if (l == 1) {rho = 0.868143258;}
        else if (l == 2) {rho = 0.759620304;}
        break;
      case 82:                   //Pb
        if (l == 0) {rho = 0.734624198;}
        else if (l == 1) {rho = 0.816789939;}
        else if (l == 2) {rho = 0.755508834;}
        break;
    }
    return rho*dist_Angstrom2aum1;
  }
  virtual double eri1Center(int atmnr, int Lbra, int Lket) {
    //function that gives back the semi-empirical 1-center eris
    //values stored in eV but returned in a.u.
    //Lbra is the sum of azimuthal quantum numbers for bra (ss = 0; pp = 2; sp = 1)
    //Lket is the sum of azimuthal quantum numbers for ket (ss = 0; pp = 2; sp = 1; p*p* = -2)
    double eri = 0.0;
    switch (atmnr) {
      case 1:            //H
        eri = 12.848;                                                                            //(ss|ss)
        break;
      case 3:           //Li
        if ((Lbra == 0)&&(Lket == 0)) {eri = 7.3;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.0;}                                          //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.42;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.52;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.83;}                                         //(sp|sp)||(ps|ps)
        break;
      case 4:          //Be
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.0;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.97;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.43;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.22;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.28;}                                         //(sp|sp)||(ps|ps)
        break;
      case 5:          //B
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.59;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.86;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 9.56;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.86;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.81;}                                         //(sp|sp)||(ps|ps)
        break;
      case 6:          //C
        if ((Lbra == 0)&&(Lket == 0)) {eri = 12.23;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.08;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 11.47;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.84;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.43;}                                         //(sp|sp)||(ps|ps)
        break;
      case 7:          //N
        if ((Lbra == 0)&&(Lket == 0)) {eri = 13.59;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 12.98;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 12.66;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 11.59;}        //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.14;}                                         //(sp|sp)||(ps|ps)
        break;
      case 8:          //O
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.42;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 14.52;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 14.48;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 12.98;}        //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.94;}                                         //(sp|sp)||(ps|ps)
        break;
      case 9:          //F
        if ((Lbra == 0)&&(Lket == 0)) {eri = 16.92;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 16.71;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 17.25;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 14.91;}        //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 4.83;}                                         //(sp|sp)||(ps|ps)
        break;
      case 13:         //Al
        if ((Lbra == 0)&&(Lket == 0)) {eri = 8.09;}                                              //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.98;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.63;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.40;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.70;}                                         //(sp|sp)||(ps|ps)
        break;
      case 14:         //Si
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.82;}                                              //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.31;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.36;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.54;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.32;}                                         //(sp|sp)||(ps|ps)
        break;
      case 15:         //P
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.56;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.64;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 10.08;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.68;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.92;}                                         //(sp|sp)||(ps|ps)
        break;
      case 16:         //S
        if ((Lbra == 0)&&(Lket == 0)) {eri = 12.88;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 9.90;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 11.26;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 8.83;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.26;}                                         //(sp|sp)||(ps|ps)
        break;
      case 17:         //Cl
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.03;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.30;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 13.16;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.97;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.42;}                                         //(sp|sp)||(ps|ps)
        break;
      case 30:         //Zn
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.8;}                                              //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 13.3;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 11.182018;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 12.93052;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.484606;}                                     //(sp|sp)||(ps|ps)
        break;
      case 32:         //Ge
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.8;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.3;}                                          //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.3;}            //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.5;}          //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.3;}                                          //(sp|sp)||(ps|ps)
        break;
      case 35:         //Br
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.03643948;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.27632539;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 13.03468242;}    //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.85442552;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.45586832;}                                   //(sp|sp)||(ps|ps)
        break;
      case 50:         //Sn
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.8;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.3;}                                          //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.3;}            //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.5;}          //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.3;}                                          //(sp|sp)||(ps|ps)
        break;
      case 53:         //I
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.04044855;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.14778369;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 13.05655798;}    //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.91409071;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.45638202;}                                   //(sp|sp)||(ps|ps)
        break;
      case 80:         //Hg
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.8;}                                              //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 14.3;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 9.3;}            //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 13.5;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.3;}                                          //(sp|sp)||(ps|ps)
        break;
      case 82:         //Pb
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.8;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.3;}                                          //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.3;}            //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.5;}          //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.3;}                                          //(sp|sp)||(ps|ps)
        break;
    }
    return eri/au2eV;
  }
  virtual double xAB(int atm1, int atm2) {
    //function returning xAB values for MNDOd, all zero
    return 0.0;
  }
  virtual double CM1ck(size_t atomicnr) {
    //function returning the atomic parameters ck for CM1
    return 0.0;
  }
  virtual double CM1dk(size_t atomicnr) {
    //function returning the atomic parameters dk for CM1
    return 0.0;
  }
  virtual double CM1ckkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters ckk' for CM1
    return 0.0;
  }
  virtual double CM1dkkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters dkk' for CM1
    return 0.0;
  }
  virtual double CM2ckkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters ckk' for CM2
    return 0.0;
  }
  virtual double CM2dkkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters dkk' for CM2
    return 0.0;
  }
  virtual double CM3ckkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters ckk' for CM3
    return 0.0;
  }
  virtual double CM3dkkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters dkk' for CM3
    return 0.0;
  }
  virtual double PAPDDG(int atmnr, int index) {return 0.0;}
  virtual double DAPDDG(int atmnr, int index) {return 0.0;}
  double ValenceElectrons(int atmA) {
    //function returning the number of valence electrons of a certain atom
    double valelectron = 0.0;
    switch (atmA) {
      case 1:      //H
        valelectron = 1.0;
        break;
      case 6:      //C
        valelectron = 4.0;
        break;
      case 7:      //N
        valelectron = 5.0;
        break;
      case 8:      //O
        valelectron = 6.0;
        break;
      case 9:      //F
        valelectron = 7.0;
        break;
      case 14:     //Si
        valelectron = 4.0;
        break;
      case 15:     //P
        valelectron = 5.0;
        break;
      case 16:     //S
        valelectron = 6.0;
        break;
      case 17:     //Cl
        valelectron = 7.0;
        break;
      case 35:     //Br
        valelectron = 7.0;
        break;
      case 53:     //I
        valelectron = 7.0;
        break;
    }
    return valelectron;
  }
};
class AM1: public MNDO {
  //this is the implementation of Dewar's AM1
  //M. J. S. Dewar, E. G. Zoebisch, E. F. Healy, J. J. P. Stewart, J. Am. Chem. Soc., 107(13), 3902, 1985
  //for charge and dipole moment model: J. W. Storer, D. J. Giesen, C. J. Cramer, D. G. Truhlar, J. Comp.-Aid. Mol. Des., 9, 87, 1995   -> CM1
  //                                    J. Li, T. Zhu, C. J. Cramer, D. G. Truhlar, J. Phys. Chem. A, 102, 1820, 1998                   -> CM2
  //bX quantities are "barred" tensors, which are used only for the open-shell case
public:
  AM1(BSet _bset, Molecule _mol, std::string _openclosed = "0", std::string _corecorrection = "0"): MNDO(_bset,_mol,_openclosed,_corecorrection) {}
  ~AM1() {}
  std::string Type() {return "AM1";}
  double gfactor(size_t iatm1, size_t iatm2, double RAB) {return 1.0;}
  virtual double gfactor_dR(size_t iatm1, size_t iatm2) {return 0.0;}
  double AM1factor(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {
    //RAB in a.u.
    double AM1factor = 0.0;
    RAB *= dist_Angstrom2au;
    for (size_t idx = 0; idx < K1.size(); ++idx) {
      AM1factor += K1[idx]*exp(-L1[idx]*(RAB - M1[idx])*(RAB - M1[idx]))/RAB;
    }
    return AM1factor;
  }
  virtual double AM1factor_dR(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {
    //RAB in a.u.
    double am1factor_dr = 0.0;
    double RR = RAB*dist_Angstrom2au;
    for (size_t idx = 0; idx < K1.size(); ++idx) {
      am1factor_dr += K1[idx]*(2.0*RR*RR*L1[idx] - 2.0*L1[idx]*M1[idx]*RR + 1.0)*exp(-L1[idx]*(RR - M1[idx])*(RR - M1[idx]))/(RAB*RR);
    }
    return am1factor_dr;  
  }
  virtual double AM1factor_dR2(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {
    //RAB in a.u.
    double am1factor_dr2 = 0.0;
    double RR = RAB*dist_Angstrom2au;
    for (size_t idx = 0; idx < K1.size(); ++idx) {
      am1factor_dr2 += 2.0*K1[idx]*exp(-L1[idx]*(RR - M1[idx])*(RR - M1[idx]))*(L1[idx]*M1[idx]*RR - 1.0 - RR*L1[idx]*(RR - M1[idx])*(2.0*RR*RR*L1[idx] - 2.0*L1[idx]*M1[idx]*RR + 1.0))/(RAB*RAB*RR);
    }
    return am1factor_dr2;  
  }
  //other auxiliary functions
  void checkAtoms() {
    //function that checks whether the atoms in the molecule are fully defined in the respective MNDO theories
    bool def;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      def = false;
      if (atoms[iatm] == 1) {def = true;}                                   //H
      else if ((atoms[iatm] > 2)&&(atoms[iatm] < 10)) {def = true;}         //Li,Be,B,C,N,O,F
      else if ((atoms[iatm] > 12)&&(atoms[iatm] < 18)) {def = true;}        //Al,Si,P,S,Cl
      else if ((atoms[iatm] == 30)||(atoms[iatm] == 32)) {def = true;}      //Zn,Ge
      else if ((atoms[iatm] == 33)||(atoms[iatm] == 34)) {def = true;}      //As,Se
      else if ((atoms[iatm] == 35)||(atoms[iatm] == 53)) {def = true;}      //Br,I
      else if ((atoms[iatm] == 51)||(atoms[iatm] == 52)) {def = true;}      //Sb,Te
      else if ((atoms[iatm] == 11)||(atoms[iatm] == 80)) {def = true;}      //Na,Hg
      if (!def) {throw("ERROR: MNDO.hpp: AM1: checkAtoms(): atom not fully specified for AM1-theory");}
    }
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapAM1(atm);}
  double ElementFormTheo(size_t atomicnr) {
    //function returning theoretical heats of formation for elements; values directly in a.u.
    double enth = 0.0;
    switch (atomicnr) {
      case 1:       //H
        enth = -0.418811010146;
        break;
      case 3:       //Li
        enth = -0.188450543318;
        break;
      case 4:       //Be
        enth = -0.889508467058;
        break;
      case 5:       //B
        enth = -2.34156653822;
        break;
      case 6:       //C
        enth = -4.41711440792;
        break;
      case 7:       //N
        enth = -7.41280733472;
        break;
      case 8:       //O
        enth = -11.5881473027;
        break;
      case 9:       //F
        enth = -17.7238538228;
        break;
      case 11:      //Na
        enth = -0.19313741269;
        break;
      case 13:      //Al
        enth = -1.70593365991;
        break;
      case 14:      //Si
        enth = -2.88911227453;
        break;
      case 15:      //P
        enth = -4.56249735046;
        break;
      case 16:      //S
        enth = -7.00454584332;
        break;
      case 17:      //Cl
        enth = -13.6780414477;
        break;
      case 30:      //Zn
        enth = -1.11277017685;
        break;
      case 32:      //Ge
        enth = -2.88509340499;
        break;
      case 33:      //As
        enth = -4.49926063551;
        break;
      case 34:      //Se
        enth = -5.25247215157;
        break;
      case 35:      //Br
        enth = -12.9473096825;
        break;
      case 51:      //Sb
        enth = -4.43488928513;
        break;
      case 52:      //Te
        enth = -5.12393240263;
        break;
      case 53:      //I
        enth = -12.7470286834;
        break;
      case 80:      //Hg
        enth = -1.06878637863;
        break;
    }
    return enth;
  }
  double betaA0(size_t atomicnr, int L) {
    //function returning the bonding parameters; values stored in eV, returned however in a.u.
    double betaa0 = 0.0;
    switch (atomicnr) {
      case 1:      //H
        betaa0 = -6.173787;
        break;
      case 3:      //Li
        betaa0 = -1.35004;
        break;
      case 4:      //Be
        betaa0 = -4.017096;
        break;
      case 5:      //B
        if (L == 0) {betaa0 = -9.599114;}
        else if (L == 1) {betaa0 = -6.273757;}
        break;
      case 6:      //C
        if (L == 0) {betaa0 = -15.715783;}
        else if (L == 1) {betaa0 = -7.719283;}
        break;
      case 7:      //N
        if (L == 0) {betaa0 = -20.299110;}
        else if (L == 1) {betaa0 = -18.238666;}
        break;
      case 8:      //O
        betaa0 = -29.272773;
        break;
      case 9:      //F
        if (L == 0) {betaa0 = -69.590277;}
        else if (L == 1) {betaa0 = -27.92236;}
        break;
      case 11:     //Na
        if (L == 0) {betaa0 = -1.4536944;}
        else if (L == 1) {betaa0 = -0.2298064;}
        break;
      case 13:     //Al
        if (L == 0) {betaa0 = -3.866822;}
        else if (L == 1) {betaa0 = -2.317146;}
        break;
      case 14:     //Si
        if (L == 0) {betaa0 = -3.784852;}
        else if (L == 1) {betaa0 = -1.968123;}
        break;
      case 15:     //P
        if (L == 0) {betaa0 = -6.353764;}
        else if (L == 1) {betaa0 = -6.590709;}
        break;
      case 16:     //S
        if (L == 0) {betaa0 = -3.920566;}
        else if (L == 1) {betaa0 = -7.905278;}
        break;
      case 17:     //Cl
        if (L == 0) {betaa0 = -24.59467;}
        else if (L == 1) {betaa0 = -14.637216;}
        break;
      case 30:     //Zn
        if (L == 0) {betaa0 = -1.997429;}
        else if (L == 1) {betaa0 = -4.758119;}
        break;
      case 32:     //Ge
        if (L == 0) {betaa0 = -4.356607;}
        else if (L == 1) {betaa0 = -0.991091;}
        break;
      case 33:     //As
        if (L == 0) {betaa0 = -5.6481504;}
        else if (L == 1) {betaa0 = -4.9979109;}
        break;
      case 34:     //Se
        if (L == 0) {betaa0 = -3.1470826;}
        else if (L == 1) {betaa0 = -6.1468406;}
        break;
      case 35:     //Br
        if (L == 0) {betaa0 = -19.39988;}
        else if (L == 1) {betaa0 = -8.957195;}
        break;
      case 51:     //Sb
        if (L == 0) {betaa0 = -7.38233;}
        else if (L == 1) {betaa0 = -3.633119;}
        break;
      case 52:     //Te
        if (L == 0) {betaa0 = -8.3897294;}
        else if (L == 1) {betaa0 = -5.1065429;}
        break;
      case 53:     //I
        if (L == 0) {betaa0 = -8.443327;}
        else if (L == 1) {betaa0 = -6.323405;}
        break;
      case 80:     //Hg
        if (L == 0) {betaa0 = -0.908657;}
        else if (L == 1) {betaa0 = -4.909384;}
        break;
    }
    return betaa0/au2eV;
  }
  double UlX(size_t atomicnr, int L) {
    //function that returns the atomic potential U
    double ulx = 0.0;
    switch (atomicnr) {
      case 1:                 //H
        ulx = -11.396427;
        break;
      case 3:                 //Li
        if (L == 0) {ulx = -5.128;}
        else if (L == 1) {ulx = -2.7212;}
        break;
      case 4:                 //Be
        if (L == 0) {ulx = -16.602378;}
        else if (L == 1) {ulx = -10.703771;}
        break;
      case 5:                 //B
        if (L == 0) {ulx = -34.492870;}
        else if (L == 1) {ulx = -22.631525;}
        break;
      case 6:                 //C
        if (L == 0) {ulx = -52.028658;}
        else if (L == 1) {ulx = -39.614239;}
        break;
      case 7:                 //N
        if (L == 0) {ulx = -71.860000;}
        else if (L == 1) {ulx = -57.167581;}
        break;
      case 8:                 //O
        if (L == 0) {ulx = -97.830000;}
        else if (L == 1) {ulx = -78.262380;}
        break;
      case 9:                 //F
        if (L == 0) {ulx = -136.105579;}
        else if (L == 1) {ulx = -104.889885;}
        break;
      case 11:                //Na
        if (L == 0) {ulx = -5.2555362;}
        else if (L == 1) {ulx = -2.0812781;}
        break;
      case 13:                //Al
        if (L == 0) {ulx = -24.353585;}
        else if (L == 1) {ulx = -18.363645;}
        break;
      case 14:                //Si
        if (L == 0) {ulx = -33.953622;}
        else if (L == 1) {ulx = -28.934749;}
        break;
      case 15:                //P
        if (L == 0) {ulx = -42.029863;}
        else if (L == 1) {ulx = -34.030709;}
        break;
      case 16:                //S
        if (L == 0) {ulx = -56.694056;}
        else if (L == 1) {ulx = -48.717049;}
        break;
      case 17:                //Cl
        if (L == 0) {ulx = -111.613948;}
        else if (L == 1) {ulx = -76.640107;}
        break;
      case 30:                //Zn
        if (L == 0) {ulx = -21.040008;}
        else if (L == 1) {ulx = -17.655574;}
        break;
      case 32:                //Ge
        if (L == 0) {ulx = -34.183889;}
        else if (L == 1) {ulx = -28.640811;}
        break;
      case 33:                //As
        if (L == 0) {ulx = -41.681751;}
        else if (L == 1) {ulx = -33.4506152;}
        break;
      case 34:                //Se
        if (L == 0) {ulx = -41.9984056;}
        else if (L == 1) {ulx = -32.8575485;}
        break;
      case 35:                //Br
        if (L == 0) {ulx = -104.656063;}
        else if (L == 1) {ulx = -74.930052;}
        break;
      case 51:                //Sb
        if (L == 0) {ulx = -44.438162;}
        else if (L == 1) {ulx = -32.389514;}
        break;
      case 52:                //Te
        if (L == 0) {ulx = -39.245423;}
        else if (L == 1) {ulx = -30.8515845;}
        break;
      case 53:                //I
        if (L == 0) {ulx = -103.589663;}
        else if (L == 1) {ulx = -74.429997;}
        break;
      case 80:                //Hg
        if (L == 0) {ulx = -19.941578;}
        else if (L == 1) {ulx = -11.11087;}
        break;
    }
    return ulx/au2eV;
  }
  double alpha(int atomicnr, int atm2) {
    //function that returns alpha values for AM1
    //alphas stored and returned in 1/Angstrom
    double alpha = 0.0;
    switch (atomicnr) {
      case 1:      //H
        alpha = 2.882324;
        break;
      case 3:      //Li
        alpha = 1.25014;
        break;
      case 4:      //Be
        alpha = 1.669434;
        break;
      case 5:      //B
        alpha = 2.446909;
        break;
      case 6:      //C
        alpha = 2.648274;
        break;
      case 7:      //N
        alpha = 2.947286;
        break;
      case 8:      //O
        alpha = 4.455371;
        break;
      case 9:      //F
        alpha = 5.5178;
        break;
      case 11:     //Na
        alpha = 2.2487164;
        break;
      case 13:     //Al
        alpha = 1.976586;
        break;
      case 14:     //Si
        alpha = 2.257816;
        break;
      case 15:     //P
        alpha = 2.455322;
        break;
      case 16:     //S
        alpha = 2.461648;
        break;
      case 17:     //Cl
        alpha = 2.919368;
        break;
      case 30:     //Zn
        alpha = 1.484563;
        break;
      case 32:     //Ge
        alpha = 2.136405;
        break;
      case 33:     //As
        alpha = 2.240538;
        break;
      case 34:     //Se
        alpha = 2.6375694;
        break;
      case 35:     //Br
        alpha = 2.576546;
        break;
      case 51:     //Sb
        alpha = 2.276331;
        break;
      case 52:     //Te          this value is totally unphysical!!!
        alpha = 6.0171167;
        break;
      case 53:     //I
        alpha = 2.299424;
        break;
      case 80:     //Hg
        alpha = 1.484734;
        break;
    }
    return alpha;
  }
  std::vector<double> AM1K(int atomicnr, int atm2 = 0) {
    std::vector<double> am1k;
    switch (atomicnr) {
      case 1:     //H
        am1k.push_back(0.122796/au2eV);
        am1k.push_back(0.005090/au2eV);
        am1k.push_back(-0.018336/au2eV);
        break;
      case 3:     //Li
        am1k.push_back(0.0);
        break;
      case 4:     //Be
        am1k.push_back(0.0);
        break;
      case 5:     //B, boron has different interactions with different atoms, as a unique set of parameters could not be fitted
        if (atm2 == 1) {                                                       //H
          am1k.push_back(0.412253/au2eV);
          am1k.push_back(-0.149917/au2eV);
        }
        else if (atm2 == 6) {                                                  //C
          am1k.push_back(0.261751/au2eV);
          am1k.push_back(0.050275/au2eV);
        }
        else if ((atm2 == 9)||(atm2 == 17)||(atm2 == 35)||(atm2 == 53)) {      //Halogen
          am1k.push_back(0.359244/au2eV);
          am1k.push_back(0.074729/au2eV);
        }
        else {                                                                 //anything else
          am1k.push_back(0.182613/au2eV);
          am1k.push_back(0.118587/au2eV);
          am1k.push_back(-0.073280/au2eV);
        }
        break;
      case 6:     //C
        am1k.push_back(0.011355/au2eV);
        am1k.push_back(0.045924/au2eV);
        am1k.push_back(-0.020061/au2eV);
        am1k.push_back(-0.001260/au2eV);
        break;
      case 7:     //N
        am1k.push_back(0.025251/au2eV);
        am1k.push_back(0.028953/au2eV);
        am1k.push_back(-0.005806/au2eV);
        break;
      case 8:     //O
        am1k.push_back(0.280962/au2eV);
        am1k.push_back(0.081430/au2eV);
        break;
      case 9:     //F
        am1k.push_back(0.242079/au2eV);
        am1k.push_back(0.003607/au2eV);
        break;
      case 11:    //Na
        am1k.push_back(0.5322668/au2eV);
        am1k.push_back(0.9223598/au2eV);
        break;
      case 13:    //Al
        am1k.push_back(0.09/au2eV);
        break;
      case 14:    //Si
        am1k.push_back(0.25/au2eV);
        am1k.push_back(0.061513/au2eV);
        am1k.push_back(0.020789/au2eV);
        break;
      case 15:    //P
        am1k.push_back(-0.031827/au2eV);
        am1k.push_back(0.01847/au2eV);
        am1k.push_back(0.03329/au2eV);
        break;
      case 16:    //S
        am1k.push_back(-0.509195/au2eV);
        am1k.push_back(-0.011863/au2eV);
        am1k.push_back(0.012334/au2eV);
        break;
      case 17:    //Cl
        am1k.push_back(0.094243/au2eV);
        am1k.push_back(0.027168/au2eV);
        break;
      case 30:    //Zn
        am1k.push_back(0.0);
        break;
      case 32:    //Ge
        am1k.push_back(0.0);
        break;
      case 33:    //As
        am1k.push_back(-0.0073614/au2eV);
        am1k.push_back(0.0437629/au2eV);
        break;
      case 34:    //Se
        am1k.push_back(0.1116681/au2eV);
        am1k.push_back(0.0396143/au2eV);
        break;
      case 35:    //Br
        am1k.push_back(0.066685/au2eV);
        am1k.push_back(0.025568/au2eV);
        break;
      case 51:    //Sb
        am1k.push_back(-0.596447/au2eV);
        am1k.push_back(0.895513/au2eV);
        break;
      case 52:    //Te
        am1k.push_back(0.4873378/au2eV);
        am1k.push_back(0.1520464/au2eV);
        break;
      case 53:    //I
        am1k.push_back(0.004361/au2eV);
        am1k.push_back(0.015706/au2eV);
        break;
      case 80:    //Ge
        am1k.push_back(0.0);
        break;
    }
    return am1k;
  }
  std::vector<double> AM1L(int atomicnr, int atm2 = 0) {
    std::vector<double> am1l;
    switch (atomicnr) {
      case 1:     //H
        am1l.push_back(5.0);
        am1l.push_back(5.0);
        am1l.push_back(2.0);
        break;
      case 3:     //Li
        am1l.push_back(0.0);
        break;
      case 4:     //Be
        am1l.push_back(0.0);
        break;
      case 5:     //B
        if (atm2 == 1) {
          am1l.push_back(10.0);
          am1l.push_back(6.0);
        }
        else if (atm2 == 6) {
          am1l.push_back(8.0);
          am1l.push_back(5.0);
        }
        else if ((atm2 == 9)||(atm2 == 17)||(atm2 == 35)||(atm2 == 53)) {
          am1l.push_back(9.0);
          am1l.push_back(9.0);
        }
        else {
          am1l.push_back(6.0);
          am1l.push_back(6.0);
          am1l.push_back(5.0);
        }
        break;
      case 6:     //C
        am1l.push_back(5.0);
        am1l.push_back(5.0);
        am1l.push_back(5.0);
        am1l.push_back(5.0);
        break;
      case 7:     //N
        am1l.push_back(5.0);
        am1l.push_back(5.0);
        am1l.push_back(2.0);
        break;
      case 8:     //O
        am1l.push_back(5.0);
        am1l.push_back(7.0);
        break;
      case 9:     //F
        am1l.push_back(4.8);
        am1l.push_back(4.6);
        break;
      case 11:    //Na
        am1l.push_back(0.4800304);
        am1l.push_back(1.9076776);
        break;
      case 13:    //Al
        am1l.push_back(12.392443);
        break;
      case 14:    //Si
        am1l.push_back(9.0);
        am1l.push_back(5.0);
        am1l.push_back(5.0);
        break;
      case 15:    //P
        am1l.push_back(6.0);
        am1l.push_back(7.0);
        am1l.push_back(9.0);
        break;
      case 16:    //S
        am1l.push_back(4.593691);
        am1l.push_back(5.865731);
        am1l.push_back(13.557336);
        break;
      case 17:    //Cl
        am1l.push_back(4.0);
        am1l.push_back(4.0);
        break;
      case 30:    //Zn
        am1l.push_back(0.0);
        break;
      case 32:    //Ge
        am1l.push_back(0.0);
        break;
      case 33:    //As
        am1l.push_back(4.9433993);
        am1l.push_back(3.1944613);
        break;
      case 34:    //Se
        am1l.push_back(6.5086644);
        am1l.push_back(6.5241228);
        break;
      case 35:    //Br
        am1l.push_back(4.0);
        am1l.push_back(4.0);
        break;
      case 51:    //Sb
        am1l.push_back(6.02795);
        am1l.push_back(3.028109);
        break;
      case 52:    //Te
        am1l.push_back(6.0519413);
        am1l.push_back(3.8304067);
        break;
      case 53:    //I
        am1l.push_back(2.3);
        am1l.push_back(3.0);
        break;
    }
    return am1l;
  }
  std::vector<double> AM1M(int atomicnr, int atm2 = 0) {
    std::vector<double> am1m;
    switch (atomicnr) {
      case 1:     //H
        am1m.push_back(1.2);
        am1m.push_back(1.8);
        am1m.push_back(2.1);
        break;
      case 3:     //Li
        am1m.push_back(0.0);
        break;
      case 4:     //Be
        am1m.push_back(0.0);
        break;
      case 5:     //B
        if (atm2 == 1) {
          am1m.push_back(0.832586);
          am1m.push_back(1.186220);
        }
        else if (atm2 == 6) {
          am1m.push_back(1.063995);
          am1m.push_back(1.936492);
        }
        else if ((atm2 == 9)||(atm2 == 17)||(atm2 == 35)||(atm2 == 53)) {
          am1m.push_back(0.819351);
          am1m.push_back(1.574414);
        }
        else {
          am1m.push_back(0.727592);
          am1m.push_back(1.466639);
          am1m.push_back(1.570975);
        }
        break;
      case 6:     //C
        am1m.push_back(1.60);
        am1m.push_back(1.85);
        am1m.push_back(2.05);
        am1m.push_back(2.65);
        break;
      case 7:     //N
        am1m.push_back(1.5);
        am1m.push_back(2.1);
        am1m.push_back(2.4);
        break;
      case 8:     //O
        am1m.push_back(0.847918);
        am1m.push_back(1.445071);
        break;
      case 9:     //F
        am1m.push_back(0.93);
        am1m.push_back(1.66);
        break;
      case 11:    //Na
        am1m.push_back(1.1681055);
        am1m.push_back(1.1537670);
        break;
      case 13:    //Al
        am1m.push_back(2.050394);
        break;
      case 14:    //Si
        am1m.push_back(0.911453);
        am1m.push_back(1.995569);
        am1m.push_back(2.990610);
        break;
      case 15:    //P
        am1m.push_back(1.474323);
        am1m.push_back(1.779354);
        am1m.push_back(3.006576);
        break;
      case 16:    //S
        am1m.push_back(0.770665);
        am1m.push_back(1.503313);
        am1m.push_back(2.009173);
        break;
      case 17:    //Cl
        am1m.push_back(1.3);
        am1m.push_back(2.1);
        break;
      case 30:    //Zn
        am1m.push_back(0.0);
        break;
      case 32:    //Ge
        am1m.push_back(0.0);
        break;
      case 33:    //As
        am1m.push_back(1.4544264);
        am1m.push_back(2.0144939);
        break;
      case 34:    //Se
        am1m.push_back(1.4981077);
        am1m.push_back(2.0751916);
        break;
      case 35:    //Br
        am1m.push_back(1.5);
        am1m.push_back(2.3);
        break;
      case 51:    //Sb
        am1m.push_back(1.710367);
        am1m.push_back(1.538318);
        break;
      case 52:    //Te
        am1m.push_back(1.3079857);
        am1m.push_back(2.0899707);
        break;
      case 53:    //I
        am1m.push_back(1.8);
        am1m.push_back(2.24);
        break;
    }
    return am1m;
  }
  double Dvalue(size_t atmnr, size_t idx) {
    //function returning the D values needed to calculate eris; values stored in Angstrom, but returned in a.u.
    //note that this function returns both D1 and D2; idx is then either 1 or 2
    double dval = 0.0;
    switch (atmnr) {
      case 1:                    //H
        dval = 0.0;
        break;
      case 3:                  //Li
        if (idx == 1) {dval = 1.0874267;}
        else if (idx == 2) {dval = 0.9227121;}
        break;
      case 4:                  //Be
        if (idx == 1) {dval = 0.7605847;}
        else if (idx == 2) {dval = 0.6453775;}
        break;
      case 5:                  //B
        if (idx == 1) {dval = 0.4819453;}
        else if (idx == 2) {dval = 0.4166779;}
        break;
      case 6:                  //C
        if (idx == 1) {dval = 0.4358609;}
        else if (idx == 2) {dval = 0.3845994;}
        break;
      case 7:                  //N
        if (idx == 1) {dval = 0.3404262;}
        else if (idx == 2) {dval = 0.3003302;}
        break;
      case 8:                  //O
        if (idx == 1) {dval = 0.2639959;}
        else if (idx == 2) {dval = 0.2567689;}
        break;
      case 9:                  //F
        if (idx == 1) {dval = 0.2193505;}
        else if (idx == 2) {dval = 0.2597917;}
        break;
      case 11:                 //Na
        if (idx == 1) {dval = 0.84137879635;}
        else if (idx == 2) {dval = 0.72756674912;}
        break;
      case 13:                 //Al
        if (idx == 1) {dval = 0.7429739;}
        else if (idx == 2) {dval = 0.6778182;}
        break;
      case 14:                 //Si
        if (idx == 1) {dval = 0.6154798;}
        else if (idx == 2) {dval = 0.6891036;}
        break;
      case 15:                 //P
        if (idx == 1) {dval = 0.5530865;}
        else if (idx == 2) {dval = 0.4722106;}
        break;
      case 16:                 //S
        if (idx == 1) {dval = 0.4764760;}
        else if (idx == 2) {dval = 0.5310894;}
        break;
      case 17:                 //Cl
        if (idx == 1) {dval = 0.2860828;}
        else if (idx == 2) {dval = 0.4263609;}
        break;
      case 30:                 //Zn
        if (idx == 1) {dval = 0.7186677;}
        else if (idx == 2) {dval = 0.8179549;}
        break;
      case 32:                 //Ge
        if (idx == 1) {dval = 0.6599821;}
        else if (idx == 2) {dval = 0.5661368;}
        break;
      case 33:                 //As
        if (idx == 1) {dval = 0.6365298;}
        else if (idx == 2) {dval = 0.6507545;}
        break;
      case 34:                 //Se
        if (idx == 1) {dval = 0.5353874;}
        else if (idx == 2) {dval = 0.5474124;}
        break;
      case 35:                 //Br
        if (idx == 1) {dval = 0.4475750;}
        else if (idx == 2) {dval = 0.5507111;}
        break;
      case 51:                 //Sb
        if (idx == 1) {dval = 0.7509809;}
        else if (idx == 2) {dval = 0.6127556;}
        break;
      case 52:                 //Te
        if (idx == 1) {dval = 0.8121000;}
        else if (idx == 2) {dval = 0.6896346;}
        break;
      case 53:                 //I
        if (idx == 1) {dval = 0.7873358;}
        else if (idx == 2) {dval = 0.6290413;}
        break;
      case 80:                 //Hg
        if (idx == 1) {dval = 0.9922320;}
        else if (idx == 2) {dval = 0.8161999;}
        break;
    }
    return dval*dist_Angstrom2aum1;
  }
  double rho(size_t atmnr, size_t l) {
    //function returning the rho values needed to calculate eris; values stored in Angstrom but returned in atomic units
    double rho = 0.0;
    switch (atmnr) {
      case 1:                  //H
        if (l == 0) {rho = 0.560345403;}
        else if (l == 1) {rho = 0.0;}
        else if (l == 2) {rho = 0.0;}
        break;
      case 3:                  //Li
        if (l == 0) {rho = 0.986207884;}
        else if (l == 1) {rho = 1.165672376;}
        else if (l == 2) {rho = 1.011953732;}
        break;
      case 4:                  //Be
        if (l == 0) {rho = 0.799924235;}
        else if (l == 1) {rho = 0.788356095;}
        else if (l == 2) {rho = 0.687877905;}
        break;
      case 5:                  //B
        if (l == 0) {rho = 0.679822279;}
        else if (l == 1) {rho = 0.524431177;}
        else if (l == 2) {rho = 0.465909859;}
        break;
      case 6:                  //C
        if (l == 0) {rho = 0.588660438;}
        else if (l == 1) {rho = 0.434959475;}
        else if (l == 2) {rho = 0.411899789;}
        break;
      case 7:                  //N
        if (l == 0) {rho = 0.529751104;}
        else if (l == 1) {rho = 0.338305732;}
        else if (l == 2) {rho = 0.335616880;}
        break;
      case 8:                  //O
        if (l == 0) {rho = 0.466881794;}
        else if (l == 1) {rho = 0.265617656;}
        else if (l == 2) {rho = 0.291866510;}
        break;
      case 9:                  //F
        if (l == 0) {rho = 0.425491557;}
        else if (l == 1) {rho = 0.218866782;}
        else if (l == 2) {rho = 0.280001651;}
        break;
      case 11:                 //Na
        if (l == 0) {rho = 0.9800619701;}
        else if (l == 1) {rho = 1.63523435444;}
        else if (l == 2) {rho = 1.41814447768;}
        break;
      case 13:                 //Al
        if (l == 0) {rho = 0.889903107;}
        else if (l == 1) {rho = 1.005933324;}
        else if (l == 2) {rho = 0.771868341;}
        break;
      case 14:                 //Si
        if (l == 0) {rho = 0.733128067;}
        else if (l == 1) {rho = 0.690852269;}
        else if (l == 2) {rho = 0.712758472;}
        break;
      case 15:                 //P
        if (l == 0) {rho = 0.622778008;}
        else if (l == 1) {rho = 0.807809865;}
        else if (l == 2) {rho = 0.603128119;}
        break;
      case 16:                 //S
        if (l == 0) {rho = 0.610819239;}
        else if (l == 1) {rho = 0.447906465;}
        else if (l == 2) {rho = 0.409892853;}
        break;
      case 17:                 //Cl
        if (l == 0) {rho = 0.478996434;}
        else if (l == 1) {rho = 0.343918655;}
        else if (l == 2) {rho = 0.431383633;}
        break;
      case 30:                 //Zn
        if (l == 0) {rho = 0.610111605;}
        else if (l == 1) {rho = 1.141714309;}
        else if (l == 2) {rho = 1.009411845;}
        break;
      case 32:                 //Ge
        if (l == 0) {rho = 0.707994522;}
        else if (l == 1) {rho = 0.831942745;}
        else if (l == 2) {rho = 0.759073299;}
        break;
      case 33:                 //As
        if (l == 0) {rho = 0.64880765;}
        else if (l == 1) {rho = 0.95995755;}
        else if (l == 2) {rho = 0.927371034;}
        break;
      case 34:                 //Se
        if (l == 0) {rho = 1.060143453;}
        else if (l == 1) {rho = 0.358473620;}
        else if (l == 2) {rho = 0.527366991;}
        break;
      case 35:                 //Br
        if (l == 0) {rho = 0.478791351;}
        else if (l == 1) {rho = 0.439172041;}
        else if (l == 2) {rho = 0.498503548;}
        break;
      case 51:                 //Sb
        if (l == 0) {rho = 0.629847642;}
        else if (l == 1) {rho = 1.086004527;}
        else if (l == 2) {rho = 0.974630761;}
        break;
      case 52:                 //Te
        if (l == 0) {rho = 1.442019673;}
        else if (l == 1) {rho = 0.457437268;}
        else if (l == 2) {rho = 0.551264090;}
        break;
      case 53:                 //I
        if (l == 0) {rho = 0.478663761;}
        else if (l == 1) {rho = 0.588287153;}
        else if (l == 2) {rho = 0.571235649;}
        break;
      case 80:                 //Hg
        if (l == 0) {rho = 0.666603429;}
        else if (l == 1) {rho = 0.904062899;}
        else if (l == 2) {rho = 0.787310536;}
        break;
    }
    return rho*dist_Angstrom2aum1;
  }
  double eri1Center(int atmnr, int Lbra, int Lket) {
    //function that gives back the semi-empirical 1-center eris
    //values stored in eV but returned in a.u.
    //Lbra is the sum of azimuthal quantum numbers for bra (ss = 0; pp = 2; sp = 1)
    //Lket is the sum of azimuthal quantum numbers for ket (ss = 0; pp = 2; sp = 1; p*p* = -2)
    double eri = 0.0;
    switch (atmnr) {
      case 1:         //H
        eri = 12.848;                                                                            //(ss|ss)
        break;
      case 3:         //Li
        if ((Lbra == 0)&&(Lket == 0)) {eri = 7.3;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.0;}                                          //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.42;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.52;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.83;}                                         //(sp|sp)||(ps|ps)
        break;
      case 4:         //Be
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.0;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.97;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.43;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.22;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.28;}                                         //(sp|sp)||(ps|ps)
        break;
      case 5:         //B
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.59;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.86;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 9.56;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.86;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.81;}                                         //(sp|sp)||(ps|ps)
        break;
      case 6:         //C
        if ((Lbra == 0)&&(Lket == 0)) {eri = 12.23;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.08;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 11.47;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.84;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.43;}                                         //(sp|sp)||(ps|ps)
        break;
      case 7:         //N
        if ((Lbra == 0)&&(Lket == 0)) {eri = 13.59;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 12.98;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 12.66;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 11.59;}        //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.14;}                                         //(sp|sp)||(ps|ps)
        break;
      case 8:         //O
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.42;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 14.52;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 14.48;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 12.98;}        //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.94;}                                         //(sp|sp)||(ps|ps)
        break;
      case 9:         //F
        if ((Lbra == 0)&&(Lket == 0)) {eri = 16.92;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 16.71;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 17.25;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 14.91;}        //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 4.83;}                                         //(sp|sp)||(ps|ps)
        break;
      case 11:        //Na
        if ((Lbra == 0)&&(Lket == 0)) {eri = 7.3459178;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 4.1130516;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.4042550;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.0370957;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.2487699;}                                    //(sp|sp)||(ps|ps)
        break;
      case 13:        //Al
        if ((Lbra == 0)&&(Lket == 0)) {eri = 8.09;}                                              //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.98;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.63;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.40;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.70;}                                         //(sp|sp)||(ps|ps)
        break;
      case 14:        //Si
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.82;}                                              //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.31;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.36;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.54;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.32;}                                         //(sp|sp)||(ps|ps)
        break;
      case 15:        //P
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.560005;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.877589;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.237449;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.307648;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.779238;}                                     //(sp|sp)||(ps|ps)
        break;
      case 16:        //S
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.786329;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 10.039308;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.663127;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.781688;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.532137;}                                     //(sp|sp)||(ps|ps)
        break;
      case 17:        //Cl
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.03;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.30;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 13.16;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.97;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.42;}                                         //(sp|sp)||(ps|ps)
        break;
      case 30:        //Zn
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.8;}                                              //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 13.3;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 11.182018;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 12.93052;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.484606;}                                     //(sp|sp)||(ps|ps)
        break;
      case 32:        //Ge
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.168605;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.671902;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.144473;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.269706;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.937093;}                                   //(sp|sp)||(ps|ps)
        break;
      case 33:        //As
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.0962258;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.8781648;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 4.9259328;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.5961088;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.6246173;}                                    //(sp|sp)||(ps|ps)
        break;
      case 34:        //Se
        if ((Lbra == 0)&&(Lket == 0)) {eri = 6.7908891;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.4769273;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.4812786;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.2796993;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 4.4548356;}                                    //(sp|sp)||(ps|ps)
        break;
      case 35:        //Br
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.03643948;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.27632539;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 13.03468242;}    //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.85442552;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.45586832;}                                   //(sp|sp)||(ps|ps)
        break;
      case 51:        //Sb
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.430251;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.424094;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.787922;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.849181;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.58834;}                                      //(sp|sp)||(ps|ps)
        break;
      case 52:        //Te
        if ((Lbra == 0)&&(Lket == 0)) {eri = 4.9925231;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.2097852;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 4.9721484;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.6211521;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 4.0071821;}                                    //(sp|sp)||(ps|ps)
        break;
      case 53:        //I
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.04044855;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.14778369;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 13.05655798;}    //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.91409071;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.45638202;}                                   //(sp|sp)||(ps|ps)
        break;
      case 80:        //Hg
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.8;}                                              //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 14.3;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 9.3;}            //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 13.5;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.3;}                                          //(sp|sp)||(ps|ps)
        break;
    }
    return eri/au2eV;
  }
  //charge model information; this is separate calculation
  double CM1ck(size_t atomicnr) {
    //function returning the atomic parameters ck for CM1
    double cm1ck = 0.0;
    if ((atomicnr == 1)||((atomicnr > 5)&&(atomicnr < 10))||(atomicnr == 14)||(atomicnr == 16)||(atomicnr == 17)||(atomicnr == 35)||(atomicnr == 53)) {
      if (atomicnr == 7) {cm1ck = 0.3846;}                 //N
      else if (atomicnr == 9) {cm1ck = 0.1468;}            //F
      else if (atomicnr == 16) {cm1ck = -0.1311;}          //S
      else if (atomicnr == 17) {cm1ck = 0.0405;}           //Cl
      else if (atomicnr == 35) {cm1ck = 0.1761;}           //Br
      else if (atomicnr == 53) {cm1ck = 0.2380;}           //I
    }
    else {std::cout << "WARNING: MNDO.hpp: AM1: CM1ck(): partial charge calculation involving non-parametrized atom";}
    return cm1ck;
  }
  double CM1dk(size_t atomicnr) {
    //function returning the atomic parameters dk for CM1
    double cm1dk = 0.0;
    if ((atomicnr == 1)||((atomicnr > 5)&&(atomicnr < 10))||(atomicnr == 14)||(atomicnr == 16)||(atomicnr == 17)||(atomicnr == 35)||(atomicnr == 53)) {
      if (atomicnr == 8) {cm1dk = -0.0283;}                //O
      else if (atomicnr == 9) {cm1dk = 0.0399;}            //F
      else if (atomicnr == 16) {cm1dk = -0.0956;}          //S
      else if (atomicnr == 17) {cm1dk = -0.0276;}          //Cl
      else if (atomicnr == 35) {cm1dk = -0.0802;}          //Br
      else if (atomicnr == 53) {cm1dk = -0.1819;}          //I
    }
    else {std::cout << "WARNING: MNDO.hpp: AM1: CM1dk(): partial charge calculation involving non-parametrized atom";}
    return cm1dk;
  }
  double CM1ckkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters ckk' for CM1
    size_t atomA = std::min(atom1,atom2);
    size_t atomB = std::max(atom1,atom2);
    double cm1ckkp = 0.0;
    if ((atomA == 6)&&(atomB == 7)) {cm1ckkp = 0.3846;}           //C-N
    return cm1ckkp;
  }
  double CM1dkkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters dkk' for CM1
    double cm1dkkp = 0.0;
    bool valid1 = false;
    bool valid2 = false;
    if ((atom1 == 1)||((atom1 > 5)&&(atom1 < 10))||(atom1 == 14)||(atom1 == 16)||(atom1 == 17)||(atom1 == 35)||(atom1 == 53)) {valid1 = true;}
    if ((atom2 == 1)||((atom2 > 5)&&(atom2 < 10))||(atom2 == 14)||(atom2 == 16)||(atom2 == 17)||(atom2 == 35)||(atom2 == 53)) {valid2 = true;}
    if ((valid1)&&(valid2)) {
      size_t atomA = std::min(atom1,atom2);
      size_t atomB = std::max(atom1,atom2);
      if ((atomA == 1)&&(atomB == 7)) {cm1dkkp = 0.0850;}           //H-N
      else if ((atomA == 1)&&(atomB == 8)) {cm1dkkp = 0.1447;}      //H-O
      else if ((atomA == 1)&&(atomB == 14)) {cm1dkkp = 0.0640;}     //H-Si
      else if ((atomA == 6)&&(atomB == 7)) {cm1dkkp = -0.0880;}     //C-N
      else if ((atomA == 7)&&(atomB == 8)) {cm1dkkp = -0.0630;}     //N-O
      else if ((atomA == 8)&&(atomB == 16)) {cm1dkkp = -0.0600;}    //O-S
    }
    else {std::cout << "WARNING: MNDO.hpp: AM1: CM1dkkp(): partial charge calculation involving non-parametrized atom";}
    return cm1dkkp;
  }
  double CM2ckkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters ckk' for CM2
    double cm2ckkp = 0.0;
    size_t atomA = std::min(atom1,atom2);
    size_t atomB = std::max(atom1,atom2);
    if ((atomA == 1)&&(atomB == 6)) {cm2ckkp = -0.020;}           //H-C
    else if ((atomA == 1)&&(atomB == 7)) {cm2ckkp = 0.207;}       //H-N
    else if ((atomA == 1)&&(atomB == 8)) {cm2ckkp = 0.177;}       //H-O
    else if ((atomA == 1)&&(atomB == 14)) {cm2ckkp = -0.083;}     //H-Si
    else if ((atomA == 1)&&(atomB == 16)) {cm2ckkp = 0.038;}      //H-S
    else if ((atomA == 6)&&(atomB == 7)) {cm2ckkp = 0.008;}       //C-N
    else if ((atomA == 6)&&(atomB == 8)) {cm2ckkp = 0.026;}       //C-O
    else if ((atomA == 6)&&(atomB == 14)) {cm2ckkp = 0.062;}      //C-Si
    else if ((atomA == 6)&&(atomB == 16)) {cm2ckkp = -0.059;}     //C-S
    else if ((atomA == 7)&&(atomB == 8)) {cm2ckkp = -0.197;}      //N-O
    if (atom1 != atomA) {cm2ckkp *= -1.0;}
    return cm2ckkp;
  }
  double CM2dkkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters dkk' for CM2
    double cm2dkkp = 0.0;
    bool valid1 = false;
    bool valid2 = false;
    if ((atom1 == 1)||((atom1 > 5)&&(atom1 < 10))||((atom1 > 13)&&(atom1 < 18))||(atom1 == 35)||(atom1 == 53)) {valid1 = true;}
    if ((atom2 == 1)||((atom2 > 5)&&(atom2 < 10))||((atom2 > 13)&&(atom2 < 18))||(atom2 == 35)||(atom2 == 53)) {valid2 = true;}
    if ((valid1)&&(valid2)) {
      size_t atomA = std::min(atom1,atom2);
      size_t atomB = std::max(atom1,atom2);
      if ((atomA == 1)&&(atomB == 15)) {cm2dkkp = 0.103;}           //H-P
      else if ((atomA == 6)&&(atomB == 7)) {cm2dkkp = 0.086;}       //C-N
      else if ((atomA == 6)&&(atomB == 8)) {cm2dkkp = 0.016;}       //C-O
      else if ((atomA == 6)&&(atomB == 9)) {cm2dkkp = 0.019;}       //C-F
      else if ((atomA == 6)&&(atomB == 15)) {cm2dkkp = -0.019;}     //C-P
      else if ((atomA == 6)&&(atomB == 16)) {cm2dkkp = 0.171;}      //C-S
      else if ((atomA == 6)&&(atomB == 17)) {cm2dkkp = 0.027;}      //C-Cl
      else if ((atomA == 6)&&(atomB == 35)) {cm2dkkp = 0.081;}      //C-Br
      else if ((atomA == 6)&&(atomB == 53)) {cm2dkkp = 0.147;}      //C-I
      else if ((atomA == 7)&&(atomB == 8)) {cm2dkkp = 0.134;}       //N-O
      else if ((atomA == 8)&&(atomB == 15)) {cm2dkkp = 0.088;}      //O-P
      else if ((atomA == 9)&&(atomB == 15)) {cm2dkkp = 0.252;}      //F-P
      else if ((atomA == 15)&&(atomB == 16)) {cm2dkkp = -0.080;}    //P-S
      if (atom1 != atomA) {cm2dkkp *= -1.0;}
    }
    else {std::cout << "WARNING: MNDO.hpp: AM1: CM2dkkp(): partial charge calculation involving non-parametrized atom";}
    return cm2dkkp;
  }
  double CM3ckkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters ckk' for CM3
    double cm3ckkp = 0.0;
    size_t atomA = std::min(atom1,atom2);
    size_t atomB = std::max(atom1,atom2);
    if ((atomA == 6)&&(atomB == 8)) {cm3ckkp = 0.029;}            //C-O
    else if ((atomA == 7)&&(atomB == 8)) {cm3ckkp = 0.439;}       //N-O
    else if ((atomA == 8)&&(atomB == 14)) {cm3ckkp = 0.201;}      //O-Si
    else if ((atomA == 8)&&(atomB == 15)) {cm3ckkp = 0.004;}      //O-P
    else if ((atomA == 15)&&(atomB == 16)) {cm3ckkp = 0.341;}     //P-S
    if (atom1 != atomA) {cm3ckkp *= -1.0;}
    return cm3ckkp;
  }
  double CM3dkkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters dkk' for CM3
    double cm3dkkp = 0.0;
    bool valid1 = false;
    bool valid2 = false;
    if ((atom1 == 1)||((atom1 > 5)&&(atom1 < 10))||((atom1 > 13)&&(atom1 < 18))||(atom1 == 35)) {valid1 = true;}
    if ((atom2 == 1)||((atom2 > 5)&&(atom2 < 10))||((atom2 > 13)&&(atom2 < 18))||(atom2 == 35)) {valid2 = true;}
    if ((valid1)&&(valid2)) {
      size_t atomA = std::min(atom1,atom2);
      size_t atomB = std::max(atom1,atom2);
      if ((atomA == 1)&&(atomB == 6)) {cm3dkkp = -0.009;}           //H-C
      else if ((atomA == 1)&&(atomB == 7)) {cm3dkkp = 0.200;}       //H-N
      else if ((atomA == 1)&&(atomB == 8)) {cm3dkkp = 0.154;}       //H-O
      else if ((atomA == 1)&&(atomB == 14)) {cm3dkkp = 0.231;}      //H-Si
      else if ((atomA == 1)&&(atomB == 15)) {cm3dkkp = 0.078;}      //H-P
      else if ((atomA == 1)&&(atomB == 16)) {cm3dkkp = 0.072;}      //H-S
      else if ((atomA == 6)&&(atomB == 7)) {cm3dkkp = 0.091;}       //C-N
      else if ((atomA == 6)&&(atomB == 8)) {cm3dkkp = -0.012;}      //C-O
      else if ((atomA == 6)&&(atomB == 9)) {cm3dkkp = 0.006;}       //C-F
      else if ((atomA == 6)&&(atomB == 14)) {cm3dkkp = 0.163;}      //C-Si
      else if ((atomA == 6)&&(atomB == 15)) {cm3dkkp = 0.106;}      //C-P
      else if ((atomA == 6)&&(atomB == 16)) {cm3dkkp = 0.103;}      //C-S
      else if ((atomA == 6)&&(atomB == 17)) {cm3dkkp = 0.013;}      //C-Cl
      else if ((atomA == 6)&&(atomB == 35)) {cm3dkkp = 0.073;}      //C-Br
      else if ((atomA == 7)&&(atomB == 8)) {cm3dkkp = -0.132;}      //N-O
      else if ((atomA == 7)&&(atomB == 15)) {cm3dkkp = -0.067;}     //N-P
      else if ((atomA == 8)&&(atomB == 14)) {cm3dkkp = 0.169;}      //O-Si
      else if ((atomA == 8)&&(atomB == 15)) {cm3dkkp = 0.033;}      //O-P
      else if ((atomA == 8)&&(atomB == 16)) {cm3dkkp = 0.103;}      //O-S
      else if ((atomA == 9)&&(atomB == 14)) {cm3dkkp = 0.228;}      //F-Si
      else if ((atomA == 9)&&(atomB == 15)) {cm3dkkp = 0.110;}      //F-P
      else if ((atomA == 14)&&(atomB == 17)) {cm3dkkp = -0.171;}    //Si-Cl
      else if ((atomA == 15)&&(atomB == 16)) {cm3dkkp = -0.685;}    //P-S
      else if ((atomA == 15)&&(atomB == 17)) {cm3dkkp = -0.252;}    //P-Cl
      if (atom1 != atomA) {cm3dkkp *= -1.0;}
    }
    else {std::cout << "WARNING: MNDO.hpp: AM1: CM3dkkp(): partial charge calculation involving non-parametrized atom";}
    return cm3dkkp;
  }
  double B0() {return 0.482;}
};
class PM3: public MNDO {
  //this is the implementation of Stewart's PM3
  //J. J. P. Stewart, J. Comput. Chem., 10(2), 209, 1989
  //J. J. P. Stewart, J. Comput. Chem., 12(3), 320, 1991
  //E. Anders, R. Koch, P. Freunscht, J. Comput. Chem., 14(11), 1301, 1993
  //for charge and dipole moment model: J. W. Storer, D. J. Giesen, C. J. Cramer, D. G. Truhlar, J. Comp.-Aid. Mol. Des., 9, 87, 1995
  //                                    J. Li, T. Zhu, C. J. Cramer, D. G. Truhlar, J. Phys. Chem. A, 102, 1820, 1998                   -> CM2
  //bX quantities are "barred" tensors, which are used only for the open-shell case
public:
  PM3(BSet _bset, Molecule _mol, std::string _openclosed = "0", std::string _corecorrection = "0"): MNDO(_bset,_mol,_openclosed,_corecorrection) {}
  ~PM3() {}
  std::string Type() {return "PM3";}
  double gfactor(size_t iatm1, size_t iatm2, double RAB) {return 1.0;}
  virtual double gfactor_dR(size_t iatm1, size_t iatm2) {return 0.0;}
  double AM1factor(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {
    //RAB in a.u.
    double AM1factor = 0.0;
    RAB *= dist_Angstrom2au;
    for (size_t idx = 0; idx < K1.size(); ++idx) {
      AM1factor += K1[idx]*exp(-L1[idx]*(RAB - M1[idx])*(RAB - M1[idx]))/RAB;
    }
    return AM1factor;
  }
  virtual double AM1factor_dR(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {
    //RAB in a.u.
    double am1factor_dr = 0.0;
    double RR = RAB*dist_Angstrom2au;
    for (size_t idx = 0; idx < K1.size(); ++idx) {
      am1factor_dr += K1[idx]*(2.0*RR*RR*L1[idx] - 2.0*L1[idx]*M1[idx]*RR + 1.0)*exp(-L1[idx]*(RR - M1[idx])*(RR - M1[idx]))/(RAB*RR);
    }
    return am1factor_dr;  
  }
  virtual double AM1factor_dR2(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {
    //RAB in a.u.
    double am1factor_dr2 = 0.0;
    double RR = RAB*dist_Angstrom2au;
    for (size_t idx = 0; idx < K1.size(); ++idx) {
      am1factor_dr2 += 2.0*K1[idx]*exp(-L1[idx]*(RR - M1[idx])*(RR - M1[idx]))*(L1[idx]*M1[idx]*RR - 1.0 - RR*L1[idx]*(RR - M1[idx])*(2.0*RR*RR*L1[idx] - 2.0*L1[idx]*M1[idx]*RR + 1.0))/(RAB*RAB*RR);
    }
    return am1factor_dr2;  
  }
  //other auxiliary functions
  void checkAtoms() {
    //function that checks whether the atoms in the molecule are fully defined in the respective MNDO theories
    bool def;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      def = false;
      if (atoms[iatm] == 1) {def = true;}                                   //H
      else if ((atoms[iatm] > 2)&&(atoms[iatm] < 10)) {def = true;}         //Li,Be,B,C,N,O,F
      else if ((atoms[iatm] > 10)&&(atoms[iatm] < 18)) {def = true;}        //Na,Mg,Al,Si,P,S,Cl
      else if ((atoms[iatm] > 29)||(atoms[iatm] < 36)) {def = true;}        //Zn,Ga,Ge,As,Se,Br
      else if ((atoms[iatm] > 47)||(atoms[iatm] < 54)) {def = true;}        //Cd-I
      else if ((atoms[iatm] > 79)||(atoms[iatm] < 84)) {def = true;}        //Hg,Tl,Pb,Bi
      if (atoms[iatm] == 5) {def = false;}                                  //remove B
      if (!def) {throw("ERROR: MNDO.hpp: PM3: checkAtoms(): atom not fully specified for PM3-theory");}
    }
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapPM3(atm);}
  double ElementFormTheo(size_t atomicnr) {
    //function returning theoretical heats of formation for elements; values directly in a.u.
    double enth = 0.0;
    switch (atomicnr) {
      case 1:      //H
        enth = -0.480435734286;
        break;
      case 3:      //Li
        enth = -0.194771427376;
        break;
      case 4:      //Be
        enth = -0.937719797485;
        break;
      case 6:      //C
        enth = -4.05540034771;
        break;
      case 7:      //N
        enth = -5.77479194152;
        break;
      case 8:      //O
        enth = -10.6102009213;
        break;
      case 9:      //F
        enth = -16.0784610392;
        break;
      case 11:     //Na
        enth = -0.19457272009;
        break;
      case 12:     //Mg
        enth = -0.828810340423;
        break;
      case 13:     //Al
        enth = -1.72224844965;
        break;
      case 14:     //Si
        enth = -2.46180717218;
        break;
      case 15:     //P
        enth = -4.32469507163;
        break;
      case 16:     //S
        enth = -6.70509021263;
        break;
      case 17:     //Cl
        enth = -11.5832018724;
        break;
      case 30:     //Zn
        enth = -1.00646113884;
        break;
      case 31:     //Ga
        enth = -2.10676627507;
        break;
      case 32:     //Ge
        enth = -3.07377982018;
        break;
      case 33:     //As
        enth = -4.50525281523;
        break;
      case 34:     //Se
        enth = -7.05046081862;
        break;
      case 35:     //Br
        enth = -12.9556035747;
        break;
      case 48:     //Cd
        enth = -0.825030010764;
        break;
      case 49:     //In
        enth = -1.91004794191;
        break;
      case 50:     //Sn
        enth = -2.89003979076;
        break;
      case 51:     //Sb
        enth = -5.47154425813;
        break;
      case 52:     //Te
        enth = -6.17695006812;
        break;
      case 53:     //I
        enth = -10.5954218289;
        break;
      case 80:     //Hg
        enth = -1.06204589076;
        break;
      case 81:     //Tl
        enth = -2.08182009765;
        break;
      case 82:     //Pb
        enth = -2.69729075229;
        break;
      case 83:     //Bi
        enth = -4.009248964;
        break;
    }
    return enth;
  }
  double betaA0(size_t atomicnr, int L) {
    //function returning the bonding parameters; values stored in eV, returned however in a.u.
    double betaa0 = 0.0;
    switch (atomicnr) {
      case 1:         //H
        betaa0 = -5.6265120;
        break;
      case 3:         //Li
        if (L == 0) {betaa0 = -0.5500000;}
        else if (L == 1) {betaa0 = -1.5000000;}
        break;
      case 4:         //Be
        if (L == 0) {betaa0 = -3.9620530;}
        else if (L == 1) {betaa0 = -2.7806840;}
        break;
      case 6:         //C
        if (L == 0) {betaa0 = -11.9100150;}
        else if (L == 1) {betaa0 = -9.8027550;}
        break;
      case 7:         //N
        if (L == 0) {betaa0 = -14.0625210;}
        else if (L == 1) {betaa0 = -20.0438480;}
        break;
      case 8:         //O
        if (L == 0) {betaa0 = -45.2026510;}
        else if (L == 1) {betaa0 = -24.7525150;}
        break;
      case 9:         //F
        if (L == 0) {betaa0 = -48.4059390;}
        else if (L == 1) {betaa0 = -27.7446600;}
        break;
      case 11:        //Na
        if (L == 0) {betaa0 = -0.1510870;}
        else if (L == 1) {betaa0 = -0.2184096;}
        break;
      case 12:        //Mg
        if (L == 0) {betaa0 = -2.0716910;}
        else if (L == 1) {betaa0 = -0.5695810;}
        break;
      case 13:        //Al
        if (L == 0) {betaa0 = -0.5943010;}
        else if (L == 1) {betaa0 = -0.9565500;}
        break;
      case 14:        //Si
        if (L == 0) {betaa0 = -2.8621450;}
        else if (L == 1) {betaa0 = -3.9331480;}
        break;
      case 15:        //P
        if (L == 0) {betaa0 = -12.6158790;}
        else if (L == 1) {betaa0 = -4.1600400;}
        break;
      case 16:        //S
        if (L == 0) {betaa0 = -8.8274650;}
        else if (L == 1) {betaa0 = -8.0914150;}
        break;
      case 17:        //Cl
        if (L == 0) {betaa0 = -27.5285600;}
        else if (L == 1) {betaa0 = -11.5939220;}
        break;
      case 30:        //Zn
        if (L == 0) {betaa0 = -0.7155780;}
        else if (L == 1) {betaa0 = -6.3518640;}
        break;
      case 31:        //Ga
        if (L == 0) {betaa0 = -4.9456180;}
        else if (L == 1) {betaa0 = -0.4070530;}
        break;
      case 32:        //Ge
        if (L == 0) {betaa0 = -5.3250024;}
        else if (L == 1) {betaa0 = -2.2501567;}
        break;
      case 33:        //As
        if (L == 0) {betaa0 = -8.2321650;}
        else if (L == 1) {betaa0 = -5.0173860;}
        break;
      case 34:        //Se
        if (L == 0) {betaa0 = -6.1578220;}
        else if (L == 1) {betaa0 = -5.4930390;}
        break;
      case 35:        //Br
        if (L == 0) {betaa0 = -31.1713420;}
        else if (L == 1) {betaa0 = -6.8140130;}
        break;
      case 48:        //Cd
        if (L == 0) {betaa0 = -8.5819440;}
        else if (L == 1) {betaa0 = -0.6010340;}
        break;
      case 49:        //In
        if (L == 0) {betaa0 = -2.9933190;}
        else if (L == 1) {betaa0 = -1.8289080;}
        break;
      case 50:        //Sn
        if (L == 0) {betaa0 = -2.7858020;}
        else if (L == 1) {betaa0 = -2.0059990;}
        break;
      case 51:        //Sb
        if (L == 0) {betaa0 = -14.7942170;}
        else if (L == 1) {betaa0 = -2.8179480;}
        break;
      case 52:        //Te
        if (L == 0) {betaa0 = -2.6651460;}
        else if (L == 1) {betaa0 = -3.8954300;}
        break;
      case 53:        //I
        if (L == 0) {betaa0 = -14.4942340;}
        else if (L == 1) {betaa0 = -5.8947030;}
        break;
      case 80:        //Hg
        if (L == 0) {betaa0 = -3.1013650;}
        else if (L == 1) {betaa0 = -3.4640310;}
        break;
      case 81:        //Tl
        if (L == 0) {betaa0 = -1.0844950;}
        else if (L == 1) {betaa0 = -7.9467990;}
        break;
      case 82:        //Pb
        if (L == 0) {betaa0 = -6.1260240;}
        else if (L == 1) {betaa0 = -1.3954300;}
        break;
      case 83:        //Bi
        if (L == 0) {betaa0 = -5.6072830;}
        else if (L == 1) {betaa0 = -5.8001520;}
        break;
    }
    return betaa0/au2eV;
  }
  double UlX(size_t atomicnr, int L) {
    //function that returns the atomic potential U
    double ulx = 0.0;
    switch (atomicnr) {
      case 1:                 //H
        ulx = -13.0733210;
        break;
      case 3:                 //Li
        if (L == 0) {ulx = -5.3000000;}
        else if (L == 1) {ulx = -3.4000000;}
        break;
      case 4:                 //Be
        if (L == 0) {ulx = -17.2647520;}
        else if (L == 1) {ulx = -11.3042430;}
        break;
      case 6:                 //C
        if (L == 0) {ulx = -47.2703200;}
        else if (L == 1) {ulx = -36.2669180;}
        break;
      case 7:                 //N
        if (L == 0) {ulx = -49.3356720;}
        else if (L == 1) {ulx = -47.5097360;}
        break;
      case 8:                 //O
        if (L == 0) {ulx = -86.9930020;}
        else if (L == 1) {ulx = -71.8795800;}
        break;
      case 9:                 //F
        if (L == 0) {ulx = -110.4353030;}
        else if (L == 1) {ulx = -105.6850470;}
        break;
      case 11:                //Na
        if (L == 0) {ulx = -5.2945929;}
        else if (L == 1) {ulx = -2.4596564;}
        break;
      case 12:                //Mg
        if (L == 0) {ulx = -14.6236880;}
        else if (L == 1) {ulx = -14.1734600;}
        break;
      case 13:                //Al
        if (L == 0) {ulx = -24.8454040;}
        else if (L == 1) {ulx = -22.2641590;}
        break;
      case 14:                //Si
        if (L == 0) {ulx = -26.7634830;}
        else if (L == 1) {ulx = -22.8136350;}
        break;
      case 15:                //P
        if (L == 0) {ulx = -40.4130960;}
        else if (L == 1) {ulx = -29.5930520;}
        break;
      case 16:                //S
        if (L == 0) {ulx = -49.8953710;}
        else if (L == 1) {ulx = -44.3925830;}
        break;
      case 17:                //Cl
        if (L == 0) {ulx = -100.6267470;}
        else if (L == 1) {ulx = -53.6143960;}
        break;
      case 30:                //Zn
        if (L == 0) {ulx = -18.5321980;}
        else if (L == 1) {ulx = -11.0474090;}
        break;
      case 31:                //Ga
        if (L == 0) {ulx = -29.8555930;}
        else if (L == 1) {ulx = -21.8753710;}
        break;
      case 32:                //Ge
        if (L == 0) {ulx = -35.4671955;}
        else if (L == 1) {ulx = -31.5863583;}
        break;
      case 33:                //As
        if (L == 0) {ulx = -38.5074240;}
        else if (L == 1) {ulx = -35.1524150;}
        break;
      case 34:                //Se
        if (L == 0) {ulx = -55.3781350;}
        else if (L == 1) {ulx = -49.8230760;}
        break;
      case 35:                //Br
        if (L == 0) {ulx = -116.6193110;}
        else if (L == 1) {ulx = -74.2271290;}
        break;
      case 48:                //Cd
        if (L == 0) {ulx = -15.8285840;}
        else if (L == 1) {ulx = 8.7497950;}
        break;
      case 49:                //In
        if (L == 0) {ulx = -26.1762050;}
        else if (L == 1) {ulx = -20.0058220;}
        break;
      case 50:                //Sn
        if (L == 0) {ulx = -34.5501920;}
        else if (L == 1) {ulx = -25.8944190;}
        break;
      case 51:                //Sb
        if (L == 0) {ulx = -56.4321960;}
        else if (L == 1) {ulx = -29.4349540;}
        break;
      case 52:                //Te
        if (L == 0) {ulx = -44.9380360;}
        else if (L == 1) {ulx = -46.3140990;}
        break;
      case 53:                //I
        if (L == 0) {ulx = -96.4540370;}
        else if (L == 1) {ulx = -61.09158209;}
        break;
      case 80:                //Hg
        if (L == 0) {ulx = -17.7622290;}
        else if (L == 1) {ulx = -18.3307510;}
        break;
      case 81:                //Tl
        if (L == 0) {ulx = -30.0531700;}
        else if (L == 1) {ulx = -26.9206370;}
        break;
      case 82:                //Pb
        if (L == 0) {ulx = -30.3227560;}
        else if (L == 1) {ulx = -24.4258340;}
        break;
      case 83:                //Bi
        if (L == 0) {ulx = -33.4959380;}
        else if (L == 1) {ulx = -35.5210260;}
        break;
    }
    return ulx/au2eV;
  }
  double alpha(int atomicnr, int atm2) {
    //function that returns alpha values for PM3
    //alphas stored and returned in 1/Angstrom
    double alpha = 0.0;
    switch (atomicnr) {
      case 1:     //H
        alpha = 3.356386;
        break;
      case 3:     //Li
        alpha = 1.255000;
        break;
      case 4:     //Be
        alpha = 1.593536;
        break;
      case 6:     //C
        alpha = 2.707807;
        break;
      case 7:     //N
        alpha = 2.830545;
        break;
      case 8:     //O
        alpha = 3.217102;
        break;
      case 9:     //F
        alpha = 3.358921;
        break;
      case 11:    //Na
        alpha = 2.3677169;
        break;
      case 12:    //Mg
        alpha = 1.329147;
        break;
      case 13:    //Al
        alpha = 1.521703;
        break;
      case 14:    //Si
        alpha = 2.135809;
        break;
      case 15:    //P
        alpha = 1.940534;
        break;
      case 16:    //S
        alpha = 2.269706;
        break;
      case 17:    //Cl
        alpha = 2.517296;
        break;
      case 30:    //Zn
        alpha = 1.350126;
        break;
      case 31:    //Ga
        alpha = 1.605115;
        break;
      case 32:    //Ge
        alpha = 1.972337;
        break;
      case 33:    //As
        alpha = 1.794477;
        break;
      case 34:    //Se
        alpha = 3.043957;
        break;
      case 35:    //Br
        alpha = 2.511842;
        break;
      case 48:    //Cd
        alpha = 1.525382;
        break;
      case 49:    //In
        alpha = 1.418385;
        break;
      case 50:    //Sn
        alpha = 1.699650;
        break;
      case 51:    //Sb
        alpha = 2.034301;
        break;
      case 52:    //Te
        alpha = 2.485019;
        break;
      case 53:    //I
        alpha = 1.990185;
        break;
      case 80:    //Hg
        alpha = 1.529377;
        break;
      case 81:    //Tl
        alpha = 1.340951;
        break;
      case 82:    //Pb
        alpha = 1.620045;
        break;
      case 83:    //Bi
        alpha = 1.857431;
        break;
    }
    return alpha;
  }
  std::vector<double> AM1K(int atomicnr, int atm2 = 0) {
    std::vector<double> am1k;
    switch (atomicnr) {
      case 1:     //H
        am1k.push_back(1.1287500/au2eV);
        am1k.push_back(-1.0603290/au2eV);
        break;
      case 3:     //Li
        am1k.push_back(-0.4500000/au2eV);
        am1k.push_back(0.8000000/au2eV);
        break;
      case 4:     //Be
        am1k.push_back(1.6315720/au2eV);
        am1k.push_back(-2.1109590/au2eV);
        break;
      case 6:     //C
        am1k.push_back(0.0501070/au2eV);
        am1k.push_back(0.0507330/au2eV);
        break;
      case 7:     //N
        am1k.push_back(1.5016740/au2eV);
        am1k.push_back(-1.5057720/au2eV);
        break;
      case 8:     //O
        am1k.push_back(-1.1311280/au2eV);
        am1k.push_back(1.1378910/au2eV);
        break;
      case 9:     //F
        am1k.push_back(-0.0121660/au2eV);
        am1k.push_back(-0.0028520/au2eV);
        break;
      case 11:    //Na
        am1k.push_back(0.6433655/au2eV);
        am1k.push_back(1.0871788/au2eV);
        break;
      case 12:    //Mg
        am1k.push_back(2.1170500/au2eV);
        am1k.push_back(-2.5477670/au2eV);
        break;
      case 13:    //Al
        am1k.push_back(-0.4730900/au2eV);
        am1k.push_back(-0.1540510/au2eV);
        break;
      case 14:    //Si
        am1k.push_back(-0.3906000/au2eV);
        am1k.push_back(0.0572590/au2eV);
        break;
      case 15:    //P
        am1k.push_back(-0.6114210/au2eV);
        am1k.push_back(-0.0939350/au2eV);
        break;
      case 16:    //S
        am1k.push_back(-0.3991910/au2eV);
        am1k.push_back(-0.0548990/au2eV);
        break;
      case 17:    //Cl
        am1k.push_back(-0.1715910/au2eV);
        am1k.push_back(-0.0134580/au2eV);
        break;
      case 30:    //Zn
        am1k.push_back(-0.1112340/au2eV);
        am1k.push_back(-0.1323700/au2eV);
        break;
      case 31:    //Ga
        am1k.push_back(-0.5601790/au2eV);
        am1k.push_back(-0.2727310/au2eV);
        break;
      case 32:    //Ge
        am1k.push_back(0.9631726/au2eV);
        am1k.push_back(-0.9593891/au2eV);
        break;
      case 33:    //As
        am1k.push_back(-0.4600950/au2eV);
        am1k.push_back(-0.0889960/au2eV);
        break;
      case 34:    //Se
        am1k.push_back(0.0478730/au2eV);
        am1k.push_back(0.1147200/au2eV);
        break;
      case 35:    //Br
        am1k.push_back(0.9604580/au2eV);
        am1k.push_back(-0.9549160/au2eV);
        break;
      case 48:    //Cd
        am1k.push_back(0.0);
        break;
      case 49:    //In
        am1k.push_back(-0.3431380/au2eV);
        am1k.push_back(-0.1095320/au2eV);
        break;
      case 50:    //Sn
        am1k.push_back(-0.1503530/au2eV);
        am1k.push_back(-0.0444170/au2eV);
        break;
      case 51:    //Sb
        am1k.push_back(3.0020280/au2eV);
        am1k.push_back(-0.0188920/au2eV);
        break;
      case 52:    //Te
        am1k.push_back(0.0333910/au2eV);
        am1k.push_back(-1.9218670/au2eV);
        break;
      case 53:    //I
        am1k.push_back(-0.1314810/au2eV);
        am1k.push_back(-0.0368970/au2eV);
        break;
      case 80:    //Hg
        am1k.push_back(1.0827200/au2eV);
        am1k.push_back(-0.0965530/au2eV);
        break;
      case 81:    //Tl
        am1k.push_back(-1.3613990/au2eV);
        am1k.push_back(-0.0454010/au2eV);
        break;
      case 82:    //Pb
        am1k.push_back(-0.1225760/au2eV);
        am1k.push_back(-0.0566480/au2eV);
        break;
      case 83:    //Bi
        am1k.push_back(2.5816930/au2eV);
        am1k.push_back(0.0603200/au2eV);
        break;
    }
    return am1k;
  }
  std::vector<double> AM1L(int atomicnr, int atm2 = 0) {
    std::vector<double> am1l;
    switch (atomicnr) {
      case 1:     //H
        am1l.push_back(5.0962820);
        am1l.push_back(6.0037880);
        break;
      case 3:     //Li
        am1l.push_back(5.0000000);
        am1l.push_back(6.5000000);
        break;
      case 4:     //Be
        am1l.push_back(2.6729620);
        am1l.push_back(1.9685940);
        break;
      case 6:     //C
        am1l.push_back(6.0031650);
        am1l.push_back(6.0029790);
        break;
      case 7:     //N
        am1l.push_back(5.9011480);
        am1l.push_back(6.0046580);
        break;
      case 8:     //O
        am1l.push_back(6.0024770);
        am1l.push_back(5.9505120);
        break;
      case 9:     //F
        am1l.push_back(6.0235740);
        am1l.push_back(6.0037170);
        break;
      case 11:    //Na
        am1l.push_back(1.5465054);
        am1l.push_back(1.4529000);
        break;
      case 12:    //Mg
        am1l.push_back(6.0094770);
        am1l.push_back(4.3953700);
        break;
      case 13:    //Al
        am1l.push_back(1.9158250);
        am1l.push_back(6.0050860);
        break;
      case 14:    //Si
        am1l.push_back(6.0000540);
        am1l.push_back(6.0071830);
        break;
      case 15:    //P
        am1l.push_back(1.9972720);
        am1l.push_back(1.9983600);
        break;
      case 16:    //S
        am1l.push_back(6.0006690);
        am1l.push_back(6.0018450);
        break;
      case 17:    //Cl
        am1l.push_back(6.0008020);
        am1l.push_back(1.9666180);
        break;
      case 30:    //Zn
        am1l.push_back(6.0014780);
        am1l.push_back(1.9958390);
        break;
      case 31:    //Ga
        am1l.push_back(5.6232730);
        am1l.push_back(1.9918430);
        break;
      case 32:    //Ge
        am1l.push_back(6.0120134);
        am1l.push_back(5.7491802);
        break;
      case 33:    //As
        am1l.push_back(1.9831150);
        am1l.push_back(1.9929440);
        break;
      case 34:    //Se
        am1l.push_back(6.0074000);
        am1l.push_back(6.0086720);
        break;
      case 35:    //Br
        am1l.push_back(5.9765080);
        am1l.push_back(5.9447030);
        break;
      case 48:    //Cd
        am1l.push_back(0.0);
        break;
      case 49:    //In
        am1l.push_back(1.9940340);
        am1l.push_back(5.6832170);
        break;
      case 50:    //Sn
        am1l.push_back(6.0056940);
        am1l.push_back(2.2573810);
        break;
      case 51:    //Sb
        am1l.push_back(6.0053420);
        am1l.push_back(6.0114780);
        break;
      case 52:    //Te
        am1l.push_back(5.9563790);
        am1l.push_back(4.9732190);
        break;
      case 53:    //I
        am1l.push_back(5.2064170);
        am1l.push_back(6.0101170);
        break;
      case 80:    //Hg
        am1l.push_back(6.4965980);
        am1l.push_back(3.9262810);
        break;
      case 81:    //Tl
        am1l.push_back(3.5572260);
        am1l.push_back(2.3069950);
        break;
      case 82:    //Pb
        am1l.push_back(6.0030620);
        am1l.push_back(4.7437050);
        break;
      case 83:    //Bi
        am1l.push_back(5.0940220);
        am1l.push_back(6.0015380);
        break;
    }
    return am1l;
  }
  std::vector<double> AM1M(int atomicnr, int atm2 = 0) {
    std::vector<double> am1m;
    switch (atomicnr) {
      case 1:     //H
        am1m.push_back(1.5374650);
        am1m.push_back(1.5701890);
        break;
      case 3:     //Li
        am1m.push_back(1.0000000);
        am1m.push_back(1.0000000);
        break;
      case 4:     //Be
        am1m.push_back(1.7916860);
        am1m.push_back(1.7558710);
        break;
      case 6:     //C
        am1m.push_back(1.6422140);
        am1m.push_back(0.8924880);
        break;
      case 7:     //N
        am1m.push_back(1.7107400);
        am1m.push_back(1.7161490);
        break;
      case 8:     //O
        am1m.push_back(1.6073110);
        am1m.push_back(1.5983950);
        break;
      case 9:     //F
        am1m.push_back(1.8568590);
        am1m.push_back(2.6361580);
        break;
      case 11:    //Na
        am1m.push_back(0.9976699);
        am1m.push_back(1.4506099);
        break;
      case 12:    //Mg
        am1m.push_back(2.0844060);
        am1m.push_back(2.0636740);
        break;
      case 13:    //Al
        am1m.push_back(1.4517280);
        am1m.push_back(2.5199970);
        break;
      case 14:    //Si
        am1m.push_back(0.6322620);
        am1m.push_back(2.0199870);
        break;
      case 15:    //P
        am1m.push_back(0.7946240);
        am1m.push_back(1.9106770);
        break;
      case 16:    //S
        am1m.push_back(0.9621230);
        am1m.push_back(1.5799440);
        break;
      case 17:    //Cl
        am1m.push_back(1.0875020);
        am1m.push_back(2.2928910);
        break;
      case 30:    //Zn
        am1m.push_back(1.5160320);
        am1m.push_back(2.5196420);
        break;
      case 31:    //Ga
        am1m.push_back(1.5317800);
        am1m.push_back(2.1838640);
        break;
      case 32:    //Ge
        am1m.push_back(2.1633655);
        am1m.push_back(2.1693724);
        break;
      case 33:    //As
        am1m.push_back(1.0867930);
        am1m.push_back(2.1400580);
        break;
      case 34:    //Se
        am1m.push_back(2.0817170);
        am1m.push_back(1.5164230);
        break;
      case 35:    //Br
        am1m.push_back(2.3216540);
        am1m.push_back(2.3281420);
        break;
      case 48:    //Cd
        am1m.push_back(0.0);
        break;
      case 49:    //In
        am1m.push_back(1.6255160);
        am1m.push_back(2.8670090);
        break;
      case 50:    //Sn
        am1m.push_back(1.7046420);
        am1m.push_back(2.4698690);
        break;
      case 51:    //Sb
        am1m.push_back(0.8530600);
        am1m.push_back(2.7933110);
        break;
      case 52:    //Te
        am1m.push_back(2.2775750);
        am1m.push_back(0.5242430);
        break;
      case 53:    //I
        am1m.push_back(1.7488240);
        am1m.push_back(2.7103730);
        break;
      case 80:    //Hg
        am1m.push_back(1.1951460);
        am1m.push_back(2.6271600);
        break;
      case 81:    //Tl
        am1m.push_back(1.0928020);
        am1m.push_back(2.9650290);
        break;
      case 82:    //Pb
        am1m.push_back(1.9015970);
        am1m.push_back(2.8618790);
        break;
      case 83:    //Bi
        am1m.push_back(0.4997870);
        am1m.push_back(2.4279700);
        break;
    }
    return am1m;
  }
  double Dvalue(size_t atmnr, size_t idx) {
    //function returning the D values needed to calculate eris; values stored in Angstrom, but returned in a.u.
    //note that this function returns both D1 and D2; idx is then either 1 or 2
    double dval = 0.0;
    switch (atmnr) {
      case 1:                  //H
        dval = 0.0;
      case 3:                  //Li
        if (idx == 1) {dval = 1.0772805453024;}
        else if (idx == 2) {dval = 0.8641427827586;}
        break;
      case 4:                  //Be
        if (idx == 1) {dval = 0.5339679024905;}
        else if (idx == 2) {dval = 0.4295641497663;}
        break;
      case 6:                  //C
        if (idx == 1) {dval = 0.4409314053780;}
        else if (idx == 2) {dval = 0.3517837786516;}
        break;
      case 7:                  //N
        if (idx == 1) {dval = 0.3480401674091;}
        else if (idx == 2) {dval = 0.2801137638434;}
        break;
      case 8:                  //O
        if (idx == 1) {dval = 0.2162309620795;}
        else if (idx == 2) {dval = 0.2712423725348;}
        break;
      case 9:                  //F
        if (idx == 1) {dval = 0.1653838587473;}
        else if (idx == 2) {dval = 0.2601608726156;}
        break;
      case 11:                 //Na
        if (idx == 1) {dval = 0.91824824182;}
        else if (idx == 2) {dval = 0.74551702212;}
        break;
      case 12:                 //Mg
        if (idx == 1) {dval = 0.6034710424661;}
        else if (idx == 2) {dval = 0.5969065462793;}
        break;
      case 13:                 //Al
        if (idx == 1) {dval = 0.6404525387508;}
        else if (idx == 2) {dval = 0.8247568110748;}
        break;
      case 14:                 //Si
        if (idx == 1) {dval = 0.6955796273438;}
        else if (idx == 2) {dval = 0.6743514719617;}
        break;
      case 15:                 //P
        if (idx == 1) {dval = 0.5633063336025;}
        else if (idx == 2) {dval = 0.5884654818765;}
        break;
      case 16:                 //S
        if (idx == 1) {dval = 0.5934358846409;}
        else if (idx == 2) {dval = 0.5337539561451;}
        break;
      case 17:                 //Cl
        if (idx == 1) {dval = 0.4855653861897;}
        else if (idx == 2) {dval = 0.4116591214170;}
        break;
      case 30:                 //Zn
        if (idx == 1) {dval = 0.7940705126955;}
        else if (idx == 2) {dval = 0.7449319638157;}
        break;
      case 31:                 //Ga
        if (idx == 1) {dval = 0.5173602579027;}
        else if (idx == 2) {dval = 1.3373119811729;}
        break;
      case 32:                 //Ge
        if (idx == 1) {dval = 0.6307953192878;}
        else if (idx == 2) {dval = 0.7049308765449;}
        break;
      case 33:                 //As
        if (idx == 1) {dval = 0.5122252810265;}
        else if (idx == 2) {dval = 0.6588189567081;}
        break;
      case 34:                 //Se
        if (idx == 1) {dval = 0.4614326300290;}
        else if (idx == 2) {dval = 0.6479255792865;}
        break;
      case 35:                 //Br
        if (idx == 1) {dval = 0.1460013147146;}
        else if (idx == 2) {dval = 0.5276178289085;}
        break;
      case 48:                 //Cd
        if (idx == 1) {dval = 0.8457670512825;}
        else if (idx == 2) {dval = 0.6578943782897;}
        break;
      case 49:                 //In
        if (idx == 1) {dval = 0.8343135397859;}
        else if (idx == 2) {dval = 0.9405893627198;}
        break;
      case 50:                 //Sn
        if (idx == 1) {dval = 0.6942825081708;}
        else if (idx == 2) {dval = 0.8298458553693;}
        break;
      case 51:                 //Sb
        if (idx == 1) {dval = 0.7457113889258;}
        else if (idx == 2) {dval = 0.7155192362237;}
        break;
      case 52:                 //Te
        if (idx == 1) {dval = 0.1843747058104;}
        else if (idx == 2) {dval = 0.8251505189177;}
        break;
      case 53:                 //I
        if (idx == 1) {dval = 0.0836877350442;}
        else if (idx == 2) {dval = 0.5539057650855;}
        break;
      case 80:                 //Hg
        if (idx == 1) {dval = 0.6518304837420;}
        else if (idx == 2) {dval = 0.6436929024681;}
        break;
      case 81:                 //Tl
        if (idx == 1) {dval = 0.0413478961836;}
        else if (idx == 2) {dval = 0.8105465509115;}
        break;
      case 82:                 //Pb
        if (idx == 1) {dval = 0.5221015798537;}
        else if (idx == 2) {dval = 0.8435382097988;}
        break;
      case 83:                 //Bi
        if (idx == 1) {dval = 0.1480960097760;}
        else if (idx == 2) {dval = 0.8250028255589;}
        break;
    }
    return dval*dist_Angstrom2aum1;
  }
  double rho(size_t atmnr, size_t l) {
    //function returning the rho values needed to calculate eris; values stored in Angstrom but returned in atomic units
    double rho = 0.0;
    switch (atmnr) {
      case 1:                  //H
        if (l == 0) {rho = 0.486640184440;}
        else if (l == 1) {rho = 0.0;}
        else if (l == 2) {rho = 0.0;}
        break;
      case 3:                  //Li
        if (l == 0) {rho = 1.599878849930;}
        else if (l == 1) {rho = 2.287368480615;}
        else if (l == 2) {rho = 0.835697784494;}
        break;
      case 4:                  //Be
        if (l == 0) {rho = 0.798799045243;}
        else if (l == 1) {rho = 0.909553000943;}
        else if (l == 2) {rho = 0.749541089292;}
        break;
      case 6:                  //C
        if (l == 0) {rho = 0.642767927834;}
        else if (l == 1) {rho = 0.449532462965;}
        else if (l == 2) {rho = 0.345972966858;}
        break;
      case 7:                  //N
        if (l == 0) {rho = 0.604753079729;}
        else if (l == 1) {rho = 0.525917048522;}
        else if (l == 2) {rho = 0.359254597636;}
        break;
      case 8:                  //O
        if (l == 0) {rho = 0.456941201522;}
        else if (l == 1) {rho = 0.499258635325;}
        else if (l == 2) {rho = 0.323472582688;}
        break;
      case 9:                  //F
        if (l == 0) {rho = 0.685880274657;}
        else if (l == 1) {rho = 0.390911556297;}
        else if (l == 2) {rho = 0.432331000322;}
        break;
      case 11:                 //Na
        if (l == 0) {rho = 1.81994308925;}
        else if (l == 1) {rho = 1.24743397132;}
        else if (l == 2) {rho = 1.02719786193;}
        break;
      case 12:                 //Mg
        if (l == 0) {rho = 1.075460694417;}
        else if (l == 1) {rho = 0.981502387101;}
        else if (l == 2) {rho = 0.956048783543;}
        break;
      case 13:                 //Al
        if (l == 0) {rho = 1.246284086598;}
        else if (l == 1) {rho = 0.412222702317;}
        else if (l == 2) {rho = 1.169277713000;}
        break;
      case 14:                 //Si
        if (l == 0) {rho = 1.426426712705;}
        else if (l == 1) {rho = 0.864466649632;}
        else if (l == 2) {rho = 0.542475229082;}
        break;
      case 15:                 //P
        if (l == 0) {rho = 0.922816001025;}
        else if (l == 1) {rho = 0.613973592320;}
        else if (l == 2) {rho = 0.708874478410;}
        break;
      case 16:                 //S
        if (l == 0) {rho = 0.803092446274;}
        else if (l == 1) {rho = 0.396143029893;}
        else if (l == 2) {rho = 0.431103562107;}
        break;
      case 17:                 //Cl
        if (l == 0) {rho = 0.449583792797;}
        else if (l == 1) {rho = 0.388271699986;}
        else if (l == 2) {rho = 0.726154842182;}
        break;
      case 30:                 //Zn
        if (l == 0) {rho = 0.743960973132;}
        else if (l == 1) {rho = 1.113734180497;}
        else if (l == 2) {rho = 0.994294413835;}
        break;
      case 31:                 //Ga
        if (l == 0) {rho = 0.851144894368;}
        else if (l == 1) {rho = 0.515831612819;}
        else if (l == 2) {rho = 1.711209644194;}
        break;
      case 32:                 //Ge
        if (l == 0) {rho = 1.338944749466;}
        else if (l == 1) {rho = 0.696619077635;}
        else if (l == 2) {rho = 0.730772694640;}
        break;
      case 33:                 //As
        if (l == 0) {rho = 0.819143788072;}
        else if (l == 1) {rho = 0.524744273631;}
        else if (l == 2) {rho = 1.027840304780;}
        break;
      case 34:                 //Se
        if (l == 0) {rho = 0.968633392541;}
        else if (l == 1) {rho = 0.352329267283;}
        else if (l == 2) {rho = 0.500760359861;}
        break;
      case 35:                 //Br
        if (l == 0) {rho = 0.451562701487;}
        else if (l == 1) {rho = 0.391670767080;}
        else if (l == 2) {rho = 0.691966653807;}
        break;
      case 48:                 //Cd
        if (l == 0) {rho = 0.781957934859;}
        else if (l == 1) {rho = 0.741084349319;}
        else if (l == 2) {rho = 0.938063861121;}
        break;
      case 49:                 //In
        if (l == 0) {rho = 1.098331941968;}
        else if (l == 1) {rho = 0.583738678882;}
        else if (l == 2) {rho = 0.717078821775;}
        break;
      case 50:                 //Sn
        if (l == 0) {rho = 0.706519361502;}
        else if (l == 1) {rho = 0.822172786627;}
        else if (l == 2) {rho = 0.934107308885;}
        break;
      case 51:                 //Sb
        if (l == 0) {rho = 0.779307247171;}
        else if (l == 1) {rho = 0.576570118943;}
        else if (l == 2) {rho = 1.091774958213;}
        break;
      case 52:                 //Te
        if (l == 0) {rho = 0.702038451270;}
        else if (l == 1) {rho = 0.221214187240;}
        else if (l == 2) {rho = 1.211050437677;}
        break;
      case 53:                 //I
        if (l == 0) {rho = 0.528131297085;}
        else if (l == 1) {rho = 0.158444790902;}
        else if (l == 2) {rho = 0.513457003310;}
        break;
      case 80:                 //Hg
        if (l == 0) {rho = 1.086756136177;}
        else if (l == 1) {rho = 0.585960015150;}
        else if (l == 2) {rho = 1.010499581625;}
        break;
      case 81:                 //Tl
        if (l == 0) {rho = 0.688257458272;}
        else if (l == 1) {rho = 0.102785527809;}
        else if (l == 2) {rho = 1.195467888842;}
        break;
      case 82:                 //Pb
        if (l == 0) {rho = 1.026734684572;}
        else if (l == 1) {rho = 0.584380347173;}
        else if (l == 2) {rho = 1.230544509879;}
        break;
      case 83:                 //Bi
        if (l == 0) {rho = 1.442927491968;}
        else if (l == 1) {rho = 0.390478300666;}
        else if (l == 2) {rho = 1.022952090658;}
        break;
    }
    return rho*dist_Angstrom2aum1;
  }
  double eri1Center(int atmnr, int Lbra, int Lket) {
    //function that gives back the semi-empirical 1-center eris
    //values stored in eV but returned in a.u.
    //Lbra is the sum of azimuthal quantum numbers for bra (ss = 0; pp = 2; sp = 1)
    //Lket is the sum of azimuthal quantum numbers for ket (ss = 0; pp = 2; sp = 1; p*p* = -2)
    double eri = 0.0;
    switch (atmnr) {
      case 1:       //H
        eri = 14.7942080;                                                                        //(ss|ss)
        break;
      case 3:       //Li
        if ((Lbra == 0)&&(Lket == 0)) {eri = 4.5000000;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.2500000;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 3.0000000;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.5000000;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.15000003;}                                   //(sp|sp)||(ps|ps)
        break;
      case 4:       //Be
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.0128510;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.0571820;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.5761990;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.0052190;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.5446790;}                                    //(sp|sp)||(ps|ps)
        break;
      case 5:       //B
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.59;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.86;}                                           //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 9.56;}             //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.86;}           //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.81;}                                           //(sp|sp)||(ps|ps)
        break;
      case 6:       //C
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.2007080;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 10.7962920;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 10.2650270;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.0425660;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.2909800;}                                    //(sp|sp)||(ps|ps)
        break;
      case 7:       //N
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.9047870;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.7546720;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.3485650;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 10.8072770;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.1367130;}                                    //(sp|sp)||(ps|ps)
        break;
      case 8:       //O
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.7557600;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 13.6540160;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 10.6211600;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 12.4060950;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.5938830;}                                    //(sp|sp)||(ps|ps)
        break;
      case 9:       //F
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.4966670;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 14.8172560;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 16.0736890;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 14.4183930;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.7277630;}                                    //(sp|sp)||(ps|ps)
        break;
      case 11:      //Na
        if ((Lbra == 0)&&(Lket == 0)) {eri = 3.9558692;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.3363963;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.1929109;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.0588074;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.5687889;}                                    //(sp|sp)||(ps|ps)
        break;
      case 12:      //Mg
        if ((Lbra == 0)&&(Lket == 0)) {eri = 6.6943000;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.9104460;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.7939950;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.0908230;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.5433000;}                                    //(sp|sp)||(ps|ps)
        break;
      case 13:      //Al
        if ((Lbra == 0)&&(Lket == 0)) {eri = 5.7767370;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.3477900;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 11.6598560;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.1210770;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 4.0062450;}                                    //(sp|sp)||(ps|ps)
        break;
      case 14:      //Si
        if ((Lbra == 0)&&(Lket == 0)) {eri = 5.0471960;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.7593670;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.9490570;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.1612970;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.9198320;}                                    //(sp|sp)||(ps|ps)
        break;
      case 15:      //P
        if ((Lbra == 0)&&(Lket == 0)) {eri = 7.8016150;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.6184780;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.1869490;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.0620020;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.5428090;}                                    //(sp|sp)||(ps|ps)
        break;
      case 16:      //S
        if ((Lbra == 0)&&(Lket == 0)) {eri = 8.9646670;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 9.9681640;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.7859360;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.9702470;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 4.0418360;}                                    //(sp|sp)||(ps|ps)
        break;
      case 17:      //Cl
        if ((Lbra == 0)&&(Lket == 0)) {eri = 16.0136010;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.5222150;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.0481150;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.5041540;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.4811530;}                                    //(sp|sp)||(ps|ps)
        break;
      case 30:      //Zn
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.6771960;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 4.9801740;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.7362040;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.6696560;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.6004130;}                                    //(sp|sp)||(ps|ps)
        break;
      case 31:      //Ga
        if ((Lbra == 0)&&(Lket == 0)) {eri = 8.4585540;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.0868550;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.9256190;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.9830450;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.0512600;}                                    //(sp|sp)||(ps|ps)
        break;
      case 32:      //Ge
        if ((Lbra == 0)&&(Lket == 0)) {eri = 5.3769635;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.6718647;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 10.2095293;}   //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.9242663;}  //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.3370204;}                                  //(sp|sp)||(ps|ps)
        break;
      case 33:      //As
        if ((Lbra == 0)&&(Lket == 0)) {eri = 8.7890010;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.2872500;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.3979830;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 8.2103460;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.9510340;}                                    //(sp|sp)||(ps|ps)
        break;
      case 34:      //Se
        if ((Lbra == 0)&&(Lket == 0)) {eri = 7.4325910;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 9.5683260;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 10.0604610;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.7242890;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 4.0165580;}                                    //(sp|sp)||(ps|ps)
        break;
      case 35:      //Br
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.9434250;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.2827630;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 16.0616800;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.8168490;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.5788690;}                                    //(sp|sp)||(ps|ps)
        break;
      case 48:      //Cd
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.2069600;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 4.9481040;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.2315390;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.6696560;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.6562340;}                                    //(sp|sp)||(ps|ps)
        break;
      case 49:      //In
        if ((Lbra == 0)&&(Lket == 0)) {eri = 6.5549000;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.2992690;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.2298730;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.9842110;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.6314610;}                                    //(sp|sp)||(ps|ps)
        break;
      case 50:      //Sn
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.1900330;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.6738100;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.2353270;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.1822140;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.0331570;}                                    //(sp|sp)||(ps|ps)
        break;
      case 51:      //Sb
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.2382770;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.3500000;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.2776800;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.2500000;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.4244640;}                                    //(sp|sp)||(ps|ps)
        break;
      case 52:      //Te
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.2550730;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.7775920;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.1691450;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.7551210;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.7724620;}                                    //(sp|sp)||(ps|ps)
        break;
      case 53:      //I
        if ((Lbra == 0)&&(Lket == 0)) {eri = 13.6319430;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.2883300;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 14.9904060;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.9664070;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.6300350;}                                    //(sp|sp)||(ps|ps)
        break;
      case 80:      //Hg
        if ((Lbra == 0)&&(Lket == 0)) {eri = 6.6247200;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 14.7092830;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 10.6392970;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 16.0007400;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.0363110;}                                    //(sp|sp)||(ps|ps)v
        break;
      case 81:      //Tl
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.4604120;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 4.9927850;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 11.2238830;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 8.9627270;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.5304060;}                                    //(sp|sp)||(ps|ps)
        break;
      case 82:      //Pb
        if ((Lbra == 0)&&(Lket == 0)) {eri = 7.0119920;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.1837800;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.7937820;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.0456510;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.5663020;}                                    //(sp|sp)||(ps|ps)
        break;
      case 83:      //Bi
        if ((Lbra == 0)&&(Lket == 0)) {eri = 4.9894800;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.6960070;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.1033080;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 8.3354470;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.5991220;}                                    //(sp|sp)||(ps|ps)
        break;
    }
    return eri/au2eV;
  }
  //charge model information; this is separate calculation
  double CM1ck(size_t atomicnr) {
    //function returning the atomic parameters ck for CM1
    double cm1ck = 0.0;
    if ((atomicnr == 1)||((atomicnr > 5)&&(atomicnr < 10))||(atomicnr == 14)||(atomicnr == 16)||(atomicnr == 17)||(atomicnr == 35)||(atomicnr == 53)) {
      if (atomicnr == 9) {cm1ck = 0.3381;}                 //F
      else if (atomicnr == 16) {cm1ck = -0.0834;}          //S
      else if (atomicnr == 17) {cm1ck = -0.1080;}          //Cl
      else if (atomicnr == 35) {cm1ck = -0.0116;}          //Br
      else if (atomicnr == 53) {cm1ck = -0.3213;}          //I
    }
    else {std::cout << "ERROR: MNDO.hpp: PM3: CM1ck(): partial charge calculation involving non-parametrized atom";}
    return cm1ck;
  }
  double CM1dk(size_t atomicnr) {
    //function returning the atomic parameters dk for CM1
    double cm1dk = 0.0;
    if ((atomicnr == 1)||((atomicnr > 5)&&(atomicnr < 10))||(atomicnr == 14)||(atomicnr == 16)||(atomicnr == 17)||(atomicnr == 35)||(atomicnr == 53)) {
      if (atomicnr == 7) {cm1dk = -0.0909;}                //N
      else if (atomicnr == 8) {cm1dk = -0.0449;}           //O
      else if (atomicnr == 9) {cm1dk = 0.0148;}            //F
      else if (atomicnr == 16) {cm1dk = -0.0848;}          //S
      else if (atomicnr == 17) {cm1dk = -0.1168;}          //Cl
      else if (atomicnr == 35) {cm1dk = -0.0338;}          //Br
      else if (atomicnr == 53) {cm1dk = -0.0636;}          //I
    }
    else {std::cout << "ERROR: MNDO.hpp: PM3: CM1dk(): partial charge calculation involving non-parametrized atom";}
    return cm1dk;
  }
  double CM1ckkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters ckk' for CM1; implemented the function just in case a refit is ever going to take place
    return 0.0;
  }
  double CM1dkkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters dkk' for CM1
    double cm1dkkp = 0.0;
    bool valid1 = false;
    bool valid2 = false;
    if ((atom1 == 1)||((atom1 > 5)&&(atom1 < 10))||(atom1 == 14)||(atom1 == 16)||(atom1 == 17)||(atom1 == 35)||(atom1 == 53)) {valid1 = true;}
    if ((atom2 == 1)||((atom2 > 5)&&(atom2 < 10))||(atom2 == 14)||(atom2 == 16)||(atom2 == 17)||(atom2 == 35)||(atom2 == 53)) {valid2 = true;}
    if ((valid1)&&(valid2)) {
      size_t atomA = std::min(atom1,atom2);
      size_t atomB = std::max(atom1,atom2);
      if ((atomA == 1)&&(atomB == 7)) {cm1dkkp = 0.1854;}           //H-N
      else if ((atomA == 1)&&(atomB == 8)) {cm1dkkp = 0.1434;}      //H-O
      else if ((atomA == 1)&&(atomB == 14)) {cm1dkkp = -0.1004;}    //H-Si
    }
    else {std::cout << "ERROR: MNDO.hpp: PM3: CM1dkkp(): partial charge calculation involving non-parametrized atom";}
    return cm1dkkp;
  }
  double CM2ckkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters ckk' for CM2
    double cm2ckkp = 0.0;
    size_t atomA = std::min(atom1,atom2);
    size_t atomB = std::max(atom1,atom2);
    if ((atomA == 1)&&(atomB == 6)) {cm2ckkp = 0.003;}            //H-C
    else if ((atomA == 1)&&(atomB == 7)) {cm2ckkp = 0.274;}       //H-N
    else if ((atomA == 1)&&(atomB == 8)) {cm2ckkp = 0.185;}       //H-O
    else if ((atomA == 1)&&(atomB == 14)) {cm2ckkp = -0.021;}     //H-Si
    else if ((atomA == 1)&&(atomB == 16)) {cm2ckkp = 0.089;}      //H-S
    else if ((atomA == 6)&&(atomB == 7)) {cm2ckkp = -0.022;}      //C-N
    else if ((atomA == 6)&&(atomB == 8)) {cm2ckkp = 0.025;}       //C-O
    else if ((atomA == 6)&&(atomB == 14)) {cm2ckkp = -0.107;}     //C-Si
    else if ((atomA == 6)&&(atomB == 16)) {cm2ckkp = -0.033;}     //C-S
    else if ((atomA == 7)&&(atomB == 8)) {cm2ckkp = -0.030;}      //N-O
    if (atom1 != atomA) {cm2ckkp *= -1.0;}
    return cm2ckkp;
  }
  double CM2dkkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters dkk' for CM2
    double cm2dkkp = 0.0;
    bool valid1 = false;
    bool valid2 = false;
    if ((atom1 == 1)||((atom1 > 5)&&(atom1 < 10))||((atom1 > 13)&&(atom1 < 18))||(atom1 == 35)||(atom1 == 53)) {valid1 = true;}
    if ((atom2 == 1)||((atom2 > 5)&&(atom2 < 10))||((atom2 > 13)&&(atom2 < 18))||(atom2 == 35)||(atom2 == 53)) {valid2 = true;}
    if ((valid1)&&(valid2)) {
      size_t atomA = std::min(atom1,atom2);
      size_t atomB = std::max(atom1,atom2);
      if ((atomA == 1)&&(atomB == 15)) {cm2dkkp = 0.253;}           //H-P
      else if ((atomA == 6)&&(atomB == 7)) {cm2dkkp = 0.156;}       //C-N
      else if ((atomA == 6)&&(atomB == 8)) {cm2dkkp = 0.016;}       //C-O
      else if ((atomA == 6)&&(atomB == 9)) {cm2dkkp = 0.025;}       //C-F
      else if ((atomA == 6)&&(atomB == 15)) {cm2dkkp = 0.082;}      //C-P
      else if ((atomA == 6)&&(atomB == 16)) {cm2dkkp = 0.112;}      //C-S
      else if ((atomA == 6)&&(atomB == 17)) {cm2dkkp = 0.117;}      //C-Cl
      else if ((atomA == 6)&&(atomB == 35)) {cm2dkkp = 0.040;}      //C-Br
      else if ((atomA == 6)&&(atomB == 53)) {cm2dkkp = -0.032;}     //C-I
      else if ((atomA == 7)&&(atomB == 8)) {cm2dkkp = -0.043;}      //N-O
      else if ((atomA == 8)&&(atomB == 15)) {cm2dkkp = 0.181;}      //O-P
      else if ((atomA == 8)&&(atomB == 16)) {cm2dkkp = 0.056;}      //O-S
      else if ((atomA == 9)&&(atomB == 15)) {cm2dkkp = 0.244;}      //F-P
      else if ((atomA == 15)&&(atomB == 16)) {cm2dkkp = -0.087;}    //P-S
      if (atom1 != atomA) {cm2dkkp *= -1.0;}
    }
    else {std::cout << "ERROR: MNDO.hpp: PM3: CM2dkkp(): partial charge calculation involving non-parametrized atom";}
    return cm2dkkp;
  }
  double CM3ckkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters ckk' for CM3
    double cm3ckkp = 0.0;
    size_t atomA = std::min(atom1,atom2);
    size_t atomB = std::max(atom1,atom2);
    if ((atomA == 6)&&(atomB == 7)) {cm3ckkp = -0.017;}           //C-N
    else if ((atomA == 6)&&(atomB == 8)) {cm3ckkp = 0.010;}       //C-O
    else if ((atomA == 7)&&(atomB == 8)) {cm3ckkp = 0.274;}       //N-O
    else if ((atomA == 8)&&(atomB == 14)) {cm3ckkp = -0.044;}     //O-Si
    else if ((atomA == 8)&&(atomB == 15)) {cm3ckkp = -0.076;}     //O-P
    else if ((atomA == 15)&&(atomB == 16)) {cm3ckkp = 0.245;}     //P-S
    if (atom1 != atomA) {cm3ckkp *= -1.0;}
    return cm3ckkp;
  }
  double CM3dkkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters dkk' for CM3
    double cm3dkkp = 0.0;
    bool valid1 = false;
    bool valid2 = false;
    if ((atom1 == 1)||(atom1 == 3)||((atom1 > 5)&&(atom1 < 10))||((atom1 > 13)&&(atom1 < 18))||(atom1 == 35)) {valid1 = true;}
    if ((atom2 == 1)||(atom2 == 3)||((atom2 > 5)&&(atom2 < 10))||((atom2 > 13)&&(atom2 < 18))||(atom2 == 35)) {valid2 = true;}
    if ((valid1)&&(valid2)) {
      size_t atomA = std::min(atom1,atom2);
      size_t atomB = std::max(atom1,atom2);
      if ((atomA == 1)&&(atomB == 6)) {cm3dkkp = 0.021;}            //H-C
      else if ((atomA == 1)&&(atomB == 7)) {cm3dkkp = 0.243;}       //H-N
      else if ((atomA == 1)&&(atomB == 8)) {cm3dkkp = 0.153;}       //H-O
      else if ((atomA == 1)&&(atomB == 14)) {cm3dkkp = 0.191;}      //H-Si
      else if ((atomA == 1)&&(atomB == 15)) {cm3dkkp = 0.233;}      //H-P
      else if ((atomA == 1)&&(atomB == 16)) {cm3dkkp = 0.118;}      //H-S
      else if ((atomA == 3)&&(atomB == 6)) {cm3dkkp = 0.192;}       //Li-C
      else if ((atomA == 3)&&(atomB == 7)) {cm3dkkp = 0.403;}       //Li-N
      else if ((atomA == 3)&&(atomB == 8)) {cm3dkkp = 0.390;}       //Li-O
      else if ((atomA == 3)&&(atomB == 9)) {cm3dkkp = 0.430;}       //Li-F
      else if ((atomA == 3)&&(atomB == 16)) {cm3dkkp = 0.268;}      //Li-S
      else if ((atomA == 3)&&(atomB == 17)) {cm3dkkp = 0.117;}      //Li-Cl
      else if ((atomA == 6)&&(atomB == 7)) {cm3dkkp = 0.126;}       //C-N
      else if ((atomA == 6)&&(atomB == 8)) {cm3dkkp = 0.018;}       //C-O
      else if ((atomA == 6)&&(atomB == 9)) {cm3dkkp = 0.019;}       //C-F
      else if ((atomA == 6)&&(atomB == 14)) {cm3dkkp = 0.069;}      //C-Si
      else if ((atomA == 6)&&(atomB == 15)) {cm3dkkp = 0.197;}      //C-P
      else if ((atomA == 6)&&(atomB == 16)) {cm3dkkp = 0.076;}      //C-S
      else if ((atomA == 6)&&(atomB == 17)) {cm3dkkp = 0.097;}      //C-Cl
      else if ((atomA == 6)&&(atomB == 35)) {cm3dkkp = 0.028;}      //C-Br
      else if ((atomA == 7)&&(atomB == 8)) {cm3dkkp = -0.085;}      //N-O
      else if ((atomA == 7)&&(atomB == 15)) {cm3dkkp = 0.018;}      //N-P
      else if ((atomA == 8)&&(atomB == 14)) {cm3dkkp = 0.141;}      //O-Si
      else if ((atomA == 8)&&(atomB == 15)) {cm3dkkp = 0.231;}      //O-P
      else if ((atomA == 8)&&(atomB == 16)) {cm3dkkp = 0.125;}      //O-S
      else if ((atomA == 9)&&(atomB == 14)) {cm3dkkp = 0.099;}      //F-Si
      else if ((atomA == 9)&&(atomB == 15)) {cm3dkkp = 0.175;}      //F-P
      else if ((atomA == 14)&&(atomB == 17)) {cm3dkkp = -0.180;}    //Si-Cl
      else if ((atomA == 15)&&(atomB == 16)) {cm3dkkp = -0.534;}    //P-S
      else if ((atomA == 15)&&(atomB == 17)) {cm3dkkp = -0.394;}    //P-Cl
      if (atom1 != atomA) {cm3dkkp *= -1.0;}
    }
    else {std::cout << "WARNING: MNDO.hpp: PM3: CM3dkkp(): partial charge calculation involving non-parametrized atom";}
    return cm3dkkp;
  }
  double B0() {return 0.460;}
};
class PM3PDDG: public MNDO {
  //this is the implementation of Jorgensen's PM3PDDG
  //M. P. Repasky, J. Chandrasekhar, W. L. Jorgensen, J. Comput. Chem., 23, 1601, 2002
  //I. Tubert-Brohman, C. R. W. Guimaraes, M. P. Repasky, W. L. Jorgensen, J. Comput. Chem., 25, 138, 2004
  //I. Tubert-Brohman, C. R. W. Guimaraes, W. L. Jorgensen, J. Chem. Theory Comput., 1, 817, 2005
  //bX quantities are "barred" tensors, which are used only for the open-shell case
public:
  PM3PDDG(BSet _bset, Molecule _mol, std::string _openclosed = "0", std::string _corecorrection = "0"): MNDO(_bset,_mol,_openclosed,_corecorrection) {
    doPDDG = true;
  }
  ~PM3PDDG() {}
  std::string Type() {return "PM3PDDG";}
  double gfactor(size_t iatm1, size_t iatm2, double RAB) {return 1.0;}
  virtual double gfactor_dR(size_t iatm1, size_t iatm2) {return 0.0;}
  double AM1factor(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {
    //RAB in a.u.
    double AM1factor = 0.0;
    RAB *= dist_Angstrom2au;
    for (size_t idx = 0; idx < K1.size(); ++idx) {
      AM1factor += K1[idx]*exp(-L1[idx]*(RAB - M1[idx])*(RAB - M1[idx]))/RAB;
    }
    return AM1factor;
  }
  virtual double AM1factor_dR(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {
    //RAB in a.u.
    double am1factor_dr = 0.0;
    double RR = RAB*dist_Angstrom2au;
    for (size_t idx = 0; idx < K1.size(); ++idx) {
      am1factor_dr += K1[idx]*(2.0*RR*RR*L1[idx] - 2.0*L1[idx]*M1[idx]*RR + 1.0)*exp(-L1[idx]*(RR - M1[idx])*(RR - M1[idx]))/(RAB*RR);
    }
    return am1factor_dr;  
  }
  virtual double AM1factor_dR2(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {
    //RAB in a.u.
    double am1factor_dr2 = 0.0;
    double RR = RAB*dist_Angstrom2au;
    for (size_t idx = 0; idx < K1.size(); ++idx) {
      am1factor_dr2 += 2.0*K1[idx]*exp(-L1[idx]*(RR - M1[idx])*(RR - M1[idx]))*(L1[idx]*M1[idx]*RR - 1.0 - RR*L1[idx]*(RR - M1[idx])*(2.0*RR*RR*L1[idx] - 2.0*L1[idx]*M1[idx]*RR + 1.0))/(RAB*RAB*RR);
    }
    return am1factor_dr2;  
  }
  //other auxiliary functions
  void checkAtoms() {
    //function that checks whether the atoms in the molecule are fully defined in the respective MNDO theories
    bool def;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      def = false;
      if (atoms[iatm] == 1) {def = true;}                                   //H
      else if ((atoms[iatm] > 5)&&(atoms[iatm] < 10)) {def = true;}         //C,N,O,F
      else if ((atoms[iatm] > 13)&&(atoms[iatm] < 18)) {def = true;}        //Si,P,S,Cl
      else if ((atoms[iatm] == 35)||(atoms[iatm] == 53)) {def = true;}      //Br,I
      if (!def) {throw("ERROR: MNDO.hpp: PM3PDDG: checkAtoms(): atom not fully specified for PM3PDDG-theory");}
    }
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapPM3PDDG(atm);}
  double ElementFormTheo(size_t atomicnr) {
    //function returning theoretical heats of formation for elements; values directly in a.u.
    double enth = 0.0;
    switch (atomicnr) {
      case 1:      //H
        enth = -13.120566;
        break;
      case 6:      //C
        enth = -113.428242;
        break;
      case 7:      //N
        enth = -158.416205;
        break;
      case 8:      //O
        enth = -292.188766;
        break;
      case 9:      //F
        enth = -442.457133;
        break;
      case 14:     //Si
        enth = -66.839000;
        break;
      case 15:     //P
        enth = -117.212854;
        break;
      case 16:     //S
        enth = -166.336554;
        break;
      case 17:     //Cl
        enth = -305.715201;
        break;
      case 35:     //Br
        enth = -351.013887;
        break;
      case 53:     //I
        enth = -291.537869;
        break;
    }
    return enth/au2eV;
  }
  double betaA0(size_t atomicnr, int L) {
    //function returning the bonding parameters; values stored in eV, returned however in a.u.
    double betaa0 = 0.0;
    switch (atomicnr) {
      case 1:         //H
        betaa0 = -6.152654;
        break;
      case 6:         //C
        if (L == 0) {betaa0 = -11.952818;}
        else if (L == 1) {betaa0 = -9.922411;}
        break;
      case 7:         //N
        if (L == 0) {betaa0 = -14.117230;}
        else if (L == 1) {betaa0 = -19.938509;}
        break;
      case 8:         //O
        if (L == 0) {betaa0 = -44.874553;}
        else if (L == 1) {betaa0 = -24.601939;}
        break;
      case 9:         //F
        if (L == 0) {betaa0 = -50.937301;}
        else if (L == 1) {betaa0 = -31.636976;}
        break;
      case 14:        //Si
        if (L == 0) {betaa0 = -3.376445;}
        else if (L == 1) {betaa0 = -3.620969;}
        break;
      case 15:        //P
        if (L == 0) {betaa0 = -12.676297;}
        else if (L == 1) {betaa0 = -7.093318;}
        break;
      case 16:        //S
        if (L == 0) {betaa0 = -2.953912;}
        else if (L == 1) {betaa0 = -8.507779;}
        break;
      case 17:        //Cl
        if (L == 0) {betaa0 = -26.913129;}
        else if (L == 1) {betaa0 = -14.991178;}
        break;
      case 35:        //Br
        if (L == 0) {betaa0 = -21.538044;}
        else if (L == 1) {betaa0 = -8.524764;}
        break;
      case 53:        //I
        if (L == 0) {betaa0 = -16.592621;}
        else if (L == 1) {betaa0 = -6.599816;}
        break;
    }
    return betaa0/au2eV;
  }
  double UlX(size_t atomicnr, int L) {
    //function that returns the atomic potential U
    double ulx = 0.0;
    switch (atomicnr) {
      case 1:                 //H
        ulx = -12.893272;
        break;
      case 6:                 //C
        if (L == 0) {ulx = -48.241241;}
        else if (L == 1) {ulx = -36.461256;}
        break;
      case 7:                 //N
        if (L == 0) {ulx = -49.454546;}
        else if (L == 1) {ulx = -47.757406;}
        break;
      case 8:                 //O
        if (L == 0) {ulx = -87.412505;}
        else if (L == 1) {ulx = -72.183070;}
        break;
      case 9:                 //F
        if (L == 0) {ulx = -111.400432;}
        else if (L == 1) {ulx = -106.395264;}
        break;
      case 14:                //Si
        if (L == 0) {ulx = -26.332522;}
        else if (L == 1) {ulx = -22.602540;}
        break;
      case 15:                //P
        if (L == 0) {ulx = -37.882113;}
        else if (L == 1) {ulx = -30.312979;}
        break;
      case 16:                //S
        if (L == 0) {ulx = -43.906366;}
        else if (L == 1) {ulx = -43.461348;}
        break;
      case 17:                //Cl
        if (L == 0) {ulx = -95.094434;}
        else if (L == 1) {ulx = -53.921651;}
        break;
      case 35:                //Br
        if (L == 0) {ulx = -115.841963;}
        else if (L == 1) {ulx = -74.205146;}
        break;
      case 53:                //I
        if (L == 0) {ulx = -97.664174;}
        else if (L == 1) {ulx = -61.167137;}
        break;
    }
    return ulx/au2eV;
  }
  double alpha(int atomicnr, int atm2) {
    //function that returns alpha values for PM3PDDG
    //alphas stored and returned in 1/Angstrom
    double alpha = 0.0;
    switch (atomicnr) {
      case 1:     //H
        alpha = 3.381686;
        break;
      case 6:     //C
        alpha = 2.725772;
        break;
      case 7:     //N
        alpha = 2.849124;
        break;
      case 8:     //O
        alpha = 3.225309;
        break;
      case 9:     //F
        alpha = 3.200571;
        break;
      case 14:    //Si
        alpha = 2.215157;
        break;
      case 15:    //P
        alpha = 2.005294;
        break;
      case 16:    //S
        alpha = 2.539751;
        break;
      case 17:    //Cl
        alpha = 2.497617;
        break;
      case 35:    //Br
        alpha = 2.424673;
        break;
      case 53:    //I
        alpha = 1.978170;
        break;
    }
    return alpha;
  }
  std::vector<double> AM1K(int atomicnr, int atm2 = 0) {
    std::vector<double> am1k;
    switch (atomicnr) {
      case 1:     //H
        am1k.push_back(1.122244/au2eV);
        am1k.push_back(-1.069737/au2eV);
        break;
      case 6:     //C
        am1k.push_back(0.048906/au2eV);
        am1k.push_back(0.047697/au2eV);
        break;
      case 7:     //N
        am1k.push_back(1.513320/au2eV);
        am1k.push_back(-1.511892/au2eV);
        break;
      case 8:     //O
        am1k.push_back(-1.138455/au2eV);
        am1k.push_back(1.146007/au2eV);
        break;
      case 9:     //F
        am1k.push_back(-0.008079/au2eV);
        am1k.push_back(-0.002659/au2eV);
        break;
      case 14:    //Si
        am1k.push_back(-0.071314/au2eV);
        am1k.push_back(0.089451/au2eV);
        break;
      case 15:    //P
        am1k.push_back(-0.398055/au2eV);
        am1k.push_back(-0.079653/au2eV);
        break;
      case 16:    //S
        am1k.push_back(-0.330692/au2eV);
        am1k.push_back(0.024171/au2eV);
        break;
      case 17:    //Cl
        am1k.push_back(-0.112222/au2eV);
        am1k.push_back(-0.013061/au2eV);
        break;
      case 35:    //Br
        am1k.push_back(0.961362/au2eV);
        am1k.push_back(-0.948834/au2eV);
        break;
      case 53:    //I
        am1k.push_back(-0.136003/au2eV);
        am1k.push_back(-0.037287/au2eV);
        break;
    }
    return am1k;
  }
  std::vector<double> AM1L(int atomicnr, int atm2 = 0) {
    std::vector<double> am1l;
    switch (atomicnr) {
      case 1:     //H
        am1l.push_back(4.707790);
        am1l.push_back(5.857995);
        break;
      case 6:     //C
        am1l.push_back(5.765340);
        am1l.push_back(5.973721);
        break;
      case 7:     //N
        am1l.push_back(5.904394);
        am1l.push_back(6.030014);
        break;
      case 8:     //O
        am1l.push_back(6.000043);
        am1l.push_back(5.963494);
        break;
      case 9:     //F
        am1l.push_back(5.938969);
        am1l.push_back(5.925105);
        break;
      case 14:    //Si
        am1l.push_back(6.000000);
        am1l.push_back(6.000000);
        break;
      case 15:    //P
        am1l.push_back(1.997272);
        am1l.push_back(1.998360);
        break;
      case 16:    //S
        am1l.push_back(6.000000);
        am1l.push_back(6.000000);
        break;
      case 17:    //Cl
        am1l.push_back(5.963719);
        am1l.push_back(1.999556);
        break;
      case 35:    //Br
        am1l.push_back(6.013600);
        am1l.push_back(5.976329);
        break;
      case 53:    //I
        am1l.push_back(3.852912);
        am1l.push_back(5.229264);
        break;
    }
    return am1l;
  }
  std::vector<double> AM1M(int atomicnr, int atm2 = 0) {
    std::vector<double> am1m;
    switch (atomicnr) {
      case 1:     //H
        am1m.push_back(1.547099);
        am1m.push_back(1.567893);
        break;
      case 6:     //C
        am1m.push_back(1.682232);
        am1m.push_back(0.894406);
        break;
      case 7:     //N
        am1m.push_back(1.728376);
        am1m.push_back(1.734108);
        break;
      case 8:     //O
        am1m.push_back(1.622362);
        am1m.push_back(1.614788);
        break;
      case 9:     //F
        am1m.push_back(1.863949);
        am1m.push_back(2.388864);
        break;
      case 14:    //Si
        am1m.push_back(0.237995);
        am1m.push_back(1.897728);
        break;
      case 15:    //P
        am1m.push_back(0.950073);
        am1m.push_back(2.336959);
        break;
      case 16:    //S
        am1m.push_back(0.823837);
        am1m.push_back(2.017756);
        break;
      case 17:    //Cl
        am1m.push_back(1.027719);
        am1m.push_back(2.286377);
        break;
      case 35:    //Br
        am1m.push_back(2.340445);
        am1m.push_back(2.348745);
        break;
      case 53:    //I
        am1m.push_back(1.697455);
        am1m.push_back(2.768669);
        break;
    }
    return am1m;
  }
  double Dvalue(size_t atmnr, size_t idx) {
    //function returning the D values needed to calculate eris; values stored and returned in a.u.
    //note that this function returns both D1 and D2; idx is then either 1 or 2
    double dval = 0.0;
    switch (atmnr) {
      case 1:                  //H
        dval = 0.0;
      case 6:                  //C
        if (idx == 1) {dval = 0.831413;}
        else if (idx == 2) {dval = 0.663222;}
        break;
      case 7:                  //N
        if (idx == 1) {dval = 0.654855;}
        else if (idx == 2) {dval = 0.526924;}
        break;
      case 8:                  //O
        if (idx == 1) {dval = 0.403741;}
        else if (idx == 2) {dval = 0.528360;}
        break;
      case 9:                  //F
        if (idx == 1) {dval = 0.246601;}
        else if (idx == 2) {dval = 0.482551;}
        break;
      case 14:                 //Si
        if (idx == 1) {dval = 1.310515;}
        else if (idx == 2) {dval = 1.126089;}
        break;
      case 15:                 //P
        if (idx == 1) {dval = 0.893978;}
        else if (idx == 2) {dval = 0.960457;}
        break;
      case 16:                 //S
        if (idx == 1) {dval = 1.006989;}
        else if (idx == 2) {dval = 0.891487;}
        break;
      case 17:                 //Cl
        if (idx == 1) {dval = 0.827561;}
        else if (idx == 2) {dval = 0.732427;}
        break;
      case 35:                 //Br
        if (idx == 1) {dval = 0.473860;}
        else if (idx == 2) {dval = 0.968214;}
        break;
      case 53:                 //I
        if (idx == 1) {dval = 0.407261;}
        else if (idx == 2) {dval = 1.062574;}
        break;
    }
    return dval;
  }
  double rho(size_t atmnr, size_t l) {
    //function returning the rho values needed to calculate eris; values stored and returned in atomic units
    double rho = 0.0;
    switch (atmnr) {
      case 1:                  //H
        if (l == 0) {rho = 0.919616;}
        else if (l == 1) {rho = 0.0;}
        else if (l == 2) {rho = 0.0;}
        break;
      case 6:                  //C
        if (l == 0) {rho = 1.214657;}
        else if (l == 1) {rho = 0.848467;}
        else if (l == 2) {rho = 0.652785;}
        break;
      case 7:                  //N
        if (l == 0) {rho = 1.142818;}
        else if (l == 1) {rho = 0.991235;}
        else if (l == 2) {rho = 0.676704;}
        break;
      case 8:                  //O
        if (l == 0) {rho = 0.863494;}
        else if (l == 1) {rho = 0.936266;}
        else if (l == 2) {rho = 0.624291;}
        break;
      case 9:                  //F
        if (l == 0) {rho = 1.296126;}
        else if (l == 1) {rho = 0.634633;}
        else if (l == 2) {rho = 0.805802;}
        break;
      case 14:                 //Si
        if (l == 0) {rho = 2.695556;}
        else if (l == 1) {rho = 1.630757;}
        else if (l == 2) {rho = 0.949200;}
        break;
      case 15:                 //P
        if (l == 0) {rho = 1.743870;}
        else if (l == 1) {rho = 1.050851;}
        else if (l == 2) {rho = 1.208907;}
        break;
      case 16:                 //S
        if (l == 0) {rho = 1.517625;}
        else if (l == 1) {rho = 0.711672;}
        else if (l == 2) {rho = 0.754336;}
        break;
      case 17:                 //Cl
        if (l == 0) {rho = 0.849590;}
        else if (l == 1) {rho = 0.696164;}
        else if (l == 2) {rho = 2.299104;}
        break;
      case 35:                 //Br
        if (l == 0) {rho = 0.853330;}
        else if (l == 1) {rho = 1.046430;}
        else if (l == 2) {rho = 1.280643;}
        break;
      case 53:                 //I
        if (l == 0) {rho = 0.998024;}
        else if (l == 1) {rho = 0.532290;}
        else if (l == 2) {rho = 0.979783;}
        break;
    }
    return rho;
  }
  double eri1Center(int atmnr, int Lbra, int Lket) {
    //function that gives back the semi-empirical 1-center eris
    //values stored in eV but returned in a.u.
    //Lbra is the sum of azimuthal quantum numbers for bra (ss = 0; pp = 2; sp = 1)
    //Lket is the sum of azimuthal quantum numbers for ket (ss = 0; pp = 2; sp = 1; p*p* = -2)
    double eri = 0.0;
    switch (atmnr) {
      case 1:       //H
        eri = 14.7942080;                                                                        //(ss|ss)
        break;
      case 6:       //C
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.2007080;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 10.7962920;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 10.2650270;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.0425660;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.2909800;}                                    //(sp|sp)||(ps|ps)
        break;
      case 7:       //N
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.9047870;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.7546720;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.3485650;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 10.8072770;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.1367130;}                                    //(sp|sp)||(ps|ps)
        break;
      case 8:       //O
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.7557600;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 13.6540160;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 10.6211600;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 12.4060950;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.5938830;}                                    //(sp|sp)||(ps|ps)
        break;
      case 9:       //F
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.4966670;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 14.8172560;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 16.0736890;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 14.4183930;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.7277630;}                                    //(sp|sp)||(ps|ps)
        break;
      case 14:      //Si
        if ((Lbra == 0)&&(Lket == 0)) {eri = 5.0471960;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.7593670;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.9490570;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.1612970;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.9198320;}                                    //(sp|sp)||(ps|ps)
        break;
      case 15:      //P
        if ((Lbra == 0)&&(Lket == 0)) {eri = 7.8016150;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.6184780;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.1869490;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.0620020;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.5428090;}                                    //(sp|sp)||(ps|ps)
        break;
      case 16:      //S
        if ((Lbra == 0)&&(Lket == 0)) {eri = 8.9646670;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 9.9681640;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.7859360;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.9702470;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 4.0418360;}                                    //(sp|sp)||(ps|ps)
        break;
      case 17:      //Cl
        if ((Lbra == 0)&&(Lket == 0)) {eri = 16.0136010;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.5222150;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.0481150;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.5041540;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.4811530;}                                    //(sp|sp)||(ps|ps)
        break;
      case 35:      //Br
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.9434250;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.2827630;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 16.0616800;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.8168490;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.5788690;}                                    //(sp|sp)||(ps|ps)
        break;
      case 53:      //I
        if ((Lbra == 0)&&(Lket == 0)) {eri = 13.6319430;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.2883300;}                                    //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 14.9904060;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.9664070;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.6300350;}                                    //(sp|sp)||(ps|ps)
        break;
    }
    return eri/au2eV;
  }
  double PAPDDG(int atmnr, int index) {
    //function returning the PDDG PAi parameters in Hartree
    double pa = 0.0;
    switch (atmnr) {
      case 1:      //H
        if (index == 1) {pa = 0.057193;}
        else if (index == 2) {pa = -0.034823;}
        break;
      case 6:      //C
        if (index == 1) {pa = -0.000743;}
        else if (index == 2) {pa = 0.000985;}
        break;
      case 7:      //N
        if (index == 1) {pa = -0.003160;}
        else if (index == 2) {pa = 0.012501;}
        break;
      case 8:      //O
        if (index == 1) {pa = -0.001000;}
        else if (index == 2) {pa = -0.001522;}
        break;
      case 9:      //F
        if (index == 1) {pa = -0.012866;}
        else if (index == 2) {pa = 0.007315;}
        break;
      case 14:     //Si
        if (index == 1) {pa = -0.091928;}
        else if (index == 2) {pa = -0.040753;}
        break;
      case 15:     //P
        if (index == 1) {pa = 0.462741;}
        else if (index == 2) {pa = -0.020444;}
        break;
      case 16:     //S
        if (index == 1) {pa = 0.120434;}
        else if (index == 2) {pa = -0.002663;}
        break;
      case 17:     //Cl
        if (index == 1) {pa = -0.016552;}
        else if (index == 2) {pa = -0.016646;}
        break;
      case 35:     //Br
        if (index == 1) {pa = -0.013772;}
        else if (index == 2) {pa = 0.008849;}
        break;
      case 53:     //I
        if (index == 1) {pa = 0.012901;}
        else if (index == 2) {pa = -0.012825;}
        break;
    }
    return pa/au2eV;
  }
  double DAPDDG(int atmnr, int index) {
    //function returning the PDDG DAi parameters in inverse Angstroem
    double da = 0.0;
    switch (atmnr) {
      case 1:      //H
        if (index == 1) {da = 0.663395;}
        else if (index == 2) {da = 1.081901;}
        break;
      case 6:      //C
        if (index == 1) {da = 0.836915;}
        else if (index == 2) {da = 1.585236;}
        break;
      case 7:      //N
        if (index == 1) {da = 1.004172;}
        else if (index == 2) {da = 1.516336;}
        break;
      case 8:      //O
        if (index == 1) {da = 1.360685;}
        else if (index == 2) {da = 1.366407;}
        break;
      case 9:      //F
        if (index == 1) {da = 1.305681;}
        else if (index == 2) {da = 1.842572;}
        break;
      case 14:     //Si
        if (index == 1) {da = 1.163190;}
        else if (index == 2) {da = 2.190526;}
        break;
      case 15:     //P
        if (index == 1) {da = 0.714296;}
        else if (index == 2) {da = 2.041209;}
        break;
      case 16:     //S
        if (index == 1) {da = 0.672870;}
        else if (index == 2) {da = 2.032340;}
        break;
      case 17:     //Cl
        if (index == 1) {da = 1.727690;}
        else if (index == 2) {da = 1.784655;}
        break;
      case 35:     //Br
        if (index == 1) {da = 1.852030;}
        else if (index == 2) {da = 2.338958;}
        break;
      case 53:     //I
        if (index == 1) {da = 1.994299;}
        else if (index == 2) {da = 2.263417;}
        break;
    }
    return da;
  }
  //charge model information; this is separate calculation
  double CM1ck(size_t atomicnr) {
    //function returning the atomic parameters ck for CM1
    double cm1ck = 0.0;
    if ((atomicnr == 1)||((atomicnr > 5)&&(atomicnr < 10))||(atomicnr == 14)||(atomicnr == 16)||(atomicnr == 17)||(atomicnr == 35)||(atomicnr == 53)) {
      if (atomicnr == 9) {cm1ck = 0.3381;}                 //F
      else if (atomicnr == 16) {cm1ck = -0.0834;}          //S
      else if (atomicnr == 17) {cm1ck = -0.1080;}          //Cl
      else if (atomicnr == 35) {cm1ck = -0.0116;}          //Br
      else if (atomicnr == 53) {cm1ck = -0.3213;}          //I
    }
    else {std::cout << "ERROR: MNDO.hpp: PM3PDDG: CM1ck(): partial charge calculation involving non-parametrized atom";}
    return cm1ck;
  }
  double CM1dk(size_t atomicnr) {
    //function returning the atomic parameters dk for CM1
    double cm1dk = 0.0;
    if ((atomicnr == 1)||((atomicnr > 5)&&(atomicnr < 10))||(atomicnr == 14)||(atomicnr == 16)||(atomicnr == 17)||(atomicnr == 35)||(atomicnr == 53)) {
      if (atomicnr == 7) {cm1dk = -0.0909;}                //N
      else if (atomicnr == 8) {cm1dk = -0.0449;}           //O
      else if (atomicnr == 9) {cm1dk = 0.0148;}            //F
      else if (atomicnr == 16) {cm1dk = -0.0848;}          //S
      else if (atomicnr == 17) {cm1dk = -0.1168;}          //Cl
      else if (atomicnr == 35) {cm1dk = -0.0338;}          //Br
      else if (atomicnr == 53) {cm1dk = -0.0636;}          //I
    }
    else {std::cout << "ERROR: MNDO.hpp: PM3PDDG: CM1dk(): partial charge calculation involving non-parametrized atom";}
    return cm1dk;
  }
  double CM1ckkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters ckk' for CM1; implemented the function just in case a refit is ever going to take place
    return 0.0;
  }
  double CM1dkkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters dkk' for CM1
    double cm1dkkp = 0.0;
    bool valid1 = false;
    bool valid2 = false;
    if ((atom1 == 1)||((atom1 > 5)&&(atom1 < 10))||(atom1 == 14)||(atom1 == 16)||(atom1 == 17)||(atom1 == 35)||(atom1 == 53)) {valid1 = true;}
    if ((atom2 == 1)||((atom2 > 5)&&(atom2 < 10))||(atom2 == 14)||(atom2 == 16)||(atom2 == 17)||(atom2 == 35)||(atom2 == 53)) {valid2 = true;}
    if ((valid1)&&(valid2)) {
      size_t atomA = std::min(atom1,atom2);
      size_t atomB = std::max(atom1,atom2);
      if ((atomA == 1)&&(atomB == 7)) {cm1dkkp = 0.1854;}           //H-N
      else if ((atomA == 1)&&(atomB == 8)) {cm1dkkp = 0.1434;}      //H-O
      else if ((atomA == 1)&&(atomB == 14)) {cm1dkkp = -0.1004;}    //H-Si
    }
    else {std::cout << "ERROR: MNDO.hpp: PM3PDDG: CM1dkkp(): partial charge calculation involving non-parametrized atom";}
    return cm1dkkp;
  }
  double CM2ckkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters ckk' for CM2
    double cm2ckkp = 0.0;
    size_t atomA = std::min(atom1,atom2);
    size_t atomB = std::max(atom1,atom2);
    if ((atomA == 1)&&(atomB == 6)) {cm2ckkp = 0.003;}            //H-C
    else if ((atomA == 1)&&(atomB == 7)) {cm2ckkp = 0.274;}       //H-N
    else if ((atomA == 1)&&(atomB == 8)) {cm2ckkp = 0.185;}       //H-O
    else if ((atomA == 1)&&(atomB == 14)) {cm2ckkp = -0.021;}     //H-Si
    else if ((atomA == 1)&&(atomB == 16)) {cm2ckkp = 0.089;}      //H-S
    else if ((atomA == 6)&&(atomB == 7)) {cm2ckkp = -0.022;}      //C-N
    else if ((atomA == 6)&&(atomB == 8)) {cm2ckkp = 0.025;}       //C-O
    else if ((atomA == 6)&&(atomB == 14)) {cm2ckkp = -0.107;}     //C-Si
    else if ((atomA == 6)&&(atomB == 16)) {cm2ckkp = -0.033;}     //C-S
    else if ((atomA == 7)&&(atomB == 8)) {cm2ckkp = -0.030;}      //N-O
    if (atom1 != atomA) {cm2ckkp *= -1.0;}
    return cm2ckkp;
  }
  double CM2dkkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters dkk' for CM2
    double cm2dkkp = 0.0;
    bool valid1 = false;
    bool valid2 = false;
    if ((atom1 == 1)||((atom1 > 5)&&(atom1 < 10))||((atom1 > 13)&&(atom1 < 18))||(atom1 == 35)||(atom1 == 53)) {valid1 = true;}
    if ((atom2 == 1)||((atom2 > 5)&&(atom2 < 10))||((atom2 > 13)&&(atom2 < 18))||(atom2 == 35)||(atom2 == 53)) {valid2 = true;}
    if ((valid1)&&(valid2)) {
      size_t atomA = std::min(atom1,atom2);
      size_t atomB = std::max(atom1,atom2);
      if ((atomA == 1)&&(atomB == 15)) {cm2dkkp = 0.253;}           //H-P
      else if ((atomA == 6)&&(atomB == 7)) {cm2dkkp = 0.156;}       //C-N
      else if ((atomA == 6)&&(atomB == 8)) {cm2dkkp = 0.016;}       //C-O
      else if ((atomA == 6)&&(atomB == 9)) {cm2dkkp = 0.025;}       //C-F
      else if ((atomA == 6)&&(atomB == 15)) {cm2dkkp = 0.082;}      //C-P
      else if ((atomA == 6)&&(atomB == 16)) {cm2dkkp = 0.112;}      //C-S
      else if ((atomA == 6)&&(atomB == 17)) {cm2dkkp = 0.117;}      //C-Cl
      else if ((atomA == 6)&&(atomB == 35)) {cm2dkkp = 0.040;}      //C-Br
      else if ((atomA == 6)&&(atomB == 53)) {cm2dkkp = -0.032;}     //C-I
      else if ((atomA == 7)&&(atomB == 8)) {cm2dkkp = -0.043;}      //N-O
      else if ((atomA == 8)&&(atomB == 15)) {cm2dkkp = 0.181;}      //O-P
      else if ((atomA == 8)&&(atomB == 16)) {cm2dkkp = 0.056;}      //O-S
      else if ((atomA == 9)&&(atomB == 15)) {cm2dkkp = 0.244;}      //F-P
      else if ((atomA == 15)&&(atomB == 16)) {cm2dkkp = -0.087;}    //P-S
      if (atom1 != atomA) {cm2dkkp *= -1.0;}
    }
    else {std::cout << "ERROR: MNDO.hpp: PM3PDDG: CM2dkkp(): partial charge calculation involving non-parametrized atom";}
    return cm2dkkp;
  }
  double CM3ckkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters ckk' for CM3
    double cm3ckkp = 0.0;
    size_t atomA = std::min(atom1,atom2);
    size_t atomB = std::max(atom1,atom2);
    if ((atomA == 6)&&(atomB == 7)) {cm3ckkp = -0.017;}           //C-N
    else if ((atomA == 6)&&(atomB == 8)) {cm3ckkp = 0.010;}       //C-O
    else if ((atomA == 7)&&(atomB == 8)) {cm3ckkp = 0.274;}       //N-O
    else if ((atomA == 8)&&(atomB == 14)) {cm3ckkp = -0.044;}     //O-Si
    else if ((atomA == 8)&&(atomB == 15)) {cm3ckkp = -0.076;}     //O-P
    else if ((atomA == 15)&&(atomB == 16)) {cm3ckkp = 0.245;}     //P-S
    if (atom1 != atomA) {cm3ckkp *= -1.0;}
    return cm3ckkp;
  }
  double CM3dkkp(size_t atom1, size_t atom2) {
    //function returning the atom-pair parameters dkk' for CM3
    double cm3dkkp = 0.0;
    bool valid1 = false;
    bool valid2 = false;
    if ((atom1 == 1)||(atom1 == 3)||((atom1 > 5)&&(atom1 < 10))||((atom1 > 13)&&(atom1 < 18))||(atom1 == 35)) {valid1 = true;}
    if ((atom2 == 1)||(atom2 == 3)||((atom2 > 5)&&(atom2 < 10))||((atom2 > 13)&&(atom2 < 18))||(atom2 == 35)) {valid2 = true;}
    if ((valid1)&&(valid2)) {
      size_t atomA = std::min(atom1,atom2);
      size_t atomB = std::max(atom1,atom2);
      if ((atomA == 1)&&(atomB == 6)) {cm3dkkp = 0.021;}            //H-C
      else if ((atomA == 1)&&(atomB == 7)) {cm3dkkp = 0.243;}       //H-N
      else if ((atomA == 1)&&(atomB == 8)) {cm3dkkp = 0.153;}       //H-O
      else if ((atomA == 1)&&(atomB == 14)) {cm3dkkp = 0.191;}      //H-Si
      else if ((atomA == 1)&&(atomB == 15)) {cm3dkkp = 0.233;}      //H-P
      else if ((atomA == 1)&&(atomB == 16)) {cm3dkkp = 0.118;}      //H-S
      else if ((atomA == 3)&&(atomB == 6)) {cm3dkkp = 0.192;}       //Li-C
      else if ((atomA == 3)&&(atomB == 7)) {cm3dkkp = 0.403;}       //Li-N
      else if ((atomA == 3)&&(atomB == 8)) {cm3dkkp = 0.390;}       //Li-O
      else if ((atomA == 3)&&(atomB == 9)) {cm3dkkp = 0.430;}       //Li-F
      else if ((atomA == 3)&&(atomB == 16)) {cm3dkkp = 0.268;}      //Li-S
      else if ((atomA == 3)&&(atomB == 17)) {cm3dkkp = 0.117;}      //Li-Cl
      else if ((atomA == 6)&&(atomB == 7)) {cm3dkkp = 0.126;}       //C-N
      else if ((atomA == 6)&&(atomB == 8)) {cm3dkkp = 0.018;}       //C-O
      else if ((atomA == 6)&&(atomB == 9)) {cm3dkkp = 0.019;}       //C-F
      else if ((atomA == 6)&&(atomB == 14)) {cm3dkkp = 0.069;}      //C-Si
      else if ((atomA == 6)&&(atomB == 15)) {cm3dkkp = 0.197;}      //C-P
      else if ((atomA == 6)&&(atomB == 16)) {cm3dkkp = 0.076;}      //C-S
      else if ((atomA == 6)&&(atomB == 17)) {cm3dkkp = 0.097;}      //C-Cl
      else if ((atomA == 6)&&(atomB == 35)) {cm3dkkp = 0.028;}      //C-Br
      else if ((atomA == 7)&&(atomB == 8)) {cm3dkkp = -0.085;}      //N-O
      else if ((atomA == 7)&&(atomB == 15)) {cm3dkkp = 0.018;}      //N-P
      else if ((atomA == 8)&&(atomB == 14)) {cm3dkkp = 0.141;}      //O-Si
      else if ((atomA == 8)&&(atomB == 15)) {cm3dkkp = 0.231;}      //O-P
      else if ((atomA == 8)&&(atomB == 16)) {cm3dkkp = 0.125;}      //O-S
      else if ((atomA == 9)&&(atomB == 14)) {cm3dkkp = 0.099;}      //F-Si
      else if ((atomA == 9)&&(atomB == 15)) {cm3dkkp = 0.175;}      //F-P
      else if ((atomA == 14)&&(atomB == 17)) {cm3dkkp = -0.180;}    //Si-Cl
      else if ((atomA == 15)&&(atomB == 16)) {cm3dkkp = -0.534;}    //P-S
      else if ((atomA == 15)&&(atomB == 17)) {cm3dkkp = -0.394;}    //P-Cl
      if (atom1 != atomA) {cm3dkkp *= -1.0;}
    }
    else {std::cout << "WARNING: MNDO.hpp: PM3PDDG: CM3dkkp(): partial charge calculation involving non-parametrized atom";}
    return cm3dkkp;
  }
};
class MNDOPDDG: public MNDO {
  //this is the implementation of Jorgensen's MNDOPDDG
  //M. P. Repasky, J. Chandrasekhar, W. L. Jorgensen, J. Comput. Chem., 23, 1601, 2002
  //I. Tubert-Brohman, C. R. W. Guimaraes, M. P. Repasky, W. L. Jorgensen, J. Comput. Chem., 25, 138, 2004
  //I. Tubert-Brohman, C. R. W. Guimaraes, W. L. Jorgensen, J. Chem. Theory Comput., 1, 817, 2005
  //bX quantities are "barred" tensors, which are used only for the open-shell case
public:
  MNDOPDDG(BSet _bset, Molecule _mol, std::string _openclosed = "0", std::string _corecorrection = "0"): MNDO(_bset,_mol,_openclosed,_corecorrection) {doPDDG = true;}
  ~MNDOPDDG() {}
  std::string Type() {return "MNDOPDDG";}
  double gfactor(size_t iatm1, size_t iatm2, double RAB) {return 1.0;}
  virtual double gfactor_dR(size_t iatm1, size_t iatm2) {return 0.0;}
  //other auxiliary functions
  void checkAtoms() {
    //function that checks whether the atoms in the molecule are fully defined in the respective MNDO theories
    bool def;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      def = false;
      if (atoms[iatm] == 1) {def = true;}                                                        //H
      else if ((atoms[iatm] > 5)&&(atoms[iatm] < 10)) {def = true;}                              //C,N,O,F
      else if ((atoms[iatm] == 17)||(atoms[iatm] == 35)||(atoms[iatm] == 53)) {def = true;}      //Cl,Br,I
      if (!def) {throw("ERROR: MNDO.hpp: MNDOPDDG: checkAtoms(): atom not fully specified for MNDOPDDG-theory");}
    }
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapMNDOPDDG(atm);}
  double ElementFormTheo(size_t atomicnr) {
    //function returning theoretical heats of formation for elements; values directly in a.u.
    double enth = 0.0;
    switch (atomicnr) {
      case 1:      //H
        enth = -12.015956;
        break;
      case 6:      //C
        enth = -123.864412;
        break;
      case 7:      //N
        enth = -206.466626;
        break;
      case 8:      //O
        enth = -310.879745;
        break;
      case 9:      //F
        enth = -488.703243;
        break;
      case 17:     //Cl
        enth = -378.909727;
        break;
      case 35:     //Br
        enth = -349.564096;
        break;
      case 53:     //I
        enth = -356.076398;
        break;
    }
    return enth/au2eV;
  }
  double betaA0(size_t atomicnr, int L) {
    //function returning the bonding parameters; values stored in eV, returned however in a.u.
    double betaa0 = 0.0;
    switch (atomicnr) {
      case 1:         //H
        betaa0 = -7.493504;
        break;
      case 6:         //C
        if (L == 0) {betaa0 = -18.841334;}
        else if (L == 1) {betaa0 = -7.922234;}
        break;
      case 7:         //N
        if (L == 0) {betaa0 = -20.375774;}
        else if (L == 1) {betaa0 = -21.085373;}
        break;
      case 8:         //O
        if (L == 0) {betaa0 = -33.606336;}
        else if (L == 1) {betaa0 = -27.984442;}
        break;
      case 9:         //F
        if (L == 0) {betaa0 = -67.827612;}
        else if (L == 1) {betaa0 = -40.924818;}
        break;
      case 17:        //Cl
        if (L == 0) {betaa0 = -15.663317;}
        else if (L == 1) {betaa0 = -15.399331;}
        break;
      case 35:        //Br
        if (L == 0) {betaa0 = -7.054170;}
        else if (L == 1) {betaa0 = -10.221030;}
        break;
      case 53:        //I
        if (L == 0) {betaa0 = -6.698375;}
        else if (L == 1) {betaa0 = -5.693814;}
        break;
    }
    return betaa0/au2eV;
  }
  double UlX(size_t atomicnr, int L) {
    //function that returns the atomic potential U
    double ulx = 0.0;
    switch (atomicnr) {
      case 1:                 //H
        ulx = -11.724114;
        break;
      case 6:                 //C
        if (L == 0) {ulx = -53.837582;}
        else if (L == 1) {ulx = -39.936409;}
        break;
      case 7:                 //N
        if (L == 0) {ulx = -71.871894;}
        else if (L == 1) {ulx = -58.216617;}
        break;
      case 8:                 //O
        if (L == 0) {ulx = -97.884970;}
        else if (L == 1) {ulx = -77.342674;}
        break;
      case 9:                 //F
        if (L == 0) {ulx = -134.220379;}
        else if (L == 1) {ulx = -107.155961;}
        break;
      case 17:                //Cl
        if (L == 0) {ulx = -111.133653;}
        else if (L == 1) {ulx = -78.062493;}
        break;
      case 35:                //Br
        if (L == 0) {ulx = -100.637007;}
        else if (L == 1) {ulx = -76.015735;}
        break;
      case 53:                //I
        if (L == 0) {ulx = -106.588422;}
        else if (L == 1) {ulx = -75.282605;}
        break;
    }
    return ulx/au2eV;
  }
  double alpha(int atomicnr, int atm2) {
    //function that returns alpha values for PM3PDDG
    //alphas stored and returned in 1/Angstrom
    double alpha = 0.0;
    switch (atomicnr) {
      case 1:     //H
        alpha = 2.491813;
        break;
      case 6:     //C
        alpha = 2.555522;
        break;
      case 7:     //N
        alpha = 2.843678;
        break;
      case 8:     //O
        alpha = 3.238842;
        break;
      case 9:     //F
        alpha = 3.322382;
        break;
      case 17:    //Cl
        alpha = 2.602846;
        break;
      case 35:    //Br
        alpha = 2.414265;
        break;
      case 53:    //I
        alpha = 2.242446;
        break;
    }
    return alpha;
  }
  double Dvalue(size_t atmnr, size_t idx) {
    //function returning the D values needed to calculate eris; values stored and returned in a.u.
    //note that this function returns both D1 and D2; idx is then either 1 or 2
    double dval = 0.0;
    switch (atmnr) {
      case 1:                  //H
        dval = 0.0;
      case 6:                  //C
        if (idx == 1) {dval = 0.794158;}
        else if (idx == 2) {dval = 0.671090;}
        break;
      case 7:                  //N
        if (idx == 1) {dval = 0.643624;}
        else if (idx == 2) {dval = 0.543495;}
        break;
      case 8:                  //O
        if (idx == 1) {dval = 0.547344;}
        else if (idx == 2) {dval = 0.454088;}
        break;
      case 9:                  //F
        if (idx == 1) {dval = 0.361556;}
        else if (idx == 2) {dval = 0.421593;}
        break;
      case 17:                 //Cl
        if (idx == 1) {dval = 0.411609;}
        else if (idx == 2) {dval = 0.821202;}
        break;
      case 35:                 //Br
        if (idx == 1) {dval = 0.574623;}
        else if (idx == 2) {dval = 0.944892;}
        break;
      case 53:                 //I
        if (idx == 1) {dval = 1.209529;}
        else if (idx == 2) {dval = 1.043559;}
        break;
    }
    return dval;
  }
  double rho(size_t atmnr, size_t l) {
    //function returning the rho values needed to calculate eris; values stored and returned in atomic units
    double rho = 0.0;
    switch (atmnr) {
      case 1:                  //H
        if (l == 0) {rho = 1.058920;}
        else if (l == 1) {rho = 0.0;}
        else if (l == 2) {rho = 0.0;}
        break;
      case 6:                  //C
        if (l == 0) {rho = 1.112429;}
        else if (l == 1) {rho = 0.805697;}
        else if (l == 2) {rho = 0.737353;}
        break;
      case 7:                  //N
        if (l == 0) {rho = 1.001103;}
        else if (l == 1) {rho = 0.639479;}
        else if (l == 2) {rho = 0.615679;}
        break;
      case 8:                  //O
        if (l == 0) {rho = 0.882296;}
        else if (l == 1) {rho = 0.527927;}
        else if (l == 2) {rho = 0.526913;}
        break;
      case 9:                  //F
        if (l == 0) {rho = 0.804078;}
        else if (l == 1) {rho = 0.383553;}
        else if (l == 2) {rho = 0.476913;}
        break;
      case 17:                 //Cl
        if (l == 0) {rho = 0.905190;}
        else if (l == 1) {rho = 0.554248;}
        else if (l == 2) {rho = 0.825604;}
        break;
      case 35:                 //Br
        if (l == 0) {rho = 0.904802;}
        else if (l == 1) {rho = 0.668867;}
        else if (l == 2) {rho = 0.884978;}
        break;
      case 53:                 //I
        if (l == 0) {rho = 0.904562;}
        else if (l == 1) {rho = 1.002309;}
        else if (l == 2) {rho = 0.992001;}
        break;
    }
    return rho;
  }
  double PAPDDG(int atmnr, int index) {
    //function returning the PDDG PAi parameters in Hartree
    double pa = 0.0;
    switch (atmnr) {
      case 1:      //H
        if (index == 1) {pa = -0.108861;}
        else if (index == 2) {pa = -0.024706;}
        break;
      case 6:      //C
        if (index == 1) {pa = -0.006889;}
        else if (index == 2) {pa = -0.027751;}
        break;
      case 7:      //N
        if (index == 1) {pa = 0.035027;}
        else if (index == 2) {pa = -0.001721;}
        break;
      case 8:      //O
        if (index == 1) {pa = 0.086344;}
        else if (index == 2) {pa = 0.030403;}
        break;
      case 9:      //F
        if (index == 1) {pa = -0.011579;}
        else if (index == 2) {pa = -0.012943;}
        break;
      case 17:     //Cl
        if (index == 1) {pa = -0.017119;}
        else if (index == 2) {pa = 0.005497;}
        break;
      case 35:     //Br
        if (index == 1) {pa = -0.017133;}
        else if (index == 2) {pa = -0.016964;}
        break;
      case 53:     //I
        if (index == 1) {pa = 0.009616;}
        else if (index == 2) {pa = -0.007505;}
        break;
    }
    return pa/au2eV;
  }
  double DAPDDG(int atmnr, int index) {
    //function returning the PDDG DAi parameters in inverse Angstroem
    double da = 0.0;
    switch (atmnr) {
      case 1:      //H
        if (index == 1) {da = 0.460721;}
        else if (index == 2) {da = 1.298731;}
        break;
      case 6:      //C
        if (index == 1) {da = 1.192456;}
        else if (index == 2) {da = 1.329522;}
        break;
      case 7:      //N
        if (index == 1) {da = 1.011630;}
        else if (index == 2) {da = 2.278423;}
        break;
      case 8:      //O
        if (index == 1) {da = 0.725408;}
        else if (index == 2) {da = 0.709728;}
        break;
      case 9:      //F
        if (index == 1) {da = 0.834606;}
        else if (index == 2) {da = 1.875603;}
        break;
      case 17:     //Cl
        if (index == 1) {da = 1.466335;}
        else if (index == 2) {da = 2.236842;}
        break;
      case 35:     //Br
        if (index == 1) {da = 2.201539;}
        else if (index == 2) {da = 2.255764;}
        break;
      case 53:     //I
        if (index == 1) {da = 2.572332;}
        else if (index == 2) {da = 2.936456;}
        break;
    }
    return da;
  }
};
class PM3BP: public MNDO {
  //this is the implementation PM3BP for DNA/RNA base pairs
  //T. J. Giese, E. C. Sherer, C. J. Cramer, D. M. York, J. Chem. Theory Comput., 1, 1275, 2005
  //bX quantities are "barred" tensors, which are used only for the open-shell case
public:
  PM3BP(BSet _bset, Molecule _mol, std::string _openclosed = "0", std::string _corecorrection = "0"): MNDO(_bset,_mol,_openclosed,_corecorrection) {}
  ~PM3BP() {}
  std::string Type() {return "PM3BP";}
  double gfactor(size_t iatm1, size_t iatm2, double RAB) {return 1.0;}
  virtual double gfactor_dR(size_t iatm1, size_t iatm2) {return 0.0;}
  double AM1factor(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {
    //RAB in a.u.
    double AM1factor = 0.0;
    RAB *= dist_Angstrom2au;
    for (size_t idx = 0; idx < K1.size(); ++idx) {
      AM1factor += K1[idx]*exp(-L1[idx]*(RAB - M1[idx])*(RAB - M1[idx]))/RAB;
    }
    return AM1factor;
  }
  virtual double AM1factor_dR(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {
    //RAB in a.u.
    double am1factor_dr = 0.0;
    double RR = RAB*dist_Angstrom2au;
    for (size_t idx = 0; idx < K1.size(); ++idx) {
      am1factor_dr += K1[idx]*(2.0*RR*RR*L1[idx] - 2.0*L1[idx]*M1[idx]*RR + 1.0)*exp(-L1[idx]*(RR - M1[idx])*(RR - M1[idx]))/(RAB*RR);
    }
    return am1factor_dr;  
  }
  virtual double AM1factor_dR2(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {
    //RAB in a.u.
    double am1factor_dr2 = 0.0;
    double RR = RAB*dist_Angstrom2au;
    for (size_t idx = 0; idx < K1.size(); ++idx) {
      am1factor_dr2 += 2.0*K1[idx]*exp(-L1[idx]*(RR - M1[idx])*(RR - M1[idx]))*(L1[idx]*M1[idx]*RR - 1.0 - RR*L1[idx]*(RR - M1[idx])*(2.0*RR*RR*L1[idx] - 2.0*L1[idx]*M1[idx]*RR + 1.0))/(RAB*RAB*RR);
    }
    return am1factor_dr2;  
  }
  //other auxiliary functions
  void checkAtoms() {
    //function that checks whether the atoms in the molecule are fully defined in PM3BP
    bool def;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      def = false;
      if (atoms[iatm] == 1) {def = true;}                            //H
      else if ((atoms[iatm] > 5)&&(atoms[iatm] < 9)) {def = true;}   //C,N,O
      if (!def) {throw("ERROR: MNDO.hpp: PM3BP: checkAtoms(): atom not fully specified for PM3BP-theory");}
    }
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapPM3(atm);}
  double ElementFormTheo(size_t atomicnr) {
    //function returning theoretical heats of formation for elements; values directly in a.u.
    double enth = 0.0;
    switch (atomicnr) {
      case 1:      //H
        enth = -0.46874201793;
        break;
      case 6:      //C
        enth = -4.05540034771;
        break;
      case 7:      //N
        enth = -5.75406193855;
        break;
      case 8:      //O
        enth = -10.6892450762;
        break;
    }
    return enth;
  }
  double betaA0(size_t atomicnr, int L) {
    //function returning the bonding parameters; values stored in eV, returned however in a.u.
    double betaa0 = 0.0;
    switch (atomicnr) {
      case 1:         //H
        betaa0 = -4.87823460;
        break;
      case 6:         //C
        if (L == 0) {betaa0 = -11.91001500;}
        else if (L == 1) {betaa0 = -9.80275500;}
        break;
      case 7:         //N
        if (L == 0) {betaa0 = -14.33884665;}
        else if (L == 1) {betaa0 = -19.30862853;}
        break;
      case 8:         //O
        if (L == 0) {betaa0 = -46.87741023;}
        else if (L == 1) {betaa0 = -24.74232518;}
        break;
    }
    return betaa0/au2eV;
  }
  double UlX(size_t atomicnr, int L) {
    //function that returns the atomic potential U
    double ulx = 0.0;
    switch (atomicnr) {
      case 1:                 //H
        ulx = -12.75511880;
        break;
      case 6:                 //C
        if (L == 0) {ulx = -47.27032000;}
        else if (L == 1) {ulx = -36.2669180;}
        break;
      case 7:                 //N
        if (L == 0) {ulx = -48.79493385;}
        else if (L == 1) {ulx = -46.57945903;}
        break;
      case 8:                 //O
        if (L == 0) {ulx = -87.38709324;}
        else if (L == 1) {ulx = -71.70268570;}
        break;
    }
    return ulx/au2eV;
  }
  double alpha(int atomicnr, int atm2) {
    //function that returns alpha values for PM3BP
    //alphas stored and returned in 1/Angstrom
    double alpha = 0.0;
    switch (atomicnr) {
      case 1:     //H
        alpha = 3.35638600;
        break;
      case 6:     //C
        alpha = 2.70780700;
        break;
      case 7:     //N
        alpha = 2.83054500;
        break;
      case 8:     //O
        alpha = 3.21710200;
        break;
    }
    return alpha;
  }
  std::vector<double> AM1K(int atomicnr, int atm2 = 0) {
    std::vector<double> am1k;
    switch (atomicnr) {
      case 1:     //H
        am1k.push_back(1.12172594/au2eV);
        am1k.push_back(-1.06492525/au2eV);
        break;
      case 6:     //C
        am1k.push_back(0.05010700/au2eV);
        am1k.push_back(0.05073300/au2eV);
        break;
      case 7:     //N
        am1k.push_back(1.50167153/au2eV);
        am1k.push_back(-1.51571618/au2eV);
        break;
      case 8:     //O
        am1k.push_back(-1.13117677/au2eV);
        am1k.push_back(1.13098909/au2eV);
        break;
    }
    return am1k;
  }
  std::vector<double> AM1L(int atomicnr, int atm2 = 0) {
    std::vector<double> am1l;
    switch (atomicnr) {
      case 1:     //H
        am1l.push_back(5.09516707);
        am1l.push_back(6.02315366);
        break;
      case 6:     //C
        am1l.push_back(6.00316500);
        am1l.push_back(6.00297900);
        break;
      case 7:     //N
        am1l.push_back(5.90399175);
        am1l.push_back(5.97579498);
        break;
      case 8:     //O
        am1l.push_back(6.00999815);
        am1l.push_back(5.87216545);
        break;
    }
    return am1l;
  }
  std::vector<double> AM1M(int atomicnr, int atm2 = 0) {
    std::vector<double> am1m;
    switch (atomicnr) {
      case 1:     //H
        am1m.push_back(1.53693700);
        am1m.push_back(1.57130732);
        break;
      case 6:     //C
        am1m.push_back(1.64221400);
        am1m.push_back(0.89248800);
        break;
      case 7:     //N
        am1m.push_back(1.71042669);
        am1m.push_back(1.71093513);
        break;
      case 8:     //O
        am1m.push_back(1.60731100);
        am1m.push_back(1.60347421);
        break;
    }
    return am1m;
  }
  double Dvalue(size_t atmnr, size_t idx) {
    //function returning the D values needed to calculate eris; values stored in Angstrom, but returned in a.u.
    //note that this function returns both D1 and D2; idx is then either 1 or 2
    double dval = 0.0;
    switch (atmnr) {
      case 1:                  //H
        dval = 0.0;
      case 6:                  //C
        if (idx == 1) {dval = 0.4409314053780;}
        else if (idx == 2) {dval = 0.3517837786516;}
        break;
      case 7:                  //N
        if (idx == 1) {dval = 0.3480401674091;}
        else if (idx == 2) {dval = 0.2801137638434;}
        break;
      case 8:                  //O
        if (idx == 1) {dval = 0.2162309620795;}
        else if (idx == 2) {dval = 0.2712423725348;}
        break;
    }
    return dval*dist_Angstrom2aum1;
  }
  double rho(size_t atmnr, size_t l) {
    //function returning the rho values needed to calculate eris; values stored in Angstrom but returned in atomic units
    double rho = 0.0;
    switch (atmnr) {
      case 1:                  //H
        if (l == 0) {rho = 0.4792425114721178;}
        else if (l == 1) {rho = 0.0;}
        else if (l == 2) {rho = 0.0;}
        break;
      case 6:                  //C
        if (l == 0) {rho = 0.6428006131604468;}
        else if (l == 1) {rho = 0.4495426970001134;}
        else if (l == 2) {rho = 0.34597896187068383;}
        break;
      case 7:                  //N
        if (l == 0) {rho = 0.5799189343842014;}
        else if (l == 1) {rho = 0.5259276440914536;}
        else if (l == 2) {rho = 0.2321791302258449;}
        break;
      case 8:                  //O
        if (l == 0) {rho = 0.4717592829251376;}
        else if (l == 1) {rho = 0.4992677822133302;}
        else if (l == 2) {rho = 0.32420165686630137;}
        break;
    }
    return rho*dist_Angstrom2aum1;
  }
  double eri1Center(int atmnr, int Lbra, int Lket) {
    //function that gives back the semi-empirical 1-center eris
    //values stored in eV but returned in a.u.
    //Lbra is the sum of azimuthal quantum numbers for bra (ss = 0; pp = 2; sp = 1)
    //Lket is the sum of azimuthal quantum numbers for ket (ss = 0; pp = 2; sp = 1; p*p* = -2)
    double eri = 0.0;
    switch (atmnr) {
      case 1:       //H
        eri = 15.02333745;                                                                        //(ss|ss)
        break;
      case 6:       //C
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.20070800;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 10.79629200;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 10.26502700;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.04256600;}    //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.29098000;}                                    //(sp|sp)||(ps|ps)
        break;
      case 7:       //N
        if ((Lbra == 0)&&(Lket == 0)) {eri = 12.41522141;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 13.96611412;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.34540226;}      //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 10.41021925;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.13671300;}                                    //(sp|sp)||(ps|ps)
        break;
      case 8:       //O
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.26164345;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 13.65930075;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 10.41332625;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 12.42051017;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.59388300;}                                    //(sp|sp)||(ps|ps)
        break;
    }
    return eri/au2eV;
  }
};
class RM1: public MNDO {
  //this is the implementation of RM1
  //G. B. Rocha, R. O. Freire, A. M. Simas, J. J. P. Stewart, J. Comp. Chem., 27(10), 1101, 2005
  //bX quantities are "barred" tensors, which are used only for the open-shell case
public:
  RM1(BSet _bset, Molecule _mol, std::string _openclosed = "0", std::string _corecorrection = "0"): MNDO(_bset,_mol,_openclosed,_corecorrection) {}
  ~RM1() {}
  std::string Type() {return "RM1";}
  double gfactor(size_t iatm1, size_t iatm2, double RAB) {return 1.0;}
  virtual double gfactor_dR(size_t iatm1, size_t iatm2) {return 0.0;}
  double AM1factor(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {
    //RAB in a.u.
    double AM1factor = 0.0;
    RAB *= dist_Angstrom2au;
    for (size_t idx = 0; idx < K1.size(); ++idx) {
      AM1factor += K1[idx]*exp(-L1[idx]*(RAB - M1[idx])*(RAB - M1[idx]))/RAB;
    }
    return AM1factor;
  }
  virtual double AM1factor_dR(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {
    //RAB in a.u.
    double am1factor_dr = 0.0;
    double RR = RAB*dist_Angstrom2au;
    for (size_t idx = 0; idx < K1.size(); ++idx) {
      am1factor_dr += K1[idx]*(2.0*RR*RR*L1[idx] - 2.0*L1[idx]*M1[idx]*RR + 1.0)*exp(-L1[idx]*(RR - M1[idx])*(RR - M1[idx]))/(RAB*RR);
    }
    return am1factor_dr;  
  }
  virtual double AM1factor_dR2(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {
    //RAB in a.u.
    double am1factor_dr2 = 0.0;
    double RR = RAB*dist_Angstrom2au;
    for (size_t idx = 0; idx < K1.size(); ++idx) {
      am1factor_dr2 += 2.0*K1[idx]*exp(-L1[idx]*(RR - M1[idx])*(RR - M1[idx]))*(L1[idx]*M1[idx]*RR - 1.0 - RR*L1[idx]*(RR - M1[idx])*(2.0*RR*RR*L1[idx] - 2.0*L1[idx]*M1[idx]*RR + 1.0))/(RAB*RAB*RR);
    }
    return am1factor_dr2;  
  }
  //other auxiliary functions
  void checkAtoms() {
    //function that checks whether the atoms in the molecule are fully defined in the respective MNDO theories
    bool def;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      def = false;
      if (atoms[iatm] == 1) {def = true;}                                   //H
      else if ((atoms[iatm] > 5)&&(atoms[iatm] < 10)) {def = true;}         //C,N,O,F
      else if ((atoms[iatm] > 14)&&(atoms[iatm] < 18)) {def = true;}        //P,S,Cl
      else if ((atoms[iatm] == 35)||(atoms[iatm] == 53)) {def = true;}      //Br,I
      if (!def) {throw("ERROR: MNDO.hpp: RM1: checkAtoms(): atom not fully specified for RM1-theory");}
    }
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapRM1(atm);}
  double ElementFormTheo(size_t atomicnr) {
    //function returning theoretical heats of formation for elements; values directly in a.u.
    double enth = 0.0;
    switch (atomicnr) {
      case 1:      //H
        enth = -0.439546816194;
        break;
      case 6:      //C
        enth = -4.30899631856;
        break;
      case 7:      //N
        enth = -7.50452619284;
        break;
      case 8:      //O
        enth = -11.440778309;
        break;
      case 9:      //F
        enth = -17.8085654062;
        break;
      case 15:     //P
        enth = -4.52307090781;
        break;
      case 16:     //S
        enth = -6.79711669201;
        break;
      case 17:     //Cl
        enth = -14.0555181302;
        break;
      case 35:     //Br
        enth = -13.1237879735;
        break;
      case 53:     //I
        enth = -9.13196215556;
        break;
    }
    return enth;
  }
  double betaA0(size_t atomicnr, int L) {
    //function returning the bonding parameters; values stored in eV, returned however in a.u.
    double betaa0 = 0.0;
    switch (atomicnr) {
      case 1:                         //H
        betaa0 = -5.76544469;
        break;
      case 6:                         //C
        if (L == 0) {betaa0 = -15.45932428;}
        else if (L == 1) {betaa0 = -8.23608638;}
        break;
      case 7:                         //N
        if (L == 0) {betaa0 = -20.87124548;}
        else if (L == 1) {betaa0 = -16.67171853;}
        break;
      case 8:                         //O
        if (L == 0) {betaa0 = -29.85101212;}
        else if (L == 1) {betaa0 = -29.15101314;}
        break;
      case 9:                         //F
        if (L == 0) {betaa0 = -70.00000512;}
        else if (L == 1) {betaa0 = -32.67982711;}
        break;
      case 15:                        //P
        if (L == 0) {betaa0 = -6.13514969;}
        else if (L == 1) {betaa0 = -5.94442127;}
        break;
      case 16:                        //S
        if (L == 0) {betaa0 = -1.95910719;}
        else if (L == 1) {betaa0 = -8.77430652;}
        break;
      case 17:                        //Cl
        if (L == 0) {betaa0 = -19.92430432;}
        else if (L == 1) {betaa0 = -11.52935197;}
        break;
      case 35:                        //Br
        if (L == 0) {betaa0 = -1.34139841;}
        else if (L == 1) {betaa0 = -8.20225991;}
        break;
      case 53:                        //I
        if (L == 0) {betaa0 = -4.19316149;}
        else if (L == 1) {betaa0 = -4.40038412;}
        break;
    }
    return betaa0/au2eV;
  }
  double UlX(size_t atomicnr, int L) {
    //function that returns the atomic potential U
    double ulx = 0.0;
    switch (atomicnr) {
      case 1:                 //H
        ulx = -11.96067697;
        break;
      case 6:                 //C
        if (L == 0) {ulx = -51.72556032;}
        else if (L == 1) {ulx = -39.40728943;}
        break;
      case 7:                 //N
        if (L == 0) {ulx = -70.85123715;}
        else if (L == 1) {ulx = -57.97730920;}
        break;
      case 8:                 //O
        if (L == 0) {ulx = -96.94948069;}
        else if (L == 1) {ulx = -77.89092978;}
        break;
      case 9:                 //F
        if (L == 0) {ulx = -134.18369591;}
        else if (L == 1) {ulx = -107.84660920;}
        break;
      case 15:                //P
        if (L == 0) {ulx = -41.81533184;}
        else if (L == 1) {ulx = -34.38342529;}
        break;
      case 16:                //S
        if (L == 0) {ulx = -55.16775121;}
        else if (L == 1) {ulx = -46.52930422;}
        break;
      case 17:                //Cl
        if (L == 0) {ulx = -118.47306918;}
        else if (L == 1) {ulx = -76.35330340;}
        break;
      case 35:                //Br
        if (L == 0) {ulx = -113.48398183;}
        else if (L == 1) {ulx = -76.18720023;}
        break;
      case 53:                //I
        if (L == 0) {ulx = -74.89997837;}
        else if (L == 1) {ulx = -51.41023805;}
        break;
    }
    return ulx/au2eV;
  }
  double alpha(int atomicnr, int atm2) {
    //function that returns alpha values for AM1
    //alphas stored and returned in 1/Angstrom
    double alpha = 0.0;
    switch (atomicnr) {
      case 1:     //H
        alpha = 3.06835947;
        break;
      case 6:     //C
        alpha = 2.79282078;
        break;
      case 7:     //N
        alpha = 2.96422542;
        break;
      case 8:     //O
        alpha = 4.17196717;
        break;
      case 9:     //F
        alpha = 6.00000062;
        break;
      case 15:    //P
        alpha = 1.90993294;
        break;
      case 16:    //S
        alpha = 2.44015636;
        break;
      case 17:    //Cl
        alpha = 3.69358828;
        break;
      case 35:    //Br
        alpha = 2.86710532;
        break;
      case 53:    //I
        alpha = 2.14157092;
        break;
    }
    return alpha;
  }
  std::vector<double> AM1K(int atomicnr, int atm2 = 0) {
    std::vector<double> am1k;
    switch (atomicnr) {
      case 1:     //H
        am1k.push_back(0.10288875/au2eV);
        am1k.push_back(0.06457449/au2eV);
        am1k.push_back(-0.03567387/au2eV);
        break;
      case 6:     //C
        am1k.push_back(0.07462271/au2eV);
        am1k.push_back(0.01177053/au2eV);
        am1k.push_back(0.03720662/au2eV);
        am1k.push_back(-0.00270657/au2eV);
        break;
      case 7:     //N
        am1k.push_back(0.06073380/au2eV);
        am1k.push_back(0.02438558/au2eV);
        am1k.push_back(-0.02283430/au2eV);
        break;
      case 8:     //O
        am1k.push_back(0.23093552/au2eV);
        am1k.push_back(0.05859873/au2eV);
        break;
      case 9:     //F
        am1k.push_back(0.40302025/au2eV);
        am1k.push_back(0.07085831/au2eV);
        break;
      case 15:    //P
        am1k.push_back(-0.41063467/au2eV);
        am1k.push_back(-0.16299288/au2eV);
        am1k.push_back(-0.04887125/au2eV);
        break;
      case 16:    //S
        am1k.push_back(-0.74601055/au2eV);
        am1k.push_back(-0.06519286/au2eV);
        am1k.push_back(-0.00655977/au2eV);
        break;
      case 17:    //Cl
        am1k.push_back(0.12947108/au2eV);
        am1k.push_back(0.00288899/au2eV);
        break;
      case 35:    //Br
        am1k.push_back(0.98689937/au2eV);
        am1k.push_back(-0.92731247/au2eV);
        break;
      case 53:    //I
        am1k.push_back(-0.08147724/au2eV);
        am1k.push_back(0.05914991/au2eV);
        break;
    }
    return am1k;
  }
  std::vector<double> AM1L(int atomicnr, int atm2 = 0) {
    std::vector<double> am1l;
    switch (atomicnr) {
      case 1:     //H
        am1l.push_back(5.90172268);
        am1l.push_back(6.41785671);
        am1l.push_back(2.80473127);
        break;
      case 6:     //C
        am1l.push_back(5.73921605);
        am1l.push_back(6.92401726);
        am1l.push_back(6.26158944);
        am1l.push_back(9.00003735);
        break;
      case 7:     //N
        am1l.push_back(4.58892946);
        am1l.push_back(4.62730519);
        am1l.push_back(2.05274659);
        break;
      case 8:     //O
        am1l.push_back(5.21828736);
        am1l.push_back(7.42932932);
        break;
      case 9:     //F
        am1l.push_back(7.20441959);
        am1l.push_back(9.00001562);
        break;
      case 15:    //P
        am1l.push_back(6.08752832);
        am1l.push_back(7.09472602);
        am1l.push_back(8.99979308);
        break;
      case 16:    //S
        am1l.push_back(4.81038002);
        am1l.push_back(7.20760864);
        am1l.push_back(9.00000180);
        break;
      case 17:    //Cl
        am1l.push_back(2.97724424);
        am1l.push_back(7.09827589);
        break;
      case 35:    //Br
        am1l.push_back(4.28484191);
        am1l.push_back(4.54005910);
        break;
      case 53:    //I
        am1l.push_back(1.56065072);
        am1l.push_back(5.76111270);
        break;
    }
    return am1l;
  }
  std::vector<double> AM1M(int atomicnr, int atm2 = 0) {
    std::vector<double> am1m;
    switch (atomicnr) {
      case 1:     //H
        am1m.push_back(1.17501185);
        am1m.push_back(1.93844484);
        am1m.push_back(1.63655241);
        break;
      case 6:     //C
        am1m.push_back(1.04396983);
        am1m.push_back(1.66159571);
        am1m.push_back(1.63158721);
        am1m.push_back(2.79557901);
        break;
      case 7:     //N
        am1m.push_back(1.37873881);
        am1m.push_back(2.08370698);
        am1m.push_back(1.86763816);
        break;
      case 8:     //O
        am1m.push_back(0.90363555);
        am1m.push_back(1.51754610);
        break;
      case 9:     //F
        am1m.push_back(0.81653013);
        am1m.push_back(1.43802381);
        break;
      case 15:    //P
        am1m.push_back(1.31650261);
        am1m.push_back(1.90721319);
        am1m.push_back(2.65857780);
        break;
      case 16:    //S
        am1m.push_back(0.59380129);
        am1m.push_back(1.29492008);
        am1m.push_back(1.80060151);
        break;
      case 17:    //Cl
        am1m.push_back(1.46749784);
        am1m.push_back(2.50002723);
        break;
      case 35:    //Br
        am1m.push_back(2.00019696);
        am1m.push_back(2.01617695);
        break;
      case 53:    //I
        am1m.push_back(2.00002063);
        am1m.push_back(2.20488800);
        break;
    }
    return am1m;
  }
  double Dvalue(size_t atmnr, size_t idx) {
    //function returning the D values needed to calculate eris; values stored in Angstrom, but returned in a.u.
    //note that this function returns both D1 and D2; idx is then either 1 or 2
    double dval = 0.0;
    switch (atmnr) {
      case 1:                   //H
        dval = 0.0;
        break;
      case 6:                   //C
        if (idx == 1) {dval = 0.421625705;}
        else if (idx == 2) {dval = 0.366514014;}
        break;
      case 7:                  //N
        if (idx == 1) {dval = 0.343733389;}
        else if (idx == 2) {dval = 0.327636952;}
        break;
      case 8:                  //O
        if (idx == 1) {dval = 0.258593059;}
        else if (idx == 2) {dval = 0.253799433;}
        break;
      case 9:                  //F
        if (idx == 1) {dval = 0.184626083;}
        else if (idx == 2) {dval = 0.244715022;}
        break;
      case 15:                 //P
        if (idx == 1) {dval = 0.534837011;}
        else if (idx == 2) {dval = 0.507940817;}
        break;
      case 16:                 //S
        if (idx == 1) {dval = 0.525839200;}
        else if (idx == 2) {dval = 0.472356644;}
        break;
      case 17:                 //Cl
        if (idx == 1) {dval = 0.240341095;}
        else if (idx == 2) {dval = 0.467043700;}
        break;
      case 35:                 //Br
        if (idx == 1) {dval = 0.111074578;}
        else if (idx == 2) {dval = 0.552580723;}
        break;
      case 53:                 //I
        if (idx == 1) {dval = 0.685994904;}
        else if (idx == 2) {dval = 0.586643919;}
        break;
    }
    return dval*dist_Angstrom2aum1;
  }
  double rho(size_t atmnr, size_t l) {
    //function returning the rho values needed to calculate eris; values stored in Angstrom but returned in atomic units
    double rho = 0.0;
    switch (atmnr) {
      case 1:                  //H
        if (l == 0) {rho = 0.51488046;}
        else if (l == 1) {rho = 0.0;}
        else if (l == 2) {rho = 0.0;}
        break;
      case 6:                  //C
        if (l == 0) {rho = 0.55156780;}
        else if (l == 1) {rho = 0.51905507;}
        else if (l == 2) {rho = 0.39998984;}
        break;
      case 7:                  //N
        if (l == 0) {rho = 0.55012484;}
        else if (l == 1) {rho = 0.27259716;}
        else if (l == 2) {rho = 0.32975880;}
        break;
      case 8:                  //O
        if (l == 0) {rho = 0.51417391;}
        else if (l == 1) {rho = 0.26285897;}
        else if (l == 2) {rho = 0.29551626;}
        break;
      case 9:                  //F
        if (l == 0) {rho = 0.43057954;}
        else if (l == 1) {rho = 0.28680176;}
        else if (l == 2) {rho = 0.42425093;}
        break;
      case 15:                 //P
        if (l == 0) {rho = 0.64975614;}
        else if (l == 1) {rho = 0.67287406;}
        else if (l == 2) {rho = 0.84709086;}
        break;
      case 16:                 //S
        if (l == 0) {rho = 0.57651500;}
        else if (l == 1) {rho = 0.38175294;}
        else if (l == 2) {rho = 0.53350588;}
        break;
      case 17:                 //Cl
        if (l == 0) {rho = 0.46872232;}
        else if (l == 1) {rho = 0.35326854;}
        else if (l == 2) {rho = 0.34329568;}
        break;
      case 35:                 //Br
        if (l == 0) {rho = 0.42064959;}
        else if (l == 1) {rho = 0.20096728;}
        else if (l == 2) {rho = 0.45124871;}
        break;
      case 53:                 //I
        if (l == 0) {rho = 0.35998881;}
        else if (l == 1) {rho = 0.71175720;}
        else if (l == 2) {rho = 0.75323162;}
        break;
    }
    return rho*dist_Angstrom2aum1;
  }
  double eri1Center(int atmnr, int Lbra, int Lket) {
    //function that gives back the semi-empirical 1-center eris
    //values stored in eV but returned in a.u.
    //Lbra is the sum of azimuthal quantum numbers for bra (ss = 0; pp = 2; sp = 1)
    //Lket is the sum of azimuthal quantum numbers for ket (ss = 0; pp = 2; sp = 1; p*p* = -2)
    double eri = 0.0;
    switch (atmnr) {
      case 1:       //H
        eri = 13.98321296;                                                                       //(ss|ss)
        break;
      case 6:       //C
        if ((Lbra == 0)&&(Lket == 0)) {eri = 13.05312440;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 10.95113739;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 11.33479389;}    //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.72395099;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.55215133;}                                   //(sp|sp)||(ps|ps)
        break;
      case 7:       //N
        if ((Lbra == 0)&&(Lket == 0)) {eri = 13.08736234;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 13.69924324;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 13.21226834;}    //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 11.94103953;}  //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 5.00000846;}                                   //(sp|sp)||(ps|ps)
        break;
      case 8:       //O
        if ((Lbra == 0)&&(Lket == 0)) {eri = 14.00242788;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 14.14515138;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 14.95625043;}    //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 12.70325497;}  //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.93217161;}                                   //(sp|sp)||(ps|ps)
        break;
      case 9:       //F
        if ((Lbra == 0)&&(Lket == 0)) {eri = 16.72091319;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 15.22581028;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 16.76142629;}    //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 14.86578679;}  //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.99766171;}                                   //(sp|sp)||(ps|ps)
        break;
      case 15:      //P
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.08059265;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.60417563;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.68339201;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.40265182;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.16181792;}                                   //(sp|sp)||(ps|ps)
        break;
      case 16:      //S
        if ((Lbra == 0)&&(Lket == 0)) {eri = 12.48828408;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.52301167;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.56910574;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.66863296;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.88978932;}                                   //(sp|sp)||(ps|ps)
        break;
      case 17:      //Cl
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.36023105;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 12.56502640;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 13.30671171;}    //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.66397083;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.76489897;}                                   //(sp|sp)||(ps|ps)
        break;
      case 35:      //Br
        if ((Lbra == 0)&&(Lket == 0)) {eri = 17.11563074;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 10.73546293;}                                  //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 15.62419251;}    //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 8.86056199;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.23512762;}                                   //(sp|sp)||(ps|ps)
        break;
      case 53:      //I
        if ((Lbra == 0)&&(Lket == 0)) {eri = 19.99974131;}                                       //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.30488343;}                                   //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.68957672;}     //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.85424614;}   //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.41602940;}                                   //(sp|sp)||(ps|ps)
        break;
    }
    return eri/au2eV;
  }
};

#endif //_MNDO_
