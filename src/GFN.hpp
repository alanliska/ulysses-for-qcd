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

#ifndef _GFN_METHODS_
#define _GFN_METHODS_
#include "QC.hpp"
#include "Hessian.hpp"
#include "Intermolecular.hpp"
#include "other/listPackage.hpp"
#include "math/SolverPackage.hpp"
#include "parameters/GFN2par.hpp"
#include "parameters/ZeroOvpar.hpp"
#include "solvation/SolvationGFN2par.hpp"

//descrition:
//xTB methods of Grimme

class GFN2: public QCbasis {
  //this is the implementation of Grimme's GFN2-xTB
  //C. Bannwarth, S. Ehrlert, S. Grimme, J. Chem. Theory Comput., 15, 1652, 2019
protected:
  size_t dispsz;
  int nopen;                              //total number of open-shell electrons
  int totalshells;
  int coreelectrons;
  int Nelectrons;
  std::vector<int> AtomShells;
  std::vector<int> refA;
  matrixE QShells;
  matrixE muAtoms;
  matrixE thetaAtoms;
  matrixE VD;                             //dipole integral proportional potential from all atoms acting on atom a certain atom
  matrixE VQ;                             //quadrupole integral proportional potential from all atoms acting on atom a certain atom
  matrixE bmatrix;                        //for Broyden damping
  matrixE dfmat;                          //for Broyden damping
  matrixE amat;                           //for Broyden damping
  matrixE umat;                           //for Broyden damping
  matrixE WDisp;                          //matrix with the weighted dispersion coefficients
  matrixE VES;                            //shell-shift required for gradients
  std::vector<double> alphavec;
  std::vector<double> betavec;
  std::vector<double> CoordNum;
  std::vector<double> CNrad;
  std::vector<double> dAtoms;
  std::vector<double> qvec;               //for Broyden damping
  std::vector<double> qvec_last;          //for Broyden damping
  std::vector<double> dqvec;              //for Broyden damping
  std::vector<double> dqvec_last;         //for Broyden damping
  std::vector<double> omega;              //for Broyden damping
  std::vector<double> gvec;               //for Broyden damping
  std::vector<double> GaussW;             //Gaussian weights for the D4 correction
  std::vector<double> zetaA;              //zeta vector
  std::vector<double> dzetaA;             //derivative of the zeta vector with respect to charge
  std::vector<double> auxv;               //auxiliary vector
  std::vector<double> VS;                 //overlap proportional potential acting on a given atom
  std::vector<double> occa;               //this is occupation number here
  std::vector<double> occb;               //this is occupation number here
  std::vector<double> rA;                 //auxiliary
  std::vector<double> rB;                 //auxiliary
  std::vector<double> muA;                //auxiliary
  std::vector<double> muB;                //auxiliary
  std::vector<double> muTemp;             //auxiliary
  std::vector<double> thetaTemp;          //auxiliary
  std::vector<double> VQA;                //auxiliary
  std::vector<double> VQB;                //auxiliary
  double Broydendamp;
  double omega0;
  double minw;
  double maxw;
  double wfac;
  double energyHueckel;                   //Hueckel-type energy
  double energyIES_IXC;                   //isotropic electrostatic and isotropic exchange contributions to energy
  double eaes;
  double eaxc;
  double ed4;
  double eabc;
  double efermia;
  double efermib;
  double gwthresh;
  bool extractC6;
public:
  GFN2(BSet _bset, Molecule _mol, const std::string & _openclosed = "closed"): QCbasis(_bset,_mol) {
    openclosed = _openclosed;
    Broydendamp = 0.4;
    omega0 = 0.01;
    minw = 1.0;
    maxw = 100000.0;
    wfac = 0.01;
    gwthresh = 1.0e-7;
    extractC6 = true;
    Telectron = 300.0;
    corecorrection = 0.0;
    rA.resize(3);
    rB.resize(3);
    rAB.resize(3);
    muA.resize(3);
    muB.resize(3);
    muTemp.resize(3);
    thetaTemp.resize(6);
    VQA.resize(6);
    VQB.resize(6);
    thresh_dens = 1.0e-4;
    thresh_en = 1.0e-6;
    nopen = 0;
    thresholdzero = 1.0e-20;
  }
  ~GFN2() {}
  //setters
  void setQShells(matrixE & _QShells) {QShells = _QShells;}
  void setmuAtoms(matrixE & _muAtoms) {muAtoms = _muAtoms;}
  void setthetaAtoms(matrixE & _thetaAtoms) {thetaAtoms = _thetaAtoms;}
  void setqvec(std::vector<double> & _qvec) {qvec = _qvec;}
  void setqvec_last(std::vector<double> & _qvec_last) {qvec_last = _qvec_last;}
  void setdqvec(std::vector<double> & _dqvec) {dqvec = _dqvec;}
  void setdqvec_last(std::vector<double> & _dqvec_last) {dqvec_last = _dqvec_last;}
  void setomega(std::vector<double> & _omega) {omega = _omega;}
  void setgvec(std::vector<double> & _gvec) {gvec = _gvec;}
  void setVS(std::vector<double> & _VS) {VS = _VS;}
  void setVD(matrixE & _VD) {VD = _VD;}
  void setVQ(matrixE & _VQ) {VQ = _VQ;}
  void setbmatrix(matrixE & _bmatrix) {bmatrix = _bmatrix;}
  void setdfmat(matrixE & _dfmat) {dfmat = _dfmat;}
  void setamat(matrixE & _amat) {amat = _amat;}
  void setumat(matrixE & _umat) {umat = _umat;}
  void setExtractC6(bool _extractC6) {extractC6 = _extractC6;}
  void setNOpen(size_t _nopen) {nopen = _nopen;}
  void setVES(matrixE & _VES) {VES = _VES;}
  //getters
  std::string Type() {return "GFN2";}
  virtual std::string TypeOverlap() {return "Non-Orthogonal";}
  matrixE getVAB() {return VAB;}
  std::vector<double> getVS() {return VS;}
  matrixE getVD() {return VD;}
  matrixE getVQ() {return VQ;}
  matrixE getQShells() {return QShells;}
  matrixE getmuAtoms() {return muAtoms;}
  matrixE getthetaAtoms() {return thetaAtoms;}
  std::vector<double> getdAtoms() {return dAtoms;}
  std::vector<double> getCoordNum() {return CoordNum;}
  bool getExtractC6() {return extractC6;}
  double getEIESIXC() {return energyIES_IXC;}
  double getEAES() {return eaes;}
  double getEAXC() {return eaxc;}
  double getED4() {return ed4;}
  double getEABC() {return eabc;}
  double EFermiA() {return efermia;}
  double EFermiB() {return efermib;}
  size_t getNOpen() {return nopen;}
  matrixE getVES() {return VES;}
  //Calculation
  void Calculate(int _print = 1, size_t maxiter = 200, bool _DIIS = true, bool _RCA = true, double d2threshold = 1600.0) {
    //check whether calculation is possible
    checkAtoms();
    matrixE geometry = mol.Geometry();
    matrixE xao = basis.XAO();
    if (usesolvation) {
      if (CoulombSolv == 1) {calcBornMatrix(geometry,GBEffectiveCoulombP16);}
      else if (CoulombSolv == 0) {calcBornMatrix(geometry,GBEffectiveCoulombStill);}
    }
    refA.resize(Natoms);
    //initialize and resize necessary quantities
    CalcAtomShells();
    Nelectrons = mol.Nelectrons();
    print = _print;
    NAOs = basis.NAO();
    sao = basis.SAO();
    //determine whether system is closed- or open-shell
    OpenClosed();
    if (print > 0) {std::cout << shell << "-shell calculation" << std::endl;}
    if (occupancy.size() != NAOs) {
      occupancy.resize(NAOs);
      occa.resize(NAOs);
      occb.resize(NAOs);
    }
    auxv.resize(7);
    if (QShells.rows() != totalshells) {QShells.resize(totalshells,1);}
    if (VES.rows() != totalshells) {VES.resize(totalshells,1);}
    if (QAtoms.size() != Natoms) {
      muAtoms.resize(Natoms,3);
      thetaAtoms.resize(Natoms,6);
      QAtoms.resize(Natoms);
      dAtoms.resize(Natoms);
      VS.resize(Natoms);
      VD.resize(Natoms,3);
      VQ.resize(Natoms,6);
    }
    if (qvec.size() != totalshells + 9*Natoms) {
      qvec.resize(totalshells + 9*Natoms);
      qvec_last.resize(totalshells + 9*Natoms);
      dqvec.resize(totalshells + 9*Natoms);
      dqvec_last.resize(totalshells + 9*Natoms);
    }
    if (restart < 1) {
      for (size_t idatm = 0; idatm < Natoms; ++idatm) {
        for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
          muAtoms(idatm + 1,idcoord) = 0.0;
          thetaAtoms(idatm + 1,idcoord) = 0.0;
          thetaAtoms(idatm + 1,idcoord + 3) = 0.0;
        }
      }
      for (size_t idx = 0; idx < totalshells + 9*Natoms; ++idx) {
        qvec[idx] = 0.0;
        qvec_last[idx] = 0.0;
        dqvec[idx] = 0.0;
        dqvec_last[idx] = 0.0;
      }
      omega.resize(25);
      dfmat.resize(25,totalshells + 9*Natoms);
      amat.resize(25,25);
      umat.resize(25,totalshells + 9*Natoms);
      dens.resize(NAOs,NAOs);
      //get the base charges, dipole moment and quadrupole moment
      GoedeckerCharges(QAtoms,atoms,geometry,mol.Charge());
      initQShell(qvec);
      if (restart == 0) {restart = 1;}
    }
    double energy = 0.0;
    double oldenergy = 0.0;
    double dqcontrol;
    efermia = 0.0;
    efermib = 0.0;
    //get dispersion things
    getWeightedGFN2(WDisp,auxMAT2,GaussW,atoms,geometry,refA,auxv,extractC6,thresholdzero,d2threshold,gwthresh);
    dispsz = GaussW.size();
    zetaA.resize(dispsz);
    dzetaA.resize(dispsz);
    //integrals and other fixed stuff
    NCoord(CoordNum,atoms,geometry,GFNcount,false,d2threshold,false);
    CNradius();
    CalcGammaAB(geometry);
    calcHcore(sao,geometry);
    //get occupation vector
    Occ();
    setOccupancy();
    if (print > 0) {
      std::cout << "occupation vector " << std::endl;
      for (size_t id = 0; id < occupancy.size(); ++id) {
        std::cout << occupancy[id] << " ";
      }
      std::cout << std::endl;
    }
    //SCF
    for (size_t iter = 1; iter < maxiter; ++iter) {
      if (print > 0) {std::cout << "iteration " << iter << std::endl;}
      if (iter%25 == 0) {
        omega.resize(iter + 25);
        dfmat.resize(iter + 25,totalshells + 9*Natoms);
        amat.resize(iter + 25,iter + 25);
        umat.resize(iter + 25,totalshells + 9*Natoms);
      }
      //get the electrostatic potential
      calcVSDQ(geometry);
      //get new MOs
      calcFock(sao);            //only VS, VD and VQ are calculated in a function that requires the geometry; rest obtained in calcFock
      transformFock(xao);
      calcMOs();
      transformCMO(xao);
      //Fermi smearing
      GFermi();
      //get new density
      Dens(dens,occupancy,CMO,DensAUX,true,1.0e-8);
      if (print == 2) {
        std::cout << "density" << std::endl;
        dens.Print(17);
      }
      //get energies
      EElectrostaticDispersion(dens,geometry);
      EHueckel(dens);
      if (usesolvation) {esolv = SolvationEnergy(QAtoms);}
      energy = eaxc + eaes + energyIES_IXC + ed4 + energyHueckel + efermia + efermib + esolv;
      if (print > 0) {
        if (print == 2) {
          std::cout << "CMO " << std::endl;
          CMO.Print(17);
        }
        std::cout << "Energy = " << energy << std::endl;
        std::cout << "orbital energies: ";
        for (size_t id = 0; id < EMOs.size(); ++id) {
          std::cout << EMOs[id] << " ";
        }
        std::cout << std::endl;
      }
      //control charges for convergence
      dqcontrol = 0.0;
      for (size_t idshell = 0; idshell < totalshells; ++idshell) {
        dqcontrol += dqvec[idshell]*dqvec[idshell];
      }
      if ((fabs(oldenergy - energy) < thresh_en)&&(dqcontrol < thresh_dens*thresh_dens)) {
        if (print > 0) {std::cout << "convergence after " << iter << " iterations" << std::endl;}
        break;
      }
      oldenergy = energy;
      //accelerate convergence
      Broyden(iter);
    }
    ErepATM(geometry,d2threshold);
    etot = energy + enuc + eabc;
    if (print > 0) {
      std::cout << "Optimized MOs" << std::endl;
      CMO.Print(17);
      std::cout << "Optimized orbital energies: ";
      for (size_t id = 0; id < EMOs.size(); ++id) {
        std::cout << EMOs[id] << " ";
      }
      std::cout << std::endl;
      std::cout << "Optimized Density" << std::endl;
      dens.Print(17);
      std::cout << "Total Energy = " << etot << std::endl;
      std::cout << "Electronic Energy = " << etot - enuc << std::endl;
      std::cout << "Nuclear Repulsion = " << enuc << std::endl;
      std::cout << "Self-Consistent Energy = " << energy - enuc << std::endl;
      if (usesolvation) {std::cout << "Gsolv = " << esolv << std::endl;}
    }
  }
  //calculate MO occupancies
  void Occ() {
    //function to calculate the occupancy vectors for a system
    int Nelec = Nelectrons - nopen - coreelectrons;    //this has to be an even number
    Nelec /= 2;                                        //total number of pairs of electrons
    for (size_t idx = 0; idx < Nelec; ++idx) {
      occa[idx] = 1;
      occb[idx] = 1;
    }
    for (size_t idx = 0; idx < nopen; ++idx) {
      occa[idx + Nelec] = 1;
    }
  }
  void setOccupancy() {
    for (size_t idorb = 0; idorb < NAOs; ++idorb) {
      occupancy[idorb] = occa[idorb] + occb[idorb];
      occa[idorb] = 0.0;
      occb[idorb] = 0.0;
    }
  }
  void ErepATM(matrixE & geometry, double d2threshold = 1600.0) {
    //function calculating the repulsion energy, equation 9 of the paper and the Axilrod-Teller-Muto term 
    enuc = 0.0;          //zero the energies
    eabc = 0.0;
    double alphaA;
    double alphaB;
    double yeffectiveA;
    double yeffectiveB;
    double krep;
    double sqrtZr4r2A;
    double sqrtZr4r2B;
    double sqrtZr4r2C;
    double c6AB;
    double c6AC;
    double c6BC;
    double cAB;
    double cAC;
    double cBC;
    double a1 = A1D4(this->Type());
    double a2 = A2D4(this->Type());
    double RAB2;
    double RAC2;
    double RBC2;
    double alpha = AlphaD4(this->Type(),9);
    int atomA;
    int atomB;
    int atomC = 0;
    int refa;
    int refb;
    int posA = 0;
    int posB = 0;
    this->calcZetas(0.0);
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      for (size_t idrefA = 0; idrefA < refA[idatm]; ++idrefA,++atomC) {
        zetaA[atomC] *= GaussW[atomC];
      }
    }
    NCoord(auxv,atoms,geometry,ERFcount,true,d2threshold,false);
    getC6D4(this->Type(),WDisp,GaussW,auxv,refA,atoms,extractC6,thresholdzero);
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      atomA = atoms[idAtm];
      refa = RefN(atomA);
      alphaA = sqrt(AlphaGFN2(atomA));
      yeffectiveA = Yeff(atomA);
      sqrtZr4r2A = R2R4(atomA);
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        rA[idcoord] = geometry(idAtm + 1,idcoord + 1);
      }
      posB = 0;
      for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
        atomB = atoms[idBtm];
        refb = RefN(atomB);
        c6AB = 0.0;
        for (size_t idrefA = 0; idrefA < refa; ++idrefA) {
          for (size_t idrefB = 0; idrefB < refb; ++idrefB) {
            c6AB += zetaA[posA + idrefA]*WDisp(posA + idrefA + 1,posB + idrefB + 1)*zetaA[posB + idrefB];
          }
        }
        WDisp(idAtm + 1,idBtm + 1) = c6AB;
        WDisp(idBtm + 1,idAtm + 1) = c6AB;
        posB += refb;
        alphaB = sqrt(AlphaGFN2(atomB));
        yeffectiveB = Yeff(atomB);
        sqrtZr4r2B = R2R4(atomB);
        RAB2 = 0.0;
        krep = 1.0;
        if ((atomA > 2)||(atomB > 2)) {krep += 0.5;}
        krep *= 0.5;
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          rB[idcoord] = geometry(idBtm + 1,idcoord + 1);
          RAB2 += (rB[idcoord] - rA[idcoord])*(rB[idcoord] - rA[idcoord]);
        }
        RAB2 *= dist_Angstrom2aum1*dist_Angstrom2aum1;                            //convert to Bohr
        if (RAB2 > d2threshold) {continue;}
        c6AB = pow(RAB2,krep);
        cAB = yeffectiveA*yeffectiveB*exp(-alphaA*alphaB*c6AB)/sqrt(RAB2);        //equation 9
        enuc += cAB;
        c6AB = WDisp(idAtm + 1,idBtm + 1);
        cAB = a1*sqrt(3.0*sqrtZr4r2A*sqrtZr4r2B) + a2;
        for (size_t idCtm = 0; idCtm < idBtm; ++idCtm) {
          atomC = atoms[idCtm];
          sqrtZr4r2C = R2R4(atomC);
          RAC2 = 0.0;
          RBC2 = 0.0;
          for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
            RAC2 += (geometry(idCtm + 1,idcoord + 1) - rA[idcoord])*(geometry(idCtm + 1,idcoord + 1) - rA[idcoord]);
            RBC2 += (geometry(idCtm + 1,idcoord + 1) - rB[idcoord])*(geometry(idCtm + 1,idcoord + 1) - rB[idcoord]);
          }
          RAC2 *= dist_Angstrom2aum1*dist_Angstrom2aum1;         //convert to Bohr
          RBC2 *= dist_Angstrom2aum1*dist_Angstrom2aum1;         //convert to Bohr
          if (RAC2 > d2threshold) {continue;}
          if (RBC2 > d2threshold) {continue;}
          c6AC = WDisp(idAtm + 1,idCtm + 1);
          c6BC = WDisp(idBtm + 1,idCtm + 1);
          cAC = a1*sqrt(3.0*sqrtZr4r2A*sqrtZr4r2C) + a2;
          cBC = a1*sqrt(3.0*sqrtZr4r2B*sqrtZr4r2C) + a2;
          eabc += EABC(c6AB,c6AC,c6BC,cAB,cAC,cBC,RAB2,RAC2,RBC2,alpha);             //part of equation 32
        }
      }                 //idBtm
      posA += refa;
    }                   //idAtm
    eabc *= s9D4(this->Type());
  }
  //energy functions
  void EElectrostaticDispersion(matrixE & dens, matrixE & geometry) {
    //calculation of several energy contributions; the D4 contribution is a simplification of function ED4 in Intermolecular.hpp
    double e01 = 0.0;
    double e02 = 0.0;
    double e11 = 0.0;
    double edouble = 0.0;
    double etriple = 0.0;
    double rAB2;
    double rAB3;
    double rAB5;
    double R0A;
    double R0B;
    double f3;
    double f5;
    double sqrtRAB2;
    double a3 = alphaMultipole(1);
    double a5 = alphaMultipole(2);
    double muA2;
    double thetaA2;
    double qA;
    double qB;
    double ed;
    double eq;
    double edd;
    int atomA;
    int kl;
    int ipos = 0;
    int nshellsA;
    int LorbA;
    int ishell = 1;
    bool transitionmetalA;
    ed4 = 0.0;
    //get the charges, dipoles and quadrupoles
    CAMM(dens);
    eaxc = 0.0;
    auxMAT3 = VAB*QShells;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      atomA = atoms[idAtm];
      //dispersion part
      for (size_t irefA = 0; irefA < refA[idAtm]; ++irefA,++ipos) {
        muA2 = 0.0;           //get the element WDisp(i,j)zetaA(j)
        for (size_t idcontr = 0; idcontr < dispsz; ++idcontr) {
          muA2 += WDisp(ipos + 1,idcontr + 1)*zetaA[idcontr];
        }
        ed4 += muA2*zetaA[ipos];
      }
      //energyIES_IXC part
      transitionmetalA = ((atomA > 20)&&(atomA < 30))||((atomA > 38)&&(atomA < 48))||((atomA > 56)&&(atomA < 80));
      nshellsA = AtomShells[idAtm];
      qA = HubbardDer(atomA);
      LorbA = 2*int(transitionmetalA);
      for (int ishellA = 0; ishellA < nshellsA; ++ishellA, ++ishell) {
        muA2 = QShells(ishell,1);
        edouble += muA2*auxMAT3(ishell,1);
        etriple += muA2*muA2*muA2*qA*Korb(LorbA);                                  //HubbardDer(atomA)*Korb(LorbA) is the {\Gamma}_{A,l} of equation 23
        ++LorbA;
        if ((transitionmetalA)&&(ishellA == 0)) {LorbA = 0;}
      }
      //AXC part
      muA2 = 0.0;
      thetaA2 = 0.0;
      for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
        rA[idcoord - 1] = geometry(idAtm + 1,idcoord);
        muA[idcoord - 1] = muAtoms(idAtm + 1,idcoord);
        muA2 += muA[idcoord - 1]*muA[idcoord - 1];
        thetaTemp[idcoord - 1] = thetaAtoms(idAtm + 1,idcoord);
        thetaA2 += thetaTemp[idcoord - 1]*thetaTemp[idcoord - 1];
        thetaA2 += thetaTemp[idcoord - 1]*thetaTemp[idcoord - 1]*double(idcoord == 2);
        thetaTemp[idcoord + 2] = thetaAtoms(idAtm + 1,idcoord + 3);
        thetaA2 += thetaTemp[idcoord + 2]*thetaTemp[idcoord + 2];
        thetaA2 += thetaTemp[idcoord + 2]*thetaTemp[idcoord + 2]*double(idcoord != 3);
      }
      eaxc += fXCmuA(atomA)*muA2 + thetaA2*fXCthetaA(atomA);                       //equation 31
      //AES part
      qA = QAtoms[idAtm];
      R0A = CNrad[idAtm];
      for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
        muB[0] = muAtoms(idBtm + 1,1);
        muB[1] = muAtoms(idBtm + 1,2);
        muB[2] = muAtoms(idBtm + 1,3);
        rAB[0] = (geometry(idBtm + 1,1) - rA[0])*dist_Angstrom2aum1;
        rAB[1] = (geometry(idBtm + 1,2) - rA[1])*dist_Angstrom2aum1;
        rAB[2] = (geometry(idBtm + 1,3) - rA[2])*dist_Angstrom2aum1;
        rAB2 = rAB[0]*rAB[0] + rAB[1]*rAB[1] + rAB[2]*rAB[2];
        sqrtRAB2 = sqrt(rAB2);
        qB = QAtoms[idBtm];
        R0B = CNrad[idBtm];
        ed = 0.0;
        eq = 0.0;
        edd = 0.0;
        for (size_t idcoord1 = 0; idcoord1 < 3; ++idcoord1) {
          ed += qB*muA[idcoord1]*rAB[idcoord1] - qA*muB[idcoord1]*rAB[idcoord1];   //equation 25b
          for (size_t idcoord2 = 0; idcoord2 < 3; ++idcoord2) {
            kl = Index(idcoord1,idcoord2);
            muA2 = rAB[idcoord2]*rAB[idcoord1];
            eq += qB*thetaTemp[kl]*muA2;
            eq += thetaAtoms(idBtm + 1,kl + 1)*qA*muA2;
            edd -= 3.0*muB[idcoord1]*muA[idcoord2]*muA2;
          }
          edd += muB[idcoord1]*muA[idcoord1]*rAB2;
        }
        rAB3 = sqrtRAB2*rAB2;
        rAB5 = rAB3*rAB2;
        f3 = MultipoleDamp(sqrtRAB2,0.5*(R0A + R0B),a3)/rAB3;                      //equation 28
        f5 = MultipoleDamp(sqrtRAB2,0.5*(R0A + R0B),a5)/rAB5;                      //equation 28
        e01 += ed*f3;
        e02 += eq*f5;
        e11 += edd*f5;
      }
    }
    ed4 *= 0.5;
    energyIES_IXC = 0.5*edouble + etriple/3.0;                                     //equation 20
    eaes = e01 + e02 + e11;                                                        //equation 25a
  }
  void EHueckel(matrixE & dens) {
    //calculation of the Hueckel type of energy, equation 15
    energyHueckel = trace(Hcore,dens);
  }
  void EnergyDispersion(matrixE & geometry, matrixE & Edisp, double d2threshold = 1600.0) {
    //calculation of the dispersion part of the energy in matrix form
    matrixE WDisp2Body;
    matrixE WDisp3Body;
    std::vector<double> zetaA2Body;
    std::vector<double> zetaA3Body;
    double ZeffA;
    double etaA;
    double ga = GA();
    double gc = GC();
    double s9d4 = s9D4(this->Type())/3.0;
    double sqrtZr4r2A;
    double sqrtZr4r2B;
    double sqrtZr4r2C;
    double c6AB;
    double c6AC;
    double c6BC;
    double cAB;
    double cAC;
    double cBC;
    double a1 = A1D4(this->Type());
    double a2 = A2D4(this->Type());
    double RAB2;
    double RAC2;
    double RBC2;
    double alpha = AlphaD4(this->Type(),9);
    double e3;
    int ipos = 0;
    int idcontr = 0;
    int ipossav;
    int atomA;
    int atomB;
    int atomC;
    int refa;
    int refb;
    int posA = 0;
    int posB = 0;
    //resize
    Edisp.resize(Natoms,Natoms);
    //zero
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      Edisp(idAtm + 1,idAtm + 1) = 0.0;
      for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
        Edisp(idAtm + 1,idBtm + 1) = 0.0;
        Edisp(idBtm + 1,idAtm + 1) = 0.0;
      }
    }
    //get the weights for 2-body
    getWeightedGFN2(WDisp2Body,auxMAT2,GaussW,atoms,geometry,refA,auxv,extractC6,thresholdzero,d2threshold,gwthresh);
    dispsz = GaussW.size();
    zetaA2Body.resize(dispsz);
    WDisp3Body.resize(dispsz,dispsz);
    zetaA3Body.resize(dispsz);
    //get the weights for 3-body
    NCoord(auxv,atoms,geometry,ERFcount,true,d2threshold,false);
    getC6D4(this->Type(),WDisp3Body,GaussW,auxv,refA,atoms,extractC6,thresholdzero);
    //calculate the zetas
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      atomA = atoms[idAtm];
      ZeffA = Zeff(atomA);
      etaA = HardnessD4(atomA);
      RefQ(atomA,auxv,1);
      for (size_t irefA = 0; irefA < refA[idAtm]; ++irefA,++ipos) {
        zetaA2Body[ipos] = 0.0;
        zetaA3Body[ipos] = 0.0;
        if (GaussW[ipos] < gwthresh) {continue;}
        zetaA2Body[ipos] = ZetaD4(ga,etaA*gc,auxv[irefA] + ZeffA,QAtoms[idAtm] + ZeffA);
        zetaA3Body[ipos] = ZetaD4(ga,etaA*gc,auxv[irefA] + ZeffA,ZeffA)*GaussW[ipos];
      }
    }
    ipos = 0;
    //dispersion part
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      atomA = atoms[idAtm];
      refa = refA[idAtm];
      ipossav = ipos;
      idcontr = 0;
      sqrtZr4r2A = R2R4(atomA);
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        rA[idcoord] = geometry(idAtm + 1,idcoord + 1);
      }
      posB = 0;
      for (size_t idBtm = 0; idBtm < Natoms; ++idBtm) {
        ipos = ipossav;
        //2-body dispersion
        for (size_t irefA = 0; irefA < refa; ++irefA,++ipos) {
          for (size_t irefB = 0; irefB < refA[idBtm]; ++irefB) {
            ZeffA = 0.5*zetaA2Body[ipos]*WDisp2Body(ipos + 1,idcontr + irefB + 1)*zetaA2Body[idcontr + irefB];
            Edisp(idAtm + 1,idBtm + 1) += ZeffA;
            Edisp(idBtm + 1,idAtm + 1) += ZeffA;
          }
        }
        idcontr += refA[idBtm];
        if (idBtm < idAtm) {
          atomB = atoms[idBtm];
          refb = RefN(atomB);
          c6AB = 0.0;
          for (size_t idrefA = 0; idrefA < refa; ++idrefA) {
            for (size_t idrefB = 0; idrefB < refb; ++idrefB) {
              c6AB += zetaA3Body[posA + idrefA]*WDisp3Body(posA + idrefA + 1,posB + idrefB + 1)*zetaA3Body[posB + idrefB];
            }
          }
          WDisp3Body(idAtm + 1,idBtm + 1) = c6AB;
          WDisp3Body(idBtm + 1,idAtm + 1) = c6AB;
          posB += refb;
          sqrtZr4r2B = R2R4(atomB);
          RAB2 = 0.0;
          for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
            rB[idcoord] = geometry(idBtm + 1,idcoord + 1);
            RAB2 += (rB[idcoord] - rA[idcoord])*(rB[idcoord] - rA[idcoord]);
          }
          RAB2 *= dist_Angstrom2aum1*dist_Angstrom2aum1;                            //convert to Bohr
          if (RAB2 > d2threshold) {continue;}
          c6AB = WDisp3Body(idAtm + 1,idBtm + 1);
          cAB = a1*sqrt(3.0*sqrtZr4r2A*sqrtZr4r2B) + a2;
          //3-body dispersion
          for (size_t idCtm = 0; idCtm < idBtm; ++idCtm) {
            atomC = atoms[idCtm];
            sqrtZr4r2C = R2R4(atomC);
            RAC2 = 0.0;
            RBC2 = 0.0;
            for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
              RAC2 += (geometry(idCtm + 1,idcoord + 1) - rA[idcoord])*(geometry(idCtm + 1,idcoord + 1) - rA[idcoord]);
              RBC2 += (geometry(idCtm + 1,idcoord + 1) - rB[idcoord])*(geometry(idCtm + 1,idcoord + 1) - rB[idcoord]);
            }
            RAC2 *= dist_Angstrom2aum1*dist_Angstrom2aum1;         //convert to Bohr
            RBC2 *= dist_Angstrom2aum1*dist_Angstrom2aum1;         //convert to Bohr
            if (RAC2 > d2threshold) {continue;}
            if (RBC2 > d2threshold) {continue;}
            c6AC = WDisp3Body(idAtm + 1,idCtm + 1);
            c6BC = WDisp3Body(idBtm + 1,idCtm + 1);
            cAC = a1*sqrt(3.0*sqrtZr4r2A*sqrtZr4r2C) + a2;
            cBC = a1*sqrt(3.0*sqrtZr4r2B*sqrtZr4r2C) + a2;
            e3 = s9d4*EABC(c6AB,c6AC,c6BC,cAB,cAC,cBC,RAB2,RAC2,RBC2,alpha);
            Edisp(idAtm + 1,idBtm + 1) += e3;
            Edisp(idBtm + 1,idAtm + 1) += e3;
            Edisp(idAtm + 1,idCtm + 1) += e3;
            Edisp(idCtm + 1,idAtm + 1) += e3;
            Edisp(idBtm + 1,idCtm + 1) += e3;
            Edisp(idCtm + 1,idBtm + 1) += e3;
          }
        }
      }                 //idBtm
      posA += refa;
    }                   //idAtm
  }
  void Enuclear(matrixE & geometry, matrixE & Emat, double d2threshold = 1600.0) {
    //function calculating the repulsion energy in matrix form
    enuc = 0.0;
    double alphaA;
    double alphaB;
    double yeffectiveA;
    double yeffectiveB;
    double krep;
    double erepAB;
    double RAB2;
    int atomA;
    int atomB;
    Emat.resize(Natoms,Natoms);
    //zero
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      Emat(idAtm + 1,idAtm + 1) = 0.0;
      for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
        Emat(idAtm + 1,idBtm + 1) = 0.0;
      }
    }
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      atomA = atoms[idAtm];
      alphaA = sqrt(AlphaGFN2(atomA));
      yeffectiveA = Yeff(atomA);
      for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
        atomB = atoms[idBtm];
        alphaB = sqrt(AlphaGFN2(atomB));
        yeffectiveB = Yeff(atomB);
        RAB2 = 0.0;
        krep = 1.0;
        if ((atomA > 2)||(atomB > 2)) {krep += 0.5;}
        krep *= 0.5;
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          RAB2 += (geometry(idBtm + 1,idcoord + 1) - geometry(idAtm + 1,idcoord + 1))*(geometry(idBtm + 1,idcoord + 1) - geometry(idAtm + 1,idcoord + 1));
        }
        RAB2 *= dist_Angstrom2aum1*dist_Angstrom2aum1;
        if (RAB2 > d2threshold) {continue;}
        erepAB = yeffectiveA*yeffectiveB*exp(-alphaA*alphaB*pow(RAB2,krep))/sqrt(RAB2);
        enuc += erepAB;
        Emat(idAtm + 1,idBtm + 1) += erepAB;
        Emat(idBtm + 1,idAtm + 1) += erepAB;
      }                 //idBtm
    }                   //idAtm
  }
  void PrepareEnergyCalculation(matrixE & geometry) {
    //function that triggers a few things to be calculated so that energy is properly evalutated
    CAMM(dens);
    getWeightedGFN2(WDisp,auxMAT2,GaussW,atoms,geometry,refA,auxv,extractC6,thresholdzero,1600.0,gwthresh);
    calcVSDQ(geometry);
  }
  void EnergyMatrix(matrixE & geometry, matrixE & Emat, bool redoCAMM = false, int cleanmatrix = -1) {
    //calculation of all energy contributions but the dispersion in atom-pair form
    int kl;
    int atomA;
    int atomB;
    int nshellsA;
    int nshellsB;
    int LorbA;
    int LorbB;
    int ishell = 1;
    int jshell = 1;
    int countA = 0;
    int countB = 0;
    bool transitionmetalA;
    bool transitionmetalB;
    double rAB2;
    double rAB3;
    double rAB5;
    double R0A;
    double R0B;
    double sqrtRAB2;
    double a3 = alphaMultipole(1);
    double a5 = alphaMultipole(2);
    double qA;
    double qB;
    double ed;
    double eq;
    double edd;
    double muA2;
    double thetaA2;
    double onethird = 1.0/3.0;
    Emat.resize(Natoms,Natoms);
    //zero
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      Emat(idAtm + 1,idAtm + 1) = 0.0;
      for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
        Emat(idAtm + 1,idBtm + 1) = 0.0;
        Emat(idBtm + 1,idAtm + 1) = 0.0;
      }
    }
    //get the charges, dipoles and quadrupoles
    if (redoCAMM) {CAMM(dens);}
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      atomA = atoms[idAtm];
      transitionmetalA = ((atomA > 20)&&(atomA < 30))||((atomA > 38)&&(atomA < 48))||((atomA > 56)&&(atomA < 80));
      nshellsA = AtomShells[idAtm];
      //IES, IXC and Hueckel parts
      qA = HubbardDer(atomA);
      LorbA = 2*int(transitionmetalA);
      for (int ishellA = 0; ishellA < nshellsA; ++ishellA, ++ishell) {
        //IES-IXC
        Emat(idAtm + 1,idAtm + 1) += onethird*QShells(ishell,1)*QShells(ishell,1)*QShells(ishell,1)*qA*Korb(LorbA);
        jshell = 1;
        countB = 0;
        for (size_t idBtm = 0; idBtm < Natoms; ++idBtm) {
          atomB = atoms[idBtm];
          transitionmetalB = ((atomB > 20)&&(atomB < 30))||((atomB > 38)&&(atomB < 48))||((atomB > 56)&&(atomB < 80));
          nshellsB = AtomShells[idBtm];
          LorbB = 2*int(transitionmetalB);
          for (int ishellB = 0; ishellB < nshellsB; ++ishellB, ++jshell) {
            //IES-IXC
            Emat(idAtm + 1,idBtm + 1) += 0.5*QShells(ishell,1)*VAB(ishell,jshell)*QShells(jshell,1);
            //Hueckel
            for (size_t idAO = 0; idAO < 2*LorbA + 1; ++idAO) {
              for (size_t idBO = 0; idBO < 2*LorbB + 1; ++idBO) {
                Emat(idAtm + 1,idBtm + 1) += Hcore(countA + idAO + 1,countB + idBO + 1)*dens(countB + idBO + 1,countA + idAO + 1);
              }
            }
            countB += 2*LorbB + 1;
            ++LorbB;
            if ((transitionmetalB)&&(ishellB == 0)) {LorbB = 0;}
          }
        }
        countA += 2*LorbA + 1;
        ++LorbA;
        if ((transitionmetalA)&&(ishellA == 0)) {LorbA = 0;}
      }
      //AXC part
      muA2 = 0.0;
      thetaA2 = 0.0;
      for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
        rA[idcoord - 1] = geometry(idAtm + 1,idcoord);
        muA[idcoord - 1] = muAtoms(idAtm + 1,idcoord);
        muA2 += muA[idcoord - 1]*muA[idcoord - 1];
        thetaTemp[idcoord - 1] = thetaAtoms(idAtm + 1,idcoord);
        thetaA2 += thetaTemp[idcoord - 1]*thetaTemp[idcoord - 1];
        thetaA2 += thetaTemp[idcoord - 1]*thetaTemp[idcoord - 1]*double(idcoord == 2);
        thetaTemp[idcoord + 2] = thetaAtoms(idAtm + 1,idcoord + 3);
        thetaA2 += thetaTemp[idcoord + 2]*thetaTemp[idcoord + 2];
        thetaA2 += thetaTemp[idcoord + 2]*thetaTemp[idcoord + 2]*double(idcoord != 3);
      }
      Emat(idAtm + 1,idAtm + 1) += fXCmuA(atomA)*muA2 + thetaA2*fXCthetaA(atomA);
      //AES part
      qA = QAtoms[idAtm];
      R0A = CNrad[idAtm];
      for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
        muB[0] = muAtoms(idBtm + 1,1);
        muB[1] = muAtoms(idBtm + 1,2);
        muB[2] = muAtoms(idBtm + 1,3);
        rAB[0] = (geometry(idBtm + 1,1) - rA[0])*dist_Angstrom2aum1;
        rAB[1] = (geometry(idBtm + 1,2) - rA[1])*dist_Angstrom2aum1;
        rAB[2] = (geometry(idBtm + 1,3) - rA[2])*dist_Angstrom2aum1;
        rAB2 = rAB[0]*rAB[0] + rAB[1]*rAB[1] + rAB[2]*rAB[2];
        sqrtRAB2 = sqrt(rAB2);
        qB = QAtoms[idBtm];
        R0B = CNrad[idBtm];
        ed = 0.0;
        eq = 0.0;
        edd = 0.0;
        for (size_t idcoord1 = 0; idcoord1 < 3; ++idcoord1) {
          ed += qB*muA[idcoord1]*rAB[idcoord1] - qA*muB[idcoord1]*rAB[idcoord1];   //equation 25b
          for (size_t idcoord2 = 0; idcoord2 < 3; ++idcoord2) {
            kl = Index(idcoord1,idcoord2);
            muA2 = rAB[idcoord2]*rAB[idcoord1];
            eq += qB*thetaTemp[kl]*muA2;
            eq += thetaAtoms(idBtm + 1,kl + 1)*qA*muA2;
            edd -= 3.0*muB[idcoord1]*muA[idcoord2]*muA2;
          }
          edd += muB[idcoord1]*muA[idcoord1]*rAB2;
        }
        rAB3 = sqrtRAB2*rAB2;
        rAB5 = rAB3*rAB2;
        Emat(idAtm + 1,idBtm + 1) += ed*MultipoleDamp(sqrtRAB2,0.5*(R0A + R0B),a3)/rAB3 + (eq + edd)*MultipoleDamp(sqrtRAB2,0.5*(R0A + R0B),a5)/rAB5;
      }
    }
    //symmetrize
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
        Emat(idAtm + 1,idBtm + 1) += Emat(idBtm + 1,idAtm + 1);
        Emat(idBtm + 1,idAtm + 1) = Emat(idAtm + 1,idBtm + 1);
      }
    }
  }
  void calcEnergy(matrixE & geom) {
    //function that triggers the calculation of energy
    EElectrostaticDispersion(dens,geom);
    EHueckel(dens);
    if (usesolvation) {esolv = SolvationEnergy(QAtoms);}
    ErepATM(geom,1600.0);
    etot = eaxc + eaes + energyIES_IXC + ed4 + energyHueckel + efermia + efermib + esolv + enuc + eabc;
  }
  //Fermi smearing
  void GFermi(double thresh = 1.0e-9, double expthresh = 50.0, double maxiter = 200) {
    if (Telectron > 0.1) {
      Occ();
      size_t homoa = 0;
      size_t homob = 0;
      for (size_t idorb = 0; idorb < NAOs; ++idorb) {
        if (occa[idorb] > 0.99) {homoa = idorb;}
        if (occb[idorb] > 0.99) {homob = idorb;}
      }
      double KBp = au2J/(Telectron*KB);                     //this is kB*T
      double occta = double(homoa + 1.0);
      double occtb = double(homob + 1.0);
      double total_numbera;
      double total_numberb;
      double total_dfermia;
      double total_dfermib;
      double exponenta;
      double exponentb;
      double change_fermia;
      double change_fermib;
      efermia = 0.5*(EMOs[homoa + 1] + EMOs[homoa]);
      efermib = 0.5*(EMOs[homob + 1] + EMOs[homob]);
      int conva = 0;
      int convb = 0;
      //tune the occupation numbers
      for (size_t iter = 0; iter < maxiter; ++iter) {
        total_numbera = 0.0;
        total_numberb = 0.0;
        total_dfermia = 0.0;
        total_dfermib = 0.0;
        for (size_t idorb = 0; idorb < NAOs; ++idorb) {
          if (conva == 0) {
            exponenta = (EMOs[idorb] - efermia)*KBp;
            if ((exponenta < expthresh)) {
              exponentb = exp(exponenta);
              exponenta = 1.0/(exponentb + 1.0);
              occa[idorb] = exponenta;                        //equation 11
              total_numbera += exponenta;
              total_dfermia += KBp*exponentb*exponenta*exponenta;
            }
            else {
              occa[idorb] = 0.0;
            }
          }
          if (convb == 0) {
            exponentb = (EMOs[idorb] - efermib)*KBp;
            if ((exponentb < expthresh)) {
              exponenta = exp(exponentb);
              exponentb = 1.0/(exponenta + 1.0);
              occb[idorb] = exponentb;                        //equation 11
              total_numberb += exponentb;
              total_dfermib += KBp*exponenta*exponentb*exponentb;
            }
            else {
              occb[idorb] = 0.0;
            }
          }
        }
        if (fabs(total_dfermia) < thresholdzero) {total_dfermia = 1.0;}
        if (fabs(total_dfermib) < thresholdzero) {total_dfermib = 1.0;}
        change_fermia = (occta - total_numbera)/total_dfermia;
        change_fermib = (occtb - total_numberb)/total_dfermib;
        efermia += change_fermia;
        efermib += change_fermib;
        if (fabs(occta - total_numbera) < thresh) {conva = 1;}
        if (fabs(occtb - total_numberb) < thresh) {convb = 1;}
        if ((conva == 1)&&(convb == 1)) {break;}
      }
      efermia = 0.0;
      efermib = 0.0;
      convb = occb.size();
      conva = occa.size() - convb;
      for (size_t idorb = 0; idorb < convb; ++idorb) {        //equation 10
        if ((occa[idorb] > thresh)&&(1.0 - occa[idorb] > thresh)) {efermia += occa[idorb]*log(occa[idorb]) + (1.0 - occa[idorb])*log(1.0 - occa[idorb]);}
        if ((occb[idorb] > thresh)&&(1.0 - occb[idorb] > thresh)) {efermib += occb[idorb]*log(occb[idorb]) + (1.0 - occb[idorb])*log(1.0 - occb[idorb]);}
      }
      for (size_t idorb = 0; idorb < conva; ++idorb) {        //equation 10
        if ((occa[idorb + convb] > thresh)&&(1.0 - occa[idorb + convb] > thresh)) {efermia += occa[idorb + convb]*log(occa[idorb + convb]) + (1.0 - occa[idorb + convb])*log(1.0 - occa[idorb + convb]);}
      }
      efermia /= KBp;
      efermib /= KBp;
      setOccupancy();
    }
  }
  //functions required to build the Hamiltonian
  void CalcAtomShells() {
    //function that calculates the number of shells per atom
    int nshells;
    int nAOs;
    AtomShells.resize(Natoms);
    AOs = basis.AtomNAOs(atoms);
    totalshells = 0;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      nshells = 1;
      nAOs = AOs[iatm];
      if ((nAOs > 1)&&(nAOs < 5)) {nshells = 2;}
      else if (nAOs > 4) {nshells = 3;}
      totalshells += nshells;
      AtomShells[iatm] = nshells;
    }
  }
  void initQShell(std::vector<double> & qvec) {
    bool transitionmetalA;
    int atomA;
    int nshellsA;
    int LorbA;
    int ishell = 1;
    double ntot;
    double zshell;
    double zcore;
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      atomA = atoms[idatm];
      transitionmetalA = ((atomA > 20)&&(atomA < 30))||((atomA > 38)&&(atomA < 48))||((atomA > 56)&&(atomA < 80));
      nshellsA = AtomShells[idatm];
      LorbA = 2*int(transitionmetalA);
      ntot = -1.0e-6;
      for (size_t idshellA = 0; idshellA < nshellsA; ++idshellA,++ishell) {
        zshell = AtomOccupancy(atomA,LorbA);
        ntot += zshell;
        zcore = double(atomA - CoreCharge[idatm]);
        if (ntot > zcore) {zshell = 0.0;}
        zshell *= QAtoms[idatm]/zcore;
        QShells(ishell,1) = zshell;
        qvec[ishell - 1] = zshell;
        ++LorbA;
        if ((transitionmetalA)&&(idshellA == 0)) {LorbA = 0;}
      }
    }
  }
  void CAMM(matrixE & dens) {
    //function that calculates the Cumulative Atomic Multipole Moments each atom
    bdens = dens;
    basis.CAMM(bdens,QShells,muAtoms,thetaAtoms);                            //equations 27
    //now add the reference populations and sum to atomic charges
    double tii;
    bool transitionmetalA;
    int atomA;
    int nshellsA;
    int LorbA;
    int ishell = 1;
    //build atomic charges and remove trace of quadrupole moment
    for (int iAtm = 0; iAtm < Natoms; ++iAtm) {
      atomA = atoms[iAtm];
      transitionmetalA = ((atomA > 20)&&(atomA < 30))||((atomA > 38)&&(atomA < 48))||((atomA > 56)&&(atomA < 80));
      nshellsA = AtomShells[iAtm];
      LorbA = 2*int(transitionmetalA);
      //QA
      QAtoms[iAtm] = 0.0;
      for (int ishellA = 0; ishellA < nshellsA; ++ishellA, ++ishell) {
        QShells(ishell,1) += AtomOccupancy(atomA,LorbA);
        QAtoms[iAtm] += QShells(ishell,1);
        ++LorbA;
        if ((transitionmetalA)&&(ishellA == 0)) {LorbA = 0;}
      }
      //traceless quadrupoles, equation 26
      tii = 0.5*(thetaAtoms(iAtm + 1,1) + thetaAtoms(iAtm + 1,3) + thetaAtoms(iAtm + 1,6));
      for (size_t idcomp = 1; idcomp < 7; ++idcomp) {
        thetaAtoms(iAtm + 1,idcomp) *= 1.50;
      }
      thetaAtoms(iAtm + 1,1) -= tii;
      thetaAtoms(iAtm + 1,3) -= tii;
      thetaAtoms(iAtm + 1,6) -= tii;
    }
  }
  void CalcGammaAB(matrixE & geometry) {
    //function calculating the distance dependent Coulomb interaction for the isotropic electrostatic and exchange term
    //this is the matrix for for equation 21, where rows and columns represent shells of the atom (each element is a shell in an atom)
    VAB.resize(totalshells,totalshells);
    double RAB;
    double RAB2;
    double gammaAB;
    double eta;
    double etaA;
    double etaB;
    double kscaleA;
    double kscaleB;
    bool transitionmetalA;
    bool transitionmetalB;
    int atomA;
    int atomB;
    int nshellsA;
    int nshellsB;
    int LorbA;
    int LorbB;
    int irow = 1;
    int icol = 1;
    for (size_t iAtm = 0; iAtm < Natoms; ++iAtm) {
      atomA = atoms[iAtm];
      transitionmetalA = ((atomA > 20)&&(atomA < 30))||((atomA > 38)&&(atomA < 48))||((atomA > 56)&&(atomA < 80));
      nshellsA = AtomShells[iAtm];
      etaA = Hardness(atomA);
      icol = 1;
      for (size_t iBtm = 0; iBtm < iAtm; ++iBtm) {
        atomB = atoms[iBtm];
        transitionmetalB = ((atomB > 20)&&(atomB < 30))||((atomB > 38)&&(atomB < 48))||((atomB > 56)&&(atomB < 80));
        nshellsB = AtomShells[iBtm];
        etaB = Hardness(atomB);
        RAB2 = 0.0;
        for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
          RAB = geometry(iAtm + 1,idcoord) - geometry(iBtm + 1,idcoord);
          RAB2 += RAB*RAB;
        }
        RAB = sqrt(RAB2)*dist_Angstrom2aum1;
        LorbA = 2*int(transitionmetalA);
        for (size_t ishellA = 0; ishellA < nshellsA; ++ishellA) {
          kscaleA = HardnessLScal(atomA,LorbA);
          LorbB = 2*int(transitionmetalB);
          for (size_t ishellB = 0; ishellB < nshellsB; ++ishellB) {
            kscaleB = HardnessLScal(atomB,LorbB);
            eta = 0.5*(etaA*kscaleA + etaB*kscaleB);                             //equation 22 of the paper, note that my kappa constants are shifted by 1
            RAB2 = 1.0/(eta*eta);
            gammaAB = 1.0/sqrt(RAB*RAB + RAB2);                                  //equation 21 of the paper
            VAB(irow + ishellA,icol + ishellB) = gammaAB;
            VAB(icol + ishellB,irow + ishellA) = gammaAB;
            ++LorbB;
            if ((transitionmetalB)&&(ishellB == 0)) {LorbB = 0;}
          }
          ++LorbA;
          if ((transitionmetalA)&&(ishellA == 0)) {LorbA = 0;}
        }
        icol += nshellsB;
      }
      LorbA = 2*int(transitionmetalA);
      for (size_t ishellA1 = 0; ishellA1 < nshellsA; ++ishellA1) {
        kscaleA = HardnessLScal(atomA,LorbA);
        LorbB = 2*int(transitionmetalA);
        for (size_t ishellA2 = 0; ishellA2 < ishellA1; ++ishellA2) {
          kscaleB = HardnessLScal(atomA,LorbB);
          eta = 0.5*etaA*(kscaleA + kscaleB);
          VAB(irow + ishellA1,irow + ishellA2) = eta;                            //equation 21 of the paper, RAB = 0
          VAB(irow + ishellA2,irow + ishellA1) = eta;                            //equation 21 of the paper, RAB = 0
          ++LorbB;
          if ((transitionmetalA)&&(ishellA2 == 0)) {LorbB = 0;}
        }
        ++LorbA;
        if ((transitionmetalA)&&(ishellA1 == 0)) {LorbA = 0;}
        VAB(irow + ishellA1,irow + ishellA1) = etaA*kscaleA;                     //equation 21 of the paper, RAB = 0
      }
      irow += nshellsA;
    }
    if (print == 2) {
      std::cout << "Gamma matrix" << std::endl;
      VAB.Print();
    }
  }
  double shellPoly(int atomA, int LorbA, std::vector<double> & rA, int atomB, int LorbB, std::vector<double> & rB) {
    //function calculating the shell polynomial for GFN2-xTB, equation 19 of the paper
    double rcovA = (rA[0] - rB[0]);
    double rcovB = (rA[1] - rB[1]);
    double aux = (rA[2] - rB[2]);
    double distance = sqrt(rcovA*rcovA + rcovB*rcovB + aux*aux);
    double kpolyA = ShellPolynomialCoeff(atomA,LorbA);
    double kpolyB = ShellPolynomialCoeff(atomB,LorbB);
    rcovA = AtmRadii(atomA);
    rcovB = AtmRadii(atomB);
    aux = sqrt(distance/(rcovA + rcovB));
    return (1.0 + 0.01*kpolyA*aux)*(1.0 + 0.01*kpolyB*aux);
  }
  double dshellPoly(std::vector<double> & rAB, std::vector<double> & drABpoly, double rAB2, double polyA, double polyB, double RcA, double RcB) {
    //derivative of the shellPoly factor
    double RAB = sqrt(rAB2);
    double rcov = (RcA + RcB)*dist_Angstrom2aum1;
    RAB /= rcov;
    polyA *= 0.01;
    polyB *= 0.01;
    RcA = sqrt(RAB);
    RcB = polyA*RcA;
    rcov = 1.0/rAB2;
    RAB = polyB*RcA + 1.0;
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      drABpoly[idcoord] = 0.5*(RcB*RAB + (1.0 + RcB)*polyB*RcA)*rcov*rAB[idcoord];
    }
    return (1.0 + RcB)*(1.0 + polyB*RcA);
  }
  void calcHcore(matrixE & sao, matrixE & geometry, double thresholdz = 1.0e-8) {
    //function that calculates the core Hamiltonian, i.e. equations 16 and 17
    Hcore.resize(NAOs,NAOs);
    //get some data
    size_t naosA;
    size_t naosB;
    int atomA;
    int atomB;
    int nshellsA;
    int nshellsB;
    int LorbA;
    int LorbB;
    int icntauxA;
    int icntauxB;
    int icntA = 0;
    int icntB = 0;
    bool transitionmetalA;
    bool transitionmetalB;
    double auxiliary;
    double kll;
    double hii;
    double hjj;
    double Sij;
    double zetaij;
    double ENA;
    double ENB;
    double den;
    double ken = kEN();
    double shellpolynomial;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      //get data for atom A
      atomA = atoms[idAtm];
      transitionmetalA = ((atomA > 20)&&(atomA < 30))||((atomA > 38)&&(atomA < 48))||((atomA > 56)&&(atomA < 80));
      nshellsA = AtomShells[idAtm];
      naosA = nshellsA*nshellsA;
      alphavec = basis.GetExponents(idAtm + 1);
      ENA = PaulingElectronegativity(atomA);
      rA[0] = geometry(idAtm + 1,1);
      rA[1] = geometry(idAtm + 1,2);
      rA[2] = geometry(idAtm + 1,3);
      //getting the main diagonal elements for atom A
      icntauxA = 0;
      LorbA = 2*int(transitionmetalA);
      for (size_t ishellA = 0; ishellA < nshellsA; ++ishellA) {
        icntauxB = 0;
        LorbB = 2*int(transitionmetalA);
        for (size_t ishellB = 0; ishellB < ishellA + 1; ++ishellB) {
          for (size_t iorb = 1; iorb < 2*LorbA + 2; ++iorb) {
            for (size_t jorb = 1; jorb < 2*LorbB + 2; ++jorb) {
              Hcore(icntA + icntauxA + iorb,icntA + icntauxB + jorb) = 0.0;
              Hcore(icntA + icntauxB + jorb,icntA + icntauxA + iorb) = 0.0;
            }
          }
          icntauxB += 2*LorbB + 1;
          ++LorbB;
          if ((transitionmetalA)&&(ishellB == 0)) {LorbB = 0;}
        }
        auxiliary = EnergyLevel(atomA,LorbA) - HCNl(atomA,LorbA)*CoordNum[idAtm];     //equation 17 of paper
        for (size_t iorb = 1; iorb < 2*LorbA + 2; ++iorb) {
          Hcore(icntA + icntauxA + iorb,icntA + icntauxA + iorb) = auxiliary;
        }
        icntauxA += 2*LorbA + 1;
        ++LorbA;
        if ((transitionmetalA)&&(ishellA == 0)) {LorbA = 0;}
      }
      icntB = 0;
      for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
        //get data for atom B
        atomB = atoms[idBtm];
        transitionmetalB = ((atomB > 20)&&(atomB < 30))||((atomB > 38)&&(atomB < 48))||((atomB > 56)&&(atomB < 80));
        nshellsB = AtomShells[idBtm];
        naosB = nshellsB*nshellsB;
        betavec = basis.GetExponents(idBtm + 1);
        ENB = PaulingElectronegativity(atomB);
        rB[0] = geometry(idBtm + 1,1);
        rB[1] = geometry(idBtm + 1,2);
        rB[2] = geometry(idBtm + 1,3);
        //pair parameters
        den = (ENA - ENB)*(ENA - ENB);
        LorbA = 2*int(transitionmetalA);
        icntauxA = 0;
        for (size_t ishellA = 0; ishellA < nshellsA; ++ishellA) {
          hii = Hcore(icntA + icntauxA + 1,icntA + icntauxA + 1);
          LorbB = 2*int(transitionmetalB);
          icntauxB = 0;
          for (size_t ishellB = 0; ishellB < nshellsB; ++ishellB) {
            hjj = Hcore(icntB + icntauxB + 1,icntB + icntauxB + 1);
            zetaij = sqrt(2.0*sqrt(alphavec[ishellA]*betavec[ishellB])/(alphavec[ishellA] + betavec[ishellB]));
            kll = kval(LorbA,LorbB);                                     //k parameter, called kScale in the xtb code
            shellpolynomial = shellPoly(atomA,LorbA,rA,atomB,LorbB,rB);
            auxiliary = 0.5*kll*(hii + hjj)*(1.0 + ken*den);
            for (size_t korb = 1; korb < 2*LorbA + 2; ++korb) {
              for (size_t lorb = 1; lorb < 2*LorbB + 2; ++lorb) {              //equation 16 of paper
                Sij = sao(icntA + icntauxA + korb,icntB + icntauxB + lorb);
                if (fabs(Sij) < thresholdz) {
                  sao(icntA + icntauxA + korb,icntB + icntauxB + lorb) = 0.0;
                  sao(icntB + icntauxB + lorb,icntA + icntauxA + korb) = 0.0;
                  Sij = 0.0;
                }
                Hcore(icntA + icntauxA + korb,icntB + icntauxB + lorb) = auxiliary*Sij*zetaij*shellpolynomial;
                Hcore(icntB + icntauxB + lorb,icntA + icntauxA + korb) = auxiliary*Sij*zetaij*shellpolynomial;
              }
            }
            icntauxB += 2*LorbB + 1;
            ++LorbB;
            if ((transitionmetalB)&&(ishellB == 0)) {LorbB = 0;}
          }                 //ishellB
          icntauxA += 2*LorbA + 1;
          ++LorbA;
          if ((transitionmetalA)&&(ishellA == 0)) {LorbA = 0;}
        }                   //ishellA
        icntB += naosB;
      }                     //idBtm
      icntA += naosA;
    }                       //idAtm
    if (print == 2) {
      std::cout << "Hcore matrix" << std::endl;
      Hcore.Print();
    }
  }
  void calcVSDQ(matrixE & geometry) {
    //function for calculating the potential terms V, proportional to S, D and Q integrals
    double STemp;
    double rA2;
    double rAB2;
    double sqrtrAB2;
    double rABn;
    double t1a;
    double t2a;
    double t3a;
    double t4a;
    double dum3a;
    double dum5a;
    double f3;
    double f5;
    double R0A;
    double R0B;
    double a3 = alphaMultipole(1);
    double a5 = alphaMultipole(2);
    double ga = GA();
    double gc = GC();
    int atomA;
    int ll;
    int ipos = 0;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      atomA = atoms[idAtm];
      //get zeta and dzeta/dq for Fock matrix
      R0A = Zeff(atomA);
      rA2 = HardnessD4(atomA);
      RefQ(atomA,auxv,1);
      for (size_t irefA = 0; irefA < refA[idAtm]; ++irefA,++ipos) {
        zetaA[ipos] = 0.0;
        dzetaA[ipos] = 0.0;
        if (GaussW[ipos] < gwthresh) {continue;}
        zetaA[ipos] = ZetaD4(ga,rA2*gc,auxv[irefA] + R0A,QAtoms[idAtm] + R0A);
        dzetaA[ipos] = dZetaD4dq(ga,rA2*gc,auxv[irefA] + R0A,QAtoms[idAtm] + R0A);
      }
      R0A = CNrad[idAtm];
      rA2 = 0.0;
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        rA[idcoord] = geometry(idAtm + 1,idcoord + 1)*dist_Angstrom2aum1;
        rA2 += rA[idcoord]*rA[idcoord];
        muA[idcoord] = muAtoms(idAtm + 1,idcoord + 1);
        muTemp[idcoord] = 0.0;
        thetaTemp[idcoord] = 0.0;
        thetaTemp[idcoord + 3] = 0.0;
      }
      STemp = 0.0;
      for (size_t idBtm = 0; idBtm < Natoms; ++idBtm) {
        if (idAtm == idBtm) {continue;}
        R0B = CNrad[idBtm];
        rAB2 = 0.0;
        for (size_t idcoord1 = 0; idcoord1 < 3; ++idcoord1) {
          rAB[idcoord1] = rA[idcoord1] - geometry(idBtm + 1,idcoord1 + 1)*dist_Angstrom2aum1;
          rAB2 += rAB[idcoord1]*rAB[idcoord1];
        }
        sqrtrAB2 = sqrt(rAB2);
        rABn = rAB2*sqrtrAB2;
        f3 = MultipoleDamp(sqrtrAB2,0.5*(R0A + R0B),a3)/rABn;           //equation 28
        rABn *= rAB2;
        f5 = MultipoleDamp(sqrtrAB2,0.5*(R0A + R0B),a5)/rABn;           //equation 28
        t1a = 0.0;
        t2a = 0.0;
        t3a = 0.0;
        t4a = 0.0;
        dum5a = 0.0;
        for (size_t idcoord1 = 0; idcoord1 < 3; ++idcoord1) {
          muB[idcoord1] = muAtoms(idBtm + 1,idcoord1 + 1);
          t1a += rA[idcoord1]*rAB[idcoord1];                                //RA.RAB, for dipole-charge and dipole-dipole terms
          t2a += muB[idcoord1]*rAB[idcoord1];                               //muA.RAB, for charge-dipole and dipole-dipole
          t3a += rA[idcoord1]*muB[idcoord1];                                //RA.muA, for diagonal dipole-dipole
          t4a += rAB[idcoord1]*rAB[idcoord1]*rA[idcoord1]*rA[idcoord1];     //square of Hadamard product for equation 42
          for (size_t idcoord2 = 0; idcoord2 < 3; ++idcoord2) {
            // potential from quadrupoles
            ll = lin(idcoord1 + 1,idcoord2 + 1);
            R0B = thetaAtoms(idBtm + 1,ll);
            dum5a -= R0B*rAB[idcoord1]*rAB[idcoord2] + 1.5*QAtoms[idBtm]*rAB[idcoord1]*rAB[idcoord2]*rA[idcoord1]*rA[idcoord2]; //second terms of second and third lines of equation 42
            if (idcoord2 >= idcoord1) {continue;}
            thetaTemp[idcoord1 + idcoord2 + 2] -= 3.0*QAtoms[idBtm]*f5*rAB[idcoord2]*rAB[idcoord1];    //first part of first line in equation 44
          }
          thetaTemp[idcoord1] -= 1.5*QAtoms[idBtm]*f5*rAB[idcoord1]*rAB[idcoord1];                     //first part of first line in equation 44
        }
        dum3a = -t1a*QAtoms[idBtm] - t2a;                                //first and last terms in second line, equation 42
        dum5a += t3a*rAB2 - 3.0*t1a*t2a + 0.5*QAtoms[idBtm]*rA2*rAB2;    //first line + first term of third line in equation 42
        STemp += dum5a*f5 + dum3a*f3;                                    //puts together the first three lines of equation 42
        for (size_t idcoord1 = 0; idcoord1 < 3; ++idcoord1) {
          dum3a = rAB[idcoord1]*QAtoms[idBtm];                                               //third term, first line of equation 43
          dum5a = 3.0*rAB[idcoord1]*t2a - rAB2*muB[idcoord1];                                //first two terms of line 1 in equation 43
          dum5a += -QAtoms[idBtm]*rAB2*rA[idcoord1] + 3.0*QAtoms[idBtm]*rAB[idcoord1]*t1a;   //line two, equation 43
          muTemp[idcoord1] += dum3a*f3 + dum5a*f5;                                           //puts together the first two lines of equation 43
          thetaTemp[idcoord1] += 0.5*rAB2*QAtoms[idBtm]*f5;                                  //second part of first line in equation 44
        }
      }
      VS[idAtm] = STemp;
      for (int idcoord = 1; idcoord < 4; ++idcoord) {
        VD(idAtm + 1,idcoord) = muTemp[idcoord - 1];
        VQ(idAtm + 1,idcoord) = thetaTemp[idcoord - 1];
        VQ(idAtm + 1,idcoord + 3) = thetaTemp[idcoord + 2];
      }
      t1a = 2.0*fXCmuA(atomA);
      t4a = 6.0*fXCthetaA(atomA);
      t2a = 0.0;
      t3a = 0.0;
      for (size_t idcoord1 = 0; idcoord1 < 3; ++idcoord1) {
        t3a += t1a*rA[idcoord1]*muA[idcoord1];                          //first term, fourth line of equation 42
        VD(idAtm + 1,idcoord1 + 1) -= t1a*muA[idcoord1];                //first term, third line, equation 43
        for (size_t idcoord2 = 0; idcoord2 < idcoord1; ++idcoord2) {
          ll = lin(idcoord1 + 1,idcoord2 + 1);
          R0B = thetaAtoms(idAtm + 1,ll);
          VQ(idAtm + 1,idcoord1 + idcoord2 + 3) -= R0B*t4a;             //first part of term in second line in equation 44
          t3a -= rA[idcoord1]*rA[idcoord2]*R0B*t4a;                     //first part, second term, fourth line of equation 42
          VD(idAtm + 1,idcoord1 + 1) += rA[idcoord2]*R0B*t4a;           //first part, second term, third line of equation 43
          VD(idAtm + 1,idcoord2 + 1) += rA[idcoord1]*R0B*t4a;           //first part, second term, third line of equation 43
        }
        ll = lin(idcoord1 + 1,idcoord1 + 1);
        R0B = thetaAtoms(idAtm + 1,ll);
        VQ(idAtm + 1,idcoord1 + 1) -= 0.5*t4a*R0B;                      //first part of term in second line in equation 44
        t3a -= 0.5*rA[idcoord1]*rA[idcoord1]*R0B*t4a;                   //first part, second term, fourth line of equation 42
        VD(idAtm + 1,idcoord1 + 1) += rA[idcoord1]*R0B*t4a;             //first part, second term, third line of equation 43
        t2a += R0B;                                                     //collect trace terms
      }
      VS[idAtm] += t3a;                                                 //puts together the last line of equation 42
      t2a *= fXCthetaA(atomA);
      for (size_t idcoord1 = 0; idcoord1 < 3; ++idcoord1) {             //trace terms in equations 42, 43 and 44
        VQ(idAtm + 1,idcoord1 + 1) += t2a;                              //equation 44
        VD(idAtm + 1,idcoord1 + 1) -= 2.0*rA[idcoord1]*t2a;             //equation 43
        VS[idAtm] += t2a*rA[idcoord1]*rA[idcoord1];                     //equation 42
      }
    }
  }
  void dSQH0(matrixE & grad, matrixE & dens, matrixE & densw, std::vector<double> & VS, matrixE & VES, matrixE & VD, matrixE & VQ, matrixE & dSEdCN, std::vector<double> & dHdCN, matrixE & geometry, double d2threshold = 1600.0, double thresholdzero = 1.0e-8) {
    //function to obtain the gradient contribution from overlap, dipole and quadrupole
    GFN2BSet bst;
    int counterA = 0;
    int counterB = 0;
    int nshellsA;
    int nshellsB;
    int naosdoneA;
    int naosdoneB;
    int norbsA;
    int norbsB;
    int atomA;
    int atomB;
    int LorbA;
    int LorbB;
    int cntA = 0;
    int cntB = 0;
    bool transitionmetalA;
    bool transitionmetalB;
    double rAB2;
    double hii;
    double hjj;
    double haverage;
    double shellpoly;
    double g_xyz[3];
    double dCN;
    double Hij;
    double Pij;
    double PEWij;
    double Sij;
    double HijPij;
    double dtmp;
    double qtmp;
    std::vector<matrixE> SDQ(10);
    matrix<matrixE> gSDQ(3,19);
    matrixE temp(6,6);
    temp.zero();
    for (size_t idx = 0; idx < 10; ++idx) {
      SDQ[idx] = temp;
      for (size_t idr = 1; idr < 4; ++idr) {
        gSDQ(idr,idx + 1) = temp;
        gSDQ(idr,idx + 10) = temp;
      }
    }
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        rA[idcoord] = geometry(idAtm + 1,idcoord + 1)*dist_Angstrom2aum1;
      }
      atomA = atoms[idAtm];
      transitionmetalA = ((atomA > 20)&&(atomA < 30))||((atomA > 38)&&(atomA < 48))||((atomA > 56)&&(atomA < 80));
      nshellsA = AtomShells[idAtm];
      norbsA = nshellsA*nshellsA;
      LorbA = 2*int(transitionmetalA);
      for (size_t idshell = 0; idshell < nshellsA; ++idshell) {
        dSEdCN(idshell + 1,idAtm + 1) = -HCNl(atomA,LorbA)*au2eV;
        ++LorbA;
        if ((transitionmetalA)&&(idshell == 0)) {LorbA = 0;}
      }
      cntB = 0;
      for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
        atomB = atoms[idBtm];
        transitionmetalB = ((atomB > 20)&&(atomB < 30))||((atomB > 38)&&(atomB < 48))||((atomB > 56)&&(atomB < 80));
        nshellsB = AtomShells[idBtm];
        norbsB = nshellsB*nshellsB;
        rAB2 = 0.0;
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          rB[idcoord] = geometry(idBtm + 1,idcoord + 1)*dist_Angstrom2aum1;
          rAB[idcoord] = rA[idcoord] - rB[idcoord];
          rAB2 += rAB[idcoord]*rAB[idcoord];
        }
        if (rAB2 > d2threshold) {
          cntB += norbsB;
          continue;
        }
        LorbA = 2*int(transitionmetalA);
        naosdoneA = 0;
        for (size_t ishellA = 0; ishellA < nshellsA; ++ishellA) {
          hii = (EnergyLevel(atomA,LorbA) - HCNl(atomA,LorbA)*CoordNum[idAtm])*au2eV;
          LorbB = 2*int(transitionmetalB);
          naosdoneB = 0;
          for (size_t ishellB = 0; ishellB < nshellsB; ++ishellB) {
            hjj = (EnergyLevel(atomB,LorbB) - HCNl(atomB,LorbB)*CoordNum[idBtm])*au2eV;
            shellpoly = dshellPoly(rAB,muTemp,rAB2,ShellPolynomialCoeff(atomA,LorbA),ShellPolynomialCoeff(atomB,LorbB),AtmRadii(atomA),AtmRadii(atomB));
            haverage = 0.5*(hii + hjj);
            basis.gOvDipQuadrGTO(bst,SDQ,gSDQ,rA,rB,rAB,rAB2,ishellA,ishellB,idAtm,idBtm,1.0e-8,2000.0,25.0);
            basis.CAO2SAO(gSDQ,temp,LorbB,LorbA);
            g_xyz[0] = 0.0;
            g_xyz[1] = 0.0;
            g_xyz[2] = 0.0;
            dCN = 0.0;
            for (size_t idAO = 0; idAO < 2*LorbA + 1; ++idAO) {
              for (size_t jdAO = 0; jdAO < 2*LorbB + 1; ++jdAO) {
                Pij = dens(cntB + naosdoneB + jdAO + 1,cntA + naosdoneA + idAO + 1);
                Sij = sao(cntB + naosdoneB + jdAO + 1,cntA + naosdoneA + idAO + 1);
                if (fabs(Sij) < thresholdzero) {Hij = 0.0;}
                else {Hij = Hcore(cntB + naosdoneB + jdAO + 1,cntA + naosdoneA + idAO + 1)/Sij;}
                HijPij = Hij*Pij;
                PEWij = densw(cntB + naosdoneB + jdAO + 1,cntA + naosdoneA + idAO + 1);
                for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                  g_xyz[idcoord] += 2.0*HijPij*Sij*muTemp[idcoord]/shellpoly;
                  g_xyz[idcoord] += gSDQ(idcoord + 1,1)(jdAO + 1,idAO + 1)*(2.0*HijPij - 2.0*PEWij - Pij*(VES(ishellA + 1,idAtm + 1) + VES(ishellB + 1,idBtm + 1)) + Pij*(VS[idAtm] + VS[idBtm]));
                }
                for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                  dtmp = 0.0;
                  qtmp = 0.0;
                  for (size_t jdcoord = 1; jdcoord < 4; ++jdcoord) {
                    dtmp += gSDQ(idcoord + 1,jdcoord + 10)(jdAO + 1,idAO + 1)*VD(idAtm + 1,jdcoord) + gSDQ(idcoord + 1,jdcoord + 1)(jdAO + 1,idAO + 1)*VD(idBtm + 1,jdcoord);
                    qtmp += gSDQ(idcoord + 1,jdcoord + 13)(jdAO + 1,idAO + 1)*VQ(idAtm + 1,jdcoord) + gSDQ(idcoord + 1,jdcoord + 16)(jdAO + 1,idAO + 1)*VQ(idAtm + 1,jdcoord + 3);
                    qtmp += gSDQ(idcoord + 1,jdcoord + 4)(jdAO + 1,idAO + 1)*VQ(idBtm + 1,jdcoord) + gSDQ(idcoord + 1,jdcoord + 7)(jdAO + 1,idAO + 1)*VQ(idBtm + 1,jdcoord + 3);
                  }
                  g_xyz[idcoord] += Pij*(dtmp + qtmp);
                }
                dCN += (HijPij/haverage)*Sij;
              }
            }
            dHdCN[idAtm] += dCN*dSEdCN(ishellA + 1,idAtm + 1);
            dHdCN[idBtm] += dCN*dSEdCN(ishellB + 1,idBtm + 1);
            for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
              grad(3*idAtm + idcoord + 1,1) += g_xyz[idcoord];
              grad(3*idBtm + idcoord + 1,1) -= g_xyz[idcoord];
            }
            naosdoneB += 2*LorbB + 1;
            ++LorbB;
            if ((transitionmetalB)&&(ishellB == 0)) {LorbB = 0;}
          }
          naosdoneA += 2*LorbA + 1;
          ++LorbA;
          if ((transitionmetalA)&&(ishellA == 0)) {LorbA = 0;}
        }
        cntB += norbsB;
      }
      naosdoneA = 0;
      LorbA = 2*int(transitionmetalA);
      for (size_t ishellA = 0; ishellA < nshellsA; ++ishellA) {
        for (size_t idAO = 0; idAO < 2*LorbA + 1; ++idAO) {
          Pij = dens(cntA + naosdoneA + idAO + 1,cntA + naosdoneA + idAO + 1);
          dHdCN[idAtm] += Pij*dSEdCN(ishellA + 1,idAtm + 1)/au2eV;
        }
        naosdoneA += 2*LorbA + 1;
        ++LorbA;
        if ((transitionmetalA)&&(ishellA == 0)) {LorbA = 0;}
      }
      cntA += norbsA;
    }
  }
  void AnalyticalGrad(matrixE & gradient, double tolerance = 1.0e-6) {
    //function calculating the gradients of the energy for GFN2 methods
    //all final derivatives come for nuclear coordinates in atomic units!
    gradient.resize(3*Natoms,1);
    double d2threshold = 1600.0;
    int maxnshells = 0;
    int counter = 1;
    int atomA;
    int atomB;
    int maxrefn;
    int szc6coeff = 0;
    int nshellsA;
    int nshellsB;
    int LorbA;
    int LorbB;
    int irow = 1;
    int icol = 1;
    bool transitionmetalA;
    bool transitionmetalB;
    double distAB;
    double yeffectiveA;
    double yeffectiveB;
    double alphaA;
    double alphaB;
    double krep;
    double STemp;
    double R0A;
    double R0B;
    double dum3a;
    double dum5a;
    double rAB2;
    double rABn;
    double f3;
    double f5;
    double a3 = alphaMultipole(1);
    double a5 = alphaMultipole(2);
    double t2a;
    double qs1;
    double qs2;
    double etaA;
    double etaB;
    std::vector<double> dHdCN(Natoms,0.0);
    matrixE geometry = mol.Geometry();
    std::vector<matrixE> dNCdR(3);
    matrixE dSEdCN(Natoms,QShells.rows());
    dSEdCN.zero();
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      dNCdR[idcoord] = dSEdCN;                            //-> for dj/dr
    }
    //block 1
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      atomA = atoms[idAtm];
      transitionmetalA = ((atomA > 20)&&(atomA < 30))||((atomA > 38)&&(atomA < 48))||((atomA > 56)&&(atomA < 80));
      nshellsA = AtomShells[idAtm];
      R0A = CNrad[idAtm];
      refA[idAtm] = RefN(atomA);
      szc6coeff += refA[idAtm];
      if (nshellsA > maxnshells) {maxnshells = nshellsA;}
      STemp = 0.0;
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        rA[idcoord] = geometry(idAtm + 1,idcoord + 1)*dist_Angstrom2aum1;
        muTemp[idcoord] = 0.0;
        thetaTemp[idcoord] = 0.0;
        thetaTemp[idcoord + 3] = 0.0;
        gradient(3*idAtm + idcoord + 1,1) = 0.0;
      }
      yeffectiveA = Yeff(atomA);
      alphaA = sqrt(AlphaGFN2(atomA));
      etaA = Hardness(atomA);
      icol = 1;
      for (size_t idBtm = 0; idBtm < Natoms; ++idBtm) {
        rAB2 = 0.0;
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          rAB[idcoord] = rA[idcoord] - geometry(idBtm + 1,idcoord + 1)*dist_Angstrom2aum1;
          rAB2 += rAB[idcoord]*rAB[idcoord];
        }
        distAB = sqrt(rAB2);
        if (idBtm < idAtm) {
          atomB = atoms[idBtm];
          transitionmetalB = ((atomB > 20)&&(atomB < 30))||((atomB > 38)&&(atomB < 48))||((atomB > 56)&&(atomB < 80));
          nshellsB = AtomShells[idBtm];
          etaB = Hardness(atomB);
          LorbA = 2*int(transitionmetalA);
          for (size_t ishellA = 0; ishellA < nshellsA; ++ishellA) {
            yeffectiveB = HardnessLScal(atomA,LorbA);
            LorbB = 2*int(transitionmetalB);
            for (size_t ishellB = 0; ishellB < nshellsB; ++ishellB) {
              alphaB = HardnessLScal(atomB,LorbB);
              t2a = 0.5*(etaA*yeffectiveB + etaB*alphaB);
              dum5a = 1.0/(t2a*t2a);
              rABn = 1.0/(rAB2 + dum5a);
              dum5a = sqrt(rABn);
              for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                dum3a = rAB[idcoord]*rABn*dum5a;
                dNCdR[idcoord](idAtm + 1,icol + ishellB) -= dum3a*QShells(irow + ishellA,1);
                dNCdR[idcoord](idBtm + 1,irow + ishellA) += dum3a*QShells(icol + ishellB,1);
              }
              ++LorbB;
              if ((transitionmetalB)&&(ishellB == 0)) {LorbB = 0;}
            }
            ++LorbA;
            if ((transitionmetalA)&&(ishellA == 0)) {LorbA = 0;}
          }
          icol += nshellsB;
          if (rAB2 <= d2threshold) {
            krep = 1.0;
            if ((atomA > 2)||(atomB > 2)) {krep += 0.5;}
            krep *= 0.5;
            distAB = sqrt(rAB2);
            yeffectiveB = Yeff(atomB);
            alphaB = sqrt(AlphaGFN2(atomB));
            rABn = pow(rAB2,krep);
            t2a = yeffectiveA*yeffectiveB*exp(-alphaA*alphaB*rABn)/distAB;
            t2a *= (1.0 + 2.0*krep*alphaA*alphaB*rABn)/rAB2;
            for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
              gradient(3*idAtm + idcoord + 1,1) -= t2a*rAB[idcoord];
              gradient(3*idBtm + idcoord + 1,1) += t2a*rAB[idcoord];
            }
          }
        }
        if (idAtm == idBtm) {continue;}
        R0B = CNrad[idBtm];
        distAB = sqrt(rAB2);
        rABn = rAB2*distAB;
        f3 = MultipoleDamp(distAB,0.5*(R0A + R0B),a3)/rABn;
        rABn *= rAB2;
        f5 = MultipoleDamp(distAB,0.5*(R0A + R0B),a5)/rABn;
        dum3a = 0.0;
        dum5a = 0.0;
        rAB2 = 0.0;
        t2a = 0.0;
        for (size_t idcoord1 = 0; idcoord1 < 3; ++idcoord1) {
          rAB2 += rAB[idcoord1]*rAB[idcoord1];
          t2a += muAtoms(idBtm + 1,idcoord1 + 1)*rAB[idcoord1];
          for (size_t idcoord2 = 0; idcoord2 < 3; ++idcoord2) {
            counter = lin(idcoord1 + 1,idcoord2 + 1);
            dum5a -= thetaAtoms(idBtm + 1,counter)*rAB[idcoord1]*rAB[idcoord2];
            if (idcoord2 >= idcoord1) {continue;}
            thetaTemp[idcoord1 + idcoord2 + 2] -= 3.0*QAtoms[idBtm]*f5*rAB[idcoord1]*rAB[idcoord2];
          }
          thetaTemp[idcoord1] -= 1.5*QAtoms[idBtm]*f5*rAB[idcoord1]*rAB[idcoord1];
        }
        dum3a = -t2a;
        STemp += dum5a*f5 + dum3a*f3;
        for (size_t idcoord1 = 0; idcoord1 < 3; ++idcoord1) {
          dum3a = rAB[idcoord1]*QAtoms[idBtm];
          dum5a = 3.0*rAB[idcoord1]*t2a - rAB2*muAtoms(idBtm + 1,idcoord1 + 1);
          muTemp[idcoord1] += dum5a*f5 + dum3a*f3;
          thetaTemp[idcoord1] += 0.5*rAB2*QAtoms[idBtm]*f5;
        }
      }
      irow += nshellsA;
      VS[idAtm] = STemp;
      for (int idcoord = 1; idcoord < 4; ++idcoord) {
        VD(idAtm + 1,idcoord) = muTemp[idcoord - 1];
        VQ(idAtm + 1,idcoord) = thetaTemp[idcoord - 1];
        VQ(idAtm + 1,idcoord + 3) = thetaTemp[idcoord + 2];
      }
      qs1 = 2.0*fXCmuA(atomA);
      qs2 = 6.0*fXCthetaA(atomA);
      t2a = 0.0;
      for (int idcoord1 = 0; idcoord1 < 3; ++idcoord1) {
        VD(idAtm + 1,idcoord1 + 1) -= qs1*muAtoms(idAtm + 1,idcoord1 + 1);
        for (size_t idcoord2 = 0; idcoord2 < idcoord1; ++idcoord2) {
          counter = lin(idcoord1 + 1,idcoord2 + 1);
          VQ(idAtm + 1,idcoord1 + idcoord2 + 3) -= thetaAtoms(idAtm + 1,counter)*qs2;
        }
        counter = lin(idcoord1 + 1,idcoord1 + 1);
        VQ(idAtm + 1,idcoord1 + 1) -= 0.5*thetaAtoms(idAtm + 1,counter)*qs2;
        t2a += thetaAtoms(idAtm + 1,counter);
      }
      t2a *= fXCthetaA(atomA);
      for (int idcoord1 = 0; idcoord1 < 3; ++idcoord1) {
        VQ(idAtm + 1,idcoord1 + 1) += t2a;
      }
    }
    dSEdCN = dNCdR[0]*QShells;
    auxMAT2 = dNCdR[1]*QShells;
    auxMAT3 = dNCdR[2]*QShells;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      gradient(3*idAtm + 1,1) += dSEdCN(idAtm + 1,1);
      gradient(3*idAtm + 2,1) += auxMAT2(idAtm + 1,1);
      gradient(3*idAtm + 3,1) += auxMAT3(idAtm + 1,1);
    }
    dSEdCN.resize(maxnshells,Natoms);
    auxMAT2 = VES;
    VES.resize(maxnshells,Natoms);
    counter = 1;
    for (size_t iAtm = 0; iAtm < Natoms; ++iAtm) {
      for (size_t ishell = 0; ishell < AtomShells[iAtm]; ++ishell, ++counter) {
        VES(ishell + 1,iAtm + 1) = auxMAT2(counter,1);
      }
    }
    boccupancy = occupancy;
    for (size_t idao = 0; idao < NAOs; ++idao) {
      boccupancy[idao] *= EMOs[idao];
    }
    if ((bdens.rows() != NAOs)&&(bdens.cols() != NAOs)) {bdens.resize(NAOs,NAOs);}
    DensW(bdens,boccupancy,CMO,1.0e-8);
    dSQH0(gradient,dens,bdens,VS,VES,VD,VQ,dSEdCN,dHdCN,geometry,d2threshold,1.0e-8);
    VES = auxMAT2;                                                            //setting VES back to what it was due to solver in geometry optimisation (if a step is skipped, then it breaks)
    //block 2
    dSEdCN.resize(Natoms,Natoms);
    dSEdCN.zero();
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      dNCdR[idcoord] = dSEdCN;
    }
    maxrefn = dCNdR(dNCdR,geometry,atoms,rA,rB,dGFNcountdx,false,d2threshold,false);
    dSEdCN = dNCdR[0]*dHdCN;
    auxMAT2 = dNCdR[1]*dHdCN;
    auxMAT3 = dNCdR[2]*dHdCN;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      gradient(3*idAtm + 1,1) += dSEdCN(idAtm + 1,1);
      gradient(3*idAtm + 2,1) += auxMAT2(idAtm + 1,1);
      gradient(3*idAtm + 3,1) += auxMAT3(idAtm + 1,1);
    }
    AnisotropicGradient(geometry,gradient,dNCdR,CoordNum);
    //block 3
    WDisp.resize(szc6coeff,szc6coeff);
    getC6D4(this->Type(),WDisp,GaussW,auxv,refA,atoms,extractC6,thresholdzero);
    NCoord(dHdCN,atoms,geometry,ERFcount,true,d2threshold,false);
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      dNCdR[idcoord].zero();
    }
    maxrefn = dCNdR(dNCdR,geometry,atoms,rA,rB,dERFcountdx,true,d2threshold,false);
    dED4GFN2(this->Type(),gradient,dNCdR,WDisp,atoms,geometry,dHdCN,QAtoms,maxrefn,true,3600.0,1600.0);
    //block 4
    if (usesolvation) {
      if (CoulombSolv == 1) {calcBornMatrix_dX(gradient,dNCdR,dSEdCN,auxMAT2,auxMAT3,geometry,GBEffectiveCoulombP16_dX,1.0e-10);}
      else if (CoulombSolv == 0) {calcBornMatrix_dX(gradient,dNCdR,dSEdCN,auxMAT2,auxMAT3,geometry,GBEffectiveCoulombStill_dX,1.0e-10);}
    }
    gradient *= dist_Angstrom2aum1;
  }
  int AvailableGradients() {return 1;}
  void calcZetas(double zerocharge = 1.0) {
    //function that calculates the zetas to be used in the calculation of the dispersion Fockian and dispersion energy
    double ga = GA();
    double gc = GC();
    double zA;
    double hardness;
    double QA;
    int atomA;
    int ipos = 0;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {                        //get zeta and dzeta/dq
      atomA = atoms[iatm];
      zA = Zeff(atomA);
      hardness = HardnessD4(atomA);
      RefQ(atomA,auxv,1);
      QA = zerocharge*QAtoms[iatm];
      for (size_t irefA = 0; irefA < refA[iatm]; ++irefA,++ipos) {
        zetaA[ipos] = 0.0;
        dzetaA[ipos] = 0.0;
        if (GaussW[ipos] < gwthresh) {continue;}
        zetaA[ipos] = ZetaD4(ga,hardness*gc,auxv[irefA] + zA,QA + zA);
        dzetaA[ipos] = dZetaD4dq(ga,hardness*gc,auxv[irefA] + zA,QA + zA);
      }
    }
  }
  void calcFock(matrixE & SAO) {
    //function that calculates the TB-Fock matrix
    int atomA;
    int atomB;
    int nshellsA;
    int nshellsB;
    int LorbA;
    int LorbB;
    int posA = 1;
    int posB = 1;
    int irow = 1;
    int icol = 1;
    bool transitionmetalA;
    bool transitionmetalB;
    double GammaA;
    double GammaB;
    double FpAl;
    double FpBl;
    double qAl;
    double qBl;
    double aux1;
    double aux2;
    double aux3;
    double solvA;
    double solvB;
    int ipos = 0;
    Fock = Hcore;
    if (usesolvation) {calcSolvationForces(solvForce,QAtoms);}
    auxMAT3 = VAB*QShells;                                     //getting the terms in the summation in the first line of equation 40
    for (int iAtm = 0; iAtm < Natoms; ++iAtm) {
      //get the dA coefficients; note that it is not worth to store zeta in matrix type to make the matrix-vector product and have an additional matrix-object; by merging I reduce the number of loops
      atomA = atoms[iAtm];
      dAtoms[iAtm] = 0.0;
      for (size_t irefA = 0; irefA < refA[iAtm]; ++irefA,++ipos) {              //equation 46
        GammaA = 0.0;                                                             //get the element WDisp(i,j)zetaA(j)
        for (size_t idcontr = 0; idcontr < dispsz; ++idcontr) {
          GammaA += WDisp(ipos + 1,idcontr + 1)*zetaA[idcontr];
        }
        dAtoms[iAtm] += GammaA*dzetaA[ipos];
      }
      transitionmetalA = ((atomA > 20)&&(atomA < 30))||((atomA > 38)&&(atomA < 48))||((atomA > 56)&&(atomA < 80));
      nshellsA = AtomShells[iAtm];
      GammaA = HubbardDer(atomA);
      solvA = solvForce[iAtm];
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        rA[idcoord] = VD(iAtm + 1,idcoord + 1);
        VQA[idcoord] = VQ(iAtm + 1,idcoord + 1);
        VQA[idcoord + 3] = VQ(iAtm + 1,idcoord + 4);
      }
      LorbA = 2*int(transitionmetalA);
      for (int ishellA = 0; ishellA < nshellsA; ++ishellA,++posA) {
        qAl = QShells(posA,1);
        FpAl = dAtoms[iAtm] + auxMAT3(posA,1) + qAl*qAl*GammaA*Korb(LorbA) - VS[iAtm] + solvA;
        VES(posA,1) = FpAl + VS[iAtm];
        posB = 1;
        icol = 1;
        for (int iBtm = 0; iBtm < iAtm; ++iBtm) {
          atomB = atoms[iBtm];
          transitionmetalB = ((atomB > 20)&&(atomB < 30))||((atomB > 38)&&(atomB < 48))||((atomB > 56)&&(atomB < 80));
          nshellsB = AtomShells[iBtm];
          GammaB = HubbardDer(atomB);
          solvB = solvForce[iBtm];
          for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
            rB[idcoord] = VD(iBtm + 1,idcoord + 1);
            VQB[idcoord] = VQ(iBtm + 1,idcoord + 1);
            VQB[idcoord + 3] = VQ(iBtm + 1,idcoord + 4);
          }
          LorbB = 2*int(transitionmetalB);
          for (int ishellB = 0; ishellB < nshellsB; ++ishellB,++posB) {
            qBl = QShells(posB,1);
            FpBl = dAtoms[iBtm] + auxMAT3(posB,1) + qBl*qBl*GammaB*Korb(LorbB) - VS[iBtm] + solvB;
            for (size_t iAorb = 0; iAorb < 2*LorbA + 1; ++iAorb) {
              for (size_t iBorb = 0; iBorb < 2*LorbB + 1; ++iBorb) {
                aux1 = 0.5*SAO(irow + iAorb,icol + iBorb)*(FpAl + FpBl);                                                                             //this builds equations 40, 41a and 45
                //get the dipole (rAB) and quadrupole (thetaTemp) moment terms
                basis.getDipQuadrEl(rAB,thetaTemp,irow + iAorb,icol + iBorb);
                aux2 = 0.0;
                aux3 = 0.0;
                for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                  aux2 += 0.5*rAB[idcoord]*(rA[idcoord] + rB[idcoord]);                                                                              //equation 41b
                  aux3 += 0.5*thetaTemp[idcoord]*(VQA[idcoord] + VQB[idcoord]) + 0.5*thetaTemp[idcoord + 3]*(VQA[idcoord + 3] + VQB[idcoord + 3]);   //equation 41c
                }
                Fock(irow + iAorb,icol + iBorb) += aux2 + aux3 - aux1;
                Fock(icol + iBorb,irow + iAorb) += aux2 + aux3 - aux1;
              }
            }
            icol += 2*LorbB + 1;
            ++LorbB;
            if ((transitionmetalB)&&(ishellB == 0)) {LorbB = 0;}
          }
        }
        LorbB = 2*int(transitionmetalA);
        for (int ishellB = 0; ishellB < ishellA + 1; ++ishellB,++posB) {
          qBl = QShells(posB,1);
          FpBl = dAtoms[iAtm] + auxMAT3(posB,1) + qAl*qAl*GammaA*Korb(LorbB) - VS[iAtm] + solvA;
          for (size_t iAorb = 0; iAorb < 2*LorbA + 1; ++iAorb) {
            for (size_t iBorb = 0; iBorb < 2*LorbB + 1; ++iBorb) {
              aux1 = 0.5*SAO(irow + iAorb,icol + iBorb)*(FpAl + FpBl);                               //this builds equations 40, 41a and 45
              //get the dipole (rAB) and quadrupole (thetaTemp) moment terms
              basis.getDipQuadrEl(rAB,thetaTemp,irow + iAorb,icol + iBorb);
              aux2 = 0.0;
              aux3 = 0.0;
              for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                aux2 += rAB[idcoord]*rA[idcoord];                                                    //equation 41b
                aux3 += thetaTemp[idcoord]*VQA[idcoord] + thetaTemp[idcoord + 3]*VQA[idcoord + 3];   //equation 41c
              }
              Fock(irow + iAorb,icol + iBorb) += aux2 + aux3 - aux1;
              if (ishellB != ishellA) {Fock(icol + iBorb,irow + iAorb) += aux2 + aux3 - aux1;}
            }
          }
          icol += 2*LorbB + 1;
          ++LorbB;
          if ((transitionmetalA)&&(ishellB == 0)) {LorbB = 0;}
        }
        irow += 2*LorbA + 1;
        ++LorbA;
        if ((transitionmetalA)&&(ishellA == 0)) {LorbA = 0;}
      }
    }
    if (print == 2) {
      std::cout << "Fock matrix" << std::endl;
      Fock.Print();
    }
  }
  void transformFock(matrixE & XAO) {
    //function that transforms the Fock matrix to the orthogonal basis
    bdens = Fock*XAO;
    Fock = XAO*bdens;
  }
  void transformCMO(matrixE & XAO) {
    //function that transforms the MO coefficients back to the non-orthogonal basis
    bdens = CMO;
    CMO = XAO*bdens;
  }
  //auxiliary functions
  int Index(int aa, int bb) {
    int idx = (aa*(aa + 1) + bb*(bb + 1))/2;
    return std::min(idx,5);
  }
  double MultipoleDamp(double RAB, double R0AB, double alpha) {
    //function that calculates the damping for multipole terms
    double damp = pow(R0AB/RAB,alpha);
    return 1.0/(1.0 + 6.0*damp);                                                               //part of equation 28
  }
  void CNradius() {
    //function that calculates the coordination number dependent atomic radius
    CNrad.resize(Natoms);
    double R0A;
    double shift;
    int atmnr;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      atmnr = atoms[idAtm];
      R0A = MultipoleRadius(atmnr);
      shift = 5.0 - R0A;
      if (shift != 0.0) {shift /= (1.0 + exp(-4.0*(CoordNum[idAtm] - Nval(atmnr) - 1.2)));}
      CNrad[idAtm] = R0A + shift;                                                              //equation 29
    }
  }
  void dZero(double rcov, double rABi, double dex, double & damp, double & ddamp) {
    damp = 1.0/(1.0 + 6.0*pow(rcov*rABi,dex));
    ddamp = -dex*rABi*damp*(damp - 1.0);
  }
  double dgAB(double exp, double rABi, double damp, double ddamp) {
    //function returning the derivative of damped Coulomb interaction with respect to rAB
    double paux = pow(rABi,exp);
    double tmp = -exp*rABi*paux;
    return tmp*damp + ddamp*paux;
  }
  void AnisotropicGradient(matrixE & geometry, matrixE & gradient, std::vector<matrixE> & dCN, std::vector<double> & CoordNum) {
    //gradient contributions from the cumulative atomic multipole moments (CAMM)
    //note that interactions are included up to quadrupole moment
    double qA;
    double qB;
    double auxiliary;
    double rAB2;
    double rABi;
    double RAB;
    double damp3;
    double damp5;
    double ddamp;
    double dgAB3;
    double dgAB5;
    double gAB3;
    double gAB5;
    double ed;
    double edd;
    double eq;
    double tt;
    double a3 = alphaMultipole(1);
    double a5 = alphaMultipole(2);
    double ddm2;
    double ddm3a;
    double ddm3b;
    double qqa;
    double qqb;
    double rCNAB;
    double t2;
    int ll;
    int atomA;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      atomA = atoms[idAtm];
      rAB2 = MultipoleRadius(atomA);
      rCNAB = exp(-4.0*(CoordNum[idAtm] - Nval(atomA) - 1.2));
      t2 = (5.0 - rAB2)/(1.0 + 2.0*rCNAB + rCNAB*rCNAB);
      t2 *= 4.0*rCNAB;
      qA = QAtoms[idAtm];
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        rA[idcoord] = geometry(idAtm + 1,idcoord + 1)*dist_Angstrom2aum1;
        muA[idcoord] = muAtoms(idAtm + 1,idcoord + 1);
        VQA[idcoord] = thetaAtoms(idAtm + 1,idcoord + 1);
        VQA[idcoord + 3] = thetaAtoms(idAtm + 1,idcoord + 4);
      }
      auxiliary = 0.0;
      for (size_t idBtm = 0; idBtm < Natoms; ++idBtm) {
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          dCN[idcoord](idBtm + 1,idAtm + 1) *= t2;
        }
        if (idAtm == idBtm) {continue;}
        rAB2 = 0.0;
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          rB[idcoord] = geometry(idBtm + 1,idcoord + 1)*dist_Angstrom2aum1;
          rAB[idcoord] = rB[idcoord] - rA[idcoord];
          rAB2 += rAB[idcoord]*rAB[idcoord];
        }
        rCNAB = 0.5*(CNrad[idAtm] + CNrad[idBtm]);
        RAB = sqrt(rAB2);
        rABi = 1.0/RAB;
        dZero(rCNAB,rABi,a3,damp3,ddamp);
        dgAB3 = dgAB(3.0,rABi,damp3,ddamp);
        dZero(rCNAB,rABi,a5,damp5,ddamp);
        dgAB5 = dgAB(5.0,rABi,damp5,ddamp);
        tt = rABi*rABi*rABi;
        gAB3 = MultipoleDamp(RAB,rCNAB,a3)*tt;
        tt *= rABi*rABi;
        gAB5 = MultipoleDamp(RAB,rCNAB,a5)*tt;
        ed = 0.0;
        edd = 0.0;
        eq = 0.0;
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          ed += QAtoms[idBtm]*muA[idcoord]*rAB[idcoord] - qA*muAtoms(idBtm + 1,idcoord + 1)*rAB[idcoord];
          tt = qA*muAtoms(idBtm + 1,idcoord + 1) - QAtoms[idBtm]*muA[idcoord];
          gradient(3*idAtm + idcoord + 1,1) += tt*gAB3;
          ddm2 = 0.0;
          ddm3a = 0.0;
          ddm3b = 0.0;
          qqa = 0.0;
          qqb = 0.0;
          for (size_t jdcoord = 0; jdcoord < 3; ++jdcoord) {
            ll = lin(idcoord + 1,jdcoord + 1);
            tt = rAB[idcoord]*rAB[jdcoord];
            eq += QAtoms[idBtm]*VQA[ll - 1]*tt + qA*thetaAtoms(idBtm + 1,ll)*tt;
            edd -= 3.0*tt*muA[jdcoord]*muAtoms(idBtm + 1,idcoord + 1);
            ddm2 += muAtoms(idBtm + 1,jdcoord + 1)*muA[jdcoord];
            ddm3a += muA[jdcoord]*rAB[jdcoord];
            ddm3b += muAtoms(idBtm + 1,jdcoord + 1)*rAB[jdcoord];
            qqa += rAB[jdcoord]*thetaAtoms(idBtm + 1,ll);
            qqb += rAB[jdcoord]*VQA[ll - 1];
          }
          edd += muAtoms(idBtm + 1,idcoord + 1)*muA[idcoord]*rAB2;
          gradient(3*idAtm + idcoord + 1,1) -= 2.0*gAB5*ddm2*rAB[idcoord];
          gradient(3*idAtm + idcoord + 1,1) += 3.0*gAB5*ddm3a*muAtoms(idBtm + 1,idcoord + 1);
          gradient(3*idAtm + idcoord + 1,1) += 3.0*gAB5*ddm3b*muA[idcoord];
          gradient(3*idAtm + idcoord + 1,1) -= 2.0*gAB5*qqa*qA;
          gradient(3*idAtm + idcoord + 1,1) -= 2.0*gAB5*qqb*QAtoms[idBtm];
        }
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          tt = rAB[idcoord]*rABi;
          gradient(3*idAtm + idcoord + 1,1) -= ed*dgAB3*tt + (eq + edd)*dgAB5*tt;
        }
        tt = ed*a3*gAB3*(damp3/rCNAB)*pow(rCNAB*rABi,a3);
        auxiliary += tt;
        tt = (eq + edd)*a5*gAB5*(damp5/rCNAB)*pow(rCNAB*rABi,a5);
        auxiliary += tt;
      }
      auxiliary *= 3.0;
      for (size_t idBtm = 0; idBtm < Natoms; ++idBtm) {
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          gradient(3*idBtm + idcoord + 1,1) -= auxiliary*dCN[idcoord](idBtm + 1,idAtm + 1);
        }
      }
    }
  }
  void Broyden(int iter) {
    //this function applies some sort of DIIS/damping to help improve convergence of the charge model
    int it1 = iter - 1;
    //get dq
    size_t idshell = 0;
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      for (size_t ishellA = 0; ishellA < AtomShells[idatm]; ++ishellA,++idshell) {
        dqvec[idshell] = QShells(idshell + 1,1) - qvec[idshell];
      }
      for (size_t idpos = 0; idpos < 3; ++idpos) {
        dqvec[9*idatm + totalshells + idpos] = muAtoms(idatm + 1,idpos + 1) - qvec[9*idatm + totalshells + idpos];
        dqvec[9*idatm + totalshells + idpos + 3] = thetaAtoms(idatm + 1,idpos + 1) - qvec[9*idatm + totalshells + idpos + 3];
        dqvec[9*idatm + totalshells + idpos + 6] = thetaAtoms(idatm + 1,idpos + 4) - qvec[9*idatm + totalshells + idpos + 6];
      }
    }
    if (it1 == 0) {       //first iteration, simple damp
      dqvec_last = dqvec;
      qvec_last = qvec;
      for (size_t idvar = 0; idvar < totalshells + 9*Natoms; ++idvar) {
        qvec[idvar] += Broydendamp*dqvec[idvar];
      }
    }
    else {
      double aux = 0.0;
      double aux2 = 0.0;
      double dfdf;
      bmatrix.resize(it1,it1);
      gvec.resize(it1);
      //get weight for current iteration (omega) and dF
      for (size_t idx = 0; idx < totalshells + 9*Natoms; ++idx) {
        aux += dqvec[idx]*dqvec[idx];
        dfmat(it1,idx + 1) = dqvec[idx] - dqvec_last[idx];
        aux2 += dfmat(it1,idx + 1)*dfmat(it1,idx + 1);
      }
      omega[it1] = sqrt(aux);
      if (omega[it1] > (wfac/maxw)) {omega[it1] = wfac/sqrt(aux);}
      else {omega[it1] = maxw;}
      if (omega[it1] < minw) {omega[it1] = minw;}
      dfdf = 1.0/sqrt(aux2);
      for (size_t idx = 0; idx < totalshells + 9*Natoms; ++idx) {
        dfmat(it1,idx + 1) *= dfdf;
      }
      for (size_t idx = 0; idx < it1; ++idx) {
        aux = 0.0;
        aux2 = 0.0;
        for (size_t idy = 0; idy < totalshells + 9*Natoms; ++idy) {
          aux += dfmat(idx + 1,idy + 1)*dfmat(it1,idy + 1);
          aux2 += dfmat(idx + 1,idy + 1)*dqvec[idy];
        }
        amat(idx + 1,it1) = aux;
        amat(it1,idx + 1) = aux;
        gvec[idx] = omega[idx + 1]*aux2;
      }
      for (size_t idx = 1; idx < it1 + 1; ++idx) {
        for (size_t idy = 1; idy < it1 + 1; ++idy) {
          bmatrix(idy,idx) = omega[idy]*amat(idy,idx)*omega[idx];
        }
        bmatrix(idx,idx) += omega0*omega0;
      }
      Solve_Ax_eq_b(bmatrix,gvec);
      //build |u>
      for (size_t idx = 0; idx < totalshells + 9*Natoms; ++idx) {
        umat(it1,idx + 1) = Broydendamp*dfmat(it1,idx + 1) + dfdf*(qvec[idx] - qvec_last[idx]);
        dqvec_last[idx] = dqvec[idx];
        qvec_last[idx] = qvec[idx];
        //get the new charges
        qvec[idx] += Broydendamp*dqvec[idx];
      }
      for (size_t idx = 0; idx < it1; ++idx) {
        for (size_t idy = 0; idy < totalshells + 9*Natoms; ++idy) {
          qvec[idy] -= omega[idx + 1]*gvec[idx]*umat(idx + 1,idy + 1);
        }
      }
    }
    //move back
    idshell = 0;
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      QAtoms[idatm] = 0.0;
      for (size_t ishellA = 0; ishellA < AtomShells[idatm]; ++ishellA,++idshell) {
        QShells(idshell + 1,1) = qvec[idshell];
        QAtoms[idatm] += QShells(idshell + 1,1);
      }
      for (size_t idpos = 0; idpos < 3; ++idpos) {
        muAtoms(idatm + 1,idpos + 1) = qvec[9*idatm + totalshells + idpos];
        thetaAtoms(idatm + 1,idpos + 1) = qvec[9*idatm + totalshells + idpos + 3];
        thetaAtoms(idatm + 1,idpos + 4) = qvec[9*idatm + totalshells + idpos + 6];
      }
    }
  }
  //other auxiliary functions
  void checkAtoms() {
    //function that checks whether the atoms in the molecule are fully defined in GFN2
    bool def;
    for (size_t iatm = 0; iatm < atoms.size(); ++iatm) {
      def = (atoms[iatm] < 87);
      if (!def) {throw("ERROR: GFN.hpp: GFN2: checkAtoms(): atom not fully specified in GFN2-xTB");}
    }
  }
  void OpenClosed() {
    //function to determine whether system is open- or closed-shell and the number of core-electrons
    if (nopen < 2) {                                       //in this case we must determine as we are working with the defaults
      shell = "closed";
      nopen = 0;
      if (mol.Nelectrons()%2 != 0) {
        shell = "open";
        nopen = mol.Nelectrons()%2;
      }
      else if (mol.Multiplicity() > 1) {
        shell = "open";
        nopen = mol.Multiplicity() - 1;
      }
      openclosed = shell;
    }
    else {shell = openclosed;}                             //nopen was defined by the user
    //get the number of core electrons (ncore function in file molecule.f90)
    if (CoreCharge.size() != Natoms) {CoreCharge.resize(Natoms);}
    coreelectrons = 0;
    for (size_t idx = 0; idx < Natoms; ++idx) {
      if (atoms[idx] < 3) {CoreCharge[idx] = 0;}                                            //H,He
      else if ((atoms[idx] > 2)&&(atoms[idx] < 11))  {CoreCharge[idx] = 2;}                 //Li-Ne
      else if ((atoms[idx] > 10)&&(atoms[idx] < 19)) {CoreCharge[idx] = 10;}                //Na-Ar
      else if ((atoms[idx] > 18)&&(atoms[idx] < 30)) {CoreCharge[idx] = 18;}                //K-Cu
      else if ((atoms[idx] > 29)&&(atoms[idx] < 37)) {CoreCharge[idx] = 28;}                //Zn-Kr
      else if ((atoms[idx] > 36)&&(atoms[idx] < 48)) {CoreCharge[idx] = 36;}                //Rb-Ag
      else if ((atoms[idx] > 47)&&(atoms[idx] < 55)) {CoreCharge[idx] = 46;}                //Cd-Xe
      else if ((atoms[idx] > 54)&&(atoms[idx] < 58)) {CoreCharge[idx] = 54;}                //Cs-La
      else if ((atoms[idx] > 57)&&(atoms[idx] < 72)) {CoreCharge[idx] = atoms[idx] - 3;}    //Ce-Lu
      else if ((atoms[idx] > 71)&&(atoms[idx] < 80)) {CoreCharge[idx] = 68;}                //Hf-Au
      else if ((atoms[idx] > 79)&&(atoms[idx] < 87)) {CoreCharge[idx] = 78;}                //Hg-Bi
      coreelectrons += CoreCharge[idx];
    }
    if (coreelectrons%2 != 0) {
      if (nopen == 0) {
        nopen += 1;
        coreelectrons -= 1;
      }
      else {
        nopen -= 1;
        Nelectrons += 1;
      }
    }
    if (print == 2) {
      std::cout << "core charge" << std::endl;
      for (size_t idx = 0; idx < CoreCharge.size(); ++idx) {
        std::cout << CoreCharge[idx] << " ";
      }
      std::cout << std::endl;
    }
  }
  virtual double ZeroOverlap(size_t atm) {return ZeroOverlapGFN2(atm);}
  //solvation model stuff
  double DeScreening(std::string solvent, size_t atmnr) {
    //descreening constant used to calculate Born-radii
    return DeScreeningGFN2(solvent,atmnr);
  }
  double FreeEnergyShift(std::string solvent_) {
    //free energy shift for ALPB model
    return FreeEnergyShiftGFN2(solvent_);
  }
  double BornRadiusScale(std::string solvent_) {
    //scaling parameter for Born radii
    return BornRadiusScaleGFN2(solvent_);
  }
  double ProbeRadius(std::string solvent_) {
    //probe radii
    return ProbeRadiusGFN2(solvent_);
  }
  double GammaScale(std::string solvent, size_t atmnr) {
    //function returning the gamma's for surface tension
    return GammaScaleGFN2(solvent,atmnr);
  }
  double HbondStrength(std::string solvent, size_t atmnr) {
    //function returning the hydrogen-bonding strength of each atom
    return HbondStrengthGFN2(solvent,atmnr);
  }
  //dispersion related
  void AtomicPolarizabilities(std::vector<double> & polarizabilities, const std::vector<double> & AtmCharge) {
    //function that returns a vector with the atom-specific polarizabilities
    matrixE geom = mol.Geometry();
    Polarizabilities(polarizabilities,atoms,geom,AtmCharge,1);
  }
  void TotalDispersion(double & mc6, double & mc8, const std::vector<double> & AtmCharge) {
    //function to calculate the molecular C6 and C8 coefficients
    double mpol;
    matrixE geom = mol.Geometry();
    MolecularDispersionData(mc6,mc8,mpol,atoms,geom,AtmCharge,1);
  }
  void TotalPolarizability(double & mpol, const std::vector<double> & AtmCharge) {
    //function to get the molecular polarizability
    double mc6;
    double mc8;
    matrixE geom = mol.Geometry();
    MolecularDispersionData(mc6,mc8,mpol,atoms,geom,AtmCharge,1);
  }
};

#endif //_GFN_METHODS_
