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

#ifndef _Basis_Set_
#define _Basis_Set_
#include <string>
#include "math/MatrixPackage.hpp"
#include "math/TensorPackage.hpp"
#include "math/IntegralPackage.hpp"
#include "Molecule.hpp"
#include "UnitConversion.hpp"
#include "basissets/burns.hpp"
#include "basissets/dewar.hpp"
#include "basissets/gfn2.hpp"
#include "basissets/pople.hpp"
#include "basissets/slater.hpp"

//description:
//main driver for basis sets

class BSet {
  size_t NAOs;                                      //total number of AOs
  int Natoms;
  int iptype;
  matrixE sao;                                      //overlap matrix
  matrixE xao;                                      //orthogonalization matrix for the SAO matrix
  matrixE geometry;
  matrixE overlaps;
  matrixE transfo;
  matrixE saux;
  matrixE s12;
  matrixE val;
  matrixE grad;
  matrixE gints;
  matrix<size_t> QNA;
  matrix<size_t> QNB;
  Molecule mol;                                     //molecule object with all its information
  std::string bset;                                 //basis set name to find the exponents
  std::vector<std::vector<double> > exponents;      //vector of exponents for each atom
  tensor<matrixE,double> dip;                       //dipole moment integrals, stored as x,y,z
  tensor<matrixE,double> quadr;                     //quadrupole moment integrals, stored as xx,yy,zz,xy,xz,yz
  std::vector<size_t> atoms;
  std::vector<double> dp;
  std::vector<double> dp2;
  std::vector<double> qp;
  std::vector<double> rA;
  std::vector<double> rB;
  std::vector<double> rP;
  std::vector<double> diffAtom;
  std::vector<double> CoeffsA;
  std::vector<double> CoeffsB;
  std::vector<double> alpha;
  std::vector<double> beta;
  std::vector<double> tk;
  std::vector<double> ints;
  std::vector<double> vi;
  std::vector<double> vj;
  std::vector<double> vv;
  std::vector<double> gi;
  std::vector<double> gg;
public:
  BSet(): mol() {
    resizeVectors();
    saux.resize(6,6);
    resizeGTOArrays();
  }
  BSet(std::string _bset): mol() {
    resizeVectors();
    bset = _bset;
    saux.resize(6,6);
    if ((bset == "gfn2")||(bset == "GFN2")) {resizeGTOArrays();}
    if ((bset == "gfn2l")||(bset == "GFN2L")) {resizeGTOArrays();}
  }
  BSet(Molecule _mol, std::string _bset, bool rotate = true): mol(_mol) {
    resizeVectors();
    bset = _bset;
    sao = Identity(1);
    xao = Identity(1);
    atoms = mol.Atoms();
    Natoms = atoms.size();
    geometry = mol.Geometry();
    saux.resize(6,6);
    bool orthogonalize = true;
    iptype = 0;
    if (bset == "slater") {
      SlaterBSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
      iptype = bst.IPType();
    }
    else if (bset == "pople") {
      PopleBSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
      iptype = bst.IPType();
    }
    else if ((bset == "pm6")||(bset == "PM6")) {
      PM6BSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
      orthogonalize = false;
      iptype = bst.IPType();
    }
    else if ((bset == "pm3")||(bset == "PM3")) {
      PM3BSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
      orthogonalize = false;
      iptype = bst.IPType();
    }
    else if ((bset == "pm3pddg")||(bset == "PM3PDDG")) {
      PM3PDDGBSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
      orthogonalize = false;
      iptype = bst.IPType();
    }
    else if ((bset == "mndopddg")||(bset == "MNDOPDDG")) {
      MNDOPDDGBSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
      orthogonalize = false;
      iptype = bst.IPType();
    }
    else if ((bset == "pm3bp")||(bset == "PM3BP")) {
      PM3BPBSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
      orthogonalize = false;
      iptype = bst.IPType();
    }
    else if ((bset == "rm1")||(bset == "RM1")) {
      RM1BSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
      orthogonalize = false;
      iptype = bst.IPType();
    }
    else if ((bset == "am1")||(bset == "AM1")) {
      AM1BSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
      orthogonalize = false;
      iptype = bst.IPType();
    }
    else if ((bset == "mndo")||(bset == "MNDO")) {
      MNDOBSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
      orthogonalize = false;
      iptype = bst.IPType();
    }
    else if ((bset == "mndod")||(bset == "MNDOd")) {
      MNDOdBSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
      orthogonalize = false;
      iptype = bst.IPType();
    }
    else if (bset == "burns") {
      BurnsBSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
      iptype = bst.IPType();
    }
    else if ((bset == "gfn2")||(bset == "GFN2")) {
      resizeGTOArrays();
      GFN2BSet bst;
      BasisSetData(bst);
      OvDipQuadrGTO(bst);
      iptype = bst.IPType();
    }
    else if ((bset == "gfn2l")||(bset == "GFN2L")) {
      resizeGTOArrays();
      GFN2BSet bst;
      BasisSetData(bst);
      OverlapGTO(bst);
      iptype = bst.IPType();
    }
    if ((rotate)&&(orthogonalize)) {OrthOverlap();}
  }
  ~BSet() {}
  //getters
  std::vector<std::vector<double> > GetExponents() {return exponents;}
  std::vector<double> GetExponents(size_t atom) {return exponents[atom - 1];}
  double elementSAO(size_t row, size_t col) {return sao(row,col);}
  matrixE SAO() {return sao;}
  tensor<matrixE,double> Dipole() {return dip;}
  tensor<matrixE,double> Quadrupole() {return quadr;}
  void getDipole(tensor<matrixE,double> & _dip) {_dip = dip;}
  void getQuadrupole(tensor<matrixE,double> & _quadr) {_quadr = quadr;}
  matrixE XAO() {return xao;}
  size_t NAO() {return NAOs;}
  int IonPotType() {return iptype;}
  //setters
  void setSAO(matrixE & newsao) {sao = newsao;}
  void setAOBasis(const std::string & _bset) {bset = _bset;}
  void setMolecule(Molecule & _mol, bool rotate = true) {
    //function receiving new molecule and restaring the whole calculation
    mol = _mol;
    atoms = mol.Atoms();
    Natoms = atoms.size();
    geometry = mol.Geometry();
    exponents.clear();
    bool orthogonalize = true;
    if (bset == "slater") {
      SlaterBSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
    }
    else if (bset == "pople") {
      PopleBSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
    }
    else if ((bset == "pm6")||(bset == "PM6")) {
      PM6BSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
      orthogonalize = false;
    }
    else if ((bset == "pm3")||(bset == "PM3")) {
      PM3BSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
      orthogonalize = false;
    }
    else if ((bset == "pm3pddg")||(bset == "PM3PDDG")) {
      PM3PDDGBSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
      orthogonalize = false;
    }
    else if ((bset == "mndopddg")||(bset == "MNDOPDDG")) {
      MNDOPDDGBSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
      orthogonalize = false;
    }
    else if ((bset == "pm3bp")||(bset == "PM3BP")) {
      PM3BPBSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
      orthogonalize = false;
    }
    else if ((bset == "rm1")||(bset == "RM1")) {
      RM1BSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
      orthogonalize = false;
    }
    else if ((bset == "am1")||(bset == "AM1")) {
      AM1BSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
      orthogonalize = false;
    }
    else if ((bset == "mndo")||(bset == "MNDO")) {
      MNDOBSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
      orthogonalize = false;
    }
    else if ((bset == "mndod")||(bset == "MNDOd")) {
      MNDOdBSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
      orthogonalize = false;
    }
    else if (bset == "burns") {
      BurnsBSet bst;
      BasisSetData(bst);
      Overlap(bst,rotate);
    }
    else if ((bset == "gfn2")||(bset == "GFN2")) {
      GFN2BSet bst;
      BasisSetData(bst);
      OvDipQuadrGTO(bst);
    }
    else if ((bset == "gfn2l")||(bset == "GFN2L")) {
      GFN2BSet bst;
      BasisSetData(bst);
      OverlapGTO(bst);
    }
    if ((rotate)&&(orthogonalize)) {OrthOverlap();}
  }
  void resizeVectors() {
    dp.resize(3);
    dp2.resize(3);
    qp.resize(6);
    rA.resize(3);
    rB.resize(3);
    diffAtom.resize(3);
  }
  void resizeGTOArrays() {
    rP.resize(3);
    tk.resize(9);
    ints.resize(10);
    vi.resize(7);
    vj.resize(7);
    vv.resize(13);
    gi.resize(7);
    gg.resize(13);
    val.resize(3,3);
    grad.resize(3,3);
    gints.resize(3,10);
  }
  template<class T>
  void BasisSetData(T & pbst) {
    //function to treat basis set information: calculate NAOs and exponents
    NAOs = 0;
    size_t atomtp;
    for (size_t idx = 0; idx < Natoms; ++idx) {
      atomtp = atoms[idx];
      NAOs += pbst.GetNAOs(atomtp);
      exponents.push_back(pbst.GetExponents(atomtp));
    }
  }
  std::vector<size_t> AtomNAOs(std::vector<size_t> atomnum) {
    //function that returns the number of AOs for a given set of atoms
    //atomnum is here a vector with the atomic numbers
    size_t Natoms = atomnum.size();
    std::vector<size_t> atomnr(Natoms,0);
    if (bset == "slater") {
      SlaterBSet bst;
      for (size_t idx = 0; idx < Natoms; ++idx) {
        atomnr[idx] = bst.GetNAOs(atomnum[idx]);
      }
    }
    else if (bset == "pople") {
      PopleBSet bst;
      for (size_t idx = 0; idx < Natoms; ++idx) {
        atomnr[idx] = bst.GetNAOs(atomnum[idx]);
      }
    }
    else if ((bset == "pm6")||(bset == "PM6")) {
      PM6BSet bst;
      for (size_t idx = 0; idx < Natoms; ++idx) {
        atomnr[idx] = bst.GetNAOs(atomnum[idx]);
      }
    }
    else if ((bset == "pm3")||(bset == "PM3")) {
      PM3BSet bst;
      for (size_t idx = 0; idx < Natoms; ++idx) {
        atomnr[idx] = bst.GetNAOs(atomnum[idx]);
      }
    }
    else if ((bset == "pm3pddg")||(bset == "PM3PDDG")) {
      PM3PDDGBSet bst;
      for (size_t idx = 0; idx < Natoms; ++idx) {
        atomnr[idx] = bst.GetNAOs(atomnum[idx]);
      }
    }
    else if ((bset == "mndopddg")||(bset == "MNDOPDDG")) {
      MNDOPDDGBSet bst;
      for (size_t idx = 0; idx < Natoms; ++idx) {
        atomnr[idx] = bst.GetNAOs(atomnum[idx]);
      }
    }
    else if ((bset == "pm3bp")||(bset == "PM3BP")) {
      PM3BPBSet bst;
      for (size_t idx = 0; idx < Natoms; ++idx) {
        atomnr[idx] = bst.GetNAOs(atomnum[idx]);
      }
    }
    else if ((bset == "rm1")||(bset == "RM1")) {
      RM1BSet bst;
      for (size_t idx = 0; idx < Natoms; ++idx) {
        atomnr[idx] = bst.GetNAOs(atomnum[idx]);
      }
    }
    else if ((bset == "am1")||(bset == "AM1")) {
      AM1BSet bst;
      for (size_t idx = 0; idx < Natoms; ++idx) {
        atomnr[idx] = bst.GetNAOs(atomnum[idx]);
      }
    }
    else if ((bset == "mndo")||(bset == "MNDO")) {
      MNDOBSet bst;
      for (size_t idx = 0; idx < Natoms; ++idx) {
        atomnr[idx] = bst.GetNAOs(atomnum[idx]);
      }
    }
    else if ((bset == "mndod")||(bset == "MNDOd")) {
      MNDOdBSet bst;
      for (size_t idx = 0; idx < Natoms; ++idx) {
        atomnr[idx] = bst.GetNAOs(atomnum[idx]);
      }
    }
    else if (bset == "burns") {
      BurnsBSet bst;
      for (size_t idx = 0; idx < Natoms; ++idx) {
        atomnr[idx] = bst.GetNAOs(atomnum[idx]);
      }
    }
    else if ((bset == "gfn2")||(bset == "GFN2")) {
      GFN2BSet bst;
      for (size_t idx = 0; idx < Natoms; ++idx) {
        atomnr[idx] = bst.GetNAOs(atomnum[idx]);
      }
    }
    else if ((bset == "gfn2l")||(bset == "GFN2L")) {
      GFN2BSet bst;
      for (size_t idx = 0; idx < Natoms; ++idx) {
        atomnr[idx] = bst.GetNAOs(atomnum[idx]);
      }
    }
    return atomnr;
  }
  template<class T>
  void fetchQN(T & basis, matrixE & NQN, matrixE & LQN, std::vector<size_t> & anumbers) {
    //where we actually fetch quantum numbers
    matrix<size_t> QNumbers(1,1);
    size_t Natoms = anumbers.size();
    size_t nshells = NQN.cols();
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      QNumbers = basis.QuantumNumbers(anumbers[idAtm]);
      for (size_t idshell = 0; idshell < nshells; ++idshell) {
        NQN(idAtm + 1,idshell + 1) = 0.0;
        LQN(idAtm + 1,idshell + 1) = 0.0;
      }
      for (size_t idshell = 0; idshell < QNumbers.rows(); ++idshell) {
        NQN(idAtm + 1,idshell + 1) = double(QNumbers(idshell + 1,1));
        LQN(idAtm + 1,idshell + 1) = double(QNumbers(idshell + 1,2));
      }
    }
  }
  void getQuantumNumbers(matrixE & NQN, matrixE & LQN, std::vector<size_t> & anumbers) {
    //function that returns matrices with the quantum numbers for each atom
    if (bset == "slater") {
      SlaterBSet bst;
      fetchQN(bst,NQN,LQN,anumbers);
    }
    else if (bset == "pople") {
      PopleBSet bst;
      fetchQN(bst,NQN,LQN,anumbers);
    }
    else if ((bset == "pm6")||(bset == "PM6")) {
      PM6BSet bst;
      fetchQN(bst,NQN,LQN,anumbers);
    }
    else if ((bset == "pm3")||(bset == "PM3")) {
      PM3BSet bst;
      fetchQN(bst,NQN,LQN,anumbers);
    }
    else if ((bset == "pm3pddg")||(bset == "PM3PDDG")) {
      PM3PDDGBSet bst;
      fetchQN(bst,NQN,LQN,anumbers);
    }
    else if ((bset == "mndopddg")||(bset == "MNDOPDDG")) {
      MNDOPDDGBSet bst;
      fetchQN(bst,NQN,LQN,anumbers);
    }
    else if ((bset == "pm3bp")||(bset == "PM3BP")) {
      PM3BPBSet bst;
      fetchQN(bst,NQN,LQN,anumbers);
    }
    else if ((bset == "rm1")||(bset == "RM1")) {
      RM1BSet bst;
      fetchQN(bst,NQN,LQN,anumbers);
    }
    else if ((bset == "am1")||(bset == "AM1")) {
      AM1BSet bst;
      fetchQN(bst,NQN,LQN,anumbers);
    }
    else if ((bset == "mndo")||(bset == "MNDO")) {
      MNDOBSet bst;
      fetchQN(bst,NQN,LQN,anumbers);
    }
    else if ((bset == "mndod")||(bset == "MNDOd")) {
      MNDOdBSet bst;
      fetchQN(bst,NQN,LQN,anumbers);
    }
    else if (bset == "burns") {
      BurnsBSet bst;
      fetchQN(bst,NQN,LQN,anumbers);
    }
    else if ((bset == "gfn2")||(bset == "GFN2")) {
      GFN2BSet bst;
      fetchQN(bst,NQN,LQN,anumbers);
    }
    else if ((bset == "gfn2l")||(bset == "GFN2L")) {
      GFN2BSet bst;
      fetchQN(bst,NQN,LQN,anumbers);
    }
  }
  template<class T>
  void Overlap(T & bst, bool rotate = true) {
    //function to calculate the overlap
    //the orbitals in the matrix follow the atom order as given in the geometry matrix; this is the vector atoms
    sao = Identity(NAOs);
    double aux = 0.0;
    double aux2;
    double RAmax;
    double Raux;
    size_t position = 1;
    size_t pos2 = 1;
    size_t pos3 = 1;
    size_t atomnaos = 0;
    size_t nao1;
    size_t nao2;
    for (size_t idx = 0; idx < Natoms; ++idx) {
      //calculate 1-center overlaps
      RAmax = bst.ZeroOverlap(atoms[idx]);
      nao1 = bst.GetNAOs(atoms[idx]);
      QNA = bst.QuantumNumbers(atoms[idx]);
      OneCenterOv(overlaps,nao1,exponents[idx],QNA);
      //copy to overlap
      atomnaos = bst.GetNAOs(atoms[idx]);
      for (size_t idx2 = 0; idx2 < atomnaos; ++idx2) {
        for (size_t idx3 = 0; idx3 < atomnaos; ++idx3) {
          sao(idx2 + position,idx3 + position) = overlaps(idx2 + 1,idx3 + 1);
        }
      }
      pos3 = position;
      position += atomnaos;
      //calculate 2-center overlaps
      pos2 = position;
      for (size_t idx2 = idx + 1; idx2 < Natoms; ++idx2) {
        nao2 = bst.GetNAOs(atoms[idx2]);
        //getting the orientation vector
        diffAtom[0] = geometry(idx2 + 1,1) - geometry(idx + 1,1);
        diffAtom[1] = geometry(idx2 + 1,2) - geometry(idx + 1,2);
        diffAtom[2] = geometry(idx2 + 1,3) - geometry(idx + 1,3);
        //normalizing it
        aux = diffAtom[0]*diffAtom[0] + diffAtom[1]*diffAtom[1] + diffAtom[2]*diffAtom[2];
        Raux = fmax(RAmax,bst.ZeroOverlap(atoms[idx2]));          //this here is an approximation, but works conservatively well
        if (aux < Raux*Raux) {
          aux2 = 1.0/sqrt(aux);
          diffAtom[0] *= aux2;
          diffAtom[1] *= aux2;
          diffAtom[2] *= aux2;
          aux2 = sqrt(aux)*dist_Angstrom2aum1;
          QNB = bst.QuantumNumbers(atoms[idx2]);
          //overlaps
          TwoCenterOv(overlaps,nao1,nao2,aux2,exponents[idx],QNA,exponents[idx2],QNB,diffAtom,rotate);
          for (size_t idx3 = 0; idx3 < overlaps.rows(); ++idx3) {
            for (size_t idx4 = 0; idx4 < overlaps.cols(); ++idx4) {
              sao(idx3 + pos3,idx4 + pos2) = overlaps(idx3 + 1,idx4 + 1);
              sao(idx4 + pos2,idx3 + pos3) = overlaps(idx3 + 1,idx4 + 1);
            }
          }
        }
        pos2 += nao2;
      }
    }
  }
  template<class T>
  void OverlapGTO(T & bst, double thresholdzero = 1.0e-8) {
    //function to calculate the overlap involving GTOs
    //the orbitals in the matrix follow the atom order as given in the geometry matrix; this is the vector atoms
    sao = Identity(NAOs);
    TransfO(transfo);                           //this is only needed for d orbitals, however, since it does not require much space nor computation time, calculated once every time
    double aux = 0.0;
    double aux2;
    size_t position = 1;
    size_t pos2 = 1;
    size_t pos3 = 1;
    size_t naosA;
    size_t naosB;
    size_t nshellsA;
    size_t nshellsB;
    size_t nrows;
    size_t ncols;
    int icntA = 1;
    int icntB = 1;
    int LshiftA;
    int LshiftB;
    int NorbA;
    int NorbB;
    int LorbA;
    int LorbB;
    int nprimLA;
    int nprimLB;
    dp[0] = 0.0;
    dp[1] = 0.0;
    dp[2] = 0.0;
    for (size_t atmA = 0; atmA < Natoms; ++atmA) {
      //calculate 1-center overlaps
      naosA = bst.GetNAOs(atoms[atmA]);
      QNA = bst.QuantumNumbers(atoms[atmA]);
      nshellsA = QNA.rows();
      for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
        rA[idcoord - 1] = geometry(atmA + 1,idcoord)*dist_Angstrom2aum1;
      }
      pos3 = position;
      position += naosA;
      //calculate 2-center overlaps
      pos2 = position;
      icntB = 1;
      for (size_t atmB = 0; atmB < atmA; ++atmB) {
        for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
          rB[idcoord - 1] = geometry(atmB + 1,idcoord)*dist_Angstrom2aum1;
        }
        naosB = bst.GetNAOs(atoms[atmB]);
        QNB = bst.QuantumNumbers(atoms[atmB]);
        nshellsB = QNB.rows();
        //now loop over the atomic shells
        LshiftA = 0;
        for (size_t ishellA = 0; ishellA < nshellsA; ++ishellA) {
          NorbA = QNA(ishellA + 1,1);
          LorbA = QNA(ishellA + 1,2);
          nprimLA = bst.getNumberPrimitives(atoms[atmA],LorbA);
          bst.slater2Gauss(nprimLA,NorbA,LorbA,exponents[atmA][ishellA],alpha,CoeffsA,true,thresholdzero);
          LshiftB = 0;
          for (size_t ishellB = 0; ishellB < nshellsB; ++ishellB) {
            NorbB = QNB(ishellB + 1,1);
            LorbB = QNB(ishellB + 1,2);
            nprimLB = bst.getNumberPrimitives(atoms[atmB],LorbB);
            bst.slater2Gauss(nprimLB,NorbB,LorbB,exponents[atmB][ishellB],beta,CoeffsB,true,thresholdzero);
            //overlaps
            TwoCenterOvBlock(overlaps,nprimLA,nprimLB,rA,alpha,CoeffsA,LorbA,rB,beta,CoeffsB,LorbB,dp,ints,val,vi,vj,tk);
            dTransf(overlaps,LorbA,LorbB,transfo,saux);
            nrows = overlaps.rows();
            ncols = overlaps.cols();
            for (size_t idpos1 = 0; idpos1 < nrows; ++idpos1) {
              for (size_t idpos2 = 0; idpos2 < ncols; ++idpos2) {
                sao(idpos1 + icntA + LshiftA,idpos2 + icntB + LshiftB) = overlaps(idpos1 + 1,idpos2 + 1);
                sao(idpos2 + icntB + LshiftB,idpos1 + icntA + LshiftA) = overlaps(idpos1 + 1,idpos2 + 1);
              }
            }
            LshiftB += 2*LorbB + 1;
          }
          LshiftA += 2*LorbA + 1;
        }
        icntB += naosB;
      }               //atmB
      icntA += naosA;
    }                 //atmA
  }
  template<class T>
  void DipQuadrGTO(T & bst, double thresholdzero = 1.0e-8) {
    //function to calculate the dipole and quadrupole moment tensors involving GTOs
    //the orbitals in the matrix follow the atom order as given in the geometry matrix; this is the vector atoms
    dip.resize(3,NAOs,NAOs);
    dip.zero();
    quadr.resize(6,NAOs,NAOs);
    quadr.zero();
    tensor<matrixE,double> auxdip(3,1,1);
    tensor<matrixE,double> auxquadr(6,1,1);
    TransfO(transfo);                           //this is only needed for d orbitals, however, since it does not require much space nor computation time, calculated once every time
    double aux = 0.0;
    double aux2;
    size_t position = 1;
    size_t pos2 = 1;
    size_t pos3 = 1;
    size_t naosA;
    size_t naosB;
    size_t nshellsA;
    size_t nshellsB;
    size_t nrows;
    size_t ncols;
    int icntA = 1;
    int icntB = 1;
    int LshiftA;
    int LshiftB;
    int NorbA;
    int NorbB;
    int LorbA;
    int LorbB;
    int nprimLA;
    int nprimLB;
    dp[0] = 0.0;
    dp[1] = 0.0;
    dp[2] = 0.0;
    for (size_t atmA = 0; atmA < Natoms; ++atmA) {
      //calculate 1-center overlaps
      naosA = bst.GetNAOs(atoms[atmA]);
      QNA = bst.QuantumNumbers(atoms[atmA]);
      nshellsA = QNA.rows();
      for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
        rA[idcoord - 1] = geometry(atmA + 1,idcoord)*dist_Angstrom2aum1;
      }
      LshiftA = 0;
      for (size_t ishellA = 0; ishellA < nshellsA; ++ishellA) {
        NorbA = QNA(ishellA + 1,1);
        LorbA = QNA(ishellA + 1,2);
        nprimLA = bst.getNumberPrimitives(atoms[atmA],LorbA);
        bst.slater2Gauss(nprimLA,NorbA,LorbA,exponents[atmA][ishellA],alpha,CoeffsA,true,thresholdzero);
        LshiftB = 0;
        for (size_t ishellB = 0; ishellB < nshellsA; ++ishellB) {
          NorbB = QNA(ishellB + 1,1);
          LorbB = QNA(ishellB + 1,2);
          nprimLB = bst.getNumberPrimitives(atoms[atmA],LorbB);
          bst.slater2Gauss(nprimLB,NorbB,LorbB,exponents[atmA][ishellB],beta,CoeffsB,true,thresholdzero);
          //integral block
          TwoCenterDipQuadrBlock(auxdip,auxquadr,nprimLA,nprimLB,rA,alpha,CoeffsA,LorbA,rA,beta,CoeffsB,LorbB,dp,ints,val,vi,vj,tk);
          dTransfT(auxdip,LorbA,LorbB,transfo,saux);
          dTransfT(auxquadr,LorbA,LorbB,transfo,saux);
          nrows = auxdip.dim2();
          ncols = auxdip.dim3();
          for (size_t idpos1 = 0; idpos1 < nrows; ++idpos1) {
            for (size_t idpos2 = 0; idpos2 < ncols; ++idpos2) {
              for (size_t intblock = 1; intblock < 4; ++intblock) {
                dip(intblock,idpos1 + icntA + LshiftA,idpos2 + icntA + LshiftB) = auxdip(intblock,idpos1 + 1,idpos2 + 1);
                dip(intblock,idpos2 + icntA + LshiftB,idpos1 + icntA + LshiftA) = auxdip(intblock,idpos1 + 1,idpos2 + 1);
              }
              for (size_t intblock = 1; intblock < 7; ++intblock) {
                quadr(intblock,idpos1 + icntA + LshiftA,idpos2 + icntA + LshiftB) = auxquadr(intblock,idpos1 + 1,idpos2 + 1);
                quadr(intblock,idpos2 + icntA + LshiftB,idpos1 + icntA + LshiftA) = auxquadr(intblock,idpos1 + 1,idpos2 + 1);
              }
            }
          }
          LshiftB += 2*LorbB + 1;
        }
        LshiftA += 2*LorbA + 1;
      }
      pos3 = position;
      position += naosA;
      //calculate 2-center integrals
      pos2 = position;
      icntB = 1;
      for (size_t atmB = 0; atmB < atmA; ++atmB) {
        for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
          rB[idcoord - 1] = geometry(atmB + 1,idcoord)*dist_Angstrom2aum1;
        }
        naosB = bst.GetNAOs(atoms[atmB]);
        QNB = bst.QuantumNumbers(atoms[atmB]);
        nshellsB = QNB.rows();
        //now loop over the atomic shells
        LshiftA = 0;
        for (size_t ishellA = 0; ishellA < nshellsA; ++ishellA) {
          NorbA = QNA(ishellA + 1,1);
          LorbA = QNA(ishellA + 1,2);
          nprimLA = bst.getNumberPrimitives(atoms[atmA],LorbA);
          bst.slater2Gauss(nprimLA,NorbA,LorbA,exponents[atmA][ishellA],alpha,CoeffsA,true,thresholdzero);
          LshiftB = 0;
          for (size_t ishellB = 0; ishellB < nshellsB; ++ishellB) {
            NorbB = QNB(ishellB + 1,1);
            LorbB = QNB(ishellB + 1,2);
            nprimLB = bst.getNumberPrimitives(atoms[atmB],LorbB);
            bst.slater2Gauss(nprimLB,NorbB,LorbB,exponents[atmB][ishellB],beta,CoeffsB,true,thresholdzero);
            TwoCenterDipQuadrBlock(auxdip,auxquadr,nprimLA,nprimLB,rA,alpha,CoeffsA,LorbA,rB,beta,CoeffsB,LorbB,dp,ints,val,vi,vj,tk);
            dTransfT(auxdip,LorbA,LorbB,transfo,saux);
            dTransfT(auxquadr,LorbA,LorbB,transfo,saux);
            nrows = auxdip.dim2();
            ncols = auxdip.dim3();
            for (size_t idpos1 = 0; idpos1 < nrows; ++idpos1) {
              for (size_t idpos2 = 0; idpos2 < ncols; ++idpos2) {
                for (size_t intblock = 1; intblock < 4; ++intblock) {
                  dip(intblock,idpos1 + icntA + LshiftA,idpos2 + icntB + LshiftB) = auxdip(intblock,idpos1 + 1,idpos2 + 1);
                  dip(intblock,idpos2 + icntB + LshiftB,idpos1 + icntA + LshiftA) = auxdip(intblock,idpos1 + 1,idpos2 + 1);
                }
                for (size_t intblock = 1; intblock < 7; ++intblock) {
                  quadr(intblock,idpos1 + icntA + LshiftA,idpos2 + icntB + LshiftB) = auxquadr(intblock,idpos1 + 1,idpos2 + 1);
                  quadr(intblock,idpos2 + icntB + LshiftB,idpos1 + icntA + LshiftA) = auxquadr(intblock,idpos1 + 1,idpos2 + 1);
                }
              }
            }
            LshiftB += 2*LorbB + 1;
          }
          LshiftA += 2*LorbA + 1;
        }
        icntB += naosB;
      }               //atmB
      icntA += naosA;
    }                 //atmA
  }
  template<class T>
  void OvDipQuadrGTO(T & bst, double thresholdzero = 1.0e-8) {
    //function to calculate the overlap, dipole and quadrupole moment tensors involving GTOs
    //the orbitals in the matrix follow the atom order as given in the geometry matrix; this is the vector atoms
    sao = Identity(NAOs);
    dip.resize(3,NAOs,NAOs);
    dip.zero();
    quadr.resize(6,NAOs,NAOs);
    quadr.zero();
    tensor<matrixE,double> auxdip(3,1,1);
    tensor<matrixE,double> auxquadr(6,1,1);
    TransfO(transfo);                           //this is only needed for d orbitals, however, since it does not require much space nor computation time, calculated once every time
    double aux = 0.0;
    double aux2;
    size_t position = 1;
    size_t pos2 = 1;
    size_t pos3 = 1;
    size_t naosA;
    size_t naosB;
    size_t nshellsA;
    size_t nshellsB;
    size_t nrows;
    size_t ncols;
    int icntA = 1;
    int icntB = 1;
    int LshiftA;
    int LshiftB;
    int NorbA;
    int NorbB;
    int LorbA;
    int LorbB;
    int nprimLA;
    int nprimLB;
    dp[0] = 0.0;
    dp[1] = 0.0;
    dp[2] = 0.0;
    for (size_t atmA = 0; atmA < Natoms; ++atmA) {
      //calculate 1-center integrals
      naosA = bst.GetNAOs(atoms[atmA]);
      QNA = bst.QuantumNumbers(atoms[atmA]);
      nshellsA = QNA.rows();
      for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
        rA[idcoord - 1] = geometry(atmA + 1,idcoord)*dist_Angstrom2aum1;
      }
      LshiftA = 0;
      for (size_t ishellA = 0; ishellA < nshellsA; ++ishellA) {
        NorbA = QNA(ishellA + 1,1);
        LorbA = QNA(ishellA + 1,2);
        nprimLA = bst.getNumberPrimitives(atoms[atmA],LorbA);
        bst.slater2Gauss(nprimLA,NorbA,LorbA,exponents[atmA][ishellA],alpha,CoeffsA,true,thresholdzero);
        LshiftB = 0;
        for (size_t ishellB = 0; ishellB < nshellsA; ++ishellB) {
          NorbB = QNA(ishellB + 1,1);
          LorbB = QNA(ishellB + 1,2);
          nprimLB = bst.getNumberPrimitives(atoms[atmA],LorbB);
          bst.slater2Gauss(nprimLB,NorbB,LorbB,exponents[atmA][ishellB],beta,CoeffsB,true,thresholdzero);
          //integral block
          TwoCenterDipQuadrBlock(auxdip,auxquadr,nprimLA,nprimLB,rA,alpha,CoeffsA,LorbA,rA,beta,CoeffsB,LorbB,dp,ints,val,vi,vj,tk);
          dTransfT(auxdip,LorbA,LorbB,transfo,saux);
          dTransfT(auxquadr,LorbA,LorbB,transfo,saux);
          nrows = auxdip.dim2();
          ncols = auxdip.dim3();
          for (size_t idpos1 = 0; idpos1 < nrows; ++idpos1) {
            for (size_t idpos2 = 0; idpos2 < ncols; ++idpos2) {
              for (size_t intblock = 1; intblock < 4; ++intblock) {
                dip(intblock,idpos1 + icntA + LshiftA,idpos2 + icntA + LshiftB) = auxdip(intblock,idpos1 + 1,idpos2 + 1);
                dip(intblock,idpos2 + icntA + LshiftB,idpos1 + icntA + LshiftA) = auxdip(intblock,idpos1 + 1,idpos2 + 1);
              }
              for (size_t intblock = 1; intblock < 7; ++intblock) {
                quadr(intblock,idpos1 + icntA + LshiftA,idpos2 + icntA + LshiftB) = auxquadr(intblock,idpos1 + 1,idpos2 + 1);
                quadr(intblock,idpos2 + icntA + LshiftB,idpos1 + icntA + LshiftA) = auxquadr(intblock,idpos1 + 1,idpos2 + 1);
              }
            }
          }
          LshiftB += 2*LorbB + 1;
        }
        LshiftA += 2*LorbA + 1;
      }
      pos3 = position;
      position += naosA;
      //calculate 2-center integrals
      pos2 = position;
      icntB = 1;
      for (size_t atmB = 0; atmB < atmA; ++atmB) {
        for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
          rB[idcoord - 1] = geometry(atmB + 1,idcoord)*dist_Angstrom2aum1;
        }
        naosB = bst.GetNAOs(atoms[atmB]);
        QNB = bst.QuantumNumbers(atoms[atmB]);
        nshellsB = QNB.rows();
        //now loop over the atomic shells
        LshiftA = 0;
        for (size_t ishellA = 0; ishellA < nshellsA; ++ishellA) {
          NorbA = QNA(ishellA + 1,1);
          LorbA = QNA(ishellA + 1,2);
          nprimLA = bst.getNumberPrimitives(atoms[atmA],LorbA);
          bst.slater2Gauss(nprimLA,NorbA,LorbA,exponents[atmA][ishellA],alpha,CoeffsA,true,thresholdzero);
          LshiftB = 0;
          for (size_t ishellB = 0; ishellB < nshellsB; ++ishellB) {
            NorbB = QNB(ishellB + 1,1);
            LorbB = QNB(ishellB + 1,2);
            nprimLB = bst.getNumberPrimitives(atoms[atmB],LorbB);
            bst.slater2Gauss(nprimLB,NorbB,LorbB,exponents[atmB][ishellB],beta,CoeffsB,true,thresholdzero);
            TwoCenterOvDipQuadrBlock(overlaps,auxdip,auxquadr,nprimLA,nprimLB,rA,alpha,CoeffsA,LorbA,rB,beta,CoeffsB,LorbB,dp,ints,val,vi,vj,tk);
            dTransf(overlaps,LorbA,LorbB,transfo,saux);
            dTransfT(auxdip,LorbA,LorbB,transfo,saux);
            dTransfT(auxquadr,LorbA,LorbB,transfo,saux);
            nrows = overlaps.rows();
            ncols = overlaps.cols();
            for (size_t idpos1 = 0; idpos1 < nrows; ++idpos1) {
              for (size_t idpos2 = 0; idpos2 < ncols; ++idpos2) {
                sao(idpos1 + icntA + LshiftA,idpos2 + icntB + LshiftB) = overlaps(idpos1 + 1,idpos2 + 1);
                sao(idpos2 + icntB + LshiftB,idpos1 + icntA + LshiftA) = overlaps(idpos1 + 1,idpos2 + 1);
                for (size_t intblock = 1; intblock < 4; ++intblock) {
                  dip(intblock,idpos1 + icntA + LshiftA,idpos2 + icntB + LshiftB) = auxdip(intblock,idpos1 + 1,idpos2 + 1);
                  dip(intblock,idpos2 + icntB + LshiftB,idpos1 + icntA + LshiftA) = auxdip(intblock,idpos1 + 1,idpos2 + 1);
                }
                for (size_t intblock = 1; intblock < 7; ++intblock) {
                  quadr(intblock,idpos1 + icntA + LshiftA,idpos2 + icntB + LshiftB) = auxquadr(intblock,idpos1 + 1,idpos2 + 1);
                  quadr(intblock,idpos2 + icntB + LshiftB,idpos1 + icntA + LshiftA) = auxquadr(intblock,idpos1 + 1,idpos2 + 1);
                }
              }
            }
            LshiftB += 2*LorbB + 1;
          }
          LshiftA += 2*LorbA + 1;
        }
        icntB += naosB;
      }               //atmB
      icntA += naosA;
    }                 //atmA
  }
  void CAO2SAO(matrix<matrixE> & gSDQ, matrixE & temp, int ishellA, int ishellB) {
    //function that converts a big block matrix from CAO to SAO
    size_t nr = gSDQ.rows() + 1;
    size_t nc = gSDQ.cols() + 1;
    size_t sztmp = temp.rows() + 1;           //temp is supposed to be square
    size_t rtmp;
    size_t ctmp;
    TransfO(transfo);
    for (size_t idk = 1; idk < nc; ++idk) {       //1 S, 2-4 D, 5-10 Q, 11-13 D, 14-19 Q
      for (size_t idcoord = 1; idcoord < nr; ++idcoord) {
        temp = gSDQ(idcoord,idk);
        dTransf(temp,ishellA,ishellB,transfo,saux);
        rtmp = temp.rows();
        ctmp = temp.cols();
        temp.resize(sztmp - 1,sztmp - 1);
        gSDQ(idcoord,idk) = temp;
      }
    }
  }
  template<class T>
  void gOvDipQuadrGTO(T & bst, std::vector<matrixE> & SDQ, matrix<matrixE> & gSDQ, std::vector<double> & rA, std::vector<double> & rB, std::vector<double> & rAB, double rAB2, size_t isA, size_t isB, int atmA, int atmB, double thresholdzero = 1.0e-8, double maxsqdist = 2000.0, double intcut = 25.0) {
    //function to calculate the dipole and quadrupole moment tensors involving GTOs
    QNA = bst.QuantumNumbers(atoms[atmA]);
    QNB = bst.QuantumNumbers(atoms[atmB]);
    int NorbA = QNA(isA + 1,1);
    int LorbA = QNA(isA + 1,2);
    int NorbB = QNB(isB + 1,1);
    int LorbB = QNB(isB + 1,2);
    int nprimLA = bst.getNumberPrimitives(atoms[atmA],LorbA);
    int nprimLB = bst.getNumberPrimitives(atoms[atmB],LorbB);
    bst.slater2Gauss(nprimLA,NorbA,LorbA,exponents[atmA][isA],alpha,CoeffsA,true,thresholdzero);
    bst.slater2Gauss(nprimLB,NorbB,LorbB,exponents[atmB][isB],beta,CoeffsB,true,thresholdzero);
    gTwoCenterOvDipQuadrBlock(SDQ,gSDQ,rA,rB,rP,rAB,rAB2,nprimLA,nprimLB,alpha,beta,CoeffsA,CoeffsB,ints,tk,gints,val,grad,vi,vj,vv,gi,gg,LorbA,LorbB,maxsqdist,intcut);
  }
  void CAMM(matrixE & dens, matrixE & QShells, matrixE & muAtoms, matrixE & thetaAtoms) {
    //function to calculate shellwise multipolar properties, to avoid passing dipole and quadrupole moment integrals
    if ((bset == "gfn2")||(bset == "GFN2")) {
      GFN2BSet bst;
      double ps;
      int totalshells = QShells.rows();
      int atomA;
      int atomB;
      int nshellsA;
      int nshellsB;
      int LorbA;
      int LorbB;
      int pos;
      int iao = 0;
      int jao = 0;
      int ishell = 1;
      int jshell = 1;
      for (int iAtm = 0; iAtm < Natoms; ++iAtm) {
        atomA = atoms[iAtm];
        QNA = bst.QuantumNumbers(atomA);
        nshellsA = QNA.rows();
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          rA[idcoord] = geometry(iAtm + 1,idcoord + 1)*dist_Angstrom2aum1;
          muAtoms(iAtm + 1,idcoord + 1) = 0.0;
          thetaAtoms(iAtm + 1,idcoord + 1) = 0.0;
          thetaAtoms(iAtm + 1,idcoord + 4) = 0.0;
        }
        for (int ishellA = 0; ishellA < nshellsA; ++ishellA) {
          QShells(ishell + ishellA,1) = 0.0;
        }
        for (int ishellA = 0; ishellA < nshellsA; ++ishellA, ++ishell) {
          LorbA = QNA(ishellA + 1,2);
          jao = 0;
          jshell = 1;
          for (int iBtm = 0; iBtm < iAtm; ++iBtm) {
            for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
              rB[idcoord] = geometry(iBtm + 1,idcoord + 1)*dist_Angstrom2aum1;
            }
            atomB = atoms[iBtm];
            QNB = bst.QuantumNumbers(atomB);
            nshellsB = QNB.rows();
            for (int ishellB = 0; ishellB < nshellsB; ++ishellB, ++jshell) {
              LorbB = QNB(ishellB + 1,2);
              //zero out where to gather data
              ps = 0.0;
              for (size_t idx = 0; idx < 3; ++idx) {
                dp[idx] = 0.0;
                qp[idx] = 0.0;
                qp[3 + idx] = 0.0;
              }
              for (int idLA = 0; idLA < 2*LorbA + 1; ++idLA) {
                for (int idLB = 0; idLB < 2*LorbB + 1; ++idLB) {
                  ps += dens(iao + idLA + 1,jao + idLB + 1)*sao(iao + idLA + 1,jao + idLB + 1);
                  for (size_t idx = 0; idx < 3; ++idx) {
                    dp[idx] += dens(iao + idLA + 1,jao + idLB + 1)*dip(idx + 1,iao + idLA + 1,jao + idLB + 1);
                    qp[idx] += dens(iao + idLA + 1,jao + idLB + 1)*quadr(idx + 1,iao + idLA + 1,jao + idLB + 1);
                    qp[3 + idx] += dens(iao + idLA + 1,jao + idLB + 1)*quadr(idx + 4,iao + idLA + 1,jao + idLB + 1);
                  }
                }
              }
              QShells(ishell,1) -= ps;
              QShells(jshell,1) -= ps;
              for (size_t idx = 0; idx < 3; ++idx) {
                pos = (idx + 1)*(idx + 2)/2;
                muAtoms(iAtm + 1,idx + 1) += rA[idx]*ps - dp[idx];
                muAtoms(iBtm + 1,idx + 1) += rB[idx]*ps - dp[idx];
                thetaAtoms(iAtm + 1,pos) += 2.0*rA[idx]*dp[idx] - rA[idx]*rA[idx]*ps - qp[idx];
                thetaAtoms(iBtm + 1,pos) += 2.0*rB[idx]*dp[idx] - rB[idx]*rB[idx]*ps - qp[idx];
                for (size_t idy = 0; idy < idx; ++idy) {
                  pos = (idx + 1)*idx/2 + idy + 1;
                  thetaAtoms(iAtm + 1,pos) += rA[idx]*dp[idy] + rA[idy]*dp[idx] - rA[idx]*rA[idy]*ps - qp[2 + idx + idy];
                  thetaAtoms(iBtm + 1,pos) += rB[idx]*dp[idy] + rB[idy]*dp[idx] - rB[idx]*rB[idy]*ps - qp[2 + idx + idy];
                }
              }
              jao += 2*LorbB + 1;
            }
          }
          //rest of A
          for (int ishellB = 0; ishellB < ishellA + 1; ++ishellB, ++jshell) {
            LorbB = QNA(ishellB + 1,2);
            ps = 0.0;
            for (size_t idx = 0; idx < 3; ++idx) {
              dp[idx] = 0.0;
              qp[idx] = 0.0;
              qp[3 + idx] = 0.0;
            }
            for (int idLA = 0; idLA < 2*LorbA + 1; ++idLA) {
              for (int idLB = 0; idLB < 2*LorbB + 1; ++idLB) {
                ps += dens(iao + idLA + 1,jao + idLB + 1)*sao(iao + idLA + 1,jao + idLB + 1);
                for (size_t idx = 0; idx < 3; ++idx) {
                  dp[idx] += dens(iao + idLA + 1,jao + idLB + 1)*dip(idx + 1,iao + idLA + 1,jao + idLB + 1);
                  qp[idx] += dens(iao + idLA + 1,jao + idLB + 1)*quadr(idx + 1,iao + idLA + 1,jao + idLB + 1);
                  qp[3 + idx] += dens(iao + idLA + 1,jao + idLB + 1)*quadr(idx + 4,iao + idLA + 1,jao + idLB + 1);
                }
              }
            }
            QShells(ishell,1) -= ps;
            QShells(jshell,1) -= ps*double(ishell != jshell);
            for (size_t idx = 0; idx < 3; ++idx) {
              pos = (idx + 1)*(idx + 2)/2;
              muAtoms(iAtm + 1,idx + 1) += rA[idx]*ps - dp[idx];
              muAtoms(iAtm + 1,idx + 1) += (rA[idx]*ps - dp[idx])*double(ishell != jshell);
              thetaAtoms(iAtm + 1,pos) += 2.0*rA[idx]*dp[idx] - rA[idx]*rA[idx]*ps - qp[idx];
              thetaAtoms(iAtm + 1,pos) += (2.0*rA[idx]*dp[idx] - rA[idx]*rA[idx]*ps - qp[idx])*double(ishell != jshell);
              for (size_t idy = 0; idy < idx; ++idy) {
                pos = (idx + 1)*idx/2 + idy + 1;
                thetaAtoms(iAtm + 1,pos) += rA[idx]*dp[idy] + rA[idy]*dp[idx] - rA[idx]*rA[idy]*ps - qp[2 + idx + idy];
                thetaAtoms(iAtm + 1,pos) += (rA[idx]*dp[idy] + rA[idy]*dp[idx] - rA[idx]*rA[idy]*ps - qp[2 + idx + idy])*double(ishell != jshell);
              }
            }
            jao += 2*LorbB + 1;
          }
          iao += 2*LorbA + 1;
        }
      }
    }
  }
  void CAMMdirect(matrixE & dens, matrixE & QShells, matrixE & muAtoms, matrixE & thetaAtoms) {
    //function to calculate shellwise multipolar properties, to avoid passing dipole and quadrupole moment integrals
    if ((bset == "gfn2l")||(bset == "GFN2L")) {
      GFN2BSet bst;
      double ps;
      int totalshells = QShells.rows();
      int atomA;
      int atomB;
      int nshellsA;
      int nshellsB;
      int NorbA;
      int NorbB;
      int nprimLA;
      int nprimLB;
      int LorbA;
      int LorbB;
      int pos;
      int iao = 0;
      int jao = 0;
      int ishell = 1;
      int jshell = 1;
      double thresholdzero = 1.0e-8;
      dip.resize(3,1,1);
      quadr.resize(6,1,1);
      dp2[0] = 0.0;
      dp2[1] = 0.0;
      dp2[2] = 0.0;
      for (int iAtm = 0; iAtm < Natoms; ++iAtm) {
        atomA = atoms[iAtm];
        QNA = bst.QuantumNumbers(atomA);
        nshellsA = QNA.rows();
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          rA[idcoord] = geometry(iAtm + 1,idcoord + 1)*dist_Angstrom2aum1;
          muAtoms(iAtm + 1,idcoord + 1) = 0.0;
          thetaAtoms(iAtm + 1,idcoord + 1) = 0.0;
          thetaAtoms(iAtm + 1,idcoord + 4) = 0.0;
        }
        for (int ishellA = 0; ishellA < nshellsA; ++ishellA) {
          QShells(ishell + ishellA,1) = 0.0;
        }
        for (int ishellA = 0; ishellA < nshellsA; ++ishellA, ++ishell) {
          NorbA = QNA(ishellA + 1,1);
          LorbA = QNA(ishellA + 1,2);
          nprimLA = bst.getNumberPrimitives(atoms[iAtm],LorbA);
          bst.slater2Gauss(nprimLA,NorbA,LorbA,exponents[iAtm][ishellA],alpha,CoeffsA,true,thresholdzero);
          jao = 0;
          jshell = 1;
          for (int iBtm = 0; iBtm < iAtm; ++iBtm) {
            for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
              rB[idcoord] = geometry(iBtm + 1,idcoord + 1)*dist_Angstrom2aum1;
            }
            atomB = atoms[iBtm];
            QNB = bst.QuantumNumbers(atomB);
            nshellsB = QNB.rows();
            for (int ishellB = 0; ishellB < nshellsB; ++ishellB, ++jshell) {
              LorbB = QNB(ishellB + 1,2);
              //integral block
              NorbB = QNB(ishellB + 1,1);
              LorbB = QNB(ishellB + 1,2);
              nprimLB = bst.getNumberPrimitives(atoms[iBtm],LorbB);
              bst.slater2Gauss(nprimLB,NorbB,LorbB,exponents[iBtm][ishellB],beta,CoeffsB,true,thresholdzero);
              TwoCenterDipQuadrBlock(dip,quadr,nprimLA,nprimLB,rA,alpha,CoeffsA,LorbA,rB,beta,CoeffsB,LorbB,dp2,ints,val,vi,vj,tk);
              dTransfT(dip,LorbA,LorbB,transfo,saux);
              dTransfT(quadr,LorbA,LorbB,transfo,saux);
              //zero out where to gather data
              ps = 0.0;
              for (size_t idx = 0; idx < 3; ++idx) {
                dp[idx] = 0.0;
                qp[idx] = 0.0;
                qp[3 + idx] = 0.0;
              }
              for (int idLA = 0; idLA < 2*LorbA + 1; ++idLA) {
                for (int idLB = 0; idLB < 2*LorbB + 1; ++idLB) {
                  ps += dens(iao + idLA + 1,jao + idLB + 1)*sao(iao + idLA + 1,jao + idLB + 1);
                  for (size_t idx = 0; idx < 3; ++idx) {
                    dp[idx] += dens(iao + idLA + 1,jao + idLB + 1)*dip(idx + 1,idLA + 1,idLB + 1);
                    qp[idx] += dens(iao + idLA + 1,jao + idLB + 1)*quadr(idx + 1,idLA + 1,idLB + 1);
                    qp[3 + idx] += dens(iao + idLA + 1,jao + idLB + 1)*quadr(idx + 4,idLA + 1,idLB + 1);
                  }
                }
              }
              QShells(ishell,1) -= ps;
              QShells(jshell,1) -= ps;
              for (size_t idx = 0; idx < 3; ++idx) {
                pos = (idx + 1)*(idx + 2)/2;
                muAtoms(iAtm + 1,idx + 1) += rA[idx]*ps - dp[idx];
                muAtoms(iBtm + 1,idx + 1) += rB[idx]*ps - dp[idx];
                thetaAtoms(iAtm + 1,pos) += 2.0*rA[idx]*dp[idx] - rA[idx]*rA[idx]*ps - qp[idx];
                thetaAtoms(iBtm + 1,pos) += 2.0*rB[idx]*dp[idx] - rB[idx]*rB[idx]*ps - qp[idx];
                for (size_t idy = 0; idy < idx; ++idy) {
                  pos = (idx + 1)*idx/2 + idy + 1;
                  thetaAtoms(iAtm + 1,pos) += rA[idx]*dp[idy] + rA[idy]*dp[idx] - rA[idx]*rA[idy]*ps - qp[2 + idx + idy];
                  thetaAtoms(iBtm + 1,pos) += rB[idx]*dp[idy] + rB[idy]*dp[idx] - rB[idx]*rB[idy]*ps - qp[2 + idx + idy];
                }
              }
              jao += 2*LorbB + 1;
            }
          }
          //rest of A
          for (int ishellB = 0; ishellB < ishellA + 1; ++ishellB, ++jshell) {
            NorbB = QNA(ishellB + 1,1);
            LorbB = QNA(ishellB + 1,2);
            nprimLB = bst.getNumberPrimitives(atoms[iAtm],LorbB);
            bst.slater2Gauss(nprimLB,NorbB,LorbB,exponents[iAtm][ishellB],beta,CoeffsB,true,thresholdzero);
            //integral block
            TwoCenterDipQuadrBlock(dip,quadr,nprimLA,nprimLB,rA,alpha,CoeffsA,LorbA,rA,beta,CoeffsB,LorbB,dp2,ints,val,vi,vj,tk);
            dTransfT(dip,LorbA,LorbB,transfo,saux);
            dTransfT(quadr,LorbA,LorbB,transfo,saux);
            ps = 0.0;
            for (size_t idx = 0; idx < 3; ++idx) {
              dp[idx] = 0.0;
              qp[idx] = 0.0;
              qp[3 + idx] = 0.0;
            }
            for (int idLA = 0; idLA < 2*LorbA + 1; ++idLA) {
              for (int idLB = 0; idLB < 2*LorbB + 1; ++idLB) {
                ps += dens(iao + idLA + 1,jao + idLB + 1)*sao(iao + idLA + 1,jao + idLB + 1);
                for (size_t idx = 0; idx < 3; ++idx) {
                  dp[idx] += dens(iao + idLA + 1,jao + idLB + 1)*dip(idx + 1,idLA + 1,idLB + 1);
                  qp[idx] += dens(iao + idLA + 1,jao + idLB + 1)*quadr(idx + 1,idLA + 1,idLB + 1);
                  qp[3 + idx] += dens(iao + idLA + 1,jao + idLB + 1)*quadr(idx + 4,idLA + 1,idLB + 1);
                }
              }
            }
            QShells(ishell,1) -= ps;
            QShells(jshell,1) -= ps*double(ishell != jshell);
            for (size_t idx = 0; idx < 3; ++idx) {
              pos = (idx + 1)*(idx + 2)/2;
              muAtoms(iAtm + 1,idx + 1) += rA[idx]*ps - dp[idx];
              muAtoms(iAtm + 1,idx + 1) += (rA[idx]*ps - dp[idx])*double(ishell != jshell);
              thetaAtoms(iAtm + 1,pos) += 2.0*rA[idx]*dp[idx] - rA[idx]*rA[idx]*ps - qp[idx];
              thetaAtoms(iAtm + 1,pos) += (2.0*rA[idx]*dp[idx] - rA[idx]*rA[idx]*ps - qp[idx])*double(ishell != jshell);
              for (size_t idy = 0; idy < idx; ++idy) {
                pos = (idx + 1)*idx/2 + idy + 1;
                thetaAtoms(iAtm + 1,pos) += rA[idx]*dp[idy] + rA[idy]*dp[idx] - rA[idx]*rA[idy]*ps - qp[2 + idx + idy];
                thetaAtoms(iAtm + 1,pos) += (rA[idx]*dp[idy] + rA[idy]*dp[idx] - rA[idx]*rA[idy]*ps - qp[2 + idx + idy])*double(ishell != jshell);
              }
            }
            jao += 2*LorbB + 1;
          }
          iao += 2*LorbA + 1;
        }
      }
    }
  }
  void getDipQuadrEldirect(tensor<matrixE,double> & auxdip, tensor<matrixE,double> & auxquadr, GFN2BSet & bst, int atmA, int shellA, int atmB, int shellB, double thresholdzero = 1.0e-8) {
    //function returning a block of elements of the dipole and quadrupole integrals
    QNA = bst.QuantumNumbers(atoms[atmA]);
    QNB = bst.QuantumNumbers(atoms[atmB]);
    int NorbA = QNA(shellA + 1,1);
    int LorbA = QNA(shellA + 1,2);
    int nprimLA = bst.getNumberPrimitives(atoms[atmA],LorbA);
    int NorbB = QNB(shellB + 1,1);
    int LorbB = QNB(shellB + 1,2);
    int nprimLB = bst.getNumberPrimitives(atoms[atmB],LorbB);
    for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
      rA[idcoord - 1] = geometry(atmA + 1,idcoord)*dist_Angstrom2aum1;
      rB[idcoord - 1] = geometry(atmB + 1,idcoord)*dist_Angstrom2aum1;
    }
    bst.slater2Gauss(nprimLA,NorbA,LorbA,exponents[atmA][shellA],alpha,CoeffsA,true,thresholdzero);
    bst.slater2Gauss(nprimLB,NorbB,LorbB,exponents[atmB][shellB],beta,CoeffsB,true,thresholdzero);
    TwoCenterDipQuadrBlock(auxdip,auxquadr,nprimLA,nprimLB,rA,alpha,CoeffsA,LorbA,rB,beta,CoeffsB,LorbB,dp2,ints,val,vi,vj,tk);
    dTransfT(auxdip,LorbA,LorbB,transfo,saux);
    dTransfT(auxquadr,LorbA,LorbB,transfo,saux);
  }
  void getDipQuadrEl(std::vector<double> & dipole, std::vector<double> & quadrupole, size_t iel, size_t jel) {
    //function returning a "row" of elements of the dipole and quadrupole integrals
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      dipole[idcoord] = dip(idcoord + 1,iel,jel);
      quadrupole[idcoord] = quadr(idcoord + 1,iel,jel);
      quadrupole[idcoord + 3] = quadr(idcoord + 4,iel,jel);
    }
  }
  void Overlap_dXCD(matrixE & sao_dxcd, std::string var, size_t atmC, size_t atmD, double RCD, std::vector<double> diffAtom, bool rotate = true) {
    //function to calculate blocks of derivatives of the overlap matrix with respect to some variable
    //the orbitals in the matrix follow the atom order as given in the geometry matrix; this is the vector atoms
    //RCD in atomic units
    if (bset == "slater") {
      SlaterBSet bst;
      if ((var == "")||(var == " ")) {bst.OV2Center(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "1")||(var == "R")||(var == "r")) {bst.OV2Center_dR(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "3")||(var == "P")||(var == "p")||(var == "PHI")||(var == "phi")||(var == "Phi")) {bst.OV2Center_dp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "2")||(var == "T")||(var == "t")||(var == "THETA")||(var == "theta")||(var == "Theta")) {bst.OV2Center_dt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "4")||(var == "RR")||(var == "rr")) {bst.OV2Center_dR2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "6")||(var == "PP")||(var == "pp")||(var == "PHIPHI")||(var == "phiphi")||(var == "PhiPhi")) {bst.OV2Center_dp2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "5")||(var == "TT")||(var == "tt")||(var == "THETATHETA")||(var == "thetatheta")||(var == "ThetaTheta")) {bst.OV2Center_dt2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "8")||(var == "RP")||(var == "rp")||(var == "PR")||(var == "pR")||(var == "phiR")||(var == "Rphi")||(var == "PHIR")||(var == "RPHI")) {bst.OV2Center_dRdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "7")||(var == "RT")||(var == "rt")||(var == "TR")||(var == "tR")||(var == "thetaR")||(var == "Rtheta")||(var == "THETAR")||(var == "RTHETA")) {bst.OV2Center_dRdt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "9")||(var == "PT")||(var == "pt")||(var == "TP")||(var == "tp")||(var == "thetaphi")||(var == "phitheta")||(var == "THETAPHI")||(var == "PHITHETA")) {bst.OV2Center_dtdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
    }
    else if (bset == "pople") {
      PopleBSet bst;
      if ((var == "")||(var == " ")) {bst.OV2Center(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "1")||(var == "R")||(var == "r")) {bst.OV2Center_dR(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "3")||(var == "P")||(var == "p")||(var == "PHI")||(var == "phi")||(var == "Phi")) {bst.OV2Center_dp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "2")||(var == "T")||(var == "t")||(var == "THETA")||(var == "theta")||(var == "Theta")) {bst.OV2Center_dt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "4")||(var == "RR")||(var == "rr")) {bst.OV2Center_dR2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "6")||(var == "PP")||(var == "pp")||(var == "PHIPHI")||(var == "phiphi")||(var == "PhiPhi")) {bst.OV2Center_dp2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "5")||(var == "TT")||(var == "tt")||(var == "THETATHETA")||(var == "thetatheta")||(var == "ThetaTheta")) {bst.OV2Center_dt2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "8")||(var == "RP")||(var == "rp")||(var == "PR")||(var == "pR")||(var == "phiR")||(var == "Rphi")||(var == "PHIR")||(var == "RPHI")) {bst.OV2Center_dRdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "7")||(var == "RT")||(var == "rt")||(var == "TR")||(var == "tR")||(var == "thetaR")||(var == "Rtheta")||(var == "THETAR")||(var == "RTHETA")) {bst.OV2Center_dRdt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "9")||(var == "PT")||(var == "pt")||(var == "TP")||(var == "tp")||(var == "thetaphi")||(var == "phitheta")||(var == "THETAPHI")||(var == "PHITHETA")) {bst.OV2Center_dtdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
    }
    else if ((bset == "pm6")||(bset == "PM6")) {
      PM6BSet bst;
      if ((var == "")||(var == " ")) {bst.OV2Center(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "1")||(var == "R")||(var == "r")) {bst.OV2Center_dR(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "3")||(var == "P")||(var == "p")||(var == "PHI")||(var == "phi")||(var == "Phi")) {bst.OV2Center_dp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "2")||(var == "T")||(var == "t")||(var == "THETA")||(var == "theta")||(var == "Theta")) {bst.OV2Center_dt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "4")||(var == "RR")||(var == "rr")) {bst.OV2Center_dR2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "6")||(var == "PP")||(var == "pp")||(var == "PHIPHI")||(var == "phiphi")||(var == "PhiPhi")) {bst.OV2Center_dp2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "5")||(var == "TT")||(var == "tt")||(var == "THETATHETA")||(var == "thetatheta")||(var == "ThetaTheta")) {bst.OV2Center_dt2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "8")||(var == "RP")||(var == "rp")||(var == "PR")||(var == "pR")||(var == "phiR")||(var == "Rphi")||(var == "PHIR")||(var == "RPHI")) {bst.OV2Center_dRdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "7")||(var == "RT")||(var == "rt")||(var == "TR")||(var == "tR")||(var == "thetaR")||(var == "Rtheta")||(var == "THETAR")||(var == "RTHETA")) {bst.OV2Center_dRdt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "9")||(var == "PT")||(var == "pt")||(var == "TP")||(var == "tp")||(var == "thetaphi")||(var == "phitheta")||(var == "THETAPHI")||(var == "PHITHETA")) {bst.OV2Center_dtdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
    }
    else if ((bset == "pm3")||(bset == "PM3")) {
      PM3BSet bst;
      if ((var == "")||(var == " ")) {bst.OV2Center(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "1")||(var == "R")||(var == "r")) {bst.OV2Center_dR(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "3")||(var == "P")||(var == "p")||(var == "PHI")||(var == "phi")||(var == "Phi")) {bst.OV2Center_dp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "2")||(var == "T")||(var == "t")||(var == "THETA")||(var == "theta")||(var == "Theta")) {bst.OV2Center_dt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "4")||(var == "RR")||(var == "rr")) {bst.OV2Center_dR2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "6")||(var == "PP")||(var == "pp")||(var == "PHIPHI")||(var == "phiphi")||(var == "PhiPhi")) {bst.OV2Center_dp2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "5")||(var == "TT")||(var == "tt")||(var == "THETATHETA")||(var == "thetatheta")||(var == "ThetaTheta")) {bst.OV2Center_dt2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "8")||(var == "RP")||(var == "rp")||(var == "PR")||(var == "pR")||(var == "phiR")||(var == "Rphi")||(var == "PHIR")||(var == "RPHI")) {bst.OV2Center_dRdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "7")||(var == "RT")||(var == "rt")||(var == "TR")||(var == "tR")||(var == "thetaR")||(var == "Rtheta")||(var == "THETAR")||(var == "RTHETA")) {bst.OV2Center_dRdt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "9")||(var == "PT")||(var == "pt")||(var == "TP")||(var == "tp")||(var == "thetaphi")||(var == "phitheta")||(var == "THETAPHI")||(var == "PHITHETA")) {bst.OV2Center_dtdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
    }
    else if ((bset == "pm3pddg")||(bset == "PM3PDDG")) {
      PM3PDDGBSet bst;
      if ((var == "")||(var == " ")) {bst.OV2Center(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "1")||(var == "R")||(var == "r")) {bst.OV2Center_dR(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "3")||(var == "P")||(var == "p")||(var == "PHI")||(var == "phi")||(var == "Phi")) {bst.OV2Center_dp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "2")||(var == "T")||(var == "t")||(var == "THETA")||(var == "theta")||(var == "Theta")) {bst.OV2Center_dt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "4")||(var == "RR")||(var == "rr")) {bst.OV2Center_dR2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "6")||(var == "PP")||(var == "pp")||(var == "PHIPHI")||(var == "phiphi")||(var == "PhiPhi")) {bst.OV2Center_dp2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "5")||(var == "TT")||(var == "tt")||(var == "THETATHETA")||(var == "thetatheta")||(var == "ThetaTheta")) {bst.OV2Center_dt2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "8")||(var == "RP")||(var == "rp")||(var == "PR")||(var == "pR")||(var == "phiR")||(var == "Rphi")||(var == "PHIR")||(var == "RPHI")) {bst.OV2Center_dRdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "7")||(var == "RT")||(var == "rt")||(var == "TR")||(var == "tR")||(var == "thetaR")||(var == "Rtheta")||(var == "THETAR")||(var == "RTHETA")) {bst.OV2Center_dRdt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "9")||(var == "PT")||(var == "pt")||(var == "TP")||(var == "tp")||(var == "thetaphi")||(var == "phitheta")||(var == "THETAPHI")||(var == "PHITHETA")) {bst.OV2Center_dtdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
    }
    else if ((bset == "mndopddg")||(bset == "MNDOPDDG")) {
      MNDOPDDGBSet bst;
      if ((var == "")||(var == " ")) {bst.OV2Center(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "1")||(var == "R")||(var == "r")) {bst.OV2Center_dR(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "3")||(var == "P")||(var == "p")||(var == "PHI")||(var == "phi")||(var == "Phi")) {bst.OV2Center_dp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "2")||(var == "T")||(var == "t")||(var == "THETA")||(var == "theta")||(var == "Theta")) {bst.OV2Center_dt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "4")||(var == "RR")||(var == "rr")) {bst.OV2Center_dR2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "6")||(var == "PP")||(var == "pp")||(var == "PHIPHI")||(var == "phiphi")||(var == "PhiPhi")) {bst.OV2Center_dp2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "5")||(var == "TT")||(var == "tt")||(var == "THETATHETA")||(var == "thetatheta")||(var == "ThetaTheta")) {bst.OV2Center_dt2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "8")||(var == "RP")||(var == "rp")||(var == "PR")||(var == "pR")||(var == "phiR")||(var == "Rphi")||(var == "PHIR")||(var == "RPHI")) {bst.OV2Center_dRdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "7")||(var == "RT")||(var == "rt")||(var == "TR")||(var == "tR")||(var == "thetaR")||(var == "Rtheta")||(var == "THETAR")||(var == "RTHETA")) {bst.OV2Center_dRdt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "9")||(var == "PT")||(var == "pt")||(var == "TP")||(var == "tp")||(var == "thetaphi")||(var == "phitheta")||(var == "THETAPHI")||(var == "PHITHETA")) {bst.OV2Center_dtdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
    }
    else if ((bset == "pm3bp")||(bset == "PM3BP")) {
      PM3BPBSet bst;
      if ((var == "")||(var == " ")) {bst.OV2Center(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "1")||(var == "R")||(var == "r")) {bst.OV2Center_dR(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "3")||(var == "P")||(var == "p")||(var == "PHI")||(var == "phi")||(var == "Phi")) {bst.OV2Center_dp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "2")||(var == "T")||(var == "t")||(var == "THETA")||(var == "theta")||(var == "Theta")) {bst.OV2Center_dt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "4")||(var == "RR")||(var == "rr")) {bst.OV2Center_dR2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "6")||(var == "PP")||(var == "pp")||(var == "PHIPHI")||(var == "phiphi")||(var == "PhiPhi")) {bst.OV2Center_dp2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "5")||(var == "TT")||(var == "tt")||(var == "THETATHETA")||(var == "thetatheta")||(var == "ThetaTheta")) {bst.OV2Center_dt2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "8")||(var == "RP")||(var == "rp")||(var == "PR")||(var == "pR")||(var == "phiR")||(var == "Rphi")||(var == "PHIR")||(var == "RPHI")) {bst.OV2Center_dRdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "7")||(var == "RT")||(var == "rt")||(var == "TR")||(var == "tR")||(var == "thetaR")||(var == "Rtheta")||(var == "THETAR")||(var == "RTHETA")) {bst.OV2Center_dRdt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "9")||(var == "PT")||(var == "pt")||(var == "TP")||(var == "tp")||(var == "thetaphi")||(var == "phitheta")||(var == "THETAPHI")||(var == "PHITHETA")) {bst.OV2Center_dtdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
    }
    else if ((bset == "rm1")||(bset == "RM1")) {
      RM1BSet bst;
      if ((var == "")||(var == " ")) {bst.OV2Center(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "1")||(var == "R")||(var == "r")) {bst.OV2Center_dR(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "3")||(var == "P")||(var == "p")||(var == "PHI")||(var == "phi")||(var == "Phi")) {bst.OV2Center_dp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "2")||(var == "T")||(var == "t")||(var == "THETA")||(var == "theta")||(var == "Theta")) {bst.OV2Center_dt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "4")||(var == "RR")||(var == "rr")) {bst.OV2Center_dR2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "6")||(var == "PP")||(var == "pp")||(var == "PHIPHI")||(var == "phiphi")||(var == "PhiPhi")) {bst.OV2Center_dp2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "5")||(var == "TT")||(var == "tt")||(var == "THETATHETA")||(var == "thetatheta")||(var == "ThetaTheta")) {bst.OV2Center_dt2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "8")||(var == "RP")||(var == "rp")||(var == "PR")||(var == "pR")||(var == "phiR")||(var == "Rphi")||(var == "PHIR")||(var == "RPHI")) {bst.OV2Center_dRdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "7")||(var == "RT")||(var == "rt")||(var == "TR")||(var == "tR")||(var == "thetaR")||(var == "Rtheta")||(var == "THETAR")||(var == "RTHETA")) {bst.OV2Center_dRdt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "9")||(var == "PT")||(var == "pt")||(var == "TP")||(var == "tp")||(var == "thetaphi")||(var == "phitheta")||(var == "THETAPHI")||(var == "PHITHETA")) {bst.OV2Center_dtdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
    }
    else if ((bset == "am1")||(bset == "AM1")) {
      AM1BSet bst;
      if ((var == "")||(var == " ")) {bst.OV2Center(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "1")||(var == "R")||(var == "r")) {bst.OV2Center_dR(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "3")||(var == "P")||(var == "p")||(var == "PHI")||(var == "phi")||(var == "Phi")) {bst.OV2Center_dp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "2")||(var == "T")||(var == "t")||(var == "THETA")||(var == "theta")||(var == "Theta")) {bst.OV2Center_dt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "4")||(var == "RR")||(var == "rr")) {bst.OV2Center_dR2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "6")||(var == "PP")||(var == "pp")||(var == "PHIPHI")||(var == "phiphi")||(var == "PhiPhi")) {bst.OV2Center_dp2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "5")||(var == "TT")||(var == "tt")||(var == "THETATHETA")||(var == "thetatheta")||(var == "ThetaTheta")) {bst.OV2Center_dt2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "8")||(var == "RP")||(var == "rp")||(var == "PR")||(var == "pR")||(var == "phiR")||(var == "Rphi")||(var == "PHIR")||(var == "RPHI")) {bst.OV2Center_dRdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "7")||(var == "RT")||(var == "rt")||(var == "TR")||(var == "tR")||(var == "thetaR")||(var == "Rtheta")||(var == "THETAR")||(var == "RTHETA")) {bst.OV2Center_dRdt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "9")||(var == "PT")||(var == "pt")||(var == "TP")||(var == "tp")||(var == "thetaphi")||(var == "phitheta")||(var == "THETAPHI")||(var == "PHITHETA")) {bst.OV2Center_dtdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
    }
    else if ((bset == "mndo")||(bset == "MNDO")) {
      MNDOBSet bst;
      if ((var == "")||(var == " ")) {bst.OV2Center(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "1")||(var == "R")||(var == "r")) {bst.OV2Center_dR(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "3")||(var == "P")||(var == "p")||(var == "PHI")||(var == "phi")||(var == "Phi")) {bst.OV2Center_dp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "2")||(var == "T")||(var == "t")||(var == "THETA")||(var == "theta")||(var == "Theta")) {bst.OV2Center_dt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "4")||(var == "RR")||(var == "rr")) {bst.OV2Center_dR2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "6")||(var == "PP")||(var == "pp")||(var == "PHIPHI")||(var == "phiphi")||(var == "PhiPhi")) {bst.OV2Center_dp2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "5")||(var == "TT")||(var == "tt")||(var == "THETATHETA")||(var == "thetatheta")||(var == "ThetaTheta")) {bst.OV2Center_dt2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "8")||(var == "RP")||(var == "rp")||(var == "PR")||(var == "pR")||(var == "phiR")||(var == "Rphi")||(var == "PHIR")||(var == "RPHI")) {bst.OV2Center_dRdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "7")||(var == "RT")||(var == "rt")||(var == "TR")||(var == "tR")||(var == "thetaR")||(var == "Rtheta")||(var == "THETAR")||(var == "RTHETA")) {bst.OV2Center_dRdt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "9")||(var == "PT")||(var == "pt")||(var == "TP")||(var == "tp")||(var == "thetaphi")||(var == "phitheta")||(var == "THETAPHI")||(var == "PHITHETA")) {bst.OV2Center_dtdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
    }
    else if ((bset == "mndod")||(bset == "MNDOd")) {
      MNDOdBSet bst;
      if ((var == "")||(var == " ")) {bst.OV2Center(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "1")||(var == "R")||(var == "r")) {bst.OV2Center_dR(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "3")||(var == "P")||(var == "p")||(var == "PHI")||(var == "phi")||(var == "Phi")) {bst.OV2Center_dp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "2")||(var == "T")||(var == "t")||(var == "THETA")||(var == "theta")||(var == "Theta")) {bst.OV2Center_dt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "4")||(var == "RR")||(var == "rr")) {bst.OV2Center_dR2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "6")||(var == "PP")||(var == "pp")||(var == "PHIPHI")||(var == "phiphi")||(var == "PhiPhi")) {bst.OV2Center_dp2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "5")||(var == "TT")||(var == "tt")||(var == "THETATHETA")||(var == "thetatheta")||(var == "ThetaTheta")) {bst.OV2Center_dt2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "8")||(var == "RP")||(var == "rp")||(var == "PR")||(var == "pR")||(var == "phiR")||(var == "Rphi")||(var == "PHIR")||(var == "RPHI")) {bst.OV2Center_dRdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "7")||(var == "RT")||(var == "rt")||(var == "TR")||(var == "tR")||(var == "thetaR")||(var == "Rtheta")||(var == "THETAR")||(var == "RTHETA")) {bst.OV2Center_dRdt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "9")||(var == "PT")||(var == "pt")||(var == "TP")||(var == "tp")||(var == "thetaphi")||(var == "phitheta")||(var == "THETAPHI")||(var == "PHITHETA")) {bst.OV2Center_dtdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
    }
    else if (bset == "burns") {
      BurnsBSet bst;
      if ((var == "")||(var == " ")) {bst.OV2Center(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "1")||(var == "R")||(var == "r")) {bst.OV2Center_dR(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "3")||(var == "P")||(var == "p")||(var == "PHI")||(var == "phi")||(var == "Phi")) {bst.OV2Center_dp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "2")||(var == "T")||(var == "t")||(var == "THETA")||(var == "theta")||(var == "Theta")) {bst.OV2Center_dt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "4")||(var == "RR")||(var == "rr")) {bst.OV2Center_dR2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "6")||(var == "PP")||(var == "pp")||(var == "PHIPHI")||(var == "phiphi")||(var == "PhiPhi")) {bst.OV2Center_dp2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "5")||(var == "TT")||(var == "tt")||(var == "THETATHETA")||(var == "thetatheta")||(var == "ThetaTheta")) {bst.OV2Center_dt2(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "8")||(var == "RP")||(var == "rp")||(var == "PR")||(var == "pR")||(var == "phiR")||(var == "Rphi")||(var == "PHIR")||(var == "RPHI")) {bst.OV2Center_dRdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "7")||(var == "RT")||(var == "rt")||(var == "TR")||(var == "tR")||(var == "thetaR")||(var == "Rtheta")||(var == "THETAR")||(var == "RTHETA")) {bst.OV2Center_dRdt(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
      else if ((var == "9")||(var == "PT")||(var == "pt")||(var == "TP")||(var == "tp")||(var == "thetaphi")||(var == "phitheta")||(var == "THETAPHI")||(var == "PHITHETA")) {bst.OV2Center_dtdp(sao_dxcd,atmC,atmD,RCD,diffAtom,rotate);}
    }
    else if ((bset == "gfn2")||(bset == "GFN2")) {
      GFN2BSet bst;
      //to implement
    }
    else if ((bset == "gfn2l")||(bset == "GFN2L")) {
      GFN2BSet bst;
      //to implement
    }
  }
  void OrthOverlap() {
    //function to obtain {S}^{-\frac{1}{2}}, where {S} is the overlap matrix
    overlaps = SAO();
    CoeffsA = MatDiag(overlaps);                 //CoeffsA = eigenvalues
    s12 = Identity(NAOs);
    for (size_t idx = 0; idx < NAOs; ++idx) {
      if (CoeffsA[idx] < 0) {std::cout << "WARNING: BSet.hpp: OrthOverlap(): negative eigenvalue!!!!!" << std::endl;}
      s12(idx + 1,idx + 1) = 1.0/sqrt(CoeffsA[idx]);
    }
    xao = overlaps*s12;
    xao *= overlaps.trans();
  }
  matrixE S12() {
    //function to obtain {S}^{\frac{1}{2}}, where {S} is the overlap matrix
    overlaps = SAO();
    CoeffsA = MatDiag(overlaps);                 //CoeffsA = eigenvalues
    s12 = Identity(NAOs);
    for (size_t idx = 0; idx < NAOs; ++idx) {
      if (CoeffsA[idx] < 0) {std::cout << "WARNING: BSet.hpp: S12(): negative eigenvalue!!!!!" << std::endl;}
      s12(idx + 1,idx + 1) = sqrt(CoeffsA[idx]);
    }
    matrixE aux = overlaps*s12;
    aux *= overlaps.trans();
    return aux;
  }
};

#endif //_Basis_Set_
