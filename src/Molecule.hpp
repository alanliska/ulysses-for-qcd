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

#ifndef _Molecule_
#define _Molecule_
#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctype.h>
#include "Global.hpp"
#include "Geometry.hpp"
#include "atoms/AtomPackage.hpp"
#include "parameters/AtomicRadiipar.hpp"
#include "parameters/BondDistancepar.hpp"
#include "ConstantsPackage.hpp"
#include "other/listPackage.hpp"
#include "other/auxPackage.hpp"
#include "other/FileFormats.hpp"
#include "other/MoleculeAnalyser.hpp"
#include "math/MatrixPackage.hpp"
#include "math/RandomPackage.hpp"
#include "math/VectorPackage.hpp"

//description:
//molecule class which is the basis for quantum chemical calculations

class Molecule {
  matrixE geometry;
  std::vector<size_t> atoms;
  std::string geomfile;
  size_t natoms;
  size_t mult;
  int nelectrons;
  int charge;
  int maxcoordinationnumber;               //the maximum coordination number allowed in getting connectivity matrix
  double MWeight;
  double mass;
  std::string symm;
public:
  Molecule(): geometry(1,1) {
    maxcoordinationnumber = 8;
  }
  Molecule(std::string _file, int _chg = 0, size_t _mult = 1, std::string _symm = "0", bool centerstructure = true) {
    charge = _chg;
    geomfile = _file;
    size_t ncharacters = _file.size();
    std::string extension = "";
    int fetch = -1;
    for (size_t idel = ncharacters - 1; idel >= 0; --idel) {
      if (std::string(1,_file[idel]) == ".") {
        fetch = idel;
        break;
      }
    }
    for (size_t idel = fetch + 1; idel < ncharacters; ++idel) {
      extension += _file[idel];
    }
    if (extension == "xyz") {ReadXYZ();}
    else if (extension == "mol2") {ReadMOL2();}
    else if (extension == "sdf") {ReadSDF();}
    /////////else if (extension == "gjf") {ReadGJF();}
    else if (extension == "pdb") {ReadPDB();}
    //if (Linear()) {linearGeom();}
    //else if (Planar()) {planarGeom();}
    masses();
    CountElectrons();
    mult = _mult;
    if ((nelectrons%2 != 0)&&(mult == 1)) {mult = 2;}
    symm = _symm;
    if (centerstructure) {CorrectedGeom();}
    maxcoordinationnumber = 8;
  }
  Molecule(const Molecule & rhs): geometry(rhs.Geometry()) {
    geometry = rhs.Geometry();
    atoms = rhs.Atoms();
    geomfile = rhs.GeomFile();
    natoms = rhs.Natoms();
    mult = rhs.Multiplicity();
    nelectrons = rhs.Nelectrons();
    charge = rhs.Charge();
    MWeight = rhs.MW();
    mass = rhs.m();
    maxcoordinationnumber = 8;
  }
  ~Molecule() {}
//--------------------------------- Counting functions
  int MoleculeCounter(std::vector<size_t> & molvector) {
    if (molvector.size() != natoms) {molvector.resize(natoms);}
    std::vector<bool> assigned(natoms,false);
    matrix<int> neighbours;
    ConnectivityMatrix(neighbours);                 //getting neighbours
    int counter = 1;
    bool neighbourassigned;
    size_t maxnrneighbours = neighbours.cols();
    for (size_t idAtm = 0; idAtm < natoms; ++idAtm) {
      molvector[idAtm] = 0;
    }
    for (size_t idAtm = 0; idAtm < natoms; ++idAtm) {
      if (!assigned[idAtm]) {
        neighbourassigned = false;
        for (size_t idneigh = 0; idneigh < maxnrneighbours; ++idneigh) {
          if (neighbours(idAtm + 1,idneigh + 1) - 1 < 0) {break;}
          if (assigned[neighbours(idAtm + 1,idneigh + 1) - 1]) {
            molvector[idAtm] = molvector[neighbours(idAtm + 1,idneigh + 1) - 1];
            neighbourassigned = true;
            break;
          }
        }
        if (!neighbourassigned) {
          molvector[idAtm] = counter;
          ++counter;
        }
        assigned[idAtm] = true;
      }
    }
    --counter;
    return counter;
  }
//--------------------------------- functions to add molecules
  void set2System(matrixE & geom, std::vector<size_t> & atm, int chrge, int multpl, std::string pg) {
    //function that sets the molecule to external variables
    this->setGeometry(geom);
    this->setAtoms(atm);
    this->calcNatoms();
    this->setCharge(chrge);
    this->setMultiplicity(multpl);
    this->masses();
    this->setPGroup(pg);
  }
  void AddMolecule(Molecule & extramol) {
    matrixE geomextra = extramol.Geometry();
    std::vector<size_t> extraatoms = extramol.Atoms();
    AddMolecule(geomextra,extraatoms,extramol.Charge(),extramol.Multiplicity());
  }
  void AddMolecule(matrixE & extraG, std::vector<size_t> & extraA, int extrachg, size_t extramult) {
    ConcatenateMR(geometry,extraG);
    ConcatenateV(atoms,extraA);
    charge += extrachg;
    mult *= extramult;
    this->calcNatoms();
    this->masses();
  }
//--------------------------------- Geometry Functions
  void Make3D(double maxshift = 0.5, bool complexanalysis = false, double thresholdzero = 1.0e-7) {
    //function that takes a potentially flat molecule (one dimension is zero) and generates distortions on the flat dimension
    bool isflat[3];
    double randomnumber;
    int randomsign;
    isflat[0] = true;
    isflat[1] = true;
    isflat[2] = true;
    //determine whether molecule is flat
    if (complexanalysis) {
      //check for planarity
      maxshift /= 3.0;
      isflat[0] = Planar();
      isflat[1] = isflat[0];
      isflat[2] = isflat[0];
    }
    else {
      //check for zeroed dimensions
      for (size_t idAtm = 0; idAtm < natoms; ++idAtm) {
        if (fabs(geometry(idAtm + 1,1)) > thresholdzero) {isflat[0] = false;}
        if (fabs(geometry(idAtm + 1,2)) > thresholdzero) {isflat[1] = false;}
        if (fabs(geometry(idAtm + 1,3)) > thresholdzero) {isflat[2] = false;}
        if (isflat[0] + isflat[1] + isflat[2] == 0) {break;}
      }
    }
    //remove flatness
    for (size_t idAtm = 0; idAtm < natoms; ++idAtm) {
      randomnumber = fRandom(0.0,1.0);
      randomsign = iRandom01();
      geometry(idAtm + 1,1) += (2*randomsign - 1)*randomnumber*maxshift;
      randomnumber = fRandom(0.0,1.0);
      randomsign = iRandom01();
      geometry(idAtm + 1,2) += (2*randomsign - 1)*randomnumber*maxshift;
      randomnumber = fRandom(0.0,1.0);
      randomsign = iRandom01();
      geometry(idAtm + 1,3) += (2*randomsign - 1)*randomnumber*maxshift;
    }
  }
  void AlignTo(Molecule & reference, matrixE & align, int firstN = -1) {
    //function that aligns this molecule to another one
    matrixE rgeom = reference.Geometry();
    AlignBonA(rgeom,geometry,align,firstN);
  }
  void AlignTo(Molecule & reference) {
    //function that aligns this molecule to another one
    matrixE rgeom = reference.Geometry();
    matrixE align;
    AlignBonA(rgeom,geometry,align,natoms);
  }
  void setAtoms(std::vector<size_t> * _atoms) {
    atoms = *_atoms;
    natoms = atoms.size();
  }
  void setAtoms(std::vector<size_t> & _atoms) {
    atoms = _atoms;
    natoms = atoms.size();
  }
  void setGeometry(matrixE * _geometry) {geometry = *_geometry;}
  void setGeometry(matrixE & _geometry) {geometry = _geometry;}
  std::vector<size_t> Atoms() {return atoms;}
  std::vector<size_t> Atoms() const {return atoms;}
  size_t Atom(size_t pos) {return atoms[pos - 1];}
  size_t Atom(size_t pos) const {return atoms[pos - 1];}
  matrixE Geometry() {return geometry;}
  matrixE Geometry() const {return geometry;}
  void transGeometry(std::vector<double> translation) {
    //function that translates geometry according to vector
    if (translation.size() != 3) {throw std::string("ERROR: Molecule.hpp: Molecule: transGeometry(): translation with vector not in R3");}
    for (size_t iatom = 0; iatom < natoms; ++iatom) {
      for (size_t icoord = 0; icoord < 3; ++icoord) {
        geometry(iatom + 1,icoord + 1) += translation[icoord];
      }
    }
  }
  void rotGeometry(size_t axis, double theta) {
    //function that rotates geometry according to certain Given's transformation matrix
    std::vector<double> cmass = this->CM();
    //recenter molecule
    for (size_t idatm = 0; idatm < geometry.rows(); ++idatm) {
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        geometry(idatm + 1,idcoord + 1) -= cmass[idcoord];
      }
    }
    matrixE rotation;
    RotateGeomMatrix(axis,theta,geometry,rotation);
    //shift molecule back to its position
    for (size_t idatm = 0; idatm < geometry.rows(); ++idatm) {
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        geometry(idatm + 1,idcoord + 1) += cmass[idcoord];
      }
    }
  }
  void rotGeometry2(size_t axis, double theta, matrixE & rotation) {
    //like function before, but quick and "dirty" as it assumes structure is centered
    RotateGeomMatrix(axis,theta,geometry,rotation);
  }
  std::string GeomFile() {return geomfile;}
  std::string GeomFile() const {return geomfile;}
  void ReadXYZ() {ReadXYZFormat(geomfile,natoms,geometry,atoms);}
  void ReadMOL2() {ReadMOL2Format(geomfile,natoms,geometry,atoms,100000000);}
  void ReadSDF() {ReadSDFFormat(geomfile,natoms,geometry,atoms,charge,100000000);}
  void ReadPDB() {ReadPDBFormat(geomfile,natoms,geometry,atoms,100000000);}
  void WriteXYZ(std::string name, int counter = -1, int prc = 7) {
    //function to write geometries into a given xyz file
    WriteXYZFormat(name,atoms,geometry,counter,prc);
  }
  std::vector<double> CM() {
    //function to calculate the center of mass
    std::vector<double> cm = CenterOfMass(geometry,atoms);
    return cm;
  }
  std::vector<double> CMg(matrixE geom) {
    //just like the previous function but it reads the geometry from somewhere else
    //the vector of atoms is the same
    std::vector<double> cm = CenterOfMass(geom,atoms);
    return cm;
  }
  void CorrectedGeom() {
    //function to correct a geometry using the center of mass, i.e., that puts geometry in the center of the referencial
    std::vector<double> cm = CM();
    for (size_t idx = 0; idx < natoms; ++idx) {
      geometry(idx+1,1) -= cm[0];
      geometry(idx+1,2) -= cm[1];
      geometry(idx+1,3) -= cm[2];
    }
  }
  bool CompGeom(matrixE MA, matrixE MB, double tolerance) {
    //function comparing 2 geometries according to a given tolerance
    bool same = true;
    std::vector<bool> symmetric(natoms,false);
    std::vector<bool> comp(3,false);
    double diff;
    double MAX1;
    double MAX2;
    double aux;
    for (size_t irw = 0; irw < MA.rows(); ++irw) {
      for (size_t irw2 = 0; irw2 < MA.rows(); ++irw2) {
        for (size_t icl = 0; icl < MA.cols(); ++icl) {
          comp[icl] = false;
          diff = fabs(MB(irw2+1,icl+1) - MA(irw+1,icl+1));
          if (diff < tolerance) {
            if (atoms[irw2] == atoms[irw]) {
              comp[icl] = true;
            }
          }
        }
        if ((comp[0] == true)&&(comp[1] == true)&&(comp[2] == true)) {
          symmetric[irw] = true;
          break;
        }
      }
      if (symmetric[irw] == false) {
        same = false;
        break;
      }
    }
    return same;
  }
  bool Linear(double tolerance = pi/180) {
    //function to check whether a molecule is Linear
    bool linear = true;
    if (natoms > 2) {
      //can just choose 2 atoms and check whether the angle is close to 0 or 180. Search by positive then
      double angle;
      linear = false;
      for (size_t atms = 3; atms < (natoms+1); ++atms) {
        angle = calcAngles(1,2,atms);
        if ((fabs(angle) < tolerance)||(fabs(angle - pi) < tolerance)) {linear = true;}
        else {
          linear = false;
          break;
        }
      }
    }
    return linear;
  }
  void linearGeom() {
    //function that rebuilds the geometry to be linear
    matrixE xgeom(natoms,3);
    if (natoms == 2) {
      xgeom(1,1) = -1e10*distance(1,2)/2;
      xgeom(1,2) = 0.0;
      xgeom(1,3) = 0.0;
      xgeom(2,1) = 1e10*distance(1,2)/2;
      xgeom(2,2) = 0.0;
      xgeom(2,3) = 0.0;
    }
    else {
      std::vector<std::vector<size_t> > connect = Connections();
      size_t begin;
      for (size_t start = 0; start < connect.size(); ++start) {
        if (connect[start].size() == 1) {
          begin = start+1;
          break;
        }
        else {
          size_t zeroes = 0;
          for (size_t idx = 0; idx < connect[start].size(); ++idx) {
            if (connect[start][idx] == 0) {++zeroes;}
          }
          if (zeroes == (connect[start].size() - 1)) {
            begin = start+1;
            break;
          }
        }
      }
      std::cout << "begin = " << begin << std::endl;
      xgeom(begin,1) = 0.0;
      xgeom(begin,2) = 0.0;
      xgeom(begin,3) = 0.0;
      size_t connected = connect[begin-1][0];
      size_t aux;
      std::cout << "B" << std::endl;
      for (size_t idx = 0; idx < natoms; ++idx) {
        xgeom(connected,1) = 1e10*distance(begin,connected)+xgeom(begin,1);
        xgeom(connected,2) = 0.0;
        xgeom(connected,3) = 0.0;
        aux = begin;
        begin = connected;
        for (size_t idx2 = 0; idx2 < connect[begin-1].size(); ++idx2) {
          if ((connect[begin-1][idx2] != aux)&&(connect[begin-1][idx2] != 0)) {
            connected = connect[begin-1][idx2];
            break;
          }
        }
      }
    }
    geometry = xgeom;
  }
  bool Planar(double tolerance = pi/180) {
    //function to check whether a molecule is planar
    bool planar = true;
    if (natoms == 4) {
      double angle = calcDihedrals(1,2,3,4);
      if (fabs(angle) < 0.5*pi) {
        if (fabs(angle) < tolerance) {planar = true;}
        else {planar = false;}
      }
      else {
        if (fabs(fabs(angle) - pi) < tolerance) {planar = true;}
        else {planar = false;}
      }
    }
    else if (natoms == 5) {
      double angle = calcDihedrals(1,2,3,4);
      if (fabs(angle) < 0.5*pi) {
        if (fabs(angle) < tolerance) {
          //planar = true, lets check the other angles
          angle = calcDihedrals(2,3,4,5);
          if (fabs(angle) < 0.5*pi) {
            if (fabs(angle) < tolerance) {planar = true;}
            else {planar = false;}
          }
          else {
            if (fabs(fabs(angle) - pi) < tolerance) {planar = true;}
            else {planar = false;}
          }
        }
        else {planar = false;}
      }
      else {
        if (fabs(fabs(angle) - pi) < tolerance) {
          //planar = true, lets check the other angles
          angle = calcDihedrals(2,3,4,5);
          if (fabs(angle) < 0.5*pi) {
            if (fabs(angle) < tolerance) {planar = true;}
            else {planar = false;}
          }
          else {
            if (fabs(fabs(angle) - pi) < tolerance) {planar = true;}
            else {planar = false;}
          }
        }
        else {planar = false;}
      }
    }
    else if (natoms > 5) {
      //loop over atoms and create lines of atoms connected
      double angle;
      std::vector<std::vector<size_t> > connect = Connections();
      size_t atm2;
      size_t atm3;
      size_t atm4;
      size_t connsz1;
      size_t connsz2;
      size_t connsz3;
      planar = true;
      for (size_t atm1 = 0; atm1 < natoms; ++atm1) {
        connsz1 = connect[atm1].size();
        if (connsz1 > 4) {
          planar = false;
          break;
        }
        for (size_t idx2 = 0; idx2 < connsz1; ++idx2) {
          atm2 = connect[atm1][idx2];
          connsz2 = connect[atm2 - 1].size();
          if (atm2 == 0) {continue;}
          for (size_t idx3 = 0; idx3 < connsz2; ++idx3) {
            atm3 = connect[atm2 - 1][idx3];
            if (((atm3 - 1) == atm1)||(atm3 == 0)) {continue;}
            connsz3 = connect[atm3 - 1].size();
            for (size_t idx4 = 0; idx4 < connsz3; ++idx4) {
              atm4 = connect[atm3 - 1][idx4];
              if (((atm4 - 1) == atm1)||(atm4 == atm2)||(atm4 == 0)) {continue;}
              angle = calcDihedrals(atm1 + 1,atm2,atm3,atm4);
              if (fabs(angle) < 0.5*pi) {
                if (fabs(angle) < tolerance) {planar = true;}
                else {
                  planar = false;
                  angle = calcDihedrals(atm1 + 1,atm3,atm4,atm2);
                  if (fabs(angle) < tolerance) {
                    planar = true;
                    continue;
                  }
                  break;
                }
              }
              else {
                if (fabs(fabs(angle) - pi) < tolerance) {planar = true;}
                else {
                  planar = false;
                  angle = calcDihedrals(atm1 + 1,atm3,atm4,atm2);
                  if ((fabs(angle) > 0.5*pi)&&(fabs(fabs(angle) - pi) < tolerance)) {
                    planar = true;
                    continue;
                  }
                  else if ((fabs(angle) < 0.5*pi)&&(fabs(angle) < tolerance)) {
                    planar = true;
                    continue;
                  }
                  break;
                }
              }
            }
            if (planar == false) {break;}
          }
          if (planar == false) {break;}
        }
        if (planar == false) {break;}
      }
    }
    return planar;
  }
  matrixE makePlanar() {
    //molecule is planar, rebuild geometry so that z component is zero for all atoms
    //this is not easily done with a couple of rotations because the molecule may have funny orientation
    //furthermore, I need to be careful with the angles, since they do not give me an absolute direction!
    matrixE ngeom = geometry;
    //defining coordinates of first 2 atoms
    for (size_t idn = 1; idn < 4; ++idn) {
      for (size_t idn2 = 1; idn2 < natoms + 1; ++idn2) {
        ngeom(idn2,idn) = 0.0;
      }
    }
    ngeom(2,1) = 1e10*distance(1,2);
    std::vector<std::vector<size_t> > connmat = Connections();
    //this is a non-elegant workaround to get correct coordinates for all atoms of the planar geometry
    size_t con1 = 0;
    size_t con2 = 0;
    for (size_t idn1 = 0; idn1 < connmat[0].size(); ++idn1) {
      if (connmat[0][idn1] > con1) {
        //check whether the atom connected to the first atom has only one connection
        size_t pos = connmat[0][idn1];
        size_t count0 = 0;
        if (connmat[pos-1].size() == 1) {con1 = connmat[0][idn1];}
        for (size_t idn2 = 0; idn2 < connmat[pos-1].size(); ++idn2) {
          if (connmat[pos-1][idn2] == 0) {++count0;}
        }
        if (connmat[pos-1].size() == (count0+1)) {con1 = connmat[0][idn1];}
      }
    }
    for (size_t idn1 = 0; idn1 < connmat[1].size(); ++idn1) {
      if (connmat[1][idn1] > con2) {
        //check whether the atom connected to the second atom has only one connection
        size_t pos = connmat[1][idn1];
        size_t count0 = 0;
        if (connmat[pos-1].size() == 1) {con2 = connmat[1][idn1];}
        for (size_t idn2 = 0; idn2 < connmat[pos-1].size(); ++idn2) {
          if (connmat[pos-1][idn2] == 0) {++count0;}
        }
        if (connmat[pos-1].size() == (count0+1)) {con2 = connmat[1][idn1];}
      }
    }
    //defining rest of coordinates
    int sign;
    double theta;
    double aux;
    for (size_t idn = 3; idn < (natoms + 1); ++idn) {
      theta = calcAngles(1,2,idn);
      //must correct theta for atoms connected to atoms 1 and 2
      aux = 1e10*distance(1,idn);
      ngeom(idn,1) = aux*cos(theta);
      if ((idn == con1)||(idn == con2)) {sign = -1;}
      else {sign = 1;}
      ngeom(idn,2) = sign*aux*sin(theta);
      //making sure there are no colisions
      for (size_t idn2 = 1; idn2 < (idn); ++idn2) {
        aux = (ngeom(idn,2) - ngeom(idn2,2))*(ngeom(idn,2) - ngeom(idn2,2));
        aux += (ngeom(idn,1) - ngeom(idn2,1))*(ngeom(idn,1) - ngeom(idn2,1));
        if (aux < 1.0) {
          ngeom(idn,2) *= -1;
          break;
        }
      }
      ngeom(idn,3) = 0.0;
    }
    //ensuring the new planar geometry is correct
    std::vector<std::vector<size_t> > connmat2 = ConnectionsG(ngeom);
    bool comparison = true;
    int position;
    size_t missing;
    if (connmat.size() != connmat2.size()) {std::cout << "something very wrong... connmat and connmat2 do not have the same size" << std::endl;}
    else {
      for (size_t idb = 0; idb < connmat.size(); ++idb) {
        comparison = true;
        //look for discrepancies
        for (size_t idc = 0; idc < connmat[idb].size(); ++idc) {
          if (connmat[idb][idc] != connmat2[idb][idc]) {
            comparison = false;
            position = idc;
            break;
          }
        }
        size_t cnt1 = 0;
        size_t cnt2 = 0;
        //where is an atom missing?
        for (size_t idc = 0; idc < connmat[idb].size(); ++idc) {
          if (connmat[idb][idc] != 0) {++cnt1;}
          if (connmat2[idb][idc] != 0) {++cnt2;}
        }
        if (cnt1 > cnt2) {missing = 2;}
        else{missing = 1;}
        //do something about it
        if (!comparison) {
          //I assume that the misplaced atoms all have single connection
          if (missing == 2) {
            //in this situation there is an atom to add to list
            double distr = 1e10*distance(1,connmat[idb][position]);
            size_t pos;
            if (position == 0) {pos = 1;}
            else {pos = 0;}
            //only do this if there is actually an atom to make an angle with
            if (connmat[idb][pos] == 0) {continue;}
            if (connmat[idb].size() == 1) {continue;}
            double angtheta = calcAngles(1,2,connmat[idb][position]);
            ngeom(connmat[idb][position],1) = distr*cos(angtheta);
            ngeom(connmat[idb][position],2) = distr*sin(angtheta);
          }
        }
      }
    }
    std::vector<double> cm = CMg(ngeom);
    for (size_t idx = 0; idx < natoms; ++idx) {
      ngeom(idx+1,1) -= cm[0];
      ngeom(idx+1,2) -= cm[1];
    }
    return ngeom;
  }
  void Cart2Zmat(matrixE & zmat, std::vector<size_t> & atomorder) {
    //function converting a geometry in cartesian coordinates into Z-matrix
    size_t atmordersz = atomorder.size();
    if (atmordersz != natoms) {               //by allowing this, I can focus the Z-matrix on some specific atoms only, and the rest is whatever
      bool present;
      //loop over atoms
      for (size_t idatm = 1; idatm < natoms + 1; ++idatm) {
        present = false;
        //check whether atom position already in atomorder vector
        for (size_t idpos = 0; idpos < atmordersz; ++idpos) {
          if (atomorder.at(idpos) == idatm) {
            present = true;
            break;
          }
        }
        if (!present) {atomorder.push_back(idatm);}
      }
    }
    zmat.resize(natoms,6);                 //reshape it; note that atomic numbers are given separately in an array
    zmat.zero();
    zmat(2,1) = 1.0;
    std::vector<double> dij = this->calcUDistances();
    zmat(2,2) = dij[LTriangMat2Array_nodiag(atomorder.at(0),atomorder.at(1),natoms)];
    if (natoms > 2) {
      zmat(3,1) = 2.0;
      zmat(3,2) = dij[LTriangMat2Array_nodiag(atomorder.at(1),atomorder.at(2),natoms)];
      zmat(3,3) = 1.0;
      zmat(3,4) = calcAngles(atomorder.at(1),atomorder.at(2),atomorder.at(0));
      if (natoms > 3) {
        for (size_t idatm = 3; idatm < natoms; ++idatm) {
          zmat(idatm + 1,1) = double(idatm);
          zmat(idatm + 1,2) = dij[LTriangMat2Array_nodiag(atomorder.at(idatm),atomorder.at(idatm - 1),natoms)];
          zmat(idatm + 1,3) = double(idatm - 1);
          zmat(idatm + 1,4) = calcAngles(atomorder.at(idatm - 1),atomorder.at(idatm),atomorder.at(idatm - 2));
          zmat(idatm + 1,5) = double(idatm - 2);
          zmat(idatm + 1,6) = calcDihedrals(atomorder.at(idatm),atomorder.at(idatm - 1),atomorder.at(idatm - 2),atomorder.at(idatm - 3));
        }
      }
    }
  }
  void Zmat2Cart(matrixE & zmat) {
    //function converting a geometry in Z-matrix into cartesian coordinates; based on TMPChem's python code
    natoms = zmat.rows();
    double rij;
    double theta;
    double factor;
    size_t atmR;
    size_t atmA;
    size_t atmD;
    std::vector<double> u21(3,0.0);
    std::vector<double> u23(3,0.0);
    std::vector<double> cp(3,0.0);
    std::vector<double> bv(3,0.0);
    matrixE cart(natoms,3);
    cart(2,3) = zmat(2,2);
    if (natoms > 1) {
      if (natoms > 2) {
        factor = 1.0;
        if (int(zmat(3,1)) == 2) {factor = -1.0;}      //if atom 3 is bonded to atom 2
        rij = zmat(3,2);
        theta = zmat(3,4);
        cart(3,2) = rij*sin(theta);
        cart(3,3) = cart(2,3) + factor*rij*cos(theta);
        for (size_t idatm = 4; idatm < natoms + 1; ++idatm) {
          //get coordinates of atoms defining the position; in my z-matrix construction this step is actually not needed; I keep this so that I am as general as possible
          atmR = size_t(zmat(idatm,1));
          atmA = size_t(zmat(idatm,3));
          atmD = size_t(zmat(idatm,5));
          GetUnitVectors(u21,u23,cart,atmA,atmR,atmD);
          GetLocalAxes(u21,u23,cp);
          GetBondVector(bv,zmat,idatm);
          //get atomic coordinates
          cart(idatm,1) = cart(atmR,1) + dotProd(bv,cp);
          cart(idatm,2) = cart(atmR,2) + dotProd(bv,u23);
          cart(idatm,3) = cart(atmR,3) + dotProd(bv,u21);
        }
      }
    }
    geometry = cart;
  }
//--------------------------------- Connectivity Functions
  int MaxCoordinationNumber() {return maxcoordinationnumber;}
  void setMaxCoordinationNumber(int _maxcdn) {maxcoordinationnumber = _maxcdn;}
  void ConnectivityMatrix(matrix<int> & ConnMat, double factor = 1.3) {
    //function that calculates connectivity between atoms as a matrix
    //factor is a scaling factor to apply on vdW radii
    double rAB2;
    double aux;
    double rABcov;
    std::vector<int> count(natoms,1);
    if ((ConnMat.rows() != natoms)&&(ConnMat.cols() != maxcoordinationnumber)) {ConnMat.resize(natoms,maxcoordinationnumber);}
    for (size_t idAtm = 0; idAtm < natoms; ++idAtm) {
      for (size_t idpos = 1; idpos < maxcoordinationnumber + 1; ++idpos) {
        ConnMat(idAtm + 1,idpos) = 0;
      }
      for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
        rAB2 = 0.0;
        for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
          aux = geometry(idAtm + 1,idcoord) - geometry(idBtm + 1,idcoord);
          rAB2 += aux*aux;
        }
        rABcov = factor*(AtmRadii(atoms[idAtm]) + AtmRadii(atoms[idBtm]));
        if (rAB2 < rABcov*rABcov) {                     //neighbour found
          ConnMat(idAtm + 1,count[idAtm]) = idBtm + 1;
          ++count[idAtm];
          ConnMat(idBtm + 1,count[idBtm]) = idAtm + 1;
          ++count[idBtm];
        }
      }
    }
  }
  bool CheckBondDistances(double minbondlen = 1.0, bool quick = false, double scalefactor = 1.3) {
    //function to check bond distances; basically whether there is a large overlap between atoms
    //boolean quick is to determine which algorithm to take
    //quick -> check only immediate neighbours
    //!quick -> loop over all atoms
    double RAB;
    bool problem = false;
    if (quick) {
      matrix<int> neighbours;
      ConnectivityMatrix(neighbours,scalefactor);
      size_t totalnghb = neighbours.cols();
      for (size_t idatm = 0; idatm < natoms; ++idatm) {
        for (size_t idneigh = 0; idneigh < totalnghb; ++idneigh) {
          if (neighbours(idatm + 1,idneigh + 1) - 1 < 0) {break;}
          RAB = Distance(idatm + 1,neighbours(idatm + 1,idneigh + 1),geometry);
          if (RAB < minbondlen) {
            problem = true;
            std::cout << "potential abnormal bond: " << idatm + 1 << "," << neighbours(idatm + 1,idneigh + 1) << "\n";
          }
        }
      }
    }
    else {
      for (size_t idatm = 0; idatm < natoms; ++idatm) {
        for (size_t idbtm = idatm + 1; idbtm < natoms; ++idbtm) {
          RAB = Distance(idatm + 1,idbtm + 1,geometry);
          if (RAB < minbondlen) {
            problem = true;
            std::cout << "potential abnormal bond: " << idatm + 1 << "," << idbtm + 1 << "\n";
          }
        }
      }
    }
    return problem;
  }
  void CheckSystem(double scalefactor = 1.25, double mindist = 1.0, bool quick = false) {
    //function that checks for particular chemical functionalities and potential inconsistencies
    matrix<int> neighbours;
    ConnectivityMatrix(neighbours,scalefactor);
    std::vector<std::string> functionalgroups = CheckGroupsByNeighborhood(atoms,neighbours,geometry);
    std::cout << "functional groups in system" << std::endl;
    for (size_t idx = 0; idx < natoms; ++idx) {
      std::cout << functionalgroups[idx] << std::endl;
    }
    std::cout << "----------------" << std::endl;
    std::cout << "running check on bond distances " << std::endl;
    bool problem = CheckBondDistances(mindist,quick,scalefactor);
    std::cout << "problems found: " << problem << std::endl;
    int totalcharge = TotalChargeEstimator(functionalgroups,atoms,neighbours);
    std::cout << "estimated total charge: " << totalcharge << std::endl;
  }
  void AnalyseNeighbours(double scalefactor = 1.3) {
    //function that analyses a system by checking it for correct number of neighbours
    matrix<int> neighbours;
    ConnectivityMatrix(neighbours,scalefactor);
    size_t totalnghb = neighbours.cols();
    size_t atomA;
    size_t nneigh;
    size_t auxie;
    //now that I have the connectivity matrix, lets check for consistency in neighbours
    for (size_t idatm = 0; idatm < natoms; ++idatm) {
      atomA = atoms[idatm];
      nneigh = 0;
      for (size_t idneigh = 0; idneigh < totalnghb; ++idneigh) {
        if (neighbours(idatm + 1,idneigh + 1) > 0) {++nneigh;}
        else {break;}
      }
      switch (atomA) {
        case 1:
          if (nneigh > 1) {std::cout << "abnormal atom H (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 2:
          break;
        case 3:
          if (nneigh > 1) {std::cout << "abnormal atom Li (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 4:
          break;
        case 5:
          if (nneigh < 2) {std::cout << "abnormal atom B (too few bonds): " << idatm + 1 << "\n";}
          else if (nneigh > 4) {std::cout << "abnormal atom B (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 6:
          if (nneigh < 2) {std::cout << "abnormal atom C (too few bonds): " << idatm + 1 << "\n";}
          else if (nneigh > 4) {std::cout << "abnormal atom C (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 7:
          if (nneigh > 4) {std::cout << "abnormal atom N (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 8:
          if (nneigh > 3) {std::cout << "abnormal atom O (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 9:
          if (nneigh > 1) {std::cout << "abnormal atom F (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 10:
          break;
        case 11:
          if (nneigh > 1) {std::cout << "abnormal atom Na (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 12:
          if ((nneigh != 2)&&(nneigh != 4)) {std::cout << "abnormal atom Mg: " << idatm + 1 << "\n";}
          break;
        case 13:
          if (nneigh != 3) {std::cout << "abnormal atom Al: " << idatm + 1 << "\n";}
          break;
        case 14:
          break;
        case 15:
          if ((nneigh != 3)&&(nneigh != 4)&&(nneigh != 5)) {std::cout << "abnormal atom P: " << idatm + 1 << "\n";}
          if (nneigh == 3) {
            //count oxygens
            auxie = 0;
            for (size_t idneigh = 0; idneigh < totalnghb; ++idneigh) {
              if (atoms[neighbours(idatm + 1,idneigh + 1) - 1] == 8) {++auxie;}
              else {break;}
            }
            if (auxie == 3) {std::cout << "potential problem on atom P: " << idatm + 1 << ". Coordination number 3 to 3 oxygens -> uncommon.\n";}
          }
          break;
        case 16:
          if (nneigh > 4) {std::cout << "abnormal atom S: " << idatm + 1 << "\n";}
          break;
        case 17:
          break;
        case 18:
          break;
        case 19:
          if (nneigh > 1) {std::cout << "abnormal atom K (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 20:
          break;
        case 21:
          break;
        case 22:
          break;
        case 23:
          break;
        case 24:
          break;
        case 25:
          break;
        case 26:
          break;
        case 27:
          break;
        case 28:
          break;
        case 29:
          break;
        case 30:
          break;
        case 31:
          break;
        case 32:
          break;
        case 33:
          break;
        case 34:
          break;
        case 35:
          break;
        case 36:
          break;
        case 37:
          if (nneigh > 1) {std::cout << "abnormal atom Rb (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 38:
          break;
        case 39:
          break;
        case 40:
          break;
        case 41:
          break;
        case 42:
          break;
        case 43:
          break;
        case 44:
          break;
        case 45:
          break;
        case 46:
          if ((nneigh != 4)&&(nneigh != 6)) {std::cout << "abnormal atom Pd: " << idatm + 1 << "\n";}
          break;
        case 47:
          break;
        case 48:
          break;
        case 49:
          break;
        case 50:
          break;
        case 51:
          break;
        case 52:
          break;
        case 53:
          break;
        case 54:
          break;
        case 55:
          if (nneigh > 1) {std::cout << "abnormal atom Cs (too many bonds): " << idatm + 1 << "\n";}
          break;
        case 56:
          break;
        case 57:
          break;
        case 58:
          break;
        case 59:
          break;
        case 60:
          break;
        case 61:
          break;
        case 62:
          break;
        case 63:
          break;
        case 64:
          break;
        case 65:
          break;
        case 66:
          break;
        case 67:
          break;
        case 68:
          break;
        case 69:
          break;
        case 70:
          break;
        case 71:
          break;
        case 72:
          break;
        case 73:
          break;
        case 74:
          break;
        case 75:
          break;
        case 76:
          break;
        case 77:
          break;
        case 78:
          if ((nneigh != 4)&&(nneigh != 6)) {std::cout << "abnormal atom Pt: " << idatm + 1 << "\n";}
          break;
        case 79:
          break;
        case 80:
          break;
        case 81:
          break;
        case 82:
          break;
        case 83:
          break;
        case 84:
          break;
        case 85:
          break;
        case 86:
          break;
        case 87:
          break;
        case 88:
          break;
        case 89:
          break;
        case 90:
          break;
        case 91:
          break;
        case 92:
          break;
        case 93:
          break;
        case 94:
          break;
      }
    }
  }
  std::vector<std::vector<size_t> > Connections() {
    //function that returns a vector of vectors
    //the outer vector contains the atoms as they appear in the geometry
    //the inner vector contains the numbers of the atoms a given atom is connected to
    matrixE connectivity = calcDistances();
    std::vector<std::vector<size_t> > connections;
    std::vector<double> dist;
    std::vector<size_t> conns;
    size_t _octet;
    double distmol;                     //to write the distance between atoms in molecule
    double maxbonddist;                 //to write the maximum expected bond distance between atoms
    //loop over all atoms
    for (size_t idx = 0; idx < natoms; ++idx) {
      //establishing the maximum number of bonds per atom
      _octet = size_t(octet(atoms[idx]));
      //std::cout << "octet atom " << atoms[idx] << " " << octet(atoms[idx]) << std::endl;
      //std::cout << "Atom " << idx+1 << " octet " << _octet << std::endl;
      dist.resize(_octet);
      conns.resize(_octet);
      //setting an upper bound for first search. Distances are in Angstrom
      for (size_t idx2 = 0; idx2 < _octet; ++idx2) {
        dist[idx2] = 100.0;          //do not really know anymore why I set this upper limit, since 
        conns[idx2] = 0;
      }
      //loop over other atoms to check which are connected to atom idx+1
      for (size_t idx2 = 0; idx2 < natoms; ++idx2) {
        if (idx2 != idx) {
          distmol = connectivity(idx+1,idx2+1)/dist_si2Angstrom;
          maxbonddist = MaxBondDistance(atoms[idx],atoms[idx2],1.30);
          //std::cout << "distance between " << AtomNr2Symbol(atoms[idx]) << " and " << AtomNr2Symbol(atoms[idx2]) << " is ";
          //std::cout << distmol << "; max dist is " << maxbonddist << std::endl;
          //std::cout << _octet << std::endl;
          if ((distmol < dist[_octet-1])&&(distmol < maxbonddist)) {
            for (size_t idx3 = 0; idx3 < _octet; ++idx3 ) {
              if (distmol <= dist[idx3]) {
                size_t aux;
                for (size_t idx4 = 0; idx4 < _octet; ++idx4) {
                  aux = _octet - idx4 - 1;
                  if (aux > idx3) {
                    dist[aux] = dist[aux-1];
                    conns[aux] = conns[aux-1];
                  }
                }
                dist[idx3] = distmol;
                conns[idx3] = idx2+1;
                break;
              }
            }
          }
        }
      }
      connections.push_back(conns);
    }
    return connections;
  }
  std::vector<std::vector<size_t> > ConnectionsG(matrixE geom) {
    //the general Connections function. Just like the previous function but also accepts another geometry.
    matrixE connectivity = calcDistancesG(geom);
    std::vector<std::vector<size_t> > connections;
    std::vector<double> dist;
    std::vector<size_t> conns;
    size_t _octet;
    double distmol;                     //to write the distance between atoms in molecule
    double maxbonddist;                 //to write the maximum expected bond distance between atoms.
    //loop over all atoms
    for (size_t idx = 0; idx < natoms; ++idx) {
      //establishing the maximum number of bonds per atom
      _octet = size_t(octet(atoms[idx]));
      //std::cout << "Atom " << idx+1 << " octet " << _octet << std::endl;
      dist.resize(_octet);
      conns.resize(_octet);
      //setting an upper bound for first search. Distances are in meters
      for (size_t idx2 = 0; idx2 < _octet; ++idx2) {
        dist[idx2] = 1.0;
        conns[idx2] = 0.0;
      }
      //loop over other atoms to check which are connected to atom idx+1
      for (size_t idx2 = 0; idx2 < natoms; ++idx2) {
        if (idx2 != idx) {
          distmol = connectivity(idx+1,idx2+1);
          maxbonddist = MaxBondDistance(atoms[idx],atoms[idx2],1.25);
          //std::cout << "distance between " << AtomNr2Symbol(atoms[idx]) << " and " << AtomNr2Symbol(atoms[idx2]) << " is ";
          //std::cout << distmol << "; max dist is " << maxbonddist << std::endl;
          //std::cout << _octet << std::endl;
          if ((distmol < dist[_octet-1])&&(distmol < maxbonddist)) {
            for (size_t idx3 = 0; idx3 < _octet; ++idx3 ) {
              if (distmol <= dist[idx3]) {
                size_t aux;
                for (size_t idx4 = 0; idx4 < _octet; ++idx4) {
                  aux = _octet - idx4 - 1;
                  if (aux > idx3) {
                    dist[aux] = dist[aux-1];
                    conns[aux] = conns[aux-1];
                  }
                }
                dist[idx3] = distmol;
                conns[idx3] = idx2+1;
                break;
              }
            }
          }
        }
      }
      connections.push_back(conns);
    }
    return connections;
  }
//--------------------------------- distance functions
  double MaxBondDistance(int atom1, int atom2, double fac = 1.2) {
    //function to determine the maximum reasonable bond distance between 2 atoms
    double rad1 = RadiusPsi4(size_t(atom1));
    double rad2 = RadiusPsi4(size_t(atom2));
    return fac*(rad1 + rad2);
  }
  std::vector<double> calcUDistances() {
    //function that determines all bond distances and stores them in an array
    std::vector<double> distances(natoms*(natoms - 1)/2,0.0);
    double distance;
    double distance2;
    for (size_t idatm1 = 1; idatm1 < natoms; ++idatm1) {
      for (size_t idatm2 = 0; idatm2 < idatm1; ++idatm2) {
        distance2 = (geometry(idatm1 + 1,1)-geometry(idatm2 + 1,1))*(geometry(idatm1 + 1,1)-geometry(idatm2 + 1,1));
        distance2 += (geometry(idatm1 + 1,2)-geometry(idatm2 + 1,2))*(geometry(idatm1 + 1,2)-geometry(idatm2 + 1,2));
        distance2 += (geometry(idatm1 + 1,3)-geometry(idatm2 + 1,3))*(geometry(idatm1 + 1,3)-geometry(idatm2 + 1,3));
        distance = sqrt(distance2);
        distances[LTriangMat2Array_nodiag(idatm1 + 1,idatm2 + 1,natoms)] = distance;
      }
    }
    return distances;
  }
  matrixE calcDistances() {
    //function that determines a matrix of all bond distances
    matrixE distances(natoms,natoms);
    double distance;
    double distance2;
    for (size_t idx = 1; idx < natoms; ++idx) {
      for (size_t idx2 = 0; idx2 < idx; ++idx2) {
        distance2 = (geometry(idx+1,1)-geometry(idx2+1,1))*(geometry(idx+1,1)-geometry(idx2+1,1));
        distance2 += (geometry(idx+1,2)-geometry(idx2+1,2))*(geometry(idx+1,2)-geometry(idx2+1,2));
        distance2 += (geometry(idx+1,3)-geometry(idx2+1,3))*(geometry(idx+1,3)-geometry(idx2+1,3));
        distance = sqrt(distance2);
        distances(idx+1,idx2+1) = distance*1e-10;
        distances(idx2+1,idx+1) = distance*1e-10;
      }
    }
    return distances;
  }
  matrixE calcDistancesG(matrixE geom) {
    //generalization of the previous function to other geometries
    matrixE distances(natoms,natoms);
    double distance;
    double distance2;
    for (size_t idx = 1; idx < natoms; ++idx) {
      for (size_t idx2 = 0; idx2 < idx; ++idx2) {
        distance2 = (geom(idx+1,1)-geom(idx2+1,1))*(geom(idx+1,1)-geom(idx2+1,1));
        distance2 += (geom(idx+1,2)-geom(idx2+1,2))*(geom(idx+1,2)-geom(idx2+1,2));
        distance2 += (geom(idx+1,3)-geom(idx2+1,3))*(geom(idx+1,3)-geom(idx2+1,3));
        distance = sqrt(distance2);
        distances(idx+1,idx2+1) = distance*1e-10;
        distances(idx2+1,idx+1) = distance*1e-10;
      }
    }
    return distances;
  }
  double Adistance(size_t el1, size_t el2) {
    //function to calculate the distance between 2 atoms (in Angstrom)
    double _dist = (geometry(el1,1)-geometry(el2,1))*(geometry(el1,1)-geometry(el2,1));
    _dist += (geometry(el1,2)-geometry(el2,2))*(geometry(el1,2)-geometry(el2,2));
    _dist += (geometry(el1,3)-geometry(el2,3))*(geometry(el1,3)-geometry(el2,3));
    return sqrt(_dist);
  }
  double AUdistance(size_t el1, size_t el2) {
    //function to calculate the distance between 2 atoms (in Bohr)
    return Adistance(el1,el2)/dist_Angstrom2au;
  }
  double distance(size_t el1, size_t el2) {
    //function to calculate the distance between 2 atoms (in SI)
    return dist_si2Angstrom*Adistance(el1,el2);
  }
  double distanceG(matrixE & geom, size_t el1, size_t el2) {
    //like the previous function but generallized to any geometry
    double _dist = (geom(el1,1)-geom(el2,1))*(geom(el1,1)-geom(el2,1));
    _dist += (geom(el1,2)-geom(el2,2))*(geom(el1,2)-geom(el2,2));
    _dist += (geom(el1,3)-geom(el2,3))*(geom(el1,3)-geom(el2,3));
    return dist_si2Angstrom*sqrt(_dist);
  }
//--------------------------------- Angle Functions
  std::vector<matrixE> calcAngles() {
    //function that calculates all molecular angles
    std::vector<matrixE> angles;
    matrixE anglesm(natoms,natoms);
    std::vector<double> aref(3,0.0);
    std::vector<double> a1(3,0.0);
    std::vector<double> a2(3,0.0);
    double intprod;
    double norm1;
    double norm2;
    //loop over reference atom. This is unique and has the simmetry of permutation a1 with a2
    for (size_t idxr = 0; idxr < natoms; ++idxr) {
      aref[0] = geometry(idxr+1,1);
      aref[1] = geometry(idxr+1,2);
      aref[2] = geometry(idxr+1,3);
      //loop over a1
      for (size_t idx1 = 0; idx1 < natoms; ++idx1) {
        if (idxr != idx1) {
          //getting a1 with respect to reference
          a1[0] = geometry(idx1+1,1) - aref[0];
          a1[1] = geometry(idx1+1,2) - aref[1];
          a1[2] = geometry(idx1+1,3) - aref[2];
          norm1 = sqrt(a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2]);
          //loop over a2
          for (size_t idx2 = idx1; idx2 < natoms; ++idx2) {
            if ((idxr != idx2)&&(idx1 != idx2)) {
              //getting a2 with respect to reference
              a2[0] = geometry(idx2+1,1) - aref[0];
              a2[1] = geometry(idx2+1,2) - aref[1];
              a2[2] = geometry(idx2+1,3) - aref[2];
              intprod = a1[0]*a2[0] + a1[1]*a2[1] + a1[2]*a2[2];
              norm2 = sqrt(a2[0]*a2[0] + a2[1]*a2[1] + a2[2]*a2[2]);
              anglesm(idx1+1,idx2+1) = acos((intprod)/(norm1*norm2));
              //anglesm(idx1+1,idx2+1) = acos((intprod)/(norm1*norm2))*180/pi;
              anglesm(idx2+1,idx1+1) = anglesm(idx1+1,idx2+1);
            }
            else {
              anglesm(idx1+1,idx2+1) = 0;
              anglesm(idx2+1,idx1+1) = 0;
            }
          }
        }
        else {
          for (size_t idx2 = 0; idx2 < natoms; ++idx2) {
            anglesm(idx2+1,idx1+1) = 0;
            anglesm(idx1+1,idx2+1) = 0;
          }
        }
      }
      angles.push_back(anglesm);
      //std::cout << "The matrix of angles for the reference atom " << idxr+1 << std::endl;
      //anglesm.Print();
    }
    return angles;
  }
  matrixE calcAngles(size_t refa) {
    //function calculating all angles for a fixed reference atom. This atom is used for the origin of referencial
    matrixE anglesm(natoms,natoms);
    std::vector<double> aref(3,0.0);
    std::vector<double> a1(3,0.0);
    std::vector<double> a2(3,0.0);
    double intprod;
    double norm1;
    double norm2;
    aref[0] = geometry(refa,1);
    aref[1] = geometry(refa,2);
    aref[2] = geometry(refa,3);
    //loop over a1
    for (size_t idx1 = 0; idx1 < natoms; ++idx1) {
      if (idx1 != (refa-1)) {
        //getting a1 with respect to reference
        a1[0] = geometry(idx1+1,1) - aref[0];
        a1[1] = geometry(idx1+1,2) - aref[1];
        a1[2] = geometry(idx1+1,3) - aref[2];
        norm1 = sqrt(a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2]);
        //loop over a2
        for (size_t idx2 = 0; idx2 < natoms; ++idx2) {
          if ((idx2 != (refa-1))&&(idx1 != idx2)) {
            //getting a2 with respect to reference
            a2[0] = geometry(idx2+1,1) - aref[0];
            a2[1] = geometry(idx2+1,2) - aref[1];
            a2[2] = geometry(idx2+1,3) - aref[2];
            intprod = a1[0]*a2[0] + a1[1]*a2[1] + a1[2]*a2[2];
            norm2 = sqrt(a2[0]*a2[0] + a2[1]*a2[1] + a2[2]*a2[2]);
            anglesm(idx1+1,idx2+1) = acos((intprod)/(norm1*norm2));
            //anglesm(idx1+1,idx2+1) = acos((intprod)/(norm1*norm2))*180/pi;
            anglesm(idx2+1,idx1+1) = anglesm(idx1+1,idx2+1);
          }
          else {
            anglesm(idx1+1,idx2+1) = 0;
            anglesm(idx2+1,idx1+1) = 0;
          }
        }
      }
      else {
        for (size_t idx2 = 0; idx2 < natoms; ++idx2) {
          anglesm(idx2+1,idx1+1) = 0;
        }
      }
    }
    return anglesm;
  }
  double calcAngles(size_t refa, size_t at1, size_t at2) {
    //function that calculates the angle made by a set of 3 atoms.
    double anglesm = 0.0;
    if ((refa != at1)&&(refa != at2)&&(at2 != at1)) {
      std::vector<double> aref(3,0.0);
      std::vector<double> a1(3,0.0);
      std::vector<double> a2(3,0.0);
      double intprod;
      double norm1;
      double norm2;
      aref[0] = geometry(refa,1);
      aref[1] = geometry(refa,2);
      aref[2] = geometry(refa,3);
      a1[0]   = geometry(at1,1) - aref[0];
      a1[1]   = geometry(at1,2) - aref[1];
      a1[2]   = geometry(at1,3) - aref[2];
      norm1   = sqrt(a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2]);
      a2[0]   = geometry(at2,1) - aref[0];
      a2[1]   = geometry(at2,2) - aref[1];
      a2[2]   = geometry(at2,3) - aref[2];
      intprod = a1[0]*a2[0] + a1[1]*a2[1] + a1[2]*a2[2];
      norm2   = sqrt(a2[0]*a2[0] + a2[1]*a2[1] + a2[2]*a2[2]);
      anglesm = acos((intprod)/(norm1*norm2));
    }
    return anglesm;
  }
  double calcDihedrals(size_t at1, size_t at2, size_t at3, size_t at4, double threshold = 1.0e-7) {
    //function calculating the dihedral angle between 4 atoms; adapted from TMPChem
    double dihedral = 0.0;
    if ((at1 != at2)&&(at1 != at3)&&(at1 != at4)&&(at2 != at3)&&(at2 != at4)&&(at3 != at4)) {
      //calculate vectors for planes
      std::vector<double> atomA(3,0.0);
      std::vector<double> atomB(3,0.0);
      std::vector<double> atomC(3,0.0);
      std::vector<double> atomD(3,0.0);
      for (size_t idx = 1; idx < 4; ++idx) {
        atomA[idx - 1] = geometry(at1,idx);
        atomB[idx - 1] = geometry(at2,idx);
        atomC[idx - 1] = geometry(at3,idx);
        atomD[idx - 1] = geometry(at4,idx);
      }
      std::vector<double> a1(3,0.0);
      std::vector<double> a2(3,0.0);
      a1[0] = (atomB[1] - atomA[1])*(atomC[2] - atomA[2]) - (atomB[2] - atomA[2])*(atomC[1] - atomA[1]);
      a1[1] = (atomB[2] - atomA[2])*(atomC[0] - atomA[0]) - (atomB[0] - atomA[0])*(atomC[2] - atomA[2]);
      a1[2] = (atomB[0] - atomA[0])*(atomC[1] - atomA[1]) - (atomB[1] - atomA[1])*(atomC[0] - atomA[0]);
      a2[0] = (atomC[1] - atomB[1])*(atomD[2] - atomB[2]) - (atomC[2] - atomB[2])*(atomD[1] - atomB[1]);
      a2[1] = (atomC[2] - atomB[2])*(atomD[0] - atomB[0]) - (atomC[0] - atomB[0])*(atomD[2] - atomB[2]);
      a2[2] = (atomC[0] - atomB[0])*(atomD[1] - atomB[1]) - (atomC[1] - atomB[1])*(atomD[0] - atomB[0]);
      //calculate angle
      double n1 = sqrt(a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2]);
      double n2 = sqrt(a2[0]*a2[0] + a2[1]*a2[1] + a2[2]*a2[2]);
      double intprod = a1[0]*a2[0] + a1[1]*a2[1] + a1[2]*a2[2];
      double theta = (intprod)/(n1*n2);
      if (theta > 1.0) {theta = 1.0;}
      else if (theta < -1.0) {theta = -1.0;}
      dihedral = acos(theta);
      //determine sign
      std::vector<double> u21(3,0.0);
      std::vector<double> * pu21 = & u21;
      std::vector<double> u23(3,0.0);
      std::vector<double> * pu23 = & u23;
      std::vector<double> u34(3,0.0);
      std::vector<double> * pu34 = & u34;
      double nu21 = 0.0;
      double nu23 = 0.0;
      double nu34 = 0.0;
      double aux;
      //get unit vectors between atoms
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        u21[idcoord] = geometry(at1,idcoord + 1) - geometry(at2,idcoord + 1);
        nu21 += u21[idcoord]*u21[idcoord];
        u23[idcoord] = geometry(at3,idcoord + 1) - geometry(at2,idcoord + 1);
        nu23 += u23[idcoord]*u23[idcoord];
        u34[idcoord] = geometry(at4,idcoord + 1) - geometry(at3,idcoord + 1);
        nu34 += u34[idcoord]*u34[idcoord];
      }
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        u21[idcoord] /= sqrt(nu21);
        u23[idcoord] /= sqrt(nu23);
        u34[idcoord] /= sqrt(nu34);
      }
      aux = dotProd(pu23,pu21);
      std::vector<double> u21xu23 = crossProd(u21,u23,sqrt(1.0 - aux*aux));
      std::vector<double> * pu21xu23 = & u21xu23;
      double sign = 1.0;
      if (dotProd(pu21xu23,pu34) > threshold) {sign = -1.0;}
      dihedral *= sign;
    }
    return dihedral;
  }
//--------------------------------- Mass Functions
  void masses() {
    //function to determine the molecular masses (mol and single molecule)
    MWeight = 0.0;
    for (size_t idx = 0; idx < natoms; ++idx) {
      MWeight += Weight(atoms[idx]);
    }
    mass = MWeight/NA;
  }
  double AtomicMasses(std::vector<double> & amass, bool convert2au = false, double scaleproton = 1.0) {
    //function that populates a vector with atomic masses
    if (amass.size() != natoms) {amass.resize(natoms);}
    double totalmass = 0.0;
    double aux;
    double unitconv = 1.0;
    if (convert2au) {unitconv = 1.0/(au2kg*NA);}         //convert to atomic units
    for (size_t idAtm = 0; idAtm < natoms; ++idAtm) {
      aux = Weight(atoms[idAtm])*unitconv;
      totalmass += aux;
      if (atoms[idAtm] == 1) {aux *= scaleproton;}
      amass[idAtm] = aux;
    }
    return totalmass;
  }
  double MW() {return MWeight;}
  double MW() const {return MWeight;}
  double m() {return mass;}
  double m() const {return mass;}
//--------------------------------- Inertia Functions
  std::vector<double> InertiaEigenvalues(bool inmeter = true) {
    //function calculating the eigenvalues of the molecular inertia matrix
    matrixE _inert = InertiaMatrix(geometry,atoms,inmeter);
    std::vector<double> eigenvalues = MatDiag(_inert);
    return eigenvalues;
  }
  matrixE Inertia(bool inmeter = true, bool scale = false) {
    //function calculating the molecular inertia matrix
    matrixE _inert = InertiaMatrix(geometry,atoms,inmeter);
    if (scale) {_inert *= 1.0e44;}
    return _inert;
  }
//--------------------------------- "particle" counting functions
  void CountElectrons() {
    //function to calculate the number of electrons
    nelectrons = 0;
    for (size_t idx = 0; idx < natoms; ++idx) {
      nelectrons += atoms[idx];
    }
    nelectrons -= charge;
  }
  size_t Natoms() {return natoms;}
  size_t Natoms() const {return natoms;}
  void calcNatoms() {natoms = atoms.size();}
  int Nelectrons() {return nelectrons;}
  int Nelectrons() const {return nelectrons;}
  size_t Multiplicity() {return mult;}
  size_t Multiplicity() const {return mult;}
  void setMultiplicity(size_t _mult) {mult = _mult;}
  int Charge() {return charge;}
  void setCharge(int _charge) {
    charge = _charge;
    CountElectrons();
  }
  int Charge() const {return charge;}
//--------------------------------- Vector Related Functions
  std::vector<double> oVector2Atom(matrixE gm, size_t atm1, size_t atm2) {
    //function that returns the vector orthogonal to the plane made by 2 atoms and the origin
    //atmi are atom numbers (start on 1)
    std::vector<double> atomB(3,0.0);
    std::vector<double> atomC(3,0.0);
    for (size_t idx = 1; idx < 4; ++idx) {
      atomB[idx-1] = gm(atm1,idx);
      atomC[idx-1] = gm(atm2,idx);
    }
    std::vector<double> a1(3,0.0);
    a1[0] = atomB[1]*atomC[2] - atomB[2]*atomC[1];
    a1[1] = atomB[2]*atomC[0] - atomB[0]*atomC[2];
    a1[2] = atomB[0]*atomC[1] - atomB[1]*atomC[0];
    double norm = sqrt(a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2]);
    a1[0] /= norm;
    a1[1] /= norm;
    a1[2] /= norm;
    return a1;
  }
  std::vector<double> oVector3Atom(matrixE gm, size_t atm1, size_t atm2, size_t atm3) {
    //function that returns the vector orthogonal to the plane made by 3 atoms
    //atmi are atom numbers (start on 1)
    std::vector<double> atomA(3,0.0);
    std::vector<double> atomB(3,0.0);
    std::vector<double> atomC(3,0.0);
    for (size_t idx = 1; idx < 4; ++idx) {
      atomA[idx-1] = gm(atm1,idx);
      atomB[idx-1] = gm(atm2,idx);
      atomC[idx-1] = gm(atm3,idx);
    }
    std::vector<double> a1(3,0.0);
    a1[0] = (atomB[1]-atomA[1])*(atomC[2]-atomA[2]) - (atomB[2]-atomA[2])*(atomC[1]-atomA[1]);
    a1[1] = (atomB[2]-atomA[2])*(atomC[0]-atomA[0]) - (atomB[0]-atomA[0])*(atomC[2]-atomA[2]);
    a1[2] = (atomB[0]-atomA[0])*(atomC[1]-atomA[1]) - (atomB[1]-atomA[1])*(atomC[0]-atomA[0]);
    double norm = sqrt(a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2]);
    a1[0] /= norm;
    a1[1] /= norm;
    a1[2] /= norm;
    return a1;
  }
  std::vector<double> oVectorPlane(double tolerance = pi/180) {
    //function that returns the vector orthogonal to the molecular plane
    double angle;
    bool planar = false;
    bool swtch = false;
    size_t atm1;
    size_t atm2;
    size_t atm3;
    size_t atm4;
    if (natoms == 3) {
      atm1 = 1;
      atm2 = 2;
      atm3 = 3;
    }
    else {
      //loop over atoms and create lines of atoms connected
      for (size_t idx1 = 0; idx1 < natoms; ++idx1) {
        atm1 = idx1+1;
        for (size_t idx2 = 0; idx2 < natoms; ++idx2) {
          atm2 = idx2+1;
          if (atm2 == atm1) {continue;}
          for (size_t idx3 = 0; idx3 < natoms; ++idx3) {
            atm3 = idx3+1;
            if ((atm3 == atm1)||(atm3 == atm2)) {continue;}
            for (size_t idx4 = 0; idx4 < natoms; ++idx4) {
              atm4 = idx4+1;
              if ((atm4 == atm1)||(atm4 == atm2)||(atm4 == atm3)) {continue;}
              angle = calcDihedrals(atm1,atm2,atm3,atm4);
              if (fabs(angle) < pi/2) {
                if (fabs(angle) < tolerance) {
                  planar = true;
                  swtch = false;
                  break;
                }
                else {
                  angle = calcDihedrals(idx1+1,atm3,atm4,atm2);
                  if (fabs(angle) < tolerance) {
                    planar = true;
                    swtch = true;
                    break;
                  }
                }
              }
              else {
                if (fabs(fabs(angle) - pi) < tolerance) {
                  planar = true;
                  swtch = false;
                  break;
                }
                else {
                  angle = calcDihedrals(idx1+1,atm3,atm4,atm2);
                  if ((fabs(angle) > pi/2)&&(fabs(fabs(angle) - pi) < tolerance)) {
                    planar = true;
                    swtch = true;
                    break;
                  }
                  else if ((fabs(angle) < pi/2)&&(fabs(angle) < tolerance)) {
                    planar = true;
                    swtch = true;
                    break;
                  }
                }
              }
            }
            if (planar == true) {break;}
          }
          if (planar == true) {break;}
        }
        if (planar == true) {break;}
      }
    }
    //getting one of the planes
    if (swtch) {
      atm2 = atm3;
      atm3 = atm4;
    }
    return oVector3Atom(geometry,atm1,atm2,atm3);
  }
  void setPGroup(const std::string & _pgroup) {symm = _pgroup;}
};
unsigned int PointGroup2Sigma(std::string pg) {
  //http://symmetry.jacobs-university.de
  //pg is point group as a string
  int _sigma = 0;
  if (pg == "DINFh") {_sigma = 2;}
  else if (pg == "CINFv") {_sigma = 1;}
  else if (pg == "C2h") {_sigma = 2;}
  else if (pg == "C3h") {_sigma = 3;}
  else if (pg == "C4h") {_sigma = 4;}
  else if (pg == "C5h") {_sigma = 5;}
  else if (pg == "C6h") {_sigma = 6;}
  else if (pg == "C2v") {_sigma = 2;}
  else if (pg == "C3v") {_sigma = 3;}
  else if (pg == "C4v") {_sigma = 4;}
  else if (pg == "C5v") {_sigma = 5;}
  else if (pg == "C6v") {_sigma = 6;}
  else if (pg == "C7v") {_sigma = 7;}
  else if (pg == "C8v") {_sigma = 8;}
  else if (pg == "C1") {_sigma = 1;}
  else if (pg == "C2") {_sigma = 2;}
  else if (pg == "C3") {_sigma = 3;}
  else if (pg == "C4") {_sigma = 4;}
  else if (pg == "C5") {_sigma = 5;}
  else if (pg == "C6") {_sigma = 6;}
  else if (pg == "C7") {_sigma = 7;}
  else if (pg == "C8") {_sigma = 8;}
  else if (pg == "Cs") {_sigma = 1;}
  else if (pg == "Ci") {_sigma = 1;}
  else if (pg == "D2h") {_sigma = 4;}
  else if (pg == "D3h") {_sigma = 6;}
  else if (pg == "D4h") {_sigma = 8;}
  else if (pg == "D5h") {_sigma = 10;}
  else if (pg == "D6h") {_sigma = 12;}
  else if (pg == "D7h") {_sigma = 14;}
  else if (pg == "D8h") {_sigma = 16;}
  else if (pg == "D2d") {_sigma = 4;}
  else if (pg == "D3d") {_sigma = 6;}
  else if (pg == "D4d") {_sigma = 8;}
  else if (pg == "D5d") {_sigma = 10;}
  else if (pg == "D6d") {_sigma = 12;}
  else if (pg == "D7d") {_sigma = 14;}
  else if (pg == "D8d") {_sigma = 16;}
  else if (pg == "D2") {_sigma = 4;}
  else if (pg == "D3") {_sigma = 6;}
  else if (pg == "D4") {_sigma = 8;}
  else if (pg == "D5") {_sigma = 10;}
  else if (pg == "D6") {_sigma = 12;}
  else if (pg == "D7") {_sigma = 14;}
  else if (pg == "D8") {_sigma = 16;}
  else if (pg == "T") {_sigma = 12;}
  else if (pg == "Th") {_sigma = 12;}    //use the same as T
  else if (pg == "Td") {_sigma = 12;}
  else if (pg == "O") {_sigma = 24;}     //use the same as Oh
  else if (pg == "Oh") {_sigma = 24;}
  else if (pg == "I") {_sigma = 60;}     //use the same as Ih
  else if (pg == "Ih") {_sigma = 60;}
  else if (pg == "S2") {_sigma = 1;}
  else if (pg == "S4") {_sigma = 2;}
  else if (pg == "S6") {_sigma = 3;}
  else if (pg == "S8") {_sigma = 4;}
  else if (pg == "S10") {_sigma = 5;}
  else if (pg == "S12") {_sigma = 6;}
  return _sigma;
}

#endif //_Molecule_
