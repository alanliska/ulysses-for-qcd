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

#ifndef _File_Format_Package_
#define _File_Format_Package_
#include <vector>
#include <fstream>
#include <iostream>
#include "BiologicalData.hpp"
#include "../atoms/AtomPackage.hpp"

//description:
//functions that read or write specific file-types associated to molecules

void ReadXYZFormat(std::string geomfile, size_t & natoms, matrixE & geometry, std::vector<size_t> & atoms) {
  //function to read geometries and atomic lists from a given xyz file
  std::ifstream gfile(geomfile,std::ios::in);
  if (!gfile.is_open()) {throw std::string("ERROR: FileFormats.hpp: ReadXYZFormat(): xyz file could not be open");}
  std::string read;
  std::getline(gfile,read);
  std::istringstream (read) >> natoms;
  std::getline(gfile,read);           //skip comment
  geometry.resize(natoms,3);
  atoms.resize(natoms);
  for (size_t idatm = 0; idatm < natoms; ++idatm) {
    gfile >> read;
    atoms[idatm] = Symbol2AtomNr(read);
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      gfile >> geometry(idatm + 1,idcoord + 1);
    }
  }
  gfile.close();
}
void WriteXYZFormat(std::string name, std::vector<size_t> & atoms, matrixE & geometry, int counter = -1, int prc = 7) {
    //function to write geometries into a given xyz file
    std::string filename = name;
    double precision = pow(10.0,-prc);
    double number;
    int natoms = atoms.size();
    if (counter >= 0) {filename += szt2str(counter);}
    filename += ".xyz";
    std::ofstream gfile(filename,std::ios::out);
    if (!gfile.is_open()) {throw std::string("ERROR: FileFormats.hpp: WriteXYZFormat(): could not open file for writing geometry");}
    gfile << std::fixed;
    gfile << std::setprecision(prc);
    gfile << natoms << "\n";
    gfile << "\n";
    for (size_t idx = 0; idx < natoms; ++idx) {
      filename = AtomNr2Symbol(atoms[idx]);
      gfile << filename << "    ";
      if (filename.length() == 1) {gfile << " ";}
      for (size_t idx2 = 0; idx2 < 3; ++idx2) {
        number = geometry(idx + 1,idx2 + 1);
        if (fabs(number) < precision) {number = 0.0;}
        if (number >= -precision) {gfile << " ";}
        gfile << number << "    ";
      }
      gfile << "\n";
    }
    gfile.close();
}
void ReadMOL2Format(std::string mol2file, size_t & Natoms, matrixE & geometry, std::vector<size_t> & atoms, int maxnelements = 100000000) {
  //function to read geometries and atomic lists from a given mol2 file
  std::ifstream gfile(mol2file,std::ios::in);
  if (!gfile.is_open()) {throw std::string("ERROR: FileFormats.hpp: ReadMOL2Format(): mol2 file could not be open");}
  std::string read;
  std::string discard;
  std::string build;
  std::string newelement;
  char element;
  bool proteinnaming = false;
  //move to where geometry begins
  for (size_t idline = 0; idline < maxnelements; ++idline) {
    std::getline(gfile,read);
    if (read == "@<TRIPOS>ATOM") {break;}
  }
  //now read geometry
  Natoms = 1;
  for (size_t idchar = 0; idchar < maxnelements; ++idchar) {
    gfile >> read;
    if (read == "@<TRIPOS>BOND") {break;}
    geometry.resize(Natoms,3);
    atoms.resize(Natoms);
    discard = read;
    gfile >> read;
    //process atomic number; mol2 adds stuff after the elements, take this out
    build = "";
    for (size_t idel = 0; idel < read.size(); ++idel) {
      element = static_cast<unsigned char>(read[idel]);
      if (isalpha(element)) {build += read[idel];}
      else {break;}
    }
    if (ValidProteinAtom(build)) {proteinnaming = true;}
    if (proteinnaming) {newelement = Convert2Element(build);}
    else {newelement = build;}
    atoms[Natoms - 1] = Symbol2AtomNr(newelement);
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      geometry(Natoms,idcoord + 1) = 0.0;
      gfile >> geometry(Natoms,idcoord + 1);
    }
    for (size_t idrepeat = 0; idrepeat < 4; ++idrepeat) {
      gfile >> read;
      discard = read;
    }
    ++Natoms;
  }
  --Natoms;
  gfile.close();
}
void FullReadMOL2Format(std::string mol2file, int & Natoms, matrixE & geometry, std::vector<size_t> & atoms, matrix<int> & Connectivity, matrixE & BondOrders, int maxneigh = 8, int maxnelements = 100000000) {
  //function to read geometries, atomic lists, connectivities and bond orders from a given mol2 file
  std::ifstream gfile(mol2file,std::ios::in);
  if (!gfile.is_open()) {throw std::string("ERROR: FullReadMOL2Format.hpp: ReadMOL2Format(): mol2 file could not be open");}
  std::string read;
  std::string discard;
  std::string build;
  std::string newelement;
  int atmA;
  int atmB;
  char element;
  double bo;
  bool proteinnaming = false;
  //move to where geometry begins
  for (size_t idline = 0; idline < maxnelements; ++idline) {
    std::getline(gfile,read);
    if (read == "@<TRIPOS>ATOM") {break;}
  }
  //now read geometry
  Natoms = 1;
  for (size_t idchar = 0; idchar < maxnelements; ++idchar) {
    gfile >> read;
    if (read == "@<TRIPOS>BOND") {break;}
    geometry.resize(Natoms,3);
    atoms.resize(Natoms);
    discard = read;
    gfile >> read;
    //process atomic number; mol2 adds stuff after the elements, take this out
    build = "";
    for (size_t idel = 0; idel < read.size(); ++idel) {
      element = static_cast<unsigned char>(read[idel]);
      if (isalpha(element)) {build += read[idel];}
      else {break;}
    }
    if (ValidProteinAtom(build)) {proteinnaming = true;}
    if (proteinnaming) {newelement = Convert2Element(build);}
    else {newelement = build;}
    atoms[Natoms - 1] = Symbol2AtomNr(newelement);
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      geometry(Natoms,idcoord + 1) = 0.0;
      gfile >> geometry(Natoms,idcoord + 1);
    }
    for (size_t idrepeat = 0; idrepeat < 4; ++idrepeat) {
      gfile >> read;
      discard = read;
    }
    ++Natoms;
  }
  --Natoms;
  //prepare connectivity and bond orders
  std::vector<int> neighbours(Natoms,0);
  Connectivity.resize(Natoms,maxneigh);           //assume a maximum of 8 neighbours per atom
  BondOrders.resize(Natoms,maxneigh);             //assume a maximum of 8 neighbours per atom
  //now read them
  for (size_t idchar = 0; idchar < maxnelements; ++idchar) {
    gfile >> read;
    if (read == "@<TRIPOS>SUBSTRUCTURE") {break;}
    discard = read;
    gfile >> atmA;
    gfile >> atmB;
    gfile >> read;
    std::cout << discard << " " << atmA << " " << atmB << " " << read << std::endl;
    if (read == "am") {bo = 1.0;}
    else if (read == "ar") {bo = 1.5;}
    else if (read == "1") {bo = 1.0;}
    else if (read == "2") {bo = 2.0;}
    else if (read == "3") {bo = 3.0;}
    ++neighbours[atmA - 1];
    ++neighbours[atmB - 1];
    Connectivity(atmA,neighbours[atmA - 1]) = atmB;
    Connectivity(atmB,neighbours[atmB - 1]) = atmA;
    BondOrders(atmA,neighbours[atmA - 1]) = bo;
    BondOrders(atmB,neighbours[atmB - 1]) = bo;
  }
  gfile.close();
}
void ReadSDFFormat(std::string geomfile, size_t & natoms, matrixE & geometry, std::vector<size_t> & atoms, int & charge, int maxnelements = 100000000) {
  //function to read geometries and atomic lists from a given sdf file
  std::ifstream gfile(geomfile,std::ios::in);
  if (!gfile.is_open()) {throw std::string("ERROR: FileFormats.hpp: ReadSDFFormat(): sdf file could not be open");}
  int sizewords = 120;          //molecule with 30 atoms
  int increment = sizewords;
  int start = -1;
  int countbetweenatoms = 0;
  int auxcount = 0;
  int firstatom = -1;
  bool e0alpha;
  bool e1alpha;
  std::vector<std::string> read(increment);
  //read the whole file to string
  for (size_t idchar = 0; idchar < maxnelements; ++idchar) {
    if (idchar + 1 == sizewords) {
      sizewords += increment;
      read.resize(sizewords);
    }
    gfile >> read[idchar];
    if (read[idchar] == "$$$$") {
      sizewords = idchar + 1;
      break;
    }
  }
  gfile.close();
  //determine where to start reading the geometry
  for (size_t idel = 0; idel < sizewords; ++idel) {
    if ((read[idel] == "V2000")||(read[idel] == "v2000")) {
      start = idel + 1;
      break;
    }
  }
  //determine what is the frequency with which atoms are declared
  for (size_t idel = start; idel < sizewords; ++idel) {
    if (read[idel].size() <= 2) {
      e0alpha = isalpha(read[idel][0]);
      e1alpha = true;
      if (read[idel].size() == 2) {e1alpha = isalpha(read[idel][1]);}
      if ((e0alpha)&&(e1alpha)) {
        ++auxcount;
        if (auxcount == 1) {firstatom = idel;}
      }
    }
    if (auxcount == 1) {++countbetweenatoms;}
    else if (auxcount == 2) {break;}
  }
  //count the number of atoms; use the fact that after coordinates there are bond orders
  natoms = 0;
  for (size_t idel = firstatom; idel < sizewords; idel += countbetweenatoms) {
    if (!isalpha(read[idel][0])) {break;}
    ++natoms;
  }
  geometry.resize(natoms,3);
  atoms.resize(natoms);
  firstatom = 0;
  //now write geometry and atom list
  for (size_t idel = start; idel < sizewords; ++idel) {
    if (isalpha(read[idel][0])) {
      atoms[firstatom] = Symbol2AtomNr(read[idel]);
      geometry(firstatom + 1,1) = stod(read[idel - 3]);
      geometry(firstatom + 1,2) = stod(read[idel - 2]);
      geometry(firstatom + 1,3) = stod(read[idel - 1]);
      ++firstatom;
    }
    if (firstatom == natoms) {
      //store last position to check whether charge was declared
      firstatom = idel + 1;
      break;
    }
  }
  //check for potentially declared charge
  for (size_t idel = firstatom; idel < sizewords; ++idel) {
    if (read[idel] == "<Formal_Charge>") {
      //charge was indeed declared
      charge = stoi(read[idel + 1]);
    }
  }
}
void WriteSDFFormat(std::string basename, std::vector<size_t> & atoms, matrixE & geometry, matrix<int> & Connectivity, matrixE & BondOrders, int counter = -1, int prc = 7) {
  //function writing an sdf file
  std::string filename = basename;
  std::string spaces = "  ";
  double number;
  double fnumber;
  double precision = pow(10.0,-prc);
  int iA;
  int iB;
  int natoms = atoms.size();
  int nbonds = 0;
  int maxnneigh = Connectivity.cols();
  bool negative;
  if (counter >= 0) {filename += szt2str(counter);}
  filename += ".sdf";
  std::ofstream gfile(filename,std::ios::out);
  if (!gfile.is_open()) {throw std::string("ERROR: FileFormats.hpp: WriteSDFFormat(): could not open file for writing geometry");}
  gfile << std::fixed;
  //write general header
  gfile << basename << "\nwritten by ULYSSES\n\n";
  //get number of atoms and number of bonds
  for (size_t idAtm = 0; idAtm < natoms; ++idAtm) {
    for (size_t idneigh = 0; idneigh < maxnneigh; ++idneigh) {
      if ((Connectivity(idAtm + 1,idneigh + 1) != 0)&&(Connectivity(idAtm + 1,idneigh + 1) - 1 < idAtm)) {++nbonds;}
      else {break;}
    }
  }
  //write molecule header
  if (natoms > 9) {spaces = " ";}
  if (natoms > 99) {spaces = "";}
  gfile << spaces << natoms;
  spaces = "  ";
  if (nbonds > 9) {spaces = " ";}
  if (nbonds > 99) {spaces = "";}
  gfile << spaces << nbonds << "  0  0  0  0  0  0  0  0999 V2000\n";
  gfile << std::setprecision(4);
  //write geometry
  for (size_t idAtm = 0; idAtm < natoms; ++idAtm) {
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      spaces = "    ";
      number = geometry(idAtm + 1,idcoord + 1);
      fnumber = fabs(number);
      if (fnumber < precision) {
        number = 0.0;
        fnumber = 0.0;
      }
      iA = int(fnumber);
      negative = (number < 0.0);
      if (iA < 10) {
        if (negative) {spaces = "   ";}
      }
      else if ((iA >= 10)&&(iA < 100)) {
        if (negative) {spaces = "  ";}
        else {spaces = "   ";}
      }
      else if ((iA >= 100)&&(iA < 1000)) {
        if (negative) {spaces = " ";}
        else {spaces = "  ";}
      }
      gfile << spaces << number;
    }
    filename = AtomNr2Symbol(atoms[idAtm]);
    spaces = "   0  0  0  0  0";
    if (filename.length() == 2) {spaces = "  0  0  0  0  0";}
    gfile << " " << filename << spaces << "\n";
  }
  //write bonds
  for (size_t idAtm = 0; idAtm < natoms; ++idAtm) {
    iA = idAtm + 1;
    for (size_t idneigh = 0; idneigh < maxnneigh; ++idneigh) {
      iB = Connectivity(idAtm + 1,idneigh + 1);
      if (iB == 0) {continue;}
      else if (iB < iA) {continue;}
      spaces = "  ";
      if ((iA > 9)&&(iA < 100)) {spaces = " ";}
      else if (iA > 99) {spaces = "";}
      gfile << spaces << iA;
      spaces = "  ";
      if ((iB > 9)&&(iB < 100)) {spaces = " ";}
      else if (iB > 99) {spaces = "";}
      gfile << spaces << iB << "  ";
      iB = int(BondOrders(idAtm + 1,idneigh + 1));
      if (BondOrders(idAtm + 1,idneigh + 1) == 1.5) {iB = 4;}
      gfile << iB << "  0  0  0\n";
    }
  }
  gfile << "M  END\n$$$$\n";
  gfile.close();
}
void ReadPDBFormat(std::string geomfile, size_t & natoms, matrixE & geometry, std::vector<size_t> & atoms, long int maxnelements = 10000000000) {
  //function that reads PDB files for geometries and atoms
  std::ifstream gfile(geomfile,std::ios::in);
  if (!gfile.is_open()) {throw std::string("ERROR: FileFormats.hpp: ReadPDBFormat(): pdb file could not be open");}
  int sizewords = 50000;
  int increment = 5000;
  int start = 0;
  int countbetweenatoms = 0;
  int skipanisotropy = 0;
  bool anisotropy = false;
  std::string rline;
  std::string line4;
  std::string aux;
  std::vector<std::string> read(sizewords);
  //move closer to where geometry begins
  for (size_t idline = 0; idline < maxnelements; ++idline) {
    std::getline(gfile,rline);
    line4 = "";
    for (size_t idchar = 0; idchar < 4; ++idchar) {
      line4 += rline[idchar];
    }
    if ((line4 == "ATOM")||(line4 == "Atom")||(line4 == "atom")) {break;}
  }
  line4 = "";
  for (size_t idchar = 0; idchar < rline.size(); ++idchar) {
    aux = rline[idchar];
    if (aux == " ") {
      if (line4.size() != 0) {
        read[start] = line4;
        ++start;
      }
      line4 = "";
      continue;
    }
    line4 += rline[idchar];
  }
  if (line4.size() != 0) {
    read[start] = line4;
    ++start;
  }
  //read the rest of file to string
  for (size_t idchar = start; idchar < maxnelements; ++idchar) {
    if (idchar + 1 == sizewords) {
      sizewords += increment;
      read.resize(sizewords);
    }
    gfile >> read[idchar];
    if ((read[idchar] == "END")||(read[idchar] == "End")||(read[idchar] == "end")) {
      sizewords = idchar + 1;
      break;
    }
  }
  gfile.close();
  //determine count between atoms
  for (int idel = 1; idel < sizewords; ++idel) {
    if ((read[idel] == "ATOM")||(read[idel] == "Atom")||(read[idel] == "atom")) {
      if (!anisotropy) {countbetweenatoms = idel;}
      else {skipanisotropy = idel - countbetweenatoms;}
      break;
    }
    else if ((read[idel] == "ANISOU")||(read[idel] == "Anisou")||(read[idel] == "anisou")) {
      countbetweenatoms = idel;
      anisotropy = true;
    }
  }
  if (countbetweenatoms < 12) {countbetweenatoms = 12;}
  //get the geometry and atom list; start by allocating excess of memory for the quantities
  geometry.resize(sizewords/countbetweenatoms,3);
  atoms.resize(sizewords/countbetweenatoms);
  natoms = 0;
  for (size_t idel = 0; idel < sizewords; idel += countbetweenatoms) {
    if ((read[idel] == "TER")||(read[idel] == "Ter")||(read[idel] == "ter")) {
      //shift further
      for (size_t idel2 = idel; idel2 < sizewords; ++idel2) {
        if ((read[idel2] == "ATOM")||(read[idel2] == "Atom")||(read[idel2] == "atom")||(read[idel2] == "HETATM")||(read[idel2] == "Hetatm")||(read[idel2] == "hetatm")) {
          idel = idel2;
          break;
        }
      }
    }
    if ((read[idel] == "ANISOU")||(read[idel] == "Anisou")||(read[idel] == "anisou")) {idel += skipanisotropy;}
    if ((read[idel] == "ATOM")||(read[idel] == "Atom")||(read[idel] == "atom")||(read[idel] == "HETATM")||(read[idel] == "Hetatm")||(read[idel] == "hetatm")) {
      atoms[natoms] = Symbol2AtomNr(read[idel + countbetweenatoms - 1]);
      geometry(natoms + 1,1) = stod(read[idel + countbetweenatoms - 6]);
      geometry(natoms + 1,2) = stod(read[idel + countbetweenatoms - 5]);
      geometry(natoms + 1,3) = stod(read[idel + countbetweenatoms - 4]);
      if (atoms[natoms] == 0) {
        //error reading atom type, search forward and fix the counter
        for (size_t iderr = 0; iderr < countbetweenatoms; ++iderr) {
          if (Symbol2AtomNr(read[idel + countbetweenatoms + iderr]) > 0) {
            atoms[natoms] = Symbol2AtomNr(read[idel + countbetweenatoms + iderr]);
            idel -= iderr - 1;
            break;
          }
        }
      }
      ++natoms;
    }
  }
  geometry.resize(natoms,3);
  atoms.resize(natoms);
}
void FullReadPDBFormat(std::string geomfile, int & natoms, matrixE & geometry, std::vector<size_t> & atoms, int & nheteroatoms, matrixE & hetgeom, std::vector<size_t> & hetatoms, std::vector<std::string> & atmtype, matrix<int> & chains, std::vector<std::string> & residue, std::vector<int> & res_start, bool & hasH, int & nres, long int maxnelements = 10000000000) {
  //function that reads PDB files for geometries and atoms
  std::ifstream gfile(geomfile,std::ios::in);
  if (!gfile.is_open()) {throw std::string("ERROR: FileFormats.hpp: FullReadPDBFormat(): pdb file could not be open");}
  int sizewords = 50000;
  int increment = 5000;
  int start = 0;
  int countbetweenatoms = 0;
  int skipanisotropy = 0;
  int ichain = 0;
  int ires = 0;
  bool Cdone = false;
  bool CAdone = false;
  bool Ndone = false;
  bool Odone = false;
  bool newresidue = false;
  bool valid_aa;
  bool valid_fr;
  bool anisotropy = false;
  std::string rline;
  std::string line4;
  std::string aux;
  std::vector<std::string> read(sizewords);
  //move closer to where geometry begins
  for (size_t idline = 0; idline < maxnelements; ++idline) {
    std::getline(gfile,rline);
    line4 = "";
    for (size_t idchar = 0; idchar < 4; ++idchar) {
      line4 += rline[idchar];
    }
    if ((line4 == "ATOM")||(line4 == "Atom")||(line4 == "atom")) {break;}
  }
  line4 = "";
  for (size_t idchar = 0; idchar < rline.size(); ++idchar) {
    aux = rline[idchar];
    if (aux == " ") {
      if (line4.size() != 0) {
        read[start] = line4;
        ++start;
      }
      line4 = "";
      continue;
    }
    line4 += rline[idchar];
  }
  if (line4.size() != 0) {
    read[start] = line4;
    ++start;
  }
  //read the rest of file to string
  for (size_t idchar = start; idchar < maxnelements; ++idchar) {
    if (idchar + 1 == sizewords) {
      sizewords += increment;
      read.resize(sizewords);
    }
    gfile >> read[idchar];
    if ((read[idchar] == "END")||(read[idchar] == "End")||(read[idchar] == "end")) {
      sizewords = idchar + 1;
      break;
    }
  }
  gfile.close();
  //determine count between atoms
  for (int idel = 1; idel < sizewords; ++idel) {
    if ((read[idel] == "ATOM")||(read[idel] == "Atom")||(read[idel] == "atom")) {
      if (!anisotropy) {countbetweenatoms = idel;}
      else {skipanisotropy = idel - countbetweenatoms;}
      break;
    }
    else if ((read[idel] == "ANISOU")||(read[idel] == "Anisou")||(read[idel] == "anisou")) {
      countbetweenatoms = idel;
      anisotropy = true;
    }
  }
  if (countbetweenatoms < 12) {countbetweenatoms = 12;}
  //get the geometry and atom list; start by allocating excess of memory for the quantities
  start = sizewords/countbetweenatoms;
  geometry.resize(start,3);
  atoms.resize(start);
  hetgeom.resize(start,3);
  hetatoms.resize(start);
  atmtype.resize(start);
  chains.resize(100,2);                          //just a large enough number
  residue.resize(start);
  natoms = 0;
  nheteroatoms = 0;
  nres = 0;
  hasH = false;
  chains(ichain + 1,1) = 1;
  for (size_t idel = 0; idel < sizewords; idel += countbetweenatoms) {
    if ((read[idel] == "TER")||(read[idel] == "Ter")||(read[idel] == "ter")) {
      //shift further
      chains(ichain + 1,2) = natoms;
      ++ichain;
      chains(ichain + 1,1) = natoms + 1;
      for (size_t idel2 = idel; idel2 < sizewords; ++idel2) {
        if ((read[idel2] == "ATOM")||(read[idel2] == "Atom")||(read[idel2] == "atom")||(read[idel2] == "HETATM")||(read[idel2] == "Hetatm")||(read[idel2] == "hetatm")) {
          idel = idel2;
          break;
        }
      }
    }
    if ((read[idel] == "ANISOU")||(read[idel] == "Anisou")||(read[idel] == "anisou")) {idel += skipanisotropy;}
    if ((read[idel] == "ATOM")||(read[idel] == "Atom")||(read[idel] == "atom")) {
      atoms[natoms] = Symbol2AtomNr(read[idel + countbetweenatoms - 1]);
      if (atoms[natoms] == 1) {hasH = true;}
      atmtype[natoms] = read[idel + 2];
      if (atmtype[natoms] == "CA") {++nres;}                //count the number of residues from the Calpha
      residue[natoms] = read[idel + 3];
      geometry(natoms + 1,1) = stod(read[idel + countbetweenatoms - 6]);
      geometry(natoms + 1,2) = stod(read[idel + countbetweenatoms - 5]);
      geometry(natoms + 1,3) = stod(read[idel + countbetweenatoms - 4]);
      if (atoms[natoms] == 0) {
        //error reading atom type, search forward and fix the counter
        for (size_t iderr = 0; iderr < countbetweenatoms; ++iderr) {
          if (Symbol2AtomNr(read[idel + countbetweenatoms + iderr]) > 0) {
            atoms[natoms] = Symbol2AtomNr(read[idel + countbetweenatoms + iderr]);
            idel -= iderr - 1;
            break;
          }
        }
      }
      ++natoms;
    }
    else if ((read[idel] == "HETATM")||(read[idel] == "Hetatm")||(read[idel] == "hetatm")) {
      hetatoms[nheteroatoms] = Symbol2AtomNr(read[idel + countbetweenatoms - 1]);
      hetgeom(nheteroatoms + 1,1) = stod(read[idel + countbetweenatoms - 6]);
      hetgeom(nheteroatoms + 1,2) = stod(read[idel + countbetweenatoms - 5]);
      hetgeom(nheteroatoms + 1,3) = stod(read[idel + countbetweenatoms - 4]);
      if (hetatoms[nheteroatoms] == 0) {
        //error reading atom type, search forward and fix the counter
        for (size_t iderr = 0; iderr < countbetweenatoms; ++iderr) {
          if (Symbol2AtomNr(read[idel + countbetweenatoms + iderr]) > 0) {
            hetatoms[nheteroatoms] = Symbol2AtomNr(read[idel + countbetweenatoms + iderr]);
            idel -= iderr - 1;
            break;
          }
        }
      }
      ++nheteroatoms;
    }
  }
  std::cout << "nres = " << nres << std::endl;
  //now lets reduce memory demands of residues and where they start
  res_start.resize(nres);
  start = 0;
  for (int idAtm = 0; idAtm < natoms; ++idAtm) {
    valid_aa = ValidAminoAcids(residue[idAtm]);
    valid_fr = ValidFragment(residue[idAtm]);
    if (valid_fr) {
      //fragments have fixed size, always
      if (residue[idAtm] == "ACE") {}
      else if (residue[idAtm] == "NME") {}
      else if (residue[idAtm] == "NHE") {}
    }
    else if (valid_aa) {
      if (atmtype[idAtm] == "C") {
        if (!Cdone) {Cdone = true;}
        else {newresidue = true;}
      }
      else if (atmtype[idAtm] == "CA") {
        if (!CAdone) {CAdone = true;}
        else {newresidue = true;}
      }
      else if (atmtype[idAtm] == "N") {
        if (!Ndone) {Ndone = true;}
        else {newresidue = true;}
      }
      else if (atmtype[idAtm] == "O") {
        if (!Odone) {Odone = true;}
        else {newresidue = true;}
      }
      if (newresidue) {
        Cdone = (atmtype[idAtm] == "C");
        CAdone = (atmtype[idAtm] == "CA");
        Ndone = (atmtype[idAtm] == "N");
        Odone = (atmtype[idAtm] == "O");
        newresidue = false;
        res_start[ires] = start;
        residue[ires] = residue[start];
        start = idAtm;
        ++ires;
      }
    }
  }
  res_start[ires] = start;
  residue.resize(nres);
  geometry.resize(natoms,3);
  atoms.resize(natoms);
  atmtype.resize(natoms);
  chains.resize(ichain,2);
  hetatoms.resize(nheteroatoms);
  hetgeom.resize(nheteroatoms,3);
}

#endif //_File_Format_Package_
