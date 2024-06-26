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

#ifndef _PES_SCAN_METHODS_
#define _PES_SCAN_METHODS_
#include <string>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include "BSet.hpp"
#include "QC.hpp"
#include "math/RandomPackage.hpp"

//descrition:
//object to perform PES scans

class PESscan {
protected:
  std::vector<Molecule> system;             //this is the container for all the molecules to be considered
  QCbasis * method;                         //the quantum chemical method to be used
  BSet basis;                               //basis set object
  bool storeGeom;
  bool center;
  std::string folder;                       //the folder where to store the geometries
  std::string systemname;                   //a name for the system, in case geometries are stored
  std::string header;
  std::string type;                         //type of method used
  double step;
  int nsteps;
  int restart;
  std::vector<double> auxiliary;
  std::vector<double> auxiliary2;
  std::vector<size_t> atomsB;
  std::vector<size_t> atomsAB;
  matrixE auxmat;
  matrixE rotmat;
  matrixE geometryA;
  matrixE geometryB;
  matrixE geometryBfix;
  matrixE geometryAB;
  Molecule B;
  Molecule AB;
public:
  PESscan(QCbasis & _method, const BSet & _basis) {
    method = & _method;
    type = method->Type();
    basis = _basis;
    storeGeom = false;
    step = -1.0;
    nsteps = -1;
    folder = "";
    systemname = "";
    restart = 0;
    method->setRestart(restart);
    center = true;
  }
  ~PESscan() {}
  //getters
  std::vector<Molecule> System() {return system;}
  bool StoreGeometry() {return storeGeom;}
  bool CenterGeometry() {return center;}
  double StepSize() {return step;}
  int NumberSteps() {return nsteps;}
  int Restart() {return restart;}
  std::string StoreFolder() {return folder;}
  std::string SystemName() {return systemname;}
  std::string Header() {return header;}
  //setters
  void setSystem(const std::vector<Molecule> & _system) {system = _system;}
  void setStoreGeometry(bool _storeGeom) {storeGeom = _storeGeom;}
  void setCenterGeometry(bool _center) {center = _center;}
  void setStepSize(double _step) {step = _step;}
  void setNumberSteps(int _nsteps) {nsteps = _nsteps;}
  void setRestart(int newval) {
    restart = newval;
    method->setRestart(restart);
  }
  void setStoreFolder(std::string _folder) {folder = _folder;}
  void setSystemName(std::string _sn) {systemname = _sn;}
  //other functions
  void AddToSystem(const Molecule & _newmol) {system.push_back(_newmol);}
  void ClearSystem() {system.clear();}
  void MoveMolecule(size_t coordinate, double shift, matrixE & geometry) {
    //function that shifts the coordinate of a molecule by a certain shift
    int Natoms = geometry.rows();
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      geometry(idatm + 1,coordinate) += shift;
    }
  }
  void RotateGeometry(int coordinate, double theta, matrixE & geometry) {
    //function that shifts the coordinate of a molecule by a certain angle
    RotateGeomMatrix(coordinate,theta,geometry,rotmat);
  }
  void SetPropertyMatrix(matrixE & PropMat, std::string property, size_t size1, size_t size2 = 1) {
    //function that defines the property matrix
    if ((property == "energy")||(property == "Energy")||(property == "ENERGY")) {
      PropMat.resize(size1,2);
      header = "RAB (A)                 Energy (Hartree)";
    }
    else if ((property == "energyfull")||(property == "EnergyFull")||(property == "ENERGYFULL")) {
      PropMat.resize(size1,size2);
      header = "RAB (A)                 Energy (Hartree)";
    }
    else if ((property == "charges")||(property == "Charges")||(property == "CHARGES")) {
      PropMat.resize(size1,size2 + 1);
      header = "RAB (A)                 Atomic Charges (e)";
    }
    else if ((property == "gradients")||(property == "grad")||(property == "Grad")||(property == "GRAD")) {
      PropMat.resize(size1,3*size2 + 1);
      header = "RAB (A)                 Gradients (Hartree/Angstroem)";
    }
  }
  void CalcProperty(matrixE & PropMat, std::string property, double RAB, size_t position, double propshift = 0.0) {
    //function that evaluates a property and stores it in the property matrix supplied
    if (property != "energyfull") {
      PropMat(position,1) = RAB;
      if (position > 1) {PropMat(position,1) += PropMat(position - 1,1);}
    }
    if ((property == "energy")||(property == "Energy")||(property == "ENERGY")) {
      double en = method->getEnergy(1) - propshift;
      PropMat(position,2) = en;
    }
    else if ((property == "energyfull")||(property == "EnergyFull")||(property == "ENERGYFULL")) {
      for (size_t idp = 0; idp < auxiliary.size(); ++idp) {
        PropMat(idp + 1,position) = auxiliary[idp];
      }
    }
    else if ((property == "charges")||(property == "Charges")||(property == "CHARGES")) {
      if (type != "GFN2") {           //in this case get the Mulliken charges
        method->getMullikenCharges();
      }
      auxiliary = method->getCharges();
      size_t Natoms = auxiliary.size();
      for (size_t idatm = 0; idatm < Natoms; ++idatm) {
        PropMat(position,idatm + 2) = auxiliary[idatm];
      }
    }
    else if ((property == "gradients")||(property == "grad")||(property == "Grad")||(property == "GRAD")) {
      if (type != "GFN2") {method->gEnergy(auxmat,1,0);}
      else {method->gEnergy(auxmat,0,0);}                     //for now only numerical
      size_t sizemat = auxmat.rows();
      for (size_t idatm = 0; idatm < sizemat; ++idatm) {
        PropMat(position,idatm + 2) = auxmat(idatm + 1,1);
      }
    }
  }
  //auxiliary functions
  void SetAndCalculate(Molecule & AB, matrixE & geom, int counter, bool store) {
    //function that sets a certain geometry and calculates the energy
    AB.setGeometry(geom);
    if (store) {AB.WriteXYZ(folder + "/" + systemname,counter);}
    method->setMolecule(AB);
    method->Calculate(0);
  }
  void CenterAndCalc(Molecule & mol) {
    if (center) {mol.CorrectedGeom();}                      //zero the coordinates
    method->setMolecule(mol);
    method->Calculate(0);
  }
  void UpdConvergence(bool minimize, double e_current, double & e_stationarypoint, size_t pos, int & pos_save) {
    //function that updates the convergence conditions
    if (minimize) {        //minimization
      if (e_current < e_stationarypoint) {
        e_stationarypoint = e_current;
        pos_save = pos;
      }
    }
    else {                 //maximization
      if (e_current > e_stationarypoint) {
        e_stationarypoint = e_current;
        pos_save = pos;
      }
    }
  }
  void PrepareMolecule(Molecule & AB, std::vector<size_t> & ABat, int charge, size_t multiplicity) {
    //function preparing the construction of molecule
    AB.setAtoms(ABat);
    AB.setCharge(charge);
    AB.setMultiplicity(multiplicity);
    AB.calcNatoms();
    AB.CountElectrons();
    AB.masses();
  }
  void PrepareBimolecularSystem(int coord, size_t mol1, size_t mol2, double & EA, double & EB, double & EAB, double Rdist, int counter, bool storegeometry) {
    //function that makes the preparation for a calculation
    //get data on mol1
    B = system[mol1 - 1];
    CenterAndCalc(B);
    EA = method->getEnergy(1);
    geometryA = B.Geometry();       //the geometries of the base structures
    //get data on mol2
    B = system[mol2 - 1];
    CenterAndCalc(B);
    EB = method->getEnergy(1);
    geometryB = B.Geometry();       //the geometries of the base structures
    geometryBfix = geometryB;
    atomsB = B.Atoms();
    //get data on megastructure
    geometryAB = geometryA;                        //this is what we are interested in, the composite matrix
    atomsAB = system[mol1 - 1].Atoms();
    ConcatenateV(atomsAB,atomsB);
    int charge = system[mol1 - 1].Charge() + B.Charge();
    size_t multiplicity = system[mol1 - 1].Multiplicity()*B.Multiplicity();
    MoveMolecule(coord,Rdist,geometryB);
    ConcatenateMR(geometryAB,geometryB);
    PrepareMolecule(AB,atomsAB,charge,multiplicity);
    SetAndCalculate(AB,geometryAB,counter,storegeometry);
    EAB = method->getEnergy(1);
  }
  void CalcStep(double delta) {
    //function that calculates step
    if ((nsteps < 0)&&(step > 0.0)) {nsteps = int(delta/step);}
    else if ((nsteps > 0)&&(step < 0.0)) {step = delta/double(nsteps);}
    else if ((nsteps < 0)&&(step < 0.0)) {
      nsteps = 100;
      step = 0.01*delta;
    }
  }
  void RotationCycle(int rcoord, int tcoord, double & EB, double sstep, double Rdist, int counter, bool printg) {
    //function applying one step in the rotation scan
    geometryAB = geometryA;
    RotateGeometry(rcoord,sstep,geometryBfix);
    geometryB = geometryBfix;
    MoveMolecule(tcoord,Rdist,geometryB);
    SetAndCalculate(B,geometryB,1,false);
    EB = method->getEnergy(1);
    ConcatenateMR(geometryAB,geometryB);
    SetAndCalculate(AB,geometryAB,counter,printg);
  }
  void ApplyTranslation(int coordinate, double sstep, double & Rdist, double & EAB, double & EA, double & EB, double & enew) {
    //function that applies a single translation step
    geometryAB = geometryA;
    MoveMolecule(coordinate,sstep,geometryB);
    Rdist += sstep;
    ConcatenateMR(geometryAB,geometryB);
    SetAndCalculate(AB,geometryAB,1,false);
    EAB = method->getEnergy(1);
    enew = EAB - EA - EB;
  }
  void CheckDirection(int coord, double & enew, double & eold, double & dir, double & sstep, double & Rdist, int & chgdir, double & EAB, double EA, double EB) {
    //function that checks the energy to see if translation goes in the direction of minimum
    if (enew > eold) {     //wrong direction
      dir *= -1.0;
      chgdir += 1;
      sstep *= -1.0;
      MoveMolecule(coord,sstep,geometryB);
      Rdist += sstep;
      geometryAB = geometryA;
      ConcatenateMR(geometryAB,geometryB);
      SetAndCalculate(AB,geometryAB,1,false);
      EAB = method->getEnergy(1);
      enew = EAB - EA - EB;
    }
    else {eold = enew;}
  }
  void CheckDirection0(int coord, double & enew, double & eold, double & sstep, double & Rdist, double & EAB, double EA, double EB) {
    //same as above but checks whether the energy goes above zero
    if ((enew < eold)||(enew > 0.0)) {     //wrong direction, so undo
      MoveMolecule(coord,-sstep,geometryB);
      Rdist -= sstep;
      geometryAB = geometryA;
      ConcatenateMR(geometryAB,geometryB);
      SetAndCalculate(AB,geometryAB,1,false);
      EAB = method->getEnergy(1);
      enew = EAB - EA - EB;
      sstep *= 0.1;
    }
    else {eold = enew;}
  }
  //PES functions
  void TwoMoleculeScan(matrixE & PropMat, size_t coordinate, std::string property, double Rmin, double Rmax, double Rfix = 0.0, size_t mol1 = 1, size_t mol2 = 2) {
    //function to scan the PES between two structures
    //coordinate is a 1 based variable, meaning
    //1 -> x
    //2 -> y
    //3 -> z
    //4 -> rotate around x
    //5 -> rotate around y
    //6 -> rotate around z
    //7 -> rotate around all angles
    //PropMat contains the property being calculated, let it be energy, charges, etc, as a function of distance
    int scalefactor = 1;
    double EA;
    double EB;
    double EAB;
    if (system.size() > 1) {
      if ((coordinate > 0)&&(coordinate < 4)) {
        //system preparation
        PrepareBimolecularSystem(coordinate,mol1,mol2,EA,EB,EAB,Rmin,1,storeGeom);
        //getting the number of steps and/or stepsize
        CalcStep(Rmax - Rmin);
        SetPropertyMatrix(PropMat,property,scalefactor*nsteps + 1,atomsAB.size());
        CalcProperty(PropMat,property,Rmin,1,EA + EB);
        for (size_t idstep = 0; idstep < nsteps; ++idstep) {
          geometryAB = geometryA;
          MoveMolecule(coordinate,step,geometryB);
          SetAndCalculate(B,geometryB,idstep + 2,false);
          EB = method->getEnergy(1);
          ConcatenateMR(geometryAB,geometryB);
          SetAndCalculate(AB,geometryAB,idstep + 2,storeGeom);
          CalcProperty(PropMat,property,step,idstep + 2,EA + EB);
        }
      }
      else if ((coordinate > 3)&&(coordinate < 7)) {
        //system preparation
        PrepareBimolecularSystem(coordinate - 3,mol1,mol2,EA,EB,EAB,Rfix,1,storeGeom);
        //getting the number of steps and/or stepsize
        CalcStep(Rmax - Rmin);
        SetPropertyMatrix(PropMat,property,scalefactor*nsteps + 1,atomsAB.size());
        CalcProperty(PropMat,property,0.0,1,EA + EB);
        for (size_t idstep = 0; idstep < nsteps; ++idstep) {
          RotationCycle(coordinate - 3,coordinate - 3,EB,step,Rfix,idstep + 2,storeGeom);
          CalcProperty(PropMat,property,step,idstep + 2,EA + EB);
        }
      }
      else if (coordinate == 7) {
        //molecule 1 is fixed, molecule 2 rotates in all directions
        int counter = 1;
        PrepareBimolecularSystem(1,mol1,mol2,EA,EB,EAB,Rfix,counter,storeGeom);
        ++counter;
        matrixE geometryBx = geometryBfix;
        matrixE geometryBy = geometryBfix;
        matrixE geometryBz = geometryBfix;
        //getting the number of steps and/or stepsize
        CalcStep(Rmax - Rmin);
        scalefactor = nsteps*nsteps;
        SetPropertyMatrix(PropMat,property,scalefactor*nsteps + 1,atomsAB.size());
        CalcProperty(PropMat,property,0.0,1,EA + EB);
        for (size_t idx = 0; idx < nsteps; ++idx) {
          RotateGeometry(1,step,geometryBx);
          geometryBy = geometryBx;
          for (size_t idy = 0; idy < nsteps; ++idy) {
            RotateGeometry(2,step,geometryBy);
            geometryBz = geometryBy;
            for (size_t idz = 0; idz < nsteps; ++idz, ++counter) {
              RotateGeometry(3,step,geometryBz);
              geometryB = geometryBz;
              MoveMolecule(1,Rfix,geometryB);
              SetAndCalculate(B,geometryB,counter,false);
              EB = method->getEnergy(1);
              geometryAB = geometryA;
              ConcatenateMR(geometryAB,geometryB);
              SetAndCalculate(AB,geometryAB,counter,storeGeom);
              CalcProperty(PropMat,property,step,counter,EA + EB);
            }
          }
        }
      }
    }
    else {std::cout << "WARNING: PES.hpp - PESscan - TwoMoleculeScan: insufficient molecules for two molecule scan" << std::endl;}
  }
  void TwoMolecule3DScan(matrixE & PropMat, size_t coordinate, double Rmin, double Rmax, double theta_min, double theta_max, size_t ntheta_step, size_t mol1 = 1, size_t mol2 = 2) {
    //function to perform a full three dimensional scan of the PES between two structures
    //it basically takes a bimolecular structure, rotates it and after each rotation the PES is scanned along intermolecular distance
    //the user must define the direction used for the scan, which is defined by coordinate, a 1 based variable, meaning
    //1 -> x
    //2 -> y
    //3 -> z
    //PropMat contains the property being calculated, which can only be energy
    if (system.size() > 1) {
      //molecule 1 is fixed, molecule 2 rotates in all directions
      int counter = 1;
      int counterg = 1;
      double EA;
      double EB;
      double EAB;
      double thetastep = (theta_max - theta_min)/double(ntheta_step - 1);
      double rotx = theta_min;
      double roty = theta_min;
      double rotz = theta_min;
      std::vector<double> rotation(ntheta_step,thetastep);
      rotation[0] = 0.0;
      PrepareBimolecularSystem(1,mol1,mol2,EA,EB,EAB,Rmin,counter,false);
      ++counter;
      matrixE geometryBx = geometryBfix;
      matrixE geometryBy = geometryBfix;
      matrixE geometryBz = geometryBfix;
      //getting the number of steps and/or stepsize
      CalcStep(Rmax - Rmin);
      auxiliary.resize(nsteps + 1);
      auxiliary2.resize(nsteps + 1);
      SetPropertyMatrix(PropMat,"energyfull",nsteps + 1,ntheta_step*ntheta_step*ntheta_step + 1);
      //CalcProperty(PropMat,"energyfull",0.0,1,EA + EB);
      for (size_t idx = 0; idx < ntheta_step; ++idx) {
        RotateGeometry(1,rotation[idx],geometryBx);
        rotx += rotation[idx];
        geometryBy = geometryBx;
        roty = theta_min;
        for (size_t idy = 0; idy < ntheta_step; ++idy) {
          RotateGeometry(2,rotation[idy],geometryBy);
          roty += rotation[idy];
          geometryBz = geometryBy;
          rotz = theta_min;
          for (size_t idz = 0; idz < ntheta_step; ++idz, ++counter) {
            RotateGeometry(3,rotation[idz],geometryBz);
            rotz += rotation[idz];
            geometryB = geometryBz;
            MoveMolecule(coordinate,Rmin,geometryB);
            SetAndCalculate(B,geometryB,counterg,false);
            EB = method->getEnergy(1);
            geometryAB = geometryA;
            ConcatenateMR(geometryAB,geometryB);
            SetAndCalculate(AB,geometryAB,counterg,storeGeom);
            ++counterg;
            auxiliary[0] = method->getEnergy(1) - EA - EB;
            auxiliary2[0] = Rmin;
            //now scan distances
            for (size_t idstep = 0; idstep < nsteps; ++idstep, ++counterg) {
              geometryAB = geometryA;
              MoveMolecule(coordinate,step,geometryB);
              auxiliary2[idstep + 1] = auxiliary2[idstep] + step;
              SetAndCalculate(B,geometryB,counterg,false);
              EB = method->getEnergy(1);
              ConcatenateMR(geometryAB,geometryB);
              SetAndCalculate(AB,geometryAB,counterg,storeGeom);
              auxiliary[idstep + 1] = method->getEnergy(1) - EA - EB;
            }
            CalcProperty(PropMat,"energyfull",step,counter,0.0);
          }
        }
      }
      //now fix the propmat
      for (size_t idp = 0; idp < nsteps + 1; ++idp) {
        PropMat(idp + 1,1) = auxiliary2[idp];
      }
    }
    else {std::cout << "WARNING: PES.hpp - PESscan - TwoMolecule3DScan: insufficient molecules for two molecule scan" << std::endl;}
  }
  void OptOrientation(matrixE & GeomA, matrixE & GeomB, double Rfix = 4.0, size_t mol1 = 1, size_t mol2 = 2, bool minimize = true, bool storeg = true, double thresholdzero = 1.0e-8, int maxiter = 200) {
    //function that finds the closest orientation between two molecules that minimizes/maximizes interaction energy
    //storeg stores the geometries in the respective system's molecules
    double EA;
    double EB;
    double EAB;
    PrepareBimolecularSystem(1,mol1,mol2,EA,EB,EAB,Rfix,1,false);
    //getting the number of steps and/or stepsize
    CalcStep(2.0*pi);
    double enew = EAB - EA - EB;
    double eold = enew + 100.0;
    double eaux = enew;
    double estationary = enew;
    int position;
    for (size_t iter = 0; iter < maxiter; ++iter) {
      if (fabs(enew - eold) < thresholdzero) {
        std::cout << "minimum in energy found: " << enew << "\n";
        break;
      }
      std::cout << "current energy minimum: " << enew << "\n";
      eold = enew;
      for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
        for (size_t idstep = 0; idstep < nsteps; ++idstep) {
          RotationCycle(idcoord,1,EB,step,Rfix,1,false);
          EAB = method->getEnergy(1);
          eaux = EAB - EA - EB;
          UpdConvergence(true,eaux,estationary,idstep + 1,position);
        }
        RotateGeometry(idcoord,(position - nsteps)*step,geometryBfix);
        enew = estationary;
        position = 0;
      }
    }
    GeomA = geometryA;
    GeomB = geometryB;
    if (storeg) {
      system[mol1 - 1].setGeometry(geometryA);
      system[mol2 - 1].setGeometry(geometryB);
    }
  }
  void OptDistance(matrixE & GeomA, matrixE & GeomB, int coordinate, double Rbegin = 2.0, size_t mol1 = 1, size_t mol2 = 2, double stepsz = 0.01, double stepmin = 0.001, bool storeg = true, double thresholdzero = 1.0e-8, int maxiter = 200) {
    //function that finds the distance between two molecules that minimizes interaction energy
    //storeg stores the geometries in the respective system's molecules
    double EA;
    double EB;
    double EAB;
    PrepareBimolecularSystem(coordinate,mol1,mol2,EA,EB,EAB,Rbegin,1,false);
    //getting the number of steps and/or stepsize
    double enew = EAB - EA - EB;
    double eold = enew;
    double direction = 1.0;
    double delta = (stepsz - stepmin)/3.0;
    double step1 = stepsz;
    double step2 = stepsz - delta;
    double step3 = step2 - delta;
    double step4 = stepmin;
    double cstep = step1;
    double auxie;
    double currentdistance = Rbegin;
    int countchgdir = 0;
    int position;
    bool converged = false;
    for (size_t iter = 0; iter < maxiter; ++iter) {
      ApplyTranslation(coordinate,cstep,currentdistance,EAB,EA,EB,enew);
      if ((fabs(enew - eold) < thresholdzero)||(countchgdir == 5)) {
        auxie = fabs(cstep);
        if (fabs(auxie - step1) < thresholdzero) {cstep = direction*step2;}
        else if (fabs(auxie - step2) < thresholdzero) {cstep = direction*step3;}
        else if (fabs(auxie - step3) < thresholdzero) {cstep = direction*step4;}
        else if (fabs(auxie - step4) < thresholdzero) {
          converged = true;
          break;
        }
        countchgdir = 0;
        continue;
      }
      CheckDirection(coordinate,enew,eold,direction,cstep,currentdistance,countchgdir,EAB,EA,EB);
      std::cout << "current step: rAB = " << currentdistance << "   EAB = " << enew << std::endl;
    }
    if (converged) {std::cout << "minimum estimation found: " << currentdistance << " " << enew << std::endl;}
    else {std::cout << "no minimum found: " << currentdistance << " " << enew << std::endl;}
    GeomA = geometryA;
    GeomB = geometryB;
    if (storeg) {
      system[mol1 - 1].setGeometry(geometryA);
      system[mol2 - 1].setGeometry(geometryB);
    }
  }
  matrixE FindGlobalMinimum(bool fullcharacterization, matrixE & GeomA, matrixE & GeomB, double Rbegin = 2.0, size_t mol1 = 1, size_t mol2 = 2, double steptrans = 0.01, bool storeg = true, double thresholdzero = 1.0e-8, int maxiter = 200) {
    //function that finds the global minimum on a two-molecule PES
    double EA;
    double EB;
    double EAB;
    PrepareBimolecularSystem(1,mol1,mol2,EA,EB,EAB,Rbegin,1,false);
    //getting the number of steps and/or stepsize
    CalcStep(2.0*pi);
    double enew = EAB - EA - EB;
    double eold = enew + 100.0;
    double eaux = enew;
    double estationary = enew;
    double direction = 1.0;
    double delta = (steptrans - 0.001*steptrans)/3.0;
    double step1 = steptrans;
    double step2 = steptrans - delta;
    double step3 = step2 - delta;
    double step4 = 0.001*steptrans;
    double cstep = steptrans;
    double auxie;
    double currentdistance = Rbegin;
    double actuaElmin = enew;
    double actuaRlmin = currentdistance;
    int countchgdir = 0;
    int position = 0;
    bool converged = false;
    std::cout << std::setprecision(10);
    for (size_t iter = 0; iter < maxiter; ++iter) {
      if (fabs(enew - eold) < thresholdzero) {
        std::cout << "minimum in energy found: " << enew << "\n";
        break;
      }
      std::cout << "iteration " << iter << "\n";
      std::cout << "  current minimum: " << currentdistance << " " << enew << "\n";
      eold = enew;
      estationary = enew;
      for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
        for (size_t idstep = 0; idstep < nsteps; ++idstep) {
          RotationCycle(idcoord,1,EB,step,currentdistance,1,false);
          EAB = method->getEnergy(1);
          eaux = EAB - EA - EB;
          UpdConvergence(true,eaux,estationary,idstep + 1,position);
        }
        RotateGeometry(idcoord,(position - nsteps)*step,geometryBfix);
        geometryB = geometryBfix;
        MoveMolecule(1,currentdistance,geometryB);
        SetAndCalculate(B,geometryB,1,false);
        EB = method->getEnergy(1);
        geometryAB = geometryA;
        ConcatenateMR(geometryAB,geometryB);
        SetAndCalculate(AB,geometryAB,1,false);
        EAB = method->getEnergy(1);
        enew = EAB - EA - EB;
        position = 0;
      }
      //with orientation optimized, now we optimize distance again
      cstep = steptrans;
      direction = 1.0;
      converged = false;
      actuaElmin = enew;
      cstep = steptrans;
      countchgdir = 0;
      for (size_t itert = 0; itert < maxiter; ++itert) {
        ApplyTranslation(1,cstep,currentdistance,EAB,EA,EB,enew);
        if (enew < actuaElmin) {
          actuaElmin = enew;
          actuaRlmin = currentdistance;
        }
        if ((fabs(enew - eaux) < thresholdzero)||(countchgdir == 5)) {
          auxie = fabs(cstep);
          if (fabs(auxie - step1) < thresholdzero) {cstep = direction*step2;}
          else if (fabs(auxie - step2) < thresholdzero) {cstep = direction*step3;}
          else if (fabs(auxie - step3) < thresholdzero) {cstep = direction*step4;}
          else if (fabs(auxie - step4) < thresholdzero) {
            converged = true;
            break;
          }
          MoveMolecule(1,-direction*auxie,geometryB);
          currentdistance -= direction*auxie;
          countchgdir = 0;
          continue;
        }
        CheckDirection(1,enew,eaux,direction,cstep,currentdistance,countchgdir,EAB,EA,EB);
      }
      if (actuaElmin < enew) {
        geometryB = geometryBfix;
        MoveMolecule(1,actuaRlmin,geometryB);
        currentdistance = actuaRlmin;
        geometryAB = geometryA;
        ConcatenateMR(geometryAB,geometryB);
        enew = actuaElmin;
      }
      SetAndCalculate(AB,geometryAB,iter + 1,true);
    }
    GeomA = geometryA;
    GeomB = geometryB;
    if (storeg) {
      system[mol1 - 1].setGeometry(geometryA);
      system[mol2 - 1].setGeometry(geometryB);
    }
    matrixE results(1 + 2*fullcharacterization,2);
    results(fullcharacterization + 1,1) = currentdistance;
    results(fullcharacterization + 1,2) = enew;
    if (fullcharacterization) {
      //find sigma
      std::cout << "\n\nfinding sigma \n";
      cstep = -steptrans;
      for (size_t iter = 0; iter < 5*maxiter; ++iter) {
        ApplyTranslation(1,cstep,currentdistance,EAB,EA,EB,enew);
        if (fabs(enew) < thresholdzero) {
          converged = true;
          break;
        }
        CheckDirection0(1,enew,eold,cstep,currentdistance,EAB,EA,EB);
        std::cout << "current step: rAB = " << currentdistance << "   EAB = " << enew << "\n";
      }
      results(1,1) = currentdistance;
      results(1,2) = enew;
      cstep = steptrans;
      //go back to minimum
      MoveMolecule(1,results(2,1) - results(1,1),geometryB);
      currentdistance += results(2,1) - results(1,1);
      std::cout << "\nfinding distance at which potential zeroes \n";
      eold = results(2,2);
      for (size_t iter = 0; iter < 5*maxiter; ++iter) {
        ApplyTranslation(1,cstep,currentdistance,EAB,EA,EB,enew);
        if (fabs(enew) < thresholdzero) {
          converged = true;
          break;
        }
        CheckDirection0(1,enew,eold,cstep,currentdistance,EAB,EA,EB);
        std::cout << "current step: rAB = " << currentdistance << "   EAB = " << enew << "\n";
      }
      results(3,1) = currentdistance;
      results(3,2) = enew;
    }
    return results;
  }
};
class ProcessStructure {
protected:
  Molecule system;
  std::vector<size_t> atomsSystem;
  std::vector<size_t> atomsSolvent;
  matrixE rotmat;
  matrixE geomSystem;
  matrixE geomSolvent;
  matrixE positions;
public:
  ProcessStructure(Molecule & initialize) {
    system = initialize;
    srand(time(NULL));
  }
  ~ProcessStructure() {}
  //getters
  Molecule getSystem() {return system;}
  //setters
  void setSystem(Molecule & newsys) {system = newsys;}
  //other functions
  bool PairCollision(double xn, double yn, double zn, double xo, double yo, double zo, double rn, double ro) {
    //function to check whether two particles are colliding
    double distance = sqrt((xn - xo)*(xn - xo) + (yn - yo)*(yn - yo) + (zn - zo)*(zn - zo));
    return (distance < rn + ro);
  }
  //preparation functions
  void AddSolventRandomly(Molecule & solvmol, int nsolvent, double boxx, double boxy, double boxz, bool ellipsoid = false, double factor = 1.4) {
    //function that adds randomly nsolvent molecules of a solvmol to the system
    solvmol.CorrectedGeom();
    system.CorrectedGeom();
    geomSystem = system.Geometry();
    geomSolvent = solvmol.Geometry();
    atomsSystem = system.Atoms();
    atomsSolvent = solvmol.Atoms();
    int chgSystem = system.Charge();
    int chgSolvent = solvmol.Charge();
    size_t multSystem = system.Multiplicity();
    size_t multSolvent = solvmol.Multiplicity();
    //get max distance between solvent atoms
    int NAtoms_solv = solvmol.Natoms();
    int atomA;
    int atomB;
    int collidingparticle;
    std::vector<double> translation(3);
    double radius_solvent = 0.0;
    double radius_solute = 0.0;
    double distance = 0.0;
    double thetax;
    double thetay;
    double thetaz;
    double xmin = -0.5*boxx;
    double ymin = -0.5*boxy;
    double zmin = -0.5*boxz;
    double xmax = 0.5*boxx;
    double ymax = 0.5*boxy;
    double zmax = 0.5*boxz;
    double auxx;
    double auxy;
    double auxz;
    bool collision;
    positions = solvmol.Geometry();
    for (size_t idAtm = 0; idAtm < NAtoms_solv; ++idAtm) {
      for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
        distance = Distance(idAtm + 1,idBtm + 1,positions);
        if (distance > radius_solvent) {
          radius_solvent = distance;
          atomA = idAtm;
          atomB = idBtm;
        }
      }
    }
    //add contribution from atomic radii to avoid collisions
    radius_solvent += factor*(AtmRadii(solvmol.Atom(atomA + 1)) + AtmRadii(solvmol.Atom(atomB + 1)));
    //now for solute
    int NAtoms_solute = system.Natoms();
    int NAtoms_solute2 = NAtoms_solute;
    positions = system.Geometry();
    for (size_t idAtm = 0; idAtm < NAtoms_solute; ++idAtm) {
      for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
        distance = Distance(idAtm + 1,idBtm + 1,positions);
        if (distance > radius_solute) {
          radius_solute = distance;
          atomA = idAtm;
          atomB = idBtm;
        }
      }
    }
    //add contribution from atomic radii to avoid collisions
    radius_solute += factor*(AtmRadii(solvmol.Atom(atomA + 1)) + AtmRadii(solvmol.Atom(atomB + 1)));
    //now the radii 
    radius_solvent *= 0.5;
    radius_solute *= 0.5;
    positions.resize(nsolvent,3);       //to add the position of the new solvent molecules
    //now add randomly the solute at random orientations
    for (size_t idsolv = 0; idsolv < nsolvent; ++idsolv) {
      translation[0] = fRandom(xmin,xmax);
      translation[1] = fRandom(ymin,ymax);
      translation[2] = fRandom(zmin,zmax);
      if (ellipsoid) {                  //in this case apply the ellipsoid conditions to determine whether to reject the step
        auxx = translation[0]/boxx;
        auxy = translation[1]/boxy;
        auxz = translation[2]/boxz;
        if (auxx*auxx + auxy*auxy + auxz*auxz >= 1.0) {       //reject step
          --idsolv;
          continue;
        }
      }
      thetax = fRandom(0.0,2.0*pi);
      thetay = fRandom(0.0,2.0*pi);
      thetaz = fRandom(0.0,2.0*pi);
      geomSolvent = solvmol.Geometry();
      RotateGeomMatrix(1,thetax,geomSolvent,rotmat);
      RotateGeomMatrix(2,thetay,geomSolvent,rotmat);
      RotateGeomMatrix(3,thetaz,geomSolvent,rotmat);
      for (size_t idAtm = 0; idAtm < NAtoms_solv; ++idAtm) {
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          geomSolvent(idAtm + 1,idcoord + 1) += translation[idcoord];
        }
      }
      collision = PairCollision(translation[0],translation[1],translation[2],0.0,0.0,0.0,radius_solvent,radius_solute);
      if (collision) {                           //collision with solute? Check atom by atom
        collision = false;
        for (size_t idAsolvent = 0; idAsolvent < NAtoms_solv; ++idAsolvent) {
          for (size_t idAsolute = 0; idAsolute < NAtoms_solute; ++idAsolute) {
            auxx = geomSystem(idAsolute + 1,1) - geomSolvent(idAsolvent + 1,1);
            auxy = geomSystem(idAsolute + 1,2) - geomSolvent(idAsolvent + 1,2);
            auxz = geomSystem(idAsolute + 1,3) - geomSolvent(idAsolvent + 1,3);
            distance = sqrt(auxx*auxx + auxy*auxy + auxz*auxz);
            if (distance < factor*(AtmRadii(atomsSolvent[idAsolvent]) + AtmRadii(atomsSystem[idAsolute]))) {
              collision = true;
              break;
            }
          }
          if (collision) {break;}
        }
      }
      if (!collision) {                                     //collision with another solvent molecule?
        collision = false;
        collidingparticle = -1;
        for (size_t idsolv2 = 0; idsolv2 < idsolv; ++idsolv2) {
          collision = PairCollision(translation[0],translation[1],translation[2],positions(idsolv2 + 1,1),positions(idsolv2 + 1,2),positions(idsolv2 + 1,3),radius_solvent,radius_solvent);
          if (collision) {
            collidingparticle = idsolv2;
            break;
          }
        }
        if (collidingparticle > -1) {
          collision = false;
          for (size_t idAsolvent = 0; idAsolvent < NAtoms_solv; ++idAsolvent) {
            for (size_t idAsolute = 0; idAsolute < NAtoms_solute2; ++idAsolute) {
              auxx = geomSystem(idAsolute + 1,1) - geomSolvent(idAsolvent + 1,1);
              auxy = geomSystem(idAsolute + 1,2) - geomSolvent(idAsolvent + 1,2);
              auxz = geomSystem(idAsolute + 1,3) - geomSolvent(idAsolvent + 1,3);
              distance = sqrt(auxx*auxx + auxy*auxy + auxz*auxz);
              if (distance < factor*(AtmRadii(atomsSolvent[idAsolvent]) + AtmRadii(atomsSystem[idAsolute]))) {
                collision = true;
                break;
              }
            }
            if (collision) {break;}
          }
          if (collision) {
            --idsolv;
            continue;
          }
        }
        //if here then placement was accepted
        positions(idsolv + 1,1) = translation[0];
        positions(idsolv + 1,2) = translation[1];
        positions(idsolv + 1,3) = translation[2];
        ConcatenateMR(geomSystem,geomSolvent);
        ConcatenateV(atomsSystem,atomsSolvent);
        chgSystem += chgSolvent;
        multSystem += multSolvent;
        system.setAtoms(atomsSystem);
        system.setGeometry(geomSystem);
        system.setCharge(chgSystem);
        system.setMultiplicity(multSystem);
        system.calcNatoms();
        system.CountElectrons();
        system.masses();
        NAtoms_solute2 += NAtoms_solv;
      }
      else {--idsolv;}
    }
    system.WriteXYZ("megastruct",0,10);
  }
};

#endif //_PES_SCAN_METHODS_
