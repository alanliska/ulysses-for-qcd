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

#ifndef _Molecular_Dynamics_
#define _Molecular_Dynamics_
#include <vector>
#include <string>
#include <math.h>
#include <time.h>
#include "QC.hpp"
#include "BSet.hpp"
#include "Molecule.hpp"
#include "Potentials.hpp"
#include "ConstantsPackage.hpp"
#include "atoms/AtomPackage.hpp"
#include "math/MatrixPackage.hpp"
#include "math/RandomPackage.hpp"
#include "math/SolverPackage.hpp"

//description:
//molecular dynamics package

std::vector<int> getNumbersFromFile(std::string filename) {
  //function reading a series of comma separated numbers from a file
  std::vector<int> numbers;
  std::ifstream file(filename);
  if (!file.is_open()) {throw std::string("ERROR: MolecularDynamics.hpp: getNumbersFromFile(): file cannot be open");}
  std::string line;
  while (std::getline(file,line)) {
      std::stringstream ss(line);
      std::string numberStr;
      while (std::getline(ss, numberStr,',')) {
          try {
            int number = std::stoi(numberStr);
            numbers.push_back(number);
          } 
          catch (std::invalid_argument & e) {std::cout << "WARNING: MolecularDynamics.hpp: getNumbersFromFile(): skipping invalid number" << std::endl;}
      }
  }
  file.close();
  return numbers;
}

class MetaDynamics {
  //class for metadynamics objects, based on
  //S. Grimme, J. Chem. Theory Comput., 15, 5, 2847, 2019
protected:
  matrixE gauxREF;
  matrixE gaux0;
  matrixE Rmat;
  matrixE Smat;
  matrixE gRMSD;                                          //RMSD gradients
  std::vector<matrixE> metaset;                           //meta geometries
  double alpha;                                           //general alpha to use everywhere
  double kpush;                                           //general k for pushing
  double kpull;                                           //general k for pulling
  double ppull;                                           //pull force
  double qi[4];                                           //auxiliary
  size_t nmtdstructs;                                     //number of meta-structures (size of metaset)
  std::vector<size_t> metaatoms;                          //list of atoms for which metadynamics applies
  std::vector<double> mtdfactors;                         //kfactors for all declared structures
  std::vector<double> mtdalphas;                          //alphas for all declared structures
  std::vector<double> auxV;
public:
  MetaDynamics(double alp = 1.2, double push = 0.003, double pull = -0.015) {
    setOptions(push,pull,0.05,alp);
    Rmat.resize(3,3);
    Smat.resize(4,4);
  }
  ~MetaDynamics() {}
  //getters
  std::vector<size_t> MetaAtoms() {return metaatoms;}
  std::vector<double> MetaFactors() {return mtdfactors;}
  std::vector<double> MetaAlphas() {return mtdalphas;}
  double Alpha() {return alpha;}
  double kPush() {return kpush;}
  double kPull() {return kpull;}
  double pPull() {return ppull;}
  std::vector<matrixE> MetaSet() {return metaset;}
  matrixE MetaGeometry(size_t istruc) {return metaset[istruc - 1];}
  size_t TotalMetaAtoms() {return metaatoms.size();}
  //setters
  void setOptions(double push = 0.003, double pull = -0.015, double p_pull = 0.05, double alp = 1.2) {
    //general object configurator
    kpush = push;
    kpull = pull;
    alpha = alp;
    ppull = p_pull;
    srand(time(NULL));
    std::cout << "temporarily using fixed seed\n";
    srand(100);
  }
  void setAlpha(double alp) {alpha = alp;}
  void setkPush(double kappa) {kpush = kappa;}
  void setkPull(double kappa) {kpull = kappa;}
  void setpPull(double pparam) {ppull = pparam;}
  void setMetaFactors(const std::vector<double> & newfactors) {mtdfactors = newfactors;}
  void setMetaAlphas(const std::vector<double> & newfactors) {mtdalphas = newfactors;}
  void setMetaAlphas2Alpha(double newalpha = 0.0) {
    if (fabs(newalpha) > 1.0e-12) {alpha = newalpha;}
    for (size_t istruct = 0; istruct < nmtdstructs; ++istruct) {
      mtdalphas[istruct] = alpha;
    }
  }
  void setMetaFactors2kPush(double newk = 0.0) {
    if (fabs(newk) > 1.0e-12) {kpush = newk;}
    for (size_t istruct = 0; istruct < nmtdstructs; ++istruct) {
      mtdfactors[istruct] = kpush;
    }
  }
  void setMetaFactors2kPull(double newk = 0.0) {
    if (fabs(newk) > 1.0e-12) {kpull = newk;}
    for (size_t istruct = 0; istruct < nmtdstructs; ++istruct) {
      mtdfactors[istruct] = kpull;
    }
  }
  //functions for setting the metaset
  void setMetaSet(const std::vector<matrixE> & newset) {
    metaset = newset;
    nmtdstructs = metaset.size();
    mtdfactors.resize(nmtdstructs);
    mtdalphas.resize(nmtdstructs);
    for (size_t idx = 0; idx < nmtdstructs; ++idx) {
      mtdfactors[idx] = kpush;
      mtdalphas[idx] = alpha;
    }
  }
  void setMetaSet2Structure(const matrixE & newgeom) {
    metaset.resize(1);
    metaset[0] = newgeom;
    nmtdstructs = 1;
    mtdfactors.resize(1);
    mtdfactors[0] = kpush;
    mtdalphas.resize(1);
    mtdalphas[0] = alpha;
  }
  void addMetaStructure(const matrixE & newgeom) {
    metaset.push_back(newgeom);
    nmtdstructs = metaset.size();
    mtdfactors.resize(nmtdstructs + 1);
    mtdalphas.resize(nmtdstructs + 1);
    mtdfactors[nmtdstructs] = kpush;
    mtdalphas[nmtdstructs] = alpha;
  }
  void clearMetaSet() {
    metaset.clear();
    mtdfactors.clear();
    mtdalphas.clear();
    nmtdstructs = 0;
  }
  //functions to define list of meta-atoms; these functions are based on index 1 (input), so internally shifted by one
  void setMetaAtoms(const std::vector<size_t> & atmmeta, bool shift = true) {
    //function that takes a vector of atomic indices to add to the meta-atom list
    metaatoms = atmmeta;
    size_t nconstr = metaatoms.size();
    gRMSD.resize(nconstr,3);
    gauxREF.resize(nconstr,3);
    gaux0.resize(nconstr,3);
    if (shift) {
      for (size_t idx = 0; idx < nconstr; ++idx) {
        --metaatoms[idx];
      }
    }
  }
  void addMetaAtom(size_t atmidx) {
    //function that adds atom atmidx - 1 to the meta-atoms list
    metaatoms.push_back(atmidx - 1);
    gRMSD.resize(metaatoms.size(),3);
    gauxREF.resize(metaatoms.size(),3);
    gaux0.resize(metaatoms.size(),3);
  }
  void addMetaElements(const std::string & atmsymbol, const std::vector<size_t> & atoms) {
    //function that adds specific elements (e.g., all protons) to meta-atoms
    int elnumber = Symbol2AtomNr(atmsymbol);
    size_t Natoms = atoms.size();
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      if (atoms[idAtm] == elnumber) {metaatoms.push_back(idAtm);}
    }
    gRMSD.resize(metaatoms.size(),3);
    gauxREF.resize(metaatoms.size(),3);
    gaux0.resize(metaatoms.size(),3);
  }
  void addMetaRange(size_t init, size_t end) {
    //adds atoms between positions init and end, inclusively, to the meta-atom list
    --init;
    int delta = int(end) - int(init);
    int currentsize = metaatoms.size();
    metaatoms.resize(currentsize + delta);
    for (size_t idpos = 0; idpos < delta; ++idpos) {
      metaatoms[currentsize + idpos] = init + idpos;
    }
    gRMSD.resize(currentsize + delta,3);
    gauxREF.resize(currentsize + delta,3);
    gaux0.resize(currentsize + delta,3);
  }
  void clearMetaAtoms() {
    //function that clears the meta-atom list
    metaatoms.clear();
    gRMSD.clear();
    gauxREF.clear();
    gaux0.clear();
  }
  //RMSD functions
  void Quaternion2Rotation(matrixE & Quat, matrixE & Rot, int pos) {
    //function that converts a quaternion in position pos of matrix Quat into a rotation matrix in Rot
    qi[0] = Quat(1,pos);
    qi[1] = Quat(2,pos);
    qi[2] = Quat(3,pos);
    qi[3] = Quat(4,pos);
    for (size_t idr = 1; idr < 5; ++idr) {
      for (size_t idc = idr; idc < 5; ++idc) {
        Quat(idr,idc) = 2.0*qi[idr - 1]*qi[idc - 1];
      }
    }
    Quat(1,1) -= 1.0;
    Rot(1,1) = Quat(1,1) + Quat(2,2);
    Rot(1,2) = Quat(2,3) - Quat(1,4);
    Rot(1,3) = Quat(2,4) + Quat(1,3);
    Rot(2,1) = Quat(2,3) + Quat(1,4);
    Rot(2,2) = Quat(1,1) + Quat(3,3);
    Rot(2,3) = Quat(3,4) - Quat(1,2);
    Rot(3,1) = Quat(2,4) - Quat(1,3);
    Rot(3,2) = Quat(3,4) + Quat(1,2);
    Rot(3,3) = Quat(1,1) + Quat(4,4);
  }
  double RMSD(matrixE & geomA, matrixE & geomB, matrixE & grad) {
    //calculation of least square root-mean-square-deviation between two sets of points, using a quaternion based-method
    double centerA[3];
    double centerB[3];
    double anorm = 0.0;
    double bnorm = 0.0;
    double ermsd;
    size_t nconstraint = geomA.rows();
    //get barycenters, centroidal coordinates
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      centerA[idcoord] = 0.0;
      centerB[idcoord] = 0.0;
      for (size_t idAtm = 0; idAtm < nconstraint; ++idAtm) {
        centerA[idcoord] += geomA(idAtm + 1,idcoord + 1);
        centerB[idcoord] += geomB(idAtm + 1,idcoord + 1);
      }
      centerA[idcoord] /= double(nconstraint);
      centerB[idcoord] /= double(nconstraint);
    }
    //norms
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      for (size_t idAtm = 0; idAtm < nconstraint; ++idAtm) {
        geomA(idAtm + 1,idcoord + 1) -= centerA[idcoord];
        anorm += geomA(idAtm + 1,idcoord + 1)*geomA(idAtm + 1,idcoord + 1);
        geomB(idAtm + 1,idcoord + 1) -= centerB[idcoord];
        bnorm += geomB(idAtm + 1,idcoord + 1)*geomB(idAtm + 1,idcoord + 1);
      }
    }
    //Rmatrix
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      for (size_t jdcoord = 0; jdcoord < 3; ++jdcoord) {
        Rmat(idcoord + 1,jdcoord + 1) = 0.0;
        for (size_t idAtm = 0; idAtm < nconstraint; ++idAtm) {
          Rmat(idcoord + 1,jdcoord + 1) += geomA(idAtm + 1,idcoord + 1)*geomB(idAtm + 1,jdcoord + 1);
        }
      }
    }
    //s matrix
    Smat(1,1) = Rmat(1,1) + Rmat(2,2) + Rmat(3,3);
    Smat(2,1) = Rmat(2,3) - Rmat(3,2);
    Smat(3,1) = Rmat(3,1) - Rmat(1,3);
    Smat(4,1) = Rmat(1,2) - Rmat(2,1);
    Smat(1,2) = Smat(2,1);
    Smat(2,2) = Rmat(1,1) - Rmat(2,2) - Rmat(3,3);
    Smat(3,2) = Rmat(1,2) + Rmat(2,1);
    Smat(4,2) = Rmat(1,3) + Rmat(3,1);
    Smat(1,3) = Smat(3,1);
    Smat(2,3) = Smat(3,2);
    Smat(3,3) = -Rmat(1,1) + Rmat(2,2) - Rmat(3,3);
    Smat(4,3) = Rmat(2,3) + Rmat(3,2);
    Smat(1,4) = Smat(4,1);
    Smat(2,4) = Smat(4,2);
    Smat(3,4) = Smat(4,3);
    Smat(4,4) = -Rmat(1,1) - Rmat(2,2) + Rmat(3,3);
    //get eigenvalues
    auxV = MatDiag(Smat);
    //convert quaternion q to rotation matrix U
    Quaternion2Rotation(Smat,Rmat,4);
    //RMSD + small number to avoid division by 0
    ermsd = sqrt(fmax(0.0,(anorm + bnorm - 2.0*auxV[3]))/double(nconstraint)) + 1.e-9;
    //now gradients
    anorm = dist_Angstrom2aum1/(ermsd*double(nconstraint));
    for (size_t idAtm = 0; idAtm < nconstraint; ++idAtm) {
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        centerA[idcoord] = 0.0;
        for (size_t jdcoord = 0; jdcoord < 3; ++jdcoord) {
          centerA[idcoord] += Rmat(jdcoord + 1,idcoord + 1)*geomB(idAtm + 1,jdcoord + 1);
        }
        grad(idAtm + 1,idcoord + 1) = anorm*(geomA(idAtm + 1,idcoord + 1) - centerA[idcoord]);
      }
    }
    return ermsd;
  }
  void MTDPrepare(const matrixE & geometry) {
    //function to prepare atom constraint if user did not apply any constraint specifically
    size_t natm = geometry.rows();
    size_t nsets = metaset.size();
    if (metaatoms.size() == 0) {
      metaatoms.resize(natm);
      for (size_t idx = 0; idx < natm; ++idx) {
        metaatoms[idx] = idx;
      }
    }
    //else -> user defined it
    gauxREF.resize(metaatoms.size(),3);
    gaux0.resize(metaatoms.size(),3);
    metaset[0] = geometry;
    for (size_t idset = 1; idset < nsets; ++idset) {
      metaset[idset].resize(1,1);
    }
    //randomize structure
    Randomize(metaset[0],1.0e-6);
    gRMSD.resize(natm,3);
  }
  //the main function
  double Metadynamic(const matrixE & geometry, matrixE & gradients) {
    double Ebias = 0.0;
    nmtdstructs = metaset.size();
    if (nmtdstructs > 0) {
      double rmsd;
      double Ei;
      size_t masz = metaatoms.size();
      size_t idAtm;
      size_t icnt;
      for (size_t istruct = 0; istruct < nmtdstructs; ++istruct) {
        //build constrained geometry and the one at stake
        if ((istruct > 0)&&(metaset[istruct].rows() != metaset[istruct - 1].rows())) {break;}
        for (size_t idconstr = 0; idconstr < masz; ++idconstr) {
          idAtm = metaatoms[idconstr];
          for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
            gauxREF(idconstr + 1,idcoord) = metaset[istruct](idAtm + 1,idcoord)*dist_Angstrom2aum1;
            gaux0(idconstr + 1,idcoord) = geometry(idAtm + 1,idcoord)*dist_Angstrom2aum1;
          }
        }
        rmsd = RMSD(gaux0,gauxREF,gRMSD);
        Ei = mtdfactors[istruct]*exp(-mtdalphas[istruct]*rmsd*rmsd);
        Ebias += Ei;
        Ei *= -2.0*mtdalphas[istruct]*rmsd;
        for (size_t idconstr = 0; idconstr < masz; ++idconstr) {
          idAtm = metaatoms[idconstr];
          icnt = 3*idAtm;
          for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
            gradients(icnt + idcoord,1) += Ei*gRMSD(idconstr + 1,idcoord);
          }
        }
      }
    }
    return Ebias;
  }
  void NumericalHessian(matrixE & hessian, const matrixE & geometry, double step = 1.0e-3, bool resetHessian = true) {
    //numerical Hessian with central differences
    int Natoms = geometry.rows();
    matrixE gplus(3*Natoms,1);
    matrixE gminus(3*Natoms,1);
    matrixE stepGeom = geometry;               //matrix containing the actual step
    double invstep = 0.5/step;
    double emtd;
    size_t icounter = 1;
    size_t jcounter = 1;
    if (resetHessian) {
      hessian.resize(3*Natoms,3*Natoms);
      hessian.zero();
    }
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      for (size_t idcoord = 0; idcoord < 3; ++idcoord,++icounter) {
        //zero gradients
        gplus.zero();
        gminus.zero();
        //do geometry + step
        stepGeom(idatm + 1,idcoord + 1) = geometry(idatm + 1,idcoord + 1) + step;
        emtd = this->Metadynamic(stepGeom,gplus);
        //do geometry - step
        stepGeom(idatm + 1,idcoord + 1) = geometry(idatm + 1,idcoord + 1) - step;
        emtd = this->Metadynamic(stepGeom,gminus);
        //reset geometry
        stepGeom(idatm + 1,idcoord + 1) = geometry(idatm + 1,idcoord + 1);
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
  void Randomize(matrixE & geom, double scale = 1.0) {
    //function that randomizes the geometry prior to metadynamics calculations in order to get the first structure for bias potential
    size_t Natoms = geom.rows();
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        geom(idAtm + 1,idcoord + 1) += fRandom(-1.0,1.0)*scale;
      }
    }
  }
};
class Dynamics: public MetaDynamics {
  //class for molecular dynamics on a certain macro-structure
  //note that though there is only one geometry, this might include several molecules
  Molecule system;
  matrixE geometry;
  matrixE geometryNew;
  matrixE inertia;
  matrixE velocity;
  matrixE acceleration;
  matrixE velocityAux;
  matrixE velocityN;
  matrixE gradients;
  matrixE block;
  matrixE property;                                       //property to be measured along time
  matrixE propauxmat;                                     //auxiliary matrix for property calculation
  matrixE constraintList;                                 //SHAKE constraints
  matrixE DRorig;
  matrixE DR;
  std::vector<double> masses;
  std::vector<double> Lmoment;
  std::vector<double> propaux;
  std::vector<double> propaux2;
  std::vector<double> auxV;
  double mtdramp;                                         //ramp parameter
  double cmass[3];
  double velcmass[3];
  double velang[3];
  double maxsimulationtime;                               //maximum simulation time for each run (enter in ps, working in a.u.)
  double minequilibrationtime;                            //minimum simulation time for equilibration (enter in ps, working in a.u.)
  double timestep;                                        //timestep (enter in ps, working in a.u.)
  double totalmass;                                       //total mass of the system
  double scaleprotons;                                    //scaling factor for protons so that they don't just dance everywhere
  double timedumpgeom;                                    //frequency at which geometry is dumped (enter in ps, working in a.u.)
  double timedumptraj;                                    //frequency at which trajectories are dumped (enter in ps, working in a.u.)
  //FMCM Add this also to MSMC
  double SCFEaccuracy;                                    //accuracy with which SCF energy is converged
  double SCFDaccuracy;                                    //accuracy with which SCF density is converged
  double temperature;                                     //temperature in Kelvin
  double Telectron;                                       //electronic temperature for fractional occupation numbers in DFT(B) methods
  double Tinit;
  double KBEh;                                            //Boltzmann constant in Hartree
  double Ekin;                                            //kinetic energy
  double Epot;                                            //potential energy
  double Etotal;                                          //total energy
  double Eatm;                                            //energy per atom to check whether something is broken
  double aux;
  std::vector<size_t> atoms;
  matrix<int> neighbours;
  int Natoms;
  int nsteps;                                             //number of steps in integration
  int ndumptraj;                                          //number of steps between trajectory dumps to file
  int ndumpgeom;                                          //number of steps between geometry dumping
  int degreesoffreedom;                                   //how many degrees of freedom allowed to play with
  int SHAKE;                                              //apply SHAKE?
                                                          // 0 -> no
                                                          // 1 -> only on protons
                                                          // 2 -> all atoms
                                                          // 3 -> user defined
  int nmeasurements;                                      //how many measurements were done?
  int measurementfreq;                                    //how frequently should property be measured?
  int blocksz;                                            //size of the block matrix
  int dumpcounter;                                        //counter for dumped geometries
  int nrowprops;
  int szprpaux;
  int mtdmaxstruct;                                       //number of structures considered in metadynamics
  int mtddumpfreq;                                        //frequency with which metadynamics structures are dumped
  bool equilibration;                                     //bool to determine whether to equilibrate
  bool optgeom;                                           //determine whether to optimize geometries at all
  bool dumpgeom;                                          //determines whether to write to file
  bool dumpvelocity;                                      //determines whether to dump velocities as well
  bool restart;                                           //control whether to restart
  bool restartF;                                          //control whether to restart from file
  bool transrotConstr;                                    //translational and rotational constraints
  bool usethermostat;                                     //using thermostat?
  bool measure;                                           //measure properties?
  bool doMetaDyn;                                         //control variable that determines whether to do metadynamics
  bool useContrainPot;                                    //use contrain potential?
  bool forceMD;                                           //this will force the MD run no matter what
  EllipsoidalConstraint ConstrPot;                        //contraining potential
  std::string ensemble;                                   //type of ensemble to be used
                                                          //NVT (canonical) -> default
                                                          //NVE (microcanonical)
  std::string restartFile;                                //file to read restart data from
  std::string restartOutput;                              //file to write at the end of the run to allow for restart from data
  std::string equilibrationFile;                          //file to write trajectories during equilibration into
  std::string trajectoryFile;                             //file to write trajectories into
  std::string geometryFile;                               //file to write geometries into
  std::string thermostat;                                 //type of thermostat to use
  std::string integration;                                //type of integration to use
  std::string propmeasure;                                //property to measure
  std::string auxstrng;
public:
  Dynamics(const Molecule & _system, double _timemax = 50, bool _equilib = true, bool optimizegeom = true, double _timestep = 1.0e-3, double tdumpg = 1.0) {
    setOptions(_timemax,_equilib,optimizegeom,_timestep,tdumpg);
    setNewSystem(_system);
  }
  Dynamics(double _timemax = 50, bool _equilib = true, bool optimizegeom = true, double _timestep = 1.0e-3, double tdumpg = 1.0) {
    setOptions(_timemax,_equilib,optimizegeom,_timestep,tdumpg);
  }
  ~Dynamics() {}
  //getters
  bool Equilibration() {return equilibration;}
  bool GeometryOptBool() {return optgeom;}
  bool GeometryDumpBool() {return dumpgeom;}
  bool VelocityDumpBool() {return dumpvelocity;}
  bool Restart() {return restart;}
  bool RestartFromFile() {return restartF;}
  bool UseThermostat() {return usethermostat;}
  bool Measure() {return measure;}
  bool DoMTD() {return doMetaDyn;}
  bool ConstraintPotential() {return useContrainPot;}
  bool ForceMD() {return forceMD;}
  int DegreesOfFreedom() {return degreesoffreedom;}
  int SHAKEType() {return SHAKE;}
  int NSteps() {return nsteps;}
  int NDumpGeometries() {return ndumpgeom;}
  int NDumpTrajectories() {return ndumptraj;}
  int DumpCounter() {return dumpcounter;}
  int MeasurementFrequency() {return measurementfreq;}
  int MaxNumberMTDStructures() {return mtdmaxstruct;}
  int MTDDumpingFrequency() {return mtddumpfreq;}
  double MaxSimulationTime() {return 0.001*maxsimulationtime*au2fs;}          //returned in ps
  double MinEquilibrationTime() {return 0.001*minequilibrationtime*au2fs;}    //returned in ps
  double TimeStep() {return 0.001*timestep*au2fs;}                            //returned in ps
  double ScaleProtons() {return scaleprotons;}
  double TimeDumpGeometry() {return 0.001*timedumpgeom*au2fs;}                //returned in ps
  double TimeDumpTrajectories() {return 0.001*timedumptraj*au2fs;}            //returned in ps
  double SCFEnergyAccuracy() {return SCFEaccuracy;}
  double SCFDensityAccuracy() {return SCFDaccuracy;}
  double Temperature() {return temperature;}
  double TElectron() {return Telectron;}
  double MTDRamp() {return mtdramp;}
  Molecule System() {return system;}
  std::string EnsembleType() {return ensemble;}
  std::string RestartFile() {return restartFile;}
  std::string RestartOutput() {return restartOutput;}
  std::string TrajectoryFile() {return trajectoryFile;}
  std::string EquilibrationFile() {return equilibrationFile;}
  std::string GeometryFile() {return geometryFile;}
  std::string Thermostat() {return thermostat;}
  std::string Integration() {return integration;}
  std::string PropertyToMeasure() {return propmeasure;}
  matrixE Property() {return property;}
  matrixE ConstraintList() {return constraintList;}
  //setters
  void setOptions(double _tmax, bool _equilib, bool optimizegeom, double _timestep, double tdumpg) {
    maxsimulationtime = 1000.0*_tmax*fs2au;          // <-> tmax
    minequilibrationtime = 5000.0*fs2au;             // <-> mintime
    equilibration = _equilib;
    optgeom = optimizegeom;
    dumpgeom = true;
    dumpvelocity = true;
    timestep = 1000.0*_timestep*fs2au;                      // <-> tstep_md,tstep0
    nsteps = this->calcNstep(maxsimulationtime,timestep);
    timedumpgeom = 1000.0*tdumpg*fs2au;                     // <-> dump_md
    timedumptraj = 1000.0*tdumpg*fs2au;                     // <-> dump_md
    ndumpgeom = this->calcNstep(timedumpgeom,timestep);     // <-> cdump0
    ndumptraj = this->calcNstep(timedumptraj,timestep);     // <-> dumpstep
    scaleprotons = 3.96847889203;
    ensemble = "NVT";                         // <-> thermostat,nvt_md
    usethermostat = true;
    thermostat = "Berendsen";
    SCFEaccuracy = 1.0e-7;                    //well converged
    SCFDaccuracy = 1.0e-4;                    //well-converged
    temperature = 298.15;
    restart = false;
    restartF = false;
    restartFile = "";
    restartOutput = "";
    SHAKE = 0;
    transrotConstr = false;
    useContrainPot = false;
    KBEh = KB/au2J;
    trajectoryFile = "trajectories.trj";
    equilibrationFile = "equilibration.trj";
    geometryFile = "geom";
    dumpcounter = 0;
    measurementfreq = 50;
    Telectron = 300.0;
    inertia.resize(3,3);
    Lmoment.resize(3);
    integration = "Verlet";
    measure = false;
    doMetaDyn = false;
    mtddumpfreq = 2*nsteps;
    //setting seed for random numbers
    srand(time(NULL));
    std::cout << "temporarily using fixed seed\n";
    srand(100);
    forceMD = false;
  }
  void setDumpCounter(int _dumpcounter) {dumpcounter = _dumpcounter;}
  void resetDumpCounter() {dumpcounter = 0;}
  void setMeasurementFrequency(int freq) {measurementfreq = freq;}
  void setEquilibration(bool _equilib) {equilibration = _equilib;}
  void setGeometryOptBool(bool optimizegeom) {optgeom = optimizegeom;}
  void setGeometryDumpBool(bool dumpgeometry) {dumpgeom = dumpgeometry;}
  void setVelocityDumpBool(bool dumpvel) {dumpvelocity = dumpvel;}
  void setRestart(bool rstrt) {restart = rstrt;}
  void setRestartFromFile(bool rstrt) {
    restartF = rstrt;
    restart = rstrt;
  }
  void ApplyContraintPotential(bool constr, double r1 = 1.0, double r2 = -1.0, double r3 = -1.0, int pottype = 0, int beta = 6.0, int alpha = 30.0) {
    //setting up a constraining potential; note that this has nothing to do with the constraints defined in the functions below
    //by default center potential with molecule
    useContrainPot = constr;
    ConstrPot.setRadii(r1,r2,r3);
    std::vector<double> cmass = system.CM();
    ConstrPot.setCenter(cmass[0],cmass[1],cmass[2]);
    ConstrPot.setTemperature(Telectron);
    ConstrPot.setTypeOfPotential(pottype);
    ConstrPot.setExponent(alpha);
    ConstrPot.setBeta(beta);
  }
  void setTransRotConstraint(bool transrot) {ApplyConstraints(transrot,SHAKE);}
  void ApplyConstraints(bool transrot, int SHAKE_ = 0) {
    transrotConstr = transrot;
    SHAKE = SHAKE_;                                  //these are applied later
    degreesoffreedom = 3*Natoms - 6*transrotConstr;
  }
  void setForceMD(bool newfmd) {forceMD = newfmd;}
  void setSHAKE(int _SHAKE) {ApplyConstraints(transrotConstr,_SHAKE);}
  void setMaxSimulationTime(double _timemax) {maxsimulationtime = _timemax*fs2au;}
  void setMinEquilibrationTime(double _timemin) {minequilibrationtime = _timemin*fs2au;}
  void setTimeStep(double _timestep) {timestep = _timestep*fs2au;}
  void setTimeDumpGeometry(double topt) {timedumpgeom = topt*fs2au;}
  void setTimeDumpTrajectory(double tdump) {timedumptraj = tdump*fs2au;}
  void setScaleProtons(double scalefactor) {scaleprotons = scalefactor;}
  void setSCFEnergyAccuracy(double acc) {SCFEaccuracy = acc;}
  void setSCFDensityAccuracy(double acc) {SCFDaccuracy = acc;}
  void setTemperature(double newT) {temperature = newT;}
  void setTElectron(double newT) {Telectron = newT;}
  void setEnsembleType(std::string newensemble) {
    ensemble = newensemble;
    usethermostat = (ensemble == "NVT");
  }
  void setRestartFile(std::string filename) {
    restartFile = filename + ".restart";
    if (filename + "blahblah" != "blahblah") {
      restartF = true;
      restart = true;
    }
    else {
      restartF = false;
      restart = false;
    }
  }
  void setRestartOutput(std::string filename) {restartOutput = filename + ".restart";}
  void setTrajectoryFile(std::string filename) {trajectoryFile = filename + ".trj";}
  void setEquilibrationFile(std::string filename) {equilibrationFile = filename + ".trj";}
  void setGeometryFile(std::string filename) {geometryFile = filename;}
  void setThermostat(std::string newthermo) {
    if ((newthermo == "Berendsen")||(newthermo == "BERENDSEN")||(newthermo == "berendsen")) {thermostat = "Berendsen";}
    else if ((newthermo == "Nose-Hoover")||(newthermo == "NOSE-HOOVER")||(newthermo == "nose-hoover")) {thermostat = "Nose-Hoover";}
    else if ((newthermo == "Maxwell-Boltzmann")||(newthermo == "MAXWELL-BOLTZMANN")||(newthermo == "maxwell-boltzmann")) {thermostat = "Maxwell-Boltzmann";}
    usethermostat = true;
    ensemble = "NVT";
  }
  void RemoveThermostat() {usethermostat = false;}
  void setIntegration(std::string method) {
    if ((method == "Verlet")||(method == "VERLET")||(method == "verlet")) {integration = "Verlet";}
    else if ((method == "LeapFrog")||(method == "LEAPFROG")||(method == "leapfrog")) {integration = "LeapFrog";}
  }
  void setPropertyToMeasure(std::string followprop) {
    if ((followprop == "harmonic-frequencies")||(followprop == "Harmonic-Frequencies")) {
      propmeasure = "harmonic-frequencies";
      nrowprops = 3*Natoms;
    }
    else if ((followprop == "atomic-polarizabilities")||(followprop == "Atomic-Polarizabilities")) {
      propmeasure = "atomic-polarizabilities";
      nrowprops = Natoms;
    }
    else if ((followprop == "partial-charges")||(followprop == "Partial-Charges")) {
      propmeasure = "partial-charges";
      nrowprops = Natoms;
    }
    else if ((followprop == "molecular-polarizabilities")||(followprop == "Molecular-Polarizabilities")) {
      propmeasure = "molecular-polarizabilities";
      propaux.resize(1);
      nrowprops = 1;
    }
    else if ((followprop == "molecular-dispersion")||(followprop == "Molecular-Dispersion")) {
      propmeasure = "molecular-dispersion";
      propaux.resize(2);
      nrowprops = 2;
    }
    else if ((followprop == "heat-of-formation")||(followprop == "Heat-Of-Formation")) {
      propmeasure = "heat-of-formation";
      propaux.resize(1);
      nrowprops = 1;
    }
    else if ((followprop == "ionization-potential")||(followprop == "IP")||(followprop == "ip")) {
      propmeasure = "ip";
      propaux.resize(1);
      nrowprops = 1;
    }
    else if ((followprop == "electron-affinity")||(followprop == "EA")||(followprop == "ea")) {
      propmeasure = "ea";
      propaux.resize(1);
      nrowprops = 1;
    }
    else if ((followprop == "HSAB")||(followprop == "hsab")) {
      propmeasure = "hsab";
      propaux.resize(2);
      nrowprops = 2;
    }
    else if ((followprop == "Fukui-Electrophiles")||(followprop == "Fukui-Indices+")||(followprop == "fukui-indices+")) {
      propmeasure = "fukui+";
      nrowprops = Natoms;
      propaux.resize(Natoms);
    }
    else if ((followprop == "Fukui-Nucleophiles")||(followprop == "Fukui-Indices-")||(followprop == "fukui-indices-")) {
      propmeasure = "fukui-";
      nrowprops = Natoms;
      propaux.resize(Natoms);
    }
    else if ((followprop == "Fukui-Radicals")||(followprop == "Fukui-Indices0")||(followprop == "fukui-indices0")) {
      propmeasure = "fukui0";
      nrowprops = Natoms;
      propaux.resize(Natoms);
    }
    else if ((followprop == "Softness-Electrophiles")||(followprop == "Softness+")||(followprop == "softness+")) {
      propmeasure = "softness+";
      nrowprops = Natoms;
      propaux.resize(Natoms);
    }
    else if ((followprop == "Softness-Nucleophiles")||(followprop == "Softness-")||(followprop == "softness-")) {
      propmeasure = "softness-";
      nrowprops = Natoms;
      propaux.resize(Natoms);
    }
    else if ((followprop == "Softness-Radicals")||(followprop == "Softness0")||(followprop == "softness0")) {
      propmeasure = "softness0";
      nrowprops = Natoms;
      propaux.resize(Natoms);
    }
    measure = true;
  }
  void setPropertyToMeasureOff() {measure = false;}
  void setMTDRamp(double _ramp) {mtdramp = _ramp;}
  void setMaxMTDStructures(int numbstruct) {mtdmaxstruct = numbstruct;}
  void setMTDDumpFrequency(int mtddpfrq) {mtddumpfreq = mtddpfrq;}
  void setMetaDynamics(bool doit, int numbstruct = 10, double k_factor = 0.02, double alpha_exponent = 1.2, int mtddpfrq = 1000, double _ramp = 0.03) {
    //function that sets up the metadynamics calculation
    doMetaDyn = doit;
    mtdmaxstruct = numbstruct;
    kpush = k_factor;
    alpha = alpha_exponent;
    mtdfactors.resize(mtdmaxstruct);
    mtdalphas.resize(mtdmaxstruct);
    for (size_t idx = 0; idx < mtdmaxstruct; ++idx) {
      mtdfactors[idx] = k_factor;
      mtdalphas[idx] = alpha;
    }
    mtdramp = _ramp;
    mtddumpfreq = mtddpfrq;
    metaset.resize(mtdmaxstruct);
  }
  void setNewSystem(const Molecule & newsystem) {
    //function that sets the system's coordinates from an already existing structure
    //it is assumed that the coordinates were already processed by the molecule class
    system = newsystem;
    Natoms = system.Natoms();
    atoms = system.Atoms();
    geometry = system.Geometry();
    system.ConnectivityMatrix(neighbours);    //getting neighbours
    masses.resize(Natoms);
    totalmass = system.AtomicMasses(masses,true,scaleprotons);
    velocity.resize(Natoms,3);
    acceleration.resize(Natoms,3);
    velocityAux.resize(Natoms,3);
    velocityN.resize(Natoms,3);
    degreesoffreedom = 3*Natoms;
  }
  //thermostat functions
  void MaxwellBoltzmann(matrixE & vel, double Tgoal) {
    //Maxwell-Boltzmann thermostat
    //T. A. Andrea, W. C. Swope, H. C. Andersen, J. Chem. Phys., 79, 4576, 1983
    double sigma;
    double velcomp;
    double totalmass = 0.0;
    double vcmass[3];
    vcmass[0] = 0.0;
    vcmass[1] = 0.0;
    vcmass[2] = 0.0;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      sigma = sqrt(KBEh*Tgoal/masses[idAtm]);
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        velcomp = GaussianRandom(0.0,sigma,1.0e-14);
        vel(idAtm + 1,idcoord + 1) = velcomp;
        vcmass[idcoord] += velcomp*masses[idAtm];
        totalmass += masses[idAtm];
      }
    }
    vcmass[0] /= totalmass;
    vcmass[1] /= totalmass;
    vcmass[2] /= totalmass;
    velcomp = 0.0;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        vel(idAtm + 1,idcoord + 1) -= vcmass[idcoord];
        sigma = vel(idAtm + 1,idcoord + 1);
        velcomp += sigma*sigma*masses[idAtm];
      }
    }
    sigma = sqrt(Tgoal*KBEh*double(degreesoffreedom)/velcomp);
    vel *= sigma;
  }
  void Berendsen(matrixE & vel, double step, double tau, double Ti, double Ta) {
    //Berendsen thermostat as declared in 
    //H. J. C. Berendsen, J. P. M. Postma, W. F. van Gunsteren, A. DiNola, J. R. Haak, J. Chem. Phys., 81(8), 3684, 1984
    double Tscal = sqrt(1.0 + (step/tau)*(Ti/Ta - 1.0));
    vel *= Tscal;
  }
  void NoseHoover() {
    //Nose'-Hoover thermostat
  }
  //SHAKE/RATTLE functions; note that RATTLE is a variant of SHAKE where velocities are also corrected
  void AddSHAKEPair(size_t atmA, size_t atmB) {
    //function that adds a particular atom pair to the SHAKE constraint list
    int ccol = constraintList.cols();
    if (ccol < 2) {
      if (constraintList.rows() == 1) {ccol = 0;}
    }
    constraintList.resize(2,ccol + 1);
    constraintList(1,ccol + 1) = atmA;
    constraintList(2,ccol + 1) = atmB;
    SHAKE = 3;
  }
  void RemoveSHAKEPair(size_t atmA, size_t atmB) {
    //function that adds a particular atom pair to the SHAKE constraint list
    size_t ccol = constraintList.cols();
    bool foundpair = false;
    for (size_t idcons = 0; idcons < ccol; ++idcons) {
      if ((atmA == constraintList(1,idcons + 1))&&(atmB == constraintList(2,idcons + 1))) {foundpair = true;}
      else if ((atmA == constraintList(2,idcons + 1))&&(atmB == constraintList(1,idcons + 1))) {foundpair = true;}
      if (foundpair) {
        constraintList(1,idcons + 1) = constraintList(1,ccol);
        constraintList(2,idcons + 1) = constraintList(2,ccol);
        --ccol;
        if (ccol > 0) {constraintList.resize(2,ccol);}
        else {constraintList.clear();}
        break;
      }
    }
    SHAKE = 3;
  }
  void SetSHAKE2Lists(std::vector<size_t> atmsA, std::vector<size_t> atmsB) {
    //function that sets the SHAKE constraint list to predefined lists supplied by the user
    size_t szA = atmsA.size();
    if (szA == atmsB.size()) {
      constraintList.resize(2,szA);
      for (size_t idcons = 0; idcons < szA; ++idcons) {
        constraintList(1,idcons + 1) = atmsA[idcons];
        constraintList(2,idcons + 1) = atmsB[idcons];
      }
      SHAKE = 3;
    }
    else {std::cout << "WARNING: MolecularDynamics.hpp: Dynamics: SetSHAKE2Lists(): user supplied lists are inconsistent\n";}
  }
  void SetSHAKE2List(matrixE & list) {
    //function that sets the SHAKE constraint list to predefined list supplied by the user
    constraintList = list.trans();
    for (size_t idcons = 0; idcons < constraintList.cols(); ++idcons) {
      constraintList(1,idcons + 1) += 1.0;
      constraintList(2,idcons + 1) += 1.0;
    }
    SHAKE = 3;
  }
  template<class T>
  void SHAKEinit(T & ElecStruct, bool simplify = true) {
    //initialization of SHAKE; this structure avoids doing anything for SHAKE = 0
    //simplify applies for SHAKE = 2, where we don't look again over all atoms and simply use the neighbour list already available
    double r2min;
    double rAB2;
    double aux;
    int Bclosest;
    int consrows = 0;
    int idAtm;
    int idBtm;
    size_t nneighbours = neighbours.cols();
    //count the number of protons
    for (idAtm = 0; idAtm < Natoms; ++idAtm) {
      if (atoms[idAtm] == 1) {++consrows;}
    }
    if (SHAKE == 2) {                     //in this case must consider other atoms as well
      Bclosest = consrows;
      //protons -> 1 neighbour; other atoms up to nneighbours
      consrows = Bclosest + (Natoms - Bclosest)*nneighbours;
    }
    else if (SHAKE == 3) {consrows = constraintList.cols();}
    constraintList.resize(2,consrows);
    consrows = 1;
    if (SHAKE == 1) {
      //protons only
      for (idAtm = 0; idAtm < Natoms; ++idAtm) {
        if (atoms[idAtm] != 1) {continue;}
        r2min = 1000.0;
        for (size_t idB = 0; idB < nneighbours; ++idB) {
          idBtm = neighbours(idAtm + 1,idB + 1) - 1;
          if (idBtm < 0) {break;}
          Distance2(idAtm + 1,idBtm + 1,geometry,aux,rAB2);
          if (rAB2 < r2min) {
            r2min = rAB2;
            Bclosest = idBtm;
          }
        }
        if (atoms[Bclosest] == 1) {
          if (Bclosest > idAtm) {
            constraintList(1,consrows) = idAtm + 1;
            constraintList(2,consrows) = Bclosest + 1;
            ++consrows;
          }
        }
        else {
          constraintList(1,consrows) = idAtm + 1;
          constraintList(2,consrows) = Bclosest + 1;
          ++consrows;
        }
      }
    }
    else if (SHAKE == 2) {
      //everyone on board
      bool metalA;
      bool metalB;
      if (!simplify) {                  //this is the way consistent with xTB
        double Rcov;
        double wbothresh;
        //get Mayer/Wiberg bond indices
        matrixE WBO;
        ElecStruct.MayerBondOrder(WBO);
        for (idAtm = 0; idAtm < Natoms; ++idAtm) {
          metalA = MetalHead(atoms[idAtm]);
          if (metalA) {
            if ((atoms[idAtm] != 3)&&(atoms[idAtm] != 4)) {continue;}
          }
          for (idBtm = idAtm + 1; idBtm < Natoms; ++idBtm) {
            metalB = MetalHead(atoms[idBtm]);
            if (metalB) {
              if ((atoms[idBtm] != 3)&&(atoms[idBtm] != 4)) {continue;}
            }
            Distance2(idAtm + 1,idBtm + 1,geometry,aux,rAB2);
            Rcov = AtmRadii(atoms[idAtm]) + AtmRadii(atoms[idBtm]);
            if (rAB2 > 1.44*Rcov*Rcov) {continue;}
            wbothresh = 0.5;
            if ((metalA)||(metalB)) {wbothresh = 0.1;}
            if (WBO(idAtm + 1,idBtm + 1) < wbothresh) {continue;}
            constraintList(1,consrows) = idAtm + 1;
            constraintList(2,consrows) = idBtm + 1;
            ++consrows;
          }
        }
      }
      else {                            //this is a faster way of doing the same
        for (idAtm = 0; idAtm < Natoms; ++idAtm) {
          for (size_t idB = 0; idB < nneighbours; ++idB) {
            idBtm = neighbours(idAtm + 1,idB + 1) - 1;
            if (idBtm < 0) {break;}
            if (idBtm < idAtm) {continue;}
            constraintList(1,consrows) = idAtm + 1;
            constraintList(2,consrows) = idBtm + 1;
            ++consrows;
          }
        }
      }
      constraintList.resize(2,consrows - 1);
    }
    consrows = constraintList.cols();
    DRorig.resize(consrows,4);
    DR.resize(consrows,4);
    for (size_t idpair = 0; idpair < consrows; ++idpair) {
      idAtm = constraintList(1,idpair + 1);
      idBtm = constraintList(2,idpair + 1);
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        DRorig(idpair + 1,idcoord + 1) = geometry(idAtm,idcoord + 1) - geometry(idBtm,idcoord + 1);
      }
      DRorig(idpair + 1,4) = DRorig(idpair + 1,1)*DRorig(idpair + 1,1) + DRorig(idpair + 1,2)*DRorig(idpair + 1,2) + DRorig(idpair + 1,3)*DRorig(idpair + 1,3);
    }
  }
  bool SHAKEit(matrixE & velocity, matrixE & geomA, matrixE & geomB, double SHAKEtol = 1.0e-7, size_t maxiter = 250) {
    //apply SHAKE geometric constraints
    bool converged = false;
    double tau1 = 1.0/timestep;
    double tau2 = tau1*tau1;
    double maxdeviation;
    double dev;
    double delta;
    double mA;
    double mB;
    double lambda;
    size_t nconstraints = constraintList.cols();
    size_t idAtm;
    size_t idBtm;
    for (size_t idpair = 0; idpair < nconstraints; ++idpair) {
      idAtm = constraintList(1,idpair + 1);
      idBtm = constraintList(2,idpair + 1);
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        DRorig(idpair + 1,idcoord + 1) = geomA(idAtm,idcoord + 1) - geomA(idBtm,idcoord + 1);
      }
    }
    geomA = geomB;
    for (size_t ishake = 0; ishake < maxiter; ++ishake) {
      maxdeviation = 0.0;
      for (size_t idpair = 0; idpair < nconstraints; ++idpair) {
        idAtm = constraintList(1,idpair + 1);
        idBtm = constraintList(2,idpair + 1);
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          DR(idpair + 1,idcoord + 1) = geomA(idAtm,idcoord + 1) - geomA(idBtm,idcoord + 1);
        }
        DR(idpair + 1,4) = DR(idpair + 1,1)*DR(idpair + 1,1) + DR(idpair + 1,2)*DR(idpair + 1,2) + DR(idpair + 1,3)*DR(idpair + 1,3);
        delta = DRorig(idpair + 1,4);
        dev = fabs(DR(idpair + 1,4) - delta)/delta;
        if (dev > maxdeviation) {maxdeviation = dev;}
      }
      //is largest deviation below tolerance?
      if (maxdeviation < SHAKEtol) {
        converged = true;
        break;
      }
      for (size_t idpair = 0; idpair < nconstraints; ++idpair) {
        idAtm = constraintList(1,idpair + 1);
        idBtm = constraintList(2,idpair + 1);
        delta = DRorig(idpair + 1,4);
        mA = 1.0/masses[idAtm - 1];
        mB = 1.0/masses[idBtm - 1];
        //get the Lagrange multiplier
        dev = 2.0*(mA + mB)*(DR(idpair + 1,1)*DRorig(idpair + 1,1) + DR(idpair + 1,2)*DRorig(idpair + 1,2) + DR(idpair + 1,3)*DRorig(idpair + 1,3));
        lambda = (delta - DR(idpair + 1,4))/dev;
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          geomA(idAtm,idcoord + 1) += mA*lambda*DRorig(idpair + 1,idcoord + 1);
          geomA(idBtm,idcoord + 1) -= mB*lambda*DRorig(idpair + 1,idcoord + 1);
        }
      }
    }
    if (converged) {
      velocity += (geomA - geomB)*tau1*dist_Angstrom2aum1;
      acceleration += (geomA - geomB)*tau2*dist_Angstrom2aum1;
      geomB = geomA;
    }
    else {std::cout << "WARNING: MolecularDynamics.hpp: Dynamics: SHAKEit(): SHAKE did not converge\n";}
    return converged;
  }
  bool RATTLEit(double tstep, matrixE & velocity, double RATTLEtol = 1.0e-7, size_t maxiter = 250) {
    //apply RATTLE velocity constraints
    //note that positions don't have to be updated, it was already done in SHAKEit (geometric part of RATTLE)
    bool converged = false;
    size_t nconstraints = constraintList.cols();
    size_t idAtm;
    size_t idBtm;
    double maxdeviation;
    double dev;
    double mA;
    double mB;
    double dRdV;
    double mu;
    double delta;
    for (size_t irattle = 0; irattle < maxiter; ++irattle) {
      maxdeviation = 0.0;
      for (size_t idpair = 0; idpair < nconstraints; ++idpair) {
        idAtm = constraintList(1,idpair + 1);
        idBtm = constraintList(2,idpair + 1);
        dRdV = 0.0;
        for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
          dRdV += DR(idpair + 1,idcoord)*(velocity(idAtm,idcoord) - velocity(idBtm,idcoord));
        }
        dev = dRdV*tstep/DR(idpair + 1,4);
        if (dev > maxdeviation) {maxdeviation = dev;}
      }
      //is largest deviation below tolerance?
      if (maxdeviation < RATTLEtol) {
        converged = true;
        break;
      }
      for (size_t idpair = 0; idpair < nconstraints; ++idpair) {
        idAtm = constraintList(1,idpair + 1);
        idBtm = constraintList(2,idpair + 1);
        delta = DRorig(idpair + 1,4);
        mA = 1.0/masses[idAtm - 1];
        mB = 1.0/masses[idBtm - 1];
        dRdV = 0.0;
        for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
          dRdV += DR(idpair + 1,idcoord)*(velocity(idAtm,idcoord) - velocity(idBtm,idcoord));
        }
        //get the Lagrange multiplier
        mu = dRdV/((mA + mB)*DR(idpair + 1,4));
        mA *= mu;
        mB *= mu;
        //update velocities
        for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
          velocity(idAtm,idcoord) -= mA*DR(idpair + 1,idcoord);
          velocity(idBtm,idcoord) += mB*DR(idpair + 1,idcoord);
        }
      }
    }
    return converged;
  }
  //other functions
  int calcNstep(double ttime, double tstep) {
    //function that calculates the number of steps
    return int(ttime/tstep);
  }
  void EstimateTimeStep() {
    //function to estimate a reasonable time step
    double minmass = 100000.0;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      if (minmass > Weight(atoms[idAtm])) {minmass = Weight(atoms[idAtm]);}
    }
    timestep = fs2au*pow(minmass/Weight(1),1.0/3.0);
  }
  void DumpTrajectory(std::ofstream & gfile, std::string comment = "", bool dumpvel = true, int prc = 7) {
    //function that dumps trajectory into file; the file must be open previously
    if (!gfile.is_open()) {throw std::string("ERROR: MolecularDynamics.hpp: Dynamics: DumpTrajectory(): file for writing not open");}
    double precision = pow(10.0,-prc);
    double number;
    std::string AtmID;
    gfile << std::fixed;
    gfile << std::setprecision(prc);
    //dump geometry
    gfile << Natoms << "\n";
    gfile << comment << "\n";
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      AtmID = AtomNr2Symbol(atoms[idAtm]);
      gfile << AtmID << "    ";
      if (AtmID.length() == 1) {gfile << " ";}
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        number = geometry(idAtm + 1,idcoord + 1);
        if (fabs(number) < precision) {number = 0.0;}
        if (number >= -precision) {gfile << " ";}
        gfile << number << "    ";
      }
      gfile << "\n";
    }
    //dump velocity
    if (dumpvel) {
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        gfile << "    ";
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          number = velocity(idAtm + 1,idcoord + 1);
          if (fabs(number) < precision) {number = 0.0;}
          if (number >= -precision) {gfile << " ";}
          gfile << number << "    ";
        }
        gfile << "\n";
      }
    }
    gfile.flush();
  }
  void ReadRestart(matrixE & vel) {
    //function to read geometries from a given xyz file
    std::ifstream gfile(restartFile, std::ios::in);
    if (!gfile.is_open()) {throw std::string("ERROR: MolecularDynamics.hpp: Dynamics: ReadRestart(): restart file could not be open");}
    std::string read;
    gfile >> Natoms;
    geometry.resize(Natoms,3);
    vel.resize(Natoms,3);
    atoms.resize(Natoms);
    //taking care of the comment
    int atmnr;
    for (size_t idx = 0; idx < 500; ++idx) {
      gfile >> read;
      atmnr = Symbol2AtomNr(read);
      if (atmnr != 0) {
        atoms[0] = atmnr;
        gfile >> geometry(1,1);
        gfile >> geometry(1,2);
        gfile >> geometry(1,3);
        break;
      }
    }
    for (size_t idAtm = 1; idAtm < Natoms; ++idAtm) {             //read geometry and atoms
      gfile >> read;
      atoms[idAtm] = Symbol2AtomNr(read);
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        gfile >> geometry(idAtm + 1,idcoord + 1);
      }
    }
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {             //read velocities
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        gfile >> vel(idAtm + 1,idcoord + 1);
      }
    }
    gfile.close();
    system.setGeometry(geometry);
    system.setAtoms(atoms);
  }
  double Ekinetic(const matrixE & velocity, const std::vector<double> & mass) {
    //function that calculates the total kinetic energy of a system
    int Natoms = velocity.rows();
    double auxV[3];
    double energy = 0.0;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      auxV[0] = velocity(idAtm + 1,1);
      auxV[1] = velocity(idAtm + 1,2);
      auxV[2] = velocity(idAtm + 1,3);
      energy += mass[idAtm]*(auxV[0]*auxV[0] + auxV[1]*auxV[1] + auxV[2]*auxV[2]);
    }
    return 0.5*energy;
  }
  void InitVel(double Ekin) {
    //function that initializes velocities
    double Eatom = 2.0*Ekin/double(3*Natoms);         //this is no temperature, but a kinetic energy per atomic coordinate
    double velA;
    double randomNum;
    double factor;
    double factor2;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      factor2 = 1.0 + double(atoms[idAtm] == 1);
      velA = sqrt(Eatom/masses[idAtm]);
      randomNum = iRandom01();
      factor = 1.0;
      if (randomNum > 0.5) {factor = -1.0;}
      velocity(idAtm + 1,1) = velA*factor*factor2;
      randomNum = iRandom01();
      factor = 1.0;
      if (randomNum > 0.5) {factor = -1.0;}
      velocity(idAtm + 1,2) = velA*factor*factor2;
      randomNum = iRandom01();
      factor = 1.0;
      if (randomNum > 0.5) {factor = -1.0;}
      velocity(idAtm + 1,3) = velA*factor*factor2;
    }
  }
  void AverageET(double & Eavg, double & Tavg, int sizebl) {
    //function that averages temperatures and potential energies from the block matrix
    Eavg = 0.0;
    Tavg = 0.0;
    for (size_t idx = 0; idx < sizebl; ++idx) {
      Eavg += block(idx + 1,1);
      Tavg += block(idx + 1,2);
    }
    Eavg /= double(sizebl);
    Tavg /= double(sizebl);
  }
  double Regression(int istart, int iend, const std::vector<double> & regressions) {
    //function that calculates the slope of energy to check whether it is drifting
    double numberelements = double(iend - istart) + 1.0;
    double sx = 0.0;
    double sy = 0.0;
    double sxx = 0.0;
    double sxy = 0.0;
    double xx = 0.0;
    for (size_t idx = istart; idx < iend; ++idx) {
      xx += 1.0;
      sx += xx;
      sy += regressions[idx];
      sxx += xx*xx;
      sxy += xx*regressions[idx];
    }
    return (numberelements*sxy - sx*sy)/(numberelements*sxx - sx*sx);
  }
  void RemoveTransRot(const std::vector<double> & masses, matrixE & velocities) {
    //remove translations and rotations from velocities
    for (size_t idx = 0; idx < 3; ++idx) {
      cmass[idx] = 0.0;
      velcmass[idx] = 0.0;
      Lmoment[idx] = 0.0;
      inertia(idx + 1,idx + 1) = 0.0;
      for (size_t idy = 0; idy < idx; ++idy) {
        inertia(idx + 1,idy + 1) = 0.0;
        inertia(idy + 1,idx + 1) = 0.0;
      }
    }
    //get center of mass and the geometry
    aux = 0.0;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      cmass[0] += geometryNew(idAtm + 1,1)*masses[idAtm];
      cmass[1] += geometryNew(idAtm + 1,2)*masses[idAtm];
      cmass[2] += geometryNew(idAtm + 1,3)*masses[idAtm];
      velcmass[0] += velocities(idAtm + 1,1)*masses[idAtm];
      velcmass[1] += velocities(idAtm + 1,2)*masses[idAtm];
      velcmass[2] += velocities(idAtm + 1,3)*masses[idAtm];
      aux += masses[idAtm];
    }
    cmass[0] /= aux;
    cmass[1] /= aux;
    cmass[2] /= aux;
    //get the angular moment and the inertia matrix
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      //recenter geometry
      geometry(idAtm + 1,1) = dist_Angstrom2aum1*(geometryNew(idAtm + 1,1) - cmass[0]);
      geometry(idAtm + 1,2) = dist_Angstrom2aum1*(geometryNew(idAtm + 1,2) - cmass[1]);
      geometry(idAtm + 1,3) = dist_Angstrom2aum1*(geometryNew(idAtm + 1,3) - cmass[2]);
      Lmoment[0] += masses[idAtm]*(geometry(idAtm + 1,2)*velocities(idAtm + 1,3) - geometry(idAtm + 1,3)*velocities(idAtm + 1,2));
      Lmoment[1] += masses[idAtm]*(geometry(idAtm + 1,3)*velocities(idAtm + 1,1) - geometry(idAtm + 1,1)*velocities(idAtm + 1,3));
      Lmoment[2] += masses[idAtm]*(geometry(idAtm + 1,1)*velocities(idAtm + 1,2) - geometry(idAtm + 1,2)*velocities(idAtm + 1,1));
      inertia(1,1) += masses[idAtm]*(geometry(idAtm + 1,2)*geometry(idAtm + 1,2) + geometry(idAtm + 1,3)*geometry(idAtm + 1,3));
      inertia(2,2) += masses[idAtm]*(geometry(idAtm + 1,3)*geometry(idAtm + 1,3) + geometry(idAtm + 1,1)*geometry(idAtm + 1,1));
      inertia(3,3) += masses[idAtm]*(geometry(idAtm + 1,1)*geometry(idAtm + 1,1) + geometry(idAtm + 1,2)*geometry(idAtm + 1,2));
      inertia(1,2) -= masses[idAtm]*geometry(idAtm + 1,1)*geometry(idAtm + 1,2);
      inertia(1,3) -= masses[idAtm]*geometry(idAtm + 1,1)*geometry(idAtm + 1,3);
      inertia(2,3) -= masses[idAtm]*geometry(idAtm + 1,2)*geometry(idAtm + 1,3);
    }
    inertia(2,1) = inertia(1,2);
    inertia(3,1) = inertia(1,3);
    inertia(3,2) = inertia(2,3);
    //get angular velocity omega
    Solve_Ax_eq_b(inertia,Lmoment);
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      velang[0] = (Lmoment[1]*geometry(idAtm + 1,3) - Lmoment[2]*geometry(idAtm + 1,2));
      velang[1] = (Lmoment[2]*geometry(idAtm + 1,1) - Lmoment[0]*geometry(idAtm + 1,3));
      velang[2] = (Lmoment[0]*geometry(idAtm + 1,2) - Lmoment[1]*geometry(idAtm + 1,1));
      velocities(idAtm + 1,1) -= velcmass[0]/aux + velang[0];
      velocities(idAtm + 1,2) -= velcmass[1]/aux + velang[1];
      velocities(idAtm + 1,3) -= velcmass[2]/aux + velang[2];
    }
  }
  //measurement related functions
  template<class T>
  void Measure(T & ElecStruct, int position) {
    //function measuring a property; always write it to propaux and then process it
    if (property.cols() <= position) {property.resize(nrowprops,position + 20);}
    if (propmeasure == "harmonic-frequencies") {propaux = ElecStruct.CalcVibrFrequencies();}
    else if (propmeasure == "atomic-polarizabilities") {
      if (ElecStruct.Type() != "GFN2") {ElecStruct.getMullikenCharges();}
      propaux2 = ElecStruct.getQAtoms();
      ElecStruct.AtomicPolarizabilities(propaux,propaux2);
    }
    else if (propmeasure == "partial-charges") {
      if (ElecStruct.Type() != "GFN2") {ElecStruct.getMullikenCharges();}
      propaux = ElecStruct.getCharges();
    }
    else if (propmeasure == "molecular-polarizabilities") {
      if (ElecStruct.Type() != "GFN2") {ElecStruct.getMullikenCharges();}
      propaux2 = ElecStruct.getQAtoms();
      ElecStruct.TotalPolarizability(propaux[0],propaux2);
    }
    else if (propmeasure == "molecular-dispersion") {
      if (ElecStruct.Type() != "GFN2") {ElecStruct.getMullikenCharges();}
      propaux2 = ElecStruct.getQAtoms();
      ElecStruct.TotalDispersion(propaux[0],propaux[1],propaux2);
    }
    else if (propmeasure == "heat-of-formation") {
      propaux[0] = ElecStruct.getHeatFormation();
    }
    else if (propmeasure == "ip") {
      propaux[0] = ElecStruct.IonizationPotential(false);
    }
    else if (propmeasure == "ea") {
      propaux[0] = ElecStruct.ElectronAffinity();
    }
    else if (propmeasure == "hsab") {
      ElecStruct.HSABdata(propaux[0],propaux[1]);
    }
    else if (propmeasure == "fukui+") {
      ElecStruct.ReactivityIndices(propauxmat,false);
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        propaux[idAtm] = propauxmat(idAtm + 1,1);
      }
    }
    else if (propmeasure == "fukui-") {
      ElecStruct.ReactivityIndices(propauxmat,false);
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        propaux[idAtm] = propauxmat(idAtm + 1,2);
      }
    }
    else if (propmeasure == "fukui0") {
      ElecStruct.ReactivityIndices(propauxmat,false);
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        propaux[idAtm] = propauxmat(idAtm + 1,3);
      }
    }
    else if (propmeasure == "softness+") {
      ElecStruct.ReactivityIndices(propauxmat,false);
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        propaux[idAtm] = propauxmat(idAtm + 1,4);
      }
    }
    else if (propmeasure == "softness-") {
      ElecStruct.ReactivityIndices(propauxmat,false);
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        propaux[idAtm] = propauxmat(idAtm + 1,5);
      }
    }
    else if (propmeasure == "softness0") {
      ElecStruct.ReactivityIndices(propauxmat,false);
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        propaux[idAtm] = propauxmat(idAtm + 1,6);
      }
    }
    szprpaux = propaux.size();
    for (size_t idcnt = 0; idcnt < szprpaux; ++idcnt) {
      property(idcnt + 1,position) = propaux[idcnt];
    }
  }
  //actual MD functions
  template<class T>
  void runMD(T & ElecStruct, double newtemp = -1.0, double driftthreshold = 2.0e-3, int printfrequency = 200) {
    //function to run the MD simulation
    if (int(newtemp) > 0) {temperature = newtemp;}
    ElecStruct.setElectronTemp(Telectron);
    //read from restart file?
    if ((restart)&&(restartF)) {                    //restart from file
      ReadRestart(velocity);
    }
    //making sure the setup is correct
    ElecStruct.setThresholdDensity(SCFDaccuracy);
    ElecStruct.setThresholdEnergy(SCFEaccuracy);
    ElecStruct.setMolecule(system);
    ElecStruct.setRestart(0);                       //use restart, otherwise this will be a killer
    ElecStruct.Calculate(0,200,true,true,1600.0);
    if (optgeom) {
      //run geometry optimization
      geometry = system.Geometry();
    }
    if (SHAKE > 0) {                                //SHAKE initialization?
      SHAKEinit(ElecStruct);
      degreesoffreedom -= constraintList.cols();
    }
    Tinit = temperature;
    if (equilibration) {
      if (!restart) {                               //initialize velocities
        Ekin = 0.5*(1.0 + double(!usethermostat))*Tinit*KBEh*double(degreesoffreedom);
        InitVel(Ekin);
      }
      //else -> already considered above
      std::cout << "Running equilibration \n";
      if (integration == "Verlet") {Verlet(ElecStruct,100.0*fs2au,true,driftthreshold,printfrequency);}
      else if (integration == "LeapFrog") {LeapFrog(ElecStruct,100.0*fs2au,true,driftthreshold,printfrequency);}
    }
    else {
      if (!restart) {            //initialize velocities
        Ekin = 0.5*(1.0 + double(!usethermostat))*Tinit*KBEh*double(degreesoffreedom);
        InitVel(Ekin);
      }
      //else -> already considered above
    }
    restart = true;
    restartF = false;
    std::cout << "Running productive calculation \n";
    if (integration == "Verlet") {Verlet(ElecStruct,500.0*fs2au,false,driftthreshold,printfrequency);}
    else if (integration == "LeapFrog") {LeapFrog(ElecStruct,500.0*fs2au,false,driftthreshold,printfrequency);}
    if (SHAKE > 0) {this->ApplyConstraints(transrotConstr,SHAKE);}        //reset SHAKE constraints
  }
  template<class T>
  void Verlet(T & ElecStruct, double Taut, bool equilibrate, double driftthreshold = 2.0e-3, int printfrequency = 200) {
    //the actual MD loop function
    //printfrequency is the frequency with which output is printed
    //equilibrate determines whether step is of equilibration
    int gradtype = ElecStruct.AvailableGradients();
    int iblock = 0;
    int nblocks = 0;
    int nregressions = -1;
    int gdump = 1;
    int imeasure = 1;
    int tdump = 1;
    int iprintoutput = 1;
    int icount;
    int mtddump = 0;
    int mtdnstruct = 1;
    auxstrng = trajectoryFile;
    if (equilibrate) {auxstrng = equilibrationFile;}
    std::ofstream trjfile(auxstrng, std::ios::out);
    double Taux = 0.0;
    double aux;
    double blockavgE;
    double blockavgT;
    double slope;
    double Epavg = 0.0;
    double Ekavg = 0.0;
    double Tavg = 0.0;
    double gnorm;
    double Edev;
    double mtdtime = 0.0;
    double Emtd;
    double Ecp;
    size_t szregr = 50;
    std::string strngaux;
    std::vector<double> regressions(szregr,0.0);
    geometry = system.Geometry();
    if (doMetaDyn) {MTDPrepare(geometry);}
    Ekin = Ekinetic(velocity,masses);
    std::cout << "initial MD velocity:\n";
    velocity.Print(12);
    std::cout << "initial MD geometry:\n";
    geometry.Print();
    ElecStruct.Calculate(0,200,true,true,1600.0);
    Epot = ElecStruct.getEnergy(1);
    ElecStruct.gEnergy(gradients,gradtype,0,1.0e-8);
    if (useContrainPot) {
      //constraint potential for the molecule
      Ecp = ConstrPot.Calculate(geometry,gradients);
      Epot += Ecp;
    }
    //get acceleration (F = m.a)
    icount = 1;
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      for (size_t idcoord = 0; idcoord < 3; ++idcoord, ++icount) {
        acceleration(idAtm + 1,idcoord + 1) = -gradients(icount,1)*dist_Angstrom2au/masses[idAtm];
      }
    }
    Etotal = Ekin + Epot;
    nmeasurements = 0;
    std::cout << "nsteps           = " << nsteps << "\n";
    std::cout << "ndumptraj        = " << ndumptraj << "\n";
    std::cout << "ndumpgeom        = " << ndumpgeom << "\n";
    std::cout << "degreesoffreedom = " << degreesoffreedom << "\n";
    std::cout << "SHAKE            = " << SHAKE << "\n";
    std::cout << "nmeasurements    = " << nmeasurements << "\n";
    std::cout << "measurementfreq  = " << measurementfreq << "\n";
    std::cout << "dumpcounter      = " << dumpcounter << "\n";
    blocksz = std::min(5000,int(5000.0*fs2au/timestep));
    block.resize(blocksz,2);
    std::cout << "blocksz          = " << blocksz << "\n";
    //main loop
    for (int idstep = 1; idstep < nsteps + 1; ++idstep) {
      if ((imeasure == measurementfreq)&&(measure)&&(!equilibrate)) {
        //measure whatever
        imeasure = 1;
        ++nmeasurements;
        this->Measure(ElecStruct,nmeasurements);
      }
      if (iblock == blocksz) {                         //average
        if (equilibrate) {                             //increase temperature?
          aux = 1.5*Tinit;
          Tinit = fmin(temperature,aux);
        }
        ++nblocks;
        iblock = 0;
        AverageET(blockavgE,blockavgT,blocksz);
        ++nregressions;
        if (nregressions > szregr) {                   //make sure there is always enough space
          szregr += 50;
          regressions.resize(szregr);
        }
        regressions[nregressions] = blockavgE;
        if (nregressions > 2) {slope = Regression(nregressions - 3,nregressions,regressions);}
        else {slope = 99.0;}
      }
      else {                                           //accumulate
        ++iblock;
        block(iblock,1) = Epot;
        block(iblock,2) = Taux;
      }
      geometry = system.Geometry();
      if ((gdump == ndumpgeom)&&(dumpgeom)) {          //dump geometry?
        system.WriteXYZ(geometryFile,dumpcounter,7);
        gdump = 1;
        ++dumpcounter;
      }
      if (iprintoutput == printfrequency) {            //write details?
        std::cout << idstep << ", " << 0.001*idstep*timestep*au2fs << ", " << Epavg/double(idstep);
        std::cout << ", " << Ekin << ", " << Tavg/double(idstep) << ", " << Taux << ", " << Ekin + Epot;
        std::cout << ", " << Etotal/double(idstep) - Ekin - Epot << "\n";
        iprintoutput = 1;
      }
      if (tdump == ndumptraj) {                        //dump trajectory?
        tdump = 0;
        gnorm = 0.0;
        for (size_t idx = 0; idx < 3*Natoms; ++idx) {
          gnorm += gradients(idx + 1,1)*gradients(idx + 1,1);
        }
        strngaux = "Energy = " + to_string(Epot) + " Eh      |g| = " + to_string(sqrt(gnorm)) + " Eh/A";
        this->DumpTrajectory(trjfile,strngaux,dumpvelocity,7);
      }
      //new positions:         drift                                kick 1/2 step
      geometryNew = geometry + velocity*timestep*dist_Angstrom2au + acceleration*0.5*timestep*timestep*dist_Angstrom2au;
      //Apply SHAKE constraints on positions
      if (SHAKE > 0) {SHAKEit(velocity,geometry,geometryNew,1.0e-7,250);}
      //upd velocities
      velocityN = velocity + acceleration*0.5*timestep;
      //update system and electronic structure
      system.setGeometry(geometryNew);
      ElecStruct.setMolecule(system);
      ElecStruct.Calculate(0,200,true,true,1600.0);
      Epot = ElecStruct.getEnergy(1);
      ElecStruct.gEnergy(gradients,gradtype,0,1.0e-8);
      if (useContrainPot) {
        //constraint potential for the molecule
        Ecp = ConstrPot.Calculate(geometryNew,gradients);
        Epot += Ecp;
      }
      if (doMetaDyn) {                  //do metadynamics
        mtdtime += 1.0;
        for (size_t idx = 0; idx < mtdnstruct; ++idx) {
          mtdfactors[idx] = kpush;
        }
        mtdfactors[mtdnstruct - 1] *= (2.0/(1.0 + exp(-mtdramp*mtdtime)) - 1.0);
        if (mtddump > mtddumpfreq) {
          mtddump = 0;
          if (mtdnstruct < mtdmaxstruct) {          //store next structure
            mtdtime = 0.0;
            ++mtdnstruct;
            metaset[mtdnstruct - 1] = geometry;
          }
          else {                                    //shift all previous structures by one
            for (size_t istruct = 1; istruct < mtdmaxstruct; ++istruct) {
              metaset[istruct - 1] = metaset[istruct];
            }
            metaset[mtdmaxstruct - 1] = geometry;
          }
        }
        Emtd = Metadynamic(geometry,gradients);
        Epot += Emtd;
      }
      Eatm = fabs(Epot/double(Natoms));
      if ((Eatm > 1.0e5)||(Taux > 10000.0)) {
        if (!forceMD) {
          std::cout << "WARNING: MolecularDynamics.hpp: Dynamics: MDcycle(): system's energy is absurd; aborting MD run\n";
          break;
        }
      }
      //get acceleration (F = m.a)
      icount = 1;
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        for (size_t idcoord = 0; idcoord < 3; ++idcoord, ++icount) {
          acceleration(idAtm + 1,idcoord + 1) = -gradients(icount,1)*dist_Angstrom2au/masses[idAtm];
        }
      }
      //kinetic energy
      Ekin = Ekinetic(velocityN,masses);
      //temperature
      Taux = 2.0*Ekin/(KBEh*double(degreesoffreedom));
      if (usethermostat) {                                      //thermostat?
        if (thermostat == "Berendsen") {this->Berendsen(velocityN,timestep,Taut,Tinit,Taux);}
        else if (thermostat == "Nose-Hoover") {}
        else if (thermostat == "Maxwell-Boltzmann") {this->MaxwellBoltzmann(velocityN,Tinit);}
      }
      //upd velocity
      velocity = velocityN + 0.5*acceleration*timestep;           //vmiddle
      //apply RATTLE constraints on velocities
      if (SHAKE > 0) {RATTLEit(timestep,velocity,1.0e-7,250);}
      //leave only internal motion
      RemoveTransRot(masses,velocity);
      //kinetic energy
      Ekin = Ekinetic(velocity,masses);
      //temperature
      Taux = 2.0*Ekin/(KBEh*double(degreesoffreedom));
      ++tdump;
      ++gdump;
      ++imeasure;
      ++iprintoutput;
      Etotal += Epot + Ekin;
      Edev = Etotal/double(idstep) - Epot - Ekin;
      Tavg += Taux;
      Epavg += Epot;
      Ekavg += Ekin;
      if ((equilibrate)&&(0.001*double(idstep)*timestep > minequilibrationtime)) {
        //can we prematurely exit?
        if ((nblocks > 1)&&(fabs(blockavgT - Tinit)/Tinit < 0.02)&&(fabs(slope) < driftthreshold)) {
          std::cout << "early equilibration of the system; prematurely exiting the equilibration\n";
          break;
        }
      }
    }
    geometry = system.Geometry();
    std::cout << "final MD geometry:\n";
    geometry.Print();
    std::cout << "final MD velocities:\n";
    velocity.Print(10);
    std::cout << "average properties " << std::endl;
    std::cout << "Epot               :" << Epavg/double(nsteps) << std::endl;
    std::cout << "Ekin               :" << Ekavg/double(nsteps) << std::endl;
    std::cout << "Etot               :" << (Ekavg + Epavg)/double(nsteps) << std::endl;
    std::cout << "T                  :" << Tavg/double(nsteps) << std::endl;
    if ((fabs(Tavg/double(nsteps) - Tinit) > 0.02*Tinit)&&(nsteps > 500)&&(usethermostat)&&(!equilibrate)) {std::cout << "WARNING: MolecularDynamics.hpp: Dynamics: MDcycle(): problem with thermostat" << std::endl;}
    trjfile.close();
    if (restartOutput + "blahblah" != "blahblah") {           //write last state to file for restart from file
      trjfile.open(restartOutput, std::ios::out);
      DumpTrajectory(trjfile,"",true,10);
      trjfile.close();
    }
  }
  template<class T>
  void LeapFrog(T & ElecStruct, double Taut, bool equilibrate, double driftthreshold = 2.0e-3, int printfrequency = 200) {
    //the actual MD loop function
    //printfrequency is the frequency with which output is printed
    //equilibrate determines whether step is of equilibration
    int gradtype = ElecStruct.AvailableGradients();
    int iblock = 0;
    int nblocks = 0;
    int nregressions = -1;
    int gdump = 1;
    int imeasure = 1;
    int tdump = 1;
    int iprintoutput = 1;
    int mtddump = 0;
    int icount;
    int mtdnstruct = 1;
    auxstrng = trajectoryFile;
    if (equilibrate) {auxstrng = equilibrationFile;}
    std::ofstream trjfile(auxstrng, std::ios::out);
    double Taux = 0.0;
    double aux;
    double blockavgE;
    double blockavgT;
    double slope;
    double Epavg = 0.0;
    double Ekavg = 0.0;
    double Tavg = 0.0;
    double gnorm;
    double Edev;
    double mtdtime = 0.0;
    double Emtd;
    double Ecp;
    size_t szregr = 50;
    std::string strngaux;
    std::vector<double> regressions(szregr,0.0);
    geometry = system.Geometry();
    Ekin = Ekinetic(velocity,masses);
    std::cout << "initial MD velocity:\n";
    velocity.Print(12);
    std::cout << "initial MD geometry:\n";
    geometry.Print();
    ElecStruct.Calculate(0,200,true,true,1600.0);
    Epot = ElecStruct.getEnergy(1);
    Etotal = Ekin + Epot;
    nmeasurements = 0;
    std::cout << "nsteps           = " << nsteps << "\n";
    std::cout << "ndumptraj        = " << ndumptraj << "\n";
    std::cout << "ndumpgeom        = " << ndumpgeom << "\n";
    std::cout << "degreesoffreedom = " << degreesoffreedom << "\n";
    std::cout << "SHAKE            = " << SHAKE << "\n";
    std::cout << "nmeasurements    = " << nmeasurements << "\n";
    std::cout << "measurementfreq  = " << measurementfreq << "\n";
    std::cout << "dumpcounter      = " << dumpcounter << "\n";
    std::cout << "mtddpfrq         = " << mtddumpfreq << "\n";
    blocksz = std::min(5000,int(5000.0*fs2au/timestep));
    block.resize(blocksz,2);
    std::cout << "blocksz          = " << blocksz << "\n";
    if (doMetaDyn) {MTDPrepare(geometry);}
    //main loop
    for (int idstep = 1; idstep < nsteps + 1; ++idstep) {
      ElecStruct.Calculate(0,200,true,true,1600.0);
      Epot = ElecStruct.getEnergy(1);
      geometry = system.Geometry();
      ElecStruct.gEnergy(gradients,gradtype,0,1.0e-8);
      if (useContrainPot) {
        //constraint potential for the molecule
        Ecp = ConstrPot.Calculate(geometry,gradients);
        Epot += Ecp;
      }
      if (doMetaDyn) {                  //do metadynamics
        mtdtime += 1.0;
        for (size_t idx = 0; idx < mtdnstruct; ++idx) {
          mtdfactors[idx] = kpush;
        }
        mtdfactors[mtdnstruct - 1] *= (2.0/(1.0 + exp(-mtdramp*mtdtime)) - 1.0);
        if (mtddump > mtddumpfreq) {
          mtddump = 0;
          if (mtdnstruct < mtdmaxstruct) {          //store next structure
            mtdtime = 0.0;
            ++mtdnstruct;
            metaset[mtdnstruct - 1] = geometry;
          }
          else {                                    //shift all previous structures by one
            for (size_t istruct = 1; istruct < mtdmaxstruct; ++istruct) {
              metaset[istruct - 1] = metaset[istruct];
            }
            metaset[mtdmaxstruct - 1] = geometry;
          }
        }
        
        Emtd = Metadynamic(geometry,gradients);
        Epot += Emtd;
      }
      if ((imeasure == measurementfreq)&&(measure)&&(!equilibrate)) {
        //measure whatever
        imeasure = 1;
        ++nmeasurements;
        this->Measure(ElecStruct,nmeasurements);
      }
      Eatm = fabs(Epot/double(Natoms));
      if ((Eatm > 1.0e5)||(Taux > 10000.0)) {
        if (!forceMD) {
          std::cout << "WARNING: MolecularDynamics.hpp: Dynamics: MDcycle(): system's energy is absurd; aborting MD run\n";
          break;
        }
      }
      if (iblock == blocksz) {                         //average
        if (equilibrate) {                             //increase temperature?
          aux = 1.5*Tinit;
          Tinit = fmin(temperature,aux);
        }
        ++nblocks;
        iblock = 0;
        AverageET(blockavgE,blockavgT,blocksz);
        ++nregressions;
        if (nregressions > szregr) {                   //make sure there is always enough space
          szregr += 50;
          regressions.resize(szregr);
        }
        regressions[nregressions] = blockavgE;
        if (nregressions > 2) {slope = Regression(nregressions - 3,nregressions,regressions);}
        else {slope = 99.0;}
      }
      else {                                           //accumulate
        ++iblock;
        block(iblock,1) = Epot;
        block(iblock,2) = Taux;
      }
      if ((gdump == ndumpgeom)&&(dumpgeom)) {          //dump geometry?
        system.WriteXYZ(geometryFile,dumpcounter,7);
        gdump = 1;
        ++dumpcounter;
      }
      if (iprintoutput == printfrequency) {            //write details?
        std::cout << idstep << ", " << 0.001*idstep*timestep*au2fs << ", " << Epavg/double(idstep);
        std::cout << ", " << Ekin << ", " << Tavg/double(idstep) << ", " << Taux << ", " << Ekin + Epot;
        std::cout << ", " << Etotal/double(idstep) - Ekin - Epot << "\n";
        iprintoutput = 1;
      }
      if (tdump == ndumptraj) {                        //dump trajectory?
        tdump = 0;
        gnorm = 0.0;
        for (size_t idx = 0; idx < 3*Natoms; ++idx) {
          gnorm += gradients(idx + 1,1)*gradients(idx + 1,1);
        }
        strngaux = "Energy = " + to_string(Epot) + " Eh      |g| = " + to_string(sqrt(gnorm)) + " Eh/A";
        this->DumpTrajectory(trjfile,strngaux,dumpvelocity,7);
      }
      //get acceleration (F = m.a)
      icount = 1;
      for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
        for (size_t idcoord = 0; idcoord < 3; ++idcoord, ++icount) {
          acceleration(idAtm + 1,idcoord + 1) = -gradients(icount,1)*dist_Angstrom2au/masses[idAtm];
        }
      }
      //upd velocities
      velocityN = velocity + acceleration*0.5*timestep;
      //kinetic energy
      Ekin = Ekinetic(velocityN,masses);
      //temperature
      Taux = 2.0*Ekin/(KBEh*double(degreesoffreedom));
      velocityAux = velocity + acceleration*timestep;
      if (usethermostat) {                                      //thermostat?
        if (thermostat == "Berendsen") {this->Berendsen(velocityAux,timestep,Taut,Tinit,Taux);}
        else if (thermostat == "Nose-Hoover") {}
        else if (thermostat == "Maxwell-Boltzmann") {this->MaxwellBoltzmann(velocityAux,Tinit);}
      }
      //positions at t + dt
      geometryNew = geometry + velocityAux*timestep*dist_Angstrom2au;
      //apply SHAKE constraints
      if (SHAKE > 0) {SHAKEit(velocityAux,geometry,geometryNew,1.0e-7,250);}
      velocity = velocityAux;
      //leave only internal motion
      RemoveTransRot(masses,velocity);
      //update system and electronic structure
      system.setGeometry(geometryNew);
      ElecStruct.setMolecule(system);
      ++tdump;
      ++gdump;
      ++imeasure;
      ++iprintoutput;
      ++mtddump;
      Etotal += Epot + Ekin;
      Edev = Etotal/double(idstep) - Epot - Ekin;
      Tavg += Taux;
      Epavg += Epot;
      Ekavg += Ekin;
      if ((equilibrate)&&(0.001*double(idstep)*timestep > minequilibrationtime)) {
        //can we exit prematurely?
        if ((nblocks > 1)&&(fabs(blockavgT - Tinit)/Tinit < 0.02)&&(fabs(slope) < driftthreshold)) {
          std::cout << "early equilibration of the system; prematurely exiting the equilibration\n";
          break;
        }
      }
    }
    geometry = system.Geometry();
    std::cout << "final MD geometry:\n";
    geometry.Print();
    std::cout << "final MD velocities:\n";
    velocity.Print(10);
    std::cout << "average properties " << std::endl;
    std::cout << "Epot               :" << Epavg/double(nsteps) << std::endl;
    std::cout << "Ekin               :" << Ekavg/double(nsteps) << std::endl;
    std::cout << "Etot               :" << (Ekavg + Epavg)/double(nsteps) << std::endl;
    std::cout << "T                  :" << Tavg/double(nsteps) << std::endl;
    if ((fabs(Tavg/double(nsteps) - Tinit) > 0.02*Tinit)&&(nsteps > 500)&&(usethermostat)&&(!equilibrate)) {std::cout << "WARNING: MolecularDynamics.hpp: Dynamics: MDcycle(): problem with thermostat" << std::endl;}
    trjfile.close();
    if (restartOutput + "blahblah" != "blahblah") {           //write last state to file for restart from file
      trjfile.open(restartOutput, std::ios::out);
      DumpTrajectory(trjfile,"",true,10);
      trjfile.close();
    }
  }
};
class CompositeMTDQM {
  //class that merges a metadynamics object with a quantum chemical method so that the main solver may be used with composite methods
  //this is not supposed to configure any object, it simply couples two objects for calculation purposes
  MetaDynamics * MTDobj;                //pointer to metadynamics class
  Method * ElectronicStr;              //pointer to quantum chemical method
  double energy;
  matrixE geometry;
  matrixE aux;
  matrixE auxHess;
public:
  CompositeMTDQM(MetaDynamics & MTD, Method & EStrct, const matrixE & geom) {
    MTDobj = & MTD;
    ElectronicStr = & EStrct;
    geometry = geom;
    aux.resize(3*geom.rows(),1);
  }
  ~CompositeMTDQM() {}
  std::string Type() {return "CompositeMTDQM";}
  std::vector<matrixE> MetaSet() {return MTDobj->MetaSet();}
  size_t NAtoms() {return ElectronicStr->NAtoms();}
  Molecule & Component() {return ElectronicStr->Component();}
  double getEnergy(bool complete = 0) {
    energy = ElectronicStr->getEnergy(complete);
    energy += MTDobj->Metadynamic(geometry,aux);
    return energy;
  }
  void gEnergy(matrixE & gx, int type = 0, bool project = 0, double threshold = 1.0e-7) {
    //first plain gradients without projection
    ElectronicStr->gEnergy(aux,type,false,threshold);
    energy = MTDobj->Metadynamic(geometry,aux);
    gx = aux;
    if (project) {ElectronicStr->ProjectGradients(gx);}
  }
  void hEnergy(matrixE & hx, int type = 0, bool project = 0, bool setrestart1 = false, double damp = 0.03, double threshold = 1.0e-7, int maxiter = 100) {
    ElectronicStr->hEnergy(hx,type,0,setrestart1,damp,threshold,maxiter);
    MTDobj->NumericalHessian(auxHess,geometry,1.0e-3,true);
    hx += auxHess;
    if (project) {
      auxHess = geometry;
      ProjectForceConstants(hx,auxHess);
      auxHess = geometry;
      ElectronicStr->WilsonBmatrix(hx,auxHess,damp);              //this modifies the Hessian as well -> diagonal matrix
    }
  }
  std::vector<double> CalcVibrFrequencies(bool restartcalcs = true) {
    //wrapper to get vibrational frequencies without headaches
    this->hEnergy(auxHess,0,0,restartcalcs,0.0,1.0e-7,50);
    std::vector<double> Vfreq = ElectronicStr->Vibrations(auxHess,1,1,true);
    return Vfreq;
  }
  void setGeometry(matrixE & _x) {
    geometry = _x;
    ElectronicStr->setGeometry(_x);
  }
  void updGeometry(matrixE & step, bool project = 0) {
    ElectronicStr->updGeometry(step,project);
    geometry = ElectronicStr->Geometry();
  }
  matrixE Geometry() {return geometry;}
  void Calculate(int _prnt = 0, size_t maxiter = 200, bool _DIIS = true, bool _RCA = true, double d2threshold = 1600.0) {
    ElectronicStr->Calculate(_prnt,maxiter,_DIIS,_RCA,d2threshold);
  }
  std::vector<size_t> Atoms() {return ElectronicStr->Atoms();}
};

#endif //_Molecular_Dynamics_
