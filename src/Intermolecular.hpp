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

#ifndef _Intermolecular_Potential_
#define _Intermolecular_Potential_
#include <vector>
#include <string>
#include <math.h>
#include "Geometry.hpp"
#include "nonbonded/DFTDfunc.hpp"
#include "nonbonded/H4func.hpp"
#include "nonbonded/Xfunc.hpp"
#include "parameters/Hplus.hpp"
#include "other/Utilities.hpp"

//description:
//set of functions that calculate contribution of intermolecular forces to total energy
//applicable to semi-empirical methods and DFT
//references:
//D3 -> S. Grimme, J. Antony, S. Ehrlich, H. Krieg, J. Chem. Phys., 132, 154104, 2010
//D4 -> E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher, C. Bannwarth, S. Grimme, J. Chem. Phys., 150, 154122, 2019
//H4 -> J. Rezac, P. Hobza, J. Chem. Theory Comput., 8, 141, 2012
//X  -> J. Rezac, P. Hobza, Chem. Phys. Lett., 506, 286, 2011

//the D3H4X correction, energies and derivatives
double D3H4X(std::vector<size_t> & atoms, matrixE & geometry, const std::string & method, double HBmin = 1.5, double HBmax = 5.5, double D3threshold = 9000.0) {
  //function that calculates the D3H4X correction to a quantum chemical method
  std::vector<double> drAH(3,0.0);
  std::vector<double> drDH(3,0.0);
  double ex = 0.0;
  double erep = 0.0;
  double eh4 = 0.0;
  double ed6 = 0.0;
  double rAB;
  double rAH;
  double rBH;
  double alphaDHA;
  double drAHdrDH;
  double fdamping = -1.0;
  double rAB2;
  double rBC2;
  double aux;
  double aux2;
  double rdH;
  double raH;
  double fo1;
  double fo2;
  double fcc;
  double erad;
  double eang;
  double ept;
  double echargeA;
  double echargeB;
  double eh2o;
  double eparam;
  int natoms = atoms.size();
  int Ccloseby;
  int Ocloseby;
  int donor;
  int acceptor;
  std::vector<double> NC;
  D3threshold *= dist_Angstrom2au*dist_Angstrom2au;
  NCoord(NC,atoms,geometry,EXPcount,false,1600.0,false);
  for (size_t idatm = 0; idatm < natoms; ++idatm) {
    for (size_t idbtm = 0; idbtm < idatm; ++idbtm) {
      //get distance between atoms
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        drAH[idcoord] = geometry(idatm + 1,idcoord + 1) - geometry(idbtm + 1,idcoord + 1);
      }
      rAB = sqrt(drAH[0]*drAH[0] + drAH[1]*drAH[1] + drAH[2]*drAH[2]);
      //dispersion correction
      if (rAB*rAB < D3threshold) {
        //get C6ij coefficients
        drAHdrDH = C6ABD3(idatm + 1,idbtm + 1,NC[idatm],NC[idbtm],atoms[idatm],atoms[idbtm]);    //C6 coefficient
        ed6 += EdispC6(rAB,drAHdrDH,atoms[idatm],atoms[idbtm],method + "D3H4X",fdamping);
      }
      //hydrogen repulsion correction
      if ((atoms[idatm] == 1)&&(atoms[idbtm] == 1)) {erep += Erep(rAB,method);}
      //X-bond
      else if (((atoms[idatm] == 17)||(atoms[idatm] == 35)||(atoms[idatm] == 53))&&((atoms[idbtm] == 7)||(atoms[idbtm] == 8))) {ex += EX(rAB,atoms[idatm],atoms[idbtm],method);}
      else if (((atoms[idbtm] == 17)||(atoms[idbtm] == 35)||(atoms[idbtm] == 53))&&((atoms[idatm] == 7)||(atoms[idatm] == 8))) {ex += EX(rAB,atoms[idbtm],atoms[idatm],method);}
      //H-bond terms
      else if (((atoms[idatm] == 7)||(atoms[idatm] == 8))&&((atoms[idbtm] == 7)||(atoms[idbtm] == 8))) {
        if ((rAB > HBmin)&&(rAB < HBmax)) {
          for (size_t idctm = 0; idctm < natoms; ++idctm) {        //find H-atoms
            if (atoms[idctm] == 1) {
              //get distances between atoms and the angle
              drAHdrDH = 0.0;
              for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                drAH[idcoord] = geometry(idatm + 1,idcoord + 1) - geometry(idctm + 1,idcoord + 1);
                drDH[idcoord] = geometry(idbtm + 1,idcoord + 1) - geometry(idctm + 1,idcoord + 1);
                drAHdrDH += drAH[idcoord]*drDH[idcoord];
              }
              rAH = sqrt(drAH[0]*drAH[0] + drAH[1]*drAH[1] + drAH[2]*drAH[2]);
              rBH = sqrt(drDH[0]*drDH[0] + drDH[1]*drDH[1] + drDH[2]*drDH[2]);
              //get H-bond angle
              aux = drAHdrDH/(rAH*rBH);
              alphaDHA = pi - acos(aux);
              if (alphaDHA < pi/2.0) {          //if you plot fsw you see that it is supposed to go tops to the value 1.0, otherwise it goes to - \infty
                //determine donor and acceptor
                if (rAH < rBH) {                //A is donor
                  donor = idatm;
                  acceptor = idbtm;
                  rdH = rAH;
                  raH = rBH;
                }
                else {                          //B is donor
                  donor = idbtm;
                  acceptor = idatm;
                  rdH = rBH;
                  raH = rAH;
                }
                //get each contribution to H-bond energy term
                erad = frad(rAB);
                aux = fsw(2.0*alphaDHA/pi);
                eang = 1.0 - aux*aux;
                ept = 1.0;
                if (rdH > 1.15) {ept -= fsw((rdH - 1.15)/(0.5*rdH + 0.5*raH - 1.15));}
                eparam = Cconst(atoms.at(donor),atoms.at(acceptor),method);
                eh2o = 1.0;
                if ((atoms.at(idatm) == 8)&&(atoms.at(idbtm) == 8)) {
                  alphaDHA = 0.0;           //counter for hydrogens in valence
                  drAHdrDH = 0.0;           //counter for other atoms in valence
                  for (size_t iddtm = 0; iddtm < natoms; ++iddtm) {
                    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                      drAH[idcoord] = geometry(donor + 1,idcoord + 1) - geometry(iddtm + 1,idcoord + 1);
                    }
                    aux = sqrt(drAH[0]*drAH[0] + drAH[1]*drAH[1] + drAH[2]*drAH[2]);
                    if (atoms.at(iddtm) == 1) {alphaDHA += ValenceCorr(aux,atoms.at(donor),1);}
                    else {drAHdrDH += ValenceCorr(aux,atoms.at(donor),atoms.at(iddtm));}
                  }
                  if (alphaDHA >= 1.0) {             //if H2O
                    aux = CH2O(method) - 1.0;
                    rAH = 0.0;
                    if (alphaDHA > 1.0 && alphaDHA <= 2.0) {rAH = alphaDHA - 1.0;}
                    if (alphaDHA > 2.0 && alphaDHA < 3.0) {rAH = 3.0 - alphaDHA;}
                    rBH = 1.0 - drAHdrDH;
                    if (rBH < 0.0) {rBH = 0.0;}
                    eh2o = 1.0 + aux*rAH*rBH;
                  }
                }
                echargeA = 1.0;
                echargeB = 1.0;
                if (atoms.at(donor) == 7) {         //look for R3NH+
                  alphaDHA = CS("R3NH+",method) - 1.0;
                  drAHdrDH = 0.0;                    //counter for valence
                  for (size_t iddtm = 0; iddtm < natoms; ++iddtm) {
                    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                      drAH[idcoord] = geometry(donor + 1,idcoord + 1) - geometry(iddtm + 1,idcoord + 1);
                    }
                    aux = sqrt(drAH[0]*drAH[0] + drAH[1]*drAH[1] + drAH[2]*drAH[2]);
                    drAHdrDH += ValenceCorr(aux,atoms.at(donor),atoms.at(iddtm));
                  }
                  if (drAHdrDH > 3.0) {drAHdrDH -= 3.0;}
                  else {drAHdrDH = 0.0;}
                  echargeA += drAHdrDH*alphaDHA;
                }
                if (atoms.at(acceptor) == 8) {         //look for RCOO-
                  alphaDHA = CS("RCOO-",method) - 1.0;
                  Ccloseby = -1;
                  Ocloseby = -1;
                  drAHdrDH = 1.0e10;                    //very large distance 
                  rAH = 0.0;
                  //closest carbon
                  for (size_t iddtm = 0; iddtm < natoms; ++iddtm) {
                    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                      drAH[idcoord] = geometry(acceptor + 1,idcoord + 1) - geometry(iddtm + 1,idcoord + 1);
                    }
                    aux = sqrt(drAH[0]*drAH[0] + drAH[1]*drAH[1] + drAH[2]*drAH[2]);
                    rBH = ValenceCorr(aux,atoms.at(acceptor),atoms.at(iddtm));
                    rAH += rBH;
                    if ((rBH > 0.0)&&(atoms.at(iddtm) == 6)&&(aux < drAHdrDH)) {
                      drAHdrDH = aux;
                      Ccloseby = iddtm;
                    }
                  }
                  fo1 = 1.0 - fabs(1.0 - rAH);
                  if (fo1 < 0.0) {fo1 = 0.0;}
                  if (Ccloseby > -1) {
                    //closest oxygen
                    drAHdrDH = 1.0e10;                  //very large distance 
                    rAH = 0.0;
                    for (size_t iddtm = 0; iddtm < natoms; ++iddtm) {
                      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                        drAH[idcoord] = geometry(Ccloseby + 1,idcoord + 1) - geometry(iddtm + 1,idcoord + 1);
                      }
                      aux = sqrt(drAH[0]*drAH[0] + drAH[1]*drAH[1] + drAH[2]*drAH[2]);
                      rBH = ValenceCorr(aux,atoms.at(Ccloseby),atoms.at(iddtm));
                      rAH += rBH;
                      if ((rBH > 0.0)&&(atoms.at(iddtm) == 8)&&(iddtm != acceptor)&&(aux < drAHdrDH)) {
                        drAHdrDH = aux;
                        Ocloseby = iddtm;
                      }
                    }
                    fcc = 1.0 - fabs(3.0 - rAH);
                    if (fcc < 0.0) {fcc = 0.0;}
                  }
                  if (Ocloseby > -1) {
                    //get valence of last oxygen
                    rAH = 0.0;
                    for (size_t iddtm = 0; iddtm < natoms; ++iddtm) {
                      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                        drAH[idcoord] = geometry(Ocloseby + 1,idcoord + 1) - geometry(iddtm + 1,idcoord + 1);
                      }
                      aux = sqrt(drAH[0]*drAH[0] + drAH[1]*drAH[1] + drAH[2]*drAH[2]);
                      rAH += ValenceCorr(aux,atoms.at(Ocloseby),atoms.at(iddtm));
                    }
                    fo2 = 1.0 - fabs(1.0 - rAH);
                    if (fo2 < 0.0) {fo2 = 0.0;}
                  }
                  echargeB += alphaDHA*fo1*fcc*fo2;
                }
                //add correction term
                eh4 += eparam*erad*eang*ept*eh2o*echargeA*echargeB;
              }
            }
          }
        }
      }
    }
  }
  erep *= sHH(method);
  ed6 *= -s6D3(method + "D3H4X");
  return ex + erep + eh4 + ed6;
}
void gD3H4X(std::vector<double> & gradient, std::vector<size_t> & atoms, matrixE & geometry, const std::string & method, double HBmin = 1.5, double HBmax = 5.5, double D3threshold = 9000.0, double thresholdzero = 1.0e-7) {
  //function that calculates the gradient of the D3H4X correction to a quantum chemical method
  //derivatives in Hartree/Angstrom
  int natoms = atoms.size();
  std::vector<double> drAH(3,0.0);
  std::vector<double> drDH(3,0.0);
  std::vector<double> dalphadrAH(3,0.0);
  std::vector<double> dalphadrDH(3,0.0);
  std::vector<double> dRABdA(3,0.0);
  std::vector<double> dRABdDon(3,0.0);
  std::vector<double> dRABdAcc(3,0.0);
  std::vector<double> dCn_terms(natoms,0.0);
  std::vector<int> carboxy(3,0);
  double gex = 0.0;
  double gerep = 0.0;
  double ged6 = 0.0;
  double rAB;
  double rAH;
  double rBH;
  double alphaDHA;
  double drAHdrDH;
  double fdamping = -1.0;
  double aux;
  double aux2;
  double rdH;
  double raH;
  double fo1;
  double fo2;
  double fcc;
  double erad;
  double gerad;
  double eang;
  double geang;
  double ept;
  double geptA;
  double geptD;
  double echargeA;
  double echargeB;
  double gechargeA;
  double gechargeB;
  double eh2o;
  double geh2o;
  double eparam;
  double donoracceptorfactor;
  double s6param = s6D3(method + "D3H4X");
  double sHHparam = sHH(method);
  double angletermderv;
  double VA;
  double gVA;
  double VH;
  double gVH;
  double signH2O;
  double cvo1;
  double cvo2;
  double cvcc;
  double dc6drab_A;
  double dc6drab_B;
  int donor;
  int acceptor;
  gradient.clear();
  gradient.resize(3*natoms);
  std::vector<double> NC;
  D3threshold *= dist_Angstrom2au*dist_Angstrom2au;
  NCoord(NC,atoms,geometry,EXPcount,false,1600.0,false);
  for (size_t idatm = 0; idatm < natoms; ++idatm) {
    for (size_t idbtm = 0; idbtm < idatm; ++idbtm) {
      //get distance between atoms
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        drAH[idcoord] = geometry(idatm + 1,idcoord + 1) - geometry(idbtm + 1,idcoord + 1);
      }
      rAB = sqrt(drAH[0]*drAH[0] + drAH[1]*drAH[1] + drAH[2]*drAH[2]);
      //getting chain rule terms
      dRABdA[0] = drAH[0]/rAB;
      dRABdA[1] = drAH[1]/rAB;
      dRABdA[2] = drAH[2]/rAB;
      ged6 = 0.0;
      gerep = 0.0;
      gex = 0.0;
      //dispersion correction
      if (rAB*rAB < D3threshold) {
        drAHdrDH = C6ABD3(idatm + 1,idbtm + 1,NC[idatm],NC[idbtm],atoms[idatm],atoms[idbtm]);      //C6 coefficient
        ged6 = s6param*gEdispC6(rAB,drAHdrDH,atoms[idatm],atoms[idbtm],method + "D3H4X",fdamping);
        dC6dRAB(dc6drab_A,dc6drab_B,atoms[idatm],atoms[idbtm],rAB,NC[idatm],NC[idbtm],EXPcount,false,1600.0,false);
        dCn_terms[idatm] += s6param*EdispC6(rAB,dc6drab_A,atoms[idatm],atoms[idbtm],method + "D3H4X",fdamping);
        dCn_terms[idbtm] += s6param*EdispC6(rAB,dc6drab_B,atoms[idatm],atoms[idbtm],method + "D3H4X",fdamping);
      }
      //hydrogen repulsion correction
      if ((atoms[idatm] == 1)&&(atoms[idbtm] == 1)) {gerep = sHHparam*gErep(rAB,method)*dist_Angstrom2au;}
      //X-bond
      else if (((atoms[idatm] == 17)||(atoms[idatm] == 35)||(atoms[idatm] == 53))&&((atoms[idbtm] == 7)||(atoms[idbtm] == 8))) {gex = gEX(rAB,atoms[idatm],atoms[idbtm],method);}
      else if (((atoms[idbtm] == 17)||(atoms[idbtm] == 35)||(atoms[idbtm] == 53))&&((atoms[idatm] == 7)||(atoms[idatm] == 8))) {gex = gEX(rAB,atoms[idbtm],atoms[idatm],method);}
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        gradient[3*idatm + idcoord] += (ged6 + gerep + gex)*dRABdA[idcoord]*dist_Angstrom2aum1;
        gradient[3*idbtm + idcoord] -= (ged6 + gerep + gex)*dRABdA[idcoord]*dist_Angstrom2aum1;
      }
      //H-bond terms
      if (((atoms[idatm] == 7)||(atoms[idatm] == 8))&&((atoms[idbtm] == 7)||(atoms[idbtm] == 8))) {
        if ((rAB > HBmin)&&(rAB < HBmax)) {
          for (size_t idctm = 0; idctm < natoms; ++idctm) {        //find H-atoms
            if (atoms[idctm] == 1) {
              //get distances between atoms and the angle
              drAHdrDH = 0.0;
              for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                drAH[idcoord] = geometry(idatm + 1,idcoord + 1) - geometry(idctm + 1,idcoord + 1);
                drDH[idcoord] = geometry(idbtm + 1,idcoord + 1) - geometry(idctm + 1,idcoord + 1);
                drAHdrDH += drAH[idcoord]*drDH[idcoord];
              }
              rAH = sqrt(drAH[0]*drAH[0] + drAH[1]*drAH[1] + drAH[2]*drAH[2]);
              rBH = sqrt(drDH[0]*drDH[0] + drDH[1]*drDH[1] + drDH[2]*drDH[2]);
              //get H-bond angle
              aux = drAHdrDH/(rAH*rBH);
              alphaDHA = pi - acos(aux);
              if (alphaDHA < pi/2.0) {          //if you plot fsw you see that it is supposed to go tops to the value 1.0, otherwise it goes to - \infty
                angletermderv = 1.0/sqrt(1.0 - aux*aux);
                //determine donor and acceptor
                if (rAH < rBH) {                //A is donor
                  donor = idatm;
                  acceptor = idbtm;
                  rdH = rAH;
                  raH = rBH;
                  rAH = 1.0/rdH;
                  rBH = 1.0/raH;
                  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                    dRABdDon[idcoord] = drAH[idcoord]*rAH;
                    dRABdAcc[idcoord] = drDH[idcoord]*rBH;
                  }
                }
                else {                          //B is donor
                  donor = idbtm;
                  acceptor = idatm;
                  rdH = rBH;
                  raH = rAH;
                  rAH = 1.0/rdH;
                  rBH = 1.0/raH;
                  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                    dRABdAcc[idcoord] = drAH[idcoord]*rBH;
                    dRABdDon[idcoord] = drDH[idcoord]*rAH;
                  }
                }
                //getting chain rule terms
                for (int idr = 0; idr < 3; ++idr) {
                  dalphadrDH[idr] = 0.0;
                  dalphadrAH[idr] = 0.0;
                  for (int idc = 0; idc < 3; ++idc) {
                    dalphadrDH[idr] -= angletermderv*dRABdDon[idc]*(dRABdDon[idc]*dRABdAcc[idr] - dRABdDon[idr]*dRABdAcc[idc])*rAH;
                    dalphadrAH[idr] -= angletermderv*dRABdAcc[idc]*(dRABdAcc[idc]*dRABdDon[idr] - dRABdAcc[idr]*dRABdDon[idc])*rBH;
                  }
                }
                //get each contribution to H-bond energy term
                erad = frad(rAB);
                gerad = gfrad(rAB);
                aux = fsw(2.0*alphaDHA/pi);
                eang = 1.0 - aux*aux;
                geang = 4.0*aux*gfsw(2.0*alphaDHA/pi)/pi;
                ept = 1.0;
                geptA = 0.0;
                geptD = 0.0;
                if (rdH > 1.15) {
                  ept -= fsw((rdH - 1.15)/(0.5*rdH + 0.5*raH - 1.15));
                  geptA = 0.5*gfsw((rdH - 1.15)/(0.5*rdH + 0.5*raH - 1.15))*(rdH - 1.15)/((0.5*rdH + 0.5*raH - 1.15)*(0.5*rdH + 0.5*raH - 1.15));    //dfPT/dRAH (acceptor)
                  geptD = -0.5*gfsw((rdH - 1.15)/(0.5*rdH + 0.5*raH - 1.15))*(raH - 1.15)/((0.5*rdH + 0.5*raH - 1.15)*(0.5*rdH + 0.5*raH - 1.15));   //dfPT/dRDH (donor)
                }
                eparam = Cconst(atoms.at(donor),atoms.at(acceptor),method);
                eh2o = 1.0;
                if ((atoms.at(idatm) == 8)&&(atoms.at(idbtm) == 8)) {
                  alphaDHA = 0.0;           //counter for hydrogens in valence
                  drAHdrDH = 0.0;           //counter for other atoms in valence
                  for (size_t iddtm = 0; iddtm < natoms; ++iddtm) {
                    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                      drAH[idcoord] = geometry(donor + 1,idcoord + 1) - geometry(iddtm + 1,idcoord + 1);
                    }
                    aux = sqrt(drAH[0]*drAH[0] + drAH[1]*drAH[1] + drAH[2]*drAH[2]);
                    if (atoms.at(iddtm) == 1) {alphaDHA += ValenceCorr(aux,atoms.at(donor),1);}
                    else {drAHdrDH += ValenceCorr(aux,atoms.at(donor),atoms.at(iddtm));}
                  }
                  if (alphaDHA >= 1.0) {             //if H2O
                    aux = CH2O(method) - 1.0;
                    rAH = 0.0;
                    signH2O = 1.0;
                    if (alphaDHA > 1.0 && alphaDHA <= 2.0) {rAH = alphaDHA - 1.0;}
                    if (alphaDHA > 2.0 && alphaDHA < 3.0) {
                      rAH = 3.0 - alphaDHA;
                      signH2O = -1.0;
                    }
                    rBH = 1.0 - drAHdrDH;
                    if (rBH < 0.0) {rBH = 0.0;}
                    eh2o = 1.0 + aux*rAH*rBH;
                  }
                }
                echargeA = 1.0;
                if (atoms.at(donor) == 7) {         //look for R3NH+
                  alphaDHA = CS("R3NH+",method) - 1.0;
                  drAHdrDH = 0.0;                    //counter for valence
                  for (size_t iddtm = 0; iddtm < natoms; ++iddtm) {
                    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                      drAH[idcoord] = geometry(donor + 1,idcoord + 1) - geometry(iddtm + 1,idcoord + 1);
                    }
                    aux = sqrt(drAH[0]*drAH[0] + drAH[1]*drAH[1] + drAH[2]*drAH[2]);
                    drAHdrDH += ValenceCorr(aux,atoms.at(donor),atoms.at(iddtm));
                  }
                  if (drAHdrDH > 3.0) {drAHdrDH -= 3.0;}
                  else {drAHdrDH = 0.0;}
                  echargeA += drAHdrDH*alphaDHA;
                }
                echargeB = 1.0;
                if (atoms.at(acceptor) == 8) {         //look for RCOO-
                  alphaDHA = CS("RCOO-",method) - 1.0;
                  carboxy[0] = -1;       //C close by
                  carboxy[1] = -1;       //O close by
                  carboxy[2] = acceptor;
                  drAHdrDH = 1.0e10;                    //very large distance 
                  cvo1 = 0.0;
                  //closest carbon
                  for (size_t iddtm = 0; iddtm < natoms; ++iddtm) {
                    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                      drAH[idcoord] = geometry(acceptor + 1,idcoord + 1) - geometry(iddtm + 1,idcoord + 1);
                    }
                    aux = sqrt(drAH[0]*drAH[0] + drAH[1]*drAH[1] + drAH[2]*drAH[2]);
                    rBH = ValenceCorr(aux,atoms.at(acceptor),atoms.at(iddtm));
                    cvo1 += rBH;
                    if ((rBH > 0.0)&&(atoms.at(iddtm) == 6)&&(aux < drAHdrDH)) {
                      drAHdrDH = aux;
                      carboxy[0] = iddtm;
                    }
                  }
                  fo1 = 1.0 - fabs(1.0 - cvo1);
                  if (fo1 < 0.0) {fo1 = 0.0;}
                  if (carboxy[0] > -1) {
                    //closest oxygen
                    drAHdrDH = 1.0e10;                  //very large distance 
                    cvcc = 0.0;
                    for (size_t iddtm = 0; iddtm < natoms; ++iddtm) {
                      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                        drAH[idcoord] = geometry(carboxy[0] + 1,idcoord + 1) - geometry(iddtm + 1,idcoord + 1);
                      }
                      aux = sqrt(drAH[0]*drAH[0] + drAH[1]*drAH[1] + drAH[2]*drAH[2]);
                      rBH = ValenceCorr(aux,atoms.at(carboxy[0]),atoms.at(iddtm));
                      cvcc += rBH;
                      if ((rBH > 0.0)&&(atoms.at(iddtm) == 8)&&(iddtm != acceptor)&&(aux < drAHdrDH)) {
                        drAHdrDH = aux;
                        carboxy[1] = iddtm;
                      }
                    }
                    fcc = 1.0 - fabs(3.0 - cvcc);
                    if (fcc < 0.0) {fcc = 0.0;}
                  }
                  if (carboxy[1] > -1) {
                    //get valence of last oxygen
                    cvo2 = 0.0;
                    for (size_t iddtm = 0; iddtm < natoms; ++iddtm) {
                      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                        drAH[idcoord] = geometry(carboxy[1] + 1,idcoord + 1) - geometry(iddtm + 1,idcoord + 1);
                      }
                      aux = sqrt(drAH[0]*drAH[0] + drAH[1]*drAH[1] + drAH[2]*drAH[2]);
                      cvo2 += ValenceCorr(aux,atoms.at(carboxy[1]),atoms.at(iddtm));
                    }
                    fo2 = 1.0 - fabs(1.0 - cvo2);
                    if (fo2 < 0.0) {fo2 = 0.0;}
                  }
                  echargeB += alphaDHA*fo1*fcc*fo2;
                }
                gerad *= eparam*1.0*eang*ept*eh2o*echargeA*echargeB;       //contribution from radial derivative
                geang *= eparam*erad*1.0*ept*eh2o*echargeA*echargeB;       //contribution from angular derivatives
                geptA *= eparam*erad*eang*1.0*eh2o*echargeA*echargeB;      //contribution from bond switch
                geptD *= eparam*erad*eang*1.0*eh2o*echargeA*echargeB;      //contribution from bond switch
                for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                  gradient[3*idatm + idcoord] += (gerad*dRABdA[idcoord])*dist_Angstrom2aum1;
                  gradient[3*idbtm + idcoord] -= (gerad*dRABdA[idcoord])*dist_Angstrom2aum1;
                  gradient[3*donor + idcoord] += (geang*dalphadrDH[idcoord] + geptD*dRABdDon[idcoord])*dist_Angstrom2aum1;
                  gradient[3*acceptor + idcoord] += (geang*dalphadrAH[idcoord] + geptA*dRABdAcc[idcoord])*dist_Angstrom2aum1;
                  gradient[3*idctm + idcoord] -= (geang*(dalphadrDH[idcoord] + dalphadrAH[idcoord]) + geptD*dRABdDon[idcoord] + geptA*dRABdAcc[idcoord])*dist_Angstrom2aum1;
                }
                //gradient contribution from water term
                if (fabs(eh2o - 1.0) > thresholdzero) {
                  aux2 = CH2O(method) - 1.0;
                  for (size_t iddtm = 0; iddtm < natoms; ++iddtm) {
                    if (iddtm == donor) {continue;}
                    //get distance
                    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                      drAH[idcoord] = geometry(donor + 1,idcoord + 1) - geometry(iddtm + 1,idcoord + 1);
                    }
                    alphaDHA = sqrt(drAH[0]*drAH[0] + drAH[1]*drAH[1] + drAH[2]*drAH[2]);
                    gechargeB = gValenceCorr(alphaDHA,atoms.at(donor),atoms.at(iddtm));
                    if (fabs(gechargeB) > thresholdzero) {
                      geh2o = eparam*erad*eang*ept*aux2*gechargeB*echargeA*echargeB;
                      if (atoms[iddtm] == 1) {geh2o *= -signH2O;}
                      drAH[0] /= alphaDHA;
                      drAH[1] /= alphaDHA;
                      drAH[2] /= alphaDHA;
                      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                        gradient[3*donor + idcoord] -= geh2o*drAH[idcoord]*dist_Angstrom2aum1;
                        gradient[3*iddtm + idcoord] += geh2o*drAH[idcoord]*dist_Angstrom2aum1;
                      }
                    }
                  }
                }
                //gradient contribution from ammonia terms
                if (fabs(echargeA - 1.0) > thresholdzero) {
                  aux2 = CS("R3NH+",method) - 1.0;
                  for (size_t iddtm = 0; iddtm < natoms; ++iddtm) {
                    if (iddtm == donor) {continue;}
                    //get distance
                    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                      drAH[idcoord] = geometry(donor + 1,idcoord + 1) - geometry(iddtm + 1,idcoord + 1);
                    }
                    drAHdrDH = sqrt(drAH[0]*drAH[0] + drAH[1]*drAH[1] + drAH[2]*drAH[2]);
                    gechargeB = gValenceCorr(drAHdrDH,atoms.at(donor),atoms.at(iddtm));
                    if (fabs(gechargeB) > thresholdzero) {
                      gechargeA = -eparam*erad*eang*ept*eh2o*aux2*gechargeB*echargeB;
                      drAH[0] /= drAHdrDH;
                      drAH[1] /= drAHdrDH;
                      drAH[2] /= drAHdrDH;
                      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                        gradient[3*donor + idcoord] -= gechargeA*drAH[idcoord]*dist_Angstrom2aum1;
                        gradient[3*iddtm + idcoord] += gechargeA*drAH[idcoord]*dist_Angstrom2aum1;
                      }
                    }
                  }
                }
                //gradient contribution from carboxy terms
                if (fabs(echargeB - 1.0) > thresholdzero) {
                  aux2 = CS("RCOO-",method) - 1.0;
                  //acceptor atom contributions
                  for (size_t iddtm = 0; iddtm < natoms; ++iddtm) {
                    if ((carboxy[0] < 0)||(carboxy[1] < 0)) {break;}         //in this case it does not really exist a carboxyl
                    if (iddtm != carboxy[2]) {
                      //get distance
                      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                        drAH[idcoord] = geometry(carboxy[2] + 1,idcoord + 1) - geometry(iddtm + 1,idcoord + 1);
                      }
                      drAHdrDH = sqrt(drAH[0]*drAH[0] + drAH[1]*drAH[1] + drAH[2]*drAH[2]);
                      gechargeA = gValenceCorr(drAHdrDH,atoms.at(carboxy[2]),atoms.at(iddtm));
                      if (fabs(gechargeA) > thresholdzero) {
                        gechargeB = -eparam*erad*eang*ept*eh2o*echargeA*aux2*gechargeA*fo2*fcc;
                        if (cvo1 > 1.0) {gechargeB *= -1.0;}
                        drAH[0] /= drAHdrDH;
                        drAH[1] /= drAHdrDH;
                        drAH[2] /= drAHdrDH;
                        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                          gradient[3*carboxy[2] + idcoord] -= gechargeB*drAH[idcoord]*dist_Angstrom2aum1;
                          gradient[3*iddtm + idcoord] += gechargeB*drAH[idcoord]*dist_Angstrom2aum1;
                        }
                      }
                    }
                    if (iddtm != carboxy[1]) {
                      //get distance
                      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                        drAH[idcoord] = geometry(carboxy[1] + 1,idcoord + 1) - geometry(iddtm + 1,idcoord + 1);
                      }
                      drAHdrDH = sqrt(drAH[0]*drAH[0] + drAH[1]*drAH[1] + drAH[2]*drAH[2]);
                      gechargeA = gValenceCorr(drAHdrDH,atoms.at(carboxy[1]),atoms.at(iddtm));
                      if (fabs(gechargeA) > thresholdzero) {
                        gechargeB = -eparam*erad*eang*ept*eh2o*echargeA*aux2*gechargeA*fo1*fcc;
                        if (cvo2 > 1.0) {gechargeB *= -1.0;}
                        drAH[0] /= drAHdrDH;
                        drAH[1] /= drAHdrDH;
                        drAH[2] /= drAHdrDH;
                        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                          gradient[3*carboxy[1] + idcoord] -= gechargeB*drAH[idcoord]*dist_Angstrom2aum1;
                          gradient[3*iddtm + idcoord] += gechargeB*drAH[idcoord]*dist_Angstrom2aum1;
                        }
                      }
                    }
                    if (iddtm != carboxy[0]) {
                      //get distance
                      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                        drAH[idcoord] = geometry(carboxy[0] + 1,idcoord + 1) - geometry(iddtm + 1,idcoord + 1);
                      }
                      drAHdrDH = sqrt(drAH[0]*drAH[0] + drAH[1]*drAH[1] + drAH[2]*drAH[2]);
                      gechargeA = gValenceCorr(drAHdrDH,atoms.at(carboxy[0]),atoms.at(iddtm));
                      if (fabs(gechargeA) > thresholdzero) {
                        gechargeB = -eparam*erad*eang*ept*eh2o*echargeA*aux2*gechargeA*fo1*fo2;
                        if (cvcc > 1.0) {gechargeB *= -1.0;}
                        drAH[0] /= drAHdrDH;
                        drAH[1] /= drAHdrDH;
                        drAH[2] /= drAHdrDH;
                        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                          gradient[3*carboxy[0] + idcoord] -= gechargeB*drAH[idcoord]*dist_Angstrom2aum1;
                          gradient[3*iddtm + idcoord] += gechargeB*drAH[idcoord]*dist_Angstrom2aum1;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  //now the missing dCn/dCN terms
  for (size_t idAtm = 0; idAtm < natoms; ++idAtm) {
    for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
      //get distance between atoms
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        drAH[idcoord] = geometry(idAtm + 1,idcoord + 1) - geometry(idBtm + 1,idcoord + 1);
      }
      rAB = sqrt(drAH[0]*drAH[0] + drAH[1]*drAH[1] + drAH[2]*drAH[2]);
      //getting chain rule terms
      gex = 1.0/rAB;
      dRABdA[0] = drAH[0]*gex;
      dRABdA[1] = drAH[1]*gex;
      dRABdA[2] = drAH[2]*gex;
      aux = dNCoorddR(atoms[idAtm],atoms[idBtm],rAB,dEXPcountdx,false,1600.0,false);   //dCN/dRAB
      ged6 = aux*(dCn_terms[idAtm] + dCn_terms[idBtm]);
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        gradient[3*idAtm + idcoord] += -ged6*dRABdA[idcoord]*dist_Angstrom2aum1;
        gradient[3*idBtm + idcoord] -= -ged6*dRABdA[idcoord]*dist_Angstrom2aum1;
      }
    }
  }
}
//the D3H+ correction, energies and derivatives, only applicable to PM6
double BondDistance(size_t atmnrA, size_t atmnrB) {
  //function that estimates the bond distance between two atoms based on covalent radii
  return covRadius(atmnrA) + covRadius(atmnrB);
}
bool HBondAngles(size_t posA, size_t posB, size_t posC, size_t posD, size_t posH, std::vector<size_t> & atoms, matrixE & geometry, int nbonds, double & cosphi, double & dcosphi, double & cospsi, double & dcospsi, bool & torsioncheck) {
  //function that calculates cos(Phi) and cos(Psi) given a geometry arrangement
  //for a group of the form R1
  //                          \
  //                       R2--X - - - H
  //                          /
  //                        R3
  //then we have posA -> X
  //             posB -> R3
  //             posC -> R2
  //             posD -> R1
  //             posH -> H
  bool check = false;                   //check whether H-bond is physical
  bool psicheck1 = false;
  torsioncheck = false;
  double psicheck = 0.0;
  double psicheckSign = 0.0;
  double psicheckFactor;
  double phishift;
  double phishift2;
  double psishift;
  double psi;
  double phi;
  double psival1;
  double psival2;
  double cosphi2;
  double cospsi2;
  psi = 1.0/180.0;
  if (atoms[posA] == 8) {
    if (nbonds == 1) {                   //this is the case of a carbonyl/carboxyl group: H - - O=C-R1
      phishift = pi;
      phishift2 = 120.0*psi*pi;
      psishift = 0.0;
      torsioncheck = true;
    }
    else {                               //alcohol or something alike
      phishift = 109.48*psi*pi;
      phishift2 = phishift;
      psishift = 54.74*psi*pi;
    }
  }
  else if (atoms[posA] == 7) {
    if (nbonds == 2) {
      phishift = 120.0*psi*pi;
      phishift2 = phishift;
      psishift = 0.0;
    }
    else {
      bool planar = false;
      //check for double bonds
      phi = 0.97*0.75;
      if ((posA != posB)&&(Distance(posA + 1,posB + 1,geometry) < BondDistance(atoms[posA],atoms[posB])*phi)) {planar = true;}
      if ((posA != posC)&&(Distance(posA + 1,posC + 1,geometry) < BondDistance(atoms[posA],atoms[posC])*phi)) {planar = true;}
      if ((posA != posD)&&(Distance(posA + 1,posD + 1,geometry) < BondDistance(atoms[posA],atoms[posD])*phi)) {planar = true;}
      if (!planar) {                  //case NR3 tetragonal
        phishift = 109.48*psi*pi;
        phishift2 = phishift;
        psishift = 54.74*psi*pi;
        psicheck1 = true;
      }
      else {                          //case NR3 planar
        phishift = 120.0*psi*pi;
        phishift2 = phishift;
        psishift = 90.0*psi*pi;
      }
    }
  }
  //extrapolation between tetragonal and planar NR3 group
  if (psicheck1) {
    psicheck = Torsion(posC + 1,posB + 1,posA + 1,posD + 1,geometry);
    if (psicheck < 0.0) {
      psicheck += pi;
      psicheck *= -1.0;
    }
    else {
      psicheck -= pi;
      psicheck *= -1.0;
    }
    psicheckSign = psicheck;
    if (psicheck < 0.0) {psicheck *= -1.0;}
    if (psicheck != 0.0) {psicheckFactor = psicheckSign/psicheck;}
    else {psicheckFactor = 0.0;}
    psicheck *= 180.0/pi;
    phishift2 = phishift;
    psicheck = psicheckSign;      //restore
  }
  //calculate phi
  phi = Angle(posB + 1,posA + 1,posH + 1,geometry);
  cosphi = cos(phishift - phi);
  cosphi2 = cos(phishift2 - phi);
  //now stuff for derivatives
  dcosphi = 2.0*cosphi*sin(phishift - phi);
  if (cosphi2 > cosphi) {
    cosphi = cosphi2;
    dcosphi = 2.0*cosphi*sin(phishift2 - phi);
  }
  if (cosphi <= 0.0) {check = true;}             //unphysical H-bond
  if (!check) {
    psicheckFactor = Torsion(posC + 1,posB + 1,posA + 1,posH + 1,geometry);
    if (torsioncheck) {psicheckFactor = Torsion(posC + 1,posD + 1,posA + 1,posH + 1,geometry);}
    if ((!torsioncheck)||(fabs(psicheckFactor*180.0/pi) > 90.0)) {
      if (psicheckFactor < 0.0) {
        psicheckFactor += pi;
        psicheckFactor *= -1.0;
      }
      else {
        psicheckFactor -= pi;
        psicheckFactor *= -1.0;
      }
    }
    //correction of NR3 torsion angle for through-bond case
    if (psicheck < 0.0) {                            //negative torsion angle occupied by -NR3
      psival1 = psishift - psicheckFactor;
      cospsi = cos(psival1);
      dcospsi = 2.0*cospsi*sin(psival1);
    }
    else if (psicheck > 0.0) {                       //positive torsion angle occupied by -NR3
      psival1 = -psishift - psicheckFactor;
      cospsi = cos(psival1);
      dcospsi = 2.0*cospsi*sin(psival1);
    }
    else {                                           //planar -NR3 or general case
      psival1 = psishift - psicheckFactor;
      psival2 = -psishift - psicheckFactor;
      cospsi = cos(psival1);
      cospsi2 = cos(psival2);
      dcospsi = 2.0*cospsi*sin(psival1);
      if (cospsi2 > cospsi) {
        cospsi = cospsi2;
        dcospsi = 2.0*cospsi*sin(psival2);
      }
    }
    //correction for H - - C=0 instead of C=0 - - H
    if ((Distance(posH + 1,posA + 1,geometry) > Distance(posH + 1,posB + 1,geometry))&&(torsioncheck)) {
      cospsi = 0.0;
      dcospsi = 0.0;
    }
    if ((posB == posA)||(posB == posH)) {
      cosphi = 1.0;
      dcosphi = 0.0;
    }
    if ((posC == posA)||(posC == posH)) {
      cospsi = 1.0;
      dcospsi = 0.0;
    }
    if (posD == posH) {
      cosphi = 1.0;
      cospsi = 1.0;
      dcosphi = 0.0;
      dcospsi = 0.0;
    }
    if (cospsi <= 0.0) {check = true;}             //unphysical H-bond
  }
  return check;
}
bool SettingHbondConfg(std::vector<int> & Hbondconfiguration, int nbonds, size_t posAtm, size_t posH, size_t shift, int Natoms, std::vector<size_t> & atoms, matrixE & geometry, const std::vector<int> & LNeighbour, double rAC, double rAD) {
  //function that sets the Hbond configuration vector
  bool hbs = false;
  int countaux;
  double distN2H[3];                             //array to store the distance of A's neighbours to proton
  double aux;
  if ((nbonds == 3)||(nbonds == 4)) {
    countaux = 0;
    for (size_t idx = 0; idx < nbonds; ++idx) {
      if ((LNeighbour[idx] == posH)&&(nbonds == 4)) {continue;}           //if there are 4 neighbours, then we want to exlude THE proton
      ++countaux;
    }
    if ((nbonds == 4)&&(countaux == 4)) {
      //in this case this is not a good Hbond, since we have an atom (N or O) with 4 neighbours that is not the H-bond donor
      //this violates the octet
      hbs = true;
    }
    else {
      countaux = 0;
      for (size_t idx = 0; idx < nbonds; ++idx) {
        if ((LNeighbour[idx] == posH)&&(nbonds == 4)) {continue;}           //if there are 4 neighbours, then we want to exlude THE proton
        distN2H[countaux] = Distance(LNeighbour[idx] + 1,posH + 1,geometry);
        Hbondconfiguration[countaux + shift] = LNeighbour[idx];
        ++countaux;
      }
      if (distN2H[0] < distN2H[1]) {                                        //re-sort
        swap(Hbondconfiguration[shift + 0],Hbondconfiguration[shift + 1]);
        swap(distN2H[0],distN2H[1]);
      }
      if (distN2H[1] < distN2H[2]) {                                        //re-sort
        swap(Hbondconfiguration[shift + 1],Hbondconfiguration[shift + 2]);
        swap(distN2H[1],distN2H[2]);
        if (distN2H[0] < distN2H[1]) {
          swap(Hbondconfiguration[shift + 0],Hbondconfiguration[shift + 1]);
          swap(distN2H[0],distN2H[1]);
        }
      }
    }
  }
  else if (nbonds == 2) {
    for (size_t idx = 0; idx < 2; ++idx) {
      distN2H[idx] = Distance(LNeighbour[idx] + 1,posH + 1,geometry);
      Hbondconfiguration[idx + shift] = LNeighbour[idx];
    }
    if (distN2H[0] < distN2H[1]) {                                        //re-sort
      swap(Hbondconfiguration[shift + 0],Hbondconfiguration[shift + 1]);
      swap(distN2H[0],distN2H[1]);
    }
    //determine whether donor or acceptor
    if (rAC < BondDistance(atoms[posAtm],atoms[posH])) {Hbondconfiguration[shift + 2] = posH;}
    else {Hbondconfiguration[shift + 2] = Hbondconfiguration[shift - 1];}
  }
  else if (nbonds == 1) {
    Hbondconfiguration[shift + 0] = LNeighbour[0];
    //look for neighbours of atom in LNeighbour[0] to get neighbours
    countaux = 0;
    int ListNeighbouraux[6];                 //this stores the list of neighbours for an auxiliary atom, which might have extended octet
    for (size_t iddtm = 0; iddtm < Natoms; ++iddtm) {
      if (iddtm == LNeighbour[0]) {continue;}
      rAD = Distance(LNeighbour[0] + 1,iddtm + 1,geometry);
      if (rAD < BondDistance(atoms[LNeighbour[0]],atoms[iddtm])) {
        distN2H[countaux] = rAD;
        ListNeighbouraux[countaux] = iddtm;
        ++countaux;
      }
      if (countaux > 5) {break;}
    }
    //find the neighbour most distant to proton
    aux = -1.0;
    for (size_t idx = 0; idx < countaux; ++idx) {
      if (distN2H[idx] > aux) {
        Hbondconfiguration[shift + 1] = ListNeighbouraux[idx];
        aux = distN2H[idx];
      }
    }
    //determine whether donor or acceptor
    if (rAC < BondDistance(atoms[posAtm],atoms[posH])) {Hbondconfiguration[shift + 2] = posH;}
    else {Hbondconfiguration[shift + 2] = Hbondconfiguration[shift - 1];}
  }
  else if (nbonds == 0) {
    //determine whether donor or acceptor
    if (rAC < BondDistance(atoms[posAtm],atoms[posH])) {
      Hbondconfiguration[shift + 0] = posH;
      Hbondconfiguration[shift + 1] = posH;
      Hbondconfiguration[shift + 2] = posH;
    }
    else {
      Hbondconfiguration[shift + 0] = Hbondconfiguration[shift - 1];
      Hbondconfiguration[shift + 1] = Hbondconfiguration[shift - 1];
      Hbondconfiguration[shift + 2] = Hbondconfiguration[shift - 1];
    }
  }
  else {hbs = true;}
  return hbs;
}
double D3Hplus(std::vector<size_t> & atoms, matrixE & geometry, const std::string & method, bool doC8 = true, double D3threshold = 9000.0, double rDAthreshold = 10.5, double rHDthreshold = 1.4) {
  //function that calculates the D3Hplus correction to a quantum chemical method, the D3 correction of Grimme with the h-bond term of Korth
  //J. C. Kromann, A. S. Christensen, C. Steinmann, M. Korth, J. H. Jensen, PeerJ., 2, e449, 2014
  //doC8 controls whether the C8 terms are also to be included
  //rDAthreshold is the maximum distance between donor and acceptor
  //rHDthreshold is the zero covalent cut between proton and the donor; this means that the proton must be covalently bound to the donor
  //nomenclature: D -> donor
  //              A -> acceptor
  double ed6 = 0.0;
  double ed8 = 0.0;
  double ehplus = 0.0;
  if (method == "PM6") {       //this was only parametrized for PM6
    int Natoms = atoms.size();
    D3threshold *= dist_Angstrom2au*dist_Angstrom2au;
    double rAB;
    double rAC;
    double rBC;
    double rAD;
    double rBD;
    double rDonorH;
    double rAC2;
    double rBC2;
    double rAD2;
    double rBD2;
    double R2R4A;
    double R2R4B;
    double C6abparam;
    double C8abparam;
    double CN = HplusCParam(7);
    double CO = HplusCParam(8);
    double scaleA;
    double scaleB;
    double scaleC;
    double costheta;
    double cosphiA;
    double cospsiA;
    double cosphiB;
    double cospsiB;
    double dcosphi;
    double dcospsi;
    double fbond;
    double dampshort;
    double damplong;
    int countA;
    int countB;
    int countaux;
    bool badHbond;
    bool hbsA;
    bool hbsB;
    bool torsioncheck;
    std::vector<double> NC;
    std::vector<int> ListNeighbourA(4);                   //this stores the list of neighbours for atom A; oxygen and nitrogen must in worst case follow octet rule, hence size 4
    std::vector<int> ListNeighbourB(4);                   //this stores the list of neighbours for atom B; oxygen and nitrogen must in worst case follow octet rule, hence size 4
    std::vector<int> Hbondconfiguration(9);               //this stores all the intervenient atoms for each H-bond
    NCoord(NC,atoms,geometry,EXPcount,false,1600.0,false);
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      R2R4A = R2R4(atoms[idatm]);
      for (size_t idbtm = idatm + 1; idbtm < Natoms; ++idbtm) {
        R2R4B = R2R4(atoms[idbtm]);
        rAB = Distance(idatm + 1,idbtm + 1,geometry);
        //dispersion correction
        if (rAB*rAB < D3threshold) {
          //get C6ij coefficients
          C6abparam = C6ABD3(idatm + 1,idbtm + 1,NC[idatm],NC[idbtm],atoms[idatm],atoms[idbtm]);   //C6 coefficient
          ed6 += EdispC6(rAB,C6abparam,atoms[idatm],atoms[idbtm],method + "D3H+",-1.0);
          if (doC8) {                                                                              //C8 coefficient
            C8abparam = C8AB(C6abparam,R2R4A,R2R4B);
            ed8 += EdispC8(rAB,C8abparam,atoms[idatm],atoms[idbtm],method + "D3H+",-1.0);
          }
        }
        //H-bond correction
        if (((atoms[idatm] == 7)||(atoms[idatm] == 8))&&((atoms[idbtm] == 7)||(atoms[idbtm] == 8))) {
          if (rAB < rDAthreshold) {
            //if atoms are correct and distance is not too large, then look for protons
            for (size_t idctm = 0; idctm < Natoms; ++idctm) {
              if (atoms[idctm] == 1) {
                rAC2 = 0.0;
                rBC2 = 0.0;
                for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                  scaleA = geometry(idatm + 1,idcoord + 1) - geometry(idctm + 1,idcoord + 1);
                  rAC2 += scaleA*scaleA;
                  scaleA = geometry(idbtm + 1,idcoord + 1) - geometry(idctm + 1,idcoord + 1);
                  rBC2 += scaleA*scaleA;
                }
                if ((rAC2 > rHDthreshold*rHDthreshold)&&(rBC2 > rHDthreshold*rHDthreshold)) {continue;}
                //if we are here we are going to check for H-bond correction
                rAC = sqrt(rAC2);
                rBC = sqrt(rBC2);
                //look for neighbours, but first clean the list and counters
                for (size_t idx = 0; idx < 4; ++idx) {
                  ListNeighbourA[idx] = -1;
                  ListNeighbourB[idx] = -1;
                }
                countA = 0;
                countB = 0;
                for (size_t iddtm = 0; iddtm < Natoms; ++iddtm) {
                  rAD2 = 0.0;
                  rBD2 = 0.0;
                  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                    scaleA = geometry(idatm + 1,idcoord + 1) - geometry(iddtm + 1,idcoord + 1);
                    rAD2 += scaleA*scaleA;
                    scaleA = geometry(idbtm + 1,idcoord + 1) - geometry(iddtm + 1,idcoord + 1);
                    rBD2 += scaleA*scaleA;
                  }
                  rAD = sqrt(rAD2);
                  rBD = sqrt(rBD2);
                  if ((rAD < BondDistance(atoms[idatm],atoms[iddtm]))&&(idatm != iddtm)) {
                    ListNeighbourA[countA] = iddtm;
                    ++countA;
                  }
                  if ((rBD < BondDistance(atoms[idbtm],atoms[iddtm]))&&(idbtm != iddtm)) {
                    ListNeighbourB[countB] = iddtm;
                    ++countB;
                  }
                  if ((countA > 3)&&(countB > 3)) {break;}          //there is no exceeding octets, this means bad geometry
                }
                //check whether H-bond is bad (1-3 and 1-4 interactions)
                badHbond = false;
                for (size_t idcntA = 0; idcntA < 4; ++idcntA) {
                  if (ListNeighbourA[idcntA] == idbtm) {badHbond = true;}     //1-3
                  if ((badHbond)||(ListNeighbourA[idcntA] == -1)) {break;}
                  for (size_t idcntB = 0; idcntB < 4; ++idcntB) {
                    if (ListNeighbourB[idcntB] == -1) {break;}
                    if (ListNeighbourA[idcntA] == ListNeighbourB[idcntB]) {
                      badHbond = true;
                      break;
                    }
                  }
                }
                //proceed only for good H-bonds
                if (badHbond) {continue;}
                //note that countA contains number of neighbours of A and countB for B
                for (size_t idx = 0; idx < 9; ++idx) {
                  Hbondconfiguration[idx] = -1;
                }
                Hbondconfiguration[0] = idatm;
                Hbondconfiguration[4] = idbtm;
                Hbondconfiguration[8] = idctm;
                hbsA = SettingHbondConfg(Hbondconfiguration,countA,idatm,idctm,1,Natoms,atoms,geometry,ListNeighbourA,rAC,rAD);
                if (hbsA) {continue;}          //unknown constellation
                hbsB = SettingHbondConfg(Hbondconfiguration,countB,idbtm,idctm,5,Natoms,atoms,geometry,ListNeighbourB,rBC,rBD);
                if (hbsB) {continue;}          //unknown constellation
                costheta = -cos(Angle(Hbondconfiguration[0] + 1,Hbondconfiguration[8] + 1,Hbondconfiguration[4] + 1,geometry));
                if (costheta < 0.0) {continue;}          //in this case skip because H-bond is not physical
                hbsA = HBondAngles(Hbondconfiguration[0],Hbondconfiguration[1],Hbondconfiguration[2],Hbondconfiguration[3],Hbondconfiguration[8],atoms,geometry,countA,cosphiA,dcosphi,cospsiA,dcospsi,torsioncheck);
                if (hbsA) {continue;}                    //if unphysical, then move on
                hbsB = HBondAngles(Hbondconfiguration[4],Hbondconfiguration[5],Hbondconfiguration[6],Hbondconfiguration[7],Hbondconfiguration[8],atoms,geometry,countB,cosphiB,dcosphi,cospsiB,dcospsi,torsioncheck);
                if (hbsB) {continue;}                    //if unphysical, then move on
                scaleA = CN*(atoms[idatm] == 7) + CO*(atoms[idatm] == 8);
                scaleB = CN*(atoms[idbtm] == 7) + CO*(atoms[idbtm] == 8);
                scaleC = 0.5*(scaleA + scaleB);
                rDonorH = fmin(rAC,rBC);
                //finally in position to calculate energies
                fbond = 1.0 - 1.0/(1.0 + exp(-60.0*(rDonorH/1.2 - 1.0)));
                dampshort = 1.0/(1.0 + exp(-100.0*(rAB/2.4 - 1.0)));
                damplong = 1.0 - 1.0/(1.0 + exp(-10.0*(rAB/7.0 - 1.0)));
                rAC = dist_Angstrom2au/rAB;
                scaleA = scaleC*costheta*costheta*cosphiA*cosphiA*cospsiA*cospsiA*cosphiB*cosphiB*cospsiB*cospsiB*fbond*dampshort*damplong*rAC*rAC;
                ehplus += scaleA;
              }
            }
          }
        }
      }
    }
    ed6 *= -s6D3(method + "D3H+");
    ed8 *= -s8D3(method + "D3H+");
  }
  return ed6 + ed8 + ehplus;
}
void gD3Hplus(std::vector<double> & gradient, std::vector<size_t> & atoms, matrixE & geometry, const std::string & method, bool doC8 = true, double D3threshold = 9000.0, double rDAthreshold = 10.5, double rHDthreshold = 1.4) {
  //function that calculates the gradient of the D3Hplus correction to a quantum chemical method, the D3 correction of Grimme with the h-bond term of Korth
  //J. C. Kromann, A. S. Christensen, C. Steinmann, M. Korth, J. H. Jensen, PeerJ., 2, e449, 2014
  //doC8 controls whether the C8 terms are also to be included
  //rDAthreshold is the maximum distance between donor and acceptor
  //rHDthreshold is the zero covalent cut between proton and the donor; this means that the proton must be covalently bound to the donor
  //nomenclature: D -> donor
  //              A -> acceptor
  //derivatives in Hartree/Angstrom
  int Natoms = atoms.size();
  gradient.clear();
  gradient.resize(3*Natoms);
  if (method == "PM6") {       //this was only parametrized for PM6
    std::vector<double> dCn_terms(Natoms,0.0);
    std::vector<double> NC;
    std::vector<int> ListNeighbourA(4);                   //this stores the list of neighbours for atom A; oxygen and nitrogen must in worst case follow octet rule, hence size 4
    std::vector<int> ListNeighbourB(4);                   //this stores the list of neighbours for atom B; oxygen and nitrogen must in worst case follow octet rule, hence size 4
    std::vector<int> Hbondconfiguration(9);               //this stores all the intervenient atoms for each H-bond
    double ged6 = 0.0;
    double ged8 = 0.0;
    double dehdr;
    double dehdtheta;
    double dehdpsiA;
    double dehdpsiB;
    double dehdphiA;
    double dehdphiB;
    double dedh_aux;
    double deds_aux;
    double dedl_aux;
    double dedh;
    double deds;
    double dedl;
    double rAB;
    double C6abparam;
    double C8abparam;
    double fdamping = -1.0;
    double s6param = s6D3(method + "D3H+");
    double s8param = s8D3(method + "D3H+");
    double dc6drab_A;
    double dc6drab_B;
    double R2R4A;
    double R2R4B;
    double CN = HplusCParam(7);
    double CO = HplusCParam(8);
    double scaleA;
    double scaleB;
    double scaleC;
    double rAC2;
    double rBC2;
    double rAD2;
    double rBD2;
    double rAC;
    double rBC;
    double rAD;
    double rBD;
    double theta;
    double costheta;
    double dcostheta;
    double cosphiA;
    double cospsiA;
    double dcosphiA;
    double dcospsiA;
    double cosphiB;
    double cospsiB;
    double dcosphiB;
    double dcospsiB;
    double rDonorH;
    double fbond;
    double dampshort;
    double damplong;
    double damp;
    int countA;
    int countB;
    bool AisDonor;
    bool badHbond;
    bool hbsA;
    bool hbsB;
    bool torsionCheckA;
    bool torsionCheckB;
    D3threshold *= dist_Angstrom2au*dist_Angstrom2au;
    NCoord(NC,atoms,geometry,EXPcount,false,1600.0,false);
    for (size_t idatm = 0; idatm < Natoms; ++idatm) {
      R2R4A = R2R4(atoms[idatm]);
      for (size_t idbtm = idatm + 1; idbtm < Natoms; ++idbtm) {
        R2R4B = R2R4(atoms[idbtm]);
        //get distance between atoms
        rAB = Distance(idatm + 1,idbtm + 1,geometry);
        ged6 = 0.0;
        ged8 = 0.0;
        //dispersion correction
        if (rAB*rAB < D3threshold) {
          C6abparam = C6ABD3(idatm + 1,idbtm + 1,NC[idatm],NC[idbtm],atoms[idatm],atoms[idbtm]);     //C6 coefficient
          ged6 = s6param*gEdispC6(rAB,C6abparam,atoms[idatm],atoms[idbtm],method + "D3H+",fdamping);
          dC6dRAB(dc6drab_A,dc6drab_B,atoms[idatm],atoms[idbtm],rAB,NC[idatm],NC[idbtm],EXPcount,false,1600.0,false);
          dCn_terms[idatm] += s6param*EdispC6(rAB,dc6drab_A,atoms[idatm],atoms[idbtm],method + "D3H+",fdamping);
          dCn_terms[idbtm] += s6param*EdispC6(rAB,dc6drab_B,atoms[idatm],atoms[idbtm],method + "D3H+",fdamping);
          if (doC8) {                                                                                //C8 coefficient
            C8abparam = C8AB(C6abparam,R2R4A,R2R4B);
            ged8 = s8param*gEdispC8(rAB,C8abparam,atoms[idatm],atoms[idbtm],method + "D3H+",fdamping);
            C8abparam = C8AB(dc6drab_A,R2R4A,R2R4B);
            dCn_terms[idatm] += s8param*EdispC8(rAB,C8abparam,atoms[idatm],atoms[idbtm],method + "D3H+",fdamping);
            C8abparam = C8AB(dc6drab_B,R2R4A,R2R4B);
            dCn_terms[idbtm] += s8param*EdispC8(rAB,C8abparam,atoms[idatm],atoms[idbtm],method + "D3H+",fdamping);
          }
        }
        Internal2Cartesian(2,ged6 + ged8,idatm + 1,idbtm + 1,0,0,geometry,gradient);
        //now the contribution from the Hbond correction
        if (((atoms[idatm] == 7)||(atoms[idatm] == 8))&&((atoms[idbtm] == 7)||(atoms[idbtm] == 8))) {
          if (rAB < rDAthreshold) {
            //if atoms are correct and distance is not too large, then look for protons
            for (size_t idctm = 0; idctm < Natoms; ++idctm) {
              if (atoms[idctm] == 1) {
                rAC2 = 0.0;
                rBC2 = 0.0;
                for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                  scaleA = geometry(idatm + 1,idcoord + 1) - geometry(idctm + 1,idcoord + 1);
                  rAC2 += scaleA*scaleA;
                  scaleA = geometry(idbtm + 1,idcoord + 1) - geometry(idctm + 1,idcoord + 1);
                  rBC2 += scaleA*scaleA;
                }
                if ((rAC2 > rHDthreshold*rHDthreshold)&&(rBC2 > rHDthreshold*rHDthreshold)) {continue;}
                //if we are here we are going to check for H-bond correction
                rAC = sqrt(rAC2);
                rBC = sqrt(rBC2);
                //look for neighbours, but first clean the list and counters
                for (size_t idx = 0; idx < 4; ++idx) {
                  ListNeighbourA[idx] = -1;
                  ListNeighbourB[idx] = -1;
                }
                countA = 0;
                countB = 0;
                for (size_t iddtm = 0; iddtm < Natoms; ++iddtm) {
                  rAD2 = 0.0;
                  rBD2 = 0.0;
                  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                    scaleA = geometry(idatm + 1,idcoord + 1) - geometry(iddtm + 1,idcoord + 1);
                    rAD2 += scaleA*scaleA;
                    scaleA = geometry(idbtm + 1,idcoord + 1) - geometry(iddtm + 1,idcoord + 1);
                    rBD2 += scaleA*scaleA;
                  }
                  rAD = sqrt(rAD2);
                  rBD = sqrt(rBD2);
                  if ((rAD < BondDistance(atoms[idatm],atoms[iddtm]))&&(idatm != iddtm)) {
                    ListNeighbourA[countA] = iddtm;
                    ++countA;
                  }
                  if ((rBD < BondDistance(atoms[idbtm],atoms[iddtm]))&&(idbtm != iddtm)) {
                    ListNeighbourB[countB] = iddtm;
                    ++countB;
                  }
                  if ((countA > 3)&&(countB > 3)) {break;}          //there is no exceeding octets, this means bad geometry
                }
                //check whether H-bond is bad (1-3 and 1-4 interactions)
                badHbond = false;
                for (size_t idcntA = 0; idcntA < 4; ++idcntA) {
                  if (ListNeighbourA[idcntA] == idbtm) {badHbond = true;}     //1-3
                  if ((badHbond)||(ListNeighbourA[idcntA] == -1)) {break;}
                  for (size_t idcntB = 0; idcntB < 4; ++idcntB) {
                    if (ListNeighbourB[idcntB] == -1) {break;}
                    if (ListNeighbourA[idcntA] == ListNeighbourB[idcntB]) {
                      badHbond = true;
                      break;
                    }
                  }
                }
                //proceed only for good H-bonds
                if (badHbond) {continue;}
                //note that countA contains number of neighbours of A and countB for B
                for (size_t idx = 0; idx < 9; ++idx) {
                  Hbondconfiguration[idx] = -1;
                }
                Hbondconfiguration[0] = idatm;
                Hbondconfiguration[4] = idbtm;
                Hbondconfiguration[8] = idctm;
                hbsA = SettingHbondConfg(Hbondconfiguration,countA,idatm,idctm,1,Natoms,atoms,geometry,ListNeighbourA,rAC,rAD);
                if (hbsA) {continue;}          //unknown constellation
                hbsB = SettingHbondConfg(Hbondconfiguration,countB,idbtm,idctm,5,Natoms,atoms,geometry,ListNeighbourB,rBC,rBD);
                if (hbsB) {continue;}          //unknown constellation
                theta = Angle(Hbondconfiguration[0] + 1,Hbondconfiguration[8] + 1,Hbondconfiguration[4] + 1,geometry);
                costheta = -cos(theta);
                if (costheta < 0.0) {continue;}          //in this case skip because H-bond is not physical
                dcostheta = 2.0*costheta*sin(theta);
                hbsA = HBondAngles(Hbondconfiguration[0],Hbondconfiguration[1],Hbondconfiguration[2],Hbondconfiguration[3],Hbondconfiguration[8],atoms,geometry,countA,cosphiA,dcosphiA,cospsiA,dcospsiA,torsionCheckA);
                if (hbsA) {continue;}                    //if unphysical, then move on
                hbsB = HBondAngles(Hbondconfiguration[4],Hbondconfiguration[5],Hbondconfiguration[6],Hbondconfiguration[7],Hbondconfiguration[8],atoms,geometry,countB,cosphiB,dcosphiB,cospsiB,dcospsiB,torsionCheckB);
                if (hbsB) {continue;}                    //if unphysical, then move on
                rAD = dist_Angstrom2au/rAB;   //1/xy_dist
                scaleA = CN*(atoms[idatm] == 7) + CO*(atoms[idatm] == 8);
                scaleB = CN*(atoms[idbtm] == 7) + CO*(atoms[idbtm] == 8);
                scaleC = 0.5*(scaleA + scaleB);
                rDonorH = fmin(rAC,rBC);
                AisDonor = (rDonorH == rAC);
                //now the damping factors
                scaleB = 1.0/1.2;
                fbond = 1.0 - 1.0/(1.0 + exp(-60.0*(rDonorH*scaleB - 1.0)));
                dampshort = 1.0/(1.0 + exp(-100.0*(rAB/2.4 - 1.0)));
                damplong = 1.0 - 1.0/(1.0 + exp(-10.0*(rAB/7.0 - 1.0)));
                damp = fbond*dampshort*damplong;
                dehdr = -2.00*scaleC*rAD*rAD*rAD*costheta*costheta*cosphiA*cosphiA*cospsiA*cospsiA*cosphiB*cosphiB*cospsiB*cospsiB*damp;
                dehdtheta = scaleC*rAD*rAD*dcostheta*cosphiA*cosphiA*cospsiA*cospsiA*cosphiB*cosphiB*cospsiB*cospsiB*damp;
                dehdphiA = scaleC*rAD*rAD*costheta*costheta*dcosphiA*cospsiA*cospsiA*cosphiB*cosphiB*cospsiB*cospsiB*damp;
                dehdphiB = scaleC*rAD*rAD*costheta*costheta*cosphiA*cosphiA*cospsiA*cospsiA*dcosphiB*cospsiB*cospsiB*damp;
                dehdpsiA = scaleC*rAD*rAD*costheta*costheta*cosphiA*cosphiA*dcospsiA*cosphiB*cosphiB*cospsiB*cospsiB*damp;
                dehdpsiB = scaleC*rAD*rAD*costheta*costheta*cosphiA*cosphiA*cospsiA*cospsiA*cosphiB*cosphiB*dcospsiB*damp;
                scaleA = exp(-60.0*(rDonorH*scaleB - 1.0));
                scaleB *= dist_Angstrom2au;
                dedh_aux = 60.0*scaleB*scaleA*(1.0/((1.0 + scaleA)*(1.0 + scaleA)));
                dedh = scaleC*rAD*rAD*costheta*costheta*cosphiA*cosphiA*cospsiA*cospsiA*cosphiB*cosphiB*cospsiB*cospsiB*dedh_aux*dampshort*damplong;
                scaleB = 1.0/2.4;
                scaleA = exp(-100.0*(rAB*scaleB - 1.0));
                scaleB *= dist_Angstrom2au;
                deds_aux = -100.0*scaleB*scaleA*(1.0/((1.0 + scaleA)*(1.0 + scaleA)));
                deds = scaleC*rAD*rAD*costheta*costheta*cosphiA*cosphiA*cospsiA*cospsiA*cosphiB*cosphiB*cospsiB*cospsiB*deds_aux*fbond*damplong;
                scaleB = 1.0/7.0;
                scaleA = exp(-10.0*(rAB*scaleB - 1.0));
                scaleB *= dist_Angstrom2au;
                dedl_aux = 10.0*scaleB*scaleA*(1.0/((1.0 + scaleA)*(1.0 + scaleA)));
                dedl = scaleC*rAD*rAD*costheta*costheta*cosphiA*cosphiA*cospsiA*cospsiA*cosphiB*cosphiB*cospsiB*cospsiB*dedl_aux*fbond*dampshort;
                dehdr += deds + dedl;
                //now convert from internal to Cartesian
                Internal2Cartesian(2,dehdr,Hbondconfiguration[0] + 1,Hbondconfiguration[4] + 1,0,0,geometry,gradient);
                Internal2Cartesian(3,dehdtheta,Hbondconfiguration[0] + 1,Hbondconfiguration[8] + 1,Hbondconfiguration[4] + 1,0,geometry,gradient);
                Internal2Cartesian(3,dehdphiA,Hbondconfiguration[1] + 1,Hbondconfiguration[0] + 1,Hbondconfiguration[8] + 1,0,geometry,gradient);
                Internal2Cartesian(3,dehdphiB, Hbondconfiguration[5] + 1,Hbondconfiguration[4] + 1,Hbondconfiguration[8] + 1,0,geometry,gradient);
                if (torsionCheckA) {Internal2Cartesian(4,-dehdpsiA,Hbondconfiguration[2] + 1,Hbondconfiguration[3] + 1,Hbondconfiguration[0] + 1,Hbondconfiguration[8] + 1,geometry,gradient);}
                else {Internal2Cartesian(4,-dehdpsiA,Hbondconfiguration[2] + 1,Hbondconfiguration[1] + 1,Hbondconfiguration[0] + 1,Hbondconfiguration[8] + 1,geometry,gradient);}
                if (torsionCheckB) {Internal2Cartesian(4,-dehdpsiB,Hbondconfiguration[6] + 1,Hbondconfiguration[7] + 1,Hbondconfiguration[4] + 1,Hbondconfiguration[8] + 1,geometry,gradient);}
                else {Internal2Cartesian(4,-dehdpsiB,Hbondconfiguration[6] + 1,Hbondconfiguration[5] + 1,Hbondconfiguration[4] + 1,Hbondconfiguration[8] + 1,geometry,gradient);}
                if (AisDonor) {Internal2Cartesian(2,dedh,Hbondconfiguration[8] + 1,Hbondconfiguration[0] + 1,0,0,geometry,gradient);}
                else {Internal2Cartesian(2,dedh,Hbondconfiguration[8] + 1,Hbondconfiguration[4] + 1,0,0,geometry,gradient);}
              }
            }
          }
        }
      }
    }
    //now the missing dCn/dCN terms
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
        //get distance between atoms
        rAB = Distance(idAtm + 1,idBtm + 1,geometry);
        fdamping = dNCoorddR(atoms[idAtm],atoms[idBtm],rAB,dEXPcountdx,false,1600.0,false);   //dCN/dRAB
        ged6 = fdamping*(dCn_terms[idAtm] + dCn_terms[idBtm]);
        Internal2Cartesian(2,-ged6,idAtm + 1,idBtm + 1,0,0,geometry,gradient);
      }
    }
  }
}

#endif //_Intermolecular_Potential_
