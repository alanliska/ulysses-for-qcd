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

#ifndef _APPROXIMATE_HESSIAN_
#define _APPROXIMATE_HESSIAN_
#include "math/FunctionPackage.hpp"
#include "math/MatrixPackage.hpp"
#include "parameters/Lindhpar.hpp"
#include "parameters/Schlegelpar.hpp"
#include "Molecule.hpp"

//description:
//methods to approximate Hessians to be used for geometry optimization and to correct the Hessian matrix by removing
//translational and rotational contributions

void LindhApproxHessian(matrixE & hessian, Molecule & mol, size_t parameterset = 1, double thresholdzero = 1.0e-14, double threshold_dist = 70.0, double thresholdr = 1.0e-10) {
  //calculation of approximated Hessian using the Lindh force field
  //R. Lindh, A. Bernhardsson, G. Karlstroem, P.-A. Malmqvist, Chem. Phys. Lett., 241(4), 423, 1995       -> parameter set 1
  //R. Lindh, personal communication - extracted from xTB                                                 -> parameter set 2
  //Hessian output in units of Eh/Angstrom^2
  std::vector<size_t> atoms = mol.Atoms();
  size_t Natoms = atoms.size();
  hessian.resize(3*Natoms,3*Natoms);
  hessian.zero();
  double rAB;
  double rAB2;
  double rCA;
  double rCA2;
  double rBC;
  double rBC2;
  double rAD2;
  double rBD2;
  double rCD2;
  double rAB0;
  double rBD0;
  double rCA0;
  double aAB;
  double aBD;
  double aCA;
  double cosf2;
  double cosf3;
  double kr = LindhRKr(parameterset);
  double kf = LindhRKf(parameterset);
  double kt = LindhRKt(parameterset);
  double zeroangle;
  double alpha;
  double RL;
  double sinp;
  double cosp;
  double rABdotrCA;
  double RAB[3];
  double RCA[3];
  double RBC[3];
  double RBD[3];
  double sA[3];
  double sB[3];
  double sC[3];
  double sD[3];
  double xx[2];
  double yy[2];
  double zz[2];
  double cosphimax = cos(35.0*pi/180.0);
  double Fstretch;
  double Fbend;
  double Ftorsion;
  double unitconv = dist_Angstrom2aum1*dist_Angstrom2aum1;
  matrixE Bmat;
  matrixE xyz(3,4);
  matrixE geometry = mol.Geometry();
  std::vector<matrixE> dBmat;
  int ACpair;
  int BDpair;
  for (size_t iAtm = 0; iAtm < Natoms; ++iAtm) {
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      xyz(idcoord + 1,2) = geometry(iAtm + 1,idcoord + 1)*dist_Angstrom2aum1;
    }
    for (size_t iBtm = 0; iBtm < Natoms; ++iBtm) {
      if (iBtm == iAtm) {continue;}
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        RAB[idcoord] = (geometry(iBtm + 1,idcoord + 1) - geometry(iAtm + 1,idcoord + 1))*dist_Angstrom2aum1;
        xyz(idcoord + 1,3) = geometry(iBtm + 1,idcoord + 1)*dist_Angstrom2aum1;
      }
      rAB2 = RAB[0]*RAB[0] + RAB[1]*RAB[1] + RAB[2]*RAB[2];
      rAB0 = LindhRAV(atoms[iAtm],atoms[iBtm],parameterset);
      aAB = LindhAAV(atoms[iAtm],atoms[iBtm],parameterset);
      if (iBtm < iAtm) {
        //stretching
        Fstretch = kr*exp(aAB*(rAB0*rAB0 - rAB2));
        //converting to cartesian Hessian
        for (size_t idcd1 = 0; idcd1 < 3; ++idcd1) {
          for (size_t idcd2 = 0; idcd2 < 3; ++idcd2) {
            hessian(3*iAtm + 1 + idcd1,3*iAtm + 1 + idcd2) += Fstretch*RAB[idcd1]*RAB[idcd2]*unitconv/rAB2;
            hessian(3*iBtm + 1 + idcd1,3*iBtm + 1 + idcd2) += Fstretch*RAB[idcd1]*RAB[idcd2]*unitconv/rAB2;
            hessian(3*iAtm + 1 + idcd1,3*iBtm + 1 + idcd2) -= Fstretch*RAB[idcd1]*RAB[idcd2]*unitconv/rAB2;
            hessian(3*iBtm + 1 + idcd1,3*iAtm + 1 + idcd2) -= Fstretch*RAB[idcd1]*RAB[idcd2]*unitconv/rAB2;
          }
        }
      }
      if (rAB2 > threshold_dist) {continue;}
      rAB = sqrt(rAB2);
      //bending
      for (size_t iCtm = 0; iCtm < iBtm; ++iCtm) {
        if (iCtm == iAtm) {continue;}
        zeroangle = 0.0;
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          RCA[idcoord] = (geometry(iCtm + 1,idcoord + 1) - geometry(iAtm + 1,idcoord + 1))*dist_Angstrom2aum1;
          RBC[idcoord] = (geometry(iBtm + 1,idcoord + 1) - geometry(iCtm + 1,idcoord + 1))*dist_Angstrom2aum1;
          zeroangle += RAB[idcoord]*RCA[idcoord];
        }
        rCA2 = RCA[0]*RCA[0] + RCA[1]*RCA[1] + RCA[2]*RCA[2];
        rBC2 = RBC[0]*RBC[0] + RBC[1]*RBC[1] + RBC[2]*RBC[2];
        if (rCA2 > threshold_dist) {continue;}
        if (rBC2 > threshold_dist) {continue;}
        rCA0 = LindhRAV(atoms[iCtm],atoms[iAtm],parameterset);
        aCA = LindhAAV(atoms[iCtm],atoms[iAtm],parameterset);
        //check for zero angle
        zeroangle /= sqrt(rAB2*rCA2);
        if (fabs(zeroangle - 1.0) < thresholdzero) {continue;}
        //move this inside the next if statement
        zeroangle = (RAB[1]*RCA[2] - RAB[2]*RCA[1])*(RAB[1]*RCA[2] - RAB[2]*RCA[1]);
        zeroangle += (RAB[2]*RCA[0] - RAB[0]*RCA[2])*(RAB[2]*RCA[0] - RAB[0]*RCA[2]);
        zeroangle += (RAB[0]*RCA[1] - RAB[1]*RCA[0])*(RAB[0]*RCA[1] - RAB[1]*RCA[0]);
        if (zeroangle < thresholdzero) {RL = 0.0;}
        else {RL = sqrt(zeroangle);}
        rCA = sqrt(rCA2);
        rBC = sqrt(rBC2);
        if ((rAB > thresholdr)&&(rCA > thresholdr)&&(rBC > thresholdr)) {
          sinp = RL/(rCA*rAB);
          rABdotrCA = RAB[0]*RCA[0] + RAB[1]*RCA[1] + RAB[2]*RCA[2];
          cosp = rABdotrCA/(rCA*rAB);
          alpha = kf*exp(aAB*rAB0*rAB0 + aCA*rCA0*rCA0);
          Fbend = alpha*exp(-(aAB*rAB2 + aCA*rCA2));
          if (sinp > thresholdr) {              //non-linear case
            for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
              sB[idcoord] = (cosp*RAB[idcoord]/rAB - RCA[idcoord]/rCA)/(rAB*sinp);
              sC[idcoord] = (cosp*RCA[idcoord]/rCA - RAB[idcoord]/rAB)/(rCA*sinp);
              sA[idcoord] = -sB[idcoord] - sC[idcoord];
            }
            for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
              hessian(3*iBtm + idcoord + 1,3*iBtm + idcoord + 1) += Fbend*sB[idcoord]*sB[idcoord]*unitconv;
              hessian(3*iCtm + idcoord + 1,3*iCtm + idcoord + 1) += Fbend*sC[idcoord]*sC[idcoord]*unitconv;
              hessian(3*iAtm + idcoord + 1,3*iAtm + idcoord + 1) += Fbend*sA[idcoord]*sA[idcoord]*unitconv;
              for (size_t jdcoord = 0; jdcoord < 3; ++jdcoord) {
                if (iAtm > iBtm) {
                  hessian(3*iAtm + idcoord + 1,3*iBtm + jdcoord + 1) += Fbend*sA[idcoord]*sB[jdcoord]*unitconv;
                  hessian(3*iBtm + jdcoord + 1,3*iAtm + idcoord + 1) += Fbend*sA[idcoord]*sB[jdcoord]*unitconv;
                }
                else {
                  hessian(3*iBtm + idcoord + 1,3*iAtm + jdcoord + 1) += Fbend*sB[idcoord]*sA[jdcoord]*unitconv;
                  hessian(3*iAtm + jdcoord + 1,3*iBtm + idcoord + 1) += Fbend*sB[idcoord]*sA[jdcoord]*unitconv;
                }
                if (iAtm > iCtm) {
                  hessian(3*iAtm + idcoord + 1,3*iCtm + jdcoord + 1) += Fbend*sA[idcoord]*sC[jdcoord]*unitconv;
                  hessian(3*iCtm + jdcoord + 1,3*iAtm + idcoord + 1) += Fbend*sA[idcoord]*sC[jdcoord]*unitconv;
                }
                else {
                  hessian(3*iCtm + idcoord + 1,3*iAtm + jdcoord + 1) += Fbend*sC[idcoord]*sA[jdcoord]*unitconv;
                  hessian(3*iAtm + jdcoord + 1,3*iCtm + idcoord + 1) += Fbend*sC[idcoord]*sA[jdcoord]*unitconv;
                }
                if (iBtm > iCtm) {
                  hessian(3*iBtm + idcoord + 1,3*iCtm + jdcoord + 1) += Fbend*sB[idcoord]*sC[jdcoord]*unitconv;
                  hessian(3*iCtm + jdcoord + 1,3*iBtm + idcoord + 1) += Fbend*sB[idcoord]*sC[jdcoord]*unitconv;
                }
                //this below should never happen!!!!!
                else {
                  hessian(3*iCtm + idcoord + 1,3*iBtm + jdcoord + 1) += Fbend*sC[idcoord]*sB[jdcoord]*unitconv;
                  hessian(3*iBtm + jdcoord + 1,3*iCtm + idcoord + 1) += Fbend*sC[idcoord]*sB[jdcoord]*unitconv;
                }
                if (jdcoord < idcoord) {
                  hessian(3*iBtm + idcoord + 1,3*iBtm + jdcoord + 1) += Fbend*sB[idcoord]*sB[jdcoord]*unitconv;
                  hessian(3*iBtm + jdcoord + 1,3*iBtm + idcoord + 1) += Fbend*sB[idcoord]*sB[jdcoord]*unitconv;
                  hessian(3*iCtm + idcoord + 1,3*iCtm + jdcoord + 1) += Fbend*sC[idcoord]*sC[jdcoord]*unitconv;
                  hessian(3*iCtm + jdcoord + 1,3*iCtm + idcoord + 1) += Fbend*sC[idcoord]*sC[jdcoord]*unitconv;
                  hessian(3*iAtm + idcoord + 1,3*iAtm + jdcoord + 1) += Fbend*sA[idcoord]*sA[jdcoord]*unitconv;
                  hessian(3*iAtm + jdcoord + 1,3*iAtm + idcoord + 1) += Fbend*sA[idcoord]*sA[jdcoord]*unitconv;
                }
              }
            }
          }
          else {                                //linear case
            if ((fabs(RAB[1]) > thresholdr)||(fabs(RAB[0]) > thresholdr)) {
              xx[0] = -RAB[1];
              yy[0] = RAB[0];
              zz[0] = 0.0;
              xx[1] = -RAB[0]*RAB[2];
              yy[1] = -RAB[1]*RAB[2];
              zz[1] = RAB[0]*RAB[0] + RAB[1]*RAB[1];
            }
            else {
              xx[0] = 1.0;
              yy[0] = 0.0;
              zz[0] = 0.0;
              xx[1] = 0.0;
              yy[1] = 1.0;
              zz[1] = 0.0;
            }
            for (size_t idx = 0; idx < 2; ++idx) {
              RL = 1.0/sqrt(xx[idx]*xx[idx] + yy[idx]*yy[idx] + zz[idx]*zz[idx]);
              sB[0] = -xx[idx]*RL/rAB;
              sC[0] = -xx[idx]*RL/rCA;
              sB[1] = -yy[idx]*RL/rAB;
              sC[1] = -yy[idx]*RL/rCA;
              sB[2] = -zz[idx]*RL/rAB;
              sC[2] = -zz[idx]*RL/rCA;
              for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                sA[idcoord] = -sB[idcoord] - sC[idcoord];
              }
              for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                hessian(3*iBtm + idcoord + 1,3*iBtm + idcoord + 1) += Fbend*sB[idcoord]*sB[idcoord]*unitconv;
                hessian(3*iCtm + idcoord + 1,3*iCtm + idcoord + 1) += Fbend*sC[idcoord]*sC[idcoord]*unitconv;
                hessian(3*iAtm + idcoord + 1,3*iAtm + idcoord + 1) += Fbend*sA[idcoord]*sA[idcoord]*unitconv;
                for (size_t jdcoord = 0; jdcoord < 3; ++jdcoord) {
                  if (iAtm > iBtm) {
                    hessian(3*iAtm + idcoord + 1,3*iBtm + jdcoord + 1) += Fbend*sA[idcoord]*sB[jdcoord]*unitconv;
                    hessian(3*iBtm + jdcoord + 1,3*iAtm + idcoord + 1) += Fbend*sA[idcoord]*sB[jdcoord]*unitconv;
                  }
                  else {
                    hessian(3*iBtm + idcoord + 1,3*iAtm + jdcoord + 1) += Fbend*sB[idcoord]*sA[jdcoord]*unitconv;
                    hessian(3*iAtm + jdcoord + 1,3*iBtm + idcoord + 1) += Fbend*sB[idcoord]*sA[jdcoord]*unitconv;
                  }
                  if (iAtm > iCtm) {
                    hessian(3*iAtm + idcoord + 1,3*iCtm + jdcoord + 1) += Fbend*sA[idcoord]*sC[jdcoord]*unitconv;
                    hessian(3*iCtm + jdcoord + 1,3*iAtm + idcoord + 1) += Fbend*sA[idcoord]*sC[jdcoord]*unitconv;
                  }
                  else {
                    hessian(3*iCtm + idcoord + 1,3*iAtm + jdcoord + 1) += Fbend*sC[idcoord]*sA[jdcoord]*unitconv;
                    hessian(3*iAtm + jdcoord + 1,3*iCtm + idcoord + 1) += Fbend*sC[idcoord]*sA[jdcoord]*unitconv;
                  }
                  if (iBtm > iCtm) {
                    hessian(3*iBtm + idcoord + 1,3*iCtm + jdcoord + 1) += Fbend*sB[idcoord]*sC[jdcoord]*unitconv;
                    hessian(3*iCtm + jdcoord + 1,3*iBtm + idcoord + 1) += Fbend*sB[idcoord]*sC[jdcoord]*unitconv;
                  }
                  //this below should never happen!!!!!
                  else {
                    hessian(3*iCtm + idcoord + 1,3*iBtm + jdcoord + 1) += Fbend*sC[idcoord]*sB[jdcoord]*unitconv;
                    hessian(3*iBtm + jdcoord + 1,3*iCtm + idcoord + 1) += Fbend*sC[idcoord]*sB[jdcoord]*unitconv;
                  }
                  if (jdcoord < idcoord) {
                    hessian(3*iBtm + idcoord + 1,3*iBtm + jdcoord + 1) += Fbend*sB[idcoord]*sB[jdcoord]*unitconv;
                    hessian(3*iBtm + jdcoord + 1,3*iBtm + idcoord + 1) += Fbend*sB[idcoord]*sB[jdcoord]*unitconv;
                    hessian(3*iCtm + idcoord + 1,3*iCtm + jdcoord + 1) += Fbend*sC[idcoord]*sC[jdcoord]*unitconv;
                    hessian(3*iCtm + jdcoord + 1,3*iCtm + idcoord + 1) += Fbend*sC[idcoord]*sC[jdcoord]*unitconv;
                    hessian(3*iAtm + idcoord + 1,3*iAtm + jdcoord + 1) += Fbend*sA[idcoord]*sA[jdcoord]*unitconv;
                    hessian(3*iAtm + jdcoord + 1,3*iAtm + idcoord + 1) += Fbend*sA[idcoord]*sA[jdcoord]*unitconv;
                  }
                }
              }
            }
          }
        }
      }
      //torsion
      for (size_t iCtm = 0; iCtm < Natoms; ++iCtm) {
        ACpair = Natoms*iAtm + iCtm;
        if (iCtm == iAtm) {continue;}
        if (iCtm == iBtm) {continue;}
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          xyz(idcoord + 1,1) = geometry(iCtm + 1,idcoord + 1)*dist_Angstrom2aum1;
          RCA[idcoord] = (geometry(iCtm + 1,idcoord + 1) - geometry(iAtm + 1,idcoord + 1))*dist_Angstrom2aum1;
          RBC[idcoord] = (geometry(iBtm + 1,idcoord + 1) - geometry(iCtm + 1,idcoord + 1))*dist_Angstrom2aum1;
        }
        rCA2 = RCA[0]*RCA[0] + RCA[1]*RCA[1] + RCA[2]*RCA[2];
        rBC2 = RBC[0]*RBC[0] + RBC[1]*RBC[1] + RBC[2]*RBC[2];
        if (rBC2 > threshold_dist) {continue;}
        if (rCA2 > threshold_dist) {continue;}
        for (size_t iDtm = 0; iDtm < Natoms; ++iDtm) {
          BDpair = Natoms*iBtm + iDtm;
          if (ACpair <= BDpair) {continue;}
          if (iDtm == iCtm) {continue;}
          if (iDtm == iAtm) {continue;}
          if (iDtm == iBtm) {continue;}
          rAD2 = 0.0;
          rCD2 = 0.0;
          for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
            xyz(idcoord + 1,4) = geometry(iDtm + 1,idcoord + 1)*dist_Angstrom2aum1;
            RBD[idcoord] = (geometry(iBtm + 1,idcoord + 1) - geometry(iDtm + 1,idcoord + 1))*dist_Angstrom2aum1;
            zeroangle = (geometry(iAtm + 1,idcoord + 1) - geometry(iDtm + 1,idcoord + 1))*dist_Angstrom2aum1;
            rAD2 += zeroangle*zeroangle;
            zeroangle = (geometry(iCtm + 1,idcoord + 1) - geometry(iDtm + 1,idcoord + 1))*dist_Angstrom2aum1;
            rCD2 += zeroangle*zeroangle;
          }
          rBD2 = RBD[0]*RBD[0] + RBD[1]*RBD[1] + RBD[2]*RBD[2];
          if (rCD2 > threshold_dist) {continue;}
          if (rAD2 > threshold_dist) {continue;}
          if (rBD2 > threshold_dist) {continue;}
          zeroangle = LindhRAV(atoms[iAtm],atoms[iBtm],parameterset);
          rAB0 = zeroangle*zeroangle;
          zeroangle = LindhRAV(atoms[iBtm],atoms[iDtm],parameterset);
          rBD0 = zeroangle*zeroangle;
          zeroangle = LindhRAV(atoms[iCtm],atoms[iAtm],parameterset);
          rCA0 = zeroangle*zeroangle;
          aAB = LindhAAV(atoms[iAtm],atoms[iBtm],parameterset);
          aBD = LindhAAV(atoms[iBtm],atoms[iDtm],parameterset);
          aCA = LindhAAV(atoms[iCtm],atoms[iAtm],parameterset);
          cosf2 = -(RCA[0]*RAB[0] + RCA[1]*RAB[1] + RCA[2]*RAB[2])/sqrt(rCA2*rAB2);
          if (fabs(cosf2) > cosphimax) {continue;}
          cosf3 = -(RBD[0]*RAB[0] + RBD[1]*RAB[1] + RBD[2]*RAB[2])/sqrt(rBD2*rAB2);
          if (fabs(cosf3) > cosphimax) {continue;}
          alpha = kt*exp(aCA*rCA0 + aAB*rAB0 + aBD*rBD0);
          Ftorsion = alpha*exp(-(aCA*rCA2 + aAB*rAB2 + aBD*rBD2));
          Torsion(xyz,alpha,Bmat,dBmat,false);
          for (size_t idx = 0; idx < 3; ++idx) {
            sC[idx] = Bmat(idx + 1,1);
            sA[idx] = Bmat(idx + 1,2);
            sB[idx] = Bmat(idx + 1,3);
            sD[idx] = Bmat(idx + 1,4);
          }
          for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
            hessian(3*iCtm + idcoord + 1,3*iCtm + idcoord + 1) += Ftorsion*sC[idcoord]*sC[idcoord]*unitconv;
            hessian(3*iAtm + idcoord + 1,3*iAtm + idcoord + 1) += Ftorsion*sA[idcoord]*sA[idcoord]*unitconv;
            hessian(3*iBtm + idcoord + 1,3*iBtm + idcoord + 1) += Ftorsion*sB[idcoord]*sB[idcoord]*unitconv;
            hessian(3*iDtm + idcoord + 1,3*iDtm + idcoord + 1) += Ftorsion*sD[idcoord]*sD[idcoord]*unitconv;
            for (size_t jdcoord = 0; jdcoord < 3; ++jdcoord) {
              hessian(3*iCtm + idcoord + 1,3*iAtm + jdcoord + 1) += Ftorsion*sC[idcoord]*sA[jdcoord]*unitconv;
              hessian(3*iAtm + jdcoord + 1,3*iCtm + idcoord + 1) += Ftorsion*sC[idcoord]*sA[jdcoord]*unitconv;
              hessian(3*iCtm + idcoord + 1,3*iBtm + jdcoord + 1) += Ftorsion*sC[idcoord]*sB[jdcoord]*unitconv;
              hessian(3*iBtm + jdcoord + 1,3*iCtm + idcoord + 1) += Ftorsion*sC[idcoord]*sB[jdcoord]*unitconv;
              hessian(3*iCtm + idcoord + 1,3*iDtm + jdcoord + 1) += Ftorsion*sC[idcoord]*sD[jdcoord]*unitconv;
              hessian(3*iDtm + jdcoord + 1,3*iCtm + idcoord + 1) += Ftorsion*sC[idcoord]*sD[jdcoord]*unitconv;
              hessian(3*iAtm + idcoord + 1,3*iBtm + jdcoord + 1) += Ftorsion*sA[idcoord]*sB[jdcoord]*unitconv;
              hessian(3*iBtm + jdcoord + 1,3*iAtm + idcoord + 1) += Ftorsion*sA[idcoord]*sB[jdcoord]*unitconv;
              hessian(3*iAtm + idcoord + 1,3*iDtm + jdcoord + 1) += Ftorsion*sA[idcoord]*sD[jdcoord]*unitconv;
              hessian(3*iDtm + jdcoord + 1,3*iAtm + idcoord + 1) += Ftorsion*sA[idcoord]*sD[jdcoord]*unitconv;
              hessian(3*iBtm + idcoord + 1,3*iDtm + jdcoord + 1) += Ftorsion*sB[idcoord]*sD[jdcoord]*unitconv;
              hessian(3*iDtm + jdcoord + 1,3*iBtm + idcoord + 1) += Ftorsion*sB[idcoord]*sD[jdcoord]*unitconv;
              if (jdcoord < idcoord) {
                hessian(3*iCtm + idcoord + 1,3*iCtm + jdcoord + 1) += Ftorsion*sC[idcoord]*sC[jdcoord]*unitconv;
                hessian(3*iCtm + jdcoord + 1,3*iCtm + idcoord + 1) += Ftorsion*sC[idcoord]*sC[jdcoord]*unitconv;
                hessian(3*iAtm + idcoord + 1,3*iAtm + jdcoord + 1) += Ftorsion*sA[idcoord]*sA[jdcoord]*unitconv;
                hessian(3*iAtm + jdcoord + 1,3*iAtm + idcoord + 1) += Ftorsion*sA[idcoord]*sA[jdcoord]*unitconv;
                hessian(3*iBtm + idcoord + 1,3*iBtm + jdcoord + 1) += Ftorsion*sB[idcoord]*sB[jdcoord]*unitconv;
                hessian(3*iBtm + jdcoord + 1,3*iBtm + idcoord + 1) += Ftorsion*sB[idcoord]*sB[jdcoord]*unitconv;
                hessian(3*iDtm + idcoord + 1,3*iDtm + jdcoord + 1) += Ftorsion*sD[idcoord]*sD[jdcoord]*unitconv;
                hessian(3*iDtm + jdcoord + 1,3*iDtm + idcoord + 1) += Ftorsion*sD[idcoord]*sD[jdcoord]*unitconv;
              }
            }
          }
        }
      }
    }
  }
  //since there might be numerical inaccuracies in the the symmetry, symmetrize
  for (size_t idr = 0; idr < 3*Natoms; ++idr) {
    for (size_t idc = 0; idc < idr; ++idc) {
      Fstretch = 0.5*(hessian(idr + 1,idc + 1) + hessian(idc + 1,idr + 1));
      hessian(idr + 1,idc + 1) = Fstretch;
      hessian(idc + 1,idr + 1) = Fstretch;
    }
  }
}
void SchlegelApproxHessian(matrixE & hessian, Molecule & mol, size_t parameterset = 1, double threshold = 1.0e-7, double thresholdzero = 1.0e-14, double threshold_dist = 70.0, double thresholdr = 1.0e-10) {
  //calculation of approximated Hessian using the Schlegel force field
  //H. B. Schlegel, Theoret. Chim. Acta (Berl.), 66, 333, 1984                      -> parameter set 1
  //J. M. Wittbrodt, H. B. Schlegel, J. Mol. Struct. (Theochem), 398, 55, 1997      -> parameter set 2
  //Hessian output in units of Eh/Angstrom^2
  std::vector<size_t> atoms = mol.Atoms();
  size_t Natoms = atoms.size();
  hessian.resize(3*Natoms,3*Natoms);
  hessian.zero();
  double aux;
  double aux2;
  double rAB;
  double rAB2;
  double rCA;
  double rCA2;
  double rBC;
  double rBC2;
  double rAD2;
  double rBD2;
  double rCD2;
  double rAB0;
  double rAC0;
  double rBD0;
  double cosf2;
  double cosf3;
  double cosf4;
  double zeroangle;
  double alpha;
  double RL;
  double sinp;
  double cosp;
  double rABdotrCA;
  double unitconv = dist_Angstrom2aum1*dist_Angstrom2aum1;
  double cosphimax = cos(35.0*pi/180.0);
  double taumax = 45.0*pi/180.0;
  double Fstretch;
  double Fbend;
  double Ftorsion;
  double Foop;
  double RAB[3];
  double RCA[3];
  double RBC[3];
  double RBD[3];
  double RAD[3];
  double sA[3];
  double sB[3];
  double sC[3];
  double sD[3];
  double xx[2];
  double yy[2];
  double zz[2];
  double cross[3];
  matrixE geometry = mol.Geometry();
  matrixE Bmat;
  matrixE xyz(3,4);
  matrixE oxyz(3,4);
  std::vector<matrixE> dBmat;
  int ACpair;
  int BDpair;
  bool cosf2fail;
  bool cosf3fail;
  bool cosf4fail;
  for (size_t iAtm = 0; iAtm < Natoms; ++iAtm) {
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      xyz(idcoord + 1,2) = geometry(iAtm + 1,idcoord + 1)*dist_Angstrom2aum1;
      oxyz(idcoord + 1,4) = geometry(iAtm + 1,idcoord + 1)*dist_Angstrom2aum1;
    }
    for (size_t iBtm = 0; iBtm < Natoms; ++iBtm) {
      if (iBtm == iAtm) {continue;}
      rAB0 = 1.30*(AtmRadii(atoms[iAtm]) + AtmRadii(atoms[iBtm]))*dist_Angstrom2aum1;
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        RAB[idcoord] = (geometry(iBtm + 1,idcoord + 1) - geometry(iAtm + 1,idcoord + 1))*dist_Angstrom2aum1;
        xyz(idcoord + 1,3) = geometry(iBtm + 1,idcoord + 1)*dist_Angstrom2aum1;
        oxyz(idcoord + 1,1) = geometry(iBtm + 1,idcoord + 1)*dist_Angstrom2aum1;
      }
      rAB2 = RAB[0]*RAB[0] + RAB[1]*RAB[1] + RAB[2]*RAB[2];
      rAB = sqrt(rAB2);
      if (rAB > rAB0) {continue;}
      if (iBtm < iAtm) {
        //stretching
        aux2 = (rAB - SchlegelB(atoms[iAtm],atoms[iBtm],parameterset));
        aux = aux2*aux2*aux2;
        Fstretch = 2.2542/aux;
        for (size_t idcd1 = 0; idcd1 < 3; ++idcd1) {
          for (size_t idcd2 = 0; idcd2 < 3; ++idcd2) {
            hessian(3*iAtm + 1 + idcd1,3*iAtm + 1 + idcd2) += Fstretch*RAB[idcd1]*RAB[idcd2]*unitconv/rAB2;
            hessian(3*iBtm + 1 + idcd1,3*iBtm + 1 + idcd2) += Fstretch*RAB[idcd1]*RAB[idcd2]*unitconv/rAB2;
            hessian(3*iAtm + 1 + idcd1,3*iBtm + 1 + idcd2) -= Fstretch*RAB[idcd1]*RAB[idcd2]*unitconv/rAB2;
            hessian(3*iBtm + 1 + idcd1,3*iAtm + 1 + idcd2) -= Fstretch*RAB[idcd1]*RAB[idcd2]*unitconv/rAB2;
          }
        }
      }
      if (rAB2 > threshold_dist) {continue;}
      //bending
      for (size_t iCtm = 0; iCtm < iBtm; ++iCtm) {
        if (iCtm == iAtm) {continue;}
        rAC0 = 1.30*(AtmRadii(atoms[iAtm]) + AtmRadii(atoms[iCtm]))*dist_Angstrom2aum1;
        zeroangle = 0.0;
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          RCA[idcoord] = (geometry(iCtm + 1,idcoord + 1) - geometry(iAtm + 1,idcoord + 1))*dist_Angstrom2aum1;
          RBC[idcoord] = (geometry(iBtm + 1,idcoord + 1) - geometry(iCtm + 1,idcoord + 1))*dist_Angstrom2aum1;
          zeroangle += RAB[idcoord]*RCA[idcoord];
        }
        rCA2 = RCA[0]*RCA[0] + RCA[1]*RCA[1] + RCA[2]*RCA[2];
        rBC2 = RBC[0]*RBC[0] + RBC[1]*RBC[1] + RBC[2]*RBC[2];
        if (rCA2 > threshold_dist) {continue;}
        if (rBC2 > threshold_dist) {continue;}
        if (rCA2 > rAC0*rAC0) {continue;}
        //check for zero angle
        zeroangle /= sqrt(rAB2*rCA2);
        if (fabs(zeroangle - 1.0) < thresholdzero) {continue;}
        zeroangle = (RAB[1]*RCA[2] - RAB[2]*RCA[1])*(RAB[1]*RCA[2] - RAB[2]*RCA[1]);
        zeroangle += (RAB[2]*RCA[0] - RAB[0]*RCA[2])*(RAB[2]*RCA[0] - RAB[0]*RCA[2]);
        zeroangle += (RAB[0]*RCA[1] - RAB[1]*RCA[0])*(RAB[0]*RCA[1] - RAB[1]*RCA[0]);
        if (zeroangle < thresholdzero) {RL = 0.0;}
        else {RL = sqrt(zeroangle);}
        rCA = sqrt(rCA2);
        rBC = sqrt(rBC2);
        if ((rAB > thresholdr)&&(rCA > thresholdr)&&(rBC > thresholdr)) {
          sinp = RL/(rCA*rAB);
          rABdotrCA = RAB[0]*RCA[0] + RAB[1]*RCA[1] + RAB[2]*RCA[2];
          cosp = rABdotrCA/(rCA*rAB);
          Fbend = 0.250;
          if (sinp > thresholdr) {              //non-linear case
            for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
              sB[idcoord] = (cosp*RAB[idcoord]/rAB - RCA[idcoord]/rCA)/(rAB*sinp);
              sC[idcoord] = (cosp*RCA[idcoord]/rCA - RAB[idcoord]/rAB)/(rCA*sinp);
              sA[idcoord] = -sB[idcoord] - sC[idcoord];
            }
            for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
              hessian(3*iBtm + idcoord + 1,3*iBtm + idcoord + 1) += Fbend*sB[idcoord]*sB[idcoord]*unitconv;
              hessian(3*iCtm + idcoord + 1,3*iCtm + idcoord + 1) += Fbend*sC[idcoord]*sC[idcoord]*unitconv;
              hessian(3*iAtm + idcoord + 1,3*iAtm + idcoord + 1) += Fbend*sA[idcoord]*sA[idcoord]*unitconv;
              for (size_t jdcoord = 0; jdcoord < 3; ++jdcoord) {
                if (iAtm > iBtm) {
                  hessian(3*iAtm + idcoord + 1,3*iBtm + jdcoord + 1) += Fbend*sA[idcoord]*sB[jdcoord]*unitconv;
                  hessian(3*iBtm + jdcoord + 1,3*iAtm + idcoord + 1) += Fbend*sA[idcoord]*sB[jdcoord]*unitconv;
                }
                else {
                  hessian(3*iBtm + idcoord + 1,3*iAtm + jdcoord + 1) += Fbend*sB[idcoord]*sA[jdcoord]*unitconv;
                  hessian(3*iAtm + jdcoord + 1,3*iBtm + idcoord + 1) += Fbend*sB[idcoord]*sA[jdcoord]*unitconv;
                }
                if (iAtm > iCtm) {
                  hessian(3*iAtm + idcoord + 1,3*iCtm + jdcoord + 1) += Fbend*sA[idcoord]*sC[jdcoord]*unitconv;
                  hessian(3*iCtm + jdcoord + 1,3*iAtm + idcoord + 1) += Fbend*sA[idcoord]*sC[jdcoord]*unitconv;
                }
                else {
                  hessian(3*iCtm + idcoord + 1,3*iAtm + jdcoord + 1) += Fbend*sC[idcoord]*sA[jdcoord]*unitconv;
                  hessian(3*iAtm + jdcoord + 1,3*iCtm + idcoord + 1) += Fbend*sC[idcoord]*sA[jdcoord]*unitconv;
                }
                if (iBtm > iCtm) {
                  hessian(3*iBtm + idcoord + 1,3*iCtm + jdcoord + 1) += Fbend*sB[idcoord]*sC[jdcoord]*unitconv;
                  hessian(3*iCtm + jdcoord + 1,3*iBtm + idcoord + 1) += Fbend*sB[idcoord]*sC[jdcoord]*unitconv;
                }
                //this below should never happen!!!!!
                else {
                  hessian(3*iCtm + idcoord + 1,3*iBtm + jdcoord + 1) += Fbend*sC[idcoord]*sB[jdcoord]*unitconv;
                  hessian(3*iBtm + jdcoord + 1,3*iCtm + idcoord + 1) += Fbend*sC[idcoord]*sB[jdcoord]*unitconv;
                }
                if (jdcoord < idcoord) {
                  hessian(3*iBtm + idcoord + 1,3*iBtm + jdcoord + 1) += Fbend*sB[idcoord]*sB[jdcoord]*unitconv;
                  hessian(3*iBtm + jdcoord + 1,3*iBtm + idcoord + 1) += Fbend*sB[idcoord]*sB[jdcoord]*unitconv;
                  hessian(3*iCtm + idcoord + 1,3*iCtm + jdcoord + 1) += Fbend*sC[idcoord]*sC[jdcoord]*unitconv;
                  hessian(3*iCtm + jdcoord + 1,3*iCtm + idcoord + 1) += Fbend*sC[idcoord]*sC[jdcoord]*unitconv;
                  hessian(3*iAtm + idcoord + 1,3*iAtm + jdcoord + 1) += Fbend*sA[idcoord]*sA[jdcoord]*unitconv;
                  hessian(3*iAtm + jdcoord + 1,3*iAtm + idcoord + 1) += Fbend*sA[idcoord]*sA[jdcoord]*unitconv;
                }
              }
            }
          }
          else {                                //linear case
            if ((fabs(RAB[1]) > thresholdr)||(fabs(RAB[0]) > thresholdr)) {
              xx[0] = -RAB[1];
              yy[0] = RAB[0];
              zz[0] = 0.0;
              xx[1] = -RAB[0]*RAB[2];
              yy[1] = -RAB[1]*RAB[2];
              zz[1] = RAB[0]*RAB[0] + RAB[1]*RAB[1];
            }
            else {
              xx[0] = 1.0;
              yy[0] = 0.0;
              zz[0] = 0.0;
              xx[1] = 0.0;
              yy[1] = 1.0;
              zz[1] = 0.0;
            }
            for (size_t idx = 0; idx < 2; ++idx) {
              RL = 1.0/sqrt(xx[idx]*xx[idx] + yy[idx]*yy[idx] + zz[idx]*zz[idx]);
              sB[0] = -xx[idx]*RL/rAB;
              sC[0] = -xx[idx]*RL/rCA;
              sB[1] = -yy[idx]*RL/rAB;
              sC[1] = -yy[idx]*RL/rCA;
              sB[2] = -zz[idx]*RL/rAB;
              sC[2] = -zz[idx]*RL/rCA;
              for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                sA[idcoord] = -sB[idcoord] - sC[idcoord];
              }
              for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
                hessian(3*iBtm + idcoord + 1,3*iBtm + idcoord + 1) += Fbend*sB[idcoord]*sB[idcoord]*unitconv;
                hessian(3*iCtm + idcoord + 1,3*iCtm + idcoord + 1) += Fbend*sC[idcoord]*sC[idcoord]*unitconv;
                hessian(3*iAtm + idcoord + 1,3*iAtm + idcoord + 1) += Fbend*sA[idcoord]*sA[idcoord]*unitconv;
                for (size_t jdcoord = 0; jdcoord < 3; ++jdcoord) {
                  if (iAtm > iBtm) {
                    hessian(3*iAtm + idcoord + 1,3*iBtm + jdcoord + 1) += Fbend*sA[idcoord]*sB[jdcoord]*unitconv;
                    hessian(3*iBtm + jdcoord + 1,3*iAtm + idcoord + 1) += Fbend*sA[idcoord]*sB[jdcoord]*unitconv;
                  }
                  else {
                    hessian(3*iBtm + idcoord + 1,3*iAtm + jdcoord + 1) += Fbend*sB[idcoord]*sA[jdcoord]*unitconv;
                    hessian(3*iAtm + jdcoord + 1,3*iBtm + idcoord + 1) += Fbend*sB[idcoord]*sA[jdcoord]*unitconv;
                  }
                  if (iAtm > iCtm) {
                    hessian(3*iAtm + idcoord + 1,3*iCtm + jdcoord + 1) += Fbend*sA[idcoord]*sC[jdcoord]*unitconv;
                    hessian(3*iCtm + jdcoord + 1,3*iAtm + idcoord + 1) += Fbend*sA[idcoord]*sC[jdcoord]*unitconv;
                  }
                  else {
                    hessian(3*iCtm + idcoord + 1,3*iAtm + jdcoord + 1) += Fbend*sC[idcoord]*sA[jdcoord]*unitconv;
                    hessian(3*iAtm + jdcoord + 1,3*iCtm + idcoord + 1) += Fbend*sC[idcoord]*sA[jdcoord]*unitconv;
                  }
                  if (iBtm > iCtm) {
                    hessian(3*iBtm + idcoord + 1,3*iCtm + jdcoord + 1) += Fbend*sB[idcoord]*sC[jdcoord]*unitconv;
                    hessian(3*iCtm + jdcoord + 1,3*iBtm + idcoord + 1) += Fbend*sB[idcoord]*sC[jdcoord]*unitconv;
                  }
                  //this below should never happen!!!!!
                  else {
                    hessian(3*iCtm + idcoord + 1,3*iBtm + jdcoord + 1) += Fbend*sC[idcoord]*sB[jdcoord]*unitconv;
                    hessian(3*iBtm + jdcoord + 1,3*iCtm + idcoord + 1) += Fbend*sC[idcoord]*sB[jdcoord]*unitconv;
                  }
                  if (jdcoord < idcoord) {
                    hessian(3*iBtm + idcoord + 1,3*iBtm + jdcoord + 1) += Fbend*sB[idcoord]*sB[jdcoord]*unitconv;
                    hessian(3*iBtm + jdcoord + 1,3*iBtm + idcoord + 1) += Fbend*sB[idcoord]*sB[jdcoord]*unitconv;
                    hessian(3*iCtm + idcoord + 1,3*iCtm + jdcoord + 1) += Fbend*sC[idcoord]*sC[jdcoord]*unitconv;
                    hessian(3*iCtm + jdcoord + 1,3*iCtm + idcoord + 1) += Fbend*sC[idcoord]*sC[jdcoord]*unitconv;
                    hessian(3*iAtm + idcoord + 1,3*iAtm + jdcoord + 1) += Fbend*sA[idcoord]*sA[jdcoord]*unitconv;
                    hessian(3*iAtm + jdcoord + 1,3*iAtm + idcoord + 1) += Fbend*sA[idcoord]*sA[jdcoord]*unitconv;
                  }
                }
              }
            }
          }
        }
      }
      for (size_t iCtm = 0; iCtm < Natoms; ++iCtm) {
        ACpair = Natoms*iAtm + iCtm;
        if (iCtm == iAtm) {continue;}
        if (iCtm == iBtm) {continue;}
        rAC0 = 1.30*(AtmRadii(atoms[iAtm]) + AtmRadii(atoms[iCtm]))*dist_Angstrom2aum1;
        rABdotrCA = 0.0;
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          xyz(idcoord + 1,1) = geometry(iCtm + 1,idcoord + 1)*dist_Angstrom2aum1;
          oxyz(idcoord + 1,2) = geometry(iCtm + 1,idcoord + 1)*dist_Angstrom2aum1;
          RCA[idcoord] = (geometry(iCtm + 1,idcoord + 1) - geometry(iAtm + 1,idcoord + 1))*dist_Angstrom2aum1;
          RBC[idcoord] = (geometry(iBtm + 1,idcoord + 1) - geometry(iCtm + 1,idcoord + 1))*dist_Angstrom2aum1;
          rABdotrCA += RAB[idcoord]*RCA[idcoord];
        }
        rCA2 = RCA[0]*RCA[0] + RCA[1]*RCA[1] + RCA[2]*RCA[2];
        rBC2 = RBC[0]*RBC[0] + RBC[1]*RBC[1] + RBC[2]*RBC[2];
        if (rBC2 > threshold_dist) {continue;}
        if (rCA2 > threshold_dist) {continue;}
        if (rCA2 > rAC0*rAC0) {continue;}
        for (size_t iDtm = 0; iDtm < Natoms; ++iDtm) {
          if (iDtm == iAtm) {continue;}
          if (iDtm == iBtm) {continue;}
          if (iDtm == iCtm) {continue;}
          rBD0 = 1.30*(AtmRadii(atoms[iBtm]) + AtmRadii(atoms[iDtm]))*dist_Angstrom2aum1;
          rCD2 = 0.0;
          alpha = 0.0;       //RAB.RAD
          RL = 0.0;          //RAC.RAD
          for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
            xyz(idcoord + 1,4) = geometry(iDtm + 1,idcoord + 1)*dist_Angstrom2aum1;
            oxyz(idcoord + 1,3) = geometry(iDtm + 1,idcoord + 1)*dist_Angstrom2aum1;
            RAD[idcoord] = (geometry(iAtm + 1,idcoord + 1) - geometry(iDtm + 1,idcoord + 1))*dist_Angstrom2aum1;
            alpha -= RAB[idcoord]*RAD[idcoord];
            RL -= RCA[idcoord]*RAD[idcoord];
            RBD[idcoord] = (geometry(iBtm + 1,idcoord + 1) - geometry(iDtm + 1,idcoord + 1))*dist_Angstrom2aum1;
            zeroangle = (geometry(iCtm + 1,idcoord + 1) - geometry(iDtm + 1,idcoord + 1))*dist_Angstrom2aum1;
            rCD2 += zeroangle*zeroangle;
          }
          rAD2 = RAD[0]*RAD[0] + RAD[1]*RAD[1] + RAD[2]*RAD[2];
          rBD2 = RBD[0]*RBD[0] + RBD[1]*RBD[1] + RBD[2]*RBD[2];
          if (rAD2 > threshold_dist) {continue;}
          if (rBD2 > threshold_dist) {continue;}
          if (rCD2 > threshold_dist) {continue;}
          if (rBD2 > rBD0*rBD0) {continue;}
          cosf2 = rABdotrCA/sqrt(rAB2*rCA2);
          cosf2fail = (fabs(fabs(cosf2) - 1.0) < 0.1);
          cosf3 = alpha/sqrt(rAB2*rAD2);
          cosf3fail = (fabs(fabs(cosf3) - 1.0) < 0.1);
          cosf4 = RL/sqrt(rCA2*rAD2);
          cosf4fail = (fabs(fabs(cosf4) - 1.0) < 0.1);
          if ((!cosf2fail)&&(!cosf3fail)&&(!cosf4fail)) {             //oop
            aux2 = mol.calcAngles(iAtm + 1,iCtm + 1,iDtm + 1);
            if ((fabs(aux2) < threshold)||(fabs(pi - aux2) < threshold)) {continue;}
            aux = sin(aux2);
            cross[0] = -(RCA[1]*RAD[2] - RCA[2]*RAD[1])/aux;
            cross[1] = -(RCA[2]*RAD[0] - RCA[0]*RAD[2])/aux;
            cross[2] = -(RCA[0]*RAD[1] - RCA[1]*RAD[0])/aux;
            aux = 0.0;
            for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
              aux += RAB[idcoord]*cross[idcoord];
            }
            aux /= sqrt(rAB2*rCA2*rAD2);
            aux2 = 1.0 - fabs(aux);
            Foop = 0.0585*aux2*aux2*aux2*aux2;
            OutOfPlane(oxyz,alpha,Bmat,1.0e-13);
            if (fabs(alpha) > taumax) {continue;}
            for (size_t idx = 0; idx < 3; ++idx) {
              sC[idx] = Bmat(idx + 1,2);
              sA[idx] = Bmat(idx + 1,4);
              sB[idx] = Bmat(idx + 1,1);
              sD[idx] = Bmat(idx + 1,3);
            }
            for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
              for (size_t jdcoord = 0; jdcoord < 3; ++jdcoord) {
                hessian(3*iAtm + idcoord + 1,3*iBtm + jdcoord + 1) += Foop*sA[idcoord]*sB[jdcoord];
                hessian(3*iAtm + idcoord + 1,3*iCtm + jdcoord + 1) += Foop*sA[idcoord]*sC[jdcoord];
                hessian(3*iAtm + idcoord + 1,3*iDtm + jdcoord + 1) += Foop*sA[idcoord]*sD[jdcoord];
                hessian(3*iBtm + jdcoord + 1,3*iAtm + idcoord + 1) += Foop*sA[idcoord]*sB[jdcoord];
                hessian(3*iCtm + jdcoord + 1,3*iAtm + idcoord + 1) += Foop*sA[idcoord]*sC[jdcoord];
                hessian(3*iDtm + jdcoord + 1,3*iAtm + idcoord + 1) += Foop*sA[idcoord]*sD[jdcoord];
                hessian(3*iBtm + idcoord + 1,3*iCtm + jdcoord + 1) += Foop*sB[idcoord]*sC[jdcoord];
                hessian(3*iBtm + idcoord + 1,3*iDtm + jdcoord + 1) += Foop*sB[idcoord]*sD[jdcoord];
                hessian(3*iCtm + idcoord + 1,3*iDtm + jdcoord + 1) += Foop*sC[idcoord]*sD[jdcoord];
                hessian(3*iCtm + jdcoord + 1,3*iBtm + idcoord + 1) += Foop*sB[idcoord]*sC[jdcoord];
                hessian(3*iDtm + jdcoord + 1,3*iBtm + idcoord + 1) += Foop*sB[idcoord]*sD[jdcoord];
                hessian(3*iDtm + jdcoord + 1,3*iCtm + idcoord + 1) += Foop*sC[idcoord]*sD[jdcoord];
              }
            }
            for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
              hessian(3*iAtm + idcoord + 1,3*iAtm + idcoord + 1) += Foop*sA[idcoord]*sA[idcoord];
              hessian(3*iBtm + idcoord + 1,3*iBtm + idcoord + 1) += Foop*sB[idcoord]*sB[idcoord];
              hessian(3*iCtm + idcoord + 1,3*iCtm + idcoord + 1) += Foop*sC[idcoord]*sC[idcoord];
              hessian(3*iDtm + idcoord + 1,3*iDtm + idcoord + 1) += Foop*sD[idcoord]*sD[idcoord];
              for (size_t jdcoord = 0; jdcoord < idcoord; ++jdcoord) {
                hessian(3*iAtm + idcoord + 1,3*iAtm + jdcoord + 1) += Foop*sA[idcoord]*sA[jdcoord];
                hessian(3*iBtm + idcoord + 1,3*iBtm + jdcoord + 1) += Foop*sB[idcoord]*sB[jdcoord];
                hessian(3*iCtm + idcoord + 1,3*iCtm + jdcoord + 1) += Foop*sC[idcoord]*sC[jdcoord];
                hessian(3*iDtm + idcoord + 1,3*iDtm + jdcoord + 1) += Foop*sD[idcoord]*sD[jdcoord];
                hessian(3*iAtm + jdcoord + 1,3*iAtm + idcoord + 1) += Foop*sA[idcoord]*sA[jdcoord];
                hessian(3*iBtm + jdcoord + 1,3*iBtm + idcoord + 1) += Foop*sB[idcoord]*sB[jdcoord];
                hessian(3*iCtm + jdcoord + 1,3*iCtm + idcoord + 1) += Foop*sC[idcoord]*sC[jdcoord];
                hessian(3*iDtm + jdcoord + 1,3*iDtm + idcoord + 1) += Foop*sD[idcoord]*sD[jdcoord];
              }
            }
          }
          //torsion
          BDpair = Natoms*iBtm + iDtm;
          if (ACpair <= BDpair) {continue;}
          cosf2 = -(RCA[0]*RAB[0] + RCA[1]*RAB[1] + RCA[2]*RAB[2])/sqrt(rCA2*rAB2);
          cosf2fail = (fabs(cosf2) > cosphimax);
          cosf3 = -(RBD[0]*RAB[0] + RBD[1]*RAB[1] + RBD[2]*RAB[2])/sqrt(rBD2*rAB2);
          cosf3fail = (fabs(cosf3) > cosphimax);
          if ((!cosf2fail)&&(!cosf3fail)) {
            alpha = (AtmRadii(atoms[iBtm]) + AtmRadii(atoms[iAtm]))*dist_Angstrom2aum1;                //I chose here to use Psi4's covalent radii
            Ftorsion = 0.00299 - 0.091*(rAB - alpha);
            if (Ftorsion < threshold) {continue;}
            Torsion(xyz,alpha,Bmat,dBmat,false);
            for (size_t idx = 0; idx < 3; ++idx) {
              sC[idx] = Bmat(idx + 1,1);
              sA[idx] = Bmat(idx + 1,2);
              sB[idx] = Bmat(idx + 1,3);
              sD[idx] = Bmat(idx + 1,4);
            }
            for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
              hessian(3*iCtm + idcoord + 1,3*iCtm + idcoord + 1) += Ftorsion*sC[idcoord]*sC[idcoord]*unitconv;
              hessian(3*iAtm + idcoord + 1,3*iAtm + idcoord + 1) += Ftorsion*sA[idcoord]*sA[idcoord]*unitconv;
              hessian(3*iBtm + idcoord + 1,3*iBtm + idcoord + 1) += Ftorsion*sB[idcoord]*sB[idcoord]*unitconv;
              hessian(3*iDtm + idcoord + 1,3*iDtm + idcoord + 1) += Ftorsion*sD[idcoord]*sD[idcoord]*unitconv;
              for (size_t jdcoord = 0; jdcoord < 3; ++jdcoord) {
                hessian(3*iCtm + idcoord + 1,3*iAtm + jdcoord + 1) += Ftorsion*sC[idcoord]*sA[jdcoord]*unitconv;
                hessian(3*iAtm + jdcoord + 1,3*iCtm + idcoord + 1) += Ftorsion*sC[idcoord]*sA[jdcoord]*unitconv;
                hessian(3*iCtm + idcoord + 1,3*iBtm + jdcoord + 1) += Ftorsion*sC[idcoord]*sB[jdcoord]*unitconv;
                hessian(3*iBtm + jdcoord + 1,3*iCtm + idcoord + 1) += Ftorsion*sC[idcoord]*sB[jdcoord]*unitconv;
                hessian(3*iCtm + idcoord + 1,3*iDtm + jdcoord + 1) += Ftorsion*sC[idcoord]*sD[jdcoord]*unitconv;
                hessian(3*iDtm + jdcoord + 1,3*iCtm + idcoord + 1) += Ftorsion*sC[idcoord]*sD[jdcoord]*unitconv;
                hessian(3*iAtm + idcoord + 1,3*iBtm + jdcoord + 1) += Ftorsion*sA[idcoord]*sB[jdcoord]*unitconv;
                hessian(3*iBtm + jdcoord + 1,3*iAtm + idcoord + 1) += Ftorsion*sA[idcoord]*sB[jdcoord]*unitconv;
                hessian(3*iAtm + idcoord + 1,3*iDtm + jdcoord + 1) += Ftorsion*sA[idcoord]*sD[jdcoord]*unitconv;
                hessian(3*iDtm + jdcoord + 1,3*iAtm + idcoord + 1) += Ftorsion*sA[idcoord]*sD[jdcoord]*unitconv;
                hessian(3*iBtm + idcoord + 1,3*iDtm + jdcoord + 1) += Ftorsion*sB[idcoord]*sD[jdcoord]*unitconv;
                hessian(3*iDtm + jdcoord + 1,3*iBtm + idcoord + 1) += Ftorsion*sB[idcoord]*sD[jdcoord]*unitconv;
                if (jdcoord < idcoord) {
                  hessian(3*iCtm + idcoord + 1,3*iCtm + jdcoord + 1) += Ftorsion*sC[idcoord]*sC[jdcoord]*unitconv;
                  hessian(3*iCtm + jdcoord + 1,3*iCtm + idcoord + 1) += Ftorsion*sC[idcoord]*sC[jdcoord]*unitconv;
                  hessian(3*iAtm + idcoord + 1,3*iAtm + jdcoord + 1) += Ftorsion*sA[idcoord]*sA[jdcoord]*unitconv;
                  hessian(3*iAtm + jdcoord + 1,3*iAtm + idcoord + 1) += Ftorsion*sA[idcoord]*sA[jdcoord]*unitconv;
                  hessian(3*iBtm + idcoord + 1,3*iBtm + jdcoord + 1) += Ftorsion*sB[idcoord]*sB[jdcoord]*unitconv;
                  hessian(3*iBtm + jdcoord + 1,3*iBtm + idcoord + 1) += Ftorsion*sB[idcoord]*sB[jdcoord]*unitconv;
                  hessian(3*iDtm + idcoord + 1,3*iDtm + jdcoord + 1) += Ftorsion*sD[idcoord]*sD[jdcoord]*unitconv;
                  hessian(3*iDtm + jdcoord + 1,3*iDtm + idcoord + 1) += Ftorsion*sD[idcoord]*sD[jdcoord]*unitconv;
                }
              }
            }
          }
        }
      }
    }
  }
  //since there might be numerical inaccuracies in the the symmetry, symmetrize
  for (size_t idr = 0; idr < 3*Natoms; ++idr) {
    for (size_t idc = 0; idc < idr; ++idc) {
      Fstretch = 0.5*(hessian(idr + 1,idc + 1) + hessian(idc + 1,idr + 1));
      hessian(idr + 1,idc + 1) = Fstretch;
      hessian(idc + 1,idr + 1) = Fstretch;
    }
  }
}
//function to project translations and rotations out of the Hessian, making their contributions exactly zero
double Tens(int xx, int yy, int zz) {
  //totally asymmetric Cartesian tensor
  double tensor = 0.0;
  if ((xx == 1)&&(yy == 2)&&(zz == 3)) {tensor = 1.0;}
  else if ((xx == 2)&&(yy == 3)&&(zz == 1)) {tensor = 1.0;}
  else if ((xx == 3)&&(yy == 1)&&(zz == 2)) {tensor = 1.0;}
  else if ((xx == 1)&&(yy == 3)&&(zz == 2)) {tensor = -1.0;}
  else if ((xx == 2)&&(yy == 1)&&(zz == 3)) {tensor = -1.0;}
  else if ((xx == 3)&&(yy == 2)&&(zz == 1)) {tensor = -1.0;}
  return tensor;
}
void ProjectForceConstants(matrixE & ForceMat, matrixE & Geometry) {
  //function that removes external degrees of freedom (zeroes them out) from the Hessian matrix
  //based now on the gtrprojm subroutine in xTB
  int Natoms = Geometry.rows();
  int dimensions = Geometry.cols();
  int Nvar = Natoms*dimensions;
  matrixE projection(Nvar,6);
  std::vector<double> PCMass(Natoms,0.0);
  if ((Nvar != ForceMat.rows())&&(Nvar != ForceMat.cols())) {throw("ERROR: Hessian.hpp: ProjectForceConstants(): force matrix and geometry inconsistent");}
  //get the point center of mass
  for (size_t idatm = 0; idatm < Natoms; ++idatm) {
    PCMass[0] += Geometry(idatm + 1,1);
    PCMass[1] += Geometry(idatm + 1,2);
    PCMass[2] += Geometry(idatm + 1,3);
  }
  PCMass[0] /= double(Natoms);
  PCMass[1] /= double(Natoms);
  PCMass[2] /= double(Natoms);
  //get basic projector
  for (size_t idatm = 0; idatm < Natoms; ++idatm) {
    for (size_t idc = 1; idc < 4; ++idc) {
      //correct geometry for point CM
      Geometry(idatm + 1,idc) -= PCMass[idc - 1];
      //translation part of the projector
      projection(3*idatm + idc,idc) = 1.0;
    }
    //rotation part of the projector
    projection(3*idatm + 2,4) = -Geometry(idatm + 1,3);
    projection(3*idatm + 3,4) = Geometry(idatm + 1,2);
    projection(3*idatm + 1,5) = Geometry(idatm + 1,3);
    projection(3*idatm + 3,5) = -Geometry(idatm + 1,1);
    projection(3*idatm + 1,6) = -Geometry(idatm + 1,2);
    projection(3*idatm + 2,6) = Geometry(idatm + 1,1);
  }
  //if different projects are to be used, then implement them here
  //possibilities include more than just translations and rotations or to have only rotations
  GramSchmidtOrthnorm(projection);
  matrixE PPT = projection*projection.trans();
  //since we multiply on both sides, it does not matter if we have (1 - PP*) or (PP* - 1)
  for (size_t iddiag = 1; iddiag < Nvar + 1; ++iddiag) {
    PPT(iddiag,iddiag) -= 1.0;
  }
  matrixE Aproj = ForceMat*PPT;
  ForceMat = PPT*Aproj;
}

#endif //_APPROXIMATE_HESSIAN_
