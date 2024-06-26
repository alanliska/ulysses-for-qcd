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

#ifndef _Geometry_
#define _Geometry_
#include "ConstantsPackage.hpp"
#include "atoms/AtomPackage.hpp"
#include "math/MatrixPackage.hpp"
#include "math/VectorPackage.hpp"
#include "other/Utilities.hpp"

//descrition:
//"geometrical" functions of molecules

std::vector<double> CenterOfMass(matrixE & geometry, std::vector<size_t> & atoms) {
  //function that calculates a center of mass
  size_t natoms = atoms.size();
  if (natoms != geometry.rows()) {throw("ERROR: Geometry.hpp: CenterOfMass(): geometry matrix and atom vector with inconsistent dimensions");}
  std::vector<double> cm(3,0.0);
  double mass_sum = 0.0;
  double mass;
  for (size_t idx = 0; idx < natoms; ++idx) {
    mass = Weight(atoms[idx]);
    cm[0] += geometry(idx + 1,1)*mass;
    cm[1] += geometry(idx + 1,2)*mass;
    cm[2] += geometry(idx + 1,3)*mass;
    mass_sum += mass;
  }
  cm[0] /= mass_sum;
  cm[1] /= mass_sum;
  cm[2] /= mass_sum;
  return cm;
}
matrixE InertiaMatrix(matrixE & geometry, std::vector<size_t> & atoms, bool inmeter = true) {
  //function calculating molecular inertia matrix
  matrixE _inertia(1,1);
  size_t natoms = atoms.size();
  double unitconvert = 1.0;
  if (inmeter) {unitconvert = 1.0e-10;}
  if (natoms != geometry.rows()) {throw("ERROR: Geometry.hpp: CenterOfMass(): geometry matrix and atom vector with inconsistent dimensions");}
  if (natoms == 2) {
    double m1 = Weight(atoms[0])/NA;                              //atomic masses in kg; convert to kg/mol
    double m2 = Weight(atoms[1])/NA;                              //atomic masses in kg; convert to kg/mol
    double rab;
    double rab2 = 0.0;
    for (size_t idx = 1; idx < 4; ++idx) {
      rab = geometry(2,idx) - geometry(1,idx);
      rab2 += rab*rab;
    }
    _inertia(1,1) = unitconvert*unitconvert*rab2*m1*m2/(m1 + m2);
  }
  else {
    _inertia.resize(3,3);
    _inertia.zero();
    double mi;
    double xi;
    double yi;
    double zi;
    for (size_t idx = 0; idx < natoms; ++idx) {
      mi = Weight(atoms[idx])/NA;                                 //atomic masses in kg; convert to kg/mol
      xi = geometry(idx + 1,1)*unitconvert;
      yi = geometry(idx + 1,2)*unitconvert;
      zi = geometry(idx + 1,3)*unitconvert;
      _inertia(1,1) += mi*(yi*yi + zi*zi);
      _inertia(2,1) -= mi*xi*yi;
      _inertia(3,1) -= mi*xi*zi;
      _inertia(2,2) += mi*(xi*xi + zi*zi);
      _inertia(3,2) -= mi*yi*zi;
      _inertia(3,3) += mi*(xi*xi + yi*yi);
    }
    _inertia(1,2) = _inertia(2,1);
    _inertia(1,3) = _inertia(3,1);
    _inertia(2,3) = _inertia(3,2);
  }
  return _inertia;
}
double Distance(size_t posA, size_t posB, const matrixE & geom) {
  //function that calculates a distance
  double aux;
  double RAB2 = 0.0;
  for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
    aux = geom(posA,idcoord) - geom(posB,idcoord);
    RAB2 += aux*aux;
  }
  return sqrt(RAB2);
}
void Distance2(size_t posA, size_t posB, const matrixE & geom, double & aux, double & RAB2) {
  //function that calculates a squared distance
  RAB2 = 0.0;
  for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
    aux = geom(posA,idcoord) - geom(posB,idcoord);
    RAB2 += aux*aux;
  }
}
double Angle(size_t posA, size_t posB, size_t posC, const matrixE & geom, double thresholdzero = 1.0e-20) {
  //function calculating the angle between three atoms
  double aux;
  double RAB2 = 0.0;
  double RAC2 = 0.0;
  double RBC2 = 0.0;
  double angle = 0.0;
  for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
    aux = geom(posA,idcoord) - geom(posB,idcoord);
    RAB2 += aux*aux;
    aux = geom(posA,idcoord) - geom(posC,idcoord);
    RAC2 += aux*aux;
    aux = geom(posC,idcoord) - geom(posB,idcoord);
    RBC2 += aux*aux;
  }
  angle = sqrt(RAB2*RBC2);
  if (angle > thresholdzero) {
    aux = 0.5*(RAB2 + RBC2 - RAC2)/angle;
    angle = fmin(1.0,aux);
    aux = fmax(-1.0,angle);
    angle = acos(aux);
  }
  return angle;
}
double AngleDer(size_t posA, size_t posB, size_t posC, const matrixE & geom, matrixE & der, std::vector<double> & rAB, std::vector<double> & rCB, std::vector<double> & rAC, double thresholdzero = 1.0e-20) {
  //function calculating the angle between three atoms and also the respective derivatives
  double RAB2 = 0.0;
  double RAC2 = 0.0;
  double RBC2 = 0.0;
  double angle = 0.0;
  double sint;
  double sintm1;
  double RABC;
  double RABC2;
  double RABm2;
  double RBCm2;
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    rAB[idcoord] = geom(posA,idcoord + 1) - geom(posB,idcoord + 1);
    RAB2 += rAB[idcoord]*rAB[idcoord];
    rAC[idcoord] = geom(posA,idcoord + 1) - geom(posC,idcoord + 1);
    RAC2 += rAC[idcoord]*rAC[idcoord];
    rCB[idcoord] = geom(posC,idcoord + 1) - geom(posB,idcoord + 1);
    RBC2 += rCB[idcoord]*rCB[idcoord];
  }
  //get the angle
  RABC = sqrt(RAB2*RBC2);
  angle = RABC;
  RABC2 = 0.5*(RAB2 + RBC2 - RAC2);
  if (angle > thresholdzero) {
    sintm1 = RABC2/angle;
    angle = fmin(1.0,sintm1);
    sintm1 = fmax(-1.0,angle);
    angle = acos(sintm1);
  }
  //now the derivatives
  sint = sin(angle);
  if (angle < 0.001) {sint = 1.0;}
  sintm1 = 1.0/(sint*RABC);
  if (RAB2 > 1.0e-5) {RABm2 = 1.0/RAB2;}
  else {RABm2 = RAB2;}
  if (RBC2 > 1.0e-5) {RBCm2 = 1.0/RBC2;}
  else {RBCm2 = RBC2;}
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    der(1,idcoord + 1) = -sintm1*(rCB[idcoord] - rAB[idcoord]*RABC2*RABm2);         //dt/dA
    der(3,idcoord + 1) = -sintm1*(rAB[idcoord] - rCB[idcoord]*RABC2*RBCm2);         //dt/dC
    der(2,idcoord + 1) = -der(1,idcoord + 1) - der(3,idcoord + 1);                  //dt/dB
  }
  return angle;
}
double DihedralAngle(double & a1, double & a2, double & b1, double & b2, double thresholdzero = 1.0e-6, double thresholdzero2 = 4.0e-5) {
  //auxiliary function for the calculation of torsion angles, returns angle in radians
  double normA;
  double normB;
  double sinth;
  double costh;
  double rcos = 0.0;
  if ((fabs(a1) > thresholdzero)||(fabs(a2) > thresholdzero)) {
    if ((fabs(b1) > thresholdzero)||(fabs(b2) > thresholdzero)) {
      normA = 1.0/sqrt(a1*a1 + a2*a2);
      normB = 1.0/sqrt(b1*b1 + b2*b2);
      a1 *= normA;
      a2 *= normA;
      b1 *= normB;
      b2 *= normB;
      sinth = a1*b2 - a2*b1;
      costh = a1*b1 + a2*b2;
      normA = fmin(1.0,costh);
      costh = fmax(-1.0,normA);
      rcos = acos(costh);
      if (fabs(rcos) > thresholdzero2) {
        if (sinth > 0.0) {rcos = 2.0*pi - rcos;}
        rcos *= -1.0;
      }
    }
  }
  return rcos;
}
double TorsionDer(size_t posA, size_t posB, size_t posC, size_t posD, const matrixE & geom, matrixE & der, std::vector<double> & rAB, std::vector<double> & rBC, std::vector<double> & rCD, matrixE & xprod, double thresholdzero = 1.0e-5) {
  //function calculating a torsion angle and its derivative, taken from openbabel, returns angle in radians
  //get vectors between atoms and the respective distances
  double RBA = 0.0;
  double RCB = 0.0;
  double RDC = 0.0;
  double RS2B = 0.0;
  double RS2C = 0.0;
  double auxA = 0.0;
  double auxB = 0.0;
  double auxC = 0.0;
  double cosABC = 0.0;
  double cosBCD = 0.0;
  double sinABC = 0.0;
  double sinBCD = 0.0;
  double aABC = 0.0;
  double aBCD = 0.0;
  double d1 = 0.0;
  double d2 = 0.0;
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    rAB[idcoord] = (geom(posB,idcoord + 1) - geom(posA,idcoord + 1));
    rBC[idcoord] = (geom(posC,idcoord + 1) - geom(posB,idcoord + 1));
    rCD[idcoord] = (geom(posD,idcoord + 1) - geom(posC,idcoord + 1));
    auxA += rAB[idcoord]*rAB[idcoord];
    auxB += rBC[idcoord]*rBC[idcoord];
    auxC += rCD[idcoord]*rCD[idcoord];
  }
  RBA = sqrt(auxA) + 1.0e-7;      //avoid division by zero
  RCB = sqrt(auxB) + 1.0e-7;      //avoid division by zero
  RDC = sqrt(auxC) + 1.0e-7;      //avoid division by zero
  auxA = 1.0/RBA;
  auxB = 1.0/RCB;
  auxC = 1.0/RDC;
  //bond vector normalization
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    rAB[idcoord] *= auxA;
    rBC[idcoord] *= auxB;
    rCD[idcoord] *= auxC;
    cosABC += rAB[idcoord]*rBC[idcoord];
    cosBCD += rBC[idcoord]*rCD[idcoord];
  }
  //get angle ABC
  if (cosABC > 1.0) {
    cosABC = 1.0;
    aABC = 0.0;
  }
  else if (cosABC < -1.0) {
    cosABC = -1.0;
    aABC = pi;
  }
  else {aABC = acos(cosABC);}
  //get angle BCD
  if (cosBCD > 1.0) {
    cosBCD = 1.0;
    aBCD = 0.0;
  }
  else if (cosBCD < -1.0) {
    cosBCD = -1.0;
    aBCD = pi;
  }
  else {aBCD = acos(cosBCD);}
  //some auxiliary stuff
  sinABC = sin(aABC);
  sinBCD = sin(aBCD);
  RS2B = 1.0/(RBA*sinABC*sinABC);
  RS2C = 1.0/(RDC*sinBCD*sinBCD);
  auxC = RDC*auxB;
  auxB *= RBA;
  auxB *= -cosABC;
  auxC *= -cosBCD;
  CrossProd(rAB,rBC,xprod,1,0);             //from vector package, the ones down in this file are to remove with time
  CrossProd(rBC,rCD,xprod,2,0);             //from vector package, the ones down in this file are to remove with time
  CrossProdM(xprod,1,2,3);
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    d1 += xprod(3,idcoord + 1)*rBC[idcoord];
    d2 += xprod(1,idcoord + 1)*xprod(2,idcoord + 1);
  }
  //get the derivatives
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    der(1,idcoord + 1) = -RS2B*xprod(1,idcoord + 1);
    der(3,idcoord + 1) = 0.0;
    der(4,idcoord + 1) = RS2C*xprod(2,idcoord + 1);
    xprod(1,idcoord + 1) = der(1,idcoord + 1)*(auxB - 1.0);
    xprod(2,idcoord + 1) = der(4,idcoord + 1)*auxC;
    der(2,idcoord + 1) = xprod(1,idcoord + 1) - xprod(2,idcoord + 1);
    xprod(1,idcoord + 1) = der(2,idcoord + 1) + der(1,idcoord + 1);
    xprod(2,idcoord + 1) = xprod(1,idcoord + 1) + der(4,idcoord + 1);
    der(3,idcoord + 1) = -xprod(2,idcoord + 1);
  }
  //return the angle
  return atan2(d1,d2);
}
double Torsion(size_t posA, size_t posB, size_t posC, size_t posD, const matrixE & geom) {
  //function calculating a torsion angle
  double AC[3];
  double BC[3];
  double DC[3];
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    AC[idcoord] = (geom(posA,idcoord + 1) - geom(posC,idcoord + 1))*dist_Angstrom2aum1;
    BC[idcoord] = (geom(posB,idcoord + 1) - geom(posC,idcoord + 1))*dist_Angstrom2aum1;
    DC[idcoord] = (geom(posD,idcoord + 1) - geom(posC,idcoord + 1))*dist_Angstrom2aum1;
  }
  double dist = sqrt(BC[0]*BC[0] + BC[1]*BC[1] + BC[2]*BC[2]);
  double cosA = 1.0;
  if (dist > 0.0) {cosA = BC[2]/dist;}
  double sinA = fmin(1.0,cosA);
  cosA = fmax(-1.0,sinA);
  sinA = 1.0 - cosA*cosA;
  double yxdist = dist*sqrt(sinA);
  double cosph;
  double sinph;
  double xi2;
  double xl2;
  double yi2;
  double yj2;
  double yl2;
  double costh;
  double sinth;
  if ((sinA > 0.0)&&(yxdist > 1.0e-6)) {
    cosph = BC[1]/yxdist;
    sinph = BC[0]/yxdist;
    xi2 = AC[0]*cosph - AC[1]*sinph;
    xl2 = DC[0]*cosph - DC[1]*sinph;
    yi2 = AC[0]*sinph + AC[1]*cosph;
    yj2 = BC[0]*sinph + BC[1]*cosph;
    yl2 = DC[0]*sinph + DC[1]*cosph;
    costh = cosA;
    sinth = yj2/dist;
  }
  else {
    xi2 = AC[0];
    xl2 = DC[0];
    yi2 = AC[1];
    yl2 = DC[1];
    costh = cosA;
    sinth = 0.0;
  }
  double yi3 = yi2*costh - AC[2]*sinth;
  double yl3 = yl2*costh - DC[2]*sinth;
  sinth = DihedralAngle(xl2,yl3,xi2,yi3);                    //the angle
  if (sinth < 0.0) {sinth += 2.0*pi;}
  if (sinth >= 2.0*pi) {sinth = 0.0;}
  return sinth;
}
double InversionDer(size_t posA, size_t posB, size_t posC, size_t posD, const matrixE & geom, matrixE & der, std::vector<double> & rAB, std::vector<double> & rCA, std::vector<double> & rDA) {
  //function calculating the inversion angle and its derivatives
  double crossABC[3];
  double crossADAC[3];
  double crossADBA[3];
  double RAD = 0.0;
  double RAD2 = 0.0;
  double nABC = 0.0;
  double nABC2 = 0.0;
  double wABCD;
  double ADABC = 0.0;
  double normcrossABD = 0.0;
  double dotABCADAC = 0.0;
  double dotABCADBA = 0.0;
  double factorB;
  double factorC;
  double factorD;
  double xval;
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    rAB[idcoord] = geom(posA,idcoord + 1) - geom(posB,idcoord + 1);
    rCA[idcoord] = geom(posC,idcoord + 1) - geom(posA,idcoord + 1);
    rDA[idcoord] = geom(posD,idcoord + 1) - geom(posA,idcoord + 1);
  }
  crossABC[0] = rAB[1]*rCA[2] - rAB[2]*rCA[1];
  crossABC[1] = rAB[2]*rCA[0] - rAB[0]*rCA[2];
  crossABC[2] = rAB[0]*rCA[1] - rAB[1]*rCA[0];
  crossADAC[0] = rDA[1]*rCA[2] - rDA[2]*rCA[1];
  crossADAC[1] = rDA[2]*rCA[0] - rDA[0]*rCA[2];
  crossADAC[2] = rDA[0]*rCA[1] - rDA[1]*rCA[0];
  crossADBA[0] = rDA[1]*rAB[2] - rDA[2]*rAB[1];
  crossADBA[1] = rDA[2]*rAB[0] - rDA[0]*rAB[2];
  crossADBA[2] = rDA[0]*rAB[1] - rDA[1]*rAB[0];
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    RAD2 += rDA[idcoord]*rDA[idcoord];
    nABC2 += crossABC[idcoord]*crossABC[idcoord];
    ADABC += rDA[idcoord]*crossABC[idcoord];
    dotABCADAC += crossABC[idcoord]*crossADAC[idcoord];
    dotABCADBA += crossABC[idcoord]*crossADBA[idcoord];
  }
  RAD = 1.0/sqrt(RAD2);
  nABC = 1.0/sqrt(nABC2);
  xval = ADABC*RAD*nABC;
  wABCD = 1.0/sqrt(1.0 - xval*xval);
  //now the derivatives
  factorB = wABCD*dotABCADAC*RAD*nABC*nABC*nABC;
  factorC = wABCD*dotABCADBA*RAD*nABC*nABC*nABC;
  factorD = wABCD*RAD*RAD*RAD*nABC;
  wABCD = acos(xval);
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    der(2,idcoord + 1) = factorB*crossABC[idcoord];                                            //terms for atom B
    der(3,idcoord + 1) = factorC*crossABC[idcoord];                                            //terms for atom C
    der(4,idcoord + 1) = factorD*(RAD2*crossABC[idcoord] - ADABC*rDA[idcoord]);                //terms for atom D
    der(1,idcoord + 1) = -der(2,idcoord + 1) - der(3,idcoord + 1) - der(4,idcoord + 1);        //terms for atom A
  }
  return 0.5*pi - wABCD;
}
void CrossProduct(size_t posA, size_t posB, size_t posC, size_t posD, std::vector<double> & crossprod, matrixE & geometry) {
  //function that calculates the cross product between differences of atom coordinates in a geometry matrix
  crossprod.resize(3);
  crossprod[0] = (geometry(posA,2) - geometry(posB,2))*(geometry(posC,3) - geometry(posD,3)) - (geometry(posA,3) - geometry(posB,3))*(geometry(posC,2) - geometry(posD,2));
  crossprod[1] = (geometry(posA,3) - geometry(posB,3))*(geometry(posC,1) - geometry(posD,1)) - (geometry(posA,1) - geometry(posB,1))*(geometry(posC,3) - geometry(posD,3));
  crossprod[2] = (geometry(posA,1) - geometry(posB,1))*(geometry(posC,2) - geometry(posD,2)) - (geometry(posA,2) - geometry(posB,2))*(geometry(posC,1) - geometry(posD,1));
}
void CrossProduct(size_t posA, size_t posB, const std::vector<double> & diff, std::vector<double> & crossprod, matrixE & geometry, double type) {
  //function that calculates the cross product between a difference of atom coordinates in a geometry matrix and a vector
  crossprod.resize(3);
  crossprod[0] = (geometry(posB,2) - geometry(posA,2))*diff[2] - (geometry(posB,3) - geometry(posA,3))*diff[1];
  crossprod[1] = (geometry(posB,3) - geometry(posA,3))*diff[0] - (geometry(posB,1) - geometry(posA,1))*diff[2];
  crossprod[2] = (geometry(posB,1) - geometry(posA,1))*diff[1] - (geometry(posB,2) - geometry(posA,2))*diff[0];
}
double DotProduct(size_t posA, size_t posB, const std::vector<double> & diff, matrixE & geometry) {
  //function that calculates the dot product between a difference of atom coordinates in a geometry matrix and a vector
  double dotprd = 0.0;
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    dotprd += (geometry(posB,idcoord + 1) - geometry(posA,idcoord + 1))*diff[idcoord];
  }
  return dotprd;
}
double diffAngle(size_t posA, size_t posB, size_t posC, matrixE & geometry) {
  //function returning the angle between difference vectors
  double rA2 = 0.0;
  double rB2 = 0.0;
  double AB = 0.0;
  double auxA;
  double auxB;
  for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
    auxA = geometry(posC,idcoord) - geometry(posA,idcoord);
    auxB = geometry(posB,idcoord) - geometry(posC,idcoord);
    rA2 += auxA*auxA;
    rB2 += auxB*auxB;
    AB += auxA*auxB;
  }
  return acos(AB/sqrt(rA2*rB2));
}
void Internal2Cartesian(size_t index, double ge, size_t posA, size_t posB, size_t posC, size_t posD, matrixE & geometry, std::vector<double> & gradient) {
  //function that transforms the gradient from internal to Cartesian coordinates (in Bohr)
  double AB[3];
  double CB[3];
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    AB[idcoord] = geometry(posA,idcoord + 1) - geometry(posB,idcoord + 1);
  }
  double rAB = 1.0/sqrt(AB[0]*AB[0] + AB[1]*AB[1] + AB[2]*AB[2]);
  double rBC = 0.0;
  if (posC != 0) {
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      CB[idcoord] = geometry(posC,idcoord + 1) - geometry(posB,idcoord + 1);
    }
    rBC = 1.0/sqrt(CB[0]*CB[0] + CB[1]*CB[1] + CB[2]*CB[2]);
  }
  if (index == 2) {                       //distances
    if (posA != posB) {
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        gradient[3*(posA - 1) + idcoord] += ge*rAB*AB[idcoord]*dist_Angstrom2aum1;
        gradient[3*(posB - 1) + idcoord] -= ge*rAB*AB[idcoord]*dist_Angstrom2aum1;
      }
    }
  }
  else if (index == 3) {                  //angles
    if ((posA != posB)&&(posA != posC)&&(posB != posC)) {
      double ABdotCB = AB[0]*CB[0] + AB[1]*CB[1] + AB[2]*CB[2];
      double tmp1;
      double tmp2;
      double aux = 1.0/sqrt(1.0 - ABdotCB*ABdotCB*rAB*rAB*rBC*rBC);
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        tmp1 = -ge*aux*(CB[idcoord]*rAB*rBC - AB[idcoord]*ABdotCB*rAB*rAB*rAB*rBC);
        tmp2 = -ge*aux*(AB[idcoord]*rAB*rBC - CB[idcoord]*ABdotCB*rAB*rBC*rBC*rBC);
        gradient[3*(posA - 1) + idcoord] += tmp1;
        gradient[3*(posB - 1) + idcoord] -= tmp1 + tmp2;
        gradient[3*(posC - 1) + idcoord] += tmp2;
      }
    }
  }
  else if (index == 4) {                  //torsions
    bool Arepeat = (posA == posB)||(posA == posC)||(posA == posD);
    bool Brepeat = (posB == posC)||(posB == posD);
    bool Crepeat = (posC == posD);
    if ((!Arepeat)&&(!Brepeat)&&(!Crepeat)) {
      double tmp1 = diffAngle(posA,posC,posB,geometry);
      double tmp2 = diffAngle(posB,posD,posC,geometry);
      if ((tmp1 != 0.0)&&(tmp1 != pi)) {
        if ((tmp2 != 0.0)&&(tmp2 != pi)) {
          std::vector<double> cd1d2;
          std::vector<double> cd2d3;
          std::vector<double> CBcd2d3;
          std::vector<double> CBcd1d2;
          std::vector<double> BADC;
          std::vector<double> BABD;
          std::vector<double> BAcd2d3;
          std::vector<double> ACcd2d3;
          std::vector<double> DCcd1d2;
          std::vector<double> BDcd1d2;
          double grad[4];
          CrossProduct(posB,posA,posC,posB,cd1d2,geometry);
          CrossProduct(posC,posB,posD,posC,cd2d3,geometry);
          CrossProduct(posB,posC,cd2d3,CBcd2d3,geometry,1.0);
          CrossProduct(posB,posC,cd1d2,CBcd1d2,geometry,1.0);
          CrossProduct(posB,posA,posD,posC,BADC,geometry);
          CrossProduct(posB,posA,posB,posD,BABD,geometry);
          CrossProduct(posA,posB,cd2d3,BAcd2d3,geometry,1.0);
          CrossProduct(posC,posA,cd2d3,ACcd2d3,geometry,1.0);
          CrossProduct(posD,posC,cd1d2,DCcd1d2,geometry,-1.0);
          CrossProduct(posB,posD,cd1d2,BDcd1d2,geometry,-1.0);
          tmp1 = 0.0;
          tmp2 = 0.0;
          for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
            tmp1 += (geometry(posB,idcoord + 1) - geometry(posA,idcoord + 1))*cd2d3[idcoord];
            tmp2 += cd1d2[idcoord]*cd2d3[idcoord];
          }
          tmp1 /= rBC;
          double tmp3 = tmp1*tmp1 + tmp2*tmp2;
          double BAdotcd2d3 = DotProduct(posA,posB,cd2d3,geometry);
          for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
            grad[0] = (CBcd2d3[idcoord]*tmp1 - cd2d3[idcoord]*tmp2/rBC)*ge/tmp3;
            grad[1] = (-BAdotcd2d3*rBC*tmp2*(geometry(posC,idcoord + 1) - geometry(posB,idcoord + 1)) + cd2d3[idcoord]*tmp2/rBC + BADC[idcoord]*tmp2/rBC + ACcd2d3[idcoord]*tmp1 - DCcd1d2[idcoord]*tmp1)*ge/tmp3;
            grad[2] = (BAdotcd2d3*rBC*tmp2*(geometry(posC,idcoord + 1) - geometry(posB,idcoord + 1)) + BABD[idcoord]/rBC*tmp2 + BAcd2d3[idcoord]*tmp1 - BDcd1d2[idcoord]*tmp1)*ge/tmp3;
            grad[3] = (CBcd1d2[idcoord]*tmp1 + cd1d2[idcoord]*tmp2/rBC)*ge/tmp3;
            gradient[3*(posA - 1) + idcoord] += grad[0];
            gradient[3*(posB - 1) + idcoord] += grad[1];
            gradient[3*(posC - 1) + idcoord] += grad[2];
            gradient[3*(posD - 1) + idcoord] += grad[3];
          }
        }
      }
    }
  }
  else {throw std::string("ERROR: Geometry.hpp - Internal2Cartesian(): strange gradient problem conversion to Cartesian");}
}
void Stretch(const matrixE & xyz, size_t idx1, size_t idx2, double & Avst, matrixE & Bmat, std::vector<matrixE> & dBmat, bool doDer = false) {
  //function to get stretching component for the Bmatrix, taken and adapted from xTB
  double Rarray[3];
  double R2;
  if ((Bmat.rows() != 3)&&(Bmat.cols() != 2)) {Bmat.resize(3,2);}
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    Rarray[idcoord] = xyz(idcoord + 1,idx1) - xyz(idcoord + 1,idx2);
  }
  R2 = Rarray[0]*Rarray[0] + Rarray[1]*Rarray[1] + Rarray[2]*Rarray[2];
  Avst = sqrt(R2);
  if (Avst > 1.0e-5) {R2 = 1.0/Avst;}
  else {R2 = Avst;}
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    Bmat(idcoord + 1,2) = Rarray[idcoord]*R2;
    Bmat(idcoord + 1,1) = -Bmat(idcoord + 1,2);        //translational invariance
  }
  if (doDer) {
    if (dBmat.size() != 4) {dBmat.resize(4);}
    for (size_t idx = 0; idx < 4; ++idx) {
      dBmat[idx].resize(3,3);
    }
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      for (size_t jdcoord = 0; jdcoord < idcoord + 1; ++jdcoord) {
        dBmat[0](idcoord + 1,jdcoord + 1) = -Bmat(jdcoord + 1,1)*Bmat(idcoord + 1,1)*R2;
        if (idcoord == jdcoord) {dBmat[0] += 1.0;}
        dBmat[0](jdcoord + 1,idcoord + 1) = dBmat[0](idcoord + 1,jdcoord + 1);
        dBmat[1](idcoord + 1,jdcoord + 1) = -dBmat[0](idcoord + 1,jdcoord + 1);
        dBmat[2](jdcoord + 1,idcoord + 1) = dBmat[1](idcoord + 1,jdcoord + 1);
        dBmat[2](idcoord + 1,jdcoord + 1) = -dBmat[0](idcoord + 1,jdcoord + 1);
        dBmat[1](jdcoord + 1,idcoord + 1) = dBmat[2](idcoord + 1,jdcoord + 1);
        dBmat[3](idcoord + 1,jdcoord + 1) = dBmat[1](idcoord + 1,jdcoord + 1);
        dBmat[3](jdcoord + 1,idcoord + 1) = dBmat[3](idcoord + 1,jdcoord + 1);
      }
    }
  }
}
void Bend(const matrixE & xyz, size_t idx1, size_t idx2, size_t idx3, double & fir, matrixE & Bmat, std::vector<matrixE> & dBmat, bool doDer = false, double threshzero1 = 1.0e-6, double threshzero2 = 1.0e-13) {
  //function to get bending component for the Bmatrix, taken and adapted from xTB
  double dfir;
  double Rij1;
  double Rjk1;
  double co = 0.0;
  double crap = 0.0;
  double sqrtcrap;
  double aux;
  double SI;
  matrixE BRij;
  matrixE BRjk;
  std::vector<matrixE> dBRij;
  std::vector<matrixE> dBRjk;
  if ((Bmat.rows() != 3)&&(Bmat.cols() != 3)) {Bmat.resize(3,3);}
  Stretch(xyz,idx1,idx2,Rij1,BRij,dBRij,doDer);
  Stretch(xyz,idx2,idx3,Rjk1,BRjk,dBRjk,doDer);
  SI = Rij1;
  if (SI > 1.0e-5) {Rij1 = 1.0/SI;}
  else {Rij1 = SI;}
  SI = Rjk1;
  if (SI > 1.0e-5) {Rjk1 = 1.0/SI;}
  else {Rjk1 = SI;}
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    co += BRij(idcoord + 1,1)*BRjk(idcoord + 1,2);
    aux = BRjk(idcoord + 1,2) + BRij(idcoord + 1,1);
    crap += aux*aux;
  }
  sqrtcrap = sqrt(crap);
  if (sqrtcrap < threshzero1) {               //special case for almost linear case
    fir = pi - asin(sqrtcrap);
    SI = sqrtcrap;
  }
  else {
    SI = fmin(co,1.0);
    aux = fmax(SI,-1.0);
    fir = acos(aux);
    SI = sqrt(1.0 - co*co);
  }
  if (fabs(fir - pi) < threshzero2) {fir = pi;}
  else {
    dfir = 180.0*fir/pi;
    //if ((fabs(dfir) > 177.5)||(fabs(dfir) < 2.5)) {std::cout << "WARNING: Geometry.hpp - Bend: angle close to end of range definition" << std::endl;}
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      Bmat(idcoord + 1,1) = (co*BRij(idcoord + 1,1) - BRjk(idcoord + 1,2))*Rij1/SI;
      Bmat(idcoord + 1,3) = (co*BRjk(idcoord + 1,2) - BRij(idcoord + 1,1))*Rjk1/SI;
      Bmat(idcoord + 1,2) = -(Bmat(idcoord + 1,1) + Bmat(idcoord + 1,3));                 //translational invariance
    }
    if (doDer) {
      if (dBmat.size() != 9) {dBmat.resize(9);}
      for (size_t idx = 0; idx < 9; ++idx) {
        dBmat[idx].resize(3,3);
      }
      //dBmat[0] -> (i,1,j,1)
      //dBmat[1] -> (i,2,j,1)
      //dBmat[2] -> (i,1,j,2)
      //dBmat[3] -> (i,2,j,2)
      //dBmat[4] -> (i,3,j,1)
      //dBmat[5] -> (i,1,j,3)
      //dBmat[6] -> (i,3,j,2)
      //dBmat[7] -> (i,2,j,3)
      //dBmat[8] -> (i,3,j,3)
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        for (size_t jdcoord = 0; jdcoord < idcoord + 1; ++jdcoord) {
          dBmat[0](idcoord + 1,jdcoord + 1) = -Bmat(idcoord + 1,1)*BRij(jdcoord + 1,1)*Rij1 + co*dBRij[0](idcoord + 1,jdcoord + 1)*Rij1/SI - Bmat(jdcoord + 1,1)*(co*Bmat(idcoord + 1,1)/SI + Rij1*BRij(idcoord + 1,1));
          dBmat[5](idcoord + 1,jdcoord + 1) = -Bmat(idcoord + 1,1)*BRjk(jdcoord + 1,2)*Rjk1 + dBRij[2](idcoord + 1,jdcoord + 1)*Rjk1/SI - Bmat(jdcoord + 1,3)*co*Bmat(idcoord + 1,1)/SI;
          dBmat[4](idcoord + 1,jdcoord + 1) = -Bmat(idcoord + 1,3)*BRij(jdcoord + 1,1)*Rij1 + dBRjk[1](idcoord + 1,jdcoord + 1)*Rij1/SI - Bmat(jdcoord + 1,1)*co*Bmat(idcoord + 1,3)/SI;
          dBmat[8](idcoord + 1,jdcoord + 1) = -Bmat(idcoord + 1,3)*BRjk(jdcoord + 1,2)*Rjk1 + co*dBRjk[2](idcoord + 1,jdcoord + 1)*Rjk1/SI - Bmat(jdcoord + 1,3)*(co*Bmat(idcoord + 1,3)/SI + Rjk1*BRjk(idcoord + 1,2));
          dBmat[0](jdcoord + 1,idcoord + 1) = dBmat[0](idcoord + 1,jdcoord + 1);
          dBmat[4](jdcoord + 1,idcoord + 1) = dBmat[5](idcoord + 1,jdcoord + 1);
          dBmat[5](jdcoord + 1,idcoord + 1) = dBmat[4](idcoord + 1,jdcoord + 1);
          dBmat[8](jdcoord + 1,idcoord + 1) = dBmat[8](idcoord + 1,jdcoord + 1);
          dBmat[2](idcoord + 1,jdcoord + 1) = -dBmat[0](idcoord + 1,jdcoord + 1) - dBmat[5](idcoord + 1,jdcoord + 1);
          dBmat[1](jdcoord + 1,idcoord + 1) = dBmat[2](idcoord + 1,jdcoord + 1);
          dBmat[2](jdcoord + 1,idcoord + 1) = -dBmat[0](jdcoord + 1,idcoord + 1) - dBmat[5](jdcoord + 1,idcoord + 1);
          dBmat[1](idcoord + 1,jdcoord + 1) = dBmat[2](jdcoord + 1,idcoord + 1);
          dBmat[6](idcoord + 1,jdcoord + 1) = -dBmat[4](idcoord + 1,jdcoord + 1) - dBmat[8](idcoord + 1,jdcoord + 1);
          dBmat[7](jdcoord + 1,idcoord + 1) = dBmat[6](idcoord + 1,jdcoord + 1);
          dBmat[6](jdcoord + 1,idcoord + 1) = -dBmat[4](jdcoord + 1,idcoord + 1) - dBmat[8](jdcoord + 1,idcoord + 1);
          dBmat[7](idcoord + 1,jdcoord + 1) = dBmat[6](jdcoord + 1,idcoord + 1);
          dBmat[3](idcoord + 1,jdcoord + 1) = -dBmat[1](idcoord + 1,jdcoord + 1) - dBmat[7](idcoord + 1,jdcoord + 1);
          dBmat[3](jdcoord + 1,idcoord + 1) = dBmat[3](idcoord + 1,jdcoord + 1);
        }
      }
    }
  }
}
void Torsion(const matrixE & xyz, double & tau, matrixE & Bmat, std::vector<matrixE> & dBmat, bool doDer = false) {
  //function to get torsion component for the Bmatrix, taken and adapted from xTB
  double Rij1;
  double Rjk1;
  double Rkl1;
  double Rijm1;
  double Rjkm1;
  double Rklm1;
  double phi2;
  double phi3;
  double sinp2;
  double cosp2;
  double sinp3;
  double cosp3;
  double cost;
  double sint;
  double aux;
  size_t jdcoord;
  size_t kdcoord;
  matrixE BRij;
  matrixE BRjk;
  matrixE BRkl;
  matrixE Bf2;
  matrixE Bf3;
  std::vector<matrixE> dBRij;
  std::vector<matrixE> dBRjk;
  std::vector<matrixE> dBRkl;
  std::vector<matrixE> dBf2;
  std::vector<matrixE> dBf3;
  if ((Bmat.rows() != 3)&&(Bmat.cols() != 4)) {Bmat.resize(3,4);}
  Stretch(xyz,2,1,Rij1,BRij,dBRij,doDer);
  Stretch(xyz,3,2,Rjk1,BRjk,dBRjk,doDer);
  Stretch(xyz,4,3,Rkl1,BRkl,dBRkl,doDer);
  if (Rij1 > 1.0e-5) {Rijm1 = 1.0/Rij1;}
  else {Rijm1 = Rij1;}
  if (Rjk1 > 1.0e-5) {Rjkm1 = 1.0/Rjk1;}
  else {Rjkm1 = Rjk1;}
  if (Rkl1 > 1.0e-5) {Rklm1 = 1.0/Rkl1;}
  else {Rklm1 = Rkl1;}
  Bend(xyz,3,2,1,phi2,Bf2,dBf2,doDer,1.0e-6,1.0e-13);
  sinp2 = sin(phi2);
  cosp2 = cos(phi2);
  Bend(xyz,4,3,2,phi3,Bf3,dBf3,doDer,1.0e-6,1.0e-13);
  sinp3 = sin(phi3);
  cosp3 = cos(phi3);
  //get the angle between the two planes using normal vectors
  cost = ((BRij(2,1)*BRjk(3,2) - BRij(3,1)*BRjk(2,2))*(BRjk(2,1)*BRkl(3,2) - BRjk(3,1)*BRkl(2,2)) + (BRij(3,1)*BRjk(1,2) - BRij(1,1)*BRjk(3,2))*(BRjk(3,1)*BRkl(1,2) - BRjk(1,1)*BRkl(3,2)) + (BRij(1,1)*BRjk(2,2) - BRij(2,1)*BRjk(1,2))*(BRjk(1,1)*BRkl(2,2) - BRjk(2,1)*BRkl(1,2)))/(sinp2*sinp3);
  //now define the sign
  sint = (BRij(1,2)*(BRjk(2,1)*BRkl(3,2) - BRjk(3,1)*BRkl(2,2)) + BRij(2,2)*(BRjk(3,1)*BRkl(1,2) - BRjk(1,1)*BRkl(3,2)) + BRij(3,2)*(BRjk(1,1)*BRkl(2,2) - BRjk(2,1)*BRkl(1,2)))/(sinp2*sinp3);
  tau = atan2(sint,cost);
  if (fabs(tau) == pi) {tau = pi;}
  aux = 180.0/pi;
  //if (fabs(tau*aux) > 177.5) {std::cout << "WARNING: Geometry.hpp - Torsion: dihedral close to end of range definition" << std::endl;}
  //if ((fabs(phi2*aux) > 177.5)||(fabs(phi2*aux) < 2.5)) {std::cout << "WARNING: Geometry.hpp - Torsion: bond angle close to end of range definition" << std::endl;}
  //if ((fabs(phi3*aux) > 177.5)||(fabs(phi3*aux) < 2.5)) {std::cout << "WARNING: Geometry.hpp - Torsion: bond angle close to end of range definition" << std::endl;}
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    jdcoord = idcoord + 1;
    if (jdcoord > 2) {jdcoord -= 3;}
    kdcoord = jdcoord + 1;
    if (kdcoord > 2) {kdcoord -= 3;}
    Bmat(idcoord + 1,1) = (BRij(jdcoord + 1,2)*BRjk(kdcoord + 1,2) - BRij(kdcoord + 1,2)*BRjk(jdcoord + 1,2))*Rijm1/(sinp2*sinp2);
    Bmat(idcoord + 1,4) = (BRkl(jdcoord + 1,1)*BRjk(kdcoord + 1,1) - BRkl(kdcoord + 1,1)*BRjk(jdcoord + 1,1))*Rklm1/(sinp3*sinp3);
    Bmat(idcoord + 1,2) = -((Rjk1 - Rij1*cosp2)*Bmat(idcoord + 1,1) + Rkl1*cosp3*Bmat(idcoord + 1,4))*Rjkm1;
    Bmat(idcoord + 1,3) = -Bmat(idcoord + 1,1) - Bmat(idcoord + 1,2) - Bmat(idcoord + 1,4);
  }
  if (doDer) {
    if (dBmat.size() != 16) {dBmat.resize(16);}
    for (size_t idx = 0; idx < 16; ++idx) {
      dBmat[idx].resize(3,4);
    }
    //dBmat[ 0] -> (i,1,j,1)
    //dBmat[ 1] -> (i,2,j,1)
    //dBmat[ 2] -> (i,1,j,2)
    //dBmat[ 3] -> (i,2,j,2)
    //dBmat[ 4] -> (i,3,j,1)
    //dBmat[ 5] -> (i,1,j,3)
    //dBmat[ 6] -> (i,3,j,2)
    //dBmat[ 7] -> (i,2,j,3)
    //dBmat[ 8] -> (i,3,j,3)
    //dBmat[ 9] -> (i,4,j,1)
    //dBmat[10] -> (i,1,j,4)
    //dBmat[11] -> (i,4,j,2)
    //dBmat[12] -> (i,2,j,4)
    //dBmat[13] -> (i,4,j,3)
    //dBmat[14] -> (i,3,j,4)
    //dBmat[15] -> (i,4,j,4)
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      jdcoord = idcoord + 1;
      for (size_t ldcoord = 0; ldcoord < idcoord + 1; ++ldcoord) {
        jdcoord = ldcoord + 1;
        if (jdcoord > 2) {jdcoord -= 3;}
        kdcoord = jdcoord + 1;
        if (kdcoord > 2) {kdcoord -= 3;}
        dBmat[0](idcoord + 1,ldcoord + 1) = (dBRij[2](idcoord + 1,jdcoord + 1)*BRjk(kdcoord + 1,2) - dBRij[2](idcoord + 1,kdcoord + 1)*BRjk(jdcoord + 1,2) - Bmat(ldcoord + 1,1)*(BRij(idcoord + 1,1)*sinp2*sinp2 + 2.0*Rij1*sinp2*cosp2*Bf2(idcoord + 1,1)))*Rijm1/(sinp2*sinp2);
        dBmat[2](idcoord + 1,ldcoord + 1) =-((-BRij(idcoord + 1,1)*cosp2 + Rij1*sinp2*Bf2(idcoord + 1,1))*Bmat(ldcoord + 1,1) + (Rjk1 - Rij1*cosp2)*dBmat[0](idcoord + 1,ldcoord + 1))*Rjkm1;
        dBmat[1](ldcoord + 1,idcoord + 1) = dBmat[2](idcoord + 1,ldcoord + 1);
        dBmat[10](idcoord + 1,ldcoord + 1) = 0.0;
        dBmat[9](ldcoord + 1,idcoord + 1) = dBmat[10](idcoord + 1,ldcoord + 1);
        dBmat[5](idcoord + 1,ldcoord + 1) = -(dBmat[0](idcoord + 1,ldcoord + 1) + dBmat[2](idcoord + 1,ldcoord + 1));
        dBmat[4](ldcoord + 1,idcoord + 1) = dBmat[5](idcoord + 1,ldcoord + 1);
        dBmat[15](idcoord + 1,ldcoord + 1) =(dBRkl[1](idcoord + 1,jdcoord + 1)*BRjk(kdcoord + 1,1) - dBRkl[1](idcoord + 1,kdcoord + 1)*BRjk(jdcoord + 1,1) - Bmat(ldcoord + 1,4)*(BRkl(idcoord + 1,2)*sinp3*sinp3 + 2.0*Rkl1*sinp3*cosp3*Bf3(idcoord + 1,3)))*Rklm1/(sinp3*sinp3);
        dBmat[13](idcoord + 1,ldcoord + 1) =-((-BRkl(idcoord + 1,2)*cosp3 + Rkl1*sinp3*Bf3(idcoord + 1,3))*Bmat(ldcoord + 1,4) + (Rjk1 - Rkl1*cosp3)*dBmat[15](idcoord + 1,ldcoord + 1))*Rjkm1;
        dBmat[14](ldcoord + 1,idcoord + 1) = dBmat[13](idcoord + 1,ldcoord + 1);
        dBmat[11](idcoord + 1,ldcoord + 1) = -(dBmat[15](idcoord + 1,ldcoord + 1) + dBmat[13](idcoord + 1,ldcoord + 1));
        dBmat[12](ldcoord + 1,idcoord + 1) = dBmat[11](idcoord + 1,ldcoord + 1);
        if (idcoord != ldcoord) {
          dBmat[0](ldcoord + 1,idcoord + 1) = dBmat[0](idcoord + 1,ldcoord + 1);
          dBmat[9](idcoord + 1,ldcoord + 1) = 0.0;
          dBmat[15](ldcoord + 1,idcoord + 1) = dBmat[15](idcoord + 1,ldcoord + 1);
          dBmat[10](ldcoord + 1,idcoord + 1) = dBmat[9](idcoord + 1,ldcoord + 1);
          dBmat[2](ldcoord + 1,idcoord + 1) = -((-BRij(ldcoord + 1,1)*cosp2 + Rij1*sinp2*Bf2(ldcoord + 1,1))*Bmat(idcoord + 1,1) + (Rjk1 - Rij1*cosp2)*dBmat[0](ldcoord + 1,idcoord + 1))*Rjkm1;
          dBmat[1](idcoord + 1,ldcoord + 1) = dBmat[2](ldcoord + 1,idcoord + 1);
          dBmat[4](idcoord + 1,ldcoord + 1) = -(dBmat[0](idcoord + 1,ldcoord + 1) + dBmat[1](idcoord + 1,ldcoord + 1) + dBmat[9](idcoord + 1,ldcoord + 1));
          dBmat[5](ldcoord + 1,idcoord + 1) = dBmat[4](idcoord + 1,ldcoord + 1);
          dBmat[13](ldcoord + 1,idcoord + 1) = -((-BRkl(ldcoord + 1,2)*cosp3 + Rkl1*sinp3*Bf3(ldcoord + 1,3))*Bmat(idcoord + 1,4) + (Rjk1 - Rkl1*cosp3)*dBmat[15](ldcoord + 1,idcoord + 1))*Rjkm1;
          dBmat[14](idcoord + 1,ldcoord + 1) = dBmat[13](ldcoord + 1,idcoord + 1);
          dBmat[12](idcoord + 1,ldcoord + 1) = -(dBmat[15](idcoord + 1,ldcoord + 1) + dBmat[14](idcoord + 1,ldcoord + 1));
          dBmat[11](ldcoord + 1,idcoord + 1) = dBmat[12](idcoord + 1,ldcoord + 1);
        }
        dBmat[7](idcoord + 1,ldcoord + 1) = -((BRjk(idcoord + 1,1) + Rkl1*sinp3*Bf3(idcoord + 1,1))*Bmat(ldcoord + 1,4) + (Rjk1 - Rkl1*cosp3)*dBmat[12](idcoord + 1,ldcoord + 1) + (BRij(idcoord + 1,2)*cosp2 - Rij1*sinp2*Bf2(idcoord + 1,2))*Bmat(ldcoord + 1,1) + Rij1*cosp2*dBmat[1](idcoord + 1,ldcoord + 1) + Bmat(ldcoord + 1,3)*BRjk(idcoord + 1,1))*Rjkm1;
        dBmat[6](ldcoord + 1,idcoord + 1) = dBmat[7](idcoord + 1,ldcoord + 1);
        dBmat[3](idcoord + 1,ldcoord + 1) = -(dBmat[1](idcoord + 1,ldcoord + 1) + dBmat[12](idcoord + 1,ldcoord + 1) + dBmat[7](idcoord + 1,ldcoord + 1));
        dBmat[8](idcoord + 1,ldcoord + 1) = -(dBmat[7](idcoord + 1,ldcoord + 1) + dBmat[5](idcoord + 1,ldcoord + 1) + dBmat[13](idcoord + 1,ldcoord + 1));
        if (idcoord != ldcoord) {
          dBmat[6](idcoord + 1,ldcoord + 1) = -(dBmat[3](idcoord + 1,ldcoord + 1) + dBmat[2](idcoord + 1,ldcoord + 1) + dBmat[11](idcoord + 1,ldcoord + 1));
          dBmat[7](ldcoord + 1,idcoord + 1) = dBmat[6](idcoord + 1,ldcoord + 1);
          dBmat[3](ldcoord + 1,idcoord + 1) = dBmat[3](idcoord + 1,ldcoord + 1);
          dBmat[8](ldcoord + 1,idcoord + 1) = dBmat[8](idcoord + 1,ldcoord + 1);
        }
      }
    }
  }
}
void OutOfPlane(const matrixE & xyz, double & tau, matrixE & Bmat, double thresholdzero = 1.0e-13) {
  //function to get out-of-plane component for the Bmatrix, taken and adapted from xTB
  double r14[3];
  double r24[3];
  double r34[3];
  double cosp = 0.0;      //angle between unit vectors in the directions r34 and r24
  double cosp2 = 0.0;     //angle between unit vectors in the directions r34 and r14
  double cosp3 = 0.0;     //angle between unit vectors in the directions r24 and r14
  double phi;
  double phi2;
  double phi3;
  double R14 = 0.0;
  double R24 = 0.0;
  double R34 = 0.0;
  double aux;
  matrixE Centers(3,3);
  if ((Bmat.rows() != 3)&&(Bmat.cols() != 4)) {Bmat.resize(3,4);}
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    r14[idcoord] = xyz(idcoord + 1,1) - xyz(idcoord + 1,4);
    r24[idcoord] = xyz(idcoord + 1,2) - xyz(idcoord + 1,4);
    r34[idcoord] = xyz(idcoord + 1,3) - xyz(idcoord + 1,4);
    R14 += r14[idcoord]*r14[idcoord];
    R24 += r24[idcoord]*r24[idcoord];
    R34 += r34[idcoord]*r34[idcoord];
    for (size_t idcol = 1; idcol < 5; ++idcol) {
      Bmat(idcoord + 1,idcol) = 0.0;
    }
  }
  aux = sqrt(R14);
  R14 = aux;
  aux = sqrt(R24);
  R24 = aux;
  aux = sqrt(R34);
  R34 = aux;
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    r14[idcoord] /= R14;
    r24[idcoord] /= R24;
    r34[idcoord] /= R34;
    cosp += r34[idcoord]*r24[idcoord];
    cosp2 += r34[idcoord]*r14[idcoord];
    cosp3 += r24[idcoord]*r14[idcoord];
  }
  phi = acos(cosp);
  if (fabs(phi - pi) < thresholdzero) {tau = 0.0;}       //bad exit, when the structure is not well defined
  else {
    phi2 = acos(cosp2);
    phi3 = acos(cosp3);
    //get the centers
    for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
      Centers(idcoord,1) = xyz(idcoord,1);
      Centers(idcoord,2) = xyz(idcoord,4);
      r24[idcoord - 1] = xyz(idcoord,2) - xyz(idcoord,4);
      r34[idcoord - 1] = xyz(idcoord,3) - xyz(idcoord,4);
    }
    Centers(1,3) = r24[1]*r34[2] - r24[2]*r34[1];
    Centers(2,3) = r24[2]*r34[0] - r24[0]*r34[2];
    Centers(3,3) = r24[0]*r34[1] - r24[1]*r34[0];
    aux = Centers(1,3)*Centers(1,3) + Centers(2,3)*Centers(2,3) + Centers(3,3)*Centers(3,3);
    if (aux < thresholdzero) {tau = 0.0;}            //check whether 2,3,4 are colinear
    else {
      for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
        Centers(idcoord,3) += xyz(idcoord,4);
      }
      matrixE BR14;
      std::vector<matrixE> dBR14;
      int jdcoord;
      int kdcoord;
      Bend(Centers,1,2,3,tau,BR14,dBR14,false,1.0e-6,1.0e-13);
      BR14 *= -1.0;
      tau -= 0.5*pi;
      for (int idcoord = 1; idcoord < 4; ++idcoord) {
        jdcoord = int(fmod(idcoord + 1,4)) + (idcoord + 1)/4;
        kdcoord = int(fmod(jdcoord + 1,4)) + (jdcoord + 1)/4;
        Bmat(idcoord,1) = BR14(idcoord,1);
        Bmat(idcoord,2) = -r34[kdcoord - 1]*BR14(jdcoord,3) + r34[jdcoord - 1]*BR14(kdcoord,3);
        Bmat(idcoord,3) = r24[kdcoord - 1]*BR14(jdcoord,3) - r24[jdcoord - 1]*BR14(kdcoord,3);
        Bmat(idcoord,4) = -(Bmat(idcoord,1) + Bmat(idcoord,2) + Bmat(idcoord,3));
      }
    }
  }
}
void RotateGeomMatrix(size_t axis, double theta, matrixE & geometry, matrixE & rotation) {
  //function that rotates geometry according to certain Given's transformation matrix
  if (axis > 3) {throw std::string("ERROR: Geometry.hpp: RotateGeomMatrix(): Ndimensions < 4");}
  size_t i = 2;
  size_t j = 1;
  if (axis == 2) {
    i = 1;
    j = 3;
  }
  else if (axis == 1) {
    i = 3;
    j = 2;
  }
  rotation = Givens(i,j,3,cos(theta),sin(theta));
  geometry *= rotation.trans();
}
void GetLocalAxes(std::vector<double> & u21, std::vector<double> & u23, std::vector<double> & cp) {
  //function that gets the local axes for the atom
  double aux = dotProd(u23,u21);
  cp = crossProd(u23,u21,sqrt(1.0 - aux*aux));
  aux = dotProd(cp,u21);
  u23 = crossProd(u21,cp,sqrt(1.0 - aux*aux));        //u23 = u21 X u23 X u21
  aux = dotProd(u23,u21);
  cp = crossProd(u23,u21,sqrt(1.0 - aux*aux));
  swap(cp[1],u23[0]);
  swap(cp[2],u21[0]);
  swap(u23[2],u21[1]);
}
void GetBondVector(std::vector<double> & bvec, matrixE & zmat, int idatm) {
  //function that gets a bond vector from the zmatrix
  double rij = zmat(idatm,2);
  double theta = zmat(idatm,4);
  double phi = zmat(idatm,6);
  bvec[0] = rij*sin(theta)*sin(phi);
  bvec[1] = rij*sin(theta)*cos(phi);
  bvec[2] = rij*cos(theta);
}
void GetUnitVectors(std::vector<double> & u21, std::vector<double> & u23, matrixE & cart, int atmA, int atmR, int atmD) {
  //getting the normalized unit vectors between relevant atoms for defining coordinates
  double nu21 = 0.0;
  double nu23 = 0.0;
  for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
    u21[idcoord - 1] = cart(atmA,idcoord) - cart(atmR,idcoord);
    nu21 += u21[idcoord - 1]*u21[idcoord - 1];
    u23[idcoord - 1] = cart(atmD,idcoord) - cart(atmA,idcoord);
    nu23 += u23[idcoord - 1]*u23[idcoord - 1];
  }
  double aux = sqrt(nu21);
  nu21 = 1.0/aux;
  aux = sqrt(nu23);
  nu23 = 1.0/aux;
  for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
    u21[idcoord] *= nu21;
    u23[idcoord] *= nu23;
  }
}
void initAlign(matrixE & align, int sz) {
  //just an initializer
  align.resize(sz,2);
  for (size_t idAtm = 0; idAtm < sz; ++idAtm) {
    align(idAtm + 1,1) = idAtm + 1;
    align(idAtm + 1,2) = idAtm + 1;
  }
}
void AlignBonA(matrixE & geomA, matrixE & geomB, matrixE & align, int firstN = -1) {
  //function that aligns the geometry of system B according to systemA, using for the alignment all the atoms given in align
  //alternatively, give a value to firstN and this will use the first firstN atoms in geometries A and B for the alignment
  //the alignment is based on SVD decomposition
  //align is given as a matrix, so that we may choose non-natural alignment combinations
  //the structure is as follows (atomA,atomB) with 1 based indices are to be aligned
  //a meaningful align matrix is kind of needed
  if ((firstN > 0)&&(align.cols() != 2)) {initAlign(align,firstN);}
  else if ((firstN < 0)&&(align.cols() != 2)) {initAlign(align,geomA.rows());}
  firstN = align.rows();
  std::vector<double> cmA(3,0.0);
  std::vector<double> cmB(3,0.0);
  matrixE newgeomA(firstN,3);
  matrixE newgeomB(firstN,3);
  //get the geometries to align and their centroids
  for (size_t idAtm = 0; idAtm < firstN; ++idAtm) {
    for (size_t icoord = 0; icoord < 3; ++icoord) {
      newgeomA(idAtm + 1,icoord + 1) = geomA(int(align(idAtm + 1,1)),icoord + 1);
      newgeomB(idAtm + 1,icoord + 1) = geomB(int(align(idAtm + 1,2)),icoord + 1);
      cmA[icoord] += newgeomA(idAtm + 1,icoord + 1);
      cmB[icoord] += newgeomB(idAtm + 1,icoord + 1);
    }
  }
  for (size_t icoord = 0; icoord < 3; ++icoord) {
    cmA[icoord] /= firstN;
    cmB[icoord] /= firstN;
  }
  //shift the geometries to the centroids
  for (size_t idAtm = 0; idAtm < firstN; ++idAtm) {
    for (size_t icoord = 0; icoord < 3; ++icoord) {
      newgeomA(idAtm + 1,icoord + 1) -= cmA[icoord];
      newgeomB(idAtm + 1,icoord + 1) -= cmB[icoord];
    }
  }
  //get covariance
  matrixE cov = newgeomA.trans()*newgeomB;
  matrixE Umat;
  matrixE Sval;
  matrixE Vmat;
  //now SVD it
  SVD(cov,Umat,Sval,Vmat);
  newgeomA = Vmat*Umat.trans();           //the rotation matrix
  double detR = newgeomA.det();
  if (detR < 0) {
    for (size_t idc = 0; idc < Vmat.cols(); ++idc) {
      Vmat(3,idc + 1) *= -1.0;
    }
    newgeomA = Vmat*Umat.trans();
  }
  newgeomB = cmA;
  Sval = cmB;
  Vmat = Sval.trans()*newgeomA;
  newgeomB -= Vmat.trans();               //the translation matrix
  //align the structure of B with respect to A
  geomB *= newgeomA;
  int NatmB = geomB.rows();
  for (size_t idBtm = 0; idBtm < NatmB; ++idBtm) {
    geomB(idBtm + 1,1) += newgeomB(1,1);
    geomB(idBtm + 1,2) += newgeomB(2,1);
    geomB(idBtm + 1,3) += newgeomB(3,1);
  }
}

#endif //_Geometry_
