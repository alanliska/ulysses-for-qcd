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

#ifndef _D3_Functions_
#define _D3_Functions_
#include <vector>
#include <string>
#include <math.h>
#include "../ConstantsPackage.hpp"
#include "../UnitConversion.hpp"
#include "../math/GaussianPackage.hpp"
#include "../math/IntegralPackage.hpp"
#include "../math/MatrixPackage.hpp"
#include "../math/TensorPackage.hpp"
#include "../parameters/D3par.hpp"
#include "../parameters/D4par.hpp"
#include "../parameters/GFN2par.hpp"

//description:
//set of functions that calculate contribution of intermolecular forces to total energy according to the D3 and D4 methods
//applicable to semi-empirical methods and DFT
//references:
//D3 -> S. Grimme, J. Antony, S. Ehrlich, H. Krieg, J. Chem. Phys., 132, 154104, 2010
//D4 -> E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher, C. Bannwarth, S. Grimme, J. Chem. Phys., 150, 154122, 2019

//counting functions
double GFNcount(double kcn, double xx, double x0) {
  //the GFN type of counting; xx is the interatomic distance (Angstrom) and x0 is a reference distance (Angstrom)
  return 1.0/((1.0 + exp(-kcn*(x0/xx - 1.0)))*(1.0 + exp(-2.0*kcn*((x0 + 2.0*dist_Angstrom2au)/xx - 1.0))));
}
double dGFNcountdx(double kcn, double xx, double x0) {
  //first derivative of the GFN type of counting
  double expterm1 = exp(-kcn*(x0/xx - 1.0));
  double x0prime = x0 + 2.0*dist_Angstrom2au;
  double expterm2 = exp(-2.0*kcn*(x0prime/xx - 1.0));
  double result = kcn*x0*expterm1/(xx*xx*(1.0 + expterm1)*(1.0 + expterm1)*(1.0 + expterm2));
  result += 2.0*kcn*x0prime*expterm2/(xx*xx*(1.0 + expterm1)*(1.0 + expterm2)*(1.0 + expterm2));
  return -result;
}
double ERFcount(double kcn, double xx, double x0) {
  //the error function type of counting; xx is the interatomic distance (Angstrom) and x0 is a reference distance (Angstrom)
  return 0.5*(1.0 + erf(-kcn*(xx/x0 - 1.0)));
}
double dERFcountdx(double kcn, double xx, double x0) {
  //first derivative of the error function type of counting
  double exponent = kcn*(xx/x0 - 1.0);
  return -kcn*exp(-exponent*exponent)/(sqrt(pi)*x0);
}
double EXPcount(double kcn, double xx, double x0) {
  //the exponential function type of counting; xx is the interatomic distance (Angstrom) and x0 is a reference distance (Angstrom)
  return 1.0/(1.0 + exp(-16.0*(x0/xx - 1.0)));
}
double dEXPcountdx(double kcn, double xx, double x0) {
  //first derivative of the exponential function type of counting
  double expterm = exp(-16.0*(x0/xx - 1.0));
  double denom = 1.0/(xx*(1.0 + expterm));
  return -16.0*x0*expterm*denom*denom;
}
//coordination numbers
void NCoordMat(matrixE & ncoord, std::vector<size_t> & atoms, matrixE & geometry, double (*count)(double, double, double), bool useden = false, double d2threshold = 1600.0, bool ERFfunctionuse = false) {
  //function calculating the coordination number; this is the same as Grimme's coordination number function for D3 or D4, but the output comes in matrix form
  //note that in comparison to paper the code has no factor 4/3; this is because the factor was included in the rcov parameters already
  int Natoms = atoms.size();
  ncoord.resize(Natoms,Natoms);
  double RAcov;
  double distance;
  double aux1;
  double aux2;
  double den;
  double ENA;
  double k4 = k4D4();
  double k5 = k5D4();
  double k6 = 1.0/k6D4();
  double kcn = 10.0*double((!useden)&&(!ERFfunctionuse)) + double((useden)||(ERFfunctionuse))*knD4();
  double shift = 2.0*double(useden);
  d2threshold *= dist_Angstrom2au*dist_Angstrom2au;
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    RAcov = covRadius(atoms[idAtm]);
    ENA = PaulingElectronegativity(atoms[idAtm]);
    for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
      ncoord(idAtm + 1,idBtm + 1) = 0.0;
      ncoord(idBtm + 1,idAtm + 1) = 0.0;
      aux1 = 0.0;
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        distance = geometry(idAtm + 1,idcoord + 1) - geometry(idBtm + 1,idcoord + 1);
        aux1 += distance*distance;
      }
      if (aux1 < d2threshold) {
        distance = sqrt(aux1);
        den = 1.0;
        if (useden) {
          aux1 = fabs(ENA - PaulingElectronegativity(atoms[idBtm])) + k5;
          den = k4*exp(-aux1*aux1*k6);
        }
        aux2 = RAcov + covRadius(atoms[idBtm]);
        aux1 = den*(*count)(kcn,distance,aux2);
        ncoord(idAtm + 1,idBtm + 1) += aux1;
        ncoord(idBtm + 1,idAtm + 1) += aux1;
      }
    }
  }
}
void NCoord(std::vector<double> & ncoord, std::vector<size_t> & atoms, matrixE & geometry, double (*count)(double, double, double), bool useden = false, double d2threshold = 1600.0, bool ERFfunctionuse = false) {
  //function calculating the coordination number
  //note that in comparison to paper the code has no factor 4/3; this is because the factor was included in the rcov parameters already
  int Natoms = atoms.size();
  ncoord.resize(Natoms);
  double RAcov;
  double distance;
  double aux1;
  double aux2;
  double den;
  double ENA;
  double k4 = k4D4();
  double k5 = k5D4();
  double k6 = 1.0/k6D4();
  double kcn = 10.0*double((!useden)&&(!ERFfunctionuse)) + double((useden)||(ERFfunctionuse))*knD4();
  double shift = 2.0*double(useden);
  d2threshold *= dist_Angstrom2au*dist_Angstrom2au;
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    RAcov = covRadius(atoms[idAtm]);
    ENA = PaulingElectronegativity(atoms[idAtm]);
    ncoord[idAtm] = 0.0;
    for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
      aux1 = 0.0;
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        distance = geometry(idAtm + 1,idcoord + 1) - geometry(idBtm + 1,idcoord + 1);
        aux1 += distance*distance;
      }
      if (aux1 < d2threshold) {
        distance = sqrt(aux1);
        den = 1.0;
        if (useden) {
          aux1 = fabs(ENA - PaulingElectronegativity(atoms[idBtm])) + k5;
          den = k4*exp(-aux1*aux1*k6);
        }
        aux2 = RAcov + covRadius(atoms[idBtm]);
        aux1 = den*(*count)(kcn,distance,aux2);
        ncoord[idAtm] += aux1;
        ncoord[idBtm] += aux1;
      }
    }
  }
}
double dNCoorddR(int AtmnrA, int AtmnrB, double RAB, double (*dcountdx)(double, double, double), bool useden = false, double d2threshold = 1600.0, bool ERFfunctionuse = false) {
  //function calculating the derivative of the coordination number with respect to internuclear distance
  double dncoorddr = 0.0;
  d2threshold *= dist_Angstrom2au*dist_Angstrom2au;
  if (RAB*RAB < d2threshold) {
    double RABcov = covRadius(AtmnrA) + covRadius(AtmnrB);
    double kcn = 10.0*double((!useden)&&(!ERFfunctionuse)) + double((useden)||(ERFfunctionuse))*knD4();
    dncoorddr = (*dcountdx)(kcn,RAB,RABcov)*dist_Angstrom2au;
  }
  return dncoorddr;
}
//functions related to the D3 dispersion correction
double C6ABD3(size_t iatm1, size_t iatm2, double nc1, double nc2, size_t iatmnr1, size_t iatmnr2, double k3 = 4.0, double threshold = 1.0e-10) {
  //function returning the C6AB parameters, i.e., for a given pair of atoms
  double c6ab;
  double c6rmin = -1.0e10;
  double rmin = 1.0e10;
  double rsum = 0.0;
  double cnsum = 0.0;
  double rr;
  double ncA;
  double ncB;
  double expfactor;
  size_t iatmA;
  size_t iatmB;
  size_t iatmnrA;
  size_t iatmnrB;
  if (iatmnr1 < iatmnr2) {
    iatmnrA = iatmnr1;
    iatmA = iatm1;
    ncA = nc1;
    iatmnrB = iatmnr2;
    iatmB = iatm2;
    ncB = nc2;
  }
  else {
    iatmnrA = iatmnr2;
    iatmA = iatm2;
    ncA = nc2;
    iatmB = iatm1;
    iatmnrB = iatmnr1;
    ncB = nc1;
  }
  //get reference coordination numbers
  std::vector<double> cnatmA = CN_D3(iatmnrA);
  std::vector<double> cnatmB = CN_D3(iatmnrB);
  //get the respective C6 reference parameters
  matrixE c6data(cnatmA.size(),cnatmB.size());
  C6parametersD3(c6data,iatmnrA,iatmnrB);
  for (size_t idA = 0; idA < cnatmA.size(); ++idA) {
    for (size_t idB = 0; idB < cnatmB.size(); ++idB) {
      rr = (cnatmA[idA] - ncA)*(cnatmA[idA] - ncA) + (cnatmB[idB] - ncB)*(cnatmB[idB] - ncB);
      if (rr < rmin) {
        rmin = rr;
        c6rmin = c6data(idA + 1,idB + 1);
      }
      expfactor = exp(-k3*rr);
      rsum += expfactor;
      cnsum += expfactor*c6data(idA + 1,idB + 1);
    }
  }
  if (rsum > threshold) {c6ab = cnsum/rsum;}
  else {c6ab = c6rmin;}
  return c6ab;
}
void dC6dCNij(size_t iatmnr1, size_t iatmnr2, double nc1, double nc2, double & dc6A, double & dc6B, double k3 = 4.0, double thresholdzero = 1.0e-10) {
  //function getting the derivative of the D3 C6 coefficients with respect to coordination number
  size_t iatmnrA;
  size_t iatmnrB;
  double ncA;
  double ncB;
  double C6reference;
  double rr;
  double expfactor;
  double rsum = 0.0;
  double drsumA = 0.0;
  double drsumB = 0.0;
  double cnsum = 0.0;
  double dcnsumA = 0.0;
  double dcnsumB = 0.0;
  if (iatmnr1 < iatmnr2) {
    iatmnrA = iatmnr1;
    ncA = nc1;
    iatmnrB = iatmnr2;
    ncB = nc2;
  }
  else {
    iatmnrA = iatmnr2;
    ncA = nc2;
    iatmnrB = iatmnr1;
    ncB = nc1;
  }
  std::vector<double> cnatmA = CN_D3(iatmnrA);
  std::vector<double> cnatmB = CN_D3(iatmnrB);
  size_t szcnA = cnatmA.size();
  size_t szcnB = cnatmB.size();
  matrixE c6data(szcnA,szcnB);
  dc6A = 0.0;
  dc6B = 0.0;
  C6parametersD3(c6data,iatmnrA,iatmnrB);
  for (size_t idelA = 0; idelA < szcnA; ++idelA) {
    for (size_t idelB = 0; idelB < szcnB; ++idelB) {
      C6reference = c6data(idelA + 1,idelB + 1);
      if (C6reference > thresholdzero) {
        rr = (cnatmA[idelA] - ncA)*(cnatmA[idelA] - ncA) + (cnatmB[idelB] - ncB)*(cnatmB[idelB] - ncB);
        expfactor = exp(-k3*rr);
        rsum += expfactor;
        cnsum += expfactor*C6reference;
        expfactor *= 2.0*k3;
        rr = -expfactor*(cnatmA[idelA] - ncA);
        dcnsumA += C6reference*rr;
        drsumA += rr;
        rr = -expfactor*(cnatmB[idelB] - ncB);
        dcnsumB += C6reference*rr;
        drsumB += rr;
      }
    }
  }
  if (rsum > thresholdzero) {
    expfactor = -1.0/(rsum*rsum);
    ncA = (dcnsumA*rsum - drsumA*cnsum)*expfactor;
    ncB = (dcnsumB*rsum - drsumB*cnsum)*expfactor;
    if (iatmnr1 < iatmnr2) {
      dc6A = ncA;
      dc6B = ncB;
    }
    else {
      dc6B = ncA;
      dc6A = ncB;
    }
  }
}
double C8AB(double C6AB, double r2r4A, double r2r4B) {
  //function calculating the C8 coefficients
  return 3.0*r2r4A*r2r4B*C6AB;
}
double fdamp(double rAB, size_t atom1, size_t atom2, const std::string & method, int order) {
  //function calculating the damping term
  return 1.0/(1.0 + 6.0*pow(rAB/(R0AB(atom1,atom2)*sR(method,order)),-AlphaD3(method,order)));
}
double gfdamp(double rAB, size_t atom1, size_t atom2, const std::string & method, int order) {
  //function calculating the first-derivative of the damping term
  double fdamping = fdamp(rAB,atom1,atom2,method,order);
  double alpha = AlphaD3(method,order);
  return 6.0*alpha*pow(sR(method,order)*R0AB(atom1,atom2)/rAB,alpha)*fdamping*fdamping/rAB;
}
double hfdamp(double rAB, size_t atom1, size_t atom2, const std::string & method, int order) {
  //function calculating the second-derivative of the damping term
  double fdamping = fdamp(rAB,atom1,atom2,method,order);
  double gfdamping = gfdamp(rAB,atom1,atom2,method,order);
  double alpha = AlphaD3(method,order);
  return 6.0*alpha*pow(sR(method,order)*R0AB(atom1,atom2)/rAB,alpha)*fdamping*(2.0*gfdamping - (alpha + 1.0)*fdamping/rAB)/rAB;
}
double EdispC6(double rAB, double c6ij, size_t atom1, size_t atom2, const std::string & method, double bjdamp = -1.0) {
  //function for the first dispersion correction D3
  //convert rAB into au
  rAB *= dist_Angstrom2aum1;
  double rab3 = rAB*rAB*rAB;
  double e6;
  if (bjdamp > 0.0) {                                               //do Becke-Johnson damping
    e6 = bjdamp*A1D3(method) + A2D3(method);
    double damp = e6*e6*e6;
    e6 = c6ij/(rab3*rab3 + damp*damp);
  }
  else {e6 = c6ij*fdamp(rAB,atom1,atom2,method,6)/(rab3*rab3);}     //zero damping
  return e6;
}
//dE(2)/dRAB = sn*[dCnAB/dRAB]*f(RAB)/(RAB^n)   (1)
//           + sn*CnAB*[df(RAB)/dRAB]/(RAB^n)   (2)
//           - n*sn*CnAB*f(RAB)/(RAB^(n + 1))   (3)
double gEdispC6(double rAB, double c6ij, size_t atom1, size_t atom2, const std::string & method, double bjdamp = -1.0) {
  //function for the first-derivatives of the first dispersion correction D3
  //this is the coordination number independent part of the derivative, i.e., terms (2) and (3)
  rAB *= dist_Angstrom2aum1;
  double rab3 = rAB*rAB*rAB;
  double ge6;
  if (bjdamp > 0.0) {                                               //do Becke-Johnson damping
    ge6 = bjdamp*A1D3(method) + A2D3(method);
    double damp = ge6*ge6*ge6;
    ge6 = 6.0*c6ij*rab3*rAB*rAB*dist_Angstrom2au/(dist_Angstrom2au*(rab3*rab3 + damp*damp)*(rab3*rab3 + damp*damp));
  }
  else {ge6 = -c6ij*(gfdamp(rAB,atom1,atom2,method,6) - 6.0*fdamp(rAB,atom1,atom2,method,6)/rAB)/(rab3*rab3);}     //zero damping
  return ge6;
}
void dC6dRAB(double & dc6A, double & dc6B, int AtmnrA, int AtmnrB, double RAB, double CNA, double CNB, double (*count)(double, double, double), bool useden = false, double d2threshold = 1600.0, bool ERFfunctionuse = false) {
  //dispersion coefficient part of the first derivative of E6, i.e., term (1); this does only contain dC6/dRAB = (dC6/dCNA)(dCNA/dRAB) + (dC6/dCNB)(dCNB/dRAB) = [(dC6/dCNA) + (dC6/dCNB)](dCNA/dRAB) since (dCNA/dRAB) = (dCNB/dRAB)
  dc6A = 0.0;
  dc6B = 0.0;
  dC6dCNij(AtmnrA,AtmnrB,CNA,CNB,dc6A,dc6B);
}
double hEdispC6(double rAB, double c6ij, size_t atom1, size_t atom2, const std::string & method, double bjdamp = -1.0) {
  //function for the second-derivatives of the first dispersion correction D3
  //convert rAB into au
  rAB *= dist_Angstrom2aum1;
  double rab2 = rAB*rAB;
  double rab3 = rab2*rAB;
  double he6;
  if (bjdamp > 0.0) {                                               //do Becke-Johnson damping
    he6 = bjdamp*A1D3(method) + A2D3(method);
    double damp = he6*he6*he6;
    double xx = rab3*rab3 + damp*damp;
    he6 = c6ij*rab2*rab2*(30.0 - 72.0*rab3*rab3/xx)/(xx*xx*dist_Angstrom2au*dist_Angstrom2au);
  }
  else {he6 = -c6ij*(hfdamp(rAB,atom1,atom2,method,6) - 12.0*gfdamp(rAB,atom1,atom2,method,6)/rAB + 42.0*fdamp(rAB,atom1,atom2,method,6)/rab2)/(rab3*rab3*dist_Angstrom2au*dist_Angstrom2au);}     //zero damping
  return he6;
}
double EdispC8(double rAB, double c8ij, size_t atom1, size_t atom2, const std::string & method, double bjdamp = -1.0) {
  //function for the second dispersion correction D3
  //convert rAB into au
  rAB *= dist_Angstrom2aum1;
  double rab4 = rAB*rAB*rAB*rAB;
  double e8;
  if (bjdamp > 0.0) {                                               //do Becke-Johnson damping
    e8 = bjdamp*A1D3(method) + A2D3(method);
    double damp = e8*e8*e8*e8;
    e8 = c8ij/(rab4*rab4 + damp*damp);
  }
  else {e8 = c8ij*fdamp(rAB,atom1,atom2,method,8)/(rab4*rab4);}     //zero damping
  return e8;
}
double gEdispC8(double rAB, double c8ij, size_t atom1, size_t atom2, const std::string & method, double bjdamp = -1.0) {
  //function for the first-derivatives of the second dispersion correction D3
  //convert rAB into au
  rAB *= dist_Angstrom2aum1;
  double rab4 = rAB*rAB*rAB*rAB;
  double ge8;
  if (bjdamp > 0.0) {                                               //do Becke-Johnson damping
    ge8 = bjdamp*A1D3(method) + A2D3(method);
    double damp = ge8*ge8*ge8*ge8;
    ge8 = 8.0*c8ij*rab4*rAB*rAB*rAB/((rab4*rab4 + damp*damp)*(rab4*rab4 + damp*damp));
  }
  else {ge8 = -c8ij*(gfdamp(rAB,atom1,atom2,method,8) - 8.0*fdamp(rAB,atom1,atom2,method,8)/rAB)/(rab4*rab4);}     //zero damping
  return ge8;
}
double hEdispC8(double rAB, double c8ij, size_t atom1, size_t atom2, const std::string & method, double bjdamp = -1.0) {
  //function for the second-derivatives of the second dispersion correction D3
  //convert rAB into au
  rAB *= dist_Angstrom2aum1;
  double rab2 = rAB*rAB;
  double rab4 = rab2*rab2;
  double he8;
  if (bjdamp > 0.0) {                                               //do Becke-Johnson damping
    he8 = bjdamp*A1D3(method) + A2D3(method);
    double damp = he8*he8*he8*he8;
    double yy = rab4*rab4 + damp*damp;
    he8 = c8ij*rab4*rab2*(56.0 - 128.0*rab4*rab4/yy)/(yy*yy*dist_Angstrom2au*dist_Angstrom2au);
  }
  else {he8 = -c8ij*(hfdamp(rAB,atom1,atom2,method,8) - 16.0*gfdamp(rAB,atom1,atom2,method,8)/rAB + 72.0*fdamp(rAB,atom1,atom2,method,8)/rab2)/(rab4*rab4*dist_Angstrom2au*dist_Angstrom2au);}     //zero damping
  return he8;
}
//functions that calculate the Axilrod-Teller-Muto terms
double EABC(double & c6AB, double & c6AC, double & c6BC, double & RcritAB, double & RcritAC, double & RcritBC, double & r2AB, double & r2AC, double & r2BC, double & alpha) {
  double c9 = sqrt(c6AB*c6AC*c6BC);                                                                                                            //equation 24 of D4 paper
  double R2ABACBC = r2AB*r2AC*r2BC;
  double Ravg = sqrt(R2ABACBC);
  double ATM = c9*(0.375*(r2AB + r2BC - r2AC)*(r2AB - r2BC + r2AC)*(-r2AB + r2BC + r2AC)/(Ravg*R2ABACBC*R2ABACBC) + 1.0/(Ravg*R2ABACBC));      //equation 22 of D4 paper
  return ATM/(1.0 + 6.0*pow(RcritAB*RcritAC*RcritBC/Ravg,alpha/3.0));                                                                          //single term of equation 25, contains implicitly also equations 26, 27
}
void calcATM(const std::string & method, std::vector<double> & eABC, std::vector<size_t> & atoms, matrixE & Geometry, matrixE & c6mat, double d2threshold = 1600.0) {
  //function that calculates the atomwise contribution to ATM thre-body term
  size_t Natoms = atoms.size();
  eABC.resize(Natoms);
  int atomA;
  int atomB;
  int atomC;
  double c6AB;
  double c6AC;
  double c6BC;
  double cAB;
  double cAC;
  double cBC;
  double sqrtZr4r2A;
  double sqrtZr4r2B;
  double sqrtZr4r2C;
  double a1 = A1D3(method);
  double a2 = A2D3(method);
  double RAB2;
  double RAC2;
  double RBC2;
  double dE;
  double alpha = AlphaD3(method,6);
  double S9 = s9D3(method)/3.0;
  std::vector<double> RA(3,0.0);
  std::vector<double> RB(3,0.0);
  d2threshold *= dist_Angstrom2au*dist_Angstrom2au;
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    atomA = atoms[idAtm];
    sqrtZr4r2A = R2R4(atomA);
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      RA[idcoord] = Geometry(idAtm + 1,idcoord + 1);
    }
    eABC[idAtm] = 0.0;
    for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
      atomB = atoms[idBtm];
      sqrtZr4r2B = R2R4(atomB);
      RAB2 = 0.0;
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        RB[idcoord] = Geometry(idBtm + 1,idcoord + 1);
        RAB2 += (RB[idcoord] - RA[idcoord])*(RB[idcoord] - RA[idcoord]);
      }
      RAB2 *= dist_Angstrom2aum1*dist_Angstrom2aum1;           //convert to Bohr
      if (RAB2 > d2threshold) {continue;}
      c6AB = c6mat(idAtm + 1,idBtm + 1);
      cAB = a1*sqrt(3.0*sqrtZr4r2A*sqrtZr4r2B) + a2;
      for (size_t idCtm = 0; idCtm < idBtm; ++idCtm) {
        atomC = atoms[idCtm];
        sqrtZr4r2C = R2R4(atomC);
        RAC2 = 0.0;
        RBC2 = 0.0;
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          RAC2 += (Geometry(idCtm + 1,idcoord + 1) - RA[idcoord])*(Geometry(idCtm + 1,idcoord + 1) - RA[idcoord]);
          RBC2 += (Geometry(idCtm + 1,idcoord + 1) - RB[idcoord])*(Geometry(idCtm + 1,idcoord + 1) - RB[idcoord]);
        }
        RAC2 *= dist_Angstrom2aum1*dist_Angstrom2aum1;         //convert to Bohr
        RBC2 *= dist_Angstrom2aum1*dist_Angstrom2aum1;         //convert to Bohr
        if (RAC2 > d2threshold) {continue;}
        if (RBC2 > d2threshold) {continue;}
        c6AC = c6mat(idAtm + 1,idCtm + 1);
        c6BC = c6mat(idBtm + 1,idCtm + 1);
        cAC = a1*sqrt(3.0*sqrtZr4r2A*sqrtZr4r2C) + a2;
        cBC = a1*sqrt(3.0*sqrtZr4r2B*sqrtZr4r2C) + a2;
        dE = EABC(c6AB,c6AC,c6BC,cAB,cAC,cBC,RAB2,RAC2,RBC2,alpha);
        eABC[idAtm] += dE*S9;
        eABC[idBtm] += dE*S9;
        eABC[idCtm] += dE*S9;
      }
    }
  }
}
//functions related to the D4 dispersion correction
double ZetaD4(double beta1, double gammaA, double ZAref, double ZA) {
  //calculation of the zeta function, equation 2 of the D4 paper
  double zeta = 0.0;
  if (ZA < 0.0) {zeta = exp(beta1);}
  else {zeta = exp(beta1*(1.0 - exp(gammaA*(1.0 - ZAref/ZA))));}
  return zeta;
}
void ZetaD4Vector(std::vector<double> & zvec, std::vector<double> & gweight, std::vector<size_t> & atoms, std::vector<double> & QAtoms, size_t reference) {
  //calculation of the zeta vector, only for testing purposes
  int Natoms = atoms.size();
  int atomA;
  int pos = 0;
  int refA;
  double ZA;
  double QA;
  double ga = GA();
  double gc = GC();
  double hardness;
  std::vector<double> auxv(7,0.0);
  zvec.resize(gweight.size());
  for (size_t iAtm = 0; iAtm < Natoms; ++iAtm) {
    atomA = atoms[iAtm];
    refA = RefN(atomA);
    RefQ(atomA,auxv,reference);
    hardness = HardnessD4(atomA)*gc;
    ZA = Zeff(atomA);
    QA = QAtoms[iAtm] + ZA;
    for (size_t idrefA = 0; idrefA < refA; ++idrefA,++pos) {
      zvec[pos] = ZetaD4(ga,hardness,auxv[idrefA] + ZA,QA)*gweight[pos];
    }
  }
}
double dZetaD4dq(double beta1, double gammaA, double ZAref, double ZA) {
  //calculation of the zeta derivative with respect to the charge (not directly in the paper)
  double dzeta = 0.0;
  if (ZA > 0.0) {dzeta = - beta1*gammaA*exp(gammaA*(1.0 - ZAref/ZA))*ZetaD4(beta1,gammaA,ZAref,ZA)*ZAref/(ZA*ZA);}
  return dzeta;
}
void getAtomicC6(matrixE & c6coeff, std::vector<size_t> & atoms, std::vector<double> & zetavec) {
  //function getting the atomic C6 coefficients
  int Natoms = atoms.size();
  int atomA;
  int atomB;
  int refA;
  int refB;
  int posA = 0;
  int posB = 0;
  double c6;
  for (size_t iAtm = 0; iAtm < Natoms; ++iAtm) {
    atomA = atoms[iAtm];
    refA = RefN(atomA);
    posB = 0;
    for (size_t iBtm = 0; iBtm < iAtm + 1; ++iBtm) {
      atomB = atoms[iBtm];
      refB = RefN(atomB);
      c6 = 0.0;
      for (size_t idrefA = 0; idrefA < refA; ++idrefA) {
        for (size_t idrefB = 0; idrefB < refB; ++idrefB) {
          c6 += zetavec[posA + idrefA]*c6coeff(posA + idrefA + 1,posB + idrefB + 1)*zetavec[posB + idrefB];
        }
      }
      c6coeff(iAtm + 1,iBtm + 1) = c6;
      c6coeff(iBtm + 1,iAtm + 1) = c6;
      posB += refB;
    }
    posA += refA;
  }
  c6coeff.resize(Natoms,Natoms);
}
void get_gAtomicC6(matrixE & refc6, matrixE & c6, matrixE & dc6dcn, std::vector<size_t> & atoms, matrixE & zvec, matrixE & dzdcn, matrixE & dzdq) {
  //function getting the atomic C6 coefficients
  int Natoms = atoms.size();
  int atomA;
  int atomB;
  int refA;
  int refB;
  int posA = 0;
  int posB = 0;
  double rc6;
  double dc6;
  double dc6dcnA;
  double dc6dcnB;
  for (size_t iAtm = 0; iAtm < Natoms; ++iAtm) {
    atomA = atoms[iAtm];
    refA = RefN(atomA);
    posB = 0;
    for (size_t iBtm = 0; iBtm < iAtm + 1; ++iBtm) {
      atomB = atoms[iBtm];
      refB = RefN(atomB);
      dc6 = 0.0;
      dc6dcnA = 0.0;
      dc6dcnB = 0.0;
      for (size_t idrefA = 0; idrefA < refA; ++idrefA) {
        for (size_t idrefB = 0; idrefB < refB; ++idrefB) {
          rc6 = refc6(posA + idrefA + 1,posB + idrefB + 1);
          dc6 += zvec(idrefA + 1,iAtm + 1)*rc6*zvec(idrefB + 1,iBtm + 1);
          dc6dcnA += dzdcn(idrefA + 1,iAtm + 1)*rc6*zvec(idrefB + 1,iBtm + 1);
          dc6dcnB += zvec(idrefA + 1,iAtm + 1)*rc6*dzdcn(idrefB + 1,iBtm + 1);
        }
      }
      c6(iAtm + 1,iBtm + 1) = dc6;
      c6(iBtm + 1,iAtm + 1) = dc6;
      dc6dcn(iAtm + 1,iBtm + 1) = dc6dcnA;
      dc6dcn(iBtm + 1,iAtm + 1) = dc6dcnB;
      posB += refB;
    }
    posA += refA;
  }
}
void initD4(std::vector<size_t> & atoms, std::vector<int> & refA, matrixE & c6coeff, std::vector<double> & gaussweight, tensor<matrixE,double> & alphapol, std::vector<double> & ncoord, double thresholdzero = 1.0e-20) {
  //gets the polarizabilities and integrates to the dispersion coefficients; gets the Gaussian weights
  //this function goes over the calculation path; parts of the calculation effort may be replaced by D4C6par.hpp
  //currently only compatible with GFN2-xTB
  size_t Natoms = atoms.size();
  int AtomA;
  int AtomB;
  int irow;
  int icol;
  int iz;
  int icn;
  int kindex = 0;
  std::vector<int> refsys(7);
  std::vector<int> cncount(19);
  std::vector<double> refcn(7);
  std::vector<double> refcovcn(7);
  std::vector<double> refh(7);
  std::vector<double> hcount(7);
  std::vector<double> secaiw(23);
  std::vector<double> alphaiw(23);
  std::vector<double> freq(23,1.0);
  matrix<int> NCount(7,Natoms);
  alphapol.resize(23,7,Natoms);
  double sec_al;
  double hardness;
  double ga = GA();
  double gc = GC();
  double sscale;
  double ascale;
  double secq;
  double zetaf;
  double wf = WF();
  double maxrefcovcn;
  bool calculated;
  for (size_t iatm = 0; iatm < Natoms; ++iatm) {
    AtomA = atoms[iatm];
    calculated = false;
    for (size_t ibtm = 0; ibtm < iatm; ++ibtm) {
      if (atoms[ibtm] == AtomA) {
        calculated = true;
        //copy data
        for (size_t idref = 0; idref < refA[iatm]; ++idref) {
          for (size_t idpos = 0; idpos < 23; ++idpos) {
            alphapol(idpos + 1,idref + 1,iatm + 1) = alphapol(idpos + 1,idref + 1,ibtm + 1);
          }
          NCount(idref + 1,iatm + 1) = NCount(idref + 1,ibtm + 1);
        }
        break;
      }
    }
    if (calculated) {continue;}
    RefSystem(AtomA,refsys);
    RefH(AtomA,refh);
    RefCovCN(AtomA,refcovcn);
    RefCN(AtomA,refcn);
    Hcount(AtomA,hcount);
    //get polarizabilities
    for (size_t idx = 0; idx < 19; ++idx) {
      cncount[idx] = 0;
    }
    cncount[0] = 1;
    for (size_t idref = 0; idref < refA[iatm]; ++idref) {
      iz = Zeff(refsys[idref]);
      sscale = SScale(refsys[idref]);
      SecaWeights(refsys[idref],secaiw);
      secq = SecQ(refsys[idref]);
      icn = round(refcn[idref]);
      ++cncount[icn];
      RefSystem(AtomA,idref + 1,alphaiw);
      ascale = AScale(AtomA,idref + 1);
      hardness = HardnessRefSys(refsys[idref]);
      zetaf = ZetaD4(ga,hardness*gc,secq + iz,refh[idref] + iz);
      for (size_t idx = 0; idx < 23; ++idx) {
        sec_al = sscale*secaiw[idx]*zetaf;
        alphapol(idx + 1,idref + 1,iatm + 1) = fmax(ascale*(alphaiw[idx] - hcount[idref]*sec_al),0.0);       //equation 5 of D4 paper
      }
    }
    for (size_t idref = 0; idref < refA[iatm]; ++idref) {
      icn = cncount[round(refcn[idref])];
      NCount(idref + 1,iatm + 1) = icn*(icn + 1)/2;
    }
  }
  //get integration step for trapezoidal rule
  freq[0] = 0.0249995;
  freq[1] = 0.0499995;
  freq[2] = 0.075;
  for (size_t idx = 3; idx < 11; ++idx) {
    freq[idx] = 0.1;
  }
  freq[11] = 0.15;
  for (size_t idx = 12; idx < 16; ++idx) {
    freq[idx] = 0.2;
  }
  freq[16] = 0.35;
  freq[17] = 0.5;
  freq[18] = 0.75;
  freq[19] = 1.0;
  freq[20] = 1.75;
  freq[21] = 2.5;
  freq[22] = 1.25;
  //integrate C6 coefficients in accordance to equation 9 of the D4 paper
  ga = 3.0/pi;
  irow = 1;
  for (size_t iatm = 0; iatm < Natoms; ++iatm) {
    AtomA = atoms[iatm];
    RefCovCN(AtomA,refcovcn);
    hardness = 0.0;                  //normalization constant for the Gaussian weights
    maxrefcovcn = 0.0;
    for (size_t idref = 0; idref < refA[iatm]; ++idref) {
      icn = NCount(idref + 1,iatm + 1);
      gaussweight[kindex + idref] = 0.0;
      if (refcovcn[idref] > maxrefcovcn) {maxrefcovcn = refcovcn[idref];}
      for (int idgauss = 0; idgauss < icn; ++idgauss) {
        ascale = wf*double(idgauss + 1);
        gc = Gauss(ascale,ncoord[iatm] - refcovcn[idref]);
        hardness += gc;
        gaussweight[kindex + idref] += gc;                 //equation 8 of D4 paper without the normalization constant
      }
    }
    iz = 0;                          //offset to gaussweight
    if (fabs(hardness) < thresholdzero) {
      for (size_t idref = 0; idref < refA[iatm]; ++idref,++iz) {
        gaussweight[kindex + idref] = double(maxrefcovcn == refcovcn[idref]);
      }
    }
    else {
      gc = 1.0/hardness;
      for (size_t idref = 0; idref < refA[iatm]; ++idref,++iz) {
        gaussweight[kindex + idref] *= gc;
      }
    }
    kindex += iz;
    icol = 1;
    for (size_t ibtm = 0; ibtm < iatm + 1; ++ibtm) {
      AtomB = atoms[ibtm];
      for (size_t irefA = 0; irefA < refA[iatm]; ++irefA) {
        for (size_t irefB = 0; irefB < refA[ibtm]; ++irefB) {
          for (size_t idx = 0; idx < 23; ++idx) {
            alphaiw[idx] = alphapol(idx + 1,irefA + 1,iatm + 1)*alphapol(idx + 1,irefB + 1,ibtm + 1);   //product of polarisabilities
          }
          sscale = ga*TrapezoidalIntegration(alphaiw,freq);                     //c6, equation 9 of D4 paper
          c6coeff(irow + irefA,icol + irefB) = sscale;
          c6coeff(icol + irefB,irow + irefA) = sscale;
        }
      }
      icol += refA[ibtm];
    }
    irow += refA[iatm];
  }
}
void getC6D4(const std::string & method, matrixE & c6coeff, std::vector<double> & gweight, std::vector<double> & ncoord, std::vector<int> & refA, std::vector<size_t> & atoms, bool extract = true, double thresholdzero = 1.0e-20) {
  //function that gets the full matrix with all the required C6 coefficients
  //if extract is true, then I just access the extracted parameters in the file D4C6par, otherwise calculation via initD4
  size_t Natoms = atoms.size();
  if (extract) {
    int AtomA;
    int AtomB;
    int kindex = 0;
    int kadd;
    int irow = 1;
    int icol = 1;
    double aux;
    double norm;
    double twf;
    double wf = WF();
    double maxrefcovcn;
    std::vector<double> refcovcn(7);
    matrixE c6aux(7,7);
    std::vector<int> ncount(7);
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      AtomA = atoms[iatm];
      getNcount(AtomA,ncount);
      norm = 0.0;
      maxrefcovcn = 0.0;
      for (size_t idref = 0; idref < refA[iatm]; ++idref) {
        RefCovCN(AtomA,refcovcn);
        gweight[kindex + idref] = 0.0;
        if (refcovcn[idref] > maxrefcovcn) {maxrefcovcn = refcovcn[idref];}
        for (int idgauss = 0; idgauss < ncount[idref]; ++idgauss) {
          twf = wf*double(idgauss + 1);
          aux = Gauss(twf,ncoord[iatm] - refcovcn[idref]);
          norm += aux;
          gweight[kindex + idref] += aux;                 //equation 8 of D4 paper without the normalization constant
        }
      }
      kadd = 0;
      if (fabs(norm) < thresholdzero) {
        for (size_t idref = 0; idref < refA[iatm]; ++idref,++kadd) {
          gweight[kindex + idref] = double(maxrefcovcn == refcovcn[idref]);
        }
      }
      else {
        aux = 1.0/norm;
        for (size_t idref = 0; idref < refA[iatm]; ++idref,++kadd) {
          gweight[kindex + idref] *= aux;
        }
      }
      kindex += kadd;
      icol = 1;
      for (size_t ibtm = 0; ibtm < iatm + 1; ++ibtm) {
        AtomB = atoms[ibtm];
        C6parametersD4(method,c6aux,AtomA,AtomB);
        if (AtomA < AtomB) {c6aux = c6aux.trans();}
        for (size_t idrefA = 0; idrefA < refA[iatm]; ++idrefA) {
          for (size_t idrefB = 0; idrefB < refA[ibtm]; ++idrefB) {
            aux = c6aux(1 + idrefA,1 + idrefB);
            c6coeff(irow + idrefA,icol + idrefB) = aux;
            c6coeff(icol + idrefB,irow + idrefA) = aux;
          }
        }
        icol += refA[ibtm];
      }
      irow += refA[iatm];
    }
  }
  else {
    tensor<matrixE,double> alphapol(1,1,1);                     //this is only needed to call the function
    initD4(atoms,refA,c6coeff,gweight,alphapol,ncoord);
  }
}
void getWeightedGFN2(matrixE & wd4, matrixE & c6coeff, std::vector<double> & gweight, std::vector<size_t> & atoms, matrixE & Geometry, std::vector<int> & refA, std::vector<double> & ncoord, bool extract = true, double thresholdzero = 1.0e-20, double d2threshold = 1600.0, double gwcutoff = 1.0e-7) {
  //calculating the weighted dispersion matrix
  int AtomA;
  int AtomB;
  int Natoms = atoms.size();
  int irow = 1;
  int icol = 1;
  int szc6coeff = 0;
  for (size_t iatm = 0; iatm < Natoms; ++iatm) {             //get size of C6coeff and the number of reference points for each atom
    AtomA = atoms[iatm];
    refA[iatm] = RefN(AtomA);                                //the number of reference Gaussians is stored here because I need the loop anyway; functions called within getWeightedD4 do not have to recall function RefN
    szc6coeff += refA[iatm];
  }
  double sqrtZr4r2A;
  double r4r2AB;
  double c10conv = 1.225;                                    // = 49.0/40.0
  double a1 = A1D4("GFN2");
  double a2 = A2D4("GFN2");
  double cutoff;
  double rAB2;
  double aux1;
  double aux2;
  double oor6;
  double oor8;
  double oor10;
  double S6 = s6D4("GFN2");
  double S8 = s8D4("GFN2");
  c10conv *= s10D4("GFN2");
  double rA[3];
  c6coeff.resize(szc6coeff,szc6coeff);
  gweight.resize(szc6coeff);
  NCoord(ncoord,atoms,Geometry,ERFcount,true,1600.0,false);  //get coordination numbers
  getC6D4("GFN2",c6coeff,gweight,ncoord,refA,atoms,extract,thresholdzero);
  wd4.resize(szc6coeff,szc6coeff);
  wd4.zero();
  for (size_t idatm = 0; idatm < Natoms; ++idatm) {
    AtomA = atoms[idatm];
    sqrtZr4r2A = R2R4(AtomA);
    rA[0] = Geometry(idatm + 1,1);
    rA[1] = Geometry(idatm + 1,2);
    rA[2] = Geometry(idatm + 1,3);
    icol = 1;
    for (size_t idbtm = 0; idbtm < idatm; ++idbtm) {
      AtomB = atoms[idbtm];
      rAB2 = 0.0;
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        aux1 = Geometry(idbtm + 1,idcoord + 1) - rA[idcoord];
        rAB2 += aux1*aux1*dist_Angstrom2aum1*dist_Angstrom2aum1;
      }
      if (rAB2 <= d2threshold) {
        r4r2AB = 3.0*sqrtZr4r2A*R2R4(AtomB);                 //equation 20 of D4 paper
        cutoff = a1*sqrt(r4r2AB) + a2;                       //equation 21 of D4 paper
        aux1 = rAB2*rAB2*rAB2;
        aux2 = cutoff*cutoff*cutoff;
        oor6  = 1.0/(aux1 + aux2*aux2);                      //equation 19 of D4 paper
        aux1 *= rAB2;
        aux2 *= cutoff;
        oor8  = 1.0/(aux1 + aux2*aux2);                      //equation 19 of D4 paper
        aux1 *= rAB2;
        aux2 *= cutoff;
        oor10 = 1.0/(aux1 + aux2*aux2);                      //equation 19 of D4 paper
        for (size_t idrefA = 0; idrefA < refA[idatm]; ++idrefA) {
          for (size_t idrefB = 0; idrefB < refA[idbtm]; ++idrefB) {
            rAB2 = gweight[irow + idrefA - 1]*gweight[icol + idrefB - 1];
            if (rAB2 < gwcutoff) {continue;}
            aux1 = c6coeff(irow + idrefA,icol + idrefB);
            cutoff = -rAB2*(S6*oor6 + S8*r4r2AB*oor8 + c10conv*r4r2AB*r4r2AB*oor10)*aux1;   //this is a modified equation 18 of D4 paper
            wd4(irow + idrefA,icol + idrefB) = cutoff;
            wd4(icol + idrefB,irow + idrefA) = cutoff;
          }
        }
      }
      icol += refA[idbtm];
    }
    irow += refA[idatm];
  }
}
double ED4(const std::string & method, matrixE & c6coeff, std::vector<double> & gweight, std::vector<size_t> & atoms, matrixE & Geometry, std::vector<int> & refA, std::vector<double> & ncoord, std::vector<double> & QAtoms, bool extract = true, double thresholdzero = 1.0e-10, double d2threshold = 1600.0, double gwcutoff = 1.0e-7) {
  //calculating the weighted dispersion matrix
  int atomA;
  int atomB;
  int atomC;
  int Natoms = atoms.size();
  int irow = 1;
  int icol = 1;
  int szc6coeff = 0;
  for (size_t iatm = 0; iatm < Natoms; ++iatm) {             //get size of C6coeff and the number of reference points for each atom
    atomA = atoms[iatm];
    refA[iatm] = RefN(atomA);                                //the number of reference Gaussians is stored here because I need the loop anyway; functions called within getWeightedD4 do not have to recall function RefN
    szc6coeff += refA[iatm];
  }
  double c6atm;
  double c6AB;
  double c6AC;
  double c6BC;
  double cAB;
  double cAC;
  double cBC;
  double sqrtZr4r2A;
  double sqrtZr4r2B;
  double sqrtZr4r2C;
  double r4r2AB;
  double c10conv = 1.225;                                    // = 49.0/40.0
  double a1 = A1D4(method);
  double a2 = A2D4(method);
  double cutoff;
  double RAB2;
  double RAC2;
  double RBC2;
  double oor6;
  double oor8;
  double oor10;
  double alpha = AlphaD4(method,9);
  double S6 = s6D4(method);
  double S8 = s8D4(method);
  double S9 = s9D4(method);
  double edftd4 = 0.0;
  c10conv *= s10D4(method);
  double rA[3];
  double rB[3];
  std::vector<double> GWsav(szc6coeff);
  std::vector<double> GWsavATM(szc6coeff);
  c6coeff.resize(szc6coeff,szc6coeff);
  gweight.resize(szc6coeff);
  d2threshold *= dist_Angstrom2au*dist_Angstrom2au;
  getC6D4(method,c6coeff,gweight,ncoord,refA,atoms,extract,thresholdzero);
  ZetaD4Vector(GWsav,gweight,atoms,QAtoms,0);
  //In D4 the ATM terms are obtained using a zero charge model. Is this actually better?
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    ncoord[idAtm] = 0.0;
  }
  ZetaD4Vector(GWsavATM,gweight,atoms,ncoord,0);
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    atomA = atoms[idAtm];
    sqrtZr4r2A = R2R4(atomA);
    rA[0] = Geometry(idAtm + 1,1);
    rA[1] = Geometry(idAtm + 1,2);
    rA[2] = Geometry(idAtm + 1,3);
    icol = 1;
    for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
      atomB = atoms[idBtm];
      sqrtZr4r2B = R2R4(atomB);
      RAB2 = 0.0;
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        rB[idcoord] = Geometry(idBtm + 1,idcoord + 1);
        sqrtZr4r2C = (rB[idcoord] - rA[idcoord])*dist_Angstrom2aum1;
        RAB2 += sqrtZr4r2C*sqrtZr4r2C;
      }
      c6AB = 0.0;
      if (RAB2 <= d2threshold) {
        r4r2AB = 3.0*sqrtZr4r2A*sqrtZr4r2B;                  //equation 20 of D4 paper
        cAB = a1*sqrt(r4r2AB) + a2;                          //equation 21 of D4 paper
        sqrtZr4r2C = RAB2*RAB2*RAB2;
        RAC2 = cAB*cAB*cAB;
        oor6  = 1.0/(sqrtZr4r2C + RAC2*RAC2);                //equation 19 of D4 paper
        sqrtZr4r2C *= RAB2;
        RAC2 *= cAB;
        oor8  = 1.0/(sqrtZr4r2C + RAC2*RAC2);                //equation 19 of D4 paper
        sqrtZr4r2C *= RAB2;
        RAC2 *= cAB;
        oor10 = 1.0/(sqrtZr4r2C + RAC2*RAC2);                //equation 19 of D4 paper
        c6atm = 0.0;
        for (size_t idrefA = 0; idrefA < refA[idAtm]; ++idrefA) {
          for (size_t idrefB = 0; idrefB < refA[idBtm]; ++idrefB) {
            RAC2 = GWsavATM[irow + idrefA - 1]*GWsavATM[icol + idrefB - 1];
            if (RAC2 >= gwcutoff) {c6AB += c6coeff(irow + idrefA,icol + idrefB)*RAC2;}
            RAC2 = GWsav[irow + idrefA - 1]*GWsav[icol + idrefB - 1];
            if (RAC2 >= gwcutoff) {c6atm += c6coeff(irow + idrefA,icol + idrefB)*RAC2;}
          }
        }
        edftd4 -= c6atm*(S6*oor6 + S8*r4r2AB*oor8 + c10conv*r4r2AB*r4r2AB*oor10);   //this is a modified equation 18 of D4 paper
        //now the atomwise contribution to ATM three-body term
        for (size_t idCtm = 0; idCtm < idBtm; ++idCtm) {
          atomC = atoms[idCtm];
          sqrtZr4r2C = R2R4(atomC);
          RAC2 = 0.0;
          RBC2 = 0.0;
          for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
            RAC2 += (Geometry(idCtm + 1,idcoord + 1) - rA[idcoord])*(Geometry(idCtm + 1,idcoord + 1) - rA[idcoord]);
            RBC2 += (Geometry(idCtm + 1,idcoord + 1) - rB[idcoord])*(Geometry(idCtm + 1,idcoord + 1) - rB[idcoord]);
          }
          RAC2 *= dist_Angstrom2aum1*dist_Angstrom2aum1;         //convert to Bohr
          RBC2 *= dist_Angstrom2aum1*dist_Angstrom2aum1;         //convert to Bohr
          if (RAC2 > d2threshold) {continue;}
          if (RBC2 > d2threshold) {continue;}
          c6AC = c6coeff(idAtm + 1,idCtm + 1);
          c6BC = c6coeff(idBtm + 1,idCtm + 1);
          cAC = a1*sqrt(3.0*sqrtZr4r2A*sqrtZr4r2C) + a2;
          cBC = a1*sqrt(3.0*sqrtZr4r2B*sqrtZr4r2C) + a2;
          edftd4 += S9*EABC(c6AB,c6AC,c6BC,cAB,cAC,cBC,RAB2,RAC2,RBC2,alpha);
        }
      }
      c6coeff(idAtm + 1,idBtm + 1) = c6AB;
      c6coeff(idBtm + 1,idAtm + 1) = c6AB;
      icol += refA[idBtm];
    }
    irow += refA[idAtm];
  }
  return edftd4;
}
void getWeightReferences(matrixE & zvec, matrixE & dz_dcn, matrixE & dz_dq, matrixE & z0vec, matrixE & dz0_dcn, std::vector<size_t> & atoms, std::vector<double> & ncoord, std::vector<double> & QAtoms, size_t reference, double thresholdzero = 1.0e-10) {
  //get the weights of the reference systems and derivatives with respect to coordination number
  int Natoms = zvec.cols();
  int nrefs = zvec.rows();
  int atomA;
  int refA;
  double zA;
  double ga = GA();
  double gc = GC();
  double hardness;
  double norm;
  double dnorm;
  double gW;
  double gWk;
  double dgWk;
  double twf;
  double wf = WF();
  double expw;
  double expd;
  double z4;
  double z04;
  double gw;
  double aux;
  double maxrefcovcn;
  std::vector<int> ncount(7);
  std::vector<double> refcovcn(7);
  std::vector<double> QrefV(7,0.0);
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    atomA = atoms[idAtm];
    refA = RefN(atomA);
    zA = Zeff(atomA);
    hardness = HardnessD4(atomA);
    getNcount(atomA,ncount);
    RefCovCN(atomA,refcovcn);
    RefQ(atomA,QrefV,reference);
    maxrefcovcn = 0.0;
    norm = 0.0;
    dnorm = 0.0;
    for (size_t idref = 0; idref < refA; ++idref) {
      if (refcovcn[idref] > maxrefcovcn) {maxrefcovcn = refcovcn[idref];}
      for (int idgauss = 0; idgauss < ncount[idref]; ++idgauss) {
        twf = wf*double(idgauss + 1);
        aux = refcovcn[idref] - ncoord[idAtm];
        gw = Gauss(twf,aux);
        norm += gw;
        dnorm += 2.0*twf*aux*gw;
      }
    }
    for (size_t idref = 0; idref < refA; ++idref) {
      expw = 0.0;
      expd = 0.0;
      for (int idgauss = 0; idgauss < ncount[idref]; ++idgauss) {
        twf = wf*double(idgauss + 1);
        aux = refcovcn[idref] - ncoord[idAtm];
        gw = Gauss(twf,aux);
        expw += gw;
        expd += 2.0*twf*aux*gw;
      }
      if (norm < thresholdzero) {
        if (maxrefcovcn == refcovcn[idref]) {gWk = 1.0;}
        else {gWk = 0.0;}
      }
      else {gWk = expw/norm;}
      z4 = ZetaD4(ga,hardness*gc,QrefV[idref] + zA,QAtoms[idAtm] + zA);
      z04 = ZetaD4(ga,hardness*gc,QrefV[idref] + zA,zA);
      zvec(idref + 1,idAtm + 1) = z4*gWk;
      z0vec(idref + 1,idAtm + 1) = z04*gWk;
      if (norm < thresholdzero) {dgWk = 0.0;}
      else {dgWk = (expd - expw*dnorm/norm)/norm;}
      dz_dcn(idref + 1,idAtm + 1) = z4*dgWk;
      dz0_dcn(idref + 1,idAtm + 1) = z04*dgWk;
      z4 = dZetaD4dq(ga,hardness*gc,QrefV[idref] + zA,QAtoms[idAtm] + zA);
      dz_dq(idref + 1,idAtm + 1) = z4*gWk;
    }
  }
}
void DispGrad(const std::string & method, matrixE & grad, matrixE & c6AB, matrixE & dc6ABdCN, matrixE & geom, matrixE & dEdCN, std::vector<size_t> & atoms, double d2thresh = 1600.0) {
  int atomA;
  int atomB;
  int Natoms = atoms.size();
  double r4r2AB;
  double sqrtZr4r2A;
  double rAB2;
  double rAB6;
  double rAB8;
  double r0;
  double r03;
  double r04;
  double t6;
  double t8;
  double t10;
  double d6;
  double d8;
  double d10;
  double a1 = A1D4(method);
  double a2 = A2D4(method);
  double S6 = s6D4(method);
  double S8 = s8D4(method);
  double S10 = 1.225*s10D4(method);
  double disp;
  double ddisp;
  double rAB[3];
  double rA[3];
  double dG;
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    atomA = atoms[idAtm];
    sqrtZr4r2A = R2R4(atomA);
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      rA[idcoord] = geom(idAtm + 1,idcoord + 1)*dist_Angstrom2aum1;
    }
    dEdCN(1,idAtm + 1) = 0.0;
    for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
      atomB = atoms[idBtm];
      rAB2 = 0.0;
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        rAB[idcoord] = rA[idcoord] - geom(idBtm + 1,idcoord + 1)*dist_Angstrom2aum1;
        rAB2 += rAB[idcoord]*rAB[idcoord];
      }
      if (rAB2 > d2thresh) {continue;}
      r4r2AB = 3.0*sqrtZr4r2A*R2R4(atomB);
      r0 = a1*sqrt(r4r2AB) + a2;
      r03 = r0*r0*r0;
      r04 = r03*r0;
      rAB6 = rAB2*rAB2*rAB2;
      rAB8 = rAB6*rAB2;
      t6 = 1.0/(rAB6 + r03*r03);
      t8 = 1.0/(rAB8 + r04*r04);
      t10 = 1.0/(rAB8*rAB2 + r04*r04*r0*r0);
      d6 = -6.0*rAB2*rAB2*t6*t6;
      d8 = -8.0*rAB6*t8*t8;
      d10 = 10.0*rAB8*t10*t10;
      disp = S6*t6 + S8*r4r2AB*t8 + S10*r4r2AB*r4r2AB*t10;
      ddisp = S6*d6 + S8*r4r2AB*d8 + S10*r4r2AB*r4r2AB*d10;
      dEdCN(1,idAtm + 1) -= dc6ABdCN(idAtm + 1, idBtm + 1)*disp;
      dEdCN(1,idBtm + 1) -= dc6ABdCN(idBtm + 1, idAtm + 1)*disp;
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        dG = -c6AB(idAtm + 1,idBtm + 1)*ddisp*rAB[idcoord];
        grad(3*idAtm + idcoord + 1,1) += dG;
        grad(3*idBtm + idcoord + 1,1) -= dG;
      }
    }
  }
}
double ATMscale(int ii, int jj, int kk) {
  //function returning the scaling factor for the Axilrod-Teller-Muto term
  double scale = 0.0;
  if (ii == jj) {
    if (ii == kk) {scale = 1.0/6.0;}                    //i i' i''
    else {scale = 0.5;}                                 //i i' k
  }
  else {
    if ((ii != kk)&&(jj != kk)) {scale = 1.0;}          //i j k
    else {scale = 0.5;}                                 //i j i' or i j j'
  }
  return scale;
}
void dATM(const std::string & method, matrixE & gradient, matrixE & c6AB, matrixE & dc6dCN, matrixE & dEdCN, std::vector<size_t> & atoms, matrixE & geometry, matrixE & dG, double d2thresh = 1600.0) {
  //function to calculate the derivatives of the Axilrod-Teller-Muto term
  int Natoms = atoms.size();
  int atomA;
  int atomB;
  int atomC;
  double c6ab;
  double c6ac;
  double c6bc;
  double c9;
  double a1 = A1D4(method);
  double a2 = A2D4(method);
  double alpha = AlphaD4(method,9);
  double cAB;
  double cAC;
  double cBC;
  double sqrtZr4r2A;
  double sqrtZr4r2B;
  double sqrtZr4r2C;
  double rAB2;
  double rAC2;
  double rBC2;
  double ccc;
  double rrr1;
  double rrr2;
  double rrr3;
  double angle;
  double dangleAB;
  double dangleAC;
  double dangleBC;
  double cubicroot;
  double fdamp;
  double dfdamp;
  double aux;
  double dGrAB;
  double dGrAC;
  double dGrBC;
  double rAB[3];
  double rAC[3];
  double rBC[3];
  double S9 = s9D4(method);
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    atomA = atoms[idAtm];
    sqrtZr4r2A = R2R4(atomA);
    for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
      atomB = atoms[idBtm];
      sqrtZr4r2B = R2R4(atomB);
      c6ab = c6AB(idAtm + 1,idBtm + 1);
      cAB = a1*sqrt(3.0*sqrtZr4r2A*sqrtZr4r2B) + a2;
      rAB2 = 0.0;
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        rAB[idcoord] = (geometry(idBtm + 1,idcoord + 1) - geometry(idAtm + 1,idcoord + 1))*dist_Angstrom2aum1;
        rAB2 += rAB[idcoord]*rAB[idcoord];
      }
      if (rAB2 > d2thresh) {continue;}
      for (size_t idCtm = 0; idCtm < idBtm; ++idCtm) {
        atomC = atoms[idCtm];
        sqrtZr4r2C = R2R4(atomC);
        c6ac = c6AB(idAtm + 1,idCtm + 1);
        c6bc = c6AB(idBtm + 1,idCtm + 1);
        cAC = a1*sqrt(3.0*sqrtZr4r2A*sqrtZr4r2C) + a2;
        cBC = a1*sqrt(3.0*sqrtZr4r2B*sqrtZr4r2C) + a2;
        rAC2 = 0.0;
        rBC2 = 0.0;
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          rAC[idcoord] = (geometry(idCtm + 1,idcoord + 1) - geometry(idAtm + 1,idcoord + 1))*dist_Angstrom2aum1;
          rBC[idcoord] = (geometry(idCtm + 1,idcoord + 1) - geometry(idBtm + 1,idcoord + 1))*dist_Angstrom2aum1;
          rAC2 += rAC[idcoord]*rAC[idcoord];
          rBC2 += rBC[idcoord]*rBC[idcoord];
        }
        if (rAC2 > d2thresh) {continue;}
        else if (rBC2 > d2thresh) {continue;}
        c9 = -sqrt(c6ab*c6ac*c6bc);
        ccc = cAB*cAC*cBC;
        rrr2 = rAB2*rAC2*rBC2;
        rrr1 = sqrt(rrr2);
        rrr3 = rrr1*rrr2;
        angle = 0.375*(rAB2 + rBC2 - rAC2)*(rAB2 - rBC2 + rAC2)*(-rAB2 + rBC2 + rAC2)/(rrr3*rrr2) + 1.0/rrr3;
        cubicroot = cbrt(ccc/rrr1);
        aux = pow(cubicroot,alpha);
        fdamp = 1.0/(1.0 + 6.0*aux);
        dfdamp = -2.0*alpha*aux*fdamp*fdamp;
        dangleAB = -0.375*(rAB2*rAB2*rAB2 + rAB2*rAB2*(rBC2 + rAC2) + rAB2*(3.0*rBC2*rBC2 + 2.0*rBC2*rAC2 + 3.0*rAC2*rAC2) -5.0*(rBC2 - rAC2)*(rBC2 - rAC2)*(rBC2 + rAC2))/(rrr3*rrr2);
        dGrAB = c9*(dfdamp*angle - dangleAB*fdamp)/rAB2;
        dangleAC = -0.375*(rAC2*rAC2*rAC2 + rAC2*rAC2*(rBC2 + rAB2) + rAC2*(3.0*rBC2*rBC2 + 2.0*rBC2*rAB2 + 3.0*rAB2*rAB2) - 5.0*(rBC2 - rAB2)*(rBC2 - rAB2)*(rBC2 + rAB2))/(rrr3*rrr2);
        dangleBC = -0.375*(rBC2*rBC2*rBC2 + rBC2*rBC2*(rAC2 + rAB2) + rBC2*(3.0*rAC2*rAC2 + 2.0*rAC2*rAB2 + 3.0*rAB2*rAB2) - 5.0*(rAC2 - rAB2)*(rAC2 - rAB2)*(rAC2 + rAB2))/(rrr3*rrr2);
        dGrAC = c9*(dfdamp*angle - dangleAC*fdamp)/rAC2;
        dGrBC = c9*(dfdamp*angle - dangleBC*fdamp)/rBC2;
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          dG(idcoord + 1,1) = -dGrAB*rAB[idcoord] - dGrAC*rAC[idcoord];
          dG(idcoord + 1,2) = dGrAB*rAB[idcoord] - dGrBC*rBC[idcoord];
          dG(idcoord + 1,3) = dGrAC*rAC[idcoord] + dGrBC*rBC[idcoord];
        }
        dangleAB = 0.5*c9*(dc6dCN(idAtm + 1,idBtm + 1)/c6ab + dc6dCN(idAtm + 1,idCtm + 1)/c6ac);
        dangleAC = 0.5*c9*(dc6dCN(idBtm + 1,idAtm + 1)/c6ab + dc6dCN(idBtm + 1,idCtm + 1)/c6bc);
        dangleBC = 0.5*c9*(dc6dCN(idCtm + 1,idAtm + 1)/c6ac + dc6dCN(idCtm + 1,idBtm + 1)/c6bc);
        aux = S9*ATMscale(idAtm + 1,idBtm + 1,idCtm + 1);
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          gradient(3*idAtm + idcoord + 1,1) += dG(idcoord + 1,1)*aux;
          gradient(3*idBtm + idcoord + 1,1) += dG(idcoord + 1,2)*aux;
          gradient(3*idCtm + idcoord + 1,1) += dG(idcoord + 1,3)*aux;
        }
        dEdCN(1,idAtm + 1) -= angle*fdamp*dangleAB*aux;
        dEdCN(1,idBtm + 1) -= angle*fdamp*dangleAC*aux;
        dEdCN(1,idCtm + 1) -= angle*fdamp*dangleBC*aux;
      }
    }
  }
}
void dED4GFN2(const std::string & method, matrixE & grad, std::vector<matrixE> & dCN, matrixE & refc6, std::vector<size_t> & atoms, matrixE & geometry, std::vector<double> & ncoord, std::vector<double> & QAtoms, int maxREFn, bool doATM, double d2thresh = 1600.0, double d3thresh = 1600.0) {
  //gradients of the D4 correction
  int Natoms = atoms.size();
  matrixE zetavec(maxREFn,Natoms);
  matrixE dzeta_dcn(maxREFn,Natoms);
  matrixE dzeta_dq(maxREFn,Natoms);
  matrixE zerovec(maxREFn,Natoms);
  matrixE dzero_dcn(maxREFn,Natoms);
  matrixE c6(Natoms,Natoms);
  matrixE dc6_dcn(Natoms,Natoms);
  getWeightReferences(zetavec,dzeta_dcn,dzeta_dq,zerovec,dzero_dcn,atoms,ncoord,QAtoms,1,1.0e-10);
  get_gAtomicC6(refc6,c6,dc6_dcn,atoms,zetavec,dzeta_dcn,dzeta_dq);
  DispGrad(method,grad,c6,dc6_dcn,geometry,zetavec,atoms,d2thresh);
  if (doATM) {
    get_gAtomicC6(refc6,c6,dc6_dcn,atoms,zerovec,dzero_dcn,dzeta_dq);
    if (maxREFn < 3) {dzero_dcn.resize(3,Natoms);}
    dATM(method,grad,c6,dc6_dcn,zetavec,atoms,geometry,dzero_dcn,d3thresh);
  }
  dzeta_dcn = dCN[0]*zetavec.trans();
  zerovec = dCN[1]*zetavec.trans();
  dzeta_dq = dCN[2]*zetavec.trans();
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    grad(3*idAtm + 1,1) += dzeta_dcn(idAtm + 1,1);
    grad(3*idAtm + 2,1) += zerovec(idAtm + 1,1);
    grad(3*idAtm + 3,1) += dzeta_dq(idAtm + 1,1);
  }
}
int dCNdR(std::vector<matrixE> & dCN, matrixE & geometry, std::vector<size_t> & atoms, std::vector<double> & rA, std::vector<double> & rB, double (*dcount)(double, double, double), bool useden = false, double d2threshold = 1600.0, bool ERFfunctionuse = false) {
  //function calculating the derivative of the coordination number with respect to interatomic distances
  int atomA;
  int atomB;
  int Natoms = atoms.size();
  double dcn;
  double distAB;
  int maxrefN = 0;
  double den;
  double aux1;
  double ENA;
  double k4 = k4D4();
  double k5 = k5D4();
  double k6 = 1.0/k6D4();
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    atomA = atoms[idAtm];
    ENA = PaulingElectronegativity(atoms[idAtm]);
    if (RefN(atomA) > maxrefN) {maxrefN = RefN(atomA);}
    for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
      rA[idcoord] = geometry(idAtm + 1,idcoord + 1);
    }
    for (size_t idBtm = 0; idBtm < idAtm; ++idBtm) {
      atomB = atoms[idBtm];
      dcn = 0.0;
      for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
        rB[idcoord] = geometry(idBtm + 1,idcoord + 1) - rA[idcoord];
        dcn += rB[idcoord]*rB[idcoord];
      }
      if (dcn < d2threshold) {
        distAB = sqrt(dcn);
        den = 1.0;
        if (useden) {
          aux1 = fabs(ENA - PaulingElectronegativity(atoms[idBtm])) + k5;
          den = k4*exp(-aux1*aux1*k6);
        }
        dcn = den*dNCoorddR(atomA,atomB,distAB,(*dcount),useden,d2threshold,ERFfunctionuse);
        dcn /= distAB;
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          rB[idcoord] *= -dcn;
        }
        for (size_t idcoord = 0; idcoord < 3; ++idcoord) {
          dCN[idcoord](idAtm + 1,idAtm + 1) += rB[idcoord];
          dCN[idcoord](idAtm + 1,idBtm + 1) += rB[idcoord];
          dCN[idcoord](idBtm + 1,idAtm + 1) -= rB[idcoord];
          dCN[idcoord](idBtm + 1,idBtm + 1) -= rB[idcoord];
        }
      }
    }
  }
  return maxrefN;
}
void MolecularDispersionData(double & mc6, double & mc8, double & mpol, std::vector<size_t> & atoms, matrixE & Geometry, const std::vector<double> & QAtoms, size_t reference) {
  //calculating dispersion properties at the molecular level
  //mc6 -> the C6 coefficients in au.(bohr)^6
  //mc8 -> the C8 coefficients in au.(bohr)^8
  //mpol -> the molecular polarizability in au
  int Natoms = atoms.size();
  int atomA;
  int iposA = 0;
  int iposB = 0;
  size_t szc6coeff = 0;
  double zA;
  double sqrtZr4r2A;
  double ga = GA();
  double gc = GC();
  double hardness;
  double c6AB;
  std::vector<int> refA(Natoms);
  for (size_t iatm = 0; iatm < Natoms; ++iatm) {             //get size of C6coeff and the number of reference points for each atom
    atomA = atoms[iatm];
    refA[iatm] = RefN(atomA);                                //the number of reference Gaussians is stored here because I need the loop anyway; functions called within getWeightedD4 do not have to recall function RefN
    szc6coeff += refA[iatm];
  }
  std::vector<double> gweight(szc6coeff,0.0);
  std::vector<double> zetaA(szc6coeff,0.0);
  std::vector<double> auxv(7,0.0);
  std::vector<double> ncoord;
  NCoord(ncoord,atoms,Geometry,ERFcount,true,1600.0,false);  //get coordination numbers
  matrixE pol(23,Natoms);
  matrixE d6(szc6coeff,szc6coeff);
  tensor<matrixE,double> alphapol(1,1,1);
  initD4(atoms,refA,d6,gweight,alphapol,ncoord);             //get polarizabilities, Gaussian weights and the dispersion coefficients
  mc6 = 0.0;
  mc8 = 0.0;
  mpol = 0.0;
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    atomA = atoms[idAtm];
    zA = Zeff(atomA);
    hardness = HardnessD4(atomA);
    RefQ(atomA,auxv,reference);
    sqrtZr4r2A = R2R4(atomA);
    for (size_t idrefA = 0; idrefA < refA[idAtm]; ++idrefA) {
      gweight[iposA + idrefA] *= ZetaD4(ga,hardness*gc,auxv[idrefA] + zA,QAtoms[idAtm] + zA);
      for (size_t idx = 0; idx < 23; ++idx) {
        pol(idx + 1,idAtm + 1) += gweight[iposA + idrefA]*alphapol(idx + 1,idrefA + 1,idAtm + 1);
      }
    }
    mpol += pol(1,idAtm + 1);
    iposB = 0;
    for (size_t idBtm = 0; idBtm < idAtm + 1; ++idBtm) {
      c6AB = 0.0;
      for (size_t idrefA = 0; idrefA < refA[idAtm]; ++idrefA) {
        for (size_t idrefB = 0; idrefB < refA[idBtm]; ++idrefB) {
          c6AB += gweight[iposA + idrefA]*d6(iposA + idrefA + 1,iposB + idrefB + 1)*gweight[iposB + idrefB];
        }
      }
      mc6 += c6AB*(1.0 + double(idAtm != idBtm));
      mc8 += 3.0*sqrtZr4r2A*c6AB*R2R4(atoms[idBtm])*(1.0 + double(idAtm != idBtm));
      iposB += refA[idBtm];
    }
    iposA += refA[idAtm];
  }
}
void Polarizabilities(std::vector<double> & polarizabilities, std::vector<size_t> & atoms, matrixE & Geometry, const std::vector<double> & QAtoms, size_t reference) {
  //calculating dispersion properties at the molecular level
  //mc6 -> the C6 coefficients in au.(bohr)^6
  //mc8 -> the C8 coefficients in au.(bohr)^8
  //mpol -> the molecular polarizability in au
  int Natoms = atoms.size();
  if (polarizabilities.size() != Natoms) {polarizabilities.resize(Natoms);}
  int atomA;
  int iposA = 0;
  size_t szc6coeff = 0;
  double zA;
  double sqrtZr4r2A;
  double ga = GA();
  double gc = GC();
  double hardness;
  std::vector<int> refA(Natoms);
  for (size_t iatm = 0; iatm < Natoms; ++iatm) {             //get size of C6coeff and the number of reference points for each atom
    atomA = atoms[iatm];
    refA[iatm] = RefN(atomA);                                //the number of reference Gaussians is stored here because I need the loop anyway; functions called within getWeightedD4 do not have to recall function RefN
    szc6coeff += refA[iatm];
  }
  std::vector<double> gweight(szc6coeff,0.0);
  std::vector<double> auxv(7,0.0);
  std::vector<double> ncoord;
  NCoord(ncoord,atoms,Geometry,ERFcount,true,1600.0,false);  //get coordination numbers
  matrixE d6(szc6coeff,szc6coeff);
  tensor<matrixE,double> alphapol(1,1,1);
  initD4(atoms,refA,d6,gweight,alphapol,ncoord);             //get polarizabilities, Gaussian weights and the dispersion coefficients
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    atomA = atoms[idAtm];
    zA = Zeff(atomA);
    hardness = HardnessD4(atomA);
    RefQ(atomA,auxv,reference);
    sqrtZr4r2A = R2R4(atomA);
    polarizabilities[idAtm] = 0.0;
    for (size_t idrefA = 0; idrefA < refA[idAtm]; ++idrefA) {
      gweight[iposA + idrefA] *= ZetaD4(ga,hardness*gc,auxv[idrefA] + zA,QAtoms[idAtm] + zA);
      polarizabilities[idAtm] += gweight[iposA + idrefA]*alphapol(1,idrefA + 1,idAtm + 1);
    }
    iposA += refA[idAtm];
  }
}

#endif //_D3_Functions_
