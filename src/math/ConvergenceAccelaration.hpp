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

#ifndef _Convergence_Acceleration_Package_
#define _Convergence_Acceleration_Package_
#include <vector>
#include <math.h>
#include "MatrixPackage.hpp"
#include "SolverPackage.hpp"

//description:
//methods for convergence accelaration in quantum chemistry

void DIISse(std::vector<matrixE> & error, std::vector<matrixE> & Fock, std::vector<matrixE> & bFock, matrixE & dens, matrixE & bdens, matrixE & Bij, std::vector<double> & occupancy, std::vector<double> & boccupancy, size_t & DIISstatus, matrixE & Fprime, std::vector<double> & avec1, std::vector<double> & avec2, matrixE & DensAUX) {
  //function performing the DIIS procedure to accelerate SCF calculations according to 
  //P. Pulay, J. Comput. Chem., 3(4), 556, 1982
  //this is the simplified version for semi-empirical methods
  //error is vector that stores error matrices
  //Fock is matrix storing Fock matrices
  //dens is density matrix
  //Bij is error matrix required in DIIS
  //DIIS status stores information about what to do with DIIS: 
  //          0 -> nothing
  //          1 -> perform DIIS
  //          2 -> DIIS converged
  double threshold1 = 0.015;          //threshold determining when to start DIIS
  double threshold2 = 1.0e-6;         //threshold determining convergence
  double maxerr = 0.0;
  int NAOs = dens.rows();
  //calculate error matrix:      E = FD - DF
  size_t posFock = Fock.size() - 1;
  bool UHF = (bFock.size() > 0);
  Fprime = Fock[posFock]*dens - dens*Fock[posFock];
  if (UHF) {Fprime += (bFock[posFock]*bdens - bdens*bFock[posFock]);}
  //look for maximum of error matrix
  for (size_t idr = 0; idr < NAOs; ++idr) {
    for (size_t idc = idr; idc < NAOs; ++idc) {
      if (fabs(Fprime(idr + 1,idc + 1)) > maxerr) maxerr = fabs(Fprime(idr + 1,idc + 1));
    }
  }
  //check conditions
  if (maxerr < threshold2) {DIISstatus = 2;}          //converge
  else {
    //start DIIS
    size_t dim = error.size() + 1;
    Bij.resize(Bij.rows() + 1,Bij.cols() + 1);
    Bij(1,dim + 1) = -1.0;
    Bij(dim + 1,1) = -1.0;
    double dtrace;
    for (size_t idx = 0; idx < dim - 1; ++idx) {
      dtrace = trace(Fprime,error[idx]);
      Bij(idx + 2,dim + 1) = dtrace;
      Bij(dim + 1,idx + 2) = dtrace;
    }
    dtrace = trace(Fprime,Fprime);
    Bij(dim + 1,dim + 1) = dtrace;
    error.push_back(Fprime);
    if (maxerr < threshold1) {
      DIISstatus = 1;
      //solve the system of equations Bx = (-1,0,0,...0) = avec
      avec1.resize(Bij.rows());
      avec1[0] = -1.0;
      for (size_t idx = 1; idx < Bij.rows(); ++idx) {
        avec1[idx] = 0.0;
      }
      Solve_Ax_eq_b(Bij,avec1);
      if (UHF) {avec2 = avec1;}
      //get DIIS Fock matrix F'
      Fprime.zero();
      for (size_t idc = 1; idc < avec1.size(); ++idc) {
        Fprime += Fock[idc - 1]*avec1[idc];
      }
      //diagonalize and get new density
      avec1 = MatDiag(Fprime);
      Dens(dens,occupancy,Fprime,DensAUX);
      if (UHF) {
        Fprime.zero();
        for (size_t idc = 1; idc < avec2.size(); ++idc) {
          Fprime += bFock[idc - 1]*avec2[idc];
        }
        //diagonalize and get new density
        avec2 = MatDiag(Fprime);
        Dens(bdens,boccupancy,Fprime,DensAUX);
      }
    }
  }
}
void RCA(std::vector<matrixE> & Fock,std::vector<matrixE> & bFock, std::vector<matrixE> & dens, std::vector<matrixE> & bdens, std::vector<double> & occupancy, std::vector<double> & boccupancy, double energy, double oldenergy, matrixE & Fprime, std::vector<double> & eigvec, matrixE & DensAUX, double threshold = 1.e-7) {
  //function performing the RCA procedure to force converge into SCF calculations according to 
  //E. Cances, C. le Bris, Int. J. Quantum Chem., Vol. 79, 82, 2000
  //Fock is the vector containing Fock matrices
  //dens is a vector containing the last two densities
  //energy and old energy are the energies of current and previous iterations
  //set pointers
  bool UHF = (bdens.size() > 0);
  int thisiter = Fock.size() - 1;
  double c = trace(Fock[thisiter - 1],dens[1]) - trace(Fock[thisiter - 1],dens[0]);
  double d = oldenergy;
  double a = trace(Fock[thisiter],dens[1]) - trace(Fock[thisiter],dens[0]);
  if (UHF) {
    a += trace(bFock[thisiter],bdens[1]) - trace(bFock[thisiter],bdens[0]);
    c += trace(bFock[thisiter - 1],bdens[1]) - trace(bFock[thisiter - 1],bdens[0]);
  }
  a += c + 2.0*(oldenergy - energy);
  double b = energy - a - c - d;
  double fmin = energy;
  double xmin = 1.0;
  if (oldenergy < energy) {
    fmin = oldenergy;
    xmin = 0.0;
  }
  double x;
  double xt1;
  double xt2;
  double fac;
  //solve cubic polynomial a*x*x*x + b*x*x + c*x + d
  if (fabs(a) < threshold) {
    //case of quadratic polynomial
    if (fabs(b) < threshold) {x = xmin;}
    else {
      xt1 = -0.5*c/b;
      if ((b > 0.0)&&((xt1 > 0.0)&&(xt1 < 1.0))) {x = xt1;}
      else {x = xmin;}
    }
  }
  else {
    fac = b*b - 3.0*a*c;
    if (fabs(fac) < threshold) {x = xmin;}
    else {
      xt1 = (-b + sqrt(fac))/(3.0*a);
      xt2 = (-b - sqrt(fac))/(3.0*a);
      if ((xt1 > 0.0)&&(xt1 < 1.0)) {
        fac = d + xt1*(c + xt1*(b + xt1*a));
        if (fac < fmin) {
          fmin = fac;
          xmin = xt1;
        }
      }
      if ((xt2 > 0.0)&&(xt2 < 1.0)) {
        fac = d + xt2*(c + xt2*(b + xt2*a));
        if (fac < fmin) {
          fmin = fac;
          xmin = xt2;
        }
      }
      x = xmin;
    }
  }
  x = fmax(0.01,x);
  //build new Fock, diagonalize and get new density
  Fprime = Fock[thisiter]*x + Fock[thisiter - 1]*(1.0 - x);
  eigvec = MatDiag(Fprime);
  Dens(dens[1],occupancy,Fprime,DensAUX);
  if (UHF) {
    Fprime = bFock[thisiter]*x + bFock[thisiter - 1]*(1.0 - x);
    eigvec = MatDiag(Fprime);
    Dens(bdens[1],boccupancy,Fprime,DensAUX);
  }
}
#endif //_Convergence_Acceleration_Package_
