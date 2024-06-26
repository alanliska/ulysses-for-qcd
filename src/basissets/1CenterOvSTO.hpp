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

#ifndef _One_Center_Overlaps_
#define _One_Center_Overlaps_
#include <vector>
#include <math.h>
#include "OvRot.hpp"

//description:
//functions for the calculation of one-center overlap functions

double FactorialFactor(size_t n, size_t N) {
  double factor = pow(2.0,double(n + N + 1));
  if (n > N) {
    for (size_t idx = (2*N + 1); idx < (2*n + 1); ++idx) {
      factor /= sqrt(double(idx));
    }
    for (size_t idx = (2*N + 1); idx < (n + N + 1); ++idx) {
      factor *= double(idx);
    }
  }
  else if (N > n) {
    for (size_t idx = (2*n + 1); idx < (2*N + 1); ++idx) {
      factor /= sqrt(double(idx));
    }
    for (size_t idx = (2*n + 1); idx < (n + N + 1); ++idx) {
      factor *= double(idx);
    }
  }
  return factor;
}
//implementation only up to shell f
void OneCenterOv(matrixE & overlaps, size_t naos, std::vector<double> & exponents, matrix<size_t> & QNumber) {
  //function that receives an atom (number) and returns a matrix with the overlap functions
  overlaps = Identity(naos);
  double aexp = 0.0;
  double bexp = 0.0;
  double aux_a = 0.0;
  double aux_b = 0.0;
  size_t aos_shellbra;
  size_t aos_shellket;
  size_t irow = 1;
  size_t icol = 1;
  size_t nrows = QNumber.rows();
  for (size_t bra_n = 1; bra_n < nrows + 1; ++bra_n) {
    aexp = exponents[bra_n - 1];
    if (aexp == 0.0) continue;
    aos_shellbra = 2*QNumber(bra_n,2) + 1;
    for (size_t imbra = 0; imbra < aos_shellbra; ++imbra, ++irow) {
      icol = 1;
      for (size_t ket_n = 1; ket_n < nrows + 1; ++ket_n) {
        bexp = exponents[ket_n - 1];
        if (bexp == 0.0) continue;
        aux_a = aexp/(aexp + bexp);
        aux_b = bexp/(aexp + bexp);
        aos_shellket = 2*QNumber(ket_n,2) + 1;
        for (size_t imket = 0; imket < aos_shellket; ++imket, ++icol) {
          if ((imket != imbra)||(aos_shellbra != aos_shellket)) continue;
          overlaps(irow,icol) = pow(aux_a,double(QNumber(bra_n,1)) + 0.5) * pow(aux_b,double(QNumber(ket_n,1)) + 0.5) * FactorialFactor(QNumber(bra_n,1),QNumber(ket_n,1));
        }
      }
    }
  }
}

#endif //_One_Center_Overlaps_
