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

#ifndef _Overlap_Functions_
#define _Overlap_Functions_
#include <vector>
#include "../math/MatrixPackage.hpp"
#include "../math/FunctionPackage.hpp"
#include "../ConstantsPackage.hpp"

//description:
//auxiliary functions for the calculation of overlap integrals

std::vector<std::vector<int> > SSOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for S-S overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  std::vector<int> Ak;
  std::vector<int> Bk;
  std::vector<int> coeffs;
  if (n > m) {
    //get the term (x^2 - y^2)^m
    coeffs = BinomialCoefficients(m);
    for (size_t idx = 0; idx < coeffs.size(); ++idx) {
      if ((coeffs.size() - idx - 1)%2 != 0) {coeffs[idx] *= -1;}
      Ak.push_back(int(idx*2));
      Bk.push_back(int((coeffs.size() - idx - 1)*2));
    }
    //get the term (x + y)^(n-m)
    std::vector<int> xk;
    std::vector<int> yk;
    std::vector<int> coeffs2 = BinomialCoefficients(n-m);
    for (size_t idx = 0; idx < coeffs2.size(); ++idx) {
      xk.push_back(int(idx));
      yk.push_back(int(coeffs2.size() - idx - 1));
    }
    ContractPolynomials(coeffs,Ak,Bk,coeffs2,xk,yk);
  }
  else if (n < m) {
    //get the term (x^2 - y^2)^n
    coeffs = BinomialCoefficients(n);
    for (size_t idx = 0; idx < coeffs.size(); ++idx) {
      if ((coeffs.size() - idx - 1)%2 != 0) {coeffs[idx] *= -1;}
      Ak.push_back(int(idx*2));
      Bk.push_back(int((coeffs.size() - idx - 1)*2));
    }
    //get the term (x - y)^(m-n)
    std::vector<int> xk;
    std::vector<int> yk;
    std::vector<int> coeffs2 = BinomialCoefficients(m-n);
    for (size_t idx = 0; idx < coeffs2.size(); ++idx) {
      if ((coeffs2.size() - idx - 1)%2 != 0) {coeffs2[idx] *= -1;}
      xk.push_back(int(idx));
      yk.push_back(int(coeffs2.size() - idx - 1));
    }
    ContractPolynomials(coeffs,Ak,Bk,coeffs2,xk,yk);
  }
  else {
    //get the term (x^2 - y^2)^n
    coeffs = BinomialCoefficients(n);
    for (size_t idx = 0; idx < coeffs.size(); ++idx) {
      if ((coeffs.size() - idx - 1)%2 != 0) {coeffs[idx] *= -1;}
      Ak.push_back(int(idx*2));
      Bk.push_back(int((coeffs.size() - idx - 1)*2));
    }
  }
  std::vector<std::vector<int> > result;
  result.push_back(coeffs);
  result.push_back(Ak);
  result.push_back(Bk);
  return result;
}
std::vector<std::vector<int> > PPsigmaOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for P-P sigma overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-1,m-1);
  //get the P-P sigma specific polynomial
  std::vector<int> xk(2,0);
  xk[1] = 2;
  std::vector<int> yk(2,0);
  yk[1] = 2;
  std::vector<int> coeffs(2,1);
  coeffs[1] = -1;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > PPpiOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for P-P pi overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-1,m-1);
  //get the P-P pi specific polynomial
  std::vector<int> xk(4,0);
  xk[0] = 2;
  xk[2] = 2;
  std::vector<int> yk(4,0);
  yk[1] = 2;
  yk[2] = 2;
  std::vector<int> coeffs(4,1);
  coeffs[2] = -1;
  coeffs[3] = -1;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > SPOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for S-P sigma overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n,m-1);
  //get the S-P sigma specific polynomial
  std::vector<int> xk(2,0);
  xk[1] = 1;
  std::vector<int> yk(2,0);
  yk[1] = 1;
  std::vector<int> coeffs(2,1);
  coeffs[1] = -1;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > PSOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for P-S sigma overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-1,m);
  //get the P-S sigma specific polynomial
  std::vector<int> xk(2,0);
  xk[1] = 1;
  std::vector<int> yk(2,0);
  yk[1] = 1;
  std::vector<int> coeffs(2,1);
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > DDsigmaOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for D-D sigma overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-2,m-2);
  //get the D-D sigma specific polynomial
  std::vector<int> xk(9,0);
  xk[1] = 2;
  xk[2] = 4;
  xk[3] = 2;
  xk[4] = 4;
  xk[5] = 2;
  xk[7] = 4;
  std::vector<int> yk(9,0);
  yk[1] = 2;
  yk[2] = 4;
  yk[3] = 4;
  yk[4] = 2;
  yk[6] = 2;
  yk[8] = 4;
  std::vector<int> coeffs(9,-6);
  coeffs[0] = 9;
  coeffs[1] = 4;
  coeffs[2] = 9;
  coeffs[7] = 1;
  coeffs[8] = 1;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > DDpiOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for D-D pi overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-2,m-2);
  //get the D-D pi specific polynomial
  std::vector<int> xk(6,2);
  xk[0] = 2;
  xk[1] = 0;
  xk[2] = 4;
  xk[3] = 2;
  xk[4] = 0;
  xk[5] = 4;
  std::vector<int> yk(6,2);
  yk[0] = 0;
  yk[1] = 2;
  yk[2] = 2;
  yk[3] = 4;
  yk[4] = 0;
  yk[5] = 4;
  std::vector<int> coeffs(6,1);
  coeffs[0] = 1;
  coeffs[1] = 1;
  coeffs[2] = -1;
  coeffs[3] = -1;
  coeffs[4] = -1;
  coeffs[5] = 1;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > DDdeltaOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for D-D delta overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-2,m-2);
  //get the D-D delta specific polynomial
  std::vector<int> xk(9,2);
  xk[0] = 4;
  xk[1] = 4;
  xk[2] = 4;
  xk[6] = 0;
  xk[7] = 0;
  xk[8] = 0;
  std::vector<int> yk(9,2);
  yk[0] = 4;
  yk[2] = 0;
  yk[3] = 4;
  yk[5] = 0;
  yk[6] = 4;
  yk[8] = 0;
  std::vector<int> coeffs(9,1);
  coeffs[1] = -2;
  coeffs[3] = -2;
  coeffs[4] = 4;
  coeffs[5] = -2;
  coeffs[7] = -2;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > SDOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for S-D sigma overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n,m-2);
  //get the S-D sigma specific polynomial
  std::vector<int> xk(5,0);
  xk[1] = 1;
  xk[2] = 2;
  xk[3] = 2;
  std::vector<int> yk(5,0);
  yk[1] = 1;
  yk[2] = 2;
  yk[4] = 2;
  std::vector<int> coeffs(5,-1);
  coeffs[0] = 3;
  coeffs[1] = -4;
  coeffs[2] = 3;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > DSOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for D-S sigma overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-2,m);
  //get the D-S sigma specific polynomial
  std::vector<int> xk(5,0);
  xk[1] = 1;
  xk[2] = 2;
  xk[3] = 2;
  std::vector<int> yk(5,0);
  yk[1] = 1;
  yk[2] = 2;
  yk[4] = 2;
  std::vector<int> coeffs(5,-1);
  coeffs[0] = 3;
  coeffs[1] = 4;
  coeffs[2] = 3;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > PDsigmaOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for P-D sigma overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-1,m-2);
  //get the P-D sigma specific polynomial
  std::vector<int> xk(8,0);
  xk[1] = 1;
  xk[2] = 2;
  xk[4] = 2;
  xk[5] = 3;
  xk[6] = 3;
  xk[7] = 1;
  std::vector<int> yk(8,0);
  yk[1] = 1;
  yk[3] = 2;
  yk[4] = 2;
  yk[5] = 3;
  yk[6] = 1;
  yk[7] = 3;
  std::vector<int> coeffs(8,-1);
  coeffs[0] = 3;
  coeffs[5] = 3;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > PDpiOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for P-D pi overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-1,m-2);
  //get the P-D pi specific polynomial
  std::vector<int> xk(8,0);
  xk[0] = 2;
  xk[3] = 1;
  xk[4] = 2;
  xk[5] = 3;
  xk[6] = 3;
  xk[7] = 1;
  std::vector<int> yk(8,0);
  yk[1] = 2;
  yk[3] = 1;
  yk[4] = 2;
  yk[5] = 3;
  yk[6] = 1;
  yk[7] = 3;
  std::vector<int> coeffs(8,1);
  coeffs[2] = -1;
  coeffs[4] = -1;
  coeffs[6] = -1;
  coeffs[7] = -1;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > DPsigmaOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for D-P sigma overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-2,m-1);
  //get the D-P sigma specific polynomial
  std::vector<int> xk(8,0);
  xk[1] = 1;
  xk[2] = 2;
  xk[4] = 2;
  xk[5] = 3;
  xk[6] = 3;
  xk[7] = 1;
  std::vector<int> yk(8,0);
  yk[1] = 1;
  yk[3] = 2;
  yk[4] = 2;
  yk[5] = 3;
  yk[6] = 1;
  yk[7] = 3;
  std::vector<int> coeffs(8,-1);
  coeffs[0] = 3;
  coeffs[1] = 1;
  coeffs[5] = -3;
  coeffs[6] = 1;
  coeffs[7] = 1;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > DPpiOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for D-P pi overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-2,m-1);
  //get the D-P pi specific polynomial
  std::vector<int> xk(8,0);
  xk[0] = 2;
  xk[1] = 2;
  xk[4] = 3;
  xk[5] = 3;
  xk[6] = 1;
  xk[7] = 1;
  std::vector<int> yk(8,0);
  yk[1] = 2;
  yk[3] = 2;
  yk[4] = 1;
  yk[5] = 3;
  yk[6] = 1;
  yk[7] = 3;
  std::vector<int> coeffs(8,-1);
  coeffs[0] = 1;
  coeffs[3] = 1;
  coeffs[4] = 1;
  coeffs[7] = 1;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > FFsigmaOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for F-F sigma overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-3,m-3);
  //get the F-F sigma specific polynomial
  std::vector<int> xk(12,2);
  xk[0] = 0;
  xk[2] = 4;
  xk[3] = 6;
  xk[5] = 0;
  xk[6] = 4;
  xk[7] = 0;
  xk[8] = 6;
  xk[9] = 4;
  xk[10] = 6;
  std::vector<int> yk(12,2);
  yk[0] = 0;
  yk[2] = 4;
  yk[3] = 6;
  yk[4] = 0;
  yk[6] = 0;
  yk[7] = 4;
  yk[8] = 4;
  yk[9] = 6;
  yk[11] = 6;
  std::vector<int> coeffs(12,9);
  coeffs[0] = 25;
  coeffs[1] = 27;
  coeffs[2] = -27;
  coeffs[3] = -25;
  coeffs[4] = -30;
  coeffs[5] = -30;
  coeffs[8] = 30;
  coeffs[9] = 30;
  coeffs[10] = -9;
  coeffs[11] = -9;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > FFpiOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for F-F pi overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-3,m-3);
  //get the F-F pi specific polynomial
  std::vector<int> xk(16,2);
  xk[0] = 6;
  xk[1] = 0;
  xk[2] = 4;
  xk[3] = 0;
  xk[5] = 0;
  xk[6] = 0;
  xk[8] = 4;
  xk[9] = 6;
  xk[10] = 4;
  xk[13] = 6;
  xk[14] = 6;
  xk[15] = 4;
  std::vector<int> yk(16,2);
  yk[0] = 0;
  yk[1] = 6;
  yk[2] = 0;
  yk[3] = 4;
  yk[4] = 0;
  yk[6] = 0;
  yk[8] = 4;
  yk[9] = 6;
  yk[11] = 4;
  yk[12] = 6;
  yk[14] = 4;
  yk[15] = 6;
  std::vector<int> coeffs(16,-11);
  coeffs[0] = 1;
  coeffs[1] = 1;
  coeffs[4] = 35;
  coeffs[5] = 35;
  coeffs[6] = -25;
  coeffs[7] = -33;
  coeffs[8] = -33;
  coeffs[9] = -25;
  coeffs[10] = 9;
  coeffs[11] = 9;
  coeffs[14] = 35;
  coeffs[15] = 35;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > FFdeltaOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for F-F delta overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-3,m-3);
  //get the F-F delta specific polynomial
  std::vector<int> xk(12,2);
  xk[0] = 4;
  xk[1] = 0;
  xk[3] = 0;
  xk[4] = 0;
  xk[6] = 4;
  xk[7] = 6;
  xk[8] = 6;
  xk[10] = 6;
  xk[11] = 4;
  std::vector<int> yk(12,2);
  yk[0] = 0;
  yk[1] = 4;
  yk[2] = 0;
  yk[4] = 0;
  yk[6] = 4;
  yk[7] = 6;
  yk[9] = 6;
  yk[10] = 4;
  yk[11] = 6;
  std::vector<int> coeffs(12,-1);
  coeffs[0] = 1;
  coeffs[1] = 1;
  coeffs[2] = -2;
  coeffs[3] = -2;
  coeffs[4] = 1;
  coeffs[5] = 3;
  coeffs[6] = -3;
  coeffs[10] = 2;
  coeffs[11] = 2;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > FFphiOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for F-F phi overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-3,m-3);
  //get the F-F phi specific polynomial
  std::vector<int> xk(16,2);
  xk[0] = 6;
  xk[1] = 0;
  xk[2] = 4;
  xk[3] = 0;
  xk[5] = 0;
  xk[6] = 0;
  xk[7] = 6;
  xk[8] = 4;
  xk[10] = 6;
  xk[12] = 6;
  xk[13] = 4;
  xk[14] = 4;
  std::vector<int> yk(16,2);
  yk[0] = 0;
  yk[1] = 6;
  yk[2] = 0;
  yk[3] = 4;
  yk[4] = 0;
  yk[6] = 0;
  yk[7] = 6;
  yk[8] = 4;
  yk[11] = 6;
  yk[12] = 4;
  yk[13] = 6;
  yk[15] = 4;
  std::vector<int> coeffs(16,-3);
  coeffs[0] = 1;
  coeffs[1] = 1;
  coeffs[4] = 3;
  coeffs[5] = 3;
  coeffs[6] = -1;
  coeffs[7] = -1;
  coeffs[8] = -9;
  coeffs[9] = -9;
  coeffs[12] = 3;
  coeffs[13] = 3;
  coeffs[14] = 9;
  coeffs[15] = 9;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > SFOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for S-F sigma overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n,m-3);
  //get the S-F sigma specific polynomial
  std::vector<int> xk(8,0);
  xk[1] = 1;
  xk[2] = 2;
  xk[3] = 3;
  xk[4] = 2;
  xk[6] = 3;
  xk[7] = 1;
  std::vector<int> yk(8,0);
  yk[1] = 1;
  yk[2] = 2;
  yk[3] = 3;
  yk[5] = 2;
  yk[6] = 1;
  yk[7] = 3;
  std::vector<int> coeffs(8,-3);
  coeffs[0] = 5;
  coeffs[1] = -9;
  coeffs[2] = 9;
  coeffs[3] = -5;
  coeffs[6] = 3;
  coeffs[7] = 3;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > FSOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for F-S sigma overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-3,m);
  //get the F-S sigma specific polynomial
  std::vector<int> xk(8,0);
  xk[1] = 1;
  xk[2] = 2;
  xk[3] = 3;
  xk[4] = 2;
  xk[6] = 3;
  xk[7] = 1;
  std::vector<int> yk(8,0);
  yk[1] = 1;
  yk[2] = 2;
  yk[3] = 3;
  yk[5] = 2;
  yk[6] = 1;
  yk[7] = 3;
  std::vector<int> coeffs(8,-3);
  coeffs[0] = 5;
  coeffs[1] = 9;
  coeffs[2] = 9;
  coeffs[3] = 5;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > PFsigmaOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for P-F sigma overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-1,m-3);
  //get the P-F sigma specific polynomial
  std::vector<int> xk(8,0);
  xk[1] = 1;
  xk[2] = 3;
  xk[3] = 4;
  xk[4] = 2;
  xk[6] = 4;
  xk[7] = 2;
  std::vector<int> yk(8,0);
  yk[1] = 1;
  yk[2] = 3;
  yk[3] = 4;
  yk[5] = 2;
  yk[6] = 2;
  yk[7] = 4;
  std::vector<int> coeffs(8,-3);
  coeffs[0] = 5;
  coeffs[1] = -4;
  coeffs[2] = 4;
  coeffs[3] = -5;
  coeffs[6] = 3;
  coeffs[7] = 3;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > PFpiOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for P-F pi overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-1,m-3);
  //get the P-F pi specific polynomial
  std::vector<int> xk(13,0);
  xk[0] = 2;
  xk[2] = 4;
  xk[5] = 1;
  xk[6] = 2;
  xk[7] = 3;
  xk[8] = 4;
  xk[9] = 3;
  xk[10] = 1;
  xk[11] = 4;
  xk[12] = 2;
  std::vector<int> yk(13,0);
  yk[1] = 2;
  yk[3] = 4;
  yk[5] = 1;
  yk[6] = 2;
  yk[7] = 3;
  yk[8] = 4;
  yk[9] = 1;
  yk[10] = 3;
  yk[11] = 2;
  yk[12] = 4;
  std::vector<int> coeffs(13,6);
  coeffs[2] = -1;
  coeffs[3] = -1;
  coeffs[4] = -5;
  coeffs[5] = 8;
  coeffs[6] = -12;
  coeffs[7] = 8;
  coeffs[8] = -5;
  coeffs[9] = -8;
  coeffs[10] = -8;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > FPsigmaOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for F-P sigma overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-3,m-1);
  //get the F-P sigma specific polynomial
  std::vector<int> xk(8,0);
  xk[1] = 1;
  xk[2] = 3;
  xk[3] = 4;
  xk[4] = 2;
  xk[6] = 4;
  xk[7] = 2;
  std::vector<int> yk(8,0);
  yk[1] = 1;
  yk[2] = 3;
  yk[3] = 4;
  yk[5] = 2;
  yk[6] = 2;
  yk[7] = 4;
  std::vector<int> coeffs(8,-3);
  coeffs[0] = 5;
  coeffs[1] = 4;
  coeffs[2] = -4;
  coeffs[3] = -5;
  coeffs[6] = 3;
  coeffs[7] = 3;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > FPpiOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for F-P pi overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-3,m-1);
  //get the F-P pi specific polynomial
  std::vector<int> xk(13,0);
  xk[0] = 2;
  xk[3] = 1;
  xk[4] = 2;
  xk[5] = 3;
  xk[6] = 4;
  xk[7] = 4;
  xk[9] = 3;
  xk[10] = 1;
  xk[11] = 4;
  xk[12] = 2;
  std::vector<int> yk(13,0);
  yk[1] = 2;
  yk[3] = 1;
  yk[4] = 2;
  yk[5] = 3;
  yk[6] = 4;
  yk[8] = 4;
  yk[9] = 1;
  yk[10] = 3;
  yk[11] = 2;
  yk[12] = 4;
  std::vector<int> coeffs(13,6);
  coeffs[2] = -5;
  coeffs[3] = -8;
  coeffs[4] = -12;
  coeffs[5] = -8;
  coeffs[6] = -5;
  coeffs[7] = -1;
  coeffs[8] = -1;
  coeffs[9] = 8;
  coeffs[10] = 8;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > DFsigmaOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for D-F sigma overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-2,m-3);
  //get the D-F sigma specific polynomial
  std::vector<int> xk(18,0);
  xk[1] = 1;
  xk[2] = 2;
  xk[3] = 3;
  xk[4] = 4;
  xk[5] = 5;
  xk[6] = 2;
  xk[8] = 4;
  xk[10] = 3;
  xk[11] = 1;
  xk[12] = 4;
  xk[13] = 2;
  xk[14] = 5;
  xk[15] = 3;
  xk[16] = 5;
  xk[17] = 1;
  std::vector<int> yk(18,0);
  yk[1] = 1;
  yk[2] = 2;
  yk[3] = 3;
  yk[4] = 4;
  yk[5] = 5;
  yk[7] = 2;
  yk[9] = 4;
  yk[10] = 1;
  yk[11] = 3;
  yk[12] = 2;
  yk[13] = 4;
  yk[14] = 3;
  yk[15] = 5;
  yk[16] = 1;
  yk[17] = 5;
  std::vector<int> coeffs(18,-14);
  coeffs[0] = 15;
  coeffs[1] = -7;
  coeffs[2] = 12;
  coeffs[3] = -12;
  coeffs[4] = 7;
  coeffs[5] = -15;
  coeffs[8] = 3;
  coeffs[9] = 3;
  coeffs[10] = 6;
  coeffs[11] = 6;
  coeffs[12] = -6;
  coeffs[13] = -6;
  coeffs[14] = 14;
  coeffs[15] = 14;
  coeffs[16] = -3;
  coeffs[17] = -3;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > DFpiOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for D-F pi overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-2,m-3);
  //get the D-F pi specific polynomial
  std::vector<int> xk(18,0);
  xk[0] = 2;
  xk[2] = 4;
  xk[5] = 2;
  xk[6] = 1;
  xk[7] = 3;
  xk[8] = 4;
  xk[9] = 5;
  xk[10] = 3;
  xk[11] = 1;
  xk[12] = 4;
  xk[13] = 2;
  xk[14] = 5;
  xk[15] = 3;
  xk[16] = 5;
  xk[17] = 1;
  std::vector<int> yk(18,0);
  yk[1] = 2;
  yk[3] = 4;
  yk[5] = 2;
  yk[6] = 1;
  yk[7] = 3;
  yk[8] = 4;
  yk[9] = 5;
  yk[10] = 1;
  yk[11] = 3;
  yk[12] = 2;
  yk[13] = 4;
  yk[14] = 3;
  yk[15] = 5;
  yk[16] = 1;
  yk[17] = 5;
  std::vector<int> coeffs(18,6);
  coeffs[2] = -1;
  coeffs[3] = -1;
  coeffs[4] = -5;
  coeffs[5] = -4;
  coeffs[6] = 3;
  coeffs[7] = -4;
  coeffs[8] = 3;
  coeffs[9] = -5;
  coeffs[10] = -2;
  coeffs[11] = -2;
  coeffs[12] = -2;
  coeffs[13] = -2;
  coeffs[16] = -1;
  coeffs[17] = -1;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > DFdeltaOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for D-F delta overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-2,m-3);
  //get the D-F delta specific polynomial
  std::vector<int> xk(18,0);
  xk[0] = 4;
  xk[2] = 2;
  xk[5] = 1;
  xk[6] = 2;
  xk[7] = 3;
  xk[8] = 4;
  xk[9] = 5;
  xk[10] = 4;
  xk[11] = 2;
  xk[12] = 5;
  xk[13] = 1;
  xk[14] = 5;
  xk[15] = 3;
  xk[16] = 3;
  xk[17] = 1;
  std::vector<int> yk(18,0);
  yk[1] = 4;
  yk[3] = 2;
  yk[5] = 1;
  yk[6] = 2;
  yk[7] = 3;
  yk[8] = 4;
  yk[9] = 5;
  yk[10] = 2;
  yk[11] = 4;
  yk[12] = 1;
  yk[13] = 5;
  yk[14] = 3;
  yk[15] = 5;
  yk[16] = 1;
  yk[17] = 3;
  std::vector<int> coeffs(18,2);
  coeffs[0] = 1;
  coeffs[1] = 1;
  coeffs[2] = -2;
  coeffs[3] = -2;
  coeffs[4] = 1;
  coeffs[5] = -1;
  coeffs[6] = 4;
  coeffs[7] = -4;
  coeffs[8] = 1;
  coeffs[9] = -1;
  coeffs[10] = -2;
  coeffs[11] = -2;
  coeffs[12] = -1;
  coeffs[13] = -1;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > FDsigmaOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for F-D sigma overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-3,m-2);
  //get the F-D sigma specific polynomial
  std::vector<int> xk(18,0);
  xk[1] = 1;
  xk[2] = 2;
  xk[3] = 3;
  xk[4] = 4;
  xk[5] = 5;
  xk[6] = 2;
  xk[8] = 4;
  xk[10] = 3;
  xk[11] = 1;
  xk[12] = 4;
  xk[13] = 2;
  xk[14] = 5;
  xk[15] = 3;
  xk[16] = 5;
  xk[17] = 1;
  std::vector<int> yk(18,0);
  yk[1] = 1;
  yk[2] = 2;
  yk[3] = 3;
  yk[4] = 4;
  yk[5] = 5;
  yk[7] = 2;
  yk[9] = 4;
  yk[10] = 1;
  yk[11] = 3;
  yk[12] = 2;
  yk[13] = 4;
  yk[14] = 3;
  yk[15] = 5;
  yk[16] = 1;
  yk[17] = 5;
  std::vector<int> coeffs(18,3);
  coeffs[0] = 15;
  coeffs[1] = 7;
  coeffs[2] = 12;
  coeffs[3] = 12;
  coeffs[4] = 7;
  coeffs[5] = 15;
  coeffs[6] = -14;
  coeffs[7] = -14;
  coeffs[10] = -6;
  coeffs[11] = -6;
  coeffs[12] = -6;
  coeffs[13] = -6;
  coeffs[14] = -14;
  coeffs[15] = -14;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > FDpiOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for F-D pi overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-3,m-2);
  //get the F-D pi specific polynomial
  std::vector<int> xk(18,0);
  xk[1] = 1;
  xk[2] = 2;
  xk[3] = 3;
  xk[4] = 4;
  xk[5] = 5;
  xk[6] = 2;
  xk[8] = 4;
  xk[10] = 3;
  xk[11] = 1;
  xk[12] = 4;
  xk[13] = 2;
  xk[14] = 5;
  xk[15] = 3;
  xk[16] = 5;
  xk[17] = 1;
  std::vector<int> yk(18,0);
  yk[1] = 1;
  yk[2] = 2;
  yk[3] = 3;
  yk[4] = 4;
  yk[5] = 5;
  yk[7] = 2;
  yk[9] = 4;
  yk[10] = 1;
  yk[11] = 3;
  yk[12] = 2;
  yk[13] = 4;
  yk[14] = 3;
  yk[15] = 5;
  yk[16] = 1;
  yk[17] = 5;
  std::vector<int> coeffs(18,6);
  coeffs[0] = -5;
  coeffs[1] = -3;
  coeffs[2] = -4;
  coeffs[3] = 4;
  coeffs[4] = 3;
  coeffs[5] = 5;
  coeffs[8] = -1;
  coeffs[9] = -1;
  coeffs[10] = 2;
  coeffs[11] = 2;
  coeffs[12] = -2;
  coeffs[13] = -2;
  coeffs[14] = -6;
  coeffs[15] = -6;
  coeffs[16] = 1;
  coeffs[17] = 1;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
std::vector<std::vector<int> > FDdeltaOvIntIndex(size_t n, size_t m) {
  //function that calculates vector indices for F-D delta overlaps
  //first vector is to use in Ak, second with Bk
  //coefficients must be mandatorily integer; that is why everything packed as vectors of ints
  //get the SS polynomial part
  std::vector<std::vector<int> > result = SSOvIntIndex(n-3,m-2);
  //get the F-D delta specific polynomial
  std::vector<int> xk(18,0);
  xk[0] = 4;
  xk[2] = 2;
  xk[5] = 1;
  xk[6] = 2;
  xk[7] = 3;
  xk[8] = 4;
  xk[9] = 5;
  xk[10] = 4;
  xk[11] = 2;
  xk[12] = 1;
  xk[13] = 5;
  xk[14] = 5;
  xk[15] = 3;
  xk[16] = 3;
  xk[17] = 1;
  std::vector<int> yk(18,0);
  yk[1] = 4;
  yk[3] = 2;
  yk[5] = 1;
  yk[6] = 2;
  yk[7] = 3;
  yk[8] = 4;
  yk[9] = 5;
  yk[10] = 2;
  yk[11] = 4;
  yk[12] = 5;
  yk[13] = 1;
  yk[14] = 3;
  yk[15] = 5;
  yk[16] = 1;
  yk[17] = 3;
  std::vector<int> coeffs(18,-2);
  coeffs[0] = 1;
  coeffs[1] = 1;
  coeffs[2] = -2;
  coeffs[3] = -2;
  coeffs[4] = 1;
  coeffs[5] = 1;
  coeffs[6] = 4;
  coeffs[7] = 4;
  coeffs[8] = 1;
  coeffs[9] = 1;
  coeffs[10] = -2;
  coeffs[11] = -2;
  coeffs[12] = 1;
  coeffs[13] = 1;
  ContractPolynomials(result[0],result[1],result[2],coeffs,xk,yk);
  return result;
}
#endif //_Overlap_Functions_