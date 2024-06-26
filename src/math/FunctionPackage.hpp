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

#ifndef _Function_Package_
#define _Function_Package_
#include <math.h>
#include <vector>
#include <string>
#include <complex>
#include "../ConstantsPackage.hpp"
#include "MatrixPackage.hpp"

//description:
//package with general functions
//check the math documentation in the folder for more details on what is in here

double FactorialRatio(size_t a, size_t b) {
  //function that calculates the ratio between 2 factorials: a!/b!
  if (a == 0) {a = 1;}
  if (b == 0) {b = 1;}
  size_t MIN = fmin(a,b);
  size_t MAX = fmax(a,b);
  ++MIN;
  ++MAX;
  double result = 1;
  for (size_t idx = MIN; idx < MAX; ++idx) {
    result *= idx;
  }
  if (b > a) {
    double aux = result;
    result = 1/aux;
  }
  return result;
}
double FactorialRatio_Stirling(size_t a, size_t b) {
  //function that calculates the ratio between 2 factorials a!/b! using the Stirling approximation: a!/b!
  double logndm = a*log(a) - b*log(b) + b - a;
  return exp(logndm);
}
double Factorial(size_t a) {
  //function that calculates factorials
  double result = 1;
  for (size_t idx = 1; idx < a + 1; ++idx) {
    result *= idx;
  }
  return result;
}
double kFactorial(int a) {
  //function that calculates factorials of negative values, according to Kurba's definition
  double result = 1;
  for (int idx = 0; idx < abs(a); ++idx) {
    result *= (idx + a);
  }
  return result;
}
double rFactorial(int a) {
  //function that calculates factorials according to Roman's definition
  double result;
  if (a >= 0) {result = Factorial(a);}
  else {result = pow(-1.0,double(-a - 1))/Factorial(-a - 1);}
  return result;
}
double DoubleFactorial(int a) {
  //function that calculates double factorials
  double result = 1;
  size_t start = 2;
  if (a%2 != 0) {--start;}
  for (size_t idx = start; idx < a+1; idx += 2) {
    result *= idx;
  }
  return result;
}
double DoubleFactorial2(int a) {
  //function that calculates double factorials with simplified formulae
  double result;
  int k;
  if (a < 0) {a = 0;}
  if (a%2 == 0) {        //a even
    k = a/2;
    result = pow(2.0,double(k))*Factorial(k);
  }
  else {                 //a odd
    k = (a + 1)/2;
    result = FactorialRatio(a,k - 1)*pow(2.0,double(1 - k));
  }
  return result;
}
double hiFactorial(double a) {
  //function that calculates factorials for half integers only!
  double n = a + 0.5;
  double result;
  if (n >= 0) {result = FactorialRatio(2*int(n),int(n))*sqrt(pi)/pow(4.0,n);}
  else {
    double aux = fabs(n);
    n = aux;
    result = sqrt(pi)*pow(-2.0,n)/DoubleFactorial2(2*int(n) - 1);
  }
  return result;
}
double Factorial_Stirling(size_t a) {
  //function that calculates factorials using Stirling's approximation
  //valid for N > 10
  return sqrt(2.0*pi*double(a))*double(pow(a,a))*exp(-double(a));
}
double DoubleFactorialS(int a) {
  //function that calculates double factorials with simplified formulae and the Stirling approximation
  double result;
  int k;
  if (a%2 == 0) {        //a even
    k = a/2;
    result = pow(2.0,double(k))*Factorial_Stirling(k);
  }
  else {                 //a odd
    k = (a + 1)/2;
    result = Factorial_Stirling(a)*pow(2.0,double(1 - k))/Factorial_Stirling(k - 1);
  }
  return result;
}
size_t Combinations(size_t n, size_t r) {
  //function that calculates C(n,r) = n!/(r!(n-r)!)
  //Code stable for n up to 60; afterwards instability at large r
  if (n < r) {throw(std::string("ERROR: FunctionPackage.hpp: Combinations(): n < r"));}
  else {
    double result = 1;
    for (size_t idx = n-r; idx < n; ++idx) {
      result *= (idx + 1);
    }
    double aux = 1;
    for (size_t idx = 0; idx < r; ++idx) {
      aux *= (idx + 1);
    }
    return size_t(result/aux);
  }
}
double dCombinations(size_t n, size_t m) {
  //function that calculates C(n,m) = n!/(m!(n-m)!)
  //this is basically like previous function but with more stability and returns type double
  //tested with n up to 100 and all binomial coefficients were calculated
  double result = 1.0;
  if (n > (m - 1)) {
    size_t MAX = std::max(m,n-m);
    size_t MIN = std::min(m,n-m);
    result *= FactorialRatio(n,MAX);
    for (size_t idx = 1; idx < MIN+1; ++idx) {
      result /= double(idx);
    }
  }
  return result;
}
double Fbinomial(int m, int n) {
  //just a wrapper...
  return dCombinations(size_t(n), size_t(m));
}
std::vector<int> BinomialCoefficients(size_t order) {
  //function that uses Newton's binomial theorem to calculate all the binomial's coefficients, returned as a vector
  //{(x+y)}^{n} = \sum_{k=0}^{n}{C(n,k){x}^{k}{y}^{n-k}}, n = order;
  //even though the coefficients have a symmetry, it is assumed that x begins with exponent 0 and y with exponent n
  std::vector<int> result(order+1,0);
  for (size_t idx = 0; idx < order+1; ++idx) {
    result[idx] = int(Combinations(order,idx));
  }
  return result;
}
void ContractPolynomials(std::vector<int> & coeff1, std::vector<int> & x1, std::vector<int> & y1, const std::vector<int> & coeff2, const std::vector<int> & x2, const std::vector<int> & y2) {
  //function that contracts two-variable polynomials described each by 3 vectors: coefficients, power first variable, power second variable
  // in p = a[0]x[0]y[0] + a[1]x[1]y[1] + a[2]x[2]y[2] + ...
  // in q = b[0]w[0]z[0] + b[1]w[1]z[1] + b[2]w[2]z[2] + ...
  // out p.q overwritten on p
  std::vector<int> fc; //final coefficients
  std::vector<int> fx; //final x
  std::vector<int> fy; //final y
  //auxiliaries
  int cc, xx, yy;
  bool found = false;
  //get size of polynomials
  size_t szp1 = coeff1.size();
  size_t szp2 = coeff2.size();
  //loop over both polynomials
  for (size_t ip1 = 0; ip1 < szp1; ++ip1) {
    for (size_t ip2 = 0; ip2 < szp2; ++ip2) {
      cc = coeff1[ip1] * coeff2[ip2];
      xx = x1[ip1] + x2[ip2];
      yy = y1[ip1] + y2[ip2];
      //loop over resulting polynomial to see whether the new term already exists
      found = false;
      for (size_t idx = 0; idx < fc.size(); ++idx) {
        if ((fx[idx] == xx)&&(fy[idx] == yy)) {
          fc[idx] += cc;
          if (fc[idx] == 0) {
            fc.erase(fc.begin() + idx);
            fx.erase(fx.begin() + idx);
            fy.erase(fy.begin() + idx);
          }
          found = true;
          break;
        }
      }
      if (!found) {
        fc.push_back(cc);
        fx.push_back(xx);
        fy.push_back(yy);
      }
    }
  }
  coeff1 = fc;
  x1 = fx;
  y1 = fy;
}
void ExpandPolynomial(std::vector<std::vector<int> > & Cmni) {
  //function that takes a polynomial of the form (x+y)^m (x-y)^n and expands to get all coefficients
  int maxx = 0;
  int maxy = 0;
  //determination of maximum power in polynomial
  for (size_t idx = 0; idx < Cmni[1].size(); ++idx) {
    if (Cmni[1][idx] > maxx) {maxx = Cmni[1][idx];}
    if (Cmni[2][idx] > maxy) {maxy = Cmni[2][idx];}
  }
  if (maxx != maxy) {throw(std::string("ERROR: FunctionPackage.hpp: ExpandPolynomial(): Cmni vector incorrectly built"));}
  //allocation of auxiliary vectors
  //+1 because zero is also a power to consider!
  std::vector<int> coeffs(maxx + 1,0);
  std::vector<int> px(maxx + 1,0);
  std::vector<int> py(maxx + 1,0);
  //assigning values to the powers and getting the exponents
  for (size_t idx = 0; idx < maxx + 1; ++idx) {
    px[idx] = idx;
    py[idx] = maxx - idx;
    for (size_t idx2 = 0; idx2 < Cmni[0].size(); ++idx2) {
      if ((px[idx] == Cmni[1][idx2])&&(py[idx] == Cmni[2][idx2])) {
        coeffs[idx] = Cmni[0][idx2];
        break;
      }
    }
  }
  Cmni[0] = coeffs;
  Cmni[1] = px;
  Cmni[2] = py;
}
matrixE LPmn(size_t l, size_t m, double x) {
  //matrix function calculating recursively all {P}_{l}^{m}(x) up to given l,m values
  ++l;
  ++m;
  matrixE Lpolynomial(l,m);
  //starting value
  //note that since matrices have > 0 indices then {P}_{0}^{0}(x) is in the position (1,1) of the matrix
  Lpolynomial(1,1) = 1.0;
  if (l - 1 > 1) {
    //get main diagonal elements
    size_t diaglimit = std::min(m,l);
    double sqrtx = 0.0;
    if (abs(x) < 1.0) {sqrtx = -sqrt(1.0 - x*x);}
    else {sqrtx = sqrt(x*x - 1.0);}
    for (size_t idx = 1; idx < diaglimit; ++idx) {
      Lpolynomial(idx + 1,idx + 1) = sqrtx*(2.0*double(idx - 1) + 1.0)*Lpolynomial(idx,idx);
    }
    //get elements below main diagonal
    if (l <= m) {--diaglimit;}
    else {diaglimit = m;}
    for (size_t idx = 0; idx < diaglimit; ++idx) {
      Lpolynomial(idx + 2,idx + 1) = x*(2.0*double(idx) + 1.0)*Lpolynomial(idx + 1,idx + 1);
    }
    //fill rest of matrix
    for (size_t idm = 0; idm < m; ++idm) {
      for (size_t idl = idm+2; idl < l; ++idl) {
        Lpolynomial(idl + 1,idm + 1) = ((2.0*double(idl) - 1.0)/(double(idl) - double(idm)))*x*Lpolynomial(idl,idm + 1);
        Lpolynomial(idl + 1,idm + 1) -= ((double(idl - 1) + double(idm))/(double(idl) - double(idm)))*Lpolynomial(idl - 1,idm + 1);
      }
    }
  }
  return Lpolynomial;
}
std::vector<double> vLPmn(size_t l, size_t m, double x) {
  //vector function calculating recursively all {P}_{l}^{m}(x) up to given l,m values
  //returns only the last column m
  matrixE Lpolynomial = LPmn(l,m,x);
  std::vector<double> lpol(l + 1,0.0);
  for (size_t id = 0; id < l + 1; ++id) {
    lpol[id] = Lpolynomial(id + 1,m + 1);
  }
  return lpol;
}
matrixE LPmn_HM(size_t l, size_t m, double x) {
  //matrix function calculating recursively all {P}_{l}^{m}(x) up to given l,m values
  //using the Harris-Michels method
  ++l;
  ++m;
  matrixE Lpolynomial(l,m);
  //starting value
  //note that since matrices have > 0 indices then {P}_{0}^{0}(x) is in the position (1,1) of the matrix
  Lpolynomial(1,1) = 1.0;
  if (l - 1 > 1) {
    Lpolynomial(2,1) = x;
    if (fabs(x) <= 1.0) {
      //get first column full, m = 0
      for (size_t idl = 1; idl < l - 1; ++idl) {
        Lpolynomial(idl + 2,1) = (1.0/double(idl + 1))*(double(2*idl + 1)*x*Lpolynomial(idl + 1,1) - double(idl)*Lpolynomial(idl,1));
      }
      //get rest of matrix
      double sqrtx = sqrt(1.0 - x*x);
      int ll;
      for (int idm = 1; idm < m; ++idm) {
        for (int idl = idm; idl < l; ++idl) {
          if (idl < l - 1) {
            Lpolynomial(idl + 1,idm + 1) = (1.0/(double(2*idl + 1)*sqrtx))*double(idl - idm + 2)*double(idl - idm + 1)*Lpolynomial(idl + 2,idm);
            Lpolynomial(idl + 1,idm + 1) -= (1.0/(double(2*idl + 1)*sqrtx))*double(idl + idm)*double(idl + idm - 1)*Lpolynomial(idl,idm);
            if (x == 1.0) {Lpolynomial(idl + 1,idm + 1) = 0.0;}
          }
          else {
            ll = idl - 1;
            Lpolynomial(ll + 2,idm + 1) = (1.0/double(ll - idm + 1))*(double(2*ll + 1)*x*Lpolynomial(ll + 1,idm + 1) - double(ll + idm)*Lpolynomial(ll,idm + 1));
          }
        }
      }
    }
    else {
      //get |P_l^m|
      //get first column full, m = 0
      for (size_t idl = 1; idl < l - 1; ++idl) {
        Lpolynomial(idl + 2,1) = (1.0/double(idl + 1))*(double(2*idl + 1)*x*Lpolynomial(idl + 1,1) - double(idl)*Lpolynomial(idl,1));
      }
      //get rest of matrix
      int ll;
      for (int idm = 1; idm < m; ++idm) {
        for (int idl = idm; idl < l; ++idl) {
          if (idl < l - 1) {
            Lpolynomial(idl + 1,idm + 1) = (1.0/double(2*idl + 1))*(Lpolynomial(idl + 2,idm) - Lpolynomial(idl,idm));
          }
          else {
            ll = idl - 1;
            Lpolynomial(ll + 2,idm + 1) = (1.0/double(ll + idm + 1))*(double(2*ll + 1)*x*Lpolynomial(ll + 1,idm + 1) - double(ll - idm)*Lpolynomial(ll,idm + 1));
          }
        }
      }
      //transform |P_l^m| to P_l^m
      double sqrtx = sqrt(x*x - 1.0);
      for (int idl = 0; idl < l; ++idl) {
        for (int idm = 0; idm < m; ++idm) {
          Lpolynomial(idl + 1,idm + 1) *= FactorialRatio(idl + idm,idl - idm)/pow(sqrtx,double(idm));
        }
      }
    }
  }
  return Lpolynomial;
}
double lqmn(size_t m, size_t n, double x) {
  //this function calculates {Q}_{m}^{n}(x) using a python function, thus brute force
  //requirements: python, scipy
  std::string result = "";
  const int max_buffer = 256;
  FILE * stream;
  char buffer[max_buffer];
  std::string command = "python other/lqmn.py "+std::to_string(n)+" "+std::to_string(m)+" "+std::to_string(x);
  stream = popen(command.c_str(), "r");
  while (!feof(stream))
    if (fgets(buffer, max_buffer, stream) != NULL) {result.append(buffer);}
  pclose(stream);
  return std::stod(result);
}
double lpmn(size_t m, size_t n, double x) {
  //test for previous function
  //this function calculates {P}_{m}^{n}(x) using a python function
  //requirements: python, scipy
  std::string result = "";
  const int max_buffer = 256;
  FILE * stream;
  char buffer[max_buffer];
  std::string command = "python other/lpmn.py "+std::to_string(n)+" "+std::to_string(m)+" "+std::to_string(x);
  stream = popen(command.c_str(), "r");
  while (!feof(stream))
    if (fgets(buffer, max_buffer, stream) != NULL) {result.append(buffer);}
  pclose(stream);
  return std::stod(result);
}
matrixE LQmn(size_t l, size_t m, double x) {
  //matrix function calculating recursively all {Q}_{l}^{m}(x) up to given l,m values
  ++l;
  ++m;
  size_t Lmax;
  if (l <= 17) {Lmax = 17;}
  else if (l <= 33) {Lmax = 33;}
  else if (l <= 65) {Lmax = 65;}
  else if (l <= 129) {Lmax = 129;}
  matrixE Lpolynomial(Lmax,m + 1);
  std::vector<double> Qm1(Lmax + 1,0.0);           //vector containing {Q}_{l}^{-1}(x)
  //starting values
  double Q1Lmax = lqmn(Lmax,1,x);
  double Q1Lmaxm1 = lqmn(Lmax-1,1,x);
  //this is required because the algorithm does not work for lower index zero
  double Q10 = lqmn(0,1,x);
  double lmlmp1 = double(Lmax)*double(Lmax + 1);
  double lmlmm1 = double(Lmax)*double(Lmax - 1);
  double sqrtx;
  if (x >= 1.0) {sqrtx = sqrt(x*x - 1.0);}
  else {sqrtx = sqrt(1.0 - x*x);}
  //calculate starting values for {Q}_{l}^{-1}(x)
  Qm1[Lmax] = sqrtx*Q1Lmax/lmlmp1;
  Qm1[Lmax-1] = sqrtx*Q1Lmaxm1/lmlmm1;
  //calculate starting values for {Q}_{l}^{0}(x)
  size_t ll;
  size_t ll1;
  size_t ll2;
  //calculate {Q}_{l}^{-1}(x) and {Q}_{l}^{0}(x)
  for (size_t idx = 2; idx < Qm1.size(); ++idx) {
    ll = Lmax - idx;
    ll1 = Lmax - idx + 1;
    ll2 = Lmax - idx + 2;
    //{Q}_{l}^{-1}(x)
    if (idx != Lmax) {Qm1[ll] = (1.0/double(ll))*(double(2*ll1 + 1)*x*Qm1[ll1] - double(ll2 + 1)*Qm1[ll2]);}
    if (idx + 1 == Qm1.size()) {Qm1[ll] = sqrtx*Q10;}
    //{Q}_{l}^{0}(x)
    if ((ll1 == Lmax - 1)||(ll1 == Lmax - 2)) {
      Lpolynomial(ll1 + 1,1) = (1.0/(double(2*ll1 + 1)*(x*x - 1.0)))*(double((ll2 + 1)*ll2)*Qm1[ll2] - double((ll + 1)*ll)*Qm1[ll]);
    }
    else {
      Lpolynomial(ll1 + 1,1) = Lpolynomial(ll1 + 3,1) - double(2*(ll1 + 1) + 1)*Qm1[ll2];
    }
  }
  Lpolynomial(1,1) = Lpolynomial(3,1) - 3.0*Qm1[1];
  //get the rest of the matrix
  double aux;
  for (size_t icm = 1; icm < m; ++icm) {
    for (size_t icl = icm; icl < Lmax - 1; ++icl) {
      aux = double(icl - icm + 1)*double(icl - icm + 2)*Lpolynomial(icl + 2,icm) - double(icl + icm - 1)*double(icl + icm)*Lpolynomial(icl,icm);
      Lpolynomial(icl + 1,icm + 1) = (1.0/(double(2*icl + 1)*sqrtx))*aux;
    }
    Lpolynomial(Lmax,icm + 1) = lqmn(Lmax-1,icm,x);
  }
  Lpolynomial.resize(l,m);
  return Lpolynomial;
}
std::vector<double> vLQmn(size_t l, size_t m, double x) {
  //vector function calculating recursively all {Q}_{l}^{m}(x) up to given l,m values
  //returns only the last column m
  matrixE Lpolynomial = LQmn(l,m,x);
  std::vector<double> lpol(l + 1,0.0);
  for (size_t id = 0; id < l + 1; ++id) {
    lpol[id] = Lpolynomial(id + 1,m + 1);
  }
  return lpol;
}
matrixE LQmn_HM(size_t l, size_t m, double x) {
  //matrix function calculating recursively all {Q}_{l}^{m}(x) up to given l,m values
  //using the Harris-Michels method
  //
  std::cout << "WARNING: FunctionPackage.hpp: LQmn_HM(): DO NOT USE THIS FUNCTION, unstable algorithm!!!!!!!!!!!!!!!" << std::endl;
  ++l;
  ++m;
  matrixE Lpolynomial(l,m);
  //add extra to the l index to increase quality of calculation
  int lmax = l + 15;
  //allocate matrix r
  matrixE rlm(lmax,m);
  //calculate asymptotic values for rlm
  for (int idm = 0; idm < m; ++idm) {
    rlm(lmax,idm + 1) = (double(lmax - idm + 1)/(double(lmax) + 1.5))*(x - sqrt(x*x - 1.0));
  }
  //back recursion to get full rlm matrix
  for (int idm = 0; idm < m; ++idm) {
    for (int idl = lmax - 2; idl >= idm; --idl) {
      rlm(idl + 1,idm + 1) = double(idl - idm)/(double(2*idl + 1)*x - double(idl + idm + 1)*rlm(idl + 2,idm + 1));
    }
  }
  std::cout << "matrix rlm" << std::endl;
  rlm.Print();
  std::cout << std::endl;
  //calculate |Q_l^m| and convert to Q_l^m
  double sqrtx = sqrt(x*x - 1.0);
  for (int idm = 0; idm < m; ++idm) {
    for (int idl = idm; idl < l; ++idl) {
      if (idl == idm) {Lpolynomial(idl + 1,idm + 1) = (pow(-1.0,double(idm))/DoubleFactorial(2*idl + 1))*(1.0/(x - rlm(idl + 2,idm + 1)))*FactorialRatio(idl + idm,idl - idm)/pow(sqrtx,double(idm));}
      else {Lpolynomial(idl + 1,idm + 1) = rlm(idl + 1,idm + 1)*Lpolynomial(idl,idm + 1)*FactorialRatio(idl + idm,idl - idm)/pow(sqrtx,double(idm));}
    }
  }
  return Lpolynomial;
}
double Pochhammer(double a, size_t n) {
  //function that calculates the Pochhammer symbol (a)_n
  double poch = 1.0;
  //(a)_0 = 1
  if (n > 0) {
    for (size_t idx = 0; idx < n; ++idx) {
      poch *= a + double(idx);
    }
  }
  return poch;
}
double Pochhammer2(double a, size_t n) {
  //function that calculates the Pochhammer symbol (a)_n using the factorial definition; valid for integer a only
  double poch = FactorialRatio(int(a) + n - 1,int(a - 1));
  return poch;
}
double hyper0F1(double l, double x, double threshold = 1e-8, size_t maxiter = 1000) {
  //calculates the hypergeometric series 0F1(l,x) = {\sum}_{k=0}^{\infty}{\frac{{x}^{k}}{{(l)_k}k!}}
  double hyper = 0.0;
  double aux;
  double old;
  double error;
  for (size_t iter = 0; iter < maxiter; ++iter) {
    //calculate x^k/k!
    aux = 1.0;
    old = hyper;
    if (iter > 0) {
      for (size_t idx = 0; idx < iter; ++idx) {
        aux *= x/double(iter - idx);
      }
    }
    hyper += aux/Pochhammer(l,iter);
    error = abs(hyper - old)/old;
    if (fabs(error) < threshold) {break;}
  }
  return hyper;
}
double hyper1F1(double a, double b, double x, double threshold = 1e-8, size_t maxiter = 1000) {
  //calculates the hypergeometric series 1F1(a,b,x) = {\sum}_{k=0}^{\infty}{\frac{{(a)_k}{x}^{k}}{{(b)_k}k!}}
  //this function has exactly the same arguments as python's
  double hyper = 0.0;
  double aux;
  double old;
  double error;
  double poch;
  if ((a == 0.0)&&(b == 0.0)) {hyper = 1.0;}
  else {
    for (size_t iter = 0; iter < maxiter; ++iter) {
      //calculate x^k/k!
      aux = 1.0;
      old = hyper;
      if (iter > 0) {
        for (size_t idx = 0; idx < iter; ++idx) {
          aux *= x/double(iter - idx);
        }
      }
      poch = Pochhammer(a,iter)/Pochhammer(b,iter);
      hyper += aux*poch;
      error = abs(hyper - old)/old;
      if (fabs(error) < threshold) {break;}
    }
  }
  return hyper;
}
std::vector<double> Bessel(double x,size_t l) {
  //function that returns a vector of Bessel functions
  std::vector<double> bessel;
  bessel.push_back(1.0);
  bessel.push_back(1.0 + 1.0/x);
  for (size_t idx = 2; idx < l + 1; ++idx) {
    bessel.push_back(bessel[idx - 2] + ((2.0*double(idx) - 1.0)/x)*bessel[idx - 1]);
  }
  return bessel;
}
std::vector<double> spheBessel_ik(double x, size_t k) {
  //spherical Bessel function ik as defined in Michels and Harris
  //definition of a very large nmax so that accuracy in calculation is high
  size_t nmax = 2*(k + 10);
  std::vector<double> rn(nmax + 1,0.0);
  for (int idx = nmax - 1; idx >= 0; --idx) {
    rn[idx] = -x/(double(2*idx + 1) - x*rn[idx + 1]);
  }
  std::vector<double> ik(k,0.0);
  double ax = fabs(x);
  double ik0 = -cosh(x)/x;
  for (size_t idx = 0; idx < ik.size(); ++idx) {
    if (idx == 0) {ik[idx] = rn[idx]*ik0;}
    else {ik[idx] = rn[idx]*ik[idx - 1];}
  }
  return ik;
}
double Gamma_rec(int alpha, double x, double threshold = 1e-8, size_t maxiter = 1000) {
  //calculation of the incomplete Gamma function of Euler
  //{\Gamma}({{\alpha},x}) = {\int}_{x}^{+\infty}{{e}^{-t}{t}^{{\alpha}-1}dt}
  //calculation via recursion or summation if alpha =< 0
  int alphap = abs(alpha);
  size_t sz = size_t(alphap);
  if (alpha == 0) {++sz;}
  std::vector<double> gamma(sz,0.0);
  gamma[0] = 1.0;
  if (x != 0.0) {
    gamma[0] = exp(-x);
    //get other values
    if (alpha > 0) {
      for (size_t idx = 1; idx < size_t(alpha); ++idx) {
        gamma[idx] = idx*gamma[idx - 1] + pow(x,double(idx))*exp(-x);
      }
    }
    else if (alpha < 0) {
      double aux = - 0.57721566490153286060651209008240243104215933593992 - log(x);
      double oldval = aux;
      for (size_t idk = 1; idk < maxiter; ++idk) {
        aux -= pow(-x,double(idk))/(double(idk)*Factorial(idk));
        if (fabs(aux - oldval) < threshold) {break;}
        oldval = aux;
      }
      gamma[gamma.size() - 1] = pow(-1.0,double(alphap))*aux/Factorial(alphap);
      aux = 0.0;
      for (size_t idk = 0; idk < alphap; ++idk) {
        aux += pow(-1.0,double(idk))*Factorial(alphap - idk - 1)*pow(x,double(idk));
      }
      gamma[gamma.size() - 1] += aux*exp(-x)/(pow(x,double(alphap))*Factorial(alphap));
    }
    else {
      gamma[0] = - 0.57721566490153286060651209008240243104215933593992 - log(x);
      double oldval = gamma[0];
      for (size_t idk = 1; idk < maxiter; ++idk) {
        gamma[0] -= pow(-x,double(idk))/(double(idk)*Factorial(idk));
        if (fabs(gamma[0] - oldval) < threshold) {break;}
        oldval = gamma[0];
      }
    }
  }
  return gamma[gamma.size() - 1];
}
double Gamma_sum(int alpha, double x, double threshold = 1e-8, size_t maxiter = 1000) {
  //calculation of the incomplete Gamma function of Euler
  //{\Gamma}({{\alpha},x}) = {\int}_{x}^{+\infty}{{e}^{-t}{t}^{{\alpha}-1}dt}
  //direct calculation using summation
  double gamma = 1.0;
  if (alpha > 1) {gamma = Factorial(size_t(alpha) - 1);}
  if (alpha != 0.0) {
    double oldval = gamma;
    for (size_t idn = 0; idn < maxiter; ++idn) {
      gamma -= (pow(-1.0,double(idn))*pow(x,double(alpha + idn)))/(Factorial(idn)*double(alpha + idn));
      if (fabs(gamma - oldval) < threshold) {break;}
      oldval = gamma;
    }
  }
  return gamma;
}
template <class type> double DeltaDirac(type x,type y) {
  //function that calculates the Dirac delta function
  return double(x == y);
}
double Clebsch_Gordan(double j1, double m1, double j2, double m2, double J, double M) {
  //function that calculates the Clebsch-Gordan coefficients
  //<j1m1;j2m2|JM>
  //note that the function can take half integer indices, if they sum up to int
  double cgcoeff = 0.0;
  if ((m1 + m2) == M) {
    int s = int(j1 + j2 + J);
    int ia = s - 2*J;
    int ib = j1 - m1;
    int ic = j2 + m2;
    int id = J - j2 + m1;
    int ie = J - j1 - m2;
    int maxv = ia;
    if (ib < maxv) {maxv = ib;}
    if (ic < maxv) {maxv = ic;}
    int minv = 0;
    if ((-id) > minv) {minv = -id;}
    if ((-ie) > minv) {minv = -ie;}
    for (int iv = minv; iv < (maxv + 1); ++iv) {
      cgcoeff += pow(-1.0,double(iv))*FactorialRatio(int(j1 - m1),ib - iv)*FactorialRatio(int(j2 + m2),ic - iv)/(Factorial(iv)*Factorial(ia - iv)*Factorial(id + iv)*Factorial(ie + iv));
    }
    cgcoeff *= sqrt(double(int(2*J + 1)*Factorial(ia)*FactorialRatio(int(s - 2*j1),s + 1)*Factorial(int(s - 2*j2))));
    cgcoeff *= sqrt(double(FactorialRatio(int(j1 + m1),int(j1 - m1))*FactorialRatio(int(j2 - m2),int(j2 + m2))*Factorial(int(J + M))*Factorial(int(J - M))));
  }
  return cgcoeff;
}
double HarmonicNumber(int n) {
  //function calculating harmonic numbers {H}_{n} = {\sum}_{k=1}^{n}{\frac{1}{k}}
  double hn = 0.0;
  for (size_t idx = 1; idx < n + 1; ++idx) {
    hn += 1.0/double(idx);
  }
  return hn;
}
double BernoulliPolynomialMatrix(double x,int n) {
  //function returning a single value for a Bernoulli polynomial
  //calculation via matrix determinant as defined in T. Agoh INTEGERS 19 (2019)
  matrixE A = Identity(n + 1);
  for (int irow = 0; irow < A.rows(); ++irow) {
    for (int icol = 0; icol < irow; ++icol) {
      A(irow + 1,icol + 1) = 1.0/Factorial(irow + 1 - icol);
    }
  }
  for (size_t idx = 0; idx < A.rows(); ++idx) {
    A(idx + 1,n + 1) = pow(x,double(idx))/Factorial(idx);
  }
  double bernoulli = Factorial(n)*A.det();
  return bernoulli;
}
double zeta(double x, size_t maxiter = 1000, double threshold = 1e-7, bool shortcut = true) {
  //function for Riemann's zeta function
  //use tabulated values from wolfram preferred; if not possible then calculate by summing series
  double rzeta = 0.0;
  bool shortcutused = false;
  if (shortcut) {
    if (x == 0.0) {
      rzeta = -0.5;
      shortcutused = true;
    }
    else if (x == 2.0) {
      rzeta = pi*pi/6.0;
      shortcutused = true;
    }
    else if (x == 3.0) {
      rzeta = 1.2020569032;
      shortcutused = true;
    }
    else if (x == 4.0) {
      rzeta = pi*pi*pi*pi/90.0;
      shortcutused = true;
    }
    else if (x == 5.0) {
      rzeta = 1.0369277551;
      shortcutused = true;
    }
    else if (x == 6.0) {
      rzeta = pi*pi*pi*pi*pi*pi/945.0;
      shortcutused = true;
    }
    else if (x == 7.0) {
      rzeta = 1.0083492774;
      shortcutused = true;
    }
    else if (x == 8.0) {
      rzeta = pi*pi*pi*pi*pi*pi*pi*pi/9450.0;
      shortcutused = true;
    }
    else if (x == 9.0) {
      rzeta = 1.0020083928;
      shortcutused = true;
    }
    else if (x == 10.0) {
      rzeta = pi*pi*pi*pi*pi*pi*pi*pi*pi*pi/93555.0;
      shortcutused = true;
    }
  }
  if ((!shortcut)||(!shortcutused)) {
    //check whether x is negative int
    double absolute = fabs(x);
    if (absolute == floor(absolute)) {
      //if here, then x is int
      if (x < 0.0) {shortcutused = true;}
    }
    if (shortcutused) {
      //shortcut for x negative integer
      rzeta = pow(-1.0,absolute)*BernoulliPolynomialMatrix(0.0,-int(x) + 1)/(-x + 1.0);
    }
    else {
      //calculate by definition
      double toadd;
      for (size_t iter = 1; iter < maxiter; ++iter) {
        toadd = pow(iter,-x);
        rzeta += toadd;
        if (toadd < threshold) {break;}
      }
    }
  }
  return rzeta;
}
double BernoulliNumber(int n) {
  //function returning Bernoulli numbers using Euler's formula
  double bernoulli = 1.0;
  if (n == 1) {bernoulli = -0.5;}
  else if (n % 2 != 0) {bernoulli = 0.0;}
  else {
    //use Euler formula
    if (n > 0) {bernoulli = 2.0*pow(-1.0,double(n)/2.0 + 1.0)*Factorial(n)*zeta(n)/pow(2.0*pi,double(n));}
  }
  return bernoulli;
}
double BernoulliNumberS(int n) {
  //function returning Bernoulli numbers using Euler's formula
  double bernoulli = 1.0;
  if (n == 1) {bernoulli = -0.5;}
  else if (n % 2 != 0) {bernoulli = 0.0;}
  else {
    //use Euler formula
    if (n > 0) {bernoulli = 2.0*pow(-1.0,double(n)/2.0 + 1.0)*Factorial_Stirling(n)*zeta(n)/pow(2.0*pi,double(n));}
  }
  return bernoulli;
}
double BernoulliPolynomial(double x,int n) {
  //function returning a single value for a Bernoulli polynomial
  double bernoulli = 0.0;
  if (n == 0) {bernoulli = 1.0;}
  else if (n == 1) {bernoulli = x - 0.5;}
  else if (n == 2) {bernoulli = x*x - x + 1.0/6.0;}
  else if (n == 3) {bernoulli = x*x*x - 1.5*x*x + 0.5*x;}
  else if (n == 4) {bernoulli = x*x*x*x - 2.0*x*x*x + x*x - 1.0/30.0;}
  else if (n == 5) {bernoulli = x*x*x*x*x - 2.5*x*x*x*x + 5.0*x*x*x/3.0 - x/6.0;}
  else if (n == 6) {bernoulli = x*x*x*x*x*x - 3.0*x*x*x*x*x + 2.5*x*x*x*x - 0.5*x*x + 1.0/42.0;}
  else {
    bernoulli = BernoulliPolynomialMatrix(x,n);
  }
  return bernoulli;
}
std::complex<double> BernoulliPolynomial(std::complex<double> z,int n) {
  //function returning a single value for a Bernoulli polynomial
  std::complex<double> bernoulli(0.0,0.0);
  if (n == 0) {bernoulli.real(1.0);}
  else if (n == 1) {bernoulli = z - 0.5;}
  else if (n == 2) {bernoulli = z*z - z + 1.0/6.0;}
  else if (n == 3) {bernoulli = z*z*z - z*z*1.5 + z*0.5;}
  else if (n == 4) {bernoulli = z*z*z*z - z*z*z*2.0 + z*z - 1.0/30.0;}
  else if (n == 5) {bernoulli = z*z*z*z*z - z*z*z*z*2.5 + z*z*z*5.0/3.0 - z/6.0;}
  else if (n == 6) {bernoulli = z*z*z*z*z*z - z*z*z*z*z*3.0 + z*z*z*z*2.5 - z*z*0.5 + 1.0/42.0;}
  return bernoulli;
}
double Lis_python(int s, double x) {
  //test for next function
  //this function calculates {Li}_{s}(x) using a python function
  //requirements: python, mpmath
  std::string result = "";
  const int max_buffer = 256;
  FILE * stream;
  char buffer[max_buffer];
  std::string command = "python other/polylog.py "+std::to_string(s)+" "+std::to_string(x);
  stream = popen(command.c_str(), "r");
  while (!feof(stream))
    if (fgets(buffer, max_buffer, stream) != NULL) {result.append(buffer);}
  pclose(stream);
  return std::stod(result);
}
double Lis(int s,double z,double threshold = 1e-7,int maxiter = 500) {
  //function that calculates polylogarithms of real numbers.
  //note that result is always real as well
  //from definition,
  //{Li}_{s}(z)={\sum}_{k=1}^{+\infty}{\frac{{z}^{k}}{{k}^{s}}}
  double lis = 0.0;
  if (fabs(z - 1.0) < threshold) {
    lis = zeta(s,maxiter,threshold);
    if (s == 1) {std::cout << std::endl << "WARNING: FunctionPackage.hpp: Lis(): polylogarithm index 1 of 1 is a pole." << std::endl;}
  }
  else if (fabs(z + 1.0) < threshold) {
    lis = -(1.0 - pow(2.0,1-s))*zeta(s,maxiter,threshold);
    if (s == 1) {lis = -log(2.0);}
  }
  else {
    if (s == 1) {lis = -log(1.0 - z);}
    else if (s == 0) {
      lis = z/(1.0 - z);
    }
    else if (s == -1) {lis = z/((1.0 - z)*(1.0 - z));}
    else if (s == -2) {lis = z*(1.0 + z)/((1.0 - z)*(1.0 - z)*(1.0 - z));}
    else if (s == -3) {lis = z*(1.0 + 4.0*z + z*z)/((1.0 - z)*(1.0 - z)*(1.0 - z)*(1.0 - z));}
    else if (s == -4) {lis = z*(1.0 + z)*(1.0 + 10.0*z + z*z)/((1.0 - z)*(1.0 - z)*(1.0 - z)*(1.0 - z)*(1.0 - z));}
    else {
      if ((z <= 0.5)&&(z >= -0.5)) {
        double num = 1.0;
        double oldval = 0.0;
        for (size_t iter = 1; iter < maxiter; ++iter) {
          num *= z;
          lis += num/pow(double(iter),double(s));
          if (fabs(oldval - lis) < threshold) {break;}
          oldval = lis;
        }
      }
      else if ((z <= -2.0)||(z >= 2.0)) {
        //the proper calculation of these terms involve either a complex Lis or, if the result is real, it requires the calculation of logarithms of negative numbers, which for the time being I cannot get
        if (z < 0) {lis = Lis_python(s,z);}
        else {throw(std::string("ERROR: FunctionPackage.hpp: Lis(): complex result for z > 0"));}
      }
      else {
        if (z <= 0.0) {throw(std::string("ERROR: FunctionPackage.hpp: Lis(): complex result for z <= 0"));}
        if (s > 0) {
          if (z >= 1.0) {throw(std::string("ERROR: FunctionPackage.hpp: Lis(): complex result for z >= 1"));}
          double num = 1.0;
          double oldval = 0.0;
          for (size_t iter = 1; iter < maxiter; ++iter) {
            num *= z;
            lis += num/pow(double(iter),double(s));
            if (fabs(oldval - lis) < threshold) {break;}
            oldval = lis;
          }
        }
        else {
          double oldval = 0.0;
          for (size_t k = 0; k < maxiter; ++k) {
            lis -= BernoulliNumber(k-s+1)*pow(log(z),double(k))/(Factorial(k)*double(k-s+1));
            if (fabs(oldval - lis) < threshold) {break;}
            oldval = lis;
          }
          lis += Factorial(-s)*pow(-log(z),double(s-1));
        }
      }
    }
  }
  return lis;
}
double trunc_en(size_t n, double x) {
  //truncated exponential function on order n
  double en = 0.0;
  double factorial;
  for (size_t idx = 0; idx < (n + 1); ++idx) {
    if (idx > 10) {factorial = Factorial_Stirling(idx);}
    else {factorial = Factorial(idx);}
    en += pow(x,double(idx))/factorial;
  }
  return en;
}
matrixE StirlingNumers2(int im, int ik) {
  //function to calculate Stirling numbers of second kind
  //m is for rows, k is for columns
  int minrc = std::min(im,ik);
  matrixE SN2 = Identity(minrc + 1);
  int maxrc = std::max(im,ik);
  if (maxrc != minrc) {SN2.resize(im + 1,ik + 1);}
  size_t nrows = SN2.rows();
  size_t ncols = SN2.cols();
  for (int idr = 1; idr < nrows; ++idr) {
    SN2(idr + 1,2) = 1.0;
    for (int idk = 2; idk < ncols; ++idk) {
      SN2(idr + 1,idk + 1) = SN2(idr,idk) + idk*SN2(idr,idk + 1);
    }
  }
  return SN2;
}
double SphericalHarmonic(std::vector<double> pos, double r, int l, int m) {
  //function evaluating unnormalized real spherical harmonics up to l = 2
  //order for quantum number m is arbitrary, just to somehow follow the "rule" for complex spherical harmonics
  double rsphharm = 0.0;
  if (l == 0) {rsphharm = 1.0;}
  else if (l == 1) {rsphharm = pos[m + 1]/r;}                                                     //x = 0, y = 1, z = 2
  else if (l == 2) {
    double r2 = r*r;
    if (m == -2) {rsphharm = pos[0]*pos[2]/r2;}                                                   //xz
    else if (m == -1) {rsphharm = pos[1]*pos[2]/r2;}                                              //yz
    else if (m == 0) {rsphharm = (2.0*pos[2]*pos[2] - pos[1]*pos[1] - pos[0]*pos[0])/r2;}         //z2
    else if (m == 1) {rsphharm = pos[0]*pos[1]/r2;}                                               //xy
    else if (m == 2) {rsphharm = (pos[0]*pos[0] - pos[1]*pos[1])/r2;}                             //x2-y2
  }
  return rsphharm;
}

#endif //_Function_Package_
