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

#ifndef _Integral_Package_
#define _Integral_Package_
#include <math.h>
#include <cstdlib>
#include <vector>
#include <string>
#include "FunctionPackage.hpp"
#include "SpecialFunctionPackage.hpp"
#include "../ConstantsPackage.hpp"

//description:
//package for the calculation of integrals

//auxiliary functions
std::vector<double> Ck(double x, size_t k, double tolerance = 1e-5) {
  //function that calculates a series relevant for the integrals
  //{\int}_{-1}^{+1}{{\eta}^{k}{e}^{-p{\eta}}d{\eta}}
  std::vector<double> result(k+1,0.0);
  if (fabs(x) > tolerance) {        //ensure that x != 0
    result[0] = -exp(x)/x;
    for (size_t idx = 1; idx < result.size(); ++idx) {
      result[idx] = -result[idx-1];
      for (size_t idx2 = 0; idx2 < idx; ++idx2) {
        result[idx] += FactorialRatio(idx-1,idx2)*result[idx2]/(pow(x,idx-idx2));
      }
    }
  }
  return result;
}
//integral functions
double int_1_infty_xpk_expmyx_dx_direct(double x, size_t k) {
  //function that calculates directly the integrals
  //{\int}_{1}^{+{\infty}}{{\xi}^{k}{e}^{-p{\xi}}d{\xi}} (a.k.a. Ak integrals)
  //all the way up to order k, returning each integral as an entry in a vector
  double xx = 1.0;
  double factor = exp(-x)/pow(x,k + 1);
  double result = 0.0;
  for (size_t idj = 0; idj < k + 1; ++idj) {
    result += xx*FactorialRatio(k,idj);
    xx *= x;
  }
  return result*factor;
}
std::vector<double> int_1_infty_xpk_expmyx_dx(double x, size_t k, double tolerance = 1e-5) {
  //function that calculates the integrals
  //{\int}_{1}^{+{\infty}}{{\xi}^{k}{e}^{-p{\xi}}d{\xi}} (a.k.a. Ak integrals)
  //all the way up to order k, returning each integral as an entry in a vector
  std::vector<double> result(k+1,0.0);
  double xx = x;
  if (x == 0.0) {xx = 0.001;}
  if (fabs(xx) > tolerance) {        //ensure that x != 0
    result[0] = exp(-xx)/xx;
    for (size_t idx = 1; idx < result.size(); ++idx) {
      if (idx < 76) {
        result[idx] = result[idx-1];
        for (size_t idx2 = 0; idx2 < idx; ++idx2) {
          result[idx] += FactorialRatio(idx-1,idx2)*result[idx2]/(pow(xx,idx-idx2));
        }
      }
      else {result[idx] = int_1_infty_xpk_expmyx_dx_direct(x,idx);}
    }
  }
  return result;
}
double int_0_1_xpk_expmyx_dx_direct(double x, size_t k, double tolerance = 1e-25, int maxiter = 500) {
  //function that calculates directly the integrals
  //{\int}_{0}^{1}{{\xi}^{k}{e}^{-p{\xi}}d{\xi}} (a.k.a. ak integrals), which are complementary 
  double result = 0.0;
  double aux;
  for (int idj = 0; idj < maxiter; ++idj) {
    aux = pow(x,idj)/FactorialRatio(k + idj + 1,k);
    result += aux;
    if (fabs(aux) < tolerance) {break;}
  }
  double expx = exp(-x);
  return result*expx;
}
std::vector<double> int_0_1_xpk_expmyx_dx(double x, size_t k, double tolerance = 1e-12) {
  //function that calculates the integrals
  //{\int}_{0}^{1}{{\xi}^{k}{e}^{-p{\xi}}d{\xi}} (a.k.a. ak integrals), which are complementary 
  //all the way up to order k, returning each integral as an entry in a vector
  std::vector<double> rk(3*k + 10,0.0);
  rk[rk.size() - 1] = 1.0;
  int nn;
  for (int idr = 1; idr < rk.size(); ++idr) {
    nn = rk.size() - 1 - idr;
    rk[nn] = double(nn + 1)/(double(nn + 2) + x*(1.0 - rk[nn + 1]));
  }
  std::vector<double> result(k + 1,0.0);
  result[0] = (1.0 - exp(-x))/x;
  for (int idr = 1; idr < result.size(); ++idr) {
    result[idr] = result[idr - 1]*rk[idr - 1];
  }
  return result;
}
double int_1_infty_xpm1_expmyx_dx(double x, int maxiter = 500, double threshold = 1.0e-16) {
  //analytical solution to {\int}_{1}^{+{\infty}}{\frac{{e}^{-px}}{x}dx} (a.k.a. E1 integrals)
  double result = -emgamma - log(x);
  double aux;
  for (size_t idj = 1; idj < maxiter; ++idj) {
    aux = pow(-x,idj)/(idj*Factorial(idj));
    result -= aux;
    if (fabs(aux) < threshold*threshold) {break;}
  }
  return result;
}
std::vector<double> int_1_infty_xpmk_expmyx_dx(double x, size_t k) {
  //function that calculates the integrals
  //{\int}_{1}^{+{\infty}}{{x}^{-k}{e}^{-px}dx} (a.k.a. Ek integrals)
  //all the way up to order k, returning each integral as an entry in a vector
  std::vector<double> en(k+1,0.0);
  en[0] = exp(-x)/x;
  en[1] = int_1_infty_xpm1_expmyx_dx(x);
  if (k > 1) {
    for (size_t idx = 2; idx < k + 1; ++idx) {
      en[idx] = (exp(-x) - x*en[idx - 1])/double(idx - 1);
    }
  }
  return en;
}
std::vector<double> int_1_infty_xpmk_expmyx_dx_2(double x, size_t k) {
  //alternative function that calculates the integrals
  //{\int}_{1}^{+{\infty}}{{x}^{-k}{e}^{-px}dx} (a.k.a. Ek integrals)
  //all the way up to order k, returning each integral as an entry in a vector
  std::vector<double> en(k+1,0.0);
  en[0] = exp(-x)/x;
  en[1] = int_1_infty_xpm1_expmyx_dx(x);
  if (k > 1) {
    for (size_t idx = 2; idx < k + 1; ++idx) {
      en[idx] = pow(-x,double(idx - 1))*en[1]/Factorial(idx - 1);
      for (size_t idj = 1; idj < idx; ++idj) {
        en[idx] += exp(-x)*pow(-x, double(idj - 1))*FactorialRatio(idx - 1 - idj,idx - 1);
      }
    }
  }
  return en;
}
std::vector<double> int_m1_1_xpk_expmyx_dx(double x, size_t k, double tolerance = 1e-5) {
  //function that calculates the integrals
  //{\int}_{1}^{+{\infty}}{{\xi}^{k}{e}^{-p{\xi}}d{\xi}} (a.k.a. Bk integrals)
  //all the way up to order k, returning each integral as an entry in a vector
  std::vector<double> result(k+1,0.0);
  std::vector<double> ak = int_1_infty_xpk_expmyx_dx(x,k,tolerance);
  std::vector<double> ck = Ck(x,k,tolerance);
  if (fabs(x) > tolerance) {        //ensure that x != 0
    for (size_t idx = 0; idx < result.size(); ++idx) {
      result[idx] = -ak[idx] - ck[idx];
    }
  }
  else {                            //if x == 0
    for (size_t idx = 0; idx < result.size(); ++idx) {
      if (idx%2 == 0) {result[idx] = 2.0/(idx+1.0);}
    }
  }
  return result;
}
std::vector<double> int_m1_1_xpk_expmyx_dx_Pople(double x, size_t k, double tolerance = 1e-5) {
  //function that calculates the integrals
  //{\int}_{1}^{+{\infty}}{{\xi}^{k}{e}^{-p{\xi}}d{\xi}} (a.k.a. Bk integrals)
  //all the way up to order k sccording to the Pople formulae gaven in his book
  // "Approximate Molecular Orbital theory". The result is a vector containing the integrals
  std::vector<double> result(k+1,0.0);
  if (fabs(x) <= tolerance) {                             //case |x| <= tolerance
    for (size_t idx = 0; idx < result.size(); ++idx) {
      if (idx%2 == 0) {result[idx] = 2.0/(idx+1.0);}
    }
  }
  else if ((fabs(x) > tolerance)&&(fabs(x) <= 0.5)) {     //case tolerance < |x| <= 0.5
    size_t last = 6;
    double y = 0.0;
    for (size_t idx = 0; idx < result.size(); ++idx) {
      y = 0.0;
      for (size_t idx2 = 0; idx2 < last; ++idx2) {
        y += (pow(-x,idx2)*(1.0 - pow(-1.0,idx+idx2+1)))/(double(Factorial(idx2)*(idx+idx2+1)));
      }
      result[idx] = y;
    }
  }
  else if ((fabs(x) > 0.5)&&(fabs(x) <= 1.0)) {           //case 0.5 < |x| <= 1.0
    if (k <= 5) {
      double expx = exp(x);
      double expmx = 1.0/expx;
      result[0] = (expx - expmx)/x;
      for (size_t idx = 0; idx < (result.size()-1); ++idx) {
        result[idx+1] = (double(idx+1)*result[idx] + pow(-1.0,idx+1)*expx - expmx)/x;
      }
    }
    else {
      size_t last = 7;
      double y = 0.0;
      for (size_t idx = 0; idx < result.size(); ++idx) {
        y = 0.0;
        for (size_t idx2 = 0; idx2 < last; ++idx2) {
          y += (pow(-x,idx2)*(1.0 - pow(-1.0,idx+idx2+1)))/(double(Factorial(idx2)*(idx+idx2+1)));
        }
        result[idx] = y;
      }
    }
  }
  else if ((fabs(x) > 1.0)&&(fabs(x) <= 2.0)) {           //case 1.0 < |x| <= 2.0
    if (k <= 7) {
      double expx = exp(x);
      double expmx = 1.0/expx;
      result[0] = (expx - expmx)/x;
      for (size_t idx = 0; idx < (result.size()-1); ++idx) {
        result[idx+1] = (double(idx+1)*result[idx] + pow(-1.0,idx+1)*expx - expmx)/x;
      }
    }
    else {
      size_t last = 12;
      double y = 0.0;
      for (size_t idx = 0; idx < result.size(); ++idx) {
        y = 0.0;
        for (size_t idx2 = 0; idx2 < last; ++idx2) {
          y += (pow(-x,idx2)*(1.0 - pow(-1.0,idx+idx2+1)))/(double(Factorial(idx2)*(idx+idx2+1)));
        }
        result[idx] = y;
      }
    }
  }
  else if ((fabs(x) > 2.0)&&(fabs(x) <= 3.0)) {           //case 2.0 < |x| <= 3.0
    if (k <= 10) {
      double expx = exp(x);
      double expmx = 1.0/expx;
      result[0] = (expx - expmx)/x;
      for (size_t idx = 0; idx < (result.size()-1); ++idx) {
        result[idx+1] = (double(idx+1)*result[idx] + pow(-1.0,idx+1)*expx - expmx)/x;
      }
    }
    else {
      size_t last = 15;
      double y = 0.0;
      for (size_t idx = 0; idx < result.size(); ++idx) {
        y = 0.0;
        for (size_t idx2 = 0; idx2 < last; ++idx2) {
          y += (pow(-x,idx2)*(1.0 - pow(-1.0,idx+idx2+1)))/(double(Factorial(idx2)*(idx+idx2+1)));
        }
        result[idx] = y;
      }
    }
  }
  else if (fabs(x) > 3.0) {                               //case 3.0 < |x|
    double expx = exp(x);
    double expmx = 1.0/expx;
    result[0] = (expx - expmx)/x;
    for (size_t idx = 0; idx < (result.size()-1); ++idx) {
      result[idx+1] = (double(idx+1)*result[idx] + pow(-1.0,idx+1)*expx - expmx)/x;
    }
  }
  return result;
}
double MullikenBExponent(double R, double a, double b, size_t na, size_t nb) {
  //function that calculates the exponent for integrals Bk according to Mulliken's rules defined in 
  //R. S. Mulliken, C. A. Rieke, D. Orloff, H. Orloff, J. Chem. Phys., 17(12), 1248, 1949
  double result = 0.0;
  if (na < nb) {result = 0.5*R*(a - b);}
  else if (na > nb) {result = 0.5*R*(b - a);}
  else {
    if (a > b) {result = 0.5*R*(a - b);}
    else {result = 0.5*R*(b - a);}
  }
  return result;
}
double GauntIntegral(size_t l, size_t u, size_t m, size_t v, size_t n, size_t w) {
  //function that calculates the Gaunt integral over 3 associated Legendre polynomials of kind 1
  //{\int}_{-1}^{1}{{P}_{l}^{u}{(x)}{P}_{m}^{v}{(x)}{P}_{n}^{w}{(x)}dx}
  //order or arguments follows above declaration
  if (int(u) != (int(v) + int(w))) {throw std::string("ERROR: IntegralPackage.hpp: GauntIntegral(): u != v+w");}
  double integral = 0.0;
  int ds = int(l) + int(m) + int(n);
  if ((ds % 2) == 0) {
    //then s is even and integral might not be zero
    int m_plus_n = int(m) + int(n);
    int m_minus_n = int(m) - int(n);
    if ((l <= m_plus_n)&&(m_minus_n <= l)) {
      //only in these conditions the integral does not zero
      //calculate variables involved in integral
      size_t s = size_t(ds/2.0);
      size_t p = 0;
      int aux = int(n) - int(m) - int(u);
      if (aux > 0) {p = size_t(aux);}
      size_t q = 0;
      aux = int(m) + int(n) - int(u);
      if (aux > 0) {q = size_t(aux);}
      aux = int(l) - int(u);
      if ((aux >= 0)&&(size_t(aux) < q)) {q = size_t(aux);}
      aux = int(n) - int(w);
      if ((aux >= 0)&&(size_t(aux) < q)) {q = size_t(aux);}
      //calculate auxiliaries
      double aux1 = 0.0;
      double aux2 = 0.0;
      double aux3 = 0.0;
      //calculate integral
      integral = 2.0*pow(-1.0,double(s - m - w))*FactorialRatio(m + v,m - v)*FactorialRatio(2*(s - n),s - n);
      size_t auxmax = std::max(s - l,s - m);
      size_t auxmin = std::min(s - l,s - m);
      integral *= FactorialRatio(s,2*s + 1)*FactorialRatio(n + w,auxmax)/Factorial(auxmin);
      double summation = 0.0;
      for (size_t t = p; t < q + 1; ++t) {
        aux1 = pow(-1.0,double(t))/(Factorial(t)*Factorial(n - w - t));
        aux2 = FactorialRatio(l + u + t,l - u - t);
        aux3 = FactorialRatio(m + n - u - t,m - n + u + t);
        summation += aux1*aux2*aux3;
      }
      integral *= summation;
    }
  }
  return integral;
}
double debye_type1(double x, size_t k, double tolerance = 1e-7, size_t maxiter = 500) {
  //function that calculates the integrals
  //{\int}_{0}^{x}{\frac{{t}^{k+1}{e}^{t}}{{({e}^{t} - 1)}^{2}}dt} 
  //these are the integrals required for the calculation of heat capacities at constant volume
  double integral = 0.0;
  double factkp1 = Factorial(k + 1);
  double powx = 1.0/x;
  double expmx = exp(-x);
  double oldval = 0.0;
  double polylog;
  for (size_t idx = 0; idx < k + 2; ++idx) {
    powx *= x;
    polylog = Lis(k + 1 - int(idx),expmx,tolerance,maxiter);
    integral -= factkp1*powx*polylog/Factorial(idx);
    if (fabs(oldval - integral) < tolerance) {break;}
    oldval = integral;
  }
  integral += factkp1*zeta(double(k) + 1.0,maxiter,tolerance);
  return integral;
}
double debye_type2(double x, size_t k, double tolerance = 1e-7, size_t maxiter = 500) {
  //function that calculates the integrals
  //{\int}_{0}^{x}{\frac{{t}^{k}{e}^{-t}}{1 - {e}^{-t}}dt} 
  //these are the integrals required for the calculation of internal energies
  double integral = pow(-1.0,double(k))*Factorial(k)*zeta(double(k) + 1.0,maxiter,tolerance);
  double expx = exp(x);
  double powx = 1.0/x;
  double factratio;
  double polylog;
  std::cout << std::endl << expx << std::endl;
  for (int m = 0; m < k + 1; ++m) {
    powx *= x;
    factratio = FactorialRatio(k,m);
    std::cout << Lis(k + 1 - m,expx,tolerance,maxiter) << std::endl;
    integral += pow(-1.0,k - m + 1)*factratio*powx*Lis(k + 1 - m,expx,tolerance,maxiter) - pow(x,double(k + 1))/double(k + 1);
  }
  return integral;
}
double int_0_a_xpmlog1mexpmx(double x, size_t m, double tolerance = 1e-7, size_t maxiter = 500) {
  //function that calculates the integrals
  //{\int}_{0}^{a}{{x}^{m}log(1.0 - {e}^{-x})dx} 
  //these are the integrals required for the calculation of heat capacities at constant volume
  double integral = 0.0;
  std::vector<double> Am;
  double factm = Factorial(m);
  double apmp1 = pow(x,double(m + 1));
  double oldval = 0.0;
  double newterm;
  for (int n = 1; n < maxiter; ++n) {
    Am = int_1_infty_xpk_expmyx_dx(double(n)*x,m,tolerance);
    newterm = pow(x,double(m + 1))*Am[m + 1] - factm/pow(double(n),double(m + 1));
    integral += newterm/double(n);
    if (fabs(integral - oldval) < tolerance) {break;}
    oldval = integral;
  }
  return integral;
}
double Shi(double x, int maxiter = 500, double threshold = 1.0e-16) {
  //function calculating the hyperbolic sine integral function by deifiniton
  //{\int}_{0}^{x}{\frac{sinht}{t}dt}
  double shix = 0.0;
  double aux;
  double factorial;
  for (size_t iter = 0; iter < maxiter; ++iter) {
    if (iter < 8) {factorial = 1.0/Factorial(2*iter + 1);}
    else {factorial = 1.0/Factorial_Stirling(2*iter + 1);}
    aux = pow(x,2*iter + 1)*factorial/double(2*iter + 1);
    shix += aux;
    if (fabs(aux) < threshold) {break;}
  }
  return shix;
}
double Chi(double x, int maxiter = 500, double threshold = 1.0e-16) {
  //function calculating the hyperbolic sine integral function by deifiniton
  //{\int}_{0}^{x}{\frac{sinht}{t}dt}
  double chix = emgamma + log(x);
  double aux;
  double factorial;
  for (size_t iter = 1; iter < maxiter; ++iter) {
    if (iter < 8) {factorial = 1.0/Factorial(2*iter);}
    else {factorial = 1.0/Factorial_Stirling(2*iter);}
    aux = pow(x,2*iter)*factorial/double(2*iter);
    chix += aux;
    if (fabs(aux) < threshold) {break;}
  }
  return chix;
}
std::vector<double> Enodd(double x, size_t k) {
  //function that only retrieves the first k odd elements in the En series of integrals, discarding En integrals with even n
  //Note that these integrals are defined in IntegralPackageSpecial
  std::vector<double> entotal = int_1_infty_xpmk_expmyx_dx(x,2*k);
  std::vector<double> enodd;
  for (size_t idx = 0; idx < k; ++idx) {
    enodd.push_back(entotal[2*idx + 1]);
  }
  return enodd;
}
double Imsqb_direct(int mu, int sigma, int qq, double beta, double threshold = 1.0e-25, int maxiter = 500) {
  //function for direct calculation of the integrals {i}^{\sigma}_{\mu}(q,{\beta}) = \frac{{(-1)}^{\sigma}}{2}\frac{({\sigma}-{\mu})!}{({\sigma}+{\mu})!}{\int}_{-1}^{+1}{{e}^{-{\beta}t}{(1-{t}^{2})}^{{\sigma}/2}{t}^{q}{P}_{\mu}^{\sigma}(t)dt}
  //calculation via series as described in
  //F. E. Harris, Int. J. Quantum Chem., 88, 701, 2002
  double imsqb = 0.0;
  if (mu >= sigma) {
    double aux;
    double maxval = 0.0;
    for (int iter = 0; iter < maxiter; ++iter) {
      if (mu - sigma - qq + 2*iter < 0) {continue;}
      aux = pow(-1.0,double(qq))*pow(beta,double(mu - sigma + 2*iter - qq))*FactorialRatio(mu - sigma + 2*iter,mu - sigma + 2*iter - qq)/(DoubleFactorial2(2*(iter + mu) + 1)*DoubleFactorial2(2*iter));
      if (fabs(aux) > maxval) {maxval = fabs(aux);}
      imsqb += aux;
      if (fabs(aux)/maxval < threshold) {break;}
    }
  }
  return imsqb;
}
void Imsqb(matrixE & imsqb, int sigma, double beta, size_t Gamma, size_t mumax = 50, double threshold = 1.0e-25, int maxiter = 500) {
  //function calculating a bunch of integrals {i}^{\sigma}_{\mu}(q,{\beta}) = \frac{{(-1)}^{\sigma}}{2}\frac{({\sigma}-{\mu})!}{({\sigma}+{\mu})!}{\int}_{-1}^{+1}{{e}^{-{\beta}t}{(1-{t}^{2})}^{{\sigma}/2}{t}^{q}{P}_{\mu}^{\sigma}(t)dt}
  //procedure as described in
  //F. E. Harris, Int. J. Quantum Chem., 88, 701, 2002
  //allocating integrals as matrix {i}^{\sigma} with elements {\mu} (rows) and q (columns); {\beta} is simply argument to integrals
  mumax += 2*sigma + 2;
  std::vector<double> imumax_0(mumax + 1, 0.0);
  std::vector<double> imumax_s(mumax + 1, 0.0);
  std::vector<double> rmu(2*mumax + 1, 0.0);
  double aux;
  double val;
  rmu[mumax] = 0.0;
  //get the rmu's
  for (int idmu = 2*mumax; idmu > 0; --idmu) {
    rmu[idmu - 1] = beta/(double(2*idmu + 1) + beta*rmu[idmu]);
  }
  imumax_s[0] = Imsqb_direct(0,0,0,beta,threshold,maxiter);
  //get {i}^{0}_{\mu}(0,{\beta})
  for (int idmu = 1; idmu < imumax_s.size(); ++idmu) {
    imumax_s[idmu] = imumax_s[idmu - 1]*rmu[idmu - 1];
  }
  imumax_0 = imumax_s;
  //get {i}^{\sigma}_{\mu}(0,{\beta})
  for (int idsigma = 1; idsigma < sigma + 1; ++idsigma) {
    for (int idmu = idsigma; idmu < mumax - 1; ++idmu) {
      imumax_s[idmu] = (imumax_0[idmu - 1] - imumax_0[idmu + 1])/double(2*idmu + 1);
    }
    imumax_0 = imumax_s;
  }
  mumax -= 2*sigma + 2;
  imsqb.resize(mumax + 1,Gamma + 1);
  imsqb.zero();
  //get first column in place
  for (size_t idrow = sigma; idrow < mumax + 1; ++idrow) {
    imsqb(idrow + 1,1) = imumax_s[idrow];
  }
  //generate rest of elements
  for (int idq = 1; idq < Gamma + 1; ++idq) {
    for (int idmu = sigma; idmu < mumax; ++idmu) {
      imsqb(idmu + 1,idq + 1) = -double(idmu + sigma + 1)*imsqb(idmu + 2,idq)/double(2*idmu + 1);
      if (idmu != sigma) {imsqb(idmu + 1,idq + 1) -= double(idmu - sigma)*imsqb(idmu,idq)/double(2*idmu + 1);}
    }
  }
}
double TrapezoidalIntegration(std::vector<double> & fx, std::vector<double> & dx) {
  //trapezoidal integration
  double trapzd = 0.0;
  size_t fxsz = fx.size();
  if (fxsz != dx.size()) {throw std::string("ERROR: IntegralPackage.hpp: TrapezoidalIntegration(): inconsistent vectors");}
  for (size_t idx = 0; idx < fxsz; ++idx) {
    trapzd += fx[idx]*dx[idx];
  }
  return trapzd;
}

#endif //_Integral_Package_