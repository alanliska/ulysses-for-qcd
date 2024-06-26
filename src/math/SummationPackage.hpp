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

#ifndef _Summation_Package_
#define _Summation_Package_
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include "xsum/xsum.hpp"

//description:
//functions for exact summation
//NOTE THAT THIS REQUIRES C++ 11

//wrappers around xsum for easier manipulation
struct xsumsmall {
  //wrapper around xsum small accumulator
  //xsum is Radford Neal's exact summation algorithm, as described in https://arxiv.org/abs/1505.05571
  xsum::xsum_small_accumulator acc;
  xsumsmall() {}
  ~xsumsmall() {}
  //getter
  double getSum() {return xsum::xsum_round(&acc);}
  //addition actions
  void Add(double val) {xsum::xsum_add(&acc,val);}
  void Add(xsumsmall xsumobj) {xsum_add(&acc, &xsumobj.acc);}
  //reset
  void Reset() {xsum::xsum_init(&acc);}
};
struct xsumlarge {
  //wrapper around xsum large accumulator
  //xsum is Radford Neal's exact summation algorithm, as described in https://arxiv.org/abs/1505.05571
  xsum::xsum_large_accumulator acc;
  xsumlarge() {}
  ~xsumlarge() {}
  //getter
  double getSum() {return xsum::xsum_round(&acc);}
  //addition actions
  void Add(double val) {xsum::xsum_add(&acc,val);}
  void Add(xsumlarge xsumobj) {xsum_add(&acc, &xsumobj.acc);}
  void Add(xsumsmall xsumobj) {xsum_add(&acc, &xsumobj.acc);}
  //reset
  void Reset() {xsum::xsum_init(&acc);}
};
struct Kahan {
  //structure performing Kahan summation to reduce truncation errors, as described in
  //W. Kahan, Comm. ACM, 8(1), 40, 1965
  //if compilation is too optimal, this won't work
  double sum;                    //main accumulator
  double error;                  //estimation of error
  size_t nelements;
  double aux1;
  double aux2;
  Kahan() {Init();}
  ~Kahan() {}
  //getters
  double getSum() {return sum;}
  //reset and initialization
  void Init() {
    sum = 0.0;
    error = 0.0;
  }
  void Reset() {Init();}
  //add elements
  void Add(double value) {
    aux1 = value - error;             //removing from next point the estimate of error up until now
    aux2 = sum + aux1;                //adding next point
    error = (aux2 - sum) - aux1;      //estimating the new error
    sum = aux2;
  }
  void Add(std::vector<double> & data) {
    nelements = data.size();
    for (size_t idx = 0; idx < nelements; ++idx) {
      aux1 = data[idx] - error;             //removing from next point the estimate of error up until now
      aux2 = sum + aux1;                    //adding next point
      error = (aux2 - sum) - aux1;          //estimating the new error
      sum = aux2;
    }
  }
};
struct Neumaier {
  //structure performing Neumaier's variant of Kahan summation, as described in 
  //A. Neumaier, Z. Ang. Math. Mech. 54(1), 39, 1974
  //if compilation is too optimal, this won't work
  double sum;                    //main accumulator
  double error;                  //estimation of error
  size_t nelements;
  double aux1;
  double aux2;
  Neumaier() {Init();}
  ~Neumaier() {}
  //getters
  double getSum() {return sum + error;}
  //reset and initialization
  void Init() {
    sum = 0.0;
    error = 0.0;
  }
  void Reset() {Init();}
  //add elements
  void Add(double value) {
    aux2 = sum + value;
    if (fabs(sum) >= fabs(value)) {error += (sum - aux2) + value;}     //in this case we lose low-order digits of value
    else {error += (value - aux2) + sum;}                              //in this case we lose low-order digits of sum
    sum = aux2;
  }
  void Add(std::vector<double> & data) {
    nelements = data.size();
    for (size_t idx = 0; idx < nelements; ++idx) {
      aux2 = sum + data[idx];
      if (fabs(sum) >= fabs(data[idx])) {error += (sum - aux2) + data[idx];}   //in this case we lose low-order digits of value
      else {error += (data[idx] - aux2) + sum;}                                //in this case we lose low-order digits of sum
      sum = aux2;
    }
  }
};
struct Naive {
  //structure performing Naive summation
  double sum;                    //main accumulator
  size_t nelements;
  Naive() {Init();}
  ~Naive() {}
  //getters
  double getSum() {return sum;}
  //reset and initialization
  void Init() {sum = 0.0;}
  void Reset() {Init();}
  //add elements
  void Add(double value) {sum += value;}
  void Add(std::vector<double> & data) {
    nelements = data.size();
    for (size_t idx = 0; idx < nelements; ++idx) {
      sum += data[idx];
    }
  }
};

#endif //_Summation_Package_
