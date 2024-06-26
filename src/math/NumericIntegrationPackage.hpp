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

#ifndef _NumIntegral_Package_
#define _NumIntegral_Package_
#include <math.h>
#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>

//description:
//functions for numerical integration

void euler(double (* function)(double, double, std::vector<double>), std::vector<double> * yn, std::vector<double> * xn, std::vector<double> alpha, double y0, double imin, double imax, double h = 1e-8, int mode = 0) {
  //function that integrates a given one dimensional ODE, i.e., it solves the equation y' = f(x,y) using Euler's method.
  //*function is a pointer to the function to integrate; assume that x is first variable, then y
  //yn is the vector to store the values of y at different points x
  //xn is the vector to store the integration grid points x
  //alpha is a vector with parameters to pass to function
  //y0 is the initial value/boundary condition
  //imin is the lower limit of integration
  //imax is the upper limit of integration
  //h is the integration step, meaning that xn = nh; note that grid points are equidistant; if omitted, machine precision taken
  //mode determines type of output:
  //    0 -> write to variable (yn,xn) only
  //    1 -> write to file and variable
  //    2 -> write graphics and to variable
  //initial values
  double x = imin;
  double y = y0;
  xn->push_back(x);
  yn->push_back(y);
  //number steps
  int nstep = round((imax - imin)/h);
  for (size_t step = 1; step < nstep; ++step) {
    x += h;
    y += h*(*function)(x,y,alpha);
    xn->push_back(x);
    yn->push_back(y);
  }
  if (mode == 1) {
    //write to file
    std::ofstream outfile;
    outfile.open("euler.csv");
    outfile << "xn, yn" << std::endl;
    for (size_t idx = 0; idx < yn->size(); ++ idx) {
      outfile << xn->at(idx) << "," << yn->at(idx) << std::endl;
    }
  }
  else if (mode == 2) {
    //generate file with plot
  }
}

#endif //_NumIntegral_Package_
