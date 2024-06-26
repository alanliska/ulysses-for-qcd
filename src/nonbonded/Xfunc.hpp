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

#ifndef _XBond_Functions_
#define _XBond_Functions_
#include <vector>
#include <string>
#include <math.h>
#include "../parameters/H4par.hpp"

//description:
//set of functions that calculate the halogen-bond contribution to intermolecular forces
//applicable to semi-empirical methods
//references:
//X  -> J. Rezac, P. Hobza, Chem. Phys. Lett., 506, 286, 2011

//functions related to the X halogen bond correction
double EX(double rij, size_t donor, size_t acceptor, const std::string & method) {
  //halogen-bond correction
  return aparameter(donor,acceptor,method)*exp(-bparameter(donor,acceptor,method)*rij);
}
double gEX(double rij, size_t donor, size_t acceptor, const std::string & method) {
  //first-derivative of halogen-bond correction
  double bparam = bparameter(donor,acceptor,method);
  return -aparameter(donor,acceptor,method)*bparam*exp(-bparam*rij);
}
double hEX(double rij, size_t donor, size_t acceptor, const std::string & method) {
  //second-derivative of halogen-bond correction
  double bparam = bparameter(donor,acceptor,method);
  return aparameter(donor,acceptor,method)*bparam*bparam*exp(-bparam*rij);
}

#endif //_XBond_Functions_
