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

#ifndef _Lindh_Model_Hessian_
#define _Lindh_Model_Hessian_
#include "../atoms/AtomPackage.hpp"

//description:
//parameterization required for Lindh's model Hessian
//references://R. Lindh, A. Bernhardsson, G. Karlstroem, P.-A. Malmqvist, Chem. Phys. Lett., 241(4), 423, 1995       -> parameter set 0
//R. Lindh, personal communication - extracted from xTB                                                              -> parameter set 1

double LindhRAV(size_t atomA, size_t atomB, size_t parameterset = 1) {
  //function returning the rAv parameter for the Lindh approximate Hessian
  //parameters in Bohr
  size_t rowA = AtomRow(atomA);
  size_t rowB = AtomRow(atomB);
  size_t row_min = std::min(rowA,rowB);
  size_t row_max = std::max(rowA,rowB);
  double rav = 0.0;
  if (parameterset == 1) {
    if ((row_min == 1)&&(row_max == 1)) {rav = 1.35;}
    else if ((row_min == 1)&&(row_max == 2)) {rav = 2.10;}
    else if ((row_min == 1)&&(row_max == 3)) {rav = 2.53;}
    else if ((row_min == 2)&&(row_max == 2)) {rav = 2.87;}
    else if ((row_min == 2)&&(row_max == 3)) {rav = 3.40;}
    else if ((row_min == 3)&&(row_max == 3)) {rav = 3.40;}
  }
  else {
    if ((row_min == 1)&&(row_max == 1)) {rav = 1.35;}
    else if ((row_min == 1)&&(row_max == 2)) {rav = 2.10;}
    else if ((row_min == 1)&&(row_max == 3)) {rav = 2.53;}
    else if ((row_min == 2)&&(row_max == 2)) {rav = 2.87;}
    else if ((row_min == 2)&&(row_max == 3)) {rav = 3.80;}
    else if ((row_min == 3)&&(row_max == 3)) {rav = 4.50;}
  }
  return rav;
}
double LindhAAV(size_t atomA, size_t atomB, size_t parameterset = 1) {
  //function returning the aAv parameter for the Lindh approximate Hessian
  //parameters in rad
  size_t rowA = AtomRow(atomA);
  size_t rowB = AtomRow(atomB);
  size_t row_min = std::min(rowA,rowB);
  size_t row_max = std::max(rowA,rowB);
  double aav = 0.0;
  if (parameterset == 1) {
    if (row_min == 1) {
      if (row_max == 1) {aav = 1.0000;}
      else {aav = 0.3949;}
    }
    else {aav = 0.2800;}
  }
  else {
    if ((row_min == 1)&&(row_max == 1)) {aav = 1.0000;}
    else if ((row_min == 1)&&(row_max == 2)) {aav = 0.3949;}
    else if ((row_min == 1)&&(row_max == 3)) {aav = 0.3949;}
    else if ((row_min == 2)&&(row_max == 2)) {aav = 0.2800;}
    else if ((row_min == 2)&&(row_max == 3)) {aav = 0.1200;}
    else if ((row_min == 3)&&(row_max == 3)) {aav = 0.0600;}
  }
  return aav;
}
double LindhDAV(size_t atomA, size_t atomB, size_t parameterset = 1) {
  //function returning the dAv parameter for the Lindh approximate Hessian
  //parameters in rad
  double dav = 0.0;
  if (parameterset == 2) {
    size_t rowA = AtomRow(atomA);
    size_t rowB = AtomRow(atomB);
    size_t row_min = std::min(rowA,rowB);
    size_t row_max = std::max(rowA,rowB);
    if (row_min == 1) {
      if (row_max == 1) {dav = 0.0;}
      else {dav = 3.6;}
    }
    else {dav = 5.3;}
  }
  return dav;
}
double LindhRKr(size_t parameterset = 1) {
  //function returning the parameter rkr (for stretching)
  double rkr = 0.4;
  return rkr;
}
double LindhRKf(size_t parameterset = 1) {
  //function returning the parameter rkf (for bending)
  double rkf = 0.13;
  return rkf;
}
double LindhRKt(size_t parameterset = 1) {
  //function returning the parameter rkt (for torsions)
  double rkt = 0.0075;
  return rkt;
}

#endif //_Schlegel_Model_Hessian_
