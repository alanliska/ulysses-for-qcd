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

#ifndef _Lists_Package_
#define _Lists_Package_
#include <vector>
#include <fstream>
#include <iostream>

//description:
//functions that aid the access to lists

size_t LTriangMat2Array_nodiag(size_t p1, size_t p2, size_t Ntotal) {
  //function that returns the position in an array of the elements of lower triangular matrix; diagonal elements excluded
  size_t pA = std::min(p1,p2);
  size_t pB = std::max(p1,p2);
  size_t pos = pB - pA - 1;
  for (size_t idx = 1; idx < pA; ++idx) {
    pos += Ntotal - idx;
  }
  return pos;
}
int lin(int i1, int i2) {
  //address in packed array
  int imin = std::min(i1,i2);
  int imax = std::max(i1,i2);
  return imin + imax*(imax - 1)/2;
}

#endif //_Lists_Package_
