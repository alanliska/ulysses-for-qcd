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

#ifndef _Search_Package_
#define _Search_Package_
#include <math.h>
#include <vector>
#include <string>
#include <fstream>

//description:
//searching functions

template <class ntype> 
void BisectionSearch(int & position, ntype value, std::vector<ntype> & data, int maxiter = 250) {
  //function performing bissection search on a vector
  size_t ndata = data.size();
  //obvious cases
  if (ndata == 0) {position = 0;}
  else if (value < data[0]) {position = 0;}
  else if (value == data[0]) {position = 1;}
  else if (value == data[ndata - 1]) {position = ndata - 1;}
  else if (value > data[ndata - 1]) {position = ndata;}
  else {
    int plower = 0;
    int pcurrent = ndata + 1;
    int pupper;
    for (size_t iter = 0; iter < maxiter; ++iter) {
      pupper = (plower + pcurrent)/2;
      if ((data[ndata - 1] >= data[0]) == (value >= data[pupper - 1])) {
        plower = pupper;
      }
      else {pcurrent = pupper;}
      if (pcurrent - plower <= 1) {break;}
    }
    position = plower;
  }
}

#endif //_Search_Package_
