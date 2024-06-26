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

#ifndef _Utilities_
#define _Utilities_
#include <iostream>
#include <sstream>

//description:
//general utilities

std::string szt2str(size_t number) {
  std::stringstream sstr;
  sstr << number;
  return sstr.str();
}
template <class T>
std::string to_string(T number) {
  std::stringstream sstr;
  sstr << number;
  return sstr.str();
}
template<typename ttype>
void swap(ttype & aval, ttype & bval) {
  //function that swaps two values
  ttype aux = aval;
  aval = bval;
  bval = aux;
}

#endif //_Utilities_
