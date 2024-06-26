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

#ifndef _Schlegel_Model_Hessian_
#define _Schlegel_Model_Hessian_

//description:
//parameterization required for Schlegel's model Hessian
//references:
//H. B. Schlegel, Theoret. Chim. Acta (Berl.), 66, 333, 1984                      -> parameter set 0
//J. M. Wittbrodt, H. B. Schlegel, J. Mol. Struct. (Theochem), 398, 55, 1997      -> parameter set 1

double SchlegelB(size_t atm1, size_t atm2, size_t parameterset = 2) {
  //auxiliary function returning the B constants to be used in Schlegel's method for approximating the Hessian
  //note that I extended the Schlegel Ansatz to rare gases as well
  double schlegelb = 0.0;
  size_t atmA = std::min(atm1,atm2);
  size_t atmB = std::max(atm1,atm2);
  if (parameterset == 1) {                       //parameter set of H. B. Schlegel, Theoret. Chim. Acta (Berl.), 66, 333, 1984
    if (atmB < 3) {schlegelb = -0.244;}                          //1st-row:1st-row
    else if ((atmB > 2)&&(atmB < 11)) {
      if (atmA < 3) {schlegelb = 0.352;}                         //1st-row:2nd-row and 2nd-row:1st-row
      else {schlegelb = 1.085;}                                  //2nd-row:2nd-row
    }
    else if ((atmB > 10)&&(atmB < 19)) {
      if (atmA < 3) {schlegelb = 0.660;}                         //1st-row:3rd-row and 3rd-row:1st-row
      else if ((atmA > 2)&&(atmA < 11)) {schlegelb = 1.522;}     //2nd-row:3rd-row and 3rd-row:2nd-row
      else {schlegelb = 2.068;}                                  //3rd-row:3rd-row
    }
  }
  else if (parameterset == 2) {                  //parameter set of J. M. Wittbrodt, H. B. Schlegel, J. Mol. Struct. (Theochem), 398, 55, 1997
    if (atmB < 3) {schlegelb = -0.2573;}                         //1st-row:1st-row
    else if ((atmB > 2)&&(atmB < 11)) {
      if (atmA < 3) {schlegelb = 0.3401;}                        //1st-row:2nd-row and 2nd-row:1st-row
      else {schlegelb = 0.9652;}                                 //2nd-row:2nd-row
    }
    else if ((atmB > 10)&&(atmB < 19)) {
      if (atmA < 3) {schlegelb = 0.6937;}                        //1st-row:3rd-row and 3rd-row:1st-row
      else if ((atmA > 2)&&(atmA < 11)) {schlegelb = 1.2843;}    //2nd-row:3rd-row and 3rd-row:2nd-row
      else {schlegelb = 1.6925;}                                 //3rd-row:3rd-row
    }
    else if ((atmB > 18)&&(atmB < 37)) {
      if (atmA < 3) {schlegelb = 0.7126;}                        //1st-row:4th-row and 4th-row:1st-row
      else if ((atmA > 2)&&(atmA < 11)) {schlegelb = 1.4725;}    //2nd-row:4th-row and 4th-row:2nd-row
      else if ((atmA > 10)&&(atmA < 19)) {schlegelb = 1.8238;}   //3rd-row:4th-row and 4th-row:3rd-row
      else {schlegelb = 2.0203;}                                 //4th-row:4th-row
    }
    else if ((atmB > 36)&&(atmB < 55)) {
      if (atmA < 3) {schlegelb = 0.8335;}                        //1st-row:5th-row and 5th-row:1st-row
      else if ((atmA > 2)&&(atmA < 11)) {schlegelb = 1.6549;}    //2nd-row:5th-row and 5th-row:2nd-row
      else if ((atmA > 10)&&(atmA < 19)) {schlegelb = 2.1164;}   //3rd-row:5th-row and 5th-row:3rd-row
      else if ((atmA > 18)&&(atmA < 37)) {schlegelb = 2.2137;}   //4th-row:5th-row and 5th-row:4th-row
      else {schlegelb = 2.3718;}                                 //5th-row:5th-row
    }
    else if ((atmB > 54)&&(atmB < 87)) {
      if (atmA < 3) {schlegelb = 0.9491;}                        //1st-row:6th-row and 6th-row:1st-row
      else if ((atmA > 2)&&(atmA < 11)) {schlegelb = 1.7190;}    //2nd-row:6th-row and 6th-row:2nd-row
      else if ((atmA > 10)&&(atmA < 19)) {schlegelb = 2.3185;}   //3rd-row:6th-row and 6th-row:3rd-row
      else if ((atmA > 18)&&(atmA < 37)) {schlegelb = 2.5206;}   //4th-row:6th-row and 6th-row:4th-row
      else if ((atmA > 36)&&(atmA < 55)) {schlegelb = 2.5110;}   //5th-row:6th-row and 6th-row:5th-row
      else {schlegelb = 2.068;}                                  //6th-row:6th-row
    }
  }
  return schlegelb;
}

#endif //_Schlegel_Model_Hessian_
