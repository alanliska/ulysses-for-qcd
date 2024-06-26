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

#ifndef _H_plus_H_Bond_Correction_
#define _H_plus_H_Bond_Correction_
#include <vector>
#include <string>
#include <math.h>
#include "../UnitConversion.hpp"

//description:
//set of functions that calculate contribution of H-bonding to total energy according to the H+ method
//applicable to PM6
//references:
//H+ -> M. Korth, J. Chem. Theory Comput., 6, 3808, 2010
//   -> J. C. Kromann, A. S. Christensen, C. Steinmann, M. Korth, J. H. Jensen, PeerJ., 2, e449, 2014
double HplusCParam(int atmnr) {
  //the C parameters for H+
  double cpar = 0.0;
  if (atmnr == 7) {cpar = -0.11;}
  else if (atmnr == 8) {cpar = -0.12;}
  return cpar;
}

#endif //_H_plus_H_Bond_Correction_
