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

#ifndef _Burns_Basis_Set_
#define _Burns_Basis_Set_
#include <vector>
#include <string>
#include "STO.hpp"

//description:
//Burns exponents for Slater type orbitals

class BurnsBSet: public STOBSet {
public:
  BurnsBSet() {}
  ~BurnsBSet() {}
  std::vector<double> GetExponents(size_t atomnum) {
    //getter for the exponents of a given atom
    std::vector<double> exponents;
    if (atomnum == 1) {exponents.push_back(1.30);}      //H
    else if (atomnum == 3) {                            //Li
      exponents.push_back(0.60   );                        //2s
      exponents.push_back(0.50   );                        //2p
    }
    else if (atomnum == 4) {                            //Be
      exponents.push_back(0.90);                           //2s
      exponents.push_back(0.75);                           //2p
    }
    else if (atomnum == 5) {                            //B
      exponents.push_back(1.1225);                         //2s
      exponents.push_back(1.0000);                         //2p
    }
    else if (atomnum == 6) {                            //C
      exponents.push_back(1.500);                          //2s
      exponents.push_back(1.325);                          //2p
    }
    else if (atomnum == 7) {                            //N
      exponents.push_back(1.79825);                        //2s
      exponents.push_back(1.56875);                        //2p
    }
    else if (atomnum == 8) {                            //O
      exponents.push_back(2.200);                          //2s
      exponents.push_back(1.975);                          //2p
    }
    else if (atomnum == 9) {                            //F
      exponents.push_back(2.525);                          //2s
      exponents.push_back(2.300);                          //2p
    }
    else if (atomnum == 11) {                           //Na
      exponents.push_back(0.900);                          //3s
      exponents.push_back(0.533);                          //3p
      exponents.push_back(0.330);                          //3d
    }
    else if (atomnum == 12) {                           //Mg
      exponents.push_back(1.100);                          //3s
      exponents.push_back(0.700);                          //3p
      exponents.push_back(0.417);                          //3d
    }
    else if (atomnum == 13) {                           //Al
      exponents.push_back(1.317);                          //3s
      exponents.push_back(0.867);                          //3p
      exponents.push_back(0.500);                          //3d
    }
    else if (atomnum == 14) {                           //Si
      exponents.push_back(1.533);                          //3s
      exponents.push_back(1.083);                          //3p
      exponents.push_back(0.667);                          //3d
    }
    else if (atomnum == 15) {                           //P
      exponents.push_back(1.750);                          //3s
      exponents.push_back(1.300);                          //3p
      exponents.push_back(0.833);                          //3d
    }
    else if (atomnum == 16) {                           //S
      exponents.push_back(1.967);                          //3s
      exponents.push_back(1.517);                          //3p
      exponents.push_back(1.000);                          //3d
    }
    else if (atomnum == 17) {                           //Cl
      exponents.push_back(2.183);                          //3s
      exponents.push_back(1.733);                          //3p
      exponents.push_back(1.167);                          //3d
    }
    else if (atomnum == 19) {                           //K
      exponents.push_back(1.500);                          //4s
      exponents.push_back(0.700);                          //4p
      exponents.push_back(0.450);                          //4d
    }
    else if (atomnum == 20) {                           //Ca
      exponents.push_back(1.783);                          //4s
      exponents.push_back(0.825);                          //4p
      exponents.push_back(0.525);                          //4d
    }
    else if (atomnum == 21) {                           //Sc
      exponents.push_back(2.067);                          //3d
      exponents.push_back(0.950);                          //4s
      exponents.push_back(0.600);                          //4p
    }
    else if (atomnum == 22) {                           //Ti
      exponents.push_back(2.283);                          //3d
      exponents.push_back(1.075);                          //4s
      exponents.push_back(0.675);                          //4p
    }
    else if (atomnum == 23) {                           //V
      exponents.push_back(2.500);                          //3d
      exponents.push_back(1.200);                          //4s
      exponents.push_back(0.750);                          //4p
    }
    else if (atomnum == 24) {                           //Cr
      exponents.push_back(2.650);                          //3d
      exponents.push_back(1.325);                          //4s
      exponents.push_back(0.825);                          //4p
    }
    else if (atomnum == 25) {                           //Mn
      exponents.push_back(2.933);                          //3d
      exponents.push_back(1.450);                          //4s
      exponents.push_back(0.900);                          //4p
    }
    else if (atomnum == 26) {                           //Fe
      exponents.push_back(3.150);                          //3d
      exponents.push_back(1.575);                          //4s
      exponents.push_back(0.975);                          //4p
    }
    else if (atomnum == 27) {                           //Co
      exponents.push_back(3.367);                          //3d
      exponents.push_back(1.700);                          //4s
      exponents.push_back(1.050);                          //4p
    }
    else if (atomnum == 28) {                           //Ni
      exponents.push_back(3.583);                          //3d
      exponents.push_back(1.825);                          //4s
      exponents.push_back(1.125);                          //4p
    }
    else if (atomnum == 29) {                           //Cu
      exponents.push_back(3.733);                          //3d
      exponents.push_back(1.950);                          //4s
      exponents.push_back(1.200);                          //4p
    }
    else if (atomnum == 30) {                           //Zn
      exponents.push_back(4.067);                          //4s
      exponents.push_back(2.075);                          //4p
      exponents.push_back(1.272);                          //4d
    }
    else if (atomnum == 31) {                           //Ga
      exponents.push_back(2.238);                          //4s
      exponents.push_back(1.400);                          //4p
    }
    else if (atomnum == 32) {                           //Ge
      exponents.push_back(2.400);                          //4s
      exponents.push_back(1.563);                          //4p
    }
    else if (atomnum == 33) {                           //As
      exponents.push_back(2.563);                          //4s
      exponents.push_back(1.725);                          //4p
    }
    else if (atomnum == 34) {                           //Se
      exponents.push_back(2.726);                          //4s
      exponents.push_back(1.888);                          //4p
    }
    else if (atomnum == 35) {                           //Br
      exponents.push_back(2.889);                          //4s
      exponents.push_back(2.050);                          //4p
    }
    else {                                              //unavailable: rare gases and all atoms beyond Br, exclusively
      throw("ERROR: basissets/burns.hpp: BurnsBSet: GetExponents(): orbital exponents not defined for atom");
    }
    return exponents;
  }
  matrix<size_t> QuantumNumbers(size_t atomnum) {
    //calculation of quantum numbers for each atom
    std::vector<double> exps = GetExponents(atomnum);
    matrix<size_t> quantumnumb(exps.size(),2);
    if (atomnum == 1) quantumnumb(1,1) = 1;         //H
    else if ((atomnum > 2)&&(atomnum < 10)) {       //Li,Be,B,C,N,O,F
      quantumnumb(1,1) = 2;    //2s
      quantumnumb(2,1) = 2;    //2p
      quantumnumb(2,2) = 1;    //2p
    }
    else if ((atomnum > 10)&&(atomnum < 18)) {      //Na,Mg,Al,Si,P,S,Cl
      quantumnumb(1,1) = 3;    //3s
      quantumnumb(2,1) = 3;    //3p
      quantumnumb(2,2) = 1;    //3p
      quantumnumb(3,1) = 3;    //3d
      quantumnumb(3,2) = 2;    //3d
    }
    else if ((atomnum > 18)&&(atomnum < 31)) {      //K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn
      quantumnumb(1,1) = 3;    //3d
      quantumnumb(1,2) = 2;    //3d
      quantumnumb(2,1) = 4;    //4s
      quantumnumb(3,1) = 4;    //4p
      quantumnumb(3,2) = 1;    //4p
    }
    else if ((atomnum > 31)&&(atomnum < 36)) {      //Ga,Ge,As,Se,Br
      quantumnumb(1,1) = 4;    //4s
      quantumnumb(2,1) = 4;    //4p
      quantumnumb(2,2) = 1;    //4p
    }
    return quantumnumb;
  }
  size_t GetNAOs(size_t atomnum) {
    //function to calculate the number of AOs used for a given atom
    std::vector<double> exps = GetExponents(atomnum);
    return exps.size()*exps.size();
  }
};
#endif //_Burns_Basis_Set_
