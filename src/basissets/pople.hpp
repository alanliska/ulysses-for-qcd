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

#ifndef _Pople_Basis_Set_
#define _Pople_Basis_Set_
#include <vector>
#include <string>
#include "STO.hpp"

//description:
//Slater exponents as modified by Pople for CNDO and INDO methods, to be used in Slater type orbitals

class PopleBSet: public STOBSet {
public:
  PopleBSet() {}
  ~PopleBSet() {}
  std::vector<double> GetExponents(size_t atomnum) {
    //getter for the exponents of a given atom
    std::vector<double> exponents;
    if (atomnum == 1) {exponents.push_back(1.20);}      //H
    else if (atomnum == 3) {                            //Li
      exponents.push_back(0.65);                            //2s
      exponents.push_back(0.65);                            //2p
    }
    else if (atomnum == 4) {                            //Be
      exponents.push_back(0.975);                           //2s
      exponents.push_back(0.975);                           //2p
    }
    else if (atomnum == 5) {                            //B
      exponents.push_back(1.30);                            //2s
      exponents.push_back(1.30);                            //2p
    }
    else if (atomnum == 6) {                            //C
      exponents.push_back(1.625);                           //2s
      exponents.push_back(1.625);                           //2p
    }
    else if (atomnum == 7) {                            //N
      exponents.push_back(1.95);                            //2s
      exponents.push_back(1.95);                            //2p
    }
    else if (atomnum == 8) {                            //O
      exponents.push_back(2.275);                           //2s
      exponents.push_back(2.275);                           //2p
    }
    else if (atomnum == 9) {                            //F
      exponents.push_back(2.60);                            //2s
      exponents.push_back(2.60);                            //2p
    }
    else if (atomnum == 11) {                           //Na
      exponents.push_back(0.733333);                        //3s
      exponents.push_back(0.733333);                        //3p
      exponents.push_back(0.306667);                        //3d
    }
    else if (atomnum == 12) {                           //Mg
      exponents.push_back(0.95);                            //3s
      exponents.push_back(0.95);                            //3p
      exponents.push_back(0.41);                            //3d
    }
    else if (atomnum == 13) {                           //Al
      exponents.push_back(1.1666667);                       //3s
      exponents.push_back(1.1666667);                       //3p
      exponents.push_back(0.52);                            //3d
    }
    else if (atomnum == 14) {                           //Si
      exponents.push_back(1.383333);                        //3s
      exponents.push_back(1.383333);                        //3p
      exponents.push_back(0.63);                            //3d
    }
    else if (atomnum == 15) {                           //P
      exponents.push_back(1.60);                            //3s
      exponents.push_back(1.60);                            //3p
      exponents.push_back(0.73);                            //3d
    }
    else if (atomnum == 16) {                           //S
      exponents.push_back(1.817);                           //3s
      exponents.push_back(1.817);                           //3p
      exponents.push_back(0.84);                            //3d
    }
    else if (atomnum == 17) {                           //Cl
      exponents.push_back(2.03333333);                      //3s
      exponents.push_back(2.03333333);                      //3p
      exponents.push_back(0.94);                            //3d
    }
    else {                                              //unavailable: He,Li, and all atoms beyond Cl
      throw("ERROR: basissets/pople.hpp: PopleBSet: GetExponents(): orbital exponents not defined for atom");
    }
    return exponents;
  }
  matrix<size_t> QuantumNumbers(size_t atomnum) {
    //calculation of quantum numbers for each atom
    std::vector<double> exps = GetExponents(atomnum);
    matrix<size_t> quantumnumb(exps.size(),2);
    if (atomnum > 86) {                             //unavailable: Fr,Ra,Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr,Rf,Db,Sg,Bh,Hs,Mt,Ds,Rg,Cn,Nh,Fl,Mc,Lv,Ts,Og
      throw("ERROR: basissets/pople.hpp: PopleBSet: QuantumNumbers(): atom not defined");
    }
    if (atomnum == 1) quantumnumb(1,1) = 1; //H
    else if ((atomnum > 2)&&(atomnum < 10)) {       //Li,Be,B,C,N,O,F,Ne
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
    return quantumnumb;
  }
  size_t GetNAOs(size_t atomnum) {
    //function to calculate the number of AOs used for a given atom
    std::vector<double> exps = GetExponents(atomnum);
    return exps.size()*exps.size();
  }
  int IPType() {return 1;}
};

#endif //_Pople_Basis_Set_
