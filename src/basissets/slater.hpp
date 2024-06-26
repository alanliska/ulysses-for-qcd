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

#ifndef _Slater_Basis_Set_
#define _Slater_Basis_Set_
#include <vector>
#include <string>
#include "STO.hpp"

//description:
//Slater exponents for Slater type orbitals

class SlaterBSet: public STOBSet {
public:
  SlaterBSet() {}
  SlaterBSet(const SlaterBSet & other) {}
  SlaterBSet & operator= (const SlaterBSet & other) {
    return *this;
  }
  ~SlaterBSet() {}
  std::vector<double> GetExponents(size_t atomnum) {
    //getter for the exponents of a given atom
    std::vector<double> exponents;
    if (atomnum == 1) {exponents.push_back(1.00);}      //H
    else if (atomnum == 2) {exponents.push_back(1.70);} //He
    else if (atomnum == 3) {                            //Li
      exponents.push_back(2.70);                            //1s
      exponents.push_back(0.65);                            //2s
      exponents.push_back(0.65);                            //2p
    }
    else if (atomnum == 4) {                            //Be
      exponents.push_back(3.70);                            //1s
      exponents.push_back(0.975);                           //2s
      exponents.push_back(0.975);                           //2p
    }
    else if (atomnum == 5) {                            //B
      exponents.push_back(4.70);                            //1s
      exponents.push_back(1.30);                            //2s
      exponents.push_back(1.30);                            //2p
    }
    else if (atomnum == 6) {                            //C
      exponents.push_back(5.70);                            //1s
      exponents.push_back(1.625);                           //2s
      exponents.push_back(1.625);                           //2p
    }
    else if (atomnum == 7) {                            //N
      exponents.push_back(6.70);                            //1s
      exponents.push_back(1.95);                            //2s
      exponents.push_back(1.95);                            //2p
    }
    else if (atomnum == 8) {                            //O
      exponents.push_back(7.70);                            //1s
      exponents.push_back(2.275);                           //2s
      exponents.push_back(2.275);                           //2p
    }
    else if (atomnum == 9) {                            //F
      exponents.push_back(8.70);                            //1s
      exponents.push_back(2.60);                            //2s
      exponents.push_back(2.60);                            //2p
    }
    else if (atomnum == 10) {                           //Ne
      exponents.push_back(9.70);                            //1s
      exponents.push_back(2.925);                           //2p
      exponents.push_back(2.925);                           //2p
    }
    else if (atomnum == 11) {                           //Na
      exponents.push_back(10.70);                           //1s
      exponents.push_back(3.425);                           //2s
      exponents.push_back(3.425);                           //2p
      exponents.push_back(0.733333);                        //3s
      exponents.push_back(0.733333);                        //3p
    }
    else if (atomnum == 12) {                           //Mg
      exponents.push_back(11.70);                           //1s
      exponents.push_back(3.925);                           //2s
      exponents.push_back(3.925);                           //2p
      exponents.push_back(0.95);                            //3s
      exponents.push_back(0.95);                            //3p
    }
    else if (atomnum == 13) {                           //Al
      exponents.push_back(12.70);                           //1s
      exponents.push_back(4.425);                           //2s
      exponents.push_back(4.425);                           //2p
      exponents.push_back(1.1666667);                       //3s
      exponents.push_back(1.1666667);                       //3p
    }
    else if (atomnum == 14) {                           //Si
      exponents.push_back(13.70);                           //1s
      exponents.push_back(4.925);                           //2s
      exponents.push_back(4.925);                           //2p
      exponents.push_back(1.383333);                        //3s
      exponents.push_back(1.383333);                        //3p
    }
    else if (atomnum == 15) {                           //P
      exponents.push_back(14.70);                           //1s
      exponents.push_back(5.425);                           //2s
      exponents.push_back(5.425);                           //2p
      exponents.push_back(1.60);                            //3s
      exponents.push_back(1.60);                            //3p
    }
    else if (atomnum == 16) {                           //S
      exponents.push_back(15.70);                           //1s
      exponents.push_back(5.925);                           //2s
      exponents.push_back(5.925);                           //2p
      exponents.push_back(1.817);                           //3s
      exponents.push_back(1.817);                           //3p
    }
    else if (atomnum == 17) {                           //Cl
      exponents.push_back(16.70);                           //1s
      exponents.push_back(6.425);                           //2s
      exponents.push_back(6.425);                           //2p
      exponents.push_back(2.03333333);                      //3s
      exponents.push_back(2.03333333);                      //3p
    }
    else if (atomnum == 18) {                           //Ar
      exponents.push_back(17.70);                           //1s
      exponents.push_back(6.925);                           //2s
      exponents.push_back(6.925);                           //2p
      exponents.push_back(2.25);                            //3s
      exponents.push_back(2.25);                            //3p
    }
    else if (atomnum == 19) {                           //K
      exponents.push_back(18.70);                           //1s
      exponents.push_back(7.425);                           //2s
      exponents.push_back(7.425);                           //2p
      exponents.push_back(2.583);                           //3s
      exponents.push_back(2.583);                           //3p
      exponents.push_back(0.00);                            //3d
      exponents.push_back(0.595);                           //4s
      exponents.push_back(0.595);                           //4p
    }
    else if (atomnum == 20) {                           //Ca
      exponents.push_back(19.70);                           //1s
      exponents.push_back(7.925);                           //2s
      exponents.push_back(7.925);                           //2p
      exponents.push_back(2.917);                           //3s
      exponents.push_back(2.917);                           //3p
      exponents.push_back(0.00);                            //3d
      exponents.push_back(0.770);                           //4s
      exponents.push_back(0.770);                           //4p
    }
    else if (atomnum == 21) {                           //Sc
      exponents.push_back(20.70);                           //1s
      exponents.push_back(8.425);                           //2s
      exponents.push_back(8.425);                           //2p
      exponents.push_back(3.25);                            //3s
      exponents.push_back(3.25);                            //3p
      exponents.push_back(1.00);                            //3d
      exponents.push_back(0.811);                           //4s
      exponents.push_back(0.811);                           //4p
    }
    else if (atomnum == 22) {                           //Ti
      exponents.push_back(21.70);                           //1s
      exponents.push_back(8.925);                           //2s
      exponents.push_back(8.925);                           //2p
      exponents.push_back(3.583);                           //3s
      exponents.push_back(3.583);                           //3p
      exponents.push_back(1.217);                           //3d
      exponents.push_back(0.851);                           //4s
      exponents.push_back(0.851);                           //4p
    }
    else if (atomnum == 23) {                           //V
      exponents.push_back(22.70);                           //1s
      exponents.push_back(9.425);                           //2s
      exponents.push_back(9.425);                           //2p
      exponents.push_back(3.917);                           //3s
      exponents.push_back(3.917);                           //3p
      exponents.push_back(1.433);                           //3d
      exponents.push_back(0.892);                           //4s
      exponents.push_back(0.892);                           //4p
    }
    else if (atomnum == 24) {                           //Cr
      exponents.push_back(23.70);                           //1s
      exponents.push_back(9.925);                           //2s
      exponents.push_back(9.925);                           //2p
      exponents.push_back(4.25);                            //3s
      exponents.push_back(4.25);                            //3p
      exponents.push_back(1.65);                            //3d
      exponents.push_back(0.932);                           //4s
      exponents.push_back(0.932);                           //4p
    }
    else if (atomnum == 25) {                           //Mn
      exponents.push_back(24.70);                           //1s
      exponents.push_back(10.425);                          //2s
      exponents.push_back(10.425);                          //2p
      exponents.push_back(4.583);                           //3s
      exponents.push_back(4.583);                           //3p
      exponents.push_back(1.867);                           //3d
      exponents.push_back(0.973);                           //4s
      exponents.push_back(0.973);                           //4p
    }
    else if (atomnum == 26) {                           //Fe
      exponents.push_back(25.70);                           //1s
      exponents.push_back(10.925);                          //2s
      exponents.push_back(10.925);                          //2p
      exponents.push_back(4.917);                           //3s
      exponents.push_back(4.917);                           //3p
      exponents.push_back(2.083);                           //3d
      exponents.push_back(1.014);                           //4s
      exponents.push_back(1.014);                           //4p
    }
    else if (atomnum == 27) {                           //Co
      exponents.push_back(26.70);                           //1s
      exponents.push_back(11.425);                          //2s
      exponents.push_back(11.425);                          //2p
      exponents.push_back(5.25);                            //3s
      exponents.push_back(5.25);                            //3p
      exponents.push_back(2.30);                            //3d
      exponents.push_back(1.054);                           //4s
      exponents.push_back(1.054);                           //4p
    }
    else if (atomnum == 28) {                           //Ni
      exponents.push_back(27.70);                           //1s
      exponents.push_back(11.925);                          //2s
      exponents.push_back(11.925);                          //2p
      exponents.push_back(5.583);                           //3s
      exponents.push_back(5.583);                           //3p
      exponents.push_back(2.517);                           //3d
      exponents.push_back(1.095);                           //4s
      exponents.push_back(1.095);                           //4p
    }
    else if (atomnum == 29) {                           //Cu
      exponents.push_back(28.70);                           //1s
      exponents.push_back(12.425);                          //2s
      exponents.push_back(12.425);                          //2p
      exponents.push_back(5.917);                           //3s
      exponents.push_back(5.917);                           //3p
      exponents.push_back(2.733);                           //3d
      exponents.push_back(1.135);                           //4s
      exponents.push_back(1.135);                           //4p
    }
    else if (atomnum == 30) {                           //Zn
      exponents.push_back(29.70);                           //1s
      exponents.push_back(12.925);                          //2s
      exponents.push_back(12.925);                          //2p
      exponents.push_back(6.25);                            //3s
      exponents.push_back(6.25);                            //3p
      exponents.push_back(2.95);                            //3d
      exponents.push_back(1.176);                           //4s
      exponents.push_back(1.176);                           //4p
    }
    else if (atomnum == 31) {                           //Ga
      exponents.push_back(30.70);                           //1s
      exponents.push_back(13.425);                          //2s
      exponents.push_back(13.425);                          //2p
      exponents.push_back(6.583);                           //3s
      exponents.push_back(6.583);                           //3p
      exponents.push_back(3.283);                           //3d
      exponents.push_back(1.351);                           //4s
      exponents.push_back(1.351);                           //4p
    }
    else if (atomnum == 32) {                           //Ge
      exponents.push_back(31.70);                           //1s
      exponents.push_back(13.925);                          //2s
      exponents.push_back(13.925);                          //2p
      exponents.push_back(6.917);                           //3s
      exponents.push_back(6.917);                           //3p
      exponents.push_back(3.617);                           //3d
      exponents.push_back(1.527);                           //4s
      exponents.push_back(1.527);                           //4p
    }
    else if (atomnum == 33) {                           //As
      exponents.push_back(32.70);                           //1s
      exponents.push_back(14.425);                          //2s
      exponents.push_back(14.425);                          //2p
      exponents.push_back(7.25);                            //3s
      exponents.push_back(7.25);                            //3p
      exponents.push_back(3.95);                            //3d
      exponents.push_back(1.703);                           //4s
      exponents.push_back(1.703);                           //4p
    }
    else if (atomnum == 34) {                           //Se
      exponents.push_back(33.70);                           //1s
      exponents.push_back(14.925);                          //2s
      exponents.push_back(14.925);                          //2p
      exponents.push_back(7.583);                           //3s
      exponents.push_back(7.583);                           //3p
      exponents.push_back(4.283);                           //3d
      exponents.push_back(1.878);                           //4s
      exponents.push_back(1.878);                           //4p
    }
    else if (atomnum == 35) {                           //Br
      exponents.push_back(34.70);                           //1s
      exponents.push_back(15.425);                          //2s
      exponents.push_back(15.425);                          //2p
      exponents.push_back(7.917);                           //3s
      exponents.push_back(7.917);                           //3p
      exponents.push_back(4.617);                           //3d
      exponents.push_back(2.054);                           //4s
      exponents.push_back(2.054);                           //4p
    }
    else if (atomnum == 36) {                           //Kr
      exponents.push_back(35.70);                           //1s
      exponents.push_back(15.925);                          //2s
      exponents.push_back(15.925);                          //2p
      exponents.push_back(8.25);                            //3s
      exponents.push_back(8.25);                            //3p
      exponents.push_back(4.95);                            //3d
      exponents.push_back(2.23);                            //4s
      exponents.push_back(2.23);                            //4p
    }
    else if (atomnum == 37) {                           //Rb
      exponents.push_back(36.70);                           //1s
      exponents.push_back(16.425);                          //2s
      exponents.push_back(16.425);                          //2p
      exponents.push_back(8.583);                           //3s
      exponents.push_back(8.583);                           //3p
      exponents.push_back(5.283);                           //3d
      exponents.push_back(2.50);                            //4s
      exponents.push_back(2.50);                            //4p
      exponents.push_back(0.00);                            //4d
      exponents.push_back(0.00);                            //4f
      exponents.push_back(0.55);                            //5s
      exponents.push_back(0.55);                            //5p
    }
    else if (atomnum == 38) {                           //Sr
      exponents.push_back(37.70);                           //1s
      exponents.push_back(16.925);                          //2s
      exponents.push_back(16.925);                          //2p
      exponents.push_back(8.917);                           //3s
      exponents.push_back(8.917);                           //3p
      exponents.push_back(5.617);                           //3d
      exponents.push_back(2.77);                            //4s
      exponents.push_back(2.77);                            //4p
      exponents.push_back(0.00);                            //4d
      exponents.push_back(0.00);                            //4f
      exponents.push_back(0.7125);                          //5s
      exponents.push_back(0.7125);                          //5p
    }
    else if (atomnum == 39) {                           //Y
      exponents.push_back(38.70);                           //1s
      exponents.push_back(17.425);                          //2s
      exponents.push_back(17.425);                          //2p
      exponents.push_back(9.25);                            //3s
      exponents.push_back(9.25);                            //3p
      exponents.push_back(5.95);                            //3d
      exponents.push_back(3.04);                            //4s
      exponents.push_back(3.04);                            //4p
      exponents.push_back(0.81);                            //4d
      exponents.push_back(0.00);                            //4f
      exponents.push_back(0.75);                            //5s
      exponents.push_back(0.75);                            //5p
    }
    else if (atomnum == 40) {                           //Zr
      exponents.push_back(39.70);                           //1s
      exponents.push_back(17.925);                          //2s
      exponents.push_back(17.925);                          //2p
      exponents.push_back(9.583);                           //3s
      exponents.push_back(9.583);                           //3p
      exponents.push_back(6.283);                           //3d
      exponents.push_back(3.31);                            //4s
      exponents.push_back(3.31);                            //4p
      exponents.push_back(0.986);                           //4d
      exponents.push_back(0.00);                            //4f
      exponents.push_back(0.7875);                          //5s
      exponents.push_back(0.7875);                          //5p
    }
    else if (atomnum == 41) {                           //Nb
      exponents.push_back(40.70);                           //1s
      exponents.push_back(18.425);                          //2s
      exponents.push_back(18.425);                          //2p
      exponents.push_back(9.917);                           //3s
      exponents.push_back(9.917);                           //3p
      exponents.push_back(6.617);                           //3d
      exponents.push_back(3.58);                            //4s
      exponents.push_back(3.58);                            //4p
      exponents.push_back(1.162);                           //4d
      exponents.push_back(0.00);                            //4f
      exponents.push_back(0.825);                           //5s
      exponents.push_back(0.825);                           //5p
    }
    else if (atomnum == 42) {                           //Mo
      exponents.push_back(41.70);                           //1s
      exponents.push_back(18.925);                          //2s
      exponents.push_back(18.925);                          //2p
      exponents.push_back(10.25);                           //3s
      exponents.push_back(10.25);                           //3p
      exponents.push_back(6.95);                            //3d
      exponents.push_back(3.851);                           //4s
      exponents.push_back(3.851);                           //4p
      exponents.push_back(1.338);                           //4d
      exponents.push_back(0.00);                            //4f
      exponents.push_back(0.8625);                          //5s
      exponents.push_back(0.8625);                          //5p
    }
    else if (atomnum == 43) {                           //Tc
      exponents.push_back(42.70);                           //1s
      exponents.push_back(19.425);                          //2s
      exponents.push_back(19.425);                          //2p
      exponents.push_back(10.583);                          //3s
      exponents.push_back(10.583);                          //3p
      exponents.push_back(7.283);                           //3d
      exponents.push_back(4.122);                           //4s
      exponents.push_back(4.122);                           //4p
      exponents.push_back(1.514);                           //4d
      exponents.push_back(0.00);                            //4f
      exponents.push_back(0.90);                            //5s
      exponents.push_back(0.90);                            //5p
    }
    else if (atomnum == 44) {                           //Ru
      exponents.push_back(43.70);                           //1s
      exponents.push_back(19.925);                          //2s
      exponents.push_back(19.925);                          //2p
      exponents.push_back(10.917);                          //3s
      exponents.push_back(10.917);                          //3p
      exponents.push_back(7.617);                           //3d
      exponents.push_back(4.392);                           //4s
      exponents.push_back(4.392);                           //4p
      exponents.push_back(1.689);                           //4d
      exponents.push_back(0.00);                            //4f
      exponents.push_back(0.9375);                          //5s
      exponents.push_back(0.9375);                          //5p
    }
    else if (atomnum == 45) {                           //Rh
      exponents.push_back(44.70);                           //1s
      exponents.push_back(20.425);                          //2s
      exponents.push_back(20.425);                          //2p
      exponents.push_back(11.25);                           //3s
      exponents.push_back(11.25);                           //3p
      exponents.push_back(7.95);                            //3d
      exponents.push_back(4.662);                           //4s
      exponents.push_back(4.662);                           //4p
      exponents.push_back(1.865);                           //4d
      exponents.push_back(0.00);                            //4f
      exponents.push_back(0.975);                           //5s
      exponents.push_back(0.975);                           //5p
    }
    else if (atomnum == 46) {                           //Pd
      exponents.push_back(45.70);                           //1s
      exponents.push_back(20.925);                          //2s
      exponents.push_back(20.925);                          //2p
      exponents.push_back(11.583);                          //3s
      exponents.push_back(11.583);                          //3p
      exponents.push_back(8.283);                           //3d
      exponents.push_back(4.932);                           //4s
      exponents.push_back(4.932);                           //4p
      exponents.push_back(2.041);                           //4d
      exponents.push_back(0.00);                            //4f
      exponents.push_back(1.0125);                          //5s
      exponents.push_back(1.0125);                          //5p
    }
    else if (atomnum == 47) {                           //Ag
      exponents.push_back(46.70);                           //1s
      exponents.push_back(21.425);                          //2s
      exponents.push_back(21.425);                          //2p
      exponents.push_back(11.917);                          //3s
      exponents.push_back(11.917);                          //3p
      exponents.push_back(8.617);                           //3d
      exponents.push_back(5.203);                           //4s
      exponents.push_back(5.203);                           //4p
      exponents.push_back(2.216);                           //4d
      exponents.push_back(0.00);                            //4f
      exponents.push_back(1.05);                            //5s
      exponents.push_back(1.05);                            //5p
    }
    else if (atomnum == 48) {                           //Cd
      exponents.push_back(47.70);                           //1s
      exponents.push_back(21.925);                          //2s
      exponents.push_back(21.925);                          //2p
      exponents.push_back(12.25);                           //3s
      exponents.push_back(12.25);                           //3p
      exponents.push_back(8.95);                            //3d
      exponents.push_back(5.473);                           //4s
      exponents.push_back(5.473);                           //4p
      exponents.push_back(2.392);                           //4d
      exponents.push_back(0.00);                            //4f
      exponents.push_back(1.0875);                          //5s
      exponents.push_back(1.0875);                          //5p
    }
    else if (atomnum == 49) {                           //In
      exponents.push_back(48.70);                           //1s
      exponents.push_back(22.425);                          //2s
      exponents.push_back(22.425);                          //2p
      exponents.push_back(12.583);                          //3s
      exponents.push_back(12.583);                          //3p
      exponents.push_back(9.283);                           //3d
      exponents.push_back(5.743);                           //4s
      exponents.push_back(5.743);                           //4p
      exponents.push_back(2.662);                           //4d
      exponents.push_back(0.00);                            //4f
      exponents.push_back(1.25);                            //5s
      exponents.push_back(1.25);                            //5p
    }
    else if (atomnum == 50) {                           //Sn
      exponents.push_back(49.70);                           //1s
      exponents.push_back(22.925);                          //2s
      exponents.push_back(22.925);                          //2p
      exponents.push_back(12.917);                          //3s
      exponents.push_back(12.917);                          //3p
      exponents.push_back(9.617);                           //3d
      exponents.push_back(6.0135);                          //4s
      exponents.push_back(6.0135);                          //4p
      exponents.push_back(2.932);                           //4d
      exponents.push_back(0.00);                            //4f
      exponents.push_back(1.4125);                          //5s
      exponents.push_back(1.4125);                          //5p
    }
    else if (atomnum == 51) {                           //Sb
      exponents.push_back(50.70);                           //1s
      exponents.push_back(23.425);                          //2s
      exponents.push_back(23.425);                          //2p
      exponents.push_back(13.25);                           //3s
      exponents.push_back(13.25);                           //3p
      exponents.push_back(9.95);                            //3d
      exponents.push_back(6.284);                           //4s
      exponents.push_back(6.284);                           //4p
      exponents.push_back(3.203);                           //4d
      exponents.push_back(0.00);                            //4f
      exponents.push_back(1.575);                           //5s
      exponents.push_back(1.575);                           //5p
    }
    else if (atomnum == 52) {                           //Te
      exponents.push_back(51.70);                           //1s
      exponents.push_back(23.925);                          //2s
      exponents.push_back(23.925);                          //2p
      exponents.push_back(13.583);                          //3s
      exponents.push_back(13.583);                          //3p
      exponents.push_back(10.283);                          //3d
      exponents.push_back(6.554);                           //4s
      exponents.push_back(6.554);                           //4p
      exponents.push_back(3.473);                           //4d
      exponents.push_back(0.00);                            //4f
      exponents.push_back(1.7375);                          //5s
      exponents.push_back(1.7375);                          //5p
    }
    else if (atomnum == 53) {                           //I
      exponents.push_back(52.70);                           //1s
      exponents.push_back(24.425);                          //2s
      exponents.push_back(24.425);                          //2p
      exponents.push_back(13.917);                          //3s
      exponents.push_back(13.917);                          //3p
      exponents.push_back(10.617);                          //3d
      exponents.push_back(6.824);                           //4s
      exponents.push_back(6.824);                           //4p
      exponents.push_back(3.743);                           //4d
      exponents.push_back(0.00);                            //4f
      exponents.push_back(1.90);                            //5s
      exponents.push_back(1.90);                            //5p
    }
    else if (atomnum == 54) {                           //Xe
      exponents.push_back(53.70);                           //1s
      exponents.push_back(24.925);                          //2s
      exponents.push_back(24.925);                          //2p
      exponents.push_back(14.25) ;                          //3s
      exponents.push_back(14.25) ;                          //3p
      exponents.push_back(10.95);                           //3d
      exponents.push_back(7.095);                           //4s
      exponents.push_back(7.095);                           //4p
      exponents.push_back(4.0135);                          //4d
      exponents.push_back(0.00);                            //4f
      exponents.push_back(2.0625);                          //5s
      exponents.push_back(2.0625);                          //5p
    }
    else if (atomnum == 55) {                           //Cs
      exponents.push_back(54.70);                           //1s
      exponents.push_back(25.425);                          //2s
      exponents.push_back(25.425);                          //2p
      exponents.push_back(14.583);                          //3s
      exponents.push_back(14.583);                          //3p
      exponents.push_back(11.283);                          //3d
      exponents.push_back(7.365);                           //4s
      exponents.push_back(7.365);                           //4p
      exponents.push_back(4.284);                           //4d
      exponents.push_back(0.00);                            //4f
      exponents.push_back(2.3125);                          //5s
      exponents.push_back(2.3125);                          //5p
      exponents.push_back(0.00);                            //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.524);                           //6s
      exponents.push_back(0.524);                           //6p
    }
    else if (atomnum == 56) {                           //Ba
      exponents.push_back(55.70);                           //1s
      exponents.push_back(25.925);                          //2s
      exponents.push_back(25.925);                          //2p
      exponents.push_back(14.917);                          //3s
      exponents.push_back(14.917);                          //3p
      exponents.push_back(11.617);                          //3d
      exponents.push_back(7.635);                           //4s
      exponents.push_back(7.635);                           //4p
      exponents.push_back(4.554);                           //4d
      exponents.push_back(0.00);                            //4f
      exponents.push_back(2.5625);                          //5s
      exponents.push_back(2.5625);                          //5p
      exponents.push_back(0.00);                            //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.679);                           //6s
      exponents.push_back(0.679);                           //6p
    }
    else if (atomnum == 57) {                           //La
      exponents.push_back(56.70);                           //1s
      exponents.push_back(26.425);                          //2s
      exponents.push_back(26.425);                          //2p
      exponents.push_back(15.25);                           //3s
      exponents.push_back(15.25);                           //3p
      exponents.push_back(11.95);                           //3d
      exponents.push_back(7.905);                           //4s
      exponents.push_back(7.905);                           //4p
      exponents.push_back(4.824);                           //4d
      exponents.push_back(15.405);                          //4f
      exponents.push_back(2.8125);                          //5s
      exponents.push_back(2.8125);                          //5p
      exponents.push_back(0.75);                            //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.714);                           //6s
      exponents.push_back(0.714);                           //6p
    }
    else if (atomnum == 58) {                           //Ce
      exponents.push_back(57.70);                           //1s
      exponents.push_back(26.925);                          //2s
      exponents.push_back(26.925);                          //2p
      exponents.push_back(15.583);                          //3s
      exponents.push_back(15.583);                          //3p
      exponents.push_back(12.283);                          //3d
      exponents.push_back(8.176);                           //4s
      exponents.push_back(8.176);                           //4p
      exponents.push_back(5.095);                           //4d
      exponents.push_back(3.243);                           //4f
      exponents.push_back(2.85);                            //5s
      exponents.push_back(2.85);                            //5p
      exponents.push_back(0.75);                            //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.714);                           //6s
      exponents.push_back(0.714);                           //6p
    }
    else if (atomnum == 59) {                           //Pr
      exponents.push_back(58.70);                           //1s
      exponents.push_back(27.425);                          //2s
      exponents.push_back(27.425);                          //2p
      exponents.push_back(15.917);                          //3s
      exponents.push_back(15.917);                          //3p
      exponents.push_back(12.617);                          //3d
      exponents.push_back(8.446);                           //4s
      exponents.push_back(8.446);                           //4p
      exponents.push_back(5.365);                           //4d
      exponents.push_back(3.324);                           //4f
      exponents.push_back(2.675);                           //5s
      exponents.push_back(2.675);                           //5p
      exponents.push_back(0.5875);                          //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.6785);                          //6s
      exponents.push_back(0.6785);                          //6p
    }
    else if (atomnum == 60) {                           //Nd
      exponents.push_back(59.70);                           //1s
      exponents.push_back(27.925);                          //2s
      exponents.push_back(27.925);                          //2p
      exponents.push_back(16.25);                           //3s
      exponents.push_back(16.25);                           //3p
      exponents.push_back(12.95);                           //3d
      exponents.push_back(8.716);                           //4s
      exponents.push_back(8.716);                           //4p
      exponents.push_back(5.635);                           //4d
      exponents.push_back(3.00);                            //4f
      exponents.push_back(2.7125);                          //5s
      exponents.push_back(2.7125);                          //5p
      exponents.push_back(0.5875);                          //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.679);                           //6s
      exponents.push_back(0.679);                           //6p
    }
    else if (atomnum == 61) {                           //Pm
      exponents.push_back(60.70);                           //1s
      exponents.push_back(28.425);                          //2s
      exponents.push_back(28.425);                          //2p
      exponents.push_back(16.583);                          //3s
      exponents.push_back(16.583);                          //3p
      exponents.push_back(13.283);                          //3d
      exponents.push_back(8.986);                           //4s
      exponents.push_back(8.986);                           //4p
      exponents.push_back(5.905);                           //4d
      exponents.push_back(3.6757);                          //4f
      exponents.push_back(2.75);                            //5s
      exponents.push_back(2.75);                            //5p
      exponents.push_back(0.5875);                          //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.679);                           //6s
      exponents.push_back(0.679);                           //6p
    }
    else if (atomnum == 62) {                           //Sm
      exponents.push_back(61.70);                           //1s
      exponents.push_back(28.925);                          //2s
      exponents.push_back(28.925);                          //2p
      exponents.push_back(16.917);                          //3s
      exponents.push_back(16.917);                          //3p
      exponents.push_back(13.617);                          //3d
      exponents.push_back(9.257);                           //4s
      exponents.push_back(9.257);                           //4p
      exponents.push_back(6.176);                           //4d
      exponents.push_back(3.851);                           //4f
      exponents.push_back(2.7875);                          //5s
      exponents.push_back(2.7875);                          //5p
      exponents.push_back(0.5875);                          //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.679);                           //6s
      exponents.push_back(0.679);                           //6p
    }
    else if (atomnum == 63) {                           //Eu
      exponents.push_back(62.70);                           //1s
      exponents.push_back(29.425);                          //2s
      exponents.push_back(29.425);                          //2p
      exponents.push_back(17.25);                           //3s
      exponents.push_back(17.25);                           //3p
      exponents.push_back(13.95);                           //3d
      exponents.push_back(9.527);                           //4s
      exponents.push_back(9.527);                           //4p
      exponents.push_back(6.446);                           //4d
      exponents.push_back(4.027);                           //4f
      exponents.push_back(2.825);                           //5s
      exponents.push_back(2.825);                           //5p
      exponents.push_back(0.5875);                          //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.679);                           //6s
      exponents.push_back(0.679);                           //6p
    }
    else if (atomnum == 64) {                           //Gd
      exponents.push_back(63.70);                           //1s
      exponents.push_back(29.925);                          //2s
      exponents.push_back(29.925);                          //2p
      exponents.push_back(17.583);                          //3s
      exponents.push_back(17.583);                          //3p
      exponents.push_back(14.283);                          //3d
      exponents.push_back(9.797);                           //4s
      exponents.push_back(9.797);                           //4p
      exponents.push_back(6.716);                           //4d
      exponents.push_back(4.297);                           //4f
      exponents.push_back(3.075);                           //5s
      exponents.push_back(3.075);                           //5p
      exponents.push_back(0.75);                            //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.714);                           //6s
      exponents.push_back(0.714);                           //6p
    }
    else if (atomnum == 65) {                           //Tb
      exponents.push_back(64.70);                           //1s
      exponents.push_back(30.425);                          //2s
      exponents.push_back(30.425);                          //2p
      exponents.push_back(17.917);                          //3s
      exponents.push_back(17.917);                          //3p
      exponents.push_back(14.617);                          //3d
      exponents.push_back(10.0676);                         //4s
      exponents.push_back(10.0676);                         //4p
      exponents.push_back(6.986);                           //4d
      exponents.push_back(4.473);                           //4f
      exponents.push_back(3.1125);                          //5s
      exponents.push_back(3.1125);                          //5p
      exponents.push_back(0.75);                            //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.714);                           //6s
      exponents.push_back(0.714);                           //6p
    }
    else if (atomnum == 66) {                           //Dy
      exponents.push_back(65.70);                           //1s
      exponents.push_back(30.925);                          //2s
      exponents.push_back(30.925);                          //2p
      exponents.push_back(18.25);                           //3s
      exponents.push_back(18.25);                           //3p
      exponents.push_back(14.95);                           //3d
      exponents.push_back(10.338);                          //4s
      exponents.push_back(10.338);                          //4p
      exponents.push_back(7.257);                           //4d
      exponents.push_back(4.554);                           //4f
      exponents.push_back(2.9375);                          //5s
      exponents.push_back(2.9375);                          //5p
      exponents.push_back(0.5875);                          //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.679);                           //6s
      exponents.push_back(0.679);                           //6p
    }
    else if (atomnum == 67) {                           //Ho
      exponents.push_back(66.70);                           //1s
      exponents.push_back(31.425);                          //2s
      exponents.push_back(31.425);                          //2p
      exponents.push_back(18.583);                          //3s
      exponents.push_back(18.583);                          //3p
      exponents.push_back(15.283);                          //3d
      exponents.push_back(10.608);                          //4s
      exponents.push_back(10.608);                          //4p
      exponents.push_back(7.527);                           //4d
      exponents.push_back(4.730);                           //4f
      exponents.push_back(2.975);                           //5s
      exponents.push_back(2.975);                           //5p
      exponents.push_back(0.5875);                          //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.679);                           //6s
      exponents.push_back(0.679);                           //6p
    }
    else if (atomnum == 68) {                           //Er
      exponents.push_back(67.70);                           //1s
      exponents.push_back(31.925);                          //2s
      exponents.push_back(31.925);                          //2p
      exponents.push_back(18.917);                          //3s
      exponents.push_back(18.917);                          //3p
      exponents.push_back(15.617);                          //3d
      exponents.push_back(10.878);                          //4s
      exponents.push_back(10.878);                          //4p
      exponents.push_back(7.797);                           //4d
      exponents.push_back(4.905);                           //4f
      exponents.push_back(3.0125);                          //5s
      exponents.push_back(3.0125);                          //5p
      exponents.push_back(0.5875);                          //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.679);                           //6s
      exponents.push_back(0.679);                           //6p
    }
    else if (atomnum == 69) {                           //Tm
      exponents.push_back(68.70);                           //1s
      exponents.push_back(32.425);                          //2s
      exponents.push_back(32.425);                          //2p
      exponents.push_back(19.25);                           //3s
      exponents.push_back(19.25);                           //3p
      exponents.push_back(15.95);                           //3d
      exponents.push_back(11.149);                          //4s
      exponents.push_back(11.149);                          //4p
      exponents.push_back(8.068);                           //4d
      exponents.push_back(5.081);                           //4f
      exponents.push_back(3.05);                            //5s
      exponents.push_back(3.05);                            //5p
      exponents.push_back(0.5875);                          //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.679);                           //6s
      exponents.push_back(0.679);                           //6p
    }
    else if (atomnum == 70) {                           //Yb
      exponents.push_back(69.70);                           //1s
      exponents.push_back(32.925);                          //2s
      exponents.push_back(32.925);                          //2p
      exponents.push_back(19.583);                          //3s
      exponents.push_back(19.583);                          //3p
      exponents.push_back(16.283);                          //3d
      exponents.push_back(11.42);                           //4s
      exponents.push_back(11.42);                           //4p
      exponents.push_back(8.338);                           //4d
      exponents.push_back(5.257);                           //4f
      exponents.push_back(3.0875);                          //5s
      exponents.push_back(3.0875);                          //5p
      exponents.push_back(0.5875);                          //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.679);                           //6s
      exponents.push_back(0.679);                           //6p
    }
    else if (atomnum == 71) {                           //Lu
      exponents.push_back(70.70);                           //1s
      exponents.push_back(33.425);                          //2s
      exponents.push_back(33.425);                          //2p
      exponents.push_back(19.917);                          //3s
      exponents.push_back(19.917);                          //3p
      exponents.push_back(16.617);                          //3d
      exponents.push_back(11.689);                          //4s
      exponents.push_back(11.689);                          //4p
      exponents.push_back(8.608);                           //4d
      exponents.push_back(5.527);                           //4f
      exponents.push_back(3.3375);                          //5s
      exponents.push_back(3.3375);                          //5p
      exponents.push_back(0.75);                            //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.714);                           //6s
      exponents.push_back(0.714);                           //6p
    }
    else if (atomnum == 72) {                           //Hf
      exponents.push_back(71.70);                           //1s
      exponents.push_back(33.925);                          //2s
      exponents.push_back(33.925);                          //2p
      exponents.push_back(20.25);                           //3s
      exponents.push_back(20.25);                           //3p
      exponents.push_back(16.95);                           //3d
      exponents.push_back(11.959);                          //4s
      exponents.push_back(11.959);                          //4p
      exponents.push_back(8.878);                           //4d
      exponents.push_back(5.797);                           //4f
      exponents.push_back(3.5875);                          //5s
      exponents.push_back(3.5875);                          //5p
      exponents.push_back(0.9125);                          //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.75);                            //6s
      exponents.push_back(0.75);                            //6p
    }
    else if (atomnum == 73) {                           //Ta
      exponents.push_back(72.70);                           //1s
      exponents.push_back(34.425);                          //2s
      exponents.push_back(34.425);                          //2p
      exponents.push_back(20.583);                          //3s
      exponents.push_back(20.583);                          //3p
      exponents.push_back(17.283);                          //3d
      exponents.push_back(12.23);                           //4s
      exponents.push_back(12.23);                           //4p
      exponents.push_back(9.149);                           //4d
      exponents.push_back(6.068);                           //4f
      exponents.push_back(3.8375);                          //5s
      exponents.push_back(3.8375);                          //5p
      exponents.push_back(1.075);                           //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.786);                           //6s
      exponents.push_back(0.786);                           //6p
    }
    else if (atomnum == 74) {                           //W
      exponents.push_back(73.70);                           //1s
      exponents.push_back(34.925);                          //2s
      exponents.push_back(34.925);                          //2p
      exponents.push_back(20.917);                          //3s
      exponents.push_back(20.917);                          //3p
      exponents.push_back(17.617);                          //3d
      exponents.push_back(12.50);                           //4s
      exponents.push_back(12.50);                           //4p
      exponents.push_back(9.419);                           //4d
      exponents.push_back(6.338);                           //4f
      exponents.push_back(4.0875);                          //5s
      exponents.push_back(4.0875);                          //5p
      exponents.push_back(1.2375);                          //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.821);                           //6s
      exponents.push_back(0.821);                           //6p
    }
    else if (atomnum == 75) {                           //Re
      exponents.push_back(74.70);                           //1s
      exponents.push_back(35.425);                          //2s
      exponents.push_back(35.425);                          //2p
      exponents.push_back(21.25);                           //3s
      exponents.push_back(21.25);                           //3p
      exponents.push_back(17.95);                           //3d
      exponents.push_back(12.77);                           //4s
      exponents.push_back(12.77);                           //4p
      exponents.push_back(9.689);                           //4d
      exponents.push_back(6.6081);                          //4f
      exponents.push_back(4.3375);                          //5s
      exponents.push_back(4.3375);                          //5p
      exponents.push_back(1.40);                            //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.857);                           //6s
      exponents.push_back(0.857);                           //6p
    }
    else if (atomnum == 76) {                           //Os
      exponents.push_back(75.70);                           //1s
      exponents.push_back(35.925);                          //2s
      exponents.push_back(35.925);                          //2p
      exponents.push_back(21.583);                          //3s
      exponents.push_back(21.583);                          //3p
      exponents.push_back(18.283);                          //3d
      exponents.push_back(13.0405);                         //4s
      exponents.push_back(13.0405);                         //4p
      exponents.push_back(9.959);                           //4d
      exponents.push_back(6.878);                           //4f
      exponents.push_back(4.5875);                          //5s
      exponents.push_back(4.5875);                          //5p
      exponents.push_back(1.5625);                          //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.893);                           //6s
      exponents.push_back(0.893);                           //6p
    }
    else if (atomnum == 77) {                           //Ir
      exponents.push_back(76.70);                           //1s
      exponents.push_back(36.425);                          //2s
      exponents.push_back(36.425);                          //2p
      exponents.push_back(21.917);                          //3s
      exponents.push_back(21.917);                          //3p
      exponents.push_back(18.617);                          //3d
      exponents.push_back(13.311);                          //4s
      exponents.push_back(13.311);                          //4p
      exponents.push_back(10.23);                           //4d
      exponents.push_back(7.149);                           //4f
      exponents.push_back(4.8375);                          //5s
      exponents.push_back(4.8375);                          //5p
      exponents.push_back(1.725);                           //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.9286);                          //6s
      exponents.push_back(0.9286);                          //6p
    }
    else if (atomnum == 78) {                           //Pt
      exponents.push_back(77.70);                           //1s
      exponents.push_back(36.925);                          //2s
      exponents.push_back(36.925);                          //2p
      exponents.push_back(22.25);                           //3s
      exponents.push_back(22.25);                           //3p
      exponents.push_back(18.95);                           //3d
      exponents.push_back(13.581);                          //4s
      exponents.push_back(13.581);                          //4p
      exponents.push_back(10.50);                           //4d
      exponents.push_back(7.419);                           //4f
      exponents.push_back(5.0875);                          //5s
      exponents.push_back(5.0875);                          //5p
      exponents.push_back(1.8875);                          //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(0.964);                           //6s
      exponents.push_back(0.964);                           //6p
    }
    else if (atomnum == 79) {                           //Au
      exponents.push_back(78.70);                           //1s
      exponents.push_back(37.425);                          //2s
      exponents.push_back(37.425);                          //2p
      exponents.push_back(22.583);                          //3s
      exponents.push_back(22.583);                          //3p
      exponents.push_back(19.283);                          //3d
      exponents.push_back(13.851);                          //4s
      exponents.push_back(13.851);                          //4p
      exponents.push_back(10.77);                           //4d
      exponents.push_back(7.689);                           //4f
      exponents.push_back(5.3375);                          //5s
      exponents.push_back(5.3375);                          //5p
      exponents.push_back(2.05);                            //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(1.00);                            //6s
      exponents.push_back(1.00);                            //6p
    }
    else if (atomnum == 80) {                           //Hg
      exponents.push_back(79.70);                           //1s
      exponents.push_back(37.925);                          //2s
      exponents.push_back(37.925);                          //2p
      exponents.push_back(22.917);                          //3s
      exponents.push_back(22.917);                          //3p
      exponents.push_back(19.617);                          //3d
      exponents.push_back(14.1216);                         //4s
      exponents.push_back(14.1216);                         //4p
      exponents.push_back(11.04);                           //4d
      exponents.push_back(7.959);                           //4f
      exponents.push_back(5.5875);                          //5s
      exponents.push_back(5.5875);                          //5p
      exponents.push_back(2.2125);                          //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(1.036);                           //6s
      exponents.push_back(1.036);                           //6p
    }
    else if (atomnum == 81) {                           //Tl
      exponents.push_back(80.70);                           //1s
      exponents.push_back(38.425);                          //2s
      exponents.push_back(38.425);                          //2p
      exponents.push_back(23.25);                           //3s
      exponents.push_back(23.25);                           //3p
      exponents.push_back(19.95);                           //3d
      exponents.push_back(14.392);                          //4s
      exponents.push_back(14.392);                          //4p
      exponents.push_back(11.31);                           //4d
      exponents.push_back(8.230);                           //4f
      exponents.push_back(5.8375);                          //5s
      exponents.push_back(5.8375);                          //5p
      exponents.push_back(2.4625);                          //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(1.19);                            //6s
      exponents.push_back(1.19);                            //6p
    }
    else if (atomnum == 82) {                           //Pb
      exponents.push_back(81.70);                           //1s
      exponents.push_back(38.925);                          //2s
      exponents.push_back(38.925);                          //2p
      exponents.push_back(23.583);                          //3s
      exponents.push_back(23.583);                          //3p
      exponents.push_back(20.283);                          //3d
      exponents.push_back(14.662);                          //4s
      exponents.push_back(14.662);                          //4p
      exponents.push_back(11.581);                          //4d
      exponents.push_back(8.50);                            //4f
      exponents.push_back(6.0875);                          //5s
      exponents.push_back(6.0875);                          //5p
      exponents.push_back(2.7125);                          //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(1.345);                           //6s
      exponents.push_back(1.345);                           //6p
    }
    else if (atomnum == 83) {                           //Bi
      exponents.push_back(82.70);                           //1s
      exponents.push_back(39.425);                          //2s
      exponents.push_back(39.425);                          //2p
      exponents.push_back(23.917);                          //3s
      exponents.push_back(23.917);                          //3p
      exponents.push_back(20.617);                          //3d
      exponents.push_back(14.932);                          //4s
      exponents.push_back(14.932);                          //4p
      exponents.push_back(11.851);                          //4d
      exponents.push_back(8.77);                            //4f
      exponents.push_back(6.3375);                          //5s
      exponents.push_back(6.3375);                          //5p
      exponents.push_back(2.9625);                          //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(1.50);                            //6s
      exponents.push_back(1.50);                            //6p
    }
    else if (atomnum == 84) {                           //Po
      exponents.push_back(83.70);                           //1s
      exponents.push_back(39.925);                          //2s
      exponents.push_back(39.925);                          //2p
      exponents.push_back(24.25);                           //3s
      exponents.push_back(24.25);                           //3p
      exponents.push_back(20.95);                           //3d
      exponents.push_back(15.203);                          //4s
      exponents.push_back(15.203);                          //4p
      exponents.push_back(12.122);                          //4d
      exponents.push_back(9.0405);                          //4f
      exponents.push_back(6.5875);                          //5s
      exponents.push_back(6.5875);                          //5p
      exponents.push_back(3.2125);                          //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(1.655);                           //6s
      exponents.push_back(1.655);                           //6p
    }
    else if (atomnum == 85) {                           //At
      exponents.push_back(84.70);                           //1s
      exponents.push_back(40.425);                          //2s
      exponents.push_back(40.425);                          //2p
      exponents.push_back(24.583);                          //3s
      exponents.push_back(24.583);                          //3p
      exponents.push_back(21.283);                          //3d
      exponents.push_back(15.473);                          //4s
      exponents.push_back(15.473);                          //4p
      exponents.push_back(12.392);                          //4d
      exponents.push_back(9.31);                            //4f
      exponents.push_back(6.8375);                          //5s
      exponents.push_back(6.8375);                          //5p
      exponents.push_back(3.4625);                          //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(1.81);                            //6s
      exponents.push_back(1.81);                            //6p
    }
    else if (atomnum == 86) {                           //Rn
      exponents.push_back(85.70);                           //1s
      exponents.push_back(40.925);                          //2s
      exponents.push_back(40.925);                          //2p
      exponents.push_back(24.917);                          //3s
      exponents.push_back(24.917);                          //3p
      exponents.push_back(21.617);                          //3d
      exponents.push_back(15.743);                          //4s
      exponents.push_back(15.743);                          //4p
      exponents.push_back(12.662);                          //4d
      exponents.push_back(9.581);                           //4f
      exponents.push_back(7.0875);                          //5s
      exponents.push_back(7.0875);                          //5p
      exponents.push_back(3.7125);                          //5d
      exponents.push_back(0.00);                            //5f
      exponents.push_back(1.964);                           //6s
      exponents.push_back(1.964);                           //6p
    }
    else if (atomnum > 86) {                            //unavailable: Fr,Ra,Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr,Rf,Db,Sg,Bh,Hs,Mt,Ds,Rg,Cn,Nh,Fl,Mc,Lv,Ts,Og
      throw("ERROR: basissets/slater.hpp: SlaterBSet: GetExponents(): orbital exponents not defined for atom");
    }
    return exponents;
  }
  matrix<size_t> QuantumNumbers(size_t atomnum) {
    //getter for quantum numbers of one atom
    matrix<size_t> quantumnumb(1,2);    //1s
    //H,He
    if (atomnum > 86) {                                 //unavailable: Fr,Ra,Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr,Rf,Db,Sg,Bh,Hs,Mt,Ds,Rg,Cn,Nh,Fl,Mc,Lv,Ts,Og
      throw("ERROR: basissets/slater.hpp: SlaterBSet: QuantumNumbers(): atom not defined");
    }
    quantumnumb(1,1) = 1;
    if (atomnum > 2) {                  //Li,Be,B,C,N,O,F,Ne
      quantumnumb.resize(3,2);
      quantumnumb(2,1) = 2;    //2s
      quantumnumb(3,1) = 2;    //2p
      quantumnumb(3,2) = 1;    //2p
    }
    if (atomnum > 10) {            //Na,Mg,Al,Si,P,S,Cl,Ar
      quantumnumb.resize(5,2);
      quantumnumb(4,1) = 3;    //3s
      quantumnumb(5,1) = 3;    //3p
      quantumnumb(5,2) = 1;    //3p
    }
    if (atomnum > 18) {            //K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr
      quantumnumb.resize(8,2);
      quantumnumb(6,1) = 3;    //3d
      quantumnumb(6,2) = 2;    //3d
      quantumnumb(7,1) = 4;    //4s
      quantumnumb(8,1) = 4;    //4p
      quantumnumb(8,2) = 1;    //4p
    }
    if (atomnum > 36) {           //Rb,Sr,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Sb,Te,I,Xe
      quantumnumb.resize(12,2);
      quantumnumb(9,1) = 4;    //4d
      quantumnumb(9,2) = 2;    //4d
      quantumnumb(10,1) = 4;   //4f
      quantumnumb(10,2) = 3;   //4f
      quantumnumb(11,1) = 5;   //5s
      quantumnumb(12,1) = 5;   //5p
      quantumnumb(12,2) = 1;   //5p
    }
    if (atomnum > 54) {           //Cs,Ba,La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi,Po,At,Rn
      quantumnumb.resize(16,2);
      quantumnumb(13,1) = 5;   //5d
      quantumnumb(13,2) = 2;   //5d
      quantumnumb(14,1) = 5;   //5f
      quantumnumb(14,3) = 3;   //5f
      quantumnumb(15,1) = 6;   //6s
      quantumnumb(16,1) = 6;   //6p
      quantumnumb(16,2) = 1;   //6p
    }
    return quantumnumb;
  }
  //direct STO stuff
  size_t GetNAOs(size_t atomnum) {
    //function to calculate the number of AOs used for a given atom
    size_t orbitals = AtomNAOs(atomnum);
    if (orbitals == 1) {
      if (GetExponents(atomnum)[0] == 0.00) {orbitals = 0;}
    }
    return orbitals;
  }
  size_t AtomNAOs(size_t atomnum) {
    //function that returns the number of AOs for a given atom
    size_t szexp = GetExponents(atomnum).size();
    size_t res = 1;                            //1s
    if (szexp == 3) {res = 5;}                 //2s,2p
    else if (szexp == 5) {res = 9;}            //3s,3p
    else if (szexp == 6) {res = 14;}           //3d
    else if (szexp == 8) {res = 18;}           //4s,4p
    else if (szexp == 9) {res = 23;}           //4d
    else if (szexp == 10) {res = 30;}          //4f
    else if (szexp == 12) {res = 34;}          //5s,5p
    else if (szexp == 13) {res = 39;}          //5d
    else if (szexp == 14) {res = 46;}          //5f
    else if (szexp == 16) {res = 50;}          //6s,6p
    else if (szexp == 17) {res = 55;}          //6d
    else if (szexp == 18) {res = 62;}          //6f
    else if (szexp == 20) {res = 66;}          //7s,7p
    else if (szexp == 21) {res = 71;}          //7d
    else if (szexp == 22) {res = 78;}          //7f
    else if (szexp == 24) {res = 82;}          //8s,8p
    else if (szexp == 25) {res = 87;}          //8d
    else if (szexp == 26) {res = 94;}          //8f
    return res;
  }
  //size_t CalcMaxAB(size_t dim1, size_t dim2) {
  //  //no longer know what this does, nor what it is good for... it is however used nowhere
  //  size_t result = 2;
  //  if (((dim1 == 2)&&(dim2 == 1))||((dim1 == 1)&&(dim2 == 2))) {result = 3;}
  //  else if (((dim1 == 3)&&(dim2 == 1))||((dim1 == 1)&&(dim2 == 3))||((dim1 == 4)&&(dim2 == 1))||((dim1 == 1)&&(dim2 == 4))) {result = 4;}
  //  else if (((dim1 == 5)&&(dim2 == 1))||((dim1 == 1)&&(dim2 == 5))) {result = 5;}
  //  else if ((dim1 == 2)&&(dim2 == 2)) {result = 4;}
  //  else if (((dim1 == 3)&&(dim2 == 2))||((dim1 == 2)&&(dim2 == 3))||((dim1 == 4)&&(dim2 == 2))||((dim1 == 2)&&(dim2 == 4))) {result = 5;}
  //  else if (((dim1 == 5)&&(dim2 == 2))||((dim1 == 2)&&(dim2 == 5))) {result = 6;}
  //  else if ((dim1 == 3)&&(dim2 == 3)) {result = 6;}
  //  else if ((dim1 == 4)&&(dim2 == 4)) {result = 6;}
  //  else if (((dim1 == 4)&&(dim2 == 3))||((dim1 == 3)&&(dim2 == 4))) {result = 6;}
  //  else if (((dim1 == 5)&&(dim2 == 3))||((dim1 == 3)&&(dim2 == 5))) {result = 7;}
  //  else if (((dim1 == 5)&&(dim2 == 4))||((dim1 == 4)&&(dim2 == 5))) {result = 7;}
  //  else if ((dim1 == 5)&&(dim2 == 5)) {result = 8;}
  //  return result;
  //}
};

#endif //_Slater_Basis_Set_
