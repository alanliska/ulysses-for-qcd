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

#ifndef _Atomic_Weight_Package_
#define _Atomic_Weight_Package_
#include <vector>
#include <sstream>
#include "../UnitConversion.hpp"

//descrition:
//atom package containing basic atomic information, like atomic weights, atomic radii, conversion between atomic numbers and atomic symbols (and vice-versa), 
//base octet rules and experimental heats of formation for the elements
/*TODO:
- convert if-else to switch, more efficient
*/

double Weight(size_t atomicnumber) {
  //function returning the atomic weights in kg/mol, averaged from the isotopic masses and compositions given by NIST
  double weight = 0.0;
  if (atomicnumber == 1)        {weight = 1.00794075;}      //H
  else if (atomicnumber == 2)   {weight = 4.00260193;}      //He
  else if (atomicnumber == 3)   {weight = 6.94003660;}      //Li
  else if (atomicnumber == 4)   {weight = 9.01218307;}      //Be
  else if (atomicnumber == 5)   {weight = 10.81102805;}     //B
  else if (atomicnumber == 6)   {weight = 12.01073590;}     //C
  else if (atomicnumber == 7)   {weight = 14.00670321;}     //N
  else if (atomicnumber == 8)   {weight = 15.99940492;}     //O
  else if (atomicnumber == 9)   {weight = 18.99840316;}     //F
  else if (atomicnumber == 10)  {weight = 20.18004638;}     //Ne
  else if (atomicnumber == 11)  {weight = 22.98976928;}     //Na
  else if (atomicnumber == 12)  {weight = 24.30505162;}     //Mg
  else if (atomicnumber == 13)  {weight = 26.98153853;}     //Al
  else if (atomicnumber == 14)  {weight = 28.08549871;}     //Si
  else if (atomicnumber == 15)  {weight = 30.97376200;}     //P
  else if (atomicnumber == 16)  {weight = 32.06478741;}     //S
  else if (atomicnumber == 17)  {weight = 35.45293758;}     //Cl
  else if (atomicnumber == 18)  {weight = 39.94779856;}     //Ar
  else if (atomicnumber == 19)  {weight = 39.09830091;}     //K
  else if (atomicnumber == 20)  {weight = 40.07802251;}     //Ca
  else if (atomicnumber == 21)  {weight = 44.95590828;}     //Sc
  else if (atomicnumber == 22)  {weight = 47.86674496;}     //Ti
  else if (atomicnumber == 23)  {weight = 50.94146504;}     //V
  else if (atomicnumber == 24)  {weight = 51.99613176;}     //Cr
  else if (atomicnumber == 25)  {weight = 54.93804391;}     //Mn
  else if (atomicnumber == 26)  {weight = 55.84514443;}     //Fe
  else if (atomicnumber == 27)  {weight = 58.93319429;}     //Co
  else if (atomicnumber == 28)  {weight = 58.69334711;}     //Ni
  else if (atomicnumber == 29)  {weight = 63.54603995;}     //Cu
  else if (atomicnumber == 30)  {weight = 65.37778253;}     //Zn
  else if (atomicnumber == 31)  {weight = 69.72306607;}     //Ga
  else if (atomicnumber == 32)  {weight = 72.62755016;}     //Ge
  else if (atomicnumber == 33)  {weight = 74.92159457;}     //As
  else if (atomicnumber == 34)  {weight = 78.95938856;}     //Se
  else if (atomicnumber == 35)  {weight = 79.90352778;}     //Br
  else if (atomicnumber == 36)  {weight = 83.79800000;}     //Kr
  else if (atomicnumber == 37)  {weight = 85.46766360;}     //Rb
  else if (atomicnumber == 38)  {weight = 87.61664447;}     //Sr
  else if (atomicnumber == 39)  {weight = 88.90584030;}     //Y
  else if (atomicnumber == 40)  {weight = 91.22364160;}     //Zr
  else if (atomicnumber == 41)  {weight = 92.90637300;}     //Nb
  else if (atomicnumber == 42)  {weight = 95.95978854;}     //Mo
  else if (atomicnumber == 43)  {weight = 97.90721240;}     //Tc
  else if (atomicnumber == 44)  {weight = 101.06494014;}    //Ru
  else if (atomicnumber == 45)  {weight = 102.90549800;}    //Rh
  else if (atomicnumber == 46)  {weight = 106.41532751;}    //Pd
  else if (atomicnumber == 47)  {weight = 107.86814963;}    //Ag
  else if (atomicnumber == 48)  {weight = 112.41155782;}    //Cd
  else if (atomicnumber == 49)  {weight = 114.81808663;}    //In
  else if (atomicnumber == 50)  {weight = 118.71011259;}    //Sn
  else if (atomicnumber == 51)  {weight = 121.75978367;}    //Sb
  else if (atomicnumber == 52)  {weight = 127.60312648;}    //Te
  else if (atomicnumber == 53)  {weight = 126.90447190;}    //I
  else if (atomicnumber == 54)  {weight = 131.29276145;}    //Xe
  else if (atomicnumber == 55)  {weight = 132.90545196;}    //Cs
  else if (atomicnumber == 56)  {weight = 137.32689163;}    //Ba
  else if (atomicnumber == 57)  {weight = 138.90546887;}    //La
  else if (atomicnumber == 58)  {weight = 140.11573074;}    //Ce
  else if (atomicnumber == 59)  {weight = 140.90765760;}    //Pr
  else if (atomicnumber == 60)  {weight = 144.24159603;}    //Nd
  else if (atomicnumber == 61)  {weight = 144.91275590;}    //Pm
  else if (atomicnumber == 62)  {weight = 150.36635571;}    //Sm
  else if (atomicnumber == 63)  {weight = 151.96437813;}    //Eu
  else if (atomicnumber == 64)  {weight = 157.25213065;}    //Gd
  else if (atomicnumber == 65)  {weight = 158.92535470;}    //Tb
  else if (atomicnumber == 66)  {weight = 162.49947282;}    //Dy
  else if (atomicnumber == 67)  {weight = 164.93032880;}    //Ho
  else if (atomicnumber == 68)  {weight = 167.25908265;}    //Er
  else if (atomicnumber == 69)  {weight = 168.93421790;}    //Tm
  else if (atomicnumber == 70)  {weight = 173.05415017;}    //Yb
  else if (atomicnumber == 71)  {weight = 174.96681496;}    //Lu
  else if (atomicnumber == 72)  {weight = 178.48497872;}    //Hf
  else if (atomicnumber == 73)  {weight = 180.94787564;}    //Ta
  else if (atomicnumber == 74)  {weight = 183.84177755;}    //W
  else if (atomicnumber == 75)  {weight = 186.20670455;}    //Re
  else if (atomicnumber == 76)  {weight = 190.22485963;}    //Os
  else if (atomicnumber == 77)  {weight = 192.21605165;}    //Ir
  else if (atomicnumber == 78)  {weight = 195.08445686;}    //Pt
  else if (atomicnumber == 79)  {weight = 196.96656879;}    //Au
  else if (atomicnumber == 80)  {weight = 200.59916703;}    //Hg
  else if (atomicnumber == 81)  {weight = 204.38341284;}    //Tl
  else if (atomicnumber == 82)  {weight = 207.21690806;}    //Pb
  else if (atomicnumber == 83)  {weight = 208.98039910;}    //Bi
  else if (atomicnumber == 84)  {weight = 208.98243080;}    //Po
  else if (atomicnumber == 85)  {weight = 209.98714790;}    //At
  else if (atomicnumber == 86)  {weight = 222.01757820;}    //Rn
  else if (atomicnumber == 87)  {weight = 223.01973600;}    //Fr
  else if (atomicnumber == 88)  {weight = 226.02541030;}    //Ra
  else if (atomicnumber == 89)  {weight = 227.02775230;}    //Ac
  else if (atomicnumber == 90)  {weight = 232.03805580;}    //Th
  else if (atomicnumber == 91)  {weight = 231.03588420;}    //Pa
  else if (atomicnumber == 92)  {weight = 238.02891046;}    //U
  else if (atomicnumber == 93)  {weight = 237.04817360;}    //Np
  else if (atomicnumber == 94)  {weight = 244.06420530;}    //Pu
  else if (atomicnumber == 95)  {weight = 243.06138130;}    //Am
  else if (atomicnumber == 96)  {weight = 247.07035410;}    //Cm
  else if (atomicnumber == 97)  {weight = 247.07030730;}    //Bk
  else if (atomicnumber == 98)  {weight = 251.07958860;}    //Cf
  else if (atomicnumber == 99)  {weight = 252.08298000;}    //Es
  else if (atomicnumber == 100) {weight = 257.09510610;}    //Fm
  else if (atomicnumber == 101) {weight = 258.09843150;}    //Md
  else if (atomicnumber == 102) {weight = 259.10103000;}    //No
  else if (atomicnumber == 103) {weight = 262.10961000;}    //Lr
  else if (atomicnumber == 104) {weight = 267.12179000;}    //Rf
  else if (atomicnumber == 105) {weight = 269.12791000;}    //Db
  else if (atomicnumber == 106) {weight = 271.13393000;}    //Sg
  else if (atomicnumber == 107) {weight = 270.13336000;}    //Bh
  else if (atomicnumber == 108) {weight = 276.14846000;}    //Hs
  else if (atomicnumber == 109) {weight = 276.15159000;}    //Mt
  else if (atomicnumber == 110) {weight = 280.16131000;}    //Ds
  else if (atomicnumber == 111) {weight = 282.16912000;}    //Rg
  else if (atomicnumber == 112) {weight = 284.17416000;}    //Cn
  else if (atomicnumber == 113) {weight = 284.17873000;}    //Nh
  else if (atomicnumber == 114) {weight = 289.19042000;}    //Fl
  else if (atomicnumber == 115) {weight = 288.19274000;}    //Mc
  else if (atomicnumber == 116) {weight = 293.20449000;}    //Lv
  else if (atomicnumber == 117) {weight = 292.20746000;}    //Ts
  else if (atomicnumber == 118) {weight = 294.21392000;}    //Og
  return weight*1.0e-3;
}
double Radius(size_t atomicnumber) {
  //values stored and returned in meters, no idea about the source anymore
  double radius = 0.0;
  if (atomicnumber == 1) {radius = 0.53e-10;}               //H
  else if (atomicnumber == 2) {radius = 0.31e-10;}          //He
  else if (atomicnumber == 3) {radius = 1.67e-10;}          //Li
  else if (atomicnumber == 4) {radius = 1.12e-10;}          //Be
  else if (atomicnumber == 5) {radius = 0.87e-10;}          //B
  else if (atomicnumber == 6) {radius = 0.67e-10;}          //C
  else if (atomicnumber == 7) {radius = 0.56e-10;}          //N
  else if (atomicnumber == 8) {radius = 0.48e-10;}          //O
  else if (atomicnumber == 9) {radius = 0.42e-10;}          //F
  else if (atomicnumber == 10) {radius = 0.38e-10;}         //Ne
  else if (atomicnumber == 11) {radius = 1.9e-10;}          //Na
  else if (atomicnumber == 12) {radius = 1.45e-10;}         //Mg
  else if (atomicnumber == 13) {radius = 1.18e-10;}         //Al
  else if (atomicnumber == 14) {radius = 1.11e-10;}         //Si
  else if (atomicnumber == 15) {radius = 0.98e-10;}         //P
  else if (atomicnumber == 16) {radius = 0.88e-10;}         //S
  else if (atomicnumber == 17) {radius = 0.79e-10;}         //Cl
  else if (atomicnumber == 18) {radius = 0.71e-10;}         //Ar
  else if (atomicnumber == 19) {radius = 2.43e-10;}         //K
  else if (atomicnumber == 20) {radius = 1.94e-10;}         //Ca
  else if (atomicnumber == 21) {radius = 1.84e-10;}         //Sc
  else if (atomicnumber == 22) {radius = 1.76e-10;}         //Ti
  else if (atomicnumber == 23) {radius = 1.71e-10;}         //V
  else if (atomicnumber == 24) {radius = 1.66e-10;}         //Cr
  else if (atomicnumber == 25) {radius = 1.61e-10;}         //Mn
  else if (atomicnumber == 26) {radius = 1.56e-10;}         //Fe
  else if (atomicnumber == 27) {radius = 1.52e-10;}         //Co
  else if (atomicnumber == 28) {radius = 1.49e-10;}         //Ni
  else if (atomicnumber == 29) {radius = 1.45e-10;}         //Cu
  else if (atomicnumber == 30) {radius = 1.42e-10;}         //Zn
  else if (atomicnumber == 31) {radius = 1.36e-10;}         //Ga
  else if (atomicnumber == 32) {radius = 1.25e-10;}         //Ge
  else if (atomicnumber == 33) {radius = 1.14e-10;}         //As
  else if (atomicnumber == 34) {radius = 1.03e-10;}         //Se
  else if (atomicnumber == 35) {radius = 0.94e-10;}         //Br
  else if (atomicnumber == 36) {radius = 0.88e-10;}         //Kr
  else if (atomicnumber == 37) {radius = 2.65e-10;}         //Rb
  else if (atomicnumber == 38) {radius = 2.19e-10;}         //Sr
  else if (atomicnumber == 39) {radius = 2.12e-10;}         //Y
  else if (atomicnumber == 40) {radius = 2.06e-10;}         //Zr
  else if (atomicnumber == 41) {radius = 1.98e-10;}         //Nb
  else if (atomicnumber == 42) {radius = 1.9e-10;}          //Mo
  else if (atomicnumber == 43) {radius = 1.83e-10;}         //Tc
  else if (atomicnumber == 44) {radius = 1.78e-10;}         //Ru
  else if (atomicnumber == 45) {radius = 1.73e-10;}         //Rh
  else if (atomicnumber == 46) {radius = 1.69e-10;}         //Pd
  else if (atomicnumber == 47) {radius = 1.65e-10;}         //Ag
  else if (atomicnumber == 48) {radius = 1.61e-10;}         //Cd
  else if (atomicnumber == 49) {radius = 1.56e-10;}         //In
  else if (atomicnumber == 50) {radius = 1.45e-10;}         //Sn
  else if (atomicnumber == 51) {radius = 1.33e-10;}         //Sb
  else if (atomicnumber == 52) {radius = 1.23e-10;}         //Te
  else if (atomicnumber == 53) {radius = 1.15e-10;}         //I
  else if (atomicnumber == 54) {radius = 1.08e-10;}         //Xe
  else if (atomicnumber == 55) {radius = 2.98e-10;}         //Cs
  else if (atomicnumber == 56) {radius = 2.53e-10;}         //Ba
  else if (atomicnumber == 57) {radius = 1.95e-10;}         //La
  else if (atomicnumber == 58) {radius = 1.85e-10;}         //Ce
  else if (atomicnumber == 59) {radius = 2.47e-10;}         //Pr
  else if (atomicnumber == 60) {radius = 2.06e-10;}         //Nd
  else if (atomicnumber == 61) {radius = 2.05e-10;}         //Pm
  else if (atomicnumber == 62) {radius = 2.38e-10;}         //Sm
  else if (atomicnumber == 63) {radius = 2.31e-10;}         //Eu
  else if (atomicnumber == 64) {radius = 2.33e-10;}         //Gd
  else if (atomicnumber == 65) {radius = 2.25e-10;}         //Tb
  else if (atomicnumber == 66) {radius = 2.28e-10;}         //Dy
  else if (atomicnumber == 67) {radius = 2.26e-10;}         //Ho
  else if (atomicnumber == 68) {radius = 2.26e-10;}         //Er
  else if (atomicnumber == 69) {radius = 2.22e-10;}         //Tm
  else if (atomicnumber == 70) {radius = 2.22e-10;}         //Yb
  else if (atomicnumber == 71) {radius = 2.17e-10;}         //Lu
  else if (atomicnumber == 72) {radius = 2.08e-10;}         //Hf
  else if (atomicnumber == 73) {radius = 2.0e-10;}          //Ta
  else if (atomicnumber == 74) {radius = 1.93e-10;}         //W
  else if (atomicnumber == 75) {radius = 1.88e-10;}         //Re
  else if (atomicnumber == 76) {radius = 1.85e-10;}         //Os
  else if (atomicnumber == 77) {radius = 1.8e-10;}          //Ir
  else if (atomicnumber == 78) {radius = 1.77e-10;}         //Pt
  else if (atomicnumber == 79) {radius = 1.74e-10;}         //Au
  else if (atomicnumber == 80) {radius = 1.71e-10;}         //Hg
  else if (atomicnumber == 81) {radius = 1.56e-10;}         //Tl
  else if (atomicnumber == 82) {radius = 1.54e-10;}         //Pb
  else if (atomicnumber == 83) {radius = 1.43e-10;}         //Bi
  else if (atomicnumber == 84) {radius = 1.35e-10;}         //Po
  else if (atomicnumber == 85) {radius = 1.27e-10;}         //At
  else if (atomicnumber == 86) {radius = 1.2e-10;}          //Rn
  else if (atomicnumber == 87) {radius = 1.94e-10;}         //Fr
  else if (atomicnumber == 88) {radius = 1.62e-10;}         //Ra
  else if (atomicnumber == 89) {radius = 1.95e-10;}         //Ac
  else if (atomicnumber == 90) {radius = 1.8e-10;}          //Th
  else if (atomicnumber == 91) {radius = 1.8e-10;}          //Pa
  else if (atomicnumber == 92) {radius = 1.75e-10;}         //U
  else if (atomicnumber == 93) {radius = 1.75e-10;}         //Np
  else if (atomicnumber == 94) {radius = 1.75e-10;}         //Pu
  else if (atomicnumber == 95) {radius = 1.75e-10;}         //Am
  else if (atomicnumber == 96) {radius = 1.11e-10;}         //Cm
  //Unavailable Radii
  //else if (atomicnumber == 97) {radius = 0.e-3;}          //Bk
  //else if (atomicnumber == 98) {radius = 0.e-3;}          //Cf
  //else if (atomicnumber == 99) {radius = 0.e-3;}          //Es
  //else if (atomicnumber == 100) {radius = 0.e-3;}         //Fm
  //else if (atomicnumber == 101) {radius = 0.e-3;}         //Md
  //else if (atomicnumber == 102) {radius = 0.e-3;}         //No
  //else if (atomicnumber == 103) {radius = 0.e-3;}         //Lr
  //else if (atomicnumber == 104) {radius = 0.e-3;}         //Rf
  //else if (atomicnumber == 105) {radius = 0.e-3;}         //Db
  //else if (atomicnumber == 106) {radius = 0.e-3;}         //Sg
  //else if (atomicnumber == 107) {radius = 0.e-3;}         //Bh
  //else if (atomicnumber == 108) {radius = 0.e-3;}         //Hs
  //else if (atomicnumber == 109) {radius = 0.e-3;}         //Mt
  //else if (atomicnumber == 110) {radius = 0.e-3;}         //Ds
  //else if (atomicnumber == 111) {radius = 0.e-3;}         //Rg
  //else if (atomicnumber == 112) {radius = 0.e-3;}         //Cn
  //else if (atomicnumber == 113) {radius = 0.e-3;}         //Nh
  //else if (atomicnumber == 114) {radius = 0.e-3;}         //Fl
  //else if (atomicnumber == 115) {radius = 0.e-3;}         //Mc
  //else if (atomicnumber == 116) {radius = 0.e-3;}         //Lv
  //else if (atomicnumber == 117) {radius = 0.e-3;}         //Ts
  //else if (atomicnumber == 118) {radius = 0.e-3;}         //Og
  return radius;
}
double RadiusPsi4(size_t atomicnumber) {
  //these are atomic radii as used in Psi4, from
  //B. Cordero, V. Gómez, A. E. Platero-Prats, M. Revés, J. Echeverría, E. Cremades, F. Barragán and S. Alvarez, Dalton Trans., 2832, 2008
  //values stored in angstrom
  double radius;
  if (atomicnumber == 1) {radius = 0.31;}               //H
  else if (atomicnumber == 2) {radius = 0.28;}          //He
  else if (atomicnumber == 3) {radius = 1.28;}          //Li
  else if (atomicnumber == 4) {radius = 0.96;}          //Be
  else if (atomicnumber == 5) {radius = 0.84;}          //B
  else if (atomicnumber == 6) {radius = 0.76;}          //C
  else if (atomicnumber == 7) {radius = 0.71;}          //N
  else if (atomicnumber == 8) {radius = 0.66;}          //O
  else if (atomicnumber == 9) {radius = 0.57;}          //F
  else if (atomicnumber == 10) {radius = 0.58;}         //Ne
  else if (atomicnumber == 11) {radius = 1.66;}         //Na
  else if (atomicnumber == 12) {radius = 1.41;}         //Mg
  else if (atomicnumber == 13) {radius = 1.21;}         //Al
  else if (atomicnumber == 14) {radius = 1.11;}         //Si
  else if (atomicnumber == 15) {radius = 1.07;}         //P
  else if (atomicnumber == 16) {radius = 1.05;}         //S
  else if (atomicnumber == 17) {radius = 1.02;}         //Cl
  else if (atomicnumber == 18) {radius = 1.06;}         //Ar
  else if (atomicnumber == 19) {radius = 2.03;}         //K
  else if (atomicnumber == 20) {radius = 1.76;}         //Ca
  else if (atomicnumber == 21) {radius = 1.70;}         //Sc
  else if (atomicnumber == 22) {radius = 1.60;}         //Ti
  else if (atomicnumber == 23) {radius = 1.53;}         //V
  else if (atomicnumber == 24) {radius = 1.39;}         //Cr
  else if (atomicnumber == 25) {radius = 1.61;}         //Mn
  else if (atomicnumber == 26) {radius = 1.52;}         //Fe
  else if (atomicnumber == 27) {radius = 1.50;}         //Co
  else if (atomicnumber == 28) {radius = 1.24;}         //Ni
  else if (atomicnumber == 29) {radius = 1.32;}         //Cu
  else if (atomicnumber == 30) {radius = 1.22;}         //Zn
  else if (atomicnumber == 31) {radius = 1.22;}         //Ga
  else if (atomicnumber == 32) {radius = 1.20;}         //Ge
  else if (atomicnumber == 33) {radius = 1.19;}         //As
  else if (atomicnumber == 34) {radius = 1.20;}         //Se
  else if (atomicnumber == 35) {radius = 1.20;}         //Br
  else if (atomicnumber == 36) {radius = 1.16;}         //Kr
  else if (atomicnumber == 37) {radius = 2.20;}         //Rb
  else if (atomicnumber == 38) {radius = 1.95;}         //Sr
  else if (atomicnumber == 39) {radius = 1.90;}         //Y
  else if (atomicnumber == 40) {radius = 1.75;}         //Zr
  else if (atomicnumber == 41) {radius = 1.64;}         //Nb
  else if (atomicnumber == 42) {radius = 1.54;}         //Mo
  else if (atomicnumber == 43) {radius = 1.47;}         //Tc
  else if (atomicnumber == 44) {radius = 1.46;}         //Ru
  else if (atomicnumber == 45) {radius = 1.42;}         //Rh
  else if (atomicnumber == 46) {radius = 1.39;}         //Pd
  else if (atomicnumber == 47) {radius = 1.45;}         //Ag
  else if (atomicnumber == 48) {radius = 1.44;}         //Cd
  else if (atomicnumber == 49) {radius = 1.42;}         //In
  else if (atomicnumber == 50) {radius = 1.39;}         //Sn
  else if (atomicnumber == 51) {radius = 1.39;}         //Sb
  else if (atomicnumber == 52) {radius = 1.38;}         //Te
  else if (atomicnumber == 53) {radius = 1.39;}         //I
  else if (atomicnumber == 54) {radius = 1.40;}         //Xe
  else if (atomicnumber == 55) {radius = 2.44;}         //Cs
  else if (atomicnumber == 56) {radius = 2.15;}         //Ba
  else if (atomicnumber == 57) {radius = 2.07;}         //La
  else if (atomicnumber == 58) {radius = 2.04;}         //Ce
  else if (atomicnumber == 59) {radius = 2.03;}         //Pr
  else if (atomicnumber == 60) {radius = 2.01;}         //Nd
  else if (atomicnumber == 61) {radius = 1.99;}         //Pm
  else if (atomicnumber == 62) {radius = 1.98;}         //Sm
  else if (atomicnumber == 63) {radius = 1.98;}         //Eu
  else if (atomicnumber == 64) {radius = 1.96;}         //Gd
  else if (atomicnumber == 65) {radius = 1.94;}         //Tb
  else if (atomicnumber == 66) {radius = 1.92;}         //Dy
  else if (atomicnumber == 67) {radius = 1.92;}         //Ho
  else if (atomicnumber == 68) {radius = 1.89;}         //Er
  else if (atomicnumber == 69) {radius = 1.90;}         //Tm
  else if (atomicnumber == 70) {radius = 1.87;}         //Yb
  else if (atomicnumber == 71) {radius = 1.87;}         //Lu
  else if (atomicnumber == 72) {radius = 1.75;}         //Hf
  else if (atomicnumber == 73) {radius = 1.70;}         //Ta
  else if (atomicnumber == 74) {radius = 1.62;}         //W
  else if (atomicnumber == 75) {radius = 1.51;}         //Re
  else if (atomicnumber == 76) {radius = 1.44;}         //Os
  else if (atomicnumber == 77) {radius = 1.41;}         //Ir
  else if (atomicnumber == 78) {radius = 1.36;}         //Pt
  else if (atomicnumber == 79) {radius = 1.36;}         //Au
  else if (atomicnumber == 80) {radius = 1.32;}         //Hg
  else if (atomicnumber == 81) {radius = 1.45;}         //Tl
  else if (atomicnumber == 82) {radius = 1.46;}         //Pb
  else if (atomicnumber == 83) {radius = 1.48;}         //Bi
  else if (atomicnumber == 84) {radius = 1.40;}         //Po
  else if (atomicnumber == 85) {radius = 1.50;}         //At
  else if (atomicnumber == 86) {radius = 1.50;}         //Rn
  else if (atomicnumber == 87) {radius = 2.60;}         //Fr
  else if (atomicnumber == 88) {radius = 2.21;}         //Ra
  else if (atomicnumber == 89) {radius = 2.15;}         //Ac
  else if (atomicnumber == 90) {radius = 2.06;}         //Th
  else if (atomicnumber == 91) {radius = 2.00;}         //Pa
  else if (atomicnumber == 92) {radius = 1.96;}         //U
  else if (atomicnumber == 93) {radius = 1.90;}         //Np
  else if (atomicnumber == 94) {radius = 1.87;}         //Pu
  else if (atomicnumber == 95) {radius = 1.80;}         //Am
  else if (atomicnumber == 96) {radius = 1.69;}         //Cm
  //Unavailable Radii
  //else if (atomicnumber == 97) {radius = 0.;}          //Bk
  //else if (atomicnumber == 98) {radius = 0.;}          //Cf
  //else if (atomicnumber == 99) {radius = 0.;}          //Es
  //else if (atomicnumber == 100) {radius = 0.;}         //Fm
  //else if (atomicnumber == 101) {radius = 0.;}         //Md
  //else if (atomicnumber == 102) {radius = 0.;}         //No
  //else if (atomicnumber == 103) {radius = 0.;}         //Lr
  //else if (atomicnumber == 104) {radius = 0.;}         //Rf
  //else if (atomicnumber == 105) {radius = 0.;}         //Db
  //else if (atomicnumber == 106) {radius = 0.;}         //Sg
  //else if (atomicnumber == 107) {radius = 0.;}         //Bh
  //else if (atomicnumber == 108) {radius = 0.;}         //Hs
  //else if (atomicnumber == 109) {radius = 0.;}         //Mt
  //else if (atomicnumber == 110) {radius = 0.;}         //Ds
  //else if (atomicnumber == 111) {radius = 0.;}         //Rg
  //else if (atomicnumber == 112) {radius = 0.;}         //Cn
  //else if (atomicnumber == 113) {radius = 0.;}         //Nh
  //else if (atomicnumber == 114) {radius = 0.;}         //Fl
  //else if (atomicnumber == 115) {radius = 0.;}         //Mc
  //else if (atomicnumber == 116) {radius = 0.;}         //Lv
  //else if (atomicnumber == 117) {radius = 0.;}         //Ts
  //else if (atomicnumber == 118) {radius = 0.;}         //Og
  return radius;
}
int Symbol2AtomNr(std::string symbol) {
  int atomicnr = 0;
  if      ((symbol == "h")||(symbol == "H"))   {atomicnr = 1;}
  else if ((symbol == "he")||(symbol == "He")) {atomicnr = 2;}
  else if ((symbol == "li")||(symbol == "Li")) {atomicnr = 3;}
  else if ((symbol == "be")||(symbol == "Be")) {atomicnr = 4;}
  else if ((symbol == "b")||(symbol == "B"))   {atomicnr = 5;}
  else if ((symbol == "c")||(symbol == "C"))   {atomicnr = 6;}
  else if ((symbol == "n")||(symbol == "N"))   {atomicnr = 7;}
  else if ((symbol == "o")||(symbol == "O"))   {atomicnr = 8;}
  else if ((symbol == "f")||(symbol == "F"))   {atomicnr = 9;}
  else if ((symbol == "ne")||(symbol == "Ne")) {atomicnr = 10;}
  else if ((symbol == "na")||(symbol == "Na")) {atomicnr = 11;}
  else if ((symbol == "mg")||(symbol == "Mg")) {atomicnr = 12;}
  else if ((symbol == "al")||(symbol == "Al")) {atomicnr = 13;}
  else if ((symbol == "si")||(symbol == "Si")) {atomicnr = 14;}
  else if ((symbol == "p")||(symbol == "P"))   {atomicnr = 15;}
  else if ((symbol == "s")||(symbol == "S"))   {atomicnr = 16;}
  else if ((symbol == "cl")||(symbol == "Cl")) {atomicnr = 17;}
  else if ((symbol == "ar")||(symbol == "Ar")) {atomicnr = 18;}
  else if ((symbol == "k")||(symbol == "K"))   {atomicnr = 19;}
  else if ((symbol == "ca")||(symbol == "Ca")) {atomicnr = 20;}
  else if ((symbol == "sc")||(symbol == "Sc")) {atomicnr = 21;}
  else if ((symbol == "ti")||(symbol == "Ti")) {atomicnr = 22;}
  else if ((symbol == "v")||(symbol == "V"))   {atomicnr = 23;}
  else if ((symbol == "cr")||(symbol == "Cr")) {atomicnr = 24;}
  else if ((symbol == "mn")||(symbol == "Mn")) {atomicnr = 25;}
  else if ((symbol == "fe")||(symbol == "Fe")) {atomicnr = 26;}
  else if ((symbol == "co")||(symbol == "Co")) {atomicnr = 27;}
  else if ((symbol == "ni")||(symbol == "Ni")) {atomicnr = 28;}
  else if ((symbol == "cu")||(symbol == "Cu")) {atomicnr = 29;}
  else if ((symbol == "zn")||(symbol == "Zn")) {atomicnr = 30;}
  else if ((symbol == "ga")||(symbol == "Ga")) {atomicnr = 31;}
  else if ((symbol == "ge")||(symbol == "Ge")) {atomicnr = 32;}
  else if ((symbol == "as")||(symbol == "As")) {atomicnr = 33;}
  else if ((symbol == "se")||(symbol == "Se")) {atomicnr = 34;}
  else if ((symbol == "br")||(symbol == "Br")) {atomicnr = 35;}
  else if ((symbol == "kr")||(symbol == "Kr")) {atomicnr = 36;}
  else if ((symbol == "rb")||(symbol == "Rb")) {atomicnr = 37;}
  else if ((symbol == "sr")||(symbol == "Sr")) {atomicnr = 38;}
  else if ((symbol == "y")||(symbol == "Y"))   {atomicnr = 39;}
  else if ((symbol == "zr")||(symbol == "Zr")) {atomicnr = 40;}
  else if ((symbol == "nb")||(symbol == "Nb")) {atomicnr = 41;}
  else if ((symbol == "mo")||(symbol == "Mo")) {atomicnr = 42;}
  else if ((symbol == "tc")||(symbol == "Tc")) {atomicnr = 43;}
  else if ((symbol == "ru")||(symbol == "Ru")) {atomicnr = 44;}
  else if ((symbol == "rh")||(symbol == "Rh")) {atomicnr = 45;}
  else if ((symbol == "pd")||(symbol == "Pd")) {atomicnr = 46;}
  else if ((symbol == "ag")||(symbol == "Ag")) {atomicnr = 47;}
  else if ((symbol == "cd")||(symbol == "Cd")) {atomicnr = 48;}
  else if ((symbol == "in")||(symbol == "In")) {atomicnr = 49;}
  else if ((symbol == "sn")||(symbol == "Sn")) {atomicnr = 50;}
  else if ((symbol == "sb")||(symbol == "Sb")) {atomicnr = 51;}
  else if ((symbol == "te")||(symbol == "Te")) {atomicnr = 52;}
  else if ((symbol == "i")||(symbol == "I"))   {atomicnr = 53;}
  else if ((symbol == "xe")||(symbol == "Xe")) {atomicnr = 54;}
  else if ((symbol == "cs")||(symbol == "Cs")) {atomicnr = 55;}
  else if ((symbol == "ba")||(symbol == "Ba")) {atomicnr = 56;}
  else if ((symbol == "la")||(symbol == "La")) {atomicnr = 57;}
  else if ((symbol == "ce")||(symbol == "Ce")) {atomicnr = 58;}
  else if ((symbol == "pr")||(symbol == "Pr")) {atomicnr = 59;}
  else if ((symbol == "nd")||(symbol == "Nd")) {atomicnr = 60;}
  else if ((symbol == "pm")||(symbol == "Pm")) {atomicnr = 61;}
  else if ((symbol == "sm")||(symbol == "Sm")) {atomicnr = 62;}
  else if ((symbol == "eu")||(symbol == "Eu")) {atomicnr = 63;}
  else if ((symbol == "gd")||(symbol == "Gd")) {atomicnr = 64;}
  else if ((symbol == "tb")||(symbol == "Tb")) {atomicnr = 65;}
  else if ((symbol == "dy")||(symbol == "Dy")) {atomicnr = 66;}
  else if ((symbol == "ho")||(symbol == "Ho")) {atomicnr = 67;}
  else if ((symbol == "er")||(symbol == "Er")) {atomicnr = 68;}
  else if ((symbol == "tm")||(symbol == "Tm")) {atomicnr = 69;}
  else if ((symbol == "yb")||(symbol == "Yb")) {atomicnr = 70;}
  else if ((symbol == "lu")||(symbol == "Lu")) {atomicnr = 71;}
  else if ((symbol == "hf")||(symbol == "Hf")) {atomicnr = 72;}
  else if ((symbol == "ta")||(symbol == "Ta")) {atomicnr = 73;}
  else if ((symbol == "w")||(symbol == "W"))   {atomicnr = 74;}
  else if ((symbol == "re")||(symbol == "Re")) {atomicnr = 75;}
  else if ((symbol == "os")||(symbol == "Os")) {atomicnr = 76;}
  else if ((symbol == "ir")||(symbol == "Ir")) {atomicnr = 77;}
  else if ((symbol == "pt")||(symbol == "Pt")) {atomicnr = 78;}
  else if ((symbol == "au")||(symbol == "Au")) {atomicnr = 79;}
  else if ((symbol == "hg")||(symbol == "Hg")) {atomicnr = 80;}
  else if ((symbol == "tl")||(symbol == "Tl")) {atomicnr = 81;}
  else if ((symbol == "pb")||(symbol == "Pb")) {atomicnr = 82;}
  else if ((symbol == "bi")||(symbol == "Bi")) {atomicnr = 83;}
  else if ((symbol == "po")||(symbol == "Po")) {atomicnr = 84;}
  else if ((symbol == "at")||(symbol == "At")) {atomicnr = 85;}
  else if ((symbol == "rn")||(symbol == "Rn")) {atomicnr = 86;}
  else if ((symbol == "fr")||(symbol == "Fr")) {atomicnr = 87;}
  else if ((symbol == "ra")||(symbol == "Ra")) {atomicnr = 88;}
  else if ((symbol == "ac")||(symbol == "Ac")) {atomicnr = 89;}
  else if ((symbol == "th")||(symbol == "Th")) {atomicnr = 90;}
  else if ((symbol == "pa")||(symbol == "Pa")) {atomicnr = 91;}
  else if ((symbol == "u")||(symbol == "U"))   {atomicnr = 92;}
  else if ((symbol == "np")||(symbol == "Np")) {atomicnr = 93;}
  else if ((symbol == "pu")||(symbol == "Pu")) {atomicnr = 94;}
  else if ((symbol == "am")||(symbol == "Am")) {atomicnr = 95;}
  else if ((symbol == "cm")||(symbol == "Cm")) {atomicnr = 96;}
  else if ((symbol == "bk")||(symbol == "Bk")) {atomicnr = 97;}
  else if ((symbol == "cf")||(symbol == "Cf")) {atomicnr = 98;}
  else if ((symbol == "es")||(symbol == "Es")) {atomicnr = 99;}
  else if ((symbol == "fm")||(symbol == "Fm")) {atomicnr = 100;}
  else if ((symbol == "md")||(symbol == "Md")) {atomicnr = 101;}
  else if ((symbol == "no")||(symbol == "No")) {atomicnr = 102;}
  else if ((symbol == "lr")||(symbol == "Lr")) {atomicnr = 103;}
  else if ((symbol == "rf")||(symbol == "Rf")) {atomicnr = 104;}
  else if ((symbol == "db")||(symbol == "Db")) {atomicnr = 105;}
  else if ((symbol == "sg")||(symbol == "Sg")) {atomicnr = 106;}
  else if ((symbol == "bh")||(symbol == "Bh")) {atomicnr = 107;}
  else if ((symbol == "hs")||(symbol == "Hs")) {atomicnr = 108;}
  else if ((symbol == "mt")||(symbol == "Mt")) {atomicnr = 109;}
  else if ((symbol == "ds")||(symbol == "Ds")) {atomicnr = 110;}
  else if ((symbol == "rg")||(symbol == "Rg")) {atomicnr = 111;}
  else if ((symbol == "cn")||(symbol == "Cn")) {atomicnr = 112;}
  else if ((symbol == "nh")||(symbol == "Nh")) {atomicnr = 113;}
  else if ((symbol == "fl")||(symbol == "Fl")) {atomicnr = 114;}
  else if ((symbol == "mc")||(symbol == "Mc")) {atomicnr = 115;}
  else if ((symbol == "lv")||(symbol == "Lv")) {atomicnr = 116;}
  else if ((symbol == "ts")||(symbol == "Ts")) {atomicnr = 117;}
  else if ((symbol == "og")||(symbol == "Og")) {atomicnr = 118;}
  return atomicnr;
}
std::string AtomNr2Symbol(int atomicnr) {
  std::string symbol;
  if      (atomicnr == 0)   {symbol = "X";}       //Pseudoatom or pseudoparticle
  else if (atomicnr == 1)   {symbol = "H";}
  else if (atomicnr == 2)   {symbol = "He";}
  else if (atomicnr == 3)   {symbol = "Li";}
  else if (atomicnr == 4)   {symbol = "Be";}
  else if (atomicnr == 5)   {symbol = "B";}
  else if (atomicnr == 6)   {symbol = "C";}
  else if (atomicnr == 7)   {symbol = "N";}
  else if (atomicnr == 8)   {symbol = "O";}
  else if (atomicnr == 9)   {symbol = "F";}
  else if (atomicnr == 10)  {symbol = "Ne";}
  else if (atomicnr == 11)  {symbol = "Na";}
  else if (atomicnr == 12)  {symbol = "Mg";}
  else if (atomicnr == 13)  {symbol = "Al";}
  else if (atomicnr == 14)  {symbol = "Si";}
  else if (atomicnr == 15)  {symbol = "P";}
  else if (atomicnr == 16)  {symbol = "S";}
  else if (atomicnr == 17)  {symbol = "Cl";}
  else if (atomicnr == 18)  {symbol = "Ar";}
  else if (atomicnr == 19)  {symbol = "K";}
  else if (atomicnr == 20)  {symbol = "Ca";}
  else if (atomicnr == 21)  {symbol = "Sc";}
  else if (atomicnr == 22)  {symbol = "Ti";}
  else if (atomicnr == 23)  {symbol = "V";}
  else if (atomicnr == 24)  {symbol = "Cr";}
  else if (atomicnr == 25)  {symbol = "Mn";}
  else if (atomicnr == 26)  {symbol = "Fe";}
  else if (atomicnr == 27)  {symbol = "Co";}
  else if (atomicnr == 28)  {symbol = "Ni";}
  else if (atomicnr == 29)  {symbol = "Cu";}
  else if (atomicnr == 30)  {symbol = "Zn";}
  else if (atomicnr == 31)  {symbol = "Ga";}
  else if (atomicnr == 32)  {symbol = "Ge";}
  else if (atomicnr == 33)  {symbol = "As";}
  else if (atomicnr == 34)  {symbol = "Se";}
  else if (atomicnr == 35)  {symbol = "Br";}
  else if (atomicnr == 36)  {symbol = "Kr";}
  else if (atomicnr == 37)  {symbol = "Rb";}
  else if (atomicnr == 38)  {symbol = "Sr";}
  else if (atomicnr == 39)  {symbol = "Y";}
  else if (atomicnr == 40)  {symbol = "Zr";}
  else if (atomicnr == 41)  {symbol = "Nb";}
  else if (atomicnr == 42)  {symbol = "Mo";}
  else if (atomicnr == 43)  {symbol = "Tc";}
  else if (atomicnr == 44)  {symbol = "Ru";}
  else if (atomicnr == 45)  {symbol = "Rh";}
  else if (atomicnr == 46)  {symbol = "Pd";}
  else if (atomicnr == 47)  {symbol = "Ag";}
  else if (atomicnr == 48)  {symbol = "Cd";}
  else if (atomicnr == 49)  {symbol = "In";}
  else if (atomicnr == 50)  {symbol = "Sn";}
  else if (atomicnr == 51)  {symbol = "Sb";}
  else if (atomicnr == 52)  {symbol = "Te";}
  else if (atomicnr == 53)  {symbol = "I";}
  else if (atomicnr == 54)  {symbol = "Xe";}
  else if (atomicnr == 55)  {symbol = "Cs";}
  else if (atomicnr == 56)  {symbol = "Ba";}
  else if (atomicnr == 57)  {symbol = "La";}
  else if (atomicnr == 58)  {symbol = "Ce";}
  else if (atomicnr == 59)  {symbol = "Pr";}
  else if (atomicnr == 60)  {symbol = "Nd";}
  else if (atomicnr == 61)  {symbol = "Pm";}
  else if (atomicnr == 62)  {symbol = "Sm";}
  else if (atomicnr == 63)  {symbol = "Eu";}
  else if (atomicnr == 64)  {symbol = "Gd";}
  else if (atomicnr == 65)  {symbol = "Tb";}
  else if (atomicnr == 66)  {symbol = "Dy";}
  else if (atomicnr == 67)  {symbol = "Ho";}
  else if (atomicnr == 68)  {symbol = "Er";}
  else if (atomicnr == 69)  {symbol = "Tm";}
  else if (atomicnr == 70)  {symbol = "Yb";}
  else if (atomicnr == 71)  {symbol = "Lu";}
  else if (atomicnr == 72)  {symbol = "Hf";}
  else if (atomicnr == 73)  {symbol = "Ta";}
  else if (atomicnr == 74)  {symbol = "W";}
  else if (atomicnr == 75)  {symbol = "Re";}
  else if (atomicnr == 76)  {symbol = "Os";}
  else if (atomicnr == 77)  {symbol = "Ir";}
  else if (atomicnr == 78)  {symbol = "Pt";}
  else if (atomicnr == 79)  {symbol = "Au";}
  else if (atomicnr == 80)  {symbol = "Hg";}
  else if (atomicnr == 81)  {symbol = "Tl";}
  else if (atomicnr == 82)  {symbol = "Pb";}
  else if (atomicnr == 83)  {symbol = "Bi";}
  else if (atomicnr == 84)  {symbol = "Po";}
  else if (atomicnr == 85)  {symbol = "At";}
  else if (atomicnr == 86)  {symbol = "Rn";}
  else if (atomicnr == 87)  {symbol = "Fr";}
  else if (atomicnr == 88)  {symbol = "Ra";}
  else if (atomicnr == 89)  {symbol = "Ac";}
  else if (atomicnr == 90)  {symbol = "Th";}
  else if (atomicnr == 91)  {symbol = "Pa";}
  else if (atomicnr == 92)  {symbol = "U";}
  else if (atomicnr == 93)  {symbol = "Np";}
  else if (atomicnr == 94)  {symbol = "Pu";}
  else if (atomicnr == 95)  {symbol = "Am";}
  else if (atomicnr == 96)  {symbol = "Cm";}
  else if (atomicnr == 97)  {symbol = "Bk";}
  else if (atomicnr == 98)  {symbol = "Cf";}
  else if (atomicnr == 99)  {symbol = "Es";}
  else if (atomicnr == 100) {symbol = "Fm";}
  else if (atomicnr == 101) {symbol = "Md";}
  else if (atomicnr == 102) {symbol = "No";}
  else if (atomicnr == 103) {symbol = "Lr";}
  else if (atomicnr == 104) {symbol = "Rf";}
  else if (atomicnr == 105) {symbol = "Db";}
  else if (atomicnr == 106) {symbol = "Sg";}
  else if (atomicnr == 107) {symbol = "Bh";}
  else if (atomicnr == 108) {symbol = "Hs";}
  else if (atomicnr == 109) {symbol = "Mt";}
  else if (atomicnr == 110) {symbol = "Ds";}
  else if (atomicnr == 111) {symbol = "Rg";}
  else if (atomicnr == 112) {symbol = "Cn";}
  else if (atomicnr == 113) {symbol = "Nh";}
  else if (atomicnr == 114) {symbol = "Fl";}
  else if (atomicnr == 115) {symbol = "Mc";}
  else if (atomicnr == 116) {symbol = "Lv";}
  else if (atomicnr == 117) {symbol = "Ts";}
  else if (atomicnr == 118) {symbol = "Og";}
  return symbol;
}
unsigned int octet(int atomicnr) {
  //atomicnr is atomic number
  unsigned int _octet = 0;
  //group 1
  if      ((atomicnr == 1)||(atomicnr == 3) ||(atomicnr == 11)||(atomicnr == 19)||(atomicnr == 37)||(atomicnr == 55)) {_octet = 1;}
  //group 2
  else if            ((atomicnr == 4) ||(atomicnr == 12)||(atomicnr == 20)||(atomicnr == 38)||(atomicnr == 56)) {_octet = 2;}
  //group 3
  else if            ((atomicnr == 5) ||(atomicnr == 13)||(atomicnr == 31)||(atomicnr == 49)||(atomicnr == 81)) {_octet = 4;}
  //group 4
  else if            ((atomicnr == 6) ||(atomicnr == 14)||(atomicnr == 32)||(atomicnr == 50)||(atomicnr == 82)) {_octet = 4;}
  //group 5, 6, 7
  //second-row
  else if            (atomicnr == 7) {_octet = 4;}
  else if            (atomicnr == 8) {_octet = 3;}
  else if            (atomicnr == 9) {_octet = 2;}
  //other
  else if            ((atomicnr == 15)||(atomicnr == 33)||(atomicnr == 51)||(atomicnr == 83)) {_octet = 6;}
  else if            ((atomicnr == 16)||(atomicnr == 34)||(atomicnr == 52)||(atomicnr == 84)) {_octet = 6;}
  else if            ((atomicnr == 17)||(atomicnr == 35)||(atomicnr == 53)||(atomicnr == 85)) {_octet = 6;}
  //group 8
  else if ((atomicnr == 2)||(atomicnr == 10)||(atomicnr == 18)||(atomicnr == 36)||(atomicnr == 54)||(atomicnr == 86)) {_octet = 0;}
  //transition metals
  else if ((atomicnr == 21)||(atomicnr == 22)||(atomicnr == 23)||(atomicnr == 24)||(atomicnr == 25)||(atomicnr == 26)||(atomicnr == 27)) {_octet = 6;}
  else if ((atomicnr == 28)||(atomicnr == 29)||(atomicnr == 30)||(atomicnr == 39)||(atomicnr == 40)||(atomicnr == 41)||(atomicnr == 42)) {_octet = 6;}
  else if ((atomicnr == 43)||(atomicnr == 44)||(atomicnr == 45)||(atomicnr == 46)||(atomicnr == 47)||(atomicnr == 48)||(atomicnr == 57)) {_octet = 6;}
  else if ((atomicnr == 71)||(atomicnr == 72)||(atomicnr == 73)||(atomicnr == 74)||(atomicnr == 75)||(atomicnr == 76)||(atomicnr == 77)) {_octet = 6;}
  else if ((atomicnr == 78)||(atomicnr == 79)||(atomicnr == 80)) {_octet = 6;}
  return _octet;
}
int ValencePairs(int atomicnr) {
  //function returning the total number of valence pairs for a certain atom
  int valence = 0;
  //transition metals
  if ((atomicnr > 20)&&(atomicnr < 31)) {valence = 18;}
  else if ((atomicnr > 38)&&(atomicnr < 49)) {valence = 18;}
  else if ((atomicnr > 56)&&(atomicnr < 81)) {valence = 18;}
  //normal dudes
  else if (atomicnr < 3) {valence = 2;}
  else {valence = 8;}
  return valence/2;
}
int OctetCalculator(int atomicnr, int charge) {
  //simple function that calculates the expected octet based on charge and atomic number
  //note that for s,p blocks this neglects extended octets
  //this is because the code currently only works for octets
  return 2*ValencePairs(atomicnr) - atomicnr + charge;
}
int Group(int atomicnr) {
  int group = 0;
  if      ((atomicnr == 1)||(atomicnr == 3) ||(atomicnr == 11)||(atomicnr == 19)||(atomicnr == 37)||(atomicnr == 55)) {group = 1;}
  else if ((atomicnr == 4)||(atomicnr == 12)||(atomicnr == 20)||(atomicnr == 38)||(atomicnr == 56)) {group = 2;}
  else if ((atomicnr == 5)||(atomicnr == 13)||(atomicnr == 31)||(atomicnr == 49)||(atomicnr == 81)) {group = 3;}
  else if ((atomicnr == 6)||(atomicnr == 14)||(atomicnr == 32)||(atomicnr == 50)||(atomicnr == 82)) {group = 4;}
  else if ((atomicnr == 7)||(atomicnr == 15)||(atomicnr == 33)||(atomicnr == 51)||(atomicnr == 83)) {group = 5;}
  else if ((atomicnr == 8)||(atomicnr == 16)||(atomicnr == 34)||(atomicnr == 52)||(atomicnr == 84)) {group = 6;}
  else if ((atomicnr == 9)||(atomicnr == 17)||(atomicnr == 35)||(atomicnr == 53)||(atomicnr == 85)) {group = 7;}
  else if ((atomicnr == 2)||(atomicnr == 10)||(atomicnr == 18)||(atomicnr == 36)||(atomicnr == 54)||(atomicnr == 86)) {group = 8;}
  //transition metals
  else if ((atomicnr == 21)||(atomicnr == 39)||(atomicnr == 57)||(atomicnr == 89)) {group = 3;}
  else if ((atomicnr == 22)||(atomicnr == 40)||(atomicnr == 72)||(atomicnr == 104)) {group = 4;}
  else if ((atomicnr == 23)||(atomicnr == 41)||(atomicnr == 73)||(atomicnr == 105)) {group = 5;}
  else if ((atomicnr == 24)||(atomicnr == 42)||(atomicnr == 74)||(atomicnr == 106)) {group = 6;}
  else if ((atomicnr == 25)||(atomicnr == 43)||(atomicnr == 75)||(atomicnr == 107)) {group = 7;}
  else if ((atomicnr == 26)||(atomicnr == 44)||(atomicnr == 76)||(atomicnr == 108)) {group = 8;}
  else if ((atomicnr == 27)||(atomicnr == 45)||(atomicnr == 77)||(atomicnr == 109)) {group = 9;}
  else if ((atomicnr == 28)||(atomicnr == 46)||(atomicnr == 78)||(atomicnr == 110)) {group = 10;}
  else if ((atomicnr == 29)||(atomicnr == 47)||(atomicnr == 79)||(atomicnr == 111)) {group = 11;}
  else if ((atomicnr == 30)||(atomicnr == 48)||(atomicnr == 80)||(atomicnr == 112)) {group = 12;}
  return group;
}
double ElementFormExp(double atomicnr) {
  //function returning experimental heats of formation for elements
  //values in kcal/mol, returned however in atomic units
  double enth = 0.0;
  if (atomicnr == 1) {enth = 52.102;}            //H
  else if (atomicnr == 2) {enth = 0.0000;}       //He
  else if (atomicnr == 3) {enth = 38.41;}        //Li
  else if (atomicnr == 4) {enth = 76.96;}        //Be
  else if (atomicnr == 5) {enth = 135.70;}       //B
  else if (atomicnr == 6) {enth = 170.89;}       //C
  else if (atomicnr == 7) {enth = 113.00;}       //N
  else if (atomicnr == 8) {enth = 59.559;}       //O
  else if (atomicnr == 9) {enth = 18.86;}        //F
  else if (atomicnr == 10) {enth = 0.000;}       //Ne
  else if (atomicnr == 11) {enth = 25.85;}       //Na
  else if (atomicnr == 12) {enth = 35.00;}       //Mg
  else if (atomicnr == 13) {enth = 79.49;}       //Al
  else if (atomicnr == 14) {enth = 108.39;}      //Si
  else if (atomicnr == 15) {enth = 75.57;}       //P
  else if (atomicnr == 16) {enth = 66.40;}       //S
  else if (atomicnr == 17) {enth = 28.99;}       //Cl
  else if (atomicnr == 18) {enth = 0.00;}        //Ar
  else if (atomicnr == 19) {enth = 21.42;}       //K
  else if (atomicnr == 20) {enth = 42.60;}       //Ca
  else if (atomicnr == 21) {enth = 90.30;}       //Sc
  else if (atomicnr == 22) {enth = 112.30;}      //Ti
  else if (atomicnr == 23) {enth = 122.90;}      //V
  else if (atomicnr == 24) {enth = 95.00;}       //Cr
  else if (atomicnr == 25) {enth = 67.70;}       //Mn
  else if (atomicnr == 26) {enth = 99.30;}       //Fe
  else if (atomicnr == 27) {enth = 102.40;}      //Co
  else if (atomicnr == 28) {enth = 102.80;}      //Ni
  else if (atomicnr == 29) {enth = 80.70;}       //Cu
  else if (atomicnr == 30) {enth = 31.17;}       //Zn
  else if (atomicnr == 31) {enth = 65.40;}       //Ga
  else if (atomicnr == 32) {enth = 89.50;}       //Ge
  else if (atomicnr == 33) {enth = 72.30;}       //As
  else if (atomicnr == 34) {enth = 54.30;}       //Se
  else if (atomicnr == 35) {enth = 26.74;}       //Br
  else if (atomicnr == 36) {enth = 0.000;}       //Kr
  else if (atomicnr == 37) {enth = 19.60;}       //Rb
  else if (atomicnr == 38) {enth = 39.10;}       //Sb
  else if (atomicnr == 39) {enth = 101.50;}      //Y
  else if (atomicnr == 40) {enth = 145.50;}      //Zr
  else if (atomicnr == 41) {enth = 172.40;}      //Nb
  else if (atomicnr == 42) {enth = 157.30;}      //Mo
  else if (atomicnr == 43) {enth = 72.30;}       //Tc
  else if (atomicnr == 44) {enth = 155.50;}      //Ru
  else if (atomicnr == 45) {enth = 133.00;}      //Rh
  else if (atomicnr == 46) {enth = 90.00;}       //Pd
  else if (atomicnr == 47) {enth = 68.10;}       //Ag
  else if (atomicnr == 48) {enth = 26.72;}       //Cd
  else if (atomicnr == 49) {enth = 58.00;}       //In
  else if (atomicnr == 50) {enth = 72.20;}       //Sn
  else if (atomicnr == 51) {enth = 63.20;}       //Sb
  else if (atomicnr == 52) {enth = 47.00;}       //Te
  else if (atomicnr == 53) {enth = 25.517;}      //I
  else if (atomicnr == 54) {enth = 0.000;}       //Xe
  else if (atomicnr == 55) {enth = 18.70;}       //Cs
  else if (atomicnr == 56) {enth = 42.50;}       //Ba
  else if (atomicnr == 57) {enth = 928.90;}      //La
  else if (atomicnr == 58) {enth = 101.30;}      //Ce
  else if (atomicnr == 62) {enth = 49.40;}       //Sm
  else if (atomicnr == 68) {enth = 75.80;}       //Er
  else if (atomicnr == 70) {enth = 36.35;}       //Yb
  else if (atomicnr == 71) {enth = 1031.20;}     //Lu
  else if (atomicnr == 72) {enth = 148.00;}      //Hf
  else if (atomicnr == 73) {enth = 186.90;}      //Ta
  else if (atomicnr == 74) {enth = 203.10;}      //W
  else if (atomicnr == 75) {enth = 185.00;}      //Re
  else if (atomicnr == 76) {enth = 188.00;}      //Os
  else if (atomicnr == 77) {enth = 160.00;}      //Ir
  else if (atomicnr == 78) {enth = 135.20;}      //Pt
  else if (atomicnr == 79) {enth = 88.00;}       //Au
  else if (atomicnr == 80) {enth = 14.69;}       //Hg
  else if (atomicnr == 81) {enth = 43.55;}       //Tl
  else if (atomicnr == 82) {enth = 46.62;}       //Pb
  else if (atomicnr == 83) {enth = 50.10;}       //Bi
  else if (atomicnr == 86) {enth = 0.000;}       //Rn
  else if (atomicnr == 102) {enth = 207.0;}      //No
  else {
    std::ostringstream ostr;
    ostr << atomicnr;
    throw("AtomPackage.hpp - ElementFormExp: Heat of formation parameter not defined for atomic number "+ostr.str());
  }
  return enth/au2kcalmol;
}
size_t AtomRow(size_t atom) {
  //function returning the row in the periodic table of a certain atom
  size_t row = 0;
  if (atom < 3) {row = 1;}
  else if ((atom > 2)&&(atom < 11)) {row = 2;}
  else if ((atom > 10)&&(atom < 19)) {row = 3;}
  else if ((atom > 18)&&(atom < 37)) {row = 4;}
  else if ((atom > 36)&&(atom < 55)) {row = 5;}
  else if ((atom > 54)&&(atom < 87)) {row = 6;}
  else if (atom > 86) {row = 7;}
  return row;
}
bool MetalHead(size_t atmnr) {
  bool metal = false;
  switch (atmnr) {
    case 1:
      metal = false;
      break;
    case 2:
      metal = false;
      break;
    case 3:
      metal = true;
      break;
    case 4:
      metal = true;
      break;
    case 5:
      metal = true;
      break;
    case 6:
      metal = false;
      break;
    case 7:
      metal = false;
      break;
    case 8:
      metal = false;
      break;
    case 9:
      metal = false;
      break;
    case 10:
      metal = false;
      break;
    case 11:
      metal = true;
      break;
    case 12:
      metal = true;
      break;
    case 13:
      metal = true;
      break;
    case 14:
      metal = false;
      break;
    case 15:
      metal = false;
      break;
    case 16:
      metal = false;
      break;
    case 17:
      metal = false;
      break;
    case 18:
      metal = false;
      break;
    case 19:
      metal = true;
      break;
    case 20:
      metal = true;
      break;
    case 21:
      metal = true;
      break;
    case 22:
      metal = true;
      break;
    case 23:
      metal = true;
      break;
    case 24:
      metal = true;
      break;
    case 25:
      metal = true;
      break;
    case 26:
      metal = true;
      break;
    case 27:
      metal = true;
      break;
    case 28:
      metal = true;
      break;
    case 29:
      metal = true;
      break;
    case 30:
      metal = true;
      break;
    case 31:
      metal = true;
      break;
    case 32:
      metal = false;
      break;
    case 33:
      metal = false;
      break;
    case 34:
      metal = false;
      break;
    case 35:
      metal = false;
      break;
    case 36:
      metal = false;
      break;
    case 37:
      metal = true;
      break;
    case 38:
      metal = true;
      break;
    case 39:
      metal = true;
      break;
    case 40:
      metal = true;
      break;
    case 41:
      metal = true;
      break;
    case 42:
      metal = true;
      break;
    case 43:
      metal = true;
      break;
    case 44:
      metal = true;
      break;
    case 45:
      metal = true;
      break;
    case 46:
      metal = true;
      break;
    case 47:
      metal = true;
      break;
    case 48:
      metal = true;
      break;
    case 49:
      metal = true;
      break;
    case 50:
      metal = false;
      break;
    case 51:
      metal = false;
      break;
    case 52:
      metal = false;
      break;
    case 53:
      metal = false;
      break;
    case 54:
      metal = false;
      break;
    case 55:
      metal = true;
      break;
    case 56:
      metal = true;
      break;
    case 57:
      metal = true;
      break;
    case 58:
      metal = true;
      break;
    case 59:
      metal = true;
      break;
    case 60:
      metal = true;
      break;
    case 61:
      metal = true;
      break;
    case 62:
      metal = true;
      break;
    case 63:
      metal = true;
      break;
    case 64:
      metal = true;
      break;
    case 65:
      metal = true;
      break;
    case 66:
      metal = true;
      break;
    case 67:
      metal = true;
      break;
    case 68:
      metal = true;
      break;
    case 69:
      metal = true;
      break;
    case 70:
      metal = true;
      break;
    case 71:
      metal = true;
      break;
    case 72:
      metal = true;
      break;
    case 73:
      metal = true;
      break;
    case 74:
      metal = true;
      break;
    case 75:
      metal = true;
      break;
    case 76:
      metal = true;
      break;
    case 77:
      metal = true;
      break;
    case 78:
      metal = true;
      break;
    case 79:
      metal = true;
      break;
    case 80:
      metal = true;
      break;
    case 81:
      metal = true;
      break;
    case 82:
      metal = false;
      break;
    case 83:
      metal = false;
      break;
    case 84:
      metal = false;
      break;
    case 85:
      metal = false;
      break;
    case 86:
      metal = false;
      break;
  }
  return metal;
}


#endif //_Atomic_Weight_Package_
