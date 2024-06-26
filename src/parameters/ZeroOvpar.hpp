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

#ifndef _Zero_Overlap_Parameters_
#define _Zero_Overlap_Parameters_

//description:
//functions returning the interatomic distance (in Angstrom) at which the overlap between two atoms of type atm is below 1.0e-10, for different methods

double ZeroOverlapMNDOd(size_t atm) {
  double zeroov = 10000000000000000.0;
  switch (atm) {
    case 1:      //H
      zeroov = 11.6;
      break;
    case 3:      //Li
      zeroov = 26.2;
      break;
    case 4:      //Be
      zeroov = 18.4;
      break;
    case 5:      //B
      zeroov = 12.3;
      break;
    case 6:      //C
      zeroov = 10.4;
      break;
    case 7:      //N
      zeroov = 8.3;
      break;
    case 8:      //O
      zeroov = 6.9;
      break;
    case 9:      //F
      zeroov = 6.6;
      break;
    case 11:     //Na
      zeroov = 23.2;
      break;
    case 12:     //Mg
      zeroov = 21.7;
      break;
    case 13:     //Al
      zeroov = 22.0;
      break;
    case 14:     //Si
      zeroov = 21.6;
      break;
    case 15:     //P
      zeroov = 19.0;
      break;
    case 16:     //S
      zeroov = 17.0;
      break;
    case 17:     //Cl
      zeroov = 16.7;
      break;
    case 30:     //Zn
      zeroov = 16.4;
      break;
    case 32:     //Ge
      zeroov = 17.2;
      break;
    case 35:     //Br
      zeroov = 17.0;
      break;
    case 48:     //Cd
      zeroov = 15.9;
      break;
    case 50:     //Sn
      zeroov = 12.9;
      break;
    case 53:     //I
      zeroov = 16.8;
      break;
    case 80:     //Hg
      zeroov = 15.7;
      break;
    case 82:     //Pb
      zeroov = 12.9;
      break;
  }
  return zeroov;
}
double ZeroOverlapPM6(size_t atm) {
  double zeroov = 10000000000000000.0;
  switch (atm) {
    case 1:      //H
      zeroov = 12.1;
      break;
    case 2:      //He
      zeroov = 5.1;
      break;
    case 3:      //Li
      zeroov = 18.2;
      break;
    case 4:      //Be
      zeroov = 14.7;
      break;
    case 5:      //B
      zeroov = 12.5;
      break;
    case 6:      //C
      zeroov = 10.9;
      break;
    case 7:      //N
      zeroov = 9.3;
      break;
    case 8:      //O
      zeroov = 8.2;
      break;
    case 9:      //F
      zeroov = 6.4;
      break;
    case 10:     //Ne
      zeroov = 4.9;
      break;
    case 11:     //Na
      zeroov = 29.2;
      break;
    case 12:     //Mg
      zeroov = 15.3;
      break;
    case 13:     //Al
      zeroov = 16.5;
      break;
    case 14:     //Si
      zeroov = 17.3;
      break;
    case 15:     //P
      zeroov = 17.0;
      break;
    case 16:     //S
      zeroov = 11.3;
      break;
    case 17:     //Cl
      zeroov = 15.8;
      break;
    case 18:     //Ar
      zeroov = 3.8;
      break;
    case 19:     //K
      zeroov = 20.3;
      break;
    case 20:     //Ca
      zeroov = 14.5;
      break;
    case 21:     //Sc
      zeroov = 17.0;
      break;
    case 22:     //Ti
      zeroov = 19.6;
      break;
    case 23:     //V
      zeroov = 21.5;
      break;
    case 24:     //Cr
      zeroov = 22.2;
      break;
    case 25:     //Mn
      zeroov = 15.0;
      break;
    case 26:     //Fe
      zeroov = 19.4;
      break;
    case 27:     //Co
      zeroov = 19.0;
      break;
    case 28:     //Ni
      zeroov = 14.0;
      break;
    case 29:     //Cu
      zeroov = 13.3;
      break;
    case 30:     //Zn
      zeroov = 14.7;
      break;
    case 31:     //Ga
      zeroov = 13.3;
      break;
    case 32:     //Ge
      zeroov = 13.4;
      break;
    case 33:     //As
      zeroov = 16.6;
      break;
    case 34:     //Se
      zeroov = 11.4;
      break;
    case 35:     //Br
      zeroov = 15.2;
      break;
    case 36:     //Kr
      zeroov = 18.4;
      break;
    case 37:     //Rb
      zeroov = 18.6;
      break;
    case 38:     //Sr
      zeroov = 14.4;
      break;
    case 39:     //Y
      zeroov = 40.6;
      break;
    case 40:     //Zr
      zeroov = 14.8;
      break;
    case 41:     //Nb
      zeroov = 17.9;
      break;
    case 42:     //Mo
      zeroov = 22.8;
      break;
    case 43:     //Tc
      zeroov = 13.1;
      break;
    case 44:     //Ru
      zeroov = 16.6;
      break;
    case 45:     //Rh
      zeroov = 18.3;
      break;
    case 46:     //Pd
      zeroov = 21.5;
      break;
    case 47:     //Ag
      zeroov = 36.4;
      break;
    case 48:     //Cd
      zeroov = 17.5;
      break;
    case 49:     //In
      zeroov = 12.0;
      break;
    case 50:     //Sn
      zeroov = 12.1;
      break;
    case 51:     //Sb
      zeroov = 14.1;
      break;
    case 52:     //Te
      zeroov = 14.4;
      break;
    case 53:     //I
      zeroov = 13.5;
      break;
    case 54:     //Xe
      zeroov = 12.6;
      break;
    case 55:     //Cs
      zeroov = 16.6;
      break;
    case 56:     //Ba
      zeroov = 18.8;
      break;
    case 57:     //La
      zeroov = 21.5;
      break;
    case 71:     //Lu
      zeroov = 15.7;
      break;
    case 72:     //Hf
      zeroov = 17.0;
      break;
    case 73:     //Ta
      zeroov = 13.7;
      break;
    case 74:     //W
      zeroov = 16.5;
      break;
    case 75:     //Re
      zeroov = 14.8;
      break;
    case 76:     //Os
      zeroov = 16.9;
      break;
    case 77:     //Ir
      zeroov = 17.4;
      break;
    case 78:     //Pt
      zeroov = 16.2;
      break;
    case 79:     //Au
      zeroov = 16.6;
      break;
    case 80:     //Hg
      zeroov = 17.7;
      break;
    case 81:     //Tl
      zeroov = 15.2;
      break;
    case 82:     //Pb
      zeroov = 15.9;
      break;
    case 83:     //Bi
      zeroov = 14.4;
      break;
  }
  return zeroov;
}
double ZeroOverlapMNDO(size_t atm) {
  double zeroov = 10000000000000000.0;
  switch (atm) {
    case 1:      //H
      zeroov = 11.6;
      break;
    case 3:      //Li
      zeroov = 26.2;
      break;
    case 4:      //Be
      zeroov = 18.4;
      break;
    case 5:      //B
      zeroov = 12.3;
      break;
    case 6:      //C
      zeroov = 10.4;
      break;
    case 7:      //N
      zeroov = 8.3;
      break;
    case 8:      //O
      zeroov = 6.9;
      break;
    case 9:      //F
      zeroov = 6.6;
      break;
    case 13:     //Al
      zeroov = 14.4;
      break;
    case 14:     //Si
      zeroov = 15.3;
      break;
    case 15:     //P
      zeroov = 11.7;
      break;
    case 16:     //S
      zeroov = 10.4;
      break;
    case 17:     //Cl
      zeroov = 10.2;
      break;
    case 30:     //Zn
      zeroov = 15.7;
      break;
    case 32:     //Ge
      zeroov = 17.2;
      break;
    case 35:     //Br
      zeroov = 10.5;
      break;
    case 50:     //Sn
      zeroov = 12.9;
      break;
    case 53:     //I
      zeroov = 11.5;
      break;
    case 80:     //Hg
      zeroov = 13.0;
      break;
    case 82:     //Pb
      zeroov = 12.9;
      break;
  }
  return zeroov;
}
double ZeroOverlapAM1(size_t atm) {
  double zeroov = 10000000000000000.0;
  switch (atm) {
    case 1:      //H
      zeroov = 13.0;
      break;
    case 3:      //Li
      zeroov = 26.2;
      break;
    case 4:      //Be
      zeroov = 18.4;
      break;
    case 5:      //B
      zeroov = 11.9;
      break;
    case 6:      //C
      zeroov = 11.0;
      break;
    case 7:      //N
      zeroov = 8.6;
      break;
    case 8:      //O
      zeroov = 7.4;
      break;
    case 9:      //F
      zeroov = 7.5;
      break;
    case 11:     //Na
      zeroov = 29.4;
      break;
    case 13:     //Al
      zeroov = 15.9;
      break;
    case 14:     //Si
      zeroov = 16.1;
      break;
    case 15:     //P
      zeroov = 11.1;
      break;
    case 16:     //S
      zeroov = 12.5;
      break;
    case 17:     //Cl
      zeroov = 10.0;
      break;
    case 30:     //Zn
      zeroov = 16.7;
      break;
    case 32:     //Ge
      zeroov = 18.2;
      break;
    case 33:     //As
      zeroov = 13.3;
      break;
    case 34:     //Se
      zeroov = 11.2;
      break;
    case 35:     //Br
      zeroov = 11.3;
      break;
    case 51:     //Sb
      zeroov = 11.3;
      break;
    case 52:     //Te
      zeroov = 12.7;
      break;
    case 53:     //I
      zeroov = 11.6;
      break;
    case 80:     //Hg
      zeroov = 13.8;
      break;
  }
  return zeroov;
}
double ZeroOverlapPM3(size_t atm) {
  double zeroov = 10000000000000000.0;
  switch (atm) {
    case 1:      //H
      zeroov = 15.9;
      break;
    case 3:      //Li
      zeroov = 27.4;
      break;
    case 4:      //Be
      zeroov = 20.3;
      break;
    case 6:      //C
      zeroov = 11.4;
      break;
    case 7:      //N
      zeroov = 8.9;
      break;
    case 8:      //O
      zeroov = 7.8;
      break;
    case 9:      //F
      zeroov = 7.5;
      break;
    case 11:     //Na
      zeroov = 17.7;
      break;
    case 12:     //Mg
      zeroov = 28.7;
      break;
    case 13:     //Al
      zeroov = 19.3;
      break;
    case 14:     //Si
      zeroov = 15.8;
      break;
    case 15:     //P
      zeroov = 13.8;
      break;
    case 16:     //S
      zeroov = 12.5;
      break;
    case 17:     //Cl
      zeroov = 9.7;
      break;
    case 30:     //Zn
      zeroov = 15.2;
      break;
    case 31:     //Ga
      zeroov = 27.2;
      break;
    case 32:     //Ge
      zeroov = 14.4;
      break;
    case 33:     //As
      zeroov = 13.5;
      break;
    case 34:     //Se
      zeroov = 13.2;
      break;
    case 35:     //Br
      zeroov = 10.8;
      break;
    case 48:     //Cd
      zeroov = 14.4;
      break;
    case 49:     //In
      zeroov = 17.2;
      break;
    case 50:     //Sn
      zeroov = 15.2;
      break;
    case 51:     //Sb
      zeroov = 13.1;
      break;
    case 52:     //Te
      zeroov = 15.1;
      break;
    case 53:     //I
      zeroov = 10.2;
      break;
    case 80:     //Hg
      zeroov = 17.7;
      break;
    case 81:     //Tl
      zeroov = 13.7;
      break;
    case 82:     //Pb
      zeroov = 14.2;
      break;
    case 83:     //Bi
      zeroov = 13.9;
      break;
  }
  return zeroov;
}
double ZeroOverlapPM3PDDG(size_t atm) {
  double zeroov = 10000000000000000.0;
  switch (atm) {
    case 1:      //H
      zeroov = 15.8;
      break;
    case 6:      //C
      zeroov = 11.4;
      break;
    case 7:      //N
      zeroov = 8.8;
      break;
    case 8:      //O
      zeroov = 8.0;
      break;
    case 9:      //F
      zeroov = 7.4;
      break;
    case 14:     //Si
      zeroov = 14.0;
      break;
    case 15:     //P
      zeroov = 11.9;
      break;
    case 16:     //S
      zeroov = 12.7;
      break;
    case 17:     //Cl
      zeroov = 9.1;
      break;
    case 35:     //Br
      zeroov = 10.5;
      break;
    case 53:     //I
      zeroov = 10.4;
      break;
  }
  return zeroov;
}
double ZeroOverlapMNDOPDDG(size_t atm) {
  double zeroov = 10000000000000000.0;
  switch (atm) {
    case 1:      //H
      zeroov = 11.7;
      break;
    case 6:      //C
      zeroov = 10.2;
      break;
    case 7:      //N
      zeroov = 8.3;
      break;
    case 8:      //O
      zeroov = 7.0;
      break;
    case 9:      //F
      zeroov = 6.4;
      break;
    case 17:     //Cl
      zeroov = 10.2;
      break;
    case 35:     //Br
      zeroov = 10.2;
      break;
    case 53:     //I
      zeroov = 10.2;
      break;
  }
  return zeroov;
}
double ZeroOverlapRM1(size_t atm) {
  double zeroov = 10000000000000000.0;
  switch (atm) {
    case 1:      //H
      zeroov = 14.2;
      break;
    case 6:      //C
      zeroov = 10.5;
      break;
    case 7:      //N
      zeroov = 9.4;
      break;
    case 8:      //O
      zeroov = 7.3;
      break;
    case 9:      //F
      zeroov = 7.1;
      break;
    case 15:     //P
      zeroov = 11.9;
      break;
    case 16:     //S
      zeroov = 11.1;
      break;
    case 17:     //Cl
      zeroov = 11.0;
      break;
    case 35:     //Br
      zeroov = 11.3;
      break;
    case 53:     //I
      zeroov = 10.8;
      break;
  }
  return zeroov;
}
double ZeroOverlapGFN2(size_t atm) {
  double zeroov = 10000000000000000.0;
  switch (atm) {
    case 1:      //H
      zeroov = 8.6;
      break;
    case 2:      //He
      zeroov = 9.7;
      break;
    case 3:      //Li
      zeroov = 23.8;
      break;
    case 4:      //Be
      zeroov = 15.3;
      break;
    case 5:      //B
      zeroov = 9.9;
      break;
    case 6:      //C
      zeroov = 8.1;
      break;
    case 7:      //N
      zeroov = 7.3;
      break;
    case 8:      //O
      zeroov = 6.9;
      break;
    case 9:      //F
      zeroov = 6.4;
      break;
    case 10:     //Ne
      zeroov = 6.4;
      break;
    case 11:     //Na
      zeroov = 23.8;
      break;
    case 12:     //Mg
      zeroov = 23.8;
      break;
    case 13:     //Al
      zeroov = 15.0;
      break;
    case 14:     //Si
      zeroov = 12.0;
      break;
    case 15:     //P
      zeroov = 12.8;
      break;
    case 16:     //S
      zeroov = 9.3;
      break;
    case 17:     //Cl
      zeroov = 8.3;
      break;
    case 18:     //Ar
      zeroov = 8.5;
      break;
    case 19:     //K
      zeroov = 23.8;
      break;
    case 20:     //Ca
      zeroov = 23.8;
      break;
    case 21:     //Sc
      zeroov = 19.2;
      break;
    case 22:     //Ti
      zeroov = 20.5;
      break;
    case 23:     //V
      zeroov = 20.9;
      break;
    case 24:     //Cr
      zeroov = 18.1;
      break;
    case 25:     //Mn
      zeroov = 16.8;
      break;
    case 26:     //Fe
      zeroov = 20.1;
      break;
    case 27:     //Co
      zeroov = 15.1;
      break;
    case 28:     //Ni
      zeroov = 14.9;
      break;
    case 29:     //Cu
      zeroov = 13.3;
      break;
    case 30:     //Zn
      zeroov = 16.7;
      break;
    case 31:     //Ga
      zeroov = 18.1;
      break;
    case 32:     //Ge
      zeroov = 17.3;
      break;
    case 33:     //As
      zeroov = 18.3;
      break;
    case 34:     //Se
      zeroov = 14.4;
      break;
    case 35:     //Br
      zeroov = 10.4;
      break;
    case 36:     //Kr
      zeroov = 10.1;
      break;
    case 37:     //Rb
      zeroov = 23.8;
      break;
    case 38:     //Sr
      zeroov = 23.8;
      break;
    case 39:     //Y
      zeroov = 19.4;
      break;
    case 40:     //Zr
      zeroov = 20.0;
      break;
    case 41:     //Nb
      zeroov = 20.0;
      break;
    case 42:     //Mo
      zeroov = 19.5;
      break;
    case 43:     //Tc
      zeroov = 16.8;
      break;
    case 44:     //Ru
      zeroov = 16.8;
      break;
    case 45:     //Rh
      zeroov = 16.2;
      break;
    case 46:     //Pd
      zeroov = 17.0;
      break;
    case 47:     //Ag
      zeroov = 17.9;
      break;
    case 48:     //Cd
      zeroov = 19.8;
      break;
    case 49:     //In
      zeroov = 18.9;
      break;
    case 50:     //Sn
      zeroov = 18.1;
      break;
    case 51:     //Sb
      zeroov = 15.9;
      break;
    case 52:     //Te
      zeroov = 14.5;
      break;
    case 53:     //I
      zeroov = 11.9;
      break;
    case 54:     //Xe
      zeroov = 10.8;
      break;
    case 55:     //Cs
      zeroov = 23.8;
      break;
    case 56:     //Ba
      zeroov = 23.8;
      break;
    case 57:     //La
      zeroov = 19.9;
      break;
    case 58:     //Ce
      zeroov = 19.8;
      break;
    case 59:     //Pr
      zeroov = 19.8;
      break;
    case 60:     //Nd
      zeroov = 19.7;
      break;
    case 61:     //Pm
      zeroov = 19.6;
      break;
    case 62:     //Sm
      zeroov = 19.5;
      break;
    case 63:     //Eu
      zeroov = 19.4;
      break;
    case 64:     //Gd
      zeroov = 19.3;
      break;
    case 65:     //Tb
      zeroov = 19.3;
      break;
    case 66:     //Dy
      zeroov = 19.2;
      break;
    case 67:     //Ho
      zeroov = 19.1;
      break;
    case 68:     //Er
      zeroov = 19.0;
      break;
    case 69:     //Tm
      zeroov = 19.0;
      break;
    case 70:     //Yb
      zeroov = 18.9;
      break;
    case 71:     //Lu
      zeroov = 18.8;
      break;
    case 72:     //Hf
      zeroov = 18.2;
      break;
    case 73:     //Ta
      zeroov = 18.5;
      break;
    case 74:     //W
      zeroov = 17.8;
      break;
    case 75:     //Re
      zeroov = 15.9;
      break;
    case 76:     //Os
      zeroov = 16.3;
      break;
    case 77:     //Ir
      zeroov = 15.5;
      break;
    case 78:     //Pt
      zeroov = 16.0;
      break;
    case 79:     //Au
      zeroov = 12.5;
      break;
    case 80:     //Hg
      zeroov = 17.7;
      break;
    case 81:     //Tl
      zeroov = 15.1;
      break;
    case 82:     //Pb
      zeroov = 13.4;
      break;
    case 83:     //Bi
      zeroov = 11.5;
      break;
    case 84:     //Po
      zeroov = 10.9;
      break;
    case 85:     //At
      zeroov = 13.3;
      break;
    case 86:     //Rn
      zeroov = 11.2;
      break;
  }
  return zeroov;
}

#endif //_Zero_Overlap_Parameters_
