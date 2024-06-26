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

#ifndef _Hbond4_Parameters_
#define _Hbond4_Parameters_

//description:
//parameters for H4 and X correction terms

double CovalentRadius(size_t atomicnr) {
  //function return the covalent radii of atoms
  double R0 = 0.0;
  if (atomicnr == 1) {R0 = 0.37;}             //H
  else if (atomicnr == 2) {R0 = 0.32;}        //He
  else if (atomicnr == 3) {R0 = 1.34;}        //Li
  else if (atomicnr == 4) {R0 = 0.90;}        //Be
  else if (atomicnr == 5) {R0 = 0.82;}        //B
  else if (atomicnr == 6) {R0 = 0.77;}        //C
  else if (atomicnr == 7) {R0 = 0.75;}        //N
  else if (atomicnr == 8) {R0 = 0.73;}        //O
  else if (atomicnr == 9) {R0 = 0.71;}        //F
  else if (atomicnr == 10) {R0 = 0.69;}       //Ne
  else if (atomicnr == 11) {R0 = 1.54;}       //Na
  else if (atomicnr == 12) {R0 = 1.30;}       //Mg
  else if (atomicnr == 13) {R0 = 1.18;}       //Al
  else if (atomicnr == 14) {R0 = 1.11;}       //Si
  else if (atomicnr == 15) {R0 = 1.06;}       //P
  else if (atomicnr == 16) {R0 = 1.02;}       //S
  else if (atomicnr == 17) {R0 = 0.99;}       //Cl
  else if (atomicnr == 18) {R0 = 0.97;}       //Ar
  else if (atomicnr == 19) {R0 = 1.96;}       //K
  else if (atomicnr == 20) {R0 = 1.74;}       //Ca
  else if (atomicnr == 21) {R0 = 1.44;}       //Sc
  else if (atomicnr == 22) {R0 = 1.36;}       //Ti
  else if (atomicnr == 23) {R0 = 1.25;}       //V
  else if (atomicnr == 24) {R0 = 1.27;}       //Cr
  else if (atomicnr == 25) {R0 = 1.39;}       //Mn
  else if (atomicnr == 26) {R0 = 1.25;}       //Fe
  else if (atomicnr == 27) {R0 = 1.26;}       //Co
  else if (atomicnr == 28) {R0 = 1.21;}       //Ni
  else if (atomicnr == 29) {R0 = 1.38;}       //Cu
  else if (atomicnr == 30) {R0 = 1.31;}       //Zn
  else if (atomicnr == 31) {R0 = 1.26;}       //Ga
  else if (atomicnr == 32) {R0 = 1.22;}       //Ge
  else if (atomicnr == 33) {R0 = 1.19;}       //As
  else if (atomicnr == 34) {R0 = 1.16;}       //Se
  else if (atomicnr == 35) {R0 = 1.14;}       //Br
  else if (atomicnr == 36) {R0 = 1.10;}       //Kr
  else if (atomicnr == 37) {R0 = 2.11;}       //Rb
  else if (atomicnr == 38) {R0 = 1.92;}       //Sr
  else if (atomicnr == 39) {R0 = 1.62;}       //Y
  else if (atomicnr == 40) {R0 = 1.48;}       //Zr
  else if (atomicnr == 41) {R0 = 1.37;}       //Nb
  else if (atomicnr == 42) {R0 = 1.45;}       //Mo
  else if (atomicnr == 43) {R0 = 1.56;}       //Tc
  else if (atomicnr == 44) {R0 = 1.26;}       //Ru
  else if (atomicnr == 45) {R0 = 1.35;}       //Rh
  else if (atomicnr == 46) {R0 = 1.31;}       //Pd
  else if (atomicnr == 47) {R0 = 1.53;}       //Ag
  else if (atomicnr == 48) {R0 = 1.48;}       //Cd
  else if (atomicnr == 49) {R0 = 1.44;}       //In
  else if (atomicnr == 50) {R0 = 1.41;}       //Sn
  else if (atomicnr == 51) {R0 = 1.38;}       //Sb
  else if (atomicnr == 52) {R0 = 1.35;}       //Te
  else if (atomicnr == 53) {R0 = 1.33;}       //I
  else if (atomicnr == 54) {R0 = 1.30;}       //Xe
  else if (atomicnr == 55) {R0 = 2.25;}       //Cs
  else if (atomicnr == 56) {R0 = 1.98;}       //Ba
  else if (atomicnr == 57) {R0 = 1.69;}       //La
  else if (atomicnr == 71) {R0 = 1.60;}       //Lu
  else if (atomicnr == 72) {R0 = 1.50;}       //Hf
  else if (atomicnr == 73) {R0 = 1.38;}       //Ta
  else if (atomicnr == 74) {R0 = 1.46;}       //W
  else if (atomicnr == 75) {R0 = 1.59;}       //Re
  else if (atomicnr == 76) {R0 = 1.28;}       //Os
  else if (atomicnr == 77) {R0 = 1.37;}       //Ir
  else if (atomicnr == 78) {R0 = 1.28;}       //Pt
  else if (atomicnr == 79) {R0 = 1.44;}       //Au
  else if (atomicnr == 80) {R0 = 1.49;}       //Hg
  else if (atomicnr == 83) {R0 = 1.46;}       //Bi
  else if (atomicnr == 86) {R0 = 1.45;}       //Rn
  return R0;
}
double Cconst(size_t donor, size_t acceptor, const std::string & method) {
  //function yielding the c constant in from of the H4 correction term
  //donor and acceptor are atomic numbers
  //values stored in kcal/mol but return in a.u.
  double cconst = 0.0;
  if ((method == "PM6")||(method == "pm6")) {
    if ((donor == 8)&&(acceptor == 8)) {cconst = 2.32;}              //O-H......O
    else if ((donor == 8)&&(acceptor == 7)) {cconst = 3.10;}         //O-H......N
    else if ((donor == 7)&&(acceptor == 8)) {cconst = 1.07;}         //N-H......O
    else if ((donor == 7)&&(acceptor == 7)) {cconst = 2.01;}         //N-H......N
  }
  else if ((method == "SCC-DFTB")||(method == "scc-dftb")) {
    if ((donor == 8)&&(acceptor == 8)) {cconst = 1.11;}              //O-H......O
    else if ((donor == 8)&&(acceptor == 7)) {cconst = 2.58;}         //O-H......N
    else if ((donor == 7)&&(acceptor == 8)) {cconst = 0.80;}         //N-H......O
    else if ((donor == 7)&&(acceptor == 7)) {cconst = 2.01;}         //N-H......N
  }
  else if ((method == "RM1")||(method == "rm1")) {
    if ((donor == 8)&&(acceptor == 8)) {cconst = 3.76;}              //O-H......O
    else if ((donor == 8)&&(acceptor == 7)) {cconst = 3.90;}         //O-H......N
    else if ((donor == 7)&&(acceptor == 8)) {cconst = 3.14;}         //N-H......O
    else if ((donor == 7)&&(acceptor == 7)) {cconst = 2.95;}         //N-H......N
  }
  else if ((method == "AM1")||(method == "am1")) {
    if ((donor == 8)&&(acceptor == 8)) {cconst = 4.89;}              //O-H......O
    else if ((donor == 8)&&(acceptor == 7)) {cconst = 6.23;}         //O-H......N
    else if ((donor == 7)&&(acceptor == 8)) {cconst = 2.54;}         //N-H......O
    else if ((donor == 7)&&(acceptor == 7)) {cconst = 4.56;}         //N-H......N
  }
  else if ((method == "PM3")||(method == "pm3")) {
    if ((donor == 8)&&(acceptor == 8)) {cconst = 2.71;}              //O-H......O
    else if ((donor == 8)&&(acceptor == 7)) {cconst = 4.37;}         //O-H......N
    else if ((donor == 7)&&(acceptor == 8)) {cconst = 2.29;}         //N-H......O
    else if ((donor == 7)&&(acceptor == 7)) {cconst = 3.86;}         //N-H......N
  }
  return cconst/au2kcalmol;
}
double CH2O(const std::string & method) {
  //function yielding the c constant for water as donor
  double ch2o = 0.0;
  if ((method == "PM6")||(method == "pm6")) {ch2o = 0.42;}
  else if ((method == "SCC-DFTB")||(method == "scc-dftb")) {ch2o = 1.32;}
  else if ((method == "RM1")||(method == "rm1")) {ch2o = 0.94;}
  else if ((method == "AM1")||(method == "am1")) {ch2o = 0.49;}
  else if ((method == "PM3")||(method == "pm3")) {ch2o = 0.91;}
  return ch2o;
}
double CS(std::string group, const std::string & method) {
  //function returning the value for the cs constant for some types of charged ions
  double cs = 0.0;
  if ((method == "PM6")||(method == "pm6")) {
    if ((group == "RCOO")||(group == "rcoo")||(group == "RCOO-")||(group == "rcoo-")) {cs = 1.41;}
    else if ((group == "R3NH")||(group == "r3nh")||(group == "R3NH+")||(group == "r3nh+")) {cs = 3.61;}
    else if ((group == "GUA")||(group == "gua")) {cs = 1.26;}
    else if ((group == "IMID")||(group == "imid")) {cs = 2.29;}
  }
  else if ((method == "SCC-DFTB")||(method == "scc-dftb")) {
    if ((group == "RCOO")||(group == "rcoo")||(group == "RCOO-")||(group == "rcoo-")) {cs = 1.22;}
    else if ((group == "R3NH")||(group == "r3nh")||(group == "R3NH+")||(group == "r3nh+")) {cs = 2.33;}
    else if ((group == "GUA")||(group == "gua")) {cs = 2.42;}
    else if ((group == "IMID")||(group == "imid")) {cs = 3.44;}
  }
  else if ((method == "RM1")||(method == "rm1")) {
    if ((group == "RCOO")||(group == "rcoo")||(group == "RCOO-")||(group == "rcoo-")) {cs = 1.10;}
    else if ((group == "R3NH")||(group == "r3nh")||(group == "R3NH+")||(group == "r3nh+")) {cs = 1.21;}
    else if ((group == "GUA")||(group == "gua")) {cs = 1.18;}
    else if ((group == "IMID")||(group == "imid")) {cs = 1.10;}
  }
  else if ((method == "AM1")||(method == "am1")) {
    if ((group == "RCOO")||(group == "rcoo")||(group == "RCOO-")||(group == "rcoo-")) {cs = 1.08;}
    else if ((group == "R3NH")||(group == "r3nh")||(group == "R3NH+")||(group == "r3nh+")) {cs = 2.78;}
    else if ((group == "GUA")||(group == "gua")) {cs = 0.86;}
    else if ((group == "IMID")||(group == "imid")) {cs = 2.11;}
  }
  else if ((method == "PM3")||(method == "pm3")) {
    if ((group == "RCOO")||(group == "rcoo")||(group == "RCOO-")||(group == "rcoo-")) {cs = 0.89;}
    else if ((group == "R3NH")||(group == "r3nh")||(group == "R3NH+")||(group == "r3nh+")) {cs = 2.54;}
    else if ((group == "GUA")||(group == "gua")) {cs = 1.54;}
    else if ((group == "IMID")||(group == "imid")) {cs = 1.84;}
  }
  return cs;
}
double sHH(const std::string & method) {
  //function returning the scaling parameter sHH
  //values stored in kcal/mol but returned in a.u.
  double shh = 0.0;
  if ((method == "PM6")||(method == "pm6")) {shh = 0.40;}
  else if ((method == "SCC-DFTB")||(method == "scc-dftb")) {shh = 0.30;}
  else if ((method == "RM1")||(method == "rm1")) {shh = 0.30;}
  else if ((method == "AM1")||(method == "am1")) {shh = 0.90;}
  else if ((method == "PM3")||(method == "pm3")) {shh = 0.90;}
  return shh/au2kcalmol;
}
double eHH(const std::string & method) {
  //function returning the exponential factor eHH
  double ehh = 0.0;
  if ((method == "PM6")||(method == "pm6")) {ehh = 12.70;}
  else if ((method == "SCC-DFTB")||(method == "scc-dftb")) {ehh = 14.31;}
  else if ((method == "RM1")||(method == "rm1")) {ehh = 4.46;}
  else if ((method == "AM1")||(method == "am1")) {ehh = 4.46;}
  else if ((method == "PM3")||(method == "pm3")) {ehh = 6.86;}
  return ehh;
}
double R0HH(const std::string & method) {
  //function returning the van der Waals distance between two hydrogen atoms R0HH
  //values stored and returned in Angstrom
  double r0hh = 0.0;
  if ((method == "PM6")||(method == "pm6")) {r0hh = 2.30;}
  else if ((method == "SCC-DFTB")||(method == "scc-dftb")) {r0hh = 2.35;}
  else if ((method == "RM1")||(method == "rm1")) {r0hh = 2.11;}
  else if ((method == "AM1")||(method == "am1")) {r0hh = 2.11;}
  else if ((method == "PM3")||(method == "pm3")) {r0hh = 2.23;}
  return r0hh;
}
//parameters for the halogen bond correction
double aparameter(size_t donor, size_t acceptor, const std::string & method) {
  //function returning the scaling parameter a
  //values stored in kcal/mol but returned in a.u.
  double aa = 0.0;
  if ((method == "PM6")||(method == "pm6")) {
    if (acceptor == 8) {
      if (donor == 17) {aa = 4.6783e8;}
      else if (donor == 35) {aa = 9.6021e3;}
      else if (donor == 53) {aa = 6.0912e5;}
    }
    else if (acceptor == 7) {
      if (donor == 17) {aa = 1.0489e12;}
      else if (donor == 35) {aa = 1.0226e5;}
      else if (donor == 53) {aa = 1.2751e12;}
    }
  }
  return aa/au2kcalmol;
}
double bparameter(size_t donor, size_t acceptor, const std::string & method) {
  //function returning the exponential parameter b
  //values stored in 1/Angstrom
  double bb = 0.0;
  if ((method == "PM6")||(method == "pm6")) {
    if (acceptor == 8) {
      if (donor == 17) {bb = 6.867;}
      else if (donor == 35) {bb = 2.900;}
      else if (donor == 53) {bb = 4.154;}
    }
    else if (acceptor == 7) {
      if (donor == 17) {bb = 9.946;}
      else if (donor == 35) {bb = 3.236;}
      else if (donor == 53) {bb = 9.534;}
    }
  }
  return bb;
}

#endif //_Hbond4_Parameters_
