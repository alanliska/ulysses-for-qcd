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

#ifndef _Dewar_Basis_Set_
#define _Dewar_Basis_Set_
#include <vector>
#include <string>
#include "STO.hpp"
#include "../parameters/ZeroOvpar.hpp"

//description:
//exponents for Slater type orbitals as optimized for each specific NDDO semi-empirical method

class PM6BSet: public STOBSet {
public:
  PM6BSet() {}
  ~PM6BSet() {}
  std::vector<double> GetExponents(size_t atomnum) {
    //getter for the exponents of a given atom
    std::vector<double> exponents;
    if (atomnum == 1) {exponents.push_back(1.268641);}  //H
    else if (atomnum == 2) {                            //He
      exponents.push_back(3.313204);                        //1s
      exponents.push_back(3.657133);                        //2p
    }
    else if (atomnum == 3) {                            //Li
      exponents.push_back(0.981041);                        //2s
      exponents.push_back(2.953445);                        //2p
    }
    else if (atomnum == 4) {                            //Be
      exponents.push_back(1.212539);                        //2s
      exponents.push_back(1.276487);                        //2p
    }
    else if (atomnum == 5) {                            //B
      exponents.push_back(1.634174);                        //2s
      exponents.push_back(1.479195);                        //2p
    }
    else if (atomnum == 6) {                            //C
      exponents.push_back(2.047558);                        //2s
      exponents.push_back(1.702841);                        //2p
    }
    else if (atomnum == 7) {                            //N
      exponents.push_back(2.380406);                        //2s
      exponents.push_back(1.999246);                        //2p
    }
    else if (atomnum == 8) {                            //O
      exponents.push_back(5.421751);                        //2s
      exponents.push_back(2.270960);                        //2p
    }
    else if (atomnum == 9) {                            //F
      exponents.push_back(6.043849);                        //2s
      exponents.push_back(2.906722);                        //2p
    }
    else if (atomnum == 10) {                           //Ne
      exponents.push_back(6.000148);                        //2p
      exponents.push_back(3.834528);                        //2p
    }
    else if (atomnum == 11) {                           //Na
      exponents.push_back(0.686327);                        //3s
      exponents.push_back(0.950068);                        //3p
    }
    else if (atomnum == 12) {                           //Mg
      exponents.push_back(1.310830);                        //3s
      exponents.push_back(1.388897);                        //3p
    }
    else if (atomnum == 13) {                           //Al
      exponents.push_back(2.364264);                        //3s
      exponents.push_back(1.749102);                        //3p
      exponents.push_back(1.269384);                        //3d
    }
    else if (atomnum == 14) {                           //Si
      exponents.push_back(1.752741);                        //3s
      exponents.push_back(1.198413);                        //3p
      exponents.push_back(2.128593);                        //3d
    }
    else if (atomnum == 15) {                           //P
      exponents.push_back(2.158033);                        //3s
      exponents.push_back(1.805343);                        //3p
      exponents.push_back(1.230358);                        //3d
    }
    else if (atomnum == 16) {                           //S
      exponents.push_back(2.192844);                        //3s
      exponents.push_back(1.841078);                        //3p
      exponents.push_back(3.109401);                        //3d
    }
    else if (atomnum == 17) {                           //Cl
      exponents.push_back(2.637050);                        //3s
      exponents.push_back(2.118146);                        //3p
      exponents.push_back(1.324033);                        //3d
    }
    else if (atomnum == 18) {                           //Ar
      exponents.push_back(6.000272);                        //3s
      exponents.push_back(5.949170);                        //3p
    }
    else if (atomnum == 19) {                           //K
      exponents.push_back(6.000478);                        //4s
      exponents.push_back(1.127503);                        //4p
    }
    else if (atomnum == 20) {                           //Ca
      exponents.push_back(1.528258);                        //4s
      exponents.push_back(2.060094);                        //4p
    }
    else if (atomnum == 21) {                           //Sc
      exponents.push_back(1.859012);                        //3d
      exponents.push_back(1.402469);                        //4s
      exponents.push_back(1.345196);                        //4p
    }
    else if (atomnum == 22) {                           //Ti
      exponents.push_back(1.418280);                        //3d
      exponents.push_back(5.324777);                        //4s
      exponents.push_back(1.164068);                        //4p
    }
    else if (atomnum == 23) {                           //V
      exponents.push_back(1.394806);                        //3d
      exponents.push_back(1.974330);                        //4s
      exponents.push_back(1.063106);                        //4p
    }
    else if (atomnum == 24) {                           //Cr
      exponents.push_back(1.623119);                        //3d
      exponents.push_back(3.283460);                        //4s
      exponents.push_back(1.029394);                        //4p
    }
    else if (atomnum == 25) {                           //Mn
      exponents.push_back(2.607800);                        //3d
      exponents.push_back(2.131680);                        //4s
      exponents.push_back(1.525880);                        //4p
    }
    else if (atomnum == 26) {                           //Fe
      exponents.push_back(1.080747);                        //3d
      exponents.push_back(1.479150);                        //4s
      exponents.push_back(6.002246);                        //4p
    }
    else if (atomnum == 27) {                           //Co
      exponents.push_back(1.860218);                        //3d
      exponents.push_back(1.166613);                        //4s
      exponents.push_back(3.000000);                        //4p
    }
    else if (atomnum == 28) {                           //Ni
      exponents.push_back(2.514761);                        //3d
      exponents.push_back(1.591828);                        //4s
      exponents.push_back(2.304739);                        //4p
    }
    else if (atomnum == 29) {                           //Cu
      exponents.push_back(2.734990);                        //3d
      exponents.push_back(1.669096);                        //4s
      exponents.push_back(3.000000);                        //4p
    }
    else if (atomnum == 30) {                           //Zn
      exponents.push_back(1.512875);                        //4s
      exponents.push_back(1.789482);                        //4p
    }
    else if (atomnum == 31) {                           //Ga
      exponents.push_back(2.339067);                        //4s
      exponents.push_back(1.729592);                        //4p
    }
    else if (atomnum == 32) {                           //Ge
      exponents.push_back(2.546073);                        //4s
      exponents.push_back(1.709130);                        //4p
    }
    else if (atomnum == 33) {                           //As
      exponents.push_back(2.926171);                        //4s
      exponents.push_back(1.765191);                        //4p
      exponents.push_back(1.392142);                        //4d
    }
    else if (atomnum == 34) {                           //Se
      exponents.push_back(2.512366);                        //4s
      exponents.push_back(2.007576);                        //4p
    }
    else if (atomnum == 35) {                           //Br
      exponents.push_back(4.670684);                        //4s
      exponents.push_back(2.035626);                        //4p
      exponents.push_back(1.521031);                        //4d
    }
    else if (atomnum == 36) {                           //Kr
      exponents.push_back(1.312248);                        //4s
      exponents.push_back(4.491371);                        //4p
    }
    else if (atomnum == 37) {                           //Rb
      exponents.push_back(5.510145);                        //5s
      exponents.push_back(1.335170);                        //5p
    }
    else if (atomnum == 38) {                           //Sr
      exponents.push_back(2.197303);                        //5s
      exponents.push_back(1.730137);                        //5p
    }
    else if (atomnum == 39) {                           //Y
      exponents.push_back(1.650893);                        //4d
      exponents.push_back(0.593368);                        //5s
      exponents.push_back(1.490422);                        //5p
    }
    else if (atomnum == 40) {                           //Zr
      exponents.push_back(1.567392);                        //4d
      exponents.push_back(1.692590);                        //5s
      exponents.push_back(1.694916);                        //5p
    }
    else if (atomnum == 41) {                           //Nb
      exponents.push_back(1.977324);                        //4d
      exponents.push_back(2.355562);                        //5s
      exponents.push_back(1.386907);                        //5p
    }
    else if (atomnum == 42) {                           //Mo
      exponents.push_back(1.827152);                        //4d
      exponents.push_back(1.060429);                        //5s
      exponents.push_back(1.350412);                        //5p
    }
    else if (atomnum == 43) {                           //Tc
      exponents.push_back(1.767360);                        //4d
      exponents.push_back(1.956245);                        //5s
      exponents.push_back(6.006299);                        //5p
    }
    else if (atomnum == 44) {                           //Ru
      exponents.push_back(2.093164);                        //4d
      exponents.push_back(1.459195);                        //5s
      exponents.push_back(5.537201);                        //5p
    }
    else if (atomnum == 45) {                           //Rh
      exponents.push_back(2.901406);                        //4d
      exponents.push_back(1.324919);                        //5s
      exponents.push_back(4.306111);                        //5p
    }
    else if (atomnum == 46) {                           //Pd
      exponents.push_back(2.219861);                        //4d
      exponents.push_back(1.658503);                        //5s
      exponents.push_back(1.156718);                        //5p
    }
    else if (atomnum == 47) {                           //Ag
      exponents.push_back(6.007328);                        //4d
      exponents.push_back(1.994004);                        //5s
      exponents.push_back(0.681817);                        //5p
    }
    else if (atomnum == 48) {                           //Cd
      exponents.push_back(1.384108);                        //5s
      exponents.push_back(1.957413);                        //5p
    }
    else if (atomnum == 49) {                           //In
      exponents.push_back(2.023087);                        //5s
      exponents.push_back(2.106618);                        //5p
    }
    else if (atomnum == 50) {                           //Sn
      exponents.push_back(2.383941);                        //5s
      exponents.push_back(2.057908);                        //5p
    }
    else if (atomnum == 51) {                           //Sb
      exponents.push_back(2.391178);                        //5s
      exponents.push_back(1.773006);                        //5p
      exponents.push_back(2.465590);                        //5d
    }
    else if (atomnum == 52) {                           //Te
      exponents.push_back(2.769862);                        //5s
      exponents.push_back(1.731319);                        //5p
    }
    else if (atomnum == 53) {                           //I
      exponents.push_back(4.498653);                        //5s
      exponents.push_back(1.917072);                        //5p
      exponents.push_back(1.875175);                        //5d
    }
    else if (atomnum == 54) {                           //Xe
      exponents.push_back(2.759787);                        //5s
      exponents.push_back(1.977446);                        //5p
    }
    else if (atomnum == 55) {                           //Cs
      exponents.push_back(5.956008);                        //6s
      exponents.push_back(1.619485);                        //6p
    }
    else if (atomnum == 56) {                           //Ba
      exponents.push_back(1.395379);                        //6s
      exponents.push_back(1.430139);                        //6p
    }
    else if (atomnum == 57) {                           //La
      exponents.push_back(1.688562);                        //5d
      exponents.push_back(2.673780);                        //6s
      exponents.push_back(1.248192);                        //6p
    }
    else if (atomnum == 71) {                           //Lu
      exponents.push_back(2.225892);                        //5d
      exponents.push_back(5.471741);                        //6s
      exponents.push_back(1.712296);                        //6p
    }
    else if (atomnum == 72) {                           //Hf
      exponents.push_back(1.840840);                        //5d
      exponents.push_back(3.085344);                        //6s
      exponents.push_back(1.575819);                        //6p
    }
    else if (atomnum == 73) {                           //Ta
      exponents.push_back(1.838249);                        //5d
      exponents.push_back(4.578087);                        //6s
      exponents.push_back(4.841244);                        //6p
    }
    else if (atomnum == 74) {                           //W
      exponents.push_back(1.794400);                        //5d
      exponents.push_back(2.664560);                        //6s
      exponents.push_back(1.624010);                        //6p
    }
    else if (atomnum == 75) {                           //Re
      exponents.push_back(2.522766);                        //5d
      exponents.push_back(2.411839);                        //6s
      exponents.push_back(1.815351);                        //6p
    }
    else if (atomnum == 76) {                           //Os
      exponents.push_back(1.775570);                        //5d
      exponents.push_back(3.031000);                        //6s
      exponents.push_back(1.593960);                        //6p
    }
    else if (atomnum == 77) {                           //Ir
      exponents.push_back(2.676047);                        //5d
      exponents.push_back(1.500907);                        //6s
      exponents.push_back(4.106373);                        //6p
    }
    else if (atomnum == 78) {                           //Pt
      exponents.push_back(3.168852);                        //5d
      exponents.push_back(2.301264);                        //6s
      exponents.push_back(1.662404);                        //6p
    }
    else if (atomnum == 79) {                           //Au
      exponents.push_back(5.053167);                        //5d
      exponents.push_back(1.814169);                        //6s
      exponents.push_back(1.618657);                        //6p
    }
    else if (atomnum == 80) {                           //Hg
      exponents.push_back(2.104896);                        //6s
      exponents.push_back(1.516293);                        //6p
    }
    else if (atomnum == 81) {                           //Tl
      exponents.push_back(3.335883);                        //6s
      exponents.push_back(1.766141);                        //6p
    }
    else if (atomnum == 82) {                           //Pb
      exponents.push_back(2.368901);                        //6s
      exponents.push_back(1.685246);                        //6p
    }
    else if (atomnum == 83) {                           //Bi
      exponents.push_back(3.702377);                        //6s
      exponents.push_back(1.872327);                        //6p
    }
    else if (((atomnum > 57)&&(atomnum < 71))||(atomnum > 83)) {
      //unavailable: Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb
      //             Po,At,Rn,Fr,Ra,Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr,Rf,Db,Sg,Bh,Hs,Mt,Ds,Rg,Cn,Nh,Fl,Mc,Lv,Ts,Og
      throw("ERROR: basissets/dewar.hpp: PM6BSet: GetExponents(): orbital exponents not defined for atom");
    }
    return exponents;
  }
  int IPType() {return 5;}
  matrix<size_t> QuantumNumbers(size_t atomnum) {
    //calculation of quantum numbers for each atom
    std::vector<double> exps = GetExponents(atomnum);
    matrix<size_t> quantumnumb(exps.size(),2);
    if (atomnum > 86) {
      //unavailable: Fr,Ra,Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr,Rf,Db,Sg,Bh,Hs,Mt,Ds,Rg,Cn,Nh,Fl,Mc,Lv,Ts,Og
      throw("ERROR: basissets/dewar.hpp: PM6BSet: QuantumNumbers(): atom not defined");
    }
    if (atomnum == 1) quantumnumb(1,1) = 1;         //H
    else if (atomnum == 2) {                        //He
      quantumnumb(1,1) = 1;    //1s
      quantumnumb(2,1) = 2;    //2p
      quantumnumb(2,2) = 1;    //2p
    }
    else if ((atomnum > 2)&&(atomnum < 10)) {       //Li,Be,B,C,N,O,F
      quantumnumb(1,1) = 2;    //2s
      quantumnumb(2,1) = 2;    //2p
      quantumnumb(2,2) = 1;    //2p
    }
    else if (atomnum == 10) {                       //Ne
      quantumnumb(1,1) = 3;    //3s
      quantumnumb(2,1) = 2;    //2p
      quantumnumb(2,2) = 1;    //2p
    }
    else if ((atomnum > 10)&&(atomnum < 18)) {      //Na,Mg,Al,Si,P,S,Cl
      quantumnumb(1,1) = 3;    //3s
      quantumnumb(2,1) = 3;    //3p
      quantumnumb(2,2) = 1;    //3p
      if (atomnum > 12) {
        quantumnumb(3,1) = 3;  //3d
        quantumnumb(3,2) = 2;  //3d
      }
    }
    else if (atomnum == 18) {                       //Ar
      quantumnumb(1,1) = 4;    //4s
      quantumnumb(2,1) = 3;    //3p
      quantumnumb(2,2) = 1;    //3p
    }
    else if (((atomnum > 18)&&(atomnum < 21))||((atomnum > 29)&&(atomnum < 36))) {//K,Ca,Zn,Ga,Ge,As,Se,Br
      quantumnumb(1,1) = 4;    //4s
      quantumnumb(2,1) = 4;    //4p
      quantumnumb(2,2) = 1;    //4p
      if ((atomnum == 33)||(atomnum == 35)) {
        quantumnumb(3,1) = 4;  //4d
        quantumnumb(3,2) = 2;  //4d
      }
    }
    else if ((atomnum > 20)&&(atomnum < 30)) {      //Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu
      quantumnumb(1,1) = 3;    //3d
      quantumnumb(1,2) = 2;    //3d
      quantumnumb(2,1) = 4;    //4s
      quantumnumb(3,1) = 4;    //4p
      quantumnumb(3,2) = 1;    //4p
    }
    else if (atomnum == 36) {                       //Kr
      quantumnumb(1,1) = 5;    //5s
      quantumnumb(2,1) = 4;    //4p
      quantumnumb(2,2) = 1;    //4p
    }
    else if (((atomnum > 36)&&(atomnum < 39))||((atomnum > 47)&&(atomnum < 54))) {//Rb,Sr,Cd,In,Sn,Sb,Te,I
      quantumnumb(1,1) = 5;    //5s
      quantumnumb(2,1) = 5;    //5p
      quantumnumb(2,2) = 1;    //5p
      if ((atomnum == 51)||(atomnum == 53)) {
        quantumnumb(3,1) = 5;  //5d
        quantumnumb(3,2) = 2;  //5d
      }
    }
    else if ((atomnum > 38)&&(atomnum < 48)) {      //Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag
      quantumnumb(1,1) = 4;    //4d
      quantumnumb(1,2) = 2;    //4d
      quantumnumb(2,1) = 5;    //5s
      quantumnumb(3,1) = 5;    //5p
      quantumnumb(3,2) = 1;    //5p
    }
    else if (atomnum == 54) {                       //Xe
      quantumnumb(1,1) = 6;    //6s
      quantumnumb(2,1) = 5;    //5p
      quantumnumb(2,2) = 1;    //5p
    }
    else if (((atomnum > 54)&&(atomnum < 57))||(atomnum > 79)) {      //Cs,Ba,Hg,Tl,Pb,Bi
      quantumnumb(1,1) = 6;    //6s
      quantumnumb(2,1) = 6;    //6p
      quantumnumb(2,2) = 1;    //6p
    }
    else if ((atomnum > 56)&&(atomnum < 80)) {      //La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au
      quantumnumb(1,1) = 5;    //5d
      quantumnumb(1,2) = 2;    //5d
      quantumnumb(2,1) = 6;    //6s
      quantumnumb(3,1) = 6;    //6p
      quantumnumb(3,2) = 1;    //6p
    }
    return quantumnumb;
  }
  size_t GetNAOs(size_t atomnum) {
    //function to calculate the number of AOs used for a given atom
    std::vector<double> exps = GetExponents(atomnum);
    return exps.size()*exps.size();
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapPM6(atm);}
};
class PM3BSet: public STOBSet {
public:
  PM3BSet() {}
  ~PM3BSet() {}
  std::vector<double> GetExponents(size_t atomnum) {
    //getter for the exponents of a given atom
    std::vector<double> exponents;
    if (atomnum == 1) {exponents.push_back(0.967807);}  //H
    else if (atomnum == 3) {                            //Li
      exponents.push_back(0.650000);                        //2s
      exponents.push_back(0.750000);                        //2p
    }
    else if (atomnum == 4) {                            //Be
      exponents.push_back(0.877439);                        //2s
      exponents.push_back(1.508755);                        //2p
    }
    else if (atomnum == 6) {                            //C
      exponents.push_back(1.565085);                        //2s
      exponents.push_back(1.842345);                        //2p
    }
    else if (atomnum == 7) {                            //N
      exponents.push_back(2.028094);                        //2s
      exponents.push_back(2.313728);                        //2p
    }
    else if (atomnum == 8) {                            //O
      exponents.push_back(3.796544);                        //2s
      exponents.push_back(2.389402);                        //2p
    }
    else if (atomnum == 9) {                            //F
      exponents.push_back(4.708555);                        //2s
      exponents.push_back(2.491178);                        //2p
    }
    else if (atomnum == 11) {                           //Na
      exponents.push_back(1.1375011);                       //3s
      exponents.push_back(1.1877433);                       //3p
    }
    else if (atomnum == 12) {                           //Mg
      exponents.push_back(0.698552);                        //3s
      exponents.push_back(1.483453);                        //3p
    }
    else if (atomnum == 13) {                           //Al
      exponents.push_back(1.702888);                        //3s
      exponents.push_back(1.073629);                        //3p
    }
    else if (atomnum == 14) {                           //Si
      exponents.push_back(1.635075);                        //3s
      exponents.push_back(1.313088);                        //3p
    }
    else if (atomnum == 15) {                           //P
      exponents.push_back(2.017563);                        //3s
      exponents.push_back(1.504732);                        //3p
    }
    else if (atomnum == 16) {                           //S
      exponents.push_back(1.891185);                        //3s
      exponents.push_back(1.658972);                        //3p
    }
    else if (atomnum == 17) {                           //Cl
      exponents.push_back(2.246210);                        //3s
      exponents.push_back(2.151010);                        //3p
    }
    else if (atomnum == 30) {                           //Zn
      exponents.push_back(1.819989);                        //4s
      exponents.push_back(1.506922);                        //4p
    }
    else if (atomnum == 31) {                           //Ga
      exponents.push_back(1.847040);                        //4s
      exponents.push_back(0.839411);                        //4p
    }
    else if (atomnum == 32) {                           //Ge
      exponents.push_back(2.237353);                        //4s
      exponents.push_back(1.592432);                        //4p
    }
    else if (atomnum == 33) {                           //As
      exponents.push_back(2.636177);                        //4s
      exponents.push_back(1.703889);                        //4p
    }
    else if (atomnum == 34) {                           //Se
      exponents.push_back(2.828051);                        //4s
      exponents.push_back(1.732536);                        //4p
    }
    else if (atomnum == 35) {                           //Br
      exponents.push_back(5.348457);                        //4s
      exponents.push_back(2.127590);                        //4p
    }
    else if (atomnum == 48) {                           //Cd
      exponents.push_back(1.679351);                        //5s
      exponents.push_back(2.066412);                        //5p
    }
    else if (atomnum == 49) {                           //In
      exponents.push_back(2.016116);                        //5s
      exponents.push_back(1.445350);                        //5p
    }
    else if (atomnum == 50) {                           //Sn
      exponents.push_back(2.373328);                        //5s
      exponents.push_back(1.638233);                        //5p
    }
    else if (atomnum == 51) {                           //Sb
      exponents.push_back(2.343039);                        //5s
      exponents.push_back(1.899992);                        //5p
    }
    else if (atomnum == 52) {                           //Te
      exponents.push_back(4.165492);                        //5s
      exponents.push_back(1.647555);                        //5p
    }
    else if (atomnum == 53) {                           //I
      exponents.push_back(7.001013);                        //5s
      exponents.push_back(2.454354);                        //5p
    }
    else if (atomnum == 80) {                           //Hg
      exponents.push_back(1.476885);                        //6s
      exponents.push_back(2.479951);                        //6p
    }
    else if (atomnum == 81) {                           //Tl
      exponents.push_back(6.867921);                        //6s
      exponents.push_back(1.969445);                        //6p
    }
    else if (atomnum == 82) {                           //Pb
      exponents.push_back(3.141289);                        //6s
      exponents.push_back(1.892418);                        //6p
    }
    else if (atomnum == 83) {                           //Bi
      exponents.push_back(4.916451);                        //6s
      exponents.push_back(1.934935);                        //6p
    }
    else {                                              //unavailable: He,B,Ne,Ar-Cu,Kr-Ag,Xe-Au,Po,At,Rn,Fr,Ra,Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr,Rf,Db,Sg,Bh,Hs,Mt,Ds,Rg,Cn,Nh,Fl,Mc,Lv,Ts,Og
      throw("ERROR: basissets/dewar.hpp: PM3BSet: GetExponents(): orbital exponents not defined for atom");
    }
    return exponents;
  }
  matrix<size_t> QuantumNumbers(size_t atomnum) {
    //calculation of quantum numbers for each atom
    std::vector<double> exps = GetExponents(atomnum);
    matrix<size_t> quantumnumb(exps.size(),2);
    if (atomnum == 1) quantumnumb(1,1) = 1;         //H
    else if ((atomnum > 2)&&(atomnum < 10)) {       //Li,Be,C,N,O,F
      quantumnumb(1,1) = 2;    //2s
      quantumnumb(2,1) = 2;    //2p
      quantumnumb(2,2) = 1;    //2p
    }
    else if ((atomnum > 10)&&(atomnum < 18)) {      //Na,Mg,Al,Si,P,S,Cl
      quantumnumb(1,1) = 3;    //3s
      quantumnumb(2,1) = 3;    //3p
      quantumnumb(2,2) = 1;    //3p
    }
    else if ((atomnum > 29)&&(atomnum < 36)) {      //Zn,Ga,Ge,As,Se,Br
      quantumnumb(1,1) = 4;    //4s
      quantumnumb(2,1) = 4;    //4p
      quantumnumb(2,2) = 1;    //4p
    }
    else if ((atomnum > 47)&&(atomnum < 54)) {      //Cd,In,Sn,Sb,Te,I
      quantumnumb(1,1) = 5;    //5s
      quantumnumb(2,1) = 5;    //5p
      quantumnumb(2,2) = 1;    //5p
    }
    else if (atomnum > 79) {                        //Hg,Tl,Pb,Bi
      quantumnumb(1,1) = 6;    //6s
      quantumnumb(2,1) = 6;    //6p
      quantumnumb(2,2) = 1;    //6p
    }
    return quantumnumb;
  }
  size_t GetNAOs(size_t atomnum) {
    //function to calculate the number of AOs used for a given atom
    std::vector<double> exps = GetExponents(atomnum);
    return exps.size()*exps.size();
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapPM3(atm);}
};
class PM3BPBSet: public STOBSet {
public:
  PM3BPBSet() {}
  ~PM3BPBSet() {}
  std::vector<double> GetExponents(size_t atomnum) {
    //getter for the exponents of a given atom
    std::vector<double> exponents;
    if (atomnum == 1) {exponents.push_back(0.967807);}  //H
    else if (atomnum == 6) {                            //C
      exponents.push_back(1.565085);                        //2s
      exponents.push_back(1.842345);                        //2p
    }
    else if (atomnum == 7) {                            //N
      exponents.push_back(2.028094);                        //2s
      exponents.push_back(2.313728);                        //2p
    }
    else if (atomnum == 8) {                            //O
      exponents.push_back(3.796544);                        //2s
      exponents.push_back(2.389402);                        //2p
    }
    else {                                              //unavailable: kind of everything else
      throw("ERROR: basissets/dewar.hpp: PM3BPBSet: GetExponents(): orbital exponents not defined for atom");
    }
    return exponents;
  }
  matrix<size_t> QuantumNumbers(size_t atomnum) {
    //calculation of quantum numbers for each atom
    std::vector<double> exps = GetExponents(atomnum);
    matrix<size_t> quantumnumb(exps.size(),2);
    if (atomnum == 1) quantumnumb(1,1) = 1;         //H
    else if ((atomnum > 5)&&(atomnum < 9)) {        //Li,Be,C,N,O,F
      quantumnumb(1,1) = 2;    //2s
      quantumnumb(2,1) = 2;    //2p
      quantumnumb(2,2) = 1;    //2p
    }
    return quantumnumb;
  }
  size_t GetNAOs(size_t atomnum) {
    //function to calculate the number of AOs used for a given atom
    std::vector<double> exps = GetExponents(atomnum);
    return exps.size()*exps.size();
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapPM3(atm);}
};
class PM3PDDGBSet: public STOBSet {
public:
  PM3PDDGBSet() {}
  ~PM3PDDGBSet() {}
  std::vector<double> GetExponents(size_t atomnum) {
    //getter for the exponents of a given atom
    std::vector<double> exponents;
    if (atomnum == 1) {exponents.push_back(0.972786);}  //H
    else if (atomnum == 6) {                            //C
      exponents.push_back(1.567864);                        //2s
      exponents.push_back(1.846659);                        //2p
    }
    else if (atomnum == 7) {                            //N
      exponents.push_back(2.035807);                        //2s
      exponents.push_back(2.324327);                        //2p
    }
    else if (atomnum == 8) {                            //O
      exponents.push_back(3.814565);                        //2s
      exponents.push_back(2.318011);                        //2p
    }
    else if (atomnum == 9) {                            //F
      exponents.push_back(5.538033);                        //2s
      exponents.push_back(2.538066);                        //2p
    }
    else if (atomnum == 14) {                           //Si
      exponents.push_back(1.586389);                        //3s
      exponents.push_back(1.485958);                        //3p
    }
    else if (atomnum == 15) {                           //P
      exponents.push_back(2.395882);                        //3s
      exponents.push_back(1.742213);                        //3p
    }
    else if (atomnum == 16) {                           //S
      exponents.push_back(1.012002);                        //3s
      exponents.push_back(1.876999);                        //3p
    }
    else if (atomnum == 17) {                           //Cl
      exponents.push_back(2.548268);                        //3s
      exponents.push_back(2.284624);                        //3p
    }
    else if (atomnum == 35) {                           //Br
      exponents.push_back(4.345079);                        //4s
      exponents.push_back(2.190961);                        //4p
    }
    else if (atomnum == 53) {                           //I
      exponents.push_back(5.062801);                        //5s
      exponents.push_back(2.417757);                        //5p
    }
    else {
      throw("ERROR: basissets/dewar.hpp: PM3PDDGBSet: GetExponents(): orbital exponents not defined for atom");
    }
    return exponents;
  }
  matrix<size_t> QuantumNumbers(size_t atomnum) {
    //calculation of quantum numbers for each atom
    std::vector<double> exps = GetExponents(atomnum);
    matrix<size_t> quantumnumb(exps.size(),2);
    if (atomnum == 1) quantumnumb(1,1) = 1;         //H
    else if ((atomnum > 5)&&(atomnum < 10)) {       //C,N,O,F
      quantumnumb(1,1) = 2;    //2s
      quantumnumb(2,1) = 2;    //2p
      quantumnumb(2,2) = 1;    //2p
    }
    else if ((atomnum > 13)&&(atomnum < 18)) {      //Si,P,S,Cl
      quantumnumb(1,1) = 3;    //3s
      quantumnumb(2,1) = 3;    //3p
      quantumnumb(2,2) = 1;    //3p
    }
    else if (atomnum == 35) {                       //Br
      quantumnumb(1,1) = 4;    //4s
      quantumnumb(2,1) = 4;    //4p
      quantumnumb(2,2) = 1;    //4p
    }
    else if (atomnum == 53) {                       //I
      quantumnumb(1,1) = 5;    //5s
      quantumnumb(2,1) = 5;    //5p
      quantumnumb(2,2) = 1;    //5p
    }
    return quantumnumb;
  }
  size_t GetNAOs(size_t atomnum) {
    //function to calculate the number of AOs used for a given atom
    std::vector<double> exps = GetExponents(atomnum);
    return exps.size()*exps.size();
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapPM3PDDG(atm);}
};
class MNDOPDDGBSet: public STOBSet {
public:
  MNDOPDDGBSet() {}
  ~MNDOPDDGBSet() {}
  std::vector<double> GetExponents(size_t atomnum) {
    //getter for the exponents of a given atom
    std::vector<double> exponents;
    if (atomnum == 1) {exponents.push_back(1.322431);}  //H
    else if (atomnum == 6) {                            //C
      exponents.push_back(1.809817);                        //2s
      exponents.push_back(1.825008);                        //2p
    }
    else if (atomnum == 7) {                            //N
      exponents.push_back(2.231424);                        //2s
      exponents.push_back(2.253460);                        //2p
    }
    else if (atomnum == 8) {                            //O
      exponents.push_back(2.569172);                        //2s
      exponents.push_back(2.697152);                        //2p
    }
    else if (atomnum == 9) {                            //F
      exponents.push_back(4.328519);                        //2s
      exponents.push_back(2.905042);                        //2p
    }
    else if (atomnum == 17) {                           //Cl
      exponents.push_back(4.212404);                        //3s
      exponents.push_back(2.037647);                        //3p
    }
    else if (atomnum == 35) {                           //Br
      exponents.push_back(3.999975);                        //4s
      exponents.push_back(2.245040);                        //4p
    }
    else if (atomnum == 53) {                           //I
      exponents.push_back(2.718404);                        //5s
      exponents.push_back(2.461813);                        //5p
    }
    else {
      throw("ERROR: basissets/dewar.hpp: MNDOPDDGBSet: GetExponents(): orbital exponents not defined for atom");
    }
    return exponents;
  }
  matrix<size_t> QuantumNumbers(size_t atomnum) {
    //calculation of quantum numbers for each atom
    std::vector<double> exps = GetExponents(atomnum);
    matrix<size_t> quantumnumb(exps.size(),2);
    if (atomnum == 1) quantumnumb(1,1) = 1;         //H
    else if ((atomnum > 5)&&(atomnum < 10)) {       //C,N,O,F
      quantumnumb(1,1) = 2;    //2s
      quantumnumb(2,1) = 2;    //2p
      quantumnumb(2,2) = 1;    //2p
    }
    else if (atomnum == 17) {                       //Cl
      quantumnumb(1,1) = 3;    //3s
      quantumnumb(2,1) = 3;    //3p
      quantumnumb(2,2) = 1;    //3p
    }
    else if (atomnum == 35) {                       //Br
      quantumnumb(1,1) = 4;    //4s
      quantumnumb(2,1) = 4;    //4p
      quantumnumb(2,2) = 1;    //4p
    }
    else if (atomnum == 53) {                       //I
      quantumnumb(1,1) = 5;    //5s
      quantumnumb(2,1) = 5;    //5p
      quantumnumb(2,2) = 1;    //5p
    }
    return quantumnumb;
  }
  size_t GetNAOs(size_t atomnum) {
    //function to calculate the number of AOs used for a given atom
    std::vector<double> exps = GetExponents(atomnum);
    return exps.size()*exps.size();
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapMNDOPDDG(atm);}
};
class RM1BSet: public STOBSet {
public:
  RM1BSet() {}
  ~RM1BSet() {}
  std::vector<double> GetExponents(size_t atomnum) {
    //getter for the exponents of a given atom
    std::vector<double> exponents;
    if (atomnum == 1) {exponents.push_back(1.08267366);}//H
    else if (atomnum == 6) {                            //C
      exponents.push_back(1.85018803);                        //2s
      exponents.push_back(1.76830093);                        //2p
    }
    else if (atomnum == 7) {                            //N
      exponents.push_back(2.37447159);                        //2s
      exponents.push_back(1.97812569);                        //2p
    }
    else if (atomnum == 8) {                            //O
      exponents.push_back(3.17936914);                        //2s
      exponents.push_back(2.55361907);                        //2p
    }
    else if (atomnum == 9) {                            //F
      exponents.push_back(4.40337913);                        //2s
      exponents.push_back(2.64841556);                        //2p
    }
    else if (atomnum == 15) {                           //P
      exponents.push_back(2.12240118);                        //3s
      exponents.push_back(1.74327954);                        //3p
    }
    else if (atomnum == 16) {                           //S
      exponents.push_back(2.13344308);                        //3s
      exponents.push_back(1.87460650);                        //3p
    }
    else if (atomnum == 17) {                           //Cl
      exponents.push_back(3.86491071);                        //3s
      exponents.push_back(1.89593144);                        //3p
    }
    else if (atomnum == 35) {                           //Br
      exponents.push_back(5.73157215);                        //4s
      exponents.push_back(2.03147582);                        //4p
    }
    else if (atomnum == 53) {                           //I
      exponents.push_back(2.53003753);                        //5s
      exponents.push_back(2.31738678);                        //5p
    }
    else {                                              //available atoms: H,C,N,O,F,P,S,Cl,Br,I
      throw("ERROR: basissets/dewar.hpp: RM1BSet: GetExponents(): orbital exponents not defined for atom");
    }
    return exponents;
  }
  matrix<size_t> QuantumNumbers(size_t atomnum) {
    //calculation of quantum numbers for each atom
    std::vector<double> exps = GetExponents(atomnum);
    matrix<size_t> quantumnumb(exps.size(),2);
    if (atomnum == 1) quantumnumb(1,1) = 1;         //H
    else if ((atomnum > 5)&&(atomnum < 10)) {       //C,N,O,F
      quantumnumb(1,1) = 2;    //2s
      quantumnumb(2,1) = 2;    //2p
      quantumnumb(2,2) = 1;    //2p
    }
    else if ((atomnum > 14)&&(atomnum < 18)) {      //P,S,Cl
      quantumnumb(1,1) = 3;    //3s
      quantumnumb(2,1) = 3;    //3p
      quantumnumb(2,2) = 1;    //3p
    }
    else if (atomnum == 35) {                       //Br
      quantumnumb(1,1) = 4;    //4s
      quantumnumb(2,1) = 4;    //4p
      quantumnumb(2,2) = 1;    //4p
    }
    else if (atomnum == 53) {                       //I
      quantumnumb(1,1) = 5;    //5s
      quantumnumb(2,1) = 5;    //5p
      quantumnumb(2,2) = 1;    //5p
    }
    return quantumnumb;
  }
  size_t GetNAOs(size_t atomnum) {
    //function to calculate the number of AOs used for a given atom
    std::vector<double> exps = GetExponents(atomnum);
    return exps.size()*exps.size();
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapRM1(atm);}
};
class AM1BSet: public STOBSet {
public:
  AM1BSet() {}
  ~AM1BSet() {}
  std::vector<double> GetExponents(size_t atomnum) {
    //getter for the exponents of a given atom
    std::vector<double> exponents;
    if (atomnum == 1) {exponents.push_back(1.188078);}  //H
    else if (atomnum == 3) {                            //Li
      exponents.push_back(0.70238);                         //2s
      exponents.push_back(0.70238);                         //2p
    }
    else if (atomnum == 4) {                            //Be
      exponents.push_back(1.00421);                         //2s
      exponents.push_back(1.00421);                         //2p
    }
    else if (atomnum == 5) {                            //B
      exponents.push_back(1.611709);                        //2s
      exponents.push_back(1.555385);                        //2p
    }
    else if (atomnum == 6) {                            //C
      exponents.push_back(1.808665);                        //2s
      exponents.push_back(1.685116);                        //2p
    }
    else if (atomnum == 7) {                            //N
      exponents.push_back(2.315410);                        //2s
      exponents.push_back(2.157940);                        //2p
    }
    else if (atomnum == 8) {                            //O
      exponents.push_back(3.108032);                        //2s
      exponents.push_back(2.524039);                        //2p
    }
    else if (atomnum == 9) {                            //F
      exponents.push_back(3.770082);                        //2s
      exponents.push_back(2.494670);                        //2p
    }
    else if (atomnum == 11) {                           //Na
      exponents.push_back(0.6797792);                       //3s
      exponents.push_back(1.2170468);                       //3p
    }
    else if (atomnum == 13) {                           //Al
      exponents.push_back(1.516593);                        //3s
      exponents.push_back(1.306347);                        //3p
    }
    else if (atomnum == 14) {                           //Si
      exponents.push_back(1.830697);                        //3s
      exponents.push_back(1.284953);                        //3p
    }
    else if (atomnum == 15) {                           //P
      exponents.push_back(1.98128);                         //3s
      exponents.push_back(1.87515);                         //3p
    }
    else if (atomnum == 16) {                           //S
      exponents.push_back(2.366515);                        //3s
      exponents.push_back(1.667263);                        //3p
    }
    else if (atomnum == 17) {                           //Cl
      exponents.push_back(3.631376);                        //3s
      exponents.push_back(2.076799);                        //3p
    }
    else if (atomnum == 30) {                           //Zn
      exponents.push_back(1.954299);                        //4s
      exponents.push_back(1.372365);                        //4p
    }
    else if (atomnum == 32) {                           //Ge
      exponents.push_back(1.219631);                        //4s
      exponents.push_back(1.982794);                        //4p
    }
    else if (atomnum == 33) {                           //As
      exponents.push_back(2.2576897);                       //4s
      exponents.push_back(1.7249710);                       //4p
    }
    else if (atomnum == 34) {                           //Se
      exponents.push_back(2.6841570);                       //4s
      exponents.push_back(2.0506164);                       //4p
    }
    else if (atomnum == 35) {                           //Br
      exponents.push_back(3.064133);                        //4s
      exponents.push_back(2.038333);                        //4p
    }
    else if (atomnum == 51) {                           //Sb
      exponents.push_back(2.254823);                        //5s
      exponents.push_back(2.218592);                        //5p
    }
    else if (atomnum == 52) {                           //Te
      exponents.push_back(2.1321165);                       //5s
      exponents.push_back(1.9712680);                       //5p
    }
    else if (atomnum == 53) {                           //I
      exponents.push_back(2.102858);                        //5s
      exponents.push_back(2.161153);                        //5p
    }
    else if (atomnum == 80) {                           //Hg
      exponents.push_back(2.036413);                        //6s
      exponents.push_back(1.955766);                        //6p
    }
    else {                                              //available: H,Li,Be,B,C,N,O,F,Na,Al,Si,P,S,Cl,Zn,Ge,As,Se,Br,Sb,Te,I,Hg
      throw("ERROR: basissets/dewar.hpp: AM1BSet: GetExponents(): orbital exponents not defined for atom");
    }
    return exponents;
  }
  matrix<size_t> QuantumNumbers(size_t atomnum) {
    //calculation of quantum numbers for each atom
    std::vector<double> exps = GetExponents(atomnum);
    matrix<size_t> quantumnumb(exps.size(),2);
    if (atomnum == 1) quantumnumb(1,1) = 1;         //H
    else if ((atomnum > 2)&&(atomnum < 10)) {       //Li,Be,C,N,O,F
      quantumnumb(1,1) = 2;    //2s
      quantumnumb(2,1) = 2;    //2p
      quantumnumb(2,2) = 1;    //2p
    }
    else if ((atomnum > 10)&&(atomnum < 18)) {      //Na,Al,Si,P,S,Cl
      quantumnumb(1,1) = 3;    //3s
      quantumnumb(2,1) = 3;    //3p
      quantumnumb(2,2) = 1;    //3p
    }
    else if ((atomnum > 29)&&(atomnum < 36)) {      //Zn,Ge,As,Se,Br
      quantumnumb(1,1) = 4;    //4s
      quantumnumb(2,1) = 4;    //4p
      quantumnumb(2,2) = 1;    //4p
    }
    else if ((atomnum > 50)&&(atomnum < 54)) {      //Sb,Te,I
      quantumnumb(1,1) = 5;    //5s
      quantumnumb(2,1) = 5;    //5p
      quantumnumb(2,2) = 1;    //5p
    }
    else if (atomnum == 80) {                       //Hg
      quantumnumb(1,1) = 6;    //6s
      quantumnumb(2,1) = 6;    //6p
      quantumnumb(2,2) = 1;    //6p
    }
    return quantumnumb;
  }
  size_t GetNAOs(size_t atomnum) {
    //function to calculate the number of AOs used for a given atom
    std::vector<double> exps = GetExponents(atomnum);
    return exps.size()*exps.size();
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapAM1(atm);}
};
class MNDOdBSet: public STOBSet {
public:
  MNDOdBSet() {}
  ~MNDOdBSet() {}
  std::vector<double> GetExponents(size_t atomnum) {
    //getter for the exponents of a given atom
    std::vector<double> exponents;
    if (atomnum == 1) {exponents.push_back(1.331967);}  //H
    else if (atomnum == 3) {                            //Li
      exponents.push_back(0.70238);                         //2s
      exponents.push_back(0.70238);                         //2p
    }
    else if (atomnum == 4) {                            //Be
      exponents.push_back(1.00421);                         //2s
      exponents.push_back(1.00421);                         //2p
    }
    else if (atomnum == 5) {                            //B
      exponents.push_back(1.506801);                        //2s
      exponents.push_back(1.506801);                        //2p
    }
    else if (atomnum == 6) {                            //C
      exponents.push_back(1.787537);                        //2s
      exponents.push_back(1.787537);                        //2p
    }
    else if (atomnum == 7) {                            //N
      exponents.push_back(2.255614);                        //2s
      exponents.push_back(2.255614);                        //2p
    }
    else if (atomnum == 8) {                            //O
      exponents.push_back(2.699905);                        //2s
      exponents.push_back(2.699905);                        //2p
    }
    else if (atomnum == 9) {                            //F
      exponents.push_back(2.848487);                        //2s
      exponents.push_back(2.848487);                        //2p
    }
    else if (atomnum == 11) {                           //Na
      exponents.push_back(0.987508);                        //3s
      exponents.push_back(0.893350);                        //3p
    }
    else if (atomnum == 12) {                           //Mg
      exponents.push_back(1.448904);                        //3s
      exponents.push_back(0.952930);                        //3p
    }
    else if (atomnum == 13) {                           //Al
      exponents.push_back(1.625160);                        //3s
      exponents.push_back(1.319936);                        //3p
      exponents.push_back(0.952667);                        //3d
    }
    else if (atomnum == 14) {                           //Si
      exponents.push_back(1.915655);                        //3s
      exponents.push_back(1.681611);                        //3p
      exponents.push_back(0.966772);                        //3d
    }
    else if (atomnum == 15) {                           //P
      exponents.push_back(2.266463);                        //3s
      exponents.push_back(1.940015);                        //3p
      exponents.push_back(1.100109);                        //3d
    }
    else if (atomnum == 16) {                           //S
      exponents.push_back(2.225851);                        //3s
      exponents.push_back(2.099706);                        //3p
      exponents.push_back(1.231472);                        //3d
    }
    else if (atomnum == 17) {                           //Cl
      exponents.push_back(2.561611);                        //3s
      exponents.push_back(2.389338);                        //3p
      exponents.push_back(1.251398);                        //3d
    }
    else if (atomnum == 30) {                           //Zn
      exponents.push_back(1.731503);                        //4s
      exponents.push_back(1.393583);                        //4p
    }
    else if (atomnum == 32) {                           //Ge
      exponents.push_back(1.293180);                        //4s
      exponents.push_back(2.020564);                        //4p
    }
    else if (atomnum == 35) {                           //Br
      exponents.push_back(2.590541);                        //4s
      exponents.push_back(2.330857);                        //4p
      exponents.push_back(1.357361);                        //4d
    }
    else if (atomnum == 48) {                           //Cd
      exponents.push_back(1.748806);                        //5s
      exponents.push_back(1.563215);                        //5p
    }
    else if (atomnum == 50) {                           //Sn
      exponents.push_back(2.080380);                        //5s
      exponents.push_back(1.937106);                        //5p
    }
    else if (atomnum == 53) {                           //I
      exponents.push_back(2.756543);                        //5s
      exponents.push_back(2.253080);                        //5p
      exponents.push_back(1.502335);                        //5d
    }
    else if (atomnum == 80) {                           //Hg
      exponents.push_back(2.333107);                        //6s
      exponents.push_back(1.708311);                        //6p
    }
    else if (atomnum == 82) {                           //Pb
      exponents.push_back(2.498286);                        //6s
      exponents.push_back(2.082071);                        //6p
    }
    else {                                              //available atoms: H,Li,Be,B,C,N,O,F,Al,Si,P,S,Cl,Zn,Ge,Br,Cd,Sn,I,Hg,Pb
      throw("ERROR: basissets/dewar.hpp: MNDOdBSet: GetExponents(): orbital exponents not defined for atom");
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
    else if ((atomnum > 10)&&(atomnum < 18)) {       //Na,Mg,Al,Si,P,S,Cl
      quantumnumb(1,1) = 3;    //3s
      quantumnumb(2,1) = 3;    //3p
      quantumnumb(2,2) = 1;    //3p
      if (atomnum > 12) {
        quantumnumb(3,1) = 3;  //3d
        quantumnumb(3,2) = 2;  //3d
      }
    }
    else if ((atomnum > 29)&&(atomnum < 36)) {       //Zn,Ge,Br
      quantumnumb(1,1) = 4;    //4s
      quantumnumb(2,1) = 4;    //4p
      quantumnumb(2,2) = 1;    //4p
      if (atomnum == 35) {
        quantumnumb(3,1) = 4;  //4d
        quantumnumb(3,2) = 2;  //4d
      }
    }
    else if ((atomnum > 47)&&(atomnum < 54)) {       //Cd,Sn,I
      quantumnumb(1,1) = 5;    //5s
      quantumnumb(2,1) = 5;    //5p
      quantumnumb(2,2) = 1;    //5p
      if (atomnum == 53) {
        quantumnumb(3,1) = 5;  //5d
        quantumnumb(3,2) = 2;  //5d
      }
    }
    else if (atomnum > 79) {                         //Hg,Pb
      quantumnumb(1,1) = 6;    //6s
      quantumnumb(2,1) = 6;    //6p
      quantumnumb(2,2) = 1;    //6p
    }
    return quantumnumb;
  }
  size_t GetNAOs(size_t atomnum) {
    //function to calculate the number of AOs used for a given atom
    std::vector<double> exps = GetExponents(atomnum);
    return exps.size()*exps.size();
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapMNDOd(atm);}
};
class MNDOBSet: public STOBSet {
public:
  MNDOBSet() {}
  ~MNDOBSet() {}
  std::vector<double> GetExponents(size_t atomnum) {
    //getter for the exponents of a given atom
    std::vector<double> exponents;
    if (atomnum == 1) {exponents.push_back(1.331967);}  //H
    else if (atomnum == 3) {                            //Li
      exponents.push_back(0.70238);                         //2s
      exponents.push_back(0.70238);                         //2p
    }
    else if (atomnum == 4) {                            //Be
      exponents.push_back(1.00421);                         //2s
      exponents.push_back(1.00421);                         //2p
    }
    else if (atomnum == 5) {                            //B
      exponents.push_back(1.506801);                        //2s
      exponents.push_back(1.506801);                        //2p
    }
    else if (atomnum == 6) {                            //C
      exponents.push_back(1.787537);                        //2s
      exponents.push_back(1.787537);                        //2p
    }
    else if (atomnum == 7) {                            //N
      exponents.push_back(2.255614);                        //2s
      exponents.push_back(2.255614);                        //2p
    }
    else if (atomnum == 8) {                            //O
      exponents.push_back(2.699905);                        //2s
      exponents.push_back(2.699905);                        //2p
    }
    else if (atomnum == 9) {                            //F
      exponents.push_back(2.848487);                        //2s
      exponents.push_back(2.848487);                        //2p
    }
    else if (atomnum == 13) {                           //Al
      exponents.push_back(1.444161);                        //3s
      exponents.push_back(1.444161);                        //3p
    }
    else if (atomnum == 14) {                           //Si
      exponents.push_back(1.315986);                        //3s
      exponents.push_back(1.709943);                        //3p
    }
    else if (atomnum == 15) {                           //P
      exponents.push_back(2.108720);                        //3s
      exponents.push_back(1.785810);                        //3p
    }
    else if (atomnum == 16) {                           //S
      exponents.push_back(2.312962);                        //3s
      exponents.push_back(2.009146);                        //3p
    }
    else if (atomnum == 17) {                           //Cl
      exponents.push_back(3.784645);                        //3s
      exponents.push_back(2.036263);                        //3p
    }
    else if (atomnum == 30) {                           //Zn
      exponents.push_back(2.047359);                        //4s
      exponents.push_back(1.460946);                        //4p
    }
    else if (atomnum == 32) {                           //Ge
      exponents.push_back(1.293180);                        //4s
      exponents.push_back(2.020564);                        //4p
    }
    else if (atomnum == 35) {                           //Br
      exponents.push_back(3.85430190);                      //4s
      exponents.push_back(2.19920914);                      //4p
    }
    else if (atomnum == 50) {                           //Sn
      exponents.push_back(2.080380);                        //5s
      exponents.push_back(1.937106);                        //5p
    }
    else if (atomnum == 53) {                           //I
      exponents.push_back(2.27296101);                      //5s
      exponents.push_back(2.16949803);                      //5p
    }
    else if (atomnum == 80) {                           //Hg
      exponents.push_back(2.218184);                        //6s
      exponents.push_back(2.065038);                        //6p
    }
    else if (atomnum == 82) {                           //Pb
      exponents.push_back(2.498286);                        //6s
      exponents.push_back(2.082071);                        //6p
    }
    else {                                              //available atoms: H,Li,Be,B,C,N,O,F,Al,Si,P,S,Cl,Zn,Ge,Br,Sn,I,Hg,Pb
      throw("ERROR: basissets/dewar.hpp: MNDOBSet: GetExponents(): orbital exponents not defined for atom");
    }
    return exponents;
  }
  matrix<size_t> QuantumNumbers(size_t atomnum) {
    //calculation of quantum numbers for each atom
    std::vector<double> exps = GetExponents(atomnum);
    matrix<size_t> quantumnumb(exps.size(),2);
    if (atomnum == 1) quantumnumb(1,1) = 1;         //H
    else if ((atomnum > 2)&&(atomnum < 10)) {       //Li,Be,C,N,O,F
      quantumnumb(1,1) = 2;    //2s
      quantumnumb(2,1) = 2;    //2p
      quantumnumb(2,2) = 1;    //2p
    }
    else if ((atomnum > 12)&&(atomnum < 18)) {      //Al,Si,P,S,Cl
      quantumnumb(1,1) = 3;    //3s
      quantumnumb(2,1) = 3;    //3p
      quantumnumb(2,2) = 1;    //3p
    }
    else if ((atomnum > 29)&&(atomnum < 36)) {      //Zn,Ge,Br
      quantumnumb(1,1) = 4;    //4s
      quantumnumb(2,1) = 4;    //4p
      quantumnumb(2,2) = 1;    //4p
    }
    else if ((atomnum > 49)&&(atomnum < 54)) {      //Sn,I
      quantumnumb(1,1) = 5;    //5s
      quantumnumb(2,1) = 5;    //5p
      quantumnumb(2,2) = 1;    //5p
    }
    else if ((atomnum > 79)&&(atomnum < 83)) {      //Hg,Pb
      quantumnumb(1,1) = 6;    //6s
      quantumnumb(2,1) = 6;    //6p
      quantumnumb(2,2) = 1;    //6p
    }
    return quantumnumb;
  }
  size_t GetNAOs(size_t atomnum) {
    //function to calculate the number of AOs used for a given atom
    std::vector<double> exps = GetExponents(atomnum);
    return exps.size()*exps.size();
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapMNDO(atm);}
};
#endif //_Dewar_Basis_Set_
