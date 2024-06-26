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

#ifndef _GFN2_Basis_Set_
#define _GFN2_Basis_Set_
#include <vector>
#include <string>
#include <math.h>
#include "STO.hpp"
#include "../parameters/ZeroOvpar.hpp"

//description:
//exponents for GFN2 basis set

class GFN2BSet: public STOBSet {
public:
  GFN2BSet() {}
  ~GFN2BSet() {}
  std::vector<double> GetExponents(size_t atomnum) {
    //getter for the exponents of a given atom
    std::vector<double> exponents;
    if (atomnum == 1) {exponents.push_back(1.230000);}  //H
    else if (atomnum == 2) {                            //He
      exponents.push_back(1.669667);
      exponents.push_back(1.500000);
    }
    else if (atomnum == 3) {                            //Li
      exponents.push_back(0.750060);
      exponents.push_back(0.557848);
    }
    else if (atomnum == 4) {                            //Be
      exponents.push_back(1.034720);
      exponents.push_back(0.949332);
    }
    else if (atomnum == 5) {                            //B
      exponents.push_back(1.479444);
      exponents.push_back(1.479805);
    }
    else if (atomnum == 6) {                            //C
      exponents.push_back(2.096432);
      exponents.push_back(1.800000);
    }
    else if (atomnum == 7) {                            //N
      exponents.push_back(2.339881);
      exponents.push_back(2.014332);
    }
    else if (atomnum == 8) {                            //O
      exponents.push_back(2.439742);
      exponents.push_back(2.137023);
    }
    else if (atomnum == 9) {                            //F
      exponents.push_back(2.416361);
      exponents.push_back(2.308399);
    }
    else if (atomnum == 10) {                           //Ne
      exponents.push_back(3.084104);
      exponents.push_back(2.312051);
      exponents.push_back(2.815609);
    }
    else if (atomnum == 11) {                           //Na
      exponents.push_back(0.763787);
      exponents.push_back(0.573553);
    }
    else if (atomnum == 12) {                           //Mg
      exponents.push_back(1.184203);
      exponents.push_back(0.717769);
      exponents.push_back(1.300000);
    }
    else if (atomnum == 13) {                           //Al
      exponents.push_back(1.352531);
      exponents.push_back(1.391201);
      exponents.push_back(1.000000);
    }
    else if (atomnum == 14) {                           //Si
      exponents.push_back(1.773917);
      exponents.push_back(1.718996);
      exponents.push_back(1.250000);
    }
    else if (atomnum == 15) {                           //P
      exponents.push_back(1.816945);
      exponents.push_back(1.903247);
      exponents.push_back(1.167533);
    }
    else if (atomnum == 16) {                           //S
      exponents.push_back(1.981333);
      exponents.push_back(2.025643);
      exponents.push_back(1.702555);
    }
    else if (atomnum == 17) {                           //Cl
      exponents.push_back(2.485265);
      exponents.push_back(2.199650);
      exponents.push_back(2.476089);
    }
    else if (atomnum == 18) {                           //Ar
      exponents.push_back(2.329679);
      exponents.push_back(2.149419);
      exponents.push_back(1.950531);
    }
    else if (atomnum == 19) {                           //K
      exponents.push_back(0.875961);
      exponents.push_back(0.631694);
    }
    else if (atomnum == 20) {                           //Ca
      exponents.push_back(1.267130);
      exponents.push_back(0.786247);
      exponents.push_back(1.380000);
    }
    else if (atomnum == 21) {                           //Sc
      exponents.push_back(2.440000);
      exponents.push_back(1.358701);
      exponents.push_back(1.019252);
    }
    else if (atomnum == 22) {                           //Ti
      exponents.push_back(1.849994);
      exponents.push_back(1.469983);
      exponents.push_back(0.957410);
    }
    else if (atomnum == 23) {                           //V
      exponents.push_back(1.673577);
      exponents.push_back(1.383176);
      exponents.push_back(0.938025);
    }
    else if (atomnum == 24) {                           //Cr
      exponents.push_back(1.568211);
      exponents.push_back(1.395427);
      exponents.push_back(1.080270);
    }
    else if (atomnum == 25) {                           //Mn
      exponents.push_back(1.839250);
      exponents.push_back(1.222190);
      exponents.push_back(1.240215);
    }
    else if (atomnum == 26) {                           //Fe
      exponents.push_back(1.911049);
      exponents.push_back(1.022393);
      exponents.push_back(1.294467);
    }
    else if (atomnum == 27) {                           //Co
      exponents.push_back(2.326507);
      exponents.push_back(1.464221);
      exponents.push_back(1.298678);
    }
    else if (atomnum == 28) {                           //Ni
      exponents.push_back(2.430756);
      exponents.push_back(1.469945);
      exponents.push_back(1.317046);
    }
    else if (atomnum == 29) {                           //Cu
      exponents.push_back(2.375425);
      exponents.push_back(1.550837);
      exponents.push_back(1.984703);
    }
    else if (atomnum == 30) {                           //Zn
      exponents.push_back(1.664847);
      exponents.push_back(1.176434);
    }
    else if (atomnum == 31) {                           //Ga
      exponents.push_back(1.720919);
      exponents.push_back(1.591570);
      exponents.push_back(1.050000);
    }
    else if (atomnum == 32) {                           //Ge
      exponents.push_back(1.990429);
      exponents.push_back(1.830340);
      exponents.push_back(1.100000);
    }
    else if (atomnum == 33) {                           //As
      exponents.push_back(2.026128);
      exponents.push_back(1.949257);
      exponents.push_back(1.040181);
    }
    else if (atomnum == 34) {                           //Se
      exponents.push_back(2.230969);
      exponents.push_back(2.150656);
      exponents.push_back(1.317549);
    }
    else if (atomnum == 35) {                           //Br
      exponents.push_back(2.077587);
      exponents.push_back(2.263120);
      exponents.push_back(1.845038);
    }
    else if (atomnum == 36) {                           //Kr
      exponents.push_back(2.445680);
      exponents.push_back(2.210494);
      exponents.push_back(1.884991);
    }
    else if (atomnum == 37) {                           //Rb
      exponents.push_back(1.017267);
      exponents.push_back(0.870130);
    }
    else if (atomnum == 38) {                           //Sr
      exponents.push_back(1.419028);
      exponents.push_back(0.928932);
      exponents.push_back(1.500000);
    }
    else if (atomnum == 39) {                           //Y
      exponents.push_back(2.670141);
      exponents.push_back(1.633876);
      exponents.push_back(1.165412);
    }
    else if (atomnum == 40) {                           //Zr
      exponents.push_back(2.238668);
      exponents.push_back(1.702480);
      exponents.push_back(1.129590);
    }
    else if (atomnum == 41) {                           //Nb
      exponents.push_back(1.706832);
      exponents.push_back(1.666463);
      exponents.push_back(1.132172);
    }
    else if (atomnum == 42) {                           //Mo
      exponents.push_back(1.777658);
      exponents.push_back(1.639917);
      exponents.push_back(1.159781);
    }
    else if (atomnum == 43) {                           //Tc
      exponents.push_back(1.918066);
      exponents.push_back(1.918167);
      exponents.push_back(1.346082);
    }
    else if (atomnum == 44) {                           //Ru
      exponents.push_back(2.102697);
      exponents.push_back(1.749643);
      exponents.push_back(1.348322);
    }
    else if (atomnum == 45) {                           //Rh
      exponents.push_back(2.458187);
      exponents.push_back(1.811796);
      exponents.push_back(1.398452);
    }
    else if (atomnum == 46) {                           //Pd
      exponents.push_back(2.353691);
      exponents.push_back(1.828354);
      exponents.push_back(1.333352);
    }
    else if (atomnum == 47) {                           //Ag
      exponents.push_back(2.843549);
      exponents.push_back(1.798462);
      exponents.push_back(1.266649);
    }
    else if (atomnum == 48) {                           //Cd
      exponents.push_back(1.846689);
      exponents.push_back(1.141823);
    }
    else if (atomnum == 49) {                           //In
      exponents.push_back(1.963283);
      exponents.push_back(1.685138);
      exponents.push_back(1.050000);
    }
    else if (atomnum == 50) {                           //Sn
      exponents.push_back(2.551510);
      exponents.push_back(1.893784);
      exponents.push_back(1.100000);
    }
    else if (atomnum == 51) {                           //Sb
      exponents.push_back(2.307407);
      exponents.push_back(2.179752);
      exponents.push_back(1.256087);
    }
    else if (atomnum == 52) {                           //Te
      exponents.push_back(2.434144);
      exponents.push_back(2.182459);
      exponents.push_back(1.373076);
    }
    else if (atomnum == 53) {                           //I
      exponents.push_back(2.159500);
      exponents.push_back(2.308379);
      exponents.push_back(1.691185);
    }
    else if (atomnum == 54) {                           //Xe
      exponents.push_back(2.715140);
      exponents.push_back(2.312510);
      exponents.push_back(1.855707);
    }
    else if (atomnum == 55) {                           //Cs
      exponents.push_back(1.225688);
      exponents.push_back(0.823818);
    }
    else if (atomnum == 56) {                           //Ba
      exponents.push_back(1.528102);
      exponents.push_back(0.991572);
      exponents.push_back(1.500000);
    }
    else if (atomnum == 57) {                           //La
      exponents.push_back(2.875048);
      exponents.push_back(1.731390);
      exponents.push_back(1.303590);
    }
    else if (atomnum == 58) {                           //Ce
      exponents.push_back(2.870000);
      exponents.push_back(1.725197);
      exponents.push_back(1.309804);
    }
    else if (atomnum == 59) {                           //Pr
      exponents.push_back(2.872308);
      exponents.push_back(1.729767);
      exponents.push_back(1.315495);
    }
    else if (atomnum == 60) {                           //Nd
      exponents.push_back(2.874615);
      exponents.push_back(1.734337);
      exponents.push_back(1.321186);
    }
    else if (atomnum == 61) {                           //Pm
      exponents.push_back(2.876923);
      exponents.push_back(1.738907);
      exponents.push_back(1.326877);
    }
    else if (atomnum == 62) {                           //Sm
      exponents.push_back(2.879231);
      exponents.push_back(1.743478);
      exponents.push_back(1.332567);
    }
    else if (atomnum == 63) {                           //Eu
      exponents.push_back(2.881538);
      exponents.push_back(1.748048);
      exponents.push_back(1.338258);
    }
    else if (atomnum == 64) {                           //Gd
      exponents.push_back(2.883846);
      exponents.push_back(1.752618);
      exponents.push_back(1.343949);
    }
    else if (atomnum == 65) {                           //Tb
      exponents.push_back(2.886154);
      exponents.push_back(1.757188);
      exponents.push_back(1.349640);
    }
    else if (atomnum == 66) {                           //Dy
      exponents.push_back(2.888462);
      exponents.push_back(1.761758);
      exponents.push_back(1.355331);
    }
    else if (atomnum == 67) {                           //Ho
      exponents.push_back(2.890769);
      exponents.push_back(1.766328);
      exponents.push_back(1.361022);
    }
    else if (atomnum == 68) {                           //Er
      exponents.push_back(2.893077);
      exponents.push_back(1.770899);
      exponents.push_back(1.366713);
    }
    else if (atomnum == 69) {                           //Tm
      exponents.push_back(2.895385);
      exponents.push_back(1.775469);
      exponents.push_back(1.372403);
    }
    else if (atomnum == 70) {                           //Yb
      exponents.push_back(2.897692);
      exponents.push_back(1.780039);
      exponents.push_back(1.378094);
    }
    else if (atomnum == 71) {                           //Lu
      exponents.push_back(2.900000);
      exponents.push_back(1.784609);
      exponents.push_back(1.383785);
    }
    else if (atomnum == 72) {                           //Hf
      exponents.push_back(2.638729);
      exponents.push_back(2.194333);
      exponents.push_back(1.427467);
    }
    else if (atomnum == 73) {                           //Ta
      exponents.push_back(2.018969);
      exponents.push_back(1.996498);
      exponents.push_back(1.407714);
    }
    else if (atomnum == 74) {                           //W
      exponents.push_back(2.155885);
      exponents.push_back(1.892022);
      exponents.push_back(1.458186);
    }
    else if (atomnum == 75) {                           //Re
      exponents.push_back(2.262783);
      exponents.push_back(2.187549);
      exponents.push_back(1.636996);
    }
    else if (atomnum == 76) {                           //Os
      exponents.push_back(2.509631);
      exponents.push_back(2.173991);
      exponents.push_back(1.597888);
    }
    else if (atomnum == 77) {                           //Ir
      exponents.push_back(2.756134);
      exponents.push_back(2.117548);
      exponents.push_back(1.680343);
    }
    else if (atomnum == 78) {                           //Pt
      exponents.push_back(2.704492);
      exponents.push_back(2.329136);
      exponents.push_back(1.623286);
    }
    else if (atomnum == 79) {                           //Au
      exponents.push_back(3.241287);
      exponents.push_back(2.183171);
      exponents.push_back(2.084484);
    }
    else if (atomnum == 80) {                           //Hg
      exponents.push_back(2.244504);
      exponents.push_back(1.470848);
    }
    else if (atomnum == 81) {                           //Tl
      exponents.push_back(2.294231);
      exponents.push_back(1.731592);
    }
    else if (atomnum == 82) {                           //Pb
      exponents.push_back(2.960592);
      exponents.push_back(1.953130);
    }
    else if (atomnum == 83) {                           //Bi
      exponents.push_back(2.788267);
      exponents.push_back(2.277039);
    }
    else if (atomnum == 84) {                           //Po
      exponents.push_back(3.314810);
      exponents.push_back(2.389456);
    }
    else if (atomnum == 85) {                           //At
      exponents.push_back(2.220421);
      exponents.push_back(2.408112);
      exponents.push_back(1.500000);
    }
    else if (atomnum == 86) {                           //Rn
      exponents.push_back(3.109394);
      exponents.push_back(2.541934);
      exponents.push_back(1.790000);
    }
    else if (atomnum > 86) {                            //unavailable: Fr,Ra,Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr,Rf,Db,Sg,Bh,Hs,Mt,Ds,Rg,Cn,Nh,Fl,Mc,Lv,Ts,Og
      throw("ERROR: basissets/slater.hpp: SlaterBSet: GetExponents(): orbital exponents not defined for atom");
    }
    return exponents;
  }
  matrix<size_t> QuantumNumbers(size_t atomnum) {
    //getter for quantum numbers of one atom
    std::vector<double> exps = GetExponents(atomnum);
    matrix<size_t> quantumnumb(exps.size(),2);
    if (atomnum > 86) {
      //unavailable: Fr,Ra,Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr,Rf,Db,Sg,Bh,Hs,Mt,Ds,Rg,Cn,Nh,Fl,Mc,Lv,Ts,Og
      throw("ERROR: basissets/dewar.hpp: PM6BSet: QuantumNumbers(): atom not defined");
    }
    if (atomnum == 1) {quantumnumb(1,1) = 1;}                                 //H
    else if (atomnum == 2) {                                                          //He
      quantumnumb(1,1) = 1;    //1s
      quantumnumb(2,1) = 2;    //2p
      quantumnumb(2,2) = 1;    //2p
    }
    else if ((atomnum > 2)&&(atomnum < 10)) {                                         //Li,Be,B,C,N,O,F
      quantumnumb(1,1) = 2;    //2s
      quantumnumb(2,1) = 2;    //2p
      quantumnumb(2,2) = 1;    //2p
    }
    else if (atomnum == 10) {                                                         //Ne
      quantumnumb(1,1) = 2;    //2s
      quantumnumb(2,1) = 2;    //2p
      quantumnumb(2,2) = 1;    //2p
      quantumnumb(3,1) = 3;    //3d
      quantumnumb(3,2) = 2;    //3d
    }
    else if ((atomnum > 10)&&(atomnum < 19)) {                                        //Na,Mg,Al,Si,P,S,Cl,Ar
      quantumnumb(1,1) = 3;    //3s
      quantumnumb(2,1) = 3;    //3p
      quantumnumb(2,2) = 1;    //3p
      if (atomnum > 11) {
        quantumnumb(3,1) = 3;  //3d
        quantumnumb(3,2) = 2;  //3d
      }
    }
    else if (((atomnum > 18)&&(atomnum < 21))||((atomnum > 29)&&(atomnum < 37))) {    //K,Ca,Zn,Ga,Ge,As,Se,Br
      quantumnumb(1,1) = 4;    //4s
      quantumnumb(2,1) = 4;    //4p
      quantumnumb(2,2) = 1;    //4p
      if (atomnum == 20) {
        quantumnumb(3,1) = 3;  //3d
        quantumnumb(3,2) = 2;  //3d
      }
      else if (atomnum > 30) {
        quantumnumb(3,1) = 4;  //4d
        quantumnumb(3,2) = 2;  //4d
      }
    }
    else if ((atomnum > 20)&&(atomnum < 30)) {                                        //Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu
      quantumnumb(1,1) = 3;    //3d
      quantumnumb(1,2) = 2;    //3d
      quantumnumb(2,1) = 4;    //4s
      quantumnumb(3,1) = 4;    //4p
      quantumnumb(3,2) = 1;    //4p
    }
    else if (((atomnum > 36)&&(atomnum < 39))||((atomnum > 47)&&(atomnum < 55))) {    //Rb,Sr,Cd,In,Sn,Sb,Te,I
      quantumnumb(1,1) = 5;    //5s
      quantumnumb(2,1) = 5;    //5p
      quantumnumb(2,2) = 1;    //5p
      if (atomnum == 38) {
        quantumnumb(3,1) = 4;  //4d
        quantumnumb(3,2) = 2;  //4d
      }
      else if (atomnum > 48) {
        quantumnumb(3,1) = 5;  //5d
        quantumnumb(3,2) = 2;  //5d
      }
    }
    else if ((atomnum > 38)&&(atomnum < 48)) {                                        //Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag
      quantumnumb(1,1) = 4;    //4d
      quantumnumb(1,2) = 2;    //4d
      quantumnumb(2,1) = 5;    //5s
      quantumnumb(3,1) = 5;    //5p
      quantumnumb(3,2) = 1;    //5p
    }
    else if (((atomnum > 54)&&(atomnum < 57))||(atomnum > 79)) {                      //Cs,Ba,Hg,Tl,Pb,Bi,Po,At,Rn
      quantumnumb(1,1) = 6;    //6s
      quantumnumb(2,1) = 6;    //6p
      quantumnumb(2,2) = 1;    //6p
      if ((atomnum == 56)||(atomnum > 84)) {
        quantumnumb(3,1) = 5;  //5d
        quantumnumb(3,2) = 2;  //5d
      }
    }
    else if ((atomnum > 56)&&(atomnum < 80)) {                                        //La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au
      quantumnumb(1,1) = 5;    //5d
      quantumnumb(1,2) = 2;    //5d
      quantumnumb(2,1) = 6;    //6s
      quantumnumb(3,1) = 6;    //6p
      quantumnumb(3,2) = 1;    //6p
    }
    return quantumnumb;
  }
  int getNumberPrimitives(size_t atomnr, int lnumb) {
    //function determining the number of Gaussian primitives used for a given element, orbital with certain azimuthal quantum number
    int nprimitives = 0;
    if (lnumb == 0) {               //s orbitals
      nprimitives = 4;
      if (atomnr < 3) {nprimitives = 3;}
      else if (atomnr > 54) {nprimitives = 6;}
    }
    else if (lnumb == 1) {          //p orbitals
      nprimitives = 4;
      if (atomnr > 54) {nprimitives = 6;}
    }
    else if (lnumb == 2) {          //d orbitals
      nprimitives = 3;
    }
    else if (lnumb == 3) {          //f orbitals
      nprimitives = 4;
    }
    return nprimitives;
  }
  size_t GetNAOs(size_t atomnum) {
    //function to calculate the number of AOs used for a given atom
    std::vector<double> exps = GetExponents(atomnum);
    return exps.size()*exps.size();
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapGFN2(atm);}
};

#endif //_GFN2_Basis_Set_
