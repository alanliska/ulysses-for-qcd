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

#ifndef _GFN2_XTB_Parameters_
#define _GFN2_XTB_Parameters_
#include "AtomicRadiipar.hpp"

//description:
//parameters for GFN2-xTB

double alphaMultipole(int multi) {
  //function returning the exponent for the damping of dipole (multi = 1) or quadrupole (multi = 2) moments
  double alpha = 0.0;
  if (multi == 1) {alpha = 3.0;}
  else if (multi == 2) {alpha = 4.0;}
  return alpha;
}
double kval(int LA, int LB) {
  //function returning the k values for orbital pairs
  int Lmin = std::min(LA,LB);
  int Lmax = std::max(LA,LB);
  double kk = 1.85;                                     //ss case
  if ((Lmin == Lmax)&&(Lmax > 0)) {kk = 2.23;}          //pp, dd cases
  else if ((Lmin == 0)&&(Lmax == 1)) {kk = 2.04;}       //sp case
  else if ((Lmin < 2)&&(Lmax == 2)) {kk = 2.0;}         //sd, pd cases
  return kk;
}
double Korb(int Lorb) {
  //function returning the K values for orbitals
  double KK = 1.0;
  if (Lorb > 0) {KK /= double(2*Lorb);}
  return KK;
}
double kEN() {
  //function returning the kEN value
  double ken = 0.02;
  return ken;
}
double EnergyLevel(int atm, int Lorb) {
  //energy levels for GFN2-xTB
  //values stored in eV, returned in Eh
  double enlev = 0.0;
  switch (atm) {
    case 1:                               //H
      enlev = -10.707211;
      break;
    case 2:                               //He
      if (Lorb == 0) {enlev = -23.716445;}
      else if (Lorb == 1) {enlev = -1.822307;}
      break;
    case 3:                               //Li
      if (Lorb == 0) {enlev = -4.900000;}
      else if (Lorb == 1) {enlev = -2.217789;}
      break;
    case 4:                               //Be
      if (Lorb == 0) {enlev = -7.743081;}
      else if (Lorb == 1) {enlev = -3.133433;}
      break;
    case 5:                               //B
      if (Lorb == 0) {enlev = -9.224376;}
      else if (Lorb == 1) {enlev = -7.419002;}
      break;
    case 6:                               //C
      if (Lorb == 0) {enlev = -13.970922;}
      else if (Lorb == 1) {enlev = -10.063292;}
      break;
    case 7:                               //N
      if (Lorb == 0) {enlev = -16.686243;}
      else if (Lorb == 1) {enlev = -12.523956;}
      break;
    case 8:                               //O
      if (Lorb == 0) {enlev = -20.229985;}
      else if (Lorb == 1) {enlev = -15.503117;}
      break;
    case 9:                               //F
      if (Lorb == 0) {enlev = -23.458179;}
      else if (Lorb == 1) {enlev = -15.746583;}
      break;
    case 10:                              //Ne
      if (Lorb == 0) {enlev = -24.500000;}
      else if (Lorb == 1) {enlev = -18.737298;}
      else if (Lorb == 2) {enlev = -5.517827;}
      break;
    case 11:                              //Na
      if (Lorb == 0) {enlev = -4.546934;}
      else if (Lorb == 1) {enlev = -1.332719;}
      break;
    case 12:                              //Mg
      if (Lorb == 0) {enlev = -6.339908;}
      else if (Lorb == 1) {enlev = -0.697688;}
      else if (Lorb == 2) {enlev = -1.458197;}
      break;
    case 13:                              //Al
      if (Lorb == 0) {enlev = -9.329017;}
      else if (Lorb == 1) {enlev = -5.927846;}
      else if (Lorb == 2) {enlev = -3.042325;}
      break;
    case 14:                              //Si
      if (Lorb == 0) {enlev = -14.360932;}
      else if (Lorb == 1) {enlev = -6.915131;}
      else if (Lorb == 2) {enlev = -1.825036;}
      break;
    case 15:                              //P
      if (Lorb == 0) {enlev = -17.518756;}
      else if (Lorb == 1) {enlev = -9.842286;}
      else if (Lorb == 2) {enlev = -0.444893;}
      break;
    case 16:                              //S
      if (Lorb == 0) {enlev = -20.029654;}
      else if (Lorb == 1) {enlev = -11.377694;}
      else if (Lorb == 2) {enlev = -0.420282;}
      break;
    case 17:                              //Cl
      if (Lorb == 0) {enlev = -29.278781;}
      else if (Lorb == 1) {enlev = -12.673758;}
      else if (Lorb == 2) {enlev = -0.240338;}
      break;
    case 18:                              //Ar
      if (Lorb == 0) {enlev = -16.487730;}
      else if (Lorb == 1) {enlev = -13.910539;}
      else if (Lorb == 2) {enlev = -1.167213;}
      break;
    case 19:                              //K
      if (Lorb == 0) {enlev = -4.510348;}
      else if (Lorb == 1) {enlev = -0.934377;}
      break;
    case 20:                              //Ca
      if (Lorb == 0) {enlev = -5.056506;}
      else if (Lorb == 1) {enlev = -1.150304;}
      else if (Lorb == 2) {enlev = -0.776883;}
      break;
    case 21:                              //Sc
      if (Lorb == 2) {enlev = -5.196187;}
      else if (Lorb == 0) {enlev = -8.877940;}
      else if (Lorb == 1) {enlev = -2.008206;}
      break;
    case 22:                              //Ti
      if (Lorb == 2) {enlev = -7.234331;}
      else if (Lorb == 0) {enlev = -10.900000;}
      else if (Lorb == 1) {enlev = -1.928783;}
      break;
    case 23:                              //V
      if (Lorb == 2) {enlev = -9.015342;}
      else if (Lorb == 0) {enlev = -9.573347;}
      else if (Lorb == 1) {enlev = -0.706647;}
      break;
    case 24:                              //Cr
      if (Lorb == 2) {enlev = -7.209794;}
      else if (Lorb == 0) {enlev = -9.201304;}
      else if (Lorb == 1) {enlev = -0.696957;}
      break;
    case 25:                              //Mn
      if (Lorb == 2) {enlev = -10.120933;}
      else if (Lorb == 0) {enlev = -5.617346;}
      else if (Lorb == 1) {enlev = -4.198724;}
      break;
    case 26:                              //Fe
      if (Lorb == 2) {enlev = -10.035473;}
      else if (Lorb == 0) {enlev = -5.402911;}
      else if (Lorb == 1) {enlev = -3.308988;}
      break;
    case 27:                              //Co
      if (Lorb == 2) {enlev = -10.580430;}
      else if (Lorb == 0) {enlev = -8.596723;}
      else if (Lorb == 1) {enlev = -2.585753;}
      break;
    case 28:                              //Ni
      if (Lorb == 2) {enlev = -12.712236;}
      else if (Lorb == 0) {enlev = -8.524281;}
      else if (Lorb == 1) {enlev = -2.878873;}
      break;
    case 29:                              //Cu
      if (Lorb == 2) {enlev = -9.506548;}
      else if (Lorb == 0) {enlev = -6.922958;}
      else if (Lorb == 1) {enlev = -2.267723;}
      break;
    case 30:                              //Zn
      if (Lorb == 0) {enlev = -7.177294;}
      else if (Lorb == 1) {enlev = -0.991895;}
      break;
    case 31:                              //Ga
      if (Lorb == 0) {enlev = -12.449656;}
      else if (Lorb == 1) {enlev = -4.469873;}
      else if (Lorb == 2) {enlev = -0.582255;}
      break;
    case 32:                              //Ge
      if (Lorb == 0) {enlev = -16.369792;}
      else if (Lorb == 1) {enlev = -8.207673;}
      else if (Lorb == 2) {enlev = -0.994226;}
      break;
    case 33:                              //As
      if (Lorb == 0) {enlev = -16.421504;}
      else if (Lorb == 1) {enlev = -9.311147;}
      else if (Lorb == 2) {enlev = -0.276830;}
      break;
    case 34:                              //Se
      if (Lorb == 0) {enlev = -20.584732;}
      else if (Lorb == 1) {enlev = -10.910799;}
      else if (Lorb == 2) {enlev = -0.110636;}
      break;
    case 35:                              //Br
      if (Lorb == 0) {enlev = -23.583718;}
      else if (Lorb == 1) {enlev = -12.588824;}
      else if (Lorb == 2) {enlev = 0.047980;}
      break;
    case 36:                              //Kr
      if (Lorb == 0) {enlev = -17.221422;}
      else if (Lorb == 1) {enlev = -13.633377;}
      else if (Lorb == 2) {enlev = -0.940657;}
      break;
    case 37:                              //Rb
      if (Lorb == 0) {enlev = -4.353793;}
      else if (Lorb == 1) {enlev = -1.392938;}
      break;
    case 38:                              //Sr
      if (Lorb == 0) {enlev = -6.291692;}
      else if (Lorb == 1) {enlev = -1.872475;}
      else if (Lorb == 2) {enlev = -0.890492;}
      break;
    case 39:                              //Y
      if (Lorb == 2) {enlev = -8.015206;}
      else if (Lorb == 0) {enlev = -12.194181;}
      else if (Lorb == 1) {enlev = -0.966195;}
      break;
    case 40:                              //Zr
      if (Lorb == 2) {enlev = -7.409832;}
      else if (Lorb == 0) {enlev = -10.199105;}
      else if (Lorb == 1) {enlev = -1.066939;}
      break;
    case 41:                              //Nb
      if (Lorb == 2) {enlev = -8.440821;}
      else if (Lorb == 0) {enlev = -11.384021;}
      else if (Lorb == 1) {enlev = -0.103760;}
      break;
    case 42:                              //Mo
      if (Lorb == 2) {enlev = -7.995133;}
      else if (Lorb == 0) {enlev = -7.336245;}
      else if (Lorb == 1) {enlev = -3.686225;}
      break;
    case 43:                              //Tc
      if (Lorb == 2) {enlev = -9.587897;}
      else if (Lorb == 0) {enlev = -6.792444;}
      else if (Lorb == 1) {enlev = -3.325525;}
      break;
    case 44:                              //Ru
      if (Lorb == 2) {enlev = -10.285405;}
      else if (Lorb == 0) {enlev = -5.332608;}
      else if (Lorb == 1) {enlev = -3.307153;}
      break;
    case 45:                              //Rh
      if (Lorb == 2) {enlev = -11.756644;}
      else if (Lorb == 0) {enlev = -7.850495;}
      else if (Lorb == 1) {enlev = -3.007906;}
      break;
    case 46:                              //Pd
      if (Lorb == 2) {enlev = -11.963518;}
      else if (Lorb == 0) {enlev = -9.714059;}
      else if (Lorb == 1) {enlev = -2.035281;}
      break;
    case 47:                              //Ag
      if (Lorb == 2) {enlev = -9.591083;}
      else if (Lorb == 0) {enlev = -8.083960;}
      else if (Lorb == 1) {enlev = -2.934333;}
      break;
    case 48:                              //Cd
      if (Lorb == 0) {enlev = -7.252341;}
      else if (Lorb == 1) {enlev = -0.744865;}
      break;
    case 49:                              //In
      if (Lorb == 0) {enlev = -13.040909;}
      else if (Lorb == 1) {enlev = -4.507143;}
      else if (Lorb == 2) {enlev = -0.805666;}
      break;
    case 50:                              //Sn
      if (Lorb == 0) {enlev = -19.970428;}
      else if (Lorb == 1) {enlev = -7.367059;}
      else if (Lorb == 2) {enlev = -2.077548;}
      break;
    case 51:                              //Sb
      if (Lorb == 0) {enlev = -18.371244;}
      else if (Lorb == 1) {enlev = -7.350148;}
      else if (Lorb == 2) {enlev = 0.909033;}
      break;
    case 52:                              //Te
      if (Lorb == 0) {enlev = -21.930653;}
      else if (Lorb == 1) {enlev = -9.480374;}
      else if (Lorb == 2) {enlev = 0.978922;}
      break;
    case 53:                              //I
      if (Lorb == 0) {enlev = -20.949407;}
      else if (Lorb == 1) {enlev = -12.180159;}
      else if (Lorb == 2) {enlev = -0.266596;}
      break;
    case 54:                              //Xe
      if (Lorb == 0) {enlev = -19.090498;}
      else if (Lorb == 1) {enlev = -11.249471;}
      else if (Lorb == 2) {enlev = -0.497097;}
      break;
    case 55:                              //Cs
      if (Lorb == 0) {enlev = -4.041706;}
      else if (Lorb == 1) {enlev = -1.394193;}
      break;
    case 56:                              //Ba
      if (Lorb == 0) {enlev = -5.900000;}
      else if (Lorb == 1) {enlev = -2.133395;}
      else if (Lorb == 2) {enlev = -1.514900;}
      break;
    case 57:                              //La
      if (Lorb == 2) {enlev = -8.958783;}
      else if (Lorb == 0) {enlev = -11.877410;}
      else if (Lorb == 1) {enlev = -0.601717;}
      break;
    case 58:                              //Ce
      if (Lorb == 2) {enlev = -7.381991;}
      else if (Lorb == 0) {enlev = -8.537781;}
      else if (Lorb == 1) {enlev = -3.017508;}
      break;
    case 59:                              //Pr
      if (Lorb == 2) {enlev = -7.280875;}
      else if (Lorb == 0) {enlev = -8.504806;}
      else if (Lorb == 1) {enlev = -2.873159;}
      break;
    case 60:                              //Nd
      if (Lorb == 2) {enlev = -7.179760;}
      else if (Lorb == 0) {enlev = -8.471830;}
      else if (Lorb == 1) {enlev = -2.728809;}
      break;
    case 61:                              //Pm
      if (Lorb == 2) {enlev = -7.078644;}
      else if (Lorb == 0) {enlev = -8.438855;}
      else if (Lorb == 1) {enlev = -2.584460;}
      break;
    case 62:                              //Sm
      if (Lorb == 2) {enlev = -6.977529;}
      else if (Lorb == 0) {enlev = -8.405879;}
      else if (Lorb == 1) {enlev = -2.440110;}
      break;
    case 63:                              //Eu
      if (Lorb == 2) {enlev = -6.876413;}
      else if (Lorb == 0) {enlev = -8.372904;}
      else if (Lorb == 1) {enlev = -2.295761;}
      break;
    case 64:                              //Gd
      if (Lorb == 2) {enlev = -6.775298;}
      else if (Lorb == 0) {enlev = -8.339929;}
      else if (Lorb == 1) {enlev = -2.151411;}
      break;
    case 65:                              //Tb
      if (Lorb == 2) {enlev = -6.674182;}
      else if (Lorb == 0) {enlev = -8.306953;}
      else if (Lorb == 1) {enlev = -2.007062;}
      break;
    case 66:                              //Dy
      if (Lorb == 2) {enlev = -6.573067;}
      else if (Lorb == 0) {enlev = -8.273978;}
      else if (Lorb == 1) {enlev = -1.862712;}
      break;
    case 67:                              //Ho
      if (Lorb == 2) {enlev = -6.471951;}
      else if (Lorb == 0) {enlev = -8.241003;}
      else if (Lorb == 1) {enlev = -1.718363;}
      break;
    case 68:                              //Er
      if (Lorb == 2) {enlev = -6.370836;}
      else if (Lorb == 0) {enlev = -8.208027;}
      else if (Lorb == 1) {enlev = -1.574013;}
      break;
    case 69:                              //Tm
      if (Lorb == 2) {enlev = -6.269720;}
      else if (Lorb == 0) {enlev = -8.175052;}
      else if (Lorb == 1) {enlev = -1.429664;}
      break;
    case 70:                              //Yb
      if (Lorb == 2) {enlev = -6.168604;}
      else if (Lorb == 0) {enlev = -8.142076;}
      else if (Lorb == 1) {enlev = -1.285314;}
      break;
    case 71:                              //Lu
      if (Lorb == 2) {enlev = -6.067489;}
      else if (Lorb == 0) {enlev = -8.109101;}
      else if (Lorb == 1) {enlev = -1.140965;}
      break;
    case 72:                              //Hf
      if (Lorb == 2) {enlev = -7.181755;}
      else if (Lorb == 0) {enlev = -10.626891;}
      else if (Lorb == 1) {enlev = -1.603430;}
      break;
    case 73:                              //Ta
      if (Lorb == 2) {enlev = -8.481353;}
      else if (Lorb == 0) {enlev = -13.073088;}
      else if (Lorb == 1) {enlev = 0.655254;}
      break;
    case 74:                              //W
      if (Lorb == 2) {enlev = -9.501505;}
      else if (Lorb == 0) {enlev = -11.093016;}
      else if (Lorb == 1) {enlev = -1.420389;}
      break;
    case 75:                              //Re
      if (Lorb == 2) {enlev = -11.189119;}
      else if (Lorb == 0) {enlev = -12.685198;}
      else if (Lorb == 1) {enlev = -3.851981;}
      break;
    case 76:                              //Os
      if (Lorb == 2) {enlev = -10.382841;}
      else if (Lorb == 0) {enlev = -8.731460;}
      else if (Lorb == 1) {enlev = -3.546379;}
      break;
    case 77:                              //Ir
      if (Lorb == 2) {enlev = -11.018475;}
      else if (Lorb == 0) {enlev = -9.349164;}
      else if (Lorb == 1) {enlev = -3.603762;}
      break;
    case 78:                              //Pt
      if (Lorb == 2) {enlev = -12.047728;}
      else if (Lorb == 0) {enlev = -10.482306;}
      else if (Lorb == 1) {enlev = -3.778297;}
      break;
    case 79:                              //Au
      if (Lorb == 2) {enlev = -9.578599;}
      else if (Lorb == 0) {enlev = -7.688552;}
      else if (Lorb == 1) {enlev = 0.883399;}
      break;
    case 80:                              //Hg
      if (Lorb == 0) {enlev = -11.538066;}
      else if (Lorb == 1) {enlev = -2.532581;}
      break;
    case 81:                              //Tl
      if (Lorb == 0) {enlev = -17.319333;}
      else if (Lorb == 1) {enlev = -4.460584;}
      break;
    case 82:                              //Pb
      if (Lorb == 0) {enlev = -24.055207;}
      else if (Lorb == 1) {enlev = -5.893816;}
      break;
    case 83:                              //Bi
      if (Lorb == 0) {enlev = -19.843840;}
      else if (Lorb == 1) {enlev = -7.297456;}
      break;
    case 84:                              //Po
      if (Lorb == 0) {enlev = -20.205380;}
      else if (Lorb == 1) {enlev = -8.476927;}
      break;
    case 85:                              //At
      if (Lorb == 0) {enlev = -17.050229;}
      else if (Lorb == 1) {enlev = -9.499822;}
      else if (Lorb == 2) {enlev = -0.096063;}
      break;
    case 86:                              //Rn
      if (Lorb == 0) {enlev = -21.000000;}
      else if (Lorb == 1) {enlev = -10.496406;}
      else if (Lorb == 2) {enlev = -1.415056;}
      break;
  }
  return enlev/au2eV;
}
double HCNl(int atm, int Lorb) {
  //HCN(l) values required for GFN2-xTB
  //values stored in eV, returned in Eh
  double hcn = 0.0;
  switch (atm) {
    case 1:                                //H
      hcn = -0.500000;
      break;
    case 2:                                //He
      if (Lorb == 0) {hcn = 2.074275;}
      break;
    case 3:                                //Li
      if (Lorb == 0) {hcn = 1.620836;}
      else if (Lorb == 1) {hcn = -0.623876;}
      break;
    case 4:                                //Be
      if (Lorb == 0) {hcn = 1.187759;}
      else if (Lorb == 1) {hcn = 0.550528;}
      break;
    case 5:                                //B
      if (Lorb == 0) {hcn = 0.120462;}
      else if (Lorb == 1) {hcn = -0.141086;}
      break;
    case 6:                                //C
      if (Lorb == 0) {hcn = -0.102144;}
      else if (Lorb == 1) {hcn = 0.161657;}
      break;
    case 7:                                //N
      if (Lorb == 0) {hcn = -1.955336;}
      else if (Lorb == 1) {hcn = 0.561076;}
      break;
    case 8:                                //O
      if (Lorb == 0) {hcn = 0.117826;}
      else if (Lorb == 1) {hcn = -0.145102;}
      break;
    case 9:                                //F
      if (Lorb == 0) {hcn = 0.394362;}
      else if (Lorb == 1) {hcn = -0.538373;}
      break;
    case 10:                               //Ne
      if (Lorb == 0) {hcn = -0.014933;}
      else if (Lorb == 1) {hcn = 0.232093;}
      else if (Lorb == 2) {hcn = 1.096714;}
      break;
    case 11:                               //Na
      if (Lorb == 0) {hcn = -0.042211;}
      else if (Lorb == 1) {hcn = -0.144323;}
      break;
    case 12:                               //Mg
      if (Lorb == 0) {hcn = 1.164444;}
      else if (Lorb == 1) {hcn = -0.079924;}
      else if (Lorb == 2) {hcn = 1.192409;}
      break;
    case 13:                               //Al
      if (Lorb == 0) {hcn = 0.715422;}
      else if (Lorb == 1) {hcn = -0.244485;}
      else if (Lorb == 2) {hcn = 0.406173;}
      break;
    case 14:                               //Si
      if (Lorb == 0) {hcn = 1.858479;}
      else if (Lorb == 1) {hcn = -1.383073;}
      else if (Lorb == 2) {hcn = -1.935494;}
      break;
    case 15:                               //P
      if (Lorb == 0) {hcn = 0.547610;}
      else if (Lorb == 1) {hcn = -0.489930;}
      else if (Lorb == 2) {hcn = 2.429507;}
      break;
    case 16:                               //S
      if (Lorb == 0) {hcn = -0.256951;}
      else if (Lorb == 1) {hcn = -0.098465;}
      else if (Lorb == 2) {hcn = 2.007690;}
      break;
    case 17:                               //Cl
      if (Lorb == 0) {hcn = 0.617972;}
      else if (Lorb == 1) {hcn = -0.181618;}
      else if (Lorb == 2) {hcn = 1.672768;}
      break;
    case 18:                               //Ar
      if (Lorb == 0) {hcn = 0.000554;}
      else if (Lorb == 1) {hcn = 0.065921;}
      else if (Lorb == 2) {hcn = -2.732165;}
      break;
    case 19:                               //K
      if (Lorb == 0) {hcn = -0.339245;}
      else if (Lorb == 1) {hcn = 0.174542;}
      break;
    case 20:                               //Ca
      if (Lorb == 0) {hcn = 0.570930;}
      else if (Lorb == 1) {hcn = -0.074926;}
      else if (Lorb == 2) {hcn = 1.013752;}
      break;
    case 21:                               //Sc
      if (Lorb == 0) {hcn = 0.991293;}
      else if (Lorb == 1) {hcn = -0.281241;}
      else if (Lorb == 2) {hcn = 2.026777;}
      break;
    case 22:                               //Ti
      if (Lorb == 0) {hcn = 1.007021;}
      else if (Lorb == 1) {hcn = -0.237074;}
      else if (Lorb == 2) {hcn = 1.028188;}
      break;
    case 23:                               //V
      if (Lorb == 0) {hcn = 0.235696;}
      else if (Lorb == 1) {hcn = -0.108232;}
      else if (Lorb == 2) {hcn = 0.164476;}
      break;
    case 24:                               //Cr
      if (Lorb == 0) {hcn = -0.232087;}
      else if (Lorb == 1) {hcn = -0.188919;}
      else if (Lorb == 2) {hcn = 0.289291;}
      break;
    case 25:                               //Mn
      if (Lorb == 0) {hcn = -0.275000;}
      else if (Lorb == 1) {hcn = -0.015839;}
      else if (Lorb == 2) {hcn = -0.195827;}
      break;
    case 26:                               //Fe
      if (Lorb == 0) {hcn = -4.049876;}
      else if (Lorb == 1) {hcn = -0.756480;}
      else if (Lorb == 2) {hcn = -0.274654;}
      break;
    case 27:                               //Co
      if (Lorb == 0) {hcn = -0.227872;}
      else if (Lorb == 1) {hcn = 0.076513;}
      else if (Lorb == 2) {hcn = 0.121980;}
      break;
    case 28:                               //Ni
      if (Lorb == 0) {hcn = 0.310301;}
      else if (Lorb == 1) {hcn = 0.226796;}
      else if (Lorb == 2) {hcn = -0.066417;}
      break;
    case 29:                               //Cu
      if (Lorb == 0) {hcn = 3.349047;}
      else if (Lorb == 1) {hcn = -2.619446;}
      else if (Lorb == 2) {hcn = -0.173684;}
      break;
    case 30:                               //Zn
      if (Lorb == 0) {hcn = 2.011910;}
      else if (Lorb == 1) {hcn = -0.055135;}
      break;
    case 31:                               //Ga
      if (Lorb == 0) {hcn = -0.234627;}
      else if (Lorb == 1) {hcn = 1.305831;}
      else if (Lorb == 2) {hcn = 0.165604;}
      break;
    case 32:                               //Ge
      if (Lorb == 0) {hcn = 0.361068;}
      else if (Lorb == 1) {hcn = -0.014474;}
      else if (Lorb == 2) {hcn = -1.042559;}
      break;
    case 33:                               //As
      if (Lorb == 0) {hcn = -0.129640;}
      else if (Lorb == 1) {hcn = -0.236470;}
      else if (Lorb == 2) {hcn = 2.330144;}
      break;
    case 34:                               //Se
      if (Lorb == 0) {hcn = -0.061654;}
      else if (Lorb == 1) {hcn = -0.435018;}
      else if (Lorb == 2) {hcn = 2.768559;}
      break;
    case 35:                               //Br
      if (Lorb == 0) {hcn = 0.006150;}
      else if (Lorb == 1) {hcn = -0.058347;}
      else if (Lorb == 2) {hcn = 2.250180;}
      break;
    case 36:                               //Kr
      if (Lorb == 0) {hcn = -0.070305;}
      else if (Lorb == 1) {hcn = 0.076023;}
      else if (Lorb == 2) {hcn = 0.349523;}
      break;
    case 37:                               //Rb
      if (Lorb == 0) {hcn = -1.516926;}
      else if (Lorb == 1) {hcn = 0.203437;}
      break;
    case 38:                               //Sr
      if (Lorb == 0) {hcn = 0.409020;}
      else if (Lorb == 1) {hcn = -0.418725;}
      else if (Lorb == 2) {hcn = 0.401255;}
      break;
    case 39:                               //Y
      if (Lorb == 0) {hcn = 1.937516;}
      else if (Lorb == 1) {hcn = -0.641897;}
      else if (Lorb == 2) {hcn = -1.270337;}
      break;
    case 40:                               //Zr
      if (Lorb == 0) {hcn = 1.266548;}
      else if (Lorb == 1) {hcn = 0.279435;}
      else if (Lorb == 2) {hcn = -0.566943;}
      break;
    case 41:                               //Nb
      if (Lorb == 0) {hcn = 2.555962;}
      else if (Lorb == 1) {hcn = -0.002341;}
      else if (Lorb == 2) {hcn = -1.356485;}
      break;
    case 42:                               //Mo
      if (Lorb == 0) {hcn = 3.008411;}
      else if (Lorb == 1) {hcn = -1.040348;}
      else if (Lorb == 2) {hcn = 0.620172;}
      break;
    case 43:                               //Tc
      if (Lorb == 0) {hcn = -0.586205;}
      else if (Lorb == 1) {hcn = -0.087319;}
      else if (Lorb == 2) {hcn = -0.066526;}
      break;
    case 44:                               //Ru
      if (Lorb == 0) {hcn = 4.471162;}
      else if (Lorb == 1) {hcn = -0.034723;}
      else if (Lorb == 2) {hcn = -0.263914;}
      break;
    case 45:                               //Rh
      if (Lorb == 0) {hcn = 0.066741;}
      else if (Lorb == 1) {hcn = -0.213308;}
      else if (Lorb == 2) {hcn = 0.104368;}
      break;
    case 46:                               //Pd
      if (Lorb == 0) {hcn = 0.266820;}
      else if (Lorb == 1) {hcn = 0.503075;}
      else if (Lorb == 2) {hcn = 0.060285;}
      break;
    case 47:                               //Ag
      if (Lorb == 0) {hcn = -0.065794;}
      else if (Lorb == 1) {hcn = 1.677171;}
      else if (Lorb == 2) {hcn = -0.062719;}
      break;
    case 48:                               //Cd
      if (Lorb == 0) {hcn = 1.418145;}
      else if (Lorb == 1) {hcn = -0.309814;}
      break;
    case 49:                               //In
      if (Lorb == 0) {hcn = -0.098312;}
      else if (Lorb == 1) {hcn = 0.994688;}
      else if (Lorb == 2) {hcn = 0.168649;}
      break;
    case 50:                               //Sn
      if (Lorb == 0) {hcn = -0.454629;}
      else if (Lorb == 1) {hcn = -0.320651;}
      else if (Lorb == 2) {hcn = -1.459407;}
      break;
    case 51:                               //Sb
      if (Lorb == 0) {hcn = -0.147626;}
      else if (Lorb == 1) {hcn = -0.091175;}
      else if (Lorb == 2) {hcn = 1.602872;}
      break;
    case 52:                               //Te
      if (Lorb == 0) {hcn = 0.115389;}
      else if (Lorb == 1) {hcn = -0.082051;}
      else if (Lorb == 2) {hcn = 3.013228;}
      break;
    case 53:                               //I
      if (Lorb == 0) {hcn = -0.506150;}
      else if (Lorb == 1) {hcn = 0.084766;}
      else if (Lorb == 2) {hcn = 3.077127;}
      break;
    case 54:                               //Xe
      if (Lorb == 0) {hcn = -0.020195;}
      else if (Lorb == 1) {hcn = 0.017246;}
      else if (Lorb == 2) {hcn = 0.327039;}
      break;
    case 55:                               //Cs
      if (Lorb == 0) {hcn = -1.312604;}
      else if (Lorb == 1) {hcn = -0.100000;}
      break;
    case 56:                               //Ba
      if (Lorb == 0) {hcn = 0.352001;}
      else if (Lorb == 1) {hcn = -0.926576;}
      else if (Lorb == 2) {hcn = 0.147995;}
      break;
    case 57:                               //La
      if (Lorb == 0) {hcn = 1.071684;}
      else if (Lorb == 1) {hcn = -0.239967;}
      else if (Lorb == 2) {hcn = -0.777542;}
      break;
    case 58:                               //Ce
      if (Lorb == 0) {hcn = 1.335148;}
      else if (Lorb == 1) {hcn = -0.198320;}
      else if (Lorb == 2) {hcn = -0.638958;}
      break;
    case 59:                               //Pr
      if (Lorb == 0) {hcn = 1.349441;}
      else if (Lorb == 1) {hcn = -0.198184;}
      else if (Lorb == 2) {hcn = -0.543909;}
      break;
    case 60:                               //Nd
      if (Lorb == 0) {hcn = 1.363734;}
      else if (Lorb == 1) {hcn = -0.198048;}
      else if (Lorb == 2) {hcn = -0.448861;}
      break;
    case 61:                               //Pm
      if (Lorb == 0) {hcn = 1.378028;}
      else if (Lorb == 1) {hcn = -0.197912;}
      else if (Lorb == 2) {hcn = -0.353812;}
      break;
    case 62:                               //Sm
      if (Lorb == 0) {hcn = 1.392321;}
      else if (Lorb == 1) {hcn = -0.197776;}
      else if (Lorb == 2) {hcn = -0.258764;}
      break;
    case 63:                               //Eu
      if (Lorb == 0) {hcn = 1.406614;}
      else if (Lorb == 1) {hcn = -0.197640;}
      else if (Lorb == 2) {hcn = -0.163715;}
      break;
    case 64:                               //Gd
      if (Lorb == 0) {hcn = 1.420907;}
      else if (Lorb == 1) {hcn = -0.197504;}
      else if (Lorb == 2) {hcn = -0.068667;}
      break;
    case 65:                               //Tb
      if (Lorb == 0) {hcn = 1.435200;}
      else if (Lorb == 1) {hcn = -0.197369;}
      else if (Lorb == 2) {hcn = 0.026382;}
      break;
    case 66:                               //Dy
      if (Lorb == 0) {hcn = 1.449494;}
      else if (Lorb == 1) {hcn = -0.197233;}
      else if (Lorb == 2) {hcn = 0.121430;}
      break;
    case 67:                               //Ho
      if (Lorb == 0) {hcn = 1.463787;}
      else if (Lorb == 1) {hcn = -0.197097;}
      else if (Lorb == 2) {hcn = 0.216479;}
      break;
    case 68:                               //Er
      if (Lorb == 0) {hcn = 1.478080;}
      else if (Lorb == 1) {hcn = -0.196961;}
      else if (Lorb == 2) {hcn = 0.311527;}
      break;
    case 69:                               //Tm
      if (Lorb == 0) {hcn = 1.492373;}
      else if (Lorb == 1) {hcn = -0.196825;}
      else if (Lorb == 2) {hcn = 0.406576;}
      break;
    case 70:                               //Yb
      if (Lorb == 0) {hcn = 1.506666;}
      else if (Lorb == 1) {hcn = -0.196689;}
      else if (Lorb == 2) {hcn = 0.501624;}
      break;
    case 71:                               //Lu
      if (Lorb == 0) {hcn = 1.520959;}
      else if (Lorb == 1) {hcn = -0.196553;}
      else if (Lorb == 2) {hcn = 0.596673;}
      break;
    case 72:                               //Hf
      if (Lorb == 0) {hcn = 2.271496;}
      else if (Lorb == 1) {hcn = -0.069771;}
      else if (Lorb == 2) {hcn = 0.176550;}
      break;
    case 73:                               //Ta
      if (Lorb == 0) {hcn = 0.988501;}
      else if (Lorb == 1) {hcn = -0.472540;}
      else if (Lorb == 2) {hcn = -0.620136;}
      break;
    case 74:                               //W
      if (Lorb == 0) {hcn = 2.543639;}
      else if (Lorb == 1) {hcn = 0.236479;}
      else if (Lorb == 2) {hcn = -0.192494;}
      break;
    case 75:                               //Re
      if (Lorb == 0) {hcn = 1.117569;}
      else if (Lorb == 1) {hcn = -1.335155;}
      else if (Lorb == 2) {hcn = -0.322139;}
      break;
    case 76:                               //Os
      if (Lorb == 0) {hcn = 0.346183;}
      else if (Lorb == 1) {hcn = -0.208758;}
      else if (Lorb == 2) {hcn = -0.095346;}
      break;
    case 77:                               //Ir
      if (Lorb == 0) {hcn = -0.123672;}
      else if (Lorb == 1) {hcn = -0.079864;}
      else if (Lorb == 2) {hcn = 0.051977;}
      break;
    case 78:                               //Pt
      if (Lorb == 0) {hcn = 1.139530;}
      else if (Lorb == 1) {hcn = 1.408029;}
      else if (Lorb == 2) {hcn = -0.204828;}
      break;
    case 79:                               //Au
      if (Lorb == 0) {hcn = 1.479337;}
      else if (Lorb == 1) {hcn = 1.048065;}
      else if (Lorb == 2) {hcn = -0.154462;}
      break;
    case 80:                               //Hg
      if (Lorb == 0) {hcn = -0.352252;}
      else if (Lorb == 1) {hcn = 0.205401;}
      break;
    case 81:                               //Tl
      if (Lorb == 0) {hcn = -0.255975;}
      else if (Lorb == 1) {hcn = 0.901364;}
      break;
    case 82:                               //Pb
      if (Lorb == 0) {hcn = -3.893459;}
      else if (Lorb == 1) {hcn = 3.437117;}
      break;
    case 83:                               //Bi
      if (Lorb == 0) {hcn = 0.160425;}
      else if (Lorb == 1) {hcn = 0.248659;}
      break;
    case 84:                               //Po
      if (Lorb == 0) {hcn = -0.046813;}
      else if (Lorb == 1) {hcn = -0.100437;}
      break;
    case 85:                               //At
      if (Lorb == 0) {hcn = -0.287369;}
      else if (Lorb == 1) {hcn = -0.007993;}
      else if (Lorb == 2) {hcn = 2.805808;}
      break;
    case 86:                               //Rn
      if (Lorb == 0) {hcn = -0.001712;}
      else if (Lorb == 1) {hcn = -0.005280;}
      else if (Lorb == 2) {hcn = -3.206020;}
  }
  return 0.1*hcn/au2eV;
}
double PaulingElectronegativity(int atm) {
  //atomic electronegativities according to Pauling
  double enpauling = 0.0;
  switch (atm) {
    case 1:
      enpauling = 2.20;
      break;
    case 2:
      enpauling = 3.00;
      break;
    case 3:
      enpauling = 0.98;
      break;
    case 4:
      enpauling = 1.57;
      break;
    case 5:
      enpauling = 2.04;
      break;
    case 6:
      enpauling = 2.55;
      break;
    case 7:
      enpauling = 3.04;
      break;
    case 8:
      enpauling = 3.44;
      break;
    case 9:
      enpauling = 3.98;
      break;
    case 10:
      enpauling = 4.50;
      break;
    case 11:
      enpauling = 0.93;
      break;
    case 12:
      enpauling = 1.31;
      break;
    case 13:
      enpauling = 1.61;
      break;
    case 14:
      enpauling = 1.90;
      break;
    case 15:
      enpauling = 2.19;
      break;
    case 16:
      enpauling = 2.58;
      break;
    case 17:
      enpauling = 3.16;
      break;
    case 18:
      enpauling = 3.50;
      break;
    case 19:
      enpauling = 0.82;
      break;
    case 20:
      enpauling = 1.00;
      break;
    case 21:
      enpauling = 1.36;
      break;
    case 22:
      enpauling = 1.54;
      break;
    case 23:
      enpauling = 1.63;
      break;
    case 24:
      enpauling = 1.66;
      break;
    case 25:
      enpauling = 1.55;
      break;
    case 26:
      enpauling = 1.83;
      break;
    case 27:
      enpauling = 1.88;
      break;
    case 28:
      enpauling = 1.91;
      break;
    case 29:
      enpauling = 1.90;
      break;
    case 30:
      enpauling = 1.65;
      break;
    case 31:
      enpauling = 1.81;
      break;
    case 32:
      enpauling = 2.01;
      break;
    case 33:
      enpauling = 2.18;
      break;
    case 34:
      enpauling = 2.55;
      break;
    case 35:
      enpauling = 2.96;
      break;
    case 36:
      enpauling = 3.00;
      break;
    case 37:
      enpauling = 0.82;
      break;
    case 38:
      enpauling = 0.95;
      break;
    case 39:
      enpauling = 1.22;
      break;
    case 40:
      enpauling = 1.33;
      break;
    case 41:
      enpauling = 1.60;
      break;
    case 42:
      enpauling = 2.16;
      break;
    case 43:
      enpauling = 1.90;
      break;
    case 44:
      enpauling = 2.20;
      break;
    case 45:
      enpauling = 2.28;
      break;
    case 46:
      enpauling = 2.20;
      break;
    case 47:
      enpauling = 1.93;
      break;
    case 48:
      enpauling = 1.69;
      break;
    case 49:
      enpauling = 1.78;
      break;
    case 50:
      enpauling = 1.96;
      break;
    case 51:
      enpauling = 2.05;
      break;
    case 52:
      enpauling = 2.10;
      break;
    case 53:
      enpauling = 2.66;
      break;
    case 54:
      enpauling = 2.60;
      break;
    case 55:
      enpauling = 0.79;
      break;
    case 56:
      enpauling = 0.89;
      break;
    case 57:
      enpauling = 1.10;
      break;
    case 58:
      enpauling = 1.12;
      break;
    case 59:
      enpauling = 1.13;
      break;
    case 60:
      enpauling = 1.14;
      break;
    case 61:
      enpauling = 1.15;
      break;
    case 62:
      enpauling = 1.17;
      break;
    case 63:
      enpauling = 1.18;
      break;
    case 64:
      enpauling = 1.20;
      break;
    case 65:
      enpauling = 1.21;
      break;
    case 66:
      enpauling = 1.22;
      break;
    case 67:
      enpauling = 1.23;
      break;
    case 68:
      enpauling = 1.24;
      break;
    case 69:
      enpauling = 1.25;
      break;
    case 70:
      enpauling = 1.26;
      break;
    case 71:
      enpauling = 1.27;
      break;
    case 72:
      enpauling = 1.30;
      break;
    case 73:
      enpauling = 1.50;
      break;
    case 74:
      enpauling = 2.36;
      break;
    case 75:
      enpauling = 1.90;
      break;
    case 76:
      enpauling = 2.20;
      break;
    case 77:
      enpauling = 2.20;
      break;
    case 78:
      enpauling = 2.28;
      break;
    case 79:
      enpauling = 2.54;
      break;
    case 80:
      enpauling = 2.00;
      break;
    case 81:
      enpauling = 1.62;
      break;
    case 82:
      enpauling = 2.33;
      break;
    case 83:
      enpauling = 2.02;
      break;
    case 84:
      enpauling = 2.00;
      break;
    case 85:
      enpauling = 2.20;
      break;
    case 86:
      enpauling = 2.20;
      break;
    case 87:
      enpauling = 1.50;
      break;
    case 88:
      enpauling = 1.50;
      break;
    case 89:
      enpauling = 1.50;
      break;
    case 90:
      enpauling = 1.50;
      break;
    case 91:
      enpauling = 1.50;
      break;
    case 92:
      enpauling = 1.50;
      break;
    case 93:
      enpauling = 1.50;
      break;
    case 94:
      enpauling = 1.50;
      break;
    case 95:
      enpauling = 1.50;
      break;
    case 96:
      enpauling = 1.50;
      break;
    case 97:
      enpauling = 1.50;
      break;
    case 98:
      enpauling = 1.50;
      break;
    case 99:
      enpauling = 1.50;
      break;
    case 100:
      enpauling = 1.50;
      break;
    case 101:
      enpauling = 1.50;
      break;
    case 102:
      enpauling = 1.50;
      break;
    case 103:
      enpauling = 1.50;
      break;
    case 104:
      enpauling = 1.50;
      break;
    case 105:
      enpauling = 1.50;
      break;
    case 106:
      enpauling = 1.50;
      break;
    case 107:
      enpauling = 1.50;
      break;
    case 108:
      enpauling = 1.50;
      break;
    case 109:
      enpauling = 1.50;
      break;
    case 110:
      enpauling = 1.50;
      break;
    case 111:
      enpauling = 1.50;
      break;
    case 112:
      enpauling = 1.50;
      break;
    case 113:
      enpauling = 1.50;
      break;
    case 114:
      enpauling = 1.50;
      break;
    case 115:
      enpauling = 1.50;
      break;
    case 116:
      enpauling = 1.50;
      break;
    case 117:
      enpauling = 1.50;
      break;
    case 118:
      enpauling = 1.50;
      break;
  }
  return enpauling;
}
double ShellPolynomialCoeff(int atm, int Lorb) {
  //{k}^{poly}_{A,l} parameters in equation 19 of paper
  double kpoly = 0.0; 
  switch (atm) {
    case 1:
      kpoly = -0.953618;
      break;
    case 2:
      if (Lorb == 0) {kpoly = -4.386816;}
      else if (Lorb == 1) {kpoly =  0.710647;}
      break;
    case 3:
      if (Lorb == 0) {kpoly = -4.750398;}
      else if (Lorb == 1) {kpoly = 20.424920;}
      break;
    case 4:
      if (Lorb == 0) {kpoly = -7.910394;}
      else if (Lorb == 1) {kpoly = -0.476438;}
      break;
    case 5:
      if (Lorb == 0) {kpoly = -5.183150;}
      else if (Lorb == 1) {kpoly = -2.453322;}
      break;
    case 6:
      if (Lorb == 0) {kpoly = -2.294321;}
      else if (Lorb == 1) {kpoly = -0.271102;}
      break;
    case 7:
      if (Lorb == 0) {kpoly = -8.506003;}
      else if (Lorb == 1) {kpoly = -2.504201;}
      break;
    case 8:
      if (Lorb == 0) {kpoly = -14.955291;}
      else if (Lorb == 1) {kpoly = -3.350819;}
      break;
    case 9:
      if (Lorb == 0) {kpoly = -13.011924;}
      else if (Lorb == 1) {kpoly = -12.300828;}
      break;
    case 10:
      if (Lorb == 0) {kpoly = -16.377759;}
      else if (Lorb == 1) {kpoly = -4.860545;}
      else if (Lorb == 2) {kpoly = -16.922333;}
      break;
    case 11:
      if (Lorb == 0) {kpoly = -4.033495;}
      else if (Lorb == 1) {kpoly = 20.873908;}
      break;
    case 12:
      if (Lorb == 0) {kpoly = -11.167374;}
      else if (Lorb == 1) {kpoly = 39.076962;}
      else if (Lorb == 2) {kpoly = 12.691061;}
      break;
    case 13:
      if (Lorb == 0) {kpoly = -10.678141;}
      else if (Lorb == 1) {kpoly = -12.442837;}
      else if (Lorb == 2) {kpoly = 16.311089;}
      break;
    case 14:
      if (Lorb == 0) {kpoly = 2.358522;}
      else if (Lorb == 1) {kpoly = -7.900406;}
      else if (Lorb == 2) {kpoly = 11.366185;}
      break;
    case 15:
      if (Lorb == 0) {kpoly = -19.831771;}
      else if (Lorb == 1) {kpoly = -5.515577;}
      else if (Lorb == 2) {kpoly = 26.397535;}
      break;
    case 16:
      if (Lorb == 0) {kpoly = -25.855520;}
      else if (Lorb == 1) {kpoly = -8.048064;}
      else if (Lorb == 2) {kpoly = 25.993857;}
      break;
    case 17:
      if (Lorb == 0) {kpoly = -16.562004;}
      else if (Lorb == 1) {kpoly = -6.986430;}
      else if (Lorb == 2) {kpoly = 38.045622;}
      break;
    case 18:
      if (Lorb == 0) {kpoly = -23.893924;}
      else if (Lorb == 1) {kpoly = -3.727322;}
      else if (Lorb == 2) {kpoly = 26.812872;}
      break;
    case 19:
      if (Lorb == 0) {kpoly = -6.076060;}
      else if (Lorb == 1) {kpoly = 21.187329;}
      break;
    case 20:
      if (Lorb == 0) {kpoly = -9.718719;}
      else if (Lorb == 1) {kpoly = 31.973372;}
      else if (Lorb == 2) {kpoly = 9.528649;}
      break;
    case 21:
      if (Lorb == 0) {kpoly =  0.686569;}
      else if (Lorb == 1) {kpoly = 38.044933;}
      else if (Lorb == 2) {kpoly = -34.502301;}
      break;
    case 22:
      if (Lorb == 0) {kpoly = 4.561234;}
      else if (Lorb == 1) {kpoly = 51.801626;}
      else if (Lorb == 2) {kpoly = -27.724389;}
      break;
    case 23:
      if (Lorb == 0) {kpoly = 9.702482;}
      else if (Lorb == 1) {kpoly = 51.178301;}
      else if (Lorb == 2) {kpoly = -29.827563;}
      break;
    case 24:
      if (Lorb == 0) {kpoly = 13.376234;}
      else if (Lorb == 1) {kpoly = 48.092152;}
      else if (Lorb == 2) {kpoly = -27.971622;}
      break;
    case 25:
      if (Lorb == 0) {kpoly = 28.519691;}
      else if (Lorb == 1) {kpoly = 26.346555;}
      else if (Lorb == 2) {kpoly = -31.255885;}
      break;
    case 26:
      if (Lorb == 0) {kpoly = 11.527794;}
      else if (Lorb == 1) {kpoly = 39.459890;}
      else if (Lorb == 2) {kpoly = -28.614961;}
      break;
    case 27:
      if (Lorb == 0) {kpoly = 9.168460;}
      else if (Lorb == 1) {kpoly = 25.424719;}
      else if (Lorb == 2) {kpoly = -22.355636;}
      break;
    case 28:
      if (Lorb == 0) {kpoly = 20.839550;}
      else if (Lorb == 1) {kpoly = 30.886445;}
      else if (Lorb == 2) {kpoly = -25.385640;}
      break;
    case 29:
      if (Lorb == 0) {kpoly = 17.798264;}
      else if (Lorb == 1) {kpoly = 14.977818;}
      else if (Lorb == 2) {kpoly = -26.508943;}
      break;
    case 30:
      if (Lorb == 0) {kpoly = -9.240315;}
      else if (Lorb == 1) {kpoly = 22.271839;}
      break;
    case 31:
      if (Lorb == 0) {kpoly = -19.018211;}
      else if (Lorb == 1) {kpoly = -1.137790;}
      else if (Lorb == 2) {kpoly = 35.401863;}
      break;
    case 32:
      if (Lorb == 0) {kpoly = -21.333673;}
      else if (Lorb == 1) {kpoly = -9.749043;}
      else if (Lorb == 2) {kpoly = 28.634699;}
      break;
    case 33:
      if (Lorb == 0) {kpoly = -23.820664;}
      else if (Lorb == 1) {kpoly = -10.644177;}
      else if (Lorb == 2) {kpoly = 30.711130;}
      break;
    case 34:
      if (Lorb == 0) {kpoly = -24.506414;}
      else if (Lorb == 1) {kpoly = -13.765750;}
      else if (Lorb == 2) {kpoly = 29.611132;}
      break;
    case 35:
      if (Lorb == 0) {kpoly = -25.005079;}
      else if (Lorb == 1) {kpoly = -14.520078;}
      else if (Lorb == 2) {kpoly = 36.614038;}
      break;
    case 36:
      if (Lorb == 0) {kpoly = -32.658744;}
      else if (Lorb == 1) {kpoly = -13.600093;}
      else if (Lorb == 2) {kpoly = 23.204749;}
      break;
    case 37:
      if (Lorb == 0) {kpoly = 4.325404;}
      else if (Lorb == 1) {kpoly = 23.255114;}
      break;
    case 38:
      if (Lorb == 0) {kpoly = -14.506782;}
      else if (Lorb == 1) {kpoly = 20.213982;}
      else if (Lorb == 2) {kpoly = 10.816244;}
      break;
    case 39:
      if (Lorb == 0) {kpoly = -2.125866;}
      else if (Lorb == 1) {kpoly = 52.161869;}
      else if (Lorb == 2) {kpoly = -39.529477;}
      break;
    case 40:
      if (Lorb == 0) {kpoly = 7.538902;}
      else if (Lorb == 1) {kpoly = 58.914069;}
      else if (Lorb == 2) {kpoly = -28.358876;}
      break;
    case 41:
      if (Lorb == 0) {kpoly = -5.141080;}
      else if (Lorb == 1) {kpoly = 55.654247;}
      else if (Lorb == 2) {kpoly = -27.963746;}
      break;
    case 42:
      if (Lorb == 0) {kpoly = -0.583137;}
      else if (Lorb == 1) {kpoly = 29.199580;}
      else if (Lorb == 2) {kpoly = -22.573688;}
      break;
    case 43:
      if (Lorb == 0) {kpoly = 36.096045;}
      else if (Lorb == 1) {kpoly = 25.095742;}
      else if (Lorb == 2) {kpoly = -27.342552;}
      break;
    case 44:
      if (Lorb == 0) {kpoly = 10.106270;}
      else if (Lorb == 1) {kpoly = 34.028722;}
      else if (Lorb == 2) {kpoly = -27.583213;}
      break;
    case 45:
      if (Lorb == 0) {kpoly = 15.413317;}
      else if (Lorb == 1) {kpoly = 31.070684;}
      else if (Lorb == 2) {kpoly = -19.656064;}
      break;
    case 46:
      if (Lorb == 0) {kpoly = 6.200145;}
      else if (Lorb == 1) {kpoly = 45.341322;}
      else if (Lorb == 2) {kpoly = -27.173113;}
      break;
    case 47:
      if (Lorb == 0) {kpoly = 1.091490;}
      else if (Lorb == 1) {kpoly = 11.561444;}
      else if (Lorb == 2) {kpoly = -16.490742;}
      break;
    case 48:
      if (Lorb == 0) {kpoly = -6.076865;}
      else if (Lorb == 1) {kpoly = 37.671924;}
      break;
    case 49:
      if (Lorb == 0) {kpoly = -21.938535;}
      else if (Lorb == 1) {kpoly = -1.949646;}
      else if (Lorb == 2) {kpoly = 31.354500;}
      break;
    case 50:
      if (Lorb == 0) {kpoly = -17.518223;}
      else if (Lorb == 1) {kpoly = -7.802870;}
      else if (Lorb == 2) {kpoly = 12.611148;}
      break;
    case 51:
      if (Lorb == 0) {kpoly = -17.543533;}
      else if (Lorb == 1) {kpoly = -12.494609;}
      else if (Lorb == 2) {kpoly = 30.872654;}
      break;
    case 52:
      if (Lorb == 0) {kpoly = -24.893875;}
      else if (Lorb == 1) {kpoly = -11.231996;}
      else if (Lorb == 2) {kpoly = 31.843215;}
      break;
    case 53:
      if (Lorb == 0) {kpoly = -26.957547;}
      else if (Lorb == 1) {kpoly = -14.183312;}
      else if (Lorb == 2) {kpoly = 28.211905;}
      break;
    case 54:
      if (Lorb == 0) {kpoly = -31.096492;}
      else if (Lorb == 1) {kpoly = -16.197851;}
      else if (Lorb == 2) {kpoly = 19.048983;}
      break;
    case 55:
      if (Lorb == 0) {kpoly = -0.713637;}
      else if (Lorb == 1) {kpoly = 20.636967;}
      break;
    case 56:
      if (Lorb == 0) {kpoly = -14.036589;}
      else if (Lorb == 1) {kpoly = 18.774072;}
      else if (Lorb == 2) {kpoly = 11.389672;}
      break;
    case 57:
      if (Lorb == 0) {kpoly = -6.732008;}
      else if (Lorb == 1) {kpoly = 54.136379;}
      else if (Lorb == 2) {kpoly = -37.820063;}
      break;
    case 58:
      if (Lorb == 0) {kpoly = -6.107742;}
      else if (Lorb == 1) {kpoly = 37.663425;}
      else if (Lorb == 2) {kpoly = -41.989203;}
      break;
    case 59:
      if (Lorb == 0) {kpoly = -6.040165;}
      else if (Lorb == 1) {kpoly = 38.194756;}
      else if (Lorb == 2) {kpoly = -41.286513;}
      break;
    case 60:
      if (Lorb == 0) {kpoly = -5.972588;}
      else if (Lorb == 1) {kpoly = 38.726087;}
      else if (Lorb == 2) {kpoly = -40.583824;}
      break;
    case 61:
      if (Lorb == 0) {kpoly = -5.905010;}
      else if (Lorb == 1) {kpoly = 39.257419;}
      else if (Lorb == 2) {kpoly = -39.881134;}
      break;
    case 62:
      if (Lorb == 0) {kpoly = -5.837433;}
      else if (Lorb == 1) {kpoly = 39.788750;}
      else if (Lorb == 2) {kpoly = -39.178445;}
      break;
    case 63:
      if (Lorb == 0) {kpoly = -5.769856;}
      else if (Lorb == 1) {kpoly = 40.320081;}
      else if (Lorb == 2) {kpoly = -38.475755;}
      break;
    case 64:
      if (Lorb == 0) {kpoly = -5.702279;}
      else if (Lorb == 1) {kpoly = 40.851412;}
      else if (Lorb == 2) {kpoly = -37.773066;}
      break;
    case 65:
      if (Lorb == 0) {kpoly = -5.634701;}
      else if (Lorb == 1) {kpoly = 41.382743;}
      else if (Lorb == 2) {kpoly = -37.070376;}
      break;
    case 66:
      if (Lorb == 0) {kpoly = -5.567124;}
      else if (Lorb == 1) {kpoly = 41.914074;}
      else if (Lorb == 2) {kpoly = -36.367686;}
      break;
    case 67:
      if (Lorb == 0) {kpoly = -5.499547;}
      else if (Lorb == 1) {kpoly = 42.445406;}
      else if (Lorb == 2) {kpoly = -35.664997;}
      break;
    case 68:
      if (Lorb == 0) {kpoly = -5.431970;}
      else if (Lorb == 1) {kpoly = 42.976737;}
      else if (Lorb == 2) {kpoly = -34.962307;}
      break;
    case 69:
      if (Lorb == 0) {kpoly = -5.364392;}
      else if (Lorb == 1) {kpoly = 43.508068;}
      else if (Lorb == 2) {kpoly = -34.259618;}
      break;
    case 70:
      if (Lorb == 0) {kpoly = -5.296815;}
      else if (Lorb == 1) {kpoly = 44.039399;}
      else if (Lorb == 2) {kpoly = -33.556928;}
      break;
    case 71:
      if (Lorb == 0) {kpoly = -5.229238;}
      else if (Lorb == 1) {kpoly = 44.570730;}
      else if (Lorb == 2) {kpoly = -32.854239;}
      break;
    case 72:
      if (Lorb == 0) {kpoly = -2.731930;}
      else if (Lorb == 1) {kpoly = 33.515028;}
      else if (Lorb == 2) {kpoly = -34.095699;}
      break;
    case 73:
      if (Lorb == 0) {kpoly = -15.707734;}
      else if (Lorb == 1) {kpoly = 60.186026;}
      else if (Lorb == 2) {kpoly = -30.396333;}
      break;
    case 74:
      if (Lorb == 0) {kpoly = 6.208982;}
      else if (Lorb == 1) {kpoly = 49.273770;}
      else if (Lorb == 2) {kpoly = -25.677080;}
      break;
    case 75:
      if (Lorb == 0) {kpoly = 13.890139;}
      else if (Lorb == 1) {kpoly = 33.973283;}
      else if (Lorb == 2) {kpoly = -31.723059;}
      break;
    case 76:
      if (Lorb == 0) {kpoly = 21.316847;}
      else if (Lorb == 1) {kpoly = 28.097231;}
      else if (Lorb == 2) {kpoly = -28.461111;}
      break;
    case 77:
      if (Lorb == 0) {kpoly = 20.733836;}
      else if (Lorb == 1) {kpoly = 18.303210;}
      else if (Lorb == 2) {kpoly = -24.693443;}
      break;
    case 78:
      if (Lorb == 0) {kpoly =  6.737556;}
      else if (Lorb == 1) {kpoly = 19.259455;}
      else if (Lorb == 2) {kpoly = -27.243898;}
      break;
    case 79:
      if (Lorb == 0) {kpoly = 4.691539;}
      else if (Lorb == 1) {kpoly = 25.250274;}
      else if (Lorb == 2) {kpoly = -6.410815;}
      break;
    case 80:
      if (Lorb == 0) {kpoly = -9.833449;}
      else if (Lorb == 1) {kpoly = 15.628858;}
      break;
    case 81:
      if (Lorb == 0) {kpoly = -22.942207;}
      else if (Lorb == 1) {kpoly = 13.109821;}
      break;
    case 82:
      if (Lorb == 0) {kpoly = -22.955084;}
      else if (Lorb == 1) {kpoly = -8.805272;}
      break;
    case 83:
      if (Lorb == 0) {kpoly = -21.750059;}
      else if (Lorb == 1) {kpoly = -10.773930;}
      break;
    case 84:
      if (Lorb == 0) {kpoly = -20.923290;}
      else if (Lorb == 1) {kpoly = -18.426448;}
      break;
    case 85:
      if (Lorb == 0) {kpoly = -30.549981;}
      else if (Lorb == 1) {kpoly = -17.108490;}
      else if (Lorb == 2) {kpoly = 23.824993;}
      break;
    case 86:
      if (Lorb == 0) {kpoly = -35.245372;}
      else if (Lorb == 1) {kpoly = -11.989735;}
      else if (Lorb == 2) {kpoly = 21.167024;}
      break;
  }
  return kpoly;
}
double Yeff(int atmnr) {
  //effective nuclear charges for atoms
  double yeff = 0.0;
  switch (atmnr) {
    case 1:
      yeff = 1.105388;
      break;
    case 2:
      yeff = 1.094283;
      break;
    case 3:
      yeff = 1.289367;
      break;
    case 4:
      yeff = 4.221216;
      break;
    case 5:
      yeff = 7.192431;
      break;
    case 6:
      yeff = 4.231078;
      break;
    case 7:
      yeff = 5.242592;
      break;
    case 8:
      yeff = 5.784415;
      break;
    case 9:
      yeff = 7.021486;
      break;
    case 10:
      yeff = 11.041068;
      break;
    case 11:
      yeff = 5.244917;
      break;
    case 12:
      yeff = 18.083164;
      break;
    case 13:
      yeff = 17.867328;
      break;
    case 14:
      yeff = 40.001111;
      break;
    case 15:
      yeff = 19.683502;
      break;
    case 16:
      yeff = 14.995090;
      break;
    case 17:
      yeff = 17.353134;
      break;
    case 18:
      yeff = 7.266606;
      break;
    case 19:
      yeff = 10.439482;
      break;
    case 20:
      yeff = 14.786701;
      break;
    case 21:
      yeff = 8.004267;
      break;
    case 22:
      yeff = 12.036336;
      break;
    case 23:
      yeff = 15.677873;
      break;
    case 24:
      yeff = 19.517914;
      break;
    case 25:
      yeff = 18.760605;
      break;
    case 26:
      yeff = 20.360089;
      break;
    case 27:
      yeff = 27.127744;
      break;
    case 28:
      yeff = 10.533269;
      break;
    case 29:
      yeff = 9.913846;
      break;
    case 30:
      yeff = 22.099503;
      break;
    case 31:
      yeff = 31.146750;
      break;
    case 32:
      yeff = 42.100144;
      break;
    case 33:
      yeff = 39.147587;
      break;
    case 34:
      yeff = 27.426779;
      break;
    case 35:
      yeff = 32.845361;
      break;
    case 36:
      yeff = 17.363803;
      break;
    case 37:
      yeff = 44.338211;
      break;
    case 38:
      yeff = 34.365525;
      break;
    case 39:
      yeff = 17.326237;
      break;
    case 40:
      yeff = 24.263093;
      break;
    case 41:
      yeff = 30.562732;
      break;
    case 42:
      yeff = 48.312796;
      break;
    case 43:
      yeff = 44.779882;
      break;
    case 44:
      yeff = 28.070247;
      break;
    case 45:
      yeff = 38.035941;
      break;
    case 46:
      yeff = 28.674700;
      break;
    case 47:
      yeff = 6.493286;
      break;
    case 48:
      yeff = 26.226628;
      break;
    case 49:
      yeff = 63.854240;
      break;
    case 50:
      yeff = 80.053438;
      break;
    case 51:
      yeff = 77.057560;
      break;
    case 52:
      yeff = 48.614745;
      break;
    case 53:
      yeff = 63.319176;
      break;
    case 54:
      yeff = 51.188398;
      break;
    case 55:
      yeff = 67.249039;
      break;
    case 56:
      yeff = 46.984607;
      break;
    case 57:
      yeff = 50.927529;
      break;
    case 58:
      yeff = 48.676714;
      break;
    case 59:
      yeff = 47.669448;
      break;
    case 60:
      yeff = 46.662183;
      break;
    case 61:
      yeff = 45.654917;
      break;
    case 62:
      yeff = 44.647651;
      break;
    case 63:
      yeff = 43.640385;
      break;
    case 64:
      yeff = 42.633120;
      break;
    case 65:
      yeff = 41.625854;
      break;
    case 66:
      yeff = 40.618588;
      break;
    case 67:
      yeff = 39.611322;
      break;
    case 68:
      yeff = 38.604057;
      break;
    case 69:
      yeff = 37.596791;
      break;
    case 70:
      yeff = 36.589525;
      break;
    case 71:
      yeff = 35.582259;
      break;
    case 72:
      yeff = 40.186772;
      break;
    case 73:
      yeff = 54.666156;
      break;
    case 74:
      yeff = 55.899801;
      break;
    case 75:
      yeff = 80.410086;
      break;
    case 76:
      yeff = 62.809871;
      break;
    case 77:
      yeff = 56.045639;
      break;
    case 78:
      yeff = 53.881425;
      break;
    case 79:
      yeff = 14.711475;
      break;
    case 80:
      yeff = 51.577544;
      break;
    case 81:
      yeff = 58.801614;
      break;
    case 82:
      yeff = 102.368258;
      break;
    case 83:
      yeff = 132.896832;
      break;
    case 84:
      yeff = 52.301232;
      break;
    case 85:
      yeff = 81.771063;
      break;
    case 86:
      yeff = 128.133580;
      break;
  }
  return yeff;
}
double AlphaGFN2(int atmnmr) {
  double alpha = 0.0;
  switch (atmnmr) {
    case 1:
      alpha = 2.213717;
      break;
    case 2:
      alpha = 3.604670;
      break;
    case 3:
      alpha = 0.475307;
      break;
    case 4:
      alpha = 0.939696;
      break;
    case 5:
      alpha = 1.373856;
      break;
    case 6:
      alpha = 1.247655;
      break;
    case 7:
      alpha = 1.682689;
      break;
    case 8:
      alpha = 2.165712;
      break;
    case 9:
      alpha = 2.421394;
      break;
    case 10:
      alpha = 3.318479;
      break;
    case 11:
      alpha = 0.572728;
      break;
    case 12:
      alpha = 0.917975;
      break;
    case 13:
      alpha = 0.876623;
      break;
    case 14:
      alpha = 1.187323;
      break;
    case 15:
      alpha = 1.143343;
      break;
    case 16:
      alpha = 1.214553;
      break;
    case 17:
      alpha = 1.577144;
      break;
    case 18:
      alpha = 0.896198;
      break;
    case 19:
      alpha = 0.482206;
      break;
    case 20:
      alpha = 0.683051;
      break;
    case 21:
      alpha = 0.574299;
      break;
    case 22:
      alpha = 0.723104;
      break;
    case 23:
      alpha = 0.928532;
      break;
    case 24:
      alpha = 0.966993;
      break;
    case 25:
      alpha = 1.071100;
      break;
    case 26:
      alpha = 1.113422;
      break;
    case 27:
      alpha = 1.241717;
      break;
    case 28:
      alpha = 1.077516;
      break;
    case 29:
      alpha = 0.998768;
      break;
    case 30:
      alpha = 1.160262;
      break;
    case 31:
      alpha = 1.122923;
      break;
    case 32:
      alpha = 1.222349;
      break;
    case 33:
      alpha = 1.249372;
      break;
    case 34:
      alpha = 1.230284;
      break;
    case 35:
      alpha = 1.296174;
      break;
    case 36:
      alpha = 0.908074;
      break;
    case 37:
      alpha = 0.574054;
      break;
    case 38:
      alpha = 0.697345;
      break;
    case 39:
      alpha = 0.706172;
      break;
    case 40:
      alpha = 0.681106;
      break;
    case 41:
      alpha = 0.865552;
      break;
    case 42:
      alpha = 1.034519;
      break;
    case 43:
      alpha = 1.019565;
      break;
    case 44:
      alpha = 1.031669;
      break;
    case 45:
      alpha = 1.094599;
      break;
    case 46:
      alpha = 1.092745;
      break;
    case 47:
      alpha = 0.678344;
      break;
    case 48:
      alpha = 0.936236;
      break;
    case 49:
      alpha = 1.024007;
      break;
    case 50:
      alpha = 1.139959;
      break;
    case 51:
      alpha = 1.122937;
      break;
    case 52:
      alpha = 1.000712;
      break;
    case 53:
      alpha = 1.017946;
      break;
    case 54:
      alpha = 1.012036;
      break;
    case 55:
      alpha = 0.585257;
      break;
    case 56:
      alpha = 0.716259;
      break;
    case 57:
      alpha = 0.737643;
      break;
    case 58:
      alpha = 0.729950;
      break;
    case 59:
      alpha = 0.734624;
      break;
    case 60:
      alpha = 0.739299;
      break;
    case 61:
      alpha = 0.743973;
      break;
    case 62:
      alpha = 0.748648;
      break;
    case 63:
      alpha = 0.753322;
      break;
    case 64:
      alpha = 0.757996;
      break;
    case 65:
      alpha = 0.762671;
      break;
    case 66:
      alpha = 0.767345;
      break;
    case 67:
      alpha = 0.772020;
      break;
    case 68:
      alpha = 0.776694;
      break;
    case 69:
      alpha = 0.781368;
      break;
    case 70:
      alpha = 0.786043;
      break;
    case 71:
      alpha = 0.790717;
      break;
    case 72:
      alpha = 0.852852;
      break;
    case 73:
      alpha = 0.990234;
      break;
    case 74:
      alpha = 1.018805;
      break;
    case 75:
      alpha = 1.170412;
      break;
    case 76:
      alpha = 1.221937;
      break;
    case 77:
      alpha = 1.197148;
      break;
    case 78:
      alpha = 1.204081;
      break;
    case 79:
      alpha = 0.919210;
      break;
    case 80:
      alpha = 1.137360;
      break;
    case 81:
      alpha = 1.399312;
      break;
    case 82:
      alpha = 1.179922;
      break;
    case 83:
      alpha = 1.130860;
      break;
    case 84:
      alpha = 0.957939;
      break;
    case 85:
      alpha = 0.963878;
      break;
    case 86:
      alpha = 0.965577;
      break;
  }
  return alpha;
}
double Hardness(int atmnmr) {
  //chemical hardness of atoms
  double eta = 0.0;
  switch (atmnmr) {
    case 1:
      eta = 0.405771;
      break;
    case 2:
      eta = 0.642029;
      break;
    case 3:
      eta = 0.245006;
      break;
    case 4:
      eta = 0.684789;
      break;
    case 5:
      eta = 0.513556;
      break;
    case 6:
      eta = 0.538015;
      break;
    case 7:
      eta = 0.461493;
      break;
    case 8:
      eta = 0.451896;
      break;
    case 9:
      eta = 0.531518;
      break;
    case 10:
      eta = 0.850000;
      break;
    case 11:
      eta = 0.271056;
      break;
    case 12:
      eta = 0.344822;
      break;
    case 13:
      eta = 0.364801;
      break;
    case 14:
      eta = 0.720000;
      break;
    case 15:
      eta = 0.297739;
      break;
    case 16:
      eta = 0.339971;
      break;
    case 17:
      eta = 0.248514;
      break;
    case 18:
      eta = 0.502376;
      break;
    case 19:
      eta = 0.247602;
      break;
    case 20:
      eta = 0.320378;
      break;
    case 21:
      eta = 0.472633;
      break;
    case 22:
      eta = 0.513586;
      break;
    case 23:
      eta = 0.589187;
      break;
    case 24:
      eta = 0.396299;
      break;
    case 25:
      eta = 0.346651;
      break;
    case 26:
      eta = 0.271594;
      break;
    case 27:
      eta = 0.477760;
      break;
    case 28:
      eta = 0.344970;
      break;
    case 29:
      eta = 0.202969;
      break;
    case 30:
      eta = 0.564152;
      break;
    case 31:
      eta = 0.432236;
      break;
    case 32:
      eta = 0.802051;
      break;
    case 33:
      eta = 0.571748;
      break;
    case 34:
      eta = 0.235052;
      break;
    case 35:
      eta = 0.261253;
      break;
    case 36:
      eta = 0.424373;
      break;
    case 37:
      eta = 0.210481;
      break;
    case 38:
      eta = 0.340000;
      break;
    case 39:
      eta = 0.711958;
      break;
    case 40:
      eta = 0.461440;
      break;
    case 41:
      eta = 0.952957;
      break;
    case 42:
      eta = 0.586134;
      break;
    case 43:
      eta = 0.368054;
      break;
    case 44:
      eta = 0.711205;
      break;
    case 45:
      eta = 0.509183;
      break;
    case 46:
      eta = 0.273310;
      break;
    case 47:
      eta = 0.263740;
      break;
    case 48:
      eta = 0.392012;
      break;
    case 49:
      eta = 0.461812;
      break;
    case 50:
      eta = 0.900000;
      break;
    case 51:
      eta = 0.942294;
      break;
    case 52:
      eta = 0.750000;
      break;
    case 53:
      eta = 0.383124;
      break;
    case 54:
      eta = 0.424164;
      break;
    case 55:
      eta = 0.236569;
      break;
    case 56:
      eta = 0.245937;
      break;
    case 57:
      eta = 0.597716;
      break;
    case 58:
      eta = 0.662889;
      break;
    case 59:
      eta = 0.660710;
      break;
    case 60:
      eta = 0.658531;
      break;
    case 61:
      eta = 0.656352;
      break;
    case 62:
      eta = 0.654173;
      break;
    case 63:
      eta = 0.651994;
      break;
    case 64:
      eta = 0.649815;
      break;
    case 65:
      eta = 0.647635;
      break;
    case 66:
      eta = 0.645456;
      break;
    case 67:
      eta = 0.643277;
      break;
    case 68:
      eta = 0.641098;
      break;
    case 69:
      eta = 0.638919;
      break;
    case 70:
      eta = 0.636740;
      break;
    case 71:
      eta = 0.634561;
      break;
    case 72:
      eta = 0.662597;
      break;
    case 73:
      eta = 0.449812;
      break;
    case 74:
      eta = 0.685426;
      break;
    case 75:
      eta = 0.224623;
      break;
    case 76:
      eta = 0.364388;
      break;
    case 77:
      eta = 0.548507;
      break;
    case 78:
      eta = 0.353574;
      break;
    case 79:
      eta = 0.438997;
      break;
    case 80:
      eta = 0.457611;
      break;
    case 81:
      eta = 0.418841;
      break;
    case 82:
      eta = 0.168152;
      break;
    case 83:
      eta = 0.900000;
      break;
    case 84:
      eta = 1.023267;
      break;
    case 85:
      eta = 0.288848;
      break;
    case 86:
      eta = 0.303400;
      break;
  }
  return eta;
}
double HardnessLScal(int atmnr, int Lorb) {
  double kscal = 0.0;
  switch (atmnr) {
    case 1:
      break;
    case 2:
      break;
    case 3:
      if (Lorb == 1) {kscal = 0.1972612;}
      break;
    case 4:
      if (Lorb == 1) {kscal = 0.9658467;}
      break;
    case 5:
      if (Lorb == 1) {kscal = 0.3994080;}
      break;
    case 6:
      if (Lorb == 1) {kscal = 0.1056358;}
      break;
    case 7:
      if (Lorb == 1) {kscal = 0.1164892;}
      break;
    case 8:
      if (Lorb == 1) {kscal = 0.1497020;}
      break;
    case 9:
      if (Lorb == 1) {kscal = 0.1677376;}
      break;
    case 10:
      if (Lorb == 1) {kscal = 0.1190576;}
      else if (Lorb == 2) {kscal = -0.3200000;}
      break;
    case 11:
      if (Lorb == 1) {kscal = 0.1018894;}
      break;
    case 12:
      if (Lorb == 1) {kscal = 1.4000000;}
      else if (Lorb == 2) {kscal = -0.0500000;}
      break;
    case 13:
      if (Lorb == 1) {kscal = -0.0603699;}
      else if (Lorb == 2) {kscal = 0.2000000;}
      break;
    case 14:
      if (Lorb == 1) {kscal = -0.5580042;}
      else if (Lorb == 2) {kscal = -0.2300000;}
      break;
    case 15:
      if (Lorb == 1) {kscal = -0.1558060;}
      else if (Lorb == 2) {kscal = -0.3500000;}
      break;
    case 16:
      if (Lorb == 1) {kscal = -0.1085866;}
      else if (Lorb == 2) {kscal = -0.2500000;}
      break;
    case 17:
      if (Lorb == 1) {kscal = 0.4989400;}
      else if (Lorb == 2) {kscal = 0.5000000;}
      break;
    case 18:
      if (Lorb == 1) {kscal = -0.0461133;}
      else if (Lorb == 2) {kscal = -0.0100000;}
      break;
    case 19:
      if (Lorb == 1) {kscal =  0.3483655;}
      break;
    case 20:
      if (Lorb == 1) {kscal = 1.5000000;}
      else if (Lorb == 2) {kscal = -0.2500000;}
      break;
    case 21:
      if (Lorb == 1) {kscal = -0.0800000;}
      else if (Lorb == 2) {kscal = -0.2046716;}
      break;
    case 22:
      if (Lorb == 1) {kscal = -0.3800000;}
      else if (Lorb == 2) {kscal = -0.4921114;}
      break;
    case 23:
      if (Lorb == 1) {kscal = -0.4500000;}
      else if (Lorb == 2) {kscal = -0.0379088;}
      break;
    case 24:
      if (Lorb == 1) {kscal = -0.4700000;}
      else if (Lorb == 2) {kscal = 0.7405872;}
      break;
    case 25:
      if (Lorb == 1) {kscal = -0.6000000;}
      else if (Lorb == 2) {kscal = 0.0545811;}
      break;
    case 26:
      if (Lorb == 1) {kscal = -0.6500000;}
      else if (Lorb == 2) {kscal = 0.4046615;}
      break;
    case 27:
      if (Lorb == 1) {kscal = -0.6500000;}
      else if (Lorb == 2) {kscal = -0.2418493;}
      break;
    case 28:
      if (Lorb == 1) {kscal = -0.6000000;}
      else if (Lorb == 2) {kscal = -0.0611188;}
      break;
    case 29:
      if (Lorb == 1) {kscal = 0.0700000;}
      else if (Lorb == 2) {kscal = 1.3333066;}
      break;
    case 30:
      if (Lorb == 1) {kscal = 0.0684343;}
      break;
    case 31:
      if (Lorb == 1) {kscal = -0.5416555;}
      else if (Lorb == 2) {kscal = -0.3000000;}
      break;
    case 32:
      if (Lorb == 1) {kscal = -0.3809089;}
      else if (Lorb == 2) {kscal = -0.1500000;}
      break;
    case 33:
      if (Lorb == 1) {kscal = -0.4104743;}
      else if (Lorb == 2) {kscal = -0.5000000;}
      break;
    case 34:
      if (Lorb == 1) {kscal = 0.1192113;}
      else if (Lorb == 2) {kscal = -0.2500000;}
      break;
    case 35:
      if (Lorb == 1) {kscal = 0.5203002;}
      else if (Lorb == 2) {kscal = 0.4000000;}
      break;
    case 36:
      if (Lorb == 1) {kscal = -0.2503223;}
      else if (Lorb == 2) {kscal = -0.0700000;}
      break;
    case 37:
      if (Lorb == 1) {kscal = 0.9386493;}
      break;
    case 38:
      if (Lorb == 1) {kscal = 1.5000000;}
      else if (Lorb == 2) {kscal = -0.2500000;}
      break;
    case 39:
      if (Lorb == 1) {kscal = -0.4500000;}
      else if (Lorb == 2) {kscal = -0.3349288;}
      break;
    case 40:
      if (Lorb == 1) {kscal = -0.1100000;}
      else if (Lorb == 2) {kscal = -0.4422630;}
      break;
    case 41:
      if (Lorb == 1) {kscal = -0.0500000;}
      else if (Lorb == 2) {kscal = -0.3562950;}
      break;
    case 42:
      if (Lorb == 1) {kscal = -0.3000000;}
      else if (Lorb == 2) {kscal = -0.4301371;}
      break;
    case 43:
      if (Lorb == 1) {kscal = -0.6000000;}
      else if (Lorb == 2) {kscal = 0.3956819;}
      break;
    case 44:
      if (Lorb == 1) {kscal = -0.6500000;}
      else if (Lorb == 2) {kscal = -0.3052305;}
      break;
    case 45:
      if (Lorb == 1) {kscal = -0.6500000;}
      else if (Lorb == 2) {kscal = -0.1881774;}
      break;
    case 46:
      if (Lorb == 1) {kscal = -0.6000000;}
      else if (Lorb == 2) {kscal = 0.0931707;}
      break;
    case 47:
      if (Lorb == 1) {kscal = -0.0300000;}
      else if (Lorb == 2) {kscal = 0.8024848;}
      break;
    case 48:
      if (Lorb == 1) {kscal = 0.2388669;}
      break;
    case 49:
      if (Lorb == 1) {kscal = -0.5867460;}
      else if (Lorb == 2) {kscal = -0.2800000;}
      break;
    case 50:
      if (Lorb == 1) {kscal = -0.5090746;}
      else if (Lorb == 2) {kscal = -0.0600000;}
      break;
    case 51:
      if (Lorb == 1) {kscal = -0.6278501;}
      else if (Lorb == 2) {kscal = -0.5500000;}
      break;
    case 52:
      if (Lorb == 1) {kscal = -0.1555334;}
      else if (Lorb == 2) {kscal = 0.0600000;}
      break;
    case 53:
      if (Lorb == 1) {kscal = -0.0338735;}
      else if (Lorb == 2) {kscal = 0.3000000;}
      break;
    case 54:
      if (Lorb == 1) {kscal = -0.2302667;}
      else if (Lorb == 2) {kscal = -0.2300000;}
      break;
    case 55:
      if (Lorb == 1) {kscal = 0.2494305;}
      break;
    case 56:
      if (Lorb == 1) {kscal =  2.2247532;}
      else if (Lorb == 2) {kscal = -0.2300000;}
      break;
    case 57:
      if (Lorb == 1) {kscal = -0.3000000;}
      else if (Lorb == 2) {kscal = -0.4699666;}
      break;
    case 58:
      if (Lorb == 1) {kscal = -0.3000000;}
      else if (Lorb == 2) {kscal = -0.5539659;}
      break;
    case 59:
      if (Lorb == 1) {kscal = -0.2769230;}
      else if (Lorb == 2) {kscal = -0.5462784;}
      break;
    case 60:
      if (Lorb == 1) {kscal = -0.2538460;}
      else if (Lorb == 2) {kscal = -0.5385909;}
      break;
    case 61:
      if (Lorb == 1) {kscal = -0.2307691;}
      else if (Lorb == 2) {kscal = -0.5309034;}
      break;
    case 62:
      if (Lorb == 1) {kscal = -0.2076921;}
      else if (Lorb == 2) {kscal = -0.5232158;}
      break;
    case 63:
      if (Lorb == 1) {kscal = -0.1846151;}
      else if (Lorb == 2) {kscal = -0.5155283;}
      break;
    case 64:
      if (Lorb == 1) {kscal = -0.1615381;}
      else if (Lorb == 2) {kscal = -0.5078408;}
      break;
    case 65:
      if (Lorb == 1) {kscal = -0.1384612;}
      else if (Lorb == 2) {kscal = -0.5001533;}
      break;
    case 66:
      if (Lorb == 1) {kscal = -0.1153842;}
      else if (Lorb == 2) {kscal = -0.4924658;}
      break;
    case 67:
      if (Lorb == 1) {kscal = -0.0923072;}
      else if (Lorb == 2) {kscal = -0.4847782;}
      break;
    case 68:
      if (Lorb == 1) {kscal = -0.0692302;}
      else if (Lorb == 2) {kscal = -0.4770907;}
      break;
    case 69:
      if (Lorb == 1) {kscal = -0.0461533;}
      else if (Lorb == 2) {kscal = -0.4694032;}
      break;
    case 70:
      if (Lorb == 1) {kscal = -0.0230763;}
      else if (Lorb == 2) {kscal = -0.4617157;}
      break;
    case 71:
      if (Lorb == 1) {kscal = 0.0000007;}
      else if (Lorb == 2) {kscal = -0.4540282;}
      break;
    case 72:
      if (Lorb == 1) {kscal = 0.1000000;}
      else if (Lorb == 2) {kscal = -0.4486165;}
      break;
    case 73:
      if (Lorb == 1) {kscal = 0.0500000;}
      else if (Lorb == 2) {kscal = -0.3394380;}
      break;
    case 74:
      if (Lorb == 1) {kscal = 0.3700000;}
      else if (Lorb == 2) {kscal = -0.3419199;}
      break;
    case 75:
      if (Lorb == 1) {kscal = -0.6000000;}
      else if (Lorb == 2) {kscal = 0.6586864;}
      break;
    case 76:
      if (Lorb == 1) {kscal = -0.6500000;}
      else if (Lorb == 2) {kscal = 0.1350223;}
      break;
    case 77:
      if (Lorb == 1) {kscal = -0.6500000;}
      else if (Lorb == 2) {kscal = -0.0977957;}
      break;
    case 78:
      if (Lorb == 1) {kscal = -0.6000000;}
      else if (Lorb == 2) {kscal = -0.0203212;}
      break;
    case 79:
      if (Lorb == 1) {kscal = -0.6000000;}
      else if (Lorb == 2) {kscal = 0.0614126;}
      break;
    case 80:
      if (Lorb == 1) {kscal = -0.5375121;}
      break;
    case 81:
      if (Lorb == 1) {kscal = -0.7133401;}
      break;
    case 82:
      if (Lorb == 1) {kscal = 0.7838251;}
      break;
    case 83:
      if (Lorb == 1) {kscal = -0.6000000;}
      break;
    case 84:
      if (Lorb == 1) {kscal = -0.8109155;}
      break;
    case 85:
      if (Lorb == 1) {kscal = -0.2532073;}
      else if (Lorb == 2) {kscal = 0.2500000;}
      break;
    case 86:
      if (Lorb == 1) {kscal = -0.0302388;}
      else if (Lorb == 2) {kscal = -0.2300000;}
      break;
  }
  return kscal + 1.0;
}
double HubbardDer(int atmnr) {
  //the gamma parameters, i.e., the Hubbard derivatives
  double gamma = 0.0;
  switch (atmnr) {
    case 1:
      gamma = 0.800000;
      break;
    case 2:
      gamma = 2.000000;
      break;
    case 3:
      gamma = 1.303821;
      break;
    case 4:
      gamma = 0.574239;
      break;
    case 5:
      gamma = 0.946104;
      break;
    case 6:
      gamma = 1.500000;
      break;
    case 7:
      gamma = -0.639780;
      break;
    case 8:
      gamma = -0.517134;
      break;
    case 9:
      gamma = 1.426212;
      break;
    case 10:
      gamma = 0.500000;
      break;
    case 11:
      gamma = 1.798727;
      break;
    case 12:
      gamma = 2.349164;
      break;
    case 13:
      gamma = 1.400000;
      break;
    case 14:
      gamma = 1.936289;
      break;
    case 15:
      gamma = 0.711291;
      break;
    case 16:
      gamma = -0.501722;
      break;
    case 17:
      gamma = 1.495483;
      break;
    case 18:
      gamma = -0.315455;
      break;
    case 19:
      gamma = 2.033085;
      break;
    case 20:
      gamma = 2.006898;
      break;
    case 21:
      gamma = 0.500000;
      break;
    case 22:
      gamma = 1.767268;
      break;
    case 23:
      gamma = 0.900000;
      break;
    case 24:
      gamma = 0.300000;
      break;
    case 25:
      gamma = 0.600000;
      break;
    case 26:
      gamma = -0.500000;
      break;
    case 27:
      gamma = 0.300000;
      break;
    case 28:
      gamma = -0.200000;
      break;
    case 29:
      gamma = 0.500000;
      break;
    case 30:
      gamma = 2.312896;
      break;
    case 31:
      gamma = 2.334269;
      break;
    case 32:
      gamma = -0.064775;
      break;
    case 33:
      gamma = 1.106041;
      break;
    case 34:
      gamma = 0.913725;
      break;
    case 35:
      gamma = 1.300000;
      break;
    case 36:
      gamma = 0.239815;
      break;
    case 37:
      gamma = 2.916203;
      break;
    case 38:
      gamma = 1.800000;
      break;
    case 39:
      gamma = 0.100000;
      break;
    case 40:
      gamma = 0.700000;
      break;
    case 41:
      gamma = 0.500000;
      break;
    case 42:
      gamma = 0.919928;
      break;
    case 43:
      gamma = 0.600000;
      break;
    case 44:
      gamma = -0.500000;
      break;
    case 45:
      gamma = 0.300000;
      break;
    case 46:
      gamma = 0.800000;
      break;
    case 47:
      gamma = 0.200000;
      break;
    case 48:
      gamma = 2.073217;
      break;
    case 49:
      gamma = 1.900000;
      break;
    case 50:
      gamma = -0.178396;
      break;
    case 51:
      gamma = 1.100000;
      break;
    case 52:
      gamma = 0.953683;
      break;
    case 53:
      gamma = 1.200000;
      break;
    case 54:
      gamma = -0.118925;
      break;
    case 55:
      gamma = 2.404185;
      break;
    case 56:
      gamma = 2.069097;
      break;
    case 57:
      gamma = 0.012793;
      break;
    case 58:
      gamma = -0.100000;
      break;
    case 59:
      gamma = -0.100002;
      break;
    case 60:
      gamma = -0.100004;
      break;
    case 61:
      gamma = -0.100006;
      break;
    case 62:
      gamma = -0.100008;
      break;
    case 63:
      gamma = -0.100010;
      break;
    case 64:
      gamma = -0.100012;
      break;
    case 65:
      gamma = -0.100013;
      break;
    case 66:
      gamma = -0.100015;
      break;
    case 67:
      gamma = -0.100017;
      break;
    case 68:
      gamma = -0.100019;
      break;
    case 69:
      gamma = -0.100021;
      break;
    case 70:
      gamma = -0.100023;
      break;
    case 71:
      gamma = -0.100025;
      break;
    case 72:
      gamma = -0.100000;
      break;
    case 73:
      gamma = 0.200000;
      break;
    case 74:
      gamma = -0.200000;
      break;
    case 75:
      gamma = 0.800000;
      break;
    case 76:
      gamma = 0.800000;
      break;
    case 77:
      gamma = -0.100000;
      break;
    case 78:
      gamma = 0.600000;
      break;
    case 79:
      gamma = 0.850000;
      break;
    case 80:
      gamma = -0.116312;
      break;
    case 81:
      gamma = -0.533933;
      break;
    case 82:
      gamma = 0.200000;
      break;
    case 83:
      gamma = -0.337508;
      break;
    case 84:
      gamma = 1.877978;
      break;
    case 85:
      gamma = 1.846485;
      break;
    case 86:
      gamma = 0.097834;
      break;
  }
  return 0.1*gamma;
}
double fXCmuA(int atmnr) {
  double fxcmua = 0.0;
  switch (atmnr) {
    case 1:
      fxcmua = 5.563889;
      break;
    case 2:
      fxcmua = -1.000000;
      break;
    case 3:
      fxcmua = -0.500000;
      break;
    case 4:
      fxcmua = -0.613341;
      break;
    case 5:
      fxcmua = -0.481186;
      break;
    case 6:
      fxcmua = -0.411674;
      break;
    case 7:
      fxcmua = 3.521273;
      break;
    case 8:
      fxcmua = -4.935670;
      break;
    case 9:
      fxcmua = -8.339183;
      break;
    case 10:
      fxcmua = 10.000000;
      break;
    case 11:
      fxcmua = 0.000000;
      break;
    case 12:
      fxcmua = -0.082005;
      break;
    case 13:
      fxcmua = 2.633341;
      break;
    case 14:
      fxcmua = -0.025750;
      break;
    case 15:
      fxcmua = 2.110225;
      break;
    case 16:
      fxcmua = -0.151117;
      break;
    case 17:
      fxcmua = -2.536958;
      break;
    case 18:
      fxcmua = -2.077329;
      break;
    case 19:
      fxcmua = -0.103383;
      break;
    case 20:
      fxcmua = -0.236675;
      break;
    case 21:
      fxcmua = -0.515177;
      break;
    case 22:
      fxcmua = -0.434506;
      break;
    case 23:
      fxcmua = -0.350000;
      break;
    case 24:
      fxcmua = 0.149669;
      break;
    case 25:
      fxcmua = -0.759168;
      break;
    case 26:
      fxcmua = 0.412929;
      break;
    case 27:
      fxcmua = -0.247938;
      break;
    case 28:
      fxcmua = -1.261887;
      break;
    case 29:
      fxcmua = -0.700000;
      break;
    case 30:
      fxcmua = -0.100000;
      break;
    case 31:
      fxcmua = 0.267219;
      break;
    case 32:
      fxcmua = 0.108460;
      break;
    case 33:
      fxcmua = -0.201294;
      break;
    case 34:
      fxcmua = -0.288648;
      break;
    case 35:
      fxcmua = -1.088586;
      break;
    case 36:
      fxcmua = -0.889357;
      break;
    case 37:
      fxcmua = -0.093328;
      break;
    case 38:
      fxcmua = -0.459925;
      break;
    case 39:
      fxcmua = -0.637291;
      break;
    case 40:
      fxcmua = -0.599615;
      break;
    case 41:
      fxcmua = -0.288729;
      break;
    case 42:
      fxcmua = 0.346327;
      break;
    case 43:
      fxcmua = -0.458416;
      break;
    case 44:
      fxcmua = -0.081922;
      break;
    case 45:
      fxcmua = 0.007016;
      break;
    case 46:
      fxcmua = -0.310361;
      break;
    case 47:
      fxcmua = -0.800314;
      break;
    case 48:
      fxcmua = -0.105364;
      break;
    case 49:
      fxcmua = 0.951079;
      break;
    case 50:
      fxcmua = 0.085029;
      break;
    case 51:
      fxcmua = -0.015519;
      break;
    case 52:
      fxcmua = -0.263414;
      break;
    case 53:
      fxcmua = -0.603648;
      break;
    case 54:
      fxcmua = -0.214447;
      break;
    case 55:
      fxcmua = -0.080000;
      break;
    case 56:
      fxcmua = -0.260000;
      break;
    case 57:
      fxcmua = -0.395198;
      break;
    case 58:
      fxcmua = -0.723806;
      break;
    case 59:
      fxcmua = -0.704819;
      break;
    case 60:
      fxcmua = -0.685832;
      break;
    case 61:
      fxcmua = -0.666845;
      break;
    case 62:
      fxcmua = -0.647858;
      break;
    case 63:
      fxcmua = -0.628871;
      break;
    case 64:
      fxcmua = -0.609884;
      break;
    case 65:
      fxcmua = -0.590897;
      break;
    case 66:
      fxcmua = -0.571910;
      break;
    case 67:
      fxcmua = -0.552923;
      break;
    case 68:
      fxcmua = -0.533936;
      break;
    case 69:
      fxcmua = -0.514949;
      break;
    case 70:
      fxcmua = -0.495961;
      break;
    case 71:
      fxcmua = -0.476974;
      break;
    case 72:
      fxcmua = -0.537685;
      break;
    case 73:
      fxcmua = -0.200343;
      break;
    case 74:
      fxcmua = 0.065886;
      break;
    case 75:
      fxcmua = -0.587636;
      break;
    case 76:
      fxcmua = -0.510090;
      break;
    case 77:
      fxcmua = -0.673822;
      break;
    case 78:
      fxcmua = -0.423684;
      break;
    case 79:
      fxcmua = 0.393418;
      break;
    case 80:
      fxcmua = -0.250000;
      break;
    case 81:
      fxcmua = 0.374018;
      break;
    case 82:
      fxcmua = 1.007016;
      break;
    case 83:
      fxcmua = -0.737252;
      break;
    case 84:
      fxcmua = -1.344854;
      break;
    case 85:
      fxcmua = -0.348123;
      break;
    case 86:
      fxcmua = -0.167597;
      break;
  }
  return fxcmua*0.01;
}
double fXCthetaA(int atmnr) {
  double fxcthetaa = 0.0;
  switch (atmnr) {
    case 1:
      fxcthetaa = 0.027431;
      break;
    case 2:
      fxcthetaa = -0.337528;
      break;
    case 3:
      fxcthetaa = 0.020000;
      break;
    case 4:
      fxcthetaa = -0.058586;
      break;
    case 5:
      fxcthetaa = -0.058228;
      break;
    case 6:
      fxcthetaa = 0.213583;
      break;
    case 7:
      fxcthetaa = 2.026786;
      break;
    case 8:
      fxcthetaa = -0.310828;
      break;
    case 9:
      fxcthetaa = -0.245955;
      break;
    case 10:
      fxcthetaa = -0.500000;
      break;
    case 11:
      fxcthetaa = 0.020000;
      break;
    case 12:
      fxcthetaa = -0.005516;
      break;
    case 13:
      fxcthetaa = -0.021887;
      break;
    case 14:
      fxcthetaa = -0.080000;
      break;
    case 15:
      fxcthetaa = 0.028679;
      break;
    case 16:
      fxcthetaa = 0.442859;
      break;
    case 17:
      fxcthetaa = 0.122783;
      break;
    case 18:
      fxcthetaa = -1.083404;
      break;
    case 19:
      fxcthetaa = 0.025000;
      break;
    case 20:
      fxcthetaa = 0.010000;
      break;
    case 21:
      fxcthetaa = -0.042004;
      break;
    case 22:
      fxcthetaa = 0.059660;
      break;
    case 23:
      fxcthetaa = 0.009764;
      break;
    case 24:
      fxcthetaa = 0.137744;
      break;
    case 25:
      fxcthetaa = 0.229903;
      break;
    case 26:
      fxcthetaa = 0.267734;
      break;
    case 27:
      fxcthetaa = 0.048237;
      break;
    case 28:
      fxcthetaa = -0.080000;
      break;
    case 29:
      fxcthetaa = -0.345631;
      break;
    case 30:
      fxcthetaa = 0.007658;
      break;
    case 31:
      fxcthetaa = -0.003616;
      break;
    case 32:
      fxcthetaa = -0.003589;
      break;
    case 33:
      fxcthetaa = 0.014149;
      break;
    case 34:
      fxcthetaa = 0.085728;
      break;
    case 35:
      fxcthetaa = 0.216935;
      break;
    case 36:
      fxcthetaa = -0.415024;
      break;
    case 37:
      fxcthetaa = 0.015000;
      break;
    case 38:
      fxcthetaa = 0.015000;
      break;
    case 39:
      fxcthetaa = 0.010460;
      break;
    case 40:
      fxcthetaa = -0.012944;
      break;
    case 41:
      fxcthetaa = 0.041491;
      break;
    case 42:
      fxcthetaa = 0.312549;
      break;
    case 43:
      fxcthetaa = 0.155242;
      break;
    case 44:
      fxcthetaa = 0.359228;
      break;
    case 45:
      fxcthetaa = 0.008570;
      break;
    case 46:
      fxcthetaa = -0.040485;
      break;
    case 47:
      fxcthetaa = -0.020810;
      break;
    case 48:
      fxcthetaa = 0.012250;
      break;
    case 49:
      fxcthetaa = -0.002031;
      break;
    case 50:
      fxcthetaa = -0.008243;
      break;
    case 51:
      fxcthetaa = -0.020630;
      break;
    case 52:
      fxcthetaa = -0.026864;
      break;
    case 53:
      fxcthetaa = 0.069660;
      break;
    case 54:
      fxcthetaa = -0.156200;
      break;
    case 55:
      fxcthetaa = 0.008000;
      break;
    case 56:
      fxcthetaa = 0.015000;
      break;
    case 57:
      fxcthetaa = -0.030000;
      break;
    case 58:
      fxcthetaa = -0.025000;
      break;
    case 59:
      fxcthetaa = -0.024615;
      break;
    case 60:
      fxcthetaa = -0.024231;
      break;
    case 61:
      fxcthetaa = -0.023846;
      break;
    case 62:
      fxcthetaa = -0.023462;
      break;
    case 63:
      fxcthetaa = -0.023077;
      break;
    case 64:
      fxcthetaa = -0.022692;
      break;
    case 65:
      fxcthetaa = -0.022308;
      break;
    case 66:
      fxcthetaa = -0.021923;
      break;
    case 67:
      fxcthetaa = -0.021538;
      break;
    case 68:
      fxcthetaa = -0.021154;
      break;
    case 69:
      fxcthetaa = -0.020769;
      break;
    case 70:
      fxcthetaa = -0.020385;
      break;
    case 71:
      fxcthetaa = -0.020000;
      break;
    case 72:
      fxcthetaa = -0.016478;
      break;
    case 73:
      fxcthetaa = 0.039599;
      break;
    case 74:
      fxcthetaa = 1.063309;
      break;
    case 75:
      fxcthetaa = 0.306870;
      break;
    case 76:
      fxcthetaa = 0.759049;
      break;
    case 77:
      fxcthetaa = 0.322935;
      break;
    case 78:
      fxcthetaa = 0.098019;
      break;
    case 79:
      fxcthetaa = -0.020320;
      break;
    case 80:
      fxcthetaa = -0.032901;
      break;
    case 81:
      fxcthetaa = -0.008506;
      break;
    case 82:
      fxcthetaa = -0.001670;
      break;
    case 83:
      fxcthetaa = 0.162529;
      break;
    case 84:
      fxcthetaa = 0.013818;
      break;
    case 85:
      fxcthetaa = 0.021624;
      break;
    case 86:
      fxcthetaa = -0.111556;
      break;
  }
  return fxcthetaa*0.01;
}
double AtomOccupancy(int atmnr, int ishell) {
  double occ = 0.0;
  switch (atmnr) {
    case 1:
      if (ishell == 0) {occ = 1.0;}
      break;
    case 2:
      if (ishell == 0) {occ = 2.0;}
      break;
    case 3:
      if (ishell == 0) {occ = 1.0;}
      break;
    case 4:
      if (ishell == 0) {occ = 2.0;}
      break;
    case 5:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 1.0;}
      break;
    case 6:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 3.0;}
      break;
    case 7:
      if (ishell == 0) {occ = 1.5;}
      else if (ishell == 1) {occ = 3.5;}
      break;
    case 8:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 4.0;}
      break;
    case 9:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 5.0;}
      break;
    case 10:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 6.0;}
      break;
    case 11:
      if (ishell == 0) {occ = 1.0;}
      break;
    case 12:
      if (ishell == 0) {occ = 2.0;}
      break;
    case 13:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 1.0;}
      break;
    case 14:
      if (ishell == 0) {occ = 1.5;}
      else if (ishell == 1) {occ = 2.5;}
      break;
    case 15:
      if (ishell == 0) {occ = 1.5;}
      else if (ishell == 1) {occ = 3.5;}
      break;
    case 16:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 4.0;}
      break;
    case 17:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 5.0;}
      break;
    case 18:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 6.0;}
      break;
    case 19:
      if (ishell == 0) {occ = 1.0;}
      break;
    case 20:
      if (ishell == 0) {occ = 2.0;}
      break;
    case 21:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 1.0;}
      break;
    case 22:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 2.0;}
      break;
    case 23:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 3.0;}
      break;
    case 24:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 4.0;}
      break;
    case 25:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 5.0;}
      break;
    case 26:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 6.0;}
      break;
    case 27:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 7.0;}
      break;
    case 28:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 8.0;}
      break;
    case 29:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 2) {occ = 10.0;}
      break;
    case 30:
      if (ishell == 0) {occ = 2.0;}
      break;
    case 31:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 1.0;}
      break;
    case 32:
      if (ishell == 0) {occ = 1.5;}
      else if (ishell == 1) {occ = 2.5;}
      break;
    case 33:
      if (ishell == 0) {occ = 1.5;}
      else if (ishell == 1) {occ = 3.5;}
      break;
    case 34:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 4.0;}
      break;
    case 35:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 5.0;}
      break;
    case 36:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 6.0;}
      break;
    case 37:
      if (ishell == 0) {occ = 1.0;}
      break;
    case 38:
      if (ishell == 0) {occ = 2.0;}
      break;
    case 39:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 1.0;}
      break;
    case 40:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 2.0;}
      break;
    case 41:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 3.0;}
      break;
    case 42:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 4.0;}
      break;
    case 43:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 5.0;}
      break;
    case 44:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 6.0;}
      break;
    case 45:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 7.0;}
      break;
    case 46:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 8.0;}
      break;
    case 47:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 2) {occ = 10.0;}
      break;
    case 48:
      if (ishell == 0) {occ = 2.0;}
      break;
    case 49:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 1.0;}
      break;
    case 50:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 2.0;}
      break;
    case 51:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 3.0;}
      break;
    case 52:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 4.0;}
      break;
    case 53:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 5.0;}
      break;
    case 54:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 6.0;}
      break;
    case 55:
      if (ishell == 0) {occ = 1.0;}
      break;
    case 56:
      if (ishell == 0) {occ = 2.0;}
      break;
    case 57:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 1.0;}
      break;
    case 58:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 1.0;}
      break;
    case 59:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 1.0;}
      break;
    case 60:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 1.0;}
      break;
    case 61:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 1.0;}
      break;
    case 62:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 1.0;}
      break;
    case 63:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 1.0;}
      break;
    case 64:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 1.0;}
      break;
    case 65:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 1.0;}
      break;
    case 66:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 1.0;}
      break;
    case 67:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 1.0;}
      break;
    case 68:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 1.0;}
      break;
    case 69:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 1.0;}
      break;
    case 70:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 1.0;}
      break;
    case 71:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 1.0;}
      break;
    case 72:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 2.0;}
      break;
    case 73:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 3.0;}
      break;
    case 74:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 4.0;}
      break;
    case 75:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 5.0;}
      break;
    case 76:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 6.0;}
      break;
    case 77:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 7.0;}
      break;
    case 78:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 1) {occ = 1.0;}
      else if (ishell == 2) {occ = 8.0;}
      break;
    case 79:
      if (ishell == 0) {occ = 1.0;}
      else if (ishell == 2) {occ = 10.0;}
      break;
    case 80:
      if (ishell == 0) {occ = 2.0;}
      break;
    case 81:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 1.0;}
      break;
    case 82:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 2.0;}
      break;
    case 83:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 3.0;}
      break;
    case 84:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 4.0;}
      break;
    case 85:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 5.0;}
      break;
    case 86:
      if (ishell == 0) {occ = 2.0;}
      else if (ishell == 1) {occ = 6.0;}
      break;
  }
  return occ;
}
double MultipoleRadius(int atmnr) {
  //atomic radii used in multipole electrostatics (in Bohr)
  double multrad = 0.0;
  switch (atmnr) {
    case 1:
      multrad = 1.4;
      break;
    case 2:
      multrad = 3.0;
      break;
    case 3:
      multrad = 5.0;
      break;
    case 4:
      multrad = 5.0;
      break;
    case 5:
      multrad = 5.0;
      break;
    case 6:
      multrad = 3.0;
      break;
    case 7:
      multrad = 1.9;
      break;
    case 8:
      multrad = 1.8;
      break;
    case 9:
      multrad = 2.4;
      break;
    case 10:
      multrad = 5.0;
      break;
    case 11:
      multrad = 5.0;
      break;
    case 12:
      multrad = 5.0;
      break;
    case 13:
      multrad = 5.0;
      break;
    case 14:
      multrad = 3.9;
      break;
    case 15:
      multrad = 2.1;
      break;
    case 16:
      multrad = 3.1;
      break;
    case 17:
      multrad = 2.5;
      break;
    case 18:
      multrad = 5.0;
      break;
    case 19:
      multrad = 5.0;
      break;
    case 20:
      multrad = 5.0;
      break;
    case 21:
      multrad = 5.0;
      break;
    case 22:
      multrad = 5.0;
      break;
    case 23:
      multrad = 5.0;
      break;
    case 24:
      multrad = 5.0;
      break;
    case 25:
      multrad = 5.0;
      break;
    case 26:
      multrad = 5.0;
      break;
    case 27:
      multrad = 5.0;
      break;
    case 28:
      multrad = 5.0;
      break;
    case 29:
      multrad = 5.0;
      break;
    case 30:
      multrad = 5.0;
      break;
    case 31:
      multrad = 5.0;
      break;
    case 32:
      multrad = 5.0;
      break;
    case 33:
      multrad = 5.0;
      break;
    case 34:
      multrad = 3.9;
      break;
    case 35:
      multrad = 4.0;
      break;
    case 36:
      multrad = 5.0;
      break;
    case 37:
      multrad = 5.0;
      break;
    case 38:
      multrad = 5.0;
      break;
    case 39:
      multrad = 5.0;
      break;
    case 40:
      multrad = 5.0;
      break;
    case 41:
      multrad = 5.0;
      break;
    case 42:
      multrad = 5.0;
      break;
    case 43:
      multrad = 5.0;
      break;
    case 44:
      multrad = 5.0;
      break;
    case 45:
      multrad = 5.0;
      break;
    case 46:
      multrad = 5.0;
      break;
    case 47:
      multrad = 5.0;
      break;
    case 48:
      multrad = 5.0;
      break;
    case 49:
      multrad = 5.0;
      break;
    case 50:
      multrad = 5.0;
      break;
    case 51:
      multrad = 5.0;
      break;
    case 52:
      multrad = 5.0;
      break;
    case 53:
      multrad = 5.0;
      break;
    case 54:
      multrad = 5.0;
      break;
    case 55:
      multrad = 5.0;
      break;
    case 56:
      multrad = 5.0;
      break;
    case 57:
      multrad = 5.0;
      break;
    case 58:
      multrad = 5.0;
      break;
    case 59:
      multrad = 5.0;
      break;
    case 60:
      multrad = 5.0;
      break;
    case 61:
      multrad = 5.0;
      break;
    case 62:
      multrad = 5.0;
      break;
    case 63:
      multrad = 5.0;
      break;
    case 64:
      multrad = 5.0;
      break;
    case 65:
      multrad = 5.0;
      break;
    case 66:
      multrad = 5.0;
      break;
    case 67:
      multrad = 5.0;
      break;
    case 68:
      multrad = 5.0;
      break;
    case 69:
      multrad = 5.0;
      break;
    case 70:
      multrad = 5.0;
      break;
    case 71:
      multrad = 5.0;
      break;
    case 72:
      multrad = 5.0;
      break;
    case 73:
      multrad = 5.0;
      break;
    case 74:
      multrad = 5.0;
      break;
    case 75:
      multrad = 5.0;
      break;
    case 76:
      multrad = 5.0;
      break;
    case 77:
      multrad = 5.0;
      break;
    case 78:
      multrad = 5.0;
      break;
    case 79:
      multrad = 5.0;
      break;
    case 80:
      multrad = 5.0;
      break;
    case 81:
      multrad = 5.0;
      break;
    case 82:
      multrad = 5.0;
      break;
    case 83:
      multrad = 5.0;
      break;
    case 84:
      multrad = 5.0;
      break;
    case 85:
      multrad = 5.0;
      break;
    case 86:
      multrad = 5.0;
      break;
    case 87:
      multrad = 5.0;
      break;
    case 88:
      multrad = 5.0;
      break;
    case 89:
      multrad = 5.0;
      break;
    case 90:
      multrad = 5.0;
      break;
    case 91:
      multrad = 5.0;
      break;
    case 92:
      multrad = 5.0;
      break;
    case 93:
      multrad = 5.0;
      break;
    case 94:
      multrad = 5.0;
      break;
    case 95:
      multrad = 5.0;
      break;
    case 96:
      multrad = 5.0;
      break;
    case 97:
      multrad = 5.0;
      break;
    case 98:
      multrad = 5.0;
      break;
    case 99:
      multrad = 5.0;
      break;
    case 100:
      multrad = 5.0;
      break;
    case 101:
      multrad = 5.0;
      break;
    case 102:
      multrad = 5.0;
      break;
    case 103:
      multrad = 5.0;
      break;
    case 104:
      multrad = 5.0;
      break;
    case 105:
      multrad = 5.0;
      break;
    case 106:
      multrad = 5.0;
      break;
    case 107:
      multrad = 5.0;
      break;
    case 108:
      multrad = 5.0;
      break;
    case 109:
      multrad = 5.0;
      break;
    case 110:
      multrad = 5.0;
      break;
    case 111:
      multrad = 5.0;
      break;
    case 112:
      multrad = 5.0;
      break;
    case 113:
      multrad = 5.0;
      break;
    case 114:
      multrad = 5.0;
      break;
    case 115:
      multrad = 5.0;
      break;
    case 116:
      multrad = 5.0;
      break;
    case 117:
      multrad = 5.0;
      break;
    case 118:
      multrad = 5.0;
      break;
  }
  return multrad;
}
double Nval(int atmnr) {
  //atomic radii used in multipole electrostatics (in Bohr)
  double nval = 0.0;
  switch (atmnr) {
    case 1:
      nval = 1.0;
      break;
    case 2:
      nval = 1.0;
      break;
    case 3:
      nval = 1.0;
      break;
    case 4:
      nval = 2.0;
      break;
    case 5:
      nval = 3.0;
      break;
    case 6:
      nval = 3.0;
      break;
    case 7:
      nval = 3.0;
      break;
    case 8:
      nval = 2.0;
      break;
    case 9:
      nval = 1.0;
      break;
    case 10:
      nval = 1.0;
      break;
    case 11:
      nval = 1.0;
      break;
    case 12:
      nval = 2.0;
      break;
    case 13:
      nval = 3.0;
      break;
    case 14:
      nval = 3.0;
      break;
    case 15:
      nval = 3.0;
      break;
    case 16:
      nval = 3.0;
      break;
    case 17:
      nval = 1.0;
      break;
    case 18:
      nval = 1.0;
      break;
    case 19:
      nval = 1.0;
      break;
    case 20:
      nval = 2.0;
      break;
    case 21:
      nval = 4.0;
      break;
    case 22:
      nval = 4.0;
      break;
    case 23:
      nval = 6.0;
      break;
    case 24:
      nval = 6.0;
      break;
    case 25:
      nval = 6.0;
      break;
    case 26:
      nval = 6.0;
      break;
    case 27:
      nval = 6.0;
      break;
    case 28:
      nval = 4.0;
      break;
    case 29:
      nval = 4.0;
      break;
    case 30:
      nval = 2.0;
      break;
    case 31:
      nval = 3.0;
      break;
    case 32:
      nval = 3.0;
      break;
    case 33:
      nval = 3.0;
      break;
    case 34:
      nval = 3.0;
      break;
    case 35:
      nval = 1.0;
      break;
    case 36:
      nval = 1.0;
      break;
    case 37:
      nval = 1.0;
      break;
    case 38:
      nval = 2.0;
      break;
    case 39:
      nval = 4.0;
      break;
    case 40:
      nval = 4.0;
      break;
    case 41:
      nval = 6.0;
      break;
    case 42:
      nval = 6.0;
      break;
    case 43:
      nval = 6.0;
      break;
    case 44:
      nval = 6.0;
      break;
    case 45:
      nval = 6.0;
      break;
    case 46:
      nval = 4.0;
      break;
    case 47:
      nval = 4.0;
      break;
    case 48:
      nval = 2.0;
      break;
    case 49:
      nval = 3.0;
      break;
    case 50:
      nval = 3.0;
      break;
    case 51:
      nval = 3.0;
      break;
    case 52:
      nval = 3.0;
      break;
    case 53:
      nval = 1.0;
      break;
    case 54:
      nval = 1.0;
      break;
    case 55:
      nval = 1.0;
      break;
    case 56:
      nval = 2.0;
      break;
    case 57:
      nval = 4.0;
      break;
    case 58:
      nval = 6.0;
      break;
    case 59:
      nval = 6.0;
      break;
    case 60:
      nval = 6.0;
      break;
    case 61:
      nval = 6.0;
      break;
    case 62:
      nval = 6.0;
      break;
    case 63:
      nval = 6.0;
      break;
    case 64:
      nval = 6.0;
      break;
    case 65:
      nval = 6.0;
      break;
    case 66:
      nval = 6.0;
      break;
    case 67:
      nval = 6.0;
      break;
    case 68:
      nval = 6.0;
      break;
    case 69:
      nval = 6.0;
      break;
    case 70:
      nval = 6.0;
      break;
    case 71:
      nval = 6.0;
      break;
    case 72:
      nval = 4.0;
      break;
    case 73:
      nval = 6.0;
      break;
    case 74:
      nval = 6.0;
      break;
    case 75:
      nval = 6.0;
      break;
    case 76:
      nval = 6.0;
      break;
    case 77:
      nval = 6.0;
      break;
    case 78:
      nval = 4.0;
      break;
    case 79:
      nval = 4.0;
      break;
    case 80:
      nval = 2.0;
      break;
    case 81:
      nval = 3.0;
      break;
    case 82:
      nval = 3.0;
      break;
    case 83:
      nval = 3.0;
      break;
    case 84:
      nval = 3.0;
      break;
    case 85:
      nval = 1.0;
      break;
    case 86:
      nval = 1.0;
      break;
  }
  return nval;
}

#endif //_GFN2_XTB_Parameters_
