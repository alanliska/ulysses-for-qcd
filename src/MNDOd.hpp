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

#ifndef _MNDOD_
#define _MNDOD_
#include <string>
#include <sstream>
#include <vector>
#include "MNDO.hpp"

//descrition:
//NDDO methods with spd bases

class MNDOd: public MNDO {
  //implementation of Thiel's and Voityuk's MNDOd
  //W. Thiel, A. A. Voityuk, J. Phys. Chem., 100, 616, 1996
  //bX quantities are "barred" tensors, which are used only for the open-shell case
public:
  MNDOd(BSet _bset, Molecule _mol, std::string _openclosed = "0", std::string _corecorrection = "0"): MNDO(_bset,_mol,_openclosed,_corecorrection) {
    Dd.resize(2);
    intn.resize(31);
  }
  ~MNDOd() {}
  //getters
  virtual std::string Type() {return "MNDOd";}
  bool AtomWithDOrbitals(size_t atomicnr) {
    //function that determines whether atom has d orbitals
    return (((atomicnr > 12)&&(atomicnr < 18))||((atomicnr > 20)&&(atomicnr < 30))||(atomicnr == 33)||(atomicnr == 35)||((atomicnr > 38)&&(atomicnr < 48))||(atomicnr == 51)||(atomicnr == 53)||(atomicnr == 57)||((atomicnr > 70)&&(atomicnr < 80)));
  }
  double SlaterCondonRadialIntegral(int atmnr, int L, int II, int JJ, double threshold = 1e-7) {
    //function returning Slater Condon radial integrals used to calculate 1-Center eris
    //II and JJ refer to the orbitals in the bra and ket respectively
    //they take the values
    //        1 -> ss; 2 -> sp; 3 -> pp; 4 -> sd; 5 -> pd; 6 -> dd
    //L is the radial number of the Slater Condon integral
    //these integrals translate into
    //    L,II,JJ = 0,6,6 -> F0dd;  2,6,6 -> F2dd;  4,6,6 -> F4dd;  0,1,6 -> F0sd;  2,4,4 -> G2sd;
    //            = 0,3,6 -> F0pd;  2,3,6 -> F2pd;  1,5,5 -> G1pd;  3,5,5 -> G3pd;
    //more combinations possible though
    bool calculate = true;
    double scint = 0.0;
    //check whether there are tabulated values so that calculation is skipped
    if ((L == 0)&&(II == 1)&&(JJ == 6)) {
      scint = F0sd(atmnr);
      if (fabs(scint) > threshold) {calculate = false;}
    }
    else if ((L == 2)&&(II == 4)&&(JJ == 4)) {
      scint = G2sd(atmnr);
      if (fabs(scint) > threshold) {calculate = false;}
    }
    if (calculate) {
      //get principal quantum numbers of orbitals in the integral (mu,nu|lambda,sigma)
      size_t l_1 = 0;                //s
      size_t l_2 = 0;                //s
      if (II == 2) {l_2 = 1;}        //sp
      else if (II == 3) {            //pp
        l_1 = 1;
        l_2 = 1;
      }
      else if (II == 4) {l_2 = 2;}   //sd
      else if (II == 5) {            //pd
        l_1 = 1;
        l_2 = 2;
      }
      else if (II == 6) {            //dd
        l_1 = 2;
        l_2 = 2;
      }
      double n_mu = PrincipalQuantNumb(atmnr,l_1);
      double n_nu = PrincipalQuantNumb(atmnr,l_2);
      double z_mu = zn(atmnr,l_1);
      double z_nu = zn(atmnr,l_2);
      l_1 = 0;                       //s
      l_2 = 0;                       //s
      if (JJ == 2) {l_2 = 1;}        //sp
      else if (JJ == 3) {            //pp
        l_1 = 1;
        l_2 = 1;
      }
      else if (JJ == 4) {l_2 = 2;}   //sd
      else if (JJ == 5) {            //pd
        l_1 = 1;
        l_2 = 2;
      }
      else if (JJ == 6) {            //dd
        l_1 = 2;
        l_2 = 2;
      }
      double n_lambda = PrincipalQuantNumb(atmnr,l_1);
      double n_sigma = PrincipalQuantNumb(atmnr,l_2);
      double z_lambda = zn(atmnr,l_1);
      double z_sigma = zn(atmnr,l_2);
      int nmn = int(n_mu + n_nu);
      int nls = int(n_lambda + n_sigma);
      int nmnls = nmn + nls;
      double lzm = log(z_mu);
      double lzn = log(z_nu);
      double lzl = log(z_lambda);
      double lzs = log(z_sigma);
      double zls = z_lambda + z_sigma;
      double zmn = z_mu + z_nu;
      double zmnls = zls + zmn;
      double lzmnls = log(zmnls);
      double lzls = log(zls);
      double lzmn = log(zmn);
      double fac_a = Factorial(nmnls - 1)/sqrt(Factorial(2*int(n_mu))*Factorial(2*int(n_nu))*Factorial(2*int(n_lambda))*Factorial(2*int(n_sigma)));
      double fac0 = fac_a*exp(n_mu*lzm + n_nu*lzn + n_lambda*lzl + n_sigma*lzs + 0.5*(lzm + lzn + lzl + lzs) + log(2.0)*double(nmnls + 2) - lzmnls*double(nmnls));
      double aux = 1.0/zmnls;
      double fac1 = 0.0;
      double fac2 = 0.0;
      double fac3;
      for (size_t i = 0; i < nls - L; i++) {
        aux = aux*zmnls/zls;
        fac1 += aux*(dCombinations(nls - L - 1,i) - dCombinations(nls + L,i))/dCombinations(nmnls - 1,i);
      }
      for (size_t i = nls - L; i < nls + L + 1; i++) {
        aux = aux*zmnls/zls;
        fac2 += aux*dCombinations(nls + L,i)/dCombinations(nmnls - 1,i);
      }
      fac3 = exp(lzmnls*nmnls - lzls*(nls + L + 1) - lzmn*(nmn - L))/dCombinations(nmnls - 1,nls + L);
      scint = fac0*(fac1 - fac2 + fac3);
    }
    return scint;
  }
  double d_eri2Center(int L1, int M1, int L2, int M2, int L3, int M3, int L4, int M4, double RAB, std::vector<double> D,size_t atmnrA,size_t atmnrB, int core = 0) {
    //function that calculates 2 electron repulsion integrals as in the NDDO approximation MNDO/d
    //AOs contains AO information on the orbitals involved in the integral, which is valid for s,p,d basis/atoms
    //RAB is distance between centers
    //Thiel, Voityuk, Theor. Chim. Acta, 81, 391, 1992
    double eri = 0.0;
    double rhoA;
    double rhoB;
    int lA;
    int lB;
    int mA;
    int mB;
    int l_pairA;
    int l_pairB;
    chg1 = NonZeroMultipoleMom_d(L1,M1,L2,M2);
    chg2 = NonZeroMultipoleMom_d(L3,M3,L4,M4);
    for (size_t mp1 = 0; mp1 < chg1.size(); ++mp1) {
      lA = chg1[mp1][0];
      mA = chg1[mp1][1];
      for (size_t mp2 = 0; mp2 < chg2.size(); ++mp2) {
        lB = chg2[mp2][0];
        mB = chg2[mp2][1];
        if (((abs(mA) % 2 == 0)&&(abs(mB) % 2 == 0))||((abs(mA) % 2 != 0)&&(abs(mB) % 2 != 0))) {
          //even-even or odd-odd only
          if (abs(mA) != abs(mB)) {
            mA = abs(mA);
            mB = abs(mB);
          }
          if ((L1 == 0)&&(L2 == 0)) {l_pairA = 0;}
          else if (((L1 == 0)&&(L2 == 1))||((L1 == 1)&&(L2 == 0))) {l_pairA = 1;}
          else if ((L1 == 1)&&(L2 == 1)) {
            l_pairA = 2;
            if (lA == 0) {l_pairA = 0;}
          }
          else if (((L1 == 0)&&(L2 == 2))||((L1 == 2)&&(L2 == 0))) {l_pairA = 3;}
          else if (((L1 == 1)&&(L2 == 2))||((L1 == 2)&&(L2 == 1))) {l_pairA = 4;}
          else if ((L1 == 2)&&(L2 == 2)) {
            l_pairA = 5;
            if (lA == 2) {l_pairA = 6;}
            //if ((lA == 2)&&(atmnrB < 11)) {l_pairA = 4;}
          }
          if ((L3 == 0)&&(L4 == 0)) {l_pairB = 0;}
          else if (((L3 == 0)&&(L4 == 1))||((L3 == 1)&&(L4 == 0))) {l_pairB = 1;}
          else if ((L3 == 1)&&(L4 == 1)) {
            l_pairB = 2;
            if (lB == 0) {l_pairB = 0;}
          }
          else if (((L3 == 0)&&(L4 == 2))||((L3 == 2)&&(L4 == 0))) {l_pairB = 3;}
          else if (((L3 == 1)&&(L4 == 2))||((L3 == 2)&&(L4 == 1))) {l_pairB = 4;}
          else if ((L3 == 2)&&(L4 == 2)) {
            l_pairB = 5;
            if (lB == 2) {l_pairB = 6;}
            //if ((lB == 2)&&(atmnrA < 11)) {l_pairB = 4;}
          }
          rhoA = rho(atmnrA,l_pairA);
          rhoB = rho(atmnrB,l_pairB);
          if (core > 0) {rhoB = rhocore(atmnrB);}                        //in this case use the correct core rho for the ss atom
          if ((core > 1)||(core < 0)) {rhoA = rhocore(atmnrA);}          //in this case use the correct core rho for the ss atom
          eri += d_se_multipole(lA,mA,lB,mB,RAB,rhoA + rhoB,D)*clm_munu(L1,M1,L2,M2,lA,mA)*clm_munu(L3,M3,L4,M4,lB,mB);
        }
      }
    }
    return eri;
  }
  double d_eri2Center_dR(int L1, int M1, int L2, int M2, int L3, int M3, int L4, int M4, double RAB, std::vector<double> D,size_t atmnrA,size_t atmnrB, int core = 0) {
    //function that calculates the first-derivatives of 2 electron repulsion integrals as in the NDDO approximation MNDO/d with respect to internuclear distance
    //AOs contains AO information on the orbitals involved in the integral, which is valid for s,p,d basis/atoms
    //RAB is distance between centers
    double eri = 0.0;
    double rhoA;
    double rhoB;
    int lA;
    int lB;
    int mA;
    int mB;
    int l_pairA;
    int l_pairB;
    chg1 = NonZeroMultipoleMom_d(L1,M1,L2,M2);
    chg2 = NonZeroMultipoleMom_d(L3,M3,L4,M4);
    for (size_t mp1 = 0; mp1 < chg1.size(); ++mp1) {
      lA = chg1[mp1][0];
      mA = chg1[mp1][1];
      for (size_t mp2 = 0; mp2 < chg2.size(); ++mp2) {
        lB = chg2[mp2][0];
        mB = chg2[mp2][1];
        if (((abs(mA) % 2 == 0)&&(abs(mB) % 2 == 0))||((abs(mA) % 2 != 0)&&(abs(mB) % 2 != 0))) {
          //even-even or odd-odd only
          if (abs(mA) != abs(mB)) {
            mA = abs(mA);
            mB = abs(mB);
          }
          if ((L1 == 0)&&(L2 == 0)) {l_pairA = 0;}
          else if (((L1 == 0)&&(L2 == 1))||((L1 == 1)&&(L2 == 0))) {l_pairA = 1;}
          else if ((L1 == 1)&&(L2 == 1)) {
            l_pairA = 2;
            if (lA == 0) {l_pairA = 0;}
          }
          else if (((L1 == 0)&&(L2 == 2))||((L1 == 2)&&(L2 == 0))) {l_pairA = 3;}
          else if (((L1 == 1)&&(L2 == 2))||((L1 == 2)&&(L2 == 1))) {l_pairA = 4;}
          else if ((L1 == 2)&&(L2 == 2)) {
            l_pairA = 5;
            if (lA == 2) {l_pairA = 6;}
          }
          if ((L3 == 0)&&(L4 == 0)) {l_pairB = 0;}
          else if (((L3 == 0)&&(L4 == 1))||((L3 == 1)&&(L4 == 0))) {l_pairB = 1;}
          else if ((L3 == 1)&&(L4 == 1)) {
            l_pairB = 2;
            if (lB == 0) {l_pairB = 0;}
          }
          else if (((L3 == 0)&&(L4 == 2))||((L3 == 2)&&(L4 == 0))) {l_pairB = 3;}
          else if (((L3 == 1)&&(L4 == 2))||((L3 == 2)&&(L4 == 1))) {l_pairB = 4;}
          else if ((L3 == 2)&&(L4 == 2)) {
            l_pairB = 5;
            if (lB == 2) {l_pairB = 6;}
          }
          rhoA = rho(atmnrA,l_pairA);
          rhoB = rho(atmnrB,l_pairB);
          if (core > 0) {rhoB = rhocore(atmnrB);}                        //in this case use the correct core rho for the ss atom
          if ((core > 1)||(core < 0)) {rhoA = rhocore(atmnrA);}          //in this case use the correct core rho for the ss atom
          eri += d_se_multipole_dR(lA,mA,lB,mB,RAB,rhoA + rhoB,D)*clm_munu(L1,M1,L2,M2,lA,mA)*clm_munu(L3,M3,L4,M4,lB,mB);
        }
      }
    }
    return eri;
  }
  double d_eri2Center_dR2(int L1, int M1, int L2, int M2, int L3, int M3, int L4, int M4, double RAB, std::vector<double> D,size_t atmnrA,size_t atmnrB, int core = 0) {
    //function that calculates the second-derivatives of 2 electron repulsion integrals as in the NDDO approximation MNDO/d with respect to internuclear distance
    //AOs contains AO information on the orbitals involved in the integral, which is valid for s,p,d basis/atoms
    //RAB is distance between centers
    double eri = 0.0;
    double rhoA;
    double rhoB;
    int lA;
    int lB;
    int mA;
    int mB;
    int l_pairA;
    int l_pairB;
    chg1 = NonZeroMultipoleMom_d(L1,M1,L2,M2);
    chg2 = NonZeroMultipoleMom_d(L3,M3,L4,M4);
    for (size_t mp1 = 0; mp1 < chg1.size(); ++mp1) {
      lA = chg1[mp1][0];
      mA = chg1[mp1][1];
      for (size_t mp2 = 0; mp2 < chg2.size(); ++mp2) {
        lB = chg2[mp2][0];
        mB = chg2[mp2][1];
        if (((abs(mA) % 2 == 0)&&(abs(mB) % 2 == 0))||((abs(mA) % 2 != 0)&&(abs(mB) % 2 != 0))) {
          //even-even or odd-odd only
          if (abs(mA) != abs(mB)) {
            mA = abs(mA);
            mB = abs(mB);
          }
          if ((L1 == 0)&&(L2 == 0)) {l_pairA = 0;}
          else if (((L1 == 0)&&(L2 == 1))||((L1 == 1)&&(L2 == 0))) {l_pairA = 1;}
          else if ((L1 == 1)&&(L2 == 1)) {
            l_pairA = 2;
            if (lA == 0) {l_pairA = 0;}
          }
          else if (((L1 == 0)&&(L2 == 2))||((L1 == 2)&&(L2 == 0))) {l_pairA = 3;}
          else if (((L1 == 1)&&(L2 == 2))||((L1 == 2)&&(L2 == 1))) {l_pairA = 4;}
          else if ((L1 == 2)&&(L2 == 2)) {
            l_pairA = 5;
            if (lA == 2) {l_pairA = 6;}
          }
          if ((L3 == 0)&&(L4 == 0)) {l_pairB = 0;}
          else if (((L3 == 0)&&(L4 == 1))||((L3 == 1)&&(L4 == 0))) {l_pairB = 1;}
          else if ((L3 == 1)&&(L4 == 1)) {
            l_pairB = 2;
            if (lB == 0) {l_pairB = 0;}
          }
          else if (((L3 == 0)&&(L4 == 2))||((L3 == 2)&&(L4 == 0))) {l_pairB = 3;}
          else if (((L3 == 1)&&(L4 == 2))||((L3 == 2)&&(L4 == 1))) {l_pairB = 4;}
          else if ((L3 == 2)&&(L4 == 2)) {
            l_pairB = 5;
            if (lB == 2) {l_pairB = 6;}
          }
          rhoA = rho(atmnrA,l_pairA);
          rhoB = rho(atmnrB,l_pairB);
          if (core > 0) {rhoB = rhocore(atmnrB);}                        //in this case use the correct core rho for the ss atom
          if ((core > 1)||(core < 0)) {rhoA = rhocore(atmnrA);}          //in this case use the correct core rho for the ss atom
          eri += d_se_multipole_dR2(lA,mA,lB,mB,RAB,rhoA + rhoB,D)*clm_munu(L1,M1,L2,M2,lA,mA)*clm_munu(L3,M3,L4,M4,lB,mB);
        }
      }
    }
    return eri;
  }
  //other auxiliary functions
  virtual void checkAtoms() {
    //function that checks whether the atoms in the molecule are defined in MNDOd
    bool def;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      def = false;
      if (atoms[iatm] == 1) {def = true;}                                                       //H
      else if ((atoms[iatm] > 2)&&(atoms[iatm] < 10)) {def = true;}                             //Li,Be,B,C,N,O,F
      else if ((atoms[iatm] > 10)&&(atoms[iatm] < 18)) {def = true;}                            //Na,Mg,Al,Si,P,S,Cl
      else if ((atoms[iatm] == 30)||(atoms[iatm] == 32)||(atoms[iatm] == 48)) {def = true;}     //Zn,Ge,Cd
      else if ((atoms[iatm] == 35)||(atoms[iatm] == 53)) {def = true;}                          //Br,I
      else if ((atoms[iatm] == 50)||(atoms[iatm] == 80)||(atoms[iatm] == 82)) {def = true;}     //Sn,Hg,Pb
      if (!def) {throw("ERROR: MNDO.hpp: MNDO: checkAtoms(): atom not fully specified for MNDO-theory");}
    }
  }
  //estimate parameters for heat of formation
  double calcEISOL(int atmnr) {
    //function that calculates the theoretical heats of formation for atoms
    //this function is not used directly by the methods, its purpose is the generation of the parameters that are then tabulated
    //integrals
    double Uss = UlX(atmnr,0);
    double Upp = UlX(atmnr,1);
    double Udd = UlX(atmnr,2);
    double Gss = eri1Center(atmnr,0,0);
    double Gpp = eri1Center(atmnr,2,2);
    double Gsp = eri1Center(atmnr,0,2);
    double Gp2 = eri1Center(atmnr,2,-2);
    double Hsp = eri1Center(atmnr,1,1);
    if (Hsp < 1.0e-7) {Hsp = 1.0e-7;}
    //Take into account constraints on the values of the integrals
    double Hpp = 0.5*(Gpp - Gp2);
    if (Hpp < 0.1) {Hpp = 0.1;}
    //counting how often integrals show up
    int nrs = NrS(atmnr);
    int nrp = NrP(atmnr);
    int nrd = NrD(atmnr);
    int ngss = (nrs - 1)*int(nrs > 1);                     //number of integrals of type <ss|ss>
    int ngsp = nrs*nrp;                                    //number of integrals of type <ss|pp>
    int iaux = std::min(nrp,6 - nrp);
    int npp = -iaux*(iaux - 1)/4;
    int ngp2 = nrp*(nrp - 1)/2 - npp;                      //number of integrals of type <pp|pp>; note that hpp is replaced by 0.5*(gpp - gp2) to insure rotational invariance
    int nhsp = -nrp;                                       //number of integrals of type <sp|sp>; note that if nrp != 0 then nrs = 2 due to Aufbau principle
    double eisol = Uss*nrs + Upp*nrp + Udd*nrd + Gss*ngss + Gpp*npp + Gsp*ngsp + Gp2*ngp2 + Hsp*nhsp;
    bool dodorbs = (((atmnr > 20)&&(atmnr < 30))||((atmnr > 38)&&(atmnr < 48))||((atmnr > 71)&&(atmnr < 80)));
    if (dodorbs) {
      //correction for atoms with partially filled d orbitals
      double r016 = SlaterCondonRadialIntegral(atmnr,0,1,6);        //F0sd
      double r066 = SlaterCondonRadialIntegral(atmnr,0,6,6);        //F0dd
      double r244 = SlaterCondonRadialIntegral(atmnr,2,4,4);        //G2sd
      double r266 = SlaterCondonRadialIntegral(atmnr,2,6,6);        //F2dd
      double r466 = SlaterCondonRadialIntegral(atmnr,4,6,6);        //F4dd
      int n016 = Nr016(atmnr);
      int n066 = Nr066(atmnr);
      int n244 = Nr244(atmnr);
      int n266 = Nr266(atmnr);
      int n466 = Nr466(atmnr);
      eisol += n016*r016 + n066*r066 - n244*r244/5.0 - n266*r266/49.0 - n466*r466/49.0;
    }
    return eisol;
  }
  //parameters
  virtual double ZeroOverlap(size_t atm) {return ZeroOverlapMNDOd(atm);}
  virtual double ElementFormTheo(size_t atomicnr) {
    //function returning theoretical heats of formation for elements; values directly in a.u.
    double enth = 0.0;
    switch (atomicnr) {
      case 1:            //H
        enth = -0.437547617217;
        break;
      case 3:            //Li
        enth = -0.188450543318;
        break;
      case 4:            //Be
        enth = -0.889508467058;
        break;
      case 5:            //B
        enth = -2.363567808404;
        break;
      case 6:            //C
        enth = -4.405531461381;
        break;
      case 7:            //N
        enth = -7.418630559411;
        break;
      case 8:            //O
        enth = -11.653156345739;
        break;
      case 9:            //F
        enth = -17.517807628747;
        break;
      case 11:           //Na
        enth = -0.191133244110;
        break;
      case 12:           //Mg
        enth = -0.838577980508;
        break;
      case 13:           //Al
        enth = -2.001448256136;
        break;
      case 14:           //Si
        enth = -2.990566579586;
        break;
      case 15:           //P
        enth = -4.811178862983;
        break;
      case 16:           //S
        enth = -6.805235065709;
        break;
      case 17:           //Cl
        enth = -9.676042537744;
        break;
      case 30:           //Zn
        enth = -1.010065292227;
        break;
      case 32:           //Ge
        enth = -2.787397563744;
        break;
      case 35:           //Br
        enth = -9.100628685109;
        break;
      case 48:           //Cd
        enth = -0.956767450850;
        break;
      case 50:           //Sn
        enth = -3.378148784301;
        break;
      case 53:           //I
        enth = -8.438500792289;
        break;
      case 80:           //Hg
        enth = -1.077330339657;
        break;
      case 82:           //Pb
        enth = -3.874646950710;
        break;
    }
    return enth;
  }
  virtual double betaA0(size_t atomicnr, int L) {
    //function returning the bonding parameters; values stored in eV, returned however in a.u.
    double betaa0 = 0.0;
    switch (atomicnr) {
      case 1:                   //H
        betaa0 = -6.989064;  
        break;
      case 3:                   //Li
        betaa0 = -1.3500400;
        break;
      case 4:                   //Be
        betaa0 = -4.0170960;
        break;
      case 5:                   //B
        betaa0 = -8.252054;
        break;
      case 6:                   //C
        if (L == 0) {betaa0 = -18.985044;}
        else if (L == 1) {betaa0 = -7.934122;}
        break;
      case 7:                   //N
        betaa0 = -20.495758;
        break;
      case 8:                   //O
        betaa0 = -32.688082;
        break;
      case 9:                   //F
        if (L == 0) {betaa0 = -48.290460;}
        else if (L == 1) {betaa0 = -36.508540;}
        break;
      case 11:                  //Na
        if (L == 0) {betaa0 = -1.087382;}
        else if (L == 1) {betaa0 = -0.486239;}
        break;
      case 12:                  //Mg
        if (L == 0) {betaa0 = -1.895884;}
        else if (L == 1) {betaa0 = -2.141089;}
        break;
      case 13:                  //Al
        if (L == 0) {betaa0 = -5.402719;}
        else if (L == 1) {betaa0 = -3.406734;}
        else if (L == 2) {betaa0 = -1.778777;}
        break;
      case 14:                  //Si
        if (L == 0) {betaa0 = -8.210734;}
        else if (L == 1) {betaa0 = -4.884620;}
        else if (L == 2) {betaa0 = -2.608011;}
        break;
      case 15:                  //P
        if (L == 0) {betaa0 = -8.902104;}
        else if (L == 1) {betaa0 = -9.386110;}
        else if (L == 2) {betaa0 = -2.091701;}
        break;
      case 16:                  //S
        if (L == 0) {betaa0 = -10.999545;}
        else if (L == 1) {betaa0 = -12.215437;}
        else if (L == 2) {betaa0 = -1.880669;}
        break;
      case 17:                  //Cl
        if (L == 0) {betaa0 = -6.037292;}
        else if (L == 1) {betaa0 = -19.183386;}
        else if (L == 2) {betaa0 = -1.877782;}
        break;
      case 30:                  //Zn
        if (L == 0) {betaa0 = -5.017261;}
        else if (L == 1) {betaa0 = -0.712060;}
        break;
      case 32:                  //Ge
        if (L == 0) {betaa0 = -4.5164790;}
        else if (L == 1) {betaa0 = -1.7555170;}
        break;
      case 35:                  //Br
        if (L == 0) {betaa0 = -8.314976;}
        else if (L == 1) {betaa0 = -10.507041;}
        else if (L == 2) {betaa0 = -0.962599;}
        break;
      case 48:                  //Cd
        if (L == 0) {betaa0 = -2.771544;}
        else if (L == 1) {betaa0 = -1.805650;}
        break;
      case 50:                  //Sn
        if (L == 0) {betaa0 = -3.2351470;}
        else if (L == 1) {betaa0 = -4.2904160;}
        break;
      case 53:                  //I
        if (L == 0) {betaa0 = -10.699487;}
        else if (L == 1) {betaa0 = -4.941178;}
        else if (L == 2) {betaa0 = -2.350461;}
        break;
      case 80:                  //Hg
        if (L == 0) {betaa0 = -2.218722;}
        else if (L == 1) {betaa0 = -2.909786;}
        break;
      case 82:                  //Pb
        if (L == 0) {betaa0 = -8.0423870;}
        else if (L == 1) {betaa0 = -3.00000;}
        break;
    }
    return betaa0/au2eV;
  }
  virtual double UlX(size_t atomicnr, int L) {
    //function that returns the atomic potential U
    double ulx = 0.0;
    switch (atomicnr) {
      case 1:                 //H
        ulx = -11.906276;
        break;
      case 3:                 //Li
        if (L == 0) {ulx = -5.1280000;}
        else if (L == 1) {ulx = -2.7212000;}
        break;
      case 4:                 //Be
        if (L == 0) {ulx = -16.6023780;}
        else if (L == 1) {ulx = -10.7037710;}
        break;
      case 5:                 //B
        if (L == 0) {ulx = -34.547130;}
        else if (L == 1) {ulx = -23.121690;}
        break;
      case 6:                 //C
        if (L == 0) {ulx = -52.279745;}
        else if (L == 1) {ulx = -39.205558;}
        break;
      case 7:                 //N
        if (L == 0) {ulx = -71.932122;}
        else if (L == 1) {ulx = -57.172319;}
        break;
      case 8:                 //O
        if (L == 0) {ulx = -99.644309;}
        else if (L == 1) {ulx = -77.797472;}
        break;
      case 9:                 //F
        if (L == 0) {ulx = -131.071548;}
        else if (L == 1) {ulx = -105.782137;}
        break;
      case 11:                //Na
        if (L == 0) {ulx = -5.201000;}
        else if (L == 1) {ulx = -2.712573;}
        break;
      case 12:                //Mg
        if (L == 0) {ulx = -15.097000;}
        else if (L == 1) {ulx = -10.650000;}
        break;
      case 13:                //Al
        if (L == 0) {ulx = -28.961830;}
        else if (L == 1) {ulx = -22.814474;}
        else if (L == 2) {ulx = -6.690879;}
        break;
      case 14:                //Si
        if (L == 0) {ulx = -36.051529;}
        else if (L == 1) {ulx = -27.535692;}
        else if (L == 2) {ulx = -14.677439;}
        break;
      case 15:                //P
        if (L == 0) {ulx = -47.055531;}
        else if (L == 1) {ulx = -38.067059;}
        else if (L == 2) {ulx = -23.691597;}
        break;
      case 16:                //S
        if (L == 0) {ulx = -56.889130;}
        else if (L == 1) {ulx = -47.274746;}
        else if (L == 2) {ulx = -25.095118;}
        break;
      case 17:                //Cl
        if (L == 0) {ulx = -69.622971;}
        else if (L == 1) {ulx = -59.100731;}
        else if (L == 2) {ulx = -36.674572;}
        break;
      case 30:                //Zn
        if (L == 0) {ulx = -18.023001;}
        else if (L == 1) {ulx = -12.242166;}
        break;
      case 32:                //Ge
        if (L == 0) {ulx = -33.9493670;}
        else if (L == 1) {ulx = -27.4251050;}
        break;
      case 35:                //Br
        if (L == 0) {ulx = -65.402779;}
        else if (L == 1) {ulx = -54.553753;}
        else if (L == 2) {ulx = -13.728099;}
        break;
      case 48:                //Cd
        if (L == 0) {ulx = -16.969700;}
        else if (L == 1) {ulx = -12.400965;}
        break;
      case 50:                //Sn
        if (L == 0) {ulx = -40.8518020;}
        else if (L == 1) {ulx = -28.5602490;}
        break;
      case 53:                //I
        if (L == 0) {ulx = -62.765354;}
        else if (L == 1) {ulx = -50.292114;}
        else if (L == 2) {ulx = -12.248305;}
        break;
      case 80:                //Hg
        if (L == 0) {ulx = -18.815649;}
        else if (L == 1) {ulx = -13.397114;}
        break;
      case 82:                //Pb
        if (L == 0) {ulx = -47.3196920;}
        else if (L == 1) {ulx = -28.8475600;}
        break;
    }
    return ulx/au2eV;
  }
  virtual double alpha(int atm1, int atm2) {
    //function that returns alpha values for MNDOd
    //alphas stored and returned in 1/Angstrom
    double alpha = 0.0;
    switch (atm1) {
      case 1:     //H
        alpha = 2.544134;
        break;
      case 3:     //Li
        alpha = 1.250140;
        break;
      case 4:     //Be
        alpha = 1.669434;
        break;
      case 5:     //B
        alpha = 2.134993;
        break;
      case 6:     //C
        alpha = 2.546380;
        break;
      case 7:     //N
        alpha = 2.861342;
        break;
      case 8:     //O
        alpha = 3.160604;
        break;
      case 9:     //F
        alpha = 3.419661;
        break;
      case 11:    //Na
        if ((atm2 == 1)||(atm2 == 6)) {alpha = 1.052252;}
        else {alpha = 1.170102;}
        break;
      case 12:    //Mg
        if (atm2 == 1) {alpha = 1.350530;}
        else if ((atm2 == 6)||(atm2 == 16)) {alpha = 1.481721;}
        else {alpha = 1.621470;}
        break;
      case 13:    //Al
        if ((atm2 == 1)||(atm2 == 6)||(atm2 == 13)) {alpha = 1.387880;}
        else {alpha = 1.439879;}
        break;
      case 14:    //Si
        alpha = 1.660069;
        break;
      case 15:    //P
        alpha = 1.852551;
        break;
      case 16:    //S
        alpha = 2.023060;
        break;
      case 17:    //Cl
        alpha = 2.180300;
        break;
      case 30:    //Zn
        alpha = 1.517637;
        break;
      case 32:    //Ge
        alpha = 1.978498;
        break;
      case 35:    //Br
        alpha = 2.091050;
        break;
      case 48:    //Cd
        alpha = 1.424613;
        break;
      case 50:    //Sn
        alpha = 1.800814;
        break;
      case 53:    //I
        alpha = 1.906174;
        break;
      case 80:    //Hg
        alpha = 1.382242;
        break;
      case 82:    //Pb
        alpha = 1.728333;
        break;
    }
    return alpha;
  }
  virtual double Dvalue(size_t atmnr, size_t idx) {
    //function returning the D values needed to calculate eris; values stored in Angstrom, but returned in a.u.
    //note that this function returns both D1 and D2; idx is then either 1 or 2
    double dval = 0.0;
    switch (atmnr) {
      case 1:                    //H
        dval = 0.0;
        break;
      case 3:                    //Li
        if (idx == 1) {dval = 1.0874267;}
        else if (idx == 2) {dval = 0.9227121;}
        else if (idx == 3) {dval = 1.3049119659857598;}
        break;
      case 4:                    //Be
        if (idx == 1) {dval = 0.7605847;}
        else if (idx == 2) {dval = 0.6453775;}
        else if (idx == 3) {dval = 0.9127016133504423;}
        break;
      case 5:                    //B
        if (idx == 1) {dval = 0.5068929;}
        else if (idx == 2) {dval = 0.4301129;}
        else if (idx == 3) {dval = 0.6082714965316228;}
        break;
      case 6:                    //C
        if (idx == 1) {dval = 0.4272845;}
        else if (idx == 2) {dval = 0.3625629;}
        else if (idx == 3) {dval = 0.5127413703933202;}
        break;
      case 7:                    //N
        if (idx == 1) {dval = 0.3386159;}
        else if (idx == 2) {dval = 0.2873251;}
        else if (idx == 3) {dval = 0.4063390532302058;}
        break;
      case 8:                    //O
        if (idx == 1) {dval = 0.2828939;}
        else if (idx == 2) {dval = 0.2400435;}
        else if (idx == 3) {dval = 0.33947277325950603;}
        break;
      case 9:                    //F
        if (idx == 1) {dval = 0.2681377;}
        else if (idx == 2) {dval = 0.2275224;}
        else if (idx == 3) {dval = 0.3217652638236763;}
        break;
      case 11:                   //Na
        if (idx == 1) {dval = 1.722047013848492;}
        else if (idx == 2) {dval = 1.568867318336547;}
        else if (idx == 3) {dval = 2.218713439155453;}
        break;
      case 12:                    //Mg
        if (idx == 1) {dval = 1.067674252498970;}
        else if (idx == 2) {dval = 0.956997330811926;}
        else if (idx == 3) {dval = 1.353398604389077;}
        break;
      case 13:                   //Al
        if (idx == 1) {dval = 0.832308566407991;}               //ss,sp
        else if (idx == 2) {dval = 0.952665777575407;}          //pp
        else if (idx == 3) {dval = 1.347272863055851;}          //D2pp
        else if (idx == 4) {dval = 1.177532394535999;}          //D2sd
        else if (idx == 5) {dval = 0.897301564955694;}          //D1pd
        else if (idx == 6) {dval = 1.154925851527145;}          //D2dd
        break;
      case 14:                   //Si
        if (idx == 1) {dval = 0.716671969127800;}               //ss,sp
        else if (idx == 2) {dval = 0.906995869298913;}          //pp
        else if (idx == 3) {dval = 1.282685859378898;}          //D2pp
        else if (idx == 4) {dval = 1.039958482223044;}          //D2sd
        else if (idx == 5) {dval = 0.864109122473992;}          //D1pd
        else if (idx == 6) {dval = 1.128112373316393;}          //D2dd
        break;
      case 15:                   //P
        if (idx == 1) {dval = 0.679314907108863;}               //ss,sp
        else if (idx == 2) {dval = 0.817686787361311;}          //pp
        else if (idx == 3) {dval = 1.156383744459651;}          //D2pp
        else if (idx == 4) {dval = 0.973232417552050;}          //D2sd
        else if (idx == 5) {dval = 0.789138952424778;}          //D1pd
        else if (idx == 6) {dval = 1.051503868894154;}          //D2dd
        break;
      case 16:                   //S
        if (idx == 1) {dval = 0.633873496222589;}               //ss,sp
        else if (idx == 2) {dval = 0.789776177529618;}          //pp
        else if (idx == 3) {dval = 1.116912181501567;}          //D2pp
        else if (idx == 4) {dval = 0.915348224929503;}          //D2sd
        else if (idx == 5) {dval = 0.759895375983441;}          //D1pd
        else if (idx == 6) {dval = 1.007144157617617;}          //D2dd
        break;
      case 17:                   //Cl
        if (idx == 1) {dval = 0.578881496107948;}               //ss,sp
        else if (idx == 2) {dval = 0.749747116762158;}          //pp
        else if (idx == 3) {dval = 1.060302540875168;}          //D2pp
        else if (idx == 4) {dval = 0.845123889350194;}          //D2sd
        else if (idx == 5) {dval = 0.712779146193005;}          //D1pd
        else if (idx == 6) {dval = 0.927879686590404;}          //D2dd
        break;
      case 30:                   //Zn
        if (idx == 1) {dval = 0.763507018483803;}
        else if (idx == 2) {dval = 1.300999463764879;}
        else if (idx == 3) {dval = 1.839891086296417;}
        break;
      case 32:                   //Ge
        if (idx == 1) {dval = 0.6644269;}
        else if (idx == 2) {dval = 0.5555542;}
        else if (idx == 3) {dval = 0.785672284273335;}
        break;
      case 35:                   //Br
        if (idx == 1) {dval = 0.600967020923898;}               //ss,sp
        else if (idx == 2) {dval = 0.783399557210443;}          //pp
        else if (idx == 3) {dval = 1.107894278564086;}          //D2pp
        else if (idx == 4) {dval = 0.856496644208800;}          //D2sd
        else if (idx == 5) {dval = 0.778100505006181;}          //D1pd
        else if (idx == 6) {dval = 1.079778579131701;}          //D2dd
        break;
      case 48:                   //Cd
        if (idx == 1) {dval = 0.992909075368672;}
        else if (idx == 2) {dval = 0.891159574668322;}
        else if (idx == 3) {dval = 1.260289956734580;}
        break;
      case 50:                   //Sn
        if (idx == 1) {dval = 0.8306740;}
        else if (idx == 2) {dval = 0.7017967;}
        else if (idx == 3) {dval = 0.9924904111686824;}
        break;
      case 53:                   //I
        if (idx == 1) {dval = 0.540129610588824;}               //ss,sp
        else if (idx == 2) {dval = 0.864645250800544;}          //pp
        else if (idx == 3) {dval = 1.222793040323616;}          //D2pp
        else if (idx == 4) {dval = 0.759451393810603;}          //D2sd
        else if (idx == 5) {dval = 0.859183859369850;}          //D1pd
        else if (idx == 6) {dval = 1.290775710431637;}          //D2dd
        break;
      case 80:                   //Hg
        if (idx == 1) {dval = 0.923531152348794;}
        else if (idx == 2) {dval = 0.936259354817651;}
        else if (idx == 3) {dval = 1.324070677481806;}
        break;
      case 82:                   //Pb
        if (idx == 1) {dval = 0.8216177;}
        else if (idx == 2) {dval = 0.7666867;}
        else if (idx == 3) {dval = 1.0842587292310726;}
        break;
    }
    return dval*dist_Angstrom2aum1;
  }
  virtual double rhocore(size_t atmnr) {
    //function returning rho values used to calculate VAB integrals; values already in atomic units
    double rho_ = 0.0;
    if (atmnr == 11) {rho_ = 1.530553;}                //Na
    else if (atmnr == 12) {rho_ = 1.350776;}           //Mg
    else if (atmnr == 29) {rho_ = 1.589234;}           //Zn
    else if (atmnr == 48) {rho_ = 1.721186;}           //Cd
    else if (atmnr == 80) {rho_ = 1.636072;}           //Hg
    else {rho_ = rho(atmnr,0);}                        //all other atoms
    return rho_;
  }
  virtual double rho(size_t atmnr, size_t l) {
    //function returning the rho values needed to calculate eris; values stored in Angstrom but returned in atomic units
    double rho = 0.0;
    switch (atmnr) {
      case 1:                    //H
        if (l == 0) {rho = 0.560345403;}
        else if (l == 1) {rho = 0.0;}
        else if (l == 2) {rho = 0.0;}
        break;
      case 3:                    //Li
        if (l == 0) {rho = 0.986207884;}
        else if (l == 1) {rho = 1.165672376;}
        else if (l == 2) {rho = 1.011953732;}
        break;
      case 4:                    //Be
        if (l == 0) {rho = 0.799924235;}
        else if (l == 1) {rho = 0.788356095;}
        else if (l == 2) {rho = 0.687877905;}
        break;
      case 5:                    //B
        if (l == 0) {rho = 0.679822279;}
        else if (l == 1) {rho = 0.539445596;}
        else if (l == 2) {rho = 0.476128306;}
        break;
      case 6:                    //C
        if (l == 0) {rho = 0.588660438;}
        else if (l == 1) {rho = 0.430253872;}
        else if (l == 2) {rho = 0.395733736;}
        break;
      case 7:                    //N
        if (l == 0) {rho = 0.529751104;}
        else if (l == 1) {rho = 0.337322211;}
        else if (l == 2) {rho = 0.325583327;}
        break;
      case 8:                    //O
        if (l == 0) {rho = 0.466881794;}
        else if (l == 1) {rho = 0.275821517;}
        else if (l == 2) {rho = 0.2786282;}
        break;
      case 9:                    //F
        if (l == 0) {rho = 0.425491557;}
        else if (l == 1) {rho = 0.243848996;}
        else if (l == 2) {rho = 0.255793341;}
        break;
      case 11:                   //Na
        if (l == 0) {rho = 1.567071093635136;}
        else if (l == 1) {rho = 1.830222533372990;}
        else if (l == 2) {rho = 1.417020172941393;}
        break;
      case 12:                   //Mg
        if (l == 0) {rho = 0.976229438953369;}
        else if (l == 1) {rho = 1.222534313856169;}
        else if (l == 2) {rho = 0.864372080385308;}
        break;
      case 13:                   //Al
        if (l == 0) {rho = 0.760157668073853;}
        else if (l == 1) {rho = 1.029192272831372;}
        else if (l == 2) {rho = 0.860459625876539;}
        else if (l == 3) {rho = 0.306374315792088;}
        else if (l == 4) {rho = 0.885024468142918;}
        else if (l == 5) {rho = 1.118515830987198;}
        else if (l == 6) {rho = 0.972795025314052;}
        break;
      case 14:                   //Si
        if (l == 0) {rho = 0.670271695839240;}
        else if (l == 1) {rho = 0.897134091517244;}
        else if (l == 2) {rho = 0.819209949824937;}
        else if (l == 3) {rho = 2.645886041500000;}
        else if (l == 4) {rho = 0.853697951546307;}
        else if (l == 5) {rho = 1.092547670500620;}
        else if (l == 6) {rho = 0.950210003881803;}
        break;
      case 15:                   //P
        if (l == 0) {rho = 0.627175684946922;}
        else if (l == 1) {rho = 0.589029003839085;}
        else if (l == 2) {rho = 0.738545109729496;}
        else if (l == 3) {rho = 0.230274672062604;}
        else if (l == 4) {rho = 0.783143597305115;}
        else if (l == 5) {rho = 1.018354314389296;}
        else if (l == 6) {rho = 0.885682597109276;}
        break;
      case 16:                   //S
        if (l == 0) {rho = 0.590328325842112;}
        else if (l == 1) {rho = 0.508801745643787;}
        else if (l == 2) {rho = 0.713335789594268;}
        else if (l == 3) {rho = 2.645886041500000;}
        else if (l == 4) {rho = 0.753084041637511;}
        else if (l == 5) {rho = 0.975393080769647;}
        else if (l == 6) {rho = 0.848318378715980;}
        break;
      case 17:                   //Cl
        if (l == 0) {rho = 0.544980795824593;}
        else if (l == 1) {rho = 0.450625897868114;}
        else if (l == 2) {rho = 0.677180931836914;}
        else if (l == 3) {rho = 1.647969914212899;}
        else if (l == 4) {rho = 0.703877431305918;}
        else if (l == 5) {rho = 0.898627489661325;}
        else if (l == 6) {rho = 0.781553842697954;}
        break;
      case 30:                   //Zn
        if (l == 0) {rho = 0.841029168406762;}
        else if (l == 1) {rho = 1.140966140980049;}
        else if (l == 2) {rho = 1.155784659246568;}
        break;
      case 32:                   //Ge
        if (l == 0) {rho = 0.734624198;}
        else if (l == 1) {rho = 0.726135254;}
        else if (l == 2) {rho = 0.608610513;}
        break;
      case 35:                   //Br
        if (l == 0) {rho = 0.589069893791056;}
        else if (l == 1) {rho = 0.485997257226682;}
        else if (l == 2) {rho = 0.695958153863445;}
        else if (l == 3) {rho = 0.364713943660608;}
        else if (l == 4) {rho = 0.778171834446568;}
        else if (l == 5) {rho = 1.059880891099237;}
        else if (l == 6) {rho = 0.896749329209849;}
        break;
      case 50:                   //Sn
        if (l == 0) {rho = 0.734624198;}
        else if (l == 1) {rho = 0.821688399;}
        else if (l == 2) {rho = 0.712428177;}
        break;
      case 53:                   //I
        if (l == 0) {rho = 0.600947580115929;}
        else if (l == 1) {rho = 0.525723746001474;}
        else if (l == 2) {rho = 0.757767568929242;}
        else if (l == 3) {rho = 0.580249398341799;}
        else if (l == 4) {rho = 0.901076613080184;}
        else if (l == 5) {rho = 1.275688097573650;}
        else if (l == 6) {rho = 1.059659488896923;}
        break;
      case 80:                   //Hg
        if (l == 0) {rho = 0.865815999544531;}
        else if (l == 1) {rho = 1.059780622389568;}
        else if (l == 2) {rho = 0.811246212763297;}
        break;
      case 82:                   //Pb
        if (l == 0) {rho = 0.734624198;}
        else if (l == 1) {rho = 0.816789939;}
        else if (l == 2) {rho = 0.755508834;}
        break;
    }
    return rho*dist_Angstrom2aum1;
  }
  virtual double PrincipalQuantNumb(size_t atmnr, size_t l) {
    //function returning principal quantum numbers for atoms for the different orbitals
    double pqn = 0.0;
    if ((atmnr < 3)&&(l == 0)) {pqn = 1.0;}
    else if ((atmnr == 2)&&(l == 1)) {pqn = 2.0;}
    else if ((atmnr > 2)&&(atmnr < 10)) {pqn = 2.0;}
    else if ((atmnr == 10)&&(l == 0)) {pqn = 3.0;}
    else if ((atmnr == 10)&&(l == 1)) {pqn = 2.0;}
    else if ((atmnr > 10)&&(atmnr < 18)) {pqn = 3.0;}
    else if ((atmnr == 18)&&(l == 0)) {pqn = 4.0;}
    else if ((atmnr == 18)&&(l == 1)) {pqn = 3.0;}
    else if ((atmnr > 18)&&(atmnr < 36)&&(l < 2)) {pqn = 4.0;}
    else if ((atmnr > 18)&&(atmnr < 30)&&(l == 2)) {pqn = 3.0;}
    else if ((atmnr > 30)&&(atmnr < 36)&&(l == 2)) {pqn = 4.0;}
    else if ((atmnr == 36)&&(l == 0)) {pqn = 5.0;}
    else if ((atmnr == 36)&&(l == 1)) {pqn = 4.0;}
    else if ((atmnr > 36)&&(atmnr < 54)&&(l < 2)) {pqn = 5.0;}
    else if ((atmnr > 36)&&(atmnr < 48)&&(l == 2)) {pqn = 4.0;}
    else if ((atmnr > 48)&&(atmnr < 54)&&(l == 2)) {pqn = 5.0;}
    else if ((atmnr == 54)&&(l == 0)) {pqn = 6.0;}
    else if ((atmnr == 54)&&(l == 1)) {pqn = 5.0;}
    else if ((atmnr > 54)&&(atmnr < 84)&&(l < 2)) {pqn = 6.0;}
    else if ((atmnr > 54)&&(atmnr < 84)&&(l == 2)) {pqn = 5.0;}
    return pqn;
  }
  virtual double zn(size_t atmnr, size_t l) {
    //function returning the zn values needed to calculate 1-center eris; values stored and returned directly in atomic units
    double zn_ = 0.0;
    if (atmnr == 11) {                                 //Na
      if (l == 0) {zn_ = 0.654113;}                        //zs
      else if (l == 1) {zn_ = 0.564409;}                   //zp
    }
    else if (atmnr == 12) {                            //Mg
      if (l == 0) {zn_ = 1.050000;}                        //zs
      else if (l == 1) {zn_ = 0.925272;}                   //zp
    }
    else if (atmnr == 13) {                            //Al
      if (l == 0) {zn_ = 1.348458;}                        //zs
      else if (l == 1) {zn_ = 0.929479;}                   //zp
      else if (l == 2) {zn_ = 0.916383;}                   //zd
    }
    else if (atmnr == 14) {                            //Si
      if (l == 0) {zn_ = 1.529292;}                        //zs
      else if (l == 1) {zn_ = 0.976281;}                   //zp
      else if (l == 2) {zn_ = 0.938164;}                   //zd
    }
    else if (atmnr == 15) {                            //P
      if (l == 0) {zn_ = 1.634376;}                        //zs
      else if (l == 1) {zn_ = 1.082912;}                   //zp
      else if (l == 2) {zn_ = 1.006515;}                   //zd
    }
    else if (atmnr == 16) {                            //S
      if (l == 0) {zn_ = 1.736391;}                        //zs
      else if (l == 1) {zn_ = 1.121182;}                   //zp
      else if (l == 2) {zn_ = 1.050847;}                   //zd
    }
    else if (atmnr == 17) {                            //Cl
      if (l == 0) {zn_ = 1.880875;}                        //zs
      else if (l == 1) {zn_ = 1.181042;}                   //zp
      else if (l == 2) {zn_ = 1.140616;}                   //zd
    }
    else if (atmnr == 30) {                            //Zn
      if (l == 0) {zn_ = 1.566000;}                        //zs
      else if (l == 1) {zn_ = 0.862840;}                   //zp
    }
    else if (atmnr == 35) {                            //Br
      if (l == 0) {zn_ = 2.235816;}                        //zs
      else if (l == 1) {zn_ = 1.432927;}                   //zp
      else if (l == 2) {zn_ = 1.242578;}                   //zd
    }
    else if (atmnr == 48) {                            //Cd
      if (l == 0) {zn_ = 1.763148;}                        //zs
      else if (l == 1) {zn_ = 1.525519;}                   //zp
    }
    else if (atmnr == 53) {                            //I
      if (l == 0) {zn_ = 2.672411;}                        //zs
      else if (l == 1) {zn_ = 1.572299;}                   //zp
      else if (l == 2) {zn_ = 1.258848;}                   //zd
    }
    else if (atmnr == 80) {                            //Hg
      if (l == 0) {zn_ = 2.186000;}                        //zs
      else if (l == 1) {zn_ = 1.705005;}                   //zp
    }
    return zn_;
  }
  virtual double F0sd(size_t atmnr) {
    //function returning the F0sd 1-center eris; values stored in eV, returned in a.u.
    return 0.0;
  }
  virtual double G2sd(size_t atmnr) {
    //function returning the G2sd 1-center eris; values stored in eV, returned in a.u.
    return 0.0;
  }
  virtual double eri1Center(int atmnr, int Lbra, int Lket) {
    //function that gives back the semi-empirical 1-center eris
    //values stored in eV but returned in a.u.
    //Lbra is the sum of azimuthal quantum numbers for bra (ss = 0; pp = 2; sp = 1)
    //Lket is the sum of azimuthal quantum numbers for ket (ss = 0; pp = 2; sp = 1; p*p* = -2)
    double eri = 0.0;
    switch (atmnr) {
      case 1:            //H
        eri = 12.848;                                                                            //(ss|ss)
        break;
      case 3:           //Li
        if ((Lbra == 0)&&(Lket == 0)) {eri = 7.3;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.0;}                                          //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.42;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.52;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.83;}                                         //(sp|sp)||(ps|ps)
        break;
      case 4:          //Be
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.0;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.97;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.43;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.22;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.28;}                                         //(sp|sp)||(ps|ps)
        break;
      case 5:          //B
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.59;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.86;}                                         //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 9.56;}           //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.86;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.81;}                                         //(sp|sp)||(ps|ps)
        break;
      case 6:          //C
        if ((Lbra == 0)&&(Lket == 0)) {eri = 12.23;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.08;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 11.47;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.84;}         //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.43;}                                         //(sp|sp)||(ps|ps)
        break;
      case 7:          //N
        if ((Lbra == 0)&&(Lket == 0)) {eri = 13.59;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 12.98;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 12.66;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 11.59;}        //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.14;}                                         //(sp|sp)||(ps|ps)
        break;
      case 8:          //O
        if ((Lbra == 0)&&(Lket == 0)) {eri = 15.42;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 14.52;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 14.48;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 12.98;}        //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.94;}                                         //(sp|sp)||(ps|ps)
        break;
      case 9:          //F
        if ((Lbra == 0)&&(Lket == 0)) {eri = 16.92;}                                             //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 16.71;}                                        //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 17.25;}          //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 14.91;}        //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 4.83;}                                         //(sp|sp)||(ps|ps)
        break;
      case 11:         //Na
        if ((Lbra == 0)&&(Lket == 0)) {eri = 4.594445;}                                          //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 4.299198;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 4.147574;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 3.796957;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.534409;}                                     //(sp|sp)||(ps|ps)
        break;
      case 12:         //Mg
        if ((Lbra == 0)&&(Lket == 0)) {eri = 7.375133;}                                          //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.047954;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.888907;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.224599;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.726734;}                                     //(sp|sp)||(ps|ps)
        break;
      case 13:         //Al
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.471485;}                                          //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.080002;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.791266;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.252903;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.778059;}                                     //(sp|sp)||(ps|ps)
        break;
      case 14:         //Si
        if ((Lbra == 0)&&(Lket == 0)) {eri = 10.741647;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.436497;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.560667;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.567751;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.877539;}                                     //(sp|sp)||(ps|ps)
        break;
      case 15:         //P
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.4797530;}                                        //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.248723;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.557570;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.285092;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.107804;}                                     //(sp|sp)||(ps|ps)
        break;
      case 16:         //S
        if ((Lbra == 0)&&(Lket == 0)) {eri = 12.196301;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.540233;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.853901;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.542547;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.646352;}                                     //(sp|sp)||(ps|ps)
        break;
      case 17:         //Cl
        if ((Lbra == 0)&&(Lket == 0)) {eri = 13.211148;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.996201;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 9.419496;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.945248;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.081499;}                                     //(sp|sp)||(ps|ps)
        break;
      case 30:         //Zn
        if ((Lbra == 0)&&(Lket == 0)) {eri = 8.560728;}                                          //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.139648;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.490036;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.505403;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.532946;}                                     //(sp|sp)||(ps|ps)
        break;
      case 32:         //Ge
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.8;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.3;}                                          //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.3;}            //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.5;}          //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.3;}                                          //(sp|sp)||(ps|ps)
        break;
      case 35:         //Br
        if ((Lbra == 0)&&(Lket == 0)) {eri = 12.222356;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.535464;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.263721;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.482167;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.749522;}                                     //(sp|sp)||(ps|ps)
        break;
      case 48:         //Cd
        if ((Lbra == 0)&&(Lket == 0)) {eri = 7.904434;}                                          //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.480000;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.515707;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.518664;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.636744;}                                     //(sp|sp)||(ps|ps)
        break;
      case 50:         //Sn
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.8;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.3;}                                          //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.3;}            //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.5;}          //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.3;}                                          //(sp|sp)||(ps|ps)
        break;
      case 53:         //I
        if ((Lbra == 0)&&(Lket == 0)) {eri = 11.980782;}                                         //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.709372;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.855902;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.718557;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.071475;}                                     //(sp|sp)||(ps|ps)
        break;
      case 80:         //Hg
        if ((Lbra == 0)&&(Lket == 0)) {eri = 8.315649;}                                          //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.115259;}                                     //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.212173;}       //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.171250;}     //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.835941;}                                     //(sp|sp)||(ps|ps)
        break;
      case 82:         //Pb
        if ((Lbra == 0)&&(Lket == 0)) {eri = 9.8;}                                               //(ss|ss)
        else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.3;}                                          //(pp|pp)
        else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.3;}            //(ss|pp)||(pp|ss)
        else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.5;}          //(pp|p*p*)||(p*p*|pp)
        else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.3;}                                          //(sp|sp)||(ps|ps)
        break;
    }
    return eri/au2eV;
  }
};
class PM6: public MNDOd {
  //implementation of Stewart's PM6
  //J. J. P Stewart, J. Mol. Model., 13, 1173, 2007
  //bX quantities are "barred" tensors, which are used only for the open-shell case
public:
  PM6(BSet _bset, Molecule _mol, std::string _openclosed = "0", std::string _corecorrection = "0"): MNDOd(_bset,_mol,_openclosed,_corecorrection) {}
  ~PM6() {}
  //getters
  std::string Type() {return "PM6";}
  double gfactor(size_t iatm1, size_t iatm2, double RAB) {
    //function determining special factors in the nuclear repulsion terms of MNDO, zero for PM6
    return 0.0;
  }
  double gfactor_dR(size_t iatm1, size_t iatm2) {
    //function determining the derivative of the special factors in the nuclear repulsion terms of MNDO, zero for PM6
    return 0.0;
  }
  double AM1factor(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {
    //RAB in a.u.
    double RR = RAB*dist_Angstrom2au;
    double AM1factor = K1[0]*exp(-L1[0]*(RR - M1[0])*(RR - M1[0]))/RAB;
    return AM1factor;
  }
  double AM1factor_dR(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {
    //RAB in a.u.
    double RR = RAB*dist_Angstrom2au;
    double am1factor_dr = exp(-L1[0]*(RR - M1[0])*(RR - M1[0]));
    am1factor_dr *= K1[0]*(2.0*L1[0]*(RR - M1[0])*RR + 1.0)/(RAB*RAB);
    return am1factor_dr;  
  }
  double AM1factor_dR2(double RAB, std::vector<double> & K1, std::vector<double> & L1, std::vector<double> & M1) {
    //RAB in a.u.
    double RR = RAB*dist_Angstrom2au;
    double am1factor_dr2 = exp(-L1[0]*(RR - M1[0])*(RR - M1[0]));
    am1factor_dr2 *= K1[0]*(L1[0]*M1[0]*RR - 1.0 - RR*L1[0]*(RR - M1[0])*(2.0*L1[0]*(RR - M1[0])*RR + 1.0))/(RAB*RAB*RAB);
    return 2.0*am1factor_dr2;  
  }
  //other auxiliary functions
  void checkAtoms() {
    //function that checks whether the atoms in the molecule are defined in PM6
    bool def;
    for (size_t iatm = 0; iatm < Natoms; ++iatm) {
      def = false;
      if ((atoms[iatm] > 0)&&(atoms[iatm] < 58)) {def = true;}                                  //H - La
      else if ((atoms[iatm] > 70)&&(atoms[iatm] < 84)) {def = true;}                            //Lu - Bi
      if (!def) {throw("ERROR: MNDOd.hpp: PM6: checkAtoms(): atom not specified in PM6");}
    }
  }
  double ZeroOverlap(size_t atm) {return ZeroOverlapPM6(atm);}
  double ElementFormTheo(size_t atomicnr) {
    //function returning theoretical heats of formation for elements; values directly in a.u.
    double enth = 0.0;
    switch (atomicnr) {
      case 1:      //H
        enth = -0.413318125151;
        break;
      case 2:      //He
        enth = -1.988015017788;
        break;
      case 3:      //Li
        enth = -0.173086091143;
        break;
      case 4:      //Be
        enth = -0.924900640372;
        break;
      case 5:      //B
        enth = -1.820439745312;
        break;
      case 6:      //C
        enth = -4.209838250663;
        break;
      case 7:      //N
        enth = -6.392702494150;
        break;
      case 8:      //O
        enth = -10.491362862110;
        break;
      case 9:      //F
        enth = -17.204475636336;
        break;
      case 10:     //Ne
        enth = -10.396582859376;
        break;
      case 11:     //Na
        enth = -0.166737314346;
        break;
      case 12:     //Mg
        enth = -0.809702454912;
        break;
      case 13:     //Al
        enth = -1.766247021100;
        break;
      case 14:     //Si
        enth = -2.507060328871;
        break;
      case 15:     //P
        enth = -5.115641699578;
        break;
      case 16:     //S
        enth = -6.295571931031;
        break;
      case 17:     //Cl
        enth = -9.294249086848;
        break;
      case 18:     //Ar
        enth = -10.718788748266;
        break;
      case 19:     //K
        enth = -0.139688156749;
        break;
      case 20:     //Ca
        enth = -0.581166445118;
        break;
      case 21:     //Sc
        enth = -1.249460624932;
        break;
      case 22:     //Ti
        enth = -2.332123707587;
        break;
      case 23:     //V
        enth = -3.697495328480;
        break;
      case 24:     //Cr
        enth = -6.825262026078;
        break;
      case 25:     //Mn
        enth = -7.195575948685;
        break;
      case 26:     //Fe
        enth = -15.685908308127;
        break;
      case 27:     //Co
        enth = -6.161267318365;
        break;
      case 28:     //Ni
        enth = -17.829507255915;
        break;
      case 29:     //Cu
        enth = -24.129433991485;
        break;
      case 30:     //Zn
        enth = -1.005987072532;
        break;
      case 31:     //Ga
        enth = -2.101569663145;
        break;
      case 32:     //Ge
        enth = -2.803685893973;
        break;
      case 33:     //As
        enth = -4.506549809187;
        break;
      case 34:     //Se
        enth = -5.094218165542;
        break;
      case 35:     //Br
        enth = -8.340660342442;
        break;
      case 36:     //Kr
        enth = -9.584257568324;
        break;
      case 37:     //Rb
        enth = -0.133639107455;
        break;
      case 38:     //Sr
        enth = -0.597238211569;
        break;
      case 39:     //Y
        enth = -0.806853345822;
        break;
      case 40:     //Zr
        enth = -1.931708769398;
        break;
      case 41:     //Nb
        enth = -3.869458292116;
        break;
      case 42:     //Mo
        enth = -6.914097351219;
        break;
      case 43:     //Tc
        enth = -7.079283134137;
        break;
      case 44:     //Ru
        enth = -6.990641236017;
        break;
      case 45:     //Rh
        enth = -7.328837778588;
        break;
      case 46:     //Pd
        enth = -17.049324692693;
        break;
      case 47:     //Ag
        enth = -9.092536335165;
        break;
      case 48:     //Cd
        enth = -0.831060281153;
        break;
      case 49:     //In
        enth = -1.933251760291;
        break;
      case 50:     //Sn
        enth = -2.608267090471;
        break;
      case 51:     //Sb
        enth = -4.533211224646;
        break;
      case 52:     //Te
        enth = -11.031220676187;
        break;
      case 53:     //I
        enth = -9.118563752632;
        break;
      case 54:     //Xe
        enth = -35.593307372443;
        break;
      case 55:     //Cs
        enth = -0.137758853888;
        break;
      case 56:     //Ba
        enth = -0.551723032188;
        break;
      case 57:     //La
        enth = -2.049453018602;
        break;
      case 71:     //Lu
        enth = -1.424045493233;
        break;
      case 72:     //Hf
        enth = -2.242740836482;
        break;
      case 73:     //Ta
        enth = -4.506186111363;
        break;
      case 74:     //W
        enth = -6.214070631262;
        break;
      case 75:     //Re
        enth = -6.721545839758;
        break;
      case 76:     //Os
        enth = -8.997800214703;
        break;
      case 77:     //Ir
        enth = -7.023877074990;
        break;
      case 78:     //Pt
        enth = -16.005543476446;
        break;
      case 79:     //Au
        enth = -20.029297219625;
        break;
      case 80:     //Hg
        enth = -1.060021149904;
        break;
      case 81:     //Tl
        enth = -2.068204031409;
        break;
      case 82:     //Pb
        enth = -3.081045808707;
        break;
      case 83:     //Bi
        enth = -4.622465267990;
        break;
    }
    return enth;
  }
  double betaA0(size_t atomicnr, int L) {
    //function returning the bonding parameters; values stored in eV, returned however in a.u.
    double betaa0 = 0.0;
    switch (atomicnr) {
      case 1:                //H
        betaa0 = -8.352984;
        break;
      case 2:                 //He
        if (L == 0) {betaa0 = -58.903774;}
        else if (L == 1) {betaa0 = -37.039974;}
        break;
      case 3:                 //Li
        if (L == 0) {betaa0 = -2.283946;}
        else if (L == 1) {betaa0 = -7.535573;}
        break;
      case 4:                 //Be
        if (L == 0) {betaa0 = -3.199549;}
        else if (L == 1) {betaa0 = -4.451920;}
        break;
      case 5:                 //B
        if (L == 0) {betaa0 = -4.959706;}
        else if (L == 1) {betaa0 = -4.656753;}
        break;
      case 6:                 //C
        if (L == 0) {betaa0 = -15.385236;}
        else if (L == 1) {betaa0 = -7.471929;}
        break;
      case 7:                 //N
        if (L == 0) {betaa0 = -17.979377;}
        else if (L == 1) {betaa0 = -15.055017;}
        break;
      case 8:                 //O
        if (L == 0) {betaa0 = -65.635137;}
        else if (L == 1) {betaa0 = -21.622604;}
        break;
      case 9:                 //F
        if (L == 0) {betaa0 = -69.922593;}
        else if (L == 1) {betaa0 = -30.448165;}
        break;
      case 10:                //Ne
        if (L == 0) {betaa0 = -69.793475;}
        else if (L == 1) {betaa0 = -33.261962;}
        break;
      case 11:                //Na
        if (L == 0) {betaa0 = 0.244853;}
        else if (L == 1) {betaa0 = 0.491998;}
        break;
      case 12:                //Mg
        if (L == 0) {betaa0 = -9.604932;}
        else if (L == 1) {betaa0 = 3.416908;}
        break;
      case 13:                //Al
        if (L == 0) {betaa0 = -18.375229;}
        else if (L == 1) {betaa0 = -9.382700;}
        else if (L == 2) {betaa0 = -20.840474;}
        break;
      case 14:                //Si
        if (L == 0) {betaa0 = -8.686909;}
        else if (L == 1) {betaa0 = -1.856482;}
        else if (L == 2) {betaa0 = -6.360627;}
        break;
      case 15:                //P
        if (L == 0) {betaa0 = -14.583780;}
        else if (L == 1) {betaa0 = -11.744725;}
        else if (L == 2) {betaa0 = -20.099893;}
        break;
      case 16:                //S
        if (L == 0) {betaa0 = -13.827440;}
        else if (L == 1) {betaa0 = -7.664613;}
        else if (L == 2) {betaa0 = -9.986172;}
        break;
      case 17:                //Cl
        if (L == 0) {betaa0 = -2.367988;}
        else if (L == 1) {betaa0 = -13.802139;}
        else if (L == 2) {betaa0 = -4.037751;}
        break;
      case 18:                //Ar
        if (L == 0) {betaa0 = -8.839842;}
        else if (L == 1) {betaa0 = -28.427303;}
        break;
      case 19:                //K
        if (L == 0) {betaa0 = -8.755195;}
        else if (L == 1) {betaa0 = -1.788061;}
        break;
      case 20:                //Ca
        if (L == 0) {betaa0 = -4.343881;}
        else if (L == 1) {betaa0 = -1.296612;}
        break;
      case 21:                //Sc
        if (L == 0) {betaa0 = -8.620944;}
        else if (L == 1) {betaa0 = 3.075948;}
        else if (L == 2) {betaa0 = -9.768661;}
        break;
      case 22:                //Ti
        if (L == 0) {betaa0 = 3.389142;}
        else if (L == 1) {betaa0 = -3.355350;}
        else if (L == 2) {betaa0 = -1.842829;}
        break;
      case 23:                //V
        if (L == 0) {betaa0 = -1.211330;}
        else if (L == 1) {betaa0 = 0.740746;}
        else if (L == 2) {betaa0 = 3.153669;}
        break;
      case 24:                //Cr
        if (L == 0) {betaa0 = -5.122615;}
        else if (L == 1) {betaa0 = 3.926711;}
        else if (L == 2) {betaa0 = -4.230550;}
        break;
      case 25:                //Mn
        if (L == 0) {betaa0 = -4.185290;}
        else if (L == 1) {betaa0 = -3.479630;}
        else if (L == 2) {betaa0 = -13.473190;}
        break;
      case 26:                //Fe
        if (L == 0) {betaa0 = 8.027621;}
        else if (L == 1) {betaa0 = -1.125760;}
        else if (L == 2) {betaa0 = -3.507531;}
        break;
      case 27:                //Co
        if (L == 0) {betaa0 = -8.992062;}
        else if (L == 1) {betaa0 = -0.100000;}
        else if (L == 2) {betaa0 = -2.481509;}
        break;
      case 28:                //Ni
        if (L == 0) {betaa0 = -9.151521;}
        else if (L == 1) {betaa0 = -8.086696;}
        else if (L == 2) {betaa0 = -8.655910;}
        break;
      case 29:                //Cu
        if (L == 0) {betaa0 = -9.369508;}
        else if (L == 1) {betaa0 = -0.100000;}
        else if (L == 2) {betaa0 = -16.982092;}
        break;
      case 30:                //Zn
        if (L == 0) {betaa0 = -13.276583;}
        else if (L == 1) {betaa0 = 1.479642;}
        break;
      case 31:                //Ga
        if (L == 0) {betaa0 = -10.808320;}
        else if (L == 1) {betaa0 = -4.185500;}
        break;
      case 32:                //Ge
        if (L == 0) {betaa0 = -14.854297;}
        else if (L == 1) {betaa0 = -2.591260;}
        break;
      case 33:                //As
        if (L == 0) {betaa0 = -11.963725;}
        else if (L == 1) {betaa0 = -7.340073;}
        else if (L == 2) {betaa0 = 3.753005;}
        break;
      case 34:                //Se
        if (L == 0) {betaa0 = 2.636001;}
        else if (L == 1) {betaa0 = -9.557700;}
        break;
      case 35:                //Br
        if (L == 0) {betaa0 = -32.131665;}
        else if (L == 1) {betaa0 = -9.514484;}
        else if (L == 2) {betaa0 = -9.839124;}
        break;
      case 36:                //Kr
        if (L == 0) {betaa0 = -2.727088;}
        else if (L == 1) {betaa0 = -16.142951;}
        break;
      case 37:                //Rb
        if (L == 0) {betaa0 = 9.998744;}
        else if (L == 1) {betaa0 = 1.343004;}
        break;
      case 38:                //Sr
        if (L == 0) {betaa0 = -6.253108;}
        else if (L == 1) {betaa0 = -9.844498;}
        break;
      case 39:                //Y
        if (L == 0) {betaa0 = 0.343336;}
        else if (L == 1) {betaa0 = -3.180807;}
        else if (L == 2) {betaa0 = -4.508957;}
        break;
      case 40:                //Zr
        if (L == 0) {betaa0 = 9.551952;}
        else if (L == 1) {betaa0 = -4.551915;}
        else if (L == 2) {betaa0 = -3.213274;}
        break;
      case 41:                //Nb
        if (L == 0) {betaa0 = -12.045244;}
        else if (L == 1) {betaa0 = 1.465762;}
        else if (L == 2) {betaa0 = -5.920160;}
        break;
      case 42:                //Mo
        if (L == 0) {betaa0 = -0.189344;}
        else if (L == 1) {betaa0 = 7.017762;}
        else if (L == 2) {betaa0 = -10.941126;}
        break;
      case 43:                //Tc
        if (L == 0) {betaa0 = -2.791024;}
        else if (L == 1) {betaa0 = -8.086697;}
        else if (L == 2) {betaa0 = -5.724335;}
        break;
      case 44:                //Ru
        if (L == 0) {betaa0 = -12.859508;}
        else if (L == 1) {betaa0 = -8.475518;}
        else if (L == 2) {betaa0 = -3.830797;}
        break;
      case 45:                //Rh
        if (L == 0) {betaa0 = -8.222141;}
        else if (L == 1) {betaa0 = -15.556691;}
        else if (L == 2) {betaa0 = -13.396182;}
        break;
      case 46:                //Pd
        if (L == 0) {betaa0 = -8.038245;}
        else if (L == 1) {betaa0 = 0.740037;}
        else if (L == 2) {betaa0 = -2.394498;}
        break;
      case 47:                //Ag
        if (L == 0) {betaa0 = -6.129623;}
        else if (L == 2) {betaa0 = 1.004115;}
        else if (L == 1) {betaa0 = -69.238347;}
        break;
      case 48:                //Cd
        if (L == 0) {betaa0 = -11.613183;}
        else if (L == 1) {betaa0 = 1.663178;}
        break;
      case 49:                //In
        if (L == 0) {betaa0 = -1.982376;}
        else if (L == 1) {betaa0 = -3.330294;}
        break;
      case 50:                //Sn
        if (L == 0) {betaa0 = -8.621087;}
        else if (L == 1) {betaa0 = -4.989752;}
        break;
      case 51:                //Sb
        if (L == 0) {betaa0 = -7.472322;}
        else if (L == 1) {betaa0 = -5.940750;}
        else if (L == 2) {betaa0 = -3.979108;}
        break;
      case 52:                //Te
        if (L == 0) {betaa0 = -70.001062;}
        else if (L == 1) {betaa0 = -6.151642;}
        break;
      case 53:                //I
        if (L == 0) {betaa0 = -30.522481;}
        else if (L == 1) {betaa0 = -5.942120;}
        else if (L == 2) {betaa0 = -7.676107;}
        break;
      case 54:                //Xe
        if (L == 0) {betaa0 = -3.980622;}
        else if (L == 1) {betaa0 = -38.822792;}
        break;
      case 55:                //Cs
        if (L == 0) {betaa0 = 2.287838;}
        else if (L == 1) {betaa0 = -5.908071;}
        break;
      case 56:                //Ba
        if (L == 0) {betaa0 = 10.003125;}
        else if (L == 1) {betaa0 = -6.335160;}
        break;
      case 57:                //La
        if (L == 0) {betaa0 = 0.796727;}
        else if (L == 1) {betaa0 = -10.856056;}
        else if (L == 2) {betaa0 = -0.484922;}
        break;
      case 71:                //Lu
        if (L == 0) {betaa0 = -5.590778;}
        else if (L == 1) {betaa0 = -0.937679;}
        else if (L == 2) {betaa0 = -7.737752;}
        break;
      case 72:                //Hf
        if (L == 0) {betaa0 = -5.366351;}
        else if (L == 1) {betaa0 = -21.550119;}
        else if (L == 2) {betaa0 = -3.884443;}
        break;
      case 73:                //Ta
        if (L == 0) {betaa0 = -17.199605;}
        else if (L == 1) {betaa0 = -5.818839;}
        else if (L == 2) {betaa0 = -9.816794;}
        break;
      case 74:                //W
        if (L == 0) {betaa0 = -16.946460;}
        else if (L == 1) {betaa0 = 5.623170;}
        else if (L == 2) {betaa0 = -2.947340;}
        break;
      case 75:                //Re
        if (L == 0) {betaa0 = 3.830075;}
        else if (L == 1) {betaa0 = -1.638530;}
        else if (L == 2) {betaa0 = -1.414411;}
        break;
      case 76:                //Os
        if (L == 0) {betaa0 = -12.508730;}
        else if (L == 1) {betaa0 = 0.846880;}
        else if (L == 2) {betaa0 = 5.164360;}
        break;
      case 77:                //Ir
        if (L == 0) {betaa0 = -10.943427;}
        else if (L == 1) {betaa0 = 2.908880;}
        else if (L == 2) {betaa0 = -3.791731;}
        break;
      case 78:                //Pt
        if (L == 0) {betaa0 = 1.151418;}
        else if (L == 1) {betaa0 = 3.298694;}
        else if (L == 2) {betaa0 = -18.044737;}
        break;
      case 79:                //Au
        if (L == 0) {betaa0 = -7.479625;}
        else if (L == 1) {betaa0 = 3.664356;}
        else if (L == 2) {betaa0 = -61.715468;}
        break;
      case 80:                //Hg
        if (L == 0) {betaa0 = -3.045239;}
        else if (L == 1) {betaa0 = -5.693556;}
        break;
      case 81:                //Tl
        if (L == 0) {betaa0 = -7.230170;}
        else if (L == 1) {betaa0 = -7.575544;}
        break;
      case 82:                //Pb
        if (L == 0) {betaa0 = -8.323792;}
        else if (L == 1) {betaa0 = -2.237891;}
        break;
      case 83:                //Bi
        if (L == 0) {betaa0 = -34.951578;}
        else if (L == 1) {betaa0 = -7.359060;}
        break;
    }
    return betaa0/au2eV;
  }
  double UlX(size_t atomicnr, int L) {
    //function that returns the atomic potential U; values stored in eV, returned in a.u.
    double ulx = 0.0;
    switch (atomicnr) {
      case 1:                 //H
        ulx = -11.246958;
        break;
      case 2:                 //He
        if (L == 0) {ulx = -31.770969;}
        else if (L == 1) {ulx = -5.856382;}
        break;
      case 3:                 //Li
        if (L == 0) {ulx = -4.709912;}
        else if (L == 1) {ulx = -2.722581;}
        break;
      case 4:                 //Be
        if (L == 0) {ulx = -16.360315;}
        else if (L == 1) {ulx = -16.339216;}
        break;
      case 5:                 //B
        if (L == 0) {ulx = -25.967679;}
        else if (L == 1) {ulx = -19.115864;}
        break;
      case 6:                 //C
        if (L == 0) {ulx = -51.089653;}
        else if (L == 1) {ulx = -39.937920;}
        break;
      case 7:                 //N
        if (L == 0) {ulx = -57.784823;}
        else if (L == 1) {ulx = -49.893036;}
        break;
      case 8:                 //O
        if (L == 0) {ulx = -91.678761;}
        else if (L == 1) {ulx = -70.460949;}
        break;
      case 9:                 //F
        if (L == 0) {ulx = -140.225626;}
        else if (L == 1) {ulx = -98.778044;}
        break;
      case 10:                //Ne
        if (L == 0) {ulx = -2.978729;}
        else if (L == 1) {ulx = -85.441118;}
        break;
      case 11:                //Na
        if (L == 0) {ulx = -4.537153;}
        else if (L == 1) {ulx = -2.433015;}
        break;
      case 12:                //Mg
        if (L == 0) {ulx = -14.574226;}
        else if (L == 1) {ulx = -7.583850;}
        break;
      case 13:                //Al
        if (L == 0) {ulx = -24.546778;}
        else if (L == 1) {ulx = -20.104434;}
        else if (L == 2) {ulx = 8.004394;}
        break;
      case 14:                //Si
        if (L == 0) {ulx = -27.358058;}
        else if (L == 1) {ulx = -20.490578;}
        else if (L == 2) {ulx = -22.751900;}
        break;
      case 15:                //P
        if (L == 0) {ulx = -48.729905;}
        else if (L == 1) {ulx = -40.354689;}
        else if (L == 2) {ulx = -7.349246;}
        break;
      case 16:                //S
        if (L == 0) {ulx = -47.530706;}
        else if (L == 1) {ulx = -39.191045;}
        else if (L == 2) {ulx = -46.306944;}
        break;
      case 17:                //Cl
        if (L == 0) {ulx = -61.389930;}
        else if (L == 1) {ulx = -54.482801;}
        else if (L == 2) {ulx = -38.258155;}
        break;
      case 18:                //Ar
        if (L == 0) {ulx = -7.797931;}
        else if (L == 1) {ulx = -83.211487;}
        break;
      case 19:                //K
        if (L == 0) {ulx = -3.801108;}
        else if (L == 1) {ulx = -3.339656;}
        break;
      case 20:                //Ca
        if (L == 0) {ulx = -10.770058;}
        else if (L == 1) {ulx = -9.754177;}
        break;
      case 21:                //Sc
        if (L == 0) {ulx = -15.544461;}
        else if (L == 1) {ulx = -18.646295;}
        else if (L == 2) {ulx = -16.069444;}
        break;
      case 22:                //Ti
        if (L == 0) {ulx = -25.507973;}
        else if (L == 1) {ulx = -17.260909;}
        else if (L == 2) {ulx = -23.809486;}
        break;
      case 23:                //V
        if (L == 0) {ulx = -32.162276;}
        else if (L == 1) {ulx = -21.572501;}
        else if (L == 2) {ulx = -34.506245;}
        break;
      case 24:                //Cr
        if (L == 0) {ulx = -34.864339;}
        else if (L == 1) {ulx = -26.978615;}
        else if (L == 2) {ulx = -54.431036;}
        break;
      case 25:                //Mn
        if (L == 0) {ulx = -51.460000;}
        else if (L == 1) {ulx = -37.543990;}
        else if (L == 2) {ulx = -47.655370;}
        break;
      case 26:                //Fe
        if (L == 0) {ulx = -70.515047;}
        else if (L == 1) {ulx = -62.963069;}
        else if (L == 2) {ulx = -103.631790;}
        break;
      case 27:                //Co
        if (L == 0) {ulx = -21.039413;}
        else if (L == 1) {ulx = 10.000000;}
        else if (L == 2) {ulx = -28.068971;}
        break;
      case 28:                //Ni
        if (L == 0) {ulx = -47.620247;}
        else if (L == 1) {ulx = -32.878408;}
        else if (L == 2) {ulx = -93.026395;}
        break;
      case 29:                //Cu
        if (L == 0) {ulx = -97.002205;}
        else if (L == 1) {ulx = -1.000000;}
        else if (L == 2) {ulx = -110.442592;}
        break;
      case 30:                //Zn
        if (L == 0) {ulx = -18.040862;}
        else if (L == 1) {ulx = -7.834895;}
        break;
      case 31:                //Ga
        if (L == 0) {ulx = -30.600226;}
        else if (L == 1) {ulx = -21.032425;}
        break;
      case 32:                //Ge
        if (L == 0) {ulx = -32.747338;}
        else if (L == 1) {ulx = -24.709016;}
        break;
      case 33:                //As
        if (L == 0) {ulx = -37.956965;}
        else if (L == 1) {ulx = -38.453701;}
        else if (L == 2) {ulx = -30.282658;}
        break;
      case 34:                //Se
        if (L == 0) {ulx = -32.671088;}
        else if (L == 1) {ulx = -32.522220;}
        break;
      case 35:                //Br
        if (L == 0) {ulx = -45.834364;}
        else if (L == 1) {ulx = -50.293675;}
        else if (L == 2) {ulx = 7.086738;}
        break;
      case 36:                //Kr
        if (L == 0) {ulx = 8.535384;}
        else if (L == 1) {ulx = -80.484321;}
        break;
      case 37:                //Rb
        if (L == 0) {ulx = -3.636505;}
        else if (L == 1) {ulx = -2.500671;}
        break;
      case 38:                //Sr
        if (L == 0) {ulx = -10.427671;}
        else if (L == 1) {ulx = -9.943751;}
        break;
      case 39:                //Y
        if (L == 0) {ulx = -14.247809;}
        else if (L == 1) {ulx = -14.817140;}
        else if (L == 2) {ulx = -16.394302;}
        break;
      case 40:                //Zr
        if (L == 0) {ulx = -20.008884;}
        else if (L == 1) {ulx = -14.559692;}
        else if (L == 2) {ulx = -21.302657;}
        break;
      case 41:                //Nb
        if (L == 0) {ulx = -31.269298;}
        else if (L == 1) {ulx = -20.151277;}
        else if (L == 2) {ulx = -35.893116;}
        break;
      case 42:                //Mo
        if (L == 0) {ulx = -53.467728;}
        else if (L == 1) {ulx = -35.291951;}
        else if (L == 2) {ulx = -55.836977;}
        break;
      case 43:                //Tc
        if (L == 0) {ulx = -41.850292;}
        else if (L == 1) {ulx = -34.910293;}
        else if (L == 2) {ulx = -45.530412;}
        break;
      case 44:                //Ru
        if (L == 0) {ulx = -44.901521;}
        else if (L == 1) {ulx = -41.424409;}
        else if (L == 2) {ulx = -37.934514;}
        break;
      case 45:                //Rh
        if (L == 0) {ulx = -20.513756;}
        else if (L == 1) {ulx = -40.045431;}
        else if (L == 2) {ulx = -35.818492;}
        break;
      case 46:                //Pd
        if (L == 0) {ulx = -76.140196;}
        else if (L == 1) {ulx = -21.073362;}
        else if (L == 2) {ulx = -85.325301;}
        break;
      case 47:                //Ag
        if (L == 0) {ulx = -25.484137;}
        else if (L == 2) {ulx = -36.116023;}
        else if (L == 1) {ulx = -35.668272;}
        break;
      case 48:                //Cd
        if (L == 0) {ulx = -14.645792;}
        else if (L == 1) {ulx = -9.318664;}
        break;
      case 49:                //In
        if (L == 0) {ulx = -28.339246;}
        else if (L == 1) {ulx = -23.373875;}
        break;
      case 50:                //Sn
        if (L == 0) {ulx = -29.888217;}
        else if (L == 1) {ulx = -22.156954;}
        break;
      case 51:                //Sb
        if (L == 0) {ulx = -41.688879;}
        else if (L == 1) {ulx = -39.541180;}
        else if (L == 2) {ulx = -6.581663;}
        break;
      case 52:                //Te
        if (L == 0) {ulx = -114.733316;}
        else if (L == 1) {ulx = -50.096389;}
        break;
      case 53:                //I
        if (L == 0) {ulx = -59.973232;}
        else if (L == 1) {ulx = -56.459835;}
        else if (L == 2) {ulx = -28.822603;}
        break;
      case 54:                //Xe
        if (L == 0) {ulx = -18.270227;}
        else if (L == 1) {ulx = -167.163063;}
        break;
      case 55:                //Cs
        if (L == 0) {ulx = -3.748609;}
        else if (L == 1) {ulx = -2.348109;}
        break;
      case 56:                //Ba
        if (L == 0) {ulx = -9.306985;}
        else if (L == 1) {ulx = -8.826713;}
        break;
      case 57:                //La
        if (L == 0) {ulx = -19.641953;}
        else if (L == 1) {ulx = -22.059431;}
        else if (L == 2) {ulx = -22.638986;}
        break;
      case 71:                //Lu
        if (L == 0) {ulx = -15.954994;}
        else if (L == 1) {ulx = -11.606213;}
        else if (L == 2) {ulx = -13.050056;}
        break;
      case 72:                //Hf
        if (L == 0) {ulx = -22.375140;}
        else if (L == 1) {ulx = -13.081670;}
        else if (L == 2) {ulx = -20.637741;}
        break;
      case 73:                //Ta
        if (L == 0) {ulx = -39.009984;}
        else if (L == 1) {ulx = 1.163975;}
        else if (L == 2) {ulx = -43.266315;}
        break;
      case 74:                //W
        if (L == 0) {ulx = -44.524950;}
        else if (L == 1) {ulx = -40.011500;}
        else if (L == 2) {ulx = -46.490410;}
        break;
      case 75:                //Re
        if (L == 0) {ulx = -41.291342;}
        else if (L == 1) {ulx = -35.089592;}
        else if (L == 2) {ulx = -44.178985;}
        break;
      case 76:                //Os
        if (L == 0) {ulx = -26.434080;}
        else if (L == 1) {ulx = -48.739500;}
        else if (L == 2) {ulx = -55.837880;}
        break;
      case 77:                //Ir
        if (L == 0) {ulx = -29.703974;}
        else if (L == 1) {ulx = -38.210924;}
        else if (L == 2) {ulx = -32.538202;}
        break;
      case 78:                //Pt
        if (L == 0) {ulx = -73.516173;}
        else if (L == 1) {ulx = -68.320056;}
        else if (L == 2) {ulx = -76.598873;}
        break;
      case 79:                //Au
        if (L == 0) {ulx = -95.041846;}
        else if (L == 1) {ulx = -63.890158;}
        else if (L == 2) {ulx = -88.066087;}
        break;
      case 80:                //Hg
        if (L == 0) {ulx = -17.608732;}
        else if (L == 1) {ulx = -18.369417;}
        break;
      case 81:                //Tl
        if (L == 0) {ulx = -29.518621;}
        else if (L == 1) {ulx = -29.826907;}
        break;
      case 82:                //Pb
        if (L == 0) {ulx = -35.038145;}
        else if (L == 1) {ulx = -25.413401;}
        break;
      case 83:                //Bi
        if (L == 0) {ulx = -42.409177;}
        else if (L == 1) {ulx = -36.393746;}
        break;
    }
    return ulx/au2eV;
  }
  double alpha(int atm1, int atm2) {
    //function that returns alpha values for PM6; values stored in 1/Angstrom, returned in 1/Bohr
    double alpha = 0.0;
    int atmA = std::min(atm1,atm2);
    int atmB = std::max(atm1,atm2);
    switch (atmA) {
      case 1:                //H
        switch (atmB) {
          case 1:                         //H
            alpha = 3.5409419537;
            break;
          case 2:                         //He
            alpha = 2.9898810387;
            break;
          case 3:                         //Li
            alpha = 2.1362650394;
            break;
          case 4:                         //Be
            alpha = 2.4754180908;
            break;
          case 5:                         //B
            alpha = 2.6152310371;
            break;
          case 6:                         //C
            alpha = 1.0278060000;
            break;
          case 7:                         //N
            alpha = 0.9694060087;
            break;
          case 8:                         //O
            alpha = 1.2609419823;
            break;
          case 9:                         //F
            alpha = 3.1367399693;
            break;
          case 10:                        //Ne
            alpha = 5.9996800423;
            break;
          case 11:                        //Na
            alpha = 0.5003259778;
            break;
          case 12:                        //Mg
            alpha = 2.6515939236;
            break;
          case 13:                        //Al
            alpha = 2.0259959698;
            break;
          case 14:                        //Si
            alpha = 1.8969500065;
            break;
          case 15:                        //P
            alpha = 1.9265370369;
            break;
          case 16:                        //S
            alpha = 2.2159750462;
            break;
          case 17:                        //Cl
            alpha = 2.4028859138;
            break;
          case 18:                        //Ar
            alpha = 4.0561671257;
            break;
          case 19:                        //K
            alpha = 0.6481729746;
            break;
          case 20:                        //Ca
            alpha = 2.1418590546;
            break;
          case 21:                        //Sc
            alpha = 1.1794849634;
            break;
          case 22:                        //Ti
            alpha = 0.8326690197;
            break;
          case 23:                        //V
            alpha = 1.2801330090;
            break;
          case 24:                        //Cr
            alpha = 0.8826609850;
            break;
          case 25:                        //Mn
            alpha = 2.3099400997;
            break;
          case 26:                        //Fe
            alpha = 0.8544880152;
            break;
          case 27:                        //Co
            alpha = 2.9665179253;
            break;
          case 28:                        //Ni
            alpha = 2.6352798939;
            break;
          case 29:                        //Cu
            alpha = 2.3353590965;
            break;
          case 30:                        //Zn
            alpha = 1.9878909588;
            break;
          case 31:                        //Ga
            alpha = 1.8473500013;
            break;
          case 32:                        //Ge
            alpha = 2.2067930698;
            break;
          case 33:                        //As
            alpha = 1.9935270548;
            break;
          case 34:                        //Se
            alpha = 2.0350680351;
            break;
          case 35:                        //Br
            alpha = 2.1928029060;
            break;
          case 36:                        //Kr
            alpha = 3.7704529762;
            break;
          case 37:                        //Rb
            alpha = 2.4435560703;
            break;
          case 38:                        //Sr
            alpha = 2.1059141159;
            break;
          case 39:                        //Y
            alpha = 1.1890530586;
            break;
          case 40:                        //Zr
            alpha = 1.3797030449;
            break;
          case 41:                        //Nb
            alpha = 2.5059120655;
            break;
          case 42:                        //Mo
            alpha = 2.0357480049;
            break;
          case 43:                        //Tc
            alpha = 2.8303449154;
            break;
          case 44:                        //Ru
            alpha = 2.8928990364;
            break;
          case 45:                        //Rh
            alpha = 3.1041650772;
            break;
          case 46:                        //Pd
            alpha = 2.1837608814;
            break;
          case 47:                        //Ag
            alpha = 2.8959360123;
            break;
          case 48:                        //Cd
            alpha = 2.6287479401;
            break;
          case 49:                        //In
            alpha = 3.0641438961;
            break;
          case 50:                        //Sn
            alpha = 2.6489100456;
            break;
          case 51:                        //Sb
            alpha = 1.5712720156;
            break;
          case 52:                        //Te
            alpha = 2.0391299725;
            break;
          case 53:                        //I
            alpha = 2.1399130821;
            break;
          case 54:                        //Xe
            alpha = 1.3568609953;
            break;
          case 55:                        //Cs
            alpha = 0.2648819983;
            break;
          case 56:                        //Ba
            alpha = 6.0001349449;
            break;
          case 57:                        //La
            alpha = 0.8336669803;
            break;
          case 71:                        //Lu
            alpha = 1.4157899618;
            break;
          case 72:                        //Hf
            alpha = 1.4237879515;
            break;
          case 73:                        //Ta
            alpha = 2.2880139351;
            break;
          case 74:                        //W
            alpha = 2.1308801174;
            break;
          case 75:                        //Re
            alpha = 1.6345000267;
            break;
          case 76:                        //Os
            alpha = 3.4041800499;
            break;
          case 77:                        //Ir
            alpha = 1.0339000225;
            break;
          case 78:                        //Pt
            alpha = 4.0011978149;
            break;
          case 79:                        //Au
            alpha = 3.3690409660;
            break;
          case 80:                        //Hg
            alpha = 1.1365870237;
            break;
          case 81:                        //Tl
            alpha = 0.6736580133;
            break;
          case 82:                        //Pb
            alpha = 1.5226759911;
            break;
          case 83:                        //Bi
            alpha = 1.6799050570;
            break;
        }
        break;
      case 2:                //He
        switch (atmB) {
          case 2:                         //He
            alpha = 3.7835590839;
            break;
          case 3:                         //Li
            alpha = 3.1124029160;
            break;
          case 4:                         //Be
            alpha = 3.3067018986;
            break;
          case 5:                         //B
            alpha = 3.1631400585;
            break;
          case 6:                         //C
            alpha = 3.0427050591;
            break;
          case 7:                         //N
            alpha = 2.8143389225;
            break;
          case 8:                         //O
            alpha = 3.6537749767;
            break;
          case 9:                         //F
            alpha = 2.8565430641;
            break;
          case 10:                        //Ne
            alpha = 3.6777579784;
            break;
          case 11:                        //Na
            alpha = 1.7030290365;
            break;
          case 12:                        //Mg
            alpha = 2.2106029987;
            break;
          case 13:                        //Al
            alpha = 2.2558300495;
            break;
          case 14:                        //Si
            alpha = 2.0404980183;
            break;
          case 15:                        //P
            alpha = 2.0931580067;
            break;
          case 16:                        //S
            alpha = 1.9591490030;
            break;
          case 17:                        //Cl
            alpha = 1.6716769934;
            break;
          case 18:                        //Ar
            alpha = 2.7165620327;
            break;
          case 19:                        //K
            alpha = 1.4185010195;
            break;
          case 20:                        //Ca
            alpha = 1.7198469639;
            break;
          case 35:                        //Br
            alpha = 2.1282749176;
            break;
          case 36:                        //Kr
            alpha = 1.9969429970;
            break;
          case 37:                        //Rb
            alpha = 1.2707409859;
            break;
          case 53:                        //I
            alpha = 2.1729838848;
            break;
          case 54:                        //Xe
            alpha = 2.4978320599;
            break;
        }
        break;
      case 3:                     //Li
        switch (atmB) {
          case 3:                         //Li
            alpha = 4.7146739960;
            break;
          case 4:                         //Be
            alpha = 2.2367279530;
            break;
          case 5:                         //B
            alpha = 3.7593970299;
            break;
          case 6:                         //C
            alpha = 3.2418739796;
            break;
          case 7:                         //N
            alpha = 2.6406230927;
            break;
          case 8:                         //O
            alpha = 2.5844419003;
            break;
          case 9:                         //F
            alpha = 3.0439009666;
            break;
          case 10:                        //Ne
            alpha = 2.1936659813;
            break;
          case 11:                        //Na
            alpha = 1.2672990561;
            break;
          case 12:                        //Mg
            alpha = 1.1843800545;
            break;
          case 13:                        //Al
            alpha = 1.5815930367;
            break;
          case 14:                        //Si
            alpha = 1.7896089554;
            break;
          case 15:                        //P
            alpha = 1.3945440054;
            break;
          case 16:                        //S
            alpha = 2.2942750454;
            break;
          case 17:                        //Cl
            alpha = 2.7830009460;
            break;
          case 18:                        //Ar
            alpha = 3.1228950024;
            break;
          case 19:                        //K
            alpha = 1.0364869833;
            break;
          case 35:                        //Br
            alpha = 2.0744409561;
            break;
          case 36:                        //Kr
            alpha = 3.3145620823;
            break;
          case 53:                        //I
            alpha = 2.1212511063;
            break;
          case 54:                        //Xe
            alpha = 2.4668951035;
            break;
          case 82:                        //Pb
            alpha = 1.0018099546;
            break;
          case 83:                        //Bi
            alpha = 0.3401400149;
            break;
        }
        break;
      case 4:                     //Be
        switch (atmB) {
          case 4:                         //Be
            alpha = 1.4999070168;
            break;
          case 5:                         //B
            alpha = 1.8889980316;
            break;
          case 6:                         //C
            alpha = 4.2128820419;
            break;
          case 7:                         //N
            alpha = 2.5808949471;
            break;
          case 8:                         //O
            alpha = 3.0518670082;
            break;
          case 9:                         //F
            alpha = 3.7269229889;
            break;
          case 10:                        //Ne
            alpha = 1.3165880442;
            break;
          case 11:                        //Na
            alpha = 1.2554800510;
            break;
          case 12:                        //Mg
            alpha = 1.5575909615;
            break;
          case 13:                        //Al
            alpha = 1.9382369518;
            break;
          case 14:                        //Si
            alpha = 1.2631319761;
            break;
          case 15:                        //P
            alpha = 1.8000700474;
            break;
          case 16:                        //S
            alpha = 2.7817358971;
            break;
          case 17:                        //Cl
            alpha = 2.8226759434;
            break;
          case 18:                        //Ar
            alpha = 3.0440070629;
            break;
          case 19:                        //K
            alpha = 1.9318879843;
            break;
          case 35:                        //Br
            alpha = 2.3671460152;
            break;
          case 36:                        //Kr
            alpha = 3.2530479431;
            break;
          case 53:                        //I
            alpha = 2.2880229950;
            break;
          case 54:                        //Xe
            alpha = 6.0000028610;
            break;
        }
        break;
      case 5:                     //B
        switch (atmB) {
          case 5:                         //B
            alpha = 3.3186240196;
            break;
          case 6:                         //C
            alpha = 2.9190070629;
            break;
          case 7:                         //N
            alpha = 2.4770040512;
            break;
          case 8:                         //O
            alpha = 2.6953508854;
            break;
          case 9:                         //F
            alpha = 2.8238370419;
            break;
          case 10:                        //Ne
            alpha = 2.7561900616;
            break;
          case 11:                        //Na
            alpha = 1.5699609518;
            break;
          case 12:                        //Mg
            alpha = 2.5274410248;
            break;
          case 13:                        //Al
            alpha = 2.0595688820;
            break;
          case 14:                        //Si
            alpha = 1.9826530218;
            break;
          case 15:                        //P
            alpha = 1.9231679440;
            break;
          case 16:                        //S
            alpha = 2.4036960602;
            break;
          case 17:                        //Cl
            alpha = 2.2593228817;
            break;
          case 18:                        //Ar
            alpha = 2.4154710770;
            break;
          case 19:                        //K
            alpha = 2.0317680836;
            break;
          case 20:                        //Ca
            alpha = 1.7000099421;
            break;
          case 22:                        //Ti
            alpha = 1.6287100315;
            break;
          case 35:                        //Br
            alpha = 2.3078899384;
            break;
          case 36:                        //Kr
            alpha = 2.3631689548;
            break;
          case 37:                        //Rb
            alpha = 5.5322389603;
            break;
          case 53:                        //I
            alpha = 2.6676049232;
            break;
          case 54:                        //Xe
            alpha = 5.0519571304;
            break;
          case 55:                        //Cs
            alpha = 1.4871100187;
            break;
          case 72:                        //Hf
            alpha = 1.6334999800;
            break;
          case 81:                        //Tl
            alpha = 1.5283470154;
            break;
          case 82:                        //Pb
            alpha = 0.9111970067;
            break;
        }
        break;
      case 6:                     //C
        switch (atmB) {
          case 6:                         //C
            alpha = 2.6137130260;
            break;
          case 7:                         //N
            alpha = 2.6861081123;
            break;
          case 8:                         //O
            alpha = 2.8896069527;
            break;
          case 9:                         //F
            alpha = 3.0276000500;
            break;
          case 10:                        //Ne
            alpha = 3.4411880970;
            break;
          case 11:                        //Na
            alpha = 2.1960499287;
            break;
          case 12:                        //Mg
            alpha = 3.0409460068;
            break;
          case 13:                        //Al
            alpha = 2.2674400806;
            break;
          case 14:                        //Si
            alpha = 1.9844980240;
            break;
          case 15:                        //P
            alpha = 1.9946529865;
            break;
          case 16:                        //S
            alpha = 2.2103049755;
            break;
          case 17:                        //Cl
            alpha = 2.1621971130;
            break;
          case 18:                        //Ar
            alpha = 1.4713089466;
            break;
          case 19:                        //K
            alpha = 2.2417569160;
            break;
          case 20:                        //Ca
            alpha = 1.0353050232;
            break;
          case 21:                        //Sc
            alpha = 2.6304900646;
            break;
          case 22:                        //Ti
            alpha = 1.5979729891;
            break;
          case 23:                        //V
            alpha = 2.7898550034;
            break;
          case 24:                        //Cr
            alpha = 3.6567540169;
            break;
          case 25:                        //Mn
            alpha = 3.0007500648;
            break;
          case 26:                        //Fe
            alpha = 3.9913430214;
            break;
          case 27:                        //Co
            alpha = 3.7162330151;
            break;
          case 28:                        //Ni
            alpha = 4.2855129242;
            break;
          case 29:                        //Cu
            alpha = 4.6387729645;
            break;
          case 30:                        //Zn
            alpha = 1.8023270369;
            break;
          case 31:                        //Ga
            alpha = 2.3254098892;
            break;
          case 32:                        //Ge
            alpha = 2.2574689388;
            break;
          case 33:                        //As
            alpha = 1.8550690413;
            break;
          case 34:                        //Se
            alpha = 2.3871181011;
            break;
          case 35:                        //Br
            alpha = 2.0150859356;
            break;
          case 36:                        //Kr
            alpha = 2.0767381191;
            break;
          case 37:                        //Rb
            alpha = 2.7658300400;
            break;
          case 38:                        //Sr
            alpha = 1.9866880178;
            break;
          case 39:                        //Y
            alpha = 1.3360940218;
            break;
          case 40:                        //Zr
            alpha = 2.0294270515;
            break;
          case 41:                        //Nb
            alpha = 2.6210119724;
            break;
          case 42:                        //Mo
            alpha = 2.1986720562;
            break;
          case 43:                        //Tc
            alpha = 3.1983261108;
            break;
          case 44:                        //Ru
            alpha = 2.7848329544;
            break;
          case 45:                        //Rh
            alpha = 3.4159910679;
            break;
          case 46:                        //Pd
            alpha = 4.7771921158;
            break;
          case 47:                        //Ag
            alpha = 4.4043359756;
            break;
          case 48:                        //Cd
            alpha = 1.4256780148;
            break;
          case 49:                        //In
            alpha = 2.1892719269;
            break;
          case 50:                        //Sn
            alpha = 2.4405379295;
            break;
          case 51:                        //Sb
            alpha = 1.6962059736;
            break;
          case 52:                        //Te
            alpha = 1.9928159714;
            break;
          case 53:                        //I
            alpha = 2.0687100887;
            break;
          case 54:                        //Xe
            alpha = 1.7044399977;
            break;
          case 55:                        //Cs
            alpha = 2.1471040249;
            break;
          case 56:                        //Ba
            alpha = 0.7706260085;
            break;
          case 57:                        //La
            alpha = 0.6048690081;
            break;
          case 71:                        //Lu
            alpha = 2.3128130436;
            break;
          case 72:                        //Hf
            alpha = 1.0021940470;
            break;
          case 73:                        //Ta
            alpha = 1.8389489651;
            break;
          case 74:                        //W
            alpha = 2.0974800587;
            break;
          case 75:                        //Re
            alpha = 2.3062849045;
            break;
          case 76:                        //Os
            alpha = 2.3364999294;
            break;
          case 77:                        //Ir
            alpha = 1.6902949810;
            break;
          case 78:                        //Pt
            alpha = 3.3067219257;
            break;
          case 79:                        //Au
            alpha = 4.5800161362;
            break;
          case 80:                        //Hg
            alpha = 0.7958160043;
            break;
          case 81:                        //Tl
            alpha = 1.3903449774;
            break;
          case 82:                        //Pb
            alpha = 1.5255930424;
            break;
          case 83:                        //Bi
            alpha = 1.5340249538;
            break;
        }
        break;
      case 7:                     //N
        switch (atmB) {
          case 7:                         //N
            alpha = 2.5745019913;
            break;
          case 8:                         //O
            alpha = 2.7842919827;
            break;
          case 9:                         //F
            alpha = 2.8566460609;
            break;
          case 10:                        //Ne
            alpha = 4.4263701439;
            break;
          case 11:                        //Na
            alpha = 2.4943840504;
            break;
          case 12:                        //Mg
            alpha = 2.0791249275;
            break;
          case 13:                        //Al
            alpha = 2.0097539425;
            break;
          case 14:                        //Si
            alpha = 1.8189879656;
            break;
          case 15:                        //P
            alpha = 2.1470420361;
            break;
          case 16:                        //S
            alpha = 2.2899899483;
            break;
          case 17:                        //Cl
            alpha = 2.1721339226;
            break;
          case 18:                        //Ar
            alpha = 2.3268051147;
            break;
          case 19:                        //K
            alpha = 2.3258590698;
            break;
          case 20:                        //Ca
            alpha = 2.3866000175;
            break;
          case 21:                        //Sc
            alpha = 2.2700040340;
            break;
          case 22:                        //Ti
            alpha = 1.6786860228;
            break;
          case 23:                        //V
            alpha = 1.6075400114;
            break;
          case 24:                        //Cr
            alpha = 3.0291860104;
            break;
          case 25:                        //Mn
            alpha = 2.9214699268;
            break;
          case 26:                        //Fe
            alpha = 2.5004858971;
            break;
          case 27:                        //Co
            alpha = 3.6186380386;
            break;
          case 28:                        //Ni
            alpha = 3.8452150822;
            break;
          case 29:                        //Cu
            alpha = 4.2143368721;
            break;
          case 30:                        //Zn
            alpha = 1.8445789814;
            break;
          case 31:                        //Ga
            alpha = 2.1218199730;
            break;
          case 32:                        //Ge
            alpha = 1.9882260561;
            break;
          case 33:                        //As
            alpha = 1.4965430498;
            break;
          case 34:                        //Se
            alpha = 1.9377640486;
            break;
          case 35:                        //Br
            alpha = 4.2249011993;
            break;
          case 36:                        //Kr
            alpha = 1.6440520287;
            break;
          case 37:                        //Rb
            alpha = 0.7610470057;
            break;
          case 38:                        //Sr
            alpha = 2.1836290359;
            break;
          case 39:                        //Y
            alpha = 1.7787959576;
            break;
          case 40:                        //Zr
            alpha = 1.7070829868;
            break;
          case 41:                        //Nb
            alpha = 2.0238630772;
            break;
          case 42:                        //Mo
            alpha = 1.8694750071;
            break;
          case 43:                        //Tc
            alpha = 2.3154170513;
            break;
          case 44:                        //Ru
            alpha = 3.0555040836;
            break;
          case 45:                        //Rh
            alpha = 3.5854620934;
            break;
          case 46:                        //Pd
            alpha = 2.3280460835;
            break;
          case 47:                        //Ag
            alpha = 4.6598711014;
            break;
          case 48:                        //Cd
            alpha = 0.9704229832;
            break;
          case 49:                        //In
            alpha = 2.4698679447;
            break;
          case 50:                        //Sn
            alpha = 2.0855889320;
            break;
          case 51:                        //Sb
            alpha = 0.6761149764;
            break;
          case 52:                        //Te
            alpha = 1.7222690582;
            break;
          case 53:                        //I
            alpha = 1.6775180101;
            break;
          case 54:                        //Xe
            alpha = 1.9329520464;
            break;
          case 55:                        //Cs
            alpha = 2.4465320110;
            break;
          case 56:                        //Ba
            alpha = 1.1482330561;
            break;
          case 57:                        //La
            alpha = 0.7588809729;
            break;
          case 71:                        //Lu
            alpha = 2.1413021088;
            break;
          case 72:                        //Hf
            alpha = 1.3324099779;
            break;
          case 73:                        //Ta
            alpha = 2.0536789894;
            break;
          case 74:                        //W
            alpha = 1.5960400105;
            break;
          case 75:                        //Re
            alpha = 1.9183319807;
            break;
          case 76:                        //Os
            alpha = 1.1430900097;
            break;
          case 77:                        //Ir
            alpha = 3.9345080853;
            break;
          case 78:                        //Pt
            alpha = 2.3079230785;
            break;
          case 79:                        //Au
            alpha = 2.1380949020;
            break;
          case 80:                        //Hg
            alpha = 0.3321520090;
            break;
          case 81:                        //Tl
            alpha = 0.9823349714;
            break;
          case 82:                        //Pb
            alpha = 1.3173940182;
            break;
          case 83:                        //Bi
            alpha = 1.1438759565;
            break;
        }
        break;
      case 8:                     //O
        switch (atmB) {
          case 8:                         //O
            alpha = 2.6239979267;
            break;
          case 9:                         //F
            alpha = 3.0154440403;
            break;
          case 10:                        //Ne
            alpha = 2.8895869255;
            break;
          case 11:                        //Na
            alpha = 1.9814490080;
            break;
          case 12:                        //Mg
            alpha = 2.2515199184;
            break;
          case 13:                        //Al
            alpha = 2.4986600876;
            break;
          case 14:                        //Si
            alpha = 1.9235999584;
            break;
          case 15:                        //P
            alpha = 2.2207679749;
            break;
          case 16:                        //S
            alpha = 2.3832890987;
            break;
          case 17:                        //Cl
            alpha = 2.3232359886;
            break;
          case 18:                        //Ar
            alpha = 2.2406730652;
            break;
          case 19:                        //K
            alpha = 1.5085710287;
            break;
          case 20:                        //Ca
            alpha = 3.2638969421;
            break;
          case 21:                        //Sc
            alpha = 2.2565159798;
            break;
          case 22:                        //Ti
            alpha = 1.7891180515;
            break;
          case 23:                        //V
            alpha = 1.6239730120;
            break;
          case 24:                        //Cr
            alpha = 2.5000000000;
            break;
          case 25:                        //Mn
            alpha = 2.5775399208;
            break;
          case 26:                        //Fe
            alpha = 1.7263129950;
            break;
          case 27:                        //Co
            alpha = 3.7269110680;
            break;
          case 28:                        //Ni
            alpha = 2.9372320175;
            break;
          case 29:                        //Cu
            alpha = 3.9599509239;
            break;
          case 30:                        //Zn
            alpha = 2.3350539207;
            break;
          case 31:                        //Ga
            alpha = 2.3483469486;
            break;
          case 32:                        //Ge
            alpha = 2.1394131184;
            break;
          case 33:                        //As
            alpha = 2.0039501190;
            break;
          case 34:                        //Se
            alpha = 2.4842629433;
            break;
          case 35:                        //Br
            alpha = 2.2830460072;
            break;
          case 36:                        //Kr
            alpha = 0.2922999859;
            break;
          case 37:                        //Rb
            alpha = 1.3349080086;
            break;
          case 38:                        //Sr
            alpha = 2.1383988857;
            break;
          case 39:                        //Y
            alpha = 1.8510299921;
            break;
          case 40:                        //Zr
            alpha = 1.7095700502;
            break;
          case 41:                        //Nb
            alpha = 2.0494890213;
            break;
          case 42:                        //Mo
            alpha = 1.7554240227;
            break;
          case 43:                        //Tc
            alpha = 2.4051899910;
            break;
          case 44:                        //Ru
            alpha = 3.1349399090;
            break;
          case 45:                        //Rh
            alpha = 3.9278299809;
            break;
          case 46:                        //Pd
            alpha = 2.1548669338;
            break;
          case 47:                        //Ag
            alpha = 1.8938740492;
            break;
          case 48:                        //Cd
            alpha = 1.6966730356;
            break;
          case 49:                        //In
            alpha = 2.6620950699;
            break;
          case 50:                        //Sn
            alpha = 2.7272601128;
            break;
          case 51:                        //Sb
            alpha = 1.8463840485;
            break;
          case 52:                        //Te
            alpha = 1.8530639410;
            break;
          case 53:                        //I
            alpha = 2.2889189720;
            break;
          case 54:                        //Xe
            alpha = 0.8392329812;
            break;
          case 55:                        //Cs
            alpha = 2.0851390362;
            break;
          case 56:                        //Ba
            alpha = 1.2830179930;
            break;
          case 57:                        //La
            alpha = 1.3183330297;
            break;
          case 71:                        //Lu
            alpha = 2.1924860477;
            break;
          case 72:                        //Hf
            alpha = 1.6332889795;
            break;
          case 73:                        //Ta
            alpha = 2.4126288891;
            break;
          case 74:                        //W
            alpha = 1.3590199947;
            break;
          case 75:                        //Re
            alpha = 1.9677469730;
            break;
          case 76:                        //Os
            alpha = 1.3503600359;
            break;
          case 77:                        //Ir
            alpha = 3.7482719421;
            break;
          case 78:                        //Pt
            alpha = 2.1105630398;
            break;
          case 79:                        //Au
            alpha = 1.5487630367;
            break;
          case 80:                        //Hg
            alpha = 1.0521450043;
            break;
          case 81:                        //Tl
            alpha = 1.5500680208;
            break;
          case 82:                        //Pb
            alpha = 1.7632100582;
            break;
          case 83:                        //Bi
            alpha = 1.5532970428;
            break;
        }
        break;
      case 9:                     //F
        switch (atmB) {
          case 9:                         //F
            alpha = 3.1757590771;
            break;
          case 10:                        //Ne
            alpha = 3.6756110191;
            break;
          case 11:                        //Na
            alpha = 2.6195509434;
            break;
          case 12:                        //Mg
            alpha = 3.3622078896;
            break;
          case 13:                        //Al
            alpha = 3.0842580795;
            break;
          case 14:                        //Si
            alpha = 2.1310279369;
            break;
          case 15:                        //P
            alpha = 2.2343559265;
            break;
          case 16:                        //S
            alpha = 2.1871860027;
            break;
          case 17:                        //Cl
            alpha = 2.3132700920;
            break;
          case 18:                        //Ar
            alpha = 3.9206581116;
            break;
          case 19:                        //K
            alpha = 3.1828169823;
            break;
          case 20:                        //Ca
            alpha = 2.6450529099;
            break;
          case 21:                        //Sc
            alpha = 3.1079850197;
            break;
          case 22:                        //Ti
            alpha = 2.3070869446;
            break;
          case 23:                        //V
            alpha = 1.8251600266;
            break;
          case 24:                        //Cr
            alpha = 2.7165210247;
            break;
          case 25:                        //Mn
            alpha = 2.7919499874;
            break;
          case 26:                        //Fe
            alpha = 4.2947068214;
            break;
          case 27:                        //Co
            alpha = 3.9563469887;
            break;
          case 28:                        //Ni
            alpha = 3.4402410984;
            break;
          case 29:                        //Cu
            alpha = 4.4788317680;
            break;
          case 30:                        //Zn
            alpha = 2.4100210667;
            break;
          case 31:                        //Ga
            alpha = 2.6798689365;
            break;
          case 32:                        //Ge
            alpha = 2.3847770691;
            break;
          case 33:                        //As
            alpha = 2.0125830173;
            break;
          case 34:                        //Se
            alpha = 2.3021800518;
            break;
          case 35:                        //Br
            alpha = 2.0317649841;
            break;
          case 36:                        //Kr
            alpha = 3.4523210526;
            break;
          case 37:                        //Rb
            alpha = 3.6381220818;
            break;
          case 38:                        //Sr
            alpha = 3.0506660938;
            break;
          case 39:                        //Y
            alpha = 2.6480460167;
            break;
          case 40:                        //Zr
            alpha = 1.9009250402;
            break;
          case 41:                        //Nb
            alpha = 3.0031569004;
            break;
          case 42:                        //Mo
            alpha = 2.2025930882;
            break;
          case 43:                        //Tc
            alpha = 3.6048150063;
            break;
          case 44:                        //Ru
            alpha = 3.8787109852;
            break;
          case 45:                        //Rh
            alpha = 4.0516538620;
            break;
          case 46:                        //Pd
            alpha = 4.2373118401;
            break;
          case 47:                        //Ag
            alpha = 4.6284232140;
            break;
          case 48:                        //Cd
            alpha = 2.3121349812;
            break;
          case 49:                        //In
            alpha = 2.9487969875;
            break;
          case 50:                        //Sn
            alpha = 3.7242860794;
            break;
          case 51:                        //Sb
            alpha = 2.1829218864;
            break;
          case 52:                        //Te
            alpha = 1.9985760450;
            break;
          case 53:                        //I
            alpha = 2.2035799026;
            break;
          case 54:                        //Xe
            alpha = 1.1288119555;
            break;
          case 55:                        //Cs
            alpha = 2.8341000080;
            break;
          case 56:                        //Ba
            alpha = 3.0006179810;
            break;
          case 57:                        //La
            alpha = 2.3793349266;
            break;
          case 72:                        //Hf
            alpha = 2.2908029556;
            break;
          case 73:                        //Ta
            alpha = 3.1073899269;
            break;
          case 74:                        //W
            alpha = 1.4460500479;
            break;
          case 75:                        //Re
            alpha = 2.1542189121;
            break;
          case 76:                        //Os
            alpha = 1.5076199770;
            break;
          case 77:                        //Ir
            alpha = 2.9827990532;
            break;
          case 78:                        //Pt
            alpha = 3.7144410610;
            break;
          case 79:                        //Au
            alpha = 4.4531450272;
            break;
          case 80:                        //Hg
            alpha = 1.2405719757;
            break;
          case 81:                        //Tl
            alpha = 1.4695160389;
            break;
          case 82:                        //Pb
            alpha = 3.2889020443;
            break;
          case 83:                        //Bi
            alpha = 2.3554000854;
            break;
        }
        break;
      case 10:                    //Ne
        switch (atmB) {
          case 10:                        //Ne
            alpha = 3.9745669365;
            break;
          case 11:                        //Na
            alpha = 1.7742359638;
            break;
          case 12:                        //Mg
            alpha = 2.0316760540;
            break;
          case 13:                        //Al
            alpha = 2.4478690624;
            break;
          case 14:                        //Si
            alpha = 2.8677840233;
            break;
          case 15:                        //P
            alpha = 2.2190361023;
            break;
          case 16:                        //S
            alpha = 2.7870581150;
            break;
          case 17:                        //Cl
            alpha = 1.7031509876;
            break;
          case 18:                        //Ar
            alpha = 2.9637470245;
            break;
          case 19:                        //K
            alpha = 1.1380209923;
            break;
          case 20:                        //Ca
            alpha = 0.9545300007;
            break;
          case 35:                        //Br
            alpha = 2.4641718864;
            break;
          case 36:                        //Kr
            alpha = 2.8136789799;
            break;
          case 37:                        //Rb
            alpha = 2.2675909996;
            break;
          case 53:                        //I
            alpha = 2.4144148827;
            break;
          case 54:                        //Xe
            alpha = 1.3302019835;
            break;
        }
        break;
      case 11:                    //Na
        switch (atmB) {
          case 11:                        //Na
            alpha = 0.4464350045;
            break;
          case 12:                        //Mg
            alpha = 1.5067729950;
            break;
          case 13:                        //Al
            alpha = 1.2028709650;
            break;
          case 14:                        //Si
            alpha = 2.0076150894;
            break;
          case 15:                        //P
            alpha = 1.5003199577;
            break;
          case 16:                        //S
            alpha = 1.4008500576;
            break;
          case 17:                        //Cl
            alpha = 1.8164290190;
            break;
          case 18:                        //Ar
            alpha = 2.1676769257;
            break;
          case 19:                        //K
            alpha = 0.8843070269;
            break;
          case 20:                        //Ca
            alpha = 3.1071040630;
            break;
          case 23:                        //V
            alpha = 2.5510098934;
            break;
          case 24:                        //Cr
            alpha = 2.2950561047;
            break;
          case 35:                        //Br
            alpha = 1.6222180128;
            break;
          case 36:                        //Kr
            alpha = 2.4805979729;
            break;
          case 41:                        //Nb
            alpha = 2.5510098934;
            break;
          case 42:                        //Mo
            alpha = 2.4407699108;
            break;
          case 53:                        //I
            alpha = 1.4030900002;
            break;
          case 54:                        //Xe
            alpha = 2.1030030251;
            break;
          case 73:                        //Ta
            alpha = 2.5511200428;
            break;
          case 74:                        //W
            alpha = 2.5510299206;
            break;
          case 76:                        //Os
            alpha = 2.5507400036;
            break;
          case 77:                        //Ir
            alpha = 2.5508201122;
            break;
        }
        break;
      case 12:                    //Mg
        switch (atmB) {
          case 12:                        //Mg
            alpha = 1.0935729742;
            break;
          case 13:                        //Al
            alpha = 1.9725300074;
            break;
          case 14:                        //Si
            alpha = 3.1397490501;
            break;
          case 15:                        //P
            alpha = 1.3837729692;
            break;
          case 16:                        //S
            alpha = 1.5001629591;
            break;
          case 17:                        //Cl
            alpha = 2.3918058872;
            break;
          case 18:                        //Ar
            alpha = 2.0926640034;
            break;
          case 19:                        //K
            alpha = 0.8848099709;
            break;
          case 20:                        //Ca
            alpha = 2.2997999191;
            break;
          case 22:                        //Ti
            alpha = 1.9113399982;
            break;
          case 35:                        //Br
            alpha = 2.1956970692;
            break;
          case 36:                        //Kr
            alpha = 1.3914870024;
            break;
          case 53:                        //I
            alpha = 2.0451369286;
            break;
          case 54:                        //Xe
            alpha = 2.6984140873;
            break;
          case 72:                        //Hf
            alpha = 1.9113500118;
            break;
        }
        break;
      case 13:                    //Al
        switch (atmB) {
          case 13:                        //Al
            alpha = 1.3877140284;
            break;
          case 14:                        //Si
            alpha = 1.8999999762;
            break;
          case 15:                        //P
            alpha = 1.9807269573;
            break;
          case 16:                        //S
            alpha = 1.9767049551;
            break;
          case 17:                        //Cl
            alpha = 2.1259388924;
            break;
          case 18:                        //Ar
            alpha = 2.6451649666;
            break;
          case 19:                        //K
            alpha = 1.9760760069;
            break;
          case 20:                        //Ca
            alpha = 1.6125650406;
            break;
          case 21:                        //Sc
            alpha = 1.0035500526;
            break;
          case 22:                        //Ti
            alpha = 1.3694859743;
            break;
          case 25:                        //Mn
            alpha = 1.7683600187;
            break;
          case 33:                        //As
            alpha = 1.1527860165;
            break;
          case 35:                        //Br
            alpha = 1.8941409588;
            break;
          case 36:                        //Kr
            alpha = 2.4671308994;
            break;
          case 37:                        //Rb
            alpha = 0.7987740040;
            break;
          case 39:                        //Y
            alpha = 1.0034999847;
            break;
          case 40:                        //Zr
            alpha = 1.2706199884;
            break;
          case 46:                        //Pd
            alpha = 1.5727200508;
            break;
          case 47:                        //Ag
            alpha = 1.9287999868;
            break;
          case 51:                        //Sb
            alpha = 1.4226410389;
            break;
          case 52:                        //Te
            alpha = 1.3875410557;
            break;
          case 53:                        //I
            alpha = 1.8160680532;
            break;
          case 54:                        //Xe
            alpha = 2.4120390415;
            break;
          case 56:                        //Ba
            alpha = 2.1059238911;
            break;
          case 57:                        //La
            alpha = 1.0035099983;
            break;
          case 72:                        //Hf
            alpha = 0.9491500258;
            break;
          case 78:                        //Pt
            alpha = 1.5723600388;
            break;
          case 79:                        //Au
            alpha = 1.5725699663;
            break;
        }
        break;
      case 14:                    //Si
        switch (atmB) {
          case 14:                        //Si
            alpha = 1.3289999962;
            break;
          case 15:                        //P
            alpha = 3.3134660721;
            break;
          case 16:                        //S
            alpha = 1.8859159946;
            break;
          case 17:                        //Cl
            alpha = 1.6849780083;
            break;
          case 18:                        //Ar
            alpha = 1.7803499699;
            break;
          case 19:                        //K
            alpha = 1.6759300232;
            break;
          case 20:                        //Ca
            alpha = 1.2187880278;
            break;
          case 21:                        //Sc
            alpha = 2.0168700218;
            break;
          case 22:                        //Ti
            alpha = 2.8560380936;
            break;
          case 24:                        //Cr
            alpha = 1.8607599735;
            break;
          case 25:                        //Mn
            alpha = 1.9379589558;
            break;
          case 27:                        //Co
            alpha = 2.4698050022;
            break;
          case 28:                        //Ni
            alpha = 2.0688810349;
            break;
          case 30:                        //Zn
            alpha = 1.8320579529;
            break;
          case 31:                        //Ga
            alpha = 1.9137799740;
            break;
          case 32:                        //Ge
            alpha = 0.2997210026;
            break;
          case 33:                        //As
            alpha = 1.9155999422;
            break;
          case 34:                        //Se
            alpha = 1.5298169851;
            break;
          case 35:                        //Br
            alpha = 1.5708249807;
            break;
          case 36:                        //Kr
            alpha = 1.7640999556;
            break;
          case 38:                        //Sr
            alpha = 2.9697799683;
            break;
          case 39:                        //Y
            alpha = 2.0168199539;
            break;
          case 40:                        //Zr
            alpha = 1.7508330345;
            break;
          case 44:                        //Ru
            alpha = 2.7759099007;
            break;
          case 45:                        //Rh
            alpha = 2.7764899731;
            break;
          case 46:                        //Pd
            alpha = 2.9481999874;
            break;
          case 48:                        //Cd
            alpha = 1.3712249994;
            break;
          case 51:                        //Sb
            alpha = 2.6865899563;
            break;
          case 53:                        //I
            alpha = 1.5595790148;
            break;
          case 54:                        //Xe
            alpha = 3.0870599747;
            break;
          case 56:                        //Ba
            alpha = 1.2404199839;
            break;
          case 57:                        //La
            alpha = 2.0168199539;
            break;
          case 72:                        //Hf
            alpha = 2.1893000603;
            break;
          case 75:                        //Re
            alpha = 2.7759299278;
            break;
          case 78:                        //Pt
            alpha = 0.9999899864;
            break;
          case 80:                        //Hg
            alpha = 2.7708599567;
            break;
        }
        break;
      case 15:                    //P
        switch (atmB) {
          case 15:                        //P
            alpha = 1.5057920218;
            break;
          case 16:                        //S
            alpha = 1.5953249931;
            break;
          case 17:                        //Cl
            alpha = 1.4683059454;
            break;
          case 18:                        //Ar
            alpha = 4.3725161552;
            break;
          case 19:                        //K
            alpha = 1.4437379837;
            break;
          case 20:                        //Ca
            alpha = 1.0241420269;
            break;
          case 21:                        //Sc
            alpha = 0.8681650162;
            break;
          case 22:                        //Ti
            alpha = 2.1519289017;
            break;
          case 23:                        //V
            alpha = 2.5491540432;
            break;
          case 24:                        //Cr
            alpha = 1.6953829527;
            break;
          case 25:                        //Mn
            alpha = 1.9470200539;
            break;
          case 26:                        //Fe
            alpha = 2.5675339699;
            break;
          case 27:                        //Co
            alpha = 1.1525050402;
            break;
          case 28:                        //Ni
            alpha = 3.2602829933;
            break;
          case 29:                        //Cu
            alpha = 0.2106399983;
            break;
          case 30:                        //Zn
            alpha = 1.2204799652;
            break;
          case 31:                        //Ga
            alpha = 2.9796500206;
            break;
          case 32:                        //Ge
            alpha = 2.4692909718;
            break;
          case 34:                        //Se
            alpha = 1.0481829643;
            break;
          case 35:                        //Br
            alpha = 1.4021389484;
            break;
          case 38:                        //Sr
            alpha = 2.7891499996;
            break;
          case 39:                        //Y
            alpha = 0.9544500113;
            break;
          case 40:                        //Zr
            alpha = 1.0918580294;
            break;
          case 41:                        //Nb
            alpha = 2.2216079235;
            break;
          case 42:                        //Mo
            alpha = 1.8504409790;
            break;
          case 44:                        //Ru
            alpha = 0.2989160120;
            break;
          case 45:                        //Rh
            alpha = 2.3346068859;
            break;
          case 46:                        //Pd
            alpha = 0.8036299944;
            break;
          case 47:                        //Ag
            alpha = 6.0000061989;
            break;
          case 52:                        //Te
            alpha = 1.4537179470;
            break;
          case 53:                        //I
            alpha = 2.1315929890;
            break;
          case 55:                        //Cs
            alpha = 2.9249529839;
            break;
          case 57:                        //La
            alpha = 0.9544500113;
            break;
          case 71:                        //Lu
            alpha = 5.6188201904;
            break;
          case 72:                        //Hf
            alpha = 1.2312200069;
            break;
          case 73:                        //Ta
            alpha = 2.5137999058;
            break;
          case 74:                        //W
            alpha = 2.3380599022;
            break;
          case 75:                        //Re
            alpha = 1.8041679859;
            break;
          case 76:                        //Os
            alpha = 2.8360900879;
            break;
          case 77:                        //Ir
            alpha = 2.7140600681;
            break;
          case 78:                        //Pt
            alpha = 1.4032390118;
            break;
          case 79:                        //Au
            alpha = 1.6187130213;
            break;
          case 80:                        //Hg
            alpha = 0.6086040139;
            break;
          case 82:                        //Pb
            alpha = 4.5167999268;
            break;
        }
        break;
      case 16:                    //S
        switch (atmB) {
          case 16:                        //S
            alpha = 1.7945560217;
            break;
          case 17:                        //Cl
            alpha = 1.7154350281;
            break;
          case 18:                        //Ar
            alpha = 2.0493979454;
            break;
          case 19:                        //K
            alpha = 2.5121560097;
            break;
          case 20:                        //Ca
            alpha = 0.9581710100;
            break;
          case 21:                        //Sc
            alpha = 0.4229390025;
            break;
          case 22:                        //Ti
            alpha = 1.8464390039;
            break;
          case 23:                        //V
            alpha = 2.7041239738;
            break;
          case 24:                        //Cr
            alpha = 2.2609779835;
            break;
          case 25:                        //Mn
            alpha = 2.4825100899;
            break;
          case 26:                        //Fe
            alpha = 0.9889910221;
            break;
          case 27:                        //Co
            alpha = 2.4292550087;
            break;
          case 28:                        //Ni
            alpha = 2.0027520657;
            break;
          case 29:                        //Cu
            alpha = 0.2731119990;
            break;
          case 30:                        //Zn
            alpha = 1.4550000429;
            break;
          case 31:                        //Ga
            alpha = 2.2321081161;
            break;
          case 32:                        //Ge
            alpha = 2.0245881081;
            break;
          case 33:                        //As
            alpha = 1.9543679953;
            break;
          case 34:                        //Se
            alpha = 1.4796060324;
            break;
          case 35:                        //Br
            alpha = 1.5098739862;
            break;
          case 37:                        //Rb
            alpha = 1.3031840324;
            break;
          case 38:                        //Sr
            alpha = 1.5981060266;
            break;
          case 39:                        //Y
            alpha = 0.9716879725;
            break;
          case 40:                        //Zr
            alpha = 2.1297609806;
            break;
          case 41:                        //Nb
            alpha = 2.2494819164;
            break;
          case 42:                        //Mo
            alpha = 1.9396580458;
            break;
          case 43:                        //Tc
            alpha = 2.4634010792;
            break;
          case 44:                        //Ru
            alpha = 2.5080759525;
            break;
          case 45:                        //Rh
            alpha = 3.1540060043;
            break;
          case 46:                        //Pd
            alpha = 2.1778008938;
            break;
          case 47:                        //Ag
            alpha = 3.6531209946;
            break;
          case 48:                        //Cd
            alpha = 1.1822019815;
            break;
          case 49:                        //In
            alpha = 2.5421309471;
            break;
          case 50:                        //Sn
            alpha = 2.1315419674;
            break;
          case 51:                        //Sb
            alpha = 1.4188369513;
            break;
          case 52:                        //Te
            alpha = 1.8301700354;
            break;
          case 53:                        //I
            alpha = 1.8551100492;
            break;
          case 55:                        //Cs
            alpha = 0.2894119918;
            break;
          case 56:                        //Ba
            alpha = 0.7051879764;
            break;
          case 57:                        //La
            alpha = 1.8341289759;
            break;
          case 72:                        //Hf
            alpha = 2.3271100521;
            break;
          case 73:                        //Ta
            alpha = 2.2467229366;
            break;
          case 74:                        //W
            alpha = 1.5425699949;
            break;
          case 75:                        //Re
            alpha = 1.0839190483;
            break;
          case 76:                        //Os
            alpha = 2.8094999790;
            break;
          case 77:                        //Ir
            alpha = 3.2048339844;
            break;
          case 78:                        //Pt
            alpha = 2.7915000916;
            break;
          case 79:                        //Au
            alpha = 4.3062381744;
            break;
          case 80:                        //Hg
            alpha = 1.0416820049;
            break;
          case 81:                        //Tl
            alpha = 0.9948509932;
            break;
          case 82:                        //Pb
            alpha = 1.0275189877;
            break;
          case 83:                        //Bi
            alpha = 1.4668790102;
            break;
        }
        break;
      case 17:                    //Cl
        switch (atmB) {
          case 17:                        //Cl
            alpha = 1.8232389688;
            break;
          case 18:                        //Ar
            alpha = 2.5544490814;
            break;
          case 19:                        //K
            alpha = 1.6221630573;
            break;
          case 20:                        //Ca
            alpha = 2.3833909035;
            break;
          case 21:                        //Sc
            alpha = 2.1414740086;
            break;
          case 22:                        //Ti
            alpha = 1.4610340595;
            break;
          case 23:                        //V
            alpha = 1.6885290146;
            break;
          case 24:                        //Cr
            alpha = 2.1526179314;
            break;
          case 25:                        //Mn
            alpha = 1.6570099592;
            break;
          case 26:                        //Fe
            alpha = 1.2297929525;
            break;
          case 27:                        //Co
            alpha = 3.2174971104;
            break;
          case 28:                        //Ni
            alpha = 2.2005119324;
            break;
          case 29:                        //Cu
            alpha = 2.7765309811;
            break;
          case 30:                        //Zn
            alpha = 1.6251759529;
            break;
          case 31:                        //Ga
            alpha = 2.0247099400;
            break;
          case 32:                        //Ge
            alpha = 1.7712279558;
            break;
          case 33:                        //As
            alpha = 1.6910699606;
            break;
          case 34:                        //Se
            alpha = 2.1288609505;
            break;
          case 35:                        //Br
            alpha = 1.7103309631;
            break;
          case 36:                        //Kr
            alpha = 1.8849740028;
            break;
          case 37:                        //Rb
            alpha = 2.2744109631;
            break;
          case 38:                        //Sr
            alpha = 1.8541899920;
            break;
          case 39:                        //Y
            alpha = 1.6301519871;
            break;
          case 40:                        //Zr
            alpha = 1.3288350105;
            break;
          case 41:                        //Nb
            alpha = 2.2152750492;
            break;
          case 42:                        //Mo
            alpha = 1.7833620310;
            break;
          case 43:                        //Tc
            alpha = 2.5720429420;
            break;
          case 44:                        //Ru
            alpha = 1.7598830462;
            break;
          case 45:                        //Rh
            alpha = 3.3001298904;
            break;
          case 46:                        //Pd
            alpha = 3.8712430000;
            break;
          case 47:                        //Ag
            alpha = 4.4411759377;
            break;
          case 48:                        //Cd
            alpha = 0.9435470104;
            break;
          case 49:                        //In
            alpha = 2.2334051132;
            break;
          case 50:                        //Sn
            alpha = 1.7715220451;
            break;
          case 51:                        //Sb
            alpha = 1.1172870398;
            break;
          case 52:                        //Te
            alpha = 1.3002599478;
            break;
          case 53:                        //I
            alpha = 1.5741610527;
            break;
          case 54:                        //Xe
            alpha = 1.5463960171;
            break;
          case 55:                        //Cs
            alpha = 1.6736630201;
            break;
          case 56:                        //Ba
            alpha = 1.0710439682;
            break;
          case 57:                        //La
            alpha = 0.9937530160;
            break;
          case 71:                        //Lu
            alpha = 2.7536358833;
            break;
          case 72:                        //Hf
            alpha = 1.2971169949;
            break;
          case 73:                        //Ta
            alpha = 1.6088049412;
            break;
          case 74:                        //W
            alpha = 1.3106900454;
            break;
          case 75:                        //Re
            alpha = 1.4338749647;
            break;
          case 76:                        //Os
            alpha = 1.8330700397;
            break;
          case 77:                        //Ir
            alpha = 2.0097699165;
            break;
          case 78:                        //Pt
            alpha = 2.1085259914;
            break;
          case 79:                        //Au
            alpha = 3.5394139290;
            break;
          case 80:                        //Hg
            alpha = 0.4307309985;
            break;
          case 81:                        //Tl
            alpha = 0.8461930156;
            break;
          case 82:                        //Pb
            alpha = 1.0941230059;
            break;
          case 83:                        //Bi
            alpha = 1.2729749680;
            break;
        }
        break;
      case 18:                    //Ar
        switch (atmB) {
          case 18:                        //Ar
            alpha = 2.3064320087;
            break;
          case 19:                        //K
            alpha = 2.3028030396;
            break;
          case 20:                        //Ca
            alpha = 1.0348809958;
            break;
          case 35:                        //Br
            alpha = 2.4508008957;
            break;
          case 36:                        //Kr
            alpha = 1.9951250553;
            break;
          case 37:                        //Rb
            alpha = 2.5109770298;
            break;
          case 53:                        //I
            alpha = 1.5765869617;
            break;
          case 54:                        //Xe
            alpha = 0.5915200114;
            break;
        }
        break;
      case 19:                    //K
        switch (atmB) {
          case 19:                        //K
            alpha = 1.4355139732;
            break;
          case 20:                        //Ca
            alpha = 1.119180;
            break;
          case 23:                        //V
            alpha = 4.5213599205;
            break;
          case 24:                        //Cr
            alpha = 2.0000000000;
            break;
          case 26:                        //Fe
            alpha = 2.0000000000;
            break;
          case 35:                        //Br
            alpha = 1.6160930395;
            break;
          case 36:                        //Kr
            alpha = 2.1824870110;
            break;
          case 41:                        //Nb
            alpha = 4.5213599205;
            break;
          case 42:                        //Mo
            alpha = 3.9394199848;
            break;
          case 53:                        //I
            alpha = 1.5397139788;
            break;
          case 54:                        //Xe
            alpha = 1.1712499857;
            break;
          case 73:                        //Ta
            alpha = 4.5214700699;
            break;
          case 74:                        //W
            alpha = 4.5213799477;
            break;
          case 76:                        //Os
            alpha = 4.5210900307;
            break;
          case 77:                        //Ir
            alpha = 4.5211701393;
            break;
        }
        break;
      case 20:                    //Ca
        switch (atmB) {
          case 20:                        //Ca
            alpha = 1.8896739483;
            break;
          case 22:                        //Ti
            alpha = 2.0000000000;
            break;
          case 25:                        //Mn
            alpha = 1.4914400578;
            break;
          case 30:                        //Zn
            alpha = 1.1191799641;
            break;
          case 35:                        //Br
            alpha = 2.0784049034;
            break;
          case 36:                        //Kr
            alpha = 1.3051970005;
            break;
          case 53:                        //I
            alpha = 2.1964900494;
            break;
          case 54:                        //Xe
            alpha = 1.5106530190;
            break;
          case 72:                        //Hf
            alpha = 2.0545001030;
            break;
        }
        break;
      case 21:                    //Sc
        switch (atmB) {
          case 21:                        //Sc
            alpha = 1.1328380108;
            break;
          case 35:                        //Br
            alpha = 1.7934859991;
            break;
          case 53:                        //I
            alpha = 1.8148839474;
            break;
        }
        break;
      case 22:                    //Ti
        switch (atmB) {
          case 22:                        //Ti
            alpha = 2.6485970020;
            break;
          case 33:                        //As
            alpha = 1.9329110384;
            break;
          case 35:                        //Br
            alpha = 1.6748470068;
            break;
          case 38:                        //Sr
            alpha = 2.8800299168;
            break;
          case 53:                        //I
            alpha = 1.9334690571;
            break;
          case 56:                        //Ba
            alpha = 2.1760399342;
            break;
          case 80:                        //Hg
            alpha = 3.4146299362;
            break;
        }
        break;
      case 23:                    //V
        switch (atmB) {
          case 23:                        //V
            alpha = 4.8323907852;
            break;
          case 35:                        //Br
            alpha = 1.9029040337;
            break;
          case 53:                        //I
            alpha = 2.6835200787;
            break;
          case 82:                        //Pb
            alpha = 1.5000000000;
            break;
        }
        break;
      case 24:                    //Cr
        switch (atmB) {
          case 24:                        //Cr
            alpha = 4.6554188728;
            break;
          case 35:                        //Br
            alpha = 1.5660279989;
            break;
          case 42:                        //Mo
            alpha = 2.6746160984;
            break;
          case 53:                        //I
            alpha = 2.6342239380;
            break;
          case 82:                        //Pb
            alpha = 1.8607599735;
            break;
        }
        break;
      case 25:                    //Mn
        switch (atmB) {
          case 25:                        //Mn
            alpha = 2.6654200554;
            break;
          case 32:                        //Ge
            alpha = 2.3828339577;
            break;
          case 34:                        //Se
            alpha = 2.6480379105;
            break;
          case 35:                        //Br
            alpha = 2.2838199139;
            break;
          case 51:                        //Sb
            alpha = 2.4003200531;
            break;
          case 53:                        //I
            alpha = 2.2665998936;
            break;
        }
        break;
      case 26:                    //Fe
        switch (atmB) {
          case 26:                        //Fe
            alpha = 2.7207849026;
            break;
          case 35:                        //Br
            alpha = 3.6417820454;
            break;
          case 53:                        //I
            alpha = 1.9128290415;
            break;
        }
        break;
      case 27:                    //Co
        switch (atmB) {
          case 27:                        //Co
            alpha = 3.2881660461;
            break;
          case 32:                        //Ge
            alpha = 2.8526101112;
            break;
          case 33:                        //As
            alpha = 3.3681399822;
            break;
          case 34:                        //Se
            alpha = 2.5234498978;
            break;
          case 35:                        //Br
            alpha = 2.6326880455;
            break;
          case 51:                        //Sb
            alpha = 2.2046298981;
            break;
          case 53:                        //I
            alpha = 3.2352039814;
            break;
        }
        break;
      case 28:                    //Ni
        switch (atmB) {
          case 28:                        //Ni
            alpha = 1.0979599953;
            break;
          case 35:                        //Br
            alpha = 2.7721359730;
            break;
          case 53:                        //I
            alpha = 1.0853430033;
            break;
        }
        break;
      case 29:                    //Cu
        switch (atmB) {
          case 29:                        //Cu
            alpha = 3.6168460846;
            break;
          case 35:                        //Br
            alpha = 5.8264069557;
            break;
          case 53:                        //I
            alpha = 0.8343049884;
            break;
        }
        break;
      case 30:                    //Zn
        switch (atmB) {
          case 30:                        //Zn
            alpha = 0.9290000200;
            break;
          case 33:                        //As
            alpha = 1.4591300488;
            break;
          case 34:                        //Se
            alpha = 1.1862419844;
            break;
          case 35:                        //Br
            alpha = 1.4161200523;
            break;
          case 52:                        //Te
            alpha = 1.2189290524;
            break;
          case 53:                        //I
            alpha = 1.3947620392;
            break;
          case 82:                        //Pb
            alpha = 1.5000000000;
            break;
        }
        break;
      case 31:                    //Ga
        switch (atmB) {
          case 31:                        //Ga
            alpha = 1.3346430063;
            break;
          case 33:                        //As
            alpha = 1.7309770584;
            break;
          case 35:                        //Br
            alpha = 1.8191050291;
            break;
          case 49:                        //In
            alpha = 1.6288700104;
            break;
          case 53:                        //I
            alpha = 1.6717289686;
            break;
        }
        break;
      case 32:                    //Ge
        switch (atmB) {
          case 32:                        //Ge
            alpha = 2.0190000534;
            break;
          case 34:                        //Se
            alpha = 2.6690568924;
            break;
          case 35:                        //Br
            alpha = 1.6023659706;
            break;
          case 43:                        //Tc
            alpha = 2.8528199196;
            break;
          case 44:                        //Ru
            alpha = 2.8523199558;
            break;
          case 45:                        //Rh
            alpha = 2.8529000282;
            break;
          case 50:                        //Sn
            alpha = 2.5246329308;
            break;
          case 52:                        //Te
            alpha = 2.3423719406;
            break;
          case 53:                        //I
            alpha = 1.8174250126;
            break;
          case 75:                        //Re
            alpha = 2.8523399830;
            break;
        }
        break;
      case 33:                    //As
        switch (atmB) {
          case 33:                        //As
            alpha = 1.5882639885;
            break;
          case 34:                        //Se
            alpha = 1.6652799845;
            break;
          case 35:                        //Br
            alpha = 1.5201699734;
            break;
          case 49:                        //In
            alpha = 2.2995519638;
            break;
          case 52:                        //Te
            alpha = 1.1892529726;
            break;
          case 53:                        //I
            alpha = 1.2452620268;
            break;
          case 72:                        //Hf
            alpha = 1.7994999886;
            break;
        }
        break;
      case 34:                    //Se
        switch (atmB) {
          case 34:                        //Se
            alpha = 1.7958940268;
            break;
          case 35:                        //Br
            alpha = 1.4837130308;
            break;
          case 43:                        //Tc
            alpha = 2.5236599445;
            break;
          case 44:                        //Ru
            alpha = 2.5231599808;
            break;
          case 45:                        //Rh
            alpha = 2.5237400532;
            break;
          case 49:                        //In
            alpha = 1.9065719843;
            break;
          case 50:                        //Sn
            alpha = 2.1273770332;
            break;
          case 52:                        //Te
            alpha = 1.5660079718;
            break;
          case 75:                        //Re
            alpha = 2.5231699944;
            break;
          case 82:                        //Pb
            alpha = 2.0000000000;
            break;
          case 83:                        //Bi
            alpha = 1.3447459936;
            break;
        }
        break;
      case 35:                    //Br
        switch (atmB) {
          case 35:                        //Br
            alpha = 1.7581460476;
            break;
          case 36:                        //Kr
            alpha = 1.5290060043;
            break;
          case 37:                        //Rb
            alpha = 1.7977659702;
            break;
          case 38:                        //Sr
            alpha = 1.5243159533;
            break;
          case 39:                        //Y
            alpha = 1.4012080431;
            break;
          case 40:                        //Zr
            alpha = 1.4468679428;
            break;
          case 41:                        //Nb
            alpha = 2.0066781044;
            break;
          case 42:                        //Mo
            alpha = 1.2833340168;
            break;
          case 43:                        //Tc
            alpha = 2.8282639980;
            break;
          case 44:                        //Ru
            alpha = 2.5847349167;
            break;
          case 45:                        //Rh
            alpha = 2.9280819893;
            break;
          case 46:                        //Pd
            alpha = 5.9948787689;
            break;
          case 47:                        //Ag
            alpha = 3.6774909496;
            break;
          case 48:                        //Cd
            alpha = 1.0014510155;
            break;
          case 49:                        //In
            alpha = 2.2579569817;
            break;
          case 50:                        //Sn
            alpha = 1.5350890160;
            break;
          case 51:                        //Sb
            alpha = 1.0639159679;
            break;
          case 52:                        //Te
            alpha = 1.2509399652;
            break;
          case 53:                        //I
            alpha = 1.5793759823;
            break;
          case 54:                        //Xe
            alpha = 1.4396179914;
            break;
          case 55:                        //Cs
            alpha = 1.1671890020;
            break;
          case 56:                        //Ba
            alpha = 1.1903460026;
            break;
          case 57:                        //La
            alpha = 0.7581840158;
            break;
          case 71:                        //Lu
            alpha = 2.3226180077;
            break;
          case 72:                        //Hf
            alpha = 1.0907590389;
            break;
          case 73:                        //Ta
            alpha = 1.6403759718;
            break;
          case 74:                        //W
            alpha = 1.2932599783;
            break;
          case 75:                        //Re
            alpha = 1.6030600071;
            break;
          case 76:                        //Os
            alpha = 1.7668800354;
            break;
          case 77:                        //Ir
            alpha = 2.0381419659;
            break;
          case 78:                        //Pt
            alpha = 2.1853070259;
            break;
          case 79:                        //Au
            alpha = 0.5819110274;
            break;
          case 80:                        //Hg
            alpha = 0.6387169957;
            break;
          case 81:                        //Tl
            alpha = 0.8744189739;
            break;
          case 82:                        //Pb
            alpha = 0.8655499816;
            break;
          case 83:                        //Bi
            alpha = 1.1462329626;
            break;
        }
        break;
      case 36:                    //Kr
        switch (atmB) {
          case 36:                        //Kr
            alpha = 1.1353189945;
            break;
          case 37:                        //Rb
            alpha = 2.2687530518;
            break;
          case 53:                        //I
            alpha = 1.2385740280;
            break;
          case 54:                        //Xe
            alpha = 0.5515609980;
            break;
        }
        break;
      case 37:                    //Rb
        switch (atmB) {
          case 37:                        //Rb
            alpha = 1.1808179617;
            break;
          case 53:                        //I
            alpha = 1.4326750040;
            break;
          case 54:                        //Xe
            alpha = 1.0878230333;
            break;
        }
        break;
      case 38:                    //Sr
        switch (atmB) {
          case 38:                        //Sr
            alpha = 1.0000400543;
            break;
          case 53:                        //I
            alpha = 1.2620420456;
            break;
        }
        break;
      case 39:                    //Y
        switch (atmB) {
          case 39:                        //Y
            alpha = 1.0126810074;
            break;
          case 53:                        //I
            alpha = 1.2791099548;
            break;
        }
        break;
      case 40:                    //Zr
        switch (atmB) {
          case 40:                        //Zr
            alpha = 3.8659679890;
            break;
          case 53:                        //I
            alpha = 1.9951820374;
            break;
        }
        break;
      case 41:                    //Nb
        switch (atmB) {
          case 41:                        //Nb
            alpha = 1.7279410362;
            break;
          case 53:                        //I
            alpha = 1.9672509432;
            break;
          case 82:                        //Pb
            alpha = 1.5000000000;
            break;
        }
        break;
      case 42:                    //Mo
        switch (atmB) {
          case 42:                        //Mo
            alpha = 2.0342540741;
            break;
          case 53:                        //I
            alpha = 0.9484609962;
            break;
          case 82:                        //Pb
            alpha = 2.0000000000;
            break;
        }
        break;
      case 43:                    //Tc
        switch (atmB) {
          case 51:                        //Sb
            alpha = 2.2048499584;
            break;
          case 53:                        //I
            alpha = 1.2923120260;
            break;
        }
        break;
      case 44:                    //Ru
        switch (atmB) {
          case 44:                        //Ru
            alpha = 0.5720559955;
            break;
          case 51:                        //Sb
            alpha = 2.2043499947;
            break;
          case 53:                        //I
            alpha = 3.9532029629;
            break;
        }
        break;
      case 45:                    //Rh
        switch (atmB) {
          case 45:                        //Rh
            alpha = 2.4973280430;
            break;
          case 51:                        //Sb
            alpha = 2.2049300671;
            break;
          case 53:                        //I
            alpha = 3.7081699371;
            break;
        }
        break;
      case 46:                    //Pd
        switch (atmB) {
          case 46:                        //Pd
            alpha = 1.0643750429;
            break;
          case 53:                        //I
            alpha = 5.1445441246;
            break;
        }
        break;
      case 47:                    //Ag
        switch (atmB) {
          case 47:                        //Ag
            alpha = 2.1276450157;
            break;
          case 53:                        //I
            alpha = 2.5931611061;
            break;
        }
        break;
      case 48:                    //Cd
        switch (atmB) {
          case 48:                        //Cd
            alpha = 1.5640439987;
            break;
          case 52:                        //Te
            alpha = 1.3072619438;
            break;
          case 53:                        //I
            alpha = 0.9962379932;
            break;
        }
        break;
      case 49:                    //In
        switch (atmB) {
          case 49:                        //In
            alpha = 2.0732409954;
            break;
          case 51:                        //Sb
            alpha = 2.1419329643;
            break;
          case 52:                        //Te
            alpha = 1.5409879684;
            break;
          case 53:                        //I
            alpha = 2.3517580032;
            break;
        }
        break;
      case 50:                    //Sn
        switch (atmB) {
          case 50:                        //Sn
            alpha = 0.9210000038;
            break;
          case 52:                        //Te
            alpha = 1.7639410496;
            break;
          case 53:                        //I
            alpha = 1.8556330204;
            break;
        }
        break;
      case 51:                    //Sb
        switch (atmB) {
          case 51:                        //Sb
            alpha = 1.3485349417;
            break;
          case 53:                        //I
            alpha = 1.1553150415;
            break;
          case 75:                        //Re
            alpha = 2.2043600082;
            break;
        }
        break;
      case 52:                    //Te
        switch (atmB) {
          case 52:                        //Te
            alpha = 1.1649780273;
            break;
          case 53:                        //I
            alpha = 1.4939509630;
            break;
          case 80:                        //Hg
            alpha = 0.2915000021;
            break;
          case 82:                        //Pb
            alpha = 1.0025589466;
            break;
        }
        break;
      case 53:                    //I
        switch (atmB) {
          case 53:                        //I
            alpha = 1.5199249983;
            break;
          case 54:                        //Xe
            alpha = 0.7991549969;
            break;
          case 55:                        //Cs
            alpha = 0.9195619822;
            break;
          case 56:                        //Ba
            alpha = 0.9825279713;
            break;
          case 57:                        //La
            alpha = 0.5926659703;
            break;
          case 71:                        //Lu
            alpha = 2.2483479977;
            break;
          case 72:                        //Hf
            alpha = 1.0140960217;
            break;
          case 73:                        //Ta
            alpha = 2.4010529518;
            break;
          case 74:                        //W
            alpha = 1.5735700130;
            break;
          case 75:                        //Re
            alpha = 2.6101191044;
            break;
          case 76:                        //Os
            alpha = 2.2037599087;
            break;
          case 77:                        //Ir
            alpha = 3.4109139442;
            break;
          case 78:                        //Pt
            alpha = 3.0773379803;
            break;
          case 79:                        //Au
            alpha = 0.5779160261;
            break;
          case 80:                        //Hg
            alpha = 0.7581620216;
            break;
          case 81:                        //Tl
            alpha = 0.9020119905;
            break;
          case 82:                        //Pb
            alpha = 0.9834740162;
            break;
          case 83:                        //Bi
            alpha = 1.3021709919;
            break;
        }
        break;
      case 54:                    //Xe
        switch (atmB) {
          case 54:                        //Xe
            alpha = 1.2447619438;
            break;
        }
        break;
      case 55:                    //Cs
        switch (atmB) {
          case 55:                        //Cs
            alpha = 1.1708430052;
            break;
        }
        break;
      case 56:                    //Ba
        switch (atmB) {
          case 56:                        //Ba
            alpha = 0.3392690122;
            break;
          case 72:                        //Hf
            alpha = 2.2648301125;
            break;
        }
        break;
      case 57:                    //La
        switch (atmB) {
          case 57:                        //La
            alpha = 4.2480669022;
            break;
        }
        break;
      case 72:                    //Hf
        switch (atmB) {
          case 72:                        //Hf
            alpha = 0.5441439748;
            break;
        }
        break;
      case 73:                    //Ta
        switch (atmB) {
          case 73:                        //Ta
            alpha = 2.0828630924;
            break;
        }
        break;
      case 74:                    //W
        switch (atmB) {
          case 74:                        //W
            alpha = 2.9408700466;
            break;
        }
        break;
      case 75:                    //Re
        switch (atmB) {
          case 75:                        //Re
            alpha = 6.0002579689;
            break;
        }
        break;
      case 76:                    //Os
        switch (atmB) {
          case 76:                        //Os
            alpha = 2.0216300488;
            break;
        }
        break;
      case 77:                    //Ir
        switch (atmB) {
          case 77:                        //Ir
            alpha = 5.7716631889;
            break;
        }
        break;
      case 78:                    //Pt
        switch (atmB) {
          case 78:                        //Pt
            alpha = 3.4042758942;
            break;
        }
        break;
      case 79:                    //Au
        switch (atmB) {
          case 79:                        //Au
            alpha = 0.9031620026;
            break;
        }
        break;
      case 80:                    //Hg
        switch (atmB) {
          case 80:                        //Hg
            alpha = 0.4744130075;
            break;
        }
        break;
      case 81:                    //Tl
        switch (atmB) {
          case 81:                        //Tl
            alpha = 1.1916840076;
            break;
        }
        break;
      case 82:                    //Pb
        switch (atmB) {
          case 82:                        //Pb
            alpha = 1.8817640543;
            break;
        }
        break;
      case 83:                    //Bi
        switch (atmB) {
          case 83:                        //Bi
            alpha = 1.0740640163;
            break;
        }
        break;
    }
    if (alpha == 0.0) {std::cout << "WARNING: MNDOd.hpp: PM6: alpha(): interaction parameter alpha between elements " << atmA << " and " << atmB << " is zero" << std::endl;}
    return alpha*dist_Angstrom2au;
  }
  double xAB(int atm1, int atm2) {
    //function that returns xAB values for PM6; values stored and returned adimensionally
    double xAB = 0.0;
    int atmA = std::min(atm1,atm2);
    int atmB = std::max(atm1,atm2);
    switch (atmA) {
      case 1:                     //H
        switch (atmB) {
          case 1:                         //H
            xAB = 2.2435870171;
            break;
          case 2:                         //He
            xAB = 2.3711988926;
            break;
          case 3:                         //Li
            xAB = 2.1919848919;
            break;
          case 4:                         //Be
            xAB = 2.5628309250;
            break;
          case 5:                         //B
            xAB = 1.3213939667;
            break;
          case 6:                         //C
            xAB = 0.2165060000;
            break;
          case 7:                         //N
            xAB = 0.1755059958;
            break;
          case 8:                         //O
            xAB = 0.1922950000;
            break;
          case 9:                         //F
            xAB = 0.8158019781;
            break;
          case 10:                        //Ne
            xAB = 5.5350208282;
            break;
          case 11:                        //Na
            xAB = 0.2078309953;
            break;
          case 12:                        //Mg
            xAB = 7.7582368851;
            break;
          case 13:                        //Al
            xAB = 2.9583790302;
            break;
          case 14:                        //Si
            xAB = 0.9241960049;
            break;
          case 15:                        //P
            xAB = 1.2349859476;
            break;
          case 16:                        //S
            xAB = 0.8497120142;
            break;
          case 17:                        //Cl
            xAB = 0.7548310161;
            break;
          case 18:                        //Ar
            xAB = 3.9334449768;
            break;
          case 19:                        //K
            xAB = 0.3693400025;
            break;
          case 20:                        //Ca
            xAB = 7.7286062241;
            break;
          case 21:                        //Sc
            xAB = 0.3511990011;
            break;
          case 22:                        //Ti
            xAB = 0.1437219977;
            break;
          case 23:                        //V
            xAB = 0.1052040011;
            break;
          case 24:                        //Cr
            xAB = 0.0444689989;
            break;
          case 25:                        //Mn
            xAB = 1.2692099810;
            break;
          case 26:                        //Fe
            xAB = 0.0251950007;
            break;
          case 27:                        //Co
            xAB = 2.4724650383;
            break;
          case 28:                        //Ni
            xAB = 1.7631239891;
            break;
          case 29:                        //Cu
            xAB = 0.6035910249;
            break;
          case 30:                        //Zn
            xAB = 3.1091930866;
            break;
          case 31:                        //Ga
            xAB = 1.3866519928;
            break;
          case 32:                        //Ge
            xAB = 1.7332259417;
            break;
          case 33:                        //As
            xAB = 1.0905890465;
            break;
          case 34:                        //Se
            xAB = 0.8479980230;
            break;
          case 35:                        //Br
            xAB = 0.8503779769;
            break;
          case 36:                        //Kr
            xAB = 5.1258969307;
            break;
          case 37:                        //Rb
            xAB = 29.861632;
            break;
          case 38:                        //Sr
            xAB = 12.9733161926;
            break;
          case 39:                        //Y
            xAB = 0.6123989820;
            break;
          case 40:                        //Zr
            xAB = 0.5937319994;
            break;
          case 41:                        //Nb
            xAB = 3.6037790775;
            break;
          case 42:                        //Mo
            xAB = 0.9346860051;
            break;
          case 43:                        //Tc
            xAB = 6.3103342056;
            break;
          case 44:                        //Ru
            xAB = 7.1379761696;
            break;
          case 45:                        //Rh
            xAB = 2.3061070442;
            break;
          case 46:                        //Pd
            xAB = 0.4432690144;
            break;
          case 47:                        //Ag
            xAB = 1.9951679707;
            break;
          case 48:                        //Cd
            xAB = 11.9142007828;
            break;
          case 49:                        //In
            xAB = 14.9752931595;
            break;
          case 50:                        //Sn
            xAB = 6.535162;
            break;
          case 51:                        //Sb
            xAB = 0.7953429818;
            break;
          case 52:                        //Te
            xAB = 1.8076790571;
            break;
          case 53:                        //I
            xAB = 0.9818980098;
            break;
          case 54:                        //Xe
            xAB = 0.7010160089;
            break;
          case 55:                        //Cs
            xAB = 0.0969009995;
            break;
          case 56:                        //Ba
            xAB = 2.0400040150;
            break;
          case 57:                        //La
            xAB = 0.6235010028;
            break;
          case 71:                        //Lu
            xAB = 0.7879199982;
            break;
          case 72:                        //Hf
            xAB = 3.4273118973;
            break;
          case 73:                        //Ta
            xAB = 2.8276689053;
            break;
          case 74:                        //W
            xAB = 1.8322700262;
            break;
          case 75:                        //Re
            xAB = 0.3458940089;
            break;
          case 76:                        //Os
            xAB = 4.3938698769;
            break;
          case 77:                        //Ir
            xAB = 0.0580470003;
            break;
          case 78:                        //Pt
            xAB = 8.9240150452;
            break;
          case 79:                        //Au
            xAB = 2.6052830219;
            break;
          case 80:                        //Hg
            xAB = 0.7993990183;
            break;
          case 81:                        //Tl
            xAB = 0.1382050067;
            break;
          case 82:                        //Pb
            xAB = 0.8400959969;
            break;
          case 83:                        //Bi
            xAB = 1.3974620104;
            break;
        }
        break;
      case 2:                //He
        switch (atmB) {
          case 2:                         //He
            xAB = 3.4509000778;
            break;
          case 3:                         //Li
            xAB = 9.2736759186;
            break;
          case 4:                         //Be
            xAB = 12.5448780060;
            break;
          case 5:                         //B
            xAB = 1.9741699696;
            break;
          case 6:                         //C
            xAB = 3.2139708996;
            break;
          case 7:                         //N
            xAB = 1.0778609514;
            break;
          case 8:                         //O
            xAB = 6.6845250130;
            break;
          case 9:                         //F
            xAB = 0.7451069951;
            break;
          case 10:                        //Ne
            xAB = 1.9609240294;
            break;
          case 11:                        //Na
            xAB = 4.2825169563;
            break;
          case 12:                        //Mg
            xAB = 3.7258501053;
            break;
          case 13:                        //Al
            xAB = 2.7014000416;
            break;
          case 14:                        //Si
            xAB = 1.8535829782;
            break;
          case 15:                        //P
            xAB = 1.4902180433;
            break;
          case 16:                        //S
            xAB = 0.4376179874;
            break;
          case 17:                        //Cl
            xAB = 0.2729640007;
            break;
          case 18:                        //Ar
            xAB = 1.1772110462;
            break;
          case 19:                        //K
            xAB = 2.8950450420;
            break;
          case 20:                        //Ca
            xAB = 2.9138519764;
            break;
          case 35:                        //Br
            xAB = 1.0620429516;
            break;
          case 36:                        //Kr
            xAB = 0.6277009845;
            break;
          case 37:                        //Rb
            xAB = 1.8625849485;
            break;
          case 53:                        //I
            xAB = 1.6307209730;
            break;
          case 54:                        //Xe
            xAB = 2.5994710922;
            break;
        }
        break;
      case 3:                     //Li
        switch (atmB) {
          case 3:                         //Li
            xAB = 16.1163845062;
            break;
          case 4:                         //Be
            xAB = 3.2871649265;
            break;
          case 5:                         //B
            xAB = 7.8860177994;
            break;
          case 6:                         //C
            xAB = 16.1800022125;
            break;
          case 7:                         //N
            xAB = 2.8234028816;
            break;
          case 8:                         //O
            xAB = 1.9685980082;
            break;
          case 9:                         //F
            xAB = 1.9759850502;
            break;
          case 10:                        //Ne
            xAB = 0.7049580216;
            break;
          case 11:                        //Na
            xAB = 0.8814820051;
            break;
          case 12:                        //Mg
            xAB = 2.4902501106;
            break;
          case 13:                        //Al
            xAB = 1.1068190336;
            break;
          case 14:                        //Si
            xAB = 3.0907909870;
            break;
          case 15:                        //P
            xAB = 1.1229499578;
            break;
          case 16:                        //S
            xAB = 2.6425020695;
            break;
          case 17:                        //Cl
            xAB = 4.2277941704;
            break;
          case 18:                        //Ar
            xAB = 3.3629100323;
            break;
          case 19:                        //K
            xAB = 4.3745670319;
            break;
          case 35:                        //Br
            xAB = 1.8588659763;
            break;
          case 36:                        //Kr
            xAB = 8.7586965561;
            break;
          case 53:                        //I
            xAB = 4.1685991287;
            break;
          case 54:                        //Xe
            xAB = 4.5820808411;
            break;
          case 82:                        //Pb
            xAB = 1.2850639820;
            break;
          case 83:                        //Bi
            xAB = 0.6953200102;
            break;
        }
        break;
      case 4:                     //Be
        switch (atmB) {
          case 4:                         //Be
            xAB = 0.2386330068;
            break;
          case 5:                         //B
            xAB = 1.1517920494;
            break;
          case 6:                         //C
            xAB = 25.0358791351;
            break;
          case 7:                         //N
            xAB = 1.7406049967;
            break;
          case 8:                         //O
            xAB = 3.2181549072;
            break;
          case 9:                         //F
            xAB = 3.8829929829;
            break;
          case 10:                        //Ne
            xAB = 0.3926280141;
            break;
          case 11:                        //Na
            xAB = 3.1216199398;
            break;
          case 12:                        //Mg
            xAB = 2.0663919449;
            break;
          case 13:                        //Al
            xAB = 5.0372138023;
            break;
          case 14:                        //Si
            xAB = 0.6234329939;
            break;
          case 15:                        //P
            xAB = 1.6848310232;
            break;
          case 16:                        //S
            xAB = 3.7915649414;
            break;
          case 17:                        //Cl
            xAB = 2.5072751045;
            break;
          case 18:                        //Ar
            xAB = 2.7554919720;
            break;
          case 19:                        //K
            xAB = 6.7322211266;
            break;
          case 35:                        //Br
            xAB = 1.9409329891;
            break;
          case 36:                        //Kr
            xAB = 10.2377958298;
            break;
          case 53:                        //I
            xAB = 2.3518979549;
            break;
          case 54:                        //Xe
            xAB = 0.6605250239;
            break;
        }
        break;
      case 5:                     //B
        switch (atmB) {
          case 5:                         //B
            xAB = 3.5936191082;
            break;
          case 6:                         //C
            xAB = 1.8748589754;
            break;
          case 7:                         //N
            xAB = 0.9528819919;
            break;
          case 8:                         //O
            xAB = 1.2698010206;
            break;
          case 9:                         //F
            xAB = 0.8627610207;
            break;
          case 10:                        //Ne
            xAB = 2.7641398907;
            break;
          case 11:                        //Na
            xAB = 3.1886079311;
            break;
          case 12:                        //Mg
            xAB = 6.1467008591;
            break;
          case 13:                        //Al
            xAB = 2.7414789200;
            break;
          case 14:                        //Si
            xAB = 1.0282870531;
            break;
          case 15:                        //P
            xAB = 1.4508860111;
            break;
          case 16:                        //S
            xAB = 1.1253939867;
            break;
          case 17:                        //Cl
            xAB = 0.8221290112;
            break;
          case 18:                        //Ar
            xAB = 1.9315860271;
            break;
          case 19:                        //K
            xAB = 8.9005413055;
            break;
          case 20:                        //Ca
            xAB = 1.7000099421;
            break;
          case 22:                        //Ti
            xAB = 0.6493600011;
            break;
          case 35:                        //Br
            xAB = 1.2264200449;
            break;
          case 36:                        //Kr
            xAB = 2.9467809200;
            break;
          case 37:                        //Rb
            xAB = 9.0404930115;
            break;
          case 53:                        //I
            xAB = 3.1613850594;
            break;
          case 54:                        //Xe
            xAB = 1.1006120443;
            break;
          case 55:                        //Cs
            xAB = 10.3926095963;
            break;
          case 72:                        //Hf
            xAB = 0.6592699885;
            break;
          case 81:                        //Tl
            xAB = 10.5043382645;
            break;
          case 82:                        //Pb
            xAB = 1.1381570101;
            break;
        }
        break;
      case 6:                     //C
        switch (atmB) {
          case 6:                         //C
            xAB = 0.8135100007;
            break;
          case 7:                         //N
            xAB = 0.8599489927;
            break;
          case 8:                         //O
            xAB = 0.9902110100;
            break;
          case 9:                         //F
            xAB = 0.7329679728;
            break;
          case 10:                        //Ne
            xAB = 5.46878004070;
            break;
          case 11:                        //Na
            xAB = 4.5204291344;
            break;
          case 12:                        //Mg
            xAB = 10.5176897049;
            break;
          case 13:                        //Al
            xAB = 2.9280560017;
            break;
          case 14:                        //Si
            xAB = 0.7857450247;
            break;
          case 15:                        //P
            xAB = 0.9795119762;
            break;
          case 16:                        //S
            xAB = 0.6668490171;
            break;
          case 17:                        //Cl
            xAB = 0.5157870054;
            break;
          case 18:                        //Ar
            xAB = 0.1223089993;
            break;
          case 19:                        //K
            xAB = 10.3179874420;
            break;
          case 20:                        //Ca
            xAB = 0.1484500021;
            break;
          case 21:                        //Sc
            xAB = 8.6080522537;
            break;
          case 22:                        //Ti
            xAB = 0.4167059958;
            break;
          case 23:                        //V
            xAB = 1.9387600422;
            break;
          case 24:                        //Cr
            xAB = 6.1101870537;
            break;
          case 25:                        //Mn
            xAB = 2.5831100941;
            break;
          case 26:                        //Fe
            xAB = 0.3668349981;
            break;
          case 27:                        //Co
            xAB = 2.1239299774;
            break;
          case 28:                        //Ni
            xAB = 7.1333241463;
            break;
          case 29:                        //Cu
            xAB = 7.0677938461;
            break;
          case 30:                        //Zn
            xAB = 0.9914649725;
            break;
          case 31:                        //Ga
            xAB = 1.9629900455;
            break;
          case 32:                        //Ge
            xAB = 1.2975100279;
            break;
          case 33:                        //As
            xAB = 0.5790979862;
            break;
          case 34:                        //Se
            xAB = 1.1147869825;
            break;
          case 35:                        //Br
            xAB = 0.5706859827;
            break;
          case 36:                        //Kr
            xAB = 0.6526229978;
            break;
          case 37:                        //Rb
            xAB = 29.9740314484;
            break;
          case 38:                        //Sr
            xAB = 6.6546568871;
            break;
          case 39:                        //Y
            xAB = 0.5043060184;
            break;
          case 40:                        //Zr
            xAB = 1.9991819859;
            break;
          case 41:                        //Nb
            xAB = 4.5754809380;
            break;
          case 42:                        //Mo
            xAB = 1.1907420158;
            break;
          case 43:                        //Tc
            xAB = 3.9724390507;
            break;
          case 44:                        //Ru
            xAB = 1.1349359751;
            break;
          case 45:                        //Rh
            xAB = 3.4880790710;
            break;
          case 46:                        //Pd
            xAB = 9.8537149429;
            break;
          case 47:                        //Ag
            xAB = 11.3354558945;
            break;
          case 48:                        //Cd
            xAB = 0.6034410000;
            break;
          case 49:                        //In
            xAB = 2.1873850822;
            break;
          case 50:                        //Sn
            xAB = 3.3743550777;
            break;
          case 51:                        //Sb
            xAB = 0.5792120099;
            break;
          case 52:                        //Te
            xAB = 0.9704939723;
            break;
          case 53:                        //I
            xAB = 0.8101559877;
            break;
          case 54:                        //Xe
            xAB = 0.8267269731;
            break;
          case 55:                        //Cs
            xAB = 24.5146236420;
            break;
          case 56:                        //Ba
            xAB = 0.1197929978;
            break;
          case 57:                        //La
            xAB = 0.1086490005;
            break;
          case 71:                        //Lu
            xAB = 4.4538249969;
            break;
          case 72:                        //Hf
            xAB = 0.3785789907;
            break;
          case 73:                        //Ta
            xAB = 0.8474389911;
            break;
          case 74:                        //W
            xAB = 1.1607700586;
            break;
          case 75:                        //Re
            xAB = 0.6906870008;
            break;
          case 76:                        //Os
            xAB = 0.4984099865;
            break;
          case 77:                        //Ir
            xAB = 0.1150470003;
            break;
          case 78:                        //Pt
            xAB = 3.4934029579;
            break;
          case 79:                        //Au
            xAB = 21.4856338501;
            break;
          case 80:                        //Hg
            xAB = 0.1471280009;
            break;
          case 81:                        //Tl
            xAB = 0.5828949809;
            break;
          case 82:                        //Pb
            xAB = 0.4046559930;
            break;
          case 83:                        //Bi
            xAB = 0.5761790276;
            break;
        }
        break;
      case 7:                     //N
        switch (atmB) {
          case 7:                         //N
            xAB = 0.6753129959;
            break;
          case 8:                         //O
            xAB = 0.7647560239;
            break;
          case 9:                         //F
            xAB = 0.6358540058;
            break;
          case 10:                        //Ne
            xAB = 29.9996089935;
            break;
          case 11:                        //Na
            xAB = 8.5863866806;
            break;
          case 12:                        //Mg
            xAB = 1.20807504655;
            break;
          case 13:                        //Al
            xAB = 1.3452019691;
            break;
          case 14:                        //Si
            xAB = 0.5929719806;
            break;
          case 15:                        //P
            xAB = 0.9721540213;
            break;
          case 16:                        //S
            xAB = 0.7387099862;
            break;
          case 17:                        //Cl
            xAB = 0.5207449794;
            break;
          case 18:                        //Ar
            xAB = 0.5625810027;
            break;
          case 19:                        //K
            xAB = 7.9777069092;
            break;
          case 20:                        //Ca
            xAB = 2.9880740643;
            break;
          case 21:                        //Sc
            xAB = 3.2318809032;
            break;
          case 22:                        //Ti
            xAB = 0.5454609990;
            break;
          case 23:                        //V
            xAB = 0.2767249942;
            break;
          case 24:                        //Cr
            xAB = 1.9203239679;
            break;
          case 25:                        //Mn
            xAB = 1.9567500353;
            break;
          case 26:                        //Fe
            xAB = 0.1553419977;
            break;
          case 27:                        //Co
            xAB = 2.6538360119;
            break;
          case 28:                        //Ni
            xAB = 4.2867999077;
            break;
          case 29:                        //Cu
            xAB = 3.2286670208;
            break;
          case 30:                        //Zn
            xAB = 0.9524760246;
            break;
          case 31:                        //Ga
            xAB = 1.1883380413;
            break;
          case 32:                        //Ge
            xAB = 0.6375060081;
            break;
          case 33:                        //As
            xAB = 0.2733370066;
            break;
          case 34:                        //Se
            xAB = 0.4828400016;
            break;
          case 35:                        //Br
            xAB = 30.0001335144;
            break;
          case 36:                        //Kr
            xAB = 0.1996060014;
            break;
          case 37:                        //Rb
            xAB = 0.0246360004;
            break;
          case 38:                        //Sr
            xAB = 6.8538661003;
            break;
          case 39:                        //Y
            xAB = 1.6279029846;
            break;
          case 40:                        //Zr
            xAB = 0.9950450063;
            break;
          case 41:                        //Nb
            xAB = 1.2135870457;
            break;
          case 42:                        //Mo
            xAB = 0.6082680225;
            break;
          case 43:                        //Tc
            xAB = 0.7271299958;
            break;
          case 44:                        //Ru
            xAB = 2.3340940475;
            break;
          case 45:                        //Rh
            xAB = 4.0009469986;
            break;
          case 46:                        //Pd
            xAB = 0.2497030050;
            break;
          case 47:                        //Ag
            xAB = 19.8037109375;
            break;
          case 48:                        //Cd
            xAB = 0.1806630045;
            break;
          case 49:                        //In
            xAB = 3.3699929714;
            break;
          case 50:                        //Sn
            xAB = 1.3918999434;
            break;
          case 51:                        //Sb
            xAB = 0.0820650011;
            break;
          case 52:                        //Te
            xAB = 0.3585929871;
            break;
          case 53:                        //I
            xAB = 0.2649030089;
            break;
          case 54:                        //Xe
            xAB = 0.9256240129;
            break;
          case 55:                        //Cs
            xAB = 29.7110767365;
            break;
          case 56:                        //Ba
            xAB = 0.2079340070;
            break;
          case 57:                        //La
            xAB = 0.1047779992;
            break;
          case 71:                        //Lu
            xAB = 2.8608279228;
            break;
          case 72:                        //Hf
            xAB = 0.6557949781;
            break;
          case 73:                        //Ta
            xAB = 1.0154609680;
            break;
          case 74:                        //W
            xAB = 0.4783500135;
            break;
          case 75:                        //Re
            xAB = 0.4452129900;
            break;
          case 76:                        //Os
            xAB = 0.0808700025;
            break;
          case 77:                        //Ir
            xAB = 8.5186395645;
            break;
          case 78:                        //Pt
            xAB = 0.5407299995;
            break;
          case 79:                        //Au
            xAB = 0.2220589966;
            break;
          case 80:                        //Hg
            xAB = 0.0502399988;
            break;
          case 81:                        //Tl
            xAB = 0.1588120013;
            break;
          case 82:                        //Pb
            xAB = 0.3357869983;
            break;
          case 83:                        //Bi
            xAB = 0.1527380049;
            break;
        }
        break;
      case 8:                     //O
        switch (atmB) {
          case 8:                         //O
            xAB = 0.5351120234;
            break;
          case 9:                         //F
            xAB = 0.6742510200;
            break;
          case 10:                        //Ne
            xAB = 0.7638990283;
            break;
          case 11:                        //Na
            xAB = 3.2700788975;
            break;
          case 12:                        //Mg
            xAB = 1.5357340574;
            break;
          case 13:                        //Al
            xAB = 2.1313960552;
            break;
          case 14:                        //Si
            xAB = 0.7510949969;
            break;
          case 15:                        //P
            xAB = 0.8787050247;
            break;
          case 16:                        //S
            xAB = 0.7472149730;
            break;
          case 17:                        //Cl
            xAB = 0.5855100155;
            break;
          case 18:                        //Ar
            xAB = 0.3557949960;
            break;
          case 19:                        //K
            xAB = 1.0122749805;
            break;
          case 20:                        //Ca
            xAB = 17.0289459229;
            break;
          case 21:                        //Sc
            xAB = 3.0586719513;
            break;
          case 22:                        //Ti
            xAB = 0.7994859815;
            break;
          case 23:                        //V
            xAB = 0.4153119922;
            break;
          case 24:                        //Cr
            xAB = 1.0555109978;
            break;
          case 25:                        //Mn
            xAB = 1.2856199741;
            break;
          case 26:                        //Fe
            xAB = 0.1364219934;
            break;
          case 27:                        //Co
            xAB = 5.2520217896;
            break;
          case 28:                        //Ni
            xAB = 0.8859419823;
            break;
          case 29:                        //Cu
            xAB = 2.0000000000;
            break;
          case 30:                        //Zn
            xAB = 2.2653129101;
            break;
          case 31:                        //Ga
            xAB = 1.5236439705;
            break;
          case 32:                        //Ge
            xAB = 0.8269640207;
            break;
          case 33:                        //As
            xAB = 0.7016140223;
            break;
          case 34:                        //Se
            xAB = 0.9551609755;
            break;
          case 35:                        //Br
            xAB = 0.7065839767;
            break;
          case 36:                        //Kr
            xAB = 0.0067329998;
            break;
          case 37:                        //Rb
            xAB = 1.1253499985;
            break;
          case 38:                        //Sr
            xAB = 3.5613958836;
            break;
          case 39:                        //Y
            xAB = 1.7429219484;
            break;
          case 40:                        //Zr
            xAB = 1.0575250387;
            break;
          case 41:                        //Nb
            xAB = 1.1847189665;
            break;
          case 42:                        //Mo
            xAB = 0.5112670064;
            break;
          case 43:                        //Tc
            xAB = 1.0246160030;
            break;
          case 44:                        //Ru
            xAB = 2.9762790203;
            break;
          case 45:                        //Rh
            xAB = 10.2986755371;
            break;
          case 46:                        //Pd
            xAB = 0.2164029926;
            break;
          case 47:                        //Ag
            xAB = 0.1656610072;
            break;
          case 48:                        //Cd
            xAB = 0.9261459708;
            break;
          case 49:                        //In
            xAB = 4.1285829544;
            break;
          case 50:                        //Sn
            xAB = 4.3740167618;
            break;
          case 51:                        //Sb
            xAB = 0.6342340112;
            break;
          case 52:                        //Te
            xAB = 0.3829259872;
            break;
          case 53:                        //I
            xAB = 0.8662040234;
            break;
          case 54:                        //Xe
            xAB = 0.0353560001;
            break;
          case 55:                        //Cs
            xAB = 8.1768426895;
            break;
          case 56:                        //Ba
            xAB = 0.3489449918;
            break;
          case 57:                        //La
            xAB = 0.5579569936;
            break;
          case 71:                        //Lu
            xAB = 2.9170761108;
            break;
          case 72:                        //Hf
            xAB = 1.0347180367;
            break;
          case 73:                        //Ta
            xAB = 1.7510830164;
            break;
          case 74:                        //W
            xAB = 0.3490099907;
            break;
          case 75:                        //Re
            xAB = 0.6359599829;
            break;
          case 76:                        //Os
            xAB = 0.1843000054;
            break;
          case 77:                        //Ir
            xAB = 9.6254024506;
            break;
          case 78:                        //Pt
            xAB = 0.4877560139;
            break;
          case 79:                        //Au
            xAB = 0.0771920010;
            break;
          case 80:                        //Hg
            xAB = 0.2407200038;
            break;
          case 81:                        //Tl
            xAB = 0.6369060278;
            break;
          case 82:                        //Pb
            xAB = 0.7825059891;
            break;
          case 83:                        //Bi
            xAB = 0.3330419958;
            break;
        }
        break;
      case 9:                     //F
        switch (atmB) {
          case 9:                         //F
            xAB = 0.6813430190;
            break;
          case 10:                        //Ne
            xAB = 2.7067539692;
            break;
          case 11:                        //Na
            xAB = 7.0473508835;
            break;
          case 12:                        //Mg
            xAB = 5.8590230942;
            break;
          case 13:                        //Al
            xAB = 1.9756350517;
            break;
          case 14:                        //Si
            xAB = 0.5435159802;
            break;
          case 15:                        //P
            xAB = 0.5145750046;
            break;
          case 16:                        //S
            xAB = 0.3752509952;
            break;
          case 17:                        //Cl
            xAB = 0.4111239910;
            break;
          case 18:                        //Ar
            xAB = 9.2697153091;
            break;
          case 19:                        //K
            xAB = 6.5929708481;
            break;
          case 20:                        //Ca
            xAB = 3.4828209877;
            break;
          case 21:                        //Sc
            xAB = 7.2523469925;
            break;
          case 22:                        //Ti
            xAB = 1.0857419968;
            break;
          case 23:                        //V
            xAB = 0.3428150117;
            break;
          case 24:                        //Cr
            xAB = 0.7376070023;
            break;
          case 25:                        //Mn
            xAB = 1.1130700111;
            break;
          case 26:                        //Fe
            xAB = 3.6573500633;
            break;
          case 27:                        //Co
            xAB = 4.5850300789;
            break;
          case 28:                        //Ni
            xAB = 1.0882079601;
            break;
          case 29:                        //Cu
            xAB = 1.2821079493;
            break;
          case 30:                        //Zn
            xAB = 1.2255450487;
            break;
          case 31:                        //Ga
            xAB = 1.4169420004;
            break;
          case 32:                        //Ge
            xAB = 0.6519770026;
            break;
          case 33:                        //As
            xAB = 0.4026280046;
            break;
          case 34:                        //Se
            xAB = 0.4448060095;
            break;
          case 35:                        //Br
            xAB = 0.2935000062;
            break;
          case 36:                        //Kr
            xAB = 4.1344070435;
            break;
          case 37:                        //Rb
            xAB = 28.8152770996;
            break;
          case 38:                        //Sr
            xAB = 10.9717054367;
            break;
          case 39:                        //Y
            xAB = 4.4338088036;
            break;
          case 40:                        //Zr
            xAB = 0.8611419797;
            break;
          case 41:                        //Nb
            xAB = 3.6636819839;
            break;
          case 42:                        //Mo
            xAB = 0.6104289889;
            break;
          case 43:                        //Tc
            xAB = 5.8117837906;
            break;
          case 44:                        //Ru
            xAB = 6.9471278191;
            break;
          case 45:                        //Rh
            xAB = 9.0653839111;
            break;
          case 46:                        //Pd
            xAB = 6.9453120232;
            break;
          case 47:                        //Ag
            xAB = 12.6958837509;
            break;
          case 48:                        //Cd
            xAB = 1.3536649942;
            break;
          case 49:                        //In
            xAB = 3.7010159492;
            break;
          case 50:                        //Sn
            xAB = 18.5986633301;
            break;
          case 51:                        //Sb
            xAB = 0.6502770185;
            break;
          case 52:                        //Te
            xAB = 0.2008219957;
            break;
          case 53:                        //I
            xAB = 0.3924250007;
            break;
          case 54:                        //Xe
            xAB = 0.0650110021;
            break;
          case 55:                        //Cs
            xAB = 22.2334156036;
            break;
          case 56:                        //Ba
            xAB = 5.5752549171;
            break;
          case 57:                        //La
            xAB = 2.4019029140;
            break;
          case 72:                        //Hf
            xAB = 1.6793349981;
            break;
          case 73:                        //Ta
            xAB = 3.1465198994;
            break;
          case 74:                        //W
            xAB = 0.2138900012;
            break;
          case 75:                        //Re
            xAB = 0.5359659791;
            break;
          case 76:                        //Os
            xAB = 0.1400499940;
            break;
          case 77:                        //Ir
            xAB = 1.4996390343;
            break;
          case 78:                        //Pt
            xAB = 5.6170139313;
            break;
          case 79:                        //Au
            xAB = 9.5943841934;
            break;
          case 80:                        //Hg
            xAB = 0.1138269976;
            break;
          case 81:                        //Tl
            xAB = 0.2261659950;
            break;
          case 82:                        //Pb
            xAB = 8.3685617447;
            break;
          case 83:                        //Bi
            xAB = 1.0353239775;
            break;
        }
        break;
      case 10:                    //Ne
        switch (atmB) {
          case 10:                        //Ne
            xAB = 2.7948300838;
            break;
          case 11:                        //Na
            xAB = 1.3430370092;
            break;
          case 12:                        //Mg
            xAB = 1.2148590088;
            break;
          case 13:                        //Al
            xAB = 1.7092000246;
            break;
          case 14:                        //Si
            xAB = 14.3786764145;
            break;
          case 15:                        //P
            xAB = 0.7749540210;
            break;
          case 16:                        //S
            xAB = 3.2961599827;
            break;
          case 17:                        //Cl
            xAB = 0.1251329929;
            break;
          case 18:                        //Ar
            xAB = 1.3046970367;
            break;
          case 19:                        //K
            xAB = 0.2339950055;
            break;
          case 20:                        //Ca
            xAB = 0.3325859904;
            break;
          case 35:                        //Br
            xAB = 1.0061589479;
            break;
          case 36:                        //Kr
            xAB = 1.4337220192;
            break;
          case 37:                        //Rb
            xAB = 7.7365632057;
            break;
          case 53:                        //I
            xAB = 1.5035680532;
            break;
          case 54:                        //Xe
            xAB = 0.2938619852;
            break;
        }
        break;
      case 11:                    //Na
        switch (atmB) {
          case 11:                        //Na
            xAB = 0.2871370018;
            break;
          case 12:                        //Mg
            xAB = 8.6756191254;
            break;
          case 13:                        //Al
            xAB = 2.0718469620;
            break;
          case 14:                        //Si
            xAB = 9.2376441956;
            break;
          case 15:                        //P
            xAB = 2.8370950222;
            break;
          case 16:                        //S
            xAB = 0.8524339795;
            break;
          case 17:                        //Cl
            xAB = 1.3578939438;
            break;
          case 18:                        //Ar
            xAB = 3.3981380463;
            break;
          case 19:                        //K
            xAB = 5.5630269051;
            break;
          case 20:                        //Ca
            xAB = 9.6575088501;
            break;
          case 23:                        //V
            xAB = 8.2760200500;
            break;
          case 24:                        //Cr
            xAB = 8.3642740250;
            break;
          case 35:                        //Br
            xAB = 1.7529369593;
            break;
          case 36:                        //Kr
            xAB = 8.3544483185;
            break;
          case 41:                        //Nb
            xAB = 8.2760200500;
            break;
          case 42:                        //Mo
            xAB = 8.2865495682;
            break;
          case 53:                        //I
            xAB = 1.9861119986;
            break;
          case 54:                        //Xe
            xAB = 8.3682041168;
            break;
          case 73:                        //Ta
            xAB = 8.2761297226;
            break;
          case 74:                        //W
            xAB = 8.2760400772;
            break;
          case 76:                        //Os
            xAB = 8.2757501602;
            break;
          case 77:                        //Ir
            xAB = 8.2758302689;
            break;
        }
        break;
      case 12:                    //Mg
        switch (atmB) {
          case 12:                        //Mg
            xAB = 0.4656449854;
            break;
          case 13:                        //Al
            xAB = 13.4724426270;
            break;
          case 14:                        //Si
            xAB = 29.9945201874;
            break;
          case 15:                        //P
            xAB = 1.1778810024;
            break;
          case 16:                        //S
            xAB = 0.5007479787;
            break;
          case 17:                        //Cl
            xAB = 2.4308559895;
            break;
          case 18:                        //Ar
            xAB = 1.9706380367;
            break;
          case 19:                        //K
            xAB = 3.2905020714;
            break;
          case 20:                        //Ca
            xAB = 8.5998001099;
            break;
          case 22:                        //Ti
            xAB = 4.3302397728;
            break;
          case 35:                        //Br
            xAB = 2.9162800312;
            break;
          case 36:                        //Kr
            xAB = 0.8884360194;
            break;
          case 53:                        //I
            xAB = 3.2769138813;
            break;
          case 54:                        //Xe
            xAB = 9.7235717773;
            break;
          case 72:                        //Hf
            xAB = 4.3302497864;
            break;
        }
        break;
      case 13:                    //Al
        switch (atmB) {
          case 13:                        //Al
            xAB = 2.1391999722;
            break;
          case 14:                        //Si
            xAB = 2.0000000000;
            break;
          case 15:                        //P
            xAB = 5.0508160591;
            break;
          case 16:                        //S
            xAB = 2.3473839760;
            break;
          case 17:                        //Cl
            xAB = 2.1534509659;
            break;
          case 18:                        //Ar
            xAB = 1.8520090580;
            break;
          case 19:                        //K
            xAB = 29.9447078705;
            break;
          case 20:                        //Ca
            xAB = 4.1885547638;
            break;
          case 21:                        //Sc
            xAB = 0.5006200075;
            break;
          case 22:                        //Ti
            xAB = 2.0918409824;
            break;
          case 25:                        //Mn
            xAB = 1.0407899618;
            break;
          case 33:                        //As
            xAB = 1.0035799742;
            break;
          case 35:                        //Br
            xAB = 2.3571300507;
            break;
          case 36:                        //Kr
            xAB = 5.0917158127;
            break;
          case 37:                        //Rb
            xAB = 2.9924569130;
            break;
          case 39:                        //Y
            xAB = 0.5006700158;
            break;
          case 40:                        //Zr
            xAB = 0.8740599751;
            break;
          case 46:                        //Pd
            xAB = 1.0572899580;
            break;
          case 47:                        //Ag
            xAB = 0.8965139985;
            break;
          case 51:                        //Sb
            xAB = 1.6166900396;
            break;
          case 52:                        //Te
            xAB = 2.1068120003;
            break;
          case 53:                        //I
            xAB = 2.9290800095;
            break;
          case 54:                        //Xe
            xAB = 7.4044651985;
            break;
          case 56:                        //Ba
            xAB = 9.5390987396;
            break;
          case 57:                        //La
            xAB = 0.5005400181;
            break;
          case 72:                        //Hf
            xAB = 0.6225200295;
            break;
          case 78:                        //Pt
            xAB = 1.0569299459;
            break;
          case 79:                        //Au
            xAB = 1.0571399927;
            break;
        }
        break;
      case 14:                    //Si
        switch (atmB) {
          case 14:                        //Si
            xAB = 0.2734769881;
            break;
          case 15:                        //P
            xAB = 13.2391214371;
            break;
          case 16:                        //S
            xAB = 0.8766580224;
            break;
          case 17:                        //Cl
            xAB = 0.5130000114;
            break;
          case 18:                        //Ar
            xAB = 1.0678900480;
            break;
          case 19:                        //K
            xAB = 8.2791996002;
            break;
          case 20:                        //Ca
            xAB = 0.3362329900;
            break;
          case 21:                        //Sc
            xAB = 3.2190699577;
            break;
          case 22:                        //Ti
            xAB = 6.7738151550;
            break;
          case 24:                        //Cr
            xAB = 1.0291099548;
            break;
          case 25:                        //Mn
            xAB = 0.9505800009;
            break;
          case 27:                        //Co
            xAB = 1.0902400017;
            break;
          case 28:                        //Ni
            xAB = 0.9386460185;
            break;
          case 30:                        //Zn
            xAB = 3.7839050293;
            break;
          case 31:                        //Ga
            xAB = 1.0022900105;
            break;
          case 32:                        //Ge
            xAB = 0.1786800027;
            break;
          case 33:                        //As
            xAB = 1.4307060242;
            break;
          case 34:                        //Se
            xAB = 0.5182269812;
            break;
          case 35:                        //Br
            xAB = 0.5895109773;
            break;
          case 36:                        //Kr
            xAB = 0.5542500019;
            break;
          case 38:                        //Sr
            xAB = 2.7647500038;
            break;
          case 39:                        //Y
            xAB = 3.2190299034;
            break;
          case 40:                        //Zr
            xAB = 1.7233430147;
            break;
          case 44:                        //Ru
            xAB = 0.8494300246;
            break;
          case 45:                        //Rh
            xAB = 0.8500099778;
            break;
          case 46:                        //Pd
            xAB = 2.2251040936;
            break;
          case 48:                        //Cd
            xAB = 2.2533459663;
            break;
          case 51:                        //Sb
            xAB = 8.7137489319;
            break;
          case 53:                        //I
            xAB = 0.7002990246;
            break;
          case 54:                        //Xe
            xAB = 16.0919990540;
            break;
          case 56:                        //Ba
            xAB = 1.2126599550;
            break;
          case 57:                        //La
            xAB = 3.2190299034;
            break;
          case 72:                        //Hf
            xAB = 3.3822999001;
            break;
          case 75:                        //Re
            xAB = 0.8494499922;
            break;
          case 78:                        //Pt
            xAB = 0.0999900028;
            break;
          case 80:                        //Hg
            xAB = 3.6807401180;
            break;
        }
        break;
      case 15:                    //P
        switch (atmB) {
          case 15:                        //P
            xAB = 0.9025009871;
            break;
          case 16:                        //S
            xAB = 0.5622659922;
            break;
          case 17:                        //Cl
            xAB = 0.3523609936;
            break;
          case 18:                        //Ar
            xAB = 0.1710139960;
            break;
          case 19:                        //K
            xAB = 4.4753842354;
            break;
          case 20:                        //Ca
            xAB = 0.4108400047;
            break;
          case 21:                        //Sc
            xAB = 0.6267489791;
            break;
          case 22:                        //Ti
            xAB = 4.1504998207;
            break;
          case 23:                        //V
            xAB = 6.2506241798;
            break;
          case 24:                        //Cr
            xAB = 0.6001769900;
            break;
          case 25:                        //Mn
            xAB = 1.1303199530;
            break;
          case 26:                        //Fe
            xAB = 0.4312910140;
            break;
          case 27:                        //Co
            xAB = 0.1059359983;
            break;
          case 28:                        //Ni
            xAB = 5.0597271919;
            break;
          case 29:                        //Cu
            xAB = 0.0201260000;
            break;
          case 30:                        //Zn
            xAB = 0.5815299749;
            break;
          case 31:                        //Ga
            xAB = 0.5000000000;
            break;
          case 32:                        //Ge
            xAB = 5.6163492203;
            break;
          case 34:                        //Se
            xAB = 0.2920520008;
            break;
          case 35:                        //Br
            xAB = 0.4565210044;
            break;
          case 38:                        //Sr
            xAB = 2.5520999432;
            break;
          case 39:                        //Y
            xAB = 0.5416600108;
            break;
          case 40:                        //Zr
            xAB = 0.7483760118;
            break;
          case 41:                        //Nb
            xAB = 6.2015070915;
            break;
          case 42:                        //Mo
            xAB = 1.5228459835;
            break;
          case 44:                        //Ru
            xAB = 0.0569740012;
            break;
          case 45:                        //Rh
            xAB = 1.0381410122;
            break;
          case 46:                        //Pd
            xAB = 0.0450170003;
            break;
          case 47:                        //Ag
            xAB = 0.0499319993;
            break;
          case 52:                        //Te
            xAB = 1.1092890501;
            break;
          case 53:                        //I
            xAB = 3.0472071171;
            break;
          case 55:                        //Cs
            xAB = 0.5065119863;
            break;
          case 57:                        //La
            xAB = 0.5416600108;
            break;
          case 71:                        //Lu
            xAB = 0.5000000000;
            break;
          case 72:                        //Hf
            xAB = 0.5055299997;
            break;
          case 73:                        //Ta
            xAB = 6.2618799210;
            break;
          case 74:                        //W
            xAB = 5.9538598061;
            break;
          case 75:                        //Re
            xAB = 0.9669420123;
            break;
          case 76:                        //Os
            xAB = 6.0583000183;
            break;
          case 77:                        //Ir
            xAB = 6.2846698761;
            break;
          case 78:                        //Pt
            xAB = 0.2337120026;
            break;
          case 79:                        //Au
            xAB = 0.0670010000;
            break;
          case 80:                        //Hg
            xAB = 0.2149509937;
            break;
          case 82:                        //Pb
            xAB = 5.0331997871;
            break;
        }
        break;
      case 16:                    //S
        switch (atmB) {
          case 16:                        //S
            xAB = 0.4738560021;
            break;
          case 17:                        //Cl
            xAB = 0.3569709957;
            break;
          case 18:                        //Ar
            xAB = 0.6537690163;
            break;
          case 19:                        //K
            xAB = 29.5289516449;
            break;
          case 20:                        //Ca
            xAB = 0.3257389963;
            break;
          case 21:                        //Sc
            xAB = 0.2118500024;
            break;
          case 22:                        //Ti
            xAB = 0.9437839985;
            break;
          case 23:                        //V
            xAB = 2.0350389481;
            break;
          case 24:                        //Cr
            xAB = 0.5503339767;
            break;
          case 25:                        //Mn
            xAB = 1.6126500368;
            break;
          case 26:                        //Fe
            xAB = 0.0334779993;
            break;
          case 27:                        //Co
            xAB = 0.4367069900;
            break;
          case 28:                        //Ni
            xAB = 0.2748520076;
            break;
          case 29:                        //Cu
            xAB = 0.0052479999;
            break;
          case 30:                        //Zn
            xAB = 0.6480000019;
            break;
          case 31:                        //Ga
            xAB = 2.4562840462;
            break;
          case 32:                        //Ge
            xAB = 1.1609569788;
            break;
          case 33:                        //As
            xAB = 1.0337840319;
            break;
          case 34:                        //Se
            xAB = 0.3917210102;
            break;
          case 35:                        //Br
            xAB = 0.2866880000;
            break;
          case 37:                        //Rb
            xAB = 0.9644110203;
            break;
          case 38:                        //Sr
            xAB = 3.1296029091;
            break;
          case 39:                        //Y
            xAB = 0.3182219863;
            break;
          case 40:                        //Zr
            xAB = 2.4293239117;
            break;
          case 41:                        //Nb
            xAB = 2.4600200653;
            break;
          case 42:                        //Mo
            xAB = 0.8304280043;
            break;
          case 43:                        //Tc
            xAB = 1.4965020418;
            break;
          case 44:                        //Ru
            xAB = 1.0066829920;
            break;
          case 45:                        //Rh
            xAB = 4.8164100647;
            break;
          case 46:                        //Pd
            xAB = 0.2552289963;
            break;
          case 47:                        //Ag
            xAB = 11.1880216599;
            break;
          case 48:                        //Cd
            xAB = 0.3613890111;
            break;
          case 49:                        //In
            xAB = 6.3411049843;
            break;
          case 50:                        //Sn
            xAB = 2.3148701191;
            break;
          case 51:                        //Sb
            xAB = 0.3969689906;
            break;
          case 52:                        //Te
            xAB = 0.9439250231;
            break;
          case 53:                        //I
            xAB = 0.7099289894;
            break;
          case 55:                        //Cs
            xAB = 0.0917429999;
            break;
          case 56:                        //Ba
            xAB = 0.2153860033;
            break;
          case 57:                        //La
            xAB = 2.6824119091;
            break;
          case 72:                        //Hf
            xAB = 1.6667599678;
            break;
          case 73:                        //Ta
            xAB = 2.9759800434;
            break;
          case 74:                        //W
            xAB = 0.4886299968;
            break;
          case 75:                        //Re
            xAB = 0.0688740015;
            break;
          case 76:                        //Os
            xAB = 4.1860499382;
            break;
          case 77:                        //Ir
            xAB = 4.1357321739;
            break;
          case 78:                        //Pt
            xAB = 2.2242629528;
            break;
          case 79:                        //Au
            xAB = 21.6191444397;
            break;
          case 80:                        //Hg
            xAB = 0.3473829925;
            break;
          case 81:                        //Tl
            xAB = 0.3034259975;
            break;
          case 82:                        //Pb
            xAB = 0.1751500070;
            break;
          case 83:                        //Bi
            xAB = 0.6209970117;
            break;
        }
        break;
      case 17:                    //Cl
        switch (atmB) {
          case 17:                        //Cl
            xAB = 0.3329190016;
            break;
          case 18:                        //Ar
            xAB = 2.2560939789;
            break;
          case 19:                        //K
            xAB = 1.2314809561;
            break;
          case 20:                        //Ca
            xAB = 5.9561438560;
            break;
          case 21:                        //Sc
            xAB = 2.9961290359;
            break;
          case 22:                        //Ti
            xAB = 0.3332970142;
            break;
          case 23:                        //V
            xAB = 0.2436569929;
            break;
          case 24:                        //Cr
            xAB = 0.3690730035;
            break;
          case 25:                        //Mn
            xAB = 0.2018499970;
            break;
          case 26:                        //Fe
            xAB = 0.0194729995;
            break;
          case 27:                        //Co
            xAB = 1.0334140062;
            break;
          case 28:                        //Ni
            xAB = 0.2023130059;
            break;
          case 29:                        //Cu
            xAB = 0.1390649974;
            break;
          case 30:                        //Zn
            xAB = 0.7213510275;
            break;
          case 31:                        //Ga
            xAB = 1.1866610050;
            break;
          case 32:                        //Ge
            xAB = 0.5452389717;
            break;
          case 33:                        //As
            xAB = 0.4544329941;
            break;
          case 34:                        //Se
            xAB = 0.9810670018;
            break;
          case 35:                        //Br
            xAB = 0.3892379999;
            break;
          case 36:                        //Kr
            xAB = 0.5202170014;
            break;
          case 37:                        //Rb
            xAB = 10.3844861984;
            break;
          case 38:                        //Sr
            xAB = 3.7839550972;
            break;
          case 39:                        //Y
            xAB = 1.1549589634;
            break;
          case 40:                        //Zr
            xAB = 0.4430989921;
            break;
          case 41:                        //Nb
            xAB = 1.8915569782;
            break;
          case 42:                        //Mo
            xAB = 0.4743250012;
            break;
          case 43:                        //Tc
            xAB = 1.6515829563;
            break;
          case 44:                        //Ru
            xAB = 0.1265860051;
            break;
          case 45:                        //Rh
            xAB = 3.5868649483;
            break;
          case 46:                        //Pd
            xAB = 2.9698910713;
            break;
          case 47:                        //Ag
            xAB = 23.7654590607;
            break;
          case 48:                        //Cd
            xAB = 0.1404239982;
            break;
          case 49:                        //In
            xAB = 2.3885519505;
            break;
          case 50:                        //Sn
            xAB = 0.8077819943;
            break;
          case 51:                        //Sb
            xAB = 0.1564749926;
            break;
          case 52:                        //Te
            xAB = 0.2854779959;
            break;
          case 53:                        //I
            xAB = 0.3104740083;
            break;
          case 54:                        //Xe
            xAB = 0.4637579918;
            break;
          case 55:                        //Cs
            xAB = 4.5319647789;
            break;
          case 56:                        //Ba
            xAB = 0.1601770073;
            break;
          case 57:                        //La
            xAB = 0.2302030027;
            break;
          case 71:                        //Lu
            xAB = 12.7570991516;
            break;
          case 72:                        //Hf
            xAB = 0.7064210176;
            break;
          case 73:                        //Ta
            xAB = 0.5164129734;
            break;
          case 74:                        //W
            xAB = 0.2779999971;
            break;
          case 75:                        //Re
            xAB = 0.1463190019;
            break;
          case 76:                        //Os
            xAB = 0.3279199898;
            break;
          case 77:                        //Ir
            xAB = 0.2589159906;
            break;
          case 78:                        //Pt
            xAB = 0.3410010040;
            break;
          case 79:                        //Au
            xAB = 2.2577021122;
            break;
          case 80:                        //Hg
            xAB = 0.0536600016;
            break;
          case 81:                        //Tl
            xAB = 0.1620370001;
            break;
          case 82:                        //Pb
            xAB = 0.1648139954;
            break;
          case 83:                        //Bi
            xAB = 0.3268710077;
            break;
        }
        break;
      case 18:                    //Ar
        switch (atmB) {
          case 18:                        //Ar
            xAB = 0.9726989865;
            break;
          case 19:                        //K
            xAB = 9.7105083466;
            break;
          case 20:                        //Ca
            xAB = 0.2910720110;
            break;
          case 35:                        //Br
            xAB = 3.2626678944;
            break;
          case 36:                        //Kr
            xAB = 0.5548740029;
            break;
          case 37:                        //Rb
            xAB = 18.4333286285;
            break;
          case 53:                        //I
            xAB = 0.3053669930;
            break;
          case 54:                        //Xe
            xAB = 0.0492659993;
            break;
        }
        break;
      case 19:                    //K
        switch (atmB) {
          case 19:                        //K
            xAB = 5.9343290329;
            break;
          case 20:                        //Ca
            xAB = 1.2403199673;
            break;
          case 23:                        //V
            xAB = 2.0265901089;
            break;
          case 24:                        //Cr
            xAB = 2.0000000000;
            break;
          case 26:                        //Fe
            xAB = 6.0000000000;
            break;
          case 35:                        //Br
            xAB = 3.3227949142;
            break;
          case 36:                        //Kr
            xAB = 8.6097822189;
            break;
          case 41:                        //Nb
            xAB = 2.0265901089;
            break;
          case 42:                        //Mo
            xAB = 2.1423900127;
            break;
          case 53:                        //I
            xAB = 4.8243532181;
            break;
          case 54:                        //Xe
            xAB = 1.2248890400;
            break;
          case 73:                        //Ta
            xAB = 2.0267000198;
            break;
          case 74:                        //W
            xAB = 2.0266098976;
            break;
          case 76:                        //Os
            xAB = 2.0263199806;
            break;
          case 77:                        //Ir
            xAB = 2.0264000893;
            break;
        }
        break;
      case 20:                    //Ca
        switch (atmB) {
          case 20:                        //Ca
            xAB = 30.0035915375;
            break;
          case 22:                        //Ti
            xAB = 4.1091408730;
            break;
          case 25:                        //Mn
            xAB = 0.6201800108;
            break;
          case 30:                        //Zn
            xAB = 1.2402900457;
            break;
          case 35:                        //Br
            xAB = 4.0529098511;
            break;
          case 36:                        //Kr
            xAB = 0.8788909912;
            break;
          case 53:                        //I
            xAB = 7.6899209023;
            break;
          case 54:                        //Xe
            xAB = 1.7171210051;
            break;
          case 72:                        //Hf
            xAB = 4.3195099831;
            break;
        }
        break;
      case 21:                    //Sc
        switch (atmB) {
          case 21:                        //Sc
            xAB = 2.5981659889;
            break;
          case 35:                        //Br
            xAB = 2.0982511044;
            break;
          case 53:                        //I
            xAB = 3.1142818928;
            break;
        }
        break;
      case 22:                    //Ti
        switch (atmB) {
          case 22:                        //Ti
            xAB = 2.0000000000;
            break;
          case 33:                        //As
            xAB = 1.5813169479;
            break;
          case 35:                        //Br
            xAB = 0.8834339976;
            break;
          case 38:                        //Sr
            xAB = 2.8172500134;
            break;
          case 53:                        //I
            xAB = 2.4267470837;
            break;
          case 56:                        //Ba
            xAB = 9.4935302734;
            break;
          case 80:                        //Hg
            xAB = 2.9572000504;
            break;
        }
        break;
      case 23:                    //V
        switch (atmB) {
          case 23:                        //V
            xAB = 10.7798919678;
            break;
          case 35:                        //Br
            xAB = 0.6126980186;
            break;
          case 53:                        //I
            xAB = 6.1981120110;
            break;
          case 82:                        //Pb
            xAB = 1.0000000000;
            break;
        }
        break;
      case 24:                    //Cr
        switch (atmB) {
          case 24:                        //Cr
            xAB = 10.3186073303;
            break;
          case 35:                        //Br
            xAB = 0.2178529948;
            break;
          case 42:                        //Mo
            xAB = 1.7419430017;
            break;
          case 53:                        //I
            xAB = 2.5985898972;
            break;
          case 82:                        //Pb
            xAB = 1.0291099548;
            break;
        }
        break;
      case 25:                    //Mn
        switch (atmB) {
          case 25:                        //Mn
            xAB = 2.4600400925;
            break;
          case 32:                        //Ge
            xAB = 2.2551510334;
            break;
          case 34:                        //Se
            xAB = 2.1807200909;
            break;
          case 35:                        //Br
            xAB = 1.1835800409;
            break;
          case 51:                        //Sb
            xAB = 2.2367100716;
            break;
          case 53:                        //I
            xAB = 1.1934100389;
            break;
        }
        break;
      case 26:                    //Fe
        switch (atmB) {
          case 26:                        //Fe
            xAB = 1.8468899727;
            break;
          case 35:                        //Br
            xAB = 6.0619211197;
            break;
          case 53:                        //I
            xAB = 0.5326219797;
            break;
        }
        break;
      case 27:                    //Co
        switch (atmB) {
          case 27:                        //Co
            xAB = 3.9196178913;
            break;
          case 32:                        //Ge
            xAB = 2.1518499851;
            break;
          case 33:                        //As
            xAB = 1.6752400398;
            break;
          case 34:                        //Se
            xAB = 2.2024099827;
            break;
          case 35:                        //Br
            xAB = 0.4251480103;
            break;
          case 51:                        //Sb
            xAB = 2.2760500908;
            break;
          case 53:                        //I
            xAB = 1.1052390337;
            break;
        }
        break;
      case 28:                    //Ni
        switch (atmB) {
          case 28:                        //Ni
            xAB = 0.0354739986;
            break;
          case 35:                        //Br
            xAB = 0.6321449876;
            break;
          case 53:                        //I
            xAB = 0.0174589995;
            break;
        }
        break;
      case 29:                    //Cu
        switch (atmB) {
          case 29:                        //Cu
            xAB = 5.1843757629;
            break;
          case 35:                        //Br
            xAB = 0.7685170174;
            break;
          case 53:                        //I
            xAB = 0.0067810002;
            break;
        }
        break;
      case 30:                    //Zn
        switch (atmB) {
          case 30:                        //Zn
            xAB = 0.4650000036;
            break;
          case 33:                        //As
            xAB = 3.1565709114;
            break;
          case 34:                        //Se
            xAB = 0.5115939975;
            break;
          case 35:                        //Br
            xAB = 0.7470269799;
            break;
          case 52:                        //Te
            xAB = 1.7560700178;
            break;
          case 53:                        //I
            xAB = 0.9766070247;
            break;
          case 82:                        //Pb
            xAB = 1.0000000000;
            break;
        }
        break;
      case 31:                    //Ga
        switch (atmB) {
          case 31:                        //Ga
            xAB = 1.1983939409;
            break;
          case 33:                        //As
            xAB = 1.6862980127;
            break;
          case 35:                        //Br
            xAB = 1.2610360384;
            break;
          case 49:                        //In
            xAB = 2.4219870567;
            break;
          case 53:                        //I
            xAB = 1.2521680593;
            break;
        }
        break;
      case 32:                    //Ge
        switch (atmB) {
          case 32:                        //Ge
            xAB = 3.0230000019;
            break;
          case 34:                        //Se
            xAB = 5.8720507622;
            break;
          case 35:                        //Br
            xAB = 0.6277369857;
            break;
          case 43:                        //Tc
            xAB = 2.1520600319;
            break;
          case 44:                        //Ru
            xAB = 2.1515600681;
            break;
          case 45:                        //Rh
            xAB = 2.1521399021;
            break;
          case 50:                        //Sn
            xAB = 12.3434114456;
            break;
          case 52:                        //Te
            xAB = 7.0190491676;
            break;
          case 53:                        //I
            xAB = 1.3232669830;
            break;
          case 75:                        //Re
            xAB = 2.1515800953;
            break;
        }
        break;
      case 33:                    //As
        switch (atmB) {
          case 33:                        //As
            xAB = 0.7373070121;
            break;
          case 34:                        //Se
            xAB = 0.7112609744;
            break;
          case 35:                        //Br
            xAB = 0.5141530037;
            break;
          case 49:                        //In
            xAB = 6.2083501816;
            break;
          case 52:                        //Te
            xAB = 0.6857740283;
            break;
          case 53:                        //I
            xAB = 0.3108240068;
            break;
          case 72:                        //Hf
            xAB = 1.2808200121;
            break;
        }
        break;
      case 34:                    //Se
        switch (atmB) {
          case 34:                        //Se
            xAB = 0.8218230009;
            break;
          case 35:                        //Br
            xAB = 0.3193419874;
            break;
          case 43:                        //Tc
            xAB = 2.2026200294;
            break;
          case 44:                        //Ru
            xAB = 2.2021200657;
            break;
          case 45:                        //Rh
            xAB = 2.2026998997;
            break;
          case 49:                        //In
            xAB = 2.3193230629;
            break;
          case 50:                        //Sn
            xAB = 3.0618851185;
            break;
          case 52:                        //Te
            xAB = 1.1878260374;
            break;
          case 75:                        //Re
            xAB = 2.2021400928;
            break;
          case 82:                        //Pb
            xAB = 0.1111949980;
            break;
          case 83:                        //Bi
            xAB = 0.6512079835;
            break;
        }
        break;
      case 35:                    //Br
        switch (atmB) {
          case 35:                        //Br
            xAB = 0.6153079867;
            break;
          case 36:                        //Kr
            xAB = 0.3080979884;
            break;
          case 37:                        //Rb
            xAB = 5.1762142181;
            break;
          case 38:                        //Sr
            xAB = 2.7665669918;
            break;
          case 39:                        //Y
            xAB = 1.0543160439;
            break;
          case 40:                        //Zr
            xAB = 0.8589090109;
            break;
          case 41:                        //Nb
            xAB = 1.9212690592;
            break;
          case 42:                        //Mo
            xAB = 0.2259179950;
            break;
          case 43:                        //Tc
            xAB = 3.8201301098;
            break;
          case 44:                        //Ru
            xAB = 0.6598809958;
            break;
          case 45:                        //Rh
            xAB = 1.5101490021;
            break;
          case 46:                        //Pd
            xAB = 4.6380510330;
            break;
          case 47:                        //Ag
            xAB = 1.7143690586;
            break;
          case 48:                        //Cd
            xAB = 0.2722670138;
            break;
          case 49:                        //In
            xAB = 3.7285981178;
            break;
          case 50:                        //Sn
            xAB = 0.6687980294;
            break;
          case 51:                        //Sb
            xAB = 0.1980440021;
            break;
          case 52:                        //Te
            xAB = 0.3942019939;
            break;
          case 53:                        //I
            xAB = 0.4830540121;
            break;
          case 54:                        //Xe
            xAB = 0.4751160145;
            break;
          case 55:                        //Cs
            xAB = 1.6584270000;
            break;
          case 56:                        //Ba
            xAB = 0.8287940025;
            break;
          case 57:                        //La
            xAB = 0.2385820001;
            break;
          case 71:                        //Lu
            xAB = 8.6482744217;
            break;
          case 72:                        //Hf
            xAB = 0.6924560070;
            break;
          case 73:                        //Ta
            xAB = 0.7914450169;
            break;
          case 74:                        //W
            xAB = 0.3723900020;
            break;
          case 75:                        //Re
            xAB = 0.2875280082;
            break;
          case 76:                        //Os
            xAB = 0.3824299872;
            break;
          case 77:                        //Ir
            xAB = 0.1718789935;
            break;
          case 78:                        //Pt
            xAB = 0.5203610063;
            break;
          case 79:                        //Au
            xAB = 0.0042369999;
            break;
          case 80:                        //Hg
            xAB = 0.1723629981;
            break;
          case 81:                        //Tl
            xAB = 0.2968359888;
            break;
          case 82:                        //Pb
            xAB = 0.1482290030;
            break;
          case 83:                        //Bi
            xAB = 0.3811700046;
            break;
        }
        break;
      case 36:                    //Kr
        switch (atmB) {
          case 36:                        //Kr
            xAB = 0.0520990007;
            break;
          case 37:                        //Rb
            xAB = 15.3075027466;
            break;
          case 53:                        //I
            xAB = 0.2011359930;
            break;
          case 54:                        //Xe
            xAB = 0.0497930013;
            break;
        }
        break;
      case 37:                    //Rb
        switch (atmB) {
          case 37:                        //Rb
            xAB = 20.1476097107;
            break;
          case 53:                        //I
            xAB = 4.0924458504;
            break;
          case 54:                        //Xe
            xAB = 0.9749649763;
            break;
        }
        break;
      case 38:                    //Sr
        switch (atmB) {
          case 38:                        //Sr
            xAB = 5.3721199036;
            break;
          case 53:                        //I
            xAB = 2.1039409637;
            break;
        }
        break;
      case 39:                    //Y
        switch (atmB) {
          case 39:                        //Y
            xAB = 1.6917250156;
            break;
          case 53:                        //I
            xAB = 1.0214020014;
            break;
        }
        break;
      case 40:                    //Zr
        switch (atmB) {
          case 40:                        //Zr
            xAB = 3.0777730942;
            break;
          case 53:                        //I
            xAB = 4.5139431953;
            break;
        }
        break;
      case 41:                    //Nb
        switch (atmB) {
          case 41:                        //Nb
            xAB = 2.1223878860;
            break;
          case 53:                        //I
            xAB = 2.3992979527;
            break;
          case 82:                        //Pb
            xAB = 1.0000000000;
            break;
        }
        break;
      case 42:                    //Mo
        switch (atmB) {
          case 42:                        //Mo
            xAB = 0.6264619827;
            break;
          case 53:                        //I
            xAB = 0.1246950030;
            break;
          case 82:                        //Pb
            xAB = 5.0000000000;
            break;
        }
        break;
      case 43:                    //Tc
        switch (atmB) {
          case 51:                        //Sb
            xAB = 2.2762598991;
            break;
          case 53:                        //I
            xAB = 0.1105939969;
            break;
        }
        break;
      case 44:                    //Ru
        switch (atmB) {
          case 44:                        //Ru
            xAB = 0.0978050008;
            break;
          case 51:                        //Sb
            xAB = 2.2757599354;
            break;
          case 53:                        //I
            xAB = 7.8377099037;
            break;
        }
        break;
      case 45:                    //Rh
        switch (atmB) {
          case 45:                        //Rh
            xAB = 2.0701138973;
            break;
          case 51:                        //Sb
            xAB = 2.2763400078;
            break;
          case 53:                        //I
            xAB = 2.3579440117;
            break;
        }
        break;
      case 46:                    //Pd
        switch (atmB) {
          case 46:                        //Pd
            xAB = 0.0519560017;
            break;
          case 53:                        //I
            xAB = 3.5220170021;
            break;
        }
        break;
      case 47:                    //Ag
        switch (atmB) {
          case 47:                        //Ag
            xAB = 0.5577419996;
            break;
          case 53:                        //I
            xAB = 0.0489040017;
            break;
        }
        break;
      case 48:                    //Cd
        switch (atmB) {
          case 48:                        //Cd
            xAB = 18.6179981232;
            break;
          case 52:                        //Te
            xAB = 1.0859190226;
            break;
          case 53:                        //I
            xAB = 0.3967840075;
            break;
        }
        break;
      case 49:                    //In
        switch (atmB) {
          case 49:                        //In
            xAB = 8.0634908676;
            break;
          case 51:                        //Sb
            xAB = 6.6608009338;
            break;
          case 52:                        //Te
            xAB = 2.0395820141;
            break;
          case 53:                        //I
            xAB = 5.9478211403;
            break;
        }
        break;
      case 50:                    //Sn
        switch (atmB) {
          case 50:                        //Sn
            xAB = 0.2870000005;
            break;
          case 52:                        //Te
            xAB = 2.9519760609;
            break;
          case 53:                        //I
            xAB = 1.7831629515;
            break;
        }
        break;
      case 51:                    //Sb
        switch (atmB) {
          case 51:                        //Sb
            xAB = 0.7248849869;
            break;
          case 53:                        //I
            xAB = 0.3181900084;
            break;
          case 75:                        //Re
            xAB = 2.2757799625;
            break;
        }
        break;
      case 52:                    //Te
        switch (atmB) {
          case 52:                        //Te
            xAB = 0.6424859762;
            break;
          case 53:                        //I
            xAB = 1.1011159420;
            break;
          case 80:                        //Hg
            xAB = 0.2127320021;
            break;
          case 82:                        //Pb
            xAB = 0.8090419769;
            break;
        }
        break;
      case 53:                    //I
        switch (atmB) {
          case 53:                        //I
            xAB = 0.5105419755;
            break;
          case 54:                        //Xe
            xAB = 0.1120899990;
            break;
          case 55:                        //Cs
            xAB = 1.0721780062;
            break;
          case 56:                        //Ba
            xAB = 0.8355969787;
            break;
          case 57:                        //La
            xAB = 0.2268829942;
            break;
          case 71:                        //Lu
            xAB = 10.0823154449;
            break;
          case 72:                        //Hf
            xAB = 0.8209480047;
            break;
          case 73:                        //Ta
            xAB = 6.5515508652;
            break;
          case 74:                        //W
            xAB = 1.0773700476;
            break;
          case 75:                        //Re
            xAB = 3.5592861176;
            break;
          case 76:                        //Os
            xAB = 2.1991899014;
            break;
          case 77:                        //Ir
            xAB = 1.4971480370;
            break;
          case 78:                        //Pt
            xAB = 4.6012477875;
            break;
          case 79:                        //Au
            xAB = 0.0088160001;
            break;
          case 80:                        //Hg
            xAB = 0.3420580029;
            break;
          case 81:                        //Tl
            xAB = 0.4300329983;
            break;
          case 82:                        //Pb
            xAB = 0.2674260139;
            break;
          case 83:                        //Bi
            xAB = 0.8623769879;
            break;
        }
        break;
      case 54:                    //Xe
        switch (atmB) {
          case 54:                        //Xe
            xAB = 0.3444739878;
            break;
        }
        break;
      case 55:                    //Cs
        switch (atmB) {
          case 55:                        //Cs
            xAB = 25.3200550079;
            break;
        }
        break;
      case 56:                    //Ba
        switch (atmB) {
          case 56:                        //Ba
            xAB = 0.3561860025;
            break;
          case 72:                        //Hf
            xAB = 9.0225200653;
            break;
        }
        break;
      case 57:                    //La
        switch (atmB) {
          case 57:                        //La
            xAB = 5.1751618385;
            break;
        }
        break;
      case 72:                    //Hf
        switch (atmB) {
          case 72:                        //Hf
            xAB = 1.0589109659;
            break;
        }
        break;
      case 73:                    //Ta
        switch (atmB) {
          case 73:                        //Ta
            xAB = 10.9870529175;
            break;
        }
        break;
      case 74:                    //W
        switch (atmB) {
          case 74:                        //W
            xAB = 7.4713897705;
            break;
        }
        break;
      case 75:                    //Re
        switch (atmB) {
          case 75:                        //Re
            xAB = 4.4888520241;
            break;
        }
        break;
      case 76:                    //Os
        switch (atmB) {
          case 76:                        //Os
            xAB = 0.8304399848;
            break;
        }
        break;
      case 77:                    //Ir
        switch (atmB) {
          case 77:                        //Ir
            xAB = 11.1751928329;
            break;
        }
        break;
      case 78:                    //Pt
        switch (atmB) {
          case 78:                        //Pt
            xAB = 9.0102519989;
            break;
        }
        break;
      case 79:                    //Au
        switch (atmB) {
          case 79:                        //Au
            xAB = 0.0130909998;
            break;
        }
        break;
      case 80:                    //Hg
        switch (atmB) {
          case 80:                        //Hg
            xAB = 0.4232760072;
            break;
        }
        break;
      case 81:                    //Tl
        switch (atmB) {
          case 81:                        //Tl
            xAB = 9.5351266861;
            break;
        }
        break;
      case 82:                    //Pb
        switch (atmB) {
          case 82:                        //Pb
            xAB = 2.3623430729;
            break;
        }
        break;
      case 83:                    //Bi
        switch (atmB) {
          case 83:                        //Bi
            xAB = 1.1682139635;
            break;
        }
        break;
    }
    if (xAB == 0.0) {std::cout << "WARNING: MNDOd.hpp: PM6: xAB(): interaction parameter xAB between elements " << atmA << " and " << atmB << " is zero" << std::endl;}
    return xAB;
  }
  std::vector<double> AM1K(int atomicnr, int atm2 = 0) {
    //function returning the a constants for PM6; parameters dimensionless
    std::vector<double> aconst(1,0.0);
    switch (atomicnr) {
      case 1:     //H
        aconst[0] = 0.024184;
        break;
      case 4:     //Be
        aconst[0] = 0.164180;
        break;
      case 6:     //C
        aconst[0] = 0.046302;
        break;
      case 7:     //N
        aconst[0] = -0.001436;
        break;
      case 8:     //O
        aconst[0] = -0.017771;
        break;
      case 9:     //F
        aconst[0] = -0.010792;
        break;
      case 11:    //Na
        aconst[0] = -1.026036;
        break;
      case 13:    //Al
        aconst[0] = 1.002222;
        break;
      case 14:    //Si
        aconst[0] = 0.208571;
        break;
      case 15:    //P
        aconst[0] = -0.034320;
        break;
      case 16:    //S
        aconst[0] = -0.036928;
        break;
      case 17:    //Cl
        aconst[0] = -0.013213;
        break;
      case 19:    //K
        aconst[0] = 0.157519;
        break;
      case 20:    //Ca
        aconst[0] = -0.025275;
        break;
      case 35:    //Br
        aconst[0] = -0.004996;
        break;
      case 38:    //Sr
        aconst[0] = -0.012948;
        break;
      case 50:    //Sn
        aconst[0] = -1.004587;
        break;
      case 53:    //I
        aconst[0] = -0.035519;
        break;
      case 82:    //Pb
        aconst[0] = -0.239463;
        break;
    }
    aconst[0] /= (au2eV*dist_Angstrom2au);
    return aconst;
  }
  std::vector<double> AM1L(int atomicnr, int atm2 = 0) {
    //function returning the B constants for PM6; parameters in Angstrom^-2
    std::vector<double> bconst(1,0.0);
    switch (atomicnr) {
      case 1:     //H
        bconst[0] = 3.055953;
        break;
      case 4:     //Be
        bconst[0] = 1.704828;
        break;
      case 6:     //C
        bconst[0] = 2.100206;
        break;
      case 7:     //N
        bconst[0] = 0.495196;
        break;
      case 8:     //O
        bconst[0] = 3.058310;
        break;
      case 9:     //F
        bconst[0] = 6.004648;
        break;
      case 11:    //Na
        bconst[0] = 2.014506;
        break;
      case 13:    //Al
        bconst[0] = 1.517400;
        break;
      case 14:    //Si
        bconst[0] = 6.000483;
        break;
      case 15:    //P
        bconst[0] = 6.001394;
        break;
      case 16:    //S
        bconst[0] = 1.795067;
        break;
      case 17:    //Cl
        bconst[0] = 3.687022;
        break;
      case 19:    //K
        bconst[0] = 6.000566;
        break;
      case 20:    //Ca
        bconst[0] = 0.500017;
        break;
      case 35:    //Br
        bconst[0] = 6.001292;
        break;
      case 38:    //Sr
        bconst[0] = 6.000126;
        break;
      case 50:    //Sn
        bconst[0] = 4.706252;
        break;
      case 53:    //I
        bconst[0] = 1.744389;
        break;
      case 82:    //Pb
        bconst[0] = 5.444338;
        break;
    }
    return bconst;
  }
  std::vector<double> AM1M(int atomicnr, int atm2 = 0) {
    //function returning the K constants for PM6; parameters in Angstrom
    std::vector<double> cconst(1,0.0);
    switch (atomicnr) {
      case 1:     //H
        cconst[0] = 1.786011;
        break;
      case 4:     //Be
        cconst[0] = 1.785591;
        break;
      case 6:     //C
        cconst[0] = 1.333959;
        break;
      case 7:     //N
        cconst[0] = 1.704857;
        break;
      case 8:     //O
        cconst[0] = 1.896435;
        break;
      case 9:     //F
        cconst[0] = 1.847724;
        break;
      case 11:    //Na
        cconst[0] = 1.271202;
        break;
      case 13:    //Al
        cconst[0] = 0.659101;
        break;
      case 14:    //Si
        cconst[0] = 1.185245;
        break;
      case 15:    //P
        cconst[0] = 2.296737;
        break;
      case 16:    //S
        cconst[0] = 2.082618;
        break;
      case 17:    //Cl
        cconst[0] = 2.544635;
        break;
      case 19:    //K
        cconst[0] = 2.047539;
        break;
      case 20:    //Ca
        cconst[0] = 2.329051;
        break;
      case 35:    //Br
        cconst[0] = 2.895153;
        break;
      case 38:    //Sr
        cconst[0] = 3.011964;
        break;
      case 50:    //Sn
        cconst[0] = 1.180218;
        break;
      case 53:    //I
        cconst[0] = 1.223844;
        break;
      case 82:    //Pb
        cconst[0] = 1.613682;
        break;
    }
    return cconst;
  }
  double Dvalue(size_t atmnr, size_t idx) {
    //function returning the D values needed to calculate eris; values stored in Angstrom, but returned in a.u.
    //note that this function returns both D1 and D2; idx is then either 1 or 2
    //idx = 1 and idx = 2 are to be used with the sp integrals;
    //idx = 1 and idx > 2 are to be used with d integrals only
    double dval = 0.0;
    switch (atmnr) {
      case 1:                             //H
        dval = 0.0;
        break;
      case 2:                             //He
          if (idx == 1) {dval = 0.131014703960855;}               //ss,sp
          else if (idx == 2) {dval = 0.112082033016238;}          //pp
          else if (idx == 3) {dval = 0.158507931189912;}          //D2pp  
          break;
      case 3:                             //Li
        if (idx == 1) {dval = 0.188309649221985;}               //ss,sp
        else if (idx == 2) {dval = 0.219441050004328;}          //pp
        else if (idx == 3) {dval = 0.310336509057513;}          //D2pp  
        break;
      case 4:                             //Be
        if (idx == 1) {dval = 0.612722978546848;}               //ss,sp
        else if (idx == 2) {dval = 0.507727122900611;}          //pp
        else if (idx == 3) {dval = 0.718034583190716;}          //D2pp  
        break;
      case 5:                             //B
        if (idx == 1) {dval = 0.487625290387018;}               //ss,sp
        else if (idx == 2) {dval = 0.438148500995496;}          //pp
        else if (idx == 3) {dval = 0.619635552441271;}          //D2pp  
        break;
      case 6:                             //C
        if (idx == 1) {dval = 0.398769026644354;}               //ss,sp
        else if (idx == 2) {dval = 0.380603398632070;}          //pp
        else if (idx == 3) {dval = 0.538254488230767;}          //D2pp  
        break;
      case 7:                             //N
        if (idx == 1) {dval = 0.342228398510806;}               //ss,sp
        else if (idx == 2) {dval = 0.324175750222850;}          //pp
        else if (idx == 3) {dval = 0.458453742557628;}          //D2pp  
        break;
      case 8:                             //O
        if (idx == 1) {dval = 0.125474829922673;}               //ss,sp
        else if (idx == 2) {dval = 0.285389030159066;}          //pp
        else if (idx == 3) {dval = 0.403601037003455;}          //D2pp 
        break;
      case 9:                             //F
        if (idx == 1) {dval = 0.122984050377085;}               //ss,sp
        else if (idx == 2) {dval = 0.222968371908298;}          //pp
        else if (idx == 3) {dval = 0.315324895532964;}          //D2pp  
        break;
      case 10:                             //Ne
        if (idx == 1) {dval = 0.137178115574348;}               //ss,sp
        else if (idx == 2) {dval = 0.169018735012505;}          //pp
        else if (idx == 3) {dval = 0.239028587349829;}          //D2pp
        break;
      case 11:                             //Na
        if (idx == 1) {dval = 1.191910304078041;}               //ss,sp
        else if (idx == 2) {dval = 0.932020480939228;}          //pp
        else if (idx == 3) {dval = 1.318076004553750;}          //D2pp
        break;
      case 12:                             //Mg
        if (idx == 1) {dval = 0.789854607260177;}               //ss,sp
        else if (idx == 2) {dval = 0.637543917428701;}          //pp
        else if (idx == 3) {dval = 0.901623254636141;}          //D2pp
        break;
      case 13:                             //Al
        if (idx == 1) {dval = 0.480350836667587;}               //ss,sp
        else if (idx == 2) {dval = 0.506249969575799;}          //pp
        else if (idx == 3) {dval = 0.715945572925063;}          //D2pp  
        else if (idx == 4) {dval = 0.662987287477118;}          //D2sd  
        else if (idx == 5) {dval = 0.501810053828850;}          //D1pd  
        else if (idx == 6) {dval = 0.833754338018994;}          //D2dd
        break;
      case 14:                             //Si
        if (idx == 1) {dval = 0.639069999484421;}               //ss,sp
        else if (idx == 2) {dval = 0.738879530082677;}          //pp
        else if (idx == 3) {dval = 1.044933452402780;}          //D2pp  
        else if (idx == 4) {dval = 0.721183686547496;}          //D2sd  
        else if (idx == 5) {dval = 0.374491926537813;}          //D1pd  
        else if (idx == 6) {dval = 0.497208445490473;}          //D2dd
        break;
      case 15:                             //P
        if (idx == 1) {dval = 0.524793731734139;}               //ss,sp
        else if (idx == 2) {dval = 0.490479002762893;}          //pp
        else if (idx == 3) {dval = 0.693642057766515;}          //D2pp  
        else if (idx == 4) {dval = 0.732787812248743;}          //D2sd  
        else if (idx == 5) {dval = 0.480199244273924;}          //D1pd  
        else if (idx == 6) {dval = 0.860200377948453;}          //D2dd
        break;
      case 16:                             //S
        if (idx == 1) {dval = 0.516188446928841;}               //ss,sp
        else if (idx == 2) {dval = 0.480958891630322;}          //pp
        else if (idx == 3) {dval = 0.680178587487533;}          //D2pp  
        else if (idx == 4) {dval = 0.508945261758166;}          //D2sd  
        else if (idx == 5) {dval = 0.263856371554635;}          //D1pd  
        else if (idx == 6) {dval = 0.340372443635254;}          //D2dd
        break;
      case 17:                             //Cl
        if (idx == 1) {dval = 0.431281700062551;}               //ss,sp
        else if (idx == 2) {dval = 0.418046175421794;}          //pp
        else if (idx == 3) {dval = 0.591206570979703;}          //D2pp  
        else if (idx == 4) {dval = 0.586016203572245;}          //D2sd  
        else if (idx == 5) {dval = 0.397418415578926;}          //D1pd  
        else if (idx == 6) {dval = 0.799341418689642;}          //D2dd
        break;
      case 18:                             //Ar
        if (idx == 1) {dval = 0.178962942795338;}               //ss,sp
        else if (idx == 2) {dval = 0.148841407168558;}          //pp
        else if (idx == 3) {dval = 0.210493536660471;}          //D2pp
        break;
      case 19:                             //K
        if (idx == 1) {dval = 0.022659672695510;}               //ss,sp
        else if (idx == 2) {dval = 0.995610989352146;}          //pp
        else if (idx == 3) {dval = 1.408006563989500;}          //D2pp
        break;
      case 20:                             //Ca
        if (idx == 1) {dval = 0.693392764275083;}               //ss,sp
        else if (idx == 2) {dval = 0.544904444810534;}          //pp
        else if (idx == 3) {dval = 0.770611256048439;}          //D2pp
        break;
      case 21:                             //Sc
        if (idx == 1) {dval = 0.998780205422073;}               //ss,sp
        else if (idx == 2) {dval = 0.834491313776961;}          //pp
        else if (idx == 3) {dval = 1.180148933625920;}          //D2pp  
        else if (idx == 4) {dval = 0.857158212122765;}          //D2sd  
        else if (idx == 5) {dval = 0.423629829188227;}          //D1pd  
        else if (idx == 6) {dval = 0.569310158628294;}          //D2dd
        break;
      case 22:                             //Ti
        if (idx == 1) {dval = 0.039096507147143;}               //ss,sp
        else if (idx == 2) {dval = 0.964337459089600;}          //pp
        else if (idx == 3) {dval = 1.363779113348922;}          //D2pp  
        else if (idx == 4) {dval = 0.284353822907186;}          //D2sd  
        else if (idx == 5) {dval = 0.284353822907186;}          //D1pd  
        else if (idx == 6) {dval = 0.746223888521239;}          //D2dd
        break;
      case 23:                             //V
        if (idx == 1) {dval = 0.592195396665819;}               //ss,sp
        else if (idx == 2) {dval = 1.055919520092552;}          //pp
        else if (idx == 3) {dval = 1.493295706089377;}          //D2pp  
        else if (idx == 4) {dval = 0.951442253215763;}          //D2sd  
        else if (idx == 5) {dval = 0.584448198350205;}          //D1pd  
        else if (idx == 6) {dval = 0.758782523599628;}          //D2dd
        break;
      case 24:                             //Cr
        if (idx == 1) {dval = 0.560965789105230;}               //ss,sp
        else if (idx == 2) {dval = 1.077160672276439;}          //pp
        else if (idx == 3) {dval = 1.523335231588261;}          //D2pp  
        else if (idx == 4) {dval = 0.973692855659483;}          //D2sd  
        else if (idx == 5) {dval = 0.746907997990843;}          //D1pd  
        else if (idx == 6) {dval = 0.905922144700335;}          //D2dd
        break;
      case 25:                             //Mn
        if (idx == 1) {dval = 0.663329241260328;}               //ss,sp
        else if (idx == 2) {dval = 0.735676709392294;}          //pp
        else if (idx == 3) {dval = 1.040403979944592;}          //D2pp  
        else if (idx == 4) {dval = 0.613525472800320;}          //D2sd  
        else if (idx == 5) {dval = 0.246728956220636;}          //D1pd  
        else if (idx == 6) {dval = 0.405841865408353;}          //D2dd
        break;
      case 26:                             //Fe
        if (idx == 1) {dval = 0.047445342043916;}               //ss,sp
        else if (idx == 2) {dval = 0.187022387507528;}          //pp
        else if (idx == 3) {dval = 0.264489596880543;}          //D2pp  
        else if (idx == 4) {dval = 1.255491159254260;}          //D2sd  
        else if (idx == 5) {dval = 0.042164309095138;}          //D1pd  
        else if (idx == 6) {dval = 0.979280457509392;}          //D2dd
        break;
      case 27:                             //Co
        if (idx == 1) {dval = 0.212175025225166;}               //ss,sp
        else if (idx == 2) {dval = 0.351315665017731;}          //pp
        else if (idx == 3) {dval = 0.496835378142199;}          //D2pp  
        else if (idx == 4) {dval = 0.824297853325991;}          //D2sd  
        else if (idx == 5) {dval = 0.343856576350667;}          //D1pd  
        else if (idx == 6) {dval = 0.568660439974243;}          //D2dd
        break;
      case 28:                             //Ni
        if (idx == 1) {dval = 0.605427298567827;}               //ss,sp
        else if (idx == 2) {dval = 0.487063557881180;}          //pp
        else if (idx == 3) {dval = 0.688811889293258;}          //D2pp  
        else if (idx == 4) {dval = 0.611098093029402;}          //D2sd  
        else if (idx == 5) {dval = 0.349115099391215;}          //D1pd  
        else if (idx == 6) {dval = 0.420856859404095;}          //D2dd
        break;
      case 29:                             //Cu
        if (idx == 1) {dval = 0.095316034880562;}               //ss,sp
        else if (idx == 2) {dval = 0.213250563320384;}          //pp
        else if (idx == 3) {dval = 0.301581838831389;}          //D2pp  
        else if (idx == 4) {dval = 0.135726037221871;}          //D2sd  
        else if (idx == 5) {dval = 0.144696523039409;}          //D1pd  
        else if (idx == 6) {dval = 0.176326691688215;}          //D2dd
        break;
      case 30:                             //Zn
        if (idx == 1) {dval = 0.806676637759663;}               //ss,sp
        else if (idx == 2) {dval = 0.627306883962796;}          //pp
        else if (idx == 3) {dval = 0.887145903070192;}          //D2pp
        break;
      case 31:                             //Ga
        if (idx == 1) {dval = 0.610208892260151;}               //ss,sp
        else if (idx == 2) {dval = 0.649028428281070;}          //pp
        else if (idx == 3) {dval = 0.917864805640783;}          //D2pp
        break;
      case 32:                             //Ge
        if (idx == 1) {dval = 0.541074302446772;}               //ss,sp
        else if (idx == 2) {dval = 0.656798708891373;}          //pp
        else if (idx == 3) {dval = 0.928853641863318;}          //D2pp
        break;
      case 33:                             //As
        if (idx == 1) {dval = 0.441036817004972;}               //ss,sp
        else if (idx == 2) {dval = 0.635939327431146;}          //pp
        else if (idx == 3) {dval = 0.899354021699552;}          //D2pp  
        else if (idx == 4) {dval = 0.616714539344592;}          //D2sd  
        else if (idx == 5) {dval = 0.633232641294486;}          //D1pd  
        else if (idx == 6) {dval = 0.963773169199262;}          //D2dd 
        break;
      case 34:                             //Se
        if (idx == 1) {dval = 0.574938532247933;}               //ss,sp
        else if (idx == 2) {dval = 0.559159094015625;}          //pp
        else if (idx == 3) {dval = 0.790770374281150;}          //D2pp 
        break;
      case 35:                             //Br
        if (idx == 1) {dval = 0.192783500671668;}               //ss,sp
        else if (idx == 2) {dval = 0.551454136136752;}          //pp
        else if (idx == 3) {dval = 0.779873918351334;}          //D2pp  
        else if (idx == 4) {dval = 0.297033902207035;}          //D2sd  
        else if (idx == 5) {dval = 0.544466297370430;}          //D1pd  
        else if (idx == 6) {dval = 0.882105037514291;}          //D2dd
        break;
      case 36:                             //Kr
        if (idx == 1) {dval = 0.095135118668742;}               //ss,sp
        else if (idx == 2) {dval = 0.249935794065445;}          //pp
        else if (idx == 3) {dval = 0.353462589689841;}          //D2pp
        break;
      case 37:                             //Rb
        if (idx == 1) {dval = 0.038009280558182;}               //ss,sp
        else if (idx == 2) {dval = 1.018208065792171;}          //pp
        else if (idx == 3) {dval = 1.439963655960966;}          //D2pp
        break;
      case 38:                             //Sr
        if (idx == 1) {dval = 0.791199250136399;}               //ss,sp
        else if (idx == 2) {dval = 0.785764863247092;}          //pp
        else if (idx == 3) {dval = 1.111239326440278;}          //D2pp
        break;
      case 39:                             //Y
        if (idx == 1) {dval = 0.522395350305567;}               //ss,sp
        else if (idx == 2) {dval = 0.912144924862712;}          //pp
        else if (idx == 3) {dval = 1.289967723590635;}          //D2pp  
        else if (idx == 4) {dval = 0.715433059296060;}          //D2sd  
        else if (idx == 5) {dval = 0.670265748745441;}          //D1pd  
        else if (idx == 6) {dval = 0.812717182346403;}          //D2dd
        break;
      case 40:                             //Zr
        if (idx == 1) {dval = 0.992092095353882;}               //ss,sp
        else if (idx == 2) {dval = 0.802093356369124;}          //pp
        else if (idx == 3) {dval = 1.134331302866574;}          //D2pp  
        else if (idx == 4) {dval = 1.210252746486354;}          //D2sd  
        else if (idx == 5) {dval = 0.710194117556750;}          //D1pd  
        else if (idx == 6) {dval = 0.856013752344914;}          //D2dd
        break;
      case 41:                             //Nb
        if (idx == 1) {dval = 0.613260425654214;}               //ss,sp
        else if (idx == 2) {dval = 0.980224963320348;}          //pp
        else if (idx == 3) {dval = 1.386247437304305;}          //D2pp  
        else if (idx == 4) {dval = 0.918907517971167;}          //D2sd  
        else if (idx == 5) {dval = 0.477973796751287;}          //D1pd  
        else if (idx == 6) {dval = 0.678547930089049;}          //D2dd
        break;
      case 42:                             //Mo
        if (idx == 1) {dval = 1.286630024468437;}               //ss,sp
        else if (idx == 2) {dval = 1.006715626937359;}          //pp
        else if (idx == 3) {dval = 1.423710893067747;}          //D2pp  
        else if (idx == 4) {dval = 0.978053525250524;}          //D2sd  
        else if (idx == 5) {dval = 0.542064048447551;}          //D1pd  
        else if (idx == 6) {dval = 0.734317181775462;}          //D2dd
        break;
      case 43:                             //Tc
        if (idx == 1) {dval = 0.081340541545497;}               //ss,sp
        else if (idx == 2) {dval = 0.226342521943003;}          //pp
        else if (idx == 3) {dval = 0.320096664273524;}          //D2pp  
        else if (idx == 4) {dval = 1.063197276025183;}          //D2sd  
        else if (idx == 5) {dval = 0.091192223734839;}          //D1pd  
        else if (idx == 6) {dval = 0.759160050762380;}          //D2dd
        break;
      case 44:                             //Ru
        if (idx == 1) {dval = 0.048976806682061;}               //ss,sp
        else if (idx == 2) {dval = 0.245517701669803;}          //pp
        else if (idx == 3) {dval = 0.347214463504107;}          //D2pp  
        else if (idx == 4) {dval = 0.921548425651970;}          //D2sd  
        else if (idx == 5) {dval = 0.153196957386564;}          //D1pd  
        else if (idx == 6) {dval = 0.640995692318136;}          //D2dd
        break;
      case 45:                             //Rh
        if (idx == 1) {dval = 0.097773028061432;}               //ss,sp
        else if (idx == 2) {dval = 0.315709665450736;}          //pp
        else if (idx == 3) {dval = 0.446480890652703;}          //D2pp  
        else if (idx == 4) {dval = 0.519532482099693;}          //D2sd  
        else if (idx == 5) {dval = 0.312681290266133;}          //D1pd  
        else if (idx == 6) {dval = 0.462434112053053;}          //D2dd
        break;
      case 46:                             //Pd
        if (idx == 1) {dval = 0.999548504038035;}               //ss,sp
        else if (idx == 2) {dval = 1.175291525854818;}          //pp
        else if (idx == 3) {dval = 1.662113215606053;}          //D2pp  
        else if (idx == 4) {dval = 0.883205139061224;}          //D2sd  
        else if (idx == 5) {dval = 0.284779357155356;}          //D1pd  
        else if (idx == 6) {dval = 0.604411315535252;}          //D2dd
        break;
      case 47:                             //Ag
        if (idx == 1) {dval = 0.276658887367092;}               //ss,sp
        else if (idx == 2) {dval = 1.993908722140594;}          //pp
        else if (idx == 3) {dval = 2.819812756985235;}          //D2pp  
        else if (idx == 4) {dval = 0.178569018592361;}          //D2sd  
        else if (idx == 5) {dval = 0.000744217229847;}          //D1pd  
        else if (idx == 6) {dval = 0.223345405364148;}          //D2dd
        break;
      case 48:                             //Cd
        if (idx == 1) {dval = 0.853338347418201;}               //ss,sp
        else if (idx == 2) {dval = 0.694529393236754;}          //pp
        else if (idx == 3) {dval = 0.982212887382173;}          //D2pp
        break;
      case 49:                             //In
        if (idx == 1) {dval = 0.811963883835729;}               //ss,sp
        else if (idx == 2) {dval = 0.645338102685790;}          //pp
        else if (idx == 3) {dval = 0.912645897134365;}          //D2pp
        break;
      case 50:                             //Sn
        if (idx == 1) {dval = 0.734455814835314;}               //ss,sp
        else if (idx == 2) {dval = 0.660613041595511;}          //pp
        else if (idx == 3) {dval = 0.934247922904913;}          //D2pp
        break;
      case 51:                             //Sb
        if (idx == 1) {dval = 0.713963189747905;}               //ss,sp
        else if (idx == 2) {dval = 0.766766081560770;}          //pp
        else if (idx == 3) {dval = 1.084370991710915;}          //D2pp  
        else if (idx == 4) {dval = 0.898985181808335;}          //D2sd  
        else if (idx == 5) {dval = 0.529231075981591;}          //D1pd  
        else if (idx == 6) {dval = 0.659027016472211;}          //D2dd
        break;
      case 52:                             //Te
        if (idx == 1) {dval = 0.552635993384303;}               //ss,sp
        else if (idx == 2) {dval = 0.785228408631646;}          //pp
        else if (idx == 3) {dval = 1.110480665047517;}          //D2pp
        break;
      case 53:                             //I
        if (idx == 1) {dval = 0.198280184739571;}               //ss,sp
        else if (idx == 2) {dval = 0.709144394787329;}          //pp
        else if (idx == 3) {dval = 1.002881620789100;}          //D2pp  
        else if (idx == 4) {dval = 0.411419500028873;}          //D2sd  
        else if (idx == 5) {dval = 0.685994486460873;}          //D1pd  
        else if (idx == 6) {dval = 0.866527348937416;}          //D2dd
        break;
      case 54:                             //Xe
        if (idx == 1) {dval = 0.741536220077310;}               //ss,sp
        else if (idx == 2) {dval = 0.687493293472354;}          //pp
        else if (idx == 3) {dval = 0.972262339669150;}          //D2pp
        break;
      case 55:                             //Cs
        if (idx == 1) {dval = 0.039699084980781;}               //ss,sp
        else if (idx == 2) {dval = 0.985700319100715;}          //pp
        else if (idx == 3) {dval = 1.393990759707719;}          //D2pp
        break;
      case 56:                             //Ba
        if (idx == 1) {dval = 1.404295537613053;}               //ss,sp
        else if (idx == 2) {dval = 1.116204006239129;}          //pp
        else if (idx == 3) {dval = 1.578550843998559;}          //D2pp
        break;
      case 57:                             //La
        if (idx == 1) {dval = 0.403140917419645;}               //ss,sp
        else if (idx == 2) {dval = 1.278911322359718;}          //pp
        else if (idx == 3) {dval = 1.808653737153623;}          //D2pp  
        else if (idx == 4) {dval = 1.021364323711469;}          //D2sd  
        else if (idx == 5) {dval = 0.694476053191406;}          //D1pd  
        else if (idx == 6) {dval = 0.962292424882071;}          //D2dd
        break;
      case 71:                             //Lu
        if (idx == 1) {dval = 0.069070466180648;}               //ss,sp
        else if (idx == 2) {dval = 0.932272738637958;}          //pp
        else if (idx == 3) {dval = 1.318432750812508;}          //D2pp  
        else if (idx == 4) {dval = 0.420176985081264;}          //D2sd  
        else if (idx == 5) {dval = 0.546316383542788;}          //D1pd  
        else if (idx == 6) {dval = 0.729995175661586;}          //D2dd
        break;
      case 72:                             //Hf
        if (idx == 1) {dval = 0.414684925090808;}               //ss,sp
        else if (idx == 2) {dval = 1.013014109665400;}          //pp
        else if (idx == 3) {dval = 1.432618292764115;}          //D2pp  
        else if (idx == 4) {dval = 0.880734147561127;}          //D2sd  
        else if (idx == 5) {dval = 0.710102707189282;}          //D1pd  
        else if (idx == 6) {dval = 0.882689653388518;}          //D2dd
        break;
      case 73:                             //Ta
        if (idx == 1) {dval = 0.419526662486526;}               //ss,sp
        else if (idx == 2) {dval = 0.329734853537401;}          //pp
        else if (idx == 3) {dval = 0.466315501859699;}          //D2pp  
        else if (idx == 4) {dval = 0.497409212442830;}          //D2sd  
        else if (idx == 5) {dval = 0.170428307876671;}          //D1pd  
        else if (idx == 6) {dval = 0.883933798709380;}          //D2dd
        break;
      case 74:                             //W
        if (idx == 1) {dval = 0.624296751090172;}               //ss,sp
        else if (idx == 2) {dval = 0.982953849593797;}          //pp
        else if (idx == 3) {dval = 1.390106665282390;}          //D2pp  
        else if (idx == 4) {dval = 1.024033740915533;}          //D2sd  
        else if (idx == 5) {dval = 0.745472395237006;}          //D1pd  
        else if (idx == 6) {dval = 0.905534118113977;}          //D2dd
        break;
      case 75:                             //Re
        if (idx == 1) {dval = 0.824436165409608;}               //ss,sp
        else if (idx == 2) {dval = 0.879348886952893;}          //pp
        else if (idx == 3) {dval = 1.243587121986467;}          //D2pp  
        else if (idx == 4) {dval = 0.929855298385377;}          //D2sd  
        else if (idx == 5) {dval = 0.452286827114394;}          //D1pd  
        else if (idx == 6) {dval = 0.644090819974472;}          //D2dd
        break;
      case 76:                             //Os
        if (idx == 1) {dval = 0.443888479145403;}               //ss,sp
        else if (idx == 2) {dval = 1.001484906320624;}          //pp
        else if (idx == 3) {dval = 1.416313537030574;}          //D2pp  
        else if (idx == 4) {dval = 0.894076717516388;}          //D2sd  
        else if (idx == 5) {dval = 0.751316453687105;}          //D1pd  
        else if (idx == 6) {dval = 0.915137348312778;}          //D2dd
        break;
      case 77:                             //Ir
        if (idx == 1) {dval = 0.145754360806309;}               //ss,sp
        else if (idx == 2) {dval = 0.388743760315690;}          //pp
        else if (idx == 3) {dval = 0.549766698126365;}          //D2pp  
        else if (idx == 4) {dval = 0.751844748459638;}          //D2sd  
        else if (idx == 5) {dval = 0.377969997910932;}          //D1pd  
        else if (idx == 6) {dval = 0.607198013167825;}          //D2dd
        break;
      case 78:                             //Pt
        if (idx == 1) {dval = 0.844464961964862;}               //ss,sp
        else if (idx == 2) {dval = 0.960252069460144;}          //pp
        else if (idx == 3) {dval = 1.358001499927367;}          //D2pp  
        else if (idx == 4) {dval = 0.726580963102600;}          //D2sd  
        else if (idx == 5) {dval = 0.220663692667707;}          //D1pd  
        else if (idx == 6) {dval = 0.512769426134045;}          //D2dd
        break;
      case 79:                             //Au
        if (idx == 1) {dval = 0.994548596311460;}               //ss,sp
        else if (idx == 2) {dval = 1.031869499942355;}          //pp
        else if (idx == 3) {dval = 1.459283841417622;}          //D2pp  
        else if (idx == 4) {dval = 0.315733007119631;}          //D2sd  
        else if (idx == 5) {dval = 0.038624555458488;}          //D1pd  
        else if (idx == 6) {dval = 0.337612769551325;}          //D2dd
        break;
      case 80:                             //Hg
        if (idx == 1) {dval = 0.921608833249811;}               //ss,sp
        else if (idx == 2) {dval = 1.052782596291628;}          //pp
        else if (idx == 3) {dval = 1.488859425905980;}          //D2pp
        break;
      case 81:                             //Tl
        if (idx == 1) {dval = 0.407862734015207;}               //ss,sp
        else if (idx == 2) {dval = 0.903850191620500;}          //pp
        else if (idx == 3) {dval = 1.278237199343232;}          //D2pp
        break;
      case 82:                             //Pb
        if (idx == 1) {dval = 0.812170704100768;}               //ss,sp
        else if (idx == 2) {dval = 0.947236712787819;}          //pp
        else if (idx == 3) {dval = 1.339595006002242;}          //D2pp
        break;
      case 83:                             //Bi
        if (idx == 1) {dval = 0.339526463166489;}               //ss,sp
        else if (idx == 2) {dval = 0.852589788684787;}          //pp
        else if (idx == 3) {dval = 1.205744042298837;}          //D2pp
        break;
    }
    return dval*dist_Angstrom2aum1;
  }
  double rhocore(size_t atmnr) {
    //function returning rho values used to calculate VAB integrals; values already in atomic units
    double rho_ = 0.0;
    if (atmnr == 21) {rho_ = 3.173734;}                //Sc
    else if (atmnr == 26) {rho_ = 1.272092;}           //Fe
    else if (atmnr == 28) {rho_ = 1.586979;}           //Ni
    else if (atmnr == 39) {rho_ = 2.773703;}           //Y
    else if (atmnr == 57) {rho_ = 2.511701;}           //La
    else if (atmnr == 71) {rho_ = 2.743262;}           //Lu
    else {rho_ = rho(atmnr,0);}                        //all other atoms
    return rho_;
  }
  double rho(size_t atmnr, size_t l) {
    //function returning the rho values needed to calculate eris; values stored in Angstrom but returned in atomic units
    double rho = 0.0;
    switch (atmnr) {
      case 1:                            //H
        rho = 0.498302887257530;
        break;
      case 2:                            //He
        if (l == 0) {rho = 0.762265117375052;}                  //ss
        else if (l == 1) {rho = 0.459664717769081;}             //sp
        else if (l == 2) {rho = 0.270856004379672;}             //pp
        break;
      case 3:                            //Li
        if (l == 0) {rho = 0.652399657851181;}                  //ss
        else if (l == 1) {rho = 0.313310125718249;}             //sp
        else if (l == 2) {rho = 0.224949823556194;}             //pp
        break;
      case 4:                            //Be
        if (l == 0) {rho = 0.953264767744198;}                  //ss
        else if (l == 1) {rho = 0.651308034778323;}             //sp
        else if (l == 2) {rho = 0.848519282029161;}             //pp
        break;
      case 5:                            //B
        if (l == 0) {rho = 0.880244747208541;}                  //ss
        else if (l == 1) {rho = 0.617848903574928;}             //sp
        else if (l == 2) {rho = 0.428829230126313;}             //pp
        break;
      case 6:                            //C
        if (l == 0) {rho = 0.539898143512633;}                  //ss
        else if (l == 1) {rho = 0.683614519401577;}             //sp
        else if (l == 2) {rho = 0.403590997141624;}             //pp
        break;
      case 7:                            //N
        if (l == 0) {rho = 0.582650060854241;}                  //ss
        else if (l == 1) {rho = 0.353250111120616;}             //sp
        else if (l == 2) {rho = 0.313774918523130;}             //pp
        break;
      case 8:                            //O
        if (l == 0) {rho = 0.636924557682770;}                  //ss
        else if (l == 1) {rho = 0.157088724234194;}             //sp
        else if (l == 2) {rho = 0.241955538461576;}             //pp
        break;
      case 9:                            //F
        if (l == 0) {rho = 0.578446792656361;}                  //ss
        else if (l == 1) {rho = 0.201847298438142;}             //sp
        else if (l == 2) {rho = 0.457966011494171;}             //pp
        break;
      case 10:                           //Ne
        if (l == 0) {rho = 0.359998765517578;}                  //ss
        else if (l == 1) {rho = 0.249803837461762;}             //sp
        else if (l == 2) {rho = 0.370982799790365;}             //pp
        break;
      case 11:                           //Na
        if (l == 0) {rho = 1.773367390434583;}                  //ss
        else if (l == 1) {rho = 1.374163035210596;}             //sp
        else if (l == 2) {rho = 1.322512003384908;}             //pp
        break;
      case 12:                           //Mg
        if (l == 0) {rho = 1.011874920014573;}                  //ss
        else if (l == 1) {rho = 0.949167694576259;}             //sp
        else if (l == 2) {rho = 1.003427815833684;}             //pp
        break;
      case 13:                           //Al
        if (l == 0) {rho = 1.082329252832721;}                  //ss
        else if (l == 1) {rho = 0.927281783420203;}             //sp
        else if (l == 2) {rho = 0.532552385374852;}             //pp
        else if (l == 3) {rho = 0.237051546320658;}             //sd
        else if (l == 4) {rho = 0.304970478680607;}             //pd
        else if (l == 5) {rho = 0.143734266923325;}             //dd0
        else if (l == 6) {rho = 0.374760792259656;}             //dd2
        break;
      case 14:                           //Si
        if (l == 0) {rho = 1.385965777517626;}                  //ss
        else if (l == 1) {rho = 0.682371382633331;}             //sp
        else if (l == 2) {rho = 0.907645849768342;}             //pp
        else if (l == 3) {rho = 2.645886041529756;}             //sd
        else if (l == 4) {rho = 0.850838326626353;}             //pd
        else if (l == 5) {rho = 1.446498578547937;}             //dd0
        else if (l == 6) {rho = 0.571202107909823;}             //dd2
        break;
      case 15:                           //P
        if (l == 0) {rho = 0.822004831553053;}                  //ss
        else if (l == 1) {rho = 0.748434355835323;}             //sp
        else if (l == 2) {rho = 0.533162648439199;}             //pp
        else if (l == 3) {rho = 0.205761257690349;}             //sd
        else if (l == 4) {rho = 0.724526081964487;}             //pd
        else if (l == 5) {rho = 0.143400215838420;}             //dd0
        else if (l == 6) {rho = 0.380911885982379;}             //dd2
        break;
      case 16:                           //S
        if (l == 0) {rho = 0.785119646565011;}                  //ss
        else if (l == 1) {rho = 0.330990693343205;}             //sp
        else if (l == 2) {rho = 0.538337898063135;}             //pp
        else if (l == 3) {rho = 1.709370452933764;}             //sd
        else if (l == 4) {rho = 1.001291902164801;}             //pd
        else if (l == 5) {rho = 0.237404983186697;}             //dd0
        else if (l == 6) {rho = 0.258667776125239;}             //dd2
        break;
      case 17:                           //Cl
        if (l == 0) {rho = 0.646149647191544;}                  //ss
        else if (l == 1) {rho = 0.304479462072066;}             //sp
        else if (l == 2) {rho = 0.416196195116528;}             //pp
        else if (l == 3) {rho = 1.239805731529637;}             //sd
        else if (l == 4) {rho = 0.548927782501082;}             //pd
        else if (l == 5) {rho = 0.159896556377262;}             //dd0
        else if (l == 6) {rho = 0.383283714391607;}             //dd2
        break;
      case 18:                           //Ar
        if (l == 0) {rho = 0.403153158473876;}                  //ss
        else if (l == 1) {rho = 0.200218599516930;}             //sp
        else if (l == 2) {rho = 0.336658539305281;}             //pp
        break;
      case 19:                           //K
        if (l == 0) {rho = 2.136920624458504;}                  //ss
        else if (l == 1) {rho = 0.144573342012369;}             //sp
        else if (l == 2) {rho = 1.386825273692823;}             //pp
        break;
      case 20:                           //Ca
        if (l == 0) {rho = 1.257441039118639;}                  //ss
        else if (l == 1) {rho = 0.759183441213040;}             //sp
        else if (l == 2) {rho = 0.893990114937423;}             //pp
        break;
      case 21:                           //Sc
        if (l == 0) {rho = 1.552282812051898;}                  //ss
        else if (l == 1) {rho = 1.983793682506516;}             //sp
        else if (l == 2) {rho = 0.587655562643245;}             //pp
        else if (l == 3) {rho = 0.452404319748588;}             //sd
        else if (l == 4) {rho = 0.846990865858787;}             //pd
        else if (l == 5) {rho = 1.298486509223875;}             //dd0
        else if (l == 6) {rho = 0.612903329194269;}             //dd2
        break;
      case 22:                           //Ti
        if (l == 0) {rho = 1.259183135051761;}                  //ss
        else if (l == 1) {rho = 0.122118184066781;}             //sp
        else if (l == 2) {rho = 0.880882715944765;}             //pp
        else if (l == 3) {rho = 0.255701049254668;}             //sd
        else if (l == 4) {rho = 0.765475174254869;}             //pd
        else if (l == 5) {rho = 1.427670493922329;}             //dd0
        else if (l == 6) {rho = 0.765657877744421;}             //dd2
        break;
      case 23:                           //V
        if (l == 0) {rho = 1.203356534333159;}                  //ss
        else if (l == 1) {rho = 0.793432123309974;}             //sp
        else if (l == 2) {rho = 1.056491002242171;}             //pp
        else if (l == 3) {rho = 0.712904075495129;}             //sd
        else if (l == 4) {rho = 0.998494237354275;}             //pd
        else if (l == 5) {rho = 0.932344606268033;}             //dd0
        else if (l == 6) {rho = 0.686298549396277;}             //dd2
        break;
      case 24:                           //Cr
        if (l == 0) {rho = 0.829185580291525;}                  //ss
        else if (l == 1) {rho = 0.767000025600732;}             //sp
        else if (l == 2) {rho = 0.993232122655000;}             //pp
        else if (l == 3) {rho = 0.795814321027551;}             //sd
        else if (l == 4) {rho = 1.100011771494041;}             //pd
        else if (l == 5) {rho = 0.719397477209990;}             //dd0
        else if (l == 6) {rho = 0.717530992991552;}             //dd2
        break;
      case 25:                           //Mn
        if (l == 0) {rho = 1.162951722718861;}                  //ss
        else if (l == 1) {rho = 0.676754315147623;}             //sp
        else if (l == 2) {rho = 0.675553175160790;}             //pp
        else if (l == 3) {rho = 0.610685207328563;}             //sd
        else if (l == 4) {rho = 0.400960969420759;}             //pd
        else if (l == 5) {rho = 1.064868207115544;}             //dd0
        else if (l == 6) {rho = 0.453664396365077;}             //dd2
        break;
      case 26:                           //Fe
        if (l == 0) {rho = 0.902568531855843;}                  //ss
        else if (l == 1) {rho = 0.124859969435664;}             //sp
        else if (l == 2) {rho = 0.263400497024990;}             //pp
        else if (l == 3) {rho = 0.897407347278094;}             //sd
        else if (l == 4) {rho = 0.155970596926469;}             //pd
        else if (l == 5) {rho = 0.474112875592390;}             //dd0
        else if (l == 6) {rho = 0.658786667194947;}             //dd2
        break;
      case 27:                           //Co
        if (l == 0) {rho = 2.691090442938072;}                  //ss
        else if (l == 1) {rho = 0.485043006926006;}             //sp
        else if (l == 2) {rho = 0.510743092362624;}             //pp
        else if (l == 3) {rho = 0.681392684079912;}             //sd
        else if (l == 4) {rho = 0.749069176475342;}             //pd
        else if (l == 5) {rho = 2.941456891431999;}             //dd0
        else if (l == 6) {rho = 0.756923712141128;}             //dd2
        break;
      case 28:                           //Ni
        if (l == 0) {rho = 1.764283466769012;}                  //ss
        else if (l == 1) {rho = 0.772090515734335;}             //sp
        else if (l == 2) {rho = 0.621710613603150;}             //pp
        else if (l == 3) {rho = 0.521712843215622;}             //sd
        else if (l == 4) {rho = 0.929794701659848;}             //pd
        else if (l == 5) {rho = 0.701401033129800;}             //dd0
        else if (l == 6) {rho = 0.415511720842509;}             //dd2
        break;
      case 29:                           //Cu
        if (l == 0) {rho = 0.673426928192856;}                  //ss
        else if (l == 1) {rho = 0.191411948872111;}             //sp
        else if (l == 2) {rho = 0.224239069093513;}             //pp
        else if (l == 3) {rho = 0.119195701231726;}             //sd
        else if (l == 4) {rho = 0.273762253800866;}             //pd
        else if (l == 5) {rho = 0.709713060079813;}             //dd0
        else if (l == 6) {rho = 0.220331309475356;}             //dd2
        break;
      case 30:                           //Zn
        if (l == 0) {rho = 0.826860154148625;}                  //ss
        else if (l == 1) {rho = 1.080450455611649;}             //sp
        else if (l == 2) {rho = 0.177492460950341;}             //pp
        break;
      case 31:                           //Ga
        if (l == 0) {rho = 0.695306799725680;}                  //ss
        else if (l == 1) {rho = 0.693004184424064;}             //sp
        else if (l == 2) {rho = 1.016644214076514;}             //pp
        break;
      case 32:                           //Ge
        if (l == 0) {rho = 0.957639492071076;}                  //ss
        else if (l == 1) {rho = 1.160242699497683;}             //sp
        else if (l == 2) {rho = 0.692844666467850;}             //pp
        break;
      case 33:                           //As
        if (l == 0) {rho = 1.080238491181203;}                  //ss
        else if (l == 1) {rho = 1.031972666586546;}             //sp
        else if (l == 2) {rho = 0.731612599005236;}             //pp
        else if (l == 3) {rho = 0.502440358466687;}             //sd
        else if (l == 4) {rho = 0.456958576085662;}             //pd
        else if (l == 5) {rho = 0.283013142517176;}             //dd0
        else if (l == 6) {rho = 0.527650113989328;}             //dd2
        break;
      case 34:                           //Se
        if (l == 0) {rho = 1.303759111306377;}                  //ss
        else if (l == 1) {rho = 0.448052879294207;}             //sp
        else if (l == 2) {rho = 0.518644989138839;}             //pp
        break;
      case 35:                           //Br
        if (l == 0) {rho = 0.945256598333530;}                  //ss
        else if (l == 1) {rho = 0.200779742100256;}             //sp
        else if (l == 2) {rho = 0.514251341643729;}             //pp
        else if (l == 3) {rho = 0.193392441234431;}             //sd
        else if (l == 4) {rho = 0.437729657484038;}             //pd
        else if (l == 5) {rho = 0.467015346407475;}             //dd0
        else if (l == 6) {rho = 0.599559161709227;}             //dd2
        break;
      case 36:                           //Kr
        if (l == 0) {rho = 0.359993671498624;}                  //ss
        else if (l == 1) {rho = 0.372795119547878;}             //sp
        else if (l == 2) {rho = 0.499303151335682;}             //pp
        break;
      case 37:                           //Rb
        if (l == 0) {rho = 1.077684721357344;}                  //ss
        else if (l == 1) {rho = 0.085737772279648;}             //sp
        else if (l == 2) {rho = 0.587092343015745;}             //pp
        break;
      case 38:                           //Sr
        if (l == 0) {rho = 1.563932978357558;}                  //ss
        else if (l == 1) {rho = 0.971880504204776;}             //sp
        else if (l == 2) {rho = 1.168833697362299;}             //pp
        break;
      case 39:                           //Y
        if (l == 0) {rho = 1.779168888263733;}                  //ss
        else if (l == 1) {rho = 0.823905421778281;}             //sp
        else if (l == 2) {rho = 0.800558363943937;}             //pp
        else if (l == 3) {rho = 0.416335233974931;}             //sd
        else if (l == 4) {rho = 0.723657858538749;}             //pd
        else if (l == 5) {rho = 0.951414914136958;}             //dd0
        else if (l == 6) {rho = 0.711941341787087;}             //dd2
        break;
      case 40:                           //Zr
        if (l == 0) {rho = 1.350504799548732;}                  //ss
        else if (l == 1) {rho = 1.160825390667300;}             //sp
        else if (l == 2) {rho = 0.902166049039207;}             //pp
        else if (l == 3) {rho = 0.706247379422544;}             //sd
        else if (l == 4) {rho = 1.210104588083458;}             //pd
        else if (l == 5) {rho = 1.106478289765514;}             //dd0
        else if (l == 6) {rho = 0.772217168567077;}             //dd2
        break;
      case 41:                           //Nb
        if (l == 0) {rho = 1.077238370770470;}                  //ss
        else if (l == 1) {rho = 0.923560797930324;}             //sp
        else if (l == 2) {rho = 1.001030448044894;}             //pp
        else if (l == 3) {rho = 0.822953963660536;}             //sd
        else if (l == 4) {rho = 0.996389658720028;}             //pd
        else if (l == 5) {rho = 0.912142214055408;}             //dd0
        else if (l == 6) {rho = 0.619200864939042;}             //dd2
        break;
      case 42:                           //Mo
        if (l == 0) {rho = 0.839467646224971;}                  //ss
        else if (l == 1) {rho = 1.027452856302049;}             //sp
        else if (l == 2) {rho = 0.881159547946004;}             //pp
        else if (l == 3) {rho = 0.826182086349990;}             //sd
        else if (l == 4) {rho = 0.797096596499414;}             //pd
        else if (l == 5) {rho = 0.701930987130646;}             //dd0
        else if (l == 6) {rho = 0.604890240104249;}             //dd2
        break;
      case 43:                           //Tc
        if (l == 0) {rho = 1.138100534032494;}                  //ss
        else if (l == 1) {rho = 0.203766225064414;}             //sp
        else if (l == 2) {rho = 0.333446639746946;}             //pp
        else if (l == 3) {rho = 0.901161403461861;}             //sd
        else if (l == 4) {rho = 0.239233423982381;}             //pd
        else if (l == 5) {rho = 1.136005387604993;}             //dd0
        else if (l == 6) {rho = 0.714657851707861;}             //dd2
        break;
      case 44:                           //Ru
        if (l == 0) {rho = 1.631265214119752;}                   //ss
        else if (l == 1) {rho = 0.811475378203632;}              //sp
        else if (l == 2) {rho = 0.230101605843808;}              //pp
        else if (l == 3) {rho = 0.455943107867301;}              //sd
        else if (l == 4) {rho = 1.552597526178067;}              //pd
        else if (l == 5) {rho = 1.720801793097790;}              //dd0
        else if (l == 6) {rho = 0.709546521555847;}              //dd2
        break;
      case 45:                           //Rh
        if (l == 0) {rho = 1.982777884408183;}                   //ss
        else if (l == 1) {rho = 2.645886041529756;}              //sp
        else if (l == 2) {rho = 0.230853456176034;}              //pp
        else if (l == 3) {rho = 0.467722216062868;}              //sd
        else if (l == 4) {rho = 2.645886041529756;}              //pd
        else if (l == 5) {rho = 2.047760929938582;}              //dd0
        else if (l == 6) {rho = 0.583142331503590;}              //dd2
        break;
      case 46:                           //Pd
        if (l == 0) {rho = 0.895106647331203;}                   //ss
        else if (l == 1) {rho = 2.645886041529756;}              //sp
        else if (l == 2) {rho = 0.481980534304428;}              //pp
        else if (l == 3) {rho = 0.602485152400958;}              //sd
        else if (l == 4) {rho = 0.866430819003188;}              //pd
        else if (l == 5) {rho = 0.807513753287343;}              //dd0
        else if (l == 6) {rho = 0.550559050315532;}              //dd2
        break;
      case 47:                           //Ag
        if (l == 0) {rho = 2.308936141170490;}                   //ss
        else if (l == 1) {rho = 2.645886041529756;}              //sp
        else if (l == 2) {rho = 0.699119333189849;}              //pp
        else if (l == 3) {rho = 0.247366033668472;}              //sd
        else if (l == 4) {rho = 0.119811056038995;}              //pd
        else if (l == 5) {rho = 2.600054248216507;}              //dd0
        else if (l == 6) {rho = 0.356555930865003;}              //dd2
        break;
      case 48:                           //Cd
        if (l == 0) {rho = 1.078256062027235;}                   //ss
        else if (l == 1) {rho = 1.018524765728737;}              //sp
        else if (l == 2) {rho = 0.239149013147166;}              //pp
        break;
      case 49:                           //In
        if (l == 0) {rho = 0.726807572318632;}                   //ss
        else if (l == 1) {rho = 0.493945855900198;}              //sp
        else if (l == 2) {rho = 1.012405300625597;}              //pp
        break;
      case 50:                           //Sn
        if (l == 0) {rho = 0.870631477477289;}                   //ss
        else if (l == 1) {rho = 1.115966388957370;}              //sp
        else if (l == 2) {rho = 0.705065254853138;}              //pp
        break;
      case 51:                           //Sb
        if (l == 0) {rho = 0.679944865578890;}                   //ss
        else if (l == 1) {rho = 0.940417578598924;}              //sp
        else if (l == 2) {rho = 0.999393923842452;}              //pp
        else if (l == 3) {rho = 0.285291664314757;}              //sd
        else if (l == 4) {rho = 0.327406769415094;}              //pd
        else if (l == 5) {rho = 0.281537286573733;}              //dd0
        else if (l == 6) {rho = 0.408755186989123;}              //dd2
        break;
      case 52:                           //Te
        if (l == 0) {rho = 1.024065559863012;}                   //ss
        else if (l == 1) {rho = 0.341470644797628;}              //sp
        else if (l == 2) {rho = 0.734070195583903;}              //pp
        break;
      case 53:                           //I
        if (l == 0) {rho = 0.995170945000027;}                   //ss
        else if (l == 1) {rho = 0.203784621857441;}              //sp
        else if (l == 2) {rho = 0.529062589367249;}              //pp
        else if (l == 3) {rho = 0.379076673667724;}              //sd
        else if (l == 4) {rho = 0.635529067008505;}              //pd
        else if (l == 5) {rho = 0.423549431263495;}              //dd0
        else if (l == 6) {rho = 0.564016393046812;}              //dd2
        break;
      case 54:                           //Xe
        if (l == 0) {rho = 0.359986561713195;}                   //ss
        else if (l == 1) {rho = 0.416047119197285;}              //sp
        else if (l == 2) {rho = 1.060393270952374;}              //pp
        break;
      case 55:                           //Cs
        if (l == 0) {rho = 1.113704449077382;}                   //ss
        else if (l == 1) {rho = 0.137522858886552;}              //sp
        else if (l == 2) {rho = 0.866138026862086;}              //pp
        break;
      case 56:                           //Ba
        if (l == 0) {rho = 1.999493435494449;}                   //ss
        else if (l == 1) {rho = 1.272475363395139;}              //sp
        else if (l == 2) {rho = 1.500036465719754;}              //pp
        break;
      case 57:                           //La
        if (l == 0) {rho = 1.169858152893597;}                   //ss
        else if (l == 1) {rho = 1.251325159865637;}              //sp
        else if (l == 2) {rho = 0.649224573824363;}              //pp
        else if (l == 3) {rho = 0.422648432123057;}              //sd
        else if (l == 4) {rho = 0.714219941146229;}              //pd
        else if (l == 5) {rho = 0.702572827801111;}              //dd0
        else if (l == 6) {rho = 0.717621013803080;}              //dd2
        break;
      case 71:                           //Lu
        if (l == 0) {rho = 1.159429775127309;}                   //ss
        else if (l == 1) {rho = 0.341615054610765;}              //sp
        else if (l == 2) {rho = 0.569981513320695;}              //pp
        else if (l == 3) {rho = 0.476328800969168;}              //sd
        else if (l == 4) {rho = 2.645886041529756;}              //pd
        else if (l == 5) {rho = 1.741756130182846;}              //dd0
        else if (l == 6) {rho = 0.773431227673866;}              //dd2
        break;
      case 72:                           //Hf
        if (l == 0) {rho = 1.998652272912038;}                   //ss
        else if (l == 1) {rho = 2.645886041529756;}              //sp
        else if (l == 2) {rho = 0.632807501860626;}              //pp
        else if (l == 3) {rho = 0.498349270337745;}              //sd
        else if (l == 4) {rho = 2.645022790314376;}              //pd
        else if (l == 5) {rho = 1.707887317193281;}              //dd0
        else if (l == 6) {rho = 0.882632625345157;}              //dd2
        break;
      case 73:                           //Ta
        if (l == 0) {rho = 1.086834541038764;}                   //ss
        else if (l == 1) {rho = 0.766226732964444;}              //sp
        else if (l == 2) {rho = 0.322644678757355;}              //pp
        else if (l == 3) {rho = 0.438614333453786;}              //sd
        else if (l == 4) {rho = 0.271267111460996;}              //pd
        else if (l == 5) {rho = 0.694833333306993;}              //dd0
        else if (l == 6) {rho = 0.675139534308140;}              //dd2
        break;
      case 74:                           //W
        if (l == 0) {rho = 1.262678170269170;}                   //ss
        else if (l == 1) {rho = 0.696473182579216;}              //sp
        else if (l == 2) {rho = 0.792917558186805;}              //pp
        else if (l == 3) {rho = 0.769749492902473;}              //sd
        else if (l == 4) {rho = 0.721573658959285;}              //pd
        else if (l == 5) {rho = 0.855816787163238;}              //dd0
        else if (l == 6) {rho = 0.733156310203074;}              //dd2
        break;
      case 75:                           //Re
        if (l == 0) {rho = 1.125982809213708;}                   //ss
        else if (l == 1) {rho = 0.841887535093940;}              //sp
        else if (l == 2) {rho = 0.850393889178076;}              //pp
        else if (l == 3) {rho = 0.643819619371760;}              //sd
        else if (l == 4) {rho = 0.669561816896736;}              //pd
        else if (l == 5) {rho = 1.077333041468201;}              //dd0
        else if (l == 6) {rho = 0.618087547802492;}              //dd2
        break;
      case 76:                           //Os
        if (l == 0) {rho = 1.025954249342538;}                   //ss
        else if (l == 1) {rho = 0.541157680262213;}              //sp
        else if (l == 2) {rho = 0.873739055787303;}              //pp
        else if (l == 3) {rho = 0.744961901238136;}              //sd
        else if (l == 4) {rho = 0.847510050388996;}              //pd
        else if (l == 5) {rho = 0.907281546585870;}              //dd0
        else if (l == 6) {rho = 0.752008379958850;}              //dd2
        break;
      case 77:                           //Ir
        if (l == 0) {rho = 2.041074108722596;}                   //ss
        else if (l == 1) {rho = 1.269651963730795;}              //sp
        else if (l == 2) {rho = 0.368466142629950;}              //pp
        else if (l == 3) {rho = 0.517161540642010;}              //sd
        else if (l == 4) {rho = 2.645886041529756;}              //pd
        else if (l == 5) {rho = 2.425801253268405;}              //dd0
        else if (l == 6) {rho = 0.737386880198948;}              //dd2
        break;
      case 78:                           //Pt
        if (l == 0) {rho = 0.833478062183219;}                   //ss
        else if (l == 1) {rho = 0.694201067212267;}              //sp
        else if (l == 2) {rho = 0.783974130827231;}              //pp
        else if (l == 3) {rho = 0.438381105910992;}              //sd
        else if (l == 4) {rho = 0.357003237589604;}              //pd
        else if (l == 5) {rho = 0.937008366089366;}              //dd0
        else if (l == 6) {rho = 0.504569020385002;}              //dd2
        break;
      case 79:                           //Au
        if (l == 0) {rho = 0.918197565701830;}                   //ss
        else if (l == 1) {rho = 1.402097594804687;}              //sp
        else if (l == 2) {rho = 0.573611108217255;}              //pp
        else if (l == 3) {rho = 0.242379105429204;}              //sd
        else if (l == 4) {rho = 0.204273322248840;}              //pd
        else if (l == 5) {rho = 1.005855385755435;}              //dd0
        else if (l == 6) {rho = 0.379722921912888;}              //dd2
        break;
      case 80:                           //Hg
        if (l == 0) {rho = 1.129769817967214;}                   //ss
        else if (l == 1) {rho = 1.012499595573335;}              //sp
        else if (l == 2) {rho = 1.443484339357354;}              //pp
        break;
      case 81:                           //Tl
        if (l == 0) {rho = 1.435623638541992;}                   //ss
        else if (l == 1) {rho = 0.965652788015512;}              //sp
        else if (l == 2) {rho = 1.293560401291496;}              //pp
        break;
      case 82:                           //Pb
        if (l == 0) {rho = 1.370317196474364;}                   //ss
        else if (l == 1) {rho = 0.904427772212870;}              //sp
        else if (l == 2) {rho = 0.742774993680170;}              //pp
        break;
      case 83:                           //Bi
        if (l == 0) {rho = 1.230359591886028;}                   //ss
        else if (l == 1) {rho = 0.856700930375696;}              //sp
        else if (l == 2) {rho = 0.862602078518111;}              //pp
        break;
    }
    return rho*dist_Angstrom2aum1;
  }
  double zn(size_t atmnr, size_t l) {
    //function returning the zn values needed to calculate 1-center eris; values stored and returned directly in atomic units
    double zn_ = 0.0;
    if (atmnr == 13) {                                 //Al
      if (l == 0) {zn_ = 4.742341;}                        //zs
      else if (l == 1) {zn_ = 4.669626;}                   //zp
      else if (l == 2) {zn_ = 7.131138;}                   //zd
    }
    else if (atmnr == 14) {                            //Si
      if (l == 0) {zn_ = 8.388111;}                        //zs
      else if (l == 1) {zn_ = 1.843048;}                   //zp
      else if (l == 2) {zn_ = 0.708600;}                   //zd
    }
    else if (atmnr == 15) {                            //P
      if (l == 0) {zn_ = 6.042706;}                        //zs
      else if (l == 1) {zn_ = 2.376473;}                   //zp
      else if (l == 2) {zn_ = 7.147750;}                   //zd
    }
    else if (atmnr == 16) {                            //S
      if (l == 0) {zn_ = 0.479722;}                        //zs
      else if (l == 1) {zn_ = 1.015507;}                   //zp
      else if (l == 2) {zn_ = 4.317470;}                   //zd
    }
    else if (atmnr == 17) {                            //Cl
      if (l == 0) {zn_ = 0.956297;}                        //zs
      else if (l == 1) {zn_ = 2.464067;}                   //zp
      else if (l == 2) {zn_ = 6.410325;}                   //zd
    }
    else if (atmnr == 21) {                            //Sc
      if (l == 0) {zn_ = 0.848418;}                        //zs
      else if (l == 1) {zn_ = 2.451729;}                   //zp
      else if (l == 2) {zn_ = 0.789372;}                   //zd
    }
    else if (atmnr == 22) {                            //Ti
      if (l == 0) {zn_ = 1.045904;}                        //zs
      else if (l == 1) {zn_ = 1.076844;}                   //zp
      else if (l == 2) {zn_ = 0.717945;}                   //zd
    }
    else if (atmnr == 23) {                            //V
      if (l == 0) {zn_ = 1.094426;}                        //zs
      else if (l == 1) {zn_ = 0.755378;}                   //zp
      else if (l == 2) {zn_ = 1.099367;}                   //zd
    }
    else if (atmnr == 24) {                            //Cr
      if (l == 0) {zn_ = 1.619853;}                        //zs
      else if (l == 1) {zn_ = 0.848266;}                   //zp
      else if (l == 2) {zn_ = 1.405015;}                   //zd
    }
    else if (atmnr == 25) {                            //Mn
      if (l == 0) {zn_ = 1.132450;}                        //zs
      else if (l == 1) {zn_ = 1.390740;}                   //zp
      else if (l == 2) {zn_ = 0.962550;}                   //zd
    }
    else if (atmnr == 26) {                            //Fe
      if (l == 0) {zn_ = 1.459152;}                        //zs
      else if (l == 1) {zn_ = 1.392614;}                   //zp
      else if (l == 2) {zn_ = 2.161909;}                   //zd
    }
    else if (atmnr == 27) {                            //Co
      if (l == 0) {zn_ = 0.519518;}                        //zs
      else if (l == 1) {zn_ = 1.000000;}                   //zp
      else if (l == 2) {zn_ = 0.352115;}                   //zd
    }
    else if (atmnr == 28) {                            //Ni
      if (l == 0) {zn_ = 0.746470;}                        //zs
      else if (l == 1) {zn_ = 0.753327;}                   //zp
      else if (l == 2) {zn_ = 1.461345;}                   //zd
    }
    else if (atmnr == 29) {                            //Cu
      if (l == 0) {zn_ = 1.899598;}                        //zs
      else if (l == 1) {zn_ = 3.000000;}                   //zp
      else if (l == 2) {zn_ = 1.484317;}                   //zd
    }
    else if (atmnr == 33) {                            //As
      if (l == 0) {zn_ = 2.006543;}                        //zs
      else if (l == 1) {zn_ = 3.316832;}                   //zp
      else if (l == 2) {zn_ = 4.653440;}                   //zd
    }
    else if (atmnr == 35) {                            //Br
      if (l == 0) {zn_ = 3.094777;}                        //zs
      else if (l == 1) {zn_ = 3.065764;}                   //zp
      else if (l == 2) {zn_ = 2.820003;}                   //zd
    }
    else if (atmnr == 39) {                            //Y
      if (l == 0) {zn_ = 0.902611;}                        //zs
      else if (l == 1) {zn_ = 1.484400;}                   //zp
      else if (l == 2) {zn_ = 1.384238;}                   //zd
    }
    else if (atmnr == 40) {                            //Zr
      if (l == 0) {zn_ = 1.189109;}                        //zs
      else if (l == 1) {zn_ = 0.809092;}                   //zp
      else if (l == 2) {zn_ = 1.190249;}                   //zd
    }
    else if (atmnr == 41) {                            //Nb
      if (l == 0) {zn_ = 1.490754;}                        //zs
      else if (l == 1) {zn_ = 0.892760;}                   //zp
      else if (l == 2) {zn_ = 1.443837;}                   //zd
    }
    else if (atmnr == 42) {                            //Mo
      if (l == 0) {zn_ = 1.912995;}                        //zs
      else if (l == 1) {zn_ = 1.355055;}                   //zp
      else if (l == 2) {zn_ = 1.876231;}                   //zd
    }
    else if (atmnr == 43) {                            //Tc
      if (l == 0) {zn_ = 1.411033;}                        //zs
      else if (l == 1) {zn_ = 1.141313;}                   //zp
      else if (l == 2) {zn_ = 1.159312;}                   //zd
    }
    else if (atmnr == 44) {                            //Ru
      if (l == 0) {zn_ = 0.984449;}                        //zs
      else if (l == 1) {zn_ = 4.586613;}                   //zp
      else if (l == 2) {zn_ = 0.765332;}                   //zd
    }
    else if (atmnr == 45) {                            //Rh
      if (l == 0) {zn_ = 0.809923;}                        //zs
      else if (l == 1) {zn_ = 6.898259;}                   //zp
      else if (l == 2) {zn_ = 0.643134;}                   //zd
    }
    else if (atmnr == 46) {                            //Pd
      if (l == 0) {zn_ = 1.794085;}                        //zs
      else if (l == 1) {zn_ = 6.158778;}                   //zp
      else if (l == 2) {zn_ = 1.630913;}                   //zd
    }
    else if (atmnr == 47) {                            //Ag
      if (l == 0) {zn_ = 0.695514;}                        //zs
      else if (l == 1) {zn_ = 4.729949;}                   //zp
      else if (l == 2) {zn_ = 0.506522;}                   //zd
    }
    else if (atmnr == 51) {                            //Sb
      if (l == 0) {zn_ = 5.993591;}                        //zs
      else if (l == 1) {zn_ = 6.145086;}                   //zp
      else if (l == 2) {zn_ = 5.704031;}                   //zd
    }
    else if (atmnr == 53) {                            //I
      if (l == 0) {zn_ = 9.135244;}                        //zs
      else if (l == 1) {zn_ = 6.888191;}                   //zp
      else if (l == 2) {zn_ = 3.791523;}                   //zd
    }
    else if (atmnr == 57) {                            //La
      if (l == 0) {zn_ = 1.617784;}                        //zs
      else if (l == 1) {zn_ = 4.331620;}                   //zp
      else if (l == 2) {zn_ = 2.285738;}                   //zd
    }
    else if (atmnr == 71) {                            //Lu
      if (l == 0) {zn_ = 1.632335;}                        //zs
      else if (l == 1) {zn_ = 4.033128;}                   //zp
      else if (l == 2) {zn_ = 0.921999;}                   //zd
    }
    else if (atmnr == 72) {                            //Hf
      if (l == 0) {zn_ = 0.946927;}                        //zs
      else if (l == 1) {zn_ = 3.538911;}                   //zp
      else if (l == 2) {zn_ = 0.940283;}                   //zd
    }
    else if (atmnr == 73) {                            //Ta
      if (l == 0) {zn_ = 1.741367;}                        //zs
      else if (l == 1) {zn_ = 3.430157;}                   //zp
      else if (l == 2) {zn_ = 2.311198;}                   //zd
    }
    else if (atmnr == 74) {                            //W
      if (l == 0) {zn_ = 1.498860;}                        //zs
      else if (l == 1) {zn_ = 1.965900;}                   //zp
      else if (l == 2) {zn_ = 1.876450;}                   //zd
    }
    else if (atmnr == 75) {                            //Re
      if (l == 0) {zn_ = 1.680823;}                        //zs
      else if (l == 1) {zn_ = 1.331218;}                   //zp
      else if (l == 2) {zn_ = 1.490623;}                   //zd
    }
    else if (atmnr == 76) {                            //Os
      if (l == 0) {zn_ = 1.844700;}                        //zs
      else if (l == 1) {zn_ = 1.564220;}                   //zp
      else if (l == 2) {zn_ = 1.770010;}                   //zd
    }
    else if (atmnr == 77) {                            //Ir
      if (l == 0) {zn_ = 0.927246;}                        //zs
      else if (l == 1) {zn_ = 3.191892;}                   //zp
      else if (l == 2) {zn_ = 0.662007;}                   //zd
    }
    else if (atmnr == 78) {                            //Pt
      if (l == 0) {zn_ = 2.270699;}                        //zs
      else if (l == 1) {zn_ = 1.949896;}                   //zp
      else if (l == 2) {zn_ = 1.713856;}                   //zd
    }
    else if (atmnr == 79) {                            //Au
      if (l == 0) {zn_ = 2.444680;}                        //zs
      else if (l == 1) {zn_ = 7.014990;}                   //zp
      else if (l == 2) {zn_ = 1.777089;}                   //zd
    }
    return zn_;
  }
  double F0sd(size_t atmnr) {
    //function returning the F0sd 1-center eris; values stored in eV, returned in a.u.
    double f0sd = 0.0;
    if (atmnr == 21) {f0sd = 4.798313;}               //Sc
    else if (atmnr == 22) {f0sd = 6.560562;}          //Ti
    else if (atmnr == 23) {f0sd = 6.810021;}          //V
    else if (atmnr == 24) {f0sd = 6.150136;}          //Cr
    else if (atmnr == 25) {f0sd = 7.690920;}          //Mn
    else if (atmnr == 26) {f0sd = 9.300165;}          //Fe
    else if (atmnr == 27) {f0sd = 1.446283;}          //Co
    else if (atmnr == 28) {f0sd = 4.651664;}          //Ni
    else if (atmnr == 29) {f0sd = 9.848807;}          //Cu
    else if (atmnr == 39) {f0sd = 4.972716;}          //Y
    else if (atmnr == 40) {f0sd = 5.010704;}          //Zr
    else if (atmnr == 41) {f0sd = 6.550674;}          //Nb
    else if (atmnr == 42) {f0sd = 10.000608;}         //Mo
    else if (atmnr == 43) {f0sd = 5.434886;}          //Tc
    else if (atmnr == 44) {f0sd = 5.917404;}          //Ru
    else if (atmnr == 45) {f0sd = 1.775497;}          //Rh
    else if (atmnr == 46) {f0sd = 8.004447;}          //Pd
    else if (atmnr == 47) {f0sd = 1.938327;}          //Ag
    else if (atmnr == 57) {f0sd = 8.856858;}          //La
    else if (atmnr == 71) {f0sd = 3.924927;}          //Lu
    else if (atmnr == 72) {f0sd = 4.842900;}          //Hf
    else if (atmnr == 73) {f0sd = 8.544427;}          //Ta
    else if (atmnr == 74) {f0sd = 7.788180;}          //W
    else if (atmnr == 75) {f0sd = 5.442818;}          //Re
    else if (atmnr == 76) {f0sd = 2.021170;}          //Os
    else if (atmnr == 77) {f0sd = 2.627170;}          //Ir
    else if (atmnr == 78) {f0sd = 7.098591;}          //Pt
    else if (atmnr == 79) {f0sd = 8.827257;}          //Au
    return f0sd/au2eV;
  }
  double G2sd(size_t atmnr) {
    //function returning the G2sd 1-center eris; values stored in eV, returned in a.u.
    double g2sd = 0.0;
    if (atmnr == 21) {g2sd = 5.380136;}                //Sc
    else if (atmnr == 22) {g2sd = 3.396235;}           //Ti
    else if (atmnr == 23) {g2sd = 1.831407;}           //V
    else if (atmnr == 24) {g2sd = 2.000300;}           //Cr
    else if (atmnr == 25) {g2sd = 1.105330;}           //Mn
    else if (atmnr == 26) {g2sd = 1.601345;}           //Fe
    else if (atmnr == 27) {g2sd = 1.680225;}           //Co
    else if (atmnr == 28) {g2sd = 1.880502;}           //Ni
    else if (atmnr == 29) {g2sd = 9.847577;}           //Cu
    else if (atmnr == 39) {g2sd = 5.016364;}           //Y
    else if (atmnr == 40) {g2sd = 2.943652;}           //Zr
    else if (atmnr == 41) {g2sd = 1.065577;}           //Nb
    else if (atmnr == 42) {g2sd = 1.216752;}           //Mo
    else if (atmnr == 43) {g2sd = 1.106875;}           //Tc
    else if (atmnr == 44) {g2sd = 5.859738;}           //Ru
    else if (atmnr == 45) {g2sd = 1.851571;}           //Rh
    else if (atmnr == 46) {g2sd = 2.613148;}           //Pd
    else if (atmnr == 47) {g2sd = 1.071901;}           //Ag
    else if (atmnr == 57) {g2sd = 7.925585;}           //La
    else if (atmnr == 71) {g2sd = 1.000946;}           //Lu
    else if (atmnr == 72) {g2sd = 4.386101;}           //Hf
    else if (atmnr == 73) {g2sd = 2.074254;}           //Ta
    else if (atmnr == 74) {g2sd = 1.684940;}           //W
    else if (atmnr == 75) {g2sd = 2.376279;}           //Re
    else if (atmnr == 76) {g2sd = 1.392130;}           //Os
    else if (atmnr == 77) {g2sd = 2.996029;}           //Ir
    else if (atmnr == 78) {g2sd = 4.484183;}           //Pt
    else if (atmnr == 79) {g2sd = 4.915625;}           //Au
    return g2sd/au2eV;
  }
  double eri1Center(int atmnr, int Lbra, int Lket) {
    //function that gives back the semi-empirical 1-center eris
    //values stored in eV but returned in a.u.
    //Lbra is the sum of azimuthal quantum numbers for bra (ss = 0; pp = 2; sp = 1)
    //Lket is the sum of azimuthal quantum numbers for ket (ss = 0; pp = 2; sp = 1; p*p* = -2)
    double eri = 0.0;
    switch (atmnr) {
      case 1:       //H
      eri = 14.448686;
      break;
    case 2:       //He
      if ((Lbra == 0)&&(Lket == 0)) {eri = 9.445299;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 9.214548;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 11.201419;}     //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 13.046115;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.299954;}                                    //(sp|sp)||(ps|ps)
      break;
    case 3:       //Li
      if ((Lbra == 0)&&(Lket == 0)) {eri = 11.035907;}                                        //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.543650;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 19.998647;}     //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.059036;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.641886;}                                    //(sp|sp)||(ps|ps)
      break;
    case 4:       //Be
      if ((Lbra == 0)&&(Lket == 0)) {eri = 7.552804;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 12.862153;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 10.203146;}     //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 13.602858;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.501452;}                                    //(sp|sp)||(ps|ps)
      break;
    case 5:       //B
      if ((Lbra == 0)&&(Lket == 0)) {eri = 8.179341;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.829395;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.294021;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.401072;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.252845;}                                    //(sp|sp)||(ps|ps)
      break;
    case 6:       //C
      if ((Lbra == 0)&&(Lket == 0)) {eri = 13.335519;}                                        //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 10.778326;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 11.528134;}     //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.486212;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.717322;}                                    //(sp|sp)||(ps|ps)
      break;
    case 7:       //N
      if ((Lbra == 0)&&(Lket == 0)) {eri = 12.357026;}                                        //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 12.570756;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 9.636190;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 10.576425;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.871545;}                                    //(sp|sp)||(ps|ps)
      break;
    case 8:       //O
      if ((Lbra == 0)&&(Lket == 0)) {eri = 11.304042;}                                        //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 13.618205;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 15.807424;}     //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 10.332765;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 5.010801;}                                    //(sp|sp)||(ps|ps)
      break;
    case 9:       //F
      if ((Lbra == 0)&&(Lket == 0)) {eri = 12.446818;}                                        //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.417366;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 18.496082;}     //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 12.179816;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.604382;}                                    //(sp|sp)||(ps|ps)
      break;
    case 10:      //Ne
      if ((Lbra == 0)&&(Lket == 0)) {eri = 19.999574;}                                        //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.963560;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 16.896951;}     //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 16.027799;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.779280;}                                    //(sp|sp)||(ps|ps)
      break;
    case 11:      //Na
      if ((Lbra == 0)&&(Lket == 0)) {eri = 4.059972;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 9.283540;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.061183;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 17.034978;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.640715;}                                    //(sp|sp)||(ps|ps)
      break;
    case 12:      //Mg
      if ((Lbra == 0)&&(Lket == 0)) {eri = 7.115328;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 4.737311;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 3.253024;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 8.428485;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.877379;}                                    //(sp|sp)||(ps|ps)
      break;
    case 13:      //Al
      if ((Lbra == 0)&&(Lket == 0)) {eri = 6.652155;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.668857;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.459435;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.673299;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.435060;}                                    //(sp|sp)||(ps|ps)
      break;
    case 14:      //Si
      if ((Lbra == 0)&&(Lket == 0)) {eri = 5.194805;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.185150;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.090534;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.769775;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.425012;}                                    //(sp|sp)||(ps|ps)
      break;
    case 15:      //P
      if ((Lbra == 0)&&(Lket == 0)) {eri = 8.758856;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.662754;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 8.483679;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.734264;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.871681;}                                    //(sp|sp)||(ps|ps)
      break;
    case 16:      //S
      if ((Lbra == 0)&&(Lket == 0)) {eri = 9.170350;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.165473;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.944296;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.301878;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 5.005404;}                                    //(sp|sp)||(ps|ps)
      break;
    case 17:      //Cl
      if ((Lbra == 0)&&(Lket == 0)) {eri = 11.142654;}                                        //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 9.551886;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.487881;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 8.128436;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 5.004267;}                                    //(sp|sp)||(ps|ps)
      break;
    case 18:      //Ar
      if ((Lbra == 0)&&(Lket == 0)) {eri = 17.858776;}                                        //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 11.852500;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 4.168451;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 15.669543;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 4.574549;}                                    //(sp|sp)||(ps|ps)
      break;
    case 19:      //K
      if ((Lbra == 0)&&(Lket == 0)) {eri = 3.369251;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 0.999505;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.129351;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 18.999148;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.300325;}                                    //(sp|sp)||(ps|ps)
      break;
    case 20:      //Ca
      if ((Lbra == 0)&&(Lket == 0)) {eri = 5.725773;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.172103;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 4.781065;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.431876;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.240572;}                                    //(sp|sp)||(ps|ps)
      break;
    case 21:      //Sc
      if ((Lbra == 0)&&(Lket == 0)) {eri = 4.638215;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 14.604872;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.739164;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 12.802595;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.193835;}                                    //(sp|sp)||(ps|ps)
      break;
    case 22:      //Ti
      if ((Lbra == 0)&&(Lket == 0)) {eri = 5.717851;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.414726;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.800015;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.623133;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.403732;}                                    //(sp|sp)||(ps|ps)
      break;
    case 23:      //V
      if ((Lbra == 0)&&(Lket == 0)) {eri = 5.983116;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 4.499763;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 4.736769;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 3.944481;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.901105;}                                    //(sp|sp)||(ps|ps)
      break;
    case 24:      //Cr
      if ((Lbra == 0)&&(Lket == 0)) {eri = 8.855572;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.053094;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.588631;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.429530;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.648039;}                                    //(sp|sp)||(ps|ps)
      break;
    case 25:      //Mn
      if ((Lbra == 0)&&(Lket == 0)) {eri = 6.190990;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.284594;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.757427;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.262255;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.520518;}                                    //(sp|sp)||(ps|ps)
      break;
    case 26:      //Fe
      if ((Lbra == 0)&&(Lket == 0)) {eri = 7.977036;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.295758;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.786867;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.272041;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.880189;}                                    //(sp|sp)||(ps|ps)
      break;
    case 27:      //Co
      if ((Lbra == 0)&&(Lket == 0)) {eri = 2.840152;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.956968;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 3.425933;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.221864;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.390087;}                                    //(sp|sp)||(ps|ps)
      break;
    case 28:      //Ni
      if ((Lbra == 0)&&(Lket == 0)) {eri = 4.080876;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 4.487545;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 4.099452;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 3.933771;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.993498;}                                    //(sp|sp)||(ps|ps)
      break;
    case 29:      //Cu
      if ((Lbra == 0)&&(Lket == 0)) {eri = 10.384910;}                                        //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 17.870905;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 12.145361;}     //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 15.665592;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 2.037394;}                                    //(sp|sp)||(ps|ps)
      break;
    case 30:      //Zn
      if ((Lbra == 0)&&(Lket == 0)) {eri = 8.707424;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 20.000041;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 3.436116;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.782785;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.662036;}                                    //(sp|sp)||(ps|ps)
      break;
    case 31:      //Ga
      if ((Lbra == 0)&&(Lket == 0)) {eri = 10.354885;}                                        //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.090184;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.993674;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.299226;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.295974;}                                    //(sp|sp)||(ps|ps)
      break;
    case 32:      //Ge
      if ((Lbra == 0)&&(Lket == 0)) {eri = 7.518301;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.066801;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.594443;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.305947;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.290742;}                                    //(sp|sp)||(ps|ps)
      break;
    case 33:      //As
      if ((Lbra == 0)&&(Lket == 0)) {eri = 6.665030;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 9.310836;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.213867;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 8.712542;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.280662;}                                    //(sp|sp)||(ps|ps)
      break;
    case 34:      //Se
      if ((Lbra == 0)&&(Lket == 0)) {eri = 5.522356;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.042391;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 2.907562;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.735106;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.095789;}                                    //(sp|sp)||(ps|ps)
      break;
    case 35:      //Br
      if ((Lbra == 0)&&(Lket == 0)) {eri = 7.616791;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 9.649216;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.010425;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 8.343792;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 4.996553;}                                    //(sp|sp)||(ps|ps)
      break;
    case 36:      //Kr
      if ((Lbra == 0)&&(Lket == 0)) {eri = 19.999857;}                                        //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 9.174784;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 1.175304;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 14.926948;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.299867;}                                    //(sp|sp)||(ps|ps)
      break;
    case 37:      //Rb
      if ((Lbra == 0)&&(Lket == 0)) {eri = 6.680824;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.068874;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 20.001098;}     //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 2.747860;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.602834;}                                    //(sp|sp)||(ps|ps)
      break;
    case 38:      //Sr
      if ((Lbra == 0)&&(Lket == 0)) {eri = 4.603664;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.334620;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.716069;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.443088;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.831527;}                                    //(sp|sp)||(ps|ps)
      break;
    case 39:      //Y
      if ((Lbra == 0)&&(Lket == 0)) {eri = 4.046733;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.278752;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 4.726277;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.343281;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.679228;}                                    //(sp|sp)||(ps|ps)
      break;
    case 40:      //Zr
      if ((Lbra == 0)&&(Lket == 0)) {eri = 5.331208;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 3.967381;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 4.150579;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 3.457490;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.743676;}                                    //(sp|sp)||(ps|ps)
      break;
    case 41:      //Nb
      if ((Lbra == 0)&&(Lket == 0)) {eri = 6.683592;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 4.377647;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 4.685339;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 3.815028;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.650679;}                                    //(sp|sp)||(ps|ps)
      break;
    case 42:      //Mo
      if ((Lbra == 0)&&(Lket == 0)) {eri = 8.576652;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.644509;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.888293;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.790552;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.317368;}                                    //(sp|sp)||(ps|ps)
      break;
    case 43:      //Tc
      if ((Lbra == 0)&&(Lket == 0)) {eri = 6.326174;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.596426;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.587138;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.877169;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.258989;}                                    //(sp|sp)||(ps|ps)
      break;
    case 44:      //Ru
      if ((Lbra == 0)&&(Lket == 0)) {eri = 4.413643;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 22.490448;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.356996;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 19.599957;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.008058;}                                    //(sp|sp)||(ps|ps)
      break;
    case 45:      //Rh
      if ((Lbra == 0)&&(Lket == 0)) {eri = 3.631179;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 33.825599;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 4.407820;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 29.478305;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.000092;}                                    //(sp|sp)||(ps|ps)
      break;
    case 46:      //Pd
      if ((Lbra == 0)&&(Lket == 0)) {eri = 8.043535;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 30.199556;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 9.755042;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 26.318284;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.086121;}                                    //(sp|sp)||(ps|ps)
      break;
    case 47:      //Ag
      if ((Lbra == 0)&&(Lket == 0)) {eri = 3.118242;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 23.193295;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 3.785152;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 20.212474;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.000432;}                                    //(sp|sp)||(ps|ps)
      break;
    case 48:      //Cd
      if ((Lbra == 0)&&(Lket == 0)) {eri = 6.677284;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 18.729843;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.953373;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 9.917452;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.825192;}                                    //(sp|sp)||(ps|ps)
      break;
    case 49:      //In
      if ((Lbra == 0)&&(Lket == 0)) {eri = 9.906091;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 4.826006;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 10.520060;}     //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.906563;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 3.500299;}                                    //(sp|sp)||(ps|ps)
      break;
    case 50:      //Sn
      if ((Lbra == 0)&&(Lket == 0)) {eri = 8.269655;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.584874;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.013349;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.855159;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.531212;}                                    //(sp|sp)||(ps|ps)
      break;
    case 51:      //Sb
      if ((Lbra == 0)&&(Lket == 0)) {eri = 10.588832;}                                        //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 9.281609;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.310023;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 8.954081;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.779112;}                                    //(sp|sp)||(ps|ps)
      break;
    case 52:      //Te
      if ((Lbra == 0)&&(Lket == 0)) {eri = 7.030626;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 7.883479;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 12.601389;}     //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 6.973163;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 5.000826;}                                    //(sp|sp)||(ps|ps)
      break;
    case 53:      //I
      if ((Lbra == 0)&&(Lket == 0)) {eri = 7.234759;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 9.877466;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 9.154406;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 8.035916;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 5.004215;}                                    //(sp|sp)||(ps|ps)
      break;
    case 54:      //Xe
      if ((Lbra == 0)&&(Lket == 0)) {eri = 20.000252;}                                        //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 2.305787;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 4.175902;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.063220;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 4.418843;}                                    //(sp|sp)||(ps|ps)
      break;
    case 55:      //Cs
      if ((Lbra == 0)&&(Lket == 0)) {eri = 6.464751;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 13.775390;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 4.004501;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 12.912537;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.026928;}                                    //(sp|sp)||(ps|ps)
      break;
    case 56:      //Ba
      if ((Lbra == 0)&&(Lket == 0)) {eri = 3.600823;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 3.345166;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 4.740579;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 3.142783;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.929429;}                                    //(sp|sp)||(ps|ps)
      break;
    case 57:      //La
      if ((Lbra == 0)&&(Lket == 0)) {eri = 6.154440;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 18.077465;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.322704;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 15.679057;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.138601;}                                    //(sp|sp)||(ps|ps)
      break;
    case 71:      //Lu
      if ((Lbra == 0)&&(Lket == 0)) {eri = 6.209796;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 16.831746;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.379102;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 14.598613;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.209008;}                                    //(sp|sp)||(ps|ps)
      break;
    case 72:      //Hf
      if ((Lbra == 0)&&(Lket == 0)) {eri = 3.602338;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 14.769194;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 4.293729;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 12.809708;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.011028;}                                    //(sp|sp)||(ps|ps)
      break;
    case 73:      //Ta
      if ((Lbra == 0)&&(Lket == 0)) {eri = 6.624580;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 14.315323;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.805321;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 12.416054;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.577263;}                                    //(sp|sp)||(ps|ps)
      break;
    case 74:      //W
      if ((Lbra == 0)&&(Lket == 0)) {eri = 5.702025;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.204433;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.323145;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.115919;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.319912;}                                    //(sp|sp)||(ps|ps)
      break;
    case 75:      //Re
      if ((Lbra == 0)&&(Lket == 0)) {eri = 6.394256;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 5.555669;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 5.555571;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 4.818577;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.220913;}                                    //(sp|sp)||(ps|ps)
      break;
    case 76:      //Os
      if ((Lbra == 0)&&(Lket == 0)) {eri = 7.017683;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.528073;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.384200;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.661968;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.508926;}                                    //(sp|sp)||(ps|ps)
      break;
    case 77:      //Ir
      if ((Lbra == 0)&&(Lket == 0)) {eri = 3.527467;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 13.320955;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 4.203820;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 11.553612;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.018501;}                                    //(sp|sp)||(ps|ps)
      break;
    case 78:      //Pt
      if ((Lbra == 0)&&(Lket == 0)) {eri = 8.638286;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.137643;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.922254;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.057990;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.892617;}                                    //(sp|sp)||(ps|ps)
      break;
    case 79:      //Au
      if ((Lbra == 0)&&(Lket == 0)) {eri = 9.300152;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 29.276168;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 11.073443;}     //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 25.391984;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.144384;}                                    //(sp|sp)||(ps|ps)
      break;
    case 80:      //Hg
      if ((Lbra == 0)&&(Lket == 0)) {eri = 6.372822;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 10.397393;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 10.143176;}     //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 14.794056;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.926128;}                                    //(sp|sp)||(ps|ps)
      break;
    case 81:      //Tl
      if ((Lbra == 0)&&(Lket == 0)) {eri = 5.015118;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 10.495551;}                                   //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 13.932049;}     //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 10.526198;}   //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.293760;}                                    //(sp|sp)||(ps|ps)
      break;
    case 82:      //Pb
      if ((Lbra == 0)&&(Lket == 0)) {eri = 5.254128;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 6.818551;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 7.061016;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 5.603019;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 1.018819;}                                    //(sp|sp)||(ps|ps)
      break;
    case 83:      //Bi
      if ((Lbra == 0)&&(Lket == 0)) {eri = 5.851803;}                                         //(ss|ss)
      else if ((Lbra == 2)&&(Lket == 2)) {eri = 8.389442;}                                    //(pp|pp)
      else if (((Lbra == 0)&&(Lket == 2))||((Lbra == 2)&&(Lket == 0))) {eri = 6.790583;}      //(ss|pp)||(pp|ss)
      else if (((Lbra == 2)&&(Lket == -2))||((Lbra == -2)&&(Lket == 2))) {eri = 7.724219;}    //(pp|p*p*)||(p*p*|pp)
      else if ((Lbra == 1)&&(Lket == 1)) {eri = 0.295606;}                                    //(sp|sp)||(ps|ps)
      break;
    }
    return eri/au2eV;
  }
};

#endif //_MNDOD_
