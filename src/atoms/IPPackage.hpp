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

#ifndef _Ionization_Potential_Package_
#define _Ionization_Potential_Package_
#include <vector>
#include <string>

//descrition:
//base function containing ionization potentials
/*TODO:
- transform the if-else into switch, for more efficiency
*/

std::vector<double> IonizationPot(size_t atom, size_t IPtype = 0, std::string hardness = "soft") {
  //function that returns ionization potentials according to different methods
  //atom is the atomic number
  std::vector<double> ip(1,0.0);
  //orbital IPs are in atomic units
  if (IPtype == 0) {                        //Slater's values calculated in SlaterExponents.xlsx
    if (atom == 1) {ip[0] = -0.500;}     //H
    else if (atom == 2) {ip[0] = -1.445;}//He
    else if (atom == 3) {                //Li
      ip[0] = -3.645;                      //1s
      ip.push_back(-0.211);                //2s
      ip.push_back(-0.211);                //2p
    }
    else if (atom == 4) {                //Be
      ip[0] = -6.845;                      //1s
      ip.push_back(-0.475);                //2s
      ip.push_back(-0.475);                //2p
    }
    else if (atom == 5) {                //B
      ip[0] = -11.045;                     //1s
      ip.push_back(-0.845);                //2s
      ip.push_back(-0.845);                //2p
    }
    else if (atom == 6) {                //C
      ip[0] = -16.245;                     //1s
      ip.push_back(-1.320);                //2s
      ip.push_back(-1.320);                //2p
    }
    else if (atom == 7) {                //N
      ip[0] = -22.445;                     //1s
      ip.push_back(-1.901);                //2s
      ip.push_back(-1.901);                //2p
    }
    else if (atom == 8) {                //O
      ip[0] = -29.645;                     //1s
      ip.push_back(-2.588);                //2s
      ip.push_back(-2.588);                //2p
    }
    else if (atom == 9) {                //F
      ip[0] = -37.845;                     //1s
      ip.push_back(-3.380);                //2s
      ip.push_back(-3.380);                //2p
    }
    else if (atom == 10) {               //Ne
      ip[0] = -47.045;                     //1s
      ip.push_back(-4.278);                //2s
      ip.push_back(-4.278);                //2p
    }
    else if (atom == 11) {               //Na
      ip[0] = -57.245;                     //1s
      ip.push_back(-5.865);                //2s
      ip.push_back(-5.865);                //2p
      ip.push_back(-0.269);                //3s
      ip.push_back(-0.269);                //3p
    }
    else if (atom == 12) {               //Mg
      ip[0] = -68.445;                     //1s
      ip.push_back(-7.703);                //2s
      ip.push_back(-7.703);                //2p
      ip.push_back(-0.451);                //3s
      ip.push_back(-0.451);                //3p
    }
    else if (atom == 13) {               //Al
      ip[0] = -80.645;                     //1s
      ip.push_back(-9.790);                //2s
      ip.push_back(-9.790);                //2p
      ip.push_back(-0.681);                //3s
      ip.push_back(-0.681);                //3p
    }
    else if (atom == 14) {               //Si
      ip[0] = -93.845;                     //1s
      ip.push_back(-12.128);               //2s
      ip.push_back(-12.128);               //2p
      ip.push_back(-0.957);                //3s
      ip.push_back(-0.957);                //3p
    }
    else if (atom == 15) {               //P
      ip[0] = -108.045;                    //1s
      ip.push_back(-14.715);               //2s
      ip.push_back(-14.715);               //2p
      ip.push_back(-1.280);                //3s
      ip.push_back(-1.280);                //3p
    }
    else if (atom == 16) {               //S
      ip[0] = -123.245;                    //1s
      ip.push_back(-17.553);               //2s
      ip.push_back(-17.553);               //2p
      ip.push_back(-1.650);                //3s
      ip.push_back(-1.650);                //3p
    }
    else if (atom == 17) {               //Cl
      ip[0] = -139.445;                    //1s
      ip.push_back(-20.640);               //2s
      ip.push_back(-20.640);               //2p
      ip.push_back(-2.067);                //3s
      ip.push_back(-2.067);                //3p
    }
    else if (atom == 18) {              //Ar
      ip[0] = -156.645;                    //1s
      ip.push_back(-23.978);               //2s
      ip.push_back(-23.978);               //2p
      ip.push_back(-2.531);                //3s
      ip.push_back(-2.531);                //3p
    }		
    else if (atom == 19) {               //K
      ip[0] = -174.845;                     //1s
      ip.push_back(-27.565);                //2s
      ip.push_back(-27.565);                //2p
      ip.push_back(-3.337);                 //3s
      ip.push_back(-3.337);                 //3p
      ip.push_back(0.00);                   //3d
      ip.push_back(-0.177);                 //4s
      ip.push_back(-0.177);                 //4p
    }
    else if (atom == 20) {               //Ca
      ip[0] = -194.045;                     //1s
      ip.push_back(-31.403);                //2s
      ip.push_back(-31.403);                //2p
      ip.push_back(-4.253);                 //3s
      ip.push_back(-4.253);                 //3p
      ip.push_back(0.00);                   //3d
      ip.push_back(-0.297);                 //4s
      ip.push_back(-0.297);                 //4p
    }
    else if (atom == 21) {               //Sc
      ip[0] = -214.245;                     //1s
      ip.push_back(-35.490);                //2s
      ip.push_back(-35.490);                //2p
      ip.push_back(-5.281);                 //3s
      ip.push_back(-5.281);                 //3p
      ip.push_back(-0.500);                 //3d
      ip.push_back(-0.329);                 //4s
      ip.push_back(-0.329);                 //4p
    }
    else if (atom == 22) {               //Ti
      ip[0] = -235.445;                     //1s
      ip.push_back(-39.828);                //2s
      ip.push_back(-39.828);                //2p
      ip.push_back(-6.420);                 //3s
      ip.push_back(-6.420);                 //3p
      ip.push_back(-0.740);                 //3d
      ip.push_back(-0.362);                 //4s
      ip.push_back(-0.362);                 //4p
    }
    else if (atom == 23) {               //V
      ip[0] = -257.645;                     //1s
      ip.push_back(-44.415);                //2s
      ip.push_back(-44.415);                //2p
      ip.push_back(-7.670);                 //3s
      ip.push_back(-7.670);                 //3p
      ip.push_back(-1.027);                 //3d
      ip.push_back(-0.398);                 //4s
      ip.push_back(-0.398);                 //4p
    }
    else if (atom == 24) {               //Cr
      ip[0] = -280.845;                     //1s
      ip.push_back(-49.253);                //2s
      ip.push_back(-49.253);                //2p
      ip.push_back(-9.031);                 //3s
      ip.push_back(-9.031);                 //3p
      ip.push_back(-1.361);                 //3d
      ip.push_back(-0.435);                 //4s
      ip.push_back(-0.435);                 //4p
    }
    else if (atom == 25) {               //Mn
      ip[0] = -305.045;                     //1s
      ip.push_back(-54.340);                //2s
      ip.push_back(-54.340);                //2p
      ip.push_back(-10.503);                //3s
      ip.push_back(-10.503);                //3p
      ip.push_back(-1.742);                 //3d
      ip.push_back(-0.473);                 //4s
      ip.push_back(-0.473);                 //4p
    }
    else if (atom == 26) {              //Fe
      ip[0] = -330.245;                     //1s
      ip.push_back(-59.678);                //2s
      ip.push_back(-59.678);                //2p
      ip.push_back(-12.087);                //3s
      ip.push_back(-12.087);                //3p
      ip.push_back(-2.170);                 //3d
      ip.push_back(-0.514);                 //4s
      ip.push_back(-0.514);                 //4p
    }
    else if (atom == 27) {               //Co
      ip[0] = -356.445;                     //1s
      ip.push_back(-65.265);                //2s
      ip.push_back(-65.265);                //2p
      ip.push_back(-13.781);                //3s
      ip.push_back(-13.781);                //3p
      ip.push_back(-2.645);                 //3d
      ip.push_back(-0.556);                 //4s
      ip.push_back(-0.556);                 //4p
    }
    else if (atom == 28) {               //Ni
      ip[0] = -383.645;                     //1s
      ip.push_back(-71.103);                //2s
      ip.push_back(-71.103);                //2p
      ip.push_back(-15.587);                //3s
      ip.push_back(-15.587);                //3p
      ip.push_back(-3.167);                 //3d
      ip.push_back(-0.599);                 //4s
      ip.push_back(-0.599);                 //4p
    }
    else if (atom == 29) {               //Cu
      ip[0] = -411.845;                     //1s
      ip.push_back(-77.190);                //2s
      ip.push_back(-77.190);                //2p
      ip.push_back(-17.503);                //3s
      ip.push_back(-17.503);                //3p
      ip.push_back(-3.736);                 //3d
      ip.push_back(-0.644);                 //4s
      ip.push_back(-0.644);                 //4p
    }
    else if (atom == 30) {               //Zn
      ip[0] = -441.045;                     //1s
      ip.push_back(-83.528);                //2s
      ip.push_back(-83.528);                //2p
      ip.push_back(-19.531);                //3s
      ip.push_back(-19.531);                //3p
      ip.push_back(-4.351);                 //3d
      ip.push_back(-0.691);                 //4s
      ip.push_back(-0.691);                 //4p
    }
    else if (atom == 31) {               //Ga
      ip[0] = -471.245;                     //1s
      ip.push_back(-90.115);                //2s
      ip.push_back(-90.115);                //2p
      ip.push_back(-21.670);                //3s
      ip.push_back(-21.670);                //3p
      ip.push_back(-5.390);                 //3d
      ip.push_back(-0.913);                 //4s
      ip.push_back(-0.913);                 //4p
    }
    else if (atom == 32) {               //Ge
      ip[0] = -502.445;                     //1s
      ip.push_back(-96.953);                //2s
      ip.push_back(-96.953);                //2p
      ip.push_back(-23.920);                //3s
      ip.push_back(-23.920);                //3p
      ip.push_back(-6.540);                 //3d
      ip.push_back(-1.166);                 //4s
      ip.push_back(-1.166);                 //4p
    }
    else if (atom == 33) {               //As
      ip[0] = -534.645;                     //1s
      ip.push_back(-104.040);               //2s
      ip.push_back(-104.040);               //2p
      ip.push_back(-26.281);                //3s
      ip.push_back(-26.281);                //3p
      ip.push_back(-7.801);                 //3d
      ip.push_back(-1.450);                 //4s
      ip.push_back(-1.450);                 //4p
    }
    else if (atom == 34) {              //Se
      ip[0] = -567.845;                     //1s
      ip.push_back(-111.378);               //2s
      ip.push_back(-111.378);               //2p
      ip.push_back(-28.753);                //3s
      ip.push_back(-28.753);                //3p
      ip.push_back(-9.173);                 //3d
      ip.push_back(-1.764);                 //4s
      ip.push_back(-1.764);                 //4p
    }
    else if (atom == 35) {               //Br
      ip[0] = -602.045;                     //1s
      ip.push_back(-118.965);               //2s
      ip.push_back(-118.965);               //2p
      ip.push_back(-31.337);                //3s
      ip.push_back(-31.337);                //3p
      ip.push_back(-10.657);                //3d
      ip.push_back(-2.110);                 //4s
      ip.push_back(-2.110);                 //4p
    }
    else if (atom == 36) {              //Kr
      ip[0] = -637.245;                     //1s
      ip.push_back(-126.803);               //2s
      ip.push_back(-126.803);               //2p
      ip.push_back(-34.031);                //3s
      ip.push_back(-34.031);                //3p
      ip.push_back(-12.251);                //3d
      ip.push_back(-2.486);                 //4s
      ip.push_back(-2.486);                 //4p
    }
    else if (atom == 37) {               //Rb
      ip[0] = -673.445;                     //1s
      ip.push_back(-134.890);               //2s
      ip.push_back(-134.890);               //2p
      ip.push_back(-36.837);                //3s
      ip.push_back(-36.837);                //3p
      ip.push_back(-13.957);                //3d
      ip.push_back(-3.125);                 //4s
      ip.push_back(-3.125);                 //4p
      ip.push_back(0.00);                   //4d
      ip.push_back(0.00);                   //4f
      ip.push_back(-0.151);                 //5s
      ip.push_back(-0.151);                 //5p
    }
    else if (atom == 38) {               //Sr
      ip[0] = -710.645;                     //1s
      ip.push_back(-143.228);               //2s
      ip.push_back(-143.228);               //2p
      ip.push_back(-39.753);                //3s
      ip.push_back(-39.753);                //3p
      ip.push_back(-15.773);                //3d
      ip.push_back(-3.837);                 //4s
      ip.push_back(-3.837);                 //4p
      ip.push_back(0.00);                   //4d
      ip.push_back(0.00);                   //4f
      ip.push_back(-0.254);                 //5s
      ip.push_back(-0.254);                 //5p
    }
    else if (atom == 39) {               //Y
      ip[0] = -748.845;                     //1s
      ip.push_back(-151.815);               //2s
      ip.push_back(-151.815);               //2p
      ip.push_back(-42.781);                //3s
      ip.push_back(-42.781);                //3p
      ip.push_back(-17.701);                //3d
      ip.push_back(-4.622);                 //4s
      ip.push_back(-4.622);                 //4p
      ip.push_back(-0.329);                 //4d
      ip.push_back(0.00);                   //4f
      ip.push_back(-0.281);                 //5s
      ip.push_back(-0.281);                 //5p
    }
    else if (atom == 40) {               //Zr
      ip[0] = -788.045;                     //1s
      ip.push_back(-160.653);               //2s
      ip.push_back(-160.653);               //2p
      ip.push_back(-45.920);                //3s
      ip.push_back(-45.920);                //3p
      ip.push_back(-19.740);                //3d
      ip.push_back(-5.481);                 //4s
      ip.push_back(-5.481);                 //4p
      ip.push_back(-0.487);                 //4d
      ip.push_back(0.00);                   //4f
      ip.push_back(-0.310);                 //5s
      ip.push_back(-0.310);                 //5p
    }
    else if (atom == 41) {               //Nb
      ip[0] = -828.245;                     //1s
      ip.push_back(-169.740);               //2s
      ip.push_back(-169.740);               //2p
      ip.push_back(-49.170);                //3s
      ip.push_back(-49.170);                //3p
      ip.push_back(-21.890);                //3d
      ip.push_back(-6.412);                 //4s
      ip.push_back(-6.412);                 //4p
      ip.push_back(-0.675);                 //4d
      ip.push_back(0.00);                   //4f
      ip.push_back(-0.340);                 //5s
      ip.push_back(-0.340);                 //5p
    }
    else if (atom == 42) {               //Mo
      ip[0] = -869.445;                     //1s
      ip.push_back(-179.078);               //2s
      ip.push_back(-179.078);               //2p
      ip.push_back(-52.531);                //3s
      ip.push_back(-52.531);                //3p
      ip.push_back(-24.151);                //3d
      ip.push_back(-7.416);                 //4s
      ip.push_back(-7.416);                 //4p
      ip.push_back(-0.895);                 //4d
      ip.push_back(0.00);                   //4f
      ip.push_back(-0.372);                 //5s
      ip.push_back(-0.372);                 //5p
    }
    else if (atom == 43) {               //Tc
      ip[0] = -911.645;                     //1s
      ip.push_back(-188.665);               //2s
      ip.push_back(-188.665);               //2p
      ip.push_back(-56.003);                //3s
      ip.push_back(-56.003);                //3p
      ip.push_back(-26.523);                //3d
      ip.push_back(-8.494);                 //4s
      ip.push_back(-8.494);                 //4p
      ip.push_back(-1.145);                 //4d
      ip.push_back(0.00);                   //4f
      ip.push_back(-0.405);                 //5s
      ip.push_back(-0.405);                 //5p
    }
    else if (atom == 44) {              //Ru
      ip[0] = -954.845;                     //1s
      ip.push_back(-198.503);               //2s
      ip.push_back(-198.503);               //2p
      ip.push_back(-59.587);                //3s
      ip.push_back(-59.587);                //3p
      ip.push_back(-29.007);                //3d
      ip.push_back(-9.644);                 //4s
      ip.push_back(-9.644);                 //4p
      ip.push_back(-1.427);                 //4d
      ip.push_back(0.00);                   //4f
      ip.push_back(-0.439);                 //5s
      ip.push_back(-0.439);                 //5p
    }
    else if (atom == 45) {               //Rh
      ip[0] = -999.045;                     //1s
      ip.push_back(-208.590);               //2s
      ip.push_back(-208.590);               //2p
      ip.push_back(-63.281);                //3s
      ip.push_back(-63.281);                //3p
      ip.push_back(-31.601);                //3d
      ip.push_back(-10.868);                //4s
      ip.push_back(-10.868);                //4p
      ip.push_back(-1.739);                 //4d
      ip.push_back(0.00);                   //4f
      ip.push_back(-0.475);                 //5s
      ip.push_back(-0.475);                 //5p
    }
    else if (atom == 46) {               //Pd
      ip[0] = -1044.245;                    //1s
      ip.push_back(-218.928);               //2s
      ip.push_back(-218.928);               //2p
      ip.push_back(-67.087);                //3s
      ip.push_back(-67.087);                //3p
      ip.push_back(-34.307);                //3d
      ip.push_back(-12.164);                //4s
      ip.push_back(-12.164);                //4p
      ip.push_back(-2.082);                 //4d
      ip.push_back(0.00);                   //4f
      ip.push_back(-0.513);                 //5s
      ip.push_back(-0.513);                 //5p
    }
    else if (atom == 47) {               //Ag
      ip[0] = -1090.445;                    //1s
      ip.push_back(-229.515);               //2s
      ip.push_back(-229.515);               //2p
      ip.push_back(-71.003);                //3s
      ip.push_back(-71.003);                //3p
      ip.push_back(-37.123);                //3d
      ip.push_back(-13.534);                //4s
      ip.push_back(-13.534);                //4p
      ip.push_back(-2.456);                 //4d
      ip.push_back(0.00);                   //4f
      ip.push_back(-0.551);                 //5s
      ip.push_back(-0.551);                 //5p
    }
    else if (atom == 48) {               //Cd
      ip[0] = -1137.645;                    //1s
      ip.push_back(-240.353);               //2s
      ip.push_back(-240.353);               //2p
      ip.push_back(-75.031);                //3s
      ip.push_back(-75.031);                //3p
      ip.push_back(-40.051);                //3d
      ip.push_back(-14.977);                //4s
      ip.push_back(-14.977);                //4p
      ip.push_back(-2.861);                 //4d
      ip.push_back(0.00);                   //4f
      ip.push_back(-0.591);                 //5s
      ip.push_back(-0.591);                 //5p
    }
    else if (atom == 49) {               //In
      ip[0] = -1185.845;                    //1s
      ip.push_back(-251.440);               //2s
      ip.push_back(-251.440);               //2p
      ip.push_back(-79.170);                //3s
      ip.push_back(-79.170);                //3p
      ip.push_back(-43.090);                //3d
      ip.push_back(-16.492);                //4s
      ip.push_back(-16.492);                //4p
      ip.push_back(-3.544);                 //4d
      ip.push_back(0.00);                   //4f
      ip.push_back(-0.781);                 //5s
      ip.push_back(-0.781);                 //5p
    }
    else if (atom == 50) {               //Sn
      ip[0] = -1235.045;                    //1s
      ip.push_back(-262.778);               //2s
      ip.push_back(-262.778);               //2p
      ip.push_back(-83.420);                //3s
      ip.push_back(-83.420);                //3p
      ip.push_back(-46.240);                //3d
      ip.push_back(-18.081);                //4s
      ip.push_back(-18.081);                //4p
      ip.push_back(-4.300);                 //4d
      ip.push_back(0.00);                   //4f
      ip.push_back(-0.998);                 //5s
      ip.push_back(-0.998);                 //5p
    }
    else if (atom == 51) {               //Sb
      ip[0] = -1285.245;                    //1s
      ip.push_back(-274.365);               //2s
      ip.push_back(-274.365);               //2p
      ip.push_back(-87.781);                //3s
      ip.push_back(-87.781);                //3p
      ip.push_back(-49.501);                //3d
      ip.push_back(-19.743);                //4s
      ip.push_back(-19.743);                //4p
      ip.push_back(-5.129);                 //4d
      ip.push_back(0.00);                   //4f
      ip.push_back(-1.240);                 //5s
      ip.push_back(-1.240);                 //5p
    }
    else if (atom == 52) {              //Te
      ip[0] = -1336.445;                    //1s
      ip.push_back(-286.203);               //2s
      ip.push_back(-286.203);               //2p
      ip.push_back(-92.253);                //3s
      ip.push_back(-92.253);                //3p
      ip.push_back(-52.873);                //3d
      ip.push_back(-21.478);                //4s
      ip.push_back(-21.478);                //4p
      ip.push_back(-6.031);                 //4d
      ip.push_back(0.00);                   //4f
      ip.push_back(-1.509);                 //5s
      ip.push_back(-1.509);                 //5p
    }
    else if (atom == 53) {               //I
      ip[0] = -1388.645;                    //1s
      ip.push_back(-298.290);               //2s
      ip.push_back(-298.290);               //2p
      ip.push_back(-96.837);                //3s
      ip.push_back(-96.837);                //3p
      ip.push_back(-56.357);                //3d
      ip.push_back(-23.286);                //4s
      ip.push_back(-23.286);                //4p
      ip.push_back(-7.006);                 //4d
      ip.push_back(0.00);                   //4f
      ip.push_back(-1.805);                 //5s
      ip.push_back(-1.805);                 //5p
    }
    else if (atom == 54) {              //Xe
      ip[0] = -1441.845;                    //1s
      ip.push_back(-310.628);               //2s
      ip.push_back(-310.628);               //2p
      ip.push_back(-101.531);               //3s
      ip.push_back(-101.531);               //3p
      ip.push_back(-59.951);                //3d
      ip.push_back(-25.167);                //4s
      ip.push_back(-25.167);                //4p
      ip.push_back(-8.054);                 //4d
      ip.push_back(0.00);                   //4f
      ip.push_back(-2.127);                 //5s
      ip.push_back(-2.127);                 //5p
    }
    else if (atom == 55) {              //Cs
      ip[0] = -1496.045;                    //1s
      ip.push_back(-323.215);               //2s
      ip.push_back(-323.215);               //2p
      ip.push_back(-106.337);               //3s
      ip.push_back(-106.337);               //3p
      ip.push_back(-63.657);                //3d
      ip.push_back(-27.121);                //4s
      ip.push_back(-27.121);                //4p
      ip.push_back(-9.175);                 //4d
      ip.push_back(0.00);                   //4f
      ip.push_back(-2.674);                 //5s
      ip.push_back(-2.674);                 //5p
      ip.push_back(0.00);                   //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.137);                 //6s
      ip.push_back(-0.137);                 //6p
    }
    else if (atom == 56) {              //Ba
      ip[0] = -1551.245;                    //1s
      ip.push_back(-336.053);               //2s
      ip.push_back(-336.053);               //2p
      ip.push_back(-111.253);               //3s
      ip.push_back(-111.253);               //3p
      ip.push_back(-67.473);                //3d
      ip.push_back(-29.148);                //4s
      ip.push_back(-29.148);                //4p
      ip.push_back(-10.370);                //4d
      ip.push_back(0.00);                   //4f
      ip.push_back(-3.283);                 //5s
      ip.push_back(-3.283);                 //5p
      ip.push_back(0.00);                   //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.230);                 //6s
      ip.push_back(-0.230);                 //6p
    }
    else if (atom == 57) {              //La
      ip[0] = -1607.445;                    //1s
      ip.push_back(-349.140);               //2s
      ip.push_back(-349.140);               //2p
      ip.push_back(-116.281);               //3s
      ip.push_back(-116.281);               //3p
      ip.push_back(-71.401);                //3d
      ip.push_back(-31.248);                //4s
      ip.push_back(-31.248);                //4p
      ip.push_back(-11.637);                //4d
      ip.push_back(0.00);                   //4f
      ip.push_back(-3.955);                 //5s
      ip.push_back(-3.955);                 //5p
      ip.push_back(-0.281);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.255);                 //6s
      ip.push_back(-0.255);                 //6p
    }
    else if (atom == 58) {              //Ce
      ip[0] = -1664.645;                    //1s
      ip.push_back(-362.478);               //2s
      ip.push_back(-362.478);               //2p
      ip.push_back(-121.420);               //3s
      ip.push_back(-121.420);               //3p
      ip.push_back(-75.440);                //3d
      ip.push_back(-33.421);                //4s
      ip.push_back(-33.421);                //4p
      ip.push_back(-12.977);                //4d
      ip.push_back(-5.259);                 //4f
      ip.push_back(-4.061);                 //5s
      ip.push_back(-4.061);                 //5p
      ip.push_back(-0.281);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.255);                 //6s
      ip.push_back(-0.255);                 //6p
    }
    else if (atom == 59) {              //Pr
      ip[0] = -1722.845;                    //1s
      ip.push_back(-376.065);               //2s
      ip.push_back(-376.065);               //2p
      ip.push_back(-126.670);               //3s
      ip.push_back(-126.670);               //3p
      ip.push_back(-79.590);                //3d
      ip.push_back(-35.667);                //4s
      ip.push_back(-35.667);                //4p
      ip.push_back(-14.391);                //4d
      ip.push_back(-5.526);                 //4f
      ip.push_back(-3.578);                 //5s
      ip.push_back(-3.578);                 //5p
      ip.push_back(-0.173);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.230);                 //6s
      ip.push_back(-0.230);                 //6p
    }
    else if (atom == 60) {              //Nd
      ip[0] = -1782.045;                    //1s
      ip.push_back(-389.903);               //2s
      ip.push_back(-389.903);               //2p
      ip.push_back(-132.031);               //3s
      ip.push_back(-132.031);               //3p
      ip.push_back(-83.851);                //3d
      ip.push_back(-37.986);                //4s
      ip.push_back(-37.986);                //4p
      ip.push_back(-15.877);                //4d
      ip.push_back(-6.125);                 //4f
      ip.push_back(-3.679);                 //5s
      ip.push_back(-3.679);                 //5p
      ip.push_back(-0.173);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.23);                  //6s
      ip.push_back(-0.23);                  //6p
    }
    else if (atom == 61) {              //Pm
      ip[0] = -1842.245;                    //1s
      ip.push_back(-403.990);               //2s
      ip.push_back(-403.990);               //2p
      ip.push_back(-137.503);               //3s
      ip.push_back(-137.503);               //3p
      ip.push_back(-88.223);                //3d
      ip.push_back(-40.378);                //4s
      ip.push_back(-40.378);                //4p
      ip.push_back(-17.437);                //4d
      ip.push_back(-6.755);                 //4f
      ip.push_back(-3.781);                 //5s
      ip.push_back(-3.781);                 //5p
      ip.push_back(-0.173);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.230);                 //6s
      ip.push_back(-0.230);                 //6p
    }
    else if (atom == 62) {              //Sm
      ip[0] = -1903.445;                    //1s
      ip.push_back(-418.328);               //2s
      ip.push_back(-418.328);               //2p
      ip.push_back(-143.087);               //3s
      ip.push_back(-143.087);               //3p
      ip.push_back(-92.707);                //3d
      ip.push_back(-42.844);                //4s
      ip.push_back(-42.844);                //4p
      ip.push_back(-19.069);                //4d
      ip.push_back(-7.416);                 //4f
      ip.push_back(-3.885);                 //5s
      ip.push_back(-3.885);                 //5p
      ip.push_back(-0.173);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.23);                  //6s
      ip.push_back(-0.23);                  //6p
    }
    else if (atom == 63) {              //Eu
      ip[0] = -1965.645;                    //1s
      ip.push_back(-432.915);               //2s
      ip.push_back(-432.915);               //2p
      ip.push_back(-148.781);               //3s
      ip.push_back(-148.781);               //3p
      ip.push_back(-97.301);                //3d
      ip.push_back(-45.382);                //4s
      ip.push_back(-45.382);                //4p
      ip.push_back(-20.775);                //4d
      ip.push_back(-8.108);                 //4f
      ip.push_back(-3.990);                 //5s
      ip.push_back(-3.990);                 //5p
      ip.push_back(-0.173);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.23);                  //6s
      ip.push_back(-0.23);                  //6p
    }
    else if (atom == 64) {              //Gd
      ip[0] = -2028.845;                    //1s
      ip.push_back(-447.753);               //2s
      ip.push_back(-447.753);               //2p
      ip.push_back(-154.587);               //3s
      ip.push_back(-154.587);               //3p
      ip.push_back(-102.007);               //3d
      ip.push_back(-47.994);                //4s
      ip.push_back(-47.994);                //4p
      ip.push_back(-22.554);                //4d
      ip.push_back(-9.233);                 //4f
      ip.push_back(-4.728);                 //5s
      ip.push_back(-4.728);                 //5p
      ip.push_back(-0.281);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.255);                 //6s
      ip.push_back(-0.255);                 //6p
    }
    else if (atom == 65) {              //Tb
      ip[0] = -2093.045;                    //1s
      ip.push_back(-462.840);               //2s
      ip.push_back(-462.840);               //2p
      ip.push_back(-160.503);               //3s
      ip.push_back(-160.503);               //3p
      ip.push_back(-106.823);               //3d
      ip.push_back(-50.678);                //4s
      ip.push_back(-50.678);                //4p
      ip.push_back(-24.405);                //4d
      ip.push_back(-10.004);                //4f
      ip.push_back(-4.844);                 //5s
      ip.push_back(-4.844);                 //5p
      ip.push_back(-0.281);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.255);                 //6s
      ip.push_back(-0.255);                 //6p
    }
    else if (atom == 66) {              //Dy
      ip[0] = -2158.245;                    //1s
      ip.push_back(-478.178);               //2s
      ip.push_back(-478.178);               //2p
      ip.push_back(-166.531);               //3s
      ip.push_back(-166.531);               //3p
      ip.push_back(-111.751);               //3d
      ip.push_back(-53.435);                //4s
      ip.push_back(-53.435);                //4p
      ip.push_back(-26.330);                //4d
      ip.push_back(-10.370);                //4f
      ip.push_back(-4.314);                 //5s
      ip.push_back(-4.314);                 //5p
      ip.push_back(-0.173);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.230);                 //6s
      ip.push_back(-0.230);                 //6p
    }
    else if (atom == 67) {              //Ho
      ip[0] = -2224.445;                    //1s
      ip.push_back(-493.765);               //2s
      ip.push_back(-493.765);               //2p
      ip.push_back(-172.670);               //3s
      ip.push_back(-172.670);               //3p
      ip.push_back(-116.790);               //3d
      ip.push_back(-56.266);                //4s
      ip.push_back(-56.266);                //4p
      ip.push_back(-28.328);                //4d
      ip.push_back(-11.185);                //4f
      ip.push_back(-4.425);                 //5s
      ip.push_back(-4.425);                 //5p
      ip.push_back(-0.173);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.23);                  //6s
      ip.push_back(-0.23);                  //6p
    }
    else if (atom == 68) {              //Er
      ip[0] = -2291.645;                    //1s
      ip.push_back(-509.603);               //2s
      ip.push_back(-509.603);               //2p
      ip.push_back(-178.920);               //3s
      ip.push_back(-178.920);               //3p
      ip.push_back(-121.940);               //3d
      ip.push_back(-59.170);                //4s
      ip.push_back(-59.170);                //4p
      ip.push_back(-30.399);                //4d
      ip.push_back(-12.032);                //4f
      ip.push_back(-4.538);                 //5s
      ip.push_back(-4.538);                 //5p
      ip.push_back(-0.173);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.230);                 //6s
      ip.push_back(-0.230);                 //6p
    }
    else if (atom == 69) {              //Tm
      ip[0] = -2359.845;                    //1s
      ip.push_back(-525.690);               //2s
      ip.push_back(-525.690);               //2p
      ip.push_back(-185.281);               //3s
      ip.push_back(-185.281);               //3p
      ip.push_back(-127.201);               //3d
      ip.push_back(-62.146);                //4s
      ip.push_back(-62.146);                //4p
      ip.push_back(-32.543);                //4d
      ip.push_back(-12.909);                //4f
      ip.push_back(-4.651);                 //5s
      ip.push_back(-4.651);                 //5p
      ip.push_back(-0.173);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.230);                 //6s
      ip.push_back(-0.230);                 //6p
    }
    else if (atom == 70) {              //Yb
      ip[0] = -2429.045;                    //1s
      ip.push_back(-542.028);               //2s
      ip.push_back(-542.028);               //2p
      ip.push_back(-191.753);               //3s
      ip.push_back(-191.753);               //3p
      ip.push_back(-132.573);               //3d
      ip.push_back(-65.196);                //4s
      ip.push_back(-65.196);                //4p
      ip.push_back(-34.760);                //4d
      ip.push_back(-13.818);                //4f
      ip.push_back(-4.766);                 //5s
      ip.push_back(-4.766);                 //5p
      ip.push_back(-0.173);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.230);                 //6s
      ip.push_back(-0.230);                 //6p
    }
    else if (atom == 71) {              //Lu
      ip[0] = -2499.245;                    //1s
      ip.push_back(-558.615);               //2s
      ip.push_back(-558.615);               //2p
      ip.push_back(-198.337);               //3s
      ip.push_back(-198.337);               //3p
      ip.push_back(-138.057);               //3d
      ip.push_back(-68.319);                //4s
      ip.push_back(-68.319);                //4p
      ip.push_back(-37.050);                //4d
      ip.push_back(-15.274);                //4f
      ip.push_back(-5.569);                 //5s
      ip.push_back(-5.569);                 //5p
      ip.push_back(-0.281);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.255);                 //6s
      ip.push_back(-0.255);                 //6p
    }
    else if (atom == 72) {              //Hf
      ip[0] = -2570.445;                    //1s
      ip.push_back(-575.453);               //2s
      ip.push_back(-575.453);               //2p
      ip.push_back(-205.031);               //3s
      ip.push_back(-205.031);               //3p
      ip.push_back(-143.651);               //3d
      ip.push_back(-71.514);                //4s
      ip.push_back(-71.514);                //4p
      ip.push_back(-39.413);                //4d
      ip.push_back(-16.804);                //4f
      ip.push_back(-6.435);                 //5s
      ip.push_back(-6.435);                 //5p
      ip.push_back(-0.416);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.281);                 //6s
      ip.push_back(-0.281);                 //6p
    }
    else if (atom == 73) {              //Ta
      ip[0] = -2642.645;                    //1s
      ip.push_back(-592.540);               //2s
      ip.push_back(-592.540);               //2p
      ip.push_back(-211.837);               //3s
      ip.push_back(-211.837);               //3p
      ip.push_back(-149.357);               //3d
      ip.push_back(-74.783);                //4s
      ip.push_back(-74.783);                //4p
      ip.push_back(-41.849);                //4d
      ip.push_back(-18.408);                //4f
      ip.push_back(-7.363);                 //5s
      ip.push_back(-7.363);                 //5p
      ip.push_back(-0.578);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.309);                 //6s
      ip.push_back(-0.309);                 //6p
    }
    else if (atom == 74) {              //W
      ip[0] = -2715.845;                    //1s
      ip.push_back(-609.878);               //2s
      ip.push_back(-609.878);               //2p
      ip.push_back(-218.753);               //3s
      ip.push_back(-218.753);               //3p
      ip.push_back(-155.173);               //3d
      ip.push_back(-78.125);                //4s
      ip.push_back(-78.125);                //4p
      ip.push_back(-44.358);                //4d
      ip.push_back(-20.084);                //4f
      ip.push_back(-8.354);                 //5s
      ip.push_back(-8.354);                 //5p
      ip.push_back(-0.766);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.337);                 //6s
      ip.push_back(-0.337);                 //6p
    }
    else if (atom == 75) {              //Re
      ip[0] = -2790.045;                    //1s
      ip.push_back(-627.465);               //2s
      ip.push_back(-627.465);               //2p
      ip.push_back(-225.781);               //3s
      ip.push_back(-225.781);               //3p
      ip.push_back(-161.101);               //3d
      ip.push_back(-81.540);                //4s
      ip.push_back(-81.540);                //4p
      ip.push_back(-46.940);                //4d
      ip.push_back(-21.834);                //4f
      ip.push_back(-9.407);                 //5s
      ip.push_back(-9.407);                 //5p
      ip.push_back(-0.980);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.367);                 //6s
      ip.push_back(-0.367);                 //6p
    }
    else if (atom == 76) {              //Os
      ip[0] = -2865.245;                    //1s
      ip.push_back(-645.303);               //2s
      ip.push_back(-645.303);               //2p
      ip.push_back(-232.920);               //3s
      ip.push_back(-232.920);               //3p
      ip.push_back(-167.140);               //3d
      ip.push_back(-85.028);                //4s
      ip.push_back(-85.028);                //4p
      ip.push_back(-49.595);                //4d
      ip.push_back(-23.656);                //4f
      ip.push_back(-10.523);                //5s
      ip.push_back(-10.523);                //5p
      ip.push_back(-1.221);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.399);                 //6s
      ip.push_back(-0.399);                 //6p
    }
    else if (atom == 77) {              //Ir
      ip[0] = -2491.445;                    //1s
      ip.push_back(-663.390);               //2s
      ip.push_back(-663.390);               //2p
      ip.push_back(-240.170);               //3s
      ip.push_back(-240.170);               //3p
      ip.push_back(-173.290);               //3d
      ip.push_back(-88.589);                //4s
      ip.push_back(-88.589);                //4p
      ip.push_back(-52.324);                //4d
      ip.push_back(-25.552);                //4f
      ip.push_back(-11.701);                //5s
      ip.push_back(-11.701);                //5p
      ip.push_back(-1.488);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.431);                 //6s
      ip.push_back(-0.431);                 //6p
    }
    else if (atom == 78) {              //Pt
      ip[0] = -3018.645;                    //1s
      ip.push_back(-681.728);               //2s
      ip.push_back(-681.728);               //2p
      ip.push_back(-247.531);               //3s
      ip.push_back(-247.531);               //3p
      ip.push_back(-179.551);               //3d
      ip.push_back(-92.223);                //4s
      ip.push_back(-92.223);                //4p
      ip.push_back(-55.125);                //4d
      ip.push_back(-27.520);                //4f
      ip.push_back(-12.941);                //5s
      ip.push_back(-12.941);                //5p
      ip.push_back(-1.781);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.465);                 //6s
      ip.push_back(-0.465);                 //6p
    }
    else if (atom == 79) {              //Au
      ip[0] = -3096.845;                    //1s
      ip.push_back(-700.315);               //2s
      ip.push_back(-700.315);               //2p
      ip.push_back(-255.003);               //3s
      ip.push_back(-255.003);               //3p
      ip.push_back(-185.923);               //3d
      ip.push_back(-95.930);                //4s
      ip.push_back(-95.930);                //4p
      ip.push_back(-57.999);                //4d
      ip.push_back(-29.562);                //4f
      ip.push_back(-14.244);                //5s
      ip.push_back(-14.244);                //5p
      ip.push_back(-2.101);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.500);                 //6s
      ip.push_back(-0.500);                 //6p
    }
    else if (atom == 80) {              //Hg
      ip[0] = -3176.045;                    //1s
      ip.push_back(-719.153);               //2s
      ip.push_back(-719.153);               //2p
      ip.push_back(-262.587);               //3s
      ip.push_back(-262.587);               //3p
      ip.push_back(-192.407);               //3d
      ip.push_back(-99.710);                //4s
      ip.push_back(-99.710);                //4p
      ip.push_back(-60.947);                //4d
      ip.push_back(-31.676);                //4f
      ip.push_back(-15.610);                //5s
      ip.push_back(-15.610);                //5p
      ip.push_back(-2.448);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.536);                 //6s
      ip.push_back(-0.536);                 //6p
    }
    else if (atom == 81) {              //Tl
      ip[0] = -3256.245;                    //1s
      ip.push_back(-738.240);               //2s
      ip.push_back(-738.240);               //2p
      ip.push_back(-270.281);               //3s
      ip.push_back(-270.281);               //3p
      ip.push_back(-199.001);               //3d
      ip.push_back(-103.563);               //4s
      ip.push_back(-103.563);               //4p
      ip.push_back(-63.967);                //4d
      ip.push_back(-33.864);                //4f
      ip.push_back(-17.038);                //5s
      ip.push_back(-17.038);                //5p
      ip.push_back(-3.032);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.709);                 //6s
      ip.push_back(-0.709);                 //6p
    }
    else if (atom == 82) {              //Pb
      ip[0] = -3337.445;                    //1s
      ip.push_back(-757.578);               //2s
      ip.push_back(-757.578);               //2p
      ip.push_back(-278.087);               //3s
      ip.push_back(-278.087);               //3p
      ip.push_back(-205.707);               //3d
      ip.push_back(-107.489);               //4s
      ip.push_back(-107.489);               //4p
      ip.push_back(-67.061);                //4d
      ip.push_back(-36.125);                //4f
      ip.push_back(-18.529);                //5s
      ip.push_back(-18.529);                //5p
      ip.push_back(-3.679);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-0.905);                 //6s
      ip.push_back(-0.905);                 //6p
    }
    else if (atom == 83) {              //Bi
      ip[0] = -3419.645;                    //1s
      ip.push_back(-777.165);               //2s
      ip.push_back(-777.165);               //2p
      ip.push_back(-286.003);               //3s
      ip.push_back(-286.003);               //3p
      ip.push_back(-212.523);               //3d
      ip.push_back(-111.489);               //4s
      ip.push_back(-111.489);               //4p
      ip.push_back(-70.227);                //4d
      ip.push_back(-38.459);                //4f
      ip.push_back(-20.082);                //5s
      ip.push_back(-20.082);                //5p
      ip.push_back(-4.388);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-1.125);                 //6s
      ip.push_back(-1.125);                 //6p
    }
    else if (atom == 84) {              //Po
      ip[0] = -3502.845;                    //1s
      ip.push_back(-797.003);               //2s
      ip.push_back(-797.003);               //2p
      ip.push_back(-294.031);               //3s
      ip.push_back(-294.031);               //3p
      ip.push_back(-219.451);               //3d
      ip.push_back(-115.561);               //4s
      ip.push_back(-115.561);               //4p
      ip.push_back(-73.467);                //4d
      ip.push_back(-40.866);                //4f
      ip.push_back(-21.698);                //5s
      ip.push_back(-21.698);                //5p
      ip.push_back(-5.160);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-1.369);                 //6s
      ip.push_back(-1.369);                 //6p
    }
    else if (atom == 85) {              //At
      ip[0] = -3587.045;                    //1s
      ip.push_back(-817.090);               //2s
      ip.push_back(-817.090);               //2p
      ip.push_back(-302.170);               //3s
      ip.push_back(-302.170);               //3p
      ip.push_back(-226.490);               //3d
      ip.push_back(-119.706);               //4s
      ip.push_back(-119.706);               //4p
      ip.push_back(-76.779);                //4d
      ip.push_back(-43.346);                //4f
      ip.push_back(-23.376);                //5s
      ip.push_back(-23.376);                //5p
      ip.push_back(-5.994);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-1.637);                 //6s
      ip.push_back(-1.637);                 //6p
    }
    else if (atom == 86) {              //Rn
      ip[0] = -3672.245;                    //1s
      ip.push_back(-837.428);               //2s
      ip.push_back(-837.428);               //2p
      ip.push_back(-310.420);               //3s
      ip.push_back(-310.420);               //3p
      ip.push_back(-233.640);               //3d
      ip.push_back(-123.925);               //4s
      ip.push_back(-123.925);               //4p
      ip.push_back(-80.165);                //4d
      ip.push_back(-45.899);                //4f
      ip.push_back(-25.116);                //5s
      ip.push_back(-25.116);                //5p
      ip.push_back(-6.891);                 //5d
      ip.push_back(0.00);                   //5f
      ip.push_back(-1.929);                 //6s
      ip.push_back(-1.929);                 //6p
    }
    else {                              //Not implemented
      ip[0] = 0.00;
    }
  }
  if (IPtype == 1) {                        //Values collected from literature for extended Hueckel
    //paper reference A: J. Phys. Chem. vol. 92(20), 5691, 1988 (Extended Hueckel Parameters from DFT)
    //paper reference B: Russ. J. In. Chem. vol. 19(11), 1591, 1974 (The Orbital Ionisation Potentials Of The Atoms And Ions Of Heavy Elements)
    if (atom == 1) {ip[0] = -0.50;}  //H       [A]
    else if (atom == 3) {            //Li      [A]
      ip[0] = -0.198;                //s
      ip.push_back(-0.129);          //p
    }
    else if (atom == 4) {            //Be      [A]
      ip[0] = -0.367;                //s
      ip.push_back(-0.220);          //p
    }
    else if (atom == 5) {            //B       [A]
      ip[0] = -0.559;                //s
      ip.push_back(-0.312);          //p
    }
    else if (atom == 6) {            //C       [A]
      ip[0] = -0.786;                //s
      ip.push_back(-0.419);          //p
    }
    else if (atom == 7) {            //N       [A]
      ip[0] = -0.955;                //s
      ip.push_back(-0.492);          //p
    }
    else if (atom == 8) {            //O       [A]
      ip[0] = -1.187;                //s
      ip.push_back(-0.544);          //p
    }
    else if (atom == 9) {            //F       [A]
      ip[0] = -1.470;                //s
      ip.push_back(-0.665);          //p
    }
    else if (atom == 11) {           //Na      [A]
      ip[0] = -0.187;                //s
      ip.push_back(-0.110);          //p
    }
    else if (atom == 13) {           //Al      [A]
      ip[0] = -0.452;                //s
      ip.push_back(-0.239);          //p
    }
    else if (atom == 14) {           //Si      [A]
      ip[0] = -0.636;                //s
      ip.push_back(-0.338);          //p
    }
    else if (atom == 15) {           //P       [A]
      ip[0] = -0.684;                //s
      ip.push_back(-0.514);          //p
    }
    else if (atom == 16) {           //S       [A]
      ip[0] = -0.735;                //s
      ip.push_back(-0.404);          //p
    }
    else if (atom == 17) {           //Cl      [A]
      ip[0] = -0.967;                //s
      ip.push_back(-0.522);          //p
    }
    else if (atom == 19) {           //K       [A]
      ip[0] = -0.159;                //s
      ip.push_back(-0.100);          //p
    }
    else if (atom == 20) {           //Ca      [A]
      ip[0] = -0.257;                //s
      ip.push_back(-0.147);            //p
    }
    else if (atom == 22) {           //Ti      [A]
      ip[0] = -0.217;                //d
      ip.push_back(-0.232);          //s
      ip.push_back(-0.118);          //p
      if (hardness == "hard") {
        ip[0] = -0.397;
        ip[1] = -0.330;
        ip[2] = -0.200;
      }
    }
    else if (atom == 23) {           //V       [A]
      ip[0] = -0.246;                //d
      ip.push_back(-0.246);          //s
      ip.push_back(-0.125);          //p
      if (hardness == "hard") {
        ip[0] = -0.404;
        ip[1] = -0.324;
        ip[2] = -0.203;
      }
    }
    else if (atom == 24) {           //Cr      [A]
      ip[0] = -0.290;                //d
      ip.push_back(-0.268);          //s
      ip.push_back(-0.132);          //p
      if (hardness == "hard") {
        ip[0] = -0.411;
        ip[1] = -0.318;
        ip[2] = -0.193;
      }
    }
    else if (atom == 25) {           //Mn      [A]
      ip[0] = -0.320;                //d
      ip.push_back(-0.276);          //s
      ip.push_back(-0.140);          //p
      if (hardness == "hard") {
        ip[0] = -0.429;
        ip[1] = -0.358;
        ip[2] = -0.216;
      }
    }
    else if (atom == 26) {           //Fe      [A]
      ip[0] = -0.338;                //d
      ip.push_back(-0.279);          //s
      ip.push_back(-0.140);          //p
      if (hardness == "hard") {
        ip[0] = -0.463;
        ip[1] = -0.334;
        ip[2] = -0.196;
      }
    }
    else if (atom == 27) {           //Co      [A]
      ip[0] = -0.356;                //d
      ip.push_back(-0.287);          //s
      ip.push_back(-0.140);          //p
      if (hardness == "hard") {
        ip[0] = -0.484;
        ip[1] = -0.338;
        ip[2] = -0.194;
      }
    }
    else if (atom == 28) {           //Ni      [A]
      ip[0] = -0.364;                //d
      ip.push_back(-0.287);          //s
      ip.push_back(-0.136);          //p
      if (hardness == "hard") {
        ip[0] = -0.496;
        ip[1] = -0.337;
        ip[2] = -0.189;
      }
    }
    else if (atom == 29) {           //Cu      [A]
      ip[0] = -0.514;                //d
      ip.push_back(-0.419);          //s
      ip.push_back(-0.223);          //p
    }
    else if (atom == 30) {           //Zn      [A]
      ip.push_back(-0.456);          //s
      ip.push_back(-0.240);          //p
    }
    else if (atom == 31) {           //Ga      [A]
      ip[0] = -0.536;                //s
      ip.push_back(-0.248);          //p
    }
    else if (atom == 32) {           //Ge      [A]
      ip[0] = -0.588;                //s
      ip.push_back(-0.331);          //p
    }
    else if (atom == 34) {           //Se      [A]
      ip[0] = -0.753;                //s
      ip.push_back(-0.529);          //p
    }
    else if (atom == 35) {           //Br      [A]
      ip[0] = -0.811;                //s
      ip.push_back(-0.481);          //p
    }
    else if (atom == 37) {           //Rb      [A]
      ip[0] = -0.154;                 //s
      ip.push_back(-0.096);           //p
    }
    else if (atom == 38) {           //Sr      [A]
      ip[0] = -0.243;                 //s
      ip.push_back(-0.144);           //p
    }
    else if (atom == 41) {           //Nb      [A]
      ip[0] = -0.445;                //d
      ip.push_back(-0.371);          //s
      ip.push_back(-0.252);          //p
    }
    else if (atom == 42) {           //Mo      [A]
      ip[0] = -0.386;                //d
      ip.push_back(-0.306);          //s
      ip.push_back(-0.193);          //p
    }
    else if (atom == 44) {           //Ru      [A]
      ip[0] = -0.413;                //d
      ip.push_back(-0.284);          //s
      ip.push_back(-0.163);          //p
      if (hardness == "hard") {
        ip[0] = -0.548;
        ip[1] = -0.382;
        ip[2] = -0.252;
      }
    }
    else if (atom == 45) {           //Rh      [A]
      ip[0] = -0.459;                //d
      ip.push_back(-0.297);          //s
      ip.push_back(-0.168);          //p
      if (hardness == "hard") {
        ip[0] = -0.468;
        ip[1] = -0.331;
        ip[2] = -0.166;
      }
    }
    else if (atom == 46) {           //Pd      [A]
      ip[0] = -0.442;                //d
      ip.push_back(-0.269);          //s
      ip.push_back(-0.138);          //p
    }
    else if (atom == 49) {           //In      [A]
      ip[0] = -0.463;                //s
      ip.push_back(-0.227);          //p
    }
    else if (atom == 50) {           //Sn      [A]
      ip[0] = -0.594;                //s
      ip.push_back(-0.306);          //p
    }
    else if (atom == 51) {           //Sb      [A]
      ip[0] = -0.691;                //s
      ip.push_back(-0.430);          //p
    }
    else if (atom == 52) {           //Te      [A]
      ip[0] = -0.764;                //s
      ip.push_back(-0.544);          //p
    }
    else if (atom == 53) {           //I       [A]
      ip[0] = -0.661;                //s
      ip.push_back(-0.467);          //p
    }
    else if (atom == 55) {           //Cs      [A]
      ip[0] = -0.143;                //s
      ip.push_back(-0.092);          //p
    }
    else if (atom == 56) {           //Ba      [B]
      ip[0] = -0.175;                //s
      ip.push_back(-0.097);          //p
    }
    else if (atom == 73) {           //Ta      [A]
      ip[0] = -0.445;                //d
      ip.push_back(-0.371);          //s
      ip.push_back(-0.252);          //p
    }
    else if (atom == 74) {           //W       [A]
      ip[0] = -0.381;                //d
      ip.push_back(-0.304);          //s
      ip.push_back(-0.190);          //p
    }
    else if (atom == 75) {           //Re      [A]
      ip[0] = -0.465;                //d
      ip.push_back(-0.344);          //s
      ip.push_back(-0.219);          //p
    }
    else if (atom == 77) {           //Ir      [A]
      ip[0] = -0.447;                //d
      ip.push_back(-0.417);          //s
      ip.push_back(-0.165);          //p
    }
    else if (atom == 78) {           //Pt      [A]
      ip[0] = -0.463;                //d
      ip.push_back(-0.334);          //s
      ip.push_back(-0.201);          //p
    }
    else if (atom == 79) {           //Au      [A]
      ip[0] = -0.554;                //d
      ip.push_back(-0.401);          //s
      ip.push_back(-0.204);          //p
    }
    else if (atom == 82) {           //Pb      [A]
      ip[0] = -0.577;                //s
      ip.push_back(-0.294);          //p
    }
    else if (atom == 83) {           //Bi      [A]
      ip[0] = -0.558;                //s
      ip.push_back(-0.286);          //p
    }
    else {                              //Not implemented
      ip[0] = 0.00;
    }
  }
  else if (IPtype == 2) {                   //Skinner and Pritchard values
    if (atom == 1) {ip[0] = -13.595;}//H
    else if (atom == 3) {            //Li
      ip[0] = -5.39;                 //s
      ip.push_back(-3.5);            //p  taken from Hoffmann
    }
    else if (atom == 4) {            //Be
      ip[0] = -9.92;                 //s
      ip.push_back(-5.96);           //p
    }
    else if (atom == 5) {            //B
      ip[0] = -15.17;                //s
      ip.push_back(-8.53);           //p
    }
    else if (atom == 6) {            //C
      ip[0] = -21.34;                //s
      ip.push_back(-11.54);          //p
    }
    else if (atom == 7) {            //N
      ip[0] = -27.5;                 //s
      ip.push_back(-14.32);          //p
    }
    else if (atom == 8) {            //O
      ip[0] = -35.3;                 //s
      ip.push_back(-17.91);          //p
    }
    else if (atom == 9) {            //F
      ip[0] = -40.0;                 //s  taken from Hoffmann
      ip.push_back(-20.88);          //p
    }
    else if (atom == 11) {           //Na
      ip[0] = -5.14;                 //s
      ip.push_back(-3.0);            //p  taken from Hoffmann
    }
    else if (atom == 12) {           //Mg
      ip[0] = -8.95;                 //s
      ip.push_back(-4.52);           //p
    }
    else if (atom == 13) {           //Al
      ip[0] = -12.48;                //s
      ip.push_back(-6.41);           //p
    }
    else if (atom == 14) {           //Si
      ip[0] = -17.93;                //s
      ip.push_back(-8.99);           //p
      ip.push_back(-6.0);            //d  taken from Hoffmann
    }
    else if (atom == 15) {           //P
      ip[0] = -18.6;                 //s  taken from Hoffmann
      ip.push_back(-12.78);          //p
      ip.push_back(-7.0);            //d  taken from Hoffmann
    }
    else if (atom == 16) {           //S
      ip[0] = -20.0;                 //s  taken from Hoffmann
      ip.push_back(-13.31);          //p
      ip.push_back(-8.0);            //d  taken from Hoffmann
    }
    else if (atom == 17) {           //Cl
      ip[0] = -30.0;                 //s  taken from Hoffmann
      ip.push_back(-15.08);          //p
      ip.push_back(-9.0);            //d  taken from Hoffmann
    }
    else {                              //Not implemented
      ip[0] = 0.00;
    }
  }
  else if (IPtype == 3) {                   //Hinze and Jaffé
    if (atom == 1) {ip[0] = -13.595;}//H
    else if (atom == 5) {            //B
      ip[0] = -14.912;               //s
      ip.push_back(-8.421);          //p
    }
    else if (atom == 6) {            //C
      ip[0] = -21.012;               //s
      ip.push_back(-11.273);         //p
    }
    else if (atom == 7) {            //N
      ip[0] = -25.588;               //s
      ip.push_back(-13.946);         //p
    }
    else if (atom == 8) {            //O
      ip[0] = -32.297;               //s
      ip.push_back(-17.274);         //p
    }
    else if (atom == 9) {            //F
      ip[0] = -39.391;               //s
      ip.push_back(-20.862);         //p
    }
    else if (atom == 14) {           //Si
      ip[0] = -17.307;               //s
      ip.push_back(-9.190);          //p
    }
    else if (atom == 15) {           //P
      ip[0] = -18.612;               //s
      ip.push_back(-10.733);         //p
    }
    else if (atom == 16) {           //S
      ip[0] = -21.135;               //s
      ip.push_back(-12.396);         //p
    }
    else if (atom == 17) {           //Cl
      ip[0] = -25.227;               //s
      ip.push_back(-15.037);         //p
    }
    else if (atom == 35) {           //Br
      ip[0] = -23.735;               //s
      ip.push_back(-13.101);         //p
    }
    else if (atom == 53) {           //I
      ip[0] = -20.833;               //s
      ip.push_back(-12.670);         //p
    }
    else {                              //Not implemented
      ip[0] = 0.00;
    }
  }
  else if (IPtype == 4) {                   //Boehm and Gleiter
    if (atom == 1) {ip[0] = -13.06;} //H
    else if (atom == 3) {            //Li
      ip[0] = -5.39;                 //s
      ip.push_back(-3.54);           //p
    }
    else if (atom == 4) {            //Be
      ip[0] = -9.32;                 //s
      ip.push_back(-5.96);           //p
    }
    else if (atom == 5) {            //B
      ip[0] = -14.05;                //s
      ip.push_back(-8.30);           //p
    }
    else if (atom == 6) {            //C
      ip[0] = -19.44;                //s
      ip.push_back(-10.67);          //p
    }
    else if (atom == 7) {            //N
      ip[0] = -25.58;                //s
      ip.push_back(-13.19);          //p
    }
    else if (atom == 8) {            //O
      ip[0] = -32.38;                //s
      ip.push_back(-15.85);          //p
    }
    else if (atom == 9) {            //F
      ip[0] = -40.42;                //s
      ip.push_back(-18.66);          //p
    }
    else if (atom == 11) {           //Na
      ip[0] = -5.14;                 //s
      ip.push_back(-3.04);           //p
      ip.push_back(-1.52);           //d
    }
    else if (atom == 12) {           //Mg
      ip[0] = -7.64;                 //s
      ip.push_back(-4.52);           //p
      ip.push_back(-1.74);           //d
    }
    else if (atom == 13) {           //Al
      ip[0] = -11.32;                //s
      ip.push_back(-5.97);           //p
      ip.push_back(-1.96);           //d
    }
    else if (atom == 14) {           //Si
      ip[0] = -15.15;                //s
      ip.push_back(-7.62);           //p
      ip.push_back(-2.04);           //d
    }
    else if (atom == 15) {           //P
      ip[0] = -19.37;                //s
      ip.push_back(-10.84);          //p
      ip.push_back(-2.86);           //d
    }
    else if (atom == 16) {           //S
      ip[0] = -20.52;                //s
      ip.push_back(-10.78);          //p
      ip.push_back(-3.06);           //d
    }
    else if (atom == 17) {           //Cl
      ip[0] = -25.29;                //s
      ip.push_back(-13.99);          //p
      ip.push_back(-3.59);           //d
    }
    else if (atom == 19) {           //K
      ip[0] = -4.44;                 //s
      ip.push_back(-2.83);           //p
      ip.push_back(-1.77);           //d
    }
    else if (atom == 20) {           //Ca
      ip[0] = -6.11;                 //s
      ip.push_back(-3.96);           //p
      ip.push_back(-3.54);           //d
    }
    else if (atom == 21) {           //Sc
      ip[0] = -5.70;                 //s
      ip.push_back(-3.22);           //p
      ip.push_back(-4.71);           //d
    }
    else if (atom == 22) {           //Ti
      ip[0] = -6.07;                 //s
      ip.push_back(-3.35);           //p
      ip.push_back(-5.58);           //d
    }
    else if (atom == 23) {           //V
      ip[0] = -6.32;                 //s
      ip.push_back(-3.47);           //p
      ip.push_back(-6.32);           //d
    }
    else if (atom == 24) {           //Cr
      ip[0] = -6.57;                 //s
      ip.push_back(-3.47);           //p
      ip.push_back(-7.19);           //d
    }
    else if (atom == 25) {           //Mn
      ip[0] = -6.82;                 //s
      ip.push_back(-3.59);           //p
      ip.push_back(-7.93);           //d
    }
    else if (atom == 26) {           //Fe
      ip[0] = -7.06;                 //s
      ip.push_back(-3.72);           //p
      ip.push_back(-8.67);           //d
    }
    else if (atom == 27) {           //Co
      ip[0] = -7.31;                 //s
      ip.push_back(-3.84);           //p
      ip.push_back(-9.42);           //d
    }
    else if (atom == 28) {           //Ni
      ip[0] = -7.56;                 //s
      ip.push_back(-3.84);           //p
      ip.push_back(-10.04);          //d
    }
    else if (atom == 29) {           //Cu
      ip[0] = -7.68;                 //s
      ip.push_back(-3.97);           //p
      ip.push_back(-10.66);          //d
    }
    else if (atom == 30) {           //Zn
      ip[0] = -7.81;                 //s
      ip.push_back(-4.09);           //p
      ip.push_back(-11.17);          //d
    }
    else if (atom == 31) {           //Ga
      ip[0] = -12.64;                //s
      ip.push_back(-5.95);           //p
    }
    else if (atom == 32) {           //Ge
      ip[0] = -15.61;                //s
      ip.push_back(-7.56);           //p
    }
    else if (atom == 33) {           //As
      ip[0] = -17.60;                //s
      ip.push_back(-9.05);           //p
    }
    else if (atom == 34) {           //Se
      ip[0] = -20.82;                //s
      ip.push_back(-10.78);          //p
    }
    else if (atom == 35) {           //Br
      ip[0] = -24.04;                //s
      ip.push_back(-12.52);          //p
    }
    else {                              //Not implemented
      ip[0] = 0.00;
    }
  }
  if (IPtype == 5) {
    if (atom == 1) {ip[0] = -0.8047249934404999;}   //H
    else if (atom == 2) {            //He
      ip[0] = -5.488660372807999;
      ip.push_back(-6.6873108898444995);
    }
    else if (atom == 3) {            //Li
      ip[0] = -0.4812207218405001;
      ip.push_back(-4.3614186840124995);
    }
    else if (atom == 4) {            //Be
      ip[0] = -0.7351254132605001;
      ip.push_back(-0.8147095305844999);
    }
    else if (atom == 5) {            //B
      ip[0] = -1.335262331138;
      ip.push_back(-1.0940089240125002);
    }
    else if (atom == 6) {            //C
      ip[0] = -2.096246881682;
      ip.push_back(-1.4498337356405);
    }
    else if (atom == 7) {            //N
      ip[0] = -2.8331663624179995;
      ip.push_back(-1.9984922842580002);
    }
    else if (atom == 8) {            //O
      ip[0] = -14.697691953000502;
      ip.push_back(-2.5786296608000003);
    }
    else if (atom == 9) {            //F
      ip[0] = -18.2640553674005;
      ip.push_back(-4.224516392641999);
    }
    else if (atom == 10) {           //Ne
      ip[0] = -18.000888010952;
      ip.push_back(-7.351802491392001);
    }
    else if (atom == 11) {           //Na
      ip[0] = -0.23552237546450003;
      ip.push_back(-0.45131460231200005);
    }
    else if (atom == 12) {           //Mg
      ip[0] = -0.85913764445;
      ip.push_back(-0.9645174383045001);
    }
    else if (atom == 13) {           //Al
      ip[0] = -2.794872130848;
      ip.push_back(-1.529678903202);
      ip.push_back(-0.8056678697280001);
    }
    else if (atom == 14) {           //Si
      ip[0] = -1.5360505065405001;
      ip.push_back(-0.7180968592845);
      ip.push_back(-2.2654540798245);
    }
    else if (atom == 15) {           //P
      ip[0] = -2.3285532145445003;
      ip.push_back(-1.6296316738244998);
      ip.push_back(-0.756890404082);
    }
    else if (atom == 16) {           //S
      ip[0] = -2.404282404168;
      ip.push_back(-1.694784101042);
      ip.push_back(-4.8341872894005);
    }
    else if (atom == 17) {           //Cl
      ip[0] = -3.4770163512499996;
      ip.push_back(-2.2432712386579996);
      ip.push_back(-0.8765316925445);
    }
    else if (atom == 18) {           //Ar
      ip[0] = -18.001632036992;
      ip.push_back(-17.696311844449998);
    }
    else if (atom == 19) {           //K
      ip[0] = -18.002868114242002;
      ip.push_back(-0.6356315075044999);
    }
    else if (atom == 20) {           //Ca
      ip[0] = -1.167786257282;
      ip.push_back(-2.1219936444179996);
    }
    else if (atom == 21) {           //Sc
      ip[0] = -1.7279628080720002;
      ip.push_back(-0.9834596479805);
      ip.push_back(-0.904776139208);
    }
    else if (atom == 22) {           //Ti
      ip[0] = -1.0057590792;
      ip.push_back(-14.1766250498645);
      ip.push_back(-0.6775271543120002);
    }
    else if (atom == 23) {           //V
      ip[0] = -0.972741888818;
      ip.push_back(-1.9489894744499998);
      ip.push_back(-0.5650971836180001);
    }
    else if (atom == 24) {           //Cr
      ip[0] = -1.3172576440805;
      ip.push_back(-5.390554785799999);
      ip.push_back(-0.5298260036179999);
    }
    else if (atom == 25) {           //Mn
      ip[0] = -3.4003104200000003;
      ip.push_back(-2.2720298111999995);
      ip.push_back(-1.1641548871999998);
    }
    else if (atom == 26) {           //Fe
      ip[0] = -0.5840070390044999;
      ip.push_back(-1.0939423612499999);
      ip.push_back(-18.013478522258);
    }
    else if (atom == 27) {           //Co
      ip[0] = -1.7302055037619999;
      ip.push_back(-0.6804929458844999);
      ip.push_back(-4.5);
    }
    else if (atom == 28) {           //Ni
      ip[0] = -3.1620114435605;
      ip.push_back(-1.266958190792);
      ip.push_back(-2.6559109290605);
    }
    else if (atom == 29) {           //Cu
      ip[0] = -3.7400851500499996;
      ip.push_back(-1.3929407286079998);
      ip.push_back(-4.5);
    }
    else if (atom == 30) {           //Zn
      ip[0] = -1.1443953828125;
      ip.push_back(-1.601122914162);
    }
    else if (atom == 31) {           //Ga
      ip[0] = -2.7356172152445;
      ip.push_back(-1.495744243232);
    }
    else if (atom == 32) {           //Ge
      ip[0] = -3.2412438606644995;
      ip.push_back(-1.46056267845);
    }
    else if (atom == 33) {           //As
      ip[0] = -4.2812383606205;
      ip.push_back(-1.5579496332404998);
      ip.push_back(-0.969029674082);
    }
    else if (atom == 34) {           //Se
      ip[0] = -3.1559914589780003;
      ip.push_back(-2.015180697888);
    }
    else if (atom == 35) {           //Br
      ip[0] = -10.907644513927998;
      ip.push_back(-2.0718866059380003);
      ip.push_back(-1.1567676514805);
    }
    else if (atom == 36) {           //Kr
      ip[0] = -0.8609974067520001;
      ip.push_back(-10.0862067298205);
    }
    else if (atom == 37) {           //Rb
      ip[0] = -15.180848960512497;
      ip.push_back(-0.89133946445);
    }
    else if (atom == 38) {           //Sr
      ip[0] = -2.4140702369044997;
      ip.push_back(-1.4966870193845);
    }
    else if (atom == 39) {           //Y
      ip[0] = -1.3627238487245;
      ip.push_back(-0.176042791712);
      ip.push_back(-1.110678869042);
    }
    else if (atom == 40) {           //Zr
      ip[0] = -1.2283588408319999;
      ip.push_back(-1.4324304540500001);
      ip.push_back(-1.4363701235280002);
    }
    else if (atom == 41) {           //Nb
      ip[0] = -1.9549051004880003;
      ip.push_back(-2.7743361679219998);
      ip.push_back(-0.9617555133244998);
    }
    else if (atom == 42) {           //Mo
      ip[0] = -1.6692422155520001;
      ip.push_back(-0.5622548320205001);
      ip.push_back(-0.9118062848719999);
    }
    else if (atom == 43) {           //Tc
      ip[0] = -1.5617806848;
      ip.push_back(-1.9134472500125);
      ip.push_back(-18.0378138387005);
    }
    else if (atom == 44) {           //Ru
      ip[0] = -2.1906677654479996;
      ip.push_back(-1.0646250240125);
      ip.push_back(-15.330297457200498);
    }
    else if (atom == 45) {           //Rh
      ip[0] = -4.209078388418001;
      ip.push_back(-0.8777051782805);
      ip.push_back(-9.271295972160498);
    }
    else if (atom == 46) {           //Pd
      ip[0] = -2.4638914296605;
      ip.push_back(-1.3753161005045);
      ip.push_back(-0.6689982657619999);
    }
    else if (atom == 47) {           //Ag
      ip[0] = -18.043994849792;
      ip.push_back(-1.9880259760080001);
      ip.push_back(-0.2324372107445);
    }
    else if (atom == 48) {           //Cd
      ip[0] = -0.9578774778319998;
      ip.push_back(-1.9157328262845001);
    }
    else if (atom == 49) {           //In
      ip[0] = -2.0464405047844996;
      ip.push_back(-2.2189196989620004);
    }
    else if (atom == 50) {           //Sn
      ip[0] = -2.8415873457405003;
      ip.push_back(-2.117492668232);
    }
    else if (atom == 51) {           //Sb
      ip[0] = -2.858866113842;
      ip.push_back(-1.571775138018);
      ip.push_back(-3.0395670240500006);
    }
    else if (atom == 52) {           //Te
      ip[0] = -3.8360677495219995;
      ip.push_back(-1.4987327398805002);
    }
    else if (atom == 53) {           //I
      ip[0] = -10.1189394072045;
      ip.push_back(-1.8375825265920003);
      ip.push_back(-1.7581406403125002);
    }
    else if (atom == 54) {           //Xe
      ip[0] = -3.8082121426845004;
      ip.push_back(-1.9551463414580001);
    }
    else if (atom == 55) {           //Cs
      ip[0] = -17.737015648032;
      ip.push_back(-1.3113658326125002);
    }
    else if (atom == 56) {           //Ba
      ip[0] = -0.9735412768204998;
      ip.push_back(-1.0226487796605002);
    }
    else if (atom == 57) {           //La
      ip[0] = -1.425620813922;
      ip.push_back(-3.5745497441999996);
      ip.push_back(-0.7789916344319999);
    }
    else if (atom == 71) {           //Lu
      ip[0] = -2.477297597832;
      ip.push_back(-14.969974785540499);
      ip.push_back(-1.4659787958080002);
    }
    else if (atom == 72) {           //Hf
      ip[0] = -1.6943459528;
      ip.push_back(-4.7596737991680005);
      ip.push_back(-1.2416027603805002);
    }
    else if (atom == 73) {           //Ta
      ip[0] = -1.6895796930005;
      ip.push_back(-10.4794402897845);
      ip.push_back(-11.718821733767998);
    }
    else if (atom == 74) {           //W
      ip[0] = -1.60993568;
      ip.push_back(-3.5499399967999996);
      ip.push_back(-1.31870424005);
    }
    else if (atom == 75) {           //Re
      ip[0] = -3.1821741453779997;
      ip.push_back(-2.9084836809605004);
      ip.push_back(-1.6477496266004998);
    }
    else if (atom == 76) {           //Os
      ip[0] = -1.5763244124500002;
      ip.push_back(-4.5934805);
      ip.push_back(-1.2703542408000001);
    }
    else if (atom == 77) {           //Ir
      ip[0] = -3.5806137731045;
      ip.push_back(-1.1263609113245);
      ip.push_back(-8.431149607564498);
    }
    else if (atom == 78) {           //Pt
      ip[0] = -5.020811498951999;
      ip.push_back(-2.6479079988480003);
      ip.push_back(-1.381793529608);
    }
    else if (atom == 79) {           //Au
      ip[0] = -12.767248364944502;
      ip.push_back(-1.6456045802805);
      ip.push_back(-1.3100252418245);
    }
    else if (atom == 80) {           //Hg
      ip[0] = -2.2152935854080003;
      ip.push_back(-1.1495722309244998);
    }
    else if (atom == 81) {           //Tl
      ip[0] = -5.5640576948445;
      ip.push_back(-1.5596270159405);
    }
    else if (atom == 82) {           //Pb
      ip[0] = -2.8058459739005004;
      ip.push_back(-1.420027040258);
    }
    else if (atom == 83) {           //Bi
      ip[0] = -6.853797725064499;
      ip.push_back(-1.7528041974645001);
    }
  }
  return ip;
}

#endif //_Ionization_Potential_Package_
