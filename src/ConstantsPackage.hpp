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

#ifndef _Constants_Package_
#define _Constants_Package_

//description:
//definition of some constants

//CCCBDB
//double h    = 6.626076e-34;                                             //Planck constant in J.s
//double KB   = 1.38066e-23;                                              //Boltzmann constant in J/K
//double NA   = 6.022137e23;                                              //Avogadro number in 1/mol

//xTB
//double KB = 3.166808578545117e-06;                                          //Boltzmann constant in Eh/K

//whenever possible from IUPAC's gold book
//https://goldbook.iupac.org/indexes/constants
double h    = 6.6260689633e-34;                                             //Planck constant in J.s
double hbar = 1.05457162853e-34;                                            //Planck constant divided by 2{\pi}
double KB   = 1.380650424e-23;                                              //Boltzmann constant in J/K
double R    = 8.31447215;                                                   //ideal gas constant in J/(K.mol)
double NA   = 6.02214076e23;                                                //Avogadro number in 1/mol
double c    = 299792458.0;                                                  //speed of light in m/s
double Mel  = 9.109389754e-31;                                              //mass of electron in kg

double pi   = 3.1415926535897932384626433832795029;
double sqrtpi = sqrt(pi);
double osqrtpi = 1.0/sqrtpi;
double emgamma = 0.57721566490153286060651209008240243104215933593992;      //Euler-Mascheroni constant

double sqrt3 = sqrt(3.0);
double sqrt5 = sqrt(5.0);
double sqrt15 = sqrt3*sqrt5;

#endif //_Constants_Package_
