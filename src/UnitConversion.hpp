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

#ifndef _Unit_Conversion_
#define _Unit_Conversion_
#include <string>
#include "ConstantsPackage.hpp"

//description:
//unit conversions
//c.f. comments to know more

//CCCBDB
//double au2cmm1 = 219474.5984588712;                            //atomic units (Hartree) to {cm}^{-1}

//MOPAC
//double dist_Angstrom2au = 0.529167;                    //Angstrom to atomic unit of length (Bohr)
//double au2eV = 27.21;                                  //atomic units (Hartree) to electron Volt

//xTB
//double au2cmm1 = 219474.63067;
//double dist_Angstrom2au = 0.52917726;                    //Angstrom to atomic unit of length (Bohr)
//double au2eV = 27.21138505;                              //atomic units (Hartree) to electron Volt
//double au2kcalmol = 627.50947428;
//double fs2au = 41.3413733365614;                       //femtosecond to atomic unit of time

//DFTD3
//double dist_Angstrom2au = 0.52917726;
//double au2kcalmol = 627.509541;
//double au2eV = 27.21138505;

//H+
//double dist_Angstrom2au = 0.52917726;
//double au2kcalmol = 627.509541;

//pDynamo
//double dist_Angstrom2au = 0.5291772083;
//double au2eV = 27.2113845;

//ULYSSES parameters
//whenever possible from IUPAC gold book and NIST
double au2J = 4.3597438100364254e-18;                  //atomic units (Hartree) to Joule
double au2cmm1 = 219474.63;                            //atomic units (Hartree) to {cm}^{-1}
double au2eV = 27.211383473452294;                     //atomic units (Hartree) to electron Volt
double au2kcalmol = 627.5094198571855;                 //atomic units (Hartree) to kcal/mol
double dist_si2au = 5.291772083e-11;                   //meter to atomic unit of length (Bohr)
double dist_si2Angstrom = 1.0e-10;                     //meter to Angstrom
double dist_Angstrom2au = 0.5291772083;                //Angstrom to atomic unit of length (Bohr)
double D2eA = 0.208194;                                //Debye to electron.Angstrom 
double au2kg = Mel;                                    //converting atomic unit to kg
double fs2au = 41.341373336493;                        //femtosecond to atomic unit of time
double kg2melec = 1.0/Mel;                             //convert kg to mass of electron
double cal2J = 4.184;

double dist_Angstrom2aum1 = 1.0/dist_Angstrom2au;
double au2fs = 1.0/fs2au;
double cmm12au = 1.0/au2cmm1;

double Energy(double hartree, std::string units) {
  //source: http://meijer.group.shef.ac.uk/local/units.html 14.07.2018
  double result = 0;
  if (units == "eV") {result = hartree*au2eV;}
  else if (units == "J") {result = hartree*au2J;}
  else if (units == "J/mol") {result = hartree*2.6254999801646E+6;}
  else if (units == "Hz") {result = hartree*6.5796838976216E+15;}
  else if (units == "1/cm") {result = hartree*au2cmm1;}
  else if (units == "kcal/mol") {result = hartree*au2kcalmol;}
  return result;
}

#endif //_Unit_Conversion_
