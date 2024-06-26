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

#ifndef _Dielectric_Constants_
#define _Dielectric_Constants_
#include <string>

//description:
//list of dielectric constants for solvation models

double SolventDiElectricConstant(std::string solvent, std::string method) {
  //list of dielectric constants taken from https://depts.washington.edu/eooptic/linkfiles/dielectric_chart%5B1%5D.pdf
  double epsilon = 1.0;
  if (solvent == "acetone") {epsilon = 20.7;}
  else if (solvent == "acetonitrile") {epsilon = 37.5;}
  else if (solvent == "aniline") {epsilon = 6.8882;}
  else if (solvent == "benzaldehyde") {epsilon = 18.22;}
  else if (solvent == "benzene") {
    if (method == "GFN2") {epsilon = 7.0;}
    else {epsilon = 2.27;}
  }
  else if (solvent == "dichloromethane") {epsilon = 8.93;}
  else if (solvent == "chloroform") {epsilon = 4.81;}
  else if (solvent == "carbon disulfide") {
    if (method == "GFN2") {epsilon = 2.64;}
    else {epsilon = 2.6;}
  }
  else if (solvent == "dioxane") {
    if (method == "GFN2") {epsilon = 2.2099;}
    else {epsilon = 2.25;}
  }
  else if (solvent == "dmf") {
    if (method == "GFN2") {epsilon = 37.0;}
    else {epsilon = 36.7;}
  }
  else if (solvent == "dmso") {epsilon = 46.68;}
  else if (solvent == "ethanol") {
    if (method == "GFN2") {epsilon = 25.3;}
    else {epsilon = 24.5;}
  }
  else if (solvent == "diethyl ether") {epsilon = 4.33;}
  else if (solvent == "ethyl acetate") {epsilon = 6.02;}
  else if (solvent == "furane") {epsilon = 2.94;}
  else if (solvent == "hexadecane") {
    if (method == "GFN2") {epsilon = 1.88;}
    else {epsilon = 2.09;}
  }
  else if (solvent == "hexane") {epsilon = 1.88;}
  else if (solvent == "methanol") {epsilon = 32.7;}
  else if (solvent == "nitromethane") {
    if (method == "GFN2") {epsilon = 36.562;}
    else {epsilon = 35.87;}
  }
  else if (solvent == "octanol") {epsilon = 9.8629;}
  else if (solvent == "phenol") {epsilon = 12.4;}
  else if (solvent == "thf") {
    if (method == "GFN2") {epsilon = 7.6;}
    else {epsilon = 7.58;}
  }
  else if (solvent == "toluene") {
    if (method == "GFN2") {epsilon = 7.0;}
    else {epsilon = 2.38;}
  }
  else if (solvent == "water") {
    if (method == "GFN2") {epsilon = 80.2;}
    else {epsilon = 80.1;}
  }
  else if (solvent == "octanol wet") {epsilon = 9.8629;}
  return epsilon;
}

#endif //_Dielectric_Constants_
