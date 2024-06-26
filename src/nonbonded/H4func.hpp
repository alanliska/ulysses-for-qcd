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

#ifndef _H4_Functions_
#define _H4_Functions_
#include <vector>
#include <string>
#include <math.h>
#include "../parameters/H4par.hpp"

//description:
//set of functions that calculate H4 contribution to intermolecular forces
//applicable to semi-empirical methods
//references:
//H4 -> J. Rezac, P. Hobza, J. Chem. Theory Comput., 8, 141, 2012

//functions related to H4 H-bond correction
double frad(double rDA) {
  //radial contribution to H4 correction
  double xx = rDA*rDA;
  double dfrad = -4.68512000000450434811 + 3.40736000001102778967*rDA + 5.31754666665572184314*xx - 7.20687407406838786983*xx*rDA + 3.25309629629461749545*xx*xx;
  dfrad += -0.70087111111082800452*xx*xx*rDA + 0.07357629629627092382*xx*xx*xx - 0.00303407407407313510*xx*xx*xx*rDA;
  return dfrad;
}
double gfrad(double rDA) {
  //gradient for radial contribution to H4 correction
  double xx = rDA*rDA;
  double dgfrad = 3.40736000001102778967 + 10.63509333331144368628*rDA - 21.62062222220516360949*xx + 13.01238518517846998179*xx*rDA - 3.50435555555413991158*xx*xx;
  dgfrad += 0.44145777777762551519*xx*xx*rDA - 0.02123851851851194655*xx*xx*xx;           
  return dgfrad;
}
double hfrad(double rDA) {
  //hessian for radial contribution to H4 correction
  double xx = rDA*rDA;
  double dhfrad = 10.63509333331144368628 - 43.24124444441033*rDA + 39.037155555535406*xx - 14.01742222221656*xx*rDA;
  dhfrad += 2.2072888888881277*xx*xx - 0.1274311111110717*xx*xx*rDA;           
  return dhfrad;
}
double fsw(double rDA) {
  //fsw auxiliary function
  double xxxx = rDA*rDA*rDA*rDA;
  double dfsw = 35.0*xxxx - 84.0*xxxx*rDA + 70.0*xxxx*rDA*rDA - 20.0*xxxx*rDA*rDA*rDA;
  return dfsw;
}
double gfsw(double rDA) {
  //gradient for fsw auxiliary function
  double xxx = rDA*rDA*rDA;
  double dgfsw = 140.0*xxx - 420.0*xxx*rDA + 420.0*xxx*rDA*rDA - 140.0*xxx*xxx;
  return dgfsw;
}
double hfsw(double rDA) {
  //hessian for fsw auxiliary function
  double xx = rDA*rDA;
  double dhfsw = 420.0*xx - 1680.0*xx*rDA + 2100.0*xx*xx - 840.0*xx*xx*rDA;
  return dhfsw;
}
double ValenceCorr(double rAB, size_t atomi, size_t atomj, double threshold = 1.0e-7) {
  //function for the valence correction for pair of atoms
  double ri = CovalentRadius(atomi);
  double rj = CovalentRadius(atomj);
  double r0 = ri + rj;
  double r1 = 1.6*r0;
  double vcorr = 0.0;
  if ((rAB > threshold)&&(rAB <= r0)) {vcorr = 1.0;}
  else if ((rAB > r0)&&(rAB < r1)) {vcorr = 1.0 - fsw((rAB - r0)/(0.6*r0));}
  return vcorr;
}
double gValenceCorr(double rAB, size_t atomi, size_t atomj, double threshold = 1.0e-7) {
  //function for the first-derivative of the valence correction for pair of atoms
  double ri = CovalentRadius(atomi);
  double rj = CovalentRadius(atomj);
  double r0 = ri + rj;
  double r1 = 1.6*r0;
  double gvcorr = 0.0;
  if ((rAB > r0)&&(rAB < r1)) {gvcorr = -gfsw((rAB - r0)/(0.6*r0))/(0.6*r0);}
  return gvcorr;
}
double hValenceCorr(double rAB, size_t atomi, size_t atomj, double threshold = 1.0e-7) {
  //function for the second-derivative of the valence correction for pair of atoms
  double ri = CovalentRadius(atomi);
  double rj = CovalentRadius(atomj);
  double r0 = ri + rj;
  double r1 = 1.6*r0;
  double hvcorr = 0.0;
  if ((rAB > r0)&&(rAB < r1)) {hvcorr = -hfsw((rAB - r0)/(0.6*r0))/(0.36*r0*r0);}
  return hvcorr;
}
//functions related to HH repulsion correction
double Erep(double rij, const std::string & method, double thresholdzero = 1.0e-8) {
  //function returning the repulsion contribution of a pair of protons in molecule
  double erep = R0HH(method);
  if (fabs(erep) < thresholdzero) {erep = rij;}
  double expfac = -eHH(method)*(rij - erep)/erep;
  erep = (1.0 - 1.0/(1.0 + exp(expfac)));
  return erep;
}
double gErep(double rij, const std::string & method, double thresholdzero = 1.0e-8) {
  //function returning the first-derivative of the repulsion contribution of a pair of protons in molecule
  double gerep = R0HH(method);
  double eHHfac = eHH(method);
  if (fabs(gerep) < thresholdzero) {gerep = rij;}
  double expterm = exp(-eHHfac*(rij - gerep)/gerep);
  eHHfac /= gerep;
  double denominator = (1.0 + expterm)*(1.0 + expterm);
  gerep = -eHHfac*expterm/denominator;
  return gerep;
}
double hErep(double rij, const std::string & method, double thresholdzero = 1.0e-8) {
  //function returning the second-derivative of the repulsion contribution of a pair of protons in molecule
  double herep = R0HH(method);
  double eHHfac = eHH(method);
  if (fabs(herep) < thresholdzero) {herep = rij;}
  double expterm = exp(-eHHfac*(rij - herep)/herep);
  eHHfac /= herep;
  double denominator = (1.0 + expterm)*(1.0 + expterm);
  herep = (eHHfac*eHHfac*expterm/denominator)*(1.0 - 2.0*expterm/(1.0 + expterm));
  return herep;
}
#endif //_H4_Functions_
