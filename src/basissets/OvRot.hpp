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

#ifndef _Overlap_Rotations_
#define _Overlap_Rotations_
#include <vector>
//#include "/Users/FilipeMM/Documents/Filipe/Thermodynamics/math/MatrixPackage.hpp"
#include "../math/MatrixPackage.hpp"

//description:
//set of functions that make the rotations of overlap matrices to get the overlap functions in cartesian coordinates and molecular frame
//also contains first- and second-derivatives with respect to theta and phi
//note that composite transformations could be easily obtained using a common function that would receive pointers to functions
//for each type of transformation this would require two functions instead, if the orbitals to rotate are the same or not
//the code would indeed be slightly shorter, but I would have a hard-time incorporating the non-rotation part, which I require for CNDO/BG
//and since I spent the time and effort implementing it, better have something that also can be used for it

bool debug = false;

matrixE SPtransf(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that returns the rotation for integrals involving s and p orbitals
  matrixE transf(3,3);
  if (rotate) {
    transf(1,1) = st*cp;       //p(sigma)->px
    transf(2,1) = st*sp;       //p(sigma)->py
    transf(3,1) = ct;          //p(sigma)->pz
    transf(1,2) = ct*cp;       //p(pi)->px
    transf(2,2) = ct*sp;       //p(pi)->py
    transf(3,2) = -st;         //p(pi)->pz
    transf(1,3) = -sp;         //p(-pi)->px
    transf(2,3) = cp;          //p(-pi)->py
    transf(3,3) = 0.0;         //p(-pi)->pz
  }
  else {
    transf(2,1) = 1.0;         //p(sigma)
    transf(3,2) = 1.0;         //p(pi)
    transf(1,3) = 1.0;         //p(-pi)
  }
  return transf;
}
matrixE SPtransf_dt(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that returns the derivative of the rotation for integrals involving s and p orbitals
  //with respect to theta
  matrixE transf(3,3);
  if (rotate) {
    transf(1,1) = ct*cp;       //p(sigma)->px
    transf(2,1) = ct*sp;       //p(sigma)->py
    transf(3,1) = -st;         //p(sigma)->pz
    transf(1,2) = -st*cp;      //p(pi)->px
    transf(2,2) = -st*sp;      //p(pi)->py
    transf(3,2) = -ct;         //p(pi)->pz
    transf(1,3) = 0.0;         //p(-pi)->px
    transf(2,3) = 0.0;         //p(-pi)->py
    transf(3,3) = 0.0;         //p(-pi)->pz
  }
  else {
    transf(2,1) = 1.0;         //p(sigma)
    transf(3,2) = 1.0;         //p(pi)
    transf(1,3) = 1.0;         //p(-pi)
  }
  return transf;
}
matrixE SPtransf_dp(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that returns the derivative of the rotation for integrals involving s and p orbitals
  //with respect to phi
  matrixE transf(3,3);
  if (rotate) {
    transf(1,1) = -st*sp;      //p(sigma)->px
    transf(2,1) = st*cp;       //p(sigma)->py
    transf(3,1) = 0.0;         //p(sigma)->pz
    transf(1,2) = -ct*sp;      //p(pi)->px
    transf(2,2) = ct*cp;       //p(pi)->py
    transf(3,2) = 0.0;         //p(pi)->pz
    transf(1,3) = -cp;         //p(-pi)->px
    transf(2,3) = -sp;         //p(-pi)->py
    transf(3,3) = 0.0;         //p(-pi)->pz
  }
  else {
    transf(2,1) = 1.0;         //p(sigma)
    transf(3,2) = 1.0;         //p(pi)
    transf(1,3) = 1.0;         //p(-pi)
  }
  return transf;
}
matrixE SPtransf_dt2(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that returns the second-derivative of the rotation for integrals involving s and p orbitals
  //with respect to theta
  matrixE transf(3,3);
  if (rotate) {
    transf(1,1) = -st*cp;      //p(sigma)->px
    transf(2,1) = -st*sp;      //p(sigma)->py
    transf(3,1) = -ct;         //p(sigma)->pz
    transf(1,2) = -ct*cp;      //p(pi)->px
    transf(2,2) = -ct*sp;      //p(pi)->py
    transf(3,2) = st;          //p(pi)->pz
    transf(1,3) = 0.0;         //p(-pi)->px
    transf(2,3) = 0.0;         //p(-pi)->py
    transf(3,3) = 0.0;         //p(-pi)->pz
  }
  else {
    transf(2,1) = 1.0;         //p(sigma)
    transf(3,2) = 1.0;         //p(pi)
    transf(1,3) = 1.0;         //p(-pi)
  }
  return transf;
}
matrixE SPtransf_dtdp(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that returns the second-derivative of the rotation for integrals involving s and p orbitals
  //with respect to theta and phi
  matrixE transf(3,3);
  if (rotate) {
    transf(1,1) = -ct*sp;      //p(sigma)->px
    transf(2,1) = ct*cp;       //p(sigma)->py
    transf(3,1) = 0.0;         //p(sigma)->pz
    transf(1,2) = st*sp;       //p(pi)->px
    transf(2,2) = -st*cp;      //p(pi)->py
    transf(3,2) = 0.0;         //p(pi)->pz
    transf(1,3) = 0.0;         //p(-pi)->px
    transf(2,3) = 0.0;         //p(-pi)->py
    transf(3,3) = 0.0;         //p(-pi)->pz
  }
  else {
    transf(2,1) = 1.0;         //p(sigma)
    transf(3,2) = 1.0;         //p(pi)
    transf(1,3) = 1.0;         //p(-pi)
  }
  return transf;
}
matrixE SPtransf_dp2(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that returns the second-derivative of the rotation for integrals involving s and p orbitals
  //with respect to phi
  matrixE transf(3,3);
  if (rotate) {
    transf(1,1) = -st*cp;      //p(sigma)->px
    transf(2,1) = -st*sp;      //p(sigma)->py
    transf(3,1) = 0.0;         //p(sigma)->pz
    transf(1,2) = -ct*cp;      //p(pi)->px
    transf(2,2) = -ct*sp;      //p(pi)->py
    transf(3,2) = 0.0;         //p(pi)->pz
    transf(1,3) = sp;          //p(-pi)->px
    transf(2,3) = -cp;         //p(-pi)->py
    transf(3,3) = 0.0;         //p(-pi)->pz
  }
  else {
    transf(2,1) = 1.0;         //p(sigma)
    transf(3,2) = 1.0;         //p(pi)
    transf(1,3) = 1.0;         //p(-pi)
  }
  return transf;
}
matrixE SDtransf(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that returns the rotation for integrals involving s and d orbitals
  matrixE transf(5,5);
  double c2t = ct*ct - st*st;
  double c2p = cp*cp - sp*sp;
  double s2t = 2.0*st*ct;
  double s2p = 2.0*sp*cp;
  if (rotate) {
    transf(1,1) = sqrt(3.0)*ct*st*cp;             //d(sigma) -> dxz
    transf(2,1) = sqrt(3.0)*ct*st*sp;             //d(sigma) -> dyz
    transf(3,1) = 1.0 - 1.5*st*st;                //d(sigma) -> dz2
    transf(4,1) = sqrt(3.0)*st*st*sp*cp;          //d(sigma) -> dxy
    transf(5,1) = 0.5*sqrt(3.0)*st*st*c2p;        //d(sigma) -> dx2-y2
    transf(1,2) = cp*c2t;                         //d(pi) -> dxz
    transf(2,2) = sp*c2t;                         //d(pi) -> dyz
    transf(3,2) = -sqrt(3.0)*s2t/2.0;             //d(pi) -> dz2
    transf(4,2) = 0.5*s2t*s2p;                    //d(pi) -> dxy
    transf(5,2) = 0.5*c2p*s2t;                    //d(pi) -> dx2-y2
    transf(1,3) = -ct*sp;                         //d(-pi) -> dxz
    transf(2,3) = ct*cp;                          //d(-pi) -> dyz
    transf(3,3) = 0.0;                            //d(-pi) -> dz2
    transf(4,3) = st*c2p;                         //d(-pi) -> dxy
    transf(5,3) = -st*s2p;                        //d(-pi) -> dx2-y2
    transf(1,4) = -0.5*s2t*cp;                    //d(delta) -> dxz
    transf(2,4) = -0.5*s2t*sp;                    //d(delta) -> dyz
    transf(3,4) = sqrt(3.0)*st*st/2.0;            //d(delta) -> dz2
    transf(4,4) = s2p*(ct*ct + 0.5*st*st);        //d(delta) -> dxy
    transf(5,4) = c2p*(ct*ct + 0.5*st*st);        //d(delta) -> dx2-y2
    transf(1,5) = st*sp;                          //d(-delta) -> dxz
    transf(2,5) = -st*cp;                         //d(-delta) -> dyz
    transf(3,5) = 0.0;                            //d(-delta) -> dz2
    transf(4,5) = ct*c2p;                         //d(-delta) -> dxy
    transf(5,5) = -ct*s2p;                        //d(-delta) -> dx2-y2
  }
  else {
    transf(3,1) = 1.0;       //d(sigma)
    transf(4,2) = 1.0;       //d(pi)
    transf(2,3) = 1.0;       //d(-pi)
    transf(5,4) = 1.0;       //d(delta)
    transf(1,5) = 1.0;       //d(-delta)
  }
  return transf;
}
matrixE SDtransf_dt(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that returns the derivatives of rotations for integrals involving s and d orbitals
  //with respect to theta
  matrixE transf(5,5);
  double c2t = ct*ct - st*st;
  double c2p = cp*cp - sp*sp;
  double s2t = 2.0*st*ct;
  double s2p = 2.0*sp*cp;
  if (rotate) {
    transf(1,1) = sqrt(3.0)*c2t*cp;               //d(sigma) -> dxz
    transf(2,1) = sqrt(3.0)*c2t*sp;               //d(sigma) -> dyz
    transf(3,1) = -1.5*s2t;                       //d(sigma) -> dz2
    transf(4,1) = 0.5*sqrt(3.0)*s2t*s2p;          //d(sigma) -> dxy
    transf(5,1) = 0.5*sqrt(3.0)*s2t*c2p;          //d(sigma) -> dx2-y2
    transf(1,2) = -2.0*cp*s2t;                    //d(pi) -> dxz
    transf(2,2) = -2.0*sp*s2t;                    //d(pi) -> dyz
    transf(3,2) = -sqrt(3.0)*c2t;                 //d(pi) -> dz2
    transf(4,2) = c2t*s2p;                        //d(pi) -> dxy
    transf(5,2) = c2p*c2t;                        //d(pi) -> dx2-y2
    transf(1,3) = st*sp;                          //d(-pi) -> dxz
    transf(2,3) = -st*cp;                         //d(-pi) -> dyz
    transf(3,3) = 0.0;                            //d(-pi) -> dz2
    transf(4,3) = ct*c2p;                         //d(-pi) -> dxy
    transf(5,3) = -ct*s2p;                        //d(-pi) -> dx2-y2
    transf(1,4) = -c2t*cp;                        //d(delta) -> dxz
    transf(2,4) = -c2t*sp;                        //d(delta) -> dyz
    transf(3,4) = 0.5*sqrt(3.0)*s2t;              //d(delta) -> dz2
    transf(4,4) = -0.5*s2p*s2t;                   //d(delta) -> dxy
    transf(5,4) = -0.5*c2p*s2t;                   //d(delta) -> dx2-y2
    transf(1,5) = ct*sp;                          //d(-delta) -> dxz
    transf(2,5) = -ct*cp;                         //d(-delta) -> dyz
    transf(3,5) = 0.0;                            //d(-delta) -> dz2
    transf(4,5) = -st*c2p;                        //d(-delta) -> dxy
    transf(5,5) = st*s2p;                         //d(-delta) -> dx2-y2
  }
  else {
    transf(3,1) = 1.0;       //d(sigma)
    transf(4,2) = 1.0;       //d(pi)
    transf(2,3) = 1.0;       //d(-pi)
    transf(5,4) = 1.0;       //d(delta)
    transf(1,5) = 1.0;       //d(-delta)
  }
  return transf;
}
matrixE SDtransf_dp(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that returns the derivative of rotations for integrals involving s and d orbitals
  //with respect to phi
  matrixE transf(5,5);
  double c2t = ct*ct - st*st;
  double c2p = cp*cp - sp*sp;
  double s2t = 2.0*st*ct;
  double s2p = 2.0*sp*cp;
  if (rotate) {
    transf(1,1) = -0.5*sqrt(3.0)*s2t*sp;          //d(sigma) -> dxz
    transf(2,1) = 0.5*sqrt(3.0)*s2t*cp;           //d(sigma) -> dyz
    transf(3,1) = 0.0;                            //d(sigma) -> dz2
    transf(4,1) = sqrt(3.0)*st*st*c2p;            //d(sigma) -> dxy
    transf(5,1) = -sqrt(3.0)*st*st*s2p;           //d(sigma) -> dx2-y2
    transf(1,2) = -sp*c2t;                        //d(pi) -> dxz
    transf(2,2) = cp*c2t;                         //d(pi) -> dyz
    transf(3,2) = 0.0;                            //d(pi) -> dz2
    transf(4,2) = s2t*c2p;                        //d(pi) -> dxy
    transf(5,2) = -s2p*s2t;                       //d(pi) -> dx2-y2
    transf(1,3) = -ct*cp;                         //d(-pi) -> dxz
    transf(2,3) = -ct*sp;                         //d(-pi) -> dyz
    transf(3,3) = 0.0;                            //d(-pi) -> dz2
    transf(4,3) = -2.0*st*s2p;                    //d(-pi) -> dxy
    transf(5,3) = -2.0*st*c2p;                    //d(-pi) -> dx2-y2
    transf(1,4) = 0.5*s2t*sp;                     //d(delta) -> dxz
    transf(2,4) = -0.5*s2t*cp;                    //d(delta) -> dyz
    transf(3,4) = 0.0;                            //d(delta) -> dz2
    transf(4,4) = c2p*(2.0*ct*ct + st*st);        //d(delta) -> dxy
    transf(5,4) = -s2p*(2.0*ct*ct + st*st);       //d(delta) -> dx2-y2
    transf(1,5) = st*cp;                          //d(-delta) -> dxz
    transf(2,5) = st*sp;                          //d(-delta) -> dyz
    transf(3,5) = 0.0;                            //d(-delta) -> dz2
    transf(4,5) = -2.0*ct*s2p;                    //d(-delta) -> dxy
    transf(5,5) = -2.0*ct*c2p;                    //d(-delta) -> dx2-y2
  }
  else {
    transf(3,1) = 1.0;       //d(sigma)
    transf(4,2) = 1.0;       //d(pi)
    transf(2,3) = 1.0;       //d(-pi)
    transf(5,4) = 1.0;       //d(delta)
    transf(1,5) = 1.0;       //d(-delta)
  }
  return transf;
}
matrixE SDtransf_dt2(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that returns the second-derivatives of rotations for integrals involving s and d orbitals
  //with respect to theta
  matrixE transf(5,5);
  double c2t = ct*ct - st*st;
  double c2p = cp*cp - sp*sp;
  double s2t = 2.0*st*ct;
  double s2p = 2.0*sp*cp;
  if (rotate) {
    transf(1,1) = -2.0*sqrt(3.0)*s2t*cp;          //d(sigma) -> dxz
    transf(2,1) = -2.0*sqrt(3.0)*s2t*sp;          //d(sigma) -> dyz
    transf(3,1) = -3.0*c2t;                       //d(sigma) -> dz2
    transf(4,1) = sqrt(3.0)*c2t*s2p;              //d(sigma) -> dxy
    transf(5,1) = sqrt(3.0)*c2t*c2p;              //d(sigma) -> dx2-y2
    transf(1,2) = -4.0*cp*c2t;                    //d(pi) -> dxz
    transf(2,2) = -4.0*sp*c2t;                    //d(pi) -> dyz
    transf(3,2) = 2.0*sqrt(3.0)*s2t;              //d(pi) -> dz2
    transf(4,2) = -2.0*s2t*s2p;                   //d(pi) -> dxy
    transf(5,2) = -2.0*s2t*c2p;                   //d(pi) -> dx2-y2
    transf(1,3) = ct*sp;                          //d(-pi) -> dxz
    transf(2,3) = -ct*cp;                         //d(-pi) -> dyz
    transf(3,3) = 0.0;                            //d(-pi) -> dz2
    transf(4,3) = -st*c2p;                        //d(-pi) -> dxy
    transf(5,3) = st*s2p;                         //d(-pi) -> dx2-y2
    transf(1,4) = 2.0*s2t*cp;                     //d(delta) -> dxz
    transf(2,4) = 2.0*s2t*sp;                     //d(delta) -> dyz
    transf(3,4) = sqrt(3.0)*c2t;                  //d(delta) -> dz2
    transf(4,4) = -s2p*c2t;                       //d(delta) -> dxy
    transf(5,4) = -c2p*c2t;                       //d(delta) -> dx2-y2
    transf(1,5) = -st*sp;                         //d(-delta) -> dxz
    transf(2,5) = st*cp;                          //d(-delta) -> dyz
    transf(3,5) = 0.0;                            //d(-delta) -> dz2
    transf(4,5) = -ct*c2p;                        //d(-delta) -> dxy
    transf(5,5) = ct*s2p;                         //d(-delta) -> dx2-y2
  }
  else {
    transf(3,1) = 1.0;       //d(sigma)
    transf(4,2) = 1.0;       //d(pi)
    transf(2,3) = 1.0;       //d(-pi)
    transf(5,4) = 1.0;       //d(delta)
    transf(1,5) = 1.0;       //d(-delta)
  }
  return transf;
}
matrixE SDtransf_dtdp(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that returns the second-derivatives of rotations for integrals involving s and d orbitals
  //with respect to theta and phi
  matrixE transf(5,5);
  double c2t = ct*ct - st*st;
  double c2p = cp*cp - sp*sp;
  double s2t = 2.0*st*ct;
  double s2p = 2.0*sp*cp;
  if (rotate) {
    transf(1,1) = -sqrt(3.0)*c2t*sp;              //d(sigma) -> dxz
    transf(2,1) = sqrt(3.0)*c2t*cp;               //d(sigma) -> dyz
    transf(3,1) = 0.0;                            //d(sigma) -> dz2
    transf(4,1) = sqrt(3.0)*s2t*c2p;              //d(sigma) -> dxy
    transf(5,1) = -sqrt(3.0)*s2t*s2p;             //d(sigma) -> dx2-y2
    transf(1,2) = 2.0*sp*s2t;                     //d(pi) -> dxz
    transf(2,2) = -2.0*cp*s2t;                    //d(pi) -> dyz
    transf(3,2) = 0.0;                            //d(pi) -> dz2
    transf(4,2) = 2.0*c2t*c2p;                    //d(pi) -> dxy
    transf(5,2) = -2.0*c2t*s2p;                   //d(pi) -> dx2-y2
    transf(1,3) = st*cp;                          //d(-pi) -> dxz
    transf(2,3) = st*sp;                          //d(-pi) -> dyz
    transf(3,3) = 0.0;                            //d(-pi) -> dz2
    transf(4,3) = -2.0*ct*s2p;                    //d(-pi) -> dxy
    transf(5,3) = -2.0*ct*c2p;                    //d(-pi) -> dx2-y2
    transf(1,4) = c2t*sp;                         //d(delta) -> dxz
    transf(2,4) = -c2t*cp;                        //d(delta) -> dyz
    transf(3,4) = 0.0;                            //d(delta) -> dz2
    transf(4,4) = -c2p*s2t;                       //d(delta) -> dxy
    transf(5,4) = s2p*s2t;                        //d(delta) -> dx2-y2
    transf(1,5) = ct*cp;                          //d(-delta) -> dxz
    transf(2,5) = ct*sp;                          //d(-delta) -> dyz
    transf(3,5) = 0.0;                            //d(-delta) -> dz2
    transf(4,5) = 2.0*st*s2p;                     //d(-delta) -> dxy
    transf(5,5) = 2.0*st*c2p;                     //d(-delta) -> dx2-y2
  }
  else {
    transf(3,1) = 1.0;       //d(sigma)
    transf(4,2) = 1.0;       //d(pi)
    transf(2,3) = 1.0;       //d(-pi)
    transf(5,4) = 1.0;       //d(delta)
    transf(1,5) = 1.0;       //d(-delta)
  }
  return transf;
}
matrixE SDtransf_dp2(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that returns the second-derivative of rotations for integrals involving s and d orbitals
  //with respect to phi
  matrixE transf(5,5);
  double c2t = ct*ct - st*st;
  double c2p = cp*cp - sp*sp;
  double s2t = 2.0*st*ct;
  double s2p = 2.0*sp*cp;
  if (rotate) {
    transf(1,1) = -0.5*sqrt(3.0)*s2t*cp;          //d(sigma) -> dxz
    transf(2,1) = -0.5*sqrt(3.0)*s2t*sp;          //d(sigma) -> dyz
    transf(3,1) = 0.0;                            //d(sigma) -> dz2
    transf(4,1) = -2.0*sqrt(3.0)*st*st*s2p;       //d(sigma) -> dxy
    transf(5,1) = -2.0*sqrt(3.0)*st*st*c2p;       //d(sigma) -> dx2-y2
    transf(1,2) = -cp*c2t;                        //d(pi) -> dxz
    transf(2,2) = -sp*c2t;                        //d(pi) -> dyz
    transf(3,2) = 0.0;                            //d(pi) -> dz2
    transf(4,2) = -2.0*s2t*s2p;                   //d(pi) -> dxy
    transf(5,2) = -2.0*c2p*s2t;                   //d(pi) -> dx2-y2
    transf(1,3) = ct*sp;                          //d(-pi) -> dxz
    transf(2,3) = -ct*cp;                         //d(-pi) -> dyz
    transf(3,3) = 0.0;                            //d(-pi) -> dz2
    transf(4,3) = -4.0*st*c2p;                    //d(-pi) -> dxy
    transf(5,3) = 4.0*st*s2p;                     //d(-pi) -> dx2-y2
    transf(1,4) = 0.5*s2t*cp;                     //d(delta) -> dxz
    transf(2,4) = 0.5*s2t*sp;                     //d(delta) -> dyz
    transf(3,4) = 0.0;                            //d(delta) -> dz2
    transf(4,4) = -2.0*s2p*(2.0*ct*ct + st*st);   //d(delta) -> dxy
    transf(5,4) = -2.0*c2p*(2.0*ct*ct + st*st);   //d(delta) -> dx2-y2
    transf(1,5) = -st*sp;                         //d(-delta) -> dxz
    transf(2,5) = st*cp;                          //d(-delta) -> dyz
    transf(3,5) = 0.0;                            //d(-delta) -> dz2
    transf(4,5) = -4.0*ct*c2p;                    //d(-delta) -> dxy
    transf(5,5) = 4.0*ct*s2p;                     //d(-delta) -> dx2-y2
  }
  else {
    transf(3,1) = 1.0;       //d(sigma)
    transf(4,2) = 1.0;       //d(pi)
    transf(2,3) = 1.0;       //d(-pi)
    transf(5,4) = 1.0;       //d(delta)
    transf(1,5) = 1.0;       //d(-delta)
  }
  return transf;
}
matrixE SFtransf(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that returns the rotation for integrals involving s and f orbitals
  matrixE transf(7,7);
  double cos2t = cost*cost - sint*sint;
  double cos2p = cosp*cosp - sinp*sinp;
  double sin2t = 2.0*sint*cost;
  double sin2p = 2.0*sinp*cosp;
  double cos3t = cos2t*cost - sin2t*sint;
  double cos3p = cos2p*cosp - sin2p*sinp;
  double sin3t = cos2t*sint + sin2t*cost;
  double sin3p = cos2p*sinp + sin2p*cosp;
  if (rotate) {
    //order of cartesian orbitals: fz(x2-y2); fx(x2-3y2); fy(3x2-y2); fxyz; fz3; fxz2; fyz2
    //sigma
    transf(1,1) = 0.5*cost*(5*cost*cost - 3.0);
    transf(2,1) = sqrt(6.0)*0.25*cosp*sint*(5*cost*cost - 1.0);
    transf(3,1) = sqrt(6.0)*0.25*sinp*sint*(5*cost*cost - 1.0);
    transf(4,1) = sqrt(15.0)*sinp*cosp*cost*sint*sint;
    transf(5,1) = sqrt(15.0)*0.5*cos2p*cost*sint*sint;
    transf(6,1) = sqrt(10.0)*0.25*cos3p*sint*sint*sint;
    transf(7,1) = sqrt(10.0)*0.25*sin3p*sint*sint*sint;
    //pi
    transf(1,1) = -sqrt(6.0)*0.25*sint*(5.0*cost*cost - 1.0);
    transf(2,2) = 0.25*cosp*cost*(15.0*cost*cost - 11.0);
    transf(3,2) = 0.25*sinp*cost*(15.0*cost*cost - 11.0);
    transf(4,2) = sqrt(10.0)*0.25*sin2p*sint*(3.0*cost*cost - 1.0);
    transf(5,2) = sqrt(10.0)*0.25*cos2p*sint*(3.0*cost*cost - 1.0);
    transf(6,2) = sqrt(15.0)*0.25*cos3p*cost*sint*sint;
    transf(7,2) = sqrt(15.0)*0.25*sin3p*cost*sint*sint;
    //-pi
    transf(1,3) = 0.0;
    transf(2,3) = -0.25*sinp*(5.0*cost*cost - 1.0);
    transf(3,3) = 0.25*cosp*(5.0*cost*cost - 1.0);
    transf(4,3) = sqrt(10.0)*0.5*cos2p*sint*cost;
    transf(5,3) = -sqrt(10.0)*sint*cost*sinp*cosp;
    transf(6,3) = -sqrt(15.0)*0.25*sin3p*sint*sint;
    transf(7,3) = sqrt(15.0)*0.25*sint*sint*cosp*(4.0*cosp*cosp - 3.0);
    //delta
    transf(1,4) = sqrt(15.0)*0.5*cost*sint*sint;
    transf(2,4) = sqrt(10.0)*0.25*cosp*sint*(1.0 - 3.0*cost*cost);
    transf(3,4) = sqrt(10.0)*0.25*sint*sinp*(1.0 - 3.0*cost*cost);
    transf(4,4) = sinp*cosp*cost*(3.0*cost*cost - 1.0);
    transf(5,4) = 0.5*cos2p*cost*(3.0*cost*cost - 1.0);
    transf(6,4) = sqrt(6.0)*0.25*cos3p*sint*(1.0 + cost*cost);
    transf(7,4) = sqrt(6.0)*0.25*sin3p*sint*(1.0 + cost*cost);
    //-delta
    transf(1,5) = 0.0;
    transf(2,5) = sqrt(10.0)*0.5*sinp*cost*sint;
    transf(3,5) = -sqrt(10.0)*0.5*cosp*cost*sint;
    transf(4,5) = cos2p*cos2t;
    transf(5,5) = -sin2p*cos2t;
    transf(6,5) = -sqrt(6.0)*0.5*sin3p*sint*cost;
    transf(7,5) = sqrt(6.0)*0.5*cost*sint*cosp*(4.0*cosp*cosp - 3.0);
    //phi
    transf(1,6) = sqrt(10.0)*0.25*sint*sint*sint;
    transf(2,6) = sqrt(15.0)*0.25*cosp*cost*sint*sint;
    transf(3,6) = sqrt(15.0)*0.25*sinp*cost*sint*sint;
    transf(4,6) = -sqrt(6.0)*0.5*sinp*cosp*sint*(1.0 + cost*cost);
    transf(5,6) = -sqrt(6.0)*0.25*cos2p*sint*(1.0 + cost*cost);
    transf(6,6) = 0.25*cos3p*cost*(3.0 + cost*cost);
    transf(7,6) = 0.25*sin3p*cost*(3.0 + cost*cost);
    //-phi
    transf(1,7) = 0.0;
    transf(2,7) = -sqrt(15.0)*0.25*sinp*sint*sint;
    transf(3,7) = sqrt(15.0)*0.25*cosp*sint*sint;
    transf(4,7) = -sqrt(6.0)*0.5*cos2p*sint*cost;
    transf(5,7) = sqrt(6.0)*sinp*cosp*sint*cost;
    transf(6,7) = -0.25*sin3p*(1.0 + 3.0*cost*cost);
    transf(7,7) = 0.25*cos3p*(1.0 + 3.0*cost*cost);
  }
  else {
    transf(4,1) = 1.0;      //f(sigma)
  }
  return transf;
}
matrixE SFtransf_dt(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that returns the derivatives of rotations for integrals involving s and f orbitals
  //with respect to theta
  matrixE transf(7,7);
  double cos2t = cost*cost - sint*sint;
  double cos2p = cosp*cosp - sinp*sinp;
  double sin2t = 2.0*sint*cost;
  double sin2p = 2.0*sinp*cosp;
  double cos3t = cos2t*cost - sin2t*sint;
  double cos3p = cos2p*cosp - sin2p*sinp;
  double sin3t = cos2t*sint + sin2t*cost;
  double sin3p = cos2p*sinp + sin2p*cosp;
  if (rotate) {
    //order of cartesian orbitals: fz(x2-y2); fx(x2-3y2); fy(3x2-y2); fxyz; fz3; fxz2; fyz2
    //sigma
    transf(1,1) = -1.5*sint*(5*cost*cost - 1.0);
    transf(2,1) = sqrt(6.0)*0.25*cosp*cost*(5*cos2t - 1.0 - 5.0*sint*sint);
    transf(3,1) = sqrt(6.0)*0.25*sinp*cost*(5*cos2t - 1.0 - 5.0*sint*sint);
    transf(4,1) = sqrt(15.0)*0.25*sin2p*(cost*sin2t + 2.0*sint*cos2t);
    transf(5,1) = sqrt(15.0)*0.25*cos2p*(cost*sin2t + 2.0*sint*cos2t);
    transf(6,1) = sqrt(10.0)*0.75*cos3p*sint*sint*cost;
    transf(7,1) = sqrt(10.0)*0.75*sin3p*sint*sint*cost;
    //pi
    transf(1,1) = -sqrt(6.0)*0.25*cost*(5.0*cos2t - 1.0 - 5.0*sint*sint);
    transf(2,2) = -0.25*cosp*sint*(45.0*cost*cost - 11.0);
    transf(3,2) = -0.25*sinp*sint*(45.0*cost*cost - 11.0);
    transf(4,2) = sqrt(10.0)*0.25*sin2p*cost*(3.0*cos2t - 1.0 - 3.0*sint*sint);
    transf(5,2) = sqrt(10.0)*0.25*cos2p*cost*(3.0*cos2t - 1.0 - 3.0*sint*sint);
    transf(6,2) = sqrt(15.0)*0.25*cos3p*sint*(cost*cost + cos2t);
    transf(7,2) = sqrt(15.0)*0.25*sin3p*sint*(cost*cost + cos2t);
    //-pi
    transf(1,3) = 0.0;
    transf(2,3) = 2.5*sinp*sint*cost;
    transf(3,3) = -2.5*cosp*sint*cost;
    transf(4,3) = sqrt(10.0)*0.5*cos2p*cos2t;
    transf(5,3) = -sqrt(10.0)*0.5*sin2p*cos2t;
    transf(6,3) = -sqrt(15.0)*0.25*sin3p*sin2t;
    transf(7,3) = sqrt(15.0)*0.25*sin2t*cosp*(4.0*cosp*cosp - 3.0);
    //delta
    transf(1,4) = sqrt(15.0)*0.5*sint*(cost*cost + cos2t);
    transf(2,4) = sqrt(10.0)*0.25*cosp*cost*(1.0 + 3.0*sint*sint - 3.0*cos2t);
    transf(3,4) = sqrt(10.0)*0.25*sinp*cost*(1.0 + 3.0*sint*sint - 3.0*cos2t);
    transf(4,4) = -0.5*sin2p*sint*(9.0*cost*cost - 1.0);
    transf(5,4) = -0.5*cos2p*sint*(9.0*cost*cost - 1.0);
    transf(6,4) = sqrt(6.0)*0.25*cos3p*cost*(1.0 + cos2t - sint*sint);
    transf(7,4) = sqrt(6.0)*0.25*sin3p*cost*(1.0 + cos2t - sint*sint);
    //-delta
    transf(1,5) = 0.0;
    transf(2,5) = sqrt(10.0)*0.5*sinp*cos2t;
    transf(3,5) = -sqrt(10.0)*0.5*cosp*cos2t;
    transf(4,5) = -2.0*cos2p*sin2t;
    transf(5,5) = 2.0*sin2p*sin2t;
    transf(6,5) = -sqrt(6.0)*0.5*sin3p*cos2t;
    transf(7,5) = sqrt(6.0)*0.5*cos3p*cos2t;
    //phi
    transf(1,6) = sqrt(10.0)*0.75*sint*sint*cost;
    transf(2,6) = sqrt(15.0)*0.25*cosp*sint*(cost*cost + cos2t);
    transf(3,6) = sqrt(15.0)*0.25*sinp*sint*(cost*cost + cos2t);
    transf(4,6) = -sqrt(6.0)*0.25*sin2p*cost*(1.0 + cos2t - sint*sint);
    transf(5,6) = -sqrt(6.0)*0.25*cos2p*cost*(1.0 + cos2t - sint*sint);
    transf(6,6) = -0.75*cos3p*sint*(1.0 + cost*cost);
    transf(7,6) = -0.75*sin3p*sint*(1.0 + cost*cost);
    //-phi
    transf(1,7) = 0.0;
    transf(2,7) = -sqrt(15.0)*0.5*sinp*sint*cost;
    transf(3,7) = sqrt(15.0)*0.5*cosp*sint*cost;
    transf(4,7) = -sqrt(6.0)*0.5*cos2p*cos2t;
    transf(5,7) = sqrt(6.0)*0.5*sin2p*cos2t;
    transf(6,7) = 0.75*sin3p*sin2t;
    transf(7,7) = -0.75*cos3p*sin2t;
  }
  else {
    transf(4,1) = 1.0;      //f(sigma)
  }
  return transf;
}
matrixE SFtransf_dp(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that returns the derivatives of rotations for integrals involving s and f orbitals
  //with respect to phi
  matrixE transf(7,7);
  double cos2t = cost*cost - sint*sint;
  double cos2p = cosp*cosp - sinp*sinp;
  double sin2t = 2.0*sint*cost;
  double sin2p = 2.0*sinp*cosp;
  double cos3t = cos2t*cost - sin2t*sint;
  double cos3p = cos2p*cosp - sin2p*sinp;
  double sin3t = cos2t*sint + sin2t*cost;
  double sin3p = cos2p*sinp + sin2p*cosp;
  if (rotate) {
    //order of cartesian orbitals: fz(x2-y2); fx(x2-3y2); fy(3x2-y2); fxyz; fz3; fxz2; fyz2
    //sigma
    transf(1,1) = 0.0;
    transf(2,1) = -sqrt(6.0)*0.25*sinp*sint*(5*cost*cost - 1.0);
    transf(3,1) = sqrt(6.0)*0.25*cosp*sint*(5*cost*cost - 1.0);
    transf(4,1) = sqrt(15.0)*0.5*sint*cos2p*sin2t;
    transf(5,1) = -sqrt(15.0)*0.5*sint*sin2p*sin2t;
    transf(6,1) = -sqrt(10.0)*0.75*sin3p*sint*sint*sint;
    transf(7,1) = sqrt(10.0)*0.75*cos3p*sint*sint*sint;
    //pi
    transf(1,1) = 0.0;
    transf(2,2) = -0.25*sinp*cost*(15.0*cost*cost - 11.0);
    transf(3,2) = 0.25*cosp*cost*(15.0*cost*cost - 11.0);
    transf(4,2) = sqrt(10.0)*0.5*cos2p*sint*(3.0*cost*cost - 1.0);
    transf(5,2) = -sqrt(10.0)*0.5*sin2p*sint*(3.0*cost*cost - 1.0);
    transf(6,2) = -sqrt(15.0)*0.75*sin3p*cost*sint*sint;
    transf(7,2) = sqrt(15.0)*0.75*cos3p*cost*sint*sint;
    //-pi
    transf(1,3) = 0.0;
    transf(2,3) = -0.25*cosp*(5.0*cost*cost - 1.0);
    transf(3,3) = -0.25*sinp*(5.0*cost*cost - 1.0);
    transf(4,3) = -sqrt(10.0)*0.5*sin2p*sin2t;
    transf(5,3) = -sqrt(10.0)*0.5*sin2t*cos2p;
    transf(6,3) = -sqrt(15.0)*0.75*cos3p*sint*sint;
    transf(7,3) = -sqrt(15.0)*0.75*sint*sint*sinp*(4.0*cosp*cosp - 1.0);
    //delta
    transf(1,4) = 0.0;
    transf(2,4) = sqrt(10.0)*0.25*sinp*sint*(3.0*cost*cost - 1.0);
    transf(3,4) = -sqrt(10.0)*0.25*cosp*sint*(3.0*cost*cost - 1.0);
    transf(4,4) = cos2p*cost*(3.0*cost*cost - 1.0);
    transf(5,4) = -sin2p*cost*(3.0*cost*cost - 1.0);
    transf(6,4) = -sqrt(6.0)*0.75*sin3p*sint*(1.0 + cost*cost);
    transf(7,4) = sqrt(6.0)*0.75*cos3p*sint*(1.0 + cost*cost);
    //-delta
    transf(1,5) = 0.0;
    transf(2,5) = sqrt(10.0)*0.25*cosp*sin2t;
    transf(3,5) = sqrt(10.0)*0.25*sinp*sin2t;
    transf(4,5) = -2.0*sin2p*cos2t;
    transf(5,5) = -2.0*cos2p*cos2t;
    transf(6,5) = -sqrt(6.0)*0.75*cos3p*sin2t;
    transf(7,5) = -sqrt(6.0)*0.75*sin2t*sin3p;
    //phi
    transf(1,6) = 0.0;
    transf(2,6) = -sqrt(15.0)*0.25*sinp*cost*sint*sint;
    transf(3,6) = sqrt(15.0)*0.25*cosp*cost*sint*sint;
    transf(4,6) = -sqrt(6.0)*0.5*cos2p*sint*(1.0 + cost*cost);
    transf(5,6) = sqrt(6.0)*0.5*sin2p*sint*(1.0 + cost*cost);
    transf(6,6) = -0.75*sin3p*cost*(3.0 + cost*cost);
    transf(7,6) = 0.75*cos3p*cost*(3.0 + cost*cost);
    //-phi
    transf(1,7) = 0.0;
    transf(2,7) = -sqrt(15.0)*0.25*cosp*sint*sint;
    transf(3,7) = -sqrt(15.0)*0.25*sinp*sint*sint;
    transf(4,7) = sqrt(6.0)*0.5*sin2p*sin2t;
    transf(5,7) = sqrt(6.0)*0.5*cos2p*sin2t;
    transf(6,7) = -0.75*cos3p*(1.0 + 3.0*cost*cost);
    transf(7,7) = -0.75*sin3p*(1.0 + 3.0*cost*cost);
  }
  else {
    transf(4,1) = 1.0;      //f(sigma)
  }
  return transf;
}
matrixE SFtransf_dt2(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that returns the second-derivatives of rotations for integrals involving s and f orbitals
  //with respect to theta
  matrixE transf(7,7);
  double cos2t = cost*cost - sint*sint;
  double cos2p = cosp*cosp - sinp*sinp;
  double sin2t = 2.0*sint*cost;
  double sin2p = 2.0*sinp*cosp;
  double cos3t = cos2t*cost - sin2t*sint;
  double cos3p = cos2p*cosp - sin2p*sinp;
  double sin3t = cos2t*sint + sin2t*cost;
  double sin3p = cos2p*sinp + sin2p*cosp;
  if (rotate) {
    //order of cartesian orbitals: fz(x2-y2); fx(x2-3y2); fy(3x2-y2); fxyz; fz3; fxz2; fyz2
    //sigma
    transf(1,1) = 15.0*cost*(sint*sint - 0.5*cost*cost) + 1.5*(cost + sint);
    transf(2,1) = sqrt(6.0)*0.25*cosp*(sint*(1.0 + 10.0*sint*sint - 5.0*cost*cost) - 15.0*cost*sin2t);
    transf(3,1) = sqrt(6.0)*0.25*sinp*(sint*(1.0 + 10.0*sint*sint - 5.0*cost*cost) - 15.0*cost*sin2t);
    transf(4,1) = sqrt(15.0)*0.25*sin2p*(4.0*cost*cos2t - 5.0*sint*sin2t);
    transf(5,1) = sqrt(15.0)*0.25*cos2p*(4.0*cost*cos2t - 5.0*sint*sin2t);
    transf(6,1) = sqrt(10.0)*0.75*cos3p*sint*(2.0*cost*cost - sint*sint);
    transf(7,1) = sqrt(10.0)*0.75*sin3p*sint*(2.0*cost*cost - sint*sint);
    //pi
    transf(1,1) = sqrt(6.0)*0.25*(sint*(5.0*cost*cost - 1.0 - 10.0*sint*sint) + 15.0*cost*sin2t);
    transf(2,2) = 0.25*cosp*cost*(90.0*sint*sint - 45.0*cost*cost + 11.0);
    transf(3,2) = 0.25*sinp*cost*(90.0*sint*sint - 45.0*cost*cost + 11.0);
    transf(4,2) = sqrt(10.0)*0.25*sin2p*sint*(-21.0*cost*cost + 1.0 + 6.0*sint*sint);
    transf(5,2) = sqrt(10.0)*0.25*cos2p*sint*(-21.0*cost*cost + 1.0 + 6.0*sint*sint);
    transf(6,2) = sqrt(15.0)*0.25*cos3p*cost*(2.0*cost*cost - 7.0*sint*sint);
    transf(7,2) = sqrt(15.0)*0.25*sin3p*cost*(2.0*cost*cost - 7.0*sint*sint);
    //-pi
    transf(1,3) = 0.0;
    transf(2,3) = 2.5*sinp*cos2t;
    transf(3,3) = -2.5*cosp*cos2t;
    transf(4,3) = -sqrt(10.0)*cos2p*sin2t;
    transf(5,3) = sqrt(10.0)*sin2p*sin2t;
    transf(6,3) = -sqrt(15.0)*0.5*sin3p*cos2t;
    transf(7,3) = sqrt(15.0)*0.5*cos2t*cosp*(4.0*cosp*cosp - 3.0);
    //delta
    transf(1,4) = sqrt(15.0)*0.5*cost*(2.0*cost*cost - 7.0*sint*sint);
    transf(2,4) = -sqrt(10.0)*0.25*cosp*sint*(1.0 + 6.0*sint*sint - 21.0*cost*cost);
    transf(3,4) = -sqrt(10.0)*0.25*sinp*sint*(1.0 + 6.0*sint*sint - 21.0*cost*cost);
    transf(4,4) = -0.5*sin2p*cost*(9.0*cos2t - 1.0 - 9.0*sint*sint);
    transf(5,4) = -0.5*cos2p*cost*(9.0*cos2t - 1.0 - 9.0*sint*sint);
    transf(6,4) = -sqrt(6.0)*0.25*cos3p*sint*(1.0 + 7.0*cost*cost - 2.0*sint*sint);
    transf(7,4) = -sqrt(6.0)*0.25*sin3p*sint*(1.0 + 7.0*cost*cost - 2.0*sint*sint);
    //-delta
    transf(1,5) = 0.0;
    transf(2,5) = -sqrt(10.0)*sinp*sin2t;
    transf(3,5) = sqrt(10.0)*cosp*sin2t;
    transf(4,5) = -4.0*cos2p*cos2t;
    transf(5,5) = 4.0*sin2p*cos2t;
    transf(6,5) = sqrt(6.0)*sin3p*sin2t;
    transf(7,5) = -sqrt(6.0)*cos3p*sin2t;
    //phi
    transf(1,6) = sqrt(10.0)*0.75*sint*(2.0*cost*cost - sint*sint);
    transf(2,6) = sqrt(15.0)*0.25*cosp*cost*(2.0*cost*cost - 7.0*sint*sint);
    transf(3,6) = sqrt(15.0)*0.25*sinp*cost*(2.0*cost*cost - 7.0*sint*sint);
    transf(4,6) = sqrt(6.0)*0.25*sin2p*sint*(1.0 + 7.0*cost*cost - 2.0*sint*sint);
    transf(5,6) = sqrt(6.0)*0.25*cos2p*sint*(1.0 + 7.0*cost*cost - 2.0*sint*sint);
    transf(6,6) = -0.75*cos3p*cost*(1.0 + cost*cost - 2.0*sint*sint);
    transf(7,6) = -0.75*sin3p*cost*(1.0 + cost*cost - 2.0*sint*sint);
    //-phi
    transf(1,7) = 0.0;
    transf(2,7) = -sqrt(15.0)*0.5*sinp*cos2t;
    transf(3,7) = sqrt(15.0)*0.5*cosp*cos2t;
    transf(4,7) = sqrt(6.0)*cos2p*sin2t;
    transf(5,7) = -sqrt(6.0)*sin2p*sin2t;
    transf(6,7) = 1.5*sin3p*cos2t;
    transf(7,7) = -1.5*cos3p*cos2t;
  }
  else {
    transf(4,1) = 1.0;      //f(sigma)
  }
  return transf;
}
matrixE SFtransf_dtdp(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that returns the second-derivatives of rotations for integrals involving s and f orbitals
  //with respect to theta and phi
  matrixE transf(7,7);
  double cos2t = cost*cost - sint*sint;
  double cos2p = cosp*cosp - sinp*sinp;
  double sin2t = 2.0*sint*cost;
  double sin2p = 2.0*sinp*cosp;
  double cos3t = cos2t*cost - sin2t*sint;
  double cos3p = cos2p*cosp - sin2p*sinp;
  double sin3t = cos2t*sint + sin2t*cost;
  double sin3p = cos2p*sinp + sin2p*cosp;
  if (rotate) {
    //order of cartesian orbitals: fz(x2-y2); fx(x2-3y2); fy(3x2-y2); fxyz; fz3; fxz2; fyz2
    //sigma
    transf(1,1) = 0.0;
    transf(2,1) = -sqrt(6.0)*0.25*sinp*cost*(5*cos2t - 1.0 - 5.0*sint*sint);
    transf(3,1) = sqrt(6.0)*0.25*cosp*cost*(5*cos2t - 1.0 - 5.0*sint*sint);
    transf(4,1) = sqrt(15.0)*0.5*cos2p*(cost*sin2t + 2.0*sint*cos2t);
    transf(5,1) = -sqrt(15.0)*0.5*sin2p*(cost*sin2t + 2.0*sint*cos2t);
    transf(6,1) = -sqrt(10.0)*2.25*sin3p*sint*sint*cost;
    transf(7,1) = sqrt(10.0)*2.25*cos3p*sint*sint*cost;
    //pi
    transf(1,1) = 0.0;
    transf(2,2) = 0.25*sinp*sint*(45.0*cost*cost - 11.0);
    transf(3,2) = -0.25*cosp*sint*(45.0*cost*cost - 11.0);
    transf(4,2) = sqrt(10.0)*0.5*cos2p*cost*(3.0*cos2t - 1.0 - 3.0*sint*sint);
    transf(5,2) = -sqrt(10.0)*0.5*sin2p*cost*(3.0*cos2t - 1.0 - 3.0*sint*sint);
    transf(6,2) = -sqrt(15.0)*0.75*sin3p*sint*(cost*cost + cos2t);
    transf(7,2) = sqrt(15.0)*0.75*cos3p*sint*(cost*cost + cos2t);
    //-pi
    transf(1,3) = 0.0;
    transf(2,3) = 2.5*cosp*sint*cost;
    transf(3,3) = 2.5*sinp*sint*cost;
    transf(4,3) = -sqrt(10.0)*sin2p*cos2t;
    transf(5,3) = -sqrt(10.0)*cos2p*cos2t;
    transf(6,3) = -sqrt(15.0)*0.75*cos3p*sin2t;
    transf(7,3) = -sqrt(15.0)*0.75*sin2t*sinp*(4.0*cosp*cosp - 1.0);
    //delta
    transf(1,4) = 0.0;
    transf(2,4) = -sqrt(10.0)*0.25*sinp*cost*(1.0 + 3.0*sint*sint - 3.0*cos2t);
    transf(3,4) = sqrt(10.0)*0.25*cosp*cost*(1.0 + 3.0*sint*sint - 3.0*cos2t);
    transf(4,4) = -cos2p*sint*(9.0*cost*cost - 1.0);
    transf(5,4) = sin2p*sint*(9.0*cost*cost - 1.0);
    transf(6,4) = -sqrt(6.0)*0.75*sin3p*cost*(1.0 + cos2t - sint*sint);
    transf(7,4) = sqrt(6.0)*0.75*cos3p*cost*(1.0 + cos2t - sint*sint);
    //-delta
    transf(1,5) = 0.0;
    transf(2,5) = sqrt(10.0)*0.5*cosp*cos2t;
    transf(3,5) = sqrt(10.0)*0.5*sinp*cos2t;
    transf(4,5) = 4.0*sin2p*sin2t;
    transf(5,5) = 4.0*cos2p*sin2t;
    transf(6,5) = -sqrt(6.0)*1.5*cos3p*cos2t;
    transf(7,5) = -sqrt(6.0)*1.5*sin3p*cos2t;
    //phi
    transf(1,6) = 0.0;
    transf(2,6) = -sqrt(15.0)*0.25*sinp*sint*(cost*cost + cos2t);
    transf(3,6) = sqrt(15.0)*0.25*cosp*sint*(cost*cost + cos2t);
    transf(4,6) = -sqrt(6.0)*0.5*cos2p*cost*(1.0 + cos2t - sint*sint);
    transf(5,6) = sqrt(6.0)*0.5*sin2p*cost*(1.0 + cos2t - sint*sint);
    transf(6,6) = 2.25*sin3p*sint*(1.0 + cost*cost);
    transf(7,6) = -2.25*cos3p*sint*(1.0 + cost*cost);
    //-phi
    transf(1,7) = 0.0;
    transf(2,7) = -sqrt(15.0)*0.5*cosp*sint*cost;
    transf(3,7) = -sqrt(15.0)*0.5*sinp*sint*cost;
    transf(4,7) = sqrt(6.0)*sin2p*cos2t;
    transf(5,7) = sqrt(6.0)*cos2p*cos2t;
    transf(6,7) = 2.25*cos3p*sin2t;
    transf(7,7) = 2.25*sin3p*sin2t;
  }
  else {
    transf(4,1) = 1.0;      //f(sigma)
  }
  return transf;
}
matrixE SFtransf_dp2(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that returns the second-derivatives of rotations for integrals involving s and f orbitals
  //with respect to phi
  matrixE transf(7,7);
  double cos2t = cost*cost - sint*sint;
  double cos2p = cosp*cosp - sinp*sinp;
  double sin2t = 2.0*sint*cost;
  double sin2p = 2.0*sinp*cosp;
  double cos3t = cos2t*cost - sin2t*sint;
  double cos3p = cos2p*cosp - sin2p*sinp;
  double sin3t = cos2t*sint + sin2t*cost;
  double sin3p = cos2p*sinp + sin2p*cosp;
  if (rotate) {
    //order of cartesian orbitals: fz(x2-y2); fx(x2-3y2); fy(3x2-y2); fxyz; fz3; fxz2; fyz2
    //sigma
    transf(1,1) = 0.0;
    transf(2,1) = -sqrt(6.0)*0.25*cosp*sint*(5*cost*cost - 1.0);
    transf(3,1) = -sqrt(6.0)*0.25*sinp*sint*(5*cost*cost - 1.0);
    transf(4,1) = -sqrt(15.0)*sint*sin2p*sin2t;
    transf(5,1) = -sqrt(15.0)*sint*cos2p*sin2t;
    transf(6,1) = -sqrt(10.0)*2.25*cos3p*sint*sint*sint;
    transf(7,1) = -sqrt(10.0)*2.25*sin3p*sint*sint*sint;
    //pi
    transf(1,1) = 0.0;
    transf(2,2) = -0.25*cosp*cost*(15.0*cost*cost - 11.0);
    transf(3,2) = -0.25*sinp*cost*(15.0*cost*cost - 11.0);
    transf(4,2) = -sqrt(10.0)*sin2p*sint*(3.0*cost*cost - 1.0);
    transf(5,2) = -sqrt(10.0)*cos2p*sint*(3.0*cost*cost - 1.0);
    transf(6,2) = -sqrt(15.0)*2.25*cos3p*cost*sint*sint;
    transf(7,2) = -sqrt(15.0)*2.25*sin3p*cost*sint*sint;
    //-pi
    transf(1,3) = 0.0;
    transf(2,3) = 0.25*sinp*(5.0*cost*cost - 1.0);
    transf(3,3) = -0.25*cosp*(5.0*cost*cost - 1.0);
    transf(4,3) = -sqrt(10.0)*sin2t*cos2p;
    transf(5,3) = sqrt(10.0)*sin2t*sin2p;
    transf(6,3) = sqrt(15.0)*2.25*sin3p*sint*sint;
    transf(7,3) = -sqrt(15.0)*0.75*sint*sint*cosp*(4.0*cos2p - 4.0*sinp*sinp - 1.0);
    //delta
    transf(1,4) = 0.0;
    transf(2,4) = sqrt(10.0)*0.25*cosp*sint*(3.0*cost*cost - 1.0);
    transf(3,4) = sqrt(10.0)*0.25*sinp*sint*(3.0*cost*cost - 1.0);
    transf(4,4) = -2.0*sin2p*cost*(3.0*cost*cost - 1.0);
    transf(5,4) = -2.0*cos2p*cost*(3.0*cost*cost - 1.0);
    transf(6,4) = -sqrt(6.0)*2.25*cos3p*sint*(1.0 + cost*cost);
    transf(7,4) = -sqrt(6.0)*2.25*sin3p*sint*(1.0 + cost*cost);
    //-delta
    transf(1,5) = 0.0;
    transf(2,5) = -sqrt(10.0)*0.25*sinp*sin2t;
    transf(3,5) = sqrt(10.0)*0.25*cosp*sin2t;
    transf(4,5) = -4.0*cos2p*cos2t;
    transf(5,5) = 4.0*sin2p*cos2t;
    transf(6,5) = sqrt(6.0)*2.25*sin2t*sin3p;
    transf(7,5) = -sqrt(6.0)*2.25*sin2t*cos3p;
    //phi
    transf(1,6) = 0.0;
    transf(2,6) = -sqrt(15.0)*0.25*cosp*cost*sint*sint;
    transf(3,6) = -sqrt(15.0)*0.25*sinp*cost*sint*sint;
    transf(4,6) = sqrt(6.0)*sin2p*sint*(1.0 + cost*cost);
    transf(5,6) = sqrt(6.0)*cos2p*sint*(1.0 + cost*cost);
    transf(6,6) = -2.25*cos3p*cost*(3.0 + cost*cost);
    transf(7,6) = -2.25*sin3p*cost*(3.0 + cost*cost);
    //-phi
    transf(1,7) = 0.0;
    transf(2,7) = sqrt(15.0)*0.25*sinp*sint*sint;
    transf(3,7) = -sqrt(15.0)*0.25*cosp*sint*sint;
    transf(4,7) = sqrt(6.0)*cos2p*sin2t;
    transf(5,7) = -sqrt(6.0)*sin2p*sin2t;
    transf(6,7) = 2.25*sin3p*(1.0 + 3.0*cost*cost);
    transf(7,7) = -2.25*cos3p*(1.0 + 3.0*cost*cost);
  }
  else {
    transf(4,1) = 1.0;      //f(sigma)
  }
  return transf;
}
matrixE PPtransf(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that calculates the rotation for integrals involving only p orbitals
  //-pi stands for pi bar
  matrixE transf(6,6);
  if (!(rotate)) {
    transf(1,2) = 1.0;      //<p(-pi)|p(-pi)> pi part
    transf(4,1) = 1.0;      //<p(sigma)|p(sigma)>,<p(sigma)|p(sigma)> sigma part
    transf(6,2) = 1.0;      //<p(pi)|p(pi)>,<p(pi)|p(pi)> pi part
  }
  else if (rotate) {
    //order of rows: px|px; px|py; px|pz; py|py; py|pz; pz|pz;
    matrixE transf_sp = SPtransf(ct,st,cp,sp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sp.rows(); ++idx1) {
      for (size_t idx2 = idx1; idx2 < transf_sp.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 4; ++idx3) {
          transf(counter + 1,idx3) = transf_sp(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx3);
          for (size_t idx4 = idx3 + 1; idx4 < 4; ++idx4, ++cnt) {
            transf(counter + 1,3 + cnt) = transf_sp(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx4) + transf_sp(idx1 + 1,idx4)*transf_sp(idx2 + 1,idx3);
          }
        }
      }
    }
  }
  return transf;
}
matrixE PPtransf_dt(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that calculates the derivatives of rotation for integrals involving only p orbitals with respect to theta
  //-pi stands for pi bar
  matrixE transf(6,6);
  if (!(rotate)) {
    transf(1,2) = 1.0;      //<p(-pi)|p(-pi)> pi part
    transf(4,1) = 1.0;      //<p(sigma)|p(sigma)>,<p(sigma)|p(sigma)> sigma part
    transf(6,2) = 1.0;      //<p(pi)|p(pi)>,<p(pi)|p(pi)> pi part
  }
  else if (rotate) {
    //order of rows: px|px; px|py; px|pz; py|py; py|pz; pz|pz;
    matrixE transf_sp = SPtransf(ct,st,cp,sp);
    matrixE transf_sp_dt = SPtransf_dt(ct,st,cp,sp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sp.rows(); ++idx1) {
      for (size_t idx2 = idx1; idx2 < transf_sp.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 4; ++idx3) {
          transf(counter + 1,idx3) = transf_sp_dt(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx3) + transf_sp(idx1 + 1,idx3)*transf_sp_dt(idx2 + 1,idx3);
          for (size_t idx4 = idx3 + 1; idx4 < 4; ++idx4, ++cnt) {
            transf(counter + 1,3 + cnt) = transf_sp_dt(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx4) + transf_sp(idx1 + 1,idx3)*transf_sp_dt(idx2 + 1,idx4);
            transf(counter + 1,3 + cnt) += transf_sp_dt(idx1 + 1,idx4)*transf_sp(idx2 + 1,idx3) + transf_sp(idx1 + 1,idx4)*transf_sp_dt(idx2 + 1,idx3);
          }
        }
      }
    }
  }
  return transf;
}
matrixE PPtransf_dp(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that calculates the derivatives of rotation for integrals involving only p orbitals with respect to phi
  //-pi stands for pi bar
  matrixE transf(6,6);
  if (!(rotate)) {
    transf(1,2) = 1.0;      //<p(-pi)|p(-pi)> pi part
    transf(4,1) = 1.0;      //<p(sigma)|p(sigma)>,<p(sigma)|p(sigma)> sigma part
    transf(6,2) = 1.0;      //<p(pi)|p(pi)>,<p(pi)|p(pi)> pi part
  }
  else if (rotate) {
    //order of rows: px|px; px|py; px|pz; py|py; py|pz; pz|pz;
    matrixE transf_sp = SPtransf(ct,st,cp,sp);
    matrixE transf_sp_dp = SPtransf_dp(ct,st,cp,sp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sp.rows(); ++idx1) {
      for (size_t idx2 = idx1; idx2 < transf_sp.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 4; ++idx3) {
          transf(counter + 1,idx3) = transf_sp_dp(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx3) + transf_sp(idx1 + 1,idx3)*transf_sp_dp(idx2 + 1,idx3);
          for (size_t idx4 = idx3 + 1; idx4 < 4; ++idx4, ++cnt) {
            transf(counter + 1,3 + cnt) = transf_sp_dp(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx4) + transf_sp(idx1 + 1,idx3)*transf_sp_dp(idx2 + 1,idx4);
            transf(counter + 1,3 + cnt) += transf_sp_dp(idx1 + 1,idx4)*transf_sp(idx2 + 1,idx3) + transf_sp(idx1 + 1,idx4)*transf_sp_dp(idx2 + 1,idx3);
          }
        }
      }
    }
  }
  return transf;
}
matrixE PPtransf_dt2(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that calculates the second-derivatives of rotation for integrals involving only p orbitals with respect to theta
  //-pi stands for pi bar
  matrixE transf(6,6);
  if (!(rotate)) {
    transf(1,2) = 1.0;      //<p(-pi)|p(-pi)> pi part
    transf(4,1) = 1.0;      //<p(sigma)|p(sigma)>,<p(sigma)|p(sigma)> sigma part
    transf(6,2) = 1.0;      //<p(pi)|p(pi)>,<p(pi)|p(pi)> pi part
  }
  else if (rotate) {
    //order of rows: px|px; px|py; px|pz; py|py; py|pz; pz|pz;
    matrixE transf_sp = SPtransf(ct,st,cp,sp);
    matrixE transf_sp_dt = SPtransf_dt(ct,st,cp,sp);
    matrixE transf_sp_dt2 = SPtransf_dt2(ct,st,cp,sp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sp.rows(); ++idx1) {
      for (size_t idx2 = idx1; idx2 < transf_sp.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 4; ++idx3) {
          transf(counter + 1,idx3) = transf_sp_dt2(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx3) + transf_sp(idx1 + 1,idx3)*transf_sp_dt2(idx2 + 1,idx3);
          transf(counter + 1,idx3) += 2.0*transf_sp_dt(idx1 + 1,idx3)*transf_sp_dt(idx2 + 1,idx3);
          for (size_t idx4 = idx3 + 1; idx4 < 4; ++idx4, ++cnt) {
            transf(counter + 1,3 + cnt) = transf_sp_dt2(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx4) + transf_sp(idx1 + 1,idx3)*transf_sp_dt2(idx2 + 1,idx4);
            transf(counter + 1,3 + cnt) += 2.0*(transf_sp_dt(idx1 + 1,idx3)*transf_sp_dt(idx2 + 1,idx4) + transf_sp_dt(idx1 + 1,idx4)*transf_sp_dt(idx2 + 1,idx3));
            transf(counter + 1,3 + cnt) += transf_sp_dt2(idx1 + 1,idx4)*transf_sp(idx2 + 1,idx3) + transf_sp(idx1 + 1,idx4)*transf_sp_dt2(idx2 + 1,idx3);
          }
        }
      }
    }
  }
  return transf;
}
matrixE PPtransf_dtdp(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that calculates the second-derivatives of rotation for integrals involving only p orbitals with respect to theta and phi
  //-pi stands for pi bar
  matrixE transf(6,6);
  if (!(rotate)) {
    transf(1,2) = 1.0;      //<p(-pi)|p(-pi)> pi part
    transf(4,1) = 1.0;      //<p(sigma)|p(sigma)>,<p(sigma)|p(sigma)> sigma part
    transf(6,2) = 1.0;      //<p(pi)|p(pi)>,<p(pi)|p(pi)> pi part
  }
  else if (rotate) {
    //order of rows: px|px; px|py; px|pz; py|py; py|pz; pz|pz;
    matrixE transf_sp = SPtransf(ct,st,cp,sp);
    matrixE transf_sp_dt = SPtransf_dt(ct,st,cp,sp);
    matrixE transf_sp_dp = SPtransf_dp(ct,st,cp,sp);
    matrixE transf_sp_dtdp = SPtransf_dtdp(ct,st,cp,sp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sp.rows(); ++idx1) {
      for (size_t idx2 = idx1; idx2 < transf_sp.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 4; ++idx3) {
          transf(counter + 1,idx3) = transf_sp_dtdp(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx3) + transf_sp(idx1 + 1,idx3)*transf_sp_dtdp(idx2 + 1,idx3);
          transf(counter + 1,idx3) += transf_sp_dt(idx1 + 1,idx3)*transf_sp_dp(idx2 + 1,idx3) + transf_sp_dp(idx1 + 1,idx3)*transf_sp_dt(idx2 + 1,idx3);
          for (size_t idx4 = idx3 + 1; idx4 < 4; ++idx4, ++cnt) {
            transf(counter + 1,3 + cnt) = transf_sp_dtdp(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx4) + transf_sp(idx1 + 1,idx3)*transf_sp_dtdp(idx2 + 1,idx4);
            transf(counter + 1,3 + cnt) += transf_sp_dt(idx1 + 1,idx3)*transf_sp_dp(idx2 + 1,idx4) + transf_sp_dp(idx1 + 1,idx3)*transf_sp_dt(idx2 + 1,idx4);
            transf(counter + 1,3 + cnt) += transf_sp_dtdp(idx1 + 1,idx4)*transf_sp(idx2 + 1,idx3) + transf_sp(idx1 + 1,idx4)*transf_sp_dtdp(idx2 + 1,idx3);
            transf(counter + 1,3 + cnt) += transf_sp_dt(idx1 + 1,idx4)*transf_sp_dp(idx2 + 1,idx3) + transf_sp_dp(idx1 + 1,idx4)*transf_sp_dt(idx2 + 1,idx3);
          }
        }
      }
    }
  }
  return transf;
}
matrixE PPtransf_dp2(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that calculates the second-derivatives of rotation for integrals involving only p orbitals with respect to phi
  //-pi stands for pi bar
  matrixE transf(6,6);
  if (!(rotate)) {
    transf(1,2) = 1.0;      //<p(-pi)|p(-pi)> pi part
    transf(4,1) = 1.0;      //<p(sigma)|p(sigma)>,<p(sigma)|p(sigma)> sigma part
    transf(6,2) = 1.0;      //<p(pi)|p(pi)>,<p(pi)|p(pi)> pi part
  }
  else if (rotate) {
    //order of rows: px|px; px|py; px|pz; py|py; py|pz; pz|pz;
    matrixE transf_sp = SPtransf(ct,st,cp,sp);
    matrixE transf_sp_dp = SPtransf_dp(ct,st,cp,sp);
    matrixE transf_sp_dp2 = SPtransf_dp2(ct,st,cp,sp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sp.rows(); ++idx1) {
      for (size_t idx2 = idx1; idx2 < transf_sp.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 4; ++idx3) {
          transf(counter + 1,idx3) = transf_sp_dp2(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx3) + transf_sp(idx1 + 1,idx3)*transf_sp_dp2(idx2 + 1,idx3);
          transf(counter + 1,idx3) += transf_sp_dp(idx1 + 1,idx3)*transf_sp_dp(idx2 + 1,idx3) + transf_sp_dp(idx1 + 1,idx3)*transf_sp_dp(idx2 + 1,idx3);
          for (size_t idx4 = idx3 + 1; idx4 < 4; ++idx4, ++cnt) {
            transf(counter + 1,3 + cnt) = transf_sp_dp2(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx4) + transf_sp(idx1 + 1,idx3)*transf_sp_dp2(idx2 + 1,idx4);
            transf(counter + 1,3 + cnt) += 2.0*(transf_sp_dp(idx1 + 1,idx3)*transf_sp_dp(idx2 + 1,idx4) + transf_sp_dp(idx1 + 1,idx4)*transf_sp_dp(idx2 + 1,idx3));
            transf(counter + 1,3 + cnt) += transf_sp_dp2(idx1 + 1,idx4)*transf_sp(idx2 + 1,idx3) + transf_sp(idx1 + 1,idx4)*transf_sp_dp2(idx2 + 1,idx3);
          }
        }
      }
    }
  }
  return transf;
}
matrixE PDtransf(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that calculates the rotation for integrals involving p and d orbitals
  matrixE transf(15,15);
  if (!(rotate)) {
    transf(4,2) = 1.0;      //<d(-pi)|p(-pi)>,<p(-pi)|d(-pi)> pi part
    transf(8,1) = 1.0;      //<d(sigma)|p(sigma)>,<p(sigma)|d(sigma)> sigma part
    transf(12,2) = 1.0;     //<d(pi)|p(pi)>,<p(pi)|d(pi)> pi part
  }
  else if (rotate) {
    //order of rows: px|dxz; py|dxz; pz|dxz; px|dyz; py|dyz; pz|dyz; px|dz2; py|dz2; pz|dz2; px|dxy; py|dxy; pz|dxy; px|dx2-y2; py|dx2-y2; pz|dx2-y2; 
    //order of columns: sigma|sigma; pi|pi; -pi|-pi; sigma|pi; sigma|-pi; sigma|delta; sigma|-delta; pi|sigma; pi|-pi; pi|delta; pi|-delta; -pi|sigma; -pi|pi; -pi|delta; -pi|-delta; 
    matrixE transf_sd = SDtransf(ct,st,cp,sp);
    matrixE transf_sp = SPtransf(ct,st,cp,sp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sd.rows(); ++idx1) {
      for (size_t idx2 = 0; idx2 < transf_sp.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 4; ++idx3) {
          transf(counter + 1,idx3) = transf_sd(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx3);
          for (size_t idx4 = 1; idx4 < 6; ++idx4) {
            if (idx4 == idx3) continue;
            transf(counter + 1,3 + cnt) = transf_sd(idx1 + 1,idx4)*transf_sp(idx2 + 1,idx3);
            ++cnt;
          }
        }
      }
    }
  }
  return transf;
}
matrixE PDtransf_dt(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that calculates the derivatives rotation for integrals involving p and d orbitals with respect to theta
  matrixE transf(15,15);
  if (!(rotate)) {
    transf(4,2) = 1.0;      //<d(-pi)|p(-pi)>,<p(-pi)|d(-pi)> pi part
    transf(8,1) = 1.0;      //<d(sigma)|p(sigma)>,<p(sigma)|d(sigma)> sigma part
    transf(12,2) = 1.0;     //<d(pi)|p(pi)>,<p(pi)|d(pi)> pi part
  }
  else if (rotate) {
    //order of rows: px|dxz; py|dxz; pz|dxz; px|dyz; py|dyz; pz|dyz; px|dz2; py|dz2; pz|dz2; px|dxy; py|dxy; pz|dxy; px|dx2-y2; py|dx2-y2; pz|dx2-y2; 
    //order of columns: sigma|sigma; pi|pi; -pi|-pi; sigma|pi; sigma|-pi; sigma|delta; sigma|-delta; pi|sigma; pi|-pi; pi|delta; pi|-delta; -pi|sigma; -pi|pi; -pi|delta; -pi|-delta; 
    matrixE transf_sd = SDtransf(ct,st,cp,sp);
    matrixE transf_sp = SPtransf(ct,st,cp,sp);
    matrixE transf_sd_dt = SDtransf_dt(ct,st,cp,sp);
    matrixE transf_sp_dt = SPtransf_dt(ct,st,cp,sp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sd.rows(); ++idx1) {
      for (size_t idx2 = 0; idx2 < transf_sp.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 4; ++idx3) {
          transf(counter + 1,idx3) = transf_sd_dt(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx3) + transf_sd(idx1 + 1,idx3)*transf_sp_dt(idx2 + 1,idx3);
          for (size_t idx4 = 1; idx4 < 6; ++idx4) {
            if (idx4 == idx3) continue;
            transf(counter + 1,3 + cnt) = transf_sd_dt(idx1 + 1,idx4)*transf_sp(idx2 + 1,idx3) + transf_sd(idx1 + 1,idx4)*transf_sp_dt(idx2 + 1,idx3);
            ++cnt;
          }
        }
      }
    }
  }
  return transf;
}
matrixE PDtransf_dp(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that calculates the derivatives rotation for integrals involving p and d orbitals with respect to phi
  matrixE transf(15,15);
  if (!(rotate)) {
    transf(4,2) = 1.0;      //<d(-pi)|p(-pi)>,<p(-pi)|d(-pi)> pi part
    transf(8,1) = 1.0;      //<d(sigma)|p(sigma)>,<p(sigma)|d(sigma)> sigma part
    transf(12,2) = 1.0;     //<d(pi)|p(pi)>,<p(pi)|d(pi)> pi part
  }
  else if (rotate) {
    //order of rows: px|dxz; py|dxz; pz|dxz; px|dyz; py|dyz; pz|dyz; px|dz2; py|dz2; pz|dz2; px|dxy; py|dxy; pz|dxy; px|dx2-y2; py|dx2-y2; pz|dx2-y2; 
    //order of columns: sigma|sigma; pi|pi; -pi|-pi; sigma|pi; sigma|-pi; sigma|delta; sigma|-delta; pi|sigma; pi|-pi; pi|delta; pi|-delta; -pi|sigma; -pi|pi; -pi|delta; -pi|-delta; 
    matrixE transf_sd = SDtransf(ct,st,cp,sp);
    matrixE transf_sp = SPtransf(ct,st,cp,sp);
    matrixE transf_sd_dp = SDtransf_dp(ct,st,cp,sp);
    matrixE transf_sp_dp = SPtransf_dp(ct,st,cp,sp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sd.rows(); ++idx1) {
      for (size_t idx2 = 0; idx2 < transf_sp.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 4; ++idx3) {
          transf(counter + 1,idx3) = transf_sd_dp(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx3) + transf_sd(idx1 + 1,idx3)*transf_sp_dp(idx2 + 1,idx3);
          for (size_t idx4 = 1; idx4 < 6; ++idx4) {
            if (idx4 == idx3) continue;
            transf(counter + 1,3 + cnt) = transf_sd_dp(idx1 + 1,idx4)*transf_sp(idx2 + 1,idx3) + transf_sd(idx1 + 1,idx4)*transf_sp_dp(idx2 + 1,idx3);
            ++cnt;
          }
        }
      }
    }
  }
  return transf;
}
matrixE PDtransf_dt2(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that calculates the second-derivatives rotation for integrals involving p and d orbitals with respect to theta
  matrixE transf(15,15);
  if (!(rotate)) {
    transf(4,2) = 1.0;      //<d(-pi)|p(-pi)>,<p(-pi)|d(-pi)> pi part
    transf(8,1) = 1.0;      //<d(sigma)|p(sigma)>,<p(sigma)|d(sigma)> sigma part
    transf(12,2) = 1.0;     //<d(pi)|p(pi)>,<p(pi)|d(pi)> pi part
  }
  else if (rotate) {
    //order of rows: px|dxz; py|dxz; pz|dxz; px|dyz; py|dyz; pz|dyz; px|dz2; py|dz2; pz|dz2; px|dxy; py|dxy; pz|dxy; px|dx2-y2; py|dx2-y2; pz|dx2-y2; 
    //order of columns: sigma|sigma; pi|pi; -pi|-pi; sigma|pi; sigma|-pi; sigma|delta; sigma|-delta; pi|sigma; pi|-pi; pi|delta; pi|-delta; -pi|sigma; -pi|pi; -pi|delta; -pi|-delta; 
    matrixE transf_sd = SDtransf(ct,st,cp,sp);
    matrixE transf_sp = SPtransf(ct,st,cp,sp);
    matrixE transf_sd_dt = SDtransf_dt(ct,st,cp,sp);
    matrixE transf_sp_dt = SPtransf_dt(ct,st,cp,sp);
    matrixE transf_sd_dt2 = SDtransf_dt2(ct,st,cp,sp);
    matrixE transf_sp_dt2 = SPtransf_dt2(ct,st,cp,sp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sd.rows(); ++idx1) {
      for (size_t idx2 = 0; idx2 < transf_sp.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 4; ++idx3) {
          transf(counter + 1,idx3) = transf_sd_dt2(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx3) + transf_sd(idx1 + 1,idx3)*transf_sp_dt2(idx2 + 1,idx3);
          transf(counter + 1,idx3) += 2.0*transf_sd_dt(idx1 + 1,idx3)*transf_sp_dt(idx2 + 1,idx3);
          for (size_t idx4 = 1; idx4 < 6; ++idx4) {
            if (idx4 == idx3) continue;
            transf(counter + 1,3 + cnt) = transf_sd_dt2(idx1 + 1,idx4)*transf_sp(idx2 + 1,idx3) + transf_sd(idx1 + 1,idx4)*transf_sp_dt2(idx2 + 1,idx3);
            transf(counter + 1,3 + cnt) += 2.0*transf_sd_dt(idx1 + 1,idx4)*transf_sp_dt(idx2 + 1,idx3);
            ++cnt;
          }
        }
      }
    }
  }
  return transf;
}
matrixE PDtransf_dtdp(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that calculates the second-derivatives rotation for integrals involving p and d orbitals with respect to theta and phi
  matrixE transf(15,15);
  if (!(rotate)) {
    transf(4,2) = 1.0;      //<d(-pi)|p(-pi)>,<p(-pi)|d(-pi)> pi part
    transf(8,1) = 1.0;      //<d(sigma)|p(sigma)>,<p(sigma)|d(sigma)> sigma part
    transf(12,2) = 1.0;     //<d(pi)|p(pi)>,<p(pi)|d(pi)> pi part
  }
  else if (rotate) {
    //order of rows: px|dxz; py|dxz; pz|dxz; px|dyz; py|dyz; pz|dyz; px|dz2; py|dz2; pz|dz2; px|dxy; py|dxy; pz|dxy; px|dx2-y2; py|dx2-y2; pz|dx2-y2; 
    //order of columns: sigma|sigma; pi|pi; -pi|-pi; sigma|pi; sigma|-pi; sigma|delta; sigma|-delta; pi|sigma; pi|-pi; pi|delta; pi|-delta; -pi|sigma; -pi|pi; -pi|delta; -pi|-delta; 
    matrixE transf_sd = SDtransf(ct,st,cp,sp);
    matrixE transf_sp = SPtransf(ct,st,cp,sp);
    matrixE transf_sd_dt = SDtransf_dt(ct,st,cp,sp);
    matrixE transf_sp_dt = SPtransf_dt(ct,st,cp,sp);
    matrixE transf_sd_dp = SDtransf_dp(ct,st,cp,sp);
    matrixE transf_sp_dp = SPtransf_dp(ct,st,cp,sp);
    matrixE transf_sd_dtdp = SDtransf_dtdp(ct,st,cp,sp);
    matrixE transf_sp_dtdp = SPtransf_dtdp(ct,st,cp,sp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sd.rows(); ++idx1) {
      for (size_t idx2 = 0; idx2 < transf_sp.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 4; ++idx3) {
          transf(counter + 1,idx3) = transf_sd_dtdp(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx3) + transf_sd(idx1 + 1,idx3)*transf_sp_dtdp(idx2 + 1,idx3);
          transf(counter + 1,idx3) += transf_sd_dt(idx1 + 1,idx3)*transf_sp_dp(idx2 + 1,idx3) + transf_sd_dp(idx1 + 1,idx3)*transf_sp_dt(idx2 + 1,idx3);
          for (size_t idx4 = 1; idx4 < 6; ++idx4) {
            if (idx4 == idx3) continue;
            transf(counter + 1,3 + cnt) = transf_sd_dtdp(idx1 + 1,idx4)*transf_sp(idx2 + 1,idx3) + transf_sd(idx1 + 1,idx4)*transf_sp_dtdp(idx2 + 1,idx3);
            transf(counter + 1,3 + cnt) += transf_sd_dt(idx1 + 1,idx4)*transf_sp_dp(idx2 + 1,idx3) + transf_sd_dp(idx1 + 1,idx4)*transf_sp_dt(idx2 + 1,idx3);
            ++cnt;
          }
        }
      }
    }
  }
  return transf;
}
matrixE PDtransf_dp2(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that calculates the second-derivatives rotation for integrals involving p and d orbitals with respect to phi
  matrixE transf(15,15);
  if (!(rotate)) {
    transf(4,2) = 1.0;      //<d(-pi)|p(-pi)>,<p(-pi)|d(-pi)> pi part
    transf(8,1) = 1.0;      //<d(sigma)|p(sigma)>,<p(sigma)|d(sigma)> sigma part
    transf(12,2) = 1.0;     //<d(pi)|p(pi)>,<p(pi)|d(pi)> pi part
  }
  else if (rotate) {
    //order of rows: px|dxz; py|dxz; pz|dxz; px|dyz; py|dyz; pz|dyz; px|dz2; py|dz2; pz|dz2; px|dxy; py|dxy; pz|dxy; px|dx2-y2; py|dx2-y2; pz|dx2-y2; 
    //order of columns: sigma|sigma; pi|pi; -pi|-pi; sigma|pi; sigma|-pi; sigma|delta; sigma|-delta; pi|sigma; pi|-pi; pi|delta; pi|-delta; -pi|sigma; -pi|pi; -pi|delta; -pi|-delta; 
    matrixE transf_sd = SDtransf(ct,st,cp,sp);
    matrixE transf_sp = SPtransf(ct,st,cp,sp);
    matrixE transf_sd_dp = SDtransf_dp(ct,st,cp,sp);
    matrixE transf_sp_dp = SPtransf_dp(ct,st,cp,sp);
    matrixE transf_sd_dp2 = SDtransf_dp2(ct,st,cp,sp);
    matrixE transf_sp_dp2 = SPtransf_dp2(ct,st,cp,sp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sd.rows(); ++idx1) {
      for (size_t idx2 = 0; idx2 < transf_sp.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 4; ++idx3) {
          transf(counter + 1,idx3) = transf_sd_dp2(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx3) + transf_sd(idx1 + 1,idx3)*transf_sp_dp2(idx2 + 1,idx3);
          transf(counter + 1,idx3) += 2.0*transf_sd_dp(idx1 + 1,idx3)*transf_sp_dp(idx2 + 1,idx3);
          for (size_t idx4 = 1; idx4 < 6; ++idx4) {
            if (idx4 == idx3) continue;
            transf(counter + 1,3 + cnt) = transf_sd_dp2(idx1 + 1,idx4)*transf_sp(idx2 + 1,idx3) + transf_sd(idx1 + 1,idx4)*transf_sp_dp2(idx2 + 1,idx3);
            transf(counter + 1,3 + cnt) += 2.0*transf_sd_dp(idx1 + 1,idx4)*transf_sp_dp(idx2 + 1,idx3);
            ++cnt;
          }
        }
      }
    }
  }
  return transf;
}
matrixE DDtransf(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that calculates the rotation for integrals involving only d orbitals
  matrixE transf(15,15);
  if (!(rotate)) {
    transf(1,4) = 1.0;       //<d(-delta)|d(-delta)> delta part
    transf(6,2) = 1.0;       //<d(-pi)|d(-pi)> pi part part
    transf(10,1) = 1.0;      //<d(sigma)|d(sigma)> sigma part
    transf(13,2) = 1.0;      //<d(pi)|d(pi)> pi part
    transf(15,4) = 1.0;      //<ddelta)|ddelta)> delta part
  }
  else if (rotate) {
    //order of rows: dxz|dxz; dxz|dyz; dxz|dz2; dxz|dxy; dxz|dx2-y2; dyz|dyz; dyz|dz2; dyz|dxy; dyz|dx2-y2; dz2|dz2; dz2|dxy; dz2|dx2-y2; dxy|dxy; dxy|dx2-y2; dx2-y2|dx2-y2;
    //order of columns: sigma|sigma; pi|pi; -pi|-pi; delta|delta; -delta|-delta; sigma|pi; sigma|-pi; sigma|delta; sigma|-delta; pi|-pi; pi|delta; pi|-delta; -pi|delta; -pi|-delta; delta|-delta; 
    matrixE transf_sd = SDtransf(ct,st,cp,sp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sd.rows(); ++idx1) {
      for (size_t idx2 = idx1; idx2 < transf_sd.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 6; ++idx3) {
          transf(counter + 1,idx3) = transf_sd(idx1 + 1,idx3)*transf_sd(idx2 + 1,idx3);
          for (size_t idx4 = idx3 + 1; idx4 < 6; ++idx4, ++cnt) {
            transf(counter + 1,5 + cnt) = transf_sd(idx1 + 1,idx3)*transf_sd(idx2 + 1,idx4) + transf_sd(idx1 + 1,idx4)*transf_sd(idx2 + 1,idx3);
          }
        }
      }
    }
  }
  return transf;
}
matrixE DDtransf_dt(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that calculates the derivatives of rotation for integrals involving only d orbitals with respect to theta
  matrixE transf(15,15);
  if (!(rotate)) {
    transf(1,4) = 1.0;       //<d(-delta)|d(-delta)> delta part
    transf(6,2) = 1.0;       //<d(-pi)|d(-pi)> pi part part
    transf(10,1) = 1.0;      //<d(sigma)|d(sigma)> sigma part
    transf(13,2) = 1.0;      //<d(pi)|d(pi)> pi part
    transf(15,4) = 1.0;      //<ddelta)|ddelta)> delta part
  }
  else if (rotate) {
    //order of rows: dxz|dxz; dxz|dyz; dxz|dz2; dxz|dxy; dxz|dx2-y2; dyz|dyz; dyz|dz2; dyz|dxy; dyz|dx2-y2; dz2|dz2; dz2|dxy; dz2|dx2-y2; dxy|dxy; dxy|dx2-y2; dx2-y2|dx2-y2;
    //order of columns: sigma|sigma; pi|pi; -pi|-pi; delta|delta; -delta|-delta; sigma|pi; sigma|-pi; sigma|delta; sigma|-delta; pi|-pi; pi|delta; pi|-delta; -pi|delta; -pi|-delta; delta|-delta; 
    matrixE transf_sd = SDtransf(ct,st,cp,sp);
    matrixE transf_sd_dt = SDtransf_dt(ct,st,cp,sp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sd.rows(); ++idx1) {
      for (size_t idx2 = idx1; idx2 < transf_sd.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 6; ++idx3) {
          transf(counter + 1,idx3) = transf_sd_dt(idx1 + 1,idx3)*transf_sd(idx2 + 1,idx3) + transf_sd(idx1 + 1,idx3)*transf_sd_dt(idx2 + 1,idx3);
          for (size_t idx4 = idx3 + 1; idx4 < 6; ++idx4, ++cnt) {
            transf(counter + 1,5 + cnt) = transf_sd_dt(idx1 + 1,idx3)*transf_sd(idx2 + 1,idx4) + transf_sd(idx1 + 1,idx3)*transf_sd_dt(idx2 + 1,idx4);
            transf(counter + 1,5 + cnt) += transf_sd_dt(idx1 + 1,idx4)*transf_sd(idx2 + 1,idx3) + transf_sd(idx1 + 1,idx4)*transf_sd_dt(idx2 + 1,idx3);
          }
        }
      }
    }
  }
  return transf;
}
matrixE DDtransf_dp(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that calculates the derivatives of rotation for integrals involving only d orbitals with respect to phi
  matrixE transf(15,15);
  if (!(rotate)) {
    transf(1,4) = 1.0;       //<d(-delta)|d(-delta)> delta part
    transf(6,2) = 1.0;       //<d(-pi)|d(-pi)> pi part part
    transf(10,1) = 1.0;      //<d(sigma)|d(sigma)> sigma part
    transf(13,2) = 1.0;      //<d(pi)|d(pi)> pi part
    transf(15,4) = 1.0;      //<ddelta)|ddelta)> delta part
  }
  else if (rotate) {
    //order of rows: dxz|dxz; dxz|dyz; dxz|dz2; dxz|dxy; dxz|dx2-y2; dyz|dyz; dyz|dz2; dyz|dxy; dyz|dx2-y2; dz2|dz2; dz2|dxy; dz2|dx2-y2; dxy|dxy; dxy|dx2-y2; dx2-y2|dx2-y2;
    //order of columns: sigma|sigma; pi|pi; -pi|-pi; delta|delta; -delta|-delta; sigma|pi; sigma|-pi; sigma|delta; sigma|-delta; pi|-pi; pi|delta; pi|-delta; -pi|delta; -pi|-delta; delta|-delta; 
    matrixE transf_sd = SDtransf(ct,st,cp,sp);
    matrixE transf_sd_dp = SDtransf_dp(ct,st,cp,sp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sd.rows(); ++idx1) {
      for (size_t idx2 = idx1; idx2 < transf_sd.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 6; ++idx3) {
          transf(counter + 1,idx3) = transf_sd_dp(idx1 + 1,idx3)*transf_sd(idx2 + 1,idx3) + transf_sd(idx1 + 1,idx3)*transf_sd_dp(idx2 + 1,idx3);
          for (size_t idx4 = idx3 + 1; idx4 < 6; ++idx4, ++cnt) {
            transf(counter + 1,5 + cnt) = transf_sd_dp(idx1 + 1,idx3)*transf_sd(idx2 + 1,idx4) + transf_sd(idx1 + 1,idx3)*transf_sd_dp(idx2 + 1,idx4);
            transf(counter + 1,5 + cnt) += transf_sd_dp(idx1 + 1,idx4)*transf_sd(idx2 + 1,idx3) + transf_sd(idx1 + 1,idx4)*transf_sd_dp(idx2 + 1,idx3);
          }
        }
      }
    }
  }
  return transf;
}
matrixE DDtransf_dt2(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that calculates the second-derivatives of rotation for integrals involving only d orbitals with respect to theta
  matrixE transf(15,15);
  if (!(rotate)) {
    transf(1,4) = 1.0;       //<d(-delta)|d(-delta)> delta part
    transf(6,2) = 1.0;       //<d(-pi)|d(-pi)> pi part part
    transf(10,1) = 1.0;      //<d(sigma)|d(sigma)> sigma part
    transf(13,2) = 1.0;      //<d(pi)|d(pi)> pi part
    transf(15,4) = 1.0;      //<ddelta)|ddelta)> delta part
  }
  else if (rotate) {
    //order of rows: dxz|dxz; dxz|dyz; dxz|dz2; dxz|dxy; dxz|dx2-y2; dyz|dyz; dyz|dz2; dyz|dxy; dyz|dx2-y2; dz2|dz2; dz2|dxy; dz2|dx2-y2; dxy|dxy; dxy|dx2-y2; dx2-y2|dx2-y2;
    //order of columns: sigma|sigma; pi|pi; -pi|-pi; delta|delta; -delta|-delta; sigma|pi; sigma|-pi; sigma|delta; sigma|-delta; pi|-pi; pi|delta; pi|-delta; -pi|delta; -pi|-delta; delta|-delta; 
    matrixE transf_sd = SDtransf(ct,st,cp,sp);
    matrixE transf_sd_dt = SDtransf_dt(ct,st,cp,sp);
    matrixE transf_sd_dt2 = SDtransf_dt2(ct,st,cp,sp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sd.rows(); ++idx1) {
      for (size_t idx2 = idx1; idx2 < transf_sd.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 6; ++idx3) {
          transf(counter + 1,idx3) = transf_sd_dt2(idx1 + 1,idx3)*transf_sd(idx2 + 1,idx3) + transf_sd(idx1 + 1,idx3)*transf_sd_dt2(idx2 + 1,idx3);
          transf(counter + 1,idx3) += 2.0*transf_sd_dt(idx1 + 1,idx3)*transf_sd_dt(idx2 + 1,idx3);
          for (size_t idx4 = idx3 + 1; idx4 < 6; ++idx4, ++cnt) {
            transf(counter + 1,5 + cnt) = transf_sd_dt2(idx1 + 1,idx3)*transf_sd(idx2 + 1,idx4) + transf_sd(idx1 + 1,idx3)*transf_sd_dt2(idx2 + 1,idx4);
            transf(counter + 1,5 + cnt) += 2.0*(transf_sd_dt(idx1 + 1,idx3)*transf_sd_dt(idx2 + 1,idx4) + transf_sd_dt(idx1 + 1,idx4)*transf_sd_dt(idx2 + 1,idx3));
            transf(counter + 1,5 + cnt) += transf_sd_dt2(idx1 + 1,idx4)*transf_sd(idx2 + 1,idx3) + transf_sd(idx1 + 1,idx4)*transf_sd_dt2(idx2 + 1,idx3);
          }
        }
      }
    }
  }
  return transf;
}
matrixE DDtransf_dtdp(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that calculates the second-derivatives of rotation for integrals involving only d orbitals with respect to theta and phi
  matrixE transf(15,15);
  if (!(rotate)) {
    transf(1,4) = 1.0;       //<d(-delta)|d(-delta)> delta part
    transf(6,2) = 1.0;       //<d(-pi)|d(-pi)> pi part part
    transf(10,1) = 1.0;      //<d(sigma)|d(sigma)> sigma part
    transf(13,2) = 1.0;      //<d(pi)|d(pi)> pi part
    transf(15,4) = 1.0;      //<ddelta)|ddelta)> delta part
  }
  else if (rotate) {
    //order of rows: dxz|dxz; dxz|dyz; dxz|dz2; dxz|dxy; dxz|dx2-y2; dyz|dyz; dyz|dz2; dyz|dxy; dyz|dx2-y2; dz2|dz2; dz2|dxy; dz2|dx2-y2; dxy|dxy; dxy|dx2-y2; dx2-y2|dx2-y2;
    //order of columns: sigma|sigma; pi|pi; -pi|-pi; delta|delta; -delta|-delta; sigma|pi; sigma|-pi; sigma|delta; sigma|-delta; pi|-pi; pi|delta; pi|-delta; -pi|delta; -pi|-delta; delta|-delta; 
    matrixE transf_sd = SDtransf(ct,st,cp,sp);
    matrixE transf_sd_dt = SDtransf_dt(ct,st,cp,sp);
    matrixE transf_sd_dp = SDtransf_dp(ct,st,cp,sp);
    matrixE transf_sd_dtdp = SDtransf_dtdp(ct,st,cp,sp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sd.rows(); ++idx1) {
      for (size_t idx2 = idx1; idx2 < transf_sd.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 6; ++idx3) {
          transf(counter + 1,idx3) = transf_sd_dtdp(idx1 + 1,idx3)*transf_sd(idx2 + 1,idx3) + transf_sd(idx1 + 1,idx3)*transf_sd_dtdp(idx2 + 1,idx3);
          transf(counter + 1,idx3) += transf_sd_dt(idx1 + 1,idx3)*transf_sd_dp(idx2 + 1,idx3) + transf_sd_dp(idx1 + 1,idx3)*transf_sd_dt(idx2 + 1,idx3);
          for (size_t idx4 = idx3 + 1; idx4 < 6; ++idx4, ++cnt) {
            transf(counter + 1,5 + cnt) = transf_sd_dtdp(idx1 + 1,idx3)*transf_sd(idx2 + 1,idx4) + transf_sd(idx1 + 1,idx3)*transf_sd_dtdp(idx2 + 1,idx4);
            transf(counter + 1,5 + cnt) += transf_sd_dt(idx1 + 1,idx3)*transf_sd_dp(idx2 + 1,idx4) + transf_sd_dp(idx1 + 1,idx3)*transf_sd_dt(idx2 + 1,idx4);
            transf(counter + 1,5 + cnt) += transf_sd_dtdp(idx1 + 1,idx4)*transf_sd(idx2 + 1,idx3) + transf_sd(idx1 + 1,idx4)*transf_sd_dtdp(idx2 + 1,idx3);
            transf(counter + 1,5 + cnt) += transf_sd_dt(idx1 + 1,idx4)*transf_sd_dp(idx2 + 1,idx3) + transf_sd_dp(idx1 + 1,idx4)*transf_sd_dt(idx2 + 1,idx3);
          }
        }
      }
    }
  }
  return transf;
}
matrixE DDtransf_dp2(double ct, double st, double cp, double sp, bool rotate = true) {
  //function that calculates the second-derivatives of rotation for integrals involving only d orbitals with respect to phi
  matrixE transf(15,15);
  if (!(rotate)) {
    transf(1,4) = 1.0;       //<d(-delta)|d(-delta)> delta part
    transf(6,2) = 1.0;       //<d(-pi)|d(-pi)> pi part part
    transf(10,1) = 1.0;      //<d(sigma)|d(sigma)> sigma part
    transf(13,2) = 1.0;      //<d(pi)|d(pi)> pi part
    transf(15,4) = 1.0;      //<ddelta)|ddelta)> delta part
  }
  else if (rotate) {
    //order of rows: dxz|dxz; dxz|dyz; dxz|dz2; dxz|dxy; dxz|dx2-y2; dyz|dyz; dyz|dz2; dyz|dxy; dyz|dx2-y2; dz2|dz2; dz2|dxy; dz2|dx2-y2; dxy|dxy; dxy|dx2-y2; dx2-y2|dx2-y2;
    //order of columns: sigma|sigma; pi|pi; -pi|-pi; delta|delta; -delta|-delta; sigma|pi; sigma|-pi; sigma|delta; sigma|-delta; pi|-pi; pi|delta; pi|-delta; -pi|delta; -pi|-delta; delta|-delta; 
    matrixE transf_sd = SDtransf(ct,st,cp,sp);
    matrixE transf_sd_dp = SDtransf_dp(ct,st,cp,sp);
    matrixE transf_sd_dp2 = SDtransf_dp2(ct,st,cp,sp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sd.rows(); ++idx1) {
      for (size_t idx2 = idx1; idx2 < transf_sd.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 6; ++idx3) {
          transf(counter + 1,idx3) = transf_sd_dp2(idx1 + 1,idx3)*transf_sd(idx2 + 1,idx3) + transf_sd(idx1 + 1,idx3)*transf_sd_dp2(idx2 + 1,idx3);
          transf(counter + 1,idx3) += 2.0*transf_sd_dp(idx1 + 1,idx3)*transf_sd_dp(idx2 + 1,idx3);
          for (size_t idx4 = idx3 + 1; idx4 < 6; ++idx4, ++cnt) {
            transf(counter + 1,5 + cnt) = transf_sd_dp2(idx1 + 1,idx3)*transf_sd(idx2 + 1,idx4) + transf_sd(idx1 + 1,idx3)*transf_sd_dp2(idx2 + 1,idx4);
            transf(counter + 1,5 + cnt) += 2.0*(transf_sd_dp(idx1 + 1,idx3)*transf_sd_dp(idx2 + 1,idx4) + transf_sd_dp(idx1 + 1,idx4)*transf_sd_dp(idx2 + 1,idx3));
            transf(counter + 1,5 + cnt) += transf_sd_dp2(idx1 + 1,idx4)*transf_sd(idx2 + 1,idx3) + transf_sd(idx1 + 1,idx4)*transf_sd_dp2(idx2 + 1,idx3);
          }
        }
      }
    }
  }
  return transf;
}
matrixE PFtransf(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that calculates the rotation for integrals involving p and f orbitals
  matrixE transf(21,21);
  if (!(rotate)) {
    transf(3,2) = 1.0;       //<p(-pi)|f(-pi)>,<f(-pi)|p(-pi)> pi part
    transf(11,1) = 1.0;      //<p(sigma)|f(sigma)>,<f(sigma)|p(sigma)> sigma part
    transf(19,2) = 1.0;      //<p(pi)|f(pi)>,<f(pi)|p(pi)> pi part
  }
  else if (rotate) {
    matrixE transf_sf = SFtransf(cost,sint,cosp,sinp);
    matrixE transf_sp = SPtransf(cost,sint,cosp,sinp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sf.rows(); ++idx1) {
      for (size_t idx2 = 0; idx2 < transf_sp.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 4; ++idx3) {
          transf(counter + 1,idx3) = transf_sf(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx3);
          for (size_t idx4 = 1; idx4 < 8; ++idx4) {
            if (idx4 == idx3) continue;
            transf(counter + 1,3 + cnt) = transf_sf(idx1 + 1,idx4)*transf_sp(idx2 + 1,idx3);
            ++cnt;
          }
        }
      }
    }
  }
  return transf;
}
matrixE PFtransf_dt(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that calculates the derivatives of rotation for integrals involving p and f orbitals with respect to theta
  matrixE transf(21,21);
  if (!(rotate)) {
    transf(3,2) = 1.0;       //<p(-pi)|f(-pi)>,<f(-pi)|p(-pi)> pi part
    transf(11,1) = 1.0;      //<p(sigma)|f(sigma)>,<f(sigma)|p(sigma)> sigma part
    transf(19,2) = 1.0;      //<p(pi)|f(pi)>,<f(pi)|p(pi)> pi part
  }
  else if (rotate) {
    matrixE transf_sf = SFtransf(cost,sint,cosp,sinp);
    matrixE transf_sp = SPtransf(cost,sint,cosp,sinp);
    matrixE transf_sf_dt = SFtransf_dt(cost,sint,cosp,sinp);
    matrixE transf_sp_dt = SPtransf_dt(cost,sint,cosp,sinp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sf.rows(); ++idx1) {
      for (size_t idx2 = 0; idx2 < transf_sp.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 4; ++idx3) {
          transf(counter + 1,idx3) = transf_sf_dt(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx3)*transf_sp_dt(idx2 + 1,idx3);
          for (size_t idx4 = 1; idx4 < 8; ++idx4) {
            if (idx4 == idx3) continue;
            transf(counter + 1,3 + cnt) = transf_sf_dt(idx1 + 1,idx4)*transf_sp(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx4)*transf_sp_dt(idx2 + 1,idx3);
            ++cnt;
          }
        }
      }
    }
  }
  return transf;
}
matrixE PFtransf_dp(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that calculates the derivatives of rotation for integrals involving p and f orbitals with respect to phi
  matrixE transf(21,21);
  if (!(rotate)) {
    transf(3,2) = 1.0;       //<p(-pi)|f(-pi)>,<f(-pi)|p(-pi)> pi part
    transf(11,1) = 1.0;      //<p(sigma)|f(sigma)>,<f(sigma)|p(sigma)> sigma part
    transf(19,2) = 1.0;      //<p(pi)|f(pi)>,<f(pi)|p(pi)> pi part
  }
  else if (rotate) {
    matrixE transf_sf = SFtransf(cost,sint,cosp,sinp);
    matrixE transf_sp = SPtransf(cost,sint,cosp,sinp);
    matrixE transf_sf_dp = SFtransf_dp(cost,sint,cosp,sinp);
    matrixE transf_sp_dp = SPtransf_dp(cost,sint,cosp,sinp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sf.rows(); ++idx1) {
      for (size_t idx2 = 0; idx2 < transf_sp.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 4; ++idx3) {
          transf(counter + 1,idx3) = transf_sf_dp(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx3)*transf_sp_dp(idx2 + 1,idx3);
          for (size_t idx4 = 1; idx4 < 8; ++idx4) {
            if (idx4 == idx3) continue;
            transf(counter + 1,3 + cnt) = transf_sf_dp(idx1 + 1,idx4)*transf_sp(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx4)*transf_sp_dp(idx2 + 1,idx3);
            ++cnt;
          }
        }
      }
    }
  }
  return transf;
}
matrixE PFtransf_dt2(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that calculates the second-derivatives of rotation for integrals involving p and f orbitals with respect to theta
  matrixE transf(21,21);
  if (!(rotate)) {
    transf(3,2) = 1.0;       //<p(-pi)|f(-pi)>,<f(-pi)|p(-pi)> pi part
    transf(11,1) = 1.0;      //<p(sigma)|f(sigma)>,<f(sigma)|p(sigma)> sigma part
    transf(19,2) = 1.0;      //<p(pi)|f(pi)>,<f(pi)|p(pi)> pi part
  }
  else if (rotate) {
    matrixE transf_sf = SFtransf(cost,sint,cosp,sinp);
    matrixE transf_sp = SPtransf(cost,sint,cosp,sinp);
    matrixE transf_sf_dt = SFtransf_dt(cost,sint,cosp,sinp);
    matrixE transf_sp_dt = SPtransf_dt(cost,sint,cosp,sinp);
    matrixE transf_sf_dt2 = SFtransf_dt2(cost,sint,cosp,sinp);
    matrixE transf_sp_dt2 = SPtransf_dt2(cost,sint,cosp,sinp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sf.rows(); ++idx1) {
      for (size_t idx2 = 0; idx2 < transf_sp.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 4; ++idx3) {
          transf(counter + 1,idx3) = transf_sf_dt2(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx3)*transf_sp_dt2(idx2 + 1,idx3);
          transf(counter + 1,idx3) += 2.0*transf_sf_dt(idx1 + 1,idx3)*transf_sp_dt(idx2 + 1,idx3);
          for (size_t idx4 = 1; idx4 < 8; ++idx4) {
            if (idx4 == idx3) continue;
            transf(counter + 1,3 + cnt) = transf_sf_dt2(idx1 + 1,idx4)*transf_sp(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx4)*transf_sp_dt2(idx2 + 1,idx3);
            transf(counter + 1,3 + cnt) += 2.0*transf_sf_dt(idx1 + 1,idx4)*transf_sp_dt(idx2 + 1,idx3);
            ++cnt;
          }
        }
      }
    }
  }
  return transf;
}
matrixE PFtransf_dtdp(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that calculates the derivatives of rotation for integrals involving p and f orbitals with respect to theta
  matrixE transf(21,21);
  if (!(rotate)) {
    transf(3,2) = 1.0;       //<p(-pi)|f(-pi)>,<f(-pi)|p(-pi)> pi part
    transf(11,1) = 1.0;      //<p(sigma)|f(sigma)>,<f(sigma)|p(sigma)> sigma part
    transf(19,2) = 1.0;      //<p(pi)|f(pi)>,<f(pi)|p(pi)> pi part
  }
  else if (rotate) {
    matrixE transf_sf = SFtransf(cost,sint,cosp,sinp);
    matrixE transf_sp = SPtransf(cost,sint,cosp,sinp);
    matrixE transf_sf_dt = SFtransf_dt(cost,sint,cosp,sinp);
    matrixE transf_sp_dt = SPtransf_dt(cost,sint,cosp,sinp);
    matrixE transf_sf_dp = SFtransf_dp(cost,sint,cosp,sinp);
    matrixE transf_sp_dp = SPtransf_dp(cost,sint,cosp,sinp);
    matrixE transf_sf_dtdp = SFtransf_dtdp(cost,sint,cosp,sinp);
    matrixE transf_sp_dtdp = SPtransf_dtdp(cost,sint,cosp,sinp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sf.rows(); ++idx1) {
      for (size_t idx2 = 0; idx2 < transf_sp.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 4; ++idx3) {
          transf(counter + 1,idx3) = transf_sf_dtdp(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx3)*transf_sp_dtdp(idx2 + 1,idx3);
          transf(counter + 1,idx3) += transf_sf_dt(idx1 + 1,idx3)*transf_sp_dp(idx2 + 1,idx3) + transf_sf_dp(idx1 + 1,idx3)*transf_sp_dt(idx2 + 1,idx3);
          for (size_t idx4 = 1; idx4 < 8; ++idx4) {
            if (idx4 == idx3) continue;
            transf(counter + 1,3 + cnt) = transf_sf_dtdp(idx1 + 1,idx4)*transf_sp(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx4)*transf_sp_dtdp(idx2 + 1,idx3);
            transf(counter + 1,3 + cnt) += transf_sf_dt(idx1 + 1,idx4)*transf_sp_dp(idx2 + 1,idx3) + transf_sf_dp(idx1 + 1,idx4)*transf_sp_dt(idx2 + 1,idx3);
            ++cnt;
          }
        }
      }
    }
  }
  return transf;
}
matrixE PFtransf_dp2(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that calculates the second-derivatives of rotation for integrals involving p and f orbitals with respect to phi
  matrixE transf(21,21);
  if (!(rotate)) {
    transf(3,2) = 1.0;       //<p(-pi)|f(-pi)>,<f(-pi)|p(-pi)> pi part
    transf(11,1) = 1.0;      //<p(sigma)|f(sigma)>,<f(sigma)|p(sigma)> sigma part
    transf(19,2) = 1.0;      //<p(pi)|f(pi)>,<f(pi)|p(pi)> pi part
  }
  else if (rotate) {
    matrixE transf_sf = SFtransf(cost,sint,cosp,sinp);
    matrixE transf_sp = SPtransf(cost,sint,cosp,sinp);
    matrixE transf_sf_dp = SFtransf_dp(cost,sint,cosp,sinp);
    matrixE transf_sp_dp = SPtransf_dp(cost,sint,cosp,sinp);
    matrixE transf_sf_dp2 = SFtransf_dp2(cost,sint,cosp,sinp);
    matrixE transf_sp_dp2 = SPtransf_dp2(cost,sint,cosp,sinp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sf.rows(); ++idx1) {
      for (size_t idx2 = 0; idx2 < transf_sp.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 4; ++idx3) {
          transf(counter + 1,idx3) = transf_sf_dp2(idx1 + 1,idx3)*transf_sp(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx3)*transf_sp_dp2(idx2 + 1,idx3);
          transf(counter + 1,idx3) += 2.0*transf_sf_dp(idx1 + 1,idx3)*transf_sp_dp(idx2 + 1,idx3);
          for (size_t idx4 = 1; idx4 < 8; ++idx4) {
            if (idx4 == idx3) continue;
            transf(counter + 1,3 + cnt) = transf_sf_dp2(idx1 + 1,idx4)*transf_sp(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx4)*transf_sp_dp2(idx2 + 1,idx3);
            transf(counter + 1,3 + cnt) += 2.0*transf_sf_dp(idx1 + 1,idx4)*transf_sp_dp(idx2 + 1,idx3);
            ++cnt;
          }
        }
      }
    }
  }
  return transf;
}
matrixE DFtransf(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that calculates the rotation for integrals involving d and f orbitals
  matrixE transf(35,35);
  if (!(rotate)) { 
    transf(2,3) = 1.0;        //<d(-delta)|f(-delta)>,<f(-delta)|d(-delta)> delta part
    transf(10,2) = 1.0;       //<d(-pi)|f(-pi)>,<f(-pi)|d(-pi)> pi part
    transf(18,1) = 1.0;       //<d(sigma)|f(sigma)>,<f(sigma)|d(sigma)> sigma part
    transf(26,2) = 1.0;       //<d(pi)|f(pi)>,<f(pi)|d(pi)> pi part
    transf(34,3) = 1.0;       //<d(delta)|f(delta)>,<f(delta)|d(delta)> delta part
  }
  else if (rotate) {
    matrixE transf_sf = SFtransf(cost,sint,cosp,sinp);
    matrixE transf_sd = SDtransf(cost,sint,cosp,sinp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sf.rows(); ++idx1) {
      for (size_t idx2 = 0; idx2 < transf_sd.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 6; ++idx3) {
          transf(counter + 1,idx3) = transf_sf(idx1 + 1,idx3)*transf_sd(idx2 + 1,idx3);
          for (size_t idx4 = 1; idx4 < 8; ++idx4) {
            if (idx4 == idx3) continue;
            transf(counter + 1,5 + cnt) = transf_sf(idx1 + 1,idx4)*transf_sd(idx2 + 1,idx3);
            ++cnt;
          }
        }
      }
    }
  }
  return transf;
}
matrixE DFtransf_dt(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that calculates the derivatives of rotation for integrals involving d and f orbitals with respect to theta
  matrixE transf(35,35);
  if (!(rotate)) { 
    transf(2,3) = 1.0;        //<d(-delta)|f(-delta)>,<f(-delta)|d(-delta)> delta part
    transf(10,2) = 1.0;       //<d(-pi)|f(-pi)>,<f(-pi)|d(-pi)> pi part
    transf(18,1) = 1.0;       //<d(sigma)|f(sigma)>,<f(sigma)|d(sigma)> sigma part
    transf(26,2) = 1.0;       //<d(pi)|f(pi)>,<f(pi)|d(pi)> pi part
    transf(34,3) = 1.0;       //<d(delta)|f(delta)>,<f(delta)|d(delta)> delta part
  }
  else if (rotate) {
    matrixE transf_sf = SFtransf(cost,sint,cosp,sinp);
    matrixE transf_sd = SDtransf(cost,sint,cosp,sinp);
    matrixE transf_sf_dt = SFtransf_dt(cost,sint,cosp,sinp);
    matrixE transf_sd_dt = SDtransf_dt(cost,sint,cosp,sinp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sf.rows(); ++idx1) {
      for (size_t idx2 = 0; idx2 < transf_sd.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 6; ++idx3) {
          transf(counter + 1,idx3) = transf_sf_dt(idx1 + 1,idx3)*transf_sd(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx3)*transf_sd_dt(idx2 + 1,idx3);
          for (size_t idx4 = 1; idx4 < 8; ++idx4) {
            if (idx4 == idx3) continue;
            transf(counter + 1,5 + cnt) = transf_sf_dt(idx1 + 1,idx4)*transf_sd(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx4)*transf_sd_dt(idx2 + 1,idx3);
            ++cnt;
          }
        }
      }
    }
  }
  return transf;
}
matrixE DFtransf_dp(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that calculates the derivatives of rotation for integrals involving d and f orbitals with respect to phi
  matrixE transf(35,35);
  if (!(rotate)) { 
    transf(2,3) = 1.0;        //<d(-delta)|f(-delta)>,<f(-delta)|d(-delta)> delta part
    transf(10,2) = 1.0;       //<d(-pi)|f(-pi)>,<f(-pi)|d(-pi)> pi part
    transf(18,1) = 1.0;       //<d(sigma)|f(sigma)>,<f(sigma)|d(sigma)> sigma part
    transf(26,2) = 1.0;       //<d(pi)|f(pi)>,<f(pi)|d(pi)> pi part
    transf(34,3) = 1.0;       //<d(delta)|f(delta)>,<f(delta)|d(delta)> delta part
  }
  else if (rotate) {
    matrixE transf_sf = SFtransf(cost,sint,cosp,sinp);
    matrixE transf_sd = SDtransf(cost,sint,cosp,sinp);
    matrixE transf_sf_dp = SFtransf_dp(cost,sint,cosp,sinp);
    matrixE transf_sd_dp = SDtransf_dp(cost,sint,cosp,sinp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sf.rows(); ++idx1) {
      for (size_t idx2 = 0; idx2 < transf_sd.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 6; ++idx3) {
          transf(counter + 1,idx3) = transf_sf_dp(idx1 + 1,idx3)*transf_sd(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx3)*transf_sd_dp(idx2 + 1,idx3);
          for (size_t idx4 = 1; idx4 < 8; ++idx4) {
            if (idx4 == idx3) continue;
            transf(counter + 1,5 + cnt) = transf_sf_dp(idx1 + 1,idx4)*transf_sd(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx4)*transf_sd_dp(idx2 + 1,idx3);
            ++cnt;
          }
        }
      }
    }
  }
  return transf;
}
matrixE DFtransf_dt2(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that calculates the second-derivatives of rotation for integrals involving d and f orbitals with respect to theta
  matrixE transf(35,35);
  if (!(rotate)) { 
    transf(2,3) = 1.0;        //<d(-delta)|f(-delta)>,<f(-delta)|d(-delta)> delta part
    transf(10,2) = 1.0;       //<d(-pi)|f(-pi)>,<f(-pi)|d(-pi)> pi part
    transf(18,1) = 1.0;       //<d(sigma)|f(sigma)>,<f(sigma)|d(sigma)> sigma part
    transf(26,2) = 1.0;       //<d(pi)|f(pi)>,<f(pi)|d(pi)> pi part
    transf(34,3) = 1.0;       //<d(delta)|f(delta)>,<f(delta)|d(delta)> delta part
  }
  else if (rotate) {
    matrixE transf_sf = SFtransf(cost,sint,cosp,sinp);
    matrixE transf_sd = SDtransf(cost,sint,cosp,sinp);
    matrixE transf_sf_dt = SFtransf_dt(cost,sint,cosp,sinp);
    matrixE transf_sd_dt = SDtransf_dt(cost,sint,cosp,sinp);
    matrixE transf_sf_dt2 = SFtransf_dt2(cost,sint,cosp,sinp);
    matrixE transf_sd_dt2 = SDtransf_dt2(cost,sint,cosp,sinp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sf.rows(); ++idx1) {
      for (size_t idx2 = 0; idx2 < transf_sd.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 6; ++idx3) {
          transf(counter + 1,idx3) = transf_sf_dt2(idx1 + 1,idx3)*transf_sd(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx3)*transf_sd_dt2(idx2 + 1,idx3);
          transf(counter + 1,idx3) += 2.0*transf_sf_dt(idx1 + 1,idx3)*transf_sd_dt(idx2 + 1,idx3);
          for (size_t idx4 = 1; idx4 < 8; ++idx4) {
            if (idx4 == idx3) continue;
            transf(counter + 1,5 + cnt) = transf_sf_dt2(idx1 + 1,idx4)*transf_sd(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx4)*transf_sd_dt2(idx2 + 1,idx3);
            transf(counter + 1,5 + cnt) += 2.0*transf_sf_dt(idx1 + 1,idx4)*transf_sd_dt(idx2 + 1,idx3);
            ++cnt;
          }
        }
      }
    }
  }
  return transf;
}
matrixE DFtransf_dtdp(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that calculates the second-derivatives of rotation for integrals involving d and f orbitals with respect to theta and phi
  matrixE transf(35,35);
  if (!(rotate)) { 
    transf(2,3) = 1.0;        //<d(-delta)|f(-delta)>,<f(-delta)|d(-delta)> delta part
    transf(10,2) = 1.0;       //<d(-pi)|f(-pi)>,<f(-pi)|d(-pi)> pi part
    transf(18,1) = 1.0;       //<d(sigma)|f(sigma)>,<f(sigma)|d(sigma)> sigma part
    transf(26,2) = 1.0;       //<d(pi)|f(pi)>,<f(pi)|d(pi)> pi part
    transf(34,3) = 1.0;       //<d(delta)|f(delta)>,<f(delta)|d(delta)> delta part
  }
  else if (rotate) {
    matrixE transf_sf = SFtransf(cost,sint,cosp,sinp);
    matrixE transf_sd = SDtransf(cost,sint,cosp,sinp);
    matrixE transf_sf_dt = SFtransf_dt(cost,sint,cosp,sinp);
    matrixE transf_sd_dt = SDtransf_dt(cost,sint,cosp,sinp);
    matrixE transf_sf_dp = SFtransf_dp(cost,sint,cosp,sinp);
    matrixE transf_sd_dp = SDtransf_dp(cost,sint,cosp,sinp);
    matrixE transf_sf_dtdp = SFtransf_dtdp(cost,sint,cosp,sinp);
    matrixE transf_sd_dtdp = SDtransf_dtdp(cost,sint,cosp,sinp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sf.rows(); ++idx1) {
      for (size_t idx2 = 0; idx2 < transf_sd.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 6; ++idx3) {
          transf(counter + 1,idx3) = transf_sf_dtdp(idx1 + 1,idx3)*transf_sd(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx3)*transf_sd_dtdp(idx2 + 1,idx3);
          transf(counter + 1,idx3) += transf_sf_dt(idx1 + 1,idx3)*transf_sd_dp(idx2 + 1,idx3) + transf_sf_dp(idx1 + 1,idx3)*transf_sd_dt(idx2 + 1,idx3);
          for (size_t idx4 = 1; idx4 < 8; ++idx4) {
            if (idx4 == idx3) continue;
            transf(counter + 1,5 + cnt) = transf_sf_dtdp(idx1 + 1,idx4)*transf_sd(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx4)*transf_sd_dtdp(idx2 + 1,idx3);
            transf(counter + 1,5 + cnt) += transf_sf_dt(idx1 + 1,idx4)*transf_sd_dp(idx2 + 1,idx3) + transf_sf_dp(idx1 + 1,idx4)*transf_sd_dt(idx2 + 1,idx3);
            ++cnt;
          }
        }
      }
    }
  }
  return transf;
}
matrixE DFtransf_dp2(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that calculates the second-derivatives of rotation for integrals involving d and f orbitals with respect to phi
  matrixE transf(35,35);
  if (!(rotate)) { 
    transf(2,3) = 1.0;        //<d(-delta)|f(-delta)>,<f(-delta)|d(-delta)> delta part
    transf(10,2) = 1.0;       //<d(-pi)|f(-pi)>,<f(-pi)|d(-pi)> pi part
    transf(18,1) = 1.0;       //<d(sigma)|f(sigma)>,<f(sigma)|d(sigma)> sigma part
    transf(26,2) = 1.0;       //<d(pi)|f(pi)>,<f(pi)|d(pi)> pi part
    transf(34,3) = 1.0;       //<d(delta)|f(delta)>,<f(delta)|d(delta)> delta part
  }
  else if (rotate) {
    matrixE transf_sf = SFtransf(cost,sint,cosp,sinp);
    matrixE transf_sd = SDtransf(cost,sint,cosp,sinp);
    matrixE transf_sf_dp = SFtransf_dp(cost,sint,cosp,sinp);
    matrixE transf_sd_dp = SDtransf_dp(cost,sint,cosp,sinp);
    matrixE transf_sf_dp2 = SFtransf_dp2(cost,sint,cosp,sinp);
    matrixE transf_sd_dp2 = SDtransf_dp2(cost,sint,cosp,sinp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sf.rows(); ++idx1) {
      for (size_t idx2 = 0; idx2 < transf_sd.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 6; ++idx3) {
          transf(counter + 1,idx3) = transf_sf_dp2(idx1 + 1,idx3)*transf_sd(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx3)*transf_sd_dp2(idx2 + 1,idx3);
          transf(counter + 1,idx3) += 2.0*transf_sf_dp(idx1 + 1,idx3)*transf_sd_dp(idx2 + 1,idx3);
          for (size_t idx4 = 1; idx4 < 8; ++idx4) {
            if (idx4 == idx3) continue;
            transf(counter + 1,5 + cnt) = transf_sf_dp2(idx1 + 1,idx4)*transf_sd(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx4)*transf_sd_dp2(idx2 + 1,idx3);
            transf(counter + 1,5 + cnt) += 2.0*transf_sf_dp(idx1 + 1,idx4)*transf_sd_dp(idx2 + 1,idx3);
            ++cnt;
          }
        }
      }
    }
  }
  return transf;
}
matrixE FFtransf(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that calculates the rotation for integrals involving f orbitals
  matrixE transf(28,28);
  if (!(rotate)) {
    transf(1,4) = 1.0;        //<f(-phi)|f(-phi)> phi part
    transf(8,3) = 1.0;        //<f(-delta)|f(-delta) delta part
    transf(14,2) = 1.0;       //<f(-pi)|f(-pi)pi part
    transf(19,1) = 1.0;       //<f(sigma)|f(sigma)> sigma part
    transf(23,2) = 1.0;       //<f(pi)|f(pi)> pi part
    transf(26,3) = 1.0;       //<f(delta)|f(delta)> delta part
    transf(28,4) = 1.0;       //<f(phi)|f(phi)> phi part
  }
  else if (rotate) {
    matrixE transf_sf = SFtransf(cost,sint,cosp,sinp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sf.rows(); ++idx1) {
      for (size_t idx2 = idx1; idx2 < transf_sf.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 8; ++idx3) {
          transf(counter + 1,idx3) = transf_sf(idx1 + 1,idx3)*transf_sf(idx2 + 1,idx3);
          for (size_t idx4 = idx3 + 1; idx4 < 8; ++idx4, ++cnt) {
            transf(counter + 1,7 + cnt) = transf_sf(idx1 + 1,idx3)*transf_sf(idx2 + 1,idx4) + transf_sf(idx1 + 1,idx4)*transf_sf(idx2 + 1,idx3);
          }
        }
      }
    }
  }
  return transf;
}
matrixE FFtransf_dt(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that calculates the derivatives of rotation for integrals involving f orbitals with respect to theta
  matrixE transf(28,28);
  if (!(rotate)) {
    transf(1,4) = 1.0;        //<f(-phi)|f(-phi)> phi part
    transf(8,3) = 1.0;        //<f(-delta)|f(-delta) delta part
    transf(14,2) = 1.0;       //<f(-pi)|f(-pi)pi part
    transf(19,1) = 1.0;       //<f(sigma)|f(sigma)> sigma part
    transf(23,2) = 1.0;       //<f(pi)|f(pi)> pi part
    transf(26,3) = 1.0;       //<f(delta)|f(delta)> delta part
    transf(28,4) = 1.0;       //<f(phi)|f(phi)> phi part
  }
  else if (rotate) {
    matrixE transf_sf = SFtransf(cost,sint,cosp,sinp);
    matrixE transf_sf_dt = SFtransf_dt(cost,sint,cosp,sinp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sf.rows(); ++idx1) {
      for (size_t idx2 = idx1; idx2 < transf_sf.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 8; ++idx3) {
          transf(counter + 1,idx3) = transf_sf_dt(idx1 + 1,idx3)*transf_sf(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx3)*transf_sf_dt(idx2 + 1,idx3);
          for (size_t idx4 = idx3 + 1; idx4 < 8; ++idx4, ++cnt) {
            transf(counter + 1,7 + cnt) = transf_sf_dt(idx1 + 1,idx3)*transf_sf(idx2 + 1,idx4) + transf_sf(idx1 + 1,idx3)*transf_sf_dt(idx2 + 1,idx4);
            transf(counter + 1,7 + cnt) += transf_sf_dt(idx1 + 1,idx4)*transf_sf(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx4)*transf_sf_dt(idx2 + 1,idx3);
          }
        }
      }
    }
  }
  return transf;
}
matrixE FFtransf_dp(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that calculates the derivatives of rotation for integrals involving f orbitals with respect to phi
  matrixE transf(28,28);
  if (!(rotate)) {
    transf(1,4) = 1.0;        //<f(-phi)|f(-phi)> phi part
    transf(8,3) = 1.0;        //<f(-delta)|f(-delta) delta part
    transf(14,2) = 1.0;       //<f(-pi)|f(-pi)pi part
    transf(19,1) = 1.0;       //<f(sigma)|f(sigma)> sigma part
    transf(23,2) = 1.0;       //<f(pi)|f(pi)> pi part
    transf(26,3) = 1.0;       //<f(delta)|f(delta)> delta part
    transf(28,4) = 1.0;       //<f(phi)|f(phi)> phi part
  }
  else if (rotate) {
    matrixE transf_sf = SFtransf(cost,sint,cosp,sinp);
    matrixE transf_sf_dp = SFtransf_dp(cost,sint,cosp,sinp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sf.rows(); ++idx1) {
      for (size_t idx2 = idx1; idx2 < transf_sf.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 8; ++idx3) {
          transf(counter + 1,idx3) = transf_sf_dp(idx1 + 1,idx3)*transf_sf(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx3)*transf_sf_dp(idx2 + 1,idx3);
          for (size_t idx4 = idx3 + 1; idx4 < 8; ++idx4, ++cnt) {
            transf(counter + 1,7 + cnt) = transf_sf_dp(idx1 + 1,idx3)*transf_sf(idx2 + 1,idx4) + transf_sf(idx1 + 1,idx3)*transf_sf_dp(idx2 + 1,idx4);
            transf(counter + 1,7 + cnt) += transf_sf_dp(idx1 + 1,idx4)*transf_sf(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx4)*transf_sf_dp(idx2 + 1,idx3);
          }
        }
      }
    }
  }
  return transf;
}
matrixE FFtransf_dt2(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that calculates the second-derivatives of rotation for integrals involving f orbitals with respect to theta
  matrixE transf(28,28);
  if (!(rotate)) {
    transf(1,4) = 1.0;        //<f(-phi)|f(-phi)> phi part
    transf(8,3) = 1.0;        //<f(-delta)|f(-delta) delta part
    transf(14,2) = 1.0;       //<f(-pi)|f(-pi)pi part
    transf(19,1) = 1.0;       //<f(sigma)|f(sigma)> sigma part
    transf(23,2) = 1.0;       //<f(pi)|f(pi)> pi part
    transf(26,3) = 1.0;       //<f(delta)|f(delta)> delta part
    transf(28,4) = 1.0;       //<f(phi)|f(phi)> phi part
  }
  else if (rotate) {
    matrixE transf_sf = SFtransf(cost,sint,cosp,sinp);
    matrixE transf_sf_dt = SFtransf_dt(cost,sint,cosp,sinp);
    matrixE transf_sf_dt2 = SFtransf_dt2(cost,sint,cosp,sinp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sf.rows(); ++idx1) {
      for (size_t idx2 = idx1; idx2 < transf_sf.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 8; ++idx3) {
          transf(counter + 1,idx3) = transf_sf_dt2(idx1 + 1,idx3)*transf_sf(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx3)*transf_sf_dt2(idx2 + 1,idx3);
          transf(counter + 1,idx3) += 2.0*transf_sf_dt(idx1 + 1,idx3)*transf_sf_dt(idx2 + 1,idx3);
          for (size_t idx4 = idx3 + 1; idx4 < 8; ++idx4, ++cnt) {
            transf(counter + 1,7 + cnt) = transf_sf_dt2(idx1 + 1,idx3)*transf_sf(idx2 + 1,idx4) + transf_sf(idx1 + 1,idx3)*transf_sf_dt2(idx2 + 1,idx4);
            transf(counter + 1,7 + cnt) += 2.0*(transf_sf_dt(idx1 + 1,idx3)*transf_sf_dt(idx2 + 1,idx4) + transf_sf_dt(idx1 + 1,idx4)*transf_sf_dt(idx2 + 1,idx3));
            transf(counter + 1,7 + cnt) += transf_sf_dt2(idx1 + 1,idx4)*transf_sf(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx4)*transf_sf_dt2(idx2 + 1,idx3);
          }
        }
      }
    }
  }
  return transf;
}
matrixE FFtransf_dtdp(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that calculates the second-derivatives of rotation for integrals involving f orbitals with respect to theta and phi
  matrixE transf(28,28);
  if (!(rotate)) {
    transf(1,4) = 1.0;        //<f(-phi)|f(-phi)> phi part
    transf(8,3) = 1.0;        //<f(-delta)|f(-delta) delta part
    transf(14,2) = 1.0;       //<f(-pi)|f(-pi)pi part
    transf(19,1) = 1.0;       //<f(sigma)|f(sigma)> sigma part
    transf(23,2) = 1.0;       //<f(pi)|f(pi)> pi part
    transf(26,3) = 1.0;       //<f(delta)|f(delta)> delta part
    transf(28,4) = 1.0;       //<f(phi)|f(phi)> phi part
  }
  else if (rotate) {
    matrixE transf_sf = SFtransf(cost,sint,cosp,sinp);
    matrixE transf_sf_dt = SFtransf_dt(cost,sint,cosp,sinp);
    matrixE transf_sf_dp = SFtransf_dp(cost,sint,cosp,sinp);
    matrixE transf_sf_dtdp = SFtransf_dtdp(cost,sint,cosp,sinp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sf.rows(); ++idx1) {
      for (size_t idx2 = idx1; idx2 < transf_sf.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 8; ++idx3) {
          transf(counter + 1,idx3) = transf_sf_dtdp(idx1 + 1,idx3)*transf_sf(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx3)*transf_sf_dtdp(idx2 + 1,idx3);
          transf(counter + 1,idx3) += transf_sf_dt(idx1 + 1,idx3)*transf_sf_dp(idx2 + 1,idx3) + transf_sf_dp(idx1 + 1,idx3)*transf_sf_dt(idx2 + 1,idx3);
          for (size_t idx4 = idx3 + 1; idx4 < 8; ++idx4, ++cnt) {
            transf(counter + 1,7 + cnt) = transf_sf_dtdp(idx1 + 1,idx3)*transf_sf(idx2 + 1,idx4) + transf_sf(idx1 + 1,idx3)*transf_sf_dtdp(idx2 + 1,idx4);
            transf(counter + 1,7 + cnt) += transf_sf_dt(idx1 + 1,idx3)*transf_sf_dp(idx2 + 1,idx4) + transf_sf_dp(idx1 + 1,idx3)*transf_sf_dt(idx2 + 1,idx4);
            transf(counter + 1,7 + cnt) += transf_sf_dtdp(idx1 + 1,idx4)*transf_sf(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx4)*transf_sf_dtdp(idx2 + 1,idx3);
            transf(counter + 1,7 + cnt) += transf_sf_dt(idx1 + 1,idx4)*transf_sf_dp(idx2 + 1,idx3) + transf_sf_dp(idx1 + 1,idx4)*transf_sf_dt(idx2 + 1,idx3);
          }
        }
      }
    }
  }
  return transf;
}
matrixE FFtransf_dp2(double cost, double sint, double cosp, double sinp, bool rotate = true) {
  //function that calculates the second-derivatives of rotation for integrals involving f orbitals with respect to phi
  matrixE transf(28,28);
  if (!(rotate)) {
    transf(1,4) = 1.0;        //<f(-phi)|f(-phi)> phi part
    transf(8,3) = 1.0;        //<f(-delta)|f(-delta) delta part
    transf(14,2) = 1.0;       //<f(-pi)|f(-pi)pi part
    transf(19,1) = 1.0;       //<f(sigma)|f(sigma)> sigma part
    transf(23,2) = 1.0;       //<f(pi)|f(pi)> pi part
    transf(26,3) = 1.0;       //<f(delta)|f(delta)> delta part
    transf(28,4) = 1.0;       //<f(phi)|f(phi)> phi part
  }
  else if (rotate) {
    matrixE transf_sf = SFtransf(cost,sint,cosp,sinp);
    matrixE transf_sf_dp = SFtransf_dp(cost,sint,cosp,sinp);
    matrixE transf_sf_dp2 = SFtransf_dp2(cost,sint,cosp,sinp);
    int counter = 0;
    int cnt;
    for (size_t idx1 = 0; idx1 < transf_sf.rows(); ++idx1) {
      for (size_t idx2 = idx1; idx2 < transf_sf.rows(); ++idx2, ++counter) {
        cnt = 1;
        for (size_t idx3 = 1; idx3 < 8; ++idx3) {
          transf(counter + 1,idx3) = transf_sf_dp2(idx1 + 1,idx3)*transf_sf(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx3)*transf_sf_dp2(idx2 + 1,idx3);
          transf(counter + 1,idx3) += 2.0*transf_sf_dp(idx1 + 1,idx3)*transf_sf_dp(idx2 + 1,idx3);
          for (size_t idx4 = idx3 + 1; idx4 < 8; ++idx4, ++cnt) {
            transf(counter + 1,7 + cnt) = transf_sf_dp2(idx1 + 1,idx3)*transf_sf(idx2 + 1,idx4) + transf_sf(idx1 + 1,idx3)*transf_sf_dp2(idx2 + 1,idx4);
            transf(counter + 1,7 + cnt) += 2.0*(transf_sf_dp(idx1 + 1,idx3)*transf_sf_dp(idx2 + 1,idx4) + transf_sf_dp(idx1 + 1,idx4)*transf_sf_dp(idx2 + 1,idx3));
            transf(counter + 1,7 + cnt) += transf_sf_dp2(idx1 + 1,idx4)*transf_sf(idx2 + 1,idx3) + transf_sf(idx1 + 1,idx4)*transf_sf_dp2(idx2 + 1,idx3);
          }
        }
      }
    }
  }
  return transf;
}

#endif //_Overlap_Rotations_