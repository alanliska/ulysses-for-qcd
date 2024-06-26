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

#ifndef _random_double_
#define _random_double_

//description:
//random number generators

double fRandom(double min, double max) {
  //function that generates a random double between limits min and max
  double frand = double(rand())/double(RAND_MAX);
  return min + frand*(max - min);
}
int iRandom(int min, int max) {
  //function generating an integer between limits min and max
  int irand = int(rand())%max;
  return min + irand;
}
int iRandom01() {
  //function generating randomly 0 or 1
  int irand = rand()%2;
  return irand;
}
double GaussianRandom(double mean, double std, double thresholdzero = 1.0e-14) {
  //calculation of Gaussian distributed values
  double w1;
  double w2;
  double r2;
resample:
  w1 = fRandom(-1.0,1.0);
  w2 = fRandom(-1.0,1.0);
  r2 = w1*w1 + w2*w2;
  if ((r2 > 1.0)||(r2 < thresholdzero)) {goto resample;}
  w2 = sqrt(-2.0*log(r2)/r2);
  return mean + w1*w2*std;
}

#endif //_random_double_
