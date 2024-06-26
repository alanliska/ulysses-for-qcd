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

#ifndef _Slater_Type_Orbitals_
#define _Slater_Type_Orbitals_
#include <vector>
#include <string>
#include "OvRot.hpp"
#include "1CenterOvSTO.hpp"
#include "2CenterOvSTO.hpp"
#include "2CenterOvGTO.hpp"
#include "../math/MatrixPackage.hpp"
#include "../math/FunctionPackage.hpp"

//description:
//base class related to STO orbitals (and GTOs based on these)
//auxliary file for the calculation of quantities associated to Slater orbitals

class STOBSet {
public:
  STOBSet() {}
  STOBSet(const STOBSet & other) {}
  STOBSet & operator= (const STOBSet & other) {
    return *this;
  }
  ~STOBSet() {}
  virtual std::vector<double> GetExponents(size_t atomnum) {
    //getter for the exponents of a given atom
    std::vector<double> exponents;
    exponents.push_back(0.0);
    return exponents;
  }
  virtual matrix<size_t> QuantumNumbers(size_t atomnum) {
    //getter for quantum numbers of one atom
    matrix<size_t> quantumnumb(1,2);
    return quantumnumb;
  }
  //direct STO stuff
  virtual size_t GetNAOs(size_t atomnum) {
    //function to calculate the number of AOs used for a given atom
    size_t orbitals = 1;
    return orbitals;
  }
  virtual double ZeroOverlap(size_t atm) {return 10000000000000000.0;}
  void OV2Center(matrixE & overlaps, size_t atom1, size_t atom2, double dist, std::vector<double> orientation, bool rotate = true) {
    //main driver for two-center overlap
    //used for blocks of the overlap
    size_t nao1 = GetNAOs(atom1);
    size_t nao2 = GetNAOs(atom2);
    std::vector<double> exponents1 = GetExponents(atom1);
    matrix<size_t> QN1 = QuantumNumbers(atom1);
    std::vector<double> exponents2 = GetExponents(atom2);
    matrix<size_t> QN2 = QuantumNumbers(atom2);
    TwoCenterOv(overlaps,nao1,nao2,dist,exponents1,QN1,exponents2,QN2,orientation,rotate);
  }
  void OV2Center_dR(matrixE & overlaps, size_t atom1, size_t atom2, double dist, std::vector<double> orientation, bool rotate = true) {
    //main driver for first-derivatives of two-center overlap with respect to internuclear distance
    //used for blocks of the overlap derivative
    size_t nao1 = GetNAOs(atom1);
    size_t nao2 = GetNAOs(atom2);
    std::vector<double> exponents1 = GetExponents(atom1);
    matrix<size_t> QN1 = QuantumNumbers(atom1);
    std::vector<double> exponents2 = GetExponents(atom2);
    matrix<size_t> QN2 = QuantumNumbers(atom2);
    TwoCenterOv_dR(overlaps,nao1,nao2,dist,exponents1,QN1,exponents2,QN2,orientation,rotate);
  }
  void OV2Center_dt(matrixE & overlaps, size_t atom1, size_t atom2, double dist, std::vector<double> orientation, bool rotate = true) {
    //main driver for first-derivatives of two-center overlap with respect to angle theta
    //used for blocks of the overlap derivative
    size_t nao1 = GetNAOs(atom1);
    size_t nao2 = GetNAOs(atom2);
    std::vector<double> exponents1 = GetExponents(atom1);
    matrix<size_t> QN1 = QuantumNumbers(atom1);
    std::vector<double> exponents2 = GetExponents(atom2);
    matrix<size_t> QN2 = QuantumNumbers(atom2);
    TwoCenterOv_dt(overlaps,nao1,nao2,dist,exponents1,QN1,exponents2,QN2,orientation,rotate);
  }
  void OV2Center_dp(matrixE & overlaps, size_t atom1, size_t atom2, double dist, std::vector<double> orientation, bool rotate = true) {
    //main driver for first-derivatives of two-center overlap with respect to angle phi
    //used for blocks of the overlap derivative
    size_t nao1 = GetNAOs(atom1);
    size_t nao2 = GetNAOs(atom2);
    std::vector<double> exponents1 = GetExponents(atom1);
    matrix<size_t> QN1 = QuantumNumbers(atom1);
    std::vector<double> exponents2 = GetExponents(atom2);
    matrix<size_t> QN2 = QuantumNumbers(atom2);
    TwoCenterOv_dp(overlaps,nao1,nao2,dist,exponents1,QN1,exponents2,QN2,orientation,rotate);
  }
  void OV2Center_dR2(matrixE & overlaps, size_t atom1, size_t atom2, double dist, std::vector<double> orientation, bool rotate = true) {
    //main driver for second-derivatives of two-center overlap with respect to internuclear distance
    //used for blocks of the overlap derivative
    size_t nao1 = GetNAOs(atom1);
    size_t nao2 = GetNAOs(atom2);
    std::vector<double> exponents1 = GetExponents(atom1);
    matrix<size_t> QN1 = QuantumNumbers(atom1);
    std::vector<double> exponents2 = GetExponents(atom2);
    matrix<size_t> QN2 = QuantumNumbers(atom2);
    TwoCenterOv_dR2(overlaps,nao1,nao2,dist,exponents1,QN1,exponents2,QN2,orientation,rotate);
  }
  void OV2Center_dRdt(matrixE & overlaps, size_t atom1, size_t atom2, double dist, std::vector<double> orientation, bool rotate = true) {
    //main driver for second-derivatives of two-center overlap with respect to internuclear distance and angle theta
    //used for blocks of the overlap derivative
    size_t nao1 = GetNAOs(atom1);
    size_t nao2 = GetNAOs(atom2);
    std::vector<double> exponents1 = GetExponents(atom1);
    matrix<size_t> QN1 = QuantumNumbers(atom1);
    std::vector<double> exponents2 = GetExponents(atom2);
    matrix<size_t> QN2 = QuantumNumbers(atom2);
    TwoCenterOv_dRdt(overlaps,nao1,nao2,dist,exponents1,QN1,exponents2,QN2,orientation,rotate);
  }
  void OV2Center_dRdp(matrixE & overlaps, size_t atom1, size_t atom2, double dist, std::vector<double> orientation, bool rotate = true) {
    //main driver for second-derivatives of two-center overlap with respect to internuclear distance and angle phi
    //used for blocks of the overlap derivative
    size_t nao1 = GetNAOs(atom1);
    size_t nao2 = GetNAOs(atom2);
    std::vector<double> exponents1 = GetExponents(atom1);
    matrix<size_t> QN1 = QuantumNumbers(atom1);
    std::vector<double> exponents2 = GetExponents(atom2);
    matrix<size_t> QN2 = QuantumNumbers(atom2);
    TwoCenterOv_dRdp(overlaps,nao1,nao2,dist,exponents1,QN1,exponents2,QN2,orientation,rotate);
  }
  void OV2Center_dt2(matrixE & overlaps, size_t atom1, size_t atom2, double dist, std::vector<double> orientation, bool rotate = true) {
    //main driver for second-derivatives of two-center overlap with respect to angle theta
    //used for blocks of the overlap derivative
    size_t nao1 = GetNAOs(atom1);
    size_t nao2 = GetNAOs(atom2);
    std::vector<double> exponents1 = GetExponents(atom1);
    matrix<size_t> QN1 = QuantumNumbers(atom1);
    std::vector<double> exponents2 = GetExponents(atom2);
    matrix<size_t> QN2 = QuantumNumbers(atom2);
    TwoCenterOv_dt2(overlaps,nao1,nao2,dist,exponents1,QN1,exponents2,QN2,orientation,rotate);
  }
  void OV2Center_dp2(matrixE & overlaps, size_t atom1, size_t atom2, double dist, std::vector<double> orientation, bool rotate = true) {
    //main driver for second-derivatives of two-center overlap with respect to angle phi
    //used for blocks of the overlap derivative
    size_t nao1 = GetNAOs(atom1);
    size_t nao2 = GetNAOs(atom2);
    std::vector<double> exponents1 = GetExponents(atom1);
    matrix<size_t> QN1 = QuantumNumbers(atom1);
    std::vector<double> exponents2 = GetExponents(atom2);
    matrix<size_t> QN2 = QuantumNumbers(atom2);
    TwoCenterOv_dp2(overlaps,nao1,nao2,dist,exponents1,QN1,exponents2,QN2,orientation,rotate);
  }
  void OV2Center_dtdp(matrixE & overlaps, size_t atom1, size_t atom2, double dist, std::vector<double> orientation, bool rotate = true) {
    //main driver for second-derivatives of two-center overlap with respect to angles theta and phi
    //used for blocks of the overlap derivative
    size_t nao1 = GetNAOs(atom1);
    size_t nao2 = GetNAOs(atom2);
    std::vector<double> exponents1 = GetExponents(atom1);
    matrix<size_t> QN1 = QuantumNumbers(atom1);
    std::vector<double> exponents2 = GetExponents(atom2);
    matrix<size_t> QN2 = QuantumNumbers(atom2);
    TwoCenterOv_dtdp(overlaps,nao1,nao2,dist,exponents1,QN1,exponents2,QN2,orientation,rotate);
  }
  //Gaussian related stuff
  std::vector<double> ExpansionExponents(int Ngauss, int nquant, int lquant) {
    //exponents for expanding STO with principal quantum number nquant and azimuthal quantum number lquant into Ngauss Gaussian functions
    //R. F. Stewart, J. Chem. Phys., 52, 431, 1970
    std::vector<double> alpha(Ngauss,0.0);
    if (Ngauss == 1) {                   // Exponents from first row Table I-V
      if (nquant == 1) {alpha[0] = 2.709498091e-1;}
      else if (nquant == 2) {
        if (lquant == 0) {alpha[0] = 1.012151084e-1;}
        else if (lquant == 1) {alpha[0] = 1.759666885e-1;}
      }
      else if (nquant == 3) {
        if (lquant == 0) {alpha[0] = 5.296881757e-2;}
        else if (lquant == 1) {alpha[0] = 9.113614253e-2;}
        else if (lquant == 2) {alpha[0] = 1.302270363e-1;}
      }
      else if (nquant == 4) {
        if (lquant == 0) {alpha[0] = 3.264600274e-2;}
        else if (lquant == 1) {alpha[0] = 5.578350235e-2;}
        else if (lquant == 2) {alpha[0] = 7.941656339e-2;}
        else if (lquant == 3) {alpha[0] = 1.033434062e-1;}
      }
      else if (nquant == 5) {
        if (lquant == 0) {alpha[0] = 2.216912938e-2;}
        else if (lquant == 1) {alpha[0] = 3.769845216e-2;}
        else if (lquant == 2) {alpha[0] = 5.352200793e-2;}
        else if (lquant == 3) {alpha[0] = 6.952785407e-2;}
        else if (lquant == 4) {alpha[0] = 8.565417784e-2;}
      }
    }
    else if (Ngauss == 2) {              // Exponents from first row Table I-V
      if (nquant == 1) {
        alpha[0] = 8.518186635e-1;
        alpha[1] = 1.516232927e-1;
      }
      else if (nquant == 2) {
        if (lquant == 0) {
          alpha[0] = 1.292278611e-1;
          alpha[1] = 4.908584205e-2;
        }
        else if (lquant == 1) {
          alpha[0] = 4.323908358e-1;
          alpha[1] = 1.069139065e-1;
        }
      }
      else if (nquant == 3) {
        if (lquant == 0) {
          alpha[0] = 6.694095822e-1;
          alpha[1] = 5.837135094e-2;
        }
        else if (lquant == 1) {
          alpha[0] = 1.458620964e-1;
          alpha[1] = 5.664210742e-2;
        }
        else if (lquant == 2) {
          alpha[0] = 2.777427345e-1;
          alpha[1] = 8.336507714e-2;
        }
      }
      else if (nquant == 4) {
        if (lquant == 0) {
          alpha[0] = 2.441785453e-1;
          alpha[1] = 4.051097664e-2;
        }
        else if (lquant == 1) {
          alpha[0] = 6.190052680e-2;
          alpha[1] = 2.648418407e-2;
        }
        else if (lquant == 2) {
          alpha[0] = 1.330958892e-1;
          alpha[1] = 5.272119659e-2;
        }
        else if (lquant == 3) {
          alpha[0] = 2.006693538e-1;
          alpha[1] = 6.865384900e-2;
        }
      }
      else if (nquant == 5) {
        if (lquant == 0) {
          alpha[0] = 1.213425654e-1;
          alpha[1] = 3.133152144e-2;
        }
        else if (lquant == 1) {
          alpha[0] = 2.691294191e-1;
          alpha[1] = 3.980805011e-2;
        }
        else if (lquant == 2) {
          alpha[0] = 6.906014388e-2;
          alpha[1] = 3.399457777e-2;
        }
        else if (lquant == 3) {
          alpha[0] = 1.156094555e-1;
          alpha[1] = 4.778940916e-2;
        }
        else if (lquant == 4) {
          alpha[0] = 1.554531559e-1;
          alpha[1] = 5.854079811e-2;
        }
      }
    }
    else if (Ngauss == 3) {
      if (nquant == 1) {
        alpha[0] = 2.227660584;
        alpha[1] = 4.057711562e-1;
        alpha[2] = 1.098175104e-1;
      }
      else if (nquant == 2) {
        if (lquant == 0) {
          alpha[0] = 2.581578398;
          alpha[1] = 1.567622104e-1;
          alpha[2] = 6.018332272e-2;
        }
        else if (lquant == 1) {
          alpha[0] = 9.192379002e-1;
          alpha[1] = 2.359194503e-1;
          alpha[2] = 8.009805746e-2;
        }
      }
      else if (nquant == 3) {
        if (lquant == 0) {
          alpha[0] = 5.641487709e-1;
          alpha[1] = 6.924421391e-2;
          alpha[2] = 3.269529097e-2;
        }
        else if (lquant == 1) {
          alpha[0] = 2.692880368;
          alpha[1] = 1.489359592e-1;
          alpha[2] = 5.739585040e-2;
        }
        else if (lquant == 2) {
          alpha[0] = 5.229112225e-1;
          alpha[1] = 1.639595876e-1;
          alpha[2] = 6.386630021e-2;
        }
      }
      else if (nquant == 4) {
        if (lquant == 0) {
          alpha[0] = 2.267938753e-1;
          alpha[1] = 4.448178019e-2;
          alpha[2] = 2.195294664e-2;
        }
        else if (lquant == 1) {
          alpha[0] = 4.859692220e-1;
          alpha[1] = 7.430216918e-2;
          alpha[2] = 3.653340923e-2;
        }
        else if (lquant == 2) {
          alpha[0] = 1.777717219e-1;
          alpha[1] = 8.040647350e-2;
          alpha[2] = 3.949855551e-2;
        }
        else if (lquant == 3) {
          alpha[0] = 3.483826963e-1;
          alpha[1] = 1.249380537e-1;
          alpha[2] = 5.349995725e-2;
        }
      }
      else if (nquant == 5) {
        if (lquant == 0) {
          alpha[0] = 1.080198458e-1;
          alpha[1] = 4.408119382e-2;
          alpha[2] = 2.610811810e-2;
        }
        else if (lquant == 1) {
          alpha[0] = 2.127482317e-1;
          alpha[1] = 4.729648620e-2;
          alpha[2] = 2.604865324e-2;
        }
        else if (lquant == 2) {
          alpha[0] = 4.913352950e-1;
          alpha[1] = 7.329090601e-2;
          alpha[2] = 3.594209290e-2;
        }
        else if (lquant == 3) {
          alpha[0] = 1.649233885e-1;
          alpha[1] = 7.487066646e-2;
          alpha[2] = 3.735787219e-2;
        }
        else if (lquant == 4) {
          alpha[0] = 2.545432122e-1;
          alpha[1] = 1.006544376e-1;
          alpha[2] = 4.624463922e-2;
        }
      }
    }
    else if (Ngauss == 4) {
      if (nquant == 1) {
        alpha[0] = 5.216844534;
        alpha[1] = 9.546182760e-1;
        alpha[2] = 2.652034102e-1;
        alpha[3] = 8.801862774e-2;
      }
      else if (nquant == 2) {
        if (lquant == 0) {
          alpha[0] = 11.61525551;
          alpha[1] = 2.000243111;
          alpha[2] = 1.607280687e-1;
          alpha[3] = 6.125744532e-2;
        }
        else if (lquant == 1) {
          alpha[0] = 1.798260992;
          alpha[1] = 4.662622228e-1;
          alpha[2] = 1.643718620e-1;
          alpha[3] = 6.543927065e-2;
        }
      }
      else if (nquant == 3) {
        if (lquant == 0) {
          alpha[0] = 1.513265591;
          alpha[1] = 4.262497508e-1;
          alpha[2] = 7.643320863e-2;
          alpha[3] = 3.760545063e-2;
        }
        else if (lquant == 1) {
          alpha[0] = 1.853180239;
          alpha[1] = 1.915075719e-1;
          alpha[2] = 8.655487938e-2;
          alpha[3] = 4.184253862e-2;
        }
        else if (lquant == 2) {
          alpha[0] = 9.185846715e-1;
          alpha[1] = 2.920461109e-1;
          alpha[2] = 1.187568890e-1;
          alpha[3] = 5.286755896e-2;
        }
      }
      else if (nquant == 4) {
        if (lquant == 0) {
          alpha[0] = 3.242212833e-1;
          alpha[1] = 1.663217177e-1;
          alpha[2] = 5.081097451e-2;
          alpha[3] = 2.829066600e-2;
        }
        else if (lquant == 1) {
          alpha[0] = 1.492607880;
          alpha[1] = 4.327619272e-1;
          alpha[2] = 7.553156064e-2;
          alpha[3] = 3.706272183e-2;
        }
        else if (lquant == 2) {
          alpha[0] = 1.995825422;
          alpha[1] = 1.823461280e-1;
          alpha[2] = 8.197240896e-2;
          alpha[3] = 4.000634951e-2;
        }
        else if (lquant == 3) {
          alpha[0] = 5.691670217e-1;
          alpha[1] = 2.074585819e-1;
          alpha[2] = 9.298346885e-2;
          alpha[3] = 4.473508853e-2;
        }
      }
      else if (nquant == 5) {
        if (lquant == 0) {
          alpha[0] = 8.602284252e-1;
          alpha[1] = 1.189050200e-1;
          alpha[2] = 3.446076176e-2;
          alpha[3] = 1.974798796e-2;
        }
        else if (lquant == 1) {
          alpha[0] = 3.962838833e-1;
          alpha[1] = 1.838858552e-1;
          alpha[2] = 4.943555157e-2;
          alpha[3] = 2.750222273e-2;
        }
        else if (lquant == 2) {
          alpha[0] = 4.230617826e-1;
          alpha[1] = 8.293863702e-2;
          alpha[2] = 4.590326388e-2;
          alpha[3] = 2.628744797e-2;
        }
        else if (lquant == 3) {
          alpha[0] = 2.017831152e-1;
          alpha[1] = 1.001952178e-1;
          alpha[2] = 5.447006630e-2;
          alpha[3] = 3.037569283e-2;
        }
        else if (lquant == 4) {
          alpha[0] = 3.945205573e-1;
          alpha[1] = 1.588100623e-1;
          alpha[2] = 7.646521729e-1;
          alpha[3] = 3.898703611e-2;
        }
      }
    }
    else if (Ngauss == 5) {
      if (nquant == 1) {
        alpha[0] = 11.30563696;
        alpha[1] = 2.071728178;
        alpha[2] = 5.786484833e-1;
        alpha[3] = 1.975724573e-1;
        alpha[4] = 7.445271746e-2;
      }
      else if (nquant == 2) {
        if (lquant == 0) {
          alpha[0] = 8.984956862;
          alpha[1] = 1.673710636;
          alpha[2] = 1.944726668e-1;
          alpha[3] = 8.806345634e-2;
          alpha[4] = 4.249068522e-2;
        }
        else if (lquant == 1) {
          alpha[0] = 3.320386533;
          alpha[1] = 8.643257633e-1;
          alpha[2] = 3.079819284e-1;
          alpha[3] = 1.273309895e-1;
          alpha[4] = 5.606243164e-2;
        }
      }
      else if (nquant == 3) {
        if (lquant == 0) {
          alpha[0] = 4.275877914;
          alpha[1] = 1.132409433;
          alpha[2] = 4.016256968e-1;
          alpha[3] = 7.732370620e-2;
          alpha[4] = 3.800708627e-2;
        }
        else if (lquant == 1) {
          alpha[0] = 6.466803859;
          alpha[1] = 1.555914802;
          alpha[2] = 1.955925255e-1;
          alpha[3] = 8.809647701e-2;
          alpha[4] = 4.234835707e-2;
        }
        else if (lquant == 2) {
          alpha[0] = 1.539033958;
          alpha[1] = 4.922090297e-1;
          alpha[2] = 2.029756928e-1;
          alpha[3] = 9.424112917e-2;
          alpha[4] = 4.569058269e-2;
        }
      }
      else if (nquant == 4) {
        if (lquant == 0) {
          alpha[0] = 2.980263783;
          alpha[1] = 3.792228833e-1;
          alpha[2] = 1.789717224e-1;
          alpha[3] = 5.002110360e-2;
          alpha[4] = 2.789361681e-2;
        }
        else if (lquant == 1) {
          alpha[0] = 1.091977298;
          alpha[1] = 3.719985051e-1;
          alpha[2] = 8.590019352e-2;
          alpha[3] = 4.786503860e-2;
          alpha[4] = 2.730479990e-2;
        }
        else if (lquant == 2) {
          alpha[0] = 1.522122079;
          alpha[1] = 2.173041823e-1;
          alpha[2] = 1.084876577e-1;
          alpha[3] = 5.836797641e-2;
          alpha[4] = 3.206682246e-2;
        }
        else if (lquant == 3) {
          alpha[0] = 8.925960415e-1;
          alpha[1] = 3.277589120e-1;
          alpha[2] = 1.492869962e-1;
          alpha[3] = 7.506099109e-2;
          alpha[4] = 3.892475795e-2;
        }
      }
      else if (nquant == 5) {
        if (lquant == 0) {
          alpha[0] = 7.403763257e-1;
          alpha[1] = 1.367990863e-1;
          alpha[2] = 9.135301779e-2;
          alpha[3] = 3.726907315e-2;
          alpha[4] = 2.241490836e-2;
        }
        else if (lquant == 1) {
          alpha[0] = 3.422168934e-1;
          alpha[1] = 1.665099900e-1;
          alpha[2] = 5.443732013e-2;
          alpha[3] = 3.367775277e-2;
          alpha[4] = 2.091949042e-2;
        }
        else if (lquant == 2) {
          alpha[0] = 9.702946470e-1;
          alpha[1] = 3.603270196e-1;
          alpha[2] = 8.668717752e-2;
          alpha[3] = 4.833708379e-2;
          alpha[4] = 2.751899341e-2;
        }
        else if (lquant == 3) {
          alpha[0] = 1.670735676;
          alpha[1] = 2.072477219e-1;
          alpha[2] = 1.024709357e-1;
          alpha[3] = 5.537913898e-2;
          alpha[4] = 3.072866652e-2;
        }
        else if (lquant == 4) {
          alpha[0] = 5.895429375e-1;
          alpha[1] = 2.393343780e-1;
          alpha[2] = 1.172646904e-1;
          alpha[3] = 6.254074479e-2;
          alpha[4] = 3.411243214e-2;
        }
      }
    }
    else if (Ngauss == 6) {
      if (nquant == 1) {
        alpha[0] = 23.10303149;
        alpha[1] = 4.235915534;
        alpha[2] = 1.185056519;
        alpha[3] = 4.070988982e-1;
        alpha[4] = 1.580884151e-1;
        alpha[5] = 6.510953954e-2;
      }
      else if (nquant == 2) {
        if (lquant == 0) {
          alpha[0] = 27.68496241;
          alpha[1] = 5.077140627;
          alpha[2] = 1.426786050;
          alpha[3] = 2.040335729e-1;
          alpha[4] = 9.260298399e-2;
          alpha[5] = 4.416183978e-2;
        }
        else if (lquant == 1) {
          alpha[0] = 5.868285913;
          alpha[1] = 1.530329631;
          alpha[2] = 5.475665231e-1;
          alpha[3] = 2.288932733e-1;
          alpha[4] = 1.046655969e-1;
          alpha[5] = 4.948220127e-2;
        }
      }
      else if (nquant == 3) {
        if (lquant == 0) {
          alpha[0] = 3.273031938;
          alpha[1] = 9.200611311e-1;
          alpha[2] = 3.593349765e-1;
          alpha[3] = 8.636686991e-2;
          alpha[4] = 4.797373812e-2;
          alpha[5] = 2.724741144e-2;
        }
        else if (lquant == 1) {
          alpha[0] = 5.077973607;
          alpha[1] = 1.340786940;
          alpha[2] = 2.248434849e-1;
          alpha[3] = 1.131741848e-1;
          alpha[4] = 6.076408893e-2;
          alpha[5] = 3.315424265e-2;
        }
        else if (lquant == 2) {
          alpha[0] = 2.488296923;
          alpha[1] = 7.981487853e-1;
          alpha[2] = 3.311327490e-1;
          alpha[3] = 1.559114463e-1;
          alpha[4] = 7.877734732e-2;
          alpha[5] = 4.058484363e-2;
        }
      }
      else if (nquant == 4) {
        if (lquant == 0) {
          alpha[0] = 1.365346;
          alpha[1] = 4.393213e-01;
          alpha[2] = 1.877069e-01;
          alpha[3] = 9.360270e-02;
          alpha[4] = 5.052263e-02;
          alpha[5] = 2.809354e-02;
          //alpha[0] = 3.232838646;
          //alpha[1] = 3.605788802e-1;
          //alpha[2] = 1.717905487e-1;
          //alpha[3] = 5.277666487e-2;
          //alpha[4] = 3.163400284e-2;
          //alpha[5] = 1.874093091e-2;
        }
        else if (lquant == 1) {
          alpha[0] = 1.365346;
          alpha[1] = 4.393213e-01;
          alpha[2] = 1.877069e-01;
          alpha[3] = 9.360270e-02;
          alpha[4] = 5.052263e-02;
          alpha[5] = 2.809354e-02;
          //alpha[0] = 2.389722618;
          //alpha[1] = 7.960947826e-1;
          //alpha[2] = 3.415541380e-1;
          //alpha[3] = 8.847434525e-2;
          //alpha[4] = 4.958248334e-2;
          //alpha[5] = 2.816929784e-2;
        }
        else if (lquant == 2) {
          alpha[0] = 4.634239420;
          alpha[1] = 1.341648295;
          alpha[2] = 2.209593028e-1;
          alpha[3] = 1.101467943e-1;
          alpha[4] = 5.904190370e-2;
          alpha[5] = 3.232628887e-2;
        }
        else if (lquant == 3) {
          alpha[0] = 1.357718039;
          alpha[1] = 5.004907278e-1;
          alpha[2] = 2.296565064e-1;
          alpha[3] = 1.173146814e-1;
          alpha[4] = 6.350097171e-2;
          alpha[5] = 3.474556673e-2;
        }
      }
      else if (nquant == 5) {
        if (lquant == 0) {
          alpha[0] = 1.410128298;
          alpha[1] = 5.077878915e-1;
          alpha[2] = 1.847926858e-1;
          alpha[3] = 1.061070594e-1;
          alpha[4] = 3.669584901e-2;
          alpha[5] = 2.213558430e-2;
        }
        else if (lquant == 1) {
          alpha[0] = 3.778623374;
          alpha[1] = 3.499121109e-1;
          alpha[2] = 1.683175469e-1;
          alpha[3] = 5.404070736e-2;
          alpha[4] = 3.328911801e-2;
          alpha[5] = 2.063815019e-2;
        }
        else if (lquant == 2) {
          alpha[0] = 8.820520428e-1;
          alpha[1] = 3.410838409e-1;
          alpha[2] = 9.204308840e-2;
          alpha[3] = 5.472831774e-2;
          alpha[4] = 3.391202830e-2;
          alpha[5] = 2.108227374e-2;
        }
        else if (lquant == 3) {
          alpha[0] = 1.334096840;
          alpha[1] = 2.372312347e-1;
          alpha[2] = 1.269485144e-1;
          alpha[3] = 7.290318381e-2;
          alpha[4] = 4.351355997e-2;
          alpha[5] = 2.598071843e-2;
        }
        else if (lquant == 4) {
          alpha[0] = 8.574668996e-1;
          alpha[1] = 3.497184772e-1;
          alpha[2] = 1.727917060e-1;
          alpha[3] = 9.373643151e-2;
          alpha[4] = 5.340032759e-2;
          alpha[5] = 3.057364464e-2;
        }
      }
      else if (nquant == 6) {
        if (lquant == 0) {
          alpha[0] = 5.800292686e-1;
          alpha[1] = 2.718262251e-1;
          alpha[2] = 7.938523262e-2;
          alpha[3] = 4.975088254e-2;
          alpha[4] = 2.983643556e-2;
          alpha[5] = 1.886067216e-2;
        }
        else if (lquant == 1) {
          alpha[0] = 6.696537714e-1;
          alpha[1] = 1.395089793e-1;
          alpha[2] = 8.163894960e-2;
          alpha[3] = 4.586329272e-2;
          alpha[4] = 2.961305556e-2;
          alpha[5] = 1.882221321e-2;
        }
      }
    }
    return alpha;
  }
  std::vector<double> ExpansionCoefficients(int Ngauss, int nquant, int lquant) {
    //exponents for expanding STO with principal quantum number nquant and azimuthal quantum number lquant into Ngauss Gaussian functions
    //R. F. Stewart, J. Chem. Phys., 52, 431, 1970
    std::vector<double> coefficient(Ngauss,0.0);
    if (Ngauss == 1) {coefficient[0] = 1.0;}
    else if (Ngauss == 2) {
      if (nquant == 1) {
        coefficient[0] = 4.301284983e-1;
        coefficient[1] = 6.789135305e-1;
      }
      else if (nquant == 2) {
        if (lquant == 0) {
          coefficient[0] = 7.470867124e-1;
          coefficient[1] = 2.855980556e-1;
        }
        else if (lquant == 1) {
          coefficient[0] = 4.522627513e-1;
          coefficient[1] = 6.713122642e-1;
        }
      }
      else if (nquant == 3) {
        if (lquant == 0) {
          coefficient[0] = -1.529645716e-1;
          coefficient[1] = 1.051370110;
        }
        else if (lquant == 1) {
          coefficient[0] = 5.349653114e-1;
          coefficient[1] = 5.299607212e-1;
        }
        else if (lquant == 2) {
          coefficient[0] = 4.666137923e-1;
          coefficient[1] = 6.644706516e-1;
        }
      }
      else if (nquant == 4) {
        if (lquant == 0) {
          coefficient[0] = -3.046656896e-1;
          coefficient[1] = 1.146877294;
        }
        else if (lquant == 1) {
          coefficient[0] = 8.743116767e-1;
          coefficient[1] = 1.513640107e-1;
        }
        else if (lquant == 2) {
          coefficient[0] = 4.932764167e-1;
          coefficient[1] = 5.918727866e-1;
        }
        else if (lquant == 3) {
          coefficient[0] = 4.769346276e-1;
          coefficient[1] = 6.587383976e-1;
        }
      }
      else if (nquant == 5) {
        if (lquant == 0) {
          coefficient[0] = -5.114756049e-1;
          coefficient[1] = 1.307377277;
        }
        else if (lquant == 1) {
          coefficient[0] = -1.034227010e-1;
          coefficient[1] = 1.033376378;
        }
        else if (lquant == 2) {
          coefficient[0] = 6.539405185e-1;
          coefficient[1] = 3.948945302e-1;
        }
        else if (lquant == 3) {
          coefficient[0] = 4.856637346e-1;
          coefficient[1] = 6.125980914e-1;
        }
        else if (lquant == 4) {
          coefficient[0] = 4.848298074e-1;
          coefficient[1] = 6.539381621e-1;
        }
      }
    }
    else if (Ngauss == 3) {
      if (nquant == 1) {
        coefficient[0] = 1.543289673e-1;
        coefficient[1] = 5.353281423e-1;
        coefficient[2] = 4.446345422e-1;
      }
      else if (nquant == 2) {
        if (lquant == 0) {
          coefficient[0] = -5.994474934e-2;
          coefficient[1] = 5.960385398e-1;
          coefficient[2] = 4.581786291e-1;
        }
        else if (lquant == 1) {
          coefficient[0] = 1.623948553e-1;
          coefficient[1] = 5.661708862e-1;
          coefficient[2] = 4.223071752e-1;
        }
      }
      else if (nquant == 3) {
        if (lquant == 0) {
          coefficient[0] = -1.782577972e-1;
          coefficient[1] = 8.612761663e-1;
          coefficient[2] = 2.261841969e-1;
        }
        else if (lquant == 1) {
          coefficient[0] = -1.061945788e-2;
          coefficient[1] = 5.218564264e-1;
          coefficient[2] = 5.450015143e-1;
        }
        else if (lquant == 2) {
          coefficient[0] = 1.686596060e-1;
          coefficient[1] = 5.847984817e-1;
          coefficient[2] = 4.056779523e-1;
        }
      }
      else if (nquant == 4) {
        if (lquant == 0) {
          coefficient[0] = -3.349048323e-1;
          coefficient[1] = 1.056744667;
          coefficient[2] = 1.256661680e-1;
        }
        else if (lquant == 1) {
          coefficient[0] = -6.147823411e-2;
          coefficient[1] = 6.604172234e-1;
          coefficient[2] = 3.932639495e-1;
        }
        else if (lquant == 2) {
          coefficient[0] = 2.308552718e-1;
          coefficient[1] = 6.042409177e-1;
          coefficient[2] = 2.595768926e-1;
        }
        else if (lquant == 3) {
          coefficient[0] = 1.737856685e-1;
          coefficient[1] = 5.973380628e-1;
          coefficient[2] = 3.929395614e-1;
        }
      }
      else if (nquant == 5) {
        if (lquant == 0) {
          coefficient[0] = -6.617401158e-1;
          coefficient[1] = 7.467595004e-1;
          coefficient[2] = 7.146490945e-1;
        }
        else if (lquant == 1) {
          coefficient[0] = -1.389529695e-1;
          coefficient[1] = 8.076691064e-1;
          coefficient[2] = 2.726029342e-1;
        }
        else if (lquant == 2) {
          coefficient[0] = -2.010175008e-2;
          coefficient[1] = 5.899310608e-1;
          coefficient[2] = 4.658445960e-1;
        }
        else if (lquant == 3) {
          coefficient[0] = 1.909729355e-1;
          coefficient[1] = 6.146060459e-1;
          coefficient[2] = 3.059611271e-1;
        }
        else if (lquant == 4) {
          coefficient[0] = 1.780980905e-1;
          coefficient[1] = 6.063757846e-1;
          coefficient[2] = 3.828552923e-1;
        }
      }
    }
    else if (Ngauss == 4) {
      if (nquant == 1) {
        coefficient[0] = 5.675242080e-2;
        coefficient[1] = 2.601413550e-1;
        coefficient[2] = 5.328461143e-1;
        coefficient[3] = 2.916254405e-1;
      }
      else if (nquant == 2) {
        if (lquant == 0) {
          coefficient[0] = -1.198411747e-2;
          coefficient[1] = -5.472052539e-2;
          coefficient[2] = 5.805587176e-1;
          coefficient[3] = 4.770079976e-1;
        }
        else if (lquant == 1) {
          coefficient[0] = 5.713170255e-2;
          coefficient[1] = 2.857455515e-1;
          coefficient[2] = 5.517873105e-1;
          coefficient[3] = 2.632314924e-1;
        }
      }
      else if (nquant == 3) {
        if (lquant == 0) {
          coefficient[0] = -3.295496352e-2;
          coefficient[1] = -1.724516959e-1;
          coefficient[2] = 7.518511194e-1;
          coefficient[3] = 3.589627317e-1;
        }
        else if (lquant == 1) {
          coefficient[0] = -1.434249391e-2;
          coefficient[1] = 2.755177589e-1;
          coefficient[2] = 5.846750879e-1;
          coefficient[3] = 2.144986514e-1;
        }
        else if (lquant == 2) {
          coefficient[0] = 5.799057705e-2;
          coefficient[1] = 3.045581349e-1;
          coefficient[2] = 5.601358038e-1;
          coefficient[3] = 2.432423313e-1;
        }
      }
      else if (nquant == 4) {
        if (lquant == 0) {
          coefficient[0] = -1.120682822e-1;
          coefficient[1] = -2.845426863e-1;
          coefficient[2] = 8.909873788e-1;
          coefficient[3] = 3.517811205e-1;
        }
        else if (lquant == 1) {
          coefficient[0] = -6.035216774e-3;
          coefficient[1] = -6.013310874e-2;
          coefficient[2] = 6.451518200e-1;
          coefficient[3] = 4.117923820e-1;
        }
        else if (lquant == 2) {
          coefficient[0] = -2.816702620e-3;
          coefficient[1] = 2.177095871e-1;
          coefficient[2] = 6.058047348e-1;
          coefficient[3] = 2.717811257e-1;
        }
        else if (lquant == 3) {
          coefficient[0] = 5.902730589e-2;
          coefficient[1] = 3.191828952e-1;
          coefficient[2] = 5.639423893e-1;
          coefficient[3] = 2.284796537e-1;
        }
      }
      else if (nquant == 5) {
        if (lquant == 0) {
          coefficient[0] = 1.103657561e-2;
          coefficient[1] = -5.606519023e-1;
          coefficient[2] = 1.179429987;
          coefficient[3] = 1.734974376e-1;
        }
        else if (lquant == 1) {
          coefficient[0] = -1.801459207e-2;
          coefficient[1] = -1.360777372e-1;
          coefficient[2] = 7.533973719e-1;
          coefficient[3] = 3.409304859e-1;
        }
        else if (lquant == 2) {
          coefficient[0] = -2.421626009e-2;
          coefficient[1] = 3.937644956e-1;
          coefficient[2] = 5.489520286e-1;
          coefficient[3] = 1.190436963e-1;
        }
        else if (lquant == 3) {
          coefficient[0] = 9.174268830e-2;
          coefficient[1] = 4.023496947e-1;
          coefficient[2] = 4.937432100e-1;
          coefficient[3] = 1.254001522e-1;
        }
        else if (lquant == 4) {
          coefficient[0] = 6.010484250e-2;
          coefficient[1] = 3.309738329e-1;
          coefficient[2] = 5.655207585e-1;
          coefficient[3] = 2.171122608e-1;
        }
      }
    }
    else if (Ngauss == 5) {
      if (nquant == 1) {
        coefficient[0] = 2.214055312e-2;
        coefficient[1] = 1.135411520e-1;
        coefficient[2] = 3.318161484e-1;
        coefficient[3] = 4.825700713e-1;
        coefficient[4] = 1.935721966e-1;
      }
      else if (nquant == 2) {
        if (lquant == 0) {
          coefficient[0] = -1.596349096e-2;
          coefficient[1] = -5.685884883e-2;
          coefficient[2] = 3.698265599e-1;
          coefficient[3] = 5.480512593e-1;
          coefficient[4] = 1.472634893e-1;
        }
        else if (lquant == 1) {
          coefficient[0] = 2.079051117e-2;
          coefficient[1] = 1.235472099e-1;
          coefficient[2] = 3.667738986e-1;
          coefficient[3] = 4.834930290e-1;
          coefficient[4] = 1.653444074e-1;
        }
      }
      else if (nquant == 3) {
        if (lquant == 0) {
          coefficient[0] = -3.920358850e-3;
          coefficient[1] = -4.168430506e-2;
          coefficient[2] = -1.637440990e-1;
          coefficient[3] = 7.419373723e-1;
          coefficient[4] = 3.724364929e-1;
        }
        else if (lquant == 1) {
          coefficient[0] = -2.329023747e-3;
          coefficient[1] = -1.357395221e-2;
          coefficient[2] = 2.632185383e-1;
          coefficient[3] = 5.880427024e-1;
          coefficient[4] = 2.242794445e-1;
        }
        else if (lquant == 2) {
          coefficient[0] = 2.020869128e-2;
          coefficient[1] = 1.321157923e-1;
          coefficient[2] = 3.911240346e-1;
          coefficient[3] = 4.779609701e-1;
          coefficient[4] = 1.463662294e-1;
        }
      }
      else if (nquant == 4) {
        if (lquant == 0) {
          coefficient[0] = 1.513948997e-3;
          coefficient[1] = -7.316801518e-2;
          coefficient[2] = -3.143703799e-1;
          coefficient[3] = 9.032615169e-1;
          coefficient[4] = 3.294210848e-1;
        }
        else if (lquant == 1) {
          coefficient[0] = -1.143929558e-2;
          coefficient[1] = -6.322651538e-2;
          coefficient[2] = 4.398907721e-1;
          coefficient[3] = 5.245859166e-1;
          coefficient[4] = 1.017072253e-1;
        }
        else if (lquant == 2) {
          coefficient[0] = -3.673711876e-3;
          coefficient[1] = 1.167122499e-1;
          coefficient[2] = 4.216476416e-1;
          coefficient[3] = 4.547673415e-1;
          coefficient[4] = 1.037803318e-1;
        }
        else if (lquant == 3) {
          coefficient[0] = 1.999839052e-2;
          coefficient[1] = 1.395427440e-1;
          coefficient[2] = 4.091508237e-1;
          coefficient[3] = 4.708252119e-1;
          coefficient[4] = 1.328082566e-1;
        }
      }
      else if (nquant == 5) {
        if (lquant == 0) {
          coefficient[0] = 1.375523371e-2;
          coefficient[1] = -3.097344179e-1;
          coefficient[2] = -3.199192259e-1;
          coefficient[3] = 1.084547038;
          coefficient[4] = 3.345288361e-1;
        }
        else if (lquant == 1) {
          coefficient[0] = -3.113958289e-2;
          coefficient[1] = -1.374007017e-1;
          coefficient[2] = 5.573881018e-1;
          coefficient[3] = 4.855428100e-1;
          coefficient[4] = 6.605423564e-2;
        }
        else if (lquant == 2) {
          coefficient[0] = -3.231527611e-3;
          coefficient[1] = -2.434931372e-2;
          coefficient[2] = 3.440817054e-1;
          coefficient[3] = 5.693674316e-1;
          coefficient[4] = 1.511340183e-1;
        }
        else if (lquant == 3) {
          coefficient[0] = -7.301193568e-4;
          coefficient[1] = 8.414991343e-2;
          coefficient[2] = 3.923683153e-1;
          coefficient[3] = 5.040033146e-1;
          coefficient[4] = 1.328979300e-1;
        }
        else if (lquant == 4) {
          coefficient[0] = 1.998085812e-2;
          coefficient[1] = 1.460384050e-1;
          coefficient[2] = 4.230565459e-1;
          coefficient[3] = 4.635699665e-1;
          coefficient[4] = 1.226411691e-1;
        }
      }
    }
    else if (Ngauss == 6) {
      if (nquant == 1) {
        coefficient[0] = 9.163596280e-3;
        coefficient[1] = 4.936149294e-2;
        coefficient[2] = 1.685383049e-1;
        coefficient[3] = 3.705627997e-1;
        coefficient[4] = 4.164915298e-1;
        coefficient[5] = 1.303340841e-1;
      }
      else if (nquant == 2) {
        if (lquant == 0) {
          coefficient[0] = -4.151277819e-3;
          coefficient[1] = -2.067024148e-2;
          coefficient[2] = -5.150303337e-2;
          coefficient[3] = 3.346271174e-1;
          coefficient[4] = 5.621061301e-1;
          coefficient[5] = 1.712994697e-1;
        }
        else if (lquant == 1) {
          coefficient[0] = 7.924233646e-3;
          coefficient[1] = 5.144104825e-2;
          coefficient[2] = 1.898400060e-1;
          coefficient[3] = 4.049863191e-1;
          coefficient[4] = 4.012362861e-1;
          coefficient[5] = 1.051855189e-1;
        }
      }
      else if (nquant == 3) {
        if (lquant == 0) {
          coefficient[0] = -6.775596947e-3;
          coefficient[1] = -5.639325779e-2;
          coefficient[2] = -1.587856086e-1;
          coefficient[3] = 5.534527651e-1;
          coefficient[4] = 5.015351020e-1;
          coefficient[5] = 7.223633674e-2;
        }
        else if (lquant == 1) {
          coefficient[0] = -3.329929840e-3;
          coefficient[1] = -1.419488340e-2;
          coefficient[2] = 1.639395770e-1;
          coefficient[3] = 4.485358256e-1;
          coefficient[4] = 3.908813050e-1;
          coefficient[5] = 7.411456232e-2;
        }
        else if (lquant == 2) {
          coefficient[0] = 7.283828112e-3;
          coefficient[1] = 5.386799363e-2;
          coefficient[2] = 2.072139149e-1;
          coefficient[3] = 4.266269092e-1;
          coefficient[4] = 3.843100204e-1;
          coefficient[5] = 8.902827546e-2;
        }
      }
      else if (nquant == 4) {
        if (lquant == 0) {
          coefficient[0] = 3.775056e-03;
          coefficient[1] = -5.585965e-02;
          coefficient[2] = -3.192946e-01;
          coefficient[3] = -2.764780e-02;
          coefficient[4] = 9.049199e-01;
          coefficient[5] = 3.406258e-01;
          //coefficient[0] = 1.374817488e-3;
          //coefficient[1] = -8.666390043e-2;
          //coefficient[2] = -3.130627309e-1;
          //coefficient[3] = 7.812787397e-1;
          //coefficient[4] = 4.389247988e-1;
          //coefficient[5] = 2.487178756e-2;
        }
        else if (lquant == 1) {
          coefficient[0] = -7.052075e-03;
          coefficient[1] = -5.259505e-02;
          coefficient[2] = -3.773450e-02;
          coefficient[3] = 3.874773e-01;
          coefficient[4] = 5.791672e-01;
          coefficient[5] = 1.221817e-01;
          //coefficient[0] = -1.665913575e-3;
          //coefficient[1] = -1.657464971e-2;
          //coefficient[2] = -5.958513378e-2;
          //coefficient[3] = 4.053115554e-1;
          //coefficient[4] = 5.433958189e-1;
          //coefficient[5] = 1.204970491e-1;
        }
        else if (lquant == 2) {
          coefficient[0] = -4.749842876e-4;
          coefficient[1] = -3.566777891e-3;
          coefficient[2] = 1.108670481e-1;
          coefficient[3] = 4.159646930e-1;
          coefficient[4] = 4.621672517e-1;
          coefficient[5] = 1.081250196e-1;
        }
        else if (lquant == 3) {
          coefficient[0] = 6.930234381e-3;
          coefficient[1] = 5.634263745e-2;
          coefficient[2] = 2.217065797e-1;
          coefficient[3] = 4.411388883e-1;
          coefficient[4] = 3.688112625e-1;
          coefficient[5] = 7.787514504e-2;
        }
      }
      else if (nquant == 5) {
        if (lquant == 0) {
          coefficient[0] = 2.695439582e-3;
          coefficient[1] = 1.850157487e-2;
          coefficient[2] = -9.588628125e-2;
          coefficient[3] = -5.200673560e-1;
          coefficient[4] = 1.087619490;
          coefficient[5] = 3.103964343e-1;
        }
        else if (lquant == 1) {
          coefficient[0] = 1.163246387e-4;
          coefficient[1] = -2.920771322e-2;
          coefficient[2] = -1.381051233e-1;
          coefficient[3] = 5.706134877e-1;
          coefficient[4] = 4.768808140e-1;
          coefficient[5] = 6.021665516e-2;
        }
        else if (lquant == 2) {
          coefficient[0] = -4.097311019e-3;
          coefficient[1] = -2.508271857e-2;
          coefficient[2] = 2.648458555e-1;
          coefficient[3] = 5.097437054e-1;
          coefficient[4] = 2.654483467e-1;
          coefficient[5] = 2.623132212e-2;
        }
        else if (lquant == 3) {
          coefficient[0] = -9.486751531e-4;
          coefficient[1] = 4.624275998e-2;
          coefficient[2] = 2.373699784e-1;
          coefficient[3] = 4.589112231e-1;
          coefficient[4] = 3.205010548e-1;
          coefficient[5] = 5.077063693e-2;
        }
        else if (lquant == 4) {
          coefficient[0] = 6.729778096e-3;
          coefficient[1] = 5.874145170e-2;
          coefficient[2] = 2.339955227e-1;
          coefficient[3] = 4.512983737e-1;
          coefficient[4] = 3.552053926e-1;
          coefficient[5] = 6.974153145e-2;
        }
      }
      else if (nquant == 6) {
        if (lquant == 0) {
          coefficient[0] = 4.554359511e-3;
          coefficient[1] = 5.286443143e-2;
          coefficient[2] = -7.561016358e-1;
          coefficient[3] = -2.269803820e-1;
          coefficient[4] = 1.332494651;
          coefficient[5] = 3.622518293e-1;
        }
        else if (lquant == 1) {
          coefficient[0] = 2.782723680e-3;
          coefficient[1] = -1.282887780e-1;
          coefficient[2] = -2.266255943e-1;
          coefficient[3] = 4.682259383e-1;
          coefficient[4] = 6.752048848e-1;
          coefficient[5] = 1.091534212e-1;
        }
      }
    }
    return coefficient;
  }
  virtual int getNumberPrimitives(size_t atomnr, int lnumb) {
    //function determining the number of Gaussian primitives used for a given element, orbital with certain azimuthal quantum number
    return 0;
  }
  void slater2Gauss(int ngauss, int nprinc, int lnumb, double zeta, std::vector<double> & alpha, std::vector<double> & coeff, bool norm, double thresholdzero = 1.0e-8) {
    //function converting one STO into a linear combination of GTOs (with expansion coefficients as well)
    //ngauss -> number of Gaussian functions for the expansion
    //nprinc -> principal quantum number (n)
    //lnumb  -> azimuthal quantum number (l)
    //zeta   -> STO exponent
    //alpha  -> vector with GTO exponents
    //coeff  -> vector with GTO contraction coefficients
    //norm   -> boolean to determine whether to normalize GTO
    //return integer saying whether successful
    //don't allow negative exponents and check that quantum numbers are consistent
    if (zeta < thresholdzero) {throw("ERROR: gfn2.hpp: GFN2BSet: slater2Gauss(): negative STO exponent not allowed");}
    if ((nprinc > 5)||(nprinc < lnumb + 1)) {
       if (!((nprinc == 6)&&(ngauss == 6))) {throw("ERROR: gfn2.hpp: GFN2BSet: slater2Gauss(): inconsistent quantum numbers");}
    }
    alpha = ExpansionExponents(ngauss,nprinc,lnumb);
    coeff = ExpansionCoefficients(ngauss,nprinc,lnumb);
    for (size_t idx = 0; idx < ngauss; ++idx) {
      alpha[idx] *= zeta*zeta;
    }
    //normalize the gaussian if requested
    //<φ|φ> = (2i-1)!!(2j-1)!!(2k-1)!!/(4α)^(i+j+k) · sqrt(π/2α)³
    //N² = (4α)^(i+j+k)/((2i-1)!!(2j-1)!!(2k-1)!!)  · sqrt(2α/π)³
    //N = (4α)^((i+j+k)/2) / sqrt((2i-1)!!(2j-1)!!(2k-1)!!) · (2α/π)^(3/4)
    double aux = 2.0/pi;
    if (norm) {
      for (size_t idx = 0; idx < ngauss; ++idx) {
        coeff[idx] *= pow(aux*alpha[idx],0.75)*pow(sqrt(4.0*alpha[idx]),lnumb)/sqrt(DoubleFactorial2(std::abs(2*lnumb - 1)));
      }
    }
  }
  virtual int IPType() {return 0;}
};

//bunch of functions necessary for STOs
double NormFactor(int na, int la, int ma, double zeta) {
  //function calculating the normalization factor for a given STO
  double ama = std::abs(ma);
  double factorials = FactorialRatio(la + ama,la - ama);
  double norm = pow(2.0*zeta,na)*sqrt(zeta*double(2*la + 1)/(factorials*Factorial(2*na)));
  return norm;
}
double RNormFactor(int na, double expa, double & rval, double r) {
  //tabulated form of the norm for radial STO functions
  double norm = 2.0*sqrt(expa*expa*expa);         //N1
  rval = 1.0;
  if (na > 1) {           //N2
    norm *= expa/sqrt(3.0);
    rval *= r;
  }
  if (na > 2) {           //N3
    norm *= expa*sqrt(2.0/15.0);
    rval *= r;
  }
  if (na > 3) {           //N4
    norm *= expa/sqrt(14.0);
    rval *= r;
  }
  if (na > 4) {           //N5
    norm *= expa*sqrt(2.0/45.0);
    rval *= r;
  }
  if (na > 5) {           //N6
    norm *= expa/sqrt(33.0);
    rval *= r;
  }
  return norm;
}
double RSlater(double r, double a, int n) {
  //function evaluating the radial part of a Slater orbital
  double rval = 1.0;
  double norm = RNormFactor(n,a,rval,r);
  return norm*rval*exp(-a*r);
}
double ASlater(int l, int ao, double x, double y, double z, double r, double r2) {
  //function evaluating the angular part of a Slater orbital
  double aval = osqrtpi;
  if (l == 0) {aval *= 0.5;}
  else if (l == 1) {aval *= 0.5*sqrt3*((ao == 0)*x + (ao == 1)*y + (ao == 2)*z)/r;}
  else if (l == 2) {
    if (ao == 2) {aval *= 0.25*sqrt5*(2.0*z*z - x*x - y*y)/r2;}
    else if (ao == 4) {aval *= 0.25*sqrt15*(x*x - y*y)/r2;}
    else {aval *= 0.5*sqrt15*((ao == 0)*x*z + (ao == 1)*y*z + (ao == 3)*x*y)/r2;}
  }
  return aval;
}

#endif //_Slater_Type_Orbitals_
