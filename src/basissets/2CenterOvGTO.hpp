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

#ifndef _Two_Center_Overlaps_GTO_
#define _Two_Center_Overlaps_GTO_
#include <vector>
#include <string>
#include <math.h>
#include "OvRot.hpp"
#include "OvFunctions.hpp"
#include "OrbitalNumeration.hpp"
#include "../math/GaussianPackage.hpp"
#include "../math/MatrixPackage.hpp"
#include "../math/TensorPackage.hpp"
#include "../math/IntegralPackage.hpp"

//description:
//functions for the calculation of two-center overlap functions using GTOs

void BuildHShift(std::vector<double> & coeff, double aa, double ee, int lorb) {
  //calculation of Cartesian prefactor for gaussians
  //this is a C++ translation of build_hshift in xTB
  //coeff -> cartesian prefactor(s)
  //a -> component of center for a single gaussian
  //a -> component of center for product gaussian
  //s case, do nothing
  if (lorb > 0) {
    ee -= aa;
    if (lorb == 1) {coeff[0] = ee*coeff[1];}     //p
    else if (lorb == 2) {                          //d
      coeff[0] = ee*ee*coeff[2];
      coeff[1] = 2.0*ee*coeff[2];
    }
    else if (lorb == 3) {                          //f
      coeff[0] = ee*ee*ee*coeff[3];
      coeff[1] = 3.0*ee*ee*coeff[3];
      coeff[2] = 3.0*ee*coeff[3];
    }
    else if (lorb == 4) {                          //g
      coeff[0] = ee*ee*ee*ee*coeff[4];
      coeff[1] = 4.0*ee*ee*ee*coeff[4];
      coeff[2] = 6.0*ee*ee*coeff[4];
      coeff[3] = 4.0*ee*coeff[4];
    }
  }
}
void HorizontalShift(std::vector<double> & cfs, double shift, int lorb) {
  switch (lorb) {
    case 0:                 //s
      //no shift
      break;
    case 1:                 //p
      cfs[0] += shift*cfs[1];
      break;
    case 2:                 //d
      cfs[0] += shift*shift*cfs[2];
      cfs[1] += 2.0*shift*cfs[2];
      break;
    case 3:                 //f
      cfs[0] += shift*shift*shift*cfs[3];
      cfs[1] += 3.0*shift*shift*cfs[3];
      cfs[2] += 3.0*shift*cfs[3];
      break;
    case 4:                 //g
      cfs[0] += shift*shift*shift*shift*cfs[4];
      cfs[1] += 4.0*shift*shift*shift*cfs[4];
      cfs[2] += 6.0*shift*shift*cfs[4];
      cfs[3] += 4.0*shift*cfs[4];
      break;
  }
}
void Prod3(std::vector<double> & RA, std::vector<double> & RB, std::vector<double> & RC, int lorbA, int lorbB) {
  //get the product of Cartesian prefactors for gaussians
  //this is a C++ translation of prod3 in xTB
  int maxL = std::max(lorbA,lorbB);
  if (maxL == 4) {                 //g
    //<s|g> = <s|*(|s> + |p> + |d> + |f> + |g>) = <s> + <p> + <d> + <f> + <g>
    RC[0] = RA[0]*RB[0];
    RC[1] = RA[0]*RB[1] + RA[1]*RB[0];
    RC[2] = RA[0]*RB[2] + RA[2]*RB[0];
    RC[3] = RA[0]*RB[3] + RA[3]*RB[0];
    RC[4] = RA[0]*RB[4] + RA[4]*RB[0];
    if ((lorbA > 0)&&(lorbB > 0)) {
      //<p|g> = (<s| + <p|)*(|s> + |p> + |d> + |f> + |g>) = <s> + <p> + <d> + <f> + <g> + <h>
      RC[2] += RA[1]*RB[1];
      RC[3] += RA[1]*RB[2] + RA[2]*RB[1];
      RC[4] += RA[1]*RB[3] + RA[3]*RB[1];
      RC[5] = RA[1]*RB[4] + RA[4]*RB[1];
      if ((lorbA > 1)&&(lorbB > 1)) {
        //<d|g> = (<s| + <p| + <d|)*(|s> + |p> + |d> + |f> + |g>) = <s> + <p> + <d> + <f> + <g> + <h> + <i>
        RC[4] += RA[2]*RB[2];
        RC[5] += RA[2]*RB[3] + RA[3]*RB[2];
        RC[6] = RA[2]*RB[4] + RA[4]*RB[2];
        if ((lorbA > 2)&&(lorbB > 2)) {
          //<f|g> = (<s| + <p| + <d| + <f|)*(|s> + |p> + |d> + |f> + |g>) = <s> + <p> + <d> + <f> + <g> + <h> + <i> + <k>
          RC[6] += RA[3]*RB[3];
          RC[7] = RA[3]*RB[4] + RA[4]*RB[3];
          if ((lorbA == 4)&&(lorbB == 4)) {
            //<g|g> = (<s| + <p| + <d| + <f| + <g|)*(|s> + |p> + |d> + |f> + |g>) = <s> + <p> + <d> + <f> + <g> + <h> + <i> + <k> + <l>
            RC[8] = RA[4]*RB[4];
          }
        }
      }
    }
  }
  else if (maxL == 3) {            //f
    //<s|f> = <s|*(|s> + |p> + |d> + |f>) = <s> + <p> + <d> + <f>
    RC[0] = RA[0]*RB[0];
    RC[1] = RA[0]*RB[1] + RA[1]*RB[0];
    RC[2] = RA[0]*RB[2] + RA[2]*RB[0];
    RC[3] = RA[0]*RB[3] + RA[3]*RB[0];
    if ((lorbA > 0)&&(lorbB > 0)) {
      //<p|f> = (<s| + <p|)*(|s> + |p> + |d> + |f>) = <s> + <p> + <d> + <f> + <g>
      RC[2] += RA[1]*RB[1];
      RC[3] += RA[1]*RB[2] + RA[2]*RB[1];
      RC[4] = RA[1]*RB[3] + RA[3]*RB[1];
      if ((lorbA > 1)&&(lorbB > 1)) {
        //<d|f> = (<s| + <p| + <d|)*(|s> + |p> + |d> + |f>) = <s> + <p> + <d> + <f> + <g> + <h>
        RC[4] += RA[2]*RB[2];
        RC[5] = RA[2]*RB[3] + RA[3]*RB[2];
        if ((lorbA == 3)&&(lorbB == 3)) {
          //<f|f> = (<s| + <p| + <d| + <f|)*(|s> + |p> + |d> + |f>) = <s> + <p> + <d> + <f> + <g> + <h> + <i>
          RC[6] = RA[3]*RB[3];
        }
      }
    }
  }
  else if (maxL == 2) {            //d
    //<s|d> = <s|*(|s> + |p> + |d>) = <s> + <p> + <d>
    RC[0] = RA[0]*RB[0];
    RC[1] = RA[0]*RB[1] + RA[1]*RB[0];
    RC[2] = RA[0]*RB[2] + RA[2]*RB[0];
    if ((lorbA > 0)&&(lorbB > 0)) {
      //<p|d> = (<s| + <p|)*(|s> + |p> + |d>) = <s> + <p> + <d> + <f>
      RC[2] += RA[1]*RB[1];
      RC[3] = RA[1]*RB[2] + RA[2]*RB[1];
      if ((lorbA == 2)&&(lorbB == 2)) {
        //<d|d> = (<s| + <p| + <d|)*(|s> + |p> + |d>) = <s> + <p> + <d> + <f> + <g>
        RC[4] = RA[2]*RB[2];
      }
    }
  }
  else {                           //s,p
    // <s|s> = <s>
    RC[0] = RA[0]*RB[0];
    if (maxL == 1) {
      //<s|p> = <s|*(|s> + |p>) = <s> + <p>
      RC[1] = RA[0]*RB[1] + RA[1]*RB[0];
      if ((lorbA == 1)&&(lorbB == 1)) {
        //<p|p> = (<s| + <p|)*(|s> + |p>) = <s> + <p> + <d>
        RC[2] = RA[1]*RB[1];
      }
    }
  }
}
void FormProduct(int LA, int LB, std::vector<double> & prod, std::vector<double> & av, std::vector<double> & bv) {
  int maxL = std::max(LA,LB);
  if (maxL < 2) {
    //<s|s> = <s>
    prod[0] = av[0]*bv[0];
    if ((LA > 0)||(LB > 0)) {
      //<s|p> = <s|*(|s>+|p>) = <s> + <p>
      prod[1] = av[0]*bv[1] + av[1]*bv[0];
      if ((LA > 0)&&(LB > 0)) {
        //<p|p> = (<s|+<p|)*(|s>+|p>) = <s> + <p> + <d>
        prod[2] = av[1]*bv[1];
      }
    }
  }
  else if (maxL == 2) {
    //<s|d> = <s|*(|s>+|p>+|d>) = <s> + <p> + <d>
    prod[0] = av[0]*bv[0];
    prod[1] = av[0]*bv[1] + av[1]*bv[0];
    prod[2] = av[0]*bv[2] + av[2]*bv[0];
    if ((LA > 0)&&(LB > 0)) {
      //<p|d> = (<s|+<p|)*(|s>+|p>+|d>) = <s> + <p> + <d> + <f>
      prod[2] += av[1]*bv[1];
      prod[3] = av[1]*bv[2] + av[2]*bv[1];
      if (LA == LB) {
        //<d|d> = (<s|+<p|+<d|)*(|s>+|p>+|d>) = <s> + <p> + <d> + <f> + <g>
        prod[4] = av[2]*bv[2];
      }
    }
  }
  else if (maxL == 3) {
    //<s|f> = <s|*(|s>+|p>+|d>+|f>) = <s> + <p> + <d> + <f>
    prod[0] = av[0]*bv[0];
    prod[1] = av[0]*bv[1] + av[1]*bv[0];
    prod[2] = av[0]*bv[2] + av[2]*bv[0];
    prod[3] = av[0]*bv[3] + av[3]*bv[0];
    if ((LA > 0)&&(LB > 0)) {
      //<p|f> = (<s|+<p|)*(|s>+|p>+|d>+|f>) = <s> + <p> + <d> + <f> + <g>
      prod[2] += av[1]*bv[1];
      prod[3] += av[1]*bv[2] + av[2]*bv[1];
      prod[4] = av[1]*bv[3] + av[3]*bv[1];
      if ((LA > 1)&&(LB > 1)) {
        //<d|f> = (<s|+<p|+<d|)*(|s>+|p>+|d>+|f>) = <s> + <p> + <d> + <f> + <g> + <h>
        prod[4] += av[2]*bv[2];
        prod[5] = av[2]*bv[3] + av[3]*bv[2];
        if (LA == LB) {
          //<f|f> = (<s|+<p|+<d|+<f|)*(|s>+|p>+|d>+|f>) = <s> + <p> + <d> + <f> + <g> + <h> + <i>
          prod[5] = av[3]*bv[3];
        }
      }
    }
  }
  else if (maxL == 4) {
    //<s|g> = <s|*(|s>+|p>+|d>+|f>+|g>) = <s> + <p> + <d> + <f> + <g>
    prod[0] = av[0]*bv[0];
    prod[1] = av[0]*bv[1] + av[1]*bv[0];
    prod[2] = av[0]*bv[2] + av[2]*bv[0];
    prod[3] = av[0]*bv[3] + av[3]*bv[0];
    prod[4] = av[0]*bv[4] + av[4]*bv[0];
    if ((LA > 0)&&(LB > 0)) {
      //<p|g> = (<s|+<p|)*(|s>+|p>+|d>+|f>+|g>) = <s> + <p> + <d> + <f> + <g> + <h>
      prod[2] += av[1]*bv[1];
      prod[3] += av[1]*bv[2] + av[2]*bv[1];
      prod[4] += av[1]*bv[3] + av[3]*bv[1];
      prod[5] = av[1]*bv[4] + av[4]*bv[1];
      if ((LA > 1)&&(LB > 1)) {
        //<d|g> = (<s|+<p|+<d|)*(|s>+|p>+|d>+|f>+|g>) = <s> + <p> + <d> + <f> + <g> + <h> + <i>
        prod[4] += av[2]*bv[2];
        prod[5] += av[2]*bv[3] + av[3]*bv[2];
        prod[6] = av[2]*bv[4] + av[4]*bv[2];
        if ((LA > 1)&&(LB > 1)) {
          //<f|g> = (<s|+<p|+<d|+<f|)*(|s>+|p>+|d>+|f>+|g>) = <s> + <p> + <d> + <f> + <g> + <h> + <i> + <k>
          prod[6] += av[3]*bv[3];
          prod[7] = av[3]*bv[4] + av[4]*bv[3];
          if (LA == LB) {
            //<g|g> = (<s|+<p|+<d|+<f|+<g|)*(|s>+|p>+|d>+|f>+|g>) = <s> + <p> + <d> + <f> + <g> + <h> + <i> + <k> + <l>
            prod[8] = av[4]*bv[4];
          }
        }
      }
    }
  }
//end subroutine form_product
}
double dftr(int Lorb) {
  double DFTR = 1.0;
  if (Lorb == 2) {DFTR = 3.0;}
  else if (Lorb == 3) {DFTR = 15.0;}
  else if (Lorb == 4) {DFTR = 105.0;}
  else if (Lorb == 5) {DFTR = 945.0;}
  else if (Lorb == 6) {DFTR = 10395.0;}
  else if (Lorb == 7) {DFTR = 135135.0;}
  return DFTR;
}
double PartialOverlapGTO(double gamma, int Lorb) {
  //function calculating the partial overlap in one Cartesian direction
  //this is a C++ translation of olapp in xTB
  double sov = 0.0;
  if (Lorb%2 == 0) {sov = pow(0.5/gamma,0.5*double(Lorb))*dftr(Lorb/2);}
  return sov;
}
void TransfO(matrixE & matrix) {
  //getting the transformation matrix for d orbitals
  matrix = Identity(6);
  double aux = 1.0/sqrt(5.0);
  matrix(1,1) = aux;
  matrix(2,1) = aux;
  matrix(3,1) = aux;
  aux = sqrt3;
  matrix(1,2) = 0.5*aux;
  matrix(2,2) = -0.5*aux;
  matrix(1,3) = 0.5;
  matrix(2,3) = 0.5;
  matrix(3,3) = -1.0;
}
void dTransf(matrixE & smat, int li, int lj, matrixE & trfo, matrixE & saux) {
  //function converting integrals involving a d-shell to the 5d basis
  double sspher;
  //transformation not needed for pure s/p overlap -> do nothing
  if ((li > 1)||(lj > 1)) {                         //do something if there are d orbitals
    if ((li == 0)&&(lj == 2)) {                     //s|d
      for (size_t idorb = 1; idorb < 7; ++idorb) {
        sspher = 0.0;
        for (size_t icontract = 1; icontract < 7; ++icontract) {
          sspher += smat(1,icontract)*trfo(icontract,idorb);
        }
        saux(idorb,1) = sspher;
      }
      for (size_t idx = 1; idx < 6; ++idx) {
        smat(1,idx) = saux(idx + 1,1);
      }
      smat.resize(1,5);
    }
    else if ((li == 1)&&(lj == 2)) {                //p|d
      for (size_t iporb = 1; iporb < 4; ++iporb) {
        for (size_t idorb = 1; idorb < 7; ++idorb) {
          sspher = 0.0;
          for (size_t icontract = 1; icontract < 7; ++icontract) {
            sspher += smat(iporb,icontract)*trfo(icontract,idorb);
          }
          saux(idorb,iporb) = sspher;
        }
        for (size_t idx = 1; idx < 6; ++idx) {
          smat(iporb,idx) = saux(idx + 1,iporb);
        }
      }
      smat.resize(3,5);
    }
    if ((li == 2)&&(lj == 0)) {                     //d|s
      for (size_t idorb = 1; idorb < 7; ++idorb) {
        sspher = 0.0;
        for (size_t icontract = 1; icontract < 7; ++icontract) {
          sspher += trfo(icontract,idorb)*smat(icontract,1);
        }
        saux(1,idorb) = sspher;
      }
      for (size_t idx = 1; idx < 6; ++idx) {
        smat(idx,1) = saux(1,idx + 1);
      }
      smat.resize(5,1);
    }
    else if ((li == 2)&&(lj == 1)) {                //d|p
      for (size_t iporb = 1; iporb < 4; ++iporb) {
        for (size_t idorb = 1; idorb < 7; ++idorb) {
          sspher = 0.0;
          for (size_t icontract = 1; icontract < 7; ++icontract) {
            sspher += trfo(icontract,idorb)*smat(icontract,iporb);
          }
          saux(iporb,idorb) = sspher;
        }
        for (size_t idx = 1; idx < 6; ++idx) {
          smat(idx,iporb) = saux(iporb,idx + 1);
        }
      }
      smat.resize(5,3);
    }
    else if ((li == 2)&&(lj == 2)) {                //d|d
      //transposing smat in first dgemm is important for integrals other than S
      saux = trfo.trans()*smat;
      smat = saux*trfo;
      saux = smat;
      smat.resize(5,5);
      for (size_t idrow = 1; idrow < 6; ++idrow) {
        for (size_t idcol = 1; idcol < 6; ++idcol) {
          smat(idrow,idcol) = saux(idrow + 1,idcol + 1);
        }
      }
    }
  }
}
void dTransfT(tensor<matrixE,double> & intblock, int li, int lj, matrixE & trfo, matrixE & saux) {
  //function converting integrals involving a d-shell to the 5d basis
  double sspher;
  size_t nintegralclass = intblock.dim1();
  //transformation not needed for pure s/p overlap -> do nothing
  if ((li > 1)||(lj > 1)) {                         //do something if there are d orbitals
    if ((li == 0)&&(lj == 2)) {                     //s|d
      for (size_t idm1 = 1; idm1 < nintegralclass + 1; ++idm1) {
        for (size_t idorb = 1; idorb < 7; ++idorb) {
          sspher = 0.0;
          for (size_t icontract = 1; icontract < 7; ++icontract) {
            sspher += intblock(idm1,1,icontract)*trfo(icontract,idorb);
          }
          saux(idorb,1) = sspher;
        }
        for (size_t idx = 1; idx < 6; ++idx) {
          intblock(idm1,1,idx) = saux(idx + 1,1);
        }
      }
      intblock.resize(nintegralclass,1,5);
    }
    else if ((li == 1)&&(lj == 2)) {                //p|d
      for (size_t idm1 = 1; idm1 < nintegralclass + 1; ++idm1) {
        for (size_t iporb = 1; iporb < 4; ++iporb) {
          for (size_t idorb = 1; idorb < 7; ++idorb) {
            sspher = 0.0;
            for (size_t icontract = 1; icontract < 7; ++icontract) {
              sspher += intblock(idm1,iporb,icontract)*trfo(icontract,idorb);
            }
            saux(idorb,iporb) = sspher;
          }
          for (size_t idx = 1; idx < 6; ++idx) {
            intblock(idm1,iporb,idx) = saux(idx + 1,iporb);
          }
        }
      }
      intblock.resize(nintegralclass,3,5);
    }
    if ((li == 2)&&(lj == 0)) {                     //d|s
      for (size_t idm1 = 1; idm1 < nintegralclass + 1; ++idm1) {
        for (size_t idorb = 1; idorb < 7; ++idorb) {
          sspher = 0.0;
          for (size_t icontract = 1; icontract < 7; ++icontract) {
            sspher += trfo(icontract,idorb)*intblock(idm1,icontract,1);
          }
          saux(1,idorb) = sspher;
        }
        for (size_t idx = 1; idx < 6; ++idx) {
          intblock(idm1,idx,1) = saux(1,idx + 1);
        }
      }
      intblock.resize(nintegralclass,5,1);
    }
    else if ((li == 2)&&(lj == 1)) {                //d|p
      for (size_t idm1 = 1; idm1 < nintegralclass + 1; ++idm1) {
        for (size_t iporb = 1; iporb < 4; ++iporb) {
          for (size_t idorb = 1; idorb < 7; ++idorb) {
            sspher = 0.0;
            for (size_t icontract = 1; icontract < 7; ++icontract) {
              sspher += trfo(icontract,idorb)*intblock(idm1,icontract,iporb);
            }
            saux(iporb,idorb) = sspher;
          }
          for (size_t idx = 1; idx < 6; ++idx) {
            intblock(idm1,idx,iporb) = saux(iporb,idx + 1);
          }
        }
      }
      intblock.resize(nintegralclass,5,3);
    }
    else if ((li == 2)&&(lj == 2)) {                //d|d
      for (size_t idm1 = 1; idm1 < nintegralclass + 1; ++idm1) {
        for (size_t idorbB = 1; idorbB < 7; ++idorbB) {
          for (size_t idorbA = 1; idorbA < 7; ++idorbA) {
            sspher = 0.0;
            for (size_t icontract = 1; icontract < 7; ++icontract) {
              sspher += trfo(icontract,idorbA)*intblock(idm1,icontract,idorbB);
            }
            saux(idorbA,idorbB) = sspher;
          }
        }
        for (size_t idx = 1; idx < 7; ++idx) {
          for (size_t idy = 1; idy < 7; ++idy) {
            intblock(idm1,idx,idy) = saux(idx,idy);
          }
        }
        for (size_t idorbA = 1; idorbA < 7; ++idorbA) {
          for (size_t idorbB = 1; idorbB < 7; ++idorbB) {
            sspher = 0.0;
            for (size_t icontract = 1; icontract < 7; ++icontract) {
              sspher += intblock(idm1,idorbA,icontract)*trfo(icontract,idorbB);
            }
            saux(idorbA,idorbB) = sspher;
          }
        }
        for (size_t idx = 1; idx < 6; ++idx) {
          for (size_t idy = 1; idy < 6; ++idy) {
            intblock(idm1,idx,idy) = saux(idx + 1,idy + 1);
          }
        }
      }
      intblock.resize(nintegralclass,5,5);
    }
  }
}
double BuildOv(double alpha, std::vector<double> & RA, double beta, std::vector<double> & RB, std::vector<double> & RC, std::vector<double> & tvec, matrixE & val, std::vector<double> & aa, std::vector<double> & bb, std::vector<double> dd, int LorbA, int LorbB) {
  //function to get the overlap between two GTOs
  //this is a partial C++ translation of build_sdq_ints in xTB
  //RA -> center of first gaussian
  //RB -> center of second gaussian
  //RC -> reference point of momentum operator
  //alpha -> exponent of RA
  //beta -> exponent of RB
  //LorbA,LorbB -> azimuthal quantum numbers
  int lmax;
  int ii[3];
  int jj[3];
  std::vector<double> RE(3,0.0);
  double gamma;
  double kab;
  ii[0] = lx(LorbA);
  ii[1] = ly(LorbA);
  ii[2] = lz(LorbA);
  jj[0] = lx(LorbB);
  jj[1] = ly(LorbB);
  jj[2] = lz(LorbB);
  //RE is center of product gaussian with exponent gamma (gaussian product theorem)
  RE = GPTheorem(alpha,RA,beta,RB);
  //RC is reference point
  BuildKAB(alpha,RA,beta,RB,gamma,kab);
  for (size_t idx = 0; idx < 3; ++idx) {
    val(idx + 1,1) = 0.0;
    val(idx + 1,2) = 0.0;
    val(idx + 1,3) = 0.0;
    for (size_t idy = 0; idy < 5; ++idy) {
      aa[idy] = 0.0;
      bb[idy] = 0.0;
    }
    aa[ii[idx]] = 1.0;
    bb[jj[idx]] = 1.0;
    //Cartesian prefactor for bra, <a|
    BuildHShift(aa,RA[idx],RE[idx],ii[idx]);
    //Cartesian prefactor for ket, |b>
    BuildHShift(bb,RB[idx],RE[idx],jj[idx]);
    //get product
    Prod3(aa,bb,dd,ii[idx],jj[idx]);
    lmax = ii[idx] + jj[idx];
    for (size_t idy = 0; idy < lmax + 3; ++idy) {
      tvec[idy] = PartialOverlapGTO(gamma,idy);
    }
    for (size_t jdx = 0; jdx < lmax + 1; ++jdx) {
      val(idx + 1,1) += dd[jdx]*tvec[jdx];
    }
  }
  gamma = kab*val(1,1)*val(2,1)*val(3,1);
  return gamma;
}
std::vector<double> BuildOvDipQuadr(double alpha, std::vector<double> & RA, double beta, std::vector<double> & RB, std::vector<double> & RC, std::vector<double> & tvec, matrixE & val, std::vector<double> & aa, std::vector<double> & bb, std::vector<double> dd, int LorbA, int LorbB) {
  //function to get the overlap, dipole and quadrupole moment integrals between two GTOs
  //this is a C++ translation of build_sdq_ints in xTB
  //RA -> center of first gaussian
  //RB -> center of second gaussian
  //RC -> reference point of momentum operator
  //alpha -> exponent of RA
  //beta -> exponent of RB
  //LorbA,LorbB -> azimuthal quantum numbers
  int lmax;
  int ii[3];
  int jj[3];
  double va[3];
  std::vector<double> RE(3);
  double gamma;
  double kab;
  double RD[3];
  ii[0] = lx(LorbA);
  ii[1] = ly(LorbA);
  ii[2] = lz(LorbA);
  jj[0] = lx(LorbB);
  jj[1] = ly(LorbB);
  jj[2] = lz(LorbB);
  //RE is center of product gaussian with exponent gamma (gaussian product theorem)
  RE = GPTheorem(alpha,RA,beta,RB);
  //RC is reference point
  BuildKAB(alpha,RA,beta,RB,gamma,kab);
  for (size_t idx = 0; idx < 3; ++idx) {
    RD[idx] = RE[idx] - RC[idx];
    val(idx + 1,1) = 0.0;
    val(idx + 1,2) = 0.0;
    val(idx + 1,3) = 0.0;
    for (size_t idy = 0; idy < 5; ++idy) {
      aa[idy] = 0.0;
      bb[idy] = 0.0;
    }
    aa[ii[idx]] = 1.0;
    bb[jj[idx]] = 1.0;
    //Cartesian prefactor for bra, <a|
    BuildHShift(aa,RA[idx],RE[idx],ii[idx]);
    //Cartesian prefactor for ket, |b>
    BuildHShift(bb,RB[idx],RE[idx],jj[idx]);
    //get product
    Prod3(aa,bb,dd,ii[idx],jj[idx]);
    lmax = ii[idx] + jj[idx];
    for (size_t idy = 0; idy < lmax + 3; ++idy) {
      tvec[idy] = PartialOverlapGTO(gamma,idy);
    }
    for (size_t jdx = 0; jdx < lmax + 1; ++jdx) {
      va[0] = tvec[jdx];                                                                  //<a|b>
      va[1] = tvec[jdx + 1] + RD[idx]*tvec[jdx];                                          //<a|x|b>
      va[2] = tvec[jdx + 2] + 2.0*RD[idx]*tvec[jdx + 1] + RD[idx]*RD[idx]*tvec[jdx];      //<a|x²|b>
      for (size_t jdy = 0; jdy < 3; ++jdy) {
        val(idx + 1,jdy + 1) += dd[jdx]*va[jdy];
      }
    }
  }
  tvec[0] = kab*val(1,1)*val(2,1)*val(3,1);
  tvec[1] = kab*val(1,2)*val(2,1)*val(3,1);
  tvec[2] = kab*val(1,1)*val(2,2)*val(3,1);
  tvec[3] = kab*val(1,1)*val(2,1)*val(3,2);
  tvec[4] = kab*val(1,3)*val(2,1)*val(3,1);
  tvec[5] = kab*val(1,1)*val(2,3)*val(3,1);
  tvec[6] = kab*val(1,1)*val(2,1)*val(3,3);
  tvec[7] = kab*val(1,2)*val(2,2)*val(3,1);
  tvec[8] = kab*val(1,2)*val(2,1)*val(3,2);
  tvec[9] = kab*val(1,1)*val(2,2)*val(3,2);
  return tvec;
}
void gBuildOvDipQuadr(std::vector<double> & s1d, std::vector<double> & s3d, matrixE & ds3d, std::vector<double> & rp, std::vector<double> & rc, std::vector<double> & ri, std::vector<double> & rj, matrixE & val, matrixE & grad, std::vector<double> & vi, std::vector<double> & vj, std::vector<double> & vv, std::vector<double> & gi, std::vector<double> & gg, int LorbA, int LorbB, double ai, double aj) {
  //gradients of the function above
  double rpc;
  int ii[3];
  int jj[3];
  ii[0] = lx(LorbA);
  ii[1] = ly(LorbA);
  ii[2] = lz(LorbA);
  jj[0] = lx(LorbB);
  jj[1] = ly(LorbB);
  jj[2] = lz(LorbB);
  for (size_t idk = 1; idk < 4; ++idk) {
    for (size_t idl = 1; idl < 4; ++idl) {
      val(idk,idl) = 0.0;
      grad(idk,idl) = 0.0;
      vi[idl - 1] = 0.0;
      vj[idl - 1] = 0.0;
      vv[idl - 1] = 0.0;
      gi[idl - 1] = 0.0;
      gg[idl - 1] = 0.0;
    }
    for (size_t idl = 3; idl < 7; ++idl) {
      vi[idl] = 0.0;
      vj[idl] = 0.0;
      vv[idl] = 0.0;
      gi[idl] = 0.0;
      gg[idl] = 0.0;
    }
    for (size_t idl = 7; idl < 13; ++idl) {
      vv[idl] = 0.0;
      gg[idl] = 0.0;
    }
    rpc = rp[idk - 1] - rc[idk - 1];
    vi[ii[idk - 1]] = 1.0;
    vj[jj[idk - 1]] = 1.0;
    gi[ii[idk - 1] + 1] = 2.0*ai;
    if (ii[idk - 1] > 0) {gi[ii[idk - 1] - 1] = -ii[idk - 1];}
    HorizontalShift(gi,rp[idk - 1] - ri[idk - 1],ii[idk - 1] - 1);
    HorizontalShift(gi,rp[idk - 1] - ri[idk - 1],ii[idk - 1] + 1);
    HorizontalShift(vi,rp[idk - 1] - ri[idk - 1],ii[idk - 1]);
    HorizontalShift(vj,rp[idk - 1] - rj[idk - 1],jj[idk - 1]);
    FormProduct(ii[idk - 1],jj[idk - 1],vv,vi,vj);
    FormProduct(ii[idk - 1] + 1,jj[idk - 1],gg,gi,vj);
    for (size_t idl = 0; idl < ii[idk - 1] + jj[idk - 1] + 2; ++idl) {
      val(idk,1) += s1d[idl]*vv[idl];
      val(idk,2) += (s1d[idl + 1] + rpc*s1d[idl])*vv[idl];
      val(idk,3) += (s1d[idl + 2] + 2.0*rpc*s1d[idl + 1] + rpc*rpc*s1d[idl])*vv[idl];
      grad(idk,1) += s1d[idl]*gg[idl];
      grad(idk,2) += (s1d[idl + 1] + rpc*s1d[idl])*gg[idl];
      grad(idk,3) += (s1d[idl + 2] + 2.0*rpc*s1d[idl + 1] + rpc*rpc*s1d[idl])*gg[idl];
    }
  }
  s3d[0] = val(1,1)*val(2,1)*val(3,1);
  s3d[1] = val(1,2)*val(2,1)*val(3,1);
  s3d[2] = val(1,1)*val(2,2)*val(3,1);
  s3d[3] = val(1,1)*val(2,1)*val(3,2);
  s3d[4] = val(1,3)*val(2,1)*val(3,1);
  s3d[5] = val(1,1)*val(2,3)*val(3,1);
  s3d[6] = val(1,1)*val(2,1)*val(3,3);
  s3d[7] = val(1,2)*val(2,2)*val(3,1);
  s3d[8] = val(1,2)*val(2,1)*val(3,2);
  s3d[9] = val(1,1)*val(2,2)*val(3,2);
  ds3d(1,1)  = grad(1,1)*val(2,1)*val(3,1);
  ds3d(2,1)  = val(1,1)*grad(2,1)*val(3,1);
  ds3d(3,1)  = val(1,1)*val(2,1)*grad(3,1);
  ds3d(1,2)  = grad(1,2)*val(2,1)*val(3,1);
  ds3d(2,2)  = val(1,2)*grad(2,1)*val(3,1);
  ds3d(3,2)  = val(1,2)*val(2,1)*grad(3,1);
  ds3d(1,3)  = grad(1,1)*val(2,2)*val(3,1);
  ds3d(2,3)  = val(1,1)*grad(2,2)*val(3,1);
  ds3d(3,3)  = val(1,1)*val(2,2)*grad(3,1);
  ds3d(1,4)  = grad(1,1)*val(2,1)*val(3,2);
  ds3d(2,4)  = val(1,1)*grad(2,1)*val(3,2);
  ds3d(3,4)  = val(1,1)*val(2,1)*grad(3,2);
  ds3d(1,5)  = grad(1,3)*val(2,1)*val(3,1);
  ds3d(2,5)  = val(1,3)*grad(2,1)*val(3,1);
  ds3d(3,5)  = val(1,3)*val(2,1)*grad(3,1);
  ds3d(1,6)  = grad(1,1)*val(2,3)*val(3,1);
  ds3d(2,6)  = val(1,1)*grad(2,3)*val(3,1);
  ds3d(3,6)  = val(1,1)*val(2,3)*grad(3,1);
  ds3d(1,7)  = grad(1,1)*val(2,1)*val(3,3);
  ds3d(2,7)  = val(1,1)*grad(2,1)*val(3,3);
  ds3d(3,7)  = val(1,1)*val(2,1)*grad(3,3);
  ds3d(1,8)  = grad(1,2)*val(2,2)*val(3,1);
  ds3d(2,8)  = val(1,2)*grad(2,2)*val(3,1);
  ds3d(3,8)  = val(1,2)*val(2,2)*grad(3,1);
  ds3d(1,9)  = grad(1,2)*val(2,1)*val(3,2);
  ds3d(2,9)  = val(1,2)*grad(2,1)*val(3,2);
  ds3d(3,9)  = val(1,2)*val(2,1)*grad(3,2);
  ds3d(1,10) = grad(1,1)*val(2,2)*val(3,2);
  ds3d(2,10) = val(1,1)*grad(2,2)*val(3,2);
  ds3d(3,10) = val(1,1)*val(2,2)*grad(3,2);
}
void TwoCenterOvBlock(matrixE & ovlp, size_t nprimA, size_t nprimB, std::vector<double> & RA, std::vector<double> & exponentsA, std::vector<double> & coeffsA, size_t LA, std::vector<double> & RB, std::vector<double> & exponentsB, std::vector<double> & coeffsB, size_t LB, std::vector<double> & RP, std::vector<double> & tvec, matrixE & val, std::vector<double> & aa, std::vector<double> & bb, std::vector<double> dd, double thresholdzero = 1e-7, double maxsqdist = 2000.0, double intcut = 25.0) {
  //function to get an overlap block between the sub-shells of two atoms described by GTOs (s,s; s,p; p,p; ...)
  //this is a partial C++ translation of get_multiints in xTB
  //nprimA,nprimB are the number of primitives for each AO subshell of each atom; because different subshells might be described by a different number of GTO primitives, I kept the function as is from xTB
  //RA,RB,... is the position vector for centers
  //exponentsA,exponentsB are vectors with exponents
  //coeffsA,coeffsB are vectors is GTO coefficients
  //contrary to Grimme's program, the vectors I pass contain exactly what I need
  double alpha;
  double beta;
  double est;
  double cA;
  double cB;
  double cAprime;
  double cBprime;
  double c2;
  double ovl;
  double RAB2 = 0.0;
  for (size_t idcoord = 0; idcoord < RA.size(); ++idcoord) {
    RAB2 += (RA[idcoord] - RB[idcoord])*(RA[idcoord] - RB[idcoord]);
  }
  int LAoffset = OrbitalOffset(LA);
  int LBoffset = OrbitalOffset(LB);
  int naoA = 2*LA + 1;
  int naoB = 2*LB + 1;
  if (LA == 2) {++naoA;}
  if (LB == 2) {++naoB;}
  ovlp.resize(naoA,naoB);
  ovlp.zero();
  if (RAB2 < maxsqdist + thresholdzero) {
    for (size_t Aprim = 0; Aprim < nprimA; ++Aprim) {
      alpha = exponentsA[Aprim];
      cA = coeffsA[Aprim];
      for (size_t Bprim = 0; Bprim < nprimB; ++Bprim) {
        beta = exponentsB[Bprim];
        cB = coeffsB[Bprim];
        est = RAB2*alpha*beta/(alpha + beta);
        if (est > intcut) {continue;}
        //get integrals for different magnetic quantum numbers
        for (size_t mlA = 1; mlA < naoA + 1; ++mlA) {
          cAprime = cA;
          if ((LAoffset + mlA > 7)&&(LAoffset + mlA < 11)) {cAprime *= sqrt3;}
          else if (LAoffset + mlA > 13) {
            cAprime *= sqrt(5.0);
            if (LAoffset + mlA == 20) {cAprime *= sqrt3;}
          }
          for (size_t mlB = 1; mlB < naoB + 1; ++mlB) {
            cBprime = cB;
            if ((LBoffset + mlB > 7)&&(LBoffset + mlB < 11)) {cBprime *= sqrt3;}
            else if (LBoffset + mlB > 13) {
              cBprime *= sqrt(5.0);
              if (LBoffset + mlB == 20) {cBprime *= sqrt3;}
            }
            ovl = BuildOv(alpha,RA,beta,RB,RP,tvec,val,aa,bb,dd,LAoffset + mlA,LBoffset + mlB);
            c2 = cBprime*cAprime;
            ovlp(mlA,mlB) += ovl*c2;
          }                      //mlB
        }                        //mlA
      }                          //Bprim
    }                            //Aprim
  }
}
void TwoCenterOvDipQuadrBlock(matrixE & ovlp, tensor<matrixE,double> & dple, tensor<matrixE,double> & quadr, size_t nprimA, size_t nprimB, std::vector<double> & RA, std::vector<double> & exponentsA, std::vector<double> & coeffsA, size_t LA, std::vector<double> & RB, std::vector<double> & exponentsB, std::vector<double> & coeffsB, size_t LB, std::vector<double> & RP, std::vector<double> & tvec, matrixE & val, std::vector<double> & aa, std::vector<double> & bb, std::vector<double> dd, double thresholdzero = 1e-7, double maxsqdist = 2000.0, double intcut = 25.0) {
  //function to get a block of overlap, dipole and quadrupole moment integrals between the sub-shells of two atoms described by GTOs (s,s; s,p; p,p; ...)
  //this is the closest translation to C++ of get_multiints in xTB
  //nprimA,nprimB are the number of primitives for each AO subshell of each atom; because different subshells might be described by a different number of GTO primitives, I kept the function as is from xTB
  //RA,RB,... is the position vector for centers
  //exponentsA,exponentsB are vectors with exponents
  //coeffsA,coeffsB are vectors is GTO coefficients
  //contrary to Grimme's program, the vectors I pass contain exactly what I need
  //tensors structured as described above
  double alpha;
  double beta;
  double est;
  double cA;
  double cB;
  double cAprime;
  double cBprime;
  double c2;
  double RAB2 = 0.0;
  std::vector<double> integrals;
  for (size_t idcoord = 0; idcoord < RA.size(); ++idcoord) {
    RAB2 += (RA[idcoord] - RB[idcoord])*(RA[idcoord] - RB[idcoord]);
  }
  int LAoffset = OrbitalOffset(LA);
  int LBoffset = OrbitalOffset(LB);
  int naoA = 2*LA + 1;
  int naoB = 2*LB + 1;
  if (LA == 2) {++naoA;}
  if (LB == 2) {++naoB;}
  ovlp.resize(naoA,naoB);
  dple.resize(3,naoA,naoB);
  quadr.resize(6,naoA,naoB);
  ovlp.zero();
  dple.zero();
  quadr.zero();
  if (RAB2 < maxsqdist + thresholdzero) {
    for (size_t Aprim = 0; Aprim < nprimA; ++Aprim) {
      alpha = exponentsA[Aprim];
      cA = coeffsA[Aprim];
      for (size_t Bprim = 0; Bprim < nprimB; ++Bprim) {
        beta = exponentsB[Bprim];
        cB = coeffsB[Bprim];
        est = RAB2*alpha*beta/(alpha + beta);
        if (est > intcut) {continue;}
        //get integrals for different magnetic quantum numbers
        for (size_t mlA = 1; mlA < naoA + 1; ++mlA) {
          cAprime = cA;
          if ((LAoffset + mlA > 7)&&(LAoffset + mlA < 11)) {cAprime *= sqrt3;}
          else if (LAoffset + mlA > 13) {
            cAprime *= sqrt(5.0);
            if (LAoffset + mlA == 20) {cAprime *= sqrt3;}
          }
          for (size_t mlB = 1; mlB < naoB + 1; ++mlB) {
            cBprime = cB;
            if ((LBoffset + mlB > 7)&&(LBoffset + mlB < 11)) {cBprime *= sqrt3;}
            else if (LBoffset + mlB > 13) {
              cBprime *= sqrt(5.0);
              if (LBoffset + mlB == 20) {cBprime *= sqrt3;}
            }
            integrals = BuildOvDipQuadr(alpha,RA,beta,RB,RP,tvec,val,aa,bb,dd,OrbitalOffset(LA) + mlA,OrbitalOffset(LB) + mlB);
            c2 = cBprime*cAprime;
            ovlp(mlA,mlB) += integrals[0]*c2;
            for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
              dple(idcoord,mlA,mlB) += integrals[idcoord]*c2;
            }
            for (size_t idcoord = 1; idcoord < 7; ++idcoord) {
              quadr(idcoord,mlA,mlB) += integrals[idcoord + 3]*c2;
            }
          }                      //mlB
        }                        //mlA
      }                          //Bprim
    }                            //Aprim
  }
}
void TwoCenterDipQuadrBlock(tensor<matrixE,double> & dple, tensor<matrixE,double> & quadr, size_t nprimA, size_t nprimB, std::vector<double> & RA, std::vector<double> & exponentsA, std::vector<double> & coeffsA, size_t LA, std::vector<double> & RB, std::vector<double> & exponentsB, std::vector<double> & coeffsB, size_t LB, std::vector<double> & RP, std::vector<double> & tvec, matrixE & val, std::vector<double> & aa, std::vector<double> & bb, std::vector<double> dd, double thresholdzero = 1e-7, double maxsqdist = 2000.0, double intcut = 25.0) {
  //function to get a block of dipole and quadrupole moment integrals between the sub-shells of two atoms described by GTOs (s,s; s,p; p,p; ...)
  //this is a partial C++ translation of get_multiints in xTB
  //nprimA,nprimB are the number of primitives for each AO subshell of each atom; because different subshells might be described by a different number of GTO primitives, I kept the function as is from xTB
  //RA,RB,... is the position vector for centers
  //exponentsA,exponentsB are vectors with exponents
  //coeffsA,coeffsB are vectors is GTO coefficients
  //contrary to Grimme's program, the vectors I pass contain exactly what I need
  //tensors structured as described above
  double alpha;
  double beta;
  double est;
  double cA;
  double cB;
  double cAprime;
  double cBprime;
  double c2;
  double RAB2 = 0.0;
  std::vector<double> integrals;
  for (size_t idcoord = 0; idcoord < RA.size(); ++idcoord) {
    RAB2 += (RA[idcoord] - RB[idcoord])*(RA[idcoord] - RB[idcoord]);
  }
  int LAoffset = OrbitalOffset(LA);
  int LBoffset = OrbitalOffset(LB);
  int naoA = 2*LA + 1;
  int naoB = 2*LB + 1;
  if (LA == 2) {++naoA;}
  if (LB == 2) {++naoB;}
  dple.resize(3,naoA,naoB);
  quadr.resize(6,naoA,naoB);
  dple.zero();
  quadr.zero();
  if (RAB2 < maxsqdist + thresholdzero) {
    for (size_t Aprim = 0; Aprim < nprimA; ++Aprim) {
      alpha = exponentsA[Aprim];
      cA = coeffsA[Aprim];
      for (size_t Bprim = 0; Bprim < nprimB; ++Bprim) {
        beta = exponentsB[Bprim];
        cB = coeffsB[Bprim];
        est = RAB2*alpha*beta/(alpha + beta);
        if (est > intcut) {continue;}
        //get integrals for different magnetic quantum numbers
        for (size_t mlA = 1; mlA < naoA + 1; ++mlA) {
          cAprime = cA;
          if ((LAoffset + mlA > 7)&&(LAoffset + mlA < 11)) {cAprime *= sqrt3;}
          else if (LAoffset + mlA > 13) {
            cAprime *= sqrt(5.0);
            if (LAoffset + mlA == 20) {cAprime *= sqrt3;}
          }
          for (size_t mlB = 1; mlB < naoB + 1; ++mlB) {
            cBprime = cB;
            if ((LBoffset + mlB > 7)&&(LBoffset + mlB < 11)) {cBprime *= sqrt3;}
            else if (LBoffset + mlB > 13) {
              cBprime *= sqrt(5.0);
              if (LBoffset + mlB == 20) {cBprime *= sqrt3;}
            }
            //not worth an implementation just for dipole and quadrupole moment integrals because the overlap block brings no extra overhead
            integrals = BuildOvDipQuadr(alpha,RA,beta,RB,RP,tvec,val,aa,bb,dd,OrbitalOffset(LA) + mlA,OrbitalOffset(LB) + mlB);
            c2 = cBprime*cAprime;
            for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
              dple(idcoord,mlA,mlB) += integrals[idcoord]*c2;
            }
            for (size_t idcoord = 1; idcoord < 7; ++idcoord) {
              quadr(idcoord,mlA,mlB) += integrals[idcoord + 3]*c2;
            }
          }                      //mlB
        }                        //mlA
      }                          //Bprim
    }                            //Aprim
  }
}
void gTwoCenterOvDipQuadrBlock(std::vector<matrixE> & SDQ, matrix<matrixE> & gSDQ, std::vector<double> & rA, std::vector<double> & rB, std::vector<double> & rP, std::vector<double> & rAB, double rAB2, size_t nprimA, size_t nprimB, std::vector<double> & exponentsA, std::vector<double> & exponentsB, std::vector<double> & coeffsA, std::vector<double> & coeffsB, std::vector<double> & ints, std::vector<double> & tk, matrixE & gints, matrixE & val, matrixE & grad, std::vector<double> & vi, std::vector<double> & vj, std::vector<double> & vv, std::vector<double> & gi, std::vector<double> & gg, size_t LA, size_t LB, double maxsqdist = 2000.0, double intcut = 25.0) {
  //function that gets the gradient for two center GTO integrals, namely overlap, dipole and quadrupole
  if (rAB2 <= maxsqdist) {
    int n1 = SDQ.size();
    int nr = gSDQ.rows();
    int nc = gSDQ.cols();
    for (size_t mlA = 0; mlA < 6; ++mlA) {
      for (size_t mlB = 0; mlB <= mlA; ++mlB) {
        for (size_t idcoord = 0; idcoord < n1; ++idcoord) {
          SDQ[idcoord](mlA + 1,mlB + 1) = 0.0;
          SDQ[idcoord](mlB + 1,mlA + 1) = 0.0;
        }
        for (size_t idr = 0; idr < nr; ++idr) {
          for (size_t idc = 0; idc < nc; ++idc) {
            gSDQ(idr + 1,idc + 1)(mlA + 1,mlB + 1) = 0.0;
            gSDQ(idr + 1,idc + 1)(mlB + 1,mlA + 1) = 0.0;
          }
        }
      }
    }
    nc = std::min(int(gSDQ.cols()),10);
    int Lsum3 = LA + LB + 3;
    int LAoffset = OrbitalOffset(LA);
    int LBoffset = OrbitalOffset(LB);
    int naoA = 2*LA + 1;
    int naoB = 2*LB + 1;
    if (LA == 2) {++naoA;}
    if (LB == 2) {++naoB;}
    double alpha;
    double beta;
    double ab;
    double est;
    double kab;
    double cA;
    double cB;
    double cAprime;
    double cBprime;
    double cc;
    //loop over primitives
    for (size_t iprimA = 0; iprimA < nprimA; ++iprimA) {
      alpha = exponentsA[iprimA];
      cA = coeffsA[iprimA];
      for (size_t iprimB = 0; iprimB < nprimB; ++iprimB) {
        beta = exponentsB[iprimB];
        cB = coeffsB[iprimB];
        ab = 1.0/(alpha + beta);
        est = alpha*beta*rAB2*ab;
        if (est > intcut) {continue;}
        kab = exp(-est)*(sqrtpi*sqrt(ab))*(sqrtpi*sqrt(ab))*(sqrtpi*sqrt(ab));
        rP[0] = (alpha*rA[0] + beta*rB[0])*ab;
        rP[1] = (alpha*rA[1] + beta*rB[1])*ab;
        rP[2] = (alpha*rA[2] + beta*rB[2])*ab;
        for (size_t idk = 0; idk < Lsum3 + 1; ++idk) {
          tk[idk] = PartialOverlapGTO(alpha + beta,idk);
        }
        //now ze gradient
        for (size_t mlA = 0; mlA < naoA; ++mlA) {
          cAprime = cA;
          if ((LAoffset + mlA > 6)&&(LAoffset + mlA < 10)) {cAprime *= sqrt3;}
          else if (LAoffset + mlA > 12) {
            cAprime *= sqrt(5.0);
            if (LAoffset + mlA == 19) {cAprime *= sqrt3;}
          }
          for (size_t mlB = 0; mlB < naoB; ++mlB) {
            cBprime = cB;
            if ((LBoffset + mlB > 6)&&(LBoffset + mlB < 10)) {cBprime *= sqrt3;}
            else if (LBoffset + mlB > 12) {
              cBprime *= sqrt(5.0);
              if (LBoffset + mlB == 19) {cBprime *= sqrt3;}
            }
            cc = cBprime*cAprime*kab;
            gBuildOvDipQuadr(tk,ints,gints,rP,rB,rA,rB,val,grad,vi,vj,vv,gi,gg,LAoffset + mlA + 1,LBoffset + mlB + 1,alpha,beta);
            for (size_t idim = 0; idim < n1; ++idim) {
              SDQ[idim](mlB + 1,mlA + 1) += cc*ints[idim];
            }
            for (size_t idr = 0; idr < nr; ++idr) {
              for (size_t idc = 0; idc < nc; ++idc) {
                gSDQ(idr + 1,idc + 1)(mlB + 1,mlA + 1) += cc*gints(idr + 1,idc + 1);
              }
            }
          }
        }
      }
    }
    //shift the integrals from center A to center B
    for (size_t mlA = 0; mlA < naoA; ++mlA) {
      for (size_t mlB = 0; mlB < naoB; ++mlB) {
        for (size_t idcoord = 1; idcoord < 4; ++idcoord) {
          gSDQ(idcoord,11)(mlB + 1,mlA + 1) = gSDQ(idcoord,2)(mlB + 1,mlA + 1) - rAB[0]*gSDQ(idcoord,1)(mlB + 1,mlA + 1);
          gSDQ(idcoord,12)(mlB + 1,mlA + 1) = gSDQ(idcoord,3)(mlB + 1,mlA + 1) - rAB[1]*gSDQ(idcoord,1)(mlB + 1,mlA + 1);
          gSDQ(idcoord,13)(mlB + 1,mlA + 1) = gSDQ(idcoord,4)(mlB + 1,mlA + 1) - rAB[2]*gSDQ(idcoord,1)(mlB + 1,mlA + 1);
          gSDQ(idcoord,14)(mlB + 1,mlA + 1) = gSDQ(idcoord,5)(mlB + 1,mlA + 1) - 2.0*rAB[0]*gSDQ(idcoord,2)(mlB + 1,mlA + 1) + rAB[0]*rAB[0]*gSDQ(idcoord,1)(mlB + 1,mlA + 1);
          gSDQ(idcoord,15)(mlB + 1,mlA + 1) = gSDQ(idcoord,6)(mlB + 1,mlA + 1) - 2.0*rAB[1]*gSDQ(idcoord,3)(mlB + 1,mlA + 1) + rAB[1]*rAB[1]*gSDQ(idcoord,1)(mlB + 1,mlA + 1);
          gSDQ(idcoord,16)(mlB + 1,mlA + 1) = gSDQ(idcoord,7)(mlB + 1,mlA + 1) - 2.0*rAB[2]*gSDQ(idcoord,4)(mlB + 1,mlA + 1) + rAB[2]*rAB[2]*gSDQ(idcoord,1)(mlB + 1,mlA + 1);
          gSDQ(idcoord,17)(mlB + 1,mlA + 1) = gSDQ(idcoord,8)(mlB + 1,mlA + 1)  - rAB[0]*gSDQ(idcoord,3)(mlB + 1,mlA + 1) - rAB[1]*gSDQ(idcoord,2)(mlB + 1,mlA + 1) + rAB[0]*rAB[1]*gSDQ(idcoord,1)(mlB + 1,mlA + 1);
          gSDQ(idcoord,18)(mlB + 1,mlA + 1) = gSDQ(idcoord,9)(mlB + 1,mlA + 1)  - rAB[0]*gSDQ(idcoord,4)(mlB + 1,mlA + 1) - rAB[2]*gSDQ(idcoord,2)(mlB + 1,mlA + 1) + rAB[0]*rAB[2]*gSDQ(idcoord,1)(mlB + 1,mlA + 1);
          gSDQ(idcoord,19)(mlB + 1,mlA + 1) = gSDQ(idcoord,10)(mlB + 1,mlA + 1) - rAB[1]*gSDQ(idcoord,4)(mlB + 1,mlA + 1) - rAB[2]*gSDQ(idcoord,3)(mlB + 1,mlA + 1) + rAB[1]*rAB[2]*gSDQ(idcoord,1)(mlB + 1,mlA + 1);
          gSDQ(idcoord,10 + idcoord)(mlB + 1,mlA + 1) -= SDQ[0](mlB + 1,mlA + 1);
        }
        gSDQ(1,14)(mlB + 1,mlA + 1) += -2.0*SDQ[1](mlB + 1,mlA + 1) + 2.0*rAB[0]*SDQ[0](mlB + 1,mlA + 1);
        gSDQ(2,15)(mlB + 1,mlA + 1) += -2.0*SDQ[2](mlB + 1,mlA + 1) + 2.0*rAB[1]*SDQ[0](mlB + 1,mlA + 1);
        gSDQ(3,16)(mlB + 1,mlA + 1) += -2.0*SDQ[3](mlB + 1,mlA + 1) + 2.0*rAB[2]*SDQ[0](mlB + 1,mlA + 1);
        gSDQ(1,17)(mlB + 1,mlA + 1) += -SDQ[2](mlB + 1,mlA + 1) + rAB[1]*SDQ[0](mlB + 1,mlA + 1);
        gSDQ(2,17)(mlB + 1,mlA + 1) += -SDQ[1](mlB + 1,mlA + 1) + rAB[0]*SDQ[0](mlB + 1,mlA + 1);
        gSDQ(1,18)(mlB + 1,mlA + 1) += -SDQ[3](mlB + 1,mlA + 1) + rAB[2]*SDQ[0](mlB + 1,mlA + 1);
        gSDQ(3,18)(mlB + 1,mlA + 1) += -SDQ[1](mlB + 1,mlA + 1) + rAB[0]*SDQ[0](mlB + 1,mlA + 1);
        gSDQ(2,19)(mlB + 1,mlA + 1) += -SDQ[3](mlB + 1,mlA + 1) + rAB[2]*SDQ[0](mlB + 1,mlA + 1);
        gSDQ(3,19)(mlB + 1,mlA + 1) += -SDQ[2](mlB + 1,mlA + 1) + rAB[1]*SDQ[0](mlB + 1,mlA + 1);
      }
    }
  }
}

#endif //_Two_Center_Overlaps_GTO_
