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

#ifndef _Density_Matrix_
#define _Density_Matrix_
#include <vector>
#include "math/MatrixPackage.hpp"

//description:
//calculation of density matrix

void Dens(matrixE & pds, std::vector<double> & occ, matrixE & CMO, matrixE & Daux, bool FermiSmear = false, double thresholdzero = 1.0e-6) {
  //function that calculates the density matrix 
  size_t nrows = CMO.rows();
  size_t ncols = CMO.cols();
  size_t occsz = occ.size();
  double nelec;
  if (!FermiSmear) {             //either RHF or UHF, i.e., occ is vector of ones or of twos
    CMO.getBlock(Daux,1,1,nrows,occsz);
    if (occsz > 0) {nelec = 1.0 + double(occ[0] == 2);}
    else {nelec = 0.0;}
    pds = Daux*Daux.trans()*nelec;
  }
  else {                         //Fermi smearing is used, meaning there is fractional occupation of orbitals
    //determine the last truly closed-shell orbital
    size_t lastcs = occsz - 1;
    size_t lastocc = occsz - 1;
    for (size_t idao = 1; idao < occsz; ++idao) {
      if ((2.0 - occ[idao] > thresholdzero)&&(lastcs > idao)) {lastcs = idao - 1;}
      if (occ[idao] < thresholdzero) {
        lastocc = idao - 1;
        break;
      }
    }
    //closed-shell part
    CMO.getBlock(Daux,1,1,nrows,lastcs + 1);
    nelec = 2.0;
    pds = Daux*Daux.trans()*nelec;
    //now the open-shell component
    if (lastcs != lastocc) {
      for (size_t irw = 0; irw < nrows; ++irw) {
        for (size_t icl = 0; icl < ncols; ++icl) {
          for (size_t idx = lastcs + 1; idx < lastocc + 1; ++idx) {
            pds(irw + 1,icl + 1) += CMO(irw + 1,idx + 1)*CMO(icl + 1,idx + 1)*occ[idx];
          }
        }
      }
    }
  }
}
void generalizedDens(matrixE & pds, std::vector<double> & occ, matrixE & CMOA, matrixE & CMOB, double thresholdzero = 1.0e-10) {
  //generalization of the function that calculates the density matrix
  //this function may be used when the MOs are not organized in two blocks (occ,vir) but in several blocks determined by the occupation
  //vector (occA,virA,occB,virB,...)
  //or two contract two different sets of "orbitals" with shared occupancies
  size_t nrows = CMOA.rows();
  size_t ncols = CMOB.cols();
  size_t ncontr = CMOA.cols();
  double nocc;
  if (ncontr != CMOB.rows()) {std::cout << "WARNING: Density.hpp: generalizedDens(): contraction dimensions not matching between matrices" << std::endl;}
  pds.resize(nrows,ncols);
  for (size_t idmu = 0; idmu < nrows; ++idmu) {
    for (size_t idalpha = 0; idalpha < ncols; ++idalpha) {
      pds(idmu + 1,idalpha + 1) = 0.0;
      for (size_t idcontr = 0; idcontr < ncontr; ++idcontr) {
        nocc = occ[idcontr];
        if (nocc < thresholdzero) {continue;}
        pds(idmu + 1,idalpha + 1) += nocc*CMOA(idmu + 1,idcontr + 1)*CMOB(idcontr + 1,idalpha + 1);
      }
    }
  }
}
void DensW(matrixE & pds, std::vector<double> & eocc, matrixE & CMO, double thresholdzero = 1.0e-8) {
  //function that calculates the density matrix 
  size_t nrows = CMO.rows();
  size_t ncols = CMO.cols();
  size_t occsz = eocc.size();
  for (size_t irw = 0; irw < nrows; ++irw) {
    for (size_t icl = 0; icl < ncols; ++icl) {
      pds(irw + 1,icl + 1) = 0.0;
      for (size_t idx = 0; idx < occsz; ++idx) {
        if (fabs(eocc[idx]) < thresholdzero) {break;}
        pds(irw + 1,icl + 1) += CMO(irw + 1,idx + 1)*CMO(icl + 1,idx + 1)*eocc[idx];
      }
    }
  }
}
void DensVirt(matrixE & pds, matrixE & CMO, size_t occsz) {
  //function that calculates a pseudo density matrix for virtual orbitals
  size_t nrows = CMO.rows();
  size_t ncols = CMO.cols();
  if ((pds.rows() != nrows)&&(pds.cols() != ncols)) {pds.resize(nrows,ncols);}
  for (size_t irw = 0; irw < nrows; ++irw) {
    for (size_t icl = 0; icl < ncols; ++icl) {
      pds(irw + 1,icl + 1) = 0.0;
      for (size_t idx = occsz; idx < ncols; ++idx) {
        pds(irw + 1,icl + 1) += CMO(irw + 1,idx + 1)*CMO(icl + 1,idx + 1);
      }
    }
  }
}
void WriteDens(std::string density_filename, tensor<matrixE,double> & grid, matrixE & geom, std::vector<size_t> & atoms, std::vector<double> & step, double xmin, double ymin, double zmin, int prc = 7) {
  //function writing the density file
  //std::cout << "the grid" << std::endl;
  //grid.Print();
  std::string filename = density_filename + ".cub";
  std::ofstream gfile(filename, std::ios::out);
  if (!gfile.is_open()) {throw std::string("ERROR: Density.hpp: WriteDens(): could not open file for writing density");}
  double precision = pow(10.0,-prc);
  double number;
  int Natoms = geom.rows();
  int icol = 1;
  int nx = grid.dim1();
  int ny = grid.dim2();
  int nz = grid.dim3();
  gfile << std::fixed;
  gfile << std::setprecision(prc);
  gfile << "Electronic Density\nCalculated by ULYSSES\n";
  gfile << Natoms << "  " << xmin << "  " << ymin << "  " << zmin << "\n";
  gfile << nx << "  " << step[0] << "  " << 0.0 << "  " << 0.0 << "\n";
  gfile << ny << "  " << 0.0 << "  " << step[1] << "  " << 0.0 << "\n";
  gfile << nz << "  " << 0.0 << "  " << 0.0 << "  " << step[2] << "\n";
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    gfile << atoms[idAtm] << "    " << 0.0 << "  " << geom(idAtm + 1,1) << "  " << geom(idAtm + 1,2) << "  " << geom(idAtm + 1,3) << "\n";
  }
  //let the density begin
  for (size_t idx = 0; idx < nx; ++idx) {
    for (size_t idy = 0; idy < ny; ++idy) {
      for (size_t idz = 0; idz < nz; ++idz) {
        gfile << grid(idx + 1,idy + 1,idz + 1) << " ";
        if (icol == 6) {
          gfile << "\n";
          icol = 0;
        }
        ++icol;
      }
    }
  }
  gfile << "\n";
  gfile.close();
}
void ElectronicDensityCalc(std::string density_filename, matrixE & dens, matrixE & geom, std::vector<size_t> & atoms, std::vector<std::vector<double> > & AOexp, matrixE & PQN, matrixE & LQN, std::vector<double> & step, double densmatthres = 0.001, double edensthres = 1.0e-4, double gridthres = 5.0) {
  //function that calculates electronic densities for GTOs
  geom *= dist_Angstrom2aum1;
  step[0] *= dist_Angstrom2aum1;
  step[1] *= dist_Angstrom2aum1;
  step[2] *= dist_Angstrom2aum1;
  int Natoms = geom.rows();
  int NAOs = dens.rows();
  int nshellsA;
  int nshellsB;
  int tshells = AOexp[0].size();
  int icnt = 0;
  int jcnt = 0;
  int kcnt = 0;
  int lcnt = 0;
  int ishell = 0;
  int jshell = 0;
  int nsteps[3];
  int nquantumA;
  int nquantumB;
  int lquantumA;
  int lquantumB;
  int xb;
  int xe;
  int yb;
  int ye;
  int zb;
  int ze;
  bool xok;
  bool yok;
  bool zok;
  double aux;
  double aexp;
  double expthres = -log(edensthres);
  double threshexp;
  double istep[3];
  double rA[3];
  double rB[3];
  double RA;
  double RA2;
  double RB;
  double RB2;
  double radialA;
  double radialB;
  double angularA;
  double angularB;
  double dmat;
  double afactor;
  //determine the limits of the grid
  double gxmin = geom(1,1);                //limits of the grid
  double gxmax = geom(1,1);                //limits of the grid
  double gymin = geom(1,2);                //limits of the grid
  double gymax = geom(1,2);                //limits of the grid
  double gzmin = geom(1,3);                //limits of the grid
  double gzmax = geom(1,3);                //limits of the grid
  for (size_t idAtm = 1; idAtm < Natoms; ++idAtm) {
    tshells += AOexp[idAtm].size();
    if (gxmin > geom(idAtm + 1,1)) {gxmin = geom(idAtm + 1,1);}
    if (gxmax < geom(idAtm + 1,1)) {gxmax = geom(idAtm + 1,1);}
    if (gymin > geom(idAtm + 1,2)) {gymin = geom(idAtm + 1,2);}
    if (gymax < geom(idAtm + 1,2)) {gymax = geom(idAtm + 1,2);}
    if (gzmin > geom(idAtm + 1,3)) {gzmin = geom(idAtm + 1,3);}
    if (gzmax < geom(idAtm + 1,3)) {gzmax = geom(idAtm + 1,3);}
  }
  //increase the grid around the extrema
  gxmin -= gridthres;
  gxmax += gridthres;
  gymin -= gridthres;
  gymax += gridthres;
  gzmin -= gridthres;
  gzmax += gridthres;
  //invert the steps to avoid constant division
  istep[0] = 1.0/step[0];
  istep[1] = 1.0/step[1];
  istep[2] = 1.0/step[2];
  //grid size and grid
  nsteps[0] = round((gxmax - gxmin)*istep[0]) + 1;
  nsteps[1] = round((gymax - gymin)*istep[1]) + 1;
  nsteps[2] = round((gzmax - gzmin)*istep[2]) + 1;
  tensor<matrixE,double> grid(nsteps[0],nsteps[1],nsteps[2],0.0);
  matrixE limits(tshells,6);              //limits for each atom,shell combination: xbegin, xend, ybegin, yend, zbegin, zend
  //precalculate the points over the grid where atomic density is not zero
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    nshellsA = AOexp[idAtm].size();
    for (size_t ishell = 0; ishell < nshellsA; ++ishell, ++icnt) {
      aexp = AOexp[idAtm][ishell];
      threshexp = expthres/aexp;
      aux = fmax(geom(idAtm + 1,1) - threshexp,gxmin) - gxmin;
      limits(icnt + 1,1) = round(aux*istep[0]);
      aux = fmin(geom(idAtm + 1,1) + threshexp,gxmax) - gxmin;
      limits(icnt + 1,2) = round(aux*istep[0]) + 1;
      aux = fmax(geom(idAtm + 1,2) - threshexp,gymin) - gymin;
      limits(icnt + 1,3) = round(aux*istep[1]);
      aux = fmin(geom(idAtm + 1,2) + threshexp,gymax) - gymin;
      limits(icnt + 1,4) = round(aux*istep[1]) + 1;
      aux = fmax(geom(idAtm + 1,3) - threshexp,gzmin) - gzmin;
      limits(icnt + 1,5) = round(aux*istep[2]);
      aux = fmin(geom(idAtm + 1,3) + threshexp,gzmax) - gzmin;
      limits(icnt + 1,6) = round(aux*istep[2]) + 1;
    }
  }
  //geom.Print();
  //std::cout << "-----------------" << std::endl;
  //limits.Print();
  //std::cout << "-----------------" << std::endl;
  matrixE newlim = limits;
  icnt = 0;
  for (size_t idx = 0; idx < Natoms; ++idx) {
    nshellsA = AOexp[idx].size();
    for (size_t iao = 0; iao < nshellsA; ++iao, ++icnt) {
      newlim(icnt + 1,1) *= step[0];
      newlim(icnt + 1,2) *= step[0];
      newlim(icnt + 1,3) *= step[1];
      newlim(icnt + 1,4) *= step[1];
      newlim(icnt + 1,5) *= step[2];
      newlim(icnt + 1,6) *= step[2];
      newlim(icnt + 1,1) += gxmin;
      newlim(icnt + 1,2) += gxmin;
      newlim(icnt + 1,3) += gymin;
      newlim(icnt + 1,4) += gymin;
      newlim(icnt + 1,5) += gzmin;
      newlim(icnt + 1,6) += gzmin;
    }
  }
  //newlim.Print();
  //std::cout << "-----------------" << std::endl;
  icnt = 0;
  //now the electronic density
  for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
    nshellsA = AOexp[idAtm].size();
    jcnt = 0;
    jshell = 0;
    for (size_t idBtm = 0; idBtm < idAtm + 1; ++idBtm) {
      nshellsB = AOexp[idBtm].size();
      kcnt = 0;
      afactor = 1.0 + (idAtm != idBtm);
      for (size_t ishellA = 0; ishellA < nshellsA; ++ishellA) {
        nquantumA = PQN(idAtm + 1,ishellA + 1);
        lquantumA = LQN(idAtm + 1,ishellA + 1);
        lcnt = 0;
        for (size_t ishellB = 0; ishellB < nshellsB; ++ishellB) {
          //see if the atom boxes overlap on each dimension
          xok = ((limits(jshell + ishellB + 1,1) >= limits(ishell + ishellA + 1,1))&&(limits(jshell + ishellB + 1,1) <= limits(ishell + ishellA + 1,2)));
          xok += ((limits(ishell + ishellA + 1,1) >= limits(jshell + ishellB + 1,1))&&(limits(ishell + ishellA + 1,1) <= limits(jshell + ishellB + 1,2)));
          yok = ((limits(jshell + ishellB + 1,3) >= limits(ishell + ishellA + 1,3))&&(limits(jshell + ishellB + 1,3) <= limits(ishell + ishellA + 1,4)));
          yok += ((limits(ishell + ishellA + 1,3) >= limits(jshell + ishellB + 1,3))&&(limits(ishell + ishellA + 1,3) <= limits(jshell + ishellB + 1,4)));
          zok = ((limits(jshell + ishellB + 1,5) >= limits(ishell + ishellA + 1,5))&&(limits(jshell + ishellB + 1,5) <= limits(ishell + ishellA + 1,6)));
          zok += ((limits(ishell + ishellA + 1,5) >= limits(jshell + ishellB + 1,5))&&(limits(ishell + ishellA + 1,5) <= limits(jshell + ishellB + 1,6)));
          nquantumB = PQN(idBtm + 1,ishellB + 1);
          lquantumB = LQN(idBtm + 1,ishellB + 1);
          if ((xok)&&(yok)&&(zok)) {
            //calculation makes sense on the overlapping regions
            xb = int(fmax(limits(jshell + ishellB + 1,1),limits(ishell + ishellA + 1,1)));
            xe = int(fmin(limits(jshell + ishellB + 1,2),limits(ishell + ishellA + 1,2)));
            yb = int(fmax(limits(jshell + ishellB + 1,3),limits(ishell + ishellA + 1,3)));
            ye = int(fmin(limits(jshell + ishellB + 1,4),limits(ishell + ishellA + 1,4)));
            zb = int(fmax(limits(jshell + ishellB + 1,5),limits(ishell + ishellA + 1,5)));
            ze = int(fmin(limits(jshell + ishellB + 1,6),limits(ishell + ishellA + 1,6)));
            for (size_t iAO = 0; iAO < 2*lquantumA + 1; ++iAO) {
              for (size_t iBO = 0; iBO < 2*lquantumB + 1; ++iBO) {
                dmat = dens(icnt + kcnt + iAO + 1,jcnt + lcnt + iBO + 1);
                if (fabs(dmat) < densmatthres) {continue;}
                dmat *= afactor;
                //loop over grid
                for (int xstep = xb; xstep < xe; ++xstep) {
                  rA[0] = geom(idAtm + 1,1) - (gxmin + xstep*step[0]);
                  rB[0] = geom(idBtm + 1,1) - (gxmin + xstep*step[0]);
                  for (int ystep = yb; ystep < ye; ++ystep) {
                    rA[1] = geom(idAtm + 1,2) - (gymin + ystep*step[1]);
                    rB[1] = geom(idBtm + 1,2) - (gymin + ystep*step[1]);
                    for (int zstep = zb; zstep < ze; ++zstep) {
                      rA[2] = geom(idAtm + 1,3) - (gzmin + zstep*step[2]);
                      rB[2] = geom(idBtm + 1,3) - (gzmin + zstep*step[2]);
                      RA2 = rA[0]*rA[0] + rA[1]*rA[1] + rA[2]*rA[2];
                      RB2 = rB[0]*rB[0] + rB[1]*rB[1] + rB[2]*rB[2];
                      RA = sqrt(RA2);
                      RB = sqrt(RB2);
                      radialA = RSlater(RA,AOexp[idAtm][ishellA],nquantumA);
                      radialB = RSlater(RB,AOexp[idBtm][ishellB],nquantumB);
                      angularA = ASlater(lquantumA,iAO,rA[0],rA[1],rA[2],RA,RA2);
                      angularB = ASlater(lquantumB,iBO,rB[0],rB[1],rB[2],RB,RB2);
                      grid(xstep + 1,ystep + 1,zstep + 1) += radialA*angularA*dmat*angularB*radialB;
                    }
                  }
                }
              }
            }
          }
          lcnt += 2*lquantumB + 1;
        }
        kcnt += 2*lquantumA + 1;
      }
      jshell += nshellsB;
      jcnt += nshellsB*nshellsB;
    }
    ishell += nshellsA;
    icnt += nshellsA*nshellsA;
  }
  WriteDens(density_filename,grid,geom,atoms,step,gxmin,gymin,gzmin,7);
}

#endif //_Density_Matrix_
