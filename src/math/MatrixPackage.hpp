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

#ifndef _Matrix_Package_
#define _Matrix_Package_
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <math.h>
#include "eigen/Dense"
#include "eigen/Cholesky"

//description:
//wrappers around eigen matrix classes for general use

struct matrixE {
  //struct allows easier code
  size_t numberrows;
  size_t numbercols;
  Eigen::MatrixXd thematrix;
//------------------------------------------------------------------
  matrixE(size_t dimr = 1, size_t dimc = 1, double initval = 0.0) {
    numberrows = dimr;
    numbercols = dimc;
    this->thematrix.resize(numberrows,numbercols);
    if (numberrows == numbercols) {
      for (size_t idr = 0; idr < numberrows; ++idr) {
        this->thematrix(idr,idr) = initval;
        for (size_t idc = 0; idc < idr; ++idc) {
          this->thematrix(idr,idc) = initval;
          this->thematrix(idc,idr) = initval;
        }
      }
    }
    else {
      for (size_t idr = 0; idr < numberrows; ++idr) {
        for (size_t idc = 0; idc < numbercols; ++idc) {
          this->thematrix(idr,idc) = initval;
        }
      }
    }
  }
  matrixE(const std::vector<double> & vec) {
    numberrows = vec.size();
    numbercols = 1;
    this->thematrix.resize(numberrows,numbercols);
    for (size_t idr = 0; idr < numberrows; ++idr) {
      this->thematrix(idr,0) = vec[idr];
    }
  }
  ~matrixE() {}
  //getters
  size_t rows() {return numberrows;}
  size_t rows() const {return numberrows;}
  size_t cols() {return numbercols;}
  size_t cols() const {return numbercols;}
  double & operator () (size_t row, size_t col) {return this->thematrix(row - 1,col - 1);}
  double operator () (size_t row, size_t col) const {return this->thematrix(row - 1,col - 1);}
  void getBlock(matrixE & blockmat, size_t aa, size_t bb, size_t cc, size_t dd) {
    //aa -> row we start fetching from
    //bb -> column we start fetching from
    //cc -> number rows to fetch
    //dd -> number columns to fetch
    --aa;
    --bb;
    blockmat.thematrix = this->thematrix.block(aa,bb,cc,dd);
    blockmat.numberrows = cc;
    blockmat.numbercols = dd;
  }
  //other functions
  void CalculateSize() {
    //set number of rows and columns to what thematrix has
    numberrows = thematrix.rows();
    numbercols = thematrix.cols();
  }
  void clear() {
    //function to clear the matrix
    numberrows = 1;
    numbercols = 1;
    this->thematrix.resize(numberrows,numbercols);
    this->thematrix(0,0) = 0.0;
  }
  void Print(int prc = 7) {
    std::cout << std::fixed;
    std::cout << std::setprecision(prc);
    double precision = pow(10.0,-prc);
    double number;
    for (size_t idr = 0; idr < numberrows; ++idr) {
      for (size_t idc = 0; idc < numbercols; ++idc) {
        number = thematrix(idr,idc);
        if (fabs(number) < precision) {number = 0.0;}
        if (number >= -precision) {std::cout << " ";}
        std::cout << number << " ";
      }
      std::cout << "\n";
    }
    std::cout.flush();
  }
  void resize(size_t dimr, size_t dimc) {
    this->thematrix.conservativeResize(dimr,dimc);
    if ((dimr > numberrows)&&(dimc <= numbercols)) {
      for (size_t idr = numberrows; idr < dimr; ++idr) {
        for (size_t idc = 0; idc < dimc; ++idc) {
          this->thematrix(idr,idc) = 0.0;
        }
      }
    }
    else if ((dimr <= numberrows)&&(dimc > numbercols)) {
      for (size_t idr = 0; idr < dimr; ++idr) {
        for (size_t idc = numbercols; idc < dimc; ++idc) {
          this->thematrix(idr,idc) = 0.0;
        }
      }
    }
    else if ((dimr > numberrows)&&(dimc > numbercols)) {
      for (size_t idr = 0; idr < dimr; ++idr) {
        for (size_t idc = 0; idc < dimc; ++idc) {
          if ((idr >= numberrows)||(idc >= numbercols)) {
            this->thematrix(idr,idc) = 0.0;
          }
        }
      }
    }
    numberrows = dimr;
    numbercols = dimc;
  }
  matrixE operator + (const matrixE & mat2) {
    matrixE res;
    res.thematrix = this->thematrix + mat2.thematrix;
    res.numberrows = this->numberrows;
    res.numbercols = this->numbercols;
    return res;
  }
  matrixE operator - (const matrixE & mat2) {
    matrixE res;
    res.thematrix = this->thematrix - mat2.thematrix;
    res.numberrows = this->numberrows;
    res.numbercols = this->numbercols;
    return res;
  }
  void operator += (const matrixE & mat2) {this->thematrix += mat2.thematrix;}
  void operator -= (const matrixE & mat2) {this->thematrix -= mat2.thematrix;}
  void operator += (double value) {
    //consistent addition of scalar to matrix
    if (numberrows == numbercols) {
      for (size_t idr = 0; idr < numberrows; ++idr) {
        this->thematrix(idr,idr) += value;
        for (size_t idc = 0; idc < idr; ++idc) {
          this->thematrix(idr,idc) += value;
          this->thematrix(idc,idr) += value;
        }
      }
    }
    else {
      for (size_t idr = 0; idr < numberrows; ++idr) {
        for (size_t idc = 0; idc < numbercols; ++idc) {
          this->thematrix(idr,idc) += value;
        }
      }
    }
  }
  void operator -= (double value) {
    //consistent subtraction of scalar to matrix
    if (numberrows == numbercols) {
      for (size_t idr = 0; idr < numberrows; ++idr) {
        this->thematrix(idr,idr) -= value;
        for (size_t idc = 0; idc < idr; ++idc) {
          this->thematrix(idr,idc) -= value;
          this->thematrix(idc,idr) -= value;
        }
      }
    }
    else {
      for (size_t idr = 0; idr < numberrows; ++idr) {
        for (size_t idc = 0; idc < numbercols; ++idc) {
          this->thematrix(idr,idc) -= value;
        }
      }
    }
  }
  void operator = (const matrixE & mat2) {
    this->thematrix = mat2.thematrix;
    numberrows = mat2.numberrows;
    numbercols = mat2.numbercols;
  }
  void operator = (const std::vector<double> & vec) {
    numberrows = vec.size();
    numbercols = 1;
    this->thematrix.conservativeResize(numberrows,numbercols);
    for (size_t idr = 0; idr < numberrows; ++idr) {
      this->thematrix(idr,0) = vec[idr];
    }
  }
  matrixE operator * (const matrixE & mat2) {
    matrixE res(numberrows,mat2.numbercols);
    res.thematrix = this->thematrix*mat2.thematrix;
    return res;
  }
  void operator *= (const matrixE & mat2) {
    matrixE res;
    res.thematrix = this->thematrix*mat2.thematrix;
    this->thematrix = res.thematrix;
    this->numbercols = mat2.numbercols;
  }
  matrixE operator * (double value) {
    matrixE res(numberrows,numbercols);
    res.thematrix = this->thematrix*value;
    return res;
  }
  matrixE operator / (double value) {
    matrixE res(numberrows,numbercols);
    res.thematrix = this->thematrix/value;
    return res;
  }
  void operator *= (double value) {this->thematrix *= value;}
  void operator /= (double value) {this->thematrix /= value;}
  matrixE operator * (const std::vector<double> & vec) {
    if (cols() != vec.size()) {throw std::string("ERROR: MatrixPackage.hpp: matrixE: operator *() with vector: contraction dimensions not matching");}
    else {
      matrixE res(numberrows,1);
      for (size_t row = 0; row < numberrows; ++row) {
        res(row + 1,1) = 0.0;
        for (size_t contr = 0; contr < numbercols; ++contr) {
          res(row + 1,1) += this->thematrix(row,contr)*vec[contr];
        }
      }
      return res;
    }
  }
  void operator *= (const std::vector<double> & vec) {
    if (cols() != vec.size()) {throw std::string("ERROR: MatrixPackage.hpp: matrixE: operator *() with vector: contraction dimensions not matching");}
    else {
      matrixE res(numberrows,1);
      for (size_t row = 0; row < numberrows; ++row) {
        res(row + 1,1) = 0.0;
        for (size_t contr = 0; contr < numbercols; ++contr) {
          res(row + 1,1) += this->thematrix(row,contr)*vec[contr];
        }
      }
      this->thematrix = res.thematrix;
      numbercols = 1;
    }
  }
  bool operator == (const matrixE & Bmat) {
    bool equal = true;
    if (numberrows != Bmat.rows()) {equal = false;} 
    if (numbercols != Bmat.cols()) {equal = false;} 
    if (equal) {
      if (numberrows == numbercols) {
        for (size_t idr = 0; idr < numberrows; ++idr) {
          if (this->thematrix(idr,idr) != Bmat(idr + 1,idr + 1)) {
            equal = false;
            break;
          }
          for (size_t idc = 0; idc < idr; ++idc) {
            if (this->thematrix(idr,idc) != Bmat(idr + 1,idc + 1)) {equal = false;}
            else if (this->thematrix(idc,idr) != Bmat(idc + 1,idr + 1)) {equal = false;}
            if (!equal) {break;}
          }
          if (!equal) {break;}
        }
      }
      else {
        for (size_t idr = 0; idr < numberrows; ++idr) {
          for (size_t idc = 0; idc < numbercols; ++idc) {
            if (this->thematrix(idr,idc) != Bmat(idr + 1,idc + 1)) {
              equal = false;
              break;
            }
          }
          if (!equal) {break;}
        }
      }
    }
    return equal;
  }
  bool operator != (const matrixE & Bmat) {
    bool unequal = false;
    if (numberrows != Bmat.rows()) {unequal = true;} 
    if (numbercols != Bmat.cols()) {unequal = true;} 
    if (!unequal) {
      unequal = false;
      if (numberrows == numbercols) {
        for (size_t idr = 0; idr < numberrows; ++idr) {
          if (this->thematrix(idr,idr) != Bmat(idr + 1,idr + 1)) {
            unequal = true;
            break;
          }
          for (size_t idc = 0; idc < idr; ++idc) {
            if (this->thematrix(idr,idc) != Bmat(idr + 1,idc + 1)) {unequal = true;}
            else if (this->thematrix(idc,idr) != Bmat(idc + 1,idr + 1)) {unequal = true;}
            if (unequal) {break;}
          }
          if (unequal) {break;}
        }
      }
      else {
        for (size_t idr = 0; idr < numberrows; ++idr) {
          for (size_t idc = 0; idc < numbercols; ++idc) {
            if (this->thematrix(idr,idc) != Bmat(idr + 1,idc + 1)) {
              unequal = true;
              break;
            }
          }
          if (unequal) {break;}
        }
      }
    }
    return unequal;
  }
  matrixE trans() {
    //function to transpose matrix
    matrixE res(numbercols,numberrows);
    res.thematrix = this->thematrix.transpose();
    return res;
  }
  matrixE adjoint() {
    //function to get adjoint of matrix
    matrixE res;
    res.thematrix = this->thematrix.adjoint();
    res.numberrows = res.thematrix.rows();
    res.numbercols = res.thematrix.cols();
    return res;
  }
  matrixE inverse() {
    //method to invert matrices
    matrixE res;
    if (numberrows == numbercols) {
      res.thematrix = this->thematrix.inverse();
      res.numberrows = res.thematrix.rows();
      res.numbercols = res.thematrix.cols();
    }
    return res;
  }
  bool isdiagonal(double tol = 1e-7) {
    bool diag = true;
    if (numberrows != numbercols) {diag = false;}
    else {
      for (size_t idr = 0; idr < numberrows; ++idr) {
        for (size_t idc = 0; idc < idr; ++idc) {
          if (fabs(this->thematrix(idr,idc)) >= tol) {diag = false;}
          else if (fabs(this->thematrix(idc,idr)) >= tol) {diag = false;}
          if (!diag) {break;}
        }
        if (!diag) {break;}
      }
    }
    return diag;
  }
  bool isidentity(double tol = 1e-7) {
    bool ident = true;
    if (numberrows != numbercols) {ident = false;}
    else {
      for (size_t idr = 0; idr < numberrows; ++idr) {
        if (fabs(this->thematrix(idr,idr) - 1.0) >= tol) {
          ident = false;
          break;
        }
        for (size_t idc = 0; idc < idr; ++idc) {
          if (fabs(this->thematrix(idr,idc)) >= tol) {ident = false;}
          else if (fabs(this->thematrix(idc,idr)) >= tol) {ident = false;}
          if (!ident) {break;}
        }
        if (!ident) {break;}
      }
    }
    return ident;
  }
  bool issymmetric(double tol = 1e-7) {
    bool symm = true;
    if (numberrows != numbercols) {symm = false;}
    else {
      for (size_t idr = 0; idr < numberrows; ++idr) {
        for (size_t idc = 0; idc < idr; ++idc) {
          if (fabs(this->thematrix(idr,idc) - this->thematrix(idc,idr)) >= tol) {
            symm = false;
            std::cout << "elements (" << idr+1 << "," << idc+1 << ") and (" <<  idc+1 << "," << idr+1 << ") contribute for asymmetry: " << this->thematrix(idr,idc) << " " << this->thematrix(idc,idr) << "\n";
            break;
          }
        }
        if (!symm) {break;}
      }
    }
    return symm;
  }
  bool isunitary(double tol = 1e-7) {
    bool uni = (numberrows == numbercols);
    if (uni) {
      matrixE INVERSE;
      INVERSE.thematrix = this->thematrix.inverse();
      matrixE ADJOINT;
      ADJOINT.thematrix = this->thematrix.adjoint();
      for (size_t idrow = 0; idrow < numberrows; ++idrow) {
        if (fabs(INVERSE(idrow + 1,idrow + 1) - ADJOINT(idrow + 1,idrow + 1)) >= tol) {
          uni = false;
          break;
        }
        for (size_t idcol = 0; idcol < idrow; ++idcol) {
          if (fabs(INVERSE(idrow + 1,idcol + 1) - ADJOINT(idrow + 1,idcol + 1)) >= tol) {uni = false;}
          else if (fabs(INVERSE(idcol + 1,idrow + 1) - ADJOINT(idcol + 1,idrow + 1)) >= tol) {uni = false;}
          if (!uni) {break;}
        }
        if (!uni) {break;}
      }
    }
    return uni;
  }
  bool isorthogonal(double tol = 1e-7) {
    bool orthog = (numberrows == numbercols);
    if (orthog) {
      matrixE INVERSE;
      INVERSE.thematrix = this->thematrix.inverse();
      matrixE TRANSPOSE;
      TRANSPOSE.thematrix = this->thematrix.transpose();
      for (size_t idrow = 0; idrow < numberrows; ++idrow) {
        if (fabs(INVERSE(idrow + 1,idrow + 1) - TRANSPOSE(idrow + 1,idrow + 1)) >= tol) {
          orthog = false;
          break;
        }
        for (size_t idcol = 0; idcol < idrow; ++idcol) {
          if (fabs(INVERSE(idrow + 1,idcol + 1) - TRANSPOSE(idrow + 1,idcol + 1)) >= tol) {orthog = false;}
          else if (fabs(INVERSE(idcol + 1,idrow + 1) - TRANSPOSE(idcol + 1,idrow + 1)) >= tol) {orthog = false;}
          if (!orthog) {break;}
        }
        if (!orthog) {break;}
      }
    }
    return orthog;
  }
  double trace() {
    //method for trace
    if (numberrows != numbercols) {throw("ERROR: MatrixPackage.hpp: matrixE: trace(): the matrix must be square to have a trace");}
    return thematrix.trace();
  }
  double det() {
    //method for determinant
    if (numberrows != numbercols) {throw("ERROR: MatrixPackage.hpp: matrixE: det(): the matrix must be square to have a determinant");}
    return thematrix.determinant();
  }
  void zero() {
    //method to zero a matrix
    if (numberrows == numbercols) {
      for (size_t idr = 0; idr < numberrows; ++idr) {
        this->thematrix(idr,idr) = 0.0;
        for (size_t idc = 0; idc < idr; ++idc) {
          this->thematrix(idr,idc) = 0.0;
          this->thematrix(idc,idr) = 0.0;
        }
      }
    }
    else {
      for (size_t idr = 0; idr < numberrows; ++idr) {
        for (size_t idc = 0; idc < numbercols; ++idc) {
          this->thematrix(idr,idc) = 0.0;
        }
      }
    }
  }
  double Max() {
    //function returning the maximum of matrix
    return thematrix.maxCoeff();
  }
  double Min() {
    //function returning the minimum of matrix
    return thematrix.minCoeff();
  }
  double Sum() {
    //function returning the sum of all elements in a matrix
    return thematrix.sum();
  }
  double Prod() {
    //function returning the product of all elements in a matrix
    return thematrix.prod();
  }
  double Avg() {
    //function returning the average of all elements in a matrix
    return thematrix.mean();
  }
};
void Vector2DiagMat(const std::vector<double> & vec, matrixE & mat) {
  //function that converts vector into matrix
  size_t sizevec = vec.size();
  if ((mat.rows() != sizevec)||(mat.cols() != sizevec)) {mat.resize(sizevec,sizevec);}
  for (size_t idr = 0; idr < sizevec; ++idr) {
    mat(idr + 1,idr + 1) = vec[idr];
    for (size_t idc = 0; idc < idr; ++idc) {
      mat(idr + 1,idc + 1) = 0.0;
      mat(idc + 1,idr + 1) = 0.0;
    }
  }
}
matrixE operator * (const std::vector<double> & lhs, const matrixE & rhs) {
  //the struct matrix constains the implementation Au, where A is a matrix and u a column vector
  //this is the implementation of u*A, where u* is a row vector and A a matrix
  //the return type is a matrix!
  size_t dimr = rhs.rows();
  size_t dimc = rhs.cols();
  if (dimr != lhs.size()) {throw std::string("ERROR: MatrixPackage.hpp: operator *(): contraction dimensions not matching");}
  else {
    matrixE res(1,dimc);
    for (size_t icol = 0; icol < dimc; ++icol) {
      res(1,icol + 1) = 0.0;
      for (size_t contr = 0; contr < dimr; ++contr) {
        res(1,icol + 1) += lhs[contr]*rhs(contr + 1,icol + 1);
      }
    }
    return res;
  }
}
matrixE operator * (double lhs, const matrixE & rhs) {
  //the struct matrix constains the implementation Aa, where A is a matrix and a is a scalar
  //this is the implementation of aA
  //the return type is a matrix!
  matrixE res = rhs;
  res *= lhs;
  return res;
}
matrixE Givens(size_t iel, size_t jel, size_t nn, double cost, double sint) {
  if (iel == jel) {throw std::string("ERROR: MatrixPackage.hpp: Givens(): i and j cannot be the same");}
  matrixE givens(nn,nn);
  for (size_t idx = 0; idx < nn; ++idx) {
    if ((idx + 1) == iel) {givens(idx + 1,idx + 1) = cost;}
    else if ((idx + 1) == jel) {givens(idx + 1,idx + 1) = cost;}
    else {givens(idx + 1,idx + 1) = 1.0;}
  }
  givens(iel,jel) = sint;
  givens(jel,iel) = -sint;
  return givens;
}
matrixE Identity(size_t nn) {
  //creation of the identity matrix
  matrixE ident(nn,nn);
  for (size_t idx = 0; idx < nn; ++idx) {
    ident(idx + 1,idx + 1) = 1.0;
  }
  return ident;
}
void CholeskyELLT(matrixE & Amatrix) {
  //wrapper around Eigen's Cholesky
  Eigen::LLT<Eigen::MatrixXd> lltOfA(Amatrix.thematrix);
  Amatrix.thematrix = lltOfA.matrixL();
}
void CholeskyELDLT(matrixE & Amatrix) {
  //wrapper around Eigen's Cholesky LDLT
  Eigen::LDLT<Eigen::MatrixXd> ldltOfA(Amatrix.thematrix);
  Amatrix.thematrix = ldltOfA.matrixL();
}
int Cholesky(matrixE & Amatrix) {
  //funciton calculating the Cholesky decomposition of a matrix, based on the LAPACK's DPOTF2
  int information = 0;
  int nrows = Amatrix.rows();
  int ncols = Amatrix.cols();
  double ajj;
  double aux;
  //check input parameters
  if (nrows < 0) {information = -2;}
  if (ncols < 1 || ncols < nrows) {information = -4;}
  if (nrows == 0) {information = -3;}                                   //this differs from the LAPACK implementation, but I want it to be an error
  if (information >= 0) {                                               //compute the Cholesky factorization A = L*L'
    for (size_t idx = 0; idx < nrows; ++idx) {
      ajj = Amatrix(idx + 1,idx + 1);
      for (size_t idy = 0; idy < idx; ++idy) {
        ajj -= Amatrix(idx + 1,idy + 1)*Amatrix(idx + 1,idy + 1);
      }
      if (ajj <= 0.0) {                                                 //test for non-positive-definiteness
        Amatrix(idx + 1,idx + 1) = ajj;
        information = idx + 1;
        break;
      }
      Amatrix(idx + 1,idx + 1) = sqrt(ajj);
      ajj = Amatrix(idx + 1,idx + 1);
      aux = 1.0/ajj;
      //get rest of elements in row idx + 1
      if (idx + 1 < nrows) {
        for (size_t idy = idx + 1; idy < nrows; ++idy) {
          for (size_t idz = 0; idz < idx; ++idz) {
            Amatrix(idy + 1,idx + 1) -= Amatrix(idy + 1,idz + 1)*Amatrix(idx + 1,idz + 1);
          }
          Amatrix(idy + 1,idx + 1) *= aux;
        }
        for (size_t idy = idx + 1; idy < nrows; ++idy) {
          Amatrix(idx + 1,idy + 1) = 0.0;
        }
      }
    }
  }
  return information;
}
std::vector<double> CholLDLt(matrixE & Amatrix) {
  //uses Cholesky decomposition to diagonalize A and replaces it by L
  CholeskyELLT(Amatrix);                                  //A contains the matrix L
  int nrows = Amatrix.rows();
  int ncols = Amatrix.cols();
  std::vector<double> diagonal(nrows,0.0);
  double aux;
  for (size_t idxc = 0; idxc < ncols; ++idxc) {
    aux = Amatrix(idxc + 1,idxc + 1);
    diagonal[idxc] = aux*aux;
    Amatrix(idxc + 1,idxc + 1) = 1.0;
    aux = 1.0/sqrt(diagonal[idxc]);
    for (size_t idxr = idxc + 1; idxr < nrows; ++idxr) {
      Amatrix(idxr + 1,idxc + 1) *= aux;
    }
  }
  return diagonal;
}
double trace(matrixE & A,matrixE & B) {
  //function calculating the trace of a product of matrices
  double trace = 0;
  size_t nrows = A.rows();
  size_t ncols = A.cols();
  if (nrows != B.cols()) {throw("ERROR: MatrixPackage.hpp: trace(): The product matrix must be square to have a trace");}
  if (ncols != B.rows()) {throw("ERROR: MatrixPackage.hpp: trace(): Incorrect dimensions for matrix contraction");}
  for (size_t idx = 0; idx < nrows; ++idx) {
    for (size_t idy = 0; idy < ncols; ++idy) {
      trace += A(idx + 1,idy + 1)*B(idy + 1,idx + 1);
    }
  }
  return trace;
}
std::vector<double> MatDiag(matrixE & eigvec) {
  //wrapper for matrix diagonalization
  size_t dim = eigvec.rows();
  std::vector<double> eigval(dim);
  //solve eigenvalue problem
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(eigvec.thematrix);
  if (eigensolver.info() != Eigen::Success) {abort();}
  //saving eigenvectors and eigenvalues
  eigvec.thematrix = eigensolver.eigenvectors();
  Eigen::VectorXd eigenval = eigensolver.eigenvalues();
  for (size_t irw = 0; irw < dim; ++irw) {
    eigval[irw] = eigenval(irw);
  }
  return eigval;
}
void SVD(matrixE & mat4svd, matrixE & Umat, matrixE & Sval, matrixE & Vmat) {
  //wrapper for SVD of small matrices
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(mat4svd.thematrix,Eigen::ComputeThinU|Eigen::ComputeThinV);
  Umat.thematrix = svd.matrixU();
  Umat.CalculateSize();
  Sval.thematrix = svd.singularValues();
  Sval.CalculateSize();
  Vmat.thematrix = svd.matrixV();
  Vmat.CalculateSize();
}
void ConcatenateMR(matrixE & MA, const matrixE & MB) {
  //function that concatenates two matrices and stores the result in MA; concatenation over rows
  size_t MArows = MA.rows();
  size_t MBrows = MB.rows();
  size_t MAcols = MA.cols();
  if (MAcols != MB.cols()) {throw std::string("ERROR: MatrixPackage.hpp: ConcatenateMR(): matrices with different number of columns");}
  MA.resize(MArows + MBrows,MAcols);
  for (size_t idr = MArows; idr < MArows + MBrows; ++idr) {
    for (size_t idc = 0; idc < MAcols; ++idc) {
      MA(idr + 1,idc + 1) = MB(idr - MArows + 1,idc + 1);
    }
  }
}
void GramSchmidtOrthnorm(matrixE & morth, double thresholdzero = 1.0e-16) {
  //function that Gram-Schmidt orthogonalizes a matrix (a vector space) according to columns
  //linearly dependent vectors are zeroed
  int nrows = morth.rows();
  int ncols = morth.cols();
  size_t idrow;
  size_t idcaux;
  double aux1;
  double aux2;
  for (size_t idcol = 0; idcol < ncols; ++idcol) {
    aux1 = 0.0;
    for (idrow = 0; idrow < nrows; ++idrow) {
      aux1 += morth(idrow + 1,idcol + 1)*morth(idrow + 1,idcol + 1);
    }
    if (aux1 < thresholdzero) {
      //making sure the zero is zero; this might come from linear dependencies
      for (idrow = 0; idrow < nrows; ++idrow) {
        morth(idrow + 1,idcol + 1) = 0.0;
      }
      continue;
    }
    else {
      aux2 = 1.0/sqrt(aux1);
      for (idrow = 0; idrow < nrows; ++idrow) {         //normalization
        morth(idrow + 1,idcol + 1) *= aux2;
      }
      for (idcaux = idcol + 1; idcaux < ncols; ++idcaux) {
        //get overlap
        aux1 = 0.0;
        for (idrow = 0; idrow < nrows; ++idrow) {
          aux1 += morth(idrow + 1,idcol + 1)*morth(idrow + 1,idcaux + 1);
        }
        //orthogonalize
        for (idrow = 0; idrow < nrows; ++idrow) {
          morth(idrow + 1,idcaux + 1) -= aux1*morth(idrow + 1,idcol + 1);
        }
      }
    }
  }
}
//a simple matrix-like container; since I won't do much algebra with this, then no great math involved
template <class tp> class matrix {
  size_t numberrows;                       //number of rows
  size_t numbercols;                       //number of columns
  std::vector<tp> thematrix;
public:
  matrix(size_t dimr_ = 1, size_t dimc_ = 1) {
    numberrows = dimr_;
    numbercols = dimc_;
    thematrix.resize(numberrows*numbercols);
  }
  matrix(size_t dimr_, size_t dimc_, tp initval) {
    numberrows = dimr_;
    numbercols = dimc_;
    thematrix.resize(numberrows*numbercols);
    for (size_t ind = 0; ind < numberrows*numbercols; ++ind) {
      thematrix[ind] = initval;
    }
  }
  ~matrix() {}
  void clear() {
    //function to clear the matrix
    thematrix.resize(1);
  }
  size_t rows() {return numberrows;}
  size_t rows() const {return numberrows;}
  size_t cols() {return numbercols;}
  size_t cols() const {return numbercols;}
  tp & element(size_t row, size_t col) {return this->thematrix[(col - 1) + numbercols*(row - 1)];}
  tp elem(size_t row, size_t col) {return this->thematrix[(col - 1) + numbercols*(row - 1)];}
  tp & operator () (size_t row, size_t col) {return this->thematrix[(col - 1) + numbercols*(row - 1)];}
  tp operator () (size_t row, size_t col) const {return this->thematrix[(col - 1) + numbercols*(row - 1)];}
  void setValue(size_t ir, size_t ic, tp value) {this->thematrix[(ic - 1) + numbercols*(ir - 1)] = value;}
  bool evalValue(size_t ir, size_t ic, tp value) {return this->thematrix[(ic - 1) + numbercols*(ir - 1)] == value;}
  void Print(size_t prc = 7) {
    std::cout << std::fixed;
    std::cout << std::setprecision(prc);
    double precision = 1.0/pow(10.0,prc);
    double number;
    size_t counter = 0;
    for (size_t ind = 0; ind < numberrows; ++ind) {
      for (size_t ind2 = 0; ind2 < numbercols; ++ind2,++counter) {
        number = double(thematrix[counter]);
        if (fabs(number) < precision) {number = 0.0;}
        if (number >= -precision) {std::cout << " ";}
        std::cout << number << " ";
      }
      std::cout << "\n";
    }
    std::cout.flush();
  }
  void resize(size_t _newr, size_t _newc) {
    matrix auxliary(_newr,_newc);
    size_t counter = 0;
    for (size_t idr = 0; idr < numberrows; ++idr) {
      for (size_t idc = 0; idc < numbercols; ++idc,++counter) {
        if ((idr < _newr)&&(idc < _newc)) {auxliary(idr + 1,idc + 1) = thematrix[counter];}
      }
    }
    numberrows = _newr;
    numbercols = _newc;
    thematrix.resize(numberrows*numbercols);
    counter = 0;
    for (size_t idr = 0; idr < numberrows; ++idr) {
      for (size_t idc = 0; idc < numbercols; ++idc,++counter) {
        if ((idr < _newr)&&(idc < _newc)) {thematrix[counter] = auxliary(idr + 1,idc + 1);}
      }
    }
  }
  matrix operator + (matrix & mat2) {
    if (numberrows != mat2.rows()) {throw std::string("ERROR: MatrixPackage.hpp: matrix: operator +(): addition of matrices with different numbers of rows");}
    else if (numbercols != mat2.cols()) {throw std::string("ERROR: MatrixPackage.hpp: matrix: operator +(): addition of matrices with different numbers of columns");}
    else {
      matrix res(numberrows,numbercols);
      size_t counter = 0;
      for (size_t irow = 0; irow < numberrows; ++irow) {
        for (size_t icol = 0; icol < numbercols; ++icol,++counter) {
          res(irow + 1,icol + 1) = thematrix[counter] + mat2(irow + 1,icol + 1);
        }
      }
      return res;
    }
  }
  void operator += (matrix & mat2) {
    if (numberrows != mat2.rows()) {throw std::string("ERROR: MatrixPackage.hpp: matrix: operator +=(): addition of matrices with different numbers of rows");}
    else if (numbercols != mat2.cols()) {throw std::string("ERROR: MatrixPackage.hpp: matrix: operator +=(): addition of matrices with different numbers of columns");}
    else {
      size_t counter = 0;
      for (size_t rs = 0; rs < numberrows; ++rs) {
        for (size_t cs = 0; cs < numbercols; ++cs,++counter) {
          thematrix[counter] += mat2(rs + 1,cs + 1);
        }
      }
    }
  }
  matrix operator - (matrix & mat2) {
    if (numberrows != mat2.rows()) {throw std::string("ERROR: MatrixPackage.hpp: matrix: operator -(): addition of matrices with different numbers of rows");}
    else if (numbercols != mat2.cols()) {throw std::string("ERROR: MatrixPackage.hpp: matrix: operator -(): addition of matrices with different numbers of columns");}
    else {
      matrix res(numberrows,numbercols);
      size_t counter = 0;
      for (size_t irow = 0; irow < numberrows; ++irow) {
        for (size_t icol = 0; icol < numbercols; ++icol,++counter) {
          res(irow + 1,icol + 1) = thematrix[counter] - mat2(irow + 1,icol + 1);
        }
      }
      return res;
    }
  }
  void operator -= (matrix & mat2) {
    if (numberrows != mat2.rows()) {throw std::string("ERROR: MatrixPackage.hpp: matrix: operator -=(): addition of matrices with different numbers of rows");}
    else if (numbercols != mat2.cols()) {throw std::string("ERROR: MatrixPackage.hpp: matrix: operator -=(): addition of matrices with different numbers of columns");}
    else {
      size_t counter = 0;
      for (size_t rs = 0; rs < numberrows; ++rs) {
        for (size_t cs = 0; cs < numbercols; ++cs,++counter) {
          thematrix[counter] -= mat2(rs + 1,cs + 1);
        }
      }
    }
  }
  bool operator == (const matrix & Bmat) {
    bool equal = true;
    if (numberrows != Bmat.rows()) {equal = false;} 
    if (numbercols != Bmat.cols()) {equal = false;} 
    if (equal) {
      size_t counter = 0;
      for (size_t idr = 0; idr < numberrows; ++idr) {
        for (size_t idc = 0; idc < numbercols; ++idc,++counter) {
          if (thematrix[counter] != Bmat(idr + 1,idc + 1)) {
            equal = false;
            break;
          }
        }
        if (!equal) {break;}
      }
    }
    return equal;
  }
  bool operator != (const matrix & Bmat) {
    bool unequal = false;
    if (numberrows != Bmat.rows()) {unequal = true;} 
    if (numbercols != Bmat.cols()) {unequal = true;} 
    if (!unequal) {
      unequal = false;
      size_t counter = 0;
      for (size_t idr = 0; idr < numberrows; ++idr) {
        for (size_t idc = 0; idc < numbercols; ++idc,++counter) {
          if (thematrix[counter] != Bmat(idr + 1,idc + 1)) {
            unequal = true;
            break;
          }
        }
        if (unequal) {break;}
      }
    }
    return unequal;
  }
  matrix trans() {
    //function to transpose matrix
    matrix res(numbercols,numberrows);
    for (size_t rs = 0; rs < numbercols; ++rs) {
      for (size_t cs = 0; cs < numberrows; ++cs) {
        res(rs + 1,cs + 1) = thematrix[cs + numbercols*rs];
      }
    }
    return res;
  }
  bool isdiagonal(double tol = 1e-5) {
    bool diag = true;
    if (numberrows != numbercols) {diag = false;}
    else {
      size_t counter = 0;
      for (size_t idx1 = 0; idx1 < numberrows; ++idx1) {
        for (size_t idx2 = 0; idx2 < numbercols; ++idx2,++counter) {
          if (idx1 != idx2) {
            if (fabs(double(thematrix[counter])) > tol) {
              diag = false;
              break;
            }
          }
        }
      }
    }
    return diag;
  }
  bool isidentity(double tol = 1e-5) {
    bool ident = true;
    if (numberrows != numbercols) {ident = false;}
    else {
      size_t counter = 0;
      for (size_t idx1 = 0; idx1 < numberrows; ++idx1) {
        for (size_t idx2 = 0; idx2 < numbercols; ++idx2,++counter) {
          if (idx1 != idx2) {
            if (fabs(double(thematrix[counter])) > tol) {
              ident = false;
              break;
            }
          }
          else {
            if (fabs(double(thematrix[counter]) - 1.0) > tol) {
              ident = false;
              break;
            }
          }
        }
      }
    }
    return ident;
  }
  bool issymmetric(double tol = 1e-5) {
    bool symm;
    if (numberrows != numbercols) {symm = false;}
    else {
      symm = true;
      for (size_t irow = 0; irow < numberrows; ++irow) {
        for (size_t icol = 0; icol < irow; ++icol) {
          if (fabs(element(irow + 1,icol + 1) - element(icol + 1,irow + 1)) > tol) {
            symm = false;
            std::cout << "elements (" << irow + 1 << "," << icol + 1 << ") and (" <<  icol + 1 << "," << irow + 1 << ") contribute for asymmetry: " << element(irow + 1,icol + 1) << " " << element(icol + 1,irow + 1) << "\n";
            break;
          }
        }
        if (!symm) {break;}
      }
    }
    return symm;
  }
  tp trace() {
    if (numberrows != numbercols) {throw("ERROR: MatrixPackage.hpp: matrix: trace(): the matrix must be square to have a trace");}
    tp trace = 0;
    for (size_t idx = 0; idx < numberrows; ++idx) {
      trace += element(idx + 1,idx + 1);
    }
    return trace;
  }
  void zero() {
    //method to zero a matrix
    for (size_t idx = 0; idx < numberrows*numbercols; ++idx) {
      thematrix[idx] = 0;
    }
  }
};

#endif //_Matrix_Package_