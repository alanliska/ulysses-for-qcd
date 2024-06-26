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

#ifndef _Vector_Package_
#define _Vector_Package_
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include "MatrixPackage.hpp"

//description:
//additional functions for vector class

template<class T>
T operator * (const std::vector<T> & lhs, const std::vector<T> & rhs) {
  if (lhs.size() != rhs.size()) {throw std::string("ERROR: VectorPackage.hpp: operator *(): vectors of different size");}
  else {
    T res = 0;
    for (size_t idx = 0; idx < lhs.size(); ++idx) {
      res += lhs.at(idx)*rhs.at(idx);
    }
    return res;
  }
} 
template<class T>
std::vector<T> operator * (const T & lhs, const std::vector<T> & rhs) {
  std::vector<T> res(rhs.size(),0);
  for (size_t idx = 0; idx < rhs.size(); ++idx) {
    res[idx] = lhs*rhs.at(idx);
  }
  return res;
} 
template<class T>
std::vector<T> operator * (matrix<T> & lhs, std::vector<T> & rhs) {
  size_t rows = lhs.rows();
  size_t cols = lhs.cols();
  if (cols != rhs.size()) {throw std::string("ERROR: VectorPackage.hpp: operator *(): matrix cannot be contracted with vector");}
  std::vector<T> res(rows,0);
  for (size_t idr = 0; idr < rows; ++idr) {
    for (size_t idc = 0; idc < cols; ++idc) {
      res[idr] += lhs(idr + 1,idc + 1)*rhs.at(idc);
    }
  }
  return res;
} 
template<class T>
std::vector<T> operator * (const std::vector<T> & lhs, const T & rhs) {
  std::vector<T> res(lhs.size(),0);
  for (size_t idx = 0; idx < lhs.size(); ++idx) {
    res[idx] = rhs*lhs.at(idx);
  }
  return res;
} 
template<class T>
std::vector<T> operator + (const std::vector<T> & lhs, const std::vector<T> & rhs) {
  if (lhs.size() != rhs.size()) {throw std::string("ERROR: VectorPackage.hpp: operator +(): vectors of different size");}
  else {
    std::vector<T> res(rhs.size(),0);
    for (size_t idx = 0; idx < lhs.size(); ++idx) {
      res[idx] = lhs.at(idx) + rhs.at(idx);
    }
    return res;
  }
} 
template<class T>
std::vector<T> operator - (const std::vector<T> & lhs, const std::vector<T> & rhs) {
  if (lhs.size() != rhs.size()) {throw std::string("ERROR: VectorPackage.hpp: operator -(): vectors of different size");}
  else {
    std::vector<T> res(rhs.size(),0);
    for (size_t idx = 0; idx < lhs.size(); ++idx) {
      res[idx] = lhs.at(idx) - rhs.at(idx);
    }
    return res;
  }
} 
template<class T>
T norm(const std::vector<T> & vec) {
  T res = 0;
  for (size_t idx = 0; idx < vec.size(); ++idx) {
    res += vec[idx]*vec[idx];
  }
  return sqrt(res);
} 
template<class T>
void ConcatenateV(std::vector<T> & VA,std::vector<T> & VB) {
  //function that concatenates two vectors and stores the result in VA
  size_t VAsize = VA.size();
  size_t VBsize = VB.size();
  VA.resize(VAsize + VBsize);
  for (size_t idv = VAsize; idv < VAsize + VBsize; ++idv) {
    VA.at(idv) = VB.at(idv - VAsize);
  }
}
double dotProd(std::vector<double> * VA,std::vector<double> * VB) {
  //function that calculates the dot product between two vectors
  if (VA->size() != VB->size()) {throw std::string("ERROR: VectorPackage.hpp: dotProd(): vectors of different size.");}
  double result = 0.0;
  for (size_t idv = 0; idv < VA->size(); ++idv) {
    result += VA->at(idv)*VB->at(idv);
  }
  return result;
}
double dotProd(const std::vector<double> & VA,const std::vector<double> & VB) {
  //function that calculates the dot product between two vectors
  if (VA.size() != VB.size()) {throw std::string("ERROR: VectorPackage.hpp: dotProd(): vectors of different size.");}
  double result = 0.0;
  for (size_t idv = 0; idv < VA.size(); ++idv) {
    result += VA.at(idv)*VB.at(idv);
  }
  return result;
}
std::vector<double> crossProd(std::vector<double> VA,std::vector<double> VB, double sint) {
  //function calculating the cross product for 3-dimensional vectors
  std::vector<double> result;
  result.push_back((VA[1]*VB[2] - VA[2]*VB[1])/sint);
  result.push_back((VA[2]*VB[0] - VA[0]*VB[2])/sint);
  result.push_back((VA[0]*VB[1] - VA[1]*VB[0])/sint);
  return result;
}
void CrossProd(std::vector<double> & v1, std::vector<double> & v2, matrixE & cross, int pos = 1, int matpos = 0) {
  //function calculating the cross product (v1 x v2) between two vectors and storing it in a matrix
  if (matpos == 0) {
    cross(pos,1) = v1[1]*v2[2] - v1[2]*v2[1];
    cross(pos,2) = v1[2]*v2[0] - v1[0]*v2[2];
    cross(pos,3) = v1[0]*v2[1] - v1[1]*v2[0];
  }
  else if (matpos > 0) {
    cross(pos,1) = v1[1]*cross(matpos,3) - v1[2]*cross(matpos,2);
    cross(pos,2) = v1[2]*cross(matpos,1) - v1[0]*cross(matpos,3);
    cross(pos,3) = v1[0]*cross(matpos,2) - v1[1]*cross(matpos,1);
  }
  else {
    matpos *= -1;
    cross(pos,1) = cross(matpos,2)*v2[2] - cross(matpos,3)*v2[1];
    cross(pos,2) = cross(matpos,3)*v2[0] - cross(matpos,1)*v2[2];
    cross(pos,3) = cross(matpos,1)*v2[1] - cross(matpos,2)*v2[0];
  }
}
void CrossProdM(matrixE & cross, int posA = 1, int posB = 2, int posC = 3) {
  //function calculating the cross product C = (A x B) between two vectors and storing it in a matrix; all vectors are in the matrix
  cross(posC,1) = cross(posA,2)*cross(posB,3) - cross(posA,3)*cross(posB,2);
  cross(posC,2) = cross(posA,3)*cross(posB,1) - cross(posA,1)*cross(posB,3);
  cross(posC,3) = cross(posA,1)*cross(posB,2) - cross(posA,2)*cross(posB,1);
}

#endif //_Vector_Package_
