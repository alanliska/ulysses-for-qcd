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

#ifndef _Tensor_Package_
#define _Tensor_Package_
#include "MatrixPackage.hpp"

//description:
//declaration of tensor classe as containers

template <class matrixtype, class valtype> class tensor {
  //matrixtype and valtype must be consistent
  std::vector<size_t> dim;
  std::vector<matrixtype> tnsr;
public:
  tensor(size_t dim1 = 1, size_t dim2 = 1, size_t dim3 = 1) {
    dim.resize(3);
    dim[0] = dim1;
    dim[1] = dim2;
    dim[2] = dim3;
    matrixtype mat(dim2,dim3);
    for (size_t ind = 0; ind < dim[0]; ++ind) {
      tnsr.push_back(mat);
    }
  }
  tensor(size_t dim1, size_t dim2, size_t dim3, valtype initval) {
    dim.resize(3);
    dim[0] = dim1;
    dim[1] = dim2;
    dim[2] = dim3;
    matrixtype mat(dim2,dim3,initval);
    for (size_t ind = 0; ind < dim[0]; ++ind) {
      tnsr.push_back(mat);
    }
  }
  ~tensor() {}
  void Print(size_t prc = 7) {
    for (size_t ind = 0; ind < dim[0]; ++ind) {
      std::cout << "dimension1: " << ind + 1 << std::endl;
      tnsr[ind].Print(prc);
    }
  }
  void resize(size_t d1, size_t d2, size_t d3) {
    size_t old0 = dim[0];
    size_t old1 = dim[1];
    size_t old2 = dim[2];
    dim[0] = d1;
    dim[1] = d2;
    dim[2] = d3;
    if (d1 >= old0) {
      for (size_t idx = 0; idx < old0; ++idx) {
        tnsr[idx].resize(d2,d3);
      }
      matrixtype mat(d2,d3);
      for (size_t idx = old0; idx < d1; ++idx) {
        tnsr.push_back(mat);
      }
    }
    else {
      tnsr.erase(tnsr.begin() + d1, tnsr.end());
      for (size_t idx = 0; idx < d1; ++idx) {
        tnsr[idx].resize(d2,d3);
      }
    }
  }
  void resize(const std::vector<int> & di) {resize(di[0],di[1],di[2]);}
  size_t dim1() {return dim[0];}
  size_t dim1() const {return dim[0];}
  size_t dim2() {return dim[1];}
  size_t dim2() const {return dim[1];}
  size_t dim3() {return dim[2];}
  size_t dim3() const {return dim[2];}
  size_t Dim(size_t ii) {return dim[ii - 1];}
  size_t Dim(size_t ii) const {return dim[ii - 1];}
  matrixtype & Block(size_t ii) {return tnsr[ii - 1];}
  valtype & element(size_t i1, size_t i2, size_t i3) {return tnsr[i1 - 1](i2,i3);}
  valtype & element(std::vector<int> i0) {return tnsr[i0[0] - 1](i0[1],i0[2]);}
  valtype elem(size_t i1, size_t i2, size_t i3) const {return tnsr[i1 - 1](i2,i3);}
  valtype & operator () (size_t i1, size_t i2, size_t i3) {return tnsr[i1 - 1](i2,i3);}
  valtype operator () (size_t i1, size_t i2, size_t i3) const {return tnsr[i1 - 1](i2,i3);}
  tensor operator + (tensor & tens2) {
    if (dim[0] != tens2.dim1()) {throw std::string("ERROR: TensorPackage.hpp: tensor: operator +(): addition of tensors with different dim1");}
    else if (dim[1] != tens2.dim2()) {throw std::string("ERROR: TensorPackage.hpp: tensor: operator +(): addition of tensors with different dim2");}
    else if (dim[2] != tens2.dim3()) {throw std::string("ERROR: TensorPackage.hpp: tensor: operator +(): addition of tensors with different dim3");}
    else {
      tensor res(dim[0],dim[1],dim[2]);
      for (size_t i1 = 0; i1 < dim[0]; ++i1) {
        res.Block(i1 + 1) = Block(i1 + 1) + tens2.Block(i1 + 1);
      }
      return res;
    }
  }
  void operator += (tensor & tens2) {
    if (dim[0] != tens2.dim1()) {throw std::string("ERROR: TensorPackage.hpp: tensor: operator +=(): addition of tensors with different dim1");}
    else if (dim[1] != tens2.dim2()) {throw std::string("ERROR: TensorPackage.hpp: tensor: operator +=(): addition of tensors with different dim2");}
    else if (dim[2] != tens2.dim3()) {throw std::string("ERROR: TensorPackage.hpp: tensor: operator +=(): addition of tensors with different dim3");}
    else {
      for (size_t i1 = 0; i1 < dim[0]; ++i1) {
        Block(i1 + 1) += tens2.Block(i1 + 1);
      }
    }
  }
  tensor operator - (tensor & tens2) {
    if (dim[0] != tens2.dim1()) {throw std::string("ERROR: TensorPackage.hpp: tensor: operator -(): addition of tensors with different dim1");}
    else if (dim[1] != tens2.dim2()) {throw std::string("ERROR: TensorPackage.hpp: tensor: operator -(): addition of tensors with different dim2");}
    else if (dim[2] != tens2.dim3()) {throw std::string("ERROR: TensorPackage.hpp: tensor: operator -(): addition of tensors with different dim3");}
    else {
      tensor res(dim[0],dim[1],dim[2]);
      for (size_t i1 = 0; i1 < dim[0]; ++i1) {
        res.Block(i1 + 1) = Block(i1 + 1) - tens2.Block(i1 + 1);
      }
      return res;
    }
  }
  void operator -= (tensor & tens2) {
    if (dim[0] != tens2.dim1()) {throw std::string("ERROR: TensorPackage.hpp: tensor: operator -=(): addition of tensors with different dim1");}
    else if (dim[1] != tens2.dim2()) {throw std::string("ERROR: TensorPackage.hpp: tensor: operator -=(): addition of tensors with different dim2");}
    else if (dim[2] != tens2.dim3()) {throw std::string("ERROR: TensorPackage.hpp: tensor: operator -=(): addition of tensors with different dim3");}
    else {
      for (size_t i1 = 0; i1 < dim[0]; ++i1) {
        Block(i1 + 1) -= tens2.Block(i1 + 1);
      }
    }
  }
  bool operator == (tensor & B) {
    bool equal = true;
    if (dim[0] != B.dim1()) {equal = false;}
    if (dim[1] != B.dim2()) {equal = false;}
    if (dim[2] != B.dim3()) {equal = false;}
    if (equal) {
      for (size_t idx = 0; idx < dim[0]; ++idx) {
        equal = (tnsr[idx] == B.Block(idx + 1));
        if (!equal) {break;}
      }
    }
    return equal;
  }
  bool operator != (tensor & B) {
    bool unequal = false;
    if (dim[0] != B.dim1()) {unequal = true;}
    if (dim[1] != B.dim2()) {unequal = true;}
    if (dim[2] != B.dim3()) {unequal = true;}
    if (!unequal) {
      unequal = true;
      for (size_t idx = 0; idx < dim[0]; ++idx) {
        unequal = (tnsr[idx] != B.Block(idx + 1));
        if (!unequal) {break;}
      }
    }
    return unequal;
  }
  void zero() {
    //method to zero a tensor
    for (size_t id1 = 0; id1 < dim[0]; ++id1) {
      tnsr[id1].zero();
    }
  }
};

#endif //_Tensor_Package_
