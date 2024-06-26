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

#ifndef _Solver_Package_
#define _Solver_Package_
#include <math.h>
#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "MatrixPackage.hpp"
#include "../Hessian.hpp"
#include "FunctionPackage.hpp"

//description:
//solvers and thresholds for geometry optimization
void xTBthresholds(std::string accuracy, double & enthresh, double & ghthresh) {
  //function returning the xTB type of thresholds for geometry optimization
  if (accuracy == "crude") {
    enthresh = 5.0e-4;          //Eh
    ghthresh = 1.0e-2;          //Eh/Bohr
  }
  else if (accuracy == "sloppy") {
    enthresh = 1.0e-4;          //Eh
    ghthresh = 6.0e-3;          //Eh/Bohr
  }
  else if (accuracy == "loose") {
    enthresh = 5.0e-5;          //Eh
    ghthresh = 4.0e-3;          //Eh/Bohr
  }
  else if (accuracy == "lax") {
    enthresh = 2.0e-5;          //Eh
    ghthresh = 2.0e-3;          //Eh/Bohr
  }
  else if (accuracy == "normal") {
    enthresh = 5.0e-6;          //Eh
    ghthresh = 1.0e-3;          //Eh/Bohr
  }
  else if (accuracy == "tight") {
    enthresh = 1.0e-6;          //Eh
    ghthresh = 8.0e-4;          //Eh/Bohr
  }
  else if (accuracy == "vtight") {
    enthresh = 1.0e-7;          //Eh
    ghthresh = 2.0e-4;          //Eh/Bohr
  }
  else if (accuracy == "extreme") {
    enthresh = 5.0e-8;          //Eh
    ghthresh = 5.0e-5;          //Eh/Bohr
  }
  ghthresh *= dist_Angstrom2aum1;
}

void Solve_Ax_eq_b(matrixE & Amat, std::vector<double> & bvec) {
  //interface to eigen to solve the system of equations Amat.x = bvec where Amat is a matrix, x and bvec vectors
  //write solution to bvec
  int ndim = bvec.size();
  int ncols = Amat.cols();
  if (Amat.rows() != ndim) {throw("ERROR: SolverPackage.hpp: Solve_Ax_eq_b: matrix A and vector b do not have consistent dimensions");}
  Eigen::VectorXd b(ndim);
  for (size_t idr = 0; idr < ndim; ++idr) {
    b(idr) = bvec[idr];
  }
  Eigen::VectorXd x = Amat.thematrix.colPivHouseholderQr().solve(b);
  for (size_t idx = 0; idx < ndim; ++idx) {
    bvec[idx] = x(idx);
  }
}
void Solve_Ax_eq_b(matrixE & Amat, matrixE & bvec) {
  //interface to eigen to solve the system of equations Amat.x = bvec where Amat is a matrix, x and bvec vectors
  //write solution to bvec
  int ndim = bvec.rows();
  int ncols = Amat.cols();
  if (Amat.rows() != ndim) {throw("ERROR: SolverPackage.hpp: Solve_Ax_eq_b: matrix A and vector b do not have consistent dimensions");}
  Eigen::VectorXd b(ndim);
  for (size_t idr = 0; idr < ndim; ++idr) {
    b(idr) = bvec(idr + 1,1);
  }
  Eigen::VectorXd x = Amat.thematrix.colPivHouseholderQr().solve(b);
  for (size_t idx = 0; idx < ndim; ++idx) {
    bvec(idx + 1,1) = x(idx);
  }
}
//solvers for polynomials
void SolveP3(std::vector<double> & xp3, double aa, double bb, double cc, int & nroots, double threshold = 1.0e-14) { 
  //solver for cubic equation x^3 + aa*x^2 + bb*x + cc = 0
  //taken and adapted from Khashin S.I. http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
  //nroots is the number of real roots; if nroots = 1, then the last two elements contain the complex number xp3[1] +- xp3[2]i
  if (xp3.size() != 3) {xp3.resize(3);}
  double a2 = aa*aa;
  double qq = (a2 - 3.0*bb)/9.0; 
  double rr = (aa*(2.0*a2 - 9.0*bb) + 27.0*cc)/54.0;
  //equation x^3 + qq*x + rr = 0
  double r2 = rr*rr;
  double q3 = qq*qq*qq;
  if ((r2 <= (q3 + threshold))&&(fabs(q3) > threshold)) {
    double tt = rr/sqrt(q3);
    if(tt < -1.0) {tt = -1.0;}
    if(tt > 1.0) {tt = 1.0;}
    tt = acos(tt);
    aa /= 3.0;
    qq = -2.0*sqrt(qq);
    xp3[0] = qq*cos(tt/3.0) - aa;
    xp3[1] = qq*cos((tt + 2.0*pi)/3.0) - aa;
    xp3[2] = qq*cos((tt - 2.0*pi)/3.0) - aa;
    nroots = 3;
  } 
  else {
    double AA = -cbrt(fabs(rr) + sqrt(r2 - q3)); 
    if (rr < 0.0) {AA = -AA;}
    double BB = 0.0;
    BB = (AA == 0? 0 : BB = qq/AA);
    aa /= 3.0;
    xp3[0] = (AA + BB) - aa;
    xp3[1] = -0.5*(AA + BB) - aa;
    xp3[2] = 0.5*sqrt(3.0)*(AA - BB);
    if (fabs(xp3[2]) < threshold) {
      xp3[2] = xp3[1];
      nroots = 2;
    }
    else {nroots = 1;}
  }
}
//minimization tools for 1 dimensional problems
double GoldenSection(double (*func)(double), double x1, double x2, bool minimize = true, double accuracy = 1.0e-6, int maxiter = 500) {
  //function applying the golden section method to determine extrema of a function
  double xmin = fmin(x1,x2);
  double xmax = fmax(x1,x2);
  double xmiddle;
  double pmiddle;
  double mmiddle;
  double fxmax;
  double fxmin;
  double fpm;
  double fmm;
  double factor = 1.0;
  double epsilon;
  if (!minimize) {factor = -1.0;}
  for (size_t iter = 0; iter < maxiter; ++iter) {
    epsilon = 0.01*(xmax - xmin);                     //define epsilon as being 1% of interval
    xmiddle = 0.5*(xmax + xmin);
    pmiddle = xmiddle + epsilon;
    mmiddle = xmiddle - epsilon;
    fxmax = factor*(*func)(xmax);
    fxmin = factor*(*func)(xmin);
    fpm = factor*(*func)(pmiddle);
    fmm = factor*(*func)(mmiddle);
    if (fpm > fmm) {xmax = pmiddle;}
    else if (fmm > fpm) {xmin = mmiddle;}
    else {
      xmax = pmiddle;
      xmin = mmiddle;
    }
    //test for convergence
    if ((xmax - xmin) < accuracy) {break;}
  }
  if ((xmax - xmin) > accuracy) {std::cout << "WARNING: SolverPackage.hpp: GoldenSection(): no convergence in " << maxiter << " iterations" << std::endl;}
  return 0.5*(xmax + xmin);
}
double NewtonRaphson1D(double (*func)(double), double (*der)(double), double x0, double accuracy = 1.0e-6, int maxiter = 500) {
  //function applying the Newton-Raphson method to find zeroes of a function func
  //note that if the function func is the first derivative and if der is the second derivative, then this function finds extrema
  double xnew = x0;
  double xold = x0;
  for (size_t iter = 0; iter < maxiter; ++iter) {
    xold = xnew;
    xnew -= (*func)(xold)/(*der)(xold);
    if ((fabs(xnew - xold) < accuracy)||(fabs((*func)(xnew)) < accuracy)) {break;}
  }
  if ((fabs(xnew - xold) > accuracy)&&(fabs((*func)(xnew)) > accuracy)) {std::cout << "WARNING: SolverPackage.hpp: NewtonRaphson1D(): no convergence in " << maxiter << " iterations" << std::endl;}
  return xnew;
}
double CubicInterpMin(double aa, double bb, double faa, double fbb, double gaa, double gbb, int & control) {
  //estimation of function's minimum by interpolation of function to cubic of the form 
  //f(x) = {c}_{3}{x}^{3} + {c}_{2}{x}^{2} + {c}_{1}{x} + {c}_{0} in the interval ]a,b] using derivatives and function values
  //based on D. P. Bersekas, Non-linear Programming, second ed. Athena Scientific, Belmont, 1999
  //aa <-> faa, gaa => function and first-derivative at point a
  //bb <-> fbb, gbb => function and first-derivative at point b
  double zz = 3.0*(faa - fbb)/(bb - aa) + gaa + gbb;
  double gamma = fmax(fabs(zz),fabs(gaa));
  double sval = fmax(gamma,fabs(gbb));
  gamma = (zz*zz - gaa*gbb)/(sval*sval);
  sval *= sqrt(fmax(0.0,gamma));
  if (aa > bb) {sval *= -1.0;}
  control = 0;
  if (((sval + zz - gaa)/(gbb - gaa + 2.0*sval) < 0.0)&&(sval != 0.0)) {control = 1;}
  return aa + (bb - aa)*(sval + zz - gaa)/(gbb - gaa + 2.0*sval);
}
double QuadrInterpMin(double aa, double bb, double faa, double fbb, double gaa) {
  //estimation of function's minimum by interpolation of function to parabola of the form 
  //f(x) = {c}_{2}{x}^{2} + {c}_{1}{x} + {c}_{0} in the interval ]a,b] using derivatives and function values
  //based on R. Fletcher, Practical Methods of Optimization, vol. 1, Unconstrained Optimization, Wiley, Chichester, 1980
  //aa <-> faa, gaa => function and first-derivative at point a
  //bb <-> fbb, gbb => function and first-derivative at point b
  return aa + 0.5*(bb - aa)*(bb - aa)*gaa/(faa - fbb + (bb - aa)*gaa);
}
double QuadrInterpGMin(double aa, double bb, double gaa, double gbb) {
  //estimation of function's minimum by interpolation of function to parabola of the form 
  //f(x) = {c}_{2}{x}^{2} + {c}_{1}{x} + {c}_{0} in the interval ]a,b] using derivatives
  //based on R. Fletcher, Practical Methods of Optimization, vol. 1, Unconstrained Optimization, Wiley, Chichester, 1980
  //aa <-> gaa => first-derivative at point a
  //bb <-> gbb => first-derivative at point b
  return bb + (bb - aa)*gbb/(gaa - gbb);
}
//minimization tools for N dimensional problems
class solver {
  //basic structure for allocation of solver
  int linesearch;             //type of line search used: 0 -> no line search
                              //                          1 -> Davidon
                              //                          2 -> Barzilai-Borwein
                              //                          3 -> Armijo
                              //                          4 -> More'-Thuente
                              //                          5 -> Fletcher's parabollic interpolation
                              //                          6 -> dogleg trust region
                              //                          7 -> trust region details for such methods
                              //                          9 -> Debug
  int bad_alphas;             //variable to count the number of bad alphas
  int last_iter;              //we want consecutive iterations with bad alphas
  int max_bad_alphas;         //variable to determine the maximum number of bad alphas
  int ijob;                   //control variable for trust-region algorithms
  double alpha_min;           //limit the minimum value for step regulator
  double alpha_max;           //limit the maximum value for step regulator
  double alphathreshzero;     //the zero for alphas
  double alpha_old;           //alpha of previous iteration
public:
  double thresholdzero;       //numerical zero
  double p5;
  double p66;
  double xtrapf;              //some stuff needed for More' and Thuente line search
  double ftol;                //this is {\mu} in the sufficient decrease strong Wolfe condition
  double gtol;                //this is {\eta} in the curvature strong Wolfe condition
  //auxiliaries to avoid constant allocation and deallocation
  double aux1;
  double aux2;
  double aux3;
  double aux4;
  double aux5;
  double aux6;
  double aux7;
  double aux8;
  double aux9;
  double ww;
  double zz;
  double stmin;
  double stmax;
  double dg;
  double ftest1;
  double fm;
  double fxm;
  double fym;
  double dgm;
  double dgxm;
  double dgym;
  double sgnd;
  double stpc;
  double stpq;
  double stpf;
  int aux;
  int xau;
  int uax;
  int axu;
  bool bracket;
  bool stage1;
  bool bound;
  matrixE Hhat;
  matrixE Hprm;
  matrixE xt;
  matrixE xn;
  matrixE cp;
  matrixE pstore;              //container for permanent storage
  solver(int _linesearch, double _alpha_min, double _alpha_max, int _max_bad_alphas = 5, double _thresholdzero = 1.0e-6, double _alphathreshzero = 1.0e-3) {
    thresholdzero = _thresholdzero;
    bad_alphas = 0;
    last_iter = -1;
    linesearch = _linesearch;
    alpha_min = _alpha_min;
    alpha_max = _alpha_max;
    if (alpha_min < thresholdzero) {alpha_min = thresholdzero;}
    if (alpha_max < thresholdzero) {throw("ERROR: SolverPackage.hpp: solver: solver(): alpha_max is below zero");}
    else if (alpha_max < alpha_min) {throw("ERROR: SolverPackage.hpp: solver: solver(): alpha_max < alpha_min");}
    max_bad_alphas = _max_bad_alphas;
    alphathreshzero = _alphathreshzero;
    if (linesearch == 4) {                         //More' and Thuente stuff
      p5 = 0.5;
      p66 = 0.66;
      xtrapf = 4.0;
      ftol = 0.001;
      gtol = 0.1;
    }
    else if (linesearch == 5) {                    //Fletcher stuff
      ftol = 0.001;
      gtol = 0.1;
      p5 = 0.05;                                   //tau in Fletcher's algorithm
      p66 = 9.0;                                   //chi in Fletcher's algorithm
    }
    ijob = 1;
  }
  ~solver() {}
  //getters
  virtual std::string Mode() {return "0";}
  int Degree() {return 0;}                                                                           //maximum order of derivatives needed
  int LineSearch() {return linesearch;}
  int MaxBadAlphas() {return max_bad_alphas;}
  int BadAlphas() {return bad_alphas;}
  double AlphaMin() {return alpha_min;}
  double AlphaMax() {return alpha_max;}
  double ThreshZero() {return thresholdzero;}
  double alphaThreshZero() {return alphathreshzero;}
  double fTolerance() {return ftol;}
  double gTolerance() {return gtol;}
  virtual std::string Type() {return "0";}
  virtual int method() {return 0;}
  virtual int Newton() {return 0;}                                                                         //some methods require a variable telling whether to use a Newton update
  //setters
  void setLineSearch(int _linesearch) {linesearch = _linesearch;}
  void setAlphaMin(double _alpha_min) {
    alpha_min = _alpha_min;
    if (alpha_min < thresholdzero) {alpha_min = thresholdzero;}
    else if (alpha_max < alpha_min) {throw("ERROR: SolverPackage.hpp: solver: solver(): alpha_max < alpha_min");}
  }
  void setAlphaMax(double _alpha_max) {
    alpha_max = _alpha_max;
    if (alpha_max < alpha_min) {throw("ERROR: SolverPackage.hpp: solver: solver(): alpha_max < alpha_min");}
  }
  void setAlphas(double _alpha_1,double _alpha_2) {
    alpha_max = fmax(_alpha_1,_alpha_2);
    alpha_min = fmin(_alpha_1,_alpha_2);
  }
  void setThreshZero(double _thresholdzero) {thresholdzero = _thresholdzero;}
  void setMaxBadAlphas(int _max_bad_alphas) {max_bad_alphas = _max_bad_alphas;}
  void setalphaThreshZero(double _alphathreshzero) {alphathreshzero = _alphathreshzero;}
  void setfTolerance(double _ftol) {ftol = _ftol;}
  void setgTolerance(double _gtol) {gtol = _gtol;}
  //calculation
  virtual void Step() {}                                                                                   //function that does the step calculation
  virtual void HessianUpdate() {}                                                                          //function that does the step calculation
  //alpha control
  void UpdateAlphaControl(double _alpha, int iteration) {
    xau = 0;
    if (bad_alphas > 0) {
      if ((fabs(_alpha) < alphathreshzero)&&(iteration - 1 == last_iter)) {xau = 1;}
    }
    else {
      if (fabs(_alpha) < alphathreshzero) {xau = 1;}
    }
    if (xau == 1) {                                            //if bad step, then add to counter; add to iteration counter as well
      ++bad_alphas;
      last_iter = iteration;
    }
    else if (xau == 0) {this->ResetAlphaControl();}            //reset counter, not consecutive
  }
  void ResetAlphaControl() {bad_alphas = 0;}                                                               //if corrective action is taken, then reset counter
  bool SteepestDescentStep() {                                                                             //determine whether corrective steepest descent step is needed
    bool res = false;
    if (bad_alphas > max_bad_alphas - 1) {res = true;}
    return res;
  }
  //line search functions
  template<class func>
  double Davidon(func & ObjFunc, matrixE & sk, matrixE & gx, matrixE & Hessian, double fnew, bool project, bool & skip) {
    //W. C. Davidon, Comp. J., 10(4), 406, 1968
    skip = false;
    xau = sk.rows();
    aux2 = 0.0;
    double stepnorm = 0.0;
    for (size_t idcoord = 0; idcoord < xau; ++idcoord, ++axu) {
      aux2 += gx(idcoord + 1,1)*sk(idcoord + 1,1);
      stepnorm += sk(idcoord + 1,1)*sk(idcoord + 1,1);
    }
    aux4 = 1.0;
    if (stepnorm > 1.0) {aux4 = 1.0/sqrt(stepnorm);}
    Hhat = aux4*sk;
    cp = ObjFunc.Geometry();                       //original geometry
    for (size_t itermicro = 0; itermicro < 100; ++itermicro) {
      ObjFunc.updGeometry(Hhat,project);
      ObjFunc.Calculate(0);
      aux3 = ObjFunc.getEnergy(1);
      if (aux3 < fnew) {break;}
      aux4 *= 0.5;
      Hhat = aux4*sk;
      ObjFunc.setGeometry(cp);
    }
    ObjFunc.gEnergy(Hprm,1,project);
    aux5 = 0.0;
    aux6 = 0.0;
    ww = 0.0;
    for (size_t idcoord = 0; idcoord < xau; ++idcoord) {
      aux5 += Hprm(idcoord + 1,1)*sk(idcoord + 1,1);
      aux6 += Hprm(idcoord + 1,1)*Hprm(idcoord + 1,1);
      ww += gx(idcoord + 1,1)*gx(idcoord + 1,1);
    }
    if (sqrt(aux6) < thresholdzero) {aux1 = aux4;}
    else {
      zz = 3.0*(fnew - aux3)/aux4 + aux2 + aux5;
      if (zz*zz < aux2*aux5) {                     //cannot really take decision, then do steepest descent on alpha
        if (aux6 < ww) {aux1 = 1.0;}
        else {aux1 = 0.0001;}
      }
      else {
        ww = sqrt(zz*zz - aux2*aux5);
        aux1 = fmin(aux4,aux4*(1.0 - (aux5 + ww - zz)/(aux5 - aux2 - 2.0*ww)));
      }
    }
    ObjFunc.setGeometry(cp);
    return aux1;
  }
  double BarzilaiBorwein(matrixE & step, matrixE & gx, matrixE & gy, size_t iter, double alpha) {
    xau = step.rows();
    aux1 = 0.0;
    aux2 = 0.0;
    aux3 = 0.0;
    for (size_t idcoord = 1; idcoord < xau + 1; ++idcoord) {
      aux1 += alpha*alpha*step(idcoord,1)*step(idcoord,1);
      aux5 = gx(idcoord,1) - gy(idcoord,1);
      aux3 += alpha*step(idcoord,1)*aux5;
      aux2 += aux5*aux5;
    }
    if (aux1 < thresholdzero) {aux1 = 1.0;}
    if (aux3 < thresholdzero) {aux3 = 1.0;}
    if (aux2 < thresholdzero) {aux2 = 1.0;}
    if ((iter == 0)&&(aux2/aux1 < 10.0)) {aux4 = aux3/aux2;}               //if gradient is too big, then do not do this step
    else {aux4 = aux1/aux3;}
    if (aux4 > 1.0) {aux4 = 1.0;}
    return aux4;
  }
  template<class func>
  double Armijo(func & ObjFunc, matrixE & step, matrixE & gnew, double fold, bool project, bool & skip, double sigma = 0.5, double beta = 0.75, int maxiter = 500) {
    skip = false;
    xau = step.rows();
    aux2 = 0.0;
    for (size_t idcoord = 1; idcoord < xau + 1; ++idcoord) {
      aux2 += gnew(idcoord,1)*step(idcoord,1);
    }
    Hhat = ObjFunc.Geometry();
    aux5 = 1.0;                         //alpha initialization
    uax = 0;
    ww = 0.0;
    zz = 0.0;
    stpf = 0.0;
    for (size_t idx = 0; idx < maxiter; ++idx, ++uax) {
      //set test value for new xk
      Hprm = step*aux5;
      ObjFunc.updGeometry(Hprm,project);
      ObjFunc.Calculate(0);
      aux1 = ObjFunc.getEnergy(1);
      ObjFunc.setGeometry(Hhat);
      stpf = aux5;
      if (aux1 <= fold - sigma*aux5*aux2) {break;}
      if (fabs(aux5 - 1.0) < thresholdzero) {aux4 = -aux2/(2.0*(aux1 - fold - aux2));}
      else {
        aux3 = aux1 - fold - aux5*aux2;
        aux6 = ww - fold - zz*aux2;
        aux7 = (aux3/(aux5*aux5) - aux6/(zz*zz))/(aux5 - zz);
        aux8 = (-zz*aux3/(aux5*aux5) + aux5*aux6/(zz*zz))/(aux5 - zz);
        if (aux7 == 0.0) {aux4 = -aux2/(2.0*aux8);}
        else {
          aux9 = aux8*aux8 - 3.0*aux7*aux2;
          if (aux9 < 0.0) {aux4 = 0.5*aux5;}
          else if (aux8 <= 0.0) {aux4 = (-aux8 + sqrt(aux9))/(3.0*aux7);}
          else {aux4 = -aux2/(aux8 + sqrt(aux9));}
        }
        if (aux4 > 0.5*aux5) {aux4 = 0.5*aux5;}                   //lambda < 0.5*lambda_1
      }
      zz = aux5;
      ww = aux1;
      aux5 = fmax(aux4,0.1*aux5);
      if (fabs(aux5) < 1.0e-12) {break;}
    }
    return aux5;
  }
  void cstep(double & stx, double & fx, double & dx, double & sty, double & fy, double & dy, double & stp, double & fp, double & dp, bool & brackt, double & stpmin, double & stpmax) {
    //auxiliary function for More' and Thuente's line search algorithm as extracted from cvsrch in MINPACK project http://www.netlib.org/minpack/
    //translation from Dianne P. O'Leary's matlab code with modifications introduced by James Melville in an R code were also considered in my adaption
    //cstep calculates a safeguarded step for a linesearch and updates interval of uncertainty for function's minimizer
    // -> stx contains the step with the least function value
    // -> stp contains the current step
    //if brackt then a minimizer has been bracketed in an interval of uncertainty with endpoints stx and sty
    // -> stx, fx, and dx are variables which specify the step, the function, and the derivative at the best step obtained so far; the derivative must be negative in the step's direction, i. e., dx and stp-stx must have opposite signs
    // -> sty, fy, and dy are variables which specify the step, the function, and the derivative at the other endpoint of the interval of uncertainty
    // -> stp, fp, and dp are variables which specify the step, the function, and the derivative at the current step; on output stp is set to the new step
    // -> brackt is a logical variable which specifies if a minimizer has been bracketed
    // -> stpmin and stpmax specify lower and upper bounds for the step
    //determine whether the derivatives have opposite sign
    sgnd = dp*dx/fabs(dx);
    if (fp > fx) {                    //case 1: higher function value (bracketed minimum); if cubic step is closer to stx than quadratic step => cubic step is taken, otherwise average
      bound = true;
      stpc = CubicInterpMin(stx,stp,fx,fp,dx,dp,uax);
      stpq = QuadrInterpMin(stx,stp,fx,fp,dx);
      if (fabs(stpc - stx) < fabs(stpq - stx)) {stpf = stpc;}
      else {stpf = stpc + 0.5*(stpq - stpc);}
      brackt = true;
    }
    else if (sgnd < 0.0) {            //case 2: lower function value and derivatives of opposite sign (bracketed minimum); if cubic step closer to stx than quadratic (secant) step => the cubic step is taken, else quadratic step
      bound = false;
      stpc = CubicInterpMin(stx,stp,fx,fp,dx,dp,uax);
      stpq = QuadrInterpGMin(stp,stx,dp,dx);
      if (fabs(stpc - stp) > fabs(stpq - stp)) {stpf = stpc;}
      else {stpf = stpq;}
      brackt = true;
    }
    else if (fabs(dp) < fabs(dx)) {   //case 3: lower function value, derivatives of same sign, and magnitude of derivative decreases; cubic step used if cubic tends to infinity in direction of step or if cubic's minimum beyond stp
                                      //otherwise cubic step is defined to be stpmin or stpmax
                                      //quadratic (secant) step also computed and if minimum is bracketed then step closest to stx is taken, else step farthest away taken
      bound = true;
      stpc = CubicInterpMin(stp,stx,fp,fx,dp,dx,uax);
      if (uax == 0) {
        if (stp > stx) {stpc = stpmax;}
        else {stpc = stpmin;}
      }
      stpq = QuadrInterpGMin(stp,stx,dp,dx);
      if (brackt) {
        if (fabs(stp - stpc) < fabs(stp - stpq)) {stpf = stpc;}
        else {stpf = stpq;}
      }
      else {
        if (fabs(stp - stpc) > fabs(stp - stpq)) {stpf = stpc;}
        else {stpf = stpq;}
      }
    }
    else {                            //case 4: lower function value, derivatives of same sign, and magnitude of derivative does not decrease; if minimum is not bracketed, step is either stpmin or stpmax, else cubic step is taken
      bound = false;
      if (brackt) {
        stpc = CubicInterpMin(stp,sty,fp,fy,dp,dy,uax);
        stpf = stpc;
      }
      else if (stp > stx) {stpf = stpmax;}
      else {stpf = stpmin;}
    }
    //update uncertainty interval; this might be independent on new step or case analysis above
    if (fp > fx) {
      sty = stp;
      fy = fp;
      dy = dp;
    }
    else {
      if (sgnd < 0.0) {
        sty = stx;
        fy = fx;
        dy = dx;
      }
      stx = stp;
      fx = fp;
      dx = dp;
    }
    //compute new step
    stpf = fmin(stpmax,stpf);
    stpf = fmax(stpmin,stpf);
    stp = stpf;
    if ((brackt)&&(bound)) {
      if (sty > stx) {stp = fmin(stx + p66*(sty - stx),stp);}
      else {stp = fmax(stx + p66*(sty - stx),stp);}
    }
  }
  template<class func>
  double MoreThuente(func & ObjFunc, matrixE & step, matrixE & gnew, double fnew, double alpha, bool & skip, bool project, int maxiter = 500) {
    //More' and Thuente line search algorithm as extracted from cvsrch in MINPACK project http://www.netlib.org/minpack/
    //translation from Dianne P. O'Leary's matlab code with modifications introduced by James Melville in an R code are also considered
    //translation of variables in my code and cvsrch is showed below, when assigning values
    //J. J. More, D. J. Thuente, Trans. Math. Soft., 20(3), 286, 1994
    //compute the initial gradient in the search direction and check that s is a descent direction
    skip = true;
    xau = step.rows();
    if (alpha < thresholdzero) {alpha = 1.0;}                   //some sort of restart
    aux1 = 0.0;                                                 // <-> dginit
    for (size_t idx = 0; idx < xau; ++idx) {
      aux1 -= gnew(idx + 1,1)*step(idx + 1,1);
    }
    if (aux1 < thresholdzero) {
      //initialization of variables
      alpha_old = 0.0;
      axu = 0;                              //control variable to check when convergence achieved
      aux2 = alpha_max - alpha_min;         // <-> width
      aux3 = 2.0*aux2;                      // <-> width1
      aux4 = fnew;                          // <-> finit
      aux5 = ftol*aux1;                     // <-> dgtest
      bracket = false;                      // <-> brackt
      stage1 = true;                        // <-> stage1
      ww = 0.0;                             // <-> stx   lowest estimate of step
      zz = 0.0;                             // <-> sty   highest estimate of step
      aux6 = fnew;                          // <-> fx    lowest estimate of f for next step
      aux7 = fnew;                          // <-> fy    highest estimate of f for next step
      aux8 = aux1;                          // <-> dgx   lowest estimate of gradients in next step
      aux9 = aux1;                          // <-> dgy   highest estimate of gradients in next step
      //other variables:     xau = n;    fnew = f;    gnew,Hprm  = g;    step = -s;    alpha = stp;    thresholdzero = xtol;    alpha_min = stpmin;    alpha_max = stpmax; maxiter <=> maxfev
      for (size_t iter = 0; iter < maxiter; ++iter) {
        if (bracket) {
          stmin = fmin(ww,zz);
          stmax = fmax(ww,zz);
        }
        else {
          stmin = ww;
          stmax = alpha + xtrapf*(alpha - ww);
        }
        //ensuring alpha is within bounds
        alpha = fmax(alpha,alpha_min);
        alpha = fmin(alpha,alpha_max);
        if ((bracket)&&((alpha < stmin)||(alpha > stmax)||((stmax - stmin) < thresholdzero*stmax))) {alpha = ww;}
        //evaluate the function and gradient at alpha and compute the directional derivative
        Hhat = step*(alpha - alpha_old);
        ObjFunc.updGeometry(Hhat,project);
        ObjFunc.Calculate(0);
        fnew = ObjFunc.getEnergy(1);
        ObjFunc.gEnergy(Hprm,1,project);
        dg = 0.0;
        for (size_t idx = 0; idx < xau; ++idx) {
          dg -= Hprm(idx + 1,1)*step(idx + 1,1);
        }
        ftest1 = aux4 + alpha*aux5;
        //check for termination
        if ((bracket)&&((alpha < alpha_min)||(alpha > alpha_max))) {
          std::cout << "WARNING: SolverPackage.hpp: solver: MoreThuente(): rounding errors prevent step optimization with sufficient decrease and curvature conditions; maybe revise tolerances" << std::endl;
          axu = -1;
        }
        if ((fabs(alpha - alpha_max) < thresholdzero)&&(fnew < ftest1 + thresholdzero)&&(dg < aux5 + thresholdzero)) {
          std::cout << "WARNING: SolverPackage.hpp: solver: MoreThuente(): step at the upper bound alpha_max" << std::endl;
          axu = -1;
        }
        if ((fabs(alpha - alpha_min) < thresholdzero)&&((fnew > ftest1 - thresholdzero)||(dg > aux5 - thresholdzero))) {
          std::cout << "WARNING: SolverPackage.hpp: solver: MoreThuente(): step at the lower bound alpha_min" << std::endl;
          axu = -1;
        }
        if ((bracket)&&((stmax - stmin) <= thresholdzero*stmax)) {
          std::cout << "WARNING: SolverPackage.hpp: solver: MoreThuente(): relative width of the interval of uncertainty is at most " << thresholdzero << std::endl;
          axu = -1;
        }
        if ((fnew < ftest1 + thresholdzero)&&(fabs(dg) < gtol*fabs(aux1))) {axu = 1;}
        if (axu != 0) {break;}
        alpha_old = alpha;
        if ((stage1)&&(fnew < ftest1 + thresholdzero)&&(dg > fmin(ftol,gtol)*aux1 - thresholdzero)) {stage1 = false;}          //first stage looks for fnew < ftest1 and dg >= 0
        if ((stage1)&&(fnew < aux6 + thresholdzero)&&(fnew > ftest1 - thresholdzero)) {
          //fnew - ftest1 used to predict step
          fm = fnew - alpha*aux5;
          fxm = aux6 - ww*aux5;
          fym = aux7 - zz*aux5;
          dgm = dg - aux5;
          dgxm = aux8 - aux5;
          dgym = aux9 - aux5;
          cstep(ww,fxm,dgxm,zz,fym,dgym,alpha,fm,dgm,bracket,stmin,stmax);
          //reset the function and gradient values for function
          aux6 = fxm + ww*aux5;
          aux7 = fym + zz*aux5;
          aux8 = dgxm + aux5;
          aux9 = dgym + aux5;
        }
        else {cstep(ww,aux6,aux8,zz,aux7,aux9,alpha,fnew,dg,bracket,stmin,stmax);}
        if (bracket) {
          //force sufficient decrease in the size of the interval of uncertainty
          if (abs(zz - ww) > p66*aux3 - thresholdzero) {alpha = ww + p5*(zz - ww);}
          aux3 = aux2;
          aux2 = fabs(zz - ww);
        }
      }
    }
    else {skip = false;}
    aux4 = alpha;
    return aux4;
  }
  template<class func>
  double Fletcher(func & ObjFunc, matrixE & step, matrixE & gnew, double fnew, double alpha, bool project, bool & skip, int maxiter = 500) {
    //Fletcher's line search as described in Seppo Pulkkinen MSc. Thesis "A Review of Methods for Unconstrained Optimization: Theory, Implementation and Testing", 2008
    //{\mu}  -> ftol (0.001)
    //{\eta} -> gtol (0.1)
    //{\tau} -> p5   (0.05)
    //{\chi} -> p66  (9.0)
    skip = false;
    xau = step.rows();
    if (alpha < thresholdzero) {alpha = 1.0;}                   //some sort of restart
    aux1 = 0.0;                                                 //{\phi}^{'}(0) = <{g}_{k}|{d}_{k}> = - <{g}_{k}|{s}_{k}>
    for (size_t idx = 0; idx < xau; ++idx) {
      aux1 -= gnew(idx + 1,1)*step(idx + 1,1);
    }
    aux4 = fnew;                                                //{f}_{initial} = {\phi}(0)
    aux2 = 0.0;                                                 //{alpha}_{l}
    aux3 = alpha_max;                                           //{alpha}_{u}
    aux8 = aux4;                                                //{\phi}_{l}
    aux9 = aux1;                                                //{\phi}^{'}_{l}
    Hhat = ObjFunc.Geometry();
    for (size_t iter_outer = 0; iter_outer < maxiter; ++iter_outer) {
      Hprm = step*alpha;
      ObjFunc.updGeometry(Hprm,project);
      ObjFunc.Calculate(0);
      fnew = ObjFunc.getEnergy(1);                                                     //{\phi}_{t}
      ObjFunc.setGeometry(Hhat);
      axu = 0;                                                                        //counter of inner iterations
      for (size_t iter_inner = 0; iter_inner < maxiter; ++iter_inner, ++axu) {        //interpolation step
        aux7 = QuadrInterpMin(aux2,alpha,aux8,fnew,aux9);                             //get {\alpha}_{new} = {\alpha}_{t}^{+}, which is stored in aux7; alpha contains {\alpha}_{t}
        aux5 = aux3 - aux2;                                                           //{\Delta}{\alpha}
        if (fabs((fnew - aux8)/(aux9*(aux7 - aux2)) - 1.0) < thresholdzero) {aux7 = 0.5*(aux2 + aux7);}        //this is to ensure sufficient decrease in the right direction
        aux7 = fmax(aux2 + p5*aux5,aux7);
        aux7 = fmin(aux3 - p5*aux5,aux7);
        aux3 = alpha;                                                                 //upd {\alpha}_{u}
        alpha = aux7;                                                                 //upd {\alpha} ({\alpha}_{t})
        Hprm = step*aux7;
        ObjFunc.updGeometry(Hprm,project);
        ObjFunc.Calculate(0);
        fnew = ObjFunc.getEnergy(1);
        //check for termination
        if ((fnew < aux4 + alpha*ftol*aux1)) {                                        //sufficient decrease condition
          ObjFunc.gEnergy(Hprm,1,project);
          ObjFunc.setGeometry(Hhat);
          break;
        }
        ObjFunc.setGeometry(Hhat);
        if (fabs(aux7) < 1.0e-14) {break;}
      }
      iter_outer += axu;
      dg = 0.0;
      for (size_t idx = 0; idx < xau; ++idx) {                                    //{\phi}^{'}_{t}
        dg -= Hprm(idx + 1,1)*step(idx + 1,1);
      }
      if (fabs(aux3 - aux2) < thresholdzero) {
        std::cout << "WARNING: SolverPackage.hpp: solver: Fletcher(): limits are too short to continue " << std::endl;
        break;
      }
      if (dg > gtol*aux1) {break;}                                                //curvature condition to leave
      else {                                                                      //extrapolation step
        aux7 = QuadrInterpGMin(aux2,alpha,aux9,dg);                                    //get {\alpha}_{new}
        if (fabs(aux9 - dg) < thresholdzero) {aux7 = 0.5*(aux3 + aux7);}
        else {
          aux6 = aux7 - aux2;                                                          //{\Delta}{\alpha}
          if (aux7 < aux7 + p5*aux6) {aux7 = aux7 + p5*aux6;}
          else if (alpha > aux7 + p66*aux6) {aux7 = aux7 + p66*aux6;}
          if (aux7 > 0.5*(aux3 + aux7)) {aux7 = 0.5*(aux7 + aux3);}
        }
        aux2 = alpha;                                                                  //upd {alpha}_{l}
        aux8 = fnew;
        aux9 = dg;
        alpha = aux7;
      }
    }
    aux4 = aux7;
    return aux4;
  }
  template<class func>
  double DogLeg(func & ObjFunc, matrixE & step, matrixE & hessc, matrixE & gradc, double funcc, double oldrad, int hessupd, bool & TROK, bool & skip) {
    //dogleg trust region algorithm as programmed by C. T. Kelley in matlab for the book (15 Dec. 1997), adapted to my infrastructure and gave robustness
    //calculation of trust region and restricting the step accordingly
    //funcc is the current function value, oldrad gives the entry trust radius, which is to be updated by the function
    //ijob determines what to do now: 1 -> fresh start; 2 -> TR radius reduction in progress; 3 -> attempt TR radius expansion; decided to make this an internal variable
    //axu (idid in original code) takes internally (functionwise) the role of ijob; additionally takes the value of 4 -> expansion step failed, use the last good step, which is not passed to the "outside"
    //pstore takes the role of jdata/sdata of Kelley's algorithm, which is the last successful Hessian stored; this is used to recalculate the Newton step
    //unlike Kelley's routine, only the step is eventually modified according to the function's results
    //bracket (bflag in original code) is internal variable (function level) that determines whether trial point is at the TR boundary and should not be Newton.
    //TROK controls whether the step is taken
    TROK = true;
    skip = true;
    xau = step.rows();
    bracket = 0;
    axu = 1;
    stpc = oldrad;
    if (hessupd == 1) {Hprm = hessc*gradc;}                              //direct methods, using directly the Hessian
    else if (hessupd == -1) {                                            //indirect methods, using the inverse Hessian
      Hprm = gradc;
      Hhat = hessc;
      Solve_Ax_eq_b(Hhat,Hprm);
    }
    if (ijob != 1) {                                                     //use old Hessian matrix if the new step is to discard
      if (hessupd == 1) {
        step = gradc;
        Hhat = pstore;
        Solve_Ax_eq_b(Hhat,step);
      }
      else if (hessupd == -1) {step = pstore*gradc;}
    }
    dgm = 0.0;
    dg = 0.0;
    stpq = 0.0;
    aux8 = 0.0;
    for (size_t idx = 0; idx < xau; ++idx) {
      dgm += gradc(idx + 1,1)*Hprm(idx + 1,1);           //<g|H|g> -> if > 0 then Cauchy point exists
      dg += gradc(idx + 1,1)*gradc(idx + 1,1);           //<g|g>
      stpq += gradc(idx + 1,1)*step(idx + 1,1);          //<g|s>
      aux8 += step(idx + 1,1)*step(idx + 1,1);           //<s|s>
    }
    aux7 = sqrt(dg);
    if (aux7 == 0.0) {aux7 = 1.0;}
    sgnd = 0.0;
    if (dgm > thresholdzero) {                                           //get Cauchy point
      sgnd = dg/dgm;                                                     //sigma in Kelley's code
      if (sgnd*aux7 > stpc) {sgnd = stpc/aux7;}
    }
    else {
      bracket = 1;
      sgnd = 1.0;
    }
    cp = gradc*sgnd;                                                     //steepest descent step
    xn = step;                                                           //Newton step
    if ((dgm > 0.0)&&(stpq > 0.0)&&(aux8 <= stpc*stpc)) {xt = xn;}       //if Cauchy point exists, Newton minimizes and the step is adequate, then Newton
    else if ((dgm > 0.0)&&(stpq <= 0.0)) {xt = cp;}                      //if Cauchy point exists, but Newton is inadequate, then do steepest descent
    else {                                                               //interpolate between Cauchy and Newton
      Hprm = gradc*sgnd;
      Hhat = Hprm - step;
      stmax = 0.0;
      stmin = 0.0;
      stpq = -stpc*stpc;
      for (size_t idx = 0; idx < xau; ++idx) {
        stmax += Hhat(idx + 1,1)*Hhat(idx + 1,1);
        stmin -= 2.0*Hprm(idx + 1,1)*Hhat(idx + 1,1);
        stpq += Hprm(idx + 1,1)*Hprm(idx + 1,1);
      }
      if (fabs(stmax) < thresholdzero) {
        stpq *= stmax;
        stmax = 1.0;
      }
      fym = 0.0;
      if (stmin*stmin - 4.0*stmax*stpq > 0.0) {fym = (sqrt(stmin*stmin - 4.0*stmax*stpq) - stmin)/(2.0*stmax);}
      xt = cp + (xn - cp)*fym;
      bracket = 1;
    }
    step = xt;
    xn = ObjFunc.Geometry();
    ObjFunc.updGeometry(xt,1);
    ObjFunc.Calculate(0);
    fxm = ObjFunc.getEnergy(1);
    ww = fxm - funcc;
    if (hessupd == 1) {Hprm = hessc*step;}            //direct methods, using directly the Hessian
    else if (hessupd == -1) {                         //indirect methods, using the inverse Hessian
      Hprm = step;
      Hhat = hessc;
      Solve_Ax_eq_b(Hhat,Hprm);
    }
    zz = 0.0;
    aux9 = 0.0;
    for (size_t idx = 0; idx < xau; ++idx) {
      zz -= gradc(idx + 1,1)*step(idx + 1,1) + 0.5*step(idx + 1,1)*Hprm(idx + 1,1);
      aux9 += step(idx + 1,1)*step(idx + 1,1);
    }
    ftest1 = ww;
    if (fabs(zz) > 1.0e-14) {ftest1 /= zz;}
    bound = false;
    if ((fabs(ww) < 1.0e-12)&&(fabs(zz) < 1.0e-12)) {bound = true;}
    if (((ftest1 < 0.05)||(ww > 0.0))&&(!bound)) {
      step.zero();
      ObjFunc.setGeometry(xn);
      stpc = 0.5*sqrt(aux9);
      axu = 2;
      if (ijob == 3) {axu = 4;}
      TROK = false;
      skip = false;
    }
    else if ((ftest1 > 0.75)&&(bracket == 1)) {
      stpc = stpc*2;
      axu = 3;
    }
    stpf = stpc;
    if (stpc < thresholdzero) {stpc = 1.0;}
    step *= 1.0/stpc;
    ijob = axu;
    if (ijob == 1) {pstore = hessc;}                                          //store last good Hessian
    else if (ijob == 4) {
      step.zero();
      ObjFunc.setGeometry(xn);
      stpf = oldrad;
      ijob = 1;
      TROK = false;
      skip = false;
    }
    this->setAlphas(stpf + thresholdzero*2.0,stpf + thresholdzero);           //just to avoid trouble or any reduction of step, that is unwanted
    return stpf;
  }
  virtual double Trustregion() {return 0.0;}
};
class SteepDesc: public solver {
  //implementation steepest descent solver
public:
  SteepDesc(int _linesearch = 4, double _alpha_min = 0.00001, double _alpha_max = 1000.0, int _max_bad_alphas = 10000000, double _thresholdzero = 1.0e-6, double _alphathreshzero = 1.0e-3): solver(_linesearch,_alpha_min,_alpha_max,_max_bad_alphas,_thresholdzero,_alphathreshzero) {
    gtol *= 0.1;
  }
  ~SteepDesc() {}
  //getters
  std::string Type() {return "Steepest-Descent";}
  int Degree() {return 1;}                //maximum order of derivatives needed
  //specific methods
  template<class func>
  double Trustregion(func & ObjFunc, matrixE & step, double funcc, double gnorm, double gmax, double fpredicted, double trustradius, double thresh_en, bool & TROK, bool & skip) {
    TROK = true;
    skip = false;
    return 0.0;
  }
  void Step(matrixE & step, matrixE & gx, matrixE & hx, double & fpredicted) {step = gx;}
  void HessianUpdate(matrixE & hessian, matrixE & sk, matrixE & dg, double alpha) {}
};
class ConjGrad: public solver {
  //implementation conjugate gradient solver
  int mthd;        //mthd == 0: Fletcher-Reeves; R. Fletcher, C. M. Reeves, Comput. J., 7, 149, 1964
                   //mthd == 1: Polak-Ribiere
                   //mthd == 2: Polak-Ribiere +, J. C. Gilbert, J. Nocedal, SIAM J. Opt., 2(1), 21, 1992
                   //mthd == 3: PR-FR, J. C. Gilbert, J. Nocedal, SIAM J. Opt., 2(1), 21, 1992
public:
  ConjGrad(int _mthod = 0, int _linesearch = 4, double _alpha_min = 0.00001, double _alpha_max = 1000.0, int _max_bad_alphas = 50, double _thresholdzero = 1.0e-6, double _alphathreshzero = 1.0e-3): solver(_linesearch,_alpha_min,_alpha_max,_max_bad_alphas,_thresholdzero,_alphathreshzero) {
    mthd = _mthod;
    gtol *= 0.5;
  }
  ~ConjGrad() {}
  //getters
  std::string Type() {return "Conjugate-Gradient";}
  int Degree() {return 1;}                //maximum order of derivatives needed
  //specific methods
  template<class func>
  double Trustregion(func & ObjFunc, matrixE & step, double funcc, double gnorm, double gmax, double fpredicted, double trustradius, double thresh_en, bool & TROK, bool & skip) {
    TROK = true;
    skip = false;
    return 0.0;
  }
  void Step(matrixE & step, matrixE & gx, matrixE & hx, double & fpredicted) {
    //for conjugate gradients gx = g(x[i]) and hx = g(x[i-1])
    Hhat = step;
    Hprm = gx;
    if (mthd > 0) {Hprm -= hx;}
    aux = step.rows();
    aux1 = 0.0;
    aux2 = 0.0;
    aux3 = 0.0;
    for (size_t idx = 0; idx < aux; ++idx) {
      aux1 += Hprm(idx + 1,1)*gx(idx + 1,1);
      aux2 += hx(idx + 1,1)*hx(idx + 1,1);
      aux3 += gx(idx + 1,1)*gx(idx + 1,1);
    }
    if (fabs(aux2) < thresholdzero) {         //do steepest descent
      aux1 = 0.0;
      aux3 = 0.0;
      aux2 = 1.0;
    }
    //calculation of beta => aux4
    if (mthd < 2) {aux4 = aux1/aux2;}                       //Fletcher-Reeves, Polak-Ribiere
    else if (mthd == 2) {aux4 = fmax(aux1/aux2,0.0);}       //Polak-Ribiere +
    else if (mthd == 3) {                                   //PR-FR
      aux5 = aux1/aux2;           //beta_PR
      aux6 = aux3/aux2;           //beta_FR
      aux4 = aux5;
      if (aux5 < -aux6) {aux4 = -aux6;}
      else if (aux5 > aux6) {aux4 = aux6;}
    }
    step = gx + Hhat*aux4;
  }
  void HessianUpdate(matrixE & hessian, matrixE & sk, matrixE & dg, double alpha) {}
};
class NewtonRaphson: public solver {
  //implementation Newton-Raphson solver
public:
  NewtonRaphson(double _thresholdzero = 1.0e-6): solver(0,0.1,1.5,100000000,_thresholdzero,1.0e-10) {}
  ~NewtonRaphson() {}
  //getters
  std::string Type() {return "Newton-Raphson";}
  int Degree() {return 2;}                //maximum order of derivatives needed
  //specific methods
  template<class func>
  double Trustregion(func & ObjFunc, matrixE & step, double funcc, double gnorm, double gmax, double fpredicted, double trustradius, double thresh_en, bool & TROK, bool & skip) {
    TROK = true;
    skip = false;
    return 0.0;
  }
  void Step(matrixE & step, matrixE & gx, matrixE & hx, double & fpredicted) {
    step = gx;
    Solve_Ax_eq_b(hx,step);
  }
  void HessianUpdate(matrixE & hessian, matrixE & sk, matrixE & dg, double alpha) {}
};
class BFGSi: public solver {
  //implementation of Broyden-Fletcher-Goldfarb-Shanno family of solvers (missing the names Davidon and Powell in the method's name)
  //with inverted Hessian update
  //mthd defines the Hessian update method:
  int mthd;        //mthd == 0: Murtagh-Sargent; B. A. Murtagh, R. W. H. Sargent, Comp. J., 13(2), 185, 1970
                   //mthd == 1: Fletcher-Powell; R. Fletcher, M. J. D. Powell, Comput. J., 6, 163, 1963
                   //mthd == 2: Shanno's optimal value; D. F. Shanno, Math. Comp., 24(111), 4899, 1970
                   //mthd == 3: Barnes-Rosen; I. G. P. Barnes, Comput J., 8, 66, 1965; E. M. Rosen, Nat. Conf. ACM, 37, 1966
                   //mthd == 4: Broyden-Fletcher-Goldfarb-Shanno (the actual BFGS update)
  double tk;
public:
  BFGSi(int _method = 4, int _linesearch = 4, double _alpha_min = 0.00001, double _alpha_max = 1.5, int _max_bad_alphas = 15, double _thresholdzero = 1.0e-6, double _alphathreshzero = 1.0e-3): solver(_linesearch,_alpha_min,_alpha_max,_max_bad_alphas,_thresholdzero,_alphathreshzero) {
    mthd = _method;
    tk = double(mthd);
    if (mthd == 3) {tk = 0.0;}
  }
  ~BFGSi() {}
  //getters
  std::string Type() {return "BFGSi";}
  int Degree() {return 2;}                //maximum order of derivatives needed
  std::string Method() {
    std::string meth;
    if (mthd == 0) {meth = "Murtagh-Sargent";}
    else if (mthd == 1) {meth = "Fletcher-Powell";}
    else if (mthd == 2) {meth = "Shanno";}
    else if (mthd == 3) {meth = "Barnes-Rosen";}
    else if (mthd == 4) {meth = "Broyden-Fletcher-Goldfarb-Shanno";}
    return meth;
  }
  int method() {return mthd;}
  //setters
  void setMethod(int _mthd) {
    mthd = _mthd;
    tk = double(mthd);
    if (mthd == 3) {tk = 0.0;}
  }
  //specific methods
  template<class func>
  double Trustregion(func & ObjFunc, matrixE & step, double funcc, double gnorm, double gmax, double fpredicted, double trustradius, double thresh_en, bool & TROK, bool & skip) {
    TROK = true;
    skip = false;
    return 0.0;
  }
  void Step(matrixE & step, matrixE & gx, matrixE & hx, double & fpredicted) {step = hx*gx;}
  void HessianUpdate(matrixE & hessian, matrixE & sk, matrixE & dg, double alpha) {
    aux = sk.rows();
    Hprm = hessian*dg;
    if (mthd < 4) {
      if (mthd == 2) {tk = (2.0*alpha - 1.0)/alpha;}                           //Shanno's optimal value
      Hprm += sk*(1.0 - tk)*alpha;
      aux1 = 0.0;
      aux2 = 0.0;
      for (size_t idx = 0; idx < aux; ++idx) {
        aux1 -= sk(idx + 1,1)*dg(idx + 1,1);
        aux2 -= Hprm(idx + 1,1)*dg(idx + 1,1);
      }
      if (aux1 == 0.0) {aux1 = 1.0;}
      if (aux2 == 0.0) {aux2 = 1.0;}
      Hhat = Hprm*Hprm.trans();
      Hhat *= (1.0/aux2);
      hessian += Hhat;
      if (mthd > 0) {
        Hprm = sk*sk.trans();
        Hprm *= (alpha*tk/aux1);
        hessian += Hprm;
      }
    }
    else if (mthd == 4) {
      aux1 = 0.0;
      aux2 = 0.0;
      for (size_t idx = 0; idx < aux; ++idx) {
        aux1 += sk(idx + 1,1)*dg(idx + 1,1);
        aux2 += dg(idx + 1,1)*Hprm(idx + 1,1);
      }
      if (aux1 == 0.0) {aux1 = 1.0;}
      Hhat = sk*sk.trans();
      hessian += Hhat*((aux2/aux1 - alpha)/aux1);
      hessian -= (sk*Hprm.trans() + Hprm*sk.trans())*(1.0/aux1);
    }
  }
};
class BFGSd: public solver {
  //implementation of Broyden-Fletcher-Goldfarb-Shanno family of solvers (missing the names Davidon and Powell in the method's name)
  //with direct Hessian update
  //mthd defines the Hessian update method:
  int mthd;        //mthd == 0: Murtagh-Sargent; B. A. Murtagh, R. W. H. Sargent, Comp. J., 13(2), 185, 1970
                   //mthd == 1: Powell-symmetric-Broyden
                   //mthd == 2: Bofill; J. M. Bofill, J. Comput. Chem., 15(1),1,1994
                   //mthd == 3: Bakken-Helgaker; V. Bakken, T. Helgaker, J. Chem. Phys., 117(20), 9160, 2002
                   //mthd == 4: Broyden-Fletcher-Goldfarb-Shanno
public:
  BFGSd(int _method = 4, int _linesearch = 4, double _alpha_min = 0.00001, double _alpha_max = 1.5, int _max_bad_alphas = 10, double _thresholdzero = 1.0e-6, double _alphathreshzero = 1.0e-3): solver(_linesearch,_alpha_min,_alpha_max,_max_bad_alphas,_thresholdzero,_alphathreshzero) {mthd = _method;}
  ~BFGSd() {}
  //getters
  virtual std::string Type() {return "BFGSd";}
  int Degree() {return 2;}                //maximum order of derivatives needed
  std::string Method() {
    std::string meth;
    if (mthd == 0) {meth = "Murtagh-Sargent";}
    else if (mthd == 1) {meth = "Powell-symmetric-Broyden";}
    else if (mthd == 2) {meth = "Bofill";}                                //this is robust for saddle-points (transition states)
    else if (mthd == 3) {meth = "Bakken-Helgaker";}
    else if (mthd == 4) {meth = "Broyden-Fletcher-Goldfarb-Shanno";}
    return meth;
  }
  int method() {return mthd;}
  //setters
  void setMethod(int _mthd) {mthd = _mthd;}
  //specific methods
  template<class func>
  double Trustregion(func & ObjFunc, matrixE & step, double funcc, double gnorm, double gmax, double fpredicted, double trustradius, double thresh_en, bool & TROK, bool & skip) {
    TROK = true;
    skip = false;
    return 0.0;
  }
  virtual void Step(matrixE & step, matrixE & gx, matrixE & hx, double & fpredicted) {
    step = gx;
    Solve_Ax_eq_b(hx,step);
  }
  void HessianUpdate(matrixE & hessian, matrixE & sk, matrixE & dg, double alpha) {
    aux = sk.rows();
    aux1 = 0.0;                                      //   phi
    aux2 = 0.0;                                      //1 - phi
    if (mthd == 0) {
      aux1 = 0.0;
      aux2 = 1.0;
    }
    else if (mthd == 1) {
      aux1 = 1.0;
      aux2 = 0.0;
    }
    else if ((mthd == 2)||(mthd == 3)) {
      Hhat = hessian*sk*alpha + dg;
      aux3 = 0.0;
      aux4 = 0.0;
      aux5 = 0.0;
      for (size_t idx = 0; idx < aux; ++idx) {
        aux3 += sk(idx + 1,1)*sk(idx + 1,1);
        aux4 += Hhat(idx + 1,1)*Hhat(idx + 1,1);
        aux5 += sk(idx + 1,1)*Hhat(idx + 1,1);
      }
      if (aux3 == 0.0) {aux3 = 1.0;}
      if (aux4 == 0.0) {aux4 = 1.0;}
      aux2 = aux5*aux5/(aux3*aux4);
      aux1 = 1 - aux2;
    }
    else if (mthd == 4) {
      aux1 = 1.0;
      aux2 = 0.0;
    }
    if ((mthd == 0)||(mthd == 2)||(mthd == 3)) {                        //Murtagh-Sargent part of update
      Hhat = dg + hessian*sk*alpha;
      aux3 = 0.0;
      for (size_t idx = 0; idx < aux; ++idx) {
        aux3 += alpha*sk(idx + 1,1)*Hhat(idx + 1,1);
      }
      if (aux3 == 0.0) {aux3 = 1.0;}
      Hprm = Hhat*Hhat.trans()*(1.0/aux3);
      hessian -= Hprm*aux2;
    }
    if ((mthd == 1)||(mthd == 2)) {                                     //Powell-Broyden part of update
      aux3 = 0.0;
      aux4 = 0.0;
      Hhat = hessian*sk*alpha + dg;
      for (size_t idx = 0; idx < aux; ++idx) {
        aux3 += alpha*sk(idx + 1,1)*sk(idx + 1,1);
        aux4 += alpha*sk(idx + 1,1)*Hhat(idx + 1,1);
      }
      if (aux3 == 0.0) {aux3 = 1.0;}
      Hprm = sk*sk.trans()*(aux4/(aux3*aux3));
      hessian += Hprm - (Hhat*sk.trans() + sk*Hhat.trans())*(1.0/aux3);
    }
    if ((mthd == 3)||(mthd == 4)) {                                     //BFGS part of update
      aux3 = 0.0;
      aux4 = 0.0;
      Hhat = hessian*sk;
      for (size_t idx = 0; idx < aux; ++idx) {
        aux3 += alpha*sk(idx + 1,1)*dg(idx + 1,1);
        aux4 += sk(idx + 1,1)*Hhat(idx + 1,1);
      }
      if (aux3 == 0.0) {aux3 = 1.0;}
      if (aux4 == 0.0) {aux4 = 1.0;}
      Hprm = dg*dg.trans()*(1.0/aux3) + Hhat*Hhat.trans()*(1.0/aux4);
      hessian -= Hprm*aux1;
    }
  }
};
class BakerRFO: public BFGSd {
  //implementation of Baker's variant of the RFO minimization algorithm
  //J. Baker, J. Comput. Chem., 7(4), 385, 1986
  //mode determines whether minimum or saddle point is to be optimized: 0 -> minimum; 1 -> TS; this is also number of negative eigenvalues allowed
  //    mthd == 0: Murtagh-Sargent; B. A. Murtagh, R. W. H. Sargent, Comp. J., 13(2), 185, 1970
  //    mthd == 1: Powell-symmetric-Broyden
  //    mthd == 2: Bofill; J. M. Bofill, J. Comput. Chem., 15(1),1,1994
  //    mthd == 3: Bakken-Helgaker; V. Bakken, T. Helgaker, J. Chem. Phys., 117(20), 9160, 2002
  //    mthd == 4: Broyden-Fletcher-Goldfarb-Shanno
  //mfollow is the mode given by user that should be followed by maximization
  //this algorithm does not require any line search, it is embedded in the step
protected:
  int mode;
  int mfollow;                        //vibrational mode to follow
  int firstiter;                      //control variable to determine which mode to follow
  int mode2follow;
  int donewton;                       //variable that controls whether Newton steps are also done
  std::vector<double> mode2vfollow;
  std::vector<double> eigval;
  std::vector<double> auxeigval;
  std::vector<double> lambda;
  double maxratio;
  double minratio;
public:
  int element;
  BakerRFO(int _mthd = 4, int _mode = 0, int _mfollow = 0, double _minratio = 0.0, double _maxratio = 2.0, double _thresholdzero = 1.0e-7): BFGSd(_mthd,7,0.0,100000.0,100000000,_thresholdzero,1.0e-10) {
    //as in MOPAC, minimum step for the Baker optimizer set to 0.005
    mode2follow = 0;
    mfollow = _mfollow;
    if ((_mode == 1)&&(mfollow == 0)) {mfollow = 1;}
    firstiter = 0;
    mode = _mode;
    maxratio = _maxratio;
    minratio = _minratio;
    lambda.resize(2);
    donewton = 0;
    ftol = 0.02;                      //maximum change in function in each step (predicted or actual); only for saddle point optimization
    gtol = 0.001;                     //criterion to determine how close we are to stationary point
    if ((mode == 1)&&(_mthd == 4)) {std::cout << "WARNING: SolverPackage.hpp: BakerRFO: BakerRFO(): search for saddle point with inadequate Hessian update"  << std::endl;}
    stage1 = true;                    //determine whether to dynamically update trust region
  }
  ~BakerRFO() {}
  //getters
  virtual std::string Type() {return "BakerRFO";}
  std::string Mode() {
    std::string modeopt;
    if (mode == 0) {modeopt = "minimization";}
    else if (mode == 1) {modeopt = "saddle-point";}
    return modeopt;
  }
  int Degree() {return 2;}                //maximum order of derivatives needed
  int CalculationMode() {return mode;}
  int ModeFollowed() {return mfollow;}
  int Newton() {return donewton;}
  double MaxRatio() {return maxratio;}
  double MinRatio() {return minratio;}
  double FTol() {return ftol;}
  double GTol() {return gtol;}
  bool DynamicAdjustTR() {return stage1;}
  //setters
  void setModeFollowed(int _mfollow) {mfollow = _mfollow;}
  void resetFirstIter() {firstiter = 0;}
  void setMode(int _mode) {mode = _mode;}
  void setMaxRatio(double _maxratio) {maxratio = _maxratio;}
  void setMinRatio(double _minratio) {minratio = _minratio;}
  void setDynamicAdjustTR(bool newvar) {stage1 = newvar;}
  void setNewton(int _donewton) {donewton = _donewton;}
  void setFTol(double value) {ftol = value;}
  void setGTol(double value) {gtol = value;}
  //specific methods
  virtual void Step(matrixE & step, matrixE & gx, matrixE & hx, double & fpredicted) {
    //fpredicted comes in as trust region, goes out as predicted decrease in f
    aux = hx.rows();
    //diagonalize Hessian
    Hhat = hx;
    eigval = MatDiag(Hhat);
    //transform gradient
    Hprm = Hhat.trans()*gx;                    //Hprm = Fbar from paper
    //count negative eigenvalues
    uax = 0;
    for (size_t idx = 0; idx < eigval.size(); ++idx) {
      if (fabs(eigval[idx]) < thresholdzero) {              //remove numerical noise
        eigval[idx] = 0.0;
        Hprm(idx + 1,1) = 0.0;
      }
      if (eigval[idx] < 0.0) {
        ++uax;                                              //number of negative eigenvalues
        axu = idx + 1;                                      //the largest eigenmode with negative eigenvalue
      }
    }
    xau = 1;
    if (mfollow != 0) {                                     //if we are interested in following a given mode, then this is where to look into
      if (firstiter == 0) {                                     //save mode to vector, if in the first iteration
        mode2vfollow.resize(aux);
        for (size_t idx = 0; idx < aux; ++idx) {                //copy eigenvector into mode2vfollow
          mode2vfollow[idx] = Hhat(idx + 1,mfollow);
        }
        firstiter = 1;
      }
      else {                                                    //search for closest mode
        xau = 0;
        aux1 = 0.0;
        for (size_t icol = 0; icol < aux; ++icol) {
          aux2 = 0.0;
          for (size_t irow = 0; irow < aux; ++irow) {
            aux2 += mode2vfollow[irow]*Hhat(irow + 1,icol + 1);
          }
          if (fabs(aux2) > aux1) {
            aux1 = fabs(aux2);
            xau = icol + 1;
          }
        }
        if (xau == 0) {std::cout << "WARNING: SolverPackage.hpp: BakerRFO: BakerRFO(): lost mode to be followed" << std::endl;}
        mfollow = xau;
      }
      mode2follow = mfollow;                                              //mode2follow is actually the mode to be followed
    }
    else {mode2follow = mfollow;}                                         //if nothing is predefined, just follow the first mode
    if (xau != 0) {                                                       //xau = 0 means that mode to be followed was lost in saddle point optimization; for minimizations it makes no sense
      if ((uax != mode)||(!donewton)) {
P_RFO:
        if ((mode == 1)||(uax != mode)) {                                 //do PRFO
          xt = Identity(aux);                                             //this is the augmented Hessian
          size_t cnt = 1;
          for (size_t idx = 0; idx < aux; ++idx) {
            if (idx + 1 == mode2follow) {continue;}
            xt(cnt,cnt) = eigval[idx];
            xt(aux,cnt) = Hprm(idx + 1,1);
            xt(cnt,aux) = Hprm(idx + 1,1);
            ++cnt;
          }
          xt(aux,aux) = 0.0;
          auxeigval = MatDiag(xt);
          if (firstiter == 1) {
            //initialize lambdas
            lambda[0] = auxeigval[aux - 1];
            lambda[1] = auxeigval[0];
            firstiter = 2;
          }
          element = mode2follow;
          if (mode2follow == 0) {element = 1;}
          aux5 = 0.5*(eigval[element - 1] + sqrt(eigval[element - 1]*eigval[element - 1] + 4.0*Hprm(element,1)*Hprm(element,1)));
          if (aux5 < lambda[0]) {lambda[0] = aux5;}
          aux5 = lambda[0];
          for (size_t idx = 0; idx < aux; ++idx) {
            if (auxeigval[idx] < aux5) {aux5 = auxeigval[idx];}
          }
          if (aux5 > lambda[1]) {lambda[1] = aux5;}
          //corrections... should never take place...
          xau = 1;
          if (mode2follow != 1) {xau = 0;}
          if (lambda[1] > eigval[xau]) {lambda[1] = eigval[xau] - 1.0;}
        }
        else {                                               //do RFO
RFO:
          xt = Identity(aux + 1);
          for (size_t idx = 0; idx < aux; ++idx) {
            xt(idx + 1,idx + 1) = eigval[idx];
            xt(aux + 1,idx + 1) = Hprm(idx + 1,1);
            xt(idx + 1,aux + 1) = Hprm(idx + 1,1);
          }
          xt(aux + 1,aux + 1) = 0.0;
          auxeigval = MatDiag(xt);
          if (firstiter == 0) {                                    //initialize lambda
            lambda[0] = auxeigval[0];
            firstiter = 1;
          }
          aux5 = eigval[0];
          for (size_t idx = 0; idx < aux + 1; ++idx) {
            if (auxeigval[idx] < aux5) {aux5 = auxeigval[idx];}
          }
          if (aux5 > lambda[0]) {lambda[0] = aux5;}
          lambda[1] = lambda[0];
        }
      }
      //modify Fbar
      aux4 = double(uax != mode) + double(!donewton);
      if (int(aux4) > 1) {aux4 = 1.0;}
      for (size_t idx = 0; idx < aux; ++idx) {
        aux5 = lambda[1];        //lambda_n
        if (mode2follow - 1 == idx) {aux5 = lambda[0];}
        if (fabs(eigval[idx] - aux4*aux5) < thresholdzero) {Hprm(idx + 1,1) = 0.0;}       //remove numerical noise
        else {Hprm(idx + 1,1) /= (eigval[idx] - aux4*aux5);}
      }
      //get step
      step = Hhat*Hprm;
      aux5 = 0.0;
      for (size_t idx = 0; idx < aux; ++idx) {
        aux5 += step(idx + 1,1)*step(idx + 1,1);
      }
      aux3 = sqrt(3.0*aux5/double(aux));
      if (aux3 > fpredicted) {
        step *= fpredicted/aux3;
        aux3 = fpredicted;
      }
      //calculate the predicted decrease in energy
      aux6 = 0.0;
      for (size_t idx = 0; idx < aux; ++idx) {
        aux7 = lambda[1];        //lambda_n
        if (mode2follow - 1 == idx) {aux7 = lambda[0];}
        aux6 += Hprm(idx + 1,1)*Hprm(idx + 1,1)*(aux7 - 0.5*eigval[idx]);
      }
      if ((aux6 > 0.0)&&(mode == 0)) {
        lambda[0] *= 2.0;
        lambda[1] *= 2.0;
        goto P_RFO;
      }
      fpredicted = aux6;
    }
    else {step.zero();}
  }
  template<class func>
  double Trustregion(func & ObjFunc, matrixE & step, double funcc, double gnorm, double gmax, double fpredicted, double trustradius, double thresh_en, bool & TROK, bool & skip) {
    //check how predicted and actual function reduction are going; only makes sense for minimization
    aux7 = trustradius;
    Hhat = step;                       //test taking of step
    cp = ObjFunc.Geometry();
    aux = cp.rows();
    ObjFunc.updGeometry(Hhat,1);
    ObjFunc.Calculate(0);
    aux1 = ObjFunc.getEnergy(1);
    aux1 -= funcc;                     //actual decrease in function
    aux2 = aux1/fpredicted;            //ratio
    ObjFunc.gEnergy(Hprm,1,1,thresholdzero);
    aux3 = 0.0;
    aux5 = 0.0;
    stpf = 0.0;
    axu = Hprm.rows();
    for (size_t idx = 0; idx < axu; ++idx) {
      if (fabs(Hprm(idx + 1,1)) > stpf) {stpf = fabs(Hprm(idx + 1,1));}
      aux3 += Hprm(idx + 1,1)*Hprm(idx + 1,1);
      aux5 += step(idx + 1,1)*step(idx + 1,1);
    }
    aux4 = sqrt(aux3/double(aux));     //gradient norm
    aux6 = sqrt(aux5/double(aux));     //step norm
    //determine whether to reject the step and reduce trust region
    bracket = false;
    if ((mode == 0)&&(aux1 > 10.0*thresh_en)) {bracket = true;}                                                                                        //if there is no decrease in function during minimization, reject the step
    else if ((mode == 0)&&(aux4 > gnorm)&&(stpf > gmax)) {bracket = true;}                                                                             //if the gradient increases, reject the step
    else if ((mode == 1)&&(aux2 > this->MaxRatio())&&((fabs(fpredicted) > ftol)||(fabs(aux1) > ftol))) {bracket = true;}                               //if the step causes too many changes, reject the step
    std::cout << bracket << " dEactual = " << aux1 << " dEpred = " << fpredicted << " |g|new = " << aux4 << " |g|old = " << gnorm << std::endl;
    if (bracket) {                            //reject step, reduce trust region
      step.zero();
      TROK = false;
      aux7 = fmin(trustradius,aux6)/2.0;
    }
    else if (stage1) {                        //dynamically adjust trust region as recommended by Fletcher
      skip = true;
      if ((aux4 > gtol)||(fabs(fpredicted) > ftol)||(fabs(aux1) > ftol)) {
        //change in trust region as in MOPAC because we focus on updating the Hessian and not on recalculating it every iteration
        if (mode == 1) {
          if ((aux2 < 0.1)||(aux2 > 3.0)) {aux7 = fmin(trustradius,aux6)/2.0;}                                             //this is an extra decrease for saddle point optimization
          else if ((aux2 > 0.75)&&(aux2 < 4.0/3.0)&&(aux6 > trustradius - thresholdzero)) {aux7 *= sqrt(2.0);}             //increase trust region
        }
        else if (mode == 0) {
          if ((aux2 > 0.5)&&(aux6 < trustradius - thresholdzero)) {
            aux7 *= sqrt(2.0);}                                                                                            //increase trust region
          if (fabs(aux2 - 1.0) < 0.1) {aux7 *= sqrt(2.0);}                                                                 //be very audacious
        }
      }
      else if ((aux4 < gtol)&&(fabs(fpredicted) < ftol)&&(fabs(aux1) < ftol)) {aux7 = fmax(trustradius,0.1);}
    }
    else {skip = true;}                       //don't adjust trust region, but accept step
    if (!skip) {ObjFunc.setGeometry(cp);}     //resetting geometry, since we reject the step
    aux6 = 1.0/aux7;
    step *= aux6;
    return aux7;
  }
};
class Wales: public BakerRFO {
  //implementation of David Wales' variant of Newton-Raphson
  //D. J. Wales, Mol. Phys., 74 (1), 1, 1991
  //D. J. Wales, J. Chem. Phys., 101 (5), 3750, 1994 (this is the main paper followed)
  bool safemode;
public:
  Wales(int _mthd = 4, int _mode = 0, int _mfollow = 0, double _minratio = 0.0, double _maxratio = 0.0, bool _safemode = 1, double _thresholdzero = 1.0e-7): BakerRFO(_mthd,_mode,_mfollow,_minratio,_maxratio,_thresholdzero) {
    safemode = _safemode;
    if ((mode == 1)&&(_mthd == 4)) {std::cout << "WARNING: SolverPackage.hpp: Wales: Wales(): search for saddle point with inadequate Hessian update"  << std::endl;}
  }
  ~Wales() {}
  //getters
  std::string Type() {return "Wales-NR";}
  bool SafeMode() {return safemode;}
  //setters
  void setSafeMode(bool _safemode) {safemode = _safemode;}
  //specific methods
  void Step(matrixE & step, matrixE & gx, matrixE & hx, double & fpredicted) {
    //fpredicted comes in as trust region, goes out as predicted decrease in f
    aux = hx.rows();
    //diagonalize Hessian
    Hhat = hx;
    eigval = MatDiag(Hhat);
    //what type of calculation?
    xau = 1;
    if (mfollow != 0) {                                     //if we are interested in following a given mode, then this is where to look into
      if (firstiter == 0) {                                     //save mode to vector, if in the first iteration
        mode2vfollow.resize(aux);
        for (size_t idx = 0; idx < aux; ++idx) {                //copy eigenvector into mode2vfollow
          mode2vfollow[idx] = Hhat(idx + 1,mfollow);
        }
        firstiter = 1;
      }
      else {                                                    //search for closest mode
        xau = 0;
        aux1 = 0.0;
        for (size_t icol = 0; icol < aux; ++icol) {
          aux2 = 0.0;
          for (size_t irow = 0; irow < aux; ++irow) {
            aux2 += mode2vfollow[irow]*Hhat(irow + 1,icol + 1);
          }
          if (fabs(aux2) > aux1) {
            aux1 = fabs(aux2);
            xau = icol + 1;
          }
        }
        if (xau == 0) {std::cout << "WARNING: SolverPackage.hpp: Wales: Wales(): lost mode to be followed" << std::endl;}
        mfollow = xau;
      }
      mode2follow = mfollow;                                              //mode2follow is actually the mode to be followed
      //Hprm(mode2follow,1) *= -1.0;                                //put it for maximization
    }
    else {mode2follow = mfollow;}                                         //if nothing is predefined, just follow mfollow
    if (xau != 0) {                                          //xau = 0 means that mode to be followed was lost in saddle point optimization; for minimizations it makes no sense to ask that question
      //transform gradient
      Hprm = Hhat.trans()*gx;                                               //Hprm = Fbar from paper
      //get some auxiliary quantities, like the {\beta}_{i} (stored in Hprm) and {\phi}
      aux4 = 0.0;
      aux8 = 0.0;
      stmin = 0.0;
      stmax = 0.0;
      for (size_t idx = 0; idx < aux; ++idx) {                                  //calculate 2Fi/|bi| = {\beta}_{i}
        if (fabs(eigval[idx]) < thresholdzero) {
          eigval[idx] = 0.0;
          Hprm(idx + 1,1) = 0.0;
        }
        else {
          fxm = eigval[idx];
          fm = Hprm(idx + 1,1);                                             //this is {\bar{F}}_{i}
          aux4 += fxm;                                                              //for the average eigenvalue
          aux6 = 2.0*fm/fabs(fxm);                                                  //{\beta}_{i}
          aux9 = -pow(-1.0,(mfollow != 0)*(mfollow == idx));                        //-1 -> minimization; +1 -> maximization
          aux7 = fxm + aux9*0.5*fabs(fxm)*(1.0 + sqrt(1.0 + aux6*aux6));            //{\lambda}_{i}
          aux8 += (aux7 - 0.5*fxm)*fm*fm/((aux7 - fxm)*(aux7 - fxm));               //predicted decrease in energy
          Hprm(idx + 1,1) = -aux9*aux6/(1.0 + sqrt(1.0 + aux6*aux6));       //-{h}_{i}
        }
        stmin += gx(idx + 1,1)*gx(idx + 1,1);                       //norm of gradient
        stmax += Hprm(idx + 1,1)*Hprm(idx + 1,1);                   //norm of step
      }
      aux6 = 2.0*double(safemode)*aux4*stmax/(double(aux)*sqrt(stmin));             //{\phi}{h}^{2}
      //calculate the shifted step
      if ((fabs(aux6) > thresholdzero)&&(1.0 + aux6 >= 0.0)) {
        aux7 = 2.0/(1.0 + sqrt(1.0 + aux6));
        for (size_t idx = 0; idx < aux; ++idx) {
          Hprm(idx + 1,1) *= aux7;                                          //{h}_{i}^{'}
        }
      }
      //get step 
      step = Hhat*Hprm;
      aux5 = 0.0;
      for (size_t idx = 0; idx < aux; ++idx) {
        aux5 += step(idx + 1,1)*step(idx + 1,1);
      }
      aux3 = sqrt(3.0*aux5/double(aux));
      if (aux3 > fpredicted) {
        step *= fpredicted/aux3;
        aux3 = fpredicted;
      }
      fpredicted = aux8;
    }
    else {step.zero();}
  }
};
class CGBFGS: public solver {
  //implementation of a mixed conjugate gradient-BFGS solver
  //this solver uses conjugate gradient until a certain gradient norm and then switches to a quasi-Newton method (with inverse Hessian update) to finalize the geometry optimization
  //for consistency reasons only More'-Thuente line search used, because this is the only line-search consistent for both solvers
  int actualsolver;  //1 -> conjugate gradient
                     //2 -> moved to quasi-Newton
  //mthd1 determines the conjugate gradient solver
  int mthd1;          //mthd1 == 0: Fletcher-Reeves; R. Fletcher, C. M. Reeves, Comput. J., 7, 149, 1964
                      //mthd1 == 1: Polak-Ribiere
                      //mthd1 == 2: Polak-Ribiere +, J. C. Gilbert, J. Nocedal, SIAM J. Opt., 2(1), 21, 1992
                      //mthd1 == 3: PR-FR, J. C. Gilbert, J. Nocedal, SIAM J. Opt., 2(1), 21, 1992
  //mthd2 determines the quasi-Newton solver
  int mthd2;          //mthd2 == 0: Murtagh-Sargent; B. A. Murtagh, R. W. H. Sargent, Comp. J., 13(2), 185, 1970
                      //mthd2 == 1: Fletcher-Powell; R. Fletcher, M. J. D. Powell, Comput. J., 6, 163, 1963
                      //mthd2 == 2: Shanno's optimal value; D. F. Shanno, Math. Comp., 24(111), 4899, 1970
                      //mthd2 == 3: Barnes-Rosen; I. G. P. Barnes, Comput J., 8, 66, 1965; E. M. Rosen, Nat. Conf. ACM, 37, 1966
                      //mthd2 == 4: Broyden-Fletcher-Goldfarb-Shanno (the actual BFGS update)
  int firstiter;      //control variable to determine whether in first iteration
  double tk;
  double gthreshold;
public:
  CGBFGS(int _mthod1 = 1, int _mthod2 = 4, double _alpha_min = 0.00001, double _alpha_max = 1000.0, int _max_bad_alphas = 50, double _thresholdzero = 1.0e-6, double _alphathreshzero = 1.0e-3): solver(4,_alpha_min,_alpha_max,_max_bad_alphas,_thresholdzero,_alphathreshzero) {
    actualsolver = 1;
    mthd1 = _mthod1;
    mthd2 = _mthod2;
    tk = double(mthd2);
    if (mthd2 == 3) {tk = 0.0;}
    gtol *= 0.5;
    gthreshold = 0.01;
    firstiter = 1;
  }
  ~CGBFGS() {}
  //getters
  std::string Type() {return "CGBFGS";}
  int Degree() {return 2;}                //maximum order of derivatives needed
  double SwitchThreshold() {return gthreshold;}
  int CurrentSolver() {return actualsolver;}
  std::string Method() {
    std::string meth;
    if (mthd2 == 0) {meth = "Murtagh-Sargent";}
    else if (mthd2 == 1) {meth = "Fletcher-Powell";}
    else if (mthd2 == 2) {meth = "Shanno";}
    else if (mthd2 == 3) {meth = "Barnes-Rosen";}
    else if (mthd2 == 4) {meth = "Broyden-Fletcher-Goldfarb-Shanno";}
    return meth;
  }
  int method() {return mthd2;}
  //setters
  void setSwitchThreshold(double _gthreshold) {gthreshold = _gthreshold;}
  void setMethod(int _mthd) {
    mthd2 = _mthd;
    tk = double(mthd2);
    if (mthd2 == 3) {tk = 0.0;}
  }
  //specific methods
  template<class func>
  double Trustregion(func & ObjFunc, matrixE & step, double funcc, double gnorm, double gmax, double fpredicted, double trustradius, double thresh_en, bool & TROK, bool & skip) {
    TROK = true;
    skip = false;
    return 0.0;
  }
  void Step(matrixE & step, matrixE & gx, matrixE & hx, double & fpredicted) {
    aux = step.rows();
    if (actualsolver == 1) {
      if (firstiter) {
        pstore.resize(aux,1);
        firstiter = 0;
      }
      aux1 = 0.0;
      for (size_t idx = 0; idx < aux; ++idx) {
        aux1 += gx(idx + 1,1)*gx(idx + 1,1);
      }
      aux2 = sqrt(3.0*aux1/double(aux));
      if (aux2 < gthreshold) {
        fpredicted = -1.0;                                //in this case use fpredicted to pass to the outside the message that solver is changing
        actualsolver = 2;
      }
    }
    if (actualsolver == 1) {
      //for conjugate gradients gx = g(x[i]) and hx = g(x[i-1])
      Hhat = step;
      Hprm = gx;
      if (mthd1 > 0) {Hprm -= pstore;}
      aux1 = 0.0;
      aux2 = 0.0;
      aux3 = 0.0;
      for (size_t idx = 0; idx < aux; ++idx) {
        aux1 += Hprm(idx + 1,1)*gx(idx + 1,1);
        aux2 += pstore(idx + 1,1)*pstore(idx + 1,1);
        aux3 += gx(idx + 1,1)*gx(idx + 1,1);
      }
      if (fabs(aux2) < thresholdzero) {         //do steepest descent
        aux1 = 0.0;
        aux3 = 0.0;
        aux2 = 1.0;
      }
      //calculation of beta => aux4
      if (mthd1 < 2) {aux4 = aux1/aux2;}                       //Fletcher-Reeves, Polak-Ribiere
      else if (mthd1 == 2) {aux4 = fmax(aux1/aux2,0.0);}       //Polak-Ribiere +
      else if (mthd1 == 3) {                                   //PR-FR
        aux5 = aux1/aux2;           //beta_PR
        aux6 = aux3/aux2;           //beta_FR
        aux4 = aux5;
        if (aux5 < -aux6) {aux4 = -aux6;}
        else if (aux5 > aux6) {aux4 = aux6;}
      }
      step = gx + Hhat*aux4;
      pstore = gx;
    }
    else {step = hx*gx;}
  }
  void HessianUpdate(matrixE & hessian, matrixE & sk, matrixE & dg, double alpha) {
    //even though there might be no quasi-Newton updates, Hessian is still accumulated
    aux = sk.rows();
    Hprm = hessian*dg;
    if (mthd2 < 4) {
      if (mthd2 == 2) {tk = (2.0*alpha - 1.0)/alpha;}                           //Shanno's optimal value
      Hprm += sk*(1.0 - tk)*alpha;
      aux1 = 0.0;
      aux2 = 0.0;
      for (size_t idx = 0; idx < aux; ++idx) {
        aux1 -= sk(idx + 1,1)*dg(idx + 1,1);
        aux2 -= Hprm(idx + 1,1)*dg(idx + 1,1);
      }
      if (aux1 == 0.0) {aux1 = 1.0;}
      if (aux2 == 0.0) {aux2 = 1.0;}
      Hhat = Hprm*Hprm.trans();
      Hhat *= (1.0/aux2);
      hessian += Hhat;
      if (mthd2 > 0) {
        Hprm = sk*sk.trans();
        Hprm *= (alpha*tk/aux1);
        hessian += Hprm;
      }
    }
    else if (mthd2 == 4) {
      aux1 = 0.0;
      aux2 = 0.0;
      for (size_t idx = 0; idx < aux; ++idx) {
        aux1 += sk(idx + 1,1)*dg(idx + 1,1);
        aux2 += dg(idx + 1,1)*Hprm(idx + 1,1);
      }
      if (aux1 == 0.0) {aux1 = 1.0;}
      Hhat = sk*sk.trans();
      hessian += Hhat*((aux2/aux1 - alpha)/aux1);
      hessian -= (sk*Hprm.trans() + Hprm*sk.trans())*(1.0/aux1);
    }
  }
};
template <class func, class slvr> 
bool SolverOpt(func & ObjFunc, slvr & Solver, int HessianType = 0, int convergencemode = 0, double thresh_energy = 1.0e-7, double thresh_grad = 1.0e-7, double thresh_step = 1.0e-7, int maxiter = 500) {
  //geometry optimization algorithm based on the main solvers available
  //convergencemode:
  //     0 -> energy + gradient
  //     1 -> gradient + displacement
  //     2 -> energy + displacement
  //     3 -> energy + gradient + displacement
  std::string type = Solver.Type();
  int degree = Solver.Degree();
  int linesearch = Solver.LineSearch();
  int doNewton = Solver.Newton();
  int szgeom = ObjFunc.Geometry().rows();
  int szgrad = 1;
  int spdim = ObjFunc.Geometry().cols();
  int failedsteps = 0;
  int totalnrsteps = 0;
  int hessupd = 1;                                                                 //algorithm with normal (1) or inverted (-1) Hessian
  matrixE sk(szgrad,1);                                                            //step; in vector form
  matrixE grad_new(1,1);                                                           //new gradient; in vector form
  matrixE grad_old(1,1);                                                           //old gradient; in vector form
  matrixE delta_g(1,1);                                                            //new gradient; in vector form
  matrixE Hessian(1,1);                                                            //hessian matrix; in vector form
  double normA;
  double normB;
  double normG;                                                                    //control variable for checking whether gradient norm is below threshold
  double maxG;                                                                     //control variable for maximum absolute value of gradient
  double normH;                                                                    //control variable for Herbert's convergence test: |x(k) - x(k-1)| < threshold
  double maxH;                                                                     //control variable for maximum absolute value of |x(k) - x(k-1)|
  double normE;                                                                    //control variable for change in energy
  double fold = 0.0;
  double fpredicted;
  double dampfactor = 0.01;
  if ((type == "BakerRFO")&&(Solver.Mode() == "saddle-point")) {dampfactor = -10000.0;}
  ObjFunc.Calculate(0);
  double fnew = ObjFunc.getEnergy(1);
  double alpha = 1.0;
  //control variables to determine when convergence obtained
  bool convergedG = false;                                      //convergence of gradient
  bool convergedS = false;                                      //convergence of step
  bool convergedE = false;                                      //convergence of energy
  bool converged = false;
  bool TROK;
  bool skip;
  bool recalculateHess;
  bool projecthessian = false;
  if ((degree == 2)&&(ObjFunc.NAtoms() > 2)) {projecthessian = true;}
  //this requires a bit more detail: matrix form means that tensor is a Nx3 matrix containing a geometry
  //in vector form means that the matrix is vectorized, as is done in gradient and hessian calculation
  if (degree == 2) {                           //get Hessian or the approximation to it
    ObjFunc.hEnergy(delta_g,HessianType,projecthessian,1,dampfactor,thresh_grad,maxiter);
    if ((type == "BFGSi")||(type == "CGBFGS")) {                                                         //this is not really the hessian but its inverse
      Hessian = delta_g.inverse();
      hessupd = -1;
    }
    else if ((type == "Newton-Raphson")||(type == "BFGSd")||(type == "BakerRFO")||(type == "Wales-NR")) {Hessian = delta_g;}
  }
  ObjFunc.gEnergy(grad_new,1,projecthessian,thresh_grad);
  szgrad = grad_new.rows();                                                        //dimensions to optimize
  grad_old.resize(szgrad,1);
  sk.resize(szgrad,1);
  if (linesearch == 4) {
    alpha = 0.1;
    if (type == "Conjugate-Gradient") {alpha *= 3.0;}
  }
  else if ((linesearch == 6)||(linesearch == 7)) {
    Solver.pstore = Hessian;
    normG = 0.0;
    for (size_t idcoord = 0; idcoord < szgrad; ++idcoord) {
      normA = grad_new(idcoord + 1,1);
      normG += normA*normA;
    }
    alpha = fmin(sqrt(normG/double(szgeom)),10.0);
  }
  else if ((type == "BakerRFO")&&(Solver.Mode() == "saddle-point")) {alpha = 0.2;}
  //begin the iteration
  for (int iter = 0; iter < maxiter; ++iter, ++totalnrsteps) {
    TROK = true;
    skip = false;
    recalculateHess = false;
    std::cout << "iter = " << iter + 1 << std::endl;
    //get search direction
    fpredicted = alpha;
    if (degree == 1) {Solver.Step(sk,grad_new,grad_old,fpredicted);}
    else if (degree == 2) {Solver.Step(sk,grad_new,Hessian,fpredicted);}           //fpredicted contains the predicted change in energy for trust region methods implemented as such
    normG = 0.0;
    maxG = 0.0;
    normB = 0.0;                 //this is to check whether there is trouble with the calculated step
    for (size_t idcoord = 0; idcoord < szgrad; ++idcoord) {
      normA = fabs(grad_new(idcoord + 1,1));
      if (maxG < normA) {maxG = normA;}
      normG += normA*normA;
      normA = fabs(sk(idcoord + 1,1));
      normB += normA*normA;
    }
    normH = sqrt(normB/double(szgeom));
    if (normH > 50.0) {                      //if the step is absurd, do something else; this something else should avoid Hessian, since we are here because it is ill-conditioned
      if ((normH < 1.0e-6)&&(Solver.Mode() != "saddle-point")) {
        recalculateHess = true;
        sk = grad_new;
      }
    }
    normA = sqrt(normG/double(szgeom));
    normG = normA;
    if ((normG/thresh_grad < 1.0 + thresh_grad)&&(0.5*maxG/thresh_grad < 1.0 + thresh_grad)) {convergedG = true;}                    //converged on gradients?
    else {convergedG = false;}
    std::cout << "        |g(k)| = " << normG << "         max|g(k)| = " << maxG << "        congG = " << convergedG << std::endl;
    if (linesearch == 1) {alpha = Solver.Davidon(ObjFunc,sk,grad_new,Hessian,fnew,projecthessian,skip);}                                //Davidon line search
    else if (linesearch == 2) {alpha = Solver.BarzilaiBorwein(sk,grad_new,grad_old,iter%2,alpha);}                                   //Barzilai-Borwein
    else if (linesearch == 3) {alpha = Solver.Armijo(ObjFunc,sk,grad_new,fold,projecthessian,skip,0.5,0.75,maxiter);}                   //Armijo rule
    else if (linesearch == 4) {alpha = Solver.MoreThuente(ObjFunc,sk,grad_new,fnew,alpha,skip,projecthessian,maxiter);}                 //More'-Thuente
    else if (linesearch == 5) {alpha = Solver.Fletcher(ObjFunc,sk,grad_new,fnew,alpha,projecthessian,skip,maxiter);}                    //Fletcher
    else if (linesearch == 6) {alpha = Solver.DogLeg(ObjFunc,sk,Hessian,grad_new,fnew,alpha,hessupd,TROK,skip);}                 //Dogleg
    else if (linesearch == 7) {alpha = Solver.Trustregion(ObjFunc,sk,fnew,normG,maxG,fpredicted,alpha,thresh_energy,TROK,skip);}     //Trust-region things for such methods
    //corrective action needed?
    Solver.UpdateAlphaControl(alpha,iter);
    if (Solver.SteepestDescentStep()) {                                                                                              //if there are too many bad steps, then do steepest descent
      std::cout << "WARNING: SolverPackage.hpp: SolverOpt(): doing steepest descent step" << std::endl;
      Solver.ResetAlphaControl();                                                                                                    //reset counter
      sk = grad_new;
      alpha = 1.0;
      skip = false;
    }
    else if (!TROK) {
      std::cout << "     -> iteration rejected" << std::endl;
      normH *= 0.5;
      maxH *= 0.5;
      normE *= 0.5;
      ++failedsteps;
      if (failedsteps > 10) {                                              //do steepest descent
        sk = grad_new;
        if (linesearch != 6) {alpha = 1.0;}
        skip = false;
        recalculateHess = true;
        failedsteps = 0;
        TROK = true;
      }
      else {--iter;}
    }
    else {                                                                                                                  //restrict alpha?
      if (alpha < Solver.AlphaMin()) {alpha = Solver.AlphaMin();}
      if (alpha > Solver.AlphaMax()) {alpha = Solver.AlphaMax();}
    }
    //update x and get norm of step
    if (TROK) {
      failedsteps = 0;
      normB = 0.0;
      maxH = 0.0;
      for (size_t idcoord = 0; idcoord < szgrad; ++idcoord) {
        normA = fabs(sk(idcoord + 1,1));
        if (maxH < normA) {maxH = normA;}
        normB += normA*normA;
      }
      normH = fabs(alpha)*sqrt(normB/double(szgeom));
      maxH *= fabs(alpha);
      //update function, gradient
      fold = fnew;
      grad_old = grad_new;
      if (!skip) {
        delta_g = sk*alpha;
        ObjFunc.updGeometry(delta_g,projecthessian);
        ObjFunc.Calculate(0);
      }
      fnew = ObjFunc.getEnergy(1);
      normE = fabs(fnew - fold);
      std::cout << "        |Delta E| = " << normE << std::endl;
      ObjFunc.gEnergy(grad_new,1,projecthessian,thresh_grad);
      if (degree == 2) {
        if (!recalculateHess) {
          recalculateHess = (int(fpredicted) == -1) + (type == "Newton-Raphson");
          if (((type == "BakerRFO")||(type == "Wales-NR"))&&(((iter + 1)%5 == 0)&&(doNewton))) {
            recalculateHess = true;
            alpha = sqrt(normG/double(szgeom));                                            //reset trust region
          }
        }
        //update Hessian or recalculate it?
        if (recalculateHess) {ObjFunc.hEnergy(Hessian,HessianType,projecthessian,1,dampfactor,thresh_grad,maxiter);}
        else {
          delta_g = grad_new - grad_old;
          Solver.HessianUpdate(Hessian,sk,delta_g,alpha);
        }
      }
    }
    std::cout << "current energy = " << fnew << " Eh" << std::endl;
    if ((normH/thresh_step < 1.0 + thresh_step)&&(maxH/thresh_step < 1.0 + thresh_step)) {convergedS = true;}       //converged?
    else if ((normH/normG > 100.0)&&(normE/thresh_energy < 1.0 + thresh_energy)) {convergedS = true;}               //can't do better than this actually
    else {convergedS = false;}
    if (normE/thresh_energy < 1.0 + thresh_energy) {convergedE = true;}                                             //converged?
    else {convergedE = false;}
    if (convergencemode > 0) {std::cout << "        |x(k) - x(k-1)| = " << normH << "         max|x(k) - x(k-1)| = " << maxH << "        congH = " << convergedS << std::endl;}
    if ((convergencemode == 0)&&(convergedG)&&(convergedE)) {converged = true;}
    else if ((convergencemode == 1)&&(convergedG)&&(convergedS)) {converged = true;}
    else if ((convergencemode == 2)&&(convergedE)&&(convergedS)) {converged = true;}
    else if ((convergencemode == 3)&&(convergedG)&&(convergedE)&&(convergedS)) {converged = true;}
    if (converged) {
      std::cout << "Convergence in " << iter + 1 + (!TROK) << " iterations." << std::endl;
      break;
    }
    if (totalnrsteps > maxiter) {break;}
  }
  if (!converged) {std::cout << "WARNING: SolverPackage.hpp: SolverOpt(): no convergence in " << maxiter << " iterations"  << std::endl;}
  return converged;
}
template <class func> 
void OptimizeTS(func & ObjFunc, int HessianType = 0, int convergencemode = 0, double thresh_energy = 1.0e-7, double thresh_grad = 1.0e-7, int maxiter = 500, int macroiter = 200, double freqmin = -1.0) {
  //freqmin is the minimum value for the first positive frequency
  BakerRFO solve(1,1,0);
  bool converged = false;
  bool vibrthresh = (freqmin < 0.0);
  for (size_t iter = 0; iter < macroiter; ++iter) {
    converged = SolverOpt(ObjFunc,solve,HessianType,convergencemode,thresh_energy,thresh_grad,1.0e-7,maxiter);
    if (freqmin > 0.0) {
      std::vector<double> vibr = ObjFunc.CalcVibrFrequencies(false);
      for (size_t idmode = 0; idmode < vibr.size(); ++idmode) {
        if (vibr[idmode] > 0.5) {
          if ((vibr[idmode] >= freqmin)&&(vibr[0] < 0.0)) {vibrthresh = true;}
          break;
        }
      }
    }
    if ((vibrthresh)&&(converged)) {break;}
  }
}

#endif //_Solver_Package_
