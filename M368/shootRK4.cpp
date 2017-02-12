/*******************************************************
Function to implement the Newton-RK4 shooting method to
find an approximate solution of a two-point BVP of the 
form
              y'' = f(x,y,y'),  a<=x<=b  
              y(a) = alpha,  y(b) = beta

Inputs:  
  feval         Function to evaluate f(x,y,y') 
  a,b           Interval params
  alpha,beta    Boundary value params
  t             Initial guess for y'(a)
  N             Number of grid points for RK4
  maxIter,tol   Control params for Newton
  eps           Finite-difference param for Newton

Outputs: 
  k       Number of Newton steps taken
  t       Approx value of y'(a) 
  x       Grid point vector: x(j)=a+jh, j=0...N
  y       Approx soln vector: y(j)=soln at x(j), j=0...N

Note 1: The function feval is assumed to be defined 
externally (e.g. by calling program).

Note 2: The parameter eps is locally-defined;  it is 
not passed as an argument.  It can be adjusted below 
if necessary.
*******************************************************/
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;

/*** Declare external function feval(x,y,yp,fval) ***/
void feval(const double&, const double&, const double&, double&) ;
/*** Define auxiliary function: RK4 solver for IVP ***/
void ivp_solve_RK4(int& N, double& a, double& b, double& alpha, 
                               double& t, vector& x, vector& y){
  vector k1(2), k2(2), k3(2), k4(2) ;
  double f, h=(b-a)/N ; 
  matrix w(2,N+1) ;

  x(0) = a ;
  for(int j=0; j<N; j++){
    x(j+1) = x(j) + h ;
  }

  w(0,0) = alpha ;
  w(1,0) = t ;
  y(0) = w(0,0) ;
  for(int j=0; j<N; j++){
    feval(x(j),w(0,j),w(1,j),f) ;
    k1(0) = h*w(1,j) ; k1(1) = h*f ;
    feval(x(j)+h/2,w(0,j)+k1(0)/2,w(1,j)+k1(1)/2,f) ;
    k2(0) = h*(w(1,j)+k1(1)/2) ; k2(1) = h*f ;
    feval(x(j)+h/2,w(0,j)+k2(0)/2,w(1,j)+k2(1)/2,f) ;
    k3(0) = h*(w(1,j)+k2(1)/2) ; k3(1) = h*f ;
    feval(x(j)+h,w(0,j)+k3(0),w(1,j)+k3(1),f) ;
    k4(0) = h*(w(1,j)+k3(1)) ; k4(1) = h*f ;
    w(0,j+1) = w(0,j) + (k1(0)+2*k2(0)+2*k3(0)+k4(0))/6 ;
    w(1,j+1) = w(1,j) + (k1(1)+2*k2(1)+2*k3(1)+k4(1))/6 ;
    y(j+1) = w(0,j+1) ;
   }
}
/*** Auxiliary function to evaluate F(t) ***/
void Feval(int& N, double& a, double& b, double& alpha, 
           double& beta, double& t, vector& x, 
                                      vector& y, double& F){
  ivp_solve_RK4(N,a,b,alpha,t,x,y) ;
  F = y(N) - beta ;
}
/*** Auxiliary function to evaluate F'(t) w/finite-difference ***/
void DFeval(int& N, double& a, double& b, double& alpha, 
            double& beta, double& t, vector& x, 
                                     vector& y, double& DF){
  double eps = 1e-05 ; 
  double Fplus, Fminus ;
  double tplus = t+eps, tminus = t-eps ;

  Feval(N,a,b,alpha,beta,tplus,x,y,Fplus) ;
  Feval(N,a,b,alpha,beta,tminus,x,y,Fminus) ;
  DF = (Fplus-Fminus)/(2.0*eps) ;
}
/*** Main function: Newton-RK4 shooting ***/
int shootRK4(int& N, double& a, double& b, double& alpha, 
             double& beta, double& t, int& maxIter, 
                         double& tol, vector& x, vector& y){

  int k=0  ;
  double error ;
  double F, DF, dt ;

  DFeval(N,a,b,alpha,beta,t,x,y,DF) ;
  Feval(N,a,b,alpha,beta,t,x,y,F) ;
  error = fabs(F) ;

  cout << endl ;
  cout << "Newton-RK4: initial t = " << t << endl ;
  cout << "Newton-RK4: |F| = " << error << endl;

  while(k<maxIter && error>=tol) {
    dt = -F/DF  ;
    t = t + dt ;
    DFeval(N,a,b,alpha,beta,t,x,y,DF) ;
    Feval(N,a,b,alpha,beta,t,x,y,F) ;
    error = fabs(F) ;
    k++ ;
    cout << "Newton-RK4: |F| = " << error << endl;
  }
  if(error < tol) {
    cout << "Newton-RK4: soln converged" << endl ;
    cout << endl ;
  }
  else {
    cout << "Newton-RK4: max iter exceeded" << endl;
    cout << endl ;
  }
  return k ;
}