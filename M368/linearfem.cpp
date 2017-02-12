/************************************************************
Function to implement the piecewise linear finite-element 
method to find an approximate solution of a two-point BVP 
of the form
              -[p(x) y']' + q(x) y = f(x),  a<=x<=b  
               y(a)=alpha,  y(b)=beta

Inputs:  
  BVPeval       Function to evaluate p,q,f and g,g' (BC func)
  a,b           Interval params
  alpha,beta    Boundary value params
  N             Number of interior pts (N+2 total pts) 
  x             Grid points: x(j), j=0...N+1

Outputs: 
  x             Grid points: x(j), j=0...N+1
  y             Approx soln: y(j), j=0...N+1


Note 1: The function BVPeval is assumed to be defined 
externally (e.g. by calling program).  
************************************************************/
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;


/*** External function: gauss_elim(A,F,x) ***/
int gauss_elim(matrix&, vector&, vector&) ;


/*** External function: BVPeval(x,p,q,f,g,g') ***/
void BVPeval(const double&, double&, 
                       double&, double&, double&, double&);


/*** Auxiliary function: evaluate phi_i(xval),phi'_i(xval) ***/
void PHIeval(int N, vector& x, int i, 
                       double xval, double& phi, double& dphi){

  /*** i=0...N+1 is actual index of grid ***/
  phi = 0 ; dphi = 0 ;
  if(i>=1 && xval>=x(i-1) && i<=N && xval<=x(i+1)){ 
    if(xval<=x(i)){
      phi = (xval-x(i-1))/(x(i)-x(i-1)) ;
      dphi = 1/(x(i)-x(i-1)) ;
    }
    else if(xval>x(i)){ 
      phi = (x(i+1)-xval)/(x(i+1)-x(i)) ;
      dphi = -1/(x(i+1)-x(i)) ;
    }
  }
}
  

/*** Auxiliary function: build A and F for FE eqns ***/
void FEMeval(int N, vector& x, matrix& A, vector& F){

  double xQ, hQ ; 
  double phi_i, dphi_i ;   //phi_{i} and phi'_{i}
  double phi_im, dphi_im ; //phi_{i-1} and phi'_{i-1}
  double phi_ip, dphi_ip ; //phi_{i+1} and phi'_{i+1}
  double p, q, f, g, dg ;
  int I ;

  for(int i=1; i<=N; i++){ //i = actual grid index
    I=i-1 ; //I = array index for A, F
    if(i>1){
      A(I,I-1) = 0 ;
      /** Use midpt quad rule for A_{I,I-1} integral **/
      for(int l=0; l<=N; l++){ //l = actual grid index
        hQ = x(l+1)-x(l) ;
        xQ = 0.5*(x(l)+x(l+1)) ;
        PHIeval(N,x,i-1,xQ,phi_im,dphi_im) ;
        PHIeval(N,x,i,xQ,phi_i,dphi_i) ;
        BVPeval(xQ,p,q,f,g,dg) ;
        A(I,I-1) = A(I,I-1) + dphi_i*dphi_im*p*hQ + phi_i*phi_im*q*hQ ;
      }
    }
    F(I) = 0 ;
    A(I,I) = 0 ;
    /** Use midpt quad rule for F_{I} and A_{I,I} integrals **/
    for(int l=0; l<=N; l++){ //l = actual grid index
      hQ = x(l+1)-x(l) ;
      xQ = 0.5*(x(l)+x(l+1)) ;
      PHIeval(N,x,i,xQ,phi_i,dphi_i) ;
      BVPeval(xQ,p,q,f,g,dg) ;
      A(I,I) = A(I,I) + dphi_i*dphi_i*p*hQ + phi_i*phi_i*q*hQ ;
      F(I) = F(I) + phi_i*f*hQ - dphi_i*p*dg*hQ - phi_i*q*g*hQ ;
    }
    if(i<N){
      A(I,I+1) = 0 ;
      /** Use midpt quad rule for A_{I,I+1} integral **/
      for(int l=0; l<=N; l++){ //l = actual grid index
        hQ = x(l+1)-x(l) ;
        xQ = 0.5*(x(l)+x(l+1)) ;
        PHIeval(N,x,i,xQ,phi_i,dphi_i) ;
        PHIeval(N,x,i+1,xQ,phi_ip,dphi_ip) ;
        BVPeval(xQ,p,q,f,g,dg) ;
        A(I,I+1) = A(I,I+1) + dphi_i*dphi_ip*p*hQ + phi_i*phi_ip*q*hQ ;
      }
    }
  }
}


/*** Main function: linear FE method ***/
int linearfem(int N, vector& x, vector& y){
  vector F(N), cvec(N) ;
  matrix A(N,N);
  int success_flag=0 ;
  double p, q, f, g, dg ;

  for (int i=0; i < N; i++) {
	  A(i,i) = 0;
	  for (int j=0; j < N; j++) {
		  A(i,j) = 0;
	  }
  }

  FEMeval(N,x,A,F) ;
  cout << "Linear-FE: arrays assembled" << endl ;
  gauss_elim(A,F,cvec) ;
  cout << "Linear-FE: equations solved" << endl ;
  
  BVPeval(x(0),p,q,f,g,dg) ;
  y(0) = g ;
  for(int i=0; i<N; i++){
    BVPeval(x(i+1),p,q,f,g,dg) ;
    y(i+1) = cvec(i) + g ; //soln at nodes 
  }
  BVPeval(x(N+1),p,q,f,g,dg) ;
  y(N+1) = g ;

  return success_flag=1 ;
}