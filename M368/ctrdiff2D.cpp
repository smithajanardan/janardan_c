/************************************************************
Function to implement the central-difference method to find 
an approximate solution of a hyperbolic IBVP in a rectangular 
domain of the form 

 utt + s ut = P uxx + Q uyy + p ux + q uy + r u + eta, 
                                   a<=x<=b, c<=y<=d, 0<=t<=T  

 u(a,y,t) = ga(y,t),   u(b,y,t) = gb(y,t),  c<=y<=d, 0<=t<=T  
 u(x,c,t) = gc(x,t),   u(x,d,t) = gd(x,t),  a<=x<=b, 0<=t<=T  

 u(x,y,0) = f(x,y),                         a<=x<=b, c<=y<=d  
 ut(x,y,0) = gamma(x,y),                    a<=x<=b, c<=y<=d  

Inputs:  
  PDEeval       Function to evaluate P,Q,p,q,r,s,eta
  BCeval        Function to evaluate ga,gb,gc,gd
  ICeval        Function to evaluate f,gamma
  a,b,c,d       Space domain parameters
  N,M           Number of interior x,y pts (N+2,M+2 total pts) 
  x             Grid point vector: x(i)=a+i*dx, i=0...N+1
  y             Grid point vector: y(j)=c+j*dy, j=0...M+1
  t,dt          Current time t and time step dt
  uold          Approx soln at time t-dt 
  u             Approx soln at time t:  u(i,j) is soln 
                at x(i),y(j), i=0...N+1, j=0...M+1

Outputs: 
  x             Grid point vector: x(i)=a+i*dx, i=0...N+1
  y             Grid point vector: y(j)=c+j*dy, j=0...M+1
  uold          Approx soln at time t
  u             Approx soln at time t+dt: u(i,j) is soln 
                at x(i),y(j), i=0...N+1, j=0...M+1

Note 1: This function is incomplete; you'll need to finish
coding the method as indicated below.  

Note 2: The functions PDEeval, BCeval and ICeval are assumed 
to be defined externally (e.g. by calling program).

Note 3: Rather than use the single-label index l=1...NM
for the interior xy-grid, it is convenient to use the 
double-label indices i=0...N+1, j=0...M+1 for the entire 
xy-grid.
************************************************************/
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;


/*** Ext fun: PDEeval(x,y,t,P,Q,p,q,r,s,eta) ***/
void PDEeval(const double&, const double&, double&, 
                     double&, double&, double&, double&, 
                                double&, double&, double&) ;

/*** Ext fun: BCeval(x,y,t,gLeft,gRight,gBottom,gTop) ***/
void BCeval(const double&, const double&, double&,
                     double&, double&, double&, double&) ;

/*** Ext fun: ICeval(x,y,f,gamma) ***/
void ICeval(const double&, const double&, double&, double&) ;


/*** Main function: central-difference method ***/
int ctrdiff2D(int N, int M, 
              double a, double b, double c, double d, 
              vector& x, vector& y, double& t, double& dt, 
                                   matrix& uold, matrix& u){

  int success_flag=0 ;
  double P, Q, p, q, r, s, eta, tn, tnn ;
  double gLeft, gRight, gBottom, gTop, f, gamma ;
  double uLeft=0, uRight=0, uBottom=0, uTop=0, uCenter=0 ;
  double dx=(b-a)/(N+1), dy=(d-c)/(M+1) ; //grid parameters
  double aij, bij, cij, dij, eij ; //fwd-diff parameters
  matrix unew(N+2,M+2) ; //temporary variable
  
  /*** Compute u^{n+1} at interior grid points ***/
  for(int i=1; i<=N; i++){ //actual i-index on grid
    for(int j=1; j<=M; j++){ //actual j-index on grid

      tn = t ;
      PDEeval(x(i),y(j),tn,P,Q,p,q,r,s,eta) ;
      BCeval(x(i),y(j),tn,gLeft,gRight,gBottom,gTop) ;
      dij = Q/(dy*dy) + q/(2.0*dy) ;
      eij = Q/(dy*dy) - q/(2.0*dy) ;
      cij = P/(dx*dx) + p/(2.0*dx) ;
      aij = P/(dx*dx) - p/(2.0*dx) ;
      bij = r - 2.0*P/(dx*dx) - 2.0*Q/(dy*dy) ;

      if( j<M ){ 
        uTop = u(i,j+1) ; 
      } 
      else {
        uTop = gTop ; 
      } 

      if( i>1 ){ 
        uLeft = u(i-1,j) ; 
      }
      else {
        uLeft = gLeft ; 
      }

      uCenter = u(i,j) ;

      if( i<N ){ 
        uRight = u(i+1,j) ; 
      }
      else {
        uRight = gRight ; 
      }

      if( j>1 ){ 
        uBottom = u(i,j-1) ; 
      }
      else {
        uBottom = gBottom; 
      }

      if(t==dt){
        /*** Compute starting value unew = u^{1} using 
        u = u^{0} = f and (du/dt)^{0} = gamma.  Upon 
        entry, the u-matrix has the f-values, but we 
        need to call ICeval to get the gamma-values.  
        This part is complete. ***/
        ICeval(x(i),y(j),f,gamma) ;
        unew(i,j) = u(i,j) 
                        + dt*gamma
                        + (dt*dt/2)*dij*uTop
                        + (dt*dt/2)*aij*uLeft 
                        + (dt*dt/2)*bij*uCenter 
                        + (dt*dt/2)*cij*uRight
                        + (dt*dt/2)*eij*uBottom
                        + (dt*dt/2)*eta 
                        - (dt*dt/2)*s*gamma ;
      } else {
        /*** Compute unew = u^{n+1} using u = u^{n} and 
        uold = u^{n-1}. THIS PART NEEDS TO BE COMPLETED.  ***/
        unew(i,j) = 0 ;
		for (int i = 0; i < N+2; i++) {
			for (int j=0; j < M+2; j++) {
				unew += 2*u(i,j)-(1-((dt/2)*s))*uold(i,j)+(pow(dt,2)*
			}
		}



      }

    }
  }
 
  /*** Compute u^{n+1} at boundary grid points ***/
  for(int i=0; i<=N+1; i++){ //actual i-index on grid
    tnn = t + dt ;
    BCeval(x(i),y(M+1),tnn,gLeft,gRight,gBottom,gTop) ;
    unew(i,M+1) = gTop ; 
    BCeval(x(i),y(0),tnn,gLeft,gRight,gBottom,gTop) ;
    unew(i,0) = gBottom ; 
  }
  for(int j=0; j<=M+1; j++){ //actual j-index on grid
    tnn = t + dt ;
    BCeval(x(N+1),y(j),tnn,gLeft,gRight,gBottom,gTop) ;
    unew(N+1,j) = gRight ; 
    BCeval(x(0),y(j),tnn,gLeft,gRight,gBottom,gTop) ;
    unew(0,j) = gLeft ; 
  }

  uold = u ; //update uold-matrix
  u = unew ; //update u-matrix

  return success_flag=1 ;
}

