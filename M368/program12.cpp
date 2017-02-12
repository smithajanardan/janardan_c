/************************************************************
Program 12.  Uses the forward-difference method to find an 
approximate solution of a parabolic IBVP in a rectangular 
domain of the form

 ut = P uxx + Q uyy + p ux + q uy + r u + eta, 
                                   a<=x<=b, c<=y<=d, 0<=t<=T  

 u(a,y,t) = ga(y,t),   u(b,y,t) = gb(y,t),  c<=y<=d, 0<=t<=T  
 u(x,c,t) = gc(x,t),   u(x,d,t) = gd(x,t),  a<=x<=b, 0<=t<=T  

 u(x,y,0) = f(x,y),                         a<=x<=b, c<=y<=d  


Inputs:  
  PDEeval       Function to evaluate P,Q,p,q,r,eta
  BCeval        Function to evaluate ga,gb,gc,gd
  ICeval        Function to evaluate f
  a,b,c,d       Space domain parameters
  N,M           Number of interior x,y pts (N+2,M+2 total pts) 
  x             Grid point vector: x(i)=a+i*dx, i=0...N+1
  y             Grid point vector: y(j)=c+j*dy, j=0...M+1
  T             Time interval parameter (final time) 
  L             Number of time steps

Outputs: 
  x             Grid point vector: x(i)=a+i*dx, i=0...N+1
  y             Grid point vector: y(j)=c+j*dy, j=0...M+1
  u             Approx soln at time t=T: u(i,j) is soln 
                at x(i),y(j), i=0...N+1, j=0...M+1


Note 1: The function file fwddiff2D.cpp is incomplete; you'll 
need to finish coding the method as indicated in that file.

Note 2: For any given problem, the functions PDEeval, BCeval 
and ICeval must be changed.

Note 3: For any given problem, the grid parameters a,b,c,d,T
and N,M,L must be specified.
************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;
/*** Define output file ***/
const char myfile[50]="C://M368/program12.out" ;
ofstream prt(myfile) ;
/*** Declare external function ***/
int fwddiff2D(int, int, double, double, double, double,
                   vector&, vector&, double&, double&, matrix&) ;
/*** Define P(x,y,t), Q(x,y,t), 
                     p(x,y,t), q(x,y,t), r(x,y,t), eta(x,y,t) ***/
void PDEeval(const double& x, const double& y, double& t,
                    double& P, double& Q, double& p, double& q,
                                          double& r, double& eta){
  vector g(3), a(3), b(3) ;
  a(0) = 0.2; b(0) = 0.2; g(0) = 0.6; a(1) = 0.6; b(1) = 0.2;
  g(1) = -0.7; a(2) = -0.5; b(2) = 0; g(2) = -0.6;
  P = 0.03 ;   Q = 0.03 ;   p = y ;   q = -x ;  r = 0 ; 
  for(int i=1; i< 4; i++){ 
	  eta += g(i-1)*exp((-30*pow(x-a(i-1),2))-(30*pow(y-b(i-1),2)));  } }
/*** Define ga(y,t), gb(y,t), gc(x,t), gd(x,t) ***/
void BCeval(const double& x, const double& y, double& t,
            double& ga, double& gb, double& gc, double& gd){
  ga = 1 ;   gb = 1 ;  gc = 1 ;   gd = 1 ;   }
/*** Define f(x,y) ***/
void ICeval(const double& x, const double& y, double& f){  f = 1 ;  }
int main() {
  /*** Define problem parameters ***/
  int N=(2*15)-1, M=(2*15)-1, L=56*40, success_flag=0 ;  
  matrix u(N+2,M+2) ;
  vector x(N+2), y(M+2) ; 
  double a=-1, b=1, c=-1, d=1 ; 
  double dx=(1.0/15.0), dy=(1.0/15.0) ;
  double T=40, dt=(1.0/56.0) ; 
  /*** Construct xy-grid ***/
  for(int i=0; i<=N+1; i++){    x(i) = a + i*dx ;  }
  for(int j=0; j<=M+1; j++){    y(j) = c + j*dy ;  } 
  /*** Load initial condition ***/
  double t=0 ;
  double gLeft, gRight, gBottom, gTop, f ;
  for(int i=0; i<=N+1; i++){ //actual i-index on grid
    for(int j=0; j<=M+1; j++){ //actual j-index on grid
      BCeval(x(i),y(j),t,gLeft,gRight,gBottom,gTop) ;
      ICeval(x(i),y(j),f) ;
      if( j==M+1 ){ u(i,j) = gTop ; }
      if( i==N+1 ){ u(i,j) = gRight ; }
      if( i==0 ){ u(i,j) = gLeft ; }
      if( j==0 ){ u(i,j) = gBottom ; }
      if((i>0)&&(i<N+1)&&(j>0)&&(j<M+1)){ u(i,j) = f ; }     }  }
  /*** Call forward-difference method ***/
  for(int n=0; n<L; n++){
    success_flag = fwddiff2D(N,M,a,b,c,d,x,y,t,dt,u) ;
    t = t + dt ;  }
  /*** Print results to output file ***/
  prt.setf(ios::fixed) ;
  prt << setprecision(5) ;
  cout << "Fwd-Diff-2D: output written to " << myfile << endl ;
  prt << "Fwd-Diff-2D results" << endl ;
  prt << "Number of interior x-grid pts: N = " << N << endl ;
  prt << "Number of interior y-grid pts: M = " << M << endl ;
  prt << "Number of time steps: L = " << L << endl ;
  prt << "Final time of simulation: t = " << t << endl ;
  prt << "Approximate solution at time t: x_i, y_j, u_ij" << endl ;
  double uavg;
  for(int i=0; i<=N+1; i++){
    for(int j=0; j<=M+1; j++){
      prt << setw(8) << x(i) ;
      prt << "   " ;
      prt << setw(8) << y(j) ;
      prt << "   " ;
      prt << setw(8) << u(i,j) ;
      prt << endl;
	  uavg = uavg + u(i,j);    }  }
  uavg = uavg/double((N+2)*(M+2));
  cout << uavg << endl;
  return 0 ;   }//terminate main program