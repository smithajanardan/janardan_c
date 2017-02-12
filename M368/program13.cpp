/************************************************************
Program 13.  Uses the central-difference method to find an 
approximate solution of a hyperbolic IBVP in a rectangular 
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
  T             Time interval parameter (final time) 
  L             Number of time steps

Outputs: 
  x             Grid point vector: x(i)=a+i*dx, i=0...N+1
  y             Grid point vector: y(j)=c+j*dy, j=0...M+1
  u             Approx soln at time t=T: u(i,j) is soln 
                at x(i),y(j), i=0...N+1, j=0...M+1

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
const char myfile[20]="program13.out" ;
ofstream prt(myfile) ;


/*** Declare external function ***/
int ctrdiff2D(int, int, double, double, double, double,
              vector&, vector&, double&, double&, matrix&, matrix&) ;

/*** Define P(x,y,t), Q(x,y,t), 
             p(x,y,t), q(x,y,t), r(x,y,t), s(x,y,t), eta(x,y,t) ***/
void PDEeval(const double& x, const double& y, double& t,
                    double& P, double& Q, double& p, double& q,
                                 double& r, double& s, double& eta){

  P = 0.1*(1 + exp(-x*y*t/10)) ;
  Q = 0.1 ;
  p = 0.02*x*y ;
  q = 0.04*x ;
  r = -0.01 ;
  s = 0.05 ;
  eta = 0.01*t - 0.02*x*sin(y) ;

}

/*** Define ga(y,t), gb(y,t), gc(x,t), gd(x,t) ***/
void BCeval(const double& x, const double& y, double& t,
            double& ga, double& gb, double& gc, double& gd){

  ga = 1 ; 
  gb = 1 ;
  gc = 1 ; 
  gd = 1 ;  

}

/*** Define f(x,y), gamma(x,y) ***/
void ICeval(const double& x, const double& y, 
                                double& f, double& gamma){

  f = 1 ;
  gamma = 0 ;

}


int main() {
  /*** Define problem parameters ***/
  int N=9, M=9, L=10*2, success_flag=0 ;  
  matrix u(N+2,M+2), uold(N+2,M+2) ;
  vector x(N+2), y(M+2) ; 
  double a=0, b=1, c=0, d=1 ; 
  double dx=(b-a)/(N+1), dy=(d-c)/(M+1) ;
  double T=2, dt=T/L ; 

  /*** Construct xy-grid ***/
  for(int i=0; i<=N+1; i++){
    x(i) = a + i*dx ;
  }
  for(int j=0; j<=M+1; j++){
    y(j) = c + j*dy ;
  } 

  /*** Load initial values into u^{0} at time t^{0} ***/
  double t=0 ;
  double gLeft, gRight, gBottom, gTop, f, gamma ;
  for(int i=0; i<=N+1; i++){ //actual i-index on grid
    for(int j=0; j<=M+1; j++){ //actual j-index on grid
      BCeval(x(i),y(j),t,gLeft,gRight,gBottom,gTop) ;
      ICeval(x(i),y(j),f,gamma) ;
      if( j==M+1 ){ uold(i,j) = gTop ; }
      if( i==N+1 ){ uold(i,j) = gRight ; }
      if( i==0 ){ uold(i,j) = gLeft ; }
      if( j==0 ){ uold(i,j) = gBottom ; }
      if((i>0)&&(i<N+1)&&(j>0)&&(j<M+1)){ uold(i,j) = f ; }
    }
  }

  /*** Load temporary values into u^{1} at time t^{1} ***/
  u = uold ;
  t = dt ;

  /*** Call central-difference method ***/
  for(int n=1; n<L; n++){
    success_flag = ctrdiff2D(N,M,a,b,c,d,x,y,t,dt,uold,u) ;
    t = t + dt ;
  }

  /*** Print results to output file ***/
  prt.setf(ios::fixed) ;
  prt << setprecision(5) ;
  cout << "Ctr-Diff-2D: output written to " << myfile << endl ;
  prt << "Ctr-Diff-2D results" << endl ;
  prt << "Number of interior x-grid pts: N = " << N << endl ;
  prt << "Number of interior y-grid pts: M = " << M << endl ;
  prt << "Number of time steps: L = " << L << endl ;
  prt << "Final time of simulation: t = " << t << endl ;
  prt << "Approximate solution at time t: x_i, y_j, u_ij" << endl ;
  for(int i=0; i<=N+1; i++){
    for(int j=0; j<=M+1; j++){
      prt << setw(8) << x(i) ;
      prt << "   " ;
      prt << setw(8) << y(j) ;
      prt << "   " ;
      prt << setw(8) << u(i,j) ;
      prt << endl;
    }
  }

  return 0 ; //terminate main program
}

