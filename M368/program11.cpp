/************************************************************
Program 11.  Uses the central-difference method to find an 
approximate solution of an elliptic BVP in a rectangular 
domain of the form

 P uxx + Q uyy + p ux + q uy + r u = f,  a<=x<=b, c<=y<=d  
 u(a,y) = ga(y),     u(b,y) = gb(y),          c<=y<=d  
 u(x,c) = gc(x),     u(x,d) = gd(x),          a<=x<=b  

Inputs:  
  PDEeval       Function to evaluate P,Q,p,q,r,f
  BCeval        Function to evaluate ga,gb,gc,gd
  a,b,c,d       Domain parameters
  N,M           Number of interior x,y pts (N+2,M+2 total pts) 
  x             Grid point vector: x(i)=a+i*dx, i=0...N+1
  y             Grid point vector: y(j)=c+j*dy, j=0...M+1

Outputs: 
  x             Grid point vector: x(i)=a+i*dx, i=0...N+1
  y             Grid point vector: y(j)=c+j*dy, j=0...M+1
  u             Approx soln: u(i,j)=soln at x(i),y(j), 
                i=0...N+1, j=0...M+1

Note 1: The function file linearcd2D.cpp is incomplete; 
you'll need to code the A-matrix and G-vector as 
indicated in that file.

Note 2: For any given problem, the functions PDEeval 
and BCeval must be changed.

Note 3: For any given problem, the grid parameters a,b,
c,d,N,M must be specified.  

Note 4: Gauss elimination is used to solve the system,
so only moderate values of N,M should be used.
************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;
/*** Define output file ***/
const char myfile[50]="C:/M368/program11.out" ;
ofstream prt(myfile) ;
/*** Declare external function ***/
int linearcd2D(int, int, double, double, double, double,
                                      vector&, vector&, matrix&) ;
/*** Define P(x,y), Q(x,y), p(x,y), q(x,y), r(x,y), f(x,y) ***/
void PDEeval(const double& x, const double& y, 
             double& P, double& Q, double& p, double& q,
                                           double& r, double& f){
  double x0 = 0.2;
  double y0 = 0.2;
  double x1 = 0.8;
  double y1 = 0.2;
  P = 0.3 + 0.05*y ;
  Q = 0.3 + 0.05*y ;
  p = 10 - 10*x ;
  q = 5*y ;
  r = 0 ;
  f = 10*exp((-30*pow((x-x0),2))+(-30*pow((y-y0),2))) + 8*exp((-30*pow((x-x1),2))+(-30*pow((y-y1),2))) ; }
/*** Define ga(y), gb(y), gc(x), gd(x) ***/
void BCeval(const double& x, const double& y,
            double& ga, double& gb, double& gc, double& gd){
  ga = 0 ;
  gb = 0 ;
  gc = 0 ;
  gd = 0 ; }
int main() {
  /*** Define problem parameters ***/
  int N=29, M=29, success_flag=0 ;  
  matrix u(N+2,M+2) ;
  vector x(N+2), y(M+2) ; 
  double a=0, b=1, c=0, d=1 ; 
  double dx=(b-a)/(N+1), dy=(d-c)/(M+1) ;
  /*** Construct grid ***/
  for(int i=0; i<=N+1; i++){
    x(i) = a + i*dx ;  }
  for(int j=0; j<=M+1; j++){
    y(j) = c + j*dy ;  } 
  /*** Call central-difference method ***/
  success_flag=linearcd2D(N,M,a,b,c,d,x,y,u) ;
  /*** Print results to output file ***/
  prt.setf(ios::fixed) ;
  prt << setprecision(5) ;
  cout << "Linear-CD-2D: output written to " << myfile << endl ;
  prt << "Linear-CD-2D results" << endl ;
  prt << "Number of interior x-grid pts: N = " << N << endl ;
  prt << "Number of interior y-grid pts: M = " << M << endl ;
  prt << "Approximate solution: x_i, y_j, u_ij" << endl ;
  for(int i=0; i<=N+1; i++){
    for(int j=0; j<=M+1; j++){
      prt << setw(8) << x(i) ;
      prt << "   " ;
      prt << setw(8) << y(j) ;
      prt << "   " ;
      prt << setw(8) << u(i,j) ;
      prt << endl;    }  }
  return 0 ; } //terminate main program