/************************************************************
Program 10.  Uses the piecewise linear finite-element method 
to find an approximate solution of a two-point BVP of the 
form
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


Note 1: For any given problem, the function BVPeval must 
be changed.  This function computes p,q,f and g,g' for 
any given x.  The function g is the BC function 
defined by

        g(x) = alpha + (x-a)[(beta-alpha)/(b-a)].
        g'(x) = (beta-alpha)/(b-a).

Note 2: For any given problem, the values of N and x(j),
j=0...N+1 must be specified.  The default choice for the
grid points is x(j) = a + jh, where h=(b-a)/(N+1).

Note 5: The program output is written to a file.
************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;
/*** Define output file ***/
const char myfile[50]="C:/M368/program10.out" ;
ofstream prt(myfile) ;
/*** Declare external function ***/
int linearfem(int, vector&, vector&) ; 
/*** Define p(x), q(x), f(x), g(x), g'(x) ***/
void BVPeval(const double& x, double& p, double& q, 
                            double& f, double& g, double& dg){
  double pi=4.0*atan(1.0) ;
  double a=-2, b=2, alphaBC=0, betaBC=0, beta=25.37, gamma=495.5863575 ;
  g = alphaBC + (x-a)*((betaBC-alphaBC)/(b-a)) ;
  dg = (betaBC-alphaBC)/(b-a) ;
  if (x < 0) {
    p = 0.1; }
  else {
    p = 0.2; }
  if ((x <= 0.5) && (x >= -0.5)) {
    q = 0;  }
  else {
	q = beta;  }
  if ((x <= 0.4) && (x >= -0.4)) {
	f = gamma;  }
  else {
	f = 0;  }}
int main() {
  /*** Define problem parameters ***/
  int N=199, success_flag ;  
  vector x(N+2), y(N+2) ;
  double a=-2, b=2, h=(b-a)/(N+1) ;
  /*** Define FE grid pts  ***/
  for(int j=0; j<=N+1; j++){
	  x(j) = a + j*h ;  }  //default
  /*** Call linear FE method ***/
  success_flag=linearfem(N,x,y) ;
  /*** Print results to output file ***/
  prt.setf(ios::fixed) ;
  prt << setprecision(5) ;
  cout << "Linear-FE: output written to " << myfile << endl ;
  prt << "Linear-FE results" << endl ;
  prt << "Number of interior grid pts: N = " << N << endl ;
  prt << "Approximate solution: x_j, y_j" << endl ;
  for(int j=0; j<=N+1; j++){
    prt << setw(8) << x(j) ;
    prt << "   " ;
    prt << setw(8) << y(j) ;
    prt << "   " ;
    prt << endl;  }
  return 0 ; }//terminate main program