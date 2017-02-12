/*******************************************************
Function to implement the Newton method to solve a 
general system of equations F(x) = 0.

Inputs:  
  Feval    Function to evaluate F(x), n vector 
  DFeval   Function to evaluate dF/dx(x), nxn matrix
  x        Initial guess, n vector  
  maxIter  Maximum number of iterations
  tol      Tolerance parameter

Outputs: 
  x        Approx soln of F(x)=0, n vector
  k        Number of Newton steps taken


Note 1: The functions Feval and DFeval are assumed 
to be defined externally (e.g. by calling program).

Note 2: Gauss elimination with partial pivoting is 
used to implement the (slightly modified) Newton 
update equations:  
                  DF^(k) dx^(k) = F^(k), 
                 x^(k+1) = x^(k) - dx^(k). 

The minus sign in the update formula has been put in 
a more convenient location.
*******************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std ;

/*** Declare user-defined functions to be used ***/
void Feval(vector&, vector&) ; //defined in program7.cpp
void DFeval(vector&, matrix&) ; //defined in program7.cpp
int gauss_elim(matrix&, vector&, vector&) ;

int newton(vector& x, int n, int maxIter, double tol){
  int k=0 ;  
  double error=1 ;
  matrix DF(n,n) ;
  vector F(n), dx(n) ;
  Feval(x,F) ;  //evaluate F(x)
  DFeval(x,DF) ; //evaluate DF(x)

  while(error>=tol && k<maxIter) {
    gauss_elim(DF,F,dx) ; // solve for dx
    x = x - dx ; 
    error = vecMaxNorm(dx) ; // error=||xnew-xold||
    Feval(x,F) ;
    DFeval(x,DF) ;
    k++ ;
    cout << "Newton: |xnew-xold|_inf = " << error << endl ;
  }

  if(error < tol) {
    cout << "Newton: solution converged" << endl ;
  }
  else {
    cout << "Newton: max iterations exceeded" << endl ;
  }
  return k ;
}

