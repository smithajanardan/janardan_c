/*******************************************************************
Function to solve Ax=b using the Jacobi method.

Inputs: matrix A(n,n), vector b(n), vector x(n),
        integer maxIter, double tol.

Outputs: vector x(n), integer iter (returned)
********************************************************************/
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;


int jacobi(matrix& A, vector& b, vector& x, 
                             int maxIter, double tol) {
  
  /*** Check input data ***/
  int n = dim(A,0);

  if(dim(A,1) != n || dim(b) != n || dim(x) != n) {
    cerr << "Jacobi: bad data -- exiting" << endl; 
    exit(EXIT_FAILURE); 
  }
  if(tol <= 0) {
    cerr << "Jacobi: bad tol -- exiting" << endl;  
    exit(EXIT_FAILURE);
  }
  if(maxIter <= 0) {
    maxIter = 1;
  }
  for(int i=0; i<n; i++) {
    if(A(i,i) == 0) {
      cerr << "Jacobi: zero on diagonal -- exiting" << endl;
      exit(EXIT_FAILURE);
    }
  }

  /*** Implement algorithm ***/
  int iter;
  double error;
  vector xold(n), r(n);

  iter = 0;
  xold = x;
  r = b - matVecMult(A,x);
  error = vecMaxNorm(r);

  while(iter<maxIter && error>=tol) {
    for(int i=0; i<n; i++) {
      double sum = 0;
      for(int j=0; j<n; j++) {
	if(j!=i) sum += A(i,j)*xold(j);
      }
      x(i) = ( -sum + b(i) ) / A(i,i);
    }
    r = b - matVecMult(A,x);
    error = vecMaxNorm(r);
    xold = x;
    iter++;
  }

  /*** Print exit message to screen ***/
  if(error<tol) {
    cout << "Jacobi: solution converged, |r|_inf = " << error << endl;
  }
  else {
    cout << "Jacobi: max iterations exceeded" << endl;
  }

  return iter;
}


