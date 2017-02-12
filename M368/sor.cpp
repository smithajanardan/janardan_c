/*******************************************************************
Function to solve Ax=b using the SOR method.

Inputs: matrix A(n,n), vector b(n), vector x(n),
        integer maxIter, double tol.

Outputs: vector x(n), integer iter (returned)
********************************************************************/
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;


int sor(matrix& A, vector& b, vector& x, 
                             int maxIter, double tol, double w) {
  
  /*** Check input data ***/
  int n = dim(A,0);

  if(dim(A,1) != n || dim(b) != n || dim(x) != n) {
    cerr << "SOR: bad data -- exiting" << endl; 
    exit(EXIT_FAILURE); 
  }
  if(tol <= 0) {
    cerr << "SOR: bad tol -- exiting" << endl;  
    exit(EXIT_FAILURE);
  }
  if(maxIter <= 0) {
    maxIter = 1;
  }
  for(int i=0; i<n; i++) {
    if(A(i,i) == 0) {
      cerr << "SOR: zero on diagonal -- exiting" << endl;
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
	if(j!=i) sum += A(i,j)*x(j);
      }
      x(i) = (w*(( -sum + b(i) ) / A(i,i))) + ((1-w)*x(i));
    }
    r = b - matVecMult(A,x);
    error = vecMaxNorm(r);
    xold = x;
    iter++;
  }

  /*** Print exit message to screen ***/
  if(error<tol) {
    cout << "SOR: solution converged, |r|_inf = " << error << endl;
  }
  else {
    cout << "SOR: max iterations exceeded" << endl;
  }

  return iter;
}


