/*******************************************************************
Function to find dominant eigenvalue and eigenvector of a square 
matrix A using the general power method.  The eigenvalue and 
eigenvector are assumed to be real, and the matrix A is assumed 
to be diagonalizable.

Inputs: matrix A(n,n), vector x(n), double lambda,
        integer n, integer maxIter, double tol.

Outputs: vector x(n), double lambda, integer iter (returned)
********************************************************************/
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;


int genpower(matrix& A, vector& x, double& lambda, 
                          int n, int maxIter, double tol) {
  
  int iter=0, pindex ;
  double error=1, ynorm ;
  vector y(n), r(n) ;

  y = x ;
  while(iter<maxIter && error>=tol) {
    ynorm = vecMaxNorm(y) ;
    x = scaleVec(1/ynorm,y) ;
    pindex=0 ;
    for(int j=0; j<n; j++) {
      if(fabs(x(j))>fabs(x(pindex))) pindex=j;
    }
    y = matVecMult(A,x) ;
    lambda = y(pindex)/x(pindex) ;
    r = scaleVec(1/lambda,y) - x ;
    error = vecMaxNorm(r) ;
    iter++ ;
  }

  if(error < tol) {
    cout << "Gen power: soln converged, |r|_inf = " << error << endl;
  }
  else {
    cout << "Gen power: max iterations exceeded" << endl;
  }
  return iter;
}


