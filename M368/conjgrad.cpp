/*******************************************************************
Function to solve Ax=b using the pre-conditioned Conjugate 
Gradient method.  A is assumed to be symmetric, pos-definite.

Inputs: matrix Cinv(n,n), matrix A(n,n), vector b(n), vector x(n),
        integer maxIter, double tol.

Outputs: vector x(n), integer iter (returned)
********************************************************************/
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;


int conjgrad(matrix& Cinv, matrix& A, vector& b, vector& x, 
                                       int maxIter, double tol) {
  
  /*** Check input data ***/
  int n = dim(A,0);

  if(dim(A,1)!=n || dim(Cinv,0)!=n || dim(Cinv,1)!=n || 
                                     dim(b)!=n || dim(x)!=n) {
    cerr << "CG: bad data -- exiting" << endl; 
    exit(EXIT_FAILURE); 
  }
  if(tol <= 0) {
    cerr << "CG: bad tol -- exiting" << endl; 
    exit(EXIT_FAILURE);
  }
  if(maxIter <= 0) {
    maxIter = 1;
  }

  /*** Implement algo ***/
  int iter ;
  double error, alpha, beta, t, s;

  matrix CinvT(n,n);
  vector r(n), u(n), v(n), w(n);

  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) {
      CinvT(i,j) = Cinv(j,i); //transpose
    }
  }

  iter = 0;

  r = b - matVecMult(A,x);
  w = matVecMult(Cinv,r);
  v = matVecMult(CinvT,w);
  alpha = vecDot(w,w);
  error = vecMaxNorm(r);

  while(iter<maxIter && error>=tol) {
    u = matVecMult(A,v);
    t = alpha/vecDot(v,u);
    x = x + scaleVec(t,v);
    r = r - scaleVec(t,u);
    w = matVecMult(Cinv,r);
    beta = vecDot(w,w);
    s = beta/alpha;
    v = matVecMult(CinvT,w) + scaleVec(s,v);
    alpha = beta; 
    error = vecMaxNorm(r);
    iter++ ;
  }

  /*** Print exit message to screen ***/
  if(error<tol) {
    cout << "CG: solution converged, |r|_inf = " << error << endl;
  }
  else {
    cout << "CG: max iterations exceeded" << endl;
  }

  return iter;
}

