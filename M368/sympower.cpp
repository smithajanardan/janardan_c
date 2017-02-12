#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;
int sympower(matrix& A, vector& x, double& lambda, 
                          int n, int maxIter, double tol) {
  int iter=0, pindex ;
  double error=1, ynorm ;
  vector y(n), r(n) ;
  y = x ;
  while(iter<maxIter && error>=tol) {
    ynorm = vecL2Norm(y) ;
    x = scaleVec(1/ynorm,y) ;
    y = matVecMult(A,x) ;
	lambda = vecDot(x,y) ;
    r = scaleVec(1/lambda,y) - x ;
    error = vecL2Norm(r) ;
//    ynorm = vecL2Norm(y) ;
//    x = scaleVec(1/ynorm,y) ;
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