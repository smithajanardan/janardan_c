#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;
int qr(matrix& A, matrix& Q, matrix& D, 
                 int n, int maxIter, double tol) {
  int iter=0 ;
  double error=1, c, s, theta ;
  double pi=4.0*atan(1.0) ; //the number pi
  matrix B(n,n), P(n,n), Pt(n,n), G(n,n), R(n,n), Eye(n,n) ;
  Eye = 0 ;
  for(int i=0; i<n; i++) {
    Eye(i,i) = 1 ; //the identity matrix of size n
  }
  B = A ; //work with copy of A to preserve original matrix
  D = 0 ;
  Q = Eye ;
  while(iter<maxIter && error>=tol) {
    G = Eye ;
    for(int i=2; i<n; i++){
      P  = Eye ; 
      Pt = Eye ;
	  double theta = atan (B(i,i-1)/B(i-1,i-1));
	  P(i,i) = cos(theta); Pt(i,i) = cos(theta);
	  P(i, i-1) = -sin(theta); Pt(i, i-1) = sin(theta);
	  P(i-1, i) = sin(theta); Pt(i-1, i) = -sin(theta);
	  P(i-1, i-1) = cos(theta); Pt(i-1, i-1) = cos(theta);
      B = matMatMult(P,B) ; //create zero at {i,i-1} 
      G = matMatMult(G,Pt) ; //accumulate Pt 
    }
    B = matMatMult(B,G) ; //finish similarity transformation
    Q = matMatMult(Q,G) ; //update current eigenvector matrix
    for(int i=0; i<n; i++) {
      D(i,i) = B(i,i) ; //extract current diagonal matrix
    }
    R = B - D ; //form residual of off-diagonal portion
    error = matMaxNorm(R) ;
    iter++ ;
  }
  if(error < tol) {
    cout << "QR method: soln converged, |R|_inf = " << error << endl;
  }
  else {
    cout << "QR method: max iterations exceeded" << endl;
  }
  return iter;
}