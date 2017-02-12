#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std ;
/*** Declare user-defined functions to be used ***/
void geval(vector&, double&) ;  //defined in program8.cpp
void dgeval(vector&, vector&) ; //defined in program8.cpp
int descent(vector& xk, int n, int maxIter, double tol){
  int k=0, kSrch=0, maxSrch=20 ;  
  double a0= 1, gk, dgkNorm ;
  double a1, a2, a3, y1, y2, y3, ah, gh, ahh, ghh ;
  vector dgk(n), xh(n), xhh(n) ;
  geval(xk,gk) ; //evaluate g(xk)
  dgeval(xk,dgk) ;  //evaluate dg(xk)
  dgkNorm = vecL2Norm(dgk) ;
  while(dgkNorm>=tol && k<maxIter) {
    cout << "fb js ka" <<k << endl;
    ah = a0 ;
    xh = xk - scaleVec(ah/dgkNorm,dgk) ;
    geval(xh,gh) ;
    kSrch = 0 ;
    while( gh>=gk && kSrch<maxSrch ){ 
      ah = ah/2 ;
      xh = xk - scaleVec(ah/dgkNorm,dgk) ;
      geval(xh,gh) ;
      kSrch++ ;    }
    if(kSrch>=maxSrch){
      cout << "Descent: max iter exceeded in alpha search" << endl;    }
    a1 = 0 ; y1 = gk ;
    a2 = ah/2 ; y2 = 0 ;
    a3 = ah ; y3 = gh ;
    ahh = ah ; ghh = gh ; xhh = xh ;
	xhh = xk - scaleVec(a2/dgkNorm,dgk) ;
	geval(xhh,y2) ;
	double h1 = (y2-y1)/a2 ;
	double h2 = (y3-y2)/(a3-a2) ;
	double h3 = (h2-h1)/a3 ;
	ahh = (a2-h1)/(2*h3) ;
    xhh = xk - scaleVec(ahh/dgkNorm,dgk) ;
	geval(xhh,ghh) ;
    if(ahh>=0 && ahh<=ah && ghh<gh){
      xk = xhh ;
      gk = ghh ;    }
    else {
      xk = xh ;
      gk = gh ;    }
    dgeval(xk,dgk) ;
	cout << dgk << endl;
    dgkNorm = vecL2Norm(dgk) ; 
    k++ ;
    cout << "Descent: |Gradg_k|_2 = " << dgkNorm << endl ;
  }
  if(dgkNorm < tol) {
    cout << "Descent: solution converged" << endl ;  } 
  else {
    cout << "Descent: max iterations exceeded" << endl; }
  return k; }