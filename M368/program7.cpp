#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std ;
/*** Declare function with Newton algorithm ***/
int newton(vector&, int, int, double) ;
/*** Define F function for problem ***/
void Feval(vector& x, vector& F){
//  double a=3, b=1.5, c=1, d=2, h=3.5, r=2.5 ; //define any constants
//  double pi=4.0*atan(1.0) ; //the number pi
  F(0) = 16280.3 + ((10973463.11)/pow(5-x(0),2))-x(1);
  F(1) = 19745.1 + ((10973463.11)/pow(6-x(0),2))-x(1);
//  F(2) = 21320.2 + ((10973463.11)/pow(7-x(0),2))-x(1);
//  F(0) = (pow(x(2),2)/(a*a))+(pow(x(3),2)/(b*b))-1;
//  F(1) = (pow(x(4),2)/(c*c))+(pow(x(5)+h,2)/(d*d))-1; ;
//  F(2) = pow(x(2)-x(0),2)+pow(x(3)-x(1),2)-pow(r,2);
//  F(3) = pow(x(4)-x(0),2)+pow(x(5)-x(1),2)-pow(r,2);
//  F(4) = (x(2)*(x(3)-x(1))*(b*b)) - (x(3)*(x(2)-x(0))*(a*a));
//  F(5) = (x(4)*(x(5)-x(1))*(d*d)) - ((x(5)+h)*(x(4)-x(0))*(c*c));
}
/*** Define DF function (Jacobian) for problem ***/
void DFeval(vector& x, matrix& DF){
//  double a=3, b=1.5, c=1, d=2, h=3.5, r=2.5; //define any constants
//  double pi=4.0*atan(1.0) ; //the number pi
  DF(0,0) = (2*10973463.11)/pow(5-x(0),3);
  DF(0,1) = 1;
//  DF(0,2) = 0;
  DF(1,0) = (2*10973463.11)/pow(6-x(0),3);
  DF(1,1) = 1;
  /*
  DF(1,2) = 0;
  DF(2,0) = (2*10973463.11)/pow(7-x(0),3);
  DF(2,1) = 1;
  DF(2,2) = 0;
  DF(0,0) = 0 ; DF(1,0) =  0 ;
  DF(2,0) =  2.0*x(0)-2.0*x(2) ; DF(3,0) = 2.0*x(0)-2.0*x(4) ;
  DF(4,0) = a*a*x(3) ; DF(5,0) = c*c*(x(5)+h) ;
  DF(0,1) = 0 ; DF(1,1) =  0 ;
  DF(2,1) =  2.0*x(1)-2.0*x(3) ; DF(3,1) = 2.0*x(1)-2.0*x(5) ;
  DF(4,1) = -b*b*x(2) ; DF(5,1) = -d*d*x(4) ;
  DF(0,2) = 2*x(2)/(a*a) ; DF(1,2) =  0 ;
  DF(2,2) =  2.0*x(2)-2.0*x(0) ; DF(3,2) = 0 ;
  DF(4,2) = (x(3)-x(1))*(b*b) - a*a*x(3) ; DF(5,2) = 0 ;
  DF(0,3) = 2*x(3)/(b*b) ; DF(1,3) =  0 ;
  DF(2,3) =  2.0*x(3)-2.0*x(1) ; DF(3,3) = 0 ;
  DF(4,3) = b*b*x(2) - (x(2)-x(0))*a*a ; DF(5,3) = 0 ;
  DF(0,4) = 0 ; DF(1,4) = 2*x(4)/(c*c) ;
  DF(2,4) =  0 ; DF(3,4) = 2.0*x(4)-2.0*x(0) ;
  DF(4,4) = 0 ; DF(5,4) = d*d*(x(5)-x(1)) - c*c*(x(5)+h) ;
  DF(0,5) = 0 ; DF(1,5) =  2*(x(4)+h)/(d*d) ;
  DF(2,5) =  0 ; DF(3,5) = 2.0*x(5)-2.0*x(1) ;
  DF(4,5) = 0 ; DF(5,5) = d*d*x(5) - c*c*(x(4)-x(0)) ;
  */
}
int main() {
  /*** Define problem parameters ***/
  int n=2, maxIter=50, iter=0 ;  
  double tol=1e-6 ;
  vector x(n) ;
  x(0) = 0 ; x(1) = 20669 ; // x(2) = 0; 
//  x(2) = 0.00000001 ; x(3) = -1.5 ;
//  x(4) = 0.00000001 ; x(5) = -1.5 ;//initial guess
  /*** Print problem info ***/
  cout << setprecision(10) ;
  cout << endl ;
  cout << "System size: n = " << n << endl ;
  cout << "Initial guess: x^(0) = " << endl ;
  cout << x << endl ;
  /*** Call Newton w/initial guess***/
  iter = newton(x,n,maxIter,tol) ;
  /*** Print results to screen ***/
  cout << endl ;
  cout << "Iteration index: k = " << iter << endl ;
  cout << "Approx solution: x = " << endl ;
  cout << x << endl ;
  return 0; //terminate main program
}