#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std ;
/*** Declare function with steepest descent algorithm ***/
int descent(vector&, int, int, double) ;
/*** Define g function for problem ***/
void geval(vector& x, double& g){
  double a=1, b=-0.1 ; //define any constants
  double pi=4.0*atan(1.0) ; //the number pi
  vector p(5) ;
  vector v(5) ;
  p(0) = 2.50; v(0) = 0.024;
  p(1) = 5.00; v(1) = 0.036;
  p(2) = 10.0; v(2) = 0.053;
  p(3) = 15.0; v(3) = 0.060;
  p(4) = 20.0; v(4) = 0.064;
  g = 0;
  for (int i = 0; i < 5; i++) {
	  g += pow((v(i) - (x(0)*p(i))/(x(1)+p(i))),2);  }}
/*** Define dg function (gradient) for problem ***/
void dgeval(vector& x, vector& dg){
  double a=1, b=-0.1 ; //define any constants
  double pi=4.0*atan(1.0) ; //the number pi
  vector p(5) ;
  vector v(5) ;
//  matrix J(5,5);
//  vector F(5);
  p(0) = 2.50; v(0) = 0.024;
  p(1) = 5.00; v(1) = 0.036;
  p(2) = 10.0; v(2) = 0.053;
  p(3) = 15.0; v(3) = 0.060;
  p(4) = 20.0; v(4) = 0.064;
  dg(0) = 0;
  dg(1) = 0;
  for (int i = 0; i < 5; i++) {
	  dg(0) += (-p(i)/(x(1)+p(i)))*2*(v(i)-(x(0)*p(i)/(x(1)+p(i))));
	  dg(1) += (x(0)*p(i)/pow(x(1)+p(i),2))*2*(v(i)-(x(0)*p(i)/(x(1)+p(i)))); 
  }}
int main() {
  /*** Define problem parameters ***/
  int n=2, maxIter=100, iter=0 ;
  double tol=1e-6 ;
  vector x(n) ;
  x(0) = 0.907 ; x(1) = 7.533 ; //initial guess
  /*** Print problem info ***/
  cout << setprecision(10) ;
  cout << endl ;
  cout << "System size: n = " << n << endl ;
  cout << "Initial guess: x^(0) = " << endl ;
  cout << x << endl ;
  /*** Call steepest descent w/initial guess***/
  iter = descent(x,n,maxIter,tol) ;
  /*** Print results to screen ***/
  cout << endl ;
  cout << "Iteration index: k = " << iter << endl ;
  cout << "Approx solution: x = " << endl ;
  cout << x << endl ;
  return 0; //terminate main program
}