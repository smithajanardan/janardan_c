#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;
/*** Declare external function ***/
int shootRK4(int&, double&, double&, double&, double&, 
                double&, int&, double&, vector&, vector&) ;
/*** Define f(x,y,yp) function ***/
void feval(const double& x, const double& y, 
                              const double& yp, double& f){
  double x1 = pow(x-1,2), y1 = pow(y-1,2);
  double x2 = pow(x+1,2), y2 = pow(y+1,2);
  double p = ((1-x)*exp(-x1-y1))-(2*(x+1)*exp(-x2-y2));
  double q = ((1-y)*exp(-x1-y1))-(2*(y+1)*exp(-x2-y2));
  double u = exp(-2*(pow(x,2)+pow(y,2)+2))*((((2*pow(x,2))+(-4*x)
	  +1)*exp(x2+y2))+(((4*pow(x,2))+(8*x)+2)*exp(x1+y1)));
  double w = exp(-2*(pow(x,2)+pow(y,2)+2))*((((2*pow(y,2))+(-4*y)
	  +1)*exp(x2+y2))+(((4*pow(y,2))+(8*y)+2)*exp(x1+y1)));
  double v = (2*(x-1)*(y-1)*exp(-x1-y1))+(4*(x+1)*(y+1)*exp(-x2-y2));
  f = (((p*yp)-q)*(u+(2*v*yp)+(w*pow(yp,2))))/(1+pow(p,2)+pow(q,2)) ; }
int main() {
  /*** Define problem parameters ***/
  double tol=1e-6 ;
  int N=20, maxIter=10, iter ;  
  double a=-3.0, b=3.0, alpha=-2.0, beta=2.0, t ; 
  vector x(N+1), y(N+1) ; 
  t= 1.2 ; // initial guess of slope 
  /*** Call Newton-RK4 method ***/
  cout << setprecision(8) ;
  iter=shootRK4(N,a,b,alpha,beta,t,maxIter,tol,x,y) ;
  /*** Print results to screen ***/
  cout << setprecision(4) ;
  cout << "Number of Newton iterations: " << iter << endl ;
  cout << "Approx solution: t = " << t << endl ;
  cout << "Approx solution: x_j, y_j =  " << endl ;
  for(int j=0; j<N+1; j++){
    cout << "x =" << setw(6) << x(j) ;
    cout << "   " ;
    cout << "y =" << setw(10) << y(j) ;
    cout << "   " ; 
    cout << endl; }
  return 0 ; //terminate main program
}