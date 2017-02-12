/*******************************************************
Note 2: If you add another function, for example to 
implement the symmetric method, then remember to include 
a header for the function in the main program below, 
and to include the filename when compiling.  Also, 
remember to re-initialize x before running each method.
*******************************************************/
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;


/*** Declare user-defined functions to be used ***/
int genpower(matrix&, vector&, double&, int, int, double) ;


int main() {
  /*** Input problem data ***/
  int n=3, maxIter=50, iter;  
  double tol=1e-4, lambda; 
  matrix A(n,n);
  vector x(n);

  A(0,0)= 1.00; A(0,1)= 0.9379; A(0,2)=-0.1863; //A(0,3)= 0.15;
  A(1,0)= 1.00; A(1,1)= 1.7103; A(1,2)=-0.8691; //A(1,3)= 0.07;
  A(2,0)= 0.00; A(2,1)=-1.1254; A(2,2)=-0.3762; //A(2,3)= 1.19;
  //A(3,0)= 0.15; A(3,1)= 0.07; A(3,2)= 1.19; A(3,3)= 6.20;

  x(0)=1; x(1)=1; x(2)=-1; //x(3)=0; //initial vec


  /*** Print data to screen ***/
  cout << endl; 
  cout << "Given: A = " << endl;
  cout << A << endl;
  cout << "Given: x^(0) = " << endl;
  cout << x << endl;

  /*** Call general power function ***/
  iter=genpower(A,x,lambda,n,maxIter,tol);
  
  /*** Print results to screen ***/
  cout << endl; 
  cout << "Iteration index: k = " << iter << endl;
  cout << "Approximate eigval: lambda^(k) = " << lambda << endl;
  cout << "Approximate eigvec: x^(k) = " << endl;
  cout << x << endl ;

  return 0; //terminate main program
}

