#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;

/*** Declare user-defined functions to be used ***/
double poly(vector&, int, double); // defined below
int gauss_elim(matrix&, vector&, vector&);

int main() {
  /*** Input (x,y) data ***/
  int m=11;
  vector x(m), y(m);

  x(0) = 0.00 ;  y(0) = 0.2817 ;
  x(1) = 0.07 ;  y(1) = 1.0845 ;
  x(2) = 0.15 ;  y(2) = 1.5160 ;
  x(3) = 0.22 ;  y(3) = 2.0952 ;
  x(4) = 0.31 ;  y(4) = 2.4731 ;
  x(5) = 0.50 ;  y(5) = 2.7769 ;
  x(6) = 0.60 ;  y(6) = 2.8647 ;
  x(7) = 0.75 ;  y(7) = 2.9361 ;
  x(8) = 0.83 ;  y(8) = 2.8972 ;
  x(9) = 0.94 ;  y(9) = 2.9753 ;
  x(10) = 1.00 ;  y(10) = 3.0817 ;

  /*** Build least-squares arrays A,b ***/
  int n; 
  cout << "Enter polynomial degree: " << flush;
  cin >> n;
  if( n<0 || n>m-1 ){
    cerr << "Least-squares: bad degree -- exiting" << endl;
    exit(EXIT_FAILURE);
  }

  matrix A(n+1,n+1);
  vector b(n+1), alpha(n+1);
  A=0 ; b=0 ; // initialize for sum
  for (int i = 0; i <= n; i++) {
	  for (int j = 0; j <= n; j++) {
		  double sum = 0;
		  for (int l = 0; l < m; l++) {
			  sum += pow(x(l),i+j);
		  }
		  A(i,j) = sum;
	  }
	  double sum  = 0;
	  for (int k = 0; k < m; k++) {
		  sum += y(k)*pow(x(k),i);
	  }
	  b(i) = sum;
  }
  /*** Print least-squares arrays to screen ***/
  cout << endl;
  cout << "Least-squares degree: n = " << endl;
  cout << n << endl << endl;
  cout << "Least-squares array: A = " << endl;
  cout << A << endl;
  cout << "Least-squares array: b = " << endl;
  cout << b << endl;

  /*** Solve normal equations and print results ***/
  int success_flag;  
  success_flag=gauss_elim(A,b,alpha);
  cout << "Least-squares solution: alpha^(k) = " << endl;
  cout << alpha << endl;

  /*** Compute fitting error and print results ***/
  double E=0;
  for(int k=0; k<m; k++) {
     E = E + pow(y(k)-poly(alpha,n,x(k)),2);
  }
  cout << "Least-squares fitting error: E = " << endl;
  cout << E << endl;

  return 0; //terminate main program
}

/*** Auxiliary function to evalute poly P_n ***/
double poly(vector& alpha, int n, double x) {
  double P=0;
  for(int i=0; i<n+1; i++) {
    P = P + alpha(i)*pow(x,i); 
  }
  return P;
}