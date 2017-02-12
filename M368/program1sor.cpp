/****************************************************
Program 1.  Solves Ax=b using Jacobi method.  

Inputs: A, b, x^(0), maxIter, tol
Outputs: x^(k), iteration count k

Here's how to get started:

1) Copy program1.cpp (this file), jacobi.cpp, 
matrix.cpp and matrix.h into your working directory.

2) Compile (and link) the files by typing
"c++ -o program1 matrix.cpp jacobi.cpp program1.cpp"
at the Linux prompt.

3) Type "program1" to run the program.

4) For any given problem, you'll need to set the 
values of the variables {n,A,b,x,maxIter,tol} below, 
and then compile and run as described above.
****************************************************/
#include <iostream> // imports librarys and shit
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;

/*** Declare user-defined functions to be used ***/
int sor(matrix&, vector&, vector&, int, double, double);


/*** Main program ***/
int main() {

  /*** Define and input problem data ***/
  int n=100, maxIter=100000, iter;  
  double tol=1e-5;
  matrix A(n,n);
  vector x(n), b(n);
  double w=1.2;

  for (int i=0; i<n; i++) {
	  b(i)=1+(i/20);
	  for (int j=0; j<n; j++) {
		  if (j==i-1) {
			  A(i,j)=-1;
		  }
		  else if (j==i) {
			  A(i,j)=2+(1/10);
		  }
		  else if (j==i+1) {
			  A(i,j)=-1;
		  }
		  else {
			  A(i,j)=0;
		  }
	  }
  }


  //A(0,0)= 4; A(0,1)= 1; A(0,2)=-1; A(0,3)= 1; b(0)=-2;
  //A(1,0)= 1; A(1,1)= 4; A(1,2)=-1; A(1,3)=-1; b(1)=-1;
  //A(2,0)=-1; A(2,1)=-1; A(2,2)= 5; A(2,3)= 1; b(2)= 0;
  //A(3,0)= 1; A(3,1)=-1; A(3,2)= 1; A(3,3)= 3; b(3)= 1;
  x=0 ; //initial x


  /*** Print data to screen ***/
  cout << endl; 
  cout << "Given: A = " << endl;
  cout << A << endl;
  cout << "Given: b = " << endl;
  cout << b << endl;
  cout << "Given: x^(0) = " << endl;
  cout << x << endl;


  /*** Call Gauss Seide function ***/
  iter=sor(A,b,x,maxIter,tol,w);
  

  /*** Print results to screen ***/
  cout << endl; 
  cout << "Iteration index: k = " << iter << endl;
  cout << endl; 
  cout << "Approximate solution: x^(k) = " << endl;
  cout << x << endl;

  return 0;
}

