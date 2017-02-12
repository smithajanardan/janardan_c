#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;

/*** Specify name of input data file ***/
//const int maxChar = 30;
//const char myfile[maxChar]="c:/M368/program4.dat" ;

int main() {

  /*** Input/read (x_j,y_j) data ***/
  double m=3;
  double pi=4.0*atan(1.0); // the number pi
  vector x(2*m), y(2*m);

//  ifstream fileRead (myfile);
//  for(int j=0; j<2*m; j++) {
//    fileRead >> x(j) >> y(j); // read xy-pair 
//  }
	for (int i = 0; i < 6; i++) {
		x(i)= -pi + ((i/m)*pi);
		y(i)= (2*pow(x(i),2))-9;
	}
  /*** Echo input file details ***/
  cout << endl;
//  cout << "xy-data read from file: " << myfile << endl;
  cout << endl;
  cout << "Number of data points 2m: " << 2*m << endl;
  cout << endl;

  /*** Compute a,b-coefficients ***/
  int n; 
  cout << "Enter trig poly degree n: " << flush;
  cin >> n;
  if( n<2 || n>= m){
    cerr << "Trig LS: bad degree -- exiting" << endl;
    exit(EXIT_FAILURE);
  }

  vector a(n+1), b(n-1);
  a = 0 ; b = 0 ; // initialize for sum

  for (int i = 0; i < n+1; i++) {
	  double sum = 0.0;
	  for (int j = 0; j < 2*m; j++) {
		  sum += y(j)*cos(i*x(j));
	  }
	  a(i) = sum/m;
  }

  for (int i = 1; i < n; i++) {
	  double sum = 0.0;
	  for (int j = 0; j < 2*m; j++) {
		  sum += y(j)*sin(i*x(j));
	  }
	  b(i-1) = sum/m;
  }

  /*** Print least-squares coeffs to screen ***/
  cout << endl;
  cout << "Trig LS results for n = " << n << ", m = " << m << endl;
  cout << endl;
  cout << "Cosine coeffs: a_0 ... a_n = " << endl;
  cout << a << endl;
  cout << "Sine coeffs: b_1 ... b_{n-1} = " << endl;
  cout << b << endl;

  /*** Compute least-squares fitting error ***/
  double E = 0;
  for (int i = 0; i < m; i++) {
	  double sm = a(0)/2;
	  for (int j = 1; j < n+1; j++) {
		  sm += a(j)*cos(j*x(i));
	  }
	  for (int k = 0; k < n-1; k++) {
		  sm += b(k)*sin((k+1)*x(i));
	  }
	  E += pow(y(i)-sm,2.0);
  }
  cout << endl;
  cout << "Fitting error: E = " << endl;
  cout << E << endl;
  return 0;
}