/*******************************************************************
Function to solve a general, square, linear system Ax=b using
the Gauss Elimination method with partial pivoting.

Inputs: matrix A(n,n), vector b(n), vector x(n)

Outputs: vector x(n), integer success_flag=1 (returned)

Note: The matrix A and vector b are not altered by this function
********************************************************************/
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;


/*** Define auxiliary function: pivot search ***/
int pivot_index(matrix& A, int k, int n){
  int i=k ;
  while( (i<n)&&(A(i,k)==0.0) ){
    i++ ;
  }
  return i ;
}

/*** Define auxiliary function: row exchange ***/
void row_exchange(matrix& A, vector& b, int k, int p, int n){
  double temp ;
  temp = b(k) ;
  b(k) = b(p) ;
  b(p) = temp ;
  for(int j=0; j<n; j++){
    temp = A(k,j) ;
    A(k,j) = A(p,j) ;
    A(p,j) = temp ;
  }
}

/*** Define auxiliary function: fwd elimination ***/
void fwd_elimination(matrix& A, vector& b, int k, int n){
  double mik ;
  for(int i=k+1; i<n; i++){
    mik = A(i,k)/A(k,k) ;
    b(i) = b(i) - mik*b(k) ;
    for(int j=k+1; j<n; j++){
      A(i,j) = A(i,j) - mik*A(k,j) ;
    }
  }
}

/*** Define auxiliary function: bwd substitution ***/
void bwd_substitution(matrix& A, vector& b, vector& x, int n){
  double temp ;
  x(n-1) = b(n-1)/A(n-1,n-1) ;
  for(int k=n-2; k>=0; k--){
    temp = 0.0 ; 
    for(int j=k+1; j<n; j++){
      temp = temp + A(k,j)*x(j) ;
    }
    x(k) = (b(k)-temp)/A(k,k) ;
  }
}

/*** Define main gauss elimination function ***/
int gauss_elim(matrix& A, vector& b, vector& x){

  /*** Check input data ***/
  int n=dim(A,0) ;
  if(dim(A,1) != n || dim(b) != n || dim(x) != n) {
    cerr << "Gauss Elim: bad data -- exiting" << endl;
    exit(EXIT_FAILURE);
  }

  matrix AA(n,n) ; 
  vector bb(n) ; 
  int p, success_flag=1 ;

  AA=A ; //copy of A-matrix (copy will be overwritten)
  bb=b ; //copy of b-vector (copy will be overwritten)

  for(int k=0; k<n-1; k++){
    p = pivot_index(AA,k,n) ;
    if( (p>k)&&(p<n) ){
      row_exchange(AA,bb,k,p,n) ;
    }
    if( p>=n ){
      cout << "Gauss Elim: no unique solution" << endl ;
      exit(EXIT_FAILURE) ;
    }
    fwd_elimination(AA,bb,k,n) ;
  }

  if( AA(n-1,n-1)==0.0 ){
    cout << "Gauss Elim: no unique solution" << endl ;
    exit(EXIT_FAILURE) ;
  }
  bwd_substitution(AA,bb,x,n) ;
  return success_flag ;
}
