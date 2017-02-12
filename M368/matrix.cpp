/*******************************************************
General purpose library: implementation file

General purpose library of vector and matrix data 
types and related functions.  See file matrix.h for 
description and usage notes.
********************************************************/

#include <math.h>
#include "matrix.h"

/****************** 
define matrix class 
*******************/
matrix::matrix(){ //default constructor 
  size[0] = 0;
  size[1] = 0;
  array = 0;
}
matrix::matrix(int m, int n){ //constructor 
  size[0] = m;
  size[1] = n;
  array = new double[size[0]*size[1]];
}
matrix::~matrix(){ //destructor 
  if(array != 0) delete []array;
}
matrix::matrix(const matrix& A) { //copy constructor
  size[0] = A.size[0];
  size[1] = A.size[1];
  array = new double[size[0]*size[1]];
  for(int i=0; i<size[0]*size[1]; i++) array[i] = A.array[i];
}
double matrix::operator() (int i, int j) const { //get-value index operator
  if(i<0 || i>=size[0] || j<0 || j>=size[1])
    cerr << "matrix index is invalid\n";
  return array[i + size[0]*j]; 
}
double& matrix::operator() (int i, int j) { //set-value index operator
  if(i<0 || i>=size[0] || j<0 || j>=size[1])
    cerr << "matrix index is invalid\n";
  return array[i + size[0]*j]; 
}
matrix& matrix::operator= (const matrix& A) { //equate-two-mats operator
  if(size[0] != A.size[0] || size[1] != A.size[1])
    cerr << "matrix equality undefined -- sizes differ\n";
  for(int i=0; i<size[0]*size[1]; i++) array[i] = A.array[i];
  return *this;
}
matrix& matrix::operator= (double x) { //fill-mat-by-scalar operator
  for(int i=0; i<size[0]*size[1]; i++) array[i] = x;
  return *this;
}
matrix matrix::operator+ (const matrix& A) { //add-two-mats operator
  if(size[0] != A.size[0] || size[1] != A.size[1])
    cerr << "matrix addition undefined -- sizes differ\n";
  matrix sum(size[0],size[1]);
  for(int i=0; i<size[0]*size[1]; i++) 
    sum.array[i] = array[i] + A.array[i];
  return sum;
}
matrix matrix::operator- (const matrix& A) { //subtract-two-mats operator
  if(size[0] != A.size[0] || size[1] != A.size[1])
    cerr << "matrix subtraction undefined -- sizes differ\n";
  matrix difference(size[0],size[1]);
  for(int i=0; i<size[0]*size[1]; i++) 
    difference.array[i] = array[i] - A.array[i];
  return difference;
}
int dim(const matrix& A, int k) { //dimension function
  if(k != 0 && k!= 1)
    cerr << "matrix dimension index must be 0 or 1\n";
  return A.size[k];
}
matrix scaleMat(const double x, const matrix& A) { //scale mat function
  matrix product(dim(A,0),dim(A,1));
  for(int i=0; i<dim(A,0)*dim(A,1); i++) 
    product.array[i] = x*A.array[i];
  return product;
}
matrix matMatMult(const matrix& A, const matrix& B) { //mat*mat function
  if(dim(A,1) != dim(B,0))
    cerr << "matrix-matrix product undefined -- size mismatch\n";
  matrix product(dim(A,0),dim(B,1));
  for(int i=0; i<dim(A,0); i++) {
  for(int j=0; j<dim(B,1); j++) {
    double sum = 0;
    for(int k=0; k<dim(A,1); k++) {
      sum += A(i,k)*B(k,j);
    }
    product(i,j) = sum;
  }
  }
  return product;
}

/****************** 
define vector class 
*******************/
vector::vector(){ //default constructor 
  size = 0;
  array = 0;
}
vector::vector(int n){ //constructor 
  size = n;
  array = new double[size];
}
vector::~vector(){ //destructor 
  if(array != 0) delete []array;
}
vector::vector(const vector& v) { //copy constructor
  size = v.size;
  array = new double[size];
  for(int i=0; i<size; i++) array[i] = v.array[i];
}
double vector::operator() (int i) const { //get-value operator
  if(i<0 || i>=size)
    cerr << "vector index is invalid\n";
  return array[i]; 
}
double& vector::operator() (int i) { //set-value operator
  if(i<0 || i>=size)
    cerr << "vector index is invalid\n";
  return array[i]; 
}
vector& vector::operator= (const vector& v) { //equate-two-vecs operator
  if(size != v.size)
    cerr << "vector equality undefined -- sizes differ\n";
  for(int i=0; i<size; i++) array[i] = v.array[i];
  return *this;
}
vector& vector::operator= (double x) { //fill-vec-by-scalar operator
  for(int i=0; i<size; i++) array[i] = x;
  return *this;
}
vector vector::operator+ (const vector& v) { //add-two-vecs operator
  if(size != v.size)
    cerr << "vector addition undefined -- sizes differ\n";
  vector sum(size);
  for(int i=0; i<size; i++) sum.array[i] = array[i] + v.array[i];
  return sum;
}
vector vector::operator- (const vector& v) { //subtract-two-vecs operator
  if(size != v.size)
    cerr << "vector subtraction undefined -- sizes differ\n";
  vector difference(size);
  for(int i=0; i<size; i++) difference.array[i] = array[i] - v.array[i];
  return difference;
}
int dim(const vector& v) { //vec dimension function
  return v.size;
}
vector scaleVec(const double x, const vector& v) { //scale vec function
  vector product(dim(v));
  for(int i=0; i<dim(v); i++) product.array[i] = x*v.array[i];
  return product;
}
vector matVecMult(const matrix& A, const vector& x) { //mat*vec function
  if(dim(A,1) != dim(x))
    cerr << "matrix-vector product undefined -- size mismatch\n";
  vector product(dim(A,0));
  for(int i=0; i<dim(A,0); i++) {
    double sum = 0;
    for(int k=0; k<dim(A,1); k++) {
      sum += A(i,k)*x(k);
    }
    product(i) = sum;
  }
  return product;
}


/***************************************************** 
define input, ouput operators for vectors and matrices 
******************************************************/
ostream& operator<< (ostream& os, const vector& v) {
  for(int i=0; i<dim(v); i++) {
    os << v(i) << "\n";
  }
  return os;
}
ostream& operator<< (ostream& os, const matrix& A) {
  for(int i=0; i<dim(A,0); i++) {
    for(int j=0; j<dim(A,1); j++) {
      os << A(i,j) << "  ";
    }
    os << "\n";
  }
  return os;
}
istream& operator>> (istream& is, vector& v) {
  for(int i=0; i<dim(v); i++) is >> v(i);
  return is;
}
istream& operator>> (istream& is, matrix& A) {
  for(int i=0; i<dim(A,0); i++)
  for(int j=0; j<dim(A,1); j++) {
    is >> A(i,j);
  }
  return is;
}

/************************************************ 
define utility functions for vectors and matrices
*************************************************/
double vecL2Norm(const vector& v) {
  double norm = 0;
  for(int i=0; i<dim(v); i++) norm += v(i)*v(i);
  return sqrt(norm);
}
double vecMaxNorm(const vector& v) {
  double norm = 0;
  for(int i=0; i<dim(v); i++) {
    double a = fabs(v(i));
    if(norm < a) norm = a;
  }
  return norm;
}
double matMaxNorm(const matrix& A) {
  double norm = 0;
  for(int i=0; i<dim(A,0); i++) {
    double sum=0;
    for(int j=0; j<dim(A,1); j++) {
      sum += fabs(A(i,j));
    }
    if(norm < sum) norm = sum;
  }
  return norm;
}
double vecDot(const vector& v1, const vector& v2) { 
  if(dim(v1) != dim(v2))
    cerr << "vector dot product undefined -- sizes differ\n";
  double dotprod = 0;
  for(int i=0; i<dim(v1); i++) dotprod += v1(i)*v2(i);
  return dotprod;
}
