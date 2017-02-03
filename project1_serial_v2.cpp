#include <iostream>
#include <iomanip>
#include <fstream>
#include <list>
#include <stdlib.h>

// Wrriten by Brandon Bakr  
// Editted by Smitha Janardan
// Date: 02/03/2014

using namespace std;

const char myfile[50]="/home/brandon/HPC/project1_serial.out" ;
ofstream prt(myfile) ;

void print(list<int> positions) {
  list<int>::iterator i;
  for(i=positions.begin(); i != positions.end(); ++i) prt << *i << " ";
  prt << endl ; }

bool check(int column, int row, list<int> positions) {
	int n = 0;
	list<int>::iterator set;
	for(set=positions.begin(); set != positions.end(); set++) {
		if (*set == row) {
			return false; }
		else if (abs(n - column) == abs(*set - row)) {
			return false; } 
		n++; }
	return true; }

void recursion(int size, int column, list<int> positions) {
	if(column == (size-1)) {
		for(int row = 0; row < size; row++) {
			if (check(column, row, positions) == true) {
				positions.push_back(row);
				print(positions);
				exit(0) ; } } }
	else {
		for(int row = 0; row < size; row++) {
			if (check(column, row, positions) == true) {
				positions.push_back(row);
				recursion(size, column + 1, positions);
				positions.pop_back(); } } } }

int main(){
  int n = 8;
  list<int> positions;
  prt.setf(ios::fixed) ;
  prt << setprecision(5) ;
  recursion(n, 0, positions);
  return 0 ; }
