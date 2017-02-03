#include <iostream>
#include <iomanip>
#include <fstream>
#include <list>
#include <stdlib.h>
#include <time.h>

using namespace std;

// start the timer right away
double t = clock();

// define output file to store solutions
const char myfile[50]="/nethome/sjanardan3/project1_serial.out" ;
ofstream prt(myfile) ;

// user defined problem specifications
const int size = 13;

// print all the solutions to the output file
void print(int positions[size]) {
  for (int n = 0; n < size; n++) {
	prt << positions[n] << " "; }
  prt << endl ; }

// function to determine if given row/column is in strike-positions of another queen 
bool check(int column, int row, int positions[size]) {
	for (int set=0; set < column; set++){
	
		// check if the new position is in the same row
		if (positions[set] == row) 	{
			return false; 	}
		
		// check if the new position is on a diagonal
		else if (abs(set - column) == abs(positions[set] - row)) {
			return false; 	} }
	return true; }

// recursively finds all solutions to problem
void recursion(int column, int positions[size]) {

	// if the processor reaches the last column, it ends the solution and looks for others
	if(column == (size-1)) {
		for(int row = 0; row < size; row++) {
			if (check(column, row, positions) == true) {
				positions[column] = row;
				print(positions); 
				} } }
	
	// for all other rows, processors looks for possible empty slots and continues
	else {
		for(int row = 0; row < size; row++) {
			if (check(column, row, positions) == true) {
				positions[column] = row;
				recursion(column + 1, positions); } } } }

int main(){

  // // open and prepare output file
  int positions[size];
  prt.setf(ios::fixed) ;
  prt << setprecision(5) ;
  
  // look for solutions
  recursion( 0, positions);
  
  // print time
  t = (clock() - t) / CLOCKS_PER_SEC ;
  cout << t << endl;
  return 0 ; }
