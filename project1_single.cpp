#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include "mpi.h"
#include <time.h>
#include <algorithm>
#include <ctime>

// Wrriten by Brandon Bakr 
// Editted by Smitha Janardan 
// Date: 02/07/2014

using namespace std;

// define output file to store solutions
const char myfile[50]="1sol.out" ;
ofstream prt(myfile) ;

// user defined problem specifications
const int SIZE = 8;
const int DEPTH = 3 ;

// solutions size is one since we only want one solutions
const int solutions_size = 1;

// global values for all functions
int MASTER = 0 ;
int solutions[solutions_size*SIZE] ;
int numsolutions = 0;

// packing of information to be send to master processor
void compress(int positions[SIZE]) 
{
	for (int i = 0; i < SIZE; i++)
	{
		solutions[(numsolutions*SIZE)+i] = positions[i] ;
	}
	numsolutions++;
}

// after everything is done, kills the processors				
void kill_processors(int numprocessors) 
{
	int kill[SIZE];
	kill[0] = SIZE*SIZE ;
	MPI_Bcast(kill, SIZE, MPI_INT, 0, MPI::COMM_WORLD) ;
}

// print solutions if it comes from the master
void master_print(int positions[SIZE]) 
{
	for (int n = 0; n < SIZE; n++) 
	{
		prt << positions[n] << " "; 
	}
	prt << endl ; 
}

// print all the solutions to the output file
int print( int numsolutions, int *rbuf, int numprocessors)
{
	for (int i = 0; i < numsolutions*SIZE; i++)
	{
		prt << rbuf[i] << " " ;
		
		// after every complete solution, add an end line
		if (((i+1)%SIZE) == 0)
		{
			prt << endl ;
			return 1;
		}
	}
return 0;
}

// function to determine if given row/column is in strike-positions of another queen 
bool check(int column, int row, int positions[SIZE]) 
{
	for (int set=0; set < column; set++)
	{
	
		// check if the new position is in the same row
		if (positions[set] == row) 
		{
			return false; 
		}
		
		// check if the new position is on a diagonal
		else if (abs(set - column) == abs(positions[set] - row)) 
		{
			return false; 
		} 
	}
	return true; 
}

// recursively finds all solutions to problem
int recursion(int column, int positions[SIZE], int rank, int numprocessors) 
{

	// if the master reaches the specified depth, it sends the rest to slaves
	if((column == DEPTH) && (rank == MASTER) && (numprocessors > 1)) 
	{
	
		// send current positions to all the slaves
		int recieve_data[solutions_size][SIZE] ;
		int n;
		MPI_Bcast(positions, SIZE, MPI_INT, 0, MPI::COMM_WORLD) ;
		
		// after processors finish calculating the data, master collects number of solutions
		int *receive_num = (int *)malloc(numprocessors*sizeof(int));
		MPI_Gather(&numsolutions, 1, MPI_INT, receive_num, 1, MPI_INT, 0, MPI::COMM_WORLD) ; 
		
		//using the number of solutions, master prepares gathering buffers/arrays
		int *count = (int *)malloc(numprocessors*sizeof(int));
		count[0] = 0 ;
		for (int i=1; i<numprocessors; ++i)
		{
			count[i] = (count[i-1]+receive_num[i-1]); 
			numsolutions = numsolutions + (receive_num[i] / SIZE) ;
		}
		int *rbuf = (int *)malloc((count[numprocessors-1]+receive_num[numprocessors-1])*sizeof(int)); 
		int *send = (int *)malloc(numprocessors*(count[numprocessors-1]+receive_num[numprocessors-1])*sizeof(int)); 
		
		// master collects all solutions from slaves
		MPI_Gatherv(send, numsolutions*SIZE, MPI_INT, rbuf, receive_num, count, MPI_INT, 0, MPI::COMM_WORLD) ;
		
		// master prints them in the output file
		int state = print(numsolutions, rbuf, numprocessors);
		if (state == 1) return 1;
	}
	
	// if the processor reaches the last column, it ends the solution and looks for others
	else if (column == (SIZE-1)) 
	{
		for(int row = 0; row < SIZE; row++) 
		{
			if (check(column, row, positions) == true) 
			{
				positions[column] = row;
				
				// if it is the master processor, print right away
				if (rank == MASTER)
				{
					master_print(positions);
				}
				
				// if it is a slave processor, store the solutions until later
				else
				{
					compress(positions);
				}
				return 1;
			} 
		} 
	}
	
	// for all other rows, processors looks for possible empty slots and continues
	else 
	{
		for(int row = 0; row < SIZE; row++) 
		{
			if (check(column, row, positions) == true) 
			{
				positions[column] = row;
				int state = recursion(column + 1, positions, rank, numprocessors); 
				
				// if there is already a solution, return 1 until you exist the recursion
				if (state == 1) return 1;
			} 
		} 
	} 
}


// function for the master to start looking for solutions
int master_task(int rank, int numprocessors)
{
	// start looking for solutions
	int positions[SIZE];
	int state = recursion(0, positions, rank, numprocessors);
	
	// after everything is done, shut down all other processors and MPI
	kill_processors(numprocessors); 
	MPI::Finalize();
	return 0;
}

// function for the slave to receive positions and find solutions
int slave_task(int rank, int numprocessors)
{
	int positions[SIZE] ;
	int column, row;
	
	// keep receiving data until kill switch is given
	while(1) 
	{
		// receive positions
		MPI_Bcast(positions, SIZE, MPI_INT, 0, MPI::COMM_WORLD);
		
		// check if the received data is a kill switch
		if (positions[0] == SIZE*SIZE) 
		{
			MPI::Finalize();
			return 0; 
		}
		
		// divide up tasks based on processor number
		for (int row = 0; row < SIZE; row++)
		{
			if ((row%(numprocessors-1))+1 == rank)
			{
				
				// if the row is available, continue finding solutions 
				if (check(DEPTH, row, positions) == true) 
				{
					positions[DEPTH] = row;
					if (DEPTH != SIZE-1)
					{
						recursion(DEPTH+1, positions, rank, numprocessors) ; 
					}
					else 
					{
						compress(positions);
					}
				}
			}
		}
		
		// send number of solutions to the master
		int *receive_num = (int *)malloc(numprocessors*sizeof(int));
		int n = numsolutions*SIZE ;
		MPI_Gather(&n, 1, MPI_INT, receive_num, 1, MPI_INT, 0, MPI::COMM_WORLD) ;
		
		// send solutions to the master
		int *count = (int *)malloc(numprocessors*sizeof(int));
		int *rbuf = (int *)malloc(numprocessors*(count[numprocessors-1]+receive_num[numprocessors-1])*sizeof(int)); 
		MPI_Gatherv(solutions, numsolutions*SIZE, MPI_INT, rbuf, receive_num, count, MPI_INT, 0, MPI::COMM_WORLD) ;
		numsolutions = 0;
	} 
}

int main()
{
	// start clock
	clock_t t = clock();
	int rank, numprocessors;

	// initialize MPI and start all processors
	MPI::Init();
	rank = MPI::COMM_WORLD.Get_rank();
	numprocessors = MPI::COMM_WORLD.Get_size();
	
	// open and prepare output file
	prt.setf(ios::fixed) ;
	prt << setprecision(5) ;
	
	// specify instructions for master processor
	if(rank == MASTER) 
	{
		// the master tasks
		int state = master_task(rank, numprocessors) ;
		
		// print the time
		t = (clock() - t) / CLOCKS_PER_SEC ;
		cout << t << endl;
	}
	
	// specify instructions for slave processor
	else
	{
		slave_task(rank, numprocessors) ;
	}
	return 0;
}
