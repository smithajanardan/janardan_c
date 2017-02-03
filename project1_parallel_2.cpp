#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include "mpi.h"
#include <time.h>
#include <algorithm>

using namespace std;

// define output file to store solutions
const char myfile[50]="/nethome/sjanardan3/project1_parallel_2.out" ;
ofstream prt(myfile) ;

// user defined problem specifications
const int SIZE = 8;
const int DEPTH = 3 ;
const int solutions_size = 99;

// global values for all functions
int MASTER = 0 ;
int solutions[solutions_size][SIZE] ;
int numsolutions ;

// packing of information to be send to master processor
void compress(int positions[SIZE]) 
{
	for (int i = 0; i < SIZE; i++)
	{
		solutions[numsolutions][i] = positions[i] ;
	}
	numsolutions++;
}

// function that will take solutions and print them in the output file
void print(int numsolutions) 
{
	if (solutions == NULL ) 
	{ 
	}
	else 
	{
		for (int i = 0; i < numsolutions; i++)
		{
			for (int j = 0; j < SIZE; j++)
			{
				prt << solutions[i][j] << " " ;
			}
			prt << endl ;
		}
	}
}

// function to determine if given row/column is in strike-positions of another queen 
bool check(int column, int row, int positions[SIZE]) 
{
	for (int set=0; set < column; set++)
	{
		if (positions[set] == row) 
		{
			return false; 
		}
		else if (abs(set - column) == abs(positions[set] - row)) 
		{
			return false; 
		} 
	}
	return true; 
}

// recursively finds all solutions to problem
void recursion(int column, int positions[SIZE], int rank, int numtasks) 
{

	// if the master reaches the specified depth, it sends the rest to slaves
	if((column == DEPTH) && (rank == MASTER)) 
	{
		MPI_Status status;
		int row = 0;
		int numsolutions;
		
		// master sends one task to each slave
		for (int r = 1; r < numtasks; ++r) 
		{
			if (row < SIZE)
			{
				cout << "Master sending " << r << " positions: " ;
				
				// master compresses information before sending
				positions[DEPTH] = row;
				for (int i = 0; i < SIZE; i++)
				{
					cout << positions[i] << " ";
				}
				cout << endl;
				MPI_Send(positions, SIZE, MPI_INT, r, 1, MPI::COMM_WORLD) ; 
				row++ ;
			}
			else
			{
				break;
			}
		}
		
		// for extra tasks, it waits until a processor is free and then sends more tasks
		while (row < SIZE) 
		{
			cout << "Master receiving 2 " ;
			MPI_Recv(&numsolutions, 1, MPI_INTEGER, MPI_ANY_SOURCE, 1, MPI::COMM_WORLD, &status) ;
			MPI_Recv(&solutions, numsolutions, MPI_INTEGER, status.MPI_SOURCE, 1, MPI::COMM_WORLD, &status) ;
			cout << "from " << status.MPI_SOURCE << " solutions " << endl;
			cout << row << " " << column ;
			print(numsolutions) ;
			positions[DEPTH] = row;
			cout << " Master sending " << rank << " send_info " << endl;
			MPI_Send(positions, SIZE, MPI_INT, status.MPI_SOURCE, 1, MPI::COMM_WORLD) ; 
			row++ ; 
		} 
		
		// after all tasks are send, master waits for all remaining solutions
		for (int r = 1; r < min(numtasks, SIZE); ++r) 
		{
			cout << "Master receiving " << endl;
			MPI_Recv(&numsolutions, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI::COMM_WORLD, &status) ;
			cout << "Master received size from " <<status.MPI_SOURCE << " sent size" << endl;
			MPI_Recv(&solutions, solutions_size*SIZE, MPI_INT, status.MPI_SOURCE, 1, MPI::COMM_WORLD, &status) ;	
			cout << "Master received solutions from " << status.MPI_SOURCE << " solutions " << endl;
			print(numsolutions) ;
		}
	}
	
	// if the processor reaches the last column, it ends the solution and looks for others
	else if (column == (SIZE-1)) 
	{
		for(int row = 0; row < SIZE; row++) 
		{
			if (check(column, row, positions) == true) 
			{
				positions[column] = row;
				compress(positions);
				positions[column] = 0; 
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
				recursion(column + 1, positions, rank, numtasks);
				positions[column] = 0; 
			} 
		} 
	} 
}

// after everything is done, kills the processors				
void kill_processors(int numtasks) 
{
	for (int n = 1; n < numtasks; n++) 
	{
		MPI_Send(0, 0, MPI::INT, n, 0, MPI::COMM_WORLD); 
	} 
}

int main()
{
	clock_t t = clock();
	int positions[SIZE];
	int rank, numtasks, rc, len, column, row;

	// initialize MPI and start all processors
	MPI::Init();
	rank = MPI::COMM_WORLD.Get_rank();
	numtasks = MPI::COMM_WORLD.Get_size();
	char name[MPI::MAX_PROCESSOR_NAME];
	memset(name,0,MPI::MAX_PROCESSOR_NAME);
	MPI::Get_processor_name(name,len);
	memset(name+len,0,MPI::MAX_PROCESSOR_NAME-len);
	
	// open and prepare output file
	prt.setf(ios::fixed) ;
	prt << setprecision(5) ;
	
	// specify instructions for master processor
	if(rank == MASTER) 
	{
		recursion(0, positions, rank, numtasks); 
		
		// after everything is done, shut down all other processors and MPI
		kill_processors(numtasks); 
		cout << t << endl;
		MPI::Finalize();
		return 0;
	}
	
	// specify instructions for slave processor
	if(rank != MASTER) 
	{
		while (1) 
		{
		
			// while possible, receive data and send to recursion
			MPI_Status stat;
			cout << "Slave " << rank << " receiving  " << endl;
			MPI_Recv(&positions, SIZE, MPI_INT, 0, MPI::ANY_TAG, MPI::COMM_WORLD, &stat);
			cout << rank << " sending worked " << endl;
			for (int i = 0; i < SIZE; i++)
			{
				cout << positions[i] << " " ;
			}
			cout << endl;
			// kill the processor if the program is done
			if (stat.MPI_TAG == 0) 
			{
				cout << "tag is 0 " << endl;
				MPI::Finalize();
				return 0; 
			}
			
			// extract the information from the packed data
			cout << rank << " almost there " << endl;
			column = DEPTH;
			cout << rank << " made it" << endl;
			row = positions[DEPTH];
			cout << "Slave rank " << rank << " " << row << " " << column << endl;
			
			// start finding all possible solutions
			if (check(column, row, positions) == true) 
			{
				recursion(column+1, positions, rank, numtasks) ; 
			}
			
			// send solutions back to master processor
			cout << "Slave " << rank << " sending solutions" << endl;
			MPI_Send(&numsolutions, 1, MPI_INT, 0, 1, MPI::COMM_WORLD);
			cout << "Slave " << rank << " sent size to master" << endl ;
			MPI_Send(&solutions, solutions_size*SIZE, MPI_INT, 0, 1, MPI::COMM_WORLD);
			cout << "Slave " << rank << " sent solutions to master" << endl ;
			
		} 
	}
}
