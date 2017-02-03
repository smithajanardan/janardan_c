#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include "mpi.h"

using namespace std;

const char myfile[50]="/nethome/sjanardan3/project1_parallel_1.out" ;
ofstream prt(myfile) ;

int SIZE = 8;
int DEPTH = 3 ;
int MASTER = 0 ;
vector<int> solutions ;

void compress(vector<int> positions) 
{
	vector<int>:: iterator i;
	for(i=positions.begin(); i != positions.end(); i++) 
	{
		solutions.push_back(*i) ; 
	} 
}

void print() 
{
	if (solutions.empty() == true) 
	{ 
	}
	else 
	{
		vector<int>::iterator i;
		for(i = solutions.begin(); i != solutions.end(); i++) 
		{
			int n = 0;
			prt << *i << " " ;
			n++ ;
			if ( n % SIZE == 0) 
			{
				prt << endl ; 
			}
		} 
	} 
	solutions.clear(); 
}


bool check(int column, int row, vector<int> positions) 
{
	int n = 0;
	vector<int>::iterator set;
	for(set=positions.begin(); set != positions.end(); set++) 
	{
		if (*set == row) 
		{
			return false; 
		}
		else if (abs(n - column) == abs(*set - row)) 
		{
			return false; 
		} 
		n++; 
	}
	return true; 
}

void recursion(int column, vector<int> positions, int rank, int numtasks) 
{
	if((column == DEPTH-1) && (rank == MASTER)) 
	{
		MPI_Status status;
		int row = 0;
		int numsolutions;
		for (rank = 1; rank < numtasks; ++rank) 
		{
			MPI_Send(&row, 1, MPI_INTEGER, rank, 0, MPI::COMM_WORLD) ;
			MPI_Send(&column, 1, MPI_INTEGER, rank, 0, MPI::COMM_WORLD) ;
			MPI_Send(&positions, DEPTH, MPI_INTEGER, rank, 0, MPI::COMM_WORLD) ; 
			row++ ;
		}
		while (row < SIZE) 
		{
			MPI_Recv(&numsolutions, 1, MPI_INTEGER, MPI_ANY_SOURCE, 0, MPI::COMM_WORLD, &status) ;
			MPI_Recv(&solutions, numsolutions, MPI_INTEGER, MPI_ANY_SOURCE, 0, MPI::COMM_WORLD, &status) ;
			print() ;
			MPI_Send(&row, 1, MPI_INTEGER, status.MPI_SOURCE, 0, MPI::COMM_WORLD) ;
			MPI_Send(&column, 1, MPI_INTEGER, status.MPI_SOURCE, 0, MPI::COMM_WORLD) ;
			MPI_Send(&positions, DEPTH, MPI_INTEGER, status.MPI_SOURCE, 0, MPI::COMM_WORLD) ; 
			row++ ; 
		} 
		for (rank = 1; rank < numtasks; ++rank) 
		{
			MPI_Recv(&numsolutions, 1, MPI_INTEGER, MPI_ANY_SOURCE, 0, MPI::COMM_WORLD, &status) ;
			MPI_Recv(&solutions, numsolutions, MPI_INTEGER, MPI_ANY_SOURCE, 0, MPI::COMM_WORLD, &status) ;	
		}
	}
	else if (column == (SIZE-1)) 
	{
		for(int row = 0; row < SIZE; row++) 
		{
			if (check(column, row, positions) == true) 
			{
				positions.push_back(row);
				compress(positions);
				positions.pop_back(); 
			} 
		} 
	}
	else 
	{
		for(int row = 0; row < SIZE; row++) 
		{
			if (check(column, row, positions) == true) 
			{
				positions.push_back(row);
				recursion(column + 1, positions, rank, numtasks);
				positions.pop_back(); 
			} 
		} 
	} 
}
				
void kill_processors(int numtasks) 
{
	for (int n = 1; n < numtasks; n++) 
	{
		MPI_Send(0, 0, MPI::INT, n, 0, MPI::COMM_WORLD); 
	} 
}

int main()
{
	int rank, numtasks, rc, len, column, row;
	MPI::Init();
	rank = MPI::COMM_WORLD.Get_rank();
	numtasks = MPI::COMM_WORLD.Get_size();
	char name[MPI::MAX_PROCESSOR_NAME];
	memset(name,0,MPI::MAX_PROCESSOR_NAME);
	MPI::Get_processor_name(name,len);
	memset(name+len,0,MPI::MAX_PROCESSOR_NAME-len);
	vector<int> positions;
	prt.setf(ios::fixed) ;
	prt << setprecision(5) ;
	if(rank == MASTER) 
	{
		recursion(0, positions, rank, numtasks); 
		kill_processors(numtasks); 
	}
	if(rank != MASTER) 
	{
		while (1) 
		{
			MPI_Status stat;
			MPI_Recv(&row, 1, MPI_INT, 0, 0, MPI::COMM_WORLD, &stat);
			if (stat.MPI_TAG == 0) 
			{
				return 0; 
			}
			MPI_Recv(&column, 1, MPI_INT, 0, 0, MPI::COMM_WORLD, &stat);
			MPI_Recv(&positions, DEPTH, MPI_INTEGER, 0, 0, MPI::COMM_WORLD, &stat) ;
			if (check(column, row, positions) == true) 
			{
				positions.push_back(rank-1) ;
				recursion(column+1, positions, rank, numtasks) ; 
			}
			int n = solutions.size();
			MPI_Send(&n, 1, MPI_INT, 0, 0, MPI::COMM_WORLD);
			MPI_Send(&solutions, solutions.size(), MPI_INT, 0, 0, MPI::COMM_WORLD);
			solutions.clear() ;	
		} 
	}
	MPI::Finalize();
	return 0;
}
