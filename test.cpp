#include "timer.hpp"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "mpi.h"
#include <math.h> 
#define DEBUG(msg) std::cerr << msg << std::endl;

using namespace std;

int count_neigh(int **arr, int i, int j, int local_m, int local_n)
{
	if (i == 0)
	{
		if (j == 0)
		{
			return arr[i][j+1] + arr[i+1][j] + arr[i+1][j+1];
		}
		else if (j == local_m-1)
		{
			return arr[i][j-1] + arr[i+1][j] + arr[i+1][j-1];
		}
		else
		{
			return arr[i][j+1] + arr[i+1][j+1] + arr[i+1][j] + arr[i+1][j-1] + arr[i][j-1];
		}
	}
	else if (i == local_n-1)
	{
		if (j == 0)
		{
			return arr[i][j+1] + arr[i-1][j+1] + arr[i-1][j];
		}
		else if (j == local_m-1)
		{
			return arr[i][j-1] + arr[i-1][j-1] + arr[i-1][j];
		}
		else
		{
			return arr[i][j+1] + arr[i-1][j+1] + arr[i-1][j] + arr[i-1][j-1] + arr[i][j-1];
		}
	}
	else
	{
		if (j == 0)
		{
			return arr[i-1][j] + arr[i-1][j+1] + arr[i][j+1] + arr[i+1][j+1] + arr[i+1][j];
		}
		else if (j == local_m-1)
		{
			return arr[i-1][j] + arr[i-1][j-1] + arr[i][j-1] + arr[i+1][j-1] + arr[i+1][j];
		}
		else
		{
			return arr[i-1][j-1] + arr[i-1][j] + arr[i-1][j+1] + arr[i][j+1] + arr[i+1][j+1] + arr[i+1][j] + arr[i+1][j-1] + arr[i][j-1];
		}
	}
}

int main(int argc, char** argv)
{

    /********************
     *  Setting up MPI  *
     ********************/
	 
	// initialize MPI and start all processors
	MPI::Init(argc, argv);
	MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
	
	// get local rank, and total number of processors
	int rank = MPI::COMM_WORLD.Get_rank();
	int numtasks = MPI::COMM_WORLD.Get_size();
	
	/*******************************
     *  Call master or worker code  *
     *******************************/
	
	MPI_Comm col_comm;
	MPI_Comm row_comm;
	MPI_Request rr_l, rr_r, rr_u, rr_d, sr_l, sr_r, sr_u, sr_d;
	MPI_Status s;
	int r, c, row_index, inner_col_index;
	int d = 2;
	int local_m = 3;
	int local_n = 3;
	int bigger_grid[5][5];
	
	int grid[2][2] = {0};
	
	if (rank == 0)
	{
		bigger_grid = {{0,0,0,0,0},{0,1,1,1,0},{0,0,0,0,0},{0,0,1,0,0},{0,0,0,0,0}};
		r = 0;
		c = 0;
	}
	else if (rank == 1)
	{
		bigger_grid = {{0,0,0,0,0},{0,0,1,0,0},{0,1,1,0,0},{0,1,0,1,0},{0,0,0,0,0}};
		r = 0;
		c = 1;
	}
	else if (rank == 2)
	{
		bigger_grid = {{0,0,0,0,0},{0,1,1,1,0},{0,1,1,0,0},{0,1,1,1,0},{0,0,0,0,0}};
		r = 1;
		c = 0;
	}
	else
	{
		bigger_grid = {{0,0,0,0,0},{0,0,1,1,0},{0,1,1,1,0},{0,1,0,0,0},{0,0,0,0,0}};
		r = 1;
		c = 1;
	}
	
	MPI_Comm_split(MPI::COMM_WORLD, r, rank, &row_comm);
	MPI_Comm_split(MPI::COMM_WORLD, c, rank, &col_comm);
	
	MPI_Comm_rank(row_comm, &row_index);
	MPI_Comm_rank(col_comm, &inner_col_index);
	
	int up_rank = inner_col_index+1;
	int down_rank = inner_col_index-1;
	
	int right_rank = row_index+1;
	int left_rank = row_index-1;
	
	MPI_Datatype col;
	MPI_Type_vector(local_n,1, local_n+2, MPI_INT, &col);
	MPI_Type_commit(&col);

	if (left_rank >= 0)
	{
		MPI_Recv_init(&bigger_grid[1][0], 1, col, left_rank, 0, row_comm, &rr_l);
		MPI_Send_init(&bigger_grid[1][1], 1, col, left_rank, 1, row_comm, &sr_r);
		MPI_Start( &rr_l );
		MPI_Start( &sr_r );
	}
	if (right_rank < d)
	{
		MPI_Recv_init(&bigger_grid[1][local_m+1], 1, col, right_rank, 1, row_comm, &rr_r);
		MPI_Send_init(&bigger_grid[1][local_m], 1, col, right_rank, 0, row_comm, &sr_l);
		MPI_Start( &sr_l );
		MPI_Start( &rr_r );
	}
	if (up_rank < d)
	{
		MPI_Recv_init(&bigger_grid[local_n+1], local_m+2, MPI_INT, up_rank, 2, col_comm, &rr_u);
		MPI_Send_init(bigger_grid[local_n], local_m+2, MPI_INT, up_rank, 3, col_comm, &sr_d);
	}
	if (down_rank >= 0)
	{
		MPI_Recv_init(&bigger_grid[0], local_m+2, MPI_INT, down_rank, 3, col_comm, &rr_d);
		MPI_Send_init(bigger_grid[1], local_m+2, MPI_INT, down_rank, 2, col_comm, &sr_u);
	}

	// do work
	
	if (right_rank < d)
	{
		MPI_Wait( &sr_l, &s );
		MPI_Wait( &rr_r, &s );
		MPI_Request_free(&sr_l);
		MPI_Request_free(&rr_r);
	}
	if (left_rank >= 0)
	{
		MPI_Wait( &sr_r, &s );
		MPI_Wait( &rr_l, &s );
		MPI_Request_free(&sr_r);
		MPI_Request_free(&rr_l);
	}

	if (up_rank < d)
	{
		MPI_Start( &rr_u );
		MPI_Start( &sr_d );
	}
	if (down_rank >= 0)
	{
		MPI_Start( &sr_u );
		MPI_Start( &rr_d );
	}

	// do work

	if (up_rank < d)
	{
		MPI_Wait( &rr_u, &s );
		MPI_Wait( &sr_d, &s );
		MPI_Request_free(&rr_u);
		MPI_Request_free(&sr_d);
	}
	if (down_rank >= 0)
	{
		MPI_Wait( &sr_u, &s );
		MPI_Wait( &rr_d, &s );
		MPI_Request_free(&sr_u);
		MPI_Request_free(&rr_d);
	}
	
	MPI_Type_free (&col);
	// do work
	
	// finish up
	MPI_Finalize();
	return 0;
}