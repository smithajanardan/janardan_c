#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<mpi.h>
#include<fstream>
#include<math.h>
#include <cstdlib>
#include <iomanip>
#define DEBUG(msg) std::cerr << msg << std::endl;

// Wrriten by Brandon Bakr  
// Editted by Smitha Janardan
// Date: 04/25/2014

using namespace std;

//function to sum arrays
int sumarray(int *buf,int size) {
	int i = 0;
	int total = 0;
	while (i<size) {
		total = total + buf[i];
		i++;
	}
	return total;
}

void slave(int numtasks, int rank)
{
	/**************************
     *  Set local parameters  *
     **************************/
	
	int dims[2], rank_coords[2], periodic[2], col_coords[2], row_coords[2], info[4], neigh[4];
	int color, my_grid_rank, size, local_m, local_n, generations;
	int row_index, col_index, up_rank, down_rank, right_rank, left_rank, d;
	int *count, *grid_decomp, *index, *part_grid, *send_buf, *info_count, *info_index;
	MPI_Comm slave_comm, grid_comm, row_comm, col_comm;
	MPI_Request *rr, *sr, rr_l, sr_l;
	MPI_Status s;
	MPI_Datatype col;
	
	color = 1;
	d = sqrt(numtasks-1);
	count = (int *)malloc(numtasks*sizeof(int));
	rr = (MPI_Request *)malloc(4 * sizeof(MPI_Request));
	sr = (MPI_Request *)malloc(4 * sizeof(MPI_Request));
	
	/******************************
     *  Create the grid topology  *
     ******************************/
	
	// take the master processor off the communicator
	MPI_Comm_split(MPI_COMM_WORLD, color, rank, &slave_comm);
	
	// Set up Cartesian Grid
	dims = {sqrt(numtasks-1), sqrt(numtasks-1)};
	periodic = {0,0};
	col_coords = {0,1};
	row_coords = {1,0};
	
	// Create the Cartesian Grid
	MPI_Cart_create(slave_comm, 2, dims, periodic, 0, &grid_comm);
	MPI_Comm_rank(grid_comm, &my_grid_rank);
	MPI_Cart_coords(grid_comm, my_grid_rank, 2, rank_coords);
	MPI_Cart_sub(grid_comm, col_coords, &col_comm);
	MPI_Cart_sub(grid_comm, row_coords, &row_comm);
	
	/******************************
     *  Retrieve Local Life Grid  *
     ******************************/

	//Get size of receive buffer and dimensions of each local grid
	MPI_Scatterv(send_buf, info_count, info_index, MPI_INT, info, 4, MPI_INT, 0, MPI_COMM_WORLD);

	//Size and local dimensions
	size = info[0];
	local_m = info[1];
	local_n = info[2];
	generations = info[3];
	
	//Get part of grid
	part_grid = (int *)malloc(size*sizeof(int));
	MPI_Scatterv(grid_decomp, count, index, MPI_INT, part_grid,size,MPI_INT, 0, MPI_COMM_WORLD);

	//Form Grid from received array
	int bigger_grid[local_m+2][local_n+2];
	int local_grid[local_m][local_n];

	int k = 0;
	for (int i = 0; i < local_m+2; i++){
		for (int j = 0; j < local_n+2; j++) {
			//if (rank == 1) cout << i << " " << j << endl;
			if ((i == 0) || (i == local_m+1)) 
			{
				//if (rank == 1) cout << " here " << endl;
				bigger_grid[i][j] = 0;
			}
			else if ((j == 0) || (j == local_n+1)) 
			{
				bigger_grid[i][j] = 0;
			}
			else 
			{
				local_grid[i-1][j-1] = part_grid[k];
				bigger_grid[i][j] = local_grid[i-1][j-1];
				k++;
			}
		}
	}
	
	/*************************************
     *  Set up persistent communication  *
     *************************************/
	
	// find the row and column index of the processor
	row_index = rank_coords[0]; 	
	col_index = rank_coords[1];

	// find the ranks of the neighbouring processors
	up_rank = row_index+1;
	down_rank = row_index-1;
	right_rank = col_index+1;
	left_rank = col_index-1;

	// create the column datatype
	MPI_Type_vector(local_n,1, local_n+2, MPI_INT, &col);
	MPI_Type_commit(&col);

	// create all persistent communications
	
	if (left_rank >= 0)
	{
		MPI_Recv_init(&bigger_grid[1][0], 1, col, left_rank, 0, col_comm, &rr[0]);
		MPI_Send_init(&bigger_grid[1][1], 1, col, left_rank, 1, col_comm, &sr[1]);
	}
	if (right_rank < d)
	{
		MPI_Recv_init(&bigger_grid[1][local_m+1], 1, col, right_rank, 1, col_comm, &rr[1]);
		MPI_Send_init(&bigger_grid[1][local_m], 1, col, right_rank, 0, col_comm, &sr[0]);
	}
	
	if (up_rank < d)
	{
		MPI_Recv_init(&bigger_grid[local_n+1], local_m+2, MPI_INT, up_rank, 2, row_comm, &rr[2]);
		MPI_Send_init(bigger_grid[local_n], local_m+2, MPI_INT, up_rank, 3, row_comm, &sr[3]);
	}
	if (down_rank >= 0)
	{
		MPI_Recv_init(&bigger_grid[0], local_m+2, MPI_INT, down_rank, 3, row_comm, &rr[3]);
		MPI_Send_init(bigger_grid[1], local_m+2, MPI_INT, down_rank, 2, row_comm, &sr[2]);
	}
	
	/*************************************
     *  Play the Game of Life  *
     *************************************/
	 
	while( generations > 0)
	{
		int neigh[local_m][local_n];
		
		if (left_rank >= 0)
		{
			MPI_Start( &rr[0] );
			MPI_Start( &sr[1] );
		}
		if (right_rank < d)
		{
			MPI_Start( &sr[0] );
			MPI_Start( &rr[1] );
		}
		
		for (int i = 2; i < (local_m/2); i++)
		{
			for (int j = 2; j < (local_n); j++)
			{
				neigh[i-1][j-i] = bigger_grid[i-1][j-1] + bigger_grid[i-1][j] + bigger_grid[i-1][j+1] + bigger_grid[i][j+1] + bigger_grid[i+1][j+1] + bigger_grid[i+1][j] + bigger_grid[i+1][j-1] + bigger_grid[i][j-1];
			}
		}
		
		if (right_rank < d)
		{
			MPI_Wait( &sr[0], &s );
			MPI_Wait( &rr[1], &s );
		}
		if (left_rank >= 0)
		{
			MPI_Wait( &sr[1], &s );
			MPI_Wait( &rr[0], &s );
		}

		if (up_rank < d)
		{
			MPI_Start( &rr[2] );
			MPI_Start( &sr[3] );
		}
		if (down_rank >= 0)
		{
			MPI_Start( &sr[2] );
			MPI_Start( &rr[3] );
		}
		
		for (int i = (local_m/2); i < (local_m); i++)
		{
			for (int j = 2; j < (local_n); j++)
			{
				neigh[i-1][j-i] = bigger_grid[i-1][j-1] + bigger_grid[i-1][j] + bigger_grid[i-1][j+1] + bigger_grid[i][j+1] + bigger_grid[i+1][j+1] + bigger_grid[i+1][j] + bigger_grid[i+1][j-1] + bigger_grid[i][j-1];
			}
		}

		if (up_rank < d)
		{
			MPI_Wait( &rr[2], &s );
			MPI_Wait( &sr[3], &s );
		}

		if (down_rank >= 0)
		{
			MPI_Wait( &sr[2], &s );
			MPI_Wait( &rr[3], &s );
		}
		
		for (int i = 1; i < (local_m+1); i=i+local_m-1)
		{
			for (int j = 1; j < (local_n+1); j++)
			{
				neigh[i-1][j-1] = bigger_grid[i-1][j-1] + bigger_grid[i-1][j] + bigger_grid[i-1][j+1] + bigger_grid[i][j+1] + bigger_grid[i+1][j+1] + bigger_grid[i+1][j] + bigger_grid[i+1][j-1] + bigger_grid[i][j-1];
			}
		}
		
		for (int i = 2; i < (local_m); i++)
		{
			for (int j = 1; j < (local_n+1); j=j+local_n-1)
			{
				neigh[i-1][j-1] = bigger_grid[i-1][j-1] + bigger_grid[i-1][j] + bigger_grid[i-1][j+1] + bigger_grid[i][j+1] + bigger_grid[i+1][j+1] + bigger_grid[i+1][j] + bigger_grid[i+1][j-1] + bigger_grid[i][j-1];
			}
		}
		
		
		for (int i = 0; i < local_m; i++)
		{
			for (int j = 0; j < (local_n); j++)
			{
				int curr_cell=2;
				if (bigger_grid[i+1][j+1]==1)
				{
					if ((neigh[i][j] < 2) || (neigh[i][j] > 3))
					{
						curr_cell=0;
					}
					else
					{
						curr_cell=1;
					}
				}
				else
				{
					if (neigh[i][j]==3)
					{
						curr_cell=1;
					}
					else
					{
						curr_cell=0;
					}
				}
				bigger_grid[i+1][j+1] = curr_cell;
			}
		}
		generations--;
	}
	
	/*********************************************
     *  Release All Communicators and Datatypes  *
     *********************************************/
	
	// release persistent communication
	if (right_rank < d)
	{
		MPI_Request_free(&sr[0]);
		MPI_Request_free(&rr[1]);
	}
	if (left_rank >= 0)
	{
		MPI_Request_free(&sr[1]);
		MPI_Request_free(&rr[0]);
	}
	if (up_rank < d)
	{
		MPI_Request_free(&rr[2]);
		MPI_Request_free(&sr[3]);
	}
	if (down_rank >= 0)
	{
		MPI_Request_free(&sr[2]);
		MPI_Request_free(&rr[3]);
	}
	
	// release datatype
	MPI_Type_free(&col);
	
	// release the grid
	MPI_Comm_free(&slave_comm);
	MPI_Comm_free(&grid_comm);
	MPI_Comm_free(&row_comm);
	MPI_Comm_free(&col_comm);
}

void master(string input, string output, int numtasks)
{
	int color = 0;
	
	MPI_Comm slave_comm;
	MPI_Comm_split(MPI_COMM_WORLD, color, 0, &slave_comm);
	
	
	int m,n,generations;
	
     /*************************
     *  Read the input file  *
     *************************/
	fstream file(input.c_str(), ios::in);
	
	// read the input size
	file >> m;
	file >> n;
	file >> generations;
	

	// Sets size of life grid
	int grid[m][n];

	// load input into array
	for(int i=0; i<m; i++)
	{
		for(int j = 0; j<n; j++) file >> grid[i][j];
		
	}
	file.close();

	// catch easy parameters
	if (generations == 0) {
		//print_2darray(grid, m, n);
		return 0;
	}

	int d = sqrt(numtasks);
	int row_base = floor(m/d);
	int col_base = floor(n/d);
	int row_extra = m % d;
	int col_extra = n % d;

	// assigns portions of grid to processors	
	int x[numtasks],y[numtasks];
	for (int i = 0; i<numtasks; i++) {
		x[i] = row_base;
		y[i] = col_base;
		if (floor(i/d) < row_extra) x[i]=x[i]++;
		if (i % d < col_extra) y[i]=y[i]++;
	}

	int grid_decomp[m*n];
	
	//Finds parts of grid each process needs
	int base_i=0;
	int base_j=0;
	int k = 0;
	
	for (int p = 0; p < numtasks; p++){
		if ((p%d == 0) && (p != 0)) base_i = base_i + x[p-1];
		for (int i = 0; i < x[p]; i++){
			for (int j = 0; j < y[p]; j++){
				grid_decomp[k] = grid[i+base_i][j+base_j];
				k++;		
			}
		}
		base_j = (base_j + y[p]) % n;
	}
	

	//Set distribution of nubers
	int *count,*info;
	count = (int *)malloc((numtasks+1)*sizeof(int));
	info = (int *)malloc(4*(numtasks+1)*sizeof(int));
	count[0] = 0;
	info[0] = 0;
	info[1] = 0;
	info[2] = 0;
	info[3] = 0;
	for (int i = 1; i < numtasks+1; i++) {
		count[i] = x[i-1]*y[i-1];
		info[4*i] = x[i-1]*y[i-1];
		info[4*i+1] = x[i-1];
		info[4*i+2] = y[i-1];
		info[4*i+3] = generations;
	}

	int info_count[numtasks+1],info_disp[numtasks+1],*info_recvbuf;
	info_count[0] = 0;
	info_disp[0] = 0;
	for (int i = 1; i < numtasks+1; i++) info_count[i] = 4;
	for (int i = 1; i < numtasks+1; i++) info_disp[i] = 4*i;
	MPI_Scatterv(info,info_count, info_disp, MPI_INT, info_recvbuf, 4*(numtasks), MPI_INT, 0, MPI_COMM_WORLD);
	int *disp,*recvbuf;
	disp = (int *)malloc((numtasks+1)*sizeof(int));
	disp[0]=0;
	for (int i = 1; i < numtasks+1; i++) disp[i] = disp[i-1] + count[i-1];
	int size = sumarray(count,numtasks+1);
	MPI_Scatterv(grid_decomp,count, disp, MPI_INT, recvbuf, m*n, MPI_INT, 0, MPI_COMM_WORLD); 
}

/**
 * Prints the usage of the program.
 */
void print_usage()
{
    cerr << "The format of the input should be './play <input_file.txt> <output_file.txt>'\n";
	cerr << "mpirun -np p ./play is the name of the program.\n";
	cerr << "p is the number of processors. This must be at least 2 and be a perfect square + 1.\n";
	cerr << "input_file.txt is file name of the input file.\n";
	cerr << "output_file.txt is the file name of the output file.\n";
}

int main(int argc, char** argv)
{
    /***************************
     *  Parse input arguments  *
     ***************************/
	 
	// check that the mandatory parameters are present
	if (argc < 3)
	{
		print_usage();
		exit(EXIT_FAILURE);
 	}
	
	// get input/output file names
	string input = string(argv[1]);
	string output = string(argv[2]);

    /********************
     *  Setting up MPI  *
     ********************/
	 
	// initialize MPI and start all processors
	MPI::Init(argc, argv);
	MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
	
	// get local rank, and total number of processors
	int rank = MPI::COMM_WORLD.Get_rank();
	int numtasks = MPI::COMM_WORLD.Get_size();

	int av_tasks = numtasks - 1;	
	// fail if number of processors is perfect square
	if (sqrt(av_tasks) != floor(sqrt(av_tasks)))
	{
		print_usage();
		exit(EXIT_FAILURE);
		return MPI_Finalize();
	}

	// fail if no slave processors exist
	if (numtasks == 1)
	{
		print_usage();
		exit(EXIT_FAILURE);
		return MPI_Finalize();
	}

	if (rank == 0)
	{
		master(input,output,av_tasks);
	}
	else
	{
		slave(numtasks,rank);
	}
	
	// finish up
	return MPI_Finalize();
}
