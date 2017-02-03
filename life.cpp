#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<mpi.h>
#include<fstream>
#include<math.h>
#include <cstdlib>
#include <iomanip>
#define DEBUG(msg) std::cerr << msg << std::endl;

using namespace std;

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

//Function to print 2-d arrays
void print_2darray(int **ptr, int row, int col)
{
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) cout << ptr[i][j] << " " ;
		cout << endl;
	}
}

int master(string input, string output, int numtasks)
{
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
	return 0;
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
	
	int color = 0;
	if (rank != 0) color = 1;
	MPI_Comm slave_comm;
	MPI_Comm_split(MPI_COMM_WORLD, color, rank, &slave_comm);

	if (rank == 0) {
		int assignment = master(input,output,av_tasks);
	}

	if (rank != 0) {
		// Set up Cartesian Grid
		MPI_Comm grid_comm;
		int *dims,*periodic;
		dims = (int *)malloc(2*sizeof(int));
		periodic = (int *)malloc(2*sizeof(int));
		
		dims[0] = sqrt(numtasks-1);
		dims[1] = dims[0];
		periodic[0] = 0;
		periodic[1] = 0;
		int rank_coords[2];
		
		MPI_Cart_create(slave_comm, 2, dims, periodic, 0, &grid_comm);
		int my_grid_rank;
		
		MPI_Comm_rank(grid_comm, &my_grid_rank);

		MPI_Cart_coords(grid_comm, my_grid_rank, 2, rank_coords);
		
		int col_coords[2] = {0,1};
		int row_coords[2] = {1,0};
		
		MPI_Comm row_comm;
		MPI_Comm col_comm;

		MPI_Cart_sub(grid_comm, col_coords, &col_comm);
		MPI_Cart_sub(grid_comm, row_coords, &row_comm);
		
		int *count, info[4];
		int *grid_decomp, *index, *part_grid;
		int *send_buf, *info_count, *info_index;
		count = (int *)malloc(numtasks*sizeof(int));

		//Get size of receive buffer and dimensions of each local grid
		MPI_Scatterv(send_buf, info_count, info_index, MPI_INT, info, 4, MPI_INT, 0, MPI_COMM_WORLD);

		//Size and local dimensions
		int size = info[0];
		int local_m = info[1];
		int local_n = info[2];
		int generations = info[3];

		//Get part of grid
		part_grid = (int *)malloc(size*sizeof(int));
		MPI_Scatterv(grid_decomp, count, index, MPI_INT, part_grid,size,MPI_INT, 0, MPI_COMM_WORLD);

		//Form Grid from received array

		int *edges;
		int bigger_grid[local_m+2][local_n+2];
		int local_grid[local_m][local_n];

		int k = 0;
		
		for (int i = 0; i < local_m+2; i++){
			for (int j = 0; j < local_n+2; j++) {
				if ((i == 0) || (i == local_m+1)) 
				{
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
		
		
		int row_index = rank_coords[0]; 	
		int col_index = rank_coords[1];

		int up_rank = col_index+1;
		int down_rank = col_index-1;
		
		int right_rank = row_index+1;
		int left_rank = row_index-1;

		MPI_Request rr_l, rr_r, rr_u, rr_d, sr_l, sr_r, sr_u, sr_d;
		MPI_Status s;

		int d = sqrt(numtasks-1);

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
		
		// do work
		
	}
	// finish up
	return MPI_Finalize();
}
