#include "timer.hpp"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "mpi.h"
#include <math.h> 
#define DEBUG(msg) std::cerr << msg << std::endl;

// Wrriten by Smitha Janardan 
// Editted by Brandon Bakr
// Date: 04/01/2014

using namespace std;

timer counter;
double t_mpi_sort_begin, t_mpi_sort_end_first, t_mpi_sort_end_all;

/**
 * @brief   Pointer structure
 *
 * @param array		The pointer.
 */
struct pointer_array 
{
    int *array;
};

int random(int seed,int probsize) {          //Generates same random number on all processors
	srand(seed);
	return rand()%probsize;
}

int compare (const void * a, const void * b) //copmare function for qsort
{
  return ( *(int*)a - *(int*)b );
}

int sumarray(int *buf,int size) {
	int i = 0;
	int total = 0;
	while (i<size) {
		total = total + buf[i];
		i++;
	}
	return total;
}

/**
 * @brief   Creates an array of the receive counts.
 *
 * @param size			The size of the selected subgroup
 * @param out_count		The distribution on the numbers.
 * @param num_p			The number of processors for the subgroup
 * @param numtasks		The number of processors on the communicator.
 * @param r				The acting rank.
 */
pointer_array count(int size, int *old_count, int num_p, int numtasks, int r)
{

	/**************************
     *  Set local parameters  *
     **************************/
	
	int *count, extra, base, addition, nums_needed, numbers_below, total_sum;
	count = (int*)malloc(sizeof(int)*numtasks);
	total_sum = 0;
	
	/**********************************
     *  Find distribution of numbers  *
     **********************************/
	 
	// use code provided by TAs to distribute number evenly among processors
	extra = size%num_p;
	base = floor((float)size/(float)num_p);
	addition = extra;
	
	// find the indices of the numbers that the processor should receive
	// the receiving numbers span from number_below to (numbers_below + nums_needed)
	nums_needed = base;
	if (r < extra)
	{
		addition = r;
		nums_needed++;
	}
	numbers_below = (r*base) + addition;
	
	/**********************
     *  Make count array  *
     *********************/
	for (int i = 0; i < numtasks; i++)
	{
		// iterate through the distribution and find what processors the numbers are coming from
		total_sum = total_sum + old_count[i];
		if ((total_sum > numbers_below))
		{
			count[i] = min(old_count[i], nums_needed);
			nums_needed = nums_needed - count[i];
		}
		else{ count[i] = 0; }
	}
	
	// return the final receiving count distribution
	return {count};
}

/**
 * @brief   Sends the old numbers and gets the new numbers.
 *
 * @param numtasks		The total number of processors.
 * @param send			The old numbers to be sent.
 * @param receive_count	The distribution of where to receive the new numbers.
 * @param send_size		The size of the old numbers.
 * @param comm			The current communicator.
 * @param rank			The current rank.
 */
pointer_array send_receive(int numtasks, int *send, int *receive_count, int send_size, MPI_Comm comm, int rank)
{
	/**************************
     *  Set local parameters  *
     **************************/
	int *send_disp,*receive_disp,*send_count, *receive;
	send_count=(int*)malloc(sizeof(int)*numtasks);
	send_disp=(int*)malloc(sizeof(int)*numtasks);
	receive_disp=(int*)malloc(sizeof(int)*numtasks);
	
	/**********************************
     *  Find send count distribution  *
     **********************************/
	
	// send the receive count distribution to other processors
	// the processor gets the send count distribution back
	MPI_Alltoall(receive_count,1,MPI_INT,send_count,1,MPI_INT,comm);
	
	/***************************************
     *  Create displacement distributions  *
     ***************************************/
	send_disp[0]=0;	
	int receive_size=0;
	receive_disp[0]=0;
	
	// find the receive size, and both distributions
	for(int i=1;i<numtasks;i++)
	{
		receive_size=receive_size+receive_count[i];
		receive_disp[i]=receive_count[i-1]+receive_disp[i-1];
		send_disp[i]=send_count[i-1]+send_disp[i-1];
	}
	
	/**************************
     *  Get new local numbers  *
     **************************/
	 
	receive = (int*)malloc(sizeof(int)*receive_size);
	MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
	
	// use all to all v to send the current local numbers and recieve the new local numbers
	int error_code = MPI_Alltoallv(send,send_count,send_disp,MPI_INT,receive,receive_count,receive_disp,MPI_INT,comm);
	if (error_code != MPI_SUCCESS) {
		char error_string[BUFSIZ];
		int length_of_error_string;
		MPI_Error_string(error_code, error_string, &length_of_error_string);
		fprintf(stderr, "%3d: %s\n", rank, error_string); }
	
	// return the new local numbers
	return {receive};
}

/**
 * @brief   Recursively sorts the numbers.
 *
 * @param numtasks	The total number of processors.
 * @param numbers	The local integers.
 * @param comm		The current communicator.
 * @param rank		The current rank.
 * @param tot_size	The total size of the subgroup.
 */
void sort(int numtasks, int size, int *numbers, MPI_Comm comm, int rank, int tot_size, int round)
{
	/***********************************
     *  Base Case (Only one processor) *
     ***********************************/
	 
	if (numtasks == 1)
	{
		// sort the local numbers serially
		qsort(numbers, size, sizeof(int),compare);
	
		// send all of the numbers back to the master processor
		int *all_sizes, *master_list, *disp;
		MPI_Gather(&size, 1, MPI_INT, all_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Gatherv(numbers,size,MPI_INT,master_list,all_sizes,disp,MPI_INT,0,MPI_COMM_WORLD);
	}
	
	/********************************************
     *  Non-Base Case (more than one processor) *
     ********************************************/
	 
	else
	{
		/**************************
		 *  Set local parameters  *
		 **************************/
		
		int proc, extra, base, k, nk, bsize, asize, tot_below, tot_above, num_b, num_a, color;
		int *parted_nums, *nums_below, *nums_above, *new_nums, *c;
		MPI_Comm next = comm;
		extra = tot_size%numtasks;
		base = floor((float)tot_size/(float)numtasks);
		parted_nums = (int *)malloc(size*sizeof(int));
		nums_below = (int*)malloc(sizeof(int)*numtasks);
		nums_above = (int*)malloc(sizeof(int)*numtasks);
		bsize = 0;
		asize = 0;
		color = 0;
		
		/************************************
		*  Step 1: Find the random number k *
		*************************************/
		
		// find a random integer 
		k = random(tot_size, numtasks);
		if (rank == 0) { cout << "round " << round << " k " << k << endl; }
		// locate the processor of the kth integer
		// after this step, k will be the negative index of the numbers
		for (int i=0; i < numtasks; i++)
		{
			k = k - base;
			if (i < extra) { k--; }
			if (k <= 0)
			{
				proc = i;
				break;
			}
		}
		
		/***********************************
		*  Step 2: Broadcast the kth digit *
		************************************/
		if (proc == rank) { nk = numbers[size+k]; }
		MPI_Bcast(&nk,1,MPI_INT,proc,comm);
		if (rank == 0) { cout << "round " << round << " nk " << nk << endl; }
		
		/**********************************************************
		*  Step 3: Partition the local numbers into two subgroups *
		***********************************************************/

		// iterate through the entire array of local numbers
		for (int i = 0; i < size; i++)
		{
		
			// if the local number is below the threshold, add to beginning of list
			if (numbers[i] <= nk)
			{
				parted_nums[bsize] = numbers[i];
				bsize++;
			}
			
			// if local number is above the threshold, add to end of list
			else
			{
				parted_nums[(size-1)-asize] = numbers[i];
				asize++;
			}
		}
		if (rank == 0) { cout << "round " << round << " bsize " << bsize << endl; }
		if (rank == 0) { cout << "round " << round << " asize " << asize << endl; }
		
		/****************************************************************
		*  Step 4: Determine the number of processors for each subgroup *
		*****************************************************************/

		// use all gather to find number distributions
        MPI_Allgather(&asize, 1, MPI_INT, nums_below, 1, MPI_INT, comm);
        MPI_Allgather(&bsize, 1, MPI_INT, nums_above, 1, MPI_INT, comm);
		
		// find the total size of each subgroup
		tot_below = sumarray(nums_below, numtasks);
        tot_above = sumarray(nums_above, numtasks);
		
		if (rank == 0) { cout << "round " << round << " tot_below " << tot_below << endl; }
		if (rank == 0) { cout << "round " << round << " tot_above " << tot_above << endl; }
		
		// if the second subgroup is empty, add one number to the above group
		if (tot_above == 0)
		{
			nums_below[proc]--;
			nums_above[proc]++;
			tot_below--;
			tot_above++;
			if (rank == proc)
			{
				bsize--;
				asize++;
			}
		}
		
		/****************************************************
		*  Step 5: Find the total size of the two subgroups *
		*****************************************************/

		if (tot_below < tot_above) 
		{
			num_a = floor((float)(tot_above*numtasks)/(float)(tot_size));
			num_b = numtasks - num_a;
		}
		else
		{
			num_b = floor((float)(tot_below*numtasks)/(float)(tot_size));
			num_a = numtasks - num_b;
		}
		if (rank == 0) { cout << "round " << round << " num_b " << num_b << endl; }
		if (rank == 0) { cout << "round " << round << " num_a " << num_a << endl; }
		
		/************************************
		*  Step 6: Switch the local numbers *
		*************************************/
		
		// Find receive count distribution
		if (rank < num_b) { c = count(tot_below, nums_below, num_b, numtasks, rank).array; }
		else { c = count(tot_above, nums_above, num_a, numtasks, (rank-num_b)).array; }
		
		// retrieve the new local numbers
		new_nums = send_receive(numtasks,parted_nums,c,size,comm,rank).array;
		
		// Get new local size
		size = sumarray(c, numtasks);
		
		// Get the new total size
		if (rank < num_b) { tot_size = tot_below; }
		else { tot_size = tot_above; }

		/****************************************
		*  Step 7: Create two new communicators *
		*****************************************/
		
		// change the colour for the second subgroup
		if (rank >= num_b) { color = 1; }

		cout << "round " << round << " rank " << rank << " colour " << color << endl;
		// split the communicators based on colour
		MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
		int error_code = MPI_Comm_split(comm, color, rank, &next);
		if (error_code != MPI_SUCCESS) {

		   char error_string[BUFSIZ];
		   int length_of_error_string;

		   MPI_Error_string(error_code, error_string, &length_of_error_string);
		   fprintf(stderr, "%3d: %s\n", rank, error_string);
		}
		cout << "round " << round << " rank " << rank << " sort 2 " << endl;
		
		// get rid of the old communicator
		MPI_Comm_free(&comm);
		
		// the new rank and numtasks
		MPI_Comm_rank(next, &rank);
		MPI_Comm_size(next, &numtasks);
		cout << "round " << round << " rank " << rank << " sort 3 " << endl;
		
		/*******************************
		*  Step 8: Repeat recursively  *
		********************************/
		
		sort(numtasks,size,new_nums,next,rank,tot_size, round+1);
	}
}

/**
 * @brief   Preforms the work for slave processors
 *
 * @param numtasks  The total number of processors.
 * @param rank		The processor rank.
 */
void slave(int numtasks, int rank)
{
	/**************************
     *  Set local parameters  *
     **************************/
	 
	int *count, *list, *index, *numbers, size, n;
	count = (int *)malloc(numtasks*sizeof(int));

	/**************************
     *  Receive local numbers *
     **************************/
	 
	// find the distribution of numbers
	MPI_Bcast(count,numtasks, MPI_INT,0,MPI::COMM_WORLD);
	
	// abstract the total size of the array 
	size = count[0];
	count[0] = 0;
	
	// get the expected number of local numbers
	n = count[rank];
	numbers = (int *)malloc(n*sizeof(int));
	
	// get the local numbers
	MPI_Scatterv(list, count, index, MPI_INT, numbers, n, MPI_INT, 0, MPI::COMM_WORLD);

	// create new communicator without master processor
	int color = 1;
	MPI_Comm comm;
	MPI_Comm_split(MPI::COMM_WORLD, color, rank, &comm);
	
	// find the new rank and numtasks
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &numtasks);
	
	// run the function that sorts the numbers
	sort(numtasks, n, numbers, comm, rank, size, 0);
}

/**
 * @brief   Prints a solution.
 *
 * @param solution  The solution to be printed.
 * @param file		The specified output file.
 * @param listSize  The size of the input list.
 */
void print_solution(int *solution, string output, int listSize)
{
	fstream file(output.c_str(), ios::out);
	
	// dump output into file
	for(int i=0; i<listSize; i++)
	{
		file << solution[i] << "\n";
	}
	file.close();
}

/**
 * @brief   Preforms the work for master processor.
 *
 * @param input		The name of the input file.
 * @param output	The name of the output file.
 * @param numtasks  The total number of processors.
 */
void master(string input, string output, int numtasks)
{
	/**************************
     *  Set local parameters  *
     **************************/
	
	int extra_elems, base_size, color, size, listSize;
	int *list, *recvbuf, *disp, *count, *all_sizes, *master_list, *send;
	MPI_Comm comm;
	count = (int *)malloc(numtasks*sizeof(int));
	disp = (int *)malloc(numtasks*sizeof(int));
	all_sizes = (int*)malloc(sizeof(int)*numtasks);
	disp = (int *)malloc(numtasks*sizeof(int));
	size = 0;
	color = 0;
	
	/*************************
     *  Read the input file  *
     *************************/
	fstream file(input.c_str(), ios::in);
	
	// read the input size
	file >> listSize;
	list = new int[listSize];
	master_list = (int*)malloc(sizeof(int)*listSize);

	// load input into array
	for(int i=0; i<listSize; i++)
	{
		file >> list[i];
	}
	file.close();
	
	// start timer
	t_mpi_sort_begin = counter.get_ms();

	// extra elements left after n/p
	extra_elems = listSize%(numtasks-1);

	base_size = floor((float)listSize/(float)(numtasks-1));

	// start filling up index
	disp[0] = 0;
	count[0] = 0;
	for (int i=1; i < numtasks; i++) 
	{
		count[i] = base_size;
		disp[i] = count[i-1] + disp[i-1];
		if (i <= extra_elems) 
		{
			count[i]++;
		}
	}
	
	/**************************
     *  Send out the list  *
     **************************/
	
	// need to include the total size of the list in the sending array
	count[0] = listSize;

	// broadcast the distribution of numbers
	MPI_Bcast(count, numtasks, MPI_INT,0,MPI::COMM_WORLD);

	// change the count array back to the original
	count[0] = 0;

	// scatter the actual numbers
	MPI_Scatterv(list, count, disp, MPI_INT, recvbuf, listSize, MPI_INT, 0, MPI::COMM_WORLD);

	// create new communicator without the master processors
	MPI_Comm_split(MPI::COMM_WORLD, color, 0, &comm);
	
	/************************
     *  Get sorted numbers  *
     ************************/
	
	// receive number distribution from slaves
	MPI_Gather(&size, 1, MPI_INT, all_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	// calculate the displacement array
	disp[0] = 0;
	for(int i=1;i<numtasks;i++)
	{
		disp[i]=all_sizes[i-1]+disp[i-1];
	}

	// receive numbers
	MPI_Gatherv(send, 0, MPI_INT, master_list, all_sizes, disp, MPI_INT, 0, MPI_COMM_WORLD);
	
	// get elapsed time
	t_mpi_sort_end_all = counter.get_ms();
	
	//print time
	cout << numtasks << "	" << listSize << "	" << t_mpi_sort_end_all - t_mpi_sort_begin << endl;
	
	// print the sorted numbers
	print_solution(master_list, output, listSize);
}

/**
 * Prints the usage of the program.
 */
void print_usage()
{
    cerr << "The format of the input should be './program2 <input_file.txt> <output_file.txt>'\n";
	cerr << "program2 is the name of the program.\n";
	cerr << "input_file.txt is file name of the input file. This should include the size or the array, and all of the numbers.\n";
	cerr << "output_file.txt is the file name of the output file. This is the location of the final sorted array.\n";
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
	
	/*******************************
     *  Call master or worker code  *
     *******************************/
	
	if (rank == 0)
	{
	
		/* I am the root process */
		
		// call the function that does the actual work
		master(input, output, numtasks);
		
	}
	else
	{
	
		/* I am a slave processor */
		
		// call the function that does the actual work
		slave(numtasks, rank);
	}
	
	// finish up
	MPI_Finalize();
	return 0;
}
