// C std lib
#include <cstdlib>

// for in/output
#include <iostream>
// for saving the solutions:
#include <vector>
// for MPI
#include <mpi.h>
// include your favorite sequential sorting implementation
//#include <sort.h>
// for timing purpose
#include "timer.hpp"

// our MPI Tags
//#define OUR_TAG_SIZE 13
//#define OUR_TAG_DATA 42

// debugging output:
//#define DEBUG(msg)
#define DEBUG(msg) std::cerr << msg << std::endl;

/* Setup timers */
timer counter;
double t_mpi_sort_begin, t_mpi_sort_end_first, t_mpi_sort_end_all;

/**
 * @brief   Performs the work.
 *
 * @param n     The size of the problem.
 * you can add more parameters if you need them.
 */
void work(unsigned int n)
{
    //Do the work that needs to happen and make the data sorted.
    //Receive some data from the processor 0  and communicate with my neighbors
    //to produce the sorted list.
}


/**
 * @brief   Prints a solution.
 *
 * @param solution  The solution to be printed.
 */
void print_solution(std::vector<unsigned int>& solution)
{
    // print out the sorted data
    // this should only run on the master because we do not have parallel IO
    // capabilities on the cluster.
}

/**
 * Prints the usage of the program.
 */
void print_usage()
{
    std::cerr << "Improve the usage text";
}


int main(int argc, char *argv[])
{
    /***************************
     *  Parse input arguments  *
     ***************************/

    // forget about first argument (which is the executable's name)
    argc--;
    argv++;

    // parse optional parameters if your program needs them
    while(argc > 0 && argv[0][0] == '-')
    {
        char option = argv[0][1];
        switch (option)
        {
            default:
                print_usage();
                exit(EXIT_FAILURE);
        }
        // iterate to next argument
        argv++;
        argc--;
    }

    // check that the mandatory parameters are present
    if (argc < 2)
    {
        print_usage();
        exit(EXIT_FAILURE);
    }
    // parse mandatory parameters
    unsigned int n = atoi(argv[0]);
    char* filepath = argv[1];

    /********************
     *  Setting up MPI  *
     ********************/

    // set up MPI
    int numprocessors, rank;
    MPI_Init(&argc, &argv);
    // get total size of processors and current rank
    MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    /*******************************
     *  Call master or worker code  *
     *******************************/

    if (rank == 0)
    {
        /* I am the root process */

        // start timer
        t_mpi_sort_begin = counter.get_ms();
        //read filepath and distribute the data to other processors

        // call the function that does the actual work
        work(n);

        // get elapsed time
        // the printing of the solution is not timed
        t_mpi_sort_end_all = counter.get_ms();

        // print output
        //print the output however you think is best
        // print timing
    }
    else
    {
        work(n)
    }

    // finish up
    MPI_Finalize();
    return 0;
}
