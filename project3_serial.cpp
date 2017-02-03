#include "timer.hpp"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h> 
#define DEBUG(msg) std::cerr << msg << std::endl;

// Wrriten by Smitha Janardan 
// Editted by Brandon Bakr
// Date: 04/20/2014

using namespace std;

timer counter;
double t_mpi_sort_begin, t_mpi_sort_end_first, t_mpi_sort_end_all;

struct data 
{
    int **data;
};


int count_neigh(int **arr, int i, int j, int xlen, int ylen)
{
	//cout << "a" << endl;
	int neigh = 0;
	if (i == 0)
	{
		//cout << "b" << endl;
		if (j == 0)
		{
			//cout << "ba" << endl;
			return arr[i][j+1] + arr[i+1][j] + arr[i+1][j+1];
		}
		else if (j == xlen-1)
		{
			//cout << xlen << " " << ylen << endl;
			//cout << i-1 << " " << j-1 << endl;
			//cout << "1 " << arr[i][j-1] << endl;
			//cout << "2 " << arr[i+1][j] << endl;
			//cout << "3 " << arr[i+1][j-1] << endl;
			return arr[i][j-1] + arr[i+1][j] + arr[i+1][j-1];
		}
		else
		{
			//cout << "bc" << endl;
			return arr[i][j+1] + arr[i+1][j+1] + arr[i+1][j] + arr[i+1][j-1] + arr[i][j-1];
		}
	}
	else if (i == ylen-1)
	{
		//cout << "c" << endl;
		if (j == 0)
		{
			return arr[i][j+1] + arr[i-1][j+1] + arr[i-1][j];
		}
		else if (j == xlen-1)
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
		//cout << "d" << endl;
		if (j == 0)
		{
			return arr[i-1][j] + arr[i-1][j+1] + arr[i][j+1] + arr[i+1][j+1] + arr[i+1][j];
		}
		else if (j == xlen-1)
		{
			return arr[i-1][j] + arr[i-1][j-1] + arr[i][j-1] + arr[i+1][j-1] + arr[i+1][j];
		}
		else
		{
			return arr[i-1][j-1] + arr[i-1][j] + arr[i-1][j+1] + arr[i][j+1] + arr[i+1][j+1] + arr[i+1][j] + arr[i+1][j-1] + arr[i][j-1];
		}
	}

}

data fill_neigh(int **arr, int xlen, int ylen)
{
	int **neigh;
	neigh = (int**)malloc(sizeof(int)*xlen*ylen);
	//cout << "51" << endl;
	for(int i = 0; i < ylen; i++)
	{
		
		int *narray;
		narray = (int*)malloc(sizeof(int)*xlen);
		//cout << xlen << endl;
		for(int j = 0; j < xlen; j++)
		{
			//cout << i << " " << j << endl;
			narray[j] = count_neigh(arr, i, j, xlen, ylen);
			//cout << j << " " << narray[j] << endl;
		}
		//cout << "53" << endl;
		neigh[i] = narray;
		//cout << "54" << endl;
	}
	//cout << "55" << endl;
	return {neigh};
}

void print_solution(int **solution, string output, int xlen, int ylen)
{
	fstream file(output.c_str(), ios::out);
	
	// dump output into file
	for(int i=0; i < xlen; i++)
	{
		for (int j = 0; j < ylen; j++)
		{
			file << solution[i][j] << " " ;
		}
		file << "\n";
	}
	file.close();
}

data game(int g, int **data, int xlen, int ylen)
{
	//cout << "4" << endl;
	int **next_data;
	next_data = (int**)malloc(sizeof(int)*xlen*ylen);
	for (int k = 0; k < g; k++)
	{
		//cout << "5" << endl;
		int **neigh_data;
		neigh_data = fill_neigh(data, xlen, ylen).data;
		//cout << "5b" << endl;
		for (int i =0; i < xlen; i++)
		{
			//cout << "6" << endl;
			int *next_row;
			next_row = (int*)malloc(sizeof(int)*ylen);
			for (int j = 0; j < ylen ; j++)
			{
				//cout << "7" << endl;
				int curr_cell=2;
				if (data[i][j]==1)
				{
					//cout << "8a" << endl;
					if (neigh_data[i][j] < 2 or neigh_data[i][j] > 3)
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
					//cout << "8b" << endl;
					if (neigh_data[i][j]==3)
					{
						curr_cell=1;
					}
					else
					{
						curr_cell=0;
					}
				}
				//cout << k << " " << curr_cell << endl;
				next_row[j] = curr_cell;
			}
			next_data[i] = next_row;
		}
		data=next_data;
	}
	return {data};
}

/**
 * Prints the usage of the program.
 */
void print_usage()
{
    cerr << "The format of the input should be './program2 <input_file.txt> <output_file.txt>'\n";
	cerr << "mpirun -np p ./program2 is the name of the program.\n";
	cerr << "p is the number of processors. This should be atleast 2.\n";
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
	int xlen, ylen, g, listSize;
	int **list, **data, *row; 
	string size;
	
	fstream file(input.c_str(), ios::in);
	
	// read the input size
	file >> xlen;
	file >> ylen;
	file >> g;
	//cout << m << " " << n << " " << g << endl;
	
	list = (int**)malloc(sizeof(int)*xlen*ylen);
	data = (int**)malloc(sizeof(int)*xlen*ylen);
	
	//cout << "1" << endl;
	for (int i = 0; i < xlen ; i++)
	{
		row = (int*)malloc(sizeof(int)*ylen);
		for (int j = 0; j < ylen ; j++)
		{
			file >> row[j];
			//cout << i << " " << j << " " << row[j] << endl;
		}
		//cout << "2" << endl;
		list[i] = row;
	}
	
	//cout << "3" << endl;
	file.close();
	
	// start timer
	t_mpi_sort_begin = counter.get_ms();
	
	data = game(g, list, xlen, ylen).data;
		
	// get elapsed time
	t_mpi_sort_end_all = counter.get_ms();
	
	//print time
	cout << 0 << " " << listSize << " " << t_mpi_sort_end_all - t_mpi_sort_begin << endl;
	
	// print the sorted numbers
	print_solution(data, output, xlen, ylen);

	return 0;
}