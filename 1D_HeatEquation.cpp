//Delsey Sabu
//1-D PDE solver using heat equation e^-(x^2)

#include <iostream>
#include<math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

const int c= 1; 
double dx = 1;
double dt = 0.1; 
double position[8] = {-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5}; 
double result; 
double r = (c*dt)/(dx*dx); 
double k = 1-(2*r); 
void printT(double left[], int ROWS); 

using namespace std; 

int main(int argc, char** argv)
{
double left_old[6]; 
double left_new[6];
double right_new[6]; 
double right_old[6];

MPI_Init(NULL, NULL);
  // Find out rank, size
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // We are assuming at least 2 processes for this task
  if (world_size < 2) {
    fprintf(stderr, "World size must be greater than 1 for %s\n", argv[0]);
    MPI_Abort(MPI_COMM_WORLD, 1);}

if (world_rank == 0) {
	 
	for (int i = 1; i<5; i++)
	{
	result = exp(-1*(position[i-1]*position[i-1])); 
	left_old[i] = result; 
	}
    	MPI_Recv(&left_old[0], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    	MPI_Recv(&left_old[5], 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    	MPI_Send(&left_old[1], 1, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD);
    	MPI_Send(&left_old[4], 1, MPI_DOUBLE, 1, 3, MPI_COMM_WORLD); 
	printT(left_old, 6); 

}//end if 
else if (world_rank == 1) {

	for (int i = 1; i<5; i++)
	{
	result = exp(-1*(position[i+3]*position[i+3])); 
	right_old[i] = result; 
	}

    	MPI_Send(&right_old[4], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    	MPI_Send(&right_old[1], 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

    	MPI_Recv(&right_old[5], 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    	MPI_Recv(&right_old[0], 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	printT(right_old, 6); 
}//end else if 

//time step
for (int R =0; R<3; R++){
	if (world_rank == 0) {
		for (int i =1; i<5; i++)
		{
		left_new[i] = (r * left_old[i-1])+(k * left_old[i])+(r * left_old[i+1]);
		}//end for loop
	MPI_Recv(&left_new[0], 1, MPI_DOUBLE, 1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    	MPI_Recv(&left_new[5], 1, MPI_DOUBLE, 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    	MPI_Send(&left_new[1], 1, MPI_DOUBLE, 1, 6, MPI_COMM_WORLD);
    	MPI_Send(&left_new[4], 1, MPI_DOUBLE, 1, 7, MPI_COMM_WORLD); 
	printT(left_new,6); 
		for (int j =0; j< 6; j++)
		{
		left_old[j] = left_new[j];  
		}//end for loop
	}//end if 
	else if (world_rank == 1) {

		for (int i =1; i<5; i++)
		{
		right_new[i] = (r * right_old[i-1])+(k * right_old[i])+(r * right_old[i+1]);
		}//end for loop
    	MPI_Send(&right_new[4], 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
    	MPI_Send(&right_new[1], 1, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);

    	MPI_Recv(&right_new[5], 1, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    	MPI_Recv(&right_new[0], 1, MPI_DOUBLE, 0, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printT(right_new, 6); 
		for (int j =0; j<6; j++)
		{
		right_old[j] = right_new[j];  
		}//end for loop 
	}//end else if
}//end for loop 




MPI_Finalize();
return 0; 
}//end main 

//--------------------------------------------------------------------
//print fucntion 
void printT(double T[], int ROWS)
{
for (int j =1; j<ROWS-1; j++)
{cout << T[j] << " "; 
}
cout << endl; 
}
//--------------------------------------------------------------------
