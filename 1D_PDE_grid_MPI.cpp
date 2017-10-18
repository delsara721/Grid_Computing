//Delsey Sabu
//PDE solver 
//given a domain, we break it in half to learn the process of communciation between grids as time steps increases
//the domain is split into right and left arrays where the computaions are done separately, 
//but left and right arrays communicate to each other at the end of each time step using MPI

#include <iostream>
#include<math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int ROWS_R =  6; //the size of right array is 4 units, the outer indexes are left for communcation
int ROWS_L = 6; //the size of left array is 4 units, the outer indexes are left for communcation
const int c= 1; //constant c for teh heat equation
int ROWS_X = 8; 
double POSITION_MIN = -3.5; 
double dx = 1;
double dt = 0.1; 
double TIME_MAX = 0.3; 
double period; 

//function prototypes 
void fillPosition(double position[]); 
void printposition(double position[]); 
void fill(double left_old[],double left_new[], double right_old[],double right_new[]); 
void fillLeft(double left_old[],double left_new[], double right_old[],double right_new[]); 
void fillRight(double left_old[],double left_new[], double right_old[],double right_new[]); 
void printT(double left[], int ROWS); 
using namespace std; 

int main(int argc, char** argv)
{
bool Lahead = false; 
bool Rahead = false; 
int tag = 0; 
	//need 4 arrays, two(right_old and left_old) for keeping values from old time step
	//and two(right_new and left_new) for saving values from new time step
double right_new[ROWS_R];
double right_old[ROWS_R];  
double left_new[ROWS_L];
double left_old[ROWS_L];
//need an array to fill the positions 
double position[ROWS_X];
fillPosition(position);
//printposition(position);   
double time[4] = {0.0,0.1,0.2,0.3}; 
double result; 
period = TIME_MAX/dt; 

// Initialize the MPI environment
  MPI_Init(NULL, NULL);
  // Find out rank, size
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // We are assuming at least 2 processes for this task
  if (world_size < 2) {
    fprintf(stderr, "World size must be greater than 1 for %s\n", argv[0]);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

	//fill left and rigt ararys with initial values with cos(x), where x is position, time =0
        //first process - left array 
	if (world_rank == 0)
  {
for (int i =1; i<ROWS_L-1; i++)
{
	//for every inner index in left array, fill with cos(x), where x is position  
result = cos(position[i-1]); 
left_old[i] = result; 
}
	    //send to right array their outer values
    MPI_Send(&left_old[1], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    MPI_Send(&left_old[ROWS_L-2], 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
	    //receive from right array the left arrays outer values
    MPI_Recv(&left_old[0], 1, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&left_old[ROWS_L-1], 1, MPI_DOUBLE, 1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
printT(left_old, ROWS_L);
}
	//for second process: right array
 else if (world_rank == 1) 
{
for (int i =1; i<ROWS_R-1; i++)
{
result = cos(position[(ROWS_R-1)+(i-1)]); 
right_old[i] = result; 
}
//receive from left array the outer values
    MPI_Recv(&right_old[ROWS_R-1], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&right_old[0], 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//send to right it's outer values
    MPI_Send(&right_old[ROWS_R-2], 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
    MPI_Send(&right_old[1], 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
printT(right_old, ROWS_R); 

}


 MPI_Finalize();

	//for every time step after initialization:
for (int i =1; i<=period+1; i++)
{
cout << i << endl; 
	//fill with new values 
fill(left_old, left_new, right_old, right_new); 

cout << "LEFT" << endl; 
printT(left_new, ROWS_L); 
cout << "RIGHT" << endl; 
printT(right_new, ROWS_R); 
}

return 0; 
}

//print array function 
void printT(double T[], int ROWS)
{
for (int j =1; j<ROWS-1; j++)
{cout << T[j] << " "; 
}
cout << endl; 
}

//fill with new values
void fill(double left_old[],double left_new[], double right_old[],double right_new[])
{ 
	double r = (c * dt)/(dx * dx); 
	double k= 1-(2 * r); 
	
	for (int i =1; i<ROWS_L-1; i++)
	{
	left_new[i] = (r * left_old[i-1])+(k * left_old[i])+(r * left_old[i+1]);
	}//end for loop

	for (int i =1; i<ROWS_R-1; i++)
	{
	right_new[i] = (r * right_old[i-1])+(k * right_old[i])+(r * right_old[i+1]);
	}//end for loop

left_new[0] = right_new[ROWS_R-2]; 
left_new[ROWS_L-1] = right_new[1]; 
right_new[ROWS_R-1] = left_new[1]; 
right_new[0] = left_new[ROWS_L-2]; 

for (int j =0; j< ROWS_L; j++)
	{
	left_old[j] = left_new[j];  
	}//end for loop 

for (int j =0; j<ROWS_R; j++)
	{
	right_old[j] = right_new[j];  
	}//end for loop 
}//end function 

void fillLeft(double left_old[],double left_new[], double right_old[],double right_new[])
{ 
	double r = (c * dt)/(dx * dx); 
	double k= 1-(2 * r); 
	
	for (int i =1; i<ROWS_L-1; i++)
	{
	left_new[i] = (r * left_old[i-1])+(k * left_old[i])+(r * left_old[i+1]);
	}//end for loop


left_new[0] = right_new[ROWS_R-2]; 
left_new[ROWS_L-1] = right_new[1]; 
 

for (int j =0; j< ROWS_L; j++)
	{
	left_old[j] = left_new[j];  
	}//end for loop 

}//end function 

void fillRight(double left_old[],double left_new[], double right_old[],double right_new[])
{ 
	double r = (c * dt)/(dx * dx); 
	double k= 1-(2 * r); 
	

for (int i =1; i<ROWS_R-1; i++)
	{
	right_new[i] = (r * right_old[i-1])+(k * right_old[i])+(r * right_old[i+1]);
	}//end for loop

right_new[ROWS_R-1] = left_new[1]; 
right_new[0] = left_new[ROWS_L-2]; 


for (int j =0; j<ROWS_R; j++)
	{
	right_old[j] = right_new[j];  
	}//end for loop 
}//end function 


void fillPosition(double position[])
{
	double i =POSITION_MIN;
	for (int counter = 0; counter<ROWS_X; counter++)
	{
	position[counter] = i; 
	i = i+dx; 	 
	}//end for loop
}//end function fillPOsition

void printposition(double position[])
{
for (int counter = 0; counter<ROWS_X; counter++)
	{
	cout << position[counter] << endl; 	 
	}//end for loop
}

