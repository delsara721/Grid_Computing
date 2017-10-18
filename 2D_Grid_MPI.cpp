//Delsey Sabu
//2-D PDE solver 
//the domain is broken into 4 different arrays where work is done, but the 4 communicate with each other after every timestep (using MPI)
#include <iostream>
#include<math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std; 

const int ROWS_X = 4; 
const int ROWS_Y = 4; 
const double c = 1.0;
const double dx = 1.0; 
const double dy = 1.0; 
const double dt = 0.1;  
const double r = (c*dt)/(dx*dx); 
const double s = (c*dt)/(dy*dy); 
const double MAX_TIME = 0.3; 
double X[4]= {-1.5,-0.5,0.5,1.5}; 
double Y[4]={-1.5,-0.5,0.5,1.5}; 
const double period = MAX_TIME/dt; 
void initialize(double T_old[ROWS_X][ROWS_Y], double T_new[ROWS_X][ROWS_Y]); 
void printT(double T[ROWS_X][ROWS_Y]); 
void fillul(double ul_new[ROWS_X][ROWS_Y], double ul_old[ROWS_X][ROWS_Y]); 
void fillur(double ur_new[ROWS_X][ROWS_Y], double ur_old[ROWS_X][ROWS_Y]); 
void filldl(double dl_new[ROWS_X][ROWS_Y], double dl_old[ROWS_X][ROWS_Y]); 
void filldr(double dr_new[ROWS_X][ROWS_Y], double dr_old[ROWS_X][ROWS_Y]); 
//////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
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
//make 8 arrays
double ul_new[ROWS_X][ROWS_Y]; 
double ul_old[ROWS_X][ROWS_Y];
double ur_new[ROWS_X][ROWS_Y]; 
double ur_old[ROWS_X][ROWS_Y];
double dl_new[ROWS_X][ROWS_Y]; 
double dl_old[ROWS_X][ROWS_Y];
double dr_new[ROWS_X][ROWS_Y]; 
double dr_old[ROWS_X][ROWS_Y];

//intialize everything to zero
initialize(ul_new,ul_old); 
initialize(ur_new,ur_old); 
initialize(dl_new,dl_old); 
initialize(dr_new,dr_old);

if (world_rank == 0) //upper left:ul
{
for (int i = 1; i <ROWS_X-1; i++)
{
for (int j =1; j<ROWS_Y-1; j++)
{
ul_old[i][j] = exp(-1*((X[i-1]*X[i-1])+(Y[j-1]*Y[j-1])));
MPI_Recv(&ul_old[0][j], 1, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //0
MPI_Recv(&ul_old[ROWS_X-1][j], 1, MPI_DOUBLE, 2, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //0
MPI_Send(&ul_old[ROWS_X-2][j], 1, MPI_DOUBLE, 2, 4, MPI_COMM_WORLD); //0
MPI_Send(&ul_old[1][j], 1, MPI_DOUBLE, 2, 5, MPI_COMM_WORLD); //0

}

MPI_Recv(&ul_old[i][0], 1, MPI_DOUBLE, 1, 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //0
MPI_Recv(&ul_old[i][ROWS_Y-1], 1, MPI_DOUBLE, 1, 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //0
MPI_Send(&ul_old[i][ROWS_Y-2], 1, MPI_DOUBLE, 1, 12, MPI_COMM_WORLD); //0
MPI_Send(&ul_old[i][1], 1, MPI_DOUBLE, 1, 13, MPI_COMM_WORLD); //0 
}
printT(ul_old);

} //end if 

else if (world_rank ==1) //upper right: ur
{
for (int i = 1; i <ROWS_X-1; i++)
{
for (int j =1; j<ROWS_Y-1; j++)
{
ur_old[i][j] = exp(-1*((X[i-1]*X[i-1])+(Y[j+1]*Y[j+1])));
MPI_Recv(&ur_old[0][j], 1, MPI_DOUBLE, 3, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //1
MPI_Recv(&ur_old[ROWS_X-1][j], 1, MPI_DOUBLE, 3, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //1
MPI_Send(&ur_old[ROWS_X-2][j], 1, MPI_DOUBLE, 3, 6, MPI_COMM_WORLD); //1
MPI_Send(&ur_old[1][j], 1, MPI_DOUBLE, 3, 7, MPI_COMM_WORLD); //1
}

MPI_Send(&ur_old[i][ROWS_Y-2], 1, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD); //1
MPI_Send(&ur_old[i][1], 1, MPI_DOUBLE, 0, 9, MPI_COMM_WORLD); //1 
MPI_Recv(&ur_old[i][0], 1, MPI_DOUBLE, 0, 12, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //1
MPI_Recv(&ur_old[i][ROWS_Y-1], 1, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //1
}
printT(ur_old);
}//end else if

else if (world_rank ==2) //down left: dl
{
for (int i = 1; i <ROWS_X-1; i++)
{
for (int j =1; j<ROWS_Y-1; j++)
{
dl_old[i][j] = exp(-1*((X[i+1]*X[i+1])+(Y[j-1]*Y[j-1])));
MPI_Send(&dl_old[ROWS_X-2][j], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); //2
MPI_Send(&dl_old[1][j], 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD); //2
MPI_Recv(&dl_old[0][j], 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //2
MPI_Recv(&dl_old[ROWS_X-1][j], 1, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //2
}

MPI_Recv(&dl_old[i][0], 1, MPI_DOUBLE, 3, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //2
MPI_Recv(&dl_old[i][ROWS_Y-1], 1, MPI_DOUBLE, 3, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //2
MPI_Send(&dl_old[i][ROWS_Y-2], 1, MPI_DOUBLE, 3, 14, MPI_COMM_WORLD); //2
MPI_Send(&dl_old[i][1], 1, MPI_DOUBLE, 3, 15, MPI_COMM_WORLD); //2
}
printT(dl_old); 
} //end else if 

else if (world_rank==3) //down right: dr
{
for (int i = 1; i <ROWS_X-1; i++)
{
for (int j =1; j<ROWS_Y-1; j++)
{
dr_old[i][j] = exp(-1*((X[i+1]*X[i+1])+(Y[j+1]*Y[j+1])));
MPI_Send(&dr_old[ROWS_X-2][j], 1, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD); //3
MPI_Send(&dr_old[1][j], 1, MPI_DOUBLE, 1, 3, MPI_COMM_WORLD); //3
MPI_Recv(&dr_old[0][j], 1, MPI_DOUBLE, 1, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //3
MPI_Recv(&dr_old[ROWS_X-1][j], 1, MPI_DOUBLE, 1, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //3
}
 
MPI_Send(&dr_old[i][ROWS_Y-2], 1, MPI_DOUBLE, 2, 10, MPI_COMM_WORLD); //3
MPI_Send(&dr_old[i][1], 1, MPI_DOUBLE, 2, 11, MPI_COMM_WORLD); //3
MPI_Recv(&dr_old[i][0], 1, MPI_DOUBLE, 2, 14, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //3
MPI_Recv(&dr_old[i][ROWS_Y-1], 1, MPI_DOUBLE, 2, 15, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //3
}
printT(dr_old);
} //end else if 


//for time steps
for (int i =0; i<period; i++)
{ 

if (world_rank == 0) //upper left:ul
{
fillul(ul_new, ul_old); 
printT(ul_new); 
}

else if (world_rank ==1) //upper right: ur
{
fillur(ur_new, ur_old); 
printT(ur_new); 
}

else if (world_rank ==2) //down left: dl
{
filldl(dl_new, dl_old); 
printT(dl_new); 
}

else if (world_rank==3) //down right: dr
{
filldr( dr_new,dr_old); 
printT(dr_new); 
}
}


MPI_Finalize();
return 0; 
}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void printT(double T[ROWS_X][ROWS_Y])
{
for (int i =1; i<ROWS_X-1; i++)
{
for (int j = 1; j<ROWS_Y-1; j++)
{
	cout << T[i][j] << " "; 
}//end inside for loop
cout << endl; 
}//end outside for loop
cout << endl; 
	
}//end printT function 
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void initialize(double T_old[ROWS_X][ROWS_Y], double T_new[ROWS_X][ROWS_Y])
{
for (int i =0; i<ROWS_X; i++)
{
for (int j = 0; j<ROWS_Y; j++)
{
	T_new[i][j] = 0.0; 
	T_old[i][j] = 0.0; 
}//end inside for loop
}//end outside for loop

} //end fucntion 
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void fillul(double ul_new[ROWS_X][ROWS_Y], double ul_old[ROWS_X][ROWS_Y])
{ 
for (int i =1; i<ROWS_X-1; i++)
{
for (int j = 1; j<ROWS_Y-1; j++)
{
ul_new[i][j] = r*ul_old[i-1][j]+r*ul_old[i+1][j]+(1-2*r-2*s)*ul_old[i][j]+s*ul_old[i][j-1]+s*ul_old[i][j+1]; 


MPI_Recv(&ul_new[0][j], 1, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //0
MPI_Recv(&ul_new[ROWS_X-1][j], 1, MPI_DOUBLE, 2, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //0
MPI_Send(&ul_new[ROWS_X-2][j], 1, MPI_DOUBLE, 2, 4, MPI_COMM_WORLD); //0
MPI_Send(&ul_new[1][j], 1, MPI_DOUBLE, 2, 5, MPI_COMM_WORLD); //0
}//end inside for loop

MPI_Recv(&ul_new[i][0], 1, MPI_DOUBLE, 1, 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //0
MPI_Recv(&ul_new[i][ROWS_Y-1], 1, MPI_DOUBLE, 1, 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //0
MPI_Send(&ul_new[i][ROWS_Y-2], 1, MPI_DOUBLE, 1, 12, MPI_COMM_WORLD); //0
MPI_Send(&ul_new[i][1], 1, MPI_DOUBLE, 1, 13, MPI_COMM_WORLD); //0 
}//end outside for loop

//copy new to old
for (int i =0; i<ROWS_X; i++)
{
for (int j = 0; j<ROWS_Y; j++)
{
	ul_old[i][j] = ul_new[i][j];  
}//end inside for loop
}//end outside for loop

}//end fillul function
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void fillur(double ur_new[ROWS_X][ROWS_Y], double ur_old[ROWS_X][ROWS_Y])
{

for (int i =1; i<ROWS_X-1; i++)
{
for (int j = 1; j<ROWS_Y-1; j++)
{

ur_new[i][j] = r*ur_old[i-1][j]+r*ur_old[i+1][j]+(1-2*r-2*s)*ur_old[i][j]+s*ur_old[i][j-1]+s*ur_old[i][j+1]; 

MPI_Recv(&ur_new[0][j], 1, MPI_DOUBLE, 3, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //1
MPI_Recv(&ur_new[ROWS_X-1][j], 1, MPI_DOUBLE, 3, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //1
MPI_Send(&ur_new[ROWS_X-2][j], 1, MPI_DOUBLE, 3, 6, MPI_COMM_WORLD); //1
MPI_Send(&ur_new[1][j], 1, MPI_DOUBLE, 3, 7, MPI_COMM_WORLD); //1
}//end inside for loop

MPI_Send(&ur_new[i][ROWS_Y-2], 1, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD); //1
MPI_Send(&ur_new[i][1], 1, MPI_DOUBLE, 0, 9, MPI_COMM_WORLD); //1 
MPI_Recv(&ur_new[i][0], 1, MPI_DOUBLE, 0, 12, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //1
MPI_Recv(&ur_new[i][ROWS_Y-1], 1, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //1
}//end outside for loop

//copy new to old
for (int i =0; i<ROWS_X; i++)
{
for (int j = 0; j<ROWS_Y; j++)
{ 
	ur_old[i][j] = ur_new[i][j]; 
	
}//end inside for loop
}//end outside for loop

}//end fillur function
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void filldl(double dl_new[ROWS_X][ROWS_Y], double dl_old[ROWS_X][ROWS_Y])
{
double r = (c*dt)/(dx*dx); 
double s = (c*dt)/(dy*dy); 
for (int i =1; i<ROWS_X-1; i++)
{
for (int j = 1; j<ROWS_Y-1; j++)
{

dl_new[i][j] = r*dl_old[i-1][j]+r*dl_old[i+1][j]+(1-2*r-2*s)*dl_old[i][j]+s*dl_old[i][j-1]+s*dl_old[i][j+1]; 

MPI_Send(&dl_new[ROWS_X-2][j], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); //2
MPI_Send(&dl_new[1][j], 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD); //2
MPI_Recv(&dl_new[0][j], 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //2
MPI_Recv(&dl_new[ROWS_X-1][j], 1, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //2
}//end inside for loop

MPI_Recv(&dl_new[i][0], 1, MPI_DOUBLE, 3, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //2
MPI_Recv(&dl_new[i][ROWS_Y-1], 1, MPI_DOUBLE, 3, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //2
MPI_Send(&dl_new[i][ROWS_Y-2], 1, MPI_DOUBLE, 3, 14, MPI_COMM_WORLD); //2
MPI_Send(&dl_new[i][1], 1, MPI_DOUBLE, 3, 15, MPI_COMM_WORLD); //2
}//end outside for loop

//copy new to old
for (int i =0; i<ROWS_X; i++)
{
for (int j = 0; j<ROWS_Y; j++)
{
	
	dl_old[i][j] = dl_new[i][j]; 

}//end inside for loop
}//end outside for loop

}//end filldl function
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void filldr(double dr_new[ROWS_X][ROWS_Y], double dr_old[ROWS_X][ROWS_Y])
{
 
for (int i =1; i<ROWS_X-1; i++)
{
for (int j = 1; j<ROWS_Y-1; j++)
{

dr_new[i][j] = r*dr_old[i-1][j]+r*dr_old[i+1][j]+(1-2*r-2*s)*dr_old[i][j]+s*dr_old[i][j-1]+s*dr_old[i][j+1]; 

MPI_Send(&dr_new[ROWS_X-2][j], 1, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD); //3
MPI_Send(&dr_new[1][j], 1, MPI_DOUBLE, 1, 3, MPI_COMM_WORLD); //3
MPI_Recv(&dr_new[0][j], 1, MPI_DOUBLE, 1, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //3
MPI_Recv(&dr_new[ROWS_X-1][j], 1, MPI_DOUBLE, 1, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //3
}//end inside for loop
MPI_Send(&dr_new[i][ROWS_Y-2], 1, MPI_DOUBLE, 2, 10, MPI_COMM_WORLD); //3
MPI_Send(&dr_new[i][1], 1, MPI_DOUBLE, 2, 11, MPI_COMM_WORLD); //3
MPI_Recv(&dr_new[i][0], 1, MPI_DOUBLE, 2, 14, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //3
MPI_Recv(&dr_new[i][ROWS_Y-1], 1, MPI_DOUBLE, 2, 15, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //3
}//end outside for loop

//copy new to old
for (int i =0; i<ROWS_X; i++)
{
for (int j = 0; j<ROWS_Y; j++)
{
	dr_old[i][j] = dr_new[i][j]; 
}//end inside for loop
}//end outside for loop

}//end filldr function
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
