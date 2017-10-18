#include<iostream>
#include <math.h>
using namespace std; 

const int ROWS_X = 4; 
const int ROWS_Y = 4; 
const double c = 1.0;
const double dx = 1.0; 
const double dy = 1.0; 
const double dt = 0.1;  
const double MAX_TIME = 0.3; 
double X[4]= {-1.5,-0.5,0.5,1.5}; 
double Y[4]={-1.5,-0.5,0.5,1.5}; 

void initialize(double T_old[ROWS_X][ROWS_Y], double T_new[ROWS_X][ROWS_Y]); 
void printT(double T[ROWS_X][ROWS_Y]); 
void fill(double ul_new[ROWS_X][ROWS_Y], double ul_old[ROWS_X][ROWS_Y], double ur_new[ROWS_X][ROWS_Y], double ur_old[ROWS_X][ROWS_Y],double dl_new[ROWS_X][ROWS_Y], double dl_old[ROWS_X][ROWS_Y], double dr_new[ROWS_X][ROWS_Y], double dr_old[ROWS_X][ROWS_Y]); 

int main()
{
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

//ul
for (int i = 1; i <ROWS_X-1; i++)
{
for (int j =1; j<ROWS_Y-1; j++)
{
ul_old[i][j] = exp(-1*((X[i-1]*X[i-1])+(Y[j-1]*Y[j-1])));
ur_old[i][j] = exp(-1*((X[i-1]*X[i-1])+(Y[j+1]*Y[j+1])));
dl_old[i][j] = exp(-1*((X[i+1]*X[i+1])+(Y[j-1]*Y[j-1])));
dr_old[i][j] = exp(-1*((X[i+1]*X[i+1])+(Y[j+1]*Y[j+1])));


ul_old[0][j] = dl_old[ROWS_X-2][j];
ul_old[ROWS_X-1][j] = dl_old[1][j]; 
ur_old[0][j] = dr_old[ROWS_X-2][j];
ur_old[ROWS_X-1][j] = dr_old[1][j]; 
dl_old[0][j] = ul_old[ROWS_X-2][j]; 
dl_old[ROWS_X-1][j] = ul_old[1][j]; 
dr_old[0][j] = ur_old[ROWS_X-2][j];
dr_old[ROWS_X-1][j] = ur_old[1][j];
}
ul_old[i][0] = ur_old[i][ROWS_Y-2];
ul_old[i][ROWS_Y-1] = ur_old[i][1]; 
dl_old[i][0] = dr_old[i][ROWS_Y-2]; 
dl_old[i][ROWS_Y-1] = dr_old[i][1]; 
ur_old[i][0] = ul_old[i][ROWS_Y-2];
ur_old[i][ROWS_Y-1] = ul_old[i][1]; 
dr_old[i][0] = dl_old[i][ROWS_Y-2]; 
dr_old[i][ROWS_Y-1] = dl_old[i][1];
}


printT(ul_old);
printT(ur_old);
printT(dl_old); 
printT(dr_old);

double period = MAX_TIME/dt; 


for (int i =0; i<period; i++)
{ 
fill(ul_new, ul_old, ur_new, ur_old, dl_new, dl_old, dr_new,dr_old); 
printT(ul_new); 
printT(ur_new); 
printT(dl_new); 
printT(dr_new); 
}

return 0; 
}

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


void fill(double ul_new[ROWS_X][ROWS_Y], double ul_old[ROWS_X][ROWS_Y], double ur_new[ROWS_X][ROWS_Y], double ur_old[ROWS_X][ROWS_Y],double dl_new[ROWS_X][ROWS_Y], double dl_old[ROWS_X][ROWS_Y], double dr_new[ROWS_X][ROWS_Y], double dr_old[ROWS_X][ROWS_Y])
{
double r = (c*dt)/(dx*dx); 
double s = (c*dt)/(dy*dy); 
for (int i =1; i<ROWS_X-1; i++)
{
for (int j = 1; j<ROWS_Y-1; j++)
{
ul_new[i][j] = r*ul_old[i-1][j]+r*ul_old[i+1][j]+(1-2*r-2*s)*ul_old[i][j]+s*ul_old[i][j-1]+s*ul_old[i][j+1]; 
ur_new[i][j] = r*ur_old[i-1][j]+r*ur_old[i+1][j]+(1-2*r-2*s)*ur_old[i][j]+s*ur_old[i][j-1]+s*ur_old[i][j+1]; 
dl_new[i][j] = r*dl_old[i-1][j]+r*dl_old[i+1][j]+(1-2*r-2*s)*dl_old[i][j]+s*dl_old[i][j-1]+s*dl_old[i][j+1]; 
dr_new[i][j] = r*dr_old[i-1][j]+r*dr_old[i+1][j]+(1-2*r-2*s)*dr_old[i][j]+s*dr_old[i][j-1]+s*dr_old[i][j+1]; 
ul_new[0][j] = dl_new[ROWS_X-2][j]; 
ul_new[ROWS_X-1][j] = dl_new[1][j]; 
ur_new[0][j] = dr_new[ROWS_X-2][j];
ur_new[ROWS_X-1][j] = dr_new[1][j]; 

dl_new[0][j] = ul_new[ROWS_X-2][j]; 
dl_new[ROWS_X-1][j] = ul_new[1][j]; 
dr_new[0][j] = ur_new[ROWS_X-2][j];
dr_new[ROWS_X-1][j] = ur_new[1][j];
}//end inside for loop
ul_new[i][0] = ur_new[i][ROWS_Y-2]; 
ul_new[i][ROWS_Y-1] = ur_new[i][1]; 
dl_new[i][0] = dr_new[i][ROWS_Y-2]; 
dl_new[i][ROWS_Y-1] = dr_new[i][1]; 

ur_new[i][0] = ul_new[i][ROWS_Y-2]; 
ur_new[i][ROWS_Y-1] = ul_new[i][1]; 
dr_new[i][0] = dl_new[i][ROWS_Y-2]; 
dr_new[i][ROWS_Y-1] = dl_new[i][1];
}//end outside for loop

//copy new to old
for (int i =0; i<ROWS_X; i++)
{
for (int j = 0; j<ROWS_Y; j++)
{
	ul_old[i][j] = ul_new[i][j]; 
	ur_old[i][j] = ur_new[i][j]; 
	dl_old[i][j] = dl_new[i][j]; 
	dr_old[i][j] = dr_new[i][j]; 
}//end inside for loop
}//end outside for loop

}//end fill function


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
