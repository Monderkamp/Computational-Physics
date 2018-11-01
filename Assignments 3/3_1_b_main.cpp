#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

int main()
{
double tmin 	= 0.0;
double tmax 	= 50.0;
const int N	= 500;
double dt 	= (double)(tmax - tmin)/(N-1);

double y[N][3];
y[0][0] = 1.0;
y[0][1] = 0.0;
y[0][2] = 0.5*(y[0][0]*y[0][0]+y[0][1]*y[0][1]);

for (int i=0;i<N-1;i++)
{
	y[i+1][0] = y[i][0] + y[i][1]*dt -0.5*y[i][0]*dt*dt-(1.0/6.0)*y[i][1]*pow(dt,3.0)+(1.0/24.0)*y[i][0]*pow(dt,4.0);
	y[i+1][1] = y[i][1] - y[i][0]*dt -0.5*y[i][1]*dt*dt+(1.0/6.0)*y[i][0]*pow(dt,3.0)+(1.0/24.0)*y[i][1]*pow(dt,4.0);
	y[i+1][2] = 0.5*(y[i+1][0]*y[i+1][0]+y[i+1][1]*y[i+1][1]);
}
ofstream outputfile;
outputfile.open("3_1_b_results.txt");

for (int i=0; i<N;i++)
{
outputfile << y[i][0] << "	" << y[i][1] << "	" << y[i][2] << endl;
}

outputfile.close();

return 0;
}