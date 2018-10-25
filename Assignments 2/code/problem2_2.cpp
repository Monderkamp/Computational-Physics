#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

const double pi = 3.14159265358979323846264338327;

int main()
{
double t0 		= 0.0;
double tmax 	= 50.0;
double dt 		= 0.2;
double x0 		= 1.0;
double v0		= 0.0;

int N = (int)(tmax - t0)/dt + 1;

double y[N][3];

y[0][0] = x0;
y[0][1] = v0;
y[0][2] = 0.5*(y[0][0]*y[0][0]+y[0][1]*y[0][1]);
ofstream out("problem_2_2.txt");

out << 0 << "	" << y[0][0] << "	" << y[0][1] << "	" << y[0][2] << endl;

for (int i = 0; i<N; i++)
	{
		y[i+1][0] = y[i][0] + dt * y[i][1];
		y[i+1][1] = y[i][1] - dt * y[i+1][0];
		y[i+1][2] = 0.5*(y[i+1][0]*y[i+1][0]+y[i+1][1]*y[i+1][1]);
		out << i+1 << "	" << y[i+1][0] << "	" << y[i+1][1] << "	" << y[i+1][2] << endl; 
		
	}

//cout << N << endl;


out.close();

return 0;	
}
