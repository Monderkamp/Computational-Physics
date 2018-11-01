#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

int main()
{
string timeinterval = "0_005";
string outputfilename = "3_2b_eulercromer_dt_" + timeinterval + ".txt";

double tmin = 0.0;
double tmax = 50.0;
double dt = 0.005;

int N = (int)(tmax - tmin)/dt + 1;

double y[N][3];

y[0][0] = (7.0/360.0)*2.0*M_PI;
y[0][1] = 0.0;
y[0][2] = 50.0*y[0][1]*y[0][1]+100.0*(1.0-cos(y[0][0]));

ofstream out;
out.open(outputfilename);
out << y[0][0] << "  " << y[0][1] << "  " << y[0][2] << endl;

for (int i=0;i<N;i++)
  {
    y[i+1][1] = y[i][1] - dt*sin(y[i][0]);
    y[i+1][0] = y[i][0] + y[i+1][1]*dt;
    y[i+1][2] = 50.0*y[i+1][1]*y[i+1][1]+100.0*(1.0-cos(y[i+1][0]));

    out << y[i+1][0] << "  " << y[i+1][1] << "  " << y[i+1][2] << endl;
  }


out.close();
return 0;
}