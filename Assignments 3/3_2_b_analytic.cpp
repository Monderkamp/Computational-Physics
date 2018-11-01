#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

int main()
{
ofstream out;
out.open("analyticsol.txt");
double dt = 0.01;

for (double t=0.0;t<=50.0;t+=dt)
  {
    double theta = (7.0/360.0)*2.0*M_PI*cos(t);
    double omega = -(7.0/360.0)*2.0*M_PI*sin(t);
    double E = 50.0*(omega*omega+theta*theta);
  out << t << "  " << theta  << "  " << omega  << "  " << E << endl;
  }
out.close();
return 0;
}