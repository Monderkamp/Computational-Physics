#include <iostream>
#include <cmath>
#include <fstream> 

using namespace std;

double f1(double x)
{	
	
	return ((2.0/sqrt(M_PI))*exp(-x*x));
}



double boolsrule(double a, double b)
{
double h = (b - a)/4.0;
double x[5];
for (int i = 0; i < 5; i++)
	{
		x[i] = a + i*h;
	}

double c = (2.0*h/45.0)*(7.0*f1(x[0])+32.0*f1(x[1])+12.0*f1(x[2])+32.0*f1(x[3])+7.0*f1(x[4]));
		
return(c);	
}
