#include <iostream>
#include <cmath>
#include <fstream> 

using namespace std;

double f(double u)
{	
	double x = 2.0;
	return ((1.0/sqrt(M_PI))*exp(-(0.5*x*(1.0+u))*(0.5*x*(1.0+u)))*x);
}



double quadrature()
{
//returns the result of the integral from -1 to 1

double a = (0.347854845137454*f(-0.861136311594053)+
		0.652145154862546*f(-0.339981043584856)+
		0.652145154862546*f(0.339981043584856) +
		0.347854845137454*f(0.861136311594053));
		
return(a);	
}
