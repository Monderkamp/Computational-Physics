#include <iostream>
#include <math.h>
#include <fstream> 

using namespace std;

double f(double x)
{
return((35.0*pow(x,4.0)+0.0*pow(x,3.0)-30.0*pow(x,2.0)-0.0*pow(x,1.0)+3.0)/8.0);	
}


double fstrich(double x)
{
return((140.0*pow(x,3.0)+0.0*pow(x,2.0)-60.0*pow(x,1.0)-0.0)/8);	
}


int main()
{
ofstream out("root1")	

double epsilon = 1.0e-12;
cout << "tolerance in f(x) if what is equal to zero = " << epsilon << endl;

		int counter = 0;
		double x0 = -0.7;		
		cout << "x0 = " << x0;	

		while (fabs(f(x0)) > epsilon)
			{
		
			x0 = x0 - f(x0)/fstrich(x0);
			counter++;
			}	
		
		cout << " root equal to " << x0;
		cout << " counter = " << counter << endl;
out.close();
getchar();
return 0;
}

