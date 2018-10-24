#include <iostream>
#include <math.h>
 

using namespace std;

double f(double x)
{
return(3.0*pow(x,4.0)+4.0*pow(x,3.0)-1.0*pow(x,2.0)-2.0*pow(x,1.0));	
}


double fstrich(double x)
{
return(12.0*pow(x,3.0)+12.0*pow(x,2.0)-2.0*pow(x,1.0)-2.0);	
}


int main()
{
double dx = 0.1;
double epsilon = 1.0e-12;
cout << "tolerance in f(x) if what is equal to zero = " << epsilon << endl;

for (double x = -4.0; x <= 4.0; x += dx)
	{
		int counter = 0;
		double x0 = x;		
		cout << "x0 = " << x0;	

		while (fabs(f(x0)) > epsilon)
			{
		
			x0 = x0 - f(x0)/fstrich(x0);
			counter++;
			}	
		
		cout << " root equal to " << x0;
		cout << " counter = " << counter << endl;
	}
getchar();
return 0;
}

/*
The observations from starting the
Newton algorithm from different nodes that are equidistantly spaced on the [-4,4] are:

- the major part of the samples in the lower half of the interval converges towards the root -1,
when taken as starting points.
- the major part of the samples in the upper half of the interval converges towards the root 2/3, 
when taken as starting points.
-although having the shortest interval of starting points that converge towards this root,
the iterations that converge towards x = 0 tend to converge faster. 
- the intervals are not exclusive such that there are several nodes within a row 
of sampled nodes that converge towards different roots than the surrounding nodes.

*/