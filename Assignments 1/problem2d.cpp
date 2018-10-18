#include <iostream>
#include <math.h>
 

using namespace std;

double f(double x)
{
return(3.0*pow(x,4.0)+4.0*pow(x,3.0)-1.0*pow(x,2.0)-2.0*pow(x,1.0));	
}

int main()
{
double dx = 0.1;
double epsilon = 1.0e-12;
double h = 0.01;

cout << "tolerance in f(x) if what is equal to zero = " << epsilon << endl;

for (double x = -4.0; x <= 4.0; x += dx)
	{
		int counter = 0;
		double x0 = x;		
		cout << "x0 = " << x0;	

		while (fabs(f(x0)) > epsilon)
			{
			double fstrich = (f(x0+h)-f(x0))/h;
			x0 = x0 - f(x0)/fstrich;
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
- the behaviour of the limits of the convergence in the interval seem to be similar to 
exercise 2c where instead of approximating the derivative via a difference quotient,
the derivative is evaluated analytically. '
- The convergence speed is for the roots 0 and 2/3 gone up by a single digit factor, 
whereas the number of iterations it takes to reach the root x = -1 is increased by 
a factor 10^3. 
- Hence the convergence is certainly not quadratic, the analytical proof of which 
can be found in the lecture notes. 


*/