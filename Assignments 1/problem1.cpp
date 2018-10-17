#include <iostream>
#include <math.h>

using namespace std;
const double pi = M_PI;

double f(double t)
	{
		return(2.0*pow(pi,-0.5)*exp(-t*t));
	}



int main()
{
	double a = 1.0;
	double b = 0.0;
	int Ntrap = 0;
	int Nsimp = 1;
	cout << "This console application computes the integral over the function f(t) = (2/sqrt(pi)) * exp(-t^2) via the numerical quadrature methods that are: 1) the trapezoidal method, 2) the Simpson rule \n\n"; 
	
	
	while (a>b)
	{
		cout << "Insert lower integral boundary and confirm with enter. The lower boundary has to be smaller than the upper boundary.\n";
		cin >> a;
		
		
		cout << "Insert upper integral boundary and confirm with enter. \n";
		cin >> b;
	}
	
	cout << "Insert the number of nodes for the trapezoidal method and confirm with enter. \n";
	cin >> Ntrap;
	
	
	
	while (Nsimp % 2 == 1)
		{
			cout << "Insert the number of nodes for the Simpson method and confirm with enter. The number of nodes has to be an even number. \n";
			cin >> Nsimp;
		}
	
	
	cout << "lower boundary = " << a << endl;
	cout << "upper boundary = " << b << endl;
	cout << "number of nodes for the trapezoidal method =" << Ntrap << endl;
	cout << "number of nodes for the Simpson method = " << Nsimp << endl;
	
	double xtrap[Ntrap];
	double xsimp[Nsimp];
	
	
	
	double htrap = (double)(b-a)/(double)(Ntrap-1);

	double hsimp = (double)(b-a)/(double)(Nsimp-1);
	
	for (int i = 0;i<Ntrap;i++)
	{
		xtrap[i] = a + i*htrap;
	}
	
		for (int i = 0;i<Nsimp;i++)
	{
		xsimp[i] = a + i*hsimp;
	}
	
	double solutiontrap = 0.0;
	
	cout << "The numerical solution to the integral with the given parameters is equal to " << solutiontrap << endl;
	
	return 0;
}
