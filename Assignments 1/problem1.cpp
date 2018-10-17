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
	
	double htrap = (double)(b-a)/(double)(Ntrap);
	double hsimp = (double)(b-a)/(double)(Nsimp);
	
	for (int i = 0;i<Ntrap;i++)
		{
			xtrap[i] = a + i*htrap;
		}
	
	for (int i = 0;i<Nsimp;i++)
		{
			xsimp[i] = a + i*hsimp;
		}
	
	double solutiontrap = 0.0;
	
	for (int i = 0; i<Ntrap;i++)
		{
			solutiontrap += (f(xtrap[i]) + f(xtrap[i+1]));
		}
		
	solutiontrap *= (0.5*htrap);
	
	cout << "The numerical solution with the trapezoidal method to the integral with the given parameters is equal to " << solutiontrap << endl;
		if (a == 0.0 && b == 0.4)
		{
			cout << "The exact value of the integral withing the given boundaries is: erf(0.4) = 0.428392" << endl;\
			double u 	= fabs(solutiontrap - 0.428392);
			double urel = fabs(solutiontrap - 0.428392)/0.428392;
			
			cout << "The corresponding absolute error is equal to " << u;
			cout << " whereas the relative error is equal to " << urel << endl;
		}
	
	double solutionsimp = 0.0;
	double term		= 0.0;
	
	for (int i = 1; i<Nsimp/2;i++)
		{
			term += f(xsimp[2*i]);
		}
	
	term	 		*= 2.0;
	
	solutionsimp 	= f(xsimp[0]) + term;
	
	term			= 0.0;
	
	for (int i = 1; i<=Nsimp/2;i++)
		{
			term += f(xsimp[2*i-1]);
		}
	
	term	 		*= 4.0;
	
	solutionsimp += (term + f(xsimp[Nsimp]));
	solutionsimp *= (hsimp/3.0);
	
	cout << "\nThe numerical solution with the Simpson method to the integral with the given parameters is equal to " << solutionsimp << endl;
	if (a == 0.0 && b == 0.4)
		{
			cout << "The exact value of the integral withing the given boundaries is: erf(0.4) = 0.428392" << endl;\
			double u 	= fabs(solutionsimp - 0.428392);
			double urel = fabs(solutionsimp - 0.428392)/0.428392;
				
			cout << "The corresponding absolute error is equal to " << u;
			cout << " whereas the relative error is equal to " << urel << endl;
		
		}
	return 0;
}
