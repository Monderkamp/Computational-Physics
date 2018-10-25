#include <iostream>
#include <cmath>
#include <fstream> 

using namespace std;

double f(double x)
{
return((35.0*pow(x,4.0)+0.0*pow(x,3.0)-30.0*pow(x,2.0)-0.0*pow(x,1.0)+3.0)/8.0);	
}


double fstrich(double x)
{
return((140.0*pow(x,3.0)+0.0*pow(x,2.0)-60.0*pow(x,1.0)-0.0)/8.0);	
}


int main()
{
ofstream out("root4.txt");	
double actualroot = 0.861136311594053;
//ddouble actualroot = 0.339981043584856;

double epsilon = 1.0e-15;
cout << "tolerance in f(x) if what is equal to zero = " << epsilon << endl;

		int counter = 0;
		double x0 = 0.7;		
			

		while (fabs(f(x0)) > epsilon)
			{
			out << counter << "	" << fabs(x0 - actualroot) << endl;
			cout << counter << "	" << fabs(x0 - actualroot) << endl;
			cout << "x0 = " << x0 << endl;
			//cout << fabs(f(x0)) << endl;
			x0 = x0 - f(x0)/fstrich(x0);
			counter++;
			}	
		
out.close();
getchar();
return 0;
}

