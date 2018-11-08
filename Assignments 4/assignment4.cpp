#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
struct particle
	{
		double x,v,F;
	};

int main()
{
double dt = 0.1;
const int N = 100;
const int Nsteps = 1e3;

particle part[N];  			// e.g. part[1].x;

for (int i=0;i < N;i++)
	{
		part[i].x = i;
		part[i].v = 0.5*(2.0*(double)rand()/RAND_MAX-1.0);
		
	}



double T[Nsteps];
double V[Nsteps];
double P[Nsteps];
double E[Nsteps];

for (int i=0;i < N-1;i++)
	{
		T[0]+= 0.5*part[i].v*part[i].v;
		V[0]+= (part[i].x-part[i+1].x)*(part[i].x-part[i+1].x);
		P[0]+=part[i].v;
	}

T[0]+= 0.5*part[N-1].v*part[N-1].v;
V[0]+= (part[N-1].x-N-part[0].x)*(part[N-1].x-N-part[0].x);
P[0]+=part[N-1].v;
E[0]=T[0]+V[0];

for (int i=1;i < Nsteps;i++)
	{
		for (int j=1;j < N-1;j++)
			{
				part[j].F = part[j+1].x-2.0*part[j].x+part[j-1].x;
			}
		
		part[0].F = part[1].x-2.0*part[0].x + part[N-1].x-N;
		part[N-1].F = part[0].x*N+ part[N-2].x-2*part[N-1].x;

		for (int j=1;j < N-1;j++)
			{
				part[j].x = part[j].x + part[j].v *dt + 0.5*part[j].F*dt*dt;
				part[j].v = part[j].v + 0.5*dt*(part[j].F+part[j+1].x+part[j-1].x-2.0*part[j].x);
				//cout << part[j].v << endl;
			}

		for (int k=0;k < N-1;k++)
			{
				T[i]+= 0.5*part[k].v*part[k].v;
				V[i]+= (part[k].x-part[k+1].x)*(part[k].x-part[k+1].x);
				P[i]+=part[k].v;
			}

		T[i]+= 0.5*part[N-1].v*part[N-1].v;
		V[i]+= (part[N-1].x-N-part[0].x)*(part[N-1].x-N-part[0].x);
		P[i]+=part[N-1].v;
		E[i]=T[0]+V[0];

		//cout << T[i] << endl;
	}

ofstream out("output.txt");

for (int i=0;i<Nsteps;i++)
	{
		out << i << "	" << T[i] << "	" << V[i] << "	" << P[i] << "	" << E[i] << endl;
	}
out.close();
getchar();
return 0;
}