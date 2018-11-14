#include <iostream>

#include <vector>

#include <cmath>
#include <fstream>
#include <string>

using namespace std;


#include "readin.cpp"


class particle 
	{
		public:
		particle()
			{
				x = 0.0;
				y = 0.0;
				vx = 0.0;
				vy = 0.0;
			}
		~particle(){}

		void setx(const double a);
		void sety(const double a);
		void setvx(const double a);
		void setvy(const double a);

		double getx();
		double gety();
		double getvx();
		double getvy();


		private:
		double x,y,vx,vy;
	};

int main()
{
vouble a,b;
int n = 0;
int m = 0;

a = readin(n,"veldat2.txt");
cout << n << endl;

b = readin(m,"veldat2.txt");
cout << m << endl;

getchar();
return 0;
}

void particle::setx(const double a){x = a;}
void particle::sety(const double a){y = a;}
void particle::setvx(const double a){vx = a;}
void particle::setvy(const double a){vy = a;}

double particle::getx(){return x;}
double particle::gety(){return y;}
double particle::getvx(){return vx;}
double particle::getvy(){return vy;}
