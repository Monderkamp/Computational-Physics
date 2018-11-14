#include <iostream>

#include <vector>

#include <cmath>
#include <fstream>

typedef vector<double> vouble; 


using namespace std;


vector<double> readin(int &counter)
{
vouble omega;


ifstream in("noisy.txt");
double number;

while (in >> number)
	{
		omega.push_back(number);
		counter++;
	}

in.close();
return omega;

}

