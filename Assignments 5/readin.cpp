#include <iostream>

#include <vector>

#include <cmath>
#include <fstream>
#include <string>

typedef vector<double> vouble; 


using namespace std;


vector<double> readin(int &counter, string filename)
{
vouble omega;


ifstream in(filename);
double number;

while (in >> number)
	{
		omega.push_back(number);
		counter++;
	}

in.close();
return omega;

}

