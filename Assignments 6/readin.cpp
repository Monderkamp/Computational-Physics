#include <iostream>

#include <vector>

#include <cmath>
#include <fstream>
#include <string>

using namespace std;


typedef vector<double> vouble;
 
vouble readin(string filename)
{
	vouble omega;


	ifstream in(filename);
	double number;

	while (in >> number)
		{
			omega.push_back(number);
		}

	in.close();
	return omega;

	}



vouble get_column(string filename,int i,int Ncol)
	{
		vouble a = readin(filename);
		int Nrows = a.size()/Ncol;
		vouble b(Nrows);
		for (int j=0;j<Nrows;j++)
			{
				b[j]=a[Ncol*j+i];
			}
		return b;
	}