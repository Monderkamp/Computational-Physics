#include<iostream>
#include<cmath>
#include <fstream>
#include "rnm.cpp"

using namespace std;

int main()
    {
        ofstream out("UvsN.txt");
        srand(1);

        for (int N = 1;N<=10000;N++)
        {
            int N_c = 0;
            for (int i=0;i<N;i++)
                {
                double x = rnm()-0.5;
                double y = rnm()-0.5;

                if (x*x+y*y < 0.25)
                    {
                        N_c++;
                    }    
                }
            double pi_approx = 4.0*N_c/N;
            cout << "N = " << N << "pi_approx = " << pi_approx << ", abs(pi_approx-pi) = " << abs(pi_approx-M_PI) << endl;
            out << N << "    " << pi_approx << "    " << abs(pi_approx-M_PI) << endl;
        }
    }