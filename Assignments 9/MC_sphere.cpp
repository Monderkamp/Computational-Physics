#include <iostream>
#include <cmath>
#include <fstream>
#include "rnm.cpp"

using namespace std;

int main()
    {        
        srand(1);
        const int d = 5;
        int N = 1e8;
        int N_c = 0;
        double* x = new double[d];
        for (int i=0;i<N;i++)
            {
                double r_squared = 0.0;
                for (int j=0;j<d;j++)
                    {
                        x[j] = rnm()-0.5;
                        r_squared += x[j]*x[j];
                    }
                


                if (r_squared < pow(0.5,2) )
                    {
                        N_c++;
                    }    
            }
        double V_d = (double)N_c/N;
        cout << "V_d = " << V_d << endl;
        getchar();
        //out << N << "    " << pi_approx << "    " << abs(pi_approx-M_PI) << endl;
        delete[] x;
        
    }