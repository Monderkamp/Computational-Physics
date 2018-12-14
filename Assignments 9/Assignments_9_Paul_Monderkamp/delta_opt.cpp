#include<iostream>
#include<cmath>
#include <fstream>
#include "rnm.cpp"

using namespace std;

int main()
    {
        ofstream out("delta_opt.txt");
        srand(1);
        int N = 1e5;
        for (double delta = 0.001;delta<=1.0;delta+=0.001)
            {

                int N_c = 0;

                double x = rnm()-0.5;
                double y = rnm()-0.5;
                if (x*x+y*y < 0.25)N_c++; 
        
                for (int i=0;i<N;i++)
                    {
                        double dx = (rnm()*2.0-1.0)*delta;
                        double dy = (rnm()*2.0-1.0)*delta;
                        x += dx;
                        y += dy;
                        
                        if (abs(x)<0.5)
                            {
                                if (abs(y)<0.5)
                                    {
                                        if (x*x+y*y < 0.25)N_c++; 
                                        continue;
                                    }            
                            }
                        if (x*x+y*y < 0.25)N_c++;
                        x -= dx;
                        y -= dy;
                    }
                double pi_approx = 4.0*N_c/N;
                out << delta << "    " << pi_approx << "    " << abs(pi_approx-M_PI) << endl;
                cout << delta << endl;

            }
    }