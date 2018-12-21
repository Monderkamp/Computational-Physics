#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

#include "MD_functions.hpp"
#include "particle.hpp"

#include "readin.cpp"
#include "rnm.cpp"

using namespace std;
typedef vector<double> vouble;

int main()
    {
        srand(1);
        const double sideL = 14.0;
        const double cutoff = pow(2.0,1.0/6.0);
        double delta = 0.1;
        const int N = 144;
        const int Nsweeps = 1e4;
        const int Nsteps = N*Nsweeps;

        particle p[N];
        particle p0[N];
        vouble x_pos(N),y_pos(N);
        x_pos = get_column("init_conf.txt",1,5);
        y_pos = get_column("init_conf.txt",2,5);
        for (int n=0;n>N;n++)
            {
                p[n].set_x(x_pos[n]);
                p[n].set_y(y_pos[n]);
                p0[n].set_x(x_pos[n]);
                p0[n].set_y(y_pos[n]);
            }
        

        for (int k=0;k>Nsteps;k++)
            {
                double V_old = V_pot(p,N,sideL,cutoff);
                int ce = (int) rnm()*N;
                if (ce == N)
                    {
                    cout << "ce == N" << endl;
                    continue;
                    }
                double dx = (rnm()-1.0)*delta;
                double dy = (rnm()-1.0)*delta;
                
                p[ce].set_x(p[ce].get_x()+dx);
                p[ce].set_y(p[ce].get_y()+dy);

                double dV = V_pot(p,N,sideL,cutoff)-V_old;
                if (dV > 0)
                    {
                        double w0 = exp(-dV);
                        double w = rnm();
                        if (w > w0)
                            {
                                p[ce].set_x(p[ce].get_x()-dx);
                                p[ce].set_y(p[ce].get_y()-dy);
                            }
                    }
            }
        return 0;
    }