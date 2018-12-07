#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>

using namespace std;

#include "readin.cpp"
#include "particle.hpp"
#include "MD_functions.hpp"
#include "rnm.cpp"

const double sideL = 14.0;

typedef vector<double> vouble;

int main()
    {
        srand(1);
        const int N = 144;
        particle* p = new particle[N];
        //const int Nsample = 1e4;
        const int Nbins = 500;
        cout << "Number of bins: " << Nbins << endl;
        vouble g(Nbins);

        for (int i=0; i<Nbins;i++)
            {
                g[i] = 0.0;
            }    

        vouble x_pos(N);
        vouble y_pos(N);
        x_pos = get_column("init_conf.txt",1,5);
        y_pos = get_column("init_conf.txt",2,5);

        for (int i=0;i<N;i++)
            {
                p[i].set_x(x_pos[i]);
                p[i].set_y(y_pos[i]);
            }
        

        g = rdf(p,N,sideL,Nbins);

        ofstream g_out("g.txt");
        
        for (int i=0;i<Nbins;i++)
            {
                g_out << i << "    " << (i+0.5)*0.5*sideL/Nbins << "    " << g[i] << endl;
                //cout << i << "    " << (i+0.5)*0.5*sideL/Nbins << "    " << g[i] << endl;
            }
        g_out.close();

        delete[] p;
        //getchar();
        return 0;
    }