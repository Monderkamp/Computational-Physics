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
        const int Nsample = 1e4;
        const int Nbins = 200;
        cout << "Number of bins: " << Nbins << endl;
        vouble g(Nbins);

        for (int i=0; i<Nbins;i++)
            {
                g[i] = 0.0;
            }    

        for (int j=0;j<Nsample;j++)
            {
                if (j % (Nsample/100) == 0) {cout << (double)j/Nsample << endl;}
                for (int i=0;i<N;i++)
                    {
                        p[i].set_x(rnm()*sideL);
                        p[i].set_y(rnm()*sideL);
                    }
        
                vouble g_imd(Nbins);
                g_imd = rdf(p,N,sideL,Nbins);
                for (int k=0;k<Nbins;k++)
                    {
                        g_imd[k] /= Nsample;
                        g[k] += g_imd[k];
                    }

            }

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