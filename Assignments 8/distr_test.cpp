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

        ofstream out("randconf.txt");        
        for (int i=0;i<N;i++)
            {
                p[i].set_x(rnm()*sideL);
                p[i].set_y(rnm()*sideL);
                out << p[i].get_x() << "    " << p[i].get_y() << endl;
                //cout << p[i].get_x() << "    " << p[i].get_y() << endl;
            }
        out.close();

        const int Nbins = 30;
        cout << Nbins << endl;
        vouble g(Nbins);

        for (int i=0; i<Nbins;i++)
            {
                g[i] = 0.0;
            }     
        //g = rdf(p,N,sideL,Nbins);

        ofstream g_out("g.txt");
        
        for (int i=0;i<Nbins;i++)
            {
                //g_out << i << "    " << g[i] << endl;
                cout << i << "    " << g[i] << endl;
            }
        g_out.close();

        delete[] p;
        getchar();
        return 0;
    }