#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>

#include "MD_functions.hpp"
#include "particle.hpp"

#include "readin.cpp"
#include "rnm.cpp"

using namespace std;
typedef vector<double> vouble;

int main()
    {
	double density = 0.3;
        srand(1);
        const double sideL = 14.0;
        const double cutoff = pow(2.0,1.0/6.0);
        double delta = 0.1;
        const int N = (int) (density*sideL*sideL);
	cout << "N = " << (int) (density*sideL*sideL) << endl;
        const int Nsweeps = 1e4;
        const int Nbins = 200;
        const int Nsample = 1e5;

        const int Nsteps = N*Nsweeps;
        int acceptence = 0;

        particle p[N];
        particle p0[N];
        vouble x_pos(N);
        vouble y_pos(N);
        x_pos = get_column("init_conf.txt",1,5);
        y_pos = get_column("init_conf.txt",2,5);

        for (int n=0;n<N;n++)
            {
                p[n].set_x(x_pos[n]);
                p[n].set_y(y_pos[n]);
                p0[n].set_x(x_pos[n]);
                p0[n].set_y(y_pos[n]);

                //cout << p[n].get_x() << "    " << p[n].get_y() << endl;
            }
        ofstream msd_out("msd_Nsweeps=" + to_string(Nsweeps) + "_delta=" + to_string(delta) + ".txt");
        

        vouble g(Nbins);
        vouble g_imd(Nbins);
        //vouble rdf(particle* p, int N, const double BoxL, const int Nbins);
        for (int k=0;k<Nsteps;k++)
            {

                if (k % (Nsteps/Nsample) == 0)
                    {
                        //cout << (double) k/Nsteps << endl;
                        msd_out << k << "    " << MSD(p0,p,N,sideL) << endl;

                        g_imd = rdf(p,N,sideL,Nbins);
                        for (int n=0;n<Nbins;n++) 
                            {
                                g_imd[n]/=Nsample;
                                g[n] += g_imd[n];
                            }
                    }
                double V_old = V_pot(p,N,sideL,cutoff);
                //cout << V_old << endl;
                int ce = rnm()*N;
                //cout << ce << endl;
                if (ce == N)
                    {
                    cout << "ce == N" << endl;
                    continue;
                    }
                double dx = (0.5*rnm()-1.0)*delta;
                double dy = (0.5*rnm()-1.0)*delta;
                
                double x_old = p[ce].get_x();
                double y_old = p[ce].get_y();

                p[ce].set_x(p[ce].get_x()+dx);
                p[ce].set_y(p[ce].get_y()+dy);

                if (p[ce].get_x() > sideL) {p[ce].set_x(p[ce].get_x()- sideL);}
                if (p[ce].get_x() < 0.0) {p[ce].set_x(p[ce].get_x()+ sideL);}

                if (p[ce].get_y() > sideL) {p[ce].set_y(p[ce].get_y()- sideL);}
                if (p[ce].get_y() < 0.0) {p[ce].set_y(p[ce].get_y()+ sideL);}

                double dV = V_pot(p,N,sideL,cutoff)-V_old;
                //cout << dV << endl;
                if (dV > 0)
                    {
                        double w0 = exp(-dV);
                        double w = rnm();
                        if (w > w0)
                            {
                                p[ce].set_x(x_old);
                                p[ce].set_y(y_old);
                                continue;
                            }
                    }
                acceptence ++;
            }
        msd_out.close();
        cout << "fertig!" << endl;
        cout << "delta = " << delta << ", acceptence rate: " << (double) acceptence/Nsteps << endl;

        ofstream final_out("finalpos_Nsweeps=" + to_string(Nsweeps) + "_delta=" + to_string(delta) + ".txt");
        for (int n=0;n<N;n++)
            {
                final_out << p[n].get_x() << "    " << p[n].get_y() << endl;
            }
        final_out.close();

    ofstream g_out("g_Nsweeps=" + to_string(Nsweeps) + "_delta=" + to_string(delta) + ".txt");
        
    for (int i=0;i<Nbins;i++)
        {
           g_out << i << "    " << (i+0.5)*0.5*sideL/Nbins << "    " << g[i] << endl;
            //cout << i << "    " << (i+0.5)*0.5*sideL/Nbins << "    " << g[i] << endl;
        }
    g_out.close(); 


        getchar();
        return 0;
    }