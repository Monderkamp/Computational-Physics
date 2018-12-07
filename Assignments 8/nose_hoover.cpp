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





typedef vector<double> vouble; 

int main()
  {
    const double sideL = 14.0;
    const double V = pow(2.0,sideL);
    const double cutoff = pow(2.0,1.0/6.0);

    vouble pos_x = get_column("init_conf.txt",1,5);
    vouble pos_y = get_column("init_conf.txt",2,5);
    vouble vel_x = get_column("init_conf.txt",3,5);
    vouble vel_y = get_column("init_conf.txt",4,5);

    const int N = pos_x.size();  

    cout << "N = " << N << endl;


    const double tmax = 10.0;
    const int Nsteps = 2e4;
    double dt = (double)tmax/Nsteps;

    const int Nbins = 500;
    const int Nsample = 1e3;
    vouble g(Nsample);
    
    double Qprime = 1.0;
    double T = 1.0;
    
    const double Q = Qprime * 2.0*N;
    cout << "dt = "<< dt << endl;
    cout << "Nsteps = "<< Nsteps << endl;
    cout << "tmax = "<< tmax << endl;

    //particle *p = new particle[N];
    particle p[N];

    for (int i=0;i<N;i++)
      {
        p[i].set_x(pos_x[i]);
        p[i].set_y(pos_y[i]);
        p[i].set_vx(vel_x[i]);
        p[i].set_vy(vel_y[i]);
      }

    double alpha = sqrt(144/T_kin(p,N));

    for (int i=0;i<N;i++)
	{
            p[i].set_vx(p[i].get_vx()*alpha);
            p[i].set_vy(p[i].get_vy()*alpha);
	}
    double eta = 0.0;
    double p_eta = 0.0;
    vouble p_schlange0(2);

    ofstream out("thst_termodyn_Nsteps=" + to_string(Nsteps) + "_dt=" + to_string(dt) + ".txt");	
    for (int k=0;k<Nsteps;k++)
      {
	//if (k*dt >= 5.0) T = 1.0;
	
        if (k>= 0.1*Nsteps)
            {
                if ((int)(k - 0.1*Nsteps) % (int)(0.9*Nsteps/Nsample) == 0)
                    {
                        vouble g_imd(Nsample);
                        g_imd = rdf(p,N,sideL,Nbins);

                        for (int k=0;k<Nbins;k++)
                            {
                                g_imd[k] /= Nsample;
                                g[k] += g_imd[k];
                            }
                    }
            }        


	double p_eta_dot_0 = p_eta_dot(p,N,T);
	double p_eta_imd = p_eta + 0.5*dt*p_eta_dot_0;
	double exp_fac = exp(-0.5*dt*(1/Q)*p_eta_imd);

        for (int i=0;i<N;i++)
          {
	    vouble f(2);
	    f = f_i(p,i,N,sideL,cutoff);

            p_schlange0[0] = p[i].get_vx()*exp_fac;
            p_schlange0[1] = p[i].get_vy()*exp_fac;

	    p[i].set_x(p[i].get_x() + dt* p_schlange0[0] + 0.5*dt*dt*f[0]);
	    p[i].set_y(p[i].get_y() + dt* p_schlange0[1] + 0.5*dt*dt*f[1]);
            
            p[i].set_vx(p[i].get_vx()*exp_fac + 0.5*dt*exp_fac*f[0]);
            p[i].set_vy(p[i].get_vy()*exp_fac + 0.5*dt*exp_fac*f[1]);

	    f = f_i(p,i,N,sideL,cutoff);

            p[i].set_vx(p[i].get_vx() + 0.5*dt*exp_fac*f[0]);
            p[i].set_vy(p[i].get_vy() + 0.5*dt*exp_fac*f[1]);

	    eta += (dt/Q)*p_eta +  0.5*dt*dt*(1/Q)*p_eta_dot_0;	 
            p_eta += 0.5*dt*(p_eta_dot_0 + p_eta_dot(p,N,T));

    	    if (p[i].get_x() > sideL) {p[i].set_x(p[i].get_x()-sideL);}
 	    if (p[i].get_x() < 0.0) {p[i].set_x(p[i].get_x()+sideL);}
    	    if (p[i].get_y() > sideL) {p[i].set_y(p[i].get_y()-sideL);}
 	    if (p[i].get_y() < 0.0) {p[i].set_y(p[i].get_y()+sideL);}
          }


	out << k*dt << "  " << 2.0*T_kin(p,N)/(3.0*N) << "  " << P(p,N,sideL,cutoff) << "  "<< V_pot(p,N,sideL,cutoff) << "  " << T_kin(p,N) << "  " << T_kin(p,N)+V_pot(p,N,sideL,cutoff) << "  " << endl;
        
	if (k % (Nsteps/100) == 0)
            {
                cout << (double)k/Nsteps << endl;
            }
           
		
      }

    out.close();
    //delete[] p;

    ofstream outpos("final_positions_Nsteps=" + to_string(Nsteps) + "_dt=" + to_string(dt) + ".txt");
    for (int i=0;i<N;i++)
        {
            outpos << p[i].get_x() << "    " << p[i].get_y() << "    " << p[i].get_vx() << "    " << p[i].get_vy()<< endl;        
        }
    outpos.close();

    ofstream g_out("g_nose_hoover.txt");
        
    for (int i=0;i<Nbins;i++)
        {
           g_out << i << "    " << (i+0.5)*0.5*sideL/Nbins << "    " << g[i] << endl;
            //cout << i << "    " << (i+0.5)*0.5*sideL/Nbins << "    " << g[i] << endl;
        }
    g_out.close();    
    cout << "100 und fertig!" << endl;
    //getchar();

    return 0;
  }

