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

const double sideL = 28.0;
const double cutoff = pow(2.0,1.0/6.0);
const double V = pow(2.0,sideL);

typedef vector<double> vouble; 



int main()
  {
    vouble pos_x = get_column("init_conf.txt",1,5);
    vouble pos_y = get_column("init_conf.txt",2,5);
    vouble vel_x = get_column("init_conf.txt",3,5);
    vouble vel_y = get_column("init_conf.txt",4,5);

    const int N = pos_x.size();  

    cout << "N = " << N << endl;
    double dt = 0.0005;
    const int Nsteps = 1e3;
    const double tmax = Nsteps * dt;
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

    for (int i;i<N;i++)
	{
            p[i].set_vx(p[i].get_vx()*alpha);
            p[i].set_vy(p[i].get_vy()*alpha);
	}


    ofstream out("termodyn_Nsteps=" + to_string(Nsteps) + "_dt=" + to_string(dt) + ".txt");	
    for (int k=0;k<Nsteps;k++)
      {
        for (int i=0;i<N;i++)
          {
	    vouble f(2);
	    f = f_i(p,i,N);
            //cout << f[0] << "    " << f[1] << endl;
            p[i].set_x(p[i].get_x() + dt*p[i].get_vx() + 0.5 * f[0] * dt*dt);
            p[i].set_y(p[i].get_y() + dt*p[i].get_vy() + 0.5 * f[1] * dt*dt);

	    p[i].set_vx(p[i].get_vx()+0.5*dt*f[0]);
	    p[i].set_vy(p[i].get_vy()+0.5*dt*f[1]);

	    f = f_i(p,i,N);

	    p[i].set_vx(p[i].get_vx()+0.5*dt*f[0]);
	    p[i].set_vy(p[i].get_vy()+0.5*dt*f[1]);

    	    if (p[i].get_x() > sideL) {p[i].set_x(p[i].get_x()-sideL);}
 	    if (p[i].get_x() < 0.0) {p[i].set_x(p[i].get_x()+sideL);}
    	    if (p[i].get_y() > sideL) {p[i].set_y(p[i].get_y()-sideL);}
 	    if (p[i].get_y() < 0.0) {p[i].set_y(p[i].get_y()+sideL);}
          }


	out << k*dt << "  " << 2.0*T_kin(p,N)/(3.0*N) << "  " << P(p,N) << "  "<< V_pot(p,N) << "  " << T_kin(p,N) << "  " << T_kin(p,N)+V_pot(p,N) << "  " << endl;
        
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

    
    cout << "100 und fertig!" << endl;
    getchar();

    return 0;
  }

