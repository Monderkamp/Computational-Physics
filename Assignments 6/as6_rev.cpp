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


const double sideL = 28.0;
const double cutoff = pow(2.0,1.0/6.0);
const double V = pow(2.0,sideL);

typedef vector<double> vouble; 

vouble f_pq(particle p, particle q); // force on particle q from particle p
vouble f_i(particle* p,int i, int N); // force on particle with index i
double V_pot(particle* p, int N);
double T_kin(particle* p, int N);
double P(particle* p, int N);

int main()
  {
    vouble pos_x = get_column("final_positions40000.txt",0,4);
    vouble pos_y = get_column("final_positions40000.txt",1,4);
    vouble vel_x = get_column("final_positions40000.txt",2,4);
    vouble vel_y = get_column("final_positions40000.txt",3,4);

    const int N = pos_x.size();  
    /*
    for (int n=0;n<N;n++)
        {
            //cout << pos_x[n] << "    " << pos_y[n] << "    " << vel_x[n] << "    " << vel_y[n] << endl;
            cout << pos_y[n] << endl;
        }
    */

 
  
    /*
    double COM_vel_x = 0.0;	
    double COM_vel_y = 0.0;

    for (int n=0;n<N;n++)
        {
            COM_vel_x += vel_x[n];
            COM_vel_y += vel_y[n];
        }
    COM_vel_x/=N;
    COM_vel_y/=N;
    cout << "COM_vel_x = " << COM_vel_x << endl;
    cout << "COM_vel_y = " << COM_vel_y << endl;

    for (int n=0;n<N;n++)
        {
        vel_x[n]-= COM_vel_x;
        vel_y[n]-= COM_vel_y;
        }
    */



    cout << "N = " << N << endl;
    double dt = -0.0005;
    const int Nsteps = 4e4;
    const double tmax = Nsteps * dt;
    cout << "dt = "<< dt << endl;
    cout << "Nsteps = "<< Nsteps << endl;
    cout << "tmax = "<< tmax << endl;

    particle *p = new particle[N];
    //particle p[N];
    for (int i=0;i<N;i++)
      {
        p[i].set_x(pos_x[i]);
        p[i].set_y(pos_y[i]);
        p[i].set_vx(vel_x[i]);
        p[i].set_vy(vel_y[i]);
      }
    
    ofstream out("termodyn_rev_Nsteps=" + to_string(Nsteps) + ".txt");	
    for (int k=0;k<Nsteps;k++)
      {
        for (int i=0;i<N;i++)
          {
	    vouble f(2);
	    f = f_i(p,i,N);
            //cout << f[0] << "    " << f[1] << endl;
            p[i].set_x(p[i].get_x() + dt*p[i].get_vx() + 0.5 * f[0] * dt*dt);
            p[i].set_y(p[i].get_y() + dt*p[i].get_vy() + 0.5 * f[1] * dt*dt);
	    
    	    if (p[i].get_x() > sideL) {p[i].set_x(p[i].get_x()-sideL);}
 	    if (p[i].get_x() < 0.0) {p[i].set_x(p[i].get_x()+sideL);}
    	    if (p[i].get_y() > sideL) {p[i].set_y(p[i].get_y()-sideL);}
 	    if (p[i].get_y() < 0.0) {p[i].set_y(p[i].get_y()+sideL);}

	    //cout << p[i].get_x() << " " << p[i].get_y() << endl;


	    p[i].set_vx(p[i].get_vx()+0.5*dt*f[0]);
	    p[i].set_vy(p[i].get_vy()+0.5*dt*f[1]);

	    f = f_i(p,i,N);

	    p[i].set_vx(p[i].get_vx()+0.5*dt*f[0]);
	    p[i].set_vy(p[i].get_vy()+0.5*dt*f[1]);
          }
	out << k*dt << "  " << 2.0*T_kin(p,N)/(3.0*N) << "  " << P(p,N) << "  "<< V_pot(p,N) << "  " << T_kin(p,N) << "  " << T_kin(p,N)+V_pot(p,N) << "  " << endl;
	//cout << k*dt << "  " << 2.0*T_kin(p,N)/(3.0*N) << "  " << P(p,N) << "  "<< V_pot(p,N) << "  " << T_kin(p,N) << "  " << T_kin(p,N)+V_pot(p,N) << "  " << endl;
	
        
	if (k % (Nsteps/100) == 0)
            {
                cout << (double)k/Nsteps << endl;
            }
           
		
      }

    out.close();
    delete[] p;

    ofstream outpos("final_positions_rev" + to_string(Nsteps) + ".txt");
    for (int i=0;i<N;i++)
        {
            outpos << p[i].get_x() << "    " << p[i].get_y() << "    " << p[i].get_vx() << "    " << p[i].get_vy()<< endl;        
        }
    outpos.close();

    
    cout << "100 und fertig!" << endl;
    getchar();

    return 0;
  }

vouble f_pq(particle p, particle q) // force on particle q from particle p
  {
    vouble force(2);
    force[0] = 0.0;
    force[1] = 0.0;

    double dx = q.get_x() - p.get_x();
    double dy = q.get_y() - p.get_y();
    if (dx > 0.5*sideL) {dx -= sideL;}
    if (dx < -0.5*sideL) {dx += sideL;}
    if (dy > 0.5*sideL) {dy -= sideL;}
    if (dy < -0.5*sideL) {dy += sideL;}

    double dr = sqrt(dx*dx+dy*dy);
    if (dr <= cutoff) 
      {
        force[0] += (dx/dr)*(48.0*pow(dr,-13.0)-24.0*pow(dr,-7.0));
        force[1] += (dy/dr)*(48.0*pow(dr,-13.0)-24.0*pow(dr,-7.0));
      }

    return force;
  }

vouble f_i(particle* p,int i, int N) // force on particle with index i
  {
    vouble force(2);
    force[0] = 0.0;
    force[1] = 0.0;

    for (int j=0; j<N;j++)
      {
        if (j != i)
	  {
	    //cout << p[i].get_x() << "    " << p[i].get_y()  << "    " << p[j].get_x()<< "    " << p[j].get_y() << endl;
            double dx = p[i].get_x() - p[j].get_x();
            double dy = p[i].get_y() - p[j].get_y();          
         
	    //cout << dx << "    " << dy << endl;

    	    if (dx > 0.5*sideL) {dx -= sideL;}
    	    if (dx < -0.5*sideL) {dx += sideL;}
    	    if (dy > 0.5*sideL) {dy -= sideL;}
    	    if (dy < -0.5*sideL) {dy += sideL;}



    	    double dr = sqrt(dx*dx+dy*dy);
	    //cout << dr << endl;
    	    if (dr <= cutoff) 
      	      {
        	force[0] += (dx/dr)*(48.0*pow(dr,-13.0)-24.0*pow(dr,-7.0));
        	force[1] += (dy/dr)*(48.0*pow(dr,-13.0)-24.0*pow(dr,-7.0));
      	      }
          }
      }
    return force;
  }

double V_pot(particle *p, int N)
  {
    double U = 0.0;
        for (int j=0; j<N;j++)
          {
            for (int i=j+1;i<N;i++)
     	      {
                double dx = p[i].get_x() - p[j].get_x();
                double dy = p[i].get_y() - p[j].get_y();
    	        if (dx > 0.5*sideL) {dx -= sideL;}
    	        if (dx < -0.5*sideL) {dx += sideL;}
    	        if (dy > 0.5*sideL) {dy -= sideL;}
    	        if (dy < -0.5*sideL) {dy += sideL;}

    	        double dr = sqrt(dx*dx+dy*dy);
		//cout << dr << endl;

    	        if (dr <= cutoff) 
      	          {
                    U += 4.0*(pow(dr,-12.0)-pow(dr,-6.0));
      	          }
              }
          }
    return U;
  }

double T_kin(particle *p, int N)
  {
    double T = 0.0;
    for (int i =0;i<N;i++)
      {
        T += 0.5*(p[i].get_vx()*p[i].get_vx()+p[i].get_vy()*p[i].get_vy());
      }
    return T;   
  }

double P(particle *p, int N)
  {
    double P = 0.0;
    
    for (int i =0;i<N;i++)
      {
        for (int j=0;j<N;j++)
	  {
	    if (i != j)
	      {
	        vouble force_ij = f_pq(p[i],p[j]);
	        P += (p[j].get_x()-p[i].get_x())*force_ij[0] + (p[j].get_y()-p[i].get_y())*force_ij[1];
	      }
	  }
      }
    P /= 6.0*sideL*sideL;
    P += (2.0/2.0*V)*T_kin(p,N);
    return P;   
  }