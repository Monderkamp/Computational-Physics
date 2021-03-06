#include <iostream>

#include <vector>

#include <cmath>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>

#include "MD_functions.hpp"
#include "particle.hpp"

using namespace std;


//const double sideL = 14.0;
//const double cutoff = pow(2.0,1.0/6.0);
//const double V = pow(2.0,sideL);


typedef vector<double> vouble; 

vector<double> f_pq(particle p, particle q,const double sideL,const double cutoff) // force on particle q from particle p
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

vouble f_i(particle* p,int i, int N,const double sideL,const double cutoff) // force on particle with index i
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

double V_pot(particle *p, int N,const double sideL,const double cutoff)
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
                    U += 4.0*(pow(dr,-12.0)-pow(dr,-6.0)) + 1.0;
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

double P(particle *p, int N,const double sideL,const double cutoff)
  {    
    double V = sideL*sideL;
    double P = 0.0;
    
    for (int i =0;i<N;i++)
      {
        for (int j=0;j<N;j++)
	  {
	    if (i != j)
	      {
	        vouble force_ij = f_pq(p[i],p[j],sideL,cutoff);
	        P += (p[j].get_x()-p[i].get_x())*force_ij[0] + (p[j].get_y()-p[i].get_y())*force_ij[1];
	      }
	  }
      }
    P /= 6.0*sideL*sideL;
    P += (2.0/2.0*V)*T_kin(p,N);
    return P;   
  }
vouble Fges(particle* p, int N,const double sideL,const double cutoff)
    {
	vouble Fges(2);
	Fges[0] = 0.0;
	Fges[1] = 0.0;

	for (int i=0;i<N;i++)
		{
			vouble kraft(2);
			kraft = f_i(p,i,N,sideL,cutoff);
			
			Fges[0] += kraft[0];
			Fges[1] += kraft[1];
		}
	//cout << Fges[0] << "    " << Fges[1] << endl;
        return Fges;
    }

double p_eta_dot(particle* p, int N, double T)
    {
        return (2*T_kin(p,N) - 2*N*T);
    }

vouble rdf(particle* p, int N, const double BoxL, const int Nbins)
    {
        vouble g(Nbins);
        for (int i=0;i<Nbins; i++)
            {
                g[i] = 0.0;
            }

        double BinL = (double)0.5*BoxL/Nbins;

        for (int i=0;i<N;i++)
            {
                for (int j=0;j<N;j++)
                    {
                        if (j==i) continue;
                        
                        double dx = p[j].get_x()-p[i].get_x();
                        if (dx > 0.5*BoxL) {dx -= BoxL;}
                        if (dx < -0.5*BoxL) {dx += BoxL;}

                        double dy = p[j].get_y()-p[i].get_y();
                        if (dy > 0.5*BoxL) {dy -= BoxL;}
                        if (dy < -0.5*BoxL) {dy += BoxL;}


                        double dr = sqrt(dx*dx+dy*dy); 
                        //cout << dr << endl;
                        if (dr > 0.5*BoxL) continue;

                        for (int k=0; k< Nbins; k++)
                            {
                                if (fabs(dr - (k+0.5)*BinL)< 0.5*BinL)
                                    {
                                        g[k] += 1.0;
                                    }
                            }
                    }
            }
            double rho = (double) N/(BoxL*BoxL);
        for (int i=0;i<Nbins;i++)
            {
                //g[i] = g[i]/((N*N)*(2.0*M_PI)*((double)((i+0.5)*BinL))/(BoxL*BoxL));
                //g[i] /= ((N*N)*(2.0*M_PI)*((double)((i+0.5)*BinL))/(BoxL*BoxL));
                //double r_i = (double) ((i+0.5)*BinL)
                double BinArea_i = M_PI*BinL*BinL*(2.0*i+1.0);
                g[i] /= N-1;
                g[i] /= (BinArea_i*rho);
                
            }

        return g;
    }

double MSD(particle* p0, particle* p, int N, const double BoxL)
    {
        double msd = 0.0;
        for (int i=0;i<N;i++)
            {
                double dx = p[i].get_x()-p0[i].get_x();
                double dy = p[i].get_y()-p0[i].get_y();
                
                if (dx > 0.5*BoxL) dx -= BoxL;
                if (dx < -0.5*BoxL) dx += BoxL;

                if (dy > 0.5*BoxL) dy -= BoxL;
                if (dy < -0.5*BoxL) dy += BoxL;

                double dr2 = dx*dx+dy*dy;
                msd += dr2/N;
            }       
        return msd;
    }











