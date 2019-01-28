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
double muexofT(double T);
double P_acc_rem(int N,double L, double mu, double T, double U_n, double U_nm);
double P_acc_ins(int N,double L, double mu, double T, double U_np, double U_n);
int main()
    {
	double density = 0.3;
        srand(1);
        //const double sideL = 14.0;
        const double cutoff = pow(2.0,1.0/6.0);
        double delta = 0.1;
        int N = 144;
	cout << "N = " << 144 << endl;
	const double sideL =  sqrt(N/(double)density);
	cout << "sideL = " << sideL << endl;

        const int Nsweeps = 1e4;
        //const int Nbins = 200;
        //const int Nsample = 1e5;
        const double T = 3.0;
        cout << "T = " << T << endl;
        const double mu_ex = muexofT(T);
        const int Nsteps = N*Nsweeps;
        int acceptence = 0;
        int distrNmax = 300;
        vector<double> N_dis(distrNmax);
        for (int n =0; n<distrNmax;n++) N_dis[n] = 0;

        const double mu_id = -T*log(sideL*sideL/(N+1)); // k_B = \lambda = 1;

        double average_exp_dU = 0.0;
        const double mu_total = mu_id + mu_ex;
        //const double mu_total = mu_ex;

        vector<particle> p(N);

        ofstream init_out("init_Nsweeps=" + to_string(Nsweeps) + "_delta=" + to_string(delta) + ".txt");

        for (int n=0;n<N;n++)
            {
                p[n].set_x(rnm()*sideL);
                p[n].set_y(rnm()*sideL);

                init_out << p[n].get_x() << "    " << p[n].get_y() << endl;
            }
        init_out.close();
;
        for (int k=0;k<Nsteps;k++)
            {
                if (k % (Nsteps/100) == 0) cout << (double) k/Nsteps << endl;
                double insdel = rnm();
                particle trial;
                double U_old = V_pot(p,N,sideL,cutoff);
                double U_new = 0.0;                

                if (insdel <= 0.5)
                    {
                        int ce = (int) rnm()*N;
                        if (ce == N) continue;
                        trial.set_x(p[ce].get_x());
                        trial.set_y(p[ce].get_y());

                        p[ce].set_x(p[N-1].get_x());
                        p[ce].set_y(p[N-1].get_y());
                        p.pop_back();
                        
                        U_new = V_pot(p,N-1,sideL,cutoff);
                        //double P_acc_rem(int N,double L, double mu, double T, double U_n, double U_nm)
                        if (rnm() > P_acc_rem(N,sideL,mu_total,T,U_old,U_new) ) 
                            {
                                p.push_back(trial);
                                
                                N_dis[N] += 1/(double)Nsteps;
                                continue;
                            }
                        else
                            {
                                N_dis[N-1] += 1/(double)Nsteps;                            
                                N--;
                            }

                    }
              
                if (insdel > 0.5)
                    {
                        particle trial;
                        trial.set_x(rnm()*sideL);
                        trial.set_y(rnm()*sideL);
                        p.push_back(trial);
                        U_new = V_pot(p,N+1,sideL,cutoff);
                        //double P_acc_ins(int N,double L, double mu, double T, double U_np, double U_n)
                        if (rnm() > P_acc_ins(N, sideL, mu_total,T,U_new, U_old) ) 
                            {
                                p.pop_back();
                                N_dis[N] += 1/(double)Nsteps;
                                continue;
                            }
                        else
                            {
                                N_dis[N+1] += 1/(double)Nsteps;                            
                                N++;
                            }


                    }
            }

        

        cout << "fertig!" << endl;
        cout << "delta = " << delta << ", acceptence rate: " << (double) acceptence/Nsteps << endl;

        ofstream final_out("finalpos_Nsweeps=" + to_string(Nsweeps) + "_delta=" + to_string(delta) + ".txt");
        for (int n=0;n<N;n++)
            {
                final_out << p[n].get_x() << "    " << p[n].get_y() << endl;
            }
        final_out.close();

        ofstream Ndis_out("Ndis_Nsweeps=" + to_string(Nsweeps) + "_T=" + to_string(T) + ".txt");
        for (int n=0;n<distrNmax;n++)
            {
                Ndis_out << n << "    " << N_dis[n] << endl;
            }
        Ndis_out.close();
        getchar();
        return 0;
    }

double min(const double a, const double b)
    {
        if (a<b) return a;
        else return b;
    }

double P_acc_ins(int N,double L, double mu, double T, double U_np, double U_n)
    {
        return  min(1.0, (double) (L*L*exp(mu/T)*exp(-(U_np-U_n)/T))/(N+1.0));
    }

double P_acc_rem(int N,double L, double mu, double T, double U_n, double U_nm)
    {
        return  min(1.0, (double) (N*exp(-mu/T)*exp(-(U_n-U_nm)/T))/(L*L));
    }
double muexofT(double T)
    {
        return  (0.8756*T+0.2860);
    }