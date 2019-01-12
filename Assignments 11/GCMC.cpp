#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>

#include "rnm.cpp"

using namespace std;

typedef vector <double> vouble;

double min(const double a, const double b);
double P_acc_ins(int N,double L, double mu); //insertion acceptence probability starting with N
double P_acc_rem(int N,double L, double mu); //removal acceptence probability starting with N
int overlap(vouble conf, double x,double sigma);
vouble density_profile(vouble conf,double L,double res, int Nbars);

int main()
    {
        srand(414);
        const int Nsteps = 1e8;
        cout << "total number of steps: " << Nsteps << endl;
        const int Nsample = 1e5;
        cout << "total number of samplings: " << Nsample << endl;
        const double L = 100.0;
        const double sigma = 1.0;
        const double mu = 1.0; //chemical potential
        cout << "chemical potential: mu = " << mu << endl;
        const double res = 0.1*sigma;

        
        vouble conf;
        double x_init = 0.5*sigma+rnm()*(L-sigma);
        //cout << x_init << endl;
        conf.push_back(x_init);
        //cout << conf.size() << endl;
        
        const int Nbars = (int)L/res;
        cout << "number of bins in density = " << Nbars << endl;
        vouble dens = vouble(Nbars);        
        vouble dens_imd = vouble(Nbars);

        ofstream N_out("N_evolution_mu=" + to_string(mu) + ".txt");
        
        for (int k=0;k<Nsteps;k++)
            {
            
                if (k % (Nsteps/Nsample) == 0)
                    {
                        //cout << "progress: " << (double)k/Nsteps << endl;
                        dens_imd = density_profile(conf,L,res,Nbars);
                        //cout << "dens_imd= " << endl;
                        for (int j=0;j<Nbars;j++)
                            {
                                //cout << dens_imd[j] << endl;

                                dens_imd[j]/=Nsample;
                                dens[j] += dens_imd[j];
                            }
                        N_out << k << "    " << conf.size() << endl;
                    }
                if (k % (Nsteps/100) == 0)
                    {
                        cout << "progress: " << (double)k/Nsteps << endl;
                    }

                int N = conf.size();
                //cout << N << endl;
                double u = rnm();
                if (u < 0.5) //removal

                    {


                        int ran_el = (int) rnm()*(N);
                        if (ran_el == N) continue;
                        double w0 = P_acc_rem(N,L,mu);  
                        double w = rnm();
                        if (w > w0) continue;

                        conf[ran_el] = conf[N-1];
                        conf.pop_back();

                        
                    }
                else        //insertion
                    {
                        double x_ins = 0.5*sigma+rnm()*(L-sigma);
                        //cout << x_ins << endl;
                        if (overlap(conf,x_ins,sigma) != 0) continue;
                        double w0 = P_acc_ins(N,L,mu);  
                        double w = rnm();
                        if (w > w0) continue;
                        conf.push_back(x_ins);
                    }
            }

        ofstream dens_out("dens_mu=" + to_string(mu) + ".txt");
        for (int i=0;i<Nbars;i++)
            {
                dens_out << (0.5+i)*res << "    " << dens[i] << endl;   
            }
        dens_out.close();
        N_out.close();
        double A = 0.0;
        for (int n=0;n<conf.size();n++)
            {
                cout << conf[n] << endl;
                A += conf[n];
            }
        A/=conf.size();
        cout << "<x> = " << A << endl;


        cout << "fertig!" << endl;
        getchar();
        return 0;
    }

double min(const double a, const double b)
    {
        if (a<b) return a;
        else return b;
    }

double P_acc_ins(int N,double L, double mu)
    {
        return  min(1.0, (double) (L*exp(mu))/(N+1.0));
    }

double P_acc_rem(int N,double L, double mu)
    {
        return  min(1.0, (double) (N*exp(-mu))/(L));
    }
int overlap(vouble conf, double x,double sigma)
    {
        int overl = 0;
        int N = conf.size();
        for (int n = 0;n<N;n++)
            {
                if (abs(conf[n]-x) < sigma) overl++;
            }
        if (overl == 0) return 0;
        else return 1;
    }
vouble density_profile(vouble conf,double L,double res, int Nbars)
    {
        int N = conf.size();
        //cout << "dens N = " << N << endl;
        //int Nbars = (int)L/res; 
        //cout << "dens L = " << L << ", res = " << res << endl;
        //cout << "L/res = " << (int)L/res << endl;
        //cout << "dens Nbars = " << Nbars << endl;
        vouble density(Nbars);

        for (int n=0; n<Nbars;n++) {density[n] =0.0;}

        for (int n = 0;n<N;n++)
            {
                for (int k=0;k<Nbars;k++)
                    {
                        if (abs(conf[n] - (0.5+k)*res)< 0.5*res ) {density[k]+= (double)1.0/N;}
                    }
            }
        return density;
    }


