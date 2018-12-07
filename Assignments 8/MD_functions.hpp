#include <iostream>

#include <vector>

#include <cmath>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>

#include "particle.hpp"

typedef std::vector<double> vouble;

vouble f_pq(particle p, particle q,const double sideL,const double cutoff); // force on particle q from particle p
vouble f_i(particle* p,int i, int N,const double sideL,const double cutoff); // force on particle with index i
vouble Fges(particle* p, int N,const double sideL,const double cutoff);
vouble rdf(particle* p, int N, const double BoxL, const int Nbins);

double V_pot(particle *p, int N,const double sideL,const double cutoff);
double T_kin(particle* p, int N);
double P(particle *p, int N,const double sideL,const double cutoff);
double p_eta_dot(particle* p, int N, double T);
