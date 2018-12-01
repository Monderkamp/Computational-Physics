#include <iostream>

#include <vector>

#include <cmath>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>

#include "particle.hpp"

typedef std::vector<double> vouble;

vouble f_pq(particle p, particle q); // force on particle q from particle p
vouble f_i(particle* p,int i, int N); // force on particle with index i
double V_pot(particle* p, int N);
double T_kin(particle* p, int N);
double P(particle* p, int N);