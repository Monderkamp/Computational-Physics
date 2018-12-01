#include <iostream>

#include <vector>

#include <cmath>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>

typedef std::vector<double> vouble; 

#ifndef PARTICLE_HPP
#define PARTICLE_HPP

class particle 
  {
    public:
    particle()
      {
        x = 0.0;
        y = 0.0;
        vx = 0.0;
        vy = 0.0;
      }
    ~particle(){}
    
    void set_x(const double a){x = a;}
    void set_y(const double a){y = a;}
    void set_vx(const double a){vx = a;}
    void set_vy(const double a){vy = a;}
    
    double get_x(){return x;}
    double get_y(){return y;}
    double get_vx(){return vx;}
    double get_vy(){return vy;}
    //friend vouble f_i(particle *p,int i, int N); // force on particle with index i

    private:
    double x,y,vx,vy;
  };

#endif