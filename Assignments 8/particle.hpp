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
        z = 0.0;
        vx = 0.0;
        vy = 0.0;
        vz = 0.0;
      }

    particle(const particle& p)
      {
        x = p.x;
        y = p.y;
        z = p.z;
        vx = p.vx;
        vy = p.vy;
        vz = p.vz;
      }
    ~particle(){}
    
    void set_x(const double a){x = a;}
    void set_y(const double a){y = a;}
    void set_z(const double a){z = a;}
    void set_vx(const double a){vx = a;}
    void set_vy(const double a){vy = a;}
    void set_vz(const double a){vz = a;}
    
    double get_x(){return x;}
    double get_y(){return y;}
    double get_z(){return z;}
    double get_vx(){return vx;}
    double get_vy(){return vy;}
    double get_vz(){return vz;}

    private:
    double x,y,z,vx,vy,vz;
  };

#endif