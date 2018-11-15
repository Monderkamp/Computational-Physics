#include <iostream>

#include <vector>

#include <cmath>
#include <fstream>
#include <string>

using namespace std;


#include "readin.cpp"
typedef vector<double> vouble; 

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

    private:
    double x,y,vx,vy;
  };

int main()
  {
    vouble pos_x = get_column("posdat2.txt",2,3);
    vouble pos_y = get_column("posdat2.txt",3,3);
    vouble vel_x = get_column("veldat2.txt",2,3);
    vouble vel_y = get_column("veldat2.txt",3,3);
    int N = pos_x.size(); 
    particle p[N];
    for (int i=0;i<N;i++)
      {
        p[i].set_x(pos_x[i]);
        p[i].set_y(pos_y[i]);
        p[i].set_vx(vel_x[i]);
        p[i].set_vy(vel_y[i]);
      }


    getchar();
    return 0;
  }


