#include <iostream>
#include <cmath>
#include <fstream> 

#include "quadrature.cpp"
#include "boolsrule.cpp"

using namespace std;

using namespace std;
int main()
{
	
	
cout << "quadrature: " <<quadrature() << endl;
cout << "boolsrule: " <<boolsrule(0.0,2.0) << endl;	
cout << "\ndeviation from exact value:" << endl;
cout << "quadrature: " <<0.995322265018953-quadrature() << endl;
cout << "boolsrule: " <<0.995322265018953-boolsrule(0.0,2.0) << endl;	
cout << "\nrelative error:" << endl;
cout << "quadrature: " <<(0.995322265018953-quadrature())/995322265018953 << endl;
cout << "boolsrule: " <<(0.995322265018953-boolsrule(0.0,2.0))/995322265018953 << endl;	

ofstream out("23d.txt");

out << "quadrature: " <<quadrature() << endl;
out << "boolsrule: " <<boolsrule(0.0,2.0) << endl;	
out << "\ndeviation from exact value:" << endl;
out << "quadrature: " <<0.995322265018953-quadrature() << endl;
out << "boolsrule: " <<0.995322265018953-boolsrule(0.0,2.0) << endl;	
out << "\nrelative error:" << endl;
out << "quadrature: " <<(0.995322265018953-quadrature())/995322265018953 << endl;
out << "boolsrule: " <<(0.995322265018953-boolsrule(0.0,2.0))/995322265018953 << endl;	

out.close();

return 0;	
}
