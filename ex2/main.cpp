#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>

#define LOG(x) cout << x << endl;
using namespace std;

const double treshold = 1e-6;

double radius(double M, double xi_bar, double phi_xi_bar)
{
    double h_bar = 1.0545718e-34;  //joule * sec
    double m = 1.674e-27; //Kg
    double k = ( pow(h_bar,2) * pow(3* pow(M_PI,2),2.0/3.0) ) / ( 5 * pow(m,8.0/3.0) );
    double hbarc = 3.16e-26; //m^3 * Kg * s^-2 
    double G = 6.67e-11; //m^3 * Kg^-1 * s^-2
    return pow( - 4 * M_PI * pow((5*k)/(8*M_PI*G),3.0) * (pow(xi_bar,5) / M) * phi_xi_bar , 1.0/3.0);
    
}
// RK4 per ODE di 2 grado 
// d(x(t))/dt = f(t,x(t),y(t)) con f = y(t)
// d(y(t))/dt = g(t,x(t),y(t)) con g = -2/t *y(t) - x(t)^n
tuple <double,double> RK4 (double ti, double xi, double yi, double h, double n)
{
  double k0,k1,k2,k3;
  double l0,l1,l2,l3;
  double xi_1,yi_1; //x(ti+1), y(ti+1)

  // calcolo predictor-corrector ordine 0
  k0 = h * yi;
  l0 = -h * (2/ti * yi + pow(xi,n));

  // cout << "k0 = "; LOG(k0);
  // cout << "l0 = "; LOG(l0);

  // calcolo predictor-corrector ordine 1
  k1 = h * (yi+l0/2);
  l1 = -h * (2/(ti+h/2) * (yi+l0/2) + pow(xi+k0/2,n)); 

  // cout << "k1 = "; LOG(k1);
  // cout << "l1 = "; LOG(l1);

  // calcolo predictor-corrector ordine 2
  k2 = h * (yi+l1/2);
  l2 = -h * (2/(ti+h/2) * (yi+l1/2) + pow(xi+k1/2,n));

  // cout << "k2 = "; LOG(k2);
  // cout << "l2 = "; LOG(l2);

  // calcolo predictor-corrector ordine 3
  k3 = h * (yi+l2);
  l3 = -h * (2/(ti+h) * (yi+l2) + pow(xi+k2,n));

  // cout << "k3 = "; LOG(k3);
  // cout << "l3 = "; LOG(l3);

  // xi_1 = x(t_i+1), yi_1 = y(t_i+1)
  xi_1 = xi + 1/6.0 * (k0 + 2*(k1+k2) + k3);
  yi_1 = yi + 1/6.0 * (l0 + 2*(l1+l2) + l3);

  return {xi_1, yi_1};
}

int main()
{
  double h = 1e-3;
  double xi_bar, phi_xi_bar;

  // n in [1.5,3]
  for (double n=1.5; n<=3; n=n+0.25)
  {
    ofstream file ("ns/" + to_string(n) + ".csv");
    file << "xi,theta,phi" << endl;
    file << "0,1,0" << endl;

    double xi_i = 1e-4;
    double theta_i = 1.0;
    double phi_i = 0;

    bool cond = (theta_i < treshold);

    for (int i=0; theta_i > treshold; ++i)
    {
      tie(theta_i, phi_i) = RK4(xi_i+i*h, theta_i, phi_i, h, n);

      if (theta_i > 0) 
      {
        file << to_string(xi_i+i*h) + ",";
        file << to_string(theta_i) + ",";
        file << to_string(phi_i) << endl;

        if (n==1.5) {
          xi_bar = xi_i+i*h;
          phi_xi_bar = phi_i;
        }
      }
    }
  }

  // n = {1,4.5}
  for (double n=1; n<=4.5; n=n+3.5)
  {
    ofstream file ("ns/" + to_string(n) + ".csv");
    file << "xi,theta,phi" << endl;
    file << "0,1,0" << endl;

    double xi_i = 1e-4;
    double theta_i = 1.0;
    double phi_i = 0;

    for (int i=0; xi_i+i*h < 30; ++i)
    {
      tie(theta_i, phi_i) = RK4(xi_i+i*h, theta_i, phi_i, h, n);

      file << to_string(xi_i+i*h) + ",";
      file << to_string(theta_i) + ",";
      file << to_string(phi_i) << endl;
    }
  }

  int size = 1e4;
  double m_sun = 1.989e30; //Kg
  double mass[size];
  mass[0] = 1.5 * m_sun;
  double raggio[size];

  double incremento_massa = (3*m_sun - 1.5 *m_sun)/size;

  ofstream massa;
  massa.open("dati_massa_neutron.csv"); //creo il csv con i vettori della massa e raggio da plottare
  massa << "massa,raggio" << endl;
  for(int j = 0; j<size; j++)
  {
    mass[j] = mass[0] + incremento_massa * j;
    raggio[j] = radius(mass[j] , xi_bar, phi_xi_bar);
    massa << mass[j] << "," << raggio[j] << endl;
  }

  return 0;
}