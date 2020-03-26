#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <tuple>
#define LOG(x) cout << x << endl;
using namespace std;

const double pi = 4*atan(1);
const double treshold = 1e-7; // soglia arbitraria approssima zero
const double gamma_ = 150; // parametro del problema 


double Xin(double E0)
{
    double a = -4;  
    double b = 4;
    double c = E0;
   
    double uin = (-b - pow(pow(b,2)-4*a*c, 0.5))/(2*a);
    double xin = pow(uin, -1.0/6);
    return xin;
}

double Xout(double E0)
{
    double a = -4;  
    double b = 4;
    double c = E0;

    double uout = (-b + pow(pow(b,2)-4*a*c, 0.5))/(2*a);
    double xout = pow(uout, -1.0/6);
    return xout;
}


// funzione del potenziale
double V (double x) {
  return 4*(1/pow(x,12) - 1/pow(x,6));
}


// derivata di V
double Vprime (double x) {
  return 24*(1/pow(x,7) - 2/(pow(x,13)));
}


// metodo di newton per soluzione di f(x) = c (default c=0)
// args: funzione, derivata, punto iniziale; return : soluzione
double newton (double(*f)(double), double(*fprime)(double), double x0, double c=0) 
{
  double xi = x0 - (f(x0)-c)/fprime(x0); // calcolo x1

  // newton's algorithm
  while (abs(f(xi) - c) > treshold) 
  {
    xi = xi - (f(xi)-c)/fprime(xi);
  }
  return xi;
}


// calcola integrale con simpson cubica su mesh valori
double simpson3 (double fi[], int N, double h) 
{
  double I = 0;

  // aggiungo componenti dispari
  for (int i=1; i<=(N/2); i++)
    {
      I = I + h/3 * 4*fi[2*i-1];
    }

  // aggiungo componenti pari
  for (int i=1; i<=(N/2-1); i++)
    {
      I = I + h/3 * 2*fi[2*i];
    }

  // aggiungo componente iniziale e finale
  I = I + (h/3)*(fi[0] + fi[N-1]);

  return I;
}

// funzione che esegue un passaggio del metodo delle secanti
// args : x0, f(x0), x1, f(x1); return: x2, x1
tuple <double,double> secant (double x0, double f0, double x1, double f1) 
{
  return {x1 - f1*(x1-x0)/(f1-f0), x1};
}


// calcola integrale d'azione usando simpon cubico
double azione (double E, double xin, double xout)
{
  int N = 1000;
  double h = (xout - xin)/N;

  double xi[N]; // creo mesh intervallo d'integrazione
  for (int i=0; i<N+1; i++)
  {
    xi[i] = xin + i*h;
  }

  double fi[N]; // creo mesh valori funzione integranda
  for (int i=0; i<N+1; i++)
  {
    fi[i] = gamma_ * sqrt( abs( E - V(xi[i]) ) );
  }

  return simpson3(fi, N, h);
}


// F = azione - pi*(n+1/2)
double F (double E, int n)
{
  double xin, xout;
  xin = Xin(E); xout = Xout(E);

  return azione(E, xin, xout) - pi*(n+0.5);
}


// cerca livello energetico vicino a E
tuple <double,double> search (double E, int n)
{
  double E0;
  double FE, FE0; // valori di F(E) e F(E0)

  E0 = E - 1e-5; // valore inziale di E0 metodo secanti

  FE0 = F(E0,n); // valore inziale di F(E0) metodo secanti
  FE = F(E,n); // valore inziale di F(E) metodo secanti

  do {
    // uso tie per leggere la tupla che ritorna secant
    tie(E, E0) = secant(E0, FE0, E, FE);

    FE0 = FE;
    FE = F(E,n);

  } while (abs(E-E0) > treshold);

  return {E, FE};
}


int main ()
{
  ofstream file ("gamma" + to_string(gamma_) + ".csv");
  file << "E,res" << endl;

  double E, Fvalue;
  E = -1 + 1e-5; // valore iniziale E
  for (int n=0; E<0; ++n)
  {
    // uso tie per leggere la tupla che ritorna search
    tie(E, Fvalue) = search(E, n);

    cout << "E" << n << " = " << E; // print livello energetico
    cout << " & F" << n << " = " << Fvalue << endl; // print valore F
    file << E << "," << Fvalue << endl;
  }
}
