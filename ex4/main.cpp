#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>

#define LOG(x) cout << x << endl;
using namespace std;

const double sigma = 3.94; // Angstrom
const double epsilon = 0.02; // eV
const double m = 2.2e-25; // Kg
const double k_B = 8.6173e-5; // eV K^-1

const double e = 1.6021e-19; // C, carica elettrone

const int n = 5; // numero particelle per lato a t=0
const int N = 125; // numero totale particelle
const double L = 50./sigma; // lato box di tutte le particelle

const double T = 300.*k_B / epsilon; 
const double delta_t = 1.e-4*sqrt(epsilon*e / m / pow(sigma, 2)); // 0.01 ps

const int steps = 1e6;

const double pi = 4*atan(1);


// potenziale di Lennard-Jones
double v (double x) {
  return 4*(1/pow(x,12) - 1/pow(x,6));
}


// derivata prima Lennard-Jones
double vprime (double x) {
  return 24/pow(x,7) - 48/pow(x,13);
}


// ritorna numero random su variabile gaussiana standardizzata
double gaussian_rand (double std)
{
  double u,v,x;
  u = (double)rand() / RAND_MAX;
  v = (double)rand() / RAND_MAX;

  x = sqrt(-2*pow(std,2)*log(1-u))*cos(2*pi*v);
  return x;
}


// calcola forza di j su i, ritorna componenti forza e potenziale
tuple <double,double,double,double> Fij (double xij, double yij, double zij)
{
  double Fxij, Fyij, Fzij;
  Fxij=0; Fyij=0; Fzij=0;

  // pbc su calcolo forze
  xij = xij - L*rint(xij/L);
  yij = yij - L*rint(yij/L);
  zij = zij - L*rint(zij/L);

  double rij = sqrt(xij*xij + yij*yij + zij*zij);

  if (rij < L/2) // cutoff distance
  {
    Fxij = -vprime(rij)/rij * xij;
    Fyij = -vprime(rij)/rij * yij;
    Fzij = -vprime(rij)/rij * zij;
  }
  return {Fxij, Fyij, Fzij, v(rij)};
}


int main ()
{ 
  ofstream xenon;
  xenon.open("xenon.csv");
  xenon << "x,y,z" << endl;

  ofstream xenon_energy;
  xenon_energy.open("xenon_energy.csv");
  xenon_energy << "K,V,E" << endl;

  ofstream xenon_cvv;
  xenon_cvv.open("xenon_cvv.csv");
  xenon_cvv << "cvv" << endl;

  double x[N], y[N], z[N];
  double vx[N], vy[N], vz[N];

  double V, K;

  // condizioni iniziali posizione
  // posiziono particella in angolo inferiore di ogni celletta
  int ind;
  for (int i=0; i<n; i++)
  {
    for (int j=0; j<n; j++)
    {
      for (int k=0; k<n; k++)
      {
        // conto in base 5
        ind = i*pow(n,2) + j*n + k;
        
        x[ind] = i*L/n;
        y[ind] = j*L/n;
        z[ind] = k*L/n;

        xenon << x[ind] << "," << y[ind] << "," << z[ind] << endl;
      }
    }
  }

  // energia potenziale iniziale
  double Fx[N], Fy[N], Fz[N];
  V = 0;
  for (int i=0; i<N; i++)
  {
    // inizializzo a 0 i vettori con le forze
    Fx[i]=0; Fy[i]=0; Fz[i]=0;

    for (int j=0; j<i; j++)
    {
      double xij = x[i] - x[j];
      double yij = y[i] - y[j];
      double zij = z[i] - z[j];

      // pbc su calcolo forze
      xij = xij - L*rint(xij/L);
      yij = yij - L*rint(yij/L);
      zij = zij - L*rint(zij/L);

      double rij = sqrt(xij*xij + yij*yij + zij*zij);
      V += v(rij);
    }
  }

  // condizioni iniziali velocità
  double std = sqrt(T);
  for (int i=0; i<N; i++)
  {
    vx[i] = gaussian_rand(std);
    vy[i] = gaussian_rand(std);
    vz[i] = gaussian_rand(std);

    // calcolo energia cinetica iniziale
    K += 0.5*(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  }

  // salvo valori energie iniziali
  xenon_energy << K << ","<< V << "," << K+V << endl;

  double Fxij, Fyij, Fzij, Vij;

  double v0x[N], v0y[N], v0z[N];
  double cvv0, cvv[5000]; int t;

  for (int i; i<5000; i++) {cvv[i] = 0;}

// passi di dinamica
  for (int step=1; step < steps; step++)
  {

    V=0; K=0;
    
    // velocity-verlet
    for (int i=0; i<N; i++)
    {
      // calcolo posizioni a t+dt
      x[i] = x[i] + vx[i]*delta_t + 0.5*Fx[i]*pow(delta_t,2);
      y[i] = y[i] + vy[i]*delta_t + 0.5*Fy[i]*pow(delta_t,2);
      z[i] = z[i] + vz[i]*delta_t + 0.5*Fz[i]*pow(delta_t,2);

      // condizioni periodiche al contorno
      x[i] = x[i] - L*round(x[i]/L);
      y[i] = y[i] - L*round(y[i]/L);
      z[i] = z[i] - L*round(z[i]/L);

      // sommo contributo F(t) a v
      vx[i] = vx[i] + 0.5*Fx[i]*delta_t;
      vy[i] = vy[i] + 0.5*Fy[i]*delta_t;
      vz[i] = vz[i] + 0.5*Fz[i]*delta_t;

      xenon << x[i] << "," << y[i] << "," << z[i] << endl;
    }

    // calcolo delle forze a t+dt
    for (int i=0; i<N; i++)
    {
      Fx[i]=0; Fy[i]=0; Fz[i]=0; Vij=0;

      for (int j=0; j<i; j++)
      {
        tie(Fxij, Fyij, Fzij, Vij) = Fij(x[i]-x[j], y[i]-y[j], z[i]-z[j]);
        // sommo a i la forza impressa da j
        Fx[i]+=Fxij; Fy[i]+=Fyij; Fz[i]+=Fzij;
        // sommo a j la forza impressa da i
        Fx[j]-=Fxij; Fy[j]-=Fyij; Fz[j]-=Fzij;
        V+=Vij;
      }
    }

    // sommo contributo F(t+dt) a v
    for (int i=0; i<N; i++)
    {
      vx[i] = vx[i] + 0.5*Fx[i]*delta_t;
      vy[i] = vy[i] + 0.5*Fy[i]*delta_t;
      vz[i] = vz[i] + 0.5*Fz[i]*delta_t;

      K += 0.5*(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    }

    xenon_energy << K << ","<< V << "," << K+V << endl;

    // calcolo cvv ogni 5000 passi
    t = 5000 - (steps - step) % 5000;

    // inizializzo punto iniziale cvv, solo dopo 100k step
    if (t==5000 and step > 1e5-1)
    {
      cvv0 = 0;
      // sommo contributo ogni particella
      for (int i=0; i<N; i++)
      {
        v0x[i] = vx[i];
        v0y[i] = vy[i];
        v0z[i] = vz[i];

        cvv0 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
      }
    }

    // cvv al passo t
    if (t > 0 and step > 1e5-1)
    {
      // sommo contributo ogni particella
      for (int i=0; i<N; i++)
      {
        cvv[t] += (v0x[i]*vx[i] + v0y[i]*vy[i] + v0z[i]*vz[i]) / cvv0;
      }
    }
  }

  // salvo valori cvv
  xenon_cvv << 1 << endl;
  for (int i=1; i<5000; i++)
  {
    // medio su numero cicli di calcolo di cvv
    cvv[i] = cvv[i] / (steps-1e5) * 5000;
    xenon_cvv << cvv[i] << endl;
  }
}