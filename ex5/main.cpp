#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>

#define LOG(x) cout << x << endl;
using namespace std;

// costante di Boltzmann in eV K^-1
const double k = 8.6173e-5;

// parametro LJ Argon in m
const double sigma = 3.4e-10;
// parametro LJ Argon in eV
const double epsilon = 120*k; 

// numero totale particelle
const int N = 125;
const int n = 5; 
// temperatura
const double T = 120*k/epsilon;

// totale passi Monte Carlo
const int M = 128e3;
// passi di equilibratura
const int E = 32e3;
// passi effettivi calcolo osservabili
const int m = M - E;


// potenziale di Lennard-Jones
double v (double x) {
  return 4*(1/pow(x,12) - 1/pow(x,6));
}


// derivata prima Lennard-Jones
double vprime (double x) {
  return 24/pow(x,7) - 48/pow(x,13);
}


// calcolo osservabili (potenziali)
tuple <double,double,double> Vtot(double x[], double y[], double z[], double L)
{
  double V, V2, F;
  double xij, yij, zij, rij;

  // calcolo potenziali cumulativi
  V = 0; V2 = 0; F = 0;
  for (int i=0; i<N; i++)
  {
    for (int j=0; j<i; j++)
    {
      xij = x[i] - x[j];
      yij = y[i] - y[j];
      zij = z[i] - z[j];

      // pbc su calcolo potenziale
      xij = xij - L*rint(xij/L);
      yij = yij - L*rint(yij/L);
      zij = zij - L*rint(zij/L);

      rij = sqrt(xij*xij + yij*yij + zij*zij);

      V += v(rij);
      V2 += pow(v(rij),2);

      if (rij < L/2) // cutoff distance per le forze
      {
        F += vprime(rij)*rij;
      }
    }
  }

  return {V, V2, F};
}


tuple <double,double,double,double> metropolis (double rho, double Vol)
{
  double L = cbrt(Vol); // cbrt: cubic root

  // coordinate configurazione
  double x[N], y[N], z[N];
  // proposta passo Monte Carlo
  double xmc[N], ymc[N], zmc[N];

  double V, V2, F;
  double Vmc, V2mc, Fmc;
  double p, eta;

  // corrego in fase di equilibratura
  double delta = L/100;

  // passi accettati
  double fw_steps;

  // osservabili
  double Pi[m], Cvi[m];

  // condizioni iniziali posizione, cellette in reticolo
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
      }
    }
  }
  
  // Monte Carlo loop
  for (int i=0; i < M; i++)
  {
    tie(V, V2, F) = Vtot(x, y, z, L);

    // proposta Monte Carlo
    for (int i=0; i<N; i++)
    {
      xmc[i] = x[i] + delta * ((double)rand() / RAND_MAX - 0.5);
      ymc[i] = y[i] + delta * ((double)rand() / RAND_MAX - 0.5);
      zmc[i] = z[i] + delta * ((double)rand() / RAND_MAX - 0.5);
    }

    for (int i=0; i<N; i++)
      {
        // condizioni periodiche al contorno
        xmc[i] = xmc[i] - L*round(xmc[i]/L);
        ymc[i] = ymc[i] - L*round(ymc[i]/L);
        zmc[i] = zmc[i] - L*round(zmc[i]/L);
      }

    // calcolo potenziale passo proposto
    tie(Vmc, V2mc, Fmc) = Vtot(xmc, ymc, zmc, L);

    eta = (double)rand() / RAND_MAX;
    p = exp(-(Vmc - V) / T);

    // controllo passo
    if (p > eta)
    {
      // passo accettato, aumento counter
      fw_steps++;
      for (int i=0; i<N; i++)
      {
        x[i] = xmc[i]; y[i] = ymc[i]; z[i] = zmc[i];
      }
    }

    // fase di equilibratura e correzione delta 
    if (i < E)
    {
      // ogni 500 passi controlla passi accettati su passi fatti
      if ((i > 0) and (i % 1000 == 0))
      {
        // cambio delta se passi accettati non in intervallo [35%,65%]
        if ((fw_steps > 650) or (fw_steps < 350))
        {
          delta = delta*(1 + 0.2 * ((double) fw_steps/1000 - 0.5));
        }
        // reset counter
        fw_steps = 0;
      }
    }

    // equilibrio
    if (i >= E)
    {
      // calcolo osservabili
      Pi[i-E] = N * T / Vol - F / N / (3*Vol);
      Cvi[i-E] = 1.5 * N + (V2/N - pow(V/N, 2)) / pow(T, 2);
    }
  }

  double P = 0, deltaP = 0;
  double Cv = 0, deltaCv = 0;

  // lunghezza massima calcolo autocorrelazione
  int l = 5000;

  double P0, Cv0;
  double corrP, corrCv;

  for (int i=0; i < m; i++)
  {
    // sommo e medio valori osservabili
    P += Pi[i] / m; Cv += Cvi[i] / m;

    // punto iniziale correlazione
    if (i % l == 0)
    {
      P0 = Pi[i]; Cv0 = Cvi[i];
    }

    // calcolo correlazione e medio
    corrP += P0 * Pi[i] / m;
    corrCv += Cv0 * Cvi[i] / m;
  }

  deltaP = sqrt(fabs(corrP - pow(P, 2))) / sqrt(m-1);
  deltaCv = sqrt(fabs(corrCv - pow(Cv, 2))) / sqrt(m-1);
  
  Cv = Cv/N; deltaCv = deltaCv/N;

  return {P, deltaP, Cv, deltaCv};
}

int main()
{
  ofstream argon;
  argon.open("argon.csv");
  argon << "rho,V,P,deltaP,Cv,deltaCv" << endl;

  // linearspace tra 0.001 e 1.2
  double rho[20] = {1.00000000e-03, 6.41052632e-02, 1.27210526e-01,
    1.90315789e-01, 2.53421053e-01, 3.16526316e-01, 3.79631579e-01,
    4.42736842e-01, 5.05842105e-01, 5.68947368e-01, 6.32052632e-01,
    6.95157895e-01, 7.58263158e-01, 8.21368421e-01, 8.84473684e-01,
    9.47578947e-01, 1.01068421, 1.07378947, 1.13689474, 1.20000000};

  double Vol, P, deltaP, Cv, deltaCv;

  for (int i=0; i < 20; i++)
  {
    Vol = N/rho[i];
    
    tie(P, deltaP, Cv, deltaCv) = metropolis(rho[i], Vol);

    argon << rho[i] << "," << Vol << ",";
    argon << P << "," << deltaP << ",";
    argon << Cv << "," << deltaCv << endl;
  }

  return 0;
}