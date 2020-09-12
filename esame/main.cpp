#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>

#define LOG(x) cout << x << endl;
using namespace std;

const double PI = 4*atan(1);

// numero particelle
const int N = 512;
// siti per lato del lattice 3-dimensionale
const int n = cbrt(N);

// totale passi Monte Carlo
const int M = 256e3;
// passi di equilibratura
const int E = 32e3;
// passi effettivi calcolo osservabili
const int m = M - E;


// implementazione operatore modulo
// diverso da operatore remainder % per numeri negativi
int mod (int i, int n)
{
  return ((i % n) + n) % n;
}


// calcolo H della configurazione s
tuple<double,double> computeH(double s[n][n][n], double smc[n][n][n], double J, double h)
{
  double H=0, Hmc=0;

  // divido J per due perchè ji = ij. Conto due volte interazione ij 
  J = J / 2;

  for (int j=0; j < N; j++)
  {
    // interazione destra e sinistra
    H -= J * cos(s[j/n/n%n][j/n%n][j%n] - s[j/n/n%n][j/n%n][mod(j%n-1, n)]);
    H -= J * cos(s[j/n/n%n][j/n%n][j%n] - s[j/n/n%n][j/n%n][mod(j%n+1, n)]);
    // interazione sopra e sotto
    H -= J * cos(s[j/n/n%n][j/n%n][j%n] - s[j/n/n%n][mod(j/n%n-1, n)][j%n]);
    H -= J * cos(s[j/n/n%n][j/n%n][j%n] - s[j/n/n%n][mod(j/n%n+1, n)][j%n]);
    // interazione dentro e fuori
    H -= J * cos(s[j/n/n%n][j/n%n][j%n] - s[mod(j/n/n%n-1, n)][j/n%n][j%n]);
    H -= J * cos(s[j/n/n%n][j/n%n][j%n] - s[mod(j/n/n%n+1, n)][j/n%n][j%n]);

    Hmc -= J * cos(smc[j/n/n%n][j/n%n][j%n] - smc[j/n/n%n][j/n%n][mod(j%n-1, n)]);
    Hmc -= J * cos(smc[j/n/n%n][j/n%n][j%n] - smc[j/n/n%n][j/n%n][mod(j%n+1, n)]);

    Hmc -= J * cos(smc[j/n/n%n][j/n%n][j%n] - smc[j/n/n%n][mod(j/n%n-1, n)][j%n]);
    Hmc -= J * cos(smc[j/n/n%n][j/n%n][j%n] - smc[j/n/n%n][mod(j/n%n+1, n)][j%n]);

    Hmc -= J * cos(smc[j/n/n%n][j/n%n][j%n] - smc[mod(j/n/n%n-1, n)][j/n%n][j%n]);
    Hmc -= J * cos(smc[j/n/n%n][j/n%n][j%n] - smc[mod(j/n/n%n+1, n)][j/n%n][j%n]);

    // interazione campo magnetico esterno
    H -= h * cos(s[j/n/n%n][j/n%n][j%n]); 
    Hmc -= h * cos(smc[j/n/n%n][j/n%n][j%n]); 
  }

  return {H,Hmc};
}


tuple<double,double,double,double,double,double> Metropolis (double T, double h)
{
  double J = 1;

  // vettore configurazione
  // assegno a ogni spin un angolo
  double s[n][n][n];
  // vettore proposta Monte Carlo 
  double smc[n][n][n];

  // grandezza tipico proposta passo mc
  // ottimizzo in fase equilibratura
  double delta = 3. / 180. * PI; // 3°

  // Hamiltoniana, Hamiltoniana proposta mc
  double H, Hmc;

  // variabili decisione passo mc
  double p, eta;
  // passi accettati
  double fw_steps;

  // osservabili
  double Ui[m], Mui[m][2];

  // init configurazione random
  for (int i=0; i<N; i++)
  {
    // numero random tra 0 e 2*pi
    s[i/n/n%n][i/n%n][i%n] = 2*PI * ((double)rand() / RAND_MAX);
  }

  // Monte Carlo loop
  for (int i=0; i < M; i++)
  {
    for (int j=0; j < N; j++)
    {
      // proposta Monte Carlo
      smc[j/n/n%n][j/n%n][j%n] = s[j/n/n%n][j/n%n][j%n];
      smc[j/n/n%n][j/n%n][j%n] += delta*((double) rand()/RAND_MAX - 0.5);
    }

    tie(H, Hmc) = computeH(s, smc, J, h);

    // estraggo numero random tra 0 e 1
    eta = (double)rand() / RAND_MAX;
    // calcolo "probabilità" transizione
    p = exp(-(Hmc - H) / T); 

    // controllo passo
    if (p > eta)
    {
      // passo accettato, aumento counter
      fw_steps++;
      // aggiorno configurazione
      for (int j=0; j < N; j++) 
      {
        s[j/n/n%n][j/n%n][j%n] = smc[j/n/n%n][j/n%n][j%n]; 
      }
      // aggiorno valore energia
      H = Hmc;
    }

    // fase di equilibratura e correzione delta 
    if (i < E)
    {
      if (fw_steps > 1000) {fw_steps = 500;}

      // ogni 1000 passi controlla passi accettati su passi fatti
      if ((i > 0) and (i % 1000 == 0))
      {
        // cambio delta se passi accettati non in intervallo [35%,65%]
        if ((fw_steps > 650) or (fw_steps < 350))
        {
          delta = delta*(1 + 0.25*(fw_steps/1000. - 0.5));
        }
        // reset counter
        fw_steps = 0;
      }
    }

    // equilibrio
    if (i >= E)
    {
      // calcolo valori osservabili passo corrente (i)
      Ui[i-E] = H;

      Mui[i-E][0] = 0; Mui[i-E][1] = 0;
      for (int j=0; j < N; j++)
      {
        Mui[i-E][0] += cos(s[j/n/n%n][j/n%n][j%n]);
        Mui[i-E][1] += sin(s[j/n/n%n][j/n%n][j%n]);
      }
    }
  }

  double U=0, Mux=0, Muy=0;
  double deltaU=0, deltaMux=0, deltaMuy=0;

  // lunghezza massima calcolo autocorrelazione
  int l = 2000;

  double U0, Mux0, Muy0;
  double corrU=0, corrMux=0, corrMuy=0;

  for (int i=0; i < m; i++)
  {
    // medio valori osservabili
    U += Ui[i] / m; 
    Mux += Mui[i][0] / m;
    Muy += Mui[i][1] / m;

    // punto iniziale correlazione
    if (i % l == 0)
    {
      U0 = Ui[i];
      Mux0 = Mui[i][0];
      Muy0 = Mui[i][1];
    }

    // calcolo correlazione e medio
    corrU += U0 * Ui[i] / m;
    corrMux += Mux0 * Mui[i][0] / m;
    corrMuy += Muy0 * Mui[i][1] / m;
  }

  deltaU = sqrt(fabs(corrU - pow(U,2))) / sqrt(m-1);
  deltaMux = sqrt(fabs(corrMux - pow(Mux,2))) / sqrt(m-1);
  deltaMuy = sqrt(fabs(corrMuy - pow(Muy,2))) / sqrt(m-1);

  double Mag = sqrt(Mux*Mux + Muy*Muy) / N;
  double deltaMag = (Mux*deltaMux + Muy*deltaMuy) / Mag / pow(N,2);
  
  return {U,deltaU,Mag,deltaMag};
}


int main ()
{
  ofstream esame;
  esame.open("esame.csv");
  esame << "Temp,h,U,deltaU,Mag,deltaMag" << endl;

  double U, deltaU;
  double Cv, deltaCv;
  double Mag, deltaMag;

  // campo magnetico esterno costante (rispetto a J)
  // range di h / J; J = 1
  double h[3] = {0, 1, 2};

  // temperatura in kB*T (rispetto a J) 
  // range di kB*T/J; J = 1
  double T[40] = {1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 
                  2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3,
                  3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4,
                  4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5};

  for (int j=0; j < 3; j++)
  {
    for (int i=0; i < 40; i++)
    {
      tie(U,deltaU,Cv,deltaCv,Mag,deltaMag) = Metropolis(T[i], h[j]);

      esame << T[i]   << "," << h[j]     << ",";
      esame << U      << "," << deltaU   << ",";
      esame << Cv     << "," << deltaCv  << ",";
      esame << Mag    << "," << deltaMag << endl;
    }
  }
}