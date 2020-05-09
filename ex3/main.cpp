#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>

#define LOG(x) cout << x << endl;
using namespace std;

const double E = 10.0;
const double r = 0.1;
const double sigma = 0.4;
const double T = 0.25;

const int N = 200; // S steps 
const int M = 2000; // time steps

// matrix A diagonal elements
double d_n (int n, double alpha, double beta)
{
  return 1 - alpha*pow(n,2) - beta;
}


// matrix A upper diagonal elements
double u_n (int n, double alpha, double beta)
{ 
  return 0.5*(alpha*pow(n-1,2) + beta*(n-1));
}


// matrix A lower diagonal elements
double l_n (int n, double alpha, double beta)
{
  return 0.5*(alpha*pow(n+1,2) - beta*(n+1));
}


int main ()
{ 
  double delta_t = T / double(M);
  double alpha = pow(sigma,2) * delta_t;
  double beta = r * delta_t;

  array<double,N+1> vm; // solution array at step m 
  array<double,N+1> vm1; // solution array at step m+1

  ofstream file;
  file.open(to_string(N)+to_string(M)+".csv");

  // init time 0 vector
  for (int n=0; n < N+1; n++)
  {
    vm[n] = max(3*E/N*n - E, 0.0);
    file << to_string(vm[n]) + ",";
  }
  file << endl;

  // boundary conditions
  vm1[0] = 0;
  vm1[N] = vm[N];

  for (int m=0; m < M-1; m++)
  {
    file << to_string(vm1[0]) + ",";
    for (int n=1; n < N; n++)
    {
      vm1[n] = l_n(n-1, alpha, beta)*vm[n-1];
      vm1[n] += d_n(n, alpha, beta)*vm[n];
      vm1[n] += u_n(n+1, alpha, beta)*vm[n+1];
      
      file << to_string(vm1[n]) + ",";
    }
    file << to_string(vm1[N]) + ",";
    file << endl;

    vm = vm1;
  }
  return 0;
}