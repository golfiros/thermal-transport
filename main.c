#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

complex double* P_MOMENTA; //dynamic pointer for allocating momentum exponentials

  //computes complex exponential corresponding to outgoing
//momentum at energy E on a band with chemical potential mu
//and hopping t
complex double find_momentum(double mu, double inv_t, double E) {
  double x = 0.5 * (mu - E) * inv_t;
  if (x > 1) { return x - sqrt(x * x - 1); }
  if (x < -1) { return x + sqrt(x * x - 1); }
  return x + I * sqrt(1 - x * x);
}

//computes the transmission probability matrix at energy E 
//between N probes coupled to a central impurity site
//out should be a pointer of length at least N*(N-1)/2

int compute_transmissions(double* out, int N, double* mu, double* inv_t, double* v, double delta, double E) {
  //we can compute all transmissions at once because the
  //expensive bit of the calculation (an absolute value)
  //is the same for all of them
  //the results are ordered in the pointer as 
  //12 13 14 ... 1N 23 24 ... 2N ... (N-1)N

  if(!out) { return 1; } //quit if bad pointer
  
  complex double gf = E - delta; //denominator of the amplitude
  for(int i=0;i<N;i++) {
    P_MOMENTA[i] = find_momentum(mu[i], inv_t[i], E); //compute and store all the momenta
    gf += v[i] * v[i] * P_MOMENTA[i] * inv_t[i]; //add contributions to the denominator
  }

  double den = 1 / (gf * conj(gf));
  for(int i=0;i<N-1;i++) {
    for(int j=i+1;j<N;j++){
      *out = 4 * v[i] * v[i] * v[j] * v[j] * inv_t[i] * inv_t[j] * den * cimag(P_MOMENTA[i]) * cimag(P_MOMENTA[j]);
      out++;
    }
  }
  return 0;
}

int main(int argc, char** argv) {
  //this program evaluates the current between two probes 
  //coupled to an impurity site, which is then coupled
  //to a third probe that functions as a bath, making the
  //process inelastic in nature.
  
  //TEST OF find_momentum 
  /* 
  double energy;
  double E_min = -6.0;
  double E_max = 6.0;
  int steps = 20;
  complex double momentum;
  for(int i=0;i<steps;i++){
    energy = E_min + (E_max - E_min) * ((double)i+0.5) / steps;
    momentum = find_momentum(0, 1, energy);
    printf("energy: %f, momentum exponential: %f+i%f\n", energy, creal(momentum), cimag(momentum));
  }
  */

  //TEST OF compute_transmissions
  /*
  #define N_CONDUCTORS 4
  P_MOMENTA = malloc(N_CONDUCTORS * sizeof(complex double));
  double transmission[N_CONDUCTORS * (N_CONDUCTORS - 1) / 2]; 
  double mu[N_CONDUCTORS] = {1.0,0.5,-0.5,-1.0};
  double inv_t[N_CONDUCTORS] = {1.2,0.8,0.4,0.6};
  double v[N_CONDUCTORS] = {0.3,0.4,0.5,0.6};
  double delta = 0.3;
  double energy;
  double E_min = -4.0;
  double E_max = 4.0;
  int steps = 20;
  for(int i=0;i<steps;i++){
    energy = E_min + (E_max - E_min) * ((double)i+0.5) / steps;
    printf("Energy: %f, Transmissions: ", energy);
    compute_transmissions(transmission, N_CONDUCTORS, mu, inv_t, v, delta, energy);
    for(int j=0;j<N_CONDUCTORS*(N_CONDUCTORS-1)/2;j++){
      printf("%f ", transmission[j]);
    }
    printf("\n");
  }
  free(P_MOMENTA);
  #undef N_CONDUCTORS
  */

  return 0;
}
