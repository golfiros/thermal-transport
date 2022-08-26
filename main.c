#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <complex.h>
#include <math.h>

  //computes complex exponential corresponding to outgoing
//momentum at energy E on a band with chemical potential mu
//and hopping t
complex double find_momentum(double mu, double t, double E) {
  double x = (mu - E) / (2 * t);
  if (x > 1) { return x - sqrt(x * x - 1); }
  if (x < -1) { return x + sqrt(x * x - 1); }
  return x + I * sqrt(1 - x * x);
}

//computes the transmission probability matrix at energy E 
//between N probes coupled to a central impurity site
//out should be a pointer of length at least N*(N-1)/2

int compute_transmissions(double* out, int N, double* mu, double* t, double* v, double delta, double E) {
  //we can compute all transmissions at once because the
  //expensive bit of the calculation (an absolute value)
  //is the same for all of them
  //the results are ordered in the pointer as 
  //12 13 14 ... 1N 23 24 ... 2N ... (N-1)N

  if(!out) { return 1; } //quit if bad pointer
  
  complex double* momenta = malloc(N * sizeof(complex double));
  if(!momenta) { return 2; } //quit if failed to allocate momenta

  complex double gf = E - delta; //denominator of the amplitude
  for(int i=0;i<N;i++) {
    momenta[i] = find_momentum(mu[i], t[i], E); //compute and store all the momenta
    gf += v[i] * v[i] * momenta[i] / t[i]; //add contributions to the denominator
  }

  complex double den = gf * conj(gf);
  for(int i=0;i<N-1;i++) {
    for(int j=i+1;j<N;j++){
      *out = (4 * v[i] * v[i] * v[j] * v[j] / (t[i] * t[j])) * cimag(momenta[i]) * cimag(momenta[j]) / den;
      out++;
    }
  }
  free(momenta);
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
  complex double momentum;
  for(int i=0; i<40;i++){
    energy = ((double)i+0.5) * 0.2 - 4.0;
    momentum = find_momentum(0, 1, energy);
    printf("energy: %f, momentum exponential: %f+i%f\n", energy, creal(momentum), cimag(momentum));
  }
  */

  //TEST OF compute_transmissions
  /*
  double transmission[3 * 2 / 2]; 
  double mu[3] = {1.0,0.0,-1.0};
  double t[3] = {1.0,1.0,1.0};
  double v[3] = {0.5,0.5,0.5};
  compute_transmissions(transmission, 3, mu, t, v, 0.3, 0.7);
  printf("%f %f %f", transmission[0], transmission[1], transmission[2]);
  */

  return 0;
}
