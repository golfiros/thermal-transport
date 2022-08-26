#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <complex.h>
#include <math.h>

complex double find_momentum(double mu, double t, double E) {
  double x = (mu - E) / (2 * t);
  if (x > 1) {
    return x - sqrt(x * x - 1);
  }
  if (x < -1) {
    return x + sqrt(x * x - 1);
  }
  return x + I * sqrt(1 - x * x);
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

  return 0;
}
