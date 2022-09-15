#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "numerical_methods.h"

//global variables defining system configuration
//don't forget to setup in main
double P_T[3]; //hopping parameters
double P_INV_T[3]; //inverse hopping parameters
double P_MU[3]; //conductor bare chemical potentials
double P_V[3]; //impurity couplings
double DELTA; //impurity site energy
double BETA; //inverse temperature

void setup(double *t, double *mu, double *v, double delta, double temp){
  for(int i=0;i<3;i++){
    P_T[i] = *t;
    P_INV_T[i] = 1 / *(t++);
    P_MU[i] = *(mu++);
    P_V[i] = *(v++);
  }
  DELTA = delta;
  BETA = 1 / temp;
}

double P_BIAS[3]; //full conductor chemical potentials
complex double P_MOMENTA[3]; //momentum exponentials
double P_TRANS_PROB[3 * (3 - 1) / 2]; //transmission probabilities

//computes complex exponentials corresponding
//to outgoing momenta at energy E
void find_momenta(double E){
  double x;
  for(int i=0;i<3;i++){
    x = 0.5 * (P_BIAS[i] - E) * P_INV_T[i];
    if (x > 1) { P_MOMENTA[i] = x - sqrt(x * x - 1); continue; }
    if (x < -1) { P_MOMENTA[i] = x + sqrt(x * x - 1); continue; }
    P_MOMENTA[i] = x + I * sqrt(1 - x * x);
  }
}

//computes the transmission probabilities at energy E 
//the results are ordered in the pointer as 12 13 23
void find_trans(double E){
  find_momenta(E);
  complex double gf = E - DELTA; //denominator of the amplitude
  for(int i=0;i<3;i++){
    gf += P_V[i] * P_V[i] * P_MOMENTA[i] * P_INV_T[i]; //add contributions to the denominator
  }
  double den = 1 / (gf * conj(gf));
  
  double *out = P_TRANS_PROB;
  for(int i=0;i<3-1;i++){
    for(int j=i+1;j<3;j++){
      *(out++) = 4 * P_V[i] * P_V[i] * P_V[j] * P_V[j] * P_INV_T[i] * P_INV_T[j] 
           * den * cimag(P_MOMENTA[i]) * cimag(P_MOMENTA[j]);
    }
  }
}

double fermi_dirac(double E){
  return 1 / (exp(BETA * E) + 1);
}

//compute excess current density flowing into the 
//bath so we can integrate it to zero
double excess_density(double E){
  find_trans(E);
  return ( P_TRANS_PROB[1] * (fermi_dirac(E - P_BIAS[0]) - fermi_dirac(E - P_BIAS[2]))
         + P_TRANS_PROB[2] * (fermi_dirac(E - P_BIAS[1]) - fermi_dirac(E - P_BIAS[2]))) / M_1_PI * 0.5;
}

//compute the current density flowing out of
//conductor 1 to get the final current
double current_density(double E){
  find_trans(E);
  return ( P_TRANS_PROB[0] * (fermi_dirac(E - P_BIAS[0]) - fermi_dirac(E - P_BIAS[1]))
         + P_TRANS_PROB[1] * (fermi_dirac(E - P_BIAS[0]) - fermi_dirac(E - P_BIAS[2]))) * M_1_PI * 0.5;
}

//total current flowing into 3
double excess_current(double mu){
  double t, min, max;
  P_BIAS[2] = mu;
  t = 2 * P_T[2];
  min = mu - t;
  max = mu + t;
  double out = integrate(min, max, excess_density);
  reset_error_numerics();
  return out;
}
//total current flowing out of 1
double total_current(){
  double x, t, min, max;
  x = P_BIAS[0];
  t = 2 * P_T[0];
  min = x - t;
  max = x + t;
  return integrate(min, max, current_density);
}

//procedure to find chemical potential of the
//bath conductor by root finding with secant
double find_bias(double voltage) {
  P_BIAS[0] = P_MU[0] + 0.5 * voltage;
  P_BIAS[1] = P_MU[1] - 0.5 * voltage;
  return find_root(P_BIAS[0], P_BIAS[1], excess_current);
}

int main(int argc, char** argv){
  //this program evaluates the current between two probes 
  //coupled to an impurity site, which is then coupled
  //to a third probe that functions as a bath, making the
  //process inelastic in nature.
  
  double t[3] = {1.0,1.0,1.0}; 
  double mu[3] = {0.0,0.0,0.0};
  double v[3] = {1.0,1.0,0.2};
  double delta = 0.6;
  double temp = 0.025;
  setup(t, mu, v, delta, temp);
 
  /*
  int steps = 100;

  FILE *ptr = fopen("out.tsv","w");
  for(int i=0;i<steps;i++){
    reset_error_numerics();
    double voltage = 4.0 * ((double)i + 0.5) / steps;
    double bias = find_bias(voltage);
    double curr = total_current();
    //fprintf(ptr,"%f\t%f\t%f\n",voltage,bias,curr);
    printf("%f\t%f\t%f\t%o\n",voltage,bias,curr,get_error_numerics());
  }
  fclose(ptr);
  */

  plot_t plot = create_plot(0, 4);
  sample_plot(plot, find_bias);
  render_plot(plot, "");
  delete_plot(plot);

  return 0;
}
