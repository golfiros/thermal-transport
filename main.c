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
  if(P_V[2] == 0){ return 0; }
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
  double v[3] = {1.0,1.0,0.0};
  double delta = 0.6;
  double temp = 0.01;
  printf("setting up\n");
  setup(t, mu, v, delta, temp);

  pyplot_import("script.py");
  pyplot_figure(1);
  pyplot_xlabel("$V$");
  pyplot_ylabel("$\\mu$");
  pyplot_figure(2);
  pyplot_xlabel("$V$");
  pyplot_ylabel("$I$");
  pyplot_figure(3);
  pyplot_xlabel("$V$");
  pyplot_ylabel("$\\mu$");
  pyplot_figure(4);
  pyplot_xlabel("$V$");
  pyplot_ylabel("$I$");
  
  plot_t plot_mu = create_plot();
  plot_t plot_rest = create_plot();
  const char* colors[5] = {"blue", "red", "green", "yellow", "purple"};
  char label[100];

  printf("doing various v\n");
  double v_cases[5] = {0.0, 0.1, 0.5, 1.0, 5.0};
  for(int i=0;i<5;i++){
    v[2] = v_cases[i];
    setup(t, mu, v, delta, temp);
    pyplot_figure(1);
    sample_plot(plot_mu, 0, 2, find_bias);
    sprintf(label,"%.1f",v_cases[i]);
    pyplot_plot(plot_mu,colors[i],label);
    pyplot_figure(2);
    int np = get_plot_points(plot_mu);
    double *x = get_plot_x(plot_mu);
    double *bias = get_plot_y(plot_mu);
    double *y = get_plot_y(plot_rest);
    set_plot_x(plot_rest, x, np);
    for(int j=0;j<np;j++){
      P_BIAS[0] = P_MU[0] + 0.5 * x[j];
      P_BIAS[1] = P_MU[1] - 0.5 * x[j];
      P_BIAS[2] = bias[j];
      y[j] = total_current();
    }
    pyplot_plot(plot_rest,colors[i],label);
  }
  
  printf("doing various T\n");
  v[2] = 0.1;
  double temp_cases[5] = {0.001, 0.01, 0.1, 1.0, 10.0};
  for(int i=0;i<5;i++){
    temp = temp_cases[i];
    setup(t, mu, v, delta, temp);
    pyplot_figure(3);
    sample_plot(plot_mu, 0, 2, find_bias);
    sprintf(label,"%.1f",temp_cases[i]);
    pyplot_plot(plot_mu,colors[i],label);
    pyplot_figure(4);
    int np = get_plot_points(plot_mu);
    double *x = get_plot_x(plot_mu);
    double *bias = get_plot_y(plot_mu);
    double *y = get_plot_y(plot_rest);
    set_plot_x(plot_rest, x, np);
    for(int j=0;j<np;j++){
      P_BIAS[0] = P_MU[0] + 0.5 * x[j];
      P_BIAS[1] = P_MU[1] - 0.5 * x[j];
      P_BIAS[2] = bias[j];
      y[j] = total_current();
    }
    pyplot_plot(plot_rest,colors[i],label);
  }
  pyplot_figure(1);
  pyplot_legend("$v_3$");
  pyplot_savefig("mu_versus_voltage_varying_v.png");
  pyplot_figure(2);
  pyplot_legend("$v_3$");
  pyplot_savefig("curr_versus_voltave_varying_v.png");
  
  pyplot_figure(3);
  pyplot_legend("$T$");
  pyplot_savefig("mu_versus_voltage_varying_T.png");
  pyplot_figure(4);
  pyplot_legend("$T$");
  pyplot_savefig("curr_versus_voltave_varying_T.png");
  delete_plot(plot_mu);
  delete_plot(plot_rest);
  pyplot_close();

  return 0;
}
