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

double P_BIAS[3]; //full conductor chemical potentials
complex double P_MOMENTA[3]; //momentum exponentials
double P_TRANS_PROB[3 * (3 - 1) / 2]; //transmission probabilities

//computes complex exponentials corresponding
//to outgoing momenta at energy E
void find_momenta(double E){
  double x;
  for(int i=0;i<3;i++){
    x = 0.5 * (P_BIAS[i] - E) * P_INV_T[i];
    if (x > 1.0) { P_MOMENTA[i] = x - sqrt(x * x - 1.0); continue; }
    if (x < -1.0) { P_MOMENTA[i] = x + sqrt(x * x - 1.0); continue; }
    P_MOMENTA[i] = x + I * sqrt(1.0 - x * x);
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
  double den = 1.0 / (gf * conj(gf));
  
  double *out = P_TRANS_PROB;
  for(int i=0;i<3-1;i++){
    for(int j=i+1;j<3;j++){
      *(out++) = 4.0 * P_V[i] * P_V[i] * P_V[j] * P_V[j] * P_INV_T[i] * P_INV_T[j] 
           * den * cimag(P_MOMENTA[i]) * cimag(P_MOMENTA[j]);
    }
  }
}

double fermi_dirac(double E){
  return 1.0 / (exp(BETA * E) + 1.0);
}

//compute excess current density flowing into the 
//bath so we can integrate it to zero
double excess_density(double E){
  find_trans(E);
  return M_1_PI * 0.5 * 
        ( P_TRANS_PROB[1] * (fermi_dirac(E - P_BIAS[0]) - fermi_dirac(E - P_BIAS[2]))
        + P_TRANS_PROB[2] * (fermi_dirac(E - P_BIAS[1]) - fermi_dirac(E - P_BIAS[2])));
}

//compute the current density flowing out of
//conductor 1 to get the final current
double current_density(double E){
  find_trans(E);
  return M_1_PI * 0.5 * 
        ( P_TRANS_PROB[0] * (fermi_dirac(E - P_BIAS[0]) - fermi_dirac(E - P_BIAS[1]))
        + P_TRANS_PROB[1] * (fermi_dirac(E - P_BIAS[0]) - fermi_dirac(E - P_BIAS[2])));
}

//compute the noise density of the current
//flowing out of conductor 1
double noise_density(double E){
  find_trans(E);
  return M_1_PI * 0.5 *
       ( P_TRANS_PROB[0] * (fermi_dirac(E - P_BIAS[0]) * (1 - fermi_dirac(E - P_BIAS[1]))
                         +  fermi_dirac(E - P_BIAS[1]) * (1 - fermi_dirac(E - P_BIAS[0])))
       + P_TRANS_PROB[1] * (fermi_dirac(E - P_BIAS[0]) * (1 - fermi_dirac(E - P_BIAS[2])) +
                         +  fermi_dirac(E - P_BIAS[2]) * (1 - fermi_dirac(E - P_BIAS[0])))
       - P_TRANS_PROB[0] * P_TRANS_PROB[0] 
       * (fermi_dirac(E - P_BIAS[0]) - fermi_dirac(E - P_BIAS[1])) 
       * (fermi_dirac(E - P_BIAS[0]) - fermi_dirac(E - P_BIAS[1])) 
       - P_TRANS_PROB[0] * P_TRANS_PROB[1] 
       * (fermi_dirac(E - P_BIAS[0]) - fermi_dirac(E - P_BIAS[1])) 
       * (fermi_dirac(E - P_BIAS[0]) - fermi_dirac(E - P_BIAS[2])) 
       - P_TRANS_PROB[1] * P_TRANS_PROB[0] 
       * (fermi_dirac(E - P_BIAS[0]) - fermi_dirac(E - P_BIAS[2])) 
       * (fermi_dirac(E - P_BIAS[0]) - fermi_dirac(E - P_BIAS[1])) 
       - P_TRANS_PROB[1] * P_TRANS_PROB[1] 
       * (fermi_dirac(E - P_BIAS[0]) - fermi_dirac(E - P_BIAS[2])) 
       * (fermi_dirac(E - P_BIAS[0]) - fermi_dirac(E - P_BIAS[2]))); 
}

//total current flowing into 3
double excess_current(double mu){
  double t, min, max;
  P_BIAS[2] = mu;
  t = 2.0 * P_T[2];
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
//total current noise on conductor 1
double total_noise(){
  double x, t, min, max;
  x = P_BIAS[0];
  t = 2 * P_T[0];
  min = x - t;
  max = x + t;
  return integrate(min, max, noise_density);
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
  double v[3] = {1.0,1.0,0.4};
  double delta = 0.6;
  double temp = 0.025;

  for(int i=0;i<3;i++){
    P_T[i] = t[i];
    P_INV_T[i] = 1 / t[i];
    P_MU[i] = mu[i];
    P_V[i] = v[i];
  }
  DELTA = delta;
  BETA = 1 / temp;
  
  enum {
    MU_VS_VOLT_V,
    CURR_VS_VOLT_V,
    NOISE_VS_VOLT_V,
    MU_VS_VOLT_T,
    CURR_VS_VOLT_T,
    NOISE_VS_VOLT_T
  };

  pyplot_import("script.py");

  pyplot_figure(MU_VS_VOLT_V);
  pyplot_xlabel("$V$");
  pyplot_ylabel("$\\mu$");

  pyplot_figure(CURR_VS_VOLT_V);
  pyplot_xlabel("$V$");
  pyplot_ylabel("$I_1$");
  
  pyplot_figure(NOISE_VS_VOLT_V);
  pyplot_xlabel("$V$");
  pyplot_ylabel("$S_{11}$");

  pyplot_figure(MU_VS_VOLT_T);
  pyplot_xlabel("$V$");
  pyplot_ylabel("$\\mu$");

  pyplot_figure(CURR_VS_VOLT_T);
  pyplot_xlabel("$V$");
  pyplot_ylabel("$I_1$");

  pyplot_figure(NOISE_VS_VOLT_T);
  pyplot_xlabel("$V$");
  pyplot_ylabel("$S_{11}$");

  plot2d_t plot_mu = create_plot2d();
  plot2d_t plot_rest = create_plot2d();
  const char* colors[5] = {"blue", "red", "green", "yellow", "purple"};
  char args[100];

  printf("doing various v\n");
  BETA = 1 / 0.025;
  double v_cases[5] = {0.0, 0.2, 0.4, 0.6, 0.8};
  for(int i=0;i<5;i++){
    P_V[2] = v_cases[i];
    sprintf(args,"c=\"%s\",label=\"%.1f\"",colors[i],v_cases[i]);

    pyplot_figure(MU_VS_VOLT_V);
    sample_plot2d(plot_mu, 0.0, 2.0, find_bias);
    pyplot_plot(plot_mu,args);

    int np = get_plot2d_points(plot_mu);
    double *x = get_plot2d_x(plot_mu);
    double *bias = get_plot2d_y(plot_mu);
    double *y = get_plot2d_y(plot_rest);
    set_plot2d_x(plot_rest, x, np);

    pyplot_figure(CURR_VS_VOLT_V);
    for(int j=0;j<np;j++){
      P_BIAS[0] = P_MU[0] + 0.5 * x[j];
      P_BIAS[1] = P_MU[1] - 0.5 * x[j];
      P_BIAS[2] = bias[j];
      y[j] = total_current();
    }
    pyplot_plot(plot_rest,args);

    pyplot_figure(NOISE_VS_VOLT_V);
    for(int j=0;j<np;j++){
      P_BIAS[0] = P_MU[0] + 0.5 * x[j];
      P_BIAS[1] = P_MU[1] - 0.5 * x[j];
      P_BIAS[2] = bias[j];
      y[j] = total_noise();
    }
    pyplot_plot(plot_rest,args);
  }

  pyplot_figure(MU_VS_VOLT_V);
  pyplot_legend("$v_3$");
  pyplot_savefig("mu_v.png");
  pyplot_figure(CURR_VS_VOLT_V);
  pyplot_legend("$v_3$");
  pyplot_savefig("curr_v.png");
  pyplot_figure(NOISE_VS_VOLT_V);
  pyplot_legend("$v_3$");
  pyplot_savefig("noise_v.png");
  

  printf("doing various T\n");
  P_V[2] = 0.4;
  double temp_cases[5] = {0.001, 0.005, 0.025, 0.125, 0.625};
  for(int i=0;i<5;i++){
    BETA = 1 / temp_cases[i];
    sprintf(args,"c=\"%s\",label=\"%.3f\"",colors[i],temp_cases[i]);

    pyplot_figure(MU_VS_VOLT_T);
    sample_plot2d(plot_mu, 0, 2, find_bias);
    pyplot_plot(plot_mu,args);

    int np = get_plot2d_points(plot_mu);
    double *x = get_plot2d_x(plot_mu);
    double *bias = get_plot2d_y(plot_mu);
    double *y = get_plot2d_y(plot_rest);
    set_plot2d_x(plot_rest, x, np);

    pyplot_figure(CURR_VS_VOLT_T);
    for(int j=0;j<np;j++){
      P_BIAS[0] = P_MU[0] + 0.5 * x[j];
      P_BIAS[1] = P_MU[1] - 0.5 * x[j];
      P_BIAS[2] = bias[j];
      y[j] = total_current();
    }
    pyplot_plot(plot_rest,args);

    pyplot_figure(NOISE_VS_VOLT_T);
    for(int j=0;j<np;j++){
      P_BIAS[0] = P_MU[0] + 0.5 * x[j];
      P_BIAS[1] = P_MU[1] - 0.5 * x[j];
      P_BIAS[2] = bias[j];
      y[j] = total_noise();
    }
    pyplot_plot(plot_rest,args);
 
  }

  pyplot_figure(MU_VS_VOLT_T);
  pyplot_legend("$T$");
  pyplot_savefig("mu_t.png");
  pyplot_figure(CURR_VS_VOLT_T);
  pyplot_legend("$T$");
  pyplot_savefig("curr_t.png");
  pyplot_figure(NOISE_VS_VOLT_T);
  pyplot_legend("$T$");
  pyplot_savefig("noise_t.png");

  delete_plot2d(plot_mu);
  delete_plot2d(plot_rest);
  pyplot_close();

  return 0;
}
