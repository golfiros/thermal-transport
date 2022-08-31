#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

//global parameters for system configuration
//remember to set up in main

double P_T[3]; //hopping parameters
double P_INV_T[3]; //inverse hopping parameters
double P_MU[3]; //conductor bare chemical potentials
double P_V[3]; //impurity couplings
double DELTA; //impurity site energy
double BETA; //inverse temperature

double P_BIAS[3]; //full conductor chemical potentials
complex double P_MOMENTA[3]; //momentum exponentials
double P_TRANS_PROB[3 * (3 - 1) / 2]; //transmission probabilities

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
      *out = 4 * P_V[i] * P_V[i] * P_V[j] * P_V[j] * P_INV_T[i] * P_INV_T[j] 
           * den * cimag(P_MOMENTA[i]) * cimag(P_MOMENTA[j]);
      out++;
    }
  }
}

//numerical integration by adaptive quadrature
//accurate within |b-a|*tol
#define MAX_REC 512 //max subdivisions for integration
#define N_INITIAL 16 //initial even subdivisions
double INTEGRAL_STACK[2 * MAX_REC];
double integrate(double (*func)(double), double a, double b, double tol){
  INTEGRAL_STACK[0] = a;
  INTEGRAL_STACK[1] = b;
  int sp;
  for(sp = 0; sp < N_INITIAL; sp++){
    INTEGRAL_STACK[2 * sp] = a + (b - a) * (double)sp / N_INITIAL;
    INTEGRAL_STACK[2 * sp + 1] = a + (b - a) * (double)(sp + 1) / N_INITIAL;
  }
  double integral = 0.0;
  double i1, i2, m;
  while(sp > 0){
    b = INTEGRAL_STACK[2 * --sp + 1];
    a = INTEGRAL_STACK[2 * sp];
    i1 = 0.5 * (b - a) * ((*func)(a) + (*func)(b));
    m = 0.5 * (a + b);
    i2 = 0.25 * (b - a) * ((*func)(a) + 2 * (*func)(m) + (*func)(b));
    if(fabs(i1-i2) < 3 * fabs(b - a) * tol || sp > MAX_REC - 2){ integral += i2; continue; }
    INTEGRAL_STACK[2 * sp] = a;
    INTEGRAL_STACK[2 * sp++ + 1] = m;
    INTEGRAL_STACK[2 * sp] = m;
    INTEGRAL_STACK[2 * sp++ + 1] = b;
  }
  return integral;
}
#undef MAX_REC
#undef N_INITIAL

double fermi_dirac(double E){
  return 1 / (exp(BETA * E) + 1);
}

//compute excess current density flowing into the 
//bath so we can integrate it to zero
double excess_density(double E){
  find_trans(E);
  return P_TRANS_PROB[1] * (fermi_dirac(E - P_BIAS[0]) - fermi_dirac(E - P_BIAS[2]))
       + P_TRANS_PROB[2] * (fermi_dirac(E - P_BIAS[1]) - fermi_dirac(E - P_BIAS[2]));
}

//compute the current density flowing out of
//conductor 1 to get the final current
double current_density(double E){
  find_trans(E);
  return P_TRANS_PROB[0] * (fermi_dirac(E - P_BIAS[0]) - fermi_dirac(E - P_BIAS[1]))
       + P_TRANS_PROB[1] * (fermi_dirac(E - P_BIAS[0]) - fermi_dirac(E - P_BIAS[2]));
}


#define TOL 0.001 //tolerance for all calculations

//total current flowing into 3
double excess_current(){
  double x, t, min, max;
  x = P_BIAS[2];
  t = 2 * P_T[2];
  min = x - t;
  max = x + t;
  return integrate(excess_density, min, max, TOL / (max - min));
}
//total current flowing out of 1
double total_current(){
  double x, t, min, max;
  x = P_BIAS[0];
  t = 2 * P_T[0];
  min = x - t;
  max = x + t;
  return integrate(current_density, min, max, TOL / (max - min));
}

//procedure to find chemical potential of the
//bath conductor by root finding with secant
double find_bias(double voltage) {
  P_BIAS[0] = P_MU[0] + voltage / 2;
  P_BIAS[1] = P_MU[1] - voltage / 2;
  double x0, x1, x2; //guesses
  double f0, f1; //evaluations
  x0 = 0;
  P_BIAS[2] = P_MU[2] + x0;
  f0 = excess_current();
  x1 = (f0 < 0) ? 10 * TOL : - 10 * TOL; //if current flowing in, raise voltage
  P_BIAS[2] = P_MU[2] + x1;
  f1 = excess_current();
  while(fabs(x1-x0) > TOL){
    //printf("%f\t%f\t%f\t%f\t%f\n",voltage,f0,f1,x0,x1);
    x2 = (x0 * f1 - x1 * f0) / (f1 - f0);
    f0 = f1;
    P_BIAS[2] = P_MU[2] + x2;
    f1 = excess_current();
    x0 = x1;
    x1 = x2;
  }
  return x1;
}

int main(int argc, char** argv){
  //this program evaluates the current between two probes 
  //coupled to an impurity site, which is then coupled
  //to a third probe that functions as a bath, making the
  //process inelastic in nature.
  
  double t[3] = {1.0,1.0,1.0}; 
  double mu[3] = {0.0,0.0,0.0};
  double v[3] = {0.2,0.2,0.2};
  double delta = 0.6;
  double temp = 1.0;
  setup(t, mu, v, delta, temp);

  int steps = 50;
  
  FILE *ptr = fopen("out.tsv","w");
  for(int i=0;i<steps;i++){
    double voltage = 4.0 * ((double)i + 0.5) / steps;
    double bias = find_bias(voltage);
    P_BIAS[0] = P_MU[0] + voltage / 2;
    P_BIAS[1] = P_MU[1] - voltage / 2;
    P_BIAS[2] = P_MU[2] + bias;
    double curr = total_current();
    fprintf(ptr,"%f\t%f\t%f\n",voltage,bias,curr);
  }
  fclose(ptr);

  return 0;
}
