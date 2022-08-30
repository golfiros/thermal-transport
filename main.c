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

double P_BIAS[3]; //full conductor chemical potentials
complex double P_MOMENTA[3]; //momentum exponentials
double P_TRANS_PROB[3 * (3 - 1) / 2]; //transmission probabilities

double E_MIN, E_MAX; //integration bounds given by bands

void setup(double *t, double *mu, double *v, double delta){
  for(int i=0;i<3;i++){
    P_T[i] = *t;
    P_INV_T[i] = 1 / *(t++);
    P_MU[i] = *(mu++);
    P_V[i] = *(v++);
  }
  DELTA = delta;
}

void set_bias(double *bias){
  double x, t, min, max;
  for(int i=0;i<3;i++){
    x = P_MU[i] + *(bias++);
    P_BIAS[i] = x;
    t = 2 * P_T[i];
    min = x - t;
    max = x + t;
    if(min < E_MIN) { E_MIN = min; }
    if(max > E_MAX) { E_MAX = max; }
  }
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

void print_momenta(){
  printf("momenta: ");
  for(int j=0;j<3;j++){
    printf("%1.2f + %1.2fi\t", creal(P_MOMENTA[j]), cimag(P_MOMENTA[j]));
  }
}

void test_find_momenta(int steps){
  double energy;
  for(int i=0;i<steps;i++){
    energy = E_MIN + (E_MAX - E_MIN) * ((double)i + 0.5) / steps;
    find_momenta(energy);
    printf("energy: %2.2f\t", energy);
    print_momenta();
    printf("\n");
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
      *out = 4 * P_V[i] * P_V[i] * P_V[j] * P_V[j] * 
        P_INV_T[i] * P_INV_T[j] * den * cimag(P_MOMENTA[i]) * cimag(P_MOMENTA[j]);
      out++;
    }
  }
}

void print_trans(){
  printf("transmissions: ");
  for(int j=0;j<3;j++){
    printf("%1.5f\t", P_TRANS_PROB[j]);
  }
}

void test_find_trans(int steps){
  double energy;
  for(int i=0;i<steps;i++){
    energy = E_MIN + (E_MAX - E_MIN) * ((double)i + 0.5) / steps;
    printf("energy: %2.2f\t", energy);
    find_trans(energy);
    print_momenta();
    print_trans();
    printf("\n");
  }
}

double fermi_dirac(double beta, double E){
  return 1 / (exp(beta * E) + 1);
}

//numerical integration by adaptive quadrature
//accurate within |b - a| * tol
#define MAX_REC 512 //max subdivisions for integration
double INTEGRAL_STACK[2 * MAX_REC];
double integrate(double (*func)(double), double a, double b, double tol){
  double integral = 0.0;
  double i1, i2, m;
  INTEGRAL_STACK[0] = a;
  INTEGRAL_STACK[1] = b;
  int sp = 1;
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

int main(int argc, char** argv){
  //this program evaluates the current between two probes 
  //coupled to an impurity site, which is then coupled
  //to a third probe that functions as a bath, making the
  //process inelastic in nature.
  
  /*
  double t[3] = {1.0,1.0,1.0}; 
  double mu[3] = {0.0,0.0,0.0};
  double v[3] = {0.2,0.3,0.1};
  double delta = 0.3;
  setup(t, mu, v, delta);

  double bias[3] = {0.5,-0.5,0.1};
  set_bias(bias);
  test_find_momenta(20);
  test_find_trans(50);
  */

  return 0;
}
