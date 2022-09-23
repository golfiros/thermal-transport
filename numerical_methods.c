#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "numerical_methods.h"

#define N_MAX_DEPTH N_INITIAL << MAX_REC_DEPTH

int num_error = 0;
void reset_error_numerics(){
  num_error = 0;
}
int get_error_numerics(){
  return num_error;
}

double integrate(double a, double b, double (*func)(double)){
  if (fabs(a - b) < TOL) { return 0; }
  int depths[N_MAX_DEPTH] = { 0 };
  double stack[N_MAX_DEPTH][4];
  double tol = TOL / fabs(b - a);
  int sp;
  for(sp = 0; sp < N_INITIAL; sp++){
    stack[sp][0] = a + (b - a) * (double)sp / N_INITIAL;
    stack[sp][1] = a + (b - a) * (double)(sp + 1) / N_INITIAL;
    stack[sp][2] = func(stack[sp][0]);
    stack[sp][3] = func(stack[sp][1]);
  }
  double integral = 0.0;
  double m;
  double fa, fb, fm;
  double i1, i2;
  while(sp > 0){
    a = stack[--sp][0];
    b = stack[sp][1];
    fb = stack[sp][3];
    fa = stack[sp][2];
    m = 0.5 * (a + b);
    fm = func(m);
    i1 = 0.5 * (b - a) * (fa + fb);
    i2 = 0.25 * (b - a) * (fa + 2.0 * fm + fb);
    if(fabs(i1-i2) < 3.0 * fabs(b - a) * tol){
      integral += i2;
      depths[sp] = 0;
      continue;
    }
    else if(depths[sp] >= MAX_REC_DEPTH){
      integral += i2;
      depths[sp] = 0;
      num_error |= ERROR_CONVERGENCE; 
      continue;
    }
    depths[sp]++;
    //stack[sp][0] = a;
    stack[sp][1] = m;
    //stack[sp][2] = fa;
    stack[sp++][3] = fm;
    depths[sp] = depths[sp - 1];
    stack[sp][0] = m;
    stack[sp][1] = b;
    stack[sp][2] = fm;
    stack[sp++][3] = fb;
  }
  return integral;
}

double find_root(double a, double b, double (*func)(double)){
  int i = 0;
  double c, d, e;
  double fa, fb, fc;
  fa = func(a);
  fb = func(b);
  if(fa * fb > 0){ num_error |= ERROR_INVALID; return b; }
  c = a;
  fc = fa;
  d = b - c;
  e = d;
  while(fabs(fb) > TOL && i++ < N_MAX_DEPTH){
    if(fa * fb > 0){
      a = c; fa = fc;
      d = b - c; e = d;
    }
    if(fabs(fa) < fabs(fb)){
      c = b; b = a; a = c;
      fc = fb; fb = fa; fa = fc;
    }
    double m = 0.5 * (a - b);
    double tol = TOL * (fabs(b) + 1 + fabs(fabs(b) - 1.0));
    if(fabs(m) < tol || fb < TOL){ break; }
    if(fabs(e) < tol || fabs(fc) < fabs(fb)){
      d = m;
      e = m;
    }
    else{
      double p, q, r, s;
      s = fb / fc;
      if(fabs(a - c) < TOL){
        p = 2.0 * m * s;
        q = 1.0 - s;
      }
      else{
        q = fc / fa;
        r = fb / fa;
        p = s * (2.0 * m * q * (q - r) - (b - c) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if(p > 0){ q = -q; }
      else{ p = -p; }
      if(2.0 * p < 3.0 * m * q - fabs(tol * q) && p < fabs(0.5 * e * q)){
        e = d;
        d = p / q;
      }
      else{
        d = m;
        e = m;
      }
    }
    c = b;
    fc = fb;
    if(fabs(d) > tol){
      b = b + d;
    }
    else{
      b = b - copysign(tol, b - a);
    }
    fb = func(b);
  }
  return b;
}

FILE *fptr;

void pyplot_open(const char* filename){
  fptr = fopen(filename,"w");
  fprintf(fptr,"import matplotlib.pyplot as plt\n");
}

void pyplot_close(){
  fclose(fptr);
}

void pyplot_figure(int figure){
  fprintf(fptr,"plt.figure(%d)\n",figure);
}

void pyplot_savefig(const char* filename){
  fprintf(fptr,"plt.savefig(\"%s\",dpi=150)\n",filename);
}

void pyplot_xlabel(const char* label){
  fprintf(fptr,"plt.xlabel(\"%s\")\n",label);
}

void pyplot_ylabel(const char* label){
  fprintf(fptr,"plt.ylabel(\"%s\")\n",label);
}

void pyplot_legend(const char* title){
  fprintf(fptr,"plt.legend(title=\"%s\")\n",title);
}

struct plot2d_s {
  int n_points;
  double x[N_MAX_DEPTH];
  double y[N_MAX_DEPTH];
};


plot2d_t create_plot2d(){
  struct plot2d_s *ptr = malloc(sizeof(struct plot2d_s));
  return ptr;
}

void delete_plot2d(plot2d_t plot){
  free(plot);
}

int comp_plot2d(const void *a, const void *b){
  const double *la = a;
  const double *lb = b;
  if(*la >= *lb){ return 1; }
  else{ return -1; }
}

void function_plot2d(plot2d_t plot, double a, double b, double (*func)(double)){
  if (fabs(a - b) < TOL) { return; }
  int depths[N_MAX_DEPTH] = { 0 };
  double intervals[N_MAX_DEPTH][4];
  double tol = 10000.0 * TOL / fabs(b - a); //we're way more lenient with plot precision
  int np;
  for(np = 0; np < N_INITIAL - 1; np++){
    intervals[np][0] = a + (b - a) * (double)np / (N_INITIAL - 1);
    intervals[np][1] = a + (b - a) * (double)(np + 1) / (N_INITIAL - 1);
    intervals[np][2] = func(intervals[np][0]);
    intervals[np][3] = func(intervals[np][1]);
  }
  intervals[N_INITIAL - 1][0] = intervals[N_INITIAL - 1][1] = b;
  intervals[N_INITIAL - 1][2] = intervals[N_INITIAL - 1][3] = func(b);
  np++;

  double m;
  double fa, fb, fm;
  double i1, i2;
  int j = 0;
  while(j < np){ //go through each interval
    a = intervals[j][0];
    b = intervals[j][1];
    fa = intervals[j][2];
    fb = intervals[j][3];
    m = 0.5 * (a + b);
    fm = func(m);
    i1 = 0.5 * (b - a) * (fa + fb);
    i2 = 0.25 * (b - a) * (fa + 2.0 * fm + fb);
    if(fabs(i1-i2) <= 3.0 * fabs(b - a) * tol){
      j++;
      continue;
    }
    else if(depths[j] >= MAX_REC_DEPTH){
      j++;
      num_error |= ERROR_CONVERGENCE; 
      continue;
    }
    depths[j]++;
    //intervals[j][0] = a;
    intervals[j][1] = m;
    //intervals[j][2] = fa;
    intervals[j][3] = fm;
    depths[np] = depths[j];
    intervals[np][0] = m;
    intervals[np][1] = b;
    intervals[np][2] = fm;
    intervals[np][3] = fb;
    np++;
  }
  //for line plots this needs to be sorted, so let's do that
  qsort(intervals, np, sizeof(*intervals), comp_plot2d);
  plot->n_points = np;
  for(j=0;j<np;j++){
    plot->x[j] = intervals[j][0];
    plot->y[j] = intervals[j][2];
  }
}

int get_plot2d_points(plot2d_t plot){
  return plot->n_points;
}

double* get_plot2d_x(plot2d_t plot){
  return plot->x;
}

void set_plot2d_x(plot2d_t plot, const double *x, int n_points){
  plot->n_points = n_points;
  memcpy(plot->x, x, n_points * sizeof(double));
}
double* get_plot2d_y(plot2d_t plot){
  return plot->y;
}

void set_plot2d_y(plot2d_t plot, const double *y){
  memcpy(plot->y, y, plot->n_points);
}

void pyplot_plot(plot2d_t plot, const char* args){
  fprintf(fptr,"plt.plot([");
  for(int i=0;i<plot->n_points;i++){
    fprintf(fptr,"%e,",plot->x[i]);
  }
  fprintf(fptr,"],[");
  for(int i=0;i<plot->n_points;i++){
    fprintf(fptr,"%e,",plot->y[i]);
  }
  fprintf(fptr,"],%s)\n",args);
}



