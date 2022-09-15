#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "numerical_methods.h"

int num_error = 0;
void reset_error_numerics(){
  num_error = 0;
}
int get_error_numerics(){
  return num_error;
}

double integrate(double a, double b, double (*func)(double)){
  if (a == b) { return 0; }
  double INTEGRAL_STACK[2 * MAX_REC];
  INTEGRAL_STACK[0] = a;
  INTEGRAL_STACK[1] = b;
  double tol = TOL / fabs(b - a);
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
    if(fabs(i1-i2) < 3 * fabs(b - a) * tol){ integral += i2; continue; }
    else if(sp > MAX_REC - 2){ integral += i2; num_error |= ERROR_CONVERGENCE; continue; }
    INTEGRAL_STACK[2 * sp] = a;
    INTEGRAL_STACK[2 * sp++ + 1] = m;
    INTEGRAL_STACK[2 * sp] = m;
    INTEGRAL_STACK[2 * sp++ + 1] = b;
  }
  return integral;
}

double find_root(double a, double b, double (*func)(double)){
  double c; //intermediate guess
  double fa, fb, fc; //function evaluations
  int np = 0; //iteration counter
  fa = (*func)(a);
  fb = (*func)(b);
  if(fa * fb > 0){ num_error |= ERROR_INVALID; return b; }
  while(fabs(a - b) > TOL && np++ < MAX_REC){
    c = (a + b) * 0.5;
    fc = (*func)(c);
    if(fa * fc > 0){ a = c; fa = fc; }
    else { b = c; fb = fc; }
  }
  if(np > MAX_REC){ num_error |= ERROR_CONVERGENCE; }
  return b;
}

struct plot_s {
  int n_points;
  double x[MAX_REC];
  double y[MAX_REC];
};

FILE *fptr;

void pyplot_import(const char* filename){
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
  fprintf(fptr,"plt.savefig(\"%s\")\n",filename);
}

void pyplot_plot(plot_t plot, const char* color,const char* label){
  fprintf(fptr,"plt.plot([");
  for(int i=0;i<plot->n_points;i++){
    fprintf(fptr,"%e,",plot->x[i]);
  }
  fprintf(fptr,"],[");
  for(int i=0;i<plot->n_points;i++){
    fprintf(fptr,"%e,",plot->y[i]);
  }
  fprintf(fptr,"],c=\"%s\",label=\"%s\")\n",color,label);
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

plot_t create_plot(){
  struct plot_s *ptr = malloc(sizeof(struct plot_s));
  return ptr;
}

void delete_plot(plot_t plot){
  free(plot);
}

void sample_plot(plot_t plot, double a, double b, double (*func)(double)){
  plot->n_points = N_INITIAL;
  for(int i=0;i<N_INITIAL;i++){
    plot->x[i] = a + (double)i * (b - a) / ((double)N_INITIAL - 1); 
    plot->y[i] = (*func)(plot->x[i]);
  }
  for(int i=0;i<plot->n_points-2 && plot->n_points < MAX_REC;i++){
    double *x = plot->x + i;
    double *y = plot->y + i;
    double ymax = y[0], ymin = y[0];
    if(y[1] > ymax) { ymax = y[1]; }
    if(y[2] > ymax) { ymax = y[2]; }
    if(y[1] < ymin) { ymin = y[1]; }
    if(y[2] < ymin) { ymin = y[2]; }
    double dx = x[2] - x[0];
    double dy = ymax - ymin;
    double area = fabs(x[0] * (y[1] - y[2]) + 
                       x[1] * (y[2] - y[0]) +
                       x[2] * (y[0] - y[1]));
    //printf("%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",x[0],x[1],x[2],area,dx*dy);
    if(area > 0.5 * dx * dy){
      memmove(x + 2, x, (plot->n_points - i) * sizeof(double));
      memmove(y + 2, y, (plot->n_points - i) * sizeof(double));
      x[2] = x[1];
      y[2] = y[1];
      x[1] = 0.5 * (x[0] + x[2]);
      x[3] = 0.5 * (x[2] + x[4]);
      y[1] = (*func)(x[1]);
      y[3] = (*func)(x[3]);
      plot->n_points += 2;
      i--;
    }
  }
}

int get_plot_points(plot_t plot){
  return plot->n_points;
}

double* get_plot_x(plot_t plot){
  return plot->x;
}

void set_plot_x(plot_t plot, const double *x, int n_points){
  plot->n_points = n_points;
  memcpy(plot->x, x, n_points * sizeof(double));
}
double* get_plot_y(plot_t plot){
  return plot->y;
}

void set_plot_y(plot_t plot, const double *y){
  memcpy(plot->y, y, plot->n_points);
}



