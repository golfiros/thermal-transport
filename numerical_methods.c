#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "numerical_methods.h"

#define N_MAX_DEPTH (N_INITIAL << MAX_REC_DEPTH)

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
  double points[N_MAX_DEPTH][2];
  int intervals[N_MAX_DEPTH][2];
  double tol = 100000.0 * TOL / fabs(b - a); //we're way more lenient with plot precision
  int np;
  for(np = 0; np < N_INITIAL; np++){
    points[np][0] = a + (b - a) * (double)np / (N_INITIAL - 1);
    points[np][1] = func(points[np][0]);
    intervals[np][0] = np;
    intervals[np][1] = np + 1;
  }
  double m;
  double fa, fb, fm;
  double i1, i2;
  int j = 0;
  while(j < np - 1){ //go through each interval
    a = points[intervals[j][0]][0];
    b = points[intervals[j][1]][0];
    fa = points[intervals[j][0]][1];
    fb = points[intervals[j][1]][1];
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
    points[np][0] = m;
    points[np][1] = fm;
    intervals[np][0] = np;
    intervals[np][1] = intervals[j][1];
    intervals[j][1] = np;
    np++;
  }
  //for line plots this needs to be sorted, so let's do that
  qsort(points, np, sizeof(*points), comp_plot2d);
  plot->n_points = np;
  for(j=0;j<np;j++){
    plot->x[j] = points[j][0];
    plot->y[j] = points[j][1];
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

struct plot3d_s {
  int n_points;
  double xy[N_MAX_DEPTH][2];
  double z[N_MAX_DEPTH];
};

plot3d_t create_plot3d(){
  struct plot3d_s *ptr = malloc(sizeof(struct plot3d_s));
  return ptr;
}

void delete_plot3d(plot3d_t plot){
  free(plot);
}

int intcmp(const void *a, const void *b){
  const int *va = a;
  const int *vb = b;
  return va - vb;
}

void triangulate(double (*points)[], int* n_entries, int n, int (*triangulation)[], int* n_triangles){
    
}

void function_plot3d(plot3d_t plot, double xa, double xb, double ya, double yb, double (*func)(double,double)){
  if (fabs(xa - xb) < TOL || fabs(ya - yb) < TOL) { return; }
  double points[N_MAX_DEPTH][3];
  int triangles[2 * N_MAX_DEPTH][3]; //more than enough memory
  int to_refine[2 * N_MAX_DEPTH];
  double tol = 100000.0 * TOL / (fabs(xb - xa) * fabs(yb - ya)); 
  const int n_initial = (int)sqrt(N_INITIAL);
  const int max_rec_depth;

  int np, nt;
  for(int j=0;j<n_initial;j++){
    for(int i=0;i<n_initial;i++){
      points[n_initial * j + i][0] = xa + (xb - xa) * (double)i / n_initial;
      points[n_initial * j + i][1] = ya + (yb - ya) * (double)j / n_initial;
      points[n_initial * j + i][2] = func(points[n_initial * j + i][0], points[n_initial * j + i][1]);
      np++;
    }
  }
  for(int j=0;j<n_initial-1;j++){
    for(int i=0;i<n_initial-1;i++){
      triangles[2 * (n_initial * j + i)][0] = n_initial * j + i;
      triangles[2 * (n_initial * j + i)][1] = n_initial * (j + 1) + i;
      triangles[2 * (n_initial * j + i)][2] = n_initial * j + i + 1;
      nt++;
      triangles[2 * (n_initial * j + i) + 1][0] = n_initial * j + i + 1;
      triangles[2 * (n_initial * j + i) + 1][1] = n_initial * (j + 1) + i;
      triangles[2 * (n_initial * j + i) + 1][2] = n_initial * (j + 1) + i + 1;
      nt++;
    }
  }

  double x0, x1, x2;
  double y0, y1, y2;
  double area;
  double xm, ym;
  double f0, f1, f2, fm;
  double i1, i2;
  int j, k, l, depth, nr;
  for(depth=0;depth<max_rec_depth;depth++){
    nr = 0;
    for(j=0;j<nt;j++){ //go through each triangle 
      x0 = points[triangles[j][0]][0];
      x1 = points[triangles[j][1]][0];
      x2 = points[triangles[j][2]][0];
      y0 = points[triangles[j][0]][1];
      y1 = points[triangles[j][1]][1];
      y2 = points[triangles[j][2]][1];
      f0 = points[triangles[j][0]][2];
      f1 = points[triangles[j][1]][2];
      f2 = points[triangles[j][2]][2];
      area = 0.5 * fabs(x0 * y1 + x1 * y2 + x2 * y0 - x0 * y2 - x1 * y0 - x2 * y1);
      xm = 1.0 / 3.0 * (x0 + x1 + x2);
      ym = 1.0 / 3.0 * (y0 + y1 + y2);
      fm = func(xm, ym);
      i1 = 0.0; //0.5 * (b - a) * (fa + fb);
      i2 = 1.0; //0.25 * (b - a) * (fa + 2.0 * fm + fb);
      if(fabs(i1-i2) > 1.0 * area * tol){ //this triangle needs refining
        points[np][0] = xm;
        points[np][1] = ym;
        points[np][2] = fm;
        to_refine[nr++] = triangles[j][0];
        to_refine[nr++] = triangles[j][1];
        to_refine[nr++] = triangles[j][2];
        to_refine[nr++] = np++;
      }
    }
    if(nr == 0){ break; }
    qsort(to_refine, nr, sizeof(int), intcmp); //sort and remove duplicate entries in the refined mesh
    for(j=0;j<nr;j++){
      for(k=j+1;k<nr;k++){
        if(to_refine[j] == to_refine[k]){
          for(l=k;k<nr-1;k++){
            to_refine[l] = to_refine[l + 1];
          }
          nr--; j--;
        }
      }
    }
    triangulate(points, to_refine, nr, triangles, &nt);
  }
  if(depth >= max_rec_depth){
    num_error |= ERROR_CONVERGENCE;
  }
}

void pyplot_tricountourf(plot3d_t plot, const char* args){

}
