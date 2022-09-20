#define TOL 0.0001 //absolute numerical tolerance for calculations
#define MAX_REC 2048 //max recursions for methods
#define N_INITIAL 16 //initial even subdivisions

//simple error handling with global variable

enum {
  ERROR_DIV_BY_ZERO = 1 << 0,
  ERROR_INVALID = 1 << 1,
  ERROR_CONVERGENCE = 1 << 2
};

void reset_error_numerics();
int get_error_numerics();

//numerical integration by adaptive quadrature
double integrate(double a, double b, double (*func)(double));

//root finding by bisection 
double find_root(double a, double b, double (*func)(double));

//opens and closes python file for pyplot
void pyplot_import(const char* filename);

void pyplot_close();

//pyplot figure manipulation
void pyplot_figure(int figure);

void pyplot_savefig(const char* filename);

void pyplot_xlabel(const char* label);

void pyplot_ylabel(const char* label);

void pyplot_legend(const char* title);

//pyplot 2d plotting
typedef struct plot2d_s *plot2d_t;

plot2d_t create_plot2d();

void delete_plot2d(plot2d_t plot);

void sample_plot2d(plot2d_t plot, double a, double b, double (*func)(double));

int get_plot2d_points(plot2d_t plot);

double* get_plot2d_x(plot2d_t plot);

void set_plot2d_x(plot2d_t plot, const double *x, int n_points);

double* get_plot2d_y(plot2d_t plot);

void set_plot2d_y(plot2d_t plot, const double *y);

void pyplot_plot(plot2d_t plot, const char* args);

//pyplot 3d plotting
typedef struct plot3d_s *plot3d_t;
