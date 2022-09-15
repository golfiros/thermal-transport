#define TOL 0.00001 //absolute numerical tolerance for calculations
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

//root finding by bracketing
double find_root(double a, double b, double (*func)(double));

typedef struct plot_s *plot_t;

void pyplot_import(const char* filename);

void pyplot_close();

void pyplot_figure(int figure);

void pyplot_savefig(const char* filename);

void pyplot_plot(plot_t plot, const char* color, const char* label);

void pyplot_xlabel(const char* label);

void pyplot_ylabel(const char* label);

void pyplot_legend(const char* title);

plot_t create_plot();

void delete_plot(plot_t plot);

void sample_plot(plot_t plot, double a, double b, double (*func)(double));

int get_plot_points(plot_t plot);

double* get_plot_x(plot_t plot);

void set_plot_x(plot_t plot, const double *x, int n_points);

double* get_plot_y(plot_t plot);

void set_plot_y(plot_t plot, const double *y);

