#include <string>
#include "vec.h"

using namespace std;

extern int n;
extern int m;
extern int dim;
extern string potential_name;
extern string band_name;
extern double mu;
extern double k_max;

extern double t;
extern double tn;
extern double U;
extern double w_D;

void init_config(double &mu, double &U, double &t, double &tn, double new_mu, double new_U, double new_t, double new_tn);
double epsilon(const Vec k);
double vp(const Vec k);
double V(const Vec k1, const Vec k2, const double T, const vector<vector<vector<double>>> &chi_cube);
