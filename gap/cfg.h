#include <string>
#include "vec.h"

using namespace std;

extern int n;
extern int s_div;
extern int m;
extern int l;
extern int dim;
extern string potential_name;
extern string band_name;
extern double mu;
extern double k_max;

extern double t;
extern double tn;
extern double tnn;
extern double U;
extern double w_D;

void init_config(double &mu, double &U, double &t, double &tn, double &w_D, double new_mu, double new_U, double new_t, double new_tn, double new_w_D);
double epsilon(const Vec k);
double e_diff(const Vec k, const Vec q);
double vp_diff(const Vec k, const Vec q);
double vp(const Vec k);
double e_base_avg(const Vec k, const Vec q);
double e_surface(const Vec k, const Vec q);
double vp_surface(const Vec k, const Vec q);
double e_base(const Vec k, const Vec q);
double e_split(const Vec k, const Vec q);
double vp_split(const Vec k, const Vec q);
double V(const Vec k1, const Vec k2, double w, const double T, const unordered_map<double, vector<vector<vector<double>>>> &chi_cube);
int get_num_points_from_delta(double &delta);

extern double *weights[5];
extern double *points[5];
