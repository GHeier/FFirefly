#include <string>
#include <unordered_map>
#include <vector>
#include "vec.h"

using namespace std;

extern int n;
extern int m;
extern int l;
extern int dim;
extern string potential_name;
extern string band_name;
extern double mu;
extern double k_max;

extern double t;
extern double tn;
extern double U;
extern double w_D;

extern float sin_arr[31416]; 
extern float cos_arr[31416]; 

void init_config(double &mu, double &U, double &t, double &tn, double &w_D, double new_mu, double new_U, double new_t, double new_tn, double new_w_D);
double epsilon(const Vec k);
double vp(const Vec k);
double V(const Vec k1, const Vec k2, double w, const double T, const unordered_map<double, vector<vector<vector<double>>>> &chi_cube);

extern double *weights[5];
extern double *points[5];

double get_sin(double x);
double get_cos(double x);
