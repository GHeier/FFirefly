#include <vector>
#include <string>
#include <unordered_map>
<<<<<<< HEAD
=======
#include <vector>
>>>>>>> origin/main
#include "vec.h"

using namespace std;

extern int n;
extern int s_div;
extern int s_pts;
extern int m;
extern int l;
extern int dim;
extern string potential_name;
extern string band_name;
extern float mu;
extern float k_max;

extern float t;
extern float tn;
extern float tnn;
extern float U;
extern float wc;

<<<<<<< HEAD
void init_config(float &mu, float &U, float &t, float &tn, float &w_D, float new_mu, float new_U, float new_t, float new_tn, float new_w_D);
void change_global_constant(float &a, float b);
float epsilon(const Vec k);
float e_diff(const Vec k, const Vec q);
float e_base_avg(const Vec k, const Vec q);
float vp_diff(const Vec k, const Vec q);
float vp(const Vec k);
float V(const Vec k1, const Vec k2, float w, const float T, const unordered_map<float, vector<vector<vector<float>>>> &chi_cube);
int get_num_points_from_delta(float &delta);
=======
extern float sin_arr[31416]; 
extern float cos_arr[31416]; 

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
>>>>>>> origin/main

extern float *weights[5];
extern float *points[5];
