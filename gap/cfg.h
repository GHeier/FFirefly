#include <vector>
#include <string>
#include <unordered_map>
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
extern int num_eigenvalues_to_save;
extern float mu;
extern float k_max;

extern float t;
extern float tn;
extern float tnn;
extern float U;
extern float wc;

void init_config(float &mu, float &U, float &t, float &tn, float &w_D, float new_mu, float new_U, float new_t, float new_tn, float new_w_D);
void change_global_constant(float &a, float b);
float epsilon(const Vec k);
float e_diff(const Vec k, const Vec q);
float e_base_avg(const Vec k, const Vec q);
float vp_diff(const Vec k, const Vec q);
float vp(const Vec k);
float V(const Vec k1, const Vec k2, float w, const float T, const unordered_map<float, vector<vector<vector<float>>>> &chi_cube);
int get_num_points_from_delta(float &delta);

extern float *weights[5];
extern float *points[5];
