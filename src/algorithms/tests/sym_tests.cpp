#include "src/algorithms/symmetry.hpp"
#include "src/config/load/c_config.h"
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include "sym_tests.hpp"

using namespace std;

float test_func(vector<float> vec) {
    float val = 0.0;
    for (float v : vec) {
        val += cos(v);
    }
    return val;
}

bool test_vector(const vector<vector<vector<int>>> &vec, vector<int> &grid) {
    for (int i = 0; i < vec.size(); i++) {
        float prev_v = 0.0, curr_v = 0.0;
        for (int j = 0; j < vec[i].size(); j++) {
            vector<float> v = ind_to_vec(vec[i][j], grid);
            curr_v = test_func(v);
            if (j > 0 && fabs(curr_v - prev_v) > 1e-5) {
                return false;
            }
            prev_v = curr_v;
        }
    }
    return true;
}

void print_vector(const vector<vector<vector<int>>> &vec, vector<int> &grid) {
    printf("Vector size: %d\n", vec.size());
    for (int i = 0; i < vec.size(); i++) {
        cout << "Set " << i << ":" << endl;
        for (int j = 0; j < vec[i].size(); j++) {
            cout << "[";
            vector<float> v = ind_to_vec(vec[i][j], grid);
            for (int k = 0; k < vec[i][j].size(); k++) {
                cout << vec[i][j][k];
                if (k < vec[i][j].size() - 1) cout << ", ";
            }
            cout << "] ";
            for (float val : v) {
                cout << val << " ";
            }
            cout << test_func(v) << endl;
        }
    }
}

bool even_small_2d() {
    vector<int> grid = {6, 6};
    string lattice = "SC";
    vector<vector<vector<int>>> reduced = get_reduced_grid(grid, lattice);
    //print_vector(reduced, grid);
    return test_vector(reduced, grid);
}

bool odd_small_2d() {
    vector<int> grid = {5, 5};
    string lattice = "SC";
    vector<vector<vector<int>>> reduced = get_reduced_grid(grid, lattice);
    //print_vector(reduced, grid);
    return test_vector(reduced, grid);
}

bool odd_big_2d() {
    vector<int> grid = {300, 300};
    string lattice = "SC";
    vector<vector<vector<int>>> reduced = get_reduced_grid(grid, lattice);
    //print_vector(reduced, grid);
    return test_vector(reduced, grid);
}

bool sym_tests() {
  int num_tests = 3;
  bool all_tests[num_tests] = {
      odd_small_2d(), 
      even_small_2d(), 
      odd_big_2d(), 
  };
  return print_test_results(all_tests, num_tests, "Symmetry tests");
}

