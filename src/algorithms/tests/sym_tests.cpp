#include "src/algorithms/symmetry.hpp"
#include "src/config/load/c_config.h"
#include <vector>
#include <string>
#include <iostream>

#include "sym_tests.hpp"

using namespace std;

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
            cout << endl;
        }
    }
}

bool even_small() {
    vector<int> grid = {5, 5};
    string lattice = "SC";
    vector<vector<vector<int>>> reduced = get_reduced_grid(grid, lattice);

    print_vector(reduced, grid);

    vector<vector<vector<int>>> expected = {
        {
            {0, 0, 0},
            {0, 0, 1},
        },
        {
            {0, 1, 0},
            {0, 1, 1},
        },
    };

    if (reduced != expected) {
        return false;
    }
    return true;
}

bool sym_tests() {
  int num_tests = 1;
  bool all_tests[num_tests] = {
      even_small(), 
  };
  return print_test_results(all_tests, num_tests, "Symmetry tests");
}

