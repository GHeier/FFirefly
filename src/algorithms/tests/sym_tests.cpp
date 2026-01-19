#include "src/algorithms/symmetry.hpp"
#include <vector>

#include "sym_tests.hpp"

void print_vector(const vector<vector<vector<int>>> &vec) {
    for (int i = 0; i < vec.size(); i++) {
        cout << "Set " << i << ":" << endl;
        for (int j = 0; j < vec[i].size(); j++) {
            cout << "[";
            for (int k = 0; k < vec[i][j].size(); k++) {
                cout << vec[i][j][k];
                if (k < vec[i][j].size() - 1) cout << ", ";
            }
            cout << "]" << endl;
        }
    }
}

bool even_small() {
    vector<int> grid = {2, 2};
    string lattice = "SC";
    vector<vector<vector<int>>> reduced = get_reduced_grid(grid, lattice);

    print_3D_vector(reduced);

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

