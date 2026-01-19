#pragma once

#include <vector>
#include <string>

using namespace std;

vector<float> ind_to_vec(vector<int> inds, vector<int> &grid);
vector<vector<vector<int>>> get_reduced_grid(vector<int> &grid, string& lattice);

