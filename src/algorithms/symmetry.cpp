#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;


vector<string> get_point_group_symmetries_from_lattice(string& lattice) {
    // Symmetry list
    //    "C" Rotation
    //    "R" Reflection
    //    "I" Inversion

    if (lattice == "SC") {
        return {"C4", "R", "I"};
    }
    else if (lattice == "BCC") {
        return {};
    }
    else if (lattice == "FCC") {
        return {};
    }
    return {}; // Default return
}

// Helper function for rotations
int get_quadrant(vector<float>& v0) {
    float x = v0[0];
    float y = v0[1];

    if (v0.size() == 2) {
        if (x >= 0 && y >= 0) return 1;
        if (x <= 0 && y >= 0) return 2;
        if (x <= 0 && y <= 0) return 3;
        if (x >= 0 && y <= 0) return 4;
    }
    // Dim = 3
    if (v0.size() > 3) {
        printf("Wrong vec size\n");
        exit(1);
    }
    float z = v0[2];
    if (x >= 0 && y >= 0 && z >= 0) return 1;
    if (x <= 0 && y >= 0 && z >= 0) return 2;
    if (x <= 0 && y <= 0 && z >= 0) return 3;
    if (x >= 0 && y <= 0 && z >= 0) return 4;
    if (x >= 0 && y >= 0 && z <= 0) return 5;
    if (x <= 0 && y >= 0 && z <= 0) return 6;
    if (x <= 0 && y <= 0 && z <= 0) return 7;
    if (x >= 0 && y <= 0 && z <= 0) return 8;
    return 0; // Default return
}

// Rotation matrix around z axis
vector<vector<float>> rot_matrix(float theta) {
    return {
        {cos(theta), -sin(theta), 0},
        {sin(theta), cos(theta), 0},
        {0, 0, 1}
    };
}

vector<float> mul(vector<vector<float>> &A, vector<float> &x) {
    vector<float> result = vector<float>(x.size());
    for (int i = 0; i < x.size(); i++) {
        for (int j = 0; j < x.size(); j++) {
            result[i] += A[i][j] * x[j];
        }
    }
    return result;
}

vector<float> action_from_sym(string& sym, vector<float> v0) {
    if (v0.size() > 3) {
        printf("Vector greater than 3D cannot have symmetries applied to it\n. Vec = ");
        for (float &x : v0) printf("%f ", x);
        exit(1);
    }
    int seed = rand();
    if (sym == "C4") { // 2D only for now
        float rot = seed % 4 + 1;
        vector<vector<float>> R = rot_matrix(rot * M_PI);
        return mul(R, v0);
    }
    else if (sym == "R") {
        int flip = seed % v0.size() + 1;
        if (flip == 1) v0[0] = -v0[0];
        if (flip == 2) v0[1] = -v0[1];
        if (flip == 3) v0[2] = -v0[2];
        return v0;
    }
    else if (sym == "I") {
        for (int i = 0; i < v0.size(); i++)
            v0[i] = -v0[i];
        return v0;
    }
    return v0; // Default return
}

// Scaled to be from -0.5 to 0.5, not -pi to pi
// k = -0.5 + (i - 1) / (nx - 1)
vector<float> ind_to_vec(vector<int> inds, vector<int> &grid) {
    vector<float> v(inds.size());
    for (int i = 0; i < inds.size(); i++) {
        v[i] = -0.5 + (float)(inds[i] - 1) / (grid[i] - 1);
    }
    return v;
}

// Scaled to be from -0.5 to 0.5, not -pi to pi
// i = 1 + (k + 0.5) * (nx - 1)
vector<int> vec_to_ind(vector<float> &v, vector<int> &grid) {
    vector<int> inds(v.size());
    for (int i = 0; i < v.size(); i++) {
        inds[i] = 1 + (v[i] + 0.5) * (grid[i] - 1);
    }
    return inds;
}

void update_inds(int &nx, int &ny, int &nz, vector<int> &grid) {
    if (grid.size() == 2) {
        ny++;
        if (ny >= grid[1]) {
            ny = 0;
            nx++;
        }
    }
    else if (grid.size() == 3) {
        nz++;
        if (nz >= grid[2]) {
            nz = 0;
            ny++;
            if (ny >= grid[1]) {
                ny = 0;
                nx++;
            }
        }
    }
}

vector<int> get_inds_from_global(int idx, vector<int> &grid) {
    int i = idx % grid[0];
    int j = (int)round(idx / grid[0]) % grid[1];
    int k = idx / (grid[0] * grid[1]);
    return {i, j, k};
}

int get_global_ind(int &nx, int &ny, int &nz, vector<int> &grid) {
    if (grid.size() == 3) {
        return nx + grid[0] * (ny + grid[1] * nz);
    }
    if (grid.size() == 2) {
        return nx + grid[0] * ny;
    }
    return 0; // Default return
}

// Grid mapping, constructs a list of equivalent points in index space
vector<int> sym_grid_map(vector<int> &grid, string& lattice) {
    int iter = 0;
    int max_iters = 10;
    vector<string> syms = get_point_group_symmetries_from_lattice(lattice);
    int prod = 1;
    for (int x : grid) prod *= x;
    vector<int> mem_list(prod);
    printf("Size of mem_list: %d\n", mem_list.size());
    int nx = 0, ny = 0, nz = 0;
    int mem_ind = 1;
    for (int i = 0; i < mem_list.size(); i++) {
        update_inds(nx, ny, nz, grid);
        vector<int> ind_set = {nx, ny, nz};
        vector<float> v0 = ind_to_vec(ind_set, grid);

        while (iter < max_iters) {
            int seed = rand() % syms.size();
            vector<float> v1 = action_from_sym(syms[seed], v0);
            vector<int> new_ind = vec_to_ind(v1, grid);
            int idx = get_global_ind(new_ind[0], new_ind[1], new_ind[2], grid);

            if (mem_list[idx] == 0) {
                if (mem_list[i] == 0) {
                    mem_list[i] = mem_ind;
                }
                mem_list[idx] = mem_list[i];
            }
            else {
                mem_list[i] = mem_list[idx];
            }
            iter++;
        }
        mem_ind++;
    }
    return mem_list;
}

vector<int> find(vector<int> &points, int p) {
    vector<int> inds;
    for (int i = 0; i < points.size(); i++) {
        if (points[i] == p) inds.push_back(i);
    }
    return inds;
}

// Returns all sets of indices that have smmetry-mapped points
// First index is the IBZ point
vector<vector<vector<int>>> get_reduced_grid(vector<int> &grid, string& lattice) {
    vector<vector<vector<int>>> reduced_grid;
    vector<int> equivalent_points = sym_grid_map(grid, lattice);
    printf("equivalent points size: %d\n", equivalent_points.size());
    vector<int> used_values;
    for (int i = equivalent_points.size() - 1; i >= 0; i--) {
        // Skip points mapped already
        vector<int> check_used = find(used_values, equivalent_points[i]);
        if (check_used.size() > 0) continue;

        // Find symmetry mapped points
        vector<int> idx_list = find(equivalent_points, equivalent_points[i]);
        vector<vector<int>> temp(idx_list.size());
        for (int j = 0; j < idx_list.size(); j++) {
            temp[j] = get_inds_from_global(idx_list[j], grid);
        }
        reduced_grid.push_back(temp);
        used_values.push_back(equivalent_points[i]);
    }
    return reduced_grid;
}
