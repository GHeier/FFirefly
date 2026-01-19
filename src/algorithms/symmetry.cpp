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
        // C4 symmetry: 0°, 90°, 180°, 270° rotations
        float rot = seed % 4;
        vector<vector<float>> R = rot_matrix(rot * M_PI / 2);
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
// k = -0.5 + i / (nx - 1) for 0-based indexing
vector<float> ind_to_vec(vector<int> inds, vector<int> &grid) {
    int dim = grid.size();
    vector<float> v(dim);
    for (int i = 0; i < dim; i++) {
        v[i] = -0.5 + (float)inds[i] / (grid[i] - 1);
    }
    return v;
}

// Scaled to be from -0.5 to 0.5, not -pi to pi
// i = (k + 0.5) * (nx - 1) for 0-based indexing
vector<int> vec_to_ind(vector<float> &v, vector<int> &grid) {
    int dim = grid.size();
    vector<int> inds(dim);
    for (int i = 0; i < dim; i++) {
        int idx = (int)round((v[i] + 0.5) * (grid[i] - 1));
        // Clamp to valid range [0, grid[i]-1]
        if (idx < 0) idx = 0;
        if (idx >= grid[i]) idx = grid[i] - 1;
        inds[i] = idx;
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
    int dim = grid.size();
    int i = idx % grid[0];

    if (dim == 1) {
        return {i};
    }

    int j = (int)round(idx / grid[0]) % grid[1];

    if (dim == 2) {
        return {i, j};
    }

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
    vector<string> syms = get_point_group_symmetries_from_lattice(lattice);
    int prod = 1;
    for (int x : grid) prod *= x;
    vector<int> mem_list(prod);
    int nx = 0, ny = 0, nz = 0;
    int mem_ind = 1;
    int dim = grid.size();

    for (int i = 0; i < mem_list.size(); i++) {
        vector<int> ind_set(dim);
        if (dim >= 1) ind_set[0] = nx;
        if (dim >= 2) ind_set[1] = ny;
        if (dim >= 3) ind_set[2] = nz;
        vector<float> v0 = ind_to_vec(ind_set, grid);

        // Apply all symmetries systematically
        for (int sym_idx = 0; sym_idx < syms.size(); sym_idx++) {
            vector<float> v1 = action_from_sym(syms[sym_idx], v0);
            vector<int> new_ind = vec_to_ind(v1, grid);
            int idx_nx = (dim >= 1) ? new_ind[0] : 0;
            int idx_ny = (dim >= 2) ? new_ind[1] : 0;
            int idx_nz = (dim >= 3) ? new_ind[2] : 0;
            int idx = get_global_ind(idx_nx, idx_ny, idx_nz, grid);

            // Bounds check
            if (idx < 0 || idx >= mem_list.size()) {
                printf("ERROR: Global index %d out of bounds [0, %d)\n", idx, (int)mem_list.size());
                continue;
            }

            if (mem_list[idx] == 0) {
                if (mem_list[i] == 0) {
                    mem_list[i] = mem_ind;
                    mem_ind++;
                }
                mem_list[idx] = mem_list[i];
            }
            else {
                if (mem_list[i] == 0) {
                    mem_list[i] = mem_list[idx];
                }
            }
        }

        // If this point wasn't mapped by any symmetry, create a new equivalence class
        if (mem_list[i] == 0) {
            mem_list[i] = mem_ind;
            mem_ind++;
        }

        update_inds(nx, ny, nz, grid);
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
