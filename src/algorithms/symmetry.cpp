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

// Apply periodic boundary conditions to wrap point into [-0.5, 0.5)
vector<float> apply_periodic_bc(vector<float> v) {
    for (int i = 0; i < v.size(); i++) {
        // Wrap to [-0.5, 0.5) using modulo arithmetic
        v[i] = fmod(v[i] + 0.5, 1.0) - 0.5;
        if (v[i] < -0.5) v[i] += 1.0;
        if (v[i] >= 0.5) v[i] -= 1.0;
    }
    return v;
}

// Generate ALL symmetry-equivalent points for a given symmetry operation
vector<vector<float>> generate_equivalent_points(string& sym, vector<float> v0) {
    if (v0.size() > 3) {
        printf("Vector greater than 3D cannot have symmetries applied to it\n. Vec = ");
        for (float &x : v0) printf("%f ", x);
        exit(1);
    }

    vector<vector<float>> equivalent_points;

    if (sym == "C4") {
        // C4 symmetry: 0°, 90°, 180°, 270° rotations around z-axis
        for (int i = 0; i < 4; i++) {
            vector<vector<float>> R = rot_matrix(i * M_PI / 2);
            vector<float> rotated = mul(R, v0);
            equivalent_points.push_back(apply_periodic_bc(rotated));
        }
    }
    else if (sym == "R") {
        // Mirror reflections across planes
        int dim = v0.size();

        // Mirror across x=0 plane
        vector<float> v1 = v0;
        v1[0] = -v1[0];
        equivalent_points.push_back(apply_periodic_bc(v1));

        // Mirror across y=0 plane
        if (dim >= 2) {
            vector<float> v2 = v0;
            v2[1] = -v2[1];
            equivalent_points.push_back(apply_periodic_bc(v2));
        }

        // Mirror across z=0 plane (for 3D)
        if (dim >= 3) {
            vector<float> v3 = v0;
            v3[2] = -v3[2];
            equivalent_points.push_back(apply_periodic_bc(v3));
        }

        // For SC lattice, also include diagonal mirrors (x=y, x=-y)
        if (dim >= 2) {
            vector<float> v4 = {v0[1], v0[0]};
            if (dim == 3) v4.push_back(v0[2]);
            equivalent_points.push_back(apply_periodic_bc(v4));

            vector<float> v5 = {-v0[1], -v0[0]};
            if (dim == 3) v5.push_back(v0[2]);
            equivalent_points.push_back(apply_periodic_bc(v5));
        }
    }
    else if (sym == "I") {
        // Inversion: (x, y, z) → (-x, -y, -z)
        vector<float> inverted = v0;
        for (int i = 0; i < v0.size(); i++)
            inverted[i] = -inverted[i];
        equivalent_points.push_back(apply_periodic_bc(inverted));
    }

    return equivalent_points;
}

// Scaled to be from -0.5 to 0.5 (half-open, [-0.5, 0.5)), not -pi to pi.
// Matches the mesh convention used by get_fractional_mesh/get_kmesh
// (response_utils.jl): k = -0.5 + i / nx for 0-based indexing.
vector<float> ind_to_vec(vector<int> inds, vector<int> &grid) {
    int dim = grid.size();
    vector<float> v(dim);
    for (int i = 0; i < dim; i++) {
        v[i] = -0.5 + (float)inds[i] / grid[i];
    }
    return v;
}

// Inverse of ind_to_vec: i = (k + 0.5) * nx for 0-based indexing.
// Grid is periodic, so wrap out-of-range indices rather than clamping.
vector<int> vec_to_ind(vector<float> &v, vector<int> &grid) {
    int dim = grid.size();
    vector<int> inds(dim);
    for (int i = 0; i < dim; i++) {
        int idx = (int)round((v[i] + 0.5) * grid[i]);
        idx = ((idx % grid[i]) + grid[i]) % grid[i];  // periodic wrap into [0, grid[i]-1]
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

// Check if two points are equal within numerical tolerance
bool points_equal(vector<float>& v1, vector<float>& v2, float tol = 1e-6) {
    if (v1.size() != v2.size()) return false;
    for (int i = 0; i < v1.size(); i++) {
        if (fabs(v1[i] - v2[i]) > tol) return false;
    }
    return true;
}

// Grid mapping, constructs a list of equivalent points in index space
vector<int> sym_grid_map(vector<int> &grid, string& lattice) {
    vector<string> syms = get_point_group_symmetries_from_lattice(lattice);
    int prod = 1;
    for (int x : grid) prod *= x;
    vector<int> mem_list(prod, 0);  // Initialize all to 0 (unassigned)
    int mem_ind = 1;
    int dim = grid.size();

    // Process each point in the grid
    for (int i = 0; i < prod; i++) {
        // Skip if already assigned to an equivalence class
        if (mem_list[i] != 0) continue;

        // Get the current point's indices and coordinates
        vector<int> current_inds = get_inds_from_global(i, grid);
        vector<float> current_vec = ind_to_vec(current_inds, grid);

        // Start a new equivalence class
        mem_list[i] = mem_ind;

        // Collect all equivalent points by applying all symmetries
        vector<vector<float>> all_equivalent_points;
        all_equivalent_points.push_back(current_vec);  // Include the original point

        // Generate equivalent points from all symmetry operations
        for (int sym_idx = 0; sym_idx < syms.size(); sym_idx++) {
            vector<vector<float>> sym_points = generate_equivalent_points(syms[sym_idx], current_vec);

            for (auto& equiv_vec : sym_points) {
                // Check if this point is already in our list
                bool already_added = false;
                for (auto& existing : all_equivalent_points) {
                    if (points_equal(equiv_vec, existing)) {
                        already_added = true;
                        break;
                    }
                }
                if (!already_added) {
                    all_equivalent_points.push_back(equiv_vec);
                }
            }
        }

        // Assign all equivalent points to the same equivalence class
        for (auto& equiv_vec : all_equivalent_points) {
            vector<int> equiv_inds = vec_to_ind(equiv_vec, grid);

            // Validate indices are in bounds
            bool valid = true;
            for (int d = 0; d < dim; d++) {
                if (equiv_inds[d] < 0 || equiv_inds[d] >= grid[d]) {
                    valid = false;
                    break;
                }
            }
            if (!valid) continue;

            int idx_nx = (dim >= 1) ? equiv_inds[0] : 0;
            int idx_ny = (dim >= 2) ? equiv_inds[1] : 0;
            int idx_nz = (dim >= 3) ? equiv_inds[2] : 0;
            int idx = get_global_ind(idx_nx, idx_ny, idx_nz, grid);

            // Bounds check
            if (idx < 0 || idx >= prod) {
                continue;
            }

            // Assign to equivalence class if not already assigned
            if (mem_list[idx] == 0) {
                mem_list[idx] = mem_ind;
            }
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
