#include "hmatrix.hpp"
#include <cstdio>
#include <cstring>
#include <algorithm>

// Static wrapper for kernel callback (using :: to disambiguate from std::real)
field HMatrix::kernel_callback(const ::real* xx, const ::real* yy, void* data) {
    HMatrix* hmat = static_cast<HMatrix*>(data);
    return hmat->kernel_func(xx, yy);
}

// Helper function to fill H-matrix from kernel matrix (adapted from test_circle.c)
void HMatrix::fill_hmatrix_kernelmatrix(pckernelmatrix km, phmatrix h) {
    uint rsons, csons;
    uint i, j;
    pamatrix f;
    prkmatrix r;
    uint rows, cols;
    uint *ridx, *cidx;

    if (h->son) {
        // Subdivided matrix - recursively fill submatrices
        rsons = h->rsons;
        csons = h->csons;

        for (j = 0; j < csons; j++) {
            for (i = 0; i < rsons; i++) {
                fill_hmatrix_kernelmatrix(km, h->son[i + j * rsons]);
            }
        }
    }
    else if (h->f) {
        // Dense block - fill with kernel evaluations
        f = h->f;
        rows = h->rc->size;
        cols = h->cc->size;
        ridx = h->rc->idx;
        cidx = h->cc->idx;

        fillN_kernelmatrix(ridx, cidx, km, f);
    }
    else if (h->r) {
        // Low-rank block - use ACA approximation
        r = h->r;
        rows = h->rc->size;
        cols = h->cc->size;
        ridx = h->rc->idx;
        cidx = h->cc->idx;

        // Create temporary dense matrix
        pamatrix tmp = new_amatrix(rows, cols);
        fillN_kernelmatrix(ridx, cidx, km, tmp);

        // Truncate to low-rank using ACA
        decomp_fullaca_rkmatrix(tmp, 1e-8, NULL, NULL, r);

        del_amatrix(tmp);
    }
}

HMatrix::HMatrix(const vector<vector<double>>& points,
                 function<double(const double*, const double*)> kernel,
                 uint interpolation_order,
                 uint leafsize,
                 double eps,
                 double eta)
    : kernel_func(kernel), n_points(points.size()), m(interpolation_order), tolerance(eps)
{
    if (points.empty()) {
        throw invalid_argument("Points vector cannot be empty");
    }

    dim = points[0].size();
    if (dim == 0) {
        throw invalid_argument("Point dimension cannot be zero");
    }

    // Initialize H2Lib (safe to call multiple times)
    int argc = 0;
    char** argv = nullptr;
    init_h2lib(&argc, &argv);

    // Create kernel matrix object
    km = new_kernelmatrix(dim, n_points, m);
    km->kernel = kernel_callback;
    km->data = this;  // Pass 'this' pointer for callback

    // Copy points to kernel matrix structure
    for (uint i = 0; i < n_points; i++) {
        if (points[i].size() != dim) {
            throw invalid_argument("All points must have the same dimension");
        }
        for (uint d = 0; d < dim; d++) {
            km->x[i][d] = points[i][d];
        }
    }

    // Create cluster geometry
    cg = creategeometry_kernelmatrix(km);

    // Create cluster tree
    uint* idx = (uint*)allocmem(sizeof(uint) * n_points);
    for (uint i = 0; i < n_points; i++) {
        idx[i] = i;
    }
    root = build_adaptive_cluster(cg, n_points, idx, leafsize);

    // Create block tree with admissibility condition
    broot = build_strict_block(root, root, &eta, admissible_2_cluster);

    // Create H-matrix structure
    hm = build_from_block_hmatrix(broot, m * m);

    // Fill H-matrix with kernel evaluations
    fill_hmatrix_kernelmatrix(km, hm);
}

HMatrix::~HMatrix() {
    // Clean up H2Lib structures in reverse order of creation
    if (hm) del_hmatrix(hm);
    if (broot) del_block(broot);
    if (root) del_cluster(root);
    if (cg) del_clustergeometry(cg);
    if (km) del_kernelmatrix(km);

    // Note: We don't call uninit_h2lib() as it's a global cleanup
    // and other HMatrix objects might still be in use
}

vector<double> HMatrix::matvec(const vector<double>& x) const {
    if (x.size() != n_points) {
        throw invalid_argument("Input vector size does not match matrix dimension");
    }

    vector<double> y(n_points, 0.0);
    matvec(x, y, 1.0, false);
    return y;
}

void HMatrix::matvec(const vector<double>& x, vector<double>& y,
                     double alpha, bool transpose) const {
    if (x.size() != n_points) {
        throw invalid_argument("Input vector size does not match matrix dimension");
    }
    if (y.size() != n_points) {
        y.resize(n_points, 0.0);
    }

    // Create H2Lib vectors
    pavector xvec = new_avector(n_points);
    pavector yvec = new_avector(n_points);

    // Copy input data
    for (uint i = 0; i < n_points; i++) {
        xvec->v[i] = x[i];
        yvec->v[i] = y[i];
    }

    // Perform matrix-vector multiplication
    mvm_hmatrix_avector(alpha, transpose, hm, xvec, yvec);

    // Copy result back
    for (uint i = 0; i < n_points; i++) {
        y[i] = yvec->v[i];
    }

    // Clean up
    del_avector(yvec);
    del_avector(xvec);
}

vector<vector<double>> HMatrix::to_dense() const {
    vector<vector<double>> dense(n_points, vector<double>(n_points, 0.0));

    // Create unit vectors and multiply
    for (uint j = 0; j < n_points; j++) {
        vector<double> ej(n_points, 0.0);
        ej[j] = 1.0;
        vector<double> col = matvec(ej);
        for (uint i = 0; i < n_points; i++) {
            dense[i][j] = col[i];
        }
    }

    return dense;
}

size_t HMatrix::memory_size() const {
    return getsize_hmatrix(hm);
}

size_t HMatrix::nearfield_size() const {
    return getnearsize_hmatrix(hm);
}

size_t HMatrix::farfield_size() const {
    return getfarsize_hmatrix(hm);
}

double HMatrix::compression_ratio() const {
    size_t dense_size = n_points * n_points * sizeof(double);
    size_t hmat_size = memory_size();
    return (double)dense_size / (double)hmat_size;
}

double HMatrix::norm() const {
    return norm2_hmatrix(hm);
}

void HMatrix::print_stats() const {
    printf("\n========================================\n");
    printf("H-Matrix Statistics\n");
    printf("========================================\n");
    printf("Matrix size:        %u x %u\n", n_points, n_points);
    printf("Spatial dimension:  %u\n", dim);
    printf("Interpolation order: %u\n", m);
    printf("Cluster depth:      %u\n", getdepth_cluster(root));
    printf("Block depth:        %u\n", getdepth_block(broot));
    printf("Number of clusters: %u\n", root->desc);
    printf("Number of blocks:   %u\n", broot->desc);
    printf("\nMemory usage:\n");
    printf("  Total:       %.2f MB (%.2f KB/DoF)\n",
           memory_size() / 1048576.0,
           memory_size() / 1024.0 / n_points);
    printf("  Nearfield:   %.2f MB\n", nearfield_size() / 1048576.0);
    printf("  Farfield:    %.2f MB\n", farfield_size() / 1048576.0);
    printf("  Dense equiv: %.2f MB\n",
           (n_points * n_points * sizeof(double)) / 1048576.0);
    printf("  Compression: %.2fx\n", compression_ratio());
    printf("========================================\n\n");
}
