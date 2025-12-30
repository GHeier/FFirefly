#include "surface_tests.hpp"
#include "src/config/load/cpp_config.hpp"
#include "src/hamiltonian/models/band_structure.hpp"
#include "src/objects/vec.hpp"
#include "src/objects/surfaces.hpp"

using namespace std;

bool surface_test_3D() {
    set_global(nbnd, 1);
    set_global(dimension, 3);
    eff_mass = 0.5;
    set_global(k_mesh, {34, 34, 34});

    float s_val = 1.0;

    auto func = [](Vec k) { return epsilon(1, k); };
    Vec k_temp(1.0, 1.0, 1.0);
    Surface surf = tetrahedron_surface(func, s_val);
    Vec temp = surf.faces[0];
    float E = epsilon(1, temp);
    return fabs(E - s_val) < 1e-2;
}

bool surface_tests() { return surface_test_3D(); }
