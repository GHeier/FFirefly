#include <jlcxx/jlcxx.hpp>
#include <jlcxx/array.hpp>
#include <jlcxx/stl.hpp>
#include "config/load/cpp_config.hpp"
#include "config/load/c_config.h"
#include "response/susceptibility.hpp"
#include "hamiltonian/interaction.hpp"
#include "hamiltonian/band_structure.hpp"
#include "objects/vec.hpp"
#include "objects/fastfield.hpp"

JLCXX_MODULE define_julia_module(jlcxx::Module& mod) {
    // Vec class
    mod.add_type<Vec>("Vec")
        .constructor<float, float, float, float, float, int, int>()
        .method("norm", &Vec::norm)
        .method("x", [](const Vec& v) { return v.x; })
        .method("y", [](const Vec& v) { return v.y; })
        .method("z", [](const Vec& v) { return v.z; })
        .method("area", [](const Vec& v) { return v.area; })
        .method("dimension", [](const Vec& v) { return v.dimension; });

    // Epsilon function
    mod.method("epsilon", [](int n, double kx, double ky, double kz) {
        Vec k_vec(kx, ky, kz);
        return epsilon(n, k_vec);
    });

    // Load C config
    mod.method("load_c_config", [](const std::string& path) {
        read_c_config_wrapper(path);
        load_cpp_config();
    });

    // Add more bindings as needed...
}
