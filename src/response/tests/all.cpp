#include <string>
#include "../../config/load/c_config.h"
#include "../../config/load/cpp_config.hpp"
#include "../../config/load/jl_interface.h"
#include "../../objects/CMField/fields.hpp"
#include "../../objects/vec.hpp"
#include "all.hpp"

using namespace std;

extern "C" bool response_tests() {
    printf("Running Response Tests\n");
    int num_tests = 1;
    bool all_tests[num_tests] = {
        sparse_ir_response_test_3D(),
    };
    return print_test_results(all_tests, num_tests, "Response Tests");
}

bool sparse_ir_response_test_3D() {
    set_global(k_mesh, {40, 40, 40});
    set_global(Temperature, 0.25);
    set_global(dimension, 3);
    set_global(fermi_energy, 0);
    set_global(dynamic, false);
    set_global(nbnd, 1);
    band.push_back("tight_binding");
    t0.push_back(1.0);

    string folder = "response/";
    string file = "sparse_ir_response";
    string module = "response_ir";
    string function = "get_ckio_ir";
    call_julia_func(folder.c_str(), file.c_str(), module.c_str(), function.c_str());
    Field_C chi(outdir + prefix + "_chi.dat");
    Vec q(3.14, 3.14, 3.14);
    if (abs(real(chi(q)) - 0.40) < 1e-4)
        return true;
    return false;
}

