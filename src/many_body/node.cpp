#include "../config/load/cpp_config.hpp"
#include "../config/load/jl_interface.h"
#include "vertex.hpp"
#include "self_energy.hpp"
#include "renormalization.hpp"
#include "node.hpp"

extern "C" void many_body_wrapper() {
    if (calculation == "vertex")
        vertex_wrapper();
    else if (calculation == "self_energy")
        self_energy_wrapper();
    else if (calculation == "renormalization")
        renormalization_wrapper();
    else
        many_body_loop();
}

void many_body_loop() {
    string folder = "many_body/";
    string filename = "many_body_loop";
    string module = "ManyBodyLoop";
    string function = "main";
    call_julia_func(folder.c_str(), filename.c_str(), module.c_str(), function.c_str());
}

void renormalization_wrapper() {
    if (method == "analytic") {
        if (interaction == "FLEX")
            FLEX_renormalization();
        else
            printf("Analytic renormalization for '%s' interaction not available\n", interaction.c_str());
    }
    else 
        self_energy_renormalization();
}

void self_energy_wrapper() {
    if (method == "sparse_ir")
        call_self_energy();
    else
        printf("Method '%s' not available\n", method.c_str());
}

