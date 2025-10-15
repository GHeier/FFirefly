#include "../config/load/py_interface.h"
#include "../config/load/jl_interface.h"
#include "../config/load/cpp_config.hpp"
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
    else if (calculation == "response")
        response_wrapper();
    else
        if (method == "triqs")
            triqs_loop();
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

void triqs_loop() {
    string folder = "many_body/";
    string filename = "many_body_triqs";
    string function = "main";
    call_python_func(folder.c_str(), filename.c_str(), function.c_str());
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


extern "C" void polarization_wrapper();

void response_wrapper() {
    if (method == "libtetrabz")
        polarization_wrapper();
    else if (method == "sparse_ir")
        ir_wrapper();
    else
        printf("Method `%s` not found\n", method.c_str());
}

void ir_wrapper() {
    string folder = "response/";
    string filename = "sparse_ir_response";
    string module = "response_ir";
    string function = "get_ckio_ir";
    call_julia_func(folder.c_str(), filename.c_str(), module.c_str(), function.c_str());
}

void vertex_wrapper() {
  if (interaction == "FLEX") {
    call_flex();
  }
  else 
      printf("Interaction '%s' not available\n", interaction.c_str());
}
