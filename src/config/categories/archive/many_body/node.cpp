#include "node.hpp"
#include "../config/load/cpp_config.hpp"
#include "../config/load/py_interface.h"
#include "../config/load/jl_interface.h"
#include "renormalization.hpp"
#include "response_tetrabz.hpp"
#include "self_energy.hpp"
#include "vertex.hpp"
#include <iostream>

using namespace std;

void loop_impl() {
    string folder = "many_body/";
    string filename = "many_body_loop";
    string module = "ManyBodyLoop";
    string function = "main";
    call_julia_func(folder.c_str(), filename.c_str(), module.c_str(), function.c_str());
}

void response_sparse_ir_impl() {
    string folder = "many_body/";
    string filename = "response_ir";
    string module = "ResponseIr";
    string function = "main";
    call_julia_func(folder.c_str(), filename.c_str(), module.c_str(), function.c_str());
}

void triqs_impl() {
    string folder = "many_body";
    string filename = "many_body_triqs";
    string function = "main";
    call_python_func(folder.c_str(), filename.c_str(), function.c_str());
}

/**
 * Wrapper function for many_body category
 * Dispatches to appropriate calculation/method based on config variables
 */
extern "C" void many_body_wrapper() {
    printv("Running many_body_wrapper\n");
    if (calculation == "loop") {
        loop_impl();
    }
    else if (calculation == "renormalization") {
        renormalization();
    }
    else if (calculation == "response") {
        if (method == "libtetrabz") {
            response_libtetrabz();
        }
        else if (method == "sparse_ir") {
            response_sparse_ir_impl();
        }
        else
            cout << "method \"" << method << "\" not recognized for calculation response" << endl;
    }
    else if (calculation == "self_energy") {
        if (method == "sparse_ir") {
            self_energy_sparse_ir();
        }
        else
            cout << "method \"" << method << "\" not recognized for calculation self_energy" << endl;
    }
    else if (calculation == "triqs") {
        triqs_impl();
    }
    else if (calculation == "vertex") {
        if (method == "FLEX") {
            vertex_FLEX();
        }
        else
            cout << "method \"" << method << "\" not recognized for calculation vertex" << endl;
    }
    else
        cout << "calculation \"" << calculation << "\" not recognized for category many_body" << endl;
}
