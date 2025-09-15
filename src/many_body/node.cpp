#include "../config/load/cpp_config.hpp"
#include "../config/load/jl_interface.h"
#include "node.hpp"

extern "C" void many_body_wrapper() {
    many_body_loop();
}

void many_body_loop() {
    string folder = "many_body/";
    string filename = "many_body_loop";
    string module = "ManyBodyLoop";
    string function = "main";
    call_julia_func(folder.c_str(), filename.c_str(), module.c_str(), function.c_str());
}

