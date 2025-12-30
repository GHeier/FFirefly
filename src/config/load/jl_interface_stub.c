#include "jl_interface.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Stub implementations for Julia interface (no julia.h dependency)
// These are deprecated - use run_julia_method2() from cpp_config.hpp instead

void load_julia() {
    fprintf(stderr, "WARNING: load_julia() is deprecated. Julia runtime is not embedded.\n");
    fprintf(stderr, "         Use run_julia_method2() instead for shell-based execution.\n");
}

bool call_julia_func(const char *folder, const char *filename, const char *module, const char *function) {
    fprintf(stderr, "ERROR: call_julia_func() is not available (julia.h not linked).\n");
    fprintf(stderr, "       Function: %s/%s::%s.%s\n", folder, filename, module, function);
    fprintf(stderr, "       Use run_julia_method2() from cpp_config.hpp instead.\n");
    exit(1);
    return false;
}

bool call_julia_func_bool(const char *folder, const char *filename, const char *module, const char *function) {
    return call_julia_func(folder, filename, module, function);
}

int call_julia_func_int(const char *folder, const char *filename, const char *module, const char *function) {
    fprintf(stderr, "ERROR: call_julia_func_int() is not available (julia.h not linked).\n");
    fprintf(stderr, "       Function: %s/%s::%s.%s\n", folder, filename, module, function);
    fprintf(stderr, "       Use run_julia_method2() from cpp_config.hpp instead.\n");
    exit(1);
    return 0;
}

float call_julia_func_float(const char *folder, const char *filename, const char *module, const char *function) {
    fprintf(stderr, "ERROR: call_julia_func_float() is not available (julia.h not linked).\n");
    fprintf(stderr, "       Function: %s/%s::%s.%s\n", folder, filename, module, function);
    fprintf(stderr, "       Use run_julia_method2() from cpp_config.hpp instead.\n");
    exit(1);
    return 0.0f;
}

double call_julia_func_double(const char *folder, const char *filename, const char *module, const char *function) {
    fprintf(stderr, "ERROR: call_julia_func_double() is not available (julia.h not linked).\n");
    fprintf(stderr, "       Function: %s/%s::%s.%s\n", folder, filename, module, function);
    fprintf(stderr, "       Use run_julia_method2() from cpp_config.hpp instead.\n");
    exit(1);
    return 0.0;
}

const char* call_julia_func_string(const char *folder, const char *filename, const char *module, const char *function) {
    fprintf(stderr, "ERROR: call_julia_func_string() is not available (julia.h not linked).\n");
    fprintf(stderr, "       Function: %s/%s::%s.%s\n", folder, filename, module, function);
    fprintf(stderr, "       Use run_julia_method2() from cpp_config.hpp instead.\n");
    exit(1);
    return strdup("");
}
