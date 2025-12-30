#include "py_interface.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Stub implementations for Python interface (no Python.h dependency)
// These are deprecated - use run_python_method2() from cpp_config.hpp instead

void start_python() {
    fprintf(stderr, "WARNING: start_python() is deprecated. Python interpreter is not embedded.\n");
    fprintf(stderr, "         Use run_python_method2() instead for shell-based execution.\n");
}

void end_python() {
    // No-op
}

void call_python_func(const char *folder, const char *filename, const char *function) {
    fprintf(stderr, "ERROR: call_python_func() is not available (Python.h not linked).\n");
    fprintf(stderr, "       Function: %s/%s::%s\n", folder, filename, function);
    fprintf(stderr, "       Use run_python_method2() from cpp_config.hpp instead.\n");
    exit(1);
}

bool call_python_func_bool(const char *folder, const char *filename, const char *function) {
    fprintf(stderr, "ERROR: call_python_func_bool() is not available (Python.h not linked).\n");
    fprintf(stderr, "       Function: %s/%s::%s\n", folder, filename, function);
    fprintf(stderr, "       Use run_python_method2() from cpp_config.hpp instead.\n");
    exit(1);
    return false;
}

int call_python_func_int(const char *folder, const char *filename, const char *function) {
    fprintf(stderr, "ERROR: call_python_func_int() is not available (Python.h not linked).\n");
    fprintf(stderr, "       Function: %s/%s::%s\n", folder, filename, function);
    fprintf(stderr, "       Use run_python_method2() from cpp_config.hpp instead.\n");
    exit(1);
    return 0;
}

float call_python_func_float(const char *folder, const char *filename, const char *function) {
    fprintf(stderr, "ERROR: call_python_func_float() is not available (Python.h not linked).\n");
    fprintf(stderr, "       Function: %s/%s::%s\n", folder, filename, function);
    fprintf(stderr, "       Use run_python_method2() from cpp_config.hpp instead.\n");
    exit(1);
    return 0.0f;
}

double call_python_func_double(const char *folder, const char *filename, const char *function) {
    fprintf(stderr, "ERROR: call_python_func_double() is not available (Python.h not linked).\n");
    fprintf(stderr, "       Function: %s/%s::%s\n", folder, filename, function);
    fprintf(stderr, "       Use run_python_method2() from cpp_config.hpp instead.\n");
    exit(1);
    return 0.0;
}

const char* call_python_func_string(const char *folder, const char *filename, const char *function) {
    fprintf(stderr, "ERROR: call_python_func_string() is not available (Python.h not linked).\n");
    fprintf(stderr, "       Function: %s/%s::%s\n", folder, filename, function);
    fprintf(stderr, "       Use run_python_method2() from cpp_config.hpp instead.\n");
    exit(1);
    return strdup("");
}
