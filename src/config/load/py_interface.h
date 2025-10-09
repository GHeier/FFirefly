#pragma once

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

void start_python();
void end_python();
void call_python_func(const char *folder, const char *filename, const char *function);
bool call_python_func_bool(const char *folder, const char *filename, const char *function);
int call_python_func_int(const char *folder, const char *filename, const char *function);
float call_python_func_float(const char *folder, const char *filename, const char *function);
double call_python_func_double(const char *folder, const char *filename, const char *function);
const char* call_python_func_string(const char *folder, const char *filename, const char *function);

#ifdef __cplusplus
}
#endif

