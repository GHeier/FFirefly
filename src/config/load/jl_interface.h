#pragma once
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif
void load_julia();

bool call_julia_func(const char *folder, const char *filename, const char* module, const char *function);
bool call_julia_func_bool(const char *folder, const char *filename, const char* module, const char *function);
int call_julia_func_int(const char *folder, const char *filename, const char* module, const char *function);
float call_julia_func_float(const char *folder, const char *filename, const char* module, const char *function);
double call_julia_func_double(const char *folder, const char *filename, const char* module, const char *function);
const char* call_julia_func_string(const char *folder, const char *filename, const char* module, const char *function);

#ifdef __cplusplus
}
#endif


