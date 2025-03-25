#pragma once
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

bool call_julia_func(const char *folder, const char *filename, const char* module, const char *function);

#ifdef __cplusplus
}
#endif


