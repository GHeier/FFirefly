#pragma once
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

void superconductor_wrapper();

#ifdef __cplusplus
}
void bcs();
void eliashberg();
void call_python_func(std::string filename, std::string function);
#endif

