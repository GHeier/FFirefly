#pragma once
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

void superconductor_wrapper();
void superconductor_tests();

#ifdef __cplusplus
}
void bcs();
void eliashberg();
void call_python_func(std::string filename, std::string function);

bool DOS_test();
#endif

