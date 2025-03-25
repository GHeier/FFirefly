#pragma once

#ifdef __cplusplus
extern "C" {
#endif

void start_python();
void end_python();
void call_python_func(const char *folder, const char *filename, const char *function);

#ifdef __cplusplus
}
#endif

