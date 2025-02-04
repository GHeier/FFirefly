#include "response.h"
#include "../config/load/c_config.h"
#include "../config/load/jl_interface.h"

extern void polarization_wrapper();

void response_wrapper() {
    if (!strcmp(c_method, "libtetrabz"))
        polarization_wrapper();
    else if (!strcmp(c_method, "sparse_ir"))
        ir_wrapper();
    else
        printf("Method not found\n");
}

void ir_wrapper() {
    char* folder = "response/";
    char* filename = "sparse_ir_response";
    char *module = "response_ir";
    char* function = "get_ckio_ir";
    call_julia_func(folder, filename, module, function);
}
