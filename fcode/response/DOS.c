#include "DOS.h"
#include "../config/load/c_config.h"
#include "custom_DOS.hpp"

extern void dos_wrapper();
extern void custom_dos_call();

void DOS_spectrum() {
    if (!strcmp(c_method, "libtetrabz")) dos_wrapper();
    else if (!strcmp(c_method, "surface_sum")) custom_dos_call();
    else {
        printf("Method '%s' not recognized, exiting\n", c_method);
        exit(1);
    }
}

