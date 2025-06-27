#include "DOS.h"
#include "../config/load/c_config.h"
#include "custom_DOS.hpp"

extern void dos_wrapper();
extern void surface_sum();

void DOS_spectrum() {
    if (!strcmp(c_method, "libtetrabz")) dos_wrapper();
    else {
        surface_sum();
    }
}

