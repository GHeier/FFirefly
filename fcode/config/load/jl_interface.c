#include <julia.h>
#include <linux/limits.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

void call_julia_func(const char *folder, const char *filename, const char* module, const char *func_name) {
    // Initialize Julia
    jl_init();

    // Set the path to the folder containing the Julia script
    char path[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
    path[len - 17] = '\0';

    strcat(path, "/fcode/");
    strcat(path, folder);

    // Add filename to path
    strcat(path, filename);

    // Add .jl to path
    char include_command[256];
    sprintf(include_command, "include(\"%s.jl\")", path);

    jl_eval_string(include_command);

    char using_command[256] = "";
    strcat(using_command, "using ");
    strcat(using_command, module);

    jl_eval_string(using_command);


    char command[256] = "";
    strcat(command, module);
    strcat(command, ".");
    strcat(command, func_name);
    strcat(command, "()");

    jl_eval_string(command);

    // Cleanup
    jl_atexit_hook(0);
}
