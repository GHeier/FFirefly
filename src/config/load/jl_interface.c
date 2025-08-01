#include <julia.h>
#include <linux/limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

bool call_julia_func(const char *folder, const char *filename,
                     const char *module, const char *func_name) {
  // Initialize Julia
  jl_init();

  // Set the path to the folder containing the Julia script
  char path[PATH_MAX];
  ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
  path[len - 16] = '\0';

  char activate_command[256];
  sprintf(activate_command, "Pkg.activate(\"%s/jlpkg/Firefly/\")", path);
  jl_eval_string("import Pkg");
  jl_eval_string(activate_command);

  strcat(path, "/src/");
  strcat(path, folder);

  // Add filename to path
  strcat(path, filename);

  // Add .jl to path
  char include_command[256];
  // sprintf(include_command, "Base.include(Main, \"%s.jl\")", path);
  sprintf(include_command, "include(\"%s.jl\")", path);

  jl_eval_string(include_command);
  // jl_value_t *ret = jl_eval_string(include_command);
  if (jl_exception_occurred()) {
    printf("Error at include\n");
    jl_call2(jl_get_function(jl_base_module, "showerror"), jl_stderr_obj(),
             jl_exception_occurred());
    fprintf(stderr, "\n");
    jl_atexit_hook(0);
    return false;
  }

  char using_command[256] = "";
  strcat(using_command, "using ");
  strcat(using_command, module);

  jl_value_t *ret = jl_eval_string(using_command);

  char command[256] = "";
  strcat(command, module);
  strcat(command, ".");
  strcat(command, func_name);
  strcat(command, "()");

  ret = jl_eval_string(command);
  if (jl_exception_occurred()) {
    printf("Julia function call error:\n");
    jl_call2(jl_get_function(jl_base_module, "showerror"), jl_stderr_obj(),
             jl_exception_occurred());
    fprintf(stderr, "\n");
    jl_atexit_hook(0);
    return false;
  }
   if (ret == jl_nothing) return true;
    jl_datatype_t *ret_type = (jl_datatype_t*)jl_typeof(ret);
    bool result = false;
    if (ret_type == jl_bool_type) {
        result = jl_unbox_bool(ret);
    }
    else {
        printf("Wrong type returned\n");
        exit(1);
    }
    return result;
  // return jl_unbox_bool(ret);

  // Cleanup
  jl_atexit_hook(0);
}
