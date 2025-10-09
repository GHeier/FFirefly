#include <julia.h>
#include <linux/limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

// Helper function to check if a Julia module is already loaded
static bool is_module_loaded(const char *module_name) {
  if (!jl_is_initialized()) {
    return false;
  }

  jl_module_t *main_module = jl_main_module;
  jl_sym_t *module_sym = jl_symbol(module_name);
  jl_binding_t *binding = jl_get_binding(main_module, module_sym);

  if (binding == NULL || binding->value == NULL) {
    return false;
  }

  // Check if the binding value is actually a module
  if (jl_is_module(binding->value)) {
    return true;
  }

  return false;
}

void load_julia() {
  jl_init();

  char path[PATH_MAX];
  ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
  path[len - 16] = '\0';

  char activate_command[256];
  sprintf(activate_command, "Pkg.activate(\"%s/jlpkg/Firefly/\")", path);
  jl_eval_string("import Pkg");
  jl_eval_string(activate_command);
}

bool call_julia_func(const char *folder, const char *filename,
                     const char *module, const char *func_name) {
  // Initialize Julia
  //jl_init();

  // Set the path to the folder containing the Julia script
  char path[PATH_MAX];
  ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
  path[len - 16] = '\0';

  //char activate_command[256];
  //sprintf(activate_command, "Pkg.activate(\"%s/jlpkg/Firefly/\")", path);
  //jl_eval_string("import Pkg");
  //jl_eval_string(activate_command);

  strcat(path, "/src/");
  strcat(path, folder);

  // Add filename to path
  strcat(path, filename);

  // Only include the file if the module isn't already loaded
  if (!is_module_loaded(module)) {
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

    jl_eval_string(using_command);
  }

  jl_value_t *ret = NULL;

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

bool call_julia_func_bool(const char *folder, const char *filename,
                          const char *module, const char *func_name) {
  return call_julia_func(folder, filename, module, func_name);
}

int call_julia_func_int(const char *folder, const char *filename,
                        const char *module, const char *func_name) {
  jl_init();

  char path[PATH_MAX];
  ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
  path[len - 16] = '\0';

  char activate_command[256];
  sprintf(activate_command, "Pkg.activate(\"%s/jlpkg/Firefly/\")", path);
  jl_eval_string("import Pkg");
  jl_eval_string(activate_command);

  strcat(path, "/src/");
  strcat(path, folder);
  strcat(path, filename);

  // Only include the file if the module isn't already loaded
  if (!is_module_loaded(module)) {
    char include_command[256];
    sprintf(include_command, "include(\"%s.jl\")", path);
    jl_eval_string(include_command);

    if (jl_exception_occurred()) {
      printf("Error at include\n");
      jl_call2(jl_get_function(jl_base_module, "showerror"), jl_stderr_obj(),
               jl_exception_occurred());
      fprintf(stderr, "\n");
      jl_atexit_hook(0);
      return 0;
    }

    char using_command[256] = "";
    strcat(using_command, "using ");
    strcat(using_command, module);
    jl_eval_string(using_command);
  }

  char command[256] = "";
  strcat(command, module);
  strcat(command, ".");
  strcat(command, func_name);
  strcat(command, "()");

  jl_value_t *ret = jl_eval_string(command);
  if (jl_exception_occurred()) {
    printf("Julia function call error:\n");
    jl_call2(jl_get_function(jl_base_module, "showerror"), jl_stderr_obj(),
             jl_exception_occurred());
    fprintf(stderr, "\n");
    jl_atexit_hook(0);
    return 0;
  }

  int result = 0;
  jl_datatype_t *ret_type = (jl_datatype_t*)jl_typeof(ret);

  if (jl_is_int32(ret) || jl_is_int64(ret)) {
    result = (int)jl_unbox_int64(ret);
  } else {
    printf("Wrong type returned (expected int)\n");
  }

  jl_atexit_hook(0);
  return result;
}

float call_julia_func_float(const char *folder, const char *filename,
                            const char *module, const char *func_name) {
  jl_init();

  char path[PATH_MAX];
  ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
  path[len - 16] = '\0';

  char activate_command[256];
  sprintf(activate_command, "Pkg.activate(\"%s/jlpkg/Firefly/\")", path);
  jl_eval_string("import Pkg");
  jl_eval_string(activate_command);

  strcat(path, "/src/");
  strcat(path, folder);
  strcat(path, filename);

  // Only include the file if the module isn't already loaded
  if (!is_module_loaded(module)) {
    char include_command[256];
    sprintf(include_command, "include(\"%s.jl\")", path);
    jl_eval_string(include_command);

    if (jl_exception_occurred()) {
      printf("Error at include\n");
      jl_call2(jl_get_function(jl_base_module, "showerror"), jl_stderr_obj(),
               jl_exception_occurred());
      fprintf(stderr, "\n");
      jl_atexit_hook(0);
      return 0.0f;
    }

    char using_command[256] = "";
    strcat(using_command, "using ");
    strcat(using_command, module);
    jl_eval_string(using_command);
  }

  char command[256] = "";
  strcat(command, module);
  strcat(command, ".");
  strcat(command, func_name);
  strcat(command, "()");

  jl_value_t *ret = jl_eval_string(command);
  if (jl_exception_occurred()) {
    printf("Julia function call error:\n");
    jl_call2(jl_get_function(jl_base_module, "showerror"), jl_stderr_obj(),
             jl_exception_occurred());
    fprintf(stderr, "\n");
    jl_atexit_hook(0);
    return 0.0f;
  }

  float result = 0.0f;

  jl_datatype_t *ret_type = (jl_datatype_t*)jl_typeof(ret);
  if (ret_type == jl_float32_type) {
    result = jl_unbox_float32(ret);
  } else if (ret_type == jl_float64_type) {
    result = (float)jl_unbox_float64(ret);
  } else {
    printf("Wrong type returned (expected float)\n");
  }

  jl_atexit_hook(0);
  return result;
}

double call_julia_func_double(const char *folder, const char *filename,
                              const char *module, const char *func_name) {
  jl_init();

  char path[PATH_MAX];
  ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
  path[len - 16] = '\0';

  char activate_command[256];
  sprintf(activate_command, "Pkg.activate(\"%s/jlpkg/Firefly/\")", path);
  jl_eval_string("import Pkg");
  jl_eval_string(activate_command);

  strcat(path, "/src/");
  strcat(path, folder);
  strcat(path, filename);

  // Only include the file if the module isn't already loaded
  if (!is_module_loaded(module)) {
    char include_command[256];
    sprintf(include_command, "include(\"%s.jl\")", path);
    jl_eval_string(include_command);

    if (jl_exception_occurred()) {
      printf("Error at include\n");
      jl_call2(jl_get_function(jl_base_module, "showerror"), jl_stderr_obj(),
               jl_exception_occurred());
      fprintf(stderr, "\n");
      jl_atexit_hook(0);
      return 0.0;
    }

    char using_command[256] = "";
    strcat(using_command, "using ");
    strcat(using_command, module);
    jl_eval_string(using_command);
  }

  char command[256] = "";
  strcat(command, module);
  strcat(command, ".");
  strcat(command, func_name);
  strcat(command, "()");

  jl_value_t *ret = jl_eval_string(command);
  if (jl_exception_occurred()) {
    printf("Julia function call error:\n");
    jl_call2(jl_get_function(jl_base_module, "showerror"), jl_stderr_obj(),
             jl_exception_occurred());
    fprintf(stderr, "\n");
    jl_atexit_hook(0);
    return 0.0;
  }

  double result = 0.0;

  jl_datatype_t *ret_type = (jl_datatype_t*)jl_typeof(ret);
  if (ret_type == jl_float64_type) {
    result = jl_unbox_float64(ret);
  } else if (ret_type == jl_float32_type) {
    result = (double)jl_unbox_float32(ret);
  } else {
    printf("Wrong type returned (expected double)\n");
  }

  jl_atexit_hook(0);
  return result;
}

const char* call_julia_func_string(const char *folder, const char *filename,
                                   const char *module, const char *func_name) {
  jl_init();

  char path[PATH_MAX];
  ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
  path[len - 16] = '\0';

  char activate_command[256];
  sprintf(activate_command, "Pkg.activate(\"%s/jlpkg/Firefly/\")", path);
  jl_eval_string("import Pkg");
  jl_eval_string(activate_command);

  strcat(path, "/src/");
  strcat(path, folder);
  strcat(path, filename);

  // Only include the file if the module isn't already loaded
  if (!is_module_loaded(module)) {
    char include_command[256];
    sprintf(include_command, "include(\"%s.jl\")", path);
    jl_eval_string(include_command);

    if (jl_exception_occurred()) {
      printf("Error at include\n");
      jl_call2(jl_get_function(jl_base_module, "showerror"), jl_stderr_obj(),
               jl_exception_occurred());
      fprintf(stderr, "\n");
      jl_atexit_hook(0);
      return strdup("");
    }

    char using_command[256] = "";
    strcat(using_command, "using ");
    strcat(using_command, module);
    jl_eval_string(using_command);
  }

  char command[256] = "";
  strcat(command, module);
  strcat(command, ".");
  strcat(command, func_name);
  strcat(command, "()");

  jl_value_t *ret = jl_eval_string(command);
  if (jl_exception_occurred()) {
    printf("Julia function call error:\n");
    jl_call2(jl_get_function(jl_base_module, "showerror"), jl_stderr_obj(),
             jl_exception_occurred());
    fprintf(stderr, "\n");
    jl_atexit_hook(0);
    return strdup("");
  }

  const char *result = strdup("");

  if (jl_is_string(ret)) {
    const char *str = jl_string_ptr(ret);
    result = strdup(str);
  } else {
    printf("Wrong type returned (expected string)\n");
  }

  jl_atexit_hook(0);
  return result;
}
