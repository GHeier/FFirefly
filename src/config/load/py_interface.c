#include <Python.h>
#include <linux/limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>

void start_python() {
    // Check if Python is already initialized (e.g., by PyCall)
    if (Py_IsInitialized()) {
        // Python already initialized, just return
        return;
    }

    // Ensure embedded Python uses the conda environment's site-packages
    // by setting PYTHONPATH before initialization
    const char* conda_prefix = getenv("CONDA_PREFIX");
    if (conda_prefix != NULL) {
        char pythonpath_env[PATH_MAX * 2];
        // Add conda environment's site-packages to PYTHONPATH
        snprintf(pythonpath_env, sizeof(pythonpath_env),
                "%s/lib/python3.11/site-packages:%s",
                conda_prefix, getenv("PYTHONPATH") ? getenv("PYTHONPATH") : "");
        setenv("PYTHONPATH", pythonpath_env, 1);
    }
    Py_Initialize();
}

void end_python() { Py_Finalize(); }
void call_python_func(const char *folder, const char *filename,
                      const char *function) {
  char path[PATH_MAX]; // Buffer to hold the executable path
  ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
  path[len - 16] = '\0';
  strcat(path, "/src/");
  strcat(path, folder);
  char pycommand[256];
   sprintf(pycommand, "import sys; sys.path.append('%s')", path);
  PyRun_SimpleString(pycommand);

  // Try loading your custom module
  const char *moduleName = filename;
  PyObject *pModule = PyImport_ImportModule(moduleName);
  if (!pModule) {
    PyErr_Print();
    fprintf(stderr, "Failed to load \"%s\"\n", moduleName);
    exit(1);
  }

  // Attempt to find and call a function from the module
  PyObject *pFunc = PyObject_GetAttrString(pModule, function);
  if (pFunc && PyCallable_Check(pFunc)) {
    PyObject *pResult = PyObject_CallObject(pFunc, NULL);
    if (pResult == NULL && PyErr_Occurred()) {
      PyErr_Print();
    }
    Py_XDECREF(pResult);
  } else {
    if (PyErr_Occurred())
      PyErr_Print();
    fprintf(stderr, "Cannot find function \"%s\"\n", function);
  }

  // Cleanup
  Py_XDECREF(pFunc);
  Py_DECREF(pModule);
}

bool call_python_func_bool(const char *folder, const char *filename,
                           const char *function) {
  char path[PATH_MAX]; // Buffer to hold the executable path
  ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
  path[len - 16] = '\0';
  strcat(path, "/src/");
  strcat(path, folder);
  char pycommand[256];
  sprintf(pycommand, "import sys; sys.path.append('%s')", path);
  PyRun_SimpleString(pycommand);

  // Try loading your custom module
  const char *moduleName = filename;
  PyObject *pModule = PyImport_ImportModule(moduleName);
  if (!pModule) {
    PyErr_Print();
    fprintf(stderr, "Failed to load \"%s\"\n", moduleName);
    return false;
  }

  // Attempt to find and call a function from the module
  PyObject *pFunc = PyObject_GetAttrString(pModule, function);
  bool result = false;

  if (pFunc && PyCallable_Check(pFunc)) {
    PyObject *pResult = PyObject_CallObject(pFunc, NULL);
    if (pResult == NULL && PyErr_Occurred()) {
      PyErr_Print();
      fprintf(stderr, "Error calling function \"%s\"\n", function);
    } else if (pResult != NULL) {
      // Check if the result is a boolean
      if (PyBool_Check(pResult)) {
        result = (pResult == Py_True);
      } else if (PyLong_Check(pResult)) {
        // Also accept integers (Python's bool is a subclass of int)
        result = (PyLong_AsLong(pResult) != 0);
      } else {
        fprintf(stderr, "Warning: Function \"%s\" did not return a boolean value\n", function);
        // Try to interpret as boolean anyway
        result = PyObject_IsTrue(pResult);
      }
    }
    Py_XDECREF(pResult);
  } else {
    if (PyErr_Occurred())
      PyErr_Print();
    fprintf(stderr, "Cannot find function \"%s\"\n", function);
  }

  // Cleanup
  Py_XDECREF(pFunc);
  Py_DECREF(pModule);

  return result;
}

int call_python_func_int(const char *folder, const char *filename,
                         const char *function) {
  char path[PATH_MAX];
  ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
  path[len - 16] = '\0';
  strcat(path, "/src/");
  strcat(path, folder);
  char pycommand[256];
  sprintf(pycommand, "import sys; sys.path.append('%s')", path);
  PyRun_SimpleString(pycommand);

  const char *moduleName = filename;
  PyObject *pModule = PyImport_ImportModule(moduleName);
  if (!pModule) {
    PyErr_Print();
    fprintf(stderr, "Failed to load \"%s\"\n", moduleName);
    return 0;
  }

  PyObject *pFunc = PyObject_GetAttrString(pModule, function);
  int result = 0;

  if (pFunc && PyCallable_Check(pFunc)) {
    PyObject *pResult = PyObject_CallObject(pFunc, NULL);
    if (pResult == NULL && PyErr_Occurred()) {
      PyErr_Print();
      fprintf(stderr, "Error calling function \"%s\"\n", function);
    } else if (pResult != NULL) {
      if (PyLong_Check(pResult)) {
        result = (int)PyLong_AsLong(pResult);
      } else {
        fprintf(stderr, "Warning: Function \"%s\" did not return an integer\n", function);
      }
    }
    Py_XDECREF(pResult);
  } else {
    if (PyErr_Occurred())
      PyErr_Print();
    fprintf(stderr, "Cannot find function \"%s\"\n", function);
  }

  Py_XDECREF(pFunc);
  Py_DECREF(pModule);
  return result;
}

float call_python_func_float(const char *folder, const char *filename,
                             const char *function) {
  char path[PATH_MAX];
  ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
  path[len - 16] = '\0';
  strcat(path, "/src/");
  strcat(path, folder);
  char pycommand[256];
  sprintf(pycommand, "import sys; sys.path.append('%s')", path);
  PyRun_SimpleString(pycommand);

  const char *moduleName = filename;
  PyObject *pModule = PyImport_ImportModule(moduleName);
  if (!pModule) {
    PyErr_Print();
    fprintf(stderr, "Failed to load \"%s\"\n", moduleName);
    return 0.0f;
  }

  PyObject *pFunc = PyObject_GetAttrString(pModule, function);
  float result = 0.0f;

  if (pFunc && PyCallable_Check(pFunc)) {
    PyObject *pResult = PyObject_CallObject(pFunc, NULL);
    if (pResult == NULL && PyErr_Occurred()) {
      PyErr_Print();
      fprintf(stderr, "Error calling function \"%s\"\n", function);
    } else if (pResult != NULL) {
      if (PyFloat_Check(pResult)) {
        result = (float)PyFloat_AsDouble(pResult);
      } else if (PyLong_Check(pResult)) {
        result = (float)PyLong_AsLong(pResult);
      } else {
        fprintf(stderr, "Warning: Function \"%s\" did not return a float\n", function);
      }
    }
    Py_XDECREF(pResult);
  } else {
    if (PyErr_Occurred())
      PyErr_Print();
    fprintf(stderr, "Cannot find function \"%s\"\n", function);
  }

  Py_XDECREF(pFunc);
  Py_DECREF(pModule);
  return result;
}

double call_python_func_double(const char *folder, const char *filename,
                               const char *function) {
  char path[PATH_MAX];
  ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
  path[len - 16] = '\0';
  strcat(path, "/src/");
  strcat(path, folder);
  char pycommand[256];
  sprintf(pycommand, "import sys; sys.path.append('%s')", path);
  PyRun_SimpleString(pycommand);

  const char *moduleName = filename;
  PyObject *pModule = PyImport_ImportModule(moduleName);
  if (!pModule) {
    PyErr_Print();
    fprintf(stderr, "Failed to load \"%s\"\n", moduleName);
    return 0.0;
  }

  PyObject *pFunc = PyObject_GetAttrString(pModule, function);
  double result = 0.0;

  if (pFunc && PyCallable_Check(pFunc)) {
    PyObject *pResult = PyObject_CallObject(pFunc, NULL);
    if (pResult == NULL && PyErr_Occurred()) {
      PyErr_Print();
      fprintf(stderr, "Error calling function \"%s\"\n", function);
    } else if (pResult != NULL) {
      if (PyFloat_Check(pResult)) {
        result = PyFloat_AsDouble(pResult);
      } else if (PyLong_Check(pResult)) {
        result = (double)PyLong_AsLong(pResult);
      } else {
        fprintf(stderr, "Warning: Function \"%s\" did not return a double\n", function);
      }
    }
    Py_XDECREF(pResult);
  } else {
    if (PyErr_Occurred())
      PyErr_Print();
    fprintf(stderr, "Cannot find function \"%s\"\n", function);
  }

  Py_XDECREF(pFunc);
  Py_DECREF(pModule);
  return result;
}

// Note: Caller must free the returned string
const char* call_python_func_string(const char *folder, const char *filename,
                                    const char *function) {
  char path[PATH_MAX];
  ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
  path[len - 16] = '\0';
  strcat(path, "/src/");
  strcat(path, folder);
  char pycommand[256];
  sprintf(pycommand, "import sys; sys.path.append('%s')", path);
  PyRun_SimpleString(pycommand);

  const char *moduleName = filename;
  PyObject *pModule = PyImport_ImportModule(moduleName);
  if (!pModule) {
    PyErr_Print();
    fprintf(stderr, "Failed to load \"%s\"\n", moduleName);
    return strdup("");
  }

  PyObject *pFunc = PyObject_GetAttrString(pModule, function);
  const char *result = strdup("");

  if (pFunc && PyCallable_Check(pFunc)) {
    PyObject *pResult = PyObject_CallObject(pFunc, NULL);
    if (pResult == NULL && PyErr_Occurred()) {
      PyErr_Print();
      fprintf(stderr, "Error calling function \"%s\"\n", function);
    } else if (pResult != NULL) {
      if (PyUnicode_Check(pResult)) {
        const char *str = PyUnicode_AsUTF8(pResult);
        if (str != NULL) {
          result = strdup(str);
        }
      } else {
        fprintf(stderr, "Warning: Function \"%s\" did not return a string\n", function);
      }
    }
    Py_XDECREF(pResult);
  } else {
    if (PyErr_Occurred())
      PyErr_Print();
    fprintf(stderr, "Cannot find function \"%s\"\n", function);
  }

  Py_XDECREF(pFunc);
  Py_DECREF(pModule);
  return result;
}
