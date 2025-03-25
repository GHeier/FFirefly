#include <Python.h>
#include <linux/limits.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

void start_python() {
    Py_Initialize();
}

void end_python() {
    Py_Finalize();
}
void call_python_func(const char *folder, const char *filename, const char *function) {
    char path[PATH_MAX]; // Buffer to hold the executable path
    ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
    path[len - 17] = '\0';
    strcat(path, "/fcode/");
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


