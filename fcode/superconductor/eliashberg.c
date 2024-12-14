#include <Python.h>
#include <stdio.h>
#include <string.h>
#include "../config/load/c_config.h"

void eliashberg() {
    Py_Initialize();
    PyRun_SimpleString("import sys; sys.path.append('/home/g/Research/bcs_diagonalization/fcode/superconductor/')");

    // Try loading your custom module
    const char *moduleName = "eliashberg1";
    PyObject *pModule = PyImport_ImportModule(moduleName);
    if (!pModule) {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n", moduleName);
        exit(1);
    }

    // Attempt to find and call a function from the module
    PyObject *pFunc = PyObject_GetAttrString(pModule, "eliashberg");
    if (pFunc && PyCallable_Check(pFunc)) {
        PyObject *pResult = PyObject_CallObject(pFunc, NULL);
        if (pResult == NULL && PyErr_Occurred()) {
            PyErr_Print();
        }
        Py_XDECREF(pResult);
    } else {
        if (PyErr_Occurred())
            PyErr_Print();
        fprintf(stderr, "Cannot find function 'eliashberg'\n");
    }

    // Cleanup
    Py_XDECREF(pFunc);
    Py_DECREF(pModule);
    Py_Finalize();
}
