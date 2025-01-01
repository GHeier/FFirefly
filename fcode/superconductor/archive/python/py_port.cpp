#include <iostream>
#include <fstream>
#include <math.h>
#include <python3.10/Python.h>

using namespace std;

// NOTE: Gets chi for cartesian momentum vector
double get_chi(double* q_vec, double T, double mu) {

    setenv("PYTHONPATH",".",1);

    Py_Initialize();
    
    PyObject* module_name = PyUnicode_FromString((char*)"py_funcs");
    PyObject* module_phys_funcs = PyImport_Import(module_name);
    PyObject* myFunction = PyObject_GetAttrString(module_phys_funcs,(char*)"chi");

    double qx = q_vec[0];
    double qy = q_vec[1];
    double qz = q_vec[2];

    PyObject* args = PyTuple_Pack(5,PyFloat_FromDouble(qx),PyFloat_FromDouble(qy), PyFloat_FromDouble(qz), PyFloat_FromDouble(T), PyFloat_FromDouble(mu));

    PyObject* myResult = PyObject_CallObject(myFunction, args);

    Py_DECREF(args);

    double result = PyFloat_AsDouble(myResult);

    //cout << result << endl;

    Py_DECREF(myFunction);
    Py_DECREF(module_phys_funcs);
    Py_DECREF(module_name);

    //Py_Finalize();

    return result;
}

double index_to_double(double i, int n) {
    return 3.1416*(i/(n-1));
}

void save_chi(double T, double mu) {
    ofstream file;
    file.open("chi_vals.dat");
    int n = 31;
    double magnitude = 0.1;
    for (double i = 0; i < n; i++) {
        cout << magnitude << endl;
        double q_vec[3] = {magnitude, 0.0, 0.0};
        double chi_val = get_chi(q_vec, T, mu);
        file << chi_val << endl;
        magnitude += 0.1;
    }
}

//int main() {
//
//    double q_vec [3] = { 2.0, 2.0, 2.0 };
//    double T = 0.25; 
//    double mu = 0.0;
//    //cout << get_chi(q_vec, T, mu);
//    save_chi(T, mu);
//
//}
