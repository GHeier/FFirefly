#include "../config/load/py_interface.h"
#include <string>

using namespace std;


extern "C" void electron_number_wrapper() {
    string folder = "algorithms/";
    string filename = "mu_to_n";
    string function = "mu_vs_n";
    call_python_func(folder.c_str(), filename.c_str(), function.c_str());
}
