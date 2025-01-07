#include <jlcxx.hpp>

int add(int a, int b) {
    return a + b;
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod) {
    mod.method("add", &add);
}
