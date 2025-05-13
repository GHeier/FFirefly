#include <iostream>
#include <cstddef>
#include "../../objects/vec.hpp"

extern "C" void print_vec_layout() {
    std::cout << "sizeof(Vec): " << sizeof(Vec) << "\n";
    std::cout << "offset x: " << offsetof(Vec, x) << "\n";
    std::cout << "offset y: " << offsetof(Vec, y) << "\n";
    std::cout << "offset z: " << offsetof(Vec, z) << "\n";
    std::cout << "offset w: " << offsetof(Vec, w) << "\n";
    std::cout << "offset area: " << offsetof(Vec, area) << "\n";
    std::cout << "offset dimension: " << offsetof(Vec, dimension) << "\n";
    std::cout << "offset n: " << offsetof(Vec, n) << "\n";
}

