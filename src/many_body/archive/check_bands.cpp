#include <iostream>
#include <cmath>

int main() {
    double t = 1.0, mu = -1.0;
    
    std::cout << "E(k) = -2t[cos(kx) + cos(ky)] - mu\n";
    std::cout << "With t=" << t << ", mu=" << mu << "\n\n";
    
    double pi = 3.14159265;
    
    std::cout << "k=(0,0):   E=" << (-2*t*(std::cos(0) + std::cos(0)) - mu) << "\n";
    std::cout << "k=(π,0):   E=" << (-2*t*(std::cos(pi) + std::cos(0)) - mu) << "\n";
    std::cout << "k=(π,π):   E=" << (-2*t*(std::cos(pi) + std::cos(pi)) - mu) << "\n";
    std::cout << "k=(π/2,π/2): E=" << (-2*t*(std::cos(pi/2) + std::cos(pi/2)) - mu) << "\n";
    
    std::cout << "\nFor Ef=0, we need states with E<0\n";
    std::cout << "Min E (at k=(0,0)): " << (-2*t*2 - mu) << "\n";
    std::cout << "Max E (at k=(π,π)): " << (-2*t*(-2) - mu) << "\n";
    
    return 0;
}
